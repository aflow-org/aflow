// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#include "aflow_apl.h"

static const double MIN_FREQ_THRESHOLD = -0.1;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////
  DOSCalculator::DOSCalculator() {
    free();
  }

  DOSCalculator::DOSCalculator(PhononCalculator& pc, const xoption &aplopts) {
    free();
    _pc = &pc;
    _pc_set = true;
    initialize(aplopts);
  }

  DOSCalculator::DOSCalculator(const DOSCalculator& that) {
    if (this != &that) free();
    copy(that);
  }

  DOSCalculator& DOSCalculator::operator=(const DOSCalculator& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void DOSCalculator::copy(const DOSCalculator& that) {
    if (this == &that) return;
    _pc = that._pc;
    _pc_set = that._pc_set;
    _bzmethod = that._bzmethod;
    _qpoints = that._qpoints;
    _freqs = that._freqs;
    _minFreq = that._minFreq;
    _maxFreq = that._maxFreq;
    _stepDOS = that._stepDOS;
    _halfStepDOS = that._halfStepDOS;
    _bins = that._bins;
    _dos = that._dos;
    _idos = that._idos;
    _eigen = that._eigen;
    _projectedDOS = that._projectedDOS;
    _projections = that._projections;
    _temperature = that._temperature;
    _system = that._system;
  }

  DOSCalculator::~DOSCalculator() {
    free();
  }

  void DOSCalculator::free() {
    _qpoints.clear();
    //_qweights.clear();  OBSOLETE ME20190423
    _freqs.clear();
    _bins.clear();
    _dos.clear();
    _idos.clear();  //ME20190614
    _eigen.clear();  //ME20190624
    _projectedDOS.clear(); //ME20190614
    _projections.clear();  //ME20190624
    _bzmethod = "";
    _temperature = 0.0;  //ME20190614
    _minFreq = AUROSTD_NAN;
    _maxFreq = AUROSTD_NAN;
    _stepDOS = 0.0;
    _halfStepDOS = 0.0;
    _system = "";
    _pc = NULL;
    _pc_set = false;
  }

  void DOSCalculator::clear(PhononCalculator& pc) {
    free();
    _pc = &pc;
    _pc_set = true;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::initialize(const xoption &aplopts){
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    //AS20200312 BEGIN: now initialization parameters are passed using xoption
    _bzmethod = aplopts.getattachedscheme("DOSMETHOD");

    // projections
    if (aplopts.flag("DOS_PROJECT")) {
      if (aplopts.flag("DOS_CART") || aplopts.flag("DOS_FRAC")) {
        string projscheme = "";
        if (aplopts.flag("DOS_CART")){
          projscheme = aplopts.getattachedscheme("DOSPROJECTIONS_CART");
        } else {
          projscheme = aplopts.getattachedscheme("DOSPROJECTIONS_FRAC");
        }
        vector<string> tokens;
        aurostd::string2tokens(projscheme, tokens, "; ");
        vector<double> proj;
        for (uint i = 0; i < tokens.size(); i++) {
          aurostd::string2tokens(tokens[i], proj, ", ");
          _projections.push_back(aurostd::vector2xvector<double>(proj));
        }
      } else {
        xvector<double> proj(3);
        _projections.push_back(proj);
      }
    }

    if ((_projections.size() > 0) && aplopts.flag("DOS_FRAC")) {
      for (uint p = 0; p < _projections.size(); p++) {
        _projections[p] = _pc->getInputCellStructure().f2c * _projections[p];
      }
    }
    //AS20200312 END
    _system = _pc->_system;

    if (!_pc->getSupercell().isConstructed()) {
      message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }

    QMesh& _qm = _pc->getQMesh();
    if (!_qm.initialized()) {
      message = "q-point mesh is not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if (_projections.size() == 0) _qm.makeIrreducible();

    calculateFrequencies();
  }

  // ///////////////////////////////////////////////////////////////////////////

  //CO START
  void DOSCalculator::calculateInOneThread(int iqp) {
    _freqs[iqp] = _pc->getFrequency(_qpoints[iqp], apl::THZ | apl::ALLOW_NEGATIVE, _eigen[iqp]);  //ME20190624
  }
  //CO END

  //////////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calculateFrequencies() {
    // Get q-points for which to calculate the frequencies
    _qpoints = _pc->getQMesh().getIrredQPointsCPOS();

    //CO START
    string message = "Calculating frequencies for the phonon DOS.";
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());

    // Prepare storage
    _freqs.clear();
    xvector<double> zero(_pc->getNumberOfBranches());
    for (uint i = 0; i < _qpoints.size(); i++)
      _freqs.push_back(zero);
    _eigen.resize(_qpoints.size(), xmatrix<xcomplex<double> >(_pc->getNumberOfBranches(), _pc->getNumberOfBranches()));  // M20190624

    int nqps = (int) _qpoints.size();
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_pc->getNCPUs());
    std::function<void(int)> fn = std::bind(&DOSCalculator::calculateInOneThread, this, std::placeholders::_1);
    xt.run(nqps, fn);
#else
    for (int i = 0; i < nqps; i++) calculateInOneThread(i);
#endif
    //CO END

    //if freq > MIN_FREQ_THRESHOLD considerd as +ve freq [PN]
    for (uint i = 0; i < _freqs.size(); i++) {
      for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
        if ((_freqs[i][j] < 0.00) && (_freqs[i][j] > MIN_FREQ_THRESHOLD)) _freqs[i][j] = 0.00;
      }
    }
    //PN END

    // Get min and max values
    _maxFreq = -1.0;
    _minFreq = 0.0;
    for (uint i = 0; i < _freqs.size(); i++) {
      for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
        if (_freqs[i](j) > _maxFreq) _maxFreq = _freqs[i](j);
        if (_freqs[i](j) < _minFreq) _minFreq = _freqs[i](j);
      }
    }
    _maxFreq += 1.0;
    if (_minFreq < MIN_FREQ_THRESHOLD) _minFreq -= 1.0;
    else _minFreq = 0.0;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190614 - added integrated DOS
  void DOSCalculator::smearWithGaussian(vector<double>& dos, vector<double>& idos, double h, double sigma) {
    // Construct table for gaussian function
    int ng = (int)(6.0 * sigma / h + 1.0);
    double fact = 1.0 / (sqrt(2.0 * M_PI) * sigma);
    vector<double> gauss, igauss;
    double gnorm = 0.0;
    for (int ig = -ng; ig <= ng; ig++) {
      double eg = ig * h;
      double arg = eg * eg / (sigma * sigma) / 2.0;
      gauss.push_back(fact * exp(-arg));
      igauss.push_back(fact * erf(-arg));
      gnorm += gauss.back();
    }

    // Norm gauss table to one
    gnorm *= h;
    for (int ig = -ng; ig <= ng; ig++) {
      gauss[ig + ng] /= gnorm;
      igauss[ig + ng] /= gnorm;
    }

    // Prepare new dos
    vector<double> newdos;
    for (uint i = 0; i < dos.size(); i++) {
      newdos.push_back(0.0);
    }
    vector<double> newidos(newdos.size(), 0.0);

    // Convolute...
    for (int ie = 0; ie < (int)dos.size(); ie++) {
      double wt = dos[ie] * h;
      double wti = idos[ie] * h;
      for (int jg = -ng; jg <= ng; jg++) {
        int je = ie + jg;
        if (je < 0) continue;
        if (je >= (int)dos.size()) continue;
        newdos[je] += gauss[jg + ng] * wt;
        newidos[je] += igauss[jg + ng] * wti;
      }
    }

    //
    dos.clear();
    for (uint i = 0; i < newdos.size(); i++) {
      dos.push_back(newdos[i]);
    }
    newdos.clear();
    gauss.clear();
  }

  //ME20200203 - DOS can now be calculated within any frequency range
  //ME20210927 - Added VERBOSE option
  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calc(int USER_DOS_NPOINTS, bool VERBOSE) {
    calc(USER_DOS_NPOINTS, 0.0, _minFreq, _maxFreq, VERBOSE);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calc(int USER_DOS_NPOINTS, double USER_DOS_SMEAR, bool VERBOSE) {
    calc(USER_DOS_NPOINTS, USER_DOS_SMEAR, _minFreq, _maxFreq, VERBOSE);
  }

  // ///////////////////////////////////////////////////////////////////////////
  void DOSCalculator::calc(int USER_DOS_NPOINTS, double USER_DOS_SMEAR,
      double fmin, double fmax, bool VERBOSE) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    // Check parameters
    if (aurostd::isequal(fmax, fmin, _FLOAT_TOL_)) {
      message = "Frequency range of phonon DOS is nearly zero.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    } else if (fmin > fmax) {
      double tmp = fmax;
      fmax = fmin;
      fmin = tmp;
    }
    // Calculate steps
    _stepDOS = (fmax - fmin) / (double)USER_DOS_NPOINTS;
    _halfStepDOS = 0.5 * _stepDOS;

    // Clear old stuff
    _dos.clear();
    _idos.clear();  //ME20190614
    _bins.clear();

    // Prepare storagearrays
    for (int k = 0; k < USER_DOS_NPOINTS; k++) {
      _dos.push_back(0);
      _idos.push_back(0);  //ME20190614
      _bins.push_back(fmin + k * _stepDOS + _halfStepDOS);
    }
    //ME20190624
    if (_projections.size() > 0)
      _projectedDOS.resize(_pc->getInputCellStructure().atoms.size(),
          vector<vector<double> >(_projections.size(), vector<double>(USER_DOS_NPOINTS)));

    // Perform the raw specific calculation by method
    //ME20190423 START
    //ME20200321 - Added logger output
    //rawCalc(USER_DOS_NPOINTS);  OBSOLETE
    if (_bzmethod == "LT") {
      if (VERBOSE) {
        message = "Calculating phonon DOS using the linear tetrahedron method.";
        pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
      }
      calcDosLT();
    } else if (_bzmethod == "RS") {
      if (VERBOSE) {
        message = "Calculating phonon DOS using the root sampling method.";
        pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
      }
      calcDosRS();
    }
    //ME20190423 END

    // Smooth DOS by gaussians
    if (USER_DOS_SMEAR > 1E-6)
      smearWithGaussian(_dos, _idos, _stepDOS, USER_DOS_SMEAR);  //ME20190614

    // Normalize to number of branches
    // ME20210927 - Only normalize when inside full frequency spectrum
    if ((fmin - _minFreq <= _ZERO_TOL_) && (_maxFreq - fmax <= _ZERO_TOL_)) {
      double sum = 0.0;
      for (int k = 0; k < USER_DOS_NPOINTS; k++)
        sum += _dos[k];
      sum /= _pc->getNumberOfBranches();

      for (int k = 0; k < USER_DOS_NPOINTS; k++) {
        _dos[k] /= (sum * _stepDOS);
        //_idos[k] /= (sum * _stepDOS);  //ME20190614  // OBSOLETE - ME20200228 - not necessary for iDOS
      }
    }
  }

  //ME20190423 START

  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calcDosRS() {
    for (uint k = 0; k < _bins.size(); k++) {
      for (uint i = 0; i < _freqs.size(); i++) {
        for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
          if ((_freqs[i][j] > (_bins[k] - _halfStepDOS)) &&
              (_freqs[i][j] < (_bins[k] + _halfStepDOS))) {
            _dos[k] += (double) _pc->getQMesh().getWeights()[i];
          }
        }
      }
    }
    //ME20190614 - calculate integrated DOS
    _idos[0] = _dos[0];
    for (uint k = 1; k < _dos.size(); k++) {
      _idos[k] = _idos[k-1] + _dos[k];
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190614 - added integrated DOS
  //ME20190625 - rearranged and added projected DOS
  //ME20200213 - added atom-projected DOS
  void DOSCalculator::calcDosLT() {
    QMesh& _qm = _pc->getQMesh();

    // Procompute projections for each q-point and branch to save time
    uint nproj = _projections.size();
    vector<xvector<double> > proj_norm(nproj, xvector<double>(3));
    for (uint p = 0; p < nproj; p++) {
      if (!aurostd::iszero(_projections[p])) proj_norm[p] = _projections[p]/aurostd::modulus(_projections[p]);
    }
    const xstructure& xstr = _pc->getInputCellStructure();
    uint natoms = xstr.atoms.size();
    vector<vector<vector<vector<double> > > > parts;
    if (nproj > 0) {
      // Precompute eigenvector projections
      xcomplex<double> eig;
      parts.assign(_qm.getnQPs(), vector<vector<vector<double> > >(_pc->getNumberOfBranches(), vector<vector<double> >(nproj, vector<double>(natoms, 0))));
      for (int q = 0; q < _qm.getnQPs(); q++) {
        for (uint br = 0; br < _pc->getNumberOfBranches(); br++) {
          int ibranch = br + _freqs[0].lrows;
          for (uint p = 0; p < nproj; p++) {
            for (uint at = 0; at < natoms; at++) {
              if (aurostd::iszero(proj_norm[p])) {  // zero-vector = atom-projected DOS
                for (int i = 1; i < 4; i++) {
                  parts[q][br][p][at] += aurostd::magsqr(_eigen[q][3*at + i][ibranch]);
                }
              } else {
                eig.re = 0.0;
                eig.im = 0.0;
                for (int i = 1; i < 4; i++) eig += proj_norm[p][i] * _eigen[q][3*at + i][ibranch];
                parts[q][br][p][at] = aurostd::magsqr(eig);
              }
            }
          }
        }
      }
    }

    vector<double> f(4);
    double max_freq = _bins.back() + _halfStepDOS;
    double min_freq = _bins.front() - _halfStepDOS;
    vector<vector<int> > tet_corners;
    _qm.generateTetrahedra();
    if (nproj == 0) {
      _qm.makeIrreducibleTetrahedra();
      tet_corners = _qm.getIrreducibleTetrahedraIbzqpt();
    } else {
      tet_corners = _qm.getTetrahedra();
    }
    for (int itet = 0; itet < _qm.getnIrredTetrahedra(); itet++) {
      double weightVolumeTetrahedron = _qm.getWeightTetrahedron(itet) * _qm.getVolumePerTetrahedron();
      const vector<int>& corners = tet_corners[itet];
      for (int ibranch = _freqs[0].lrows; ibranch <= _freqs[0].urows; ibranch++) {
        for (int icorner = 0; icorner < 4; icorner++) {
          f[icorner] = _freqs[corners[icorner]][ibranch];
        }
        std::sort(f.begin(), f.end());
        double fmin = f[0];
        double fmax = f[3];
        if (fmax > max_freq) continue;
        if (fmin < min_freq) continue;

        int kstart = (int) ((fmin - min_freq)/_stepDOS) + 1;
        if (kstart < 0) kstart = 0;
        if (kstart > (int) _bins.size()) kstart = _bins.size() - 1;
        int kstop = (int) ((fmax - min_freq)/_stepDOS) + 1;
        if (kstop < 0) kstop = 0;
        if (kstop < (int) _bins.size()) kstop = _bins.size() - 1;

        double f21 = f[1]  - f[0];
        double f31 = f[2]  - f[0];
        double f32 = f[2]  - f[1];
        double f41 = f[3]  - f[0];
        double f42 = f[3]  - f[1];
        double f43 = f[3]  - f[2];
        double cc12 = 3.0 * weightVolumeTetrahedron/(f21 * f31 * f41);
        double cc23 = weightVolumeTetrahedron/(f31 * f41);
        double cc23a = 3.0 * f21 * cc23;
        double cc23b = 6.0 * cc23;
        double cc23c = -3.0 * cc23 * (f31 + f42)/(f32 * f42);
        double cc34 = 3.0 * weightVolumeTetrahedron/(f41 * f42 * f43);

        double fbin, dos, part;
        int br = ibranch - _freqs[0].lrows;
        for (int k = kstart; k <= kstop; k++) {
          //ME20200203 - Use bins to accommodate different frequency range
          fbin = _bins[k]; // _minFreq + k * _stepDOS + _halfStepDOS;
          dos = 0.0;
          if ((f[0] <= fbin) && (fbin <= f[1])) {
            double df = fbin - f[0];
            dos = cc12 * df * df;
            _idos[k] += cc12 * (df * df * df)/3.0;
          } else if ((f[1] < fbin) && (fbin <= f[2])) {
            double df = fbin - f[1];
            dos = cc23a + cc23b * df + cc23c * df * df;
            _idos[k] += cc23a * f21/3.0 + 3.0 * cc23 * f21 * df + 3.0 * cc23 * df * df + cc23c * (df * df * df)/3.0;
          } else if ((f[2] < fbin) && (fbin <= f[3])) {
            double df = f[3] - fbin;
            dos = cc34 * df * df;
            _idos[k] += weightVolumeTetrahedron + cc34 * (df * df * df)/3.0;
          } else if (f[3] < fbin) {
            _idos[k] += weightVolumeTetrahedron;
          }
          _dos[k] += dos;
          for (uint p = 0; p < nproj; p++) {
            for (uint at = 0; at < natoms; at++) {
              part = 0.0;
              //ME20200320 - due to the big nesting level, loop unrolling
              // leads to a huge speed-up (x2 or more)
              //[OBSOLETE] for (int icorner = 0; icorner < 4; icorner++) {
              //[OBSOLETE]   part += parts[corners[icorner]][br][p][at];
              //[OBSOLETE] }
              part += parts[corners[0]][br][p][at];
              part += parts[corners[1]][br][p][at];
              part += parts[corners[2]][br][p][at];
              part += parts[corners[3]][br][p][at];
              _projectedDOS[at][p][k] += 0.25 * dos * part;
            }
          }
        }
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190423 END

  void DOSCalculator::writePDOS(const string& directory) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    // Write PHDOS file
    //CO START
    //ofstream outfile("PDOS",ios_base::out);
    stringstream outfile;
    //if( !outfile.is_open() )
    //{
    //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
    //}
    //CO END

    string filename = directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDOS_FILE;
    double factorTHz2Raw = _pc->getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2rcm = _pc->getFrequencyConversionFactor(apl::RAW, apl::RECIPROCAL_CM);
    double factorRaw2meV = _pc->getFrequencyConversionFactor(apl::RAW, apl::MEV);

    message = "Writing phonon density of states into file " + aurostd::CleanFileName(filename) + "."; //ME20181226
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    //outfile << "############### ############### ############### ###############" << std::endl;
    outfile << "#    f(THz)      1/lambda(cm-1)      E(meV)          pDOS      " << std::endl;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    for (uint k = 0; k < _dos.size(); k++) {
      outfile << setw(15) << _bins[k] << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2rcm << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2meV << " "
        << setw(15) << _dos[k] << std::endl;
    }

    //CO START
    aurostd::stringstream2file(outfile, filename); //ME20181226
    if (!aurostd::FileExist(filename)) { //ME20181226
      message = "Cannot open output file " + filename + "."; //ME20181226
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
    }
    //outfile.clear();
    //outfile.close();
    //CO END
  }

  //ME20190614 - writes phonon DOS in DOSCAR format
  void DOSCalculator::writePHDOSCAR(const string& directory) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHDOSCAR_FILE);
    message = "Writing phonon density of states into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    stringstream doscar;
    xDOSCAR xdos = createDOSCAR();
    doscar << xdos;
    aurostd::stringstream2file(doscar, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    // OBSOLETE ME20191219 - PHPOSCAR is already written in KBIN::RunPhonons_APL
    // if (xdos.partial) {  // Write PHPOSCAR if there are projected DOS
    //   filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHPOSCAR_FILE);
    //   xstructure xstr = _pc->getInputCellStructure();
    //   xstr.is_vasp5_poscar_format = true;
    //   stringstream poscar;
    //   poscar << xstr;
    //   aurostd::stringstream2file(poscar, filename);
    //   if (!aurostd::FileExist(filename)) {
    //     string message = "Cannot open output file " + filename + ".";
    //     throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    //   }
    // }
  }

  xDOSCAR DOSCalculator::createDOSCAR() const {
    if (!_pc_set) {
      string message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    xDOSCAR xdos;
    xdos.spin = 0;
    // Header values
    xdos.number_atoms = _pc->getInputCellStructure().atoms.size();
    xdos.partial = (_projectedDOS.size() > 0);
    xdos.Vol = GetVolume(_pc->getInputCellStructure())/xdos.number_atoms;
    xvector<double> lattice;
    lattice[1] = _pc->getInputCellStructure().a * 1E-10;
    lattice[2] = _pc->getInputCellStructure().b * 1E-10;
    lattice[3] = _pc->getInputCellStructure().c * 1E-10;
    xdos.lattice = lattice;
    xdos.POTIM = 0.5E-15;
    xdos.temperature = _temperature;
    xdos.carstring = "PHON";
    xdos.title = _system;

    // Data
    double factorTHz2Raw = _pc->getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2meV = _pc->getFrequencyConversionFactor(apl::RAW, apl::MEV);
    double conv = factorTHz2Raw * factorRaw2meV/1000;
    //ME20200203 - use _bins instead of _minFreq in case the DOS was calculated
    // using different frequency ranges
    xdos.energy_max = _bins.back() * conv;
    xdos.energy_min = _bins[0] * conv;
    xdos.number_energies = _dos.size();
    xdos.Efermi = 0.0;  // phonon DOS have no Fermi energy
    xdos.venergy = aurostd::vector2deque(_bins);
    xdos.venergyEf = xdos.venergy;  //ME20200324
    for (uint i = 0; i < xdos.number_energies; i++) xdos.venergy[i] *= conv;
    xdos.viDOS.resize(1);
    xdos.viDOS[0] = aurostd::vector2deque(_idos);
    deque<deque<deque<deque<double> > > > vDOS;
    //ME20190625
    if (_projections.size() > 0) {
      vDOS.resize(xdos.number_atoms + 1, deque<deque<deque<double> > >(_projections.size() + 1, deque<deque<double> >(1, deque<double>(xdos.number_energies, 0.0))));
    } else {
      vDOS.resize(1, deque<deque<deque<double> > >(1, deque<deque<double> >(1, deque<double>(xdos.number_energies, 0.0))));
    }

    vDOS[0][0][0] = aurostd::vector2deque<double>(_dos);

    //ME20190624 - projected DOS
    if (_projections.size() > 0) {
      for (uint at = 0; at < xdos.number_atoms; at++) {
        for (uint p = 0; p < _projections.size(); p++) {
          vDOS[at + 1][p + 1][0] = aurostd::vector2deque(_projectedDOS[at][p]);
          //ME20200321 - Add to totals for consistency
          for (uint e = 0; e < xdos.number_energies; e++) {
            vDOS[0][p + 1][0][e] += vDOS[at + 1][p + 1][0][e];
            vDOS[at + 1][0][0][e] += vDOS[at + 1][p + 1][0][e];
          }
        }
      }
    }
    xdos.vDOS = vDOS;

    return xdos;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20200108 - added const
  const vector<double>& DOSCalculator::getBins() const {
    return _bins;
  }

  const vector<double>& DOSCalculator::getDOS() const {
    return _dos;
  }

  //ME20200210
  const vector<double>& DOSCalculator::getIDOS() const {
    return _idos;
  }

  //AS20200512
  const vector<xvector<double> >& DOSCalculator::getFreqs() const {
    return _freqs;
  }

  bool DOSCalculator::hasImaginaryFrequencies() const {
    return (_minFreq < MIN_FREQ_THRESHOLD ? true : false);
  }

  // ME20210927
  double DOSCalculator::getMinFreq() const {
    return _minFreq;
  }

  // ME20210927
  double DOSCalculator::getMaxFreq() const {
    return _maxFreq;
  }

  // ME20210927
  const xstructure& DOSCalculator::getInputStructure() const {
    return _pc->getInputCellStructure();
  }

  // ME20210927
  uint DOSCalculator::getNumberOfBranches() const {
    return _pc->getNumberOfBranches();
  }

  // ///////////////////////////////////////////////////////////////////////////
  //PN START
  void DOSCalculator::writePDOS(string path, string ex)  //[PN]
  {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    //CO START
    // Write PHDOS file
    stringstream outfile;
    //CO END

    double factorTHz2Raw = _pc->getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2rcm = _pc->getFrequencyConversionFactor(apl::RAW, apl::RECIPROCAL_CM);
    double factorRaw2meV = _pc->getFrequencyConversionFactor(apl::RAW, apl::MEV);

    string filename = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDOS_FILE; //ME20181226
    message = "Writing phonon density of states into file " + filename + "."; //ME20181226
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    //outfile << "############### ############### ############### ###############" << std::endl;
    outfile << "#    f(THz)      1/lambda(cm-1)      E(meV)          pDOS      " << std::endl;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    for (uint k = 0; k < _dos.size(); k++) {
      outfile << setw(15) << _bins[k] << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2rcm << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2meV << " "
        << setw(15) << _dos[k] << std::endl;
    }

    //CO START
    string file = path + "/" + filename + "." + ex; //ME20181226
    aurostd::stringstream2file(outfile, file);
    if (!aurostd::FileExist(file)) {
      message = "Cannot open output file " + filename + "."; //ME20181226
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    //CO END
  }
  //PN END
  // ///////////////////////////////////////////////////////////////////////////

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
