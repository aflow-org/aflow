// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *                 Aflow MARCO ESTERS - Duke University 2020               *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters based on work by Pinku Nath

#include "aflow_apl.h"

static const xcomplex<double> iONE(0.0, 1.0);

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  AtomicDisplacements::AtomicDisplacements() {
    free();
  }

  AtomicDisplacements::AtomicDisplacements(PhononCalculator& pc) {
    free();
    _pc = &pc;
    _pc_set = true;
  }

  AtomicDisplacements::AtomicDisplacements(const AtomicDisplacements& that) {
    if (this != &that) free();
    copy(that);
  }

  AtomicDisplacements& AtomicDisplacements::operator=(const AtomicDisplacements& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  AtomicDisplacements::~AtomicDisplacements() {
    free();
  }

  void AtomicDisplacements::copy(const AtomicDisplacements& that) {
    if (this == &that) return;
    _eigenvectors = that._eigenvectors;
    _frequencies = that._frequencies;
    _displacement_matrices = that._displacement_matrices;
    _displacement_modes = that._displacement_modes;
    _pc = that._pc;
    _pc_set = that._pc_set;
    _qpoints = that._qpoints;
    _temperatures = that._temperatures;
  }

  void AtomicDisplacements::free() {
    _eigenvectors.clear();
    _frequencies.clear();
    _displacement_matrices.clear();
    _displacement_modes.clear();
    _pc = NULL;
    _pc_set = false;
    _qpoints.clear();
    _temperatures.clear();
  }

  void AtomicDisplacements::clear(PhononCalculator& pc) {
    free();
    _pc = &pc;
    _pc_set = true;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             EIGENVECTORS                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //calculateEigenvectors/////////////////////////////////////////////////////
  // Calculates eigenvectors and rearranges them from the APL matrix format
  // into atom-resolved vectors.
  void AtomicDisplacements::calculateEigenvectors() {
    _eigenvectors.clear();
    _frequencies.clear();
    int nq = (int) _qpoints.size();
    if (nq == 0) return;
    uint natoms = _pc->getInputCellStructure().atoms.size();
    uint nbranches = _pc->getNumberOfBranches();
    _eigenvectors.resize(nq, vector<vector<xvector<xcomplex<double> > > >(nbranches, vector<xvector<xcomplex<double> > >(natoms, xvector<xcomplex<double> >(3))));
    _frequencies.resize(nq, vector<double> (nbranches, 0.0));
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_pc->getNCPUs());
    std::function<void(int)> fn = std::bind(&AtomicDisplacements::calculateEigenvectorsInThread, this, std::placeholders::_1);
    xt.run(nq, fn);
#else
    for (int i = 0; i < nq; ++i) calculateEigenvectorsInThread(i);
#endif
  }

  void AtomicDisplacements::calculateEigenvectorsInThread(int i) {
    uint nbranches = _pc->getNumberOfBranches();
    uint natoms = _pc->getInputCellStructure().atoms.size();
    xvector<double> freq(nbranches);
    xmatrix<xcomplex<double> > eig(nbranches, nbranches, 1, 1);
    freq = _pc->getFrequency(_qpoints[i].cpos, apl::THZ, eig);
    for (uint br = 0; br < nbranches; br++) {
      _frequencies[i][br] = freq[br + 1];
      for (uint at = 0; at < natoms; at++) {
        for (int j = 1; j < 4; j++) {
          _eigenvectors[i][br][at][j] = eig[3 * at + j][br + 1];
        }
      }
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             DISPLACEMENTS                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //calculateMeanSquareDisplacements//////////////////////////////////////////
  // Calculates the complex mean square displacement matrices for a range of
  // temperatures. For the formula, see e.g. A. A. Mardudin "Theory of
  // Lattice Dynamics in the Harmonic Approximation", eq. 2.4.23 and 2.4.24.
  // Units are Angstrom^2.
  void AtomicDisplacements::calculateMeanSquareDisplacements(double Tstart, double Tend, double Tstep) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    _qpoints.clear();
    _temperatures.clear();

    QMesh& _qm = _pc->getQMesh();
    if (!_qm.initialized()) {
      message = "q-point mesh is not initialized.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    _qpoints = _pc->getQMesh().getPoints();

    if (Tstart > Tend) {
      message = "Tstart cannot be higher than Tend.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    for (double T = Tstart; T <= Tend; T += Tstep) _temperatures.push_back(T);

    calculateMeanSquareDisplacementMatrices();
  }

  void AtomicDisplacements::calculateMeanSquareDisplacementMatrices() {
    string message = "Calculating mean square displacement matrices.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    _displacement_matrices.clear();
    _displacement_modes.clear();
    calculateEigenvectors();

    uint ntemps = _temperatures.size();
    uint natoms = _pc->getInputCellStructure().atoms.size();
    uint nq = _qpoints.size();
    uint nbranches = _pc->getNumberOfBranches();

    _displacement_matrices.resize(ntemps, vector<xmatrix<xcomplex<double> > >(natoms, xmatrix<xcomplex<double> >(3, 3)));
    vector<xmatrix<xcomplex<double> > > outer(natoms, xmatrix<xcomplex<double> >(3, 3));
    vector<double> masses(natoms);
    for (uint at = 0; at < natoms; at++) masses[at] = AMU2KILOGRAM * _pc->getSupercell().getAtomMass(_pc->getSupercell().pc2scMap(at));

    // Factor 2pi necessary because frequencies are raw
    double prefactor = (PLANCKSCONSTANT_hbar * Hz2THz * std::pow(1e10, 2)/(2 * PI * (double) _qpoints.size()));
    for (uint q = 0; q < nq; q++) {
      for (uint br = 0; br < nbranches; br++) {
        if (_frequencies[q][br] > _FLOAT_TOL_) {
          for (uint at = 0; at < natoms; at++) {
            outer[at] = aurostd::outer_product(_eigenvectors[q][br][at], conj(_eigenvectors[q][br][at]));
          }
          for (uint t = 0; t < ntemps; t++) {
            double prefactor_T = prefactor * ((0.5 + getOccupationNumber(_temperatures[t], _frequencies[q][br]))/_frequencies[q][br]);
            for (uint at = 0; at < natoms; at++) {
              // Add element-wise (much faster)
              for (int i = 1; i < 4; i++) {
                for (int j = 1; j < 4; j++) {
                  _displacement_matrices[t][at][i][j].re += (prefactor_T/masses[at]) * outer[at][i][j].re;
                  _displacement_matrices[t][at][i][j].im += (prefactor_T/masses[at]) * outer[at][i][j].im;
                }
              }
            }
          }
        }
      }
    }
  }

  //calculateModeDisplacements////////////////////////////////////////////////
  // Calculates the displacements for phonon modes along specific q-points.
  // Units are 1/sqrt(amu).
  void AtomicDisplacements::calculateModeDisplacements(const vector<xvector<double> >& qpts, bool coords_are_fractional) {
    if (!_pc_set) {
      string message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    _qpoints.clear();
    uint nq = qpts.size();
    _qpoints.resize(nq);
    if (coords_are_fractional) {
      xmatrix<double> f2c = trasp(ReciprocalLattice(_pc->getInputCellStructure().lattice));
      for (uint q = 0; q < nq; q++) {
        _qpoints[q].fpos = qpts[q];
        _qpoints[q].cpos = f2c * qpts[q];
      }
    } else {
      xmatrix<double> c2f = inverse(trasp(ReciprocalLattice(_pc->getInputCellStructure().lattice)));
      for (uint q = 0; q < nq; q++) {
        _qpoints[q].cpos = qpts[q];
        _qpoints[q].fpos = c2f * qpts[q];
      }
    }
    calculateModeDisplacements();
  }

  void AtomicDisplacements::calculateModeDisplacements() {
    _displacement_matrices.clear();
    _displacement_modes.clear();
    _temperatures.clear();
    calculateEigenvectors();

    uint nq = _qpoints.size();
    uint nbranches = _pc->getNumberOfBranches();
    uint natoms = _pc->getInputCellStructure().atoms.size();
    _displacement_modes.resize(nq, vector<vector<xvector<xcomplex<double> > > >(nbranches, vector<xvector<xcomplex<double> > >(natoms)));

    vector<double> masses(natoms);
    for (uint at = 0; at < natoms; at++) masses[at] = _pc->getSupercell().getAtomMass(_pc->getSupercell().pc2scMap(at));

    for (uint q = 0; q < nq; q++) {
      for (uint br = 0; br < nbranches; br++) {
        for (uint at = 0; at < natoms; at++) {
          _displacement_modes[q][br][at] = _eigenvectors[q][br][at]/sqrt(masses[at]);
        }
      }
    }
  }

  //getOccupationNumber///////////////////////////////////////////////////////
  // Calculates the phonon numbers for a specific frequency and temperature
  // using Bose-Einstein statistics.
  double AtomicDisplacements::getOccupationNumber(double T, double f) {
    if (T < _FLOAT_TOL_) return 0.0;
    else return (1.0/(exp(BEfactor_h_THz * f/T) - 1));
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            GETTER FUNCTIONS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  const vector<double>& AtomicDisplacements::getTemperatures() const {
    return _temperatures;
  }

  const vector<vector<xmatrix<xcomplex<double> > > >& AtomicDisplacements::getDisplacementMatrices() const {
    return _displacement_matrices;
  }

  // The diagonal of the displacement matrix are the x-, y- and z-components.
  // They are always real.
  vector<vector<xvector<double> > > AtomicDisplacements::getDisplacementVectors() const {
    vector<vector<xvector<double> > >  disp_vec;
    uint ntemps = _displacement_matrices.size();
    if (ntemps > 0) {
      uint natoms = _displacement_matrices[0].size();
      disp_vec.resize(ntemps, vector<xvector<double> >(natoms, xvector<double>(3)));
      for (uint t = 0; t < ntemps; t++) {
        for (uint at = 0; at < natoms; at++) {
          for (int i = 1; i < 4; i++) {
            disp_vec[t][at][i] = _displacement_matrices[t][at][i][i].re;
          }
        }
      }
    }
    return disp_vec;
  }

  const vector<vector<vector<xvector<xcomplex<double> > > > >& AtomicDisplacements::getModeDisplacements() const {
    return _displacement_modes;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                FILE I/O                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //writeMeanSquareDisplacementsToFile////////////////////////////////////////
  // Writes the mean square displacement vectors to a file.
  void AtomicDisplacements::writeMeanSquareDisplacementsToFile(string filename) {
    if (!_pc_set) {
      string message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    filename = aurostd::CleanFileName(filename);
    string message = "Writing mean square displacements into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    vector<vector<xvector<double> > > disp_vec = getDisplacementVectors();
    stringstream output;
    string tag = "[APL_DISPLACEMENTS]";

    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << tag << "SYSTEM=" << _pc->_system << std::endl;
    output << tag << "START" << std::endl;
    output << "#" << std::setw(9) << "T (K)" << setw(15) << "Species"
      << std::setw(15) << "x (A^2)" << std::setw(15) << "y (A^2)" << std::setw(15) << "z (A^2)" << std::endl;
    output << std::fixed << std::showpoint;
    for (uint t = 0; t < _temperatures.size(); t++) {
      output << std::setw(10) << std::setprecision(2) << _temperatures[t];
      for (uint at = 0; at < disp_vec[t].size(); at++) {
        output << (at==0?std::setw(15):std::setw(25)) << _pc->getInputCellStructure().atoms[at].cleanname;
        for (int i = 1; i < 4; i++) output << std::setw(15) << std::setprecision(8) << disp_vec[t][at][i];
        output << std::endl;
      }
    }
    output << tag << "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not write to file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //writeSceneFileXcrysden////////////////////////////////////////////////////
  // Writes an animated XCRYSDEN structure file that can be used to create a
  // gif or mpeg of a phonon mode displacement.
  void AtomicDisplacements::writeSceneFileXcrysden(string filename, const xstructure& scell, const vector<vector<vector<double> > >& disp, int nperiods) {
    if (!_pc_set) {
      string message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    filename = aurostd::CleanFileName(filename);
    string message = "Writing atomic displacements in XCRYSDEN format into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(),*_pc->getOSS());

    uint nsteps = disp.size();
    uint natoms = scell.atoms.size();

    stringstream output;
    output << "ANIMSTEPS " << (nperiods * nsteps) << std::endl;
    output << "CRYSTAL" << std::endl;
    output << "PRIMVEC" << std::endl;
    output << scell.lattice << std::endl;

    int step = 1;
    for (int i = 0; i < nperiods; i++) {
      for (uint j = 0; j < nsteps; j++) {
        output << "PRIMCOORD " << step << std::endl;
        output << std::setw(4) << natoms << " 1" << std::endl;
        for (uint at = 0; at < natoms; at++) {
          output << std::setw(5) << scell.atoms[at].atomic_number;
          for (uint k = 0; k < 6; k++) {
            output << std::setw(15) << std::fixed << std::setprecision(8) << disp[j][at][k];
          }
          output << std::endl;
        }
        step++;
      }
    }

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not write to file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //writeSceneFileVsim////////////////////////////////////////////////////////
  // Writes displacements into a V_sim-formatted file that can be visualized
  // by V_sim or ASCII-phonons.
  void AtomicDisplacements::writeSceneFileVsim(string filename, const xstructure& xstr_projected,
      const vector<vector<vector<xvector<xcomplex<double> > > > >& displacements) {
    if (!_pc_set) {
      string message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    filename = aurostd::CleanFileName(filename);
    string message = "Writing atomic displacements in V_SIM format into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());

    stringstream output;
    // Lattice
    output << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[1][1]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[2][1]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[2][2] << std::endl;
    output << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[3][1]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[3][2]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[3][3] << std::endl;
    // Atoms
    uint natoms = xstr_projected.atoms.size();
    for (uint at = 0; at < natoms; at++) {
      for (int i = 1; i < 4; i++) output << std::setw(15) << std::setprecision(8) << xstr_projected.atoms[at].cpos[i];
      output << std::setw(4) << xstr_projected.atoms[at].cleanname << std::endl;
    }

    for (uint q = 0; q < displacements.size(); q++) {
      for (uint br = 0; br < displacements[q].size(); br++) {
        output << "#metaData: qpt=[";
        for (int i = 1; i < 4; i++) output << _qpoints[q].fpos[i] << ";";
        output << _frequencies[q][br] << "\\" << std::endl;
        for (uint at = 0; at < natoms; at++) {
          output << "#;";
          for (int i = 1; i < 4; i++) output << displacements[q][br][at][i].re << ";";
          for (int i = 1; i < 4; i++) {
            output << displacements[q][br][at][i].im;
            if (i < 3) output << ";";
            else output << "\\" << std::endl;
          }
        }
        output << "# ]" << std::endl;
      }
    }

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not write to file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INTERFACE                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //createAtomicDisplacementSceneFile/////////////////////////////////////////
  // Interfaces with the command line to create a phonon visualization file.
  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow, ostream& oss) {
    ofstream mf("/dev/null");
    createAtomicDisplacementSceneFile(vpflow, mf, oss);
  }

  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow, ofstream& mf, ostream& oss) {
    string message = "";

    // Parse command line options
    string directory = vpflow.getattachedscheme("ADISP::DIRECTORY");
    if (directory.empty()) directory = "./";
    else directory = aurostd::CleanFileName(directory + "/");

    // Format
    string format = aurostd::toupper(vpflow.getattachedscheme("ADISP::FORMAT"));
    if (format.empty()) format = aurostd::toupper(DEFAULT_APL_ADISP_SCENE_FORMAT);
    string allowed_formats_str = "XCRYSDEN,V_SIM";
    vector<string> allowed_formats;
    aurostd::string2tokens(allowed_formats_str, allowed_formats, ",");
    if (!aurostd::WithinList(allowed_formats, format)) {
      message = "Unrecognized format " + format + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // Amplitude
    string amplitude_str = vpflow.getattachedscheme("ADISP::AMPLITUDE");
    double amplitude = 0.0;
    if (amplitude_str.empty()) amplitude = DEFAULT_APL_ADISP_AMPLITUDE;
    else amplitude = aurostd::string2utype<double>(amplitude_str);
    if (amplitude < _FLOAT_TOL_) {
      message = "Amplitude must be positive.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // Number of steps per period
    string nsteps_str = vpflow.getattachedscheme("ADISP::STEPS");
    int nsteps = 0;
    if (format != "V_SIM") {
      if (nsteps_str.empty()) nsteps = DEFAULT_APL_ADISP_NSTEPS;
      else nsteps = aurostd::string2utype<int>(nsteps_str);
      if (nsteps < 1) {
        message = "Number of steps must be a positive integer";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
    }

    // Number of periods
    string nperiods_str = vpflow.getattachedscheme("ADISP::PERIODS");
    int nperiods = 0;
    if (format != "V_SIM") {
      if (nperiods_str.empty()) nperiods = DEFAULT_APL_ADISP_NPERIODS;
      else nperiods = aurostd::string2utype<int>(nperiods_str);
      if (nperiods < 1) {
        message = "Number of periods must be a positive integer";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
    }

    // Range/supercell
    string supercell_str = vpflow.getattachedscheme("ADISP::SUPERCELL");
    xvector<int> sc_dim(3); sc_dim.set(1);
    if (format != "V_SIM") {
      if (supercell_str.empty()) {
        supercell_str = "1x1x1";
      } else {
        vector<int> tokens;
        aurostd::string2tokens(supercell_str, tokens, "xX");
        if (tokens.size() == 3) {
          sc_dim = aurostd::vector2xvector(tokens);
        } else {
          message = "Broken supercell format.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        }
      }
    }

    // Branches
    string branches_str = vpflow.getattachedscheme("ADISP::BRANCHES");
    vector<int> branches;
    if (format != "V_SIM") {
      if (branches_str.empty()) {
        message = "No branches selected. Displacements will be calculated for all modes.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, directory, mf, oss);
      } else {
        aurostd::string2tokens(branches_str, branches, ",");
        // Branch index is based on 1, not 0
        for (uint br = 0; br < branches.size(); br++) branches[br] -= 1;
      }
    }

    // q-points
    string qpoints_str = vpflow.getattachedscheme("ADISP::QPOINTS");
    vector<xvector<double> > qpoints;
    if (qpoints_str.empty()) {
      message = "No q-points given.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
    } else {
      vector<string> tokens;
      aurostd::string2tokens(qpoints_str, tokens, ",");
      xvector<double> q(3);
      if (tokens.size() % 3 == 0) {
        for (uint i = 0; i < tokens.size(); i += 3) {
          for (int j = 0; j < 3; j++) q[j + 1] = aurostd::string2utype<double>(tokens[i + j]);
          qpoints.push_back(q);
        }
      } else {
        message = "Broken q-points format.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
    }

    // Initialize phonon calculator
    string statefile = directory + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_STATE_FILE;
    PhononCalculator pc(mf, oss);
    pc.setDirectory(directory);
    pc.initialize_supercell(statefile);
    pc.awake();
    apl::Supercell& sc_pcalc = pc.getSupercell();
    // Must project to primitive or the vibrations will be incorrect
    if (!sc_pcalc.projectToPrimitive()) {
      message = "Could not project to the primitive structure.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    // Check branches
    int nbr = (int) branches.size();
    if (nbr == 0) {
      nbr = (int) pc.getNumberOfBranches();
      for (int br = 0; br < nbr; br++) branches.push_back(br);
    } else {
      int nbranches = pc.getNumberOfBranches();
      for (int br = 0; br < nbr; br++) {
        if (branches[br] >= nbranches || branches[br] < 0) {
          message = "Index " + aurostd::utype2string<int>(branches[br] + 1) + " out of range.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
        }
      }
    }

    // Done with setup - calculate displacements
    AtomicDisplacements ad(pc);
    ad.calculateModeDisplacements(qpoints);

    if (format == "XCRYSDEN") {
      // Create supercell for the scene file
      Supercell scell(mf);
      scell.initialize(pc.getSupercell().getOriginalStructure(), false);
      scell.build(sc_dim, false);
      if (!scell.projectToPrimitive()) {
        message = "Could not project to primitive structure.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      vector<vector<vector<double> > > disp;
      for (uint q = 0; q < qpoints.size(); q++) {
        for (uint br = 0; br < branches.size(); br++) {
          disp = ad.createDisplacementsXcrysden(scell, amplitude, q, br, nsteps);
          stringstream filename;
          filename << directory <<  DEFAULT_APL_FILE_PREFIX << "displacements_q";
          for (int i = 1; i < 4; i++) filename << "_" << qpoints[q][i];
          filename << "_b_" << (branches[br] + 1) << ".axsf";
          ad.writeSceneFileXcrysden(filename.str(), scell.getSupercellStructure(), disp, nperiods);
        }
      }
    } else if (format == "V_SIM") {
      xstructure xstr_oriented;
      vector<vector<vector<xvector<xcomplex<double> > > > > displacements_oriented;
      ad.getOrientedDisplacementsVsim(xstr_oriented, displacements_oriented, amplitude);
      string filename = directory + DEFAULT_APL_FILE_PREFIX + "displacements.ascii";
      ad.writeSceneFileVsim(filename, xstr_oriented, displacements_oriented);
    }
  }

  //createDisplacementsXcrysden///////////////////////////////////////////////
  // Creates the displacements in a format that can be used to output into an
  // XCRYSDEN animated structure file. Since XCRYSDEN cannot perform the time
  // propagation itself, this function needs to create each time step.
  vector<vector<vector<double> > > AtomicDisplacements::createDisplacementsXcrysden(const Supercell& scell,
      double amplitude, int q, int br, int nsteps) {
    const xvector<double>& qpt = _qpoints[q].cpos;
    const vector<xvector<xcomplex<double> > >& disp = _displacement_modes[q][br];

    // Generate a list of atoms that are inside the desired sphere
    const xstructure& scell_str = scell.getSupercellStructure();
    uint natoms = scell_str.atoms.size();

    // Calculate original displacements
    const vector<int>& sc2pcMap = scell._sc2pcMap;
    const vector<int>& pc2scMap = scell._pc2scMap;

    // Calculate original displacements
    vector<vector<vector<double> > > displacements(nsteps, vector<vector<double> >(natoms, vector<double>(6, 0.0)));
    vector<xvector<xcomplex<double> > > displacements_orig(natoms, xvector<xcomplex<double> >(3));
    int at_pc = -1, at_eq_sc = -1;
    xvector<double> dist_scell(3);
    for (uint at = 0; at < natoms; at++) {
      at_pc = sc2pcMap[at];
      at_eq_sc = pc2scMap[at_pc];
      dist_scell = scell_str.atoms[at].cpos - scell_str.atoms[at_eq_sc].cpos;
      displacements_orig[at] = amplitude * exp(iONE * scalar_product(qpt, dist_scell)) * disp[at_pc];
      for (int i = 0; i < 3; i++) {
        displacements[0][at][i] = scell_str.atoms[at].cpos[i + 1];
        displacements[0][at][i + 3] = displacements_orig[at][i + 1].re;
      }
    }

    // Calculate displacements for the rest of the period
    xcomplex<double> phase;
    xvector<xcomplex<double> > disp_step;
    for (int s = 1; s < nsteps; s++) {
      phase = exp(-(2.0 * PI * (double) s/(double) nsteps) * iONE);
      for (uint at = 0; at < natoms; at++) {
        disp_step = phase * displacements_orig[at];
        for (int i = 0; i < 3; i++) {
          displacements[s][at][i] = displacements[s - 1][at][i] + displacements[s - 1][at][i + 3];
          displacements[s][at][i + 3] = disp_step[i + 1].re;
        }
      }
    }

    return displacements;
  }

  //getOrientedDisplacementsVsim//////////////////////////////////////////////
  // Creates displacements in the V_sim format. V_sim expects structures to
  // be projected such that a points along the x-axis and b points along the
  // x-y-plane, so the structure and the displacements needs to be proejcted.
  void AtomicDisplacements::getOrientedDisplacementsVsim(xstructure& xstr_oriented,
      vector<vector<vector<xvector<xcomplex<double> > > > >& displacements, double amplitude) {
    // Project the structure as required by V_sim
    const xstructure& xstr_orig = _pc->getInputCellStructure();
    xstr_oriented.clear();

    // Project lattice
    xvector<double> params = Getabc_angles(xstr_orig.lattice, RADIANS);
    xmatrix<double> lattice(3, 3);
    double cosalpha = cos(params[4]);
    double cosbeta = cos(params[5]);
    double cosgamma = cos(params[6]);
    double singamma = sin(params[6]);
    double l32 = (2* cosalpha - 2 * cosbeta * cosgamma)/(2 * singamma);
    lattice[1][1] = params[1];
    lattice[2][1] = params[2] * cosgamma;
    lattice[2][2] = params[2] * singamma;
    lattice[3][1] = params[3] * cosbeta;
    lattice[3][2] = params[3] * l32;
    lattice[3][3] = params[3] * sqrt(1 - std::pow(cosbeta, 2) - std::pow(l32, 2));
    xstr_oriented.lattice = lattice;
    xstr_oriented.f2c = trasp(lattice);
    xstr_oriented.c2f = inverse(xstr_oriented.f2c);

    // Project positions
    uint natoms = xstr_orig.atoms.size();
    for (uint i = 0; i < natoms; i++) {
      _atom at = xstr_orig.atoms[i];
      at.cpos = xstr_oriented.f2c * at.fpos;
      xstr_oriented.atoms.push_back(at);
    }

    // Project displacements
    uint nq = _qpoints.size();
    uint nbr = _pc->getNumberOfBranches();
    displacements.clear();
    displacements.resize(nq, vector<vector<xvector<xcomplex<double> > > >(nbr, vector<xvector<xcomplex<double> > >(natoms)));

    xmatrix<double> transf = amplitude * xstr_oriented.f2c * xstr_orig.c2f;
    for (uint q = 0; q < nq; q++) {
      for (uint br = 0; br < nbr; br++) {
        for (uint at = 0; at < natoms; at++) {
          displacements[q][br][at] = transf * _displacement_modes[q][br][at];
        }
      }
    }
  }

}  // namespace apl


// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow PINKU NATH - Duke University 2014-2016             *
// *            Aflow MARCO ESTERS - Duke University 2020                    *
// *                                                                         *
// ***************************************************************************
