// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// This file contains the classes PhononDispersionCalculator and PathBuilder

#include "aflow_apl.h"

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                        PhononDispersionCalculator                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  PhononDispersionCalculator::PhononDispersionCalculator() {
    free();
  }

  PhononDispersionCalculator::PhononDispersionCalculator(PhononCalculator& pc) {
    free();
    _pc = &pc;
    _pc_set = true;
    _system = _pc->_system;  //ME20190614
  }

  PhononDispersionCalculator::PhononDispersionCalculator(const PhononDispersionCalculator& that) {
    if (this != &that) free();
    copy(that);
  }

  PhononDispersionCalculator& PhononDispersionCalculator::operator=(const PhononDispersionCalculator& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void PhononDispersionCalculator::copy(const PhononDispersionCalculator& that) {
    if (this == &that) return;
    _frequencyFormat = that._frequencyFormat;
    _freqs = that._freqs;
    _pc = that._pc;
    _pc_set = that._pc_set;
    _pb = that._pb;
    _qpoints = that._qpoints;
    _system = that._system;
    _temperature = that._temperature;
  }

  PhononDispersionCalculator::~PhononDispersionCalculator() {
    free();
  }

  void PhononDispersionCalculator::free() {
    _qpoints.clear();
    _freqs.clear();
    _frequencyFormat = apl::NONE;
    _pc = NULL;
    _pc_set = false;
    _pb.clear();
    _system = "";
    _temperature = 0.0;  //ME20190614
  }

  void PhononDispersionCalculator::clear(PhononCalculator& pc) {
    free();
    _pc = &pc;
    _pc_set = true;
  }

  //////////////////////////////////////////////////////////////////////////////

  void PhononDispersionCalculator::initPathCoords(  //CO20180406
      const string& USER_DC_INITCOORDS,
      const string& USER_DC_INITLABELS,
      int USER_DC_NPOINTS, 
      bool CARTESIAN_COORDS) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    if(USER_DC_INITCOORDS.empty() || USER_DC_INITLABELS.empty()) {
      string message = "Inputs are empty.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }
    _pb.defineCustomPoints(USER_DC_INITCOORDS,USER_DC_INITLABELS,_pc->getSupercell(),CARTESIAN_COORDS);
    _pb.setDensity(USER_DC_NPOINTS);
    //_qpoints = _pb.getPath(); // Get points // OBSOLETE ME20190429 - this function should just define points; there is no path to set or get
  }

  void PhononDispersionCalculator::initPathLattice(const string& USER_DC_INITLATTICE,int USER_DC_NPOINTS){
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    string lattice = USER_DC_INITLATTICE;
    if (lattice.empty()) {
      xstructure a(_pc->getInputCellStructure());
      //CO START
      if (a.bravais_lattice_variation_type == "") {
        if (a.spacegroup == "") {
          if(a.space_group_ITC<1 || a.space_group_ITC>230){a.space_group_ITC = a.SpaceGroup_ITC();} //if (a.space_group_ITC == 0) { //CO20180214 - if not set then it could be 32767 //[CO20200106 - close bracket for indenting]}
        a.spacegroup = GetSpaceGroupName(a.space_group_ITC) + " #" + aurostd::utype2string(a.space_group_ITC);  //will break here if spacegroup is bad
        }
        // Use PLATON to get space group number if user did not get it...
        //a.platon2sg(_PLATON_P_EQUAL_DEFAULT,    //CO
        //	      _PLATON_P_EXACT_DEFAULT,
        //	      _PLATON_P_ANG_DEFAULT,
        //	      _PLATON_P_D1_DEFAULT,
        //	      _PLATON_P_D2_DEFAULT,
        //            _PLATON_P_D3_DEFAULT);
        //if( a.spacegroup.empty() )
        //throw apl::APLRuntimeError("apl::PhononDispersionCalculator::initPath(); The PLATON call to get spacegroup number failed. You have to specify it by DCINITSG in "+_AFLOWIN_);

        //vector<string> tokens;
        //aurostd::string2tokens(a.spacegroup,tokens,"#");
        //int spacegroupNumber = aurostd::string2utype<int>(tokens[1]);
        //tokens.clear();

        //lattice = LATTICE_Lattice_Variation_SpaceGroup(spacegroupNumber,_pc->getInputCellStructure());
        //lattice = LATTICE::SpaceGroup2LatticeVariation(spacegroupNumber,_pc->getInputCellStructure());
        //AS+DX20210602 BEGIN
        // lattice = LATTICE::SpaceGroup2LatticeVariation(a.space_group_ITC, a);
        // Sometimes SpaceGroup2LatticeVariation might not detect a correct
        // lattice variation, that is consistent with the space group number:
        // AS encountered a problem when MCLC is returned instead of MCLC1
        a.GetExtendedCrystallographicData();
        lattice = a.bravais_lattice_variation_type;
        //AS+DX20210602 END
      } else {
        lattice = a.bravais_lattice_variation_type;
      }
      message = "The phonon dispersion curves will be generated for lattice variation " + lattice + ".";
      pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(),*_pc->getOSS());
    }
    //CO END

    // cerr << "LATTICE=" << lattice << std::endl;
    // Suck point definition from the electronic structure part of AFLOW...
    _pb.takeAflowElectronicPath(lattice,_pc->getSupercell());             //CO20180406
    //_pc->getInputCellStructure(),        //CO20180406
    //_pc->getSuperCellStructure());       //CO20180406

    _pb.setDensity(USER_DC_NPOINTS);
    _qpoints = _pb.getPath(); // Get points
  }

  //////////////////////////////////////////////////////////////////////////////

  void PhononDispersionCalculator::setPath(const string& USER_DC_OWNPATH) {
    // Get user's path...
    if (!USER_DC_OWNPATH.empty()) {
      if (USER_DC_OWNPATH.find('|') != string::npos) {
        //ME20190614 START
        // This breaks "mixed" paths such as G-X-W-L|K-U (interprets as G-X|W-L|K-U
        // _qpoints = _pb.getPath(apl::PathBuilder::COUPLE_POINT_MODE, USER_DC_OWNPATH);
        vector<string> tokens, tokens_pt;
        aurostd::string2tokens(USER_DC_OWNPATH, tokens, "-");
        string path;
        if (tokens[0].find('|') != string::npos) {
          string message = "Cannot have | in the first path coordinate";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
        } else {
          path = tokens[0];
        }
        for (uint i = 1; i < tokens.size(); i++) {
          if ((tokens[i].find('|') != string::npos) || (i == tokens.size() - 1)) {
            path += '-' + tokens[i];
          } else {
            path += '-' + tokens[i] + '|' + tokens[i];
          }
        }
        _qpoints = _pb.getPath(apl::PathBuilder::COUPLE_POINT_MODE, path);
        //ME20190614 END
      } else {
        _qpoints = _pb.getPath(apl::PathBuilder::SINGLE_POINT_MODE, USER_DC_OWNPATH);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  void PhononDispersionCalculator::calculateInOneThread(int iqp) {
    //ME20200206 - get direction for q-points near Gamma for non-analytical correction
    // or the discontinuity due to LO-TO splitting is not accurately captured.
    if (_pc->isPolarMaterial() && (aurostd::modulus(_qpoints[iqp]) < 0.005)) {
      int npts = _pb.getDensity() + 1;
      int i = iqp/npts;
      xvector<double> qpoint_nac = _qpoints[i * npts] - _qpoints[(i + 1) * npts - 1];
      _freqs[iqp] = _pc->getFrequency(_qpoints[iqp], qpoint_nac, _frequencyFormat);
    } else {
      _freqs[iqp] = _pc->getFrequency(_qpoints[iqp], _frequencyFormat);
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  void PhononDispersionCalculator::calc(const IPCFreqFlags frequencyFormat) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    // Save
    _frequencyFormat = frequencyFormat;

    // Maybe there was some error and the list of q-points is empty, hence bye-bye...
    if (_qpoints.empty()) {
      //throw apl::APLRuntimeError("There are no points for calculation.");
      message = "There are no points for the calculation.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Compute frequencies for each q-point
    message = "Calculating frequencies for the phonon dispersion.";
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());

    // Prepare storage
    _freqs.clear();
    xvector<double> zero(_pc->getNumberOfBranches());
    for (uint i = 0; i < _qpoints.size(); i++)
      _freqs.push_back(zero);

    int nqps = (int) _qpoints.size();
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(_pc->getNCPUs());
    std::function<void(int)> fn = std::bind(&PhononDispersionCalculator::calculateInOneThread, this, std::placeholders::_1);
    xt.run(nqps, fn);
#else
    for (int i = 0; i < nqps; i++) calculateInOneThread(i);
#endif
  }

  //////////////////////////////////////////////////////////////////////////////

  void PhononDispersionCalculator::writePDIS(const string& directory) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDIS_FILE); //ME20181226
    message = "Writing dispersion curves into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());

    //CO START
    stringstream outfile;
    //CO END

    // Write header
    outfile << "# Phonon dispersion curves calculated by Aflow" << std::endl;
    outfile << "#" << std::endl;
    outfile << "# <system>    \"" << _system << "\"" << std::endl;  //ME20190614 - use system name, not structure title
    outfile << "#" << std::endl;
    outfile << "# <units>     " << _frequencyFormat << std::endl;
    outfile << "# <nbranches> " << _pc->getNumberOfBranches() << std::endl;
    outfile << "# <npoints>   " << _freqs.size() << std::endl;
    outfile << "# <nsubpathp> " << _pb.getDensity() + 1 << std::endl;
    outfile << "#" << std::endl;

    // Write table of label points
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    double x = 0.0;
    double wholePathLength = _pb.getPathLength();
    std::map<double, string> labelMap;
    for (uint i = 1; i < _pb.getPointSize(); i++) {
      outfile << "# <label>     " << x << " "
        << setw(5) << _pb.getPointLabel(i)
        << std::endl;
      labelMap.insert(std::pair<double, string>(x, _pb.getPointLabel(i)));
      x += _pb.getPathLength(i) / wholePathLength;
    }
    outfile << "# <label>     " << 1.0 << " " << setw(5) << _pb.getPointLabel(_pb.getPointSize()) << std::endl;
    labelMap.insert(std::pair<double, string>(1.0, _pb.getPointLabel(_pb.getPointSize())));
    outfile << "#" << std::endl;

    //writing high-symmetry qpoints [PN] //PN20180705
    stringstream ouths;  //PN20180705
    ouths << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);  //PN20180705
    ouths << setprecision(8);  //PN20180705

    // Write table of exact _qpoints + use label map to identify its labels
    x = 0.0;
    int subpath = 0;
    double xstep = 0.0;
    int p = 0;
    vector<double> exactPointPositions;
    for (uint i = 0; i < _qpoints.size(); i++) {
      // Check it
      if (isExactQPoint(_qpoints[i], _pc->getSuperCellStructure().lattice)) {
        // Is it new exact points
        uint j = 0;
        for (; j < exactPointPositions.size(); j++)
          if (exactPointPositions[j] == x) break;

        // If yes, add it....
        if (j == exactPointPositions.size()) {
          exactPointPositions.push_back(x);
          string name = "-";
          std::map<double, string>::iterator iter = labelMap.begin();
          for (; iter != labelMap.end(); iter++)
            if (aurostd::abs(iter->first - x) < _FLOAT_TOL_) break;
          if (iter != labelMap.end())
            name = iter->second;
          outfile << "# <exact>     " << x << " "
            << setw(5) << name
            << std::endl;
          ouths << "# <exact>     " << x << " "  //PN20180705
            << setw(5) << name  //PN20180705
            << setw(15) << _qpoints[i][1]<< setw(15) << _qpoints[i][2]<< setw(15) << _qpoints[i][3]  //PN20180705
            << '\n';  //PN20180705
        }
      }

      // Step of x will change in each subpath
      if (i % (_pb.getDensity() + 1) == 0) {
        if (i + 1 != _freqs.size())
          xstep = _pb.getPathLength(++subpath) / wholePathLength / (_pb.getDensity());
      }
      x += xstep;
      if (p != 0 && (p % _pb.getDensity()) == 0) {
        x -= xstep;
        p = -1;
      }
      p++;
    }
    outfile << "#" << std::endl;
    labelMap.clear();
    exactPointPositions.clear();

    // Write frequencies
    //[OBSOLETE PN20180705]path_segment.clear();  //[PN]
    //[OBSOLETE PN20180705]path.clear();          //[PN]
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    x = 0.0;
    subpath = 0;
    xstep = 0.0;
    p = 0;
    for (uint i = 0; i < _freqs.size(); i++) {
      ouths<<setw(15)<<_qpoints[i][1]<<setw(15) //PN20180705
        <<_qpoints[i][2]<<setw(15)<<_qpoints[i][3]<<setw(15)<<p<<setw(15)<<x<<"\n"; //PN20180705

      outfile << setw(4) << p << " ";
      //[OBSOLETE PN20180705]path_segment.push_back(p);  //[PN]
      outfile << setw(15) << x << " ";
      //[OBSOLETE PN20180705]path.push_back(x);  //[PN]
      for (uint j = 1; j <= _pc->getNumberOfBranches(); j++)
        outfile << setw(15) << _freqs[i](j) << " ";
      outfile << std::endl;

      // Step of x will change in each subpath
      if (i % (_pb.getDensity() + 1) == 0) {
        if (i + 1 != _freqs.size())
          xstep = _pb.getPathLength(++subpath) / wholePathLength / (_pb.getDensity());
      }
      x += xstep;
      if (p != 0 && (p % _pb.getDensity()) == 0) {
        x -= xstep;
        p = -1;
      }
      p++;
    }

    aurostd::stringstream2file(outfile, filename); //ME20181226
    if (!aurostd::FileExist(filename)) { //ME20181226
      message = "Cannot open output file " + filename + "."; //ME20181226
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    //PN //PN20180705
    string hskptsfile = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HSKPTS_FILE); //ME20181226
    aurostd::stringstream2file(ouths, hskptsfile); //ME20181226
    if (!aurostd::FileExist(hskptsfile)) { //ME20181226
      message = "Cannot open output file " + hskptsfile + "."; //ME20181226
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    //PN
  }

  //////////////////////////////////////////////////////////////////////////////

  bool PhononDispersionCalculator::isExactQPoint(const xvector<double>& qpoint,
      const xmatrix<double>& lattice) {
    xcomplex<double> iONE(0.0, 1.0);

    bool isExact = false;
    for (int i = 1; i <= 1; i++) {
      for (int j = 1; j <= 1; j++) {
        for (int k = 1; k <= 1; k++) {
          xvector<double> L = (((double)i) * lattice(1) +
              ((double)j) * lattice(2) +
              ((double)k) * lattice(3));
          xcomplex<double> p = exp(iONE * scalar_product(qpoint, L));
          if ((aurostd::abs(p.imag()) < _FLOAT_TOL_) &&
              (aurostd::abs(p.real() - 1.0) < _FLOAT_TOL_)) {
            isExact = true;
            break;
          }
        }
        if (isExact) break;
      }
      if (isExact) break;
    }

    return (isExact);
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190614 START
  // Write the eigenvalues into a VASP EIGENVAL-formatted file
  void PhononDispersionCalculator::writePHEIGENVAL(const string& directory) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHEIGENVAL_FILE);
    message = "Writing phonon eigenvalues into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, "APL", message, _pc->getDirectory(), *_pc->getOFStream(), *_pc->getOSS());
    stringstream eigenval;
    eigenval << createEIGENVAL();
    aurostd::stringstream2file(eigenval, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    // Also write PHKPOINTS and PHPOSCAR file
    writePHKPOINTS(directory);
    // OBSOLETE ME20191219 - PHPOSCAR is already written in KBIN::RunPhonons_APL
    // filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHPOSCAR_FILE);
    // xstructure xstr = _pc->getInputCellStructure();
    // xstr.is_vasp5_poscar_format = true;
    // stringstream poscar;
    // poscar << xstr;
    // aurostd::stringstream2file(poscar, filename);
    // if (!aurostd::FileExist(filename)) {
    //   string message = "Cannot open output file " + filename + ".";
    //   throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    // }
  }

  xEIGENVAL PhononDispersionCalculator::createEIGENVAL() {
    if (!_pc_set) {
      string message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    xEIGENVAL xeigen;
    stringstream outfile;
    // Header values
    xeigen.number_atoms = _pc->getInputCellStructure().atoms.size();
    xeigen.number_loops = 0;
    xeigen.spin = 0;
    xeigen.Vol = GetVolume(_pc->getInputCellStructure())/xeigen.number_atoms;
    xvector<double> lattice(3);
    lattice[1] = _pc->getInputCellStructure().a * 1E-10;
    lattice[2] = _pc->getInputCellStructure().b * 1E-10;
    lattice[3] = _pc->getInputCellStructure().c * 1E-10;
    xeigen.lattice = lattice;
    xeigen.POTIM = 0.5E-15;
    xeigen.temperature = _temperature;
    xeigen.carstring = "PHON";
    xeigen.title = _pc->_system;
    xeigen.number_electrons = 0;
    for (uint at = 0; at < xeigen.number_atoms; at++) {
      xeigen.number_electrons += _pc->getInputCellStructure().species_pp_ZVAL[at];
    }
    xeigen.number_kpoints = _freqs.size();
    xeigen.number_bands = _pc->getNumberOfBranches();

    // Data
    double weight = 1.0/_freqs.size();
    double factorTHz2Raw = _pc->getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2meV = _pc->getFrequencyConversionFactor(apl::RAW, apl::MEV);
    double conv = factorTHz2Raw * factorRaw2meV/1000;
    apl::PathBuilder::StoreEnumType store = _pb.getStore();
    xmatrix<double> c2f = inverse(trasp(ReciprocalLattice(_pc->getInputCellStructure().lattice)));

    xeigen.vweight.assign(xeigen.number_kpoints, weight);
    xeigen.vkpoint.resize(xeigen.number_kpoints, xvector<double>(3));
    xeigen.venergy.resize(xeigen.number_kpoints, deque<deque<double> >(xeigen.number_bands, deque<double>(1)));

    for (uint q = 0; q < xeigen.number_kpoints; q++) {
      if (store == apl::PathBuilder::CARTESIAN_LATTICE) xeigen.vkpoint[q] = c2f * _qpoints[q];
      else xeigen.vkpoint[q] = _qpoints[q];
      for (uint br = 0; br < xeigen.number_bands; br++) {
        xeigen.venergy[q][br][0] = conv * _freqs[q][br + 1];
      }
    }

    return xeigen;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // Write the k-point path into a VASP KPOINTS-formatted file
  void PhononDispersionCalculator::writePHKPOINTS(const string& directory) {
    string message = "";
    if (!_pc_set) {
      message = "PhononCalculator pointer not set.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_INIT_);
    }
    string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHKPOINTS_FILE);
    stringstream kpoints;
    kpoints << _pb.createKPOINTS(_pc->getSupercell());
    aurostd::stringstream2file(kpoints, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////
  //ME20190614 END

  //AS20201110 BEGIN
  /// Return q-path (as xKPOINTS) used to calculate phonon dispersion relations.
  xKPOINTS PhononDispersionCalculator::getPHKPOINTS()
  {
    return _pb.createKPOINTS(_pc->getSupercell());
  }
  //AS20201110 END

}  // namespace apl


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               PathBuilder                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  PathBuilder::PathBuilder() {
    free();
  }

  PathBuilder::PathBuilder(ModeEnumType mode) {
    free();
    setMode(mode);
  }

  PathBuilder::PathBuilder(const PathBuilder& that) {
    if (this != &that) free();
    copy(that);
  }

  PathBuilder& PathBuilder::operator=(const PathBuilder& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  void PathBuilder::copy(const PathBuilder& that) {
    if (this == &that) return;
    _mode = that._mode;
    _store = that._store;
    _path = that._path;
    _points = that._points;
    _labels = that._labels;
    reciprocalLattice = that.reciprocalLattice;
    cartesianLattice = that.cartesianLattice;
    _pointsVectorDimension = that._pointsVectorDimension;
    _pointsVectorStartingIndex = that._pointsVectorStartingIndex;
    _nPointsPerSubPath = that._nPointsPerSubPath;
  }

  PathBuilder::~PathBuilder() {
    free();
  }

  void PathBuilder::free() {
    _mode = SINGLE_POINT_MODE;
    _store = CARTESIAN_LATTICE;
    _path.clear();
    _points.clear();
    _labels.clear();
    reciprocalLattice.clear();
    cartesianLattice.clear();
    _pointsVectorDimension = 0;
    _pointsVectorStartingIndex = 0;
    _nPointsPerSubPath = 0;
  }

  void PathBuilder::clear() {
    free();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setMode(ModeEnumType mode) {
    _mode = mode;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setStore(StoreEnumType store) {
    _store = store;
  }

  // ///////////////////////////////////////////////////////////////////////////
  //ME20190614
  const apl::PathBuilder::StoreEnumType& PathBuilder::getStore() const {
    return _store;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::addPoint(const string& l, int dim, ...) {
    va_list arguments;
    xvector<double> point(dim,1);

    va_start(arguments, dim);
    for(int i = 0; i < dim; i++) {
      point(i+1) = va_arg(arguments, double);
    }
    va_end(arguments);
    //
    addPoint(l, point);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::addPoint(const string& l, const xvector<double>& p) {
    if( _points.empty() ) {
      _pointsVectorDimension = p.rows;
      _pointsVectorStartingIndex = p.lrows;
    }

    if( p.rows != _pointsVectorDimension ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::addPoint(); Wrong dimension of the point.");
      string message = "Wrong dimension of the point.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    if( p.lrows != _pointsVectorStartingIndex ) {
      xvector<double> pp(_pointsVectorDimension+_pointsVectorStartingIndex,_pointsVectorStartingIndex);
      for(int i = 0; i < _pointsVectorDimension; i++){
        pp(_pointsVectorStartingIndex+i) = p(p.lrows+i);
      }
      _points.push_back(pp);
    } else {
      _points.push_back(p);
    }

    _labels.push_back(l);
  }

  // ///////////////////////////////////////////////////////////////////////////

  int PathBuilder::getDensity() {
    return _nPointsPerSubPath;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setDensity(int n) {
    if( n < 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::setDensity(); The density should be >= 0.");
      string message = "The density should be >= 0.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    _nPointsPerSubPath = n;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::buildPath() {
    string message = "";
    // Remove the old path
    _path.clear();

    // Quick solution...
    if( _points.empty() ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::buildPath; There are no points.");
      message ="There are no points.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    };

    // Create path in the SINGLE_POINT_MODE
    if( _mode == SINGLE_POINT_MODE ) {
      for(uint i = 1; i < _points.size(); i++) {
        xvector<double> startPoint = _points[i-1];
        _path.push_back(startPoint);
        if( _nPointsPerSubPath != 0 ) {
          xvector<double> dPoint = ( _points[i] - _points[i-1] ) * ( 1.0 / (_nPointsPerSubPath) );
          for(int j = 1; j <= (int)_nPointsPerSubPath; j++)
          {
            xvector<double> p = startPoint + j*dPoint;
            _path.push_back(p);
          }
        }
      }
      _path.push_back(_points[_points.size()-1]);
    }

    // Create path in the COUPLE_POINT_MODE
    if( _mode == COUPLE_POINT_MODE ) {
      // If the number of points is odd -> escape
      if( _points.size() % 2 != 0 ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::buildPath(); The number of points is odd.");
        message = "The number of points is odd.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      xvector<double> startPoint(3), endPoint(3);
      endPoint = _points[0];
      for(uint i = 0; i < _points.size(); i+=2) {
        startPoint = _points[i];
        //if( aurostd::modulus( startPoint - endPoint ) > _FLOAT_TOL_ )
        //    _path.push_back(endPoint);
        _path.push_back(startPoint);
        endPoint = _points[i+1];

        if( _nPointsPerSubPath != 0 ) {
          xvector<double> dPoint = ( endPoint - startPoint ) * ( 1.0 / (_nPointsPerSubPath) );
          for(int j = 1; j <= (int)_nPointsPerSubPath; j++) {
            xvector<double> p = startPoint + j*dPoint;
            _path.push_back(p);
          }
        }
      }
      //_path.push_back(endPoint);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  uint PathBuilder::getPathSize() {
    return _path.size();
  }

  // ///////////////////////////////////////////////////////////////////////////

  uint PathBuilder::getPointSize() {
    if( _mode == SINGLE_POINT_MODE ) {
      return( _points.size() );
    }

    if( _mode == COUPLE_POINT_MODE ) {
      return( ( _points.size() / 2 ) + 1 );
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPointSize(); Unknown mode.");
    string message = "Unknown mode.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PathBuilder::getPathLength() {
    // Quick solution...
    if( _points.empty() ) return 0;

    //
    uint npaths;
    if( _mode == SINGLE_POINT_MODE )
      npaths = _points.size()-1;
    else if( _mode == COUPLE_POINT_MODE )
      npaths = _points.size()/2;
    else {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPointLength(); Unknown mode.");
      string message = "Unknown mode.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }

    //
    double length = 0.0;

    for(uint i = 1; i <= npaths; i++) {
      length += getPathLength(i);
    }

    return length;
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PathBuilder::getPathLength(uint i) {
    string message = "";
    if( i <= 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index. The index has to start from 1.");
      message = "Wrong index. The index has to start from 1.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }

    // Quick solution 1...
    if( _points.empty() ) return 0;

    double length;
    if( _mode == SINGLE_POINT_MODE ) {
      // Quick solution 2...
      if( i > _points.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index.");
        message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if( _store == RECIPROCAL_LATTICE ) {
        length = aurostd::modulus( F2C(trasp(reciprocalLattice),_points[i]) - F2C(trasp(reciprocalLattice),_points[i-1]) );
      } else {
        length = aurostd::modulus( _points[i] - _points[i-1] );
      }
      return length;
    }
    else if( _mode == COUPLE_POINT_MODE ) {
      // Quick solution 2...
      if( i > _points.size()/2 ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index.");
        message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if( _store == RECIPROCAL_LATTICE ) {
        length = aurostd::modulus( F2C(trasp(reciprocalLattice),_points[(i*2)-1]) - F2C(trasp(reciprocalLattice),_points[(i-1)*2]) );
      } else {
        length = aurostd::modulus( _points[(i*2)-1] - _points[(i-1)*2] );
      }
      return length;
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Unknown mode.");
    message = "Unknown mode.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  xvector<double> PathBuilder::getPoint(uint i) {
    string message = "";
    if( i <= 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index. The index has to start from 1.");
      message = "Wrong index. The index has to start from 1.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    if( _mode == SINGLE_POINT_MODE ) {
      if( i > _points.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index.");
        message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      return _points[i-1];
    } else if( _mode == COUPLE_POINT_MODE ) {
      if( i > (_points.size()/2)+1 ) {
        //ME20191031 - use xerro
        //throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index.");
        message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if( i == 1 ) return _points[0];
      if( i == (_points.size()/2)+1 ) return _points[_points.size()-1];
      return _points[2*i-3];
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPoint(); Unknown mode.");
    message = "Unknown mode.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  string PathBuilder::getPointLabel(uint i) {
    string message = "";
    if( i <= 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index. The index has to start from 1.");
      message = "Wrong index. The index has to start from 1.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }

    if( _mode == SINGLE_POINT_MODE ) {
      if( i > _labels.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index.");
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      return _labels[i-1];
    } else if( _mode == COUPLE_POINT_MODE ) {
      if( i > (_labels.size()/2)+1 ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index.");
        message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
      }
      if( i == 1 ) return _labels[0];
      if( i == (_labels.size()/2)+1 ) return _labels[_labels.size()-1];

      int h = 2*i-2;
      int l = 2*i-3;

      if( _labels[l] == _labels[h] )
        return _labels[l];
      else
        return(_labels[l]+"|"+_labels[h]);
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Unknown mode.");
    message = "Unknown mode.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  vector< aurostd::xvector<double> > PathBuilder::getPath(ModeEnumType mode, const string& userPath) {
    vector< xvector<double> > new_points;
    vector<string> new_labels;
    vector<string> tokens;

    // Generate users path by points
    aurostd::string2tokens(userPath,tokens,"|-");
    for(uint i = 0; i < tokens.size(); i++) {
      uint j = 0;
      for( ; j < _labels.size(); j++) {
        // This works fine, but notorious changes of the Gamma point label
        // in aflow_kpoints.cpp lead me to do it like: G, Gamma, of \Gamma or
        // or anything like GAMMA \G or GaMmA or I do not know what else is
        // still a gamma point labeled as G!
        //if( _labels[j] == tokens.at(i) ) break;
        if( _labels[j].find(tokens.at(i)) != string::npos ) break;
      }
      if( j == _labels.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPath(); Undefined label of the point.");
        string message = "Undefined label of the point.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      new_points.push_back(_points[j]);      
      //new_labels.push_back(_labels[j]);
      // I will use my labels...
      new_labels.push_back(tokens.at(i));
    }

    // update
    _points.clear();
    _labels.clear();
    for(uint i = 0; i < new_points.size(); i++) {
      _points.push_back(new_points[i]);
      _labels.push_back(new_labels[i]);
    }
    new_points.clear();
    new_labels.clear();

    // Build full path
    _mode = mode;
    buildPath();

    //
    return _path;
  }

  // ///////////////////////////////////////////////////////////////////////////

  std::vector< aurostd::xvector<double> > PathBuilder::getPath() {
    buildPath();
    return _path;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::transform(const aurostd::xmatrix<double>& m) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG){
      cerr << "PathBuilder::transform():" << " matrix:" << std::endl;
      cerr << m;
    }
    for(uint i = 0; i < _points.size(); i++) {
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " IN: " << _labels[i] << ": " << _points[i] << std::endl;}
      _points[i] = m * _points[i];
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " OUT: " << _labels[i] << ": " << _points[i] << std::endl;}
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::pointsAreDefinedFor(const xstructure& primitiveStructure, StoreEnumType store) {
    // Transform from the reciprocal lattice of the primitive cell to cartesian coordinates
    if( store == RECIPROCAL_LATTICE ) {
      transform( trasp(ReciprocalLattice(primitiveStructure.lattice)) );
    }

    //
    _store = store;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //void PathBuilder::transformPointsFor(const xstructure& supercellStructure, StoreEnumType store) {
  //  if( store == CARTESIAN_LATTICE ) {
  //    if( _store == RECIPROCAL_LATTICE )
  //      transform( trasp(ReciprocalLattice(supercellStructure.lattice)) );
  //    else
  //  }
  //  else if( store == RECIPROCAL_LATTICE ) {
  //    transform( inverse(trasp(ReciprocalLattice(supercellStructure.lattice))) );
  //  }

  //  _store = store;
  //  cartesianLattice = supercellStructure.lattice;
  //  reciprocalLattice = ReciprocalLattice(supercellStructure.lattice);

  //}

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::defineCustomPoints(const string& coords,const string& labels,const Supercell& sc,bool CARTESIAN_COORDS){ //CO20180409
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string message = "";
    vector<string> vcoords,vlabels;
    vector<double> coordinate;
    aurostd::string2tokens(coords,vcoords,";");
    aurostd::string2tokens(labels,vlabels,",;");

    if(vcoords.empty()){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); No input coordinates found");
      message = "No input coordinates found";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }
    if(vlabels.empty()){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); No input labels found");
      message = "No input labels found";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }
    if(vcoords.size()!=vlabels.size()){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); Size of input coordinates does not match size of input labels");
      message = "Size of input coordinates does not match size of input labels";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }
    if(vcoords.size()<2){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); Path size should include at least two points");
      message = "Path size should include at least two points";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
    }

    for(uint i=0;i<vcoords.size();i++){
      aurostd::string2tokens<double>(vcoords[i],coordinate,",");
      if(coordinate.size()!=3){
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); Input coordinate["+aurostd::utype2string(i)+"] is not dimension==3");
        message = "Input coordinate["+aurostd::utype2string(i)+"] is not dimension==3";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " found point " << vlabels[i] << ": (" << (CARTESIAN_COORDS?"cartesian":"fractional") << ") " << coordinate[0] << "," << coordinate[1] << "," << coordinate[2] << std::endl;}
      addPoint(vlabels[i],3,coordinate[0],coordinate[1],coordinate[2]);
    }

    if(!CARTESIAN_COORDS){
      //ME20200203 - these are custom points, i.e. they are user generated based
      // on a user-defined input structure. AFLOW should not switch the basis
      // behind the scenes. Assume that the user knows what they are doing.
      //transform( trasp(ReciprocalLattice(sc.getPrimitiveStructure().lattice)) ); //must be reciprocal
      transform( trasp(ReciprocalLattice(sc.getInputStructure().lattice)) );
    }

    //set lattices
    cartesianLattice = sc.getSupercellStructure().lattice;
    reciprocalLattice = ReciprocalLattice(sc.getSupercellStructure().lattice);

    setMode(SINGLE_POINT_MODE); //unless modified later
    _store = CARTESIAN_LATTICE;
  }

  void PathBuilder::takeAflowElectronicPath(const string& latticeName,const Supercell& sc) //CO20180409
    //const xstructure& inStructure,
    //const xstructure& scStructure) 
  {
    // Get path from electronic structure...
    stringstream fileKPOINTS;
    // input is not the reciprocal lattice!
    bool foundBZ;
    fileKPOINTS << LATTICE::KPOINTS_Directions(latticeName,sc.getInputStructure().lattice,10,sc.getInputStructure().iomode,foundBZ);
    fileKPOINTS.flush();

    if( !foundBZ ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::takeAflowElectronicPath(); The BZ not found for this lattice.");
      string message = "The Brillouin zone was not found for this lattice.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    //   cerr << fileKPOINTS.str() << std::endl;
    string line;
    vector<string> tokens;
    int nLine = 0;
    while( getline( fileKPOINTS, line ) ) {
      if( nLine++ < 4 ) continue;
      if( line.empty() ) continue;
      if( line.size() == 1 ) continue;
      //cout << line << std::endl;
      aurostd::string2tokens(line,tokens,string(" !"));
      addPoint(tokens[3],3,atof(tokens[0].c_str()),atof(tokens[1].c_str()),atof(tokens[2].c_str()));
      tokens.clear();
      line.clear();
    }
    fileKPOINTS.clear();

    //CO20180406 - with help from Mike Mehl
    //these coordinates are wrt PRIMITIVE lattice, so if input is not PRIMITIVE, the conversion will not work
    // Transform from the reciprocal lattice of the primitive cell to
    // cartesian coordinates (note; in klattice vectors are rows! )
    transform( trasp(ReciprocalLattice(sc.getPrimitiveStructure().lattice)) );  //inStructure.lattice)) );  //CO20180406

    //set lattices
    cartesianLattice = sc.getSupercellStructure().lattice;
    reciprocalLattice = ReciprocalLattice(sc.getSupercellStructure().lattice);

    // set mode
    setMode(COUPLE_POINT_MODE);
    _store = CARTESIAN_LATTICE;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190614
  // Writes the k-points path into a VASP KPOINTS-formatted file
  xKPOINTS PathBuilder::createKPOINTS(const Supercell& sc) {
    vector<xvector<double> > points(2 * (_points.size() - 1), xvector<double>(3));
    vector<string> labels(2 * (_points.size() - 1));
    if (_mode == SINGLE_POINT_MODE) {
      points[0] = _points[0];
      labels[0] = _labels[0];
      for (uint i = 1; i < _points.size() - 1; i++) {
        points[2 * i - 1] = _points[i];
        points[2 * i] = _points[i];
        labels[2 * i - 1] = _labels[i];
        labels[2 * i] = _labels[i];
      }
      points.back() = _points.back();
      labels.back() = _labels.back();
    } else {
      points = _points;
      labels = _labels;
    }
    //ME20200117 - Convert to reciprocal coordinates of the
    // original structure or the distances will be wrong
    if (_store == CARTESIAN_LATTICE) {
      xmatrix<double> c2f = inverse(trasp(ReciprocalLattice(sc.getOriginalStructure().lattice)));
      for (uint i = 0; i < points.size(); i++) points[i] = c2f * points[i];
    } else {
      xmatrix<double> f2c = trasp(ReciprocalLattice(sc.getInputStructure().lattice));
      xmatrix<double> c2f_orig = inverse(trasp(ReciprocalLattice(sc.getOriginalStructure().lattice)));
      for (uint i = 0; i < points.size(); i++) points[i] = c2f_orig * f2c * points[i];
    }

    xKPOINTS xkpts;
    xkpts.is_KPOINTS_PATH = true;
    xkpts.is_KPOINTS_NNN = false;
    xkpts.vpath = labels;
    xkpts.vkpoints = points;
    xkpts.title = xkpts.createStandardTitlePath(sc.getPrimitiveStructure());
    xkpts.path_grid = _path.size()/(_points.size()/2);
    xkpts.grid_type = "Line-mode";
    xkpts.path_mode = "reciprocal";

    return xkpts;
  }

} // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
