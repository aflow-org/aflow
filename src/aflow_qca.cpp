//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_QCA_CPP_
#define _AFLOW_QCA_CPP_

#include "aflow_qca.h"
using namespace std::placeholders;

// ###############################################################################
//            AFLOW Quasi-Chemical Approximation (QCA) (2022-)
// ###############################################################################

namespace qca {
  /// @class QuasiChemApproxCalculator
  /// 
  /// @brief calculator to perform the quasi-chemical approximation
  ///
  /// @authors
  /// @mod{SD,20220204,created class}
  ///
  /// @see
  /// @doi{10.1016/j.actamat.2018.07.042}
  QuasiChemApproxCalculator::QuasiChemApproxCalculator(ostream& oss) : xStream(oss), m_initialized(false) {initialize();}
  QuasiChemApproxCalculator::QuasiChemApproxCalculator(ofstream& FileMESSAGE, ostream& oss) : xStream(FileMESSAGE, oss), m_initialized(false) {initialize();}
  QuasiChemApproxCalculator::QuasiChemApproxCalculator(const aurostd::xoption& qca_flags, ostream& oss) : xStream(oss), m_initialized(false) {initialize(qca_flags);}
  QuasiChemApproxCalculator::QuasiChemApproxCalculator(const aurostd::xoption& qca_flags, ofstream& FileMESSAGE, ostream& oss) : xStream(FileMESSAGE, oss), m_initialized(false) {initialize(qca_flags);}
  QuasiChemApproxCalculator::QuasiChemApproxCalculator(const QuasiChemApproxCalculator& b) : xStream(*b.getOFStream(), *b.getOSS()) {copy(b);}

  QuasiChemApproxCalculator::~QuasiChemApproxCalculator() {xStream::free(); free();}

  const QuasiChemApproxCalculator& QuasiChemApproxCalculator::operator=(const QuasiChemApproxCalculator& b) {
    if (this != &b) {copy(b);}
    return *this;
  }
  void QuasiChemApproxCalculator::clear() {QuasiChemApproxCalculator a; copy(a);}

  void QuasiChemApproxCalculator::free() {
    m_initialized = false;
    m_aflags.clear();
    min_sleep = 0;
    print = "";
    screen_only = false;
    image_only = false;
    calc_binodal = false;
    use_sg = false;
    rootdirpath = "";
    aflowlibpath = "";
    plattice = "";
    elements.clear();
    aflow_max_num_atoms = 0;
    max_num_atoms = 0;
    cv_cut = 0.0;
    conc_npts = 0;
    conc_curve = false;
    conc_curve_range.clear();
    conc_macro.clear();
    temp_npts = 0;
    temp_range.clear();
    temp.clear();
    beta.clear();
    alloyname = "";
    rundirpath = "";
    vstr_aflow.clear();
    lat_atat = "";
    vstr_ce.clear();
    mapstr.clear();
    skipstr.clear();
    nelem = 0;
    ncluster = 0;
    nconc = 0;
    cv_cluster = 0.0;
    num_atom_cluster.clear();
    degeneracy_cluster.clear();
    conc_cluster.clear();
    excess_energy_cluster.clear();
    prob_ideal_cluster.clear();
    soln0.clear();
    prob_cluster.clear();
    param_ec.first = 0.0; param_ec.second = 0.0;
    rel_s.clear();
    binodal_curve.clear();
  }

  void QuasiChemApproxCalculator::copy(const QuasiChemApproxCalculator& b) {
    xStream::copy(b);
    m_initialized = b.m_initialized;
    m_aflags = b.m_aflags;
    min_sleep = b.min_sleep;
    print = b.print;
    screen_only = b.screen_only;
    image_only = b.image_only;
    calc_binodal = b.calc_binodal;
    use_sg = b.use_sg;
    rootdirpath = b.rootdirpath;
    aflowlibpath = b.aflowlibpath;
    plattice = b.plattice;
    elements = b.elements;
    aflow_max_num_atoms = b.aflow_max_num_atoms;
    max_num_atoms = b.max_num_atoms;
    cv_cut = b.cv_cut;
    conc_npts = b.conc_npts;
    conc_curve = b.conc_curve;
    conc_curve_range = b.conc_curve_range;
    conc_macro = b.conc_macro;
    temp_npts = b.temp_npts;
    temp_range = b.temp_range;
    temp = b.temp;
    beta = b.beta;
    alloyname = b.alloyname;
    rundirpath = b.rundirpath;
    vstr_aflow = b.vstr_aflow;
    lat_atat = b.lat_atat;
    vstr_ce = b.vstr_ce;
    mapstr = b.mapstr;
    skipstr = b.skipstr;
    nelem = b.nelem;
    ncluster = b.ncluster;
    nconc = b.nconc;
    cv_cluster = b.cv_cluster;
    num_atom_cluster = b.num_atom_cluster;
    degeneracy_cluster = b.degeneracy_cluster;
    conc_cluster = b.conc_cluster;
    excess_energy_cluster = b.excess_energy_cluster;
    prob_ideal_cluster = b.prob_ideal_cluster;
    soln0 = b.soln0;
    prob_cluster = b.prob_cluster;
    param_ec = b.param_ec;
    rel_s = b.rel_s;
    binodal_curve = b.binodal_curve;
  }
  
  bool QuasiChemApproxCalculator::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }
  bool QuasiChemApproxCalculator::initialize(ofstream& FileMESSAGE, ostream& oss) {
    xStream::initialize(FileMESSAGE, oss);
    return initialize();
  }
  bool QuasiChemApproxCalculator::initialize() {
    free();
    return m_initialized;
  }
  bool QuasiChemApproxCalculator::initialize(const aurostd::xoption& qca_flags, ostream& oss) {
    xStream::initialize(oss);
    return initialize(qca_flags);
  }
  bool QuasiChemApproxCalculator::initialize(const aurostd::xoption& qca_flags, ofstream& FileMESSAGE, ostream& oss) {
    xStream::initialize(FileMESSAGE, oss);
    return initialize(qca_flags);
  }
  /// @brief initializes the QCA variables
  ///
  /// @param qca_flags QCA xoptions object
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  bool QuasiChemApproxCalculator::initialize(const aurostd::xoption& qca_flags) {
    free();
    min_sleep = DEFAULT_QCA_MIN_SLEEP_SECONDS;
    rootdirpath = aurostd::getPWD();
    aflow_max_num_atoms = DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS;
    max_num_atoms = DEFAULT_QCA_MAX_NUM_ATOMS;
    cv_cut = DEFAULT_QCA_CV_CUTOFF;
    conc_npts = DEFAULT_QCA_CONC_NPTS;
    temp_range = {DEFAULT_QCA_TEMP_MIN, DEFAULT_QCA_TEMP_MAX};
    temp_npts = DEFAULT_QCA_TEMP_NPTS;
    print = DEFAULT_QCA_PRINT;
    readQCAFlags(qca_flags);
    m_initialized = true;
    return m_initialized;
  }

  /// @brief reads the QCA flags
  ///
  /// @param qca_flags QCA xoptions object
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::readQCAFlags(const aurostd::xoption& qca_flags) {
    if (qca_flags.getattachedscheme("QCA::PLATTICE").empty() || qca_flags.getattachedscheme("QCA::ELEMENTS").empty()) {
      string message = "QCA::PLATTICE and QCA::ELEMENTS flags must be defined";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
    }
    plattice = qca_flags.getattachedscheme("QCA::PLATTICE");
    aurostd::string2tokens(qca_flags.getattachedscheme("QCA::ELEMENTS"), elements, ",");
    if (!qca_flags.getattachedscheme("QCA::DIRECTORY").empty()) {rootdirpath = qca_flags.getattachedscheme("QCA::DIRECTORY");}
    if (!qca_flags.getattachedscheme("QCA::AFLOWLIB_DIRECTORY").empty()) {aflowlibpath = qca_flags.getattachedscheme("QCA::AFLOWLIB_DIRECTORY");}
    if (!qca_flags.getattachedscheme("QCA::AFLOW_MAX_NUM_ATOMS").empty()) {aflow_max_num_atoms = aurostd::string2utype<int>(qca_flags.getattachedscheme("QCA::AFLOW_MAX_NUM_ATOMS"));}
    if (!qca_flags.getattachedscheme("QCA::CV_CUTOFF").empty()) {cv_cut = aurostd::string2utype<double>(qca_flags.getattachedscheme("QCA::CV_CUTOFF"));}
    if (!qca_flags.getattachedscheme("QCA::CONC_NPTS").empty()) {conc_npts = aurostd::string2utype<int>(qca_flags.getattachedscheme("QCA::CONC_NPTS"));}
    if (!qca_flags.getattachedscheme("QCA::TEMP_RANGE").empty()) {aurostd::string2tokens(qca_flags.getattachedscheme("QCA::TEMP_RANGE"), temp_range, ",");}
    if (!qca_flags.getattachedscheme("QCA::TEMP_NPTS").empty()) {temp_npts = aurostd::string2utype<int>(qca_flags.getattachedscheme("QCA::TEMP_NPTS"));}
    if (!qca_flags.getattachedscheme("QCA::PRINT").empty()) {print = qca_flags.getattachedscheme("QCA::PRINT");}
    if (qca_flags.flag("QCA::SCREEN_ONLY")) {screen_only = true;}
    if (qca_flags.flag("QCA::IMAGE_ONLY")) {image_only = true;}
    if (qca_flags.flag("QCA::BINODAL")) {calc_binodal = true;}
    if (qca_flags.flag("QCA::USE_SG")) {use_sg = true;}
  }

  /// @brief reads the QCA parameters
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::printParams() {
    stringstream message;
    message << "Root directory = " << rootdirpath;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "AFLOWLIB directory = " << aflowlibpath;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Maximum number of atoms in AFLOW = " << aflow_max_num_atoms;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Parent lattice = " << plattice;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Alloy elements = " << aurostd::joinWDelimiter(elements, ",");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Maximum number of atoms in cluster expansion = " << max_num_atoms;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Cross validation cutoff = " << cv_cut;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Concentration curve range = " << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < conc_curve_range.size(); i++) {message << conc_curve_range[i] << " ";}
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Number of points (concentration) = " << conc_npts;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Temperature range (K) = " <<  std::fixed << std::setprecision(2) << temp_range[0] << " " << temp_range[1];
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Number of points (temperature) = " << temp_npts;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Output format = " << print;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  /// @brief checks and fixes errors in the input variables
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::errorFix() {
    // Clean-up names and remove duplicates
    rootdirpath = aurostd::CleanFileName(rootdirpath);
    plattice = aurostd::tolower(plattice);
    print = aurostd::tolower(print);
    aurostd::sort_remove_duplicates(elements);
    nelem = elements.size();
    // Check if min sleep is at least 1 sec
    if (min_sleep < 1) {min_sleep = 1;}
    // Check if max number of atoms in cluster is valid
    if (aflow_max_num_atoms < 1 || max_num_atoms < 1) {
      string message = "Maximum number of atoms per cluster must be at least 1";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if directory is writable
    if (!aurostd::DirectoryMake(rootdirpath)) {
      string message = "Cannot create directory";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    // Check if parent lattice is valid
    if (plattice != "fcc" && plattice != "bcc" && plattice != "hcp") {
      string message = "Parent lattice \"" + plattice + "\" is invalid";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    // Check if alloy is at least binary
    if (nelem < 2) {
      string message = "Alloy must be at least binary";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if HCP parent lattice is used for greater than binary
    if (plattice == "hcp" && nelem != 2) {
      string message = "HCP parent lattice of alloys greater than binary not supported";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    // Check if elements are valid, construct alloy name
    for (uint i = 0; i < nelem; i++) {
      if (!xelement::xelement::isElement(elements[i])) {
        string message = "Element \"" + elements[i] + "\" is invalid";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      alloyname += elements[i];
    }
    rundirpath += rootdirpath + "/" + pflow::arity_string(nelem, false, false) + "/" + plattice + "/" + alloyname;
    // Check if concentration range has enough points
    if (conc_npts < 2) {
      string message = "Number of points for the concentration range must be at least 2";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if concentration curve format is valid
    if (conc_curve) {
      if (conc_curve_range.size() != 2 * nelem) {
        string message = "Concentration curve must have format [X1_start, X1_end, X2_start, X2_end,...X(K)_start, X(K)_end]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      // Check if concentration curve is within [0,1] and sums to 1
      double totconc_init = 0.0, totconc_final = 0.0;
      vector<double> vconc_init, vconc_final;
      vector<size_t> vzeros_init, vzeros_final;
      for (size_t i = 0; i < nelem && totconc_init <= 1.0 && totconc_final <= 1.0; i++) {
        if (conc_curve_range[2 * i] < 0 || conc_curve_range[2 * i] > 1 ||
            conc_curve_range[2 * i + 1] < 0 || conc_curve_range[2 * i + 1] > 1) {
          string message = "Concentration curve must be defined on [0,1]";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
        totconc_init += conc_curve_range[2 * i]; totconc_final += conc_curve_range[2 * i + 1];
        vconc_init.push_back(conc_curve_range[2 * i]); vconc_final.push_back(conc_curve_range[2 * i + 1]);
      }
      if (!aurostd::isequal(totconc_init, 1.0) || !aurostd::isequal(totconc_final, 1.0)) {
        string message = "Total concentration must sum to 1";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
      // Redefine concentration curve from [0,1] to (0,1)
      aurostd::WithinList(vconc_init, 0.0, vzeros_init); aurostd::WithinList(vconc_final, 0.0, vzeros_final);
      double cdelta_init = CONC_SHIFT * (double)vzeros_init.size() / ((double)vconc_init.size() - (double)vzeros_init.size());
      double cdelta_final = CONC_SHIFT * (double)vzeros_final.size() / ((double)vconc_final.size() - (double)vzeros_final.size());
      for (size_t i = 0; i < conc_curve_range.size(); i++) {
        if (aurostd::isequal(conc_curve_range[i], 0.0)) {
          conc_curve_range[i] += CONC_SHIFT;
        }
        else {
          conc_curve_range[i] -= (i % 2) ? cdelta_final : cdelta_init;
        }
      }
    }
    // Check if temperature range format is valid
    if (temp_range.size() != 2) {
      string message = "Temperature range must have format [T_start T_end]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    // Check if temperature range has enough points
    if (temp_npts < 2) {
      string message = "Number of points for the temperature range must be at least 2";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if temperature values are valid
    if (temp_range[0] < 0 || temp_range[1] < 0) {
      string message = "Temperature cannot be below 0K";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if output format is valid
    if (print != "txt" && print != "json" &&
        print != "pdf" && print != "eps" && print != "png") {
      string message = "Format \"" + print + "\" is invalid";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
  }

  /// @brief calculates the binodal curve and the quantities needed for the calculation
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  ///
  /// @see
  /// @doi{10.1016/j.actamat.2018.07.042}
  void QuasiChemApproxCalculator::calculateBinodal() {
    aurostd::DirectoryMake(rundirpath);
    if (aurostd::FileExist(rundirpath + "/fit.out") && aurostd::FileExist(rundirpath + "/predstr.out")) { // read ATAT data
      vstr_ce = getATATXstructures(0, true);
    }
    else { // run ATAT
      readAFLOWXstructures();
      lat_atat = getLatForATAT(); // unscaled
      vstr_ce = getATATXstructures(max_num_atoms);
      calculateMapForXstructures(getATATXstructures(aflow_max_num_atoms), vstr_aflow); // map ATAT xstrs to AFLOW xstrs because ATAT cannot identify AFLOW xstrs
      lat_atat = getLatForATAT(true); // scaled
      generateFilesForATAT();
      runATAT();
    }
    readCVCluster();
    calculateNumAtomCluster();
    calculateConcentrationCluster();
    readExcessEnergyCluster();
    setCongruentClusters();
    calculateDegeneracyCluster();
    calculateConcentrationMacro();
    calculateTemperatureRange();
    calculateRelativeEntropyEC();
    if (calc_binodal) {
      calculateProbabilityIdealCluster();
      checkProbabilityIdeal();
      calculateProbabilityCluster();
      try {
        checkProbabilityEquilibrium();
      }
      catch (aurostd::xerror& err) {
        stringstream message;
        message << "Equilibrium probabilities do not satisfy the concentration constraints, skipping binodal curve calculation";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        return;
      }
      calculateRelativeEntropy();
      calculateBinodalCurve();
    }
    return;
  }

  /// @brief reads the AFLOW xstructures from the AFLOW database
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::readAFLOWXstructures() {
    stringstream message;
    vector<string> vstr_aflowurl, elements_xc;
    string alloyname_PBE = "", aflowlib = "", aflowurl = "";
    stringstream oss;
    aurostd::xcombos xc;
    vector<int> index_xc;
    for (uint ialloy = 2; ialloy <= 3; ialloy++) { // AFLOW database only contains data for binaries and ternaries
      xc = aurostd::xcombos(elements.size(), ialloy, 'C');
      while (xc.increment()) {
        elements_xc.clear();
        index_xc = xc.getCombo();
        for (size_t i = 0; i < index_xc.size(); i++) {
          if (index_xc[i] == 1) {elements_xc.push_back(elements[i]);}
        }
        aurostd::sort(elements_xc);
        alloyname_PBE = "";
        for (size_t i = 0; i < elements_xc.size(); i++) {alloyname_PBE += AVASP_Get_PseudoPotential_PAW_PBE(elements_xc[i]);}
        aflowurl = "aflowlib.duke.edu:AFLOWDATA/LIB" + aurostd::utype2string<uint>(ialloy) + "_RAW/" + alloyname_PBE;
        if (plattice == "fcc") {
          if (ialloy == 2) {
            for (uint i = 1; i < 30; i++) {
              vstr_aflowurl.push_back(aflowurl + "/" + aurostd::utype2string<uint>(i));
            }
          }
          if (ialloy == 3) {
            vstr_aflowurl.push_back(aflowurl + "/" + "TFCC001.ABC");
            vstr_aflowurl.push_back(aflowurl + "/" + "TFCC002.ABC");
            vstr_aflowurl.push_back(aflowurl + "/" + "TFCC003.ABC");
            for (uint i = 4; i < 17; i++) {
              vstr_aflowurl.push_back(aflowurl + "/" + "TFCC00" + aurostd::utype2string<uint>(i) + ".ABC");
              vstr_aflowurl.push_back(aflowurl + "/" + "TFCC00" + aurostd::utype2string<uint>(i) + ".BCA");
              vstr_aflowurl.push_back(aflowurl + "/" + "TFCC00" + aurostd::utype2string<uint>(i) + ".CAB");
            }
          }
        }
        else if (plattice == "bcc") {
          if (ialloy == 2) {
            for (uint i = 58; i < 87; i++) {
              vstr_aflowurl.push_back(aflowurl + "/" + aurostd::utype2string<uint>(i));
            }
          }
          if (ialloy == 3) {
            vstr_aflowurl.push_back(aflowurl + "/" + "TBCC001.ABC");
            vstr_aflowurl.push_back(aflowurl + "/" + "TBCC002.ABC");
            vstr_aflowurl.push_back(aflowurl + "/" + "TBCC003.ABC");
            for (uint i = 4; i < 17; i++) {
              vstr_aflowurl.push_back(aflowurl + "/" + "TBCC00" + aurostd::utype2string<uint>(i) + ".ABC");
              vstr_aflowurl.push_back(aflowurl + "/" + "TBCC00" + aurostd::utype2string<uint>(i) + ".BCA");
              vstr_aflowurl.push_back(aflowurl + "/" + "TBCC00" + aurostd::utype2string<uint>(i) + ".CAB");
            }
          }
        }
        else if (plattice == "hcp") {
          if (ialloy == 2) {
            for (uint i = 115; i < 178; i++) {
              vstr_aflowurl.push_back(aflowurl + "/" + aurostd::utype2string<uint>(i));
            }
          }
        }
      }
    }
    message << "Reading AFLOW xstructures from AFLOW database";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    aflowlib::EntryLoader el;
    el.m_xstructure_original = true;
    el.m_xstructure_relaxed = true;
    el.loadAURL(vstr_aflowurl);
    if (!aflowlibpath.empty()) {el.loadFolders(aflowlibpath);}
    size_t ie = 0;
    for (std::shared_ptr<aflowlib::_aflowlib_entry> entry : *el.m_entries_flat) {
      ie++;
      entry->vstr[0].qm_E_cell = entry->enthalpy_cell; // ATAT needs energy per cell
      if (use_sg && (entry->spacegroup_orig == entry->spacegroup_relax)) {
        vstr_aflow.push_back(entry->vstr[0]);
      }
      else {
        vstr_aflow.push_back(entry->vstr[0]);
      }
      pflow::updateProgressBar(ie, (*el.m_entries_flat).size(), *p_oss);
    }
  }

  /// @brief calculates the map between two groups of xstructures
  ///
  /// @param vstr1 first group of xstructures
  /// @param vstr2 second group of xstructures
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    message << "Mapping xstructures between AFLOW and ATAT";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<xstructure> vstr3 = vstr1;
    vstr3.insert(vstr3.end(), vstr2.begin(), vstr2.end());
    XtalFinderCalculator xtal_calc;
    vector<vector<uint>> gsx = xtal_calc.groupSimilarXstructures(vstr3);
    for (size_t i = 0; i < gsx.size(); i++) {
      std::sort(gsx[i].begin(), gsx[i].end());
      if (gsx[i][0] > vstr1.size() -  1) { // index is not part of vstr1
        continue;
      }
      else if (gsx[i].size() == 1) { // no match for xstr in vstr1
        mapstr.push_back(-1);
      }
      else { // take first match of xstr in vstr2 to xstr in vstr1
        for (size_t j = 1; j < gsx[i].size(); j++) {
          if (gsx[i][j] > vstr1.size() - 1) {
            mapstr.push_back((int)gsx[i][j] - (int)vstr1.size());
            break;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " i=" << i << " | ";
        for (size_t j = 0; j < gsx[i].size(); j++) {cerr << gsx[i][j] << " ";}
        cerr << endl;
      }
    }
    if (mapstr.size() != vstr1.size()) {
      message << "Something went wrong in the mapping of the xstructures; vstr1.size()=" << vstr1.size() << " | mapstr.size()=" << mapstr.size() << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
  }

  /// @brief gets the lattice file for ATAT
  ///
  /// @param scale scale the lattice vectors by the lattice constant
  ///
  /// @return lattice string for ATAT
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  string QuasiChemApproxCalculator::getLatForATAT(bool scale) const {
    stringstream oss;
    double alat = 0.0;
    for (size_t i = 0; i < nelem; i++) {alat += GetAtomRadiusCovalent(elements[i]);}
    alat /= nelem;
    xmatrix<double> lattice(3, 3);
    xvector<double> angles = 90.0 * aurostd::ones_xv<double>(3);
    xvector<double> coorsys = aurostd::ones_xv<double>(3);
    oss.precision(_DOUBLE_WRITE_PRECISION_MAX_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    if (plattice == "fcc") {
      alat *= 2.0 * std::sqrt(2.0);
      lattice(1, 1) = 0.0; lattice(2, 1) = 0.5; lattice(3, 1) = 0.5;
      lattice(1, 2) = 0.5; lattice(2, 2) = 0.0; lattice(3, 2) = 0.5;
      lattice(1, 3) = 0.5; lattice(2, 3) = 0.5; lattice(3, 3) = 0.0;
    }
    else if (plattice == "bcc") {
      alat *= 4.0 / std::sqrt(3.0);
      lattice(1, 1) = -0.5; lattice(2, 1) = 0.5; lattice(3, 1) = 0.5;
      lattice(1, 2) = 0.5; lattice(2, 2) = -0.5; lattice(3, 2) = 0.5;
      lattice(1, 3) = 0.5; lattice(2, 3) = 0.5; lattice(3, 3) = -0.5;
    }
    else if (plattice == "hcp") {
      alat *= 2.0;
      lattice = aurostd::eye<double>(3, 3);
      coorsys(3) = std::sqrt(8.0 / 3.0);
      angles(3) = 120.0;
    }
    alat = (scale) ? aurostd::round(alat, 3) : 1.0;
    for (uint i = 1; i <= 3; i++) {oss << alat * coorsys(i) << " ";}
    for (uint i = 1; i <= 3; i++) {oss << angles(i) << " ";}
    oss << endl;
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        oss << lattice(i, j) << " ";
      }
      oss << endl;
    }
    oss << 0.0 << " " << 0.0 << " " << 0.0 << " ";
    for (size_t i = 0; i < nelem; i++) {
      oss << elements[i] << ",";
    }
    if (plattice == "hcp") {
      oss << endl << 2.0 / 3.0 << " " << 1.0 / 3.0 << " " << 0.5 << " ";
      for (size_t i = 0; i < nelem; i++) {
        oss << elements[i] << ",";
      }
    }
    oss << endl;
    return oss.str();
  }
  
  /// @brief gets the ATAT xstructures
  ///
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  /// @param fromfile read ATAT xstructures from file
  ///
  /// @return xstructures from ATAT runs
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  vector<xstructure> QuasiChemApproxCalculator::getATATXstructures(const int max_num_atoms, bool fromfile) const {
    stringstream message;
    vector<xstructure> vstr;
    stringstream oss;
    if (fromfile) {
      message << "Reading ATAT xstructures from directory = " << rundirpath;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      vector<string> files;
      vector<xstructure> vstr_tmp;
      vector<uint> index;
      aurostd::DirectoryLS(rundirpath, files);
      for (size_t i = 0; i < files.size(); i++) {
        if (aurostd::FileExist(rundirpath + "/" + files[i] + "/str.out")) {
          aurostd::file2stringstream(rundirpath + "/" + files[i] + "/str.out", oss);
          vstr_tmp.push_back(xstructure(oss, IOATAT_STR));
          index.push_back(aurostd::string2utype<uint>(files[i]));
          aurostd::StringstreamClean(oss);
        }
        pflow::updateProgressBar(i, files.size() - 1, *p_oss);
      }
      vstr.resize(index.size());
      for (size_t i = 0; i < index.size(); i++) {vstr[index[i]] = vstr_tmp[i];}
      return vstr;
    }
    double scale = 0.0;
    for (size_t i = 0; i < nelem; i++) {scale += GetAtomRadiusCovalent(elements[i]);}
    scale /= nelem;
    if (plattice == "fcc") {
      scale *= 2.0 * std::sqrt(2.0);
    }
    else if (plattice == "bcc") {
      scale *= 4.0 / std::sqrt(3.0);
    }
    else if (plattice == "hcp") {
      scale *= 2.0;
    }
    scale = aurostd::round(scale, 3);
    vector<string> vinput, tokens;
    string tmpfile = aurostd::TmpStrCreate();
    aurostd::string2file(lat_atat, tmpfile);
    string sstr = aurostd::execute2string("genstr -n " + aurostd::utype2string<int>(max_num_atoms) + " -l " + tmpfile, stdouterr_fsio);
    aurostd::RemoveFile(tmpfile);
    if (sstr.size() == 0 || aurostd::substring2bool(sstr, "Unable to open lattice file")) {
      message << "Invalid lat.in file";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    else if (aurostd::substring2bool(sstr, "command not found")) {
      message << "Missing genstr program";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::string2vectorstring(sstr, vinput);
    message << "Reading ATAT xstructures from genstr";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    for (size_t line = 0; line < vinput.size(); line++) {
      if (aurostd::substring2bool(vinput[line], "end")) {
        vstr.push_back(xstructure(oss, IOATAT_STR));
        vstr.back().scale = scale;
        vstr.back().ReScale(1.0);
        aurostd::StringstreamClean(oss);
      }
      else {
        oss << vinput[line] << endl;
      }
      pflow::updateProgressBar(line, vinput.size() - 1, *p_oss);
    }
    return vstr;
  }

  /// @brief generates the files for the ATAT program
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::generateFilesForATAT() {
    stringstream message, oss;
    aurostd::string2file(lat_atat, rundirpath + "/lat.in"); // generate lat.in file
    uint msize = mapstr.size();
    message << "Generating ATAT files";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    for (size_t i = 0; i < vstr_ce.size(); i++) {
      aurostd::DirectoryMake(rundirpath + "/" + aurostd::utype2string<uint>(i));
      oss << vstr_ce[i];
      aurostd::string2file(oss.str(), rundirpath + "/" + aurostd::utype2string<uint>(i) + "/str.out");
      aurostd::StringstreamClean(oss);
      if (i < msize && mapstr[i] != -1) {
        aurostd::string2file(aurostd::utype2string<double>(vstr_aflow[mapstr[i]].qm_E_cell) + "\n", rundirpath + "/" + aurostd::utype2string<uint>(i) + "/energy");
      }
      pflow::updateProgressBar(i, vstr_ce.size() - 1, *p_oss);
    }
  }

  /// @brief runs the ATAT program
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::runATAT() {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    string cdirpath = aurostd::getPWD();
    uint iter = 0;
    if (aurostd::substring2bool(aurostd::execute2string("mmaps", stdouterr_fsio), "command not found")) {
      message << "Missing mmaps program";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::execute("touch " + rundirpath + "/stop"); // pre-kill mmaps gracefully
    aurostd::RemoveFile(rundirpath + "/maps_is_running");
    aurostd::RemoveFile(rundirpath + "/maps.log");
    if(chdir(rundirpath.c_str())!=0){ //CO20230507 - compiler warning on gitrunner
      message << "chdir failed";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    string tmpfile = aurostd::TmpStrCreate(), logstring = "";
    aurostd::execute("mmaps -d > " + tmpfile + " 2>&1 &");
    message << "Running ATAT cluster expansion";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    while (!aurostd::substring2bool(logstring, "true and predicted ground states agree") && !aurostd::substring2bool(logstring, "No other ground states")) {
      iter++;
      if (LDEBUG) {cerr << __AFLOW_FUNC__ << " Sleeping, iter=" << iter << endl;}
      if (iter > 60) { // wait 60 times the minimum sleep
        message << "mmaps is taking too long to predict structures, directory = " << rundirpath;
        aurostd::RemoveFile(tmpfile);
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      aurostd::Sleep(min_sleep);
      aurostd::file2string(rundirpath + "/maps.log", logstring);
    }
    aurostd::RemoveFile(tmpfile);
    if(chdir(cdirpath.c_str())!=0){ //CO20230507 - compiler warning on gitrunner
      message << "chdir failed";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  /// @brief reads the coefficient of variation
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::readCVCluster() {
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/maps.log", vinput);
    aurostd::string2tokens(vinput.back(), tokens, " ");
    if (tokens[0] != "Crossvalidation") {
      string message = "ATAT run did not complete";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    cv_cluster = aurostd::string2utype<double>(tokens[2]);
    if (cv_cluster > cv_cut) {
      stringstream message;
      message << "Cluster cross-validation score is above the cut-off, CV_cluster=" << cv_cluster << ", CV_cutoff=" << cv_cut;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  /// @brief calculates the cluster number of atoms
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateNumAtomCluster() {
    num_atom_cluster = xvector<int>(vstr_ce.size());
    int natom = 0;
    for (size_t i = 0; i < vstr_ce.size(); i++) {
      natom = 0;
      for (size_t j = 0; j < vstr_ce[i].num_each_type.size(); j++) {natom += vstr_ce[i].num_each_type[j];}
      num_atom_cluster(i + 1) = natom;
    }
  }

  /// @brief calculates the cluster concentration
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateConcentrationCluster() {
    conc_cluster = xmatrix<double>(vstr_ce.size(), nelem);
    int ie = 0;
    xvector<double> stoich(nelem);
    vector<string> str_elements;
    for (size_t i = 0; i < vstr_ce.size(); i++) {
      if (nelem != vstr_ce[i].stoich_each_type.size()) {
        str_elements = vstr_ce[i].GetElements(true, true);
        stoich.reset();
        for (size_t j = 0; j < nelem; j++) {
          if (aurostd::WithinList(str_elements, elements[j], ie)) {stoich(j + 1) = vstr_ce[i].stoich_each_type[ie];}
        }
      }
      else {
        stoich = aurostd::vector2xvector(aurostd::deque2vector(vstr_ce[i].stoich_each_type));
      }
      for (size_t j = 0; j < nelem; j++) {conc_cluster(i + 1, j + 1) = stoich(j + 1);}
    }
    num_elem_cluster = (double)max_num_atoms * conc_cluster;
  }

  /// @brief reads the cluster excess energy
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::readExcessEnergyCluster() {
    excess_energy_cluster = xvector<double>(conc_cluster.rows);
    int ind;
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/ref_energy.out", vinput);
    xvector<double> nrg_ref = aurostd::vector2xvector(aurostd::vectorstring2vectorutype<double>(vinput));
    // Read energies per atom
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    aurostd::string2tokens(vinput[0], tokens, " ");
    if (nelem > tokens.size()) {
      string message = "fit.out has the wrong output format";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    for (size_t line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ind = aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1;
      excess_energy_cluster(ind) = aurostd::string2utype<double>(tokens[nelem]);
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    aurostd::string2tokens(vinput[0], tokens, " ");
    if (nelem > tokens.size()) {
      string message = "predstr.out has the wrong output format";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    for (size_t line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ind = aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1;
      excess_energy_cluster(ind) = aurostd::string2utype<double>(tokens[nelem]);
    }
    excess_energy_cluster *= max_num_atoms; // per cell
  }

  /// @brief prune the clusters to only include clusters congruent with the maximum number of atoms in the cluster expansion
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::setCongruentClusters() {
    vector<int> index_cluster;
    for (int i = num_atom_cluster.lrows; i <= num_atom_cluster.urows; i++) {
      if (!(max_num_atoms % num_atom_cluster(i))) {index_cluster.push_back(i);}
    }
    ncluster = index_cluster.size();
    vector<xstructure> _vstr_ce(ncluster);
    xmatrix<double> m1(ncluster, nelem), m2(ncluster, nelem);
    xvector<int> v1(ncluster);
    xvector<double> v2(ncluster);
    for (size_t i = 0; i < index_cluster.size(); i++) {
      _vstr_ce[i] = vstr_ce[index_cluster[i] - 1];
      m1.setmat(conc_cluster.getxmat(index_cluster[i], index_cluster[i], 1, nelem), i + 1, 1);
      m2.setmat(num_elem_cluster.getxmat(index_cluster[i], index_cluster[i], 1, nelem), i + 1, 1);
      v1(i + 1) = num_atom_cluster(index_cluster[i]);
      v2(i + 1) = excess_energy_cluster(index_cluster[i]);
    }
    vstr_ce = _vstr_ce;
    conc_cluster = m1;
    num_elem_cluster = m2;
    num_atom_cluster = v1;
    excess_energy_cluster = v2;
  }

  /// @brief calculates the cluster degeneracy
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateDegeneracyCluster() {
    degeneracy_cluster = xvector<unsigned long long int>(ncluster);
    string filepath = "";
    if (!rundirpath.empty()) {
      vector<string> tokens;
      aurostd::string2tokens(rundirpath, tokens, "/");
      filepath = "/" + QCA_FILE_PREFIX + "degen_" + aurostd::utype2string<int>(max_num_atoms) + "_atom.txt";
      for (int i = tokens.size() - 2; i >= 0; i--) {filepath = "/" + tokens[i] + filepath;} // well defined path due to how we define rundirpath
      if (aurostd::file2vectorstring(filepath, tokens)) {
        degeneracy_cluster = aurostd::vector2xvector(aurostd::vectorstring2vectorutype<unsigned long long int>(tokens));
        return;
      }
    }
    stringstream message;
    message << "Generating derivative structures";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    xstructure str_pocc, str;
    deque<_atom> patoms, atoms;
    _atom atom;
    uint itype;
    xmatrix<double> lattice(3,3);
    // Define primitive parent lattice and pocc atoms
    if (plattice == "fcc") {
      lattice(1, 1) = 0.0; lattice(1, 2) = 0.5; lattice(1, 3) = 0.5;
      lattice(2, 1) = 0.5; lattice(2, 2) = 0.0; lattice(2, 3) = 0.5;
      lattice(3, 1) = 0.5; lattice(3, 2) = 0.5; lattice(3, 3) = 0.0;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      patoms.push_back(atom);
    }
    else if (plattice == "bcc") {
      lattice(1, 1) = -0.5; lattice(1, 2) = 0.5; lattice(1, 3) = 0.5;
      lattice(2, 1) = 0.5; lattice(2, 2) = -0.5; lattice(2, 3) = 0.5;
      lattice(3, 1) = 0.5; lattice(3, 2) = 0.5; lattice(3, 3) = -0.5;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      patoms.push_back(atom);
    }
    else if (plattice == "hcp") {
      double sqrt3 = std::sqrt(3.0), sqrt8 = std::sqrt(8.0);
      lattice(1, 1) = 1.0; lattice(1, 2) = 0.0; lattice(1, 3) = 0.0;
      lattice(2, 1) = -0.5; lattice(2, 2) = 0.5 * sqrt3; lattice(2, 3) = 0.0;
      lattice(3, 1) = 0.0; lattice(3, 2) = 0.0; lattice(3, 3) = sqrt8 / sqrt3;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      patoms.push_back(atom);
      atom.clear();
      atom.fpos(1) = 2.0 / 3.0; atom.fpos(2) = 1.0 / 3.0; atom.fpos(3) = 0.5;
      patoms.push_back(atom);
    }
    aurostd::xcombos xcpocc(max_num_atoms + 1, elements.size(), 'E', false);
    vector<xvector<double>> vpocc;
    while (xcpocc.increment()) {
      if (aurostd::sum(xcpocc.getCombo()) == max_num_atoms) {
        vpocc.push_back(aurostd::xvector2utype<int, double>(aurostd::vector2xvector(xcpocc.getCombo())) / (double)max_num_atoms);
      }
    }
    vector<xstructure> vstr_ds, vstr_u;
    long int sum_calculated = 0, sum_accepted = 0;
    for (size_t i = 0; i < vpocc.size(); i++) {
      str_pocc.clear();
      str_pocc.lattice = lattice;
      str_pocc.partial_occupation_HNF = max_num_atoms;
      str_pocc.partial_occupation_flag = true;
      str_pocc.neg_scale_second = true;
      atoms.clear();
      for (size_t j = 0; j < elements.size(); j++) {
        for (size_t k = 0; k < patoms.size(); k++) {
          atom.clear();
          atom.name = atom.cleanname = elements[j];
          atom.name_is_given = true;
          atom.partial_occupation_flag = true;
          atom.partial_occupation_value = vpocc[i](j + 1);
          atom.fpos = patoms[k].fpos;
          atom.cpos = str_pocc.f2c * patoms[k].fpos;
          atoms.push_back(atom);
        }
      }
      std::stable_sort(atoms.begin(), atoms.end(), sortAtomsNames);
      //CO20230319 - DX, I believe this code is obsolete now that I patch AddAtom() with add_species, keep for safety
      itype = 0;
      atoms[0].type = itype;
      for (size_t j = 1; j < atoms.size(); j++) {
        if (atoms[j].name != atoms[j - 1].name) {itype++;}
        atoms[j].type = itype;
      }
      str_pocc.AddAtom(atoms,true); //CO20230319 - add by species
      // Enumerate unique superlattices
      pocc::POccCalculator pcalc(str_pocc);
      pcalc.calculateHNF();
      pcalc.getTotalPermutationsCount();
      pcalc.calculate();
      vstr_u = pcalc.getUniqueDerivativeStructures();
      for (vector<xstructure>::iterator it = vstr_u.begin(); it != vstr_u.end(); it++) {vstr_ds.push_back(*it);}
      sum_accepted += pcalc.total_permutations_count;
    }
    message << "Generated all unique derivative structures, count = " << vstr_ds.size();
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    // Find degenerate structures
    message << "Finding the degeneracy of clusters";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vstr_ds.insert(vstr_ds.begin(), vstr_ce.begin(), vstr_ce.end()); // concatenate xstructures
    XtalFinderCalculator xtal_calc;
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
    comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS", true);
    vector<vector<uint>> vindex = xtal_calc.groupSimilarXstructures(vstr_ds, comparison_options); // costly part of the function
    for (size_t i = 0; i < vindex.size(); i++) {
      std::sort(vindex[i].begin(), vindex[i].end()); // first index is the cluster index
      if (vindex[i][0] < ncluster) {
        for (size_t j = 1; j < vindex[i].size(); j++) {
          degeneracy_cluster(vindex[i][0] + 1) += pocc::getDGFromXStructureTitle(vstr_ds[vindex[i][j]].title);
        }
      }
    }
    sum_calculated = aurostd::sum(degeneracy_cluster);
    if (!aurostd::isequal(sum_calculated, sum_accepted)) {
      message << "Degeneracies do not satisfy the sum rule, the sum is " << sum_calculated << " but should be " << sum_accepted;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    if (!filepath.empty()) {
      stringstream ss;
      for (int i = degeneracy_cluster.lrows; i <= degeneracy_cluster.urows; i++) {ss << degeneracy_cluster(i) << endl;}
      aurostd::stringstream2file(ss, filepath);
    }
  }

  /// @brief calculates the macroscopic concentration
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateConcentrationMacro() {
    if (!conc_curve_range.empty()) { // curve in concentration space
      xmatrix<double> _conc_macro(conc_npts, nelem);
      for (uint i = 0; i < nelem; i++) {
        _conc_macro.setcol(aurostd::linspace(conc_curve_range[2 * i], conc_curve_range[2 * i + 1], conc_npts), i + 1);
      }
      conc_macro = _conc_macro;
    }
    else {
      if (nelem == 2) { // trivial
        xmatrix<double> _conc_macro(conc_npts, nelem);
        xvector<double> x = aurostd::linspace(CONC_SHIFT, 1.0 - CONC_SHIFT, conc_npts);
        _conc_macro.setcol(x, 1);
        _conc_macro.setcol(1.0 - x, 2);
        conc_macro = _conc_macro;
      }
      else if (nelem == 3) { // barycentric coordinates
        xmatrix<double> _conc_macro((conc_npts * conc_npts - conc_npts) / 2.0, nelem); // always even
        xvector<double> x = aurostd::linspace(CONC_SHIFT, 1.0 - CONC_SHIFT, conc_npts);
        xvector<double> p1(2), p2(2), p3(2), p4(2);
        p1(1) = 0.0; p1(2) = 0.0;
        p2(1) = 1.0; p2(2) = 0.0;
        p3(1) = 1.0; p3(2) = 1.0;
        int ii = 0;
        for (int i = x.lrows; i <= x.rows; i++) {
          for (int j = x.lrows + i; j <= x.rows; j++) {
            ii++;
            p4(1) = x(j); p4(2) = x(i);
            _conc_macro(ii, 1) = (p1(1) - p3(1)) * (p4(2) - p1(2)) - (p1(1) - p4(1)) * (p3(2) - p1(2));
            _conc_macro(ii, 2) = (p1(1) - p4(1)) * (p2(2) - p1(2)) - (p1(1) - p2(1)) * (p4(2) - p1(2));
            _conc_macro(ii, 3) = (p4(1) - p3(1)) * (p2(2) - p4(2)) - (p4(1) - p2(1)) * (p3(2) - p4(2));
          }
        }
        conc_macro = _conc_macro;
      }
      else { // for greater than ternary only evaluate at equi-concentration
        xmatrix<double> _conc_macro(1, nelem);
        _conc_macro.set(1.0 / (double)nelem);
        conc_macro = _conc_macro;
      }
    }
    nconc = conc_macro.rows;
  }

  /// @brief calculates the temperature profile
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateTemperatureRange() {
    temp = aurostd::linspace(temp_range[0], temp_range[1], temp_npts);
    beta = aurostd::pow(KBOLTZEV * temp, -1.0);
  }

  /// @brief calculates the ideal (high-T) cluster probability
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateProbabilityIdealCluster() {
    prob_ideal_cluster = xmatrix<double>(nconc, ncluster);
    for (int ix = conc_macro.lrows; ix <= conc_macro.urows; ix++) {
      for (int ic = conc_cluster.lrows; ic <= conc_cluster.urows; ic++) {
        prob_ideal_cluster(ix, ic) = (double)degeneracy_cluster(ic);
        for (int ie = conc_macro.lcols; ie <= conc_macro.ucols; ie++) {
          prob_ideal_cluster(ix, ic) *= std::pow(conc_macro(ix, ie), conc_cluster(ic, ie) * max_num_atoms);
        }
      }
      prob_ideal_cluster.setmat(prob_ideal_cluster.getxmat(ix, ix, conc_cluster.lrows, conc_cluster.urows) /
                                aurostd::sum(prob_ideal_cluster.getxmat(ix, ix, conc_cluster.lrows, conc_cluster.urows)), ix, 1); // normalize sum to 1
    }
  }

  /// @brief checks the ideal (high-T) probability for errors
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::checkProbabilityIdeal() const {
    for (int ix = prob_ideal_cluster.lrows; ix <= prob_ideal_cluster.urows; ix++) {
      if (!aurostd::isequal(aurostd::sum(prob_ideal_cluster(ix)), 1.0)) { // unnormalized
        stringstream message;
        message << "Ideal solution (high-T) probability is unnormalized for ix=" << ix << " | SUM[P_cluster]=" << aurostd::sum(prob_ideal_cluster(ix));
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
      else if (!aurostd::isequal(prob_ideal_cluster(ix) * conc_cluster, conc_macro(ix))) { // does not satisfy concentration constraints
        stringstream message;
        message << "Ideal solution (high-T) probability does not satisfy concentration contraint X=[" << conc_macro(ix) << " ], X_calc=[" << prob_ideal_cluster(ix) * conc_cluster << "]";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
    }
  }

  /// @brief calculates the equilibrium cluster probability of a 1-D system at a particular concentration
  ///
  /// @param iix shifted concentration index
  /// @param it temperature index
  ///
  /// @note helper function for parallelization
  ///
  /// @authors
  /// @mod{SD,20220812,created function}
  void QuasiChemApproxCalculator::calculateProbabilityCluster1D(int iix, const int it) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    int ix = iix + prob_ideal_cluster.lrows;
    xvector<double> soln(nelem - 1), coeff(max_num_atoms + 1), rootr(max_num_atoms), rooti(max_num_atoms);
    bool found_soln = false;
    for (int ic = prob_ideal_cluster.lcols; ic <= prob_ideal_cluster.ucols; ic++) {
      coeff((int)num_elem_cluster(ic, 1) + 1) += prob_ideal_cluster(ix, ic) * std::exp(-beta(it) * excess_energy_cluster(ic)) * (conc_cluster(ic, 1) - conc_macro(ix, 1));
    }
    if (LDEBUG) {cerr << __AFLOW_FUNC__ << " it=" << it << " ix=" << ix << " | coeff=" << coeff << endl;}
    aurostd::polynomialFindRoots(coeff, rootr, rooti);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "   real roots=" << rootr << endl;
      cerr << __AFLOW_FUNC__ << "   imag roots=" << rooti << endl;
    }
    for (int isol = rootr.lrows; isol <= rootr.urows; isol++) {
      if (rootr(isol) > soln(1) &&
          aurostd::isequal(rooti(isol), 0.0) &&
          max_num_atoms * std::pow(rootr(isol), max_num_atoms) != INFINITY) { // solution must be positive, real, and finite
        soln(1) = rootr(isol);
        found_soln = true;
      }
    }
    if (!found_soln) { // physical solution does not exist
      stringstream message;
      message << "Physical equilibrium probability does not exist for T=" << temp(it) << "K, X=[" << conc_macro(ix) << " ]" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    for (int ic = prob_ideal_cluster.lcols; ic <= prob_ideal_cluster.ucols; ic++) {
      prob_cluster[it - 1](ix, ic) = prob_ideal_cluster(ix, ic) * std::exp(-beta(it) * excess_energy_cluster(ic)) * std::pow(soln(1), num_elem_cluster(ic, 1));
    }
    prob_cluster[it - 1].setmat(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols) /
                                aurostd::sum(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols)), ix, 1);
  }

  /// @brief calculates the equilibrium cluster probability of a N-D system at a particular concentration
  ///
  /// @param iix shifted concentration index
  /// @param it temperature index
  ///
  /// @note helper function for parallelization
  ///
  /// @authors
  /// @mod{SD,20220812,created function}
  void QuasiChemApproxCalculator::calculateProbabilityClusterND(int iix, const int it) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    int ix = iix + prob_ideal_cluster.lrows;
    xvector<double> soln(nelem - 1);
    xmatrix<double> msoln;
    vector<std::function<double(xvector<double>)>> vpoly, vdpoly;
    vector<vector<std::function<double(xvector<double>)>>> jac;
    bool found_soln = false;
    for (int ieq = 1; ieq <= (int)nelem - 1; ieq++) {
      vpoly.push_back([this, it, ix, ieq](xvector<double> xvar) {return getProbabilityConstraint(it, ix, ieq, 0, xvar);});
      for (int ideq = 1; ideq <= (int)nelem - 1; ideq++) {
        vdpoly.push_back([this, it, ix, ieq, ideq](xvector<double> xvar) {return getProbabilityConstraint(it, ix, ieq, ideq, xvar);});
      }
      jac.push_back(vdpoly);
      vdpoly.clear();
    }
    if (LDEBUG) {cerr << __AFLOW_FUNC__ << " it=" << it << " ix=" << ix << endl;}
    found_soln = aurostd::findZeroDeflation(soln0.getcol(ix), vpoly, jac, msoln);
    if (LDEBUG) {cerr << __AFLOW_FUNC__ << "   exist=" << found_soln << " | real roots=" << endl << msoln << endl;}
    if (found_soln) { // if solution is found, it must also be positive, real, and finite
      for (int isol = msoln.lcols; isol <= msoln.ucols; isol++) {
        found_soln = true;
        for (int ieq = 1; ieq <= (int)nelem - 1 && found_soln; ieq++) {
          found_soln = !(msoln.getcol(isol)(ieq) <= 0.0);
        }
        if (found_soln) {
          soln = msoln.getcol(isol);
          break;
        }
      }
    }
    if (!found_soln) {
      stringstream message;
      message << "Physical equilibrium probability does not exist for T=" << temp(it) << "K, X=[" << conc_macro(ix) << " ]" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    soln0.setcol(soln, ix); // use the higher temperature solution as an initial guess for lower temperature solution
    for (int ic = prob_ideal_cluster.lcols; ic <= prob_ideal_cluster.ucols; ic++) {
      prob_cluster[it - 1](ix, ic) = prob_ideal_cluster(ix, ic) * std::exp(-beta(it) * excess_energy_cluster(ic)) *
                                    aurostd::elements_product(aurostd::pow(soln, num_elem_cluster.getxvec(ic, ic, 1, nelem - 1)));
    }
    prob_cluster[it - 1].setmat(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols) /
                                aurostd::sum(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols)), ix, 1); // normalize sum to 1
  }

  /// @brief calculates the equilibrium cluster probability
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateProbabilityCluster() {
    stringstream message;
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt;
    std::function<void(int, const int&)> calculateProbabilityCluster1D_MT = std::bind(&QuasiChemApproxCalculator::calculateProbabilityCluster1D, this, _1, _2);
    std::function<void(int, const int&)> calculateProbabilityClusterND_MT = std::bind(&QuasiChemApproxCalculator::calculateProbabilityClusterND, this, _1, _2);
#endif
    prob_cluster.clear();
    prob_cluster.assign(temp.rows, xmatrix<double>(nconc, ncluster)); // initialize
    message << "Calculating cluster probabilities";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    if (nelem - 1 == 1) {
      for (int it = temp.lrows; it <= temp.urows; it++) {
#ifdef AFLOW_MULTITHREADS_ENABLE
        xt.run(prob_ideal_cluster.rows, calculateProbabilityCluster1D_MT, it);
#else
        for (int iix = 0; iix < prob_ideal_cluster.rows; iix++) {calculateProbabilityCluster1D(iix, it);}
#endif
        pflow::updateProgressBar(it - temp.lrows + 1, temp.rows, *p_oss);
      }
    }
    else {
      soln0 = aurostd::ones_xm<double>(nelem - 1, nconc);
      for (int it = temp.urows; it >= temp.lrows; it--) { // go backwards
#ifdef AFLOW_MULTITHREADS_ENABLE
        xt.run(prob_ideal_cluster.rows, calculateProbabilityClusterND_MT, it);
#else
        for (int iix = 0; iix < prob_ideal_cluster.rows; iix++) {calculateProbabilityClusterND(iix, it);}
#endif
        pflow::updateProgressBar(temp.urows - it + 1, temp.rows, *p_oss);
      }
    }
  }

  /// @brief gets the value of the probability constraint for a particular set of indices
  ///
  /// @param it temperature index
  /// @param ix macroscopic concentration index
  /// @param ie element index
  /// @param ideq index of the equation that is differentiated
  /// @param xvar concentration variable
  ///
  /// @return value of the constraint for a particular set of indices
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  double QuasiChemApproxCalculator::getProbabilityConstraint(const int it, const int ix, const int ie, const int ideq, const xvector<double>& xvar) const {
    double totsum = 0.0, prodx = 1.0;
    for (int ic = prob_ideal_cluster.lcols; ic <= prob_ideal_cluster.ucols; ic++) {
      prodx = 1.0;
      for (int ieq = 1; ieq <= (int)nelem - 1; ieq++) {
        if (ieq == ideq) {
          prodx *= num_elem_cluster(ic, ieq) * std::pow(xvar(ieq), num_elem_cluster(ic, ieq) - 1);
        }
        else {
          prodx *= std::pow(xvar(ieq), num_elem_cluster(ic, ieq));
        }
      }
      totsum += prob_ideal_cluster(ix, ic) * std::exp(-beta(it) * excess_energy_cluster(ic)) * (conc_cluster(ic, ie) - conc_macro(ix, ie)) * prodx;
    }
    return totsum;
  }


  /// @brief checks the equilibrium probability for errors
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::checkProbabilityEquilibrium() const {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    double diff = AUROSTD_MAX_DOUBLE, diff_old;
    for (int it = temp.lrows; it <= temp.urows; it++) {
      for (int ix = prob_ideal_cluster.lrows; ix <= prob_ideal_cluster.urows; ix++) {
        if (!aurostd::isequal(aurostd::sum(prob_cluster[it - 1](ix)), 1.0)) { // unnormalized
          stringstream message;
          message << "Equilibrium probability is unnormalized for T=" << temp(it) << "K, ix=" << ix << " | SUM[P_cluster]=" << aurostd::sum(prob_cluster[it - 1](ix));
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
        else if (!aurostd::isequal(prob_cluster[it - 1](ix) * conc_cluster, conc_macro(ix))) { // does not satisfy concentration constraints
          stringstream message;
          message << "Equilibrium probability does not satisfy concentration contraint for T=" << temp(it) << "K, X=[" << conc_macro(ix) << " ], X_calc=[" << prob_cluster[it - 1](ix) * conc_cluster << "]";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
      }
      diff_old = diff;
      diff = aurostd::sum(aurostd::abs(prob_cluster[it - 1] - prob_ideal_cluster));
      if (LDEBUG) {cerr << __AFLOW_FUNC__ << " |P - P0| = " << diff << endl;}
      if (diff > diff_old) { // P(T_inf) does not equal P0(T_inf)
        stringstream message;
        message << "Equilibrium probability does not converge to the high-T limit for T=" << temp(it) << "K";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
  }

  /// @brief calculates the relative entropy evaluated at equi-concentration and the transition temperature
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateRelativeEntropyEC() {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    int n_fit = 8;
    nconc = 1;
    xmatrix<double> conc_macro_orig = conc_macro; // store original macroscopic concentration
    conc_macro = xmatrix<double>(1, nelem);
    conc_macro.set(1.0 / (double)nelem);
    calculateProbabilityIdealCluster();
    // Check that the probability is physical, if not, shift the temperature range upward
    double dtemp = temp(2) - temp(1);
    message << "Checking for physical probabilities at equi-composition";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    while (aurostd::min(temp) < DEFAULT_QCA_TEMP_MIN_LIMIT) {
      try {
        calculateProbabilityCluster();
        break;
      }
      catch (aurostd::xerror& err) {
        message << "Could not find physical solution at equi-concentration, increasing the temperature range";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        temp += dtemp;
        beta = aurostd::pow(KBOLTZEV * temp, -1.0);
      }
    }
    if (aurostd::min(temp) > DEFAULT_QCA_TEMP_MIN_LIMIT) {
      message << "Could not find a temperature at equi-concentration that leads to a physical solution";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    xvector<double> temp_orig = temp; // store original temperature
    xvector<double> beta_orig = beta; // store original inverse temperature
    message << "Equi-composition temperature range = [" << temp(1) << "K, " << temp(temp.urows) << "K]";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    xvector<double> order_param(temp.rows);
    xmatrix<double> m1, m2, m3 = prob_ideal_cluster * aurostd::trasp(prob_ideal_cluster);
    for (int it = temp.lrows; it <= temp.urows; it++) {
      m1 = prob_cluster[it - 1] * aurostd::trasp(prob_ideal_cluster);
      m2 = prob_cluster[it - 1] * aurostd::trasp(prob_cluster[it - 1]);
      order_param(it) = m1(1, 1) / (aurostd::sqrt(m2(1, 1)) * aurostd::sqrt(m3(1, 1)) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
    }
    // Check for divergence of the order parameter at low temperature
    bool low_t_div = false;
    if (order_param.rows > 3) {
      double d1 = (order_param(3) - order_param(1)) / (temp(3) - temp(1)), d2 = (order_param(4) - order_param(2)) / (temp(4) - temp(2));
      if (std::abs(d1) > std::abs(d2)) {low_t_div = true;}
    }
    double temp_mean = aurostd::mean(temp), temp_std = aurostd::stddev(temp);
    xvector<double> temp_scaled = (temp - temp_mean) / temp_std; // scale for numerical stability
    xvector<double> wts = aurostd::ones_xv<double>(order_param.rows);
    xvector<double> coeff = aurostd::polynomialCurveFit(temp_scaled, order_param, n_fit, wts);
    xvector<double> coeff1 = aurostd::evalPolynomialCoeff(coeff, 1), coeff2 = aurostd::evalPolynomialCoeff(coeff, 2);
    // If the order parameter diverges, then we can interpolate
    if (low_t_div) {
      message << "Order parameter diverges at low temperature, using interpolation";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      xvector<double> dadt = aurostd::evalPolynomial_xv(temp_scaled, coeff1), temp_interp(temp.rows + 1), dadt_interp(temp.rows + 1);
      for (int i = temp.lrows; i <= temp.urows; i++) {temp_interp(i + 1) = temp(i);} // add T = 0K
      for (int i = dadt.lrows; i <= dadt.urows; i++) {dadt_interp(i + 1) = dadt(i);} // add Da/DT(0K) = 0
      temp_mean = aurostd::mean(temp_interp), temp_std = aurostd::stddev(temp_interp);
      temp_scaled = (temp_interp - temp_mean) / temp_std;
      wts = aurostd::exp(temp_scaled); // scale weights by temperature, since high-T is more accurate
      wts(1) = wts(wts.rows); // force convergence at T = 0K
      coeff1 = aurostd::polynomialCurveFit(temp_scaled, dadt_interp, n_fit - 1, wts);
      coeff2 = aurostd::evalPolynomialCoeff(coeff1, 1);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " alpha_orig=" << order_param << endl;
      cerr << __AFLOW_FUNC__ << " D[alpha, 0]=" << aurostd::evalPolynomial_xv(temp_scaled, coeff) << endl;
      cerr << __AFLOW_FUNC__ << " D[alpha, 1]=" << aurostd::evalPolynomial_xv(temp_scaled, coeff1) << endl;
      cerr << __AFLOW_FUNC__ << " D[alpha, 2]=" << aurostd::evalPolynomial_xv(temp_scaled, coeff2) << endl;
    }
    xvector<double> rootr(n_fit - 2), rooti(n_fit - 2);
    aurostd::polynomialFindRoots(coeff2, rootr, rooti);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " coeff2=" << coeff2 << endl;
      cerr << __AFLOW_FUNC__ << "   Real roots=" << rootr << endl;
      cerr << __AFLOW_FUNC__ << "   Imag roots=" << rooti << endl;
    }
    bool unset = true;
    for (int j = 1; j <= n_fit; j++) {
      if (rootr(j) >= temp_scaled(1) && rootr(j) <= temp_scaled(temp.rows) && aurostd::isequal(rooti(j), 0.0)) { // solution must be real and within temp range
        if (unset) {
          param_ec.second = rootr(j);
          unset = false;
        }
        else if (aurostd::evalPolynomial(param_ec.second, coeff1) < aurostd::evalPolynomial(rootr(j), coeff1)) { // keep largest gradient
          param_ec.second = rootr(j);
        }
      }
    }
    param_ec.second = temp_std * param_ec.second + temp_mean;
    message << "Equi-composition transition temperature = " << param_ec.second << "K";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    temp = xvector<double> {param_ec.second};
    beta = aurostd::pow(KBOLTZEV * temp, -1.0);
    calculateProbabilityCluster();
    for (int ic = prob_ideal_cluster.lcols; ic <= prob_ideal_cluster.ucols; ic++) {
      param_ec.first += prob_cluster[0](1, ic) * aurostd::log((prob_cluster[0](1, ic) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) / (prob_ideal_cluster(1, ic) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_));
    }
    // Revert variables to the original
    nconc = conc_macro_orig.rows; 
    conc_macro = conc_macro_orig;
    temp = temp_orig;
    beta = beta_orig;
  }

  /// @brief calculates the relative entropy.
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateRelativeEntropy() {
    rel_s = xmatrix<double>(nconc, temp_npts);
    for (int ix = prob_ideal_cluster.lrows; ix <= prob_ideal_cluster.urows; ix++) {
      for (int it = temp.lrows; it <= temp.urows; it++) {
        for (int ic = prob_ideal_cluster.lcols; ic <= prob_ideal_cluster.ucols; ic++) {
          rel_s(ix, it) += prob_cluster[it - 1](ix, ic) * aurostd::log((prob_cluster[it - 1](ix, ic) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) / (prob_ideal_cluster(ix, ic) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_));
        }
      }
    }
  }

  /// @brief calculates the binodal curve
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::calculateBinodalCurve() {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    int n_fit = 8;
    binodal_curve = xvector<double>(nconc);
    xvector<double> coeff, rootr(n_fit), rooti(n_fit);
    xvector<double> wts = aurostd::ones_xv<double>(temp_npts);
    double temp_mean = aurostd::mean(temp), temp_std = aurostd::stddev(temp);
    xvector<double> temp_scaled = (temp - temp_mean) / temp_std; // scale for numerical stability
    message << "Calculating binodal curve";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    for (int ix = rel_s.lrows; ix <= rel_s.urows; ix++) {
      coeff = aurostd::polynomialCurveFit(temp_scaled, rel_s(ix) - param_ec.first, n_fit, wts);
      aurostd::polynomialFindRoots(coeff, rootr, rooti);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ix=" << ix << " | coeff=" << coeff << endl;
        cerr << __AFLOW_FUNC__ << "   Real roots=" << rootr << endl;
        cerr << __AFLOW_FUNC__ << "   Imag roots=" << rooti << endl;
      }
      rootr = temp_std * rootr + temp_mean;
      for (int isol = 1; isol <= n_fit; isol++) {
        if (rootr(isol) >= temp(1) && rootr(isol) <= temp(temp.rows) && aurostd::isequal(rooti(isol), 0.0) && binodal_curve(ix) < rootr(isol)) { // largest solution must be real and within temp range
          binodal_curve(ix) = rootr(isol);
        }
      }
      pflow::updateProgressBar(ix - rel_s.lrows + 1, rel_s.rows, *p_oss);
    }
  }

  /// @brief write the QCA variables to screen or to file
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::writeData() {
    string filepath = rundirpath + "/" + QCA_FILE_PREFIX + "output." + print;
    stringstream output;
    if (print == "txt") {
      output << AFLOWIN_SEPARATION_LINE << endl;
      output << QCA_AFLOW_TAG << "START" << endl;
      output << "ALLOY=" << alloyname << endl;
      output << "PLATTICE=" << plattice << endl;
      output << "CONC_MACRO=" << endl << conc_macro << endl;
      output << "TEMP_RANGE=" << temp_range[0] << "," << temp_range[1] << endl;
      output << "MAX_ATOMS_CELL=" << max_num_atoms << endl;
      output << "CV=" << cv_cluster << endl;
      output << "TEMP_EC=" << param_ec.second << endl;
      output << "BINODAL=" << endl << trasp(binodal_curve) << endl;
      output << QCA_AFLOW_TAG << "END" << endl;
      output << AFLOWIN_SEPARATION_LINE << endl;
      if (!screen_only) {
        aurostd::stringstream2file(output, filepath);
        return;
      }
    }
    else if (print == "json") {
      aurostd::JSON::object jo(aurostd::JSON::object_types::DICTIONARY);
      jo["Alloy name"] = alloyname;
      jo["Elements"] = elements;
      jo["Parent lattice"] = plattice;
      jo["Concentration curve"] = conc_curve;
      jo["Macroscopic concentration"] = conc_macro;
      jo["Temperature range (K)"] = temp;
      jo["Max atoms per cell"] = max_num_atoms;
      jo["Cluster CV score (eV)"] = cv_cluster;
      jo["Cluster degeneracy"] = degeneracy_cluster;
      jo["Cluster concentration"] = conc_cluster;
      jo["Cluster excess energy (eV)"] = excess_energy_cluster;
      jo["EC transition temperature (K)"] = param_ec.second;
      jo["Binodal curve (K)"] = binodal_curve;
      jo.saveFile(filepath);
      return;
    }
    cout << output.str() << endl;
  }

  /// @brief reads the QCA variables from file
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::readData() {
    string filepath = rundirpath + "/" + QCA_FILE_PREFIX + "output.json";
    if (!aurostd::FileExist(filepath)) {
      string message = "The JSON file does not exist, filepath=" + filepath; 
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    aurostd::JSON::object jo = aurostd::JSON::loadFile(filepath);
    alloyname = (string)jo["Alloy name"];
    conc_curve = bool(jo["Concentration curve"]);
    conc_macro = jo["Macroscopic concentration"];
    temp = jo["Temperature range (K)"];
    binodal_curve = jo["Binodal curve (K)"];
  }

  /// @brief plots the QCA data
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void QuasiChemApproxCalculator::plotData() {
    if (print == "txt" || print == "json") {return;}
    stringstream gpfile;
    aurostd::xoption plotoptions;
    vector<vector<double>> data;
    // BEGIN - set standard
    plotoptions.push_attached("BACKGROUND_COLOR", "#FFFFFF");
    plotoptions.push_attached("GRID_COLOR", "#808080");
    plotoptions.push_attached("GRID_WIDTH", "1");
    plotoptions.push_attached("GRID_LINE_TYPE", "0");
    plotoptions.push_attached("PLOT_SIZE", "5.333, 3");
    // END - set standard
    // START - set default
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    plotoptions.push_attached("FILE_NAME", rundirpath + "/" + QCA_FILE_PREFIX + "plot");
    plotoptions.push_attached("FILE_NAME_LATEX", rundirpath + "/" + QCA_FILE_PREFIX + "plot");
    plotoptions.push_attached("IMAGE_FORMAT", print);
    plotoptions.push_attached("PLOT_TITLE", alloyname);
    plotoptions.push_attached("XLABEL", "Concentration");
    plotoptions.push_attached("XMIN", "0.0");
    plotoptions.push_attached("XMAX", "1.0");
    plotoptions.push_attached("XUNIT", "");
    plotoptions.push_attached("YLABEL", "Temperature");
    plotoptions.push_attached("YMIN", "1.0");
    plotoptions.push_attached("YMAX", aurostd::utype2string(1.2 * aurostd::max(binodal_curve)));
    plotoptions.push_attached("YUNIT", "K");
    plotoptions.push_attached("TITLES", "Binodal");
    plotoptions.push_attached("COLORS", "#000000");
    plotoptions.push_attached("LINETYPES", "-1");
    // END - set default
    plotter::generateHeader(gpfile, plotoptions, false);
    if (conc_curve) { // generic x-axis for the concentration curve
      xvector<double> lsc = aurostd::linspace(0.0, 1.0, conc_macro.rows);
      for (int ix = conc_macro.lrows; ix <= conc_macro.urows; ix++) {data.push_back({lsc(ix), binodal_curve(ix)});}
        // START - custom
        gpfile << "# Custom" << endl;
        gpfile << "set label '\\shortstack{c}{$" << elements[0] << "$=" << conc_macro(conc_macro.lrows, 1);
        for (size_t i = 0; i < nelem; i++) {
          gpfile << "\\\\$" << elements[i] << "$=" << conc_macro(conc_macro.lrows, i + 1);
        }
        gpfile << "}' at 0.0,0.0 center offset 0.0,-2.5" << endl;
        gpfile << "set label '\\shortstack{c}{$" << elements[0] << "$=" << conc_macro(conc_macro.urows, 1);
        for (size_t i = 0; i < nelem; i++) {
          gpfile << "\\\\$" << elements[i] << "$=" << conc_macro(conc_macro.urows, i + 1);
        }
        gpfile << "}' at 0.0,0.0 center offset 0.0,-2.5" << endl;
        // END - custom
    }
    else if (elements.size() == 2) {
      for (int ix = conc_macro.lrows; ix <= conc_macro.urows; ix++) {data.push_back({conc_macro(ix, 1), binodal_curve(ix)});}
      // START - custom
      gpfile << "# Custom" << endl;
      gpfile << "set label '$" << elements[0] << "$' at 0.0,0.0 center offset 0.0,-2.5" << endl;
      gpfile << "set label '$" << elements[1] << "$' at 1.0,0.0 center offset 0.0,-2.5" << endl;
      // END - custom
    }
    else {
      return; // no plotting support for alloys higher than binary
    }
    plotter::generatePlotGNUPLOT(gpfile, plotoptions, data);
    if (screen_only) {
      cout << gpfile.str() << endl;
    }
    else {
      plotter::savePlotGNUPLOT(plotoptions, gpfile);
    }
  }
}

namespace qca {
  /// @brief modifies the QCA variables based on the input flags and starts the module
  ///
  /// @param vpflow command line options
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void quasiChemicalApprox(const aurostd::xoption& vpflow) {
    if (vpflow.flag("QCA::USAGE")) {
      if (!vpflow.flag("QCA::SCREEN_ONLY")) {displayUsage();}
      return;
    }
    qca::QuasiChemApproxCalculator qca_calc;
    qca_calc.initialize(vpflow);
    qca_calc.errorFix();
    qca_calc.printParams();
    // Only plot data from JSON file. This is useful when the plotting routine cannot be run on the machine
    // running the calculation.
    if (vpflow.flag("QCA::IMAGE_ONLY")) {
      qca_calc.readData();
      qca_calc.plotData();
      return;
    }
    qca_calc.calculateBinodal();
    // Write and plot results
    if (vpflow.flag("QCA::BINODAL")) {
      qca_calc.writeData();
      qca_calc.plotData();
    }
  }

  /// @brief displays the usage commands and options
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void displayUsage(void) {
    vector<string> usage_options;
    usage_options.push_back("aflow --quasi_chem_approx|--qca --plattice=|--plat=fcc --elements=|--elem=Au,Pt[,Zn] [qca_options] [--directory=[DIRECTORY]]");
    usage_options.push_back(" ");
    usage_options.push_back("qca_options:");
    usage_options.push_back(" ");
    usage_options.push_back("GENERAL OPTIONS:");
    usage_options.push_back("--usage");
    usage_options.push_back("--screen_only");
    usage_options.push_back("--image_only|--image");
    usage_options.push_back("--aflowlib_directory=|--aflowlib_dir=...");
    usage_options.push_back("--print=|--p=|--output=|--o=txt");
    usage_options.push_back(" ");
    usage_options.push_back("BINODAL OPTIONS:");
    usage_options.push_back("--binodal");
    usage_options.push_back("--use_sg");
    usage_options.push_back("--aflow_max_num_atoms=4");
    usage_options.push_back("--max_num_atoms=|--mna=8");
    usage_options.push_back("--cv_cutoff=|--cv_cut=0.05");
    usage_options.push_back("--conc_curve_range=|--conc_curve=0,1,1,0");
    usage_options.push_back("--conc_npts=20");
    usage_options.push_back("--temp_range=|--temp=300,5000");
    usage_options.push_back("--temp_npts=150");
    usage_options.push_back(" ");
    init::MessageOption("--usage", "QCA()", usage_options);
  }
}

#endif
