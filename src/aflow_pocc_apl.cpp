//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                  Marco Esters - Duke University 2021                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters, 2021.
//
// This file provides a framework to calculate phonon properties for
// disordered materials modeled using the POCC algorithm.

#include "aflow.h"
#include "aflow_pocc.h"
#include "APL/aflow_apl.h"

#define _DEBUG_POCC_APL_ false

using std::string;
using std::vector;
using std::deque;

static const string _POCC_APL_MODULE_ = "POCC-APL";
static const double _FREQ_WARNING_THRESHOLD_ = 1e-3;

namespace pocc {

  void POccCalculator::calculatePhononPropertiesAPL(const vector<double>& v_temperatures) {
    stringstream message;
    // Make sure that everything is consistent
    uint nruns = m_ARUN_directories.size();
    if (nruns == 0) {
      message << "Number of ARUN directories is zero.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if (l_supercell_sets.size() != nruns) {
      message << "Number of directories and number of supercells are different.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    vector<int> vexclude;
    // Store arun2skip so we can restore later
    string aruns2skip_backup = "";
    if (m_kflags.KBIN_POCC_EXCLUDE_UNSTABLE) {
      // vflag_control (command line) has priority
      if (XHOST.vflag_control.flag("ARUNS2SKIP")) {
        aruns2skip_backup = XHOST.vflag_control.getattachedscheme("ARUNS2SKIP");
      } else if (!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) {
        aruns2skip_backup = m_kflags.KBIN_POCC_ARUNS2SKIP_STRING;
      }
    }

    // Parse aflow.in options
    const vector<aurostd::xoption>& vxopts = m_kflags.KBIN_MODULE_OPTIONS.aplflags;
    aurostd::xoption aplopts;
    for (uint i = 0; i < vxopts.size(); i++) {
      const string& key = vxopts[i].keyword;
      aplopts.push_attached(key, vxopts[i].xscheme);
      aplopts.flag(key, vxopts[i].option);
    }
    apl::validateParametersDosAPL(aplopts, m_aflags, *p_FileMESSAGE);

    vector<apl::PhononCalculator> vphcalc = initializePhononCalculators();

    // Get phonon DOS from each run
    vector<xDOSCAR> vxdos = getPhononDoscars(vphcalc, aplopts, vexclude);
    if (vxdos.size() == 0) {
      if (m_kflags.KBIN_POCC_EXCLUDE_UNSTABLE) {
        message << "No structures left after excluding dynamically unstable representatives.";
        pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        return;
      } else {
        message << "No phonon DOSCARs calculated.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
    vphcalc.clear();

    // Compress new files in subdirectories
    for (uint i = 0; i < nruns; i++) {
      if (m_kflags.KZIP_COMPRESS) KBIN::CompressDirectory(m_ARUN_directories[i], m_kflags);
      aurostd::DirectoryChmod("777", m_ARUN_directories[i]);
    }

    // Store PHPOSCAR for plotting
    stringstream phposcar;
    xstructure xstr_phposcar = xstr_pocc;
    xstr_phposcar.is_vasp4_poscar_format = false;
    xstr_phposcar.is_vasp5_poscar_format = true;
    phposcar << xstr_phposcar;
    string phposcar_file = aurostd::CleanFileName(m_aflags.Directory + "/" + DEFAULT_APL_PHPOSCAR_FILE);
    aurostd::stringstream2file(phposcar, phposcar_file);
    if (!aurostd::FileExist(phposcar_file)) {
      message << "Could not write file " << phposcar_file << ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    uint nexclude = vexclude.size();
    if (nexclude > 0) {
      string aruns2skip = "";
      for (uint i = 0; i < nexclude; i++) {
        aruns2skip += m_ARUN_directories[vexclude[i]];
        if (i < nexclude - 1) aruns2skip += ",";
      }
      if (XHOST.vflag_control.flag("ARUNS2SKIP")) {
        aruns2skip = XHOST.vflag_control.getattachedscheme("ARUNS2SKIP") + "," + aruns2skip;
        XHOST.vflag_control.push_attached("ARUNS2SKIP", aruns2skip);
      } else if (!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) {
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING += "," + aruns2skip;
      } else { // Nothing set in aflow.in or command line, so add to kflags
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING = aruns2skip;
      }
      loadDataIntoCalculator();
      setDFTEnergies();
    }

    // Calculate vibrational properties
    stringstream aplout;

    apl::ThermalPropertiesCalculator tpc(*p_FileMESSAGE, *p_oss);
    double tpt_start = aurostd::string2utype<double>(aplopts.getattachedscheme("TSTART"));
    double tpt_end = aurostd::string2utype<double>(aplopts.getattachedscheme("TEND"));
    double tpt_step = aurostd::string2utype<double>(aplopts.getattachedscheme("TSTEP"));
    xDOSCAR xdos_T;
    double T = 0;
    for (uint t = 0; t < v_temperatures.size(); t++) {
      T = v_temperatures[t];
      xdos_T = getAveragePhononDos(T, vxdos);
      uint i = 0;
      for ( ; i < xdos_T.number_energies; i++) {
        if (xdos_T.venergy[i] > -_ZERO_TOL_) break;
      }
      if (i > 0) {
        double rel_idos = xdos_T.viDOS[0][i-1]/xdos_T.viDOS[0].back();
        if (rel_idos > _FREQ_WARNING_THRESHOLD_) {
          double percent = 100 * xdos_T.viDOS[0][i-1]/xdos_T.viDOS[0].back();
          message << "There are imaginary frequencies in the DOS at " << T << " K, covering "
            << std::dec << setprecision(3) << percent <<  "\% of the integrated DOS. These frequencies were"
            << " omitted in the calculation of thermodynamic properties.";
          pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }
      }

      string tstring = getTemperatureString(T);
      string filename = aurostd::CleanFileName(m_aflags.Directory + "/" + POCC_PHDOSCAR_FILE + "_T" + tstring + "K");
      message << "Writing out " << filename << ".";
      pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, *p_FileMESSAGE, *p_oss);
      stringstream phdoscar;
      phdoscar << xdos_T;
      aurostd::stringstream2file(phdoscar, filename);
      if (!aurostd::FileExist(filename)) {
        message << "Could not write file " << filename << ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }

      tpc.clear();
      tpc._directory = m_aflags.Directory;
      tpc.initialize(xdos_T, *p_FileMESSAGE, *p_oss);
      tpc.natoms = (uint) aurostd::round(aurostd::sum(xstr_pocc.comp_each_type));  // Set to number of atoms in parent structure
      tpc.calculateThermalProperties(tpt_start, tpt_end, tpt_step);
      aplout << "[POCC_APL_RESULTS]START_TEMPERATURE=" << tstring << "_K" << endl;
      tpc.addToAPLOut(aplout);
      aplout << "[POCC_APL_RESULTS]STOP_TEMPERATURE=" << tstring << "_K" << endl;
      tpc.writePropertiesToFile(m_aflags.Directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_FILE);
    }

    if (!aplout.str().empty()) {
      string filename = aurostd::CleanFileName(m_aflags.Directory + "/" + POCC_FILE_PREFIX + POCC_APL_OUT_FILE);
      aurostd::stringstream2file(aplout, filename);
      if (!aurostd::FileExist(filename)) {
        message << "Cannot open output file " << filename << ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
    }

    // Restore ARUNS2SKIP
    if (nexclude > 0) {
      if (XHOST.vflag_control.flag("ARUNS2SKIP")) {
        XHOST.vflag_control.push_attached("ARUNS2SKIP", aruns2skip_backup);
      } else if (!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) {
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING = aruns2skip_backup;
      } else {
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING = "";
      }
      loadDataIntoCalculator();
      setDFTEnergies();
    }
  }

  vector<apl::PhononCalculator> POccCalculator::initializePhononCalculators() {
    string message = "Initializing phonon calculators.";
    pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

    vector<apl::PhononCalculator> vphcalc;
    unsigned long long int isupercell = 0;
    int imax = 0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);
      imax = std::max((int) (*it).getHNFIndex(), imax);
      string directory = aurostd::CleanFileName(m_aflags.Directory + "/" + m_ARUN_directories[isupercell]);
      message = "Descending into directory " + directory + ".";
      pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

      string statefile = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_STATE_FILE);
      if (!aurostd::EFileExist(statefile)) {
        message = "Cannot find state file in directory " + directory + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }

      // Initialize phonon calculator
      apl::PhononCalculator phcalc(*p_FileMESSAGE, *p_oss);
      phcalc.initialize_supercell(statefile);
      phcalc.setDirectory(directory);
      phcalc.setNCPUs(m_kflags);
      phcalc._system = m_vflags.AFLOW_SYSTEM.content_string + "." + m_ARUN_directories[isupercell];
      phcalc.setPolarMaterial(aurostd::EFileExist(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_POLAR_FILE));

      // Get force constants
      string hibfile = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HARMIFC_FILE);
      bool awakeHarmIFCs = aurostd::EFileExist(hibfile);
      if (awakeHarmIFCs) {
        try {
          message = "Reading force constant from hibernate file.";
          pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, directory, *p_FileMESSAGE, *p_oss);
          phcalc.awake();
        } catch (aurostd::xerror& e) {
          message = e.buildMessageString() + " Recalculating force constants.";
          pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, directory, *p_FileMESSAGE, *p_oss);
          awakeHarmIFCs = false;
        }
      }
      if (!awakeHarmIFCs) {
        // Calculate force constants
        _xvasp xvasp;
        xvasp.Directory = directory;
        _aflags _aflowFlags = m_aflags;
        _aflowFlags.Directory = directory;
        _xinput _xInput(xvasp);
        _xflags _xFlags(m_vflags);
        string AflowIn = "";  // Dummy for AIMS, which is not supported by POCC
        apl::ForceConstantCalculator fccalc(phcalc.getSupercell(), *p_FileMESSAGE, *p_oss);

        // Read parameters and calculate IFCs
        fccalc.readFromStateFile(statefile);
        if (fccalc.run()) {
          fccalc.hibernate();
          phcalc.setHarmonicForceConstants(fccalc);
        } else {
          message = "Could not calculate harmonic force constants for directory " + directory + ".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
      }

      // Convert to primitive cell for projections
      phcalc.getSupercell().projectToPrimitive();
      vphcalc.push_back(phcalc);
    }
    return vphcalc;
  }

  vector<xDOSCAR> POccCalculator::getPhononDoscars(vector<apl::PhononCalculator>& vphcalc, xoption& aplopts, vector<int>& vexclude) {
    string message = "Calculating phonon densities of states.";
    pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

    bool projected = aplopts.flag("DOS_PROJECT");

    // q-point mesh
    vector<int> qpt_mesh;
    aurostd::string2tokens(aplopts.getattachedscheme("DOSMESH"), qpt_mesh, " xX");

    // Calculate phonon DOS
    uint nruns = vphcalc.size();
    vector<apl::DOSCalculator> vphdos(nruns);
    double minfreq = 0.0, maxfreq = 0.0;
    unsigned long long int isupercell = 0;
    vector<uint> vcalc;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);

      string directory = aurostd::CleanFileName(m_aflags.Directory + "/" + m_ARUN_directories[isupercell]);
      vphcalc[isupercell].initialize_qmesh(qpt_mesh, true, true);
      if (!projected) vphcalc[isupercell].getQMesh().makeIrreducible();

      apl::DOSCalculator dosc(vphcalc[isupercell], aplopts);
      if (m_kflags.KBIN_POCC_EXCLUDE_UNSTABLE && dosc.hasImaginaryFrequencies()) {
        vexclude.push_back(isupercell);
      } else {
        if (dosc.getMinFreq() < minfreq) minfreq = dosc.getMinFreq();
        if (dosc.getMaxFreq() > maxfreq) maxfreq = dosc.getMaxFreq();
        vphdos[isupercell] = dosc;
        vcalc.push_back((uint) isupercell);
      }
    }

    nruns = vcalc.size();
    vector<xDOSCAR> vxdos(nruns);

    aplopts.push_attached("MINFREQ", aurostd::utype2string<double>(minfreq));
    aplopts.push_attached("MAXFREQ", aurostd::utype2string<double>(maxfreq));

    if (aplopts.getattachedscheme("DOSMETHOD") == "LT") message = "Calculating phonon DOS using the linear tetrahedron method.";
    else message = "Calculating phonon DOS using the root sampling method.";
    pflow::logger(__AFLOW_FILE__, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex m;
    xthread::xThread xt(KBIN::get_NCPUS(m_kflags), 1);
    std::function<void(uint, const vector<uint>&, const aurostd::xoption&, vector<apl::DOSCalculator>&,
        vector<xDOSCAR>&, std::mutex&)> fn = std::bind(&POccCalculator::calculatePhononDOSThread, this,
          std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
          std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
    xt.run(nruns, fn, vcalc, aplopts, vphdos, vxdos, m);
#else
    for (uint i = 0; i < nruns; i++) calculatePhononDOSThread(i, vcalc, aplopts, vphdos, vxdos);
#endif

    return vxdos;
  }

#ifdef AFLOW_MULTITHREADS_ENABLE
  void POccCalculator::calculatePhononDOSThread(uint i, const vector<uint>& vcalc,
      const aurostd::xoption& aplopts, vector<apl::DOSCalculator>& vphdos, vector<xDOSCAR>& vxdos, std::mutex& m)
#else
  void POccCalculator::calculatePhononDOSThread(uint i, const vector<uint>& vcalc,
      const aurostd::xoption& aplopts, vector<apl::DOSCalculator>& vphdos, vector<xDOSCAR>& vxdos)
#endif
  {

    // Normalize to number of branches in the parent structure
    double pocc_sum = aurostd::sum(xstr_pocc.comp_each_type);  // Will be needed for projections
    double nbranches = 3 * pocc_sum;

    // DOS options
    int dos_npoints = aurostd::string2utype<int>(aplopts.getattachedscheme("DOSPOINTS"));
    double dos_smear = aurostd::string2utype<double>(aplopts.getattachedscheme("DOSSMEAR"));
    double minfreq = aurostd::string2utype<double>(aplopts.getattachedscheme("MINFREQ"));
    double maxfreq = aurostd::string2utype<double>(aplopts.getattachedscheme("MAXFREQ"));

    uint icalc = vcalc[i];
    vphdos[icalc].calc(dos_npoints, dos_smear, minfreq, maxfreq, false);
    xDOSCAR phdos = vphdos[icalc].createDOSCAR();

    // Normalize
    double norm_factor = nbranches/((double) vphdos[icalc].getNumberOfBranches());
    for (uint e = 0; e < phdos.number_energies; e++) phdos.viDOS[0][e] *= norm_factor;

    const xstructure& xstr_phcalc = vphdos[icalc].getInputStructure();
    uint natoms_phcalc = xstr_phcalc.atoms.size();
    uint nproj = phdos.vDOS[0].size();
    bool projected = (nproj > 1);
    for (uint at = 0; at < (projected?(natoms_phcalc + 1):1); at++) {  // +1 for totals
      for (uint p = 0; p < nproj; p++) {
        for (uint e = 0; e < phdos.number_energies; e++) {
          phdos.vDOS[at][p][0][e] *= norm_factor;
        }
      }
    }

    uint natoms_pocc = xstr_pocc.atoms.size();
    phdos.number_atoms = natoms_pocc;

    // If there are projected DOS, add the DOS contributions that belong to
    // each site in the parent structure.
    if (projected) {
      vector<uint> map = getMapToPARTCAR((unsigned long long int) icalc, xstr_phcalc);

      deque<deque<deque<deque<double> > > > vDOS_pocc(natoms_pocc + 1, deque<deque<deque<double> > >(nproj, deque<deque<double> >(1, deque<double>(phdos.number_energies, 0.0))));
      vDOS_pocc[0] = phdos.vDOS[0];  // Totals stay the same

      uint at_mapped = 0;
      for (uint at = 0; at < natoms_phcalc; at++) {
        at_mapped = map[at];
        for (uint p = 0; p < nproj; p++) {
          for (uint e = 0; e < phdos.number_energies; e++) {
            vDOS_pocc[at_mapped + 1][p][0][e] += phdos.vDOS[at + 1][p][0][e];
          }
        }
      }
      phdos.vDOS = vDOS_pocc;
    }
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::lock_guard<std::mutex> lk(m);
#endif
    vxdos[i] = std::move(phdos);
  }

  xDOSCAR POccCalculator::getAveragePhononDos(double T, const vector<xDOSCAR>& vxdos) {
    stringstream message;
    xDOSCAR xdos;

    setPOccStructureProbabilities(T);
    unsigned long long int isupercell = 0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);
      if (isupercell == 0) { // reset xDOSCAR
        xdos = vxdos[isupercell];
        std::fill(xdos.viDOS[0].begin(), xdos.viDOS[0].end(), 0.0);
        for (uint at = 0; at < xdos.vDOS.size(); at++) {
          for (uint p = 0; p < xdos.vDOS[at].size(); p++) {
            std::fill(xdos.vDOS[at][p][0].begin(), xdos.vDOS[at][p][0].end(), 0.0);
          }
        }
        xdos.temperature = T;
        xdos.Efermi = 0.0;
        xdos.spin = 0;
        xdos.spinF = 0.0;
        if (m_vflags.AFLOW_SYSTEM.content_string.empty()) {
          xdos.title = "Unknown system";
        } else {
          xdos.title = m_vflags.AFLOW_SYSTEM.content_string;
        }
      } else { // check that all dimensions are the same
        bool mismatch = false;
        // Check frequencies
        if (!mismatch) {
          if ((xdos.energy_max != vxdos[isupercell].energy_max)
              || (xdos.energy_min != vxdos[isupercell].energy_min)
              || (xdos.number_energies != vxdos[isupercell].venergy.size())) {
            message << "Frequencies do not match for phonon DOS in " << m_ARUN_directories[isupercell] << ".";
            mismatch = true;
          }
        }
        if (!mismatch) {
          if (xdos.vDOS.size() != vxdos[isupercell].vDOS.size()) {
            message << "Phonon DOS in " << m_ARUN_directories[isupercell] << " has different number of atoms.";
            mismatch = true;
          }
        }
        if (!mismatch) {
          for (uint i = 0; i < xdos.vDOS.size() && !mismatch; i++) {
            if (xdos.vDOS[i].size() != vxdos[isupercell].vDOS[i].size()) {
              message << "Phonon DOS in " << m_ARUN_directories[isupercell] << " has different number of projections.";
              mismatch = true;
            }
          }
        }
        if (mismatch) throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      uint nat = xdos.vDOS.size();
      uint nproj = xdos.vDOS[0].size();
      for (uint e = 0; e < xdos.number_energies; e++) {
        xdos.viDOS[0][e] += (*it).m_probability * vxdos[isupercell].viDOS[0][e];
        for (uint at = 0; at < nat; at++) {
          for (uint p = 0; p < nproj; p++) {
            xdos.vDOS[at][p][0][e] += (*it).m_probability * vxdos[isupercell].vDOS[at][p][0][e];
          }
        }
      }
    }
    return xdos;
  }

  bool POccCalculator::inputFilesFoundAnywhereAPL() {
    if (m_ARUN_directories.size() == 0) loadDataIntoCalculator();
    string directory = "";
    for (uint i = 0; i < m_ARUN_directories.size(); i++) {
      directory = m_aflags.Directory + "/" + m_ARUN_directories[i];
      if (aurostd::EFileExist(directory + "/" + DEFAULT_APL_PHPOSCAR_FILE)
          || aurostd::EFileExist(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_STATE_FILE)) {
        return true;
      }
    }
    return false;
  }

}  // namespace pocc

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                  Marco Esters - Duke University 2020                    *
// *                                                                         *
//****************************************************************************
