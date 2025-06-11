// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_DEBYE_CPP
#define _AFLOW_AGL_DEBYE_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"
#include "aflow_ael_elasticity.h"

// ###############################################################################
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2021)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. Materials 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//
// ***************************************************************************
//
// Constructor for _AGL_data class, which is a container for variables required for the AGL method
// Initializes values of variables
// User controlled variables are passed to the constructor and initialized to user selected values or default values as appropriate
// Other variables are initialized to zero
//
_AGL_data::_AGL_data() {
  i_eqn_of_state = 0;
  i_optimize_beta = 0;
  i_debye = 0;
  poissonratio = 0.0;
  energy_infinity = 0.0;
  maxloops = 0;
  maxfit = 0;
  maxpolycoeffs = 0;
  birchfitorder_iG = 0;
  fittype = 0;
  poissonratiosource="";
  precaccalgonorm = false;
  relax_static = false;
  static_only = false;
  relax_only = false;
  hugoniotrun = false;
  hugoniotextrapolate = false;
  ael_pressure_calc = false;
  postprocess = false;
  natoms = 0.0;
  cellmass = 0.0;
  dirpathname = "";
  sysname = "";
  pressure_external.clear();
  temperature_external.clear();
  ntemperature = 0;
  stemperature = 0.0;
  npressure = 0;
  spressure = 0.0;
  energyinput.clear();
  volumeinput.clear();
  energyinput_orig.clear();
  volumeinput_orig.clear();
  pressurecalculated.clear();
  stresscalculated.clear();
  structurecalculated.clear();
  pressurecalculatedmax = 0.0;
  tdebye.clear();
  gaussxm_debug = false;
  // Data calculated within AGL
  d2EnergydVolume2_static.clear();
  poissonratiofunction = 0.0;
  // Record of noise in E-V data
  EV_noise = false;
  // Record difference between minimum of E-V data and center
  itdiff = 0;
  // Equation of state data
  voleqmin.clear();
  bulkmodulus.clear();
  alpha.clear();
  d2EnergydVolume2_dynamic.clear();
  gamma_G.clear();
  gamma_poly.clear();
  pressure_static.clear();
  rms = 0.0;
  bulkmodulus_0pressure = 0.0;
  dbulkmodulusdpV_0pressure = 0.0;
  d2bulkmodulusdpV2_0pressure = 0.0; 
  // Saved data
  bcnt_beta = 0.0;
  x_K_opt = 0.0;
  x_m_opt = 0.0;
  bcntvolume0pressure = 0.0;
  optimize_beta = false;
  pfit.clear();
  IntEnergStatic.clear();
  bcnt_beta_statcalc = 0.0;
  x_K_opt_statcalc = 0.0;
  x_m_opt_statcalc = 0.0;
  volumestatcalc_0pressure = 0.0;
  gibbsenergystatcalc_0pressure = 0.0;
  x_Press_sp_statcalc = 0.0;
  bulkmodulusstatcalc_0pressure = 0.0;
  Vol_sp_statcalc = 0.0;
  astatic.clear();
  Avinetstatcalc_0pressure = 0.0;
  xsup_K_final = 0.0;
  Press_sp_final = 0.0;
  Vol_sp_final = 0.0;
  bcnt_beta_final = 0.0;
  // Highest temperature reached 
  max_temperature = 0.0;
  // zero pressure output data
  InternalEnergy0pressurekjmol.clear();
  InternalEnergy0pressuremeV.clear();
  Entropy0pressurekjmol.clear(); 
  Entropy0pressuremeV.clear();
  Entropy0pressureunitkB.clear();
  Cvkjmol0pressure.clear(); 
  Cpkjmol0pressure.clear();
  CvunitkB0pressure.clear(); 
  CpunitkB0pressure.clear(); 
  DebyeTemperature0pressure.clear(); 
  GruneisenParameter0pressure.clear();
  HelmholtzEnergy0pressurekjmol.clear();
  HelmholtzEnergy0pressuremeV.clear();
  GibbsFreeEnergy0pressurekjmol.clear();
  GibbsFreeEnergy0pressureeV.clear();
  ThermalExpansion0pressure.clear(); 
  bulkmodulusstatic_0pressure.clear();
  bulkmodulusisothermal_0pressure.clear();
  // Thermodynamic properties as a function of temperature and pressure
  AGL_pressure_temperature_energy_list.clear();
  AGL_pressure_enthalpy_list.clear();
  // Avoid truncating pressure or temperature range if no minimum energy is found 
  run_all_pressure_temperature = false;
  // Electronic properties as a function of pressure
  AGL_edos_properties.clear();
  edosgap_pressure = 0.0;
  dosvalEF_pressure = 0.0;
  egap_min = 0.0;
  egap_max = 0.0;
  egap_min_pressure = 0.0;
  egap_max_pressure = 0.0;
  dosvalEF_min = 0.0;
  dosvalEF_max = 0.0;
  dosvalEF_min_pressure = 0.0;
  dosvalEF_max_pressure = 0.0;
  // output data for all pressure values (optional)
  savedatapressure = false;
  InternalEnergyPressurekjmol.clear();
  EntropyPressurekjmol.clear();
  CvkjmolPressure.clear();
  DebyeTemperaturePressure.clear();
  GruneisenParameterPressure.clear();
  HelmholtzEnergyPressurekjmol.clear();
  InternalEnergyPressuremeV.clear();
  EntropyPressuremeV.clear();
  HelmholtzEnergyPressuremeV.clear();
  CvunitkBpressure.clear();
  EntropyunitkBpressure.clear();
  GibbsFreeEnergyPressureeV.clear();
  EnergyDFT_UIntVib.clear();
  // [OBSOLETE] EnergyDFT_UIntVibeV.clear();
  EnthalpyPressureeV.clear();
  xminsav.clear();
  VolumeEquilibrium.clear();
  mass_density_gcm3.clear();
  VolumeStaticPressure.clear();
  StaticPressure.clear();
  VolumeFactors.clear();
  failed_arun_list.clear();
  autoskipfailedaruns = false;
  skiparunsmax = 7;
  nstructsinit = 28;
}

//
// Destructor for AGL _AGL_data class
//
_AGL_data::~_AGL_data() {
  free();
}

void _AGL_data::free() {
}

//
// Copy constructor for AGL _AGL_data class
//
//_AGL_data::_AGL_data(const _AGL_data & b)
const _AGL_data& _AGL_data::operator=(const _AGL_data& b) {       // operator=
  if(this != &b) {
    i_eqn_of_state = b.i_eqn_of_state;
    i_optimize_beta = b.i_optimize_beta;
    i_debye = b.i_debye;
    poissonratio = b.poissonratio;
    energy_infinity = b.energy_infinity;
    maxloops = b.maxloops;
    maxfit = b.maxfit;
    maxpolycoeffs = b.maxpolycoeffs;
    birchfitorder_iG = b.birchfitorder_iG;
    fittype = b.fittype;
    poissonratiosource = b.poissonratiosource;
    precaccalgonorm = b.precaccalgonorm;
    relax_static = b.relax_static;
    static_only = b.static_only;
    relax_only = b.relax_only;
    hugoniotrun = b.hugoniotrun;
    hugoniotextrapolate = b.hugoniotextrapolate;
    ael_pressure_calc = b.ael_pressure_calc;
    postprocess = b.postprocess;
    natoms = b.natoms;
    cellmass = b.cellmass;
    dirpathname = b.dirpathname;
    sysname = b.sysname;
    pressure_external = b.pressure_external;
    temperature_external = b.temperature_external;
    ntemperature = b.ntemperature;
    stemperature = b.stemperature;
    npressure = b.npressure;
    spressure = b.spressure;
    energyinput = b.energyinput;
    volumeinput = b.volumeinput;
    energyinput_orig = b.energyinput_orig;
    volumeinput_orig = b.volumeinput_orig;
    pressurecalculated = b.pressurecalculated;
    stresscalculated = b.stresscalculated;
    structurecalculated = b.structurecalculated;    
    pressurecalculatedmax = b.pressurecalculatedmax;
    tdebye = b.tdebye;
    gaussxm_debug = b.gaussxm_debug;
    // Record of noise in E-V data
    EV_noise = b.EV_noise;
    // Record difference between minimum of E-V data and center
    itdiff = b.itdiff;
    // Data calculated within GIBBS
    d2EnergydVolume2_static = b.d2EnergydVolume2_static;
    poissonratiofunction = b.poissonratiofunction;
    // Equation of state data
    voleqmin = b.voleqmin;
    bulkmodulus = b.bulkmodulus;
    alpha = b.alpha;
    d2EnergydVolume2_dynamic = b.d2EnergydVolume2_dynamic;
    gamma_G = b.gamma_G;
    gamma_poly = b.gamma_poly;
    pressure_static = b.pressure_static;
    rms = b.rms;
    bulkmodulus_0pressure = b.bulkmodulus_0pressure;
    dbulkmodulusdpV_0pressure = b.dbulkmodulusdpV_0pressure;
    d2bulkmodulusdpV2_0pressure = b.d2bulkmodulusdpV2_0pressure; 
    // Saved data
    bcnt_beta = b.bcnt_beta;
    x_K_opt = b.x_K_opt;
    x_m_opt = b.x_m_opt;
    bcntvolume0pressure = b.bcntvolume0pressure;
    optimize_beta = b.optimize_beta;
    pfit = b.pfit;
    IntEnergStatic = b.IntEnergStatic;
    bcnt_beta_statcalc = b.bcnt_beta_statcalc;
    x_K_opt_statcalc = b.x_K_opt_statcalc;
    x_m_opt_statcalc = b.x_m_opt_statcalc;
    volumestatcalc_0pressure = b.volumestatcalc_0pressure;
    gibbsenergystatcalc_0pressure = b.gibbsenergystatcalc_0pressure;
    x_Press_sp_statcalc = b.x_Press_sp_statcalc;
    bulkmodulusstatcalc_0pressure = b.bulkmodulusstatcalc_0pressure;
    Vol_sp_statcalc = b.Vol_sp_statcalc;
    astatic = b.astatic;
    Avinetstatcalc_0pressure = b.Avinetstatcalc_0pressure;
    xsup_K_final = b.xsup_K_final;
    Press_sp_final = b.Press_sp_final;
    Vol_sp_final = b.Vol_sp_final;
    bcnt_beta_final = b.bcnt_beta_final;
    // Highest temperature reached 
    max_temperature = b.max_temperature;
    // zero pressure output data
    InternalEnergy0pressurekjmol = b.InternalEnergy0pressurekjmol;
    InternalEnergy0pressuremeV = b.InternalEnergy0pressuremeV;
    Entropy0pressurekjmol = b.Entropy0pressurekjmol;
    Entropy0pressuremeV = b.Entropy0pressuremeV;
    Entropy0pressureunitkB = b.Entropy0pressureunitkB;
    Cvkjmol0pressure = b.Cvkjmol0pressure;
    Cpkjmol0pressure = b.Cpkjmol0pressure;
    CvunitkB0pressure = b.CvunitkB0pressure; 
    CpunitkB0pressure = b.CpunitkB0pressure; 
    DebyeTemperature0pressure = b.DebyeTemperature0pressure; 
    GruneisenParameter0pressure = b.GruneisenParameter0pressure;
    HelmholtzEnergy0pressurekjmol = b.HelmholtzEnergy0pressurekjmol;
    HelmholtzEnergy0pressuremeV = b.HelmholtzEnergy0pressuremeV;
    GibbsFreeEnergy0pressurekjmol = b.GibbsFreeEnergy0pressurekjmol;
    GibbsFreeEnergy0pressureeV = b.GibbsFreeEnergy0pressureeV;
    ThermalExpansion0pressure = b.ThermalExpansion0pressure; 
    bulkmodulusstatic_0pressure = b.bulkmodulusstatic_0pressure;
    bulkmodulusisothermal_0pressure = b.bulkmodulusisothermal_0pressure;
    // Thermodynamic properties as a function of temperature and pressure
    AGL_pressure_temperature_energy_list = b.AGL_pressure_temperature_energy_list;
    AGL_pressure_enthalpy_list = b.AGL_pressure_enthalpy_list;
    //  Avoid truncating pressure or temperature range if no minimum energy is found 
    run_all_pressure_temperature = b.run_all_pressure_temperature;
    // Electronic properties as a function of pressure
    AGL_edos_properties = b.AGL_edos_properties;
    edosgap_pressure = b.edosgap_pressure;
    dosvalEF_pressure = b.dosvalEF_pressure;
    egap_min = b.egap_min;
    egap_max = b.egap_max;
    egap_min_pressure = b.egap_min_pressure;
    egap_max_pressure = b.egap_max_pressure;
    dosvalEF_min = b.dosvalEF_min;
    dosvalEF_max = b.dosvalEF_max;
    dosvalEF_min_pressure = b.dosvalEF_min_pressure;
    dosvalEF_max_pressure = b.dosvalEF_max_pressure;
    // output data for all pressure values (optional)
    savedatapressure = b.savedatapressure;
    InternalEnergyPressurekjmol = b.InternalEnergyPressurekjmol;
    EntropyPressurekjmol = b.EntropyPressurekjmol;
    CvkjmolPressure = b.CvkjmolPressure;
    DebyeTemperaturePressure = b.DebyeTemperaturePressure;
    GruneisenParameterPressure = b.GruneisenParameterPressure;
    HelmholtzEnergyPressurekjmol = b.HelmholtzEnergyPressurekjmol;
    InternalEnergyPressuremeV = b.InternalEnergyPressuremeV;
    EntropyPressuremeV = b.EntropyPressuremeV;
    HelmholtzEnergyPressuremeV = b.HelmholtzEnergyPressuremeV;
    CvunitkBpressure = b.CvunitkBpressure;
    EntropyunitkBpressure = b.EntropyunitkBpressure;
    GibbsFreeEnergyPressureeV = b.GibbsFreeEnergyPressureeV;
    EnergyDFT_UIntVib = b.EnergyDFT_UIntVib;
    // [OBSOLETE] EnergyDFT_UIntVibeV = b.EnergyDFT_UIntVibeV;
    EnthalpyPressureeV = b.EnthalpyPressureeV;
    xminsav = b.xminsav;
    VolumeEquilibrium = b.VolumeEquilibrium;
    mass_density_gcm3 = b.mass_density_gcm3;
    VolumeStaticPressure = b.VolumeStaticPressure;
    StaticPressure = b.StaticPressure;
    VolumeFactors = b.VolumeFactors;
    failed_arun_list = b.failed_arun_list;
    autoskipfailedaruns = b.autoskipfailedaruns;
    skiparunsmax = b.skiparunsmax;
    nstructsinit = b.nstructsinit;
  }
  return *this;
}

// ***************************************************************************
// KBIN::VASP_RunPhonons_AGL
// ***************************************************************************
namespace KBIN {
  //
  // Run AGL method: uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
  // See Computer Physics Communications 158, 57-72 (2004), Journal of Molecular Structure (Theochem) 368, 245-255 (1996), Phys. Rev. B 90, 174107 (2014) and Phys. Rev. Materials 1, 015401 (2017) for details
  // This function reads in user selections, creates set of strained structures to be run with VASP, reads VASP output and calls gibbsrun function to calculate thermodynamic properties
  //
  void VASP_RunPhonons_AGL(  _xvasp&  xvasp,
      string  AflowIn,
      _aflags& aflags,
      _kflags& kflags,
      _vflags& vflags, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    uint aglerror = 1;
    ostringstream aus;
    bool USER_RELAX=true;

    if (XHOST.GENERATE_AFLOWIN_ONLY && USER_RELAX) {
      USER_RELAX = false;
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "RELAX will be switched OFF for generate_aflowin_only." << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }

    if (USER_RELAX) {
      aglerror = relaxStructureAGL_VASP(AflowIn, xvasp, aflags, kflags, vflags, FileMESSAGE);
      if (aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failed to relax structure prior to running AGL." << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
        return;
      }
    }


    // Call RunDebye_AGL to run AGL
    aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if (aglerror == 9) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye failure to generate Hugoniot data" << endl;
      aus << _AGLSTR_MESSAGE_ << "Main AGL Debye completed successfully" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
  }
} // namespace KBIN

namespace KBIN {
  //CT20200323 - Relax structure if it hasn't relaxed yet before running AGL
  uint relaxStructureAGL_VASP(const string& AflowIn,
      _xvasp& xvasp,
      _aflags& aflags,
      _kflags& kflags,
      _vflags& vflags,
      ofstream& FileMESSAGE) {
    // bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    // uint num_relax = xvasp.NRELAX;
    uint num_relax = 2;
    bool relax_complete = true;
    bool krun = false;

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ + "Number of requested relaxation runs = " << num_relax << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // First check if relaxation has already been performed
    for (uint i = 1; i <= num_relax; i++) {
      string filename = aurostd::CleanFileName(xvasp.Directory + "/CONTCAR.relax" + aurostd::utype2string<int>(i));  //CO20200624 - patching --run vs --run=1 bug
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Relaxation CONTCAR filename = " << filename << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if (aurostd::FileExist(filename) || aurostd::FileExist(filename + ".xz")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Relaxation " << i << " has already completed"  << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Relaxation " << i << " has not completed"  << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        relax_complete = false;
      }
    }
    if (relax_complete) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Relaxation has already completed"  << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    } else {
      for (uint i = 1; i <= num_relax; i++) {
        krun = VASP_Produce_and_Modify_INPUT(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
        krun = (krun && VASP_Write_INPUT(xvasp, vflags));
        krun = VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string<int>(i), true, FileMESSAGE);
        XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp, aflags, vflags, i, (i<num_relax), FileMESSAGE);  //CO20210315 - write_incar only if (i<num_relax), no static afterwards
        if (i < num_relax) {
          XVASP_KPOINTS_IBZKPT_UPDATE(xvasp, aflags, vflags, i, true, FileMESSAGE); //CO20210315 - always true because (i < num_relax)
        }
      }
    }
    if (krun) {
      return 0;
    } else {
      return 1;
    }
  }
}   

// ***************************************************************************
// KBIN::VASP_RunPhonons_AGL
// ***************************************************************************
namespace KBIN {
  //
  // Run AGL postprocessing: calls AGL from other parts of AFLOW, to run AGL postprocessing on previously run calculations
  // See Computer Physics Communications 158, 57-72 (2004), Journal of Molecular Structure (Theochem) 368, 245-255 (1996), Phys. Rev. B 90, 174107 (2014) and Phys. Rev. Materials 1, 015401 (2017) for details
  //  void VASP_RunPhonons_AGL_postprocess(  _xvasp&  xvasp, string  AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE)
  void VASP_RunPhonons_AGL_postprocess(const string& directory_LIB, string& AflowInName, string& FileLockName) {  
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror = 0;
    ostringstream aus;
    _xvasp xvasp;
    string AflowIn = "";
    _aflags aflags;
    _kflags kflags;
    _vflags vflags;
    ofstream FileMESSAGE;

    // Call AGL_xvasp_flags_populate to populate xvasp, aflags, kflags and vflags classes
    uint aglerror = AGL_functions::AGL_xvasp_flags_populate(xvasp, AflowIn, AflowInName, FileLockName, directory_LIB, aflags, kflags, vflags, FileMESSAGE);
    if (aglerror != 0) { //CT20200624: Check for error, aglerror=2 implies not an AGL main directory
      if (aglerror == 2) {
        cerr << _AGLSTR_MESSAGE_ << "Not AGL main directory" << endl;
        cerr << _AGLSTR_MESSAGE_ << AflowInName << endl;
        return;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "AGL set xvasp flags failed" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return;
      }
    }

    // Set AGL postprocess flag to true
    AGL_data.postprocess = true;

    // Call RunDebye_AGL to run AGL
    aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
  }
} // namespace KBIN

// ***************************************************************************
// AGL_functions::Get_ThermalProperties_AGL_postprocess
// ***************************************************************************
namespace AGL_functions {
  //
  // Run AGL postprocessing: calls AGL from other parts of AFLOW, to run AGL postprocessing on previously run calculations
  // Returns thermal properties: Debye temperature, Gruneisen parameter, Gibbs free energy per atom, vibrational free energy per atom
  // See Computer Physics Communications 158, 57-72 (2004), Journal of Molecular Structure (Theochem) 368, 245-255 (1996), Phys. Rev. B 90, 174107 (2014) and Phys. Rev. Materials 1, 015401 (2017) for details
  uint Get_ThermalProperties_AGL_postprocess(const string& directory, uint ntemperature, double stemperature, uint npressure, double spressure, vector<double>& agl_temperature, vector<double>& agl_gibbs_energy_atom, vector<double>& agl_vibrational_energy_atom) {  
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror = 0;
    ostringstream aus;
    _xvasp xvasp;
    string AflowIn = "";
    _aflags aflags;
    _kflags kflags;
    _vflags vflags;
    ofstream FileMESSAGE;
    // string AflowInName = _AFLOWIN_AGL_DEFAULT_;
    // string FileLockName = "agl.LOCK";
    // string AflowInName = _AFLOWIN_;
    // string FileLockName = _AFLOWLOCK_;    
    string AflowInName = "", FileLockName = "";

    // string filedirname = aurostd::CleanFileName(directory) + "agl.LOCK";
    if (aurostd::FileExist(aurostd::CleanFileName(directory) + "agl.LOCK")) {
      FileLockName = "agl.LOCK";
    } else if (aurostd::FileExist(aurostd::CleanFileName(directory) + "ael.LOCK")) {
      FileLockName = "ael.LOCK";
    } else {
      FileLockName = _AFLOWLOCK_;
    }

    if (aurostd::FileExist(aurostd::CleanFileName(directory) + _AFLOWIN_AGL_DEFAULT_)) {
      AflowInName = _AFLOWIN_AGL_DEFAULT_;
    } else if (aurostd::FileExist(aurostd::CleanFileName(directory) + _AFLOWIN_AEL_DEFAULT_)) {
      AflowInName = _AFLOWIN_AEL_DEFAULT_;
    } else {
      AflowInName = _AFLOWIN_;
    }

    // Call AGL_xvasp_flags_populate to populate xvasp, aflags, kflags and vflags classes
    uint aglerror = AGL_functions::AGL_xvasp_flags_populate(xvasp, AflowIn, AflowInName, FileLockName, directory, aflags, kflags, vflags, FileMESSAGE);

    // Set AGL postprocess flag to true
    AGL_data.postprocess = true;

    // Set temperature and pressure number/size of steps
    AGL_data.ntemperature = ntemperature;
    AGL_data.stemperature = stemperature;
    AGL_data.npressure = npressure;
    AGL_data.spressure = spressure;

    // Call RunDebye_AGL to run AGL
    aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      for (uint i = 0; i < AGL_data.GibbsFreeEnergy0pressureeV.size(); i++) {
        agl_temperature.push_back(AGL_data.temperature_external.at(i));
        agl_gibbs_energy_atom.push_back(AGL_data.GibbsFreeEnergy0pressureeV.at(i) / AGL_data.natoms);
        agl_vibrational_energy_atom.push_back(AGL_data.HelmholtzEnergy0pressuremeV.at(i) / AGL_data.natoms);
      }
      return aglerror;
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aglerror;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aglerror;
    }
  }
} // namespace AGL_functions

// ***************************************************************************
//
// Set of functions to allow AGL to be called from other parts of AFLOW to obtain equilibrium volume, bulk modulus, etc. as a function of temperature
// See Computer Physics Communications 158, 57-72 (2004), Journal of Molecular Structure (Theochem) 368, 245-255 (1996), Phys. Rev. B 90, 174107 (2014) and Phys. Rev. Materials 1, 015401 (2017) for details
// Note that AGL will need to be run twice if the E(V) calculations are not already available - see the README file for details. 
// If AGL needs to be rerun after performing the E(V) calculations, these functions will return a value of 8. If AGL has run and the data is available, they will return 0.
// If there is an error, they will return an unsigned integer value other than 0 or 8.
//
// ***************************************************************************
// AGL_functions::Get_EquilibriumVolumeTemperature
// ***************************************************************************
namespace AGL_functions {
  uint Get_EquilibriumVolumeTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror;
    ostringstream aus;

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
        Temperature.push_back(AGL_data.temperature_external.at(i));
        EquilibriumVolume.push_back(AGL_data.VolumeEquilibrium.at(i).at(0));
      }
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::Get_EquilibriumVolumeAngstromTemperature
// ***************************************************************************
namespace AGL_functions {
  uint Get_EquilibriumVolumeAngstromTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror;
    ostringstream aus;

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
        Temperature.push_back(AGL_data.temperature_external.at(i));
        // [OBSOLETE] EquilibriumVolume.push_back((AGL_data.VolumeEquilibrium.at(i).at(0)/ pow(angstrom2bohr, 3.0)));
        EquilibriumVolume.push_back((AGL_data.VolumeEquilibrium.at(i).at(0)));
      }
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::Get_BulkModulusStaticTemperature
// ***************************************************************************
namespace AGL_functions {
  uint Get_BulkModulusStaticTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusStatic, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror;
    ostringstream aus;

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
        Temperature.push_back(AGL_data.temperature_external.at(i));
        BulkModulusStatic.push_back(AGL_data.bulkmodulusstatic_0pressure.at(i));
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::Get_BulkModulusIsothermalTemperature
// ***************************************************************************
namespace AGL_functions {
  uint Get_BulkModulusIsothermalTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusIsothermal, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror;
    ostringstream aus;

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
        Temperature.push_back(AGL_data.temperature_external.at(i));
        BulkModulusIsothermal.push_back(AGL_data.bulkmodulusisothermal_0pressure.at(i));
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::Get_BulkModulusVolumeTemperature
// ***************************************************************************
namespace AGL_functions {
  // Function to be specifically called from APL 2.0 which calculates Gruneisen parameter 
  // Runs AGL to obtain Equilibrium Volume in cubic Bohr and Isothermal Bulk Modulus in GPa as a function of temperature for use in calculating thermal expansion 
  // Used in APL 2.0 to get thermal expansion and heat capacity at constant pressure 
  uint Get_BulkModulusVolumeTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusIsothermal, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror = 0;
    ostringstream aus;

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
        Temperature.push_back(AGL_data.temperature_external.at(i));
        BulkModulusIsothermal.push_back(AGL_data.bulkmodulusisothermal_0pressure.at(i));
        EquilibriumVolume.push_back(AGL_data.VolumeEquilibrium.at(i).at(0));
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::Get_BulkModulusVolumeAngstromTemperature
// ***************************************************************************
namespace AGL_functions {
  // Function to be specifically called from APL 2.0 which calculates Gruneisen parameter 
  // Runs AGL to obtain Equilibrium Volume in cubic Angstrom and Isothermal Bulk Modulus in GPa as a function of temperature for use in calculating thermal expansion 
  // Used in APL 2.0 to get thermal expansion and heat capacity at constant pressure 
  uint Get_BulkModulusVolumeAngstromTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusIsothermal, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror = 0;
    ostringstream aus;

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      for (uint i = 0; i < AGL_data.temperature_external.size(); i++) {
        Temperature.push_back(AGL_data.temperature_external.at(i));
        BulkModulusIsothermal.push_back(AGL_data.bulkmodulusisothermal_0pressure.at(i));
        // [OBSOLETE] EquilibriumVolume.push_back((AGL_data.VolumeEquilibrium.at(i).at(0)/ pow(angstrom2bohr, 3.0)));
        EquilibriumVolume.push_back((AGL_data.VolumeEquilibrium.at(i).at(0)));
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::Get_VolumeStaticPressure
// ***************************************************************************
namespace AGL_functions {
  // Function to be called from APL or AEL to calculate thermal or mechanical properties as a function of pressure
  // Runs AGL to obtain Volume in cubic Angstrom and Volume scaling factors as a function of static pressure in GPa 
  // Required pressure values should be specified in the aflow.in file
  uint Get_VolumeStaticPressure(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Pressure, vector<double>& PressureVolumes, vector<double>& VolumeScaleFactors, bool postprocess, ofstream& FileMESSAGE) {
    // Class to contain AGL input and output data
    _AGL_data AGL_data;
    // [OBSOLETE] uint aglerror = 0;
    ostringstream aus;

    // Set ael_pressure_calc to true
    AGL_data.ael_pressure_calc = true;

    // Set AGL postprocess flag
    AGL_data.postprocess = postprocess;    

    // Call RunDebye_AGL to run AGL
    uint aglerror = AGL_functions::RunDebye_AGL(xvasp, AflowIn, aflags, kflags, vflags, AGL_data, FileMESSAGE);
    if(aglerror == 0) {
      if(AGL_data.VolumeStaticPressure.size() != AGL_data.StaticPressure.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Mismatch in pressure-volume vectors sizes" << endl;  
        aus << _AGLSTR_ERROR_ + "AGL_data.StaticPressure.size() = " << AGL_data.StaticPressure.size() << endl;  
        aus << _AGLSTR_ERROR_ + "AGL_data.VolumeStaticPressure.size() = " << AGL_data.VolumeStaticPressure.size() << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = 3;
        return aglerror;
      } else if(AGL_data.VolumeStaticPressure.size() != AGL_data.VolumeFactors.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Mismatch in pressure-volume vectors sizes" << endl;  
        aus << _AGLSTR_ERROR_ + "AGL_data.StaticPressure.size() = " << AGL_data.StaticPressure.size() << endl;  
        aus << _AGLSTR_ERROR_ + "AGL_data.VolumeFactors.size() = " << AGL_data.VolumeFactors.size() << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = 3;
        return aglerror;
      } 
      for (uint i = 0; i < AGL_data.StaticPressure.size(); i++) {
        Pressure.push_back(AGL_data.StaticPressure.at(i));
        PressureVolumes.push_back(AGL_data.VolumeStaticPressure.at(i));
        VolumeScaleFactors.push_back(AGL_data.VolumeFactors.at(i));
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aglerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AGL Debye run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "AGL Debye run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::RunDebye_AGL
// ***************************************************************************
namespace AGL_functions {
  // Function to actually run the AGL method
  uint RunDebye_AGL(_xvasp& xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    // [OBSOLETE] bool LDEBUG=(FALSE || XHOST.DEBUG);
    // User's control of parameters for GIBBS calculation; setting defaults

    // Options to decide what to calculate and write
    aurostd::xoption USER_WRITE_GIBBS_INPUT;             USER_WRITE_GIBBS_INPUT.option = false;  
    aurostd::xoption USER_WRITE_FULL_RESULTS;            USER_WRITE_FULL_RESULTS.option = false;
    aurostd::xoption USER_SYSOUT_FILENAME;               USER_SYSOUT_FILENAME.option = false;
    aurostd::xoption USER_SET_PRESSURE_VALUE;            USER_SET_PRESSURE_VALUE.option = false;
    aurostd::xoption USER_SAVE_ALL_PRESSURES;            USER_SAVE_ALL_PRESSURES.option = false;
    aurostd::xoption USER_WRITE_ALL_PRESSURES;           USER_WRITE_ALL_PRESSURES.option = false;
    aurostd::xoption USER_CHECK_EV_CONCAVITY;            USER_CHECK_EV_CONCAVITY.option = false;  
    aurostd::xoption USER_CHECK_EV_MIN;                  USER_CHECK_EV_MIN.option = false;
    aurostd::xoption USER_THETA_GAMMA_KAPPA;             USER_THETA_GAMMA_KAPPA.option = false;  
    aurostd::xoption USER_GRUNEISEN;                     USER_GRUNEISEN.option = false;
    aurostd::xoption USER_THETA_COND;                    USER_THETA_COND.option = false;
    aurostd::xoption USER_PLOT_RESULTS;                  USER_PLOT_RESULTS.option = false;
    aurostd::xoption USER_KAPPA_VOLUME;                  USER_KAPPA_VOLUME.option = false;
    aurostd::xoption USER_AEL_POISSON_RATIO;             USER_AEL_POISSON_RATIO.option = false;
    aurostd::xoption USER_SPECIES_MASS;                  USER_SPECIES_MASS.option = false;
    aurostd::xoption USER_DIRNAME_ARUN;                  USER_DIRNAME_ARUN.option = false;
    aurostd::xoption USER_GAUSSXM_DEBUG;                 USER_GAUSSXM_DEBUG.option = false;
    aurostd::xoption USER_PRECACC_ALGONORM;              USER_PRECACC_ALGONORM.option = false;
    aurostd::xoption USER_RELAX_STATIC;                  USER_RELAX_STATIC.option = false;
    aurostd::xoption USER_STATIC;                        USER_STATIC.option = false;
    aurostd::xoption USER_AUTOSKIP_FAILED_ARUNS;         USER_AUTOSKIP_FAILED_ARUNS.option = false;
    aurostd::xoption USER_WRITE_HUGONIOT_INPUT;          USER_WRITE_HUGONIOT_INPUT.option = false;
    aurostd::xoption USER_HUGONIOT_CALC;                 USER_HUGONIOT_CALC.option = true;
    aurostd::xoption USER_HUGONIOT_EXTRAPOLATE;          USER_HUGONIOT_EXTRAPOLATE.option = false;
    aurostd::xoption USER_RUN_ALL_PRESSURE_TEMPERATURE;  USER_RUN_ALL_PRESSURE_TEMPERATURE.option = false;

    // Parameters for determining calculation details
    string USER_SYSTEM_NAME                 = "AGL";
    int    USER_NSTRUCTURES                 = 28;
    double USER_STRAIN_STEP                 = 0.01;
    int    USER_IEOS                        = 0;
    int    USER_IOPTIMIZE_BETA              = 0;
    int    USER_IDEBYE                      = 0;
    double USER_POISSON_RATIO               = 0.25;
    uint   USER_NPRESSURE                   = 101;
    double USER_SPRESSURE                   = 1.0;
    uint   USER_NTEMPERATURE                = 201;
    double USER_STEMPERATURE                = 10.0;
    double USER_EINF                        = 0.0;
    int    USER_MAX_LOOPS                   = 250;
    int    USER_MAX_FITS                    = 500;
    int    USER_MAX_PAR                     = 50;
    int    USER_BIRCH_FIT_ORDER_IG          = 2;
    int    USER_MAX_CC_ITER                 = 1;  
    int    USER_MAX_CM_ITER                 = 5;  
    int    USER_FIT_TYPE                    = 0;
    string USER_SKIP_FAILED_ARUNS           = "";
    int    USER_SKIP_ARUNS_MAX              = 7;
    double USER_DEBYE_TEMP;
    double USER_PRES_VAL;
    double tdmin, tdmax, DEB_min, DEB_max, Tmin, Tmax, CV_min, CV_max, cvmin, cvmax, CP_min, CP_max, cpmin, cpmax, GA_min, GA_max, gamin, gamax, VF_min, VF_max, vfmin, vfmax;
    double Ecellmin, Ecellmax, Eatommin, Eatommax, Vcellmin, Vcellmax, Vatommin, Vatommax, ecmin, ecmax, eamin, eamax, vcmin, vcmax, vamin, vamax;
    double tdbest, cellmass_grams, energy_kJg, volume_cm3, massdensity_gcm3;
    ostringstream aus;
    stringstream oss;
    // Boltzmann's constant / number of atoms per cell pre-factor for heat capacity; set to 1.0 so heat capacity is in units of (kB / cell)
    double nkb = 1.0;
    // Error return
    uint aglerror = 0;
    // Tolerance for zero value; if a double is less than this it will be considered to be less than or equal to zero (which could generate an error)
    double tolzero = 1e-12;
    // Record position of first non-NaN value in vector (used to ensure only non-NaN results are written in JSON files)
    uint firstnonnan = 0;
    // Conversion from eV to Hartree
    // [OBSOLETE] double ev2hart1=1.0/27.211383;
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);

    aurostd::StringstreamClean(aus);
    aus << "00000  MESSAGE Starting GIBBS RUN" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Get user's parameters from _AFLOWIN_ 
    // Get user's values of what to calculate and write

    // Check's if the system name is present in the _AFLOWIN_ file
    // If it exists, it is used to write the output filenames
    if( aurostd::substring2bool(AflowIn,_AFSTROPT_+"SYSTEM=",TRUE) ) {
      USER_SYSTEM_NAME = aurostd::substring2string(AflowIn,_AFSTROPT_+"SYSTEM=",1,TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "System name = " << USER_SYSTEM_NAME.c_str() << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      

    // Get the user's selection of whether to print the input file for the original version of the GIBBS program (useful for comparison and debugging).
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITE_GIBBS_INPUT=",TRUE) ) {
      USER_WRITE_GIBBS_INPUT.options2entry(AflowIn,_AGLSTROPT_+"WRITE_GIBBS_INPUT=",USER_WRITE_GIBBS_INPUT.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write GIBBS input = " << USER_WRITE_GIBBS_INPUT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"WRITE_GIBBS_INPUT=",TRUE) ) {
      USER_WRITE_GIBBS_INPUT.options2entry(AflowIn,_AGIBBSTROPT_+"WRITE_GIBBS_INPUT=",USER_WRITE_GIBBS_INPUT.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write GIBBS input = " << USER_WRITE_GIBBS_INPUT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITEGIBBSINPUT=",TRUE) ) {
      USER_WRITE_GIBBS_INPUT.options2entry(AflowIn,_AGLSTROPT_+"WRITEGIBBSINPUT=",USER_WRITE_GIBBS_INPUT.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write GIBBS input = " << USER_WRITE_GIBBS_INPUT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"WRITEGIBBSINPUT=",TRUE) ) {
      USER_WRITE_GIBBS_INPUT.options2entry(AflowIn,_AGIBBSTROPT_+"WRITEGIBBSINPUT=",USER_WRITE_GIBBS_INPUT.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write GIBBS input = " << USER_WRITE_GIBBS_INPUT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }           
    // Get the user's selection of whether to write out the results for the Debye temperature, heat capacity, Gruneisen parameter, vibrational free energy and other thermal properties as a function of temperature. 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITE_FULL_RESULTS=",TRUE) ) {
      USER_WRITE_FULL_RESULTS.options2entry(AflowIn,_AGLSTROPT_+"WRITE_FULL_RESULTS=",USER_WRITE_FULL_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_WRITE_FULL_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"WRITE_FULL_RESULTS=",TRUE) ) {
      USER_WRITE_FULL_RESULTS.options2entry(AflowIn,_AGIBBSTROPT_+"WRITE_FULL_RESULTS=",USER_WRITE_FULL_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_WRITE_FULL_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITEFULLRESULTS=",TRUE) ) {
      USER_WRITE_FULL_RESULTS.options2entry(AflowIn,_AGLSTROPT_+"WRITEFULLRESULTS=",USER_WRITE_FULL_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_WRITE_FULL_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"WRITEFULLRESULTS=",TRUE) ) {
      USER_WRITE_FULL_RESULTS.options2entry(AflowIn,_AGIBBSTROPT_+"WRITEFULLRESULTS=",USER_WRITE_FULL_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_WRITE_FULL_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to use the system name to name the file containing the GIBBS output; default filename is AGL.out.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SYSOUT_FILENAME=",TRUE) ) {
      USER_SYSOUT_FILENAME.options2entry(AflowIn,_AGLSTROPT_+"SYSOUT_FILENAME=",USER_SYSOUT_FILENAME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_SYSOUT_FILENAME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SYSOUT_FILENAME=",TRUE) ) {
      USER_SYSOUT_FILENAME.options2entry(AflowIn,_AGIBBSTROPT_+"SYSOUT_FILENAME=",USER_SYSOUT_FILENAME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_SYSOUT_FILENAME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SYSOUTFILENAME=",TRUE) ) {
      USER_SYSOUT_FILENAME.options2entry(AflowIn,_AGLSTROPT_+"SYSOUTFILENAME=",USER_SYSOUT_FILENAME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_SYSOUT_FILENAME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SYSOUTFILENAME=",TRUE) ) {
      USER_SYSOUT_FILENAME.options2entry(AflowIn,_AGIBBSTROPT_+"SYSOUTFILENAME=",USER_SYSOUT_FILENAME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write full results = " << USER_SYSOUT_FILENAME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to check the concavity of the (E, V) data before performing the GIBBS method.
    // If the data is not concave, then it reruns the VASP calculations with extra k-points.
    // If this option is switched off, the data will be checked for concavity anyway within the GIBBS method, and, depending on the fitting method chosen, only the suitable points are chosen
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"CHECK_EV_CONCAVITY=",TRUE) ) {
      USER_CHECK_EV_CONCAVITY.options2entry(AflowIn,_AGLSTROPT_+"CHECK_EV_CONCAVITY=",USER_CHECK_EV_CONCAVITY.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) concavity = " << USER_CHECK_EV_CONCAVITY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"CHECK_EV_CONCAVITY=",TRUE) ) {
      USER_CHECK_EV_CONCAVITY.options2entry(AflowIn,_AGIBBSTROPT_+"CHECK_EV_CONCAVITY=",USER_CHECK_EV_CONCAVITY.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) concavity = " << USER_CHECK_EV_CONCAVITY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"CHECKEVCONCAVITY=",TRUE) ) {
      USER_CHECK_EV_CONCAVITY.options2entry(AflowIn,_AGLSTROPT_+"CHECKEVCONCAVITY=",USER_CHECK_EV_CONCAVITY.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) concavity = " << USER_CHECK_EV_CONCAVITY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"CHECKEVCONCAVITY=",TRUE) ) {
      USER_CHECK_EV_CONCAVITY.options2entry(AflowIn,_AGIBBSTROPT_+"CHECKEVCONCAVITY=",USER_CHECK_EV_CONCAVITY.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) concavity = " << USER_CHECK_EV_CONCAVITY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }        
    // Get the user's selection of whether to check the position of the energy minimum. 
    // If the minimum is for the maximum or minimum volumes, then it runs extra VASP calculations for smaller or larger structures until it brackets the minimum.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"CHECK_EV_MIN=",TRUE) ) {
      USER_CHECK_EV_MIN.options2entry(AflowIn,_AGLSTROPT_+"CHECK_EV_MIN=",USER_CHECK_EV_MIN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) minimum = " << USER_CHECK_EV_MIN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"CHECK_EV_MIN=",TRUE) ) {
      USER_CHECK_EV_MIN.options2entry(AflowIn,_AGIBBSTROPT_+"CHECK_EV_MIN=",USER_CHECK_EV_MIN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) minimum = " << USER_CHECK_EV_MIN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"CHECKMIN=",TRUE) ) {
      USER_CHECK_EV_MIN.options2entry(AflowIn,_AGLSTROPT_+"CHECKMIN=",USER_CHECK_EV_MIN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) minimum = " << USER_CHECK_EV_MIN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"CHECKMIN=",TRUE) ) {
      USER_CHECK_EV_MIN.options2entry(AflowIn,_AGIBBSTROPT_+"CHECKMIN=",USER_CHECK_EV_MIN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Check (E, V) minimum = " << USER_CHECK_EV_MIN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's choice for whether or not to explicitly specify the pressure values to be used
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SET_PRESSURE_VALUE=",TRUE) ) {
      USER_SET_PRESSURE_VALUE.options2entry(AflowIn,_AGLSTROPT_+"SET_PRESSURE_VALUE=",USER_SET_PRESSURE_VALUE.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "User defined pressure values option = " << USER_SET_PRESSURE_VALUE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SET_PRESSURE_VALUE=",TRUE) ) {
      USER_SET_PRESSURE_VALUE.options2entry(AflowIn,_AGIBBSTROPT_+"SET_PRESSURE_VALUE=",USER_SET_PRESSURE_VALUE.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "User defined pressure values option = " << USER_SET_PRESSURE_VALUE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SETPRESSUREVALUE=",TRUE) ) {
      USER_SET_PRESSURE_VALUE.options2entry(AflowIn,_AGLSTROPT_+"SETPRESSUREVALUE=",USER_SET_PRESSURE_VALUE.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "User defined pressure values option = " << USER_SET_PRESSURE_VALUE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SETPRESSUREVALUE=",TRUE) ) {
      USER_SET_PRESSURE_VALUE.options2entry(AflowIn,_AGIBBSTROPT_+"SETPRESSUREVALUE=",USER_SET_PRESSURE_VALUE.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "User defined pressure values option = " << USER_SET_PRESSURE_VALUE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }          
    // Get the user's selection of whether to save the values of thermal properties such as the Debye temperature, heat capacity and vibrational free energy for all temperatures and pressures
    // This option will require a lot of memory, but is useful for generating (p, T) phase diagrams 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SAVE_ALL_PRESSURES=",TRUE) ) {
      USER_SAVE_ALL_PRESSURES.options2entry(AflowIn,_AGLSTROPT_+"SAVE_ALL_PRESSURES=",USER_SAVE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Save results for all pressure values = " << USER_SAVE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SAVE_ALL_PRESSURES=",TRUE) ) {
      USER_SAVE_ALL_PRESSURES.options2entry(AflowIn,_AGIBBSTROPT_+"SAVE_ALL_PRESSURES=",USER_SAVE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Save results for all pressure values = " << USER_SAVE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SAVEALLPRESSURES=",TRUE) ) {
      USER_SAVE_ALL_PRESSURES.options2entry(AflowIn,_AGLSTROPT_+"SAVEALLPRESSURES=",USER_SAVE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Save results for all pressure values = " << USER_SAVE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SAVEALLPRESSURES=",TRUE) ) {
      USER_SAVE_ALL_PRESSURES.options2entry(AflowIn,_AGIBBSTROPT_+"SAVEALLPRESSURES=",USER_SAVE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Save results for all pressure values = " << USER_SAVE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }        	
    // Get the user's selection of whether to write the values of thermal properties such as the Debye temperature, heat capacity and vibrational free energy for all temperatures and pressures
    // This option will require a lot of memory, but is useful for generating (p, T) phase diagrams 
    // This option will only work if the SAVE_ALL_PRESSURES option is switched on
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITE_ALL_PRESSURES=",TRUE) ) {
      USER_WRITE_ALL_PRESSURES.options2entry(AflowIn,_AGLSTROPT_+"WRITE_ALL_PRESSURES=",USER_WRITE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write results for all pressure values = " << USER_WRITE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"WRITE_ALL_PRESSURES=",TRUE) ) {
      USER_WRITE_ALL_PRESSURES.options2entry(AflowIn,_AGIBBSTROPT_+"WRITE_ALL_PRESSURES=",USER_WRITE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write results for all pressure values = " << USER_WRITE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITEALLPRESSURES=",TRUE) ) {
      USER_WRITE_ALL_PRESSURES.options2entry(AflowIn,_AGLSTROPT_+"WRITEALLPRESSURES=",USER_WRITE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write results for all pressure values = " << USER_WRITE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"WRITEALLPRESSURES=",TRUE) ) {
      USER_WRITE_ALL_PRESSURES.options2entry(AflowIn,_AGIBBSTROPT_+"WRITEALLPRESSURES=",USER_WRITE_ALL_PRESSURES.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write results for all pressure values = " << USER_WRITE_ALL_PRESSURES.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }        	
    // Get the user's choice for whether or not to call AEL to calculate the Poisson ratio of the material
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"AEL_POISSON_RATIO=",TRUE) ) {
      USER_AEL_POISSON_RATIO.options2entry(AflowIn,_AGLSTROPT_+"AEL_POISSON_RATIO=",USER_AEL_POISSON_RATIO.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"AEL_POISSON_RATIO=",TRUE) ) {
      USER_AEL_POISSON_RATIO.options2entry(AflowIn,_AGIBBSTROPT_+"AEL_POISSON_RATIO=",USER_AEL_POISSON_RATIO.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"AELPOISSONRATIO=",TRUE) ) {
      USER_AEL_POISSON_RATIO.options2entry(AflowIn,_AGLSTROPT_+"AELPOISSONRATIO=",USER_AEL_POISSON_RATIO.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"AELPOISSONRATIO=",TRUE) ) {
      USER_AEL_POISSON_RATIO.options2entry(AflowIn,_AGIBBSTROPT_+"AELPOISSONRATIO=",USER_AEL_POISSON_RATIO.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }        	


    // Get user's values of parameters for AGL calculation

    // Get the user's number of strained structures to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"NSTRUCTURES=",TRUE) ) {
      USER_NSTRUCTURES = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"NSTRUCTURES=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Number of structures = " << USER_NSTRUCTURES << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"NSTRUCTURES=",TRUE) ) {
      USER_NSTRUCTURES = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"NSTRUCTURES=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Number of structures = " << USER_NSTRUCTURES << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_NSTRUCTURES < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Number of structures = " << USER_NSTRUCTURES << " < 1" << endl;  
      USER_NSTRUCTURES = 28;
      aus << _AGLSTR_WARNING_ + "Increasing number of structures to default value of " << USER_NSTRUCTURES << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's strain step size to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"STRAIN_STEP=",TRUE) ) {
      USER_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"STRAIN_STEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Strain step size = " << USER_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"STRAIN_STEP=",TRUE) ) {
      USER_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+"STRAIN_STEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Strain step size = " << USER_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"STRAINSTEP=",TRUE) ) {
      USER_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"STRAINSTEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Strain step size = " << USER_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"STRAINSTEP=",TRUE) ) {
      USER_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+"STRAINSTEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Strain step size = " << USER_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }   
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_STRAIN_STEP < tolzero) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Strain step size = " << USER_STRAIN_STEP << " <= 0.0" << endl;  
      USER_STRAIN_STEP = 0.01;
      aus << _AGLSTR_WARNING_ + "Increasing strain step size to default value of " << USER_STRAIN_STEP << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }  
    // Get the user's selection of equation of state to be used. Options are:
    // -1: Minimal output only
    // 0: Numerical equation of state
    // 1: Vinet equation of state
    // 2: Birch-Murnaghan equation of state
    // 3: Vinet equation of state using numerical evaluation of thermal properties
    // 4: Birch-Murnaghan equation of state using numerical evaluation of thermal properties
    // 5: Spinodal equation of state
    // 6: Spinodal equation of state using numerical evaluation of thermal properties
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"IEOS=",TRUE) ) {
      USER_IEOS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"IEOS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Equation of state option = " << USER_IEOS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"IEOS=",TRUE) ) {
      USER_IEOS = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"IEOS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Equation of state option = " << USER_IEOS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Get the user's selection of whether or not to optimize g in BCNT EOS. Only used for IEOS=5 or IEOS=6. Options are:
    // 0: No optimization
    // 1: Optimization
    // 2: Optimization setting initial bcnt_beta value to saved value
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"IOPT_G=",TRUE) ) {
      USER_IOPTIMIZE_BETA = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"IOPT_G=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "BCNT optimization option = " << USER_IOPTIMIZE_BETA << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"IOPT_G=",TRUE) ) {
      USER_IOPTIMIZE_BETA = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"IOPT_G=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "BCNT optimization option = " << USER_IOPTIMIZE_BETA << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's selection of whether or not to use self-consistent calculation of Debye temperature. Options are:
    // -1: Static calculation only (no finite temperature calculations)
    // 0: Non self-consistent, Debye temperature from static bulk modulus				       
    // 1: User supplied Debye temperatures
    // 2: Self-consistent calculation of Debye temperature
    // 3: Debye temperature calculated from static bulk modulus at equilibrium volume
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"IDEBYE=",TRUE) ) {
      USER_IDEBYE = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"IDEBYE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Debye temperature option = " << USER_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"IDEBYE=",TRUE) ) {
      USER_IDEBYE = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"IDEBYE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Debye temperature option = " << USER_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's value of Poisson ratio to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"POISSON_RATIO=",TRUE) ) {
      USER_POISSON_RATIO = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"POISSON_RATIO=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Poisson ratio = " << USER_POISSON_RATIO << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Warn user if USER_AEL_POISSON_RATIO is switched on as this will override the user's value
      if(USER_AEL_POISSON_RATIO.option) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
        aus << _AGLSTR_MESSAGE_ << "User supplied value of Poisson ratio = " << USER_POISSON_RATIO << " will not be used unless AEL run fails" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"POISSON_RATIO=",TRUE) ) {
      USER_POISSON_RATIO = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+"POISSON_RATIO=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Poisson ratio = " << USER_POISSON_RATIO << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Warn user if USER_AEL_POISSON_RATIO is switched on as this will override the user's value
      if(USER_AEL_POISSON_RATIO.option) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
        aus << _AGLSTR_MESSAGE_ << "User supplied value of Poisson ratio = " << USER_POISSON_RATIO << " will not be used unless AEL run fails" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"POISSONRATIO=",TRUE) ) {
      USER_POISSON_RATIO = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"POISSONRATIO=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Poisson ratio = " << USER_POISSON_RATIO << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Warn user if USER_AEL_POISSON_RATIO is switched on as this will override the user's value
      if(USER_AEL_POISSON_RATIO.option) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
        aus << _AGLSTR_MESSAGE_ << "User supplied value of Poisson ratio = " << USER_POISSON_RATIO << " will not be used unless AEL run fails" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"POISSONRATIO=",TRUE) ) {
      USER_POISSON_RATIO = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+"POISSONRATIO=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Poisson ratio = " << USER_POISSON_RATIO << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Warn user if USER_AEL_POISSON_RATIO is switched on as this will override the user's value
      if(USER_AEL_POISSON_RATIO.option) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Use AEL to calculate Poisson ratio = " << USER_AEL_POISSON_RATIO.option << endl;  
        aus << _AGLSTR_MESSAGE_ << "User supplied value of Poisson ratio = " << USER_POISSON_RATIO << " will not be used unless AEL run fails" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } 
    // Warn user if the Poisson ratio given is zero or negative
    if(USER_POISSON_RATIO < tolzero) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Poisson ratio = " << USER_POISSON_RATIO << " <= 0.0" << endl;  
      // [OBSOLETE] USER_POISSON_RATIO = 0.25;
      // [OBSOLETE] aus << _AGLSTR_WARNING_ + "Increasing Poisson ratio to default value of " << USER_POISSON_RATIO << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Warn user if the Poisson ratio given is outside the range from -1.0 to 0.5
    if(USER_POISSON_RATIO < -1.0 || USER_POISSON_RATIO > 0.5) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Poisson ratio = " << USER_POISSON_RATIO << " is outside range from -1.0 to 0.5" << endl;  
      // [OBSOLETE] USER_POISSON_RATIO = 0.25;
      // [OBSOLETE] aus << _AGLSTR_WARNING_ + "Increasing Poisson ratio to default value of " << USER_POISSON_RATIO << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Get the user's value of energy for infinite volume. This is just used to set the zero of the energy scale.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"EINF=",TRUE) ) {
      USER_EINF = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"EINF=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Energy at infinite volume = " << USER_EINF << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"EINF=",TRUE) ) {
      USER_EINF = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+"EINF=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Energy at infinite volume = " << USER_EINF << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's number of pressure points to be used in the calculation
    if(AGL_data.postprocess && (AGL_data.npressure > 0)) {
      USER_NPRESSURE = AGL_data.npressure;
    } else {
      // Otherwise, check in aflow.in file
      if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"NPRESSURE=",TRUE) ) {
        USER_NPRESSURE = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"NPRESSURE=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Number of pressure points = " << USER_NPRESSURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"NPRESSURE=",TRUE) ) {
        USER_NPRESSURE = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"NPRESSURE=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Number of pressure points = " << USER_NPRESSURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      AGL_data.npressure = USER_NPRESSURE;
    }    
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_NPRESSURE < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Number of pressure points = " << USER_NPRESSURE << " < 1" << endl;  
      USER_NPRESSURE = 11;
      aus << _AGLSTR_WARNING_ + "Increasing number of pressure points default value of " << USER_NPRESSURE << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Get the user's step size for the pressure values to be used in the calculation
    // First check if a different step size has been requested by a postprocessing command line option
    if(AGL_data.postprocess && (AGL_data.spressure > 0.0)) {
      USER_SPRESSURE = AGL_data.spressure;
    } else {
      // Otherwise, check in aflow.in file    
      if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SPRESSURE=",TRUE) ) {
        USER_SPRESSURE = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"SPRESSURE=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Pressure interval = " << USER_SPRESSURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SPRESSURE=",TRUE) ) {
        USER_SPRESSURE = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"SPRESSURE=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Pressure interval = " << USER_SPRESSURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      AGL_data.spressure = USER_SPRESSURE;
    }
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_SPRESSURE < tolzero) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Pressure interval = " << USER_SPRESSURE << " <= 0.0" << endl;  
      USER_SPRESSURE = 2.0;
      aus << _AGLSTR_WARNING_ + "Increasing pressure interval to default value of " << USER_SPRESSURE << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's number of temperature points to be used in the calculation
    // First check if a different number of points have been requested by a postprocessing command line option
    if(AGL_data.postprocess && (AGL_data.ntemperature > 0)) {
      USER_NTEMPERATURE = AGL_data.ntemperature;
    } else {
      // Otherwise, check in aflow.in file    
      if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"NTEMP=",TRUE) ) {
        USER_NTEMPERATURE = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"NTEMP=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Number of temperature points = " << USER_NTEMPERATURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"NTEMP=",TRUE) ) {
        USER_NTEMPERATURE = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"NTEMP=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Number of temperature points = " << USER_NTEMPERATURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      AGL_data.ntemperature = USER_NTEMPERATURE;
    }
    if(USER_NTEMPERATURE < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Number of temperature points = " << USER_NTEMPERATURE << " < 1" << endl;  
      aus << _AGLSTR_WARNING_ + "This number will be used, and will result in only a static zero-temperature GIBBS calculation being run" << endl;
      aus << _AGLSTR_WARNING_ + "To calculate the finite temperature properties please provide a value of NTEMP >= 1" << endl;
      aus << _AGLSTR_WARNING_ + "Alternatively, you can remove the NTEMP line from the " << _AFLOWIN_ << " file to use the default value" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's step size for the temperature values to be used in the calculation
    // First check if a different step size has been requested by a postprocessing command line option
    if(AGL_data.postprocess && (AGL_data.stemperature > 0.0)) {
      USER_STEMPERATURE = AGL_data.stemperature;
    } else {
      // Otherwise, check in aflow.in file    
      if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"STEMP=",TRUE) ) {
        USER_STEMPERATURE = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+"STEMP=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Temperature interval = " << USER_STEMPERATURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"STEMP=",TRUE) ) {
        USER_STEMPERATURE = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+"STEMP=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Temperature interval = " << USER_STEMPERATURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      AGL_data.stemperature = USER_STEMPERATURE;
    }      
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_STEMPERATURE < tolzero) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Temperature interval = " << USER_STEMPERATURE << " <= 0.0" << endl;  
      USER_STEMPERATURE = 10.0;
      aus << _AGLSTR_WARNING_ + "Increasing temperature interval to default value of " << USER_STEMPERATURE << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }   
    // Get the user's maximum limit on the number of covergence iterations to be performed in the functions scdebye and gauleg. 
    // This limit is to prevent these functions from getting stuck in an infinite loop. Default value is 250.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAX_LOOPS=",TRUE) ) {
      USER_MAX_LOOPS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAX_LOOPS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of convergence loops = " << USER_MAX_LOOPS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAX_LOOPS=",TRUE) ) {
      USER_MAX_LOOPS = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAX_LOOPS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of convergence loops = " << USER_MAX_LOOPS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAXLOOPS=",TRUE) ) {
      USER_MAX_LOOPS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAXLOOPS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of convergence loops = " << USER_MAX_LOOPS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAXLOOPS=",TRUE) ) {
      USER_MAX_LOOPS = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAXLOOPS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of convergence loops = " << USER_MAX_LOOPS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_MAX_LOOPS < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Maximum number of convergence loops = " << USER_MAX_LOOPS << " < 1" << endl;  
      USER_MAX_LOOPS = 250;
      aus << _AGLSTR_WARNING_ + "Increasing maximum number of covergence loops to default value of " << USER_MAX_LOOPS << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Get the user's maximum limit on the number of fitting iterations to be performed in the functions fitt and debfitt. 
    // Default value is 500.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAX_FITS=",TRUE) ) {
      USER_MAX_FITS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAX_FITS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of fitting iterations = " << USER_MAX_FITS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAX_FITS=",TRUE) ) {
      USER_MAX_FITS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAX_FITS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of fitting iterations = " << USER_MAX_FITS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAXFITS=",TRUE) ) {
      USER_MAX_FITS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAXFITS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of fitting iterations = " << USER_MAX_FITS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAXFITS=",TRUE) ) {
      USER_MAX_FITS = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAXFITS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of fitting iterations = " << USER_MAX_FITS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_MAX_FITS < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Maximum number of fitting iterations = " << USER_MAX_FITS << " < 1" << endl;  
      USER_MAX_FITS = 500;
      aus << _AGLSTR_WARNING_ + "Increasing maximum number of fitting iterations to default value of " << USER_MAX_FITS << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's maximum limit on the number of terms in the polynomial to be fitted to the (E, V) data. 
    // Default value is 50.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAX_PAR=",TRUE) ) {
      USER_MAX_PAR = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAX_PAR=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of terms in the polynomial = " << USER_MAX_PAR << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAX_PAR=",TRUE) ) {
      USER_MAX_PAR = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAX_PAR=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of terms in the polynomial = " << USER_MAX_PAR << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAXPAR=",TRUE) ) {
      USER_MAX_PAR = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAXPAR=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of terms in the polynomial = " << USER_MAX_PAR << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAXPAR=",TRUE) ) {
      USER_MAX_PAR = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAXPAR=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum number of terms in the polynomial = " << USER_MAX_PAR << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_MAX_PAR < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Maximum number of terms in the polynomial = " << USER_MAX_PAR << " < 1" << endl;  
      USER_MAX_PAR = 50;
      aus << _AGLSTR_WARNING_ + "Increasing maximum number of terms in the polynomial to default value of " << USER_MAX_PAR << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's maximum limit on the fitting order for the Birch-Murnaghan equation of state. 
    // Default value is 2.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"BIRCH_FIT_ORDER_IG=",TRUE) ) {
      USER_BIRCH_FIT_ORDER_IG = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"BIRCH_FIT_ORDER_IG=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum fitting order for Birch-Murnaghan equation of state = " << USER_BIRCH_FIT_ORDER_IG << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"BIRCH_FIT_ORDER_IG=",TRUE) ) {
      USER_BIRCH_FIT_ORDER_IG = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"BIRCH_FIT_ORDER_IG=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum fitting order for Birch-Murnaghan equation of state = " << USER_BIRCH_FIT_ORDER_IG << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"BIRCHFITORDER_IG=",TRUE) ) {
      USER_BIRCH_FIT_ORDER_IG = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"BIRCHFITORDER_IG=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum fitting order for Birch-Murnaghan equation of state = " << USER_BIRCH_FIT_ORDER_IG << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"BIRCHFITORDER_IG=",TRUE) ) {
      USER_BIRCH_FIT_ORDER_IG = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"BIRCHFITORDER_IG=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum fitting order for Birch-Murnaghan equation of state = " << USER_BIRCH_FIT_ORDER_IG << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_BIRCH_FIT_ORDER_IG < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Maximum number of fitting iterations = " << USER_BIRCH_FIT_ORDER_IG << " < 1" << endl;  
      USER_BIRCH_FIT_ORDER_IG = 2;
      aus << _AGLSTR_WARNING_ + "Increasing maximum number of fitting iterations to default value of " << USER_BIRCH_FIT_ORDER_IG << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }   
    // Get the user's value of the maximum number of concavity check iterations. The default value is "1".
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAX_CC_ITER=",TRUE) ) {
      USER_MAX_CC_ITER = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAX_CC_ITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) concavity check iterations = " << USER_MAX_CC_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAX_CC_ITER=",TRUE) ) {
      USER_MAX_CC_ITER = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAX_CC_ITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) concavity check iterations = " << USER_MAX_CC_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAXCCITER=",TRUE) ) {
      USER_MAX_CC_ITER = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAXCCITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) concavity check iterations = " << USER_MAX_CC_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAXCCITER=",TRUE) ) {
      USER_MAX_CC_ITER = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAXCCITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) concavity check iterations = " << USER_MAX_CC_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }         
    // Check that value is not negative; warn user and reset to default values if it is zero or negative
    if(USER_MAX_CC_ITER < 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Maximum number of concavity check iterations = " << USER_MAX_CC_ITER << " < 0" << endl;  
      USER_MAX_CC_ITER = 1;
      aus << _AGLSTR_WARNING_ + "Increasing maximum number of concavity check iterations to default value of " << USER_MAX_CC_ITER << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Get the user's value of the maximum number of minimum check iterations. The default value is "5".
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAX_CM_ITER=",TRUE) ) {
      USER_MAX_CM_ITER = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAX_CM_ITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) minimum check iterations = " << USER_MAX_CM_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAX_CM_ITER=",TRUE) ) {
      USER_MAX_CM_ITER = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAX_CM_ITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) minimum check iterations = " << USER_MAX_CM_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"MAXCMITER=",TRUE) ) {
      USER_MAX_CM_ITER = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"MAXCMITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) minimum check iterations = " << USER_MAX_CM_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"MAXCMITER=",TRUE) ) {
      USER_MAX_CM_ITER = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"MAXCMITER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Maximum (E, V) minimum check iterations = " << USER_MAX_CM_ITER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      
    // Check that value is not negative; warn user and reset to default values if it is zero or negative
    if(USER_MAX_CM_ITER < 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Maximum number of minimum check iterations = " << USER_MAX_CM_ITER << " < 0" << endl;  
      USER_MAX_CM_ITER = 5;
      aus << _AGLSTR_WARNING_ + "Increasing maximum number of minimum check iterations to default value of " << USER_MAX_CM_ITER << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }          	
    // Get the user's selection of which type of fit to perform to the (E, V) data
    // A value of "0" selects the standard GIBBS fit, skipping all convex points (default setting)
    // A value of "1" results in the fitting algorithm only skipping points where the curve changes direction twice in a row, except if this point is the global minimum
    // A value of "2" results in the fitting algorithm fitting all points 
    // A value of "3" results in the fitting algorithm selecting the same points as a value of "1", and then fitting these by a polynomial to generate a new concave set of points
    // A value of "4" results in the fitting algorithm selecting all points, and then fitting these by a polynomial to generate a new concave set of points
    // A value of "5" results (E, V) data 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"FIT_TYPE=",TRUE) ) {
      USER_FIT_TYPE = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"FIT_TYPE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Fitting algorithm selection = " << USER_FIT_TYPE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"FIT_TYPE=",TRUE) ) {
      USER_FIT_TYPE = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"FIT_TYPE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Fitting algorithm selection = " << USER_FIT_TYPE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"FITTYPE=",TRUE) ) {
      USER_FIT_TYPE = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"FITTYPE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Fitting algorithm selection = " << USER_FIT_TYPE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"FITTYPE=",TRUE) ) {
      USER_FIT_TYPE = aurostd::substring2utype<int>(AflowIn,_AGIBBSTROPT_+"FITTYPE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Fitting algorithm selection = " << USER_FIT_TYPE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Check that value is not negative; warn user and reset to default values if it is zero or negative
    if((USER_FIT_TYPE < 0) || (USER_FIT_TYPE > 5)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Fit type = " << USER_FIT_TYPE << "is outside of range 0 -> 5" << endl;  
      USER_FIT_TYPE = 0;
      aus << _AGLSTR_WARNING_ + "Setting fit type to default value of " << USER_FIT_TYPE << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Get the user's value of Gruneisen parameter to be used in the thermal conductivity calculation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"GRUNEISEN=",TRUE) ) {
      USER_GRUNEISEN.args2addattachedscheme(vAflowIn,"GRUNEISEN_EQUAL",_AGLSTROPT_+"GRUNEISEN=","");
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "Gruneisen parameter on/off = " << USER_GRUNEISEN.isflag("GRUNEISEN_EQUAL") << endl; 
      aus << _AGLSTR_MESSAGE_ << "Gruneisen parameter on/off = " << USER_GRUNEISEN.flag("GRUNEISEN_EQUAL") << endl;;
      aus << _AGLSTR_MESSAGE_ << "Gruneisen parameter = " << USER_GRUNEISEN.getattachedscheme("GRUNEISEN_EQUAL") << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"GRUNEISEN=",TRUE) ) {
      USER_GRUNEISEN.args2addattachedscheme(vAflowIn,"GRUNEISEN_EQUAL",_AGIBBSTROPT_+"GRUNEISEN=","");
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "Gruneisen parameter on/off = " << USER_GRUNEISEN.isflag("GRUNEISEN_EQUAL") << endl; 
      aus << _AGLSTR_MESSAGE_ << "Gruneisen parameter on/off = " << USER_GRUNEISEN.flag("GRUNEISEN_EQUAL") << endl;;
      aus << _AGLSTR_MESSAGE_ << "Gruneisen parameter = " << USER_GRUNEISEN.getattachedscheme("GRUNEISEN_EQUAL") << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Get the user's value of Debye temperature to be used in the thermal conductivity calculation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"THETA_COND=",TRUE) ) {
      USER_THETA_COND.args2addattachedscheme(vAflowIn,"THETA_EQUAL",_AGLSTROPT_+"THETA_COND=","");
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "Debye temperature on/off = " << USER_THETA_COND.isflag("THETA_EQUAL") << endl; 
      aus << _AGLSTR_MESSAGE_ << "Debye temperature on/off = " << USER_THETA_COND.flag("THETA_EQUAL") << endl; 
      aus << _AGLSTR_MESSAGE_ << "Debye temperature for thermal conductivity check = " << USER_THETA_COND.getattachedscheme("THETA_EQUAL") << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"THETA_COND=",TRUE) ) {
      USER_THETA_COND.args2addattachedscheme(vAflowIn,"THETA_EQUAL",_AGIBBSTROPT_+"THETA_COND=","");
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "Debye temperature on/off = " << USER_THETA_COND.isflag("THETA_EQUAL") << endl; 
      aus << _AGLSTR_MESSAGE_ << "Debye temperature on/off = " << USER_THETA_COND.flag("THETA_EQUAL") << endl; 
      aus << _AGLSTR_MESSAGE_ << "Debye temperature for thermal conductivity check = " << USER_THETA_COND.getattachedscheme("THETA_EQUAL") << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to generate plots of the Debye temperature, heat capacity, Gruneisen parameter and vibrational free energy as a function of temperature. 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"PLOT_RESULTS=",TRUE) ) {
      USER_PLOT_RESULTS.options2entry(AflowIn,_AGLSTROPT_+"PLOT_RESULTS=",USER_PLOT_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_PLOT_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"PLOT_RESULTS=",TRUE) ) {
      USER_PLOT_RESULTS.options2entry(AflowIn,_AGIBBSTROPT_+"PLOT_RESULTS=",USER_PLOT_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_PLOT_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"PLOTRESULTS=",TRUE) ) {
      USER_PLOT_RESULTS.options2entry(AflowIn,_AGLSTROPT_+"PLOTRESULTS=",USER_PLOT_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_PLOT_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"PLOTRESULTS=",TRUE) ) {
      USER_PLOT_RESULTS.options2entry(AflowIn,_AGIBBSTROPT_+"PLOTRESULTS=",USER_PLOT_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_PLOT_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      	
    // Get the user's selection of whether to calculate the thermal conductivity as a function of temperature using the volume as a function of temperature
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"KAPPA_VOLUME=",TRUE) ) {
      USER_KAPPA_VOLUME.options2entry(AflowIn,_AGLSTROPT_+"KAPPA_VOLUME=",USER_KAPPA_VOLUME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_KAPPA_VOLUME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"KAPPA_VOLUME=",TRUE) ) {
      USER_KAPPA_VOLUME.options2entry(AflowIn,_AGIBBSTROPT_+"KAPPA_VOLUME=",USER_KAPPA_VOLUME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_KAPPA_VOLUME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"KAPPAVOLUME=",TRUE) ) {
      USER_KAPPA_VOLUME.options2entry(AflowIn,_AGLSTROPT_+"KAPPAVOLUME=",USER_KAPPA_VOLUME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_KAPPA_VOLUME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"KAPPAVOLUME=",TRUE) ) {
      USER_KAPPA_VOLUME.options2entry(AflowIn,_AGIBBSTROPT_+"KAPPAVOLUME=",USER_KAPPA_VOLUME.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Plot AGL results = " << USER_KAPPA_VOLUME.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's selection of whether to get the atomic species from the potential name in order to calculate the unit cell mass
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SPECIES_MASS=",TRUE) ) {
      USER_SPECIES_MASS.options2entry(AflowIn,_AGLSTROPT_+"SPECIES_MASS=",USER_SPECIES_MASS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Obtain atomic species from potential name = " << USER_SPECIES_MASS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SPECIES_MASS=",TRUE) ) {
      USER_SPECIES_MASS.options2entry(AflowIn,_AGIBBSTROPT_+"SPECIES_MASS=",USER_SPECIES_MASS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Obtain atomic species from potential name = " << USER_SPECIES_MASS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SPECIESMASS=",TRUE) ) {
      USER_SPECIES_MASS.options2entry(AflowIn,_AGLSTROPT_+"SPECIESMASS=",USER_SPECIES_MASS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Obtain atomic species from potential name = " << USER_SPECIES_MASS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"SPECIESMASS=",TRUE) ) {
      USER_SPECIES_MASS.options2entry(AflowIn,_AGIBBSTROPT_+"SPECIESMASS=",USER_SPECIES_MASS.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Obtain atomic species from potential name = " << USER_SPECIES_MASS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }       		
    // Get the user's selection of whether to use the prefix "ARUN" before the individual directory names
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"DIRNAME_ARUN=",TRUE) ) {
      USER_DIRNAME_ARUN.options2entry(AflowIn,_AGLSTROPT_+"DIRNAME_ARUN=",USER_DIRNAME_ARUN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use ARUN prefix in directory names = " << USER_DIRNAME_ARUN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"DIRNAME_ARUN=",TRUE) ) {
      USER_DIRNAME_ARUN.options2entry(AflowIn,_AGIBBSTROPT_+"DIRNAME_ARUN=",USER_DIRNAME_ARUN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use ARUN prefix in directory names = " << USER_DIRNAME_ARUN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"DIRNAMEARUN=",TRUE) ) {
      USER_DIRNAME_ARUN.options2entry(AflowIn,_AGLSTROPT_+"DIRNAMEARUN=",USER_DIRNAME_ARUN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use ARUN prefix in directory names = " << USER_DIRNAME_ARUN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"DIRNAMEARUN=",TRUE) ) {
      USER_DIRNAME_ARUN.options2entry(AflowIn,_AGIBBSTROPT_+"DIRNAMEARUN=",USER_DIRNAME_ARUN.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Use ARUN prefix in directory names = " << USER_DIRNAME_ARUN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }          		
    // Get the user's selection of whether to write the debugging information for the function AGL_functions::gaussxm
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"GAUSSXM_DEBUG=",TRUE) ) {
      USER_GAUSSXM_DEBUG.options2entry(AflowIn,_AGLSTROPT_+"GAUSSXM_DEBUG=",USER_GAUSSXM_DEBUG.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write AGL_functions::gaussxm debugging information = " << USER_GAUSSXM_DEBUG.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"GAUSSXM_DEBUG=",TRUE) ) {
      USER_GAUSSXM_DEBUG.options2entry(AflowIn,_AGIBBSTROPT_+"GAUSSXM_DEBUG=",USER_GAUSSXM_DEBUG.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write AGL_functions::gaussxm debugging information = " << USER_GAUSSXM_DEBUG.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"GAUSSXMDEBUG=",TRUE) ) {
      USER_GAUSSXM_DEBUG.options2entry(AflowIn,_AGLSTROPT_+"GAUSSXMDEBUG=",USER_GAUSSXM_DEBUG.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write AGL_functions::gaussxm debugging information = " << USER_GAUSSXM_DEBUG.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGIBBSTROPT_+"GAUSSXMDEBUG=",TRUE) ) {
      USER_GAUSSXM_DEBUG.options2entry(AflowIn,_AGIBBSTROPT_+"GAUSSXMDEBUG=",USER_GAUSSXM_DEBUG.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Write AGL_functions::gaussxm debugging information = " << USER_GAUSSXM_DEBUG.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which values to use for PREC and ALGO.
    // By setting this variable to ON, the values for PREC=ACCURATE and ALGO=NORMAL will be used for AGL, overriding other settings in the aflow.in file.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"PRECACCALGONORM=",TRUE) ) {
      USER_PRECACC_ALGONORM.options2entry(AflowIn,_AGLSTROPT_+"PRECACCALGONORM=",USER_PRECACC_ALGONORM.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Use PREC=ACCURATE and ALGO=NORMAL = " << USER_PRECACC_ALGONORM.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"PRECACC_ALGONORM=",TRUE) ) {
      USER_PRECACC_ALGONORM.options2entry(AflowIn,_AGLSTROPT_+"PRECACC_ALGONORM=",USER_PRECACC_ALGONORM.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Use PREC=ACCURATE and ALGO=NORMAL = " << USER_PRECACC_ALGONORM.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which type of runs to do to obtain the thermal properties
    // The default is to use STATIC, independently of what is used in these lines in the aflow.in file. 
    // By setting RELAX_STATIC to ON, the run type RELAX_STATIC=2 will be used.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"RELAX_STATIC=",TRUE) ) {
      USER_RELAX_STATIC.options2entry(AflowIn,_AGLSTROPT_+"RELAX_STATIC=",USER_RELAX_STATIC.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Use RELAX_STATIC=2 = " << USER_RELAX_STATIC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which type of runs to do to obtain the elastic constants
    // The default is to use STATIC, independently of what is used in these lines in the aflow.in file. 
    // By setting STATIC to ON, the run type STATIC will be used.
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"STATIC=",TRUE) ) {
      USER_STATIC.options2entry(AflowIn,_AGLSTROPT_+"STATIC=",USER_STATIC.option);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Use STATIC = " << USER_STATIC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's list of failed run subdirectories to skip when calculating the thermal properties
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SKIP_FAILED_ARUNS=",TRUE) ) {
      USER_SKIP_FAILED_ARUNS = aurostd::substring2string(AflowIn,_AGLSTROPT_+"SKIP_FAILED_ARUNS=",1,TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Failed ARUNS to skip = " << USER_SKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SKIPFAILEDARUNS=",TRUE) ) {
      USER_SKIP_FAILED_ARUNS = aurostd::substring2string(AflowIn,_AGLSTROPT_+"SKIPFAILEDARUNS=",1,TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Failed ARUNS to skip = " << USER_SKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }	
    // Get the user's option of whether to automatically skip failed run subdirectories 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"AUTOSKIP_FAILED_ARUNS=",TRUE) ) {
      USER_AUTOSKIP_FAILED_ARUNS.options2entry(AflowIn,_AGLSTROPT_+"AUTOSKIP_FAILED_ARUNS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Automatically detect and skip failed ARUNS = " << USER_AUTOSKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"AUTOSKIPFAILEDARUNS=",TRUE) ) {
      USER_AUTOSKIP_FAILED_ARUNS.options2entry(AflowIn,_AGLSTROPT_+"AUTOSKIPFAILEDARUNS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Automatically detect and skip failed ARUNS = " << USER_AUTOSKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }	
    // Get the user's option of maximum number of failed run subdirectories to skip 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SKIP_ARUNS_MAX=",TRUE) ) {
      USER_SKIP_ARUNS_MAX = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"SKIP_ARUNS_MAX=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Maximum number of ARUNS to skip = " << USER_SKIP_ARUNS_MAX << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"SKIPARUNSMAX=",TRUE) ) {
      USER_SKIP_ARUNS_MAX = aurostd::substring2utype<int>(AflowIn,_AGLSTROPT_+"SKIPARUNSMAX=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Maximum number of ARUNS to skip = " << USER_SKIP_ARUNS_MAX << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 	
    // Get the user's option of whether to write the input data for calculating the Hugoniot relation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITE_HUGONIOT_INPUT=",TRUE) ) {
      USER_WRITE_HUGONIOT_INPUT.options2entry(AflowIn,_AGLSTROPT_+"WRITE_HUGONIOT_INPUT=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Write Hugoniot input = " << USER_WRITE_HUGONIOT_INPUT << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"WRITEHUGONIOTINPUT=",TRUE) ) {
      USER_WRITE_HUGONIOT_INPUT.options2entry(AflowIn,_AGLSTROPT_+"WRITEHUGONIOTINPUT=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Write Hugoniot input = " << USER_WRITE_HUGONIOT_INPUT << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's option of whether to write the input data for calculating the Hugoniot relation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"HUGONIOT_CALC=",TRUE) ) {
      USER_HUGONIOT_CALC.options2entry(AflowIn,_AGLSTROPT_+"HUGONIOT_CALC=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Calculate Hugoniot relation = " << USER_HUGONIOT_CALC << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"HUGONIOTCALC=",TRUE) ) {
      USER_HUGONIOT_CALC.options2entry(AflowIn,_AGLSTROPT_+"HUGONIOTCALC=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Calculate Hugoniot relation = " << USER_HUGONIOT_CALC << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's option of whether to write the input data for calculating the Hugoniot relation
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"HUGONIOT_EXTRAPOLATE=",TRUE) ) {
      USER_HUGONIOT_EXTRAPOLATE.options2entry(AflowIn,_AGLSTROPT_+"HUGONIOT_EXTRAPOLATE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Extrapolate data for Hugoniot calculation = " << USER_HUGONIOT_EXTRAPOLATE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"HUGONIOTEXTRAPOLATE=",TRUE) ) {
      USER_HUGONIOT_EXTRAPOLATE.options2entry(AflowIn,_AGLSTROPT_+"HUGONIOTEXTRAPOLATE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Extrapolate data for Hugoniot calculation = " << USER_HUGONIOT_EXTRAPOLATE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's option of whether to write the input data for whether to avoid truncating pressure or temperature range if no minimum energy is found 
    if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"RUN_ALL_PRESSURE_TEMPERATURE=",TRUE) ) {
      USER_RUN_ALL_PRESSURE_TEMPERATURE.options2entry(AflowIn,_AGLSTROPT_+"RUN_ALL_PRESSURE_TEMPERATURE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Run all pressures and temperatures = " << USER_RUN_ALL_PRESSURE_TEMPERATURE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"RUNALLPRESSURETEMPERATURE=",TRUE) ) {
      USER_RUN_ALL_PRESSURE_TEMPERATURE.options2entry(AflowIn,_AGLSTROPT_+"RUNALLPRESSURETEMPERATURE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Run all pressures and temperatures = " << USER_RUN_ALL_PRESSURE_TEMPERATURE << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }	

    // If this is an AGL run to get the pressure-volume data for finite pressure AEL runs, then get user's options for AEL pressures
    if (AGL_data.ael_pressure_calc) {
      // Get the user's number of pressure points to be used in the calculation
      if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NPRESSURE=",TRUE) ) {
        USER_NPRESSURE = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NPRESSURE=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Number of pressure points = " << USER_NPRESSURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }  
      // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
      if(USER_NPRESSURE < 1) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "Number of pressure points = " << USER_NPRESSURE << " < 1" << endl;  
        USER_NPRESSURE = 11;
        aus << _AGLSTR_WARNING_ + "Increasing number of pressure points default value of " << USER_NPRESSURE << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } 
      // Get the user's step size for the pressure values to be used in the calculation
      if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SPRESSURE=",TRUE) ) {
        USER_SPRESSURE = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"SPRESSURE=",TRUE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Pressure interval = " << USER_SPRESSURE << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } 
      // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
      if(USER_SPRESSURE < tolzero) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "Pressure interval = " << USER_SPRESSURE << " <= 0.0" << endl;  
        USER_SPRESSURE = 2.0;
        aus << _AGLSTR_WARNING_ + "Increasing pressure interval to default value of " << USER_SPRESSURE << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    double presi = 0.0;
    double tempi = 0.0;
    vector<double> pres_orig;
    pres_orig.clear();

    // Set variables in class AGL_data to user-defined values
    AGL_data.i_eqn_of_state = USER_IEOS;
    AGL_data.i_optimize_beta = USER_IOPTIMIZE_BETA;
    AGL_data.i_debye = USER_IDEBYE;
    // [OBSOLETE] AGL_data.poissonratio = USER_POISSON_RATIO;
    AGL_data.energy_infinity = USER_EINF;
    AGL_data.maxloops = USER_MAX_LOOPS;
    AGL_data.maxfit = USER_MAX_FITS;
    AGL_data.maxpolycoeffs = USER_MAX_PAR;
    AGL_data.birchfitorder_iG = USER_BIRCH_FIT_ORDER_IG;
    AGL_data.fittype = USER_FIT_TYPE;
    AGL_data.gaussxm_debug = USER_GAUSSXM_DEBUG.option;
    AGL_data.precaccalgonorm = USER_PRECACC_ALGONORM.option;
    AGL_data.relax_static = USER_RELAX_STATIC.option;
    AGL_data.static_only = USER_STATIC.option;
    AGL_data.autoskipfailedaruns = USER_AUTOSKIP_FAILED_ARUNS.option;
    AGL_data.hugoniotrun = USER_HUGONIOT_CALC.option;
    AGL_data.hugoniotextrapolate = USER_HUGONIOT_EXTRAPOLATE.option;

    // STATIC is the default; if the RELAXSTATIC option is true it should be false; otherwise it should be true
    if(AGL_data.relax_static) {
      AGL_data.static_only = false;
    } else {
      AGL_data.static_only = true;
    }

    // Insert user's values for maximum number of failed runs to skip
    AGL_data.skiparunsmax = USER_SKIP_ARUNS_MAX;
    // Insert user's values for total number of runs
    // Will be used to check number eliminated in error correction does not exceed skiparunsmax
    AGL_data.nstructsinit = USER_NSTRUCTURES;

    // Check that user's option for whether to avoid truncating pressures or temperature ranges if no minimum energy is found is consistent with the equation of state being used
    // The semi-empirical equations of state fit to pressure-volume data, and so the fitting algorithms require all pressure values including zero
    // Therefore, running all pressure/temperature values, skipping low pressure values, requires the use of the numerical equation of state
    if (USER_RUN_ALL_PRESSURE_TEMPERATURE.option) {
      if (AGL_data.i_eqn_of_state != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "Run all pressures and temperatures = " << USER_RUN_ALL_PRESSURE_TEMPERATURE << endl;
        aus << _AGLSTR_WARNING_ + "This requires use of the numerical equation of state (IEOS=0)" << endl;
        aus << _AGLSTR_WARNING_ + "Currently, IEOS=" << AGL_data.i_eqn_of_state << endl;
        aus << _AGLSTR_WARNING_ + "Running standard AGL run, including truncating pressure/temperature ranges when fitting fails" << endl;
        aus << _AGLSTR_WARNING_ + "To avoid truncation, set IEOS=0 in input" << endl;	
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        AGL_data.run_all_pressure_temperature = false;
      } else {
        AGL_data.run_all_pressure_temperature = true;
      }
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "AGL_data class constructed" << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // If this is an AGL run to get the pressure-volume data for finite pressure AEL runs, switches off AEL Poisson ratio
    // Otherwise, AEL-AGL will could end up calling each other in an infinite loop
    if (AGL_data.ael_pressure_calc) {
      USER_AEL_POISSON_RATIO.option = false;
    }

    // Calls AFLOW AEL to calculate the Poisson ratio if selected by user
    if(USER_AEL_POISSON_RATIO.option) {
      aglerror = AEL_functions::Get_PoissonRatio(xvasp, AflowIn, aflags, kflags, vflags, AGL_data.poissonratio, AGL_data.postprocess, FileMESSAGE);
      if(aglerror == 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "AEL successfully used to calculate Poisson ratio = " << AGL_data.poissonratio << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if((AGL_data.poissonratio < -1.0) || (AGL_data.poissonratio > 0.5)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Poisson ratio = " << AGL_data.poissonratio << " is outside range from -1.0 to 0.5" << endl;
          aus << _AGLSTR_ERROR_ + "Poisson ratios outside of this range will not give real values for the Poisson ratio function in AGL" << endl;
          aus << _AGLSTR_ERROR_ + "Check AEL calculations to ensure that they have run correctly and that the material is not unstable" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          aglerror = 1;
          return aglerror;
        }
        // Record source and value of Poisson ratio used
        stringstream prss;
        prss << AGL_data.poissonratio;
        string prs = prss.str();
        AGL_data.poissonratiosource = "ael_poisson_ratio_" + prs;	
      } else if(aglerror == 8) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "AEL waiting for other calculations to run" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        //throw AGLStageBreak();
        return aglerror;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ << "AEL calculation failed" << endl;
        aus << _AGLSTR_WARNING_ << "Using default / user provided value of Poisson ratio = " << USER_POISSON_RATIO << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        AGL_data.poissonratio = USER_POISSON_RATIO;
        // Record source and value of Poisson ratio used
        if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"POISSONRATIO=",TRUE) ) {
          stringstream prss;
          prss << AGL_data.poissonratio;
          string prs = prss.str();
          AGL_data.poissonratiosource = "user_defined_value_" + prs;
        } else {
          AGL_data.poissonratiosource = "default_value_Cauchy_solid_0.25";
        }
      }	
    } else {
      // Use default or user supplied value of Poisson ratio if AEL Poisson ratio is not selected by user
      AGL_data.poissonratio = USER_POISSON_RATIO;
      // Record source and value of Poisson ratio used
      if( aurostd::substring2bool(AflowIn,_AGLSTROPT_+"POISSONRATIO=",TRUE) || (aurostd::substring2bool(AflowIn,_AGLSTROPT_+"POISSON_RATIO=",TRUE)) ) {
        stringstream prss;
        prss << AGL_data.poissonratio;
        string prs = prss.str();
        AGL_data.poissonratiosource = "user_defined_value_" + prs;
      } else {
        AGL_data.poissonratiosource = "default_value_Cauchy_solid_0.25";
      }
    }

    // Set user's selection of whether or not to save thermal properties for all pressures
    if(USER_SAVE_ALL_PRESSURES.option) {
      AGL_data.savedatapressure = true;
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Saving results for all pressure values" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // Convert units of energy value at infinity from eV to Hartree for GIBBS input
    // [OBSOLETE] AGL_data.energy_infinity = AGL_data.energy_infinity * ev2hart1;
    AGL_data.energy_infinity = AGL_data.energy_infinity / hart2ev;
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "AGL_data.energy_infinity = " << AGL_data.energy_infinity << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Calculate values of pressure and temperature points.

    // Read in user supplied pressure values (if present)
    // Activated by including [AFLOW_AGL]SETPRESSUREVALUE=ON or [AFLOW_GIBBS]SETPRESSUREVALUE=ON in _AFLOWIN_ file
    if(USER_SET_PRESSURE_VALUE.option) {
      for (uint ids = 0; ids < USER_NPRESSURE; ids++) {
        string strids;
        ostringstream cnvids;
        cnvids << ids;
        strids = cnvids.str();
        string stridsafl = _AGLSTROPT_ + "PRESSURE_VALUE_" + strids + "=";
        string stridsagbfl = _AGIBBSTROPT_ + "PRESSURE_VALUE_" + strids + "=";
        string stridsaflt = "PRESSURE_VALUE_" + strids + "=";
        if( aurostd::substring2bool(AflowIn,stridsafl,TRUE) || aurostd::substring2bool(AflowIn,stridsagbfl,TRUE) ) {
          if( aurostd::substring2bool(AflowIn,stridsafl,TRUE) ) {
            USER_PRES_VAL = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+stridsaflt,TRUE);
          } else if( aurostd::substring2bool(AflowIn,stridsagbfl,TRUE) ) {
            USER_PRES_VAL = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+stridsaflt,TRUE);
          }
          if(USER_PRES_VAL < 0.0) {
            // Skip negative pressure values
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Error in user supplied values of pressures" << endl;
            aus << _AGLSTR_ERROR_ + "Pressure value = " << USER_PRES_VAL << endl;
            aus << _AGLSTR_ERROR_ + "Pressure values should not be negative" << endl;
            aus << _AGLSTR_ERROR_ + "Check section in AFLOW manual on GIBBS for details" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          } else if(ids == 0 && USER_PRES_VAL > 0.0) {
            // Check that first pressure value = 0.0; if not, adds in extra pressure point at start with p = 0.0
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Error in user supplied values of pressures from " << _AFLOWIN_ << "" << endl;
            aus << _AGLSTR_WARNING_ + "First value of SETPRESVAL = " << USER_PRES_VAL << endl;
            aus << _AGLSTR_WARNING_ + "First pressure value should be equal to zero, setting first pressure value to zero" << endl;
            aus << _AGLSTR_WARNING_ + "Check section in AFLOW manual on GIBBS for details" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            AGL_data.pressure_external.push_back(0.0);
            AGL_data.pressure_external.push_back(USER_PRES_VAL);
          } else {	  
            // Otherwise, simply add pressure value to list
            AGL_data.pressure_external.push_back(USER_PRES_VAL);
          }
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error reading in user supplied values of pressures from " << _AFLOWIN_ << "" << endl;
          aus << _AGLSTR_ERROR_ + "Value of SETPRESVAL = " << USER_SET_PRESSURE_VALUE.option << endl;
          aus << _AGLSTR_ERROR_ + "This requires user to supply NPRESSURE pressure values, where NPRESSURE = " << USER_NPRESSURE << endl;
          aus << _AGLSTR_ERROR_ + "Check section in AFLOW manual on GIBBS for details" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
      }
    } else {
      // Generate pressure values using stepsize and number of values
      // Resize vector in advance to allocate memory to store pressure values 
      AGL_data.pressure_external.resize(USER_NPRESSURE);
      for(uint i = 0; i < USER_NPRESSURE; i++) {
        AGL_data.pressure_external.at(i) = presi; 
        presi = presi + USER_SPRESSURE;
      }
    }
    for (uint i = 0; i < AGL_data.pressure_external.size(); i++) {
      pres_orig.push_back(AGL_data.pressure_external.at(i));
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Pressure values generated successfully" << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Generate temperature values using stepsize and number of values
    // Resize vector in advance to allocate memory to store temperature values
    AGL_data.temperature_external.resize(USER_NTEMPERATURE);
    for(uint i = 0; i < USER_NTEMPERATURE; i++) {
      AGL_data.temperature_external.at(i) = tempi; 
      tempi = tempi + USER_STEMPERATURE;
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Temperature values generated successfully" << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Read in user supplied Debye temperatures (if present)
    // Activated by including [AFLOW_AGL]IDEBYE=1 or [AFLOW_GIBBS]IDEBYE=1 in _AFLOWIN_ file
    if(USER_IDEBYE == 1) {
      for (int ids = 0; ids < USER_NSTRUCTURES; ids++) {
        string strids;
        ostringstream cnvids;
        cnvids << ids;
        strids = cnvids.str();
        string stridsafl = _AGLSTROPT_ + "DEBYE_TEMP_" + strids + "=";
        string stridsagbfl = _AGIBBSTROPT_ + "DEBYE_TEMP_" + strids + "=";
        string stridsaflt = "DEBYE_TEMP_" + strids + "=";
        if( aurostd::substring2bool(AflowIn,stridsafl,TRUE) ) {
          USER_DEBYE_TEMP = aurostd::substring2utype<double>(AflowIn,_AGLSTROPT_+stridsaflt,TRUE);
          AGL_data.tdebye.push_back(USER_DEBYE_TEMP);
        } else if( aurostd::substring2bool(AflowIn,stridsagbfl,TRUE) ) {
          USER_DEBYE_TEMP = aurostd::substring2utype<double>(AflowIn,_AGIBBSTROPT_+stridsaflt,TRUE);
          AGL_data.tdebye.push_back(USER_DEBYE_TEMP);
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error reading in values of Debye temperatures from " << _AFLOWIN_ << "" << endl;
          aus << _AGLSTR_ERROR_ + "Value of IDEBYE = " << USER_IDEBYE << endl;
          aus << _AGLSTR_ERROR_ + "This requires user to supply Debye temperature for each structure" << endl;
          aus << _AGLSTR_ERROR_ + "Check section in AFLOW manual on GIBBS for details" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
      }
    }

    // Set up name of calculation and output filename and directory path
    AGL_data.dirpathname = aurostd::CleanFileName(xvasp.Directory);
    AGL_data.sysname = aurostd::CleanStringASCII(USER_SYSTEM_NAME);
    AGL_data.sysname = aurostd::RemoveWhiteSpaces(AGL_data.sysname);

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Output directory name = " << AGL_data.dirpathname << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Check if there are failed ARUN subdirectories to skip; if so tokenize into vector of strings
    if(USER_SKIP_FAILED_ARUNS != "") {
      vector<string> tokens;
      string failed_arun;
      aurostd::string2tokens(USER_SKIP_FAILED_ARUNS, tokens, ",");
      for(uint i = 0; i < tokens.size(); i++) {
        failed_arun = aurostd::RemoveWhiteSpaces(tokens.at(i));
        AGL_data.failed_arun_list.push_back(failed_arun);
      }
      for (uint i = 0; i < AGL_data.failed_arun_list.size(); i++) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Failed ARUN directory = " << AGL_data.failed_arun_list.at(i) << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // [OBSOLETE] return 1;
    }

    // Calculate number of atoms in unit cell of system
    AGL_data.natoms = xvasp.str.atoms.size();
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Number of atoms = " << AGL_data.natoms << endl;  
    aus << _AGLSTR_MESSAGE_ << "xvasp.str.atoms.size() = " << xvasp.str.atoms.size() << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Calculate total mass of unit cell of system 
    if(USER_SPECIES_MASS.option) {
      // Use species pp names 
      for(uint i = 0; i < xvasp.str.species.size(); i++) {
        AGL_data.cellmass = AGL_data.cellmass + xvasp.str.num_each_type.at(i) * GetAtomMass(xvasp.str.species.at(i));
      }

      // Write out the names of the atoms used in the calculation
      uint ij = 0;
      for(uint i = 0; i < xvasp.str.species.size(); i++) {
        for (int j = 0; j < xvasp.str.num_each_type.at(i); j++) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species.at(i) << endl; 
          aus << _AGLSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species_pp.at(i) << endl; 
          aus << _AGLSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species_pp_type.at(i) << endl; 
          aus << _AGLSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species_pp_version.at(i) << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          ij++;
        }
      }
    } else {
      // Use atom names
      for(uint i = 0; i < xvasp.str.atoms.size(); i++) {
        AGL_data.cellmass = AGL_data.cellmass + GetAtomMass(xvasp.str.atoms.at(i).cleanname);
      }

      // Write out the names of the atoms used in the calculation
      for(uint i = 0; i < xvasp.str.atoms.size(); i++) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Atom " << i << " = " << xvasp.str.atoms.at(i).cleanname << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    // [OBSOLETE] Calculate conversion factor to convert atomic mass units (AMU) to atomic units
    // [OBSOLETE] double amu2au = physconstamu/physconstme;

    // Get average mass in kg for calculation of thermal conductivity at end of AGL run
    double avmasskg = AGL_data.cellmass / AGL_data.natoms;
    // Get total cell mass in kg for calculation of mass density
    double cellmasskg = AGL_data.cellmass;

    // [OBSOLETE] Convert mass of cell from kg to AMU and then to atomic units for GIBBS input
    // Convert mass of cell from kg to AMU
    AGL_data.cellmass /= AMU2KILOGRAM;
    // [OBSOLETE] AGL_data.cellmass *= amu2au;

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Mass of cell = " << AGL_data.cellmass << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "starting try loop to create strained structures and run GIBBS" << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Create set of strained structures, write _AFLOWIN_ files to directories, run VASP, collect results, and pass to GIBBS method
    // This "try" loop is adapted from that in AFLOW APL function DirectMethodPC::runVASPCalculations()
    try {
      int mid = USER_NSTRUCTURES / 2;
      int imid;
      double dimid, strainfactor;
      xstructure initialstructure = xvasp.str;
      double initialvolume = xvasp.str.Volume();

      xstructure strainedstructure;
      _xvasp vaspRun;
      vector<_xvasp> vaspRuns;
      _vflags vaspFlags;
      vaspFlags = vflags;
      _kflags kbinFlags;
      kbinFlags = kflags;
      _aflags aflowFlags;
      aflowFlags = aflags;
      vector<string> runname;
      vector<string> dirrunname;
      vector<string> arunname;
      vector<double> strainfactorlist;
      bool avasp_pop_xvasp_success = true;
      bool aflowin_success = true;
      stringstream arunaflowin;
      bool write = true;      

      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "creating strained structures" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      bool skipdir = false;
      string dfilename;      
      // Create set of strained structures
      for(int i = 0; i < USER_NSTRUCTURES; i++) {
        imid = i - mid;
        dimid = imid;
        strainedstructure = initialstructure;;
        strainfactor = (dimid * USER_STRAIN_STEP) + 1.0;
        strainedstructure.InflateLattice(strainfactor);
        strainfactorlist.push_back(strainfactor);

        // Create vector of classes of type _xvasp
        // Each element will contain the input and output data for the VASP run for one strained structure
        vaspRun = xvasp;
        vaspRuns.push_back(vaspRun);
        int idVaspRun = vaspRuns.size()-1;

        // Set up name of separate directory for AFLOW run for each strained structure
        // [OBSOLETE] string stridVaspRun;
        ostringstream cnvidVaspRun;
        cnvidVaspRun << idVaspRun;
        string stridVaspRun = cnvidVaspRun.str();

        // [OBSOLETE] string strstrainfactor;
        ostringstream cnvstrainfactor;
        cnvstrainfactor << strainfactor;
        string strstrainfactor = cnvstrainfactor.str();

        // [OBSOLETE] string runname;
        if(USER_DIRNAME_ARUN.option) {
          // [OBSOLETE] runname = _ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor;
          runname.push_back(_ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor);
        } else {
          // [OBSOLETE] runname = _AGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor;
          runname.push_back(_AGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor);
        }
        arunname.push_back(stridVaspRun + "_SF_" + strstrainfactor);
        dirrunname.push_back(AGL_data.dirpathname + "/" + runname.at(idVaspRun));

        // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AGL_data.dirpathname + "/" + runname;
        // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AGL_data.dirpathname + "/" + runname.at(idVaspRun);
        vaspRuns.at(idVaspRun).str = strainedstructure;
        vaspRuns.at(idVaspRun).AVASP_arun_runname = arunname.at(idVaspRun);
        vaspRuns.at(idVaspRun).AVASP_arun = true;
        vaspRuns.at(idVaspRun).AVASP_arun_mode = "AGL";

        // Assign the values of the flags provided by the user in the _AFLOWIN_ file to the class containing the input data for the VASP run
        aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), vaspFlags, kbinFlags, dirrunname.at(idVaspRun), AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to assign values of flags from " << _AFLOWIN_ << " file" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }	
        avasp_pop_xvasp_success = AVASP_populateXVASP(aflowFlags, kbinFlags, vaspFlags, vaspRuns.at(idVaspRun));
        if (!avasp_pop_xvasp_success) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to set AFLOW flags for ARUN " << arunname.at(idVaspRun) << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }

        // If there are already LOCK or OUTCAR.static files in this directory, it means this directory was already generated and computed.
        // Therefore do not touch, but store this structure in the list, so that it can be used in the next part of code.
        if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/"+_AFLOWLOCK_ ) ||
            aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
            aurostd::EFileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
            aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) ||
            aurostd::FileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ||
            aurostd::EFileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ||
            ((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) &&
             ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) ||
              (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
              (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
              (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK")))) ) continue;

        // Check if structure is on list of failed runs to be skipped
        // If so, then skip generation and continue to next structure
        for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
          // [OBSOLETE] if(aurostd::substring2bool(runname,AGL_data.failed_arun_list.at(ij),TRUE))
          if(aurostd::substring2bool(runname.at(idVaspRun),AGL_data.failed_arun_list.at(ij),TRUE)) {
            skipdir = true;
          }
        }	
        if(skipdir) {
          aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << runname << endl;
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << runname.at(idVaspRun) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = false;
          continue;
        }

        // If files do not exist, and the postprocess flag is not set, continue on to prepare generation of _AFLOWIN_ ...
        if (!(XHOST.POSTPROCESS || XHOST.ARUN_POSTPROCESS || AGL_data.postprocess)) { //CO20210701
          // Assign the values of the flags provided by the user in the _AFLOWIN_ file to the class containing the input data for the VASP run
          // [OBSOLETE] aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname, FileMESSAGE);
          // [OBSOLETE] aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname.at(idVaspRun), FileMESSAGE);
          // [OBSOLETE] if(aglerror != 0) {
          // [OBSOLETE]   aurostd::StringstreamClean(aus);
          // [OBSOLETE]   aus << _AGLSTR_ERROR_ + "Failed to assign values of flags from " << _AFLOWIN_ << " file" << endl;  
          // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE]   return aglerror;
          // [OBSOLETE] }

          // Create _AFLOWIN_
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "writing " << _AFLOWIN_ << "" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

          // Writes _AFLOWIN_ file to appropriate directory for each VASP run
          // [OBSOLETE] aglerror = AGL_functions::createAFLOWIN(vaspRuns.at(idVaspRun), xvasp, _kbinFlags, _vaspFlags, AGL_data, FileMESSAGE);
          aflowin_success = AVASP_MakeSingleAFLOWIN(vaspRuns.at(idVaspRun), arunaflowin, write);
          if(!aflowin_success) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to create " << _AFLOWIN_ << " files" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }
        }
      }

      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "extract energies" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // Extract final energies from results of VASP calculations for strained structures
      aglerror = AGL_functions::extractenerg(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failed to extract total energies" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aglerror;
      }
      // Extract electronic DOS from results of VASP calculations for strained structures
      aglerror = AGL_functions::extractedos(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failed to extract electronic DOS" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aglerror;
      }

      // Checks if the (E, V) data sets contains a minimum in the energy, which is required by the GIBBS algorithm
      // If there is no minimum in the data passed to it, the GIBBS algorithm will have problems fitting the data by a polynomial, and will exit with an error
      // If the data does not have a minimum, then this section of the code will add additional strained structures until a minimum is found
      // If the structure has been properly relaxed prior to running AFLOW AGL, then this problem should not occur
      // If this part of the code produces warnings, you should check the relaxation of the initial structure
      if(USER_CHECK_EV_MIN.option && !(XHOST.ARUN_POSTPROCESS || AGL_data.postprocess)) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Checking location of energy minimum" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // The function "checkmin" checks for a minimum
        // If a minimum is found, the value of cm is set to 0
        // If the lowest energy corresponds to the smallest volume, the value of cm is set to 1
        // If the lowest energy corresponds to the largest volume, the value of cm is set to 2
        // Default to assuming that minimum is present to avoid running unnecessary additional calculations
        int cm = 0;
        aglerror = AGL_functions::checkmin(AGL_data, cm, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to find energy minimum" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        int i = USER_NSTRUCTURES - vaspRuns.size();
        int ncm = 0;
        // While the lowest energy corresponds to the smallest volume, additional structures with smaller volumes are created until a minimum is found or the volume is less than or equal to zero
        while (cm == 1 &&  ncm <= USER_MAX_CM_ITER) {
          ncm++;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Strained structure with the smallest volume has the lowest energy" << endl;  
          aus << _AGLSTR_WARNING_ + "Extending number of strained structures to smaller volumes to find minimum" << endl;  
          aus << _AGLSTR_WARNING_ + "Initial structure should be fully relaxed prior to running AFLOW AGL" << endl;
          aus << _AGLSTR_WARNING_ + "Check relaxation of initial structure provided in " << _AFLOWIN_ << " file" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          i--;
          imid = i - mid;
          dimid = imid;
          strainedstructure = initialstructure;
          strainfactor = (dimid * USER_STRAIN_STEP) + 1.0;
          strainfactorlist.push_back(strainfactor);
          if(strainfactor <= 0.0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Strain factor = " << strainfactor << endl;
            aus << _AGLSTR_ERROR_ + "Cannot compress structure any further" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "New strain factor = " << strainfactor << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          strainedstructure.InflateLattice(strainfactor);
          vaspRuns.push_back(vaspRun);
          int idVaspRun = vaspRuns.size()-1;
          // [OBSOLETE] string stridVaspRun;
          ostringstream cnvidVaspRun;
          cnvidVaspRun << idVaspRun;
          string stridVaspRun = cnvidVaspRun.str();
          // [OBSOLETE] string strstrainfactor;
          ostringstream cnvstrainfactor;
          cnvstrainfactor << strainfactor;
          string strstrainfactor = cnvstrainfactor.str();
          // [OBSOLETE] string runname;
          if(USER_DIRNAME_ARUN.option) {
            // [OBSOLETE] runname = _ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor;
            runname.push_back(_ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor);
          } else {
            // [OBSOLETE] runname = _AGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor;
            runname.push_back(_ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor);
          }	  
          arunname.push_back(stridVaspRun + "_SF_" + strstrainfactor);
          dirrunname.push_back(AGL_data.dirpathname + "/" + runname.at(idVaspRun));

          // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AGL_data.dirpathname + "/" + runname;
          // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AGL_data.dirpathname + "/" + runname.at(idVaspRun);
          vaspRuns.at(idVaspRun).str = strainedstructure;
          vaspRuns.at(idVaspRun).AVASP_arun_runname = arunname.at(idVaspRun);
          vaspRuns.at(idVaspRun).AVASP_arun = true;
          vaspRuns.at(idVaspRun).AVASP_arun_mode = "AGL";

          // Assign the values of the flags provided by the user in the _AFLOWIN_ file to the class containing the input data for the VASP run
          aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), vaspFlags, kbinFlags, dirrunname.at(idVaspRun), AGL_data, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to assign values of flags from " << _AFLOWIN_ << " file" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }	  
          avasp_pop_xvasp_success = AVASP_populateXVASP(aflowFlags, kbinFlags, vaspFlags, vaspRuns.at(idVaspRun));
          if (!avasp_pop_xvasp_success) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to set AFLOW flags for ARUN " << arunname.at(idVaspRun) << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }

          // If there are already LOCK or OUTCAR.static files, it means this directory was already generated and computed.
          // Therefore do not touch, but store this structure in the
          // list, so that it can be used in next part of the code.
          if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/"+_AFLOWLOCK_ ) ||
              aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
              aurostd::EFileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ||
              aurostd::EFileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") )) {
            aglerror = AGL_functions::extractenerg(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Failed to extract total energies" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }
            aglerror = AGL_functions::extractedos(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Failed to extract electronic DOS" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }
            aglerror = AGL_functions::checkmin(AGL_data, cm, FileMESSAGE);	 
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Failed to find energy minimum" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }
            continue;
          } 	  
          // Check if structure is on list of failed runs to be skipped
          // If so, then skip generation and continue to next structure
          for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
            // [OBSOLETE] if(aurostd::substring2bool(runname,AGL_data.failed_arun_list.at(ij),TRUE))
            if(aurostd::substring2bool(runname.at(idVaspRun),AGL_data.failed_arun_list.at(ij),TRUE)) {
              skipdir = true;
            }
          }	
          if(skipdir) {
            aurostd::StringstreamClean(aus);
            // [OBSOLETE] aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << runname << endl;
            aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << runname.at(idVaspRun) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            skipdir = false;
            continue;
          }
          // If files do not exist, and the postprocess flag is not set, continue on to prepare generation of _AFLOWIN_ ...
          if (!(XHOST.POSTPROCESS || XHOST.ARUN_POSTPROCESS || AGL_data.postprocess)) { //CO20210701
            // [OBSOLETE] aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname, FileMESSAGE);
            // [OBSOLETE] aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname.at(idVaspRun), FileMESSAGE);
            // [OBSOLETE] if(aglerror != 0) {
            // [OBSOLETE]   aurostd::StringstreamClean(aus);
            // [OBSOLETE]   aus << _AGLSTR_ERROR_ + "Failed to assign flags" << endl;  
            // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            // [OBSOLETE]   return aglerror;
            // [OBSOLETE] }
            // [OBSOLETE]  aglerror = AGL_functions::createAFLOWIN(vaspRuns.at(idVaspRun), xvasp, _kbinFlags, _vaspFlags, AGL_data, FileMESSAGE);
            // [OBSOLETE] if(aglerror != 0) {
            // [OBSOLETE]   aurostd::StringstreamClean(aus);
            // [OBSOLETE]   aus << _AGLSTR_ERROR_ + "Failed to create " << _AFLOWIN_ << " files" << endl;  
            // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            // [OBSOLETE]   return aglerror;
            // [OBSOLETE] }
            aflowin_success = AVASP_MakeSingleAFLOWIN(vaspRuns.at(idVaspRun), arunaflowin, write);
            if(!aflowin_success) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Failed to create " << _AFLOWIN_ << " files" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }
          }
          aglerror = AGL_functions::extractenerg(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to extract total energies" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          aglerror = AGL_functions::extractedos(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to extract electronic DOS" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          aglerror = AGL_functions::checkmin(AGL_data, cm, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to find energy minimum" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;	  
          }	  
        }
        i = vaspRuns.size() - 1;
        // While the lowest energy corresponds to the largest volume, additional structures with larger volumes are created until a minimum is found
        while (cm == 2) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Strained structure with the largest volume has the lowest energy" << endl;  
          aus << _AGLSTR_WARNING_ + "Extending number of strained structures to larger volumes to find minimum" << endl;  
          aus << _AGLSTR_WARNING_ + "Initial structure should be fully relaxed prior to running AFLOW AGL" << endl;
          aus << _AGLSTR_WARNING_ + "Check relaxation of initial structure provided in " << _AFLOWIN_ << " file" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          i++;
          imid = i - mid;
          dimid = imid;
          strainedstructure = initialstructure;
          strainfactor = (dimid * USER_STRAIN_STEP) + 1.0;
          strainfactorlist.push_back(strainfactor);

          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "New strain factor = " << strainfactor << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          strainedstructure.InflateLattice(strainfactor);
          vaspRuns.push_back(vaspRun);
          int idVaspRun = vaspRuns.size()-1;
          // [OBSOLETE] string stridVaspRun;
          ostringstream cnvidVaspRun;
          cnvidVaspRun << idVaspRun;
          string stridVaspRun = cnvidVaspRun.str();
          // [OBSOLETE] string strstrainfactor;
          ostringstream cnvstrainfactor;
          cnvstrainfactor << strainfactor;
          string strstrainfactor = cnvstrainfactor.str();
          // [OBSOLETE] string runname;
          if(USER_DIRNAME_ARUN.option) {
            // [OBSOLETE] runname = _ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor;
            runname.push_back(_ARAGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor);
          } else {
            //OBSOLETE runname = _AGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor;
            runname.push_back(_AGLSTR_DIRNAME_ + stridVaspRun + "_SF_" + strstrainfactor);
          }
          arunname.push_back(stridVaspRun + "_SF_" + strstrainfactor);
          dirrunname.push_back(AGL_data.dirpathname + "/" + runname.at(idVaspRun));

          // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AGL_data.dirpathname + "/" + runname;
          // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AGL_data.dirpathname + "/" + runname.at(idVaspRun);
          vaspRuns.at(idVaspRun).str = strainedstructure;
          vaspRuns.at(idVaspRun).AVASP_arun_runname = arunname.at(idVaspRun);
          vaspRuns.at(idVaspRun).AVASP_arun = true;
          vaspRuns.at(idVaspRun).AVASP_arun_mode = "AGL";

          // Assign the values of the flags provided by the user in the _AFLOWIN_ file to the class containing the input data for the VASP run
          aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), vaspFlags, kbinFlags, dirrunname.at(idVaspRun), AGL_data, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to assign values of flags from " << _AFLOWIN_ << " file" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }	  
          avasp_pop_xvasp_success = AVASP_populateXVASP(aflowFlags, kbinFlags, vaspFlags, vaspRuns.at(idVaspRun));
          if (!avasp_pop_xvasp_success) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to set AFLOW flags for ARUN " << arunname.at(idVaspRun) << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }

          // If there are already LOCK or OUTCAR.static files, it means this directory was already generated and computed.
          // Therefore do not touch, but store this structure in the
          // list, so that it can be used in next part of the code.
          if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/"+_AFLOWLOCK_ ) ||
              aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
              aurostd::EFileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ||
              aurostd::EFileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ) {
            AGL_functions::extractenerg(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
            AGL_functions::extractedos(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
            AGL_functions::checkmin(AGL_data, cm, FileMESSAGE);	 
            continue;
          }
          // Check if structure is on list of failed runs to be skipped
          // If so, then skip generation and continue to next structure
          for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
            // [OBSOLETE] if(aurostd::substring2bool(runname,AGL_data.failed_arun_list.at(ij),TRUE))
            if(aurostd::substring2bool(runname.at(idVaspRun),AGL_data.failed_arun_list.at(ij),TRUE)) {
              skipdir = true;
            }
          }	
          if(skipdir) {
            aurostd::StringstreamClean(aus);
            // [OBSOLETE] aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << runname << endl;
            aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << runname.at(idVaspRun) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            skipdir = false;
            continue;
          } 	  	  
          // If files do not exist, and the postprocess flag is not set, continue on to prepare generation of _AFLOWIN_ ...
          if (!(XHOST.POSTPROCESS || XHOST.ARUN_POSTPROCESS || AGL_data.postprocess)) { //CO20210701
            // [OBSOLETE] aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname, FileMESSAGE);
            // [OBSOLETE] aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname.at(idVaspRun), FileMESSAGE);
            // [OBSOLETE] if(aglerror != 0) {
            // [OBSOLETE]   aurostd::StringstreamClean(aus);
            // [OBSOLETE]   aus << _AGLSTR_ERROR_ + "Failed to assign flags" << endl;  
            // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            // [OBSOLETE]   return aglerror;
            // [OBSOLETE] }
            // [OBSOLETE] aglerror = AGL_functions::createAFLOWIN(vaspRuns.at(idVaspRun), xvasp, _kbinFlags, _vaspFlags, AGL_data, FileMESSAGE);
            // [OBSOLETE] if(aglerror != 0) {
            // [OBSOLETE]   aurostd::StringstreamClean(aus);
            // [OBSOLETE]   aus << _AGLSTR_ERROR_ + "Failed to create " << _AFLOWIN_ << " files" << endl;  
            // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            // [OBSOLETE]   return aglerror;
            // [OBSOLETE] }
            // Writes _AFLOWIN_ file to appropriate directory for each VASP run
            aflowin_success = AVASP_MakeSingleAFLOWIN(vaspRuns.at(idVaspRun), arunaflowin, write);
            if(!aflowin_success) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Failed to create " << _AFLOWIN_ << " files" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }
          }
          aglerror = AGL_functions::extractenerg(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to extract total energies" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          aglerror = AGL_functions::extractedos(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to extract electronic DOS" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          aglerror = AGL_functions::checkmin(AGL_data, cm, FileMESSAGE);	  
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to find energy minimum" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
        }
      }

      // Checks if the (E, V) data is concave, which is required by the GIBBS algorithm
      // If there is no concave part around the minimum of the data passed to it, the GIBBS algorithm will have problems fitting the data by a polynomial, and will exit with an error
      // Problems with the concavity of the (E, V) data are usually caused by an insufficient basis set or k-point mesh
      // If the data is not concave, then this section of the code will increase the density of the k-point mesh and rerun the VASP calculations
      // If this part of the code produces warnings, you should check the k-point mesh and the basis set
      if(USER_CHECK_EV_CONCAVITY.option && !(XHOST.ARUN_POSTPROCESS || AGL_data.postprocess)) {	  
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Checking concavity of (E, V) data" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // The function "checkconcav" checks the concavity of the (E, V) data
        // If the data is concave, it sets cc = 0
        // If the data is not concave, it sets cc = 1
        // Initialize to 0 so that default assumption is that data is concave to avoid running unnecessary additional calculations
        int cc = 0;
        aglerror = AGL_functions::checkconcav(AGL_data, cc, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failure checking energy concavity" << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        int ncc = 0;
        while (cc == 1 && ncc <= USER_MAX_CC_ITER) {
          ncc++;
          for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
            string runnamedir = vaspRuns.at(idVaspRun).Directory;
            aglerror = AGL_functions::aglvaspflags(vaspRuns.at(idVaspRun), vaspFlags, kbinFlags, runnamedir, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Failed to assign flags" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }
            vaspRuns.at(idVaspRun).AVASP_value_KPPRA = vaspRuns.at(idVaspRun).AVASP_value_KPPRA * 2;
            // [OBSOLETE] string strkppra;
            ostringstream cnvkppra;
            cnvkppra << vaspRuns.at(idVaspRun).AVASP_value_KPPRA;
            string strkppra = cnvkppra.str();
            vaspRuns.at(idVaspRun).Directory = vaspRuns.at(idVaspRun).Directory + "_K_" + strkppra;	    
          }
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "(E, V) data is not concave" << endl;
          aus << _AGLSTR_WARNING_ + "Increasing the number of k-points and rerunning VASP calculations" << endl;
          aus << _AGLSTR_WARNING_ + "Resetting value of KPPRA from " << vaspFlags.KBIN_VASP_KPOINTS_KPPRA.content_int << " to " << vaspRuns.at(0).AVASP_value_KPPRA << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE]for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
          // [OBSOLETE]  aglerror = AGL_functions::createAFLOWIN(vaspRuns.at(idVaspRun), xvasp, _kbinFlags, _vaspFlags, AGL_data, FileMESSAGE);
          // [OBSOLETE]  if(aglerror != 0) {
          // [OBSOLETE]    aurostd::StringstreamClean(aus);
          // [OBSOLETE]    aus << _AGLSTR_ERROR_ + "Failed to create " << _AFLOWIN_ << " files" << endl;  
          // [OBSOLETE]    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE]    return aglerror;
          // [OBSOLETE]  }
          // [OBSOLETE]}
          aglerror = AGL_functions::extractenerg(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to extract total energies" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          aglerror = AGL_functions::extractedos(vaspRuns, AGL_data, dirrunname, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failed to extract electronic DOS" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          aglerror = AGL_functions::checkconcav(AGL_data, cc, FileMESSAGE);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failure checking energy concavity" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
        }
        if(cc == 1 && ncc >= USER_MAX_CC_ITER) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Maximum number of concavity check iterations exceeded without producing concave data" << endl;
          aus << _AGLSTR_WARNING_ + "Check settings for DFT calculations such as k-point mesh and basis set" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }

      // Checks that number of energy values and volume values is the same
      if(AGL_data.volumeinput.size() != AGL_data.energyinput.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Number of energy and volume points are different" << endl; 
        aus << _AGLSTR_ERROR_ + "Number of volume points = " << AGL_data.volumeinput.size() << endl; 
        aus << _AGLSTR_ERROR_ + "Number of energy points = " << AGL_data.energyinput.size() << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }

      // Checks that number of pressure values and volume values is the same
      if(AGL_data.pressurecalculated.size() != AGL_data.volumeinput.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Number of pressure and volume points are different" << endl; 
        aus << _AGLSTR_ERROR_ + "Number of volume points = " << AGL_data.volumeinput.size() << endl; 
        aus << _AGLSTR_ERROR_ + "Number of pressure points = " << AGL_data.pressurecalculated.size() << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }

      // Checks that number of skipped points is less than the maximum permitted value
      if((vaspRuns.size() - AGL_data.volumeinput.size()) > AGL_data.skiparunsmax) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Number of skipped runs exceeds maximum permitted value of " << AGL_data.skiparunsmax << endl; 
        aus << _AGLSTR_ERROR_ + "Number of accepted runs = " << AGL_data.volumeinput.size() << endl; 
        aus << _AGLSTR_ERROR_ + "Number of initial runs = " << vaspRuns.size() << endl; 
        aus << _AGLSTR_ERROR_ + "Number of skipped runs = " << (vaspRuns.size() - AGL_data.volumeinput.size()) << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;	
      }

      // Adds (E, V) data to AGL_data class
      // [OBSOLETE] Converts units of energy values calculated by VASP from eV to Hartree
      // [OBSOLETE] Converts units of volume from cubic Angstrom to cubic Bohr
      // Converts units of pressure from kB to GPa
      // [OBSOLETE] AGL_data.volumeinput.resize(vaspRuns.size());
      // [OBSOLETE] for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++)
      // Vectors to store original E(V) data in case GIBBS needs to be re-run
      AGL_data.energyinput_orig.resize(AGL_data.volumeinput.size());
      AGL_data.volumeinput_orig.resize(AGL_data.volumeinput.size());
      for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
        // [OBSOLETE] AGL_data.volumeinput.at(idVaspRun) = vaspRuns.at(idVaspRun).str.Volume() * pow(angstrom2bohr, 3.0);
        // [OBSOLETE] AGL_data.energyinput.at(idVaspRun) *= ev2hart1;
        // [OBSOLETE] AGL_data.energyinput.at(idVaspRun) = AGL_data.energyinput.at(idVaspRun) / hart2ev;
        // [OBSOLETE] AGL_data.volumeinput.at(i) = AGL_data.volumeinput.at(i) * pow(angstrom2bohr, 3.0);
        // [OBSOLETE] AGL_data.energyinput.at(i) = AGL_data.energyinput.at(i) / hart2ev;
        // Store original E(V) data in case GIBBS needs to be re-run
        AGL_data.volumeinput_orig.at(i) = AGL_data.volumeinput.at(i);
        AGL_data.energyinput_orig.at(i) = AGL_data.energyinput.at(i);
        // Convert pressure units from kB to GPa
        AGL_data.pressurecalculated.at(i) = AGL_data.pressurecalculated.at(i) / 10.0;
        aurostd::StringstreamClean(aus);
        // Print out total energy, volume and calculated residual pressure
        aus << _AGLSTR_MESSAGE_ << "Energy (eV) = " << AGL_data.energyinput.at(i) << endl;
        aus << _AGLSTR_MESSAGE_ << "Volume (Ang^3) = " << AGL_data.volumeinput.at(i) << endl;
        aus << _AGLSTR_MESSAGE_ << "Pressure (GPa) = " << AGL_data.pressurecalculated.at(i) << endl;      
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }

      // Find maximum pressure for VASP calculations
      AGL_data.pressurecalculatedmax = AGL_data.pressurecalculated.at(0);
      for (uint i = 0; i < AGL_data.pressurecalculated.size(); i++) {
        if (AGL_data.pressurecalculatedmax < AGL_data.pressurecalculated.at(i)) {
          AGL_data.pressurecalculatedmax = AGL_data.pressurecalculated.at(i);
        }
      }

      // Writes input file for original, non-AFLOW version of GIBBS (useful for testing and debugging)
      if(USER_WRITE_GIBBS_INPUT.option) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Writing input file for original, non-AFLOW version of GIBBS" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = AGL_functions::gibbsinpwrite(AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to write GIBBS input file" << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
      }

      if(USER_WRITE_FULL_RESULTS.option) {
        // Write structures to JSON format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_structures.json" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        string strjson = xstructure2json(xvasp.str);
        oss << "{\"AGL_structures\":[";
        oss << "{\"distorted_structure\":[{\"distortion\":null},"; 
        oss << "{\"structure\":" << strjson << "}]}";
        for (uint i = 0; i < vaspRuns.size(); i++) {
          strjson = xstructure2json(vaspRuns.at(i).str);
          oss << ",{\"distorted_structure\":[{\"distortion\":\"" << runname.at(i) << "\"},";
          oss << "{\"structure\":" << strjson << "}]}";
        }
        oss << "]}";
        string ofilestrjson = AGL_data.dirpathname + "/AGL_structures.json";
        if(!aurostd::stringstream2file(oss, ofilestrjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_structure.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
      }

      // Write structures and corresponding energies, pressures and stress tensors to JSON format
      aurostd::xmatrix<double> stress_tensor(3, 3);      
      // First check that there is something to write
      if (AGL_data.structurecalculated.size() < 1) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "No runs to write out yet" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 8;
      } else {
        // [OBSOLETE] aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "Opening file AGL_energy_structures.json" <<  endl;
        // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // Only write structures that have been successfully calculated
        uint idSuccessRun = AGL_data.structurecalculated.at(0); 
        string strjson = xstructure2json(vaspRuns.at(idSuccessRun).str);
        oss << "{\"AGL_energy_structures\":[";
        oss << "{\"distortion\":\"" << runname.at(idSuccessRun) << "\",";
        oss << "\"strainfactor\":" << strainfactorlist.at(idSuccessRun) << ",";
        oss << "\"energy\":" << AGL_data.energyinput.at(0) / AGL_data.natoms << ",";
        oss << "\"pressure\":" << AGL_data.pressurecalculated.at(0) << ",";
        stress_tensor = AGL_data.stresscalculated.at(0);
        // [OBSOLETE] aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "stress_tensor = " << stress_tensor << endl;
        // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        oss << "\"stress_tensor\":[";
        for (uint j = 1; j < 4; j++) {
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "j = " << j << endl;
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if (j == 1) {
            oss << "[";
          } else {
            oss << ",[";
          }
          for (uint k = 1; k < 4; k++) {
            // [OBSOLETE] aurostd::StringstreamClean(aus);
            // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "k = " << k << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if (k == 1) {
              oss << 0.1 * stress_tensor[j][k];
            } else {
              oss << "," << 0.1 * stress_tensor[j][k];
            }
          }
          oss << "]";	
        }
        oss << "],";
        oss << "\"structure\":" << strjson << "}";
        for (uint i = 1; i < AGL_data.structurecalculated.size(); i++) {
          // Only write structures that have been successfully calculated	  
          idSuccessRun = AGL_data.structurecalculated.at(i); 
          strjson = xstructure2json(vaspRuns.at(idSuccessRun).str);
          oss << ",{\"distortion\":\"" << runname.at(idSuccessRun) << "\",";
          oss << "\"strainfactor\":" << strainfactorlist.at(idSuccessRun) << ",";	
          oss << "\"energy\":" << AGL_data.energyinput.at(i) / AGL_data.natoms << ",";
          oss << "\"pressure\":" << AGL_data.pressurecalculated.at(i) << ",";
          stress_tensor = AGL_data.stresscalculated.at(i);
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "stress_tensor = " << stress_tensor << endl;
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          oss << "\"stress_tensor\":[";
          for (uint j = 1; j < 4; j++) {
            if (j == 1) {
              oss << "[";
            } else {
              oss << ",[";
            }
            for (uint k = 1; k < 4; k++) {
              if (k == 1) {
                oss << 0.1 * stress_tensor[j][k];
              } else {
                oss << "," << 0.1 * stress_tensor[j][k];
              }
            }
            oss << "]";	
          }
          oss << "],";
          oss << "\"structure\":" << strjson << "}";
        }
        oss << "]}";
        string ofileenergstrjson = AGL_data.dirpathname + "/AGL_energy_structures.json";
        if(!aurostd::stringstream2file(oss, ofileenergstrjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energy_structure.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
      }

      // Clear vaspRuns vector to free up memory
      vaspRuns.clear();

      // Runs GIBBS algorithm within AFLOW to calculate thermal properties
      aglerror = AGL_functions::gibbsrun(AGL_data, FileMESSAGE);

      // Opening file to write thermal properties in plottable format and write stringstream to it
      // Similar format to the "THERMO" file written by APL
      string outfiletname = AGL_data.dirpathname + "/AGL_THERMO";
      if(AGL_data.i_eqn_of_state >= 0) {
        if(!aurostd::stringstream2file(AGL_data.outfiletss, outfiletname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_THERMO.dat" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
      }
      aurostd::StringstreamClean(AGL_data.outfiletss);

      // Open output file for GIBBS output data for all p, T values and write stringstream to it
      // This file is in the same format as the original Fortran version of the GIBBS program output file
      string outputfilename = AGL_data.dirpathname + "/AGL.out";
      if(USER_SYSOUT_FILENAME.option) {
        string outputfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + ".out";
      } 
      if(!aurostd::stringstream2file(AGL_data.outfiless, outputfilename, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << outputfilename.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(AGL_data.outfiless);

      // Checks to see if GIBBS algorithm ran successfully
      // A value of 0 indicates successful completion
      // A value of 1 indicates an error opening a file
      // A value of 2 indicates an error finding the minimum of the (E, V) (perhaps at non-zero temperature or pressure, see previously written warnings)
      // A value of 3 indicates an error with the vector sizes
      // A value of 4 indicates an error in inverting a matrix to obtain a polynomial fit 
      // Another, non-zero value indicates some other error, see previously written warnings
      if(aglerror != 0) {	  
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "gibbsrun failed" << endl;
        if(aglerror == 1) {
          aus << _AGLSTR_ERROR_ + "error opening a file" << endl;
        }
        if(aglerror == 2) {
          aus << _AGLSTR_ERROR_ + "Problem finding minimum of (E, V) data calculated with VASP" << endl;
          aus << _AGLSTR_ERROR_ + "Try increasing the number of k-points or the number of strained structures and rerunning VASP calculations" << endl;
          if(AGL_data.EV_noise) {
            aus << _AGLSTR_ERROR_ + "Noise in E-V data caused run failure: re-run with increased STATIC_KPPRA" << endl;
          } else if(abs(AGL_data.itdiff) > 1)  {
            aus << _AGLSTR_ERROR_ + "Minimum in E-V data not at center of E-V data: check relaxation of initial structure" << endl;	    
          }
        }
        if(aglerror == 3) {
          aus << _AGLSTR_ERROR_ + "Error in vector sizes" << endl;
        }
        if(aglerror == 4) {
          aus << _AGLSTR_ERROR_ + "Problem inverting matrix to find polynomial" << endl;
        }
        aus << _AGLSTR_ERROR_ + "Exiting GIBBS" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aglerror;
      }

      // Calculate volume and energy per atom for calculated E(V) data
      // Include minimum energy structure from EoS fit
      vector<double> volumetotal;
      vector<double> energytotal;
      vector<double> volume_atom;
      vector<double> energy_atom;
      volumetotal = AGL_data.volumeinput;
      energytotal = AGL_data.energyinput;
      volumetotal.push_back(AGL_data.volume_equilibrium_0p0T);
      energytotal.push_back(AGL_data.energy_equilibrium_0p0T);
      aglerror = AGL_functions::qcksortev(volumetotal, energytotal, FileMESSAGE);
      for (uint i = 0; i < volumetotal.size(); i++) {
        volume_atom.push_back(volumetotal.at(i) / AGL_data.natoms);
        energy_atom.push_back(energytotal.at(i) / AGL_data.natoms);
      }

      // Write calculated E(V) data including equilibrium volume determined by fitting calculated E(V) data
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Opening file AGL_energy_volume.out" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      oss << "# E(V) data including equilibrium volume determined by fitting calculated E(V) data" << endl;
      oss << "# Equilibrium volume at zero temperature and pressure = " << AGL_data.volume_equilibrium_0p0T << " Angstrom^3/cell" << endl;
      oss << "# Energy at equilibrium volume at zero temperature and pressure = " << AGL_data.energy_equilibrium_0p0T << " eV/cell" << endl;
      oss << "# Equilibrium volume at zero temperature and pressure = " << AGL_data.volume_equilibrium_0p0T / AGL_data.natoms << " Angstrom^3/atom" << endl;
      oss << "# Energy at equilibrium volume at zero temperature and pressure = " << AGL_data.energy_equilibrium_0p0T / AGL_data.natoms << " eV/atom" << endl;
      oss << "# Volume (Ang^3/atom)        Energy(eV/atom)        Volume (Ang^3/cell)        Energy(eV/cell)" << endl;
      for (uint i = 0; i < volumetotal.size(); i++) {
        oss << setw(21) << setprecision(6) << volume_atom.at(i) << "\t" << setw(20) << setprecision(6) << energy_atom.at(i) << "\t" << setw(23) << setprecision(6) << volumetotal.at(i) << "\t" << setw(22) << setprecision(6) << energytotal.at(i) << endl;
      }
      string ofileenergyvolume = AGL_data.dirpathname + "/AGL_energy_volume.out";
      if(!aurostd::stringstream2file(oss, ofileenergyvolume, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energy_volume.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      oss.clear();
      oss.str(std::string());

      // Calculate and write electronic properties as a function of pressure in JSON and plain text format
      aglerror = AGL_functions::edosbandgap(AGL_data, FileMESSAGE);
      if(aglerror != 0 && aglerror != 2) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failure in calculation of band gap as a function of pressure" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aglerror;
      } else {
        // Write calculated band gap as a function of pressure
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_edos_gap_pressure.json" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        oss << "{\"edos_gap_pressure_derivative\": " << AGL_data.edosgap_pressure;
        oss << ",\"edos_value_EF_pressure_derivative\": " << AGL_data.dosvalEF_pressure;
        oss << ",\"edos_gap_minimum\": " << AGL_data.egap_min;
        oss << ",\"edos_gap_maximum\": " << AGL_data.egap_max;
        oss << ",\"edos_gap_minimum_pressure\": " << AGL_data.egap_min_pressure;
        oss << ",\"edos_gap_maximum_pressure\": " << AGL_data.egap_max_pressure;
        oss << ",\"edos_value_EF_minimum\": " << AGL_data.dosvalEF_min;
        oss << ",\"edos_value_EF_maximum\": " << AGL_data.dosvalEF_max;
        oss << ",\"edos_value_EF_minimum_pressure\": " << AGL_data.dosvalEF_min_pressure;
        oss << ",\"edos_value_EF_maximum_pressure\": " << AGL_data.dosvalEF_max_pressure;
        oss << ",\"edos_gap_pressure\": [";
        bool firstpres = true;
        for (uint i = 0; i < AGL_data.AGL_edos_properties.size(); i++) {
          if (AGL_data.AGL_edos_properties.at(i).edos_gap_set && AGL_data.AGL_edos_properties.at(i).dosval_EF_set && AGL_data.AGL_edos_properties.at(i).edos_gap_poly_set && AGL_data.AGL_edos_properties.at(i).dosval_EF_poly_set) {
            if (firstpres) {
              if (!std::isnan(AGL_data.AGL_edos_properties.at(i).edos_band_gap) && !std::isnan(AGL_data.AGL_edos_properties.at(i).dosval_EF) && !std::isnan(AGL_data.AGL_edos_properties.at(i).edos_band_gap_polyval) && !std::isnan(AGL_data.AGL_edos_properties.at(i).dosval_EF_polyval)) {
                oss << "{\"pressure\":" << AGL_data.AGL_edos_properties.at(i).pressure_external << ",\"edos_band_gap\":" << AGL_data.AGL_edos_properties.at(i).edos_band_gap << ",\"edos_value_EF\":" << AGL_data.AGL_edos_properties.at(i).dosval_EF << ",\"edos_band_gap_polyval\":" << AGL_data.AGL_edos_properties.at(i).edos_band_gap_polyval << ",\"edos_value_EF_polyval\":" << AGL_data.AGL_edos_properties.at(i).dosval_EF_polyval << "}";
                firstpres = false;
              }
            } else {
              if (!std::isnan(AGL_data.AGL_edos_properties.at(i).edos_band_gap) && !std::isnan(AGL_data.AGL_edos_properties.at(i).dosval_EF) && !std::isnan(AGL_data.AGL_edos_properties.at(i).edos_band_gap_polyval) && !std::isnan(AGL_data.AGL_edos_properties.at(i).dosval_EF_polyval)) {
                oss << ",{\"pressure\":" << AGL_data.AGL_edos_properties.at(i).pressure_external << ",\"edos_band_gap\":" << AGL_data.AGL_edos_properties.at(i).edos_band_gap << ",\"edos_value_EF\":" << AGL_data.AGL_edos_properties.at(i).dosval_EF << ",\"edos_band_gap_polyval\":" << AGL_data.AGL_edos_properties.at(i).edos_band_gap_polyval << ",\"edos_value_EF_polyval\":" << AGL_data.AGL_edos_properties.at(i).dosval_EF_polyval << "}";
              }
            }
          }
        }      	
        oss << "],\"units\":{\"edos_band_gap\":\"eV\",\"pressure\":\"GPa\",\"edos_value_EF\":\"states/(cell*eV)\",\"edos_gap_pressure_derivative\":\"eV/GPa\",\"edos_value_EF_pressure_derivative\":\"(states/(cell*eV))/GPa\"}";
        oss << "}" << endl;
        string ofilepressurebandgapjson = AGL_data.dirpathname + "/AGL_edos_gap_pressure.json";
        if(!aurostd::stringstream2file(oss, ofilepressurebandgapjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_edos_gap_pressure.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());	
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_edos_gap_pressure.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        oss << "# Calculated electronic DOS band gap as a function of pressure" << endl;
        oss << "# Derivative of pressure-band gap at zero pressure = " << AGL_data.edosgap_pressure << " eV/GPa" << endl;
        oss << "# Derivative of pressure-DOS(EF) at zero pressure = " << AGL_data.dosvalEF_pressure << " (states/(cell*eV))/GPa" << endl;
        oss << "# Minimum value of electronic band gap = " << AGL_data.egap_min << " eV" << endl;
        oss << "# Maximum value of electronic band gap = " << AGL_data.egap_max << " eV" << endl;
        oss << "# Minimum value of electronic band gap is at pressure = " << AGL_data.egap_min_pressure << " GPa" << endl;
        oss << "# Maximum value of electronic band gap is at pressure = " << AGL_data.egap_max_pressure << " GPa" << endl;
        oss << "# Minimum value of electronic DOS at Fermi energy = " << AGL_data.dosvalEF_min << " eV" << endl;
        oss << "# Maximum value of electronic DOS at Fermi energy = " << AGL_data.dosvalEF_max << " eV" << endl;
        oss << "# Minimum value of electronic DOS at Fermi energy is at pressure = " << AGL_data.dosvalEF_min_pressure << " GPa" << endl;
        oss << "# Maximum value of electronic DOS at Fermi energy is at pressure = " << AGL_data.dosvalEF_max_pressure << " GPa" << endl;
        oss << "# Pressure (GPa)          Electronic Band Gap (eV)          DOS at EF (states/(cell*eV))          Electronic Band Gap Polynomial Fit (eV)          DOS at EF Polynomial Fit (states/(cell*eV))" << endl;
        for (uint i = 0; i < AGL_data.AGL_edos_properties.size(); i++) {
          if (AGL_data.AGL_edos_properties.at(i).edos_gap_set && AGL_data.AGL_edos_properties.at(i).dosval_EF_set && AGL_data.AGL_edos_properties.at(i).edos_gap_poly_set && AGL_data.AGL_edos_properties.at(i).dosval_EF_poly_set) {
            if (!std::isnan(AGL_data.AGL_edos_properties.at(i).edos_band_gap) && !std::isnan(AGL_data.AGL_edos_properties.at(i).dosval_EF) && !std::isnan(AGL_data.AGL_edos_properties.at(i).edos_band_gap_polyval) && !std::isnan(AGL_data.AGL_edos_properties.at(i).dosval_EF_polyval)) {
              oss << setw(16) << setprecision(6) << AGL_data.AGL_edos_properties.at(i).pressure_external << "\t" << setw(26) << setprecision(6) << AGL_data.AGL_edos_properties.at(i).edos_band_gap  << "\t" << setw(32) << setprecision(6) << AGL_data.AGL_edos_properties.at(i).dosval_EF << "\t" << setw(41) << setprecision(6) << AGL_data.AGL_edos_properties.at(i).edos_band_gap_polyval  << "\t" << setw(46) << setprecision(6) << AGL_data.AGL_edos_properties.at(i).dosval_EF_polyval << endl;
            }
          }
        }
        string ofilepressurebandgap = AGL_data.dirpathname + "/AGL_edos_gap_pressure.out";
        if(!aurostd::stringstream2file(oss, ofilepressurebandgap, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_edos_gap_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
      }

      // Container vector for Hugoniot output
      vector<vector<double> > hugoniot_output;
      // If AGL_data.run_all_pressure_temperature is selected, writes out pressure and temperature properties from list of structs
      // Exits without writing other output since this output might not be available for all pressures and temperatures in non-truncated ranges
      if (AGL_data.run_all_pressure_temperature) {
        // Run Hugoniot calculation
        if(USER_HUGONIOT_CALC.option) {
          // [OBSOLETE] vector<vector<double> > hugoniot_output;
          double cellmass_grams = cellmasskg * 1000.0;
          aglerror = AGL_functions::runHugoniotAllTemperaturesAndPressures(AGL_data.AGL_pressure_temperature_energy_list, cellmass_grams, hugoniot_output, AGL_data.hugoniotextrapolate, FileMESSAGE);
          if(aglerror != 0 && aglerror != 2) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failure in Hugoniot calculation" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          } else {
            oss << std::setw(15) << "#     rho(g/cc)" << std::setw(15) << "T(K)" << std::setw(15) << "E(kJ/g)" << std::setw(15) << "P(GPa)" << endl;
            for (uint i = 0; i < hugoniot_output.size(); i++) {
              for (uint j = 0; j < hugoniot_output.at(i).size(); j++) {
                // [OBSOLETE] oss << "\t" << std::left << std::setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
                oss << setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
              }
              oss << endl;
            }
            string ofilehugoniotoutput = AGL_data.dirpathname + "/AGL_Hugoniot.out";
            if(!aurostd::stringstream2file(oss, ofilehugoniotoutput, "WRITE")) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Hugoniot.out" <<  endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }
            // Reset aglerror to 0 on successful output of Hugoniot data
            aglerror = 0;
            // [OBSOLETE] aurostd::StringstreamClean(oss);
          }
          oss.clear();
          oss.str(std::string());
        }
        // Write calculated energies to JSON format:
        // E(V) data including equilibrium volume determined by fitting calculated E(V) data
        // Gibbs free energy as a function of temperature and pressure
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_energy.json" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        oss << "{\"energy_volume_atom\": [";
        for (uint i = 0; i < volume_atom.size(); i++) {
          if (i==0) {
            oss << "{\"energy\":" << energy_atom.at(i) << ",\"volume\":" << volume_atom.at(i) << "}";
          } else {
            oss << ",{\"energy\":" << energy_atom.at(i) << ",\"volume\":" << volume_atom.at(i) << "}";
          }
        }
        oss << "],\"enthalpy_pressure_atom\": [";
        for (uint k = 0; k < AGL_data.AGL_pressure_enthalpy_list.size(); k++) {
          if (!std::isnan(AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy)) {
            firstnonnan = k;
            break;
          }
        }
        oss << "{\"pressure\":" << AGL_data.AGL_pressure_enthalpy_list.at(firstnonnan).pressure_external << ",\"enthalpy\":" << AGL_data.AGL_pressure_enthalpy_list.at(firstnonnan).enthalpy / AGL_data.natoms << "}";
        for (uint k = firstnonnan+1; k < AGL_data.AGL_pressure_enthalpy_list.size(); k++) {	
          if (!std::isnan(AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy)) {
            oss << ",{\"pressure\":" << AGL_data.AGL_pressure_enthalpy_list.at(k).pressure_external << ",\"enthalpy\":" << AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy / AGL_data.natoms << "}";
          }
        }      	
        oss << "],\"gibbs_free_energy_atom\": [";
        for (uint i = 0; i < AGL_data.AGL_pressure_temperature_energy_list.size(); i++) {
          if (!std::isnan(AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy)) {
            firstnonnan = i;
            break;
          }
        }
        oss << "{\"temperature\":" << AGL_data.AGL_pressure_temperature_energy_list.at(firstnonnan).temperature_external << ",\"pressure\":" << AGL_data.AGL_pressure_temperature_energy_list.at(firstnonnan).pressure_external << ",\"gibbs_energy\":" << AGL_data.AGL_pressure_temperature_energy_list.at(firstnonnan).gibbs_free_energy / AGL_data.natoms << "}";
        for (uint i = firstnonnan+1; i < AGL_data.AGL_pressure_temperature_energy_list.size(); i++) {
          if (!std::isnan(AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy)) {
            oss << ",{\"temperature\":" << AGL_data.AGL_pressure_temperature_energy_list.at(i).temperature_external << ",\"pressure\":" << AGL_data.AGL_pressure_temperature_energy_list.at(i).pressure_external << ",\"gibbs_energy\":" << AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy / AGL_data.natoms << "}";
          }
        }
        oss << "],\"temperature_stepnumber\":" << USER_NTEMPERATURE << ",\"temperature_stepsize\":" << USER_STEMPERATURE << ",\"pressure_stepnumber\":" << USER_NPRESSURE << ",\"pressure_stepsize\":" << USER_SPRESSURE;
        oss << ",\"units\":{\"energy\":\"eV/atom\",\"volume\":\"Angstrom^3/atom\",\"temperature\":\"K\",\"pressure\":\"GPa\"}";
        if (USER_HUGONIOT_CALC.option) {
          oss << ",\"hugoniot_output\": [";
          for (uint i = 0; i < hugoniot_output.size(); i++) {
            if ((hugoniot_output.at(i).size() < 4) || std::isnan(hugoniot_output.at(i).at(0)) || std::isnan(hugoniot_output.at(i).at(1)) || std::isnan(hugoniot_output.at(i).at(2)) || std::isnan(hugoniot_output.at(i).at(3))) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ << "Hugoniot entry " << i << " only has " << hugoniot_output.at(i).size() << " components" <<  endl;
              aus << _AGLSTR_WARNING_ << "It should have 4 components, not writing this entry to JSON" <<  endl;	      
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              continue;
            } else {
              firstnonnan = i;
              break;
            }
          }
          oss << "{\"mass_density\":" << hugoniot_output.at(firstnonnan).at(0) << ",\"temperature\":" << hugoniot_output.at(firstnonnan).at(1) << ",\"energy\":" << hugoniot_output.at(firstnonnan).at(2) << ",\"pressure\":"  << hugoniot_output.at(firstnonnan).at(3) << "}";
          for (uint i = firstnonnan + 1; i < hugoniot_output.size(); i++) {
            oss << ",{\"mass_density\":" << hugoniot_output.at(i).at(0) << ",\"temperature\":" << hugoniot_output.at(i).at(1) << ",\"energy\":" << hugoniot_output.at(i).at(2) << ",\"pressure\":"  << hugoniot_output.at(i).at(3) << "}";
          }
          oss << "],\"hugoniot_units\":{\"mass_density\":\"g/cc\",\"energy\":\"kJ/g\",\"temperature\":\"K\",\"pressure\":\"GPa\"}";
        }
        oss << "}" << endl;
        string ofileenergyjson = AGL_data.dirpathname + "/AGL_energy.json";
        if(!aurostd::stringstream2file(oss, ofileenergyjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energy.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());	
        oss << "# temperature_stepnumber=" << USER_NTEMPERATURE << endl;
        oss << "# temperature_stepsize=" << USER_STEMPERATURE << endl;
        oss << "# pressure_stepnumber=" << USER_NPRESSURE << endl;
        oss << "# pressure_stepsize=" << USER_SPRESSURE << endl;
        oss << "#  T (K) " << "        " << aurostd::PaddedPRE("P (GPa)",12," ") << "        " << aurostd::PaddedPRE("G (eV/atom)",13," ") << endl;
        for (uint i = 0; i < AGL_data.AGL_pressure_temperature_energy_list.size(); i++) {
          if (!std::isnan(AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy)) {
            oss << setw(8) << setprecision(2) << fixed << AGL_data.AGL_pressure_temperature_energy_list.at(i).temperature_external << "          " << setw(11) << setprecision(2) << AGL_data.AGL_pressure_temperature_energy_list.at(i).pressure_external << "          " << setw(11) << setprecision(6) << AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy / AGL_data.natoms << endl;
          }
        }	
        string ofileenergies = AGL_data.dirpathname + "/AGL_gibbs_energy_pT.out";
        if(!aurostd::stringstream2file(oss, ofileenergies, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_gibbs_energy_pT.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        oss.clear();
        oss.str(std::string());
        oss << "# pressure_stepnumber=" << USER_NPRESSURE << endl;
        oss << "# pressure_stepsize=" << USER_SPRESSURE << endl;
        oss << "#  P (GPa) " << "        " << aurostd::PaddedPRE("H (eV/atom)",13," ") << endl;
        for (uint k = 0; k < AGL_data.AGL_pressure_enthalpy_list.size(); k++) {
          if (!std::isnan(AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy)) {
            oss << setw(10) << setprecision(2) << fixed << setprecision(2) << AGL_data.AGL_pressure_enthalpy_list.at(k).pressure_external << "          " << setw(12) << setprecision(6) << AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy / AGL_data.natoms << endl;
          }
        }	
        string ofileenthalpies = AGL_data.dirpathname + "/AGL_enthalpy_pressure.out";
        if(!aurostd::stringstream2file(oss, ofileenthalpies, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_enthalpy_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        oss.clear();
        oss.str(std::string());		
        // Exits AGL without writing other output since this output might not be available for all pressures and temperatures in non-truncated ranges
        return 0;
      }		

      // [OBSOLETE] // Converts units of volume as a function of pressure from cubic Bohr to cubic Angstrom 
      // [OBSOLETE] for (uint i = 0; i < AGL_data.VolumeStaticPressure.size(); i++) {
      // [OBSOLETE] AGL_data.VolumeStaticPressure.at(i) = AGL_data.VolumeStaticPressure.at(i) * pow(bohr2angstrom, 3.0);
      // [OBSOLETE] }

      // Checks that number of volumes as a function of pressure is equal to number of pressures
      if(AGL_data.VolumeStaticPressure.size() == AGL_data.StaticPressure.size()) {
        // Prints volume as a function of pressure
        double volumefactor = 1.0;
        for (uint i = 0; i < AGL_data.VolumeStaticPressure.size(); i++) {
          volumefactor = AGL_data.VolumeStaticPressure.at(i) / initialvolume;
          AGL_data.VolumeFactors.push_back(volumefactor);
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "Pressure = " << AGL_data.StaticPressure.at(i) << ", Volume = " << AGL_data.VolumeStaticPressure.at(i) << ", Volume scale factor = " << AGL_data.VolumeFactors.at(i) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else {
        aus << _AGLSTR_ERROR_ + "Number of volumes as a function of pressure is not equal to the number of pressures" << endl;
      }

      // If USER_IDEBYE is negative, only a static calculation has been performed, and no finite temperature calculations have been done
      // AGL exits at this point, writing a message for the user to let them know that only a static calculation was performed
      if(USER_IDEBYE < 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Static GIBBS calculation only requested in " << _AFLOWIN_ << " file" << endl;  
        aus << _AGLSTR_MESSAGE_ << "[AFLOW_GIBBS]IDEBYE=" << USER_IDEBYE << endl;
        aus << _AGLSTR_MESSAGE_ << "To calculate thermal properties, increase value of IDEBYE to greater than or equal to zero in " << _AFLOWIN_ << " file" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 0;
      }

      // If USER_NTEMPERATURE < 1, only a static calculation has been performed, and no finite temperature calculations have been done
      // AFLOW AGL exits at this point, writing a message for the user to let them know that only a static calculation was performed
      if(USER_NTEMPERATURE < 1) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Static GIBBS calculation only requested in " << _AFLOWIN_ << " file" << endl;  
        aus << _AGLSTR_MESSAGE_ << "[AFLOW_GIBBS]NTEMP=" << USER_NTEMPERATURE << endl;
        aus << _AGLSTR_MESSAGE_ << "To calculate thermal properties, increase value of NTEMP to greater than or equal to one in " << _AFLOWIN_ << " file" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 0;
      }

      // If USER_IEOS is negative, only minimal output has been provided, and several thermal properties have not been calculated
      // AFLOW AGL exits at this point, writing a message for the user to let them know that minimal output has been written
      if(USER_IEOS < 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Minimal GIBBS output requested in " << _AFLOWIN_ << " file" << endl;  
        aus << _AGLSTR_MESSAGE_ << "[AFLOW_GIBBS]IEOS=" << USER_IEOS << endl;
        aus << _AGLSTR_MESSAGE_ << "To get more output, increase value of IEOS to greater than or equal to zero in " << _AFLOWIN_ << " file" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 0;
      }

      // Calculates mass density as a function of pressure and temperature in units of grams/cm^3
      AGL_data.mass_density_gcm3.resize(AGL_data.temperature_external.size());
      // First get mass of cell in units of grams
      cellmass_grams = cellmasskg * 1000.0;
      for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          // Next get equilibrium cell volume in units of cm^3
          // [OBSOLETE] volume_cm3 = (AGL_data.VolumeEquilibrium.at(j).at(k) / pow(angstrom2bohr, 3.0)) * 1e-24;
          volume_cm3 = AGL_data.VolumeEquilibrium.at(j).at(k) * 1e-24;
          // Calculate mass density and save for each temperature and pressure
          massdensity_gcm3 = cellmass_grams / volume_cm3;
          AGL_data.mass_density_gcm3.at(j).push_back(massdensity_gcm3);
        }
      }

      // Initializes minimum and maximum values of calculated Debye temperature, heat capacity, Gruneisen parameter, and vibrational free energy
      // Used to determine range for fitting Debye temperature to heat capacity, and to determine x and y ranges for plotting graphs
      tdmin = AGL_data.DebyeTemperature0pressure.at(0);
      tdmax = AGL_data.DebyeTemperature0pressure.at(0);

      cvmin = AGL_data.CvunitkB0pressure.at(0);
      cvmax = AGL_data.CvunitkB0pressure.at(0);

      cpmin = AGL_data.CpunitkB0pressure.at(0);
      cpmax = AGL_data.CpunitkB0pressure.at(0);

      gamin = AGL_data.GruneisenParameter0pressure.at(0);
      gamax = AGL_data.GruneisenParameter0pressure.at(0);

      vfmin = AGL_data.HelmholtzEnergy0pressuremeV.at(0);
      vfmax = AGL_data.HelmholtzEnergy0pressuremeV.at(0);

      // Prefactor for Cv: number of atoms per cell * kB
      // Setting kB = 1 so Cv is in units of kB/cell
      nkb = nkb * AGL_data.natoms;

      // Dulong-Petit value for heat capacity (high-temperature saturation limit)
      // [OBSOLETE] double cvdp = 3.0 * nkb;

      // Half Dulong-Petit value - used for Debye temperatures given in Table 23.3, Solid State Physics, Ashcroft & Mermin
      // [OBSOLETE] double cvdph = cvdp / 2.0;

      int ntdpoints = 0;
      // [OBSOLETE] int jtdcvh = 0;
      int jtdbest = 0;
      int jtdmax = 0;
      int jtdmin = 0;
      int jtd300K = 0;

      // Determines minimum and maximum values of calculated Debye temperature, heat capacity, Gruneisen parameter, and vibrational free energy
      // Used to determine range for fitting Debye temperature to heat capacity, and to determine x and y ranges for plotting graphs
      //for (int j = 1; j < AGL_data.ntemp; j++)
      for (uint j = 1; j < AGL_data.temperature_external.size(); j++) {
        if(AGL_data.DebyeTemperature0pressure.at(j) > tdmax) {
          tdmax = AGL_data.DebyeTemperature0pressure.at(j);
          jtdmax = j;
        }
        if(AGL_data.DebyeTemperature0pressure.at(j) < tdmin) {
          tdmin = AGL_data.DebyeTemperature0pressure.at(j);
          jtdmin = j;
        }
        if(AGL_data.CvunitkB0pressure.at(j) > cvmax) {
          cvmax = AGL_data.CvunitkB0pressure.at(j);
        }
        if(AGL_data.CvunitkB0pressure.at(j) < cvmin) {
          cvmin = AGL_data.CvunitkB0pressure.at(j);
        }
        if(AGL_data.CpunitkB0pressure.at(j) > cpmax) {
          cpmax = AGL_data.CpunitkB0pressure.at(j);
        }
        if(AGL_data.CpunitkB0pressure.at(j) < cpmin) {
          cpmin = AGL_data.CpunitkB0pressure.at(j);
        }
        if(AGL_data.GruneisenParameter0pressure.at(j) > gamax) {
          gamax = AGL_data.GruneisenParameter0pressure.at(j);
        }
        if(AGL_data.GruneisenParameter0pressure.at(j) < gamin) {
          gamin = AGL_data.GruneisenParameter0pressure.at(j);
        }
        if(AGL_data.HelmholtzEnergy0pressuremeV.at(j) > vfmax) {
          vfmax = AGL_data.HelmholtzEnergy0pressuremeV.at(j);
        }
        if(AGL_data.HelmholtzEnergy0pressuremeV.at(j) < vfmin) {
          vfmin = AGL_data.HelmholtzEnergy0pressuremeV.at(j);
        }
        // Checks if heat capacity is NaN at any point
        // If so, gives a warning and skips the remaining points
        // At the end of the loop, ntdpoints will be the number of plottable points
        if(std::isnan(AGL_data.CvunitkB0pressure.at(j))) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Heat capacity at T = " << AGL_data.temperature_external.at(j) << "K is NaN" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          break;
        } else {
          ntdpoints++;
        }	
        // Finds temperature point corresponding to 300K (~room temperature)
        if((AGL_data.temperature_external.at(j) > (300.0 - USER_STEMPERATURE)) && (AGL_data.temperature_external.at(j) < (300.0 + USER_STEMPERATURE))) {
          jtd300K = j;
        }
        // Finds temperature point corresponding to half Dulong-Petit value for heat capacity
        // [OBSOLETE] if(AGL_data.CvunitkB0pressure.at(j) < cvdph) {
        // [OBSOLETE] jtdcvh = j;
        // [OBSOLETE] }	  
      }
      AGL_data.max_temperature = AGL_data.temperature_external.at(ntdpoints);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Maximum fitted temperature = " << AGL_data.max_temperature << "K" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // Checks that the maximum temperature calculated in the AGL method exceeds 300K
      // Gives a warning if the maximum calculated temperature is less than 300K
      // Also sets jtd300K equal to ntdpoints (as this is the closest fitted value to 300K)
      if(AGL_data.max_temperature < 300.0) {
        aurostd::StringstreamClean(aus);	
        aus << _AGLSTR_WARNING_ + "Maximum temperature calculated in AGL has value of " << AGL_data.max_temperature << endl;  
        aus << _AGLSTR_WARNING_ + "Maximum temperature value calculated using AGL is less than 300K" << endl;  
        if(AGL_data.EV_noise) {
          aus << _AGLSTR_ERROR_ + "Noise in E-V data caused run failure: re-run with increased STATIC_KPPRA" << endl;
        } else if(abs(AGL_data.itdiff) > 1)  {
          aus << _AGLSTR_ERROR_ + "Minimum in E-V data not at center of E-V data: check relaxation of initial structure" << endl;	 
        }
        jtd300K = ntdpoints;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
      }

      // Sets x-y range for plotting Debye temperature, heat capacity, Gruneisen parameter, and vibrational free energy
      DEB_min = tdmin - fmod(tdmin,100.0);
      DEB_max = tdmax - fmod(tdmax,100.0) + 100.0;

      CV_min = cvmin - fmod(cvmin,1.0);
      CV_max = cvmax - fmod(cvmax,1.0) + 1.0;

      CP_min = cpmin - fmod(cpmin,1.0);
      CP_max = cpmax - fmod(cpmax,1.0) + 1.0;

      GA_min = gamin - fmod(gamin,1.0);
      GA_max = gamax - fmod(gamax,1.0) + 1.0;

      if(vfmin < 0.0) { 
        VF_min = vfmin - fmod(vfmin,1.0) - 1.0;
      } else {
        VF_min = vfmin - fmod(vfmin,1.0);
      }
      if(vfmax < 0.0) {
        VF_max = vfmax - fmod(vfmax,1.0);
      } else {
        VF_max = vfmax - fmod(vfmax,1.0) + 1.0;
      }

      Tmin = AGL_data.temperature_external.at(0);
      Tmax = AGL_data.temperature_external.at(AGL_data.temperature_external.size()-1);

      // Finds Debye temperature which produces best fit to heat capacity data
      // [OBSOLETE] aglerror = AGL_functions::cvdebfit(AGL_data.CvunitkB0pressure, tdmin, tdmax, nkb, AGL_data.temperature_external, ntdpoints, tdbest, _AGL_data& AGL_data, FileMESSAGE);
      aglerror = AGL_functions::cvdebfit(tdmin, tdmax, nkb, ntdpoints, tdbest, AGL_data, FileMESSAGE);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failed to fit Debye temperature to heat capacity" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aglerror;
      }

      // Finds temperature value which produces Debye temperature which is closest to the value which produces the best fit to the heat capacity data
      // [OBSOLETE] double difftd, difftdmin;
      double difftdmin = aurostd::abs(tdbest - AGL_data.DebyeTemperature0pressure.at(0));
      double difftd = difftdmin;
      for (uint j = 1; j < AGL_data.temperature_external.size(); j++) {
        difftd = aurostd::abs(tdbest - AGL_data.DebyeTemperature0pressure.at(j));
        if(difftd < difftdmin) {
          jtdbest = j;
          difftdmin = difftd;
        }
      }

      // Checks Gruneisen parameter value to ensure it is in the range from 0.0 to 4.0
      // Values outside of this range could indicate a problem in the calculation and result in AGL issueing a warning to the user
      if((AGL_data.GruneisenParameter0pressure.at(jtdbest) < 0.0) || (AGL_data.GruneisenParameter0pressure.at(jtdbest) > 4.0)) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "Gruneisen parameter has value of " << AGL_data.GruneisenParameter0pressure.at(jtdbest) << endl;  
        aus << _AGLSTR_WARNING_ + "Gruneisen parameter is outside the range from 0.0 to 4.0 and may indicate a problem with the calculation" << endl;
        aus << _AGLSTR_WARNING_ + "It is recommended that the user check the E(V) data used for noise" << endl;
        aus << _AGLSTR_WARNING_ + "It may be necessary to re-run the E(V) calculations with a denser k-point mesh or higher basis set energy cut-off" << endl;	
        if(AGL_data.EV_noise) {
          aus << _AGLSTR_ERROR_ + "Noise in E-V data caused run failure: re-run with increased STATIC_KPPRA" << endl;
        } else if(abs(AGL_data.itdiff) > 1)  {
          aus << _AGLSTR_ERROR_ + "Minimum in E-V data not at center of E-V data: check relaxation of initial structure" << endl;	 
        }
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
      }

      // Calculate thermal conductivity at Debye temperature using method of Julian / Slack (Liebfried-Schlomann equation)

      // Determine acoustic Debye temperature
      double acousticthetaD = tdbest / (pow(AGL_data.natoms, third));

      // Determine temperature value which is closest to the acoustic Debye temperature
      double diffatd, diffatdmin;
      uint jtdatd = 0;
      diffatdmin = aurostd::abs(acousticthetaD - AGL_data.temperature_external.at(0));
      for (uint j = 1; j < AGL_data.temperature_external.size(); j++) {
        diffatd = aurostd::abs(acousticthetaD - AGL_data.temperature_external.at(j));
        if(diffatd < diffatdmin) {
          jtdatd = j;
          diffatdmin = diffatd;
        }
      }    

      // Equilibrium volume at acoustic Debye temperature in units of m^3
      // [OBSOLETE] double voleqDm = (AGL_data.VolumeEquilibrium.at(jtdatd).at(0) / pow(angstrom2bohr, 3.0)) * 1e-30;
      // [OBSOLETE] double voleq300K = (AGL_data.VolumeEquilibrium.at(jtd300K).at(0) / pow(angstrom2bohr, 3.0)) * 1e-30;
      double voleqDm = AGL_data.VolumeEquilibrium.at(jtdatd).at(0) * 1e-30;
      double voleq300K = AGL_data.VolumeEquilibrium.at(jtd300K).at(0) * 1e-30;

      // Lattice thermal conductivity as a function of temperature
      double kappaD;
      vector<double> kappaT;

      aglerror = AGL_functions::thermalconductD(acousticthetaD, AGL_data.temperature_external, kappaD, kappaT, AGL_data.GruneisenParameter0pressure.at(jtdbest), voleqDm, avmasskg);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failed to calculate thermal conductivity" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aglerror;
      }

      // Checks that fitted Debye temperature is less than the maximum temperature calculated in the AGL method
      // Gives a warning if the maximum calculated temperature is less than the Debye temperature
      if(AGL_data.max_temperature < AGL_data.DebyeTemperature0pressure.at(jtdbest)) {
        aurostd::StringstreamClean(aus);	
        aus << _AGLSTR_WARNING_ + "Maximum temperature calculated in AGL has value of " << AGL_data.max_temperature << endl;  
        aus << _AGLSTR_WARNING_ + "Debye temperature has value of " << AGL_data.DebyeTemperature0pressure.at(jtdbest) << endl;  
        aus << _AGLSTR_WARNING_ + "Debye temperature exceeds maximum temperature value calculated using AGL" << endl;  
        if(AGL_data.EV_noise) {
          aus << _AGLSTR_ERROR_ + "Noise in E-V data caused run failure: re-run with increased STATIC_KPPRA" << endl;
        } else if(abs(AGL_data.itdiff) > 1)  {
          aus << _AGLSTR_ERROR_ + "Minimum in E-V data not at center of E-V data: check relaxation of initial structure" << endl;	 
        }
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
      }

      // Writes thermal property values at 300K in file for REST-API
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << "[AGL_RESULTS]START" << endl;
      oss << "agl_thermal_conductivity_300K=" << kappaT.at(jtd300K) << "  (W/m*K)" << endl;
      oss << "agl_debye=" << AGL_data.DebyeTemperature0pressure.at(jtdbest) << "  (K)" << endl;
      oss << "agl_acoustic_debye=" << acousticthetaD << "  (K)" << endl;
      oss << "agl_gruneisen=" << AGL_data.GruneisenParameter0pressure.at(jtdbest) << " " << endl;
      oss << "agl_heat_capacity_Cv_300K=" << AGL_data.CvunitkB0pressure.at(jtd300K) << "  (kB/cell)" << endl;
      oss << "agl_heat_capacity_Cp_300K=" << AGL_data.CpunitkB0pressure.at(jtd300K) << "  (kB/cell)" << endl;
      oss << "agl_thermal_expansion_300K=" << AGL_data.ThermalExpansion0pressure.at(jtd300K) << "  (1/K)" << endl;
      oss << "agl_bulk_modulus_static_300K=" << AGL_data.bulkmodulusstatic_0pressure.at(jtd300K) << "  (GPa)" << endl;
      oss << "agl_bulk_modulus_isothermal_300K=" << AGL_data.bulkmodulusisothermal_0pressure.at(jtd300K) << "  (GPa)" << endl;
      oss << "agl_poisson_ratio_source=" << AGL_data.poissonratiosource << endl;
      oss << "agl_vibrational_free_energy_300K_cell=" << AGL_data.HelmholtzEnergy0pressuremeV.at(jtd300K) << "  (meV/cell)" << endl;
      oss << "agl_vibrational_free_energy_300K_atom=" << AGL_data.HelmholtzEnergy0pressuremeV.at(jtd300K) / AGL_data.natoms << "  (meV/atom)" << endl;
      oss << "agl_vibrational_entropy_300K_cell=" << AGL_data.Entropy0pressuremeV.at(jtd300K) << "  (meV/cell*K)" << endl;
      oss << "agl_vibrational_entropy_300K_atom=" << AGL_data.Entropy0pressuremeV.at(jtd300K) / AGL_data.natoms << "  (meV/atom*K)" << endl;
      oss << "[AGL_RESULTS]STOP" << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << "[AGL_THERMAL_PROPERTIES_TEMPERATURE]START" << endl;
      // [OBSOLETE] oss << "#  T (K)" << "\t" << "Kappa (W/(m*K)) " << "\t" << "Debye temp. (K) " << "\t" << "Gruneisen Parameter" << "\t" << "Cv (kB/cell)" << "\t" << "Cp (kB/cell)" << "\t" << "Thermal Expan. (10^5/K)" << "\t" << "            B_static (Gpa)" << "\t" << "    B_isothermal (GPa)" << endl; 
      oss << "#  T (K)" << "        " << aurostd::PaddedPRE("Kappa (W/(m*K))",8," ") << "         " << aurostd::PaddedPRE("Debye temp. (K)",9," ") << "         " << aurostd::PaddedPRE("Gruneisen Parameter",9," ") << "     " << aurostd::PaddedPRE("Cv (kB/cell)",5," ") << "    " << aurostd::PaddedPRE("Cp (kB/cell)",4," ") << "    " << aurostd::PaddedPRE("Thermal Expan. (10^{-5}/K)",4," ") << "             " << aurostd::PaddedPRE("B_static (Gpa)",13," ") << "          " << aurostd::PaddedPRE("B_isothermal (GPa)",10," ") << endl;
      for (int i = 0; i < ntdpoints; i++) {
        // [OBSOLETE] oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "\t" << setw(15) << setprecision(6) << kappaT.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.DebyeTemperature0pressure.at(i) << "\t" << setw(27) << setprecision(6) << AGL_data.GruneisenParameter0pressure.at(i) << "\t" << setw(12) << setprecision(6) << AGL_data.CvunitkB0pressure.at(i) << "\t" << setw(12) << setprecision(6) << AGL_data.CpunitkB0pressure.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.ThermalExpansion0pressure.at(i) << "\t" << setw(12) << setprecision(6) << "\t" << AGL_data.bulkmodulusstatic_0pressure.at(i) << "\t" << setw(22) << setprecision(6) << AGL_data.bulkmodulusisothermal_0pressure.at(i)  << endl;
        oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "        " << setw(15) << setprecision(6) << kappaT.at(i) << " " << setw(23) << setprecision(6) << AGL_data.DebyeTemperature0pressure.at(i) << "   " << setw(25) << setprecision(6) << AGL_data.GruneisenParameter0pressure.at(i) << "      " << setw(12) << setprecision(6) << AGL_data.CvunitkB0pressure.at(i) << "   " << setw(12) << setprecision(6) << AGL_data.CpunitkB0pressure.at(i) << "    " << setw(26) << setprecision(6) << AGL_data.ThermalExpansion0pressure.at(i) * 100000.0 << "               " << setw(12) << setprecision(6) << AGL_data.bulkmodulusstatic_0pressure.at(i) << "      " << setw(22) << setprecision(6) << AGL_data.bulkmodulusisothermal_0pressure.at(i) << endl;
      }
      oss << "[AGL_THERMAL_PROPERTIES_TEMPERATURE]STOP" << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << "[AGL_ENERGIES_TEMPERATURE]START" << endl;
      // [OBSOLETE] oss << "#  T (K)" << "\t" << "G (eV/cell) " << "\t" << "Fvib (meV/cell) " << "\t" << "Uvib (meV/cell) "  << "\t" << "Svib (meV/cell) " << "\t" << "G (eV/atom) " << "\t" << "Fvib (meV/atom) " << "\t" << "Uvib (meV/atom) "  << "\t" << "Svib (meV/atom) " << endl;
      oss << "#  T (K)" << "        " << aurostd::PaddedPRE("G (eV/cell)",8," ") << "     " << aurostd::PaddedPRE("Fvib (meV/cell)",5," ") << "         " << aurostd::PaddedPRE("Uvib (meV/cell)",9," ")  << "         " << aurostd::PaddedPRE("Svib (meV/cell*K)",9," ") << "         " << aurostd::PaddedPRE("G (eV/atom)",8," ") << "     " << aurostd::PaddedPRE("Fvib (meV/atom)",5," ") << "         " << aurostd::PaddedPRE("Uvib (meV/atom)",9," ")  << "         " << aurostd::PaddedPRE("Svib (meV/atom*K)",9," ") << endl;
      for (int i = 0; i < ntdpoints; i++) {
        // [OBSOLETE] oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "\t" << setw(11) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) << "\t" << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) << "\t" << setw(19) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) / AGL_data.natoms << "\t" << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) / AGL_data.natoms << "\t" << setw(23) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) / AGL_data.natoms << "\t" << setw(23) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) / AGL_data.natoms << endl;
        oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "        " << setw(11) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) << "     " << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) << "    " << setw(20) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) << "    " << setw(22) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) << "     " << setw(15) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) / AGL_data.natoms << "     " << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) / AGL_data.natoms << "    " << setw(20) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) / AGL_data.natoms << "    " << setw(22) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) / AGL_data.natoms << endl;
      }
      oss << "[AGL_ENERGIES_TEMPERATURE]STOP" << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
      string ofileafaglname = AGL_data.dirpathname + "/aflow.agl.out";
      if(!aurostd::stringstream2file(oss, ofileafaglname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file aflow.agl.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      oss.clear();
      oss.str(std::string());

      // Calculates thermal conductivity using user supplied values for Debye temperature and Gruneisen parameter if they are supplied in _AFLOWIN_
      // This is useful for testing and debugging where experimental values are available
      // [OBSOLETE] if(USER_THETA_COND.isflag("THETA_EQUAL") && USER_GRUNEISEN.isflag("GRUNEISEN_EQUAL"))
      if(USER_THETA_COND.flag("THETA_EQUAL") && USER_GRUNEISEN.flag("GRUNEISEN_EQUAL")) {  
        double kappaU = 0.0;
        vector<double> kappaTU;
        double gammaU = USER_GRUNEISEN.getattachedutype<double>("GRUNEISEN_EQUAL");
        double thetaU = USER_THETA_COND.getattachedutype<double>("THETA_EQUAL");
        // Check that values are not zero or negative; warn user and reset to default values if they are zero or negative
        if(gammaU < tolzero) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "User Gruneisen parameter = " << gammaU << " <= 0.0" << endl;  
          gammaU = 2.0;
          aus << _AGLSTR_WARNING_ + "Increasing Gruneisen parameter to default value of " << gammaU << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(thetaU < tolzero) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "User Debye temperature = " << thetaU << " <= 0.0" << endl;  
          thetaU = 300.0;
          aus << _AGLSTR_WARNING_ + "Increasing Debye temperature to default value of " << thetaU << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }

        // Finds temperature value which is closest to the user's value of the Debye temperature to get the equilibrium volume
        double difftdU = 0.0, difftdminU = 0.0;
        int jtdbestU = 0;
        difftdminU = aurostd::abs(thetaU - AGL_data.temperature_external.at(0));
        for (uint j = 1; j < AGL_data.temperature_external.size(); j++) {
          difftdU = aurostd::abs(thetaU - AGL_data.temperature_external.at(j));
          if(difftdU < difftdminU) {
            jtdbestU = j;
            difftdminU = difftdU;
          }
        }

        // Equilibrium volume at Debye temperature in units of m^3
        // [OBSOLETE] double voleqDmU = (AGL_data.VolumeEquilibrium.at(jtdbestU).at(0) / pow(angstrom2bohr, 3.0)) * 1e-30;
        double voleqDmU = AGL_data.VolumeEquilibrium.at(jtdbestU).at(0) * 1e-30;

        aglerror = AGL_functions::thermalconductD(thetaU, AGL_data.temperature_external, kappaU, kappaTU, gammaU, voleqDmU, avmasskg);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to calculate thermal conductivity with user-provided Debye temperature and Gruneisen parameter" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // Writes user defined Debye temperature, Gruneisen parameter, and thus-calculated thermal conductivity values in script-readable format
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "material=" << AGL_data.sysname << " | ";
        oss << "nspecies=" << xvasp.str.num_each_type.size() << " | ";
        oss << "natoms=" << AGL_data.natoms << " | ";
        oss << "average_atomic_mass=" << avmasskg << " | ";
        oss << "debye_temperature_user=" << thetaU << " | ";
        oss << "gruneisen_parameter_user=" << gammaU << " | ";
        oss << "thermal_conductivity_user_300K=" << kappaTU.at(jtd300K) << " | ";
        oss << "thermal_conductivity_user_tdebye=" << kappaU << " | ";
        oss << "equilibrium_unit_cell_volume_bestfit=" << voleqDm << " | ";
        oss << "equilibrium_unit_cell_volume_user=" << voleqDmU << " | ";
        string ofileutclibname = AGL_data.dirpathname + "/AGL_thermal_conductivity_user_lib.dat";;
        if(!aurostd::stringstream2file(oss, ofileutclibname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_thermal_conductivity_user_lib.dat" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());
      }

      // Writes thermal properties as a function of temperature to a file in an easy-to-plot format
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss << "# Thermal property values obtained from GIBBS" << endl;
      oss << "# Thermal conductivity at acoustic Debye temperature = " << kappaD << "W/(m*K) " << endl;
      oss << "# Thermal conductivity at " << AGL_data.temperature_external.at(jtd300K) << "K = " << kappaT.at(jtd300K) << "W/(m*K) " << endl;
      oss << "# Debye temperature giving best fit to Cv data = " << tdbest << "K " << endl;
      oss << "# Minimum value of Debye temperature = " << tdmin << "K " << "at T = " << AGL_data.temperature_external.at(jtdmin) << "K " << endl;
      oss << "# Maximum value of Debye temperature = " << tdmax << "K " << "at T = " << AGL_data.temperature_external.at(jtdmax) << "K " << endl;
      oss << "# Debye temperature at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.DebyeTemperature0pressure.at(jtd300K) << "K " << endl;
      oss << "# Acoustic Debye temperature = " << acousticthetaD << "K " << endl;
      oss << "# Gruneisen parameter at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.GruneisenParameter0pressure.at(jtd300K) << endl;
      oss << "# Gruneisen parameter from polynomial at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.gamma_poly.at(jtd300K) << endl;
      oss << "# Equilibrium volume at " <<  AGL_data.temperature_external.at(jtd300K) << "K = " << voleq300K << "m^3" << endl;
      oss << "# Equilibrium volume at " <<  AGL_data.temperature_external.at(jtdatd) << "K = " << voleqDm << "m^3" << endl;
      oss << "# Heat capacity at constant volume (Cv) at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.CvunitkB0pressure.at(jtd300K) << "kB/cell " << endl;
      oss << "# Heat capacity at constant pressure (Cp) at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.CpunitkB0pressure.at(jtd300K) << "kB/cell " << endl;
      oss << "# Thermal expansion at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.ThermalExpansion0pressure.at(jtd300K) * 100000.0 << "10^{-5}/K " << endl;
      oss << "# Static bulk modulus at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.bulkmodulusstatic_0pressure.at(jtd300K) << "GPa " << endl;
      oss << "# Isothermal bulk modulus at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.bulkmodulusisothermal_0pressure.at(jtd300K) << "GPa " << endl;
      // [OBSOLETE] oss << "#  T (K)" << "\t" << "Kappa (W/(m*K)) " << "\t" << "Debye temp. (K) " << "\t" << "Gruneisen Parameter" << "\t" << "Cv (kB/cell)" << "\t" << "Cp (kB/cell)" << "\t" << "Thermal Expan. (10^5/K)" << "\t" << "            B_static (Gpa)" << "\t" << "    B_isothermal (GPa)" << endl; 
      oss << "#  T (K)" << "        " << aurostd::PaddedPRE("Kappa (W/(m*K))",8," ") << "         " << aurostd::PaddedPRE("Debye temp. (K)",9," ") << "         " << aurostd::PaddedPRE("Gruneisen Parameter",9," ") << "     " << aurostd::PaddedPRE("Cv (kB/cell)",5," ") << "    " << aurostd::PaddedPRE("Cp (kB/cell)",4," ") << "    " << aurostd::PaddedPRE("Thermal Expan. (10^{-5}/K)",4," ") << "             " << aurostd::PaddedPRE("B_static (Gpa)",13," ") << "          " << aurostd::PaddedPRE("B_isothermal (GPa)",10," ") << endl;
      for (int i = 0; i < ntdpoints; i++) {
        // [OBSOLETE] oss << setw(8) << setprecision(2) << fixed << aurostd::PaddedPOST(AGL_data.temperature_external.at(i),8," ") << "        " << setw(15) << setprecision(6) << aurostd::PaddedPOST(kappaT.at(i),4," ") << "    " << setw(23) << setprecision(6) << aurostd::PaddedPOST(AGL_data.DebyeTemperature0pressure.at(i),2," ") << "  " << setw(25) << setprecision(6) << aurostd::PaddedPOST(AGL_data.GruneisenParameter0pressure.at(i),4," ") << "    " << setw(12) << setprecision(6) << aurostd::PaddedPOST(AGL_data.CvunitkB0pressure.at(i),3," ") << "   " << setw(12) << setprecision(6) << aurostd::PaddedPOST(AGL_data.CpunitkB0pressure.at(i),3," ") << "   " << setw(23) << setprecision(6) << aurostd::PaddedPOST(AGL_data.ThermalExpansion0pressure.at(i),14," ") << "              " << setw(12) << setprecision(6) << aurostd::PaddedPOST(AGL_data.bulkmodulusstatic_0pressure.at(i),4," ") << "     " << setw(22) << setprecision(6) << aurostd::PaddedPOST(AGL_data.bulkmodulusisothermal_0pressure.at(i),2," ") << "  " << endl;
        oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "        " << setw(15) << setprecision(6) << kappaT.at(i) << " " << setw(23) << setprecision(6) << AGL_data.DebyeTemperature0pressure.at(i) << "   " << setw(25) << setprecision(6) << AGL_data.GruneisenParameter0pressure.at(i) << "      " << setw(12) << setprecision(6) << AGL_data.CvunitkB0pressure.at(i) << "   " << setw(12) << setprecision(6) << AGL_data.CpunitkB0pressure.at(i) << "    " << setw(26) << setprecision(6) << AGL_data.ThermalExpansion0pressure.at(i) * 100000.0 << "               " << setw(12) << setprecision(6) << AGL_data.bulkmodulusstatic_0pressure.at(i) << "      " << setw(22) << setprecision(6) << AGL_data.bulkmodulusisothermal_0pressure.at(i) << endl;
      }
      string ofiletemppropsname = AGL_data.dirpathname + "/AGL_thermal_properties_temperature.out";
      if(!aurostd::stringstream2file(oss, ofiletemppropsname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file AGL_thermal_properties_temperature.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss.clear();
      oss.str(std::string());

      // Writes energy values (Gibbs, Helmholtz, etc.) in units of both eV/cell and eV/atom in a plottable format
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss << "# Energy values obtained from AGL" << endl;
      // [OBSOLETE] oss << "#  T (K)" << "\t" << "G (eV/cell) " << "\t" << "Fvib (meV/cell) " << "\t" << "Uvib (meV/cell) "  << "\t" << "Svib (meV/cell) " << "\t" << "G (kJ/mol)" << "\t" << "Fvib (kJ/mol)" << "\t" << "Uvib (kJ/mol)" << "\t" << "Svib (kJ/mol)" << endl;
      oss << "#  T (K)" << "        " << aurostd::PaddedPRE("G (eV/cell)",8," ") << "     " << aurostd::PaddedPRE("Fvib (meV/cell)",5," ") << "         " << aurostd::PaddedPRE("Uvib (meV/cell)",9," ")  << "         " << aurostd::PaddedPRE("Svib (meV/cell*K)",9," ") << "         " << aurostd::PaddedPRE("G (eV/atom)",8," ") << "     " << aurostd::PaddedPRE("Fvib (meV/atom)",5," ") << "         " << aurostd::PaddedPRE("Uvib (meV/atom)",9," ")  << "         " << aurostd::PaddedPRE("Svib (meV/atom*K)",9," ") << endl;
      for (int i = 0; i < ntdpoints; i++) {
        // [OBSOLETE] oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "\t" << setw(11) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) << "\t" << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) << "\t" << setw(18) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressure.at(i) << "\t" << setw(13) << setprecision(6) << AGL_data.HelmholtzEnergy0pressure.at(i) << "\t" << setw(13) << setprecision(6) << AGL_data.InternalEnergy0pressure.at(i) << "\t" << setw(13) << setprecision(6) << AGL_data.Entropy0pressure.at(i) << endl;
        oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "        " << setw(11) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) << "     " << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) << "    " << setw(20) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) << "    " << setw(22) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) << "     " << setw(15) << setprecision(6) << AGL_data.GibbsFreeEnergy0pressureeV.at(i) / AGL_data.natoms << "     " << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergy0pressuremeV.at(i) / AGL_data.natoms << "    " << setw(20) << setprecision(6) << AGL_data.InternalEnergy0pressuremeV.at(i) / AGL_data.natoms << "    " << setw(22) << setprecision(6) << AGL_data.Entropy0pressuremeV.at(i) / AGL_data.natoms << endl;
      }
      string ofileenergyname = AGL_data.dirpathname + "/AGL_energies_temperature.out";
      if(!aurostd::stringstream2file(oss, ofileenergyname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energies_temperature.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss.clear();
      oss.str(std::string());

      // Write out values for all temperatures and pressures for several different properties
      if(USER_KAPPA_VOLUME.option) {
        vector<double> kappaTV;
        aglerror = AGL_functions::thermalconductDvol(acousticthetaD, AGL_data.temperature_external, kappaTV, AGL_data.GruneisenParameter0pressure.at(jtdbest), AGL_data.VolumeEquilibrium, avmasskg);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to calculate thermal conductivity using temperature-dependent volume" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        // Writes thermal properties as a function of temperature to a file in an easy-to-plot format
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Thermal property values obtained from GIBBS" << endl;
        oss << "# Thermal conductivity at acoustic Debye temperature = " << kappaD << "W/(m*K) " << endl;
        oss << "# Thermal conductivity at " << AGL_data.temperature_external.at(jtd300K) << "K = " << kappaTV.at(jtd300K) << "W/(m*K) " << endl;
        oss << "# Debye temperature giving best fit to Cv data = " << tdbest << "K " << endl;
        oss << "# Minimum value of Debye temperature = " << tdmin << "K " << "at T = " << AGL_data.temperature_external.at(jtdmin) << "K " << endl;
        oss << "# Maximum value of Debye temperature = " << tdmax << "K " << "at T = " << AGL_data.temperature_external.at(jtdmax) << "K " << endl;
        oss << "# Debye temperature at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.DebyeTemperature0pressure.at(jtd300K) << "K " << endl;
        oss << "# Acoustic Debye temperature = " << acousticthetaD << "K " << endl;
        oss << "# Gruneisen parameter at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.GruneisenParameter0pressure.at(jtd300K) << endl;
        oss << "# Gruneisen parameter from polynomial at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.gamma_poly.at(jtd300K) << endl;
        oss << "# Heat capacity at constant volume (Cv) at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.CvunitkB0pressure.at(jtd300K) << "kB/cell " << endl;
        oss << "# Heat capacity at constant pressure (Cp) at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.Cpkjmol0pressure.at(jtd300K) << "J/(mol*K) " << endl;
        oss << "# Thermal expansion at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.ThermalExpansion0pressure.at(jtd300K)  * 100000.0 << "10^{-5}/K " << endl;
        oss << "# Static bulk modulus at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.bulkmodulusstatic_0pressure.at(jtd300K) << "GPa " << endl;
        oss << "# Isothermal bulk modulus at " << AGL_data.temperature_external.at(jtd300K) << "K = " << AGL_data.bulkmodulusisothermal_0pressure.at(jtd300K) << "GPa " << endl;
        oss << "#  T (K)" << "\t" << "Kappa (W/(m*K)) " << "\t" << "Debye temp. (K) " << "\t" << "Gruneisen Parameter" << "\t" << "Cv (kB/cell)" << "\t" << "Cp (kB/cell)" << "\t" << "Thermal Expan. (10^{-5}/K)" << "\t" << "            B_static (Gpa)" << "\t" << "    B_isothermal (GPa)" << endl; 
        for (int i = 0; i < ntdpoints; i++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(i) << "\t" << setw(15) << setprecision(6) << kappaTV.at(i) << "\t" << setw(23) << setprecision(6) << AGL_data.DebyeTemperature0pressure.at(i) << "\t" << setw(27) << setprecision(6) << AGL_data.GruneisenParameter0pressure.at(i) << "\t" << setw(12) << setprecision(6) << AGL_data.CvunitkB0pressure.at(i) << "\t" << setw(14) << setprecision(6) << AGL_data.CpunitkB0pressure.at(i) << "\t" << setw(26) << setprecision(6) << AGL_data.ThermalExpansion0pressure.at(i) * 100000.0 << "\t" << setw(12) << setprecision(6) << "\t" << AGL_data.bulkmodulusstatic_0pressure.at(i) << "\t" << setw(22) << setprecision(6) << AGL_data.bulkmodulusisothermal_0pressure.at(i)  << endl;
        }
        string ofiletempvpropsname = AGL_data.dirpathname + "/AGL_thermal_properties_temperature_vol.out";
        if(!aurostd::stringstream2file(oss, ofiletempvpropsname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_thermal_properties_temperature_vol.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());
      }

      // Write calculated E(V) data including equilibrium volume determined by fitting calculated E(V) data
      // [OBSOLETE] aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "Opening file AGL_energy_volume.out" <<  endl;
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE] oss << "# E(V) data including equilibrium volume determined by fitting calculated E(V) data" << endl;
      // [OBSOLETE] oss << "# Equilibrium volume at zero temperature and pressure = " << AGL_data.volume_equilibrium_0p0T << " Angstrom^3/cell" << endl;
      // [OBSOLETE] oss << "# Energy at equilibrium volume at zero temperature and pressure = " << AGL_data.energy_equilibrium_0p0T << " eV/cell" << endl;
      // [OBSOLETE] oss << "# Equilibrium volume at zero temperature and pressure = " << AGL_data.volume_equilibrium_0p0T / AGL_data.natoms << " Angstrom^3/atom" << endl;
      // [OBSOLETE] oss << "# Energy at equilibrium volume at zero temperature and pressure = " << AGL_data.energy_equilibrium_0p0T / AGL_data.natoms << " eV/atom" << endl;
      // [OBSOLETE] oss << "# Volume (Ang^3/atom)        Energy(eV/atom)        Volume (Ang^3/cell)        Energy(eV/cell)" << endl;
      // [OBSOLETE] vector<double> volumetotal;
      // [OBSOLETE] vector<double> energytotal;
      // [OBSOLETE] vector<double> volume_atom;
      // [OBSOLETE] vector<double> energy_atom;
      // [OBSOLETE] volumetotal = AGL_data.volumeinput;
      // [OBSOLETE] energytotal = AGL_data.energyinput;
      // [OBSOLETE] volumetotal.push_back(AGL_data.volume_equilibrium_0p0T);
      // [OBSOLETE] energytotal.push_back(AGL_data.energy_equilibrium_0p0T);
      // [OBSOLETE] aglerror = AGL_functions::qcksortev(volumetotal, energytotal, FileMESSAGE);
      // [OBSOLETE] for (uint i = 0; i < volumetotal.size(); i++) {
      // [OBSOLETE] volume_atom.push_back(volumetotal.at(i) / AGL_data.natoms);
      // [OBSOLETE] energy_atom.push_back(energytotal.at(i) / AGL_data.natoms);
      // [OBSOLETE] }
      // [OBSOLETE] for (uint i = 0; i < volumetotal.size(); i++) {
      // [OBSOLETE] oss << setw(21) << setprecision(6) << volume_atom.at(i) << "\t" << setw(20) << setprecision(6) << energy_atom.at(i) << "\t" << setw(23) << setprecision(6) << volumetotal.at(i) << "\t" << setw(22) << setprecision(6) << energytotal.at(i) << endl;
      // [OBSOLETE] }
      // [OBSOLETE] string ofileenergyvolume = AGL_data.dirpathname + "/AGL_energy_volume.out";
      // [OBSOLETE] if(!aurostd::stringstream2file(oss, ofileenergyvolume, "WRITE")) {
      // [OBSOLETE] aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energy_volume.out" <<  endl;
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE] return 1;
      // [OBSOLETE] }	
      // [OBSOLETE] oss.clear();
      // [OBSOLETE] oss.str(std::string());

      // Determine minimum and maximum energy and volume values to set plotting boundaries
      ecmin = energytotal.at(0);
      ecmax = energytotal.at(0);
      vcmin = volumetotal.at(0);
      vcmax = volumetotal.at(0);
      eamin = energy_atom.at(0);
      eamax = energy_atom.at(0);
      vamin = volume_atom.at(0);
      vamax = volume_atom.at(0);

      int nevpoints = 0;
      for (uint i = 0; i < volumetotal.size(); i++) {
        if(energytotal.at(i) < ecmin) {
          ecmin = energytotal.at(i);
        }
        if(energytotal.at(i) > ecmax) {
          ecmax = energytotal.at(i);
        }
        if(volumetotal.at(i) < vcmin) {
          vcmin = volumetotal.at(i);
        }
        if(volumetotal.at(i) > vcmax) {
          vcmax = volumetotal.at(i);
        }
        if(energy_atom.at(i) < eamin) {
          eamin = energy_atom.at(i);
        }
        if(energy_atom.at(i) > eamax) {
          eamax = energy_atom.at(i);
        }
        if(volume_atom.at(i) < vamin) {
          vamin = volume_atom.at(i);
        }
        if(volume_atom.at(i) > vamax) {
          vamax = volume_atom.at(i);
        }
        nevpoints++;
      }

      // Set plotting boundaries to nearest integer outside of value ranges
      if(ecmin < 0.0) {
        Ecellmin = ecmin - fmod(ecmin, 1.0) - 1.0;
      } else {
        Ecellmin = ecmin - fmod(ecmin, 1.0);
      }
      if(ecmax < 0.0) { 
        Ecellmax = ecmax - fmod(ecmax, 1.0);
      } else {
        Ecellmax = ecmax - fmod(ecmax, 1.0) + 1.0;
      }

      Vcellmin = vcmin - fmod(vcmin, 1.0);
      Vcellmax = vcmax - fmod(vcmax, 1.0) + 1.0;

      if(eamin < 0.0) {
        Eatommin = eamin - fmod(eamin, 1.0) - 1.0;
      } else {
        Eatommin = eamin - fmod(eamin, 1.0);
      }
      if(eamax < 0.0) {
        Eatommax = eamax - fmod(eamax, 1.0);
      } else {
        Eatommax = eamax - fmod(eamax, 1.0) + 1.0;
      }

      Vatommin = vamin - fmod(vamin, 1.0);
      Vatommax = vamax - fmod(vamax, 1.0) + 1.0;      

      // Write out values for all temperatures and pressures for several different properties
      if(USER_WRITE_FULL_RESULTS.option) {
        // Writes thermal property values at 300K in script-readable format
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());
        oss << "material=" << AGL_data.sysname << " | ";
        oss << "nspecies=" << xvasp.str.num_each_type.size() << " | ";
        oss << "natoms=" << AGL_data.natoms << " | ";
        oss << "average_atomic_mass=" << avmasskg << " | ";
        oss << "unit_cell_volume_m3_acoustic_debye_temperature=" << voleqDm << " | ";
        oss << "unit_cell_volume_m3_300K=" << voleq300K << " | ";
        oss << "thermal_conductivity_300K=" << kappaT.at(jtd300K) << " | ";
        oss << "thermal_conductivity_acoustic_debye_temperature=" << kappaD << " | ";
        oss << "debye_temperature=" << AGL_data.DebyeTemperature0pressure.at(jtdbest) << " | ";
        oss << "acoustic_debye_temperature=" << acousticthetaD << " | ";
        oss << "gruneisen_parameter=" << AGL_data.GruneisenParameter0pressure.at(jtdbest) << " | ";
        oss << "gruneisen_parameter_polynomial=" << AGL_data.gamma_poly.at(jtdbest) << " | ";
        // [OBSOLETE] oss << "heat_capacity_Cv_300K=" << AGL_data.Cv0pressure.at(jtd300K) << " | ";
        // [OBSOLETE] oss << "heat_capacity_Cp_300K=" << AGL_data.Cp0pressure.at(jtd300K) << " | ";
        oss << "heat_capacity_Cv_300K=" << AGL_data.CvunitkB0pressure.at(jtd300K) << " | ";
        oss << "heat_capacity_Cp_300K=" << AGL_data.CpunitkB0pressure.at(jtd300K) << " | ";
        oss << "thermal_expansion_300K=" << AGL_data.ThermalExpansion0pressure.at(jtd300K) * 100000.0 << " | ";
        oss << "bulk_modulus_static_300K=" << AGL_data.bulkmodulusstatic_0pressure.at(jtd300K) << " | ";
        oss << "bulk_modulus_isothermal_300K=" << AGL_data.bulkmodulusisothermal_0pressure.at(jtd300K) << " | ";
        string ofiletplibname = AGL_data.dirpathname + "/AGL_thermal_properties_lib.dat";
        if(!aurostd::stringstream2file(oss, ofiletplibname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_thermal_properties_lib.dat" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        //aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes thermal property values at 300K in file in parseable format
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "agl_thermal_conductivity_300K=" << kappaT.at(jtd300K) << " | ";
        oss << "agl_debye=" << AGL_data.DebyeTemperature0pressure.at(jtdbest) << " | ";
        oss << "agl_acoustic_debye=" << acousticthetaD << " | ";
        oss << "agl_gruneisen=" << AGL_data.GruneisenParameter0pressure.at(jtdbest) << " | ";
        // [OBSOLETE] oss << "agl_heat_capacity_Cv_300K=" << AGL_data.Cv0pressure.at(jtd300K) << " | ";
        // [OBSOLETE] oss << "agl_heat_capacity_Cp_300K=" << AGL_data.Cp0pressure.at(jtd300K) << " | ";
        oss << "agl_heat_capacity_Cv_300K=" << AGL_data.CvunitkB0pressure.at(jtd300K) << " | ";
        oss << "agl_heat_capacity_Cp_300K=" << AGL_data.CpunitkB0pressure.at(jtd300K) << " | ";
        oss << "agl_thermal_expansion_300K=" << AGL_data.ThermalExpansion0pressure.at(jtd300K)  * 100000.0 << " | ";
        oss << "agl_bulk_modulus_static_300K=" << AGL_data.bulkmodulusstatic_0pressure.at(jtd300K) << " | ";
        // [OBSOLETE] oss << "agl_bulk_modulus_isothermal_300K=" << AGL_data.bulkmodulusisothermal_0pressure.at(jtd300K) << " | ";
        oss << "agl_bulk_modulus_isothermal_300K=" << AGL_data.bulkmodulusisothermal_0pressure.at(jtd300K);
        string ofileagllibname = AGL_data.dirpathname + "/agllib.out";
        if(!aurostd::stringstream2file(oss, ofileagllibname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file agllib.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());

        // Writes equilibrium volume in units of Angstrom^3 for all temperatures and pressures in a plottable format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_Volume_equilibrium_p.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Equilibrium volumes for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "V (Ang^3) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.VolumeEquilibrium.at(j).at(k));
            oss << setw(16) << setprecision(6) << (AGL_data.VolumeEquilibrium.at(j).at(k));
          }
          oss << endl;
        }
        string ofilevolminpname = AGL_data.dirpathname + "/AGL_Volume_equilibrium_pressure.out";
        if(!aurostd::stringstream2file(oss, ofilevolminpname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Volume_equilibrium_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes equilibrium volume in terms of x = (V/V0)^(1/3) for all temperatures and pressures in a plottable format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_X_equilibrium_p.dat" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Equilibrium lattice scaling factors for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "x = (V/V0)^(1/3) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.xminsav.at(j).at(k));
            oss << setw(16) << setprecision(6) << (AGL_data.xminsav.at(j).at(k));
          }
          oss << endl;
        }
        string ofilexminpname = AGL_data.dirpathname + "/AGL_X_equilibrium_pressure.out";
        if(!aurostd::stringstream2file(oss, ofilexminpname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_X_equilibrium_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes Gibbs free energy in units of eV/atom for all temperatures and pressures in a plottable format
        // This data can be used to generate (p, T) phase diagrams
        // For a given value of temperature and pressure, the phase with the lowest Gibbs free energy will be the equilibrium phase
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_Gibbs_free_energy_p.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Gibbs free energy values for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "G (eV/atom) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k) / AGL_data.natoms);
            oss << setw(16) << setprecision(6) << (AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k) / AGL_data.natoms);
          }
          oss << endl;
        }
        string ofilegfrenergpname = AGL_data.dirpathname + "/AGL_Gibbs_free_energy_pressure.out";
        if(!aurostd::stringstream2file(oss, ofilegfrenergpname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Gibbs_free_energy_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Write POSCAR with equilibrium volume determined by fitting calculated E(V) data
        xstructure equilibriumstructure = xvasp.str;
        equilibriumstructure.InflateLattice(AGL_data.xlattice_equilibrium_0p0T);
        equilibriumstructure.is_vasp4_poscar_format=TRUE;
        equilibriumstructure.is_vasp5_poscar_format=FALSE;
        oss << equilibriumstructure << endl;
        string ofileposcarequilibrium = AGL_data.dirpathname + "/POSCAR.aglmin";
        if(!aurostd::stringstream2file(oss, ofileposcarequilibrium, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file POSCAR.aglmin" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());      
      }

      // Write out values for all temperatures and pressures for a wide range of thermal properties
      // This option uses a lot of memory but is useful for creating (p, T) phase diagrams
      if(USER_WRITE_ALL_PRESSURES.option && USER_SAVE_ALL_PRESSURES.option) {     
        // Write Helmholtz free energy for all temperatures and pressures
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Writing free energy for all temperatures and pressures" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // Writes vibrational Helmholtz free energy in units of meV/cell for all temperatures and pressures in a plottable format
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Helmholtz free energy values for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "Fvib (meV/cell) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << AGL_data.HelmholtzEnergyPressuremeV.at(j).at(k);
            oss << setw(16) << setprecision(6) << AGL_data.HelmholtzEnergyPressuremeV.at(j).at(k);
          }
          oss << endl;
        }
        string ofilefrenergpname = AGL_data.dirpathname + "/AGL_Helmholtz_free_energy_pressure.out";
        if(!aurostd::stringstream2file(oss, ofilefrenergpname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Helmholtz_free_energy_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes Entropy in units of eV/atom for all temperatures and pressures in a plottable format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file Entropy_pressure.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Vibrational entropy values for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "S (eV/atom*K) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        //for (int k = 0; k < AGL_data.npres; k++)
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          //for (int k = 0; k < AGL_data.npres; k++)
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.EntropyPressuremeV.at(j).at(k) / AGL_data.natoms);
            oss << setw(16) << setprecision(6) << (AGL_data.EntropyPressuremeV.at(j).at(k) / AGL_data.natoms);
          }
          oss << endl;
        }
        string ofileentropypname = AGL_data.dirpathname + "/AGL_Entropy_pressure.out";
        if(!aurostd::stringstream2file(oss, ofileentropypname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Entropy_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes internal energy in units of eV/atom for all temperatures and pressures in a plottable format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file Internal_energy_pressure.dat" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Vibrational internal energy values for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "S (eV/atom) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.InternalEnergyPressuremeV.at(j).at(k) / AGL_data.natoms);
            oss << setw(16) << setprecision(6) << (AGL_data.InternalEnergyPressuremeV.at(j).at(k) / AGL_data.natoms);
          }
          oss << endl;
        }
        string ofileintenergpname = AGL_data.dirpathname + "/AGL_Internal_energy_pressure.out";
        if(!aurostd::stringstream2file(oss, ofileintenergpname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Internal_energy_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes Debye temperature in units of K for all temperatures and pressures in a plottable format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file Debye_temperature_pressure.dat" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Debye temperature values for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "Theta (K) for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.DebyeTemperaturePressure.at(j).at(k));
            oss << setw(16) << setprecision(6) << (AGL_data.DebyeTemperaturePressure.at(j).at(k));
          }
          oss << endl;
        }
        string ofiletdebyepname = AGL_data.dirpathname + "/AGL_Debye_temperature_pressure.out";
        if(!aurostd::stringstream2file(oss, ofiletdebyepname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Debye_temperature_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Writes Gruneisen parameter for all temperatures and pressures in a plottable format
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file Gruneisen_parameter_pressure.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Debye temperature values for all temperatures and pressures obtained from GIBBS" << endl;
        oss << "#  T (K)" << "\t" << "Gamma for Pressure (GPa)" << endl;
        oss << "#  T (K)";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "\t" << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k);
        }
        oss << endl;
        for (int j = 0; j < ntdpoints; j++) {
          oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j);
          //for (int k = 0; k < AGL_data.npres; k++)
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] oss << "\t" << setw(15) << setprecision(6) << (AGL_data.GruneisenParameterPressure.at(j).at(k));
            oss << setw(16) << setprecision(6) << (AGL_data.GruneisenParameterPressure.at(j).at(k));
          }
          oss << endl;
        }
        string ofilegammapname = AGL_data.dirpathname + "/AGL_Gruneisen_parameter_pressure.out";
        if(!aurostd::stringstream2file(oss, ofilegammapname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Gruneisen_parameter_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());
      }

      // Write out data for Hugoniot plots: moved to after plot function calls to allow re-run of GIBBS fitting without truncating pressure and temperature
      // Mass density (g/cm^3); Temperature (K); E_DFT + U_int_vib (kJ/g); Pressure (GPa)
      // Write out data for all temperatures for each pressure
      // [OBSOLETE] if(USER_WRITE_HUGONIOT_INPUT.option) {
      // [OBSOLETE] for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
      // [OBSOLETE]   oss << "# Data for pressure = " << AGL_data.pressure_external.at(k) << endl;
      // [OBSOLETE]   oss << "# rho(g/cc)     T(K)            E(kJ/g)         P(GPa)" << endl;
      // [OBSOLETE]   for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
      // [OBSOLETE] Convert energy from eV/cell to kJ/g
      // [OBSOLETE] energy_kJg = (AGL_data.EnergyDFT_UIntVib.at(j).at(k) * 1.6e-22) / cellmass_grams;
      // [OBSOLETE]  Write data for this (p, T) point
      // [OBSOLETE] oss << AGL_data.mass_density_gcm3.at(j).at(k) << "\t" << AGL_data.temperature_external.at(j) << "\t" << energy_kJg << "\t" << AGL_data.pressure_external.at(k) << endl;
      // [OBSOLETE]  }
      // [OBSOLETE] }
      // [OBSOLETE] string ofilehugoniotinput = AGL_data.dirpathname + "/AGL_Hugoniot_input_data.out";
      // [OBSOLETE] if(!aurostd::stringstream2file(oss, ofilehugoniotinput, "WRITE")) {
      // [OBSOLETE]   aurostd::StringstreamClean(aus);
      // [OBSOLETE]   aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Hugoniot_input_data.out" <<  endl;
      // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE]   return 1;
      // [OBSOLETE] }	
      // [OBSOLETE] oss.clear();
      // [OBSOLETE] oss.str(std::string());
      // [OBSOLETE] }

      // Run Hugoniot calculation: moved to after plot function calls to allow re-run of GIBBS fitting without truncating pressure and temperature
      // Mass density (g/cm^3); Temperature (K); E_DFT + U_int_vib (kJ/g); Pressure (GPa)
      // Pass data for all temperatures for each pressure
      // [OBSOLETE] vector<vector<double> > energy_pT_kJg;
      // [OBSOLETE] vector<vector<double> > hugoniot_output;
      // [OBSOLETE] energy_pT_kJg.resize(AGL_data.pressure_external.size());
      // [OBSOLETE] hugoniot_output.resize(AGL_data.pressure_external.size());
      // [OBSOLETE] if(USER_HUGONIOT_CALC.option) {
      // [OBSOLETE] for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
      // [OBSOLETE] for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
      // Convert energy from eV/cell to kJ/g
      // [OBSOLETE] energy_kJg = (AGL_data.EnergyDFT_UIntVib.at(j).at(k) * 1.6e-22) / cellmass_grams;
      // [OBSOLETE] energy_pT_kJg.at(k).push_back(energy_kJg);
      // [OBSOLETE] }
      // [OBSOLETE] }
      // [OBSOLETE] aglerror = AGL_functions::runHugoniot(AGL_data.pressure_external, AGL_data.temperature_external, AGL_data.mass_density_gcm3, energy_pT_kJg, hugoniot_output, FileMESSAGE);
      // [OBSOLETE] if(aglerror != 0 && aglerror != 2) {
      // [OBSOLETE] aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_ERROR_ + "Failure in Hugoniot calculation" <<  endl;
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE] return aglerror;
      // [OBSOLETE] } else {
      // [OBSOLETE] oss << std::setw(15) << "rho(g/cc)" << std::setw(15) << "T(K)" << std::setw(15) << "E(kJ/g)" << std::setw(15) << "P(GPa)" << endl;
      // [OBSOLETE] for (uint i = 0; i < hugoniot_output.size(); i++) {
      // [OBSOLETE] for (uint j = 0; j < hugoniot_output.at(i).size(); j++) {
      // [OBSOLETE] oss << setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
      // [OBSOLETE] }
      // [OBSOLETE] oss << endl;
      // [OBSOLETE] }
      // [OBSOLETE] string ofilehugoniotoutput = AGL_data.dirpathname + "/AGL_Hugoniot.out";
      // [OBSOLETE] if(!aurostd::stringstream2file(oss, ofilehugoniotoutput, "WRITE")) {
      // [OBSOLETE] aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Hugoniot.out" <<  endl;
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE] return 1;
      // [OBSOLETE] }	
      // [OBSOLETE] }
      // [OBSOLETE] oss.clear();
      // [OBSOLETE] oss.str(std::string());
      // [OBSOLETE] }

      // Plot Debye temperature, heat capacity, Gruneisen parameter, and vibrational free energy as a function of temperature; and energy as a function of cell volume
      if(USER_PLOT_RESULTS.option) {

        // ****************** Plot Debye temperature as a function of temperature **********************
        string debyedatalabel = "debye";
        string xlabel = "Temperature (K)";
        string ylabel = "Debye Temperature (K)";
        string plotname = "Debye Temperature (K)";
        string plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_debye_temperature";
        // Call plotting function
        aglerror = AGL_functions::plotaglresults (AGL_data.temperature_external, AGL_data.DebyeTemperature0pressure, Tmin, Tmax, DEB_min, DEB_max, ntdpoints, debyedatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot Debye temperature" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // ****************** Plot heat capacity (constant volume) as a function of temperature **********************
        string cvdatalabel = "cv";
        xlabel = "Temperature (K)";
        ylabel = "Heat Capacity (kB/cell)";
        plotname = "Heat Capacity (kB/cell)";
        plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_heat_capacity_Cv";
        // Call plotting function	
        aglerror = AGL_functions::plotaglresults (AGL_data.temperature_external, AGL_data.CvunitkB0pressure, Tmin, Tmax, CV_min, CV_max, ntdpoints, cvdatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot heat capacity at constant volume" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // ****************** Plot heat capacity (constant pressure) as a function of temperature **********************
        string cpdatalabel = "cp";
        xlabel = "Temperature (K)";
        ylabel = "Heat Capacity (kB/cell)";
        plotname = "Heat Capacity at constant pressure (kB/cell)";
        plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_heat_capacity_Cp";
        // Call plotting function	
        aglerror = AGL_functions::plotaglresults (AGL_data.temperature_external, AGL_data.CpunitkB0pressure, Tmin, Tmax, CP_min, CP_max, ntdpoints, cpdatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot heat capacity at constant pressure" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // ****************** Plot Gruneisen parameter as a function of temperature **********************
        string gadatalabel = "ga";
        xlabel = "Temperature (K)";
        ylabel = "Gruneisen Parameter";
        plotname = "Gruneisen Parameter";
        plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_gruneisen_parameter";
        // Call plotting function	
        aglerror = AGL_functions::plotaglresults (AGL_data.temperature_external, AGL_data.GruneisenParameter0pressure, Tmin, Tmax, GA_min, GA_max, ntdpoints, gadatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot Gruneisen parameter" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // ****************** Plot vibrational free energy as a function of temperature **********************
        string hedatalabel = "he";
        xlabel = "Temperature (K)";
        ylabel = "Vibrational Free Energy (meV/cell)";
        plotname = "Vibrational Free Energy (meV/cell)";
        plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_vibrational_free_energy";
        // Call plotting function	
        aglerror = AGL_functions::plotaglresults (AGL_data.temperature_external, AGL_data.HelmholtzEnergy0pressuremeV, Tmin, Tmax, VF_min, VF_max, ntdpoints, hedatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot vibrational free energy" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // ****************** Plot calculated energy/cell as a function of volume/cell ***********************
        string evcdatalabel = "evc";
        xlabel = "Volume (Ang^3/cell)";
        ylabel = "Energy (eV/cell)";
        plotname = "Energy as a function of volume (eV/cell)";
        plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_energy_volume_cell";
        // Call plotting function	
        aglerror = AGL_functions::plotaglresults (volumetotal, energytotal, Vcellmin, Vcellmax, Ecellmin, Ecellmax, nevpoints, evcdatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot energy(volume) per cell" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

        // ****************** Plot calculated energy/atom as a function of volume/atom ***********************
        string evadatalabel = "eva";
        xlabel = "Volume (Ang^3/atom)";
        ylabel = "Energy (eV/atom)";
        plotname = "Energy as a function of volume (eV/atom)";
        plotfilename = AGL_data.dirpathname + "/" + AGL_data.sysname + "_energy_volume_atom";
        // Call plotting function	
        aglerror = AGL_functions::plotaglresults (volume_atom, energy_atom, Vatommin, Vatommax, Eatommin, Eatommax, nevpoints, evadatalabel, xlabel, ylabel, plotname, plotfilename, AGL_data.sysname, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to plot energy(volume) per atom" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }

      }

      // Write out data for Hugoniot plots
      // Mass density (g/cm^3); Temperature (K); E_DFT + U_int_vib (kJ/g); Pressure (GPa)
      // Write out data for all temperatures for each pressure
      if(USER_WRITE_HUGONIOT_INPUT.option && !AGL_data.run_all_pressure_temperature) {
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          oss << "# Data for pressure = " << AGL_data.pressure_external.at(k) << endl;
          oss << "# rho(g/cc)     T(K)            E(kJ/g)         P(GPa)" << endl;
          for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
            // Convert energy from eV/cell to kJ/g
            energy_kJg = (AGL_data.EnergyDFT_UIntVib.at(j).at(k) * 1.6e-22) / cellmass_grams;
            // Write data for this (p, T) point
            oss << AGL_data.mass_density_gcm3.at(j).at(k) << "\t" << AGL_data.temperature_external.at(j) << "\t" << energy_kJg << "\t" << AGL_data.pressure_external.at(k) << endl;
          }
        }
        string ofilehugoniotinput = AGL_data.dirpathname + "/AGL_Hugoniot_input_data.out";
        if(!aurostd::stringstream2file(oss, ofilehugoniotinput, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Hugoniot_input_data.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());
      }

      // Run Hugoniot calculation
      // Mass density (g/cm^3); Temperature (K); E_DFT + U_int_vib (kJ/g); Pressure (GPa)
      // Pass data for all temperatures for each pressure
      vector<vector<double> > energy_pT_kJg;
      // [OBSOLETE] vector<vector<double> > hugoniot_output;
      energy_pT_kJg.resize(AGL_data.pressure_external.size());
      hugoniot_output.resize(AGL_data.pressure_external.size());
      bool json_written = false;
      if(USER_HUGONIOT_CALC.option && !AGL_data.run_all_pressure_temperature) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Calculating Hugoniot" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if ((AGL_data.pressure_external.size() < USER_NPRESSURE) || (AGL_data.temperature_external.size() < USER_NTEMPERATURE)) {
          AGL_data.run_all_pressure_temperature = true;
          AGL_data.tdebye.clear();
          AGL_data.pressure_external.resize(USER_NPRESSURE);
          for(uint i = 0; i < AGL_data.pressure_external.size(); i++) {
            AGL_data.pressure_external.at(i) = pres_orig.at(i); 
          }
          // Regenerate temperature values using stepsize and number of values
          // Resize vector in advance to allocate memory to store temperature values
          AGL_data.temperature_external.resize(USER_NTEMPERATURE);
          tempi = 0.0;
          for(uint i = 0; i < USER_NTEMPERATURE; i++) {
            AGL_data.temperature_external.at(i) = tempi; 
            tempi = tempi + USER_STEMPERATURE;
          }

          // Reset E(V) data to initial values
          AGL_data.energyinput.resize(AGL_data.volumeinput_orig.size());
          AGL_data.volumeinput.resize(AGL_data.volumeinput_orig.size());
          for (uint i = 0; i < AGL_data.volumeinput_orig.size(); i++) {
            AGL_data.volumeinput.at(i) = AGL_data.volumeinput_orig.at(i);
            AGL_data.energyinput.at(i) = AGL_data.energyinput_orig.at(i);
            // Print out total energy and volume
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ << "Energy (eV) = " << AGL_data.energyinput.at(i) << endl;
            aus << _AGLSTR_MESSAGE_ << "Volume (Ang^3) = " << AGL_data.volumeinput.at(i) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }

          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Running GIBBS" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

          // Runs GIBBS algorithm within AFLOW to calculate thermal properties
          aglerror = AGL_functions::gibbsrun(AGL_data, FileMESSAGE);
          if (aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Failure in re-run of GIBBS to get all pressure-temperature data to calculate Hugoniot" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 9;
          }

          // Run Hugoniot
          aglerror = AGL_functions::runHugoniotAllTemperaturesAndPressures(AGL_data.AGL_pressure_temperature_energy_list, cellmass_grams, hugoniot_output, AGL_data.hugoniotextrapolate, FileMESSAGE);
          if(aglerror != 0 && aglerror != 2) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Failure in Hugoniot calculation" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 9;
          } else {
            oss << std::setw(15) << "#     rho(g/cc)" << std::setw(15) << "T(K)" << std::setw(15) << "E(kJ/g)" << std::setw(15) << "P(GPa)" << endl;
            for (uint i = 0; i < hugoniot_output.size(); i++) {
              for (uint j = 0; j < hugoniot_output.at(i).size(); j++) {
                // [OBSOLETE] oss << "\t" << std::left << std::setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
                oss << setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
              }
              oss << endl;
            }
            string ofilehugoniotoutput = AGL_data.dirpathname + "/AGL_Hugoniot.out";
            if(!aurostd::stringstream2file(oss, ofilehugoniotoutput, "WRITE")) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Hugoniot.out" <<  endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }
            // Reset aglerror to 0 on successful output of Hugoniot data
            aglerror = 0;
            // [OBSOLETE] aurostd::StringstreamClean(oss);
          }
          oss.clear();
          oss.str(std::string());
          // Write calculated energies to JSON format:
          // E(V) data including equilibrium volume determined by fitting calculated E(V) data
          // Gibbs free energy as a function of temperature and pressure
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "Opening file AGL_energy.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          oss << "{\"energy_volume_atom\": [";
          for (uint i = 0; i < volume_atom.size(); i++) {
            if (i==0) {
              oss << "{\"energy\":" << energy_atom.at(i) << ",\"volume\":" << volume_atom.at(i) << "}";
            } else {
              oss << ",{\"energy\":" << energy_atom.at(i) << ",\"volume\":" << volume_atom.at(i) << "}";
            }
          }
          oss << "],\"enthalpy_pressure_atom\": [";
          for (uint k = 0; k < AGL_data.AGL_pressure_enthalpy_list.size(); k++) {
            if (!std::isnan(AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy)) {
              firstnonnan = k;
              break;
            }
          }
          oss << "{\"pressure\":" << AGL_data.AGL_pressure_enthalpy_list.at(firstnonnan).pressure_external << ",\"enthalpy\":" << AGL_data.AGL_pressure_enthalpy_list.at(firstnonnan).enthalpy / AGL_data.natoms << "}";
          for (uint k = firstnonnan+1; k < AGL_data.AGL_pressure_enthalpy_list.size(); k++) {
            if (!std::isnan(AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy)) {
              oss << ",{\"pressure\":" << AGL_data.AGL_pressure_enthalpy_list.at(k).pressure_external << ",\"enthalpy\":" << AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy / AGL_data.natoms << "}";
            }
          }      	
          oss << "],\"gibbs_free_energy_atom\": [";
          for (uint i = 0; i < AGL_data.AGL_pressure_temperature_energy_list.size(); i++) {
            if (!std::isnan(AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy)) {
              firstnonnan = i;
              break;
            }
          }
          oss << "{\"temperature\":" << AGL_data.AGL_pressure_temperature_energy_list.at(firstnonnan).temperature_external << ",\"pressure\":" << AGL_data.AGL_pressure_temperature_energy_list.at(firstnonnan).pressure_external << ",\"gibbs_energy\":" << AGL_data.AGL_pressure_temperature_energy_list.at(firstnonnan).gibbs_free_energy / AGL_data.natoms << "}";
          for (uint i = firstnonnan+1; i < AGL_data.AGL_pressure_temperature_energy_list.size(); i++) {
            if (!std::isnan(AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy)) {
              oss << ",{\"temperature\":" << AGL_data.AGL_pressure_temperature_energy_list.at(i).temperature_external << ",\"pressure\":" << AGL_data.AGL_pressure_temperature_energy_list.at(i).pressure_external << ",\"gibbs_energy\":" << AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy / AGL_data.natoms << "}";
            }
          }
          oss << "],\"temperature_stepnumber\":" << USER_NTEMPERATURE << ",\"temperature_stepsize\":" << USER_STEMPERATURE << ",\"pressure_stepnumber\":" << USER_NPRESSURE << ",\"pressure_stepsize\":" << USER_SPRESSURE;
          oss << ",\"units\":{\"energy\":\"eV/atom\",\"volume\":\"Angstrom^3/atom\",\"temperature\":\"K\",\"pressure\":\"GPa\"}";
          if (USER_HUGONIOT_CALC.option) {
            oss << ",\"hugoniot_output\": [";
            for (uint i = 0; i < hugoniot_output.size(); i++) {
              if ((hugoniot_output.at(i).size() < 4) || std::isnan(hugoniot_output.at(i).at(0)) || std::isnan(hugoniot_output.at(i).at(1)) || std::isnan(hugoniot_output.at(i).at(2)) || std::isnan(hugoniot_output.at(i).at(3))) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_WARNING_ << "Hugoniot entry " << i << " only has " << hugoniot_output.at(i).size() << " components" <<  endl;
                aus << _AGLSTR_WARNING_ << "It should have 4 components, not writing this entry to JSON" <<  endl;	      
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                continue;
              } else {
                firstnonnan = i;
                break;
              }
            }
            oss << "{\"mass_density\":" << hugoniot_output.at(firstnonnan).at(0) << ",\"temperature\":" << hugoniot_output.at(firstnonnan).at(1) << ",\"energy\":" << hugoniot_output.at(firstnonnan).at(2) << ",\"pressure\":"  << hugoniot_output.at(firstnonnan).at(3) << "}";
            for (uint i = firstnonnan+1; i < hugoniot_output.size(); i++) {
              if ((hugoniot_output.at(i).size() < 4) || std::isnan(hugoniot_output.at(i).at(0)) || std::isnan(hugoniot_output.at(i).at(1)) || std::isnan(hugoniot_output.at(i).at(2)) || std::isnan(hugoniot_output.at(i).at(3))) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_WARNING_ << "Hugoniot entry " << i << " only has " << hugoniot_output.at(i).size() << " components" <<  endl;
                aus << _AGLSTR_WARNING_ << "It should have 4 components, not writing this entry to JSON" <<  endl;	      
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                continue;
              } else {
                oss << ",{\"mass_density\":" << hugoniot_output.at(i).at(0) << ",\"temperature\":" << hugoniot_output.at(i).at(1) << ",\"energy\":" << hugoniot_output.at(i).at(2) << ",\"pressure\":"  << hugoniot_output.at(i).at(3) << "}";
              }
            }
            oss << "],\"hugoniot_units\":{\"mass_density\":\"g/cc\",\"energy\":\"kJ/g\",\"temperature\":\"K\",\"pressure\":\"GPa\"}";
          }
          oss << "}" << endl;
          string ofileenergyjson = AGL_data.dirpathname + "/AGL_energy.json";
          if(!aurostd::stringstream2file(oss, ofileenergyjson, "WRITE")) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energy.json" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }	
          oss.clear();
          oss.str(std::string());
          oss << "# temperature_stepnumber=" << USER_NTEMPERATURE << endl;
          oss << "# temperature_stepsize=" << USER_STEMPERATURE << endl;
          oss << "# pressure_stepnumber=" << USER_NPRESSURE << endl;
          oss << "# pressure_stepsize=" << USER_SPRESSURE << endl;
          oss << "#  T (K) " << "        " << aurostd::PaddedPRE("P (GPa)",12," ") << "        " << aurostd::PaddedPRE("G (eV/atom)",13," ") << endl;
          for (uint i = 0; i < AGL_data.AGL_pressure_temperature_energy_list.size(); i++) {
            if (!std::isnan(AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy)) {
              oss << setw(8) << setprecision(2) << fixed << AGL_data.AGL_pressure_temperature_energy_list.at(i).temperature_external << "          " << setw(11) << setprecision(2) << AGL_data.AGL_pressure_temperature_energy_list.at(i).pressure_external << "          " << setw(11) << setprecision(6) << AGL_data.AGL_pressure_temperature_energy_list.at(i).gibbs_free_energy / AGL_data.natoms << endl;
            }
          }	
          string ofileenergies = AGL_data.dirpathname + "/AGL_gibbs_energy_pT.out";
          if(!aurostd::stringstream2file(oss, ofileenergies, "WRITE")) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Unable to open file AGL_gibbs_energy_pT.out" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }	
          oss.clear();
          oss.str(std::string());
          oss << "# pressure_stepnumber=" << USER_NPRESSURE << endl;
          oss << "# pressure_stepsize=" << USER_SPRESSURE << endl;
          oss << "#  P (GPa) " << "        " << aurostd::PaddedPRE("H (eV/atom)",13," ") << endl;
          for (uint k = 0; k < AGL_data.AGL_pressure_enthalpy_list.size(); k++) {
            if (!std::isnan(AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy)) {
              oss << setw(10) << setprecision(2) << fixed << setprecision(2) << AGL_data.AGL_pressure_enthalpy_list.at(k).pressure_external << "          " << setw(12) << setprecision(6) << AGL_data.AGL_pressure_enthalpy_list.at(k).enthalpy / AGL_data.natoms << endl;
            }
          }	
          string ofileenthalpies = AGL_data.dirpathname + "/AGL_enthalpy_pressure.out";
          if(!aurostd::stringstream2file(oss, ofileenthalpies, "WRITE")) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Unable to open file AGL_enthalpy_pressure.out" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }
          oss.clear();
          oss.str(std::string());			  
          json_written = true;
        } else {
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] energy_pT_kJg.resize(AGL_data.pressure_external.size());
            // [OBSOLETE] oss << "# Data for pressure = " << AGL_data.pressure_external.at(k) << endl;
            // [OBSOLETE] oss << "# rho(g/cc)     T(K)            E(kJ/g)         P(GPa)" << endl;
            for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
              // Convert energy from eV/cell to kJ/g
              energy_kJg = (AGL_data.EnergyDFT_UIntVib.at(j).at(k) * 1.6e-22) / cellmass_grams;
              energy_pT_kJg.at(k).push_back(energy_kJg);
            }
          }
          aglerror = AGL_functions::runHugoniot(AGL_data.pressure_external, AGL_data.temperature_external, AGL_data.mass_density_gcm3, energy_pT_kJg, hugoniot_output, AGL_data.hugoniotextrapolate, FileMESSAGE);
          if(aglerror != 0 && aglerror != 2) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Failure in Hugoniot calculation" <<  endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          } else {
            oss << std::setw(15) << "#     rho(g/cc)" << std::setw(15) << "T(K)" << std::setw(15) << "E(kJ/g)" << std::setw(15) << "P(GPa)" << endl;
            for (uint i = 0; i < hugoniot_output.size(); i++) {
              for (uint j = 0; j < hugoniot_output.at(i).size(); j++) {
                // [OBSOLETE] oss << "\t" << std::left << std::setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
                oss << setw(15) << setprecision(6) << hugoniot_output.at(i).at(j);
              }
              oss << endl;
            }
            string ofilehugoniotoutput = AGL_data.dirpathname + "/AGL_Hugoniot.out";
            if(!aurostd::stringstream2file(oss, ofilehugoniotoutput, "WRITE")) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Unable to open file AGL_Hugoniot.out" <<  endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }	
            // Reset aglerror to 0 on successful output of Hugoniot data
            aglerror = 0;
            // [OBSOLETE] aurostd::StringstreamClean(oss);
          }
          oss.clear();
          oss.str(std::string());
        }
      }
      // If the JSON files hasn't already been written out, then writes calculated energies to JSON format
      if (!json_written) {
        // Write calculated energies to JSON format:
        // E(V) data including equilibrium volume determined by fitting calculated E(V) data
        // Enthalpy as a function of pressure (H = E + pV)
        // Gibbs free energy as a function of temperature and pressure
        // Write calculated energies to JSON format:
        // E(V) data including equilibrium volume determined by fitting calculated E(V) data
        // Gibbs free energy as a function of temperature and pressure
        bool firstgibbs = true;
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "Opening file AGL_energy.json" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        oss << "{\"energy_volume_atom\": [";
        for (uint i = 0; i < volume_atom.size(); i++) {
          if (i==0) {
            oss << "{\"energy\":" << energy_atom.at(i) << ",\"volume\":" << volume_atom.at(i) << "}";
          } else {
            oss << ",{\"energy\":" << energy_atom.at(i) << ",\"volume\":" << volume_atom.at(i) << "}";
          }
        }
        oss << "],\"enthalpy_pressure_atom\": [";
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          if (!std::isnan(AGL_data.EnthalpyPressureeV.at(k))) {
            firstnonnan = k;
            break;
          }
        }
        oss << "{\"pressure\":" << AGL_data.pressure_external.at(firstnonnan) << ",\"enthalpy\":" << AGL_data.EnthalpyPressureeV.at(firstnonnan) / AGL_data.natoms << "}";
        for (uint k = firstnonnan+1; k < AGL_data.pressure_external.size(); k++) {
          if (!std::isnan(AGL_data.EnthalpyPressureeV.at(k))) {
            oss << ",{\"pressure\":" << AGL_data.pressure_external.at(k) << ",\"enthalpy\":" << AGL_data.EnthalpyPressureeV.at(k) / AGL_data.natoms << "}";
          }
        }      
        oss << "],\"gibbs_free_energy_atom\": [";
        for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // [OBSOLETE] if ((j==0) && (k==0))
            if (firstgibbs) {
              if (!std::isnan(AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k))) {
                oss << "{\"temperature\":" << AGL_data.temperature_external.at(j) << ",\"pressure\":" << AGL_data.pressure_external.at(k) << ",\"gibbs_energy\":" << AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k) / AGL_data.natoms << "}";
                firstgibbs = false;
              }
            } else {
              if (!std::isnan(AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k))) {
                oss << ",{\"temperature\":" << AGL_data.temperature_external.at(j) << ",\"pressure\":" << AGL_data.pressure_external.at(k) << ",\"gibbs_energy\":" << AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k) / AGL_data.natoms << "}";
              }
            }
          }
        }
        oss << "],\"temperature_stepnumber\":" << USER_NTEMPERATURE << ",\"temperature_stepsize\":" << USER_STEMPERATURE << ",\"pressure_stepnumber\":" << USER_NPRESSURE << ",\"pressure_stepsize\":" << USER_SPRESSURE;
        oss << ",\"units\":{\"energy\":\"eV/atom\",\"volume\":\"Angstrom^3/atom\",\"temperature\":\"K\",\"pressure\":\"GPa\"}";
        if (USER_HUGONIOT_CALC.option) {
          oss << ",\"hugoniot_output\": [";
          for (uint i = 0; i < hugoniot_output.size(); i++) {
            if ((hugoniot_output.at(i).size() < 4) || std::isnan(hugoniot_output.at(i).at(0)) || std::isnan(hugoniot_output.at(i).at(1)) || std::isnan(hugoniot_output.at(i).at(2)) || std::isnan(hugoniot_output.at(i).at(3))) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ << "Hugoniot entry " << i << " only has " << hugoniot_output.at(i).size() << " components" <<  endl;
              aus << _AGLSTR_WARNING_ << "It should have 4 components, not writing this entry to JSON" <<  endl;	      
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              continue;
            } else {
              firstnonnan = i;
              break;
            }
          }
          oss << "{\"mass_density\":" << hugoniot_output.at(firstnonnan).at(0) << ",\"temperature\":" << hugoniot_output.at(firstnonnan).at(1) << ",\"energy\":" << hugoniot_output.at(firstnonnan).at(2) << ",\"pressure\":"  << hugoniot_output.at(firstnonnan).at(3) << "}";
          for (uint i = 0; i < hugoniot_output.size(); i++) {	      
            if ((hugoniot_output.at(i).size() < 4) || std::isnan(hugoniot_output.at(i).at(0)) || std::isnan(hugoniot_output.at(i).at(1)) || std::isnan(hugoniot_output.at(i).at(2)) || std::isnan(hugoniot_output.at(i).at(3))) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ << "Hugoniot entry " << i << " only has " << hugoniot_output.at(i).size() << " components" <<  endl;
              aus << _AGLSTR_WARNING_ << "It should have 4 components, not writing this entry to JSON" <<  endl;	      
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              continue;
            } else {
              oss << ",{\"mass_density\":" << hugoniot_output.at(i).at(0) << ",\"temperature\":" << hugoniot_output.at(i).at(1) << ",\"energy\":" << hugoniot_output.at(i).at(2) << ",\"pressure\":"  << hugoniot_output.at(i).at(3) << "}";
            }
          }
          oss << "],\"hugoniot_units\":{\"mass_density\":\"g/cc\",\"energy\":\"kJ/g\",\"temperature\":\"K\",\"pressure\":\"GPa\"}";
        }
        oss << "}" << endl;
        string ofileenergyjson = AGL_data.dirpathname + "/AGL_energy.json";
        if(!aurostd::stringstream2file(oss, ofileenergyjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_energy.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
        oss << "# temperature_stepnumber=" << USER_NTEMPERATURE << endl;
        oss << "# temperature_stepsize=" << USER_STEMPERATURE << endl;
        oss << "# pressure_stepnumber=" << USER_NPRESSURE << endl;
        oss << "# pressure_stepsize=" << USER_SPRESSURE << endl;
        oss << "#  T (K) " << "        " << aurostd::PaddedPRE("P (GPa)",12," ") << "        " << aurostd::PaddedPRE("G (eV/atom)",13," ") << endl;
        for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            if (!std::isnan(AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k))) {
              oss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << "          "  << setw(11) << setprecision(2) << AGL_data.pressure_external.at(k) << "          "  << setw(11) << setprecision(6) << AGL_data.GibbsFreeEnergyPressureeV.at(j).at(k) / AGL_data.natoms << endl;
            }
          }
        }	
        string ofileenergies = AGL_data.dirpathname + "/AGL_gibbs_energy_pT.out";
        if(!aurostd::stringstream2file(oss, ofileenergies, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_gibbs_energy_pT.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
        oss << "# pressure_stepnumber=" << USER_NPRESSURE << endl;
        oss << "# pressure_stepsize=" << USER_SPRESSURE << endl;
        oss << "#  P (GPa) " << "        " << aurostd::PaddedPRE("H (eV/atom)",13," ") << endl;
        for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
          if (!std::isnan(AGL_data.EnthalpyPressureeV.at(k))) {
            oss << setw(10) << setprecision(2) << fixed << setprecision(2) << AGL_data.pressure_external.at(k) << "          " << setw(12) << setprecision(6) << AGL_data.EnthalpyPressureeV.at(k) / AGL_data.natoms << endl;
          }
        }	
        string ofileenthalpies = AGL_data.dirpathname + "/AGL_enthalpy_pressure.out";
        if(!aurostd::stringstream2file(oss, ofileenthalpies, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file AGL_enthalpy_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        oss.clear();
        oss.str(std::string());		
      }
    }
    catch (AGLStageBreak& e) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_NOTICE_ + "Stopped. Waiting for completion of required calculations..." << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aglerror = 8;
    }

    return aglerror;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL
// **************************************************************************


#endif  // _AFLOW_AGL_DEBYE_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
