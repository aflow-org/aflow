// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// aflow_agl_debye.h
// Functions written by Cormac Toher
// 2013-2018: cormac.toher@duke.edu
// Hugoniot functions written by Patrick Avery
// 2016-2018: psavery@buffalo.edu

#ifndef _AFLOW_AGL_DEBYE_H_
#define _AFLOW_AGL_DEBYE_H_

// ***************************************************************************
// Numerical constants for use in AGL
#define zero 0.0
#define one 1.0
#define three 3.0
#define third one/three

// Physical constants for use in AGL
// [OBSOLETE] #define hartree2kjmol 2625.5             // hartree -> kJ/mol
// [OBSOLETE] #define au2GPa 29421.4                   // atomic units (pressure) -> GPa
// [OBSOLETE] #define physconstkbau 3.166830e-6        // Boltzmann constant [Hartree/K]
// [OBSOLETE] #define physconstme 9.1093897e-28        // electron mass [g]
// [OBSOLETE] #define physconstamu 1.6605402e-24       // atomic mass unit [g]
#define kj2unit 0.1202707236;                 // kJ/mol -> kB/cell
#define hart2ev 27.211383                     // Hartree -> eV
#define eV_Ang3_to_GPa 160.21766              // eV/Ang^3 -> GPa
#define eV_Ang3_to_amu_Ang_s 9.648533e27      // eV/Ang^3 -> atomic mass unit / (Ang s^2)
#define hbar_eV 6.5821195e-16                 // Planck's constant (angular frequency) hbar = h/2pi in units of eV*s/radian
#define eV2kjmol 96.48533                     // eV/cell to kJ/mol

// Strings for use in read/write operations
#define _AGLSTROPT_ string("[AFLOW_AGL]")
#define _AGIBBSTROPT_ string("[AFLOW_GIBBS]")
#define _AFSTROPT_ string("[AFLOW]")
#define _AELSTROPT_ string("[AFLOW_AEL]")
#define _AGLSTR_DIRNAME_ string("AGL_")
#define _ARAGLSTR_DIRNAME_ string("ARUN.AGL_")
#define _AGLSTR_MESSAGE_ string("00000  MESSAGE AGL: ")
#define _AGLSTR_NOTICE_ string("00000  NOTICE AGL: ")
#define _AGLSTR_WARNING_ string("WWWWW  WARNING AGL: ")
#define _AGLSTR_ERROR_ string("EEEEE  ERROR AGL: ")

using std::string;
using std::vector;
using std::fixed;

struct AGL_pressure_temperature_energies {
  double pressure_external;
  double temperature_external;
  double gibbs_free_energy;
  double volume_equilibrium;
  double internal_energy;
  double E_DFT_internal_vib_energy;
  double helmholtz_vib_energy;
  double vib_entropy;
};

struct AGL_pressure_enthalpies {
  double pressure_external;
  double enthalpy;
};

struct AGL_dos_pressures {
  vector<double> energy;
  vector<double> dosval;
  vector<double> idosval;
  double pressure_external;
  double Fermi_energy;
  double edos_band_gap;
  double dosval_EF;
  double edos_band_gap_polyval;
  double dosval_EF_polyval;
  bool edos_gap_set;
  bool dosval_EF_set;
  bool edos_gap_poly_set;
  bool dosval_EF_poly_set;
};

class AGLStageBreak : public std::exception
{
  public:
    AGLStageBreak() { }
};

// Class to hold input and output data for AGL
class _AGL_data {
  public:
    // constructors/destructors
    _AGL_data();
    ~_AGL_data();
    const _AGL_data& operator=(const _AGL_data &b);    

    // Input title and filename
    string dirpathname;
    string sysname;

    // Stringstreams for output data
    stringstream outfiless, outfiletss;

    // User input data
    double natoms;
    double energy_infinity;
    int i_eqn_of_state;
    int i_debye;
    int i_optimize_beta;
    double poissonratio;
    vector<double> pressure_external;
    vector<double> temperature_external;
    int maxloops;
    int maxfit;
    int maxpolycoeffs;
    int birchfitorder_iG;
    int fittype;
    bool gaussxm_debug;
    // Variable to record origin of Poisson ratio used
    string poissonratiosource;
    // Variable to control whether PREC and ALGO values are set to ACCURATE and NORMAL, or are left as defaults/read from aflow.in file
    bool precaccalgonorm;
    // Variable to control what aflow run type is used (relax, static or relax_static)
    bool relax_static;
    bool static_only;
    bool relax_only;
    // Variable to control running of Hugoniot calculation
    bool hugoniotrun;
    bool hugoniotextrapolate;
    // Number/size of temperature and pressure steps to be used for re-run of postprocessing from command line options
    uint ntemperature;
    double stemperature;
    uint npressure;
    double spressure;


    // Variable to control generation of pressure-volume data for AEL calculations
    bool ael_pressure_calc;

    // Variable to control if this is a postprocessing run (suppresses creation of new ARUN directories)
    bool postprocess;

    // List of failed ARUN directories to be skipped by AGL fitting procedure
    vector<string> failed_arun_list;
    bool autoskipfailedaruns;
    uint skiparunsmax;
    // Initial number of ARUN directories set-up
    // Will be used to check number eliminated in error correction does not exceed skiparunsmax
    uint nstructsinit;

    // Calculated input data
    double cellmass;
    vector<double> energyinput;
    vector<double> volumeinput;
    vector<double> energyinput_orig;
    vector<double> volumeinput_orig;
    vector<double> pressurecalculated;
    vector<xmatrix<double> > stresscalculated;
    vector<int> structurecalculated;
    vector<double> tdebye;
    double pressurecalculatedmax;

    // Data calculated within GIBBS
    vector<double> d2EnergydVolume2_static;
    double poissonratiofunction;
    // Equation of state data
    vector<double> voleqmin;
    vector<double> bulkmodulus;
    vector<double> alpha;
    vector<double> d2EnergydVolume2_dynamic;
    vector<double> gamma_G;
    vector<double> gamma_poly;
    vector<double> pressure_static;
    double rms;
    double bulkmodulus_0pressure;
    double dbulkmodulusdpV_0pressure;
    double d2bulkmodulusdpV2_0pressure;
    // Equilibrium volume and lattice strain factor at zero pressure and temperature
    // Can be used to create a fully relaxed volume structure
    double volume_equilibrium_0p0T;
    double xlattice_equilibrium_0p0T;
    double energy_equilibrium_0p0T;

    // Vector of structs with Gibbs free energy and other quantities as a function of pressure and temperature
    vector<AGL_pressure_temperature_energies> AGL_pressure_temperature_energy_list;
    vector<AGL_pressure_enthalpies> AGL_pressure_enthalpy_list;  
    // Variable to control whether to avoid truncating pressure or temperature range if no minimum energy is found
    // This results in thermodynamic quantities being stored in the vector of structs defined above
    bool run_all_pressure_temperature;

    // Vector of structs with data on electronic properties as a function of pressure
    vector<AGL_dos_pressures> AGL_edos_properties;
    double edosgap_pressure;
    double dosvalEF_pressure;
    double egap_min;
    double egap_max;
    double egap_min_pressure;
    double egap_max_pressure;
    double dosvalEF_min;
    double dosvalEF_max;
    double dosvalEF_min_pressure;
    double dosvalEF_max_pressure;

    // Variable to record whether noise was detected in the E-V data
    bool EV_noise;
    // Variable to record the distance of the minimum of the E-V data from the center
    int itdiff;

    // Saved data for equations for state
    double bcnt_beta;
    double x_K_opt;
    double x_m_opt;
    double bcntvolume0pressure;
    bool optimize_beta;
    vector<double> pfit;
    vector<double> IntEnergStatic;
    double bcnt_beta_statcalc;
    double x_K_opt_statcalc;
    double x_m_opt_statcalc;
    double Vol_sp_statcalc;
    double x_Press_sp_statcalc;
    double volumestatcalc_0pressure;
    double gibbsenergystatcalc_0pressure;
    double bulkmodulusstatcalc_0pressure;
    vector<double> astatic;
    double Avinetstatcalc_0pressure;
    double xsup_K_final;
    double Press_sp_final;
    double Vol_sp_final;
    double bcnt_beta_final;

    // Highest temperature reached
    double max_temperature;

    // zero pressure output data
    vector<double> InternalEnergy0pressurekjmol, Entropy0pressurekjmol, Cvkjmol0pressure, Cpkjmol0pressure, CvunitkB0pressure, CpunitkB0pressure; 
    vector<double> DebyeTemperature0pressure, GruneisenParameter0pressure, HelmholtzEnergy0pressurekjmol, GibbsFreeEnergy0pressurekjmol;
    vector<double> InternalEnergy0pressuremeV, Entropy0pressuremeV, HelmholtzEnergy0pressuremeV, Entropy0pressureunitkB, GibbsFreeEnergy0pressureeV;
    vector<double> ThermalExpansion0pressure, bulkmodulusstatic_0pressure, bulkmodulusisothermal_0pressure;

    // output data for all pressures
    vector<vector<double> > InternalEnergyPressurekjmol, EntropyPressurekjmol, CvkjmolPressure, DebyeTemperaturePressure, GruneisenParameterPressure, HelmholtzEnergyPressurekjmol;
    vector<vector<double> > InternalEnergyPressuremeV, EntropyPressuremeV, HelmholtzEnergyPressuremeV, CvunitkBpressure, EntropyunitkBpressure;
    // [OBSOLETE] vector<vector<double> > GibbsFreeEnergyPressureeV, EnergyDFT_UIntVib, EnergyDFT_UIntVibeV;
    vector<vector<double> > GibbsFreeEnergyPressureeV, EnergyDFT_UIntVib;
    vector<vector<double> > xminsav, VolumeEquilibrium, mass_density_gcm3;

    // Enthalpy at finite pressure (no vibrational energy)
    vector<double> EnthalpyPressureeV;

    // Volume and pressure data for static pressure
    vector<double> VolumeStaticPressure, StaticPressure, VolumeFactors;
    bool savedatapressure;
  private:                                                //
    void free();                                           // free space
};

// Declaration of functions in AGL

// Namespace for functions used by AGL method, to avoid potential naming conflicts with other parts of AFLOW
namespace AGL_functions {
  // Functions to actually run AGL, either directly or from another part of AFLOW
  uint RunDebye_AGL(_xvasp& xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint AGL_xvasp_flags_populate(_xvasp& xvasp, string& AflowIn, string& AflowInName, string& FileLockName, const string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE);
  uint AGL_Get_AflowInName(string& AflowInName, const string& directory_LIB, bool& agl_aflowin_found); //CT20200713 Find aflow.in filename
  uint Get_ThermalProperties_AGL_postprocess(const string& directory, uint ntemperature, double stemperature, uint npressure, double spressure, vector<double>& agl_temperature, vector<double>& agl_gibbs_energy_atom, vector<double>& agl_vibrational_energy_atom);
  uint Get_EquilibriumVolumeTemperature(_xvasp& xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE);
  uint Get_EquilibriumVolumeAngstromTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE);
  uint Get_BulkModulusStaticTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusStatic, ofstream& FileMESSAGE);
  uint Get_BulkModulusIsothermalTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusIsothermal, ofstream& FileMESSAGE);
  uint Get_BulkModulusVolumeTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusIsothermal, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE);
  uint Get_BulkModulusVolumeAngstromTemperature(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Temperature, vector<double>& BulkModulusIsothermal, vector<double>& EquilibriumVolume, ofstream& FileMESSAGE);
  uint Get_VolumeStaticPressure(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<double>& Pressure, vector<double>& PressureVolumes, vector<double>& VolumeScaleFactors, bool postprocess, ofstream& FileMESSAGE);
  // Functions for generating _AFLOWIN_ input files for strained structures and checking (E,V) data calculated with VASP
  uint aglvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, string& dirrunname, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  // [OBSOLETE] uint createAFLOWIN(_xvasp& vaspRun, _xvasp& xvasp, _kflags& kbinFlags, _vflags& vaspFlags, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint extractenerg(vector<_xvasp>& vaspRuns, _AGL_data& AGL_data, vector<string>& dirrunname, ofstream& FileMESSAGE);
  uint checkmin(_AGL_data& AGL_data, int& cmerr, ofstream& FileMESSAGE);
  uint checkconcav(_AGL_data& AGL_data, int& ccerr, ofstream& FileMESSAGE);
  uint gibbsinpwrite(_AGL_data& AGL_data, ofstream& FileMESSAGE);
  // Functions for implementing GIBBS method for calculating thermal properties
  uint gibbsrun(_AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint qcksortev(vector<double>& vol, vector<double>& energ, ofstream& FileMESSAGE);
  uint qcksortevt(vector<double>& vol, vector<double>& energ, vector<double>& tdebye, ofstream& FileMESSAGE);
  uint polynom_eval (double& xval, vector<double>& polynomialcoeffs, double& yval, uint n);
  uint gaussxm (vector<vector<double> >& cmatrix, int& n, vector<double>& xdata, bool gxmdebug, ofstream& FileMESSAGE);
  uint polynom_fit (uint& first_entry_tofit, uint& last_entry_tofit, vector<double>& xdata, vector<double>& ydata, vector<double>& weight, double& rms, int& npolycoeffwork, vector<double>& polycoeffwork, bool gxmdebug, ofstream& FileMESSAGE);
  uint polynomial_fit_weight_ave (uint& imin, vector<double>& polynomialcoeffs, vector<double>& polynomialerrors, vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint polynomial_fit_weight_ave_debug (uint& imin, vector<double>& polynomialcoeffs, vector<double>& polynomialerrors, vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint bracket_minimum (uint& imin, vector<double>& ydata_to_fit, ofstream& FileMESSAGE);
  uint bracket_minimum_global (uint& imin, vector<double>& ydata_to_fit, ofstream& FileMESSAGE);
  uint polynom_bracket_minimum (uint& imin, vector<double>& x, vector<double>& polynomialcoeffs, ofstream& FileMESSAGE);
  bool val_bracketed_check (double x_a, double x_b);
  bool zero_bracketed_check (double& f_a, double& f_b);
  uint polynom_minimum (double& xini, double& xmin, double& xmax, vector<double>& polynomialcoeffs, double& pxmin, ofstream& FileMESSAGE);
  uint autocorrect_polynom_minimum (uint& itry, vector<double>& xvalues, vector<double>& polynomialcoeffs, double& xmin, ofstream& FileMESSAGE);
  uint numerical_eos (double& volumereference, vector<double>& epol, vector<double>& fpol,  vector<double>& xconfigurationvector, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint numerical_eos_run_all_pT (double& volumereference, vector<double>& energypolynomialcoeffs, vector<double>& helmholtzenergypolynomialcoeffs, double& voleq, double& bulkmodpres, double& pres_static, double& d2EnergydVolume2_dynamic_pres, double& Gruneisen_pres);
  uint vinet_eos (double& volume_0pressure, double& gibbsenergy_0pressure, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint birch_murnaghan_eos (double& volume_0pressure, double& gibbsenergy_0pressure, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint bcnt_optimize_beta (double& x_Press_sp, _AGL_data& AGL_data, double& lB_lK_lbppsp, ofstream& FileMESSAGE);
  uint bcnt_bracket_minimum (double& x_a, double& x_b, double& x_c, double& f_a, double& f_b, double& f_c, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint bcnt_brent_minimum (double& x_a, double& x_b, double& x_c, double& tol, double& xmin, _AGL_data& AGL_data, double& brentx, ofstream& FileMESSAGE);
  uint gauss_legendre (double& xupper, double& xlower, vector<double>& xabscissa, vector<double>& weight, int& n, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint bcnt_eos (double& volume_0pressure, double& gibbsenergy_0pressure, double& bulkmod_0pressure, bool statcalc,  _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint debye_polynom_fit (vector<double>& polynomialcoeffs, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint self_consistent_debye (double& Temperature, vector<double>& polynomialcoeffs, vector<double>& polynomialerror, vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, uint& imin, bool firsttime, vector<double>& press_stat, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  double debye_function (double& z);
  uint thermal_properties (double& ThetaD, double& Temperature, double& DebyeIntegral, double& DebyeIntegralerror, double& internalenergy, double& Cv, double& helmholtzenergy, double& entropy, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  // Functions for analysing and plotting results generated by GIBBS method
  uint cvdebfit (double& tdmin, double& tdmax, double& nkb, int& npoints, double& tdbest, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  double fdebyeint (double z);
  uint debyeint (double& yval, double& Debye_int_val, _AGL_data& AGL_data, ofstream& FileMESSAGE);
  uint plotaglresults(vector<double>& xdata, vector<double>& ydata, double& xmin, double& xmax, double& ymin, double& ymax, int& nplotpoints, string& datalabel, string& xlabel, string& ylabel, string& plotname, string& plotfilename, string& sysname, ofstream& FileMESSAGE);
  uint thermalconductD (double& acousticthetaD, vector<double>& temperature, double& kappaD, vector<double>& kappaT, double& gammaD, double& voleqD, double& avmass);
  uint thermalconductDvol (double& acousticthetaD, vector<double>& temperature, vector<double>& kappaTV, double& gammaD, vector< vector<double> >& volumeq, double& avmass);
  // Functions for calculating Hugoniot relation
  // Written by Patrick Avery: psavery@buffalo.edu
  uint runHugoniot(const std::vector<double>& pressures_external, const std::vector<double>& temperatures_external, const std::vector<std::vector<double> >& mass_densities_gcm3, const std::vector<std::vector<double> >& energies_pT_kJg, std::vector<std::vector<double> >& hugoniot_output, bool hugoniotextrapolate, ofstream& FileMESSAGE, double desired_initial_pressure_external = 0.0, double desired_initial_temperature_external = 300.0);
  uint runHugoniotAllTemperaturesAndPressures(const std::vector<AGL_pressure_temperature_energies>& AGL_pressure_temperature_energy_list, double cellmass_grams, std::vector<std::vector<double> >& results, bool hugoniotextrapolate, ofstream& FileMESSAGE, double desired_initial_pressure_external = 0.0, double desired_initial_temperature_external = 300.0);
  // [OBSOLETE] uint calculateHugoniot(const std::vector<std::vector<std::vector<double> > >& data, std::vector<std::vector<double> >& hugoniotData, ofstream& FileMESSAGE, double desired_initial_pressure_external = 0.0, double desired_initial_temperature_external = 300.0);
  // OSBOLETE uint calculateHugoniotDataPoint(double initial_mass_density_gcm3, double initial_energyDFT_UIntVib_kJg, double initial_pressure_external, const std::vector<double>& mass_densities_gcm3, const std::vector<double>& temperatures_external, const std::vector<double>& energiesDFT_UIntVib_kJg, double pressure_external, double& hugoniotDensity, double& hugoniotTemperature, double& hugoniotEnergy, ofstream& FileMESSAGE);
  uint calculateHugoniotDataPoint(double initial_mass_density_gcm3, double initial_energyDFT_UIntVib_kJg, double initial_pressure_external, const std::vector<double>& mass_densities_gcm3, const std::vector<double>& temperatures_external, const std::vector<double>& energiesDFT_UIntVib_kJg, double pressure_external, double& hugoniotDensity, double& hugoniotTemperature, double& hugoniotEnergy, bool hugoniotextrapolate, ofstream& FileMESSAGE);
  // [OBSOLETE] uint findBestInitialConditions(const std::vector<std::vector<std::vector<double> > >& data, double desired_initial_pressure_external, double desired_initial_temperature_external, double& initial_mass_density_gcm3, double& initial_temperature_external, double& initial_energyDFT_UIntVib_kJg, double& initial_pressure_external, ofstream& FileMESSAGE);
  uint findBestInitialConditions(const std::vector<double>& pressures_external, const std::vector<double>& temperatures_external, const std::vector<std::vector<double> >& mass_densities_gcm3, const std::vector<std::vector<double> >& energies_pT_kJg, double desired_initial_pressure_external, double desired_initial_temperature_external, double& initial_mass_density_gcm3, double& initial_temperature_external, double& initial_energyDFT_UIntVib_kJg, double& initial_pressure_external, ofstream& FileMESSAGE);
  uint findBestInitialConditionsAllTemperaturesAndPressures(const std::vector<AGL_pressure_temperature_energies>& AGL_pressure_temperature_energy_list, const std::vector<double>& mass_densities_gcm3, const std::vector<double>& energies_pT_kJg, double desired_initial_pressure_external, double desired_initial_temperature_external, double& initial_mass_density_gcm3, double& initial_temperature_external, double& initial_energyDFT_UIntVib_kJg, double& initial_pressure_external, ofstream& FileMESSAGE);
  void interpolateToFindHugoniotEnergyAndDensity(double initial_mass_density_gcm3, double initial_energyDFT_UIntVib_kJg, double initial_pressure_external, double pressure_external, double energy1, double energy2, double density1, double density2, double& hugoniot_energy, double& hugoniot_density);
  // Functions for extracting electronic properties as a function of pressure
  uint extractedos(vector<_xvasp>& vaspRuns, _AGL_data& AGL_data, vector<string>& dirrunname, ofstream& FileMESSAGE);
  uint edosbandgap(_AGL_data& AGL_data, ofstream& FileMESSAGE); 
  uint edosfit(vector<double>& energtofit, vector<double>& dostofit, vector<double>& energtoeval, vector<double>& dospolyeval, uint& nminlimit, uint& nmaxlimit, bool gxmdebug, ofstream& FileMESSAGE);
  uint edosgap_pressure_fit(vector<double>& pressure, vector<double>& edosgap, double& edosgap_pressure, vector<double>& edosgapfit, uint& nminlimit, uint& nmaxlimit, bool gxmdebug, ofstream& FileMESSAGE);
} // end of namespace

// ***************************************************************************

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *              AFlow CORMAC TOHER - Duke University 2013-2021             *
// *                                                                         *
// ***************************************************************************

