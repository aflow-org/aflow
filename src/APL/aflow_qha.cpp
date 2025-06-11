// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Andriy Smolyanyuk (andriy.smolyanyuk@duke.edu)
//
//    This file is part of Aflow.
//
//    Aflow is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Aflow is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Aflow.  If not, see <https://www.gnu.org/licenses/>.

#include "aflow_apl.h"

// Field width modifiers: used to format text during the output via std::setw function
#define SW 5  // width of columns with blank space separator
#define TW 15 // width of columns containing label/number
#define SWTM 2  // width of columns with blank space separator for the THERMO_MESH block
#define TWTM 18 // width of columns containing label/number for the THERMO_MESH block
#define DCOEFF 1e-2 // this coefficient is used in numerical differentiation
// (should not be too small). Usage dT = DCOEFF*T

// to avoid division by zero in the self-consistent loop,
// the initial volume is taken to be 10% bigger
#define SCQHA_INITIAL_VOLUME_FACTOR 1.1

#define DEBUG_QHA false // toggles QHA-related debug output
#define DEBUG_QHA_GP_FIT false // toggles debug output related to the fit functionality
// in the calcGrueneisen function

// files required by static DFT calculation
enum ST_DATA_FILE {ST_DF_DIRECTORY, ST_DF_OUTCAR, ST_DF_EIGENVAL, ST_DF_IBZKPT};
// files required by APL calculation
enum PH_DATA_FILE {PH_DF_DIRECTORY, PH_DF_HARMIFC, PH_DF_PHPOSCAR};

#define EOS_METHOD_FILE_SJ "stabilized-jellium."
#define EOS_METHOD_FILE_BIRCH_MURNAGHAN2 "birch-murnaghan2."
#define EOS_METHOD_FILE_BIRCH_MURNAGHAN3 "birch-murnaghan3."
#define EOS_METHOD_FILE_BIRCH_MURNAGHAN4 "birch-murnaghan4."
#define EOS_METHOD_FILE_MURNAGHAN "murnaghan."

// labels to specify what contributions to the free energy are included
#define F_ELEC 1
#define F_VIB  2

//================================================================================
//                    EOS related

/// Murnaghan equation of state
/// https://en.wikipedia.org/wiki/Murnaghan_equation_of_state
/// Ref: Proceedings of the National Academy of Sciences of the United States of 
/// America, 30 (9): 244–247
/// 
double Murnaghan(const double x, const xvector<double> &p, xvector<double> &dydp)
{
  double Eeq = p[1];
  double Veq = p[2];
  double B   = p[3];
  double Bp  = p[4];
  double V = x;

  dydp[1] = 1.0;
  dydp[2] = (B*(Veq-V*pow(Veq/V,Bp)))/(Veq-Bp*Veq);
  dydp[3] = -(Veq/(Bp-1.0))+(V*(1+pow(Veq/V,Bp)/(Bp-1.0)))/Bp;
  dydp[4] = (B*(-V+2.0*Bp*V-pow(Bp,2)*V+pow(Bp,2)*Veq+V*pow(Veq/V,Bp)*(1.0-2.0*Bp
          +(Bp-1.0)*Bp*log(Veq/V))))/(pow(Bp-1.0,2)*pow(Bp,2));

  return Eeq-(B*Veq)/(Bp-1.0)+(B*V*(1+pow(Veq/V,Bp)/(Bp-1)))/Bp;
}

/// Birch-Murnaghan equation of state
/// https://en.wikipedia.org/wiki/Birch%E2%80%93Murnaghan_equation_of_state
/// Ref: Physical Review. 71 (11): 809–824
/// 
double BirchMurnaghan(const double x, const xvector<double> &p, xvector<double> &dydp)
{
  double Eeq = p[1];
  double Veq = p[2];
  double B   = p[3];
  double Bp  = p[4];
  double V = x;

  double kappa = pow(Veq/V, 2.0/3.0);

  dydp[1] = 1.0;
  dydp[2] = 3.0*B*kappa*(kappa-1)*(3.0*(Bp-6.0)*V*pow(Veq/V,1.0/3.0)+
      Veq*(62.0-36.0*kappa+3*Bp*(3*kappa-4)))/(16.0*Veq);
  dydp[3] = 9.0/16.0*Veq*pow(kappa-1,2.0)*(6.0-4.0*kappa+Bp*(kappa-1));
  dydp[4] = 9.0/16.0*B*Veq*pow(kappa-1,3.0);

  return Eeq + 9.0/16.0*B*Veq*(pow(kappa-1.0,3.0)*Bp + pow(kappa-1.0,2.0)*(6.0-4.0*kappa));
}

/// Checks if there is a minimum within a given data set 
/// (at least one internal point should be lower than the edge points)
bool isMinimumWithinBounds(const xvector<double> &y){
  for (int i=y.lrows+1; i<y.urows; i++){
    if (y[y.lrows] > y[i] && y[y.urows] > y[i]) return true;
  }

  return false;
}

/// Calculates the bulk modulus for the Birch-Murnaghan EOS model.
double calcBulkModulus_BM(double x_eq, const xvector<double> &dEdp)
{
  double B = 4.0/9.0*std::pow(x_eq,5)*dEdp[3]+10.0/9.0*std::pow(x_eq,4)*dEdp[2];
  B *= eV2GPa*std::pow(x_eq,-3.0/2.0);
  return B;
}

/// Calculates the pressure derivative of the bulk modulus for the Birch-Murnaghan EOS
/// model.
double Bprime_BM(double x, const xvector<double> &dEdp)
{
  // -(V*d^3E/dV^3/d^2E/dV^2 + 1)
  double d3EdV3 = -8.0/27.0*std::pow(x,15.0/2.0)*dEdp[4];
  d3EdV3 -= 20.0/9.0*std::pow(x, 13.0/2.0)*dEdp[3];
  d3EdV3 -= 80.0/28.0*std::pow(x, 11.0/2.0)*dEdp[2];

  double d2EdV2 = 4.0/9.0*std::pow(x, 5)*dEdp[3];
  d2EdV2 += 10.0/9.0*std::pow(x, 4)*dEdp[2];

  return -(std::pow(x,-3.0/2.0)*d3EdV3/d2EdV2 + 1);
}

/// Calculates the bulk modulus for the stabilized jellium EOS model.
/// Stabilized jellium model: https://doi.org/10.1103/PhysRevB.63.224115
/// It is equivalent to the model introduced in:
/// https://doi.org/10.1103/PhysRevB.52.8064
double calcBulkModulus_SJ(double x, const xvector<double> &dEdp)
{
  double B = 1.0/9.0*std::pow(x,8)*dEdp[3]+4.0/9.0*std::pow(x,7)*dEdp[2];
  B *= eV2GPa*std::pow(x,-3);
  return B;
}

/// Calculates the pressure derivative of the bulk modulus for the stabilized jellium EOS
/// model.
/// Stabilized jellium model: https://doi.org/10.1103/PhysRevB.63.224115
/// It is equivalent to the model introduced in:
/// https://doi.org/10.1103/PhysRevB.52.8064
double Bprime_SJ(double x, const xvector<double> &dEdp)
{
  // -(V*d^3E/dV^3/d^2E/dV^2 + 1)
  double d3EdV3 = -1.0/27.0*std::pow(x,12)*dEdp[4];
  d3EdV3 -= 4.0/9.0*std::pow(x,11)*dEdp[3];
  d3EdV3 -= 28.0/27.0*std::pow(x,10)*dEdp[2];

  double d2EdV2 = 1.0/9.0*std::pow(x, 8)*dEdp[3];
  d2EdV2 += 4.0/9.0*std::pow(x, 7)*dEdp[2];

  return -(std::pow(x,-3)*d3EdV3/d2EdV2 + 1);
}

//=============================================================================
//                         Definitions of the QHA class members
namespace apl
{
  QHA::QHA(ostream& oss) : xStream(oss) { free(); }
  QHA::QHA(const QHA &qha) : xStream(*qha.getOFStream(),*qha.getOSS()) {
    free(); copy(qha);
  }
  QHA::~QHA() { xStream::free(); free(); }
  const QHA& QHA::operator=(const QHA &qha){
    copy(qha);
    return *this;
  }

  void QHA::clear()
  {
    free();
  }

  /// Initializes all values to "zero" and attempts to clear all containers.
  ///
  void QHA::free()
  {
    apl_options.clear();
    qha_options.clear();
    system_title = "";
    isEOS = false; isGP_FD = false;
    ignore_imaginary = false;
    isQHA   = false; isQHA3P = false; isSCQHA = false; isQHANP = false;
    isInitialized = false;
    includeElectronicContribution = false;
    doSommerfeldExpansion = false;
    Ntemperatures = 0;
    N_GPvolumes = 3;
    N_EOSvolumes = 0;
    N_QHANPvolumes = 0;
    Nbranches = 0;
    NatomsOrigCell = 0;
    Nelectrons = 0;
    TaylorExpansionOrder = 0;
    gp_distortion = 0.0;
    origStructure.clear();
    Temperatures.clear();
    ph_disp_temperatures.clear();
    GPvolumes.clear();
    EOSvolumes.clear();
    QHANPvolumes.clear();
    coefGPVolumes.clear();
    coefEOSVolumes.clear();
    coefQHANPVolumes.clear();
    DOS_Ef.clear();
    Efermi_V.clear();
    E0_V.clear();
    static_eigvals.clear();
    static_ibzkpts.clear();
    qpWeights.clear();
    qPoints.clear();
    gp_fit_matrix.clear();
    omegaV_mesh.clear();
    omegaV_mesh_EOS.clear();
    gp_ph_dispersions.clear();
    eos_ph_dispersions.clear();
    eos_vib_thermal_properties.clear();
    subdirectories_apl_gp.clear();
    subdirectories_apl_eos.clear();
    subdirectories_apl_qhanp.clear();
    subdirectories_static.clear();
    arun_runnames_apl_eos.clear();
    arun_runnames_apl_gp.clear();
    arun_runnames_apl_qhanp.clear();
    arun_runnames_static.clear();
    xinput.clear();
    currentDirectory = ".";
  }

  void QHA::copy(const QHA &qha)
  {
    if (this==&qha) return;

    apl_options       = qha.apl_options;
    qha_options       = qha.qha_options;
    system_title      = qha.system_title;
    isEOS             = qha.isEOS;
    isGP_FD           = qha.isGP_FD;
    ignore_imaginary  = qha.ignore_imaginary;
    isQHA             = qha.isQHA;
    isQHA3P           = qha.isQHA3P;
    isSCQHA           = qha.isSCQHA;
    isQHANP           = qha.isQHANP;
    isInitialized     = qha.isInitialized;
    includeElectronicContribution = qha.includeElectronicContribution;
    doSommerfeldExpansion = qha.doSommerfeldExpansion;
    Ntemperatures     = qha.Ntemperatures;
    N_GPvolumes       = qha.N_GPvolumes;
    N_EOSvolumes      = qha.N_EOSvolumes;
    N_QHANPvolumes    = qha.N_QHANPvolumes;
    Nbranches         = qha.Nbranches;
    NatomsOrigCell    = qha.NatomsOrigCell;
    Nelectrons        = qha.Nelectrons;
    TaylorExpansionOrder = qha.TaylorExpansionOrder;
    gp_distortion     = qha.gp_distortion;
    origStructure     = qha.origStructure;
    Temperatures      = qha.Temperatures;
    ph_disp_temperatures = qha.ph_disp_temperatures;
    GPvolumes         = qha.GPvolumes;
    EOSvolumes        = qha.EOSvolumes;
    QHANPvolumes      = qha.QHANPvolumes;
    coefGPVolumes     = qha.coefGPVolumes;
    coefEOSVolumes    = qha.coefEOSVolumes;
    coefQHANPVolumes  = qha.coefQHANPVolumes;
    DOS_Ef            = qha.DOS_Ef;
    Efermi_V          = qha.Efermi_V;
    E0_V              = qha.E0_V;
    static_eigvals    = qha.static_eigvals;
    static_ibzkpts    = qha.static_ibzkpts;
    qpWeights         = qha.qpWeights;
    qPoints           = qha.qPoints;
    gp_fit_matrix     = qha.gp_fit_matrix;
    omegaV_mesh       = qha.omegaV_mesh;
    omegaV_mesh_EOS   = qha.omegaV_mesh_EOS;
    gp_ph_dispersions = qha.gp_ph_dispersions;
    eos_ph_dispersions = qha.eos_ph_dispersions;
    eos_vib_thermal_properties = qha.eos_vib_thermal_properties;
    subdirectories_apl_eos  = qha.subdirectories_apl_eos;
    subdirectories_apl_gp   = qha.subdirectories_apl_gp;
    subdirectories_apl_qhanp   = qha.subdirectories_apl_qhanp;
    subdirectories_static   = qha.subdirectories_static;
    arun_runnames_apl_eos   = qha.arun_runnames_apl_eos;
    arun_runnames_apl_gp    = qha.arun_runnames_apl_gp;
    arun_runnames_apl_qhanp = qha.arun_runnames_apl_qhanp;
    arun_runnames_static    = qha.arun_runnames_static;
    xinput            = qha.xinput;
    currentDirectory  = qha.currentDirectory;

    xStream::copy(qha);
  }

  ///////////////////////////////////////////////////////////////////////////////

  QHA::QHA(const xstructure &in_structure, _xinput &xinput, xoption &qha_options,
      xoption &apl_options, ofstream &FileMESSAGE, ostream &oss)
  {
    initialize(in_structure, xinput, qha_options, apl_options, FileMESSAGE, oss);
  }

  //xStream initializers
  void QHA::initialize(ostream& oss) {
    xStream::initialize(oss);
  }

  void QHA::initialize(ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
  }

  /// Initializes the QHA class with all the necessary data.
  ///
  void QHA::initialize(const xstructure &in_structure, _xinput &xinput,
      xoption &qha_options, xoption &apl_options, ofstream &FileMESSAGE, ostream &oss)
  {
    static const int precision_format = 4;

    isInitialized = false;

    xStream::initialize(FileMESSAGE, oss);

    string msg  = "Initializing QHA.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    free();

    this->xinput = xinput;
    this->apl_options = apl_options;
    this->qha_options = qha_options;

    currentDirectory = xinput.xvasp.Directory; // remember the current directory

    // initialization of the data that is derived from the structure
    origStructure = in_structure;
    origStructure.ReScale(1.0);
    NatomsOrigCell = origStructure.atoms.size();
    Nbranches = NatomsOrigCell * 3;
    double Volume = origStructure.GetVolume();

    // parse QHA-related aflow.in options
    isEOS = qha_options.flag("EOS");
    vector<double> eosrange(3);
    aurostd::string2tokens(qha_options.getattachedscheme("EOS_DISTORTION_RANGE"),
        eosrange, " :");
    gp_distortion = aurostd::string2utype<double>(qha_options.getattachedscheme("GP_DISTORTION"));
    gp_distortion /= 100.0;
    isGP_FD = qha_options.flag("GP_FINITE_DIFF");
    ignore_imaginary = qha_options.flag("IGNORE_IMAGINARY");
    includeElectronicContribution = qha_options.flag("INCLUDE_ELEC_CONTRIB");
    isQHA = qha_options.flag("MODE:QHA");
    isQHA3P = qha_options.flag("MODE:QHA3P");
    isQHANP = qha_options.flag("MODE:QHANP");
    isSCQHA = qha_options.flag("MODE:SCQHA");
    aurostd::string2tokens(qha_options.getattachedscheme("PDIS_T"), ph_disp_temperatures, ",");
    doSommerfeldExpansion = qha_options.flag("SOMMERFELD_EXPANSION");
    TaylorExpansionOrder = aurostd::string2utype<int>(qha_options.getattachedscheme("TAYLOR_EXPANSION_ORDER"));

    // retrieve the minimum number of parameters for the EOS calculation with a
    // given set of requested EOS models
    int REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT = 5;
    if (qha_options.flag("EOS_MODEL:BM2")){
      REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT = 3;
    }

    if (qha_options.flag("EOS_MODEL:BM3") ||
        qha_options.flag("EOS_MODEL:SJ")  ||
        qha_options.flag("EOS_MODEL:M")){
      REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT = 4;
    }

    if (qha_options.flag("EOS_MODEL:BM4")){
      REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT = 5;
    }

    // output with what parameters QHA will be run
    stringstream message;
    message << "QHA will run with the following parameters:" << std::endl;
    message << "Methods:" << (isQHA ? " QHA" : "") << (isQHA3P ? " QHA3P" : "");
    message << (isQHANP ? " QHANP" : "") << (isSCQHA ? " SCQHA" : "") << "." << std::endl;

    message << "Finite differences calculation of Grueneisen parameters will ";
    message << (isGP_FD ? "" : "NOT ") << "be run." << std::endl;

    message << "Equation of state will " << (isEOS ? "" : "NOT ") << "be calculated ";
    message << "when QHA method is used." << std::endl;
    message << "EOS range of volumes and volume increment is ";
    message << qha_options.getattachedscheme("EOS_DISTORTION_RANGE") << "." << std::endl;

    message << "Electronic contribution to the free energy will ";
    message << (includeElectronicContribution ? "" : "NOT ") << "be included." << std::endl;
    message << "Sommerfeld expansion will " << (doSommerfeldExpansion ? "" : "NOT ");
    message << "be employed." << std::endl;

    message << "GP_DISTORTION = " << gp_distortion << "." << std::endl;
    message << "TAYLOR_EXPANSION_ORDER = " << TaylorExpansionOrder << "." << std::endl;
    message << "Temperature-dependent phonon dispersions will be calculated for the";
    message << " following list of temperatures: ";
    message << qha_options.getattachedscheme("PDIS_T") << "." << std::endl;

    message << "If there are unstable phonon modes, they will ";
    message << (ignore_imaginary ? "" : "NOT ") << "be ignored.";

    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, message, currentDirectory, *p_FileMESSAGE,
        *p_oss, _LOGGER_MESSAGE_);

    // determine the names for the directories for the calculation of the Grueneisen parameter
    // (calculated using finite differences method)
    vector<double> gprange(3);
    gprange[0] = 1.0-gp_distortion; gprange[1] = 1.0; gprange[2] = 1.0+gp_distortion;
    N_GPvolumes = 3;
    for (int j=0; j<N_GPvolumes; j++){
      arun_runnames_apl_gp.push_back("PHONON_"
          + aurostd::utype2string(gprange[j],precision_format,false,FIXED_STREAM));

      // if APL state file for V=1.0 exists, reuse it and do not create PHONON_1.000 directory
      if (aurostd::isequal(gprange[j], 1.0) &&
          aurostd::EFileExist(currentDirectory + "/" +
            DEFAULT_APL_FILE_PREFIX+DEFAULT_APL_STATE_FILE)){
        subdirectories_apl_gp.push_back(".");
      }
      else{
        subdirectories_apl_gp.push_back(ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
            arun_runnames_apl_gp.back());
      }

      coefGPVolumes.push_back(gprange[j]);
      GPvolumes.push_back(gprange[j]*Volume/NatomsOrigCell);
    }

    // determine the names for the directories used for the EOS calculation
    if (isEOS || isQHA3P || isSCQHA || isQHANP){
      eosrange[0] = 1.0 + eosrange[0]/100.0;
      eosrange[1] = 1.0 + eosrange[1]/100.0;
      eosrange[2] = eosrange[2]/100.0;

      // get a set of volumes that would be used for the QHA-EOS calculation
      string dirname = "";
      for (double i=eosrange[0]; i<=(eosrange[1]+AUROSTD_IDENTITY_TOL); i+=eosrange[2]){
        arun_runnames_apl_eos.push_back( "PHONON_" +
            aurostd::utype2string(i,precision_format,false,FIXED_STREAM));

        // if APL state file for V=1.0 exists, reuse it and do not create PHONON_1.000 directory
        if (aurostd::isequal(i, 1.0) &&
            aurostd::EFileExist(currentDirectory + "/" +
              DEFAULT_APL_FILE_PREFIX+DEFAULT_APL_STATE_FILE)){
          subdirectories_apl_eos.push_back(".");
        }
        else{
          subdirectories_apl_eos.push_back(ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
              arun_runnames_apl_eos.back());
        }

        arun_runnames_static.push_back("STATIC_" + 
            aurostd::utype2string(i, precision_format, false, FIXED_STREAM));
        dirname = ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
          arun_runnames_static.back();
        subdirectories_static.push_back(dirname);

        coefEOSVolumes.push_back(i);
        EOSvolumes.push_back(i*Volume/NatomsOrigCell);
      }
      N_EOSvolumes = EOSvolumes.size();

      // check that we have enough datapoints to do the requested set of EOS
      // models. If not print an error and swith the EOS calculation off
      if (N_EOSvolumes < REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT){
        isEOS = false;

        stringstream msg;
        msg << "QHA EOS calculation requires at least ";
        msg << aurostd::utype2string<int>(REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT);
        msg << " APL calculations for the requested set of EOS models." << std::endl;
        msg << "The current choice of volume range and increment produces ";
        msg << aurostd::utype2string<int>(N_EOSvolumes);
        msg << " APL calculations." << std::endl;
        msg << "QHA EOS calculation will be skipped!";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_ERROR_);
      }
    }

    // determine the names for the directories used for the QHANP calculation
    if (isQHANP){
      N_QHANPvolumes = 2*TaylorExpansionOrder + 1;

      // get a set of volumes that would be used for the QHANP calculation
      double coef = 0.0;
      for (int i=0; i<N_QHANPvolumes; i++){
        coef = 1.0 + (-TaylorExpansionOrder + i)*gp_distortion; //dV = 2*gp_distortion
        arun_runnames_apl_qhanp.push_back("PHONON_"+
            aurostd::utype2string(coef,precision_format,false,FIXED_STREAM));

        // if APL state file for V=1.0 exists, reuse it and do not create PHONON_1.000 directory
        if (aurostd::isequal(coef, 1.0) &&
            aurostd::EFileExist(currentDirectory + "/" +
              DEFAULT_APL_FILE_PREFIX+DEFAULT_APL_STATE_FILE)){
          subdirectories_apl_qhanp.push_back(".");
        }
        else{
          subdirectories_apl_qhanp.push_back(ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
              arun_runnames_apl_qhanp.back());
        }

        coefQHANPVolumes.push_back(coef);
        QHANPvolumes.push_back(coef*Volume/NatomsOrigCell);
      }
    }


    vector<double> tokens;
    aurostd::string2tokens(apl_options.getattachedscheme("TPT"), tokens, string (" :"));
    if (tokens.size() != 3){
      stringstream msg;
      msg << "Wrong setting in ";
      msg << "[AFLOW_APL]TPT.";
      msg << "Specify as TPT="+AFLOWRC_DEFAULT_APL_TPT+"." << std::endl;
      msg << "See README_AFLOW_APL.TXT for the details.";
      throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg,
          _INPUT_NUMBER_);
    }
    double tp_start = tokens[0];
    double tp_end   = tokens[1];
    double tp_step  = tokens[2];


    // define a set of temperatures for thermodynamic calculations
    Ntemperatures = floor((tp_end - tp_start)/tp_step) + 1;

    for (int T=tp_start; T<=tp_end; T+=tp_step) Temperatures.push_back(T);

    // this matrix will be used in the fit of the frequency-volume dependency
    // w(V) = a + b*V + c*V^2 + d*V^3
    gp_fit_matrix = xmatrix<double> (N_EOSvolumes, 4);
    for (int Vid=1; Vid<=N_EOSvolumes; Vid++){
      gp_fit_matrix[Vid][1] = 1;
      gp_fit_matrix[Vid][2] = EOSvolumes[Vid-1];
      gp_fit_matrix[Vid][3] = pow(EOSvolumes[Vid-1],2);
      gp_fit_matrix[Vid][4] = pow(EOSvolumes[Vid-1],3);
    }

    isInitialized = true;
  }

  /// Performs a QHA calculation.
  /// 
  /// For a regular QHA calculation, there are two possible choices:
  /// 1) calculate the Grueneisen parameter using the finite difference method.
  ///    This calculation requires 3 phonon calculations (for volumes V, V-dV and V+dV);
  /// 2) calculate temperature-dependent parameters (such as equilibrium volume,
  /// free energy, bulk modulus, thermal expansion, isochoric and isobaric specific heat,
  /// average Grueneisen parameter) employing a set of phonon calculations and
  /// making a fit to some model equation of state.
  /// It requires N_EOSvolumes static and N_EOSvolumes phonon calculations
  /// (N_EOSvolumes is determined by the EOS_DISTORTION_RANGE option).
  /// 
  /// QHA3P and SCQHA calculations require 3 phonon calculations and N_EOSvolumes static
  /// calculations.
  /// 
  /// QHANP calculation requires 2*TaylorExpansionOrder+1 phonon calculations and 
  /// N_EOSvolumes static calculations.
  ///
  void QHA::run(_xflags &xflags, _aflags &aflags, _kflags &kflags)
  {
    bool LDEBUG = (FALSE || DEBUG_QHA || XHOST.DEBUG);

    // ME20220427 - must be one or the other
    if (!isEOS && !isGP_FD) {
      string message = "Nothing to be done. Selector EOS=ON or GP_FINITE_DIFF=ON";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      return;
    }
    string msg = "Performing a QHA calculation.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return;
    }

    try{
      bool eos_static_data_available = false;
      bool eos_apl_data_available = false;
      bool gp_data_available = false;

      // we need to clean THERMO file from previous calculations, since calculated data from
      // each QHA/EOS model is appended to the THERMO file and we do not want to mix
      // results of calculations with potentially different parameters
      if (eos_static_data_available){
        aurostd::RemoveFile(currentDirectory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_THERMO_FILE);
        aurostd::RemoveFile(currentDirectory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_COEFF_FILE);
        aurostd::RemoveFile(currentDirectory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_IMAG_FILE);
      }

      // In a QHA calculation, the EOS flag performs APL calculations for a set of volumes.
      // This flag is used when one is interested in T-dependent properties.
      if (isEOS && isQHA){
        eos_apl_data_available = runAPL(xflags, aflags, kflags, QHA_EOS);

        if (eos_apl_data_available && eos_static_data_available){
          if (includeElectronicContribution && doSommerfeldExpansion) calcDOSatEf();
          if (LDEBUG) writeFrequencies();
          writeFVT(currentDirectory);

          if (qha_options.flag("EOS_MODEL:SJ")){
            writeQHAresults(currentDirectory);
            writeThermalProperties(EOS_SJ, QHA_CALC, currentDirectory);
            writeTphononDispersions(EOS_SJ, QHA_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:BM2")){
            writeThermalProperties(EOS_BIRCH_MURNAGHAN2, QHA_CALC, currentDirectory);
            writeTphononDispersions(EOS_BIRCH_MURNAGHAN2, QHA_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:BM3")){
            writeThermalProperties(EOS_BIRCH_MURNAGHAN3, QHA_CALC, currentDirectory);
            writeTphononDispersions(EOS_BIRCH_MURNAGHAN3, QHA_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:BM4")){
            writeThermalProperties(EOS_BIRCH_MURNAGHAN4, QHA_CALC, currentDirectory);
            writeTphononDispersions(EOS_BIRCH_MURNAGHAN4, QHA_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:M")){
            writeThermalProperties(EOS_MURNAGHAN, QHA_CALC, currentDirectory);
            writeTphononDispersions(EOS_MURNAGHAN, QHA_CALC, currentDirectory);
          }
        }
      }

      // ME20220427 - Do static last. There is so much phonon output that
      // the NOTICE won't be easy to spot otherwise.
      // QHA3P, QHANP and SCQHA require a set of static EOS calculations.
      // But for a QHA calculation, the EOS flag is used to toggle these types of
      // calculations.
      if ((isQHA && isEOS) || isQHA3P || isSCQHA || isQHANP){
        msg = "Checking if all required files from static DFT calculations exist.";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_MESSAGE_);
        vector<vector<bool> > file_is_present(subdirectories_static.size(),
            vector<bool>(4));

        uint n_static_calcs = checkStaticCalculations(file_is_present);
        // read data from a set of static DFT calculations only when all required
        // output files are present.
        // Post-processing that requires a set of static DFT calculations (QHA+EOS,
        // QHA3P, QHANP and SCQHA) will happen only when all required data from this
        // set was read successfully.
        if (n_static_calcs == subdirectories_static.size()){
          eos_static_data_available = readStaticCalculationsData();
        }
        else{
          // if there exists data for at least one completed static DFT calculation,
          // print an error listing missing files/directories.
          // Thus, errors are not printed if it is a first QHA run or QHA was run by
          // mistake after the first run before the static DFT calculations were executed
          if (n_static_calcs > 0) {
            printMissingStaticFiles(file_is_present, subdirectories_static);
          } else {
            // ME20220427
            // Add NOTICE when no static calculations have been performed, or QHA will
            // finish as DONE without telling the user why no output has been produed.
            string message = "Waiting for required STATIC calculations.";
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
          }

          createSubdirectoriesStaticRun(xflags, aflags, kflags, file_is_present);
        }
      }

      bool qhanp_data_available = false;
      if (isQHANP){
        qhanp_data_available = runAPL(xflags, aflags, kflags, QHA_TE);
        if (qhanp_data_available){
          if (qha_options.flag("EOS_MODEL:SJ")){
            writeThermalProperties(EOS_SJ, QHANP_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:BM2")){
            writeThermalProperties(EOS_BIRCH_MURNAGHAN2, QHANP_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:BM3")){
            writeThermalProperties(EOS_BIRCH_MURNAGHAN3, QHANP_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:BM4")){
            writeThermalProperties(EOS_BIRCH_MURNAGHAN4, QHANP_CALC, currentDirectory);
          }

          if (qha_options.flag("EOS_MODEL:M")){
            writeThermalProperties(EOS_MURNAGHAN, QHANP_CALC, currentDirectory);
          }
        }
      }

      // QHA3P and SCQHA require that mode-dependent Grueneisen parameters are
      // calculated via a finite difference numerical derivative.
      // For a QHA calculation, use the GP_FD flag to turn on this type of calculation.
      if (isGP_FD || isQHA3P || isSCQHA){
        gp_data_available = runAPL(xflags, aflags, kflags, QHA_FD);

        if (gp_data_available){
          if (isQHA || isQHA3P){
            double V0 = origStructure.GetVolume()/NatomsOrigCell;
            writeGPpath(V0, currentDirectory);
            writeGPmeshFD(currentDirectory);
            writeAverageGPfiniteDifferences(currentDirectory);

            if (isQHA3P && eos_static_data_available){
              if (qha_options.flag("EOS_MODEL:SJ")){
                writeThermalProperties(EOS_SJ, QHA3P_CALC, currentDirectory);
                writeTphononDispersions(EOS_SJ, QHA3P_CALC, currentDirectory);
              }

              if (qha_options.flag("EOS_MODEL:BM2")){
                writeThermalProperties(EOS_BIRCH_MURNAGHAN2, QHA3P_CALC,
                    currentDirectory);
                writeTphononDispersions(EOS_BIRCH_MURNAGHAN2, QHA3P_CALC,
                    currentDirectory);
              }

              if (qha_options.flag("EOS_MODEL:BM3")){
                writeThermalProperties(EOS_BIRCH_MURNAGHAN3, QHA3P_CALC,
                    currentDirectory);
                writeTphononDispersions(EOS_BIRCH_MURNAGHAN3, QHA3P_CALC,
                    currentDirectory);
              }

              if (qha_options.flag("EOS_MODEL:BM4")){
                writeThermalProperties(EOS_BIRCH_MURNAGHAN4, QHA3P_CALC,
                    currentDirectory);
                writeTphononDispersions(EOS_BIRCH_MURNAGHAN4, QHA3P_CALC,
                    currentDirectory);
              }

              if (qha_options.flag("EOS_MODEL:M")){
                writeThermalProperties(EOS_MURNAGHAN, QHA3P_CALC, currentDirectory);
                writeTphononDispersions(EOS_MURNAGHAN, QHA3P_CALC, currentDirectory);
              }
            }
          }

          if (isSCQHA && eos_static_data_available){
            if (qha_options.flag("EOS_MODEL:SJ")){
              runSCQHA(EOS_SJ, true, currentDirectory);
              writeTphononDispersions(EOS_SJ, SCQHA_CALC, currentDirectory);
            }
            if (qha_options.flag("EOS_MODEL:BM2")){
              runSCQHA(EOS_BIRCH_MURNAGHAN2, true, currentDirectory);
              writeTphononDispersions(EOS_BIRCH_MURNAGHAN2, SCQHA_CALC,
                  currentDirectory);
            }
            if (qha_options.flag("EOS_MODEL:BM3")){
              runSCQHA(EOS_BIRCH_MURNAGHAN3, true, currentDirectory);
              writeTphononDispersions(EOS_BIRCH_MURNAGHAN3, SCQHA_CALC,
                  currentDirectory);
            }
            if (qha_options.flag("EOS_MODEL:BM4")){
              runSCQHA(EOS_BIRCH_MURNAGHAN4, true, currentDirectory);
              writeTphononDispersions(EOS_BIRCH_MURNAGHAN4, SCQHA_CALC,
                  currentDirectory);
            }
            if (qha_options.flag("EOS_MODEL:M")){
              runSCQHA(EOS_MURNAGHAN, true, currentDirectory);
              writeTphononDispersions(EOS_MURNAGHAN, SCQHA_CALC,
                  currentDirectory);
            }
          }
        }
      }
    }
    catch(aurostd::xerror &e){
      pflow::logger(e.whereFileName(),e.whereFunction(), e.what(), currentDirectory,
          *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      return;
    }
  }

  /// Creates subdirectories with aflow.in for a set of static DFT calculations.
  ///
  void QHA::createSubdirectoriesStaticRun(const _xflags &xflags, const _aflags &aflags,
      const _kflags &kflags, const vector<vector<bool> > &file_is_present)
  {
    string msg = "";
    // use static_bands calculations to get a reasonable electronic DOS
    if (qha_options.flag("RELAX_IONS_CELL")){ 
      xinput.xvasp.AVASP_flag_RUN_STATIC_BANDS       = false;
      xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS = true;
    }
    else {
      xinput.xvasp.AVASP_flag_RUN_STATIC_BANDS       = true;
      xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS = false;
    }
    xinput.xvasp.AVASP_flag_RUN_STATIC             = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC       = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX              = false;
    xinput.xvasp.AVASP_flag_GENERATE               = false;
    xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::MODULE");

    stringstream aflow;
    string incar_explicit = xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP.str();

    // create aflow.in if there are no outputs with results in the directory:
    // OUTCAR, EIGENVAL and IBZKPT. The last two are only required when electronic
    // contribution is included
    // The existing aflow.in file will not be overwritten.
    for (uint i=0; i<subdirectories_static.size(); i++){
      if (includeElectronicContribution){
        if (file_is_present[i][ST_DF_EIGENVAL] && file_is_present[i][ST_DF_IBZKPT] &&
            file_is_present[i][ST_DF_OUTCAR]) continue;
      }
      else if (file_is_present[i][ST_DF_OUTCAR]) continue;

      if (aurostd::FileExist(currentDirectory+'/'+subdirectories_static[i]+'/'
            +_AFLOWIN_)) continue;

      msg = "Generate aflow.in file in " + subdirectories_static[i] + " directory.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_MESSAGE_);

      xinput.xvasp.str = origStructure;
      xinput.xvasp.str.InflateVolume(coefEOSVolumes[i]);
      xinput.xvasp.str.is_vasp5_poscar_format = true;

      xinput.setDirectory(currentDirectory);
      xinput.xvasp.AVASP_arun_mode = QHA_ARUN_MODE;
      xinput.xvasp.AVASP_arun_runname = arun_runnames_static[i];

      AVASP_populateXVASP(aflags,kflags,xflags.vflags, xinput.xvasp);
      if (xinput.xvasp.AVASP_path_BANDS.empty()){
        if (xflags.vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string.empty())
          xinput.xvasp.AVASP_path_BANDS = AFLOWRC_DEFAULT_BANDS_LATTICE;
        else
          xinput.xvasp.AVASP_path_BANDS = xflags.vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string;
      }
      if (!xinput.xvasp.AVASP_value_BANDS_GRID)
        xinput.xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;

      // set CONVERT_UNIT_CELL to PRES
      xinput.xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL",true);
      xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL");
      xinput.xvasp.aopts.push_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL", "PRES");
      if (qha_options.flag("RELAX_IONS_CELL")){
        xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::RELAX_TYPE");
        xinput.xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE","IONS_CELL_SHAPE");
      }

      AVASP_MakeSingleAFLOWIN(xinput.xvasp, aflow, true);

      xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP.str(std::string());
      xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP.clear();
      xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP << incar_explicit;
    }
  }

  /// Checks if all required static calculations exist and returns the number of
  /// finished calculations.
  /// 
  int QHA::checkStaticCalculations(vector<vector<bool> > &file_is_present)
  {
    string msg = "";

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return 0;
    }

    int count = 0;
    string outcarfile = "", eigenvfile = "", ibzkptfile = "", directory = "";
    bool all_files_are_present = true;
    for (uint i=0; i<subdirectories_static.size(); i++){
      all_files_are_present = true;

      directory = currentDirectory + '/' + subdirectories_static[i];
      if (aurostd::IsDirectory(directory))
        file_is_present[i][ST_DF_DIRECTORY] = true;
      else
        all_files_are_present = false;

      outcarfile = directory+"/OUTCAR.static";
      if (aurostd::EFileExist(outcarfile))
        file_is_present[i][ST_DF_OUTCAR] = true;
      else
        all_files_are_present = false;

      if (includeElectronicContribution){
        eigenvfile = directory+"/EIGENVAL.static";
        ibzkptfile = directory+"/IBZKPT.static";

        if (aurostd::EFileExist(eigenvfile))
          file_is_present[i][ST_DF_EIGENVAL] = true;
        else
          all_files_are_present = false;

        if (aurostd::EFileExist(ibzkptfile))
          file_is_present[i][ST_DF_IBZKPT] = true;
        else
          all_files_are_present = false;
      }
      if (all_files_are_present) count++;
    }

    return count;
  }

  /// Prints what data output files from static DFT calculations are missing.
  void QHA::printMissingStaticFiles(const vector<vector<bool> > & list,
      const vector<string> &subdirectories)
  {
    string msg = "", missing_files = "";
    for (uint i=0; i<subdirectories.size(); i++){
      if (!list[i][ST_DF_DIRECTORY]){
        msg += "Directory " + subdirectories[i] + " is missing.\n";
      }
      else{
        missing_files = "";
        if (!list[i][ST_DF_OUTCAR]){
          missing_files = "OUTCAR.static";
        }
        if (includeElectronicContribution){
          if(!list[i][ST_DF_EIGENVAL]){
            if (!missing_files.empty()) missing_files += ",";
            missing_files += "EIGENVAL.static";
          }

          if(!list[i][ST_DF_IBZKPT]){
            if (!missing_files.empty()) missing_files += ",";
            missing_files += "IBZKPT.static";
          }
        }
        if (!missing_files.empty()){
          msg += "File(s)   "+subdirectories[i]+"/"+missing_files+" is (are) missing.\n";
        }
      }
    }

    // No exception is thrown since QHA is handling the problem of missing files
    // in a way that a method that depends on missing data will not be evaluated.
    // Since the user might request a number of methods to be run in one aflow run,
    // we want to allow the ones that do not depend on the missing data to finish
    // successfully. For example, a finite difference calculation of Grueneisen
    // parameters does not depend on any data from any static DFT calculation.
    if (!msg.empty()){
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////

  /// Creates subdirectories with aflow.in for a set of APL calculations.
  ///
  void QHA::createSubdirectoriesAPLRun(const _xflags &xflags, const _aflags &aflags,
      const _kflags &kflags, const vector<vector<bool> > &file_is_present,
      QHAtype qhatype)
  {
    string msg = "";

    xinput.xvasp.aopts.opattachedscheme("AFLOWIN_FLAG::MODULE","APL",true);
    xinput.xvasp.aopts.flag("FLAG::VOLUME_PRESERVED",true);
    if (qha_options.flag("RELAX_IONS_CELL")){
      xinput.xvasp.aplopts.opattachedscheme("AFLOWIN_FLAG::APL_RELAX","ON",true);
    }
    else{
      xinput.xvasp.aplopts.opattachedscheme("AFLOWIN_FLAG::APL_RELAX","OFF",true);
    }
    xinput.xvasp.aplopts.opattachedscheme("AFLOWIN_FLAG::APL_HIBERNATE","ON",true);

    stringstream aflow;
    string incar_explicit = xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP.str();

    vector<string> &subdirectories = (QHA_FD==qhatype) ? subdirectories_apl_gp :
      (QHA_EOS==qhatype) ? subdirectories_apl_eos : subdirectories_apl_qhanp;
    vector<string> &arun_runnames = (QHA_FD==qhatype) ? arun_runnames_apl_gp :
      (QHA_EOS==qhatype) ? arun_runnames_apl_eos : arun_runnames_apl_qhanp;
    vector<double> &coefVolumes = (QHA_FD==qhatype) ? coefGPVolumes :
      (QHA_EOS==qhatype) ? coefEOSVolumes : coefQHANPVolumes;

    // create aflow.in if there are no outputs with results in the directory.
    // The existing aflow.in file will not be overwritten.
    for (uint i=0; i<subdirectories.size(); i++){
      if (file_is_present[i][PH_DF_HARMIFC]) continue;

      if (aurostd::FileExist(currentDirectory+'/'+subdirectories[i]+'/'+_AFLOWIN_)) 
        continue;

      msg = "Generating aflow.in file in " + subdirectories[i] + " directory.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_MESSAGE_);

      xinput.xvasp.str = origStructure;
      xinput.xvasp.str.InflateVolume(coefVolumes[i]);
      xinput.xvasp.str.is_vasp5_poscar_format = true;

      xinput.setDirectory(currentDirectory);
      xinput.xvasp.AVASP_arun_mode = QHA_ARUN_MODE;
      xinput.xvasp.AVASP_arun_runname = arun_runnames[i];

      AVASP_populateXVASP(aflags,kflags,xflags.vflags, xinput.xvasp);
      if (xinput.xvasp.AVASP_path_BANDS.empty()){
        if (xflags.vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string.empty())
          xinput.xvasp.AVASP_path_BANDS = AFLOWRC_DEFAULT_BANDS_LATTICE;
        else
          xinput.xvasp.AVASP_path_BANDS = xflags.vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string;
      }
      if (!xinput.xvasp.AVASP_value_BANDS_GRID)
        xinput.xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;

      // set CONVERT_UNIT_CELL to PRES
      xinput.xvasp.aopts.flag("AFLOWIN_FLAG::CONVERT_UNIT_CELL",true);
      xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL");
      xinput.xvasp.aopts.push_attached("AFLOWIN_FLAG::CONVERT_UNIT_CELL", "PRES");
      if (qha_options.flag("RELAX_IONS_CELL")){
        xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::RELAX_TYPE");
        xinput.xvasp.aopts.push_attached("AFLOWIN_FLAG::RELAX_TYPE","IONS_CELL_SHAPE");
      }

      AVASP_MakeSingleAFLOWIN(xinput.xvasp, aflow, true);

      xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP.str(std::string());
      xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP.clear();
      xinput.xvasp.AVASP_INCAR_EXPLICIT_START_STOP << incar_explicit;
    }
  }

  /// Checks if all required APL calculations exist and returns the number of
  /// finished calculations.
  ///
  int QHA::checkAPLCalculations(vector<vector<bool> > &file_is_present, QHAtype qhatype)
  {
    string msg = "";

    if (!isInitialized){
      msg = "QHA was not initialized properly.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return 0;
    }

    vector<string> &subdirectories = (QHA_FD==qhatype) ? subdirectories_apl_gp :
      (QHA_EOS==qhatype) ? subdirectories_apl_eos : subdirectories_apl_qhanp;

    int count = 0;
    string harmifcfile = "", phposcarfile = "", directory = "";
    bool all_files_are_present = true;
    for (uint i=0; i<subdirectories.size(); i++){
      all_files_are_present = true;

      directory = currentDirectory + '/' + subdirectories[i];
      if (aurostd::IsDirectory(directory))
        file_is_present[i][PH_DF_DIRECTORY] = true;
      else
        all_files_are_present = false;

      harmifcfile=directory+'/'+DEFAULT_APL_FILE_PREFIX+DEFAULT_APL_HARMIFC_FILE;
      if (aurostd::EFileExist(harmifcfile))
        file_is_present[i][PH_DF_HARMIFC] = true;
      else
        all_files_are_present = false;

      phposcarfile=directory+'/'+DEFAULT_APL_PHPOSCAR_FILE;
      if (aurostd::EFileExist(phposcarfile))
        file_is_present[i][PH_DF_PHPOSCAR] = true;
      else
        all_files_are_present = false;

      if (all_files_are_present) count++;
    }

    return count;
  }

  /// Prints what data output files from APL calculations are missing.
  void QHA::printMissingAPLFiles(const vector<vector<bool> > & list, QHAtype qhatype)
  {

    vector<string> &subdirectories = (QHA_FD==qhatype) ? subdirectories_apl_gp :
      (QHA_EOS==qhatype) ? subdirectories_apl_eos : subdirectories_apl_qhanp;

    string msg = "", missing_files = "";
    for (uint i=0; i<subdirectories.size(); i++){
      if (!list[i][PH_DF_DIRECTORY]){
        msg += "Directory " + subdirectories[i] + " is missing.\n";
      }
      else{
        missing_files = "";
        if (!list[i][PH_DF_HARMIFC]){
          missing_files = DEFAULT_APL_FILE_PREFIX+DEFAULT_APL_HARMIFC_FILE;
        }
        if (!list[i][PH_DF_PHPOSCAR]){
          if (!missing_files.empty()) missing_files += ",";
          missing_files +=  DEFAULT_APL_PHPOSCAR_FILE;
        }
        if (!missing_files.empty()){
          msg += "File(s)   "+subdirectories[i]+"/"+missing_files+" is (are) missing.\n";
        }
      }
    }

    // No exception is thrown since QHA is handling the problem of missing files
    // in a way that a method that depends on missing data will not be evaluated.
    // Since the user might request a number of methods to be run in one aflow run,
    // we want to allow the ones that do not depend on the missing data to finish
    // successfully. For example, a finite difference calculation of Grueneisen
    // parameters does not depend on any data from any static DFT calculation.
    if (!msg.empty()){
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
    }
  }

  /// Runs APL-related processing and post-processing logic.
  ///
  /// Reads data from a set of APL calculations only when all required output files are
  /// present.
  /// Post-processing will happen only when all required data from this set was read
  /// successfully.
  bool QHA::runAPL(_xflags &xflags, _aflags &aflags, _kflags &kflags, QHAtype qhatype)
  {
    string msg = "Checking if all required files from ";
    switch(qhatype){
      case (QHA_EOS):
        msg += "EOS";
        break;
      case (QHA_FD):
        msg += "FD";
        break;
      case (QHA_TE):
        msg += "TE";
        break;
    }
    msg+=" APL calculations exist.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
        *p_oss, _LOGGER_MESSAGE_);

    vector<string> &subdirectories = (QHA_FD==qhatype) ? subdirectories_apl_gp :
      (QHA_EOS==qhatype) ? subdirectories_apl_eos : subdirectories_apl_qhanp;

    vector<vector<bool> > file_is_present(subdirectories.size(),vector<bool>(2));

    bool apl_data_available = false;
    uint n_apl_calcs = checkAPLCalculations(file_is_present, qhatype);
    if (n_apl_calcs == subdirectories.size()){
      apl_data_available = readAPLCalculationData(subdirectories, kflags, qhatype);
    }
    else{
      /// if there exists data for at least one completed APL calculation, an error is
      /// printed listing missing files/directories.
      /// Thus, errors are not printed if it is a first QHA run or QHA was run by
      /// mistake after the first run before the APL calculations were executed.
      if (n_apl_calcs > 0) {
        printMissingAPLFiles(file_is_present, qhatype);
      } else {
        // ME20220427
        // Add NOTICE when no phonon calculations have been performed, or QHA will
        // finish as DONE without telling the user why no output has been produed.
        string message = "Waiting for required PHONON calculations.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
      }
      createSubdirectoriesAPLRun(xflags, aflags, kflags, file_is_present, qhatype);
    }
    return apl_data_available;
  }

  /// Gathers and processes data from finished APL calculations.
  /// 
  bool QHA::readAPLCalculationData(const vector<string> &subdirectories, _kflags &kflags,
      QHAtype type)
  {
    string msg = "Reading phonon DOS and dispersion relations.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return false;
    }

    int Nqpoints = 0;
    bool apl_data_read_successfully = true;
    bool has_negative_frequencies   = false;

    string phposcarfile = "", directory = "";
    for (uint i=0; i<subdirectories.size(); i++){
      directory = currentDirectory + '/' + subdirectories[i];
      phposcarfile = directory + '/' + DEFAULT_APL_PHPOSCAR_FILE;

      apl::PhononCalculator phcalc(*p_FileMESSAGE, *p_oss);
      phcalc.initialize_supercell(xstructure(phposcarfile, IOVASP_POSCAR), false);// verbose=false: do not write files related to symmetry-analysis
      phcalc.getSupercell().build(apl_options);  // ME20200518
      phcalc.setDirectory(directory);
      phcalc.setNCPUs(kflags);
      phcalc._system = system_title;

      string hibernation_file = directory+'/'+DEFAULT_APL_FILE_PREFIX + 
        DEFAULT_APL_HARMIFC_FILE;
      if (aurostd::EFileExist(hibernation_file)){
        msg = "Reading hibernation file...";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_MESSAGE_);
        try{
          phcalc.awake();
        } catch (aurostd::xerror& e){
          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, e.buildMessageString(), currentDirectory,
              *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
          msg = "Reading data from the hibernation file failed.";
          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_ERROR_);
          return false;
        }
        msg = "Hibernation file was read successfully.";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_MESSAGE_);
      }

      // calculate all phonon-related data: DOS, frequencies along the q-mesh and
      // phonon dispersions
      vector<xvector<double> > dummy_dos_projections; // do not need for QHA
      vector<int> dos_mesh(3);

      vector<string> tokens;

      //ME20200518 - Changed keywords to correspond to rest of APL
      aurostd::string2tokens(apl_options.getattachedscheme("DOSMESH"), tokens,
          string(" xX"));
      for (uint j=0; j<tokens.size(); j++){
        dos_mesh[j] = aurostd::string2utype<int>(tokens[j]);
      }

      phcalc.initialize_qmesh(dos_mesh);

      apl_options.flag("DOS_PROJECT", false);// do not need to calculated DOS projections for QHA
      apl::DOSCalculator dosc(phcalc, apl_options);
      dosc.calc(aurostd::string2utype<double>(apl_options.getattachedscheme("DOSPOINTS")),
          aurostd::string2utype<double>(apl_options.getattachedscheme("DOSSMEAR")));
      // AS20200824 we do not want to overwrite the PHDOSCAR created by the independent APL
      // calculation.
      // A.S. was not sure if it is necessary to write the file at QHA level at all (maybe for
      // debug purposes).
      // In the case if it would be necessary to output the file, for convenience the
      // following line is left commented instead of being deleted.
      //      dosc.writePHDOSCAR(subdirectories[i]);//AS20200824

      if (dosc.hasImaginaryFrequencies()){
        stringstream msg;
        msg << "Phonon dispersions of the APL calculation in the " << subdirectories[i];
        msg << " directory contain imaginary frequencies." << std::endl;
        if (ignore_imaginary){
          msg << "The imaginary parts of the phonon dispersions and phonon DOS will be ignored.";
          msg << " Check if the results are still meaningful!" << std::endl;
          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_WARNING_);
        }
        else{
          msg << "QHA calculation will be stopped after checking all available APL calculations." << std::endl;
          msg << "Please check the phonon DOS and phonon dispersions." << std::endl;
          msg << "Workaround: adjust the EOS_DISTORTION_RANGE option to exclude ";
          msg << "the problematic calculations or ignore this error by setting "; 
          msg << "IGNORE_IMAGINARY=ON." << std::endl;

          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_ERROR_);
          apl_data_read_successfully = false;
          has_negative_frequencies   = true;
        }
      }

      if (phcalc.getSupercell().projectToPrimitive()){
        // if projection to primitive was successful update origStructure on the last
        // iteration
        if (i==subdirectories.size()-1){
          double Volume = origStructure.GetVolume();
          origStructure = phcalc.getSupercell().getInputStructure();
          origStructure.InflateVolume(Volume/origStructure.GetVolume());
          NatomsOrigCell = origStructure.atoms.size();
        }
      }
      else{
        msg = "Could not map the AFLOW standard primitive cell to the supercell. ";
        msg += "Phonon dispersions will be calculated using the original structure instead.";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_WARNING_);
      }
      Nbranches = phcalc.getNumberOfBranches();

      string USER_DC_INITLATTICE="";
      //ME20200518 - Changed keyword to common APL keywords
      int USER_DC_NPOINTS = aurostd::string2utype<int>(
          apl_options.getattachedscheme("DCPOINTS"));
      apl::PhononDispersionCalculator pdisc(phcalc);
      pdisc.initPathLattice(USER_DC_INITLATTICE, USER_DC_NPOINTS);
      pdisc.calc(apl::THZ | apl::ALLOW_NEGATIVE);
      // AS20200824 we do not want to overwrite the PHEIGENVAL created by the independent APL
      // calculation.
      // A.S. was not sure if it is necessary to write the file at QHA level at all (maybe for
      // debug purposes).
      // In the case if it would be necessary to output the file, for convenience the
      // following line is left commented instead of being deleted.
      //      pdisc.writePHEIGENVAL(subdirectories[i]);//AS20200824

      // save all the data that is necessary for QHA calculations
      if (i==0){// qmesh data is the same for all volumes: need to store only once
        qpWeights = phcalc.getQMesh().getWeights();
        qPoints   = phcalc.getQMesh().getIrredQPointsFPOS();
        Nqpoints  = qPoints.size();
        stringstream qpath_stream;
        qpath_stream <<  pdisc.getPHKPOINTS();
        string filename = currentDirectory+"/";
        filename += DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_KPOINTS_FILE;
        if (!aurostd::stringstream2file(qpath_stream, filename)){
          msg = "Error writing to " + filename + "file.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
        }
      }
      else{
        if (qpWeights.size() != phcalc.getQMesh().getWeights().size()){
          msg = "Inconsistent size of q-points weights for calculations at different volumes.";
          throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INDEX_MISMATCH_);
        }

        if (qPoints.size() != phcalc.getQMesh().getIrredQPointsFPOS().size()){
          msg = "Inconsistent number of irreducible q-points for calculations at different volumes.";
          throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INDEX_MISMATCH_);
        }
      }

      // we have two different sets of data: one for the finite difference Grueneisen
      // parameter calculation and one for the EOS APL calculation
      vector<xvector<double> > freqs = dosc.getFreqs();
      switch (type){
        case(QHA_FD):
          gp_ph_dispersions.push_back(pdisc.createEIGENVAL());

          // allocate memory at the first run (entire chunk of memory is required because
          // data will be reordered)
          if (i==0){
            omegaV_mesh = vector<vector<vector<double> > > (Nqpoints,
                vector<vector<double> >(Nbranches, vector<double>(N_GPvolumes)));
          }

          for (int q=0; q<Nqpoints; q++){
            for (int branch=0; branch<Nbranches; branch++){
              // clean imaginary frequencies in case if user wants to ignore them
              omegaV_mesh[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
            }
          }
          break;
        case (QHA_EOS):
          eos_ph_dispersions.push_back(pdisc.createEIGENVAL());
          eos_vib_thermal_properties.push_back(ThermalPropertiesCalculator(dosc,
                *p_FileMESSAGE));

          // allocate memory at the first run (an entire chunk of memory is required because
          // data will be reordered)
          if (i==0){
            omegaV_mesh_EOS = vector<vector<vector<double> > > (Nqpoints,
                vector<vector<double> >(Nbranches, vector<double>(N_EOSvolumes)));
          }

          for (int q=0; q<Nqpoints; q++){
            for (int branch=0; branch<Nbranches; branch++){
              // clean imaginary frequencies in case if user wants to ignore them
              omegaV_mesh_EOS[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
            }
          }
          break;
        case (QHA_TE):
          // allocate memory at the first run (entire chunk of memory is required because
          // data will be reordered)
          if (i==0){
            omegaV_mesh_QHANP = vector<vector<vector<double> > > (Nqpoints,
                vector<vector<double> >(Nbranches, vector<double>(N_QHANPvolumes)));
          }

          for (int q=0; q<Nqpoints; q++){
            for (int branch=0; branch<Nbranches; branch++){
              // clean imaginary frequencies in case if user wants to ignore them
              omegaV_mesh_QHANP[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
            }
          }
          break;
      }
    }

    // write the flag that tells if a specific type of QHA calculation contains
    // imaginary frequencies
    switch(type){
      case(QHA_EOS):
        msg = (std::string)"[QHA]IMAG=" + (has_negative_frequencies ? "YES" : "NO") + "\n";
        break;
      case(QHA_FD):
        msg = (std::string)"[QHA3P]IMAG=" + (has_negative_frequencies ? "YES" : "NO") + "\n";
        msg = (std::string)"[SCQHA]IMAG=" + (has_negative_frequencies ? "YES" : "NO") + "\n";
        break;
      case(QHA_TE):
        msg = (std::string)"[QHANP]IMAG=" + (has_negative_frequencies ? "YES" : "NO") + "\n";
        break;
    }

    string filename = currentDirectory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_IMAG_FILE;
    if (aurostd::FileExist(filename)){
      if (!aurostd::string2file(msg, filename, "APPEND")){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    }
    else{
      if (!aurostd::string2file(msg, filename)){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    }

    return apl_data_read_successfully;
  }

  /// Reads data from a set of static DFT calculations.
  bool QHA::readStaticCalculationsData()
  {
    string msg = "";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return false;
    }

    bool data_read_success = true;

    xOUTCAR outcar;
    string outcarfile = "", directory = "";
    for (uint i=0; i<subdirectories_static.size(); i++){
      msg = "Reading data from the static DFT calculation in the ";
      msg += subdirectories_static[i] + " directory.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_MESSAGE_);

      directory = currentDirectory + '/' + subdirectories_static[i];
      outcarfile = directory+'/'+"OUTCAR.static";
      if (!outcar.GetPropertiesFile(outcarfile)){
        msg = "Could not read the " + outcarfile + " file";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_ERROR_);
        data_read_success = false;
      }

      Efermi_V.push_back(outcar.Efermi);
      E0_V.push_back(outcar.energy_cell/outcar.natoms);

      Nelectrons = outcar.nelectrons;

      if (includeElectronicContribution){
        static_eigvals.push_back(xEIGENVAL(directory+"/EIGENVAL.static"));
        static_ibzkpts.push_back(xIBZKPT(directory+"/IBZKPT.static"));

        if (!static_eigvals.back().m_initialized){
          msg = "Could not read the " + subdirectories_static[i]+"/EIGENVAL.static file.";
          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_ERROR_);
          data_read_success = false;
        }
      }
    }

    return data_read_success;
  }

  /// Returns a frequency obtained by approximation of frequency-volume dependence
  /// to a polynomial.
  ///
  /// The following polynomial is used:
  /// w = a + b*V  + c*V**2 + d*V**3
  ///
  double QHA::calcFrequencyFit(double V, xvector<double> &xomega)
  {
    // set all weight in fit to 1
    xvector<double> s(xomega.rows); for (int i=s.lrows;i<=s.urows;i++) s[i]=1;
    xvector<double> Vpoly(4); for (int i=1; i<=4; i++) Vpoly[i] = pow(V,i-1);

    aurostd::cematrix lsfit(gp_fit_matrix);
    lsfit.LeastSquare(xomega,s);

    return scalar_product(Vpoly, lsfit.GetFitVector());
  }

  /// Calculates the Grueneisen parameter of an individual vibrational mode for a
  /// given volume.
  /// 
  /// gamma_i = -V/w*dw/dV,
  /// where w is a frequency of a given mode at a given volume
  /// 
  /// Volume dependence of the frequency is approximated by the following polynomial:
  /// w = a + b*V  + c*V**2 + d*V**3
  /// 
  /// @param V is a volume at which the Grueneisen parameter is calculated.
  /// @param xomega is an array of the frequency-volume dependence for a given
  /// phonon branch.
  /// 
  double QHA::calcGrueneisen(double V, xvector<double> &xomega, double &w)
  {
    // set all weight in fit to 1
    xvector<double> s(xomega.rows); for (int i=s.lrows;i<=s.urows;i++) s[i]=1;
    aurostd::cematrix lsfit(gp_fit_matrix);

    double a=0.0, b=0.0, c=0.0, d=0.0; // fit parameters
    double gamma = 0; // mode-dependent Grueneisen parameter
    w = 0; // frequency of a given branch that corresponds to the input volume V

    /** Workaround: for w=0, gamma_i is assumed to be zero, since
      / the actual value depends on the path used to approach w->0. */
    bool freqs_are_nonzero = true;
    for (int i=xomega.lrows; i<=xomega.urows; i++){
      if (!(xomega[i] > AUROSTD_IDENTITY_TOL)){
        freqs_are_nonzero = false;
        break;
      }
    }

    // calculate gamma_i only if all frequencies in a given branch are nonzero
    if (freqs_are_nonzero){
      lsfit.LeastSquare(xomega,s);
      a = lsfit.AVec()[0];
      b = lsfit.AVec()[1];
      c = lsfit.AVec()[2];
      d = lsfit.AVec()[3];

#if DEBUG_QHA_GP_FIT
      double err = 0;
      xvector<double> tmp(4);
      tmp[1] = a; tmp[2] = b; tmp[3] = c; tmp[4] = d;

      // check the fit error
      for (int Vid=0; Vid<N_EOSvolumes; Vid++){
        err = std::abs(a+b*EOSvolumes[Vid]+c*pow(EOSvolumes[Vid],2)+d*pow(EOSvolumes[Vid],3)
            -xomega[Vid+1])/xomega[Vid+1];
        err *= 100.0;
        if (err>=10.0){
          string msg="Relative error of the log(w)=f(V) fit (used to ";
          msg+="determine mode-decomposed Grueneisen parameter) is larger than ";
          msg+="10\% for V="+aurostd::utype2string<double>(EOSvolumes[Vid]);

          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory,
              *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
          cerr << __AFLOW_FUNC__ << " original frequencies:" << std::endl;
          cerr << __AFLOW_FUNC__ << " " << xomega << std::endl;
          cerr << __AFLOW_FUNC__ << " fit to frequencies:" << std::endl;
          cerr << __AFLOW_FUNC__ << " " << gp_fit_matrix * tmp << std::endl;
          break;
        }
      }
#endif

      w = a + b*V + c*pow(V,2) + d*pow(V,3);
      gamma = -V/w * (b + 2*c*V + 3*d*pow(V,2));
    }

    return gamma;
  }

  double QHA::calcGrueneisen(double V, xvector<double> &xomega)
  {
    double w = 0.0;
    return calcGrueneisen(V, xomega, w);
  }

  /// Calculates the Grueneisen parameter of an individual vibrational mode using
  /// the central finite difference method.
  /// 
  /// gamma = -V0/w0*dw/dV |V->V0 ~ -V0/w0*(w(V0+dV)-w(V0-dV))/(2*dV)
  /// 
  double QHA::calcGrueneisenFD(const xvector<double> &xomega)
  {
    string msg = "";
    if (xomega.rows != 3){
      msg = "Wrong size of the xomega array passed to calcGrueneisenFD function.";
      msg += "Expected size: 3. Actual size: ";
      msg += aurostd::utype2string<int>(xomega.rows);
      throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INDEX_BOUNDS_);
    }

    if (xomega[2] > _ZERO_TOL_){
      return -GPvolumes[1]/xomega[2]*(xomega[3]-xomega[1])/
        (GPvolumes[2]-GPvolumes[0]);
    }
    else{
      return 0.0;
    }
  }


  /// Calculates the free energy as a function of volume and temperature.
  /// The volume dependency is obtained via a fit to the model equation of state.
  /// @param qha_method defines what kind of QHA calculation is performed.
  ///
  double QHA::calcFreeEnergyFit(double T, double V, EOSmethod eos_method,
      QHAmethod qha_method, uint contrib)
  {
    string msg = "";
    if (T<_ZERO_TOL_) T = 0.0;

    xvector<double> F = calcFreeEnergy(T, qha_method, contrib);
    xvector<double> p = fitToEOSmodel(F, eos_method);
    return evalEOSmodel(V, p, eos_method);
  }

  /// Calculates the free energy as a function of temperature.
  xvector<double> QHA::calcFreeEnergy(double T, QHAmethod qha_method,
      uint contrib)
  {
    string msg = "";

    xvector<double> F(aurostd::vector2xvector(E0_V));
    if (T<_ZERO_TOL_) T = 0.0;

    // include electronic contribution to the free energy
    if (contrib & F_ELEC){
      if (doSommerfeldExpansion) F += calcElectronicFreeEnergySommerfeld(T);
      else{
        for (int id=0; id<N_EOSvolumes; id++){
          F[id+1]+=calcElectronicFreeEnergy(T, id);
        }
      }
    }

    // include vibrational contribution to the free energy
    if (contrib & F_VIB){
      switch(qha_method){
        // here QHA3P and QHANP share the same code
        case (QHA3P_CALC):
        case (QHANP_CALC):
          for (int i=F.lrows; i<=F.urows; i++){
            F[i] += calcVibFreeEnergyTaylorExpansion(T, i-1, qha_method);
          }
          break;
        case(QHA_CALC):
          for (int i=F.lrows; i<=F.urows; i++){
            F[i] += calcVibFreeEnergy(T, i-1);
          }
          break;
        default:
          msg = "Nonexistent QHA method was passed to " + string(__func__);
          throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg,
              _INPUT_UNKNOWN_);
          break;
      }
    }

    return F;
  }

  /// Calculates the bulk modulus for a given EOS model.
  double QHA::calcBulkModulus(double V, const xvector<double> &parameters, EOSmethod method)
  {
    double B = 0.0;
    switch(method){
      case(EOS_SJ):
        {
          double x = std::pow(V, -1.0/3.0);
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x, parameters, 2);
          B = calcBulkModulus_SJ(x, dEdp);
        }
        break;
      case(EOS_BIRCH_MURNAGHAN2):
      case(EOS_BIRCH_MURNAGHAN3):
      case(EOS_BIRCH_MURNAGHAN4):
        {
          double x = std::pow(V, -2.0/3.0);
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x, parameters, 2);
          B = calcBulkModulus_BM(x, dEdp);
        }
        break;
      case(EOS_MURNAGHAN):
        // B = B0*(V0/V)^Bp0 | parameters = {E, V0, B0, Bp0}
        B = eV2GPa*parameters[3]*std::pow(parameters[2]/V, parameters[4]);
        break;
    }

    return B;
  }

  /// Calculates the pressure derivative of the bulk modulus for a given EOS model.
  double QHA::calcBprime(double V, const xvector<double> &parameters, EOSmethod method)
  {
    double Bp = 0.0;
    switch(method){
      case(EOS_SJ):
        {
          double x = std::pow(V, -1.0/3.0);
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x, parameters, 3);
          Bp = Bprime_SJ(x, dEdp);
        }
        break;
      case(EOS_BIRCH_MURNAGHAN2):
      case(EOS_BIRCH_MURNAGHAN3):
      case(EOS_BIRCH_MURNAGHAN4):
        {
          double x = std::pow(V, -2.0/3.0);
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x, parameters, 3);
          Bp = Bprime_BM(x, dEdp);
        }
        break;
      case(EOS_MURNAGHAN):
        Bp = parameters[4]; // for Murnaghan EOS Bp does not depend on volume
        break;
    }
    return Bp;
  }

  /// For a given E(V) relations calculates a pressure the corresponds to a given
  /// volume V.
  double QHA::calcEOS2Pressure(double V, const xvector<double> &parameters, EOSmethod method)
  {
    double P = 0.0;
    switch(method){
      case(EOS_SJ):
        {
          double x = std::pow(V, -1.0/3.0);
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x, parameters, 1);
          P = -std::pow(x,4)/3.0*dEdp[2];
        }
        break;
      case(EOS_BIRCH_MURNAGHAN2):
      case(EOS_BIRCH_MURNAGHAN3):
      case(EOS_BIRCH_MURNAGHAN4):
        {
          double x = std::pow(V, -2.0/3.0);
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x, parameters, 1);
          P = -2.0/3.0*std::pow(x,5.0/2.0)*dEdp[2];
        }
        break;
      case(EOS_MURNAGHAN):
        P = parameters[3]*(1-std::pow(parameters[2]/V,parameters[4]))/parameters[4];
        break;
    }
    return -P;
  }

  /// Calculates an equilibrium volume for a given EOS model.
  double QHA::calcEquilibriumVolume(const xvector<double> &parameters, EOSmethod method)
  {
    double Veq = 0.0;
    switch(method){
      case(EOS_SJ):
        {
          double xmin = std::pow(max(EOSvolumes), -1.0/3.0);
          double xmax = std::pow(min(EOSvolumes), -1.0/3.0);
          Veq = std::pow(aurostd::polynomialFindExtremum(parameters, xmin, xmax), -3);
        }
        break;
      case(EOS_BIRCH_MURNAGHAN2):
      case(EOS_BIRCH_MURNAGHAN3):
      case(EOS_BIRCH_MURNAGHAN4):
        {
          double xmin = std::pow(max(EOSvolumes), -2.0/3.0);
          double xmax = std::pow(min(EOSvolumes), -2.0/3.0);
          Veq = std::pow(aurostd::polynomialFindExtremum(parameters, xmin, xmax), -3.0/2.0);
        }
        break;
      case(EOS_MURNAGHAN):
        Veq = parameters[2]; // parameters = {E, V, B, Bp}
    break;
    }
    return Veq;
  }


  /// Calculates the internal energy as a function of volume and temperature.
  ///
  double QHA::calcInternalEnergyFit(double T, double V, EOSmethod method)
  {
    static xvector<double> U(N_EOSvolumes);
    static xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);

    if (T<0) return 0;

    for (int i=0; i<N_EOSvolumes; i++){
      U[i+1] = eos_vib_thermal_properties[i].getInternalEnergy(T, apl::eV);
      U[i+1] /= NatomsOrigCell;
    }

    return evalEOSmodel(V, fitToEOSmodel(U, method), method);
  }

  /// Calculates the vibrational contribution to the free energy for a
  /// calculation at a specific volume (given by its id in the EOSvolumes array).
  ///
  double QHA::calcVibFreeEnergy(double T, int id)
  {
    return eos_vib_thermal_properties[id].getVibrationalFreeEnergy(T, apl::eV)
      /NatomsOrigCell;
  }

  /// Fits the (free) energy-volume dependency to one of the following equation of state
  /// models:
  /// stabilized jellium: https://doi.org/10.1103/PhysRevB.63.224115
  /// Murnaghan: Proceedings of the National Academy of Sciences of the United States of 
  /// America, 30 (9): 244–247
  /// Birch-Murnaghan: Physical Review. 71 (11): 809–824
  ///
  xvector<double> QHA::fitToEOSmodel(xvector<double> &V, xvector<double> &E,
      EOSmethod method)
  {
    string msg = "";

    if (!isMinimumWithinBounds(E)){
      msg = "Calculation is stopped since there is no free energy minimum";
      msg += " within a given volume range.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_VALUE_RANGE_);
    }

    xvector<double> fit_params;
    xvector<double> guess(4);
    switch(method){
      case(EOS_SJ):
        {
          xvector<double> x(V.rows);
          for (int i=V.lrows; i<=V.urows; i++) x[i] = std::pow(V[i],-1.0/3.0);
          aurostd::cematrix M(aurostd::Vandermonde_matrix(x, 4));
          M.LeastSquare(E);
          fit_params = M.GetFitVector();

          double x_eq = aurostd::polynomialFindExtremum(fit_params, min(x), max(x));
          EOS_volume_at_equilibrium = std::pow(x_eq, -3);
          // index 1 - value, 2 - 1st derivative, 3 - 2nd derivative and so on
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x_eq, fit_params, 3);
          EOS_energy_at_equilibrium = dEdp[1];
          EOS_bulk_modulus_at_equilibrium = calcBulkModulus_SJ(x_eq, dEdp);
          EOS_Bprime_at_equilibrium  = Bprime_SJ(x_eq, dEdp);
        }
        break;
      case(EOS_BIRCH_MURNAGHAN2):
        {
          xvector<double> x(V.rows);
          for (int i=V.lrows; i<=V.urows; i++) x[i] = std::pow(V[i],-2.0/3.0);
          aurostd::cematrix M(aurostd::Vandermonde_matrix(x, 3));
          M.LeastSquare(E);
          fit_params = M.GetFitVector();

          double x_eq = aurostd::polynomialFindExtremum(fit_params, min(x), max(x));
          EOS_volume_at_equilibrium = std::pow(x_eq, -3.0/2.0);
          // index 1 - value, 2 - 1st derivative, 3 - 2nd derivative and so on
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x_eq, fit_params, 3);
          EOS_energy_at_equilibrium = dEdp[1];
          EOS_bulk_modulus_at_equilibrium = calcBulkModulus_BM(x_eq, dEdp);
          EOS_Bprime_at_equilibrium  = Bprime_BM(x_eq, dEdp);
        }
        break;
      case(EOS_BIRCH_MURNAGHAN3):
        {
          xvector<double> x(V.rows);
          for (int i=V.lrows; i<=V.urows; i++) x[i] = std::pow(V[i],-2.0/3.0);
          aurostd::cematrix M(aurostd::Vandermonde_matrix(x, 4));
          M.LeastSquare(E);
          fit_params = M.GetFitVector();

          double x_eq = aurostd::polynomialFindExtremum(fit_params, min(x), max(x));
          EOS_volume_at_equilibrium = std::pow(x_eq, -3.0/2.0);
          // index 1 - value, 2 - 1st derivative, 3 - 2nd derivative and so on
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x_eq, fit_params, 3);
          EOS_energy_at_equilibrium = dEdp[1];
          EOS_bulk_modulus_at_equilibrium = calcBulkModulus_BM(x_eq, dEdp);
          EOS_Bprime_at_equilibrium  = Bprime_BM(x_eq, dEdp);
        }
        break;
      case(EOS_BIRCH_MURNAGHAN4):
        {
          xvector<double> x(V.rows);
          for (int i=V.lrows; i<=V.urows; i++) x[i] = std::pow(V[i],-2.0/3.0);
          aurostd::cematrix M(aurostd::Vandermonde_matrix(x, 5));
          M.LeastSquare(E);
          fit_params = M.GetFitVector();

          double x_eq = aurostd::polynomialFindExtremum(fit_params, min(x), max(x));
          EOS_volume_at_equilibrium = std::pow(x_eq, -3.0/2.0);
          // index 1 - value, 2 - 1st derivative, 3 - 2nd derivative and so on
          xvector<double> dEdp = aurostd::evalPolynomialDeriv(x_eq, fit_params, 3);
          EOS_energy_at_equilibrium = dEdp[1];
          EOS_bulk_modulus_at_equilibrium = calcBulkModulus_BM(x_eq, dEdp);
          EOS_Bprime_at_equilibrium  = Bprime_BM(x_eq, dEdp);
        }
        break;
      case(EOS_MURNAGHAN):
        {
          guess[1] = min(E);
          guess[2] = (max(V)+min(V))/2;
          guess[3] = V[1]*(E[3]-2*E[2]+E[1])/pow(V[2]-V[1],2); // B from central differences 
          guess[4] = 3.5; // a reasonable initial value for most materials

          aurostd::NonlinearFit nlfit(V,E,guess,Murnaghan);
          nlfit.fitLevenbergMarquardt();

          fit_params = nlfit.p;
          EOS_energy_at_equilibrium = nlfit.p[1];
          EOS_volume_at_equilibrium = nlfit.p[2];
          EOS_bulk_modulus_at_equilibrium   = nlfit.p[3]*eV2GPa;
          EOS_Bprime_at_equilibrium  = nlfit.p[4];
        }
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + string(__func__);
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    return fit_params;
  }
  xvector<double> QHA::fitToEOSmodel(xvector<double> &E, EOSmethod method)
  {
    xvector<double> V = aurostd::vector2xvector(EOSvolumes);
    return fitToEOSmodel(V, E, method);
  }

  /// Returns the (free) energy at a given volume for a chosen EOS model.
  double QHA::evalEOSmodel(double V, const xvector<double> &p, EOSmethod method)
  {
    string msg = "";
    double energy = 0;

    static xvector<double> dydp(4);
    switch(method){
      case(EOS_SJ):
        energy = aurostd::evalPolynomial(std::pow(V, -1.0/3.0), p);
        break;
      case(EOS_BIRCH_MURNAGHAN2):
      case(EOS_BIRCH_MURNAGHAN3):
      case(EOS_BIRCH_MURNAGHAN4):
        energy = aurostd::evalPolynomial(std::pow(V, -2.0/3.0), p);
        break;
      case(EOS_MURNAGHAN):
        energy = Murnaghan(V, p, dydp);
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + string(__func__);
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    return energy;
  }


  /// Calculates the equilibrium volume at a given temperature.
  /// @param qha_method defines what kind of QHA calculation is performed.
  /// 
  double QHA::getEqVolumeT(double T, EOSmethod eos_method, QHAmethod qha_method,
      uint contrib)
  {
    string msg = "";
    xvector<double> F = calcFreeEnergy(T, qha_method, contrib);
    fitToEOSmodel(F, eos_method);

    return EOS_volume_at_equilibrium;
  }


  /// Calculates the average Grueneisen parameter (GP) and the isochoric specific 
  /// heat (CV) at a given temperature using a weighted sum over the k-points mesh.
  /// This function is used when the Grueneisen parameter is calculated using the finite
  /// difference method.
  ///
  void QHA::calcCVandGP(double T, double &CV, double &GP)
  {
    CV = 0; GP = 0;
    if (T<_ZERO_TOL_) return;

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh.size();
    int NQpoints = 0; //  total number of q-points

    double Cvi  = 0;   // mode-dependent specific heat at V=const
    double w    = 0;   // frequency for a given volume
    double expx = 0;   // temperature-dependent exponential factor
    double beta = 1.0/KBOLTZEV/T;

    for (uint q=0; q<NIrQpoints; q++){
      for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = xomega[2]*THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
        if (w>0){
          expx = exp(beta*w);

          Cvi = pow(w,2)*expx/pow(expx-1.0,2) * qpWeights[q];

          GP += calcGrueneisenFD(xomega) * Cvi;
          CV += Cvi;
        }
      }
      NQpoints += qpWeights[q];
    }

    GP /= CV;
    CV /= NQpoints; CV /= Nbranches;
    CV *= 3*pow(beta,2);
  }

  /// Calculates the average Grueneisen parameter (GP) and the isochoric specific 
  /// heat (CV) at a given volume and a given temperature.
  /// This function is used when the Grueneisen parameter is calculated using a fit to the
  /// w(V) relation, that is obtained from a set of EOS calculations.
  ///
  void QHA::calcCVandGPfit(double T, double V, double &CV, double &GP)
  {
    CV = 0; GP = 0;
    if (T<_ZERO_TOL_) return;

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh_EOS.size();
    int NQpoints = 0; // the total number of q-points

    double Cvi  = 0;  // mode-dependent specific heat at V=const
    double w    = 0;  // frequency for a given volume
    double expx = 0;  // temperature-dependent exponential factor
    double gamma = 0;
    double beta = 1.0/KBOLTZEV/T;

    for (uint q=0; q<NIrQpoints; q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh_EOS[q][branch]);

        gamma = calcGrueneisen(V, xomega, w);
        w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
        if (w > AUROSTD_IDENTITY_TOL){
          expx = exp(w*beta);

          Cvi = pow(w,2)*expx/pow(expx-1.0,2) * qpWeights[q];

          GP += gamma * Cvi;
          CV += Cvi;
        }
      }
      NQpoints += qpWeights[q];
    }

    GP /= CV;
    CV /= NQpoints; CV /= NatomsOrigCell;
    CV *= pow(beta,2);
  }

  /// Calculates average Grueneisen parameter at infinite temperature.
  /// It corresponds to an unweighted average over mode-dependent Grueneisen parameters.
  double QHA::calcGPinfFit(double V)
  {
    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh_EOS.size();
    int NQpoints = 0; // the total number of q-points
    double w    = 0;  // frequency for a given volume
    double gamma = 0;
    double GP = 0;

    for (uint q=0; q<NIrQpoints; q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh_EOS[q][branch]);

        gamma = calcGrueneisen(V, xomega, w);
        w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
        if (w > AUROSTD_IDENTITY_TOL) GP += gamma * qpWeights[q];
      }
      NQpoints += qpWeights[q];
    }

    GP /= NQpoints; GP/= Nbranches;
    return GP;
  }

  /// Calculates the volumetric thermal expansion coefficient (beta).
  /// 
  /// Central finite differences are used to calculate beta = 1/V dV/dT.
  /// Values of V(T+dT) and V(T-dT) are taken from the EOS fit to the corresponding
  /// free energies.
  /// 
  /// For T=0 K, the function returns 0. Negative temperatures are ignored.
  /// 
  /// @param T is the temperature at which the calculation is done.
  /// @param eos_method defines which model is used for the EOS fit
  /// @return volumetric thermal expansion coefficient.
  /// 
  double QHA::calcThermalExpansion(double T, EOSmethod eos_method,
      QHAmethod qha_method, uint contrib)
  {
    if (!(T>0)) return 0;

    double dT = DCOEFF*T;
    return 0.5*(getEqVolumeT(T+dT,eos_method,qha_method,contrib)
        -getEqVolumeT(T-dT,eos_method,qha_method,contrib))/
      dT/getEqVolumeT(T,eos_method,qha_method,contrib);
  }

  /// Calculates the isochoric specific heat as a second temperature derivative
  /// of the free energy using the central finite differences method.
  /// 
  /// @param eos_method defines which model is used for the EOS fit
  /// 
  double QHA::calcIsochoricSpecificHeat(double T, double V, EOSmethod eos_method,
      QHAmethod qha_method, uint contrib)
  {
    double dT = DCOEFF*T;
    double CV = 0;
    if (T>0){
      CV = -(calcFreeEnergyFit(T+dT, V, eos_method, qha_method, contrib)
          -2*calcFreeEnergyFit(T, V, eos_method, qha_method, contrib)
          +calcFreeEnergyFit(T-dT ,V, eos_method, qha_method, contrib));
      CV *= T/pow(dT,2);
    }

    return CV;
  }

  /// Calculates the entropy as a temperature derivative of the free energy using
  /// the central finite differences method.
  /// 
  /// @param eos_method defines which model is used for EOS fit
  /// 
  double QHA::calcEntropy(double T, double V, EOSmethod eos_method,
      QHAmethod qha_method, uint contrib)
  {
    double dT = DCOEFF*T;
    return -0.5*(calcFreeEnergyFit(T+dT,V,eos_method,qha_method,contrib)
        -calcFreeEnergyFit(T-dT,V,eos_method,qha_method,contrib))/dT;
  }

  /// Calculates the DOS value at the Fermi level using the linear tetrahedron method. 
  /// For the details check: https://doi.org/10.1103/PhysRevB.49.16223
  ///
  xvector<double> QHA::calcDOSatEf()
  {
    xvector<double> energy_tetrahedron(4);//energy eigenvalues at the corners of the tetrahedron
    double weighted_volume = 0; // weighted volume of tetrahedron
    double DEf = 0.0; // DOS at Fermi level
    double Ef = 0.0;  // Fermi energy
    // Eij = energy_tetrahedron[i] - energy_tetrahedron[j]
    double E21 = 0.0, E31 = 0.0, E41 = 0.0, E32 = 0.0, E42 = 0.0, E43 = 0.0; 

    DOS_Ef =  xvector<double>(N_EOSvolumes);

    for (int id=0; id<N_EOSvolumes; id++){
      DEf = 0;
      for (uint band=0; band<static_eigvals[id].number_bands; band++){
        for (uint i=0; i<static_ibzkpts[id].ntetrahedra; i++){
          weighted_volume = static_ibzkpts[id].wtetrahedra * static_ibzkpts[id].vtetrahedra[i][1];

          for (uint s=0; s<static_eigvals[id].spin+1; s++){
            for (int j=1; j<=4; j++)
              energy_tetrahedron[j] = static_eigvals[id].venergy[static_ibzkpts[id].vtetrahedra[i][j+1]-1][band][s];

            aurostd::sort(energy_tetrahedron);
            Ef = Efermi_V[id];

            if (Ef < energy_tetrahedron[1] && energy_tetrahedron[4] < Ef) DEf += 0;

            if (energy_tetrahedron[1] < Ef && Ef <= energy_tetrahedron[2]){
              E21 = energy_tetrahedron[2] - energy_tetrahedron[1];
              E31 = energy_tetrahedron[3] - energy_tetrahedron[1];
              E41 = energy_tetrahedron[4] - energy_tetrahedron[1];

              DEf += weighted_volume*3*pow(Ef-energy_tetrahedron[1],2)/(E21*E31*E41);
            }
            else if (energy_tetrahedron[2] < Efermi_V[id] && Efermi_V[id] <= energy_tetrahedron[3]){
              E21 = energy_tetrahedron[2] - energy_tetrahedron[1];
              E31 = energy_tetrahedron[3] - energy_tetrahedron[1];
              E41 = energy_tetrahedron[4] - energy_tetrahedron[1];

              E32 = energy_tetrahedron[3] - energy_tetrahedron[2];
              E42 = energy_tetrahedron[4] - energy_tetrahedron[2];

              DEf += weighted_volume/(E31*E41)*(3*E21+6*(Ef-energy_tetrahedron[2])
                  -3*(E31+E42)*pow(Ef-energy_tetrahedron[2],2)/(E32*E42));
            }
            else if (energy_tetrahedron[3] < Ef && Ef <= energy_tetrahedron[4]){
              E41 = energy_tetrahedron[4] - energy_tetrahedron[1];
              E42 = energy_tetrahedron[4] - energy_tetrahedron[2];
              E43 = energy_tetrahedron[4] - energy_tetrahedron[3];

              DEf += 3*weighted_volume*pow(energy_tetrahedron[4]-Ef,2)/(E41*E42*E43);
            }
          }
        }
      }
      DOS_Ef[id+1] = (2-static_eigvals[id].spin)*DEf; // factor 2 if non-magnetic and 1 otherwise
    }
    return DOS_Ef;
  }

  /// Calculates integrated DOS as a function of energy and temperature.
  double QHA::calcIDOS(double e, double T, xEIGENVAL &eig)
  {
    double res = 0.0;
    double weight = 0.0;
    for (uint k=0; k<eig.number_kpoints; k++){
      weight = (2-eig.spin)*eig.vweight[k];// factor 2 if non-magnetic and 1 otherwise
      for (uint b=0; b<eig.number_bands; b++){
        for (uint s=0; s<=eig.spin; s++){
          res += weight*aurostd::FermiDirac(eig.venergy[k][b][s], e, T);
        }
      }
    }
    return res;
  }

  /// Returns the chemical potential at a given temperature for a given volume Vid.
  ///
  /// The chemical potential is found as a solution of the equation IDOS(T)-Nelectrons = 0,
  /// where IDOS(T) = integral over the energy region of DOS(E)*FermiDirac(E, mu, T),
  /// where mu is a chemical potential.
  /// The equation is solved using the bisection method and the bracketing interval is
  /// automatically determined starting with the Fermi energy at the one end of the
  /// interval.
  double QHA::calcChemicalPotential(double T, int Vid)
  {
    // step taken to determine the bracketing interval
    static const double dE = max(0.01, KBOLTZEV*T);

    bool LDEBUG = (FALSE || DEBUG_QHA || XHOST.DEBUG);
    if (LDEBUG) cerr << __AFLOW_FUNC__ << " begin" << std::endl;

    // Fermi energy is used as a starting guess for the value of chemical potential
    double guess = Efermi_V[Vid];
    double f = calcIDOS(guess, T, static_eigvals[Vid]) - Nelectrons;
    if (std::abs(f) < AUROSTD_IDENTITY_TOL) return guess;

    double left_end = 0.0, middle = 0.0, right_end = 0.0;
    double f_at_left_end = 0.0, f_at_middle = 0.0, f_at_right_end = 0.0;

    // The value of the initial guess is used to bracket the root.
    // Since IDOS is a monotonically increasing function another end of the bracketing 
    // interval is found by stepping in the direction of increase, if the IDOS value at
    // the "guess" energy is lower than the number of electrons, or in the direction of
    // decrease in the opposite case until the bracketing interval is found (opposite
    // signs of IDOS-number_of_electrons at bracketing interval ends).
    left_end = right_end = guess;
    if (sign(f) > 0){
      f_at_right_end = f;
      do {
        left_end -= dE;
        f_at_left_end = calcIDOS(left_end, T, static_eigvals[Vid]) - Nelectrons;
      } while (sign(f_at_left_end) == sign(f_at_right_end));
    }
    else{
      f_at_left_end = f;
      do {
        right_end += dE;
        f_at_right_end = calcIDOS(right_end, T, static_eigvals[Vid]) - Nelectrons;
      } while (sign(f_at_left_end) == sign(f_at_right_end));
    }

    middle = (left_end + right_end)/2;
    f_at_middle = calcIDOS(middle, T, static_eigvals[Vid]) - Nelectrons;

    if (LDEBUG){
      cerr << __AFLOW_FUNC__ << " Bracketing interval was determined." << std::endl;
      cerr << __AFLOW_FUNC__ << " left_end= "  << left_end  << "middle= ";
      cerr << middle  << "right_end= "  << right_end  << std::endl;
      cerr << __AFLOW_FUNC__ << " f_left= " << f_at_left_end;
      cerr << "f_middle= " << f_at_middle << "f_right= ";
      cerr << f_at_right_end << std::endl;
    }

    // Iterate until the convergence criterion is reached:
    // f(middle) is sufficiently close to zero or the bracketing interval is sufficiently
    // small. The latter is used to avoid an infinite loop.
    // For T->0 it is likely that the IDOS is not smooth but step-like due to numerical
    // discretization.
    // Meanwhile, do a sanity check that the function has opposite signs at the interval 
    // ends.
    while ((sign(f_at_left_end) != sign(f_at_right_end)) && 
        (std::abs(f_at_middle)>_ZERO_TOL_) && (std::abs(left_end-right_end) > _ZERO_TOL_)){
      if (sign(f_at_left_end) == sign(f_at_middle)){
        std::swap(left_end, middle);
        std::swap(f_at_left_end, f_at_middle);
      }

      if (sign(f_at_right_end) == sign(f_at_middle)){
        std::swap(right_end, middle);
        std::swap(f_at_right_end, f_at_middle);
      }

      middle = 0.5*(left_end + right_end);
      f_at_middle = calcIDOS(middle, T, static_eigvals[Vid]) - Nelectrons;

      if (LDEBUG){
        cerr << __AFLOW_FUNC__ << " left_end= "  << left_end  << "middle= ";
        cerr << middle  << "right_end= "  << right_end  << std::endl;
        cerr << __AFLOW_FUNC__ << " f_left= " << f_at_left_end;
        cerr << "f_middle= " << f_at_middle << "f_right= ";
        cerr << f_at_right_end << std::endl;
      }
    }

    if (LDEBUG) cerr << __AFLOW_FUNC__ << " end" << std::endl;
    return middle;
  }

  /// Calculates the electronic free energy at a given temperature for a given volume Vid.
  ///
  /// The electronic free energy is calculated as a weighted sum over energy eigenvalues
  /// with occupancies as defined by Fermi-Dirac statistics.
  ///
  /// Notice: there might be a problem for non-magnetic systems when magnetic calculation
  /// is turned on. Even though eigenvalues for spin up and spin down should be the same,
  /// A.S. encountered a situation when the small differences in the values of
  /// eigenvalues for spin up and down down lead to wrong electronic free energies for 
  /// spin down while the calculation was correct for spin up.
  double QHA::calcElectronicFreeEnergy(double T, int Vid)
  {
    static const double Tmin = 0.1;
    if (T < Tmin) return 0.0;
    double U = 0.0, S = 0.0;
    xEIGENVAL eig = static_eigvals[Vid];
    double mu  = calcChemicalPotential(T, Vid);
    double mu0 = calcChemicalPotential(Tmin, Vid);
    double E = 0.0;
    double f = 0.0, f0 = 0.0;
    double weight = 0.0;
    for (uint k=0; k<eig.number_kpoints; k++){
      weight = (2-eig.spin)*eig.vweight[k];// factor 2 if non-magnetic and 1 otherwise
      for (uint s=0; s<=eig.spin; s++){
        for (uint b=0; b<eig.number_bands; b++){
          E = eig.venergy[k][b][s];
          f = aurostd::FermiDirac(E, mu, T); f0 = aurostd::FermiDirac(E, mu0, Tmin);
          U += E * weight*(f - f0);

          if (f>0 && f<1) S -= (f*log(f) + (1-f)*log(1-f)) * weight;
        }
      }
    }
    S *= KBOLTZEV;

    return (U - T*S)/NatomsOrigCell;
  }

  /// Calculates the electronic free energy using the Sommerfeld expansion.
  xvector<double> QHA::calcElectronicFreeEnergySommerfeld(double T)
  {
    xvector<double> F_Som = DOS_Ef;
    for (int i=F_Som.lrows; i<=F_Som.urows; i++){
      F_Som[i] *= -pow(M_PI*KBOLTZEV*T,2)/6.0/NatomsOrigCell;
    }
    return F_Som;
  }

  /// Calculates the electronic specific heat at a fixed volume using the
  /// Sommerfeld expansion. Units are eV/atom.
  xvector<double> QHA::calcElectronicSpecificHeatSommerfeld(double T)
  {
    xvector<double> Cv_Som = DOS_Ef;
    for (int i=Cv_Som.lrows; i<=Cv_Som.urows; i++){
      Cv_Som[i] *= pow(M_PI*KBOLTZEV,2)*T/3.0/NatomsOrigCell;
    }
    return Cv_Som;
  }

  // QHA3P-related functions

  /// Returns the frequency extrapolated by a Taylor expansion.
  ///
  double QHA::extrapolateFrequency(double V, const xvector<double> &xomega,
      QHAmethod qha_method)
  {
    string msg = "";
    double result = 0.0;
    switch(qha_method){
      // here QHA3P and SCQHA share the same code
      case(QHA3P_CALC):
      case(SCQHA_CALC):
        {
          double dwdV = (xomega[1]-xomega[3])/(GPvolumes[0]-GPvolumes[2]);
          double d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
            pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);

          double dV = (V-GPvolumes[1]);
          result = xomega[2] + dwdV*dV + 0.5*d2wdV2*pow(dV,2);
        }
        break;
      case(QHANP_CALC):
        {
          // volume-derivative of frequency is calculated using central finite 
          // differences:
          // d^n w/dV^2 = Sum_i=0^n (-1)^i binomial(n,i) w(V0 + (n/2-i)*dV)/dV^n
          //
          // Frequency itself is extrapolated using Taylor expansion
          double V0 = QHANPvolumes[TaylorExpansionOrder];
          int order_begin = 0, order_end = 0;
          double deriv = 0.0;
          // iterate to calculate derivatives up to a given order
          for (int order=1; order<=TaylorExpansionOrder; order++){
            deriv = 0.0; // derivative of a given order
            // example for TE of 3rd order: we need the following data:
            // data:            |0|1|2|3|4|5|6|
            // value itself     | | | |*| | | |
            // 1st derivative   | | |*| |*| | |
            // 2nd derivative   | |*| |*| |*| |
            // 3rd derivative   |*| |*| |*| |*|
            // example: for 2nd derivative order_begin = 1 and order_end = 5
            order_begin = TaylorExpansionOrder - order;
            order_end = order_begin + 2*order;
            for (int i=order_begin, j=order; i<=order_end; i+=2, j--){
              deriv += std::pow(-1,j) * aurostd::combinations(order, j) * xomega[i+1];
            }
            deriv /= pow(2*gp_distortion*V0, order);
            result += deriv * pow(V-V0,order)/aurostd::factorial(order);
          }
          result += xomega[TaylorExpansionOrder+1];
        }
        break;
      case (QHA_CALC):
        msg = "QHA::extrapolateFrequency() function is designed to be used only with QHA3P, QHANP or SCQHA method.";
        msg += " However, QHA method was requested.";
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg,
            _INPUT_ILLEGAL_);
        break;
    }
    return result;
  }

  /// Returns the Grueneisen parameter extrapolated by a Taylor expansion.
  ///
  double QHA::extrapolateGrueneisen(double V, const xvector<double> &xomega,
      QHAmethod qha_method)
  {
    double gamma = 0.0;
    double w  = extrapolateFrequency(V, xomega, qha_method);

    if (w>1e-6){
      switch (qha_method){
        // here QHA3P and SCQHA share the same code
        case (QHA3P_CALC):
        case (SCQHA_CALC):
          {
            double dwdV = 0, d2wdV2 = 0;
            double V0 = GPvolumes[1];
            dwdV = (xomega[1]-xomega[3])/(GPvolumes[0]-GPvolumes[2]);
            d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
              pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);

            gamma = -V/w*(dwdV + d2wdV2 * (V-V0));
          }
          break;
        case (QHANP_CALC):
          {
            // volume-derivative of frequency is calculated using central finite 
            // differences:
            // d^n w/dV^2 = Sum_i=0^n (-1)^i binomial(n,i) w(V0 + (n/2-i)*dV)/dV^n
            //
            // Frequency itself is extrapolated using Taylor expansion
            double V0 = QHANPvolumes[TaylorExpansionOrder];
            int order_begin = 0, order_end = 0;
            double deriv = 0.0;
            // iterate to calculate derivatives up to a given order
            for (int order=1; order<=TaylorExpansionOrder; order++){
              deriv = 0.0; // derivative of a given order
              order_begin = TaylorExpansionOrder - order;
              order_end = order_begin + 2*order;
              for (int i=order_begin, j=order; i<=order_end; i+=2, j--){
                deriv += std::pow(-1,j) * aurostd::combinations(order, j) * xomega[i+1];
              }
              deriv /= pow(2*gp_distortion*V0, order);
              // gamma = -V/w dw/dV: reduce power in (V-V0)^order expression
              gamma += deriv * pow(V-V0,order-1)/aurostd::factorial(order-1);
            }
            gamma *= -V/w;
          }
          break;
        default:
          string message = "QHA::extrapolateGrueneisen() function is designed to be used only with QHA3P, QHANP or SCQHA method.";
          message += " However, QHA method was requested.";
          throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, message,
              _INPUT_ILLEGAL_);
          break;
      }
    }

    return gamma;
  }

  /// Calculates the vibrational contribution to the free energy obtained
  /// using a Taylor expansion of the frequencies.
  /// This function is used in the QHA3P and SCQHA methods.
  /// 
  double QHA::calcVibFreeEnergyTaylorExpansion(double T, int Vid, QHAmethod qha_method)
  {
    double w = 0.0; // extrapolated frequency at V_id volume
    double F = 0.0; // Free energy
    double fi = 0.0;
    double beta = 1.0/KBOLTZEV/T;
    int NQpoints = 0;
    xvector<double> xomega;

    vector<vector<vector<double> > > &omegaV = qha_method == QHA3P_CALC ? omegaV_mesh : omegaV_mesh_QHANP;

    for (uint q=0; q<omegaV.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV[q][branch]);

        w = extrapolateFrequency(EOSvolumes[Vid], xomega, qha_method) * THz2Hz *
          PLANCKSCONSTANTEV_h;
        fi = 0.5*w;

        if (w> AUROSTD_IDENTITY_TOL && T>_ZERO_TOL_) fi += KBOLTZEV*T*log(1-exp(-w*beta));

        fi *= qpWeights[q];
        F += fi;
      }
      NQpoints += qpWeights[q];
    }

    F /= NQpoints;
    F /= NatomsOrigCell;

    return F;
  }

  /// Calculates the vibrational internal energy obtained using a Taylor expansion of
  /// the frequencies.
  /// This function is used in the QHA3P and SCQHA methods.
  double QHA::calcInternalEnergyTaylorExpansion(double T, double V, QHAmethod qha_method)
  {
    if (T<_ZERO_TOL_) return 0.0;

    double U = 0.0, ui = 0.0,  w = 0.0; 
    double beta = 1.0/KBOLTZEV/T;
    int NQpoints = 0;
    xvector<double> xomega;

    vector<vector<vector<double> > > &omegaV = qha_method == QHA3P_CALC ? omegaV_mesh : omegaV_mesh_QHANP;

    for (uint q=0; q<omegaV.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV[q][branch]);

        w = extrapolateFrequency(V, xomega, qha_method) * THz2Hz * PLANCKSCONSTANTEV_h;
        ui = 0.5*w;

        if (w>_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ && T>_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) ui += w/(exp(w*beta)-1.0);

        ui *= qpWeights[q];
        U += ui;
      }
      NQpoints += qpWeights[q];
    }

    U /= NQpoints;
    U /= NatomsOrigCell;

    return U;
  }

  // Definition of the functions used by SCQHA.
  // Implementation is based on http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  // and https://doi.org/10.1016/j.commatsci.2016.04.012

  /// Calculates the phononic pressure multiplied by volume for a given temperature and
  /// a given volume.
  /// Check http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  /// for more details
  /// 
  double QHA::calcVPgamma(double T, double V)
  {
    double calcVPgamma = 0.0, ui = 0.0,  w = 0.0; 
    double beta = 1.0/KBOLTZEV/T;
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(V, xomega, SCQHA_CALC) * THz2Hz * PLANCKSCONSTANTEV_h;
        ui = 0.5*w;

        if (w>_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_ && T>_AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) ui += w/(exp(w*beta)-1.0);

        ui *= qpWeights[q];
        calcVPgamma += ui * extrapolateGrueneisen(V, xomega, SCQHA_CALC);
      }
      NQpoints += qpWeights[q];
    }

    calcVPgamma /= NQpoints;
    calcVPgamma /= NatomsOrigCell;

    return calcVPgamma;
  }

  /// Calculates the equilibrium volume for a given temperatures using SCQHA self-consistent
  /// loop procedure.
  /// Check for details:
  /// http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  /// and https://doi.org/10.1016/j.commatsci.2016.04.012
  double QHA::calcSCQHAequilibriumVolume(double T, double Vguess, xvector<double> &fit_params, EOSmethod method)
  {
    const static int max_scqha_iteration = 10000;
    const static double Vtol = 1e-5;
    const static double dV = 1e-3;

    double Pe = 0.0, VPg = 0.0; // electronic pressure and volume multiplied by phononic pressure
    double Vnew = 0.0, V = Vguess;
    int iter = 0;
    while (iter++ < max_scqha_iteration){
      Pe   = calcEOS2Pressure(V, fit_params, method);
      VPg  = calcVPgamma(T, V);
      Vnew = -VPg/Pe; // be careful: note that sign in this expression depends on 
      // whether the pressure or the volume derivative of energy is used (dE/dV = -Pe)
      if (std::abs(V - Vnew)/V > Vtol) V += (Vnew - V) * dV; else break;
    }

    if (iter == max_scqha_iteration){
      string msg="Maximum number of iterations in self consistent loop is reached";
      msg += " at T="+aurostd::utype2string<double>(T)+"K.";
      pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_MESSAGE_);
    }

    return V;
  }

  double QHA::calcSCQHAequilibriumVolume(double T, EOSmethod method)
  {
    xvector<double> E = aurostd::vector2xvector<double>(E0_V);
    xvector<double> fit_params = fitToEOSmodel(E, method);
    double Vguess = SCQHA_INITIAL_VOLUME_FACTOR * calcEquilibriumVolume(fit_params, method);

    return calcSCQHAequilibriumVolume(T, Vguess, fit_params, method);
  }

  /// Performs SCQHA calculations.
  /// There are two implementations:
  /// 1) perform a self-consistent loop for the initial nonzero temperature and extrapolate
  /// the volume at the next temperature step using V(T+dT) ~ (1+beta dT)*V.
  /// Expect it to be inaccurate at high temperatures.
  /// 2) perform a self-consistent loop for each temperature.
  ///
  /// 1) corresponds to the original implementation as described by the following papers:
  /// http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  /// and https://doi.org/10.1016/j.commatsci.2016.04.012
  void QHA::runSCQHA(EOSmethod method, bool all_iterations_self_consistent,
      const string &directory)
  {
    string msg = "Running SCQHA ";
    if (all_iterations_self_consistent)
      msg += "with all temperature steps computed self-consistenly.";
    else
      msg += "with temperature steps computed using the V *= (1 + beta*dT) approximation.";

    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    // get the equilibrium volume from the fit to the EOS model fitted to the set of 
    // energies obtained from static DFT calculations
    xvector<double> E = aurostd::vector2xvector<double>(E0_V);
    xvector<double> fit_params = fitToEOSmodel(E, method);

    // self-consistent loop to determine equilibrium volume at T=0K
    // to avoid division by zero in the self-consistent loop,
    // the initial volume is taken to be 10% bigger
    double V0K = calcSCQHAequilibriumVolume(0,SCQHA_INITIAL_VOLUME_FACTOR *
        EOS_volume_at_equilibrium, fit_params, method);

    // the name of the output file depends on the EOS fit method
    stringstream file;
    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_THERMO_FILE;
    string sc = all_iterations_self_consistent ? "SC_" : "";
    string blockname = "[SCQHA_" + sc;
    switch(method){
      case(EOS_SJ):
        blockname += "SJ_THERMO]";
        break;
      case(EOS_BIRCH_MURNAGHAN2):
        blockname += "BM2_THERMO]";
        break;
      case(EOS_BIRCH_MURNAGHAN3):
        blockname += "BM3_THERMO]";
        break;
      case(EOS_BIRCH_MURNAGHAN4):
        blockname += "BM4_THERMO]";
        break;
      case(EOS_MURNAGHAN):
        blockname += "M_THERMO]";
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + string(__func__);
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    file.precision(10);

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << blockname + "SYSTEM=" << system_title << std::endl;
    file << blockname + "START" << std::endl;

    // print header
    file << setw(5)  << "# T[K]"               << setw(SW) << ' ' <<
      setw(TW) << "Veq[ev/atom]"         << setw(SW) << ' ' <<
      setw(TW) << "F(Veq)[eV/atom]"      << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-5/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv[kB/atom]"         << setw(SW) << ' ' <<
      setw(TW) << "Cp[kB/atom]"         << setw(SW) << ' ' <<
      setw(TW) << "gamma"                << setw(SW) << ' ' <<
      std::endl;

    double T = Temperatures[0];

    // to avoid division by zero in the self-consistent loop,
    // the initial volume is taken to be 10% bigger
    double V = calcSCQHAequilibriumVolume(T, SCQHA_INITIAL_VOLUME_FACTOR * 
        EOS_volume_at_equilibrium, fit_params, method);

    double dT = (Temperatures[Ntemperatures-1]-Temperatures[0])/(Ntemperatures-1);

    // a set of variables predefined here
    double Cvi     = 0.0;  // mode-dependent specific heat at V=const
    double w       = 0.0;  // frequency for a given volume
    double w0K     = 0.0;  // frequency for a volume at T=0K
    double expx    = 0.0;  // temperature-dependent exponential factor
    double expx0   = 0.0;  // temperature-dependent exponential factor
    double ui      = 0.0;  // mode-dependent internal energy
    double Belec   = 0.0;  // "electronic" bulk modulus
    double Bphononic = 0.0;// "phononic" bulk modulus
    double B       = 0.0;  // total bulk modulus
    double gamma   = 0.0;  // Grueneisen parameter
    double fi      = 0.0;  // mode-dependent free energy
    double Feq     = 0.0;  // total free energy for equilibrium volume at given T
    double CP      = 0.0;  // isobaric specific heat
    double CV      = 0.0;  // isochoric specific heat
    double GP      = 0.0;  // average Grueneisen parameter
    double beta    = 0.0;  // coefficient of thermal expansion
    double d2wdV2  = 0.0;  // second derivative of frequency w.r.t volume
    double betaT   = 1.0/KBOLTZEV/T;

    uint NIrQpoints = omegaV_mesh.size(); // number of irreducible q-points
    int NQpoints = 0; // the total number of q-points (to be determined in the loop)

    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];
      if (T > _ZERO_TOL_) betaT = 1.0/KBOLTZEV/T;
      else betaT = 0.0;

      // calculate the next equilibrium volume using a self-consistent loop if beta=0 or
      // if the user wants so
      if (all_iterations_self_consistent || !(std::abs(beta)>0)){
        V = calcSCQHAequilibriumVolume(T, V, fit_params, method);
      }
      else{
        V *= (1 + beta*dT);
      }

      // calculate thermodynamic properties
      Cvi = 0.0; CP = 0.0; CV = 0.0; 
      GP = 0.0; gamma = 0.0;
      w = 0.0; w0K = 0.0; expx = 0.0; expx0 = 0.0;
      ui = 0.0;
      Belec = 0.0; Bphononic = 0.0; B = 0.0;
      fi = 0.0; Feq = 0.0;
      NQpoints = 0;

      for (uint q=0; q<NIrQpoints; q++){
        for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
          xvector<double> xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);
          w = extrapolateFrequency(V, xomega, SCQHA_CALC);
          w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]

          w0K = extrapolateFrequency(V0K, xomega, SCQHA_CALC);
          w0K *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]

          ui = 0.5*w;
          fi = 0.5*w;

          Cvi = 0.0;

          if (T > _ZERO_TOL_){
            expx = exp(w*betaT);
            // use volume at T=0K for calculation of isochoric specific heat
            expx0 = exp(w0K*betaT);
            if (w0K > _ZERO_TOL_) Cvi = pow(w0K,2)*expx0/pow(expx0-1.0,2) * qpWeights[q];

            if (w > _ZERO_TOL_){
              ui += w/(expx - 1.0);
              fi += KBOLTZEV*T*log(1-exp(-w*betaT));
            }
          }
          ui *= qpWeights[q];
          fi *= qpWeights[q];

          d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
            pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);
          d2wdV2 *= THz2Hz*PLANCKSCONSTANTEV_h;

          gamma = extrapolateGrueneisen(V, xomega, SCQHA_CALC);

          if (w > _ZERO_TOL_) Bphononic += ui*pow(V,2)/w*d2wdV2;
          Bphononic -= T*Cvi*pow(betaT,2)*KBOLTZEV*pow(gamma, 2);

          GP += extrapolateGrueneisen(V, xomega, SCQHA_CALC) * Cvi;
          CV += Cvi;
          Feq  += fi;
        }
        NQpoints += qpWeights[q];
      }

      Feq /= NQpoints;
      Feq /= NatomsOrigCell;
      Feq += evalEOSmodel(V, fit_params, method);

      Belec  = calcBulkModulus(V, fit_params, method);
      Bphononic *= eV2GPa/(NQpoints*NatomsOrigCell*V);

      B = Belec + Bphononic;

      GP /= CV;
      CV /= NQpoints; CV /= Nbranches;
      CV *= 3*pow(betaT,2); // [kB/atom]

      if (T > _ZERO_TOL_){
        beta = KBOLTZEV*CV*GP/V/(B/eV2GPa); // [K^-1]
        CP = CV + V*T*B*pow(beta,2)/eV2GPa/KBOLTZEV; // [kB/atom]
      }

      file << setw(5)  << T                   << setw(SW) << ' ' <<
        setw(TW) << V                   << setw(SW) << ' ' <<
        setw(TW) << Feq                 << setw(SW) << ' ' <<
        setw(TW) << B                   << setw(SW) << ' ' <<
        setw(TW) << beta * 1e5          << setw(SW) << ' ' << //[10^-5/K]
        setw(TW) << CV                  << setw(SW) << ' ' <<
        setw(TW) << CP                  << setw(SW) << ' ' <<
        setw(TW) << GP                  << setw(SW) << ' ' <<
        std::endl;
    }
    file << blockname + "STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename, "APPEND")){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // output-related functions

  /// Writes temperature-dependent properties to file.
  ///
  /// List of properties: temperature, volume, free energy, bulk modulus,
  /// thermal expansion coefficient, specific heat at const V, specific heat at
  /// const P, average Grueneisen parameter calculated from mode-dependent
  /// Grueneisen parameters and average Grueneisen calculated from the thermal
  /// expansion coefficient.
  ///
  void QHA::writeThermalProperties(EOSmethod eos_method, QHAmethod qha_method,
      const string &directory)
  {
    static int PRECISION = 10;
    string msg = "";

    // type of qha calculation
    string blockname = "[";
    switch(qha_method){
      case (QHANP_CALC):
        blockname += "QHANP_";
        break;
      case (QHA3P_CALC):
        blockname += "QHA3P_";
        break;
      case(QHA_CALC):
        blockname += "QHA_";
        break;
      default:
        msg = "Nonexistent QHA method was passed to " + string(__func__);
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    // the name of the output file depends on the EOS fit method and on the type of
    // QHA calculation
    switch(eos_method){
      case(EOS_SJ):
        blockname += "SJ";
        break;
      case(EOS_BIRCH_MURNAGHAN2):
        blockname += "BM2";
        break;
      case(EOS_BIRCH_MURNAGHAN3):
        blockname += "BM3";
        break;
      case(EOS_BIRCH_MURNAGHAN4):
        blockname += "BM4";
        break;
      case(EOS_MURNAGHAN):
        blockname += "M";
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + string(__func__);
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    string blockname_coeff  = blockname + "_COEFF]";
    string blockname_mesh   = blockname + "_THERMO_MESH]";
    blockname += "_THERMO]";

    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_THERMO_FILE;
    msg = "Writing T-dependent properties to "+filename;
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file_coeff;
    file_coeff.precision(20);// it's better to have as many significant digits as possible
    stringstream thermo_block, thermo_mesh_block;

    // write header
    thermo_block << setw(5)  << "#T[K]"  << setw(SW) << ' ' <<
      setw(TW) << "V[A^3/atom]"          << setw(SW) << ' ' <<
      setw(TW) << "F(V)[eV/atom]"        << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-5/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "Cp(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "gamma"                << setw(SW) << ' ' <<
      setw(TW) << "Bprime";
    // the following properties are calculated only with a regular QHA calculation
    if (qha_method==QHA_CALC){
      thermo_mesh_block  << setw(5)  << "#T[K]"  << setw(SWTM) << ' ' <<
        setw(TWTM) << "gamma"                    << setw(SWTM) << ' ' <<
        setw(TWTM) << "beta[10^-5/K]"            << setw(SWTM) << ' ' <<
        setw(TWTM) << "Cv[kB/atom]"              << setw(SWTM) << ' ' <<
        setw(TWTM) << "Cp[kB/atom]"              << setw(SWTM) << ' ' <<
        setw(TWTM) << "gamma(V0)"                << setw(SWTM) << ' ' <<
        setw(TWTM) << "beta(V0)[10^-5/K]"        << setw(SWTM) << ' ' <<
        setw(TWTM) << "Cv(V0)[kB/atom]"          << setw(SWTM) << ' ' <<
        setw(TWTM) << "Cp(V0)[kB/atom]";
    }
    thermo_block      << std::endl;
    thermo_mesh_block << std::endl;

    xvector<double> F(N_EOSvolumes); // free energy
    xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);
    aurostd::cematrix VM(aurostd::Vandermonde_matrix(xvolumes, 3));

    F = calcFreeEnergy(0, qha_method, F_VIB);

    fitToEOSmodel(F, eos_method);
    double V0K = EOS_volume_at_equilibrium; // equilibrium volume at 0K
    double B0K = EOS_bulk_modulus_at_equilibrium; // bulk modulus at 0K

    double T = 0.0, Veq = 0.0, Feq = 0.0, B = 0.0, Bp = 0.0, beta = 0.0, CV = 0.0, CP = 0.0, GP = 0.0;
    // the following properties are calculated by weighted sum over q-points mesh
    double CV_mesh_V0K = 0.0, GP_mesh_V0K = 0.0, CP_mesh_V0K = 0.0, beta_mesh_V0K = 0.0;
    double CV_mesh_V   = 0.0, GP_mesh_V   = 0.0, CP_mesh_V   = 0.0, beta_mesh_V   = 0.0;

    // if electronic contributions are requested
    double beta_elec_V0K = 0.0, beta_elec_V = 0.0;
    double CV_elec_V0K   = 0.0, CV_elec_V   = 0.0;

    uint f_contrib = includeElectronicContribution ? F_ELEC | F_VIB : F_VIB;

    xvector<double> p; // EOS fit coefficients

    file_coeff << AFLOWIN_SEPARATION_LINE << std::endl;
    file_coeff << blockname_coeff + "VMIN=" << min(xvolumes) << std::endl;
    file_coeff << blockname_coeff + "VMAX=" << max(xvolumes) << std::endl;
    file_coeff << AFLOWIN_SEPARATION_LINE << std::endl;
    file_coeff << blockname_coeff + "SYSTEM=" << system_title << std::endl;
    file_coeff << blockname_coeff + "START" << std::endl;

    try{
      for (int Tid=0; Tid<Ntemperatures; Tid++){
        T = Temperatures[Tid];

        F = calcFreeEnergy(T, qha_method, f_contrib);
        // stop if energy minimum is no longer within a given set of volumes
        if (!isMinimumWithinBounds(F)){
          msg = "Calculation is stopped at T=" + aurostd::utype2string<double>(T) + " [K]";
          msg+= " since there is no free energy minimum within a given volume range.";
          pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_WARNING_);
          break;
        }

        p = fitToEOSmodel(F, eos_method);

        // output fit coefficients to the file
        file_coeff << setw(5)  << T;
        for (int i=p.lrows; i<=p.urows; i++){
          file_coeff << setw(SW) << ' ' << setw(TW) << p[i];
        }
        file_coeff << std::endl;

        // note that the state of the EOS fit is changed by the calcThermalExpansion and/or
        // calcIsochoricSpecificHeat functions, so save Veq, Feq and B for future use
        Veq = EOS_volume_at_equilibrium;
        Feq = EOS_energy_at_equilibrium;
        B   = EOS_bulk_modulus_at_equilibrium;  // [GPa]
        Bp  = EOS_Bprime_at_equilibrium;
        beta = calcThermalExpansion(T, eos_method, qha_method, f_contrib); // [K^-1]
        CV   = calcIsochoricSpecificHeat(T, Veq, eos_method, qha_method, f_contrib)/KBOLTZEV; // [kB/atom]
        CP   = CV + Veq*T*B*pow(beta,2)/eV2GPa/KBOLTZEV; // [kB/atom]
        GP   = (beta/CV)*B*Veq/eV2GPa/KBOLTZEV;

        // the following properties are calculated only with a regular QHA calculation
        if (qha_method==QHA_CALC){
          if (includeElectronicContribution){
            if (doSommerfeldExpansion){
              xvector<double> Cv_elec = calcElectronicSpecificHeatSommerfeld(T);
              VM.LeastSquare(Cv_elec);
              CV_elec_V0K = aurostd::evalPolynomial(V0K, VM.GetFitVector())/KBOLTZEV; // [kB/atom]
              CV_elec_V   = aurostd::evalPolynomial(Veq, VM.GetFitVector())/KBOLTZEV; // [kB/atom]

              beta_elec_V0K = 2.0*CV_elec_V0K*KBOLTZEV*eV2GPa/(3.0*B0K*V0K);
              beta_elec_V   = 2.0*CV_elec_V *KBOLTZEV*eV2GPa/(3.0*B*Veq);
            }
            else{
              CV_elec_V0K = calcIsochoricSpecificHeat(T, V0K, eos_method,
                  qha_method, F_ELEC)/KBOLTZEV; // [kB/atom]
              CV_elec_V   = calcIsochoricSpecificHeat(T, Veq, eos_method,
                  qha_method, F_ELEC)/KBOLTZEV; // [kB/atom]
              beta_elec_V   = calcThermalExpansion(T, eos_method, qha_method, F_ELEC);
              beta_elec_V0K = beta_elec_V;
            }
          }

          calcCVandGPfit(T, V0K, CV_mesh_V0K, GP_mesh_V0K);
          CV_mesh_V0K += CV_elec_V0K;

          beta_mesh_V0K = KBOLTZEV*CV_mesh_V0K*GP_mesh_V0K/V0K/(B0K/eV2GPa); // [K^-1]
          beta_mesh_V0K += beta_elec_V0K;

          CP_mesh_V0K = CV_mesh_V0K + V0K*T*B0K*pow(beta_mesh_V0K,2)/eV2GPa/KBOLTZEV;//[kB/atom]

          calcCVandGPfit(T, Veq, CV_mesh_V,  GP_mesh_V);
          CV_mesh_V += CV_elec_V;

          beta_mesh_V = KBOLTZEV*CV_mesh_V*GP_mesh_V/Veq/(B/eV2GPa); // [K^-1]
          beta_mesh_V += beta_elec_V;

          CP_mesh_V   = CV_mesh_V + Veq*T*B*pow(beta_mesh_V,2)/eV2GPa/KBOLTZEV; //[kB/atom]
        }

        // write values to file
        thermo_block.unsetf(ios_base::floatfield);
        thermo_block << setw(5)  << T     << setw(SW) << ' ' <<
          std::fixed << std::setprecision(PRECISION) <<
          setw(TW) << Veq                 << setw(SW) << ' ' <<
          setw(TW) << Feq                 << setw(SW) << ' ' << //[eV/atom]
          setw(TW) << B                   << setw(SW) << ' ' <<
          setw(TW) << beta * 1e5          << setw(SW) << ' ' << //[10^-5/K]
          setw(TW) << CV                  << setw(SW) << ' ' << //[kB/atom]
          setw(TW) << CP                  << setw(SW) << ' ' << //[kB/atom]
          setw(TW) << GP                  << setw(SW) << ' ' <<
          setw(TW) << Bp;
        // the following properties are calculated only with a regular QHA calculation
        if (qha_method==QHA_CALC){
          thermo_mesh_block.unsetf(ios_base::floatfield);
          thermo_mesh_block << setw(5)  << T     << setw(SWTM) << ' ' <<
            std::fixed << std::setprecision(PRECISION) <<
            setw(TWTM) << GP_mesh_V              << setw(SWTM) << ' ' <<
            setw(TWTM) << beta_mesh_V * 1e5      << setw(SWTM) << ' ' << //[10^-5/K]
            setw(TWTM) << CV_mesh_V              << setw(SWTM) << ' ' << //[kB/atom]
            setw(TWTM) << CP_mesh_V              << setw(SWTM) << ' ' << //[kB/atom]
            setw(TWTM) << GP_mesh_V0K            << setw(SWTM) << ' ' <<
            setw(TWTM) << beta_mesh_V0K * 1e5    << setw(SWTM) << ' ' << //[10^-5/K]
            setw(TWTM) << CV_mesh_V0K            << setw(SWTM) << ' ' << //[kB/atom]
            setw(TWTM) << CP_mesh_V0K            << setw(SWTM); //[kB/atom]
        }
        thermo_block      << std::endl;
        thermo_mesh_block << std::endl;
      }
    } catch (aurostd::xerror e){
      // QHA throws _VALUE_RANGE_ exception only when there is no minimum in
      // the energy-volume relation: at this point the calculation of 
      // thermodynamic properties should be stopped and a warning should be
      // printed, and all calculated data should be saved to the file
      if (e.whatCode() == _VALUE_RANGE_){
        pflow::logger(e.whereFileName(), e.whereFunction(), e.buildMessageString(),
            currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
      else{
        throw;
      }
    }

    // output the data to the file
    stringstream file;
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << blockname + "SYSTEM=" << system_title << std::endl;
    file << blockname + "START" << std::endl;
    file << thermo_block.rdbuf();
    file << blockname + "STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (qha_method==QHA_CALC){
      file << AFLOWIN_SEPARATION_LINE << std::endl;
      file << blockname_mesh + "SYSTEM=" << system_title << std::endl;
      file << blockname_mesh + "START" << std::endl;
      file << thermo_mesh_block.rdbuf();
      file << blockname_mesh + "STOP" << std::endl;
      file << AFLOWIN_SEPARATION_LINE << std::endl;
    }

    if (aurostd::FileExist(filename)){
      if (!aurostd::stringstream2file(file, filename, "APPEND")){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    }
    else{
      if (!aurostd::stringstream2file(file, filename)){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    }

    file_coeff << blockname_coeff + "STOP" << std::endl;
    file_coeff << AFLOWIN_SEPARATION_LINE << std::endl;

    // write EOS coefficients to file
    filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_COEFF_FILE;
    if (aurostd::FileExist(filename)){
      if (!aurostd::stringstream2file(file_coeff, filename, "APPEND")){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    }
    else{
      if (!aurostd::stringstream2file(file_coeff, filename)){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    }
  }

  /// Writes the F(V,T) data to aflow.qha.FVT.out file
  ///
  void QHA::writeFVT(const string &directory)
  {
    string msg = "Writing F(V,T) relations to file.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file;
    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_FVT_FILE;
    file.precision(10);

    // output the calculation parameters
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_FVT_PARAMETERS]START" << std::endl;
    file << "N_VOLUMES="      << N_EOSvolumes  << std::endl;
    file << "N_TEMPERATURES=" << Ntemperatures << std::endl;
    file << "[QHA_FVT_PARAMETERS]STOP"  << std::endl;

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_FVT]START" << std::endl;

    // output the FVT data
    xvector<double> Felec(N_EOSvolumes);
    double T = 0.0, Ftot = 0.0;

    for (int Tid = 0; Tid < Ntemperatures; Tid++){
      T = Temperatures[Tid];
      file << "# T = " << T << " K" << std::endl;
      if (includeElectronicContribution && doSommerfeldExpansion){
        Felec = calcElectronicFreeEnergySommerfeld(T);
      }
      for (int Vid = 0; Vid < N_EOSvolumes; Vid++){
        if (includeElectronicContribution && !doSommerfeldExpansion){
          Felec[Vid+1]=calcElectronicFreeEnergy(T,Vid);
        }

        Ftot = E0_V[Vid] + calcVibFreeEnergy(T, Vid) + Felec[Vid+1];
        file << setw(TW) << EOSvolumes[Vid]    << setw(SW) << ' '
          << setw(TW) << Ftot << setw(SW) << ' '
          << setw(TW) << calcVibFreeEnergy(T, Vid) << setw(SW) << ' '
          << setw(TW) << Felec[Vid+1]       << setw(SW) << ' '
          << setw(TW) << E0_V[Vid] <<
          std::endl;
      }
      file << std::endl << std::endl;
    }
    file << "[QHA_FVT]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }

  /// Writes mode-dependent Grueneisen parameter along a given path in k-space.
  ///
  void QHA::writeGPpath(double V, const string &directory)
  {
    string msg = "Calculating Grueneisen parameters along the path in";
    msg += " that was used to for the phonon dispersion.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);
    // we will save bands-projected Grueneisen parameter in xEIGENVAL
    xEIGENVAL GPpath(gp_ph_dispersions.front());
    GPpath.Vol = V;
    GPpath.temperature = 0;

    xstructure struc = origStructure;
    struc.SetVolume(V);
    xvector<double> lattice(3);
    lattice[1] = struc.a * 1E-10;
    lattice[2] = struc.b * 1E-10;
    lattice[3] = struc.c * 1E-10;
    GPpath.lattice = lattice;
    GPpath.carstring = "GRUEN_PATH";

    xvector<double> xomega(N_GPvolumes);
    //venergy.at(kpoint number).at(band number).at(spin number)
    for (uint q=0; q<GPpath.venergy.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        for (int Vid=0; Vid<N_GPvolumes; Vid++){
          xomega[Vid+1] = gp_ph_dispersions[Vid].venergy.at(q).at(branch).at(0);
        }

        //GPpath.venergy.at(q).at(branch).at(0) = calcGrueneisen(V, xomega);
        GPpath.venergy.at(q).at(branch).at(0) = calcGrueneisenFD(xomega);
      }
    }

    stringstream eigenval;
    eigenval << GPpath;

    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_GP_PATH_FILE;
    aurostd::stringstream2file(eigenval, filename);
    if (!aurostd::FileExist(filename)){
      msg = "Cannot open "+filename+" file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }


  /// Writes the average Grueneisen parameter, which is calculated using the finite
  /// difference method and isochoric specific heat.
  /// Output file: "aflow.qha.gp.avg.out"
  ///
  void QHA::writeAverageGPfiniteDifferences(const string &directory)
  {
    string msg = "Writing T-dependence of the average Grueneisen parameter.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file;
    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_GP_AVG_FILE;
    file.precision(10);

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_AVG_GP]START" << std::endl;

    file << setw(5)  << "# T[K]"           << setw(SW) << ' ' <<
      setw(TW) << "gamma"            << setw(SW) << ' ' <<
      setw(TW) << "CV"               << setw(SW) << ' ' <<
      std::endl;

    double T = 0.0, CV =0.0, GP = 0.0;
    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];

      calcCVandGP(T, CV, GP);

      file << setw(5)  << T                   << setw(SW) << ' ' <<
        setw(TW) << GP                  << setw(SW) << ' ' <<
        setw(TW) << CV                  << setw(SW) << ' ' <<
        std::endl;
    }

    file << "[QHA_AVG_GP]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }

  /// Writes the mode-dependent Grueneisen parameters calculated at each q-point in the IBZ.
  ///
  void QHA::writeGPmeshFD(const string &directory)
  {
    string msg = "Writing Grueneisen parameters calculated on the mesh of q-points.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file;
    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_GP_MESH_FILE;
    file.precision(6);

    xvector<double> xomega;
    uint NIrQpoints = omegaV_mesh.size();

    xEIGENVAL xeigen;
    xeigen.number_atoms = NatomsOrigCell;
    xeigen.number_loops = 0;
    xeigen.spin = 1;
    xeigen.Vol = origStructure.GetVolume()/NatomsOrigCell;
    xvector<double> lattice(3);
    lattice[1] = origStructure.a * 1e-10;
    lattice[2] = origStructure.b * 1e-10;
    lattice[3] = origStructure.c * 1e-10;
    xeigen.lattice = lattice;
    xeigen.POTIM = 0;
    xeigen.temperature = 0;
    xeigen.carstring = "GRUEN_MESH";
    xeigen.title = system_title;
    xeigen.number_electrons = 0;
    xeigen.number_kpoints = NIrQpoints;
    xeigen.number_bands = Nbranches;

    xeigen.vweight.resize(NIrQpoints);
    xeigen.vkpoint.resize(NIrQpoints, xvector<double>(3));
    xeigen.venergy.resize(NIrQpoints, deque<deque<double> >(xeigen.number_bands,
          deque<double>(2)));

    for (uint q=0; q<NIrQpoints; q++){
      xeigen.vkpoint[q] = qPoints[q];
      xeigen.vweight[q] = qpWeights[q];
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);
        xeigen.venergy[q][branch][0] = xomega[2];
        xeigen.venergy[q][branch][1] = calcGrueneisenFD(xomega);
      }
    }

    file << xeigen;
    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }

  /// Writes the volume-dependent phonon frequencies obtained from the series of EOS APL
  /// calculations.
  ///
  void QHA::writeFrequencies(const string &directory)
  {
    string msg = "";
    stringstream file;
    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_FREQS_FILE;

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_FREQUENCIES]START" << std::endl;
    for (uint q=0; q<omegaV_mesh_EOS.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        for (int Vid=0; Vid<N_EOSvolumes; Vid++){
          file << EOSvolumes[Vid] << " " << omegaV_mesh_EOS[q][branch][Vid] << std::endl;
        }
        file << std::endl << std::endl;
      }
    }
    file << "[QHA_FREQUENCIES]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
    }
  }

  void QHA::writeTphononDispersions(EOSmethod eos_method, QHAmethod qha_method,
      const string &directory)
  {
    string msg = "";
    double V = 0.0; int T = 0;
    xvector<double> xomega;

    int ndigits = aurostd::getZeroPadding(max(ph_disp_temperatures));
    uint f_contrib = includeElectronicContribution ? F_ELEC | F_VIB : F_VIB;

    for (uint i=0; i<ph_disp_temperatures.size(); i++){
      try{
        T = ph_disp_temperatures[i];
        switch(qha_method){
          case (QHA_CALC):
            xomega = xvector<double>(N_EOSvolumes);
            V = getEqVolumeT(T, eos_method, qha_method, f_contrib);
            break;
          case (QHA3P_CALC):
            xomega = xvector<double>(N_GPvolumes);
            V = getEqVolumeT(T, eos_method, qha_method, f_contrib);
            break;
          case (SCQHA_CALC):
            xomega = xvector<double>(N_GPvolumes);
            V = calcSCQHAequilibriumVolume(T, eos_method);
            break;
          case (QHANP_CALC):
            msg = "T-dependent phonon dispersion calculation is not supported for QHANP method";
            pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
            return;
            break;
        }

        msg = "Writing phonon dispersions corresponding to a ";
        msg += "temperature of " + aurostd::utype2string<double>(T) + " K.";
        pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
            _LOGGER_MESSAGE_);

        // we will save T-dependent phonon bands in xEIGENVAL
        xEIGENVAL eig;
        switch (qha_method){
          case(QHA_CALC):
            eig = eos_ph_dispersions.front();
            break;
            // here QHA3P and SCQHA share the same code
          case(QHA3P_CALC):
          case(SCQHA_CALC):
            eig = gp_ph_dispersions.front();
            break;
          case (QHANP_CALC):
            // not supported, this case is handled in an earlier switch statement
            break;
        }
        eig.Vol = V;
        eig.temperature = T;

        xstructure struc = origStructure;
        struc.InflateVolume(V/struc.GetVolume());
        xvector<double> lattice(3);
        lattice[1] = struc.a * 1E-10;
        lattice[2] = struc.b * 1E-10;
        lattice[3] = struc.c * 1E-10;
        eig.lattice = lattice;
        eig.carstring = "PHON";
        eig.title = system_title;

        //venergy.at(kpoint number).at(band number).at(spin number)
        for (uint q=0; q<eig.venergy.size(); q++){
          for (int branch=0; branch<Nbranches; branch++){
            switch (qha_method){
              case (QHA_CALC):
                for (int Vid=0; Vid<N_EOSvolumes; Vid++){
                  xomega[Vid+1] = eos_ph_dispersions[Vid].venergy[q][branch][0];
                }
                eig.venergy[q][branch][0] = calcFrequencyFit(V, xomega);
                break;
                // here QHA3P and SCQHA share the same code
              case (QHA3P_CALC):
              case (SCQHA_CALC):
                for (int Vid=0; Vid<N_GPvolumes; Vid++){
                  xomega[Vid+1] = gp_ph_dispersions[Vid].venergy[q][branch][0];
                }
                eig.venergy[q][branch][0] = extrapolateFrequency(V, xomega, SCQHA_CALC);
                break;
              case (QHANP_CALC):
                // not supported, this case is handled in an earlier switch statement
                break;
            }
          }
        }

        stringstream eig_stream;
        eig_stream << eig;

        string filename = directory + '/';
        switch (qha_method){
          case (QHA_CALC):
            filename += DEFAULT_QHA_FILE_PREFIX;
            break;
          case (QHA3P_CALC):
            filename += DEFAULT_QHA3P_FILE_PREFIX;
            break;
          case (SCQHA_CALC):
            filename += DEFAULT_SCQHA_FILE_PREFIX;
            break;
          case (QHANP_CALC):
            // not supported, this case is handled in an earlier switch statement
            break;
        }
        filename += DEFAULT_QHA_PDIS_FILE;
        filename += ".T"+aurostd::PaddedNumString(T, ndigits)+"K.out";
        if (!aurostd::stringstream2file(eig_stream, filename)){
          msg = "An error occurred when attempted to write "+filename+" file.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
        }
      } catch (aurostd::xerror e){
        // QHA throws _VALUE_RANGE_ exception only when there is no minimum in
        // the energy-volume relation: at this point the calculation of
        // thermodynamic properties should be stopped and a warning should be
        // printed, and all calculated data should be saved to the file
        if (e.whatCode() == _VALUE_RANGE_){
          pflow::logger(e.whereFileName(), e.whereFunction(),
              "called by " + __AFLOW_FUNC__ + " " + e.what(),
              currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }
        else{
          throw;
        }
      }
    }
  }

  void QHA::writeQHAresults(const string &directory)
  {
    string msg = "Writing properties calculated by QHA into aflow.qha.out file.";
    pflow::logger(QHA_ARUN_MODE, __AFLOW_FUNC__, msg, currentDirectory, *p_FileMESSAGE,
        *p_oss, _LOGGER_MESSAGE_);

    try{
      xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);

      uint f_contrib = includeElectronicContribution ? F_ELEC | F_VIB : F_VIB;
      xvector<double> F = calcFreeEnergy(0, QHA_CALC, f_contrib);

      fitToEOSmodel(F, EOS_SJ);
      double V0K = EOS_volume_at_equilibrium; // equilibrium volume at 0K

      double T = 300.0;
      F = calcFreeEnergy(T, QHA_CALC, f_contrib);
      fitToEOSmodel(F, EOS_SJ);

      double volume = EOS_volume_at_equilibrium;
      double free_energy = EOS_energy_at_equilibrium;
      double bulk_modulus = EOS_bulk_modulus_at_equilibrium;
      double bprime = EOS_Bprime_at_equilibrium;
      double thermal_expansion = calcThermalExpansion(T, EOS_SJ, QHA_CALC, f_contrib);

      double CV = 0.0, grueneisen_300K = 0.0;
      calcCVandGPfit(T, V0K, CV, grueneisen_300K);
      double grueneisen = calcGPinfFit(V0K);

      CV = calcIsochoricSpecificHeat(T, volume, EOS_SJ, QHA_CALC, f_contrib)/KBOLTZEV; // [kB/atom]
      double CP = CV + volume*T*bulk_modulus*pow(thermal_expansion,2)/eV2GPa/KBOLTZEV; // [kB/atom]

      stringstream aflow_qha_out;
      aflow_qha_out << AFLOWIN_SEPARATION_LINE << endl;
      aflow_qha_out << "[QHA_RESULTS]START" << endl;
      aflow_qha_out << "gruneisen_qha = " << grueneisen << endl;
      aflow_qha_out << "gruneisen_qha_300K = " << grueneisen_300K << endl;
      aflow_qha_out << "thermal_expansion_qha_300K = " << thermal_expansion;
      aflow_qha_out << " (1/K)" << endl;
      aflow_qha_out << "modulus_bulk_qha_300K = " << bulk_modulus;
      aflow_qha_out << " (GPa)" << endl;
      aflow_qha_out << "modulus_bulk_derivative_pressure_qha_300K = " << bprime << endl;
      aflow_qha_out << "heat_capacity_Cv_atom_qha_300K = " << CV;
      aflow_qha_out << " (kB/atom)" << endl;
      aflow_qha_out << "heat_capacity_Cv_cell_qha_300K = " << CV * NatomsOrigCell;
      aflow_qha_out << " (kB/cell)" << endl;
      aflow_qha_out << "heat_capacity_Cp_atom_qha_300K = " << CP;
      aflow_qha_out << " (kB/atom)" << endl;
      aflow_qha_out << "heat_capacity_Cp_cell_qha_300K = " << CP * NatomsOrigCell;
      aflow_qha_out << " (kB/cell)" << endl;
      aflow_qha_out << "volume_atom_qha_300K = " << volume;
      aflow_qha_out << " (A^3/atom)" << endl;
      aflow_qha_out << "energy_free_atom_qha_300K = " << free_energy;
      aflow_qha_out << " (eV/atom)" << endl;
      aflow_qha_out << "energy_free_cell_qha_300K = " << free_energy * NatomsOrigCell;
      aflow_qha_out << " (eV/cell)" << endl;
      aflow_qha_out << "[QHA_RESULTS]STOP" << endl;
      aflow_qha_out << AFLOWIN_SEPARATION_LINE << endl;

      string filename = directory + '/' + DEFAULT_QHA_FILE_PREFIX + "out";
      if (!aurostd::stringstream2file(aflow_qha_out, filename)){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,msg,_FILE_ERROR_);
      }
    } catch (aurostd::xerror e){
      // QHA throws _VALUE_RANGE_ exception only when there is no minimum in
      // the energy-volume relation: at this point the calculation of
      // thermodynamic properties should be stopped and a warning should be
      // printed, and all calculated data should be saved to the file
      if (e.whatCode() == _VALUE_RANGE_){
        pflow::logger(e.whereFileName(), e.whereFunction(),
            "called by " + __AFLOW_FUNC__ + " " + e.what(),
            currentDirectory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
      else{
        throw;
      }
    }
  }
}

//=============================================================================
//                         Auxiliary functions
namespace apl{
  /// Checks if there is an aflow.in-type file in directory_LIB with a directive to
  /// run QHA and returns its name in AflowInName parameter.
  /// The following variables/filenames are used for the check: AflowInName, _AFLOWIN_
  /// and aflow_qha.in.
  /// If there are a few suitable files, the topmost in the list is returned.
  bool QHA_Get_AflowInName(string &AflowInName, const string &directory_LIB)
  {
    bool qha_aflowin_is_found = false;
    vector<string> vaflowins;
    string stmp="";
    if (!AflowInName.empty()) vaflowins.push_back(AflowInName);
    if (!_AFLOWIN_.empty()) vaflowins.push_back(_AFLOWIN_);
    vaflowins.push_back(_AFLOWIN_QHA_DEFAULT_);

    string aflowin_name = "", aflowin = "", fullpath_aflowin_name = "";
    for (uint i=0; i<vaflowins.size(); i++){
      aflowin_name = vaflowins[i];
      fullpath_aflowin_name = directory_LIB + '/' + aflowin_name;
      if(aurostd::EFileExist(fullpath_aflowin_name,stmp)&&aurostd::IsCompressed(stmp)){aurostd::UncompressFile(stmp);}  //CO20210204 - fix aflow.in.xz
      if (aurostd::FileExist(fullpath_aflowin_name)){
        aflowin = aurostd::file2string(fullpath_aflowin_name);
        if (aflowin.empty()) continue;
      }

      aflowin = aurostd::RemoveComments(aflowin);
      if (aurostd::substring2bool(aflowin,"[AFLOW_QHA]CALC",TRUE) ||
          aurostd::substring2bool(aflowin,"[VASP_QHA]CALC",TRUE)){
        AflowInName = aflowin_name;
        qha_aflowin_is_found = true;
        break;
      }
    }

    return qha_aflowin_is_found;
  }

  string EOSmethod2label(EOSmethod eos_method)
  {
    string function = "EOSmethod2label():", msg = "";

    switch(eos_method){
      case(EOS_SJ):
        return "SJ";
        break;
      case(EOS_BIRCH_MURNAGHAN2):
        return "BM2";
        break;
      case(EOS_BIRCH_MURNAGHAN3):
        return "BM3";
        break;
      case(EOS_BIRCH_MURNAGHAN4):
        return "BM4";
        break;
      case(EOS_MURNAGHAN):
        return "M";
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + function;
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }
  }

  string QHAmethod2label(QHAmethod qha_method)
  {
    string function = "QHAmethod2label():", msg = "";

    switch(qha_method){
      case(QHA_CALC):
        return "QHA";
        break;
      case(QHA3P_CALC):
        return "QHA3P";
        break;
      case(SCQHA_CALC):
        return "SCQHA";
        break;
      case(QHANP_CALC):
        return "QHANP";
        break;
      default:
        msg = "Nonexistent QHA method was passed to " + function;
        throw aurostd::xerror(__AFLOW_FILE__, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }
  }

  /// Checks if the given calculation contains imaginary frequencies by
  /// reading the corresponding flag from filename.
  /// The flag is IMAG and is set to YES if it contains imaginary frequencies, i.e:
  /// [QHA_method]IMAG=NO
  bool hasImaginary(const string& filename, const string &QHA_method)
  {
    string function = "hasImaginary():", msg = "";

    vector<string> vlines;
    bool has_imaginary = false;
    if (!aurostd::efile2vectorstring(filename, vlines)){
      msg = "File " + filename + " does not exist.";
      throw aurostd::xerror(__AFLOW_FILE__, function, msg, _FILE_NOT_FOUND_);
    }

    vector<string> tokens;
    for (uint i=0; i<vlines.size(); i++){
      if (vlines[i].find("["+QHA_method+"]") != std::string::npos){
        if (vlines[i].find("IMAG") != std::string::npos){
          aurostd::string2tokens(vlines[i],tokens,"=");
          if (tokens.size() != 2){
            msg = "Incorrect number of tokens: should be 2 instead of ";
            msg += tokens.size();
            throw aurostd::xerror(__AFLOW_FILE__, function, msg, _FILE_CORRUPT_);
          }

          has_imaginary = tokens[1].find("YES") != std::string::npos;
          break;
        }
      }
    }
    return has_imaginary;
  }

  // AS20201211 OBSOLETE, but might be useful
  /// Links existing APL calculation at the input volume to a corresponding QHA directory
  /// (ARUN.QHA_PHONONS_1.0000).
  /// It is assumed that linking should be performed before the first QHA run.
  // void linkAPLtoQHA()
  // {
  //    string qha_directory_default = ARUN_DIRECTORY_PREFIX+QHA_ARUN_MODE+"_PHONON_1.0000";
  //    string currentDirectory = aurostd::getPWD();
  //    string qha_directory = currentDirectory + '/' + qha_directory_default;
  //    if (!aurostd::IsDirectory(qha_directory)) aurostd::DirectoryMake(qha_directory);

  //    vector<string> vfiles;
  //    aurostd::DirectoryLS(currentDirectory, vfiles);
  //    string dirname = "";
  //    for (uint i=0; i<vfiles.size(); i++){
  //      if (aurostd::substring2bool(vfiles[i], ARUN_DIRECTORY_PREFIX+"APL")){
  //        dirname = currentDirectory+'/'+vfiles[i];
  //        cout << "Linking "+dirname+" to "+qha_directory << endl;
  //        aurostd::LinkFile(dirname, qha_directory);
  //      }
  //    }
  // }
}
