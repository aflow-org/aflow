// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *              AFlow CORMAC TOHER - Duke University 2013-2021             *
// *                                                                         *
// ***************************************************************************
// aflow_ael_elasticity.h
// functions written by
// 2013-2019: cormac.toher@duke.edu

#ifndef _AFLOW_AEL_ELASTICITY_H_
#define _AFLOW_AEL_ELASTICITY_H_

// ***************************************************************************
// Strings for use in read/write operations
#define _AELSTROPT_ string("[AFLOW_AEL]")
#define _AFSTROPT_ string("[AFLOW]")
#define _AELSTR_DIRNAME_ string("AEL_")
#define _ARAELSTR_DIRNAME_ string("ARUN.AEL_")
#define _AELSTR_MESSAGE_ string("00000  MESSAGE AEL: ")
#define _AELSTR_NOTICE_ string("00000  NOTICE AEL: ")
#define _AELSTR_WARNING_ string("WWWWW  WARNING AEL: ")
#define _AELSTR_ERROR_ string("EEEEE  ERROR AEL: ")

using std::string;
using std::vector;
using std::fixed;

class AELStageBreak : public std::exception
{
  public:
    AELStageBreak() { }
};

// Class to hold input and output data for Elastic constants
class _AEL_data {
  public:
    // constructors/destructors
    _AEL_data();
    ~_AEL_data();
    const _AEL_data& operator=(const _AEL_data &b);    

    // Input title and filename
    string dirpathname;
    string sysname;

    // User input data
    double symtolfrac;

    // Variables to control what is calculated within AEL
    bool calcstrainorigin;
    bool fitstrainorigin;
    bool fitrelaxedstruct;
    bool negstrain;
    bool gaussxm_debug;

    // Variable to control where stress tensor is read from vasprun.xml file
    bool vasprunxmlstress;

    // Variable to control whether PREC and ALGO values are set to ACCURATE and NORMAL, or are left as defaults/read from aflow.in file
    bool precaccalgonorm;

    // Variable to control what aflow run type is used (relax, static or relax_static)
    bool relax_static;
    bool static_only;
    bool relax_only;
    bool aflowin_overwrite;

    // Variable to control if this is a postprocessing run (suppresses creation of new ARUN directories)
    bool postprocess;

    // Variable to let other parts of AFLOW know where or not the material is mechanically stable
    bool mechanically_stable;

    // Variable to record applied and external pressure for AEL calculation
    double applied_pressure;
    double average_external_pressure;

    // Variables to control fitting order for polynomial for stress-strain points
    int normalpolyfitorder;
    int shearpolyfitorder;

    // List of failed ARUN directories to be skipped by AEL fitting procedure
    vector<string> failed_arun_list;
    bool autoskipfailedaruns;
    uint skiparunsmax;

    // Variable to control whether elastic stiffness tensor is symmetrized
    bool symmetrize_elastic_tensor;

    // Variable to switch off VASP symmetrization
    bool vasp_symmetry;

    // Variables for cell mass, volume and density used to calculate speed of sound
    double cellmasskg;
    double cellvolumem3;
    double mass_density;

    // Data calculated within AEL
    double poisson_ratio;
    double bulkmodulus_voigt;
    double bulkmodulus_voigt_half;
    double bulkmodulus_reuss;
    double bulkmodulus_reuss_half;
    double shearmodulus_voigt;
    double shearmodulus_voigt_half;
    double shearmodulus_reuss;
    double shearmodulus_reuss_half;
    double bulkmodulus_vrh;
    double bulkmodulus_vrh_half;
    double shearmodulus_vrh;
    double shearmodulus_vrh_half;
    double elastic_anisotropy;
    // Elastic modulus or Young's modulus, using definition E = 9BG/(3B+G)
    double youngsmodulus_vrh;
    // Speed of sound: see J.-P. Poirer, Introduction to the Physics of the Earth's Interior
    double speed_sound_transverse;
    double speed_sound_longitudinal;
    double speed_sound_average;
    // Pugh's modulus = G/B: see Phil. Mag. S7, 45, 367; PRB 84, 121405 and Intermetallics 19, 1275
    double pughs_modulus_ratio;
    // Debye temperature at zero temperature from speed of sound: see J.-P. Poirer, Introduction to the Physics of the Earth's Interior
    double debye_temperature;
    double natomscell;

    // Elastic constants / compliance tensors
    vector<vector<double> > elastic_tensor;
    vector<vector<double> > elastic_tensor_half;
    vector<vector<double> > compliance_tensor;
    vector<vector<double> > compliance_tensor_half;
    vector<double> elastic_eigenvalues;
    vector<double> youngsmodulus_directional;
    vector<double> shearmodulus_directional;
    vector<double> poissonratio_directional;
    // [OBSOLETE] aurostd::xmatrix<double> elastic_tensor;
    // [OBSOLETE] aurostd::xmatrix<double> compliance_tensor;
    // [OBSOLETE] aurostd::xmatrix<double> elastic_tensor(6, 6);
    // [OBSOLETE] aurostd::xmatrix<double> compliance_tensor(6, 6);
    // [OBSOLETE] xmatrix<double> elastic_tensor(6, 6);
    // [OBSOLETE] xmatrix<double> compliance_tensor(6, 6);

    // Normal strain deformations
    vector<double> normal_deformations;

    // Shear strain deformations
    vector<double> shear_deformations;

    // Normal strain matrix
    vector<vector<xmatrix<double> > > normal_strain;

    // Shear strain matrix
    vector<vector<xmatrix<double> > > shear_strain;

    // Normal strain matrix
    vector<vector<double> > normal_deformations_complete;

    // Shear strain matrix
    vector<vector<double> > shear_deformations_complete;

    // Normal strain deformations: only smallest deformations (for linearity check)
    vector<double> normal_deformations_half;

    // Shear strain deformations: only smallest deformations (for linearity check)
    vector<double> shear_deformations_half;

    // Stress matrices
    vector<vector<xmatrix<double> > > normal_stress;
    vector<vector<xmatrix<double> > > shear_stress;
    // vector<vector<xmatrix<double> > > normal_stress_complete;
    // vector<vector<xmatrix<double> > > shear_stress_complete;
    vector<xmatrix<double> > origin_stress;

    // Save energies, pressures and stresses
    vector<double> energycalculated; 
    vector<double> pressurecalculated;
    vector<xmatrix<double> > stresscalculated;
    vector<int> structurecalculated;
    vector<xmatrix<double> > strain_matrix_list;

  private:                                                //
    void free();                                           // free space
};

// Declaration of functions in AEL

// Namespace for functions used by AEL method, to avoid potential naming conflicts with other parts of AFLOW
namespace AEL_functions {
  // Functions to actually run AEL, either directly or from another part of AFLOW
  uint RunElastic_AEL(_xvasp& xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, _AEL_data& AEL_data, ofstream& FileMESSAGE);
  uint AEL_xvasp_flags_populate(_xvasp& xvasp, string& AflowIn, string& AflowInName, string& FileLockName, const string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE);
  uint AEL_Get_AflowInName(string& AflowInName, const string& directory_LIB, bool& ael_aflowin_found); //CT20200715 Function to find aflow.in filename
  uint Get_ElasticProperties_AEL_postprocess(const string& directory, double& ael_bulk_modulus_voigt, double& ael_bulk_modulus_reuss, double& ael_bulk_modulus_vrh, double& ael_shear_modulus_voigt, double& ael_shear_modulus_reuss, double& ael_shear_modulus_vrh, double& ael_poisson_ratio, vector<vector<double> >& ael_elastic_tensor, vector<vector<double> >& ael_compliance_tensor);
  uint Get_PoissonRatio(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& Poissonratio, bool postprocess, ofstream& FileMESSAGE);
  uint Get_BulkModulus(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& BulkModulus, ofstream& FileMESSAGE);
  uint Get_ShearModulus(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& ShearModulus, ofstream& FileMESSAGE);
  // Functions for generating aflow.in input files for strained structures and extracting stress tensor data calculated with VASP
  uint aelvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, string& dirrunname, _AEL_data& AEL_data, ofstream& FileMESSAGE);
  // [OBSOLETE] uint createAFLOWIN(_xvasp& vaspRun, _xvasp& xvasp, _kflags& kbinFlags, _vflags& vaspFlags, _AEL_data& AEL_data, ofstream& FileMESSAGE);
  uint extractstress(vector<_xvasp>& vaspRuns, _AEL_data& AEL_data, vector<string>& dirrunname, ofstream& FileMESSAGE);
  // Functions for implementing the elastic constants method 
  uint elasticityfit(_AEL_data& AEL_data, ofstream& FileMESSAGE);
  uint elasticitycheck(_AEL_data& AEL_data, _xvasp& xvasp, ofstream& FileMESSAGE);
  uint cij_fit(vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, double& Cij, int& npolycoeffwork, bool gxmdebug, ofstream& FileMESSAGE);
  uint elasticmoduli(_AEL_data& AEL_data, ofstream& FileMESSAGE);
} // end of namespace

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *              AFlow CORMAC TOHER - Duke University 2013-2021             *
// *                                                                         *
// ***************************************************************************

