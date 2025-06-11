// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *              AFlow CORMAC TOHER - Duke University 2013-2021             *
// *                                                                         *
// ***************************************************************************
// aflow_ael_elasticity.h
// functions written by
// 2013-2019: cormac.toher@duke.edu

#ifndef _AFLOW_AEL_ELASTICITY_H_
#define _AFLOW_AEL_ELASTICITY_H_

#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xmatrix.h"

#include "flow/aflow_xclasses.h"

// ***************************************************************************
// Strings for use in read/write operations
#define _AELSTROPT_ std::string("[AFLOW_AEL]")
#define _AFSTROPT_ std::string("[AFLOW]")
#define _AELSTR_DIRNAME_ std::string("AEL_")
#define _ARAELSTR_DIRNAME_ std::string("ARUN.AEL_")
#define _AELSTR_MESSAGE_ std::string("00000  MESSAGE AEL: ")
#define _AELSTR_NOTICE_ std::string("00000  NOTICE AEL: ")
#define _AELSTR_WARNING_ std::string("WWWWW  WARNING AEL: ")
#define _AELSTR_ERROR_ std::string("EEEEE  ERROR AEL: ")

class AELStageBreak : public std::exception {
public:
  AELStageBreak() {}
};

// Class to hold input and output data for Elastic constants
class _AEL_data {
public:
    // constructors/destructors
  _AEL_data();
  ~_AEL_data();
  const _AEL_data& operator=(const _AEL_data& b);

    // Input title and filename
  std::string dirpathname;
  std::string sysname;

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
  std::vector<std::string> failed_arun_list;
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
  std::vector<std::vector<double>> elastic_tensor;
  std::vector<std::vector<double>> elastic_tensor_half;
  std::vector<std::vector<double>> compliance_tensor;
  std::vector<std::vector<double>> compliance_tensor_half;
  std::vector<double> elastic_eigenvalues;
  std::vector<double> youngsmodulus_directional;
  std::vector<double> shearmodulus_directional;
  std::vector<double> poissonratio_directional;

    // Normal strain deformations
  std::vector<double> normal_deformations;

    // Shear strain deformations
  std::vector<double> shear_deformations;

    // Normal strain matrix
  std::vector<std::vector<aurostd::xmatrix<double>>> normal_strain;

    // Shear strain matrix
  std::vector<std::vector<aurostd::xmatrix<double>>> shear_strain;

    // Normal strain matrix
  std::vector<std::vector<double>> normal_deformations_complete;

    // Shear strain matrix
  std::vector<std::vector<double>> shear_deformations_complete;

    // Normal strain deformations: only smallest deformations (for linearity check)
  std::vector<double> normal_deformations_half;

    // Shear strain deformations: only smallest deformations (for linearity check)
  std::vector<double> shear_deformations_half;

    // Stress matrices
  std::vector<std::vector<aurostd::xmatrix<double>>> normal_stress;
  std::vector<std::vector<aurostd::xmatrix<double>>> shear_stress;
    // vector<vector<xmatrix<double> > > normal_stress_complete;
    // vector<vector<xmatrix<double> > > shear_stress_complete;
  std::vector<aurostd::xmatrix<double>> origin_stress;

    // Save energies, pressures and stresses
  std::vector<double> energycalculated;
  std::vector<double> pressurecalculated;
  std::vector<aurostd::xmatrix<double>> stresscalculated;
  std::vector<int> structurecalculated;
  std::vector<aurostd::xmatrix<double>> strain_matrix_list;

private:                                                //
  void free();                                           // free space
};

// Declaration of functions in AEL

// Namespace for functions used by AEL method, to avoid potential naming conflicts with other parts of AFLOW
namespace AEL_functions {
  // Functions to actually run AEL, either directly or from another part of AFLOW
  uint RunElastic_AEL(_xvasp& xvasp, const std::string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, _AEL_data& AEL_data, std::ofstream& FileMESSAGE);
  uint AEL_xvasp_flags_populate(_xvasp& xvasp, std::string& AflowIn, std::string& AflowInName, std::string& FileLockName, const std::string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMESSAGE);
  uint AEL_Get_AflowInName(std::string& AflowInName, const std::string& directory_LIB, bool& ael_aflowin_found); // CT20200715 Function to find aflow.in filename
  uint Get_ElasticProperties_AEL_postprocess(const std::string& directory,
                                             double& ael_bulk_modulus_voigt,
                                             double& ael_bulk_modulus_reuss,
                                             double& ael_bulk_modulus_vrh,
                                             double& ael_shear_modulus_voigt,
                                             double& ael_shear_modulus_reuss,
                                             double& ael_shear_modulus_vrh,
                                             double& ael_poisson_ratio,
                                             std::vector<std::vector<double>>& ael_elastic_tensor,
                                             std::vector<std::vector<double>>& ael_compliance_tensor);
  uint Get_PoissonRatio(_xvasp& xvasp, const std::string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& Poissonratio, bool postprocess, std::ofstream& FileMESSAGE);
  uint Get_BulkModulus(_xvasp& xvasp, const std::string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& BulkModulus, std::ofstream& FileMESSAGE);
  uint Get_ShearModulus(_xvasp& xvasp, const std::string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& ShearModulus, std::ofstream& FileMESSAGE);
  // Functions for generating aflow.in input files for strained structures and extracting stress tensor data calculated with VASP
  uint aelvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, std::string& dirrunname, _AEL_data& AEL_data, std::ofstream& FileMESSAGE);
  uint extractstress(std::vector<_xvasp>& vaspRuns, _AEL_data& AEL_data, std::vector<std::string>& dirrunname, std::ofstream& FileMESSAGE);
  // Functions for implementing the elastic constants method
  uint elasticityfit(_AEL_data& AEL_data, std::ofstream& FileMESSAGE);
  uint elasticitycheck(_AEL_data& AEL_data, _xvasp& xvasp, std::ofstream& FileMESSAGE);
  uint cij_fit(std::vector<double>& xdata_to_fit, std::vector<double>& ydata_to_fit, double& Cij, int& npolycoeffwork, bool gxmdebug, std::ofstream& FileMESSAGE);
  uint elasticmoduli(_AEL_data& AEL_data, std::ofstream& FileMESSAGE);
} // namespace AEL_functions

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *              AFlow CORMAC TOHER - Duke University 2013-2021             *
// *                                                                         *
// ***************************************************************************
