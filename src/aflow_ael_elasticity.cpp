// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023*
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AEL_ELASTICITY_CPP
#define _AFLOW_AEL_ELASTICITY_CPP
#include "aflow.h"
#include "aflow_ael_elasticity.h"
#include "aflow_agl_debye.h"


// [OBSOLETE] #define LIB2RAW_FIX_THIS_LATER FALSE // Not needed anymore due to "--quiet" option in lib2lib

// ###############################################################################
//                  AFLOW Automatic Elasticity Library (AEL) (2014-2021)
// ###############################################################################
//
// Uses strain-stress calculations to obtain elastic constants of materials
// Based on original Python program written by M. de Jong et al.
// See Scientific Data 2, 150009 (2015) for details of original program
// See Phys. Rev. Materials 1, 015401 (2017) for details of this implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AEL
//
// ***************************************************************************
//
// Constructor for _AEL_data class, which is a container for variables required for the AEL method
// Initializes values of variables
//
_AEL_data::_AEL_data() {
  dirpathname = "";
  sysname = "";
  symtolfrac = 0.0;
  fitstrainorigin = false;
  calcstrainorigin = false;
  fitrelaxedstruct = false;
  vasprunxmlstress = false;
  precaccalgonorm = false;
  relax_static = false;
  static_only = false;
  relax_only = false;
  postprocess = false;
  aflowin_overwrite = false;
  mechanically_stable = true; 
  applied_pressure = 0.0;
  average_external_pressure = 0.0;
  negstrain = true;
  gaussxm_debug = false;
  normalpolyfitorder = 4;
  shearpolyfitorder = 4;
  cellmasskg = 0.0;
  cellvolumem3 = 0.0;
  mass_density = 0.0;
  poisson_ratio = 0.0;
  bulkmodulus_voigt = 0.0;
  bulkmodulus_voigt_half = 0.0;
  bulkmodulus_reuss = 0.0;
  bulkmodulus_reuss_half = 0.0;
  shearmodulus_voigt = 0.0;
  shearmodulus_voigt_half = 0.0;
  shearmodulus_reuss = 0.0;
  shearmodulus_reuss_half = 0.0;
  bulkmodulus_vrh = 0.0;
  bulkmodulus_vrh_half = 0.0;
  shearmodulus_vrh = 0.0;
  shearmodulus_vrh_half = 0.0;
  elastic_anisotropy = 0.0;
  youngsmodulus_vrh = 0.0;
  speed_sound_transverse = 0.0;
  speed_sound_longitudinal = 0.0;
  speed_sound_average = 0.0;
  debye_temperature = 0.0;
  natomscell = 0.0;
  pughs_modulus_ratio = 0.0;
  elastic_tensor.clear();
  elastic_tensor_half.clear();
  compliance_tensor.clear();
  compliance_tensor_half.clear();
  elastic_eigenvalues.clear();
  normal_strain.clear();
  shear_strain.clear();
  normal_stress.clear();
  shear_stress.clear();
  origin_stress.clear();
  energycalculated.clear();
  pressurecalculated.clear();
  stresscalculated.clear();
  structurecalculated.clear();
  strain_matrix_list.clear();
  normal_deformations.clear();
  shear_deformations.clear();
  normal_deformations_complete.clear();
  shear_deformations_complete.clear();  
  failed_arun_list.clear();
  autoskipfailedaruns = false;
  skiparunsmax = 1;
  symmetrize_elastic_tensor = false;
  vasp_symmetry = false;
}

//
// Destructor for AEL _AEL_data class
//
_AEL_data::~_AEL_data() {
  free();
}

void _AEL_data::free() {
}

//
// Copy constructor for AEL _AEL_data class
//
//_AGL_data::_AGL_data(const _AGL_data & b)
const _AEL_data& _AEL_data::operator=(const _AEL_data& b) {       // operator=
  if(this != &b) {
    dirpathname = b.dirpathname;
    sysname = b.sysname;
    symtolfrac = b.symtolfrac;
    fitstrainorigin = b.fitstrainorigin;
    calcstrainorigin = b.calcstrainorigin;
    fitrelaxedstruct = b.fitrelaxedstruct;
    vasprunxmlstress = b.vasprunxmlstress;
    precaccalgonorm = b.precaccalgonorm;
    relax_static = b.relax_static;
    static_only = b.static_only;
    relax_only = b.relax_only;
    postprocess = b.postprocess;
    aflowin_overwrite = b.aflowin_overwrite;
    mechanically_stable = b.mechanically_stable;
    applied_pressure = b.applied_pressure;
    average_external_pressure = b.average_external_pressure;    
    negstrain = b.negstrain;
    gaussxm_debug = b.gaussxm_debug;
    normalpolyfitorder = b.normalpolyfitorder;
    shearpolyfitorder = b.shearpolyfitorder;
    cellmasskg = b.cellmasskg;
    cellvolumem3 = b.cellvolumem3;
    mass_density = b.mass_density;
    poisson_ratio = b.poisson_ratio;
    bulkmodulus_voigt = b.bulkmodulus_voigt;
    bulkmodulus_voigt_half = b.bulkmodulus_voigt_half;
    bulkmodulus_reuss = b.bulkmodulus_reuss;
    bulkmodulus_reuss_half = b.bulkmodulus_reuss_half;
    shearmodulus_voigt = b.shearmodulus_voigt;
    shearmodulus_voigt_half = b.shearmodulus_voigt_half;
    shearmodulus_reuss = b.shearmodulus_reuss;
    shearmodulus_reuss_half = b.shearmodulus_reuss_half;
    bulkmodulus_vrh = b.bulkmodulus_vrh;
    bulkmodulus_vrh_half = b.bulkmodulus_vrh_half;
    shearmodulus_vrh = b.shearmodulus_vrh;
    shearmodulus_vrh_half = b.shearmodulus_vrh_half;
    elastic_anisotropy = b.elastic_anisotropy;
    youngsmodulus_vrh = b.youngsmodulus_vrh;
    speed_sound_transverse = b.speed_sound_transverse;
    speed_sound_longitudinal = b.speed_sound_longitudinal;
    speed_sound_average = b.speed_sound_average;
    debye_temperature = b.debye_temperature;
    natomscell = b.natomscell;
    pughs_modulus_ratio = b.pughs_modulus_ratio;
    elastic_tensor = b.elastic_tensor;
    elastic_tensor_half = b.elastic_tensor_half;
    compliance_tensor = b.compliance_tensor;
    compliance_tensor_half = b.compliance_tensor_half;
    elastic_eigenvalues= b.elastic_eigenvalues;
    normal_strain = b.normal_strain;
    shear_strain = b.shear_strain;
    normal_stress = b.normal_stress;
    shear_stress = b.shear_stress;
    origin_stress = b.origin_stress;
    energycalculated = b.energycalculated;
    pressurecalculated = b.pressurecalculated;
    stresscalculated = b.stresscalculated;
    structurecalculated = b.structurecalculated;
    strain_matrix_list = b.strain_matrix_list;
    normal_deformations = b.normal_deformations;
    shear_deformations = b.shear_deformations;
    normal_deformations_complete = b.normal_deformations_complete;
    shear_deformations_complete = b.shear_deformations_complete;  
    failed_arun_list = b.failed_arun_list;
    autoskipfailedaruns = b.autoskipfailedaruns;
    skiparunsmax = b.skiparunsmax;
    symmetrize_elastic_tensor = b.symmetrize_elastic_tensor;
    vasp_symmetry = b.vasp_symmetry;
  }
  return *this;
}
// ***************************************************************************
// KBIN::VASP_RunPhonons_AEL
// ***************************************************************************
namespace KBIN {
  //
  // Run AEL method: uses strain-stress calculations to obtain elastic constants of materials
  // See Phys. Rev. Materials 1, 015401 (2017) and Scientific Data 2, 150009 (2015) for details
  // This function reads in user selections, creates set of strained structures to be run with VASP, reads VASP output and calls gibbsrun function to calculate elastic constants
  //
  void VASP_RunPhonons_AEL(  _xvasp&  xvasp,
      string  AflowIn,
      _aflags& aflags,
      _kflags& kflags,
      _vflags& vflags, ofstream& FileMESSAGE) {
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    uint aelerror = 1;
    bool USER_RELAX=true;
    ostringstream aus;

    if (XHOST.GENERATE_AFLOWIN_ONLY && USER_RELAX) {
      USER_RELAX = false;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "RELAX will be switched OFF for generate_aflowin_only." << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }

    if (USER_RELAX) {
      aelerror = relaxStructureAEL_VASP(AflowIn, xvasp, aflags, kflags, vflags, FileMESSAGE);
      if (aelerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Failed to relax structure prior to running AEL." << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
        return;
      }
    }

    // Call RunElastic_AEL to run AEL
    aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic constants run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic constants run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
  }
} // namespace KBIN

namespace KBIN {
  //CT20200323 - Relax structure if it hasn't relaxed yet before running AEL
  uint relaxStructureAEL_VASP(const string& AflowIn,
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

    // First check if relaxation has already been performed
    for (uint i = 1; i <= num_relax; i++) {
      string filename = aurostd::CleanFileName(xvasp.Directory + "/CONTCAR.relax" + aurostd::utype2string<int>(i)); //CO20200624 - patching --run vs --run=1 bug
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Relaxation CONTCAR filename = " << filename << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if (aurostd::FileExist(filename) || aurostd::FileExist(filename + ".xz")) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Relaxation " << i << " has already completed"  << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Relaxation " << i << " has not completed"  << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        relax_complete = false;
      }
    }
    if (relax_complete) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Relaxation has already completed"  << endl;  
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
// KBIN::VASP_RunPhonons_AEL
// ***************************************************************************
namespace KBIN {
  //
  // Run AEL postprocessing: calls AEL from other parts of AFLOW, to run AEL postprocessing on previously run calculations
  // See Computer Physics Communications 158, 57-72 (2004), Journal of Molecular Structure (Theochem) 368, 245-255 (1996), Phys. Rev. B 90, 174107 (2014) and Phys. Rev. Materials 1, 015401 (2017) for details
  void VASP_RunPhonons_AEL_postprocess(const string& directory_LIB, string& AflowInName, string& FileLockName) {  
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // OBSOLETE uint aelerror = 0;
    ostringstream aus;
    _xvasp xvasp;
    string AflowIn = "";
    _aflags aflags;
    _kflags kflags;
    _vflags vflags;
    ofstream FileMESSAGE;

    // Call AEL_xvasp_flags_populate to populate xvasp, aflags, kflags and vflags classes
    uint aelerror = AEL_functions::AEL_xvasp_flags_populate(xvasp, AflowIn, AflowInName, FileLockName, directory_LIB, aflags, kflags, vflags, FileMESSAGE);
    if (aelerror != 0) { //CT20200624: Check for error, aelerror=2 implies not an AEL main directory
      if (aelerror == 2) {
        cerr << _AELSTR_MESSAGE_ << "Not AEL main directory" << endl;
        return;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "AEL set xvasp flags failed" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return;
      }
    }

    // Set AEL postprocess flag to true
    AEL_data.postprocess = true;

    // Call RunDebye_AEL to run AEL
    aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic constants run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic constants run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
  }
} // namespace KBIN

// ***************************************************************************
// AEL_functions::Get_ElasticProperties_AEL_postprocess
// ***************************************************************************
namespace AEL_functions {
  //
  // Run AEL postprocessing: calls AEL from other parts of AFLOW, to run AEL postprocessing on previously run calculations
  // Returns elastic properties: bulk and shear moduli, Poisson ratio, elastic stiffness tensor
  // See Computer Physics Communications 158, 57-72 (2004), Journal of Molecular Structure (Theochem) 368, 245-255 (1996), Phys. Rev. B 90, 174107 (2014) and Phys. Rev. Materials 1, 015401 (2017) for details
  uint Get_ElasticProperties_AEL_postprocess(const string& directory, double& ael_bulk_modulus_voigt, double& ael_bulk_modulus_reuss, double& ael_bulk_modulus_vrh, double& ael_shear_modulus_voigt, double& ael_shear_modulus_reuss, double& ael_shear_modulus_vrh, double& ael_poisson_ratio, vector<vector<double> >& ael_elastic_tensor, vector<vector<double> >& ael_compliance_tensor) {  
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // [OBSOLETE] uint aelerror = 0;
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
    vector<double> elastic_tensor_row;
    string AflowInName, FileLockName;

    // string filedirname = aurostd::CleanFileName(directory) + "agl.LOCK";
    if (aurostd::FileExist(aurostd::CleanFileName(directory) + "ael.LOCK")) {
      FileLockName = "ael.LOCK";
    } else if (aurostd::FileExist(aurostd::CleanFileName(directory) + "agl.LOCK")) {
      FileLockName = "agl.LOCK";
    } else {
      FileLockName = _AFLOWLOCK_;
    }

    if (aurostd::FileExist(aurostd::CleanFileName(directory) + _AFLOWIN_AEL_DEFAULT_)) {
      AflowInName = _AFLOWIN_AEL_DEFAULT_;
    } else if (aurostd::FileExist(aurostd::CleanFileName(directory) + _AFLOWIN_AGL_DEFAULT_)) {
      AflowInName = _AFLOWIN_AGL_DEFAULT_;
    } else {
      AflowInName = _AFLOWIN_;
    }

    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Directory path = " << directory << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Call AEL_xvasp_flags_populate to populate xvasp, aflags, kflags and vflags classes
    uint aelerror = AEL_functions::AEL_xvasp_flags_populate(xvasp, AflowIn, AflowInName, FileLockName, directory, aflags, kflags, vflags, FileMESSAGE);

    // Set AEL postprocess flag to true
    AEL_data.postprocess = true;

    // Call RunDebye_AEL to run AEL
    aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aelerror;
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic constants run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ael_bulk_modulus_voigt = AEL_data.bulkmodulus_voigt;
      ael_bulk_modulus_reuss = AEL_data.bulkmodulus_reuss;
      ael_bulk_modulus_vrh = AEL_data.bulkmodulus_vrh;
      ael_shear_modulus_voigt = AEL_data.shearmodulus_voigt;
      ael_shear_modulus_reuss = AEL_data.shearmodulus_reuss;
      ael_shear_modulus_vrh = AEL_data.shearmodulus_vrh;
      ael_poisson_ratio = AEL_data.poisson_ratio;
      ael_elastic_tensor = AEL_data.elastic_tensor;
      ael_compliance_tensor = AEL_data.compliance_tensor;
      return aelerror;
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aelerror;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic constants run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aelerror;
    }
  }
} // namespace KBIN

// ***************************************************************************
//
// Set of functions to allow AEL to be called from other parts of AFLOW to obtain Poisson ratio, bulk modulus, etc.
// See Phys. Rev. Materials 1, 015401 (2017) and Scientific Data 2, 150009 (2015) for details
// Note that AEL will need to be run twice if the strain-stress calculations are not already available - see the README file for details. 
// If AEL needs to be rerun after performing the strain-stress calculations, these functions will return a value of 8. If AEL has run and the data is available, they will return 0.
// If there is an error, they will return an unsigned integer value other than 0 or 8.
// For the function Get_ElasticTensor_xmatrix, the data structure passed to ElasticTensor must be a 6x6 xmatrix<double>
//
// ***************************************************************************
// AEL_functions::Get_PoissonRatio
// ***************************************************************************
namespace AEL_functions {
  uint Get_PoissonRatio(_xvasp& xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& Poissonratio, bool postprocess, ofstream& FileMESSAGE) {
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // [OBSOLETE] uint aelerror;
    ostringstream aus;

    // Set AEL postprocess flag
    AEL_data.postprocess = postprocess;

    // Call RunElastic_AEL to run AEL
    uint aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      Poissonratio = AEL_data.poisson_ratio;
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aelerror;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::Get_BulkModulus
// ***************************************************************************
namespace AEL_functions {
  uint Get_BulkModulus(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& BulkModulus, ofstream& FileMESSAGE) {
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // [OBSOLETE] uint aelerror;
    ostringstream aus;

    // Call RunElastic_AEL to run AEL
    uint aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      BulkModulus = AEL_data.bulkmodulus_vrh;
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aelerror;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::Get_ShearModulus
// ***************************************************************************
namespace AEL_functions {
  uint Get_ShearModulus(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, double& ShearModulus, ofstream& FileMESSAGE) {
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // [OBSOLETE] uint aelerror;
    ostringstream aus;

    // Call RunElastic_AEL to run AEL
    uint aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ShearModulus = AEL_data.shearmodulus_vrh;
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aelerror;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::Get_ElasticTensor
// ***************************************************************************
namespace AEL_functions {
  uint Get_ElasticTensor(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, vector<vector<double> > ElasticTensor, ofstream& FileMESSAGE) {
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // uint aelerror;
    ostringstream aus;

    // Call RunElastic_AEL to run AEL
    uint aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ElasticTensor = AEL_data.elastic_tensor;
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aelerror;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::Get_ElasticTensor_xmatrix
// ***************************************************************************
namespace AEL_functions {
  uint Get_ElasticTensor_xmatrix(_xvasp&  xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, xmatrix<double> ElasticTensor, ofstream& FileMESSAGE) {
    // Class to contain AEL input and output data
    _AEL_data AEL_data;
    // [OBSOLETE] uint aelerror;
    ostringstream aus;

    // Call RunElastic_AEL to run AEL
    uint aelerror = AEL_functions::RunElastic_AEL(xvasp, AflowIn, aflags, kflags, vflags, AEL_data, FileMESSAGE);
    if(!AEL_data.mechanically_stable) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Negative stiffness tensor eigenvalues indicate mechanical instability" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }
    // Copy values from AEL_data.elastic_tensor (vector of vectors) to ElasticTensor (xmatrix)
    // Note that the data structure passed in ElasticTensor must be defined as a 6x6 xmatrix<double> 
    if(aelerror == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run completed successfully!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for(uint i = 0; i < AEL_data.elastic_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.elastic_tensor.at(i).size(); j++) {
          ElasticTensor[i+1][j+1] = AEL_data.elastic_tensor.at(i).at(j);
        }
      } 
    } else if(aelerror == 8) { 
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "AEL Elastic run waiting for other calculations!" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "AEL Elastic run failed" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return aelerror;
  }
} // namespace AEL_functions

// ***************************************************************************
// AGL_functions::RunElastic_AEL
// ***************************************************************************
namespace AEL_functions {
  // Function to actually run the AEL method
  uint RunElastic_AEL(_xvasp& xvasp, const string& AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
    // User's control of parameters for elastic constants calculation; setting defaults

    // Options to decide what to calculate and write
    aurostd::xoption USER_WRITE_FULL_RESULTS;             USER_WRITE_FULL_RESULTS.option = false;
    aurostd::xoption USER_NEG_STRAINS;                    USER_NEG_STRAINS.option = true;
    aurostd::xoption USER_STRAIN_SYMMETRY;                USER_STRAIN_SYMMETRY.option = false;
    aurostd::xoption USER_CHECK_ELASTIC_SYMMETRY;         USER_CHECK_ELASTIC_SYMMETRY.option = true;
    aurostd::xoption USER_ORIGIN_STRAIN_FIT;              USER_ORIGIN_STRAIN_FIT.option = false;
    aurostd::xoption USER_ORIGIN_STRAIN_CALC;             USER_ORIGIN_STRAIN_CALC.option = false;
    aurostd::xoption USER_RELAXED_STRUCT_FIT;             USER_RELAXED_STRUCT_FIT.option = false;
    aurostd::xoption USER_VASPRUNXML_STRESS;              USER_VASPRUNXML_STRESS.option = false;
    aurostd::xoption USER_PRECACC_ALGONORM;               USER_PRECACC_ALGONORM.option = false;
    aurostd::xoption USER_RELAX_STATIC;                   USER_RELAX_STATIC.option = false;
    aurostd::xoption USER_STATIC;                         USER_STATIC.option = false;
    aurostd::xoption USER_RELAX;                          USER_RELAX.option = false;
    aurostd::xoption USER_SPECIES_MASS;                   USER_SPECIES_MASS.option = false;
    aurostd::xoption USER_DIRNAME_ARUN;                   USER_DIRNAME_ARUN.option = false;
    aurostd::xoption USER_GAUSSXM_DEBUG;                  USER_GAUSSXM_DEBUG.option = false;
    aurostd::xoption USER_AFLOWIN_OVERWRITE;              USER_AFLOWIN_OVERWRITE.option = false;
    aurostd::xoption USER_AUTOSKIP_FAILED_ARUNS;          USER_AUTOSKIP_FAILED_ARUNS.option = false;
    aurostd::xoption USER_PRESSURE_CALC;                  USER_PRESSURE_CALC.option = false;
    aurostd::xoption USER_USE_STABLE_PRESSURE;            USER_USE_STABLE_PRESSURE.option = false;        
    aurostd::xoption USER_SYMMETRIZE;                     USER_SYMMETRIZE.option = false;
    aurostd::xoption USER_VASPSYM;                        USER_VASPSYM.option = false;        

    // Parameters for determining calculation details
    string USER_SYSTEM_NAME                 = "AEL";
    int    USER_NNORMAL_STRAINS             = 4;
    int    USER_NSHEAR_STRAINS              = 4;
    double USER_NORMAL_STRAIN_STEP           = 0.005;
    double USER_SHEAR_STRAIN_STEP            = 0.005;   
    // [OBSOLETE] double USER_SHEAR_STRAIN_STEP             = 0.04; 
    // [OBSOLETE] double USER_NORMAL_STRAIN_STEP            = 0.02;
    // [OBSOLETE] double USER_SHEAR_STRAIN_STEP             = 0.02;
    double USER_NIND_STRAIN_DIRS           = 3;
    int    USER_NORMAL_POLYFIT_ORDER       = USER_NNORMAL_STRAINS;
    int    USER_SHEAR_POLYFIT_ORDER        = USER_NSHEAR_STRAINS;
    double USER_SYM_TOL_FRAC               = 0.05;
    string USER_SKIP_FAILED_ARUNS          = "";
    int    USER_SKIP_ARUNS_MAX             = 1;
    ostringstream aus;
    stringstream oss;
    // Error return
    uint aelerror = 0;
    // Tolerance for zero value; if a double is less than this it will be considered to be less than or equal to zero (which could generate an error)
    double tolzero = 1e-12;
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
    vector<double> Pressure, PressureVolumes, VolumeScaleFactors;
    uint stablePressureStructure = 0;
    bool stablePressureFound = false;

    aurostd::StringstreamClean(aus);
    aus << "00000  MESSAGE Starting Elastic Constants RUN" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Get user's parameters from aflow.in 

    // Get user's values of what to calculate and write

    // Check's if the system name is present in the aflow.in file
    // If it exists, it is used to write the output filenames
    if( aurostd::substring2bool(AflowIn,_AFSTROPT_+"SYSTEM=",TRUE) ) {
      USER_SYSTEM_NAME = aurostd::substring2string(AflowIn,_AFSTROPT_+"SYSTEM=",1,TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System name = " << USER_SYSTEM_NAME.c_str() << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      
    // Get the user's selection of whether to write out the results for the elastic constants tensor and elastic moduli. 
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"WRITE_FULL_RESULTS=",TRUE) ) {
      USER_WRITE_FULL_RESULTS.options2entry(AflowIn,_AELSTROPT_+"WRITE_FULL_RESULTS=",USER_WRITE_FULL_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Write full results = " << USER_WRITE_FULL_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"WRITEFULLRESULTS=",TRUE) ) {
      USER_WRITE_FULL_RESULTS.options2entry(AflowIn,_AELSTROPT_+"WRITEFULLRESULTS=",USER_WRITE_FULL_RESULTS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Write full results = " << USER_WRITE_FULL_RESULTS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to use the lattice symmetry to reduce the number of strains which need to be performed.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"STRAIN_SYMMETRY=",TRUE) ) {
      USER_STRAIN_SYMMETRY.options2entry(AflowIn,_AELSTROPT_+"STRAIN_SYMMETRY=",USER_STRAIN_SYMMETRY.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use lattice symmetry to reduce number of strains = " << USER_STRAIN_SYMMETRY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"STRAINSYMMETRY=",TRUE) ) {
      USER_STRAIN_SYMMETRY.options2entry(AflowIn,_AELSTROPT_+"STRAINSYMMETRY=",USER_STRAIN_SYMMETRY.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use lattice symmetry to reduce number of strains = " << USER_STRAIN_SYMMETRY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to use the lattice symmetry to check the symmetry and consistency of the elastic constants.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"CHECK_ELASTIC_SYMMETRY=",TRUE) ) {
      USER_CHECK_ELASTIC_SYMMETRY.options2entry(AflowIn,_AELSTROPT_+"CHECK_ELASTIC_SYMMETRY=",USER_CHECK_ELASTIC_SYMMETRY.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use lattice symmetry to reduce number of strains = " << USER_CHECK_ELASTIC_SYMMETRY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"CHECKELASTICSYMMETRY=",TRUE) ) {
      USER_CHECK_ELASTIC_SYMMETRY.options2entry(AflowIn,_AELSTROPT_+"CHECKELASTICSYMMETRY=",USER_CHECK_ELASTIC_SYMMETRY.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use lattice symmetry to reduce number of strains = " << USER_CHECK_ELASTIC_SYMMETRY.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to fit the origin of the strain-stress curve.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"ORIGIN_STRAIN_FIT=",TRUE) ) {
      USER_ORIGIN_STRAIN_FIT.options2entry(AflowIn,_AELSTROPT_+"ORIGIN_STRAIN_FIT=",USER_ORIGIN_STRAIN_FIT.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Fit origin of stress-strain curve = " << USER_ORIGIN_STRAIN_FIT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"FITSTRAINORIGIN=",TRUE) ) {
      USER_ORIGIN_STRAIN_FIT.options2entry(AflowIn,_AELSTROPT_+"FITSTRAINORIGIN=",USER_ORIGIN_STRAIN_FIT.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Fit origin of stress-strain curve = " << USER_ORIGIN_STRAIN_FIT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to calculate the origin of the strain-stress curve.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"ORIGIN_STRAIN_CALC=",TRUE) ) {
      USER_ORIGIN_STRAIN_CALC.options2entry(AflowIn,_AELSTROPT_+"ORIGIN_STRAIN_CALC=",USER_ORIGIN_STRAIN_CALC.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Calculate origin of stress-strain curve = " << USER_ORIGIN_STRAIN_CALC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"CALCSTRAINORIGIN=",TRUE) ) {
      USER_ORIGIN_STRAIN_CALC.options2entry(AflowIn,_AELSTROPT_+"CALCSTRAINORIGIN=",USER_ORIGIN_STRAIN_CALC.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Calculate origin of stress-strain curve = " << USER_ORIGIN_STRAIN_CALC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to use the existing relaxed structure to fit the origin of the strain-stress curve.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"RELAXED_STRUCT_FIT=",TRUE) ) {
      USER_RELAXED_STRUCT_FIT.options2entry(AflowIn,_AELSTROPT_+"RELAXED_STRUCT_FIT=",USER_RELAXED_STRUCT_FIT.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use existing relaxed structure to fit origin of stress-strain curve = " << USER_RELAXED_STRUCT_FIT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"FITRELAXEDSTRUCT=",TRUE) ) {
      USER_RELAXED_STRUCT_FIT.options2entry(AflowIn,_AELSTROPT_+"FITRELAXEDSTRUCT=",USER_RELAXED_STRUCT_FIT.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use existing relaxed structure to fit origin of stress-strain curve = " << USER_RELAXED_STRUCT_FIT.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to calculate the origin of the strain-stress curve.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"VASPRUNXML_STRESS=",TRUE) ) {
      USER_VASPRUNXML_STRESS.options2entry(AflowIn,_AELSTROPT_+"VASPRUNXML_STRESS=",USER_VASPRUNXML_STRESS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use stress from vasprun.xml file = " << USER_VASPRUNXML_STRESS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"VASPRUNXMLSTRESS=",TRUE) ) {
      USER_VASPRUNXML_STRESS.options2entry(AflowIn,_AELSTROPT_+"VASPRUNXMLSTRESS=",USER_VASPRUNXML_STRESS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use stress from vasprun.xml file = " << USER_VASPRUNXML_STRESS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which values to use for PREC and ALGO.
    // By setting this variable to ON, the values for PREC=ACCURATE and ALGO=NORMAL will be used for AEL, overriding other settings in the aflow.in file.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"PRECACCALGONORM=",TRUE) ) {
      USER_PRECACC_ALGONORM.options2entry(AflowIn,_AELSTROPT_+"PRECACCALGONORM=",USER_PRECACC_ALGONORM.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use PREC=ACCURATE and ALGO=NORMAL = " << USER_PRECACC_ALGONORM.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"PRECACC_ALGONORM=",TRUE) ) {
      USER_PRECACC_ALGONORM.options2entry(AflowIn,_AELSTROPT_+"PRECACC_ALGONORM=",USER_PRECACC_ALGONORM.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use PREC=ACCURATE and ALGO=NORMAL = " << USER_PRECACC_ALGONORM.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which type of runs to do to obtain the elastic constants
    // The default is to use RELAX_STATIC=2, independently of what is used in these lines in the aflow.in file. 
    // By setting RELAX_STATIC to ON, the run type RELAX_STATIC=2 will be used.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"RELAX_STATIC=",TRUE) ) {
      USER_RELAX_STATIC.options2entry(AflowIn,_AELSTROPT_+"RELAX_STATIC=",USER_RELAX_STATIC.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use RELAX_STATIC = " << USER_RELAX_STATIC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which type of runs to do to obtain the elastic constants
    // The default is to use RELAX_STATIC=2, independently of what is used in these lines in the aflow.in file. 
    // By setting STATIC to ON, the run type STATIC will be used.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"STATIC=",TRUE) ) {
      USER_STATIC.options2entry(AflowIn,_AELSTROPT_+"STATIC=",USER_STATIC.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use STATIC = " << USER_STATIC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of which type of runs to do to obtain the elastic constants
    // The default is to use RELAX_STATIC=2, independently of what is used in these lines in the aflow.in file. 
    // By setting RELAX to ON, the run type RELAX=2 will be used.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"RELAX=",TRUE) ) {
      USER_RELAX.options2entry(AflowIn,_AELSTROPT_+"RELAX=",USER_RELAX.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use RELAX = " << USER_RELAX.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to overwrite an existing aflow.in file in subdirectory if no LOCK file is present.
    // The default is to set AFLOWIN_OVERWRITE=OFF, so that aflow.in files in the subdirectories are protected from being automatically overwritten. 
    // By setting AFLOWIN_OVERWRITE to ON, the aflow.in files in the subdirectories will be overwritten if no LOCK file is present.
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"AFLOWIN_OVERWRITE=",TRUE) ) {
      USER_AFLOWIN_OVERWRITE.options2entry(AflowIn,_AELSTROPT_+"AFLOWIN_OVERWRITE=",USER_AFLOWIN_OVERWRITE.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use AFLOWIN_OVERWRITE = " << USER_AFLOWIN_OVERWRITE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"AFLOWINOVERWRITE=",TRUE) ) {
      USER_AFLOWIN_OVERWRITE.options2entry(AflowIn,_AELSTROPT_+"AFLOWINOVERWRITE=",USER_AFLOWIN_OVERWRITE.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use AFLOWINOVERWRITE = " << USER_AFLOWIN_OVERWRITE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // Get user's values of parameters for AGL calculation

    // Get the user's number of normal strains to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NNORMAL_STRAINS=",TRUE) ) {
      USER_NNORMAL_STRAINS = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NNORMAL_STRAINS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Number of normal strains = " << USER_NNORMAL_STRAINS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NNORMALSTRAINS=",TRUE) ) {
      USER_NNORMAL_STRAINS = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NNORMALSTRAINS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Number of normal strains = " << USER_NNORMAL_STRAINS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_NNORMAL_STRAINS < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Number of normal strains = " << USER_NNORMAL_STRAINS << " < 1" << endl;  
      USER_NNORMAL_STRAINS = 4;
      aus << _AELSTR_WARNING_ + "Increasing number of normal strains to default value of " << USER_NNORMAL_STRAINS << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's number of shear strains to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NSHEAR_STRAINS=",TRUE) ) {
      USER_NSHEAR_STRAINS = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NSHEAR_STRAINS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Number of shear strains = " << USER_NSHEAR_STRAINS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NSHEARSTRAINS=",TRUE) ) {
      USER_NSHEAR_STRAINS = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NSHEARSTRAINS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Number of shear strains = " << USER_NSHEAR_STRAINS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_NSHEAR_STRAINS < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Number of shear strains = " << USER_NSHEAR_STRAINS << " < 1" << endl;  
      USER_NSHEAR_STRAINS = 4;
      aus << _AELSTR_WARNING_ + "Increasing number of shear strains to default value of " << USER_NSHEAR_STRAINS << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Get the user's number of polynomial coefficients to be used to fit the normal strain-stress data
    // This includes the zeroth order coefficient, so actual polynomial degree will be one less than this number
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NORMAL_POLYFIT_ORDER=",TRUE) ) {
      USER_NORMAL_POLYFIT_ORDER = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NORMAL_POLYFIT_ORDER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Normal strain polynomial fitting order = " << USER_NORMAL_POLYFIT_ORDER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NORMALPOLYFITORDER=",TRUE) ) {
      USER_NORMAL_POLYFIT_ORDER = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NORMALPOLYFITORDER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Normal strain polynomial fitting order = " << USER_NORMAL_POLYFIT_ORDER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(USER_ORIGIN_STRAIN_FIT.option || USER_RELAXED_STRUCT_FIT.option) {
      USER_NORMAL_POLYFIT_ORDER = USER_NNORMAL_STRAINS + 1;
    }
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_NORMAL_POLYFIT_ORDER < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Normal strain polynomial fitting order = " << USER_NORMAL_POLYFIT_ORDER << " < 1" << endl;  
      if(USER_ORIGIN_STRAIN_FIT.option || USER_RELAXED_STRUCT_FIT.option) {
        USER_NORMAL_POLYFIT_ORDER = USER_NNORMAL_STRAINS + 1;
      } else {
        USER_NORMAL_POLYFIT_ORDER = USER_NNORMAL_STRAINS;
      }	
      aus << _AELSTR_WARNING_ + "Increasing polynomial fitting order to default value of " << USER_NORMAL_POLYFIT_ORDER << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }     
    // Get the user's number of polynomial coefficients to be used to fit the shear strain-stress data
    // This includes the zeroth order coefficient, so actual polynomial degree will be one less than this number
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SHEAR_POLYFIT_ORDER=",TRUE) ) {
      USER_SHEAR_POLYFIT_ORDER = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"SHEAR_POLYFIT_ORDER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Shear strain polynomial fitting order = " << USER_SHEAR_POLYFIT_ORDER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SHEARPOLYFITORDER=",TRUE) ) {
      USER_SHEAR_POLYFIT_ORDER = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"SHEARPOLYFITORDER=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Shear strain polynomial fitting order = " << USER_SHEAR_POLYFIT_ORDER << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if(USER_ORIGIN_STRAIN_FIT.option || USER_RELAXED_STRUCT_FIT.option) {
      USER_SHEAR_POLYFIT_ORDER = USER_NSHEAR_STRAINS + 1;
    }
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_SHEAR_POLYFIT_ORDER < 1) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Shear strain polynomial fitting order = " << USER_SHEAR_POLYFIT_ORDER << " < 1" << endl;  
      if(USER_ORIGIN_STRAIN_FIT.option || USER_RELAXED_STRUCT_FIT.option) {
        USER_SHEAR_POLYFIT_ORDER = USER_NSHEAR_STRAINS + 1;
      } else {
        USER_SHEAR_POLYFIT_ORDER = USER_NSHEAR_STRAINS;
      }
      aus << _AELSTR_WARNING_ + "Increasing shear strain polynomial fitting order to default value of " << USER_SHEAR_POLYFIT_ORDER << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }    
    // Get the user's number of independent strain directions to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NIND_STRAIN_DIRS=",TRUE) ) {
      USER_NIND_STRAIN_DIRS = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NIND_STRAIN_DIRS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Number of independent strain directions = " << USER_NIND_STRAIN_DIRS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NINDSTRAINDIRS=",TRUE) ) {
      USER_NIND_STRAIN_DIRS = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"NINDSTRAINDIRS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Number of independent strain directions = " << USER_NIND_STRAIN_DIRS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Check that value is not zero, negative or greater than 3; warn user and reset to default values if it is zero or negative or greater than 3
    if(USER_NIND_STRAIN_DIRS < 1 || USER_NIND_STRAIN_DIRS > 3) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Number of independent strain directions = " << USER_NNORMAL_STRAINS << " < 1 or > 3" << endl;  
      USER_NIND_STRAIN_DIRS = 3;
      aus << _AELSTR_WARNING_ + "Resetting number of independent strain directions to default value of " << USER_NIND_STRAIN_DIRS << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      
    // Get the user's fractional tolerance for checking symmetry and consistency of elastic tensor
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SYM_TOL_FRAC=",TRUE) ) {
      USER_SYM_TOL_FRAC = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"SYM_TOL_FRAC=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Fractional tolerance for checking symmetry and consistency of elastic tensor = " << USER_SYM_TOL_FRAC << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SYMTOLFRAC=",TRUE) ) {
      USER_SYM_TOL_FRAC = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"SYMTOLFRAC=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Fractional tolerance for checking symmetry and consistency of elastic tensor = " << USER_SYM_TOL_FRAC << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Get the user's normal strain step size to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NORMAL_STRAIN_STEP=",TRUE) ) {
      USER_NORMAL_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"NORMAL_STRAIN_STEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Normal strain step size = " << USER_NORMAL_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NORMALSTRAINSTEP=",TRUE) ) {
      USER_NORMAL_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"NORMALSTRAINSTEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Normal strain step size = " << USER_NORMAL_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_NORMAL_STRAIN_STEP < tolzero) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Strain step size = " << USER_NORMAL_STRAIN_STEP << " <= 0.0" << endl;  
      USER_NORMAL_STRAIN_STEP = 0.005;
      aus << _AELSTR_WARNING_ + "Increasing normal strain step size to default value of " << USER_NORMAL_STRAIN_STEP << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }  
    // Get the user's shear strain step size to be used in the calculation
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SHEAR_STRAIN_STEP=",TRUE) ) {
      USER_SHEAR_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"SHEAR_STRAIN_STEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Shear strain step size = " << USER_SHEAR_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SHEARSTRAINSTEP=",TRUE) ) {
      USER_SHEAR_STRAIN_STEP = aurostd::substring2utype<double>(AflowIn,_AELSTROPT_+"SHEARSTRAINSTEP=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Shear strain step size = " << USER_SHEAR_STRAIN_STEP << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Check that value is not zero or negative; warn user and reset to default values if it is zero or negative
    if(USER_SHEAR_STRAIN_STEP < tolzero) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Shear strain step size = " << USER_SHEAR_STRAIN_STEP << " <= 0.0" << endl;  
      //USER_SHEAR_STRAIN_STEP = 0.04;
      USER_SHEAR_STRAIN_STEP = 0.005;
      aus << _AELSTR_WARNING_ + "Increasing shear strain step size to default value of " << USER_SHEAR_STRAIN_STEP << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }  
    // Get the user's selection of whether to use negative strains. 
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NEG_STRAINS=",TRUE) ) {
      USER_NEG_STRAINS.options2entry(AflowIn,_AELSTROPT_+"NEG_STRAINS=",USER_NEG_STRAINS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Plot AEL results = " << USER_NEG_STRAINS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"NEGSTRAINS=",TRUE) ) {
      USER_NEG_STRAINS.options2entry(AflowIn,_AELSTROPT_+"NEGSTRAINS=",USER_NEG_STRAINS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Plot AEL results = " << USER_NEG_STRAINS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }  	
    // Get the user's selection of whether to get the atomic species from the potential name in order to calculate the unit cell mass
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SPECIES_MASS=",TRUE) ) {
      USER_SPECIES_MASS.options2entry(AflowIn,_AELSTROPT_+"SPECIES_MASS=",USER_SPECIES_MASS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Obtain atomic species from potential name = " << USER_SPECIES_MASS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SPECIESMASS=",TRUE) ) {
      USER_SPECIES_MASS.options2entry(AflowIn,_AELSTROPT_+"SPECIESMASS=",USER_SPECIES_MASS.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Obtain atomic species from potential name = " << USER_SPECIES_MASS.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Get the user's selection of whether to use the prefix "ARUN" before the individual directory names
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"DIRNAME_ARUN=",TRUE) ) {
      USER_DIRNAME_ARUN.options2entry(AflowIn,_AELSTROPT_+"DIRNAME_ARUN=",USER_DIRNAME_ARUN.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Use ARUN prefix in directory names = " << USER_DIRNAME_ARUN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"DIRNAMEARUN=",TRUE) ) {
      USER_DIRNAME_ARUN.options2entry(AflowIn,_AELSTROPT_+"DIRNAMEARUN=",USER_DIRNAME_ARUN.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Use ARUN prefix in directory names = " << USER_DIRNAME_ARUN.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's selection of whether to write the debugging information for the function AGL_functions::gaussxm
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"GAUSSXM_DEBUG=",TRUE) ) {
      USER_GAUSSXM_DEBUG.options2entry(AflowIn,_AELSTROPT_+"GAUSSXM_DEBUG=",USER_GAUSSXM_DEBUG.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Write AGL_functions::gaussxm debugging information = " << USER_GAUSSXM_DEBUG.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"GAUSSXMDEBUG=",TRUE) ) {
      USER_GAUSSXM_DEBUG.options2entry(AflowIn,_AELSTROPT_+"GAUSSXMDEBUG=",USER_GAUSSXM_DEBUG.option);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Write AGL_functions::gaussxm debugging information = " << USER_GAUSSXM_DEBUG.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 	
    // Get the user's list of failed run subdirectories to skip when calculating the elastic constants
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SKIP_FAILED_ARUNS=",TRUE) ) {
      USER_SKIP_FAILED_ARUNS = aurostd::substring2string(AflowIn,_AELSTROPT_+"SKIP_FAILED_ARUNS=",1,TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Failed ARUNS to skip = " << USER_SKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SKIPFAILEDARUNS=",TRUE) ) {
      USER_SKIP_FAILED_ARUNS = aurostd::substring2string(AflowIn,_AELSTROPT_+"SKIPFAILEDARUNS=",1,TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Failed ARUNS to skip = " << USER_SKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 		
    // Get the user's option of whether to automatically skip failed run subdirectories 
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"AUTOSKIP_FAILED_ARUNS=",TRUE) ) {
      USER_AUTOSKIP_FAILED_ARUNS.options2entry(AflowIn,_AELSTROPT_+"AUTOSKIP_FAILED_ARUNS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Automatically detect and skip failed ARUNS = " << USER_AUTOSKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"AUTOSKIPFAILEDARUNS=",TRUE) ) {
      USER_AUTOSKIP_FAILED_ARUNS.options2entry(AflowIn,_AELSTROPT_+"AUTOSKIPFAILEDARUNS=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Automatically detect and skip failed ARUNS = " << USER_AUTOSKIP_FAILED_ARUNS << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 		
    // Get the user's option of maximum number of failed run subdirectories to skip in each independent direction
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SKIP_ARUNS_MAX=",TRUE) ) {
      USER_SKIP_ARUNS_MAX = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"SKIP_ARUNS_MAX=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Maximum number of ARUNS to skip = " << USER_SKIP_ARUNS_MAX << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SKIPARUNSMAX=",TRUE) ) {
      USER_SKIP_ARUNS_MAX = aurostd::substring2utype<int>(AflowIn,_AELSTROPT_+"SKIPARUNSMAX=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Maximum number of ARUNS to skip = " << USER_SKIP_ARUNS_MAX << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 		
    // Get the user's option of whether to calculate the elastic constants under pressure
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"PRESSURE_CALC=",TRUE) ) {
      USER_PRESSURE_CALC.options2entry(AflowIn,_AELSTROPT_+"PRESSURE_CALC=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Calculate elastic constants at pressure = " << USER_PRESSURE_CALC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"PRESSURECALC=",TRUE) ) {
      USER_PRESSURE_CALC.options2entry(AflowIn,_AELSTROPT_+"PRESSURECALC=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Calculate elastic constants at pressure = " << USER_PRESSURE_CALC.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's option of whether to use the elastic constants for a pressure where the material is elastically stable
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"USE_STABLE_PRESSURE=",TRUE) ) {
      USER_USE_STABLE_PRESSURE.options2entry(AflowIn,_AELSTROPT_+"USE_STABLE_PRESSURE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use elastic constants at pressure where material is stable = " << USER_USE_STABLE_PRESSURE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } else if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"USESTABLEPRESSURE=",TRUE) ) {
      USER_USE_STABLE_PRESSURE.options2entry(AflowIn,_AELSTROPT_+"USESTABLEPRESSURE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Use elastic constants at pressure where material is stable = " << USER_USE_STABLE_PRESSURE.option << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 	   	        
    // Get the user's option of whether to symmetrize the elastic stiffness tensor
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"SYMMETRIZE=",TRUE) ) {
      USER_SYMMETRIZE.options2entry(AflowIn,_AELSTROPT_+"SYMMETRIZE=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Symmetrize elastic stiffness tensor = " << USER_SYMMETRIZE.option << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Get the user's option of whether to use the VASP symmetry routines to reduce the cell size
    if( aurostd::substring2bool(AflowIn,_AELSTROPT_+"VASPSYM=",TRUE) ) {
      USER_VASPSYM.options2entry(AflowIn,_AELSTROPT_+"VASPSYM=",TRUE);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Switch off VASP symmetry = " << USER_VASPSYM.option << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }        

    // Set up name of calculation and output filename and directory path
    AEL_data.dirpathname = aurostd::CleanFileName(xvasp.Directory);
    AEL_data.sysname = aurostd::CleanStringASCII(USER_SYSTEM_NAME);
    AEL_data.sysname = aurostd::RemoveWhiteSpaces(AEL_data.sysname);

    // Set up options on what to calculate
    AEL_data.fitstrainorigin = USER_ORIGIN_STRAIN_FIT.option;
    AEL_data.calcstrainorigin = USER_ORIGIN_STRAIN_CALC.option;
    AEL_data.fitrelaxedstruct = USER_RELAXED_STRUCT_FIT.option;
    AEL_data.vasprunxmlstress = USER_VASPRUNXML_STRESS.option;
    AEL_data.precaccalgonorm = USER_PRECACC_ALGONORM.option;
    AEL_data.relax_static = USER_RELAX_STATIC.option;
    AEL_data.static_only = USER_STATIC.option;
    AEL_data.relax_only = USER_RELAX.option;
    AEL_data.negstrain = USER_NEG_STRAINS.option;
    AEL_data.gaussxm_debug = USER_GAUSSXM_DEBUG.option;
    AEL_data.aflowin_overwrite = USER_AFLOWIN_OVERWRITE.option;
    AEL_data.autoskipfailedaruns = USER_AUTOSKIP_FAILED_ARUNS.option;
    AEL_data.symmetrize_elastic_tensor = USER_SYMMETRIZE.option;
    AEL_data.vasp_symmetry = USER_VASPSYM.option;

    // RELAX_STATIC=2 is the default; if either of the other two options are true it should be false; otherwise it should be true
    if(AEL_data.static_only || AEL_data.relax_only) {
      AEL_data.relax_static = false;
    } else {
      AEL_data.relax_static = true;
    }

    // Insert user's value of fractional tolerance for checking symmetry and consistency of elastic tensor into AEL_data class
    AEL_data.symtolfrac = USER_SYM_TOL_FRAC;

    // Insert user's values of polynomial fitting order for stress-strain data into AEL_data class
    AEL_data.normalpolyfitorder = USER_NORMAL_POLYFIT_ORDER;
    AEL_data.shearpolyfitorder = USER_SHEAR_POLYFIT_ORDER;

    // Insert user's values for maximum number of failed runs to skip in each independent direction
    AEL_data.skiparunsmax = USER_SKIP_ARUNS_MAX;

    // Check if user has requested that the elastic constants be calculated as a function of pressure
    if(USER_PRESSURE_CALC.option) {
      // Call AGL to calculate pressure vs. volume values
      aelerror = AGL_functions::Get_VolumeStaticPressure(xvasp, AflowIn, aflags, kflags, vflags, Pressure, PressureVolumes, VolumeScaleFactors, AEL_data.postprocess, FileMESSAGE);
      // Check if calculation of pressure vs. volume values has finished properly
      if(aelerror == 0) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Pressure vs. volume calculation run successfully" << endl;        
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        for (uint k = 0; k < Pressure.size(); k++) {
          aurostd::StringstreamClean(aus);
          // [OBSOLETE] if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
          aus << _AELSTR_MESSAGE_ + "Pressure = " << Pressure.at(k) << "GPa, Volume = " << PressureVolumes.at(k) << "Ang^3" << endl;        
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE] }
        }
      } else if(aelerror == 8) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ << "AGL waiting for other calculations to run" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aelerror;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ << "AGL calculation failed" << endl;
        aus << _AELSTR_WARNING_ << "Cannot calculate elastic properties at finite pressure" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aelerror;
      }
    }

    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Normal polynomial fitting order = " << AEL_data.normalpolyfitorder << endl;        
    aus << _AELSTR_MESSAGE_ + "Shear polynomial fitting order = " << AEL_data.shearpolyfitorder << endl;       
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Output directory name = " << AEL_data.dirpathname << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Check if there are failed ARUN subdirectories to skip; if so tokenize into vector of strings
    if(USER_SKIP_FAILED_ARUNS != "") {
      vector<string> tokens;
      string failed_arun;
      aurostd::string2tokens(USER_SKIP_FAILED_ARUNS, tokens, ",");
      for(uint i = 0; i < tokens.size(); i++) {
        failed_arun = aurostd::RemoveWhiteSpaces(tokens.at(i));
        AEL_data.failed_arun_list.push_back(failed_arun);
      }
      for (uint i = 0; i < AEL_data.failed_arun_list.size(); i++) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ << "Failed ARUN directory = " << AEL_data.failed_arun_list.at(i) << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // [OBSOLETE] return 1;
    }

    // Calculate total mass of unit cell of system in units of kg
    if(USER_SPECIES_MASS.option) {
      // Use species pp names 
      for(uint i = 0; i < xvasp.str.species.size(); i++) {
        AEL_data.cellmasskg = AEL_data.cellmasskg + xvasp.str.num_each_type.at(i) * GetAtomMass(xvasp.str.species.at(i));
      }

      // Write out the names of the atoms used in the calculation
      uint ij = 0;
      for(uint i = 0; i < xvasp.str.species.size(); i++) {
        for (int j = 0; j < xvasp.str.num_each_type.at(i); j++) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species.at(i) << endl; 
          aus << _AELSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species_pp.at(i) << endl; 
          aus << _AELSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species_pp_type.at(i) << endl; 
          aus << _AELSTR_MESSAGE_ << "Atom " << ij << " = " << xvasp.str.species_pp_version.at(i) << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          ij++;
        }
      }
    } else {
      // Use atom names
      for(uint i = 0; i < xvasp.str.atoms.size(); i++) {
        AEL_data.cellmasskg = AEL_data.cellmasskg + GetAtomMass(xvasp.str.atoms.at(i).cleanname);
      }
      // Write out the names of the atoms used in the calculation
      for(uint i = 0; i < xvasp.str.atoms.size(); i++) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ << "Atom " << i << " = " << xvasp.str.atoms.at(i).cleanname << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    // Calculate volume of unit cell of system in units of m^3
    AEL_data.cellvolumem3 = xvasp.str.Volume() * 1e-30;

    // Calculate mass density of material in units of kg / m^3
    AEL_data.mass_density = AEL_data.cellmasskg / AEL_data.cellvolumem3;

    // Save number of atoms in cell
    AEL_data.natomscell = 1.0 * xvasp.str.atoms.size();

    // Write out values for mass, volume and density
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ << "Mass of calculated cell = " <<  AEL_data.cellmasskg << " kg" << endl;  
    aus << _AELSTR_MESSAGE_ << "Volume of calculated cell = " <<  AEL_data.cellvolumem3 << " m^3" << endl;  
    aus << _AELSTR_MESSAGE_ << "Density of material = " <<  AEL_data.mass_density << " kg / m^3" << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Check which Bravais lattice type the structure has
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Bravais lattice type = " << xvasp.str.bravais_lattice_type << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Use Bravais lattice symmetry to set the number of strains if requested by user
    // So far, only cubic symmetry is used as there are issues with identifying the unique directions for other symmetries
    if(USER_STRAIN_SYMMETRY.option) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Using Bravais lattice type to reduce the number of required strains" << endl;        
      aus << _AELSTR_MESSAGE_ + "Bravais lattice type = " << xvasp.str.bravais_lattice_type << endl;  
      aus << _AELSTR_WARNING_ + "Number of independent strains will be reset automatically from initial value of " << USER_NIND_STRAIN_DIRS << endl;
      aus << _AELSTR_WARNING_ + "If you want to specify a different number of strains, please use the line " + _AELSTROPT_ + "STRAINSYMMETRY=OFF in the aflow.in file" << endl;      
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if((xvasp.str.bravais_lattice_type == "FCC") || (xvasp.str.bravais_lattice_type == "BCC") || (xvasp.str.bravais_lattice_type == "CUB")) {
        USER_NIND_STRAIN_DIRS = 1;
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice has cubic symmetry" << endl;        
        aus << _AELSTR_MESSAGE_ + "Number of independent strains set to " << USER_NIND_STRAIN_DIRS << endl;        
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      } else {
        USER_NIND_STRAIN_DIRS = 3;
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice symmetry requires full set of strains " << endl;        
        aus << _AELSTR_MESSAGE_ + "Number of independent strains set to " << USER_NIND_STRAIN_DIRS << endl;        
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else {
      // If not using lattice symmetry to set strains, checks that number of strains are sufficient for that lattice symmetry
      if((xvasp.str.bravais_lattice_type == "FCC") || (xvasp.str.bravais_lattice_type == "BCC") || (xvasp.str.bravais_lattice_type == "CUB")) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice has cubic symmetry" << endl;        
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(USER_NIND_STRAIN_DIRS >= 1) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Number of independent strains = " << USER_NIND_STRAIN_DIRS << " which is sufficient for cubic symmetry" << endl;        
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          USER_NIND_STRAIN_DIRS = 1;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Number of independent strains = " << USER_NIND_STRAIN_DIRS << " which is not sufficient for cubic symmetry" << endl;        
          aus << _AELSTR_MESSAGE_ + "Number of independent strains reset to " << USER_NIND_STRAIN_DIRS << endl;        
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } 
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice symmetry requires full set of strains " << endl;        
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(USER_NIND_STRAIN_DIRS >= 3) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Number of independent strains = " << USER_NIND_STRAIN_DIRS << " which is sufficient for this symmetry" << endl;        
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          USER_NIND_STRAIN_DIRS = 3;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Number of independent strains = " << USER_NIND_STRAIN_DIRS << " which is not sufficient for this symmetry" << endl;        
          aus << _AELSTR_MESSAGE_ + "Number of independent strains reset to " << USER_NIND_STRAIN_DIRS << endl;        
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
    }

    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "starting try loop to create strained structures and run elastic constant calculation" << endl;  
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Create set of strained structures, write aflow.in files to directories, run VASP, collect results, and pass to elastic constants method
    // "try" loop adapted from that in AFLOW APL function DirectMethodPC::runVASPCalculations()
    try {
      // Create strain sets
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "creating strain sets" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(USER_NEG_STRAINS.option) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "negative strains" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        int nstrainsteps = USER_NNORMAL_STRAINS / 2;
        int nstrainstepshalf = nstrainsteps / 2;
        double dnstrainsteps = nstrainsteps;
        double deformation, di;
        for (int i = 0; i < nstrainsteps; i++) {
          di = i;
          deformation = (di - dnstrainsteps) * USER_NORMAL_STRAIN_STEP;
          AEL_data.normal_deformations.push_back(deformation);
        }
        for (int i = 0; i < nstrainsteps; i++) {
          di = i + 1.0;
          deformation = di * USER_NORMAL_STRAIN_STEP;
          AEL_data.normal_deformations.push_back(deformation);
        }
        for (int i = nstrainstepshalf; i < (USER_NNORMAL_STRAINS - nstrainstepshalf); i++) {
          AEL_data.normal_deformations_half.push_back(AEL_data.normal_deformations.at(i));
        }	  
        nstrainsteps = USER_NSHEAR_STRAINS / 2;
        nstrainstepshalf = nstrainsteps / 2;
        dnstrainsteps = nstrainsteps;
        for (int i = 0; i < nstrainsteps; i++) {
          di = i;
          deformation = (di - dnstrainsteps) * USER_SHEAR_STRAIN_STEP;
          AEL_data.shear_deformations.push_back(deformation);
        }
        for (int i = 0; i < nstrainsteps; i++) {
          di = i + 1.0;
          deformation = di * USER_SHEAR_STRAIN_STEP;
          AEL_data.shear_deformations.push_back(deformation);
        }	
        for (int i = nstrainstepshalf; i < (USER_NSHEAR_STRAINS - nstrainstepshalf); i++) {
          AEL_data.shear_deformations_half.push_back(AEL_data.shear_deformations.at(i));
        }
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "positive strains only" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        int nstrainsteps = USER_NNORMAL_STRAINS;
        int nstrainstepshalf = nstrainsteps / 2;
        double deformation, di;
        for (int i = 0; i < nstrainsteps; i++) {
          di = i + 1.0;
          deformation = di * USER_NORMAL_STRAIN_STEP;
          AEL_data.normal_deformations.push_back(deformation);
        }
        for (int i = 0; i < nstrainstepshalf; i++) {
          AEL_data.normal_deformations_half.push_back(AEL_data.normal_deformations.at(i));
        }	  	  
        nstrainsteps = USER_NSHEAR_STRAINS;
        nstrainstepshalf = nstrainsteps / 2;
        for (int i = 0; i < nstrainsteps; i++) {
          di = i + 1.0;
          deformation = di * USER_SHEAR_STRAIN_STEP;
          AEL_data.shear_deformations.push_back(deformation);
        }
        for (int i = 0; i < nstrainstepshalf; i++) {
          AEL_data.shear_deformations_half.push_back(AEL_data.shear_deformations.at(i));
        }
      }

      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "normal_deformations.size() = " << AEL_data.normal_deformations.size() << endl;  
      aus << _AELSTR_MESSAGE_ + "shear_deformations.size() = " << AEL_data.shear_deformations.size() << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // Create strain matrices
      aurostd::xmatrix<double> strain_matrix(3, 3);

      AEL_data.normal_strain.resize(USER_NIND_STRAIN_DIRS);
      AEL_data.shear_strain.resize(USER_NIND_STRAIN_DIRS);

      // Normal strain deformation     
      for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
        for (uint j = 0; j < AEL_data.normal_deformations.size(); j++) {
          strain_matrix = aurostd::identity(strain_matrix);
          strain_matrix[i][i] = strain_matrix[i][i] + AEL_data.normal_deformations.at(j);
          AEL_data.normal_strain.at(i-1).push_back(strain_matrix);
        }
      }

      // Shear strain deformation
      // Note that here the strain is applied asymmetrically, so that epsilon[k][l] = deformation, while epsilon[l][k] = 0
      // Therefore the shear strain elements of the strain tensor in Voigt notation become (epsilon[k][l] + epsilon[l][k]) = deformation
      for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
        uint k = 3 - i;
        uint l = 3;
        if(i == 3) {
          k = 1;
          l = 2;
        }
        for (uint j = 0; j < AEL_data.shear_deformations.size(); j++) {
          strain_matrix = aurostd::identity(strain_matrix);
          strain_matrix[k][l] = strain_matrix[k][l] + AEL_data.shear_deformations.at(j);
          AEL_data.shear_strain.at(i-1).push_back(strain_matrix);
        }
      }

      // Create strained structures
      xstructure initialstructure = xvasp.str;
      xstructure strainedstructure;
      xstructure compressedstructure;
      vector<xstructure> pressurestructures;
      _xvasp vaspRun;
      vector<_xvasp> vaspRuns;
      vector<vector <_xvasp> > vaspRunsPressures;
      vector<_xvasp> vaspRunsPressure;
      vector<_AEL_data> AEL_data_pressures;
      _AEL_data AEL_data_orig;
      _vflags vaspFlags;
      vaspFlags = vflags;
      _kflags kbinFlags;
      kbinFlags = kflags;
      _aflags aflowFlags;
      aflowFlags = aflags;
      vector<string> runname;
      vector<string> runnamepressure;
      vector<vector<string> > runnamepressures;
      vector<string> dirrunname;
      vector<string> dirrunnamepressure;
      vector<vector<string> > dirrunnamepressures;
      vector<string> arunname;
      vector<string> arunnamepressure;
      vector<vector<string> > arunnamepressures;
      bool avasp_pop_xvasp_success = true;
      bool aflowin_success = true;
      stringstream arunaflowin;
      bool write = true;

      // Initialize vaspRun to get values of directory path, etc.
      vaspRun = xvasp;

      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Initial structure lattice:" << endl; 
      aus << initialstructure.lattice << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "creating strained structures" << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // Generates structures as a function of pressure
      if(USER_PRESSURE_CALC.option) {
        for (uint i = 0; i < VolumeScaleFactors.size(); i++) {
          compressedstructure = initialstructure;
          compressedstructure.InflateVolume(VolumeScaleFactors.at(i));
          pressurestructures.push_back(compressedstructure);
        }
        for (uint i = 0; i < VolumeScaleFactors.size(); i++) {
          aurostd::StringstreamClean(aus);
          // [OBSOLETE] if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
          aus << _AELSTR_MESSAGE_ + "Pressure = " << Pressure.at(i) << ", Volume = " << PressureVolumes.at(i) << endl; 
          aus << _AELSTR_MESSAGE_ + "VolumeScaleFactor = " << VolumeScaleFactors.at(i) << ", Rescaled Volume = " << pressurestructures.at(i).Volume() << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE] }
        }
      }

      // Applies normal strain to structures
      for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
        for (uint j = 0; j < AEL_data.normal_deformations.size(); j++) {
          strainedstructure = initialstructure;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "AEL_data.normal_strain.at(i-1).at(j) = " << endl; 
          aus << AEL_data.normal_strain.at(i-1).at(j) << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          strainedstructure.lattice = strainedstructure.lattice * AEL_data.normal_strain.at(i-1).at(j);
          vaspRuns.push_back(vaspRun);
          uint idVaspRun = vaspRuns.size()-1;
          vaspRuns.at(idVaspRun).str = strainedstructure;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "idVaspRun = " << idVaspRun << endl; 
          aus << _AELSTR_MESSAGE_ + "vaspRuns.at(idVaspRun).str.lattice = " << endl; 
          aus << vaspRuns.at(idVaspRun).str.lattice << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // Set up name of separate directory for AFLOW run for each strained structure
          string stridVaspRun;
          ostringstream cnvidVaspRun;
          cnvidVaspRun << idVaspRun;
          stridVaspRun = cnvidVaspRun.str();
          double strainfactor = 1.0 + AEL_data.normal_deformations.at(j);
          string strstrainfactor;
          ostringstream cnvstrainfactor;
          cnvstrainfactor << strainfactor;
          strstrainfactor = cnvstrainfactor.str();
          string strnormindstrain;
          ostringstream cnvnormindstrain;
          cnvnormindstrain << i;
          strnormindstrain = cnvnormindstrain.str();
          if(USER_DIRNAME_ARUN.option) {
            runname.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
          } else {
            runname.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
          }
          arunname.push_back(stridVaspRun + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
          dirrunname.push_back(AEL_data.dirpathname + "/" + runname.at(idVaspRun));	  	  
          // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AEL_data.dirpathname + "/" + runname.at(idVaspRun);
          vaspRuns.at(idVaspRun).str = strainedstructure;
          AEL_data.strain_matrix_list.push_back(AEL_data.normal_strain.at(i-1).at(j));
          vaspRuns.at(idVaspRun).AVASP_arun_runname = arunname.at(idVaspRun);
          vaspRuns.at(idVaspRun).AVASP_arun = true;
          vaspRuns.at(idVaspRun).AVASP_arun_mode = "AEL";	  
        }
      }

      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "creating shear strained structures" << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // Applies shear strain to structures
      for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
        for (uint j = 0; j < AEL_data.shear_deformations.size(); j++) {
          strainedstructure = initialstructure;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "i = " << i << ", j = " << j << endl;
          aus << _AELSTR_MESSAGE_ + "AEL_data.shear_strain.at(i-1).at(j) = " << endl;
          aus << AEL_data.shear_strain.at(i-1).at(j) << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          strainedstructure.lattice = strainedstructure.lattice * AEL_data.shear_strain.at(i-1).at(j);
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "strainedstructure.lattice = " << endl; 
          aus << strainedstructure.lattice << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          vaspRuns.push_back(vaspRun);
          uint idVaspRun = vaspRuns.size()-1;
          vaspRuns.at(idVaspRun).str = strainedstructure;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "idVaspRun = " << idVaspRun << endl; 
          aus << _AELSTR_MESSAGE_ + "vaspRuns.at(idVaspRun).str.lattice = " << endl; 
          aus << vaspRuns.at(idVaspRun).str.lattice << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // Set up name of separate directory for AFLOW run for each strained structure
          string stridVaspRun;
          ostringstream cnvidVaspRun;
          cnvidVaspRun << idVaspRun;
          stridVaspRun = cnvidVaspRun.str();
          double strainfactor = 1.0 + AEL_data.shear_deformations.at(j);
          string strstrainfactor;
          ostringstream cnvstrainfactor;
          cnvstrainfactor << strainfactor;
          strstrainfactor = cnvstrainfactor.str();
          string strshearindstrain;
          ostringstream cnvshearindstrain;
          cnvshearindstrain << i;
          strshearindstrain = cnvshearindstrain.str();
          if(USER_DIRNAME_ARUN.option) {
            runname.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
          } else {
            runname.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
          }
          arunname.push_back(stridVaspRun + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
          dirrunname.push_back(AEL_data.dirpathname + "/" + runname.at(idVaspRun));
          // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AEL_data.dirpathname + "/" + runname.at(idVaspRun);
          vaspRuns.at(idVaspRun).str = strainedstructure;
          AEL_data.strain_matrix_list.push_back(AEL_data.shear_strain.at(i-1).at(j));
          vaspRuns.at(idVaspRun).AVASP_arun_runname = arunname.at(idVaspRun);
          vaspRuns.at(idVaspRun).AVASP_arun = true;
          vaspRuns.at(idVaspRun).AVASP_arun_mode = "AEL";
        }
      }

      // Sets up calculation for unstrained structure
      if(AEL_data.calcstrainorigin) {
        strainedstructure = initialstructure;
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "strainedstructure.lattice = " << endl; 
        aus << strainedstructure.lattice << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        vaspRuns.push_back(vaspRun);
        uint idVaspRun = vaspRuns.size()-1;
        vaspRuns.at(idVaspRun).str = strainedstructure;
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "idVaspRun = " << idVaspRun << endl; 
        aus << _AELSTR_MESSAGE_ + "vaspRuns.at(idVaspRun).str.lattice = " << endl; 
        aus << vaspRuns.at(idVaspRun).str.lattice << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // Set up name of separate directory for AFLOW run for each strained structure
        string stridVaspRun;
        ostringstream cnvidVaspRun;
        cnvidVaspRun << idVaspRun;
        stridVaspRun = cnvidVaspRun.str();
        double strainfactor = 1.0;
        string strstrainfactor;
        ostringstream cnvstrainfactor;
        cnvstrainfactor << strainfactor;
        strstrainfactor = cnvstrainfactor.str();
        string strorgindstrain;
        ostringstream cnvorgindstrain;
        cnvorgindstrain << 0;
        strorgindstrain = cnvorgindstrain.str();
        if(USER_DIRNAME_ARUN.option) {
          runname.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
        } else {
          runname.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
        }
        arunname.push_back(stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
        dirrunname.push_back(AEL_data.dirpathname + "/" + runname.at(idVaspRun));
        // [OBSOLETE] vaspRuns.at(idVaspRun).Directory = AEL_data.dirpathname + "/" + runname.at(idVaspRun);
        vaspRuns.at(idVaspRun).str = strainedstructure;
        AEL_data.strain_matrix_list.push_back(aurostd::identity(strain_matrix));
        vaspRuns.at(idVaspRun).AVASP_arun_runname = arunname.at(idVaspRun);
        vaspRuns.at(idVaspRun).AVASP_arun = true;
        vaspRuns.at(idVaspRun).AVASP_arun_mode = "AEL";
      }

      // Set-up normal and shear strained structures for finite pressures
      if(USER_PRESSURE_CALC.option) {
        if(aurostd::abs(Pressure.at(0)) > tolzero) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Calculation of elastic constants at finite pressure" << endl; 
          aus << _AELSTR_ERROR_ + "First pressure value = " << Pressure.at(0) << " > 0.0" << endl; 
          aus << _AELSTR_ERROR_ + "First pressure value needs to be equal to 0.0" << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aelerror = 1;
        }
        // [OBSOLETE] for (uint k = 1; k < Pressure.size(); k++)
        for (uint k = 0; k < Pressure.size(); k++) {
          vaspRunsPressure.clear();
          runnamepressure.clear();
          arunnamepressure.clear();
          dirrunnamepressure.clear();
          compressedstructure = pressurestructures.at(k);
          string strpressure;
          ostringstream cnvstrpressure;
          cnvstrpressure << Pressure.at(k);
          strpressure = cnvstrpressure.str();	
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "creating normal strained structures at P = " << Pressure.at(k) << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // Applies normal strain to structures
          for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
            for (uint j = 0; j < AEL_data.normal_deformations.size(); j++) {
              strainedstructure = compressedstructure;
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "AEL_data.normal_strain.at(i-1).at(j) = " << endl; 
              aus << AEL_data.normal_strain.at(i-1).at(j) << endl; 
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              strainedstructure.lattice = strainedstructure.lattice * AEL_data.normal_strain.at(i-1).at(j);
              vaspRunsPressure.push_back(vaspRun);
              uint idVaspRun = vaspRunsPressure.size()-1;
              vaspRunsPressure.at(idVaspRun).str = strainedstructure;
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "idVaspRun = " << idVaspRun << endl; 
              aus << _AELSTR_MESSAGE_ + "vaspRunsPressure.at(idVaspRun).str.lattice = " << endl; 
              aus << vaspRunsPressure.at(idVaspRun).str.lattice << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // Set up name of separate directory for AFLOW run for each strained structure
              string stridVaspRun;
              ostringstream cnvidVaspRun;
              cnvidVaspRun << idVaspRun;
              stridVaspRun = cnvidVaspRun.str();
              double strainfactor = 1.0 + AEL_data.normal_deformations.at(j);
              string strstrainfactor;
              ostringstream cnvstrainfactor;
              cnvstrainfactor << strainfactor;
              strstrainfactor = cnvstrainfactor.str();
              string strnormindstrain;
              ostringstream cnvnormindstrain;
              cnvnormindstrain << i;
              strnormindstrain = cnvnormindstrain.str();
              if(USER_DIRNAME_ARUN.option) {
                if(aurostd::abs(Pressure.at(k)) < tolzero) {
                  runnamepressure.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
                } else {
                  runnamepressure.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_P_" + strpressure + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
                }
              } else {
                if(aurostd::abs(Pressure.at(k)) < tolzero) {
                  runnamepressure.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
                } else {
                  runnamepressure.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_P_" + strpressure + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
                }
              }
              if(aurostd::abs(Pressure.at(k)) < tolzero) {
                arunnamepressure.push_back(stridVaspRun + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
              } else {
                arunnamepressure.push_back(stridVaspRun + "_P_" + strpressure + "_SF_N_" + strnormindstrain + "_" + strstrainfactor);
              }
              dirrunnamepressure.push_back(AEL_data.dirpathname + "/" + runnamepressure.at(idVaspRun));
              // [OBSOLETE] vaspRunsPressure.at(idVaspRun).Directory = AEL_data.dirpathname + "/" + runnamepressure.at(idVaspRun);
              vaspRunsPressure.at(idVaspRun).str = strainedstructure;
              vaspRunsPressure.at(idVaspRun).AVASP_arun_runname = arunnamepressure.at(idVaspRun);
              vaspRunsPressure.at(idVaspRun).AVASP_arun = true;
              vaspRunsPressure.at(idVaspRun).AVASP_arun_mode = "AEL";
            }
          }

          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "creating shear strained structures at P = " << Pressure.at(k) << endl; 
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

          // Applies shear strain to structures
          for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
            for (uint j = 0; j < AEL_data.shear_deformations.size(); j++) {
              strainedstructure = compressedstructure;
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "i = " << i << ", j = " << j << endl;
              aus << _AELSTR_MESSAGE_ + "AEL_data.shear_strain.at(i-1).at(j) = " << endl;
              aus << AEL_data.shear_strain.at(i-1).at(j) << endl; 
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              strainedstructure.lattice = strainedstructure.lattice * AEL_data.shear_strain.at(i-1).at(j);
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "strainedstructure.lattice = " << endl; 
              aus << strainedstructure.lattice << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              vaspRunsPressure.push_back(vaspRun);
              uint idVaspRun = vaspRunsPressure.size()-1;
              vaspRunsPressure.at(idVaspRun).str = strainedstructure;
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "idVaspRun = " << idVaspRun << endl; 
              aus << _AELSTR_MESSAGE_ + "vaspRunsPressure.at(idVaspRun).str.lattice = " << endl; 
              aus << vaspRunsPressure.at(idVaspRun).str.lattice << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // Set up name of separate directory for AFLOW run for each strained structure
              string stridVaspRun;
              ostringstream cnvidVaspRun;
              cnvidVaspRun << idVaspRun;
              stridVaspRun = cnvidVaspRun.str();
              double strainfactor = 1.0 + AEL_data.shear_deformations.at(j);
              string strstrainfactor;
              ostringstream cnvstrainfactor;
              cnvstrainfactor << strainfactor;
              strstrainfactor = cnvstrainfactor.str();
              string strshearindstrain;
              ostringstream cnvshearindstrain;
              cnvshearindstrain << i;
              strshearindstrain = cnvshearindstrain.str();
              if(USER_DIRNAME_ARUN.option) {
                if(aurostd::abs(Pressure.at(k)) < tolzero) {
                  runnamepressure.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
                } else {
                  runnamepressure.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_P_" + strpressure + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
                }
              } else {
                if(aurostd::abs(Pressure.at(k)) < tolzero) {
                  runnamepressure.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
                } else {
                  runnamepressure.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_P_" + strpressure + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
                }
              }
              if(aurostd::abs(Pressure.at(k)) < tolzero) {
                arunnamepressure.push_back(stridVaspRun + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
              } else {
                arunnamepressure.push_back(stridVaspRun + "_P_" + strpressure + "_SF_S_" + strshearindstrain + "_" + strstrainfactor);
              }
              dirrunnamepressure.push_back(AEL_data.dirpathname + "/" + runnamepressure.at(idVaspRun));
              // [OBSOLETE] vaspRunsPressure.at(idVaspRun).Directory = AEL_data.dirpathname + "/" + runnamepressure.at(idVaspRun);
              vaspRunsPressure.at(idVaspRun).str = strainedstructure;
              vaspRunsPressure.at(idVaspRun).AVASP_arun_runname = arunnamepressure.at(idVaspRun);
              vaspRunsPressure.at(idVaspRun).AVASP_arun = true;
              vaspRunsPressure.at(idVaspRun).AVASP_arun_mode = "AEL";
            }
          }

          // Sets up calculation for unstrained structure
          if(AEL_data.calcstrainorigin) {
            strainedstructure = compressedstructure;
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "strainedstructure.lattice = " << endl; 
            aus << strainedstructure.lattice << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            vaspRunsPressure.push_back(vaspRun);
            uint idVaspRun = vaspRunsPressure.size()-1;
            vaspRunsPressure.at(idVaspRun).str = strainedstructure;
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "idVaspRun = " << idVaspRun << endl; 
            aus << _AELSTR_MESSAGE_ + "vaspRunsPressure.at(idVaspRun).str.lattice = " << endl; 
            aus << vaspRunsPressure.at(idVaspRun).str.lattice << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            // Set up name of separate directory for AFLOW run for each strained structure
            string stridVaspRun;
            ostringstream cnvidVaspRun;
            cnvidVaspRun << idVaspRun;
            stridVaspRun = cnvidVaspRun.str();
            double strainfactor = 1.0;
            string strstrainfactor;
            ostringstream cnvstrainfactor;
            cnvstrainfactor << strainfactor;
            strstrainfactor = cnvstrainfactor.str();
            string strorgindstrain;
            ostringstream cnvorgindstrain;
            cnvorgindstrain << 0;
            strorgindstrain = cnvorgindstrain.str();
            if(USER_DIRNAME_ARUN.option) {
              if(aurostd::abs(Pressure.at(k)) < tolzero) {
                runnamepressure.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
              } else {
                runnamepressure.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_P_" + strpressure + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
              }
            } else {
              if(aurostd::abs(Pressure.at(k)) < tolzero) {
                runnamepressure.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
              } else {
                runnamepressure.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_P_" + strpressure + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
              }
            }
            if(aurostd::abs(Pressure.at(k)) < tolzero) {
              arunnamepressure.push_back(stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
            } else {
              arunnamepressure.push_back(stridVaspRun + "_P_" + strpressure + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
            }
            dirrunnamepressure.push_back(AEL_data.dirpathname + "/" + runnamepressure.at(idVaspRun));
            // [OBSOLETE] vaspRunsPressure.at(idVaspRun).Directory = AEL_data.dirpathname + "/" + runnamepressure.at(idVaspRun);
            vaspRunsPressure.at(idVaspRun).str = strainedstructure;
            vaspRunsPressure.at(idVaspRun).AVASP_arun_runname = arunnamepressure.at(idVaspRun);
            vaspRunsPressure.at(idVaspRun).AVASP_arun = true;
            vaspRunsPressure.at(idVaspRun).AVASP_arun_mode = "AEL";	    
          }
          vaspRunsPressures.push_back(vaspRunsPressure);
          runnamepressures.push_back(runnamepressure);
          arunnamepressures.push_back(arunnamepressure);
          dirrunnamepressures.push_back(dirrunnamepressure);
          AEL_data_pressures.push_back(AEL_data);
          if (AEL_data_pressures.size() < (k + 1)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_ERROR_ + "Failed to create correct number of pressure configurations" << endl;
            aus << _AELSTR_ERROR_ + "Size of pressure configurations vector = " << AEL_data_pressures.size() << endl;
            aus << _AELSTR_ERROR_ + "Number of pressure configurations generated = " << k + 1 << endl;	    
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          } else {
            AEL_data_pressures.at(k).applied_pressure = Pressure.at(k);
          }	  
        }
      }

      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "setting up AFLOW runs for strained structures" << endl; 
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      bool skipdir = false;
      string dfilename;
      // Set up AFLOW runs for strained structures
      for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
        // Assign the values of the flags provided by the user in the aflow.in file to the class containing the input data for the VASP run
        aelerror = AEL_functions::aelvaspflags(vaspRuns.at(idVaspRun), vaspFlags, kbinFlags, dirrunname.at(idVaspRun), AEL_data, FileMESSAGE);
        if(aelerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Failed to assign values of flags from aflow.in file" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aelerror;
        }
        // [OBSOLETE] aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry = " << vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry << endl;
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme = " << vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme << endl;  
        // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        avasp_pop_xvasp_success = AVASP_populateXVASP(aflowFlags, kbinFlags, vaspFlags, vaspRuns.at(idVaspRun));
        if (!avasp_pop_xvasp_success) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Failed to set AFLOW flags for ARUN " << arunname.at(idVaspRun) << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        // [OBSOLETE] aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry = " << vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry << endl;
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme = " << vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme << endl;  
        // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		

        // If there are already LOCK or OUTCAR.static files in this directory, it means this directory was already generated and computed.
        // Therefore do not touch, but store this structure in the list, so that it can be used in the next part of code.
        // [OBSOLETE] if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/LOCK") ) ||
        if(AEL_data.relax_static || AEL_data.static_only) {
          if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ||
              aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
              aurostd::EFileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ||
              aurostd::EFileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.static") ) ||
              ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
               ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) ||
                (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/ael.LOCK")) ||
                (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
                (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
                (aurostd::FileExist( dirrunname.at(idVaspRun) + "/ael.LOCK")) ||
                (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK")))) ) continue; 	 
        } else {
          if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ||
              aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.relax2") ) ||
              aurostd::EFileExist( vaspRuns.at(idVaspRun).Directory + string("/OUTCAR.relax2") )||
              aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) ||
              aurostd::FileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.relax2") ) ||
              aurostd::EFileExist( dirrunname.at(idVaspRun) + string("/OUTCAR.relax2") ) ||
              ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
               ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) ||
                (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/ael.LOCK")) ||
                (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
                (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
                (aurostd::FileExist( dirrunname.at(idVaspRun) + "/ael.LOCK")) ||
                (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK"))))  ) continue; 	 
        }

        // Check if structure is on list of failed runs to be skipped
        // If so, then skip generation and continue to next structure
        for (uint ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
          if(aurostd::substring2bool(runname.at(idVaspRun),AEL_data.failed_arun_list.at(ij),TRUE)) {
            skipdir = true;
          }
        }	
        if(skipdir) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Skipping directory: " << runname.at(idVaspRun) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = false;
          continue;
        }

        // If files do not exist, and the postprocess flag is not set, continue on to prepare generation of _AFLOWIN_ ...
        if (!(XHOST.POSTPROCESS || XHOST.ARUN_POSTPROCESS || AEL_data.postprocess)) { //CO20210701
          // Assign the values of the flags provided by the user in the aflow.in file to the class containing the input data for the VASP run
          // [OBSOLETE] aelerror = AEL_functions::aelvaspflags(vaspRuns.at(idVaspRun), _vaspFlags, _kbinFlags, runname.at(idVaspRun), AEL_data, FileMESSAGE);
          // [OBSOLETE] if(aelerror != 0) {
          // [OBSOLETE]   aurostd::StringstreamClean(aus);
          // [OBSOLETE]   aus << _AELSTR_ERROR_ + "Failed to assign values of flags from aflow.in file" << endl;  
          // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE]   return aelerror;
          // [OBSOLETE] }

          // Create aflow.in
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "writing aflow.in" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

          // Writes aflow.in file to appropriate directory for each VASP run
          // [OBSOLETE] aelerror = AEL_functions::createAFLOWIN(vaspRuns.at(idVaspRun), xvasp, _kbinFlags, _vaspFlags, AEL_data, FileMESSAGE);
          aflowin_success = AVASP_MakeSingleAFLOWIN(vaspRuns.at(idVaspRun), arunaflowin, write);
          if(!aflowin_success) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_ERROR_ + "Failed to create aflow.in files" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 1;
          }
        }
      }

      // Set up AFLOW runs for normal and shear strained structures for different applied pressures
      if(USER_PRESSURE_CALC.option) {
        for (uint k = 0; k < vaspRunsPressures.size(); k++) {
          // Set up AFLOW runs for strained structures
          for(uint idVaspRun = 0; idVaspRun < vaspRunsPressures.at(k).size(); idVaspRun++) {
            // Assign the values of the flags provided by the user in the aflow.in file to the class containing the input data for the VASP run
            aelerror = AEL_functions::aelvaspflags(vaspRunsPressures.at(k).at(idVaspRun), vaspFlags, kbinFlags, dirrunnamepressures.at(k).at(idVaspRun), AEL_data, FileMESSAGE);
            if(aelerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_ERROR_ + "Failed to assign values of flags from aflow.in file" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aelerror;
            }
            // If there are already LOCK or OUTCAR.static files in this directory, it means this directory was already generated and computed.
            // Therefore do not touch, but store this structure in the list, so that it can be used in the next part of code.
            // [OBSOLETE] if( aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/LOCK") ) ||
            if(AEL_data.relax_static || AEL_data.static_only) {
              if( aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ||
                  aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
                  aurostd::EFileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + string("/OUTCAR.static") ) ||
                  aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/" + _AFLOWLOCK_ ) ||
                  aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + string("/OUTCAR.static") ) ||
                  aurostd::EFileExist( dirrunnamepressures.at(k).at(idVaspRun) + string("/OUTCAR.static") ) ||
                  ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
                   ((aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/agl.LOCK")) ||
                    (aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/ael.LOCK")) ||
                    (aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/LOCK")) ||
                    (aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/agl.LOCK")) ||
                    (aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/ael.LOCK")) ||
                    (aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/LOCK")))) ) continue; 	 
            } else {
              if( aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) ||
                  aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + string("/OUTCAR.relax2") ) ||
                  aurostd::EFileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + string("/OUTCAR.relax2") ) ||
                  aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/" + _AFLOWLOCK_ ) ||
                  aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + string("/OUTCAR.relax2") ) ||
                  aurostd::EFileExist( dirrunnamepressures.at(k).at(idVaspRun) + string("/OUTCAR.relax2") ) ||
                  ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
                   ((aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/agl.LOCK")) ||
                    (aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/ael.LOCK")) ||
                    (aurostd::FileExist( vaspRunsPressures.at(k).at(idVaspRun).Directory + "/LOCK")) ||
                    (aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/agl.LOCK")) ||
                    (aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/ael.LOCK")) ||
                    (aurostd::FileExist( dirrunnamepressures.at(k).at(idVaspRun) + "/LOCK")))) ) continue; 	 
            }

            // Check if structure is on list of failed runs to be skipped
            // If so, then skip generation and continue to next structure
            for (uint ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
              if(aurostd::substring2bool(runnamepressures.at(k).at(idVaspRun),AEL_data.failed_arun_list.at(ij),TRUE)) {
                skipdir = true;
              }
            }	
            if(skipdir) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << runnamepressures.at(k).at(idVaspRun) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              skipdir = false;
              continue;
            }

            // If files do not exist, and the postprocess flag is not set, continue on to prepare generation of _AFLOWIN_ ...
            if (!(XHOST.POSTPROCESS || XHOST.ARUN_POSTPROCESS || AEL_data.postprocess)) { //CO20210701
              // Assign the values of the flags provided by the user in the aflow.in file to the class containing the input data for the VASP run
              // [OBSOLETE] aelerror = AEL_functions::aelvaspflags(vaspRunsPressures.at(k).at(idVaspRun), _vaspFlags, _kbinFlags, runnamepressures.at(k).at(idVaspRun), AEL_data, FileMESSAGE);
              // [OBSOLETE] if(aelerror != 0) {
              // [OBSOLETE]   aurostd::StringstreamClean(aus);
              // [OBSOLETE]   aus << _AELSTR_ERROR_ + "Failed to assign values of flags from aflow.in file" << endl;  
              // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE]   return aelerror;
              // [OBSOLETE] }

              // Create aflow.in
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "writing aflow.in" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

              // Writes aflow.in file to appropriate directory for each VASP run
              // [OBSOLETE] aelerror = AEL_functions::createAFLOWIN(vaspRunsPressures.at(k).at(idVaspRun), xvasp, _kbinFlags, _vaspFlags, AEL_data, FileMESSAGE);
              aflowin_success = AVASP_MakeSingleAFLOWIN(vaspRuns.at(idVaspRun), arunaflowin, write);
              if(!aflowin_success) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_ERROR_ + "Failed to create aflow.in files" << endl;  
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                return 1;
              }
            }
          }
        }
      }	  

      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "extract stress tensors" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // Extract stress tensors from results of VASP calculations for strained structures
      aelerror = AEL_functions::extractstress(vaspRuns, AEL_data, dirrunname, FileMESSAGE);
      if(aelerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Failed to extract stress tensors" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aelerror;
      }
      // Check number of skipped runs in each direction does not exceed limit
      uint numskippedruns;
      for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
        numskippedruns = AEL_data.normal_deformations.size() - AEL_data.normal_deformations_complete.at(i-1).size();
        if(numskippedruns > AEL_data.skiparunsmax) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Number of failed runs in one independent directions = " << numskippedruns << endl;
          aus << _AELSTR_ERROR_ + "Number of failed runs exceeds limit in one independent direction of " << AEL_data.skiparunsmax << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
      }
      for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
        numskippedruns = AEL_data.shear_deformations.size() - AEL_data.shear_deformations_complete.at(i-1).size();
        if(numskippedruns > AEL_data.skiparunsmax) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Number of failed runs in one independent directions = " << numskippedruns << endl;
          aus << _AELSTR_ERROR_ + "Number of failed runs exceeds limit per independent directions of " << AEL_data.skiparunsmax << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
      }

      if(USER_WRITE_FULL_RESULTS.option) {
        // Write structures to JSON format
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ << "Opening file AEL_structures.json" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        string strjson = xstructure2json(xvasp.str);
        oss << "{\"AEL_structures\":[";
        oss << "{\"distorted_structure\":[{\"distortion\":null},"; 
        oss << "{\"structure\":" << strjson << "}]}";
        for (uint i = 0; i < vaspRuns.size(); i++) {
          strjson = xstructure2json(vaspRuns.at(i).str);
          oss << ",{\"distorted_structure\":[{\"distortion\":\"" << runname.at(i) << "\"},";
          oss << "{\"structure\":" << strjson << "}]}";
        }
        oss << "]}";
        string ofilestrjson = AEL_data.dirpathname + "/AEL_structures.json";
        if(!aurostd::stringstream2file(oss, ofilestrjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Unable to open file AEL_structure.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
      }

      if (AEL_data.fitrelaxedstruct) {
        AEL_data.strain_matrix_list.push_back(aurostd::identity(strain_matrix));
        vaspRuns.push_back(vaspRun);
        uint idVaspRun = vaspRuns.size()-1;
        vaspRuns.at(idVaspRun).str = initialstructure;
        string stridVaspRun;
        ostringstream cnvidVaspRun;
        cnvidVaspRun << idVaspRun;
        stridVaspRun = cnvidVaspRun.str();
        double strainfactor = 1.0;
        string strstrainfactor;
        ostringstream cnvstrainfactor;
        cnvstrainfactor << strainfactor;
        strstrainfactor = cnvstrainfactor.str();
        string strorgindstrain;
        ostringstream cnvorgindstrain;
        cnvorgindstrain << 0;
        strorgindstrain = cnvorgindstrain.str();
        if(USER_DIRNAME_ARUN.option) {
          runname.push_back(_ARAELSTR_DIRNAME_ + stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
        } else {
          runname.push_back(_AELSTR_DIRNAME_ + stridVaspRun + "_SF_O_" + strorgindstrain + "_" + strstrainfactor);
        }
      }

      // Write structures and corresponding energies, pressures and stress tensors to JSON format
      // First check that the number of strains is equal to the number of runs
      aurostd::xmatrix<double> stress_tensor(3, 3);
      if (AEL_data.structurecalculated.size() < 1) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "No runs to write out yet" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 8;
      } else {
        uint idSuccessRun = AEL_data.structurecalculated.at(0); 
        string strjson = xstructure2json(vaspRuns.at(idSuccessRun).str);
        oss << "{\"AEL_energy_structures\":[";
        oss << "{\"distortion\":\"" << runname.at(idSuccessRun) << "\",";
        strain_matrix = AEL_data.strain_matrix_list.at(idSuccessRun);
        oss << "\"strain_matrix\":[";
        for (uint j = 1; j < 4; j++) {
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "j = " << j << endl;
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if (j == 1) {
            oss << "[";
          } else {
            oss << ",[";
          }
          for (uint k = 1; k < 4; k++) {
            // [OBSOLETE] aurostd::StringstreamClean(aus);
            // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "k = " << k << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            // Write stress tensor components, converting from kB to GPa
            if (k == 1) {
              oss << 0.1 * strain_matrix[j][k];
            } else {
              oss << "," << 0.1 * strain_matrix[j][k];
            }
          }
          oss << "]";	
        }
        oss << "],";	
        oss << "\"energy\":" << AEL_data.energycalculated.at(0) / AEL_data.natomscell << ",";
        // Convert pressure from kB to GPa before writing
        oss << "\"pressure\":" << 0.1 * AEL_data.pressurecalculated.at(0) << ",";
        stress_tensor = AEL_data.stresscalculated.at(0);
        oss << "\"stress_tensor\":[";
        for (uint j = 1; j < 4; j++) {
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "j = " << j << endl;
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if (j == 1) {
            oss << "[";
          } else {
            oss << ",[";
          }
          for (uint k = 1; k < 4; k++) {
            // [OBSOLETE] aurostd::StringstreamClean(aus);
            // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "k = " << k << endl;
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
        for (uint i = 1; i < AEL_data.structurecalculated.size(); i++) {
          idSuccessRun = AEL_data.structurecalculated.at(i);
          strjson = xstructure2json(vaspRuns.at(idSuccessRun).str);
          oss << ",{\"distortion\":\"" << runname.at(idSuccessRun) << "\",";
          strain_matrix = AEL_data.strain_matrix_list.at(idSuccessRun);
          oss << "\"strain_matrix\":[";
          for (uint j = 1; j < 4; j++) {
            // [OBSOLETE] aurostd::StringstreamClean(aus);
            // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "j = " << j << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if (j == 1) {
              oss << "[";
            } else {
              oss << ",[";
            }
            for (uint k = 1; k < 4; k++) {
              // [OBSOLETE] aurostd::StringstreamClean(aus);
              // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "k = " << k << endl;
              // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              if (k == 1) {
                oss << strain_matrix[j][k];
              } else {
                oss << "," << strain_matrix[j][k];
              }
            }
            oss << "]";	
          }
          oss << "],";		
          oss << "\"energy\":" << AEL_data.energycalculated.at(i) / AEL_data.natomscell << ",";
          // Convert pressure from kB to GPa before writing
          oss << "\"pressure\":" << 0.1 * AEL_data.pressurecalculated.at(i) << ",";
          stress_tensor = AEL_data.stresscalculated.at(i);
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "stress_tensor = " << stress_tensor << endl;
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
        string ofileenergstrjson = AEL_data.dirpathname + "/AEL_energy_structures.json";
        if(!aurostd::stringstream2file(oss, ofileenergstrjson, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Unable to open file AEL_energy_structure.json" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());                       
      }

      // Clear vaspRuns vector to free up memory
      vaspRuns.clear();

      // Runs elastic constants algorithm within AFLOW
      aelerror = AEL_functions::elasticityfit(AEL_data, FileMESSAGE);
      if(aelerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Failed to calculate elastic constants" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aelerror;
      }

      // Checks that Poisson ratio is within the range -1.0 to 0.5
      // Materials with a Poisson ratio outside of this range are typically unstable
      if((AEL_data.poisson_ratio < -1.0) || (AEL_data.poisson_ratio > 0.5)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Poisson ratio = " << AEL_data.poisson_ratio << " is outside range from -1.0 to 0.5" << endl; ;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }

      // Checks that eigenvalues of elastic stiffness tensor are all greater than zero
      // Negative eigenvalues imply that the material is mechanically unstable
      for (uint i = 0; i < AEL_data.elastic_eigenvalues.size(); i++) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Elastic tensor eigenvalue " << i << " = " << AEL_data.elastic_eigenvalues.at(i) << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(AEL_data.elastic_eigenvalues.at(i) < 0.0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Elastic eigenvalue " << i << " = " << AEL_data.elastic_eigenvalues.at(i) << " is negative" << endl;
          aus << _AELSTR_WARNING_ + "Negative eigenvalue implies material is unstable" << endl;	  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
          AEL_data.mechanically_stable = false;
        }
      }

      // Uses Bravais lattice symmetry to check symmetry and consistency of elastic constant tensor
      if(USER_CHECK_ELASTIC_SYMMETRY.option) {
        aelerror = AEL_functions::elasticitycheck(AEL_data, xvasp, FileMESSAGE);
        if(aelerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Error checking elastic constants symmetry" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aelerror;
        }
      }      

      // Calculates elastic moduli using elastic constant tensor
      aelerror = AEL_functions::elasticmoduli(AEL_data, FileMESSAGE);
      if(aelerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Failed to calculate elastic moduli" << endl;  
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return aelerror;
      }

      // Compare bulk and shear moduli calculated using full strain set to those calculated using half strain set
      // If values differ by more than 15%, it is likely that the stress-strain curve is no longer linear
      double bulkvoigthalfratio = AEL_data.bulkmodulus_voigt_half / AEL_data.bulkmodulus_voigt;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Voigt bulk modulus = " << AEL_data.bulkmodulus_voigt << endl; 
      aus << _AELSTR_MESSAGE_ + "Voigt bulk modulus (half strain) = " << AEL_data.bulkmodulus_voigt_half << endl; 
      aus << _AELSTR_MESSAGE_ + "Voigt bulk modulus ratio = " << bulkvoigthalfratio << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
      if((bulkvoigthalfratio < 0.85) || (bulkvoigthalfratio > 1.15)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Voigt bulk modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
      }
      double bulkreusshalfratio = AEL_data.bulkmodulus_reuss_half / AEL_data.bulkmodulus_reuss;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Reuss bulk modulus = " << AEL_data.bulkmodulus_reuss << endl; 
      aus << _AELSTR_MESSAGE_ + "Reuss bulk modulus (half strain) = " << AEL_data.bulkmodulus_reuss_half << endl; 
      aus << _AELSTR_MESSAGE_ + "Reuss bulk modulus ratio = " << bulkreusshalfratio << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
      if((bulkreusshalfratio < 0.85) || (bulkreusshalfratio > 1.15)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Reuss bulk modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
      }
      double shearvoigthalfratio = AEL_data.shearmodulus_voigt_half / AEL_data.shearmodulus_voigt;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Voigt shear modulus = " << AEL_data.shearmodulus_voigt << endl; 
      aus << _AELSTR_MESSAGE_ + "Voigt shear modulus (half strain) = " << AEL_data.shearmodulus_voigt_half << endl; 
      aus << _AELSTR_MESSAGE_ + "Voigt shear modulus ratio = " << shearvoigthalfratio << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
      if((shearvoigthalfratio < 0.85) || (shearvoigthalfratio > 1.15)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Voigt shear modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
      }
      double shearreusshalfratio = AEL_data.shearmodulus_reuss_half / AEL_data.shearmodulus_reuss;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Reuss shear modulus = " << AEL_data.shearmodulus_reuss << endl; 
      aus << _AELSTR_MESSAGE_ + "Reuss shear modulus (half strain) = " << AEL_data.shearmodulus_reuss_half << endl; 
      aus << _AELSTR_MESSAGE_ + "Reuss shear modulus ratio = " << shearreusshalfratio << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
      if((shearreusshalfratio < 0.85) || (shearreusshalfratio > 1.15)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Reuss shear modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
      }
      double bulkvrhhalfratio = AEL_data.bulkmodulus_vrh_half / AEL_data.bulkmodulus_vrh;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "VRH bulk modulus = " << AEL_data.bulkmodulus_vrh << endl; 
      aus << _AELSTR_MESSAGE_ + "VRH bulk modulus (half strain) = " << AEL_data.bulkmodulus_vrh_half << endl; 
      aus << _AELSTR_MESSAGE_ + "VRH bulk modulus ratio = " << bulkvrhhalfratio << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
      if((bulkvrhhalfratio < 0.85) || (bulkvrhhalfratio > 1.15)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "VRH bulk modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
      }
      double shearvrhhalfratio = AEL_data.shearmodulus_vrh_half / AEL_data.shearmodulus_vrh;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "VRH shear modulus = " << AEL_data.shearmodulus_vrh << endl; 
      aus << _AELSTR_MESSAGE_ + "VRH shear modulus (half strain) = " << AEL_data.shearmodulus_vrh_half << endl; 
      aus << _AELSTR_MESSAGE_ + "VRH shear modulus ratio = " << shearvrhhalfratio << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
      if((shearvrhhalfratio < 0.85) || (shearvrhhalfratio > 1.15)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "VRH shear modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
      }

      // Extracts stress tensors for normal and shear strains at different pressures
      // Fits stress-strain data to obtain elastic constants at different pressures
      if(USER_PRESSURE_CALC.option) {
        for (uint k = 0; k < vaspRunsPressures.size(); k++) {
          aurostd::StringstreamClean(aus);
          // [OBSOLETE] if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "extract stress tensors for Pressure = " << Pressure.at(k+1) << "GPa" << endl;  
          aus << _AELSTR_MESSAGE_ + "extract stress tensors for Pressure = " << Pressure.at(k) << "GPa" << endl;  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // [OBSOLETE] }
          // Extract stress tensors from results of VASP calculations for strained structures
          aelerror = AEL_functions::extractstress(vaspRunsPressures.at(k), AEL_data_pressures.at(k), dirrunnamepressures.at(k), FileMESSAGE);
          if(aelerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_ERROR_ + "Failed to extract stress tensors" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aelerror;
          }
          // Check number of skipped runs in each direction does not exceed limit
          uint numskippedruns;
          for (uint i = 1; i <= AEL_data_pressures.at(k).normal_strain.size(); i++) {
            numskippedruns = AEL_data_pressures.at(k).normal_deformations.size() - AEL_data_pressures.at(k).normal_deformations_complete.at(i-1).size();
            if(numskippedruns > AEL_data_pressures.at(k).skiparunsmax) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_ERROR_ + "Number of failed runs in one independent directions = " << numskippedruns << endl;
              aus << _AELSTR_ERROR_ + "Number of failed runs exceeds limit in one independent direction of " << AEL_data.skiparunsmax << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }
          }
          for (uint i = 1; i <= AEL_data_pressures.at(k).shear_strain.size(); i++) {
            numskippedruns = AEL_data_pressures.at(k).shear_deformations.size() - AEL_data_pressures.at(k).shear_deformations_complete.at(i-1).size();
            if(numskippedruns > AEL_data_pressures.at(k).skiparunsmax) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_ERROR_ + "Number of failed runs in one independent directions = " << numskippedruns << endl;
              aus << _AELSTR_ERROR_ + "Number of failed runs exceeds limit per independent directions of " << AEL_data_pressures.at(k).skiparunsmax << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 1;
            }
          }

          // Clear vaspRuns vector to free up memory
          vaspRuns.clear();

          // Runs elastic constants algorithm within AFLOW
          aelerror = AEL_functions::elasticityfit(AEL_data_pressures.at(k), FileMESSAGE);
          if(aelerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_ERROR_ + "Failed to calculate elastic constants" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aelerror;
          }

          // Checks that Poisson ratio is within the range -1.0 to 0.5
          // Materials with a Poisson ratio outside of this range are typically unstable
          if((AEL_data_pressures.at(k).poisson_ratio < -1.0) || (AEL_data_pressures.at(k).poisson_ratio > 0.5)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Poisson ratio = " << AEL_data.poisson_ratio << " is outside range from -1.0 to 0.5" << endl; ;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }

          // Checks that eigenvalues of elastic stiffness tensor are all greater than zero
          // Negative eigenvalues imply that the material is mechanically unstable
          for (uint i = 0; i < AEL_data_pressures.at(k).elastic_eigenvalues.size(); i++) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Elastic tensor eigenvalue " << i << " = " << AEL_data_pressures.at(k).elastic_eigenvalues.at(i) << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(AEL_data_pressures.at(k).elastic_eigenvalues.at(i) < 0.0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "Elastic eigenvalue " << i << " = " << AEL_data_pressures.at(k).elastic_eigenvalues.at(i) << " is negative" << endl;
              aus << _AELSTR_WARNING_ + "Negative eigenvalue implies material is unstable" << endl;	  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
              AEL_data_pressures.at(k).mechanically_stable = false;
            }
          }

          // Uses Bravais lattice symmetry to check symmetry and consistency of elastic constant tensor
          if(USER_CHECK_ELASTIC_SYMMETRY.option) {
            aelerror = AEL_functions::elasticitycheck(AEL_data_pressures.at(k), xvasp, FileMESSAGE);
            if(aelerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_ERROR_ + "Error checking elastic constants symmetry" << endl;  
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aelerror;
            }
          }      

          // Calculates elastic moduli using elastic constant tensor
          aelerror = AEL_functions::elasticmoduli(AEL_data_pressures.at(k), FileMESSAGE);
          if(aelerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_ERROR_ + "Failed to calculate elastic moduli" << endl;  
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aelerror;
          }

          // Compare bulk and shear moduli calculated using full strain set to those calculated using half strain set
          // If values differ by more than 15%, it is likely that the stress-strain curve is no longer linear
          double bulkvoigthalfratio = AEL_data_pressures.at(k).bulkmodulus_voigt_half / AEL_data_pressures.at(k).bulkmodulus_voigt;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Voigt bulk modulus = " << AEL_data_pressures.at(k).bulkmodulus_voigt << endl; 
          aus << _AELSTR_MESSAGE_ + "Voigt bulk modulus (half strain) = " << AEL_data_pressures.at(k).bulkmodulus_voigt_half << endl; 
          aus << _AELSTR_MESSAGE_ + "Voigt bulk modulus ratio = " << bulkvoigthalfratio << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
          if((bulkvoigthalfratio < 0.85) || (bulkvoigthalfratio > 1.15)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Voigt bulk modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
          }
          double bulkreusshalfratio = AEL_data_pressures.at(k).bulkmodulus_reuss_half / AEL_data_pressures.at(k).bulkmodulus_reuss;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Reuss bulk modulus = " << AEL_data_pressures.at(k).bulkmodulus_reuss << endl; 
          aus << _AELSTR_MESSAGE_ + "Reuss bulk modulus (half strain) = " << AEL_data_pressures.at(k).bulkmodulus_reuss_half << endl; 
          aus << _AELSTR_MESSAGE_ + "Reuss bulk modulus ratio = " << bulkreusshalfratio << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
          if((bulkreusshalfratio < 0.85) || (bulkreusshalfratio > 1.15)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Reuss bulk modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
          }
          double shearvoigthalfratio = AEL_data_pressures.at(k).shearmodulus_voigt_half / AEL_data_pressures.at(k).shearmodulus_voigt;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Voigt shear modulus = " << AEL_data_pressures.at(k).shearmodulus_voigt << endl; 
          aus << _AELSTR_MESSAGE_ + "Voigt shear modulus (half strain) = " << AEL_data_pressures.at(k).shearmodulus_voigt_half << endl; 
          aus << _AELSTR_MESSAGE_ + "Voigt shear modulus ratio = " << shearvoigthalfratio << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
          if((shearvoigthalfratio < 0.85) || (shearvoigthalfratio > 1.15)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Voigt shear modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
          }
          double shearreusshalfratio = AEL_data_pressures.at(k).shearmodulus_reuss_half / AEL_data_pressures.at(k).shearmodulus_reuss;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Reuss shear modulus = " << AEL_data_pressures.at(k).shearmodulus_reuss << endl; 
          aus << _AELSTR_MESSAGE_ + "Reuss shear modulus (half strain) = " << AEL_data_pressures.at(k).shearmodulus_reuss_half << endl; 
          aus << _AELSTR_MESSAGE_ + "Reuss shear modulus ratio = " << shearreusshalfratio << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
          if((shearreusshalfratio < 0.85) || (shearreusshalfratio > 1.15)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Reuss shear modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
          }
          double bulkvrhhalfratio = AEL_data_pressures.at(k).bulkmodulus_vrh_half / AEL_data_pressures.at(k).bulkmodulus_vrh;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "VRH bulk modulus = " << AEL_data_pressures.at(k).bulkmodulus_vrh << endl; 
          aus << _AELSTR_MESSAGE_ + "VRH bulk modulus (half strain) = " << AEL_data_pressures.at(k).bulkmodulus_vrh_half << endl; 
          aus << _AELSTR_MESSAGE_ + "VRH bulk modulus ratio = " << bulkvrhhalfratio << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
          if((bulkvrhhalfratio < 0.85) || (bulkvrhhalfratio > 1.15)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "VRH bulk modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
          }
          double shearvrhhalfratio = AEL_data_pressures.at(k).shearmodulus_vrh_half / AEL_data_pressures.at(k).shearmodulus_vrh;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "VRH shear modulus = " << AEL_data_pressures.at(k).shearmodulus_vrh << endl; 
          aus << _AELSTR_MESSAGE_ + "VRH shear modulus (half strain) = " << AEL_data_pressures.at(k).shearmodulus_vrh_half << endl; 
          aus << _AELSTR_MESSAGE_ + "VRH shear modulus ratio = " << shearvrhhalfratio << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 
          if((shearvrhhalfratio < 0.85) || (shearvrhhalfratio > 1.15)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "VRH shear modulus difference exceeds 15%: stress-strain curve might be non-linear" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	
          }
        }
      }

      // Checks if results at a pressure where the material is stable should be written as main results
      // If yes, finds lowest absolute value pressure where material is stable
      // Save a copy of the original AEL_data class first, just in case
      AEL_data_orig = AEL_data;
      if(USER_PRESSURE_CALC.option && USER_USE_STABLE_PRESSURE.option) {
        // Find highest pressure value
        double maxPressure = 0.0;
        for (uint k = 0; k < Pressure.size(); k++) {
          if (std::abs(Pressure.at(k)) > maxPressure) {
            maxPressure = Pressure.at(k);
          }
        }	
        double lowestStablePressure = maxPressure + 1.0;
        stablePressureStructure = 0;
        for (uint k = 0; k < AEL_data_pressures.size(); k++) {
          if (AEL_data_pressures.at(k).mechanically_stable) {
            stablePressureFound = true;
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Structure at pressure = " << AEL_data_pressures.at(k).applied_pressure << "GPa is elastically stable" << endl; 
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	    
            if (std::abs(AEL_data_pressures.at(k).applied_pressure) < lowestStablePressure) {
              stablePressureStructure = k;
              lowestStablePressure = AEL_data_pressures.at(k).applied_pressure;
            }
          }
        }
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lowest pressure with stable structure = " << lowestStablePressure << "GPa" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	    
        if (stablePressureFound) {
          AEL_data = AEL_data_pressures.at(stablePressureStructure);
        }
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Pressure for structure used for results  = " << AEL_data.applied_pressure << "GPa" << endl; 
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	 	    
      }

      // Writes elastic moduli values in file for REST-API
      oss.clear();
      oss.str(std::string());
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << "[AEL_RESULTS]START" << endl;
      oss << "ael_poisson_ratio=" << AEL_data.poisson_ratio << "  " << endl;
      oss << "ael_bulk_modulus_voigt=" << AEL_data.bulkmodulus_voigt << "  (GPa)" << endl;
      oss << "ael_bulk_modulus_reuss=" << AEL_data.bulkmodulus_reuss  << "  (GPa)" << endl;
      oss << "ael_shear_modulus_voigt=" << AEL_data.shearmodulus_voigt << "  (GPa)" << endl;
      oss << "ael_shear_modulus_reuss=" << AEL_data.shearmodulus_reuss << "  (GPa)" << endl;
      oss << "ael_bulk_modulus_vrh=" << AEL_data.bulkmodulus_vrh << "  (GPa)" << endl;
      oss << "ael_shear_modulus_vrh=" << AEL_data.shearmodulus_vrh << "  (GPa)" << endl;
      oss << "ael_elastic_anisotropy=" << AEL_data.elastic_anisotropy << "  " << endl;
      oss << "ael_youngs_modulus_vrh=" << AEL_data.youngsmodulus_vrh << "  (GPa)" << endl;
      oss << "ael_speed_sound_transverse=" << AEL_data.speed_sound_transverse << " (m/s)" << endl;
      oss << "ael_speed_sound_longitudinal=" << AEL_data.speed_sound_longitudinal << " (m/s)" << endl;
      oss << "ael_speed_sound_average=" << AEL_data.speed_sound_average << " (m/s)" << endl;
      oss << "ael_pughs_modulus_ratio=" << AEL_data.pughs_modulus_ratio << "  " << endl;
      oss << "ael_debye_temperature=" << AEL_data.debye_temperature << "  (K)" << endl;
      oss << "ael_applied_pressure=" << AEL_data.applied_pressure << "  (GPa)" << endl;
      oss << "ael_average_external_pressure=" << AEL_data.average_external_pressure << "  (GPa)" << endl;   
      oss << "[AEL_RESULTS]STOP" << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << "[AEL_STIFFNESS_TENSOR]START" << endl;
      for (uint i = 0; i < AEL_data.elastic_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.elastic_tensor.size(); j++) {
          oss << setw(15) << setprecision(6) << (AEL_data.elastic_tensor.at(i).at(j));
        }
        oss << endl;
      }
      oss << "[AEL_STIFFNESS_TENSOR]STOP" << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
      oss << "[AEL_COMPLIANCE_TENSOR]START" << endl;
      for (uint i = 0; i < AEL_data.compliance_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.compliance_tensor.size(); j++) {
          oss << setw(15) << setprecision(6) << (AEL_data.compliance_tensor.at(i).at(j));
        }
        oss << endl;
      }
      oss << "[AEL_COMPLIANCE_TENSOR]STOP" << endl;
      oss << AFLOWIN_SEPARATION_LINE << endl;
      string ofileafaelname = AEL_data.dirpathname + "/aflow.ael.out";
      if(!aurostd::stringstream2file(oss, ofileafaelname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Unable to open file aflow.ael.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      oss.clear();
      oss.str(std::string());

      // Elastic tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Opening file AEL_Elastic_constants.out" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss << "# Elastic constants: stiffness tensor" << endl;
      for (uint i = 0; i < AEL_data.elastic_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.elastic_tensor.size(); j++) {
          oss << setw(15) << setprecision(6) << (AEL_data.elastic_tensor.at(i).at(j));
        }
        oss << endl;
      }
      string ofileelconstname = AEL_data.dirpathname + "/AEL_Elastic_constants.out";
      if(!aurostd::stringstream2file(oss, ofileelconstname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Unable to open file AEL_Elastic_constants.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss.clear();
      oss.str(std::string());

      // Compliance tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Opening file AEL_Compliance_tensor.out" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss << "# Compliance tensor" << endl;
      for (uint i = 0; i < AEL_data.compliance_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.compliance_tensor.size(); j++) {
          oss << setw(15) << setprecision(6) << (AEL_data.compliance_tensor.at(i).at(j));
        }
        oss << endl;
      }
      string ofilecomptensname = AEL_data.dirpathname + "/AEL_Compliance_tensor.out";
      if(!aurostd::stringstream2file(oss, ofilecomptensname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Unable to open file AEL_Compliance_tensor.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      // [OBSOLETE] aurostd::StringstreamClean(oss);
      oss.clear();
      oss.str(std::string());

      // Write calculated elastic tensors to JSON format:
      // Elastic stiffness tensor
      // Elastic compliance tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "Opening file AEL_elastic_tensor.json" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      uint elasticindexi, elasticindexj;
      string strelasticindexi, strelasticindexj;
      string elastic_constant_label;
      ostringstream cnvstrelasticindexi, cnvstrelasticindexj;
      // [OBSOLETE] oss << "{\"elastic_stiffness_tensor\": [";
      oss << "{\"elastic_stiffness_tensor\": {";
      for (uint i = 0; i < AEL_data.elastic_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.elastic_tensor.size(); j++) {
          elasticindexi = i + 1;
          elasticindexj = j + 1;
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index i = " << elasticindexi <<  endl;
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index j = " << elasticindexj <<  endl;	  
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
          cnvstrelasticindexi.clear();
          cnvstrelasticindexi.str(std::string());
          cnvstrelasticindexj.clear();
          cnvstrelasticindexj.str(std::string());
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index i ss = " << cnvstrelasticindexi.str() <<  endl;
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index j ss = " << cnvstrelasticindexj.str() <<  endl;	  
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
          cnvstrelasticindexi << elasticindexi;
          cnvstrelasticindexj << elasticindexj;
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index i ss = " << cnvstrelasticindexi.str() <<  endl;
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index j ss = " << cnvstrelasticindexj.str() <<  endl;	  
          // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
          strelasticindexi = cnvstrelasticindexi.str();
          strelasticindexj = cnvstrelasticindexj.str();
          // [OBSOLETE] aurostd::StringstreamClean(aus);
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index i str = " << strelasticindexi <<  endl;
          // [OBSOLETE] aus << _AELSTR_MESSAGE_ << "Elastic index j str = " << strelasticindexj <<  endl;	  
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
          elastic_constant_label = "c_" + strelasticindexi + strelasticindexj;
          if ((i==0) && (j==0)) {
            // [OBSOLETE] oss << "{\"" + elastic_constant_label + "\":" << AEL_data.elastic_tensor.at(i).at(j) << "}";
            oss << "\"" + elastic_constant_label + "\":" << AEL_data.elastic_tensor.at(i).at(j);
          } else {
            // [OBSOLETE] oss << ",{\"" + elastic_constant_label + "\":" << AEL_data.elastic_tensor.at(i).at(j) << "}";
            oss << ",\"" + elastic_constant_label + "\":" << AEL_data.elastic_tensor.at(i).at(j);
          }
        }
      }
      // [OBSOLETE] oss << "],\"elastic_compliance_tensor\": [";
      oss << "},\"elastic_compliance_tensor\": {";
      for (uint i = 0; i < AEL_data.compliance_tensor.size(); i++) {
        for (uint j = 0; j < AEL_data.compliance_tensor.size(); j++) {
          elasticindexi = i + 1;
          elasticindexj = j + 1;
          cnvstrelasticindexi.clear();
          cnvstrelasticindexi.str(std::string());
          cnvstrelasticindexj.clear();
          cnvstrelasticindexj.str(std::string());
          cnvstrelasticindexi << elasticindexi;
          cnvstrelasticindexj << elasticindexj;
          strelasticindexi = cnvstrelasticindexi.str();
          strelasticindexj = cnvstrelasticindexj.str();
          elastic_constant_label = "s_" + strelasticindexi + strelasticindexj;
          if ((i==0) && (j==0)) {
            // [OBSOLETE] oss << "{\"" + elastic_constant_label + "\":" << AEL_data.compliance_tensor.at(i).at(j) << "}";
            oss << "\"" + elastic_constant_label + "\":" << AEL_data.compliance_tensor.at(i).at(j);
          } else {
            // [OBSOLETE] oss << ",{\"" + elastic_constant_label + "\":" << AEL_data.compliance_tensor.at(i).at(j) << "}";
            oss << ",\"" + elastic_constant_label + "\":" << AEL_data.compliance_tensor.at(i).at(j);
          }
        }
      }
      // [OBSOLETE] oss << "],\"units\":\"GPa\"}" << endl;
      oss << "},\"units\":\"GPa\"}" << endl;
      string ofileelastictensorjson = AEL_data.dirpathname + "/AEL_elastic_tensor.json";
      if(!aurostd::stringstream2file(oss, ofileelastictensorjson, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Unable to open file AEL_elastic_tensor.json" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      oss.clear();
      oss.str(std::string());

      // Write out values for full elastic constants and compliance tensors
      if(USER_WRITE_FULL_RESULTS.option) {
        // Writes elastic moduli values in file in parseable format
        oss << "ael_poisson_ratio=" << AEL_data.poisson_ratio << " | ";
        oss << "ael_bulk_modulus_voigt=" << AEL_data.bulkmodulus_voigt << " | ";
        oss << "ael_bulk_modulus_reuss=" << AEL_data.bulkmodulus_reuss << " | ";
        oss << "ael_shear_modulus_voigt=" << AEL_data.shearmodulus_voigt << " | ";
        oss << "ael_shear_modulus_reuss=" << AEL_data.shearmodulus_reuss << " | ";
        oss << "ael_bulk_modulus_vrh=" << AEL_data.bulkmodulus_vrh << " | ";
        oss << "ael_shear_modulus_vrh=" << AEL_data.shearmodulus_vrh << " | ";
        // [OBSOLETE] oss << "ael_elastic_anisotropy=" << AEL_data.elastic_anisotropy << " | ";
        oss << "ael_elastic_anisotropy=" << AEL_data.elastic_anisotropy;
        string ofileaellibname = AEL_data.dirpathname + "/aellib.out";
        if(!aurostd::stringstream2file(oss, ofileaellibname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Unable to open file aellib.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());

        // Directional elastic moduli
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Opening file AEL_elastic_moduli_directional.out" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss << "# Elastic moduli and Poisson ratio for different stress/strain directions" << endl;
        oss << "# Valid for orthotropic materials" << endl;
        oss << "#########################" << endl;
        oss << "# Young's modulus X Y Z #" << endl;
        for (uint i = 0; i < AEL_data.youngsmodulus_directional.size(); i++) {
          oss << setw(15) << setprecision(6) << (AEL_data.youngsmodulus_directional.at(i));
        }
        oss << endl;
        oss << "##########################" << endl;
        oss << "# Shear modulus YZ XZ XY #" << endl;
        for (uint i = 0; i < AEL_data.shearmodulus_directional.size(); i++) {
          oss << setw(15) << setprecision(6) << (AEL_data.shearmodulus_directional.at(i));
        }
        oss << endl;
        oss << "###################################" << endl;
        oss << "# Poisson ratio XY XZ YX YZ ZX ZY #" << endl;
        for (uint i = 0; i < AEL_data.poissonratio_directional.size(); i++) {
          oss << setw(15) << setprecision(6) << (AEL_data.poissonratio_directional.at(i));
        }
        oss << endl;
        string ofileelmodirname = AEL_data.dirpathname + "/AEL_elastic_moduli_directional.out";
        if(!aurostd::stringstream2file(oss, ofileelmodirname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Unable to open file AEL_elastic_moduli_directional.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        // [OBSOLETE] aurostd::StringstreamClean(oss);
        oss.clear();
        oss.str(std::string());
      }

      // Write out elastic properties for different pressures
      if(USER_PRESSURE_CALC.option) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Writing elastic constants as a function of pressure" << endl;
        aus << _AELSTR_MESSAGE_ + "Number of pressure values = " << vaspRunsPressures.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        oss.clear();
        oss.str(std::string());
        for (uint k = 0; k < vaspRunsPressures.size(); k++) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "k = " << k << endl;
          aus << _AELSTR_MESSAGE_ + "Pressure = " << Pressure.at(k) << "GPa" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // Writes elastic moduli values in file for REST-API
          oss << AFLOWIN_SEPARATION_LINE << endl;
          oss << "[AEL_PRESSURE]PRESSURE = " << Pressure.at(k) << "GPa" << endl;
          if(AEL_data_pressures.at(k).mechanically_stable) {
            oss << "[AEL_PRESSURE]Structure is mechanically stable at this pressure" << endl;
          } else {
            oss << "[AEL_PRESSURE]Structure is not mechanically stable at this pressure" << endl;
          }
          oss << "[AEL_RESULTS]START" << endl;
          oss << "ael_poisson_ratio=" << AEL_data_pressures.at(k).poisson_ratio << "  " << endl;
          oss << "ael_bulk_modulus_voigt=" << AEL_data_pressures.at(k).bulkmodulus_voigt << "  (GPa)" << endl;
          oss << "ael_bulk_modulus_reuss=" << AEL_data_pressures.at(k).bulkmodulus_reuss  << "  (GPa)" << endl;
          oss << "ael_shear_modulus_voigt=" << AEL_data_pressures.at(k).shearmodulus_voigt << "  (GPa)" << endl;
          oss << "ael_shear_modulus_reuss=" << AEL_data_pressures.at(k).shearmodulus_reuss << "  (GPa)" << endl;
          oss << "ael_bulk_modulus_vrh=" << AEL_data_pressures.at(k).bulkmodulus_vrh << "  (GPa)" << endl;
          oss << "ael_shear_modulus_vrh=" << AEL_data_pressures.at(k).shearmodulus_vrh << "  (GPa)" << endl;
          oss << "ael_elastic_anisotropy=" << AEL_data_pressures.at(k).elastic_anisotropy << "  " << endl;
          oss << "ael_speed_sound_transverse=" << AEL_data_pressures.at(k).speed_sound_transverse << " (m/s)" << endl;
          oss << "ael_speed_sound_longitudinal=" << AEL_data_pressures.at(k).speed_sound_longitudinal << " (m/s)" << endl;
          oss << "ael_speed_sound_average=" << AEL_data_pressures.at(k).speed_sound_average << " (m/s)" << endl;
          // [OBSOLETE] oss << "ael_pughs_modulus_ratio=" << AEL_data_pressures.at(k).pughs_modulus_ratio << "  " << endl;
          oss << "ael_debye_temperature=" << AEL_data_pressures.at(k).debye_temperature << "  " << endl;
          oss << "ael_applied_pressure=" << AEL_data_pressures.at(k).applied_pressure << "  (GPa)" << endl;
          oss << "ael_average_external_pressure=" << AEL_data_pressures.at(k).average_external_pressure << "  (GPa)" << endl;  
          oss << "[AEL_RESULTS]STOP" << endl;
          oss << AFLOWIN_SEPARATION_LINE << endl;
          oss << "[AEL_STIFFNESS_TENSOR]START" << endl;
          for (uint i = 0; i < AEL_data_pressures.at(k).elastic_tensor.size(); i++) {
            for (uint j = 0; j < AEL_data_pressures.at(k).elastic_tensor.size(); j++) {
              oss << setw(15) << setprecision(6) << (AEL_data_pressures.at(k).elastic_tensor.at(i).at(j));
            }
            oss << endl;
          }
          oss << "[AEL_STIFFNESS_TENSOR]STOP" << endl;
          oss << AFLOWIN_SEPARATION_LINE << endl;
          oss << "[AEL_COMPLIANCE_TENSOR]START" << endl;
          for (uint i = 0; i < AEL_data_pressures.at(k).compliance_tensor.size(); i++) {
            for (uint j = 0; j < AEL_data_pressures.at(k).compliance_tensor.size(); j++) {
              oss << setw(15) << setprecision(6) << (AEL_data_pressures.at(k).compliance_tensor.at(i).at(j));
            }
            oss << endl;
          }
          oss << "[AEL_COMPLIANCE_TENSOR]STOP" << endl;
          oss << AFLOWIN_SEPARATION_LINE << endl;
        }
        string ofileafaelname = AEL_data.dirpathname + "/aflow.ael_pressure.out";
        if(!aurostd::stringstream2file(oss, ofileafaelname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Unable to open file aflow.ael_pressure.out" <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }	
        oss.clear();
        oss.str(std::string());
      }
    }
    catch (AELStageBreak& e) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_NOTICE_ + "Stopped. Waiting for completion of required calculations..." << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aelerror = 8;
    }

    return aelerror;
  }
} // namespace AEL_functions

// **************************************************************************
//  End of AFLOW AEL
// **************************************************************************


#endif  // _AFLOW_AEL_ELASTICITY_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
