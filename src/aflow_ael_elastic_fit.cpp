// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AEL_ELASTIC_FIT_CPP
#define _AFLOW_AEL_ELASTIC_FIT_CPP
#include "aflow.h"
#include "aflow_ael_elasticity.h"
#include "aflow_agl_debye.h"

// ###############################################################################
//                  AFLOW Automatic Elasticity Library (AEL) (2014-2019)
// ###############################################################################
//
// Uses strain-stress calculations to obtain elastic constants of materials
// Based on original Python program written by M. de Jong et al.
// See Scientific Data 2, 150009 (2015) for details of original program
// See Phys. Rev. Materials 1, 015401 (2017) for details of this implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AEL
//
// ************************************************************************************************
// This set of functions calculate elastic constants and elastic moduli from stress-strain data
// ************************************************************************************************

// ***************************************************************************
// AEL_functions::elasticityfit
// ***************************************************************************
namespace AEL_functions {
  //
  // elasticityfit: Fit stress-strain data to obtain elastic constants
  // See Phys. Rev. Materials 1, 015401 (2017) and Scientific Data 2, 150009 (2015) for details
  //
  uint elasticityfit(_AEL_data& AEL_data, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Running elasticity fit" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    xmatrix<double> strain_components(3, 3);
    xmatrix<double> stress_tensor(3, 3);
    xmatrix<double> str_identity(3, 3);
    xmatrix<double> elastic_const(6, 6);
    xmatrix<double> elastic_const_half(6, 6);
    xmatrix<double> elastic_const_inv(6, 6);
    xmatrix<double> elastic_const_half_inv(6, 6);
    xvector<double> elastic_eigenvals(6);
    xvector<double> elastic_eigenval_real(6);
    xvector<double> elastic_eigenval_imag(6);
    xmatrix<double> zero_stress(3, 3);
    vector<vector<double> > stress_components;
    vector<vector<double> > stress_components_half;
    vector<double> shear_deformations_doubled;
    vector<double> normal_deformations_origin;
    vector<double> shear_deformations_origin;
    vector<vector<double> > normal_deformations_tofit;
    vector<vector<double> > shear_deformations_tofit;
    vector<vector<double> > normal_deformations_half_tofit;
    vector<vector<double> > shear_deformations_half_tofit;
    vector<vector<xmatrix<double> > > normal_stress_tofit;
    vector<vector<xmatrix<double> > > shear_stress_tofit;
    vector<vector<xmatrix<double> > > normal_stress_half_tofit;
    vector<vector<xmatrix<double> > > shear_stress_half_tofit;
    double Cij = 0.0;
    double Cijmax = 0.0;
    double Cijtol = 0.0;
    double Cijhalf = 0.0;
    double tol = 1e-12;
    uint aelerror = 0;
    // OBSOLETE uint halffirstpointnormal = 0;
    // OBSOLETE uint halffirstpointshear = 0;
    // OBSOLETE uint halfnumberruns = 0;
    vector<uint> halfstartpointnormal;
    vector<uint> halfstartpointshear;
    uint halfstartpoint = 0;
    uint midpointnormal = 0;
    uint midpointshear = 0;
    int imidpointnormal = 0;
    int imidpointshear = 0;
    uint mid_midpointnormal = 0;
    uint mid_midpointshear = 0;
    uint firsthalfcompletenormal = 0;
    uint firsthalfcompleteshear = 0;
    uint firsthalfincompletenormal = 0;
    uint firsthalfincompleteshear = 0;
    int available_points = 0;
    int halfstartpointshift = 0;
    // OBSOLETE int nlinfit = 2;
    int npolyfitorderhalf;
    int npolyfitorder;
    uint nupolyfitorder;
    midpointnormal = AEL_data.normal_deformations.size() / 2;
    midpointshear = AEL_data.shear_deformations.size() / 2;
    mid_midpointnormal = AEL_data.normal_deformations.size() / 4;
    mid_midpointshear = AEL_data.shear_deformations.size() / 4;
    // OBSOLETE if(AEL_data.negstrain) {
    // OBSOLETE   halffirstpointnormal = AEL_data.normal_deformations.size() / 4;
    // OBSOLETE } 
    // OBSOLETE if(AEL_data.negstrain) {
    // OBSOLETE   halffirstpointshear = AEL_data.shear_deformations.size() / 4;
    // OBSOLETE }
    // OBSOLETE halfnumberruns = AEL_data.normal_deformations.size() / 2;
    normal_deformations_tofit.resize(AEL_data.normal_strain.size());
    shear_deformations_tofit.resize(AEL_data.shear_strain.size());
    normal_deformations_half_tofit.resize(AEL_data.normal_strain.size());
    shear_deformations_half_tofit.resize(AEL_data.shear_strain.size());
    normal_stress_tofit.resize(AEL_data.normal_strain.size());
    shear_stress_tofit.resize(AEL_data.shear_strain.size());
    normal_stress_half_tofit.resize(AEL_data.normal_strain.size());
    shear_stress_half_tofit.resize(AEL_data.shear_strain.size());
    // Initialize "zero stress" tensor
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Initializing zero stress tensor" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        zero_stress[i][j] = 0.0;
      }
    }
    for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
      if(AEL_data.normal_deformations_complete.at(i-1).size() < midpointnormal) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Not enough complete runs to fit half number of points" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
    }
    for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
      if(AEL_data.shear_deformations_complete.at(i-1).size() < midpointnormal) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_ERROR_ + "Not enough complete runs to fit half number of points" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
    }
    // Set up normal deformations to be fitted
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Setting up normal deformations to be fitted" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
      // Adds completed deformation runs to set of deformations to fit for each strain direction
      for (uint j = 0; j < AEL_data.normal_deformations_complete.at(i-1).size(); j++) {   
        normal_deformations_tofit.at(i-1).push_back(AEL_data.normal_deformations_complete.at(i-1).at(j));
        normal_stress_tofit.at(i-1).push_back(AEL_data.normal_stress.at(i-1).at(j));
      }
      // Checks if number of complete deformation runs is less than number of initial runs
      if(AEL_data.normal_deformations_complete.at(i-1).size() < AEL_data.normal_deformations.size()) {
        firsthalfcompletenormal = 0;
        firsthalfincompletenormal = 0;
        // Counts number of runs in first half of deformation set which were completed
        for (uint j = 0; j < midpointnormal; j++) {
          for (uint ij = 0; ij < AEL_data.normal_deformations_complete.at(i-1).size(); ij++) {
            if((AEL_data.normal_deformations_complete.at(i-1).at(ij) - AEL_data.normal_deformations.at(j)) < tol) {
              // OBSOLETE halfstartpointnormal = halfstartpointnormal + 1;
              firsthalfcompletenormal = firsthalfcompletenormal + 1;
            }
          }
        }
        // Calculates number of runs in first half of deformation set which were skipped
        firsthalfincompletenormal = midpointnormal - firsthalfcompletenormal;
        // Moves initial point for half runs fit back to account for skipped runs
        halfstartpointshift = mid_midpointnormal - firsthalfincompletenormal;
        // Makes sure that there are still sufficient remaining runs for half fitting
        // OBSOLETE uint jmax = AEL_data.normal_deformations_complete.at(i-1).size() - 1;
        // OBSOLETE while ((AEL_data.normal_deformations_complete.at(i-1).at(jmax) - halfstartpointnormal) < halfnumberruns) {
        // OBSOLETE   if(halfstartpointnormal < 1) {
        // OBSOLETE     aurostd::StringstreamClean(aus);
        // OBSOLETE     aus << _AELSTR_ERROR_ + "Insufficient successful runs to fit half of initial points" << endl;
        // OBSOLETE     aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
        // OBSOLETE     return 1;
        // OBSOLETE   }
        // OBSOLETE   halfstartpointnormal = halfstartpointnormal - 1;
        // OBSOLETE }
        available_points = AEL_data.normal_deformations_complete.at(i-1).size() - halfstartpointshift;
        imidpointnormal = midpointnormal;
        if(available_points < imidpointnormal) {
          halfstartpointshift = AEL_data.normal_deformations_complete.at(i-1).size() - midpointnormal;
        }
        if(halfstartpointshift < 0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Not enough complete runs to fit half number of points" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        } else {
          halfstartpoint = halfstartpointshift;
        }
        if(AEL_data.negstrain) {
          halfstartpointnormal.push_back(halfstartpoint);
        } else {
          halfstartpoint = 0;
          halfstartpointnormal.push_back(halfstartpoint);
        }
        // Appends half of runs to initial
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "halfstartpoint = " << halfstartpoint << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        for (uint j = halfstartpoint; j < (halfstartpoint + midpointnormal); j++) {
          normal_deformations_half_tofit.at(i-1).push_back(AEL_data.normal_deformations_complete.at(i-1).at(j));
          normal_stress_half_tofit.at(i-1).push_back(AEL_data.normal_stress.at(i-1).at(j));
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "j = " << j << endl;
          aus << _AELSTR_MESSAGE_ + "normal_deformations_half_tofit.at(i-1).at(j - halfstartpoint) = " << normal_deformations_half_tofit.at(i-1).at(j - halfstartpoint) << endl;
          aus << _AELSTR_MESSAGE_ + "normal_stress_half_tofit.at(i-1).at(j - halfstartpoint) = " << normal_stress_half_tofit.at(i-1).at(j - halfstartpoint) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
        }
      } else {
        if(AEL_data.negstrain) {
          halfstartpoint = mid_midpointnormal;
          halfstartpointnormal.push_back(halfstartpoint);
        } else {
          halfstartpoint = 0;
          halfstartpointnormal.push_back(halfstartpoint);
        }
        // If all runs are successful, appends original set of half runs
        for (uint j = 0; j < AEL_data.normal_deformations_half.size(); j++) {
          normal_deformations_half_tofit.at(i-1).push_back(AEL_data.normal_deformations_half.at(j));
          normal_stress_half_tofit.at(i-1).push_back(AEL_data.normal_stress.at(i-1).at(j + halfstartpoint));	  
        }
      }
    }
    // Set up shear deformations to be fitted
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Setting up shear deformations to be fitted" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // OBSOLETE for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
    // OBSOLETE   for (uint j = 0; j < AEL_data.shear_deformations_complete.at(i-1).size(); j++) {   
    // OBSOLETE 	shear_deformations_tofit.at(i-1).push_back(AEL_data.shear_deformations_complete.at(i-1).at(j));
    // OBSOLETE 	shear_stress_tofit.at(i-1).push_back(AEL_data.shear_stress.at(i-1).at(j));
    // OBSOLETE   }
    // OBSOLETE }
    for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
      // Adds completed deformation runs to set of deformations to fit for each strain direction
      for (uint j = 0; j < AEL_data.shear_deformations_complete.at(i-1).size(); j++) {   
        shear_deformations_tofit.at(i-1).push_back(AEL_data.shear_deformations_complete.at(i-1).at(j));
        shear_stress_tofit.at(i-1).push_back(AEL_data.shear_stress.at(i-1).at(j));
      }
      // Checks if number of complete deformation runs is less than number of initial runs
      if(AEL_data.shear_deformations_complete.at(i-1).size() < AEL_data.shear_deformations.size()) {
        firsthalfcompleteshear = 0;
        firsthalfincompleteshear = 0;
        // Counts number of runs in first half of deformation set which were completed
        for (uint j = 0; j < midpointshear; j++) {
          for (uint ij = 0; ij < AEL_data.shear_deformations_complete.at(i-1).size(); ij++) {
            if((AEL_data.shear_deformations_complete.at(i-1).at(ij) - AEL_data.shear_deformations.at(j)) < tol) {
              // OBSOLETE halfstartpointshear = halfstartpointshear + 1;
              firsthalfcompleteshear = firsthalfcompleteshear + 1;
            }
          }
        }
        // Calculates number of runs in first half of deformation set which were skipped
        firsthalfincompleteshear = midpointshear - firsthalfcompleteshear;
        // Moves initial point for half runs fit back to account for skipped runs
        halfstartpointshift = mid_midpointshear - firsthalfincompleteshear;
        // Makes sure that there are still sufficient remaining runs for half fitting
        // OBSOLETE uint jmax = AEL_data.shear_deformations_complete.at(i-1).size() - 1;
        // OBSOLETE while (AEL_data.shear_deformations_complete.at(i-1).at(jmax) < (halfstartpointshear + halfnumberruns)) {
        // OBSOLETE   if(halfstartpointshear < 1) {
        // OBSOLETE     aurostd::StringstreamClean(aus);
        // OBSOLETE     aus << _AELSTR_ERROR_ + "Insufficient successful runs to fit half of initial points" << endl;
        // OBSOLETE     aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
        // OBSOLETE     return 1;
        // OBSOLETE   }
        // OBSOLETE   halfstartpointshear = halfstartpointshear - 1;
        // OBSOLETE }
        available_points = AEL_data.shear_deformations_complete.at(i-1).size() - halfstartpointshift;
        imidpointshear = midpointshear;
        if(available_points < imidpointshear) {
          halfstartpointshift = AEL_data.shear_deformations_complete.at(i-1).size() - midpointshear;
        }
        if(halfstartpointshift < 0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_ERROR_ + "Not enough complete runs to fit half number of points" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        } else {
          halfstartpoint = halfstartpointshift;
        }
        if(AEL_data.negstrain) {
          halfstartpointshear.push_back(halfstartpoint);
        } else {
          halfstartpoint = 0;
          halfstartpointshear.push_back(halfstartpoint);
        }
        // Appends half of runs to initial
        for (uint j = halfstartpoint; j < (halfstartpoint + midpointshear); j++) {
          shear_deformations_half_tofit.at(i-1).push_back(AEL_data.shear_deformations_complete.at(i-1).at(j));
          shear_stress_half_tofit.at(i-1).push_back(AEL_data.shear_stress.at(i-1).at(j));	  
        }
      } else {
        if(AEL_data.negstrain) {
          halfstartpoint = mid_midpointshear;
          halfstartpointshear.push_back(halfstartpoint);
        } else {
          halfstartpoint = 0;
          halfstartpointshear.push_back(halfstartpoint);
        }
        // If all runs are successful, appends original set of half runs
        for (uint j = 0; j < AEL_data.shear_deformations_half.size(); j++) {
          shear_deformations_half_tofit.at(i-1).push_back(AEL_data.shear_deformations_half.at(j));
          shear_stress_half_tofit.at(i-1).push_back(AEL_data.shear_stress.at(i-1).at(j + halfstartpoint));	  
        }
      }
    }
    stress_components.resize(6);
    // OBSOLETE for (uint k = 0; k < stress_components.size(); k++) {
    // OBSOLETE   stress_components.at(k).resize(AEL_data.normal_deformations.size());
    // OBSOLETE }
    stress_components_half.resize(6);
    // OBSOLETE for (uint k = 0; k < stress_components_half.size(); k++) {
    // OBSOLETE   stress_components_half.at(k).resize(AEL_data.normal_deformations_half.size());
    // OBSOLETE }
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Checking origin fit" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(AEL_data.fitstrainorigin || AEL_data.fitrelaxedstruct) {
      for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
        normal_deformations_tofit.at(i-1).push_back(0.0);	
        if(AEL_data.calcstrainorigin || AEL_data.fitrelaxedstruct) {
          normal_stress_tofit.at(i-1).push_back(AEL_data.origin_stress.at(0));
        } else {
          normal_stress_tofit.at(i-1).push_back(zero_stress);
        }
      }
      // OBSOLETE for (uint k = 0; k < stress_components.size(); k++) {
      // OBSOLETE   stress_components.at(k).resize(normal_deformations_origin.size());
      // OBSOLETE }
      for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
        shear_deformations_tofit.at(i-1).push_back(0.0);
        if(AEL_data.calcstrainorigin || AEL_data.fitrelaxedstruct) {
          shear_stress_tofit.at(i-1).push_back(AEL_data.origin_stress.at(0));
        } else {
          shear_stress_tofit.at(i-1).push_back(zero_stress);
        }
      }
      // OBSOLETE shear_deformations_origin.push_back(0.0);
    }
    aurostd::identity(str_identity);
    // Fit normal stress-strain
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Fitting normal stress-strain data" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "i = " << i << ", AEL_data.normal_strain.size() = " << AEL_data.normal_strain.size() << ", normal_deformations_tofit.size() = " << normal_deformations_tofit.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (uint k = 0; k < stress_components.size(); k++) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "k = " << k << ", stress_components.size() = " << stress_components.size() << ", normal_deformations_tofit.at(i-1).size() = " << normal_deformations_tofit.at(i-1).size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        stress_components.at(k).resize(normal_deformations_tofit.at(i-1).size());
      }
      for (uint j = 0; j < normal_deformations_tofit.at(i-1).size(); j++) {   
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "i = " << i << ", j = " << j << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Normal strain tensor = " << endl;
          // OBSOLETE aus << AEL_data.normal_strain.at(i-1).at(j) << endl;
          aus << _AELSTR_MESSAGE_ + "Normal strain deformations = " << endl;
          aus << normal_deformations_tofit.at(i-1).at(j) << endl;
          aus << _AELSTR_MESSAGE_ + "Normal stress tensor = " << endl; 
          aus << AEL_data.normal_stress.at(i-1).at(j) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // OBSOLETE strain_components = AEL_data.normal_strain.at(i-1).at(j) - str_identity;
        // OBSOLETE stress_tensor = AEL_data.normal_stress.at(i-1).at(j);
        stress_tensor = normal_stress_tofit.at(i-1).at(j);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Strain tensor = " << endl;
          // OBSOLETE aus << strain_components << endl;
          aus << _AELSTR_MESSAGE_ + "Stress tensor = " << endl; 
          aus << stress_tensor << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Stress components for stress-strain fitting; see Eq. 1, M. de Jong et al. 
        stress_components.at(0).at(j) = stress_tensor[1][1];
        stress_components.at(1).at(j) = stress_tensor[2][2];
        stress_components.at(2).at(j) = stress_tensor[3][3];
        stress_components.at(3).at(j) = (stress_tensor[2][3] + stress_tensor[3][2]) / 2.0;
        stress_components.at(4).at(j) = (stress_tensor[1][3] + stress_tensor[3][1]) / 2.0;
        stress_components.at(5).at(j) = (stress_tensor[1][2] + stress_tensor[2][1]) / 2.0;
      }
      // Stress components for smallest deformations only; used to check linearity
      for (uint k = 0; k < stress_components_half.size(); k++) {
        stress_components_half.at(k).resize(normal_deformations_half_tofit.at(i-1).size());
      }
      for (uint j = 0; j < normal_deformations_half_tofit.at(i-1).size(); j++) {   
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "i = " << i << ", j = " << j << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Normal strain tensor = " << endl;
          // OBSOLETE aus << AEL_data.normal_strain.at(i-1).at(j + halffirstpoint) << endl;
          aus << _AELSTR_MESSAGE_ + "Normal strain deformations = " << endl;
          aus << normal_deformations_half_tofit.at(i-1).at(j) << endl;
          aus << _AELSTR_MESSAGE_ + "Normal stress tensor = " << endl; 
          aus << normal_stress_half_tofit.at(i-1).at(j) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // OBSOLETE strain_components = AEL_data.normal_strain.at(i-1).at(j + halffirstpoint) - str_identity;
        // OBSOLETE stress_tensor = AEL_data.normal_stress.at(i-1).at(j + halffirstpoint);
        stress_tensor = normal_stress_half_tofit.at(i-1).at(j);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Strain tensor = " << endl;
          // OBSOLETE aus << strain_components << endl;
          aus << _AELSTR_MESSAGE_ + "Stress tensor = " << endl; 
          aus << stress_tensor << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Stress components for stress-strain fitting; see Eq. 1, Phys. Rev. Materials 1, 015401 (2017) and Scientific Data 2, 150009 (2015)
        stress_components_half.at(0).at(j) = stress_tensor[1][1];
        stress_components_half.at(1).at(j) = stress_tensor[2][2];
        stress_components_half.at(2).at(j) = stress_tensor[3][3];
        stress_components_half.at(3).at(j) = (stress_tensor[2][3] + stress_tensor[3][2]) / 2.0;
        stress_components_half.at(4).at(j) = (stress_tensor[1][3] + stress_tensor[3][1]) / 2.0;
        stress_components_half.at(5).at(j) = (stress_tensor[1][2] + stress_tensor[2][1]) / 2.0;
      }    
      // OBSOLETE if(AEL_data.fitstrainorigin || AEL_data.fitrelaxedstruct) {
      // OBSOLETE uint j = normal_deformations_origin.size() - 1;
      // OBSOLETE if(AEL_data.calcstrainorigin || AEL_data.fitrelaxedstruct) {
      // OBSOLETE   stress_tensor = AEL_data.origin_stress.at(0);
      // OBSOLETE   stress_components.at(0).at(j) = stress_tensor[1][1];
      // OBSOLETE   stress_components.at(1).at(j) = stress_tensor[2][2];
      // OBSOLETE   stress_components.at(2).at(j) = stress_tensor[3][3];
      // OBSOLETE   stress_components.at(3).at(j) = (stress_tensor[2][3] + stress_tensor[3][2]) / 2.0;
      // OBSOLETE   stress_components.at(4).at(j) = (stress_tensor[1][3] + stress_tensor[3][1]) / 2.0;
      // OBSOLETE   stress_components.at(5).at(j) = (stress_tensor[1][2] + stress_tensor[2][1]) / 2.0;
      // OBSOLETE } else {
      // OBSOLETE   stress_components.at(0).at(j) = 0.0;
      // OBSOLETE   stress_components.at(1).at(j) = 0.0;
      // OBSOLETE   stress_components.at(2).at(j) = 0.0;
      // OBSOLETE   stress_components.at(3).at(j) = 0.0;
      // OBSOLETE   stress_components.at(4).at(j) = 0.0;
      // OBSOLETE   stress_components.at(5).at(j) = 0.0;
      // OBSOLETE }	  
      // OBSOLETE for (uint k = 0; k < 6; k++) {
      // OBSOLETE   aelerror =  cij_fit(normal_deformations_origin, stress_components.at(k), Cij, AEL_data.normalpolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
      // Convert units from kB to GPa
      // OBSOLETE Cij = 0.1 * Cij;
      // OBSOLETE // OBSOLETE Remove any values less than 1.0; these are likely due to noise
      // OBSOLETE // OBSOLETE if(aurostd::abs(Cij) < 1.0) {
      // OBSOLETE   // OBSOLETE Cij = 0.0;
      // OBSOLETE // OBSOLETE }
      // OBSOLETE elastic_const[k+1][i] = Cij;
      // OBSOLETE }      
      // OBSOLETE } else
      // Checks if number of complete deformation runs is less than number of initial runs
      npolyfitorder = AEL_data.normalpolyfitorder;
      npolyfitorderhalf = midpointnormal;
      nupolyfitorder = npolyfitorder;
      if(AEL_data.normal_deformations_complete.at(i-1).size() < AEL_data.normal_deformations.size()) {
        if(AEL_data.normal_deformations_complete.at(i-1).size() < nupolyfitorder) {
          npolyfitorder = AEL_data.normal_deformations_complete.at(i-1).size();
          // OBSOLETE if(normal_deformations_tofit.at(i-1).size() < AEL_data.normal_deformations.size())
          // OBSOLETE if(normal_deformations_tofit.at(i-1).size() < nupolyfitorder)
          // OBSOLETE   npolyfitorder = normal_deformations_tofit.at(i-1).size();
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Polynomial order for fitting less than number of available points" << endl; 
          aus << _AELSTR_WARNING_ + "Resetting polynomial order for fitting to " << npolyfitorder << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } 
      }
      for (uint k = 0; k < 6; k++) {
        // OBSOLETE aelerror =  cij_fit(AEL_data.normal_deformations, stress_components.at(k), Cij, AEL_data.normalpolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
        // OBSOLETE aelerror =  cij_fit(normal_deformations_tofit.at(i-1), stress_components.at(k), Cij, AEL_data.normalpolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
        aelerror =  cij_fit(normal_deformations_tofit.at(i-1), stress_components.at(k), Cij, npolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
        // Convert units from kB to GPa
        Cij = 0.1 * Cij;
        // OBSOLETE Remove any values less than 1.0; these are likely due to noise
        // OBSOLETE if(aurostd::abs(Cij) < 1.0) {
        // OBSOLETE Cij = 0.0;
        // OBSOLETE }
        elastic_const[k+1][i] = Cij;
      }
      // Fit small deformations only to check stress-strain linearity
      for (uint k = 0; k < 6; k++) {
        // OBSOLETE aelerror =  cij_fit(AEL_data.normal_deformations_half, stress_components_half.at(k), Cijhalf, nlinfit, AEL_data.gaussxm_debug, FileMESSAGE);
        // OBSOLETE aelerror =  cij_fit(normal_deformations_half_tofit.at(i-1), stress_components_half.at(k), Cijhalf, nlinfit, AEL_data.gaussxm_debug, FileMESSAGE);
        aelerror =  cij_fit(normal_deformations_half_tofit.at(i-1), stress_components_half.at(k), Cijhalf, npolyfitorderhalf, AEL_data.gaussxm_debug, FileMESSAGE);
        // Convert units from kB to GPa
        Cijhalf = 0.1 * Cijhalf;
        elastic_const_half[k+1][i] = Cijhalf;
      }     
    }
    if(AEL_data.normal_strain.size() < 3) {
      if(AEL_data.normal_strain.size() == 1) {
        for (uint i = 2; i <= 3; i++) {
          for (uint k = 1; k <= 3; k++) {
            if(k == i) {
              elastic_const[k][i] = elastic_const[1][1];
              elastic_const_half[k][i] = elastic_const_half[1][1];
            } else {
              elastic_const[k][i] = elastic_const[2][1];
              elastic_const_half[k][i] = elastic_const_half[2][1];
            }
          }
          for (uint k = 4; k <= 6; k++) {
            elastic_const[k][i] = elastic_const[4][1];
            elastic_const_half[k][i] = elastic_const_half[4][1];
          }
        }
      } else if(AEL_data.normal_strain.size() == 2) {
        for (uint k = 1; k <= 3; k++) {
          if(k == 2) {
            elastic_const[k][2] = elastic_const[1][1];
            elastic_const_half[k][2] = elastic_const_half[1][1];
          } else {
            elastic_const[k][2] = elastic_const[2][1];
            elastic_const_half[k][2] = elastic_const_half[2][1];
          }
        }
      }
    }	    
    // OBSOLETE for (uint k = 0; k < stress_components.size(); k++) {
    // OBSOLETE   stress_components.at(k).resize(AEL_data.shear_deformations.size());
    // OBSOLETE }
    // OBSOLETE for (uint k = 0; k < stress_components_half.size(); k++) {
    // OBSOLETE   stress_components_half.at(k).resize(AEL_data.shear_deformations_half.size());
    // OBSOLETE }
    // OBSOLETE if(AEL_data.negstrain) {
    // OBSOLETE   halffirstpoint = AEL_data.shear_deformations.size() / 4;
    // OBSOLETE }
    // OBSOLETE if(AEL_data.fitstrainorigin || AEL_data.fitrelaxedstruct) {
    // OBSOLETE   for (uint k = 0; k < stress_components.size(); k++) {
    // OBSOLETE stress_components.at(k).resize(shear_deformations_origin.size());
    // OBSOLETE   }      
    // OBSOLETE }
    // OBSOLETE for (uint j = 0; j < AEL_data.shear_deformations.size(); j++) {
    // OBSOLETE  shear_deformations_doubled.push_back(2.0*AEL_data.shear_deformations.at(j));
    // OBSOLETE }
    // Fit shear stress-strains
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Fitting shear stress-strain data" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
      for (uint k = 0; k < stress_components.size(); k++) {
        stress_components.at(k).resize(shear_deformations_tofit.at(i-1).size());
      }
      if(LVERBOSE) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "shear_deformations_tofit.at(i-1).size() = " << shear_deformations_tofit.at(i-1).size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
      }
      for (uint j = 0; j < shear_deformations_tofit.at(i-1).size(); j++) {   
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "i = " << i << ", j = " << j << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Shear strain tensor = " << endl;
          // OBSOLETE aus << AEL_data.shear_strain.at(i-1).at(j) << endl;
          aus << _AELSTR_MESSAGE_ + "Shear strain deformations = " << endl;
          aus << shear_deformations_tofit.at(i-1).at(j) << endl;
          aus << _AELSTR_MESSAGE_ + "Shear stress tensor = " << endl; 
          aus << AEL_data.shear_stress.at(i-1).at(j) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // OBSOLETE strain_components = AEL_data.shear_strain.at(i-1).at(j) - str_identity;
        // OBSOLETE stress_tensor = AEL_data.shear_stress.at(i-1).at(j);
        stress_tensor = shear_stress_tofit.at(i-1).at(j);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Strain tensor = " << endl;
          // OBSOLETE aus << strain_components << endl;
          aus << _AELSTR_MESSAGE_ + "Stress tensor = " << endl; 
          aus << stress_tensor << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Stress components for stress-strain fitting; see Eq. 1, M. de Jong et al. 
        stress_components.at(0).at(j) = stress_tensor[1][1];
        stress_components.at(1).at(j) = stress_tensor[2][2];
        stress_components.at(2).at(j) = stress_tensor[3][3];
        stress_components.at(3).at(j) = (stress_tensor[2][3] + stress_tensor[3][2]) / 2.0;
        stress_components.at(4).at(j) = (stress_tensor[1][3] + stress_tensor[3][1]) / 2.0;
        stress_components.at(5).at(j) = (stress_tensor[1][2] + stress_tensor[2][1]) / 2.0;
      }
      // Stress components for smallest deformations only; used to check linearity
      for (uint k = 0; k < stress_components_half.size(); k++) {
        stress_components_half.at(k).resize(shear_deformations_half_tofit.at(i-1).size());
      }
      for (uint j = 0; j < shear_deformations_half_tofit.at(i-1).size(); j++) {   
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "i = " << i << ", j = " << j << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Shear strain tensor = " << endl;
          // OBSOLETE aus << AEL_data.shear_strain.at(i-1).at(j + halffirstpoint) << endl;
          aus << _AELSTR_MESSAGE_ + "Shear strain deformations = " << endl;
          aus << shear_deformations_half_tofit.at(i-1).at(j) << endl;
          aus << _AELSTR_MESSAGE_ + "Shear stress tensor = " << endl; 
          aus << shear_stress_half_tofit.at(i-1).at(j) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // OBSOLETE strain_components = AEL_data.shear_strain.at(i-1).at(j + halffirstpoint) - str_identity;
        // OBSOLETE stress_tensor = AEL_data.shear_stress.at(i-1).at(j + halffirstpoint);
        stress_tensor = shear_stress_half_tofit.at(i-1).at(j);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          // OBSOLETE aus << _AELSTR_MESSAGE_ + "Strain tensor = " << endl;
          // OBSOLETE aus << strain_components << endl;
          aus << _AELSTR_MESSAGE_ + "Stress tensor = " << endl; 
          aus << stress_tensor << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Stress components for stress-strain fitting; see Eq. 1, M. de Jong et al. 
        stress_components_half.at(0).at(j) = stress_tensor[1][1];
        stress_components_half.at(1).at(j) = stress_tensor[2][2];
        stress_components_half.at(2).at(j) = stress_tensor[3][3];
        stress_components_half.at(3).at(j) = (stress_tensor[2][3] + stress_tensor[3][2]) / 2.0;
        stress_components_half.at(4).at(j) = (stress_tensor[1][3] + stress_tensor[3][1]) / 2.0;
        stress_components_half.at(5).at(j) = (stress_tensor[1][2] + stress_tensor[2][1]) / 2.0;
      }
      // OBSOLETE if(AEL_data.fitstrainorigin || AEL_data.fitrelaxedstruct) {
      // OBSOLETE uint j = shear_deformations_origin.size() - 1;
      // OBSOLETE if(AEL_data.calcstrainorigin || AEL_data.fitrelaxedstruct) {
      // OBSOLETE   stress_tensor = AEL_data.origin_stress.at(0);
      // OBSOLETE   stress_components.at(0).at(j) = stress_tensor[1][1];
      // OBSOLETE   stress_components.at(1).at(j) = stress_tensor[2][2];
      // OBSOLETE   stress_components.at(2).at(j) = stress_tensor[3][3];
      // OBSOLETE   stress_components.at(3).at(j) = (stress_tensor[2][3] + stress_tensor[3][2]) / 2.0;
      // OBSOLETE   stress_components.at(4).at(j) = (stress_tensor[1][3] + stress_tensor[3][1]) / 2.0;
      // OBSOLETE   stress_components.at(5).at(j) = (stress_tensor[1][2] + stress_tensor[2][1]) / 2.0;
      // OBSOLETE } else {
      // OBSOLETE   stress_components.at(0).at(j) = 0.0;
      // OBSOLETE   stress_components.at(1).at(j) = 0.0;
      // OBSOLETE   stress_components.at(2).at(j) = 0.0;
      // OBSOLETE   stress_components.at(3).at(j) = 0.0;
      // OBSOLETE   stress_components.at(4).at(j) = 0.0;
      // OBSOLETE   stress_components.at(5).at(j) = 0.0;
      // OBSOLETE }	  
      // OBSOLETE for (uint k = 0; k < 6; k++) {
      // OBSOLETE   aelerror =  cij_fit(shear_deformations_origin, stress_components.at(k), Cij, AEL_data.shearpolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
      // Convert units from kB to GPa
      // OBSOLETE Cij = 0.1 * Cij;
      // OBSOLETE // OBSOLETE Remove any values less than 1.0; these are likely due to noise
      // OBSOLETE // OBSOLETE if(aurostd::abs(Cij) < 1.0) {
      // OBSOLETE // OBSOLETE Cij = 0.0;
      // OBSOLETE // OBSOLETE }
      // OBSOLETE  elastic_const[k+1][i+3] = Cij;
      // OBSOLETE }      
      // OBSOLETE } else
      npolyfitorder = AEL_data.shearpolyfitorder;
      npolyfitorderhalf = midpointshear;
      nupolyfitorder = npolyfitorder;
      if(AEL_data.shear_deformations_complete.at(i-1).size() < AEL_data.shear_deformations.size()) {
        if(AEL_data.shear_deformations_complete.at(i-1).size() < nupolyfitorder) {
          npolyfitorder = AEL_data.shear_deformations_complete.at(i-1).size();
          // OBSOLETE if(shear_deformations_tofit.at(i-1).size() < AEL_data.shear_deformations.size())
          // OBSOLETE 	if(shear_deformations_tofit.at(i-1).size() < nupolyfitorder)
          // OBSOLETE   npolyfitorder = shear_deformations_tofit.at(i-1).size();
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Polynomial order for fitting less than number of available points" << endl; 
          aus << _AELSTR_WARNING_ + "Resetting polynomial order for fitting to " << npolyfitorder << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } 
      }
      for (uint k = 0; k < 6; k++) {
        // OBSOLETE aelerror =  cij_fit(AEL_data.shear_deformations, stress_components.at(k), Cij, AEL_data.shearpolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
        // OBSOLETE aelerror =  cij_fit(shear_deformations_tofit.at(i-1), stress_components.at(k), Cij, AEL_data.shearpolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
        aelerror =  cij_fit(shear_deformations_tofit.at(i-1), stress_components.at(k), Cij, npolyfitorder, AEL_data.gaussxm_debug, FileMESSAGE);
        // Convert units from kB to GPa
        Cij = 0.1 * Cij;
        // OBSOLETE if(aurostd::abs(Cij) < 1.0) {
        // OBSOLETE Cij = 0.0;
        // OBSOLETE }
        elastic_const[k+1][i+3] = Cij;
      }      
      for (uint k = 0; k < 6; k++) {
        // OBSOLETE aelerror =  cij_fit(AEL_data.shear_deformations_half, stress_components_half.at(k), Cijhalf, nlinfit, AEL_data.gaussxm_debug, FileMESSAGE);
        // OBSOLETE aelerror =  cij_fit(shear_deformations_half_tofit.at(i-1), stress_components_half.at(k), Cijhalf, nlinfit, AEL_data.gaussxm_debug, FileMESSAGE);
        aelerror =  cij_fit(shear_deformations_half_tofit.at(i-1), stress_components_half.at(k), Cijhalf, npolyfitorderhalf, AEL_data.gaussxm_debug, FileMESSAGE);
        // Convert units from kB to GPa
        Cijhalf = 0.1 * Cijhalf;
        elastic_const_half[k+1][i+3] = Cijhalf;
      }      
    }
    if(AEL_data.shear_strain.size() < 3) {
      if(AEL_data.shear_strain.size() == 1) {
        for (uint i = 5; i <= 6; i++) {
          for (uint k = 4; k <= 6; k++) {
            if(k == i) {
              elastic_const[k][i] = elastic_const[4][4];
              elastic_const_half[k][i] = elastic_const_half[4][4];
            } else {
              elastic_const[k][i] = elastic_const[5][4];
              elastic_const_half[k][i] = elastic_const_half[5][4];
            }
          }
          for (uint k = 1; k <= 3; k++) {
            elastic_const[k][i] = elastic_const[1][4];
            elastic_const_half[k][i] = elastic_const_half[1][4];
          }
        }
      } else if(AEL_data.shear_strain.size() == 2) {
        for (uint k = 4; k <= 6; k++) {
          if(k == 2) {
            elastic_const[k][5] = elastic_const[4][4];
            elastic_const_half[k][5] = elastic_const_half[4][4];
          } else {
            elastic_const[k][5] = elastic_const[5][4];
            elastic_const_half[k][5] = elastic_const_half[5][4];
          }
        }
      }
    }	   
    // Remove any off-diagonal values which are less than 2% of the largest value; these are likely due to noise
    // First find the matrix element with the largest value
    for (uint i = 1; i <= 6; i++) {
      for (uint j = 1; j <= 6; j++) {
        if(Cijmax < aurostd::abs(elastic_const[i][j])) {
          Cijmax = elastic_const[i][j];
        }
      }
    }
    // Next set tolerance to 0.02 times Cijmax or 1.0; whichever is larger
    Cijtol = aurostd::abs(0.02 * Cijmax);
    if(Cijtol < 1.0) {
      Cijtol = 1.0;
    }
    // Next set off-diagonal elements with an absolute value less than than the tolerance to zero
    for (uint i = 1; i <= 6; i++) {
      for (uint j = 1; j <= 6; j++) {
        if((i != j) && (aurostd::abs(elastic_const[i][j]) < Cijtol)) {
          elastic_const[i][j] = 0.0;
        }
        if((i != j) && (aurostd::abs(elastic_const_half[i][j]) < Cijtol)) {
          elastic_const_half[i][j] = 0.0;
        }
      }
    }     
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Elastic constants = " << endl;
    aus << elastic_const << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Elastic constants (half strain) = " << endl;
    aus << elastic_const_half << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // Symmetrize elastic constant matrix and add to AEL_data class
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Allocating tensor storage" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    AEL_data.elastic_tensor.resize(6);
    AEL_data.elastic_tensor_half.resize(6);
    AEL_data.compliance_tensor.resize(6);
    AEL_data.compliance_tensor_half.resize(6);
    AEL_data.elastic_eigenvalues.resize(6);
    for (uint i = 0; i < 6; i++) {
      AEL_data.elastic_tensor.at(i).resize(6);
      AEL_data.elastic_tensor_half.at(i).resize(6);
      AEL_data.compliance_tensor.at(i).resize(6);
      AEL_data.compliance_tensor_half.at(i).resize(6);
    }      
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Storing elastic tensor" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    for (uint i = 0; i < 6; i++) {
      for (uint j = 0; j < 6; j++) {
        AEL_data.elastic_tensor.at(i).at(j) = elastic_const[i+1][j+1];
        AEL_data.elastic_tensor_half.at(i).at(j) = elastic_const_half[i+1][j+1];
      }
    }
    // Symmetrize elastic stiffness tensor if requested by user
    if (AEL_data.symmetrize_elastic_tensor) {
      for (uint i = 0; i < 6; i++) {
        for (uint j = 0; j < 6; j++) {
          AEL_data.elastic_tensor.at(i).at(j) = (elastic_const[i+1][j+1] + elastic_const[j+1][i+1]) / 2.0;
          AEL_data.elastic_tensor_half.at(i).at(j) = (elastic_const_half[i+1][j+1] + elastic_const_half[j+1][i+1]) / 2.0;
        }
      }
    }    
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Elastic constant tensor = " << endl;
    for (uint i = 0; i < 6; i++) {
      for (uint j = 0; j < 6; j++) {
        aus <<  AEL_data.elastic_tensor.at(i).at(j) << "\t";
      }
      aus << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aurostd::cematrix eltens(elastic_const);
    aurostd::cematrix eltenshalf(elastic_const_half);
    // OBSOLETE aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AELSTR_MESSAGE_ + "Elastic tensor (cematrix) = " << endl;
    // OBSOLETE aus << eltens << endl;
    // OBSOLETE aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // Calculate elastic stiffness tensor eigenvalues to test mechanical stability
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Computing eigenvalues of elastic tensor" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OBSOLETE elastic_eigenvals = eltens.EigenValues();
    aurostd::eigen(elastic_const, elastic_eigenval_real, elastic_eigenval_imag);
    elastic_eigenvals = elastic_eigenval_real;
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Elastic tensor eigenvalues = " << endl;
    aus << elastic_eigenvals << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    for (uint i = 0; i < 6; i++) {
      AEL_data.elastic_eigenvalues.at(i) = elastic_eigenvals[i+1];
    }    
    // Invert elastic stiffness tensor to obtain compliance tensor
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Inverting elastic tensor" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    elastic_const_inv = eltens.InverseMatrix();
    elastic_const_half_inv = eltenshalf.InverseMatrix();
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Storing compliance tensor" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    for (uint i = 0; i < 6; i++) {
      for (uint j = 0; j < 6; j++) {
        AEL_data.compliance_tensor.at(i).at(j) = elastic_const_inv[i+1][j+1];
        AEL_data.compliance_tensor_half.at(i).at(j) = elastic_const_half_inv[i+1][j+1];
      }
    }
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Compliance tensor = " << endl;
    for (uint i = 0; i < 6; i++) {
      for (uint j = 0; j < 6; j++) {
        aus <<  AEL_data.compliance_tensor.at(i).at(j) << "\t";
      }
      aus << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
    return aelerror;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::cij_fit
// ***************************************************************************
namespace AEL_functions {
  //
  // cij_fit: Calls AGL::polynom_fit to fit x-y obtain elastic constant element
  // Elastic constant element is first derivative of polynomial at x=0
  // See Phys. Rev. Materials 1, 015401 (2017) and Scientific Data 2, 150009 (2015) for details
  //
  uint cij_fit(vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, double& Cij, int& npolycoeffwork, bool gxmdebug, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    uint first_entry_tofit = 0;
    uint last_entry_tofit = ydata_to_fit.size() - 1;
    double rms = 0;
    vector<double> polycoeffwork, weight;
    // OBSOLETE polycoeffwork.resize(ydata_to_fit.size());
    polycoeffwork.resize(npolycoeffwork);
    for (uint i = 0; i < polycoeffwork.size(); i++) {
      polycoeffwork.at(i) = 0.0;
    }
    // OBSOLETE int npolycoeffwork = ydata_to_fit.size();
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Fitting x-y data" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    if(LVERBOSE) {
      aurostd::StringstreamClean(aus);
      for (uint i = 0; i < ydata_to_fit.size(); i++) {
        aus << xdata_to_fit.at(i) << "\t" << ydata_to_fit.at(i) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    } 
    // Initialize the weights
    weight.resize(ydata_to_fit.size());
    for (uint i = 0; i < weight.size(); i++) {
      weight.at(i) = 1.0;
    }
    // Fit stress-strain data to obtain polynomial
    uint pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, xdata_to_fit, ydata_to_fit, weight, rms, npolycoeffwork, polycoeffwork, gxmdebug, FileMESSAGE);
    if(pferr != 0) {
      return 2;
    }
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "npolycoeffwork = " << npolycoeffwork << endl;
    for (uint i = 0; i < polycoeffwork.size(); i++) {
      aus << _AELSTR_MESSAGE_ + "polycoeffwork.at(" << i << ") = " << polycoeffwork.at(i) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // Evaluate slope of polynomial at x=0 to obtain linear-response elastic constant element
    // Polynomial is of form: y = p[0] + p[1]x + p[2]x^2 + ....
    // First derivative at x=0 is: dy/dx = p[1]
    Cij = polycoeffwork.at(1);
    return 0;    
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::elasticitycheck
// ***************************************************************************
namespace AEL_functions {
  //
  // elasticitycheck: Checks symmetry and consistency of elastic constants
  // See F. Mouhat and F.-X. Coudert, Phys. Rev. B 90, 224104 (2014) for details
  //
  uint elasticitycheck(_AEL_data& AEL_data, _xvasp& xvasp, ofstream& FileMESSAGE) {
    // OBSOLETE bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    double lowerbound = 1.0 - AEL_data.symtolfrac;
    double upperbound = 1.0 + AEL_data.symtolfrac;
    double tol = 1.0e-8;
    bool tetI = false;
    bool tetII = false;
    bool rhlI = false;
    bool rhlII = false;
    // Independent direction in crystal lattice; take z-direction as default
    int inddir = 2;
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Checking symmetry and consistency of elasticity tensor" << endl;
    aus << _AELSTR_MESSAGE_ + "Elastic constant tensor = " << endl; 
    for (uint i = 0; i < 6; i++) {
      for (uint j = 0; j < 6; j++) {
        aus <<  AEL_data.elastic_tensor.at(i).at(j) << "\t";
      }
      aus << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);   
    // First check symmetry about diagonal - tensor should be symmetric for all lattice types
    for (uint i = 0; i < 6; i++) {
      for (uint j = 0; j < 6; j++) {
        if(aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) < tol) {
          if(aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > tol) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_" << i+1 << j+1 << endl;
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << endl;
            aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << endl;
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " / C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(j).at(i) << endl;	    
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: asymmetry in elastic tensor: increase strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        } else if((AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(j).at(i) > upperbound) || (AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(j).at(i) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_" << i+1 << j+1 << endl;
          aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << endl;
          aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << endl;
          aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " / C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(j).at(i) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: asymmetry in elastic tensor: increase strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
      }
    }    
    // Next check symmetry relations for cubic lattices
    if((xvasp.str.bravais_lattice_type == "FCC") || (xvasp.str.bravais_lattice_type == "BCC") || (xvasp.str.bravais_lattice_type == "CUB")) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice has cubic symmetry" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (uint i = 0; i < 3; i++) {
        // Check diagonal elements in upper left quadrant are all equal
        if((AEL_data.elastic_tensor.at(i).at(i) / AEL_data.elastic_tensor.at(0).at(0) > upperbound) || (AEL_data.elastic_tensor.at(i).at(i) / AEL_data.elastic_tensor.at(0).at(0) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_11" << endl;
          aus << _AELSTR_WARNING_ + "C_" << i+1 << i+1 << " / C_11 = " << AEL_data.elastic_tensor.at(i).at(i) / AEL_data.elastic_tensor.at(0).at(0) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Check off-diagonal elements in upper left quadrant are all equal
        for (uint j = 0; j < 3; j++) {
          if((i != j) && ((AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(0).at(1) > upperbound) || (AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(0).at(1) < lowerbound))) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_12" << endl;
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " / C_12 = " << AEL_data.elastic_tensor.at(i).at(j) / AEL_data.elastic_tensor.at(0).at(0) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
      }
      for (uint i = 3; i < 6; i++) {
        // Check diagonal elements in lower right quadrant are all equal
        if((AEL_data.elastic_tensor.at(i).at(i) / AEL_data.elastic_tensor.at(3).at(3) > upperbound) || (AEL_data.elastic_tensor.at(i).at(i) / AEL_data.elastic_tensor.at(3).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
          aus << _AELSTR_WARNING_ + "C_" << i+1 << i+1 << " / C_44 = " << AEL_data.elastic_tensor.at(i).at(i) / AEL_data.elastic_tensor.at(3).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Check off-diagonal elements in lower right quadrant are zero
        for (uint j = 3; j < 6; j++) {
          if((i != j) && (aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
      }
      // Check upper right and lower left quadrants are zero
      for (uint i = 3; i < 6; i++) {
        for (uint j = 0; j < 3; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
      }
      // Next check symmetry relations for hexagonal and tetragonal lattices
    } else if((xvasp.str.bravais_lattice_type == "HEX") || (xvasp.str.bravais_lattice_type == "TET") || (xvasp.str.bravais_lattice_type == "BCT")) {  
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice has hexagonal or tetragonal symmetry" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // For tetragonal lattice types, check whether lattice is of type I (space group # 89 to 142) or of type II (space group # 75 to 88) 
      // This changes the number of independent elastic constants (see Phys. Rev. B 90, 224104 (2014) and http://en.wikipedia.org/wiki/List_of_space_groups for details)  
      if((xvasp.str.bravais_lattice_type == "TET") || (xvasp.str.bravais_lattice_type == "BCT")) {
        // OBSOLETE bool sg_manual_toggle = false;
        // OBSOLETE vector<string> sg_args;
        // OBSOLETE sg_args.push_back("dummyfiller");
        // OBSOLETE uint ael_spacegroupnumber = xvasp.str.SpaceGroup_ITC(sg_manual_toggle, sg_args);
        uint ael_spacegroupnumber = xvasp.str.SpaceGroup_ITC();
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice space group number = " << ael_spacegroupnumber << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
        // OBSOLETE if((xvasp.str.spacegroupnumber > 74) && (xvasp.str.spacegroupnumber < 89)) {
        if((ael_spacegroupnumber > 74) && (ael_spacegroupnumber < 89)) {
          tetII = true;
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Lattice is of type tetragonal II" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
          // OBSOLETE } else if((xvasp.str.spacegroupnumber > 88) && (xvasp.str.spacegroupnumber < 143))
      } else if((ael_spacegroupnumber > 88) && (ael_spacegroupnumber < 143)) {
        tetI = true;
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice is of type tetragonal I" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Lattice type = " << xvasp.str.bravais_lattice_type << " but space group # = " <<  xvasp.str.spacegroupnumber << endl;
        aus << _AELSTR_WARNING_ + "This is outside the range of 89 to 142 for tetragonal space groups" << endl;
        aus << _AELSTR_WARNING_ + "Possible mismatch between lattice type and space group" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
      }
      } 
      // First, determine which direction is the independent one
      if((AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) < upperbound) && (AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) > lowerbound)) {
        inddir = 2;
      }	else if((AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(0).at(0) < upperbound) && (AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(0).at(0) > lowerbound)) {
        inddir = 1;
      }	else if((AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(1).at(1) < upperbound) && (AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(1).at(1) > lowerbound)) {
        inddir = 0;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_11: no two diagonal elements are equal for this material" << endl;
        aus << _AELSTR_WARNING_ + "C_22 / C_11 = " << AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) << endl;
        aus << _AELSTR_WARNING_ + "C_33 / C_11 = " << AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(0).at(0) << endl;
        aus << _AELSTR_WARNING_ + "C_33 / C_22 = " << AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(1).at(1) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // Next check off-diagonal elements in lower right quadrant are all zero
      for (uint i = 3; i < 6; i ++) {
        for (uint j = 3; j < 6; j++) {
          if((i != j) && (aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
      }    
      // OBSOLETE if((AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) > upperbound) || (AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) < lowerbound)) {
      // OBSOLETE aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_11" << endl;
      // OBSOLETE aus << _AELSTR_WARNING_ + "C_22 / C_11 = " << AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) << endl;
      // OBSOLETE aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE }
      // OBSOLETE if((AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(0).at(1) > upperbound) || (AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(0).at(1) < lowerbound)) {
      // OBSOLETE aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_12" << endl;
      // OBSOLETE aus << _AELSTR_WARNING_ + "C_21 / C_12 = " << AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(0).at(1) << endl;
      // OBSOLETE aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
      // OBSOLETE aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE }
      // OBSOLETE for (uint i = 0; i < 2; i++)
      // Next check off-diagonal elements in upper left quadrant and diagonal elements in lower right quadrant   
      if(inddir == 2) {
        if((AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) > upperbound) || (AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
          aus << _AELSTR_WARNING_ + "C_23 / C_13 = " << AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) > upperbound) || (AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
          aus << _AELSTR_WARNING_ + "C_32 / C_31 = " << AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(4).at(4) / AEL_data.elastic_tensor.at(3).at(3) > upperbound) || (AEL_data.elastic_tensor.at(4).at(4) / AEL_data.elastic_tensor.at(3).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
          aus << _AELSTR_WARNING_ + "C_55 / C_44 = " << AEL_data.elastic_tensor.at(4).at(4) / AEL_data.elastic_tensor.at(3).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else if(inddir == 1) { 
        if((AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(1).at(0) > upperbound) || (AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(1).at(0) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
          aus << _AELSTR_WARNING_ + "C_23 / C_13 = " << AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(0).at(1) > upperbound) || (AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(0).at(1) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
          aus << _AELSTR_WARNING_ + "C_32 / C_31 = " << AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(3).at(3) > upperbound) || (AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(3).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
          aus << _AELSTR_WARNING_ + "C_66 / C_44 = " << AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(3).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else if(inddir == 0) { 
        if((AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(2).at(0) > upperbound) || (AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(2).at(0) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
          aus << _AELSTR_WARNING_ + "C_21 / C_31 = " << AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(2).at(0) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(0).at(1) / AEL_data.elastic_tensor.at(0).at(2) > upperbound) || (AEL_data.elastic_tensor.at(0).at(1) / AEL_data.elastic_tensor.at(0).at(2) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
          aus << _AELSTR_WARNING_ + "C_12 / C_13 = " << AEL_data.elastic_tensor.at(0).at(1) / AEL_data.elastic_tensor.at(0).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(4).at(4) > upperbound) || (AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(4).at(4) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
          aus << _AELSTR_WARNING_ + "C_66 / C_55 = " << AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(4).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      // Check that off-diagonal quadrants are all zero for hexagonal and tetragonal (I) class
      if((xvasp.str.bravais_lattice_type == "HEX") || (tetI)) {
        for (uint i = 3; i < 6; i++) {
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
        }
      } else if(tetII) {
        if(inddir == 2) {
          if((AEL_data.elastic_tensor.at(0).at(5) / (-AEL_data.elastic_tensor.at(1).at(5)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(5) / (-AEL_data.elastic_tensor.at(1).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_16" << endl;
            aus << _AELSTR_WARNING_ + "C_16 / C_26 = " << AEL_data.elastic_tensor.at(0).at(5) / AEL_data.elastic_tensor.at(1).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(5).at(0) / (-AEL_data.elastic_tensor.at(5).at(1)) > upperbound) || (AEL_data.elastic_tensor.at(5).at(0) / (-AEL_data.elastic_tensor.at(5).at(1)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_61" << endl;
            aus << _AELSTR_WARNING_ + "C_61 / C_62 = " << AEL_data.elastic_tensor.at(5).at(0) / AEL_data.elastic_tensor.at(5).at(1) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint i = 3; i < 5; i++) {
            for (uint j = 0; j < 3; j++) {
              if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }	  
            }
          }
          if(AEL_data.elastic_tensor.at(5).at(2) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_63 = " << AEL_data.elastic_tensor.at(5).at(2) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(2).at(5) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_36 = " << AEL_data.elastic_tensor.at(2).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else if(inddir == 1) {
          if((AEL_data.elastic_tensor.at(0).at(4) / (-AEL_data.elastic_tensor.at(2).at(4)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(4) / (-AEL_data.elastic_tensor.at(2).at(4)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_16" << endl;
            aus << _AELSTR_WARNING_ + "C_16 / C_26 = " << AEL_data.elastic_tensor.at(0).at(4) / AEL_data.elastic_tensor.at(2).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(4).at(0) / (-AEL_data.elastic_tensor.at(4).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(4).at(0) / (-AEL_data.elastic_tensor.at(4).at(2)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_61" << endl;
            aus << _AELSTR_WARNING_ + "C_61 / C_62 = " << AEL_data.elastic_tensor.at(4).at(0) / AEL_data.elastic_tensor.at(4).at(2) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_3" << j+1 << " = " << AEL_data.elastic_tensor.at(3).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(3)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "3 = " << AEL_data.elastic_tensor.at(j).at(3) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_5" << j+1 << " = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(5)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "5 = " << AEL_data.elastic_tensor.at(j).at(5) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          if(AEL_data.elastic_tensor.at(4).at(1) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_52 = " << AEL_data.elastic_tensor.at(4).at(1) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(1).at(4) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_25 = " << AEL_data.elastic_tensor.at(1).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else if(inddir == 0) {
          if((AEL_data.elastic_tensor.at(1).at(3) / (-AEL_data.elastic_tensor.at(2).at(3)) > upperbound) || (AEL_data.elastic_tensor.at(1).at(3) / (-AEL_data.elastic_tensor.at(2).at(3)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_16" << endl;
            aus << _AELSTR_WARNING_ + "C_16 / C_26 = " << AEL_data.elastic_tensor.at(1).at(3) / AEL_data.elastic_tensor.at(1).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(3).at(1) / (-AEL_data.elastic_tensor.at(3).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(3).at(1) / (-AEL_data.elastic_tensor.at(3).at(2)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_61" << endl;
            aus << _AELSTR_WARNING_ + "C_61 / C_62 = " << AEL_data.elastic_tensor.at(3).at(1) / AEL_data.elastic_tensor.at(3).at(2) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint i = 4; i < 6; i++) {
            for (uint j = 0; j < 3; j++) {
              if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }	  
            }
          }
          if(AEL_data.elastic_tensor.at(3).at(0) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_41 = " << AEL_data.elastic_tensor.at(3).at(0) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(0).at(3) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_14 = " << AEL_data.elastic_tensor.at(0).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
      }
    }
    // Check symmetry relations for rhombohedral lattices
    else if(xvasp.str.bravais_lattice_type == "RHL") {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice has rhombohedral symmetry" << endl;
      aus << _AELSTR_WARNING_ + "The rhombohedral lattice type has some ambiguity in the choice of the coordinate frame which can affect the form of the elastic constants tensor" << endl;	  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // For rhombohedral lattice types, check whether lattice is of type I (space group # 155, 160, 161, 166, 167) or of type II (space group # 146 and 148) 
      // This changes the number of independent elastic constants (see Phys. Rev. B 90, 224104 (2014) and http://en.wikipedia.org/wiki/List_of_space_groups for details)  
      // OBSOLETE bool sg_manual_toggle = false;
      // OBSOLETE vector<string> sg_args;
      // OBSOLETE sg_args.push_back("dummyfiller");
      // OBSOLETE uint ael_spacegroupnumber = xvasp.str.SpaceGroup_ITC(sg_manual_toggle, sg_args);
      uint ael_spacegroupnumber = xvasp.str.SpaceGroup_ITC();
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice space group number = " << ael_spacegroupnumber << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      // OBSOLETE if((xvasp.str.spacegroupnumber == 146) || (xvasp.str.spacegroupnumber == 148)) {
      if((ael_spacegroupnumber == 146) || (ael_spacegroupnumber == 148)) {
        rhlII = true;
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Lattice is of type RHL II" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
        // OBSOLETE } else if((xvasp.str.spacegroupnumber == 155) || (xvasp.str.spacegroupnumber == 160) || (xvasp.str.spacegroupnumber == 161) || (xvasp.str.spacegroupnumber == 166) || (xvasp.str.spacegroupnumber == 167))
    } else if((ael_spacegroupnumber == 155) || (ael_spacegroupnumber == 160) || (ael_spacegroupnumber == 161) || (ael_spacegroupnumber == 166) || (ael_spacegroupnumber == 167)) {
      rhlI = true;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice is of type RHL I" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    } else {
      aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AELSTR_WARNING_ + "Lattice type = " << xvasp.str.bravais_lattice_type << " but space group # = " <<  xvasp.str.spacegroupnumber << endl;
      aus << _AELSTR_WARNING_ + "Lattice type = " << xvasp.str.bravais_lattice_type << " but space group # = " <<  ael_spacegroupnumber << endl;
      aus << _AELSTR_WARNING_ + "This is outside the range of 146, 148, 155, 160, 161, 166, 167 for rhombohedral space groups" << endl;
      aus << _AELSTR_WARNING_ + "Possible mismatch between lattice type and space group" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	  
    }
    // First, determine which direction is the independent one
    if((AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) < upperbound) && (AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) > lowerbound)) {
      inddir = 2;
    }	else if((AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(0).at(0) < upperbound) && (AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(0).at(0) > lowerbound)) {
      inddir = 1;
    }	else if((AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(1).at(1) < upperbound) && (AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(1).at(1) > lowerbound)) {
      inddir = 0;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_11: no two diagonal elements are equal for this material" << endl;
      aus << _AELSTR_WARNING_ + "C_22 / C_11 = " << AEL_data.elastic_tensor.at(1).at(1) / AEL_data.elastic_tensor.at(0).at(0) << endl;
      aus << _AELSTR_WARNING_ + "C_33 / C_11 = " << AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(0).at(0) << endl;
      aus << _AELSTR_WARNING_ + "C_33 / C_22 = " << AEL_data.elastic_tensor.at(2).at(2) / AEL_data.elastic_tensor.at(1).at(1) << endl;
      aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
      aus << _AELSTR_WARNING_ + "Cannot identify independent elastic tensor direction" << endl;
      aus << _AELSTR_WARNING_ + "Taking z-direction as independent direction by default" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      inddir = 2;
    }     
    // Next check off-diagonal elements in upper left quadrant and diagonal elements in lower right quadrant
    if(inddir == 2) {
      if((AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) > upperbound) || (AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
        aus << _AELSTR_WARNING_ + "C_23 / C_13 = " << AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if((AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) > upperbound) || (AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_31" << endl;
        aus << _AELSTR_WARNING_ + "C_32 / C_31 = " << AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(2).at(0) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if((AEL_data.elastic_tensor.at(4).at(4) / AEL_data.elastic_tensor.at(3).at(3) > upperbound) || (AEL_data.elastic_tensor.at(4).at(4) / AEL_data.elastic_tensor.at(3).at(3) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
        aus << _AELSTR_WARNING_ + "C_55 / C_44 = " << AEL_data.elastic_tensor.at(4).at(4) / AEL_data.elastic_tensor.at(3).at(3) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else if(inddir == 1) { 
      if((AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(1).at(0) > upperbound) || (AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(1).at(0) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
        aus << _AELSTR_WARNING_ + "C_23 / C_21 = " << AEL_data.elastic_tensor.at(1).at(2) / AEL_data.elastic_tensor.at(0).at(2) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if((AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(0).at(1) > upperbound) || (AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(0).at(1) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_31" << endl;
        aus << _AELSTR_WARNING_ + "C_32 / C_12 = " << AEL_data.elastic_tensor.at(2).at(1) / AEL_data.elastic_tensor.at(0).at(1) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if((AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(3).at(3) > upperbound) || (AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(3).at(3) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
        aus << _AELSTR_WARNING_ + "C_66 / C_44 = " << AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(3).at(3) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else if(inddir == 0) { 
      if((AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(2).at(0) > upperbound) || (AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(2).at(0) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_31" << endl;
        aus << _AELSTR_WARNING_ + "C_21 / C_31 = " << AEL_data.elastic_tensor.at(1).at(0) / AEL_data.elastic_tensor.at(2).at(0) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if((AEL_data.elastic_tensor.at(0).at(1) / AEL_data.elastic_tensor.at(0).at(2) > upperbound) || (AEL_data.elastic_tensor.at(0).at(1) / AEL_data.elastic_tensor.at(0).at(2) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_13" << endl;
        aus << _AELSTR_WARNING_ + "C_12 / C_13 = " << AEL_data.elastic_tensor.at(0).at(1) / AEL_data.elastic_tensor.at(0).at(2) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }	
      if((AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(4).at(4) > upperbound) || (AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(4).at(4) < lowerbound)) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_44" << endl;
        aus << _AELSTR_WARNING_ + "C_66 / C_55 = " << AEL_data.elastic_tensor.at(5).at(5) / AEL_data.elastic_tensor.at(4).at(4) << endl;
        aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // Next check off-diagonal elements 
    if(rhlI) {
      if(inddir == 2) {
        if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(3)) > 1.0) && (aurostd::abs(AEL_data.elastic_tensor.at(0).at(4)) < 1.0)) {
          if((AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(1).at(3)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(1).at(3)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_14 / C_24 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(1).at(3) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (AEL_data.elastic_tensor.at(4).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_14 / C_56 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(1)) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(1)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_41" << endl;
            aus << _AELSTR_WARNING_ + "C_41 / C_42 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(3).at(1) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_41" << endl;
            aus << _AELSTR_WARNING_ + "C_41 / C_65 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint i = 4; i < 6; i++) {
            for (uint j = 0; j < 4; j++) {
              if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }	  
            }
          }
          if(AEL_data.elastic_tensor.at(3).at(2) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_43 = " << AEL_data.elastic_tensor.at(3).at(2) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(2).at(3) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_34 = " << AEL_data.elastic_tensor.at(2).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(4)) > 1.0) && (aurostd::abs(AEL_data.elastic_tensor.at(0).at(3)) < 1.0)) {
          if((AEL_data.elastic_tensor.at(0).at(4) / (-AEL_data.elastic_tensor.at(1).at(4)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(4) / (-AEL_data.elastic_tensor.at(1).at(4)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_15 / C_25 = " << AEL_data.elastic_tensor.at(0).at(4) / AEL_data.elastic_tensor.at(1).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(0).at(4) / AEL_data.elastic_tensor.at(3).at(5) > upperbound) || (AEL_data.elastic_tensor.at(0).at(4) / (AEL_data.elastic_tensor.at(3).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_15 / C_46 = " << AEL_data.elastic_tensor.at(0).at(4) / AEL_data.elastic_tensor.at(3).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(4).at(0) / (-AEL_data.elastic_tensor.at(4).at(1)) > upperbound) || (AEL_data.elastic_tensor.at(4).at(0) / (-AEL_data.elastic_tensor.at(4).at(1)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_51 / C_52 = " << AEL_data.elastic_tensor.at(4).at(0) / AEL_data.elastic_tensor.at(4).at(1) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((AEL_data.elastic_tensor.at(4).at(0) / AEL_data.elastic_tensor.at(5).at(3) > upperbound) || (AEL_data.elastic_tensor.at(4).at(0) / AEL_data.elastic_tensor.at(5).at(3) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_51 / C_64 = " << AEL_data.elastic_tensor.at(4).at(0) / AEL_data.elastic_tensor.at(5).at(3) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_4" << j+1 << " = " << AEL_data.elastic_tensor.at(3).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(3)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "4 = " << AEL_data.elastic_tensor.at(j).at(3) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_6" << j+1 << " = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(5)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "6 = " << AEL_data.elastic_tensor.at(j).at(5) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          for (uint j = 2; j < 4; j++) {
            if(AEL_data.elastic_tensor.at(4).at(j) > 1.0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_5" << j+1 << " = " << AEL_data.elastic_tensor.at(4).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
            if(AEL_data.elastic_tensor.at(j).at(4) > 1.0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "5 = " << AEL_data.elastic_tensor.at(j).at(4) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          if(AEL_data.elastic_tensor.at(5).at(4) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_65 = " << AEL_data.elastic_tensor.at(5).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(4).at(5) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_56 = " << AEL_data.elastic_tensor.at(4).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: asymmetry in elastic tensor: increase strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else if(inddir == 1) {
        if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(3)) > 1.0) && (aurostd::abs(AEL_data.elastic_tensor.at(0).at(5)) < 1.0)) {
          if((AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(2).at(3)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(2).at(3)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_14 / C_34 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(2).at(3) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (AEL_data.elastic_tensor.at(4).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_14 / C_56 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(2)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_41 / C_43 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(3).at(2) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_41 / C_65 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint i = 4; i < 6; i++) {
            for (uint j = 0; j < 4; j++) {
              if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }	  
            }
          }
          if(AEL_data.elastic_tensor.at(3).at(1) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_42 = " << AEL_data.elastic_tensor.at(3).at(1) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(1).at(3) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_24 = " << AEL_data.elastic_tensor.at(1).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(5)) > 1.0) && (aurostd::abs(AEL_data.elastic_tensor.at(0).at(3)) < 1.0)) {
          if((AEL_data.elastic_tensor.at(0).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_16 / C_36 = " << AEL_data.elastic_tensor.at(0).at(5) / AEL_data.elastic_tensor.at(2).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(0).at(5) / AEL_data.elastic_tensor.at(3).at(4) > upperbound) || (AEL_data.elastic_tensor.at(0).at(5) / (AEL_data.elastic_tensor.at(3).at(4)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_16 / C_45 = " << AEL_data.elastic_tensor.at(0).at(5) / AEL_data.elastic_tensor.at(3).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(5).at(0) / (-AEL_data.elastic_tensor.at(5).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(5).at(0) / (-AEL_data.elastic_tensor.at(5).at(2)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_61 / C_62 = " << AEL_data.elastic_tensor.at(5).at(0) / AEL_data.elastic_tensor.at(5).at(1) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((AEL_data.elastic_tensor.at(5).at(0) / AEL_data.elastic_tensor.at(4).at(3) > upperbound) || (AEL_data.elastic_tensor.at(5).at(0) / AEL_data.elastic_tensor.at(4).at(3) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_61 / C_54 = " << AEL_data.elastic_tensor.at(5).at(0) / AEL_data.elastic_tensor.at(4).at(3) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint i = 3; i < 5; i++) {
            for (uint j = 0; j < 3; j++) {
              if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }	  
            }
          }
          if(AEL_data.elastic_tensor.at(5).at(1) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_62 = " << AEL_data.elastic_tensor.at(5).at(1) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(1).at(5) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_26 = " << AEL_data.elastic_tensor.at(1).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(5).at(3) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_64 = " << AEL_data.elastic_tensor.at(5).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(3).at(5) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_46 = " << AEL_data.elastic_tensor.at(3).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: asymmetry in elastic tensor: increase strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }  
      } else if(inddir == 0) {
        if((aurostd::abs(AEL_data.elastic_tensor.at(1).at(4)) > 1.0) && (aurostd::abs(AEL_data.elastic_tensor.at(1).at(5)) < 1.0)) {
          if((AEL_data.elastic_tensor.at(1).at(4) / (-AEL_data.elastic_tensor.at(2).at(4)) > upperbound) || (AEL_data.elastic_tensor.at(1).at(4) / (-AEL_data.elastic_tensor.at(2).at(4)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_25 / C_35 = " << AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(2).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(3).at(5) > upperbound) || (AEL_data.elastic_tensor.at(1).at(4) / (AEL_data.elastic_tensor.at(3).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_25 / C_46 = " << AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(3).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(4).at(1) / (-AEL_data.elastic_tensor.at(4).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(4).at(1) / (-AEL_data.elastic_tensor.at(4).at(2)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_52 / C_53 = " << AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(4).at(2) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) > upperbound) || (AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_52 / C_64 = " << AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_3" << j+1 << " = " << AEL_data.elastic_tensor.at(3).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(3)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "3 = " << AEL_data.elastic_tensor.at(j).at(3) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          for (uint j = 0; j < 3; j++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(j)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_5" << j+1 << " = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(5)) > 1.0) ) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "5 = " << AEL_data.elastic_tensor.at(j).at(5) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(0)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_51 = " << AEL_data.elastic_tensor.at(4).at(0) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(4)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_15 = " << AEL_data.elastic_tensor.at(0).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(3)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_54 = " << AEL_data.elastic_tensor.at(4).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(4)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_45 = " << AEL_data.elastic_tensor.at(3).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(4)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_65 = " << AEL_data.elastic_tensor.at(5).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(5)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_56 = " << AEL_data.elastic_tensor.at(4).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        } else if((aurostd::abs(AEL_data.elastic_tensor.at(1).at(5)) > 1.0) && (aurostd::abs(AEL_data.elastic_tensor.at(1).at(4)) < 1.0)) {
          if((AEL_data.elastic_tensor.at(1).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) > upperbound) || (AEL_data.elastic_tensor.at(1).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_26 / C_36 = " << AEL_data.elastic_tensor.at(1).at(5) / AEL_data.elastic_tensor.at(2).at(5) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(1).at(5) / AEL_data.elastic_tensor.at(3).at(4) > upperbound) || (AEL_data.elastic_tensor.at(1).at(5) / (AEL_data.elastic_tensor.at(3).at(4)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_26 / C_45 = " << AEL_data.elastic_tensor.at(0).at(5) / AEL_data.elastic_tensor.at(3).at(4) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	
          if((AEL_data.elastic_tensor.at(5).at(1) / (-AEL_data.elastic_tensor.at(5).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(5).at(1) / (-AEL_data.elastic_tensor.at(5).at(2)) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_62 / C_63 = " << AEL_data.elastic_tensor.at(5).at(1) / AEL_data.elastic_tensor.at(5).at(2) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((AEL_data.elastic_tensor.at(5).at(1) / AEL_data.elastic_tensor.at(4).at(3) > upperbound) || (AEL_data.elastic_tensor.at(5).at(1) / AEL_data.elastic_tensor.at(4).at(3) < lowerbound)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
            aus << _AELSTR_WARNING_ + "C_62 / C_54 = " << AEL_data.elastic_tensor.at(5).at(1) / AEL_data.elastic_tensor.at(4).at(3) << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          for (uint i = 3; i < 5; i++) {
            for (uint j = 0; j < 3; j++) {
              if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
                aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }	  
            }
          }
          if(AEL_data.elastic_tensor.at(5).at(0) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_61 = " << AEL_data.elastic_tensor.at(5).at(0) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if(AEL_data.elastic_tensor.at(0).at(5) > 1.0) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_16 = " << AEL_data.elastic_tensor.at(0).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          for (uint j = 3; j < 5; j++) {
            if(AEL_data.elastic_tensor.at(5).at(j) > 1.0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_6" << j+1 << " = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
            if(AEL_data.elastic_tensor.at(j).at(5) > 1.0) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << "6 = " << AEL_data.elastic_tensor.at(j).at(5) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }	  
          }
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: asymmetry in elastic tensor: increase strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }  	  
      }
    } else if(rhlII) {	  
      if(inddir == 2) {
        if((AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(1).at(3)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(1).at(3)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_14 / C_24 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(1).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (AEL_data.elastic_tensor.at(4).at(5)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_14 / C_56 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(1)) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(1)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_41 / C_42 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(3).at(1) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_41 / C_65 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(0).at(4) / (-AEL_data.elastic_tensor.at(1).at(4)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(4) / (-AEL_data.elastic_tensor.at(1).at(4)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_15 / C_25 = " << AEL_data.elastic_tensor.at(0).at(4) / AEL_data.elastic_tensor.at(1).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(3).at(5) > upperbound) || (AEL_data.elastic_tensor.at(1).at(4) / (AEL_data.elastic_tensor.at(3).at(5)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_25 / C_46 = " << AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(3).at(5) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(4).at(0) / (-AEL_data.elastic_tensor.at(4).at(1)) > upperbound) || (AEL_data.elastic_tensor.at(4).at(0) / (-AEL_data.elastic_tensor.at(4).at(1)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_51 / C_52 = " << AEL_data.elastic_tensor.at(4).at(0) / AEL_data.elastic_tensor.at(4).at(1) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) > upperbound) || (AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_51 / C_64 = " << AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        for (uint j = 2; j < 4; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(j)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_5" << j+1 << " = " << AEL_data.elastic_tensor.at(4).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(4)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << "5 = " << AEL_data.elastic_tensor.at(j).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
        for (uint j = 0; j < 3; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(j)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_6" << j+1 << " = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(5)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << "6 = " << AEL_data.elastic_tensor.at(j).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }	  
        if(AEL_data.elastic_tensor.at(3).at(2) > 1.0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_43 = " << AEL_data.elastic_tensor.at(3).at(2) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if(AEL_data.elastic_tensor.at(2).at(3) > 1.0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_34 = " << AEL_data.elastic_tensor.at(2).at(3) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
      } else if(inddir == 1) {
        if((AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(2).at(3)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (-AEL_data.elastic_tensor.at(2).at(3)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_14 / C_24 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(2).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) > upperbound) || (AEL_data.elastic_tensor.at(0).at(3) / (AEL_data.elastic_tensor.at(4).at(5)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_14 / C_56 = " << AEL_data.elastic_tensor.at(0).at(3) / AEL_data.elastic_tensor.at(4).at(5) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / (-AEL_data.elastic_tensor.at(3).at(2)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_41 / C_42 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(3).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) > upperbound) || (AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_41 / C_65 = " << AEL_data.elastic_tensor.at(3).at(0) / AEL_data.elastic_tensor.at(5).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(0).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) > upperbound) || (AEL_data.elastic_tensor.at(0).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_16 / C_36 = " << AEL_data.elastic_tensor.at(0).at(5) / AEL_data.elastic_tensor.at(2).at(5) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(2).at(5) / AEL_data.elastic_tensor.at(3).at(4) > upperbound) || (AEL_data.elastic_tensor.at(2).at(5) / (AEL_data.elastic_tensor.at(3).at(4)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_36 / C_45 = " << AEL_data.elastic_tensor.at(2).at(5) / AEL_data.elastic_tensor.at(3).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(5).at(0) / (-AEL_data.elastic_tensor.at(5).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(5).at(0) / (-AEL_data.elastic_tensor.at(5).at(2)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_61 / C_63 = " << AEL_data.elastic_tensor.at(5).at(0) / AEL_data.elastic_tensor.at(5).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(5).at(2) / AEL_data.elastic_tensor.at(4).at(3) > upperbound) || (AEL_data.elastic_tensor.at(5).at(2) / AEL_data.elastic_tensor.at(4).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_63 / C_54 = " << AEL_data.elastic_tensor.at(5).at(2) / AEL_data.elastic_tensor.at(4).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        for (uint j = 0; j < 3; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(j)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_5" << j+1 << " = " << AEL_data.elastic_tensor.at(4).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(4)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << "5 = " << AEL_data.elastic_tensor.at(j).at(4) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(1)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_62 = " << AEL_data.elastic_tensor.at(5).at(1) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(1).at(5)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_26 = " << AEL_data.elastic_tensor.at(1).at(5) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(3)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_64 = " << AEL_data.elastic_tensor.at(5).at(3) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(5)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_46 = " << AEL_data.elastic_tensor.at(3).at(5) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if(AEL_data.elastic_tensor.at(3).at(2) > 1.0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_43 = " << AEL_data.elastic_tensor.at(3).at(2) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if(AEL_data.elastic_tensor.at(2).at(3) > 1.0) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_34 = " << AEL_data.elastic_tensor.at(2).at(3) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else if(inddir == 0) {
        if((AEL_data.elastic_tensor.at(1).at(4) / (-AEL_data.elastic_tensor.at(2).at(4)) > upperbound) || (AEL_data.elastic_tensor.at(1).at(4) / (-AEL_data.elastic_tensor.at(2).at(4)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_25 / C_35 = " << AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(2).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(3).at(5) > upperbound) || (AEL_data.elastic_tensor.at(1).at(4) / (AEL_data.elastic_tensor.at(3).at(5)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_25 / C_46 = " << AEL_data.elastic_tensor.at(1).at(4) / AEL_data.elastic_tensor.at(3).at(5) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(4).at(1) / (-AEL_data.elastic_tensor.at(4).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(4).at(1) / (-AEL_data.elastic_tensor.at(4).at(2)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_52 / C_53 = " << AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(4).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) > upperbound) || (AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_14" << endl;
          aus << _AELSTR_WARNING_ + "C_52 / C_64 = " << AEL_data.elastic_tensor.at(4).at(1) / AEL_data.elastic_tensor.at(5).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(1).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) > upperbound) || (AEL_data.elastic_tensor.at(1).at(5) / (-AEL_data.elastic_tensor.at(2).at(5)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_26 / C_36 = " << AEL_data.elastic_tensor.at(1).at(5) / AEL_data.elastic_tensor.at(2).at(5) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(2).at(5) / AEL_data.elastic_tensor.at(3).at(4) > upperbound) || (AEL_data.elastic_tensor.at(2).at(5) / (AEL_data.elastic_tensor.at(3).at(4)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_36 / C_45 = " << AEL_data.elastic_tensor.at(2).at(5) / AEL_data.elastic_tensor.at(3).at(4) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	
        if((AEL_data.elastic_tensor.at(5).at(1) / (-AEL_data.elastic_tensor.at(5).at(2)) > upperbound) || (AEL_data.elastic_tensor.at(5).at(1) / (-AEL_data.elastic_tensor.at(5).at(2)) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_62 / C_63 = " << AEL_data.elastic_tensor.at(5).at(1) / AEL_data.elastic_tensor.at(5).at(2) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((AEL_data.elastic_tensor.at(5).at(2) / AEL_data.elastic_tensor.at(4).at(3) > upperbound) || (AEL_data.elastic_tensor.at(5).at(2) / AEL_data.elastic_tensor.at(4).at(3) < lowerbound)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant C_15" << endl;
          aus << _AELSTR_WARNING_ + "C_63 / C_54 = " << AEL_data.elastic_tensor.at(5).at(2) / AEL_data.elastic_tensor.at(4).at(3) << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        for (uint j = 0; j < 3; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(j)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_4" << j+1 << " = " << AEL_data.elastic_tensor.at(3).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(3)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << "4 = " << AEL_data.elastic_tensor.at(j).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(0)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_51 = " << AEL_data.elastic_tensor.at(4).at(0) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(4)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_15 = " << AEL_data.elastic_tensor.at(0).at(4) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(0)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_61 = " << AEL_data.elastic_tensor.at(5).at(0) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(0).at(5)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_16 = " << AEL_data.elastic_tensor.at(0).at(5) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  

        if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(4)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_65 = " << AEL_data.elastic_tensor.at(5).at(4) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(5)) > 1.0) ) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_56 = " << AEL_data.elastic_tensor.at(4).at(5) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
      }
    }
    }  
    // Check symmetry relations for orthorhombic lattices
    else if((xvasp.str.bravais_lattice_type == "ORC") || (xvasp.str.bravais_lattice_type == "ORCC") || (xvasp.str.bravais_lattice_type == "ORCF") || (xvasp.str.bravais_lattice_type == "ORCI")) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice has orthorhombic symmetry" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
      // Check off-diagonal elements in lower right quadrant are zero     
      for (uint i = 3; i < 6; i++) {
        for (uint j = 3; j < 6; j++) {
          if((i != j) && (aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
      }
      // Check upper right and lower left quadrants are zero
      for (uint i = 3; i < 6; i++) {
        for (uint j = 0; j < 3; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0) ) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
      }
    }
    // Check symmetry relations for monoclinic lattices
    else if((xvasp.str.bravais_lattice_type == "MCL") || (xvasp.str.bravais_lattice_type == "MCLC")) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice has monoclinic symmetry" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);       
      if(AEL_data.elastic_tensor.at(4).at(0) > 1.0) {
        for (uint j = 0; j < 3; j++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(j)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_4" << j+1 << " = " << AEL_data.elastic_tensor.at(3).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(j)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_6" << j+1 << " = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(3)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << "4 = " << AEL_data.elastic_tensor.at(j).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(5)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << j+1 << "6 = " << AEL_data.elastic_tensor.at(5).at(j) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
        if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(3)) > 1.0)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_54 = " << AEL_data.elastic_tensor.at(4).at(3) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if((aurostd::abs(AEL_data.elastic_tensor.at(4).at(5)) > 1.0)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_56 = " << AEL_data.elastic_tensor.at(4).at(5) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(4)) > 1.0)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_45 = " << AEL_data.elastic_tensor.at(3).at(4) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(4)) > 1.0)) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "C_65 = " << AEL_data.elastic_tensor.at(5).at(4) << " != 0.0" << endl;
          aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
      } else if(AEL_data.elastic_tensor.at(3).at(0) > 1.0) {
        for (uint j = 0; j < 3; j++) {
          for (uint i = 4; i < 6; i++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0)) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0)) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
          }
        }
        for (uint i = 4; i < 6; i++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(3).at(i)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_4" << i+1 << " = " << AEL_data.elastic_tensor.at(3).at(i) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(3)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << "4 = " << AEL_data.elastic_tensor.at(i).at(3) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
      } else if(AEL_data.elastic_tensor.at(5).at(0) > 1.0) {
        for (uint j = 0; j < 3; j++) {
          for (uint i = 3; i < 5; i++) {
            if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(j)) > 1.0)) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << i+1 << j+1 << " = " << AEL_data.elastic_tensor.at(i).at(j) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in normal stress: Increase normal strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            if((aurostd::abs(AEL_data.elastic_tensor.at(j).at(i)) > 1.0)) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_WARNING_ + "C_" << j+1 << i+1 << " = " << AEL_data.elastic_tensor.at(j).at(i) << " != 0.0" << endl;
              aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
          }
        }
        for (uint i = 3; i < 5; i++) {
          if((aurostd::abs(AEL_data.elastic_tensor.at(5).at(i)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_6" << i+1 << " = " << AEL_data.elastic_tensor.at(5).at(i) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
          if((aurostd::abs(AEL_data.elastic_tensor.at(i).at(5)) > 1.0)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "C_" << i+1 << "6 = " << AEL_data.elastic_tensor.at(i).at(5) << " != 0.0" << endl;
            aus << _AELSTR_WARNING_ + "Inconsistency in elastic constant: Noise in shear stress: Increase shear strain" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }	  
        }
      }
    }
    // Checks if lattice is triclinic - no structure other than symmetry about the diagonal
    else if(xvasp.str.bravais_lattice_type == "TRI") {    
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Lattice has triclinic symmetry" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);             
    }
    // Warns user that lattice symmetry type was not found
    else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "Lattice symmetry type not found!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    }  
    return 0;
  }
} // namespace AEL_functions

// ***************************************************************************
// AEL_functions::elasticmoduli
// ***************************************************************************
namespace AEL_functions {
  //
  // elasticmoduli: Calculate elastic moduli from elastic constants
  // See Phys. Rev. Materials 1, 015401 (2017) and Scientific Data 2, 150009 (2015) for details
  // Also calculates average calculated external pressure in strained structures
  //
  uint elasticmoduli(_AEL_data& AEL_data, ofstream& FileMESSAGE) {
    ostringstream aus;
    double sumpressure = 0.0;
    double dnpressure = 0.0;
    // First calculate average external pressure
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Calculating elastic average external pressure" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    for (uint i = 0; i < AEL_data.pressurecalculated.size(); i++) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Pressure = " << AEL_data.pressurecalculated.at(i) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      sumpressure = sumpressure + AEL_data.pressurecalculated.at(i);
      dnpressure = dnpressure + 1.0;
    }
    // Take average and convert from kB to GPa
    AEL_data.average_external_pressure = 0.1 * (sumpressure / dnpressure);
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Average external pressure = " << AEL_data.average_external_pressure << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // Calculate elastic moduli   
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Calculating elastic moduli" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    AEL_data.bulkmodulus_voigt = ((AEL_data.elastic_tensor.at(0).at(0) + AEL_data.elastic_tensor.at(1).at(1) + AEL_data.elastic_tensor.at(2).at(2)) + 2.0 * (AEL_data.elastic_tensor.at(0).at(1) + AEL_data.elastic_tensor.at(1).at(2) + AEL_data.elastic_tensor.at(2).at(0))) / 9.0; 
    AEL_data.bulkmodulus_reuss = 1.0 / ((AEL_data.compliance_tensor.at(0).at(0) + AEL_data.compliance_tensor.at(1).at(1) + AEL_data.compliance_tensor.at(2).at(2)) + 2.0 * (AEL_data.compliance_tensor.at(0).at(1) + AEL_data.compliance_tensor.at(1).at(2) + AEL_data.compliance_tensor.at(2).at(0))); 
    AEL_data.shearmodulus_voigt = ((AEL_data.elastic_tensor.at(0).at(0) + AEL_data.elastic_tensor.at(1).at(1) + AEL_data.elastic_tensor.at(2).at(2)) - (AEL_data.elastic_tensor.at(0).at(1) + AEL_data.elastic_tensor.at(1).at(2) + AEL_data.elastic_tensor.at(2).at(0)) + 3.0 * (AEL_data.elastic_tensor.at(3).at(3) + AEL_data.elastic_tensor.at(4).at(4) + AEL_data.elastic_tensor.at(5).at(5))) / 15.0; 
    AEL_data.shearmodulus_reuss = 15.0 / (4.0 * (AEL_data.compliance_tensor.at(0).at(0) + AEL_data.compliance_tensor.at(1).at(1) + AEL_data.compliance_tensor.at(2).at(2)) - 4.0 * (AEL_data.compliance_tensor.at(0).at(1) + AEL_data.compliance_tensor.at(1).at(2) + AEL_data.compliance_tensor.at(2).at(0)) + 3.0 * (AEL_data.compliance_tensor.at(3).at(3) + AEL_data.compliance_tensor.at(4).at(4) + AEL_data.compliance_tensor.at(5).at(5))); 
    AEL_data.bulkmodulus_vrh = (AEL_data.bulkmodulus_voigt + AEL_data.bulkmodulus_reuss) / 2.0;
    AEL_data.shearmodulus_vrh = (AEL_data.shearmodulus_voigt + AEL_data.shearmodulus_reuss) / 2.0;
    AEL_data.elastic_anisotropy = 5.0 * (AEL_data.shearmodulus_voigt / AEL_data.shearmodulus_reuss) + (AEL_data.bulkmodulus_voigt / AEL_data.bulkmodulus_reuss) - 6.0;
    if(AEL_data.elastic_anisotropy < 0.0) {
      AEL_data.elastic_anisotropy = 0.0;
    }
    AEL_data.poisson_ratio = (3.0 * AEL_data.bulkmodulus_vrh - 2.0 * AEL_data.shearmodulus_vrh) / (6.0 * AEL_data.bulkmodulus_vrh + 2.0 * AEL_data.shearmodulus_vrh);
    AEL_data.bulkmodulus_voigt_half = ((AEL_data.elastic_tensor_half.at(0).at(0) + AEL_data.elastic_tensor_half.at(1).at(1) + AEL_data.elastic_tensor_half.at(2).at(2)) + 2.0 * (AEL_data.elastic_tensor_half.at(0).at(1) + AEL_data.elastic_tensor_half.at(1).at(2) + AEL_data.elastic_tensor_half.at(2).at(0))) / 9.0; 
    AEL_data.bulkmodulus_reuss_half = 1.0 / ((AEL_data.compliance_tensor_half.at(0).at(0) + AEL_data.compliance_tensor_half.at(1).at(1) + AEL_data.compliance_tensor_half.at(2).at(2)) + 2.0 * (AEL_data.compliance_tensor_half.at(0).at(1) + AEL_data.compliance_tensor_half.at(1).at(2) + AEL_data.compliance_tensor_half.at(2).at(0))); 
    AEL_data.shearmodulus_voigt_half = ((AEL_data.elastic_tensor_half.at(0).at(0) + AEL_data.elastic_tensor_half.at(1).at(1) + AEL_data.elastic_tensor_half.at(2).at(2)) - (AEL_data.elastic_tensor_half.at(0).at(1) + AEL_data.elastic_tensor_half.at(1).at(2) + AEL_data.elastic_tensor_half.at(2).at(0)) + 3.0 * (AEL_data.elastic_tensor_half.at(3).at(3) + AEL_data.elastic_tensor_half.at(4).at(4) + AEL_data.elastic_tensor_half.at(5).at(5))) / 15.0; 
    AEL_data.shearmodulus_reuss_half = 15.0 / (4.0 * (AEL_data.compliance_tensor_half.at(0).at(0) + AEL_data.compliance_tensor_half.at(1).at(1) + AEL_data.compliance_tensor_half.at(2).at(2)) - 4.0 * (AEL_data.compliance_tensor_half.at(0).at(1) + AEL_data.compliance_tensor_half.at(1).at(2) + AEL_data.compliance_tensor_half.at(2).at(0)) + 3.0 * (AEL_data.compliance_tensor_half.at(3).at(3) + AEL_data.compliance_tensor_half.at(4).at(4) + AEL_data.compliance_tensor_half.at(5).at(5))); 
    AEL_data.bulkmodulus_vrh_half = (AEL_data.bulkmodulus_voigt_half + AEL_data.bulkmodulus_reuss_half) / 2.0;
    AEL_data.shearmodulus_vrh_half = (AEL_data.shearmodulus_voigt_half + AEL_data.shearmodulus_reuss_half) / 2.0;
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "Elastic anisotropy = " << AEL_data.elastic_anisotropy << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // Calculates Young's modulus along principle strain directions
    // This expression is only correct for elastically isotropic (i.e. cubic) systems 
    for (uint i = 0; i < 3; i++) {
      AEL_data.youngsmodulus_directional.push_back(1.0 / AEL_data.compliance_tensor.at(i).at(i));
    }
    // Calculates shear modulus for different face/strain combinations
    for (uint i = 3; i < 6; i++) {
      AEL_data.shearmodulus_directional.push_back(1.0 / AEL_data.compliance_tensor.at(i).at(i));
    }
    // Calculates Poisson ratio for different stress/strain combinations
    for (uint i = 0; i < 3; i++) {
      for (uint j = 0; j < 3; j++) {
        if(i != j) {
          AEL_data.poissonratio_directional.push_back(- AEL_data.compliance_tensor.at(j).at(i) * AEL_data.youngsmodulus_directional.at(i));
        }
      }
    }
    // Calculates Young's modulus using expression E = 9BG/(3B+G)
    AEL_data.youngsmodulus_vrh = (9.0 * AEL_data.bulkmodulus_vrh * AEL_data.shearmodulus_vrh) / ((3.0 * AEL_data.bulkmodulus_vrh) + AEL_data.shearmodulus_vrh);
    // Calculates tranverse, longitudinal and average speed of sound
    // See J.-P. Poirer, Introduction to the Physics of the Earth's Interior
    // First convert elastic moduli from GPa to Pa
    double shearmod_Pa = AEL_data.shearmodulus_vrh * (pow(10, 9));
    double bulkmod_Pa = AEL_data.shearmodulus_vrh * (pow(10, 9));
    // AEL_data.speed_sound_transverse = pow(AEL_data.shearmodulus_vrh / AEL_data.mass_density, 0.5);
    // AEL_data.speed_sound_longitudinal = pow((AEL_data.bulkmodulus_vrh + ((4.0 / 3.0) * AEL_data.shearmodulus_vrh)) / AEL_data.mass_density, 0.5);
    AEL_data.speed_sound_transverse = pow(shearmod_Pa / AEL_data.mass_density, 0.5);
    AEL_data.speed_sound_longitudinal = pow((bulkmod_Pa + ((4.0 / 3.0) * shearmod_Pa)) / AEL_data.mass_density, 0.5);
    double vtcubed = pow(AEL_data.speed_sound_transverse, 3.0);
    double vlcubed = pow(AEL_data.speed_sound_longitudinal, 3.0);
    AEL_data.speed_sound_average = pow(((2.0 / vtcubed) + (1.0 / vlcubed)) / 3.0, -1.0 / 3.0);
    // Calculates Pugh's modulus ratio
    // Low values imply high ductility and malleability; high values imply high hardness and brittleness
    // See Philosophical Magazine Series 7, Vol. 45, page 367; PRB 84, 121405 and Intermetallics 19, 1275
    AEL_data.pughs_modulus_ratio = AEL_data.shearmodulus_vrh / AEL_data.bulkmodulus_vrh;
    // Calculates Debye temperature based on speed of sound: see J.-P. Poirer, Introduction to the Physics of the Earth's Interior
    // OBSOLETE double hbar = PLANKSCONSTANT_h / (2 * PI);
    double hbar = PLANCKSCONSTANT_hbar;  //ME20181020
    double sqrtvol = pow(AEL_data.cellvolumem3, 0.5);
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ + "hbar = " << hbar << endl;
    aus << _AELSTR_MESSAGE_ + "KBOLTZ = " << KBOLTZ << endl;
    aus << _AELSTR_MESSAGE_ + "sqrtvol = " << sqrtvol << endl;
    aus << _AELSTR_MESSAGE_ + "PI = " << PI << endl;
    aus << _AELSTR_MESSAGE_ + "No. atoms / cell = " << AEL_data.natomscell << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    AEL_data.debye_temperature = (hbar / KBOLTZ) * (pow(((6 * PI * PI * AEL_data.natomscell) / AEL_data.cellvolumem3), 1.0/3.0)) * AEL_data.speed_sound_average;
    return 0;
  }
} // namespace AEL_functions

// ***************************************************************************************
//  End of AFLOW AEL functions to fit stress-strain and calculate elastic properties
// ***************************************************************************************

#endif // _AFLOW_AEL_ELASTIC_FIT_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
