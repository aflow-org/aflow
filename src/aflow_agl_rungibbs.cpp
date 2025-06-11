// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_RUNGIBBS_CPP
#define _AFLOW_AGL_RUNGIBBS_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"

#define LIB2RAW_FIX_THIS_LATER TRUE // (XHOST.vflag_control.flag("BEEP"))

// ###############################################################################
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2018)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. Materials 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//

// ******************************************************************
// This set of functions implement the GIBBS algorithm
// ******************************************************************

// ***************************************************************************
// AGL_functions::gibbsrun
// ***************************************************************************
namespace AGL_functions {
  //
  // gibbsrun: runs GIBBS algorithm
  // Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  // Main function for GIBBS method
  //
  uint gibbsrun(_AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    uint imin, itry;
    uint ndatapoints, nskippedpoints;
    uint ilow;
    double energyreference, volumereference, tdebyevol, logvol, gammaalphaT, Bstatic, DebyeIntegral, DebyeIntegralerror, UIntEnergyVib, HelmholtzEnergy, Entropy, Cvt, UIntEnergy;
    double Emin, Vmin;
    double Vmin_Ang, Vmin_Bohr;
    double volume_0pressure, gibbsenergykjmol_0pressure, gibbsenergyeV_0pressure, binp_bcnt;
    double energyvalue, energyerror;
    double tdebyemin, xmin, plnv;
    double IntEnergyVibkjmol, CvVibkjmol, CvVibUnitkB, HelmholtzVibkjmol, EntropyVibkjmol, EdftUintVib;
    // OBSOLETE double IntEnergyVibmeV, HelmholtzVibmeV, EntropyVibmeV, EntropyVibUnitkB, EdftUintVibeV;
    double IntEnergyVibmeV, HelmholtzVibmeV, EntropyVibmeV, EntropyVibUnitkB;
    // OBSOLETE double gibbsenergy0pressureeV;
    double S_Entropy, Pbeta, Cpkjmol, CpUnitkB;
    // OBSOLETE double g0dhy2kjm, dzero;
    // OBSOLETE double ge0peV, dzero;
    double dzero;
    vector<int> idx;
    vector<double> xconfigurationvector, gibbsenergydatapoints, helmholtzenergydatapoints;
    vector<double> energypolynomialcoeffs, gibbsenergypolynomialcoeffs, helmholtzenergypolynomialcoeffs, tdebyepolynomialcoeffs;
    vector<double> energypolynomialerror, helmholtzpolynomialerror, gibbspolynomialerror, gibbsenergykjmol, gibbsenergyeV, relative_error, press_stat, F_Helmholtzvinp;
    vector<double> bcnt_K_sp, bcnt_P_sp, bcnt_V_sp, bcnt_beta_sp;
    vector<double> volumetemperature0pressure, gibbsenergytemperature0pressurekjmol, dbulkmodulusdpt0pressure, d2bulkmodulusdp2t0pressure;
    vector<double> theta, Cv;
    int pmerr = 0;
    ostringstream aus;
    uint aglerror = 0, minpolyerror = 0;
    // Declare and initialize structs to hold thermodynamic quantities as a function of pressure and temperature
    // Only used if AGL_data.run_all_pressure_temperature is set to true
    AGL_pressure_temperature_energies AGL_pressure_temperature_energy;
    AGL_pressure_enthalpies AGL_pressure_enthalpy;
    AGL_pressure_temperature_energy.pressure_external = 0.0;
    AGL_pressure_temperature_energy.temperature_external = 0.0;
    AGL_pressure_temperature_energy.gibbs_free_energy = 0.0;
    AGL_pressure_temperature_energy.volume_equilibrium = 0.0;
    AGL_pressure_temperature_energy.E_DFT_internal_vib_energy = 0.0;    
    AGL_pressure_temperature_energy.internal_energy = 0.0;
    AGL_pressure_temperature_energy.helmholtz_vib_energy = 0.0;    
    AGL_pressure_temperature_energy.vib_entropy = 0.0;
    AGL_pressure_enthalpy.pressure_external = 0.0;
    AGL_pressure_enthalpy.enthalpy = 0.0;
    // Variables to record thermal and elastic properties at specific pressure values
    // Only used if AGL_data.run_all_pressure_temperature is set to true
    double bulkmodpres, pres_static, d2EnergydVolume2_dynamic_pres, Gruneisen_pres, theta_pres, Cv_pres;

    if(AGL_data.i_optimize_beta != 0) {
      AGL_data.optimize_beta = true;
    }
    aurostd::StringstreamClean(AGL_data.outfiless);
    // Compute the function of the Poisson ratio which appears in the expression for the Debye temperature.
    // This function comes from averaging over the transverse and longitudinal speeds of sounds
    // See, for example and more details, Miguel A. Blanco, PhD Thesis; 
    // or J-P Poirier, Introduction to the Physics of the Earth's Interior
    double poissonratiofunc1 = (2.0/3.0) * ((1.0 + AGL_data.poissonratio)/(1.0 - 2.0 * AGL_data.poissonratio));
    double poissonratiofunc2 = (1.0/3.0) * ((1.0 + AGL_data.poissonratio)/(1.0 - AGL_data.poissonratio));
    double poissonratiofunc3 = 2.0 * sqrt(pow(poissonratiofunc1, 3)) + sqrt(pow(poissonratiofunc2, 3));
    AGL_data.poissonratiofunction = pow((3.0 / poissonratiofunc3), 1.0 / 3.0);

    // Write appropriate header for GIBBS output data for all p, T values
    if(AGL_data.i_eqn_of_state >= 0) {
      AGL_data.outfiless << "AFLOW AGL - (P,T) thermodynamics of crystals from (E,V) data" << endl;
      AGL_data.outfiless << endl;
      AGL_data.outfiless << AGL_data.sysname << endl;
      AGL_data.outfiless << "Number of data points: " << AGL_data.volumeinput.size() << endl;
      AGL_data.outfiless << endl;
    } else {
      AGL_data.outfiless << "# AFLOW AGL - (P,T) thermodynamics of crystals from (E,V) data" << endl;
      AGL_data.outfiless << "# " << endl;
      AGL_data.outfiless << "# " << AGL_data.sysname << endl;
      AGL_data.outfiless << "# " << endl;
      AGL_data.outfiless << "# P(GPa)              T(K)       V(bohr^3)     G(kJ/mol)        Bt(GPa)" << endl;
    }

    // Reads E-V data in from a file
    if(AGL_data.fittype == 5) {
      ifstream evtest;
      vector<double> energy, volume;
      string evdataname = AGL_data.dirpathname + "/evdata.in";
      evtest.open(evdataname.c_str());
      evtest >> ndatapoints;
      energy.resize(ndatapoints);
      volume.resize(ndatapoints);
      AGL_data.volumeinput.resize(ndatapoints);
      AGL_data.energyinput.resize(ndatapoints);
      for (uint i = 0; i < ndatapoints; i++) {
        evtest >> volume.at(i) >> energy.at(i);
        AGL_data.volumeinput.at(i) = volume.at(i);
        AGL_data.energyinput.at(i) = energy.at(i);
      }
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Number of data points = " << ndatapoints << endl;
      for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
        aus << _AGLSTR_MESSAGE_ << "(V, E) = " << AGL_data.volumeinput.at(i) << "\t" << AGL_data.energyinput.at(i) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // Sort data points in order of increasing volume
    if(AGL_data.i_debye != 1) {
      aglerror = AGL_functions::qcksortev (AGL_data.volumeinput, AGL_data.energyinput, FileMESSAGE);
    } else {
      aglerror = AGL_functions::qcksortevt (AGL_data.volumeinput, AGL_data.energyinput, AGL_data.tdebye, FileMESSAGE);
    }
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V) data" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 5;
    }

    // Find global minimum of energy data points
    double etref = AGL_data.energyinput.at(0);
    uint itmin = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      if(AGL_data.energyinput.at(i) < etref) {
        etref = AGL_data.energyinput.at(i);
        itmin = i;
      }
    }
    // Find midpoint of energy data points
    uint itmidpoint = (AGL_data.energyinput.size() / 2) - 1;
    // Check difference between midpoint and minimum; give warning if it is greater than 1
    AGL_data.itdiff = itmidpoint - itmin;
    if(abs(AGL_data.itdiff) > 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Global energy minimum is not at center of E(V) data" << endl;
      aus << _AGLSTR_WARNING_ + "Global energy minimum is at point " << itmin << endl;
      aus << _AGLSTR_WARNING_ + "Center of E(V) data is at point " << itmidpoint << endl;
      aus << _AGLSTR_WARNING_ + "Number of E(V) points between global energy minimum and center of E(V) data = " << AGL_data.itdiff << endl;    
      aus << _AGLSTR_WARNING_ + "Check for noise in E(V) data and / or relaxation of initial structure" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // Search for the minimum of the accepted data
    double vtref = AGL_data.volumeinput.at(itmin);

    xconfigurationvector.resize(AGL_data.volumeinput.size());
    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      xconfigurationvector.at(i) = pow(AGL_data.volumeinput.at(i)/vtref, third);
    }

    // Check that all points show concave patterns
    // First check where the concave pattern begins
    if(AGL_data.energyinput.size() != AGL_data.volumeinput.size()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error in vector sizes" << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
      aus << _AGLSTR_ERROR_ + "Error in vector sizes" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 3;
    }    
    uint evj = 1;
    while ( ( (AGL_data.energyinput.at(evj)-AGL_data.energyinput.at(evj-1))/(AGL_data.volumeinput.at(evj)-AGL_data.volumeinput.at(evj-1)) >= (AGL_data.energyinput.at(evj)-AGL_data.energyinput.at(evj+1))/(AGL_data.volumeinput.at(evj)-AGL_data.volumeinput.at(evj+1)) ) && evj < AGL_data.volumeinput.size()-2) {
      evj = evj + 1;
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Concavity check: evj = " << evj << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // If only the last three points are concave, then GIBBS will not be able to use them to fit a polynomial
    // GIBBS exits giving an error message
    if(evj == AGL_data.volumeinput.size()-2) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "All points show convex patterns" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    }
    // If evj > 1, then some E-V points do not form part of the concave pattern
    // If this is the case, then AGL gives a warning and records that there is noise in the E-V data
    // If the AGL runs fails, this information will be useful for suggesting potential corrections
    if(evj > 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Some points show convex patterns and will be removed from the data set" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      AGL_data.EV_noise = true;
    }
    AGL_data.energyinput.at(0) = AGL_data.energyinput.at(evj-1);
    AGL_data.volumeinput.at(0) = AGL_data.volumeinput.at(evj-1);
    AGL_data.energyinput.at(1) = AGL_data.energyinput.at(evj);
    AGL_data.volumeinput.at(1) = AGL_data.volumeinput.at(evj);
    AGL_data.energyinput.at(2) = AGL_data.energyinput.at(evj+1);
    AGL_data.volumeinput.at(2) = AGL_data.volumeinput.at(evj+1);
    if(AGL_data.i_debye == 1) {
      AGL_data.tdebye.at(0) = AGL_data.tdebye.at(evj-1);
      AGL_data.tdebye.at(1) = AGL_data.tdebye.at(evj);
      AGL_data.tdebye.at(2) = AGL_data.tdebye.at(evj+1);
    }
    evj = evj + 1;
    uint jtmax = 2;
    // evj marks the last accepted point, i the new trial point
    if(LDEBUG) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // OBSOLETE for (uint i=evj+1; i < AGL_data.volumeinput.size()-1; i++)
    for (uint i=evj+1; i < AGL_data.volumeinput.size(); i++) {
      if(LDEBUG) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "i = " << i << endl;
        aus << _AGLSTR_MESSAGE_ + "evj = " << evj << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(AGL_data.fittype == 0) {
        if( (AGL_data.energyinput.at(evj)-AGL_data.energyinput.at(evj-1))/(AGL_data.volumeinput.at(evj)-AGL_data.volumeinput.at(evj-1)) < (AGL_data.energyinput.at(evj)-AGL_data.energyinput.at(i))/(AGL_data.volumeinput.at(evj)-AGL_data.volumeinput.at(i)) ) {
          evj = evj + 1;
          jtmax = i;
          AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(i);
          AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(i);
          if(AGL_data.i_debye == 1) AGL_data.tdebye.at(evj) = AGL_data.tdebye.at(i);
        }
      } else if(AGL_data.fittype == 1 || AGL_data.fittype == 3) {
        if( i < itmin ) {
          if((AGL_data.energyinput.at(i) < AGL_data.energyinput.at(evj)) && (AGL_data.energyinput.at(i) > AGL_data.energyinput.at(i+1))) {
            evj = evj + 1;
            jtmax = i;
            AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(i);
            AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(i);
            if(AGL_data.i_debye == 1) AGL_data.tdebye.at(evj) = AGL_data.tdebye.at(i);
          }
        } else if(i == itmin) {
          evj = evj + 1;
          jtmax = i;
          AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(i);
          AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(i);
          if(AGL_data.i_debye == 1) AGL_data.tdebye.at(evj) = AGL_data.tdebye.at(i);
        } else if(i > itmin) {
          if((AGL_data.energyinput.at(i) > AGL_data.energyinput.at(evj)) && (AGL_data.energyinput.at(i) < AGL_data.energyinput.at(i+1))) {
            evj = evj + 1;
            jtmax = i;
            AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(i);
            AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(i);
            if(AGL_data.i_debye == 1) AGL_data.tdebye.at(evj) = AGL_data.tdebye.at(i);
          }
        }	
      } else if(AGL_data.fittype == 2 || AGL_data.fittype == 4 || AGL_data.fittype == 5) {
        evj = evj + 1;
        jtmax = i;
        AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(i);
        AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(i);
        if(AGL_data.i_debye == 1) AGL_data.tdebye.at(evj) = AGL_data.tdebye.at(i);
      } 
    }
    // If the global minimum lies outside of the range of the accepted points, a warning is given and the concavity check is restarted at the global minimum
    // This is only done for AGL_data.fittype = 0, as in the other cases the global minimum is always included
    // Problems with noise in the concavity data are often caused by insufficient k-points or basis set
    if(jtmax < itmin && (AGL_data.fittype == 0)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Global minimum of (E, V) data lies outside of initial range of points accepted for concavity" << endl;
      aus << _AGLSTR_WARNING_ + "Restarting concavity check at global minimum" << endl;
      aus << _AGLSTR_WARNING_ + "Consider checking (E, V) data and increasing the number of k-points or size of the basis set" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if( (AGL_data.energyinput.at(itmin)-AGL_data.energyinput.at(itmin-1))/(AGL_data.volumeinput.at(itmin)-AGL_data.volumeinput.at(itmin-1)) < (AGL_data.energyinput.at(itmin)-AGL_data.energyinput.at(itmin+1))/(AGL_data.volumeinput.at(itmin)-AGL_data.volumeinput.at(itmin+1)) ) {
        AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(itmin);
        AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(itmin);
        for (uint i=itmin+1; i < AGL_data.energyinput.size()-1; i++) {
          if( (AGL_data.energyinput.at(evj)-AGL_data.energyinput.at(evj-1))/(AGL_data.volumeinput.at(evj)-AGL_data.volumeinput.at(evj-1)) < (AGL_data.energyinput.at(evj)-AGL_data.energyinput.at(i))/(AGL_data.volumeinput.at(evj)-AGL_data.volumeinput.at(i)) ) {
            evj = evj + 1;
            AGL_data.energyinput.at(evj) = AGL_data.energyinput.at(i);
            AGL_data.volumeinput.at(evj) = AGL_data.volumeinput.at(i);
            if(AGL_data.i_debye == 1) AGL_data.tdebye.at(evj) = AGL_data.tdebye.at(i);
          }
        }
      }
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Concavity check: evj = " << evj << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    ndatapoints = evj + 1;
    if(ndatapoints < AGL_data.nstructsinit) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Some points showed convex patterns and have been removed from the data set" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      AGL_data.EV_noise = true;
      // OBSOLETE int nskippedpoints = AGL_data.nstructsinit - ndatapoints;
      nskippedpoints = AGL_data.nstructsinit - ndatapoints;
      if(nskippedpoints > AGL_data.skiparunsmax) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Total number of skipped points exceeds maximum limit" << endl;
        aus << _AGLSTR_ERROR_ + "Total number of skipped points = " << nskippedpoints << endl;
        aus << _AGLSTR_ERROR_ + "Maximum limit of skipped points = " << AGL_data.skiparunsmax << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = 2;
        return aglerror;
      }
    }
    AGL_data.energyinput.resize(ndatapoints);
    AGL_data.volumeinput.resize(ndatapoints);

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "(E, V) resized" << endl;
    aus << _AGLSTR_MESSAGE_ << "(E, V) data points to be used in AGL:" << endl;
    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      aus <<_AGLSTR_MESSAGE_ << "V = " << AGL_data.volumeinput.at(i) << "; E = " << AGL_data.energyinput.at(i) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // For fit types 3 and 4, fit the (E, V) data by a polynomial, and then use this polynomial to generate a new (E, V) data set
    // This can be useful for smoothing out slightly noisy (E, V) data sets, making the GIBBS algorithm more robust
    if(AGL_data.fittype == 3 || AGL_data.fittype == 4) {
      imin = 0;
      energyreference = AGL_data.energyinput.at(0);
      for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
        if(AGL_data.energyinput.at(i) < energyreference) {
          energyreference = AGL_data.energyinput.at(i);
          imin = i;
        }
      }
      if(!LDEBUG) {
        aglerror = AGL_functions::polynomial_fit_weight_ave (imin, energypolynomialcoeffs, energypolynomialerror, xconfigurationvector, AGL_data.energyinput, AGL_data, FileMESSAGE);
      } else {
        aglerror = AGL_functions::polynomial_fit_weight_ave_debug (imin, energypolynomialcoeffs, energypolynomialerror, xconfigurationvector, AGL_data.energyinput, AGL_data, FileMESSAGE);
      }
      if(aglerror != 0) {
        if(aglerror == 2) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Problem inverting matrix to fit polynomial" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 4;
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 2;
        }
      }
      AGL_data.energyinput.resize(xconfigurationvector.size());
      AGL_data.volumeinput.resize(xconfigurationvector.size());
      for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, plnv, 0);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Failed to evaluate polynomial" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 5;
        }
        AGL_data.energyinput.at(i) = plnv;
        AGL_data.volumeinput.at(i) = (pow(xconfigurationvector.at(i), 3.0)) * vtref;
      }
    }

    // Search for the minimum energy point of the accepted data points
    volumereference = AGL_data.volumeinput.at(0);
    energyreference = AGL_data.energyinput.at(0);
    imin = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      if(AGL_data.energyinput.at(i) < energyreference) {
        volumereference = AGL_data.volumeinput.at(i);
        energyreference = AGL_data.energyinput.at(i);
        imin = i;
      }
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Minimum of accepted data: imin = " << imin << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    xconfigurationvector.resize(AGL_data.volumeinput.size());
    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      xconfigurationvector.at(i) = pow(AGL_data.volumeinput.at(i) / volumereference, third);
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "x configuration vector assigned" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // If the lowest energy corresponds to the smallest or largest volume, then the minimum lies outside of the accepted input data
    // GIBBS exits giving an error message
    if(imin == 0 || imin == (AGL_data.volumeinput.size() - 1)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "static minimum outside input data" << endl;
      aus << _AGLSTR_ERROR_ + "structure with lowest energy is " << imin << endl;
      aus << _AGLSTR_ERROR_ + "energy of this structure = " << energyreference << endl;
      aus << _AGLSTR_ERROR_ + "total number of structures = " << AGL_data.volumeinput.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    }

    // Obtain the polynomial fit of E(static) vs. x
    // xconfigurationvector = V^(1/3)
    // energypolynomialcoeffs contains coefficients of fitted polynomial
    if(AGL_data.energyinput.size() != xconfigurationvector.size()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Vectors have different size" << endl;
      aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size()" << AGL_data.volumeinput.size() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 3;
    }
    if(!LDEBUG) {
      aglerror = AGL_functions::polynomial_fit_weight_ave (imin, energypolynomialcoeffs, energypolynomialerror, xconfigurationvector, AGL_data.energyinput, AGL_data, FileMESSAGE);
    } else {
      aglerror = AGL_functions::polynomial_fit_weight_ave_debug (imin, energypolynomialcoeffs, energypolynomialerror, xconfigurationvector, AGL_data.energyinput, AGL_data, FileMESSAGE);
    }
    if(LDEBUG) {
      string polyerepname = AGL_data.dirpathname + "/AGL_Plot_adjusted_energyinput.dat";
      stringstream polyerepss;
      for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
        polyerepss << xconfigurationvector.at(i) << "\t" << AGL_data.energyinput.at(i) << endl;
      }
      if(!aurostd::stringstream2file(polyerepss, polyerepname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << polyerepname.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }	
      aurostd::StringstreamClean(polyerepss);
    }
    if(aglerror != 0) {
      if(aglerror == 2) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Problem inverting matrix to fit polynomial" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = 4;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = 2;
      }
      return aglerror;
    }
    if(energypolynomialcoeffs.size() == 3) {
      energypolynomialcoeffs.push_back(0.0);
    }

    if(LDEBUG) {
      // Plot polynomial fitted to (E, V) data
      string einppolyfitname = AGL_data.dirpathname + "/AGL_Plot_fitted_energyinput.dat";
      stringstream einppolyfitss;
      if(AGL_data.volumeinput.size() != xconfigurationvector.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors have different size" << einppolyfitname.c_str() <<  endl;
        aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size()" << AGL_data.volumeinput.size() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      } else {
        for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
          aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, plnv, 0);
          einppolyfitss << AGL_data.volumeinput.at(i) << "\t" << plnv << endl;
        }
      }
      if(!aurostd::stringstream2file(einppolyfitss, einppolyfitname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << einppolyfitname.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(einppolyfitss);	
    }

    // Find minimum of polynomial epol to find minimum of (E, V) curve
    // First bracket (i.e. find points on each side of) the minimum of (E, V) data
    itry = imin;
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Calling bracket_minimum to find imin" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    if(AGL_data.fittype == 0) {
      aglerror = AGL_functions::bracket_minimum (itry, AGL_data.energyinput, FileMESSAGE);
    } else {
      aglerror = AGL_functions::bracket_minimum_global (itry, AGL_data.energyinput, FileMESSAGE);
    }
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Cannot find minimum of (E, V) data" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    } else {
      aurostd::StringstreamClean(aus);
      if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
        aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is at point imin = " << itry << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // Find minimum of polynomial described by energypolynomialcoeffs, using itry as initial guess
    // Explicitly creating upper and lower bounds for minimum search as using "min" and "max" functions could be problematic with unsigned integers
    aglerror = AGL_functions::autocorrect_polynom_minimum (itry, xconfigurationvector, energypolynomialcoeffs, xmin, FileMESSAGE);
    if(aglerror == 0) {
      aurostd::StringstreamClean(aus);
      // OBSOLETE aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is at point xmin = " << xmin << " Bohr" << endl;
      if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
        aus << _AGLSTR_MESSAGE_ << "Minimum (relative) of (E, V) data is at point xmin = " << xmin << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Cannot find minimum of (E, V) data" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    }   
    Vmin = pow(xmin, 3);
    Vmin_Ang = Vmin * volumereference;
    Vmin_Bohr = Vmin_Ang * pow(angstrom2bohr, 3.0);
    aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_MESSAGE_ << "Volume at minimum of (E, V) data is V = " << Vmin << " Bohr^3" << endl;
    aus << _AGLSTR_MESSAGE_ << "Relative volume at minimum of (E, V) data is V = " << Vmin << endl;
    aus << _AGLSTR_MESSAGE_ << "Volume at minimum of (E, V) data is V = " << Vmin_Ang << " Angstrom^3" << endl;
    aus << _AGLSTR_MESSAGE_ << "Volume at minimum of (E, V) data is V = " << Vmin_Bohr << " Bohr^3" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // Evaluate polynomial at minimum point xmin to get Emin 
    AGL_functions::polynom_eval (xmin, energypolynomialcoeffs, Emin, 0);
    // Write out Emin in different units
    aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is Emin = " << Emin << " Hartree/cell" << endl;
    if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
      aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is Emin = " << Emin << " eV/cell" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is Emin = " << Emin * hartree2kjmol << " kJ/mol" << endl;
    if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
      aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is Emin = " << Emin * eV2kjmol << " kJ/mol" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is Emin = " << Emin * hart2ev << " eV/cell" << endl;
    if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
      aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is Emin = " << Emin / hart2ev << " Hartree/cell" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Record values of energy, volume and lattice strain factor at zero temperature and pressure for use elsewhere in AFLOW code
    AGL_data.volume_equilibrium_0p0T = Vmin_Ang;
    AGL_data.xlattice_equilibrium_0p0T = xmin;
    AGL_data.energy_equilibrium_0p0T = Emin;

    // Allocate vectors for storing static pressure properties if option not to truncate pressure range is selected
    if (AGL_data.run_all_pressure_temperature) {
      AGL_data.voleqmin.resize(AGL_data.pressure_external.size());
      AGL_data.VolumeStaticPressure.resize(AGL_data.pressure_external.size());
      AGL_data.StaticPressure.resize(AGL_data.pressure_external.size());
      gibbsenergyeV.resize(AGL_data.pressure_external.size());
      gibbsenergykjmol.resize(AGL_data.pressure_external.size());
      AGL_data.bulkmodulus.resize(AGL_data.pressure_external.size());
      relative_error.resize(AGL_data.pressure_external.size());
    }

    // Minimize G(static) as a function of V; G = Gibbs free energy
    gibbsenergypolynomialcoeffs.resize(energypolynomialcoeffs.size());
    for (uint i = 0; i < gibbsenergypolynomialcoeffs.size(); i++) {
      gibbsenergypolynomialcoeffs.at(i) = energypolynomialcoeffs.at(i);
    }
    gibbsenergydatapoints.resize(AGL_data.volumeinput.size());
    itry = imin;
    for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
      // Bracket the minimum with the numerical function
      if((gibbsenergydatapoints.size() != AGL_data.energyinput.size()) || (gibbsenergydatapoints.size() != AGL_data.volumeinput.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vector sizes are not equal" << endl;
        aus << _AGLSTR_ERROR_ + "gibbsenergydatapoints.size() = " << gibbsenergydatapoints.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      } else {     
        for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
          // OBSOLETE gibbsenergydatapoints.at(i) = AGL_data.energyinput.at(i) + AGL_data.volumeinput.at(i) * AGL_data.pressure_external.at(k) / au2GPa;
          gibbsenergydatapoints.at(i) = AGL_data.energyinput.at(i) + AGL_data.volumeinput.at(i) * (AGL_data.pressure_external.at(k) / eV_Ang3_to_GPa);
        }
      }
      if(AGL_data.fittype == 0) {
        aglerror = AGL_functions::bracket_minimum (itry, gibbsenergydatapoints, FileMESSAGE);
      } else {
        aglerror = AGL_functions::bracket_minimum_global (itry, gibbsenergydatapoints, FileMESSAGE);
      }
      // If bracket_minimum returns an error for zero pressure, AFLOW AGL exits giving an error
      // If bracket_minimum returns an error for p > zero, AFLOW AGL resets the maximum pressure to the previous value and skips the rest of the pressure loop
      // It then continues to complete the remainder of the GIBBS algorithm
      if(aglerror != 0) {
        // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of the iterations of this loop and continues with the next temperature value
        if (AGL_data.run_all_pressure_temperature) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + pV at p = " << AGL_data.pressure_external.at(k) << endl;
          aus << _AGLSTR_WARNING_ + "Skipping this pressure point" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {	
          if(k == 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Static: p = " << AGL_data.pressure_external.at(k) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", gibbsenergydatapoints.size() = " << gibbsenergydatapoints.size() << endl;
            aus << _AGLSTR_ERROR_ + "gibbsenergydatapoints = ";
            for(uint i = 0; i < gibbsenergydatapoints.size(); i++) {
              aus << gibbsenergydatapoints.at(i) << "\t";
            }
            aus << endl;
            aus << _AGLSTR_ERROR_ + "Cannot find minimum for pressure = " << AGL_data.pressure_external.at(k) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 2;
          } else {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Static: " << AGL_data.pressure_external.at(k) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", gibbsenergydatapoints.size() = " << gibbsenergydatapoints.size() << endl;
            aus << _AGLSTR_WARNING_ + "gibbsenergydatapoints = ";
            for(uint i = 0; i < gibbsenergydatapoints.size(); i++) {
              aus << gibbsenergydatapoints.at(i) << "\t";
            }
            aus << endl;
            aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + pV at p = " << AGL_data.pressure_external.at(k) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting maximum pressure value to p = " << AGL_data.pressure_external.at(k-1) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting number of pressure points to " << k << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            AGL_data.pressure_external.resize(k);
            break;
          }
        }
      }
      // Obtain the minimum of the fitted function
      // Adds pV term onto coefficient of x^3; x^3 ~ V
      // OBSOLETE gibbsenergypolynomialcoeffs.at(3) = energypolynomialcoeffs.at(3) + AGL_data.pressure_external.at(k) / au2GPa * volumereference;
      gibbsenergypolynomialcoeffs.at(3) = energypolynomialcoeffs.at(3) + (AGL_data.pressure_external.at(k) / eV_Ang3_to_GPa) * volumereference;
      aglerror = AGL_functions::autocorrect_polynom_minimum (itry, xconfigurationvector, gibbsenergypolynomialcoeffs, xmin, FileMESSAGE);
      // If autocorrect_polynom_minimum returns an error for zero pressure, AFLOW AGL exits giving an error
      // If autocorrect_polynom_minimum returns an error for p > zero, AFLOW AGL resets the maximum pressure to the previous value and skips the rest of the pressure loop
      // It then continues to complete the remainder of the GIBBS algorithm
      if(aglerror == 0) {
        aurostd::StringstreamClean(aus);
        if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
          aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is at point xmin = " << xmin << " Bohr" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else {
        if (AGL_data.run_all_pressure_temperature) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + pV at p = " << AGL_data.pressure_external.at(k) << endl;
          aus << _AGLSTR_WARNING_ + "Skipping this pressure point" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {
          if(k == 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Cannot find minimum of polynomial fit for E + pV at p = " << AGL_data.pressure_external.at(k) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 2;
          } else {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Cannot find minimum of polynomial fit for E + pV at p = " << AGL_data.pressure_external.at(k) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting maximum pressure value to p = " << AGL_data.pressure_external.at(k-1) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting number of pressure points to " << k << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            AGL_data.pressure_external.resize(k);
            break;
          }
        }
      }
      // Evaluates polynomial for (G, V) at minimum to get equilibrium values for V, G, and B
      // If option not to truncate pressure range is selected, writes values to appropriate vector element
      // Otherwise, populates vector dynamically up to highest pressure for which an enthalpy minimum can be obtained
      if (AGL_data.run_all_pressure_temperature) {
        AGL_data.voleqmin.at(k) = pow(xmin, 3.0) * volumereference;
        AGL_data.VolumeStaticPressure.at(k) = pow(xmin, 3.0) * volumereference;
        AGL_data.StaticPressure.at(k) = AGL_data.pressure_external.at(k);
        aglerror = AGL_functions::polynom_eval (xmin, gibbsenergypolynomialcoeffs, plnv, 0);
        // OBSOLETE gibbsenergy.push_back(plnv * hartree2kjmol);
        gibbsenergyeV.at(k) = plnv;
        gibbsenergykjmol.at(k) = plnv * eV2kjmol;
        // Save static Gibbs free energy (i.e. enthalpy: H = E + pV)
        AGL_pressure_enthalpy.pressure_external = AGL_data.pressure_external.at(k);
        AGL_pressure_enthalpy.enthalpy = plnv;
        AGL_data.AGL_pressure_enthalpy_list.push_back(AGL_pressure_enthalpy);
        aglerror = AGL_functions::polynom_eval (xmin, gibbsenergypolynomialcoeffs, plnv, 2);
        // OBSOLETE AGL_data.bulkmodulus.push_back(plnv * au2GPa * xmin * xmin / (9.0 * AGL_data.voleqmin.at(k)));
        AGL_data.bulkmodulus.at(k) = plnv * eV_Ang3_to_GPa * xmin * xmin / (9.0 * AGL_data.voleqmin.at(k));
        aglerror = AGL_functions::polynom_eval (xmin, energypolynomialcoeffs, energyvalue, 0);
        aglerror = AGL_functions::polynom_eval (xmin, energypolynomialerror, plnv, 0);
        energyerror = sqrt (aurostd::abs( plnv - energyvalue*energyvalue ));
        relative_error.at(k) = energyerror / max (aurostd::abs(energyvalue), aurostd::abs(energyvalue) + (energyerror/2.0));
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Error evaluating polynomials for pressure = " << AGL_data.pressure_external.at(k) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        }
      } else {
        AGL_data.voleqmin.push_back(pow(xmin, 3.0) * volumereference);
        AGL_data.VolumeStaticPressure.push_back(pow(xmin, 3.0) * volumereference);
        AGL_data.StaticPressure.push_back(AGL_data.pressure_external.at(k));
        aglerror = AGL_functions::polynom_eval (xmin, gibbsenergypolynomialcoeffs, plnv, 0);
        // OBSOLETE gibbsenergy.push_back(plnv * hartree2kjmol);
        gibbsenergyeV.push_back(plnv);
        gibbsenergykjmol.push_back(plnv * eV2kjmol);
        // Save static Gibbs free energy (i.e. enthalpy: H = E + pV)
        AGL_data.EnthalpyPressureeV.push_back(plnv);
        // Calculate bulk modulus from second derivative of energy with respect to volume
        aglerror = AGL_functions::polynom_eval (xmin, gibbsenergypolynomialcoeffs, plnv, 2);
        // OBSOLETE AGL_data.bulkmodulus.push_back(plnv * au2GPa * xmin * xmin / (9.0 * AGL_data.voleqmin.at(k)));
        AGL_data.bulkmodulus.push_back(plnv * eV_Ang3_to_GPa * xmin * xmin / (9.0 * AGL_data.voleqmin.at(k)));
        aglerror = AGL_functions::polynom_eval (xmin, energypolynomialcoeffs, energyvalue, 0);
        aglerror = AGL_functions::polynom_eval (xmin, energypolynomialerror, plnv, 0);
        energyerror = sqrt (aurostd::abs( plnv - energyvalue*energyvalue ));
        relative_error.push_back(energyerror / max (aurostd::abs(energyvalue), aurostd::abs(energyvalue)+energyerror/2.0));
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error evaluating polynomials" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 5;
        }
      }
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Finished initial pressure loop, starting EoS" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Write static Equation of State results to stringstream
    volume_0pressure = AGL_data.voleqmin.at(0);
    gibbsenergykjmol_0pressure = gibbsenergykjmol.at(0);
    binp_bcnt = AGL_data.bulkmodulus.at(0);
    if(AGL_data.i_eqn_of_state >= 0) {
      AGL_data.outfiless << "Static EOS calculation - Numerical results" << endl;
      // OBSOLETE AGL_data.outfiless << "Vmin(static;  P=0)    = " << volume_0pressure << " bohr^3" << endl;
      AGL_data.outfiless << "Vmin(static;  P=0)    = " << volume_0pressure * pow(angstrom2bohr, 3.0) << " bohr^3" << endl;
      AGL_data.outfiless << "Gmin(static;  P=0)    = " << gibbsenergykjmol_0pressure << " kJ/mol" << endl;
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "NUMERICAL EQUILIBRIUM PROPERTIES" << endl;
      AGL_data.outfiless << "================================" << endl; 
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "  P(GPa)         G(kJ/mol)       V(bohr^3)         V/V0          B(GPa)          rel.err." << endl;
      AGL_data.outfiless << " ----------------------------------------------------------------------------------------- " << endl;
      for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
        // OBSOLETE AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "        " << setw(10) << gibbsenergy.at(k) << "      " << setw(10) << AGL_data.voleqmin.at(k) << "      " << setw(6) << setprecision(5) << AGL_data.voleqmin.at(k)/volume_0pressure << "    " << setw(12) << setprecision(2) << AGL_data.bulkmodulus.at(k) << "   " << setw(15) << setprecision(6) << relative_error.at(k) << endl;
        AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "        " << setw(10) << gibbsenergykjmol.at(k) << "      " << setw(10) << AGL_data.voleqmin.at(k) * pow(angstrom2bohr, 3.0) << "      " << setw(6) << setprecision(5) << AGL_data.voleqmin.at(k)/volume_0pressure << "    " << setw(12) << setprecision(2) << AGL_data.bulkmodulus.at(k) << "   " << setw(15) << setprecision(6) << relative_error.at(k) << endl;

      }
    }
    else if(AGL_data.i_debye == -1) {
      for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
        // OBSOLETE AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "\t" << setw(6) << 0.0 << "\t" << setw(10) << AGL_data.voleqmin.at(k) << "\t" << setw(11) << setprecision(5) << gibbsenergy.at(k) << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k)  << endl;
        AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "\t" << setw(6) << 0.0 << "\t" << setw(10) << AGL_data.voleqmin.at(k) * pow(angstrom2bohr, 3.0) << "\t" << setw(11) << setprecision(5) << gibbsenergykjmol.at(k) << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k)  << endl;
      }
    } else {
      for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
        // OBSOLETE AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "\t" << setw(6) << "static" << "\t" << setw(10) << AGL_data.voleqmin.at(k) << "\t" << setw(11) << setprecision(5) << gibbsenergy.at(k) << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k)  << endl;
        AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "\t" << setw(6) << "static" << "\t" << setw(10) << AGL_data.voleqmin.at(k) * pow(angstrom2bohr, 3.0) << "\t" << setw(11) << setprecision(5) << gibbsenergykjmol.at(k) << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k)  << endl;
      }
    }
    // Allocating arrays to save EOS variables
    AGL_data.d2EnergydVolume2_static.resize(AGL_data.volumeinput.size());
    AGL_data.IntEnergStatic.resize(AGL_data.volumeinput.size());
    F_Helmholtzvinp.resize(AGL_data.volumeinput.size());
    AGL_data.d2EnergydVolume2_dynamic.resize(AGL_data.pressure_external.size());
    AGL_data.pressure_static.resize(AGL_data.pressure_external.size());
    AGL_data.gamma_G.resize(AGL_data.pressure_external.size());
    AGL_data.pfit.resize(AGL_data.pressure_external.size());
    AGL_data.alpha.resize(AGL_data.pressure_external.size());
    theta.resize(AGL_data.pressure_external.size());
    Cv.resize(AGL_data.pressure_external.size());
    AGL_data.astatic.resize(AGL_data.birchfitorder_iG+1);
    press_stat.resize(AGL_data.volumeinput.size());

    if (AGL_data.run_all_pressure_temperature) {
      if((AGL_data.d2EnergydVolume2_static.size() != xconfigurationvector.size()) || (AGL_data.energyinput.size() != xconfigurationvector.size()) || (AGL_data.volumeinput.size() != xconfigurationvector.size()) ) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
        aus << _AGLSTR_ERROR_ + "xconfigurationvector.size() = " << xconfigurationvector.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      }
      for (uint i = 0; i < xconfigurationvector.size(); i++) {
        double Ex, dEdx, d2Edx2, volumex3;
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, Ex, 0);
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, dEdx, 1);
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, d2Edx2, 2);
        double dEdx2_xd2Edx2 = xconfigurationvector.at(i) * d2Edx2 - 2.0 * dEdx;
        volumex3 = 3.0 * AGL_data.volumeinput.at(i);
        AGL_data.d2EnergydVolume2_static.at(i) = dEdx2_xd2Edx2 * xconfigurationvector.at(i) / (volumex3 * volumex3);
      }
    } else {     
      // Debye temperature: Numerical evaluation of thermal properties
      if(AGL_data.i_eqn_of_state == 0) {
        aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, energypolynomialcoeffs, xconfigurationvector, true, AGL_data, FileMESSAGE);
        AGL_data.outfiless << endl;
        AGL_data.outfiless <<  "Debye temperature - numerical derivatives" << endl;
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
      }
      else if(AGL_data.i_eqn_of_state < 0) {
        aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, energypolynomialcoeffs, xconfigurationvector, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
      }
      // Vinet equation of state evaluation of thermal properties
      else if(AGL_data.i_eqn_of_state == 1) {
        // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
        // OBSOLETE ge0peV = gibbsenergykjmol_0pressure/eV2kjmol;
        // OBSOLETE aglerror = AGL_functions::vinet_eos (volume_0pressure, g0dhy2kjm, true, AGL_data, FileMESSAGE);
        gibbsenergyeV_0pressure = gibbsenergyeV.at(0);
        aglerror = AGL_functions::vinet_eos (volume_0pressure, gibbsenergyeV_0pressure, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "Debye temperature - Vinet EOS derivatives" << endl;
      }
      // Second order (fourth order in energy) Birch-Murnaghan equation of state evaluation of thermal properties
      else if(AGL_data.i_eqn_of_state == 2) {
        // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
        // OBSOLETE ge0peV = gibbsenergykjmol_0pressure/eV2kjmol;
        // OBSOLETE aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, g0dhy2kjm, true, AGL_data, FileMESSAGE);
        gibbsenergyeV_0pressure = gibbsenergyeV.at(0);
        aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, gibbsenergyeV_0pressure, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error in Birch-Murnaghan fit" << endl;
          aus << _AGLSTR_ERROR_ + "Check input pressure values (first pressure point must be p = 0)" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "Debye temperature - BIRCH-MURNAGHAN EOS derivatives" << endl;
      }
      // Vinet equation fitting, but numerical calculation to obtain all properties.
      else if(AGL_data.i_eqn_of_state == 3) {
        // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
        // OBSOLETE ge0peV = gibbsenergykjmol_0pressure/eV2kjmol;
        // OBSOLETE aglerror = AGL_functions::vinet_eos (volume_0pressure, g0dhy2kjm, true, AGL_data, FileMESSAGE);
        gibbsenergyeV_0pressure = gibbsenergyeV.at(0);
        aglerror = AGL_functions::vinet_eos (volume_0pressure, gibbsenergyeV_0pressure, true, AGL_data, FileMESSAGE);
        aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, energypolynomialcoeffs, xconfigurationvector, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "Debye temperature - numerical derivatives" << endl;
      }
      // Second order (fourth order in energy) Birch equation fitting, but numerical calculation to obtain all properties.
      else if(AGL_data.i_eqn_of_state == 4) {
        // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
        // OBSOLETE ge0peV = gibbsenergykjmol_0pressure/eV2kjmol;      
        // OBSOLETE aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, g0dhy2kjm, true, AGL_data, FileMESSAGE);
        gibbsenergyeV_0pressure = gibbsenergyeV.at(0);
        aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, gibbsenergyeV_0pressure, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error in Birch-Murnaghan fit" << endl;
          aus << _AGLSTR_ERROR_ + "Check input pressure values (first pressure point must be p = 0)" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, energypolynomialcoeffs, xconfigurationvector, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "Debye temperature - numerical derivatives" << endl;
      }
      // Equation of state of Baonza-Caceres-Nuñez-Taravillo to obtain thermal properties
      else if(AGL_data.i_eqn_of_state == 5) {
        // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
        // OBSOLETE ge0peV = gibbsenergykjmol_0pressure/eV2kjmol;    
        // OBSOLETE aglerror = AGL_functions::bcnt_eos (volume_0pressure, g0dhy2kjm, binp_bcnt, true, AGL_data, FileMESSAGE);
        gibbsenergyeV_0pressure = gibbsenergyeV.at(0);
        aglerror = AGL_functions::bcnt_eos (volume_0pressure, gibbsenergyeV_0pressure, binp_bcnt, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);;
          return aglerror;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "Debye temperature - BCNT EOS derivatives" << endl;
      }
      // Equation of state of Baonza-Caceres-Nuñez-Taravillo, but numerical calculation to obtain all properties.
      else if(AGL_data.i_eqn_of_state == 6) {
        // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
        // OBSOLETE ge0peV = gibbsenergykjmol_0pressure/eV2kjmol;   
        // OBSOLETE aglerror = AGL_functions::bcnt_eos (volume_0pressure, g0dhy2kjm, binp_bcnt, true, AGL_data, FileMESSAGE);
        gibbsenergyeV_0pressure = gibbsenergyeV.at(0);
        aglerror = AGL_functions::bcnt_eos (volume_0pressure, gibbsenergyeV_0pressure, binp_bcnt, true, AGL_data, FileMESSAGE);
        aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, energypolynomialcoeffs, xconfigurationvector, true, AGL_data, FileMESSAGE);
        if(aglerror != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Error equation of state fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return aglerror;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "Debye temperature - numerical derivatives" << endl;
      }
      // Write Poisson coefficient and poisson ratio function
      if(AGL_data.i_eqn_of_state >= 0 && (AGL_data.i_debye == 0 || AGL_data.i_debye >= 2)) {
        AGL_data.outfiless << "Poisson coefficient: " << AGL_data.poissonratio << ", Poisson ratio function: " << AGL_data.poissonratiofunction << endl;
        AGL_data.outfiless << endl;
      }
      // Calculate Debye temperatures at each volume
      if(AGL_data.dbulkmodulusdpV_0pressure < 0) {
        AGL_data.outfiless << "gibbs: Warning! B'<0, will use i_debye=3" << endl;
        AGL_data.i_debye = 3;
      }
    }
    if(AGL_data.i_debye == 1) {
      if(AGL_data.i_eqn_of_state >= 0) {
        AGL_data.outfiless << "  V(bohr^3)      TDebye(K)      Computed(K)" << endl;
        AGL_data.outfiless << "-----------    -----------    -------------" << endl;
      }
      aglerror = AGL_functions::debye_polynom_fit (tdebyepolynomialcoeffs, AGL_data, FileMESSAGE);
      if(aglerror != 0) {
        if(aglerror == 2) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Problem inverting matrix to fit polynomial" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 4;
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 2;
        }
      }
    }
    else if(AGL_data.i_debye == 3) {
      if(AGL_data.volumeinput.size() != AGL_data.d2EnergydVolume2_static.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      } else {
        // OBSOLETE tdebyemin = pow(6*pi*pi*AGL_data.natoms*AGL_data.volumeinput.at(imin)*AGL_data.volumeinput.at(imin), third) / physconstkbau * AGL_data.poissonratiofunction * sqrt(aurostd::abs(AGL_data.d2EnergydVolume2_static.at(imin))/AGL_data.cellmass);
        // OBSOLETE tdebyemin = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * (sqrt(AGL_data.volumeinput.at(imin))), third) * AGL_data.poissonratiofunction * sqrt((aurostd::abs(AGL_data.d2EnergydVolume2_static.at(imin) * AGL_data.volumeinput.at(imin))/AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
        tdebyemin = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * AGL_data.volumeinput.at(imin) * AGL_data.volumeinput.at(imin), third) * AGL_data.poissonratiofunction * sqrt((aurostd::abs(AGL_data.d2EnergydVolume2_static.at(imin)) / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
      }
    }
    else {
      if(AGL_data.i_eqn_of_state >= 0) {
        AGL_data.outfiless << "   V(bohr^3)      TDebye(K)" << endl;
        AGL_data.outfiless << " -----------     ----------" << endl;
      }
    }

    // Test: print out polynomials for dE / dx, etc.
    if(LDEBUG) {
      double tf0, tf1, tf2;
      string ofiletf0name = AGL_data.dirpathname + "/AGL_Poly_Ex0.dat";
      stringstream ofiletf0ss;
      string ofiletf1name = AGL_data.dirpathname + "/AGL_Poly_Ex1.dat";
      stringstream ofiletf1ss;
      string ofiletf2name = AGL_data.dirpathname + "/AGL_Poly_Ex2.dat";
      stringstream ofiletf2ss;
      for (uint i = 0; i < xconfigurationvector.size(); i++) {
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, tf0, 0);
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, tf1, 1);
        aglerror = AGL_functions::polynom_eval (xconfigurationvector.at(i), energypolynomialcoeffs, tf2, 2);
        ofiletf0ss << xconfigurationvector.at(i) << "\t" << tf0 << endl;
        ofiletf1ss << xconfigurationvector.at(i) << "\t" << tf1 << endl;
        ofiletf2ss << xconfigurationvector.at(i) << "\t" << tf2 << endl;
      }
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Failed to evaluate polynomial" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 5;
      }
      if(!aurostd::stringstream2file(ofiletf0ss, ofiletf0name, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofiletf0name.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(ofiletf0ss);
      if(!aurostd::stringstream2file(ofiletf0ss, ofiletf0name, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofiletf0name.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(ofiletf1ss);
      if(!aurostd::stringstream2file(ofiletf0ss, ofiletf0name, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofiletf0name.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(ofiletf2ss);
    }

    uint ij = 0;
    // Checks second derivative is positive at all points
    // If not, the function is convex at that point and that point is skipped
    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      if(AGL_data.i_debye != 1 && (AGL_data.d2EnergydVolume2_static.at(i) <= 0.0 && (AGL_data.fittype < 2 || AGL_data.fittype == 3))) {
        if(i > itry) {
          // OBSOLETE AGL_data.volumeinput.resize(i-1);
          // OBSOLETE AGL_data.energyinput.resize(i-1);
          // OBSOLETE xconfigurationvector.resize(i-1);
          // OBSOLETE gibbsenergydatapoints.resize(i-1);
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Warning! convex function at i = " << i << endl;
          aus << _AGLSTR_WARNING_ + "Warning! d2EnergydVolume2_static = " << AGL_data.d2EnergydVolume2_static.at(i) << endl;
          aus << _AGLSTR_WARNING_ + "The points following this one will be discarded" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          AGL_data.volumeinput.resize(i);
          AGL_data.energyinput.resize(AGL_data.volumeinput.size());
          xconfigurationvector.resize(AGL_data.volumeinput.size());
          gibbsenergydatapoints.resize(AGL_data.volumeinput.size());
          AGL_data.IntEnergStatic.resize(AGL_data.volumeinput.size());
          F_Helmholtzvinp.resize(AGL_data.volumeinput.size());
          AGL_data.d2EnergydVolume2_static.resize(AGL_data.volumeinput.size());
        } else {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Warning! convex function at i = " << i << endl;
          aus << _AGLSTR_WARNING_ + "Warning! d2EnergydVolume2_static = " << AGL_data.d2EnergydVolume2_static.at(i) << endl;
          aus << _AGLSTR_WARNING_ + "This point will be skipped and the next concave one taken as the next point" << endl;
          aus << _AGLSTR_WARNING_ + "Recommend increasing the number of k-points and rerunning AGL" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else {
        if(i != ij) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "i = " << i << endl;
          aus << _AGLSTR_WARNING_ + "ij = " << ij << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          AGL_data.volumeinput.at(ij) = AGL_data.volumeinput.at(i);
          AGL_data.energyinput.at(ij) = AGL_data.energyinput.at(i);
          xconfigurationvector.at(ij) = xconfigurationvector.at(i);
          gibbsenergydatapoints.at(ij) = gibbsenergydatapoints.at(i);
          AGL_data.IntEnergStatic.at(ij) = AGL_data.IntEnergStatic.at(i);
          F_Helmholtzvinp.at(ij) = F_Helmholtzvinp.at(ij);
          AGL_data.d2EnergydVolume2_static.at(ij) = AGL_data.d2EnergydVolume2_static.at(i);
          if(i == imin) {
            imin = ij;
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "imin = " << imin << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
        // tdebyevol is the Debye temperature for the structure with volume AGL_data.volumeinput.at(ij)
        if(AGL_data.volumeinput.size() != AGL_data.d2EnergydVolume2_static.size()) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
          aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
          aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 3;
        } else {
          // OBSOLETE tdebyevol = pow(6*pi*pi*AGL_data.natoms*AGL_data.volumeinput.at(ij)*AGL_data.volumeinput.at(ij), third) / physconstkbau * AGL_data.poissonratiofunction * sqrt(aurostd::abs(AGL_data.d2EnergydVolume2_static.at(ij))/AGL_data.cellmass);
          // OBSOLETE tdebyevol = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * (sqrt(AGL_data.volumeinput.at(ij))), third) * AGL_data.poissonratiofunction * sqrt((aurostd::abs(AGL_data.d2EnergydVolume2_static.at(ij) * AGL_data.volumeinput.at(imin)) / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
          tdebyevol = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * AGL_data.volumeinput.at(ij) * AGL_data.volumeinput.at(ij), third) * AGL_data.poissonratiofunction * sqrt((aurostd::abs(AGL_data.d2EnergydVolume2_static.at(ij)) / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
          if(LDEBUG) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_MESSAGE_ + "tdebyevol = " << tdebyevol << endl;
            aus << _AGLSTR_MESSAGE_ + "AGL_data.d2EnergydVolume2_static.at(ij) = " << AGL_data.d2EnergydVolume2_static.at(ij) << endl;
            aus << _AGLSTR_MESSAGE_ + "AGL_data.volumeinput.at(ij) = " << AGL_data.volumeinput.at(ij) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        }
        if(AGL_data.i_debye == 3) {
          AGL_data.tdebye.push_back(tdebyemin);
          if(AGL_data.i_eqn_of_state >= 0) {
            // OBSOLETE AGL_data.outfiless << setw(12) << setprecision(6) << fixed << AGL_data.volumeinput.at(ij) << "    " << setw(12) << AGL_data.tdebye.at(ij) << "    " << setw(12) << tdebyevol << endl;
            AGL_data.outfiless << setw(12) << setprecision(6) << fixed << (AGL_data.volumeinput.at(ij) * pow(angstrom2bohr, 3.0)) << "    " << setw(12) << AGL_data.tdebye.at(ij) << "    " << setw(12) << tdebyevol << endl;
          }
        } else if(AGL_data.i_debye == 1) {
          if(AGL_data.i_eqn_of_state >= 0) { 
            // OBSOETE AGL_data.outfiless << setw(12) << setprecision(6) << fixed << AGL_data.volumeinput.at(ij) << "    " << setw(12) << AGL_data.tdebye.at(ij) << endl;
            AGL_data.outfiless << setw(12) << setprecision(6) << fixed << (AGL_data.volumeinput.at(ij) * pow(angstrom2bohr, 3.0)) << "    " << setw(12) << AGL_data.tdebye.at(ij) << endl;
          }
        } else {
          AGL_data.tdebye.push_back(tdebyevol);
        }
        if(AGL_data.i_eqn_of_state >= 0) { 
          // OBSOLETE AGL_data.outfiless << setw(12) << setprecision(6) << fixed << AGL_data.volumeinput.at(ij) << "    " << setw(12) << AGL_data.tdebye.at(ij) << endl;
          AGL_data.outfiless << setw(12) << setprecision(6) << fixed << (AGL_data.volumeinput.at(ij) * pow(angstrom2bohr, 3.0)) << "    " << setw(12) << AGL_data.tdebye.at(ij) << endl;
        }
        ij++;
      }
    }
    if(LDEBUG) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "ij = " << ij << endl;
      aus << _AGLSTR_WARNING_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
      aus << _AGLSTR_WARNING_ + "AGL_data.volumeinput = ";
      for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
        aus << AGL_data.volumeinput.at(i) << "\t";	
      }
      aus << endl;
      aus << _AGLSTR_WARNING_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
      aus << _AGLSTR_WARNING_ + "AGL_data.energyinput = ";
      for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
        aus << AGL_data.energyinput.at(i) << "\t";	
      }
      aus << endl;
      aus << _AGLSTR_WARNING_ + "AGL_data.tdebye.size() = " << AGL_data.tdebye.size() << endl;
      aus << _AGLSTR_WARNING_ + "AGL_data.tdebye = ";
      for (uint i = 0; i < AGL_data.tdebye.size(); i++) {
        aus << AGL_data.tdebye.at(i) << "\t";	
      }
      aus << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(AGL_data.volumeinput.size() != ij) {
      AGL_data.volumeinput.resize(ij);
      AGL_data.energyinput.resize(ij);
      xconfigurationvector.resize(ij);
      gibbsenergydatapoints.resize(ij);
      AGL_data.IntEnergStatic.resize(ij);
      F_Helmholtzvinp.resize(ij);
      AGL_data.d2EnergydVolume2_static.resize(ij);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Number of (E, V) points set to = " << ij << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(LDEBUG) {
      string ofiled2EnergydVolume2_staticname = AGL_data.dirpathname + "/AGL_Poly_d2EnergydVolume2_static.dat";
      stringstream ofiled2EnergydVolume2_staticss;
      if(AGL_data.d2EnergydVolume2_static.size() != xconfigurationvector.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors have different size" << endl;
        aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size()" << AGL_data.d2EnergydVolume2_static.size() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      }
      for (uint i = 0; i < xconfigurationvector.size(); i++) {
        ofiled2EnergydVolume2_staticss << xconfigurationvector.at(i) << "\t" << AGL_data.d2EnergydVolume2_static.at(i) << endl;
      }	
      if(!aurostd::stringstream2file(ofiled2EnergydVolume2_staticss, ofiled2EnergydVolume2_staticname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofiled2EnergydVolume2_staticname.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(ofiled2EnergydVolume2_staticss);
      string ofiletdebiname = AGL_data.dirpathname + "/AGL_Plot_tdebye_initial.dat";
      stringstream ofiletdebiss;
      for (uint i = 0; i < xconfigurationvector.size(); i++) {
        ofiletdebiss << xconfigurationvector.at(i) << "\t" << AGL_data.tdebye.at(i) << endl;
      }	
      if(!aurostd::stringstream2file(ofiletdebiss, ofiletdebiname, "WRITE")) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Unable to open file " << ofiletdebiname.c_str() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(ofiletdebiss);
    }
    if(AGL_data.volumeinput.size() < AGL_data.nstructsinit) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Some points showed convex patterns and have been removed from the data set" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      AGL_data.EV_noise = true;
      // int nskippedpoints = AGL_data.nstructsinit - AGL_data.volumeinput.size();
      nskippedpoints = AGL_data.nstructsinit - AGL_data.volumeinput.size();
      if(nskippedpoints > AGL_data.skiparunsmax) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Total number of skipped points exceeds maximum limit" << endl;
        aus << _AGLSTR_ERROR_ + "Total number of skipped points = " << nskippedpoints << endl;
        aus << _AGLSTR_ERROR_ + "Maximum limit of skipped points = " << AGL_data.skiparunsmax << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aglerror = 2;
        return aglerror;
      }
    }
    // End of static calculation: exit if IDEBYE = -1
    if(AGL_data.i_debye == -1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "End of static run ok!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    }
    // Allocate memory for calculating E + F (Helmholtz free energy)
    helmholtzenergydatapoints.resize(AGL_data.volumeinput.size());
    // Initialize self-consistent Debye variables
    // Debye temperature is calculated self-consistently if IDEBYE = 2
    if(AGL_data.i_debye == 2) {
      if(helmholtzenergydatapoints.size() != xconfigurationvector.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors have different size" << endl;
        aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
        aus << _AGLSTR_ERROR_ + "helmholtzenergydatapoints.size()" << helmholtzenergydatapoints.size() <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      }
      aglerror = AGL_functions::self_consistent_debye (dzero, energypolynomialcoeffs, energypolynomialerror, xconfigurationvector, helmholtzenergydatapoints, imin, true, press_stat, AGL_data, FileMESSAGE);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Problem in self-consistent Debye!" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 4;
      }
    }
    // If i_optimize_beta = 2 the variable optimize_beta is changed to false
    // Optimization of bcnt_beta is only performed for the static calculation
    // Optimized value of bcnt_beta from the static calculation used for the finite temperature calculations
    if(AGL_data.i_eqn_of_state == 5 || AGL_data.i_eqn_of_state == 6) {
      if(AGL_data.i_optimize_beta == 2 ) AGL_data.optimize_beta = false;
    }

    // Stringstream to save thermal properties in a plottable format
    // Similar format to "THERMO" file written by APL
    aurostd::StringstreamClean(AGL_data.outfiletss);
    if(AGL_data.i_eqn_of_state >= 0) {
      AGL_data.outfiletss << "#   T(K)    U(meV/cell)     F(meV/cell)      S(kB/cell)     Cv(kB/cell)      Theta_D(K)     Gruneisen parameter" << std::endl;
    }

    // Allocate vectors to save thermal properties at all (p, T) if requested by user
    if(AGL_data.savedatapressure) {
      AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
      AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
      AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
      AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
      AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
      AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
      AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
      AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
      AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
      AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
      AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
    }
    AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
    AGL_data.xminsav.resize(AGL_data.temperature_external.size());
    AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
    AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
    // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());

    // Loop over temperatures
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Starting temperature loop" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
      aurostd::StringstreamClean(aus);
      if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
        aus << _AGLSTR_MESSAGE_ << "Temperature = " << AGL_data.temperature_external.at(j) << "K" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      AGL_data.max_temperature = AGL_data.temperature_external.at(j);
      // Self consistent Debye temperatures
      if(AGL_data.i_debye == 2) {
        if(helmholtzenergydatapoints.size() != xconfigurationvector.size()) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Vectors have different size" << endl;
          aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
          aus << _AGLSTR_ERROR_ + "helmholtzenergydatapoints.size()" << helmholtzenergydatapoints.size() <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 3;
        }
        aglerror = AGL_functions::self_consistent_debye(AGL_data.temperature_external.at(j), helmholtzenergypolynomialcoeffs, helmholtzpolynomialerror, xconfigurationvector, helmholtzenergydatapoints, imin, false, press_stat, AGL_data, FileMESSAGE);
        // If self_consistent_debye returns an error for zero temperature, AFLOW AGL exits giving an error
        // If self_consistent_debye returns an error for T > zero, AFLOW AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop
        // It then continues to complete the remainder of the GIBBS algorithm
        if(aglerror != 0) {
          // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of the iterations of this loop and continues with the next temperature value
          if (AGL_data.run_all_pressure_temperature) {
            continue;
          } else {
            // Otherwise, truncates temperature range and skips remainder of temperature loop
            if(j == 0) {      
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Self-consistent polynomial fit: T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
              aus << _AGLSTR_ERROR_ + "helmholtzenergydatapoints = ";
              for (uint i = 0; i < helmholtzenergydatapoints.size(); i++) {
                aus << helmholtzenergydatapoints.at(i) << "\t";
              }
              aus << endl;
              aus << _AGLSTR_ERROR_ + "Cannot find minimum for temperature = " << AGL_data.temperature_external.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 2;
            } else {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ + "Self-consistent polynomial fit: T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
              aus << _AGLSTR_WARNING_ + "helmholtzenergydatapoints = ";
              for (uint i = 0; i < helmholtzenergydatapoints.size(); i++) {
                aus << helmholtzenergydatapoints.at(i) << "\t";
              }
              aus << endl;
              aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + F at T = " << AGL_data.temperature_external.at(j) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting maximum value of T to T = " << AGL_data.temperature_external.at(j) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting number of temperature points to " << j << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              AGL_data.temperature_external.resize(j);
              if(AGL_data.savedatapressure) {
                AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
                AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
                AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
                AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
                AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
                AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
                AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
              }
              AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
              AGL_data.xminsav.resize(AGL_data.temperature_external.size());
              AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
              AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
              // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());
              break;
            }
          }
        }
        aglerror = AGL_functions::debye_polynom_fit (tdebyepolynomialcoeffs, AGL_data, FileMESSAGE);
        // If debye_polynom_fit returns an error for zero temperature, AFLOW AGL exits giving an error
        // If debye_polynom_fit returns an error for T > zero, AFLOW AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop 
        // It then continues to complete the remainder of the GIBBS algorithm
        if(aglerror != 0) {
          // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of the iterations of this loop and continues with the next temperature value
          if (AGL_data.run_all_pressure_temperature) {
            continue;
          } else {
            if(j == 0) {      
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
              aus << _AGLSTR_ERROR_ + "Debye polynomial fit: T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
              aus << _AGLSTR_ERROR_ + "Cannot find minimum for temperature = " << AGL_data.temperature_external.at(j) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              if(aglerror == 2) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_ERROR_ + "Problem inverting matrix to fit polynomial" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                return 4;
              } else {
                return 2;
              }
            } else {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
              aus << _AGLSTR_WARNING_ + "Debye polynomial fit: T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
              aus << _AGLSTR_WARNING_ + "Resetting maximum value of T to T = " << AGL_data.temperature_external.at(j) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting number of temperature points to " << j << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              AGL_data.temperature_external.resize(j);
              if(AGL_data.savedatapressure) {
                AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
                AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
                AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
                AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
                AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
                AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
                AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
              }
              AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
              AGL_data.xminsav.resize(AGL_data.temperature_external.size());
              AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
              AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
              // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());
              break;
            }
          }    
        }
      }
      // End of self-consistent temperature loop

      // Obtain vibrational Helmholtz function fit (Debye model) for this value of the temperature, T
      // At each volume, use the Debye temperature to get the vibrational Helmholtz energy (F = U - TS) 
      // Add the value of F at each volume to the DFT final energy at that volume to AGL_data.F_Helmholtzvinp
      // Fit (F + E, V) data to get the polynomial helmholtzenergypolynomialcoeffs (zero pressure, constant volume)
      if(LDEBUG && AGL_data.fittype == 5) {
        string ofiletdebname = AGL_data.dirpathname + "/AGL_Plot_tdebye.dat";
        stringstream ofiletdebss;
        if(AGL_data.tdebye.size() != xconfigurationvector.size()) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Vectors have different size" << endl;
          aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
          aus << _AGLSTR_ERROR_ + "AGL_data.tdebye.size()" << AGL_data.tdebye.size() <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 3;
        }
        for (uint i = 0; i < xconfigurationvector.size(); i++) {
          ofiletdebss << xconfigurationvector.at(i) << "\t" << AGL_data.tdebye.at(i) << endl;
        }
        if(!aurostd::stringstream2file(ofiletdebss, ofiletdebname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file " << ofiletdebname.c_str() <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        aurostd::StringstreamClean(ofiletdebss);
      }
      if((helmholtzenergydatapoints.size() != AGL_data.energyinput.size()) || (helmholtzenergydatapoints.size() != F_Helmholtzvinp.size()) || (helmholtzenergydatapoints.size() != AGL_data.tdebye.size()) ) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
        aus << _AGLSTR_ERROR_ + "F_Helmholtzvinp.size() = " << F_Helmholtzvinp.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.tdebye.size() = " << AGL_data.tdebye.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      } else {
        for(uint i = 0; i < helmholtzenergydatapoints.size(); i++) {
          aglerror = AGL_functions::thermal_properties (AGL_data.tdebye.at(i), AGL_data.temperature_external.at(j), DebyeIntegral, DebyeIntegralerror, UIntEnergy, Cvt, F_Helmholtzvinp.at(i), S_Entropy, AGL_data, FileMESSAGE);
          helmholtzenergydatapoints.at(i) = AGL_data.energyinput.at(i) + F_Helmholtzvinp.at(i);
        }
      }
      if(LDEBUG && AGL_data.fittype == 5) {
        string ofilehelmname = AGL_data.dirpathname + "/AGL_Plot_Helmholtz_energy.dat";
        stringstream ofilehelmss;
        if(F_Helmholtzvinp.size() != xconfigurationvector.size()) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Vectors have different size" << endl;
          aus << _AGLSTR_ERROR_ + "xconfigurationvector.size()" << xconfigurationvector.size() <<  endl;
          aus << _AGLSTR_ERROR_ + "F_Helmholtzvinp.size()" << F_Helmholtzvinp.size() <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 3;
        }
        for (uint i = 0; i < xconfigurationvector.size(); i++) {
          ofilehelmss << xconfigurationvector.at(i) << "\t" << F_Helmholtzvinp.at(i) << endl;
        }
        if(!aurostd::stringstream2file(ofilehelmss, ofilehelmname, "WRITE")) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Unable to open file " << ofilehelmname.c_str() <<  endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        aurostd::StringstreamClean(ofilehelmss);
      } else if(LDEBUG) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "Polynomial fit: T = " << AGL_data.temperature_external.at(j) << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
        aus << _AGLSTR_WARNING_ + "helmholtzenergydatapoints = ";
        for (uint i = 0; i < helmholtzenergydatapoints.size(); i++) {
          aus << helmholtzenergydatapoints.at(i) << "\t";
        }
        aus << endl;
        aus << _AGLSTR_WARNING_ + "F_Helmholtzvinp.size() = " << F_Helmholtzvinp.size() << endl;
        aus << _AGLSTR_WARNING_ + "F_Helmholtzvinp = ";
        for (uint i = 0; i < F_Helmholtzvinp.size(); i++) {
          aus << F_Helmholtzvinp.at(i) << "\t";
        }
        aus << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
      }
      // Brackets minimum of E + F
      if(AGL_data.fittype == 0) {
        aglerror = AGL_functions::bracket_minimum (imin, helmholtzenergydatapoints, FileMESSAGE);
      } else {
        aglerror = AGL_functions::bracket_minimum_global (imin, helmholtzenergydatapoints, FileMESSAGE);
      }
      // If bracket_minimum returns an error for zero temperature, AFLOW AGL exits giving an error
      // If bracket_minimum returns an error for T > zero, AFLOW AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop
      // It then continues to complete the remainder of the GIBBS algorithm
      if(aglerror != 0) {
        // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of the iterations of this loop and continues with the next temperature value
        if (AGL_data.run_all_pressure_temperature) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + F at T = " << AGL_data.temperature_external.at(j) << endl;
          aus << _AGLSTR_WARNING_ + "Skipping this temperature point" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // OBSOLETE continue;
          minpolyerror = aglerror;
          aglerror = 0;
        } else {
          if(j == 0) {      
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Bracket minimum: T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
            aus << _AGLSTR_ERROR_ + "helmholtzenergydatapoints = ";
            for (uint i = 0; i < helmholtzenergydatapoints.size(); i++) {
              aus << helmholtzenergydatapoints.at(i) << "\t";
            }
            aus << endl;
            aus << _AGLSTR_ERROR_ + "AGL_data.energyinput = ";
            for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
              aus << AGL_data.energyinput.at(i) << "\t";
            }
            aus << endl;
            aus << _AGLSTR_ERROR_ + "F_helmholtzvinp = ";
            for (uint i = 0; i < F_Helmholtzvinp.size(); i++) {
              aus << F_Helmholtzvinp.at(i) << "\t";
            }
            aus << endl;
            aus << _AGLSTR_ERROR_ + "Cannot find minimum for temperature = " << AGL_data.temperature_external.at(j) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 2;
          } else {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Bracket minimum: T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
            aus << _AGLSTR_WARNING_ + "AGL_data.functiontofit = ";
            for (uint i = 0; i < helmholtzenergydatapoints.size(); i++) {
              aus << helmholtzenergydatapoints.at(i) << "\t";
            }
            aus << endl;
            aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + F at T = " << AGL_data.temperature_external.at(j) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting maximum value of T to T = " << AGL_data.temperature_external.at(j) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting number of temperature points to " << j << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            AGL_data.temperature_external.resize(j);
            if(AGL_data.savedatapressure) {
              AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
              AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
              AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
              AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
              AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
              AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
              AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
              AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
              AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
              AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
              AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
            }
            AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
            AGL_data.xminsav.resize(AGL_data.temperature_external.size());
            AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
            AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
            // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());
            break;
          }
        }
      }
      itry = imin;
      // Fits polynomial to E + F
      if(!LDEBUG) {
        aglerror = AGL_functions::polynomial_fit_weight_ave (imin, helmholtzenergypolynomialcoeffs, helmholtzpolynomialerror, xconfigurationvector, helmholtzenergydatapoints, AGL_data, FileMESSAGE);
      } else {
        aglerror = AGL_functions::polynomial_fit_weight_ave_debug (imin, helmholtzenergypolynomialcoeffs, helmholtzpolynomialerror, xconfigurationvector, helmholtzenergydatapoints, AGL_data, FileMESSAGE);
      }
      // If polynomial fit returns an error for zero temperature, AFLOW AGL exits giving an error
      // If polynomial fit returns an error for T > zero, AFLOW AGL resets the maximum temperature to the previous value and skips the rest of the temperature loop
      // It then continues to complete the remainder of the GIBBS algorithm
      if(aglerror != 0) {
        // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of the iterations of this loop and continues with the next temperature value
        if (AGL_data.run_all_pressure_temperature) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "Problem fitting polynomial to E + F at T = " << AGL_data.temperature_external.at(j) << endl;
          aus << _AGLSTR_WARNING_ + "Skipping this temperature point" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          // OBSOLETE continue;
          minpolyerror = aglerror;
          aglerror = 0;
        } else {
          if(j == 0) {      
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
            aus << _AGLSTR_ERROR_ + "Polynomial fit:  T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
            aus << _AGLSTR_ERROR_ + "Cannot find minimum for temperature = " << AGL_data.temperature_external.at(j) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(aglerror == 2) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Problem inverting matrix to fit polynomial" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return 4;
            } else {
              return 2;
            }
          } else {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
            aus << _AGLSTR_WARNING_ + "Polynomial fit:  T = " << AGL_data.temperature_external.at(j) << "\t" << ", itry = " << itry << "\t" << ", imin = " << imin << "\t" << ", helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
            aus << _AGLSTR_WARNING_ + "Resetting maximum value of T to T = " << AGL_data.temperature_external.at(j-1) << endl;
            aus << _AGLSTR_WARNING_ + "Resetting number of temperature points to " << j << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            AGL_data.temperature_external.resize(j);
            if(AGL_data.savedatapressure) {
              AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
              AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
              AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
              AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
              AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
              AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
              AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
              AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
              AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
              AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
              AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
            }
            AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
            AGL_data.xminsav.resize(AGL_data.temperature_external.size());
            AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
            AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
            // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());
            break;
          }
        }
      }
      if(helmholtzenergypolynomialcoeffs.size() == 3) {
        helmholtzenergypolynomialcoeffs.push_back(0.0);
        helmholtzpolynomialerror.push_back(0.0);
        helmholtzpolynomialerror.push_back(0.0);
      }
      // Loop over pressures
      gibbsenergypolynomialcoeffs.resize(helmholtzenergypolynomialcoeffs.size());
      for (uint i = 0; i < helmholtzenergypolynomialcoeffs.size(); i++) {
        gibbsenergypolynomialcoeffs.at(i) = helmholtzenergypolynomialcoeffs.at(i);
      }
      for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
        aurostd::StringstreamClean(aus);
        if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
          aus << _AGLSTR_MESSAGE_ << "Pressure = " << AGL_data.pressure_external.at(k) << "GPa" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Calculate the value of the Gibbs free energy E + F + pV for each volume AGL_data.volumeinput.at(i) and store in gibbsenergydatapoints
        // Bracket the minimum of E + F + pV  with the numerical function
        // Gibbs free energy = A + pV; constant pressure, variable volume
        if((gibbsenergydatapoints.size() != AGL_data.energyinput.size()) || (gibbsenergydatapoints.size() != F_Helmholtzvinp.size()) || (helmholtzenergydatapoints.size() != AGL_data.volumeinput.size()) ) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
          aus << _AGLSTR_ERROR_ + "gibbsenergydatapoints.size() = " << gibbsenergydatapoints.size() << endl;
          aus << _AGLSTR_ERROR_ + "F_Helmholtzvinp.size() = " << F_Helmholtzvinp.size() << endl;
          aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
          aus << _AGLSTR_ERROR_ + "AGL_data.tdebye.size() = " << AGL_data.tdebye.size() << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 3;
        } else {
          for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
            // OBSOLETE gibbsenergydatapoints.at(i) = AGL_data.energyinput.at(i) + F_Helmholtzvinp.at(i) + AGL_data.volumeinput.at(i) * AGL_data.pressure_external.at(k) / au2GPa;
            gibbsenergydatapoints.at(i) = AGL_data.energyinput.at(i) + F_Helmholtzvinp.at(i) + (AGL_data.volumeinput.at(i) * AGL_data.pressure_external.at(k) / eV_Ang3_to_GPa);
          }
        }
        if(AGL_data.fittype == 0) {
          if (minpolyerror > 0) {
            itry = gibbsenergydatapoints.size() / 2;
          }
          aglerror = AGL_functions::bracket_minimum (itry, gibbsenergydatapoints, FileMESSAGE);
        } else {
          aglerror = AGL_functions::bracket_minimum_global (itry, gibbsenergydatapoints, FileMESSAGE);
        }
        // If bracket_minimum returns an error for zero pressure and zero temperature, AFLOW AGL exits giving an error
        // If bracket_minimum returns an error for p = 0, T > 0, AFLOW AGL resets the maximum temperature to the previous value and skips the rest of the loop
        // If bracket_minimum returns an error for p > zero, AFLOW AGL resets the maximum pressure to the previous value and skips the rest of the pressure loop
        // It then continues to complete the remainder of the GIBBS algorithm
        if(aglerror != 0) {
          // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of the iterations of this loop and continues with the next temperature value
          if (AGL_data.run_all_pressure_temperature) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + pV + F at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
            aus << _AGLSTR_WARNING_ + "Skipping this pressure and temperature point" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            continue;
          } else {	  
            if(k == 0) {
              if(j == 0) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_ERROR_ + "gmin, T = " << AGL_data.temperature_external.at(j) << ", P = " << AGL_data.pressure_external.at(k) << ", trial point = " << itry << ", minimum point = " << imin << ", total points = " << gibbsenergydatapoints.size() << endl;
                aus << _AGLSTR_ERROR_ + "gibbsenergydatapoints = ";
                for (uint i = 0; i < gibbsenergydatapoints.size(); i++) {
                  aus << gibbsenergydatapoints.at(i) << "\t"; 
                }
                aus << endl;
                aus << _AGLSTR_ERROR_ + "Cannot find minimum for pressure = " << AGL_data.pressure_external.at(k)  << " and temperature = " << AGL_data.temperature_external.at(j) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                return 2;
              } else {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_WARNING_ + "gmin, T = " << AGL_data.temperature_external.at(j) << ", P = " << AGL_data.pressure_external.at(k) << ", trial point = " << itry << ", minimum point = " << imin << ", total points = " << gibbsenergydatapoints.size() << endl;
                aus << _AGLSTR_WARNING_ + "gibbsenergydatapoints = ";
                for (uint i = 0; i < gibbsenergydatapoints.size(); i++) {
                  aus << gibbsenergydatapoints.at(i) << "\t";
                }
                aus << endl;
                aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + pV + F at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
                aus << _AGLSTR_WARNING_ + "Resetting maximum value of T to T = " << AGL_data.temperature_external.at(j-1) << endl;
                aus << _AGLSTR_WARNING_ + "Resetting number of pressure points to " << j << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                AGL_data.temperature_external.resize(j);
                if(AGL_data.savedatapressure) {
                  AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                  AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
                  AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
                  AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
                  AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
                  AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                  AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                  AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
                  AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                  AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
                  AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
                }
                AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
                AGL_data.xminsav.resize(AGL_data.temperature_external.size());
                AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
                AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
                // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());
                break;
              }	 
            } else {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ + "gmin, T = " << AGL_data.temperature_external.at(j) << ", P = " << AGL_data.pressure_external.at(k) << ", trial point = " << itry << ", minimum point = " << imin << ", total points = " << gibbsenergydatapoints.size() << endl;
              aus << _AGLSTR_WARNING_ + "gibbsenergydatapoints = ";
              for (uint i = 0; i < gibbsenergydatapoints.size(); i++) {
                aus << gibbsenergydatapoints.at(i) << "\t";
              }
              aus << endl;
              aus << _AGLSTR_WARNING_ + "Cannot find minimum of E + pV + F at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting maximum value of p to p = " << AGL_data.pressure_external.at(k-1) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting number of pressure points to " << k << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              AGL_data.pressure_external.resize(k);
              AGL_data.d2EnergydVolume2_dynamic.resize(AGL_data.pressure_external.size());
              AGL_data.pressure_static.resize(AGL_data.pressure_external.size());
              AGL_data.gamma_G.resize(AGL_data.pressure_external.size());
              AGL_data.pfit.resize(AGL_data.pressure_external.size());
              AGL_data.alpha.resize(AGL_data.pressure_external.size());
              AGL_data.bulkmodulus.resize(AGL_data.pressure_external.size());
              AGL_data.voleqmin.resize(AGL_data.pressure_external.size());
              gibbsenergykjmol.resize(AGL_data.pressure_external.size());
              gibbsenergyeV.resize(AGL_data.pressure_external.size());
              theta.resize(AGL_data.pressure_external.size());
              Cv.resize(AGL_data.pressure_external.size());
              break;	 
            }
          }
        }
        if(LDEBUG) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "bracket_minimum: itry = " << itry << endl;
          aus << _AGLSTR_MESSAGE_ << "bracket_minimum: var(itry) = " << xconfigurationvector.at(itry) << endl;	
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Find the volume which minimizes the fitted G function at (p, T) 
        // Coefficients of fitted polynomial are stored in gibbsenergypolynomialcoeffs)
        // G function is Gibbs free energy: G = E + F + pV; E = DFT energy
        // For a given temperature and pressure, the equilibrium system is the one which minimizes the GIBBS free energy
        // OBSOLETE gibbsenergypolynomialcoeffs.at(3) = helmholtzenergypolynomialcoeffs.at(3) + AGL_data.pressure_external.at(k)/au2GPa * volumereference;
        if (AGL_data.run_all_pressure_temperature) {
          if(!LDEBUG) {
            aglerror = AGL_functions::polynomial_fit_weight_ave (itry, gibbsenergypolynomialcoeffs, gibbspolynomialerror, xconfigurationvector, gibbsenergydatapoints, AGL_data, FileMESSAGE);
          } else {
            aglerror = AGL_functions::polynomial_fit_weight_ave_debug (itry, gibbsenergypolynomialcoeffs, gibbspolynomialerror, xconfigurationvector, gibbsenergydatapoints, AGL_data, FileMESSAGE);
          }
          if (aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Cannot fit polynomial to E + pV + F at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
            aus << _AGLSTR_WARNING_ + "Skipping this pressure and temperature point" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            continue;
          }
        } else {
          gibbsenergypolynomialcoeffs.at(3) = helmholtzenergypolynomialcoeffs.at(3) + (AGL_data.pressure_external.at(k) / eV_Ang3_to_GPa) * volumereference;
        }
        aglerror = AGL_functions::autocorrect_polynom_minimum (itry, xconfigurationvector, gibbsenergypolynomialcoeffs, xmin, FileMESSAGE);
        // If autocorrect_polynom_minimum returns an error for zero temperature and pressure, AFLOW AGL exits giving an error
        // If autocorrect_polynom_minimum returns an error for T > zero, AFLOW AGL resets the maximum temperature or pressure to the previous value and skips the rest of that loop
        // It then continues to complete the remainder of the GIBBS algorithm
        if(pmerr == 0) {
          aurostd::StringstreamClean(aus);
          if(LIB2RAW_FIX_THIS_LATER) {  // to avoid too much garbage in lib2raw
            aus << _AGLSTR_MESSAGE_ << "Minimum of (E, V) data is at point xmin = " << xmin << " Bohr" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
        } else {
          // If the option to avoid truncating the pressure and temperature ranges is set, then skips the rest of this iteration of this loop and continues with the next pressure value
          if (AGL_data.run_all_pressure_temperature) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "Cannot find minimum of polynomial fit for E + pV + F at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
            aus << _AGLSTR_WARNING_ + "Skipping this pressure and temperature point" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            continue;
          } else {
            if(k == 0) {
              if(j == 0) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_ERROR_ + "Cannot find minimum of polynomial fit for E + pV + F for p = " << AGL_data.pressure_external.at(k)  << " and T = " << AGL_data.temperature_external.at(j) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                return 2;
              } else { 
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_WARNING_ + "Cannot find minimum of polynomial fit for E + pV + F at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
                aus << _AGLSTR_WARNING_ + "Resetting maximum temperature value to T = " << AGL_data.temperature_external.at(j-1) << endl;
                aus << _AGLSTR_WARNING_ + "Resetting number of temperature points to " << j << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                AGL_data.temperature_external.resize(j);
                if(AGL_data.savedatapressure) {
                  AGL_data.InternalEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                  AGL_data.EntropyPressurekjmol.resize(AGL_data.temperature_external.size());
                  AGL_data.CvkjmolPressure.resize(AGL_data.temperature_external.size());
                  AGL_data.DebyeTemperaturePressure.resize(AGL_data.temperature_external.size());
                  AGL_data.GruneisenParameterPressure.resize(AGL_data.temperature_external.size());
                  AGL_data.HelmholtzEnergyPressurekjmol.resize(AGL_data.temperature_external.size());
                  AGL_data.InternalEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                  AGL_data.EntropyPressuremeV.resize(AGL_data.temperature_external.size());
                  AGL_data.HelmholtzEnergyPressuremeV.resize(AGL_data.temperature_external.size());
                  AGL_data.CvunitkBpressure.resize(AGL_data.temperature_external.size());
                  AGL_data.EntropyunitkBpressure.resize(AGL_data.temperature_external.size());
                }
                AGL_data.GibbsFreeEnergyPressureeV.resize(AGL_data.temperature_external.size());
                AGL_data.xminsav.resize(AGL_data.temperature_external.size());
                AGL_data.VolumeEquilibrium.resize(AGL_data.temperature_external.size());
                AGL_data.EnergyDFT_UIntVib.resize(AGL_data.temperature_external.size());
                // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.resize(AGL_data.temperature_external.size());
                break;
              }
            } else {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ + "Cannot find minimum of polynomial fit for E + F + pV at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting maximum pressure value to p = " << AGL_data.pressure_external.at(k-1) << endl;
              aus << _AGLSTR_WARNING_ + "Resetting number of pressure points to " << k << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              AGL_data.pressure_external.resize(k);
              break;
            }
          }
        }
        if(LDEBUG) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "xmin = " << xmin << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Evaluates polynomial at minimum to get equilibrium V, G, and B
        AGL_data.voleqmin.at(k) = pow(xmin,3.0) * volumereference;
        AGL_functions::polynom_eval (xmin, gibbsenergypolynomialcoeffs, plnv, 0);
        // OBSOLETE gibbsenergy.at(k) = plnv * hartree2kjmol;
        gibbsenergyeV.at(k) = plnv;
        gibbsenergykjmol.at(k) = plnv * eV2kjmol;
        // OBSOLETE AGL_data.GibbsFreeEnergyPressureeV.at(j).push_back(plnv*hart2ev);
        // If AGL_data.run_all_pressure_temperature is selected, then saves G(p, T), p, T in structure, otherwise saves G(p, T) in a vector of vectors
        if (AGL_data.run_all_pressure_temperature) {
          AGL_pressure_temperature_energy.gibbs_free_energy = plnv;
          AGL_pressure_temperature_energy.pressure_external = AGL_data.pressure_external.at(k);
          AGL_pressure_temperature_energy.temperature_external = AGL_data.temperature_external.at(j);
        } else {
          AGL_data.GibbsFreeEnergyPressureeV.at(j).push_back(plnv);
        }	
        AGL_functions::polynom_eval (xmin, gibbsenergypolynomialcoeffs, plnv, 2);
        // OBSOLETE AGL_data.bulkmodulus.at(k) = plnv * au2GPa * xmin*xmin/(9.0*AGL_data.voleqmin.at(k));
        AGL_data.bulkmodulus.at(k) = plnv * eV_Ang3_to_GPa * xmin * xmin / (9.0 * AGL_data.voleqmin.at(k));
        AGL_functions::polynom_eval (xmin, helmholtzenergypolynomialcoeffs, energyvalue, 0);
        AGL_functions::polynom_eval (xmin, helmholtzpolynomialerror, plnv, 0);
        energyerror = sqrt(aurostd::abs(plnv - energyvalue * energyvalue));
        relative_error.at(k) = energyerror / max (aurostd::abs(energyvalue), aurostd::abs(energyvalue) + energyerror/2.0);     
        // If AGL_data.run_all_pressure_temperature is selected, calculates remaining thermal properties for this temperature and pressure and saves them in structure
        if (AGL_data.run_all_pressure_temperature) {
          AGL_pressure_temperature_energy.volume_equilibrium = AGL_data.voleqmin.at(k);
          aglerror = AGL_functions::numerical_eos_run_all_pT (volumereference, energypolynomialcoeffs, helmholtzenergypolynomialcoeffs, AGL_pressure_temperature_energy.volume_equilibrium, bulkmodpres, pres_static, d2EnergydVolume2_dynamic_pres, Gruneisen_pres);
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            return aglerror;
          }
          // If d2EnergydVolume2_dynamic is not consistent, skip this pressure and go to the next one
          if(d2EnergydVolume2_dynamic_pres <= 0.0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "d^2 E / d V^2 is inconsistent at T = " << AGL_data.temperature_external.at(j) << ", p = " << AGL_data.pressure_external.at(k) << endl;
            aus << _AGLSTR_WARNING_ + "Skipping this temperature and pressure point" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            continue;
          }
          theta_pres = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * AGL_pressure_temperature_energy.volume_equilibrium * AGL_pressure_temperature_energy.volume_equilibrium, third) * AGL_data.poissonratiofunction * sqrt((d2EnergydVolume2_dynamic_pres / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
          AGL_functions::thermal_properties (theta_pres, AGL_data.temperature_external.at(j), DebyeIntegral, DebyeIntegralerror, UIntEnergyVib, Cv_pres, HelmholtzEnergy, Entropy, AGL_data, FileMESSAGE);
          IntEnergyVibmeV = UIntEnergyVib * 1000;
          CvVibUnitkB = Cv_pres * eV2kjmol * 1000 * kj2unit;
          HelmholtzVibmeV = HelmholtzEnergy * 1000;
          EntropyVibkjmol = Entropy * eV2kjmol * 1000;
          EntropyVibmeV = Entropy * 1000;
          EntropyVibUnitkB =  Entropy * eV2kjmol * 1000 * kj2unit;
          // Evaluates polynomial energypolynomialcoeffs to get DFT energy at equilibrium V
          AGL_functions::polynom_eval (xmin, energypolynomialcoeffs, energyvalue, 0);	
          EdftUintVib = energyvalue + UIntEnergyVib;
          AGL_pressure_temperature_energy.internal_energy = UIntEnergyVib;
          AGL_pressure_temperature_energy.E_DFT_internal_vib_energy = EdftUintVib;  
          AGL_pressure_temperature_energy.helmholtz_vib_energy = HelmholtzVibmeV;
          AGL_pressure_temperature_energy.vib_entropy = EntropyVibmeV;
          AGL_data.AGL_pressure_temperature_energy_list.push_back(AGL_pressure_temperature_energy);
        } else {
          AGL_data.xminsav.at(j).push_back(xmin);
          AGL_data.VolumeEquilibrium.at(j).push_back(AGL_data.voleqmin.at(k));
        }
      }
      //  If AGL_data.run_all_pressure_temperature is selected, skips rest of temperature loop for this iteration and moves on to next temperature
      if (AGL_data.run_all_pressure_temperature) {
        continue;
      } else {
        // i_eqn_of_state < 0 means minimum output
        if(AGL_data.i_eqn_of_state < 0) {
          if((AGL_data.pressure_external.size() != AGL_data.voleqmin.size()) || (AGL_data.pressure_external.size() != gibbsenergykjmol.size()) || (AGL_data.pressure_external.size() != AGL_data.bulkmodulus.size()) ) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
            aus << _AGLSTR_ERROR_ + "helmholtzenergydatapoints.size() = " << helmholtzenergydatapoints.size() << endl;
            aus << _AGLSTR_ERROR_ + "F_Helmholtzvinp.size() = " << F_Helmholtzvinp.size() << endl;
            aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
            aus << _AGLSTR_ERROR_ + "AGL_data.tdebye.size() = " << AGL_data.tdebye.size() << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return 3;
          } else {
            for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
              // OBSOLETE AGL_data.outfiless << setw(8) << setprecision(2) << AGL_data.pressure_external.at(k) << "\t" << setw(6) << AGL_data.temperature_external.at(j) << "\t" << setw(10) << AGL_data.voleqmin.at(k) << "\t" << setw(11) << setprecision(5) << gibbsenergy.at(k) << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k) << endl;
              AGL_data.outfiless << setw(8) << setprecision(2) << AGL_data.pressure_external.at(k) << "\t" << setw(6) << AGL_data.temperature_external.at(j) << "\t" << setw(10) << AGL_data.voleqmin.at(k) << "\t" << setw(11) << setprecision(5) << gibbsenergykjmol.at(k) << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k) << endl;
            }
          }
        } else {
          // Numerical energetic, geometric and elastic properties
          volume_0pressure = AGL_data.voleqmin.at(0);
          gibbsenergykjmol_0pressure = gibbsenergykjmol.at(0);
          gibbsenergyeV_0pressure = gibbsenergyeV.at(0);	
          binp_bcnt = AGL_data.bulkmodulus.at(0);
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "Temperature:  T = " << AGL_data.temperature_external.at(j) << "K" << endl;
          // OBSOLETE AGL_data.outfiless << "Vmin(T; P=0) = " << volume_0pressure << "bohr^3" << endl;
          AGL_data.outfiless << "Vmin(T; P=0) = " << volume_0pressure * pow(angstrom2bohr, 3.0) << "bohr^3" << endl;
          AGL_data.outfiless << "Gmin(T; P=0) = " << gibbsenergykjmol_0pressure << "kJ/mol" << endl;
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "NUMERICAL EQUILIBRIUM PROPERTIES" << endl;
          AGL_data.outfiless << "================================" << endl;
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "  P(GPa)         G(kJ/mol)       V(bohr^3)          V/V0          B(GPa)         rel.err. " << endl;
          AGL_data.outfiless << " ----------------------------------------------------------------------------------------- " << endl;
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // OBSOLETE AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "        " << setw(10) << gibbsenergy.at(k) << "      " << setw(10) << AGL_data.voleqmin.at(k) << "      " << setw(8) << setprecision(6) << AGL_data.voleqmin.at(k)/volume_0pressure << "        " << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k) << "        " << setw(9) << setprecision(6) << relative_error.at(k) << endl;
            AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "        " << setw(10) << gibbsenergykjmol.at(k) << "      " << setw(10) << AGL_data.voleqmin.at(k) * pow(angstrom2bohr, 3.0) << "      " << setw(8) << setprecision(6) << AGL_data.voleqmin.at(k)/volume_0pressure << "        " << setw(8) << setprecision(2) << AGL_data.bulkmodulus.at(k) << "        " << setw(9) << setprecision(6) << relative_error.at(k) << endl;
          }
          // Finite temperature results using numerical method of calculating bulk modulus from E(V) data.
          if(AGL_data.i_eqn_of_state == 0) {
            aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, helmholtzenergypolynomialcoeffs, xconfigurationvector, false, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              return aglerror;
            }
          }
          // Finite temperature results from Vinet equation of state.
          else if(AGL_data.i_eqn_of_state == 1) {
            // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
            // OBSOLETE ge0peV = gibbsenergykjmol_0pressure / eV2kjmol; 
            // OBSOLETE aglerror = AGL_functions::vinet_eos (volume_0pressure, g0dhy2kjm, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::vinet_eos (volume_0pressure, gibbsenergyeV_0pressure, false, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              return aglerror;
            }
          }
          // Finite temperature results from Birch-Murnaghan equation of state.
          else if(AGL_data.i_eqn_of_state == 2) {
            // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
            // OBSOLETE ge0peV = gibbsenergykjmol_0pressure / eV2kjmol; 
            // OBSOLETE aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, g0dhy2kjm, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, gibbsenergyeV_0pressure, false, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in Birch-Murnaghan fit" << endl;
              aus << _AGLSTR_ERROR_ + "Check input pressure values (first pressure point must be p = 0)" << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              return aglerror;
            }
          }
          // Finite temperature results using Vinet equation of state, but calculating other properties numerically.
          else if(AGL_data.i_eqn_of_state == 3) {
            // OBSOLETE g0dhy2kjm = g0dhy2kjm;
            // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
            // OBSOLETE ge0peV = gibbsenergykjmol_0pressure / eV2kjmol; 
            // OBSOLETE aglerror = AGL_functions::vinet_eos (volume_0pressure, g0dhy2kjm, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::vinet_eos (volume_0pressure, gibbsenergyeV_0pressure, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, helmholtzenergypolynomialcoeffs, xconfigurationvector, false, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              return aglerror;
            }
          }
          // Finite temperature results from Birch-Murnaghan equation of state, but calculating other properties numerically.
          else if(AGL_data.i_eqn_of_state == 4) {
            // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
            // OBSOLETE ge0peV = gibbsenergykjmol_0pressure / eV2kjmol; 
            // OBSOLETE aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, g0dhy2kjm, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::birch_murnaghan_eos (volume_0pressure, gibbsenergyeV_0pressure, false, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in Birch-Murnaghan fit" << endl;
              aus << _AGLSTR_ERROR_ + "Check input pressure values (first pressure point must be p = 0)" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }
            aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, helmholtzenergypolynomialcoeffs, xconfigurationvector, false, AGL_data, FileMESSAGE);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }
          }
          // Equation of state of Baonza-Caceres-Nuñez-Taravillo.
          else if(AGL_data.i_eqn_of_state == 5) {
            // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
            // OBSOLETE ge0peV = gibbsenergykjmol_0pressure / eV2kjmol; 
            // OBSOLETE aglerror = AGL_functions::bcnt_eos (volume_0pressure, g0dhy2kjm, binp_bcnt, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::bcnt_eos (volume_0pressure, gibbsenergyeV_0pressure, binp_bcnt, false, AGL_data, FileMESSAGE);
            bcnt_K_sp.push_back(AGL_data.xsup_K_final);
            bcnt_P_sp.push_back(AGL_data.Press_sp_final);
            bcnt_V_sp.push_back(AGL_data.Vol_sp_final);
            bcnt_beta_sp.push_back(AGL_data.bcnt_beta_final);
            if(aglerror != 0) {
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              return aglerror;
            }		  
          }
          // Equation of state of Baonza-Caceres-Nuñez-Taravillo, but numerical calculation to obtain all properties.
          else if(AGL_data.i_eqn_of_state == 6) {
            // OBSOLETE g0dhy2kjm = gibbsenergy_0pressure/hartree2kjmol;
            // OBSOLETE ge0peV = gibbsenergykjmol_0pressure / eV2kjmol; 
            // OBSOLETE aglerror = AGL_functions::bcnt_eos (volume_0pressure, g0dhy2kjm, binp_bcnt, false, AGL_data, FileMESSAGE);
            aglerror = AGL_functions::bcnt_eos (volume_0pressure, gibbsenergyeV_0pressure, binp_bcnt, false, AGL_data, FileMESSAGE);
            bcnt_K_sp.push_back(AGL_data.xsup_K_final);
            bcnt_P_sp.push_back(AGL_data.Press_sp_final);
            bcnt_V_sp.push_back(AGL_data.Vol_sp_final);
            bcnt_beta_sp.push_back(AGL_data.bcnt_beta_final);		  
            aglerror = AGL_functions::numerical_eos (volumereference, energypolynomialcoeffs, helmholtzenergypolynomialcoeffs, xconfigurationvector, false, AGL_data, FileMESSAGE);
          }
          if(aglerror != 0) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_ERROR_ + "Error in equation of state fit" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          // Save properties at P=0 for all of the temperatures.
          volumetemperature0pressure.push_back(AGL_data.voleqmin.at(0));
          gibbsenergytemperature0pressurekjmol.push_back(gibbsenergykjmol.at(0));
          AGL_data.bulkmodulusisothermal_0pressure.push_back(AGL_data.bulkmodulus_0pressure);
          dbulkmodulusdpt0pressure.push_back(AGL_data.dbulkmodulusdpV_0pressure);
          d2bulkmodulusdp2t0pressure.push_back(AGL_data.d2bulkmodulusdpV2_0pressure);
          AGL_data.GibbsFreeEnergy0pressurekjmol.push_back(gibbsenergykjmol.at(0));
          // OBSOLETE gibbsenergy0pressureeV = (gibbsenergy.at(0) / hartree2kjmol) * hart2ev;
          // OBSOLETE gibbsenergy0pressureeV = gibbsenergykjmol.at(0) / eV2kjmol;
          // OBSOLETE AGL_data.GibbsFreeEnergy0pressureeV.push_back(gibbsenergy0pressureeV);gibbsenergyeV_0pressure
          AGL_data.GibbsFreeEnergy0pressureeV.push_back(gibbsenergyeV.at(0));
          // Compute vibrational properties at the equilibrium volumes.
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "VIBRATIONAL PROPERTIES" << endl;
          AGL_data.outfiless << "======================" << endl;
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "  P(GPa)          U(kJ/mol)      Cv(J/mol*K)     A(kJ/mol)       S(J/mol*K)      Theta(K)        gamma" << endl;
          AGL_data.outfiless << " ------------------------------------------------------------------------------------------------------ " << endl;
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            if(AGL_data.i_debye == 1 || AGL_data.i_debye == 2) {
              if((theta.size() != AGL_data.pressure_external.size()) || (Cv.size() != AGL_data.pressure_external.size()) || (AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) ) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
                aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
                aus << _AGLSTR_ERROR_ + "theta.size() = " << theta.size() << endl;
                aus << _AGLSTR_ERROR_ + "Cv.size() = " << Cv.size() << endl;
                aus << _AGLSTR_ERROR_ + "AGL_data.voleqmin.size() = " << AGL_data.voleqmin.size() << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                return 3;
              } 	    
              logvol = log (AGL_data.voleqmin.at(k));
              aglerror = AGL_functions::polynom_eval (logvol, tdebyepolynomialcoeffs, plnv, 0);
              theta.at(k) = exp ( plnv );
              aglerror = AGL_functions::polynom_eval (logvol, tdebyepolynomialcoeffs, plnv, 1);
              AGL_data.gamma_G.at(k) = - plnv;
              aglerror = AGL_functions::thermal_properties (theta.at(k), AGL_data.temperature_external.at(j), DebyeIntegral, DebyeIntegralerror, UIntEnergyVib, Cv.at(k), HelmholtzEnergy, Entropy, AGL_data, FileMESSAGE);
              // OBSOLETE IntEnergyVib = UIntEnergyVib * hartree2kjmol;
              IntEnergyVibkjmol = UIntEnergyVib * eV2kjmol;
              // OBSOLETE IntEnergyVibmeV = UIntEnergyVib * hart2ev * 1000;
              IntEnergyVibmeV = UIntEnergyVib * 1000;
              // OBSOLETE CvVib = Cv.at(k) * hartree2kjmol * 1000;
              CvVibkjmol = Cv.at(k) * eV2kjmol * 1000;
              CvVibUnitkB = CvVibkjmol * kj2unit;
              // OBSOLETE HelmholtzVib = HelmholtzEnergy * hartree2kjmol;
              HelmholtzVibkjmol = HelmholtzEnergy * eV2kjmol;
              // OBSOLETE HelmholtzVibmeV = HelmholtzEnergy * hart2ev * 1000;
              HelmholtzVibmeV = HelmholtzEnergy * 1000;
              // OBSOLETE EntropyVib = Entropy * hartree2kjmol * 1000;
              EntropyVibkjmol = Entropy * eV2kjmol * 1000;
              // OBSOLETE EntropyVibmeV = Entropy * hart2ev * 1000;
              EntropyVibmeV = Entropy * 1000;
              EntropyVibUnitkB = EntropyVibkjmol * kj2unit;
              // Evaluates polynomial energypolynomialcoeffs to get DFT energy at equilibrium V
              AGL_functions::polynom_eval (AGL_data.xminsav.at(j).at(k), energypolynomialcoeffs, energyvalue, 0);	
              EdftUintVib = energyvalue + UIntEnergyVib;
              // OBSOLETE EdftUintVibeV = (energyvalue + UIntEnergyVib) * hart2ev;
              AGL_data.EnergyDFT_UIntVib.at(j).push_back(EdftUintVib);
              // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.at(j).push_back(EdftUintVibeV);
            } else {
              if(k == 0) {
                aglerror = AGL_functions::debye_polynom_fit (tdebyepolynomialcoeffs, AGL_data, FileMESSAGE);
                if(aglerror != 0) {
                  if(aglerror == 2) {
                    aurostd::StringstreamClean(aus);
                    aus << _AGLSTR_ERROR_ + "Problem inverting matrix to fit polynomial" << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    return 4;
                  } else {
                    aurostd::StringstreamClean(aus);
                    aus << _AGLSTR_ERROR_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    return 2;
                  }
                }
                logvol = log (AGL_data.voleqmin.at(k));
                aglerror = AGL_functions::polynom_eval (logvol, tdebyepolynomialcoeffs, plnv, 1);
                AGL_data.gamma_poly.push_back(- plnv);
              }
              if((theta.size() != AGL_data.pressure_external.size()) || (Cv.size() != AGL_data.pressure_external.size()) || (AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) || (AGL_data.d2EnergydVolume2_dynamic.size() != AGL_data.pressure_external.size()) ) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
                aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
                aus << _AGLSTR_ERROR_ + "theta.size() = " << theta.size() << endl;
                aus << _AGLSTR_ERROR_ + "Cv.size() = " << Cv.size() << endl;
                aus << _AGLSTR_ERROR_ + "AGL_data.voleqmin.size() = " << AGL_data.voleqmin.size() << endl;
                aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_dynamic.size() = " << AGL_data.d2EnergydVolume2_dynamic.size() << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                return 3;
              } 	    
              // If d2EnergydVolume2_dynamic is not consistent, make a Gruneisen interpolation
              if(AGL_data.d2EnergydVolume2_dynamic.at(k) <= 0.0) {
                aurostd::StringstreamClean(aus);
                aus << _AGLSTR_WARNING_ + "inconsistent derivative, v[k] = " << AGL_data.voleqmin.at(k) << ", p[k] = " << AGL_data.pressure_external.at(k) << endl;
                aus << _AGLSTR_WARNING_ + "k = " << k << ", T = " << AGL_data.temperature_external.at(j) << ", d2EnergydVolume2_dynamic[k] = " << AGL_data.d2EnergydVolume2_dynamic.at(k) << endl; 
                aus << _AGLSTR_WARNING_ + "making a Gruneisen interpolation" << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                ilow = 0;
                while (AGL_data.voleqmin.at(k) >= AGL_data.volumeinput.at(ilow) && ilow < (AGL_data.volumeinput.size()-1)) {
                  ilow = ilow + 1;
                }
                if((ilow == (AGL_data.volumeinput.size()-1) && AGL_data.voleqmin.at(k) >= AGL_data.volumeinput.at(ilow)) || ilow == 0) {
                  if(k == 0) {
                    aurostd::StringstreamClean(aus);
                    aus << _AGLSTR_ERROR_ + "inconsistent volume, v[k] = " << AGL_data.voleqmin.at(k) << ", p[k] = " << AGL_data.pressure_external.at(k) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    return 5;
                  } else {
                    aurostd::StringstreamClean(aus);
                    aus << _AGLSTR_WARNING_ + "Inconsistent volume, v[k] = " << AGL_data.voleqmin.at(k) << ", p[k] = " << AGL_data.pressure_external.at(k) << ", k = " << k << endl;
                    aus << _AGLSTR_WARNING_ + "Resetting maximum pressure value to p = " << AGL_data.pressure_external.at(k-1) << endl;
                    aus << _AGLSTR_WARNING_ + "Resetting number of pressure points to " << k << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    AGL_data.pressure_external.resize(k);
                    AGL_data.d2EnergydVolume2_dynamic.resize(AGL_data.pressure_external.size());
                    AGL_data.pressure_static.resize(AGL_data.pressure_external.size());
                    AGL_data.gamma_G.resize(AGL_data.pressure_external.size());
                    AGL_data.pfit.resize(AGL_data.pressure_external.size());
                    AGL_data.alpha.resize(AGL_data.pressure_external.size());
                    AGL_data.bulkmodulus.resize(AGL_data.pressure_external.size());
                    AGL_data.voleqmin.resize(AGL_data.pressure_external.size());
                    gibbsenergykjmol.resize(AGL_data.pressure_external.size());
                    theta.resize(AGL_data.pressure_external.size());
                    Cv.resize(AGL_data.pressure_external.size());
                    break;
                  }
                }
                ilow = ilow - 1;
                AGL_data.d2EnergydVolume2_dynamic.at(k) = AGL_data.d2EnergydVolume2_static.at(ilow) * pow(AGL_data.voleqmin.at(k)/AGL_data.volumeinput.at(ilow), log(AGL_data.d2EnergydVolume2_static.at(ilow+1)/AGL_data.d2EnergydVolume2_static.at(ilow)) / log(AGL_data.volumeinput.at(ilow+1)/AGL_data.volumeinput.at(ilow)) );
              }
              // Isotropic Debye model properties
              // OBSOLETE theta.at(k) = pow(6 * pi * pi * AGL_data.natoms * AGL_data.voleqmin.at(k) * AGL_data.voleqmin.at(k), third) / physconstkbau * AGL_data.poissonratiofunction * sqrt(AGL_data.d2EnergydVolume2_dynamic.at(k)/AGL_data.cellmass);
              // OBSOLETE theta.at(k) = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * sqrt(AGL_data.voleqmin.at(k)), third) * AGL_data.poissonratiofunction * sqrt((aurostd::abs(AGL_data.d2EnergydVolume2_dynamic.at(k) * AGL_data.voleqmin.at(k)) / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
              theta.at(k) = (hbar_eV / KBOLTZEV) * pow(6 * pi * pi * AGL_data.natoms * AGL_data.voleqmin.at(k) * AGL_data.voleqmin.at(k), third) * AGL_data.poissonratiofunction * sqrt((AGL_data.d2EnergydVolume2_dynamic.at(k) / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
              AGL_functions::thermal_properties (theta.at(k), AGL_data.temperature_external.at(j), DebyeIntegral, DebyeIntegralerror, UIntEnergyVib, Cv.at(k), HelmholtzEnergy, Entropy, AGL_data, FileMESSAGE);
              // OBSOLETE IntEnergyVib = UIntEnergyVib * hartree2kjmol;
              IntEnergyVibkjmol = UIntEnergyVib * eV2kjmol;
              // OBSOLETE IntEnergyVibmeV = UIntEnergyVib * hart2ev * 1000;
              IntEnergyVibmeV = UIntEnergyVib * 1000;
              // OBSOLETE CvVib = Cv.at(k) * hartree2kjmol * 1000;
              CvVibkjmol = Cv.at(k) * eV2kjmol * 1000;
              CvVibUnitkB = CvVibkjmol * kj2unit;
              // OBSOLETE HelmholtzVib = HelmholtzEnergy * hartree2kjmol;
              HelmholtzVibkjmol = HelmholtzEnergy * eV2kjmol;
              // OBSOLETE HelmholtzVibmeV = HelmholtzEnergy * hart2ev * 1000;
              HelmholtzVibmeV = HelmholtzEnergy * 1000;
              // OBSOLETE EntropyVib = Entropy * hartree2kjmol * 1000;
              EntropyVibkjmol = Entropy * eV2kjmol * 1000;
              // OBSOLETE EntropyVibmeV = Entropy * hart2ev * 1000;
              EntropyVibmeV = Entropy * 1000;
              EntropyVibUnitkB = EntropyVibkjmol * kj2unit;
              // Evaluates polynomial energypolynomialcoeffs to get DFT energy at equilibrium V
              AGL_functions::polynom_eval (AGL_data.xminsav.at(j).at(k), energypolynomialcoeffs, energyvalue, 0);	
              EdftUintVib = energyvalue + UIntEnergyVib;
              // OBSOLETE EdftUintVibeV = (energyvalue + UIntEnergyVib) * hart2ev;
              AGL_data.EnergyDFT_UIntVib.at(j).push_back(EdftUintVib);
              // OBSOLETE AGL_data.EnergyDFT_UIntVibeV.at(j).push_back(EdftUintVibeV);
            }
            AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "        " << setw(11) << setprecision(5) << IntEnergyVibkjmol << "     " << setw(12) << setprecision(6) << CvVibkjmol << "    " << setw(10) << setprecision(4) << HelmholtzVibkjmol << "      " << setw(11) << setprecision(6) << EntropyVibkjmol << "     " << setw(9) << setprecision(2) << theta.at(k) << "       " << setw(6) << setprecision(3) << AGL_data.gamma_G.at(k) << endl;
            if(k == 0) {
              // Save heat capacity, internal energy, entropy, Debye temperature, Gruneisen parameter, Helmholtz free energy and Gibbs free energy at zero pressure for analysis and plotting later on
              AGL_data.Cvkjmol0pressure.push_back(CvVibkjmol);
              AGL_data.CvunitkB0pressure.push_back(CvVibUnitkB);
              AGL_data.InternalEnergy0pressurekjmol.push_back(IntEnergyVibkjmol);
              AGL_data.InternalEnergy0pressuremeV.push_back(IntEnergyVibmeV);
              AGL_data.Entropy0pressurekjmol.push_back(EntropyVibkjmol);
              AGL_data.Entropy0pressuremeV.push_back(EntropyVibmeV);
              AGL_data.Entropy0pressureunitkB.push_back(EntropyVibUnitkB);
              AGL_data.DebyeTemperature0pressure.push_back(theta.at(k));
              AGL_data.GruneisenParameter0pressure.push_back(AGL_data.gamma_G.at(k));
              AGL_data.HelmholtzEnergy0pressurekjmol.push_back(HelmholtzVibkjmol);
              AGL_data.HelmholtzEnergy0pressuremeV.push_back(HelmholtzVibmeV);
            }
            if(AGL_data.savedatapressure) {
              // Save heat capacity, internal energy, entropy, Debye temperature, Gruneisen parameter, and Helmholtz free energy at all pressures for analysis and plotting later on
              // Saving these values is optional and can be activated by setting [AFLOW_GIBBS]SAVEALLPRESSURES=ON in the _AFLOWIN_ file
              // Note that this can use a lot of memory, but is useful if you want to generate (p, T) phase diagrams
              AGL_data.CvkjmolPressure.at(j).push_back(CvVibkjmol);
              AGL_data.CvunitkBpressure.at(j).push_back(CvVibUnitkB);
              AGL_data.InternalEnergyPressurekjmol.at(j).push_back(IntEnergyVibkjmol);
              AGL_data.InternalEnergyPressuremeV.at(j).push_back(IntEnergyVibmeV);
              AGL_data.EntropyPressurekjmol.at(j).push_back(EntropyVibkjmol);
              AGL_data.EntropyPressuremeV.at(j).push_back(EntropyVibmeV);
              AGL_data.EntropyunitkBpressure.at(j).push_back(EntropyVibUnitkB);
              AGL_data.DebyeTemperaturePressure.at(j).push_back(theta.at(k));
              AGL_data.GruneisenParameterPressure.at(j).push_back(AGL_data.gamma_G.at(k));
              AGL_data.HelmholtzEnergyPressurekjmol.at(j).push_back(HelmholtzVibkjmol);
              AGL_data.HelmholtzEnergyPressuremeV.at(j).push_back(HelmholtzVibmeV);
            }
          }
          // Compute Equation of State derivatives
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "THERMAL EOS DERIVATIVES" << endl;
          AGL_data.outfiless << "=======================" << endl;
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "  P(GPa)         alpha(10^{-5}/K)   dp/dt(GPa/K)    Bs(GPa)     Cp(J/mol*K) " << endl;
          AGL_data.outfiless << " --------------------------------------------------------------------------- " << endl;
          // OBSOLETE AGL_data.outfiless << endl;
          for (uint k = 0; k < AGL_data.pressure_external.size(); k++) {
            // OBSOLETE Pbeta = Cv.at(k) * AGL_data.gamma_G.at(k) / AGL_data.voleqmin.at(k) * au2GPa;
            Pbeta = Cv.at(k) * AGL_data.gamma_G.at(k) / AGL_data.voleqmin.at(k) * eV_Ang3_to_GPa;
            AGL_data.alpha.at(k) = Pbeta / AGL_data.bulkmodulus.at(k);
            gammaalphaT = 1.0 + AGL_data.gamma_G.at(k) * AGL_data.alpha.at(k) * AGL_data.temperature_external.at(j);
            // OBSOLETE Cp = Cv.at(k) * gammaalphaT * hartree2kjmol * 1000;
            Cpkjmol = Cv.at(k) * gammaalphaT * eV2kjmol * 1000;
            CpUnitkB = Cpkjmol * kj2unit;
            Bstatic = AGL_data.bulkmodulus.at(k) * gammaalphaT;
            AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.pressure_external.at(k) << "        " << setw(17) << setprecision(8) << AGL_data.alpha.at(k)*1e5 << "  " << setw(13) << setprecision(9) << Pbeta << "   " << setw(8) << setprecision(2) << Bstatic << "    " << setw(12) << setprecision(6) << Cpkjmol << endl;
            if(k == 0) {
              AGL_data.ThermalExpansion0pressure.push_back(AGL_data.alpha.at(k));
              AGL_data.bulkmodulusstatic_0pressure.push_back(Bstatic);
              AGL_data.Cpkjmol0pressure.push_back(Cpkjmol);
              AGL_data.CpunitkB0pressure.push_back(CpUnitkB);
            }
          }   
          AGL_data.outfiletss << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << setprecision(8) << setw(15) << AGL_data.InternalEnergy0pressuremeV.at(j) << " " << setw(15) << AGL_data.HelmholtzEnergy0pressuremeV.at(j) << " " << setw(15) << AGL_data.Entropy0pressureunitkB.at(j) << " " << setw(15) << AGL_data.CvunitkB0pressure.at(j) << " " << setw(15) << AGL_data.DebyeTemperature0pressure.at(j) << " " << setw(23) << AGL_data.GruneisenParameter0pressure.at(j) << endl;
        }
      }
      // Close temperature loop
    }
    // If AGL_data.run_all_pressure_temperature is selected, skips writing out zero pressure properties, since these are not necessarily present for all temperatures
    if (AGL_data.run_all_pressure_temperature) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "End of run ok!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    } else { 
      // Write all of the results at P=0 for all of the temperatures to the stringstream for the main GIBBS output file.
      if(AGL_data.i_eqn_of_state >= 0) {
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "RESULTS AT P=0 FOR ALL TEMPERATURES" << endl;
        AGL_data.outfiless << "===================================" << endl;
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "    T(K)         V(bohr^3)       G(kJ/mol)       U(kJ/mol)       S(J/mol K)      Cv(J/mol K) " << endl;
        AGL_data.outfiless << " -------------------------------------------------------------------------------------------- " << endl;
        for(uint j = 0; j < AGL_data.temperature_external.size(); j++) {
          // OBSOLETE AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << "        " << setw(10) << setprecision(4) << volumetemperature0pressure.at(j) << "      " << setw(10) << setprecision(4) << gibbsenergytemperature0pressure.at(j) << "      " << setw(10) << setprecision(4) << AGL_data.InternalEnergy0pressuremeV.at(j) << "      " << setw(11) << setprecision(5) << AGL_data.Entropy0pressure.at(j) << "     " << setw(12) << setprecision(6) << AGL_data.Cv0pressure.at(j) << endl;
          AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << "        " << setw(10) << setprecision(4) << (volumetemperature0pressure.at(j) * pow(angstrom2bohr, 3.0)) << "      " << setw(10) << setprecision(4) << gibbsenergytemperature0pressurekjmol.at(j) << "      " << setw(10) << setprecision(4) << AGL_data.InternalEnergy0pressurekjmol.at(j) << "      " << setw(11) << setprecision(5) << AGL_data.Entropy0pressurekjmol.at(j) << "     " << setw(12) << setprecision(6) << AGL_data.Cvkjmol0pressure.at(j) << endl;
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "OTHER THERMODYNAMIC PROPERTIES AT P=0" << endl;
        AGL_data.outfiless << "=====================================" << endl;
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "    T(K)         B0(GPa)             B0'         B0''(GPa-1)     Bs(GPa)         alpha(10^{-5}/K) " << endl;  
        AGL_data.outfiless << " ------------------------------------------------------------------------------------------------- " << endl;
        for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
          AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << "        " << setw(8) << setprecision(2) << AGL_data.bulkmodulusisothermal_0pressure.at(j) << "        " << setw(8) << setprecision(5) << dbulkmodulusdpt0pressure.at(j) << "        " << setw(12) << setprecision(6) << d2bulkmodulusdp2t0pressure.at(j) << "    " << setw(8) << setprecision(2) << AGL_data.bulkmodulusstatic_0pressure.at(j) << "        " << setw(17) << setprecision(5) << AGL_data.ThermalExpansion0pressure.at(j)*100000.0 << endl;
        }
        if(AGL_data.i_eqn_of_state == 5 || AGL_data.i_eqn_of_state == 6) {
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "SPINODAL PROPERTIES" << endl;
          AGL_data.outfiless << "===================" << endl;
          AGL_data.outfiless << endl;
          AGL_data.outfiless << "    T(K) \t     K* \t Psp(GPa) \t Vsp(bohr^3) \t    beta " << endl;
          AGL_data.outfiless << " ------------------------------------------------------------------------ " << endl;
          for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
            AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << "\t" << setw(8) << setprecision(6) << bcnt_K_sp.at(j) << "\t" << setw(9) << setprecision(2) << -bcnt_P_sp.at(j) << "\t" << setw(12) << setprecision(6) << bcnt_V_sp.at(j) << "\t" << setw(8) << setprecision(5) << bcnt_beta_sp.at(j) << endl;
          }
        }
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "DEBYE MODEL RELATED PROPERTIES AT P=0" << endl;
        AGL_data.outfiless << "=====================================" << endl;
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "    T(K)         Theta(K)        gamma " << endl;
        AGL_data.outfiless << " -------------------------------------- " << endl;
        for (uint j = 0; j < AGL_data.temperature_external.size(); j++) {
          AGL_data.outfiless << setw(8) << setprecision(2) << fixed << AGL_data.temperature_external.at(j) << "        " << setw(9) << setprecision(2) << AGL_data.DebyeTemperature0pressure.at(j) << "       " << setw(6) << setprecision(2) << AGL_data.GruneisenParameter0pressure.at(j) << endl;
        }
      }
    }
    // End of the GIBBS program
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "End of run ok!" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    return 0;
  }  
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL run GIBBS
// **************************************************************************

#endif // _AFLOW_AGL_RUNGIBBS_CPP


// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
