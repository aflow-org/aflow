// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_THERMAL_CPP
#define _AFLOW_AGL_THERMAL_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"

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


// **************************************************************************************
//  This set of functions calculate thermal properties 
//  and determine the Debye temperature self-consistently
// **************************************************************************************

// ***************************************************************************
// AGL_functions::self_consistent_debye
// ***************************************************************************
namespace AGL_functions {
  //
  // self_consistent_debye: calculates the table of self consistent Debye
  //     temperatures at T for all the input volumes.
  //
  //     If firsttime is true, only calculates the static pressures table
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint self_consistent_debye (double& Temperature, vector<double>& polynomialcoeffs, vector<double>& polynomialerror, vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, uint& imin, bool firsttime, vector<double>& press_stat, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    double eps = 5.0;
    int mloops = 25;
    bool converged = false;
    double theta0 = 0.0, theta = 0.0;
    double DebyeIntegral = 0.0, DebyeIntegralerror = 0.0, UIntEnergy = 0.0, Cvt = 0.0, F_Helmholtz = 0.0, S_Entropy = 0.0;
    uint itry = 0;
    uint ierr = 0;
    uint aglerror = 0;
    double dFdx = 0.0, d2Fdx2 = 0.0, plnv = 0.0;
    double pressure_F = 0.0, bulkmodT = 0.0, mie_gruneisen = 0.0, bulkmodstV = 0.0, delta_theta = 0.0;
    ostringstream aus;
    // If in static part of algorithm (before temperature loop), then calculates static pressures
    if(firsttime) {
      for (uint i = 0; i < xdata_to_fit.size(); i++) {
        aglerror = polynom_eval (xdata_to_fit.at(i), polynomialcoeffs, plnv, 1);
        press_stat.at(i) = -xdata_to_fit.at(i) * plnv / (3.0 * AGL_data.volumeinput.at(i));
      }
      return aglerror;
    }
    // Iterate until convergence is achieved
    converged = false;
    int iloops = 0;
    itry = imin;
    if(ydata_to_fit.size() != AGL_data.energyinput.size()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
      aus << _AGLSTR_ERROR_ + "ydata_to_fit.size() = " << ydata_to_fit.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 3;
    }
    while ((!converged || (iloops < mloops)) && (iloops < AGL_data.maxloops)) {
      iloops = iloops + 1;
      for (uint i = 0; i < ydata_to_fit.size(); i++) {
        thermal_properties (AGL_data.tdebye.at(i), Temperature, DebyeIntegral, DebyeIntegralerror, UIntEnergy, Cvt, F_Helmholtz, S_Entropy, AGL_data, FileMESSAGE);
        ydata_to_fit.at(i) = AGL_data.energyinput.at(i) + F_Helmholtz;
      }
      if(AGL_data.fittype == 0) {
        ierr = AGL_functions::bracket_minimum (imin, ydata_to_fit, FileMESSAGE);
      } else {
        ierr = AGL_functions::bracket_minimum_global (imin, ydata_to_fit, FileMESSAGE);
      }
      if(ierr != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "self_consistent_debye: T = " << Temperature << ", minimum point = " << imin << ", trial point = " << itry << ", total points = " << ydata_to_fit.size() << endl;
        aus << _AGLSTR_WARNING_ + "self_consistent_debye: Helmholtz energy points to fit = ";
        for(uint i = 0; i < ydata_to_fit.size(); i++) {
          aus << ydata_to_fit.at(i) << "\t";
        }
        aus << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
        return 1;
      }
      aglerror = polynomial_fit_weight_ave (imin, polynomialcoeffs, polynomialerror, xdata_to_fit, ydata_to_fit, AGL_data, FileMESSAGE);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "No polynomial fit found with a minimum within bounds of input data" << endl;
        aus << _AGLSTR_WARNING_ + "self_consistent_debye: T = " << Temperature << ", minimum point = " << imin << ", trial point = " << itry << ", total points = " << ydata_to_fit.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
        return 1;
      }
      // Calculate the new set of debye temperatures and its convergence
      converged = true;
      theta0 = AGL_data.tdebye.at(0);
      if((AGL_data.volumeinput.size() != xdata_to_fit.size()) || (AGL_data.tdebye.size() != xdata_to_fit.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.tdebye.size() = " << AGL_data.tdebye.size() << endl;
        aus << _AGLSTR_ERROR_ + "xdata_to_fit.size() = " << xdata_to_fit.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 3;
      }
      for (uint i=0; i < xdata_to_fit.size(); i++) {
        aglerror = AGL_functions::thermal_properties (AGL_data.tdebye.at(i), Temperature, DebyeIntegral, DebyeIntegralerror, UIntEnergy, Cvt, F_Helmholtz, S_Entropy, AGL_data, FileMESSAGE);
        aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(i), polynomialcoeffs, dFdx, 1);
        aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(i), polynomialcoeffs, d2Fdx2, 2);
        pressure_F = -xdata_to_fit.at(i) * dFdx / (3.0 * AGL_data.volumeinput.at(i));
        bulkmodT = -xdata_to_fit.at(i) / (9.0 * AGL_data.volumeinput.at(i)) * (2.0 * dFdx - xdata_to_fit.at(i) * d2Fdx2);
        mie_gruneisen = (pressure_F - press_stat.at(i)) * AGL_data.volumeinput.at(i) / UIntEnergy;
        bulkmodstV = AGL_data.volumeinput.at(i) * bulkmodT + Temperature * mie_gruneisen * mie_gruneisen * Cvt;
        if(bulkmodstV < 0.0) {
          if(i <= imin) {
            aurostd::StringstreamClean(aus);
            aus << _AGLSTR_WARNING_ + "self_consistent_debye: Bstatic < 0 for equilibrium V!" << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            return aglerror;
          }
          if(i > 0) {
            theta = AGL_data.tdebye.at(i) / theta0 * AGL_data.tdebye.at(i-1);
          }
          delta_theta = 0.0;
        } else {
          // OBSOLETE theta = pow((6.0 * pi * pi * AGL_data.natoms / AGL_data.volumeinput.at(i)), third) / physconstkbau * AGL_data.poissonratiofunction * sqrt(bulkmodstV / AGL_data.cellmass);
          // OBSOLETE theta = (hbar_eV / KBOLTZEV) * pow((6.0 * pi * pi * AGL_data.natoms * sqrt(AGL_data.volumeinput.at(i))), third) * AGL_data.poissonratiofunction * sqrt((aurostd::abs(bulkmodstV / AGL_data.volumeinput.at(i)) / AGL_data.cellmass)  * eV_Ang3_to_amu_Ang_s);
          theta = (hbar_eV / KBOLTZEV) * pow((6.0 * pi * pi * AGL_data.natoms / AGL_data.volumeinput.at(i)), third) * AGL_data.poissonratiofunction * sqrt((bulkmodstV / AGL_data.cellmass) * eV_Ang3_to_amu_Ang_s);
          delta_theta = theta - AGL_data.tdebye.at(i);
          if(i > 0) {
            if(theta > AGL_data.tdebye.at(i-1)) {
              theta = AGL_data.tdebye.at(i)/theta0 * AGL_data.tdebye.at(i-1);
              aurostd::StringstreamClean(aus);
              aus << _AGLSTR_WARNING_ + "self_consistent_debye: Warning! gamma < 0, i = " << i << ", T = " << Temperature << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              delta_theta = 0.0;
            }
          }
        }
        theta0 = AGL_data.tdebye.at(i);
        AGL_data.tdebye.at(i) = theta;
        if((AGL_data.volumeinput.at(i) - AGL_data.volumeinput.at(imin) < AGL_data.volumeinput.at(imin) * 0.1) && (i >= 3 && i >= (AGL_data.volumeinput.size()/10))) {
          converged = converged && (aurostd::abs(delta_theta) < eps);
        }
      }	
    }

    // Warn of convergence failure
    if(iloops >= AGL_data.maxloops && !converged) {    
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "self_consistent_debye: Maximum number of convergence iterations exceeded" << endl;
      aus << _AGLSTR_WARNING_ + "self_consistent_debye: Number of convergence loops = " << iloops << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
    }

    if(!converged) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "self_consistent_debye: Warning! convergence not achieved" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // End of self_consistent_debye function
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::debye_polynom_fit
// ***************************************************************************
namespace AGL_functions {
  //
  // debye_polynom_fit: fits (log_e (theta_debye), log_e (V)) data by polynomials.
  //     It then takes a weighted average of these polynomials.
  //
  //     It returns the averaged polynomial coefficients in polynomialcoeffs.
  //
  //
  // Adapted from original Fortran version written by V. Luana and M. A. Blanco, Departamento de Quimica Fisica y Analitica, Universidad de Oviedo, 33006-Oviedo, Spain.  
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint debye_polynom_fit (vector<double>& polynomialcoeffs, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<double> polycoeffwork, weight, xdata, ydata;
    double rms;
    int npolycoeffwork, nfit;
    // Store parameters from each fitting procedure iteration
    vector<int> npolycoeffworkfit(AGL_data.maxfit+1), ndatafit(AGL_data.maxfit+1);
    // Service variables
    int npolycoeffworkmax, ifit;
    double weight_norm, weight_rms, weight_rms_gaussian, rmsmin;
    // Limit in the dimensions of polynomials (# of points - limit)
    int limit = 4;
    // Number of pairs of data points to delete in the fitting procedure
    int ndel = 3;
    int npolycoeffworkmin;
    int ndatamin, ndataiterfit;
    uint pferr = 0;
    ostringstream aus;
    if(LDEBUG) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Starting debye_polynom_fit" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(LDEBUG) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "AGL_data.volumeinput = ";
      for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
        aus << AGL_data.volumeinput.at(i) << "\t";
      }
      aus << endl;
      aus << _AGLSTR_MESSAGE_ + "AGL_data.tdebye = ";
      for (uint i = 0; i < AGL_data.tdebye.size(); i++) {
        aus << AGL_data.tdebye.at(i) << "\t";
      }
      aus << endl;      
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // Initialize the weights and the variables
    xdata.resize(AGL_data.volumeinput.size());
    ydata.resize(AGL_data.volumeinput.size());
    weight.resize(AGL_data.volumeinput.size());
    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      xdata.at(i) = log (AGL_data.volumeinput.at(i));
      ydata.at(i) = log (AGL_data.tdebye.at(i));
      weight.at(i) = 1.0;
    }
    if(LDEBUG) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "xdata = ";
      for (uint i = 0; i < xdata.size(); i++) {
        aus << xdata.at(i) << "\t";
      }
      aus << endl;
      aus << _AGLSTR_MESSAGE_ + "ydata = ";
      for (uint i = 0; i < ydata.size(); i++) {
        aus << ydata.at(i) << "\t";
      }
      aus << endl;      
      aus << _AGLSTR_MESSAGE_ + "weight = ";
      for (uint i = 0; i < weight.size(); i++) {
        aus << weight.at(i) << "\t";
      }
      aus << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }      
    // OBSOLETE polycoeffwork.resize(AGL_data.maxpolycoeffs+1);

    // Fitting algorithm loops over different polynomials
    nfit = 0;
    int npolynomialcoeffs = 0;
    npolycoeffworkmin = 2;
    uint last_entry_tofit = 0;
    uint first_entry_tofit = AGL_data.volumeinput.size() - 1;
    ndataiterfit = AGL_data.volumeinput.size();
    ndatamin = max(ndataiterfit - 2 * ndel, npolycoeffworkmin + 1);
    rmsmin = 1e33;
    // Eliminate pairs of outermost elements ndel times:
    while (ndataiterfit >= ndatamin) {
      // Increasing the number of parameters until the limit is reached
      npolycoeffworkmax = min (ndataiterfit - limit - 1, AGL_data.maxpolycoeffs);
      for (npolycoeffwork = npolycoeffworkmin; npolycoeffwork <= npolycoeffworkmax; npolycoeffwork++) {
        nfit = nfit + 1;
        polycoeffwork.resize(npolycoeffwork+1);
        pferr = AGL_functions::polynom_fit (last_entry_tofit, first_entry_tofit, xdata, ydata, weight, rms, npolycoeffwork, polycoeffwork, AGL_data.gaussxm_debug, FileMESSAGE);
        if(pferr != 0) {
          return 2;
        }
        // Check if number of fitted polynomials exceeds maximum number allowed
        // Discard excess fitted polynomials
        if(nfit > AGL_data.maxfit) {
          nfit = AGL_data.maxfit;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "debye_polynom_fit: maximum number of fits exceeded" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Save fit parameters
        else {
          npolynomialcoeffs = max (npolynomialcoeffs,npolycoeffwork);
          npolycoeffworkfit.at(nfit) = npolycoeffwork;
          ndatafit.at(nfit) = ndataiterfit;
          rmsmin = min(rmsmin, rms * npolycoeffwork / ndataiterfit);
        }
      }
      last_entry_tofit = last_entry_tofit + 1;
      first_entry_tofit = first_entry_tofit - 1;
      ndataiterfit = ndataiterfit - 2;
    }
    // Check that there are valid fitted polynomials to average over	
    if(nfit == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "debye_polynom_fit: no fits to average!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    // Take weighted average of fitted polynomials (repeating the fitting procedure)
    weight_norm = 0.0;
    polynomialcoeffs.resize(npolynomialcoeffs+1);
    for (int i = 0; i <= npolynomialcoeffs; i++) {
      polynomialcoeffs.at(i) = 0.0;
    }
    for ( ifit = 1; ifit <= nfit; ifit++) {
      ndataiterfit = ndatafit.at(ifit);
      npolycoeffwork = npolycoeffworkfit.at(ifit);
      polycoeffwork.resize(npolycoeffwork + 1);
      last_entry_tofit = 1 + ((AGL_data.volumeinput.size() - ndataiterfit) / 2) - 1;
      first_entry_tofit = AGL_data.volumeinput.size() - last_entry_tofit - 1;
      pferr = AGL_functions::polynom_fit (last_entry_tofit, first_entry_tofit, xdata, ydata, weight, rms, npolycoeffwork, polycoeffwork, AGL_data.gaussxm_debug, FileMESSAGE);
      if(pferr != 0) {
        return 2;
      }
      weight_rms = rms * npolycoeffwork / (rmsmin * ndataiterfit);
      weight_rms_gaussian = exp(-weight_rms * weight_rms);
      for (int i = 0; i <= npolycoeffwork; i++) {
        polynomialcoeffs.at(i) = polynomialcoeffs.at(i) + weight_rms_gaussian * polycoeffwork.at(i);
      }
      weight_norm = weight_norm + weight_rms_gaussian;
    }
    // Apply the proper weight to each of the polynomials
    for (int i = 0; i <= npolynomialcoeffs; i++) {
      polynomialcoeffs.at(i) = polynomialcoeffs.at(i) / weight_norm;
    }
    // End of  routine
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************************
//  This set of functions compute the Debye model vibrational properties
// **************************************************************************************

// ***************************************************************************
// AGL_functions::thermal_properties
// ***************************************************************************
namespace AGL_functions {
  //
  // thermal_properties: compute Debye model vibrational properties.
  //
  //     This routine obtains the vibrational properties per cell of a given
  //     crystal using the Debye model: internal energy (U), heat
  //     capacity at constant volume (Cv), Helmholtz free energy (F),
  //     and vibrational entropy (S).
  //
  //     To evaluate these properties, the following integral is required:
  //                         
  //                                      |    x^3     |
  //     Debye (y) = 3*y^(-3) * INT (0,y) | ---------- | dx
  //                                      | exp(x) - 1 |
  //
  //    where y = ThetaD / T, where ThetaD is the Debye temperature (K) and
  //    T is the absolute (thermodynamic) temperature. The integral is evaluated
  //    using a Gauss-Legendre quadrature.
  //
  // Input: 
  //       ThetaD : Debye temperature (K).     
  //            T : Absolute temperature (K).
  //       natoms : Number of atoms in the unit cell.
  // Output: 
  //  internalenergy : Vibrational internal energy, U (eV/cell).
  //              Cv : Constant V heat capacity, Cv (eV/(K*cell)).
  // helmholtzenergy : Helmholtz's free energy (eV/cell).
  //         entropy : Entropy (eV/(K*cell)).
  //           Debye : Debye integral.
  //            xabs : Maximum error in Debye integral evaluation.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  double debye_function(double& z) {
    return (pow(z, 3.0)) / (exp(z) - 1);
  }

  uint thermal_properties (double& ThetaD, double& Temperature, double& DebyeIntegral, double& DebyeIntegralerror, double& internalenergy, double& Cv, double& helmholtzenergy, double& entropy, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    double eps=1e-12;
    double cero=0.0;
    int maxnl=100;
    vector<double> xpoints(maxnl), weight(maxnl);
    double yval, debye_int_prev; 
    int i, nl;
    double sumdebyeint;
    double tol2=1e-8;
    ostringstream aus;
    // Check that temperature and Debye temperature are above zero
    DebyeIntegral = cero;
    DebyeIntegralerror = cero;
    if(LDEBUG) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "thermal_properties: ThetaD = " << ThetaD << endl;
      aus << _AGLSTR_WARNING_ + "thermal_properties: Temperature = " << Temperature << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(ThetaD < 0.0 || Temperature < 0.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "thermal_properties: Bad ThetaD or T value" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    // For zero temperature, calculate zero-point vibrational properties
    if(Temperature < tol2) {
      // OBSOLETE internalenergy = AGL_data.natoms * 9.0 * physconstkbau * ThetaD / 8.0;
      internalenergy = (9.0 / 8.0) * KBOLTZEV * AGL_data.natoms * ThetaD;
      Cv = cero;
      // OBSOLETE helmholtzenergy = AGL_data.natoms * 9.0 * physconstkbau * ThetaD / 8.0;
      helmholtzenergy = (9.0 / 8.0) * KBOLTZEV * AGL_data.natoms * ThetaD;
      entropy = cero;
      return 0;
    }
    if(ThetaD < tol2) {
      // OBSOLETE internalenergy = AGL_data.natoms * 3.0 * physconstkbau * Temperature;
      internalenergy = 3.0 * KBOLTZEV * AGL_data.natoms * Temperature;
      // OBSOLETE Cv = AGL_data.natoms * 3.0 * physconstkbau;
      Cv = 3.0 * KBOLTZEV * AGL_data.natoms;
      entropy = 1e100;
      helmholtzenergy = internalenergy - Temperature * entropy;
      return 0;
    }
    yval = ThetaD / Temperature;
    DebyeIntegral = (3.0 * pow(PI, 4.0)) / (15.0 * pow(yval, 3.0));
    if(yval <= 250) {
      // Iterate with increasing number of Legendre points to evaluate the Debye integral.
      debye_int_prev = 1e30;
      nl = 5;
      DebyeIntegralerror = aurostd::abs(DebyeIntegral - debye_int_prev);
      while ((nl <= maxnl) && (DebyeIntegralerror >= eps)) {
        AGL_functions::gauss_legendre (cero, yval, xpoints, weight, nl, AGL_data, FileMESSAGE);
        sumdebyeint = 0.0;
        for (i = 0; i < nl; i++) {
          sumdebyeint = sumdebyeint + weight.at(i) * debye_function (xpoints.at(i));
        }
        DebyeIntegral = sumdebyeint * (3.0 / pow(yval, 3.0));
        DebyeIntegralerror = aurostd::abs(DebyeIntegral - debye_int_prev);
        debye_int_prev = DebyeIntegral;
        nl = nl + 5;
      }
    }

    // Thermodynamic vibrational properties
    // internalenergy = vibrational internal energy, cv = heat capacity, entropy = vibrational entropy, helmholtzenergy = Helmholtz vibrational free energy
    // OBSOLETE internalenergy = AGL_data.natoms * 3.0 * physconstkbau * (ThetaD * 3.0/8.0 + Temperature * DebyeIntegral);
    internalenergy = 3.0 * KBOLTZEV * AGL_data.natoms * (ThetaD * 3.0/8.0 + Temperature * DebyeIntegral);
    // OBSOLETE Cv = AGL_data.natoms * 3.0 * physconstkbau * (4.0 * DebyeIntegral - 3.0 * yval/(exp(yval)-1.0));
    Cv = 3.0 * KBOLTZEV * AGL_data.natoms * (4.0 * DebyeIntegral - 3.0 * yval/(exp(yval)-1.0));
    // OBSOLETE entropy = AGL_data.natoms * 3.0 * physconstkbau * (DebyeIntegral * 4.0/3.0 - log(1.0 - exp(-yval)));
    entropy = 3.0 * KBOLTZEV * AGL_data.natoms * (DebyeIntegral * 4.0/3.0 - log(1.0 - exp(-yval)));
    helmholtzenergy = internalenergy - Temperature * entropy;
    //
    return 0;
  }
} // namespace AGL_functions


// ***************************************************************************
// AGL_functions::thermalconductD
// ***************************************************************************
namespace AGL_functions {
  // Calculate thermal conductivity at acoustic Debye temperature and as a function of temperature
  // See Glen A. Slack, Solid State Physics 34, 1-75 (1979) for original method
  // See C. Toher et al., Phys. Rev. B 90, 174107 (2014) and references therein for explanation of method used here
  // Please cite these works and other relevant works cited in Toher et al. if you use results from this subroutine
  uint thermalconductD(double& acousticthetaD, vector<double>& temperature, double& kappaD, vector<double>& kappaT, double& gammaD, double& voleqD, double& avmass) {
    // OBSOLETE double kboltz = 1.3807e-23;
    // OBSOLETE double hbar = 1.05459e-34;
    double kth, kths, vDcr;
    double tol = 1e-12;

    // OBSOLETE double hbar = PLANKSCONSTANT_h / (2 * pi);
    double hbar = PLANCKSCONSTANT_hbar;

    // OBSOLETE kth = kboltz * acousticthetaD / hbar;
    kth = KBOLTZ * acousticthetaD / hbar;
    kths = kth * kth;

    vDcr = pow(voleqD, third);

    // ONSOLETE kappaD = ((0.849 * 3 * pow(4, third)) / (20 * pow(PI, 3) * (1.0 - (0.514 / gammaD) + (0.228 / (gammaD * gammaD)) ))) * kths * kboltz * vDcr * avmass / (hbar * gammaD * gammaD);
    kappaD = ((0.849 * 3 * pow(4, third)) / (20 * pow(PI, 3) * (1.0 - (0.514 / gammaD) + (0.228 / (gammaD * gammaD)) ))) * kths * KBOLTZ * vDcr * avmass / (hbar * gammaD * gammaD);
    kappaT.resize(temperature.size());
    for(uint i = 0; i < temperature.size(); i++) {
      if(aurostd::abs(temperature.at(i)) < tol) {
        kappaT.at(i) = 0.0;
      } else {
        kappaT.at(i) =  kappaD * acousticthetaD / temperature.at(i);
      }
    }
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::thermalconductDvol
// ***************************************************************************
namespace AGL_functions {
  // Calculate thermal conductivity using temperature dependent volume
  uint thermalconductDvol(double& acousticthetaD, vector<double>& temperature, vector<double>& kappaTV, double& gammaD, vector< vector<double> >& volumeq, double& avmass) {
    // OBSOLETE double kboltz = 1.3807e-23;
    // OBSOLETE double hbar = 1.05459e-34;
    double kth, kths, vDcr;
    double tol = 1e-12;

    // OBSOLETE kth = kboltz * acousticthetaD / hbar;

    // OBSOLETE double hbar = PLANKSCONSTANT_h / (2 * pi);
    double hbar = PLANCKSCONSTANT_hbar;  //ME20181020

    kth = KBOLTZ * acousticthetaD / hbar;
    kths = kth * kth;
    double gammafunc = ((0.849 * 3 * pow(4, third)) / (20 * pow(PI, 3) * (1.0 - (0.514 / gammaD) + (0.228 / (gammaD * gammaD)) )));
    kappaTV.resize(temperature.size());
    for(uint i = 0; i < temperature.size(); i++) {
      if(aurostd::abs(temperature.at(i)) < tol) {
        kappaTV.at(i) = 0.0;
      } else {
        // OBSOLETE vDcr = (pow(volumeq.at(i).at(0), third) / angstrom2bohr) * 1e-10;
        vDcr = (pow(volumeq.at(i).at(0), third)) * 1e-10;
        // OBSOLETE kappaTV.at(i) =  gammafunc * kths * kboltz * vDcr * avmass * acousticthetaD / (hbar * gammaD * gammaD *  temperature.at(i));
        kappaTV.at(i) =  gammafunc * kths * KBOLTZ * vDcr * avmass * acousticthetaD / (hbar * gammaD * gammaD *  temperature.at(i));
      }
    }
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL thermal conductivity
// **************************************************************************

#endif // _AFLOW_AGL_THERMAL_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
