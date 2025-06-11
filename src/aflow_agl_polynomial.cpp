// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_POLYNOMIAL_CPP
#define _AFLOW_AGL_POLYNOMIAL_CPP
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

// ************************************************************************************************
//  This set of functions fits (E, V) by a polynomial and finds the minimum
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::polynomial_fit
// ***************************************************************************
namespace AGL_functions {
  //
  // polynomial_fit: fits polynomials to (variable, function) data and averages
  //     them, weighted by its chi-square test probabilities. It returns
  //     the averaged polynomial coefficients in polynomialcoeffs, and the
  //     coefficients of its square in polynomialerrors.
  //
  // Adapted from original Fortran version written by V. Luana and M. A. Blanco, Departamento de Quimica Fisica y Analitica, Universidad de Oviedo, 33006-Oviedo, Spain.  
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint polynomial_fit_weight_ave (uint& imin, vector<double>& polynomialcoeffs, vector<double>& polynomialerrors, vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    // Data variables
    vector<double> polycoeffwork, weight;
    vector<vector <double> > allfits;
    int npolycoeffwork, nfit;
    // Store parameters from each fitting procedure iteration
    vector<int> npolycoeffworkfit(AGL_data.maxfit), ndatafit(AGL_data.maxfit);
    // Service variables:
    int npolycoeffworkmax, ifit;
    double weight_norm, weight_rms, weight_rms_gaussian, rmsmin, deriv_prod, plnvp, plnvm;
    int ndatapoints = ydata_to_fit.size();
    // Upper limit on the dimensions of the fitted polynomials (# of points - nlimit)
    // Restrict to be less than the number of available E(V) points
    // This helps to prevent overfitting
    // OBSOLETE int nlimit = 4;
    int nlimit;
    nlimit = min (ndatapoints / 2, 4);
    // Number of pairs of data points to delete in the fitting procedure
    // OBSOLETE int ndel = 3;
    int ndel;
    ndel = min (ndatapoints / 2, 3);
    ostringstream aus;
    // Polynomial fit error warning
    uint pferr = 0;
    uint aglerror = 0;

    // Check that there are enough data points available for fitting algorithm to actually be able to obtain a minimum
    int max_poly_coeffs = ndatapoints - nlimit - 1;
    if(max_poly_coeffs < 2) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "polynomial_fit_weight_ave: insufficient points to find minimum" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    // Initialize the weights
    weight.resize(ydata_to_fit.size());
    for (uint i = 0; i < weight.size(); i++) {
      weight.at(i) = 1.0;
    }

    // OBSOLETE polycoeffwork.resize(AGL_data.maxpolycoeffs+1);
    // OBSOLETE for (int i = 0; i <= AGL_data.maxpolycoeffs; i++) {
    // OBSOLETE   polycoeffwork.at(i) = 0.0;
    // OBSOLETE }

    // Fitting algorithm first loops over different order polynomials
    nfit = 0;
    int npolynomialcoeffs = 0;
    int npolycoeffworkmin = 2;
    uint first_entry_tofit = 0;
    uint last_entry_tofit = ydata_to_fit.size() - 1;
    int ndataiterfit = ydata_to_fit.size();
    int ndatamin = max (ndataiterfit-2*ndel, npolycoeffworkmin+1);  
    rmsmin = 1e33;
    // Eliminate pairs of outermost elements ndel times
    while ((ndataiterfit >= ndatamin) && ((first_entry_tofit < imin) && (last_entry_tofit > imin))) {
      // Increase the number of parameters until limit is reached:
      npolycoeffworkmax = min (ndataiterfit-nlimit-1, AGL_data.maxpolycoeffs);
      for (npolycoeffwork = npolycoeffworkmin; npolycoeffwork <= npolycoeffworkmax; npolycoeffwork++) {
        nfit = nfit + 1;
        polycoeffwork.resize(npolycoeffwork+1);
        pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, xdata_to_fit, ydata_to_fit, weight, AGL_data.rms, npolycoeffwork, polycoeffwork, AGL_data.gaussxm_debug, FileMESSAGE);
        if(pferr != 0) {
          return 2;
        }
        // Discard polynomial fits that don't give a minimum inside the bracket around the input data minimum
        aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(imin-1), polycoeffwork, plnvm, 1);
        aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(imin+1), polycoeffwork, plnvp, 1);
        deriv_prod = plnvm * plnvp;
        if(deriv_prod >= 0.0) {
          nfit = nfit - 1;
        }
        // Check that the number of fitted polynomials does not exceed the maximum number requested
        // Discard fitted polynomials which exceed the maximum number allowed
        else if(nfit > (AGL_data.maxfit - 1)) {
          nfit = AGL_data.maxfit - 1;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "polynomial_fit_weight_ave: Warning! maximum number of fits exceeded" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Save parameters obtained from fitting procedure
        else {
          npolynomialcoeffs = max (npolynomialcoeffs, npolycoeffwork);
          npolycoeffworkfit.at(nfit) = npolycoeffwork;
          ndatafit.at(nfit) = ndataiterfit;
          rmsmin = min(rmsmin, AGL_data.rms * npolycoeffwork / ndataiterfit);
        }
      }
      first_entry_tofit = first_entry_tofit + 1;
      last_entry_tofit = last_entry_tofit - 1;
      ndataiterfit = ndataiterfit - 2;
    }
    // Check that valid fitted polynomials have been obtained
    if(nfit == 0) {   
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "polynomial_fit_weight_ave: no fits to average!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    // Take a weighted average of the polynomial coefficients (repeating the fitting procedure)
    allfits.resize(nfit);
    weight_norm = 0.0;
    polynomialcoeffs.resize(npolynomialcoeffs + 1);
    polynomialerrors.resize(2 * npolynomialcoeffs + 1);
    for (int i = 0; i <= npolynomialcoeffs; i++) {
      polynomialcoeffs.at(i) = 0.0;
    }
    for (int i = 0; i <= 2 * npolynomialcoeffs; i++) {
      polynomialerrors.at(i) = 0.0;
    }
    for (ifit = 1; ifit <= nfit; ifit++) {
      allfits.at(ifit - 1).resize(ydata_to_fit.size());
      ndataiterfit = ndatafit.at(ifit);
      npolycoeffwork = npolycoeffworkfit.at(ifit);
      polycoeffwork.resize(npolycoeffwork + 1);
      first_entry_tofit = 1 + ((ydata_to_fit.size() - ndataiterfit) / 2) - 1;
      last_entry_tofit = ydata_to_fit.size() - first_entry_tofit - 1;
      pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, xdata_to_fit, ydata_to_fit, weight, AGL_data.rms, npolycoeffwork, polycoeffwork, AGL_data.gaussxm_debug, FileMESSAGE);
      if(pferr != 0) {
        return 2;
      }
      // OBSOLETE weight_rms = AGL_data.rms * (npolycoeffwork + 1) / (rmsmin * (ndataiterfit+1));
      weight_rms = (AGL_data.rms * npolycoeffwork) / (rmsmin * ndataiterfit);
      weight_rms_gaussian = exp(-weight_rms * weight_rms);
      for (int i = 0; i <= npolycoeffwork; i++) {
        polynomialcoeffs.at(i) = polynomialcoeffs.at(i) + weight_rms_gaussian * polycoeffwork.at(i);
        polynomialerrors.at(2*i) = polynomialerrors.at(2*i) + weight_rms_gaussian * polycoeffwork.at(i) * polycoeffwork.at(i);
        for (int j = 0; j < i; j++) {
          polynomialerrors.at(i+j) = polynomialerrors.at(i+j) + 2.0 * weight_rms_gaussian * polycoeffwork.at(j) * polycoeffwork.at(i);
        }
      }
      weight_norm = weight_norm + weight_rms_gaussian;
    }
    // Apply the proper weight to each of the fitted polynomials
    for (int i = 0; i <= npolynomialcoeffs; i++) {
      polynomialcoeffs.at(i) = polynomialcoeffs.at(i) / weight_norm;
      polynomialerrors.at(i) = polynomialerrors.at(i) / weight_norm;
    }
    for (int i = npolynomialcoeffs+1; i <= 2*npolynomialcoeffs; i++) {
      polynomialerrors.at(i) = polynomialerrors.at(i) / weight_norm;
    }
    // End of the polynomial fitting routine
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::polynomial_fit_weight_ave_debug
// ***************************************************************************
namespace AGL_functions {
  //
  // polynomial_fit_weight_ave_debug: fits polynomials to (var,function) data and averages
  //     them, weighted by its chi-square test probabilities. It returns
  //     the averaged polynomial coefficients in polynomialcoeffs, and the
  //     coefficients of its square in polynomialerrors.
  //     This routine differs from polynomial_fit_weight_ave in that it writes additional data 
  //     to files so that the user can investigate the behavior of the fitting algorithm in detail.
  //     This is useful for debugging and for understanding problematic systems.
  //     This routine is used when aflow is run with the "--debug" option.
  //
  // Adapted from original Fortran version written by V. Luana and M. A. Blanco, Departamento de Quimica Fisica y Analitica, Universidad de Oviedo, 33006-Oviedo, Spain.  
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint polynomial_fit_weight_ave_debug (uint& imin, vector<double>& polynomialcoeffs, vector<double>& polynomialerrors, vector<double>& xdata_to_fit, vector<double>& ydata_to_fit, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    // Data variables
    vector<double> polycoeffwork, weight;
    vector<vector <double> > allfits;
    int npolycoeffwork, nfit;
    // Save data parameters for each fitting procedure
    vector<int> npolycoeffworkfit(AGL_data.maxfit), ndatafit(AGL_data.maxfit);
    // Service variables
    int npolycoeffworkmax, ifit;
    double weight_norm, weight_rms, weight_rms_gaussian, rmsmin, deriv_prod, plnvp, plnvm, plnvf;
    int ndatapoints = ydata_to_fit.size();
    // Upper limit the dimensions of the fitted polynomials (# of points - limit)
    // Restrict to be less than the number of available E(V) points
    // This helps to prevent overfitting
    // OBSOLETE int nlimit = 4;
    int nlimit;
    nlimit = min (ndatapoints / 2, 4);
    // Number of pairs of data points to delete in the fitting procedure
    // OBSOLETE int ndel = 3;
    int ndel;
    ndel = min (ndatapoints / 2, 3);
    ostringstream aus;
    // Polynomial fit error warning
    uint pferr = 0;
    uint aglerror = 0;

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ + "Starting polynomial_fit_weight_ave_debug" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // Check that there are enough data points available for fitting to actually be able to obtain a minimum
    int max_poly_coeffs = ndatapoints - nlimit - 1;
    if(max_poly_coeffs < 2) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "polynomial_fit_weight_ave: insufficient points to find minimum" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ + "max_poly_coeffs = " << max_poly_coeffs << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    string ofilefittname = AGL_data.dirpathname + "/AGL_Plot_Fitting_coefficients.dat";
    stringstream ofilefittss;
    string ofileplotxname = AGL_data.dirpathname + "/AGL_Polynomial_plot_x.dat";
    stringstream ofileplotxss;
    string ofileplotvname = AGL_data.dirpathname + "/AGL_Polynomial_plot_v.dat";
    stringstream ofileplotvss;
    string ofileplotfitname = AGL_data.dirpathname + "/AGL_Polynomial_plot_fit.dat";
    stringstream ofileplotfitss;
    string ofileplotafitname = AGL_data.dirpathname + "/AGL_Polynomial_plot_all_fits.dat";
    stringstream ofileplotafitss;
    string ofilermsname = AGL_data.dirpathname + "/AGL_Polynomial_plot_rms.dat";
    stringstream ofilermsss;

    for (uint i = 0; i < xdata_to_fit.size(); i++) {
      ofileplotxss << xdata_to_fit.at(i) << "\t" << ydata_to_fit.at(i) << endl;
    }

    for (uint i = 0; i < AGL_data.volumeinput.size(); i++) {
      ofileplotvss << AGL_data.volumeinput.at(i) << "\t" << AGL_data.energyinput.at(i) << endl;
    }

    // Initialize the weights
    weight.resize(ydata_to_fit.size());
    for (uint i = 0; i < weight.size(); i++) {
      weight.at(i) = 1.0;
    }

    // OBSOLETE polycoeffwork.resize(AGL_data.maxpolycoeffs+1);
    // OBSOLETE for (int i = 0; i <= AGL_data.maxpolycoeffs; i++) {
    // OBSOLETE   polycoeffwork.at(i) = 0.0;
    // OBSOLETE }

    ofilefittss << "# Raw fitts" << endl;
    ofilefittss << endl;

    // Fitting algorithm first loops over different order polynomials
    nfit = 0;
    int npolynomialcoeffs = 0;
    int npolycoeffworkmin = 2;
    uint first_entry_tofit = 0;
    uint last_entry_tofit = ydata_to_fit.size() - 1;
    int ndataiterfit = ydata_to_fit.size();
    int ndatamin = max (ndataiterfit-2*ndel, npolycoeffworkmin+1); 
    rmsmin = 1e33;
    // Eliminate outermost pairs of elements ndel times:
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: ndataiterfit = " << ndataiterfit << endl;
    aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: ndatamin = " << ndatamin << endl;
    aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: imin = " << imin << endl;
    aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: first_entry_tofit = " << first_entry_tofit << endl;
    aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: last_entry_tofit = " << last_entry_tofit << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    while ((ndataiterfit >= ndatamin) && ((first_entry_tofit < imin) && (last_entry_tofit > imin))) {
      // Increase the number of parameters until the limit is reached
      npolycoeffworkmax = min (ndataiterfit-nlimit-1, AGL_data.maxpolycoeffs);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: npolycoeffworkmax = " << npolycoeffworkmax << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (npolycoeffwork = npolycoeffworkmin; npolycoeffwork <= npolycoeffworkmax; npolycoeffwork++) {
        nfit = nfit + 1;
        ofilefittss << "npolycoeffwork = " << npolycoeffwork << ", nfit = " << nfit << endl;
        polycoeffwork.resize(npolycoeffwork+1);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: AGL_data.rms = " << AGL_data.rms << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: npolycoeffwork = " << npolycoeffwork << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: nfit = " << nfit << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: first_entry_tofit = " << first_entry_tofit << ", last_entry_tofit = " << last_entry_tofit << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: xdata_to_fit = ";
        for (uint ij = 0; ij < xdata_to_fit.size(); ij++) {
          aus << xdata_to_fit.at(ij) << "\t";
        }
        aus << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: ydata_to_fit = ";
        for (uint ij = 0; ij < ydata_to_fit.size(); ij++) {
          aus << ydata_to_fit.at(ij) << "\t";
        }
        aus << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: weight = ";
        for (uint ij = 0; ij < weight.size(); ij++) {
          aus << weight.at(ij) << "\t";
        }
        aus << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, xdata_to_fit, ydata_to_fit, weight, AGL_data.rms, npolycoeffwork, polycoeffwork, AGL_data.gaussxm_debug, FileMESSAGE);
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: AGL_data.rms = " << AGL_data.rms << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: npolycoeffwork = " << npolycoeffwork << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: polycoeffwork = ";
        for (uint ij = 0; ij < polycoeffwork.size(); ij++) {
          aus << polycoeffwork.at(ij) << "\t";
        }
        aus << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: nfit = " << nfit << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: first_entry_tofit = " << first_entry_tofit << ", last_entry_tofit = " << last_entry_tofit << endl;
        aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: pferr = " << pferr << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        for (int ij = 0; ij < npolycoeffwork; ij++) {
          ofilefittss << "polycoeffwork[" << ij << "] = " << polycoeffwork.at(ij) << endl;
        }
        if(pferr != 0) {
          return 2;
        }
        // Discard polynomials fits that don't give a minimum inside the bracket around the input data minimum
        aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(imin-1), polycoeffwork, plnvm, 1);
        aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(imin+1), polycoeffwork, plnvp, 1);
        deriv_prod = plnvm * plnvp;
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Product of derivatives on either side of expected minimum = " << deriv_prod << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(deriv_prod >= 0.0) {
          nfit = nfit - 1;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "polynomial_fit_weight_ave_debug: no minimum in input bracket: removing fit" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Check that the number of fitted polynomials does not exceed the maximum number requested
        // Discard fitted polynomials which exceed the maximum number allowed
        else if(nfit > AGL_data.maxfit) {
          nfit = AGL_data.maxfit;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "polynomial_fit_weight_ave_debug: Warning! maximum number of fits exceeded" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // Save fitted parameters
        else {
          npolynomialcoeffs = max (npolynomialcoeffs,npolycoeffwork);
          npolycoeffworkfit.at(nfit) = npolycoeffwork;
          ndatafit.at(nfit) = ndataiterfit;
          rmsmin = min(rmsmin,AGL_data.rms*npolycoeffwork/ndataiterfit);
        }
      }
      first_entry_tofit = first_entry_tofit + 1;
      last_entry_tofit = last_entry_tofit - 1;
      ndataiterfit = ndataiterfit - 2;
    }
    // Check that valid fitted polynomials have been obtained
    if(nfit == 0) {   
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "polynomial_fit_weight_ave_debug: no fits to average!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    ofilefittss << endl;
    ofilefittss << "# Repeated fits" << endl;

    // Take a weighted average of the polynomial coefficients (repeating the fitting procedure)
    allfits.resize(nfit);
    weight_norm = 0.0;
    polynomialcoeffs.resize(npolynomialcoeffs + 1);
    polynomialerrors.resize(2*npolynomialcoeffs + 1);
    for (int i = 0; i <= npolynomialcoeffs; i++) {
      polynomialcoeffs.at(i) = 0.0;
    }
    for (int i = 0; i <= 2*npolynomialcoeffs; i++) {
      polynomialerrors.at(i) = 0.0;
    }
    for (ifit = 1; ifit <= nfit; ifit++) {
      allfits.at(ifit-1).resize(ydata_to_fit.size());
      ndataiterfit = ndatafit.at(ifit);

      npolycoeffwork = npolycoeffworkfit.at(ifit);
      first_entry_tofit = 1 + ((ydata_to_fit.size()-ndataiterfit) / 2) - 1;
      last_entry_tofit = ydata_to_fit.size() - first_entry_tofit - 1;
      ofilefittss << "npolycoeffwork = " << npolycoeffwork << ", ifit = " << ifit << endl;
      polycoeffwork.resize(npolycoeffwork+1);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: AGL_data.rms = " << AGL_data.rms << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: npolycoeffwork = " << npolycoeffwork << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: ifit = " << ifit << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: first_entry_tofit = " << first_entry_tofit << ", last_entry_tofit = " << last_entry_tofit << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: xdata_to_fit = ";
      for (uint ij = 0; ij < xdata_to_fit.size(); ij++) {
        aus << xdata_to_fit.at(ij) << "\t";
      }
      aus << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: ydata_to_fit = ";
      for (uint ij = 0; ij < ydata_to_fit.size(); ij++) {
        aus << ydata_to_fit.at(ij) << "\t";
      }
      aus << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: weight = ";
      for (uint ij = 0; ij < weight.size(); ij++) {
        aus << weight.at(ij) << "\t";
      }
      aus << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      pferr = AGL_functions::polynom_fit (first_entry_tofit, last_entry_tofit, xdata_to_fit, ydata_to_fit, weight, AGL_data.rms, npolycoeffwork, polycoeffwork, AGL_data.gaussxm_debug, FileMESSAGE);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: AGL_data.rms = " << AGL_data.rms << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: npolycoeffwork = " << npolycoeffwork << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: polycoeffwork = ";
      for (uint ij = 0; ij < polycoeffwork.size(); ij++) {
        aus << polycoeffwork.at(ij) << "\t";
      }
      aus << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: ifit = " << ifit << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: first_entry_tofit = " << first_entry_tofit << ", last_entry_tofit = " << last_entry_tofit << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (int ij = 0; ij < npolycoeffwork; ij++) {
        ofilefittss << "polycoeffwork[" << ij << "] = " << polycoeffwork.at(ij) << endl;
        ofilermsss << ifit << "\t" << AGL_data.rms << endl; 
        for (uint i = 0; i < ydata_to_fit.size(); i++) {
          aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(i), polynomialcoeffs, plnvf, 0); 
          allfits.at(ifit-1).at(i) = plnvf;
        }
      }
      if(pferr != 0) {
        return 2;
      }
      // OBSOLETE weight_rms = AGL_data.rms*(npolycoeffwork+1)/(rmsmin*(ndataiterfit+1));
      weight_rms = (AGL_data.rms * npolycoeffwork)/(rmsmin * ndataiterfit);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: weight_rms_gaussian = " << weight_rms << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: AGL_data.rms = " << AGL_data.rms << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: rmsmin = " << rmsmin << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: npolycoeffwork = " << npolycoeffwork << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      weight_rms_gaussian = exp(-weight_rms * weight_rms);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: weight_rms_gaussian = " << weight_rms_gaussian << endl;
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: weight_norm = " << weight_norm << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for (int i = 0; i <= npolycoeffwork; i++) {
        polynomialcoeffs.at(i) = polynomialcoeffs.at(i) + weight_rms_gaussian * polycoeffwork.at(i);
        polynomialerrors.at(2*i) = polynomialerrors.at(2*i) + weight_rms_gaussian * polycoeffwork.at(i) * polycoeffwork.at(i);
        for (int j = 0; j < i; j++) {
          polynomialerrors.at(i+j) = polynomialerrors.at(i+j) + 2.0 * weight_rms_gaussian * polycoeffwork.at(j) * polycoeffwork.at(i);
        }
      }
      weight_norm = weight_norm + weight_rms_gaussian;
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "polynomial_fit_weight_ave_debug: weight_norm = " << weight_norm << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    ofilefittss << endl;
    ofilefittss << "# Averaged fit" << endl;
    for (int ij = 0; ij < npolynomialcoeffs; ij++) {
      ofilefittss << "polynomialcoeffs[" << ij << "] = " << polynomialcoeffs.at(ij) << endl;
    }

    // Apply the proper weight to each of the fitted polynomials
    for (int i = 0; i <= npolynomialcoeffs; i++) {
      polynomialcoeffs.at(i) = polynomialcoeffs.at(i) / weight_norm;
      polynomialerrors.at(i) = polynomialerrors.at(i) / weight_norm;
    }
    for (int i = npolynomialcoeffs+1; i <= 2*npolynomialcoeffs; i++) {
      polynomialerrors.at(i) = polynomialerrors.at(i) / weight_norm;
    }

    ofilefittss << endl;
    ofilefittss << "# Weighted fit" << endl;
    for (int ij = 0; ij < npolynomialcoeffs; ij++) {
      ofilefittss << "polynomialcoeffs[" << ij << "] = " << polynomialcoeffs.at(ij) << endl;
    }

    for (uint i = 0; i < xdata_to_fit.size(); i++) {
      aglerror = AGL_functions::polynom_eval (xdata_to_fit.at(i), polynomialcoeffs, plnvf, 0);   
      ofileplotfitss << xdata_to_fit.at(i) << "\t" << plnvf << endl;
    }

    for (uint i = 0; i < xdata_to_fit.size(); i++) {
      ofileplotafitss << xdata_to_fit.at(i);
      for (ifit = 1; ifit <= nfit; ifit++) {
        ofileplotafitss << "\t" << allfits.at(ifit-1).at(i);
      }
      ofileplotafitss << endl;
    }

    ofilefittss << "End of fit" << endl;
    ofilefittss << endl;
    if(!aurostd::stringstream2file(ofilefittss, ofilefittname, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << ofilefittname.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    aurostd::StringstreamClean(ofilefittss);
    if(!aurostd::stringstream2file(ofileplotxss, ofileplotxname, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << ofileplotxname.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    aurostd::StringstreamClean(ofileplotxss);
    if(!aurostd::stringstream2file(ofileplotvss, ofileplotvname, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << ofileplotvname.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    aurostd::StringstreamClean(ofileplotvss);
    if(!aurostd::stringstream2file(ofileplotfitss, ofileplotfitname, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << ofileplotfitname.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    aurostd::StringstreamClean(ofileplotfitss);
    if(!aurostd::stringstream2file(ofileplotafitss, ofileplotafitname, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << ofileplotafitname.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    aurostd::StringstreamClean(ofileplotafitss);
    if(!aurostd::stringstream2file(ofilermsss, ofilermsname, "WRITE")) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Unable to open file " << ofilermsname.c_str() <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    aurostd::StringstreamClean(ofilermsss);
    //
    // End of polynomial fitting routine
    //
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::minbrack
// ***************************************************************************
namespace AGL_functions {
  //
  //  bracket_mininum: searches for the minimum value on a list, using a downhill algorithm.
  //     Starts from a provided initial index imin.
  //     The routine searches for the minimum on the list by a downhill
  //     algorithm, assuming the index as an abscissa.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  // 
  uint bracket_minimum (uint& imin, vector<double>& ydata_to_fit, ofstream& FileMESSAGE) {
    int istep;
    ostringstream aus;
    double tol = 1e-12;
    // Checks if value of imin is at or outside the limits of the provided data
    if(imin <= 0 || imin >= (ydata_to_fit.size()-1)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "bracket_minimum: trial point is outside limits" << endl;
      aus << _AGLSTR_WARNING_ + "bracket_minimum: imin = " << imin << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }
    // Decide which direction to start searching in
    istep = 0;
    if(ydata_to_fit.at(imin) < ydata_to_fit.at(imin+1) && ydata_to_fit.at(imin) < ydata_to_fit.at(imin-1)) {
      return 0;
    } 
    else if(aurostd::abs(ydata_to_fit.at(imin) - ydata_to_fit.at(imin+1)) < tol && aurostd::abs(ydata_to_fit.at(imin) - ydata_to_fit.at(imin-1)) < tol) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "bracket_minimum: Flat function, impossible to search" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    }
    else if(ydata_to_fit.at(imin) >= ydata_to_fit.at(imin+1) && ydata_to_fit.at(imin) <= ydata_to_fit.at(imin-1)) {
      istep = 1;
    }
    else if(ydata_to_fit.at(imin) <= ydata_to_fit.at(imin+1) && ydata_to_fit.at(imin) >= ydata_to_fit.at(imin-1)) {
      istep = -1;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "bracket_minimum: trial point is a maximum" << endl;
      aus << _AGLSTR_WARNING_ + "bracket_minimum: imin = " << imin << endl;
      aus << _AGLSTR_WARNING_ + "bracket_minimum: ydata_to_fit[imin - 1] = " << ydata_to_fit.at(imin-1) << endl;
      aus << _AGLSTR_WARNING_ + "bracket_minimum: ydata_to_fit[imin] = " << ydata_to_fit.at(imin) << endl;
      aus << _AGLSTR_WARNING_ + "bracket_minimum: ydata_to_fit[imin + 1] = " << ydata_to_fit.at(imin+1) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 3; 
    }
    imin = imin + istep;
    // Search for a minimum pattern
    while (imin > 0 && imin < (ydata_to_fit.size()-1)) {
      if(ydata_to_fit.at(imin) > ydata_to_fit.at(imin+istep)) {
        imin = imin + istep;
      }
      // Minimum pattern found
      else {
        return 0;
      }
    }
    // No minimum pattern
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_WARNING_ + "bracket_minimum: monotonic function, there's no minimum" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OBSOLETE ierr = 4;
    return 4;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::bracket_minimum_global
// ***************************************************************************
namespace AGL_functions {
  //
  // bracket_minimum_global: given a table of values, and an initial index imin,
  //     the routine searches for the global minimum on the list by checking
  //     all points on the list to see which one is lowest
  // 
  uint bracket_minimum_global (uint& imin, vector<double>& ydata_to_fit, ofstream& FileMESSAGE) {
    ostringstream aus;
    double funcmin;

    imin = 0;
    funcmin = ydata_to_fit.at(imin);
    for (uint istep = 1; istep < ydata_to_fit.size(); istep++) {
      if(ydata_to_fit.at(istep) < funcmin) {
        imin = istep;
        funcmin = ydata_to_fit.at(imin);
      }
    }
    // imin outside limits
    if(imin <= 0 || imin >= (ydata_to_fit.size()-1)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "minbrack: trial point outside limits" << endl;
      aus << _AGLSTR_WARNING_ + "minbrack: imin = " << imin << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    } else {
      return 0;
    }
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::polynom_bracket_minimum
// ***************************************************************************
namespace AGL_functions {
  //
  // polynom_bracket_minimum: given a polynomial and a set of x values, evaluates the polynomial at each x-point,
  //     then searches for the global minimum of the polynomial by checking all of the calculated values
  //     at all points on the list
  // 
  uint polynom_bracket_minimum (uint& imin, vector<double>& x, vector<double>& polynomialcoeffs, ofstream& FileMESSAGE) {
    ostringstream aus;
    double funcmin, plnv;
    vector<double> y(x.size());

    for(uint istep = 0; istep < x.size(); istep++) {
      // OBSOLETE AGL_functions::polin0(x.at(istep), polynomialcoeffs, plnv);
      AGL_functions::polynom_eval (x.at(istep), polynomialcoeffs, plnv, 0);
      y.at(istep) = plnv;
    }

    imin = 0;
    funcmin = y.at(imin);
    for (uint istep = 1; istep < x.size(); istep++) {
      if(y.at(istep) < funcmin) {
        imin = istep;
        funcmin = y.at(imin);
      }
    }
    // imin outside limits
    if(imin <= 0 || imin >= (x.size()-1)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "polynom_bracket_minimum: trial point outside limits" << endl;
      aus << _AGLSTR_WARNING_ + "polynom_bracket_minimum: imin = " << imin << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    } else {
      return 0;
    }
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::polynom_minimum
// ***************************************************************************
namespace AGL_functions {
  //
  // polynom_minimum: gets the minimum of a polymomial, using the Newton-Raphson
  //     method to find where its first derivative goes to zero.
  //
  //     The minimum must be in the interval [xmin, xmax], and whenever
  //     Newton's method doesn't converge, a bisection step will be given
  //     for this interval
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //

  // Checks if zero lies between the values f_a and f_b
  bool zero_bracketed_check (double& f_a, double& f_b) { 
    if((f_a >= 0.0  && f_b <= 0.0) || (f_b >= 0.0  && f_a <= 0.0)) {
      return true;
    } else {
      return false;
    }
  }

  // Checks if zero lies between the values x_a and x_b
  bool val_bracketed_check (double x_a, double x_b) { 
    if((x_a >= 0.0  && x_b <= 0.0)) {
      return true;
    } else {
      return false;
    }
  }


  uint polynom_minimum (double& xinitial, double& xmin, double& xmax, vector<double>& polynomialcoeffs, double& pxmin, ofstream& FileMESSAGE) {
    double x, x_0, x_a, x_b, dx, tol;
    double f_a, f_b, f_x;
    ostringstream aus;
    uint aglerror = 0;
    x = xinitial;
    x_a = xmin;
    x_b = xmax;
    aglerror = AGL_functions::polynom_eval (x, polynomialcoeffs, f_x, 1);
    aglerror = AGL_functions::polynom_eval (x_a, polynomialcoeffs, f_a, 1);
    aglerror = AGL_functions::polynom_eval (x_b, polynomialcoeffs, f_b, 1);
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "polynom_minimum: error evaluating polynomial" << endl;
      aurostd::StringstreamClean(aus);
      return 7;
    }
    dx = x_b - x_a;
    tol = dx * 1e-6;
    uint pmerr = 0;
    // Checks that x_a and x_b correctly bracket minimum
    if(x_a > x_b || !(AGL_functions::val_bracketed_check(x - x_a, x - x_b)) || !(AGL_functions::zero_bracketed_check(f_a, f_b))) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "polynom_minimum: wrong input data" << endl;
      aus << _AGLSTR_WARNING_ + "polynom_minimum: x = " << x << ", x_a = " << x_a << ", x_b = " << x_b << endl;
      aus << _AGLSTR_WARNING_ + "polynom_minimum: f(x) = " << f_x << ", f(x_a) = " << f_a << ", f(x_b) = " << f_b << endl;
      aus << _AGLSTR_WARNING_ + "polynom_minimum: polynomialcoeffs.size() = " << polynomialcoeffs.size() << endl;
      aus << _AGLSTR_WARNING_ + "polynom_minimum: polynomialcoeffs = ";
      for (uint i = 0; i < polynomialcoeffs.size(); i++) {
        aus << polynomialcoeffs.at(i) << "\t";
      }
      aus << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aurostd::StringstreamClean(aus);
      // Checks error type and sets pmerr appropriately to tell gibbsrun function which correction steps to attempt
      // Errors are checked in reverse order of severity (i.e. most severe error is checked last)
      // This sets pmerr signal to most severe error warning so this can be corrected first (if possible)
      if(!(AGL_functions::zero_bracketed_check(f_a, f_b))) {
        if(f_a > 0.0 && f_b > 0.0) {
          aus << _AGLSTR_WARNING_ + "polynom_minimum: f(x_a) and f(x_b) are both greater than 0" << endl;
          if(f_a > f_b) {
            pmerr = 1;
          } else {
            pmerr = 2;
          }
        } else if(f_a < 0.0 && f_b < 0.0) {
          aus << _AGLSTR_WARNING_ + "polynom_minimum: f(x_a) and f(x_b) are both less than 0" << endl;
          if(f_a > f_b) {
            pmerr = 2;
          } else {
            pmerr = 1;
          }
        }
      }
      if(!(AGL_functions::val_bracketed_check(x - x_a, x - x_b))) {
        if(x_a > x) {
          aus << _AGLSTR_WARNING_ + "polynom_minimum: x_a > x" << endl;
          pmerr = 3;
        } else if(x_b < x) {
          aus << _AGLSTR_WARNING_ + "polynom_minimum: x_b < x" << endl;
          pmerr = 4;
        }
      }
      if(x_a > x_b) {
        aus << _AGLSTR_WARNING_ + "polynom_minimum: x_a > x_b" << endl;
        pmerr = 5;
      }	
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      pxmin = x;
      return pmerr;
    }
    while (aurostd::abs(x_b - x_a) > tol && aurostd::abs(dx) > tol) {
      // Reset bracket range to previous root estimate on one side
      if(zero_bracketed_check(f_a, f_x)) {
        x_b = x;
        f_b = f_x;
      } else {
        x_a = x;
        f_a = f_x;
      }
      aglerror = AGL_functions::polynom_eval (x, polynomialcoeffs, dx, 2);
      aglerror = AGL_functions::polynom_eval (x, polynomialcoeffs, f_x, 1);
      if(aglerror != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "polynom_minimum: error evaluating polynomial" << endl;
        aurostd::StringstreamClean(aus);
        return 7;
      }
      // Estimate axis-intersection from slope of function
      x_0 = x;
      if( !(AGL_functions::val_bracketed_check(dx * (x_0 - x_a) - f_x, dx * (x_0 - x_b) - f_x))) {
        x = 0.5 * (x_b + x_a);
        dx = x - x_0;
      } else {
        dx = - f_x / dx;
        x = x_0 + dx;
      }
    }
    pxmin = x;
    return pmerr;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::autocorrect_polynom_minimum
// ***************************************************************************
namespace AGL_functions {
  //
  // autocorrect_polynom_minimum - wrapper function for polynom_minimum which checks if polynom_minimum has run correctly.
  //     If polynom_minimum returns an error, this function resets the boundaries of the search to try to fix
  //     the problem in polynom_minimum.
  //     It also checks that polynom_minimum has returned a minimum and not a maximum.
  //
  uint autocorrect_polynom_minimum (uint& itry, vector<double>& xvalues, vector<double>& polynomialcoeffs, double& xmin, ofstream& FileMESSAGE) {
    uint upperbound, lowerbound;
    ostringstream aus;
    double plnv;
    uint ierr;
    uint pmerr;
    uint aglerror;
    // Need to explicitly create upper and lower bounds for minimum search as using "min" and "max" functions could be problematic with unsigned integers
    if(itry < 2) {
      lowerbound = 0;
    } else {
      lowerbound = itry - 2;
    }
    if(itry + 2 > xvalues.size() - 1) {
      upperbound = xvalues.size() - 1;
    } else {
      upperbound = itry + 2;
    }
    pmerr = AGL_functions::polynom_minimum (xvalues.at(itry), xvalues.at(lowerbound), xvalues.at(upperbound), polynomialcoeffs, xmin, FileMESSAGE);
    // Check that polynom_minimum has actually found a minimum and not a maximum
    aglerror = AGL_functions::polynom_eval (xmin, polynomialcoeffs, plnv, 2);
    if((aglerror != 0) || (pmerr == 7)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "autocorrect_polynom_minimum: error evaluating polynomial" << endl;
      aurostd::StringstreamClean(aus);
      return 7;
    }
    if(plnv < 0.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "autocorrect_polynom_minimum: Minimum polynomial fitted to (E, V) data is actually a maximum" << endl;
      aus << _AGLSTR_WARNING_ + "autocorrect_polynom_minimum: Rebracketing minimum of polynomial points" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ierr = AGL_functions::polynom_bracket_minimum (itry, xvalues, polynomialcoeffs, FileMESSAGE);
      if(ierr != 0) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "autocorrect_polynom_minimum: Cannot find minimum of (E, V) data" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // OBSOLETE pmerr = 6;
        return 6;
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ << "autocorrect_polynom_minimum: Minimum of (E, V) data is at point imin = " << itry << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(itry < 2) {
        lowerbound = 0;
      } else {
        lowerbound = itry - 2;
      }
      if(itry + 2 > xvalues.size() - 1) {
        upperbound = xvalues.size() - 1;
      } else {
        upperbound = itry + 2;
      }
      pmerr = AGL_functions::polynom_minimum (xvalues.at(itry), xvalues.at(lowerbound), xvalues.at(upperbound), polynomialcoeffs, xmin, FileMESSAGE);
    }
    // If polynom_minimum returns an error, then AFLOW AGL tries shifting the trial point to try to correct the error
    // If this does not correct the error, then AFLOW AGL exits giving an error message
    if(pmerr != 0) {
      if(pmerr == 1) {
        while ((pmerr == 1) && ((itry + 2) < (xvalues.size() - 1))) {
          itry = itry + 1;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "autocorrect_polynom_minimum: Resetting itry to itry = " << itry << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if(itry < 2) {
            lowerbound = 0;
          } else {
            lowerbound = itry - 2;
          }
          if(itry + 2 > xvalues.size() - 1) {
            upperbound = xvalues.size() - 1;
          } else {
            upperbound = itry + 2;
          }
          pmerr = AGL_functions::polynom_minimum (xvalues.at(itry), xvalues.at(lowerbound), xvalues.at(upperbound), polynomialcoeffs, xmin, FileMESSAGE);
        }
        // If error indicator has changed from 1 to 2, then the bracket has shifted from one side of the minimum to the other
        // Need to expand the size of the bracket to incorporate the minimum
        if(pmerr == 2) {
          if(itry < 4) {
            lowerbound = 0;
          } else {
            lowerbound = itry - 4;
          }
          if(itry + 4 > xvalues.size() - 1) {
            upperbound = xvalues.size() - 1;
          } else {
            upperbound = itry + 4;
          }
          pmerr = AGL_functions::polynom_minimum (xvalues.at(itry), xvalues.at(lowerbound), xvalues.at(upperbound), polynomialcoeffs, xmin, FileMESSAGE);
        }
        if(pmerr != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "autocorrect_polynom_minimum: Cannot find minimum of (E, V) data" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 6;
        }
      } else if(pmerr == 2) {
        while ((pmerr == 2) && (itry >= 2)) {
          itry = itry - 1;
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ << "autocorrect_polynom_minimum: Resetting itry to itry = " << itry << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if(itry < 2) {
            lowerbound = 0;
          } else {
            lowerbound = itry - 2;
          }
          if(itry + 2 > xvalues.size() - 1) {
            upperbound = xvalues.size() - 1;
          } else {
            upperbound = itry + 2;
          }
          pmerr = AGL_functions::polynom_minimum (xvalues.at(itry), xvalues.at(lowerbound), xvalues.at(upperbound), polynomialcoeffs, xmin, FileMESSAGE);
        }
        // If error indicator has changed from 2 to 1, then the bracket has shifted from one side of the minimum to the other
        // Need to expand the size of the bracket to incorporate the minimum
        if(pmerr == 1) {
          if(itry < 2) {
            lowerbound = 0;
          } else {
            lowerbound = itry - 2;
          }
          if(itry + 2 > xvalues.size() - 1) {
            upperbound = xvalues.size() - 1;
          } else {
            upperbound = itry + 2;
          }
          pmerr = AGL_functions::polynom_minimum (xvalues.at(itry), xvalues.at(lowerbound), xvalues.at(upperbound), polynomialcoeffs, xmin, FileMESSAGE);
        }
        if(pmerr != 0) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_ERROR_ + "autocorrect_polynom_minimum: Cannot find minimum of (E, V) data" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 6;
        }
      } else {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "autocorrect_polynom_minimum: Cannot find minimum of (E, V) data" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 6;
      }
    } 
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************************
//  This set of functions implement numerical routines for integration and data fitting
// **************************************************************************************

// ***************************************************************************
// AGL_functions::gauss_legendre
// ***************************************************************************
namespace AGL_functions {
  //
  // gauss_legendre: Gauss-Legendre quadrature coefficients.
  //
  // Given the lower and upper limits of integration xlower and xupper, 
  //     and given n, this routine returns arrays xabscissa and weight of length n, 
  //     containing the abscissas and weights of the Gauss-Legendre 
  //     n-point quadrature formula.
  //-----------------------------------------------------------------------
  // High precision is a good idea for this routine
  //
  // Adapted from original Fortran version written by M. A. Blanco et al., which was based on the Numerical Recipes routine
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //    
  uint gauss_legendre(double& xlower, double& xupper, vector<double>& xabscissa, vector<double>& weight, int& n, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    double eps=3.0e-16;
    int i = 0, j = 0, mid = 0, id = 0;
    double xm = 0.0, xl = 0.0;
    double polynom1 = 0.0, polynom2 = 0.0, polynom3 = 0.0;
    double polynomderiv = 0.0, z = 0.0, z1 = 0.0;
    double zdiff = 0.0;
    int iloops = 0;
    ostringstream aus;
    // The roots are symmetric in the interval, so we only have to find half of them.
    mid = (n + 1) / 2;
    xm = 0.5 * (xupper + xlower);
    xl = 0.5 * (xupper - xlower);
    // Iterate over the desired roots.
    for ( i = 0; i < mid; i++) {
      id = i + 1;
      z = cos(PI * (id - 0.25)/(n + 0.5));
      // Initialize convergence test
      zdiff = 1.0;
      iloops = 0;
      while ((zdiff > eps) && (iloops < AGL_data.maxloops)) {
        iloops = iloops + 1;
        polynom1 = 1.0;
        polynom2 = 0.0;
        // Iterate over the recurrence relation to get the Legendre polynomial evaluated at z.
        for ( j = 1; j <= n; j++) {
          polynom3 = polynom2;
          polynom2 = polynom1;
          polynom1 = ((2.0*j-1.0)*z*polynom2-(j-1.0)*polynom3)/j;
        }
        // polynom1 is the desired Legendre polynomial
        // It's derivative, polynomderiv, is computed using a standard relation involving polynom2, the polynomial of one lower order.     
        polynomderiv = n * (z * polynom1 - polynom2)/(z * z - 1.0);
        z1 = z;
        // Newton's method.
        z = z1 - polynom1 / polynomderiv;
        // Convergence test
        zdiff = aurostd::abs(z-z1);
      }
      if((iloops >= AGL_data.maxloops) && (zdiff > eps)) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "gauss_legendre: Maximum number of convergence loops exceeded" << endl;
        aus << _AGLSTR_WARNING_ + "gauss_legendre: Number of loops for convergence = " << iloops << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // Scale the root to the desired interval.
      xabscissa.at(i) = xm - xl * z;  
      // Insert its symmetric counterpart inside the array.
      xabscissa.at(n-id) = xm + xl * z;
      // Compute the weights.
      weight.at(i) = 2.0 * xl / ((1.0 - z * z) * polynomderiv * polynomderiv);
      // Insert the weight for its symmetric counterpart.
      weight.at(n-id) = weight.at(i);
    }
    return 0;
  }
} // namespace AGL_functions



// ***************************************************************************
// AGL_functions::polynom_fit
// ***************************************************************************
namespace AGL_functions {
  //
  // polynom_fit: fit data by a polynomial of order npolycoeffwork 
  //
  //      Calculates the coefficients of a polynomial of order npolycoeffwork to fit a x-y data set
  //      Generates linear matrix system cmatrix from data set to be fitted
  //      Calls Gaussian elimination function to solve this system and obtain coefficients
  //      Calculates RMS deviation for the fit between the polynomial and the input data
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  // 
  uint polynom_fit(uint& first_entry_tofit, uint& last_entry_tofit, vector<double>& xdata, vector<double>& ydata, vector<double>& weight, double& rms, int& npolycoeffwork, vector<double>& polycoeffwork, bool gxmdebug, ofstream& FileMESSAGE) {
    vector<vector<double> > cmatrix;
    double weightnorm;
    int ij;
    uint aglerror = 0;
    double variance, plnv;
    ostringstream aus;
    // Calculate the normalization of the weights
    weightnorm = 0.0;
    for (uint k = first_entry_tofit; k <= last_entry_tofit; k++) {
      weightnorm = weightnorm + weight.at(k);
    }
    weightnorm = 1.0 / weightnorm;
    cmatrix.resize(npolycoeffwork+1);
    for (int i = 0; i < npolycoeffwork + 1; i++) {
      cmatrix.at(i).resize(npolycoeffwork + 2);
    }
    // Construct the linear system matrix cmatrix:
    for (int j = 0; j <= npolycoeffwork; j++) {
      cmatrix.at(j).at(npolycoeffwork+1) = 0.0;
      if(j > 0) {
        for (uint k = first_entry_tofit; k <= last_entry_tofit; k++) {
          cmatrix.at(j).at(npolycoeffwork+1) = cmatrix.at(j).at(npolycoeffwork+1) + weight.at(k) * ydata.at(k) * (pow(xdata.at(k), j));
        }
      } else {
        for (uint k = first_entry_tofit; k <= last_entry_tofit; k++) {
          cmatrix.at(j).at(npolycoeffwork+1) = cmatrix.at(j).at(npolycoeffwork+1) + weight.at(k) * ydata.at(k);
        }
      }
      cmatrix.at(j).at(npolycoeffwork+1) = weightnorm * cmatrix.at(j).at(npolycoeffwork+1);
      for (int i = j; i <= npolycoeffwork; i++) {
        cmatrix.at(i).at(j) = 0.0;
        ij = i + j;
        if(ij > 0) {
          for (uint k = first_entry_tofit; k <= last_entry_tofit; k++) {
            cmatrix.at(i).at(j) = cmatrix.at(i).at(j) + weight.at(k) * pow(xdata.at(k), ij);
          }
          cmatrix.at(i).at(j) = weightnorm * cmatrix.at(i).at(j);
        } else {
          cmatrix.at(i).at(j) = 1.0;
        }
        cmatrix.at(j).at(i) = cmatrix.at(i).at(j);
      }
    }
    // Solve the linear system for the best A()'s:
    aglerror = AGL_functions::gaussxm (cmatrix, npolycoeffwork, polycoeffwork, gxmdebug, FileMESSAGE);
    if(aglerror != 0) {
      return aglerror;
    }
    // Compute the Root Mean Square (RMS) deviation
    variance = 0.0;
    for (uint k = first_entry_tofit; k <= last_entry_tofit; k++) {
      aglerror = AGL_functions::polynom_eval (xdata.at(k), polycoeffwork, plnv, 0);
      variance = variance + weight.at(k) * pow((ydata.at(k) - plnv), 2);
    }
    rms = sqrt(variance * weightnorm);

    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::gaussxm
// ***************************************************************************
// 
// gaussxm: interface function for aurostd::GaussJordan
//      
//      Converts data sets of type vector to type xmatrix and passes to aurostd::GaussJordan 
//	  Also contains optional debugging to print out values of input matrices		
// 
namespace AGL_functions {
  uint gaussxm (vector<vector<double> >& cmatrix, int& n, vector<double>& xdata, bool gxmdebug, ofstream& FileMESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    xmatrix<double> cxm(n, n);
    xmatrix<double> bxm(n, 1);
    ostringstream aus;

    for (int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        cxm[i+1][j+1] = cmatrix.at(i).at(j);
        if(LDEBUG && gxmdebug) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_WARNING_ + "gaussxm: cxm[" << i+1 << "][" << j+1 << "] = " << cxm[i+1][j+1] << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
    }

    for (int i = 0; i < n; i++) {
      bxm[i+1][1] = cmatrix.at(i).at(n+1);
      if(LDEBUG && gxmdebug) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "gaussxm: bxm(i+1)(1) = " << bxm[i+1][1] << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    aurostd::GaussJordan(cxm, bxm);

    for (int i = 0; i < n; i++) {
      xdata.at(i) = bxm[i+1][1];
    }
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************************
//  This function evaluates polynomial and their derivatives
// **************************************************************************************

namespace AGL_functions {
  //
  // polynom_eval: Evaluates the polynomial y(x) or its derivatives at a specific point using Horner's evaluation of a polynomial.
  //
  //     The polynomial is given by:
  //
  //       y(x) = SUM(i=0, polynomialcoeffs.size()) polynomialcoeffs(i) * x**i = a_m * x^m + a_(m-1) * x^(m-1) + ... + a_2 * x^2 + a_1 * x + a_0
  //
  //     The order of the derivative of the polynomial is set by the value of norder, so that d^n y / d x^n is returned.
  //
  //     A value of norder = 0 causes the value of the polynomial itself to be returned. 
  //
  //     It is assumed that polynomialcoeffs.size() >= norder. If not, then a value of zero is returned.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint polynom_eval (double& xval, vector<double>& polynomialcoeffs, double& yval, uint norder) {
    // If the order of the polynomial is less than the order of the derivative, then return zero
    if((polynomialcoeffs.size() - norder) < 1) {
      yval = 0.0;
      return 0;
    }

    // If the order of the polynomial is greater than the order of the derivative, then evaluate the polynomial or the appropriate derivative
    // Variable to store coefficient of next term in polynomial
    double next_coeff;
    // Initialize value of y with value of coefficient of highest order term
    yval = polynomialcoeffs.at(polynomialcoeffs.size()-1);
    // Multiply value of y by index values of its coefficient for each derivative
    for (uint i = 1; i < (norder+1); i++) {
      yval = yval * (polynomialcoeffs.size()-i);
    }
    // Multiplies value of y by values of x for each term
    for (uint i = (polynomialcoeffs.size()-1); i > norder; i--) {
      next_coeff = polynomialcoeffs.at(i-1);
      for (uint j = 1; j < (norder+1); j++) {
        next_coeff = next_coeff * (i - j);
      }
      yval = yval * xval + next_coeff;
    }
    return 0;
  }
} // namespace AGL_functions



// ***************************************************************************
// AGL_functions::cvdebfit
// ***************************************************************************
namespace AGL_functions {
  //
  // cvdebfit: determines which value of the Debye temperature produces best fit to heat capacity data
  // 
  // Runs through all integer values of Debye temperature in range between minimum and maximum values
  // Generates heat capacity curve for each Debye temperature using Debye model (see eqn 23.26, Solid State Physics, Ashcroft and Mermin):
  //
  // c_V = 9 n k_B (T / Theta_D)^3 \int_0^{Theta_D / T} (x^4 e^x) / (e^x - 1)^2 dx
  // 
  uint cvdebfit(double& tdmin, double& tdmax, double& nkb, int& npoints, double& tdbest, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    int itmin = 0, itmax = 0, id = 0;
    double rmsmin = 0.0, tdtrial = 0.0, smsq = 0.0, rms = 0.0, yval = 0.0, cvt = 0.0, Debye_int_val = 0.0;
    uint aglerror = 0;
    tdmin = tdmin - fmod(tdmin, 1.0);
    tdmax = tdmax - fmod(tdmax, 1.0);

    itmin = tdmin;
    itmax = tdmax;
    rmsmin = 1e30;

    for (int i = itmin; i <= itmax; i++) {
      id = i;
      tdtrial = id;

      smsq = 0.0;

      for (int j = 1; j < npoints; j++) {
        yval = tdtrial / AGL_data.temperature_external.at(j);
        aglerror = AGL_functions::debyeint(yval, Debye_int_val, AGL_data, FileMESSAGE);
        cvt = 9.0 * nkb * Debye_int_val;
        smsq = smsq + pow((cvt - AGL_data.CvunitkB0pressure.at(j)), 2);
      }

      rms = sqrt(smsq);

      if(rms < rmsmin) {
        rmsmin = rms;
        tdbest = tdtrial;
      }
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::debyeint
// ***************************************************************************
namespace AGL_functions {
  // debyeint: Evaluate Debye integral
  double fdebyeint(double z) {
    return (pow(z, 4) * exp(z)) / (pow((exp(z) - 1), 2));
  }

  uint debyeint (double& yval, double& Debye_int_val, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    // Evaluate Debye integral
    double eps = 1e-12;
    double cero = 0.0;
    int maxnl = 100;
    vector<double> xpoints, weight;
    double debye_integral, debye_int_prev, diff_int_abs, sumdebyeint;
    int i, nl;
    xpoints.resize(maxnl);
    weight.resize(maxnl);
    // Error condition controls
    debye_integral = (3.0 * pow(pi, 4.0)) / (15.0 * pow(yval, 3.0));
    if(yval <= 250) {
      // Iterate increasing the number of Legendre points.
      debye_int_prev = 1e30;
      for (nl = 5; nl <= maxnl; nl = nl+5) {
        AGL_functions::gauss_legendre (cero, yval, xpoints, weight, nl, AGL_data, FileMESSAGE);
        sumdebyeint = 0.0;
        for (i = 0; i < nl; i++) {
          sumdebyeint = sumdebyeint + weight.at(i) * fdebyeint(xpoints.at(i));
        }
        debye_integral = sumdebyeint / (pow(yval, 3.0));
        // Convergence check
        diff_int_abs = aurostd::abs(debye_integral - debye_int_prev);
        if(diff_int_abs < eps) {
          break;
        }
        else {
          debye_int_prev = debye_integral;
        }
      }
    }
    Debye_int_val = debye_integral;
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL polynomial fitting and processing
// **************************************************************************

#endif // _AFLOW_AGL_POLYNOMIAL_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
