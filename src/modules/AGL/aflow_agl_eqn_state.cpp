// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xscalar.h"

#include "aflow_xhost.h"
#include "modules/AGL/aflow_agl_debye.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::stringstream;
using std::vector;

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
//  These functions calculate the thermal properties using different equations of state
// **************************************************************************************

// ***************************************************************************
// AGL_functions::numerical_eos
// ***************************************************************************
namespace AGL_functions {
  //
  // numerical_eos: numerical EOS calculation.
  //
  // numerical_eos computes the derivatives of the Helmholtz function and the
  //     static energy needed to obtain the Debye temperature, the static
  //     pressure, and succesive derivatives of the bulk modulus.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint numerical_eos(double& volumereference, vector<double>& energypolynomialcoeffs, vector<double>& helmholtzenergypolynomialcoeffs, vector<double>& xconfigurationvector, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    double xeqmin;
    double dFdx;
    double d2Fdx2;
    double d3Fdx3;
    double d4Fdx4;
    double Ex;
    double dEdx;
    double d2Edx2;
    double d3Edx3;
    double d4Edx4;
    double volumex3;
    double pressure_fit;
    double dFdx2_xd2Fdx2;
    double dEdx2_xd2Edx2;
    double d2Fdx2_xd3Fdx3;
    double d2Edx2_xd3Edx3;
    double dbdp;
    double d2bdp2;
    uint aglerror = 0;
    ostringstream aus;
    // Compute Pfit(p), and bulk modulus and pressure derivatives B(p), dB/dp (p), and d^2B/dp^2 (p)
    if (AGL_data.i_eqn_of_state >= 0) {
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "NUMERICAL EOS PRESSURE DERIVATIVES" << endl;
      AGL_data.outfiless << "==================================" << endl;
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "  P(GPa)         V(bohr^3)          V/V0         Pfit(GPa)       B(GPa)          B'      B''(GPa-1) " << endl;
      AGL_data.outfiless << " --------------------------------------------------------------------------------------------------- " << endl;
    }
    for (size_t k = 0; k < AGL_data.pressure_external.size(); k++) {
      xeqmin = pow(AGL_data.voleqmin[k] / volumereference, THIRD);
      aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, dFdx, 1);
      aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, d2Fdx2, 2);
      aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, d3Fdx3, 3);
      aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, d4Fdx4, 4);
      pressure_fit = -xeqmin * dFdx / (3.0 * AGL_data.voleqmin[k]) * eV_Ang3_to_GPa;
      dFdx2_xd2Fdx2 = 2.0 * dFdx - xeqmin * d2Fdx2;
      // Bulk modulus
      AGL_data.bulkmodulus[k] = -xeqmin / (9.0 * AGL_data.voleqmin[k]) * dFdx2_xd2Fdx2 * eV_Ang3_to_GPa;
      d2Fdx2_xd3Fdx3 = (d2Fdx2 - xeqmin * d3Fdx3) / dFdx2_xd2Fdx2;
      // First and second derivatives of bulk modulus with respect to pressure
      dbdp = THIRD * (2.0 - xeqmin * d2Fdx2_xd3Fdx3);
      d2bdp2 = -AGL_data.voleqmin[k] * (d2Fdx2_xd3Fdx3 * (1.0 - xeqmin * d2Fdx2_xd3Fdx3) - xeqmin * xeqmin * d4Fdx4 / dFdx2_xd2Fdx2) / (eV_Ang3_to_GPa * dFdx2_xd2Fdx2);
      if (k == 0) {
        AGL_data.bulkmodulus_0pressure = AGL_data.bulkmodulus[k];
        AGL_data.dbulkmodulusdpV_0pressure = dbdp;
        AGL_data.d2bulkmodulusdpV2_0pressure = d2bdp2;
      }
      if (AGL_data.i_eqn_of_state >= 0) {
        AGL_data.outfiless << setw(8) << setprecision(2) << std::fixed << AGL_data.pressure_external[k] << "        " << setw(10) << (AGL_data.voleqmin[k] * pow(angstrom2bohr, 3.0)) << "      " << setw(7)
                           << setprecision(6) << AGL_data.voleqmin[k] / AGL_data.voleqmin[0] << "        " << setw(10) << setprecision(2) << pressure_fit << "      " << setw(7) << AGL_data.bulkmodulus[k] << " "
                           << setw(11) << dbdp << "     " << setw(11) << setprecision(6) << d2bdp2 << endl;
      }
    }
    // Static calculation: get second derivative of the static energy with respect to volume
    if (statcalc) {
      if (AGL_data.i_eqn_of_state >= 0) {
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "INPUT AND FITTED VALUES OF THE LATTICE ENERGY" << endl;
        AGL_data.outfiless << "=============================================" << endl;
        AGL_data.outfiless << endl;
        AGL_data.outfiless << "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) " << endl;
        AGL_data.outfiless << " -------------------------------------------------- " << endl;
      }
      // Check vectors are the same size
      if ((AGL_data.d2EnergydVolume2_static.size() != xconfigurationvector.size()) || (AGL_data.energyinput.size() != xconfigurationvector.size()) || (AGL_data.volumeinput.size() != xconfigurationvector.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
        aus << _AGLSTR_ERROR_ + "xconfigurationvector.size() = " << xconfigurationvector.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < xconfigurationvector.size(); i++) {
        aglerror = AGL_functions::polynom_eval(xconfigurationvector[i], energypolynomialcoeffs, Ex, 0);
        aglerror = AGL_functions::polynom_eval(xconfigurationvector[i], energypolynomialcoeffs, dEdx, 1);
        aglerror = AGL_functions::polynom_eval(xconfigurationvector[i], energypolynomialcoeffs, d2Edx2, 2);
        dEdx2_xd2Edx2 = xconfigurationvector[i] * d2Edx2 - 2.0 * dEdx;
        volumex3 = 3.0 * AGL_data.volumeinput[i];
        AGL_data.d2EnergydVolume2_static[i] = dEdx2_xd2Edx2 * xconfigurationvector[i] / (volumex3 * volumex3);
        if (AGL_data.i_eqn_of_state >= 0) {
          AGL_data.outfiless << setw(12) << setprecision(6) << std::fixed << (AGL_data.volumeinput[i] * pow(angstrom2bohr, 3.0)) << "    " << setw(15) << (AGL_data.energyinput[i] / hart2ev) << " " << setw(18)
                             << (Ex / hart2ev) << endl;
        }
      }
      // Dynamic calculation: get static pressure and second derivative of the energy with respect to volume
    } else {
      // Check vectors are the same size
      if ((AGL_data.d2EnergydVolume2_dynamic.size() != AGL_data.pressure_external.size()) || (AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) ||
          (AGL_data.pressure_static.size() != AGL_data.pressure_external.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_dynamic.size() = " << AGL_data.d2EnergydVolume2_dynamic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.pressure_static.size() = " << AGL_data.pressure_static.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.voleqmin.size() = " << AGL_data.voleqmin.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t k = 0; k < AGL_data.pressure_external.size(); k++) {
        xeqmin = pow(AGL_data.voleqmin[k] / volumereference, THIRD);
        aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, dEdx, 1);
        aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, d2Edx2, 2);
        aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, d3Edx3, 3);
        aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, d4Edx4, 4);
        volumex3 = 3.0 * AGL_data.voleqmin[k];
        AGL_data.pressure_static[k] = (-dEdx * xeqmin / volumex3) * eV_Ang3_to_GPa;
        dEdx2_xd2Edx2 = xeqmin * d2Edx2 - 2.0 * dEdx;
        // Second derivative of energy with respect to volume
        AGL_data.d2EnergydVolume2_dynamic[k] = dEdx2_xd2Edx2 * xeqmin / (volumex3 * volumex3);
        d2Edx2_xd3Edx3 = d2Edx2 - xeqmin * d3Edx3;
        // Evaluation of Gruneisen parameter using Slater gamma approximation
        AGL_data.gamma_G[k] = (1.0 + xeqmin * d2Edx2_xd3Edx3 / dEdx2_xd2Edx2) / 6.0;
      }
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::numerical_eos_run_all_pT
// ***************************************************************************
namespace AGL_functions {
  //
  // numerical_eos_run_all_pT: numerical EOS calculation.
  //
  // numerical_eos_run_all_pT computes the derivatives of the Helmholtz function and the
  //     static energy needed to obtain the Debye temperature, the static
  //     pressure, and succesive derivatives of the bulk modulus.
  //
  // Similar to numerical_eos, except that it only calculates properties for a pressure value
  // This enables the AGL algorithm to skip certain pressure values, and thus avoid truncating
  // the pressure-temperature range.
  // Only used when AGL_data.run_all_pressure_temperature is set to true
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint numerical_eos_run_all_pT(
      double& volumereference, vector<double>& energypolynomialcoeffs, vector<double>& helmholtzenergypolynomialcoeffs, double& voleq, double& bulkmodpres, double& pres_static, double& d2EnergydVolume2_dynamic_pres, double& Gruneisen_pres) {
    double xeqmin;
    double dFdx;
    double d2Fdx2;
    double d3Fdx3;
    double d4Fdx4;
    double dEdx;
    double d2Edx2;
    double d3Edx3;
    double d4Edx4;
    double volumex3;
    double dFdx2_xd2Fdx2;
    double dEdx2_xd2Edx2;
    double d2Edx2_xd3Edx3;
    uint aglerror = 0;
    const ostringstream aus;
    xeqmin = pow(voleq / volumereference, THIRD);
    aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, dFdx, 1);
    aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, d2Fdx2, 2);
    aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, d3Fdx3, 3);
    aglerror = AGL_functions::polynom_eval(xeqmin, helmholtzenergypolynomialcoeffs, d4Fdx4, 4);
    dFdx2_xd2Fdx2 = 2.0 * dFdx - xeqmin * d2Fdx2;
    // Bulk modulus
    bulkmodpres = -xeqmin / (9.0 * voleq) * dFdx2_xd2Fdx2 * eV_Ang3_to_GPa;
    // Dynamic calculation: get static pressure and second derivative of the energy with respect to volume
    xeqmin = pow(voleq / volumereference, THIRD);
    aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, dEdx, 1);
    aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, d2Edx2, 2);
    aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, d3Edx3, 3);
    aglerror = AGL_functions::polynom_eval(xeqmin, energypolynomialcoeffs, d4Edx4, 4);
    volumex3 = 3.0 * voleq;
    pres_static = (-dEdx * xeqmin / volumex3) * eV_Ang3_to_GPa;
    dEdx2_xd2Edx2 = xeqmin * d2Edx2 - 2.0 * dEdx;
    // Second derivative of energy with respect to volume
    d2EnergydVolume2_dynamic_pres = dEdx2_xd2Edx2 * xeqmin / (volumex3 * volumex3);
    d2Edx2_xd3Edx3 = d2Edx2 - xeqmin * d3Edx3;
    // Evaluation of Gruneisen parameter using Slater gamma approximation
    Gruneisen_pres = (1.0 + xeqmin * d2Edx2_xd3Edx3 / dEdx2_xd2Edx2) / 6.0;
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::vinet_eos
// ***************************************************************************
namespace AGL_functions {
  //
  // vinet_eos: computes Vinet equation of state from (P,V) data.
  //
  // vinet_eos computes the equation of state from the (P,V) data.
  //	The equation of state has the following form:
  //     log_e (H) = A + B(1 - x)
  //     where H = P x^2 / (3 (1 - x) )
  //           A = log_e (B_0)
  //           B = 3/2 ((B_0)' - 1)
  //           x = (V / V_0)^(1/3)
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint vinet_eos(double& volume_0pressure, double& gibbsenergy_0pressure, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    vector<double> dbdp;
    vector<double> d2bdp2;
    vector<double> xconfigvector;
    vector<double> log_H;
    double log_B0;
    double sum_omxconfig;
    double sum_logH;
    double sum_1mx_logH;
    double sum_omxconfig2;
    double sum_logH2;
    double H;
    double omx_lH_diffsqsm;
    double sr_omx_diffsqsm_lH;
    double omx_lH_dsqsm_sr;
    double f_x;
    double d1fdx1;
    double d2fdx2;
    double dfdx_fx;
    double d2fdx2_fx;
    double omxdfdx_fx;
    double x2inv;
    double b0exp;
    double omx_lH_dsqsm_omx;
    double omx_lH_dsqsm_xpo;
    double omxconfig;
    double xstatcalc;
    const uint aglerror = 0;
    ostringstream aus;

    // Fit Log_e H vs (1-x)
    log_H.resize(AGL_data.pressure_external.size());
    xconfigvector.resize(AGL_data.pressure_external.size());

    sum_omxconfig = 0.0;
    sum_logH = 0.0;
    sum_1mx_logH = 0.0;
    sum_omxconfig2 = 0.0;
    sum_logH2 = 0.0;
    int n = 0;
    xconfigvector[0] = 1.0;
    for (size_t i = 1; i < AGL_data.pressure_external.size(); i++) {
      xconfigvector[i] = pow((AGL_data.voleqmin[i] / volume_0pressure), THIRD);
      H = AGL_data.pressure_external[i] * xconfigvector[i] * xconfigvector[i] / (3.0 * (1.0 - xconfigvector[i]));
      log_H[i] = log(H);
      omxconfig = 1.0 - xconfigvector[i];
      n = n + 1;
      sum_omxconfig = sum_omxconfig + omxconfig;
      sum_logH = sum_logH + log_H[i];
      sum_1mx_logH = sum_1mx_logH + omxconfig * log_H[i];
      sum_omxconfig2 = sum_omxconfig2 + omxconfig * omxconfig;
      sum_logH2 = sum_logH2 + log_H[i] * log_H[i];
    }
    log_B0 = (sum_logH * sum_omxconfig2 - sum_1mx_logH * sum_omxconfig) / (n * sum_omxconfig2 - sum_omxconfig * sum_omxconfig);
    omx_lH_diffsqsm = (n * sum_1mx_logH - sum_omxconfig * sum_logH) / (n * sum_omxconfig2 - sum_omxconfig * sum_omxconfig);
    sr_omx_diffsqsm_lH = sqrt((sum_omxconfig2 - sum_omxconfig * sum_omxconfig / n) * (sum_logH2 - sum_logH * sum_logH / n));
    omx_lH_dsqsm_sr = (sum_1mx_logH - sum_omxconfig * sum_logH / n) / sr_omx_diffsqsm_lH;
    log_H[0] = log_B0;
    // Obtain bulk modulus and its volume derivatives at zero pressure: B (p=0), dB/dp (p=0), d^2B/dp^2 (p=0)
    AGL_data.bulkmodulus_0pressure = exp(log_B0);
    AGL_data.dbulkmodulusdpV_0pressure = 2.0 * omx_lH_diffsqsm * THIRD + 1.0;
    AGL_data.d2bulkmodulusdpV2_0pressure = -(2.0 + omx_lH_diffsqsm * (omx_lH_diffsqsm + 6.0)) / (9.0 * AGL_data.bulkmodulus_0pressure);
    // Save static values
    if (statcalc) {
      AGL_data.gibbsenergystatcalc_0pressure = gibbsenergy_0pressure;
      AGL_data.volumestatcalc_0pressure = volume_0pressure;
      AGL_data.bulkmodulusstatcalc_0pressure = AGL_data.bulkmodulus_0pressure / eV_Ang3_to_GPa;
      AGL_data.Avinetstatcalc_0pressure = omx_lH_diffsqsm;
    }
    // Compute Pfit(p) and bulk modulus and pressure derivatives: B(p), dB/dp (p), and d^2B/dp^2 (p)
    dbdp.resize(AGL_data.pressure_external.size());
    d2bdp2.resize(AGL_data.pressure_external.size());
    AGL_data.bulkmodulus[0] = AGL_data.bulkmodulus_0pressure;
    dbdp[0] = AGL_data.dbulkmodulusdpV_0pressure;
    d2bdp2[0] = AGL_data.d2bulkmodulusdpV2_0pressure;
    AGL_data.pfit[0] = 0.0;
    // Check vectors are the same size
    if ((AGL_data.bulkmodulus.size() != AGL_data.pressure_external.size()) || (AGL_data.pfit.size() != AGL_data.pressure_external.size())) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.bulkmodulus.size() = " << AGL_data.bulkmodulus.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.pfit.size() = " << AGL_data.pfit.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 3;
    }
    for (size_t i = 1; i < AGL_data.pressure_external.size(); i++) {
      omx_lH_dsqsm_omx = omx_lH_diffsqsm * (1.0 - xconfigvector[i]);
      omx_lH_dsqsm_xpo = omx_lH_diffsqsm * xconfigvector[i] + 1.0;
      f_x = xconfigvector[i] * (1.0 - omx_lH_dsqsm_omx) - 2.0;
      d1fdx1 = omx_lH_dsqsm_xpo - omx_lH_dsqsm_omx;
      d2fdx2 = 2.0 * omx_lH_diffsqsm;
      dfdx_fx = d1fdx1 / f_x;
      d2fdx2_fx = d2fdx2 / f_x;
      omxdfdx_fx = 1.0 - xconfigvector[i] * dfdx_fx;
      x2inv = 1.0 / (xconfigvector[i] * xconfigvector[i]);
      b0exp = AGL_data.bulkmodulus_0pressure * exp(omx_lH_dsqsm_omx);
      AGL_data.bulkmodulus[i] = -b0exp * f_x * x2inv;
      dbdp[i] = THIRD * (omx_lH_dsqsm_xpo + omxdfdx_fx);
      d2bdp2[i] = xconfigvector[i] / (9.0 * AGL_data.bulkmodulus[i]) * (xconfigvector[i] * d2fdx2_fx - omx_lH_diffsqsm + dfdx_fx * omxdfdx_fx);
      AGL_data.pfit[i] = 3.0 * (1.0 - xconfigvector[i]) * x2inv * b0exp;
    }
    // Write output to file AGL.out
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "VINET EOS PRESSURE DERIVATIVES" << endl;
    AGL_data.outfiless << "==============================" << endl;
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "  1-V/V0 \t Vinet-Func \t P(GPa) \t Pfit(GPa) \t  B(GPa) \t      B' \t  B''(GPa-1)" << endl;
    AGL_data.outfiless << " ------------------------------------------------------------------------------------------------------------ " << endl;
    for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
      AGL_data.outfiless << setw(8) << setprecision(5) << std::fixed << 1.0 - xconfigvector[i] << "\t" << setw(11) << setprecision(6) << log_H[i] << "\t" << setw(7) << setprecision(2)
                         << AGL_data.pressure_external[i] << "\t" << setw(18) << setprecision(2) << AGL_data.pfit[i] << "\t" << setw(8) << setprecision(2) << AGL_data.bulkmodulus[i] << "\t" << setw(8)
                         << setprecision(4) << dbdp[i] << "\t" << setw(12) << setprecision(6) << d2bdp2[i] << endl;
    }
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "B0 = " << AGL_data.bulkmodulus_0pressure << ", B0' = " << AGL_data.dbulkmodulusdpV_0pressure << ", B0'' = " << AGL_data.d2bulkmodulusdpV2_0pressure << " reg.coef = " << omx_lH_dsqsm_sr << endl;
    AGL_data.outfiless << endl;

    // Static calculation: get static energy and its second derivative with respect to volume
    if (statcalc) {
      // Check vectors are the same size
      if (AGL_data.d2EnergydVolume2_static.size() != AGL_data.volumeinput.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < AGL_data.volumeinput.size(); i++) {
        xstatcalc = pow(AGL_data.volumeinput[i] / AGL_data.volumestatcalc_0pressure, THIRD);
        omx_lH_dsqsm_omx = AGL_data.Avinetstatcalc_0pressure * (1.0 - xstatcalc);
        f_x = xstatcalc * (1.0 - omx_lH_dsqsm_omx) - 2.0;
        b0exp = AGL_data.bulkmodulusstatcalc_0pressure * exp(omx_lH_dsqsm_omx);
        AGL_data.IntEnergStatic[i] = AGL_data.gibbsenergystatcalc_0pressure + 9.0 * AGL_data.volumestatcalc_0pressure / (AGL_data.Avinetstatcalc_0pressure * AGL_data.Avinetstatcalc_0pressure) *
                                                                                  (b0exp * (omx_lH_dsqsm_omx - 1.0) + AGL_data.bulkmodulusstatcalc_0pressure);
        AGL_data.d2EnergydVolume2_static[i] = -f_x / (xstatcalc * xstatcalc * AGL_data.volumeinput[i]) * b0exp;
      }
      // Print input and fitted values of the lattice energy
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "INPUT AND FITTED VALUES OF THE LATTICE ENERGY" << endl;
      AGL_data.outfiless << "=============================================" << endl;
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) " << endl;
      AGL_data.outfiless << " -------------------------------------------------- " << endl;
      // Check vectors are the same size
      if ((AGL_data.energyinput.size() != AGL_data.volumeinput.size()) || (AGL_data.IntEnergStatic.size() != AGL_data.volumeinput.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.IntEnergStatic.size() = " << AGL_data.IntEnergStatic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.energyinput.size() = " << AGL_data.energyinput.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < AGL_data.volumeinput.size(); i++) {
        AGL_data.outfiless << setw(12) << setprecision(6) << std::fixed << AGL_data.volumeinput[i] << "\t" << setw(15) << AGL_data.energyinput[i] << "\t" << setw(18) << AGL_data.IntEnergStatic[i] << endl;
      }
      // Dynamic calculation: get static pressure and second derivative of the energy
    } else {
      // Check vectors are the same size
      if ((AGL_data.d2EnergydVolume2_dynamic.size() != AGL_data.pressure_external.size()) || (AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) ||
          (AGL_data.pressure_static.size() != AGL_data.pressure_external.size()) || (AGL_data.gamma_G.size() != AGL_data.pressure_external.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_dynamic.size() = " << AGL_data.d2EnergydVolume2_dynamic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
        xstatcalc = pow(AGL_data.voleqmin[i] / AGL_data.volumestatcalc_0pressure, THIRD);
        omx_lH_dsqsm_omx = AGL_data.Avinetstatcalc_0pressure * (1.0 - xstatcalc);
        omx_lH_dsqsm_xpo = AGL_data.Avinetstatcalc_0pressure * xstatcalc + 1.0;
        f_x = xstatcalc * (1.0 - omx_lH_dsqsm_omx) - 2.0;
        d1fdx1 = omx_lH_dsqsm_xpo - omx_lH_dsqsm_omx;
        d2fdx2 = 2.0 * AGL_data.Avinetstatcalc_0pressure;
        dfdx_fx = d1fdx1 / f_x;
        omxdfdx_fx = 1.0 - xstatcalc * dfdx_fx;
        x2inv = 1.0 / (xstatcalc * xstatcalc);
        b0exp = AGL_data.bulkmodulusstatcalc_0pressure * exp(omx_lH_dsqsm_omx);
        AGL_data.pressure_static[i] = 3.0 * (1.0 - xstatcalc) * x2inv * b0exp * eV_Ang3_to_GPa;
        AGL_data.d2EnergydVolume2_dynamic[i] = -f_x * x2inv * x2inv / (xstatcalc * AGL_data.volumestatcalc_0pressure) * b0exp;
        AGL_data.gamma_G[i] = (omx_lH_dsqsm_xpo + omxdfdx_fx - 1.0) / 6.0;
      }
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::birch_murnaghan_eos
// ***************************************************************************
namespace AGL_functions {
  //
  // birch_murnaghan_eos - computes the Birch-Murnaghan equation of state of order birchfitorder_iG from the (P,V) data.
  //
  //     The equation of state has the following expression:
  //
  //     F = Sum (i = 0, birchfitorder_iG) a(i) * f^i
  //
  //     where : F = P / [3f(1 + 2f)^(5/2)]
  //             f = [x^(-2) - 1]/2
  //             x = [V(i) / V(1)]^(1/3)
  //
  // Input:
  //  volume_0pressure      : Cell volume (Ang^3/cell) at P = 0.
  //  gibbsenergy_0pressure : Gibbs energy (or 0k static energy) at volume_0pressure (in Hartree).
  //  birchfitorder_iG      : order of Birch-Murnaghan fit.
  //  pressure()    : Pressure values (GPa).
  //  volumeinput() : Initial values of the volume (Ang^3/cell).
  //  statcalc      : Logical variable that determines if the calculation is
  //            static or dynamic. In the first case the second derivative
  //            of the static energy (d2EnergydVolume2_static) is computed for all the input
  //            values of the volume. In the second case the second
  //            derivative of the static energy (d2EnergydVolume2_dynamic) is computed for
  //            the equilibrium volumes at the different pressures.
  //
  // Output:
  //  pressure_static()          : Static pressures in GPa (only on dynamic calculations).
  //  d2EnergydVolume2_static()  : Second derivative of IntEnergStatic(k) for each volumeinput() (eV/Ang^6).
  //  d2EnergydVolume2_dynamic() : Second derivative of IntEnergStatic(k) for each voleqmin() (eV/Ang^6).
  //  rms                        : Root mean square deviation.
  //  bulkmodulus_0pressure,dbulkmodulusdpV_0pressure,d2bulkmodulusdVp2_0pressure : Bulk modulus and their derivatives at P = 0.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint birch_murnaghan_eos(double& volume_0pressure, double& gibbsenergy_0pressure, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    vector<double> weight;
    vector<double> dbdp;
    vector<double> d2bdp2;
    vector<double> a_poly_coeff;
    vector<double> ybirch;
    vector<double> fbirch;
    double x0birch;
    double fbirch_im1;
    double fbirchk;
    double fb2_p1_sq;
    double fb2_p1;
    double fb2_p1_32;
    double fb2_p1_52;
    double fb2_p1_2;
    double d2bdf2;
    double d2pdf2;
    double pvd0_pvd1;
    double pvd0_pvd1_inv;
    double pvd1_pvd2;
    double volume_3;
    double volume_9;
    double pvd1_pvd0;
    double pvd2_pvd1;
    double fb52_pvd1_pvd0;
    double polval_d0;
    double polval_d1;
    double polval_d2;
    double polval_d3;
    const double tol = 1e-12;
    uint pferr = 0;
    uint aglerror = 0;
    ostringstream aus;

    uint npresm2 = AGL_data.pressure_external.size() - 2;
    uint izero = 0;

    if (aurostd::abs(AGL_data.pressure_external[0]) > tol) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "birch : P(0) must be 0.0" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 1;
    }
    a_poly_coeff.resize(AGL_data.birchfitorder_iG + 1);
    // Compute the Birch function F and strain variable f.
    weight.resize(AGL_data.pressure_external.size());
    weight[0] = 0.0;
    for (size_t i = 1; i < AGL_data.pressure_external.size(); i++) {
      x0birch = pow(AGL_data.voleqmin[i] / volume_0pressure, 1.0 / 3.0);
      fbirch.push_back((pow(x0birch, -2) - 1) / 2.0);
      ybirch.push_back(AGL_data.pressure_external[i] / (3 * fbirch[i - 1] * pow(1 + 2 * fbirch[i - 1], 2.5)) / eV_Ang3_to_GPa);
      weight[i] = 1.0;
    }
    // Fit to a polynomial of order iG.
    pferr = AGL_functions::polynom_fit(izero, npresm2, fbirch, ybirch, weight, AGL_data.rms, AGL_data.birchfitorder_iG, a_poly_coeff, AGL_data.gaussxm_debug, FileMESSAGE);
    if (pferr != 0) {
      return 4;
    }
    // Compute bulk modulus and volume derivatives at zero pressure: B(p=0), dB/dV (p=0), d^2B/dV^2 (p=0)
    // AGL_data.bulkmodulus_0pressure = a_poly_coeff[0]*au2GPa;
    AGL_data.bulkmodulus_0pressure = a_poly_coeff[0] * eV_Ang3_to_GPa;
    if (AGL_data.birchfitorder_iG == 0) {
      AGL_data.dbulkmodulusdpV_0pressure = 4.0;
      AGL_data.d2bulkmodulusdpV2_0pressure = -35.0 / (9.0 * AGL_data.bulkmodulus_0pressure);
    } else if (AGL_data.birchfitorder_iG == 1) {
      AGL_data.dbulkmodulusdpV_0pressure = 4.0 + 2.0 * a_poly_coeff[1] * eV_Ang3_to_GPa / (3.0 * AGL_data.bulkmodulus_0pressure);
      AGL_data.d2bulkmodulusdpV2_0pressure = (-AGL_data.dbulkmodulusdpV_0pressure * (AGL_data.dbulkmodulusdpV_0pressure - 7.0) - 143.0 / 9.0) / AGL_data.bulkmodulus_0pressure;
    } else if (AGL_data.birchfitorder_iG >= 2) {
      AGL_data.dbulkmodulusdpV_0pressure = 4.0 + 2.0 * a_poly_coeff[1] * eV_Ang3_to_GPa / (3.0 * AGL_data.bulkmodulus_0pressure);
      AGL_data.d2bulkmodulusdpV2_0pressure = (2.0 * a_poly_coeff[2] / (AGL_data.bulkmodulus_0pressure * 3.0) - AGL_data.dbulkmodulusdpV_0pressure * (AGL_data.dbulkmodulusdpV_0pressure - 7.0) - 143.0 / 9.0) /
                                             AGL_data.bulkmodulus_0pressure;
    }
    // Compute bulk modulus and derivatives with respect to pressure: B(p), dB/dp (p) and d^B/dp^2 (p) (b(), dbdp(), and d2bdp2())
    dbdp.resize(AGL_data.pressure_external.size());
    d2bdp2.resize(AGL_data.pressure_external.size());
    // Check vectors are the same size
    if ((AGL_data.bulkmodulus.size() != AGL_data.pressure_external.size()) || (AGL_data.pfit.size() != AGL_data.pressure_external.size())) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.bulkmodulus.size() = " << AGL_data.bulkmodulus.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.pfit.size() = " << AGL_data.pfit.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 3;
    }
    for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
      if (i == 0) {
        AGL_data.pfit[i] = 0.0;
        AGL_data.bulkmodulus[i] = AGL_data.bulkmodulus_0pressure;
        dbdp[i] = AGL_data.dbulkmodulusdpV_0pressure;
        d2bdp2[i] = AGL_data.d2bulkmodulusdpV2_0pressure;
      } else {
        fbirch_im1 = fbirch[i - 1];
        fb2_p1 = 1.0 + 2.0 * fbirch_im1;
        fb2_p1_sq = sqrt(fb2_p1);
        fb2_p1_32 = fb2_p1_sq * fb2_p1;
        fb2_p1_52 = fb2_p1_32 * fb2_p1;
        aglerror = AGL_functions::polynom_eval(fbirch_im1, a_poly_coeff, polval_d0, 0);
        aglerror = AGL_functions::polynom_eval(fbirch_im1, a_poly_coeff, polval_d1, 1);
        aglerror = AGL_functions::polynom_eval(fbirch_im1, a_poly_coeff, polval_d2, 2);
        polval_d3 = 0.0;
        if (AGL_data.birchfitorder_iG > 2) {
          aglerror = AGL_functions::polynom_eval(fbirch_im1, a_poly_coeff, polval_d3, 3);
        }
        // Fitted pressure and bulk modulus B(p)
        AGL_data.pfit[i] = 3.0 * fbirch_im1 * fb2_p1_52 * polval_d0 * eV_Ang3_to_GPa;
        pvd1_pvd0 = fb2_p1_32 * (fbirch_im1 * fb2_p1 * polval_d1 + (1.0 + 7.0 * fbirch_im1) * polval_d0);
        AGL_data.bulkmodulus[i] = fb2_p1 * pvd1_pvd0;
        pvd2_pvd1 = fb2_p1_52 * (fb2_p1 * polval_d1 + 2 * fbirch_im1 * polval_d1 + fbirch_im1 * fb2_p1 * polval_d2 + 7 * polval_d0 + (1 + 7 * fbirch_im1) * polval_d1);
        fb52_pvd1_pvd0 = 3 * fbirch_im1 * fb2_p1_52 * polval_d1 + (3.0 * fb2_p1_52 + 15.0 * fbirch_im1 * fb2_p1_32) * polval_d0;
        // dB/dp (p)
        dbdp[i] = (5 * pvd1_pvd0 + pvd2_pvd1) / fb52_pvd1_pvd0;
        d2bdf2 = 25 * fb2_p1_sq * (fbirch_im1 * fb2_p1 * polval_d1 + (1.0 + 7.0 * fbirch_im1) * polval_d0);
        d2bdf2 = d2bdf2 + 10.0 * fb2_p1_32 * ((2.0 + 11.0 * fbirch_im1) * polval_d1 + 7 * polval_d0 + fbirch_im1 * fb2_p1 * polval_d2);
        d2bdf2 = d2bdf2 + fb2_p1_52 * ((3.0 + 15.0 * fbirch_im1) * polval_d2 + 18 * polval_d1 + fbirch_im1 * fb2_p1 * polval_d3);
        d2pdf2 = 3 * fb2_p1_52 * polval_d1 + 15 * fbirch_im1 * fb2_p1_32 * polval_d1 + 3 * fbirch_im1 * fb2_p1_52 * polval_d2;
        d2pdf2 = d2pdf2 + (30 * fb2_p1_32 + 45 * fbirch_im1 * fb2_p1_sq) * polval_d0;
        d2pdf2 = d2pdf2 + (3 * fb2_p1_52 + 15 * fbirch_im1 * fb2_p1_32) * polval_d1;
        // d^2B/dp^2 (p)
        d2bdp2[i] = (fb52_pvd1_pvd0 * d2bdf2 - (5 * pvd1_pvd0 + pvd2_pvd1) * d2pdf2) / (pow(fb52_pvd1_pvd0, 3));
        AGL_data.bulkmodulus[i] = AGL_data.bulkmodulus[i] * eV_Ang3_to_GPa;
        d2bdp2[i] = d2bdp2[i] / eV_Ang3_to_GPa;
      }
    }
    // Write output to file AGL.out
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "BIRCH-MURNAGHAN EOS PRESSURE DERIVATIVES" << endl;
    AGL_data.outfiless << "========================================" << endl;
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "   Strain \t Birch-Func \t P(GPa) \t Pfit(GPa) \t  B(GPa) \t      B' \t B''(GPa-1) " << endl;
    AGL_data.outfiless << " ----------------------------------------------------------------------------------------------------------- " << endl;
    for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
      if (i == 0) {
        AGL_data.outfiless << setw(9) << setprecision(5) << std::fixed << 0.0 << "\t" << setw(11) << setprecision(5) << AGL_data.bulkmodulus_0pressure / eV_Ang3_to_GPa << "\t" << setw(7) << setprecision(4)
                           << AGL_data.pressure_external[i] << "\t" << setw(18) << setprecision(5) << AGL_data.pfit[i] << "\t" << setw(8) << setprecision(3) << AGL_data.bulkmodulus[i] << "\t" << setw(8)
                           << setprecision(5) << dbdp[i] << "\t" << setw(11) << setprecision(5) << d2bdp2[i] << endl;
      } else {
        AGL_data.outfiless << setw(9) << setprecision(5) << std::fixed << fbirch[i - 1] << "\t" << setw(11) << setprecision(5) << ybirch[i - 1] << "\t" << setw(7) << setprecision(4)
                           << AGL_data.pressure_external[i] << "\t" << setw(18) << setprecision(5) << AGL_data.pfit[i] << "\t" << setw(8) << setprecision(3) << AGL_data.bulkmodulus[i] << "\t" << setw(8)
                           << setprecision(5) << dbdp[i] << "\t" << setw(11) << setprecision(5) << d2bdp2[i] << endl;
      }
    }
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "B0 = " << AGL_data.bulkmodulus_0pressure << ", B0' = " << AGL_data.dbulkmodulusdpV_0pressure << ", B0'' = " << AGL_data.d2bulkmodulusdpV2_0pressure << ", reg.coef = " << AGL_data.rms << endl;
    AGL_data.outfiless << endl;

    if (statcalc) {
      // Compute the static potential energy U(V) and its second derivative d^2U/dV^2 (V) with respect to V for all the input values of the volume.
      AGL_data.volumestatcalc_0pressure = volume_0pressure;
      AGL_data.gibbsenergystatcalc_0pressure = gibbsenergy_0pressure;
      // Check vectors are the same size
      if (AGL_data.astatic.size() != a_poly_coeff.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.astatic.size() = " << AGL_data.astatic.size() << endl;
        aus << _AGLSTR_ERROR_ + "a_poly_coeff.size() = " << a_poly_coeff.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t k = 0; k < AGL_data.astatic.size(); k++) {
        AGL_data.astatic[k] = a_poly_coeff[k];
      }
      // Check vectors are the same size
      if ((AGL_data.IntEnergStatic.size() != AGL_data.volumeinput.size()) || (AGL_data.d2EnergydVolume2_static.size() != AGL_data.volumeinput.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.IntEnergStatic.size() = " << AGL_data.IntEnergStatic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t k = 0; k < AGL_data.volumeinput.size(); k++) {
        AGL_data.IntEnergStatic[k] = AGL_data.gibbsenergystatcalc_0pressure;
        AGL_data.d2EnergydVolume2_static[k] = 0.0;
        fbirchk = pow(AGL_data.volumeinput[k] / AGL_data.volumestatcalc_0pressure, 1.0 / 3.0);
        fbirchk = (pow(fbirchk, -2) - 1) / 2.0;
        fb2_p1 = (1.0 + 2.0 * fbirchk);
        aglerror = AGL_functions::polynom_eval(fbirchk, AGL_data.astatic, polval_d0, 0);
        aglerror = AGL_functions::polynom_eval(fbirchk, AGL_data.astatic, polval_d1, 1);
        volume_9 = 9.0 * AGL_data.volumestatcalc_0pressure;
        for (int j = 0; j <= AGL_data.birchfitorder_iG; j++) {
          AGL_data.IntEnergStatic[k] = AGL_data.IntEnergStatic[k] + volume_9 * AGL_data.astatic[j] / (j + 2) * pow(fbirchk, j + 2);
        }
        AGL_data.d2EnergydVolume2_static[k] = fb2_p1 * fb2_p1 * fb2_p1 * fb2_p1 / AGL_data.volumestatcalc_0pressure * (fbirchk * fb2_p1 * polval_d1 + (1.0 + 7.0 * fbirchk) * polval_d0);
      }
      // Print input and fitted values of the lattice energy to file AGL.out
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "INPUT AND FITTED VALUES OF THE LATTICE ENERGY" << endl;
      AGL_data.outfiless << "=============================================" << endl;
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) " << endl;
      AGL_data.outfiless << " -------------------------------------------------- " << endl;
      for (size_t i = 0; i < AGL_data.volumeinput.size(); i++) {
        AGL_data.outfiless << setw(12) << setprecision(6) << std::fixed << AGL_data.volumeinput[i] << "\t" << setw(15) << AGL_data.energyinput[i] << "\t" << setw(18) << AGL_data.IntEnergStatic[i] << endl;
      }
      return aglerror;
    } else {
      // Check vectors are the same size
      if ((AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) || (AGL_data.gamma_G.size() != AGL_data.pressure_external.size()) ||
          (AGL_data.d2EnergydVolume2_dynamic.size() != AGL_data.pressure_external.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.voleqmin.size() = " << AGL_data.voleqmin.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.gamma_G.size() = " << AGL_data.gamma_G.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_dynamic.size() = " << AGL_data.d2EnergydVolume2_dynamic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      // Compute the second derivative d^2U/dV^2(V) with respect to V for all the equilibrium values of the volume at the different pressures
      for (size_t k = 0; k < AGL_data.pressure_external.size(); k++) {
        fbirchk = pow(AGL_data.voleqmin[k] / AGL_data.volumestatcalc_0pressure, THIRD);
        fbirchk = (pow(fbirchk, -2) - 1) / 2.0;
        fb2_p1 = (1.0 + 2.0 * fbirchk);
        fb2_p1_2 = fb2_p1 * fb2_p1;
        aglerror = AGL_functions::polynom_eval(fbirchk, AGL_data.astatic, polval_d0, 0);
        aglerror = AGL_functions::polynom_eval(fbirchk, AGL_data.astatic, polval_d1, 1);
        aglerror = AGL_functions::polynom_eval(fbirchk, AGL_data.astatic, polval_d2, 2);
        polval_d3 = 0.0;
        if (AGL_data.birchfitorder_iG > 2) {
          aglerror = AGL_functions::polynom_eval(fbirchk, AGL_data.astatic, polval_d3, 3);
        }
        AGL_data.pressure_static[k] = 3.0 * fbirchk * pow(fb2_p1, 2.5) * polval_d0 * eV_Ang3_to_GPa;
        pvd0_pvd1 = (1.0 + 7.0 * fbirchk) * polval_d0 + fbirchk * fb2_p1 * polval_d1;
        AGL_data.d2EnergydVolume2_dynamic[k] = fb2_p1_2 * fb2_p1_2 / AGL_data.volumestatcalc_0pressure * pvd0_pvd1;
        pvd0_pvd1_inv = 1.0 / pvd0_pvd1;
        pvd1_pvd2 = fb2_p1 * pvd0_pvd1_inv * (7.0 * polval_d0 + (2.0 + 11.0 * fbirchk) * polval_d1 + fbirchk * fb2_p1 * polval_d2);
        volume_3 = AGL_data.voleqmin[k] / (3.0 * AGL_data.volumestatcalc_0pressure);
        AGL_data.gamma_G[k] = -2.0 * THIRD + 0.5 * fb2_p1 * sqrt(fb2_p1) * volume_3 * (8.0 + pvd1_pvd2);
      }
    }
    // End of Birch-Murnaghan routine
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::bcnt_eos
// ***************************************************************************
namespace AGL_functions {
  //
  // bcnt_eos: compute the Spinodal (BCNT) equation of state from (B, p) data.
  //
  //     The EOS has the following expresion:
  //
  //     B(p) = ( p - p_sp )^beta / K
  //
  //     where//c
  //     beta = 0.85  (If optimize_beta = true => beta is optimized)
  //     (-p_sp) and K are the parameters to optimize.
  //
  //     These parameters bear the following relation with B_0 and B_0'.
  //
  //     B_0  = (-p_sp)^{beta} / K
  //
  //     dB_0/dp = beta B_0 (-p_sp)^{-1}
  //
  // Input parameters:
  //     volume_0pressure      : Zero pressure volume, either static or dynamic.
  //     gibbsenergy_0pressure : Zero pressure Gibbs function.
  //     bulkmod_0pressure     : Bulk modulus used to compute the initial value of  -p_sp (GPa).
  //
  //     optimize_beta  : if true  => beta is optimized.
  //     statcalc : if true  => static calculation.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint bcnt_eos(double& volume_0pressure, double& gibbsenergy_0pressure, double& bulkmod_0pressure, bool statcalc, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    const double eps = 1e-10;
    const int maxnl = 100;
    vector<double> dbdp;
    vector<double> d2bdp2;
    vector<double> xpoints(maxnl);
    vector<double> weights(maxnl);
    double x_a;
    double x_b;
    double x_c;
    double x_Press_sp;
    double f_a;
    double f_b;
    double f_c;
    double xbm1;
    double xmlV;
    double inv1mbeta;
    double pxbm1;
    double xmstlVi;
    double inv1mbst;
    double bst_1mbst;
    double xswap;
    double xinf;
    double xsup;
    double signfactor;
    double sum_xm;
    double sum_prev_iter;
    double xmlxp;
    double xabs;
    double tol = 1e-12;
    double f_x_Press_sp;
    double Vol_sp;
    int nl;
    uint aglerror = 0;
    ostringstream aus;

    AGL_data.bcntvolume0pressure = volume_0pressure;

    // Initial values of properties to optimize.
    AGL_data.bcnt_beta = 0.85;
    if (!statcalc) {
      if (AGL_data.i_optimize_beta == 2) {
        AGL_data.bcnt_beta = AGL_data.bcnt_beta_statcalc;
      }
    }
    x_Press_sp = AGL_data.bcnt_beta * bulkmod_0pressure / 4.0;
    // Optimize beta and x_Press_sp.
    x_a = x_Press_sp * 0.5;
    x_b = x_Press_sp;
    x_c = x_Press_sp * 2.0;

    aglerror = AGL_functions::bcnt_bracket_minimum(x_a, x_b, x_c, f_a, f_b, f_c, AGL_data, FileMESSAGE);
    aglerror = AGL_functions::bcnt_brent_minimum(x_a, x_b, x_c, tol, x_Press_sp, AGL_data, f_x_Press_sp, FileMESSAGE);
    if (aglerror != 0) {
      return aglerror;
    }

    // Final properties
    xbm1 = pow(x_Press_sp, (AGL_data.bcnt_beta - 1));
    AGL_data.x_m_opt = AGL_data.x_K_opt / xbm1 / (1.0 - AGL_data.bcnt_beta);
    AGL_data.bulkmodulus_0pressure = xbm1 * x_Press_sp / AGL_data.x_K_opt;
    AGL_data.dbulkmodulusdpV_0pressure = AGL_data.bcnt_beta * AGL_data.bulkmodulus_0pressure / x_Press_sp;
    AGL_data.d2bulkmodulusdpV2_0pressure = AGL_data.bcnt_beta * (AGL_data.bcnt_beta - 1.0) * AGL_data.bulkmodulus_0pressure / (x_Press_sp * x_Press_sp);
    Vol_sp = volume_0pressure * exp(AGL_data.bcnt_beta / (1 - AGL_data.bcnt_beta) / AGL_data.dbulkmodulusdpV_0pressure);
    AGL_data.Press_sp_final = x_Press_sp;
    AGL_data.xsup_K_final = AGL_data.x_K_opt;
    AGL_data.Vol_sp_final = Vol_sp;
    AGL_data.bcnt_beta_final = AGL_data.bcnt_beta;

    // Save static values
    if (statcalc) {
      AGL_data.gibbsenergystatcalc_0pressure = gibbsenergy_0pressure;
      AGL_data.bulkmodulusstatcalc_0pressure = AGL_data.bulkmodulus_0pressure / eV_Ang3_to_GPa;
      AGL_data.volumestatcalc_0pressure = volume_0pressure;
      AGL_data.Vol_sp_statcalc = Vol_sp;
      AGL_data.x_K_opt_statcalc = AGL_data.x_K_opt;
      AGL_data.x_m_opt_statcalc = AGL_data.x_m_opt;
      AGL_data.x_Press_sp_statcalc = x_Press_sp;
      AGL_data.bcnt_beta_statcalc = AGL_data.bcnt_beta;
    }
    // Compute Pfit(p) and bulk modulus and its pressure derivatives B(p), dB/dp (p), and d^2B/dp^2 (p)
    dbdp.resize(AGL_data.pressure_external.size());
    d2bdp2.resize(AGL_data.pressure_external.size());
    // Check vectors are the same size
    if ((AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) || (AGL_data.pfit.size() != AGL_data.pressure_external.size()) || (AGL_data.bulkmodulus.size() != AGL_data.pressure_external.size())) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.bulkmodulus.size() = " << AGL_data.bulkmodulus.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.pfit.size() = " << AGL_data.pfit.size() << endl;
      aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 3;
    }
    for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
      xmlV = (AGL_data.x_m_opt + log(volume_0pressure / AGL_data.voleqmin[i])) / AGL_data.x_m_opt;
      inv1mbeta = 1.0 / (1.0 - AGL_data.bcnt_beta);
      AGL_data.pfit[i] = x_Press_sp * (exp(inv1mbeta * log(xmlV)) - 1.0);
      pxbm1 = pow((AGL_data.pressure_external[i] + x_Press_sp), (AGL_data.bcnt_beta - 1));
      AGL_data.bulkmodulus[i] = pxbm1 * (AGL_data.pressure_external[i] + x_Press_sp) / AGL_data.x_K_opt;
      dbdp[i] = AGL_data.bcnt_beta * pxbm1 / AGL_data.x_K_opt;
      d2bdp2[i] = AGL_data.bcnt_beta * (AGL_data.bcnt_beta - 1) * pxbm1 / AGL_data.x_K_opt / (AGL_data.pressure_external[i] + x_Press_sp);
    }
    // Write output to file AGL.out
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "SPINODAL EOS PRESSURE DERIVATIVES" << endl;
    AGL_data.outfiless << "=================================" << endl;
    AGL_data.outfiless << endl;
    AGL_data.outfiless << "   P(GPa) \t Pfit(GPa) \t B(GPa) \t B' \t B''(GPa-1)" << endl;
    AGL_data.outfiless << " ------------------------------------------------------------------- " << endl;
    for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
      AGL_data.outfiless << setw(9) << setprecision(2) << std::fixed << AGL_data.pressure_external[i] << "\t" << setw(10) << setprecision(2) << AGL_data.pfit[i] << "\t" << setw(7) << setprecision(2)
                         << AGL_data.bulkmodulus[i] << "\t" << setw(11) << setprecision(4) << dbdp[i] << "\t" << setw(11) << setprecision(6) << d2bdp2[i] << endl;
    }
    AGL_data.outfiless << endl;
    if (AGL_data.optimize_beta) {
      AGL_data.outfiless << "B0 = " << AGL_data.bulkmodulus_0pressure << ", B0' = " << dbdp[0] << ", B0'' = " << d2bdp2[0] << ", reg.coef = " << f_x_Press_sp << ", Vsp = " << Vol_sp << ", -Psp = " << x_Press_sp
                         << ", K* = " << AGL_data.x_K_opt << ", gamma = " << AGL_data.bcnt_beta << endl;
    } else {
      AGL_data.outfiless << "B0 = " << AGL_data.bulkmodulus[0] << ", B0' = " << dbdp[0] << ", B0'' = " << d2bdp2[0] << ", reg.coef = " << f_x_Press_sp << ", Vsp = " << Vol_sp << ", -Psp = " << x_Press_sp
                         << ", K* = " << AGL_data.x_K_opt << ", gamma = " << AGL_data.bcnt_beta << endl;
    }
    // Static calculation: get static energy and its second derivative
    if (statcalc) {
      // Check vectors are the same size
      if ((AGL_data.d2EnergydVolume2_static.size() != AGL_data.volumeinput.size()) || (AGL_data.IntEnergStatic.size() != AGL_data.volumeinput.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_static.size() = " << AGL_data.d2EnergydVolume2_static.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.IntEnergStatic.size() = " << AGL_data.IntEnergStatic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.volumeinput.size() = " << AGL_data.volumeinput.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < AGL_data.volumeinput.size(); i++) {
        xmstlVi = (AGL_data.x_m_opt_statcalc + log(AGL_data.volumestatcalc_0pressure / AGL_data.volumeinput[i])) / AGL_data.x_m_opt_statcalc;
        inv1mbst = 1.0 / (1.0 - AGL_data.bcnt_beta_statcalc);
        bst_1mbst = AGL_data.bcnt_beta_statcalc / (1.0 - AGL_data.bcnt_beta_statcalc);
        // Numerically integrate the Helmholtz function by means of a loop with increasing number of Legendre points (obtained from Gauss-Legendre quadrature)
        xinf = 1.0;
        xsup = xmstlVi;
        signfactor = 1.0;
        if (xsup < xinf) {
          xswap = xinf;
          xinf = xsup;
          xsup = xswap;
          signfactor = -1.0;
        }
        // Iterative loop with increasing number of Legendre points
        sum_prev_iter = 1e30;
        nl = 5;
        xabs = 1.0;
        while ((nl <= maxnl) && (xabs >= eps)) {
          gauss_legendre(xinf, xsup, xpoints, weights, nl, AGL_data, FileMESSAGE);
          sum_xm = 0.0;
          for (int ii = 0; ii < nl; ii++) {
            xmlxp = exp(AGL_data.x_m_opt_statcalc * (1.0 - xpoints[ii])) * (exp(inv1mbst * log(xpoints[ii])) - 1.0);
            sum_xm = sum_xm + weights[ii] * signfactor * xmlxp * AGL_data.x_m_opt_statcalc * AGL_data.volumestatcalc_0pressure * AGL_data.x_Press_sp_statcalc;
          }
          xabs = aurostd::abs(sum_xm - sum_prev_iter);
          sum_prev_iter = sum_xm;
          nl = nl + 5;
        }
        AGL_data.IntEnergStatic[i] = AGL_data.gibbsenergystatcalc_0pressure + sum_xm / eV_Ang3_to_GPa;
        AGL_data.d2EnergydVolume2_static[i] = AGL_data.x_Press_sp_statcalc / (AGL_data.x_m_opt_statcalc * AGL_data.volumeinput[i] * (1.0 - AGL_data.bcnt_beta_statcalc)) * exp(bst_1mbst * log(xmstlVi)) / eV_Ang3_to_GPa;
      }
      // Print input and fitted values of the lattice energy to file AGL.out
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "INPUT AND FITTED VALUES OF THE LATTICE ENERGY" << endl;
      AGL_data.outfiless << "=============================================" << endl;
      AGL_data.outfiless << endl;
      AGL_data.outfiless << "   V(bohr^3)     E_inp(hartree)     E_fit(hartree) " << endl;
      AGL_data.outfiless << " -------------------------------------------------- " << endl;
      for (size_t i = 0; i < AGL_data.volumeinput.size(); i++) {
        AGL_data.outfiless << setw(12) << setprecision(6) << std::fixed << (AGL_data.volumeinput[i] * pow(angstrom2bohr, 3.0)) << "\t" << setw(15) << (AGL_data.energyinput[i] / hart2ev) << "\t" << setw(18)
                           << (AGL_data.IntEnergStatic[i] / hart2ev) << endl;
      }
    }
    // Dynamic calculation: get static pressure and second derivative of the energy with respect to volume
    else {
      if ((AGL_data.voleqmin.size() != AGL_data.pressure_external.size()) || (AGL_data.gamma_G.size() != AGL_data.pressure_external.size()) ||
          (AGL_data.d2EnergydVolume2_dynamic.size() != AGL_data.pressure_external.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.voleqmin.size() = " << AGL_data.voleqmin.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.gamma_G.size() = " << AGL_data.gamma_G.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.d2EnergydVolume2_dynamic.size() = " << AGL_data.d2EnergydVolume2_dynamic.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
        xmlV = (AGL_data.x_m_opt_statcalc + log(AGL_data.volumestatcalc_0pressure / AGL_data.voleqmin[i])) / AGL_data.x_m_opt_statcalc;
        inv1mbeta = 1.0 / (1.0 - AGL_data.bcnt_beta_statcalc);
        bst_1mbst = AGL_data.bcnt_beta_statcalc / (1.0 - AGL_data.bcnt_beta_statcalc);
        AGL_data.pressure_static[i] = AGL_data.x_Press_sp_statcalc * (exp(inv1mbeta * log(xmlV)) - 1.0);
        AGL_data.d2EnergydVolume2_dynamic[i] = AGL_data.x_Press_sp_statcalc / (AGL_data.x_m_opt_statcalc * AGL_data.volumestatcalc_0pressure * (1.0 - AGL_data.bcnt_beta_statcalc)) *
                                               exp(-AGL_data.x_m_opt_statcalc * (1.0 - xmlV)) * exp(bst_1mbst * log(xmlV)) / eV_Ang3_to_GPa;
        AGL_data.gamma_G[i] = -1.0 / 6.0 + AGL_data.bcnt_beta_statcalc * inv1mbeta / 2.0 / AGL_data.x_m_opt_statcalc / xmlV;
      }
    }
    // End of function

    return aglerror;
  }
} // namespace AGL_functions

// **************************************************************************************
//  This set of functions implement routines required for the BCNT EOS
// **************************************************************************************

// ***************************************************************************
// AGL_functions::bcnt_bracket_minimum
// ***************************************************************************
namespace AGL_functions {
  //
  // bcnt_bracket_minimum - brackets a minimum of the function f.
  //
  //     Given a function, and two distinct initial points x_a and x_b,
  //     this routine searches in the downhill direction (defined by the
  //     function as evaluated at the initial points) and returns new
  //     points x_a, x_b, and x_c which bracket a minimum of the function.
  //     Also returned are the function values at the three points: f_a, f_b,
  //     and f_c.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint bcnt_bracket_minimum(double& x_a, double& x_b, double& x_c, double& f_a, double& f_b, double& f_c, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    const double golden_ratio = 1.618034;
    const double glimit = 100.0;
    const double tiny = 1e-20;
    double x_ab_swap;
    double xba_fbc;
    double xbc_fba;
    double x_u;
    bool endpart;
    double f_u;
    double ulim;
    const ofstream mbrkfile;
    uint aglerror = 0;

    aglerror = AGL_functions::bcnt_optimize_beta(x_a, AGL_data, f_a, FileMESSAGE);
    aglerror = AGL_functions::bcnt_optimize_beta(x_b, AGL_data, f_b, FileMESSAGE);
    if (aglerror != 0) {
      return aglerror;
    }
    // Ensure set-up is such that f_a > f_b
    if (f_b > f_a) {
      x_ab_swap = x_a;
      x_a = x_b;
      x_b = x_ab_swap;
      x_ab_swap = f_b;
      f_b = f_a;
      f_a = x_ab_swap;
    }
    // Choose new point using golden ratio to select point between x_a and x_b
    x_c = x_b + golden_ratio * (x_b - x_a);
    aglerror = AGL_functions::bcnt_optimize_beta(x_c, AGL_data, f_c, FileMESSAGE);
    if (aglerror != 0) {
      return aglerror;
    }
    while (f_b >= f_c) {
      endpart = true;
      xba_fbc = (x_b - x_a) * (f_b - f_c);
      xbc_fba = (x_b - x_c) * (f_b - f_a);
      x_u = x_b - ((x_b - x_c) * xbc_fba - (x_b - x_a) * xba_fbc) / (2.0 * copysign(std::max(aurostd::abs(xbc_fba - xba_fbc), tiny), xbc_fba - xba_fbc));
      ulim = x_b + glimit * (x_c - x_b);
      if ((x_b - x_u) * (x_u - x_c) > 0.0) {
        aglerror = AGL_functions::bcnt_optimize_beta(x_u, AGL_data, f_u, FileMESSAGE);
        if (aglerror != 0) {
          return aglerror;
        }
        if (f_u < f_c) {
          x_a = x_b;
          f_a = f_b;
          x_b = x_u;
          f_b = f_u;
          endpart = false;
        } else if (f_u > f_b) {
          x_c = x_u;
          f_c = f_u;
          endpart = false;
        } else {
          x_u = x_c + golden_ratio * (x_c - x_b);
          aglerror = bcnt_optimize_beta(x_u, AGL_data, f_u, FileMESSAGE);
          if (aglerror != 0) {
            return aglerror;
          }
        }
      } else if ((x_c - x_u) * (x_u - ulim) > 0.0) {
        aglerror = AGL_functions::bcnt_optimize_beta(x_u, AGL_data, f_u, FileMESSAGE);
        if (aglerror != 0) {
          return aglerror;
        }
        if (f_u < f_c) {
          x_b = x_c;
          x_c = x_u;
          x_u = x_c + golden_ratio * (x_c - x_b);
          f_b = f_c;
          f_c = f_u;
          aglerror = AGL_functions::bcnt_optimize_beta(x_u, AGL_data, f_u, FileMESSAGE);
          if (aglerror != 0) {
            return aglerror;
          }
        }
      } else if ((x_u - ulim) * (ulim - x_c) >= 0.0) {
        x_u = ulim;
        aglerror = AGL_functions::bcnt_optimize_beta(x_u, AGL_data, f_u, FileMESSAGE);
        if (aglerror != 0) {
          return aglerror;
        }
      } else {
        x_u = x_c + golden_ratio * (x_c - x_b);
        aglerror = AGL_functions::bcnt_optimize_beta(x_u, AGL_data, f_u, FileMESSAGE);
        if (aglerror != 0) {
          return aglerror;
        }
      }
      if (endpart) {
        x_a = x_b;
        x_b = x_c;
        x_c = x_u;
        f_a = f_b;
        f_b = f_c;
        f_c = f_u;
      }
    }
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::bcnt_brent_minimum
// ***************************************************************************
namespace AGL_functions {
  //
  // bcnt_brent_minimum - unidimensional minimization of f in the range [x_a, x_c].
  //
  //    Given a function f, and a bracketing triplet of abscissas this
  //    routine isolates the minimum to within a given tolerance
  //    using Brent's method. The bracketing triplet must be such that x_b
  //    is between x_a and x_c, and that f(x_b) is less than both f(x_a) and
  //    f(x_c). The abscissa of the minimum is returned as xmin, and the
  //    minimum function value as brent_x.
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint bcnt_brent_minimum(double& x_a, double& x_b, double& x_c, double& tol, double& xmin, _AGL_data& AGL_data, double& brent_x, ofstream& FileMESSAGE) {
    const int itmax = 100;
    const double goldenratio_diff2 = 0.3819660;
    const double zeps = 1.0e-10;
    double xlower;
    double xupper;
    double xwim1;
    double xwim2;
    double xwork;
    double xbound_diff;
    double xbdg;
    double f_xwork;
    double f_xwim1;
    double f_xwim2;
    double xwdg;
    double f_xwdg;
    double tolx;
    double tolx2;
    const double tole12 = 1e-12;
    double xmid;
    int iter;
    double xw1_fxw2;
    double xw2_fxw1;
    double xw12_fxw21;
    double xbd_tmp;
    bool xbdg_reset = true;
    ostringstream aus;
    uint aglerror = 0;

    xlower = std::min(x_a, x_c);
    xupper = std::max(x_a, x_c);
    xwim2 = x_b;
    xwim1 = x_b;
    xwork = x_b;
    xbound_diff = 0.0;
    xbdg = 0.0;
    aglerror = AGL_functions::bcnt_optimize_beta(xwork, AGL_data, f_xwork, FileMESSAGE);
    if (aglerror != 0) {
      return aglerror;
    }
    f_xwim2 = f_xwork;
    f_xwim1 = f_xwork;
    for (iter = 1; iter <= itmax; iter++) {
      xmid = 0.5 * (xlower + xupper);
      tolx = tol * aurostd::abs(xwork) + zeps;
      tolx2 = 2.0 * tolx;
      // Check if xwork is within tolerance range
      if (aurostd::abs(xwork - xmid) <= (tolx2 - 0.5 * (xupper - xlower))) {
        xmin = xwork;
        brent_x = f_xwork;
        return aglerror;
      } else {
        if (aurostd::abs(xbound_diff) > tolx) {
          xw1_fxw2 = (xwork - xwim1) * (f_xwork - f_xwim2);
          xw2_fxw1 = (xwork - xwim2) * (f_xwork - f_xwim1);
          xw12_fxw21 = (xwork - xwim2) * xw2_fxw1 - (xwork - xwim1) * xw1_fxw2;
          xw2_fxw1 = 2.0 * (xw2_fxw1 - xw1_fxw2);
          if (xw2_fxw1 > 0.0) {
            xw12_fxw21 = -xw12_fxw21;
          }
          xw2_fxw1 = aurostd::abs(xw2_fxw1);
          xbd_tmp = xbound_diff;
          xbound_diff = xbdg;
          if (aurostd::abs(xw12_fxw21) >= aurostd::abs(0.5 * xw2_fxw1 * xbd_tmp) || xw12_fxw21 <= xw2_fxw1 * (xlower - xwork) || xw12_fxw21 >= xw2_fxw1 * (xupper - xwork)) {
            xbdg_reset = true;
          } else {
            xbdg = xw12_fxw21 / xw2_fxw1;
            xwdg = xwork + xbdg;
            if (xwdg - xlower < tolx2 || xupper - xwdg < tolx2) {
              xbdg = copysign(tolx, xmid - xwork);
              if (aurostd::abs(xbdg) >= tolx) {
                xwdg = xwork + xbdg;
              } else {
                xwdg = xwork + copysign(tolx, xbdg);
              }
              xbdg_reset = false;
            }
          }
        }
        if (xbdg_reset) {
          if (xwork >= xmid) {
            xbound_diff = xlower - xwork;
          } else {
            xbound_diff = xupper - xwork;
          }
          xbdg = goldenratio_diff2 * xbound_diff;
        }
        if (aurostd::abs(xbdg) >= tolx) {
          xwdg = xwork + xbdg;
        } else {
          xwdg = xwork + copysign(tolx, xbdg);
        }
        aglerror = AGL_functions::bcnt_optimize_beta(xwdg, AGL_data, f_xwdg, FileMESSAGE);
        if (aglerror != 0) {
          return aglerror;
        }
        if (f_xwdg <= f_xwork) {
          if (xwdg >= xwork) {
            xlower = xwork;
          } else {
            xupper = xwork;
          }
          xwim2 = xwim1;
          f_xwim2 = f_xwim1;
          xwim1 = xwork;
          f_xwim1 = f_xwork;
          xwork = xwdg;
          f_xwork = f_xwdg;
        } else {
          if (xwdg < xwork) {
            xlower = xwdg;
          } else {
            xupper = xwdg;
          }
          if (f_xwdg <= f_xwim1 || aurostd::abs(xwim1 - xwork) < tole12) {
            xwim2 = xwim1;
            f_xwim2 = f_xwim1;
            xwim1 = xwdg;
            f_xwim1 = f_xwdg;
          } else if (f_xwdg <= f_xwim2 || aurostd::abs(xwim2 - xwork) < tole12 || aurostd::abs(xwim2 - xwork) < tole12) {
            xwim2 = xwdg;
            f_xwim2 = f_xwdg;
          }
        }
      }
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_WARNING_ + "bcnt_brent_minimum: exceeded maximum iterations." << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    xmin = xwork;
    brent_x = f_xwork;
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::bcnt_optimize_beta
// ***************************************************************************
namespace AGL_functions {
  //
  // bcnt_optimize_beta - optimization of exponent parameter "beta" required for BCNT EOS.
  //
  // Evaluates the function expression [log B + log K* - beta log (p - p_sp)]^2,
  // summed over p and B for all values for all pressures p.
  // Minimizing this expression (i.e. when this expression is equal to zero)
  // which is always positive since it is a sum of squares,
  // gives the set of values for which the equation of state is satisfied.
  // Uses least squares fit of this expression to optimize the value of "beta".
  //
  // Adapted from original Fortran version written by M. A. Blanco et al.
  // See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details
  //
  uint bcnt_optimize_beta(double& x_Press_sp, _AGL_data& AGL_data, double& lB_lK_lbppsp, ofstream& FileMESSAGE) {
    vector<double> xfunc(AGL_data.pressure_external.size());
    vector<double> yfunc(AGL_data.pressure_external.size());
    double mat_a_11;
    double mat_a_12;
    double mat_a_21;
    double mat_a_22;
    double z_1;
    double z_2;
    double det_mat_a;
    double log_K;
    double beta_fit;
    ostringstream aus;

    if (x_Press_sp < 0.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "bcnt_optimize_beta: Warning: Spinodal pressure is negative" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      lB_lK_lbppsp = 1e30;
      return 0;
    }
    if (AGL_data.optimize_beta) {
      mat_a_11 = 0.0;
      mat_a_12 = 0.0;
      mat_a_21 = 0.0;
      mat_a_22 = 0.0;
      z_1 = 0.0;
      z_2 = 0.0;
      if ((xfunc.size() != AGL_data.pressure_external.size()) || (yfunc.size() != AGL_data.pressure_external.size()) || (AGL_data.bulkmodulus.size() != AGL_data.pressure_external.size())) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Vectors are different size" << endl;
        aus << _AGLSTR_ERROR_ + "xfunc.size() = " << xfunc.size() << endl;
        aus << _AGLSTR_ERROR_ + "yfunc.size() = " << yfunc.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.bulkmodulus.size() = " << AGL_data.bulkmodulus.size() << endl;
        aus << _AGLSTR_ERROR_ + "AGL_data.pressure_external.size() = " << AGL_data.pressure_external.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 3;
      }
      for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
        xfunc[i] = log(AGL_data.pressure_external[i] + x_Press_sp);
        yfunc[i] = log(AGL_data.bulkmodulus[i]);
        mat_a_11 = mat_a_11 + 1.0;
        mat_a_12 = mat_a_12 + xfunc[i];
        mat_a_21 = mat_a_12;
        mat_a_22 = mat_a_22 + xfunc[i] * xfunc[i];
        z_1 = z_1 + yfunc[i];
        z_2 = z_2 + xfunc[i] * yfunc[i];
      }
      det_mat_a = mat_a_11 * mat_a_22 - mat_a_12 * mat_a_21;
      log_K = (z_1 * mat_a_22 - z_2 * mat_a_12) / det_mat_a;
      beta_fit = (mat_a_11 * z_2 - z_1 * mat_a_21) / det_mat_a;
      lB_lK_lbppsp = 0.0;
      for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
        lB_lK_lbppsp = lB_lK_lbppsp + pow(yfunc[i] - log_K - beta_fit * xfunc[i], 2);
      }
      AGL_data.x_K_opt = exp(-log_K);
      AGL_data.bcnt_beta = beta_fit;
      return 0;
    } else {
      mat_a_12 = 0.0;
      z_1 = 0.0;
      for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
        xfunc[i] = log(AGL_data.pressure_external[i] + x_Press_sp);
        yfunc[i] = log(AGL_data.bulkmodulus[i]);
        mat_a_12 = mat_a_12 + xfunc[i];
        z_1 = z_1 + yfunc[i];
      }
      log_K = (z_1 - AGL_data.bcnt_beta * mat_a_12) / AGL_data.pressure_external.size();
      AGL_data.x_K_opt = exp(-log_K);
      lB_lK_lbppsp = 0.0;
      for (size_t i = 0; i < AGL_data.pressure_external.size(); i++) {
        lB_lK_lbppsp = lB_lK_lbppsp + pow(yfunc[i] - log_K - AGL_data.bcnt_beta * xfunc[i], 2);
      }
    }

    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL Eqn of State
// **************************************************************************

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
