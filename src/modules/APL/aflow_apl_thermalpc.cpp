// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
//
// Calculates thermal properties from phonon densities of states. Originally
// written by Jahnatek and adapted/rewritten by Marco Esters for use with POCC.
//
// The old code assumed that there was one DOS for the entire temperature range.
// However, POCC has different DOS for different temperatures.
//
// Properties (all units per cell):
// U0:    zero point energy (in meV)
// U:     internal energy (in meV)
// Fvib:  vibrational free energy (in meV)
// Cv:    isochoric heat capacity (in kB)
// Svib:  vibrational entropy (in kB)
//
// For the calculations, frequencies are assumed to be in THz divided by 2pi (as
// is standard for the APL DOS calculator).

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"

#include "aflow.h"
#include "aflow_apl.h"
#include "aflow_defs.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"

using std::ofstream;
using std::ostream;
using std::string;
using std::stringstream;
using std::vector;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  ThermalPropertiesCalculator::ThermalPropertiesCalculator(const DOSCalculator& dosc, ofstream& mf, const string& directory, ostream& oss) {
    _directory = aurostd::CleanFileName(directory);
    initialize(dosc.createDOSCAR(), mf, oss);
  }

  ThermalPropertiesCalculator::ThermalPropertiesCalculator(const xDOSCAR& xdos, ofstream& mf, const string& directory, ostream& oss) {
    _directory = aurostd::CleanFileName(directory);
    initialize(xdos, mf, oss);
  }

  void ThermalPropertiesCalculator::clear() {
    _freqs_0K.clear();
    _dos_0K.clear();
    _directory = "";
    natoms = 0;
    system = "";
    temperatures.clear();
    Cv.clear();
    Fvib.clear();
    Svib.clear();
    U.clear();
    U0 = 0.0;
  }
}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                 OVERHEAD                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // initialize////////////////////////////////////////////////////////////////
  //  Initializes the thermal properties calculator with a 0 K solution.
  void ThermalPropertiesCalculator::initialize(const xDOSCAR& xdos, ofstream& mf, ostream& oss) {
    natoms = xdos.number_atoms;
    vector<double> freq = aurostd::deque2vector(xdos.venergy);
    // Convert to THz
    for (size_t i = 0; i < freq.size(); i++) {
      freq[i] *= eV2Hz * Hz2THz;
    }
    const vector<double> dos = aurostd::deque2vector(xdos.vDOS[0][0][0]);
    initialize(freq, dos, mf, xdos.title, oss);
  }

  void ThermalPropertiesCalculator::initialize(const vector<double>& freqs, const vector<double>& dos, ofstream& mf, const string& _system, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(freqs, dos, _system);
  }

  void ThermalPropertiesCalculator::initialize(const vector<double>& freqs, const vector<double>& dos, const string& _system) {
    _freqs_0K = freqs;
    _dos_0K = dos;
    system = _system;
    U0 = getZeroPointEnergy();
  }

  // calculateThermalProperties////////////////////////////////////////////////
  //  Calculates thermal properties within a desired temperature range using
  //  the 0 K density of states.
  void ThermalPropertiesCalculator::calculateThermalProperties(double Tstart, double Tend, double Tstep) {
    string message = "Calculating thermal properties.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);
    if (Tstart > Tend) {
      message = "Tstart cannot be higher than Tend.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }

    temperatures.clear();
    U.clear();
    Fvib.clear();
    Svib.clear();
    Cv.clear();
    for (double T = Tstart; T <= Tend; T += Tstep) {
      addPoint(T, _freqs_0K, _dos_0K);
    }
  }

  // addPoint//////////////////////////////////////////////////////////////////
  //  Adds a temperature data point to the thermal properties. This function
  //  is especially useful when each temperature has a different DOS.
  void ThermalPropertiesCalculator::addPoint(double T, const xDOSCAR& xdos) {
    vector<double> freq = aurostd::deque2vector(xdos.venergy);
    // Convert to THz
    for (size_t i = 0; i < freq.size(); i++) {
      freq[i] *= eV2Hz * Hz2THz;
    }
    const vector<double> dos = aurostd::deque2vector(xdos.vDOS[0][0][0]);
    addPoint(T, freq, dos);
  }

  void ThermalPropertiesCalculator::addPoint(double T, const vector<double>& freq, const vector<double>& dos) {
    temperatures.push_back(T);
    U.push_back(getInternalEnergy(T, freq, dos));
    Fvib.push_back(getVibrationalFreeEnergy(T, freq, dos));
    Svib.push_back(getVibrationalEntropy(T, U.back(), Fvib.back()));
    Cv.push_back(getIsochoricSpecificHeat(T, freq, dos));
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                PROPERTIES                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// Except for the zero point energy, each function is overloaded to use another
// DOS as the input. This is useful for cases that have different DOS for each
// temperature.

namespace apl {

  // getZeroPointEnergy////////////////////////////////////////////////////////
  //  Calculates the zero point (internal) energy from the 0 K DOS.
  double ThermalPropertiesCalculator::getZeroPointEnergy() {
    double zpe = 0.0;
    const double stepDOS = getStepDOS(_freqs_0K);
    for (size_t i = 0; i < _freqs_0K.size(); i++) {
      if (_freqs_0K[i] < _FLOAT_TOL_) {
        continue;
      }
      zpe += _freqs_0K[i] * THz2Hz * _dos_0K[i];
    }
    zpe *= 0.5 * 1000 * PLANCKSCONSTANTEV_h * stepDOS;  // Convert to meV
    return zpe;
  }

  // getInternalEnergy/////////////////////////////////////////////////////////
  //  Calculates the internal energy.
  double ThermalPropertiesCalculator::getInternalEnergy(double T, ThermalPropertiesUnits unit) {
    return getInternalEnergy(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getInternalEnergy(double T, const vector<double>& freq, const vector<double>& dos, ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) {
      return getScalingFactor(unit) * U0; // AS20200508
    }

    const double stepDOS = getStepDOS(freq);
    const double beta = 1.0 / (KBOLTZEV * T);  // beta = 1/kBT (in 1/eV)

    double E = 0.0;
    double hni = 0;
    for (size_t i = 0; i < freq.size(); i++) {
      if (freq[i] < _FLOAT_TOL_) {
        continue;
      }
      hni = PLANCKSCONSTANTEV_h * freq[i] * THz2Hz;  // h * freq in eV
      E += dos[i] * hni / (exp(beta * hni) - 1.0);
    }
    E *= 1000 * stepDOS;  // Convert to meV
    E += U0;
    return getScalingFactor(unit) * E;
  }

  // getVibrationalFreeEnergy//////////////////////////////////////////////////
  //  Calculates the vibrational free energy using the zero point energy, which
  //  is faster than calculating from scratch.
  double ThermalPropertiesCalculator::getVibrationalFreeEnergy(double T, ThermalPropertiesUnits unit) {
    return getVibrationalFreeEnergy(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getVibrationalFreeEnergy(double T, const vector<double>& freq, const vector<double>& dos, ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) {
      return getScalingFactor(unit) * U0; // AS20200508
    }

    const double stepDOS = getStepDOS(freq);
    const double beta = 1.0 / (KBOLTZEV * T);  // beta = 1/kBT (in 1/eV)

    double F = 0.0;
    double hni = 0.0;
    for (size_t i = 0; i < freq.size(); i++) {
      if (freq[i] < _FLOAT_TOL_) {
        continue;
      }
      hni = PLANCKSCONSTANTEV_h * freq[i] * THz2Hz;  // h * freq in eV
      F += dos[i] * aurostd::ln(1.0 - exp(-beta * hni)) / beta;
    }
    F *= 1000 * stepDOS;  // Convert to meV
    F += U0;
    return getScalingFactor(unit) * F;
  }

  // getVibrationalEntropy/////////////////////////////////////////////////////
  //  Calculates the vibrational entropy using the free energy and the internal
  //  energy (faster than calculating directly if they are available).
  double ThermalPropertiesCalculator::getVibrationalEntropy(double T, ThermalPropertiesUnits unit) {
    return getVibrationalEntropy(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getVibrationalEntropy(double T, const vector<double>& freq, const vector<double>& dos, ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) {
      return 0.0;
    }

    const double E = getInternalEnergy(T, freq, dos);
    const double F = getVibrationalFreeEnergy(T, freq, dos);

    return getVibrationalEntropy(T, E, F, unit);
  }

  double ThermalPropertiesCalculator::getVibrationalEntropy(double T, double E, double F, ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) {
      return 0.0;
    }

    const double S = (E - F) / T;
    return getScalingFactor(unit) * S;
  }

  // getIsochoricSpecificHeat//////////////////////////////////////////////////
  //  Calculates the isochoric heat capacity.
  double ThermalPropertiesCalculator::getIsochoricSpecificHeat(double T, ThermalPropertiesUnits unit) {
    return getIsochoricSpecificHeat(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getIsochoricSpecificHeat(double T, const vector<double>& freq, const vector<double>& dos, ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) {
      return 0.0;
    }

    const double stepDOS = getStepDOS(freq);
    const double beta = 1.0 / (KBOLTZEV * T);  // beta = 1/kBT (in 1/eV)

    double cv = 0.0;
    double bhni = 0.0;
    double ebhni = 0.0;
    for (size_t i = 0; i < freq.size(); i++) {
      if (freq[i] < _FLOAT_TOL_) {
        continue;
      }
      bhni = beta * PLANCKSCONSTANTEV_h * freq[i] * THz2Hz;  // h * freq/(kB * T)
      ebhni = exp(bhni);
      cv += dos[i] * (KBOLTZEV * bhni * bhni / ((1.0 - 1.0 / ebhni) * (ebhni - 1.0)));
    }
    if (std::isnan(cv)) {
      return 0.0;  // ME20200428 - the only problem here is when T = 0 (checked above), phase transitions not captured by the model
    }
    cv *= 1000.0 * stepDOS;  // Convert to meV
    return getScalingFactor(unit) * cv;
  }

  // getStepDOS////////////////////////////////////////////////////////////////
  //  Calculates the step size of the DOS. While it may be redundant to do
  //  this for every temperature when the 0 K DOS is used, it is essential when
  //  using different DOS for each temperature.
  double ThermalPropertiesCalculator::getStepDOS(const vector<double>& freq) {
    double stepDOS = 0.0;
    if (freq.size() > 2) {
      stepDOS = freq[1] - freq[0];
    } else {
      const string message = "Not enough DOS points (need at least two).";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return stepDOS;
  }

  // getScalingFactor//////////////////////////////////////////////////////////
  //  Used to convert meV (default unit) into another unit.
  double ThermalPropertiesCalculator::getScalingFactor(const ThermalPropertiesUnits& units) {
    switch (units) {
      case eV:
      case eVK: return 0.001;
      case meV:
      case meVK: return 1.0;
      case ueV:
      case ueVK: return 1000.0;
      case kB: return 1.0 / (1000 * KBOLTZEV);
    }
    return 1.0;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                FILE I/O                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // writePropertiesToFile/////////////////////////////////////////////////////
  //  Outputs the thermal properties into a file that can be plotted using the
  //  AFLOW plotter.
  void ThermalPropertiesCalculator::writePropertiesToFile(string filename, filetype ft) {
    filename = aurostd::CleanFileName(filename);
    string message = "Writing thermal properties into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    stringstream outfile;
    if (ft == json_ft) {
      aurostd::JSON::object json = getPropertiesJSON();
      if (!system.empty()) {
        json["system"] = system;
      }
      outfile << json.toString();
    } else {
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      if (!system.empty()) {
        outfile << "[APL_THERMO]SYSTEM=" << system << std::endl;
      }
      outfile << getPropertiesFileString(ft);
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    }
    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  // ME20210927
  void ThermalPropertiesCalculator::addToAPLOut(stringstream& apl_outfile) {
    apl_outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    apl_outfile << "[APL_THERMO_RESULTS]START" << std::endl;
    apl_outfile << "energy_zero_point_cell_apl=" << std::setprecision(8) << U0 << " (meV/cell)" << std::endl;
    if (natoms > 0) {
      apl_outfile << "energy_zero_point_atom_apl=" << std::setprecision(8) << (U0 / natoms) << " (meV/atom)" << std::endl;
    }
    for (size_t t = 0; t < temperatures.size(); t++) {
      if (aurostd::isequal(temperatures[t], 300.0)) {
        apl_outfile << "energy_free_vibrational_cell_apl_300K=" << std::setprecision(8) << Fvib[t] << " (meV/cell)" << std::endl;
        if (natoms > 0) {
          apl_outfile << "energy_free_vibrational_atom_apl_300K=" << std::setprecision(8) << (Fvib[t] / (double) natoms) << " (meV/atom)" << std::endl;
        }
        apl_outfile << "entropy_vibrational_cell_apl_300K=" << std::setprecision(8) << Svib[t] << " (kB/cell)" << std::endl;
        if (natoms > 0) {
          apl_outfile << "entropy_vibrational_atom_apl_300K=" << std::setprecision(8) << (Svib[t] / (double) natoms) << " (kB/atom)" << std::endl;
        }
        apl_outfile << "energy_internal_vibrational_cell_apl_300K=" << std::setprecision(8) << U[t] << " (meV/cell)" << std::endl;
        if (natoms > 0) {
          apl_outfile << "energy_internal_vibrational_atom_apl_300K=" << std::setprecision(8) << (U[t] / (double) natoms) << " (meV/atom)" << std::endl;
        }
        apl_outfile << "heat_capacity_Cv_cell_apl_300K=" << std::setprecision(5) << Cv[t] << " (kB/cell)" << std::endl;
        if (natoms > 0) {
          apl_outfile << "heat_capacity_Cv_atom_apl_300K=" << std::setprecision(5) << (Cv[t] / (double) natoms) << " (kB/atom)" << std::endl;
        }
        break;
      }
    }
    apl_outfile << "[APL_THERMO_RESULTS]STOP" << std::endl;
    apl_outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if (!system.empty()) {
      apl_outfile << "[APL_THERMO]SYSTEM=" << system << std::endl;
    }
    apl_outfile << getPropertiesFileString();
    apl_outfile << AFLOWIN_SEPARATION_LINE << std::endl;
  }

  aurostd::JSON::object ThermalPropertiesCalculator::getPropertiesJSON() {
    using namespace aurostd::JSON;
    object json(object_types::DICTIONARY);
    json["T"] = temperatures;
    json["U0_cell"] = object(object_types::DICTIONARY);
    json["U0_cell"]["unit"] = "meV";
    json["U0_cell"]["unit_latex"] = "meV";
    json["U0_cell"]["unit_html"] = "meV";
    json["U0_cell"]["value"] = U0;
    json["U_cell"] = object(object_types::DICTIONARY);
    json["U_cell"]["unit"] = "meV";
    json["U_cell"]["unit_latex"] = "meV";
    json["U_cell"]["unit_html"] = "meV";
    json["U_cell"]["value"] = U;
    json["Fvib_cell"] = object(object_types::DICTIONARY);
    json["Fvib_cell"]["unit"] = "meV";
    json["Fvib_cell"]["unit_latex"] = "meV";
    json["Fvib_cell"]["unit_html"] = "meV";
    json["Fvib_cell"]["value"] = Fvib;
    json["Svib_cell"] = object(object_types::DICTIONARY);
    json["Svib_cell"]["unit"] = "kB";
    json["Svib_cell"]["unit_latex"] = "$k_\\\\textnormal{B}$";
    json["Svib_cell"]["unit_html"] = "<i>k</i><sub>B</sub>";
    json["Svib_cell"]["value"] = Svib;
    json["Cv_cell"] = object(object_types::DICTIONARY);
    json["Cv_cell"]["unit"] = "kB";
    json["Cv_cell"]["unit_latex"] = "$k_\\\\textnormal{B}$";
    json["Cv_cell"]["unit_html"] = "<i>k</i><sub>B</sub>";
    json["Cv_cell"]["value"] = Svib;

    if (natoms > 0) {
      const uint ntemps = temperatures.size();
      vector<double> val_per_atom(ntemps);

      json["U0_atom"] = object(object_types::DICTIONARY);
      json["U0_atom"]["unit"] = "meV";
      json["U0_atom"]["unit_latex"] = "meV";
      json["U0_atom"]["unit_html"] = "meV";
      json["U0_atom"]["value"] = U0 / ((double) natoms);

      json["U_atom"] = object(object_types::DICTIONARY);
      json["U_atom"]["unit"] = "meV";
      json["U_atom"]["unit_latex"] = "meV";
      json["U_atom"]["unit_html"] = "meV";
      for (uint i = 0; i < ntemps; i++) {
        val_per_atom[i] = U[i] / ((double) natoms);
      }
      json["U_atom"]["value"] = val_per_atom;

      json["Fvib_atom"] = object(object_types::DICTIONARY);
      json["Fvib_atom"]["unit"] = "meV";
      json["Fvib_atom"]["unit_latex"] = "meV";
      json["Fvib_atom"]["unit_html"] = "meV";
      for (uint i = 0; i < ntemps; i++) {
        val_per_atom[i] = Fvib[i] / ((double) natoms);
      }
      json["Fvib_atom"]["value"] = val_per_atom;

      json["Svib_atom"] = object(object_types::DICTIONARY);
      json["Svib_atom"]["unit"] = "kB";
      json["Svib_atom"]["unit_latex"] = "$k_\\\\textnormal{B}$";
      json["Svib_atom"]["unit_html"] = "<i>k</i><sub>B</sub>";
      for (uint i = 0; i < ntemps; i++) {
        val_per_atom[i] = Svib[i] / ((double) natoms);
      }
      json["Svib_atom"]["value"] = val_per_atom;

      json["Cv_atom"] = object(object_types::DICTIONARY);
      json["Cv_atom"]["unit"] = "kB";
      json["Cv_atom"]["unit_latex"] = "$k_\\\\textnormal{B}$";
      json["Cv_atom"]["unit_html"] = "<i>k</i><sub>B</sub>";
      for (uint i = 0; i < ntemps; i++) {
        val_per_atom[i] = Cv[i] / ((double) natoms);
      }
      json["Cv_atom"]["value"] = val_per_atom;
    }
    return json;
  }

  string ThermalPropertiesCalculator::getPropertiesFileString(filetype ft) {
    stringstream props;

    if (ft == json_ft) {
      props << getPropertiesJSON().toString();
    } else {
      // Header
      props << "[APL_THERMO]START" << std::endl;
      props << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      props << "#" << std::setw(7) << "T(K)" << std::setw(15) << "U0 (meV/cell)" << "   " << std::setw(15) << "U (meV/cell)" << "   " << std::setw(15) << "F (meV/cell)" << "   " << std::setw(15)
            << "S (kB/cell)" << "   " << std::setw(15) << "Cv (kB/cell)";
      if (natoms > 0) {
        props << std::setw(15) << "U0 (meV/atom)" << "   " << std::setw(15) << "U (meV/atom)" << "   " << std::setw(15) << "F (meV/atom)" << "   " << std::setw(15) << "S (kB/atom)" << "   " << std::setw(15)
              << "Cv (kB/atom)";
      }
      props << std::endl;

      for (size_t t = 0; t < temperatures.size(); t++) {
        props << std::setw(8) << std::setprecision(2) << temperatures[t] << std::setprecision(8) << std::setw(15) << U0 << "   " << std::setw(15) << U[t] << "   " << std::setw(15) << Fvib[t] << "   "
              << std::setw(15) << Svib[t] << "   " << std::setw(15) << Cv[t];
        if (natoms > 0) {
          props << std::setw(15) << (U0 / (double) natoms) << "   " << std::setw(15) << (U[t] / (double) natoms) << "   " << std::setw(15) << (Fvib[t] / (double) natoms) << "   " << std::setw(15)
                << (Svib[t] / (double) natoms) << "   " << std::setw(15) << (Cv[t] / (double) natoms);
        }
        props << std::endl;
      }

      // Footer
      props << "[APL_THERMO]STOP" << std::endl;
    }

    return props.str();
  }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
