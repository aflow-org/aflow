// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
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

#include "aflow_apl.h"

using std::vector;
using std::string;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Default Constructor
  ThermalPropertiesCalculator::ThermalPropertiesCalculator(ostream& oss): xStream(oss) {
    free();
  }

  ThermalPropertiesCalculator::ThermalPropertiesCalculator(ofstream& mf, ostream& oss) : xStream(mf,oss) {
    free();
  }

  ThermalPropertiesCalculator::ThermalPropertiesCalculator(const DOSCalculator& dosc, ofstream& mf, const string& directory, ostream& oss) {
    free();
    _directory = aurostd::CleanFileName(directory);
    initialize(dosc.createDOSCAR(), mf, oss);
  }

  ThermalPropertiesCalculator::ThermalPropertiesCalculator(const xDOSCAR& xdos, ofstream& mf, const string& directory, ostream& oss) {
    free();
    _directory = aurostd::CleanFileName(directory);
    initialize(xdos, mf, oss);
  }

  // Copy constructors
  ThermalPropertiesCalculator::ThermalPropertiesCalculator(const ThermalPropertiesCalculator& that) : xStream(*that.getOFStream(),*that.getOSS()) {
    if (this != &that) free();
    copy(that);
  }

  ThermalPropertiesCalculator& ThermalPropertiesCalculator::operator=(const ThermalPropertiesCalculator& that) {
    if (this != &that) free();
    copy(that);
    return *this;
  }

  // Destructor
  ThermalPropertiesCalculator::~ThermalPropertiesCalculator() {
    xStream::free();
    free();
  }

  void ThermalPropertiesCalculator::copy(const ThermalPropertiesCalculator& that) {
    if (this == &that) return;
    xStream::copy(that);
    _freqs_0K = that._freqs_0K;
    _dos_0K = that._dos_0K;
    _directory = that._directory;
    natoms = that.natoms;
    system = that.system;
    temperatures = that.temperatures;
    Cv = that.Cv;
    Fvib = that.Fvib;
    Svib = that.Svib;
    U = that.U;
    U0 = that.U0;
  }


  void ThermalPropertiesCalculator::free() {
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

  void ThermalPropertiesCalculator::clear() {
    free();
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                 OVERHEAD                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //initialize////////////////////////////////////////////////////////////////
  // Initializes the thermal properties calculator with a 0 K solution.
  void ThermalPropertiesCalculator::initialize(const xDOSCAR& xdos, ofstream& mf, ostream& oss) {
    natoms = xdos.number_atoms;
    vector<double> freq = aurostd::deque2vector(xdos.venergy);
    // Convert to THz
    for (uint i = 0; i < freq.size(); i++) freq[i] *= eV2Hz * Hz2THz;
    vector<double> dos = aurostd::deque2vector(xdos.vDOS[0][0][0]);
    initialize(freq, dos, mf, xdos.title, oss);
  }

  void ThermalPropertiesCalculator::initialize(const vector<double>& freqs,
      const vector<double>& dos, ofstream& mf, const string& _system, ostream& oss) {
    xStream::initialize(mf, oss);
    initialize(freqs, dos, _system);
  }

  void ThermalPropertiesCalculator::initialize(const vector<double>& freqs,
      const vector<double>& dos, const string& _system) {
    _freqs_0K = freqs;
    _dos_0K = dos;
    system = _system;
    U0 = getZeroPointEnergy();
  }

  //xStream initializers
  void ThermalPropertiesCalculator::initialize(ostream& oss) {
    xStream::initialize(oss);
  }

  void ThermalPropertiesCalculator::initialize(ofstream& mf, ostream& oss) {
    xStream::initialize(mf, oss);
  }

  //calculateThermalProperties////////////////////////////////////////////////
  // Calculates thermal properties within a desired temperature range using
  // the 0 K density of states.
  void ThermalPropertiesCalculator::calculateThermalProperties(double Tstart,
      double Tend,
      double Tstep) {
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
    for (double T = Tstart; T <= Tend; T += Tstep) addPoint(T, _freqs_0K, _dos_0K);
  }

  //addPoint//////////////////////////////////////////////////////////////////
  // Adds a temperature data point to the thermal properties. This function
  // is especially useful when each temperature has a different DOS.
  void ThermalPropertiesCalculator::addPoint(double T,
      const xDOSCAR& xdos) {
    vector<double> freq = aurostd::deque2vector(xdos.venergy);
    // Convert to THz
    for (uint i = 0; i < freq.size(); i++) freq[i] *= eV2Hz * Hz2THz;
    vector<double> dos = aurostd::deque2vector(xdos.vDOS[0][0][0]);
    addPoint(T, freq, dos);
  }

  void ThermalPropertiesCalculator::addPoint(double T,
      const vector<double>& freq,
      const vector<double>& dos) {
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

  //getZeroPointEnergy////////////////////////////////////////////////////////
  // Calculates the zero point (internal) energy from the 0 K DOS.
  double ThermalPropertiesCalculator::getZeroPointEnergy() {
    double zpe = 0.0;
    double stepDOS = getStepDOS(_freqs_0K);
    for (uint i = 0; i < _freqs_0K.size(); i++) {
      if (_freqs_0K[i] < _FLOAT_TOL_) continue;
      zpe += _freqs_0K[i] * THz2Hz * _dos_0K[i];
    }
    zpe *= 0.5 * 1000 * PLANCKSCONSTANTEV_h * stepDOS;  // Convert to meV
    return zpe;
  }

  //getInternalEnergy/////////////////////////////////////////////////////////
  // Calculates the internal energy.
  double ThermalPropertiesCalculator::getInternalEnergy(double T,
      ThermalPropertiesUnits unit) {
    return getInternalEnergy(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getInternalEnergy(double T,
      const vector<double>& freq,
      const vector<double>& dos,
      ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) return getScalingFactor(unit) * U0; //AS20200508

    double stepDOS = getStepDOS(freq);
    double beta = 1.0/(KBOLTZEV * T);  // beta = 1/kBT (in 1/eV)

    double E = 0.0, hni = 0;
    for (uint i = 0; i < freq.size(); i++) {
      if (freq[i] < _FLOAT_TOL_) continue;
      hni = PLANCKSCONSTANTEV_h * freq[i] * THz2Hz;  // h * freq in eV
      E += dos[i] * hni / (exp(beta * hni) - 1.0);
    }
    E *= 1000 * stepDOS;  // Convert to meV
    E += U0;
    return getScalingFactor(unit) * E;
  }

  //getVibrationalFreeEnergy//////////////////////////////////////////////////
  // Calculates the vibrational free energy using the zero point energy, which
  // is faster than calculating from scratch.
  double ThermalPropertiesCalculator::getVibrationalFreeEnergy(double T,
      ThermalPropertiesUnits unit) {
    return getVibrationalFreeEnergy(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getVibrationalFreeEnergy(double T,
      const vector<double>& freq,
      const vector<double>& dos,
      ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) return getScalingFactor(unit) * U0; //AS20200508

    double stepDOS = getStepDOS(freq);
    double beta = 1.0/(KBOLTZEV * T);  // beta = 1/kBT (in 1/eV)

    double F = 0.0, hni = 0.0;
    for (uint i = 0; i < freq.size(); i++) {
      if (freq[i] < _FLOAT_TOL_) continue;
      hni = PLANCKSCONSTANTEV_h * freq[i] * THz2Hz;  // h * freq in eV
      F += dos[i] * aurostd::ln(1.0 - exp(-beta * hni)) / beta;
    }
    F *= 1000 * stepDOS;  // Convert to meV
    F += U0;
    return getScalingFactor(unit) * F;
  }


  //getVibrationalEntropy/////////////////////////////////////////////////////
  // Calculates the vibrational entropy using the free energy and the internal
  // energy (faster than calculating directly if they are available).
  double ThermalPropertiesCalculator::getVibrationalEntropy(double T,
      ThermalPropertiesUnits unit) {
    return getVibrationalEntropy(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getVibrationalEntropy(double T,
      const vector<double>& freq,
      const vector<double>& dos,
      ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) return 0.0;

    double E = getInternalEnergy(T, freq, dos);
    double F = getVibrationalFreeEnergy(T, freq, dos);

    return getVibrationalEntropy(T, E, F, unit);
  }

  double ThermalPropertiesCalculator::getVibrationalEntropy(double T, double E,
      double F, ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) return 0.0;

    double S = (E - F)/T;
    return getScalingFactor(unit) * S;
  }

  //getIsochoricSpecificHeat//////////////////////////////////////////////////
  // Calculates the isochoric heat capacity.
  double ThermalPropertiesCalculator::getIsochoricSpecificHeat(double T,
      ThermalPropertiesUnits unit) {
    return getIsochoricSpecificHeat(T, _freqs_0K, _dos_0K, unit);
  }

  double ThermalPropertiesCalculator::getIsochoricSpecificHeat(double T,
      const vector<double>& freq,
      const vector<double>& dos,
      ThermalPropertiesUnits unit) {
    if (T < _FLOAT_TOL_) return 0.0;

    double stepDOS = getStepDOS(freq);
    double beta = 1.0/(KBOLTZEV * T);  // beta = 1/kBT (in 1/eV)

    double cv = 0.0, bhni = 0.0, ebhni = 0.0;
    for (uint i = 0; i < freq.size(); i++) {
      if (freq[i] < _FLOAT_TOL_) continue;
      bhni = beta * PLANCKSCONSTANTEV_h * freq[i] * THz2Hz;  // h * freq/(kB * T)
      ebhni = exp(bhni);
      cv += dos[i] * (KBOLTZEV * bhni * bhni / ((1.0 - 1.0 / ebhni) * (ebhni - 1.0)));
    }
    if (std::isnan(cv)) return 0.0;  //ME20200428 - the only problem here is when T = 0 (checked above), phase transitions not captured by the model
    cv *= 1000.0 * stepDOS;  // Convert to meV
    return getScalingFactor(unit) * cv;
  }

  //getStepDOS////////////////////////////////////////////////////////////////
  // Calculates the step size of the DOS. While it may be redundant to do
  // this for every temperature when the 0 K DOS is used, it is essential when
  // using different DOS for each temperature.
  double ThermalPropertiesCalculator::getStepDOS(const vector<double>& freq) {
    double stepDOS = 0.0;
    if (freq.size() > 2) {
      stepDOS = freq[1] - freq[0];
    } else {
      string message = "Not enough DOS points (need at least two).";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return stepDOS;
  }

  //getScalingFactor//////////////////////////////////////////////////////////
  // Used to convert meV (default unit) into another unit.
  double ThermalPropertiesCalculator::getScalingFactor(const ThermalPropertiesUnits& units) {
    switch (units) {
      case eV:
      case eVK:
        return 0.001;
      case meV:
      case meVK:
        return 1.0;
      case ueV:
      case ueVK:
        return 1000.0;
      case kB:
        return 1.0/(1000 * KBOLTZEV);
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

  //writePropertiesToFile/////////////////////////////////////////////////////
  // Outputs the thermal properties into a file that can be plotted using the
  // AFLOW plotter.
  void ThermalPropertiesCalculator::writePropertiesToFile(string filename, filetype ft) {
    filename = aurostd::CleanFileName(filename);
    string message = "Writing thermal properties into file " + filename + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _directory, *p_FileMESSAGE, *p_oss);

    stringstream outfile;
    if (ft == json_ft) {
      aurostd::JSONwriter json;
      if (!system.empty()) json.addString("system", system);
      json.mergeRawJSON(getPropertiesFileString(ft));
      outfile << json.toString();
    } else {
      outfile << AFLOWIN_SEPARATION_LINE << std::endl;
      if (!system.empty()) outfile << "[APL_THERMO]SYSTEM=" << system << std::endl;
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
    if (natoms > 0) apl_outfile << "energy_zero_point_atom_apl=" << std::setprecision(8) << (U0/natoms) << " (meV/atom)" << std::endl;
    for (uint t = 0; t < temperatures.size(); t++) {
      if (aurostd::isequal(temperatures[t], 300.0)) {
        apl_outfile << "energy_free_vibrational_cell_apl_300K=" << std::setprecision(8) << Fvib[t] << " (meV/cell)" << std::endl;
        if (natoms > 0) apl_outfile << "energy_free_vibrational_atom_apl_300K=" << std::setprecision(8) << (Fvib[t]/(double) natoms) << " (meV/atom)" << std::endl;
        apl_outfile << "entropy_vibrational_cell_apl_300K=" << std::setprecision(8) << Svib[t] << " (kB/cell)" << std::endl;
        if (natoms > 0) apl_outfile << "entropy_vibrational_atom_apl_300K=" << std::setprecision(8) << (Svib[t]/(double) natoms) << " (kB/atom)" << std::endl;
        apl_outfile << "energy_internal_vibrational_cell_apl_300K=" << std::setprecision(8) << U[t] << " (meV/cell)" << std::endl;
        if (natoms > 0) apl_outfile << "energy_internal_vibrational_atom_apl_300K=" << std::setprecision(8) << (U[t]/(double) natoms) << " (meV/atom)" << std::endl;
        apl_outfile << "heat_capacity_Cv_cell_apl_300K=" << std::setprecision(5) << Cv[t] << " (kB/cell)" << std::endl;
        if (natoms > 0) apl_outfile << "heat_capacity_Cv_atom_apl_300K=" << std::setprecision(5) << (Cv[t]/(double) natoms) << " (kB/atom)" << std::endl;
        break;
      }
    }
    apl_outfile << "[APL_THERMO_RESULTS]STOP" << std::endl;
    apl_outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    if (!system.empty()) apl_outfile << "[APL_THERMO]SYSTEM=" << system << std::endl;
    apl_outfile << getPropertiesFileString();
    apl_outfile << AFLOWIN_SEPARATION_LINE << std::endl;
  }

  string ThermalPropertiesCalculator::getPropertiesFileString(filetype ft) {
    stringstream props;

    if (ft == json_ft) {
      aurostd::JSONwriter json,property;
      json.addVector("T", temperatures);
      property.clear();
      property.addString("unit", "meV");
      property.addString("unit_latex", "meV");
      property.addString("unit_html", "meV");
      property.addNumber("value", U0);
      json.addJSON("U0_cell", property);
      property.clear();
      property.addString("unit", "meV");
      property.addString("unit_latex", "meV");
      property.addString("unit_html", "meV");
      property.addVector("value", U);
      json.addJSON("U_cell", property);
      property.clear();
      property.addString("unit", "meV");
      property.addString("unit_latex", "meV");
      property.addString("unit_html", "meV");
      property.addVector("value", Fvib);
      json.addJSON("Fvib_cell", property);
      property.clear();
      property.addString("unit", "kB");
      property.addString("unit_latex", "$k_\\\\textnormal{B}$");
      property.addString("unit_html", "<i>k</i><sub>B</sub>");
      property.addVector("value", Svib);
      json.addJSON("Svib_cell", property);
      property.clear();
      property.addString("unit", "kB");
      property.addString("unit_latex", "$k_\\\\textnormal{B}$");
      property.addString("unit_html", "<i>k</i><sub>B</sub>");
      property.addVector("value", Cv);
      json.addJSON("Cv_cell", property);
      if (natoms > 0) {
        property.clear();
        property.addString("unit", "meV");
        property.addString("unit_latex", "meV");
        property.addString("unit_html", "meV");
        property.addNumber("value", (U0/(double) natoms));
        json.addJSON("U0_atom", property);
        uint ntemps = temperatures.size();
        vector<double> val_per_atom(ntemps);
        property.clear();
        property.addString("unit", "meV");
        property.addString("unit_latex", "meV");
        property.addString("unit_html", "meV");
        for (uint i = 0; i < ntemps; i++) val_per_atom[i] = U[i]/((double) natoms);
        property.addVector("value", val_per_atom);
        property.addJSON("U_atom", property);
        property.clear();
        property.addString("unit", "meV");
        property.addString("unit_latex", "meV");
        property.addString("unit_html", "meV");
        for (uint i = 0; i < ntemps; i++) val_per_atom[i] = Fvib[i]/((double) natoms);
        property.addVector("value", val_per_atom);
        property.addJSON("Fvib_atom", property);
        property.clear();
        property.addString("unit", "kB");
        property.addString("unit_latex", "$k_\\\\textnormal{B}$");
        property.addString("unit_html", "<i>k</i><sub>B</sub>");
        for (uint i = 0; i < ntemps; i++) val_per_atom[i] = Svib[i]/((double) natoms);
        property.addVector("value", val_per_atom);
        property.addJSON("Svib_atom", property);
        property.clear();
        property.addString("unit", "kB");
        property.addString("unit_latex", "$k_\\\\textnormal{B}$");
        property.addString("unit_html", "<i>k</i><sub>B</sub>");
        for (uint i = 0; i < ntemps; i++) val_per_atom[i] = Cv[i]/((double) natoms);
        property.addVector("value", val_per_atom);
        property.addJSON("Cv_atom", property);
      }
      props << json.toString(false);
    } else {
      // Header
      props << "[APL_THERMO]START" << std::endl;
      props << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      props << "#"
        << std::setw(7) << "T(K)"
        << std::setw(15) << "U0 (meV/cell)" << "   "
        << std::setw(15) << "U (meV/cell)" << "   "
        << std::setw(15) << "F (meV/cell)" << "   "
        << std::setw(15) << "S (kB/cell)" << "   "
        << std::setw(15) << "Cv (kB/cell)";
      if (natoms > 0) {
        props << std::setw(15) << "U0 (meV/atom)" << "   "
          << std::setw(15) << "U (meV/atom)" << "   "
          << std::setw(15) << "F (meV/atom)" << "   "
          << std::setw(15) << "S (kB/atom)" << "   "
          << std::setw(15) << "Cv (kB/atom)";
      }
      props << std::endl;

      for (uint t = 0; t < temperatures.size(); t++) {
        props << std::setw(8) << std::setprecision(2) << temperatures[t]
          << std::setprecision(8)
          << std::setw(15) << U0 << "   "
          << std::setw(15) << U[t] << "   "
          << std::setw(15) << Fvib[t] << "   "
          << std::setw(15) << Svib[t] << "   "
          << std::setw(15) << Cv[t];
        if (natoms > 0) {
          props  << std::setw(15) << (U0/(double) natoms) << "   "
            << std::setw(15) << (U[t]/(double) natoms) << "   "
            << std::setw(15) << (Fvib[t]/(double) natoms) << "   "
            << std::setw(15) << (Svib[t]/(double) natoms) << "   "
            << std::setw(15) << (Cv[t]/(double) natoms);
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
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
