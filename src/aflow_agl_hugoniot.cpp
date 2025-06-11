// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *               Aflow PATRICK AVERY - Duke University 2016-2018           *
// *                                                                         *
// ***************************************************************************
// Written by Patrick Avery
// psavery@buffalo.edu
#ifndef _AFLOW_AGL_HUGONIOT_CPP
#define _AFLOW_AGL_HUGONIOT_CPP

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <climits>

#include "aflow.h"
#include "aflow_agl_debye.h"

namespace AGL_functions {
  // CHANGES BY PA20171224
  //
  // This function finds energy and mass density so that the following Hugoniot
  // expression is true: E - E0 = 0.5 * (P + P0)*(1.0/rho0 - 1.0/rho). If the
  // values that satisfy this expression are found between given data points,
  // interpolation will be performed to get the correct value.
  //
  // The initial conditions will be extracted automatically from the data.
  // The initial conditions are usually set to match experimental conditions,
  // and they are by default set to be 0 GPa and 300 K. At the start of the
  // function, the pressures will be investigated to see if there are any
  // that match initial conditions (the closest one will be chosen - a
  // warning will be printed if the closest one is more than 1 GPa away).
  // Then, the temperatures in that middle vector will be investigated to
  // see if there are any that match initial conditions (once again, the
  // closest one will be chosen - a warning will be printed if the closest
  // one is more than 10 K away). If the caller does not want to use the
  // default pressure of 0 GPa or the default temperature of 300 K, they can
  // change it by changing the parameters @p initial_pressure_external and
  // @p initial_temperature_external.
  // END CHANGES BY PA20171224
  //
  // Input for this function is organized as follows:
  // @param pressures_external The different pressures we will be using in GPa.
  // @param temperatures_external The different temperatures we will be using
  //                              in K.
  // @param mass_densities_gcm3 The different mass densities in grams per
  //                            cubic centimeter. The outer vector should be
  //                            varying the temperature (using the exact same
  //                            temperatures in the exact same order as in the
  //                            temperatures_external parameter), and the inner
  //                            vector should be varying the pressures (using the
  //                            exact same pressures in the exact same order as
  //                            in the pressures_external parameter).
  // @param energies_pT_kJg The different DFT energy + temperature-adjusted
  //                        internal energy values. The outer vector should be
  //                        (opposite to the mass_densities_gcm3 parameter)
  //                        varying the pressure (using the exact same pressures
  //                        in the exact same order as in the pressures_external
  //                        parameter), and the inner vector should be varying
  //                        the temperature (using the exact same temperatures in
  //                        the exact same order as in the temperatures_external
  //                        parameter).
  // @param desired_initial_pressure_external This is simply forwarded to
  //                                          calculateHugoniot(). See its
  //                                          description in the
  //                                          calculateHugoniot() section.
  // @param desired_initial_temperature_external This is simply forwarded to
  //                                             calculateHugoniot(). See its
  //                                             description in the
  //                                             calculateHugoniot() section.
  //
  // The output is set to @p results. These provide the calculated
  // solutions to the Hugoniot equation from the input data.
  // The output is organized as follows:
  //
  // inner vector: mass density (g/cc)
  //               external temperature (K)
  //               DFT energy + temperature-adjusted internal energy (kJ/g)
  //               pressure (GPa)
  //
  // outer vector: a series of inner vectors for the different pressures
  //               given in the input.
  //
  //
  // CHANGES BY PA20170603
  // @return 0 on success and 1 on critical failure. Returns 2 if one of the
  //         data points was found to extrapolate too far. All of the data
  //         stored in the output variables will be valid, but not every
  //         requested data point will be there. As soon as an invalid
  //         data point is found, the function returns. So there won't be
  //         any "gaps" in the valid data.
  // END CHANGES BY PA20170603
  //
  uint runHugoniot(const std::vector<double>& pressures_external,
      const std::vector<double>& temperatures_external,
      const std::vector<std::vector<double> >& mass_densities_gcm3,
      const std::vector<std::vector<double> >& energies_pT_kJg,
      std::vector<std::vector<double> >& hugoniotData, bool hugoniotextrapolate, ofstream& FileMESSAGE,
      double desired_initial_pressure_external,
      double desired_initial_temperature_external) {
    // bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    if(pressures_external.size() != energies_pT_kJg.size()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + " Error in " << __FUNCTION__ << ": pressures_external and "
        << "energies_pT_kJg should be the same size, but they are not!" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": pressures_external and "
      // OBSOLETE    << "energies_pT_kJg should be the same size, but they are "
      // OBSOLETE    << "not!\n";
      return 1;
    }
    if(temperatures_external.size() != mass_densities_gcm3.size()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + " Error in " << __FUNCTION__ << ": temperatures_external and "
        << "mass_densities_gcm3 should be the same size, but they are not!" <<  endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": temperatures_external and "
      // OBSOLETE           << "mass_densities_gcm3 should be the same size, but they are "
      // OBSOLETE           << "not!\n";
      return 1;
    }

    // First, find the best initial conditions
    double initial_mass_density_gcm3 = 0.0, initial_temperature_external = 0.0, initial_energyDFT_UIntVib_kJg = 0.0, initial_pressure_external = 0.0;

    // CHANGES BY PA20171224
    if(0 != findBestInitialConditions(pressures_external, temperatures_external, mass_densities_gcm3, energies_pT_kJg, desired_initial_pressure_external, desired_initial_temperature_external, initial_mass_density_gcm3, initial_temperature_external, initial_energyDFT_UIntVib_kJg, initial_pressure_external, FileMESSAGE)) {
      return 1;
    }
    // END CHANGES BY PA20171224

    // Clear the output if anything was written to it
    hugoniotData.clear();

    // CHANGES BY PA20171224
    // Now let's get all the data
    for(size_t i = 0; i < pressures_external.size(); ++i) {
      double pressure_external = pressures_external.at(i);
      // Create the input vectors
      const std::vector<double>& energiesDFT_UIntVib_kJg = energies_pT_kJg.at(i);
      std::vector<double> tmp_mass_densities_gcm3;
      for(size_t j = 0; j < temperatures_external.size(); ++j)
        tmp_mass_densities_gcm3.push_back(mass_densities_gcm3.at(j).at(i));
      // END CHANGES BY PA20171224

      // Now calculate the Hugoniot!
      double hugoniot_density = 0.0, hugoniot_temperature = 0.0, hugoniot_energy = 0.0;
      // CHANGES BY PA20170603
      uint ret = calculateHugoniotDataPoint(initial_mass_density_gcm3,
          initial_energyDFT_UIntVib_kJg,
          initial_pressure_external,
          // OBSOLETE mass_densities_gcm3,
          tmp_mass_densities_gcm3,
          temperatures_external,
          energiesDFT_UIntVib_kJg,
          pressure_external,
          hugoniot_density,
          hugoniot_temperature,
          hugoniot_energy, hugoniotextrapolate, FileMESSAGE);
      if(ret == 1) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + " calculateHugoniotDataPoint returned with an error!" <<  endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": calculateHugoniotDataPoint "
        // OBSOLETE          << "returned with an error!\n";	
        return ret;
      }
      else if(ret == 2) {
        // This means that we are extrapolating too far and should stop
        // collecting data. Just return. All data so far is valid.
        return ret;
      }
      // END CHANGES BY PA20170603

      // Store the results
      std::vector<double> results;
      results.clear();
      results.push_back(hugoniot_density);
      results.push_back(hugoniot_temperature);
      results.push_back(hugoniot_energy);
      results.push_back(pressure_external);
      hugoniotData.push_back(results);
    }

    // Success!
    return 0;
  }
} // namespace AGL_functions

namespace AGL_functions {
  // Convenience function to check if a double is an element of a vector of doubles
  bool inVector(double item, const std::vector<double>& vec, double tol = 1.e-5) {
    for(size_t i = 0; i < vec.size(); ++i) {
      if(aurostd::abs(item - vec.at(i)) < tol)
        return true;
    }
    return false;
  }
}

namespace AGL_functions {
  //
  // Same as runHugoniot() except it runs over all combinations of
  // temperatures and pressures. It needs the cell mass in grams as an
  // additional input.
  //
  // It probably takes longer to run this function than the above, but it has
  // more flexibility in terms of possible combinations of temperatures and
  // pressures.
  //
  // The definition of struct 'AGL_pressure_temperature_energies' can be
  // found in aflow_agl_debye.h
  //
  uint runHugoniotAllTemperaturesAndPressures(
      const std::vector<AGL_pressure_temperature_energies>& AGL_pressure_temperature_energy_list,
      double cellmass_grams,
      std::vector<std::vector<double> >& hugoniotData, bool hugoniotextrapolate,
      ofstream& FileMESSAGE,
      double desired_initial_pressure_external,
      double desired_initial_temperature_external) {
    ostringstream aus;
    // This tolerance will be used when comparing pressure values
    double pressureTol = 1.e-5;

    // Indices of these two vectors will correspond to indices
    // in AGL_pressure_temperature_energy_list
    std::vector<double> energies_pT_kJg;
    std::vector<double> mass_densities_gcm3;
    energies_pT_kJg.reserve(AGL_pressure_temperature_energy_list.size());
    mass_densities_gcm3.reserve(AGL_pressure_temperature_energy_list.size());

    // Cache all the different pressures that we have
    std::vector<double> pressures_external;
    pressures_external.reserve(100);
    for(size_t i = 0; i < AGL_pressure_temperature_energy_list.size(); ++i) {
      // Convert energy from eV/cell to kJ/g
      double energy_kJg = (AGL_pressure_temperature_energy_list.at(i).E_DFT_internal_vib_energy * 1.6e-22) / cellmass_grams;
      energies_pT_kJg.push_back(energy_kJg);

      // Next get equilibrium cell volume in units of cm^3
      double volume_cm3 = AGL_pressure_temperature_energy_list.at(i).volume_equilibrium * 1.e-24;
      // Calculate mass density and save for each temperature and pressure
      double massdensity_gcm3 = cellmass_grams / volume_cm3;
      mass_densities_gcm3.push_back(massdensity_gcm3);

      // Add the pressure if it isn't already in the pressures vector
      double pressure = AGL_pressure_temperature_energy_list.at(i).pressure_external;
      if(!inVector(pressure, pressures_external, pressureTol))
        pressures_external.push_back(pressure);
    }

    // First, find the best initial conditions
    double initial_mass_density_gcm3 = 0.0,
           initial_temperature_external = 0.0,
           initial_energyDFT_UIntVib_kJg = 0.0,
           initial_pressure_external = 0.0;

    if(0 != findBestInitialConditionsAllTemperaturesAndPressures(
          AGL_pressure_temperature_energy_list,
          mass_densities_gcm3, energies_pT_kJg,
          desired_initial_pressure_external,
          desired_initial_temperature_external,
          initial_mass_density_gcm3,
          initial_temperature_external,
          initial_energyDFT_UIntVib_kJg,
          initial_pressure_external, FileMESSAGE))
      return 1;

    // Clear the output if anything was written to it
    hugoniotData.clear();

    // Now let's get all the data
    for(size_t i = 0; i < pressures_external.size(); ++i) {
      double currentPressure = pressures_external.at(i);

      // Create the input vectors
      std::vector<double> temperatures_external;
      std::vector<double> energiesDFT_UIntVib_kJg;
      std::vector<double> tmp_mass_densities_gcm3;

      // '300' is probably a typical size
      temperatures_external.reserve(300);
      energiesDFT_UIntVib_kJg.reserve(300);
      tmp_mass_densities_gcm3.reserve(300);
      for(size_t j = 0; j < AGL_pressure_temperature_energy_list.size(); ++j) {
        double pressure = AGL_pressure_temperature_energy_list.at(j).pressure_external;
        if(aurostd::abs(pressure - currentPressure) < pressureTol) {
          temperatures_external.push_back(AGL_pressure_temperature_energy_list.at(j).temperature_external);
          energiesDFT_UIntVib_kJg.push_back(energies_pT_kJg.at(j));
          tmp_mass_densities_gcm3.push_back(mass_densities_gcm3.at(j));
        }
      }

      // Now calculate the Hugoniot!
      double hugoniot_density = 0.0, hugoniot_temperature = 0.0, hugoniot_energy = 0.0;
      uint ret = calculateHugoniotDataPoint(initial_mass_density_gcm3,
          initial_energyDFT_UIntVib_kJg,
          initial_pressure_external,
          tmp_mass_densities_gcm3,
          temperatures_external,
          energiesDFT_UIntVib_kJg,
          currentPressure,
          hugoniot_density,
          hugoniot_temperature,
          hugoniot_energy, hugoniotextrapolate, FileMESSAGE);
      if(ret == 1) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__
          << ": calculateHugoniotDataPoint returned with an error!" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return ret;
      }
      else if(ret == 2) {
        // This means that we are extrapolating too far and should stop
        // collecting data. Just return. All data so far is valid.
        return ret;
      }

      // Store the results
      std::vector<double> results;
      results.push_back(hugoniot_density);
      results.push_back(hugoniot_temperature);
      results.push_back(hugoniot_energy);
      results.push_back(currentPressure);
      hugoniotData.push_back(results);
    }

    // Success!
    return 0;
  }
} // namespace AGL_functions

// OBSOLETE namespace AGL_functions {
//
// This function finds energy and mass density so that the following Hugoniot
// expression is true: E - E0 = 0.5// (P + P0)*(1.0/rho0 - 1.0/rho). If the
// values that satisfy this expression are found between given data points,
// interpolation will be performed to get the correct value. The @p data
// vector should be organized as follows:
//
// inner vector: mass density (g/cc)
//               external temperature (K)
//               DFT energy + temperature-adjusted internal energy (kJ/g)
//               pressure (GPa)
//
// middle vector: a series of inner vectors each with a different
//                temperature. Usually the temperature is in series like:
//                0 K, 10 K, 20 K, 30 K.... The pressure should be constant
//                in every middle vector. It is ideal for the Hugoniot
//                solution to lie between temperatures. Otherwise,
//                extrapolation will have to be performed. Also, if the
//                spacing between the temperature values is smaller, the
//                solution will likely be more accurate.
//
// outer vector: a series of middle vectors each with a different
//               pressure.
//
// The initial conditions will be extracted automatically from the data.
// The initial conditions are usually set to match experimental conditions,
// and they are by default set to be 0 GPa and 300 K. At the start of the
// function, the pressures will be investigated to see if there are any
// that match initial conditions (the closest one will be chosen - a
// warning will be printed if the closest one is more than 1 GPa away).
// Then, the temperatures in that middle vector will be investigated to
// see if there are any that match initial conditions (once again, the
// closest one will be chosen - a warning will be printed if the closest
// one is more than 10 K away). If the caller does not want to use the
// default pressure of 0 GPa or the default temperature of 300 K, they can
// change it by changing the parameters @p initial_pressure_external and
// @p initial_temperature_external.
//
// The output is set to @p hugoniotData. These provide the calculated
// solutions to the Hugoniot equation from the input data.
// The output is organized as follows:
//
// inner vector: mass density (g/cc)
//               external temperature (K)
//               DFT energy + temperature-adjusted internal energy (kJ/g)
//               pressure (GPa)
//
// outer vector: a series of inner vectors for the different pressures
//               given in the input.
//
//
// @param data The input data. It is described in the text above.
// @param hugoniotData The output data. it is described in the text above.
// @param initial_pressure_external The pressure at the initial conditions.
//                                  Default is 0 GPa.
// @param initial_temperature_external The temperature at the initial
//                                     conditions. Default is 300 K.
//
// CHANGES BY PA20170603
// @return 0 on success and 1 on critical failure. Returns 2 if one of the
//         data points was found to extrapolate too far. All of the data
//         stored in the output variables will be valid, but not every
//         requested data point will be there. As soon as an invalid
//         data point is found, the function returns. So there won't be
//         any "gaps" in the valid data.
// END CHANGES BY PA20170603
//
// OBSOLETE uint calculateHugoniot(const std::vector<std::vector<std::vector<double> > >& data,
// OBSOLETE                      std::vector<std::vector<double> >& hugoniotData, ofstream& FileMESSAGE,
// OBSOLETE                      double desired_initial_pressure_external,
// OBSOLETE                      double desired_initial_temperature_external) {
// OBSOLETE  ostringstream aus;
// First, find the best initial conditions
// OBSOLETE  double initial_mass_density_gcm3 = 0.0,
// OBSOLETE         initial_temperature_external = 0.0,
// OBSOLETE         initial_energyDFT_UIntVib_kJg = 0.0,
// OBSOLETE         initial_pressure_external = 0.0;
// OBSOLETE
// OBSOLETE if(0 != findBestInitialConditions(data,
// OBSOLETE                                    desired_initial_pressure_external,
// OBSOLETE                                    desired_initial_temperature_external,
// OBSOLETE                                    initial_mass_density_gcm3,
// OBSOLETE                                    initial_temperature_external,
// OBSOLETE                                    initial_energyDFT_UIntVib_kJg,
// OBSOLETE                                     initial_pressure_external, FileMESSAGE))
// OBSOLETE   return 1;
// OBSOLETE
// Clear the output if anything was written to it
// OBSOLETE hugoniotData.clear();
// OBSOLETE
// Now let's get all the data
// OBSOLETE for(size_t i = 0; i < data.size(); ++i) {
// OBSOLETE  std::vector<double> mass_densities_gcm3, temperatures_external, energiesDFT_UIntVib_kJg;
// OBSOLETE  double pressure_external = -1.0;
// Create the input vectors
// OBSOLETE  for(size_t j = 0; j < data.at(i).size(); ++j) {
// OBSOLETE    mass_densities_gcm3.push_back(data.at(i).at(j).at(0));
// OBSOLETE    temperatures_external.push_back(data.at(i).at(j).at(1));
// OBSOLETE    energiesDFT_UIntVib_kJg.push_back(data.at(i).at(j).at(2));
// OBSOLETE    pressure_external = data.at(i).at(j).at(3);
// All of the pressures in this vector should be the same. Double check
// that and print an error if they are not.
// OBSOLETE    if(j != 0 && aurostd::abs(pressure_external - data.at(i).at(j - 1).at(3)) > 1e-5) {
// OBSOLETE	aurostd::StringstreamClean(aus);
// OBSOLETE	  aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": pressures should all be "
// OBSOLETE            << "the same in each middle vector! In middle vector number "
// OBSOLETE            << i << ", the last pressure was " << data.at(i).at(j - 1).at(3)
// OBSOLETE            << " GPa, and the current pressure is " << pressure_external
// OBSOLETE	      << " GPa!" << endl;
// OBSOLETE	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
// OBSOLETE      return 1;
// OBSOLETE    }
// OBSOLETE  }
// OBSOLETE
// Now calculate the Hugoniot!
// OBSOLETE  double hugoniot_density = 0.0, hugoniot_temperature = 0.0, hugoniot_energy = 0.0;
// CHANGES BY PA20170603
// OBSOLETE uint ret = calculateHugoniotDataPoint(initial_mass_density_gcm3,
// OBSOLETE                                      initial_energyDFT_UIntVib_kJg,
// OBSOLETE                                      initial_pressure_external,
// OBSOLETE                                      mass_densities_gcm3,
// OBSOLETE                                      temperatures_external,
// OBSOLETE                                      energiesDFT_UIntVib_kJg,
// OBSOLETE                                      pressure_external,
// OBSOLETE                                      hugoniot_density,
// OBSOLETE                                      hugoniot_temperature,
// OBSOLETE                                      hugoniot_energy, FileMESSAGE);
// OBSOLETE if(ret == 1) {
// OBSOLETE   aurostd::StringstreamClean(aus);
// OBSOLETE   aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": calculateHugoniotDataPoint returned with an error!" << endl;
// OBSOLETE   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
// OBSOLETE  return ret;
// OBSOLETE }
// OBSOLETE else if(ret == 2) {
// This means that we are extrapolating too far and should stop
// collecting data. Just return. All data so far is valid.
// OBSOLETE return ret;
// OBSOLETE }
// END CHANGES BY PA20170603

// Store the results
// OBSOLETE std::vector<double> results;
// OBSOLETE results.push_back(hugoniot_density);
// OBSOLETE results.push_back(hugoniot_temperature);
// OBSOLETE results.push_back(hugoniot_energy);
// OBSOLETE results.push_back(pressure_external);
// OBSOLETE hugoniotData.push_back(results);
// OBSOLETE }

// Success!
// OBSOLETE return 0;
// OBSOLETE }
// OBSOLETE } // namespace AGL_functions

namespace AGL_functions {
  //
  // This function finds energy and mass density so that the following Hugoniot
  // expression is true: E - E0 = 0.5 * (P + P0)*(1.0/rho0 - 1.0/rho). If the
  // values that satisfy this expression are found between data points given in
  // the density and energy vectors, interpolation will be performed to get the
  // correct value. The more points there are in the vector, the more accurate
  // the solution will be. It is very important that each index in each
  // vector correspond with the same index in the other vector. For example,
  // the mass density at index i should correspond to the temperature at
  // index i and the energy at index i. In addition, the vectors should
  // be sorted so that temperature is constantly increasing.
  //
  // The initial conditions are usually at 300K and at 0GPa (approximate
  // atmospheric conditions). The initial condition values have to be
  // extracted from the 0 GPa (or whatever the initial pressure is)
  // data and input as parameters here.
  //
  // @param initial_mass_density_gcm3 The initial conditions mass density of
  //                                  the material in g/cc.
  // @param initial_energyDFT_UIntVib_kJg The initial conditions energy of the
  //                                      material in kJ/g. This energy is
  //                                      usually the 0 K DFT energy +
  //                                      the vibrational internal energy
  //                                      at the initial conditions temperature.
  // @param initial_pressure_external The initial pressure of the system to
  //                                  be investigated (usually set to match
  //                                  experimental conditions) in GPa. This is
  //                                  typically 0 GPa.
  // @param mass_densities_gcm3 The mass densities (g/cc) corresponding to the
  //                            temperatures in @p temperatures_external and the
  //                            energies in @p energiesDFT_UIntVib_kJg.
  // @param temperatures_external The temperatures (K) corresponding to the
  //                              mass densities in @p mass_densities_gcm3 and
  //                              the energies in @p energiesDFT_UIntVib_kJg.
  // @param energiesDFT_UIntVib_kJg The energies (kJ/g) corresponding to the
  //                                temperatures in @p temperatures_external and
  //                                the mass densities in @p mass_densities_gcm3.
  //                                This energy is usually the 0 K DFT energy +
  //                                the vibrational internal energy at the
  //                                specified temperature.
  // @param pressure_external The constant external pressure (GPa) for this
  //                          dataset of densities, temperatures, and energies.
  // @param hugoniotDensity The output mass density (g/cc) solution to the
  //                        Hugoniot equation. This is an out variable - its
  //                        value will be overwritten if the function succeeds.
  // @param hugoniotTemperature The output temperature (K) solution to the
  //                            Hugoniot equation. This is an out variable - its
  //                            value will be overwritten if the function
  //                            succeeds.
  // @param hugoniotEnergy The output energy (kJ/g) solution to the
  //                       Hugoniot equation. This is an out variable - its
  //                       value will be overwritten if the function succeeds.
  //
  // CHANGES BY PA20170603
  // @return 0 on success. 1 on a critical failure. 2 if the data point is
  //         extrapolated too much too be trusted.
  // END CHANGES BY PA20170603
  //
  uint calculateHugoniotDataPoint(double initial_mass_density_gcm3,
      double initial_energyDFT_UIntVib_kJg,
      double initial_pressure_external,
      const std::vector<double>& mass_densities_gcm3,
      const std::vector<double>& temperatures_external,
      const std::vector<double>& energiesDFT_UIntVib_kJg,
      double pressure_external,
      double& hugoniot_density,
      double& hugoniot_temperature,
      double& hugoniot_energy, bool hugoniotextrapolate, ofstream& FileMESSAGE) {
    ostringstream aus;
    // We are finding E and P so that this will be equal to zero
    // E - E0 - 0.5 * (P + P0)(1.0/rho0 - 1.0/rho)

    // First do some initial checks
    if(mass_densities_gcm3.size() < 2) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error: " << __FUNCTION__ << " was called with too few points!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE std::cerr << "Error: " << __FUNCTION__
      // OBSOLETE           << " was called with too few points!\n";
      return 1;
    }

    // All of the sizes of the input vectors must be equal
    if(mass_densities_gcm3.size() != temperatures_external.size() ||
        mass_densities_gcm3.size() != energiesDFT_UIntVib_kJg.size()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + __FUNCTION__ << " was called, but the vectors do not have the same size!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE std::cerr << "Error: " << __FUNCTION__
      // OBSOLETE           << " was called, but the vectors do not have the same size!\n";
      return 1;
    }

    // It is VERY preferred that the solution be found in between two
    // temperatures in the TDE vector. Interpolation can be performed if
    // that is the case. If it can't be, it will print a warning to the
    // user. For now, try to find the two indices on either side of the
    // solution. ind2 will always be ind1 + 1 unless we have not yet found
    // the indices.
    size_t ind1 = 0, ind2 = 0;
    double lastVal = 0.0;
    for(size_t i = 0; i < mass_densities_gcm3.size(); ++i) {
      if(i != 0 && temperatures_external.at(i) < temperatures_external.at(i - 1)) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": Temperatures are not sorted in "
          << "increasing order! They must be sorted in increasing order to call this function." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // OBSOLETE  std::cerr << "Error in " << __FUNCTION__ << ": Temperatures are not sorted in "
        // OBSOLETE            << "increasing order! They must be sorted in increasing order "
        // OBSOLETE            << "to call this function.\n";
        return 1;
      }
      // We are trying to zero this expression.
      double val = energiesDFT_UIntVib_kJg.at(i) - initial_energyDFT_UIntVib_kJg - 0.5 *
        (pressure_external + initial_pressure_external) * (1.0/initial_mass_density_gcm3 - 1.0/mass_densities_gcm3.at(i));
      //std::cout << "E is " << energiesDFT_UIntVib_kJg.at(i) << "; D is " << mass_densities_gcm3.at(i) << "\n";
      // Figure out when we switch signs. The solution is between these
      // two points.
      if((val == 0.0 || aurostd::abs(val + lastVal) < aurostd::abs(val) + aurostd::abs(lastVal)) &&
          i != 0) {
        ind1 = i - 1;
        ind2 = i;
        break;
      }
      lastVal = val;
    }

    // If a solution was not found, find which end is closest to the
    // solution, and just use those two points. This SHOULD NOT
    // happen often, and it may result in an inaccurate solution.
    if(ind1 == 0 && ind2 == 0) {
      // Perform Hugoniot extrapolation if requested
      if (hugoniotextrapolate) {
        size_t lastInd = mass_densities_gcm3.size() - 1;
        // If the first point is closer to the solution, use the first two
        // points. Otherwise, use the last two points.
        if(aurostd::abs(energiesDFT_UIntVib_kJg.at(0) - initial_energyDFT_UIntVib_kJg -
              0.5 * (pressure_external + initial_pressure_external) * (1.0/initial_mass_density_gcm3 - 1.0/mass_densities_gcm3.at(0))) <
            aurostd::abs(energiesDFT_UIntVib_kJg.at(lastInd) - initial_energyDFT_UIntVib_kJg - 0.5 *
              (pressure_external + initial_pressure_external) * (1.0/initial_mass_density_gcm3 - 1.0/mass_densities_gcm3.at(lastInd)))) {
          ind1 = 0;
          ind2 = 1;
        }
        else {
          ind1 = lastInd - 1;
          ind2 = lastInd;
        }

        // Print a warning - the data may be bad.
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "WARNING: in " << __FUNCTION__ << ", the solution was found to be outside the temperature range." << endl;
        aus << _AGLSTR_WARNING_ + "Data should be used cautiously. Extrapolation from the two nearest points was used." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // OBSOLETE  std::cerr << "WARNING: in " << __FUNCTION__ << ", the solution was found "
        // OBSOLETE            << "to be outside the temperature range. Data should be used "
        // OBSOLETE            << "cautiously. Extrapolation from the two nearest points "
        // OBSOLETE            << "was used.\n";
      } else {
        // Exit loop for this pressure 
        return 2;
      }
    }

    // We may replace this function in the future if we decide it is necessary
    // to use a different function to perform better fitting.
    interpolateToFindHugoniotEnergyAndDensity(initial_mass_density_gcm3,
        initial_energyDFT_UIntVib_kJg,
        initial_pressure_external,
        pressure_external,
        energiesDFT_UIntVib_kJg.at(ind1),
        energiesDFT_UIntVib_kJg.at(ind2),
        mass_densities_gcm3.at(ind1),
        mass_densities_gcm3.at(ind2),
        hugoniot_energy,
        hugoniot_density);

    // To find the temperature, first find what fraction of a distance we are
    // from the ind1 energy. This fraction should also correspond to the distance
    // we are from the ind1 temperature. This should always be positive unless
    // the solution is less than temperatures_external[0]. Then (and only then) should hugoniot_energy
    // be less than energiesDFT_UIntVib_kJg.at(ind1). In addition, fraction may be greater than 1.0
    // if and only if the solution is greater than temperatures_external.at(lastInd - 1).
    double fraction = (hugoniot_energy - energiesDFT_UIntVib_kJg.at(ind1)) / aurostd::abs(energiesDFT_UIntVib_kJg.at(ind2) - energiesDFT_UIntVib_kJg.at(ind1));

    // The temperature should be the same fraction away from the ind1
    // temperature.
    hugoniot_temperature = aurostd::abs(temperatures_external.at(ind1) - temperatures_external.at(ind2)) * fraction +
      temperatures_external.at(ind1);

    // CHANGES BY PA20170603
    // If this temperature is more than 1.5 times greater than the nearest temperature
    // in the temperatures_external vector, we should assume that the data is invalid and
    // return 2. The answer will still be stored one function up for further analysis,
    // though...
    if(hugoniot_temperature > 1.5 * temperatures_external.back()) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "In " << __FUNCTION__ << ", the solution was found to be greater than 1.5 times the greatest temperature value in the temperatures_external vector." << endl;
      aus << _AGLSTR_WARNING_ + "Data should beyond this temperature should be considered invalid since too great an extrapolation was used." << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE  std::cerr << "Error: in " << __FUNCTION__ << ", the solution was found "
      // OBSOLETE            << "to be greater than 1.5 times the greatest temperature "
      // OBSOLETE            << "value in the temperatures_external vector. Data should "
      // OBSOLETE            << "be considered invalid since too great an extrapolation "
      // OBSOLETE            << "was used.";
      return 2;
    }
    // END CHANGES BY PA20170603

    return 0;
  }
} // namespace AGL_functions

namespace AGL_functions {
  //
  // NOT TO CALLED (only to be used by calculateHugoniot()).
  //
  // Looks through the data given to calculateHugoniot() and tries to
  // find data matching the desired_initial_pressure_external and the
  // desired_initial_temperature_external. It then sets the initial
  // conditions and returns 0.
  //
  // CHANGES BY PA20171224  
  uint findBestInitialConditions(const std::vector<double>& pressures_external,
      const std::vector<double>& temperatures_external,
      const std::vector<std::vector<double> >& mass_densities_gcm3,
      const std::vector<std::vector<double> >& energies_pT_kJg,
      double desired_initial_pressure_external,
      double desired_initial_temperature_external,
      double& initial_mass_density_gcm3,
      double& initial_temperature_external,
      double& initial_energyDFT_UIntVib_kJg,
      double& initial_pressure_external, ofstream& FileMESSAGE) {
    ostringstream aus;
    // OBSOLETE if(data.size() == 0) {
    // OBSOLETE   aurostd::StringstreamClean(aus);
    // OBSOLETE   aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": data is empty!" << endl;
    // OBSOLETE   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    // OBSOLETE   return 1;
    // OBSOLETE }

    // Find a matching pressure first. Then find a matching temperature in that pressure.
    int bestInitialPressureInd = -1, bestInitialTemperatureInd = -1;
    double bestInitialPressure = -1.0, bestInitialTemperature = -1.0;
    // OBSOLETE bestInitialEnergy = 0.0, bestInitialDensity = -1.0;
    // OBSOLETE for(size_t i = 0; i < data.size(); ++i)
    // Must not be empty. Skip over it if it is.
    // OBSOLETE if(data.at(i).size() == 0) {
    // OBSOLETE aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_WARNING_ + "Warning in " << __FUNCTION__ << ": data at index " << i << " is empty! Skipping over it." << endl;
    // OBSOLETE aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OBSOLETE std::cerr << "Warning in " << __FUNCTION__ << ": data at index "
    // OBSOLETE           << i << " is empty! Skipping over it\n";
    // OBSOLETE continue;
    // OBSOLETE }
    // The inner vector must be of size 4. Return with an error if it isn't.
    // OBSOLETE if(data.at(i).at(0).size() != 4) {
    // OBSOLETE aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": data.at("<< i << ").at(0) should have a size of 4, "
    // OBSOLETE    << "but instead it has a size of " << data.at(i).at(0).size() << "! Data is invalid." << endl;
    // OBSOLETE aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);     
    // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": data.at("
    // OBSOLETE           << i << ").at(0) should have a size of 4, but instead it "
    // OBSOLETE           << "has a size of " << data.at(i).at(0).size() << "! Data "
    // OBSOLETE           << "is invalid.\n";
    // OBSOLETE return 1;
    // OBSOLETE }
    // Otherwise, check the pressure. If it hasn't been set or it is
    // closer to the desired initial, keep it.
    // OBSOLETE double currentPressure = data.at(i).at(0).at(3);
    // CHANGES BY PA20171224
    for(size_t i = 0; i < pressures_external.size(); ++i) {
      // Check the pressure. If it hasn't been set or it is
      // closer to the desired initial, keep it.
      double currentPressure = pressures_external.at(i);
      // END CHANGES BY PA20171224
      if(bestInitialPressureInd == -1 || aurostd::abs(currentPressure - desired_initial_pressure_external) < aurostd::abs(bestInitialPressure - desired_initial_pressure_external)) {
        bestInitialPressureInd = i;
        bestInitialPressure = currentPressure;
      }
    }

    // If a best initial pressure wasn't found, return false.
    if(bestInitialPressureInd == -1 || bestInitialPressure < 0.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": best initial pressure index was not found!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);           
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": best initial pressure "
      // OBSOLETE          << "index was not found!\n";
      return 1;
    }

    // If the best initial pressure is greater than a difference of
    // 1.0 GPa, print a warning.
    if(aurostd::abs(bestInitialPressure - desired_initial_pressure_external) > 1.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Warning in " << __FUNCTION__ << ": the closest pressure to the desired initial pressure of "
        << desired_initial_pressure_external << " GPa was found to be " << bestInitialPressure << " GPa." << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);           
      // OBSOLETE std::cerr << "Warning in " << __FUNCTION__ << ": the closest pressure "
      // OBSOLETE          << "to the desired initial pressure of "
      // OBSOLETE          << desired_initial_pressure_external << " GPa was found to be "
      // OBSOLETE          << bestInitialPressure << " GPa.\n";
    }

    // Now let's find a best temperature at this pressure.
    // OBSOLETE for(size_t i = 0; i < data.at(bestInitialPressureInd).size(); ++i)
    // The inner vector must be of size 4 again. Return an error if it isn't.
    // OBSOLETE if(data.at(bestInitialPressureInd).at(i).size() != 4) {
    // OBSOLETE aurostd::StringstreamClean(aus);
    // OBSOLETE aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": data.at(" << bestInitialPressureInd << ").at(" << i << ") should have a size of 4,"
    // OBSOLETE   << " but instead it has a size of " << data.at(bestInitialPressureInd).at(i).size() << "! Data " << "is invalid." << endl;
    // OBSOLETE aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);           	
    // OBSOLETE  std::cerr << "Error in " << __FUNCTION__ << ": data.at("
    // OBSOLETE            << bestInitialPressureInd << ").at(" << i
    // OBSOLETE            << ") should have a size of 4, but instead it "
    // OBSOLETE            << "has a size of "
    // OBSOLETE            << data.at(bestInitialPressureInd).at(i).size()
    // OBSOLETE            << "! Data " << "is invalid.\n";
    // OBSOLETE return 1;
    // OBSOLETE }
    // Otherwise, check the temperature. If it hasn't been set or it is
    // closer to the desired initial, keep it.
    // OBSOLETE double currentTemperature = data.at(bestInitialPressureInd).at(i).at(1);
    // CHANGES BY PA20171224
    for(size_t i = 0; i < temperatures_external.size(); ++i) {
      // Check the temperature. If it hasn't been set or it is
      // closer to the desired initial, keep it.
      double currentTemperature = temperatures_external.at(i);
      // END CHANGES BY PA20171224   
      if(bestInitialTemperatureInd == -1 ||
          aurostd::abs(currentTemperature - desired_initial_temperature_external) <
          aurostd::abs(bestInitialTemperature - desired_initial_temperature_external)) {
        bestInitialTemperatureInd = i;
        bestInitialTemperature = currentTemperature;
      }
    }

    // If a best initial temperature wasn't found, return false.
    if(bestInitialTemperatureInd == -1 || bestInitialTemperature < 0.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": best initial temperature index was not found!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": best initial temperature "
      // OBSOLETE           << "index was not found!\n";
      return 1;
    }

    // If the best initial temperature is greater than a difference of
    // 10.0, print a warning.
    if(aurostd::abs(bestInitialTemperature - desired_initial_temperature_external) > 10.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Warning in " << __FUNCTION__ << ": the closest temperature to the desired initial temperature of "
        << desired_initial_temperature_external << " K was found to be " << bestInitialTemperature << " K." << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);       
      // OBSOLETE std::cerr << "Warning in " << __FUNCTION__ << ": the closest temperature "
      // OBSOLETE           << "to the desired initial temperature of "
      // OBSOLETE           << desired_initial_temperature_external << " K was found to be "
      // OBSOLETE           << bestInitialTemperature << " K.\n";
    }

    // Now that we have the indices, extract the data for the initial conditions.
    // OBSOLETE initial_energyDFT_UIntVib_kJg = data.at(bestInitialPressureInd).at(bestInitialTemperatureInd).at(2);
    initial_energyDFT_UIntVib_kJg = energies_pT_kJg.at(bestInitialPressureInd).at(bestInitialTemperatureInd);
    initial_temperature_external = bestInitialTemperature;
    // OBSOLETE initial_mass_density_gcm3 = data.at(bestInitialPressureInd).at(bestInitialTemperatureInd).at(0);
    initial_mass_density_gcm3 = mass_densities_gcm3.at(bestInitialTemperatureInd).at(bestInitialPressureInd);
    initial_pressure_external = bestInitialPressure;

    return 0;
  }
} // namespace AGL_functions

namespace AGL_functions {
  //
  // NOT TO CALLED (only to be used by calculateHugoniot()).
  //
  // Same as findBestInitialConditions() except it can be run for all
  // combinations of temperatures and pressures.
  //
  // The definition of the struct 'AGL_pressure_temperature_energies' can be
  // found in aflow_agl_debye.h
  //
  uint findBestInitialConditionsAllTemperaturesAndPressures(
      const std::vector<AGL_pressure_temperature_energies>& AGL_pressure_temperature_energy_list,
      const std::vector<double>& mass_densities_gcm3,
      const std::vector<double>& energies_pT_kJg,
      double desired_initial_pressure_external,
      double desired_initial_temperature_external,
      double& initial_mass_density_gcm3,
      double& initial_temperature_external,
      double& initial_energyDFT_UIntVib_kJg,
      double& initial_pressure_external, ofstream& FileMESSAGE) {
    ostringstream aus;
    // Find a matching pressure first. Then find a matching temperature that has that pressure.
    double bestInitialPressure = 1.e300;
    for(size_t i = 0; i < AGL_pressure_temperature_energy_list.size(); ++i) {
      const double& currentPressure = AGL_pressure_temperature_energy_list.at(i).pressure_external;
      if(aurostd::abs(currentPressure - desired_initial_pressure_external) <
          aurostd::abs(bestInitialPressure - desired_initial_pressure_external)) {
        bestInitialPressure = currentPressure;
      }
    }

    if(bestInitialPressure > 1.e299) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": best initial pressure was not found!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);           
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": best initial pressure "
      // OBSOLETE          << "was not found!\n";
      return 1;
    }

    double pressureTol = 1.e-5;
    double bestInitialTemperature = 1.e300;
    int bestConditionsInd = -1;
    for(size_t i = 0; i < AGL_pressure_temperature_energy_list.size(); ++i) {
      const double& currentPressure = AGL_pressure_temperature_energy_list.at(i).pressure_external;
      if(aurostd::abs(currentPressure - bestInitialPressure) < pressureTol) {
        const double& currentTemperature = AGL_pressure_temperature_energy_list.at(i).temperature_external;
        if(aurostd::abs(currentTemperature - desired_initial_temperature_external) <
            aurostd::abs(bestInitialTemperature - desired_initial_temperature_external)) {
          bestInitialTemperature = currentTemperature;
          bestConditionsInd = i;
        }
      }
    }

    if(bestInitialTemperature > 1.e299) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": best initial temperature was not found!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);           
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": best initial temperature "
      // OBSOLETE          << "was not found!\n";
      return 1;
    }

    if(bestConditionsInd == -1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Error in " << __FUNCTION__ << ": best initial conditions were not found!" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // OBSOLETE std::cerr << "Error in " << __FUNCTION__ << ": best initial conditions "
      // OBSOLETE          << "were not found!\n";
      return 1;
    }

    // If the best initial pressure is greater than a difference of
    // 1.0 GPa, print a warning.
    if(aurostd::abs(bestInitialPressure - desired_initial_pressure_external) > 1.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Warning in " << __FUNCTION__ << ": the closest pressure to the desired initial pressure of "
        << desired_initial_pressure_external << " GPa was found to be " << bestInitialPressure << " GPa" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);           
      // OBSOLETE std::cerr << "Warning in " << __FUNCTION__ << ": the closest pressure "
      // OBSOLETE           << "to the desired initial pressure of "
      // OBSOLETE           << desired_initial_pressure_external << " GPa was found to be "
      // OBSOLETE           << bestInitialPressure << " GPa.\n";
    }

    // If the best initial temperature is greater than a difference of
    // 10.0, print a warning.
    if(aurostd::abs(bestInitialTemperature - desired_initial_temperature_external) > 10.0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Warning in " << __FUNCTION__ << ": the closest temperature at pressure " << bestInitialPressure << " GPa to the desired initial temperature of "
        << desired_initial_temperature_external << " K was found to be " << bestInitialTemperature << " K." << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);                 
      // OBSOLETE std::cerr << "Warning in " << __FUNCTION__ << ": the closest temperature "
      // OBSOLETE          << "at pressure " << bestInitialPressure
      // OBSOLETE          << " to the desired initial temperature of "
      // OBSOLETE          << desired_initial_temperature_external << " K was found to be "
      // OBSOLETE          << bestInitialTemperature << " K.\n";
    }

    // Now that we have the indices, extract the data for the initial conditions.
    initial_energyDFT_UIntVib_kJg = energies_pT_kJg.at(bestConditionsInd);
    initial_temperature_external = bestInitialTemperature;
    initial_mass_density_gcm3 = mass_densities_gcm3.at(bestConditionsInd);
    initial_pressure_external = bestInitialPressure;

    return 0;
  }
} // namespace AGL_functions

namespace AGL_functions {
  //
  // NOT TO CALLED (only to be used by calculateHugoniotDataPoint()).
  //
  // Creates a density vs. energy line using two energy points and two
  // density points. It then interpolates using the equation of this line
  // in order to find a density and energy answer for the Hugoniot equation.
  //
  // If we find out in the future that an interpolation is not an adequate
  // approximation, we can use some other kind of fit in order to find the
  // answer.
  //
  // This function is used to calculate hugoniot_energy and hugoniot_density by
  // creating an equation for a line with two points. It then analyically
  // calculates the solution through interpolation with those two points.
  // It can and will perform extrapolation, but it is highly preferred that
  // the solution be between energy1 and energy2.
  // This function can be replaced in the future using some other fitting
  // functions if we find that fitting in this way is unsatisfactory.
  //
  void interpolateToFindHugoniotEnergyAndDensity(double initial_mass_density_gcm3,
      double initial_energyDFT_UIntVib_kJg,
      double initial_pressure_external,
      double pressure_external,
      double energy1,
      double energy2,
      double density1,
      double density2,
      double& hugoniot_energy,
      double& hugoniot_density) {
    // Find the equation of the interpolation line: D vs. E.
    // This only uses two points...
    double slope = (density2 - density1) / (energy2 - energy1);
    double slope_intercept = density1 - slope * energy1;

    // Now calculate the point at which the Hugoniot is satisfied.
    // This expression is cached only because it will be used several times
    double cached_expression = pressure_external + initial_pressure_external + 2.0 * initial_mass_density_gcm3 * initial_energyDFT_UIntVib_kJg;
    // The following equation is ugly, but it is an analytical solution that results from the
    // interpolation.
    hugoniot_energy = (sqrt(4.0 * pow(slope_intercept, 2.0) * pow(initial_mass_density_gcm3, 2.0) +
          4.0 * slope_intercept * slope * cached_expression * initial_mass_density_gcm3 +
          slope * (slope * pow(cached_expression, 2.0) -
            8.0 * (pressure_external + initial_pressure_external) * pow(initial_mass_density_gcm3, 2.0))) -
        2.0 * slope_intercept * initial_mass_density_gcm3 + slope * cached_expression) /
      (4.0 * slope * initial_mass_density_gcm3);

    // Plug this into our equation to find density
    hugoniot_density = slope * hugoniot_energy + slope_intercept;
  }
} // namespace AGL_functions

#endif // _AFLOW_AGL_HUGONIOT_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *               Aflow PATRICK AVERY - Duke University 2016-2018           *
// *                                                                         *
// ***************************************************************************
