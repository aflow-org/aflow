// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu

#include "aflow_nhull_stats.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_xhost.h"
#include "modules/HULL/aflow_nhull_entry.h"

#define NHULL_STATS_DEBUG false

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::ios_base;
using std::istream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;
using std::vector;

namespace nhull {
  //***************** STATISTICS FUNCTIONS

  ///@brief takes a vector of pair<energy, degeneracy> values and adds up all degeneracies
  ///@param v vector of pair<energy, degeneracy>
  ///@return double of the sum of all degeneracies
  /// @authors
  /// mod{NHA,20260206,created}
  double degenSize(vector<std::pair<double, int>>& v) {
    double summation = 0;
    for (auto element : v) {
      summation += element.second;
    }
    return summation;
  }

  ///@brief takes a vector of pair<energy, degeneracy> and takes weighted mean of energies
  ///@param v vector of pair<energy, degeneracy>\
  /// @authors
  /// mod{NHA,20260206,created}
  double weightedMean(vector<std::pair<double, int>>& v) {
    if (v.empty()) {
      stringstream message;
      message << "ERROR IN FUNCTION median(): empty input!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
    }
    double total_elements = 0;
    double total_energy = 0;
    for (auto element : v) {
      total_elements += element.second;
      total_energy += element.first * element.second;
    }

    return total_energy / total_elements;
  }

  ///@brief takes a vector of pair<energy, degeneracy> values and computes the variance with respect to compound degeneracies
  ///@param v vector of pair<energy, degeneracy>
  ///@note used in DEED and EFA
  ///@return double of the sum of all degeneracies
  /// @authors
  /// mod{NHA,20260206,created}
  double degenVariance(vector<std::pair<double, int>>& v) {
    int n = degenSize(v);
    double squared_sum = 0;
    double variance = 0;
    double this_mean = weightedMean(v);

    if (n < 2) {
      stringstream message;
      message << "ERROR IN FUNCTION variance(): input too small!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
    }

    for (auto element : v) {
      double dif = element.first - this_mean;
      squared_sum += element.second * pow(dif, 2);
    }

    variance = squared_sum / (n - 1);

    return variance;
  }

  /// @brief statistics constructor for POCC entries. Degeneracies OF ENTRIES HAVE BEEN INITIALIZED BY THIS POINT
  /// @param members Entries that pertain to this POCC folder
  /// @authors
  /// @mod{NHA,20251005,created}
  compoundStatStore::compoundStatStore(const vector<Entry>& members) : members(members) {
    for (const auto& entry : members) {
      std::pair<double, int> enthalpy_store = {entry.enthalpy_formation_atom, entry.pocc_degeneracy};
      std::pair<double, int> dhull_store = {entry.distance_hull_enthalpy_formation_atom, entry.pocc_degeneracy};
      this->enthalpy_degeneracy_store.push_back(enthalpy_store);
      this->delta_hull_degeneracy_store.push_back(dhull_store);
    }
  }

  /// @brief function used to compile general, non-POCC statistics
  /// @authors
  /// mod{NHA,20260206,created}
  void compoundStatStore::compileStats() {
    for (const auto& element : this->members) {
      vector<double> tmp_enthalpies;
      vector<double> tmp_dhulls;
      if (element.flags.extraneous_energy_entry == false) {
        tmp_enthalpies.push_back(element.enthalpy_formation_atom);
        tmp_dhulls.push_back(element.distance_hull_enthalpy_formation_atom);
        total_members++;
      }
      this->enthalpies = aurostd::vector2xvector(tmp_enthalpies);
      this->delta_hulls = aurostd::vector2xvector(tmp_dhulls);
    }
    if (total_members == 0) {
      compute_statistics = false;
    }

    if (compute_statistics) {
      try {
        this->enthalpy_mean = aurostd::mean(this->enthalpies);
        this->enthalpy_median = aurostd::median(this->enthalpies);
        this->enthalpy_variance = aurostd::var(this->enthalpies);
        this->enthalpy_stdDev = sqrt(this->enthalpy_variance);
        this->delta_hull_mean = mean(this->delta_hulls);
        this->delta_hull_median = aurostd::median(this->delta_hulls);
        this->delta_hull_variance = aurostd::var(this->delta_hulls);
        this->delta_hull_stdDev = sqrt(this->delta_hull_variance);
        this->EFA = 1 / (this->enthalpy_stdDev);
        this->DEED = sqrt(this->EFA / (this->delta_hull_mean));
      } catch (aurostd::xerror& e) {
        if (XHOST.DEBUG || NHULL_STATS_DEBUG) {
          cout << "Non-fatal error in aflow_nhull_stats: either empty or out-of-bounds input in statistic functions." << endl;
        }
      }
    }
  }

  ///@brief function compiles statistics for a POCC directory
  /// @authors
  /// mod{NHA,20260206,created}
  void compoundStatStore::poccStats() {
    try {
      this->enthalpy_mean = weightedMean(this->enthalpy_degeneracy_store);
      this->enthalpy_variance = degenVariance(this->enthalpy_degeneracy_store);
      this->enthalpy_stdDev = sqrt(this->enthalpy_variance);
      this->delta_hull_mean = weightedMean(this->delta_hull_degeneracy_store);
      this->delta_hull_variance = degenVariance(this->delta_hull_degeneracy_store);
      this->delta_hull_stdDev = sqrt(this->delta_hull_variance);
      this->EFA = 1 / (this->enthalpy_stdDev);
      this->DEED = sqrt(this->EFA / (this->delta_hull_mean));
    } catch (aurostd::xerror& e) {
      if (XHOST.DEBUG || NHULL_STATS_DEBUG) {
        cout << "Non-fatal error in aflow_nhull_stats: either empty or out-of-bounds input in statistic functions." << endl;
      }
    }
  }

  /// @brief statistics done for single pocc compound using degeneracy.
  /// @param entries
  /// @param pocc_store list of degeneracies associated with a POCC compound
  /// @authors
  /// @mod{NHA,20251005,created}
  void singlePoccStats(vector<Entry>& entries, vector<pair<string, unsigned long long int>> pocc_store) {
    // updates degeneracy values of relevant entries
    // each pocc represented as aflowID (string), Degeneracy (int) value pair
    vector<Entry> members;
    for (Entry& entry : entries) {
      for (auto pocc : pocc_store) {
        if (pocc.first == entry.auid) {
          entry.set_pocc_degeneracy(pocc.second);
          members.push_back(entry);
        }
      }
    }

    // create new stats object and do stats:
    compoundStatStore pocc_stats(members);
    pocc_stats.poccStats();

    cout << endl;
    cout << "STATISTICS FOR POCC COMPOUND: " << endl;
    cout << "EFA = " << pocc_stats.EFA << endl;
    cout << "DEED = " << pocc_stats.DEED << endl << endl;
  }

  /// @brief function used in compound stats for all entries in run
  /// @param entries
  /// @param compounds vector of data structures for statistics of 1 compound; includes EFA and DEED
  /// @authors
  /// @mod {NHA,20251005,created}
  void statCalcLoop(vector<Entry>& entries, vector<compoundStatStore>& compounds) {
    vector<Entry> remaining_entries;
    compoundStatStore curr_cmpd_store;

    curr_cmpd_store.compound_map = entries[0].compound_map;
    for (auto& entry : entries) {
      if (entry.compound_map == curr_cmpd_store.compound_map) {
        curr_cmpd_store.members.push_back(entry);
      } else {
        remaining_entries.push_back(entry);
      }
    }
    // compile statistics for this compound
    curr_cmpd_store.compileStats();
    compounds.push_back(curr_cmpd_store);

    // reset for next loop
    entries = remaining_entries;
  }

  ///@brief compiles and outputs enthalpy statistics for each compound
  ///@param root_path path to save the statistics file
  /// @authors
  /// mod{NHA,20260206,created}
  void compoundStats(vector<Entry> entries, const string& root_path) {
    vector<compoundStatStore> compounds;

    while (!entries.empty()) {
      statCalcLoop(entries, compounds);
    }
    string save_path = root_path + "_stats.json";
    aurostd::JSON::object(compounds).saveFile(save_path);
  }
} // namespace nhull
