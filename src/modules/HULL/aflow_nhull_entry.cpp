#include "aflow_nhull_entry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_nhull.h"
#include "aflow_nhull_util.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "extern/QHULL/libqhullcpp/RoadError.h"

using std::string;
using std::vector;

//********************** BEGIN CLASS ENTRY

namespace nhull {

  ///@brief constructor for thermodynamic hull calculations
  ///@param aflow_lib_entry
  ///@param  global_specie_vector
  ///@param nhull_index sets the integer index of this entry for ease of manipulation later on
  /// @authors
  /// mod{NHA,20260206,created}
  Entry::Entry(const aflowlib::_aflowlib_entry& aflow_lib_entry, const vector<string>& global_specie_vector, int nhull_index, bool neglect_cce) :
      compound(aflow_lib_entry.compound),
      vcomposition(aflow_lib_entry.vcomposition),
      vspecies(aflow_lib_entry.vspecies),
      global_vspecies(global_specie_vector),
      prototype(aflow_lib_entry.prototype),
      auid(aflow_lib_entry.auid),
      aurl(aflow_lib_entry.aurl),
      nspecies(aflow_lib_entry.nspecies),
      spin_atom(aflow_lib_entry.spin_atom),
      entropic_temperature(aflow_lib_entry.entropic_temperature),
      ldau_TLUJ(aflow_lib_entry.ldau_TLUJ),
      nh_index(nhull_index),
      aflow_lib_entry(aflow_lib_entry.entry),
      lib_entry_direct(aflow_lib_entry),
      // necessary constant initialization
      is_hull_point(false),
      lowest_energy_alloy(false),
      flagged_entry(false),
      m_is_artificial(false) {
    bool cce_used = false;
    enthalpy_formation_atom = H_f_atom(aflow_lib_entry, neglect_cce, cce_used);
    m_uses_cce = cce_used;
    this->compoundMapInit(); // generate compound map for this alloy (reduced stoichiometry)
    this->init_mCoord(global_specie_vector); //generate highest dimension qhull point for distance calcs
    if (this->nspecies == 1) {
      // sets pure elements to decompose into themselves
      this->nhull_phase_decomp = {
          {this->compound, 1}
      };
      //initialize all flags:
      this->flags = Flags();
    }
  }

  ///@brief lightweight constructor: instantiates point coordinates and whether the point is an artificial unary
  /// @authors
  /// mod{NHA,20260206,created}
  Entry::Entry(const aurostd::xvector<double>& coord, bool is_artificial) {
    m_coord = coord;
    nspecies = coord.rows;
    m_is_artificial = is_artificial;
    m_initialized = true;
  }

  ///@brief constructor used only for GFA
  ///@param has_stoich_coords legacy parameter
  ///@param  formation_energy_coord legacay parameter
  /// @authors
  /// mod{NHA,20260206,created}
  Entry::Entry(const aurostd::xvector<double>& coord, std::ostream& oss, bool has_stoich_coords, bool formation_energy_coord, bool is_artificial, const uint ch_index) {
    m_coord = coord;
    nspecies = coord.rows;
    nh_index = ch_index;
    m_is_artificial = is_artificial;
    m_initialized = true;
  }

  ///@brief constructor used to retrieve test data for unit tests
  /// @authors
  /// mod{NHA,20260206,created}
  Entry::Entry(const std::map<std::string, double>& nhull_phase_decomp, const std::string& auid, const double dhull) :
      nhull_phase_decomp(nhull_phase_decomp), auid(auid), distance_hull_enthalpy_formation_atom(dhull) {}

  ///@brief takes a list of point coordinates and creates a list of entries
  ///@return list of entries
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> Entry::pointsToEntries(const vector<aurostd::xvector<double>>& points) {
    vector<Entry> entries;
    for (const auto& point : points) {
      Entry new_entry(point, false);
      entries.push_back(new_entry);
    }

    return entries;
  }

  ///@brief given a point, finds the associated entry, updates found_entry with the found entry or returns false if it could not find one
  ///@note assumes point is non-stoich (first element index removed)
  ///@return true if function was able to find an entry in entries
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::getEntryFromPoint(const vector<Entry>& entries, Entry& found_entry, const aurostd::xvector<double>& point) {
    bool found_entry_flag = false;
    for (const Entry& entry : entries) {
      if (entry.entryEqualsPoint(point)) {
        found_entry = entry;
        found_entry_flag = true;
        break;
      }
    }
    return found_entry_flag;
  }

  ///@brief returns true if the current entry shares the same coordinates as point
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::entryEqualsPoint(const aurostd::xvector<double>& point) const {
    return aurostd::identical(this->m_coord, point, QHULL_POINT_TOL);
  }

  /// @brief returns true if a given coordinate represents a unary: assumes point is not a stoich point (first index removed and with energy included)
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::isPointUnary(const aurostd::xvector<double>& point) {
    uint non_zero_row_count = 0;
    //ignore energy coordinate
    for (int i = point.lrows; i < point.urows; i++) {
      if (fabs(point[i]) > QHULL_POINT_TOL) {
        non_zero_row_count++;
      }
    }
    if (non_zero_row_count == 1) {
      return true;
    }
    return false;
  }

  /// @brief converts an entry to a qhull point given a sorted specie_vector
  /// @param specie_vector
  /// @param full_point full_point = true returns a full stoich point (ratios of all species included), full_point = false returns a qhull point (first specie removed)
  /// @return returns a vector (point) in composition-enthalpy space
  /// @authors
  /// @mod{NHA,20251005,created}
  aurostd::xvector<double> Entry::compoundToPoint(const vector<string>& specie_vector, bool full_point) const {
    const int dimension = specie_vector.size(); // dimension of qhull system = (n-1)
    vector<double> stoich_point(dimension, 0);

    // use the following standard: (specie_2, specie_3, .... specie_n, energy)
    for (const auto& [symbol, multiplicity] : this->compound_map) {
      for (int i = 0; i < dimension; i++) {
        // This if statement ensures consistent indexing of stoichiometries among qhull points
        if (specie_vector[i] == symbol) {
          stoich_point[i] = multiplicity;
        }
      }
    }

    vector<double> return_point;
    if (full_point) {
      return_point = fullCoords(stoichCoord(stoich_point), false);
    } else {
      return_point = stoichCoord(stoich_point);
    }
    return_point.push_back(this->enthalpy_formation_atom); // add back enthalpy here

    return aurostd::vector2xvector(return_point);
  }

  /// @brief initializes the qhull point representation of this compound in the given specie vector
  /// @authors
  /// mod{NHA,20260206,created}
  void Entry::init_mCoord(vector<string> specie_vector) {
    //ensure specie vector is sorted
    std::sort(specie_vector.begin(), specie_vector.end());
    vector<double> mcoord;
    std::map<string, int> tmp_map;
    double total_comp = 0;

    for (const double i : vcomposition) {
      total_comp += i;
    }

    if (this->vcomposition.size() != this->vspecies.size()) {
      std::cerr << "ERROR IN FUNCTION compoundMap(): vcomposition and vspecies are not same dimension!" << std::endl;
    } else {
      for (size_t i = 0; i < this->vcomposition.size(); i++) {
        tmp_map[this->vspecies[i]] = this->vcomposition[i];
      }
    }

    //skip first element
    for (size_t i = 1; i < specie_vector.size(); i++) {
      auto it = tmp_map.find(specie_vector[i]);
      if (it != tmp_map.end()) {
        mcoord.push_back(it->second / total_comp);
      } else {
        mcoord.push_back(0);
      }
    }

    //add energy coordinate:
    mcoord.push_back(enthalpy_formation_atom);

    m_coord = aurostd::vector2xvector(mcoord);
  }

  ///@brief takes a list of entries and returns a list of their point representations
  /// @authors
  /// mod{NHA,20260206,created}
  vector<aurostd::xvector<double>> Entry::entriesToPoints(vector<Entry> entries) {
    vector<aurostd::xvector<double>> b;
    for (const Entry& entry : entries) {
      b.push_back(entry.m_coord);
    }
    return b;
  }

  /// @brief takes in compound string and returns sorted specie vector
  /// @authors
  /// mod{NHA,20260206,created}
  vector<string> Entry::toSpecieVector() const {
    string specie_string;

    // first convert current entry to specie vector
    for (const auto chr : this->compound) {
      if (isdigit(chr) == 0) {
        specie_string.push_back(chr);
      }
    }
    auto elements = aurostd::getElements(specie_string);
    std::sort(elements.begin(), elements.end());
    return elements;
  }

  ///@brief converts all energy quantities from eVs to meVs for this entry
  /// @authors
  /// mod{NHA,20260206,created}
  void Entry::enthalpyUnitConversion() {
    // converts an entry from eV to meV per convention
    this->enthalpy_formation_atom *= 1000;
    this->distance_hull_enthalpy_formation_atom *= 1000;
    this->m_nminus1_enthalpy_gain *= 1000;
    this->m_stability_criterion *= 1000;
  }

  /// @brief checks ldau_TLUJ is empty; if not +U was used for this entry and should be ignored
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::plusUCheck() {
    return !this->ldau_TLUJ.empty();
  }

  /// @brief checks whether given compound contains species not present in alloy_name_store
  /// @param alloy_name_store
  /// @return false if found element not in current entry and true otherwise
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::isCompoundCompatible(const vector<string>& alloy_name_store) {
    // now check if this compound contains elements not found in alloy_name store:

    for (const auto element : this->toSpecieVector()) {
      if (count(alloy_name_store.begin(), alloy_name_store.end(), element) == 0) {
        return false;
      }
    }

    return true;
  }

  ///@brief returns true if this entry contains only elements of alloy_index_store
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::isCompoundCompatible(const aurostd::xvector<int>& alloy_index_store) {
    // now check if this compound contains elements not found in alloy_name store:
    vector<int> alloy_conv = aurostd::xvector2vector(alloy_index_store);

    aurostd::xvector<int> elements_present = getElementsPresent();

    if (alloy_index_store.rows < elements_present.rows) {
      return false;
    }

    for (int i = alloy_index_store.lrows; i < alloy_index_store.urows; i++) {
      if (alloy_index_store[i] == 0 && elements_present[i] != 0) {
        return false;
      }
    }

    return true;
  }

  /// @brief finds and returns to the highest dimension found among alloys in entries vector
  /// @param entries vector of entries to analyze
  /// @return max_dim
  /// @authors
  /// mod{NHA,20260206,created}
  uint Entry::maxDimension(const vector<Entry>& entries) {
    int max_dim = 2;
    for (const auto& entry : entries) {
      max_dim = std::max(max_dim, entry.nspecies);
    }
    return max_dim;
  }

  ///@brief returns the lowest energy from a vector of entries
  /// @authors
  /// mod{NHA,20260206,created}
  double Entry::minEnergy(const vector<Entry>& entries) {
    double min_energy = 0;
    for (const auto& entry : entries) {
      min_energy = std::min(min_energy, entry.enthalpy_formation_atom);
    }
    return min_energy;
  }

  ///@brief return true if a given auid is in a list of entries
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::findAUID(const string auid, const vector<Entry>& entries) {
    for (const Entry& entry : entries) {
      if (entry.auid == auid) {
        return true;
      }
    }
    return false;
  }

  ///@brief sets the enthalpy_formation_atom of an entry and updates is_hull_point if it is close enough to hull
  /// @authors
  /// mod{NHA,20260206,created}
  void Entry::set_distance_hull_enthalpy_formation_atom(double distance_hull_enthalpy_formation_atom) {
    const double cutoff = NHULL_HULL_POINT_CUTOFF; // accounts for numerical error
    this->distance_hull_enthalpy_formation_atom = distance_hull_enthalpy_formation_atom;
    if (distance_hull_enthalpy_formation_atom < cutoff && this->nspecies != 1) {
      this->is_hull_point = true;
    }
  }

  void Entry::set_relative_stability_criterion(double distance_hull_enthalpy_formation_atom, double enthalpy_formation_atom) {
    this->relative_stability_criterion = std::abs(distance_hull_enthalpy_formation_atom / enthalpy_formation_atom);
  }

  ///@brief returns true if entries enthalpy_formation_atom lies outside of cutoff value
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::is_outside_IRQ(const double inner_quartile_cutoff) {
    bool outside_iqr = inner_quartile_cutoff != 0 && this->enthalpy_formation_atom < inner_quartile_cutoff;
    return outside_iqr;
  }

  /// @brief intitializes the map of compounds for an entry
  /// @authors
  /// mod{NHA,20260206,created}
  void Entry::compoundMapInit() {
    std::map<string, int> compound_map; // an individual compound

    if (this->vcomposition.size() != this->vspecies.size()) {
      std::cerr << "ERROR IN FUNCTION compoundMap(): vcomposition and vspecies are not same dimension!" << std::endl;
    } else {
      for (int i = 0; i < this->vcomposition.size(); i++) {
        compound_map[this->vspecies[i]] = this->vcomposition[i];
      }
    }

    // reduce stoichiometry of this compound
    compound_map = reducedStoich(compound_map);

    if (compound_map.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Stability criterion cannot be done: Chosen auid has no compound map! Check to make sure all entries have been properly initialized.");
    }

    this->compound_map = compound_map;
  }

  /// @brief converts a string of a compounds into a map
  /// @param compound
  /// @return compound_map
  /// @authors
  /// mod{NHA,20260206,created}
  std::map<string, int> Entry::compoundMap(string compound) {
    string compound_store; // Stores only chars of compound entry
    string stoich_store;
    std::map<string, int> compound_map; // an individual compound

    for (int i = 0; i < compound.size(); i++) {
      // take in chars
      if (std::isdigit(compound[i]) == 0) {
        compound_store += compound[i];
      } else {
        // ensure that we capture all digits of stoichiometry
        stoich_store += compound[i];

        // if last entry, set stoich to stoich_store
        if (i == compound.size() - 1) {
          compound_map.insert({compound_store, stoi(stoich_store)});
          compound_store = "";
          stoich_store = "";
        } else if (isdigit(compound[i + 1]) == 0) {
          compound_map.insert({compound_store, stoi(stoich_store)});
          compound_store = "";
          stoich_store = "";
        }
      }
    }

    // reduce stoichiometry of this compound
    compound_map = reducedStoich(compound_map);

    return compound_map;
  }

  // Json serialization of Flags
  aurostd::JSON::object Flags::serialize() const {
    return aurostd::JSON::object({AST_JSON_GETTER(JSON_Flags_MEMBERS)});
  }

  Flags Flags::deserialize(const aurostd::JSON::object& jo) {
    AST_JSON_SETTER(JSON_Flags_MEMBERS)
    return *this;
  }

  // Json serialization of Entry
  aurostd::JSON::object Entry::serialize() const {
    return aurostd::JSON::object({AST_JSON_GETTER(JSON_Entry_MEMBERS)});
  }

  Entry Entry::deserialize(const aurostd::JSON::object& jo) {
    AST_JSON_SETTER(JSON_Entry_MEMBERS)
    return *this;
  }

  /// @brief takes a vector of entries and a specie vector and converts them to stoichiometry-energy coordinates
  /// @return a list of hull points
  /// @authors
  /// mod{NHA,20260206,created}
  vector<aurostd::xvector<double>> entriesToPoints(const vector<Entry>& entries, const vector<string>& specie_vector) {
    vector<aurostd::xvector<double>> qhull_many_points(entries.size());

    // use the following standard: (specie_2, .... specie_n, energy)
    for (int i = 0; i < entries.size(); i++) {
      const aurostd::xvector<double> qhull_point = entries[i].compoundToPoint(specie_vector, false);
      qhull_many_points[i] = (qhull_point);
    }

    return qhull_many_points;
  }

  /// @brief takes a vector of entries and a specie vector and converts them to stoichiometry-energy coordinates
  /// @return a list of hull points
  /// @authors
  /// mod{NHA,20260206,created}
  vector<vector<double>> entriesToPoints(const std::map<std::map<string, int>, double>& entry_map, const vector<string>& specie_vector) {
    const int dimension = specie_vector.size(); // dimension of qhull system = (n-1)
    vector<vector<double>> qhull_many_points(entry_map.size());
    int i = 0; // iterator

    // use the following standard: (specie_1, specie_2, .... specie_n, energy)
    for (const auto& entry : entry_map) {
      vector<double> stoich_point(dimension, 0);
      for (const auto& specie : entry.first) {
        for (int i = 0; i < specie_vector.size(); i++) {
          // This if statement ensures consistent indexing of stoichiometries
          if (specie_vector[i] == specie.first) {
            stoich_point[i] = specie.second;
          }
        }
      }
      vector<double> qhull_point = stoichCoord(stoich_point);
      qhull_point.push_back(entry.second);
      qhull_many_points[i] = qhull_point;
      i++;
    }

    return qhull_many_points;
  }

  ///@brief returns true if entry contains 1 element only
  bool Entry::isUnary() const {
    return nspecies == 1;
  }

  ///@brief legacy function used for GFA: returns a hull point with the first element included
  aurostd::xvector<double> Entry::getStoichiometricCoords() const {
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (m_coord.rows == 1) { // trivial unary hull
      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " m_coord.rows==1: generating null xvector" << endl;
      }
      if (!isUnary()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Found non-unary point on unary hull", _INPUT_UNKNOWN_);
      }
      const aurostd::xvector<double> r_coords = aurostd::null_xv<double>(); // get null vector //(m_coord.lrows,m_coord.lrows);  //just return vector length 1 with 0 inside
      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " r_coords=" << r_coords << endl;
      }
      return r_coords;
    }
    aurostd::xvector<double> r_coords(m_coord.lrows, m_coord.urows - 1);
    for (int i = m_coord.lrows; i <= m_coord.urows - 1; i++) {
      r_coords[i] = m_coord[i];
    }
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " r_coords=" << r_coords << endl;
    }
    return r_coords;
  }

  /// @brief returns last index of m_coord. For thermodynamic calculations this is the energy
  double Entry::getLastCoord() const {
    return m_coord[m_coord.urows];
  }

  ///@returns the coordinate dimensions of this entry
  uint Entry::getCoordDimension() const {
    return m_coord.rows;
  }

  /// @brief returns a list of the phase decomposition ratios in the same order as the global specie_vector
  /// @authors
  /// mod{NHA,20260206,created}
  aurostd::xvector<double> Entry::findPhaseDecompRatios() const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Entry has not been initialized!");
    }
    auto map = m_nhull_phase_decomp_entries;
    vector<double> phase_list;
    //maps are already sorted by key implicitly, just pull out values here
    std::transform(map.begin(), map.end(), std::back_inserter(phase_list), [](const auto& pair) { return pair.second; });
    return aurostd::vector2xvector(phase_list);
  }

  ///@brief returns list of entries with the lowest energies in their respective compositions
  ///@param filter_out_below_dim filters out points below this dimension
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> Entry::getLowestEnergyEntries(const vector<Entry>& entries, int filter_out_below_dim) {
    vector<Entry> lowest_energy_entries;
    for (auto entry : entries) {
      if (entry.lowest_energy_alloy && entry.nspecies >= filter_out_below_dim) {
        lowest_energy_entries.push_back(entry);
      }
    }
    return lowest_energy_entries;
  }

  ///@brief returns list of entries with the lowest energies in their respective compositions and that are below a set dimension
  ///@param filter_out_above_dim filters out points above this dimension
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> Entry::getLowestEnergyNminus1Entries(const vector<Entry>& entries, int filter_out_above_dim) {
    vector<Entry> lowest_energy_entries;
    for (auto entry : entries) {
      if (entry.lowest_energy_alloy && entry.nspecies < filter_out_above_dim) {
        lowest_energy_entries.push_back(entry);
      }
    }
    return lowest_energy_entries;
  }

  ///@brief returns list of entries with the lowest energies in their respective compositions, avoiding compositions associated with filter_out_entry
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> Entry::getLowestEnergySCEntries(const vector<Entry>& entries, const Entry& filter_out_entry) {
    vector<Entry> lowest_energy_entries;
    for (auto entry : entries) {
      if (entry.lowest_energy_alloy && entry.compound_map != filter_out_entry.compound_map) {
        lowest_energy_entries.push_back(entry);
      }
    }
    return lowest_energy_entries;
  }

  /// @brief evaluates whether point lies inside the convex hull
  /// @param lower_hull true if this is a thermodynamic calculation with energies, and at least one point lies below 0-energy tie-line
  /// @return true if the point is below the zero energy tie line
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::isWithinHalfHull(bool lower_hull) const {
    return (lower_hull ? aurostd::lessEqualZero(getLastCoord()) : aurostd::greaterEqualZero(getLastCoord()));
  }

  ///@brief converts m_nhull_phase_decomp_entries to a string/double map for ease of displaying
  /// @authors
  /// mod{NHA,20260206,created}
  std::map<std::string, double> Entry::phaseDecompToMap() const {
    std::map<std::string, double> phase_decomp;
    vector<double> decomp_ratios;
    vector<aurostd::xvector<double>> vertices;

    for (const auto& [entry, ratio] : m_nhull_phase_decomp_entries) {
      phase_decomp[reducedCompoundString(entry.compound)] = ratio;
    }

    return phase_decomp;
  }

  ///@brief returns and entry from a list given an auid
  /// @authors
  /// mod{NHA,20260206,created}
  Entry Entry::getEntryFromAuid(const vector<Entry>& entries, const string auid) {
    for (const Entry& entry : entries) {
      if (entry.auid == auid) {
        return entry;
      }
    }
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No entry found in list!");
  }

  ///@brief legacy function for use in GFA: returns is_hull_point
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::isGState() const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Uninitialized point");
    }
    return is_hull_point; //|| m_is_g_state || m_is_equivalent_g_state;
  }

  ///@brief utility method used to flag duplicate entries in filterEntries()
  ///@param  unique_entries list of entry indices of unique entries
  ///@param other entry to test
  ///@param entries global filtered list of entries
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::entryUnique(const vector<uint>& unique_entries, Entry& other, const vector<Entry>& entries) {
    // points have already been created, determined to be unique
    // hack, go backwards, as the way entries are ordered, duplicates occur near each other
    for (size_t fl_size_i = unique_entries.size(), i = fl_size_i - 1; i < fl_size_i; i--) {
      const Entry& point = entries[unique_entries[i]];
      if (point.entryIdentical(other)) {
        //flag this point here:
        other.flags.duplicate_entry = true;
        return false;
      }
    }
    return true;
  }

  /// @brief returns true if other is the same as current entry
  /// @authors
  /// mod{NHA,20260206,created}
  bool Entry::entryIdentical(const Entry& other) const {
    return (compound == other.compound)
        && (prototype == other.prototype)
        && // believe it or not, compound + prototype is PROBABLY enough, but let's be sure
           (aurostd::identical(enthalpy_formation_atom, other.enthalpy_formation_atom, ENERGY_TOL));
  }

  ///@brief utility method used for sorting points in IQR filter
  ///@return list of integer indices corresponding to elements present in this entry
  /// @authors
  /// mod{NHA,20260206,created}
  aurostd::xvector<int> Entry::getElementsPresent() const {
    aurostd::xvector<double> stoich_coords = getStoichCoords(m_coord);
    for (int i = stoich_coords.lrows; i <= stoich_coords.urows; i++) {
      if (std::signbit(stoich_coords[i]) || stoich_coords[i] > 1.0) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Coord(" + aurostd::utype2string(i) + ") is outside of [0,1] range of a generalized stoichiometry coordinate");
      }
    }
    aurostd::xvector<int> elements_present(stoich_coords.lrows, stoich_coords.urows);
    for (int i = stoich_coords.lrows; i <= stoich_coords.urows; i++) {
      if (aurostd::nonZeroWithinTol(stoich_coords[i], ZERO_TOL)) {
        elements_present[i] = 1;
      }
    }
    return elements_present;
  }

  ///@brief converts a qhull coordinate (first specie removed, energy coord included), to a stoich coordinate (all species included, no energy coord)
  ///@param qhull_coord coordinate used in qhull calculation: no first (element) index, energy included as last index
  ///@return stoichiometric coordinates of this entry
  /// @authors
  /// mod{NHA,20260206,created}
  aurostd::xvector<double> Entry::getStoichCoords(aurostd::xvector<double> qhull_coord) {
    // The ordering of species here is the same as in the original compound vector
    // returns list of all elements present in fractional ratios specified in their compositions
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    double c_sum = 0.0; // concentration sum
    aurostd::xvector<double> stoich(qhull_coord.urows, qhull_coord.lrows);
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " qhull_coord=" << qhull_coord << endl;
    }
    for (int j = qhull_coord.lrows; j < qhull_coord.urows; j++) {
      if (std::signbit(qhull_coord[j])) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Negative stoich coordinate found");
      } // no negative numbers in stoich coordinates, only energy
      stoich[j + 1] = qhull_coord[j];
      c_sum += qhull_coord[j];
    }
    stoich[stoich.lrows] = (1.0 - c_sum); // hidden dimension set to first index
    if (std::signbit(stoich[stoich.lrows])) {
      // necessary check now because POCC entries have no composition, only stoich, so write-out errors will be prevalent
      if (aurostd::zeroWithinTol(stoich[stoich.lrows], ZERO_TOL)) {
        stoich[stoich.urows] = 0.0;
      } // only zero out for the last coord, as it's derived from the subtraction of the others
      else {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Negative stoich coordinate found");
      }
    } // no negative numbers
    return stoich;
  }
}// namespace nhull

template <> inline aurostd::JSON::object JsonSerializable<nhull::Entry>::dumpJsonMeta() const {
  return aurostd::JSON::object(aurostd::JSON::Dictionary{{}});
}

template <> inline nhull::Entry JsonSerializable<nhull::Entry>::loadFromJson(const aurostd::JSON::object& jo) {
  nhull::Entry tr{};
  return tr.deserialize(jo);
}
