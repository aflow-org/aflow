// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu
// Previous versions also written by Eric Perim, Eric Gossett, and Corey Oses

#include "aflow_nhull.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <extern/SYMBOLICCPLUSPLUS/symbolicc++.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib_entry_loader.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "extern/QHULL/libqhull_r/libqhull_r.h"
#include "extern/QHULL/libqhullcpp/Qhull.h"
#include "extern/QHULL/libqhullcpp/QhullPoint.h"
#include "extern/QHULL/libqhullcpp/QhullRidge.h"
#include "extern/QHULL/libqhullcpp/QhullVertexSet.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/HULL/aflow_nhull_entry.h"
#include "modules/HULL/aflow_nhull_facet.h"
#include "modules/HULL/aflow_nhull_latex.h"
#include "modules/HULL/aflow_nhull_stats.h"
#include "modules/HULL/aflow_nhull_util.h"

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::get;
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
using std::tuple;
using std::vector;

using aurostd::xvector;

namespace nhull {
  //************* MAIN PROGRAM LOOP FUNCTIONS IN CHRONOLOGICAL ORDER for entryloader stuff

  ///@brief flags entries that have either internal inconsistencies, inapropriate DFT run options, or have bad convergence
  ///@param entry the loaded entry to be evaluated and potentially flagged
  ///@param specie_vector the global user-submitted compound
  /// @authors
  /// mod{NHA,20260206,created}
  void entryValid(Entry& entry, const vector<string>& specie_vector) {
    //checks whether given compound contains species not present in specie_vector
    if (!entry.isCompoundCompatible(specie_vector)) {
      entry.flags.incompatible_entry = true;
    }

    //checks whether +U was used to calculate this entry
    if (entry.uses_LDAU()) {
      entry.flags.plus_U_used = true;
    }

    //Test is needed to ensure we have not loaded a pocc parent directory
    if (entry.vspecies.size() != entry.vcomposition.size()) {
      entry.flags.vspecie_vcomposition_size_error = true;
      if (entry.prototype.find("POCC") == string::npos) {
        entry.flags.pocc_prototype = true;
      }
    }

    if (!aurostd::WithinList(ALLOWED_DFT_TYPES, entry.lib_entry_direct.dft_type)) {
      entry.flags.uses_disallowed_dft_type = true;
    }

    if (aurostd::substring2bool(entry.aflow_lib_entry, "NUPDOWN")) {
      entry.flags.uses_NUPDOWN = true;
    }

    if (entry.aurl.find("ARUN.AEL_") != string::npos || entry.aurl.find("ARUN.AGL_") != string::npos || entry.aurl.find("ARUN.APL_") != string::npos || entry.aurl.find("ARUN.QHA_") != string::npos || false) {
      entry.flags.arun_entry = true;
    }

    //remove entries with very high energies (>1000 eVs)
    if (entry.is_energy_extraneous()) {
      entry.flags.extraneous_energy_entry = true;
    }

    //remove bad database entries:
    string bad_database_error_reason;
    if (entry.lib_entry_direct.ignoreBadDatabase(bad_database_error_reason)) {
      entry.flags.bad_database = true;
    }
  }

  /// @brief takes in vector of entries and determines which entry has lowest energy in respective compound
  /// @note updates value of global entries and returns smaller vector<Entry> object of sorted elements
  /// @note this represents the last level of filtering and determines all the points ingested into qhull
  /// @param filtered_entries ingested entries must have no extraneous entries to avoid inclusion into hull
  /// @param m_points
  /// @return list of entries of lowest energy of their compound
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> lowestEnergyAlloys(vector<Entry>& filtered_entries, vector<Entry>& m_points) {
    if (filtered_entries.empty()) {
      return vector<Entry>{};
    }
    vector<Entry> lowest_energy_points;

    for (Entry& entry : filtered_entries) {
      if (entry.nspecies != 1) {
        bool in_lowest_points = false;
        for (Entry& lowest_energy_point : lowest_energy_points) {
          if (lowest_energy_point.compound_map == entry.compound_map) {
            in_lowest_points = true;

            if (entry.enthalpy_formation_atom < lowest_energy_point.enthalpy_formation_atom) {
              lowest_energy_point = entry;
            }
          }
        }
        if (!in_lowest_points) {
          lowest_energy_points.push_back(entry);
        }
      } else if (entry.m_is_artificial) {
        lowest_energy_points.push_back(entry); //add artificial unaries here
      }
    }

    // now update parameters in both vectors of entries to reflect which alloys are lowest energy
    for (Entry& lowest_entry : lowest_energy_points) {
      lowest_entry.lowest_energy_alloy = true;
    }

    // now update entries with same info as in lowestEnergyAlloys
    for (const Entry& lowest_entry : lowest_energy_points) {
      for (Entry& entry : m_points) {
        if (lowest_entry.auid == entry.auid) {
          entry.lowest_energy_alloy = true;
        }
      }
    }

    return lowest_energy_points;
  }

  ///@brief Takes global list of all load entries and applies all filters that determine acceptable hull points. Global entry flags are only updated in this function.
  ///@param entries global list of entries
  ///@param specie_vector the global user-submitted compound
  ///@param do_IQR if true, applies interquartile range filtering to entries
  ///@param pocc_run true if user submitted --pocc= option at runtime. If true, adds pocc_prototype identification to flags.
  ///@param dimension the global dimension of this hull run: i.e a ternary "MnPdPt" would have dimension = 3.
  ///@param nminus1_calc
  ///@param sc_calc
  ///@return list of filtered entries
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> filterEntries(vector<Entry>& entries,
                              const vector<string>& specie_vector,
                              const bool do_IQR,
                              const bool pocc_run,
                              const int dimension,
                              const std::pair<bool, string> nminus1_calc = {false, ""},
                              const std::pair<bool, string>& sc_calc = {false, ""}) {
    vector<Entry> filtered_entries;
    bool lower_half_hull = false; //indicates whether points sit below 0 energy tie line.
    int lower_half_hull_count = 0;

    //check to make sure we aren't filtering both nminus1 and sc: this would produce nonsensical results
    if (nminus1_calc.first && sc_calc.first) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Tried to filter for nminus1 and stability criterion simultaneously!");
    }

    //Test whether given entry is suitable for hull calculations:
    vector<uint> unique_entries;
    for (int i = 0; i < entries.size(); i++) {
      Entry& entry = entries[i];
      if (!entry.m_is_artificial) { //only filter non-artifical entries
        entryValid(entry, specie_vector);
        if (entry.flags.is_flagged(pocc_run)) {
          continue;
        }
        if (!Entry::entryUnique(unique_entries, entry, entries)) {
          continue;
        }
      }
      unique_entries.push_back(i);
    }

    bool half_hull = ConvexHull::isHalfHull(entries, pocc_run);

    if (do_IQR) {
      IQRFilter(entries, pocc_run, half_hull, dimension);
    }

    if (sc_calc.first) {
      stabilityCriteriaFilter(entries, sc_calc.second);
    }

    if (nminus1_calc.first) {
      uint entry_dimension = Entry::getEntryFromAuid(entries, nminus1_calc.second).getCoordDimension();
      Nminus1Filter(entries, entry_dimension);
    }

    //only return entries with no flags thrown:
    //UPDATE ALL GLOBAL Entry FLAGS HERE:
    //After this point no global flags are updated
    for (auto& entry : entries) {
      if (!entry.flags.is_flagged(pocc_run)) {
        filtered_entries.push_back(entry);
        if (entry.enthalpy_formation_atom < 0) {
          lower_half_hull = true;
          lower_half_hull_count++;
        }
      } else {
        entry.flagged_entry = true;
      }
    }

    // Now only take lowest energy alloy of each compound
    vector<Entry> lowest_energy_alloys = lowestEnergyAlloys(filtered_entries, entries);
    return lowest_energy_alloys;
  }

  ///@brief returns all points except prototypes corresponding to excluded entry
  ///@param entries list of all entries to potentially flag
  ///@param excluded_auid the excluded entry
  /// @authors
  /// mod{NHA,20260206,created}
  void stabilityCriteriaFilter(vector<Entry>& entries, const string& excluded_auid) {
    std::string proto = Entry::getEntryFromAuid(entries, excluded_auid).lib_entry_direct.aflow_prototype_label_relax;
    for (Entry& entry : entries) {
      if (entry.lib_entry_direct.aflow_prototype_label_relax == proto) {
        entry.flags.excluded_stability_criterion = true;
      }
    }
  }

  ///@brief returns all points except compositions of dimension > N_compound
  ///@param entries list of all entries to potentially flag
  ///@param N_compound dimension to start excluding entries
  /// @authors
  /// mod{NHA,20260206,created}
  void Nminus1Filter(vector<Entry>& entries, uint N_compound) {
    for (Entry& entry : entries) {
      if (entry.vspecies.size() >= N_compound) {
        entry.flags.excluded_nminus1 = true;
      }
    }
  }

  /// @brief filter used before statistics: checks if entry was flagged
  /// @param entries vector of nhull_entries
  /// @return vector of sorted entries
  /// @authors
  /// @mod{NHA,20251005,created}
  vector<Entry> lightFilter(const vector<Entry>& entries) {
    vector<Entry> filtered_entries;

    for (const auto& entry : entries) {
      if (!entry.flagged_entry) {
        filtered_entries.push_back(entry);
      }
    }

    return filtered_entries;
  }

  ///@brief logging for user-submitted flags
  ///@param velements list of elements
  /// @authors
  /// mod{NHA,20260206,created}
  void flagCheck(aurostd::xoption& vpflow, const vector<string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss, bool silent) {
    const string directory = aurostd::getPWD();
    _aflags aflags;
    aflags.Directory = directory;
    if (vpflow.flag("NHULL::INIT")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Initiallizing NHULL", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if (vpflow.flag("NHULL::PRINT")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Printing plots and tables", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if (vpflow.flag("NHULL::STATS")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Compiling compound statistics including EFA and DEED", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if (vpflow.flag("NHULL::N-1_ENTHALPY_GAIN")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Creating N-1 Pseudohull", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
    if (vpflow.flag("NHULL::POCC")) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Running NHULL on POCC directory", aflags, FileMESSAGE, oss, _LOGGER_OPTION_, silent);
    }
  }

  //*************** QHULL SECTION *******************

  //************* QHULL SPECIFIC UTILITY METHODS

  ///@brief To run an n-dimension hull, it is necessary to run hulls of lower dimension first. This function returns the full n-dimension representation of a point after running a lower-dimension hull.
  ///@param alloy_name name of compound of the current lower-dimension hull run
  ///@param specie_vector global compound at dimension n
  ///@param facet_points points determined to be on the hull from the current hull run
  ///@return vector of n-dimension hull points
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> unsquashVec(const vector<string>& alloy_name, const vector<string>& specie_vector, const vector<vector<double>>& facet_points) {
    const int original_dim = specie_vector.size();
    vector<vector<double>> squashed_vectors;
    vector<xvector<double>> new_qhull_vectors;
    vector<vector<double>> unsquashed_vectors;
    const vector<double> stoich_coords;
    double energy_coord;
    const vector<double> coord_store;
    const vector<double> flattened_vector;
    const vector<double> energy_coords;
    vector<double> unsquashed_store;

    vector<bool> deleted_indices(original_dim);
    // first determine the indices we are restoring
    for (int i = 0; i < original_dim; i++) {
      if (count(alloy_name.begin(), alloy_name.end(), specie_vector[i]) == 0) {
        deleted_indices[i] = true;
      } else {
        deleted_indices[i] = false;
      }
    }

    // then convert from qhull coords to full stoich coords
    for (auto coord : facet_points) {
      energy_coord = coord.back();
      coord.pop_back(); // remove energy coord here
      vector<double> true_stoich_coord = fullCoords(coord, false);

      true_stoich_coord.push_back(energy_coord);
      squashed_vectors.push_back(true_stoich_coord);
    }

    // using our deleted indices vector, strategically reinsert where zeroes should be in full dimension vector
    int squashed_index = 0; // index in the original squashed vector
    for (const auto& element : squashed_vectors) {
      unsquashed_store.clear();
      squashed_index = 0;
      for (int i = 0; i < deleted_indices.size(); i++) {
        if (deleted_indices[i] == false) {
          unsquashed_store.push_back(element[squashed_index]);
          squashed_index++;
        } else {
          unsquashed_store.push_back(0);
        }
      }
      // make sure to reintroduce energy coordinate:
      unsquashed_store.push_back(element.back());
      unsquashed_vectors.push_back(unsquashed_store);
    }

    // convert back to qhull coords by removing first index
    for (auto element : unsquashed_vectors) {
      element.erase(element.begin());
      const xvector<double> xvector_coord = aurostd::vector2xvector(element, 0);
      new_qhull_vectors.push_back(xvector_coord);
    }

    return new_qhull_vectors;
  }

  ///@brief returns all hull points from a qhull run and returns them as qhull coordinates which are used as input to higher-dimension hulls in iterateRun
  ///@param qhull an instance of Qhull which contains a list of facets
  ///@return a vector of points
  /// @authors
  /// mod{NHA,20260206,created}
  vector<vector<double>> getFacetPoints(const orgQhull::Qhull& qhull) {
    const orgQhull::QhullFacetList facet_list = qhull.facetList();
    const orgQhull::QhullFacet first_facet = qhull.beginFacet();
    const int dimension = first_facet.dimension();
    orgQhull::QhullFacet curr_facet = first_facet;
    orgQhull::QhullVertexSet curr_vertexset = first_facet.vertices();
    vector<double> all_facet_points_store;
    vector<vector<double>> all_facet_points;

    for (int i = 0; i < facet_list.count(); i++) {
      curr_vertexset = curr_facet.vertices();
      const vector<orgQhull::QhullVertex> curr_vertexvector = curr_vertexset.toStdVector();
      // filter out facets that have normal that point up, these are not lower half plane facets!
      if (getNormEnergy(getQhullFacetNormal(curr_facet)) > 0) {
        for (const auto& curr_vertex : curr_vertexvector) {
          all_facet_points_store.clear();
          const orgQhull::QhullPoint point = curr_vertex.point();
          coordT* point_coordinates = curr_vertex.point().coordinates();
          for (int i = 0; i < dimension; i++) {
            // current point in list for facet
            all_facet_points_store.push_back(*point_coordinates);
            point_coordinates++;
          }
          all_facet_points.push_back(all_facet_points_store);
        }
      }
      curr_facet = curr_facet.next();
    }

    return all_facet_points;
  }

  ///@brief To run an n-dimension hull, it is necessary to run hulls of lower dimension first. This function converts n-dimension points to m<n dimension points.
  ///@param converted_dim the dimension (m) being translated to from the highest dimension (n); (m<n)
  ///@param alloy_name name of compound of the current lower-dimension hull run
  ///@param specie_vector global compound at dimension n
  ///@param qhull_coords List of n-dimension points to be ingested into qhull
  ///@return tuple where first element is the converted dimension, the second is a flattened vector of converted points
  /// @authors
  /// mod{NHA,20260206,created}
  tuple<int, vector<double>> squashQhull(const int converted_dim, const vector<string>& alloy_name, const vector<string>& specie_vector, const vector<vector<double>>& qhull_coords) {
    const int original_dim = specie_vector.size();
    vector<vector<double>> new_qhull_vectors;
    const vector<vector<double>> unsquashed_vector;
    const vector<double> stoich_coords;
    vector<double> flattened_vector;
    const vector<vector<double>> facet_points;

    // get significant indices here:
    vector<bool> deleted_indices(original_dim);
    // first determine the indices we are removing
    for (int i = 0; i < original_dim; i++) {
      if (count(alloy_name.begin(), alloy_name.end(), specie_vector[i]) == 0) {
        deleted_indices[i] = true;
      } else {
        deleted_indices[i] = false;
      }
    }

    // only convert to squashed if there is a difference in coordinate size
    if (converted_dim != original_dim) {
      vector<double> coord_store;
      vector<vector<double>> squashed_vectors;
      for (auto point : qhull_coords) {
        const double energy_coord = point.back();
        point.pop_back(); // remove energy point here
        vector<double> true_stoich_coord = fullCoords(point, false);
        coord_store.clear();

        // only add significant indices to vector
        for (size_t i = 0; i < true_stoich_coord.size(); i++) {
          if (deleted_indices[i] == false) {
            coord_store.push_back(true_stoich_coord[i]);
          }
        }
        // add back energy coordinate
        coord_store.push_back(energy_coord);
        squashed_vectors.push_back(coord_store);
      }

      // by convention we remove the first index specie and define this as our qhull vector
      for (auto element : squashed_vectors) {
        element.erase(element.begin());
        new_qhull_vectors.push_back(element);
      }
    } else {
      new_qhull_vectors = qhull_coords;
    }

    flattened_vector = flatten_vector(new_qhull_vectors);

    // return qhull input for later use in run
    return {converted_dim, flattened_vector};
  }

  ///@brief function ensures that no duplicate coordinates are added to unsquashed vector (list of points carried over between hull calcs)
  ///@param unsquashed_facet_points list of all found facet points. These are updated here, and are selectively added to lower_dim_point_list before other hull runs
  ///@param unsquashed_points list of all unsquashed_facet_points from previous, lower dimension hulls
  ///@return list of qhull points at highest dimension in hull run
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> unsquashedInsert(vector<xvector<double>>& unsquashed_facet_points, const vector<xvector<double>>& unsquashed_points) {
    bool add = true;
    for (const auto& element_point : unsquashed_points) {
      add = true;
      for (const auto& element_facet : unsquashed_facet_points) {
        if (element_facet == element_point) {
          add = false;
          break;
        }
      }
      if (add) {
        unsquashed_facet_points.push_back(element_point);
      }
    }
    return unsquashed_facet_points;
  }

  ///@brief Utility function used in lower_dim_insert. Determines whether given points are permissible for an alloy
  ///@param deleted_indices list that indicates which indices have been removed in accordance to this m<n hull calc
  ///@return list of qhull points that can be used for the current composition hull run
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> lowDimPermissible(const int original_dim, const vector<xvector<double>>& points, const vector<bool>& deleted_indices) {
    vector<xvector<double>> unsquashed_permissable_points;
    vector<bool> found_deleted_indices(original_dim);

    for (size_t i = 0; i < points.size(); i++) {
      bool out_of_bounds_flag = false;
      xvector<double> true_stoich_coord = Entry::getStoichCoords(points[i]);
      for (size_t i = true_stoich_coord.lrows; i <= true_stoich_coord.urows; i++) {
        if (true_stoich_coord[i] == 0) {
          found_deleted_indices[i] = true;
        } else {
          found_deleted_indices[i] = false;
        }
      }
      for (size_t i = 0; i < found_deleted_indices.size(); i++) {
        if (found_deleted_indices[i] == false && deleted_indices[i] == true) {
          out_of_bounds_flag = true;
          break;
        }
      }
      // only add point if it does not violate "not found in alloy" condition
      if (out_of_bounds_flag == false) {
        unsquashed_permissable_points.push_back(points[i]);
      }
    }
    return unsquashed_permissable_points;
  }

  ///@brief tests whether given qhull points are permissible for this alloy, and inserts these points into lower_dim point array which is used to carry lower dimension points over in calculations
  ///@param specie_vector global compound at dimension n
  ///@param alloy_name_store allowed elements for the current iteration of the hull calculation
  ///@param unsquashed_facet_points list of all found facet points. These are updated here, and are selectively added to lower_dim_point_list before other hull runs
  ///@return list of permissible hull points for the current hull iteration
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> lower_dim_insert(const vector<string>& specie_vector, const vector<string>& alloy_name_store, const vector<xvector<double>>& unsquashed_facet_points) {
    const size_t original_dim = specie_vector.size();
    const vector<xvector<double>> lower_dim_permissable_points;
    vector<bool> deleted_indices(original_dim);
    const vector<bool> found_deleted_indices(original_dim); // recorded indices of given qhull vector

    // first determine the significant indices
    for (size_t i = 0; i < original_dim; i++) {
      if (count(alloy_name_store.begin(), alloy_name_store.end(), specie_vector[i]) == 0) {
        deleted_indices[i] = true;
      } else {
        deleted_indices[i] = false;
      }
    }

    // now test qhull points:
    vector<xvector<double>> unsquashed_permissable_points = lowDimPermissible(original_dim, unsquashed_facet_points, deleted_indices);

    return unsquashed_permissable_points;
  }

  ///@brief returns all lower-dimension points from iterateRunNhullFacet function; ensures unaries are present aswell
  ///@param unsquashed_facet_points list of all found facet points. These are updated here, and are selectively added to lower_dim_point_list before other hull runs
  /// @return list of qhull points (no first compound, energy included)
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> getLowerDimPoints(const vector<string>& specie_vector, const vector<xvector<double>>& unsquashed_facet_points) {
    size_t dimension = specie_vector.size();
    vector<xvector<double>> all_lower_dim_points = unsquashed_facet_points;

    //now make sure all unaries are placed:
    for (size_t i = 0; i < dimension; i++) {
      bool missing_pure = true;
      xvector<double> coordinates(dimension - 1, 0);
      if (i != dimension - 1) {
        coordinates[i + coordinates.lrows] = 1;
      }

      for (xvector<double> point : all_lower_dim_points) {
        if (aurostd::identical(point, coordinates)) {
          missing_pure = false;
          break;
        }
      }
      if (missing_pure) {
        all_lower_dim_points.push_back(coordinates);
      }
    }
    return all_lower_dim_points;
  }

  //************************* MAIN RUN FUNCTIONS

  ///@brief function provides and interface for calling Qhull within the context of a thermodynamic hull calculation
  ///@param qhull instance of qhull to be run
  ///@param alloy_store map of allowed compounds for the current iteration of the hull calculation
  ///@param alloy_name_store allowed elements for the current iteration of the hull calculation
  ///@param specie_vector global compound at dimension n
  ///@param lower_dim_point_list list of points from previous hull calculations done at a lower dimension than the current
  ///@param unsquashed_facet_points list of all found facet points. These are updated here, and are selectively added to lower_dim_point_list before other hull runs
  /// @authors
  /// mod{NHA,20260206,created}
  void qhullCalcSimple(orgQhull::Qhull& qhull,
                       const std::map<std::map<string, int>, double>& alloy_store,
                       const vector<string>& alloy_name_store,
                       const vector<string>& specie_vector,
                       const vector<xvector<double>>& lower_dim_point_list,
                       vector<xvector<double>>& unsquashed_facet_points) {
    const char* qhull_command2 = "Qt"; // triangulated output: necessary for correct phase decomposition

    // now calculate qhull for these points and save to array
    // WE ARE IN DIMENSION i HERE
    vector<vector<double>> point_list = entriesToPoints(alloy_store, specie_vector);
    // point list coords is unsquashed always
    auto conv_lower_dim_point_list = compoundXvectorConversion(lower_dim_point_list);
    point_list.insert(point_list.end(), conv_lower_dim_point_list.begin(), conv_lower_dim_point_list.end());

    // Run qhull here
    tuple<int, vector<double>> squash_qhull = squashQhull(alloy_name_store.size(), alloy_name_store, specie_vector, point_list);
    const int squashed_dim = get<0>(squash_qhull);
    vector<double> qhull_point_list = get<1>(squash_qhull);
    const int point_count = qhull_point_list.size() / squashed_dim;
    const double* qhull_input_points = &qhull_point_list[0];
    qhull.runQhull("run", squashed_dim, point_count, qhull_input_points, qhull_command2);

    const vector<vector<double>> facet_points = getFacetPoints(qhull);
    const vector<xvector<double>> unsquashed_points = unsquashVec(alloy_name_store, specie_vector, facet_points);
    unsquashed_facet_points = unsquashedInsert(unsquashed_facet_points, unsquashed_points);
  }

  /// @brief main hull calculation loop: goes through binary, ternary, quaternary, etc and calculates hulls progressively before returning a list of all facets at the highest dimension
  /// @param entries vector of nhull_entries
  /// @param specie_vector global compound at dimension n
  /// @param dimension the global dimension of this hull run: i.e a ternary "MnPdPt" would have dimension = 3.
  /// @param filtered_entries list of all entries deemed acceptable to use for hull generation
  /// @return list of facets
  /// @authors
  /// mod{NHA,20260206,created}
  vector<NhullFacet> iterateRunFacetsOnly(vector<Entry>& entries, const vector<string> specie_vector, const int dimension, const vector<Entry>& filtered_entries) {
    map<map<string, int>, double> alloy_dim_store;
    map<map<string, int>, double> alloy_store;
    map<map<string, int>, double> pure_stoich_store; // storage var for pure elements
    vector<string> curr_alloy_name;
    vector<string> alloy_name_store;
    vector<xvector<double>> lower_dim_point_list;
    vector<xvector<double>> unsquashed_facet_points;
    vector<NhullFacet> nhull_facets; // facet store for internal functionality: Entry object and pdf plot

    if (dimension < 2) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Tried to initialize a hull with less than 2 dimensions!");
    }

    // LOOP START DIMENSION: 2
    // loop goes up to dimension-1 (we calculate max dimension after all loops are done to determine facets as a last step)
    // this must be done to avoid degerate planes of points being placed into higher dimension hulls. See: DOI: 10.1021/acs.jcim.8b00393
    for (int i = 2; i < dimension; i++) {
      bool missing_alloys = false;
      // first determine relevant points:
      for (const Entry& element : filtered_entries) {
        // only keep ith dimensional and 1 dimension (pure) alloys
        if (element.nspecies == i) {
          alloy_dim_store[element.compound_map] = element.enthalpy_formation_atom;
        }
        // On first go-around make sure to add pure elements that would otherwise be ignored.
        if (element.nspecies == 1 && i == 2) {
          pure_stoich_store[element.compound_map] = element.enthalpy_formation_atom;
        }
      }

      // now sort ith dimension alloys by their species, (each specie combo is unique; Ex: Mn10Pd2 sorted with Mn3Pd1)
      //add additional check: if no compounds of specie n exist in filtered points, check if they exist in entries
      while (alloy_dim_store.empty() == false) {
        alloy_name_store.clear();
        alloy_store.clear();
        for (const auto& alloy : alloy_dim_store) {
          curr_alloy_name.clear();
          for (const auto& specie : alloy.first) {
            curr_alloy_name.push_back(specie.first);
          }
          if (alloy_name_store.empty()) {
            alloy_name_store = curr_alloy_name; // ENSURE THAT NAME STORE IS SORTED PROPERLY
            alloy_store.insert(alloy);
          } else {
            if (curr_alloy_name == alloy_name_store) {
              alloy_store.insert(alloy);
            }
          }
        }

        // now add any relevent pure stoich coordinates:
        // These correspond to single species that are contained in alloy_name_store
        // to prevent duplicate pure coordinates, only add these at the lowest dimension level (2)
        if (i == 2) {
          for (const string& element : alloy_name_store) {
            for (auto compound : pure_stoich_store) {
              if (compound.first.count(element) > 0) {
                alloy_store.insert(compound);
              }
            }
          }
        }

        // if we are above 2 dimensions here, incorporate hull points from previous calculations
        // Ensure that we add only relevent points: Ex: no binaries "Cr, P" in ternary "Cr, Mn, Ni"
        if (i > 2) {
          // lower dim_point_list is unsquashed always
          lower_dim_point_list = lower_dim_insert(specie_vector, alloy_name_store, unsquashed_facet_points);
        }

        //*************QHULL CALCULATION DONE HERE
        orgQhull::Qhull qhull;
        std::sort(alloy_name_store.begin(), alloy_name_store.end());
        qhullCalcSimple(qhull, alloy_store, alloy_name_store, specie_vector, lower_dim_point_list, unsquashed_facet_points);

        // remove found elements from alloy_dim_store before we loop again
        for (const auto& alloy : alloy_store) {
          alloy_dim_store.erase(alloy.first);
        }
      }
    }

    vector<xvector<double>> lower_dim_facet_points = getLowerDimPoints(specie_vector, unsquashed_facet_points);
    nhull_facets = NhullFacet::getFacetData(specie_vector, lower_dim_facet_points, entries, dimension);
    return nhull_facets;
  }

  ///@brief wrapper function for the main function: translates all user input into parameters and runs nhull
  ///@param vpflow All user input flags
  ///@return ConvexHull instance
  ///@author
  ///@mod{NHA,20251213,created}
  ConvexHull convexHull(const aurostd::xoption& _vpflow) {
    string alloy;
    vector<string> velements;
    ConvexHull run;
    aurostd::xoption vpflow(_vpflow);
    _aflags aflags;

    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    ostream& oss = cout;
    ofstream FileMESSAGE;
    stringstream message;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (vpflow.flag("NHULL::POCC")) {
      string pocc_root = vpflow.getattachedscheme("NHULL::POCC");
      if (!aurostd::IsDirectory(pocc_root)) {
        message << "POCC directory entered does not exist!";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
      }
    }

    string input = vpflow.getattachedscheme("PFLOW::ALLOY");
    velements = aurostd::getElements(input, pp_string, FileMESSAGE, true, true, false, oss);
    alloy = aurostd::joinWDelimiter(velements, "");

    if (velements.empty()) {
      message << "Invalid input (" << input << "), please capitalize element symbols";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }
    if (velements.size() < 2) {
      message << "Trivial input (" << input << "), enter binaries or higher";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_ERROR_);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }

    if (!vpflow.flag("NHULL::POCC")) {
      if (vpflow.flag(convex_flags.m_stability_criterion_flag)) {
        string auid = vpflow.getattachedscheme(convex_flags.m_stability_criterion_flag);
        if (auid.empty()) {
          message << "Empty auid found";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
        }
      }
    }

    // spit out banner for only the first request
    message << aflow::Banner("BANNER_NORMAL");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_RAW_, true);
    // i //no screen, first to be logged
    message << "Starting " << alloy << " " << pflow::arity_string(velements.size(), false, false) << " convex hull";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    //getPath(vpflow, FileMESSAGE, oss, false); // CO20180220 - directory stuff for logging
    nhull::flagCheck(vpflow, velements, FileMESSAGE, oss, false); // spit out all flag options

    run = ConvexHull(vpflow, aflags, velements);

    if (vpflow.flag("NHULL::POCC")) {
      run.calculatePOCC();
    } else {
      run.createThermoHull();
    }
    return run;
  }

  //****************************************** Class ConvexHull implementations:

  /// @brief generates list of facets, if two neighboring facets are coplanar join them
  /// @param facet_collection output vector containing lists of vertex indexes
  /// @param angle_threshold max angle between two facts in radian to be still considered coplanar
  /// @note the facets contain only vertices
  /// @authors
  /// mod{HE,20210510,created}
  void ConvexHull::getJoinedFacets(vector<vector<uint>>& facet_collection, const double angle_threshold) {
    // HE20210510
    const bool LDEBUG = (false || XHOST.DEBUG);

    if (m_dimension != 3) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "facet joining is just available in 3D", _VALUE_RANGE_);
    }
    vector<vector<uint>> raw_facets;
    vector<xvector<double>> normals;
    std::map<uint, vector<uint>> point_neighbors;
    std::list<std::pair<uint, uint>> raw_join_list;
    vector<vector<uint>> join_list;
    vector<uint> remove_facet;
    vector<vector<uint>> raw_facet_collection;
    std::map<uint, uint> corner_check;

    // Collect information on each facet
    for (std::vector<NhullFacet>::const_iterator facet = m_facets.begin(); facet != m_facets.end(); ++facet) {
      vector<uint> vertices;

      // Represent facet based on the vertices indexes
      for (std::vector<nhull::Entry>::const_iterator vert = facet->entries.begin(); vert != facet->entries.end(); ++vert) {
        vertices.push_back(vert->nh_index);
      }
      raw_facets.push_back(vertices);
      normals.push_back(facet->m_normal);
    }

    // Build lookup for neighboring points
    // (Base point is included to make a check easier)
    for (std::vector<vector<uint>>::const_iterator facet = raw_facets.begin(); facet != raw_facets.end(); ++facet) {
      for (std::vector<uint>::const_iterator ind = facet->begin(); ind != facet->end(); ++ind) {
        std::copy(facet->begin(), facet->end(), std::inserter(point_neighbors[*ind], point_neighbors[*ind].end()));
      }
    }

    // Ensure that point_neighbors is unique
    for (std::map<uint, vector<uint>>::iterator point_neighbors_entry = point_neighbors.begin(); point_neighbors_entry != point_neighbors.end(); ++point_neighbors_entry) {
      std::sort(point_neighbors_entry->second.begin(), point_neighbors_entry->second.end());
      point_neighbors_entry->second.erase(std::unique(point_neighbors_entry->second.begin(), point_neighbors_entry->second.end()), point_neighbors_entry->second.end());
    }

    const uint raw_facets_size = raw_facets.size();
    // Check for each facet, if their neighbors have an equivalent normal vector
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " coplanar | angle | facets | n1 | n2" << endl;
    }
    for (size_t i1 = 0; i1 < raw_facets_size; i1++) {
      const vector<uint> base_facet = raw_facets[i1];
      for (size_t i2 = i1 + 1; i2 < raw_facets_size; i2++) {
        uint check = 0;
        const vector<uint> compare_facet = raw_facets[i2];
        for (uint k = 0; k < 3; k++) {
          for (uint j = 0; j < 3; j++) {
            if (base_facet[k] == compare_facet[j]) {
              check++;
            }
          }
        }
        // neighbor share two vertices
        if (check != 2) {
          continue;
        }
        const double check_angle = aurostd::angle(normals[i1], normals[i2]);
        if (check_angle < angle_threshold) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " YES | ";
          }
          raw_join_list.emplace_back(i1, i2);
          remove_facet.push_back(i1);
          remove_facet.push_back(i2);
        } else {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " NO  | ";
          }
        }
        if (LDEBUG) {
          cerr << check_angle << " | ";
          cerr << i1 << ", " << i2 << " | ";
          for (uint k = 1; k < 4; k++) {
            cerr << normals[i1][k] << ", ";
          }
          cerr << "| ";
          for (uint k = 1; k < 4; k++) {
            cerr << normals[i2][k] << ", ";
          }
          cerr << endl;
        }
      }
    }

    // Combine the joined pairs into complete facets
    while (!raw_join_list.empty()) {
      const std::pair<uint, uint> start = *raw_join_list.begin();
      vector<uint> new_facet;
      new_facet.push_back(start.first);
      new_facet.push_back(start.second);
      raw_join_list.erase(std::find(raw_join_list.begin(), raw_join_list.end(), start));
      std::vector<std::pair<uint, uint>> to_delete;
      uint found = 1;
      while (found) {
        found = 0;
        to_delete.clear();
        for (std::list<std::pair<uint, uint>>::const_iterator next_ptr = raw_join_list.begin(); next_ptr != raw_join_list.end(); ++next_ptr) {
          const std::pair<uint, uint> next = *next_ptr;
          if (std::find(new_facet.begin(), new_facet.end(), next.first) != new_facet.end()) {
            found++;
            new_facet.push_back(next.second);
            to_delete.push_back(next);
          } else if (std::find(new_facet.begin(), new_facet.end(), next.second) != new_facet.end()) {
            found++;
            new_facet.push_back(next.first);
            to_delete.push_back(next);
          }
        }
        for (std::vector<std::pair<uint, uint>>::const_iterator next = to_delete.begin(); next != to_delete.end(); ++next) {
          raw_join_list.erase(std::find(raw_join_list.begin(), raw_join_list.end(), *next));
        }
      }
      join_list.push_back(new_facet);
    }

    // Build the raw facet collection
    // HE20220319 sorting needed for consistent 3D rendering (the order defines the facet normal)
    for (size_t i = 0; i < raw_facets_size; i++) {
      sortFacetVertices(raw_facets[i], i);
      if (std::find(remove_facet.begin(), remove_facet.end(), i) == remove_facet.end()) {
        raw_facet_collection.push_back(raw_facets[i]);
      }
    }
    // Add joined facets
    for (std::vector<vector<uint>>::const_iterator to_join = join_list.begin(); to_join != join_list.end(); ++to_join) {
      vector<uint> vertices;
      for (vector<uint>::const_iterator f_id = to_join->begin(); f_id != to_join->end(); ++f_id) {
        for (std::vector<uint>::const_iterator p_id = raw_facets[*f_id].begin(); p_id != raw_facets[*f_id].end(); ++p_id) {
          if (std::find(vertices.begin(), vertices.end(), *p_id) == vertices.end()) {
            vertices.push_back(*p_id); // ensure vertices vector is unique
          }
        }
      }
      // sort the vertices for std::set_difference(); point_neighbors are already sorted
      std::sort(vertices.begin(), vertices.end());

      // If a vertex has no outside neighbor it is removed
      vector<uint> vertices_to_remove;
      for (vector<uint>::const_iterator p_id = vertices.begin(); p_id != vertices.end(); ++p_id) {
        std::vector<uint> diff_result;
        std::set_difference(point_neighbors[*p_id].begin(), point_neighbors[*p_id].end(), vertices.begin(), vertices.end(), std::inserter(diff_result, diff_result.end()));
        if (diff_result.empty()) {
          vertices_to_remove.push_back(*p_id);
        }
      }
      for (std::vector<uint>::const_iterator p_id = vertices_to_remove.begin(); p_id != vertices_to_remove.end(); ++p_id) {
        vertices.erase(std::find(vertices.begin(), vertices.end(), *p_id));
      }
      // HE20220319 never add facet with less than three vertices
      if (vertices.size() >= 3) {
        sortFacetVertices(vertices, *to_join->begin());
        raw_facet_collection.push_back(vertices);
      }
    }

    // remove points that are on an edge and not a corner
    for (size_t i_facet = 0; i_facet < raw_facet_collection.size(); i_facet++) {
      for (size_t i_vert = 0; i_vert < raw_facet_collection[i_facet].size(); i_vert++) {
        corner_check[raw_facet_collection[i_facet][i_vert]]++;
      }
    }

    for (size_t i_facet = 0; i_facet < raw_facet_collection.size(); i_facet++) {
      vector<uint> vec_vertices;
      for (size_t i_vert = 0; i_vert < raw_facet_collection[i_facet].size(); i_vert++) {
        if (corner_check[raw_facet_collection[i_facet][i_vert]] > 2) {
          vec_vertices.push_back(raw_facet_collection[i_facet][i_vert]);
        }
      }
      facet_collection.push_back(vec_vertices);
    }
    // HE20220413 START
    //  Check for ghost facets with less than 3 vertices
    //  (they can be created when joining facets that are not perfectly coplanar)
    //  start with the largest index when removing, as vector will shrink and indexes would change
    for (size_t i_facet = facet_collection.size() - 1; i_facet > 0; i_facet--) {
      if (facet_collection[i_facet].size() <= 2) {
        facet_collection.erase(facet_collection.begin() + i_facet);
      }
    }
    // HE20220413 END
  }

  ///@brief returns the hull coordinates of an Entry at the highest dimension level (first specie coordinate missing)
  /// @authors
  /// mod{NHA,20260206,created}
  vector<aurostd::xvector<double>> ConvexHull::getAllCoordinates() {
    vector<aurostd::xvector<double>> all_coords;
    for (const auto& entry : m_points) {
      all_coords.push_back(entry.m_coord);
    }
    return all_coords;
  }

  ///@brief point-intializer for non-thermodynamic hull calculations
  ///@param vcoords list of coordinates to initialize entries with
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::initializePointsSimple(const vector<aurostd::xvector<double>>& vcoords) {
    m_points.clear();
    m_dimension = vcoords[0].rows;
    size_t i = 0;
    for (const auto& coord : vcoords) {
      if (coord.rows != m_dimension) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Input points do not share dimension!");
      } else {
        Entry new_entry(coord, std::cout, false, false, false, i);
        m_points.push_back(new_entry);
        i++;
      }
    }
  }

  /// @brief sort the vertex indexes based on their angle around a central normal vector
  /// @param facet list of vertex indexes forming the facet
  /// @param facet_id id of one of the original facets (to use its already calculated normal vector)
  /// @authors
  /// mod{HE,20210510,created}
  void ConvexHull::sortFacetVertices(vector<uint>& facet, const uint& facet_id) {
    //Sorting the vertex indexes avoid crossing lines when drawing or calculating the area or volume later.
    xvector<double> center(3, 1);
    const uint num_points = facet.size();
    xvector<uint> index_list(num_points, 1);
    xvector<double> angle_list(num_points, 1);
    for (std::vector<uint>::const_iterator p_id = facet.begin(); p_id != facet.end(); ++p_id) {
      center += m_points[*p_id].m_coord;
    }

    center /= num_points;
    xvector<double> start_vector = m_points[facet[0]].m_coord - center;
    xvector<double> normal = m_facets[facet_id].m_normal;

    // first index is used for the start_vector, therefore the angle is set to 0.0
    angle_list[1] = 0;
    index_list[1] = facet[0];
    for (uint i = 1; i < num_points; i++) {
      xvector<double> next_vector = m_points[facet[i]].m_coord - center;
      const double dot = aurostd::scalar_product(start_vector, next_vector);
      const double det = start_vector[1] * next_vector[2] * normal[3]
                       + next_vector[1] * normal[2] * start_vector[3]
                       + normal[1] * start_vector[2] * next_vector[3]
                       - start_vector[3] * next_vector[2] * normal[1]
                       - next_vector[3] * normal[2] * start_vector[1]
                       - normal[3] * start_vector[2] * next_vector[1];
      angle_list[i + 1] = std::atan2(det, dot);
      index_list[i + 1] = facet[i];
    }
    aurostd::quicksort2(num_points, angle_list, index_list);
    for (uint i = 0; i < num_points; i++) {
      facet[i] = index_list[i + 1];
    }
  }

  ///@brief sets the directory to save output in for a ConvexHull run
  ///@param vpflow list of user inputs
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::setDirectory(aurostd::xoption vpflow) {
    bool get_path = vpflow.flag(convex_flags.m_nhull_path, true);
    string directory = vpflow.getattachedscheme(convex_flags.m_nhull_path);

    if (directory.empty() || !aurostd::IsDirectory(directory)) {
      if (get_path) {
        std::cerr << "WARNING: Destination folder was not specified or does not exist! Outputting to working directory." << std::endl;
      }
      m_aflags.Directory = aurostd::getPWD(); //+ (false ? string("/") : string(""));
    } else {
      m_aflags.Directory = directory;
    }
  }

  ///@brief simple path tool: returns working directory
  string getPath(bool add_backslash) {
    return aurostd::getPWD() + (add_backslash ? string("/") : string(""));
  }

  ///@brief returns desired directory from user input
  string getPath(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& oss, bool silent) {
    // main function
    stringstream message;
    string path = vpflow.getattachedscheme("NHULL::PATH");
    if (!vpflow.flag("NHULL::PATH") && std::filesystem::is_directory(path)) {
      string pwd = getPath();
      if (!silent) {
        message << "Directing output to current directory: " << pwd;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_OPTION_);
      }
      return pwd;
    }
    return vpflow.getattachedscheme("NHULL::PATH");
  }

  ///@brief setter saves all user input to member var
  ///@param vpflow user options
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::setHullFlags(const aurostd::xoption& vpflow) {
    m_nflags = vpflow;
  }

  ///@brief simplest type of qhull run: just points, no thermodynamic analysis
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::calculateSimpleHull() {
    const char* qhull_command2 = "Qt";
    orgQhull::Qhull qhull;
    // Run qhull here
    const int dimension = m_dimension;
    vector<double> qhull_point_list = flatten_vector(this->getAllCoordinates());
    if (dimension == 0) {
      stringstream message;
      message << "Dimension is zero. This will cause divide by zero error!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }
    const int point_count = qhull_point_list.size() / dimension;
    const double* qhull_input_points = &qhull_point_list[0];

    //check to make sure we can calculate hull (must have n+1 points for n-dimension hull)
    if (point_count > dimension) {
      qhull.runQhull("run", dimension, point_count, qhull_input_points, qhull_command2);
    } else {
      //if we can't, see if we can calculate a single facet from n-1 dimension points
      stringstream message;
      message << "Hull dimension exceeds the number of input points!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }

    //update facets here:
    this->m_facets.clear();
    NhullFacet::qhullRunToNhullFacets(m_facets, qhull, m_points);
  }

  ///@brief calculates a thermodynamic hull: takes a list of aflowlib entries, determines facets, and updates thermo quantitities
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::calculateThermoHull() {
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    stringstream message;
    vector<Entry> pre_flag_points;

    //retrieve all flags that modify run
    const bool do_IQR = true;
    const bool read_cache = m_nflags.flag(convex_flags.read_cache);
    const bool write_cache = m_nflags.flag(convex_flags.write_cache);
    const bool calc_nminus1 = m_nflags.flag(convex_flags.m_nminus1_flag);
    const string nminus1_auid = m_nflags.getattachedscheme(convex_flags.m_nminus1_flag);
    const bool calc_stability_criterion = m_nflags.flag(convex_flags.m_stability_criterion_flag);
    const string stability_auid = m_nflags.getattachedscheme(convex_flags.m_stability_criterion_flag);
    const bool calc_stats = m_nflags.flag(convex_flags.m_stats);
    const bool print = m_nflags.flag(convex_flags.m_nhull_print);
    const bool print_plots = m_nflags.flag(convex_flags.m_latex_doc);
    const bool print_JSON = m_nflags.flag(convex_flags.m_print_JSON);
    const bool get_dhull = m_nflags.flag(convex_flags.m_dhull);
    const string dhull_auid = m_nflags.getattachedscheme(convex_flags.m_dhull);

    //if user entered auid for stability calc check that this exists in hull first
    if (calc_stability_criterion && !Entry::findAUID(stability_auid, m_points)) {
      stringstream message;
      message << "Specified auid not found on hull (auid=" << stability_auid << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }

    //if we are doing sc or nminus, save a copy of m_points before hull point flags are changed
    if (calc_stability_criterion || calc_nminus1) {
      pre_flag_points = m_points;
    }

    //filter points
    const vector<Entry> filtered_entries = filterEntries(m_points, m_velements, do_IQR, false, m_dimension);

    //check filtered entries
    if (filtered_entries.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No entries could pass filtering for hull calculation!");
    }

    //qhull requires dimension+1 points ALWAYS:
    if (filtered_entries.size() < m_dimension + 1) {
      stringstream message;
      message << "Not enough entries filtered to create a " << std::to_string(m_dimension) << " dimensional qhull!" << "\n";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }

    //run hull calc, update points, and return facets
    m_facets = iterateRunFacetsOnly(m_points, m_velements, m_dimension, filtered_entries);

    //now initialize all distance to hulls here:
    postIterateRunInitializeEntries(m_points, m_facets);

    //run other hulls if they were asked for
    if (calc_nminus1) {
      //updates nminus1 value in ENTRY
      std::pair nminus1_calc = {calc_nminus1, nminus1_auid};
      const vector<Entry> nminus1_filtered_entries = filterEntries(pre_flag_points, m_velements, do_IQR, false, m_dimension, nminus1_calc, {false, ""});
      vector<NhullFacet> nminus1_facets = calculateNminus1(pre_flag_points, m_velements, m_dimension, nminus1_filtered_entries, nminus1_auid);
      if (print_plots && print) {
        string nminus1_plot_name = specieString(m_velements) + "_nminus1";
        savePlotAndTable(nminus1_plot_name, m_points, m_velements, nminus1_facets, m_aflags, m_nflags);
      }
    }

    if (calc_stability_criterion) {
      std::pair sc_calc = {calc_stability_criterion, stability_auid};
      const vector<Entry> sc_filtered_entries = filterEntries(pre_flag_points, m_velements, do_IQR, false, m_dimension, {false, ""}, sc_calc);
      vector<NhullFacet> sc_facets = calculateStabilityCriterion(pre_flag_points, m_velements, m_dimension, sc_filtered_entries, stability_auid);
      if (print_plots && print) {
        string sc_plot_name = specieString(m_velements) + "_stability_criterion";
        savePlotAndTable(sc_plot_name, m_points, m_velements, sc_facets, m_aflags, m_nflags);
      }
    }

    if (get_dhull) {
      distanceToHullPrintToConsole(dhull_auid);
    }

    //save full-dimension hull, table, and plots here:
    if (print) {
      if (print_JSON) {
        saveHullJSON(specieString(m_velements));
      }
      if (print_plots) {
        savePlotAndTable(specieString(m_velements), m_points, m_velements, m_facets, m_aflags, m_nflags);
      }
    }

    if (calc_stats) {
      compoundStats(lightFilter(m_points), m_aflags.Directory);
    }
  }

  ///@brief after calculateThermoHull has been run, this function determines distances to hull and phase decompositions
  ///@param entries global list of entries
  ///@param facets hull facets
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::postIterateRunInitializeEntries(vector<Entry>& entries, vector<NhullFacet>& facets) {
    vector<Entry> hull_entries;
    for (const Entry& entry : entries) {
      if (entry.is_hull_point) {
        hull_entries.push_back(entry);
      }
    }

    for (Entry& entry : entries) {
      //hull points have already been initialized
      if (!entry.is_hull_point) {
        NhullFacet::calc_vertical_distance_and_phase_decomposition(entry, hull_entries, facets);
      }
    }
  }

  ///@brief determines the distance to hull of a point at dimension n by removing all points with dimension (m<n)
  ///@param pre_flag_entries a clean copy of the global entries list before we start modifying entries here
  ///@param filtered_entries list of all entries deemed acceptable to use for hull generation
  ///@param nminus1_auid
  ///@param verbose
  ///@return list of facets
  /// @authors
  /// mod{NHA,20260206,created}
  vector<NhullFacet> ConvexHull::calculateNminus1(vector<Entry> pre_flag_entries, const vector<string> specie_vector, const int dimension, const vector<Entry>& filtered_entries, const string& nminus1_auid, bool verbose) {
    //assumes all hull points have been initialized and the full dimension hull calculation has completed
    //make a copy of the global entries (entries) and update here with n-1 values
    vector<NhullFacet> nhull_facets;
    uint dim_nminus1 = m_dimension - 1;

    Entry nminus1_target = Entry::getEntryFromAuid(pre_flag_entries, nminus1_auid);

    //first check trivial case
    if (m_dimension == 2) {
      for (Entry& entry : m_points) {
        if (entry.nspecies) {
          entry.m_nminus1_enthalpy_gain = entry.distance_hull_enthalpy_formation_atom;
        }
      }
    } else {
      nhull_facets = iterateRunFacetsOnly(pre_flag_entries, specie_vector, dimension, filtered_entries);
      postIterateRunInitializeEntries(pre_flag_entries, nhull_facets);

      //now update m_nminus1_enthalpy_gain globally
      for (const auto& entry : pre_flag_entries) {
        m_points[entry.nh_index].m_nminus1_enthalpy_gain = entry.distance_hull_enthalpy_formation_atom;
        if (entry.auid == nminus1_auid) {
          nminus1_target = m_points[entry.nh_index];
        }
      }
    }

    if (verbose) {
      stringstream return_stream;
      return_stream << "Nminus1 analysis for " << nminus1_target.auid << "= " << nminus1_target.m_nminus1_enthalpy_gain * 1000;
      std::cout << return_stream.str() << std::endl;
    }

    //return facets if we want to plot n-1 hulls
    return nhull_facets;
  }

  ///@brief determines the distance to hull of a point at dimension n by removing all points with the same composition as the given auid
  ///@param pre_flag_entries a clean copy of the global entries list before we start modifying entries here
  ///@param filtered_entries list of all entries deemed acceptable to use for hull generation
  ///@param auid auid of the compound to be excluded from hull creation
  ///@param verbose if true, outputs results to console
  ///@return list of facets
  /// @authors
  /// mod{NHA,20260206,created}
  vector<NhullFacet> ConvexHull::calculateStabilityCriterion(vector<Entry> preflag_entries, const vector<string>& specie_vector, const int dimension, const vector<Entry>& filtered_entries, const string& auid, bool verbose) {
    //assumes all hull points have been initialized and the full dimension hull calculation has completed
    vector<NhullFacet> nhull_facets;
    Entry entropy_stabilization_target = Entry::getEntryFromAuid(preflag_entries, auid);
    auto compound_map = entropy_stabilization_target.compound_map;
    nhull_facets = iterateRunFacetsOnly(preflag_entries, specie_vector, dimension, filtered_entries);
    postIterateRunInitializeEntries(preflag_entries, nhull_facets);

    //now update the stability criterion for the chosen auid and all others of the same compound:
    for (Entry& entry : preflag_entries) {
      if (entry.compound_map == compound_map) {
        m_points[entry.nh_index].m_stability_criterion = entry.distance_hull_enthalpy_formation_atom;
        m_points[entry.nh_index].relative_stability_criterion = entry.distance_hull_enthalpy_formation_atom / fabs(entry.enthalpy_formation_atom);
      }
      if (entry.auid == auid) {
        entropy_stabilization_target = m_points[entry.nh_index];
      }
    }

    if (verbose) {
      stringstream return_stream;
      return_stream << "Stability Criterion analysis for " << entropy_stabilization_target.auid << "= " << entropy_stabilization_target.m_stability_criterion * 1000;
      std::cout << return_stream.str() << std::endl;
    }

    //return facets if we want to plot stability_criterion hulls
    return nhull_facets;
  }

  /// @brief Hull calculation routine for POCC directories. Calculates EFA and DEED values
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::calculatePOCC() {
    //assumes all hull points have been initialized and the full dimension hull calculation has completed
    stringstream message;
    vector<Entry> pre_flag_points;

    //retrieve all flags
    const bool do_IQR = true;
    const bool calc_stats = m_nflags.flag(convex_flags.m_stats);
    const string pocc_path = m_nflags.getattachedscheme(convex_flags.m_pocc);
    const bool print = m_nflags.flag(convex_flags.m_nhull_print);
    const bool print_plots = m_nflags.flag(convex_flags.m_latex_doc);
    const bool print_JSON = m_nflags.flag(convex_flags.m_print_JSON);

    // get pocc stuff: DG's and associated AUID's
    vector<pair<string, unsigned long long int>> dg_list = loadPocc(pocc_path);

    //filter points
    const vector<Entry> filtered_entries = filterEntries(m_points, m_velements, do_IQR, true, m_dimension);

    //check filtered entries
    if (filtered_entries.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No entries could pass filtering for hull calculation!");
    }

    if (filtered_entries.size() < m_dimension + 1) {
      stringstream message;
      message << "Not enough entries filtered to create a " << std::to_string(m_dimension) << " dimensional hull!" << "\n";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }

    //run hull calc, update points, and return facets
    m_facets = iterateRunFacetsOnly(m_points, m_velements, m_dimension, filtered_entries);

    //now initialize all distance to hulls here:
    postIterateRunInitializeEntries(m_points, m_facets);

    //save full-dimension hull, table, and plots here:
    if (print) {
      if (print_JSON) {
        saveHullJSON(specieString(m_velements));
      }
      if (print_plots) {
        savePlotAndTable(specieString(m_velements), m_points, m_velements, m_facets, m_aflags, m_nflags);
      }
    }

    //compile statistics
    singlePoccStats(m_points, dg_list);
  }

  ///@brief Simplest convex hull generator for general non-thermodynamic calculations
  ///@param vcoords list of input points
  ///@return true if no errors were thrown
  /// @authors
  /// mod{NHA,20260206,created}
  bool ConvexHull::createHull(const vector<xvector<double>>& vcoords) {
    try {
      initializePointsSimple(vcoords);
      calculateSimpleHull();
      m_initialized = true;
      // hull must be initialized for these analyses
      //thermodynamicsPostProcessing(); // will return if not m_thermo_hull
    } catch (aurostd::xerror& err) {
      ostream& oss = cout;
      ofstream FileMESSAGE;
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
    }
    return m_initialized;
  }

  /// @brief convex hull generator for thermodynamic calculations (energy distance to hull, phase decomp, n-1, stability criterion)
  /// @return true if no errors thrown
  /// @authors
  /// mod{NHA,20260206,created}
  bool ConvexHull::createThermoHull() {
    try {
      calculateThermoHull();
      m_initialized = true;
      // hull must be initialized for these analyses
    } catch (aurostd::xerror& err) {
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    return m_initialized;
  }

  ///@brief Initializes points for thermodynamic hull calcs (dhulls, phase decomp, n-1 enthalpy, sc): Assumes points are non-stoichiometric (first index removed)
  ///@param vpoints list of entries to evaluate
  ///@param add_artificial_unaries if true, adds n points for n-dimension hull representing pure elements at zero enthalpy
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::initializePointsThermo(const std::vector<Entry>& vpoints, bool add_artificial_unaries) {
    //NHULL:: this function was copied from loadPoints() from chull (ln# 5134)
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    stringstream message;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " initializing Entrys, count=" << vpoints.size() << endl;
    }
    for (size_t i = 0, fl_size_i = vpoints.size(); i < fl_size_i; i++) {
      m_points.push_back(vpoints[i]);
    };
    if (m_points.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No points loaded, no way to determine dimensionality of hull");
    }

    m_dimension = m_points[0].m_coord.rows;
    if (m_dimension < 2) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "1D hulls are trivial");
    } // this MUST be true: nhull cannot find facets for 1D hulls
    // test of stupidity
    for (size_t i = 0, fl_size_i = m_points.size(); i < fl_size_i; i++) {
      if (m_points[i].getCoordDimension() != m_dimension) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Dimension mismatch among points");
      }
    }

    // ensures proper construction of hull
    // duplicates DO NOT AFFECT performance/accuracy of algorithm
    // we ignore these points after hull construction anyway
    if (add_artificial_unaries) {
      addArtificialUnaries();
    }
  }

  /// @brief Function meant to work with legacy CHULL datastructures for running GFA
  /// @note to avoid rewriting GFA module we needed to import a chull datastructure: CoordGroups
  /// @note CoordGroups HAS to be initialized after finding outliers (IQR)
  /// @note code taken from chull.cpp ln# 5629 structurePoints(). This file exists in versions of AFLOW <4.1
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::initializeCoordGroups() {
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    stringstream message;
    uint i_coord_group_sort; // so it doesn't conflict with i_coord_group in for-loops
    CoordGroup cg;
    xvector<double> r_coords;
    bool perform_outliers_analysis = m_nflags.flag(convex_flags.m_interquartile_range);

    for (size_t i = 0, fl_size_i = m_points.size(); i < fl_size_i; i++) {
      const Entry& point = m_points[i];
      r_coords = point.getStoichiometricCoords();
      if (!setCoordGroupIndex(r_coords, i_coord_group_sort)) {
        cg.setCoords(r_coords, true);
        m_coord_groups.push_back(cg);
        i_coord_group_sort = m_coord_groups.size() - 1;
      }
      m_coord_groups[i_coord_group_sort].m_points.push_back(i);
      if (m_coord_groups[i_coord_group_sort].m_has_stoich_coords && false) {
        message << "Mismatch among coord types (stoich vs. non-stoich coords), assuming non-stoich coords";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        m_coord_groups[i_coord_group_sort].m_has_stoich_coords = false;
      }
      if (point.isUnary() && point.m_is_artificial) {
        m_coord_groups[i_coord_group_sort].m_has_artificial_unary = true;
      }
    }

    // remove outliers before sort (BY MATCHING COORDS)
    // we do this analysis EVEN if we keep outliers in the end, simply print out useful warnings for user
    uint i_point = AUROSTD_MAX_UINT;
    vector<uint> outliers;

    if (perform_outliers_analysis) {
      vector<size_t> allowed_entry_indices;
      for (int i = 0; i < m_points.size(); i++) {
        if (!m_points[i].flags.is_flagged(false)) {
          allowed_entry_indices.push_back(i);
        }
      }

      xvector<int> elements_present;
      vector<uint> tmp_store;
      vector<uint> outliers_found;
      //get all binary combos:
      aurostd::xcombos xc(m_dimension, 2, 'C');
      while (xc.increment()) {
        elements_present = aurostd::vector2xvector<int>(xc.getCombo());
        tmp_store = getOutliers(m_points, allowed_entry_indices, elements_present, ConvexHull::isHalfHull(m_points, false));
        outliers_found.insert(outliers_found.end(), tmp_store.begin(), tmp_store.end());
      }
      outliers = outliers_found;
    }

    //remove outliers from coord group:
    bool found_outlier = false;
    vector<uint> points_to_remove;
    uint valid_count = 0;
    for (size_t i_coord_group = 0, fl_size_i_coord_group = m_coord_groups.size(); i_coord_group < fl_size_i_coord_group; i_coord_group++) {
      points_to_remove.clear();
      for (size_t i = 0, fl_size_i = m_coord_groups[i_coord_group].m_points.size(); i < fl_size_i; i++) {
        i_point = m_coord_groups[i_coord_group].m_points[i];

        found_outlier = false;
        for (size_t j = 0, fl_size_j = outliers.size(); j < fl_size_j && !found_outlier; j++) {
          if (i_point == outliers[j]) {
            found_outlier = true;
          }
        }
        if (found_outlier) {
          points_to_remove.push_back(i);
        } // not i_point, so I can remove this index
        else {
          if (!m_points[i_point].m_is_artificial) {
            valid_count++;
          }
        }
      }
      if (!points_to_remove.empty()) {
        std::sort(points_to_remove.rbegin(), points_to_remove.rend()); // descending
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " before outlier removal, count = " << m_coord_groups[i_coord_group].m_points.size() << endl;
        }
        for (size_t i = 0, fl_size_i = points_to_remove.size(); i < fl_size_i; i++) {
          i_point = m_coord_groups[i_coord_group].m_points[points_to_remove[i]];
          if (!m_points[i_point].m_initialized) {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Point[" + aurostd::utype2string(i_point) + "] is not initialized");
          }
          message << "Removing outlier ";
          if (m_points[i_point].m_has_entry) {
            message << "auid=" << m_points[i_point].auid;
          } else {
            message << "m_coord=" << m_points[i_point].m_coord;
          }
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_OPTION_);
          m_coord_groups[i_coord_group].m_points.erase(m_coord_groups[i_coord_group].m_points.begin() + points_to_remove[i]);
        }
      }
    }
    message << "Employing " << valid_count << " total entries for " << pflow::arity_string(m_dimension, false, false) << " convex hull analysis";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);

    // remove empty m_coord_groups
    vector<uint> empty_coord_groups;
    for (size_t i_coord_group = 0, fl_size_i_coord_group = m_coord_groups.size(); i_coord_group < fl_size_i_coord_group; i_coord_group++) {
      if (m_coord_groups[i_coord_group].m_points.empty()) {
        empty_coord_groups.push_back(i_coord_group);
      }
    }
    std::sort(empty_coord_groups.rbegin(), empty_coord_groups.rend()); // descending
    for (size_t i = 0, fl_size_i = empty_coord_groups.size(); i < fl_size_i; i++) {
      m_coord_groups.erase(m_coord_groups.begin() + empty_coord_groups[i]);
    }

    // sort
    for (size_t i_coord_group = 0, fl_size_i_coord_group = m_coord_groups.size(); i_coord_group < fl_size_i_coord_group; i_coord_group++) {
      std::sort(m_coord_groups[i_coord_group].m_points.begin(), m_coord_groups[i_coord_group].m_points.end(), sortWithinCoordGroup(m_points, true)); // ascending order
    }
    std::sort(m_coord_groups.rbegin(), m_coord_groups.rend()); // descending order for alphabetic print out later

    // assign coord group indices to points, useful later
    for (size_t i_coord_group = 0, fl_size_i_coord_group = m_coord_groups.size(); i_coord_group < fl_size_i_coord_group; i_coord_group++) {
      // set ref state
      if (!m_coord_groups[i_coord_group].m_points.empty()) {
        m_coord_groups[i_coord_group].m_ref_state = artificialMap(m_coord_groups[i_coord_group].m_points[0]);
      }
      for (size_t i = 0, fl_size_i = m_coord_groups[i_coord_group].m_points.size(); i < fl_size_i; i++) {
        i_point = m_coord_groups[i_coord_group].m_points[i];
        if (!m_points[i_point].m_initialized) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Point[" + aurostd::utype2string(i_point) + "] is not initialized");
        }
        m_points[i_point].m_i_coord_group = i_coord_group;
      }
    }
  }

  /// @brief Utility function used to initialize CoordGroups member var
  /// @brief Takes the first point (in order) of a point in m_coord_groups that is not artificial
  /// @brief If only artificial points exist, return the first one
  /// @return returns an index corresponding to the global entries vector (m_points)
  /// @authors
  /// mod{NHA,20260206,created}
  uint ConvexHull::artificialMap(uint i_point) const {
    if (i_point > m_points.size() - 1) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Invalid index within m_points");
    }
    if (!m_points[i_point].m_is_artificial) {
      return i_point;
    }

    const Entry& art_point = m_points[i_point];
    uint i_coord_group = AUROSTD_MAX_UINT;
    if (!setCoordGroupIndex(art_point.m_coord, i_coord_group)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Coordgroup index not set");
    }
    if (!m_coord_groups[i_coord_group].m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Uninitialized coordgroup");
    }
    if (m_coord_groups[i_coord_group].m_points.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Coordgroup[" + aurostd::utype2string(i_coord_group) + "] has no points");
    }

    uint i_point_new = AUROSTD_MAX_UINT;
    for (size_t i = 0, fl_size_i = m_coord_groups[i_coord_group].m_points.size(); i < fl_size_i; i++) {
      i_point_new = m_coord_groups[i_coord_group].m_points[i];
      if (!m_points[i_point_new].m_is_artificial) {
        return i_point_new;
      }
    }

    // if we get here, there are no viable replacements, simply return artificial point
    if (m_coord_groups[i_coord_group].m_points.size() == 1) {
      if (m_coord_groups[i_coord_group].m_points[0] == i_point) {
        return i_point;
      }
    }

    // really bad if we get here
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot determine artificial point mapping");
  }

  ///@brief default constructor for the hull
  /// @authors
  /// mod{NHA,20260206,created}
  ConvexHull::ConvexHull() {
    m_initialized = false; // no points
  }

  ///@brief simplest hull constructor for non-thermodynamic general hull calcs
  ///@param hull_options contains user-specified settings for running this hull
  ///@param points list of points to create hull with
  /// @authors
  /// mod{NHA,20260206,created}
  ConvexHull::ConvexHull(aurostd::xoption hull_options, vector<aurostd::xvector<double>> points) {
    try {
      setHullFlags(hull_options);
      setDirectory(hull_options);
      m_initialized = createHull(points);
    } catch (aurostd::xerror& err) {
      ostream& oss = cout;
      ofstream FileMESSAGE;
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
    }
  }

  ///@brief thermodynamic hull constructor used only for GFA
  ///@param vpflow run options
  ///@param vpoints hull points
  ///@param oss
  ///@param formation_energy_hull legacy argument needed for GFA (indicates that this is a thermodynamic hull)
  ///@param add_artificial_unaries if true, adds n points for n-dimension hull representing pure elements at zero enthalpy
  /// @authors
  /// mod{NHA,20260206,created}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow, const vector<Entry>& vpoints, ostream& oss, bool formation_energy_hull, bool add_artificial_unaries) : xStream(oss), m_initialized(false) {
    try {
      setHullFlags(vpflow);
      setDirectory(vpflow);
      initializePointsThermo(vpoints, add_artificial_unaries);
      m_initialized = createThermoHull();
      initializeCoordGroups();
    } catch (aurostd::xerror& err) {
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
  }

  ///@brief hull constructor used for pocc or for a regular thermodynamic hull calculation
  ///@param vpflow run options
  ///@param aflags directory information
  ///@param specie_vector global compound at dimension n
  ///@param add_artificial_unaries if true, adds n points for n-dimension hull representing pure elements at zero enthalpy
  /// @authors
  /// mod{NHA,20260206,created}
  ConvexHull::ConvexHull(const aurostd::xoption& vpflow, const _aflags& aflags, vector<string> specie_vector, ostream& oss, bool add_artificial_unaries) :
      xStream(oss), m_initialized(false), m_aflags(aflags), m_velements(specie_vector) {
    try {
      setHullFlags(vpflow);
      setDirectory(vpflow);
      vector<Entry> vpoints = getEntriesFromAFLOWlib(specie_vector, vpflow);
      initializePointsThermo(vpoints, add_artificial_unaries);
      m_initialized = true;
    } catch (aurostd::xerror& err) {
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
  }

  //--------------------------------------------------------------------------------
  // Coord Group Implementations
  //--------------------------------------------------------------------------------

  /// @brief default initializer for CoordGroup class
  /// @note this is a legacy function taken from chull for use in GFA function
  /// @authors
  /// mod{NHA,20260206,created}
  CoordGroup::CoordGroup() :
      m_initialized(false),
      m_has_stoich_coords(false),
      m_has_artificial_unary(false),
      m_is_on_hull(false),
      m_hull_member(0),
      m_ref_state(0),
      m_i_nary(0),
      m_i_alloy(0),
      m_nearest_facet(0),
      m_nearest_distance(0),
      m_calculated_equivalent_g_states(false),
      m_stability_criterion(0),
      m_n_plus_1_enthalpy_gain(0),
      m_found_icsd_g_state(false),
      m_i_icsd_g_state(0) {};

  ///@brief operator defined to allow sorting of CoordGroup
  bool CoordGroup::operator<(const CoordGroup& other) const {
    // safety, so it doesn't break, but it's outside scope of function
    if (m_coord.rows != other.m_coord.rows) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Dimension mismatch between stoichiometries");
    } // return (m_coord.rows<other.m_coord.rows);
    for (int i = m_coord.lrows; i <= m_coord.urows; i++) {
      if (m_coord(i) != other.m_coord(i)) {
        return (m_coord(i) < other.m_coord(i));
      }
    }
    return false;
  }

  ///@brief assigns the element coordinates of a CoordGroup object
  ///@param has_stoich_coords if true, then m_coord represents a coordinate with the first element index included
  /// @authors
  /// mod{NHA,20260206,created}
  void CoordGroup::setCoords(const xvector<double>& coord, bool has_stoich_coords) {
    m_coord = coord;
    m_has_stoich_coords = has_stoich_coords;
    m_has_artificial_unary = false;
    m_initialized = true;
  }

  //***************** END OF CLASS COORDGROUP

  /// @brief substitutes pure elements in calculation with formation_enthalpy_atom = 0 "dummy" entries
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::addArtificialUnaries() {
    for (uint i = 0; i < m_dimension; i++) {
      // points are really dimension+1 dimensional (hidden dimension)
      xvector<double> coordinates(m_dimension);
      if (i != 0) {
        coordinates[i + coordinates.lrows - 1] = 1;
      }

      Entry entry(coordinates, true);
      entry.enthalpy_formation_atom = 0;
      entry.set_distance_hull_enthalpy_formation_atom(0);
      entry.is_hull_point = true;
      entry.m_is_artificial = true;
      entry.nspecies = 1;
      entry.lowest_energy_alloy = true;
      entry.compound = m_velements[i] + "1";
      entry.compound_map = {
          {m_velements[i], 1}
      };
      entry.nh_index = m_points.size();
      entry.flagged_entry = false;
      m_points.push_back(entry);
    }
  }

  /// @brief returns distance to hull (vertical) of a point. Used only for GFA
  /// @authors
  /// mod{NHA,20260206,created}
  double ConvexHull::getDistanceToHull(Entry& entry, bool redo, bool get_signed_distance) const {
    if (!entry.m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "entry has not been initialized!");
    }
    return entry.distance_hull_enthalpy_formation_atom;
  }

  /// @brief returns distance to hull (vertical) of a point initialized outside of current hull calculation.
  /// @returns vertical distance to hull facet
  /// @authors
  /// mod{NHA,20260206,created}
  double ConvexHull::getDistanceToHull(const string auid) const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Hull has not been initialized!");
    }
    Entry dhull_target = Entry::getEntryFromAuid(m_points, auid);
    return dhull_target.distance_hull_enthalpy_formation_atom;
  }

  /// @brief returns m_nminus1_enthalpy_gain given an auid
  /// @authors
  /// mod{NHA,20260206,created}
  double ConvexHull::getNminus1EnthalpyGain(const std::string& auid) const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Hull has not been initialized!");
    }

    Entry nminus1_target = Entry::getEntryFromAuid(m_points, auid);
    return nminus1_target.m_nminus1_enthalpy_gain;
  }

  ///@brief returns total number of ingested hull points (flagged entries included)
  /// @authors
  /// mod{NHA,20260206,created}
  uint ConvexHull::getEntriesCount() const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Hull has not been initialized!");
    }
    return m_points.size();
  }

  ///@brief prints distance to hull of selected auid to console
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::distanceToHullPrintToConsole(const string auid) const {
    stringstream message;
    message << "Distance to hull for " << auid << " is: ";
    message << getDistanceToHull(auid) * 1000 << " meV.";
    std::cout << message.str() << std::endl;
  }

  /// @brief returns list of ch_indices of entries associated with the decomposition of the given entry
  /// @authors
  /// mod{NHA,20260206,created}
  std::vector<uint> ConvexHull::getDecompositionPhases(Entry& entry) const {
    if (!entry.m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "entry has not been initialized!");
    }
    //convert to nh_index here:
    vector<uint> indices;
    for (const auto& [entry, dhull] : entry.m_nhull_phase_decomp_entries) {
      indices.push_back(entry.nh_index);
    }

    return indices;
  }

  /// @brief returns stability criterion of an auid
  /// @authors
  /// mod{NHA,20260206,created}
  double ConvexHull::getStabilityCriterion(const std::string& auid) const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Hull has not been initialized!");
    }

    Entry target = Entry::getEntryFromAuid(m_points, auid);
    return target.m_stability_criterion;
  }

  /// @brief returns decomposition ratios of global entry with index i_point
  /// @authors
  /// mod{NHA,20260206,created}
  aurostd::xvector<double> ConvexHull::getDecompositionCoefficients(uint i_point) const {
    if (!m_initialized) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Hull has not been initialized!");
    }

    aurostd::xvector<double> ratios = m_points[i_point].findPhaseDecompRatios();
    return ratios;
  }

  /// @brief converts enthalpy_formation_atom units from eVs to meVs
  /// @return list of converted entries
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> ConvexHull::convertUnitsAllEntries(const vector<Entry>& entries) {
    vector<Entry> converted_entries;
    for (Entry entry : entries) {
      entry.enthalpyUnitConversion();
      converted_entries.push_back(entry);
    }
    if (converted_entries.size() != entries.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Failed to convert all entries!");
    }
    return converted_entries;
  }

  /// @brief converts enthalpy_formation_atom units from eVs to meVs for given facets
  /// @returns list of converted facets
  /// @authors
  /// mod{NHA,20260206,created}
  vector<NhullFacet> ConvexHull::convertUnitsAllFacets(const vector<NhullFacet>& facets) {
    vector<NhullFacet> converted_facets;
    for (NhullFacet facet : facets) {
      facet.enthalpyUnitConversion();
      converted_facets.push_back(facet);
    }
    if (converted_facets.size() != facets.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Failed to convert all facets!");
    }
    return converted_facets;
  }

  /// @brief saves all hull points to JSON
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::saveHullJSON(const string& file_name) const {
    if (m_initialized) {
      string directory = m_aflags.Directory;
      const string json_path = directory + "/" + file_name + "_hull.json";
      //convert units before saving:
      std::vector<Entry> converted_points = convertUnitsAllEntries(m_points);
      aurostd::JSON::object(converted_points).saveFile(json_path);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Hull not initialized");
    }
  }

  /// @brief saves a pdf plot and table from a list of entries
  /// @param aflags directory information
  /// @param nflags user-specified run flags
  /// @authors
  /// mod{NHA,20260206,created}
  void ConvexHull::savePlotAndTable(const string& file_name, const vector<Entry>& entries, const vector<string>& specie_vector, const vector<NhullFacet>& facets, const _aflags aflags, const aurostd::xoption nflags) {
    string directory = aflags.Directory;
    const string save_root = directory + "/" + file_name;
    bool plot_3D = nflags.flag("NHULL::PRINT3D");
    //convert units before saving:
    vector<Entry> converted_points = convertUnitsAllEntries(entries);
    vector<NhullFacet> converted_facets = convertUnitsAllFacets(facets);
    nhullxplotter::toTexXplotter(converted_points, specie_vector, converted_facets, save_root, plot_3D);
  }

  ///brief returns all entries associated with a string alloy: i.e all m <= n dimensional alloys with same species
  /// @brief brief returns all entries associated with a string alloy: i.e all m <= n dimensional alloys with same species
  /// @note points are loaded and sorted as to optimize outlier identification later on
  /// @return list of Entry objects
  /// @authors
  /// mod{NHA,20260206,created}
  vector<Entry> ConvexHull::getEntriesFromAFLOWlib(vector<string> specie_vector, const aurostd::xoption& vpflow) {
    //Loadpoints overload on 4985 (#1)
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);

    //extract whether or not to do cce
    bool neglect_cce = vpflow.flag(convex_flags.m_neglect_cce);

    //make sure in same order (alphabetic) each time:
    std::sort(specie_vector.begin(), specie_vector.end());

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " initializing velements, velements=" << aurostd::joinWDelimiter(specie_vector, ",") << endl;
    }

    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Loading raw entries from the database", m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<vector<vector<aflowlib::_aflowlib_entry>>> entries;
    // HE new method
    {
      // ensures that all resources are released as soon as possible
      aflowlib::EntryLoader el;
      el.m_out_debug = LDEBUG; // leave this on for testing
      el.loadAlloy(specie_vector);
      el.getEntriesThreeLayer(entries); // the alloy will be loaded recursive by default, (elements will be cleaned and sorted)

      if (entries.empty()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No entries could be loaded from the database.");
      }
    }
    const uint i_nary = 0;
    uint i_alloy = 0;
    const bool sortLIB2 = true;
    if (sortLIB2) {
      uint i = 0;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " reshuffling unaries to optimize duplicate removal" << endl;
      }
      vector<_aflowlib_entry_LIB2sorting> vaes;
      uint i_entry = 0;
      for (i_alloy = 0; i_alloy < entries[i_nary].size(); i_alloy++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " i_alloy=" << i_alloy << endl;
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " INITIAL ORDER=" << endl;
          for (i = 0; i < entries[i_nary][i_alloy].size() - 1; i++) {
            cerr << __AFLOW_FUNC__ << " i_entry.auid=" << entries[i_nary][i_alloy][i].auid << " i_entry.aurl=" << entries[i_nary][i_alloy][i].aurl << endl;
          }
          cerr << endl;
        }
        vaes.clear();
        for (i_entry = 0; i_entry < entries[i_nary][i_alloy].size(); i_entry++) {
          vaes.emplace_back(entries[i_nary][i_alloy][i_entry], i_entry, specie_vector, *p_FileMESSAGE, *p_oss);
        }
        std::stable_sort(vaes.begin(), vaes.end());
        vector<uint> order_new;
        for (i_entry = 0; i_entry < vaes.size(); i_entry++) {
          order_new.push_back(vaes[i_entry].m_index);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " order_new=" << aurostd::joinWDelimiter(order_new, ",") << endl;
        }
        aurostd::reorder(entries[i_nary][i_alloy], order_new, 1); // use mode==1
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " FINAL ORDER=" << endl;
          for (i = 0; i < entries[i_nary][i_alloy].size() - 1; i++) {
            cerr << __AFLOW_FUNC__ << " i_entry.auid=" << entries[i_nary][i_alloy][i].auid << " i_entry.aurl=" << entries[i_nary][i_alloy][i].aurl << endl;
          }
          cerr << endl;
        }
      }
    }
    size_t count = 0;
    for (size_t i = 0; i < entries.size(); i++) {
      for (size_t j = 0; j < entries[i].size(); j++) {
        count += entries[i][j].size();
      }
    }
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Loaded " + aurostd::utype2string(count) + " raw entries from the database", m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " initializing velements WITH entries, velements=" << aurostd::joinWDelimiter(specie_vector, ",") << endl;
    }

    vector<Entry> points;
    ofstream FileMESSAGE;
    int index = 0;
    for (size_t i = 0, fl_size_i = entries.size(); i < fl_size_i; i++) {
      for (size_t j = 0, fl_size_j = entries[i].size(); j < fl_size_j; j++) {
        for (size_t k = 0, fl_size_k = entries[i][j].size(); k < fl_size_k; k++) {
          Entry cp(entries[i][j][k], specie_vector, index, neglect_cce);
          index++;
          points.push_back(cp);
        }
      }
    }
    if (points.empty()) {
      // simply always die here, we cannot grab dimensionality of hull without ANY points
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No entries loaded");
    }
    return points;
  }

  /// @brief sets the indexing of a coord group given a coordinate (r_coord) that it shares
  /// @param r_coords coordinate to do matching with
  /// @param i_coord_group index given to a coordgroup
  /// @return true if r_coords is in this coord_group, false otherwise
  /// @authors
  /// mod{NHA,20260206,created}
  bool ConvexHull::setCoordGroupIndex(const aurostd::xvector<double>& r_coords, uint& i_coord_group) const {
    const bool LDEBUG = (false || NHULL_DEBUG || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " r_coords=" << r_coords << endl;
      cerr << __AFLOW_FUNC__ << " m_coord_groups.size()=" << m_coord_groups.size() << endl;
      for (size_t i = 0, fl_size_i = m_coord_groups.size(); i < fl_size_i; i++) {
        cerr << __AFLOW_FUNC__ << " m_coord_groups[i=" << i << "].m_coord=" << m_coord_groups[i].m_coord << endl;
      }
    }
    for (size_t i = 0, fl_size_i = m_coord_groups.size(); i < fl_size_i; i++) {
      if (!m_coord_groups[i].m_initialized) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Uninitialized coordgroup");
      }
      if (aurostd::identical(m_coord_groups[i].m_coord, r_coords)) {
        i_coord_group = i;
        return true;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " m_coord_groups[i=" << i << "].m_coord=" << m_coord_groups[i].m_coord << " != " << "r_coords=" << r_coords << endl;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " no coord_group index found for r_coords=" << r_coords << endl;
    }
    return false;
  }

  ///@brief sorting operator for CoordGroups
  ///@returns true if the CoordGroup at index j is "greater than" the one at i
  bool sortWithinCoordGroup::sortWithinCoordGroup::operator()(uint i, uint j) {
    // ascending order
    if ((i > m_points.size() - 1) || (j > m_points.size() - 1)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Invalid index within CoordGroup");
    } // safety
    const Entry& ci = m_points[i];
    const Entry& cj = m_points[j];
    if (ci.isGState() != cj.isGState()) {
      return ci.isGState() > cj.isGState();
    }
    // always sort by last coord first
    const bool energy_sort = (m_sort_energy_ascending ? (ci.getLastCoord() < cj.getLastCoord()) : (ci.getLastCoord() > cj.getLastCoord()));
    if (ci.m_has_entry && cj.m_has_entry) {
      // if entry, then we also sort by proto, compound, and then aurl to make final print out pretty
      return (energy_sort
              || ((ci.getLastCoord() == cj.getLastCoord()) && (ci.prototype < cj.prototype))
              || ((ci.getLastCoord() == cj.getLastCoord()) && (ci.prototype == cj.prototype) && (ci.compound < cj.compound))
              || ((ci.getLastCoord() == cj.getLastCoord()) && (ci.prototype == cj.prototype) && (ci.compound == cj.compound) && (ci.aurl < cj.aurl)));
      // aurl is guaranteed to be unique (more so than auid)
    } else {
      return energy_sort;
    }
  }

  ///@brief returns true if there are candidate hull points that lie below zero-energy tie line
  ///@param true if this hull is being generated for a POCC directory
  /// @authors
  /// mod{NHA,20260206,created}
  bool ConvexHull::isHalfHull(const vector<Entry>& entries, bool pocc_run) {
    for (const Entry entry : entries) {
      if (!entry.flags.is_flagged(pocc_run)) {
        if (entry.getLastCoord() < 0) {
          return true;
        }
      }
    }
    return false;
  }
} // namespace nhull
