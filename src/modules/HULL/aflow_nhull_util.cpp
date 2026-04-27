// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu

#include "aflow_nhull_util.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <AUROSTD/aurostd_xscalar.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_math.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_aflowrc.h"
#include "aflow_nhull_entry.h"
#include "aflow_nhull_facet.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "extern/QHULL/libqhullcpp/QhullFacet.h"
#include "extern/QHULL/libqhullcpp/QhullFacetList.h"
#include "flow/aflow_xclasses.h"
#include "modules/POCC/aflow_pocc.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xvector;

namespace nhull {

  //***************** UTILITY FUNCTIONS

  ///@brief converts string of vectors to single string. used for creating keys in maps for caching functionality for each hull
  /// @authors
  /// mod{NHA,20260206,created}
  string specieVecToString(vector<string> specie_vector) {
    string specie_string;
    for (auto element : specie_vector) {
      specie_string += element;
    }
    return specie_string;
  }

  /// @brief determines interquartile range cutoff value given a list of entries. This is a filter utility.
  /// @param loaded_entries list of entries
  /// @return double
  /// @authors
  /// @mod {NHA,20251005,created}
  double IQRCutoff(vector<Entry>& loaded_entries) {
    // std::map<string, double> loaded_entries
    vector<double> sorted_entries;
    const vector<double> Q1_vector = {0};
    const vector<double> Q3_vector = {0};
    const double IQR_multiple = 3.25;

    //load entries, only keep those below 0 enthalpy
    //Duplicates have been ignored before this point
    for (const auto& entry : loaded_entries) {
      double enthalpy_formation_atom = entry.enthalpy_formation_atom;
      if (enthalpy_formation_atom <= 0) {
        sorted_entries.push_back(entry.enthalpy_formation_atom);
      }
    }
    //convert to xvector
    xvector<double> sorted_entries_xvector = aurostd::vector2xvector<double>(sorted_entries);

    double q1;
    double q2;
    double q3;

    aurostd::getQuartiles(sorted_entries_xvector, q1, q2, q3); // we sort in here
    // cout<<"AUROSTD Quartiles:"<<endl;
    // cout<<"Q1 = "<<q1*1000<<endl;
    // cout<<"Q3 = "<<q3*1000<<endl<<endl;

    // //only use IQR if we have more than 4 entries
    if (sorted_entries.size() > 4) {
      double inner_quartile_range = q3 - q1;
      // cout<<inner_quartile_range*1000<<endl;
      // cout<<(q1 - (IQR_multiple * inner_quartile_range))*1000<<endl;
      return (q1 - (IQR_multiple * inner_quartile_range));
    } else {
      return 0;
    }
  }

  /// @brief takes in compound stoichiometries and returns a reduced stoich map (greatest common factors)
  /// @param compound_map
  /// @return a compound map that has been reduced down by a greatest common factor
  /// @authors
  /// @mod{NHA,20251005,created}
  std::map<string, int> reducedStoich(std::map<string, int> compound_map) {
    vector<int> stoichs;
    int itts = 0;

    for (const auto& specie : compound_map) {
      stoichs.push_back(specie.second);
      itts++;
    }

    // find minimum stoichiometry
    int min = stoichs[0];
    for (int i = 1; i < stoichs.size(); i++) {
      if (min > stoichs[i]) {
        min = stoichs[i];
      };
    }

    // now find greatest common factor and divide each stoich by this
    bool gcf_clear = false;
    int final_min = 1;
    while (!gcf_clear) {
      for (int i = 0; i < stoichs.size(); i++) {
        if (stoichs[i] % min != 0) {
          break;
        };
        if (i == stoichs.size() - 1) {
          gcf_clear = true;
          final_min = min;
        }
      }
      min--;
    }
    // now apply that common factor to all entries in compound:
    for (const auto& specie : compound_map) {
      compound_map.at(specie.first) = specie.second / (final_min);
    }

    return compound_map;
  }

  /// @brief takes integer stoich coords and converts them to relative ratios for qhull input
  /// @param input_coords
  /// @return list of points ready for orgQhull input
  /// @authors
  /// @mod {NHA,20251005,created}
  vector<double> stoichCoord(const vector<double>& input_coords) {
    // takes n-dimensional vector and returns n-1 dimensional vector
    // used in points()
    // first coordinate corresponds to central specie and is not "directly used" in output
    const int dimension = input_coords.size();
    // vector<double> stoich_coords(dimension,0);
    vector<double> stoich_coords(dimension);
    const double central_specie = input_coords[0];
    double net_stoich = 0;

    // add all stoichs together to determine relative stoichs in next step
    for (int i = 0; i < input_coords.size(); i++) {
      net_stoich += input_coords[i];
    }
    // first find ratio of each stoich to the whole:
    for (int i = 0; i < input_coords.size(); i++) {
      if (central_specie != 0 || input_coords[i] != 0) {
        stoich_coords[i] = input_coords[i] / net_stoich;
      } else {
        stoich_coords[i] = 0;
      }
    }

    // assume compounds are loaded by atomic number ordering to ensure correct ordering
    stoich_coords.erase(stoich_coords.begin()); // remove central compound, as this can be found later with (1-sum)
    vector<double> remaining_species = stoich_coords;

    return remaining_species;
  }

  /// @brief takes in a qhull point (no first index), and returns its full relative stoich coords (adds central compound back)
  /// @param qhull_coord qhull point (no first index), either has energy or not for specific cases
  /// @param has_enthalpy true if stoich coord has enthalpy coordinate, false otherwise
  /// @return a single point for orgQhull ingestion
  /// @authors
  /// @mod{NHA,20251005,created}
  vector<double> fullCoords(const vector<double>& qhull_coord, const bool has_enthalpy) {
    // this is a utility function for the squashing function
    // The ordering of species here is the same as in the original compound vector
    // NO ENTHALPY COORD USED HERE
    double total_sum = 0;
    double central_compound = 0;
    vector<double> total_vector;
    size_t max_index = qhull_coord.size();

    if (has_enthalpy) {
      max_index -= 1;
    }

    for (size_t i = 0; i < max_index; i++) {
      total_sum += qhull_coord[i];
    }

    central_compound = 1 - total_sum;
    total_vector.push_back(central_compound);
    for (double element : qhull_coord) {
      total_vector.push_back(element);
    }

    return total_vector;
  }

  /// @brief takes in a stoich point (relative stoich vector, no energy) and returns the integer values of compounds instead of relative fractions
  /// @param stoich_coord relative coordinate with no energy
  /// @return list of integers representing the reduced integer stoichiometry of compound
  /// @authors
  /// @mod  {NHA,20251005,created}
  vector<int> fullIntCoords(const xvector<double>& stoich_coord) {
    // The ordering of species here is the same as in the original compound vector
    double total_sum = 0;
    double central_compound = 0;
    vector<double> total_vector;
    vector<int> itts_array;
    vector<int> true_int_vector;
    double residual;
    int itts = 0;
    int lcm_var;
    int lcm_store;
    int greatest_non_factor = 0;
    bool all_multiple = false;
    int itts2 = 0;

    for (int i = 0; i < stoich_coord.rows; i++) {
      total_sum += stoich_coord[i];
    }

    central_compound = 1 - total_sum;
    total_vector.push_back(central_compound);
    for (int i = 0; i < stoich_coord.rows; i++) {
      total_vector.push_back(stoich_coord[i]);
    }

    // now convert to integers
    // first find least integer that produces an integer when mulitplied by a stoich coord
    for (auto element : total_vector) {
      itts = 0;
      do {
        itts++;
        residual = fabs(round(element * itts) - element * itts);
      } while (residual > .001);

      itts_array.push_back(itts);
    }

    // now find least common multiple of elements in itts_array
    if (itts_array.size() > 1) {
      // first find the greatest lcm_var out of all combinations
      lcm_var = itts_array[0];
      for (int i = 0; i < itts_array.size(); i++) {
        for (int j = 0; j < itts_array.size(); j++) {
          lcm_store = std::lcm(itts_array[i], itts_array[j]);
          if (lcm_store > lcm_var) {
            lcm_var = lcm_store;
          }
        }
      }
      // now test if this lcm_var is a multiple of all entries
      // if not true, multiply by greatest non-factor in list and try again until true

      while (all_multiple == false) {
        // if we went through all array elements without error, we have found lcm_var
        if (all_multiple == false && itts2 == itts_array.size() - 1) {
          break;
        }
        if (lcm_var % itts_array[itts2] != 0 && itts_array[itts2] > greatest_non_factor) {
          greatest_non_factor = itts_array[itts2];
        }
        // if we went through all array elements with error, update lcm_var and loop again:
        if (all_multiple == false && itts2 == itts_array.size() - 1) {
          lcm_var = lcm_var * greatest_non_factor;
          itts2 = -1;
        }
        itts2++;
      }

      // final vector is then each element of total_vector times lcm_var.
      // Round off to nearest int to prevent num error
      for (auto element : total_vector) {
        true_int_vector.push_back((int) round(element * lcm_var));
      }
    } else {
      // if itts_array.size() = 1 this is a pure compound
      true_int_vector = {1};
    }

    return true_int_vector;
  }

  /// @brief Utility method used in decomp ratios: Calculates geometric term used later for calculating relative phase decomp using G-S orthogonalization
  /// @param v vertices of the nearest facet
  /// @param vdif qhull_point - vertex_of_nearest_facet
  /// @return list of the unormalized decomposition ratios
  /// @authors
  /// mod{NHA,20260206,created}
  vector<double> getAlphas(const int dimension, const vector<xvector<double>>& v, const vector<xvector<double>>& vdif) {
    double alpha;
    vector<double> alphas;
    double dist_h;
    xvector<double> k_vec;

    // calculate phase for each vertice:
    for (int i = 0; i < v.size(); i++) {
      // for nth dimension vector populate gs_vectors with n-1 vectors:
      // first vectors are difference in any 2 vectors not related to central vector
      // last vector is difference in central vector to any other vector
      vector<xvector<double>> gs_vectors;
      // pick central index (arbitrary):
      int central_index;
      if (i == 0) {
        central_index = 1;
      } else {
        central_index = 0;
      }
      // create gs vectors here:
      for (int j = 0; j < dimension; j++) {
        if (j != central_index && j != i) {
          const xvector<double> xvector_store = v[central_index] - v[j];
          gs_vectors.push_back(xvector_store);
        }
      }
      // add in last vector that is dependent on current vertice:
      const xvector<double> xvector_store = v[central_index] - v[i];
      gs_vectors.push_back(xvector_store);

      const xvector<double> normal_vec = getHyperplaneNorm(dimension, gs_vectors);

      // check if norm zero, if so this vertice has no contributions to decomp
      if (aurostd::modulus(normal_vec) == 0) {
        alpha = 0;
      } else {
        dist_h = std::abs(scalar_product(normal_vec, vdif[central_index]));
        k_vec = vdif[i] / aurostd::modulus(vdif[i]);
        alpha = dist_h / std::abs(scalar_product(k_vec, normal_vec));
      }
      alphas.push_back(alpha);
    }
    return alphas;
  }

  /// @brief Takes a qhull point and vertices of facet and returns the decomposition percents for each component
  /// @param dimension dimension of the hull (number of stoichiometric coordinates + 1 energy coord)
  /// @param vertices vertices of the nearest facet
  /// @param qhull_point point of interest to find decomposition of
  /// @return list of decompostion ammounts
  /// @authors
  /// @mod {NHA,20251005,created}
  vector<double> decompRatios(const int dimension, const vector<xvector<double>>& vertices, xvector<double> qhull_point) {
    // The ordering of return vector is identical to sortAtomic(); return vector same length as dimension
    bool zero_distance_flag = false; // bool used to flag if compound is directly above a hull point
    size_t zero_distance_index;
    vector<double> return_ratios; // return variable

    // remove energy coordinate
    auto vector_qhull_point = aurostd::xvector2vector(qhull_point);
    vector_qhull_point.pop_back();
    auto xvector_qhull_point = aurostd::vector2xvector(vector_qhull_point, 0);

    // now find distances from vertices to xvector_qhull_point

    vector<double> vdist;
    for (size_t i = 0; i < vertices.size(); i++) {
      const double distance = aurostd::modulus(xvector_qhull_point - vertices[i]);

      if (distance == 0) {
        // check if point is directly above hull point
        zero_distance_flag = true;
        zero_distance_index = i;
        break;
      }
      vdist.push_back(distance);
    }

    if (!zero_distance_flag) {
      // find the vectors from vertices to xvector_qhull_point
      vector<xvector<double>> vdif;
      for (auto vertex : vertices) {
        const xvector<double> xvector_store = xvector_qhull_point - vertex;
        vdif.push_back(xvector_store);
      }

      // we use gsNorm for dimension > 2 and a special case for dimension = 2 due to its unique simplicity
      switch (dimension) {
        case 2: {
          const double phase_0 = vdist[1] / (vdist[0] + vdist[1]);
          const double phase_1 = vdist[0] / (vdist[0] + vdist[1]);
          return_ratios = {phase_0, phase_1};
          break;
        }
        default: {
          vector<double> alphas = getAlphas(dimension, vertices, vdif);

          for (int i = 0; i < dimension; i++) {
            return_ratios.push_back(alphas[i] / (alphas[i] + vdist[i]));
          }
          break;
        }
      }
      // clean return ratios from very small non-zero values:
      return_ratios = aurostd::absorbResidualsToMax(return_ratios);
    } else {
      for (int i = 0; i < dimension; i++) {
        if (i == zero_distance_index) {
          return_ratios.push_back(1);
        } else {
          return_ratios.push_back(0);
        }
      }
    }

    return return_ratios;
  }

  ///@brief returns the surface normal coordinates of a qhull facet
  ///@param current_facet qhull facet
  ///@return surface normal of a qhull facet as a xvector<double>
  /// @authors
  /// mod{NHA,20260206,created}
  xvector<double> getQhullFacetNormal(const orgQhull::QhullFacet& current_facet) {
    const double* surface_normal = current_facet.hyperplane().coordinates();
    const int dimension = current_facet.dimension();
    vector<double> surface_normal_coordinates;
    for (int i = 0; i < dimension; i++) {
      surface_normal_coordinates.push_back(-*surface_normal);
      // include negative sign here since qhull assumes surface_normals face outwards of hull
      surface_normal++;
    }

    return aurostd::vector2xvector(surface_normal_coordinates);
  }

  /// @brief returns the hyperplane offset from a qhull facet
  /// @authors
  /// mod{NHA,20260206,created}
  double getQhullFacetOffset(const orgQhull::QhullFacet& current_facet) {
    return -current_facet.hyperplane().offset();
  }

  ///@brief Takes in surface_normal vector of facet and returns its energy (last) coordinate
  /// @authors
  /// mod{NHA,20260206,created}
  double getNormEnergy(const aurostd::xvector<double>& surface_normal) {
    if (fabs(aurostd::modulus(surface_normal) - 1) > NHULL_NORMAL_VECTOR_TOL) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "surface normal is not normalized!");
    }
    return surface_normal[surface_normal.urows]; // last entry of surface_normal_coord
  }

  ///@brief calculates vertical distance to hull in energy coordinates
  /// @authors
  /// mod{NHA,20260206,created}
  double disttoHull(const double norm_energy, const double distance) {
    return fabs(distance / norm_energy); // since norm_energy and distance can be negative (and we onl  y care about distance to hull), take fabs() here
  }

  //************* OUTPUT UTILITY METHODS

  ///@brief takes a list of coordinates and displays them to console
  /// @authors
  /// mod{NHA,20260206,created}
  void printPoints(vector<aurostd::xvector<double>> qhull_many_points) {
    for (int i = 0; i < qhull_many_points.size(); i++) {
      for (int j = qhull_many_points[i].lrows; j <= qhull_many_points[i].urows; j++) {
        cout << qhull_many_points[i][j] << " ";
      }
      cout << endl;
    }
  }

  //************* Utility methods associated with qhull runs:

  ///@brief convert lower_dim_point_list to vector<xvector<double>> for easy concatenation
  vector<vector<double>> compoundXvectorConversion(const vector<xvector<double>>& v) {
    vector<vector<double>> vector_lower_dim_point_list;
    for (auto point : v) {
      vector_lower_dim_point_list.push_back(aurostd::xvector2vector(point));
    }
    return vector_lower_dim_point_list;
  }

  ///@brief turns a 2d vector into a 1d vector. Used for qhull ingestion.
  vector<double> flatten_vector(const vector<vector<double>>& v) {
    vector<double> flattened_vector;
    // now flatten entries of new_qhull_vectors for qhull ingestion
    for (const auto& element : v) {
      for (const auto subelement : element) {
        flattened_vector.push_back(subelement);
      }
    }
    return flattened_vector;
  }

  ///@brief turns a 2d vector into a 1d vector. Used for qhull ingestion.
  vector<double> flatten_vector(const vector<aurostd::xvector<double>>& v) {
    vector<double> flattened_vector;
    // now flatten entries of new_qhull_vectors for qhull ingestion
    for (const xvector<double>& element : v) {
      for (int i = element.lrows; i <= element.urows; i++) {
        flattened_vector.push_back(element[i]);
      }
    }

    return flattened_vector;
  }

  ///@brief finds nearest facet to an entry and calculates distance to hull and phase decomp
  ///@param entries_on_hull
  ///@param qhull_point point to determine distance from hull
  ///@param facet_list list of all hull facets
  ///@return tuple of the distance to hull, and the decomposition map in terms of entries and their ratios
  /// @authors
  /// mod{NHA,20260206,created}
  std::tuple<double, std::map<Entry, double>> nearestDistDecomp(const vector<Entry>& entries_on_hull, aurostd::xvector<double> qhull_point, const vector<NhullFacet>& facet_list) {
    NhullFacet closest_facet;

    // shortest travel distance to facet of hull (not vertical)
    // Qhull will use negative distances so to prevent ambiguity absolute value used here
    double shortest_distance = -1; // set shortest distance to -1; should be replaced by finite value in main loop
    vector<double> test;
    for (const auto& facet : facet_list) {
      const double norm_energy = getNormEnergy(facet.m_normal);

      // If norm_energy <= 0, the surface normal is upward or perfectly horizontal, and the point will never intersect so ignore
      // also ignore zero-energy tie line as this isn't part of hull
      if (norm_energy > 0) {
        const double distance = facet.distance(qhull_point);
        //if (distance > 0) //special case that occurs for points below the hull
        if (checkPointBelowHull(facet, qhull_point)) {
          closest_facet = facet;
          shortest_distance = disttoHull(norm_energy, fabs(distance));
          break;
        }
      }
    }

    if (shortest_distance < 0) {
      stringstream message;
      message << "No facet found for point!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }

    // now find decomposition compounds and relative amounts:
    std::map<Entry, double> decomp_map;
    try {
      decomp_map = decompReactionEntries(entries_on_hull, closest_facet, qhull_point);
    } catch (aurostd::xerror& e) {
      //don't crash run, but make a note of it
      std::cerr << "Error in function nearestDistDecomp: failure to produce decomposition map! Assigned a trivial map." << std::endl;
      Entry null_entry;
      std::map<Entry, double> trivial_map = {
          {null_entry, 1}
      };
      decomp_map = trivial_map;
    }

    return {shortest_distance, decomp_map};
  }

  /// @brief For n-1 and stability criterion distance to hull needs to be calc for points below energy hull
  /// @note Uses barycentric coordinates to determine whether point lies within simplex: https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Determining_location_with_respect_to_a_triangle
  /// @returns true if point is vertically below/above a given facet
  /// @authors
  /// mod{NHA,20260206,created}
  bool checkPointBelowHull(const NhullFacet& facet, const aurostd::xvector<double>& qhull_point) {
    vector<double> qhull_point_vector = aurostd::xvector2vector(qhull_point);
    //all points have energy coord removed
    if (facet.m_dimension != qhull_point.urows) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Facet and point do not share dimension!");
    }

    vector<xvector<double>> vertices = NhullFacet::nhullFacetToVertices(facet, true);
    vector<vector<double>> vertices_test = nhull::compoundXvectorConversion(vertices);
    qhull_point_vector.pop_back();
    xvector<double> point = aurostd::vector2xvector(qhull_point_vector);
    size_t dimension = vertices.size();
    aurostd::xmatrix<double> barycentric_matrix(dimension - 1, dimension - 1, 1, 1);
    aurostd::xvector<double> rdif = point - vertices.back();

    //set barycentric_matrix:
    for (size_t i = 0; i < dimension - 1; i++) {
      barycentric_matrix.setcol(vertices[i] - vertices.back(), i + 1);
    }

    aurostd::xmatrix<double> barycentric_matrix_inverse = inverse(barycentric_matrix);
    aurostd::xvector<double> barycentric_coordinates = barycentric_matrix_inverse * rdif;

    if (aurostd::min(barycentric_coordinates) + NHULL_HULL_CHECK_TOL < 0) {
      return false;
    }

    double last_barycentric_coordinate = 1 - aurostd::sum(barycentric_coordinates);
    if (last_barycentric_coordinate + NHULL_HULL_CHECK_TOL < 0) {
      return false;
    }

    return true;
  }

  ///@brief takes a compound string ex: "Mn3Pd6" and returns a reduced stoichiometric string ex: "Mn1Pd2"
  /// @authors
  /// mod{NHA,20260206,created}
  string reducedCompoundString(const string& compound) {
    vector<string> ele;
    vector<double> comp;
    vector<int> int_comp;
    bool integer_coords = true;
    aurostd::elementsFromCompositionString(compound, ele, comp);

    //ensure in reduced form if integers:
    for (const auto value : comp) {
      if (std::floor(value) != std::ceil(value)) {
        integer_coords = false;
        break;
      }
    }

    if (integer_coords) {
      //now reduce common factors:
      std::transform(comp.begin(), comp.end(), std::back_inserter(int_comp), [](double d) { return static_cast<int>(d); });
      int gcd = 1;
      aurostd::GCD(int_comp, gcd);
      std::transform(int_comp.begin(), int_comp.end(), int_comp.begin(), [gcd](int d) { return d / gcd; });

      //create return string here:
      string formatted_compound;
      for (size_t i = 0; i < int_comp.size(); i++) {
        if (int_comp[i] != 0) {
          formatted_compound += ele[i] + std::to_string(int_comp[i]);
        }
      }
      return formatted_compound;
    }
    return compound;
  }

  /// @brief Returns list of auid, degeneracy value pairs from a directory of given pocc compound. User provides function with directory.
  /// @param root_path User-provided directory
  /// @return vector of auid, degeneracy value pairs
  /// @authors
  /// @mod {NHA,20251013,created}
  vector<std::pair<string, unsigned long long int>> loadPocc(const string& root_path) {
    //root path example: /common/LIB6/RAW/CHf_pvNb_svTa_pvTi_svZr_sv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF/
    vector<std::pair<string, unsigned long long int>> dg_list;

    for (const auto& dir_pocc_entry : fs::directory_iterator(root_path)) {
      if (dir_pocc_entry.is_directory()) {
        string aflowin_path = (string) dir_pocc_entry.path() + "/aflowlib.json";
        string poscar_path = (string) dir_pocc_entry.path() + "/POSCAR.orig";

        //process aflowlib.json for auid here:
        aurostd::JSON::object aflowjson = aurostd::JSON::loadFile(aflowin_path);
        string auid = (string) aflowjson["auid"];

        //process POSCAR here:
        //need only first line from POSCAR to get DG
        std::ifstream poscar_file(poscar_path);
        string first_line;
        getline(poscar_file, first_line);

        unsigned long long int dg = pocc::getDGFromXStructureTitle(first_line);

        std::pair<string, unsigned long long int> curr_pair = {auid, dg};
        dg_list.push_back(curr_pair);
      }
    }
    return dg_list;
  }

  /// @brief applies an interquartile filter to points corresponding to a binary calculation
  /// @param m_half_hull true if there are candidate hull points that lie below zero-energy tie line
  /// @authors
  /// mod{NHA,20260206,created}
  void IQRFilter(vector<Entry>& entries, bool pocc_run, bool m_half_hull, const int dimension) {
    //Take only points that are not flagged, then organize into binaries and perform outlier detection binary-by-binary
    vector<size_t> allowed_entry_indices;
    for (int i = 0; i < entries.size(); i++) {
      if (!entries[i].flags.is_flagged(pocc_run)) {
        allowed_entry_indices.push_back(i);
      }
    }

    xvector<int> elements_present;
    vector<uint> tmp_store;
    vector<uint> outliers_found;
    //get all binary combos:
    aurostd::xcombos xc(dimension, 2, 'C');
    while (xc.increment()) {
      elements_present = aurostd::vector2xvector<int>(xc.getCombo());
      tmp_store = getOutliers(entries, allowed_entry_indices, elements_present, m_half_hull);
      outliers_found.insert(outliers_found.end(), tmp_store.begin(), tmp_store.end());
    }

    //now flag all outliers:
    for (const uint found_index : outliers_found) {
      entries[found_index].flags.outside_IQR = true;
    }
  }

  ///@brief utility function for IQRfilter(): calculates a threshold value used to exclude entries
  ///@param energies list of energies from a filtered entries list
  ///@param upper_threshold largest allowed energy (passed by reference)
  ///@param lower_threshold smallest allowed energy (passed by reference)
  /// @authors
  /// mod{NHA,20260206,created}
  void calculateOutlierThreshold(const xvector<double>& energies, double& upper_threshold, double& lower_threshold) {
    stringstream message;

    upper_threshold = AUROSTD_MAX_DOUBLE; // effectively NOT a threshold
    lower_threshold = -AUROSTD_MAX_DOUBLE; // effectively NOT a threshold

    const uint iqr_count_threshold = 4; // 3 results in degenerate quartile indices
    if ((uint) energies.rows < iqr_count_threshold) { // ALWAYS not enough points to do statistics (need 3 quartiles)
      return;
    }

    double q1;
    double q2;
    double q3;
    // so we sort full anyway to be completely robust, should be easy considering how we sorted before
    aurostd::getQuartiles(energies, q1, q2, q3); // we sort in here

    double lower_anchor = q1;
    double upper_anchor = q3;
    double range = q3 - q1;
    double multiplier = DEFAULT_NHULL_OUTLIER_MULTIPLIER;

    upper_threshold = upper_anchor + (multiplier * range);
    lower_threshold = lower_anchor - (multiplier * range);
  }

  /// @brief returns indices of entries that fall outside of the IQR thresholds
  /// @param half_hull true if there are candidate hull points that lie below zero-energy tie line
  /// @param allowed_entry_indices list of unflagged entries
  /// @return list of entry indexes of flagged outliers corresponding to the global list: m_points[]
  /// @authors
  /// mod{NHA,20260206,created}
  vector<uint> getOutliers(vector<Entry>& entries, const vector<size_t>& allowed_entry_indices, const xvector<int>& elements_present, bool half_hull) {
    stringstream message;
    vector<uint> points_to_consider;

    if (half_hull) {
      // we only care about points above/below hull
      //keep points if they share elements present
      //keep points if they are within half hull (inside hull)
      for (const size_t index : allowed_entry_indices) {
        if (entries[index].getElementsPresent() != elements_present) {
          continue;
        }

        if (entries[index].isWithinHalfHull(half_hull)) {
          points_to_consider.push_back(index);
        }
      }
    } else {
      //keep points if they share elements present
      for (const size_t index : allowed_entry_indices) {
        if (entries[index].getElementsPresent() != elements_present) {
          continue;
        }
        points_to_consider.push_back(index);
      }
    }

    if (half_hull && aurostd::sum(elements_present) == 2) {
      const uint binaries_half_hull_threshold = DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES;
      if (points_to_consider.size() < binaries_half_hull_threshold) {
        vector<uint> outliers;
        return outliers; //if under threshold, we return no outliers for this binary
      }
    }

    return calculateOutliers(entries, points_to_consider, half_hull);
  }

  ///@brief determines upper and lower IQR filter bounds and returns entries that fall outside
  ///@param half_hull true if there are candidate hull points that lie below zero-energy tie line
  ///@return list of entry indexes of flagged outliers corresponding to the global list: m_points[]
  /// @authors
  /// mod{NHA,20260206,created}
  vector<uint> calculateOutliers(vector<Entry>& entries, const vector<uint>& points_to_consider, bool half_hull) {
    // it is very important that we do not define outliers using std/mean, as these are
    // very sensitive to outliers
    // instead, use median
    // see discussion here:  doi=10.1016/j.jesp.2013.03.013
    vector<uint> outliers;

    // get vector of last coords (we want to find outliers in this dimension)
    vector<double> energies; // vector and not xvector because we need push_back()
    uint i_point = AUROSTD_MAX_UINT;
    for (size_t i = 0, fl_size_i = points_to_consider.size(); i < fl_size_i; i++) {
      i_point = points_to_consider[i];
      energies.push_back(entries[i_point].getLastCoord());
    }

    double upper_threshold;
    double lower_threshold;
    const xvector<double> energies_converted = aurostd::vector2xvector<double>(energies);
    calculateOutlierThreshold(energies_converted, upper_threshold, lower_threshold);

    if (half_hull) {
      // look at lower range
      for (size_t i = 0, fl_size_i = points_to_consider.size(); i < fl_size_i; i++) {
        i_point = points_to_consider[i];
        if (entries[i_point].getLastCoord() < lower_threshold) {
          outliers.push_back(i_point);
        }
      }
    } else {
      // look at upper range
      for (size_t i = 0, fl_size_i = points_to_consider.size(); i < fl_size_i; i++) {
        i_point = points_to_consider[i];
        if (entries[i_point].getLastCoord() > upper_threshold) {
          outliers.push_back(i_point);
        }
      }
    }
    return outliers;
  }

  //********************** aflowlib sorting functions:

  /// @brief constructor for a LIB2 sorting object: copies certain member vars from aflowlib_entry that are used to order aflowlib entries
  /// @param entry an aflowlib entry
  /// @param index object gets index according to where it was instantiated in for loop. Index used later for ease of manipulation.
  /// @param specie_vector global list of elements in this hull
  /// @authors
  /// @mod {NHA,20260210,created}
  _aflowlib_entry_LIB2sorting::_aflowlib_entry_LIB2sorting(aflowlib::_aflowlib_entry& entry, uint index, vector<string>& specie_vector, std::ofstream& FileMESSAGE, std::ostream& oss) :
      xStream(FileMESSAGE, oss), m_index(index) {
    m_auid = entry.auid;
    m_aurl = entry.aurl;
    m_catalog = entry.catalog;
    m_prototype = entry.prototype;
    m_nspecies = entry.nspecies;
    m_velements_nhull = specie_vector;
    m_species_AURL = entry.getSpeciesAURL(*p_FileMESSAGE, *p_oss);
  }

  ///@brief operator used for sorting aflowlib_entry objects
  ///@return true if other is "greater than" this entry
  /// @authors
  /// mod{NHA,20260206,created}
  bool _aflowlib_entry_LIB2sorting::operator<(const _aflowlib_entry_LIB2sorting& other) const {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (m_catalog != other.m_catalog) {
      return m_catalog < other.m_catalog;
    }
    if (m_prototype != other.m_prototype) {
      return m_prototype < other.m_prototype;
    }
    if (m_catalog == "LIB2" && m_nspecies == 1 && other.m_nspecies == 1) {
      if (m_species_AURL.size() != other.m_species_AURL.size()) {
        return m_species_AURL.size() < other.m_species_AURL.size();
      } // effectively sorting by catalog
      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " m_auid=" << m_auid << " m_aurl=" << m_aurl << " m_species_AURL=" << aurostd::joinWDelimiter(m_species_AURL, ",") << endl;
        std::cerr << __AFLOW_FUNC__ << " other.m_auid=" << other.m_auid << " other.m_aurl=" << other.m_aurl << " other.m_species_AURL=" << aurostd::joinWDelimiter(other.m_species_AURL, ",") << endl;
      }
      bool i_within_hull = true;
      for (size_t k = 0; k < m_species_AURL.size() && i_within_hull; k++) {
        if (!aurostd::WithinList(m_velements_nhull, m_species_AURL[k])) {
          i_within_hull = false;
        }
      }
      bool j_within_hull = true;
      for (size_t k = 0; k < other.m_species_AURL.size() && j_within_hull; k++) {
        if (!aurostd::WithinList(m_velements_nhull, other.m_species_AURL[k])) {
          j_within_hull = false;
        }
      }
      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " i_within_hull=" << i_within_hull << " j_within_hull=" << j_within_hull << endl;
      }
      if (i_within_hull != j_within_hull) {
        return i_within_hull;
      }
      for (size_t k = 0; k < m_species_AURL.size(); k++) {
        if (m_species_AURL[k] != other.m_species_AURL[k]) {
          return m_species_AURL[k] < other.m_species_AURL[k];
        }
      }
    }
    //if we exhausted all sort methods, sort directories lexicographically as they are guaranteed to be different.
    return m_aurl < other.m_aurl;
  }

  /// @brief utility function for initializing the formation enthalpies of nhull::Entry objects
  /// @note uses CCE if available
  /// @return formation enthalpy per atom
  /// @authors
  /// @mod {NHA,20251005,created}
  double H_f_atom(const aflowlib::_aflowlib_entry& entry, bool neglect_cce, bool& cce_used) {
    if (neglect_cce == false) {
      const double d_tmp = entry.enthalpyFormationAtom(cce_used, 0); // entry.enthalpy_formation_atom - ADDING CCE @ 0K
      return d_tmp;
    }
    return entry.enthalpy_formation_atom;
  }

  ///@brief determines # of elements present in a qhull point: first specie removed, energy coordinate included
  ///@return number of atomic species present
  /// @authors
  /// @mod {NHA,20251005,created}
  int dimensionOfQhullPoint(aurostd::xvector<double> qhull_point) {
    double coord_sum = 0;
    int number_of_empty_coords = 0;
    //ignore energy coordinate
    for (int i = qhull_point.lrows; i <= qhull_point.urows - 1; i++) {
      if (aurostd::zeroWithinTol(qhull_point[i], ZERO_TOL)) {
        number_of_empty_coords++;
      } else {
        coord_sum += qhull_point[i];
      }
    }

    if (aurostd::zeroWithinTol(coord_sum - 1, ZERO_TOL)) {
      return qhull_point.rows - 1 - number_of_empty_coords;
    }

    return qhull_point.rows - number_of_empty_coords;
  }

} // namespace nhull
