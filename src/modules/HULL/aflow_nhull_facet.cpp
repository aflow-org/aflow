#include "aflow_nhull_facet.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_nhull_util.h"
#include "extern/QHULL/libqhull_r/libqhull_r.h"
#include "extern/QHULL/libqhullcpp/Qhull.h"
#include "extern/QHULL/libqhullcpp/QhullFacet.h"
#include "extern/QHULL/libqhullcpp/QhullPoint.h"
#include "extern/QHULL/libqhullcpp/QhullVertex.h"
#include "extern/QHULL/libqhullcpp/QhullVertexSet.h"
#include "modules/HULL/aflow_nhull_entry.h"

using aurostd::xvector;
using std::string;
using std::vector;

namespace nhull {

  /// @brief utility function assigns distance to hull values after running qhull and determining facet points
  /// @param hull_entry entry to modify
  /// @authors
  /// mod{NHA,20260206,created}
  void NhullFacet::initializeHullEntry(Entry& hull_entry) {
    if (hull_entry.is_hull_point) {
      hull_entry.distance_hull_enthalpy_formation_atom = 0;
      std::map<nhull::Entry, double> trivial_map = {
          {hull_entry, 1}
      };
      hull_entry.m_nhull_phase_decomp_entries = trivial_map;
      hull_entry.m_initialized = true;
    } else {
      std::stringstream message;
      message << "Tried to initialize a non-hull point" << std::endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }
  }

  /// @brief Unit conversion of facet vertices from eV to meV
  /// @authors
  /// mod{NHA,20260206,created}
  void NhullFacet::enthalpyUnitConversion() {
    for (auto& entry : entries) {
      entry.enthalpyUnitConversion();
    }
    for (auto& facetpoint : vertices) {
      facetpoint[facetpoint.urows] *= 1000;
    }
  }

  /// @brief takes output from QHULL and stores simple facet data in NhullFacet
  /// @param nhull_facets
  /// @authors
  /// mod{NHA,20260206,created}
  void NhullFacet::qhullRunToNhullFacets(std::vector<NhullFacet>& nhull_facets, orgQhull::Qhull& qhull, const vector<Entry>& entries) {
    nhull_facets.clear();
    const orgQhull::QhullFacet first_facet = qhull.beginFacet();
    const int dimension = first_facet.dimension();
    const orgQhull::QhullFacetList& facet_list = qhull.facetList();

    for (const auto& facet : facet_list) {
      qhullFacetToNhullFacet(nhull_facets, entries, facet, dimension);
    }
  }

  ///@brief returns list of vertices of a nhull facet in vector form
  ///@param remove_energy if true, removes the last (energy) coordinate
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> NhullFacet::nhullFacetToVertices(const NhullFacet& facet, const bool remove_energy) {
    std::vector<aurostd::xvector<double>> vertices = facet.vertices;
    std::vector<aurostd::xvector<double>> ret_verts;
    if (remove_energy) {
      for (auto& vertice : vertices) {
        auto tmp = aurostd::xvector2vector(vertice);
        tmp.pop_back();
        ret_verts.push_back(aurostd::vector2xvector(tmp));
      }
    } else {
      ret_verts = vertices;
    }
    return ret_verts;
  }

  ///@brief returns list of vertices of a qhull facet in vector form
  ///@param remove_energy if true, removes the last (energy) coordinate
  /// @authors
  /// mod{NHA,20260206,created}
  vector<xvector<double>> NhullFacet::qhullFacetToVertices(const orgQhull::QhullFacet& facet, const bool remove_energy) {
    int facet_dim = facet.dimension();
    if (remove_energy) {
      facet_dim = facet_dim - 1;
    }
    vector<double> vertice;
    vector<xvector<double>> vertices;
    vector<Entry> found_entries;

    // first ensure that facets and points share dimension
    // first get species from facet
    const orgQhull::QhullVertexSet curr_vertexset = facet.vertices();
    const vector<orgQhull::QhullVertex> curr_vertexvector = curr_vertexset.toStdVector();

    // retrieve vertices
    for (const auto& curr_vertex : curr_vertexvector) {
      orgQhull::QhullPoint point = curr_vertex.point();
      coordT* point_coordinates = point.coordinates();
      vertice.clear();
      for (int i = 0; i < facet_dim; i++) {
        vertice.push_back(*point_coordinates);
        point_coordinates++;
      }
      // make sure point is below zero-energy tie-line:
      vertices.push_back(aurostd::vector2xvector(vertice));
    }
    return vertices;
  }

  /// @brief Takes a facet of type orgQhull and converts to an nhull native facet for ease of calculation
  /// @param nhull_facets vector of uninitialized nhullFacets
  /// @param filtered_entries vector of entries that are on the hull
  /// @param facet an orgQhull facet
  /// @param dimension dimension of facets
  /// @param specie_vector vector of atomic species that make up points in this space
  /// @author
  /// mod{NHA,20251005,created}
  void NhullFacet::qhullFacetToNhullFacet(vector<NhullFacet>& nhull_facets, vector<Entry>& entries, const orgQhull::QhullFacet& facet, const int dimension, const vector<string>& specie_vector) {
    const int facet_dim = facet.dimension();
    vector<aurostd::xvector<double>> vertices;
    vector<Entry> found_entries;

    // first ensure that facets and points share dimension
    if (facet_dim == dimension) {
      // retrieve vertices
      vertices = qhullFacetToVertices(facet, false);
    } else {
      const std::string message = "Facet and qhull coordinate do not share dimensions!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    // now find associated hull entries:
    for (Entry& entry : entries) {
      if (entry.lowest_energy_alloy) {
        auto entry_point = entry.compoundToPoint(specie_vector, false);
        for (auto vertice : vertices) {
          if (aurostd::identical(entry_point, vertice, QHULL_POINT_TOL)) {
            //update entry flags here:
            entry.is_hull_point = true;
            initializeHullEntry(entry);
            found_entries.push_back(entry);
          }
        }
      }
    }

    NhullFacet nhull_facet(specie_vector, vertices, found_entries, getQhullFacetNormal(facet), getQhullFacetOffset(facet));
    nhull_facets.push_back(nhull_facet);
  }

  /// @brief overload function: lightweight facet initializer used for simple hull runs (no thermodynamic analysis)
  /// @authors
  /// mod{NHA,20260206,created}
  void NhullFacet::qhullFacetToNhullFacet(vector<NhullFacet>& nhull_facets, const std::vector<Entry>& hull_points, const orgQhull::QhullFacet& facet, const int dimension) {
    const int facet_dim = facet.dimension();
    vector<aurostd::xvector<double>> vertices;
    vector<Entry> found_entries;

    // first ensure that facets and points share dimension
    if (facet_dim == dimension) {
      // retrieve vertices
      vertices = qhullFacetToVertices(facet, false);
    } else {
      const std::string message = "Facet and qhull coordinate do not share dimensions!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    // now find associated hull entries:
    for (const Entry& entry : hull_points) {
      for (const auto& vertice : vertices) {
        if (entry.entryEqualsPoint(vertice)) {
          found_entries.push_back(entry);
        }
      }
    }

    const NhullFacet nhull_facet(found_entries, getQhullFacetNormal(facet));
    nhull_facets.push_back(nhull_facet);
  }

  /// @brief performs a hull calculation given lists of lower dimensional points and entries
  /// @return vector of NhullFacets
  /// author
  /// mod{NHA,20251005,created}
  vector<NhullFacet> NhullFacet::getFacetData(const vector<string>& specie_vector, const vector<xvector<double>>& lower_dim_point_list, vector<Entry>& entries, const int dimension) {
    const char* qhull_command2 = "Qt"; // triangulated output: necessary for correct phase decomposition

    // convert lower_dim_point_list to vector<xvector<double>> for easy concatenation
    vector<vector<double>> vector_lower_dim_point_list = compoundXvectorConversion(lower_dim_point_list);
    vector<Entry> filtered_entries = Entry::getLowestEnergyEntries(entries, dimension); //filter out lower dimension points since we have already calculated those

    // done at max dimension here:
    vector<vector<double>> point_list = compoundXvectorConversion(entriesToPoints(filtered_entries, specie_vector));
    point_list.insert(point_list.end(), vector_lower_dim_point_list.begin(), vector_lower_dim_point_list.end());
    vector<double> flattened_vector;

    if (point_list.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No points found to ingest into qhull!");
    }

    for (const auto& element : point_list) {
      // flatten vector
      for (const auto subelement : element) {
        flattened_vector.push_back(subelement);
      }
    }
    const int qhull_dim = dimension;
    vector<double> qhull_point_list = flattened_vector;
    const int point_count = qhull_point_list.size() / qhull_dim;
    const double* qhull_input_points = &qhull_point_list[0];
    const orgQhull::Qhull qhull("test", qhull_dim, point_count, qhull_input_points, qhull_command2);

    const orgQhull::QhullFacet first_facet = qhull.beginFacet();
    const int facet_dimension = first_facet.dimension();
    orgQhull::QhullFacet curr_facet = first_facet;
    vector<NhullFacet> nhull_facets;
    bool break_condition = false;

    do {
      // check to make sure current facet does not have vertical or perfectly horizontal normal
      if (getNormEnergy(getQhullFacetNormal(curr_facet)) > QHULL_FACET_NORM_CUTOFF) {
        qhullFacetToNhullFacet(nhull_facets, entries, curr_facet, facet_dimension, specie_vector);
      }

      if (curr_facet.hasNext()) {
        curr_facet = curr_facet.next();
      } else {
        break_condition = true;
      }
    } while (break_condition == false);

    return nhull_facets;
  }

  //

  /// @brief takes a facet and a qhull coordinate and determines its decompositions; returns a vector of those decomps
  /// @param alloy_name_store stores strings of each element in compound
  /// @param closest_facet the nearest facet to qhull_point
  /// @param qhull_point the point of interest in stoichiometry-energy space
  /// @return a compound map: key is compound string and value is the percent decomposition
  /// author
  /// mod{NHA,20251005,created}
  std::map<string, double> decompReaction(const vector<string>& alloy_name_store, const orgQhull::QhullFacet& closest_facet, const aurostd::xvector<double>& qhull_point) {
    const int closest_facet_dimension = closest_facet.dimension();
    vector<double> vertice;
    vector<xvector<double>> vertices;
    orgQhull::QhullPoint point;
    coordT* point_coordinates;
    vector<double> decomp_ratios;
    std::map<string, double> decomp_map;

    // first ensure that facets and points share dimension
    if (closest_facet_dimension == qhull_point.rows) {
      // first get species from facet
      const orgQhull::QhullVertexSet curr_vertexset = closest_facet.vertices();
      const vector<orgQhull::QhullVertex> curr_vertexvector = curr_vertexset.toStdVector();

      for (const auto& curr_vertex : curr_vertexvector) {
        point = curr_vertex.point();
        point_coordinates = point.coordinates();
        vertice.clear();
        for (int i = 0; i < closest_facet_dimension - 1; i++) {
          // LOOP IGNORES ENERGY COORDINATE BY DESIGN THROUGH UPPER BOUND closest_facet_dimension-1
          vertice.push_back(*point_coordinates);
          point_coordinates++;
        }
        const xvector<double> xvector_coord = aurostd::vector2xvector(vertice, 0);
        vertices.push_back(xvector_coord);
      }

      // determine relative amounts of each compound:
      decomp_ratios = decompRatios(closest_facet_dimension, vertices, qhull_point);

      // Now convert to map:
      // first find each compound by taking vertices and converting to string
      // alloy_name_store has already been sorted up to this point
      const std::map<string, double> vertice_compounds;
      for (size_t i = 0; i < vertices.size(); i++) {
        vector<int> int_coords = fullIntCoords(vertices[i]);
        string compound_string;
        for (size_t j = 0; j < alloy_name_store.size(); j++) {
          if (int_coords[j] != 0) {
            compound_string += alloy_name_store[j] + std::to_string(int_coords[j]);
          }
        }

        decomp_map[compound_string] = decomp_ratios[i];
      }
    } else {
      const std::string message = "Closest_facet_dimension and qhull_point do not share size!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    return decomp_map;
  }

  /// @brief determines the entries (phases) and ratios of the nearest facet (nearest vertical distance)
  /// @param closest_facet  the nearest facet to given qhull_point from which to determine phase decomp reaction
  /// @return map of the entries associated with the phase decomp and their ratios
  /// @authors
  /// mod{NHA,20260206,created}
  std::map<Entry, double> decompReactionEntries(const vector<Entry>& entries_on_hull, const NhullFacet& closest_facet, const aurostd::xvector<double>& qhull_point) {
    const int closest_facet_dimension = closest_facet.m_dimension;
    vector<double> vertice;
    vector<xvector<double>> vertices_no_energy;
    vector<double> decomp_ratios;
    std::map<Entry, double> decomp_map;

    // first ensure that facets and points share dimension
    if (closest_facet_dimension == qhull_point.rows) {
      for (auto vertice : closest_facet.vertices) {
        // LOOP IGNORES ENERGY COORDINATE BY DESIGN
        auto tmp = aurostd::xvector2vector(vertice);
        tmp.pop_back();
        const xvector<double> xvector_coord = aurostd::vector2xvector(tmp, 0);
        vertices_no_energy.push_back(xvector_coord);
      }

      // determine relative amounts of each compound:
      decomp_ratios = decompRatios(closest_facet_dimension, vertices_no_energy, qhull_point);

      // Now convert to map:
      const std::map<string, double> vertice_compounds;
      for (size_t i = 0; i < closest_facet.vertices.size(); i++) {
        Entry found_entry;
        bool found_point = Entry::getEntryFromPoint(entries_on_hull, found_entry, closest_facet.vertices[i]); //fullIntCoords(vertices[i]);
        if (found_point) {
          decomp_map[found_entry] = decomp_ratios[i];
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No vertice found");
        }
      }
    } else {
      const std::string message = "Closest_facet_dimension and qhull_point do not share size!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    return decomp_map;
  }

  /// @brief default constructor
  NhullFacet::NhullFacet() {}

  /// @brief lightweight constructor for general non-thermodynamic hulls
  NhullFacet::NhullFacet(const std::vector<aurostd::xvector<double>>& vertices, const aurostd::xvector<double> normal, const double offset, const int dimension) :
      vertices(vertices), m_normal(normal), hyperplane_offset(offset), m_dimension(dimension) {}

  /// @brief constructor used initialize from a qhull facet during thermodynamic calcs
  NhullFacet::NhullFacet(const std::vector<std::string>& specie_vector, const std::vector<aurostd::xvector<double>>& vertices, const std::vector<Entry>& found_entries, aurostd::xvector<double> normal, const double offset) :
      specie_vector(specie_vector), entries(found_entries), vertices(vertices), m_normal(normal), hyperplane_offset(offset) {
    this->m_dimension = specie_vector.size();
  }

  /// @brief lightweight constructor used initialize from a qhull facet during non-thermodynamic calcs with entry objects
  NhullFacet::NhullFacet(const std::vector<Entry>& found_entries, aurostd::xvector<double> normal) : entries(found_entries), m_normal(normal) {
    this->m_dimension = m_normal.rows;
  }

  /// @brief distance method from qhull copied here to allow easy post-processing distance calcs
  /// @note original qhull method can be found in extern/QHULL/libqhullcpp/QHULLHyperplane.cpp ln# 91
  /// @note   If greater than zero, the point is above the facet (i.e., outside).
  /// @return distance from point to hyperplane.
  /// @authors
  /// mod{NHA,20260206,created}
  double NhullFacet::distance(const aurostd::xvector<double>& p) const {
    int dimension = p.rows;

    if (dimension != m_dimension) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Input dimension mismatch!");
    }
    const double offset = hyperplane_offset;
    aurostd::xvector<double> point = p;
    point.shift(0);
    aurostd::xvector<double> normal = m_normal;
    normal.shift(0);

    double dist;

    switch (dimension) {
      case 2: dist = offset + point[0] * normal[0] + point[1] * normal[1]; break;
      case 3: dist = offset + point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2]; break;
      case 4: dist = offset + point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2] + point[3] * normal[3]; break;
      case 5: dist = offset + point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2] + point[3] * normal[3] + point[4] * normal[4]; break;
      case 6: dist = offset + point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2] + point[3] * normal[3] + point[4] * normal[4] + point[5] * normal[5]; break;
      case 7: dist = offset + point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2] + point[3] * normal[3] + point[4] * normal[4] + point[5] * normal[5] + point[6] * normal[6]; break;
      case 8:
        dist = offset + point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2] + point[3] * normal[3] + point[4] * normal[4] + point[5] * normal[5] + point[6] * normal[6] + point[7] * normal[7];
        break;
      default:
        dist = offset;
        for (int k = dimension; k--;) {
          //dist += *point++ * *normal++;
          dist += point[k] * normal[k];
        }
        break;
    }
    return dist;
  }//distance

  ///@brief calculates distance to hull (vertical) and phase decompositions for a thermodynamic hull calc
  ///@param entry
  ///@param hull_entries list of all entries on the hull
  ///@param nhull_facet_list list of all hull facets
  /// @authors
  /// mod{NHA,20260206,created}
  void NhullFacet::calc_vertical_distance_and_phase_decomposition(Entry& entry, const std::vector<Entry>& hull_entries, const vector<NhullFacet>& nhull_facet_list) {
    aurostd::xvector point = entry.m_coord;

    std::tuple<double, std::map<Entry, double>> nearest_dist_store = nearestDistDecomp(hull_entries, point, nhull_facet_list);
    const double NRG_store = std::get<0>(nearest_dist_store);
    std::map<Entry, double> decomp_map_store = std::get<1>(nearest_dist_store);

    entry.set_distance_hull_enthalpy_formation_atom(NRG_store);
    entry.m_nhull_phase_decomp_entries = decomp_map_store;
    entry.nhull_phase_decomp = entry.phaseDecompToMap();
  }

  ///@brief returns the hyperplane offset of a facet
  /// @authors
  /// mod{NHA,20260206,created}
  double NhullFacet::getOffset() {
    return -aurostd::scalar_product(vertices[0], m_normal);
  }

} // namespace nhull
