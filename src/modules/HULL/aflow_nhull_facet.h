#ifndef AFLOW_AFLOW_NHULL_FACET_H
#define AFLOW_AFLOW_NHULL_FACET_H

#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_nhull_entry.h"
#include "extern/QHULL/libqhullcpp/QhullFacet.h"

namespace nhull {

  // general struct for storing entries of 1 facet
  class NhullFacet {
  public:
      //vars:
    uint m_dimension;
    bool m_initialized;
    aurostd::xvector<double> m_normal;
    std::vector<std::string> specie_vector;
    std::vector<Entry> entries;
    std::vector<aurostd::xvector<double>> vertices; // full qhull coords here: includes enthalpy (but no central compound)
    std::unordered_set<std::string> auid_identifier; // stores all three auid's which uniquely define this facet
    uint f_index;   //indexing variable for performance and ease of calling: CORRESPONDS TO QHULL FACET LIST STORED IN CONVEXHULL class

      //QHULL specific vars (hyperplane parameters associated with this facet):
    double hyperplane_offset;

      //constructors:
    NhullFacet();
    NhullFacet(const std::vector<aurostd::xvector<double>>& vertices, aurostd::xvector<double> normal, const double offset, const int dimension);
    NhullFacet(const std::vector<Entry>& found_entries, aurostd::xvector<double> normal);
    NhullFacet(const std::vector<std::string>& specie_vector, const std::vector<aurostd::xvector<double>>& vertices, const std::vector<Entry>& found_entries, aurostd::xvector<double> normal, const double offset);
      //operators:
    bool operator==(const NhullFacet& other_facet) const { return auid_identifier == other_facet.auid_identifier; }

      //conversion functions:
    static void qhullRunToNhullFacets(std::vector<NhullFacet>& nhull_facets, orgQhull::Qhull& qhull, const std::vector<Entry>& entries);
    static void qhullFacetToNhullFacet(std::vector<NhullFacet>& nhull_facets, std::vector<Entry>& filtered_entries, const orgQhull::QhullFacet& facet, int dimension, const std::vector<std::string>& specie_vector);
    static void qhullFacetToNhullFacet(std::vector<NhullFacet>& nhull_facets, const std::vector<Entry>& hull_points, const orgQhull::QhullFacet& facet, int dimension);
    static std::vector<NhullFacet> getFacetData(const std::vector<std::string>& specie_vector, const std::vector<aurostd::xvector<double>>& lower_dim_point_list, std::vector<Entry>& filtered_entries, int dimension);

      //utility functions:
    static std::vector<aurostd::xvector<double>> qhullFacetToVertices(const orgQhull::QhullFacet& facet, bool remove_energy = false);
    static std::vector<aurostd::xvector<double>> nhullFacetToVertices(const NhullFacet& facet, bool remove_energy);
    void enthalpyUnitConversion();
    static void initializeHullEntry(Entry& hull_entry);

      //math functions
    static void calc_vertical_distance_and_phase_decomposition(Entry& entry, const std::vector<Entry>& hull_entries, const std::vector<NhullFacet>& nhull_facet_list);
    [[nodiscard]] double distance(const aurostd::xvector<double>& p) const; //CAREFUL: this is not vertical distance to facet! Just regular distance calc.
    double getOffset();

      //getters
    [[nodiscard]] aurostd::xvector<double> get_m_normal() const { return m_normal; }
    [[nodiscard]] std::vector<std::string> get_specie_vector() const { return specie_vector; }
    [[nodiscard]] std::vector<Entry> get_entries() const { return entries; }
    [[nodiscard]] std::vector<aurostd::xvector<double>> get_vertices() const { return vertices; }
    [[nodiscard]] std::unordered_set<std::string> get_auid_identifier() const { return auid_identifier; }
  };

  std::map<Entry, double> decompReactionEntries(const std::vector<Entry>& entries_on_hull, const NhullFacet& closest_facet, const aurostd::xvector<double>& qhull_point);
  std::map<std::string, double> decompReaction(const std::vector<std::string>& alloy_name_store, const orgQhull::QhullFacet& closest_facet, const aurostd::xvector<double>& qhull_point);
} // namespace nhull

#endif // AFLOW_AFLOW_NHULL_FACET_H
