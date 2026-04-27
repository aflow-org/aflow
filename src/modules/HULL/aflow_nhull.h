// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                                                                         *
// ***************************************************************************
// Written by Nick Anderson
// nicholas.anderson@duke.edu
// Previous versions also written by Eric Perim, Eric Gossett, and Corey Oses

#ifndef AFLOW_NHULL_H
#define AFLOW_NHULL_H

#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_nhull_entry.h"
#include "aflow_nhull_facet.h"
#include "flow/aflow_xclasses.h"

#define NHULL_DEBUG false

// non-member functions in nhull.cpp that need to be used elsewhere:
namespace nhull {
  // Just another group of Entries
  class CoordGroup {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    CoordGroup();
    // constructors - STOP
    bool operator<(const CoordGroup& other) const;
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // initializer
    void setCoords(const aurostd::xvector<double>& coord, bool has_stoich_coords);

    // attributes
    bool m_initialized;
    aurostd::xvector<double> m_coord;
    std::vector<uint> m_points; // points to nhull::Entry
    bool m_has_stoich_coords; //if true, then m_coord represents a coordinate with the first element index included
    bool m_has_artificial_unary;
    bool m_is_on_hull; // is ground state
    uint m_hull_member; // points to nhull::Entry
    uint m_ref_state; // first in order (artificial map)
    std::vector<uint> m_candidate_hull_points; // also points to nhull::Entry, refers to extremes of CoordGroup (half_hull's only has extrema)

    // for organization of points
    uint m_i_nary; // stoich_coords only
    uint m_i_alloy; // stoich_coords_only

    // the follow properties are ONLY found if stoich_coords + half_hull + artificial_points
    // assumes ONE gstate per coordgroup
    uint m_nearest_facet; // which facet is directly below/above me?
    double m_nearest_distance; // how close is the nearest facet?
    std::vector<uint> m_decomp_phases;
    aurostd::xvector<double> m_decomp_coefs;
    std::vector<std::vector<uint>> m_equilibrium_phases;
    bool m_calculated_equivalent_g_states;
    std::vector<uint> m_equivalent_g_states; // structure comparison
    std::vector<uint> m_sym_equivalent_g_states; // structure comparison
    double m_stability_criterion; // g-states only
    double m_n_plus_1_enthalpy_gain; // g-states only
    bool m_found_icsd_g_state; // whether icsd exists among equivalent states
    uint m_i_icsd_g_state; // canonical icsd entry (lowest number)

    friend class ConvexHull; // ConvexHull defines everything!
  };
  class ConvexHull : public xStream {
  private:
      // attributes
    bool m_initialized;
    std::vector<std::string> m_velements;
    std::vector<Entry> m_points;
    std::vector<CoordGroup> m_coord_groups; // index to m_points
    uint m_dimension; // full dimensionality
    std::vector<NhullFacet> m_facets; //nhull native facets for thermodynamic properties

      // flags
    aurostd::xoption m_nflags; //general nhull run arguments
    _aflags m_aflags; // used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags
  public:
      //constructors:
    ConvexHull();
    ConvexHull(aurostd::xoption hull_options, std::vector<aurostd::xvector<double>> points);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<Entry>& vpoints, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow, const _aflags& aflags, std::vector<std::string> specie_vector, std::ostream& _oss = std::cout, bool add_artificial_unaries = true);

      //functions for running everything from user input parameters
    void mainRun();

      // wrappers for try/catch's
    bool createHull(const std::vector<aurostd::xvector<double>>& vcoords);
    bool createThermoHull();

      //hull calculators
    void calculateSimpleHull(); //a general qhull calc given a set of points
    void calculateThermoHull();
    std::vector<NhullFacet> calculateNminus1(std::vector<Entry> pre_flag_entries,
                                             const std::vector<std::string> specie_vector,
                                             const int dimension,
                                             const std::vector<Entry>& filtered_entries,
                                             const std::string& nminus1_auid,
                                             bool verbose = true);
    std::vector<NhullFacet> calculateStabilityCriterion(std::vector<Entry> entries, const std::vector<std::string>& specie_vector, int dimension, const std::vector<Entry>& filtered_entries, const std::string& auid, bool verbose = true);
    void calculatePOCC();

      //post-hull calc point initializers
    static void postIterateRunInitializeEntries(std::vector<Entry>& entries, std::vector<NhullFacet>& facets);
    void initializeCoordGroups(); //used only for GFA right now

      //utilities:
    std::vector<Entry> getEntriesFromAFLOWlib(std::vector<std::string> specie_vector, const aurostd::xoption& vpflow);
    void setDirectory(std::string directory = "");
    void setDirectory(aurostd::xoption vpflow);
    void setHullFlags(const aurostd::xoption& vpflow);
    void clear();
    void saveHullJSON(const std::string& file_name) const;
    static void savePlotAndTable(const std::string& file_name, const std::vector<Entry>& entries, const std::vector<std::string>& specie_vector, const std::vector<NhullFacet>& facets, const _aflags aflags, aurostd::xoption nflags);
    static std::vector<Entry> convertUnitsAllEntries(const std::vector<Entry>& entries);
    static std::vector<NhullFacet> convertUnitsAllFacets(const std::vector<NhullFacet>& facets);
    [[nodiscard]] std::vector<std::string> alloyToElements(uint i_nary, uint i_alloy) const;

    // getters for member vars:
    [[nodiscard]] bool isHullInitialized() const { return m_initialized; }
    [[nodiscard]] std::vector<std::string> getVelements() const { return m_velements; }
    [[nodiscard]] std::vector<Entry> getPoints() const { return m_points; }
    [[nodiscard]] std::vector<CoordGroup> getCoordGroups() const { return m_coord_groups; }
    [[nodiscard]] uint getDimension() const { return m_dimension; }
    [[nodiscard]] std::vector<NhullFacet> getFacets() const { return m_facets; }
    [[nodiscard]] aurostd::xoption getNflags() const { return m_nflags; }
    [[nodiscard]] _aflags getAflags() const { return m_aflags; }

    //getters for specific entries in m_points:
    [[nodiscard]] double getDistanceToHull(Entry& entry, bool redo = false, bool get_signed_distance = false) const;
    [[nodiscard]] double getDistanceToHull(const std::string auid) const;
    void distanceToHullPrintToConsole(const std::string auid) const;
    [[nodiscard]] std::vector<uint> getDecompositionPhases(Entry& entry) const;
    [[nodiscard]] static bool isHalfHull(const std::vector<Entry>& filtered_entries, bool pocc_run);
    [[nodiscard]] double getStabilityCriterion(const std::string& cauid) const;
    [[nodiscard]] aurostd::xvector<double> getDecompositionCoefficients(uint i_point) const;
    [[nodiscard]] double getNminus1EnthalpyGain(const std::string& auid) const;
    [[nodiscard]] uint getEntriesCount() const;
    [[nodiscard]] bool setCoordGroupIndex(const aurostd::xvector<double>& r_coords, uint& i_coord_group) const;

      //facet functions
    void sortFacetVertices(std::vector<uint>& facet, const uint& facet_id); // HE20210510
    void getJoinedFacets(std::vector<std::vector<uint>>& facet_collection, const double angle_threshold = PI / 180); // HE20210510 - 0.018 radians ~ 1 degree
    [[nodiscard]] std::vector<aurostd::xvector<double>> getAllCoordinates();

      //methods associated with points:
    void addArtificialUnaries();
    void initializePointsSimple(const std::vector<aurostd::xvector<double>>& vcoords);
    void initializePointsThermo(const std::vector<Entry>& vpoints, bool add_artificial_unaries = false);

    [[nodiscard]] uint artificialMap(uint i_point) const;
  };

  /// @brief struct containing all run flags that modify hull calculations
  const struct ConvexFlags {
    std::string m_pocc = "NHULL::POCC";
    std::string m_stability_criterion_flag = "NHULL::STABILITY_CRITERION";
    std::string m_alloy = "PFLOW::ALLOY";
    std::string m_nminus1_flag = "NHULL::N-1_ENTHALPY_GAIN";
    std::string m_print_3d = "NHULL::PRINT3D";
    std::string m_stats = "NHULL::STATS";
    std::string m_include_outliers = "NHULL::INCLUDE_OUTLIERS";
    std::string m_no_output = "NHULL:NO_OUTPUT";
    std::string m_nhull_print = "NHULL::OUTPUT";
    std::string m_multi_output = "NHULL::MULTI_OUTPUT";
    std::string m_print_JSON = "NHULL::JSON_DOC";
    std::string m_latex_doc = "NHULL::LATEX_DOC";
    std::string m_dhull = "NHULL::DIST2HULL";

    std::string m_interquartile_range = "NHULL::IQR";
    std::string read_cache = "NHULL::READ_CACHE";
    std::string write_cache = "NHULL::READ_CACHE";
    std::string save_qhull_facets = "NHULL::SAVE_QHULL_FACETS"; //explicitly save facets from calculations for distance calcs later
    std::string m_nhull_path = "NHULL::PATH";
    std::string m_screen_only = "NHULL::SCREEN_ONLY"; //TODO see if we need screen only
    std::string m_neglect_cce = "NHULL::NEGLECT_CCE";
  } convex_flags;

  // sorting structures
  struct sortWithinCoordGroup {
    sortWithinCoordGroup(const std::vector<Entry>& vpoints, bool sort_energy_ascending = true) : m_points(vpoints), m_sort_energy_ascending(sort_energy_ascending) {};
    const std::vector<Entry>& m_points;
    bool m_sort_energy_ascending; // good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
    bool operator()(uint ci, uint cj);
  };

  ConvexHull convexHull(const aurostd::xoption& vpflow);
  std::vector<Entry> lightFilter(const std::vector<Entry>& entries);
  void stabilityCriteriaFilter(std::vector<Entry>& entries, const std::string& excluded_auid);
  void Nminus1Filter(std::vector<Entry>& entries, uint N_compound);
  void flagCheck(aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss, bool silent);
  std::string getPath(bool add_backslash = true);
  std::string getPath(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool silent = true);  // CO20180220

} // namespace nhull

#endif
