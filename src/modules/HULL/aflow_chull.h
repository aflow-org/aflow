// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu
// Previous versions also written by Eric Perim and Eric Gossett

#ifndef _AFLOW_CHULL_H_
#define _AFLOW_CHULL_H_

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/types.h>

#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"

// UNITS
const char _std_ = 'S';  // standard aflow/vasp units
const char _m_ = 'm';    // convert to milli-

// FORMATS
//[ME20190628 - moved to aflow.h] const char _apool_ = 'a';  // apool
//[ME20190628 - moved to aflow.h] const char _json_ = 'j';   // standard json
//[ME20190628 - moved to aflow.h] const char _pdf_ = 'p';    // pdf
//[ME20190628 - moved to aflow.h] const char _txt_ = 't';    // plain text
//[ME20190628 - moved to aflow.h] const char _web_ = 'w';    // web json
//[ME20190628 - moved to aflow.h] const char _latex_ = 'l';    // latex
//[ME20190628 - moved to aflow.h] const char _gnuplot_ = 'g';  // gnuplot
//[ME20190628 - moved to aflow.h] const char _jupyterthree_ = 'y';  // jupyter python 3
//[ME20190628 - moved to aflow.h] const char _jupytertwo_ = 'z'; // jupyter python 2

// REDUCTION MODES
//[ME20190628 - moved to aflow.h] const char _frac_ = 'f';  //fractional
//[ME20190628 - moved to aflow.h] const char _gcd_ = 'g';   //gcd
//[ME20190628 - moved to aflow.h] const char _none_ = 'n';  //none

// DEFAULTS
const int CHULL_PRECISION = 8;                          // must be less than _precision_ in aflow_xatom.cpp, which is currently set to 14
const int FULL_PRECISION = 15;                          // max printing precision
const int COEF_PRECISION = 4;
const int MEV_PRECISION = 3;                            // precision to within 1 meV
const double ZERO_TOL = pow(10, -CHULL_PRECISION);       // lower bound for absolute resolution of floats, significant differences among floats should be well above this threshold
const double ROUNDOFF_TOL = pow(10, -CHULL_PRECISION + 2); // make less stringent so we don't get 1e-6
const double ZERO_FULL_TOL = pow(10, -FULL_PRECISION);
const double ZERO_COEF_TOL = pow(10, -COEF_PRECISION);
const double ZERO_MEV_TOL = pow(10, -MEV_PRECISION);
const double ENERGY_TOL = 0.015;                        // eV, CO NOTES - structures within this threshold may be equivalent, I've seen as large as 5meV, keep at 15 to be safe
const int ZERO_RANGE_TOL = 1;
//[CO20180316 - moved to aflowrc]const uint BINARY_ENTRIES_THRESHOLD = 200;

// CO20180419 - moved to AFLOWRuntimeError and AFLOWLogicError
// namespace chull {
//   class CHullRuntimeError : public std::runtime_error {
//     public:
//       CHullRuntimeError(const std::string& function,const std::string& message);
//       CHullRuntimeError(const std::string& function,std::stringstream& message);
//       string where();
//       ~CHullRuntimeError() throw() {};
//     private:
//       string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//   };
//   class CHullLogicError : public std::logic_error {
//     public:
//       CHullLogicError(const std::string& function,const std::string& message);
//       CHullLogicError(const std::string& function,std::stringstream& message);
//       string where();
//       ~CHullLogicError() throw() {};
//     private:
//       string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//   };
// } // namespace chull

namespace chull {
  class ChullPoint; // forward declaration
  class ConvexHull; // forward declaration
  bool convexHull(const aurostd::xoption& vpflow);
  bool convexHull(const std::string& input, const aurostd::xoption& vpflow, const _aflags& aflags, std::ostream& oss = std::cout, bool silence_flag_check = false);
  ////////////////////////////////////////////////////////////////////////////////
  // gets path to redirect output
  std::string getPath(bool add_backslash = true);
  std::string getPath(const aurostd::xoption& vpflow, std::ostream& oss = std::cout, bool silent = true); // CO20180220
  std::string getPath(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool silent = true);  // CO20180220
  std::string getPath(std::string _path, std::ostream& oss = std::cout, bool silent = true); // CO20180220
  std::string getPath(std::string _path, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool silent = true);  // CO20180220
  ////////////////////////////////////////////////////////////////////////////////
  // logs which flags are on
  void flagCheck(aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ostream& oss = std::cout, bool silent = false);
  void flagCheck(aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool silent = false);
  ////////////////////////////////////////////////////////////////////////////////
  // stability criterion calculation
  ////////////////////////////////////////////////////////////////////////////////
  // returns value in desired units
  double convertUnits(double value, char units = _std_);
  double H_f_atom(const ChullPoint& point, char units = _std_);
  double H_f_atom(const aflowlib::_aflowlib_entry& entry, char units = _std_);
  double T_S(const ChullPoint& point);
  double T_S(const aflowlib::_aflowlib_entry& entry);
  double EFA(const ChullPoint& point, char units = _std_);
  double EFA(const aflowlib::_aflowlib_entry& entry, char units = _std_);
  double isoMaxLatentHeat(const ChullPoint& point, double x, char units = _std_);
  double isoMaxLatentHeat(const aflowlib::_aflowlib_entry& entry, double x, char units = _std_);
  ////////////////////////////////////////////////////////////////////////////////
  bool subspaceBelongs(const aurostd::xvector<int>& space, const aurostd::xvector<int>& subspace);
  bool correctSignVerticalDistance(double dist_2_hull, bool should_be_positive);
  aurostd::xvector<double> getTruncatedCoords(const aurostd::xvector<double>& coords, const aurostd::xvector<int>& elements_present); // truncated arbitrary coords
  std::vector<uint> getRelevantIndices(const aurostd::xvector<int>& elements_present);
  bool coordsIdentical(const aurostd::xvector<double>& coords1, const aurostd::xvector<double>& coords2);  // CO20210315
} // namespace chull

// CO20180420 - moved to xStream (aflow.h)
// namespace chull {
//   class ChullClassTemplate {
//     public:
//       //NECESSARY PUBLIC CLASS METHODS - START
//       //constructors - START
//       ChullClassTemplate();
//       //constructors - STOP
//       ~ChullClassTemplate();
//       //NECESSARY PUBLIC CLASS METHODS - END
//
//       ////flags
//       //aurostd::xoption m_cflags;
//       //_aflags m_aflags;                  //used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags
//
//     protected:
//       //NECESSARY private CLASS METHODS - START
//       void free();
//       void freeStreams();
//       void freeAll();
//       //NECESSARY END CLASS METHODS - END
//
//       //logger variables
//       ostream* m_p_oss;
//       ofstream* m_p_FileMESSAGE;
//       bool m_new_ofstream;  //for deletion later
//
//       //general setters
//       //void setDefaultCFlags();
//       //void setCFlags(const aurostd::xoption& vpflow);
//       void setOFStream(ofstream& FileMESSAGE);
//       void setOSS(ostream& oss);
//       //void setDirectory();
//   };
// } // namespace chull

namespace chull {
  class ChullPointLight : virtual public xStream {
  public:
      // NECESSARY PUBLIC CLASS METHODS - START
      // constructors - START
    ChullPointLight(std::ostream& oss = std::cout);
    ChullPointLight(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    ChullPointLight(const ChullPointLight& b);  // upcasting is allowed, works for ChullPointLight and ChullPoint
      // constructors - STOP
    ~ChullPointLight();
    const ChullPointLight& operator=(const ChullPointLight& other); // upcasting is allowed, works for ChullPointLight and ChullPoint
    bool operator<(const ChullPointLight& other) const;
    void clear();
      // NECESSARY PUBLIC CLASS METHODS - STOP

      // general attributes
    bool m_initialized;
    aurostd::xvector<double> m_coords;           // most general hull coordinates
    bool m_has_stoich_coords;
    bool m_formation_energy_coord;
    bool m_is_artificial;
    bool m_has_entry;
    uint m_i_nary;    // stoich_coords only
    uint m_i_alloy;   // stoich_coords_only

      // temporary variables that come with every new calculate() command
    aurostd::xvector<double> h_coords;  // dummy coord that changes per convex hull iteration, mere projection

      // getters
    [[nodiscard]] double getLastCoord() const;

  protected:
      // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const ChullPointLight& b);  // upcasting is allowed, works for ChullPointLight and ChullPoint
      // NECESSARY END CLASS METHODS - END
  };
  class ChullPoint : public ChullPointLight {
  public:
      // NECESSARY PUBLIC CLASS METHODS - START
      // constructors - START
    ChullPoint(std::ostream& oss = std::cout, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    ChullPoint(const aurostd::xvector<double>& coord, std::ostream& oss = std::cout, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    ChullPoint(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, std::ostream& oss = std::cout, bool formation_energy_coord = true);
    ChullPoint(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    ChullPoint(const aurostd::xvector<double>& coord, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    ChullPoint(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool formation_energy_coord = true);
    ChullPoint(const ChullPoint& b);
      // constructors - STOP
    ~ChullPoint();
    const ChullPoint& operator=(const ChullPoint& other);
    void clear();
      // NECESSARY PUBLIC CLASS METHODS - STOP

      // general attributes
    aflowlib::_aflowlib_entry m_entry;

      // for organization of points
    uint m_i_coord_group;
    bool m_found_icsd;
    uint m_i_icsd;                    // index of equivalent icsd in m_points

      // stoich_coords only!
    aurostd::xvector<double> s_coords;         // stoich_coords, m_coords[0:rows-1]+hidden_dimension
    aurostd::xvector<double> c_coords;         // composition coords, similar to stoich_coords, but with integers (except for POCC)
    aurostd::xvector<int> m_elements_present;  // 1 if s_coords[i]>ZERO_TOL, zero otherwise, that way, nary=sum(m_elements_present)

      // calculate per hull
    bool m_calculated_equivalent_entries; // equivalent entries have been calculated
    std::vector<uint> m_equivalent_entries;    // equivalent entries
    bool m_is_on_hull;  // one max per coordgroup
    bool m_is_g_state;  // one max per coordgroup, must have m_entry (artificialMap())
    bool m_is_equivalent_g_state; // can be many, includes original g_state
    bool m_is_sym_equivalent_g_state; // can be many, includes original g_state
    double m_dist_2_hull; // warning, this is not true dist to hull (facet), this is vertical distance
    double m_stability_criterion;
    double m_n_plus_1_enthalpy_gain;

      // initialization methods
    bool initialize(std::ostream& oss, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    bool initialize(const aurostd::xvector<double>& coord, std::ostream& oss, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    bool initialize(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, std::ostream& oss, bool formation_energy_coord = true);
    bool initialize(std::ofstream& FileMESSAGE, std::ostream& oss, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    bool initialize(const aurostd::xvector<double>& coord, std::ofstream& FileMESSAGE, std::ostream& oss, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    bool initialize(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, std::ofstream& FileMESSAGE, std::ostream& oss, bool formation_energy_coord = true);
    bool initialize(bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    bool initialize(const aurostd::xvector<double>& coord, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false);
    bool initialize(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, bool formation_energy_coord = true);

      // getters
    [[nodiscard]] bool isWithinHalfHull(bool lower_hull = true) const;
    [[nodiscard]] bool isGState() const;
    [[nodiscard]] aurostd::xvector<double> getStoichiometricCoords() const; // get stoichiometric coordinates (sans energetic coordinate)
    [[nodiscard]] aurostd::xvector<double> getTruncatedReducedCoords(const aurostd::xvector<int>& elements_present, vector_reduction_type vred = frac_vrt) const;
    [[nodiscard]] aurostd::xvector<double> getTruncatedSCoords(const aurostd::xvector<int>& elements_present) const; // truncate stoichiometry
    [[nodiscard]] aurostd::xvector<double> getTruncatedCCoords(const aurostd::xvector<int>& elements_present, bool reduce = true) const; // similar to truncated stoichiometry, but in integer form (if not POCC)
    [[nodiscard]] aurostd::xvector<double> getReducedCCoords() const; // reduce by gcd()
    uint loadXstructures(bool relaxed_only);
    bool getMostRelaxedXstructure(xstructure& xstr) const;
    [[nodiscard]] uint getDim() const;
    [[nodiscard]] bool isUnary() const;
    [[nodiscard]] double getFormationEnthalpy() const;
    [[nodiscard]] double getEntropicTemperature() const;
    [[nodiscard]] double getEntropyFormingAbility() const;
    [[nodiscard]] const std::vector<std::string>& getVSG() const;
    [[nodiscard]] const std::string& getSG() const;
    [[nodiscard]] double getDist2Hull(char units = _std_) const;
    [[nodiscard]] double getStabilityCriterion(char units = _std_) const;
    [[nodiscard]] double getRelativeStabilityCriterion() const;
    [[nodiscard]] double getNPlus1EnthalpyGain(char units = _std_) const;
    [[nodiscard]] double getEntropyStabilizationCoefficient(char units = _std_) const;

    // setters
    void setHullCoords();
    void setHullCoords(const aurostd::xvector<double>& coords);
    void setHullCoords(const aurostd::xvector<int>& elements_present);
    void reduceCoords(const aurostd::xvector<int>& elements_present);

    // general methods
    [[nodiscard]] bool entryIdentical(const aflowlib::_aflowlib_entry& other) const;

    friend class ChullFacet; // ConvexHull defines everything!
    friend class ConvexHull; // ConvexHull defines everything!
  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const ChullPoint& b);
    // NECESSARY END CLASS METHODS - END

    // initialization methods
    void initializeCoords(const aurostd::xvector<double>& coord, bool formation_energy_coord = false);
    void initializeCoords(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, bool formation_energy_coord = true);
    void addEntry(const aflowlib::_aflowlib_entry& entry);
    void setGenCoords(const aurostd::xvector<double>& coord, bool formation_energy_coord = false);
    void setGenCoords(const std::vector<std::string>& velements, const aflowlib::_aflowlib_entry& entry, bool formation_energy_coord = true);

    // hull specific methods
    [[nodiscard]] std::vector<uint> getRelevantIndices(const aurostd::xvector<int>& elements_present) const;
    void setStoichCoords();
    void cleanPointForHullCalc(); // clean hull state of point
    void cleanPointForHullTransfer();
  };
} // namespace chull

namespace chull {
  // allows chullfacet to be (mostly) INDEPENDENT of convex hull
  // there is a HUGE advantage to considering points by index, rather than coordinate
  // comparisons/groupings are easier
  // make sure that such comparisons are being made per facets of the same hull, never among different hulls!
  class FacetPoint {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    FacetPoint();
    FacetPoint(const ChullPointLight& point, uint index);
    FacetPoint(const FacetPoint& b);
    // constructors - STOP
    ~FacetPoint();
    const FacetPoint& operator=(const FacetPoint& other);
    bool operator<(const FacetPoint& other) const;
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // attributes
    bool m_initialized;
    uint ch_index; // belongs to chull
    ChullPointLight ch_point; // take only what you need from chullpoint, don't copy the whole thing (has entry which is large)

    // initializer
    void initialize(const ChullPointLight& point, uint index);

  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const FacetPoint& b);
    // NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

// nice sorting for points we know have energy vs. stoich coords
namespace chull {
  struct sortThermoPoints {
    sortThermoPoints(bool sort_stoich_ascending = true, bool sort_energy_ascending = true) : m_sort_stoich_ascending(sort_stoich_ascending), m_sort_energy_ascending(sort_energy_ascending) {};
    bool m_sort_stoich_ascending; // good for sorting points in facets (lower vs. upper hemisphere), if facet in lower hemisphere, then sort ascending order (CLOCKWISE, graphing)
    bool m_sort_energy_ascending; // good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
    bool operator()(const FacetPoint& fpi, const FacetPoint& fpj) const;
    bool operator()(const ChullPointLight& ci, const ChullPointLight& cj) const; // upcasting is allowed, works for ChullPointLight and ChullPoint
  };
  //[CO20221111 - VERY SLOW]struct sortLIB2Entries: public xStream { //fixes issue with AFLUX
  //[CO20221111 - VERY SLOW]  sortLIB2Entries(const vector<string>& velements,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_velements(velements) {;}
  //[CO20221111 - VERY SLOW]  ~sortLIB2Entries() {xStream::free();}
  //[CO20221111 - VERY SLOW]  const vector<string>& m_velements;
  //[CO20221111 - VERY SLOW]  bool operator() (const aflowlib::_aflowlib_entry i_entry,const aflowlib::_aflowlib_entry j_entry) const;
  //[CO20221111 - VERY SLOW]};
  struct _aflowlib_entry_LIB2sorting : public xStream { // for sorting LIB2 unaries
    _aflowlib_entry_LIB2sorting(aflowlib::_aflowlib_entry& entry, uint index, std::vector<std::string>& velements_chull, std::ofstream& FileMESSAGE, std::ostream& oss);
    ~_aflowlib_entry_LIB2sorting();
    //
    //[CO20221112 - no pointers]aflowlib::_aflowlib_entry* m_entry; //BE CAREFUL: this is a pointer, do not use out of scope
    std::string m_auid;
    std::string m_aurl;
    std::string m_catalog;
    std::string m_prototype;
    uint m_nspecies;
    uint m_index;
    std::vector<std::string> m_velements_chull; // BE CAREFUL: this is a pointer, do not use out of scope
    std::vector<std::string> m_species_AURL;
    bool operator<(const _aflowlib_entry_LIB2sorting& other) const;
  };
} // namespace chull

namespace chull {
  class ChullFacet : public xStream {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    ChullFacet(std::ostream& oss = std::cout);
    ChullFacet(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    ChullFacet(const ChullFacet& b);
    // constructors - STOP
    ~ChullFacet();
    const ChullFacet& operator=(const ChullFacet& other);
    bool operator<(const ChullFacet& other) const;
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // general attributes
    bool m_initialized;
    std::vector<FacetPoint> m_vertices;
    uint m_dim;
    bool m_has_stoich_coords;
    bool m_formation_energy_coord;
    double m_content; // length of line, area of triangle, vol of tetrahedron...
    std::vector<aurostd::xvector<double>> m_directive_vectors;
    aurostd::xvector<double> m_normal;
    double m_offset;
    aurostd::xvector<double> m_facet_centroid;
    aurostd::xvector<double> m_hull_reference;
    bool m_is_hypercollinear;
    bool m_is_vertical;
    bool m_is_artificial;
    bool m_in_lower_hemisphere; // this determines how we sort wrt stoich (descending in lower_hemisphere)
    std::vector<ChullFacet> m_ridges;

    // flags - MOVED TO xStream
    // aurostd::xoption m_cflags;
    //_aflags m_aflags;                  //used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags

    // getters
    [[nodiscard]] std::vector<uint> getCHIndices() const;

    // setters
    void addVertex(const ChullPointLight& point, uint index = AUROSTD_MAX_UINT);
    void addVertex(const FacetPoint& fp);
    void initialize(const aurostd::xvector<double>& hull_centroid, uint hull_dim, bool check_validity = true); // needs some hull data (context) to align correctly

    // methods using facets to build hull
    [[nodiscard]] bool shareRidge(const ChullFacet& other) const;
    [[nodiscard]] bool isPointOnFacet(const FacetPoint& fp) const;
    [[nodiscard]] bool isPointOnFacet(uint i_point) const;
    [[nodiscard]] bool isPointOutside(const FacetPoint& fp) const;
    [[nodiscard]] bool isPointOutside(const ChullPointLight& point) const;
    [[nodiscard]] double getSignedPointPlaneDistance(const ChullPointLight& point) const;
    [[nodiscard]] double getSignedPointPlaneDistance(const aurostd::xvector<double>& point) const;
    [[nodiscard]] double getSignedVerticalDistanceToZero(const ChullPointLight& point) const;
    [[nodiscard]] double getSignedVerticalDistanceToZero(const aurostd::xvector<double>& point) const;
    [[nodiscard]] double getSignedVerticalDistance(const ChullPointLight& point) const;
    [[nodiscard]] double getSignedVerticalDistance(const aurostd::xvector<double>& point) const;

    friend class ConvexHull; // ConvexHull defines everything!
  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const ChullFacet& b);
    // NECESSARY END CLASS METHODS - END

    // logger variables - MOVED TO xStream
    // ostream* m_p_oss;
    // ofstream* m_p_FileMESSAGE;
    // bool m_new_ofstream;  //for deletion later

    // temporary variables that come with every new calculate() command
    bool f_visited; // temporary per calculation(), has it been visited?
    std::vector<FacetPoint> f_outside_set;
    FacetPoint f_furthest_point;
    std::vector<uint> f_neighbors;

    // general setters
    // MOVED TO xStream
    // void setCFlags(const aurostd::xoption& vpflow);
    // void setOFStream(ofstream& FileMESSAGE);
    // void setOSS(ostream& oss);

    // validation methods - we split up as we build facets in pieces, and not all work
    // sometimes we find a (d-1)flat facet, so we remove
    // we want to return false vs. throw exception appropriately
    bool hasValidPoints(std::string& error);
    void setContent();
    void setDirectiveVectors(bool check_validity = true);
    bool pointsMatchDirectiveVectors(std::string& error);
    bool hasValidDirectiveVectors(std::string& error);
    bool hasCollinearVectors(bool check_validity = true);
    bool isValid(std::string& error);

    // initialization methods, again split up because we build in pieces
    void setNormal(bool check_validity = true);
    void setOffset();
    void setCentroid();
    void setVertical();
    void setArtificial();
    void alignNormalInward();
    void setHemisphere();
    void setFurthestPoint();
    void setRidges();

    // hull specific
    void cleanFacet(); // clean state of facet
  };
} // namespace chull

// Just another group of ChullPoints
namespace chull {
  class CoordGroup {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    CoordGroup();
    CoordGroup(const aurostd::xvector<double>& coord, bool has_stoich_coords);
    CoordGroup(const CoordGroup& b);
    // constructors - STOP
    ~CoordGroup();
    const CoordGroup& operator=(const CoordGroup& other);
    bool operator<(const CoordGroup& other) const;
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // initializer
    void initialize(const aurostd::xvector<double>& coord, bool has_stoich_coords);

    // getters
    [[nodiscard]] aurostd::xvector<int> getElementsPresent() const;
    [[nodiscard]] uint getDim() const;

    // attributes
    bool m_initialized;
    aurostd::xvector<double> m_coords; // ChullPoint.getStoichiometricCoords()
    std::vector<uint> m_points; // points to ChullPoints
    bool m_has_stoich_coords;
    bool m_has_artificial_unary;
    bool m_is_on_hull; // is ground state
    uint m_hull_member; // points to ChullPoints
    uint m_ref_state; // first in order (artificial map)
    std::vector<uint> m_candidate_hull_points; // also points to ChullPoints, refers to extremes of CoordGroup (half_hull's only has extrema)

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
  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const CoordGroup& b);
    // NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

// simple, nice organization of points
// usually, we further organize points into CoordGroups EXCEPT for hull_points_unary
// hull_points_unary === single CoordGroup, trivial
namespace chull {
  class Alloy {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    Alloy();
    Alloy(const aurostd::xvector<int>& elements_present);
    Alloy(const Alloy& b);
    // constructors - STOP
    ~Alloy();
    const Alloy& operator=(const Alloy& other);
    bool operator<(const Alloy& other) const;
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // initializer
    void initialize(const aurostd::xvector<int>& elements_present);
    [[nodiscard]] bool belongs2Hull(const aurostd::xvector<int>& elements_present_hull) const; // checks if alloy is a member of the hull (via elements_present_hull)

    // attributes
    bool m_initialized;
    aurostd::xvector<int> m_elements_present;
    uint m_dim; // projected hull dimensionality
    std::vector<uint> m_coord_groups;
    std::vector<uint> m_facets;

    friend class ConvexHull; // ConvexHull defines everything!
  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const Alloy& b);
    // NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

// simple, nice organization of points
namespace chull {
  class Nary {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    Nary();
    Nary(uint nary);
    Nary(const Nary& b);
    // constructors - STOP
    ~Nary();
    const Nary& operator=(const Nary& other);
    bool operator<(const Nary& other) const;
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // initializer
    void initialize(uint nary);

    // attributes
    bool m_initialized;
    uint nary;
    std::vector<Alloy> m_alloys;

    friend class ConvexHull; // ConvexHull defines everything!
  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const Nary& b);
    // NECESSARY private CLASS METHODS - STOP
  };
} // namespace chull

namespace chull {
  class ConvexHull : public xStream {
  public:
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    ConvexHull(std::ostream& _oss = std::cout);
    ConvexHull(const std::string& alloy, std::ostream& _oss = std::cout);
    ConvexHull(const std::vector<std::string>& velements, std::ostream& _oss = std::cout);
    ConvexHull(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& _oss = std::cout);
    ConvexHull(const std::vector<aurostd::xvector<double>>& vcoords, std::ostream& _oss = std::cout, bool has_stoich_coords = false, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const std::vector<ChullPoint>& vpoints, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const std::string& alloy, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const std::vector<aurostd::xvector<double>>& vcoords, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout, bool has_stoich_coords = false, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const std::vector<ChullPoint>& vpoints, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::string& alloy, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<aurostd::xvector<double>>& vcoords, std::ostream& _oss = std::cout, bool has_stoich_coords = false, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::string& alloy, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout);
    ConvexHull(const aurostd::xoption& vpflow,
               const std::vector<aurostd::xvector<double>>& vcoords,
               std::ofstream& FileMESSAGE,
               std::ostream& _oss = std::cout,
               bool has_stoich_coords = false,
               bool formation_energy_hull = false,
               bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, std::ofstream& FileMESSAGE, std::ostream& _oss = std::cout, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    ConvexHull(const aurostd::xoption& vpflow,
               const std::vector<ChullPoint>& vpoints,
               const std::vector<std::string>& velements,
               std::ofstream& FileMESSAGE,
               std::ostream& _oss = std::cout,
               bool formation_energy_hull = false,
               bool add_artificial_unaries = false);
    ConvexHull(const ConvexHull& b);
    // constructors - STOP
    ~ConvexHull();
    const ConvexHull& operator=(const ConvexHull& other);
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // attributes
    bool m_initialized;
    std::vector<std::string> m_velements;
    std::vector<uint> m_icsd_entries; // save these for checks for experimental validation
    std::vector<ChullPoint> m_points;
    std::vector<Nary> m_naries; // one really nice container for indicies for points, naries, then alloys, then CoordGroups, then points
    std::vector<CoordGroup> m_coord_groups; // index to m_points
    uint m_dim; // full dimensionality
    bool m_half_hull; // huge speed up
    bool m_lower_hull; // true == lower half only, false == upper half only
    bool m_has_stoich_coords; //[0,1]
    bool m_add_artificial_unaries; // force unaries to 0
    bool m_thermo_hull; // enthalpy_formation/entropic_temperature with stoich_coords and artificial unaries
    bool m_formation_energy_hull;
    std::vector<ChullFacet> m_facets;
    std::vector<uint> m_i_facets; // contains indices of most recently calculated facets, points to m_facets
    bool m_sort_energy_ascending; // true for lower half hulls

    // flags
    aurostd::xoption m_cflags;
    _aflags m_aflags; // used PURELY for the logger (path), so no need to pass into constructor, pull from m_cflags

    // initialization methods
    bool initialize(std::ostream& oss);
    bool initialize(const std::string& alloy, std::ostream& oss);
    bool initialize(const std::vector<std::string>& velements, std::ostream& oss);
    bool initialize(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& oss);
    bool initialize(const std::vector<aurostd::xvector<double>>& vcoords, std::ostream& oss, bool has_stoich_coords, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const std::vector<ChullPoint>& vpoints, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const std::string& alloy, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const std::vector<aurostd::xvector<double>>& vcoords, std::ofstream& FileMESSAGE, std::ostream& oss, bool has_stoich_coords, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const std::vector<ChullPoint>& vpoints, std::ofstream& FileMESSAGE, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize();
    bool initialize(const std::string& alloy);
    bool initialize(const std::vector<std::string>& velements);
    bool initialize(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries);
    bool initialize(const std::vector<aurostd::xvector<double>>& vcoords, bool has_stoich_coords, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const std::vector<ChullPoint>& vpoints, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::string& alloy, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<aurostd::xvector<double>>& vcoords, std::ostream& oss, bool has_stoich_coords, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::string& alloy, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<aurostd::xvector<double>>& vcoords, std::ofstream& FileMESSAGE, std::ostream& oss, bool has_stoich_coords, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, std::ofstream& FileMESSAGE, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, std::ofstream& FileMESSAGE, std::ostream& oss, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow);
    bool initialize(const aurostd::xoption& vpflow, const std::string& alloy);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<std::string>& velements);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<aurostd::xvector<double>>& vcoords, bool has_stoich_coords, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, bool formation_enthalpy_hull, bool add_artificial_unaries);
    bool initialize(const aurostd::xoption& vpflow, const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, bool formation_enthalpy_hull, bool add_artificial_unaries);

    // initialize points ONLY
    void initializePoints(const std::string& alloy);
    void initializePoints(const std::vector<std::string>& velements);
    void initializePoints(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries);
    void initializePoints(const std::vector<aurostd::xvector<double>>& vcoords, bool has_stoich_coords = false, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    void initializePoints(const std::vector<ChullPoint>& vpoints, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    void initializePoints(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, bool formation_energy_hull = false, bool add_artificial_unaries = false);

    // getters
    [[nodiscard]] uint getDim() const;
    [[nodiscard]] uint getEntriesCount(bool only_within_half_hull = true) const;
    [[nodiscard]] uint getEntriesCount(uint i_nary, bool only_within_half_hull = true) const;
    [[nodiscard]] uint getEntriesCount(uint i_nary, uint i_alloy, bool only_within_half_hull = true) const;
    [[nodiscard]] std::vector<std::vector<uint>> getHullSizes(bool only_within_half_hull = true) const;
    [[nodiscard]] uint getGStateCount() const;
    [[nodiscard]] uint getGStateCount(uint i_nary) const;
    [[nodiscard]] std::vector<uint> getHullPoints(bool sort_stoich_ascending = true) const;
    [[nodiscard]] std::vector<uint> getGStates(bool include_unaries = true, bool sort_stoich_ascending = true) const;
    [[nodiscard]] uint getUnaryGState(uint i_alloy) const;
    [[nodiscard]] bool isViablePoint(uint i_point) const;
    [[nodiscard]] bool isViableGState(uint g_state) const;
    bool findPoint(const std::string& auid, uint& i_point) const;
    bool findPoint(const aurostd::xvector<double>& coords, uint& i_point) const;
    bool getNariesIndex(uint i_point, uint& i_nary, uint& i_alloy, uint& i_coord_group, bool redo = false) const;
    bool getNariesIndex(const ChullPoint& point, uint& i_nary, uint& i_alloy, uint& i_coord_group, bool redo = false) const;
    bool getCoordGroupIndex(uint i_point, uint& i_coord_group, bool redo = false) const;
    bool getCoordGroupIndex(const ChullPoint& point, uint& i_coord_group, bool redo = false) const;
    bool getCoordGroupIndex(const aurostd::xvector<double>& r_coords, uint& i_coord_group) const;
    bool getAlloyIndex(const ChullPoint& point, uint& i_nary, uint& i_alloy, bool redo = false) const;
    bool getAlloyIndex(const CoordGroup& cg, uint& i_nary, uint& i_alloy, bool redo = false) const;
    bool getAlloyIndex(const aurostd::xvector<int>& elements_present, uint& i_nary, uint& i_alloy) const;
    [[nodiscard]] uint artificialMap(uint i_point) const;
    [[nodiscard]] uint getNearestFacetVertically(const std::vector<uint>& i_facets, const ChullPoint& point) const;
    [[nodiscard]] uint getNearestFacetVertically(const std::vector<uint>& i_facets, const aurostd::xvector<double>& point) const;
    [[nodiscard]] double getSignedVerticalDistanceWithinCoordGroup(uint i_coord_group, uint i_point) const;
    [[nodiscard]] double getSignedVerticalDistanceWithinCoordGroup(uint i_coord_group, const ChullPoint& point) const;
    [[nodiscard]] double getDistanceToHull(uint i_point, bool redo = false, bool get_signed_distance = false) const; // CO20190808
    [[nodiscard]] double getDistanceToHull(const ChullPoint& point, bool redo = false, bool get_signed_distance = false) const; // CO20190808
    [[nodiscard]] std::vector<double> getDistancesToHull(const std::vector<std::string>& vauid, bool redo = false) const;
    [[nodiscard]] std::vector<uint> extractDecompositionPhases(const ChullFacet& facet) const;
    [[nodiscard]] std::vector<uint> getDecompositionPhases(uint i_point) const;
    [[nodiscard]] std::vector<uint> getDecompositionPhases(const ChullPoint& point) const;
    [[nodiscard]] aurostd::xvector<double> getDecompositionCoefficients(uint i_point, vector_reduction_type vred = frac_vrt) const;
    [[nodiscard]] aurostd::xvector<double> getDecompositionCoefficients(const ChullPoint& point, vector_reduction_type vred = frac_vrt) const;
    [[nodiscard]] aurostd::xvector<double> getDecompositionCoefficients(uint i_point, const std::vector<uint>& decomp_phases, vector_reduction_type vred = frac_vrt) const;
    [[nodiscard]] aurostd::xvector<double> getDecompositionCoefficients(const ChullPoint& point, const std::vector<uint>& decomp_phases, vector_reduction_type vred = frac_vrt) const;
    [[nodiscard]] std::vector<uint> getAdjacentFacets(uint hull_member, bool ignore_hypercollinear = true, bool ignore_vertical = true, bool ignore_artificial = true) const;
    [[nodiscard]] std::vector<std::vector<uint>> getEquilibriumPhases(uint hull_member) const;
    [[nodiscard]] std::vector<uint> getEquivalentEntries(uint i_point) const;
    [[nodiscard]] std::vector<uint> getSymEquivalentGStates(uint g_state) const;
    [[nodiscard]] double getStabilityCriterion(const std::string& cauid) const;
    [[nodiscard]] double getStabilityCriterion(uint cpoint) const;
    [[nodiscard]] std::vector<double> getStabilityCriterion(const std::vector<std::string>& vcauid) const;
    [[nodiscard]] std::vector<double> getStabilityCriterion(const std::vector<uint>& vcpoints) const;
    void getFakeHull(const std::vector<uint>& vcpoints, ConvexHull& fake_hull) const;
    void getFakeHull(const std::vector<uint>& vcpoints, const aurostd::xvector<int>& elements_present, ConvexHull& fake_hull) const;
    void getFakeHull(const std::vector<uint>& vcpoints, const aurostd::xvector<int>& elements_present_hull, const aurostd::xvector<int>& elements_present_points, ConvexHull& fake_hull) const;
    [[nodiscard]] double getNPlus1EnthalpyGain(const std::string& cauid) const;
    double getNPlus1EnthalpyGain(const std::string& cauid, ConvexHull& fake_hull, bool hull_set) const;
    [[nodiscard]] std::vector<double> getNPlus1EnthalpyGain(const std::vector<std::string>& vcauid) const;
    std::vector<double> getNPlus1EnthalpyGain(const std::vector<std::string>& vcauid, ConvexHull& fake_hull, bool hull_set) const;
    [[nodiscard]] double getNPlus1EnthalpyGain(uint cpoint) const;
    double getNPlus1EnthalpyGain(uint cpoint, ConvexHull& fake_hull, bool hull_set) const;
    [[nodiscard]] std::vector<double> getNPlus1EnthalpyGain(const std::vector<uint>& vcpoints) const;
    std::vector<double> getNPlus1EnthalpyGain(const std::vector<uint>& vcpoint, ConvexHull& fake_hull, bool hull_set) const;
    void getJoinedFacets(std::vector<std::vector<uint>>& facet_collection, const double angle_threshold = PI / 180); // HE20210510 - 0.018 radians ~ 1 degree

    // writer
    [[nodiscard]] bool write(filetype ftype = latex_ft) const;

  private:
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const ConvexHull& b);
    // NECESSARY END CLASS METHODS - END

    // logger variables - MOVED TO xStream
    // ostream* m_p_oss;
    // ofstream* m_p_FileMESSAGE;
    // bool m_new_ofstream;  //for deletion later

    // aflowrc definitions we don't want to redefine over and over again, just once at initialization
    bool m_allow_all_formation_energies;
    std::vector<std::string> m_allowed_dft_types;

    // temporary variables that come with every new calculate() command
    uint h_dim; // projected hull dimensionality
    aurostd::xvector<int> m_elements_present; // stoich_coords only
    std::vector<uint> h_points; // current hull points, subset of m_points
    aurostd::xvector<double> h_centroid; // centroid of ALL these points
    aurostd::xvector<double> h_reference; // centroid of INITIAL points on hull (more stable reference)
    std::vector<ChullFacet> h_facets; // current hull facets, real facets are stored in Alloys (m_naries)
    std::vector<uint> h_visible_facets;
    std::vector<ChullFacet> h_horizon_ridges;

    // general setters
    void setDefaultCFlags();
    void setCFlags(const aurostd::xoption& vpflow);
    void setDirectory();
    // MOVED TO xStream
    // void setOFStream(ofstream& FileMESSAGE);
    // void setOSS(ostream& oss);

    // wrapper for try/catch's
    bool createHull(const std::string& alloy);
    bool createHull(const std::vector<std::string>& velements);
    bool createHull(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries);
    bool createHull(const std::vector<aurostd::xvector<double>>& vcoords, bool has_stoich_coords = false, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    bool createHull(const std::vector<ChullPoint>& vpoints, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    bool createHull(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, bool formation_energy_hull = false, bool add_artificial_unaries = false);

    // methods associated with points
    [[nodiscard]] bool entryValid(const aflowlib::_aflowlib_entry& entry, bool ignore_bad_database = true) const;
    bool entryValid(const aflowlib::_aflowlib_entry& entry, std::string& reason, bool ignore_bad_database = true) const;
    bool entryValid(const aflowlib::_aflowlib_entry& entry, std::string& reason, char& LOGGER_TYPE, bool ignore_bad_database = true) const;
    bool entryUnique(const std::vector<uint>& unique_entries, const aflowlib::_aflowlib_entry& entry2, std::string& canonical_auid) const;
    void addArtificialUnaries(uint dim);
    void addArtificialUnaries();
    void loadPoints(const std::string& alloy);
    void loadPoints(const std::vector<std::string>& velements);
    void loadPoints(const std::vector<std::string>& velements, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries);
    void loadPoints(const std::vector<aurostd::xvector<double>>& vcoords, bool has_stoich_coords = false, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    void loadPoints(const std::vector<ChullPoint>& vpoints, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    void loadPoints(const std::vector<ChullPoint>& vpoints, const std::vector<std::string>& velements, bool formation_energy_hull = false, bool add_artificial_unaries = false);
    void calculateOutlierThreshold(const std::vector<double>& _energies, double& upper_threshold, double& lower_threshold);
    void calculateOutlierThreshold(const aurostd::xvector<double>& energies, double& upper_threshold, double& lower_threshold);
    std::vector<uint> calculateOutliers(const std::vector<uint>& points_to_consider);
    std::vector<uint> getOutliers();
    std::vector<uint> getOutliers(const aurostd::xvector<int>& elements_present);
    std::vector<uint> findArtificialPoints(uint i_coord_group);
    uint findArtificialUnary(uint i_coord_group);
    void organizeHullPoints(uint i_coord_group);
    void organizeHullPoints();
    void initializeNaries();
    void structurePoints();
    [[nodiscard]] std::vector<std::string> alloyToElements(const ChullPoint& point) const;
    [[nodiscard]] std::vector<std::string> alloyToElements(uint i_nary, uint i_alloy) const;
    void checkStructurePoints();
    void sortFacetVertices(std::vector<uint>& facet, const uint& facet_id); // HE20210510

    // other methods associated with calculating hull
    void addPointToFacet(ChullFacet& facet, uint i_point);
    void initializeFacet(ChullFacet& facet, bool check_validity = true);
    bool isPointOutsideFacet(ChullFacet& facet, uint i_point);
    uint getExtremePoint(uint dim);
    uint getExtremePoint(uint dim, const std::vector<FacetPoint>& points_to_avoid);
    void setCentroid();
    std::vector<FacetPoint> getInitialExtremePoints();
    void setNeighbors();
    void createInitializeSimplex();
    void setVisibleFacets(uint i_facet);
    void setHorizonRidges();
    uint createNewFacets(FacetPoint furthest_point); // returns true if new facets were, in fact created (not necessary at each iteration)
    void updateOutsideSet(uint new_facet_count);
    void deleteVisibleFacets();
    void removeDuplicateHullPoints();
    void calculateFacets();
    [[nodiscard]] const aurostd::xvector<int>& getElementsPresent(uint i_nary, uint i_alloy) const;
    [[nodiscard]] const aurostd::xvector<int>& getElementsPresent(uint ipoint) const;
    [[nodiscard]] aurostd::xvector<int> getElementsPresent(const std::vector<uint>& vcpoints) const;
    void setElementsPresent(uint i_nary, uint i_alloy);
    void addRelevantUnariesToHullCalculation(uint i_nary, uint i_alloy);
    void addRelevantUnariesToHullCalculation(aurostd::xvector<int>& elements_present);
    void addLowerDimensionPointsToHullCalculation(uint i_nary_max);
    void addPointToHullCalculation(uint i_point);
    void addPointToHullCalculation(uint i_point, aurostd::xvector<int>& elements_present);
    void preparePointsForHullCalculation(uint i_nary, uint i_alloy);
    void preparePointsForHullCalculation();
    [[nodiscard]] const std::vector<uint>& getRelevantFacets(uint i_nary, uint i_alloy) const;
    void setHullMembers();
    void setHullMembers(uint i_nary, uint i_alloy);
    void setHullMembers(const std::vector<uint>& i_facets);
    void setNearestFacet(uint i_nary, uint i_alloy, uint i_coord_group);
    [[nodiscard]] bool hullDistanceExtractionRequired() const;
    [[nodiscard]] bool thermoPropertiesExtractionRequired() const;
    [[nodiscard]] bool thermoPostProcessingExtractionRequired() const;
    void setDistancesToHull(uint i_nary, uint i_alloy);
    void setDistancesToHull(uint i_nary, uint i_alloy, uint i_coord_group);
    void setDecompositionPhases(uint i_nary, uint i_alloy, uint i_coord_group);
    void setDecompositionCoefficients(uint i_nary, uint i_alloy, uint i_coord_group);
    void setOffHullProperties(uint i_nary, uint i_alloy);
    void setEquilibriumPhases(uint i_nary, uint i_alloy, uint i_coord_group);
    [[nodiscard]] bool phasesEquivalent(uint i_point1, uint i_point2, bool perform_structure_comparison) const;
    [[nodiscard]] bool energiesDiffer(uint i_point1, uint i_point2, bool strict = true) const;
    [[nodiscard]] bool spacegroupsDiffer(uint i_point1, uint i_point2, bool strict = true) const;
    [[nodiscard]] bool structuresEquivalent(uint i_point1, uint i_point2) const;
    bool isICSD(uint i_point, uint& i_point_icsd) const;
    void setEquivalentGStates(uint i_nary, uint i_alloy, uint i_coord_group);
    void setEquivalentStates(uint i_nary, uint i_alloy, uint i_coord_group);
    void setSymEquivalentGStates(uint i_nary, uint i_alloy, uint i_coord_group);
    void setOnHullProperties(uint i_nary, uint i_alloy);
    void storeHullData(uint i_nary, uint i_alloy);
    void storeHullData();
    void extractThermodynamicProperties(uint i_nary, uint i_alloy);
    void thermodynamicsPostProcessing();
    void calculate();
    void setStabilityCriterion();
    void setNPlus1EnthalpyGain(uint i_point);
    void setNPlus1EnthalpyGain(uint i_point, ConvexHull& fake_hull, bool hull_set);
    void setNPlus1EnthalpyGain();
    void cleanHull(); // clean state of hull

    // writer functions
    [[nodiscard]] std::string prettyPrintCompound(const ChullPoint& point, vector_reduction_type vred = gcd_vrt, bool exclude1 = true, filetype ftype = latex_ft) const;
    [[nodiscard]] std::string prettyPrintCompound(const aflowlib::_aflowlib_entry& entry, vector_reduction_type vred = gcd_vrt, bool exclude1 = true, filetype ftype = latex_ft) const;
    //[ME20190628 - moved to pflow.h] string prettyPrintCompound(const vector<string>& vspecies,const vector<double>& vcomposition,vector_reduction_type vred=gcd_vrt,bool exclude1=true,filetype ftype=latex_ft) const;
    //[ME20190628 - moved to pflow.h] string prettyPrintCompound(const vector<string>& vspecies,const aurostd::xvector<double>& vcomposition,vector_reduction_type vred=gcd_vrt,bool exclude1=true,filetype ftype=latex_ft) const;
    [[nodiscard]] std::string getICSDNumber(uint i_point, bool remove_suffix = true) const;
    [[nodiscard]] std::string getICSDNumber(const ChullPoint& point, bool remove_suffix = true) const;
    [[nodiscard]] std::string getICSDNumber(const aflowlib::_aflowlib_entry& entry, bool remove_suffix = true) const;
    [[nodiscard]] std::string prettyPrintPrototype(const ChullPoint& point, bool double_back_slash, bool icsd_label_skim = false) const;
    [[nodiscard]] std::string prettyPrintPrototype(const aflowlib::_aflowlib_entry& entry, bool double_back_slash, bool icsd_label_skim = false) const;
    //[CO20190419 - moved to aurostd_main.cpp]string fixStringLatex(const string& input, bool double_back_slash,bool symmetry_string) const;
    [[nodiscard]] std::string getPlotHeaderPDF(char function_mode, const std::string& column_header, bool points_color_gradient = DEFAULT_CHULL_LATEX_COLOR_GRADIENT) const;
    [[nodiscard]] std::string getPlotPointContentPDF(const ChullPoint& point, bool zero_end_point = true, bool zero_dist_2_hull = false) const;
    [[nodiscard]] std::string getNodeCoordPosition(const ChullPoint& point) const;
    [[nodiscard]] std::string getNodeCoordPosition(const aflowlib::_aflowlib_entry& entry, const aurostd::xvector<double>& coord) const;
    std::string nodeCreator(std::stringstream& option, std::stringstream& position, std::stringstream& content) const;
    [[nodiscard]] std::string nodeCreator(const std::string& option, const std::string& position, const std::string& content) const;
    [[nodiscard]] bool unwantedFacetLine(uint vi, uint vj, bool check_border = true) const;
    bool unwantedFacetLine(uint vi, uint vj, std::vector<std::vector<uint>>& facet_lines, bool check_border = true) const;
    [[nodiscard]] std::string getPointsPropertyHeaderList(filetype ftype) const;
    [[nodiscard]] std::string getDelta(bool helvetica_font) const;
    [[nodiscard]] std::string getStabilityCriterionSymbol(bool helvetica_font) const;
    [[nodiscard]] std::string getSnapshotTableHeader(std::string headers, bool designate_HEADER = false) const;
    [[nodiscard]] bool addInternalHyperlinks(bool internal_links_graph2report = true, bool internal_links_withinreport = true) const;
    [[nodiscard]] bool addExternalHyperlinks() const;
    [[nodiscard]] double getRoundToValue(double point_range) const;
    [[nodiscard]] double getYTickDistance(double y_range, int approx_num_ticks, double round_to_value) const;
    [[nodiscard]] std::vector<std::string> grabAcceptableLatexColors(bool replace_pranab_standard = true, bool allow_dvips_colors = true, uint count = 10) const;
    [[nodiscard]] std::vector<std::string> grabAcceptableLatexColors(const std::string& banned_colors_str, bool replace_pranab_standard = true, bool allow_dvips_colors = true, uint count = 10) const;
    [[nodiscard]] std::vector<std::string> grabAcceptableLatexColors(const std::vector<std::string>& banned_colors, bool replace_pranab_standard = true, bool allow_dvips_colors = true, uint count = 10) const;
    [[nodiscard]] aurostd::xoption resolvePlotLabelSettings() const;
    void writeLatex() const;
    [[nodiscard]] std::string getPlainTextHeader() const;
    [[nodiscard]] std::string getJSONHeader() const;
    [[nodiscard]] std::string grabCHPointProperty(const ChullPoint& point, const std::string& property, filetype ftype = txt_ft) const;
    [[nodiscard]] std::string grabCHFacetProperty(const ChullFacet& facet, const std::string& property, filetype ftype = txt_ft) const;
    std::vector<std::vector<std::string>> getPointsData(const std::string& properties_str, std::vector<std::string>& headers, filetype ftype = txt_ft) const;
    std::vector<std::vector<std::vector<std::vector<std::string>>>> getFacetsData(const std::string& properties_str, std::vector<std::string>& headers, filetype ftype = txt_ft) const;
    void getPlainTextColumnSizes(const std::vector<std::string>& headers, const std::vector<std::vector<std::string>>& ventry, std::vector<uint>& sizes) const;
    void getPlainTextColumnSizesPoints(const std::vector<std::string>& headers, const std::vector<std::vector<std::string>>& ventries, std::vector<uint>& sizes) const;
    void getPlainTextColumnSizesFacets(const std::vector<std::string>& headers, const std::vector<std::vector<std::vector<std::vector<std::string>>>>& ventries, std::vector<uint>& sizes) const;
    [[nodiscard]] std::string getPlainTextTable(const std::vector<std::string>& headers, const std::vector<std::vector<std::string>>& ventries, const std::vector<uint>& sizes) const;
    [[nodiscard]] std::string getJSONTable(const std::vector<std::string>& headers, const std::vector<std::vector<std::string>>& ventries) const;
    void writeText(filetype ftype = txt_ft) const;
    void writeWebApp() const;
    void writeAPool() const;

    // sorting structures
    struct sortWithinCoordGroup {
      sortWithinCoordGroup(const std::vector<ChullPoint>& vpoints, bool sort_energy_ascending = true) : m_points(vpoints), m_sort_energy_ascending(sort_energy_ascending) {};
      const std::vector<ChullPoint>& m_points;
      bool m_sort_energy_ascending; // good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
      bool operator()(uint ci, uint cj);
    };
    struct sortCHullPoints {
      sortCHullPoints(const std::vector<ChullPoint>& vpoints, bool sort_stoich_ascending = true, bool sort_energy_ascending = true) :
          m_points(vpoints), m_sort_stoich_ascending(sort_stoich_ascending), m_sort_energy_ascending(sort_energy_ascending) {};
      const std::vector<ChullPoint>& m_points;
      bool m_sort_stoich_ascending; // good for sorting points in facets (lower vs. upper hemisphere), if facet in lower hemisphere, then sort ascending order (CLOCKWISE, graphing)
      bool m_sort_energy_ascending; // good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
      bool operator()(uint i, uint j) const;
    };
    struct sortFacetsByPoints {
      sortFacetsByPoints(const std::vector<ChullPoint>& vpoints, bool auto_sort_stoich = true, bool sort_stoich_ascending = true, bool auto_sort_energy = true, bool sort_energy_ascending = true) :
          m_points(vpoints), m_auto_sort_stoich(auto_sort_stoich), m_sort_stoich_ascending(sort_stoich_ascending), m_auto_sort_energy(auto_sort_energy), m_sort_energy_ascending(sort_energy_ascending) {};
      const std::vector<ChullPoint>& m_points;
      bool m_auto_sort_stoich;
      bool m_sort_stoich_ascending; // good for sorting points in facets (lower vs. upper hemisphere), if facet in lower hemisphere, then sort ascending order (CLOCKWISE, graphing)
      bool m_auto_sort_energy;
      bool m_sort_energy_ascending; // good for sorting points in facets (lower vs. upper hulls), if lower hull, then sort ascending (ground state is always first)
      bool operator()(const ChullFacet& fi, const ChullFacet& fj) const;
    };
  };
} // namespace chull

#endif // _AFLOW_CHULL_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
