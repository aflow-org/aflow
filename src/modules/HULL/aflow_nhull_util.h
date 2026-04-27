// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Nicholas Anderson
// nicholas.anderson@duke.edu

#ifndef AFLOW_NHULL_UTIL_H
#define AFLOW_NHULL_UTIL_H

#include <cstddef>
#include <fstream>
#include <map>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_nhull_entry.h"
#include "aflow_nhull_facet.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "extern/QHULL/libqhullcpp/QhullFacet.h"
#include "flow/aflow_xclasses.h"

//define all constants used in nhull
constexpr double ZERO_TOL = 1e-8;
constexpr double QHULL_POINT_TOL = 1e-8;
constexpr double NHULL_HULL_CHECK_TOL = 1e-10;
constexpr double NHULL_HULL_POINT_CUTOFF = 1e-6; //if dhull less than, entry is a hull point
constexpr double QHULL_FACET_NORM_CUTOFF = 1e-2;
constexpr double NHULL_UNIT_TEST_TOL = 1e-4;
constexpr double NHULL_NORMAL_VECTOR_TOL = 1e-6;
const double ENERGY_TOL = 0.015; //tolerance used to define identical entries
const double ENERGY_CORRECTION_CUTOFF = 1e-4; //cut-off energy (meV) for database enthalpy corrections
const std::vector<std::string> ALLOWED_DFT_TYPES = {"PAW_PBE", "PAW_PBE_KIN"};

namespace nhull {

  // math/general utitilies
  double IQRCutoff(std::vector<Entry>& loaded_entries);
  std::map<std::string, int> reducedStoich(std::map<std::string, int> compound_map);
  std::vector<double> stoichCoord(const std::vector<double>& input_coords);
  std::vector<int> fullIntCoords(const aurostd::xvector<double>& stoich_coord);
  std::vector<double> fullCoords(const std::vector<double>& qhull_coord, const bool has_enthalpy);
  std::vector<std::vector<double>> compoundXvectorConversion(const std::vector<aurostd::xvector<double>>& v);
  int dimensionOfQhullPoint(aurostd::xvector<double> qhull_point);

  //QHULL specific
  aurostd::xvector<double> getQhullFacetNormal(const orgQhull::QhullFacet& current_facet);
  double getQhullFacetOffset(const orgQhull::QhullFacet& current_facet);
  std::vector<double> flatten_vector(const std::vector<std::vector<double>>& v);
  std::vector<double> flatten_vector(const std::vector<aurostd::xvector<double>>& v);

  //output
  void printPoints(std::vector<aurostd::xvector<double>> qhull_many_points);
  std::string specieVecToString(std::vector<std::string> specie_vector);
  std::string reducedCompoundString(const std::string& compound);

  //Thermodynamic post-processing
  double getNormEnergy(const aurostd::xvector<double>& normal);
  double disttoHull(double norm_energy, double distance);
  std::vector<double> decompRatios(int dimension, const std::vector<aurostd::xvector<double>>& vertices, aurostd::xvector<double> qhull_point);
  std::tuple<double, std::map<Entry, double>> nearestDistDecomp(const std::vector<Entry>& entries_on_hull, aurostd::xvector<double> qhull_point, const std::vector<NhullFacet>& facet_list);
  bool checkPointBelowHull(const NhullFacet& facet, const aurostd::xvector<double>& qhull_point);

  std::vector<std::pair<std::string, unsigned long long int>> loadPocc(const std::string& root_path);

  //IQR filtering
  void IQRFilter(std::vector<Entry>& entries, bool pocc_run, bool m_half_hull, const int dimension);
  void calculateOutlierThreshold(const aurostd::xvector<double>& energies, double& upper_threshold, double& lower_threshold);
  std::vector<uint> getOutliers(std::vector<Entry>& entries, const std::vector<size_t>& allowed_entry_indices, const aurostd::xvector<int>& elements_present, bool m_half_hull);
  [[nodiscard]] std::vector<uint> calculateOutliers(std::vector<Entry>& entries, const std::vector<uint>& points_to_consider, bool m_half_hull);

  //formation_enthalpy getter
  double H_f_atom(const aflowlib::_aflowlib_entry& entry, bool neglect_cce, bool& cce_used);

  struct _aflowlib_entry_LIB2sorting : public xStream { // for sorting LIB2 unaries
    _aflowlib_entry_LIB2sorting(aflowlib::_aflowlib_entry& entry, uint index, std::vector<std::string>& specie_vector, std::ofstream& FileMESSAGE, std::ostream& oss);
    // List of relevant aflowlib_entry member vars that we need for sorting:
    std::string m_auid;
    std::string m_aurl;
    std::string m_catalog;
    std::string m_prototype;
    uint m_nspecies;
    uint m_index;
    std::vector<std::string> m_velements_nhull;
    std::vector<std::string> m_species_AURL;
    bool operator<(const _aflowlib_entry_LIB2sorting& other) const;
  };

} // namespace nhull

#endif // AFLOW_NHULL_UTIL_H
