
#ifndef AFLOW_SURFACE_H
#define AFLOW_SURFACE_H

#include <iosfwd>
#include <iostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"
#include "flow/aflow_xclasses.h"
#include "structure/aflow_xstructure.h"

namespace surface {
  double PointInTriangleContribution(const aurostd::xvector<double>& point, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
  double PointInRhombusContribution(const aurostd::xvector<double>& point, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3, const aurostd::xvector<double>& v4);
  double TriangleArea(const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
  bool PlaneGetABCD(double& a, double& b, double& c, double& d, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
  double PlaneDistance(const aurostd::xvector<double>& r, const double& a, const double& b, const double& c, const double& d);
  double PlaneDistance(const aurostd::xvector<double>& r, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
  aurostd::xvector<double> PlaneGetProjection(const aurostd::xvector<double>& r, const double& a, const double& b, const double& c, const double& d);
  aurostd::xvector<double> PlaneGetProjection(const aurostd::xvector<double>& r, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
  aurostd::xvector<double> PlaneGetHKL(
      const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3, const aurostd::xvector<double>& a1, const aurostd::xvector<double>& a2, const aurostd::xvector<double>& a3);
  bool PlaneGetVVV(const aurostd::xvector<double>& hkl,
                   double& area,
                   aurostd::xvector<double>& v1,
                   aurostd::xvector<double>& v2,
                   aurostd::xvector<double>& v3,
                   aurostd::xvector<double>& v4,
                   const aurostd::xvector<double>& a1,
                   const aurostd::xvector<double>& a2,
                   const aurostd::xvector<double>& a3);
  bool PlaneGetVVV_V2(const aurostd::xvector<double>& hkl,
                      double& area,
                      aurostd::xvector<double>& v1,
                      aurostd::xvector<double>& v2,
                      aurostd::xvector<double>& v3,
                      aurostd::xvector<double>& v4,
                      const aurostd::xvector<double>& a1,
                      const aurostd::xvector<double>& a2,
                      const aurostd::xvector<double>& a3);
  double GetPlaneDensityAtoms(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness, const int& type_at);
  double GetPlaneDensityAtoms(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness);
  double GetPlaneDensityBBonds(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness, const double& bbdistance, const int& type_at1, const int& type_at2);
  double GetPlaneDensityBBonds(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness, const double& bbdistance);
  double GetNNeighbors(const xstructure& _str, const int& type_at1, const int& type_at2);
  double GetNNeighbors(const xstructure& _str);
  std::string PrintHKLSigma(int num_types, int num_types_combinations);
  std::string PrintHKLSigmaBB(int num_types, int num_types_combinations, const double& bbfrac, const double& bbdistance, const aurostd::xmatrix<double>& bbdistances);
  bool GetSurfaceHKL(const xstructure& _str, _aflags& aflags, const aurostd::xvector<double>& _hkl, std::vector<std::vector<double>>& planesreducible, std::vector<std::vector<double>>& planesirreducible, std::ostream& oss);
  bool GetSurfaceHKLSearch(const xstructure& _str,
                           _aflags& aflags,
                           const aurostd::xvector<double>& iparams,
                           std::vector<std::vector<double>>& planesreducible,
                           std::vector<std::vector<double>>& planesirreducible,
                           std::vector<std::vector<uint>>& planesirreducible_images,
                           std::ostream& oss,
                           const std::string& smode);
} // namespace surface

namespace slab { // ROMAN CHEPULSKYY
  xstructure MAKE_SLAB(std::string options, std::istream& cin);
  xstructure MAKE_SLAB(std::string options, xstructure& str_in);
  void POSCAR_reading(std::istream& cin);
  double VectorAbsValue(int Layer, int NinLayer, const aurostd::xmatrix<double>& UnitCellVector, const std::vector<std::vector<std::vector<int>>>& LayerSitesDirCoords);
  double VectorScalarMult(int Layer1, int NinLayer1, int Layer2, int NinLayer2, const aurostd::xmatrix<double>& UnitCellVector, const std::vector<std::vector<std::vector<int>>>& LayerSitesDirCoords);
  double CosAngle(int Layer1, int NinLayer1, int Layer2, int NinLayer2, const aurostd::xmatrix<double>& UnitCellVector, const std::vector<std::vector<std::vector<int>>>& LayerSitesDirCoords);
  double hkl_CartCoord_Length(aurostd::xvector<double>& hkl_CartCoord, const aurostd::xvector<double>& hkl);
} // namespace slab

namespace slab { // CO20190601
  aurostd::xvector<double> HKLPlane2Normal(const xstructure& xstr_in, int h, int k, int l);  // CO20190321
  aurostd::xvector<double> HKLPlane2Normal(const aurostd::xmatrix<double>& lattice, int h, int k, int l); // CO20190321
  aurostd::xvector<double> HKLPlane2Normal(const xstructure& xstr_in, const aurostd::xvector<int>& hkl);  // CO20190321
  aurostd::xvector<double> HKLPlane2Normal(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl); // CO20190321
  bool Normal2HKLPlane(const xstructure& xstr_in, const aurostd::xvector<double>& n, aurostd::xvector<int>& hkl);  // CO20190321
  bool Normal2HKLPlane(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<double>& n, aurostd::xvector<int>& hkl); // CO20190321
  std::vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const xstructure& xstr_in, int h, int k, int l); // CO20190321
  std::vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const aurostd::xmatrix<double>& lattice, int h, int k, int l);  // CO20190321
  std::vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const xstructure& xstr_in, const aurostd::xvector<int>& hkl); // CO20190321
  std::vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl);  // CO20190321
  double getSpacingHKLPlane(const xstructure& xstr_in, int h, int k, int l);  // CO20190321
  double getSpacingHKLPlane(const aurostd::xmatrix<double>& lattice, int h, int k, int l); // CO20190321
  double getSpacingHKLPlane(const xstructure& xstr_in, const aurostd::xvector<int>& hkl);  // CO20190321
  double getSpacingHKLPlane(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl); // CO20190321
  double getAngleHKLPlanes(const xstructure& xstr_in, int h1, int k1, int l1, int h2, int k2, int l2);  // CO20190321
  double getAngleHKLPlanes(const aurostd::xmatrix<double>& lattice, int h1, int k1, int l1, int h2, int k2, int l2); // CO20190321
  double getAngleHKLPlanes(const xstructure& xstr_in, const aurostd::xvector<int>& hkl1, const aurostd::xvector<int>& hkl2);  // CO20190321
  double getAngleHKLPlanes(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl1, const aurostd::xvector<int>& hkl2); // CO20190321
  void BringInBoundary(aurostd::xvector<double>& vec, double padding = 0.0);
  //[CO20190520 - plugged into BringInBoundary() with no padding]void Bring2OppositeBoundary(xvector<double>& vec);
  aurostd::xvector<double> getNextAtomInPath(
      const xstructure& xstr_in, const aurostd::xvector<double>& l_cpos, const aurostd::xvector<double>& cpos_starting, std::vector<uint>& atoms2skip, uint& loop_iteration, bool outside_current_cell = false); // CO20190321
  double getDistanceBetweenImages(const xstructure& xstr_in, const aurostd::xvector<double>& n, bool outside_cell = false); // CO20190321
  bool distanceBetweenImages_HKL(const xstructure& xstr_in, const aurostd::xvector<double>& n, double& distance_between_images, bool outside_cell = false); // CO20190321
  bool distanceBetweenImages_Tracing(const xstructure& xstr_in, const aurostd::xvector<double>& n, double& distance_between_images, bool outside_cell = false); // CO20190321

  aurostd::xmatrix<double> getSlabLattice(std::istream& input, const aurostd::xvector<int>& hkl, aurostd::xmatrix<double>& lattice_slab_origbasis, double ang_dev = DEFAULT_V3_ANGLE_DEVIATION, double vlen_max_strict = AUROSTD_MAX_DOUBLE); // CO20190321
  aurostd::xmatrix<double> getSlabLattice(
      const xstructure& _xstr_in, const aurostd::xvector<int>& hkl, aurostd::xmatrix<double>& lattice_slab_origbasis, double ang_dev = DEFAULT_V3_ANGLE_DEVIATION, double vlen_max_strict = AUROSTD_MAX_DOUBLE); // CO20190321
  // easy inputs
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, double v3len_max_strict = AUROSTD_MAX_DOUBLE, std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(
      const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, const _aflags& aflags, double v3len_max_strict = AUROSTD_MAX_DOUBLE, std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(
      const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, std::ofstream& FileMESSAGE, double v3len_max_strict = AUROSTD_MAX_DOUBLE, std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(
      const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, const _aflags& aflags, std::ofstream& FileMESSAGE, double v3len_max_strict = AUROSTD_MAX_DOUBLE, std::ostream& oss = std::cout); // CO20190321
  // load from xoptions
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       std::istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       std::istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       std::istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       std::istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  // input directly
  xstructure CreateSlab_SurfaceLattice(std::istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(std::istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       std::ofstream& FileMESSAGE,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(std::istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(std::istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       std::vector<int>& sc2pcMap_slab,
                                       std::vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict = AUROSTD_MAX_DOUBLE,
                                       std::ostream& oss = std::cout); // CO20190321

} // namespace slab

#endif // AFLOW_SURFACE_H
