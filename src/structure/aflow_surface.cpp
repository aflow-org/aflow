// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - Sept/Oct/Nov 2007, fixes May 2013
// Corey Oses contributed new slab functionality (see CreateSlab_RigidRotation(), CreateSlab_SurfaceLattice(), and
// related functions) - May/June 2019

#include "structure/aflow_surface.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xtensor.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/SYM/aflow_symmetry.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using std::cerr;
using std::deque;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::string;
using std::stringstream;
using std::vector;

namespace surface {

  constexpr double _eps_ = 0.005;
  constexpr double BBFRAC = 1.3;
  constexpr double RRFRAC = 1.3;
  constexpr int HKLDEF = 4;
  constexpr int oss_short_precision_aflow_surface = 6;

  /// @brief Returns how much angle of a point is in the given triangle measured by an arbitrary circle. Points on
  /// vertices contribute the angle of the vertex, edges contribute 180 deg or pi, inside contributes 360 deg or 2pi.
  /// @param point The point to consider
  /// @param v1 first triangle vertex
  /// @param v2 second triangle vertex
  /// @param v3 third triangle vertex
  /// @return How much of the point is in the triangle (as if the point is a circle snapped to the edges and vertices)
  /// @authors
  /// @mod{ST,20250514,added doxy and cleanup}
  double PointInTriangleContribution(const xvector<double>& point, const xvector<double>& v1, const xvector<double>& v2, const xvector<double>& v3) {
    constexpr double eps = 1.1 * _eps_; // relax a little bit
    const aurostd::xvector<double> proj_point = surface::PlaneGetProjection(point, v1, v2, v3);
    if (const double d = distance(proj_point, point); d > eps) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "point too far from triangle plane = " + aurostd::utype2string(d), _INPUT_ILLEGAL_);
    } // CO20200624

    // check if triangle vertices are too close to each other
    vector<string> problems;
    if (distance(v1, v2) < eps) {
      problems.emplace_back("v1-v2 distance too close");
    }
    if (distance(v2, v3) < eps) {
      problems.emplace_back("v2-v3 distance too close");
    }
    if (distance(v3, v1) < eps) {
      problems.emplace_back("v3-v1 distance too close");
    }
    if (!problems.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aurostd::joinWDelimiter(problems, ',', _INPUT_ILLEGAL_));
    }

    // is in vertices? then contribution is the angle of that vertex
    if (distance(proj_point, v1) < eps) {
      return angle(v1, v2, v3);
    }
    if (distance(proj_point, v2) < eps) {
      return angle(v2, v3, v1);
    }
    if (distance(proj_point, v3) < eps) {
      return angle(v3, v1, v2);
    }

    // is in edges? then contribution is pi if the point is actually between the vertices of the edge, 0 if not
    if (aurostd::abs(sin(v1 - proj_point, v2 - proj_point)) < eps) {
      return angle(proj_point, v1, v2) > eps ? pi : 0.0;
    }
    if (aurostd::abs(sin(v2 - proj_point, v3 - proj_point)) < eps) {
      return angle(proj_point, v2, v3) > eps ? pi : 0.0;
    }
    if (aurostd::abs(sin(v3 - proj_point, v1 - proj_point)) < eps) {
      return angle(proj_point, v3, v1) > eps ? pi : 0.0;
    }

    // is inside, then the angles are 360... otherwise less !
    if (aurostd::abs(angle(proj_point, v1, v2) + angle(proj_point, v2, v3) + angle(proj_point, v3, v1) - 2 * pi) < eps) {
      return 2.0 * pi; // inside, return 2pi
    }
    return 0.0; // not inside
  }

  /// @brief Returns how much angle of a point is in the given rhombus measured by an arbitrary circle. Same rules as
  /// @c PointInTriangleContribution. Delegates to @c PointInTriangleContribution.
  /// @param point The point to consider
  /// @param v1 of rhombus vertex
  /// @param v2 of rhombus vertex
  /// @param v3 of rhombus vertex
  /// @param v4 of rhombus vertex
  /// @return How much of the point is in the rhombus (as if the point is a circle snapped to the edges and vertices)
  /// @authors
  /// @mod{ST,20250514,added doxy and cleanup}
  double PointInRhombusContribution(const aurostd::xvector<double>& point, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3, const aurostd::xvector<double>& v4) {
    double out = 0.0;
    out += surface::PointInTriangleContribution(point, v2, v3, v4); // double count because pick twice area in rhombi
    out += surface::PointInTriangleContribution(point, v3, v4, v1); // double count because pick twice area in rhombi
    out += surface::PointInTriangleContribution(point, v4, v1, v2); // double count because pick twice area in rhombi
    out += surface::PointInTriangleContribution(point, v1, v2, v3); // double count because pick twice area in rhombi
    return out / 2.0;
  }

  /// @brief Calculate the area of the triangle bounded by the 3 points in 3D space.
  /// @param v1 First point as a 3-vec
  /// @param v2 Second point as a 3-vec
  /// @param v3 Third point as a 3-vec
  /// @return The area of the triangle
  /// @authors
  /// @mod{ST,20241205,added doxy\, clean}
  ///
  /// see "https://mathworld.wolfram.com/TriangleArea.html"
  double TriangleArea(const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3) {
    // area of ABC triangle = 1/2 * sqrt( |AB|^2 * |AC|^2 - (AB dot AC)^2 )
    const aurostd::xvector<double> v12 = v2 - v1;
    const aurostd::xvector<double> v31 = v1 - v3;
    const double v12dotv31 = scalar_product(v12, v31);
    return 0.5 * sqrt(scalar_product(v12, v12) * scalar_product(v31, v31) - v12dotv31 * v12dotv31);
  }

  /// @brief Calculate the general equation of the plane defined by three points in 3D space.
  /// @param a[out] Coefficient a of the general equation
  /// @param b[out] Coefficient b of the general equation
  /// @param c[out] Coefficient c of the general equation
  /// @param d[out] Constant d of the general equation
  /// @param v1[in] First point defining the plane
  /// @param v2[in] Second point defining the plane
  /// @param v3[in] Third point defining the plane
  /// @return true
  /// @authors
  /// @mod{ST,20241205,added doxy\, clean}
  ///
  /// The general equation of the plane is of the form \f$ ax + by + cz + d = 0 \f$
  /// where \f$ d = -ax_0 - by_0 - cz_0 \f$ with a nonzero normal vector \f$ n = (a, b, c) \f$
  /// through the point \f$ (x_0, y_0, z_0) \f$.
  ///
  /// @note Because there are infinitely many solutions, the normal vector defined by a,b,c is not a unit vector.
  /// In this implementation, its magnitude is the area of the parallelogram formed by the three given points.
  /// You will often want to normalize the vector before/during use.
  ///
  /// see "https://mathworld.wolfram.com/Plane.html"
  bool PlaneGetABCD(double& a, double& b, double& c, double& d, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3) {
    const aurostd::xvector<double> v21 = v2 - v1;
    const aurostd::xvector<double> v31 = v3 - v1;

    a = v21(2) * v31(3) - v21(3) * v31(2);
    b = v21(3) * v31(1) - v21(1) * v31(3);
    c = v21(1) * v31(2) - v21(2) * v31(1);
    d = -a * v1(1) - b * v1(2) - c * v1(3);
    return true;
  }

  /// @brief Gives the absolute valued distance of a point from a plane defined by its general equation.
  /// @param r The point as a 3-vec
  /// @param a Coefficient a of the general equation
  /// @param b Coefficient b of the general equation
  /// @param c Coefficient c of the general equation
  /// @param d Constant d of the general equation
  /// @return Absolute valued distance of point r from the plane
  /// @authors
  /// @mod{ST,20241205,added doxy}
  ///
  /// The general equation of the plane is of the form \f$ ax + by + cz + d = 0 \f$
  ///
  /// see "https://mathworld.wolfram.com/Point-PlaneDistance.html"
  double PlaneDistance(const aurostd::xvector<double>& r, const double& a, const double& b, const double& c, const double& d) {
    return aurostd::abs(a * r(1) + b * r(2) + c * r(3) + d) / sqrt(a * a + b * b + c * c);
  }

  /// @brief Gives the absolute valued distance of a point from a plane defined by three points.
  /// @param r The point distant from the plane
  /// @param v1 First point defining the plane
  /// @param v2 Second point defining the plane
  /// @param v3 Third point defining the plane
  /// @return Absolute valued distance of point r from the plane
  /// @authors
  /// @mod{ST,20241205,added doxy}
  double PlaneDistance(const aurostd::xvector<double>& r, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3) {
    double a;
    double b;
    double c;
    double d;
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    return surface::PlaneDistance(r, a, b, c, d);
  }

  /// @brief Get the projection of a given distant point onto a plane defined by its general equation
  /// @param r The point distant from the plane
  /// @param a Coefficient a of the general equation
  /// @param b Coefficient b of the general equation
  /// @param c Coefficient c of the general equation
  /// @param d Constant d of the general equation
  /// @return The projected point from the projection of the distant point onto the plane
  /// @authors
  /// @mod{ST,20241205,added doxy\, clean}
  ///
  /// @note The returned vector is from the plane to the point. ( | --> . )
  ///
  /// Calculation: projected = given - dist * unit_normal
  aurostd::xvector<double> PlaneGetProjection(const aurostd::xvector<double>& r, const double& a, const double& b, const double& c, const double& d) {
    const aurostd::xvector<double> rorth{a, b, c};
    const double dist = surface::PlaneDistance(r, a, b, c, d);
    const aurostd::xvector<double> rproj = r - dist * rorth / sqrt(a * a + b * b + c * c);
    return rproj;
  }

  /// @brief Get the projection of a given distant point onto a plane defined by three points
  /// @param r The point distant from the plane
  /// @param v1 First point defining the plane
  /// @param v2 Second point defining the plane
  /// @param v3 Third point defining the plane
  /// @return The projected point from the projection of the distant point onto the plane
  /// @authors
  /// @mod{ST,20241205,added doxy}
  ///
  /// Calculation: projected = given - dist * unit_normal
  aurostd::xvector<double> PlaneGetProjection(const aurostd::xvector<double>& r, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3) {
    double a;
    double b;
    double c;
    double d;
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    return surface::PlaneGetProjection(r, a, b, c, d);
  }

  /// @brief Get the hkl indices for a given plane intersecting the lattice vectors of a given lattice
  /// @param v1 First point defining the plane
  /// @param v2 Second point defining the plane
  /// @param v3 Third point defining the plane
  /// @param a1 First row of the cartesian lattice matrix
  /// @param a2 Second row of the cartesian lattice matrix
  /// @param a3 Third row of the cartesian lattice matrix
  /// @return The hkl indices
  /// @authors
  /// @mod{ST,20241205,added doxy\, clean}
  aurostd::xvector<double> PlaneGetHKL(
      const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3, const aurostd::xvector<double>& a1, const aurostd::xvector<double>& a2, const aurostd::xvector<double>& a3) {
    const aurostd::xvector<double> hkl(3);
    constexpr double eps = _eps_;
    double a;
    double b;
    double c;
    double d;
    surface::PlaneGetABCD(a, b, c, d, v1, v2, v3);
    // hkl index = xyz component of ( plane normal dot lattice vector) / d
    hkl[1] = -(a * a1(1) + b * a1(2) + c * a1(3)) / d;
    hkl[2] = -(a * a2(1) + b * a2(2) + c * a2(3)) / d;
    hkl[3] = -(a * a3(1) + b * a3(2) + c * a3(3)) / d;
    if (std::abs(hkl[1]) < eps) {
      hkl[1] = 0.0;
    }
    if (std::abs(hkl[2]) < eps) {
      hkl[2] = 0.0;
    }
    if (std::abs(hkl[3]) < eps) {
      hkl[3] = 0.0;
    }
    return hkl;
  }

  /// @brief Gets the vertices and area of the intersection of hkl plane with lattice vectors
  /// @param hkl[in] the hkl indices of a lattice plane
  /// @param area[out] the area of the plane within the lattice
  /// @param v1[out] First intersection of hkl plane with lattice vectors
  /// @param v2[out] Second intersection of hkl plane with lattice vectors
  /// @param v3[out] Third intersection of hkl plane with lattice vectors
  /// @param v4[out] Fourth intersection of hkl plane with lattice vectors if exists
  /// @param a1[in] First row of the cartesian lattice matrix
  /// @param a2[in] Second row of the cartesian lattice matrix
  /// @param a3[in] Third row of the cartesian lattice matrix
  /// @return true if rhombus
  /// @authors
  /// @mod{ST,20241205,added doxy\, clean}
  bool PlaneGetVVV(const aurostd::xvector<double>& hkl,
                   double& area,
                   aurostd::xvector<double>& v1,
                   aurostd::xvector<double>& v2,
                   aurostd::xvector<double>& v3,
                   aurostd::xvector<double>& v4,
                   const aurostd::xvector<double>& a1,
                   const aurostd::xvector<double>& a2,
                   const aurostd::xvector<double>& a3) {
    const bool LDEBUG = XHOST.DEBUG;
    bool isrhombus = true;
    const double h = hkl(1);
    const double k = hkl(2);
    const double l = hkl(3);
    constexpr double eps = _eps_;
    v1.clear();
    v2.clear();
    v3.clear();
    v4.clear();

    // *****************************************
    // get the 4 vertices of rhombus
    if (std::abs(h) > eps && std::abs(k) > eps && std::abs(l) > eps) { // XYZ axis DEFINED
      if (LDEBUG) {
        cerr << "XYZ axis DEFINED" << endl;
      }
      const aurostd::xvector<double> a1h = a1 * (1.0 / h);
      const aurostd::xvector<double> a2k = a2 * (1.0 / k);
      const aurostd::xvector<double> a3l = a3 * (1.0 / l);
      v1 = -a1h + a2k + a3l;
      v2 = +a1h - a2k + a3l;
      v3 = +a1h + a2k - a3l;
      isrhombus = false;
      area = surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) < eps && std::abs(k) > eps && std::abs(l) > eps) { // X axis INFINITE
      if (LDEBUG) {
        cerr << "X axis INFINITE - 0kl" << endl;
      }
      v1 = a2 * (1.0 / k);
      v2 = a3 * (1.0 / l);
      v3 = v1 + a1;
      v4 = v2 + a1;
      if (hkl != surface::PlaneGetHKL(v1, v2, v3, a1, a2, a3)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "hkl problem in \"[X] axis INFINITE\"", _INPUT_ILLEGAL_);
      }
      isrhombus = true;
      area = 2.0 * surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) > eps && std::abs(k) < eps && std::abs(l) > eps) { // Y axis INFINITE
      if (LDEBUG) {
        cerr << "Y axis INFINITE - h0l" << endl;
      }
      v1 = a1 * (1.0 / h);
      v2 = a3 * (1.0 / l);
      v3 = v1 + a2;
      v4 = v2 + a2;
      if (hkl != surface::PlaneGetHKL(v1, v2, v3, a1, a2, a3)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "hkl problem in \"[Y] axis INFINITE\"", _INPUT_ILLEGAL_);
      } // CO20200624
      isrhombus = true;
      area = 2.0 * surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) > eps && std::abs(k) > eps && std::abs(l) < eps) { // Z axis INFINITE
      if (LDEBUG) {
        cerr << "Z axis INFINITE - hk0" << endl;
      }
      v1 = a1 * (1.0 / h);
      v2 = a2 * (1.0 / k);
      v3 = v1 + a3;
      v4 = v2 + a3;
      if (hkl != surface::PlaneGetHKL(v1, v2, v3, a1, a2, a3)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "hkl problem in \"[Z] axis INFINITE\"", _INPUT_ILLEGAL_);
      } // CO20200624
      isrhombus = true;
      area = 2.0 * surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) > eps && std::abs(k) < eps && std::abs(l) < eps) { // YZ axis INFINITE
      if (LDEBUG) {
        cerr << "YZ axis INFINITE - h00" << endl;
      }
      v1 = a1 * (1.0 / h);
      v2 = v1 + a2;
      v3 = v1 + a3;
      v4 = v1 + a2 + a3;
      isrhombus = true;
      area = 2.0 * surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) < eps && std::abs(k) > eps && std::abs(l) < eps) { // XZ axis INFINITE
      if (LDEBUG) {
        cerr << "XZ axis INFINITE - 0k0" << endl;
      }
      v2 = a2 * (1.0 / k);
      v1 = v2 + a1;
      v3 = v2 + a3;
      v4 = v2 + a1 + a3;
      isrhombus = true;
      area = 2.0 * surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) < eps && std::abs(k) < eps && std::abs(l) > eps) { // XY axis INFINITE
      if (LDEBUG) {
        cerr << "XY axis INFINITE - 00l" << endl;
      }
      v3 = a3 * (1.0 / l);
      v1 = v3 + a1;
      v2 = v3 + a2;
      v4 = v3 + a1 + a2;
      isrhombus = true;
      area = 2.0 * surface::TriangleArea(v1, v2, v3);
    }
    if (std::abs(h) < eps && std::abs(k) < eps && std::abs(l) < eps) { // XYZ axis INFINITE
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "h,k,l cannot be 0 0 0", _INPUT_ILLEGAL_);
    } // CO20200624
    for (int i = 1; i <= 3; i++) {
      if (std::abs(v1[i]) < eps) {
        v1[i] = 0.0;
      }
      if (std::abs(v2[i]) < eps) {
        v2[i] = 0.0;
      }
      if (std::abs(v3[i]) < eps) {
        v3[i] = 0.0;
      }
      if (std::abs(v4[i]) < eps) {
        v4[i] = 0.0;
      }
    }
    return isrhombus;
  }

  /// @brief Gets the vertices and area of the intersection of hkl plane with lattice vectors
  /// @param hkl[in] the hkl indices of a lattice plane
  /// @param area[out] the area of the plane within the lattice
  /// @param v1[out] First intersection of hkl plane with lattice vectors
  /// @param v2[out] Second intersection of hkl plane with lattice vectors
  /// @param v3[out] Third intersection of hkl plane with lattice vectors
  /// @param v4[out] Fourth intersection of hkl plane with lattice vectors if exists
  /// @param a1[in] First row of the cartesian lattice matrix
  /// @param a2[in] Second row of the cartesian lattice matrix
  /// @param a3[in] Third row of the cartesian lattice matrix
  /// @return true if the intersection is a rhombus, false if it is a triangle
  ///
  /// @note This version condenses the math. It should be faster and mathematically equivalent.
  /// @authors
  /// @mod{ST,20241205,added doxy\, clean}
  bool PlaneGetVVV_V2(const aurostd::xvector<double>& hkl,
                      double& area,
                      aurostd::xvector<double>& v1,
                      aurostd::xvector<double>& v2,
                      aurostd::xvector<double>& v3,
                      aurostd::xvector<double>& v4,
                      const aurostd::xvector<double>& a1,
                      const aurostd::xvector<double>& a2,
                      const aurostd::xvector<double>& a3) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    constexpr double eps = _eps_;
    const double h = hkl(1);
    const double k = hkl(2);
    const double l = hkl(3);
    bool isrhombus = true;

    if (std::abs(h) < eps && std::abs(k) < eps && std::abs(l) < eps) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "h,k,l cannot be 0 0 0", _INPUT_ILLEGAL_);
    }

    const aurostd::xvector<double> a1h = std::abs(h) > eps ? a1 * (1.0 / h) : a1;
    const aurostd::xvector<double> a2k = std::abs(k) > eps ? a2 * (1.0 / k) : a2;
    const aurostd::xvector<double> a3l = std::abs(l) > eps ? a3 * (1.0 / l) : a3;

    if (std::abs(h) > eps && std::abs(k) > eps && std::abs(l) > eps) { // XYZ axis DEFINED
      if (LDEBUG) {
        cerr << "XYZ axis DEFINED" << endl;
      }
      v1 = -a1h + a2k + a3l;
      v2 = +a1h - a2k + a3l;
      v3 = +a1h + a2k - a3l;
      v4.clear();
      isrhombus = false;
      area = surface::TriangleArea(v1, v2, v3);
    } else {
      v1 = a1h;
      v2 = a2k;
      v3 = a3l;
      v4 = a2k;

      if (std::abs(h) <= eps) {
        v1 += std::abs(l) > eps ? a3l : a2k;
        v4 += a1h;
      }
      if (std::abs(k) <= eps) {
        v2 += std::abs(h) > eps ? a1h : a3l;
        v4 += a3l;
      }
      if (std::abs(l) <= eps) {
        v3 += std::abs(h) > eps ? a1h : a2k;
        v4 += std::abs(k) > eps ? a3l : a1h;
      }

      isrhombus = true;
      area = 2 * surface::TriangleArea(v1, v2, v3);
    }

    for (int i = 1; i <= 3; i++) {
      if (std::abs(v1[i]) < eps) {
        v1[i] = 0.0;
      }
      if (std::abs(v2[i]) < eps) {
        v2[i] = 0.0;
      }
      if (std::abs(v3[i]) < eps) {
        v3[i] = 0.0;
      }
      if (std::abs(v4[i]) < eps) {
        v4[i] = 0.0;
      }
    }

    return isrhombus;
  }

  /// @brief Gets the density of atoms in the plane
  /// @param _str The atomic structure
  /// @param hkl the hkl indices
  /// @param roughness
  /// @param type_at index for the atom type to consider for density, use -1 to consider all types
  /// @return the density
  /// @authors
  /// @mod{ST,20241205,added doxy\, cleanup}
  double GetPlaneDensityAtoms(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness, const int& type_at) {
    if (type_at > 0 && type_at >= _str.num_each_type.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "type_at greater than size", _INPUT_ILLEGAL_);
    }

    const aurostd::xvector<double> a1 = _str.lattice(1);
    const aurostd::xvector<double> a2 = _str.lattice(2);
    const aurostd::xvector<double> a3 = _str.lattice(3);
    const xstructure str = BringInCell(_str);
    aurostd::xvector<double> v1(3);
    aurostd::xvector<double> v2(3);
    aurostd::xvector<double> v3(3);
    aurostd::xvector<double> v4(3);
    double area;
    bool isrhombus;

    vector<double> atom_densities(_str.num_each_type.size(), 0.0);

    isrhombus = PlaneGetVVV(hkl, area, v1, v2, v3, v4, a1, a2, a3);

    double a;
    double b;
    double c;
    double d;
    PlaneGetABCD(a, b, c, d, v1, v2, v3);

    const double radius = 1.1 * aurostd::max(aurostd::modulus(v1), aurostd::modulus(v2), aurostd::modulus(v3), aurostd::modulus(v4)) + aurostd::max(aurostd::modulus(a1), aurostd::modulus(a2), aurostd::modulus(a3));
    const aurostd::xvector<int> dims = LatticeDimensionSphere(_str.lattice, radius);

    for (int ii = dims(1); ii >= -dims(1); ii--) {
      for (int jj = dims(2); jj >= -dims(2); jj--) {
        for (int kk = dims(3); kk >= -dims(3); kk--) {
          for (size_t iat = 0; iat < str.atoms.size(); iat++) {
            const aurostd::xvector<double> rrr = ii * a1 + jj * a2 + kk * a3 + str.atoms[iat].cpos;
            const double dist = PlaneDistance(rrr, a, b, c, d);
            if (dist < roughness) {
              if (!isrhombus) {
                atom_densities.at(str.atoms[iat].type) += surface::PointInTriangleContribution(rrr, v1, v2, v3); // 1st triangle
              }
              if (isrhombus) {
                atom_densities.at(str.atoms[iat].type) += surface::PointInRhombusContribution(rrr, v1, v2, v3, v4); // rhombus
              }
            }
          }
        }
      }
    }
    double out;
    if (type_at < 0) {
      out = aurostd::sum(atom_densities);
    } else {
      out = atom_densities[type_at];
    }
    return out / (2 * pi * area * str.scale * str.scale);
  }

  double GetPlaneDensityAtoms(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness) {
    return GetPlaneDensityAtoms(_str, hkl, roughness, -1);
  }

  double GetPlaneDensityBBonds(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness, const double& bbdistance, const int& type_at1, const int& type_at2) {
    aurostd::xmatrix<double> lattice(3, 3);
    lattice = (_str.lattice);
    aurostd::xvector<double> a1(3);
    aurostd::xvector<double> a2(3);
    aurostd::xvector<double> a3(3); // lattice vectors
    a1 = lattice(1);
    a2 = lattice(2);
    a3 = lattice(3); // a1,a2,a3 are the rows of the lattice matrix
    aurostd::xvector<double> v1(3);
    aurostd::xvector<double> v2(3);
    aurostd::xvector<double> v3(3);
    aurostd::xvector<double> v4(3);
    double area;
    xstructure str(_str);
    str = BringInCell(str);
    bool isrhombus;
    double a;
    double b;
    double c;
    double d;
    double afound = 0.0;
    const double eps = _eps_;

    if (roughness) {
      ;
    } // phony, just to use

    vector<aurostd::xvector<double>*> grid_far_atoms_cpos;
    vector<aurostd::xvector<double>*> grid_close_atoms_cpos;
    vector<int> grid_far_atoms_type;
    vector<int> grid_close_atoms_type;
    aurostd::xvector<double>* grid_atoms_cpos_ptr;

    isrhombus = PlaneGetVVV(hkl, area, v1, v2, v3, v4, a1, a2, a3);
    PlaneGetABCD(a, b, c, d, v1, v2, v3);
    //  cerr << "a=" << a << " b=" << b << " c=" << c << " d=" << d << endl;

    // *****************************************
    // create search over all atoms.
    // distance http://en.wikipedia.org/wiki/Plane_(mathematics)
    aurostd::xvector<int> dims(3);
    double radius;
    double dist;
    double num;
    double den;
    double u;
    aurostd::xvector<double> rrr(3);
    aurostd::xvector<double> rrr1(3);
    aurostd::xvector<double> rrr2(3);
    radius = 1.1 * max(modulus(v1), modulus(v2), modulus(v3), modulus(v4)) + max(modulus(a1), modulus(a2), modulus(a3));
    dims = LatticeDimensionSphere(lattice, radius);

    for (int ii = dims(1); ii >= -dims(1); ii--) {
      for (int jj = dims(2); jj >= -dims(2); jj--) {
        for (int kk = dims(3); kk >= -dims(3); kk--) {
          for (size_t iat = 0; iat < str.atoms.size(); iat++) {
            rrr = ((double) ii) * a1 + ((double) jj) * a2 + ((double) kk) * a3 + str.atoms[iat].cpos;
            dist = aurostd::abs(a * rrr(1) + b * rrr(2) + c * rrr(3) + d) / sqrt(a * a + b * b + c * c);
            if (dist <= bbdistance) {
              if (dist > 0.01) { // FAR
                grid_atoms_cpos_ptr = new aurostd::xvector<double>(3);
                *grid_atoms_cpos_ptr = rrr;
                grid_far_atoms_cpos.push_back(grid_atoms_cpos_ptr);
                grid_far_atoms_type.push_back(str.atoms[iat].type);
              } else { // CLOSE
                grid_atoms_cpos_ptr = new aurostd::xvector<double>(3);
                *grid_atoms_cpos_ptr = rrr;
                grid_close_atoms_cpos.push_back(grid_atoms_cpos_ptr);
                grid_close_atoms_type.push_back(str.atoms[iat].type);
              }
            }
          }
        }
      }
    }
    //  cout << grid_close_atoms_cpos.size() << " " << grid_far_atoms_cpos.size() << endl;
    int matches = 0;
    afound = 0.0;
    for (size_t iat_close = 0; iat_close < grid_close_atoms_cpos.size(); iat_close++) {
      rrr1 = *grid_close_atoms_cpos[iat_close];
      for (size_t iat_far = 0; iat_far < grid_far_atoms_cpos.size(); iat_far++) {
        rrr2 = *grid_far_atoms_cpos[iat_far];
        if ((type_at1 < 0 && type_at2 < 0) || (grid_close_atoms_type.at(iat_close) == type_at1 && grid_far_atoms_type.at(iat_far) == type_at2) ||
            (grid_close_atoms_type.at(iat_close) == type_at2 && grid_far_atoms_type.at(iat_far) == type_at1)) {
          dist = distance(rrr2, rrr1);
          if (dist <= bbdistance && dist > 0.2) { // they are close... but not the same
            // intersection http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/
            // P = P1 + u (P2 - P1)    rrr = rrr1 + u (rrr2 - rrr1)
            num = a * rrr1(1) + b * rrr1[2] + c * rrr1[3] + d;
            den = a * (rrr1(1) - rrr2(1)) + b * (rrr1[2] - rrr2[2]) + c * (rrr1[3] - rrr2[3]);
            // if(aurostd::abs(den)<eps/10.0 && aurostd::abs(num)>eps/10.0) {throw
            // aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,aurostd::utype2string(num)+"
            // "+aurostd::utype2string(den),_INPUT_ILLEGAL_);}  //CO20200624 if(aurostd::abs(den)<eps/10.0 &&
            // aurostd::abs(num)<=aurostd::abs(den)) {num=0.0;den=1.0;}
            if (aurostd::abs(den) > eps / 10.0) // to avoid rrr1,rrr2,rrr coplanar with the plane
            {                                   // found
              u = num / den;
              if (u < -2 * eps || u > 1.0 + 2 * eps) {
                cerr << "u=" << u << " num=" << num << " den=" << den << endl;
              }
              // if(u<0.0 && u>=-eps) u=0.0;
              // if(u>1.0 && u<=1.0+eps) u=1.0;
              // if(u<-eps || u>1+eps) {throw
              // aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"u="aurostd::utype2string(u),_INPUT_ILLEGAL_);}
              // //CO20200624
              if (u < 0.0) {
                u = 0.0;
              }
              if (u > 1.0) {
                u = 1.0;
              }
              rrr = rrr1 + u * (rrr2 - rrr1);
              //  point=PlaneGetProjection(_point,v1,v2,v3);
              if (PlaneDistance(rrr, a, b, c, d) < eps) {
                if (!isrhombus) {
                  afound += surface::PointInTriangleContribution(rrr, v1, v2, v3);   // 1st triangle
                }
                if (isrhombus) {
                  afound += surface::PointInRhombusContribution(rrr, v1, v2, v3, v4); // rhombus
                }
                matches++;
              }
            }
          }
        }
      }
    }
    // EXIT
    // cerr << matches << " " << afound/(2*pi*area*str.scale*str.scale) << endl;
    for (size_t i = 0; i < grid_far_atoms_cpos.size(); i++) {
      delete grid_far_atoms_cpos[i];
    }
    grid_far_atoms_cpos.clear();
    for (size_t i = 0; i < grid_close_atoms_cpos.size(); i++) {
      delete grid_close_atoms_cpos[i];
    }
    grid_close_atoms_cpos.clear();
    return afound / (2 * pi * area * str.scale * str.scale);
  }

  double GetPlaneDensityBBonds(const xstructure& _str, const aurostd::xvector<double>& hkl, const double& roughness, const double& bbdistance) {
    return surface::GetPlaneDensityBBonds(_str, hkl, roughness, bbdistance, -1, -1);
  }

  /// @brief Calculates the shortest pairwise distance of atoms in a structure
  /// @param _str The input structure
  /// @param type_at1 The element type to consider for first atom, use -1 for any
  /// @param type_at2 The element type to consider for second atom, use -1 for any
  /// @return The shortest distance between the atom types
  /// @authors
  /// @mod{ST,20241217,added doxy\, clean declarations}
  double GetNNeighbors(const xstructure& _str, const int& type_at1, const int& type_at2) {
    const xstructure str = ReScale(BringInCell(_str), 1.0);
    const aurostd::xvector<double> a1 = str.lattice(1);
    const aurostd::xvector<double> a2 = str.lattice(2);
    const aurostd::xvector<double> a3 = str.lattice(3);
    const double radius = 1.5 * max(modulus(a1), modulus(a2), modulus(a3));
    const aurostd::xvector<int> dims = LatticeDimensionSphere(str.lattice, radius);
    double nndistance = 10.0 * radius;

    for (size_t iat1 = 0; iat1 < str.atoms.size(); iat1++) {
      for (int i = -dims(1); i <= dims(1); i++) {
        for (int j = -dims(2); j <= dims(2); j++) {
          for (int k = -dims(3); k <= dims(3); k++) {
            for (size_t iat2 = 0; iat2 < str.atoms.size(); iat2++) {
              if ((type_at1 < 0 && type_at2 < 0) || (str.atoms.at(iat1).type == type_at1 && str.atoms.at(iat2).type == type_at2) || (str.atoms.at(iat1).type == type_at2 && str.atoms.at(iat2).type == type_at1)) {
                const aurostd::xvector<double> rrr = i * a1 + j * a2 + k * a3 + str.atoms.at(iat2).cpos;
                const double r = modulus(str.atoms.at(iat1).cpos - rrr);
                if (r <= nndistance && r > _eps_) {
                  nndistance = r;
                }
              }
            }
          }
        }
      }
    }
    return nndistance;
  }

  /// @brief Calculates the shortest pairwise distance of any atoms in a structure
  /// @param _str The input structure
  /// @return The shortest distance between the atoms
  /// @authors
  /// @mod{ST,20241217,added doxy}
  double GetNNeighbors(const xstructure& _str) {
    return surface::GetNNeighbors(_str, -1, -1);
  }

  /// @brief Outputs headers for a hkl info output list
  /// @param num_types number of atom types used to control number of headers
  /// @param num_types_combinations number of atom type combinations to control number of headers
  /// @return lines for the header of hkl output info
  string PrintHKLSigma(const int num_types, const int num_types_combinations) {
    ostringstream aus;
    if (num_types == 1) {
      aus << "     h           k           l        sigma(#/AA)  Nb(#/AA)    " << endl;
    } else {
      // 1st line
      aus << "     h           k           l        ";
      for (int j = 0; j <= num_types; j++) {
        aus << "sigma(#/AA) ";
      }
      for (int j = 0; j <= num_types_combinations; j++) {
        aus << " Nb(#/AA)   ";
      }
      aus << endl;
      // 2nd line
      aus << "                                      ";
      aus << "   T=(*)    ";
      for (int j = 0; j < num_types; j++) {
        aus << "   T=(" << j << ")    ";
      }
      aus << " TT=(*-*)   ";
      for (int it1 = 0; it1 < num_types; it1++) {
        for (int it2 = it1; it2 < num_types; it2++) {
          aus << " TT=(" << it1 << "-" << it2 << ")   ";
        }
      }
      aus << endl;
      // 3rd line
    }
    return aus.str();
  }

  string PrintHKLSigmaBB(const int num_types, const int num_types_combinations, const double& /*bbfrac*/, const double& bbdistance, const aurostd::xmatrix<double>& bbdistances) {
    ostringstream aus;
    aus.setf(std::ios::fixed, std::ios::floatfield);
    aus.precision(4);
    aus << surface::PrintHKLSigma(num_types, num_types_combinations);
    if (num_types == 1) {
    } else {
      // 3rd line
      aus << "                                       ";
      for (int j = 0; j < num_types + 1; j++) {
        aus << "            ";
      }
      aus << "b=" << bbdistance << "    ";
      for (int it1 = 0; it1 < num_types; it1++) {
        for (int it2 = it1; it2 < num_types; it2++) {
          aus << "b=" << bbdistances(it1, it2) << "    ";
        }
      }
      aus << endl;
    }
    return aus.str();
  }

  string PrintNNdists(const int num_types, const int /*num_types_combinations*/, const double& /*rrdist*/, const double& nndist, const aurostd::xmatrix<double>& nndists) {
    ostringstream aus;
    aus.setf(std::ios::fixed, std::ios::floatfield);
    aus.precision(5);
    aus << "nndist(*-*)=" << nndist << endl;
    if (num_types > 1) {
      for (int it1 = 0; it1 < num_types; it1++) {
        for (int it2 = it1; it2 < num_types; it2++) {
          aus << "nndists(" << it1 << "-" << it2 << ")=" << nndists(it1, it2) << " ";
        }
      }
      aus << endl;
    }
    return aus.str();
  }

  bool AddPlaneHKL(const xstructure& str,
                   const double& h,
                   const double& k,
                   const double& l,
                   const double& roughness,
                   const double& bbdistance,
                   const aurostd::xmatrix<double>& bbdistances,
                   const double& eps,
                   const int& jossmax,
                   vector<double>& plane,
                   vector<vector<double>>& planes,
                   const bool& osswrite1,
                   ostream& oss1,
                   const bool& osswrite2,
                   ostream& oss2) {
    if (std::abs(h) < eps && std::abs(k) < eps && std::abs(l) < eps) {
      return false;
    }
    bool plane_added = false;
    const aurostd::xvector<double> hkl = {h, k, l};
    const int num_types = str.num_each_type.size();
    if ((std::abs(h) >= 1.0 || std::abs(h) < eps) && (std::abs(k) >= 1.0 || std::abs(k) < eps) && (std::abs(l) >= 1.0 || std::abs(l) < eps)) {
      bool hkl_found = false;
      for (size_t i = 0; i < planes.size() && !hkl_found; i++) {
        hkl_found = std::abs(h - planes[i].at(0)) < eps && std::abs(k - planes[i].at(1)) < eps && std::abs(l - planes[i].at(2)) < eps;
      }
      if (hkl_found == false) {
        // new operation, generate and save it
        // ---------------------- hkl ----------------------
        plane.at(0) = h;
        plane.at(1) = k;
        plane.at(2) = l;
        // ---------------------- density ------------------
        const double density = surface::GetPlaneDensityAtoms(str, hkl, roughness);
        int i = 3;
        plane.at(i++) = density;
        if (num_types > 1) {
          for (int it1 = 0; it1 < num_types; it1++) {
            plane.at(i++) = surface::GetPlaneDensityAtoms(str, hkl, roughness, it1);
          }
        }
        // ---------------------- bbonds -------------------
        const double bbonds = surface::GetPlaneDensityBBonds(str, hkl, roughness, bbdistance);
        plane.at(i++) = bbonds;
        if (num_types > 1) {
          for (int it1 = 0; it1 < num_types; it1++) {
            for (int it2 = it1; it2 < num_types; it2++) {
              plane.at(i++) = surface::GetPlaneDensityBBonds(str, hkl, roughness, bbdistances(it1, it2), it1, it2);
            }
          }
        }
        if (density > 0.0) {
          planes.push_back(plane); // save them
          plane_added = true;
          //    oss1 << "*";oss1.flush();
        }
      }
      if (plane_added) {
        for (int j = 0; j <= jossmax; j++) {
          if (osswrite1) {
            oss1 << (plane[j] >= 0 ? "  " : " ") << (std::abs(plane[j]) < 10.0 ? " " : "") << plane[j] << " ";
          }
          if (osswrite2) {
            oss2 << (plane[j] >= 0 ? "  " : " ") << (std::abs(plane[j]) < 10.0 ? " " : "") << plane[j] << " ";
          }
        }
        if (osswrite1) {
          oss1 << planes.size() << " ";
        }
        if (osswrite2) {
          oss2 << planes.size() << " ";
        }
        if (osswrite1) {
          oss1 << endl;
        }
        if (osswrite1) {
          oss1.flush();
        }
        if (osswrite2) {
          oss2 << endl;
        }
        if (osswrite2) {
          oss2.flush();
        }
      }
    }
    return plane_added;
  }

  bool AddPlaneHKL(const xstructure& str,
                   const aurostd::xvector<double>& hkl,
                   const double& roughness,
                   const double& bbdistance,
                   const aurostd::xmatrix<double>& bbdistances,
                   const double& eps,
                   const int& jossmax,
                   vector<double>& plane,
                   vector<vector<double>>& planes,
                   const bool& osswrite1,
                   ostream& oss1,
                   const bool& osswrite2,
                   ostream& oss2) {
    return surface::AddPlaneHKL(str, hkl(1), hkl(2), hkl(3), roughness, bbdistance, bbdistances, eps, jossmax, plane, planes, osswrite1, oss1, osswrite2, oss2);
  }

  bool ReducePrintSurfaces(const xstructure& str,
                           const double& eps,
                           vector<vector<double>>& planesreducible,
                           vector<vector<double>>& planesirreducible,
                           vector<vector<uint>>& planesirreducible_images,
                           const bool& osswrite1,
                           ostream& oss1,
                           const bool& osswrite2,
                           ostream& oss2) {
    const uint num_types = str.num_each_type.size();
    const uint num_types_combinations = num_types * (num_types + 1) / 2;
    uint jossmax = 4;
    if (num_types_combinations > 1) {
      jossmax += num_types + num_types_combinations;
    }
    const vector<double> plane(3 + (1 + num_types) + (1 + num_types_combinations)); // 3 for hkl, 1 for all types and all combinations

    sort(planesreducible.begin(), planesreducible.end(), aurostd::_sort_double_value012());
    sort(planesreducible.begin(), planesreducible.end(), aurostd::_isort_double_value3());

    if (osswrite1) {
      oss1 << surface::PrintHKLSigma(num_types, num_types_combinations);
    }
    if (osswrite2) {
      oss2 << surface::PrintHKLSigma(num_types, num_types_combinations);
    }
    // THIS ROUTINE SHOULD BE CHANCED IN MULTITHREADS TO SPEED UP THE REDUCTION
    for (size_t i = 0; i < planesreducible.size(); i++) {
      const aurostd::xvector<double> rhkl(3);
      rhkl(1) = planesreducible[i][0];
      rhkl(2) = planesreducible[i][1];
      rhkl(3) = planesreducible[i][2];
      const double rdensity = planesreducible[i][3];
      const double rbbonds = planesreducible[i][4];
      if (rdensity > 0.0) {
        // must make sense
        bool hkl_found = false;
        for (size_t ii = 0; ii < planesirreducible.size() && !hkl_found; ii++) {
          // check if triplet of vectors are equivalent !
          const aurostd::xvector<double> ihkl(3);
          ihkl(1) = planesirreducible[ii][0];
          ihkl(2) = planesirreducible[ii][1];
          ihkl(3) = planesirreducible[ii][2];
          const double idensity = planesirreducible[ii][3];
          const double ibbonds = planesirreducible[ii][4];
          if (std::abs(rdensity - idensity) < eps && std::abs(rbbonds - ibbonds) < eps) {
            // check that the density and bbonds makes sense
            for (size_t sg = 0; sg < str.fgroup.size() && !hkl_found; sg++) {
              if (aurostd::modulus(str.fgroup[sg].ftau) < eps) {
                // only for non shift, otherwise no meaning
                if (aurostd::modulus(str.fgroup[sg].Uf * ihkl - rhkl) < eps) {
                  // same but faster than the other one
                  hkl_found = true;
                }
              }
            }
          }
        }
        if (hkl_found == false) {
          // new irreducible operation, generate and save it
          planesirreducible.push_back(planesreducible[i]);
          planesirreducible_images.emplace_back(0);
        }
        planesirreducible_images.back().push_back(i);
      }
    }

    for (size_t i = 0; i < planesirreducible.size(); i++) {
      // new irreducible operation, generate and save it
      for (uint j = 0; j <= jossmax; j++) {
        if (osswrite1) {
          oss1 << (planesirreducible[i].at(j) >= 0 ? "  " : " ") << (std::abs(planesirreducible[i].at(j)) < 10.0 ? " " : "") << planesirreducible[i].at(j) << " ";
        }
      }
      for (uint j = 0; j <= 4; j++) {
        if (osswrite2) {
          oss2 << (planesirreducible[i].at(j) >= 0 ? "  " : " ") << (std::abs(planesirreducible[i].at(j)) < 10.0 ? " " : "") << planesirreducible[i].at(j) << " ";
        }
      }
      if (osswrite1) {
        oss1 << planesirreducible.size() << " ";
      }
      if (osswrite2) {
        oss2 << planesirreducible.size() << " ";
      }
      if (osswrite1) {
        oss1 << endl;
      }
      if (osswrite2) {
        oss2 << endl;
      }
      if (osswrite1) {
        oss1 << "           equivalent family " << endl;
        for (size_t k = 0; k < planesirreducible_images.at(i).size(); k++) {
          for (uint j = 0; j < 3; j++) {
            oss1 << "" << (planesreducible.at(planesirreducible_images.at(i).at(k)).at(j) >= 0 ? "  " : " ") << (std::abs(planesreducible.at(planesirreducible_images.at(i).at(k)).at(j)) < 10.0 ? " " : "")
                 << planesreducible.at(planesirreducible_images.at(i).at(k)).at(j) << " ";
          }
          oss1 << endl;
        }
      }
    }
    // routine to make the planes as POSITIVE as POSSIBLE
    // ------------------------------------------------------
    // CODE FOR NUM_TYPES AND NUM_TYPES_COMBINATIONS
    return true;
  }

  // ------------------------------------------------------------------------------------------------------------------
  // only one HKL
  bool GetSurfaceHKL(const xstructure& _str, _aflags& aflags, const aurostd::xvector<double>& iparams, vector<vector<double>>& planesreducible, vector<vector<double>>& planesirreducible, ostream& oss) {
    const int num_types = _str.num_each_type.size();
    const int num_types_combinations = num_types * (num_types + 1) / 2;
    int jossmax = 4;
    const int mode = iparams.rows;
    if (num_types_combinations > 1) {
      jossmax += num_types + num_types_combinations;
    }
    vector<double> plane(3 + (1 + num_types) + (1 + num_types_combinations));
    // 3 for hkl, 1 for all types and all combinations
    const xstructure str = ReScale(BringInCell(_str), 1.0);
    const double roughness = _eps_ / 2.0;

    // vectors for symmetry search

    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(oss_short_precision_aflow_surface);

    const double nndist = surface::GetNNeighbors(str);
    const aurostd::xmatrix<double> nndists(num_types - 1, num_types - 1, 0, 0);
    aurostd::xmatrix<double> bbdistances(num_types - 1, num_types - 1, 0, 0);
    for (int it1 = 0; it1 < num_types; it1++) {
      for (int it2 = it1; it2 < num_types; it2++) {
        nndists(it1, it2) = surface::GetNNeighbors(str, it1, it2);
      }
    }

    //  if(mode==3 || mode==4) { // all three are given
    if (mode != 3 && mode != 4) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "only mode 3 and 4 are defined", _INPUT_ILLEGAL_);
    } // CO20200624
    if (mode == 3 || mode == 4) {
      double bbfrac = BBFRAC;
      constexpr bool Krun = true;
      if (mode == 3) {
        bbfrac = BBFRAC;
      }
      if (mode == 4) {
        bbfrac = iparams(4);
      }
      const double bbdistance = bbfrac * nndist;
      bbdistances = bbfrac * nndists;

      oss << "HKL CALCULATION" << endl;
      oss << surface::PrintNNdists(num_types, num_types_combinations, bbfrac, nndist, nndists);
      oss << surface::PrintHKLSigmaBB(num_types, num_types_combinations, bbfrac, bbdistance, bbdistances);
      // double h,k,l;
      const aurostd::xvector<double> hkl = iparams;
      // hkl
      plane.at(0) = hkl(1);
      plane.at(1) = hkl(2);
      plane.at(2) = hkl(3);
      // density
      const double density = surface::GetPlaneDensityAtoms(str, hkl, roughness);
      int i = 3;
      plane.at(i++) = density;
      if (num_types > 1) {
        for (int it1 = 0; it1 < num_types; it1++) {
          plane.at(i++) = surface::GetPlaneDensityAtoms(str, hkl, roughness, it1);
        }
      }
      // bbonds
      const double bbonds = surface::GetPlaneDensityBBonds(str, hkl, roughness, bbdistance);
      plane.at(i++) = bbonds;
      if (num_types > 1) {
        for (int it1 = 0; it1 < num_types; it1++) {
          for (int it2 = it1; it2 < num_types; it2++) {
            plane.at(i++) = surface::GetPlaneDensityBBonds(str, hkl, roughness, bbdistances(it1, it2), it1, it2);
          }
        }
      }
      for (int j = 0; j <= jossmax; j++) {
        oss << (plane[j] >= 0 ? "  " : " ") << (std::abs(plane[j]) < 10.0 ? " " : "") << plane[j] << " ";
      }
      oss << endl;
      planesreducible.push_back(plane);   // save them
      planesirreducible.push_back(plane); // save them
      return Krun;
    }
    return false;
  }

  // ------------------------------------------------------------------------------------------------------------------
  // Search HKL the trivial/simple/complete
  bool GetSurfaceHKLSearch(const xstructure& _str,
                           _aflags& aflags,
                           const aurostd::xvector<double>& iparams,
                           vector<vector<double>>& planesreducible,
                           vector<vector<double>>& planesirreducible,
                           vector<vector<uint>>& planesirreducible_images,
                           ostream& oss,
                           const string& smode) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    bool search_trivial = false;
    bool search_simple = false;
    bool search_complete = false;

    if (smode != "HKL_SEARCH_TRIVIAL" && smode != "HKL_SEARCH_SIMPLE" && smode != "HKL_SEARCH_COMPLETE") {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "unknown mode", _INPUT_ILLEGAL_);
    } // CO20200624
    if (smode == "HKL_SEARCH_TRIVIAL") {
      search_trivial = true;
    }
    if (smode == "HKL_SEARCH_SIMPLE") {
      search_simple = true;
    }
    if (smode == "HKL_SEARCH_COMPLETE") {
      search_complete = true;
    }
    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: SURFACE begin " << endl;
    }
    const double eps = _eps_;
    const double roughness = _eps_ / 2.0;
    double bbdistance;
    double step = 0;
    double hklmax = 0;
    aurostd::xvector<double> dims_hklmax(3);
    int jossmax = 4;
    const int mode = iparams.rows;
    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: [1] mode=" << mode << endl;
    }
    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: [1] iparams=" << iparams << endl;
    }

    const int num_types = _str.num_each_type.size();
    const int num_types_combinations = num_types * (num_types + 1) / 2;
    if (num_types_combinations > 1) {
      jossmax += num_types + num_types_combinations;
    }
    vector<double> plane(3 + (1 + num_types) + (1 + num_types_combinations));
    // 3 for hkl, 1 for all types and all combinations

    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: [1] " << endl;
    }

    string banner = "--------------------------------------------------------------";
    for (int i = 1; i <= num_types + num_types_combinations && num_types > 1; i++) {
      banner += "------------";
    }

    xstructure str = ReScale(BringInCell(_str), 1.0);
    double bbfrac = BBFRAC;

    // a1,a2,a3 are the rows of the lattice matrix
    const aurostd::xvector<double> a1 = str.lattice(1);
    const aurostd::xvector<double> a2 = str.lattice(2);
    const aurostd::xvector<double> a3 = str.lattice(3);

    vector<aurostd::xvector<double>*> grid_atoms_cpos;
    ofstream FileDevNull("/dev/null");
    constexpr bool PFSWRITE = true;
    constexpr bool FFFflag = true;
    bool hkl_found;

    aflags.QUIET = true;
    constexpr bool OSSWRITE = false; // true;

    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(oss_short_precision_aflow_surface);

    ofstream FFF;
    const string FileNameSURFACE = FFFflag ? DEFAULT_AFLOW_SURFACE_OUT : "/dev/null";

    FFF.open(FileNameSURFACE.c_str(), std::ios::out);
    FFF.setf(std::ios::fixed, std::ios::floatfield);
    FFF.precision(oss_short_precision_aflow_surface);

    const double nndist = surface::GetNNeighbors(str);
    const aurostd::xmatrix<double> nndists(num_types - 1, num_types - 1, 0, 0);
    aurostd::xmatrix<double> bbdistances(num_types - 1, num_types - 1, 0, 0);
    for (int it1 = 0; it1 < num_types; it1++) {
      for (int it2 = it1; it2 < num_types; it2++) {
        nndists(it1, it2) = surface::GetNNeighbors(str, it1, it2);
      }
    }

    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: [2] " << endl;
    }

    if (search_trivial) {
      if (mode == 0 || mode == 1 || mode == 2 || mode == 3) { // seek
        switch (mode) {
          case 0: {
            hklmax = HKLDEF;
            bbfrac = BBFRAC;
            step = 1.0;
          }
          case 1: {
            hklmax = std::abs(iparams(1));
            bbfrac = BBFRAC;
            step = 1.0;
          }
          case 2: {
            hklmax = std::abs(iparams(1));
            bbfrac = iparams(2);
            step = std::abs(iparams(3));
          }
          case 3: {
            hklmax = std::abs(iparams(1));
            bbfrac = iparams(2);
            step = std::abs(iparams(3));
          }
          default:;
        }
        dims_hklmax(1) = hklmax;
        dims_hklmax(2) = hklmax;
        dims_hklmax(3) = hklmax;
        bbdistance = bbfrac * nndist;
        bbdistances = bbfrac * nndists;

        oss << banner << endl; // ----------------------------------------------------------------
        oss << aflow::Banner("BANNER_TINY") << endl;
        oss << banner << endl; // ----------------------------------------------------------------
        oss << "HKL CALCULATION" << endl;
        if (search_simple) {
          oss << "SIMPLE SEARCH" << endl;
        }
        if (search_trivial) {
          oss << "TRIVIAL SEARCH" << endl;
        }
        if (search_complete) {
          oss << "COMPLETE SEARCH" << endl;
        }
        str.LatticeReduction_avoid = true;                           // DOES NOT DO LATTICE REDUCTION // NIGGLI and MINK
        str.sgroup_radius = 1.05 * RadiusSphereLattice(str.lattice); // CO20171024 - new sym framework
        _kflags kflags;
        pflow::defaultKFlags4SymCalc(kflags, true); // CO20171024 - new sym framework
        pflow::defaultKFlags4SymWrite(kflags, PFSWRITE);
        kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;                                   // CO20171024 - new sym framework
        pflow::PerformFullSymmetry(str, FileDevNull, aflags, kflags, OSSWRITE, oss); // CO20171024 - new sym framework
        //   oss << banner << endl; // ----------------------------------------------------------------
        oss << surface::PrintNNdists(num_types, num_types_combinations, bbfrac, nndist, nndists);
        oss << "hklmax=" << hklmax << endl;
        oss << "bbdistance=" << bbdistance / nndist << "  surface::GetNNeighbors=" << surface::GetNNeighbors(str) << " real_bbdistance=" << bbdistance << endl;
        oss << "step=" << step << endl;
        oss << "SCANNING TRIVIAL PLANES" << endl;
        oss << surface::PrintHKLSigmaBB(num_types, num_types_combinations, bbfrac, bbdistance, bbdistances);
        // FFF
        FFF << banner << endl; // ----------------------------------------------------------------
        FFF << aflow::Banner("BANNER_TINY") << endl;
        FFF << banner << endl; // ----------------------------------------------------------------
        FFF << "HKL CALCULATION" << endl;
        if (search_trivial) {
          FFF << "TRIVIAL SEARCH" << endl;
        }
        if (search_simple) {
          FFF << "SIMPLE SEARCH" << endl;
        }
        if (search_complete) {
          FFF << "COMPLETE SEARCH" << endl;
        }
        FFF << surface::PrintNNdists(num_types, num_types_combinations, bbfrac, nndist, nndists);
        FFF << "hklmax=" << hklmax << endl;
        FFF << "bbdistance=" << bbdistance / nndist << "  surface::GetNNeighbors=" << surface::GetNNeighbors(str) << " real_bbdistance=" << bbdistance << endl;
        FFF << "step=" << step << endl;
        FFF << "SCANNING TRIVIAL PLANES" << endl;
        FFF << surface::PrintHKLSigmaBB(num_types, num_types_combinations, bbfrac, bbdistance, bbdistances);
        //
        // THIS ROUTINE IS VERY SLOW AND SHOULD BE MADE MULTI-THREADS
        vector<aurostd::xvector<double>> vhkl;
        for (double h = dims_hklmax[1]; h >= -dims_hklmax[1]; h -= step) {
          for (double k = dims_hklmax[2]; k >= -dims_hklmax[2]; k -= step) {
            for (double l = dims_hklmax[3]; l >= -dims_hklmax[3]; l -= step) {
              if (std::abs(h) < eps) {
                h = 0.0;
              }
              if (std::abs(k) < eps) {
                k = 0.0;
              }
              if (std::abs(l) < eps) {
                l = 0.0;
              }
              vhkl.emplace_back(aurostd::xvector{h, k, l});
            }
          }
        }
        for (size_t i = 0; i < vhkl.size(); i++) {
          hkl_found = surface::AddPlaneHKL(str, vhkl[i], roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
        }
      }
    }

    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: [3] " << endl;
    }

    if (search_simple || search_complete) {
      if (mode == 0 || mode == 1 || mode == 2 || mode == 3 || mode == 4) {
        double rrfrac = RRFRAC;
        // seek
        if (LDEBUG) {
          cerr << "surface::GetSurfaceHKLSearch: [3] mode=" << mode << endl;
        }
        if (search_simple) {
          step = 1.0;
        };
        if (search_complete) {
          step = 1 / 3.0 / 4.0;
        };

        switch (mode) {
          case 0: {
            rrfrac = RRFRAC;
            bbfrac = BBFRAC;
            hklmax = -1;
            break;
          }
          case 1: {
            rrfrac = std::abs(iparams(1));
            bbfrac = BBFRAC;
            hklmax = -1;
            break;
          }
          case 2: {
            rrfrac = std::abs(iparams(1));
            bbfrac = std::abs(iparams(2));
            hklmax = -1;
            break;
          }
          case 3: {
            rrfrac = std::abs(iparams(1));
            bbfrac = std::abs(iparams(2));
            hklmax = std::abs(iparams(3));
            break;
          }
          case 4: {
            rrfrac = std::abs(iparams(1));
            bbfrac = std::abs(iparams(2));
            hklmax = std::abs(iparams(3));
            step = std::abs(iparams(4));
            break;
          }
          default: break;
        }
        bbdistance = bbfrac * nndist;
        bbdistances = bbfrac * nndists;
        const double radius = rrfrac * modulus(a1 + a2 + a3) + bbdistance;
        const aurostd::xvector<int> dims = LatticeDimensionSphere(str.lattice, radius);
        int imin;
        int imax;
        int jmin;
        int jmax;
        int kmin;
        int kmax;
        if (search_simple) {
          imin = 0;
          jmin = 0;
          kmin = 0;
          imax = dims(1);
          jmax = dims(2);
          kmax = dims(3);
        }
        if (search_complete) {
          imin = -dims(1);
          jmin = -dims(2);
          kmin = -dims(3);
          imax = dims(1);
          jmax = dims(2);
          kmax = dims(3);
        }
        if (hklmax < 0) {
          // xvector needs an iterator so we can assign
          dims_hklmax(1) = dims(1);
          dims_hklmax(2) = dims(2);
          dims_hklmax(3) = dims(3);
        } else {
          dims_hklmax = {hklmax, hklmax, hklmax};
        }

        oss << banner << endl; // ----------------------------------------------------------------
        oss << aflow::Banner("BANNER_TINY") << endl;
        oss << banner << endl; // ----------------------------------------------------------------
        if (search_simple) {
          oss << "SIMPLE SEARCH" << endl;
        }
        if (search_trivial) {
          oss << "TRIVIAL SEARCH" << endl;
        }
        if (search_complete) {
          oss << "COMPLETE SEARCH" << endl;
        }
        str.LatticeReduction_avoid = true;
        str.sgroup_radius = 1.05 * RadiusSphereLattice(str.lattice);
        _kflags kflags;
        pflow::defaultKFlags4SymCalc(kflags, true); // CO20171024 - new sym framework
        pflow::defaultKFlags4SymWrite(kflags, PFSWRITE);
        kflags.KBIN_SYMMETRY_SGROUP_WRITE = false;                                   // CO20171024 - new sym framework
        pflow::PerformFullSymmetry(str, FileDevNull, aflags, kflags, OSSWRITE, oss); // CO20171024 - new sym framework
        // oss
        oss << banner << endl; // ----------------------------------------------------------------
        oss << "HKL CALCULATION" << endl;
        if (search_trivial) {
          oss << "TRIVIAL SEARCH" << endl;
        }
        if (search_simple) {
          oss << "SIMPLE SEARCH" << endl;
        }
        if (search_complete) {
          oss << "COMPLETE SEARCH" << endl;
        }
        oss << surface::PrintNNdists(num_types, num_types_combinations, bbfrac, nndist, nndists);
        oss << "CUTOFF rrfrac=" << rrfrac << endl;
        oss << "BOND   bbfrac=" << bbfrac << endl;
        oss << "HKLMAX hklmax=" << hklmax << endl;
        oss << "STEP   step=  " << step << endl;
        oss << " radius=" << radius << endl;
        oss << " dims=(" << dims(1) << "," << dims(2) << "," << dims(3) << ")" << endl;
        oss << " dims_hklmax=(" << dims_hklmax[1] << "," << dims_hklmax[2] << "," << dims_hklmax[3] << ")" << endl;
        oss << "SCANNING " << endl;
        // FFF
        FFF << banner << endl; // ----------------------------------------------------------------
        FFF << aflow::Banner("BANNER_TINY") << endl;
        FFF << banner << endl; // ----------------------------------------------------------------
        FFF << "HKL CALCULATION" << endl;
        if (search_trivial) {
          FFF << "TRIVIAL SEARCH" << endl;
        }
        if (search_simple) {
          FFF << "SIMPLE SEARCH" << endl;
        }
        if (search_complete) {
          FFF << "COMPLETE SEARCH" << endl;
        }
        FFF << surface::PrintNNdists(num_types, num_types_combinations, bbfrac, nndist, nndists);
        FFF << "CUTOFF rrfrac=" << rrfrac << endl;
        FFF << "BOND   bbfrac=" << bbfrac << endl;
        FFF << "HKLMAX hklmax=" << hklmax << endl;
        FFF << "STEP   step=  " << step << endl;
        FFF << " radius=" << radius << endl;
        FFF << " dims=(" << dims(1) << "," << dims(2) << "," << dims(3) << ")" << endl;
        FFF << " dims_hklmax=(" << dims_hklmax[1] << "," << dims_hklmax[2] << "," << dims_hklmax[3] << ")" << endl;
        FFF << "SCANNING " << endl;
        // ----

        vector<int> grid_atoms_number;
        for (int i = imin; i <= imax; i++) {
          for (int j = jmin; j <= jmax; j++) {
            for (int k = kmin; k <= kmax; k++) {
              for (size_t iat = 0; iat < str.atoms.size(); iat++) {
                const aurostd::xvector<double> rrr = i * a1 + j * a2 + k * a3 + str.atoms.at(iat).cpos;
                if (modulus(rrr) <= radius && modulus(rrr) > eps) {
                  grid_atoms_cpos.emplace_back(new aurostd::xvector<double>(rrr));
                  grid_atoms_number.push_back(str.atoms.at(iat).basis);
                }
              }
            }
          }
        }
        const uint grid_atoms_size = grid_atoms_cpos.size();
        oss << "grid_atoms_size=" << grid_atoms_size << endl;
        FFF << "grid_atoms_size=" << grid_atoms_size << endl;
        oss << banner << endl; // ----------------------------------------------------------------
        FFF << banner << endl; // ----------------------------------------------------------------
        oss << surface::PrintHKLSigmaBB(num_types, num_types_combinations, bbfrac, bbdistance, bbdistances);
        FFF << surface::PrintHKLSigmaBB(num_types, num_types_combinations, bbfrac, bbdistance, bbdistances);
        //
        // extra juice !
        // only if complete
        oss << "SCANNING TRIVIAL PLANES" << endl;
        FFF << "SCANNING TRIVIAL PLANES" << endl;
        for (double h = -dims_hklmax[1]; h <= dims_hklmax[1]; h += step) {
          for (double k = -dims_hklmax[2]; k <= dims_hklmax[2]; k += step) {
            for (double l = -dims_hklmax[3]; l <= dims_hklmax[3]; l += step) {
              if (std::abs(h) < eps) {
                h = 0.0;
              }
              if (std::abs(k) < eps) {
                k = 0.0;
              }
              if (std::abs(l) < eps) {
                l = 0.0;
              }
              hkl_found = surface::AddPlaneHKL(str, h, k, l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
            }
          }
        }
        oss << "SCANNING COMPLICATE PLANES" << endl;
        FFF << "SCANNING COMPLICATE PLANES" << endl;
        for (size_t iat1 = 0; iat1 < grid_atoms_size; iat1++) {
          const int number1 = grid_atoms_number[iat1];
          for (size_t iat2 = iat1 + 1; iat2 < grid_atoms_size; iat2++) {
            if (iat2 != iat1) {
              const int number2 = grid_atoms_number[iat2];
              if (number2 == number1) {
                for (size_t iat3 = iat2 + 1; iat3 < grid_atoms_size; iat3++) {
                  if (iat3 != iat1 && iat3 != iat2) {
                    const int number3 = grid_atoms_number[iat3];
                    if (number3 == number2 && number2 == number1) {
                      // all threee numbers identical
                      const aurostd::xvector<double> rrr1 = *grid_atoms_cpos[iat1];
                      const aurostd::xvector<double> rrr2 = *grid_atoms_cpos[iat2];
                      const aurostd::xvector<double> rrr3 = *grid_atoms_cpos[iat3];
                      const double determinant = det(rrr1, rrr2, rrr3); // is zero if origin goes through the triangle (better avoid)
                      if (std::abs(determinant) > eps) {
                        const double area = surface::TriangleArea(rrr1, rrr2, rrr3);
                        if (area > eps) {
                          double a;
                          double b;
                          double c;
                          double d;
                          PlaneGetABCD(a, b, c, d, rrr1, rrr2, rrr3); //  abs(d/sqrt(a*a+b*b+c*c));
                          if (std::abs(d / std::sqrt(a * a + b * b + c * c)) > eps) {
                            const xvector<double> hkl = PlaneGetHKL(rrr1, rrr2, rrr3, a1, a2, a3);
                            double h = hkl(1);
                            double k = hkl(2);
                            double l = hkl(3);
                            if ((std::abs(h) >= 1.0 || std::abs(h) < eps) && (std::abs(k) >= 1.0 || std::abs(k) < eps) && (std::abs(l) >= 1.0 || std::abs(l) < eps)) {
                              if (std::abs(h) < eps) {
                                h = 0.0;
                              }
                              if (std::abs(k) < eps) {
                                k = 0.0;
                              }
                              if (std::abs(l) < eps) {
                                l = 0.0;
                              }
                              if (search_simple) {
                                hkl_found = surface::AddPlaneHKL(str, h, k, l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, h, k, l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, -h, k, l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, h, -k, l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, h, k, -l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, h, -k, -l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, -h, k, -l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                              if (search_complete) {
                                hkl_found = surface::AddPlaneHKL(str, -h, -k, l, roughness, bbdistance, bbdistances, eps, jossmax, plane, planesreducible, true, oss, FFFflag, FFF);
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: [4] " << endl;
    }

    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    oss << "REDUCIBLE: " << planesreducible.size() << endl;
    FFF << "REDUCIBLE: " << planesreducible.size() << endl;
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    // reduce --------------------------------------------------------------------
    surface::ReducePrintSurfaces(str, eps, planesreducible, planesirreducible, planesirreducible_images, true, oss, FFFflag, FFF);
    // print end  ----------------------------------------------------------------
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    oss << "IRREDUCIBLE= " << planesirreducible.size() << "  REDUCIBLE= " << planesreducible.size() << endl;
    FFF << "IRREDUCIBLE= " << planesirreducible.size() << "  REDUCIBLE= " << planesreducible.size() << endl;
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    oss << aflow::Banner("BANNER_TINY") << endl;
    FFF << aflow::Banner("BANNER_TINY") << endl;
    oss << banner << endl; // ----------------------------------------------------------------
    FFF << banner << endl; // ----------------------------------------------------------------
    // END CLEAN EVERYTHING ------------------------------------------------------------------
    FFF.close();
    for (size_t i = 0; i < grid_atoms_cpos.size(); i++) {
      delete grid_atoms_cpos[i];
    }
    grid_atoms_cpos.clear();

    if (LDEBUG) {
      cerr << "surface::GetSurfaceHKLSearch: END " << endl;
    }

    return hkl_found;
  }

} // namespace surface

// ********************************************************************************************************/
// ********************************************************************************************************/
// RC2009: Depending on the initial input POSCAR, the obtained 2 basis vectors in the slab layer may be not as minimum as possible.
//        It should be corrected afterwards, by applying the general methods of unit cell minimization to the obtained here final slab POSCAR.
// RC2009: The slab uc basis vectors coordinates are presented wrt initial POSCAR Cartesian system. So it is easy to find the position of slab uc basis vectors wrt initial POSCAR sites.
// RC2009: If (a) parent latice is undistorted cubic-like and (b) Cartesian axises (wrt which the uc basis vectors of initial POSCAR are determined) are directed along cubic edges
//        then the Cartesian coordinates of third (normal to layers) basis bector of slab uc (S3) determine the slab Miller indices wrt parent CUBIC unit cell.
//        E.g. for L10 with 4at/uc and 2at/us POSCAR definitions: (111)of4/uc=(101)of2/uc, (010)of4/uc=(110)of2/uc, (001)of4/uc=(001)of2/uc, (110)of4/uc=(100)of2/uc

#define _slab_file_OUT string("POSCAR_slab_hkl")
#define _slab_file2_OUT string("Plane.dat")
#define _slab_epsilon_ double(1e-4)        // Small nonzero number

// ****************************************************************************************************************/
namespace slab {
  double VectorAbsValue(int Layer, int NinLayer /*IN*/, const aurostd::xmatrix<double>& UnitCellVector, const vector<vector<vector<int>>>& LayerSitesDirCoords) {
    int i1;
    int j1;
    double AbsValue = 0;
    for (j1 = 1; j1 <= 3; j1++) {
      double x = 0;
      for (i1 = 1; i1 <= 3; i1++) {
        x += LayerSitesDirCoords[Layer][NinLayer][i1] * UnitCellVector[i1][j1];
      }
      AbsValue += x * x;
    }
    AbsValue = sqrt(AbsValue);
    return AbsValue;
  }
} // namespace slab

// ****************************************************************************************************************/
namespace slab {
  double VectorScalarMult(int Layer1, int NinLayer1, int Layer2, int NinLayer2 /*IN*/, const aurostd::xmatrix<double>& UnitCellVector, const vector<vector<vector<int>>>& LayerSitesDirCoords) {
    double Help1;
    double Help2;
    int i1;
    int j1;
    double ScalarMult = 0;
    for (j1 = 1; j1 <= 3; j1++) {
      Help1 = 0;
      for (i1 = 1; i1 <= 3; i1++) {
        Help1 += LayerSitesDirCoords[Layer1][NinLayer1][i1] * UnitCellVector[i1][j1];
      }
      Help2 = 0;
      for (i1 = 1; i1 <= 3; i1++) {
        Help2 += LayerSitesDirCoords[Layer2][NinLayer2][i1] * UnitCellVector[i1][j1];
      }
      ScalarMult += Help1 * Help2;
    }
    return ScalarMult;
  }
} // namespace slab

// ****************************************************************************************************************/
namespace slab {
  double CosAngle(int Layer1, int NinLayer1, int Layer2, int NinLayer2 /*IN*/, const aurostd::xmatrix<double>& UnitCellVector, const vector<vector<vector<int>>>& LayerSitesDirCoords) {
    double CosAngleOUT;
    double ScalarMult;
    double AbsValue1;
    double AbsValue2;
    ScalarMult = slab::VectorScalarMult(Layer1, NinLayer1, Layer2, NinLayer2, UnitCellVector, LayerSitesDirCoords);
    AbsValue1 = slab::VectorAbsValue(Layer1, NinLayer1, UnitCellVector, LayerSitesDirCoords);
    AbsValue2 = slab::VectorAbsValue(Layer2, NinLayer2, UnitCellVector, LayerSitesDirCoords);
    CosAngleOUT = ScalarMult / (AbsValue1 * AbsValue2);
    return CosAngleOUT;
  }
} // namespace slab

// ****************************************************************************************************************/
namespace slab {
  double hkl_CartCoord_Length(aurostd::xvector<double>& hkl_CartCoord, const aurostd::xmatrix<double>& UnitCellVector, const aurostd::xvector<double>& hkl) {
    // ReciprUnitCellVector,hkl_CartCoord_Length are in 2PI/LattPar[0] (VolumeReciprUC in LattPar[0]^3) units
    int i;
    int j;
    double VolumeReciprUC;
    const aurostd::xmatrix<double> ReciprUnitCellVector(3, 3);
    double hkl_Length;
    VolumeReciprUC = UnitCellVector(1, 1) * (UnitCellVector(2, 2) * UnitCellVector(3, 3) - UnitCellVector(3, 2) * UnitCellVector(2, 3)) -
                     UnitCellVector(1, 2) * (UnitCellVector(2, 1) * UnitCellVector(3, 3) - UnitCellVector(3, 1) * UnitCellVector(2, 3)) +
                     UnitCellVector(1, 3) * (UnitCellVector(2, 1) * UnitCellVector(3, 2) - UnitCellVector(3, 1) * UnitCellVector(2, 2));
    ReciprUnitCellVector(1, 1) = UnitCellVector(2, 2) * UnitCellVector(3, 3) - UnitCellVector(3, 2) * UnitCellVector(2, 3);
    ReciprUnitCellVector(1, 2) = -(UnitCellVector(2, 1) * UnitCellVector(3, 3) - UnitCellVector(3, 1) * UnitCellVector(2, 3));
    ReciprUnitCellVector(1, 3) = UnitCellVector(2, 1) * UnitCellVector(3, 2) - UnitCellVector(3, 1) * UnitCellVector(2, 2);

    ReciprUnitCellVector(2, 1) = UnitCellVector(3, 2) * UnitCellVector(1, 3) - UnitCellVector(1, 2) * UnitCellVector(3, 3);
    ReciprUnitCellVector(2, 2) = -(UnitCellVector(3, 1) * UnitCellVector(1, 3) - UnitCellVector(1, 1) * UnitCellVector(3, 3));
    ReciprUnitCellVector(2, 3) = UnitCellVector(3, 1) * UnitCellVector(1, 2) - UnitCellVector(1, 1) * UnitCellVector(3, 2);

    ReciprUnitCellVector(3, 1) = UnitCellVector(1, 2) * UnitCellVector(2, 3) - UnitCellVector(2, 2) * UnitCellVector(1, 3);
    ReciprUnitCellVector(3, 2) = -(UnitCellVector(1, 1) * UnitCellVector(2, 3) - UnitCellVector(2, 1) * UnitCellVector(1, 3));
    ReciprUnitCellVector(3, 3) = UnitCellVector(1, 1) * UnitCellVector(2, 2) - UnitCellVector(2, 1) * UnitCellVector(1, 2);
    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        ReciprUnitCellVector(i, j) /= VolumeReciprUC;
      }
    }
    // cerr << "Reciprocal Basis: "; for(i=1;i<=3;i++) { cerr << "("; for(j=1;j<=3;j++) {cerr << ReciprUnitCellVector(i,j) << ","; } cerr << ")"; }; cerr << endl << endl;
    for (j = 1; j <= 3; j++) {
      hkl_CartCoord[j] = 0;
      for (i = 1; i <= 3; i++) {
        hkl_CartCoord[j] += hkl(i) * ReciprUnitCellVector(i, j);
      }
    }
    hkl_Length = 0;
    for (j = 1; j <= 3; j++) {
      hkl_Length += hkl_CartCoord[j] * hkl_CartCoord[j];
    }
    hkl_Length = sqrt(hkl_Length);
    return hkl_Length;
  } // END Reciprocal_Lattice_Calc //
} // namespace slab

// ********************************************************************************************************/
namespace slab {
  xstructure MAKE_SLAB(string options, istream& cin) {
    xstructure str_in(cin, IOAFLOW_AUTO);
    xstructure str_out("");
    str_out = MAKE_SLAB(options, str_in);
    // MAKE_SLAB(options,str_in);
    //   cout << str_out;
    return str_out;
  }
} // namespace slab

namespace slab {
  xstructure MAKE_SLAB(string options, xstructure& _str_in) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (tokens.size() < 3 || tokens.size() > 5) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --slab=h,k,l[,#filled_layers[,#vacuum layers]] < POSCAR");
    }
    int i = 0;
    int j = 0;
    int k = 0;
    //  options[0-2]="h" "k" "l" options[3-4]="NumFilledLayers" "NumEmptyLayers"

    // CO20180724 START - AddAtom() requires atom types AND names, so assign fake names if necessary, then remove later
    string specie;
    xstructure str_in = _str_in;
    bool assigning_fake_names = false;
    for (size_t i = 0; i < str_in.num_each_type.size(); i++) {
      specie = str_in.SpeciesLabel(i);
      if (specie.empty() || aurostd::substring2bool(specie, "name not given")) {
        assigning_fake_names = true;
        break;
      }
    }
    if (assigning_fake_names) {
      int iiatom = 0;
      for (size_t i = 0; i < str_in.num_each_type.size(); i++) {
        for (int j = 0; j < str_in.num_each_type[i]; j++) {
          str_in.atoms[iiatom].name = char('A' + i);
          str_in.atoms[iiatom].name_is_given = true;
          iiatom++;
        }
      }
    }
    // CO20180724 STOP - AddAtom() requires atom types AND names, so assign fake names if necessary, then remove later

    const aurostd::xvector<double> hkl(3);
    hkl(1) = 1;
    hkl(2) = 1;
    hkl(3) = 1;
    if (!tokens.empty()) {
      hkl(1) = aurostd::string2utype<double>(tokens.at(0));
    }
    if (tokens.size() >= 2) {
      hkl(2) = aurostd::string2utype<double>(tokens.at(1));
    }
    if (tokens.size() >= 3) {
      hkl(3) = aurostd::string2utype<double>(tokens.at(2));
    }

    string file_OUT = _slab_file_OUT;
    stringstream temp;
    temp << "POSCAR_slab_" << hkl(1) << hkl(2) << hkl(3);
    temp >> file_OUT;

    int NumFilledLayers = 2;
    int NumEmptyLayers = 3;  // number of filled/empty layers in slab
    const int SearchMax = 20 + NumFilledLayers + NumEmptyLayers;
    if (tokens.size() >= 4) {
      NumFilledLayers = aurostd::string2utype<int>(tokens.at(3));
    }
    if (tokens.size() >= 5) {
      NumEmptyLayers = aurostd::string2utype<int>(tokens.at(4));
    }

    const int In_Plane_Multiplication[3] = {0, 1, 1};    // In_Plane_Multiplication^2 = number of plane unit cells in slab unit cell
    vector<int> NumAtomsForElementUC(2, 0);
    vector<vector<vector<double>>> AtomDirCoords;
    vector<vector<vector<int>>> LayerSitesDirCoords;

    //   int const SearchMax=20+ NumFilledLayers + NumEmptyLayers;
    // [-SearchMax to SearchMax]^3 number of tested sites in lattice=12 + NumFilledLayers + NumEmptyLayers

    int NumSitesInPlane[NumFilledLayers + 1];
    int LayerBasis[3];
    int Layer0Basis[3];
    int t1 = 0;
    int t2 = 0;
    int VECTOR[3];
    int NumLayers;
    int NumBasisSitesInSlab;
    int Layer;
    int l[4];
    int Sum_Int;
    int Element;
    int ElementSite;
    double DET = 0.0;
    double s[3];
    double fractpart[3];
    double intpart;
    double Layer0BasisCosAngle = 1;
    double Layer0BasisCosAngleCurrent;
    double LastCandidateCoord3;
    double CurrCandidateCoord3;
    double hkl_Length;
    const aurostd::xvector<double> SiteCartCoord(3);
    const aurostd::xvector<double> RS(3);
    const aurostd::xvector<double> AbsValueSlab_Basis(3);
    const aurostd::xvector<double> SiteDirectCoordWRTslab(3);
    aurostd::xvector<double> hkl_CartCoord(3);
    const aurostd::xvector<double> SiteDirectCoord(3);
    const aurostd::xmatrix<double> MATRIX(3, 3);
    const aurostd::xmatrix<double> INVERSE_MATRIX(3, 3);
    const aurostd::xmatrix<double> Slab_Basis_CartCoord(3, 3);

    bool Continue;
    bool BasisFound;

    const int NumberElements = str_in.num_each_type.size();
    const string Title = str_in.title;
    const double LattPar = str_in.scale;

    const aurostd::xmatrix<double> UnitCellVector(3, 3);

    for (i = 1; i < 4; i++) {
      for (j = 1; j < 4; j++) {
        UnitCellVector(i, j) = str_in.lattice(i, j);
      }
    }

    for (i = 1; i < NumberElements + 1; i++) {
      NumAtomsForElementUC[i] = str_in.num_each_type[i - 1];
    }

    // i=1;NumberElements=0;
    // while (iss >> NumAtomsForElementUC[i]) {i++; NumberElements++; NumAtomsForElementUC.push_back(0);};
    // here

    int iatom = -1;
    AtomDirCoords.resize(NumberElements + 1);
    for (k = 1; k <= NumberElements; k++) {
      AtomDirCoords[k].resize(NumAtomsForElementUC[k] + 1);
    }
    for (k = 1; k <= NumberElements; k++) {
      for (i = 1; i <= NumAtomsForElementUC[k]; i++) {
        AtomDirCoords[k][i].resize(4);
      }
    }
    for (k = 1; k <= NumberElements; k++) {
      for (i = 1; i <= NumAtomsForElementUC[k]; i++) {
        iatom++;
        for (j = 1; j <= 3; j++) {
          AtomDirCoords[k][i][j] = str_in.atoms.at(iatom).fpos[j];
        }
      }
    }

    //  int itmp=0;
    //  AtomDirCoords.resize(NumberElements+1);
    //  for(i=1;i<=NumberElements+1;i++) {
    // itmp=NumAtomsForElementUC[i];
    // AtomDirCoords[i].resize(itmp);
    // for(j=1;j<=itmp;j++)
    //  AtomDirCoords[i][j].resize(4);
    //  }
    //  int iatom=0;
    //  for(i=1;i<=NumberElements;i++) {
    // for(j=1;j<=NumAtomsForElementUC[i];j++) {
    //  iatom++;
    //  for(k=1;k<=3;k++) {
    // AtomDirCoords[i][j][k]=str_in.atoms.at(iatom).fpos[k];
    //  }} }

    vector<int> NumSites(NumberElements + 1, 0);

    //---------------- Sorting (by layer number) of all sites in [-SearchMax to SearchMax]^3 box -----------------------------------------------//
    //  from equation Sum_i (l_i*n_i)=Layer number; {l_i}=sites coordinates wrt direct-lattice basis; {n_i}=hkl

    LayerSitesDirCoords.resize(NumFilledLayers + 1);

    for (Layer = 0; Layer <= NumFilledLayers; Layer++) {
      NumSitesInPlane[Layer] = 0;
    }

    for (l[1] = SearchMax; l[1] >= -SearchMax; l[1]--) {
      for (l[2] = SearchMax; l[2] >= -SearchMax; l[2]--) {
        for (l[3] = SearchMax; l[3] >= -SearchMax; l[3]--) {
          Sum_Int = l[1] * hkl(1) + l[2] * hkl(2) + l[3] * hkl(3);
          if ((Sum_Int > -1) && (Sum_Int < NumFilledLayers + 1)) {
            NumSitesInPlane[Sum_Int]++;
            LayerSitesDirCoords[Sum_Int].resize(NumSitesInPlane[Sum_Int] + 1);
            LayerSitesDirCoords[Sum_Int][NumSitesInPlane[Sum_Int]].resize(4);
            for (i = 1; i <= 3; i++) {
              LayerSitesDirCoords[Sum_Int][NumSitesInPlane[Sum_Int]][i] = l[i];
            }
          }
        }
      }
    }
    if (LDEBUG) {
      for (Layer = 0; Layer <= 0; Layer++) {
        cerr << "Layer= " << Layer << endl;
        for (j = 1; j <= NumSitesInPlane[Layer]; j++) {
          cerr << j << " " << LayerSitesDirCoords[Layer][j][1] << " " << LayerSitesDirCoords[Layer][j][2] << " " << LayerSitesDirCoords[Layer][j][3] << endl;
        }
      }
    }

    stringstream f2out;
    double x;
    for (k = 1; k <= NumSitesInPlane[0]; k++) {
      for (j = 1; j <= 3; j++) {
        x = 0.0;
        for (i = 1; i <= 3; i++) {
          x += LayerSitesDirCoords[0][k][i] * UnitCellVector(i, j);
        }
        f2out << x << " ";
      }
      f2out << k << endl;
    }
    aurostd::stringstream2file(f2out, _slab_file2_OUT);
    // return 1;
    // cin >> Help_string;
    //----------------  Searching for basis vectors in m=0 layer ------------------------------------------------------//
    //   Condition: any layer-site-vector (k-numbered) is a integer-linear superposition of two layer basis vectors (LayerBasis[1,2]-numbered in a search)
    //   The angle between basis vectors in layer is seek to be closest to 90 degree

    // cerr << "Searching for basis vectors in m=0 layer with closest to 90 degree inter-angle" << endl;

    BasisFound = false;
    LayerBasis[1] = 1;
    while (LayerBasis[1] <= NumSitesInPlane[0]) {// cerr << LayerBasis[1] << " / " << NumSitesInPlane[0] << endl;
      LayerBasis[2] = LayerBasis[1] + 1;
      while (LayerBasis[2] <= NumSitesInPlane[0]) {
        for (i = 1; i <= 2; i++) {
          for (j = 1; j <= 3; j++) {
            MATRIX[j][i] = LayerSitesDirCoords[0][LayerBasis[i]][j];
          }
        }
        Continue = true;
        i = 1;
        while (i <= 2 && Continue == true) {
          j = i + 1;
          while (j <= 3 && Continue == true) {
            DET = MATRIX[i][1] * MATRIX[j][2] - MATRIX[i][2] * MATRIX[j][1];
            if (aurostd::abs(DET) > _slab_epsilon_) {
              t1 = i;
              t2 = j;
              Continue = false;
            }
            j++;
          }
          i++;
        }

        // if(LayerBasis[1]==48 && LayerBasis[2]==49)
        //{  cerr << " DET= " << DET << endl;
        // for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][LayerBasis[i]][j] << " ";} cerr << endl;}
        // cerr << MATRIX[1][1] << " " << MATRIX[1][2] << endl << MATRIX[2][1] << " " << MATRIX[2][2] << endl;
        // }

        if (aurostd::abs(DET) > _slab_epsilon_) {
          INVERSE_MATRIX[1][1] = MATRIX[t2][2] / DET;
          INVERSE_MATRIX[1][2] = -MATRIX[t1][2] / DET;
          INVERSE_MATRIX[2][1] = -MATRIX[t2][1] / DET;
          INVERSE_MATRIX[2][2] = MATRIX[t1][1] / DET;
          VECTOR[1] = t1;
          VECTOR[2] = t2;
          Continue = true;
          k = 1;
          while (k <= NumSitesInPlane[0] && Continue == true) {
            for (i = 1; i <= 2; i++) {
              s[i] = 0;
              for (j = 1; j <= 2; j++) {
                s[i] += INVERSE_MATRIX(i, j) * LayerSitesDirCoords[0][k][VECTOR[j]];
              }
            }
            for (i = 1; i <= 2; i++) {
              fractpart[i] = modf(s[i], &intpart);
            }

            // if(aurostd::abs(fractpart[1]-1.0)<_slab_epsilon_ || aurostd::abs(fractpart[2]-1.0)<_slab_epsilon_) {for(i=1;i<=2;i++) { cerr << s[i] << " (" << fractpart[i] << ") " << endl; };}
            // cerr << s[1] << " (" << fractpart[1] << "); "  << s[2] << " (" << fractpart[2] << ") " << endl; cin >> Help_string;
            // cerr << aurostd::abs(fractpart[1]) << "absFrac-_slab_epsilon_" << _slab_epsilon_ << endl;
            //   if((aurostd::abs(fractpart[1])<_slab_epsilon_) &&  (aurostd::abs(fractpart[2])<_slab_epsilon_)) { cerr << LayerBasis[1] << " =LayerBasis= " << LayerBasis[2] << endl; cerr << s[1] << " " << s[2] << endl; /*cin >> Help_double;*/}
            if ((aurostd::abs(fractpart[1]) > _slab_epsilon_ && aurostd::abs(aurostd::abs(fractpart[1]) - 1.0) > _slab_epsilon_) ||
                (aurostd::abs(fractpart[2]) > _slab_epsilon_ && aurostd::abs(aurostd::abs(fractpart[2]) - 1.0) > _slab_epsilon_)) {
              Continue = false;
              // if(LayerBasis[1]==48 && LayerBasis[2]==49) {
              // cerr << LayerBasis[1] << " =LayerBasis= " << LayerBasis[2] << endl;
              // cerr << "WRONG-Basis:" << endl; for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][LayerBasis[i]][j] << " ";} cerr << endl;}
              // cerr << "WRONG-Plane-point:" << endl; for(j=1;j<=3;j++) {cerr << LayerSitesDirCoords[0][k][j] << " ";} ; cerr << endl;
              // cerr << s[1] << " (" << fractpart[1] << "); "  << s[2] << " (" << fractpart[2] << ") " << endl;
              // cin >> Help_string;}
            }
            k++;
          }
        } else {
          Continue = false;
        }
        if (Continue == true) {
          // cerr << "Basis:" << endl; for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][LayerBasis[i]][j] << " ";} cerr << endl;}
          Layer0BasisCosAngleCurrent = CosAngle(0, LayerBasis[1], 0, LayerBasis[2], UnitCellVector, LayerSitesDirCoords);
          // cerr << "CosAngle=" << Layer0BasisCosAngleCurrent << endl;
          if (Layer0BasisCosAngleCurrent >= 0 && Layer0BasisCosAngleCurrent < Layer0BasisCosAngle && aurostd::abs(Layer0BasisCosAngleCurrent - Layer0BasisCosAngle) > _slab_epsilon_) {
            for (i = 1; i <= 2; i++) {
              Layer0Basis[i] = LayerBasis[i];
            }
            Layer0BasisCosAngle = Layer0BasisCosAngleCurrent;
            BasisFound = true;
            // cerr << "Cos0Angle=" << Layer0BasisCosAngle << endl; //cin >> Help_double;
          }
        }
        LayerBasis[2]++;
      }
      LayerBasis[1]++;
    }
    if (BasisFound == false || 180.0 / PI * acos(Layer0BasisCosAngle) < 10.0) {
      cerr << __AFLOW_FUNC__ << " Basis in plane was not found" << endl; // return 1;
    }
    for (i = 1; i <= 2; i++) {
      AbsValueSlab_Basis[i] = slab::VectorAbsValue(0, Layer0Basis[i], UnitCellVector, LayerSitesDirCoords);
    }
    // cerr << "Primitive Basis in plane:" << endl; for(i=1;i<=2;i++) {for(j=1;j<=3;j++) { cerr << LayerSitesDirCoords[0][Layer0Basis[i]][j] << " ";} cerr << endl;}
    // cerr << "------------------------------" << endl;
    //----------------  Build "vertical" basis vector perpendicular to slab  ------------------------------------------------------//
    hkl_Length = slab::hkl_CartCoord_Length(hkl_CartCoord, UnitCellVector, hkl);
    // cerr << "hkl_CartCoord="; for(i=1;i<=3;i++) { cerr << hkl_CartCoord[i] << " ";}; cerr << endl;
    // cerr << "hkl_Length=" << hkl_Length << endl;
    double Help_double = (NumFilledLayers + NumEmptyLayers) / (hkl_Length * hkl_Length);
    for (j = 1; j <= 3; j++) {
      Slab_Basis_CartCoord[3][j] = Help_double * hkl_CartCoord[j];
    }
    AbsValueSlab_Basis[3] = 0;
    for (j = 1; j <= 3; j++) {
      AbsValueSlab_Basis[3] += Slab_Basis_CartCoord[3][j] * Slab_Basis_CartCoord[3][j];
    };
    AbsValueSlab_Basis[3] = sqrt(AbsValueSlab_Basis[3]);

    // cerr << "AbsValueSlab_Basis[3]=" << AbsValueSlab_Basis[3] << endl;

    //----------------  Build two "horisontal" slab basis vectors  ------------------------------------------------------//
    for (k = 1; k <= 2; k++) {
      for (j = 1; j <= 3; j++) {
        Slab_Basis_CartCoord(k, j) = 0;
        for (i = 1; i <= 3; i++) {
          Slab_Basis_CartCoord(k, j) += LayerSitesDirCoords[0][Layer0Basis[k]][i] * UnitCellVector(i, j);
        }
      }
    }
    for (k = 1; k <= 2; k++) {
      for (j = 1; j <= 3; j++) {
        Slab_Basis_CartCoord(k, j) *= In_Plane_Multiplication[k];
      }
    }
    //----------------  Making Slab_Basis to be right-handed: if 1[23]<0 then 1<->2 ------------------------------------------------------//
    Help_double = Slab_Basis_CartCoord(1, 1) * (Slab_Basis_CartCoord(2, 2) * Slab_Basis_CartCoord(3, 3) - Slab_Basis_CartCoord(3, 2) * Slab_Basis_CartCoord(2, 3)) -
                  Slab_Basis_CartCoord(1, 2) * (Slab_Basis_CartCoord(2, 1) * Slab_Basis_CartCoord(3, 3) - Slab_Basis_CartCoord(3, 1) * Slab_Basis_CartCoord(2, 3)) +
                  Slab_Basis_CartCoord(1, 3) * (Slab_Basis_CartCoord(2, 1) * Slab_Basis_CartCoord(3, 2) - Slab_Basis_CartCoord(3, 1) * Slab_Basis_CartCoord(2, 2));
    if (Help_double < 0) {
      for (j = 1; j <= 3; j++) {
        Help_double = Slab_Basis_CartCoord(1, j);
        Slab_Basis_CartCoord(1, j) = Slab_Basis_CartCoord(2, j);
        Slab_Basis_CartCoord(2, j) = Help_double;
      }
    }
    // for(i=1;i<=3;i++) {for(j=1;j<=3;j++) {cerr << Slab_Basis_CartCoord(i,j) << " ";} cerr << endl;}
    //----------------  Searching for sites (in filled layers) that are within the slab unit cell  ------------------------------------------------------//
    MATRIX(1, 1) = 0;
    for (j = 1; j <= 3; j++) {
      MATRIX(1, 1) += Slab_Basis_CartCoord(1, j) * Slab_Basis_CartCoord(1, j);
    }
    MATRIX(2, 2) = 0;
    for (j = 1; j <= 3; j++) {
      MATRIX(2, 2) += Slab_Basis_CartCoord(2, j) * Slab_Basis_CartCoord(2, j);
    }
    MATRIX(1, 2) = 0;
    for (j = 1; j <= 3; j++) {
      MATRIX(1, 2) += Slab_Basis_CartCoord(1, j) * Slab_Basis_CartCoord(2, j);
    }
    MATRIX(2, 1) = MATRIX(1, 2);
    DET = MATRIX(1, 1) * MATRIX(2, 2) - MATRIX(1, 2) * MATRIX(2, 1);

    // cerr << MATRIX(1,1) << " " << MATRIX(1,2) << endl << MATRIX(2,1) << " " << MATRIX(2,2) << endl << " DET= " << DET << endl;
    INVERSE_MATRIX(1, 1) = MATRIX(2, 2) / DET;
    INVERSE_MATRIX(1, 2) = -MATRIX(1, 2) / DET;
    INVERSE_MATRIX(2, 1) = -MATRIX(2, 1) / DET;
    INVERSE_MATRIX(2, 2) = MATRIX(1, 1) / DET;
    vector<vector<double>> BasisSitesInSlabDirectCoord;
    NumBasisSitesInSlab = 0;
    for (Layer = 0; Layer <= (NumFilledLayers - 1); Layer++) {
      for (k = 1; k <= NumSitesInPlane[Layer]; k++) {
        for (j = 1; j <= 3; j++) {
          SiteCartCoord[j] = 0;
          for (i = 1; i <= 3; i++) {
            SiteCartCoord[j] += LayerSitesDirCoords[Layer][k][i] * UnitCellVector(i, j);
          }
        }
        for (i = 1; i <= 2; i++) {
          RS[i] = 0;
          for (j = 1; j <= 3; j++) {
            RS[i] += SiteCartCoord[j] * Slab_Basis_CartCoord(i, j);
          }
        }
        for (i = 1; i <= 2; i++) {
          SiteDirectCoord[i] = 0;
          for (j = 1; j <= 2; j++) {
            SiteDirectCoord[i] += INVERSE_MATRIX(i, j) * RS[j];
          }
        }
        // cerr << SiteDirectCoord[1] << " " << SiteDirectCoord[2] << endl; // cin >> Help_double;
        if (SiteDirectCoord[1] >= 0 && SiteDirectCoord[1] < 1 && aurostd::abs(SiteDirectCoord[1] - 1) > _slab_epsilon_ && SiteDirectCoord[2] >= 0 && SiteDirectCoord[2] < 1 &&
            aurostd::abs(SiteDirectCoord[2] - 1) > _slab_epsilon_) {
          // SiteDirectCoord[3]=Layer/hkl_Length/AbsValueSlab_Basis[3];
          SiteDirectCoord[3] = 1.0 * Layer / (NumFilledLayers + NumEmptyLayers);
          // cerr << SiteDirectCoord[1] << " " << SiteDirectCoord[2] << " " << SiteDirectCoord[3] << endl;
          NumBasisSitesInSlab = NumBasisSitesInSlab + 1;
          BasisSitesInSlabDirectCoord.resize(NumBasisSitesInSlab + 1);
          BasisSitesInSlabDirectCoord[NumBasisSitesInSlab].resize(4);
          for (j = 1; j <= 3; j++) {
            BasisSitesInSlabDirectCoord[NumBasisSitesInSlab][j] = SiteDirectCoord[j];
          }
        }
      }
    }
    //---  Coordinates of initial slab unit cell sites (maybe with actual number of layers to be larger than necessary because only layers of one Bravais lattice were layer-enumerated) ------------//
    vector<vector<vector<double>>> ListSiteDirectCoordWRTslab;
    vector<vector<vector<double>>> ListSiteCartCoord;
    ListSiteDirectCoordWRTslab.resize(NumberElements + 1);
    ListSiteCartCoord.resize(NumberElements + 1);
    for (Element = 1; Element <= NumberElements; Element++) {
      for (ElementSite = 1; ElementSite <= NumAtomsForElementUC[Element]; ElementSite++) {
        for (l[1] = SearchMax; l[1] >= -SearchMax; l[1]--) {
          for (l[2] = SearchMax; l[2] >= -SearchMax; l[2]--) {
            for (l[3] = SearchMax; l[3] >= -SearchMax; l[3]--) {
              for (j = 1; j <= 3; j++) {
                SiteCartCoord[j] = 0;
                for (i = 1; i <= 3; i++) {
                  SiteCartCoord[j] += (AtomDirCoords[Element][ElementSite][i] + l[i]) * UnitCellVector(i, j);
                }
              }
              for (i = 1; i <= 3; i++) {
                RS[i] = 0;
                for (j = 1; j <= 3; j++) {
                  RS[i] += SiteCartCoord[j] * Slab_Basis_CartCoord(i, j);
                }
              }
              // cerr << Element << endl;
              // for(i=1;i<=3;i++) { cerr << l[i] << " ";}; cerr << " l" << endl;
              // for(i=1;i<=3;i++) { cerr << AtomDirCoords[Element][ElementSite][i]+l[i] << " ";}; cerr << " AtomDirCoords" << endl;
              // for(i=1;i<=3;i++) { cerr << SiteCartCoord[i] << " ";}; cerr << " SiteCartCoord" << endl;

              /*
                 for(i=1;i<=3;i++) { cerr << RS[i] << " ";}; cerr << " RS" << endl;
                 cerr << MATRIX[1][1] << " " << MATRIX[1][2] << endl;
                 cerr << MATRIX[2][1] << " " << MATRIX[2][2] << endl;
                 cerr  << " DET= " << DET << endl;
                 cerr << INVERSE_MATRIX[1][1] << " " << INVERSE_MATRIX[1][2] << endl;
                 cerr << INVERSE_MATRIX[2][1] << " " << INVERSE_MATRIX[2][2] << endl;
                 */

              for (i = 1; i <= 2; i++) {
                SiteDirectCoordWRTslab[i] = 0;
                for (j = 1; j <= 2; j++) {
                  SiteDirectCoordWRTslab[i] += INVERSE_MATRIX(i, j) * RS[j];
                }
              }
              SiteDirectCoordWRTslab[3] = RS[3] / (AbsValueSlab_Basis[3] * AbsValueSlab_Basis[3]);
              for (i = 1; i <= 3; i++) {
                if (aurostd::abs(SiteDirectCoordWRTslab[i] - 0.0) < _slab_epsilon_) {
                  SiteDirectCoordWRTslab[i] = 0.0;
                }
              }
              for (i = 1; i <= 3; i++) {
                if (aurostd::abs(SiteDirectCoordWRTslab[i] - 1.0) < _slab_epsilon_) {
                  SiteDirectCoordWRTslab[i] = 1.0;
                }
              }
              // for(i=1;i<=3;i++) { cerr << SiteDirectCoordWRTslab[i] << " ";}; cerr << " SiteDirectCoordWRTslab" << endl;
              // if(SiteCartCoord[1]==0.5 && Element==2) {cin >> Help_double;}

              if (SiteDirectCoordWRTslab[1] >= 0.0 && SiteDirectCoordWRTslab[1] < 1.0 && SiteDirectCoordWRTslab[2] >= 0.0 && SiteDirectCoordWRTslab[2] < 1.0 && SiteDirectCoordWRTslab[3] >= 0.0 &&
                  SiteDirectCoordWRTslab[3] < 1.0) {
                NumSites[Element]++;
                ListSiteDirectCoordWRTslab[Element].resize(NumSites[Element] + 1);
                ListSiteCartCoord[Element].resize(NumSites[Element] + 1);
                ListSiteDirectCoordWRTslab[Element][NumSites[Element]].resize(4);
                ListSiteCartCoord[Element][NumSites[Element]].resize(4);
                for (i = 1; i <= 3; i++) {
                  ListSiteDirectCoordWRTslab[Element][NumSites[Element]][i] = SiteDirectCoordWRTslab[i];
                }
                for (i = 1; i <= 3; i++) {
                  ListSiteCartCoord[Element][NumSites[Element]][i] = SiteCartCoord[i];
                }
                // for(i=1;i<=3;i++) { cerr << SiteDirectCoordWRTslab[i] << " ";} cerr << "  Direct" << Element << endl;
                // for(i=1;i<=3;i++) { cerr << SiteCartCoord[i] << " ";} cerr << "  Cart" << Element << endl;
              }
            }
          }
        }
      }
    }
    //-----------  SORTING atoms inside slab u.c. by Layers -------------------------------------//
    vector<int> NumInLayer(1, 0);
    vector<double> LayerDirectCoordWRTslab3(1, 0);
    vector<vector<vector<int>>> AtomInLayer;
    NumLayers = 0;
    for (k = 1; k <= NumberElements; k++) {
      for (i = 1; i <= NumSites[k]; i++) {
        Layer = 1;
        Continue = true;
        while (Layer <= NumLayers && Continue == true) {
          if (aurostd::abs(ListSiteDirectCoordWRTslab[k][i][3] - LayerDirectCoordWRTslab3[Layer]) < _slab_epsilon_) {
            NumInLayer[Layer]++;
            AtomInLayer[Layer].resize(NumInLayer[Layer] + 1);
            AtomInLayer[Layer][NumInLayer[Layer]].resize(3);
            AtomInLayer[Layer][NumInLayer[Layer]][1] = k;
            AtomInLayer[Layer][NumInLayer[Layer]][2] = i;
            Continue = false;
          }
          Layer = Layer + 1;
        }
        if (Continue == true) {
          NumLayers++;
          LayerDirectCoordWRTslab3.push_back(ListSiteDirectCoordWRTslab[k][i][3]);
          NumInLayer.push_back(1);
          AtomInLayer.resize(NumLayers + 1);
          AtomInLayer[NumLayers].resize(2);
          AtomInLayer[NumLayers][1].resize(3);
          AtomInLayer[NumLayers][1][1] = k;
          AtomInLayer[NumLayers][1][2] = i;
        }
      }
    }
    /*
       for(Layer=1; Layer<=NumLayers; Layer++) {
       cerr << "------ " << LayerDirectCoordWRTslab3[Layer] << endl;
       for(k=1;k<=NumInLayer[Layer];k++) {
       cerr << AtomInLayer[Layer][k][1] << "  " << AtomInLayer[Layer][k][2] << endl;
       for(i=1;i<=3;i++) {
       cerr << ListSiteDirectCoordWRTslab[AtomInLayer[Layer][k][1]][AtomInLayer[Layer][k][2]][i] << " ";
       }      cerr  << endl;      }  }
       cerr << "2--------------------" << endl;
       */
    //-----------  Order Layers by third coordinate wrt S3 ------------------------------------//
    vector<int> OrderLayers(NumLayers + 1);
    vector<double> LayerCoord3(NumLayers + 1);
    LayerCoord3[0] = -100;
    for (Layer = 1; Layer <= NumLayers; Layer++) {
      LastCandidateCoord3 = 1000;
      for (k = 1; k <= NumLayers; k++) {
        CurrCandidateCoord3 = LayerDirectCoordWRTslab3[k];
        if (CurrCandidateCoord3 > LayerCoord3[Layer - 1] && CurrCandidateCoord3 < LastCandidateCoord3) {
          LastCandidateCoord3 = CurrCandidateCoord3;
          OrderLayers[Layer] = k;
        }
      }
      LayerCoord3[Layer] = LayerDirectCoordWRTslab3[OrderLayers[Layer]];
    }

    // for(Layer=1; Layer<=NumLayers; Layer++) {cerr << OrderLayers[Layer]<< "  S3=" << LayerDirectCoordWRTslab3[OrderLayers[Layer]]  << endl;}
    // cerr << "--------------------" << endl;

    //-----------  Shorten S3 (third, "normal to layers" basis vector of slab uc) and project sites coordinates wrt new S3 ------------------------------------//

    if ((NumFilledLayers + NumEmptyLayers) < NumLayers) {
      for (j = 1; j <= 3; j++) {
        Slab_Basis_CartCoord[3][j] *= LayerDirectCoordWRTslab3[OrderLayers[NumFilledLayers + NumEmptyLayers + 1]];
      }

      for (k = 1; k <= NumberElements; k++) {
        for (i = 1; i <= NumSites[k]; i++) {
          ListSiteDirectCoordWRTslab[k][i][3] /= LayerDirectCoordWRTslab3[OrderLayers[NumFilledLayers + NumEmptyLayers + 1]];
        }
      }
    }

    // Adjust the number of atoms in u.c. for true number of layers
    // for(k=1;k<=NumberElements;k++) {NumSites[k]=0;}
    // for(k=1;k<=NumberElements;k++) {
    //   for(Layer=NumFilledLayers; Layer>=1; Layer--) {
    //   for(i=1;i<=NumInLayer[OrderLayers[Layer]];i++) {
    //   if(AtomInLayer[OrderLayers[Layer]][i][1]==k)
    //   {NumSites[k]++;}
    // } } }
    for (k = 1; k <= NumberElements; k++) {
      NumSites[k] = 0;
    } // here
    for (k = 1; k <= NumberElements; k++) {
      for (Layer = NumFilledLayers; Layer >= 1; Layer--) {
        for (i = 1; i <= NumInLayer[OrderLayers[Layer]]; i++) {
          if (AtomInLayer[OrderLayers[Layer]][i][1] == k) {
            NumSites[k]++;
          }
        }
      }
    }

    //-----------  FILE OUTPUT ------------------------------------//

    // Output xstructure for slab
    xstructure str_out("");

    str_out.is_vasp4_poscar_format = str_in.is_vasp4_poscar_format;
    str_out.is_vasp5_poscar_format = str_in.is_vasp5_poscar_format;

    stringstream sout;

    sout << str_in.title << ", Slab " << hkl(1) << hkl(2) << hkl(3) << ": (" << NumFilledLayers << " full + " << NumEmptyLayers << " empty)*(" << In_Plane_Multiplication[1] << "x" << In_Plane_Multiplication[2]
         << " InPlaneMultipl)";
    str_out.title = sout.str();
    str_out.scale = LattPar;

    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        str_out.lattice(i, j) = Slab_Basis_CartCoord(i, j);
      }
    }
    str_out.FixLattices(); // CO20180202

    // str_out.num_each_type.clear();
    // str_out.comp_each_type.clear();

    _atom newatom;
    iatom = 0;
    for (k = 1; k <= NumberElements; k++) {
      for (Layer = NumFilledLayers; Layer >= 1; Layer--) {
        for (i = 1; i <= NumInLayer[OrderLayers[Layer]]; i++) {
          if (AtomInLayer[OrderLayers[Layer]][i][1] == k) {
            iatom++;
            newatom.type = k - 1; // CO20180724 - assign type needed for AddAtom()
            specie = str_in.SpeciesLabel(k - 1); // CO20180724 - this should be NON-empty as per assigning_fake_names above
            if (specie.empty() || aurostd::substring2bool(specie, "name not given")) {
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "ERROR! Species label not found!", _FILE_WRONG_FORMAT_); // CO20180724
            }
            newatom.name = specie; // CO20180724 - assign name needed for AddAtom()
            newatom.name_is_given = true; // CO20180724 - assign name needed for AddAtom()
            for (j = 1; j <= 3; j++) {
              newatom.fpos[j] = ListSiteDirectCoordWRTslab[k][AtomInLayer[OrderLayers[Layer]][i][2]][j];
            }
            newatom.cpos = F2C(str_out.scale, str_out.lattice, newatom.fpos); // CO20180202
            str_out.AddAtom(newatom, false); // CO20230319 - add by type
          }
        }
      }
    }
    if (assigning_fake_names) {
      for (size_t i = 0; i < str_out.atoms.size(); i++) {
        str_out.atoms[i].name_is_given = false;
      }
    } // CO20180724 - since these are fake names, don't print out them out

    // CO20180202
    if (LDEBUG) {
      cerr << "PRINTING OUT STRUCTURE ATTRIBUTES" << endl;
      cerr << "str_out.atoms.size()=" << str_out.atoms.size() << endl;
      cerr << "str_out.num_each_type.size()=" << str_out.num_each_type.size() << endl;
      cerr << "str_out.comp_each_type.size()=" << str_out.comp_each_type.size() << endl;
      for (size_t i = 0; i < str_out.num_each_type.size(); i++) {
        cerr << "str_out.num_each_type[i]=" << str_out.num_each_type[i] << endl;
      }
      cerr << str_out << endl;
    }

    ////////////////////////////////////////////////////////////
    // CO+DU20180705 START
    // everything here (below AddAtom()) is a HACK and needs to be
    // fixed
    // AddAtom() must handle names/num_each_type/comp_each_type
    // otherwise we get uneven species parameter vectors
    //(e.g., specices_pp) and hence CORE files
    // to fix, we need to perform full copies of atoms inside
    // xstructure (newatom=str_in.atoms[XX]), otherwise
    // we lose type + name which is absolutely critical for AddAtom()

    // CO20180202 - this is obsolete, it is done INSIDE AddAtom()
    // DU20180705 - putting back as it doesn't work without it
    //  we need to add type + name before AddAtom()
    //  this is a temporary patch, fix later
    // CO+DU20180705 STOP
    ////////////////////////////////////////////////////////////

    // for(i=0;i<NumberElements;i++)
    // str_out.num_each_type[i]=(NumSites[i+1]);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return str_out;
  }
} // namespace slab

// CO20190601 START
namespace slab {
  //[CO20190520 - this is wrong, do not convert to real space]#define HKL_DUAL_TEST 0 //CO20190520 - this is WRONG, so keep 0: do NOT convert to real space
  aurostd::xvector<double> HKLPlane2Normal(const xstructure& xstr_in, int h, int k, int l) {
    return HKLPlane2Normal(xstr_in.scale * xstr_in.lattice, h, k, l);
  } // CO20190320
  aurostd::xvector<double> HKLPlane2Normal(const aurostd::xmatrix<double>& lattice, int h, int k, int l) {
    const xvector<int> hkl;
    hkl[1] = h;
    hkl[2] = k;
    hkl[3] = l;
    return HKLPlane2Normal(lattice, hkl);
  } // CO20190320
  aurostd::xvector<double> HKLPlane2Normal(const xstructure& xstr_in, const aurostd::xvector<int>& hkl) {
    return HKLPlane2Normal(xstr_in.scale * xstr_in.lattice, hkl);
  } // CO20190320
  aurostd::xvector<double> HKLPlane2Normal(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);

    // http://www.mse.mtu.edu/~drjohn/my3200/stereo/sg5.html
    // use metric tensor of reciprocal lattice to write
    // also here: http://ssd.phys.strath.ac.uk/resources/crystallography/crystallographic-direction-calculator/
    // hkl is nothing more than fractional coordinates in reciprocal space
    // useful relationship: kM=(2*PI)^2*inverse(M)
    // https://it.iucr.org/Ba/ch1o1v0001/ - metric tensors of the covariant (direct) and contravariant (reciprocal) bases
    // http://physastro-msci.tripod.com/webonmediacontents/notes1.pdf
    const aurostd::xvector<double> dhkl = aurostd::xvector2utype<double>(hkl); // need double for operations
    const aurostd::xmatrix<double> klattice = ReciprocalLattice(lattice);
    const aurostd::xmatrix<double> kf2c = trasp(klattice); // convert fractional to cartesian
    //[CO20190520 - this is wrong, do not convert to real space]#if !HKL_DUAL_TEST
    aurostd::xvector<double> n = kf2c * dhkl; // h*b1+k*b2+l*b3
    n /= aurostd::modulus(n); // normalize
    //[CO20190520 - this is wrong, do not convert to real space]#else
    //[CO20190520 - this is wrong, do not convert to real space]  //[CO20190528 - do not convert vector from reciprocal to direct, kn is the direction]
    //[CO20190520 - this is wrong, do not convert to real space]  aurostd::xvector<double> kn=kf2c*dhkl;               //h*b1+k*b2+l*b3
    //[CO20190520 - this is wrong, do not convert to real space]  aurostd::xmatrix<double> M=MetricTensor(lattice);    //metric tensor of direct space
    //[CO20190520 - this is wrong, do not convert to real space]  aurostd::xvector<double> n=M*kn;                     //convert from reciprocal (contravariant) to direct (covariant): direct = metric(direct) *
    // reciprocal [CO20190520 - this is wrong, do not convert to real space]  n/=aurostd::modulus(n);                     //normalize [CO20190520 - this is wrong, do not convert to real space]#endif

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl=" << hkl << endl;
      cerr << __AFLOW_FUNC__ << " lattice=" << endl;
      cerr << lattice << endl;
      cerr << __AFLOW_FUNC__ << " klattice=" << endl;
      cerr << klattice << endl;
      //[CO20190520 - this is wrong, do not convert to real space]#if HKL_DUAL_TEST
      //[CO20190520 - this is wrong, do not convert to real space]    //[CO20190528 - do not convert vector from reciprocal to direct, kn is the direction]
      //[CO20190520 - this is wrong, do not convert to real space]    cerr << __AFLOW_FUNC__ << " M=" << endl;cerr << M << endl;
      //[CO20190520 - this is wrong, do not convert to real space]    cerr << __AFLOW_FUNC__ << " kn=" << kn/aurostd::modulus(kn) << endl;
      //[CO20190520 - this is wrong, do not convert to real space]#endif
      cerr << __AFLOW_FUNC__ << " n=" << n << endl;
    }

    return n;
  }
  bool Normal2HKLPlane(const xstructure& xstr_in, const aurostd::xvector<double>& n, aurostd::xvector<int>& hkl) {
    return Normal2HKLPlane(xstr_in.scale * xstr_in.lattice, n, hkl);
  } // CO20190320
  bool Normal2HKLPlane(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<double>& n, aurostd::xvector<int>& hkl) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);

    // http://www.mse.mtu.edu/~drjohn/my3200/stereo/sg5.html
    // use metric tensor of reciprocal lattice to write
    // also here: http://ssd.phys.strath.ac.uk/resources/crystallography/crystallographic-direction-calculator/
    // useful relationship: kM=(2*PI)^2*inverse(M)
    // https://it.iucr.org/Ba/ch1o1v0001/ - metric tensors of the covariant (direct) and contravariant (reciprocal) bases
    // http://physastro-msci.tripod.com/webonmediacontents/notes1.pdf
    // https://physcourses.lums.edu.pk/wp-content/uploads/2012/09/Reciprocal-lattices.pdf
    const aurostd::xmatrix<double> klattice = ReciprocalLattice(lattice);
    const aurostd::xmatrix<double> kc2f = inverse(trasp(klattice)); // convert cartesian to fractional
    //[CO20190520 - this is wrong, do not convert to real space]#if !HKL_DUAL_TEST
    aurostd::xvector<double> dhkl_frac = kc2f * n; // hkl (double) in fractional form (need to convert to integers)
    //[CO20190520 - this is wrong, do not convert to real space]#else
    //[CO20190520 - this is wrong, do not convert to real space]  //[CO20190528 - do not convert vector from reciprocal to direct, kn is the direction]
    //[CO20190520 - this is wrong, do not convert to real space]  aurostd::xmatrix<double> kM=MetricTensor(klattice);        //metric tensor of reciprocal space
    //[CO20190520 - this is wrong, do not convert to real space]  aurostd::xvector<double> kn=kM*n;                          //convert from direct (covariant) to reciprocal (contravariant): reciprocal =
    // metric(reciprocal) * direct [CO20190520 - this is wrong, do not convert to real space]  aurostd::xvector<double> dhkl_frac=kc2f*kn;                //hkl (double) in fractional form (need to convert to
    // integers) [CO20190520 - this is wrong, do not convert to real space]#endif
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " dhkl_frac=" << dhkl_frac << endl;
    }

    // define tolerance based on how far we explore (up to max_multiple in hkl)
    const uint max_multiple = 1e4;
    const double zero_tol = pow(10, (int) ceil(log10(1.0 / max_multiple)));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " max_multiple=" << max_multiple << endl;
      cerr << __AFLOW_FUNC__ << " zero_tol=" << zero_tol << endl;
    }
    // find minimum of hkl, be careful of 0s
    // dhkl_frac/=aurostd::min(dhkl_frac); //what if it's 0, we need to be more careful
    double min = AUROSTD_MAX_DOUBLE;
    for (int i = dhkl_frac.lrows; i <= dhkl_frac.urows; i++) {
      if (std::abs(dhkl_frac[i]) > zero_tol && std::abs(dhkl_frac[i]) < std::abs(min)) {
        min = dhkl_frac[i];
      }
    }
    if (min == AUROSTD_MAX_DOUBLE) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Could not find minimum value of dhkl_frac", _VALUE_ERROR_);
    }
    dhkl_frac /= min;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " dhkl_frac=" << dhkl_frac << endl;
    }

    // explore multiples of hkl up to tolerance allows
    aurostd::xvector<double> dhkl = dhkl_frac;
    bool found = false;
    for (uint i = 1; i <= max_multiple && !found; i++) { // BRUTE (stupid) force, there's probably an algorithm out there for it... //start with 1, we may already have the solution
      dhkl = dhkl_frac * (double) i;
      // if(LDEBUG) {cerr << __AFLOW_FUNC__ << " dhkl_frac*" << i << "=" << dhkl << endl;}
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " dhkl_frac*" << i << "=" << setprecision(15) << dhkl[1] << " " << dhkl[2] << " " << dhkl[3] << endl;
      }
      if (aurostd::isinteger(dhkl, zero_tol)) {
        found = true;
        break;
      }
    }
    if (!found) {
      // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not find valid hkl",_VALUE_ERROR_);
      return false;
    }

    // convert dhkl to hkl (integer)
    for (int i = dhkl.lrows; i <= dhkl.urows; i++) {
      hkl[i] = (int) nint(dhkl[i]);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n=" << n << endl;
      //[CO20190520 - this is wrong, do not convert to real space]#if HKL_DUAL_TEST
      //[CO20190520 - this is wrong, do not convert to real space]    //[CO20190528 - do not convert vector from reciprocal to direct, kn is the direction]
      //[CO20190520 - this is wrong, do not convert to real space]    cerr << __AFLOW_FUNC__ << " kM=" << endl;cerr << kM << endl;
      //[CO20190520 - this is wrong, do not convert to real space]#endif
      cerr << __AFLOW_FUNC__ << " hkl=" << hkl << endl;
    }

    return true;
  }

  vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const xstructure& xstr_in, int h, int k, int l) {
    return getHKLPlaneIntercepts(xstr_in.scale * xstr_in.lattice, h, k, l);
  } // CO20190320
  vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const aurostd::xmatrix<double>& lattice, int h, int k, int l) {
    const xvector<int> hkl;
    hkl[1] = h;
    hkl[2] = k;
    hkl[3] = l;
    return getHKLPlaneIntercepts(lattice, hkl);
  } // CO20190320
  vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const xstructure& xstr_in, const aurostd::xvector<int>& hkl) {
    return getHKLPlaneIntercepts(xstr_in.scale * xstr_in.lattice, hkl);
  } // CO20190320
  vector<aurostd::xvector<double>> getHKLPlaneIntercepts(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;

    if (hkl[1] == 0 && hkl[2] == 0 && hkl[3] == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "hkl=(0,0,0)", _INPUT_ERROR_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl=" << hkl << endl;
    }

    const aurostd::xmatrix<double> f2c = trasp(lattice);
    const aurostd::xmatrix<double> c2f = inverse(f2c);

    // http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/G_is_orthogonal_to_hkl_plane.html
    const aurostd::xvector<double>& a1 = lattice(1);
    const aurostd::xvector<double>& a2 = lattice(2);
    const aurostd::xvector<double>& a3 = lattice(3);
    const int multiple = aurostd::LCM(hkl); // 1.0 //See W. Sun and G. Ceder, Surface Science 617 (2013) 53-59
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " multiple=" << multiple << endl;
    }

    vector<aurostd::xvector<double>> intercepts;

    // get zeros
    vector<uint> zero_indices;
    int count_zeros = 0;
    for (int i = hkl.lrows; i <= hkl.urows; i++) {
      if (hkl[i] == 0) {
        zero_indices.push_back(i);
        count_zeros++;
      }
    }

    // See W. Sun and G. Ceder, Surface Science 617 (2013) 53-59
    if (count_zeros == 0) {
      intercepts.push_back(((double) multiple / (double) hkl[1]) * a1);
      intercepts.push_back(((double) multiple / (double) hkl[2]) * a2);
      intercepts.push_back(((double) multiple / (double) hkl[3]) * a3);
    } else if (count_zeros == 1) {
      if (aurostd::WithinList(zero_indices, (uint) 1)) {
        intercepts.push_back(((double) multiple / (double) hkl[2]) * a2);
        intercepts.push_back(((double) multiple / (double) hkl[3]) * a3);
        // intercepts.push_back( intercepts[0] + a1 ); //consistent order, but it really doesn't matter
        const aurostd::xvector<double> tmp = (intercepts[0] + a1); // consistent order, but it really doesn't matter
        intercepts.insert(intercepts.begin(), tmp); // consistent order, but it really doesn't matter
      } else if (aurostd::WithinList(zero_indices, (uint) 2)) {
        intercepts.push_back(((double) multiple / (double) hkl[1]) * a1);
        intercepts.push_back(intercepts[0] + a2);
        intercepts.push_back(((double) multiple / (double) hkl[3]) * a3);
      } else { // aurostd::WithinList(zero_indices,(uint)3)
        intercepts.push_back(((double) multiple / (double) hkl[1]) * a1);
        intercepts.push_back(((double) multiple / (double) hkl[2]) * a2);
        intercepts.push_back(intercepts[0] + a3);
      }
    } else { // count_zeros==2
      const aurostd::xvector<double> tmp; // 0,0,0
      if (!aurostd::WithinList(zero_indices, (uint) 1)) {
        intercepts.push_back(tmp);
        intercepts.push_back(a2);
        intercepts.push_back(a3);
      } else if (!aurostd::WithinList(zero_indices, (uint) 2)) {
        intercepts.push_back(a1);
        intercepts.push_back(tmp);
        intercepts.push_back(a3);
      } else { //! aurostd::WithinList(zero_indices,(uint)3)
        intercepts.push_back(a1);
        intercepts.push_back(a2);
        intercepts.push_back(tmp);
      }
    }

    if (LDEBUG) {
      for (size_t i = 0; i < intercepts.size(); i++) {
        cerr << __AFLOW_FUNC__ << " intercepts[" << i << "]=" << intercepts[i] << endl;
      }
    }

    // check that vectors created with intercepts are Bravais lattice vectors
    const aurostd::xvector<double> v1 = intercepts[1] - intercepts[0];
    const aurostd::xvector<double> v2 = intercepts[2] - intercepts[0];
    const aurostd::xvector<double> v3 = intercepts[2] - intercepts[1];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " v1=" << v1 << endl;
      cerr << __AFLOW_FUNC__ << " v2=" << v2 << endl;
      cerr << __AFLOW_FUNC__ << " v3=" << v3 << endl;
    }
    aurostd::xvector<double> fpos;
    fpos = c2f * v1;
    if (!aurostd::isinteger(fpos)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v1 is not a Bravais lattice", _INPUT_ERROR_);
    }
    fpos = c2f * v2;
    if (!aurostd::isinteger(fpos)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v2 is not a Bravais lattice", _INPUT_ERROR_);
    }
    fpos = c2f * v3;
    if (!aurostd::isinteger(fpos)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v3 is not a Bravais lattice", _INPUT_ERROR_);
    }

    // check that normal is orthogonal with vectors created
    const aurostd::xvector<double> n = HKLPlane2Normal(lattice, hkl);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n=" << n << endl;
    }
    if (!aurostd::isequal(aurostd::scalar_product(n, v1), 0.0, _ZERO_TOL_)) {
      message << "n[" << n << "] is not orthogonal to v1[" << v1 << "]: scalar_product=" << aurostd::scalar_product(n, v1) << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    if (!aurostd::isequal(aurostd::scalar_product(n, v2), 0.0, _ZERO_TOL_)) {
      message << "n[" << n << "] is not orthogonal to v2[" << v2 << "]: scalar_product=" << aurostd::scalar_product(n, v2) << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    if (!aurostd::isequal(aurostd::scalar_product(n, v3), 0.0, _ZERO_TOL_)) {
      message << "n[" << n << "] is not orthogonal to v3[" << v3 << "]: scalar_product=" << aurostd::scalar_product(n, v3) << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }

    return intercepts;
  }

  double getSpacingHKLPlane(const xstructure& xstr_in, int h, int k, int l) {
    return getSpacingHKLPlane(xstr_in.scale * xstr_in.lattice, h, k, l);
  } // CO20190320
  double getSpacingHKLPlane(const aurostd::xmatrix<double>& lattice, int h, int k, int l) {
    const xvector<int> hkl;
    hkl[1] = h;
    hkl[2] = k;
    hkl[3] = l;
    return getSpacingHKLPlane(lattice, hkl);
  } // CO20190320
  double getSpacingHKLPlane(const xstructure& xstr_in, const aurostd::xvector<int>& hkl) {
    return getSpacingHKLPlane(xstr_in.scale * xstr_in.lattice, hkl);
  } // CO20190320
  double getSpacingHKLPlane(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);

    // http://lafactoria.lec.csic.es/mcc/attachments/article/12/Introduction%20to%20Reciprocal%20Space.pdf
    // https://web.stanford.edu/group/glam/xlab/MatSci162_172/LectureNotes/02_Geometry,%20RecLattice.pdf
    // useful relationship: kM=(2*PI)^2*inverse(M)
    const aurostd::xvector<double> dhkl = aurostd::xvector2utype<double>(hkl); // need double for operations
    const aurostd::xmatrix<double> klattice = ReciprocalLattice(lattice);
    const aurostd::xmatrix<double> kM = MetricTensor(klattice);
    const double d_spacing = 2.0 * PI / sqrt(aurostd::scalar_product(dhkl, kM * dhkl)); // 2*pi factor here is very important (counters the one in ReciprocalLattice())

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl=" << hkl << endl;
      cerr << __AFLOW_FUNC__ << " kM=" << endl;
      cerr << kM << endl;
      cerr << __AFLOW_FUNC__ << " d_spacing=" << d_spacing << endl;
    }

    return d_spacing;
  }
  double getAngleHKLPlanes(const xstructure& xstr_in, int h1, int k1, int l1, int h2, int k2, int l2) {
    return getAngleHKLPlanes(xstr_in.scale * xstr_in.lattice, h1, k1, l1, h2, k2, l2);
  } // CO20190320
  double getAngleHKLPlanes(const aurostd::xmatrix<double>& lattice, int h1, int k1, int l1, int h2, int k2, int l2) {
    const xvector<int> hkl1;
    const xvector<int> hkl2;
    hkl1[1] = h1;
    hkl1[2] = k1;
    hkl1[3] = l1;
    hkl2[1] = h2;
    hkl2[2] = k2;
    hkl2[3] = l2;
    return getAngleHKLPlanes(lattice, hkl1, hkl2);
  } // CO20190320
  double getAngleHKLPlanes(const xstructure& xstr_in, const aurostd::xvector<int>& hkl1, const aurostd::xvector<int>& hkl2) {
    return getAngleHKLPlanes(xstr_in.scale * xstr_in.lattice, hkl1, hkl2);
  } // CO20190320
  double getAngleHKLPlanes(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<int>& hkl1, const aurostd::xvector<int>& hkl2) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);

    // http://lafactoria.lec.csic.es/mcc/attachments/article/12/Introduction%20to%20Reciprocal%20Space.pdf
    // https://web.stanford.edu/group/glam/xlab/MatSci162_172/LectureNotes/02_Geometry,%20RecLattice.pdf
    // useful relationship: kM=(2*PI)^2*inverse(M)
    const aurostd::xvector<double> dhkl1 = aurostd::xvector2utype<double>(hkl1); // need double for operations
    const aurostd::xvector<double> dhkl2 = aurostd::xvector2utype<double>(hkl2); // need double for operations
    const aurostd::xmatrix<double> klattice = ReciprocalLattice(lattice);
    const aurostd::xmatrix<double> kf2c = trasp(klattice); // convert fractional to cartesian
    const aurostd::xvector<double> n1 = kf2c * dhkl1; // h*b1+k*b2+l*b3
    const aurostd::xvector<double> n2 = kf2c * dhkl2; // h*b1+k*b2+l*b3
    const aurostd::xmatrix<double> kM = MetricTensor(klattice);
    const double angle = std::acos(aurostd::scalar_product(dhkl1, kM * dhkl2) / ((2.0 * PI) * (2.0 * PI) * aurostd::modulus(n1) * aurostd::modulus(n2))); //(2*pi)^2 factor here is very important (counters ReciprocalLattice())

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl1=" << hkl1 << endl;
      cerr << __AFLOW_FUNC__ << " hkl2=" << hkl2 << endl;
      cerr << __AFLOW_FUNC__ << " n1=" << n1 << endl;
      cerr << __AFLOW_FUNC__ << " n2=" << n2 << endl;
      cerr << __AFLOW_FUNC__ << " kM=" << endl;
      cerr << kM << endl;
      cerr << __AFLOW_FUNC__ << " angle=" << angle << endl;
    }

    return angle;
  }

  void BringInBoundary(aurostd::xvector<double>& vec, double padding) { // different than BringInCell()
    for (int i = vec.lrows; i <= vec.urows; i++) { // specialized BringInCell() for our purposes, no preference for origin (1.0 wall is fine)
      while (vec[i] > 1.0 + padding) {
        vec[i] -= 1.0;
      } // BringInCell() has preference for origin, but we don't here, so no >=
      while (vec[i] < -padding) {
        vec[i] += 1.0;
      }
    }
  }

  //[CO20190520 - plugged into BringInBoundary() with no padding]void Bring2OppositeBoundary(aurostd::xvector<double>& vec){  //bounces position to the opposite boundary
  //[CO20190520 - plugged into BringInBoundary() with no padding]  for(int i=vec.lrows;i<=vec.urows;i++){
  //[CO20190520 - plugged into BringInBoundary() with no padding]    //if(aurostd::isequal(vec[i],1.0) || vec[i]>1.0){vec[i]-=1.0;}      //_ZERO_TOL_ is too strict
  //[CO20190520 - plugged into BringInBoundary() with no padding]    //else if(aurostd::isequal(vec[i],0.0) || vec[i]<0.0){vec[i]+=1.0;} //_ZERO_TOL_ is too strict
  //[CO20190520 - plugged into BringInBoundary() with no padding]    if(vec[i]>1.0){vec[i]-=1.0;}      //_ZERO_TOL_ is too strict
  //[CO20190520 - plugged into BringInBoundary() with no padding]    else if(vec[i]<0.0){vec[i]+=1.0;} //_ZERO_TOL_ is too strict
  //[CO20190520 - plugged into BringInBoundary() with no padding]  }
  //[CO20190520 - plugged into BringInBoundary() with no padding]}

#define LOOP_ITERATION_MAX 1e3
  aurostd::xvector<double> getNextAtomInPath(const xstructure& xstr_in, const aurostd::xvector<double>& _l_cpos, const aurostd::xvector<double>& cpos_starting, vector<uint>& atoms2skip, uint& loop_iteration, bool outside_current_cell) {
    // this function takes inputs direction l_cpos and cpos_starting and finds next atom (cpos_final) in that direction
    // good for hkl calculations
    // it also returns loop_iteration (how many times we loop unit cell in that direction)
    const bool LDEBUG = (false || XHOST.DEBUG);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " starting" << endl;
    }

    // assume NO changes in xstr_in (too slow)
    const aurostd::xmatrix<double>& lattice = xstr_in.lattice;
    const aurostd::xmatrix<double>& f2c = xstr_in.f2c;
    const aurostd::xmatrix<double>& c2f = xstr_in.c2f;
    const deque<_atom>& atoms = xstr_in.atoms;
    double min_dist = xstr_in.dist_nn_min;
    if (min_dist == AUROSTD_NAN) {
      min_dist = SYM::minimumDistance(xstr_in);
    }
    double sym_eps = xstr_in.sym_eps;
    if (sym_eps == AUROSTD_NAN) {
      sym_eps = SYM::defaultTolerance(xstr_in);
    }
    const bool skew = SYM::isLatticeSkewed(lattice, min_dist, sym_eps);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " min_dist=" << min_dist << endl;
      cerr << __AFLOW_FUNC__ << " sym_eps=" << sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }

    const aurostd::xvector<double> l_cpos = _l_cpos / aurostd::modulus(_l_cpos); // ensure it is normal!
    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> l_fpos=c2f*l_cpos;l_fpos/=aurostd::modulus(l_fpos);
    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> l_cpos=c2f*l_fpos;l_cpos/=aurostd::modulus(l_cpos);
    const uint loop_iteration_starting = loop_iteration;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " l_cpos=" << l_cpos << endl;
      //[CO20190520 - safe-guard for working in fractional space (skew)]cerr << __AFLOW_FUNC__ << " l_fpos=" << l_fpos << endl;
      cerr << __AFLOW_FUNC__ << " loop_iteration=" << loop_iteration << endl;
    }

    // define unit box (in fpos FIRST, then convert to cpos)
    // 2 points
    aurostd::xvector<double> p_origin;
    aurostd::xvector<double> p_top; // p_origin [0,0,0]
    p_top[1] = p_top[2] = p_top[3] = 1.0;
    p_origin = f2c * p_origin; // convert to cpos
    p_top = f2c * p_top; // convert to cpos
    // 6 normals all pointing outward
    vector<aurostd::xvector<double>> v_n_fpos;
    vector<aurostd::xvector<double>> v_n_cpos;
    vector<aurostd::xvector<double>> v_line_plane_intersection_cpos; // v_line_plane_intersection_fpos
    // xvector<double> line_plane_intersection_cpos;
    vector<double> v_d;
    vector<bool> v_line_plane_intersect;
    for (uint i = 0; i < 6; i++) {
      v_n_fpos.emplace_back();
      v_n_cpos.emplace_back();
      v_line_plane_intersection_cpos.emplace_back(); // v_line_plane_intersection_fpos.push_back(aurostd::xvector<double>());
      v_d.push_back(0.0);
      v_line_plane_intersect.push_back(false);
    }
    // borrow instead from LatticeDimensionSphere()
    const aurostd::xmatrix<double> normals;
    for (int m = 1; m <= 3; m++) {
      for (int n = 1; n <= 3; n++) {
        for (int l = 1; l <= 3; l++) {
          normals(1, l) += aurostd::eijk(l, m, n) * lattice(2, m) * lattice(3, n);
          normals(2, l) += aurostd::eijk(l, m, n) * lattice(3, m) * lattice(1, n);
          normals(3, l) += aurostd::eijk(l, m, n) * lattice(1, m) * lattice(2, n);
        }
      }
    }
    double length;
    for (int i = 1; i <= 3; i++) {
      length = aurostd::modulus(normals(i));
      for (int j = 1; j <= 3; j++) {
        normals(i, j) /= length;
      }
    }

    v_n_cpos[0] = -normals(1);
    v_n_cpos[1] = -normals(2);
    v_n_cpos[2] = -normals(3);
    v_n_cpos[3] = normals(1);
    v_n_cpos[4] = normals(2);
    v_n_cpos[5] = normals(3);
    //[CO20190520 - this does NOT work]i}

    if (LDEBUG) {
      for (uint i = 0; i < 6; i++) {
        cerr << __AFLOW_FUNC__ << " plane[i=" << i << "]: n=" << v_n_cpos[i] << ", p=" << (i < 3 ? p_origin : p_top) << endl;
      }
    }

    // find what atom fpos_starting corresponds to
    const aurostd::xvector<double> fpos_starting = BringInCell(c2f * cpos_starting);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cpos_starting=" << cpos_starting << endl;
      cerr << __AFLOW_FUNC__ << " fpos_starting=" << fpos_starting << endl;
    }
    uint starting_atom = AUROSTD_MAX_UINT;
    for (size_t i = 0; i < atoms.size(); i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " checking atoms[i=" << i << "].fpos=" << atoms[i].fpos << endl;
      }
      if (SYM::FPOSMatch(fpos_starting, atoms[i].fpos, lattice, f2c, skew, sym_eps)) {
        starting_atom = i;
        break;
      } // DX20190619 - lattice and f2c as input
    }
    if (starting_atom == AUROSTD_MAX_UINT) {
      if (!atoms2skip.empty()) {
        starting_atom = atoms2skip.back();
      } else {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot find starting atom", _INPUT_ERROR_);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " fpos_starting matches atom[" << starting_atom << "]" << endl;
    }
    if (starting_atom != AUROSTD_MAX_UINT && !aurostd::WithinList(atoms2skip, starting_atom)) {
      atoms2skip.push_back(starting_atom);
    } // no duplicates is better

    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> fpos_current=fpos_starting;
    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> cpos_current=f2c*fpos_current;
    aurostd::xvector<double> cpos_current;
    aurostd::xvector<double> cpos_final;
    aurostd::xvector<double> fpos_final;
    cpos_current = cpos_final = cpos_starting;
    aurostd::xvector<double> fpos_current = BringInCell(c2f * cpos_current);
    aurostd::xvector<double> fpos_current_prev;
    cpos_current = f2c * fpos_current; // bring current inside
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " fpos_current(start)=" << fpos_current << endl;
      cerr << __AFLOW_FUNC__ << " cpos_current(start)=" << cpos_current << endl;
    }

    aurostd::xvector<double> cdiff; // fdiff
    double dist_line = 0.0;
    double dist_rorigin = 0.0;
    double dist_rorigin_new = 0.0;
    aurostd::xvector<double> point_line_intersection_cpos;
    aurostd::xvector<double> line_plane_intersection_fpos;
    aurostd::xvector<double> line_plane_intersection_fpos_BIC;
    aurostd::xvector<double> line_plane_intersection_fpos_BIB; // point_line_intersection_fpos
    const double loop_shift = sym_eps / 2.0; // THIS IS CRITICAL //0.0; //introduces numerical inaccuracies, just get to unit cell boundary for BringInCell() //_ZERO_TOL_; //sym_eps //small bump into next loop
    double dist_rorigin_min = AUROSTD_MAX_DOUBLE;
    double dist_line_min = AUROSTD_MAX_DOUBLE; // must be positive
    uint ind_min = AUROSTD_MAX_UINT;
    //[CO20190520 - fixed more robustly with loop_shift]bool bring_2_opposite_boundary=false;
    while (loop_iteration < LOOP_ITERATION_MAX) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " loop_iteration=" << loop_iteration << endl;
        cerr << __AFLOW_FUNC__ << " fpos_current=" << fpos_current << endl;
        cerr << __AFLOW_FUNC__ << " cpos_current=" << cpos_current << endl;
        cerr << __AFLOW_FUNC__ << " atoms2skip=" << aurostd::joinWDelimiter(atoms2skip, ",") << endl;
      }
      // look for an atom in the path
      dist_rorigin_min = AUROSTD_MAX_DOUBLE;
      dist_line_min = AUROSTD_MAX_DOUBLE;
      ind_min = AUROSTD_MAX_UINT;
      if (!(loop_iteration == loop_iteration_starting && outside_current_cell)) { // go OUTSIDE current cell
        for (size_t i = 0; i < atoms.size(); i++) { // loop through all atoms, find nearest in line of sight
          if (loop_iteration == loop_iteration_starting && i == starting_atom) {
            continue;
          } // only for the first
          if (aurostd::WithinList(atoms2skip, (uint) i)) {
            continue;
          }
          //[CO20190520 - safe-guard for working in fractional space (skew)]point_line_intersection_fpos=aurostd::pointLineIntersection(fpos_current,l_fpos,atoms[i].fpos);
          //[CO20190520 - safe-guard for working in fractional space (skew)]point_line_intersection_cpos=f2c*point_line_intersection_fpos;
          point_line_intersection_cpos = aurostd::pointLineIntersection(cpos_current, l_cpos, atoms[i].cpos);
          cdiff = point_line_intersection_cpos - atoms[i].cpos; // distance from line (should be practically 0)
          dist_line = aurostd::modulus(cdiff); // distance from line (should be practical 0)
          if (dist_line < dist_line_min) {
            dist_line_min = dist_line;
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " atoms[i=" << i << "].cpos=" << atoms[i].cpos << ", atoms[i=" << i << "].fpos=" << atoms[i].fpos << endl;
            //[CO20190520 - safe-guard for working in fractional space (skew)]cerr << __AFLOW_FUNC__ << " point_line_intersection_fpos=" << point_line_intersection_fpos << endl;
            cerr << __AFLOW_FUNC__ << " point_line_intersection_cpos=" << point_line_intersection_cpos << endl;
            cerr << __AFLOW_FUNC__ << " dist_line=" << dist_line << endl;
          }
          cdiff = atoms[i].cpos - cpos_current; // distance from cpos_current (relative origin)
          dist_rorigin = aurostd::modulus(cdiff); // distance from cpos_current (relative origin)
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " dist_rorigin=" << dist_rorigin << endl;
          }
          if (dist_line < sym_eps && dist_rorigin < dist_rorigin_min) { // find point on line (dist_line) that is nearest to cpos_current/rorigin (smallest dist_rorigin)
            dist_rorigin_min = dist_rorigin;
            ind_min = i;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " dist_line_min=" << dist_line_min << endl;
      }
      if (ind_min != AUROSTD_MAX_UINT) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " cpos_final(pre)=" << cpos_final << endl;
        }
        //[CO20190520 - fractional space considerations]fdiff=fpos_current-atoms[ind_min].fpos; (within cell)
        //[CO20190520 - safe-guard for working in fractional space (skew)]fdiff=SYM::FPOSDistance(fpos_current,atoms[ind_min].fpos,lattice,c2f,f2c,skew); //CO20190520 - fractional space considerations (within
        // cell) [CO20190520 - safe-guard for working in fractional space (skew)]cdiff=f2c*fdiff;
        cdiff = atoms[ind_min].cpos - cpos_current; // distance from cpos_current (relative origin)
        dist_rorigin = aurostd::modulus(cdiff); // distance from cpos_current (relative origin)
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " dist_rorigin=" << dist_rorigin << endl;
        }
        cpos_final += cdiff;
        fpos_final = BringInCell(c2f * cpos_final);
        dist_rorigin_new = aurostd::modulus(f2c * SYM::FPOSDistFromFPOS(fpos_final, atoms[ind_min].fpos, lattice, c2f, f2c, skew)); // using fpos_final as a check //DX20190620 - changed function name
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " cpos_final(post)=" << cpos_final << endl;
          cerr << __AFLOW_FUNC__ << " fpos_final(post)=" << fpos_final << endl;
          cerr << __AFLOW_FUNC__ << " FOUND IMAGE! atom[i_atom=" << ind_min << ",fpos=" << atoms[ind_min].fpos << ",";
          cerr << "dist=" << dist_rorigin_new << "], ";
          cerr << "cpos_final=" << cpos_final << endl;
        }
        if (!aurostd::isequal(dist_rorigin_new, 0.0, sym_eps)) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Mismatch between atom.fpos and fpos_final", _RUNTIME_ERROR_);
        }
        atoms2skip.push_back(ind_min);
        //[CO20190520 - safe-guard for working in fractional space (skew)]return atoms[ind_min].fpos;
        return cpos_final;
      }
      // no atom found, so go to the next loop (ijk++)
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " no atoms found, loop_iteration++" << endl;
        cerr << __AFLOW_FUNC__ << " cpos_current=" << cpos_current << endl;
        cerr << __AFLOW_FUNC__ << " fpos_current=" << fpos_current << endl;
      }
      dist_rorigin_min = AUROSTD_MAX_DOUBLE;
      ind_min = AUROSTD_MAX_UINT;
      for (uint i = 0; i < 6; i++) {
        // this must be done in fractional, as
        //[CO20190520 - safe-guard for working in fractional space (skew)]v_line_plane_intersect[i]=aurostd::linePlaneIntersect( (i<3?p_origin:p_top),v_n_cpos[i],fpos_current,l_fpos,v_d[i],v_line_plane_intersection_fpos[i]);
        v_line_plane_intersect[i] = aurostd::linePlaneIntersect((i < 3 ? p_origin : p_top), v_n_cpos[i], cpos_current, l_cpos, v_d[i], v_line_plane_intersection_cpos[i]);
        line_plane_intersection_fpos = line_plane_intersection_fpos_BIB = c2f * v_line_plane_intersection_cpos[i];
        line_plane_intersection_fpos_BIC = BringInCell(line_plane_intersection_fpos);
        BringInBoundary(line_plane_intersection_fpos_BIB, _ZERO_TOL_); // small padding === no shift to opposite boundary
        if (!aurostd::isequal(line_plane_intersection_fpos, line_plane_intersection_fpos_BIB)) {
          v_line_plane_intersect[i] = false;
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " line-plane intersection is outside cell: line_plane_intersection_fpos=" << line_plane_intersection_fpos << endl;
          }
        }
        if (v_line_plane_intersect[i] &&
            !(aurostd::isequal(line_plane_intersection_fpos_BIC[1], 0.0) || aurostd::isequal(line_plane_intersection_fpos_BIC[2], 0.0) || aurostd::isequal(line_plane_intersection_fpos_BIC[3], 0.0))) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Line did not intersect plane", _RUNTIME_ERROR_);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " v_line_plane_intersect[i=" << i << "]=" << v_line_plane_intersect[i];
          cerr << ", v_d[i=" << i << "]=" << v_d[i] << endl;
          cerr << __AFLOW_FUNC__ << " line_plane_intersection[i=" << i << "].cpos=" << v_line_plane_intersection_cpos[i];
          cerr << ", fpos[i=" << i << "]=" << line_plane_intersection_fpos << endl;
        }
        if (v_line_plane_intersect[i] && v_d[i] > sym_eps) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " valid possible intersection at plane[i=" << i << "]" << endl;
          }
          if (v_d[i] < dist_rorigin_min) {
            dist_rorigin_min = v_d[i];
            ind_min = i;
          }
        }
      }

      if (ind_min == AUROSTD_MAX_UINT) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot find intersecting plane", _RUNTIME_ERROR_);
      }

      //[CO20190520 - do NOT take FPOSDistance(), we do not want to minimize vector]cpos_final += aurostd::modulus(f2c*SYM::FPOSDistance(fpos_current,v_line_plane_intersection_fpos[ind_min],lattice,c2f,f2c,skew)) + sym_eps;

      //[CO20190520 - safe-guard for working in fractional space (skew)]line_plane_intersection_cpos=f2c*v_line_plane_intersection_fpos[ind_min];
      //[CO20190520 - safe-guard for working in fractional space (skew)]line_plane_intersection_cpos=f2c*v_line_plane_intersection_fpos[ind_min];
      //[CO20190520 - safe-guard for working in fractional space (skew)]cdiff=cpos_current-line_plane_intersection_cpos;  //distance between cpos_current and unit cell boundary
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " found intersection plane = " << ind_min << endl;
        cerr << __AFLOW_FUNC__ << " cpos_current=" << cpos_current << endl;
        cerr << __AFLOW_FUNC__ << " fpos_current=" << fpos_current << endl;
      }
      cdiff = v_line_plane_intersection_cpos[ind_min] - cpos_current; // distance between cpos_current and unit cell boundary
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " cpos_final(pre)=" << cpos_final << endl;
      }
      cpos_final += cdiff + loop_shift * l_cpos; //(loop_shift * l_cpos / aurostd::modulus(l_cpos) ); //loop_shift is to ensure BringInCell() gets us to next ijk
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " cpos_final(post)=" << cpos_final << endl;
      }
      //[CO20190520 - safe-guard for working in fractional space (skew)]cpos_current=line_plane_intersection_cpos + (loop_shift * l_cpos / aurostd::modulus(l_cpos) );  //loop_shift is to ensure BringInCell() gets us to next ijk
      cpos_current = v_line_plane_intersection_cpos[ind_min] + loop_shift * l_cpos; //(loop_shift * l_cpos / aurostd::modulus(l_cpos) );  //loop_shift is to ensure BringInCell() gets us to next ijk
      fpos_current_prev = fpos_current;
      fpos_current = c2f * cpos_current;
      //[CO20190520 - need a smarter BringInCell() approach]fpos_current=BringInCell(c2f*cpos_current); //bring to boundary and BringInCell() to shift back to origin //BringInCell(): need to recalculate cpos (stay inside cell)
      //[NOT GOOD ENOUGH]fpos_current+=(-v_n_fpos[ind_min]); //a trick! add negative of boundary normal
      //[NOT GOOD ENOUGH]fpos_current=BringInCell(c2f*cpos_current); //we still need to BringInCell()
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " fpos_current(pre-boundary)=" << fpos_current << endl;
      }
      //[CO20190520 - loop shift prevents us from having to play with tols too much]bring_2_opposite_boundary=!(
      //[CO20190520 - loop shift prevents us from having to play with tols too much]    fpos_current_prev[1]<-_ZERO_TOL_ || fpos_current_prev[1]>1.0+_ZERO_TOL_ ||
      //[CO20190520 - loop shift prevents us from having to play with tols too much]    fpos_current_prev[2]<-_ZERO_TOL_ || fpos_current_prev[2]>1.0+_ZERO_TOL_ ||
      //[CO20190520 - loop shift prevents us from having to play with tols too much]    fpos_current_prev[3]<-_ZERO_TOL_ || fpos_current_prev[3]>1.0+_ZERO_TOL_ );
      //[CO20190520 - loop shift prevents us from having to play with tols too much]if(LDEBUG){cerr << __AFLOW_FUNC__ << " bring_2_opposite_boundary=" << bring_2_opposite_boundary << endl;}
      //[CO20190520 - fixed more robustly with loop_shift]if(bring_2_opposite_boundary){Bring2OppositeBoundary(fpos_current);}
      //[CO20190520 - plugged into BringInBoundary() with no padding]Bring2OppositeBoundary(fpos_current);
      BringInBoundary(fpos_current); // no padding === shift to opposite boundary
      cpos_current = f2c * fpos_current;

      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " v_line_plane_intersection_cpos[i=" << ind_min << "]=" << v_line_plane_intersection_cpos[ind_min] << endl;
        cerr << __AFLOW_FUNC__ << " v_line_plane_intersection_fpos[i=" << ind_min << "]=" << c2f * v_line_plane_intersection_cpos[ind_min] << endl;
        cerr << __AFLOW_FUNC__ << " fpos_current(post-boundary)=" << fpos_current << endl;
        cerr << __AFLOW_FUNC__ << " cpos_current(post-boundary)=" << cpos_current << endl;
      }

      atoms2skip.clear();
      loop_iteration++;
    }

    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot find next atom", _RUNTIME_ERROR_);
    //[CO20190520 - safe-guard for working in fractional space (skew)]return fpos_current;
    return cpos_final;
  }

  // returns back how many times you need to go in hkl direction before you return back to equivalent site
  // differs from getSpacingHKLPlane(), which considers ONLY lattice
  // getDistanceBetweenImages() considers both lattice + basis
  // outside_cell makes sure you loop outside the cell at least once
  double getDistanceBetweenImages(const xstructure& xstr_in, const aurostd::xvector<double>& n_cpos, bool outside_cell) { // CO20190320
    double dist = 0.0;
    if (distanceBetweenImages_HKL(xstr_in, n_cpos, dist, outside_cell)) {
      return dist;
    }
    if (distanceBetweenImages_Tracing(xstr_in, n_cpos, dist, outside_cell)) {
      return dist;
    }

    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot find distance between images", _INPUT_ERROR_);
    return dist;
  }
  bool distanceBetweenImages_HKL(const xstructure& xstr_in, const aurostd::xvector<double>& n_cpos, double& distance_between_images, bool outside_cell) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);
    distance_between_images = 0.0;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " starting" << endl;
      cerr << __AFLOW_FUNC__ << " xstr_in" << endl;
      cerr << xstr_in << endl;
      cerr << __AFLOW_FUNC__ << " n_cpos" << n_cpos << endl;
    }
    if (xstr_in.atoms.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No atoms found in xstructure", _INPUT_ERROR_);
    }

    // assume NO changes in xstr_in (too slow)
    const aurostd::xmatrix<double>& lattice = xstr_in.lattice;
    const aurostd::xmatrix<double>& f2c = xstr_in.f2c;
    const aurostd::xmatrix<double>& c2f = xstr_in.c2f;
    double min_dist = xstr_in.dist_nn_min;
    if (min_dist == AUROSTD_NAN) {
      min_dist = SYM::minimumDistance(xstr_in);
    }
    double sym_eps = xstr_in.sym_eps;
    if (sym_eps == AUROSTD_NAN) {
      sym_eps = SYM::defaultTolerance(xstr_in);
    }
    const bool skew = SYM::isLatticeSkewed(lattice, min_dist, sym_eps);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " min_dist=" << min_dist << endl;
      cerr << __AFLOW_FUNC__ << " sym_eps=" << sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }

    aurostd::xvector<int> hkl;
    if (!Normal2HKLPlane(lattice, n_cpos, hkl)) {
      return false;
    }
    const double d_spacing = getSpacingHKLPlane(lattice, hkl);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl=" << hkl << endl;
      cerr << __AFLOW_FUNC__ << " d_spacing=" << d_spacing << endl;
    }

    int count_d_spacings = 0;

    // use getFullSymBasis() to map atoms to same type
    _sym_op symop;
    symop.is_fgroup = true;
    symop.Uc = symop.Uf = aurostd::eye<double>(); // no rotation

    if (xstr_in.atoms.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No atoms found in xstructure", _INPUT_ERROR_);
    }
    vector<int> basis_atoms_map;
    vector<int> basis_types_map; // dummies
    uint loop_iteration = 0;
    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> fpos_prev;
    //[CO20190520 - safe-guard for working in fractional space (skew)]double dist_prev;
    bool found_map = false;
    aurostd::xvector<double> ftau_BIC;

    while (!found_map && loop_iteration < LOOP_ITERATION_MAX) {
      symop.ctau = ((double) ++count_d_spacings * d_spacing) * n_cpos;
      symop.ftau = c2f * symop.ctau;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " symop.ftau=" << symop.ftau << endl;
        cerr << __AFLOW_FUNC__ << " symop.ctau=" << symop.ctau << endl;
      }
      loop_iteration++;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " symop=" << endl;
        cerr << symop << endl;
      }
      if (outside_cell) {
        ftau_BIC = BringInCell(symop.ftau);
        if (aurostd::isequal(symop.ftau, ftau_BIC)) {
          continue;
        }
      }
      if (SYM::getFullSymBasis(xstr_in.atoms, lattice, c2f, f2c, symop, true, skew, sym_eps, basis_atoms_map, basis_types_map)) {
        found_map = true;
        break;
      }
    }

    if (loop_iteration == LOOP_ITERATION_MAX) {
      // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find distance between images",_RUNTIME_ERROR_);
      return false;
    }

    distance_between_images = (double) count_d_spacings++ * d_spacing;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " distance_between_images=" << distance_between_images << endl;
    }
    return distance_between_images;
  }

  bool distanceBetweenImages_Tracing(const xstructure& xstr_in, const aurostd::xvector<double>& n_cpos, double& distance_between_images, bool outside_cell) { // CO20190320
    const bool LDEBUG = (false || XHOST.DEBUG);
    distance_between_images = 0.0;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " starting" << endl;
      cerr << __AFLOW_FUNC__ << " xstr_in" << endl;
      cerr << xstr_in << endl;
      cerr << __AFLOW_FUNC__ << " n_cpos" << n_cpos << endl;
    }
    if (xstr_in.atoms.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No atoms found in xstructure", _INPUT_ERROR_);
    }

    // assume NO changes in xstr_in (too slow)
    const aurostd::xmatrix<double>& lattice = xstr_in.lattice;
    const aurostd::xmatrix<double>& f2c = xstr_in.f2c;
    const aurostd::xmatrix<double>& c2f = xstr_in.c2f;
    double min_dist = xstr_in.dist_nn_min;
    if (min_dist == AUROSTD_NAN) {
      min_dist = SYM::minimumDistance(xstr_in);
    }
    double sym_eps = xstr_in.sym_eps;
    if (sym_eps == AUROSTD_NAN) {
      sym_eps = SYM::defaultTolerance(xstr_in);
    }
    const bool skew = SYM::isLatticeSkewed(lattice, min_dist, sym_eps);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " min_dist=" << min_dist << endl;
      cerr << __AFLOW_FUNC__ << " sym_eps=" << sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }

    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> n_fpos=c2f*n_cpos;n_fpos/=aurostd::modulus(n_fpos);
    //[CO20190520 - safe-guard for working in fractional space (skew)]if(LDEBUG){
    //[CO20190520 - safe-guard for working in fractional space (skew)]  cerr << __AFLOW_FUNC__ << " n_cpos=" << n_cpos << endl;
    //[CO20190520 - safe-guard for working in fractional space (skew)]  cerr << __AFLOW_FUNC__ << " n_fpos=" << n_fpos << endl;
    //[CO20190520 - safe-guard for working in fractional space (skew)]}

    double cpos_diff = 0.0;

    // use getFullSymBasis() to map atoms to same type
    aurostd::xvector<double> cpos_prev;
    aurostd::xvector<double> cpos_new;
    aurostd::xvector<double> cpos_orig;
    aurostd::xvector<double> cpos_direct;
    _sym_op symop;
    symop.is_fgroup = true;
    symop.Uc = symop.Uf = aurostd::eye<double>(); // no rotation

    cpos_prev = xstr_in.atoms.front().cpos;
    vector<int> basis_atoms_map;
    vector<int> basis_types_map; // dummies
    uint loop_iteration = 0;
    uint loop_iteration_prev = 0;
    //[CO20190520 - safe-guard for working in fractional space (skew)]xvector<double> fpos_prev;
    //[CO20190520 - safe-guard for working in fractional space (skew)]double dist_prev;
    cpos_orig = cpos_prev = cpos_new = xstr_in.atoms.front().cpos; // initialize
    vector<uint> atoms2skip;
    atoms2skip.push_back(0); // first atom
    bool found_map = false;

    while (!found_map && loop_iteration < LOOP_ITERATION_MAX) {
      // there is some noise associated with moving atom to atom (of order of sym_eps)
      // it is better to increment by loop all at once
      /*
         symop.ftau=getNextAtomInPath(xstr_in,n_fpos,symop.ftau,atoms2skip,distance_between_images,loop_iteration,outside_cell);
         symop.ctau=f2c*symop.ftau;
         if(LDEBUG){
         cerr << __AFLOW_FUNC__ << " symop.ftau=" << symop.ftau << endl;
         cerr << __AFLOW_FUNC__ << " symop.ctau=" << symop.ctau << endl;
         }
         if(LDEBUG){cerr << __AFLOW_FUNC__ << " symop=" << endl;cerr << symop << endl;}
         if(SYM::getFullSymBasis(xstr_in.atoms,lattice,c2f,f2c,symop,true,skew,sym_eps,basis_atoms_map,basis_types_map)){break;}
         */
      cpos_prev = cpos_new;
      loop_iteration_prev = loop_iteration;
      while (!found_map && loop_iteration == loop_iteration_prev) {
        //[CO20190520 - safe-guard for working in fractional space (skew)]symop.ftau=getNextAtomInPath(xstr_in,n_fpos,fpos_prev,atoms2skip,distance_between_images,loop_iteration,outside_cell);
        //[CO20190520 - safe-guard for working in fractional space (skew)]symop.ctau=f2c*symop.ftau;
        cpos_new = getNextAtomInPath(xstr_in, n_cpos, cpos_prev, atoms2skip, loop_iteration, outside_cell);
        symop.ctau = cpos_new - cpos_orig;
        symop.ftau = c2f * symop.ctau;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " atoms2skip=" << aurostd::joinWDelimiter(atoms2skip, ",") << endl;
          cerr << __AFLOW_FUNC__ << " symop.ctau=" << symop.ctau << endl;
          cerr << __AFLOW_FUNC__ << " symop.ftau=" << symop.ftau << endl;
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " symop=" << endl;
          cerr << symop << endl;
        }
        if (SYM::getFullSymBasis(xstr_in.atoms, lattice, c2f, f2c, symop, true, skew, sym_eps, basis_atoms_map, basis_types_map)) {
          found_map = true;
          break;
        }
      }
      distance_between_images += aurostd::modulus(cpos_new - cpos_prev);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " cpos_prev=" << cpos_prev << endl;
        cerr << __AFLOW_FUNC__ << " cpos_new=" << cpos_new << endl;
        cerr << __AFLOW_FUNC__ << " distance_between_images=" << distance_between_images << endl;
        cerr << __AFLOW_FUNC__ << " cpos_new(pre)=" << cpos_new[1] << "," << cpos_new[2] << "," << cpos_new[3] << endl;
      }
      cpos_direct = cpos_orig + distance_between_images * n_cpos; // rectify for atoms close to line that did not yield getFullSymBasis()
      cpos_diff = aurostd::modulus(cpos_new - cpos_direct);
      cpos_new = cpos_direct;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " cpos_new(post)=" << cpos_new[1] << "," << cpos_new[2] << "," << cpos_new[3] << endl;
        cerr << __AFLOW_FUNC__ << " cpos_diff=" << cpos_diff << endl;
      }
      if (!aurostd::isequal(cpos_diff, 0.0, sym_eps)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "cpos_diff!=0 (cpos_diff=" + aurostd::utype2string(cpos_diff) + ")", _RUNTIME_ERROR_);
      }
    }

    if (loop_iteration == LOOP_ITERATION_MAX) {
      // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find distance between images",_RUNTIME_ERROR_);
      return false;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " distance_between_images=" << distance_between_images << endl;
    }

    return distance_between_images;
  }

  aurostd::xmatrix<double> getSlabLattice(istream& input, const aurostd::xvector<int>& hkl, aurostd::xmatrix<double>& lattice_slab_origbasis, double ang_dev, double vlen_max_strict) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return getSlabLattice(xstr_in, hkl, lattice_slab_origbasis, ang_dev, vlen_max_strict);
  }
  aurostd::xmatrix<double> getSlabLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, aurostd::xmatrix<double>& lattice_slab_origbasis, double ang_dev, double vlen_max_strict) {
    // ang_dev is acceptable angle deviation from v1Xv2 in degrees (5 degrees, W. Sun and G. Ceder, Surface Science 617 (2013) 53-59)
    // vlen_max_strict is the absolute limit for vlen (default -> infinity)
    // restricting vlen is NOT a priority (unless you want a particular v3), it is better to restrict angle via ang_dev
    // hence, the order of the default variables
    const bool LDEBUG = (false || XHOST.DEBUG);

    // assume NO changes in xstr_in (too slow)
    const aurostd::xmatrix<double>& lattice = xstr_in.lattice;
    const aurostd::xvector<double>& a1 = lattice(1);
    const aurostd::xvector<double>& a2 = lattice(2);
    const aurostd::xvector<double>& a3 = lattice(3);

    vector<aurostd::xvector<double>> intercepts; // plane intercepts
    intercepts = getHKLPlaneIntercepts(lattice, hkl);
    const aurostd::xvector<double> v1 = intercepts[1] - intercepts[0];
    const aurostd::xvector<double> v2 = intercepts[2] - intercepts[0];
    const aurostd::xvector<double> v1Xv2 = aurostd::vector_product(v1, v2); // pseudo v3, NOT a basis vector

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " v1=" << v1 << endl;
      cerr << __AFLOW_FUNC__ << " v2=" << v2 << endl;
      cerr << __AFLOW_FUNC__ << " v1Xv2=" << v1Xv2 << endl;
    }

    // need to search for v3
    int dim = 1; // initial search

    aurostd::xvector<double> v3_test;
    aurostd::xvector<double> v3;
    double adiff_v1;
    double adiff_v2;
    double adiff_v1Xv2;
    double adiff_max = AUROSTD_MAX_DOUBLE;
    double vlen = AUROSTD_MAX_DOUBLE;
    const double vlen_maxv1v2 = max(aurostd::modulus(v1), aurostd::modulus(v2));
    double vlen_max = min(vlen_maxv1v2, vlen_max_strict);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vlen_max=" << vlen_max << endl;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " searching for v3 with search constraint" << endl;
    } // tight search, this should be the best
    for (int i = -dim; i <= dim; i++) {
      for (int j = -dim; j <= dim; j++) {
        for (int k = -dim; k <= dim; k++) {
          if (!i && !j && !k) {
            continue;
          } // no vector
          v3_test = (double) i * a1 + (double) j * a2 + (double) k * a3;
          adiff_v1 = aurostd::angle(v1, v3_test);
          adiff_v2 = aurostd::angle(v2, v3_test);
          adiff_v1Xv2 = aurostd::angle(v1Xv2, v3_test);
          vlen = aurostd::modulus(v3_test);
          if (false && LDEBUG) {
            cerr << __AFLOW_FUNC__ << " i=" << i << ",j=" << j << ",k=" << k << endl;
            cerr << __AFLOW_FUNC__ << " v3_test=" << v3_test << endl;
            cerr << __AFLOW_FUNC__ << " adiff_v1=" << adiff_v1 << endl;
            cerr << __AFLOW_FUNC__ << " adiff_v2=" << adiff_v2 << endl;
            cerr << __AFLOW_FUNC__ << " adiff_v1Xv2=" << adiff_v1Xv2 << endl;
            cerr << __AFLOW_FUNC__ << " vlen=" << vlen << " ?<= " << vlen_max << " == " << bool(vlen <= vlen_max) << endl;
          }
          if ((adiff_v1 >= 0.0 && adiff_v1 <= PI / 2.0) && // angle with v1 must be within acceptable range
              (adiff_v2 >= 0.0 && adiff_v2 <= PI / 2.0) && // angle with v2 must be within acceptable range
              (adiff_v1Xv2 >= 0.0 && adiff_v1Xv2 <= PI / 2.0) && // angle with v1Xv2 must be within acceptable range
              (adiff_v1Xv2 < adiff_max) && // angle with v1Xv2 is minimized
              (vlen <= vlen_max) // look for vectors same length or smaller than max(||v1||,||v2||) (initial constraint to keep cell size from EXPLODING)
          ) {
            adiff_max = adiff_v1Xv2;
            v3 = v3_test;
            vlen_max = vlen;
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " adiff_max[i=" << i << ",j=" << j << ",k=" << k << "](degrees)=" << rad2deg * adiff_max << endl;
              cerr << __AFLOW_FUNC__ << " v3[i=" << i << ",j=" << j << ",k=" << k << "]=" << v3 << endl;
            }
          }
        }
      }
    }

    if (adiff_max > deg2rad * ang_dev) { // if not found or it's greater than acceptable angle deviation, relax length constraint
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " searching for v3 WITHOUT search constraint" << endl;
      }
      const double radius = RadiusSphereLattice(lattice);
      const aurostd::xvector<int> dims = LatticeDimensionSphere(lattice, radius);
      dim = max(dims); //+1  //do not go too far out
      const int dim_found = dim;
      vlen_max = vlen_max_strict;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " dim=" << dim << endl;
      }
      for (int i = -dim; i <= dim && i <= dim_found; i++) {
        for (int j = -dim; j <= dim && j <= dim_found; j++) {
          for (int k = -dim; k <= dim && k <= dim_found; k++) {
            if (!i && !j && !k) {
              continue;
            } // no vector
            v3_test = (double) i * a1 + (double) j * a2 + (double) k * a3;
            adiff_v1 = aurostd::angle(v1, v3_test);
            adiff_v2 = aurostd::angle(v2, v3_test);
            adiff_v1Xv2 = aurostd::angle(v1Xv2, v3_test);
            vlen = aurostd::modulus(v3_test);
            if (false && LDEBUG) {
              cerr << __AFLOW_FUNC__ << " i=" << i << ",j=" << j << ",k=" << k << endl;
              cerr << __AFLOW_FUNC__ << " v3_test=" << v3_test << endl;
              cerr << __AFLOW_FUNC__ << " adiff_v1=" << adiff_v1 << endl;
              cerr << __AFLOW_FUNC__ << " adiff_v2=" << adiff_v2 << endl;
              cerr << __AFLOW_FUNC__ << " adiff_v1Xv2=" << adiff_v1Xv2 << endl;
              cerr << __AFLOW_FUNC__ << " vlen=" << vlen << " ?<= " << vlen_max << " == " << bool(vlen <= vlen_max) << endl;
            }
            if ((adiff_v1 >= 0.0 && adiff_v1 <= PI / 2.0) && // angle with v1 must be within acceptable range
                (adiff_v2 >= 0.0 && adiff_v2 <= PI / 2.0) && // angle with v2 must be within acceptable range
                (adiff_v1Xv2 >= 0.0 && adiff_v1Xv2 <= PI / 2.0) && // angle with v1Xv2 must be within acceptable range
                (adiff_v1Xv2 < adiff_max) && // angle with v1Xv2 is minimized
                (vlen <= vlen_max) // look for vectors smaller than current vlen_max
            ) {
              adiff_max = adiff_v1Xv2;
              v3 = v3_test;
              if (adiff_max <= deg2rad * ang_dev) { // constrain search here, we found a good one
                vlen_max = vlen;
                //[not needed, add if you want]if(dim_found==dim){dim_found=max(abs(i),max(abs(j),abs(k)));}  //max of abs(i),abs(j),abs(k)
              }
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " adiff_max[i=" << i << ",j=" << j << ",k=" << k << "](degrees)=" << rad2deg * adiff_max << endl;
                cerr << __AFLOW_FUNC__ << " v3[i=" << i << ",j=" << j << ",k=" << k << "]=" << v3 << endl;
              }
            }
          }
        }
      }
    }

    if (adiff_max == AUROSTD_MAX_DOUBLE) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot find acceptable v3", _INPUT_ERROR_);
    }

    // try linear combinations with v1, v2
    bool try_once = true;
    int multiplier = 1;
    vlen_max = min(vlen_maxv1v2, vlen_max_strict);
    while (try_once || (double) multiplier * aurostd::modulus(v3) <= vlen_max) {
      v3 *= (double) multiplier++;
      adiff_v1Xv2 = aurostd::angle(v1Xv2, v3); // shouldn't change
      while (true) {
        if (rad2deg * aurostd::angle(v1Xv2, v3) < 5.0) {
          break;
        } else if (aurostd::modulus(v3 - v1) <= vlen_max && rad2deg * aurostd::angle(v1Xv2, v3 - v1) < adiff_v1Xv2) {
          v3 -= v1;
          adiff_v1Xv2 = aurostd::angle(v1Xv2, v3);
        } else if (aurostd::modulus(v3 + v1) <= vlen_max && rad2deg * aurostd::angle(v1Xv2, v3 + v1) < adiff_v1Xv2) {
          v3 += v1;
          adiff_v1Xv2 = aurostd::angle(v1Xv2, v3);
        } else if (aurostd::modulus(v3 - v2) <= vlen_max && rad2deg * aurostd::angle(v1Xv2, v3 - v2) < adiff_v1Xv2) {
          v3 -= v2;
          adiff_v1Xv2 = aurostd::angle(v1Xv2, v3);
        } else if (aurostd::modulus(v3 + v2) <= vlen_max && rad2deg * aurostd::angle(v1Xv2, v3 + v2) < adiff_v1Xv2) {
          v3 += v2;
          adiff_v1Xv2 = aurostd::angle(v1Xv2, v3);
        } else {
          break;
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " found viable v3=" << v3;
        cerr << ", len(v3)=" << aurostd::modulus(v3);
        cerr << ", adiff_v1Xv2=" << aurostd::angle(v1Xv2, v3) << endl;
      }
      try_once = false;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " v1=" << v1 << endl;
      cerr << __AFLOW_FUNC__ << " v2=" << v2 << endl;
      cerr << __AFLOW_FUNC__ << " v3=" << v3 << endl;
      cerr << __AFLOW_FUNC__ << " ||v3||=" << aurostd::modulus(v3) << endl;
      cerr << __AFLOW_FUNC__ << " adiff_v3=" << adiff_max << endl;
    }

    // load in vectors as columns of matrix (as opposed to usual rows) for householder, we will transpose later
    lattice_slab_origbasis(1, 1) = v1[1];
    lattice_slab_origbasis(1, 2) = v2[1];
    lattice_slab_origbasis(1, 3) = v3[1];
    lattice_slab_origbasis(2, 1) = v1[2];
    lattice_slab_origbasis(2, 2) = v2[2];
    lattice_slab_origbasis(2, 3) = v3[2];
    lattice_slab_origbasis(3, 1) = v1[3];
    lattice_slab_origbasis(3, 2) = v2[3];
    lattice_slab_origbasis(3, 3) = v3[3];

    // orthogonalize as much as possible (rotate)
    aurostd::xmatrix<double> Q;
    aurostd::xmatrix<double> lattice_slab_newbasis;
    QRDecomposition_HouseHolder(lattice_slab_origbasis, Q, lattice_slab_newbasis);

    // immediate check for NEGATIVE determinant, VASP has issues here (see aflow_ivasp.cpp for negative triple product)
    // fix with S diagonal matrix with +-1: http://www.math.purdue.edu/~kkloste/cs515fa14/qr-uniqueness.pdf
    double trip_prod = det(lattice_slab_newbasis);
    if (trip_prod < 0.0) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " applying triple product correction" << endl;
      }
      int n_neg = 1;
      aurostd::xcombos xc;
      vector<int> v_combo;
      aurostd::xmatrix<double> S = aurostd::eye<double>();
      aurostd::xmatrix<double> Q_tmp;
      aurostd::xmatrix<double> lat_tmp;
      int index;
      int row;
      int col;
      while (trip_prod < 0.0) {
        if (n_neg > S.rows) {
          break;
        }
        xc.reset(lattice_slab_newbasis.cols, n_neg++, 'C');
        while (xc.increment()) {
          v_combo = xc.getCombo();
          S = aurostd::eye<double>();
          for (size_t i = 0; i < v_combo.size(); i++) {
            if (v_combo[i] == 1) {
              index = i + (S.urows - 1); // shift index to start by negating z first (keep x-y as calculated)
              row = aurostd::boundary_conditions_periodic(S.lrows, S.urows, index + S.lrows);
              col = aurostd::boundary_conditions_periodic(S.lcols, S.ucols, index + S.lcols);
              S(row, col) = -1;
            }
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " S=" << endl;
            cerr << S << endl;
          }
          Q_tmp = Q * S;
          lat_tmp = S * lattice_slab_newbasis;
          if (!aurostd::isequal(lattice_slab_origbasis, Q_tmp * lat_tmp)) {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "S matrix is not a viable triple-product correction", _RUNTIME_ERROR_);
          }
          trip_prod = det(lat_tmp);
          if (trip_prod >= 0.0) {
            Q = Q_tmp;
            lattice_slab_newbasis = lat_tmp;
            break;
          }
        }
      }
      if (det(lattice_slab_newbasis) < 0.0) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Triple product of lattice remains negative despite attempts to rectify", _INPUT_ERROR_);
      }
    }

    lattice_slab_newbasis = aurostd::roundoff(lattice_slab_newbasis, 1e-12); // even smaller than _ZERO_TOL_  //do round-off - here it is important because householder introduces some noise (order of machine epsilon)
    lattice_slab_newbasis = trasp(lattice_slab_newbasis);
    lattice_slab_origbasis = trasp(lattice_slab_origbasis); // do AFTER triple product correction

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Q=" << endl;
      cerr << Q << endl;
      cerr << __AFLOW_FUNC__ << " lattice_slab_origbasis=" << endl;
      cerr << lattice_slab_origbasis << endl;
      cerr << __AFLOW_FUNC__ << " lattice_slab_newbasis=" << endl;
      cerr << lattice_slab_newbasis << endl;
    }

    return lattice_slab_newbasis;
  }

  // follows procedure outlined in: W. Sun and G. Ceder, Surface Science 617 (2013) 53-59
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, double v3len_max_strict, ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, aflags, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, const _aflags& aflags, double v3len_max_strict, ostream& oss) {
    ofstream FileMESSAGE;
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, aflags, FileMESSAGE, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, ofstream& FileMESSAGE, double v3len_max_strict, ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, aflags, FileMESSAGE, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, const _aflags& aflags, ofstream& FileMESSAGE, double v3len_max_strict, ostream& oss) {
    aurostd::xmatrix<double> rotation;
    xstructure xstr_slab_newbasis; // xstr_rotated
    vector<int> sc2pcMap_slab;
    vector<int> pc2scMap_slab;
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, vector<int>& sc2pcMap_slab, vector<int>& pc2scMap_slab, double v3len_max_strict, ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, sc2pcMap_slab, pc2scMap_slab, aflags, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(
      const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, vector<int>& sc2pcMap_slab, vector<int>& pc2scMap_slab, const _aflags& aflags, double v3len_max_strict, ostream& oss) {
    ofstream FileMESSAGE;
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(
      const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, vector<int>& sc2pcMap_slab, vector<int>& pc2scMap_slab, ofstream& FileMESSAGE, double v3len_max_strict, ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  } // CO20190321
  xstructure CreateSlab_SurfaceLattice(
      const xstructure& xstr_in, const aurostd::xvector<int>& hkl, int total_layers, double vacuum, vector<int>& sc2pcMap_slab, vector<int>& pc2scMap_slab, const _aflags& aflags, ofstream& FileMESSAGE, double v3len_max_strict, ostream& oss) {
    aurostd::xmatrix<double> rotation;
    xstructure xstr_slab_newbasis;
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  } // CO20190321
  // STOP - EASY INPUTS
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       double v3len_max_strict,
                                       ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict,
                                       ostream& oss) {
    ofstream FileMESSAGE;
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       istream& input,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(vpflow, xstr_in, hkl, total_layers, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,
                                       const xstructure& xstr_in,
                                       aurostd::xvector<int>& hkl,
                                       int& total_layers,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    const std::streamsize prec = 8;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - read flags
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " reading flags" << endl;
    }

    int h_i = 1;
    int k_i = 1;
    int l_i = 1; // hkl of interest
    total_layers = DEFAULT_TOTAL_LAYERS; // size of supercell (~ 2x layers)
    double vacuum = 15; // vacuum in Angstroms

    vector<string> tokens;
    aurostd::string2tokens(vpflow.getattachedscheme("CREATE_SLAB::PLANE_INTEREST"), tokens, ",");
    if (tokens.size() == 3) {
      h_i = aurostd::string2utype<int>(tokens[0]);
      k_i = aurostd::string2utype<int>(tokens[1]);
      l_i = aurostd::string2utype<int>(tokens[2]);
    }
    hkl[1] = h_i;
    hkl[2] = k_i;
    hkl[3] = l_i;
    if (hkl[1] == 0 && hkl[2] == 0 && hkl[3] == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "hkl=(0,0,0)", _INPUT_ERROR_);
    }
    const string total_layers_string = vpflow.getattachedscheme("CREATE_SLAB::TOTAL_LAYERS");
    if (aurostd::isfloat(total_layers_string)) {
      const int _total_layers = aurostd::string2utype<int>(total_layers_string);
      if (_total_layers > 0) {
        total_layers = _total_layers;
      }
    }
    const string vacuum_string = vpflow.getattachedscheme("CREATE_SLAB::VACUUM");
    if (aurostd::isfloat(vacuum_string)) {
      const double _vacuum = aurostd::string2utype<double>(vacuum_string);
      if (_vacuum > 0) {
        vacuum = _vacuum;
      }
    }

    const std::streamsize prec_original = message.precision(); // original
    const std::ios_base::fmtflags ff_original = message.flags(); // original
    message.precision(prec);
    message.unsetf(std::ios_base::floatfield);

    message << "plane_interest=" << hkl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    message << "total_layers=" << total_layers;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    message << "vacuum=" << vacuum;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    message.precision(prec_original); // set back
    message.flags(ff_original); // set back

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - read flags
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }

  xstructure CreateSlab_SurfaceLattice(istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       double v3len_max_strict,
                                       ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       double v3len_max_strict,
                                       ostream& oss) {
    ofstream FileMESSAGE;
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    _aflags aflags;
    aflags.Directory = ".";
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(istream& input,
                                       const aurostd::xvector<int>& hkl,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const xstructure xstr_in(input, IOAFLOW_AUTO);
    return CreateSlab_SurfaceLattice(xstr_in, hkl, total_layers, vacuum, rotation, xstr_slab_newbasis, sc2pcMap_slab, pc2scMap_slab, aflags, FileMESSAGE, v3len_max_strict, oss);
  }
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,
                                       const aurostd::xvector<int>& hkl_i,
                                       int total_layers,
                                       double vacuum,
                                       aurostd::xmatrix<double>& rotation,
                                       xstructure& xstr_slab_newbasis,
                                       vector<int>& sc2pcMap_slab,
                                       vector<int>& pc2scMap_slab,
                                       const _aflags& aflags,
                                       ofstream& FileMESSAGE,
                                       double v3len_max_strict,
                                       ostream& oss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    const bool check_min_dist = true; // turn off if it gets too slow
    int count_check_min_dist = 0;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " starting" << endl;
    }

    const int xy_dims = 1; // dimensions of supercell in x-y dimensions
    const aurostd::xvector<double> zero_xvector; // zero aurostd::xvector //DX20201124

    xstructure xstr_bulk(xstr_in);
    xstr_bulk.ReScale(1.0); // do NOT modify further
    double min_dist = xstr_bulk.dist_nn_min;
    if (min_dist == AUROSTD_NAN) {
      min_dist = xstr_bulk.dist_nn_min = SYM::minimumDistance(xstr_bulk);
    }
    const double min_dist_orig = min_dist;
    double sym_eps = xstr_bulk.sym_eps;
    if (sym_eps == AUROSTD_NAN) {
      sym_eps = xstr_bulk.sym_eps = SYM::defaultTolerance(xstr_bulk);
    }
    const bool skew = SYM::isLatticeSkewed(xstr_bulk.lattice, min_dist, sym_eps);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_in=" << endl;
      cerr << xstr_in << endl;
      cerr << __AFLOW_FUNC__ << " min_dist=" << min_dist << endl;
      cerr << __AFLOW_FUNC__ << " sym_eps=" << sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }

    message << "Constructing slab (surface lattice) along (" << aurostd::joinWDelimiter(hkl_i, ",") << ")";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    if (check_min_dist) { // sanity check as we rotate structure/atoms
      min_dist = xstr_bulk.MinDist();
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;
      }
      if (!aurostd::isequal(min_dist_orig, min_dist)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Minimum distance changed", _VALUE_ERROR_);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - defining hkl normals
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " defining HKL normals" << endl;
    }

    const aurostd::xvector<double> n_i = HKLPlane2Normal(xstr_bulk.lattice, hkl_i);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n_i[h=" << hkl_i << "]=" << n_i << endl;
    }

    // quick test to make sure everything works
    aurostd::xvector<int> hkl;
    vector<aurostd::xvector<double>> intercepts; // plane intercepts
    const aurostd::xvector<double> v1;
    const aurostd::xvector<double> v2;
    const aurostd::xvector<double> v3; // plane-defining vectors (need only two)
    // test that we can go back and forth between n and hkl
    if (!Normal2HKLPlane(xstr_bulk.lattice, n_i, hkl)) {
      message << "Cannot convert normal -> (hkl): normal=" << n_i;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl_i=" << hkl_i << endl;
      cerr << __AFLOW_FUNC__ << " n(hkl_i)=" << n_i << endl;
      cerr << __AFLOW_FUNC__ << " hkl_i(test)=" << hkl << endl;
    }
    if (!aurostd::isequal(hkl_i, hkl)) {
      message << "Normal2HKLPlane() function failed on hkl_i=" << hkl_i << " (Normal2HKLPlane(n_i=" << n_i << ")=" << hkl << ")" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // test that hkl plane is orthogonal to n
    intercepts = getHKLPlaneIntercepts(xstr_bulk.lattice, hkl_i);
    if (LDEBUG) {
      for (size_t i = 0; i < intercepts.size(); i++) {
        cerr << __AFLOW_FUNC__ << " intercepts[" << i << "]=" << intercepts[i] << endl;
      }
    }

    const aurostd::xvector<int> hkl_test;
    aurostd::xvector<double> n_test;
    // test that we can go back and forth between n and hkl
    hkl_test[1] = 7;
    hkl_test[2] = 3;
    hkl_test[3] = 2;
    n_test = HKLPlane2Normal(xstr_bulk.lattice, hkl_test);
    if (!Normal2HKLPlane(xstr_bulk.lattice, n_test, hkl)) {
      message << "Cannot convert normal -> (hkl): normal=" << n_test;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hkl_test=" << hkl_test << endl;
      cerr << __AFLOW_FUNC__ << " n(hkl_test)=" << n_test << endl;
      cerr << __AFLOW_FUNC__ << " hkl_test(test)=" << hkl << endl;
    }
    if (!aurostd::isequal(hkl_test, hkl)) {
      message << "Normal2HKLPlane() function failed on hkl=" << hkl_test << " (Normal2HKLPlane(n_i=" << n_test << ")=" << hkl << ")" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // test that hkl plane is orthogonal to n
    intercepts = getHKLPlaneIntercepts(xstr_bulk.lattice, hkl_test);
    if (LDEBUG) {
      for (size_t i = 0; i < intercepts.size(); i++) {
        cerr << __AFLOW_FUNC__ << " intercepts[" << i << "]=" << intercepts[i] << endl;
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - defining hkl normals
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - create slab structure
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (check_min_dist) { // sanity check as we rotate structure/atoms
      min_dist = xstr_bulk.MinDist();
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;
      }
      if (!aurostd::isequal(min_dist_orig, min_dist)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Minimum distance changed", _VALUE_ERROR_);
      }
    }

    // test whether fpos/cpos work
    if (LDEBUG) {
      cerr << xstr_bulk << endl;
    }
    aurostd::xvector<double> fpos;
    aurostd::xvector<double> cpos;
    for (size_t i = 0; i < xstr_bulk.atoms.size(); i++) {
      fpos = xstr_bulk.c2f * xstr_bulk.atoms[i].cpos;
      cpos = xstr_bulk.f2c * xstr_bulk.atoms[i].fpos;
      if (!aurostd::isequal(xstr_bulk.atoms[i].fpos, fpos)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "atoms[i=" + aurostd::utype2string(i) + "].fpos mismatch", _INPUT_ERROR_);
      }
      if (!aurostd::isequal(xstr_bulk.atoms[i].cpos, cpos)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "atoms[i=" + aurostd::utype2string(i) + "].cpos mismatch", _INPUT_ERROR_);
      }
    }

    stringstream title;
    xstr_slab_newbasis.clear(); // DX20191220 - uppercase to lowercase clear
    xstructure xstr_slab_origbasis;
    xstr_slab_newbasis.lattice = getSlabLattice(xstr_bulk, hkl_i, xstr_slab_origbasis.lattice, DEFAULT_V3_ANGLE_DEVIATION, v3len_max_strict);
    xstr_slab_newbasis.FixLattices(); // ang_dev==5.0 is standard (DEFAULT_V3_ANGLE_DEVIATION), the vlen_max_strict is very important here, as the test from Sun et al. takes a shortcut here
    rotation = trasp(xstr_slab_newbasis.lattice) * inverse(trasp(xstr_slab_origbasis.lattice));
    //  inverse(xstr_slab_origbasis.lattice)*xstr_slab_newbasis.lattice; //equivalent to trasp( trasp(new) * inv(trasp(orig)) )  //new = rotation * orig

    // quick check
    // xstr_slab_newbasis.lattice[1][1]=4.736233;xstr_slab_newbasis.lattice[1][2]=0.0;xstr_slab_newbasis.lattice[1][3]=0.0;
    // xstr_slab_newbasis.lattice[2][1]=9.472473;xstr_slab_newbasis.lattice[2][2]=21.242108;xstr_slab_newbasis.lattice[2][3]=0.0;
    // xstr_slab_newbasis.lattice[3][1]=2.368118;xstr_slab_newbasis.lattice[3][2]=3.168035;xstr_slab_newbasis.lattice[3][3]=2.605281;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_slab_origbasis.lattice=" << endl;
      cerr << xstr_slab_origbasis.lattice << endl;
      cerr << __AFLOW_FUNC__ << " xstr_slab_newbasis.lattice=" << endl;
      cerr << xstr_slab_newbasis.lattice << endl;
      cerr << __AFLOW_FUNC__ << " rotation=" << endl;
      cerr << rotation << endl;
      cerr << __AFLOW_FUNC__ << " xstr_slab_newbasis.c2f=" << endl;
      cerr << xstr_slab_newbasis.c2f << endl;
      cerr << __AFLOW_FUNC__ << " xstr_slab_newbasis.f2c=" << endl;
      cerr << xstr_slab_newbasis.f2c << endl;
    }

    const double volume_original = std::abs(aurostd::det(xstr_bulk.lattice));
    const double volume_new = std::abs(aurostd::det(xstr_slab_origbasis.lattice));
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " volume_original=" << volume_original << endl;
      cerr << __AFLOW_FUNC__ << " volume_new=" << volume_new << endl;
      cerr << __AFLOW_FUNC__ << " sym_eps=" << sym_eps << endl;
      cerr << __AFLOW_FUNC__ << " skew=" << skew << endl;
    }
    // fold_in_only=false;
    // VERY IMPORTANT, use xstr_slab_origbasis.lattice and not xstr_slab_newbasis.lattice
    // xstr_slab_newbasis.lattice is rotated relative to original lattice
    deque<_atom> atoms = foldAtomsInCell(xstr_bulk, xstr_slab_origbasis.lattice, skew, sym_eps); // do NOT fold_in_only
    xstr_slab_origbasis.ReplaceAtoms(atoms);

    // clean up structure
    xstr_slab_origbasis.ReScale(1.0);
    xstr_slab_origbasis.ShiftOriginToAtom(0);
    xstr_slab_origbasis.origin = zero_xvector; // reset origin //DX20201124
    xstr_slab_origbasis.BringInCell();
    xstr_slab_origbasis.clean(); // DX20191220 - uppercase to lowercase clean

    // set title
    title.str("");
    title << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr_bulk.title) << " (SLAB surface lattice original basis: ";
    title << "hkl=(" << aurostd::joinWDelimiter(hkl_i, ",") << "), ";
    title << "total_layers=" << total_layers << ", ";
    title << "vacuum=" << vacuum << ")";
    xstr_slab_newbasis.title = title.str();

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_slab_origbasis=" << endl;
      cerr << xstr_slab_origbasis << endl;
    }

    // now we can convert cpos for xstr_slab_newbasis.lattice
    for (size_t i = 0; i < atoms.size(); i++) {
      atoms[i].cpos = xstr_slab_newbasis.f2c * atoms[i].fpos;
    } // rotate to new lattice
    xstr_slab_newbasis.ReplaceAtoms(atoms);

    // clean up structure
    xstr_slab_newbasis.ReScale(1.0);
    xstr_slab_newbasis.ShiftOriginToAtom(0);
    xstr_slab_newbasis.origin = zero_xvector; // reset origin
    xstr_slab_newbasis.BringInCell();
    xstr_slab_newbasis.clean(); // DX20191220 - uppercase to lowercase clean

    // set title
    title.str("");
    title << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr_bulk.title) << " (SLAB surface lattice new basis: ";
    title << "hkl=(" << aurostd::joinWDelimiter(hkl_i, ",") << "), ";
    title << "total_layers=" << total_layers << ", ";
    title << "vacuum=" << vacuum << ")";
    xstr_slab_newbasis.title = title.str();

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_slab_newbasis=" << endl;
      cerr << xstr_slab_newbasis << endl;
    }

    if (check_min_dist) { // sanity check as we rotate structure/atoms
      min_dist = xstr_slab_newbasis.MinDist();
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;
      }
      if (!aurostd::isequal(min_dist_orig, min_dist)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Minimum distance changed", _VALUE_ERROR_);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - create slab structure
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - resolve layers count
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " resolving layers count" << endl;
    }

    const double d_spacing = slab::getSpacingHKLPlane(xstr_bulk, hkl_i); // aurostd::modulus(xstr_slab.lattice(1))/sqrt(h_s*h_s+k_s*k_s+l_s*l_s);
    const double d_layers = slab::getDistanceBetweenImages(xstr_bulk, n_i, false); // this depends on UN-ROTATED lattice
    const double d_cells = slab::getDistanceBetweenImages(xstr_bulk, n_i, true); // go outside cell
    const int layers_per_cell = (int) (d_cells / d_layers); // floor
    const int supercell_layers = (total_layers + layers_per_cell - 1) / layers_per_cell; // ceil //(double)total_layers;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " n_i[h=" << hkl_i << "]=" << n_i << endl;
      cerr << __AFLOW_FUNC__ << " d_spacing=" << d_spacing << endl;
      cerr << __AFLOW_FUNC__ << " d_layers=" << d_layers << endl;
      cerr << __AFLOW_FUNC__ << " d_cells=" << d_cells << endl;
      cerr << __AFLOW_FUNC__ << " abs(d_layers-d_cells)=" << std::abs(d_layers - d_cells) << endl;
      cerr << __AFLOW_FUNC__ << " layers_per_cell=" << layers_per_cell << endl;
      cerr << __AFLOW_FUNC__ << " supercell_layers=" << supercell_layers << endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - resolve layers count
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - create supercell
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " creating supercell" << endl;
    }

    // now create a supercell
    const aurostd::xmatrix<double> supercell_mat;
    supercell_mat(1, 1) = (double) xy_dims;
    supercell_mat(2, 2) = (double) xy_dims;
    supercell_mat(3, 3) = supercell_layers;
    xstructure xstr_slab = GetSuperCell(xstr_slab_newbasis, supercell_mat, sc2pcMap_slab, pc2scMap_slab, false, false, false);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_slab=" << endl;
      cerr << xstr_slab << endl;
    }
    if (check_min_dist) { // sanity check as we rotate structure/atoms
      min_dist = xstr_slab.MinDist();
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;
      }
      if (!aurostd::isequal(min_dist_orig, min_dist)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Minimum distance changed", _VALUE_ERROR_);
      }
    }

    // clean up structure
    xstr_slab.ReScale(1.0);
    xstr_slab.ShiftOriginToAtom(0);
    xstr_slab.origin = zero_xvector; // reset origin //DX20201124
    xstr_slab.BringInCell();
    // xstr_slab.clean();  //clear origin! //do not clear ijk! origin is okay here, only a problem for Rotate() //DX20191220 - uppercase to lowercase clean

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - create supercell
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START - add vacuum
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " adding vacuum" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " old_c_lattice=" << endl;
      cerr << xstr_slab.lattice << endl;
    }
    aurostd::xvector<double> new_c_lattice = xstr_slab.lattice(3);
    new_c_lattice += vacuum * new_c_lattice / aurostd::modulus(new_c_lattice);
    xstr_slab.lattice[3][1] = new_c_lattice(1);
    xstr_slab.lattice[3][2] = new_c_lattice(2);
    xstr_slab.lattice[3][3] = new_c_lattice(3);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " new_c_lattice=" << endl;
      cerr << xstr_slab.lattice << endl;
    }

    // fix fpos
    xstr_slab.FixLattices();
    const aurostd::xmatrix<double>& c2f = xstr_slab.c2f;
    for (size_t i = 0; i < xstr_slab.atoms.size(); i++) {
      xstr_slab.atoms[i].fpos = c2f * xstr_slab.atoms[i].cpos;
    }

    if (check_min_dist) { // sanity check as we rotate structure/atoms
      min_dist = xstr_slab.MinDist();
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " mindist[" << count_check_min_dist++ << "]=" << min_dist << endl;
      }
      if (!aurostd::isequal(min_dist_orig, min_dist)) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Minimum distance changed", _VALUE_ERROR_);
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // STOP - add vacuum
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    message << "Slab (surface lattice) along (" << aurostd::joinWDelimiter(hkl_i, ",") << ") constructed";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

    return xstr_slab;
  }

} // namespace slab
// CO20190601 STOP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
