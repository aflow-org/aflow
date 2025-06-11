
#ifndef AFLOW_DEFECTS_H
#define AFLOW_DEFECTS_H

#include <deque>
#include <ostream>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xvector.h"

#include "flow/aflow_xclasses.h"
#include "structure/aflow_xatom.h"

class acage {
public:
    // constructors/destructors                                   // --------------------------------------
  acage();                                                      // constructor default
  acage(const acage& b);                                        // constructor copy
  ~acage();                                                     // destructor
    // CONTENT                                                    // --------------------------------------
  void clear();                                             // clear everything
  aurostd::xvector<double> origin_fpos;// default 3D
  aurostd::xvector<double> origin_cpos;// default 3D
  double radius;
  uint coordination_position;
  int cages_position;
  int cages_irrtype;
  std::deque<_atom> atoms;

private:                                                        // ---------------------------------------
  void free();                                                  // to free everything
  void copy(const acage& b);                                    // the flag is necessary because sometimes you need to allocate the space.
};
class _isort_acage_radius {                   // sorting through reference
public:
  bool operator()(const acage& a1, const acage& a2) const { return (bool) (a1.radius > a2.radius); }
};

bool GetSphereFromFourPoints(aurostd::xvector<double>& orig, double& radius, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3, const aurostd::xvector<double>& v4);
bool GetCircumCircleeFromThreePoints(aurostd::xvector<double>& orig, double& radius, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
bool GetCircleFromTwoPoints(aurostd::xvector<double>& orig, double& radius, const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2);
bool FPositionInsideCell(const aurostd::xvector<double>& r);
bool EmptySphere(const std::deque<_atom>& grid_atoms, const aurostd::xvector<double>& origin_cpos, const double& radius);
bool EmptySphere(const xstructure& str, const aurostd::xvector<double>& origin_cpos, const double& radius);
uint CoordinationPoint(const std::deque<_atom>& atoms, std::deque<_atom>& ratoms, const aurostd::xvector<double>& point, const double& rmin, const double& rmax);
uint CoordinationPoint(const std::deque<_atom>& atoms, std::deque<_atom>& ratoms, const aurostd::xvector<double>& point, const double& rmin);
uint CoordinationPoint(const xstructure& str, std::deque<_atom>& ratoms, const aurostd::xvector<double>& point, const double& rmin, const double& rmax);
uint CoordinationPoint(const xstructure& str, std::deque<_atom>& ratoms, const aurostd::xvector<double>& point, const double& rmin);
bool AddCageToCages(const xstructure& str,
                    const aurostd::xvector<double>& origin_cpos,
                    const aurostd::xvector<double>& origin_fpos,
                    const double& radius,
                    const int& cage_points_type,
                    const double& roughness,
                    std::vector<acage>& cages,
                    std::vector<acage>& cagesX,
                    const bool& osswrite1,
                    std::ostream& oss1,
                    const bool& osswrite2,
                    std::ostream& oss2,
                    int ithread);
uint GetCages4(const xstructure& str, const double& roughness, std::vector<acage>& cages, std::vector<acage>& cages4, const bool& osswrite1, std::ostream& oss1, const bool& osswrite2, std::ostream& oss2);
uint GetCages3(const xstructure& str, const double& roughness, std::vector<acage>& cages, std::vector<acage>& cages3, const bool& osswrite1, std::ostream& oss1, const bool& osswrite2, std::ostream& oss2);
uint GetCages2(const xstructure& str, const double& roughness, std::vector<acage>& cages, std::vector<acage>& cages2, const bool& osswrite1, std::ostream& oss1, const bool& osswrite2, std::ostream& oss2);
bool GetCages(const xstructure& _str,
              _aflags& aflags,
              std::vector<acage>& cagesirreducible,
              std::vector<acage>& cagesreducible,
              std::vector<acage>& cages4,
              std::vector<acage>& cages3,
              std::vector<acage>& cages2,
              double _roughness,
              bool FFFflag,
              bool osswrite,
              std::ostream& oss);

#endif // AFLOW_DEFECTS_H
