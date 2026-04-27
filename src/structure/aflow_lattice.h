
#ifndef AFLOW_LATTICE_H
#define AFLOW_LATTICE_H

#include <istream>
#include <map>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xvector.h"

#include "structure/aflow_xstructure.h"

namespace LATTICE {
  bool lattice_is_working(std::string lat);
  std::string Lattice2TypeAndCentering(const std::string& lattice_type); // DX20191031
  std::string SpaceGroup2Lattice(uint sg);
  std::string SpaceGroup2LatticeTypeAndCentering(uint sg); // DX20191031
  uint Conventional2PrimitiveRatio(char& lattice_centering); // DX20200427
  uint Lattice2SpaceGroup(const std::string& lattice, std::vector<uint>& vsg);
  std::string SpaceGroup2LatticeVariation(uint sg, const xstructure& str);
  std::string ConventionalLattice_SpaceGroup(uint sg, double a, double b, double c);
  std::string ConventionalLattice_SpaceGroup(uint sg, const xstructure& str);
  aurostd::xvector<double> Getabc_angles_Conventional(const aurostd::xmatrix<double>& rlattice, std::string lattice, int mode);
  void findLattices(const std::vector<aurostd::xvector<double>>& translation_vectors,
                    const aurostd::xmatrix<double>& lattice_original,
                    std::vector<aurostd::xmatrix<double>>& lattices,
                    std::vector<aurostd::xmatrix<double>>& lattices_aaa,
                    const std::string& crystal_system,
                    double eps); // DX20210316
  bool fix_sts_sp(xstructure& str_sp, aurostd::xmatrix<double>& rlattice, aurostd::xmatrix<double>& plattice);
  bool Standard_Lattice_Structure(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc, bool full_sym = true);
  bool Standard_Lattice_StructureDefault(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc, bool full_sym = true);
  bool Standard_Lattice_StructureCoarse(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc);
  bool Standard_Lattice_StructureNormal(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc);
  bool Standard_Lattice_StructureMedium(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc);
  bool Standard_Lattice_StructurePrecise(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc);
  bool Standard_Lattice_StructureUltra(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc);
  // bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang); //SC OLD VERSION
  bool Standard_Lattice_Structure(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc, double eps, double epsang, int& time, double symeps);
  bool Standard_Lattice_Structure(const xstructure& str_in, xstructure& str_sp, xstructure& str_sc, double eps, double epsang, int& time, double symeps, bool histogram);
  bool Bravais_Lattice_Structure(xstructure& str_in, xstructure& str_sp, xstructure& str_sc, double eps, double epsang); // calculate everything
  bool Bravais_Lattice_StructureDefault(xstructure& str_in, xstructure& str_sp, xstructure& str_sc, bool full_sym = true); // calculate everything

  // AZ 20231027 START
  /// @brief map of lattice names used in aflow_lattice.cpp
  static const std::map<std::string, std::string> lattice_names = {
      {  "CUB",                                                         "simple cubic"},
      {  "FCC",                                                  "face-centered cubic"},
      {  "BCC",                                                  "body-centered cubic"},
      {  "TET",                                                           "tetragonal"},
      { "BCT1",                                       "body-centered tetragonal c < a"},
      { "BCT2",                                       "body-centered tetragonal a < c"},
      {  "ORC",                                                         "orthorhombic"},
      {"ORCF1",                       "face-centered orthorhombic 1/a^2 > 1/b^2+1/c^2"},
      {"ORCF2",                       "face-centered orthorhombic 1/a^2 < 1/b^2+1/c^2"},
      {"ORCF3",                       "face-centered orthorhombic 1/a^2 = 1/b^2+1/c^2"},
      { "ORCI",                                 "body-centered orthorhombic a < b < c"},
      { "ORCC",                                        "C-centered orthorhombic a < b"},
      {  "HEX",                                                            "hexagonal"},
      { "RHL1",                                              "rhombohedral alpha < 90"},
      { "RHL2",                                              "rhombohedral alpha > 90"},
      {  "MCL",                                                           "monoclinic"},
      {"MCLC1",                                    "C-centered monoclinic kgamma > 90"},
      {"MCLC2",                                    "C-centered monoclinic kgamma = 90"},
      {"MCLC3", "C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2<1"},
      {"MCLC4", "C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2=1"},
      {"MCLC5", "C-centered monoclinic kgamma < 90, bcos(alpha)/c+(bsin(alpha)/a)^2>1"},
      {"TRI1A",                                                            "triclinic"},
      {"TRI1B",                                                            "triclinic"},
      {"TRI2A",                                                            "triclinic"},
      {"TRI2B",                                                            "triclinic"}
  };
  // AZ20231027 END
  bool Lattice(const aurostd::xmatrix<double>& lattice,
               aurostd::xmatrix<double>& lattice_sp,
               aurostd::xmatrix<double>& lattice_sc,
               std::string& bravais_lattice_type,
               std::string& bravais_lattice_variation_type,
               std::string& bravais_lattice_system,
               double eps,
               double epsang);
  std::string Bravais_Lattice_Type(const aurostd::xmatrix<double>& lattice, aurostd::xmatrix<double>& lattice_sp, aurostd::xmatrix<double>& lattice_sc, double eps, double epsang);
  std::string Bravais_Lattice_Type(const aurostd::xmatrix<double>& lattice, aurostd::xmatrix<double>& lattice_sp, aurostd::xmatrix<double>& lattice_sc);
  std::string Bravais_Lattice_Type(const aurostd::xmatrix<double>& lattice, double eps, double epsang);
  std::string Bravais_Lattice_Type(const aurostd::xmatrix<double>& lattice);
  std::string Bravais_Lattice_Variation_Type(const aurostd::xmatrix<double>& lattice, aurostd::xmatrix<double>& lattice_sp, aurostd::xmatrix<double>& lattice_sc, double eps, double epsang);
  std::string Bravais_Lattice_Variation_Type(const aurostd::xmatrix<double>& lattice, aurostd::xmatrix<double>& lattice_sp, aurostd::xmatrix<double>& lattice_sc);
  std::string Bravais_Lattice_Variation_Type(const aurostd::xmatrix<double>& lattice, double eps, double epsang);
  std::string Bravais_Lattice_Variation_Type(const aurostd::xmatrix<double>& lattice);
  std::string Bravais_Lattice_System(const aurostd::xmatrix<double>& lattice, aurostd::xmatrix<double>& lattice_sp, aurostd::xmatrix<double>& lattice_sc, double eps, double epsang);
  std::string Bravais_Lattice_System(const aurostd::xmatrix<double>& lattice, aurostd::xmatrix<double>& lattice_sp, aurostd::xmatrix<double>& lattice_sc);
  std::string Bravais_Lattice_System(const aurostd::xmatrix<double>& lattice, double eps, double epsang);
  std::string Bravais_Lattice_System(const aurostd::xmatrix<double>& lattice);
  std::string Primitive_Lattice_Type(const xstructure& str);
  std::string Bravais_Lattice_System(const xstructure& str);
  std::string Conventional_Lattice_Type(const xstructure& str);
  xstructure Standard_Primitive_Lattice_Structure(const xstructure& str);
  xstructure Standard_Conventional_Lattice_Structure(const xstructure& str);
  std::string Get_Primitive_Lattice_Structure(const xstructure& str);
  aurostd::xmatrix<double> sc2sp(const aurostd::xmatrix<double>& rlattice, std::string lattice, bool inverseflag);
  aurostd::xmatrix<double> sp2sc(const aurostd::xmatrix<double>& rlattice, std::string lattice, bool inverseflag);
  void BZPLOTDATA(std::string options, std::istream& poscar, int mode);
} // namespace LATTICE
aurostd::xvector<double> Vrotate(aurostd::xvector<double> v, aurostd::xvector<double> axisrot, double theta);
void CheckLatticeHistogram();
namespace LATTICE {
  // kpoints and brillouin zones
  std::string KPOINTS_Directions(xstructure str_in, double grid, bool& foundBZ);
  std::string KPOINTS_Directions(std::string lattice_type, aurostd::xmatrix<double> sp, double _grid, int iomode, bool& foundBZ);
  std::string KPOINTS_Directions(std::string lattice_type, aurostd::xmatrix<double> sp, aurostd::xmatrix<double> transformation_matrix, double _grid, int iomode, bool& foundBZ); // DX20181101
} // namespace LATTICE

#endif // AFLOW_LATTICE_H
