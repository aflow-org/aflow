
#ifndef AFLOW_SYMMETRY_H
#define AFLOW_SYMMETRY_H

#include <deque>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xcomplex.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xvector.h"

#include "flow/aflow_xclasses.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

namespace SYM {
  // DX+CO START
  bool ApplyAtomValidate(const _atom& atom_in, _atom& atom_out, const _sym_op& symop, const xstructure& a); // DX+CO
  bool ApplyAtomValidate(const _atom& atom_in, _atom& atom_out, const _sym_op& symop, const xstructure& a, bool _incell_, bool roff); // DX+CO
  bool ApplyAtomValidate(const _atom& atom_in, _atom& atom_out, const _sym_op& symop, const xstructure& a, bool skew, bool _incell_, bool roff, double _eps_); // DX+CO
  bool ApplyAtomValidate(
      const _atom& atom_in, _atom& atom_out, const _sym_op& symop, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& c2f, const aurostd::xmatrix<double>& f2c, bool skew, bool _incell_, bool roff, double _eps_); // DX+CO
  _atom ApplyAtom(const _atom& atom_in, const _sym_op& symop, const xstructure& str); // DX
  // DX+CO END
  _atom ApplyAtom(const _atom& atom_in, const _sym_op& symop, const xstructure& str, bool _incell_);
  // DX+CO START
  _atom ApplyAtom(const _atom& atom_in, const _sym_op& symop, const xstructure& str, bool _incell_, bool roff); // DX
  _atom ApplyAtom(const _atom& atom_in, const _sym_op& symop, const xstructure& str, bool _incell_, bool roff, bool validatePosition); // DX
  _atom ApplyAtom(const _atom& atom_in, const _sym_op& symop, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& c2f, const aurostd::xmatrix<double>& f2c, bool skew, bool _incell_, bool roff, bool validatePosition, double eps); // DX
  // DX+CO END
  aurostd::xvector<double> ApplyCpos(const aurostd::xvector<double>& cpos_in, const _sym_op& symop, const xstructure& str, bool _incell_);
  aurostd::xvector<double> ApplyCpos(const aurostd::xvector<double>& cpos_in, const _sym_op& symop, const xstructure& str);
  aurostd::xvector<double> ApplyFpos(const aurostd::xvector<double>& fpos_in, const _sym_op& symop, const xstructure& str, bool _incell_);
  aurostd::xvector<double> ApplyFpos(const aurostd::xvector<double>& fpos_in, const _sym_op& symop, const xstructure& str);
  aurostd::xvector<int> ApplyIJK(const aurostd::xvector<int>& ijk_in, const _sym_op& symop, const xstructure& str);
  int ApplyL(const int& l_in, const _sym_op& symop, const xstructure& str);
  xstructure ApplyXstructure(const _sym_op& symop, const xstructure& str);
  xstructure ApplyXstructure(const _sym_op& symop, const xstructure& str, bool _incell_);
  // DX+CO START
  bool AtomsEquivalent(xstructure& str, _atom& atom1, _atom& atom2); // DX
  bool AtomsEquivalent(xstructure& str, _atom& atom1, _atom& atom2, double& eps); // DX
  bool AtomsEquivalent(xstructure& str, _atom& a, _atom& b, bool skew, double tol); // DX //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool AtomsEquivalent_Basis(xstructure& str, int atom1_indx, int atom2_indx);
  // DX+CO END
  bool CposEquivalent(const xstructure& str, const aurostd::xvector<double>& cpos1, const aurostd::xvector<double>& cpos2, const double& eps);
  bool FposEquivalent(const xstructure& str, const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2, const double& eps);
  bool CposEquivalent(const xstructure& str, const aurostd::xvector<double>& cpos1, const aurostd::xvector<double>& cpos2);
  bool FposEquivalent(const xstructure& str, const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2);
  bool TypePointGroupOperation(const aurostd::xmatrix<double>& Uc,
                               const aurostd::xmatrix<double>& Uf,
                               std::string& _string,
                               bool& _inversion,
                               double& _angle,
                               aurostd::xvector<double>& _axis,
                               aurostd::xmatrix<double>& _generator,
                               aurostd::xvector<double>& _generator_coefficients,
                               aurostd::xmatrix<aurostd::xcomplex<double>>& _SU2_matrix,
                               aurostd::xvector<aurostd::xcomplex<double>>& _su2_coefficients,
                               double _eps_); // calculate the symmetry inversion,type,axis,generator //DX20171206 - Added generator coefficients //DX20171207 - Added Uf //DX20180117 - Added SU2 and su2 coefficients
  bool TypePointGroupOperationInternational(const aurostd::xmatrix<double>& Uc,
                                            std::string& _stringHM,
                                            std::string& _stringSC,
                                            const bool& _inversion,
                                            const double& _angle,
                                            const aurostd::xvector<double>& _axis,
                                            const aurostd::xmatrix<double>& _generator,
                                            aurostd::xvector<double>& _generator_coefficients,
                                            aurostd::xmatrix<aurostd::xcomplex<double>>& _SU2_matrix,
                                            aurostd::xvector<aurostd::xcomplex<double>>& _su2_coefficients,
                                            double _eps_); // International symbol = Hermann-Mauguin notation & Schonflies notation //DX20171206 - Added generator coefficients //DX20180117 - Added SU2 and su2 coefficients
  // DX+CO START
  uint AddSymmetryToStructure(xstructure& a,
                              const uint& iat,
                              const aurostd::xmatrix<double>& Uc,
                              const aurostd::xmatrix<double>& Uf,
                              const aurostd::xvector<double>& ctau,
                              const aurostd::xvector<double>& ftau,
                              const aurostd::xvector<double>& ctrasl,
                              const aurostd::xvector<double>& ftrasl,
                              const std::vector<int>& basis_atoms_map,
                              const std::vector<int>& basis_types_map,
                              bool basis_map_calculated,
                              char group);
  uint AddSymmetryToStructure(xstructure& a,
                              const uint& iat,
                              const aurostd::xmatrix<double>& Uc,
                              const aurostd::xmatrix<double>& Uf,
                              const aurostd::xvector<double>& ctau,
                              const aurostd::xvector<double>& ftau,
                              const aurostd::xvector<double>& ctrasl,
                              const aurostd::xvector<double>& ftrasl,
                              const std::vector<int>& basis_atoms_map,
                              const std::vector<int>& basis_types_map,
                              bool basis_map_calculated,
                              char group,
                              bool roff); // DX
  uint AddSymmetryToStructure(xstructure& a,
                              const aurostd::xmatrix<double>& Uc,
                              const aurostd::xmatrix<double>& Uf,
                              const aurostd::xvector<double>& ctau,
                              const aurostd::xvector<double>& ftau,
                              const aurostd::xvector<double>& ctrasl,
                              const aurostd::xvector<double>& ftrasl,
                              const std::vector<int>& basis_atoms_map,
                              const std::vector<int>& basis_types_map,
                              bool basis_map_calculated,
                              char group);
  uint AddSymmetryToStructure(xstructure& a,
                              const aurostd::xmatrix<double>& Uc,
                              const aurostd::xmatrix<double>& Uf,
                              const aurostd::xvector<double>& ctau,
                              const aurostd::xvector<double>& ftau,
                              const aurostd::xvector<double>& ctrasl,
                              const aurostd::xvector<double>& ftrasl,
                              const std::vector<int>& basis_atoms_map,
                              const std::vector<int>& basis_types_map,
                              bool basis_map_calculated,
                              char group,
                              bool roff); // DX
  uint AddSymmetryToStructure(
      xstructure& a, const uint& iat, const aurostd::xmatrix<double>& Uc, const aurostd::xmatrix<double>& Uf, const std::vector<int>& basis_atoms_map, const std::vector<int>& basis_types_map, bool basis_map_calculated, char group);
  uint AddSymmetryToStructure(xstructure& a,
                              const uint& iat,
                              const aurostd::xmatrix<double>& Uc,
                              const aurostd::xmatrix<double>& Uf,
                              const std::vector<int>& basis_atoms_map,
                              const std::vector<int>& basis_types_map,
                              bool basis_map_calculated,
                              char group,
                              bool roff); // DX
  uint AddSymmetryToStructure(xstructure& a, const aurostd::xmatrix<double>& Uc, const aurostd::xmatrix<double>& Uf, const std::vector<int>& basis_atoms_map, const std::vector<int>& basis_types_map, bool basis_map_calculated, char group); // DX
  uint AddSymmetryToStructure(
      xstructure& a, const aurostd::xmatrix<double>& Uc, const aurostd::xmatrix<double>& Uf, const std::vector<int>& basis_atoms_map, const std::vector<int>& basis_types_map, bool basis_map_calculated, char group, bool roff); // DX
  bool PointGroupsIdentical(const std::vector<_sym_op>& vpg1, const std::vector<_sym_op>& vpg2, double eps, bool is_same_lattice = false); // DX20171207 - added is_same_lattice
  // GG START
  bool CalculateQuaternion(_sym_op& a);
  // GG STOP
  bool ComplexSU2Rotations(aurostd::xmatrix<aurostd::xcomplex<double>>& _SU2_matrix, aurostd::xvector<aurostd::xcomplex<double>>& _su2_coefficients, double& theta, aurostd::xvector<double>& _axis); // DX20180117 - add SU(2) and su(2) coefficients
  // DX+CO END
  bool CalculatePointGroup(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const aurostd::xmatrix<double>& lattice, std::vector<_sym_op> pgroup, std::ofstream& FileMESSAGE, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_);
  uint CalculatePointGroup(const aurostd::xmatrix<double>& lattice, std::vector<_sym_op> pgroup, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_); // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const aurostd::xmatrix<double>& lattice, double _eps_); // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const aurostd::xmatrix<double>& lattice); // POINT GROUP      _PGROUP_
  bool CalculatePointGroup(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_, std::string format = "txt"); // POINT GROUP _PGROUP_
  bool CalculatePointGroupKLattice(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // POINT GROUP KLATTICE _PGROUPK_
  bool CalculatePointGroupKCrystal(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // POINT GROUP KCRYSTAL _PGROUPK_XTAL_ //DX20171205 - New group: reciprocal space counterpart of pgroup_xtal
  bool TransformSymmetryFromRealToReciprocal(std::ofstream& FileMESSAGE, xstructure& real_space_crystal, xstructure& reciprocal_space, _aflags& aflags, const bool& osswrite, std::ostream& oss, std::string& pgroup_type); // DX20170808 - New klattice routine //DX20171205 - Added pgroup_type option to account for pgroupk_xtal
  bool CalculateSitePointGroup(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // SITE POINT GROUP _AGROUP_
  // DX+CO START
  bool CalculateSitePointGroup(std::ofstream& FileMESSAGE, xstructure& a, int CALCULATION_MODE, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // SITE POINT GROUP _AGROUP_
  bool CalculateSitePointGroup_EquivalentSites(xstructure& a, double _eps_); // DX
  bool CalculateSitePointGroup_EquivalentSites(xstructure& a, bool get_full_basis, double _eps_); // DX
  // DX+CO END
  bool CalculatePointGroupCrystal(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // POINT GROUP      _PGROUP_
  bool CalculatePointGroupCrystal(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_, std::string format = "txt"); // POINT GROUP _PGROUP_
  // DX+CO START
  bool CalculatePointGroupKPatterson(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // POINT GROUP PATTERSON _PGROUPK_PATTERSON_ //DX20200129
  bool CalculatePointGroupKPatterson(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_, std::string format = "txt"); // POINT GROUP PATTERSON     _PGROUPK_PATTERSON_ //DX20200129
  bool PointGroupMap(xstructure& a, std::string& pgname, std::string& operations, char group); // DX20170906
  bool PointGroupLookUpTable(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format);
  // DX+CO END
  void CalculateSitePointGroup2(xstructure& a, bool ComMidss); // for --agroup2 and --agroup2m
  // DX START
  // xstructure and _sym_op
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp, bool map_types, std::vector<int>& basis_atoms_map, std::vector<int>& basis_types_map);
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp, bool map_types, double tolerance, std::vector<int>& basis_atoms_map, std::vector<int>& basis_types_map); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp, bool map_types, bool skew, double tolerance, std::vector<int>& basis_atoms_map, std::vector<int>& basis_types_map); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp, bool map_types, bool skew, double tolerance, std::vector<int>& basis_atoms_map, std::vector<int>& basis_types_map); // CO20190520 - removed pointers for bools and doubles, added const where possible
  // atoms, c2f, f2c and _sym_op
  bool getFullSymBasis(const std::deque<_atom>& atoms,
                       const aurostd::xmatrix<double>& lattice,
                       const aurostd::xmatrix<double>& c2f,
                       const aurostd::xmatrix<double>& f2c,
                       _sym_op& symOp,
                       bool map_types,
                       bool skew,
                       double tolerance,
                       std::vector<int>& basis_atoms_map,
                       std::vector<int>& basis_types_map); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis(const std::deque<_atom>& atoms,
                       const aurostd::xmatrix<double>& lattice,
                       const aurostd::xmatrix<double>& c2f,
                       const aurostd::xmatrix<double>& f2c,
                       _sym_op& symOp,
                       bool map_types,
                       bool skew,
                       double tolerance,
                       std::vector<int>& basis_atoms_map,
                       std::vector<int>& basis_types_map); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool CalculateFactorGroup(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // FACTOR GROUP     _FGROUP_
  bool CalculateFactorGroup(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_, std::string format = "txt"); // FACTOR GROUP _FGROUP_
  // DX START
  bool AtomsMapped(const _atom& a, const _atom& b, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // DX20190620
  bool AtomsMapped(const _atom& a, const _atom& b, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  aurostd::xvector<double> minimizeDistanceCartesianMethod(const aurostd::xvector<double>& cpos1, const aurostd::xvector<double>& cpos2, const aurostd::xmatrix<double>& lattice); // DX20190613
  aurostd::xvector<double> minimizeDistanceCartesianMethod(const aurostd::xvector<double>& cpos1, const aurostd::xvector<double>& cpos2, const aurostd::xmatrix<double>& lattice, aurostd::xvector<int>& ijk); // DX20190613
  aurostd::xvector<double> minimizeDistanceFractionalMethod(const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2); // DX20190613
  aurostd::xvector<double> minimizeDistanceFractionalMethod(const aurostd::xvector<double>& fdiff); // DX20190613
  aurostd::xvector<double> minimizeDistanceFractionalMethod(const aurostd::xvector<double>& fdiff, aurostd::xvector<int>& ijk); // DX20190613
  // DX20190613 [OBOSLETE] bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, const xmatrix<double>& c2f, const xmatrix<double>& f2c, double tol);
  // //CO20190520 - removed pointers for bools and doubles, added const where possible DX20190613 [OBOSLETE] bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2,
  // xvector<double>& out, const xmatrix<double>& c2f, const xmatrix<double>& f2c, xvector<int>& ijk, bool& restriction, double tol); //CO20190520 - removed pointers for bools and doubles, added const where
  // possible DX20190613 [OBOSLETE] double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice); DX20190613 [OBOSLETE] double
  // minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice,xvector<double>& min_vec,xvector<int>& ijk); DX20190613 [OBOSLETE] xvector<double>
  // minimumCartesianVector(const xvector<double>&, const xvector<double>&, DX20190613 [OBOSLETE]                                        const xmatrix<double>&);  //ME20180730 DX20190613 [OBOSLETE]
  // xvector<double> minimumCartesianVector(const xvector<double>&, const xvector<double>&, DX20190613 [OBOSLETE]                                        const xmatrix<double>&, xvector<int>&);  //ME20180730
  // DX20190613 [OBOSLETE] bool PBC(xvector<double>& v_in, xvector<int>& ijk, bool& restriction); DX20190613 [OBOSLETE] bool PBC(xvector<double>& v_in);
  aurostd::xvector<double> FPOSDistFromFPOS(const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2, const aurostd::xmatrix<double>& lattice, bool skew = false); // CO20190525
  aurostd::xvector<double> FPOSDistFromFPOS(
      const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& c2f, const aurostd::xmatrix<double>& f2c, bool skew = false); // CO20190525
  aurostd::xvector<double> CPOSDistFromFPOS(const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2, const aurostd::xmatrix<double>& lattice, bool skew = false); // DX20190620
  aurostd::xvector<double> CPOSDistFromFPOS(const aurostd::xvector<double>& fpos1, const aurostd::xvector<double>& fpos2, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew = false); // DX20190620

  bool FPOSMatch(const std::deque<_atom>& atom_set, const _atom& atom2, uint& match_type, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool FPOSMatch(const std::deque<_atom>& atom_set, const _atom& atom2, uint& match_type, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // DX20190620 - overload
  bool FPOSMatch(const _atom& atom1, const _atom& atom2, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool FPOSMatch(const _atom& atom1, const _atom& atom2, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // DX20190620 - overload
  bool FPOSMatch(const aurostd::xvector<double>& atom1, const aurostd::xvector<double>& atom2, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool FPOSMatch(const aurostd::xvector<double>& atom1, const aurostd::xvector<double>& atom2, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // DX20190620 - overload
  bool validateAtomPosition(const _atom& atom, const aurostd::xmatrix<double>& c2f, const aurostd::xmatrix<double>& f2c, bool skew, double& _eps_); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool validateAtomPosition(const aurostd::xvector<double>& cpos, const aurostd::xvector<double>& fpos, const aurostd::xmatrix<double>& c2f, const aurostd::xmatrix<double>& f2c, bool skew, double& _eps_); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool MapAtom(const std::deque<_atom>& a_vec, const _atom& b, bool map_types, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // DX20190620
  bool MapAtom(const std::deque<_atom>& a_deq, const _atom& b, bool map_types, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool MapAtom(const std::vector<_atom> a_vec, const _atom b, bool map_types, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // DX20190620
  bool MapAtom(const std::vector<_atom> a_vec, const _atom b, bool map_types, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool MapAtom(const aurostd::xvector<double>& a, const aurostd::xvector<double>& b, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // DX20190620
  bool MapAtom(const aurostd::xvector<double>& a, const aurostd::xvector<double>& b, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool MapAtom(const _atom& a, const _atom& b, bool map_types, const aurostd::xmatrix<double>& lattice, bool skew, double tol); // DX20190620
  bool MapAtom(const _atom& a, const _atom& b, bool map_types, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& f2c, bool skew, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
  bool BringInCellTolABC(xstructure& a, aurostd::xvector<double> tol_abc_res);
  //[CO20190515 - not needed and is ambiguous with overload]bool MapAtomWithBasis(vector<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool skew, double tol,bool fast=true);    //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool MapAtomWithBasis(const std::vector<_atom>& vec, const _atom& a, bool map_types, std::deque<uint>& index_to_check, const aurostd::xmatrix<double>& lattice, bool skew, double tol, uint& mapped_index, bool fast = true); // DX20190620
  bool MapAtomWithBasis(const std::vector<_atom>& vec,
                        const _atom& a,
                        bool map_types,
                        std::deque<uint>& index_to_check,
                        const aurostd::xmatrix<double>& lattice,
                        const aurostd::xmatrix<double>& f2c,
                        bool skew,
                        double tol,
                        uint& mapped_index,
                        bool fast = true); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  //[CO20190515 - not needed and is ambiguous with overload]bool MapAtomWithBasis(deque<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool skew, double tol,bool fast=true); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool MapAtomWithBasis(const std::deque<_atom>& vec, const _atom& a, bool map_types, std::deque<uint>& index_to_check, const aurostd::xmatrix<double>& lattice, bool skew, double tol, uint& mapped_index, bool fast = true); // DX20190620
  bool MapAtomWithBasis(const std::deque<_atom>& vec,
                        const _atom& a,
                        bool map_types,
                        std::deque<uint>& index_to_check,
                        const aurostd::xmatrix<double>& lattice,
                        const aurostd::xmatrix<double>& f2c,
                        bool skew,
                        double tol,
                        uint& mapped_index,
                        bool fast = true); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
  bool isLatticeSkewed(const aurostd::xmatrix<double>& lattice, double& min_dist, double tol); // CO20190520 - removed pointers for bools and doubles, added const where possible
  double minimumDistance(const xstructure& xstr);
  double minimumDistance(const std::deque<_atom>& atoms); // CO20190808 - for NON periodic systems
  double minimumDistance(const std::deque<_atom>& atoms, const aurostd::xmatrix<double>& lattice, double scale = 1.0);
  double defaultTolerance(const xstructure& xstr);
  bool checkAngle(aurostd::xvector<double>& v1, aurostd::xvector<double>& v2, double input_angle, double tolerance); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool checkAngle(aurostd::xvector<double>& v1, aurostd::xvector<double>& v2, double input_angle, bool& is_deg, double tolerance); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, double tolerance); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, bool& is_deg, double tolerance); // CO20190520 - removed pointers for bools and doubles, added const where possible
  bool change_tolerance(xstructure& xstr, double& tolerance, double& min_dist, bool& no_scan); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190524 - need pointer for tolerance, otherwise it will not update
  std::deque<std::deque<_atom>> break_up_by_type(std::deque<_atom>& expanded_crystal);
  std::vector<std::vector<_atom>> break_up_by_type(std::vector<_atom> expanded_crystal);
  bool CheckForIdentity(const xstructure& xstr); // DX
  bool checkSuperCellLatticePoints(xstructure& xstr, int& num_lattice_points, char& centering, uint& expand_size); // DX
  bool ComparePointGroupAndSpaceGroupString(xstructure& xstr, int& multiplicity_of_primitive, bool& derivative_structure); // DX
  // DX END
  bool CalculateSpaceGroup(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // SPACE GROUP      _SGROUP_

  bool CalculateInequivalentAtoms(xstructure&);
  bool CalculateInequivalentAtoms(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // EQUIVALENT ATOMS _IATOMS_
  // DX+CO START
  bool CalculateInequivalentAtoms(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_, std::string format = "txt"); // EQUIVALENT ATOMS _IATOMS_ //DX
  bool CalculateInequivalentAtoms(std::ofstream& FileMESSAGE, xstructure& a, bool rely_on_basis, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, std::string format = "txt"); // EQUIVALENT ATOMS _IATOMS_ //DX
  bool CalculateInequivalentAtoms(std::ofstream& FileMESSAGE, xstructure& a, bool rely_on_basis, _aflags& aflags, bool _write_, const bool& osswrite, std::ostream& oss, double _eps_, std::string format = "txt"); // EQUIVALENT ATOMS _IATOMS_ //DX
  // DX+CO END

  void writePythonScript(std::ostream& oss); // DX20201228
} // namespace SYM
std::string AgroupSymmetryToJson(std::vector<std::vector<_sym_op>>& group, char& mode); // DX20170803 - For Python wrapper
std::string EquivalentAtomsToJson(std::vector<std::vector<int>>& iatoms); // DX20170803 - For Python wrapper
std::string SymmetryToJson(std::vector<_sym_op>& group, char& mode); // DX20170803 - For Python wrapper
bool KBIN_SymmetryWrite(std::ofstream& FileMESSAGE, xstructure& a, _aflags& aflags, char group, const bool& osswrite, std::ostream& oss, const std::string& format = "txt");
// bool KBIN_SymmetryToScreen(xstructure& a, std::string& format, std::ostream& oss); //DX20170803 - For Python wrapper
bool KBIN_SymmetryToScreen(xstructure& a, const std::string& format, std::ostream& oss, char mode = '\0'); // DX20170822 - For Python wrapper
bool KBIN_SymmetryToScreenWeb(xstructure& a, std::ostream& oss, char mode); // ME20210402
bool KBIN_StepSymmetryPerform(xstructure& a, std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, const bool& osswrite, std::ostream& oss);
std::vector<double> PointGroupHistogramCheck(xstructure& a);
// bool SYM_CalculatePointGroup(std::ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,std::ostream& oss);      // MIKNOWSKI BASIS REDUCTION

#endif // AFLOW_SYMMETRY_H
