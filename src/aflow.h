// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_H_
#define _AFLOW_H_

#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <istream>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xcomplex.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xtensor.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"
#include "flow/aflow_xclasses.h"
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

namespace aflowlib // forward declaration
{
  class _aflowlib_entry;
}

// CO20200731 START
static const std::string SEP_TAG1 = ":";
static const std::string SEP_TAG2 = "_";
static const std::string TAG_POCC = "POCC";
static const std::string TAG_TOL = "TOL";
static const std::string TAG_ARUN = "ARUN";
static const std::string TAG_TITLE_POCC = SEP_TAG1 + TAG_POCC + SEP_TAG2;
static const std::string TAG_TITLE_POCC_TOL = SEP_TAG1 + TAG_TOL + SEP_TAG2;
static const std::string TAG_TITLE_ARUN = SEP_TAG1 + TAG_ARUN + ".";
static const std::string TAG_TITLE_POCC_ARUN = TAG_TITLE_ARUN + TAG_POCC + SEP_TAG2;
static const std::string POCC_DOSCAR_PREFIX = "DOSCAR.pocc_T";
static const std::string POCC_PHDOSCAR_PREFIX = "PHDOSCAR.pocc_T";  // ME20210927
// CO20200731 END
static const std::vector<std::string> BRAVAIS_LATTICES = {"BCC", "BCT", "CUB", "FCC", "HEX", "MCL", "MCLC", "ORC", "ORCC", "ORCF", "ORCI", "RHL", "TET", "TRI"}; // HE20220420

extern std::string _AFLOWIN_;
extern std::string _AFLOWLOCK_;

const std::string VASP_KEYWORD_EXECUTION = " Executing: ";

// --------------------------------------------------------------------------
// definitions for projects
extern uint LIBRARY_AUID, LIBRARY_ICSD, LIBRARY_LIB0, LIBRARY_LIB1, LIBRARY_LIB2, LIBRARY_LIB3, LIBRARY_LIB4, LIBRARY_LIB5, LIBRARY_LIB6, LIBRARY_LIB7, LIBRARY_LIB8, LIBRARY_LIB9; // not in order.. will be nailed by init.cpp

// --------------------------------------------------------------------------
// aflow_arguments
uint PflowARGs(std::vector<std::string>& argv, std::vector<std::string>& cmds, aurostd::xoption& vpflow); // called inside Init::InitMachine coded in aflow_pflow_main.cpp

namespace soliquidy {
  void run_cmd(aurostd::xoption& vpflow, std::istream& input);
  aurostd::JSON::object CalculateEntryLoader(const std::shared_ptr<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>& aflux_result, const std::string& result_folder = "", bool create_x3d = false);
  aurostd::JSON::object CalculateWorklist(const std::vector<std::string>& worklist, const std::string& result_folder = "", bool create_x3d = false);
  aurostd::JSON::object Calculate(const xstructure& work_structure, const std::filesystem::path& out_folder = "", bool create_x3d = false);
  template <class utype> aurostd::JSON::object CalculateAUID(const utype& auid, const std::string& result_folder = "", bool create_x3d = false);
} // namespace soliquidy

std::vector<uint> getAtomIndicesByType(const xstructure& xstr, int type); // DX20210322
std::vector<uint> getAtomIndicesByName(const xstructure& xstr, const std::string& name); // DX20210322
std::vector<uint> getLeastFrequentAtomTypes(const xstructure& xstr); // DX20210322
std::vector<std::string> getLeastFrequentAtomSpecies(const xstructure& xstr, bool clean = true); // DX20201230 - moved from XtalFinder

xstructure GetStructure(const int& iomode, std::ifstream& input); // plug from cin
xstructure GetStructure(const int& iomode, const std::string& Directory); // plug from a directory
xstructure SetSDNumbers(const xstructure& a, const std::vector<std::string>& in_sd);
xstructure SetSDTypes(const xstructure& a, const std::vector<std::string>& in_sd);
std::vector<int> GetTypes(const xstructure& a);
std::vector<std::string> GetNames(const xstructure& a);
std::vector<std::string> GetCleanNames(const xstructure& a);
std::vector<double> GetSpins(const xstructure& a);
std::string GetElementName(std::string stringin);
std::string GetSpaceGroupName(int spacegroupnumber, const std::string& directory = ""); // DX20180526 - add directory
int GetSpaceGroupNumber(const std::string& spacegroupsymbol, const std::string& directory = ""); // DX20190708
std::string GetSpaceGroupLabel(int spacegroupnumber);
std::string GetSpaceGroupSchoenflies(int spacegroupnumber, const std::string& directory = ""); // DX20170901 //DX20180526 - add directory
std::string GetSpaceGroupHall(int spacegroupnumber, int setting = SG_SETTING_1, const std::string& directory = ""); // DX20170901 //DX20180526 - add directory //DX20180806 - added setting
std::string GetLaueLabel(std::string& point_group); // DX20170901 //DX20180526 - add directory

aurostd::xmatrix<double> MetricTensor(const xstructure& a); // CO20180409
aurostd::xmatrix<double> MetricTensor(const aurostd::xmatrix<double>& lattice, double scale = 1.0); // CO20180409
aurostd::xmatrix<double> ReciprocalLattice(const xstructure& a); // CO20180409
aurostd::xmatrix<double> ReciprocalLattice(const aurostd::xmatrix<double>& rlattice, double scale = 1.0); // CO20180409
std::string KPPRA(int& k1, int& k2, int& k3, const aurostd::xmatrix<double>& rlattice, const int& NK);
std::string KPPRA(xstructure& str, const int& _NK);
std::string KPPRA_DELTA(int& k1, int& k2, int& k3, const aurostd::xmatrix<double>& rlattice, const double& DK);
std::string KPPRA_DELTA(xstructure& str, const double& DK);
[[deprecated("Use GetNBANDS_VASP_SERIAL for a more appropriate estimate.")]]
int GetNBANDS_AFLOW3(int electrons, int nions, int spineach, bool ispin = true, int NPAR = 1); // CO20210315 - spin==true is safer, added NPAR
int GetNBANDS_VASP_SERIAL_spineach(int electrons, int nions, int spineach, bool ispin);
int GetNBANDS_VASP_SERIAL(int electrons, int nions, const std::vector<int>& spins, bool ispin);
int GetNBANDS_VASP(int electrons, int nions, int spineach, bool ispin, int NPAR = 1);
double GetZVAL(const std::stringstream& sss, std::vector<double>& vZVAL);
double GetZVAL(const _xvasp& xvasp, std::vector<double>& vZVAL);
double GetZVAL(const std::string& directory, std::vector<double>& vZVAL);
double GetCellAtomZVAL(const std::stringstream& sss, std::vector<double>& vZVAL, const std::stringstream& sstr, std::vector<double>& sZVAL, std::string mode); // sss sstr returns ZVAL cell, VAL and sZVAL
double GetCellAtomZVAL(const std::string& directory, std::vector<double>& vZVAL, std::vector<double>& sZVAL, std::string mode); // from directory POT/POS returns total ZVAL cell, vZVAL and sZVAL
double GetPOMASS(const std::stringstream& sss, std::vector<double>& vPOMASS);
double GetPOMASS(const _xvasp& xvasp, std::vector<double>& vPOMASS);
double GetPOMASS(const std::string& directory, std::vector<double>& vPOMASS);
double GetCellAtomPOMASS(const std::stringstream& sss, std::vector<double>& vPOMASS, const std::stringstream& sstr, std::vector<double>& sPOMASS, std::string mode); // sss sstr returns POMASS cell, VAL and sPOMASS
double GetCellAtomPOMASS(const std::string& directory, std::vector<double>& vPOMASS, std::vector<double>& sPOMASS, std::string mode); // from directory POT/POS returns total POMASS cell, vPOMASS and sPOMASS
double GetVol(const aurostd::xmatrix<double>& lat);
double det(const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
double GetVol(const aurostd::xvector<double>& v1, const aurostd::xvector<double>& v2, const aurostd::xvector<double>& v3);
double det(const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
aurostd::xvector<double> Getabc_angles(const aurostd::xmatrix<double>& lat, int mode);
aurostd::xvector<long double> Getabc_angles(const aurostd::xmatrix<long double>& lat, int mode);
aurostd::xvector<double> Getabc_angles(const aurostd::xmatrix<double>& lat, const aurostd::xvector<int>& permut, int mode);
aurostd::xvector<double> Getabc_angles(const aurostd::xvector<double>& r1, const aurostd::xvector<double>& r2, const aurostd::xvector<double>& r3, int mode);
aurostd::xvector<double> Getabc_angles(const aurostd::xvector<double>& r1, const aurostd::xvector<double>& r2, const aurostd::xvector<double>& r3, const aurostd::xvector<int>& permut, int mode);
#define _Getabc_angles Getabc_angles
// #define _Getabc_angles __NO_USE_Sortabc_angles
aurostd::xvector<double> Sortabc_angles(const aurostd::xmatrix<double>& lat, int mode);
aurostd::xmatrix<double> GetClat(const aurostd::xvector<double>& abc_angles);
aurostd::xmatrix<double> GetClat(const double& a, const double& b, const double& c, const double& alpha, const double& beta, const double& gamma);
xstructure GetIntpolStr(xstructure strA, xstructure strB, const double& f, const std::string& path_flag);
double RadiusSphereLattice(const aurostd::xmatrix<double>& lattice, double scale = 1.0); // CO20180409
aurostd::xvector<int> LatticeDimensionSphere(const aurostd::xmatrix<double>& lattice, double radius, double scale = 1.0); // CO20180409
aurostd::xvector<int> LatticeDimensionSphere(const xstructure& str, double radius);
void resetLatticeDimensions(const aurostd::xmatrix<double>& lattice,
                            double radius,
                            aurostd::xvector<int>& dims,
                            std::vector<aurostd::xvector<double>>& l1,
                            std::vector<aurostd::xvector<double>>& l2,
                            std::vector<aurostd::xvector<double>>& l3,
                            std::vector<int>& a_index,
                            std::vector<int>& b_index,
                            std::vector<int>& c_index); // DX20191122
aurostd::xvector<double> F2C(const double& scale, const aurostd::xmatrix<double>& lattice, const aurostd::xvector<double>& fpos); // fpos are F components per COLUMS !
aurostd::xvector<double> F2C(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<double>& fpos); // fpos are F components per COLUMS !
aurostd::xvector<double> C2F(const double& scale, const aurostd::xmatrix<double>& lattice, const aurostd::xvector<double>& cpos); // cpos are C components per COLUMS !
aurostd::xvector<double> C2F(const aurostd::xmatrix<double>& lattice, const aurostd::xvector<double>& cpos); // cpos are C components per COLUMS !
aurostd::xmatrix<double> F2C(const double& scale, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& fpos); // fpos are F components per COLUMS !
aurostd::xmatrix<double> F2C(const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& fpos); // fpos are F components per COLUMS !
aurostd::xmatrix<double> C2F(const double& scale, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& cpos); // cpos are C components per COLUMS !
aurostd::xmatrix<double> C2F(const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& cpos); // cpos are C components per COLUMS !
_atom F2C(const double& scale, const aurostd::xmatrix<double>& lattice, const _atom& iatom); // atom.fpos are F components per COLUMS !
_atom F2C(const xstructure& str, const _atom& iatom); // atom.fpos are F components per COLUMS !
_atom C2F(const double& scale, const aurostd::xmatrix<double>& lattice, const _atom& iatom); // atom.cpos are C components per COLUMS !
_atom C2F(const aurostd::xmatrix<double>& lattice, const _atom& iatom); // atom.cpos are C components per COLUMS !
_atom F2C(const double& scale, const xstructure& str, const _atom& iatom); // atom.fpos are F components per COLUMS !
_atom F2C(const xstructure& str, const _atom& iatom); // atom.fpos are F components per COLUMS !
_atom C2F(const double& scale, const xstructure& str, const _atom& iatom); // atom.fpos are F components per COLUMS !
_atom C2F(const xstructure& str, const _atom& iatom); // atom.cpos are C components per COLUMS !
aurostd::xmatrix<double> FF2CC(const double& scale, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& fmat); // fmat is an operation in F coordinates
aurostd::xmatrix<double> FF2CC(const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& fmat); // fmat is an operation in F coordinates
aurostd::xmatrix<double> CC2FF(const double& scale, const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& cmat); // cmat is an operation in C coordinates
aurostd::xmatrix<double> CC2FF(const aurostd::xmatrix<double>& lattice, const aurostd::xmatrix<double>& cmat); // cmat is an operation in C coordinates
// DX20190905 START
// BringInCellInPlace() overloads
void BringInCellInPlace(double&, double = _ZERO_TOL_, double = 1.0, double = 0.0); // ME+DX20190409
void BringInCellInPlace(aurostd::xvector<double>&, double = _ZERO_TOL_, double = 1.0, double = 0.0); // ME+DX20190409
void BringInCellInPlace(_atom& atom_in, const aurostd::xmatrix<double>& lattice, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904
void BringInCellInPlace(xstructure& xstr, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904

// BringInCell() overloads
double BringInCell(double, double = _ZERO_TOL_, double = 1.0, double = 0.0); // ME+DX20190409
aurostd::xvector<double> BringInCell(const aurostd::xvector<double>& fpos_in, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904
_atom BringInCell(const _atom& atom_in, const aurostd::xmatrix<double>& lattice, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904
xstructure BringInCell(const xstructure& xstr_in, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904

// BringInCellFPOS overloads
void BringInCellInPlaceFPOS(_atom& atom_in, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904
_atom BringInCellFPOS(const _atom& atom_in, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904
// DX20190905 END
// DX+CO START
// DX+CO END
xstructure IdenticalAtoms(const xstructure& a); // Make identical atoms
// xstructure SwapSpecies(const xstructure& a,const uint& A,const uint& B);       // Permute Species A with B (safe for species C).
// xstructure SwapCoordinates(const xstructure& str,const uint& i,const uint& j); // Permute Coordinates i with j
// string SpeciesLabel(const xstructure& a,const uint& A);                        // Returns the Label of the specie A (if available)
// string SpeciesString(const xstructure& a);                                           // Gives a string with the list of all the species
bool GetNiggliCell(const aurostd::xmatrix<double>& in_lat, aurostd::xmatrix<double>& niggli_lat, aurostd::xmatrix<double>& P, aurostd::xmatrix<double>& Q);
// standard lattice reduction and type
std::string GetLatticeType(aurostd::xmatrix<double> lattice);
std::string GetLatticeType(aurostd::xvector<double> data);
xstructure Standard_Primitive_UnitCellForm(const xstructure& a);
xstructure GetStandardPrimitive(const xstructure& a);
aurostd::xmatrix<double> GetStandardPrimitive(aurostd::xmatrix<double> lattice);
aurostd::xvector<double> GetStandardPrimitive(aurostd::xvector<double> data);
xstructure Standard_Conventional_UnitCellForm(const xstructure& a);
xstructure GetStandardConventional(const xstructure& a);
aurostd::xmatrix<double> GetStandardConventional(aurostd::xmatrix<double> lattice);
aurostd::xvector<double> GetStandardConventional(aurostd::xvector<double> data);
// niggli
xstructure GetNiggliStr(const xstructure& in_str);
aurostd::xmatrix<double> GetNiggliStr(const aurostd::xmatrix<double>& lattice);
xstructure NiggliUnitCellForm(const xstructure& a);
aurostd::xmatrix<double> NiggliUnitCellForm(const aurostd::xmatrix<double>& lattice);
// minkowsky
xstructure MinkowskiBasisReduction(const xstructure& a);
aurostd::xmatrix<double> MinkowskiBasisReduction(const aurostd::xmatrix<double>& lattice);
// optimal lattice reduction
xstructure LatticeReduction(const xstructure& a);
aurostd::xmatrix<double> LatticeReduction(const aurostd::xmatrix<double>& lattice);
// CO20170807 START
std::deque<_atom> foldAtomsInCell(const xstructure& a, const aurostd::xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists = true); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - added check_min_dists bool
std::deque<_atom> foldAtomsInCell(const std::deque<_atom>& atoms, const aurostd::xmatrix<double>& lattice_orig, const aurostd::xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists = true); // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 = added check_min_dists bool
xstructure GetPrimitiveVASP(const xstructure& a);
xstructure GetPrimitiveVASP(const xstructure& a, double tol);
// CO20170807 STOP
//  bring cell in,compact, wigner seitz
xstructure BringInCompact(const xstructure& a);
xstructure BringInWignerSeitz(const xstructure& a);
// primitive stuff
xstructure GetPrimitive_20210322(const xstructure& a, double eps = AUROSTD_MAX_DOUBLE); // DX20210323
xstructure GetPrimitive(const xstructure& a);
xstructure GetPrimitive(const xstructure& a, double tol);
xstructure GetPrimitive1(const xstructure& a);
xstructure GetPrimitive2(const xstructure& a);
xstructure GetPrimitive3(const xstructure& a);
bool isTranslationVector(const xstructure& xstr, const aurostd::xvector<double>& vec, double tolerance = 0.5, bool is_frac = false); // DX20210316
bool IsTranslationFVector(const xstructure& a, const aurostd::xvector<double>& ftvec);
bool IsTranslationCVector(const xstructure& a, const aurostd::xvector<double>& ctvec);
// other eggs
xstructure ReScale(const xstructure& a, const double& in_scale);
xstructure SetScale(const xstructure& a, const double& in_scale);
xstructure SetVolume(const xstructure& a, const double& in_volume);
xstructure InflateLattice(const xstructure& a, const double& coefficient);
xstructure InflateVolume(const xstructure& a, const double& coefficient);
double GetVolume(const xstructure& a);
double Volume(const xstructure& a);
// DX20180726 START
_atom BringCloseToOrigin(_atom& atom, aurostd::xmatrix<double>& f2c);
bool uniqueAtomInCell(_atom& atom, std::deque<_atom>& atoms);
bool alreadyInCell(_atom& atom, std::deque<_atom> atoms);
// DX20180726 END
// DX+CO START
bool atomInCell(const _atom& atom, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20191125 //ME+DX20210203 - added bounds
bool inCell(const aurostd::xvector<double>& pos_vec, double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20191125 - added tolerance  // ME20210128 - added bounds
// DX+CO END
xstructure GetSuperCell(const xstructure& a, const aurostd::xmatrix<double>& sc);
xstructure GetSuperCell(const xstructure& a, const aurostd::xvector<double>& sc);
xstructure GetSuperCell(const xstructure& a, const aurostd::xvector<int>& sc);
xstructure GetSuperCell(const xstructure& a, const int& sc11, const int& sc12, const int& sc13, const int& sc21, const int& sc22, const int& sc23, const int& sc31, const int& sc32, const int& sc33);
xstructure GetSuperCell(const xstructure& a, const int& sc1, const int& sc2, const int& sc3);
// CO START
xstructure GetSuperCell(const xstructure& a,
                        const aurostd::xmatrix<double>& sc,
                        std::vector<int>& sc2pcMap,
                        std::vector<int>& pc2scMap,
                        bool get_symmetry,
                        bool get_full_basis,
                        bool force_supercell_matrix = false,
                        bool force_strict_pc2scMap = false); // DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,
                        const aurostd::xvector<double>& sc,
                        std::vector<int>& sc2pcMap,
                        std::vector<int>& pc2scMap,
                        bool get_symmetry,
                        bool get_full_basis,
                        bool force_supercell_matrix = false,
                        bool force_strict_pc2scMap = false); // DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,
                        const aurostd::xvector<int>& sc,
                        std::vector<int>& sc2pcMap,
                        std::vector<int>& pc2scMap,
                        bool get_symmetry,
                        bool get_full_basis,
                        bool force_supercell_matrix = false,
                        bool force_strict_pc2scMap = false); // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,
                        const int& sc11,
                        const int& sc12,
                        const int& sc13,
                        const int& sc21,
                        const int& sc22,
                        const int& sc23,
                        const int& sc31,
                        const int& sc32,
                        const int& sc33,
                        std::vector<int>& sc2pcMap,
                        std::vector<int>& pc2scMap,
                        bool get_symmetry,
                        bool get_full_basis,
                        bool force_supercell_matrix = false,
                        bool force_strict_pc2scMap = false); // DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,
                        const int& sc1,
                        const int& sc2,
                        const int& sc3,
                        std::vector<int>& sc2pcMap,
                        std::vector<int>& pc2scMap,
                        bool get_symmetry,
                        bool get_full_basis,
                        bool force_supercell_matrix = false,
                        bool force_strict_pc2scMap = false); // DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
// CO END
bool CalculateSymmetry(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose, double radius);
bool CalculateSymmetry(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose);
bool CalculateSymmetry(xstructure& str, bool ossverbose, std::ostream& oss, double radius);
bool CalculateSymmetry(xstructure& str, bool ossverbose, double radius);
bool CalculateSymmetry(xstructure& str, double radius);
bool CalculateSymmetry(xstructure& str, bool ossverbose);
bool CalculateSymmetry(xstructure& str);
void CalculateSymmetryPointGroup(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose);
void CalculateSymmetryPointGroup(xstructure& str, bool ossverbose, std::ostream& oss);
void CalculateSymmetryPointGroup(xstructure& str, bool ossverbose);
void CalculateSymmetryPointGroup(xstructure& str);
void CalculateSymmetryPointGroupCrystal(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose);
void CalculateSymmetryPointGroupCrystal(xstructure& str, bool ossverbose, std::ostream& oss);
void CalculateSymmetryPointGroupCrystal(xstructure& str, bool ossverbose);
void CalculateSymmetryPointGroupCrystal(xstructure& str);
void CalculateSymmetryFactorGroup(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose);
void CalculateSymmetryFactorGroup(xstructure& str, bool ossverbose, std::ostream& oss);
void CalculateSymmetryFactorGroup(xstructure& str, bool ossverbose);
void CalculateSymmetryFactorGroup(xstructure& str);
void CalculateSymmetryPointGroupKLattice(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose);
void CalculateSymmetryPointGroupKLattice(xstructure& str, bool ossverbose, std::ostream& oss);
void CalculateSymmetryPointGroupKLattice(xstructure& str, bool ossverbose);
void CalculateSymmetryPointGroupKLattice(xstructure& str);
void CalculateSymmetryPointGroupKCrystal(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose); // ME20200114
void CalculateSymmetryPointGroupKCrystal(xstructure& str, bool ossverbose, std::ostream& oss); // ME20200114
void CalculateSymmetryPointGroupKCrystal(xstructure& str, bool ossverbose); // ME20200114
void CalculateSymmetryPointGroupKCrystal(xstructure& str); // ME20200114
void CalculateSymmetryPointGroupKPatterson(xstructure& str, bool ossverbose, std::ostream& oss, bool fffverbose); // ME20200129
void CalculateSymmetryPointGroupKPatterson(xstructure& str, bool ossverbose, std::ostream& oss); // ME20200129
void CalculateSymmetryPointGroupKPatterson(xstructure& str, bool ossverbose); // ME20200129
void CalculateSymmetryPointGroupKPatterson(xstructure& str); // ME20200129
xstructure Rotate(const xstructure& a, const aurostd::xmatrix<double>& rm);
xstructure GetLTCell(const aurostd::xmatrix<double>& lt, const xstructure& str);
xstructure GetLTFVCell(const aurostd::xvector<double>& nvec, const double phi, const xstructure& str);
xstructure ShiftPos(const xstructure& a, const aurostd::xvector<double>& shift, bool is_frac); // DX20210111
xstructure ShiftCPos(const xstructure& a, const aurostd::xvector<double>& shift);
xstructure ShiftFPos(const xstructure& a, const aurostd::xvector<double>& shift);
double MaxStructureLattice(const xstructure& str);
double MinStructureLattice(const xstructure& str);
double AtomDist(const xstructure& str, const _atom& atom1, const _atom& atom2);
bool SameAtom(const xstructure& str, const _atom& atom1, const _atom& atom2);
bool SameAtom(const _atom& atom1, const _atom& atom2);
bool DifferentAtom(const xstructure& str, const _atom& atom1, const _atom& atom2);
aurostd::xmatrix<double> GetDistMatrix(const xstructure& a); // CO20171025
std::vector<double> GetNBONDXX(const xstructure& a);
int GenerateGridAtoms(xstructure& str, int i1, int i2, int j1, int j2, int k1, int k2); // DX20191218 [ORIG]
int GenerateGridAtoms(xstructure& str, double radius); // CO20200912 - double
int GenerateGridAtoms(xstructure& str, int d);
int GenerateGridAtoms(xstructure& str, int d1, int d2, int d3);
int GenerateGridAtoms(xstructure& str, const aurostd::xvector<int>& dims);
int GenerateGridAtoms(xstructure& str);

void l2ijk(const xstructure& str, const int& l, int& i, int& j, int& k);
void l2ijk(const xstructure& str, const int& l, aurostd::xvector<int>& ijk);
aurostd::xvector<int> l2ijk(const xstructure& str, const int& l);
void ijk2l(const xstructure& str, int& l, const int& i, const int& j, const int& k);
void ijk2l(const xstructure& str, int& l, const aurostd::xvector<int>& ijk);
int ijk2l(const xstructure& str, const int& i, const int& j, const int& k);
int ijk2l(const xstructure& str, const aurostd::xvector<int>& ijk);
aurostd::xvector<double> r_lattice(const xstructure& str, const int& l);
aurostd::xvector<double> r_lattice(const xstructure& str, const int& i, const int& j, const int& k);
aurostd::xvector<double> r_lattice(const xstructure& str, const aurostd::xvector<int>& ijk);
xstructure input2AIMSxstr(std::istream& input);
xstructure input2ABINITxstr(std::istream& input);
xstructure input2QExstr(std::istream& input);
xstructure input2VASPxstr(std::istream& input, bool vasp5 = false);
xstructure input2ITCxstr(std::istream& input); // CO20220613
xstructure input2ELKxstr(std::istream& input); // DX20200313
xstructure input2ATATxstr(std::istream& input); // SD20220123
xstructure input2LMPxstr(std::istream& input); // SD20240111

// ----------------------------------------------------------------------------
// centroid functions for structures //DX20200728
aurostd::xvector<double> getCentroidOfStructure(const xstructure& xstr, bool use_cpos = true, bool use_atom_mass = false);
aurostd::xvector<double> getCentroidOfStructure(const std::deque<_atom>& atoms, bool use_cpos = true, bool use_atom_mass = false);
aurostd::xvector<double> getCentroidOfStructurePBC(const xstructure& xstr, bool use_cpos = true, bool use_atom_mass = false);
aurostd::xvector<double> getCentroidOfStructurePBC(const std::deque<_atom>& atoms, aurostd::xmatrix<double> lattice, bool use_cpos = true, bool use_atom_mass = false);

// ----------------------------------------------------------------------------
// functions related to AtomEnvironment - DX20191122
std::vector<AtomEnvironment> getAtomEnvironments(const xstructure& xstr, uint mode = ATOM_ENVIRONMENT_MODE_1);
void writeAtomEnvironments(const std::vector<AtomEnvironment> environments, const std::map<std::string, std::string> meta_data = std::map<std::string, std::string>()); // HE20210723 add separate write function
std::vector<AtomEnvironment> getLFAAtomEnvironments(const xstructure& xstr, const std::string& lfa, const std::vector<std::string>& LFAs, uint mode = ATOM_ENVIRONMENT_MODE_1);

void minimumCoordinationShellLatticeOnly(const aurostd::xmatrix<double>& lattice, double& min_dist, uint& frequency, std::vector<aurostd::xvector<double>>& coordinates); // DX20191122
void minimumCoordinationShellLatticeOnly(const aurostd::xmatrix<double>& lattice, double& min_dist, uint& frequency, std::vector<aurostd::xvector<double>>& coordinates, double radius); // DX20191122
void minimumCoordinationShellLatticeOnly(const aurostd::xmatrix<double>& lattice,
                                         aurostd::xvector<int>& dims,
                                         std::vector<aurostd::xvector<double>>& l1,
                                         std::vector<aurostd::xvector<double>>& l2,
                                         std::vector<aurostd::xvector<double>>& l3,
                                         std::vector<int>& a_index,
                                         std::vector<int>& b_index,
                                         std::vector<int>& c_index,
                                         double& min_dist,
                                         uint& frequency,
                                         std::vector<aurostd::xvector<double>>& coordinates,
                                         double radius); // DX20191122
void minimumCoordinationShell(const xstructure& xstr, uint center_index, double& min_dist, uint& frequency, std::vector<aurostd::xvector<double>>& coordinates); // DX20191122
void minimumCoordinationShell(const xstructure& xstr, uint center_index, double& min_dist, uint& frequency, std::vector<aurostd::xvector<double>>& coordinates, const std::string& type); // DX20191122

// makefile tests
bool smithTest(std::ostream& oss = std::cout);
bool smithTest(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
bool PrototypeGeneratorTest(std::ostream& oss = std::cout, bool check_symmetry = false, bool check_uniqueness = false); // DX20200928
bool PrototypeGeneratorTest(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool check_symmetry = false, bool check_uniqueness = false); // DX20200928
// ----------------------------------------------------------------------------
// Structure Prototypes

// aflow_xproto.cpp
namespace aflowlib {
  std::string PrototypeCleanLatticeString(const std::string& latticeIN);
}
double NearestNeighbor(const xstructure& str_in);
std::vector<double> NearestNeighbors(const xstructure& xstr); // DX20201230 - moved from XtalFinder
double NearestNeighborToAtom(const xstructure& xstr, uint k); // DX20201230 - moved from XtalFinder

namespace aflowlib {
  struct _PROTO_PARAMS {
    std::string label;
    std::string parameters;
    std::deque<std::string> vatomX;
    std::deque<double> vvolumeX;
    double volume_in;
    int mode;
    bool flip_option;
  };

  uint PrototypeLibrariesSpeciesNumber(const std::string& label); // CO20181226
  std::set<std::string> GetPrototypesByStoichiometry(const std::vector<uint>& stoichiometry, const std::string& library = "all"); // DX20181009
  std::set<std::string> GetPrototypesBySymmetry(const std::vector<uint>& stoichiometry,
                                                uint space_group_number,
                                                const std::vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions,
                                                uint setting,
                                                const std::string& library = "all"); // DX20181010
  xstructure PrototypeLibraries(std::ostream& oss, const std::string& label, const std::string& parameters, int mode);
  xstructure PrototypeLibraries(std::ostream& oss, const std::string& label, const std::string& parameters, std::deque<std::string>& vatomX, int mode);
  xstructure PrototypeLibraries(std::ostream& oss, const std::string& label, const std::string& parameters, std::deque<std::string>& vatomX, std::deque<double>& vvolumeX, double volume_in, int mode); //=LIBRARY_MODE_HTQC);
  xstructure PrototypeLibraries(std::ostream& oss, std::string label, std::string parameters, std::deque<std::string>& vatomX, std::deque<double>& vvolumeX, double volume_in, int mode, bool flip_option);
  xstructure PrototypeLibraries(std::ostream& oss, _PROTO_PARAMS* PARAMS);

  std::string PrototypesHelp();
  std::string PrototypesIcsdHelp(const std::string& options);
  std::string CALCULATED();
  std::string CALCULATED_ICSD_RANDOM();
} // namespace aflowlib

extern std::string PrototypeBinaryGUS_Cache_Library[];

// ----------------------------------------------------------------------------
// Various prototypes to be moved somewhere sometime
// PROTOTYPES
// uint argsprint(std::vector<std::string> argv);
// ----------------------------------------------------------------------------
// aflow.cpp
//[CO20200502 - DUPLICATE?]string aflow_get_time_string(void);
//[CO20200502 - DUPLICATE?]string aflow_get_time_string_short(void);
//[CO20200502 - DUPLICATE?]string strPID(void);
//[CO20200502 - DUPLICATE?]string strTID(void);  //CO20200502 - threadID
int AFLOW_main(std::vector<std::string>& argv);
namespace aflow {
  std::string License_Preamble_aflow();
  std::string Intro_aflow(std::string x);
  std::string Intro_sflow(std::string x);
  std::string Intro_HELP(std::string x);
  std::string Banner(std::string type);
} // namespace aflow
int VASP_Main(std::vector<std::string> argv);
int GRND_Main(std::vector<std::string> argv);
namespace KBIN {
  int KBIN_Main(std::vector<std::string> argv);
}
std::string MessageTime();
std::string MessageHostTime(const _aflags& aflags);
std::string MessageDir(const _aflags& aflags);
std::string MessageDirTime(const _aflags& aflags);
std::string MessageDirHostTime(const _aflags& aflags);
//[CO20200624 - REDUNDANT]bool AFLOW_BlackList(string hostname);

// interfaces
namespace KBIN {
  bool MoveRun2NewDirectory(_aflags& aflags, const std::string& subdirectory_orig, const std::string& subdirectory_new); // DX20210901 //SD20220319 - return bool
  void RUN_Directory_PTHREADS(_aflags& aflags);
  void* _threaded_interface_RUN_Directory(void* ptr);
} // namespace KBIN

// ----------------------------------------------------------------------------
// aflow_kbin.cpp
// int KbinCheckInputFiles(string Directory,ofstream& FileERROR);
namespace KBIN {
  void MPI_Extract(std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  void RUN_Directory(_aflags& aflags);
  void AFLOW_RUN_Directory(const _aflags& aflags);
  void RUN_DirectoryScript(const _aflags& aflags, const std::string& script, const std::string& output);
  bool CompressDirectory(const std::string& directory);
  bool CompressDirectory(const _aflags& aflags);
  bool DecompressDirectory(const std::string& directory);
  bool DecompressDirectory(const _aflags& aflags);
  void Clean(const _aflags& aflags);
  void Clean(const std::string& directory);
  void Clean(const _aflags& aflags, const aurostd::xoption& opts_clean);
  void Clean(const std::string& directory, const aurostd::xoption& opts_clean);
  void XClean(std::string options);
  void GenerateAflowinFromVASPDirectory(_aflags& aflags);
  void StartStopCheck(const std::string& AflowIn, std::string str1, std::string str2, bool& flag, bool& flagS);
  void StartStopCheck(const std::string& AflowIn, std::string str1, bool& flag, bool& flagS);
  bool Legitimate_krun(const _aflags& aflags, const bool osswrite, std::ostringstream& oss); // SD20220224
  bool Legitimate_krun(const _aflags& aflags); // SD20220224
  bool Legitimate_aflowin(const std::string& aflowindir, const bool osswrite, std::ostringstream& oss); // SD20220224 - made aflowindir const, removed reference from bool
  bool Legitimate_aflowin(const std::string& aflowindir); // SD20220224 - made aflowindir const
  bool Legitimate_aflowdir(const std::string& aflowindir, const _aflags& aflags, const bool osswrite, std::ostringstream& oss); // SD20220224
  bool Legitimate_aflowdir(const std::string& aflowindir, const _aflags& aflags); // SD20220224
  void getAflowInFromAFlags(const _aflags& aflags, std::string& AflowIn_file, std::string& AflowIn, std::ostream& oss = std::cout); // CO20191110
  void getAflowInFromAFlags(const _aflags& aflags, std::string& AflowIn_file, std::string& AflowIn, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20191110
  void getAflowInFromDirectory(const std::string& directory, std::string& AflowIn_file, std::string& AflowIn, std::ostream& oss = std::cout); // CO20191110
  void getAflowInFromDirectory(const std::string& directory, std::string& AflowIn_file, std::string& AflowIn, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20191110
  int get_NCPUS(); // ME20200219
  int get_NCPUS(const _kflags&); // ME20200219

  // ----------------------------------------------------------------------------
  // aflow_modules.cpp
  // ME20181027
  void setModules(_xvasp&);
  void setModules(_xinput&);
  void readModulesFromAflowIn(const std::string&, _kflags&, _xvasp&);
  void readModulesFromAflowIn(const std::string&, _kflags&, _xinput&);
  std::vector<aurostd::xoption> loadDefaultsAPL();
  bool writeFlagAPL(const std::string& key, const aurostd::xoption& xopt); // CO20181226  //ME20190113
  void readParametersAPL(const std::string&, _moduleOptions&, _xinput&);
  std::vector<aurostd::xoption> loadDefaultsQHA(); // AS20200302
  void readParametersQHA(const std::string&, _moduleOptions&, _xinput&); // AS20200302
  std::vector<aurostd::xoption> loadDefaultsAAPL();
  bool writeFlagAAPL(const std::string& key, const aurostd::xoption& xopt); // CO20181226  //ME20190113
  void readParametersAAPL(const std::string&, _moduleOptions&, _xinput&);
  std::vector<aurostd::xoption> loadDefaultsAEL();
  bool writeFlagAEL(const std::string& key, const aurostd::xoption& xopt);
  std::vector<aurostd::xoption> loadDefaultsAGL();
  bool writeFlagAGL(const std::string& key, const aurostd::xoption& xopt);

  // ----------------------------------------------------------------------------
  // aflow_qsub.cpp
  bool QSUB_Extract(_xqsub& xqsub, std::string AflowIn, std::ifstream& FileAFLOWIN, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  bool QSUB_RunFinished(_aflags& aflags, std::ofstream& FileMESSAGE, bool = false);
  void QSUB_WaitFinished(_aflags& aflags, std::ofstream& FileMESSAGE, bool = false);
  bool QSUB_Extract_Mode1(_xqsub& xqsub, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  bool QSUB_Extract_Mode2(_xqsub& xqsub, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  bool QSUB_Extract_Mode3(_xqsub& xqsub, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
} // namespace KBIN

// ----------------------------------------------------------------------------
// aflow_ialien.cpp
namespace ALIEN {
  bool Produce_INPUT(_xalien& xalien, std::string AflowIn, std::ifstream& FileAFLOWIN, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _alienflags& alienflags);
  bool Modify_INPUT(_xalien& xalien, std::ofstream& FileMESSAGE, _aflags& aflags, _alienflags& alienflags);
  bool Write_INPUT(_xalien& xalien);
  bool Produce_INPUT_FILE(_xalien& xalien, std::string AflowIn, std::ifstream& FileAFLOWIN, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _alienflags& alienflags);
  bool Modify_INPUT_FILE(_xalien& xalien, std::ofstream& FileMESSAGE, _aflags& aflags, _alienflags& alienflags);

  // ----------------------------------------------------------------------------
  // aflow_kalien.cpp
  _alienflags Get_Alienflags_from_AflowIN(std::string& AflowIn);
  bool Run_Directory(std::ofstream& FileERROR, _aflags& aflags, _kflags& kflags);
} // namespace ALIEN

// ----------------------------------------------------------------------------
// aflow_matlab.cpp aflow_matlab_funcs.cpp
bool KBIN_MATLAB_Extract(std::string AflowIn, std::ifstream& FileAFLOWIN, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
bool KBIN_MATLAB_Run(_kflags& kflags);
_kflags KBIN_MATLAB_Get_Matlabflags_from_AflowIN(std::string& AflowIn);
bool KBIN_MATLAB_Directory(std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);

// ----------------------------------------------------------------------------
// aflow_ifrozsl.cpp

namespace KBIN {
  void VASP_RunPhonons_FROZSL(_xvasp& xvasp, std::string AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMESSAGE);
}

namespace FROZSL {
  bool Extract_INPUT(const std::string& AflowIn, std::ofstream& FileMESSAGE, std::stringstream& input_file, _aflags& aflags, _kflags& kflags);
  bool Setup_frozsl_init_input(const std::string& AflowIn, std::ofstream& FileMESSAGE, std::stringstream& input_file, _aflags& aflags, _kflags& kflags);
  bool Already_Calculated_Input(const std::string& AflowIn);
  bool WGET_INPUT(std::ofstream& FileMESSAGE, std::string AflowIn, _aflags& aflags, _kflags& kflags);
  bool WGET_OUTPUT(std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  bool input_TO_poscar(std::ofstream& FileMESSAGE, std::stringstream& input_file, _aflags& aflags, _kflags& kflags);
  std::string Generate_Input_file(std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  bool File_INPUT(const std::string& AflowIn, std::ofstream& FileMESSAGE, std::stringstream& input_file, _aflags& aflags, _kflags& kflags);
  bool Write(const std::string& filename, const std::string& directory);
  bool Delete(const std::string& filename, const std::string& directory);
} // namespace FROZSL
namespace FINDSYM {
  bool Write(const std::string& filename, const std::string& directory);
}

// -------------------------------------------------------------------------------------------------
// aflow_ovasp.cpp
// TODO should move to a new ovasp.h
class xOUTCAR;
class xDOSCAR;
class xEIGENVAL;
class xPOTCAR;
class xVASPRUNXML;
class xIBZKPT;
class xKPOINTS;
class xCHGCAR;
class xVASPOUT;
class xQMVASP;
class xPLASMONICS; // CO20190803
namespace aflowlib {
  class _aflowlib_entry;
}

namespace aurostd {
  // REGEX expressions for quick finding/replacements in strings
  /// REGEX to find all chemical elements in a string
  const std::regex regex_elements{
      "(A[cglmrstu]|B[aehikr]?|C[adeflmnorsu]?|D[bsy]|E[rsu]|F[elmr]?|G[ade]|H[efgos]?|I[nr]?|Kr?|L[airuv]|M[dgnot]|N[abdeiop]?|Os?|P[abdmortu]?|R[abefghnu]|S[bcegimnr]?|T[abcehilm]|U(u[opst])?|V|W|Xe|Yb?|Z["
      "nr])"};
  /// REGEX to find all pseudo potentials that contain uppercase letters from a string (could be mistaken for a chemical element)
  const std::regex regex_ppclean{"(" + std::regex_replace(CAPITAL_LETTERS_PP_LIST, std::regex(","), "|") + ")"};
  /// @brief REGEX to help change a AURL into a file path
  /// @note the content of the group `((?:(?:LIB\d{1,})|(?:ICSD)))` can be used in the replacement  with `$1`;
  ///       the second group `(?:(?:RAW)|(?:LIB)|(?:WEB))` is there to select the full substring to be replaced
  const std::regex regex_aurl2file{"((?:(?:LIB\\d{1,})|(?:ICSD)))_(?:(?:RAW)|(?:LIB)|(?:WEB))\\/"};
} // namespace aurostd

// -------------------------------------------------------------------------------------------------
/// @brief Abstract base class for the vasp data containers such as @c xOUTCAR, @c xIBZKPT, etc.
/// Implements some convenience static methods to construct derived classes from files. Derived classes must
/// implement @c GetProperties to read from string, @c GetPropertiesFile to read from file, and @c clear to reset the object to default.
/// Derived class may use @code using xOVASP::xOVASP; @endcode to inherit constructors.
template <class Derived> class xOVASP : public xStream {
public:
  xOVASP(std::ostream& oss = std::cout) : xStream(oss) {}
  xOVASP(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss) {}

  static Derived fromFile(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) {
    Derived derived(oss);
    derived.m_initialized = derived.initialize(fileIN, QUIET);
    return derived;
  }
  static Derived fromFile(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) {
    Derived derived(FileMESSAGE, oss);
    derived.m_initialized = derived.initialize(fileIN, QUIET);
    return derived;
  }

  bool initialize(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) {
    if (aurostd::CompressFileExist(fileIN, filename)) {
      return GetPropertiesFile(filename, QUIET);
    }
    return false;
  }

  virtual bool GetProperties(const std::string& stringIN, bool QUIET = true) = 0;
  virtual bool GetPropertiesFile(const std::string& fileIN, bool QUIET = true) = 0;

  virtual void clear() = 0;

  bool m_initialized = false;

  // CONTENT
  std::string content;
  std::vector<std::string> vcontent;
  std::string filename;
};

class xOUTCAR : public xOVASP<xOUTCAR>, public JsonSerializable<xOUTCAR> {
public:
  using xOVASP::xOVASP;

  xOUTCAR(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xOUTCAR(fromFile(fileIN, QUIET, oss)) {}
  xOUTCAR(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xOUTCAR(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xOUTCAR(); }

  std::string SYSTEM;
  int NELM = 0;
  int NIONS = 0;
  double Efermi = 0.0;
  bool isLSCOUPLING = false;
  aurostd::xvector<double> efield_pead;
  int nelectrons = 0; // AS20200528
  std::vector<aurostd::xvector<double>> vmag_noncoll; // DX20171205 - non-collinear

  // for aflowlib_libraries.cpp
  std::vector<double> vmag;
  std::vector<aurostd::xvector<double>> vforces;
  std::vector<aurostd::xvector<double>> vpositions_cartesian;
  aurostd::xmatrix<double> stress;
  double natoms = 0.0;
  double energy_cell = 0.0;
  double energy_atom = 0.0;
  double enthalpy_cell = 0.0;
  double enthalpy_atom = 0.0;
  double eentropy_cell = 0.0;
  double eentropy_atom = 0.0;
  double PV_cell = 0.0;
  double PV_atom = 0.0;
  double mag_cell = 0.0;
  double mag_atom = 0.0;
  double volume_cell = 0.0;
  double volume_atom = 0.0;
  double pressure = 0.0; // SAME AS PSTRESS
  double pressure_residual = 0.0;
  double Pulay_stress = 0.0;
  double ENCUT = 0.0;
  double EDIFF = 0.0;
  double EDIFFG = 0.0;
  double POTIM = 0.0;
  double TEIN = 0.0;
  double TEBEG = 0.0;
  double TEEND = 0.0;
  double SMASS = 0.0;
  double NPACO = 0.0;
  double APACO = 0.0;
  double PSTRESS = 0.0;
  int NBANDS = 0;
  int NKPTS = 0;
  int NSW = 0;
  int NBLOCK = 0;
  int KBLOCK = 0;
  int IBRION = 0;
  int NFREE = 0;
  int ISIF = 0;
  int IWAVPR = 0;
  int ISYM = 0;
  int ISPIN = 0;
  double total_energy_change = 0.0;
  // DOS related values
  double EMIN = 0.0;
  double EMAX = 0.0;
  double SIGMA = 0.0; // eV - energy-range for DOS
  int ISMEAR = 0; // broadening in eV -4-tet -1-fermi 0-gaus
  //  Electronic relaxation
  int IALGO = 0; // algorithm
  std::string LDIAG; // sub-space diagonalisation
  int IMIX = 0;
  int INIMIX = 0;
  int MIXPRE = 0; // mixing-type and parameters
  double AMIX = 0.0;
  double BMIX = 0.0;
  double AMIX_MAG = 0.0;
  double BMIX_MAG = 0.0;
  double AMIN = 0.0;
  double WC = 0.0; // parameters
  // Intra band minimization
  double WEIMIN = 0.0;
  double EBREAK = 0.0;
  double DEPER = 0.0;
  double TIME = 0.0;
  // begin shared xPOTCAR
  double ENMAX = 0.0;
  std::vector<double> vENMAX; // eV
  double ENMIN = 0.0;
  std::vector<double> vENMIN; // eV
  double POMASS_sum = 0.0;
  double POMASS_min = 0.0;
  double POMASS_max = 0.0;
  std::vector<double> vPOMASS; // mass
  double ZVAL_sum = 0.0;
  double ZVAL_min = 0.0;
  double ZVAL_max = 0.0;
  std::vector<double> vZVAL; // valence
  double EATOM_min = 0.0;
  double EATOM_max = 0.0;
  std::vector<double> vEATOM; // eV
  double RCORE_min = 0.0;
  double RCORE_max = 0.0;
  std::vector<double> vRCORE; // outmost cutoff radius
  double RWIGS_min = 0.0;
  double RWIGS_max = 0.0;
  std::vector<double> vRWIGS; // wigner-seitz radius (au A)
  double EAUG_min = 0.0;
  double EAUG_max = 0.0;
  std::vector<double> vEAUG; // augmentation
  double RAUG_min = 0.0;
  double RAUG_max = 0.0;
  std::vector<double> vRAUG; // augmentation
  double RMAX_min = 0.0;
  double RMAX_max = 0.0;
  std::vector<double> vRMAX; // unicity
  std::vector<std::string> vTITEL; // unicity
  std::vector<std::string> vLEXCH; // unicity
  // end shared xPOTCAR
  std::string pp_type;
  std::vector<std::string> species; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<int> species_Z; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_type; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_version; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_AUID; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_AUID_collisions; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<double> species_pp_groundstate_energy; // meV/atom
  std::vector<std::string> species_pp_groundstate_structure; // name that we have, maybe ANRL
  std::deque<std::deque<double>> species_pp_vLDAU; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  bool isKIN = false; // METAGGA
  bool isMETAGGA = false;
  std::string METAGGA; // METAGGA
  std::string string_LDAU;
  uint nweights = 0;
  uint nkpoints_irreducible = 0; // kpoints reading
  std::vector<aurostd::xvector<double>> vkpoint_reciprocal; // kpoints reading
  std::vector<aurostd::xvector<double>> vkpoint_cartesian; // kpoints reading
  std::vector<double> vweights; // kpoints reading
  double calculation_time = 0.0;
  double calculation_memory = 0.0;
  uint calculation_cores = 1;
  xstructure xstr; // for GetBandGap()
  std::vector<xstructure> vxstr_ionic; // for all ionic steps  //CO20211106
  std::vector<double> venergy_ionic; // for all ionic steps  //CO20211106
  std::vector<aurostd::xvector<double>> vstresses_ionic; // for all ionic steps  //CO20211106
  std::vector<std::string> GetCorrectEntriesFromLine(const std::string& line, const uint expected_count); // CO20170725 - vasp issues with lattice spacing (negative sign)
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, uint natoms_check, bool = true); // get everything QUIET  //CO20200404 - added default for bool
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
  std::vector<int> band_index;
  std::vector<int> carrier_spin;
  std::vector<std::string> carrier_type;
  std::vector<std::vector<double>> extrema_cart_coord;
  std::vector<std::vector<double>> effective_mass_axes;
  std::vector<int> equivalent_valley;
  std::vector<double> effective_mass_DOS;
  std::vector<double> effective_mass_COND;
  std::vector<double> mass_elec_dos;
  std::vector<double> mass_hole_dos;
  std::vector<double> mass_elec_conduction;
  std::vector<double> mass_hole_conduction;
  // BAND GAPS
  bool GetXStructure();
  int isKPointLine(uint iline, aurostd::xvector<double>& kpoint); // if returns 0 if not KPointLine, -1 means it gave *** for kpoint
  int isKPointLine(uint iline); // if returns 0 if not KPointLine, -1 means it gave *** for kpoint
  bool GetStartingKPointLines(std::vector<uint>& ilines);
  bool GetNextKPointLine(uint& iline);
  bool ProcessKPoint(uint iline, double EFERMI, std::vector<double>& b_energies, std::vector<double>& b_occs);
  bool GetBandEdge(std::vector<double>& b_energies, std::vector<double>& b_occs, double EFERMI, uint& iedge, double efermi_tol = AUROSTD_NAN, double energy_tol = 1e-4, double occ_tol = 1e-5);
  bool identicalKPoints(std::vector<aurostd::xvector<double>>& vkpoints, uint kpt1, uint kpt2, double tol = 1e-12);
  bool identicalKPoints(aurostd::xvector<double>& kpoint1, aurostd::xvector<double>& kpoint2, double tol = 1e-12);
  bool removeDuplicateKPoints(std::vector<aurostd::xvector<double>>& vkpoints, std::vector<uint>& vikpt);
  bool removeDuplicateKPoints(std::vector<std::vector<aurostd::xvector<double>>>& vkpoints, std::vector<uint>& vikpt, std::vector<uint>& vispin);
  double minimumDistanceKPoints(std::vector<aurostd::xvector<double>>& vkpoints, uint ikp1, uint ikp2);
  double minimumDistanceKPoints(aurostd::xvector<double>& kpoint1, aurostd::xvector<double>& kpoint2);
  struct bandEnergyOcc {
    double energy;
    double occ;
  };
  struct bandEnergyOccCompare {
    bandEnergyOccCompare(double _energy_tol) : energy_tol(_energy_tol) {};
    double energy_tol;
    bool operator()(const bandEnergyOcc& a, const bandEnergyOcc b) const;
  };
  bool orderBands(std::vector<double>& b_energies, std::vector<double>& b_occs, double energy_tol = 1e-4);
  enum BROAD_TYPES { empty, metal, insulator }; // bandgap types
  enum EMPTY_TYPES { empty_all, empty_partial }; // bandgap types
  enum INSULATOR_TYPES { insulator_direct, insulator_indirect }; // bandgap types
  enum GAP_TYPES { zero_gap, non_zero_gap }; // bandgap types
  bool GetBandGap(double EFERMI = AUROSTD_NAN, double efermi_tol = AUROSTD_NAN, double energy_tol = 1e-4, double occ_tol = 1e-5);
  std::vector<double> conduction_band_min;
  double conduction_band_min_net = 0.0;
  std::vector<double> valence_band_max;
  double valence_band_max_net = 0.0;
  std::vector<double> Egap;
  double Egap_net = 0.0;
  std::vector<double> Egap_fit;
  double Egap_fit_net = 0.0;
  std::vector<std::string> Egap_type;
  std::string Egap_type_net;
  // DIELECTRIC
  aurostd::xmatrix<double> freq_plasma; //extracted from OUTCAR
  std::vector<double> freq_grid;
  aurostd::xmatrix<double> dielectric_static;
  std::vector<aurostd::xmatrix<double>> dielectric_interband_real; //extracted from OUTCAR
  std::vector<aurostd::xmatrix<double>> dielectric_interband_imag; //extracted from OUTCAR
  double freq_plasma_iso = 0.0;
  std::vector<double> dielectric_interband_iso_real;
  std::vector<double> dielectric_interband_iso_imag;
  std::vector<aurostd::xmatrix<double>> dielectric_drude_real;
  std::vector<aurostd::xmatrix<double>> dielectric_drude_imag;
  std::vector<double> dielectric_drude_iso_real;
  std::vector<double> dielectric_drude_iso_imag;
  std::vector<double> energy_loss_function_iso;
  std::vector<double> reflectivity_iso;
  std::vector<double> dielectric_full_iso_real;
  std::vector<double> dielectric_full_iso_imag;
  std::vector<aurostd::xmatrix<double>> dielectric_full_real;
  std::vector<aurostd::xmatrix<double>> dielectric_full_imag;
  void GetDielectricData();
  bool GetOptical(double freq_relax = 0.2);
  // EXPORT
  bool GetIonicStepsData();
  void AddStepsIAPCFG(aurostd::JSON::object& jo, aflowlib::_aflowlib_entry& entry);
  friend std::ostream& operator<<(std::ostream&, const xOUTCAR&); // ME20190623

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xOUTCAR deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xOUTCAR"; }

private: //
  void copy(const xOUTCAR& b); //

  // SERIALIZATION MEMBERS
#define JSON_xOUTCAR_MEMBERS                                                                                                                                                                                     \
  m_initialized, content, vcontent, filename, SYSTEM, NELM, NIONS, Efermi, isLSCOUPLING, efield_pead, nelectrons, natoms, energy_cell, energy_atom, enthalpy_cell, enthalpy_atom, eentropy_cell, eentropy_atom,  \
      PV_cell, PV_atom, stress, mag_cell, mag_atom, vmag, vmag_noncoll, volume_cell, volume_atom, pressure, pressure_residual, Pulay_stress, vforces, vpositions_cartesian, ENCUT, EDIFF, EDIFFG, POTIM, TEIN,   \
      TEBEG, TEEND, SMASS, NPACO, APACO, PSTRESS, NBANDS, NKPTS, NSW, NBLOCK, KBLOCK, IBRION, NFREE, ISIF, IWAVPR, ISYM, ISPIN, total_energy_change, EMIN, EMAX, SIGMA, ISMEAR, IALGO, LDIAG, IMIX, INIMIX,      \
      MIXPRE, AMIX, BMIX, AMIX_MAG, BMIX_MAG, AMIN, WC, WEIMIN, EBREAK, DEPER, TIME, ENMAX, vENMAX, ENMIN, vENMIN, POMASS_sum, POMASS_min, POMASS_max, vPOMASS, ZVAL_sum, ZVAL_min, ZVAL_max, vZVAL, EATOM_min,  \
      EATOM_max, vEATOM, RCORE_min, RCORE_max, vRCORE, RWIGS_min, RWIGS_max, vRWIGS, EAUG_min, EAUG_max, vEAUG, RAUG_min, RAUG_max, vRAUG, RMAX_min, RMAX_max, vRMAX, vTITEL, vLEXCH, pp_type, species,          \
      species_Z, species_pp, species_pp_type, species_pp_version, species_pp_AUID, species_pp_AUID_collisions, species_pp_groundstate_energy, species_pp_groundstate_structure, species_pp_vLDAU, isKIN,         \
      isMETAGGA, METAGGA, string_LDAU, nweights, nkpoints_irreducible, vkpoint_reciprocal, vkpoint_cartesian, vweights, calculation_time, calculation_memory, calculation_cores, xstr, vxstr_ionic,              \
      venergy_ionic, vstresses_ionic, band_index, carrier_spin, carrier_type, extrema_cart_coord, effective_mass_axes, equivalent_valley, effective_mass_DOS, effective_mass_COND, mass_elec_dos, mass_hole_dos, \
      mass_elec_conduction, mass_hole_conduction
};

// EFFECTIVE MASSES //CO20200404 - moved from "friend" of xOUTCAR
bool GetEffectiveMass(xOUTCAR& outcar, xDOSCAR& doscar, xEIGENVAL& eigenval, xstructure xstr, std::ostream& oss = std::cout); // CO20200404
bool GetEffectiveMass(xOUTCAR& outcar, xDOSCAR& doscar, xEIGENVAL& eigenval, xstructure xstr, std::ofstream& FileMeSSAGE, std::ostream& oss = std::cout); // CO20200404

//-------------------------------------------------------------------------------------------------
class xDOSCAR : public xOVASP<xDOSCAR>, public JsonSerializable<xDOSCAR> {
public:
  using xOVASP::xOVASP;

  xDOSCAR(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xDOSCAR(fromFile(fileIN, QUIET, oss)) {}
  xDOSCAR(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xDOSCAR(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xDOSCAR(); }

  std::string title;
  uint spin = 0;
  double Vol = 0.0;
  double POTIM = 0.0;
  aurostd::xvector<double> lattice; // CO20200922 - an xvector in the style of Getabc_angles(), only the abc are printed/read, must be in meters: https://www.vasp.at/wiki/index.php/DOSCAR
  double temperature = 0.0;
  bool RWIGS = false;
  double Efermi = 0.0;
  double spinF = AUROSTD_NAN;
  double energy_max = 0.0;
  double energy_min = 0.0;
  uint number_energies = 0;
  uint number_atoms = 0; // ME20190614
  bool partial = false; // ME20190614
  double denergy = 0.0;
  std::deque<double> venergy; // venergy.at(energy_number)
  std::deque<double> venergyEf; // venergyEf.at(energy_number)
  // ME20190614 BEGIN
  std::deque<std::deque<double>> viDOS; // viDOS.at(spin).at(energy_number)
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS; // vDOS.at(atom).at(orbital).at(spin).at(energy_number); 0 = total for atoms and orbitals
  // ME20190614 END
  // ME20190620 BEGIN
  bool isLSCOUPLING = false; // Contains spin-orbit coupling
  bool lmResolved = false; // Is it lm-resolved?
  std::string carstring; // The fourth line of the DOSCAR
  // ME20190620 END
  std::vector<double> conduction_band_min; // CO20191004
  double conduction_band_min_net = AUROSTD_NAN; // CO20191004
  std::vector<double> valence_band_max; // CO20191004
  double valence_band_max_net = AUROSTD_NAN; // CO20191004
  std::vector<double> Egap; // CO20191004
  double Egap_net = AUROSTD_NAN; // CO20191004
  std::vector<double> Egap_fit; // CO20191004
  double Egap_fit_net = AUROSTD_NAN; // CO20191004
  std::vector<std::string> Egap_type; // CO20191004
  std::string Egap_type_net; // CO20191004
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
  void convertSpinOFF2ON(); // CO20191217 - copies everything from spin channel 1 to spin channel 2
  void addAtomChannel(); // CO20211124 - creates another atom channel, mimicking size of orbital, spin, and energy
  void addOrbitalChannel(); // CO20211124 - creates another orbital channel, mimicking sizes of spin and energy
  void resetVDOS(); // CO20211124 - set all vDOS to 0
  bool checkDOS(std::string& ERROR_out) const; // CO20191010
  bool GetBandGap(double EFERMI = AUROSTD_NAN, double efermi_tol = AUROSTD_NAN, double energy_tol = 1e-3, double occ_tol = 1e-4); // CO20191110
  [[nodiscard]] std::deque<std::deque<std::deque<std::deque<double>>>> GetVDOSSpecies(const xstructure& xstr) const; // vDOS.at(species).at(spin).at(energy_number)  //CO20191110
  [[nodiscard]] std::deque<std::deque<std::deque<std::deque<double>>>> GetVDOSSpecies(std::deque<int> num_each_type) const; // vDOS.at(species).at(spin).at(energy_number)  //CO20191110
  friend std::ostream& operator<<(std::ostream&, const xDOSCAR&); // ME20190623

  //SD20230214 BEGIN
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS_atom; // vDOS.at(atom).at(orbital).at(spin).at(energy_number); 0 = total for atoms and orbitals
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS_lm_atom; // vDOS.at(atom).at(orbital_lm).at(spin).at(energy_number); 0 = total for atoms and orbitals
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS_species; // vDOS.at(species).at(orbital).at(spin).at(energy_number); 0 = total for atoms and orbitals
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS_lm_species; // vDOS.at(species).at(orbital_lm).at(spin).at(energy_number); 0 = total for atoms and orbitals
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS_iatom; // vDOS.at(iatom).at(orbital).at(spin).at(energy_number); 0 = total for atoms and orbitals
  std::deque<std::deque<std::deque<std::deque<double>>>> vDOS_lm_iatom; // vDOS.at(iatom).at(orbital_lm).at(spin).at(energy_number); 0 = total for atoms and orbitals
  //SD20230214 END
  void GetVDOSSpecies(const xstructure& xstr); //vDOS.at(species).at(spin).at(energy_number)  //CO20191110 //SD20230213 - void function
  void GetVDOSSpecies(const std::deque<int>& num_each_type); //vDOS.at(species).at(spin).at(energy_number)  //CO20191110 //SD20230213 - void function
  void GetVDOSIAtom(const xstructure& xstr); //vDOS.at(iatom).at(spin).at(energy_number)  //SD20230213
  void GetVDOSIAtom(const std::vector<std::vector<int>>& iatoms_index); //vDOS.at(iatom).at(spin).at(energy_number)  //SD20230213
  friend std::ostream& operator<<(std::ostream&, const xDOSCAR&); //ME20190623

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xDOSCAR deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xDOSCAR"; }

public:
  // debugging methods:
  void printVidos() const { std::cout << viDOS.size() << std::endl; }

// SERIALIZATION MEMBERS, ignoring vcontent and content
#define JSON_xDOSCAR_MEMBERS                                                                                                                                                                              \
  m_initialized, filename, title, spin, Vol, POTIM, lattice, temperature, RWIGS, Efermi, spinF, energy_max, energy_min, number_energies, number_atoms, partial, denergy, venergy, venergyEf, viDOS, vDOS, \
      isLSCOUPLING, lmResolved, carstring, conduction_band_min, conduction_band_min_net, valence_band_max, valence_band_max_net, Egap, Egap_net, Egap_fit, Egap_fit_net, Egap_type, Egap_type_net
};
//-------------------------------------------------------------------------------------------------
class xEIGENVAL : public xOVASP<xEIGENVAL>, public JsonSerializable<xEIGENVAL> {
public:
  using xOVASP::xOVASP;

  xEIGENVAL(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xEIGENVAL(fromFile(fileIN, QUIET, oss)) {}
  xEIGENVAL(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xEIGENVAL(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xEIGENVAL(); }

  std::string title;
  uint number_atoms = 0;
  uint number_loops = 0;
  uint spin = 0;
  double Vol = 0.0;
  double POTIM = 0.0;
  aurostd::xvector<double> lattice;
  double temperature = 0.0;
  uint number_electrons = 0;
  uint number_kpoints = 0;
  uint number_bands = 0;
  std::deque<double> vweight; // vweight.at(kpoint number)
  std::deque<aurostd::xvector<double>> vkpoint; // vkpoint.at(kpoint number)[1,2,3]=xyz.
  std::deque<std::deque<std::deque<double>>> venergy; // venergy.at(kpoint number).at(band number).at(spin number)
  std::string carstring; // ME20190620 - the fourth line of the EIGENVAL file
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
  double energy_max = 0.0;
  double energy_min = 0.0;
  friend std::ostream& operator<<(std::ostream&, const xEIGENVAL&); // ME20190623

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xEIGENVAL deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xEIGENVAL"; }

// SERIALIZATION MEMBERS
#define JSON_xEIGENVAL_MEMBERS \
  m_initialized, filename, title, number_atoms, number_loops, spin, Vol, POTIM, lattice, temperature, number_electrons, number_kpoints, number_bands, vweight, vkpoint, venergy, carstring, energy_max, energy_min
};
//-------------------------------------------------------------------------------------------------
class xPOTCAR : public xOVASP<xPOTCAR>, public JsonSerializable<xPOTCAR> {
public:
  using xOVASP::xOVASP;

  xPOTCAR(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xPOTCAR(fromFile(fileIN, QUIET, oss)) {}
  xPOTCAR(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xPOTCAR(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xPOTCAR(); }; // clear

  std::string title;
  bool POTCAR_PAW = false;
  std::string POTCAR_TYPE;
  bool POTCAR_KINETIC = false;
  bool POTCAR_GW = false;
  bool POTCAR_AE = false;
  double ENMAX = 0.0;
  std::vector<double> vENMAX; // eV
  double ENMIN = 0.0;
  std::vector<double> vENMIN; // eV
  double POMASS_sum = 0.0;
  double POMASS_min = 0.0;
  double POMASS_max = 0.0;
  std::vector<double> vPOMASS; // mass
  double ZVAL_sum = 0.0;
  double ZVAL_min = 0.0;
  double ZVAL_max = 0.0;
  std::vector<double> vZVAL; // valence
  double EATOM_min = 0.0;
  double EATOM_max = 0.0;
  std::vector<double> vEATOM; // eV
  double RCORE_min = 0.0;
  double RCORE_max = 0.0;
  std::vector<double> vRCORE; // outmost cutoff radius
  double RWIGS_min = 0.0;
  double RWIGS_max = 0.0;
  std::vector<double> vRWIGS; // wigner-seitz radius (au A)
  double EAUG_min = 0.0;
  double EAUG_max = 0.0;
  std::vector<double> vEAUG; // augmentation
  double RAUG_min = 0.0;
  double RAUG_max = 0.0;
  std::vector<double> vRAUG; // augmentation
  double RMAX_min = 0.0;
  double RMAX_max = 0.0;
  std::vector<double> vRMAX; // unicity
  std::vector<std::string> vTITEL; // unicity
  std::vector<std::string> vLEXCH; // unicity
  std::string pp_type;
  std::vector<std::string> species; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<int> species_Z; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_type; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_version; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_AUID; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<std::string> species_pp_AUID_collisions; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::vector<double> species_pp_groundstate_energy; // meV/atom
  std::vector<std::string> species_pp_groundstate_structure; // name that we have, maybe ANRL
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
  // objects/functions for references energies defined only with one specie
  std::string AUID; // crc32 - THIS IS A GLOBAL PROPERTY OF THE WHOLE POTCAR
  friend std::ostream& operator<<(std::ostream&, const xPOTCAR&); // print // FIX COPY CONSTRUCTOR
  // xPOTCAR xPOTCAR_initialize(uint Z);                         // function to clean up the name // FIX COPY CONSTRUCTOR

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xPOTCAR deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xPOTCAR"; }

  // SERIALIZATION MEMBERS
#define JSON_xPOTCAR_MEMBERS                                                                                                                                                                                    \
  m_initialized, filename, title, POTCAR_PAW, POTCAR_TYPE, POTCAR_KINETIC, POTCAR_GW, POTCAR_AE, vENMAX, vENMIN, vPOMASS, vZVAL, vEATOM, vRCORE, vRWIGS, vEAUG, vRAUG, vRMAX, vTITEL, vLEXCH, pp_type, species, \
      species_Z, species_pp, species_pp_type, species_pp_version, species_pp_AUID, species_pp_AUID_collisions, species_pp_groundstate_energy, species_pp_groundstate_structure
};

std::vector<xPOTCAR> get_pseudopotential_data();
aurostd::JSON::object get_pseudopotential_enthalpy_references();
bool xPOTCAR_PURE_Printer(xPOTCAR& xPOT, std::ostream& oss, bool LVERBOSE = false);
xPOTCAR xPOTCAR_Finder(std::vector<std::string>& species_pp_AUID, std::vector<std::string>& species_pp_AUID_collisions, const std::string& TITEL, const std::string& LEXCH, const double& EATOM, const double& RMAX, bool LVERBOSE = false);
xPOTCAR xPOTCAR_Finder(const std::string& AUID, bool LVERBOSE = false);
bool xPOTCAR_EnthalpyReference_AUID(std::string AUID, std::string METAGGA = ""); // returns if available
bool xPOTCAR_EnthalpyReference_AUID(std::string AUID, std::string METAGGA, std::string& gs, double& enthalpy_atom, double& volume_atom, double& spin_atom);

// -------------------------------------------------------------------------------------------------
class xVASPRUNXML : public xOVASP<xVASPRUNXML>, public JsonSerializable<xVASPRUNXML> { // CO20200404 - xStream integration for logging
public:
  using xOVASP::xOVASP;

  xVASPRUNXML(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xVASPRUNXML(fromFile(fileIN, QUIET, oss)) {}
  xVASPRUNXML(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xVASPRUNXML(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xVASPRUNXML(); }; // clear

  double natoms = 0.0; // for aflowlib_libraries.cpp
  aurostd::xmatrix<double> stress; // for aflowlib_libraries.cpp
  std::vector<aurostd::xvector<double>> vkpoint; // for aflowlib_libraries.cpp
  std::vector<aurostd::xvector<double>> vweights; // for aflowlib_libraries.cpp
  std::vector<aurostd::xvector<double>> vforces; // for aflowlib_libraries.cpp
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
  bool GetForces(const std::string&, bool = true); // ME20190204
  bool GetForcesFile(const std::string&, bool = true); // ME20190204
  bool GetForces(std::stringstream&, bool = true); // ME20190204

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xVASPRUNXML deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xVASPRUNXML"; }

// SERIALIZATION MEMBERS
#define JSON_xVASPRUNXML_MEMBERS m_initialized, filename, natoms, stress, vkpoint, vweights, vforces
};
// -------------------------------------------------------------------------------------------------
class xIBZKPT : public xOVASP<xIBZKPT>, public JsonSerializable<xIBZKPT> { // CO20200404 - xStream integration for logging
public:
  using xOVASP::xOVASP;

  xIBZKPT(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xIBZKPT(fromFile(fileIN, QUIET, oss)) {}
  xIBZKPT(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xIBZKPT(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xIBZKPT(); }; // clear

  uint nweights = 0; // for aflowlib_libraries.cpp
  uint nkpoints_irreducible = 0; // for aflowlib_libraries.cpp
  std::vector<aurostd::xvector<double>> vkpoint; // for aflowlib_libraries.cpp
  std::vector<uint> vweights; // for aflowlib_libraries.cpp
  uint ntetrahedra = 0; // for aflowlib_libraries.cpp
  double wtetrahedra = 0.0; // for aflowlib_libraries.cpp
  std::vector<aurostd::xvector<int>> vtetrahedra; // for aflowlib_libraries.cpp
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xIBZKPT deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xIBZKPT"; }

// SERIALIZATION MEMBERS
#define JSON_xIBZKPT_MEMBERS m_initialized, filename, nweights, nkpoints_irreducible, vkpoint, vweights, ntetrahedra, wtetrahedra, vtetrahedra
};
// -------------------------------------------------------------------------------------------------
class xKPOINTS : public xOVASP<xKPOINTS>, public JsonSerializable<xKPOINTS> { // CO20200404 - xStream integration for logging
public:
  using xOVASP::xOVASP;

  xKPOINTS(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xKPOINTS(fromFile(fileIN, QUIET, oss)) {}
  xKPOINTS(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xKPOINTS(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xKPOINTS(); }; // clear

  std::string title; // first line
  int mode; // sort of mode
  std::string grid_type; // if grid specified
  bool is_KPOINTS_NNN, is_KPOINTS_PATH; // control parameters
  aurostd::xvector<int> nnn_kpoints; // N*N*N                          // triplet of kpoints
  aurostd::xvector<double> ooo_kpoints; // ORIGIN                         // triplet of origin
  int nkpoints; // total kpoints
  std::string path_mode, path;
  std::vector<std::string> vpath;
  int path_grid; // path if any
  std::vector<aurostd::xvector<double>> vkpoints; // ME20190614 - k-point coordinates of the path

  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
  friend std::ostream& operator<<(std::ostream&, const xKPOINTS&); // ME20190623
  std::string createStandardTitlePath(const xstructure&); // ME20190623

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xKPOINTS deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xKPOINTS"; }

// SERIALIZATION MEMBERS
#define JSON_xKPOINTS_MEMBERS title, filename, mode, grid_type, is_KPOINTS_NNN, is_KPOINTS_PATH, nnn_kpoints, ooo_kpoints, nkpoints, path_mode, path, vpath, path_grid, vkpoints, m_initialized
};
// -------------------------------------------------------------------------------------------------
class xCHGCAR : public xOVASP<xCHGCAR> { // CO20200404 - xStream integration for logging
public:
  using xOVASP::xOVASP;

  xCHGCAR(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xCHGCAR(fromFile(fileIN, QUIET, oss)) {}
  xCHGCAR(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xCHGCAR(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xCHGCAR(); }; // clear

  aurostd::xvector<int> grid; // N*N*N                                // triplet of grid
  std::vector<std::string> vstring; // ORIGIN                            // std::string of values
  aurostd::xvector<double> vvalues; // ORIGIN                            // xvector of values
  aurostd::xtensor<double> tvalues; // ORIGIN                             // xtensor of values ME20180705
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
};
// -------------------------------------------------------------------------------------------------
class xQMVASP : public xOVASP<xQMVASP>, public JsonSerializable<xQMVASP> { // CO20191110 //CO20200404 - xStream integration for logging
public:
  using xOVASP::xOVASP;

  xQMVASP(const std::string& fileIN, bool QUIET = true, std::ostream& oss = std::cout) : xQMVASP(fromFile(fileIN, QUIET, oss)) {}
  xQMVASP(const std::string& fileIN, std::ofstream& FileMESSAGE, bool QUIET = true, std::ostream& oss = std::cout) : xQMVASP(fromFile(fileIN, FileMESSAGE, QUIET, oss)) {}

  void clear() override { *this = xQMVASP(); }; // clear

  double H_atom_relax = AUROSTD_NAN;
  double H_atom_static = AUROSTD_NAN;
  std::vector<aurostd::xvector<double>> vforces; // for APL - only one (no relax vs. static), get most relaxed forces  //CO20191112
  xstructure xstr_final;
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xQMVASP deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xQMVASP"; }

// SERIALIZATION MEMBERS
#define JSON_xQMVASP_MEMBERS m_initialized, filename, H_atom_relax, H_atom_static, vforces
};
// -------------------------------------------------------------------------------------------------
// aflow_kaims.cpp
namespace KBIN {
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(std::string& AflowIn, _aflags& aflags, _kflags& kflags);
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(std::string& AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
  bool AIMS_Directory(std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);
} // namespace KBIN
// -------------------------------------------------------------------------------------------------
// aflow_iaims.cpp
namespace KBIN {
  bool AIMS_Produce_INPUT(_xaims& xaims, std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _aimsflags& aimsflags);
  bool AIMS_Modify_INPUT(_xaims& xaims, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _aimsflags& aimsflags);
  bool AIMS_Write_INPUT(_xaims& xaims, _aimsflags& aimsflags);
  bool AIMS_Write_CONTROL(_xaims& xaims, _aimsflags& aimsflags);
  bool AIMS_Write_GEOM(_xaims& xaims, _aimsflags& aimsflags);
  bool AIMS_Produce_CONTROL(_xaims& xaims, std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _aimsflags& aimsflags);
  bool AIMS_Modify_CONTROL(_xaims& xaims, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _aimsflags& aimsflags);
  bool AIMS_Reread_CONTROL(_xaims& xaims, std::ofstream& FileMESSAGE, _aflags& aflags);
  bool AIMS_Produce_GEOM(_xaims& xaims, std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, _aimsflags& aimsflags);
  bool AIMS_Produce_GEOM(_xaims& xaims);
  bool AIMS_Modify_GEOM(_xaims& xaims, std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _aimsflags& aimsflags);
  bool AIMS_Reread_GEOM(_xaims& xaims, std::ofstream& FileMESSAGE, _aflags& aflags);
  bool XAIMS_CONTROL_PREPARE_GENERIC(std::string command, _xaims& xaims, _aimsflags& aimsflags, std::string svalue, int ivalue, double dvalue, bool OPTION);
  void XAIMS_CONTROL_REMOVE_ENTRY(_xaims& xaims, std::string ENTRY, std::string COMMENT, bool VERBOSE);
} // namespace KBIN
// -------------------------------------------------------------------------------------------------
// aflow_oaims.cpp
class xAIMSOUT;
/// @brief Abstract base class for the AIMS data containers such as @c xAIMSOUT.
/// Implements some convenience static methods to construct derived classes from files. Derived classes must
/// implement @c GetProperties to read from string, @c GetPropertiesFile to read from file, and @c clear to reset the object to default.
/// Derived class may use @code using xOAIMS::xOAIMS; @endcode to inherit constructors.
template <class Derived> class xOAIMS {
public:
  virtual ~xOAIMS() = default;
  xOAIMS() = default;

  static Derived fromFile(const std::string& fileIN, bool QUIET = true) {
    Derived derived;
    derived.filename = fileIN;
    derived.GetPropertiesFile(fileIN, QUIET);
    return derived;
  }

  virtual bool GetProperties(const std::string& stringIN, bool QUIET = true) = 0;
  virtual bool GetPropertiesFile(const std::string& fileIN, bool QUIET = true) = 0;

  virtual void clear() = 0;

  // CONTENT
  std::string content;
  std::vector<std::string> vcontent;
  std::string filename;
};
class xAIMSOUT : public xOAIMS<xAIMSOUT> {
public:
  using xOAIMS::xOAIMS;

  xAIMSOUT(const std::string& fileIN, bool QUIET = true) : xAIMSOUT(fromFile(fileIN, QUIET)) {}

  void clear() override { *this = xAIMSOUT(); }

  std::vector<aurostd::xvector<double>> vforces; // for aflowlib_libraries.cpp
  double natoms = 0;
  bool GetProperties(const std::stringstream& stringstreamIN, bool = true); // get everything QUIET
  bool GetProperties(const std::string& stringIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, bool = true) override; // get everything QUIET
  bool GetPropertiesFile(const std::string& fileIN, uint natoms_check, bool); // get everything QUIET
  bool GetPropertiesUrlFile(const std::string& url, const std::string& file, bool = true); // get everything from an aflowlib entry
};
// -----------------------------------------------------------------------------------------------
bool PrintBandGap(std::string& WorkDir, std::ostream& oss);
bool PrintBandGap_DOS(std::string& WorkDir, std::ostream& oss); // CO20191110
bool PrintEffectiveMass(std::string& WorkDir, std::ostream& oss);
bool PrintEigCurv(std::string& WorkDir, std::ostream& oss);
// -----------------------------------------------------------------------------------------------
bool ParseKPOINTS(std::stringstream& File_Kpoints, int& GRIDS, std::vector<aurostd::xvector<double>>& special_kpts, std::vector<aurostd::xvector<double>>& unique_kpts, std::vector<int>& repeat_kpts_num);
bool AdjacencyList_KPT(std::vector<aurostd::xvector<double>>& special_kpts, std::vector<aurostd::xvector<double>>& unique_kpts, std::vector<aurostd::xvector<int>>& connect_kpts, std::vector<int>& connect_kpts_num);
bool AdjacencyList_EIG(std::vector<aurostd::xvector<double>>& unique_kpts,
                       std::vector<aurostd::xvector<int>>& connect_kpts,
                       std::vector<int>& connect_kpts_num,
                       xEIGENVAL& xeigenval,
                       std::vector<aurostd::xvector<double>>& unique_kpts_EIG,
                       std::vector<aurostd::xvector<int>>& connect_kpts_EIG,
                       std::vector<aurostd::xvector<double>>& vkpoint_eig);
bool RepeatsList(std::vector<aurostd::xvector<double>>& unique_kpts_EIG, std::vector<int>& repeat_kpts_num, std::vector<aurostd::xvector<double>>& vkpoint_eig, std::vector<aurostd::xvector<int>>& repeat_kpts_EIG);
bool VertexPaths(std::vector<aurostd::xvector<int>>& repeat_kpts_EIG, std::vector<aurostd::xvector<int>>& connect_kpts_EIG, std::vector<int>& repeat_kpts_num, int& GRIDS, std::vector<aurostd::xvector<int>>& vrtx_path);
bool RepeatedEdges(std::vector<aurostd::xvector<int>>& vrtx_path, std::vector<aurostd::xvector<int>>& repeat_kpts_EIG, std::vector<int>& repeat_kpts_num, std::vector<aurostd::xvector<int>>& ndx_edges);
bool VertexBranches(std::vector<aurostd::xvector<int>>& ndx_edges, std::vector<int>& repeat_kpts_num, std::vector<aurostd::xvector<int>>& repeat_kpts_EIG, std::vector<std::vector<aurostd::xvector<int>>>& branches);
bool PathDataStuct(xEIGENVAL& xeigenval,
                   std::vector<aurostd::xvector<double>>& vkpoint_eig,
                   std::vector<std::vector<aurostd::xvector<int>>>& branches,
                   std::vector<std::vector<std::vector<int>>>& branches_indx,
                   std::vector<std::vector<std::vector<aurostd::xvector<double>>>>& branches_kpts,
                   std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& branches_bnds);
bool IBZextrema(xEIGENVAL& xeigenval, std::vector<aurostd::xvector<double>>& vkpoint_eig, std::vector<std::vector<aurostd::xvector<int>>>& branches);
void CompareDoublesChar(bool& MATCH, double& number1, double& number2);
void CompareEdges(std::vector<std::vector<aurostd::xvector<int>>>& branches, std::vector<aurostd::xvector<int>>& vertex_edges, aurostd::xvector<int>& test_edge, bool& MATCH);
void NaiveCurvatures(aurostd::xvector<double>& eigvec, std::vector<aurostd::xvector<double>>& posvec, std::vector<double>& curvature);
double StencilLinear1D(std::vector<aurostd::xvector<double>>& positions, aurostd::xvector<double>& eigenvals);
//-------------------------------------------------------------------------------------------------
struct kEn_st {
  aurostd::xvector<double> kpoint;
  double energy[2];
  int band_index;
  int band_type; // 0 -- valence band; 1 -- conduction band
};
#define _SIGMA 1.0 // default standard deviation of input data
// range of energy point to fit the ellipse curve
const double _FIT_ENERGY_RANGE = 0.026; // eV range of band
const int _FIT_POINTS_NUMBER = 8; // minimum fit points in Irreducible BZ
// range of band extremes to determine the number of bands for effective mass calculations
const double _BANDS_ENERGY_RANGE = 0.026; // eV
// used to determine cluster of points can be changed to other values
const double _BANDS_PARAMETER_MIN_RATIO = 0.2;
// factor unit
// mass is in unit of electron mass
const double _MASS_FACTOR = 3.80998; // hbar^2*10^{20}/(2.0*me*eV)
bool comparison_kEn_str_up(const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_dn(const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_position(const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_band_type_up(const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_band_type_dn(const kEn_st& k1, const kEn_st& k2);
bool is_equal_position_kEn_str(const kEn_st& k1, const kEn_st& k2);
bool near_to(const aurostd::xvector<double>& k1, const aurostd::xvector<double>& k2, const std::vector<double>& max_distance);
namespace aurostd {
  class JSONwriter; // forward-declaration of JSONwriter class: later in plotter
  // namespace JSONwriter class defined in aurostd.h is not visible; dependencies race?
} // namespace aurostd
//-------------------------------------------------------------------------------------------------
// ME20190614 - plotter functions
namespace plotter {
  // Plot setup --------------------------------------------------------------
  // Plot options
  aurostd::xoption getPlotOptions(const aurostd::xoption&, const std::string&, bool = false);
  aurostd::xoption getPlotOptionsEStructure(const aurostd::xoption&, const std::string&, bool = false);
  aurostd::xoption getPlotOptionsPhonons(const aurostd::xoption&, const std::string&);
  aurostd::xoption getPlotOptionsQHAthermo(const aurostd::xoption& xopt, const std::string& key); // AS20210705

  // Plot functions
  void generateHeader(std::stringstream&, const aurostd::xoption&, bool = false);
  void savePlotGNUPLOT(const aurostd::xoption&, const std::stringstream&);
  void setFileName(aurostd::xoption&, std::string = "");
  void setTitle(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void setTitle(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  std::string formatDefaultPlotTitle(const aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  std::string formatDefaultPlotTitle(const aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  std::vector<double> getCompositionFromHTQCPrototype(const std::string&, const std::string&); // ME20190813
  std::vector<double> getCompositionFromANRLPrototype(const std::string&);
  std::string formatDefaultTitlePOCC(const aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  std::string formatDefaultTitlePOCC(const aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  std::vector<double> getCompositionFromPoccString(const std::string&, bool&);

  // Electronic structure ----------------------------------------------------
  void patchDefaultTitleAFLOWIN(aurostd::xoption& plotoptions); // CO20191110
  // Plot functions
  void PLOT_DOS(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_DOS(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_DOS(aurostd::xoption&, const xDOSCAR&, std::ostream& oss = std::cout); // CO20191110 //CO20200404
  void PLOT_DOS(aurostd::xoption&, const xDOSCAR&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20191110 //CO20200404
  void PLOT_DOS(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_DOS(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_DOS(aurostd::xoption&, std::stringstream&, const xDOSCAR&, std::ostream& oss = std::cout); // CO20191110  //CO20200404
  void PLOT_DOS(aurostd::xoption&, std::stringstream&, const xDOSCAR&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20191110  //CO20200404

  void PLOT_PDOS(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PDOS(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PDOS(aurostd::xoption&, const xDOSCAR&, std::ostream& oss = std::cout); // CO20191110 //CO20200404
  void PLOT_PDOS(aurostd::xoption&, const xDOSCAR&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20191110 //CO20200404
  void PLOT_PDOS(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PDOS(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PDOS(aurostd::xoption&, std::stringstream&, const xDOSCAR&, std::ostream& oss = std::cout); // CO20191110 //CO20200404
  void PLOT_PDOS(aurostd::xoption&, std::stringstream&, const xDOSCAR&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20191110 //CO20200404

  void PLOT_BAND(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_BAND(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_BAND(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_BAND(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void BANDDOS2JSON(std::ostream&, std::string);
  void PLOT_BANDDOS(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_BANDDOS(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_BANDDOS(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_BANDDOS(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404

  // Helper functions
  xstructure getStructureWithNames(const aurostd::xoption&, const std::string& carstring = "CAR", std::ostream& oss = std::cout); // CO20191110 //CO20200404
  xstructure getStructureWithNames(const aurostd::xoption&, std::ofstream& FileMESSAGE, const std::string& carstring = "CAR", std::ostream& oss = std::cout); // CO20191110 //CO20200404
  std::string getLatticeFromKpointsTitle(const std::string&);
  void shiftEfermiToZero(xEIGENVAL&, double);
  void setEMinMax(aurostd::xoption&, double, double);
  aurostd::JSON::object DOS2JSON(aurostd::xoption& xopt, const xDOSCAR& xdos, std::ofstream& FileMESSAGE,
                                 std::ostream& oss); // AS20201102
  aurostd::JSON::object bands2JSON(const xEIGENVAL& xeigen,
                                   const xKPOINTS& xkpts,
                                   const std::vector<double>& distances,
                                   const std::vector<double>& segment_points,
                                   const aurostd::xoption& plotoptions); // AS2021102
  aurostd::JSON::object bandsDOS2JSON(const xDOSCAR& xdos, const xEIGENVAL& xeigen, const xKPOINTS& xkpts, aurostd::xoption& xopt, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // AS20201102  //ME20211014 - added default for oss

  // DOS
  bool dosDataAvailable(const std::deque<std::deque<std::deque<std::deque<double>>>>& vdos, int pdos); // ME20200305
  void generateDosPlot(std::stringstream&, const xDOSCAR&, aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void generateDosPlot(std::stringstream&, const xDOSCAR&, aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404

  // Bands
  void generateBandPlot(std::stringstream&, const xEIGENVAL&, const xKPOINTS&, const xstructure&, const aurostd::xoption&);
  void generateBandPlot(std::stringstream&, aurostd::JSON::object&, const xEIGENVAL&, const xKPOINTS&, const xstructure&, const aurostd::xoption&);
  std::string convertKPointLabel(const std::string&, const std::string&);
  std::string convertKPointLetter(std::string, const std::string&);

  // Gnuplot
  void generateDosPlotGNUPLOT(std::stringstream&, const xDOSCAR&, const std::deque<double>&, const std::deque<std::deque<std::deque<double>>>&, const std::vector<std::string>&, const aurostd::xoption&);
  double getDosLimits(const aurostd::xoption&, const xDOSCAR&, const std::deque<std::deque<std::deque<double>>>&, const std::deque<double>&);
  void generateBandPlotGNUPLOT(std::stringstream&, const xEIGENVAL&, const std::vector<double>&, const std::vector<double>&, const std::vector<std::string>&, const aurostd::xoption&);
  std::string getFormattedUnit(const std::string&);

  // Phonons -----------------------------------------------------------------
  void PLOT_PHDOS(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDOS(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDOS(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDOS(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDOS(aurostd::xoption&, const xDOSCAR&, std::ostream& oss = std::cout); // ME20210927
  void PLOT_PHDOS(aurostd::xoption&, const xDOSCAR&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // ME20210927
  void PLOT_PHDOS(aurostd::xoption&, std::stringstream& out, xDOSCAR, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // ME20210927

  void PLOT_PHDISP(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDISP(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDISP(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDISPDOS(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDISPDOS(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_PHDISPDOS(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20204004
  void PLOT_PHDISPDOS(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20204004

  void convertEnergies(xEIGENVAL&, const std::string&);
  void convertEnergies(xDOSCAR&, const std::string&);
  double getEnergyConversionFactor(const std::string&);

  // Properties plotter ------------------------------------------------------
  void PLOT_THERMO(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_THERMO(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_THERMO(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_THERMO(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_TCOND(aurostd::xoption&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_TCOND(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void PLOT_TCOND(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void PLOT_TCOND(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404

  // QHA properties plotter -------------------------------------------------
  void PLOT_THERMO_QHA(aurostd::xoption&, std::ostream& oss = std::cout); // AS20200909
  void PLOT_THERMO_QHA(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // AS20200909
  void PLOT_THERMO_QHA(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // AS20200909
  void PLOT_THERMO_QHA(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // AS20200909
  void PLOT_GRUENEISEN_DISPERSION(aurostd::xoption&, std::ostream& oss = std::cout); // AS20210701
  void PLOT_GRUENEISEN_DISPERSION(aurostd::xoption&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // AS20210701
  void PLOT_GRUENEISEN_DISPERSION(aurostd::xoption&, std::stringstream&, std::ostream& oss = std::cout); // AS20210701
  void PLOT_GRUENEISEN_DISPERSION(aurostd::xoption&, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // AS20210701

  // General plots -----------------------------------------------------------
  void plotSingleFromSet(aurostd::xoption&, std::stringstream&, const std::vector<std::vector<double>>&, int, std::ostream& oss = std::cout); // CO20200404
  void plotSingleFromSet(aurostd::xoption&, std::stringstream&, const std::vector<std::vector<double>>&, int, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void plotMatrix(aurostd::xoption& plotoptions, std::stringstream&, std::ostream& oss = std::cout); // CO20200404
  void plotMatrix(aurostd::xoption& plotoptions, std::stringstream&, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200404
  void setPlotLabels(aurostd::xoption&, const std::string&, const std::string&, const std::string&, const std::string&);
  std::vector<std::vector<double>> readAflowDataFile(aurostd::xoption&);
  void generatePlotGNUPLOT(std::stringstream&, const aurostd::xoption&, const std::vector<std::vector<double>>&);
} // namespace plotter

//-------------------------------------------------------------------------------------------------
// aflow_estructure_dos.cpp

namespace estructure {
  std::string PEDOS_GENERATE_GNUPLOTSCRIPT(const std::string&, const std::string&, const double&, const double&, const double&, const double&, const int&, const std::vector<std::vector<std::vector<double>>>&, const std::string&);
  bool isSpecialKPOINT(std::string kpoint); // CO20170830
  std::string fixSpecialKPOINT_GNUPLOT(std::string kpoint, bool json = false); // CO20170830
  std::string fixSpecialKPOINT_HTML(std::string kpoint); // CO20170830
  std::string fixSpecialKPOINT_LATEX(std::string kpoint); // CO20170830
  std::string fixKPOINT_GNUPLOT(std::string kpoint, bool json = false); // CO20170830
  std::string fixKPOINT_HTML(std::string kpoint); // CO20170830
  std::string fixKPOINT_LATEX(std::string kpoint); // CO20170830
  std::string fixKPOINT_SPECIALONLY(std::string kpoint); // CO20170830
  void PLOT_BANDDOS(std::string options);
  void PLOT_BAND(std::string options);
  void PLOT_DOS(std::string options);
  void PLOT_PEDOS(std::string options);
  void PLOT_PEDOSALL(std::string options);
  void PLOT_PEDOSALL_AFLOWLIB(std::string options, _aflags& aflags);
  void PLOT_BAND2(std::string options);
  void PLOT_BAND_SPINSPLIT(std::string options);
  void PLOT_DOSWEB(std::string options);
  // manipulation
  std::string changeICSDNameGunplot(std::string ICSDName);
  void CombineTDOSAndTOTALPDOS(const std::vector<std::vector<double>>& TDOS, const std::vector<std::vector<double>>& TOTALPDOS, std::vector<std::vector<double>>& vvDOS);
  double GET_TDOSDATA(const std::string& str_dir, std::vector<std::vector<double>>& TDOS);
  double GET_TDOSDATA(std::stringstream& ss_dosfile, std::stringstream& ss_outcarfile, std::vector<std::vector<double>>& TDOS);
  double GET_TOTALPDOSDATA(const std::string& str_dir, std::vector<std::vector<double>>& TOTALPDOS);
  double GET_TOTALPDOSDATA(std::stringstream& ss_dosfile, std::stringstream& ss_outfile, std::vector<std::vector<double>>& TOTALPDOS);
  double GET_PDOSDATA(const std::string& str_dir, std::vector<std::vector<std::vector<double>>>& PDOS);
  double GET_PDOSDATA(std::stringstream& ss_dosfile, std::stringstream& ss_outfile, std::vector<std::vector<std::vector<double>>>& PDOS);
  bool GET_DOS_DATA(std::stringstream& ss_dosfile, std::stringstream& ss_outfile, double& Efermi, std::vector<std::vector<double>>& TDOS, std::vector<std::vector<double>>& TOTALPDOS); // CO20180216
  bool GET_DOS_DATA(const std::string& str_dir, double& Efermi, std::vector<std::vector<double>>& TDOS, std::vector<std::vector<double>>& TOTALPDOS, std::vector<std::vector<std::vector<double>>>& PDOS); // CO20180216
  bool GET_DOS_DATA(std::stringstream& ss_dosfile,
                    std::stringstream& ss_outfile,
                    double& Efermi,
                    std::vector<std::vector<double>>& TDOS,
                    std::vector<std::vector<double>>& TOTALPDOS,
                    std::vector<std::vector<std::vector<double>>>& PDOS); // CO20180216
  void FormatSpinofPDOS(std::vector<std::vector<std::vector<double>>>& vvva);

  // Functions for serializing bands data to JSON
  // Added by EG
  bool DOSDATA_JSON(aurostd::xoption& vpflow, std::ostream& oss = std::cout);
  bool DOSDATA_JSON(aurostd::xoption& vpflow, std::string directory, std::stringstream& json, bool wrapping_brackets = true);
  bool BANDSDATA_JSON(aurostd::xoption& vpflow, std::ostream& oss = std::cout);
  bool BANDSDATA_JSON(aurostd::xoption& vpflow, std::string directory, std::stringstream& json, bool wrapping_brackets = true);
  // uint DOSDATA_JSON(std::string options);
  // uint DOSDATA_JSON(std::string options, std::ostream& json);
  // uint BANDSDATA_JSON(std::string options);
  // uint BANDSDATA_JSON(std::string options, std::string json_dir);
  // uint BANDSDATA_JSON(std::string options, std::ostream& json);
  std::string linelabel2HTML(std::string linelabel);
  uint inequivalentAtomsJSON(std::vector<std::vector<std::vector<double>>>& PDOS, std::vector<int>& iatoms, std::vector<double>& numbers, std::vector<std::string>& vspecies, std::ostream& json);
  uint constructInequivalentAtomPDOSJSON(std::vector<std::vector<std::vector<double>>>& PDOS, int iatom, std::ostream& json);
  // End of bands data JSON serializers

} // namespace estructure

// ----------------------------------------------------------------------------
// aflow_poccupation_*.cpp
//  #include "modules/POCC/aflow_pocc.h"

// aflow_poccupation_params.cpp
namespace pocc {
  bool poccInput(); // CO20170805

  std::string ReturnAtomSpecies(std::string atom);
  std::string ReturnAtomSpeciesPotential(std::string atom);
  std::string ReturnUFFParameters(std::string atom);
  class UFFPara {
  public:
    UFFPara(); // constructor
    ~UFFPara(); // destructor
    std::string symbol;
    double r1, theta0, x1, D1, zeta, Z1, Vi, Uj, Xi, hard, radius;
    void GetUFFParameters(std::string);

  private:
    void free(); // free space
  };
  std::string ReturnAtomProperties(std::string atom);
  // Atomic Properties Database
  class Atom {
  public:
    Atom();
    ~Atom();
    std::string name, symbol;
    int number; // atomic number
    double mass, radius, Xi; // atomic, weight radius /pauling electronegativity
    void GetAtomicProperties(std::string);

  private:
    void free();
  };

  // aflow_poccupation_forcefield.cpp
  class Bond {
  public:
    Bond();
    Bond(const Bond& b);
    ~Bond();
    _atom bgn, end;
    double length;
    void Set(xstructure, _atom, _atom);
    const Bond& operator=(const Bond& other);
    bool operator==(const Bond& other) const;
    bool operator!=(const Bond& other) const;
    friend std::ostream& operator<<(std::ostream&, const Bond&);

  private:
    void free();
    void copy(const Bond& b);
  };
  void SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij);
  double CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj);
  double CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj);
  double CalculateUFFEnergy(xstructure xstr);
  void RemoveSameBond(std::vector<Bond>& Bonds_orig, std::vector<Bond>& Bonds_new);
  void ExtractBonds(const xstructure& xstr, std::deque<std::deque<_atom>>& neigh_mat_bonded, std::deque<std::deque<_atom>>& neigh_mat_nonbonded);
  void AnalyzeBonds(const xstructure& xstr, std::vector<Bond>& Bonds, std::vector<Bond>& NonBonds);
  void UFFENERGY(std::istream& input);
} // namespace pocc

// ----------------------------------------------------------------------------
// aflow_mix.cpp  aflow_nomix.cpp   aflow_mix_pauling.cpp

int MiscibilityCheck(int speciesA, int speciesB); // aflow_mix.cpp
int MiscibilityCheck(std::string speciesA, std::string speciesB); // aflow_mix.cpp
int MiscibilityExperimentsCheck(int speciesA, int speciesB); // aflow_mix.cpp
int MiscibilityExperimentsCheck(std::string speciesA, std::string speciesB); // aflow_mix.cpp
int MiscibilityMiedemaCheck(int speciesA, int speciesB); // aflow_mix.cpp
int MiscibilityMiedemaCheck(std::string speciesA, std::string speciesB); // aflow_mix.cpp
int MiscibilityMiedemaCheck(std::string system_in); // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(int speciesA, int speciesB); // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(std::string speciesA, std::string speciesB); // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(std::string system_in); // aflow_mix.cpp
int MiscibilityCheck(std::string system_in); // aflow_nomix.cpp
int MiscibilityExperimentsCheck(std::string system_in); // aflow_mix_pauling.cpp

// ----------------------------------------------------------------------------
// neighbors prototypes
// aflow_neighbors.cpp
bool StepNeighborsPerform(xstructure& a, std::string AflowIn, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags);

// ----------------------------------------------------------------------------
// aflow_pocc //CO20180502
namespace KBIN {
  void VASP_RunPOCC(const std::string& directory, std::ostream& oss = std::cout); // CO20200624
  void VASP_RunPOCC(const std::string& directory, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20200624
  void VASP_RunPOCC(const _xvasp& xvasp, const std::string& AflowIn, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
} // namespace KBIN
// ----------------------------------------------------------------------------
// aflow_phonons.cpp
namespace KBIN {
  bool relaxStructureAPL_VASP(int, const std::string&, aurostd::xoption&, const aurostd::xvector<int>&, bool, _xvasp&, _aflags&, _kflags&, _vflags&, std::ofstream&, std::ostream& oss = std::cout); // ME20181107
  bool runRelaxationsAPL_VASP(int, const std::string&, _xvasp&, _aflags&, _kflags&, _vflags&, std::ofstream&); // ME20200427
  void VASP_RunPhonons_APL(_xvasp& xvasp, std::string AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  void RunPhonons_APL(_xinput& xinput, std::string AflowIn, _aflags& aflags, _kflags& kflags, _xflags& xflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // now it's general
  // ----------------------------------------------------------------------------
  // aflow_agl_debye.cpp
  uint relaxStructureAGL_VASP(const std::string& AflowIn, _xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMessage); // CT20200501
  void VASP_RunPhonons_AGL(_xvasp& xvasp, std::string AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMESSAGE);
  void VASP_RunPhonons_AGL_postprocess(const std::string& directory_LIB, std::string& AflowInName, std::string& FileLockName); // CT20200624
  // ----------------------------------------------------------------------------
  // aflow_ael_elasticity.cpp
  uint relaxStructureAEL_VASP(const std::string& AflowIn, _xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMessage); // CT20200501
  void VASP_RunPhonons_AEL(_xvasp& xvasp, std::string AflowIn, _aflags& aflags, _kflags& kflags, _vflags& vflags, std::ofstream& FileMESSAGE);
  void VASP_RunPhonons_AEL_postprocess(const std::string& directory_LIB, std::string& AflowInName, std::string& FileLockName); // CT20200624
} // namespace KBIN

// --------------------------------------------------------------------------------------------------------------------------------------------------------
// aconvasp_aflow.cpp

xstructure PutInCell(const xstructure& a); // Bring all atoms in the cell (will be moved to external function)
xstructure PutInCompact(const xstructure& a); // Bring all atoms in a compact shape (will be moved to external function)
xstructure GetPrim(const xstructure& a);
bool IsTranslationFVector(const xstructure& a, const aurostd::xvector<double>& ftvec);
bool IsTranslationCVector(const xstructure& a, const aurostd::xvector<double>& ctvec);
aurostd::xvector<double> GetMom1(const xstructure& a); // get moment_1 position of the atoms
xstructure SetMom1(const xstructure& a, const aurostd::xvector<double>& mom1_in); // set moment_1 position of atoms
aurostd::xvector<double> AtomCDisp(const _atom& at1, const _atom& at2);
double AtomDist(const xstructure& str, const _atom& atom1, const _atom& atom2); // with structure
double AtomDist(const _atom& at1, const _atom& at2); // without structure
aurostd::xvector<double> GetCDispFromOrigin(const _atom& atom);
double GetDistFromOrigin(const _atom& atom);
_atom ConvertAtomToLat(const _atom& in_at, const aurostd::xmatrix<double>& lattice);
double GetXrayScattFactor(const std::string& name, double lambda = XRAY_RADIATION_COPPER_Kalpha, bool clean = true); // CO20190322
aurostd::xmatrix<double> RecipLat(const aurostd::xmatrix<double>& lat);
double Normal(const double& x, const double& mu, const double& sigma);
xstructure SetLat(const xstructure& a, const aurostd::xmatrix<double>& in_lat);
aurostd::xmatrix<double> GetLat(const xstructure& a);

namespace pflow {
  double GetVol(const aurostd::xmatrix<double>& lat);
  double GetVol(const aurostd::matrix<double>& lat); // CO20200404 pflow::matrix()->aurostd::matrix()
  double GetSignedVol(const aurostd::xmatrix<double>& lat);
  double GetSignedVol(const aurostd::matrix<double>& lat); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::xmatrix<double> RecipLat(const aurostd::xmatrix<double>& lat);
  aurostd::matrix<double> RecipLat(const aurostd::matrix<double>& lat); // CO20200404 pflow::matrix()->aurostd::matrix()
  _atom SetCpos(const _atom& a, const std::vector<double>& in_cpos);
  _atom SetFpos(const _atom& a, const std::vector<double>& in_fpos);
  std::vector<double> vecF2C(const aurostd::matrix<double>& lat, const std::vector<double>& vf); // CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<double> vecC2F(const aurostd::matrix<double>& lat, const std::vector<double>& vc); // CO20200404 pflow::matrix()->aurostd::matrix()
  _atom SetName(const _atom& a, const std::string& in_name);
  _atom SetType(const _atom& a, const int in_type);
  _atom SetNum(const _atom& a, const int in_num);
  // [RF20200415 - duplicate from xatom]std::vector<int> GetTypes(const xstructure& a);
  // [RF20200415 - duplicate from xatom]std::vector<std::string> GetNames(const xstructure& a);
  // [RF20200415 - duplicate from xatom]std::vector<std::string> GetCleanNames(const xstructure& a);
  // [RF20200415 - duplicate from xatom]std::vector<double> GetSpins(const xstructure& a);
  aurostd::matrix<double> GetFpos(const xstructure& str); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> GetCpos(const xstructure& str); // CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure SetLat(const xstructure& a, const aurostd::matrix<double>& in_lat); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> GetLat(const xstructure& a); // CO20200404 pflow::matrix()->aurostd::matrix()
  double GetScale(const xstructure& a);
  aurostd::matrix<double> GetScaledLat(const xstructure& a); // CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure AddAllAtomPos(const xstructure& a, const aurostd::matrix<double>& in_pos, const int in_coord_flag); // CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure SetAllAtomPos(const xstructure& a, const aurostd::matrix<double>& in_pos, const int in_coord_flag); // CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure SetAllAtomNames(const xstructure& a, const std::vector<std::string>& in_names);
  xstructure SetNamesWereGiven(const xstructure& a, const std::vector<int>& in_names_were_given);
  xstructure SetOrigin(const xstructure& a, const std::vector<double>& in_origin);
  xstructure SetOrigin(const xstructure& a, const aurostd::xvector<double>& in_origin);
  bool VVequal(const std::vector<double>& a, const std::vector<double>& b);
  bool VVequal(const std::vector<int>& a, const std::vector<int>& b);
  bool VVequal(const std::deque<double>& a, const std::deque<double>& b);
  bool VVequal(const std::deque<int>& a, const std::deque<int>& b);
  std::vector<double> SmoothFunc(const std::vector<double>& func, const double& sigma);
  void VVset(aurostd::matrix<double>& mat, const double& value); // CO20200404 pflow::matrix()->aurostd::matrix()
  void VVset(std::vector<std::vector<int>>& mat, const int& value);
  double norm(const std::vector<double>& v);
  double getcos(const std::vector<double>& a, const std::vector<double>& b);
  //  vector<double> Getabc_angles(const aurostd::matrix<double>& lat);   // confuses namespace  //CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<double> Sort_abc_angles(const std::vector<double>& abc_angles);
  void Vout(const std::vector<double>& a, std::ostream& out);
  void Vout(const std::vector<int>& a, std::ostream& out);
  void Vout(const std::vector<std::string>& a, std::ostream& out);
  void Mout(const aurostd::matrix<double>& m, std::ostream& out); // CO20200404 pflow::matrix()->aurostd::matrix()
  void Mout(const std::vector<std::vector<double>>& m, std::ostream& out);
  std::vector<double> SVprod(const double& s, const std::vector<double>& b);
  std::vector<int> SVprod(const int& s, const std::vector<int>& b);
  std::vector<double> VVsum(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> VVsum(const std::vector<double>& a, const std::vector<int>& b);
  std::vector<double> VVdiff(const std::vector<double>& a, const std::vector<double>& b);
  double VVprod(const std::vector<double>& a, const std::vector<double>& b);
  double VVprod(const std::vector<double>& a, const std::vector<int>& b);
  aurostd::matrix<double> MMmult(const aurostd::matrix<double>& a, const aurostd::matrix<double>& b); // CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<double> MVmult(const aurostd::matrix<double>& A, const std::vector<double>& v); // CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<double> VMmult(const std::vector<double>& v, const aurostd::matrix<double>& A); // CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<double> VMmult(const std::vector<int>& v, const aurostd::matrix<double>& A); // CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<double> VVcross(const std::vector<double>& a, const std::vector<double>& b);
  double VVdot(const std::vector<double>& a, const std::vector<double>& b);
  int GetNumAtoms(const xstructure& a);
  void SetSpline(const std::vector<double>& x, const std::vector<double>& y, const double& yp1, const double& ypn, std::vector<double>& y2);
  void GetSplineInt(const std::vector<double>& xa, const std::vector<double>& ya, std::vector<double>& y2a, const double& x, double& y);
  void PrintSpline(const std::vector<double>& x, const std::vector<double>& y, const int& npts, std::ostream& outf);
} // namespace pflow

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// aflow_xelement.h stuff
namespace xelement {
  class xelement { // simple class.. nothing fancy
  public:
    // constructor destructor                                       // constructor/destructor
    xelement(); // default, just allocate
    xelement(uint Z, int oxidation_state = AUROSTD_MAX_INT); // look at it by Z
    xelement(const std::string&, int oxidation_state = AUROSTD_MAX_INT); // look at it by symbol or name  //CO20200520
    xelement(const xelement& b); // CO20210201 copy
    ~xelement(); // kill everything
    const xelement& operator=(const xelement& b); // copy
    void clear();
    static uint isElement(const std::string& element); // CO20201220 //SD20220223 - made static
    void loadDefaultUnits(); // CO20201111
    void populate(const std::string& element, int oxidation_state = AUROSTD_MAX_INT); // CO20200520
    void populate(uint ZZ, int oxidation_state = AUROSTD_MAX_INT); // CO20200520
    [[nodiscard]] std::string getPropertyStringVector(const std::string& property, const std::string& delim = ",", uint ncols = AUROSTD_MAX_UINT) const; // CO20201111
    [[nodiscard]] std::string getPropertyString(const std::string& property, const std::string& delim = ",", uint ncols = AUROSTD_MAX_UINT) const; // CO20201111
    [[nodiscard]] double getPropertyDouble(const std::string& property, int index = AUROSTD_MAX_INT) const;
    [[nodiscard]] const aurostd::xvector<double>& getPropertyXVectorDouble(const std::string& property) const;
    [[nodiscard]] const std::vector<double>& getPropertyVectorDouble(const std::string& property) const;
    [[nodiscard]] std::string getType(const std::string& property) const; // CO20201111
    [[nodiscard]] std::string getUnits(const std::string& property) const; // CO20201111
    void convertUnits(const std::string& property = "ALL", const std::string& units_new = "SI"); // CO20201111

    // content                                             // content
    bool verbose;

    // [AFLOW]START=DECLARATION
    uint Z; // Z
    std::string symbol; // http://periodictable.com      //DU20190517   // DONE SC20190524
    std::string name; // http://periodictable.com      //DU20190517   // DONE SC20190524
    uint period; // http://periodictable.com      //DU20190517
    uint group; // http://periodictable.com      //DU20190517
    std::string series; // http://periodictable.com For Nh,Fl,Mc,Lv,Ts Value is a guess based on periodic table trend.      //DU20190517
    std::string block; // http://periodictable.com      //DU20190517
    //
    double mass; // (kg)     // DONE SC20190524
    double volume_molar; // (m^3/mol) http://periodictable.com      //DU20190517
    double volume; // atomic volume in A^3 from the FCC vasp table and/or successive calculations // DONE SC20190524
    double area_molar_Miedema; // (V_m^{2/3} in (cm^2)) (molar volume)^{2/3} surface area Miedema Rule Table 1a Physica 100B (1980) 1-28 10.1016/0378-4363(80)90054-6
    // for lanthines from J.A. Alonso and N.H. March. Electrons in Metals and Alloys, Academic Press, London (1989) (except La)
    double valence_std; // http://en.wikipedia.org/wiki/Valence_(chemistry) standard: number electrons minus closed shell at leff (noble gas)
    double valence_iupac; // http://en.wikipedia.org/wiki/Valence_(chemistry) IUPAC Maximum number of univalent atoms that may combine with an atom of the element under consideration, or with a fragment, or for which an atom of this element can be substituted.
    double valence_PT; //           http://periodictable.com      //DU20190517
    double valence_s; // number of valence s electrons (http://periodictable.com) //CO20201111
    double valence_p; // number of valence p electrons (http://periodictable.com) //CO20201111
    double valence_d; // number of valence f electrons (http://periodictable.com) //CO20201111
    double valence_f; // number of valence f electrons (http://periodictable.com) //CO20201111
    double density_PT; // (g/cm^3)  http://periodictable.com      //DU20190517
    std::string crystal; // Ashcroft-Mermin
    std::string crystal_structure_PT; // http://periodictable.com      //DU20190517
    std::string spacegroup; // http://periodictable.com      //DU20190517
    uint spacegroup_number; // http://periodictable.com      //DU20190517
    double variance_parameter_mass; // Pearson mass deviation coefficient: the square deviation of the isotope masses (weighted by occurrence): 10.1103/PhysRevB.27.858 (isotope corrections), 10.1351/PAC-REP-10-06-02 (isotope distributions) //ME20181020
    aurostd::xvector<double> lattice_constants; // (pm) http://periodictable.com      //DU20190517
    aurostd::xvector<double> lattice_angles; // (rad) http://periodictable.com      //DU20190517
    std::string phase; //      http://periodictable.com      //DU20190517
    double radius_Saxena; // Saxena (nm)
    double radius_PT; // (pm)       http://periodictable.com      //DU20190517
    double radius_covalent_PT; // (pm)       http://periodictable.com      //DU20190517
    double radius_covalent; // (Angstrom) Dalton Trans. 2836, 2832-2838 (2008) //DX+CO20170904
    double radius_VanDerWaals_PT; // (pm)       http://periodictable.com      //DU20190517
    double radii_Ghosh08; // (Angstrom) Journal of Molecular Structure: THEOCHEM 865, 60–67 (2008)      //DU20190517
    double radii_Slatter; // (Angstrom) J. of Chem. Phys. 41, 3199 (1964)      //DU20190517
    double radii_Pyykko; // (Angstrom) single bond covalent radii  Chem. Eur. J. 15, 186-197 (2009)      //DU20190517
    //
    double conductivity_electrical; // (S/m)  http://periodictable.com  Value given for graphite. Diamond electrical conductivity is approximately 0.001.      //DU20190517
    double electronegativity_Pauling; // Saxena
    double hardness_chemical_Ghosh; // (eV) Int. J. Quantum Chem 110, 1206-1213 (2010) Table III       //DU20190517
    double electronegativity_Pearson; // (eV) Inorg. Chem., 27(4), 734–740 (1988)      //DU20190517
    double electronegativity_Ghosh; // (eV) Journal of Theoretical and Computational Chemistry, 4, 21-33 (2005)      //DU20190517

    // RF+SK20200410 START
    //  Allen electronegativities were chosen for CCE since the IUPAC definition of oxidation states seems to use Allen electronegativities and since they also gave the best results
    //  https://en.wikipedia.org/wiki/Oxidation_state#Determination
    //  since there were no Allen electronegativities available for f-elements besides Lu but these elements are usually very similar,
    //  the Lu electronegativity was also used for the other f-elements listed (e.g. La)
    //  this is confirmed by the Allred and Rochow electronegativities that are all very similar for all lanthanides
    double electronegativity_Allen; // https://pubs.acs.org/doi/abs/10.1021/ja00207a003; https://pubs.acs.org/doi/10.1021/ja992866e; https://pubs.acs.org/doi/10.1021/ja9928677
    // preferred and all oxidation states of the elements according to the periodic table of the elements from Wiley-VCH, 5th edition (2012) with some modifications (e. g. for Cr, Cu, Fe, Ti)
    std::vector<double> oxidation_states_preferred;
    std::vector<double> oxidation_states;
    // RF+SK20200410 END

    double electron_affinity_PT; // (kJ/mol)  http://periodictable.com       //DU20190517
    std::vector<double> energies_ionization; // (kJ/mol) http://periodictable.com //CO20201111
    double work_function_Miedema; // (V)        (phi^{\star} empirically-adjusted work function   Miedema Rule Table 1a Physica 100B 1-28 (1980) 10.1016/0378-4363(80)90054-6
    double density_line_electron_WS_Miedema; // (d.u.)^1/3 n_{ws}^{1/3} (averaged electron density at the boundary of the Wigner-Seitz cell)^{1/3}  Miedema Rule Table 1a Physica 100B 1-28 (1980) 10.1016/0378-4363(80)90054-6
    double energy_surface_0K_Miedema; // (mJ/m^2)   \gamma_s^0 surface energy at T=0   Miedema Rule Table 1a Physica 100B 1-28 (1980) 10.1016/0378-4363(80)90054-6
    double chemical_scale_Pettifor; // Chemical Scale Pettifor Solid State Communications 51 31-34 (1984) //updated with D.G. Pettifor 1986 J. Phys. C: Solid State Phys. 19 285  10.1088/0022-3719/19/3/002 //CO20201111
    uint Mendeleev_number; // D.G. Pettifor 1986 J. Phys. C: Solid State Phys. 19 285  10.1088/0022-3719/19/3/002 //CO20201111
    //
    double temperature_boiling; // (Celsius), http://periodictable.com C:diamond, P:"YELLOW" Phosphorus, As:sublimates at this T.      //DU20190517
    double temperature_melting; // (Celsius), http://periodictable.com He does not solidify at standard pressure,C: Value given for diamond form, P : Value given for "YELLOW" phosphorus form, S : Value given for monoclinic, beta form, Se: Value given for hexagonal, gray form, Bk: Value given for alpha form.           //DU20190517
    double enthalpy_fusion; // (kJ/mol)   http://periodictable.com primarily, also https://www.webelements.com/     //CO20201111
    double enthalpy_vaporization; // (kJ/mol)   http://periodictable.com primarily, also https://www.webelements.com/     //DU20190517  //CO20201111
    double enthalpy_atomization_WE; // (kJ/mol)   https://www.webelements.com   //CO20201111
    double energy_cohesive; // (kJ/mol)   http://www.knowledgedoor.com/2/elements_handbook/cohesive_energy.html pulled mostly from Kittel pg 50  //CO20201111
    double specific_heat_PT; // (J/(kg K)) http://periodictable.com Gas_Phase:H(H2),He,N(N2),O(O2),F(F2),Ne,Cl(Cl2),Ar,Kr,Tc,Xe,Rn,Ra,Pa -- Liquid_Phase:Br,Hg -- Solid Phase: B(rhombic),C(graphite),S(rhombic),P(phase of P.4),As(alpha),Se(hexagonal),Cd(gamma),Sn(gray),Li,In,Be,Na,Mg,Al,Si,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,Rb,Sr,Y,Zr,Nb,Mo,Ru,Rh,Pd,Ag,Sb,Te,I,Cs,Ba,La,Ce,Pr,Nd,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Tl,Pb,Bi,Ac,Th,U.      //DU20190517
    double critical_pressure; // (Atm)      http://periodictable.com Li,Na,K,Rb: Value estimated based on extrapolation.      //DU20190517
    double critical_temperature_PT; // (K)        http://periodictable.com Li,Na,K,Rb: Value estimated based on extrapolation.      //DU20190517
    double thermal_expansion; // (K^{-1})   http://periodictable.com C:graphite      //DU20190517
    double conductivity_thermal; // (W/(m K))   http://periodictable.com      //DU20190517
    //
    double hardness_mechanical_Brinell; // (MPa)  http://periodictable.com For Ge value is converted from Mohs scale      //DU20190517
    double hardness_mechanical_Mohs; //        http://periodictable.com For C, value given for graphite. Diamond value is 10.0; For Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Lu converted from Vickers scale. //DU20190517
    double hardness_mechanical_Vickers; // (MPa)  http://periodictable.com For Si,Ge,As,Ru,Os converted from Brinell scale.      //DU20190517
    double hardness_chemical_Pearson; // (eV)   Inorg. Chem. 27(4) 734-740 (1988).      //DU20190517
    double hardness_chemical_Putz; // (eV/atom) International Journal of Quantum Chemistry, Vol 106, 361–389 (2006), TABLE-V. 10.1002/qua.20787      //DU20190517
    double hardness_chemical_RB; // (eV)   Robles and Bartolotti, J. Am. Chem. Soc. 106, 3723-3727 (1984).  10.1021/ja00325a003 using Gunnarsson-Lundqvist (GL) for XC functional    //DU20190517
    double modulus_shear; // (GPa)  http://periodictable.com      //DU20190517
    double modulus_Young; // (GPa)  http://periodictable.com      //DU20190517
    double modulus_bulk; // (GPa)  http://periodictable.com      //DU20190517
    double Poisson_ratio_PT; // (--)   http://periodictable.com      //DU20190517
    double modulus_bulk_x_volume_molar_Miedema; // (kJ/mol) B*V_m Miedema Rule Table 1a Physica 100B 1-28 (1980) 10.1016/0378-4363(80)90054-6
    //
    std::string magnetic_type_PT; //           http://periodictable.com  //DU20190517
    double susceptibility_magnetic_mass; // (m^3/K)   http://periodictable.com //DU20190517
    double susceptibility_magnetic_volume; //           http://periodictable.com //DU20190517
    double susceptibility_magnetic_molar; // (m^3/mol) http://periodictable.com //DU20190517
    double temperature_Curie; // (K)       http://periodictable.com   //DU20190517
    //
    double refractive_index; // http://periodictable.com C:diamond      //DU20190517
    std::string color_PT; // http://periodictable.com      //DU20190517
    //
    double HHIP; // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIP: for elemental production, Uncertinities in HHI_P: C,O,F,Cl,Sc,Ga,Rb,Ru,Rh,Cs,Hf,Os,Ir,Tl.      //DU20190517
    double HHIR; // Chem. Mater. 25(15), 2911–2920 (2013) Herfindahl–Hirschman Index (HHI), HHIR: for elemental reserves,   Uncertinities in HHI_R: Be,C,N,O,F,Na,Mg,Al,Si,S,Cl,Ca,Sc,Ga,Ge,As,Rb,Sr,Ru,Rh,Pd,In,Cs,Hf,Os,Ir,Pt,Tl. //DU20190517
    double xray_scatt; // e-/atom //shift+1 // All data collected from the NIST online tables: http://physics.nist.gov/PhysRefData/FFast/html/form.html  //CO20201111 - another good source: https://henke.lbl.gov/optical_constants/asf.html

    // Xray_scatt_vector All data collected from the NIST online tables
    // http://physics.nist.gov/PhysRefData/FFast/html/form.html
    // All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
    // These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).

    // All data collected from the online tables:
    // http://www-cxro.lbl.gov/optical_constants/pert_form.html
    // All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV].

    // [AFLOW]STOP=DECLARATION

    // UNITS
    std::string units_Z;
    std::string units_symbol;
    std::string units_name;
    std::string units_period;
    std::string units_group;
    std::string units_series;
    std::string units_block;
    //
    std::string units_mass;
    std::string units_volume_molar;
    std::string units_volume;
    std::string units_area_molar_Miedema;
    //
    std::string units_valence_std;
    std::string units_valence_iupac;
    std::string units_valence_PT;
    std::string units_valence_s; // CO20201111
    std::string units_valence_p; // CO20201111
    std::string units_valence_d; // CO20201111
    std::string units_valence_f; // CO20201111
    std::string units_density_PT;
    std::string units_crystal;
    std::string units_crystal_structure_PT;
    std::string units_spacegroup;
    std::string units_spacegroup_number;
    std::string units_variance_parameter_mass;
    std::string units_lattice_constants;
    std::string units_lattice_angles;
    std::string units_phase;
    std::string units_radius_Saxena;
    std::string units_radius_PT;
    std::string units_radius_covalent_PT;
    std::string units_radius_covalent;
    std::string units_radius_VanDerWaals_PT;
    std::string units_radii_Ghosh08;
    std::string units_radii_Slatter;
    std::string units_radii_Pyykko;
    //
    std::string units_conductivity_electrical;
    std::string units_electronegativity_Pauling;
    std::string units_hardness_chemical_Ghosh;
    std::string units_electronegativity_Pearson;
    std::string units_electronegativity_Ghosh;
    std::string units_electronegativity_Allen;
    std::string units_oxidation_states;
    std::string units_oxidation_states_preferred;
    std::string units_electron_affinity_PT;
    std::string units_energies_ionization;
    std::string units_work_function_Miedema;
    std::string units_density_line_electron_WS_Miedema;
    std::string units_energy_surface_0K_Miedema;
    //
    std::string units_chemical_scale_Pettifor;
    std::string units_Mendeleev_number; // CO20201111
    //
    std::string units_temperature_boiling;
    std::string units_temperature_melting;
    std::string units_enthalpy_fusion; // CO20201111
    std::string units_enthalpy_vaporization;
    std::string units_enthalpy_atomization_WE; // CO20201111
    std::string units_energy_cohesive; // CO20201111
    std::string units_specific_heat_PT;
    std::string units_critical_pressure;
    std::string units_critical_temperature_PT;
    std::string units_thermal_expansion;
    std::string units_conductivity_thermal;
    //
    std::string units_hardness_mechanical_Brinell;
    std::string units_hardness_mechanical_Mohs;
    std::string units_hardness_mechanical_Vickers;
    std::string units_hardness_chemical_Pearson;
    std::string units_hardness_chemical_Putz;
    std::string units_hardness_chemical_RB;
    std::string units_modulus_shear;
    std::string units_modulus_Young;
    std::string units_modulus_bulk;
    std::string units_Poisson_ratio_PT;
    std::string units_modulus_bulk_x_volume_molar_Miedema;
    //
    std::string units_magnetic_type_PT;
    std::string units_susceptibility_magnetic_mass;
    std::string units_susceptibility_magnetic_volume;
    std::string units_susceptibility_magnetic_molar;
    std::string units_temperature_Curie;
    //
    std::string units_refractive_index;
    std::string units_color_PT;
    //
    std::string units_HHIP;
    std::string units_HHIR;
    std::string units_xray_scatt;

    // operators/functions                                    // operator/functions
    friend std::ostream& operator<<(std::ostream&, const xelement&); // print
    xelement Initialize(uint Z); // function to clean up the name

  private: //
    void free(); // free space
    void copy(const xelement& b); // copy space //CO20200520
  };

  void Initialize();
  std::string symbol2name(const std::string& symbol);
  std::string name2symbol(const std::string& name);
  int symbol2Z(const std::string& symbol);
  std::string Z2symbol(const int& Z);
  std::string Z2name(const int& Z);
  int name2Z(const std::string& name);

} // namespace xelement

extern std::vector<xelement::xelement> velement; // store starting from ONE

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// CO20201111 - START

namespace aflowMachL { // CO20211111
  void insertElementalProperties(const std::vector<std::string>& vproperties, const xelement::xelement& xel, std::vector<std::string>& vitems);
  void insertElementalPropertiesCoordCE(const std::vector<std::string>& vproperties, const xelement::xelement& xel, double M_X_bonds, double natoms_per_fu, std::vector<std::string>& vitems);
  void insertCrystalProperties(const std::string& structure_path, const std::string& anion, const std::vector<std::string>& vheaders, std::vector<std::string>& vitems, const std::string& e_props = _AFLOW_XELEMENT_PROPERTIES_ALL_);
  double getStatistic(const aurostd::xvector<double>& xvec, const std::string& stat);
  void insertElementalCombinations(const std::vector<std::string>& vproperties, std::vector<std::string>& vheaders);
  void insertElementalCombinations(const std::vector<std::string>& vproperties,
                                   const xelement::xelement& xel_cation,
                                   const xelement::xelement& xel_anion,
                                   const aflowlib::_aflowlib_entry& entry,
                                   double M_X_bonds,
                                   double natoms_per_fu_cation,
                                   double natoms_per_fu_anion,
                                   std::vector<std::string>& vheaders,
                                   std::vector<double>& vfeatures,
                                   bool vheaders_only = false,
                                   uint count_vcols = AUROSTD_MAX_UINT);
  void getColumn(const std::vector<std::vector<std::string>>& table, uint icol, std::vector<std::string>& column, bool& isfloat, bool& isinteger, bool include_header = false);
  void delColumn(std::vector<std::vector<std::string>>& table, uint icol);
  void oneHotFeatures(std::vector<std::vector<std::string>>& table, const std::string& features_categories);
  void removeNaN(const aurostd::xvector<double>& xvec, aurostd::xvector<double>& xvec_new);
  void replaceNaN(aurostd::xvector<double>& xvec, double val = 0.0);
  void MinMaxScale(aurostd::xvector<double>& xvec);
  void reduceFeatures(std::vector<std::vector<std::string>>& table, const std::string& yheader, double var_threshold = _VAR_THRESHOLD_STD_, double ycorr_threshold = _Y_CORR_THRESHOLD_STD_, double selfcorr_threshold = _SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(std::vector<std::vector<std::string>>& table,
                      const std::string& yheader,
                      const std::string& header2skip,
                      double var_threshold = _VAR_THRESHOLD_STD_,
                      double ycorr_threshold = _Y_CORR_THRESHOLD_STD_,
                      double selfcorr_threshold = _SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(std::vector<std::vector<std::string>>& table,
                      const std::string& yheader,
                      const std::vector<std::string>& vheaders2skip,
                      double var_threshold = _VAR_THRESHOLD_STD_,
                      double ycorr_threshold = _Y_CORR_THRESHOLD_STD_,
                      double selfcorr_threshold = _SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(std::vector<std::vector<std::string>>& table,
                      const std::string& yheader,
                      uint icol2skip,
                      double var_threshold = _VAR_THRESHOLD_STD_,
                      double ycorr_threshold = _Y_CORR_THRESHOLD_STD_,
                      double selfcorr_threshold = _SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(std::vector<std::vector<std::string>>& table,
                      const std::string& yheader,
                      const std::vector<uint>& vicol2skip,
                      double var_threshold = _VAR_THRESHOLD_STD_,
                      double ycorr_threshold = _Y_CORR_THRESHOLD_STD_,
                      double selfcorr_threshold = _SELF_CORR_THRESHOLD_STD_);
  std::string reduceEProperties(double var_threshold = _VAR_THRESHOLD_STD_, double selfcorr_threshold = _SELF_CORR_THRESHOLD_STD_);
  void writeCoordCECSV();
  void WriteFileIAPCFG(const aurostd::xoption& vpflow); // CO20211111 //SD20221207 - rewritten using EntryLoader and JSON
  void WriteFilesXYZF(const std::string& outfile, const aurostd::JSON::object& jo); // SD20230919
} // namespace aflowMachL
// CO20201111 - END
//  ----------------------------------------------------------------------------
//  ----------------------------------------------------------------------------
//  aflow_xprototype.h stuff by DAVID

namespace xprototype {
  class xprototype { // stuff in aflow_xprototype.cpp
  public:
    // constructor destructor                          // constructor/destructor
    xprototype(); // default, just allocate
    xprototype(const std::string&); // look at it by symbol or name IN ANRL database
    ~xprototype(); // kill everything
    const xprototype& operator=(const xprototype& b); // copy
    void clear();
    void populate(const std::string& prototype);
    //    void populate(uint ZZ);
    // content                                         // content
    bool verbose;
    // label/params info
    std::string catalog; // prototype catalog 'anrl' or 'htqc'
    uint volume; // volume/part of Encyclopedia
    std::string label; // label (e.g., 201 or AB_cF8_225_a_b)
    std::vector<std::string> parameter_list; // list of degrees of freedom (a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,...)
    std::vector<double> parameter_values; // values for degrees of freedom
    std::string parameter_set_id; // parameter std::set enumeration (e.g., 001, 002, 003, etc.)
    std::string weblink; // link to the corresponding CrystalDatabase web page
    std::vector<uint> stoichiometry; // reduced stoichiometry for prototype (e.g., equicompositional ternary=1:1:1)
    // symmetry
    std::string Pearson_symbol; // Pearson symbol
    uint space_group_number; // space group number
    std::string space_group_symbol_H_M; // space group symbol Hermann-Mauguin (std::optional or use AFLOW lookup table)
    std::string space_group_symbol_Hall; // space group symbol Hall (std::optional or use AFLOW lookup table)
    std::string space_group_symbol_Schoenflies; // space group symbol Schoenflies (std::optional or use AFLOW lookup table)
    std::vector<std::vector<std::string>> Wyckoff_letters; // list of Wyckoff letters grouped by species ([[a,b],[c,d,e],[f,g,h,i],...])
    std::vector<std::vector<std::string>> Wyckoff_site_symmetries; // list of Wyckoff site symmetries grouped by species ([mmm],[2mm,m2m],[mm2],...]) (std::optional, I can grab from look-up table)
    std::vector<std::vector<uint>> Wyckoff_multiplicities; // list of Wyckoff multiplicities grouped by species ([48],[24,24],[12,12,12][4,4,4,4],...]) (std::optional, I can grab from look-up table)
    // designations
    std::string prototype_material; // common prototype material, e.g., NaCl
    std::string common_name; // common prototype name, e.g., half-Heusler
    std::string mineral_name; // mineral name, e.g., corundum
    std::string phase; // compound phase designation (alpha, beta, gamma, delta, etc.) (if applicable)
    std::string strukturbericht; // Strukturbericht designation (if applicable)
    std::vector<std::string> similar_materials; // list of similar compounds (if in same order as stoichiometry we can easily decorate prototypes)
    std::vector<std::string> comments; // noteworthy comments (included in ANRL document and webpage)
    std::string title; // title (for ANRL document/webpage)
    // operators/functions                                    // operator/functions
    friend std::ostream& operator<<(std::ostream&, const xprototype&); // print
    xprototype Iinitialize(uint Z); // function to clean up the name
  private: //
    void free(); // free space
    void copy(const xprototype& b); // copy space
  };
} // namespace xprototype

#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
