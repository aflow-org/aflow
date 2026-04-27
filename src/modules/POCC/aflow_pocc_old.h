// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2011              *
// *                                                                         *
// ***************************************************************************
// aflow_contrib_kesong.h
// functions written by
// 2010-2010: kesong.yang@gmail.com

#ifndef _AFLOW_CONTRIB_KESONG_H_
#define _AFLOW_CONTRIB_KESONG_H_

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd_xmatrix.h"

#include "flow/aflow_xclasses.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

// ***************************************************************************

// pocc definitions
const std::string _ELEMENTS_STRING_
    = "Ac Ag Al Am Ar As At Au B Ba Be Bi Bk Br C Ca Cd Ce Cf Cl Cm Co Cr Cs Cu Dy Er Es Eu F Fe Fm Fr Ga Gd Ge H He Hf Hg Ho I In Ir K Kr La Li Lr Lu Md Mg Mn Mo N Na Nb Nd Ne Ni No Np O Os P Pa Pb Pd Pm Po "
      "Pr "
      "Pt Pu Ra Rb Re Rh Rn Ru S Sb Sc Se Si Sm Sn Sr Ta Tb Tc Te Th Ti Tl Tm U V W Xe Y Yb Zn Zr";

// ***************************************************************************
// aflow_contrib_kesong_multienum.h
// ***************************************************************************
namespace pocc {
  void POSCAR2ENUM(std::istream& input);
  std::string POSCAR2ENUM(xstructure& a);
  bool MultienumPrintAllXstr(std::istream& input);
  bool MultienumPrintSortedXstr(std::istream& input);
  void POSCAR2ENUM(xstructure& a, std::stringstream& oss, std::ofstream& FileMESSAGE, _aflags& aflags);
  std::vector<xstructure> MultienumGenerateXstr(xstructure& xstr, std::ofstream& FileMESSAGE, _aflags& aflags);
} // namespace pocc

// ***************************************************************************
// aflow_contrib_kesong_pocc_basic.h
// ***************************************************************************
// Subroutine for partial occupation
const std::string name_vacancy = "ZZ";
const double epsilon_etot = 1E-6;  // tollerance to remove equivalent structures

void string_replaceAll(std::string& str, const std::string& from, const std::string& to);
bool is_number(const std::string& s);
int gcd(int a, int b);
int lcm(int a, int b);
int CalculateLcmofVector(std::vector<int> Num);

struct str_num_data {
  double decimal;
  int fac;
  int denom;
};
str_num_data OptimizePoccValue(double dvalue, double tol);
void OptimizeXstr(xstructure& a, std::ofstream& FileMESSAGE, _aflags& aflags);
void UpdateXstr_comp_each_type(xstructure& b);
void UpdateXstr(xstructure& xstr_orig);
void UpdateXstr(xstructure& xstr_orig, std::ofstream& FileMESSAGE, _aflags& aflags);

str_num_data double2str_num_data(double a);
bool sortdenom(str_num_data str_i, str_num_data str_j);
std::vector<int> double2fraction(double a);
std::vector<int> Decimal2Fraction(double Num);
std::vector<int> GetFraction(std::vector<int> IntList);
std::vector<std::vector<int>> NumberNormalisedXstructure(xstructure& xstr);
std::vector<std::vector<int>> NormalisedNumberXstructure(xstructure& xstr);
std::vector<std::vector<int>> CalculateLableXstructure(xstructure& xstr);
std::vector<std::vector<int>> CalculatecRange(xstructure& xstr);
double CalculatePartialValueOfVacancy(xstructure& xstr, unsigned int k);
double CalculateDistanceXstructure(xstructure& xstr, int i, int j);
bool CheckVacancyOnOnesite(xstructure& xstr, unsigned int k);
bool CheckVacancy(xstructure& xstr_in);
bool CheckPartialOccupation(xstructure& xstr);
bool CheckOneSite(xstructure& xstr, int i, int j);
bool CheckMultiOccupied(xstructure& xstr);
xstructure CleanVacancy(xstructure& xstr);

struct strno_energy {
  double energy;
  int number;
};

struct xstr_energy {
  xstructure xstr;
  double energy;
};
void PrintGulpEnergies(std::vector<xstructure>& vxstr);
bool sortenergy(strno_energy str_i, strno_energy str_j);
namespace pocc {
  std::string POSCAR2GulpInput(xstructure& xstr);
  std::string POSCAR2GulpInput(xstructure& xstr, std::vector<std::string> AtomSpecies);
  void POSCAR2GULP(std::istream& input);
  std::vector<double> CalculateEnergyUsingGulp(std::vector<xstructure>& vxstr);
  double CalculateEnergyUsingGulp(xstructure& xstr);
  double CalculateEnergyUsingGulp(xstructure& xstr, std::vector<std::string> AtomSpecies);
  std::vector<xstructure> SortGroupXstrUFFEnergy(std::vector<xstructure> groupxstr);
  std::vector<xstructure> SortGroupXstr(std::vector<xstructure> groupxstr, std::vector<std::string> AtomSpecies);
  bool MultienumPrintSortedXstr(std::istream& input);
} // namespace pocc

// ***************************************************************************
// aflow_contrib_kesong_hnfcell.h
// ***************************************************************************
struct atom_number_name {
  int number;
  std::string name;
};
struct xstr_atom_number {
  int number;
  std::vector<int> vec_atom_number;
};
void combine(std::vector<int>& range, std::vector<int>& cur, std::vector<std::vector<int>>& final_result, int start, int depth);
void combine(std::vector<int>& range, std::vector<std::vector<int>>& final_result, int n);
std::vector<aurostd::xmatrix<double>> CalculateHNF(xstructure a, int n);
xstructure XstrSubstitute(xstructure& a, int n, std::string b);
bool sort_function_struct_num_name(atom_number_name i, atom_number_name j);
std::vector<atom_number_name> SortMax2Min(std::vector<atom_number_name>& a);
bool myfunction_int(int i, int j);
std::vector<int> SortMax2Min(std::vector<int>& a);
bool myfunction_double(double i, double j);
std::vector<double> SortMax2Min(std::vector<double>& a);
xstructure XstrSubstitute(xstructure& a, std::vector<int> vec_n, std::string b);
xstructure XstrSubstitute(xstructure& a, std::vector<int> vec_n, std::vector<std::string> b);
xstructure NormalizeXstructure(xstructure a);
std::vector<std::vector<int>> CalculateXstrNumberSupercell(xstructure a, int n);
int hnf_double2int(double a);
bool CheckDoubleOccupied(xstructure xstr_orig);
bool CheckTripleOccupied(xstructure xstr_orig);
bool CheckFourfoldOccupied(xstructure xstr_orig);
bool CheckFivefoldOccupied(xstructure xstr_orig);
void CombineAll(const std::vector<std::vector<std::vector<int>>>& allVecs, size_t vecIndex, std::vector<std::vector<int>>& intSoFar, std::vector<std::vector<std::vector<int>>>& results);
std::vector<std::string> CalculateSecondAtomicNameSupercell(xstructure xstr_orig, int n);
std::vector<std::vector<int>> GenerateSupercellAtomNumber(xstructure xstr_orig, int n);
std::vector<xstructure> AssignPartialOccupation(xstructure xstr_supercell, xstructure xstr_orig, int n);
int CalculateSupercellSize(xstructure xstr_orig);
xstructure AssignAlphabeticNameXstr(xstructure& xstr);
xstructure CleanNameXstr(xstructure& xstr);
xstructure AssignNameXstr(xstructure& xstr, std::vector<std::string> names);
std::vector<xstructure> AssignNameXstr(std::vector<xstructure>& vxstr, std::vector<std::string> names);
namespace pocc {
  void HNFCELL(std::istream& input);
}
std::vector<xstructure> Partial2Supercell(xstructure xstr_ori);
std::vector<xstructure> CalculateInitialSupercell(xstructure xstr, int n, std::ofstream& FileMESSAGE, _aflags& aflags);
int InitializeXstr(xstructure& xstr, std::vector<std::string> vxstr_species_ori, std::ofstream& FileMESSAGE, _aflags& aflags);
std::vector<std::vector<xstructure>> Partial2Xstr_Fivefold_Occupied(xstructure xstr, int nHNF, std::ofstream& FileMESSAGE, _aflags& aflags);
std::vector<std::vector<xstructure>> Partial2Xstr_Fourfold_Occupied(xstructure xstr, int nHNF, std::ofstream& FileMESSAGE, _aflags& aflags);
std::vector<std::vector<xstructure>> Partial2Xstr_Triple_Occupied(xstructure xstr, int nHNF, std::ofstream& FileMESSAGE, _aflags& aflags);
std::vector<std::vector<xstructure>> Partial2Xstr(xstructure xstr, int nHNF, std::ofstream& FileMESSAGE, _aflags& aflags);
std::vector<xstructure> CalculatePrimitiveCell(std::vector<xstructure>& vec_xstr);
bool comparison_atom_fpos(_atom atom1, _atom atom2);
namespace pflow {
  void Sort_atom_fpos(xstructure& xstr);
}
bool comparison_atom_cpos(_atom atom1, _atom atom2);
namespace pflow {
  void Sort_atom_cpos(xstructure& xstr);
}
bool LatticeCompare(xstructure xstr1, xstructure xstr2);
bool CoordCompare(xstructure xstr1, xstructure xstr2);
std::vector<xstr_atom_number> CalculateXstrNumEachType(xstructure xstr);
std::vector<xstr_energy> xstructure2xstr_energy(std::vector<xstructure> vec_xstr);
std::vector<xstructure> RemoveEquivalentXstr(std::vector<xstructure> vec_xstr, std::ofstream& FileMESSAGE, _aflags& aflags);
namespace pocc {
  bool DIFF(xstructure xstr1, xstructure xstr2);
  bool DIFF_LEVEL1(xstructure xstr1, xstructure xstr2);
  bool DIFF_LEVEL4(xstructure xstr1, xstructure xstr2);
  std::vector<xstructure> GenerateRotatedXstr(xstructure xstr);
  void DIFF(std::string options);
  bool CompareLattice(xstructure xstr1, xstructure xstr2);
} // namespace pocc

// ***************************************************************************
// aflow_contrib_kesong_ipocc.h
// ***************************************************************************
namespace pflow {
  void POCC_INPUT();
}
bool POCC_GENERATE_INPUT(std::ofstream& FileMESSAGE, _aflags& aflags);

// ***************************************************************************
// aflow_contrib_kesong_std.h
// ***************************************************************************

// ***************************************************************************
// aflow_contrib_kesong_poccdos.h
// ***************************************************************************
void ExtracAllPOSCARSFromAflowin(std::vector<xstructure>& vxstr, const std::string& str_aflowin);
void GetDegeneracyFromVPOSCAR(const std::vector<xstructure>& vxstr, std::vector<int>& vDE);

namespace pocc {
  void POCC_SetOptions(std::string options, std::string& directory, double& T, double& DOS_Emin, double& DOS_Emax, double& DOSSCALE);
  void POCC_DOS(std::ostream& oss, std::string options);
  bool POCC_Already_Calculated_Input(const std::string& str_AflowIn, const std::string& dir = ".");
  void POCC_CalculateDistribution(std::vector<double>& vdelta_toten, std::vector<double>& vprob, const double& T, const std::string& NameDist, const std::vector<int>& vDE);
  void POCC_GENERATE_DOSDATA(
      const std::string& str_dir, const double& T, std::vector<std::vector<double>>& TDOS_ONLY, std::vector<std::vector<double>>& PDOS_ONLY, std::vector<double>& POCC_Efermi, double& POCC_mag, std::vector<double>& vprob);
  void POCC_COMBINE_TDOS_PDOS_ONEDOS(const std::vector<std::vector<double>>& TDOS, const std::vector<std::vector<double>>& PDOS, std::vector<std::vector<double>>& DOS, std::vector<std::vector<double>>& DOS_IDOS); // CO20190808 - one without IDOS and one with
  std::string POCC_GENERATE_GNUPLOTSCRIPT(
      std::vector<std::vector<double>>& DOS, const std::string& SystemName, const std::string& dosdatafile, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE, const double& DOSMAX);
  void POCC_GENERATE_OUTPUT(const std::string& str_dir, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE);
  void POCC_BANDGAP(std::string options);
  void POCC_MAG(std::string options);
  std::vector<std::vector<double>> SpinFlipDOS(const std::vector<std::vector<double>>& vva);
  std::vector<std::vector<double>> SpinSplitDOS(const std::vector<std::vector<double>>& vva);
  std::vector<std::vector<double>> POCC_Formalise(const bool& SPIN_FLAG, const double& weight, const double& mag, const std::vector<std::vector<double>>& TDOS);
} // namespace pocc

// ***************************************************************************
// aflow_contrib_kesong.h
// ***************************************************************************

int StringCrop(std::string directory, std::vector<std::string>& vstring);

// ***************************************************************************

#endif

// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2011              *
// *                                                                         *
// ***************************************************************************
