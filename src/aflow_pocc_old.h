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

#include "aflow.h"
#include "AUROSTD/aurostd.h"
#include "aflow_pflow.h"

// ***************************************************************************
using std::setfill;
//using aurostd::file2stringstream;
//using aurostd::utype2string;
//using aurostd::PaddedPRE;

//pocc definitions
const std::string _ELEMENTS_STRING_="Ac Ag Al Am Ar As At Au B Ba Be Bi Bk Br C Ca Cd Ce Cf Cl Cm Co Cr Cs Cu Dy Er Es Eu F Fe Fm Fr Ga Gd Ge H He Hf Hg Ho I In Ir K Kr La Li Lr Lu Md Mg Mn Mo N Na Nb Nd Ne Ni No Np O Os P Pa Pb Pd Pm Po Pr Pt Pu Ra Rb Re Rh Rn Ru S Sb Sc Se Si Sm Sn Sr Ta Tb Tc Te Th Ti Tl Tm U V W Xe Y Yb Zn Zr";

// ***************************************************************************
// aflow_contrib_kesong_multienum.h
// ***************************************************************************
namespace pocc {
  void POSCAR2ENUM(istream& input);
  string POSCAR2ENUM(xstructure &a);
  bool MultienumPrintAllXstr(istream& input);
  bool MultienumPrintSortedXstr(istream& input);
  void POSCAR2ENUM(xstructure &a, stringstream &oss, ofstream &FileMESSAGE, _aflags &aflags);
  vector<xstructure> MultienumGenerateXstr(xstructure& xstr, ofstream &FileMESSAGE, _aflags &aflags);
}

// ***************************************************************************
// aflow_contrib_kesong_pocc_basic.h
// ***************************************************************************
//Subroutine for partial occupation
const string name_vacancy="ZZ";
const double epsilon_etot = 1E-6;  //tollerance to remove equivalent structures

void string_replaceAll(std::string& str, const std::string& from, const std::string& to);
bool is_number(const std::string& s);
int gcd(int a, int b);
int lcm (int a, int b);
int CalculateLcmofVector(vector<int> Num);

struct str_num_data {
  double decimal;
  int fac;
  int denom;
};
str_num_data OptimizePoccValue(double dvalue, double tol);
void OptimizeXstr(xstructure &a, ofstream &FileMESSAGE, _aflags &aflags);
void UpdateXstr_comp_each_type(xstructure &b);
void UpdateXstr(xstructure &xstr_orig);
void UpdateXstr(xstructure &xstr_orig, ofstream &FileMESSAGE, _aflags &aflags);

str_num_data double2str_num_data(double a);
bool sortdenom (str_num_data str_i, str_num_data str_j);
vector<int> double2fraction(double a);
vector<int> Decimal2Fraction(double Num);
vector<int> GetFraction(vector<int> IntList);
vector<vector<int> > NumberNormalisedXstructure(xstructure& xstr);
vector<vector<int> > NormalisedNumberXstructure(xstructure& xstr);
vector<vector<int> > CalculateLableXstructure(xstructure& xstr);
vector<vector<int> > CalculatecRange(xstructure& xstr);
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
void PrintGulpEnergies(vector<xstructure>& vxstr);
bool sortenergy (strno_energy str_i, strno_energy str_j);
namespace pocc {
  string POSCAR2GulpInput(xstructure& xstr);
  string POSCAR2GulpInput(xstructure& xstr, vector<string> AtomSpecies);
  void POSCAR2GULP(istream& input);
  vector<double> CalculateEnergyUsingGulp(vector<xstructure>& vxstr);
  double CalculateEnergyUsingGulp(xstructure& xstr);
  double CalculateEnergyUsingGulp(xstructure & xstr, vector<string> AtomSpecies);
  vector<xstructure> SortGroupXstrUFFEnergy(vector<xstructure> groupxstr);
  vector<xstructure> SortGroupXstr(vector<xstructure> groupxstr, vector<string> AtomSpecies);
  // [OBSOLETE]  bool MULTIENUM(vector<string> argv,istream& input);
  bool MultienumPrintSortedXstr(istream& input);
}

// ***************************************************************************
// aflow_contrib_kesong_hnfcell.h
// ***************************************************************************
struct atom_number_name{
  int number;
  string name;
};
struct xstr_atom_number {
  int number;
  vector<int> vec_atom_number;
};
void combine(vector<int> &range, vector<int> &cur, vector<vector<int> > &final_result, int start, int depth);
void combine(vector<int> &range, vector<vector<int> > &final_result, int n);
vector<xmatrix<double> > CalculateHNF(xstructure a, int n);
xstructure XstrSubstitute(xstructure &a, int n, string b);
bool sort_function_struct_num_name (atom_number_name i, atom_number_name j);
vector<atom_number_name> SortMax2Min(vector<atom_number_name>& a);
bool myfunction_int (int i,int j) ;
vector<int> SortMax2Min(vector<int>& a);
bool myfunction_double (double i,double j) ;
vector<double> SortMax2Min(vector<double>& a);
xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, string b);
xstructure XstrSubstitute(xstructure &a, vector<int> vec_n, vector<string> b);
xstructure NormalizeXstructure(xstructure a);
vector<vector<int> > CalculateXstrNumberSupercell(xstructure a, int n);
int hnf_double2int(double a);
bool CheckDoubleOccupied(xstructure xstr_orig);
bool CheckTripleOccupied(xstructure xstr_orig);
bool CheckFourfoldOccupied(xstructure xstr_orig);
bool CheckFivefoldOccupied(xstructure xstr_orig);
void CombineAll(const vector<vector<vector<int> > > &allVecs, size_t vecIndex, vector<vector<int> > &intSoFar, vector<vector<vector<int> > > &results);
vector<string> CalculateSecondAtomicNameSupercell(xstructure xstr_orig, int n);
vector<vector<int> > GenerateSupercellAtomNumber(xstructure xstr_orig, int n);
vector<xstructure> AssignPartialOccupation(xstructure xstr_supercell, xstructure xstr_orig, int n);
int CalculateSupercellSize(xstructure xstr_orig);
xstructure AssignAlphabeticNameXstr(xstructure& xstr);
xstructure CleanNameXstr(xstructure& xstr);
xstructure AssignNameXstr(xstructure& xstr, vector<string> names);
vector<xstructure> AssignNameXstr(vector<xstructure>& vxstr, vector<string> names);
namespace pocc {
  void HNFCELL(istream& input);
}
vector<xstructure> Partial2Supercell(xstructure xstr_ori);
vector<xstructure> CalculateInitialSupercell(xstructure xstr, int n, ofstream &FileMESSAGE, _aflags &aflags);
int InitializeXstr(xstructure &xstr, vector<string> vxstr_species_ori, ofstream &FileMESSAGE, _aflags &aflags);
vector<vector<xstructure> > Partial2Xstr_Fivefold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags);
vector<vector<xstructure> > Partial2Xstr_Fourfold_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags);
vector<vector<xstructure> > Partial2Xstr_Triple_Occupied(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags);
vector<vector<xstructure> > Partial2Xstr(xstructure xstr, int nHNF, ofstream &FileMESSAGE, _aflags &aflags);
vector<xstructure> CalculatePrimitiveCell(vector<xstructure> &vec_xstr);
bool comparison_atom_fpos(_atom atom1, _atom atom2);
namespace pflow {
  void Sort_atom_fpos(xstructure &xstr);
}
bool comparison_atom_cpos(_atom atom1, _atom atom2);
namespace pflow {
  void Sort_atom_cpos(xstructure &xstr);
}
bool LatticeCompare(xstructure xstr1, xstructure xstr2);
bool CoordCompare(xstructure xstr1, xstructure xstr2);
vector<xstr_atom_number> CalculateXstrNumEachType(xstructure xstr);
vector<xstr_energy> xstructure2xstr_energy(vector<xstructure> vec_xstr);
vector<xstructure> RemoveEquivalentXstr(vector<xstructure> vec_xstr, ofstream &FileMESSAGE, _aflags &aflags);
namespace pocc {
  bool DIFF(xstructure xstr1, xstructure xstr2);
  bool DIFF_LEVEL1(xstructure xstr1, xstructure xstr2);
  bool DIFF_LEVEL4(xstructure xstr1, xstructure xstr2);
  vector<xstructure> GenerateRotatedXstr(xstructure xstr);
  void DIFF(string options);
  bool CompareLattice(xstructure xstr1, xstructure xstr2);
}

// ***************************************************************************
// aflow_contrib_kesong_ipocc.h
// ***************************************************************************
namespace pflow {
  void POCC_INPUT(void);
}
bool POCC_GENERATE_INPUT(ofstream &FileMESSAGE,_aflags &aflags);

// ***************************************************************************
// aflow_contrib_kesong_std.h
// ***************************************************************************

// ***************************************************************************
// aflow_contrib_kesong_poccdos.h
// ***************************************************************************
void ExtracAllPOSCARSFromAflowin(vector<xstructure>& vxstr, const string& str_aflowin);
void GetDegeneracyFromVPOSCAR(const vector<xstructure>& vxstr, vector<int>& vDE);

namespace pocc {  
  void POCC_SetOptions(string options, string& directory, double& T, double& DOS_Emin, double& DOS_Emax, double& DOSSCALE);
  void POCC_DOS(ostream &oss,string options);
  bool POCC_Already_Calculated_Input(const string& str_AflowIn,const string& dir=".");
  void POCC_CalculateDistribution(vector<double>& vdelta_toten, vector<double>& vprob, const double& T, const string& NameDist, const vector<int>& vDE);
  void POCC_GENERATE_DOSDATA(const string& str_dir, const double& T, vector<vector<double> >& TDOS_ONLY, vector<vector<double> >& PDOS_ONLY, vector<double>& POCC_Efermi, double& POCC_mag, vector<double>& vprob);
  void POCC_COMBINE_TDOS_PDOS_ONEDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& PDOS, vector<vector<double> >& DOS, vector<vector<double> >& DOS_IDOS);  //CO20190808 - one without IDOS and one with
  string POCC_GENERATE_GNUPLOTSCRIPT(vector<vector<double> >& DOS, const string& SystemName, const string& dosdatafile, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE, const double& DOSMAX);
  void POCC_GENERATE_OUTPUT(const string& str_dir, const double& T, const double& DOS_Emin, double& DOS_Emax, const double& DOSSCALE);
  void POCC_BANDGAP(string options);
  void POCC_MAG(string options);
  vector<vector<double> > SpinFlipDOS(const vector<vector<double> >& vva);
  vector<vector<double> > SpinSplitDOS(const vector<vector<double> >& vva);
  vector<vector<double> > POCC_Formalise(const bool& SPIN_FLAG, const double& weight, const double& mag, const vector<vector<double> >& TDOS);
}

// ***************************************************************************
// aflow_contrib_kesong.h
// ***************************************************************************

int StringCrop(string directory, vector<string>& vstring);

// ***************************************************************************

#endif

// ***************************************************************************
// *                                                                         *
// *              AFlow KESONG YANG - Duke University 2010-2011              *
// *                                                                         *
// ***************************************************************************

