// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Richard H. Taylor
// UPDATED by David Hicks
// d.hicks@duke.edu

#ifndef _AFLOW_SYMMETRY_SPACEGROUP_H_
#define _AFLOW_SYMMETRY_SPACEGROUP_H_
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xvector.h"

#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_xatom.h"

const double Pi_r = 3.141592653589793;

namespace SYM {
  // Symop
  struct symop {
    friend std::ostream& operator<<(std::ostream& output, const symop& a);
    std::string symbol;
    aurostd::xvector<double> direction;
    aurostd::xvector<double> screwglide;  // will also store inversion pnts for rotoinversion
    aurostd::xvector<double> shift;

    void clear();
  };
} // namespace SYM

class SymmetryInformationITC {
public:
  SymmetryInformationITC();
  ~SymmetryInformationITC();
  friend std::ostream& operator<<(std::ostream& oss, const SymmetryInformationITC& SymmetryInformationITC);
  const SymmetryInformationITC& operator=(const SymmetryInformationITC& b);
  SymmetryInformationITC(const SymmetryInformationITC& b);

    // glides
  std::vector<aurostd::xvector<double>> glideplanes;
  std::vector<aurostd::xvector<double>> glideplanes_hex;
  std::vector<aurostd::xvector<double>> glidetrans;
  std::vector<aurostd::xvector<double>> glidetrans_hex;
  std::vector<std::string> glidesymbols;
  std::vector<std::string> glidesymbols_hex;

    // symmetry matrices
  std::vector<int> index_cubic;
  std::vector<int> index_hex;  // To be used with the *_hex vectors
  std::vector<int> index_rhom;
  std::vector<int> index_tetr;
  std::vector<int> index_ortho;
  std::vector<int> index_mono_b;
  std::vector<int> index_mono_c;
  std::vector<int> index_tric;
  std::vector<aurostd::xmatrix<double>> sym_mats;
  std::vector<std::string> symbol;
  std::vector<std::string> dirparam;
  std::map<int, aurostd::xvector<double>> sym_mats_direction;
  std::vector<aurostd::xmatrix<double>> sym_mats_hex;
  std::vector<std::string> symbol_hex;
  std::vector<std::string> dirparam_hex;
  std::map<int, aurostd::xvector<double>> sym_mats_direction_hex;

    // sym_ops
  std::vector<std::string> sym_ops;

    // generators
  std::vector<std::vector<SYM::symop>> generators;
  std::vector<int> sgindex;

  std::vector<std::string> gl_sgs;

    // functions
  bool initsymmats();
  bool initglides();
  bool initsymops();
  bool initgenerators(std::string axis_cell);
  bool initsgs(std::string axis_cell);

private:
  void free();
  void copy(const SymmetryInformationITC& b);
};

// ******************************************************************************
// xstructure
class xstructure;  // forward define xstructure for compilation.

std::vector<int> AllCombination41(int num, int total_num, int index);
std::vector<int> AllCombination42(int num, int total_num, std::vector<int>& str_in);
unsigned long int CombinationNr(int num, int total_num);

namespace SYM {
  // stringdouble (Used for Wyckoff position algebra)
  class stringdouble {
    friend std::ostream& operator<<(std::ostream& output, const stringdouble& a);

  public:
    std::string type;
    double distance;
    aurostd::xvector<double> coord;
    stringdouble() {
      type = "XX";
      distance = 0;
    }
    ~stringdouble() {};
  };

  // wyckoff site
  class wyckoffsite  // Also for wyckoff sites
  {
    // friend bool operator==(const wyckoffsite& a, const wyckoffsite& b);
    friend std::ostream& operator<<(std::ostream& output, const wyckoffsite& a);

  public:
    wyckoffsite() { wyckoffSymbol = " "; };
    ~wyckoffsite() {};
    aurostd::xvector<double> coord;
    std::string type;  // chemical label etc
    std::string wyckoffSymbol;
  };

  // atom class manipulation std::function

  // Screw (Screw Operation)
  class Screw {
    friend aurostd::xvector<double> operator*(Screw& S, const aurostd::xvector<double>& P);
    friend aurostd::xvector<double> operator*(const aurostd::xvector<double>& P, Screw& S);
    friend std::ostream& operator<<(std::ostream& output, const Screw& S);

  private:
    aurostd::xmatrix<double> A;  // operator matrix
    void get_A();
    std::string linestring;          // parametric line
    aurostd::xvector<double> T;          // translation component of screw
    double angle;               // angle of rotation
    aurostd::xvector<double> one_point;  // a point on the axis (in get_A)
    aurostd::xvector<double> direction_vector;

  public:
    Screw() {};
    ~Screw() {};
    void get_screw_direct(aurostd::xvector<double> axis_dir, aurostd::xvector<double> point, double order);
    void get_line_string(std::string str);
    void get_trans_comp(aurostd::xvector<double> trans);
    void get_angle(double theta);
    aurostd::xvector<double> return_direction();
    aurostd::xvector<double> return_point();
    [[nodiscard]] double return_order() const;
  };

  // Glide (Glide Operation)
  class Glide {  // for a plane defined on the plane ax+by+cz=d with translation
    //(t1,t2,t3) the glide operation works as follows:
    // A.P + d*a + |a.P - d|*a + (t1,t2,t3)^T
  private:
    bool HEX;
    bool DIRECT;
    std::string planestring;           // parametric plane
    aurostd::xmatrix<double> A;            // reflecting matrix
    aurostd::xvector<double> a;            // normal std::vector or fixed point (if HEX)
    aurostd::xvector<double> T;            // translation std::vector
    aurostd::xvector<double> plane_point;  // point in plane for get_glide_direct
    double d;                     // d in ax+by+cz+d
    void get_A(aurostd::xvector<double> n);

  public:
    Glide() {
      HEX = false;
      DIRECT = false;
      d = 0.0;
    };
    ~Glide() {};
      // overload operator * for Glide
      //  cartesian
    void get_glide_direct(aurostd::xvector<double> n, aurostd::xvector<double> p);
    void get_glide_direct(aurostd::xvector<double> n, aurostd::xvector<double> p, aurostd::xvector<double> trans);
    aurostd::xvector<double> return_point();      // returns point on plane
    aurostd::xvector<double> return_direction();  // returns normal to plane
    friend aurostd::xvector<double> operator*(Glide& G, const aurostd::xvector<double>& P);
    friend aurostd::xvector<double> operator*(const aurostd::xvector<double>& P, Glide& G);
    friend std::ostream& operator<<(std::ostream& output, const Glide& S);
  };

  // Translation
  class Translation {
  private:
    aurostd::xvector<double> translation_vec;

  public:
    Translation() {};
    ~Translation() {};
    void get_translation(std::string ITCstring);
    friend aurostd::xvector<double> operator*(Translation& T, const aurostd::xvector<double>& P);
    friend aurostd::xvector<double> operator*(const aurostd::xvector<double>& P, Translation& T);
    friend aurostd::xvector<double> operator+(Translation& T, const aurostd::xvector<double>& P);
    friend aurostd::xvector<double> operator+(const aurostd::xvector<double>& P, Translation& T);
  };

  // Inversion
  class Inversion {
    friend aurostd::xvector<double> operator*(Inversion& I, const aurostd::xvector<double>& P);
    friend aurostd::xvector<double> operator*(const aurostd::xvector<double>& P, Inversion& I);

  private:
  public:
    Inversion() {};
    ~Inversion() {};
    // overload operator * for Glide
    aurostd::xvector<double> inversion_point;
    void get_inversion(std::string ITCstring);
  };
} // namespace SYM

// End Classes
//  ******************************************************************************

// ******************************************************************************
// Structure Declarations
// ******************************************************************************

namespace SYM {
  // Symfolder (Stores all symmetry elements/operators)
  struct symfolder {
    std::vector<Screw> twofold_ops;
    std::vector<Screw> rot_ops;
    std::vector<Glide> mirror_ops;
    std::vector<aurostd::xvector<double>> twofold_lattice_vectors;
    std::vector<aurostd::xvector<double>> rot_lattice_vectors;
    std::vector<aurostd::xvector<double>> mirror_lattice_vectors;
    bool commensurate;
    std::string crystalsystem;
    std::string latticesystem;
    std::vector<int> lattice_pgroups;                // DX NEW
    std::vector<aurostd::xmatrix<double>> lattice_sym_mats;  // DX NEW
    std::vector<aurostd::xmatrix<double>> crystal_sym_mats;  // DX NEW
    std::vector<int> insym;
    std::vector<aurostd::xvector<double>> translations;
    std::vector<std::vector<int>> all_atom_maps;
    std::vector<std::vector<int>> all_type_maps;
    std::vector<aurostd::xvector<double>> centeringops;
    std::vector<std::string> pointgroupops;
  };

  ////Symop
  // struct symop {
  //   friend std::ostream& operator<<(std::ostream& output, const symop& a);
  //   std::string symbol;
  //   xvector<double> direction;
  //   xvector<double> screwglide;  //will also store inversion pnts for rotoinversion
  //   xvector<double> shift;
  //
  //     void clear();
  //   };

  // Eqatoms
  struct eqatoms {
    std::string label;
    std::vector<std::vector<double>> atoms;
  };

  // sdouble (String and double)
  struct sdouble {
    char chr;
    double dbl;
    sdouble() {
      chr = '\0';
      dbl = 0;
    }
  };

  // enum_alpha (enumerate alphabet: for Wyckoff letters)
  struct enum_alpha {
    std::string letter;
    int number;
  };

  // End Structures
  //  ******************************************************************************

  // ******************************************************************************
  // Template Declarations
  // ******************************************************************************

  template <class dmmy> bool allsame(std::vector<dmmy> vec) {
    bool all = true;
    for (uint i = 0; i < vec.size(); i++) {
      if (vec[i] != vec[0]) {
        all = false;
      }
    }
    return all;
  }

} // namespace SYM

// End Templates
//  ******************************************************************************

// ******************************************************************************
// Function Declarations
// ******************************************************************************
// MAIN FUNCITONS
namespace SYM {
  void calculateSpaceGroups(std::vector<xstructure>& vxstrs, uint start_index = 0, uint end_index = AUROSTD_MAX_UINT, uint setting = 0); // DX20191230 add setting option
  std::string OrthoDefect(std::istream& cin);
  xstructure SpaceGroup(std::istream& cin);
  void rgcd(std::vector<std::string> num);
  // DX void AtomicEnvironment(std::istream & cin, std::vector<std::string> av);
  std::string ReverseSpaceGroup(std::vector<std::string> num);
  // END MAIN FUNCTIONS

  // FUNCTION THAT TAKES WYCCAR FROM XSTRUCTURE AND PUTS IT IN OSTREAM
  void printWyccar(std::ofstream& FileMESSAGE, xstructure& str, const bool& osswrite, std::ostream& oss);

  // LINEAR ALGEBRA FUNCTIONS
  bool solve_overdetermined_system(std::vector<aurostd::xvector<double>>& LHS, std::vector<double>& RHS, aurostd::xvector<double>& SOL, aurostd::xmatrix<double>& lattice, double& min_dist, double& tol); // DX20190215
  void ReducedRowEchelonForm(aurostd::xmatrix<double>& M, double& tol); // DX20190215
  bool find_solution_UD(aurostd::xmatrix<double> M, aurostd::xvector<double>& SOL, double& tol); // check if underdetermined system has a solution, and get a particular solution. //DX20190215
  bool checkLinearSystem(std::vector<aurostd::xvector<double>>& LHS, std::vector<double>& RHS, aurostd::xmatrix<double>& lattice, double& tol); // DX20190215
  // WYCKOFF FUNCITONS
  std::vector<double> system_solve(std::vector<double> numeric, std::vector<std::vector<sdouble>> variable);
  std::vector<std::vector<double>> wyckoff_solve(std::vector<std::vector<std::vector<sdouble>>> win, std::vector<eqatoms> pin);
  std::vector<std::vector<double>> wyckoff_solve_2(std::vector<std::vector<std::vector<sdouble>>> win, eqatoms pin);
  std::vector<std::vector<double>> wyckoff_sites(std::vector<std::vector<std::vector<sdouble>>> win, std::vector<_atom> atomgroup);

  // SPACE GROUP LIBRARY FUNCTIONS
  bool initsgs(std::string axis_cell);
  bool initsymops();
  bool initsymmats();
  bool initglides();
  bool initgenerators(std::string axis_cell);
  std::vector<int> generatorindex(int spacegroup);
  std::vector<aurostd::xvector<double>> ReturnITCGenShift(int spacegroupnum, std::string axis_cell);

  std::string point_group_library_search(std::vector<std::string> symmetryoperations, std::string crystalsystem, int centeringops);
  std::string crystal_system_library_search(std::vector<std::string> symmetryoperations);
  std::vector<int> spacegroups_CS_LS(std::string crystalsystem, char latticesystem);
  char discern_sym_op(std::string ITCstring);
  int* PointGroup_SpaceGroup(std::string pgroup);
  std::vector<int> PointGroup_SpaceGroup(std::string pgroup, char cl);
  aurostd::xmatrix<double> spacegroup_to_lattice(int spacegroup);

  std::vector<int> get_multiplicities(std::string sg);
  std::vector<std::string> get_symmetry_symbols(std::string sg);

  std::vector<std::string> get_wyckoff_equation(std::string spaceg, int mult); // DX20170830
  std::vector<std::vector<std::vector<std::string>>> get_wyckoff_pos(std::string spaceg, int mult);
  std::vector<std::vector<std::vector<std::string>>> get_wyckoff_pos(std::string spaceg, int mult, bool getcentering);
  std::vector<std::vector<std::string>> get_wyckoff_pos(std::string spaceg, int Wyckoff_multiplicity, std::string Wyckoff_letter); // DX20190129
  std::vector<std::string> get_minimum_enumerated_Wyckoff_letters(std::string spacegroupstring, std::vector<int>& multiplicities, std::vector<std::string> site_symmetries);
  int enumerate_wyckoff_letter(std::string& wyckoff_letter);
  std::vector<int> enumerate_wyckoff_letters(std::vector<std::string>& wyckoff_letters); // DX20180927
  void get_all_wyckoff_for_site_symmetry(std::string spaceg, int mult, std::string site_symmetry, std::vector<std::vector<std::string>>& all_positions);
  void get_Wyckoff_from_letter(uint space_group_number, std::string& space_group_setting, std::string& Wyckoff_letter, uint& Wyckoff_multiplicity, std::string& site_symmetry, std::vector<std::string>& positions); // DX20191029
  void get_Wyckoff_from_letter(std::string& spaceg, std::string& Wyckoff_letter, uint& Wyckoff_multiplicity, std::string& site_symmetry, std::vector<std::string>& positions);
  aurostd::xvector<double> Wyckoff_position_string2xvector(std::string& string_position);
  std::vector<aurostd::xvector<double>> get_possible_origin_shifts(std::string spacegroupstring, int multiplicity, std::string site_symmetry);
  void get_certain_wyckoff_pos(std::string spaceg, int mult, std::string site_symmetry, std::vector<std::string>& site_symmetries, std::vector<std::string>& letters, std::vector<std::string>& positions);
  void getWyckoffAutomorphismSets(const std::string& spaceg,
                                  int mult,
                                  std::string site_symmetry,
                                  std::vector<std::string>& site_symmetries,
                                  std::vector<std::string>& letters,
                                  std::vector<std::vector<std::string>>& all_positions,
                                  bool origin_shifts_only = false); // DX20220828
  void getGeneralWyckoffMultiplicityAndPosition(uint space_group_number, std::string& space_group_setting, int& general_wyckoff_multiplicity, std::vector<std::string>& general_wyckoff_position);
  std::vector<std::string> findGeneralWyckoffPosition(std::string& spacegroupstring, int& general_wyckoff_multiplicity);
  uint numberOccupiedSitesInConventionalCell(const std::vector<wyckoffsite_ITC>& Wyckoff_sites); // DX20200512
  std::vector<uint> numberEachTypeFromWyckoff(const std::vector<wyckoffsite_ITC>& Wyckoff_sites); // DX20200512
  std::vector<std::string> findWyckoffEquations(uint space_group_number, std::string& space_group_setting, std::string& Wyckoff_letter, uint Wyckoff_multiplicity, bool convert2frac = true); // DX 20191029 //DX20200423 - add convert2frac
  std::vector<std::string> findWyckoffEquations(std::string& spacegroupstring, std::string& Wyckoff_letter, uint Wyckoff_multplicity, bool convert2frac = true, bool keep_multiplication_symbol = false); // DX 20190128 //DX20200423 - add convert2frac, keep_multiplication_symbol
  std::string formatWyckoffPosition(const std::vector<sdouble>& sd_coordinate,
                                    bool convert2frac = true,
                                    bool keep_multiplication_symbol = false,
                                    int precision = AUROSTD_DEFAULT_PRECISION); // DX 20190723 //DX20200423 - add convert2frac, keep_multiplication_symbol  //CO20220607 - added precision
  std::string reorderWyckoffPosition(const std::string& orig_position); // DX 20190708
  bool shiftWyckoffPositions(std::deque<std::deque<_atom>>& equivalent_atoms_shifted, aurostd::xvector<double>& previous_shift, aurostd::xvector<double>& new_shift);
  bool findWyckoffPositions(xstructure& CCell,
                            std::deque<_atom>& atomicbasis,
                            std::vector<std::vector<std::vector<std::string>>>& tmpvvvstring,
                            std::deque<std::deque<_atom>>& equivalent_atoms,
                            std::deque<std::deque<_atom>>& equivalent_atoms_shifted,
                            bool& foundspacegroup,
                            std::string& spacegroupstring,
                            bool& orig_origin_shift,
                            aurostd::xvector<double>& OriginShift,
                            std::vector<int>& wyckoffmult,
                            std::vector<std::string>& wyckoffsymbols,
                            std::vector<wyckoffsite_ITC>& wyckoffVariables,
                            std::deque<_atom>& wyckoffPositionsVector,
                            std::vector<std::string>& wyckoffSymbols,
                            std::ostringstream& woss,
                            bool& obverse_force_transformed);

  std::vector<std::vector<std::string>> getWyckoffEquations(const uint space_group_number, const std::string& space_group_setting, const std::string& Wyckoff_letter); // DX20191030
  std::vector<std::vector<std::string>> getWyckoffEquations(const std::string& Wyckoff_string, const std::string& Wyckoff_letter); // DX20191030
  uint getWyckoffMultiplicity(const uint space_group_number, const std::string& space_group_setting, const std::string& Wyckoff_letter); // DX20191030
  uint getWyckoffMultiplicity(const std::string& Wyckoff_string, const std::string& Wyckoff_letter); // DX20191030
  std::string getWyckoffSiteSymmetry(const uint space_group_number, const std::string& space_group_setting, const std::string& Wyckoff_letter); // DX20191030
  std::string getWyckoffSiteSymmetry(const std::string& Wyckoff_string, const std::string& Wyckoff_letter); // DX20191030
  void getWyckoffInformation(const uint space_group_number,
                             const std::string& space_group_setting,
                             const std::string& Wyckoff_letter,
                             uint& Wyckoff_multiplicity,
                             std::string& site_symmetry,
                             std::vector<std::vector<std::string>>& all_positions); // DX20191030
  void getWyckoffInformation(const std::string& Wyckoff_string, const std::string& Wyckoff_letter, uint& Wyckoff_multiplicity, std::string& site_symmetry, std::vector<std::vector<std::string>>& all_positions); // DX20191030

  std::vector<std::vector<std::vector<std::string>>> GetSameSymmetryWyckoffLetters(uint space_group_number, std::vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint setting);
  void print_wyckoff_pos(std::vector<std::vector<std::vector<std::string>>> wyckoff_positions);
  std::vector<std::vector<std::vector<std::vector<sdouble>>>> convert_wyckoff_pos_sd(std::vector<std::vector<std::vector<std::string>>> wyckoff_positions);
  void convert_wyckoff_pos(std::vector<std::vector<std::vector<std::string>>> wyckoff_positions);
  std::vector<std::vector<std::string>> get_centering(std::string spaceg);

  std::vector<double> ExtractLatticeParametersFromWyccar(const std::vector<std::string>& wyccar_ITC); // DX20191030 - added const
  std::string ExtractWyckoffAttributesString(const std::vector<std::string>& wyccar_ITC, uint attribute_index); // DX201780823 //DX20191030 - added const
  std::string ExtractWyckoffLettersString(const std::vector<std::string>& wyccar_ITC); // DX201780823 //DX20191030 - added const
  std::string ExtractWyckoffMultiplicitiesString(const std::vector<std::string>& wyccar_ITC); // DX201780823 //DX20191030 - added const
  std::string ExtractWyckoffSiteSymmetriesString(const std::vector<std::string>& wyccar_ITC); // DX201780823 //DX20191030 - added const
  std::vector<std::vector<std::vector<std::string>>> getWyckoffLettersWithSameMultiplcityAndSiteSymmetry(uint& space_group_number, std::vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint& cell_choice); // DX20190201
  std::vector<std::string> splitSiteSymmetry(const std::string& site_symmetry); // DX20190219 //DX20190730 - added const

  // TOPOLOGY FUNCTIONS
  std::vector<aurostd::xvector<double>> find_vectors_inplane(const std::vector<aurostd::xvector<double>>& big_expanded, const aurostd::xvector<double>& perp_to_vec, double& tol); // DX20190215

  // check if three distances can define a screw translation
  // shortest vectors. "m" specifies how many you want (e.g., 2 --> the smallest 2)
  std::vector<int> shortest_vec(std::vector<aurostd::xvector<double>> vecs, int m);
  // Closest point in terms of lattice:
  aurostd::xmatrix<double> get_mod_angle_tensor(std::vector<aurostd::xvector<double>>& points); // Points is a std::vector of points (ordered std::pair, triplet, etc.);
  // Check if two points are equivalent under the lattice L
  bool points_equivalent(aurostd::xmatrix<double>& c2f, aurostd::xmatrix<double>& f2c, aurostd::xvector<double> P1, aurostd::xvector<double> P2, aurostd::xvector<double>& lattice_vector, double& radius, bool& skew, double& tol); // DX20190215
  // Check if the candidate operation is a symmetry of the crystal C
  bool symmetry_axis_equivalent(aurostd::xmatrix<double> L, aurostd::xmatrix<double> Linv, aurostd::xvector<double> P1, aurostd::xvector<double> N1, aurostd::xvector<double> P2, aurostd::xvector<double> N2, double& tol); // DX20190215
  bool screw_equivalent(Screw S1, Screw S2);
  bool mirror_plane_equivalent(aurostd::xvector<double> normal);
  aurostd::xvector<double> next_point_on_line(aurostd::xvector<double> P, aurostd::xmatrix<double> L);
  void add_3d_point(std::vector<aurostd::xvector<double>>& points, double x, double y, double z);
  char discern_rot_sym(aurostd::xmatrix<double> m, double& tol); // DX20190215
  bool is_lattice_point(aurostd::xmatrix<double> L, aurostd::xvector<double> point, aurostd::xvector<double>& lattice_vector, double& radius, bool& skew, double& tol); // DX20190215

  double get_angle(aurostd::xvector<double> a, aurostd::xvector<double> b, std::string c);
  aurostd::xvector<double> get_random_point_in_plane(std::string plane);
  aurostd::xvector<double> get_point_in_line(std::string line, double param);
  aurostd::xvector<double> random_point();

  // STRING FUNCTIONS
  bool blank(std::string str_in);
  bool containschar(std::string str_in);
  bool iselem(std::string str_in);
  bool havechar(std::string str_in, char c); // see template std::function "invec"
  int whereischar(std::string str, char c);
  char whichchar(const std::string& str_in);
  double whichnum(std::string str_in);
  // DX20200313 [MOVED TO AUROSTD] double frac2dbl(std::string str);  //expand to cover case when input is e.g., ".5"
  // DX20190724 [MOVED TO AUROSTD] std::string dbl2frac(double a, bool sign_prefix=true);
  void multiply(std::vector<std::string> A, std::vector<std::string> B);
  void cleanupstring(std::string& str); // eliminates blank spaces before and after std::string
  std::vector<std::string> splitstring(std::string str, char c); // c is delimeter
  std::vector<std::string> splitstringbyspaces(std::string str);
  std::vector<std::string> splitstring(std::string str);
  std::vector<std::string> splitstring_2(std::string str);
  std::vector<sdouble> simplify(const std::string& str);
  aurostd::xvector<double> get_triplet(std::string str);

  // Used to check atomic basis
  bool GCD_conventional_atomic_basis(std::deque<_atom>& conventional_basis_atoms, std::deque<std::deque<_atom>>& prim_split_atom_types, int& prim_GCD);
  //[CO20180409 - moved to xatom]std::deque<_atom> foldAtomsInCell(std::deque<_atom>& atoms, aurostd::xmatrix<double>& c2f_new, aurostd::xmatrix<double>& f2c_new, bool& skew);
  //[CO20180409 - moved to xatom]std::deque<_atom> foldAtomsInCell(std::deque<_atom>& atoms, aurostd::xmatrix<double>& c2f_new, aurostd::xmatrix<double>& f2c_new, bool& skew, double& tol);
  bool MapAtomsInNewCell(_atom& a, _atom& b, aurostd::xmatrix<double>& lattice_new, bool& skew, double& tol); // DX20190619 - changed c2f_orig and f2c_new to lattice_new
  bool MapAtomsInNewCell(aurostd::xvector<double>& a, aurostd::xvector<double>& b, aurostd::xmatrix<double>& lattice_new, bool& skew, double& tol); // DX20190619 - changed c2f_orig and f2c_new to lattice_new
  std::deque<std::deque<_atom>> groupSymmetryEquivalentAtoms(std::deque<_atom>& atoms,
                                                             aurostd::xmatrix<double>& lattice,
                                                             std::vector<aurostd::xmatrix<double>>& sym_ops,
                                                             std::vector<aurostd::xvector<double>>& translations,
                                                             double& min_dist,
                                                             double& tol); // DX20190215
  std::deque<std::deque<_atom>> shiftSymmetryEquivalentAtoms(std::deque<std::deque<_atom>>& equivalent_atoms, aurostd::xmatrix<double>& lattice, aurostd::xvector<double>& translation, double& min_dist, double& tol);

  // ******************************************************************************
  // RSTD Namespace Functions
  // ******************************************************************************
  //  namespace rstd {
  typedef std::map<int, aurostd::xvector<double>> hash;
  aurostd::xmatrix<double> xvec2xmat(aurostd::xvector<double> a, aurostd::xvector<double> b, aurostd::xvector<double> c);
  aurostd::xmatrix<double> xvec2xmat(std::vector<aurostd::xvector<double>> V);
  aurostd::xmatrix<double> xvec2xmat(std::vector<aurostd::xvector<double>> V, std::vector<double> R);

  aurostd::xvector<double> extract_row(aurostd::xmatrix<double> a, int row);
  aurostd::xvector<double> dir(double a, double b, double c);
  aurostd::xmatrix<double> a2m3x3(double* array);

  double get_random_double(double min, double max);

  // SYMMETRY OPERATIONS FUNCTIONS
  symfolder check_ccell(xstructure& xstr, SymmetryInformationITC& ITC_sym_info); // DX20190215
  std::vector<aurostd::xvector<double>> expand_space_group_on_point(int sg, aurostd::xvector<double> point);
  std::vector<aurostd::xvector<double>> grid(double t);
  aurostd::xvector<double> find_inversion_point(double tol, std::vector<eqatoms> poscar_atoms);
  std::vector<aurostd::xvector<double>> symmetry_directions(char lattice_type);

  std::vector<Glide> mirror_operations(std::vector<aurostd::xvector<double>> expanded_lattice,
                                       std::vector<aurostd::xvector<double>> expanded_cell,
                                       aurostd::xmatrix<double> L,
                                       aurostd::xmatrix<double> Linv,
                                       std::vector<aurostd::xvector<double>>& lattice_vectors,
                                       double& radius,
                                       bool& skew,
                                       double& tol); // DX20190215 - added tol
  std::vector<Screw> triplet_operations(std::vector<aurostd::xvector<double>> expanded_lattice,
                                        std::vector<aurostd::xvector<double>> expanded_cell,
                                        aurostd::xmatrix<double> L,
                                        aurostd::xmatrix<double> Linv,
                                        std::vector<aurostd::xvector<double>>& lattice_vectors,
                                        double& radius,
                                        bool& skew,
                                        double& tol); // DX20190215 - added tol
  std::vector<Screw> twofold_operations(std::vector<aurostd::xvector<double>> expanded_lattice,
                                        std::vector<aurostd::xvector<double>> expanded_cell,
                                        aurostd::xmatrix<double> L,
                                        aurostd::xmatrix<double> Linv,
                                        std::vector<aurostd::xvector<double>>& lattice_vectors,
                                        double& radius,
                                        bool& skew,
                                        double& tol); // DX20190215 - added tol
  std::vector<aurostd::xvector<double>> getLatticeVectorsFromOriginalMirrorOperations(std::vector<Glide>& old_mirrors, std::vector<Glide>& new_mirrors, std::vector<aurostd::xvector<double>>& lattice_vectors, bool& all_matched);
  std::vector<aurostd::xvector<double>> getLatticeVectorsFromOriginalRotationOperations(std::vector<Screw>& old_rotations_twofold,
                                                                                        std::vector<Screw>& old_rotations_higher,
                                                                                        std::vector<Screw>& new_rotations,
                                                                                        std::vector<aurostd::xvector<double>>& twofold_lattice_vectors,
                                                                                        std::vector<aurostd::xvector<double>>& rot_lattice_vectors,
                                                                                        bool& all_matched);

  // COMBINATORICS FUNCTIONS
  void reduce_atom_deques(std::deque<_atom>& expanded, aurostd::xmatrix<double>& lattice, double& min_dist, double& sym_tol); // DX20190215
  // NOT IN SYM  std::vector<int> AllCombination41(int num, int total_num, int index);
  // NOT IN SYM  std::vector<int> AllCombination42(int num, int total_num, std::vector<int>& str_in);
  // NOT IN SYM  unsigned long int CombinationNr(int num, int total_num);
  std::vector<std::vector<int>> permute(int n);

  // LATTICE FUNCTIONS
  double length_lattice_vectors(aurostd::xmatrix<double> lattice);
  double orthogonality_defect(aurostd::xmatrix<double> xmat);

  // NUMBER THEORY FUNCTIONS
  double smallest_gt_min(double min, std::vector<double> vec);
  int smallest_gt_min_index(double min, int not_index1, int not_index2, std::vector<double> vec);

  // VISUALIZATION FUNCTIONS
  std::stringstream* latex_plot_lattice(aurostd::xmatrix<double> L, std::string color);
  std::string plot_lattice(std::vector<std::string> files);

  // ATOM CLASS ELEMENT MANIPULATION FUNCTIONS
  // Do (1-atm.coord) for the atomic basis, for the column associated with the lattice std::vector row
  void minus_one(std::deque<_atom>& atomicBasis, int lvec);
  // when lattice vectors are permuted, you must swap columns in basis (in direct)
  void swap_columns(std::deque<_atom>& atomicBasis, int col1, int col2);
  void rearrange_columns(std::deque<_atom>& atomicBasis, int c1, int c2, int c3);
} // namespace SYM
////A structure to store "decomposition" of points in plane--plane locations and distances from planes:
// struct Proj {
// std::vector<_atom> inplane_locations;
// std::vector<double> distances_from_plane;
// aurostd::xvector<double> plane_normal;
// };
void xb(); // print a break
template <class d> void print(std::vector<d> vec) {
  for (uint i = 0; i < vec.size(); i++) {
    std::cout << vec[i] << std::endl;
  }
}
void xb(); // print a break
void print(std::vector<int> vec);
void print(const std::vector<aurostd::xvector<double>>& points);
void print(const std::vector<std::vector<double>>& points);
void print(const std::vector<aurostd::xvector<double>>& points, aurostd::xmatrix<double> T);
void print(const std::vector<aurostd::xmatrix<double>>& mats);
void print(const std::vector<_atom>& atoms);
void print(const std::deque<_atom>& atoms);

// Function to eliminate duplicate projections:
// DX TEST void eliminate_duplicates(Proj& P);

namespace SYM {
  // Operations on class-atoms
  std::vector<int> count_types(std::deque<_atom>& vatom);

  std::vector<int> countWyckoffTypes(const std::vector<wyckoffsite_ITC>& Wyckoff_positions); // DX20210526

  // EXPAND LATTICE/CRYSTAL FUNCTIONS
  std::vector<aurostd::xvector<double>> expand_lattice_positive_only(int& a, int& b, int& c, aurostd::xmatrix<double>& L);
  std::vector<aurostd::xvector<double>> expand_lattice(int& a, int& b, int& c, aurostd::xmatrix<double>& L);
  std::vector<aurostd::xvector<double>> expand_cell(aurostd::xmatrix<double>& L);
  std::deque<_atom> add_basis(std::vector<aurostd::xvector<double>>& expanded_lattice_points, aurostd::xmatrix<double>& L, xstructure& xstr);

  // CONVENTIONAL LATTICE VECTOR FUNCTIONS
  void orient(xstructure& xstr, bool update_atom_positions = true);
  xstructure ConventionalCell(xstructure& xstr,
                              int& IT,
                              int& cell_choice,
                              bool& last_orientation,
                              std::string& crystalsystem_prev,
                              xstructure& CrystOut_prev,
                              std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors_prev,
                              std::vector<char>& candidate_lattice_chars_prev,
                              symfolder& checkops,
                              SymmetryInformationITC& ITC_sym_info,
                              bool& lattice_reformed,
                              std::vector<int>& lattice_pgroups,
                              std::vector<aurostd::xmatrix<double>>& lattice_sym_mats,
                              std::vector<aurostd::xmatrix<double>>& crystal_sym_mats,
                              bool& symmetry_found);
  bool findCubicLattice(std::vector<aurostd::xvector<double>>& rot_lattice_vectors,
                        std::vector<Screw>& rot_ops_vec,
                        std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors,
                        std::vector<char>& candidate_lattice_chars,
                        double& tol); // DX20190215 - added tol
  bool findTrigonalLattice(std::vector<aurostd::xvector<double>>& rot_lattice_vectors,
                           std::vector<Screw>& rot_ops_vec,
                           std::vector<aurostd::xvector<double>>& big_expanded,
                           std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors,
                           std::vector<char>& candidate_lattice_chars,
                           double& tol); // DX20190215 - added tol
  bool findTetragonalLattice(std::vector<aurostd::xvector<double>>& rot_lattice_vectors,
                             std::vector<aurostd::xvector<double>>& twofold_lattice_vectors,
                             std::vector<Screw>& rot_ops_vec,
                             std::vector<Screw>& twofold_ops_vec,
                             std::vector<aurostd::xvector<double>>& big_expanded,
                             std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors,
                             std::vector<char>& candidate_lattice_chars,
                             double& tol); // DX20190215
  bool findMonoclinicLattice(std::vector<aurostd::xvector<double>>& mirror_lattice_vectors,
                             std::vector<aurostd::xvector<double>>& twofold_lattice_vectors,
                             std::vector<aurostd::xvector<double>>& big_expanded,
                             std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors,
                             std::vector<char>& candidate_lattice_chars,
                             int& cell_choice,
                             double& tol); // DX20180816 - added cell_choice //DX20190215 - added tol
  bool findTriclinicLattice(aurostd::xmatrix<double>& lattice, std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors, std::vector<char>& candidate_lattice_chars);
  bool findOrthorhombicLattice(std::vector<aurostd::xvector<double>>& twofold_lattice_vectors,
                               std::vector<aurostd::xvector<double>>& mirror_lattice_vectors,
                               std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors,
                               std::vector<char>& candidate_lattice_chars,
                               double& tol); // DX20190215 - added tol
  bool findRhombohedralLattice(std::vector<aurostd::xvector<double>>& rot_lattice_vectors,
                               std::vector<Screw>& rot_ops_vec,
                               std::vector<aurostd::xvector<double>>& big_expanded,
                               std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors,
                               std::vector<char>& candidate_lattice_chars,
                               double& tol); // DX20190215 - added tol
  bool findRhombohedralSetting(std::vector<aurostd::xvector<double>>& big_expanded, std::vector<aurostd::xmatrix<double>>& candidate_lattice_vectors, std::vector<char>& candidate_lattice_chars, double& tol); // DX20190215 - added tol
  bool determineLatticeCentering(std::vector<aurostd::xvector<double>>& bravais_basis,
                                 int& bravais_count,
                                 aurostd::xmatrix<double>& c2f,
                                 aurostd::xmatrix<double>& f2c,
                                 bool& skew,
                                 std::vector<aurostd::xvector<double>>& big_expanded,
                                 std::string& crystalsystem,
                                 std::vector<char>& candidate_lattice_chars,
                                 double& tol); // DX20190215 - added tol
  std::string getPearsonSymbol(char& centering, char& lattice_char, std::deque<_atom> atoms);
  std::string spacegroup2latticeAndCentering(uint space_group_number); // DX20190418
  uint getEnantiomorphSpaceGroupNumber(uint space_group_number); // DX20181010
  bool getAtomGCD(std::deque<_atom>& atomic_basis, std::deque<std::deque<_atom>>& split_atom_types, int& gcd_num);
  void updateAtomPositions(std::deque<_atom>& atoms, Screw& S, aurostd::xmatrix<double>& lattice); // DX20190805 - return to void

  // RHOMBOHEDRAL OBVERSE/REVERSE FUNCTIONS
  bool isObverseSetting(xstructure& xstr, double& tolerance);
  bool isObverseSetting(aurostd::xmatrix<double>& lattice, std::deque<_atom>& atomic_basis_, double& dist_nn_min, double& tolerance);
  bool isObverseSetting(aurostd::xmatrix<double>& lattice, std::deque<std::deque<_atom>>& equivalent_atoms, double& dist_nn_min, double& tolerance);
  bool transformToObverse(aurostd::xmatrix<double>& lattice, std::deque<_atom>& atoms);
  bool transformToObverse(aurostd::xmatrix<double>& lattice, std::deque<std::deque<_atom>>& equivalent_atoms);

  // LATTICE VECTOR CHECK FUNCTIONS
  bool latticeVectorsSame(aurostd::xvector<double>& a, aurostd::xvector<double>& b, aurostd::xvector<double>& c, double& tol);
  bool orientVectorsPositiveQuadrant(std::vector<aurostd::xvector<double>>& lattice_vectors, double& tol);
  bool orientVectorsPositiveQuadrant(aurostd::xvector<double>& vec, double& tol);
  bool alignLatticeWithXYZ(aurostd::xvector<double>& a, aurostd::xvector<double>& b, aurostd::xvector<double>& c, double& tol);
  bool anyNullVectors(std::vector<aurostd::xvector<double>>& vecs, double& tol);
  bool nullVector(aurostd::xvector<double>& vec, double& tol);
} // namespace SYM

#endif
