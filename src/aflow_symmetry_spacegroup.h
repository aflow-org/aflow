// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Richard H. Taylor
// UPDATED by David Hicks
// d.hicks@duke.edu

#ifndef _AFLOW_SYMMETRY_SPACEGROUP_H_
#define _AFLOW_SYMMETRY_SPACEGROUP_H_
#include <cstring>
#include <map>
#include "aflow.h"
const double Pi_r = 3.141592653589793;

namespace SYM{
  //Symop
  struct symop {
    friend ostream& operator<<(ostream& output, const symop& a);
    string symbol;
    xvector<double> direction;
    xvector<double> screwglide;  //will also store inversion pnts for rotoinversion
    xvector<double> shift;

    void clear();
  };
}

class SymmetryInformationITC {
  public:
    SymmetryInformationITC();
    ~SymmetryInformationITC();
    friend ostream& operator<<(ostream& oss, const SymmetryInformationITC& SymmetryInformationITC);
    const SymmetryInformationITC& operator=(const SymmetryInformationITC& b);
    SymmetryInformationITC(const SymmetryInformationITC& b);

    //glides
    vector<xvector<double> > glideplanes;
    vector<xvector<double> > glideplanes_hex;
    vector<xvector<double> > glidetrans;
    vector<xvector<double> > glidetrans_hex;
    vector<string> glidesymbols;
    vector<string> glidesymbols_hex;

    //symmetry matrices
    vector<int> index_cubic;
    vector<int> index_hex;  //To be used with the *_hex vectors
    vector<int> index_rhom;
    vector<int> index_tetr;
    vector<int> index_ortho;
    vector<int> index_mono_b;
    vector<int> index_mono_c;
    vector<int> index_tric;
    vector<xmatrix<double> > sym_mats;
    vector<string> symbol;
    vector<string> dirparam;
    std::map<int, xvector<double> > sym_mats_direction;
    vector<xmatrix<double> > sym_mats_hex;
    vector<string> symbol_hex;
    vector<string> dirparam_hex;
    std::map<int, xvector<double> > sym_mats_direction_hex;

    //sym_ops
    vector<string> sym_ops;

    //generators
    vector<vector<SYM::symop> > generators;
    vector<int> sgindex;

    vector<string> gl_sgs;

    //functions
    bool initsymmats();
    bool initglides();
    bool initsymops();
    bool initgenerators(string axis_cell);
    bool initsgs(string axis_cell);
  private:
    void free();
    void copy(const SymmetryInformationITC& b);
};

//DX20190215 [OBSOLETE] // ******************************************************************************
//DX20190215 [OBSOLETE] // Tolerance Functions
//DX20190215 [OBSOLETE] // ******************************************************************************
//DX20190215 [OBSOLETE] namespace SYM{
//DX20190215 [OBSOLETE]   void SetTolerance(const double& starting_tol);
//DX20190215 [OBSOLETE]}

// ******************************************************************************
// Class Declarations
// ******************************************************************************
//xstructure
class xstructure;  //forward define xstructure for compilation.

vector<int> AllCombination41(int num, int total_num, int index);
vector<int> AllCombination42(int num, int total_num, vector<int>& str_in);
unsigned long int CombinationNr(int num, int total_num);

namespace SYM {
  //stringdouble (Used for Wyckoff position algebra)
  class stringdouble {
    friend ostream& operator<<(ostream& output, const stringdouble& a);

    public:
    string type;
    double distance;
    xvector<double> coord;
    stringdouble() {
      type = "XX";
      distance = 0;
    }
    ~stringdouble(){};
  };

  //wyckoff site
  class wyckoffsite  //Also for wyckoff sites
  {
    //friend bool operator==(const wyckoffsite& a, const wyckoffsite& b);
    friend ostream& operator<<(ostream& output, const wyckoffsite& a);

    public:
    wyckoffsite() {
      wyckoffSymbol = " ";
    };
    ~wyckoffsite(){};
    xvector<double> coord;
    string type;  //chemical label etc
    string wyckoffSymbol;
  };

  //atom class manipulation function
  //DX20210118 [OBSOLETE - now GetNumEachType in xatom] deque<int> arrange_atoms(deque<_atom>& atoms);

  //Screw (Screw Operation)
  class Screw {
    friend xvector<double> operator*(Screw& S, const xvector<double>& P);
    friend xvector<double> operator*(const xvector<double>& P, Screw& S);
    friend ostream& operator<<(ostream& output, const Screw& S);

    private:
    xmatrix<double> A;  //operator matrix
    void get_A();
    string linestring;          //parametric line
    xvector<double> T;          //translation component of screw
    double angle;               //angle of rotation
    xvector<double> one_point;  //a point on the axis (in get_A)
    xvector<double> direction_vector;

    public:
    Screw(){};
    ~Screw(){};
    void get_screw_direct(xvector<double> axis_dir, xvector<double> point, double order);
    void get_line_string(string str);
    void get_trans_comp(xvector<double> trans);
    void get_angle(double theta);
    xvector<double> return_direction();
    xvector<double> return_point();
    double return_order();
  };

  //Glide (Glide Operation)
  class Glide {  //for a plane defined on the plane ax+by+cz=d with translation
    //(t1,t2,t3) the glide operation works as follows:
    // A.P + d*a + |a.P - d|*a + (t1,t2,t3)^T
    private:
      bool HEX;
      bool DIRECT;
      string planestring;           //parametric plane
      xmatrix<double> A;            //reflecting matrix
      xvector<double> a;            //normal vector or fixed point (if HEX)
      xvector<double> T;            //translation vector
      xvector<double> plane_point;  //point in plane for get_glide_direct
      double d;                     // d in ax+by+cz+d
      void get_A(xvector<double> n);

    public:
      Glide() {
        HEX = false;
        DIRECT = false;
        d = 0.0;
      };
      ~Glide(){};
      //overload operator * for Glide
      // cartesian
      void get_glide_direct(xvector<double> n, xvector<double> p);
      void get_glide_direct(xvector<double> n, xvector<double> p, xvector<double> trans);
      xvector<double> return_point();      //returns point on plane
      xvector<double> return_direction();  //returns normal to plane
      friend xvector<double> operator*(Glide& G, const xvector<double>& P);
      friend xvector<double> operator*(const xvector<double>& P, Glide& G);
      friend ostream& operator<<(ostream& output, const Glide& S);
  };

  //Translation
  class Translation {
    private:
      xvector<double> translation_vec;

    public:
      Translation(){};
      ~Translation(){};
      void get_translation(string ITCstring);
      friend xvector<double> operator*(Translation& T, const xvector<double>& P);
      friend xvector<double> operator*(const xvector<double>& P, Translation& T);
      friend xvector<double> operator+(Translation& T, const xvector<double>& P);
      friend xvector<double> operator+(const xvector<double>& P, Translation& T);
  };

  //Inversion
  class Inversion {
    friend xvector<double> operator*(Inversion& I, const xvector<double>& P);
    friend xvector<double> operator*(const xvector<double>& P, Inversion& I);

    private:
    public:
    Inversion(){};
    ~Inversion(){};
    //overload operator * for Glide
    xvector<double> inversion_point;
    void get_inversion(string ITCstring);
  };
} //namespace SYM

//End Classes
// ******************************************************************************

// ******************************************************************************
// Structure Declarations
// ******************************************************************************

namespace SYM {
  //Symfolder (Stores all symmetry elements/operators)
  struct symfolder {
    vector<Screw> twofold_ops;
    vector<Screw> rot_ops;
    vector<Glide> mirror_ops;
    vector<xvector<double> > twofold_lattice_vectors;
    vector<xvector<double> > rot_lattice_vectors;
    vector<xvector<double> > mirror_lattice_vectors;
    bool commensurate;
    string crystalsystem;
    string latticesystem;
    vector<int> lattice_pgroups;                //DX NEW
    vector<xmatrix<double> > lattice_sym_mats;  //DX NEW
    vector<xmatrix<double> > crystal_sym_mats;  //DX NEW
    vector<int> insym;
    vector<xvector<double> > translations;
    vector<vector<int> > all_atom_maps;
    vector<vector<int> > all_type_maps;
    vector<xvector<double> > centeringops;
    vector<string> pointgroupops;
  };

  ////Symop
  //struct symop {
  //  friend ostream& operator<<(ostream& output, const symop& a);
  //  string symbol;
  //  xvector<double> direction;
  //  xvector<double> screwglide;  //will also store inversion pnts for rotoinversion
  //  xvector<double> shift;
  //
  //    void clear();
  //  };

  //Eqatoms
  struct eqatoms {
    string label;
    vector<vector<double> > atoms;
  };

  //sdouble (String and double)
  struct sdouble {
    char chr;
    double dbl;
    sdouble() {
      chr = '\0';
      dbl = 0;
    }
  };

  //enum_alpha (enumerate alphabet: for Wyckoff letters)
  struct enum_alpha {
    string letter;
    int number;
  };


  //End Structures
  // ******************************************************************************

  // ******************************************************************************
  // Template Declarations
  // ******************************************************************************

  template <class dmmy>
    bool allsame(vector<dmmy> vec) {
      bool all = true;
      for (uint i = 0; i < vec.size(); i++) {
        if (vec[i] != vec[0])
          all = false;
      }
      return all;
    }

} //namespace SYM

//End Templates
// ******************************************************************************

// ******************************************************************************
// Function Declarations
// ******************************************************************************
//MAIN FUNCITONS
namespace SYM {
  void calculateSpaceGroups(vector<xstructure>& vxstrs, uint start_index=0, uint end_index=AUROSTD_MAX_UINT, uint setting=0); //DX20191230 add setting option
  string OrthoDefect(istream& cin);
  xstructure SpaceGroup(istream& cin);
  void rgcd(vector<string> num);
  //DX void AtomicEnvironment(istream & cin, vector<string> av);
  string ReverseSpaceGroup(vector<string> num);
  //END MAIN FUNCTIONS

  //FUNCTION THAT TAKES WYCCAR FROM XSTRUCTURE AND PUTS IT IN OSTREAM
  void printWyccar(ofstream& FileMESSAGE, xstructure& str, const bool& osswrite, ostream& oss);

  //LINEAR ALGEBRA FUNCTIONS
  bool solve_overdetermined_system(vector<xvector<double> >& LHS, vector<double>& RHS, xvector<double>& SOL, xmatrix<double>& lattice, double& min_dist, double& tol); //DX20190215
  void ReducedRowEchelonForm(xmatrix<double>& M, double& tol); //DX20190215
  bool find_solution_UD(xmatrix<double> M, xvector<double>& SOL, double& tol);  //check if underdetermined system has a solution, and get a particular solution. //DX20190215
  bool checkLinearSystem(vector<xvector<double> >& LHS, vector<double>& RHS, xmatrix<double>& lattice, double& tol); //DX20190215
  //WYCKOFF FUNCITONS
  vector<double> system_solve(vector<double> numeric, vector<vector<sdouble> > variable);
  vector<vector<double> > wyckoff_solve(vector<vector<vector<sdouble> > > win, vector<eqatoms> pin);
  vector<vector<double> > wyckoff_solve_2(vector<vector<vector<sdouble> > > win, eqatoms pin);
  vector<vector<double> > wyckoff_sites(vector<vector<vector<sdouble> > > win, vector<_atom> atomgroup);

  //SPACE GROUP LIBRARY FUNCTIONS
  bool initsgs(string axis_cell);
  bool initsymops();
  bool initsymmats();
  bool initglides();
  bool initgenerators(string axis_cell);
  vector<int> generatorindex(int spacegroup);
  vector<xvector<double> > ReturnITCGenShift(int spacegroupnum, string axis_cell);

  string point_group_library_search(vector<string> symmetryoperations, string crystalsystem, int centeringops);
  string crystal_system_library_search(vector<string> symmetryoperations);
  vector<int> spacegroups_CS_LS(string crystalsystem, char latticesystem);
  char discern_sym_op(string ITCstring);
  int* PointGroup_SpaceGroup(string pgroup);
  vector<int> PointGroup_SpaceGroup(string pgroup, char cl);
  xmatrix<double> spacegroup_to_lattice(int spacegroup);

  vector<int> get_multiplicities(string sg);
  vector<string> get_symmetry_symbols(string sg);

  vector<string> get_wyckoff_equation(string spaceg, int mult);  //DX20170830
  vector<vector<vector<string> > > get_wyckoff_pos(string spaceg, int mult);
  vector<vector<vector<string> > > get_wyckoff_pos(string spaceg, int mult, bool getcentering);
  vector<vector<string> > get_wyckoff_pos(string spaceg, int Wyckoff_multiplicity, string Wyckoff_letter); //DX20190129
  vector<string> get_minimum_enumerated_Wyckoff_letters(string spacegroupstring, vector<int>& multiplicities, vector<string> site_symmetries);
  int enumerate_wyckoff_letter(string& wyckoff_letter);
  vector<int> enumerate_wyckoff_letters(vector<string>& wyckoff_letters); //DX20180927
  void get_all_wyckoff_for_site_symmetry(string spaceg, int mult, string site_symmetry, vector<vector<string> >& all_positions);
  void get_Wyckoff_from_letter(uint space_group_number, string& space_group_setting,
      string& Wyckoff_letter, uint& Wyckoff_multiplicity, string& site_symmetry, vector<string>& positions); //DX20191029
  void get_Wyckoff_from_letter(string& spaceg, string& Wyckoff_letter, 
      uint& Wyckoff_multiplicity, string& site_symmetry, vector<string>& positions);
  xvector<double> Wyckoff_position_string2xvector(string& string_position);
  vector<xvector<double> > get_possible_origin_shifts(string spacegroupstring, int multiplicity, string site_symmetry);
  void get_certain_wyckoff_pos(string spaceg, int mult, string site_symmetry, vector<string>& site_symmetries, vector<string>& letters, vector<string>& positions);
  void getWyckoffAutomorphismSets(const string& spaceg, int mult, string site_symmetry, vector<string>& site_symmetries, vector<string>& letters, vector<vector<string> >& all_positions, bool origin_shifts_only=false); //DX20220828
  void getGeneralWyckoffMultiplicityAndPosition(uint space_group_number, string& space_group_setting, int& general_wyckoff_multiplicity, vector<string>& general_wyckoff_position);
  vector<string> findGeneralWyckoffPosition(string& spacegroupstring, int& general_wyckoff_multiplicity);
  uint numberOccupiedSitesInConventionalCell(const vector<wyckoffsite_ITC>& Wyckoff_sites); //DX20200512
  vector<uint> numberEachTypeFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_sites); //DX20200512
  vector<string> findWyckoffEquations(uint space_group_number, string& space_group_setting,
      string& Wyckoff_letter, uint Wyckoff_multiplicity, bool convert2frac=true); //DX 20191029 //DX20200423 - add convert2frac
  vector<string> findWyckoffEquations(string& spacegroupstring, string& Wyckoff_letter, uint Wyckoff_multplicity, bool convert2frac=true, bool keep_multiplication_symbol=false); //DX 20190128 //DX20200423 - add convert2frac, keep_multiplication_symbol
  string formatWyckoffPosition(const vector<sdouble>& sd_coordinate, bool convert2frac=true, bool keep_multiplication_symbol=false,int precision=AUROSTD_DEFAULT_PRECISION); //DX 20190723 //DX20200423 - add convert2frac, keep_multiplication_symbol  //CO20220607 - added precision
  string reorderWyckoffPosition(const string& orig_position); //DX 20190708
  bool shiftWyckoffPositions(deque<deque<_atom> >& equivalent_atoms_shifted, xvector<double>& previous_shift, xvector<double>& new_shift);
  bool findWyckoffPositions(xstructure& CCell, deque<_atom>& atomicbasis, vector<vector<vector<string> > >& tmpvvvstring,
      deque<deque<_atom> >& equivalent_atoms, deque<deque<_atom> >& equivalent_atoms_shifted,
      bool& foundspacegroup, string& spacegroupstring, bool& orig_origin_shift, xvector<double>& OriginShift,
      vector<int>& wyckoffmult, vector<string>& wyckoffsymbols, vector<wyckoffsite_ITC>& wyckoffVariables,
      deque<_atom>& wyckoffPositionsVector, vector<string>& wyckoffSymbols, ostringstream& woss,
      bool& obverse_force_transformed);

  vector<vector<string> > getWyckoffEquations(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter); //DX20191030
  vector<vector<string> > getWyckoffEquations(const string& Wyckoff_string, const string& Wyckoff_letter); //DX20191030
  uint getWyckoffMultiplicity(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter); //DX20191030
  uint getWyckoffMultiplicity(const string& Wyckoff_string, const string& Wyckoff_letter); //DX20191030
  string getWyckoffSiteSymmetry(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter); //DX20191030
  string getWyckoffSiteSymmetry(const string& Wyckoff_string, const string& Wyckoff_letter); //DX20191030
  void getWyckoffInformation(const uint space_group_number, const string& space_group_setting, const string& Wyckoff_letter,
      uint& Wyckoff_multiplicity, string& site_symmetry, vector<vector<string> >& all_positions); //DX20191030
  void getWyckoffInformation(const string& Wyckoff_string, const string& Wyckoff_letter,
      uint& Wyckoff_multiplicity, string& site_symmetry, vector<vector<string> >& all_positions); //DX20191030

  vector<vector<vector<string> > > GetSameSymmetryWyckoffLetters(uint space_group_number, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint setting);
  void print_wyckoff_pos(vector<vector<vector<string> > > wyckoff_positions);
  vector<vector<vector<vector<sdouble> > > > convert_wyckoff_pos_sd(vector<vector<vector<string> > > wyckoff_positions);
  void convert_wyckoff_pos(vector<vector<vector<string> > > wyckoff_positions);
  vector<vector<string> > get_centering(string spaceg);

  vector<double> ExtractLatticeParametersFromWyccar(const vector<string>& wyccar_ITC); //DX20191030 - added const
  string ExtractWyckoffAttributesString(const vector<string>& wyccar_ITC, uint attribute_index); //DX201780823 //DX20191030 - added const 
  string ExtractWyckoffLettersString(const vector<string>& wyccar_ITC); //DX201780823 //DX20191030 - added const
  string ExtractWyckoffMultiplicitiesString(const vector<string>& wyccar_ITC); //DX201780823 //DX20191030 - added const
  string ExtractWyckoffSiteSymmetriesString(const vector<string>& wyccar_ITC); //DX201780823 //DX20191030 - added const
  vector<vector<vector<string> > > getWyckoffLettersWithSameMultiplcityAndSiteSymmetry(uint& space_group_number, 
      vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint& cell_choice); //DX20190201  
  vector<string> splitSiteSymmetry(const string& site_symmetry); //DX20190219 //DX20190730 - added const

  //TOPOLOGY FUNCTIONS
  vector<xvector<double> > find_vectors_inplane(const vector<xvector<double> >& big_expanded, const xvector<double>& perp_to_vec, double& tol); //DX20190215

  //check if three distances can define a screw translation
  //shortest vectors. "m" specifies how many you want (e.g., 2 --> the smallest 2)
  vector<int> shortest_vec(vector<xvector<double> > vecs, int m);
  //Closest point in terms of lattice:
  //[OBSOLETE] xvector<double> closest_point(xvector<double> norm, xmatrix<double> Linv, xmatrix<double> L);
  xmatrix<double> get_mod_angle_tensor(vector<xvector<double> >& points);  //Points is a vector of points (ordered pair, triplet, etc.);
  //Check if two points are equivalent under the lattice L
  bool points_equivalent(xmatrix<double>& c2f, xmatrix<double>& f2c, xvector<double> P1, xvector<double> P2, xvector<double>& lattice_vector, double& radius, bool& skew, double& tol); //DX20190215
  //Check if the candidate operation is a symmetry of the crystal C
  bool symmetry_axis_equivalent(xmatrix<double> L, xmatrix<double> Linv, xvector<double> P1, xvector<double> N1, xvector<double> P2, xvector<double> N2, double& tol); //DX20190215
  bool screw_equivalent(Screw S1, Screw S2);
  bool mirror_plane_equivalent(xvector<double> normal);
  xvector<double> next_point_on_line(xvector<double> P, xmatrix<double> L);
  void add_3d_point(vector<xvector<double> >& points, double x, double y, double z);
  char discern_rot_sym(xmatrix<double> m, double& tol); //DX20190215
  //DX20210422 [OBSOLETE] bool allsame(vector<double> v);
  bool is_lattice_point(xmatrix<double> L, xvector<double> point, xvector<double>& lattice_vector, double& radius, bool& skew, double& tol); //DX20190215
  //DX20210422 [OBSOLETE] double distance_between_points(const xvector<double>& a, const xvector<double>& b);

  double get_angle(xvector<double> a, xvector<double> b, string c);
  xvector<double> get_random_point_in_plane(string plane);
  xvector<double> get_point_in_line(string line, double param);
  xvector<double> random_point();

  //STRING FUNCTIONS
  bool blank(string str_in);
  bool containschar(string str_in);
  bool iselem(string str_in);
  bool havechar(string str_in, char c);  //see template function "invec"
  //DX20210421 [OBSOLETE] bool intinvec(vector<int> vint, int n);
  int whereischar(string str, char c);
  char whichchar(const string& str_in);
  double whichnum(string str_in);
  //DX20200313 [MOVED TO AUROSTD] double frac2dbl(string str);  //expand to cover case when input is e.g., ".5"
  //DX20190724 [MOVED TO AUROSTD] string dbl2frac(double a, bool sign_prefix=true);
  void multiply(vector<string> A, vector<string> B);
  //DX20210422 [OBSOLETE] void xstring(ostream& output, xmatrix<double> a);
  void cleanupstring(string& str);  //eliminates blank spaces before and after string
  vector<string> splitstring(string str, char c);  //c is delimeter
  vector<string> splitstringbyspaces(string str);
  vector<string> splitstring(string str);
  vector<string> splitstring_2(string str);
  vector<sdouble> simplify(const string& str);
  xvector<double> get_triplet(string str);

  // Used to check atomic basis
  bool GCD_conventional_atomic_basis(deque<_atom>& conventional_basis_atoms, deque<deque<_atom> >& prim_split_atom_types, int& prim_GCD);
  //[CO20180409 - moved to xatom]deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew);
  //[CO20180409 - moved to xatom]deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool& skew, double& tol);
  bool MapAtomsInNewCell(_atom& a, _atom& b, xmatrix<double>& lattice_new, bool& skew, double& tol); //DX20190619 - changed c2f_orig and f2c_new to lattice_new
  bool MapAtomsInNewCell(xvector<double>& a, xvector<double>& b, xmatrix<double>& lattice_new, bool& skew, double& tol); //DX20190619 - changed c2f_orig and f2c_new to lattice_new
  deque<deque<_atom> > groupSymmetryEquivalentAtoms(deque<_atom>& atoms, xmatrix<double>& lattice, vector<xmatrix<double> >& sym_ops,
      vector<xvector<double> >& translations, double& min_dist, double& tol); //DX20190215
  deque<deque<_atom> > shiftSymmetryEquivalentAtoms(deque<deque<_atom> >& equivalent_atoms, xmatrix<double>& lattice, xvector<double>& translation, double& min_dist, double& tol);

  // ******************************************************************************
  // RSTD Namespace Functions
  // ******************************************************************************
  //  namespace rstd {
  typedef std::map<int, xvector<double> > hash;
  xmatrix<double> xvec2xmat(xvector<double> a, xvector<double> b, xvector<double> c);
  xmatrix<double> xvec2xmat(vector<xvector<double> > V);
  xmatrix<double> xvec2xmat(vector<xvector<double> > V, vector<double> R);

  xvector<double> extract_row(xmatrix<double> a, int row);
  xvector<double> dir(double a, double b, double c);
  xmatrix<double> a2m3x3(double* array);

  double get_random_double(double min, double max);

  //SYMMETRY OPERATIONS FUNCTIONS
  //DX20190215 [OBSOLETE] symfolder check_ccell(xstructure& xstr);
  symfolder check_ccell(xstructure& xstr, SymmetryInformationITC& ITC_sym_info); //DX20190215
  vector<xvector<double> > expand_space_group_on_point(int sg, xvector<double> point);
  vector<xvector<double> > grid(double t);
  xvector<double> find_inversion_point(double tol, vector<eqatoms> poscar_atoms);
  vector<xvector<double> > symmetry_directions(char lattice_type);

  vector<Glide> mirror_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew, double& tol); //DX20190215 - added tol
  vector<Screw> triplet_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew, double& tol); //DX20190215 - added tol
  vector<Screw> twofold_operations(vector<xvector<double> > expanded_lattice, vector<xvector<double> > expanded_cell, xmatrix<double> L, xmatrix<double> Linv, vector<xvector<double> >& lattice_vectors, double& radius, bool& skew, double& tol); //DX20190215 - added tol
  vector<xvector<double> > getLatticeVectorsFromOriginalMirrorOperations(vector<Glide>& old_mirrors, vector<Glide>& new_mirrors,
      vector<xvector<double> >& lattice_vectors, bool& all_matched);
  vector<xvector<double> > getLatticeVectorsFromOriginalRotationOperations(vector<Screw>& old_rotations_twofold,
      vector<Screw>& old_rotations_higher, vector<Screw>& new_rotations,
      vector<xvector<double> >& twofold_lattice_vectors,
      vector<xvector<double> >& rot_lattice_vectors, bool& all_matched);

  // COMBINATORICS FUNCTIONS
  void reduce_atom_deques(deque<_atom>& expanded, xmatrix<double>& lattice, double& min_dist, double& sym_tol); //DX20190215
  //NOT IN SYM  vector<int> AllCombination41(int num, int total_num, int index);
  //NOT IN SYM  vector<int> AllCombination42(int num, int total_num, vector<int>& str_in);
  //NOT IN SYM  unsigned long int CombinationNr(int num, int total_num);
  vector<vector<int> > permute(int n);

  //LATTICE FUNCTIONS
  double length_lattice_vectors(xmatrix<double> lattice);
  double orthogonality_defect(xmatrix<double> xmat);

  //NUMBER THEORY FUNCTIONS
  //DX20191202 [OBSOLETE] long long int gcd(long long int u, long long int v);                             //Euclid's
  //DX20191202 [OBSOLETE] unsigned long long int gcd(unsigned long long int u, unsigned long long int v);  //Dijkstra's GCD Algorithm
  //DX20191202 [OBSOLETE] int gcdD(int u, int v);
  //[OBSOLETE]long long int cast2int(double d, long long int prec);
  double smallest_gt_min(double min, vector<double> vec);
  int smallest_gt_min_index(double min, int not_index1, int not_index2, vector<double> vec);

  //VISUALIZATION FUNCTIONS
  stringstream* latex_plot_lattice(xmatrix<double> L, string color);
  string plot_lattice(vector<string> files);

  //ATOM CLASS ELEMENT MANIPULATION FUNCTIONS
  //Do (1-atm.coord) for the atomic basis, for the column associated with the lattice vector row
  void minus_one(deque<_atom>& atomicBasis, int lvec);
  //when lattice vectors are permuted, you must swap columns in basis (in direct)
  void swap_columns(deque<_atom>& atomicBasis, int col1, int col2);
  void rearrange_columns(deque<_atom>& atomicBasis, int c1, int c2, int c3);
} //namespace SYM
////A structure to store "decomposition" of points in plane--plane locations and distances from planes:
//struct Proj {
//vector<_atom> inplane_locations;
//vector<double> distances_from_plane;
//xvector<double> plane_normal;
//};
void xb();  //print a break
template <class d>
void print(vector<d> vec) {
  for (uint i = 0; i < vec.size(); i++) {
    cout << vec[i] << endl;
  }
}
void xb();  //print a break
void print(vector<int> vec);
void print(const vector<xvector<double> >& points);
void print(const vector<vector<double> >& points);
void print(const vector<xvector<double> >& points, xmatrix<double> T);
void print(const vector<xmatrix<double> >& mats);
void print(const vector<_atom>& atoms);
void print(const deque<_atom>& atoms);

//Function to eliminate duplicate projections:
//DX TEST void eliminate_duplicates(Proj& P);

namespace SYM {
  //Operations on class-atoms
  vector<int> count_types(deque<_atom>& vatom);

  vector<int> countWyckoffTypes(const vector<wyckoffsite_ITC>& Wyckoff_positions); //DX20210526

  //EXPAND LATTICE/CRYSTAL FUNCTIONS
  vector<xvector<double> > expand_lattice_positive_only(int& a, int& b, int& c, xmatrix<double>& L);
  vector<xvector<double> > expand_lattice(int& a, int& b, int& c, xmatrix<double>& L);
  vector<xvector<double> > expand_cell(xmatrix<double>& L);
  deque<_atom> add_basis(vector<xvector<double> >& expanded_lattice_points, xmatrix<double>& L, xstructure& xstr);

  // CONVENTIONAL LATTICE VECTOR FUNCTIONS
  void orient(xstructure& xstr, bool update_atom_positions=true);
  xstructure ConventionalCell(xstructure& xstr, int& IT, int& cell_choice, bool& last_orientation, string& crystalsystem_prev,
      xstructure& CrystOut_prev, vector<xmatrix<double> >& candidate_lattice_vectors_prev,
      vector<char>& candidate_lattice_chars_prev, symfolder& checkops, SymmetryInformationITC& ITC_sym_info, bool& lattice_reformed,
      vector<int>& lattice_pgroups, vector<xmatrix<double> >& lattice_sym_mats,
      vector<xmatrix<double> >& crystal_sym_mats, bool& symmetry_found);
  bool findCubicLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol); //DX20190215 - added tol
  bool findTrigonalLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec, vector<xvector<double> >& big_expanded,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol); //DX20190215 - added tol
  bool findTetragonalLattice(vector<xvector<double> >& rot_lattice_vectors, vector<xvector<double> >& twofold_lattice_vectors,
      vector<Screw>& rot_ops_vec, vector<Screw>& twofold_ops_vec, vector<xvector<double> >& big_expanded,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol); //DX20190215
  bool findMonoclinicLattice(vector<xvector<double> >& mirror_lattice_vectors, vector<xvector<double> >& twofold_lattice_vectors,
      vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars, int& cell_choice, double& tol); //DX20180816 - added cell_choice //DX20190215 - added tol
  bool findTriclinicLattice(xmatrix<double>& lattice, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars);
  bool findOrthorhombicLattice(vector<xvector<double> >& twofold_lattice_vectors, vector<xvector<double> >& mirror_lattice_vectors,
      vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol); //DX20190215 - added tol
  bool findRhombohedralLattice(vector<xvector<double> >& rot_lattice_vectors, vector<Screw>& rot_ops_vec,
      vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors, vector<char>& candidate_lattice_chars, double& tol); //DX20190215 - added tol
  bool findRhombohedralSetting(vector<xvector<double> >& big_expanded, vector<xmatrix<double> >& candidate_lattice_vectors,
      vector<char>& candidate_lattice_chars, double& tol); //DX20190215 - added tol
  bool determineLatticeCentering(vector<xvector<double> >& bravais_basis, int& bravais_count, xmatrix<double>& c2f, xmatrix<double>& f2c, bool& skew, vector<xvector<double> >& big_expanded, string& crystalsystem, vector<char>& candidate_lattice_chars, double& tol); //DX20190215 - added tol
  string getPearsonSymbol(char& centering, char& lattice_char, deque<_atom> atoms);
  string spacegroup2latticeAndCentering(uint space_group_number); //DX20190418
  uint getEnantiomorphSpaceGroupNumber(uint space_group_number); //DX20181010
  bool getAtomGCD(deque<_atom>& atomic_basis, deque<deque<_atom> >& split_atom_types, int& gcd_num);
  void updateAtomPositions(deque<_atom>& atoms, Screw& S, xmatrix<double>& lattice); //DX20190805 - return to void

  //RHOMBOHEDRAL OBVERSE/REVERSE FUNCTIONS
  bool isObverseSetting(xstructure& xstr, double& tolerance);
  bool isObverseSetting(xmatrix<double>& lattice, deque<_atom>& atomic_basis_, double& dist_nn_min, double& tolerance);
  bool isObverseSetting(xmatrix<double>& lattice, deque<deque<_atom> >& equivalent_atoms, double& dist_nn_min, double& tolerance);
  bool transformToObverse(xmatrix<double>& lattice, deque<_atom>& atoms);
  bool transformToObverse(xmatrix<double>& lattice, deque<deque<_atom> >& equivalent_atoms);

  //LATTICE VECTOR CHECK FUNCTIONS
  bool latticeVectorsSame(xvector<double>& a, xvector<double>& b, xvector<double>& c, double& tol);
  bool orientVectorsPositiveQuadrant(vector<xvector<double> >& lattice_vectors, double& tol);
  bool orientVectorsPositiveQuadrant(xvector<double>& vec, double& tol);
  bool alignLatticeWithXYZ(xvector<double>& a, xvector<double>& b, xvector<double>& c, double& tol);
  bool anyNullVectors(vector<xvector<double> >& vecs, double& tol);
  bool nullVector(xvector<double>& vec, double& tol);
} //namespace SYM

//DX20190215 [OBSOLETE] // external variables
//DX20190215 [OBSOLETE] namespace SYM {
//DX20190215 [OBSOLETE] #ifdef AFLOW_SYMMETRY_MULTITHREADS_ENABLE
//DX20190215 [OBSOLETE]   extern thread_local vector<xvector<double> > glideplanes;
//DX20190215 [OBSOLETE]   extern thread_local vector<xvector<double> > glidetrans;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> glidesymbols;
//DX20190215 [OBSOLETE]   extern thread_local vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> symbol;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> dirparam;
//DX20190215 [OBSOLETE]   extern thread_local hash sym_mats_direction;
//DX20190215 [OBSOLETE]   extern thread_local hash sym_mats_direction_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<xmatrix<double> > sym_mats_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> symbol_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> dirparam_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<xvector<double> > glideplanes_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<xvector<double> > glidetrans_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> glidesymbols_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_cubic;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_hex;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_rhom;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_tetr;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_ortho;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_mono_b;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_mono_c;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> index_tric;
//DX20190215 [OBSOLETE]   extern thread_local vector<vector<symop> > generators;
//DX20190215 [OBSOLETE]   extern thread_local vector<int> sgindex;
//DX20190215 [OBSOLETE]   extern thread_local vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> symbol;
//DX20190215 [OBSOLETE]   extern thread_local hash sym_mats_direction;
//DX20190215 [OBSOLETE]   extern thread_local vector<string> gl_sgs;
//DX20190215 [OBSOLETE] #else
//DX20190215 [OBSOLETE]   extern vector<xvector<double> > glideplanes;
//DX20190215 [OBSOLETE]   extern vector<xvector<double> > glidetrans;
//DX20190215 [OBSOLETE]   extern vector<string> glidesymbols;
//DX20190215 [OBSOLETE]   extern vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE]   extern vector<string> symbol;
//DX20190215 [OBSOLETE]   extern vector<string> dirparam;
//DX20190215 [OBSOLETE]   extern hash sym_mats_direction;
//DX20190215 [OBSOLETE]   extern hash sym_mats_direction_hex;
//DX20190215 [OBSOLETE]   extern vector<xmatrix<double> > sym_mats_hex;
//DX20190215 [OBSOLETE]   extern vector<string> symbol_hex;
//DX20190215 [OBSOLETE]   extern vector<string> dirparam_hex;
//DX20190215 [OBSOLETE]   extern vector<xvector<double> > glideplanes_hex;
//DX20190215 [OBSOLETE]   extern vector<xvector<double> > glidetrans_hex;
//DX20190215 [OBSOLETE]   extern vector<string> glidesymbols_hex;
//DX20190215 [OBSOLETE]   extern vector<int> index_cubic;
//DX20190215 [OBSOLETE]   extern vector<int> index_hex;
//DX20190215 [OBSOLETE]   extern vector<int> index_rhom;
//DX20190215 [OBSOLETE]   extern vector<int> index_tetr;
//DX20190215 [OBSOLETE]   extern vector<int> index_ortho;
//DX20190215 [OBSOLETE]   extern vector<int> index_mono_b;
//DX20190215 [OBSOLETE]   extern vector<int> index_mono_c;
//DX20190215 [OBSOLETE]   extern vector<int> index_tric;
//DX20190215 [OBSOLETE]   extern vector<vector<symop> > generators;
//DX20190215 [OBSOLETE]   extern vector<int> sgindex;
//DX20190215 [OBSOLETE]   extern vector<xmatrix<double> > sym_mats;
//DX20190215 [OBSOLETE]   extern vector<string> symbol;
//DX20190215 [OBSOLETE]   extern hash sym_mats_direction;
//DX20190215 [OBSOLETE]   extern vector<string> gl_sgs;
//DX20190215 [OBSOLETE] 
//DX20190215 [OBSOLETE] #endif
//DX20190215 [OBSOLETE] }

#endif
