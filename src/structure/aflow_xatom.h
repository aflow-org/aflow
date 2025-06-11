#ifndef AFLOW_XATOM_H
#define AFLOW_XATOM_H

#include <deque>
#include <iosfwd>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xcomplex.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"

class _atom : public JsonSerializable<_atom> {
  // simple class... nothing fancy
public:
  // constructor destructor                              // constructor/destructor
  _atom(); // default, just allocate
  ~_atom(); // kill everything
  _atom(const _atom& b); // constructor copy
  const _atom& operator=(const _atom& b); // copy
  void clear();
  // content                                             // content
  aurostd::xvector<double> fpos; // positions are with respect to ijk lattice cell
  aurostd::xvector<double> cpos; // so if fpos/cpos is outside a cell, you can shift
  aurostd::xvector<double> corigin; // origin for convasp purposes
  aurostd::xvector<double> coord; // general coordinate for symmetry routines (RHT)
  std::vector<std::string> fpos_equation; // DX20180607 - lattice equation for atomic position
  std::vector<std::string> cpos_equation; // DX20180607 - Cartesian equation for atomic position
  double spin; // spin along z in VASP MODE
  bool spin_is_given; // true if spin has been set //DX20170921
  aurostd::xvector<double> noncoll_spin; // non-collinear spin                //DX20171205
  bool noncoll_spin_is_given; // true if noncoll_spin has been set //DX20171205
  double mass; // mass
  int type; // with bringincell, which adjust cpos/fpos and ijk as well
  std::string name; // the name read from the INPUT
  bool name_is_given; // true is atom name has been given
  std::string cleanname; // a chemical clean version of the name
  int info; // container for misc. information  //RHT
  int atomic_number; // 0 by defauls
  //[CO20200130 - number->basis]int    number;                                         // atom number reference for convasp, from zero to the sky
  std::string sd; // ?
  aurostd::xvector<int> ijk; // xvector identifier of the lattice (but you must give str)
  bool isincell; // is in cell ? (i==j==k==0 ?)
  int basis; // identifier of position in the basis, from zero to the sky
  double reference; // reference/measure for ordering
  int ireference; // sort of number in the list
  // for symmetry
  int equivalent; // points to the equivalent atom in the cell (-1 if irreducible)
  bool is_inequivalent; // if atom is irreducible
  uint num_equivalents; // say how many they are (only for is_inequivalent)
  uint index_iatoms; // if calculated on the xstructure, the index within iatoms for the identical atoms
  // for order parameter                                 // order parameter
  int order_parameter_value; // order parameter
  bool order_parameter_atom; // order parameter
  // for partial occupation                              // partial occupation
  double partial_occupation_value; // partial occupation
  bool partial_occupation_flag; // partial occupation
  int shell; // neighbor shell number
  // for xOUTCAR
  aurostd::xvector<double> force; // force components from OUTCAR  //CO20211106
  // printing
  bool verbose; // verbose in printing
  bool print_RHT; // a printer for coord and name (general position)   //RHT
  bool print_cartesian; // print frac or cartesian
  // operators/functions                                 // operator/functions
  friend std::ostream& operator<<(std::ostream&, const _atom&); // print
  void CleanName(); // function to clean up the name
  void CleanSpin(); // function to clean up the spin from EZ vasp script
  void ClearSymmetry(); // clears symmetry //CO20190219

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  _atom deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "_atom"; }

private: //
  void free(); // free space
  void copy(const _atom& b); //

  // SERIALIZATION MEMBERS
#define JSON_atom_MEMBERS                                                                                                                                                                                        \
  fpos, cpos, corigin, coord, fpos_equation, cpos_equation, spin, spin_is_given, noncoll_spin, noncoll_spin_is_given, mass, type, name, name_is_given, cleanname, info, atomic_number, sd, ijk, isincell, basis, \
      reference, ireference, equivalent, is_inequivalent, num_equivalents, index_iatoms, order_parameter_value, order_parameter_atom, partial_occupation_value, partial_occupation_flag, shell, force, verbose,  \
      print_RHT, print_cartesian
};

class _atom_reference_cmp {
  // sorting through reference
public:
  bool operator()(const _atom& atom1, const _atom& atom2) const { return (bool) (atom1.reference < atom2.reference); }
};

class _atom_type_cmp {
  // sorting through type
public:
  bool operator()(const _atom& atom1, const _atom& atom2) const { return (bool) (atom1.type < atom2.type); }
};

void atoms_initialize();
uint GetAtomNumber(const std::string& symbol);
std::string GetAtomName(const std::string& symbol);
std::string GetAtomName(const uint& atnum);
std::string GetAtomSymbol(const std::string& symbol);
std::string GetAtomSymbol(const uint& atnum);
double GetAtomMass(const std::string& symbol, bool clean = true); // in Kg //CO20181129
double GetAtomMass(const uint& atnum); // in Kg
double GetAtomComptonCrossSection(const std::string& symbol); // barn (1 barn = 1e-28 m^2)
double GetAtomComptonCrossSection(const uint& atnum); // barn (1 barn = 1e-28 m^2)
double GetAtomPhotoelectricCrossSection(const std::string& symbol); // barn (1 barn = 1e-28 m^2)
double GetAtomPhotoelectricCrossSection(const uint& atnum); // barn (1 barn = 1e-28 m^2)
double GetAtomVolume(const std::string& symbol, bool clean = true); // CO20181129
double GetAtomVolume(const uint& atnum);
int GetAtomValenceIupac(const std::string& symbol);
int GetAtomValenceIupac(const uint& atnum);
int GetAtomValenceStd(const std::string& symbol);
int GetAtomValenceStd(const uint& atnum);
double GetAtomRadius(const std::string& symbol);
double GetAtomRadius(const uint& atnum);
double GetAtomRadiusCovalent(const std::string& symbol); // DX+CO20170904
double GetAtomRadiusCovalent(const uint& atnum); // DX+CO20170904
double GetAtomElectronegativity(const std::string& symbol);
double GetAtomElectronegativity(const uint& atnum);
std::string GetAtomCrystal(const std::string& symbol);
std::string GetAtomCrystal(const uint& atnum);
double GetAtomPettiforScale(const std::string& symbol);
double GetAtomPettiforScale(const uint& atnum);
bool GetAtomPettiforScale(const std::vector<std::string>& vsymbol, std::vector<double>& vvalue);
bool GetAtomPettiforScale(const std::vector<uint>& vatnum, std::vector<double>& vvalue);
bool GetAtomPettiforScale(const std::vector<std::string>& vsymbol, aurostd::xvector<double>& vvalue);
bool GetAtomPettiforScale(const std::vector<uint>& vatnum, aurostd::xvector<double>& vvalue);
bool SortAtomsPettiforScale(std::vector<std::string>& vsymbols, aurostd::xvector<int>& vorders, aurostd::xvector<double>& vvalues);
bool SortAtomsPettiforScale(std::vector<std::string>& vsymbols, std::vector<int>& vorders, std::vector<double>& vvalues);
bool SortAtomsPettiforScale(std::vector<std::string>& vsymbol, std::vector<int>& vorder);
bool SortAtomsPettiforScale(std::vector<std::string>& vsymbol, std::vector<double>& vvalue);
bool SortAtomsPettiforScale(std::vector<std::string>& vsymbol);
double GetPearsonCoefficient(const std::string&);
double GetPearsonCoefficient(const int&);
double GetAtomXrayScatt(const std::string& symbol);
double GetAtomXrayScatt(const uint& atnum);
std::vector<std::string> GetGroupOfAtoms(std::string& group_name); // DX20181220
double GetCompoundAttenuationLength(const std::vector<std::string>& species, const std::vector<double>& composition,
                                    const double& density); // density in g/cm^3, return in cm
double GetCompoundAttenuationLength(const std::deque<std::string>& _species, const std::deque<int>& _composition,
                                    const double& density); // density in g/cm^3, return in cm
// DX+CO START
// DX+CO END
//  routines of general use
std::string XATOM_AlphabetizationSpecies(const std::string& speciesA, const std::string& speciesB);
std::string XATOM_AlphabetizationSpecies(const std::vector<std::string>& vspecies);
std::string XATOM_AlphabetizationSpecies(const std::vector<std::string>& vspecies, const std::vector<double>& vnumbers);
void XATOM_AlphabetizationSpecies(std::string& system, std::vector<std::string>& vspecies, std::vector<double>& vnumbers);
void XATOM_AlphabetizationCompound(std::string& system, std::vector<std::string>& vspecies, std::vector<double>& vnumbers);
void XATOM_AlphabetizationSpecies(std::string& system, std::vector<std::string>& vspecies);
void XATOM_AlphabetizationSpecies(std::string& system);
void XATOM_AlphabetizationCompound(std::string& system);
uint XATOM_SplitAlloySpecies(const std::string& alloy_in, std::vector<std::string>& speciesX);
uint XATOM_SplitAlloySpecies(const std::string& alloy_in, std::vector<std::string>& speciesX, std::vector<double>& natomsX);
uint XATOM_SplitAlloyPseudoPotentials(const std::string& alloy_in, std::vector<std::string>& species_ppX);
uint XATOM_SplitAlloyPseudoPotentials(const std::string& alloy_in, std::vector<std::string>& species_ppX, std::vector<double>& natomsX);
// neighbor things
void GetUnitCellRep(const aurostd::xvector<double>& ppos, aurostd::xvector<double>& p_cell0, aurostd::xvector<int>& ijk, const aurostd::xmatrix<double>& lattice, const bool coord_flag);

aurostd::JSON::object atom2json(const _atom& atom, int coord_flag, int poccupation); // DX20170831 - atom2json

// --------------------------------------------------------------------------
class _sym_op : public JsonSerializable<_sym_op> {
public:
  // constructor destructor
  _sym_op(); // default, just allocate
  _sym_op(const _sym_op& b); // constructor copy
  ~_sym_op(); // kill everything
  // content
  // for _PGROUP_
  aurostd::xmatrix<double> Uc; // 3x3                        // uniques (not irreducible) operations on positions (Uc cartesian)
  aurostd::xmatrix<double> Uf; // 3x3                        // uniques (not irreducible) operations on indices   (Uf fractional)
  aurostd::xmatrix<double> generator; // 3x3                        // generator A, U=exp(A*theta)
  aurostd::xvector<double> generator_coefficients; // generator coefficients on Lx, Ly, Lz basis //DX20171206
  aurostd::xmatrix<aurostd::xcomplex<double>> SU2_matrix; // 2x2                 // SU(2) 2x2 complex matrix //DX20180115
  aurostd::xvector<aurostd::xcomplex<double>> su2_coefficients;
  // su(2) coefficients on sigma_1, sigma_2, sigma_3 basis (Pauli matrices) //DX20180115
  double angle; // angle axis
  aurostd::xvector<double> axis; // 3                          // (1,2,3)=axis
  aurostd::xvector<double> quaternion_vector; // GG
  aurostd::xmatrix<double> quaternion_matrix; // GG
  std::string str_type; // generic type of the operation
  std::string str_Hermann_Mauguin; // Hermann_Mauguin notation
  std::string str_Schoenflies; // Schoenflies notation
  bool flag_inversion; // flag if inversion
  bool is_pgroup; // bool is_pgroup
  // for _PGROUP_XTAL_
  bool is_pgroup_xtal; // bool is_pgroup_xtal
  // for _PGROUPK_PATTERSON_
  bool is_pgroupk_Patterson; // bool is_pgroup_Patterson //DX20200129
  // for _PGROUPK_
  bool is_pgroupk; // bool is_pgroupk
  // for _PGROUPK_XTAL_
  bool is_pgroupk_xtal; // bool is_pgroupk_xtal //DX20171205 - Added pgroupk_xtal
  // for _FGROUP_
  aurostd::xvector<double> ctau; // 3                          // translation in CARTESIAN       // FACTOR GROUP only, [0,1[
  aurostd::xvector<double> ftau; // 3                          // translation in FRACTIONAL      // FACTOR GROUP only, [0,1[
  std::vector<int> basis_atoms_map; // this is the vector that tell where the basis atom gets mapped by the operation
  std::vector<int> basis_types_map; // this is the vector that tell where the basis species gets mapped by the operation
  bool basis_map_calculated; // have we've calculated it?
  bool is_fgroup; // bool is_fgroup
  // for _SGROUP_
  aurostd::xvector<double> ctrasl;
  // 3                          // translation in CARTESIAN       // SPACE GROUP only, [integers]
  aurostd::xvector<double> ftrasl;
  // 3                          // translation in FRACTIONAL      // SPACE GROUP only, [ingegers]
  bool is_sgroup; // bool is_sgroup
  // operators
  // for _AGROUP_
  uint site; // uint site          // site index //DX20170803
  bool is_agroup; // bool is_agroup     // for site operation point group

  const _sym_op& operator=(const _sym_op& b);
  friend std::ostream& operator<<(std::ostream&, const _sym_op&);

  void setUc(const aurostd::xmatrix<double>& Uc, const aurostd::xmatrix<double>& lattice); // CO20190321
  void setUf(const aurostd::xmatrix<double>& Uf, const aurostd::xmatrix<double>& lattice); // CO20190321
  void setctau(const aurostd::xmatrix<double>& Uc, const aurostd::xmatrix<double>& lattice); // CO20190321
  void setftau(const aurostd::xmatrix<double>& Uf, const aurostd::xmatrix<double>& lattice); // CO20190321

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  _sym_op deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "_sym_op"; }

private:
  void free();
  void copy(const _sym_op& b);

  // SERIALIZATION MEMBERS
#define JSON_symop_MEMBERS                                                                                                                                                                               \
  Uc, Uf, generator, generator_coefficients, SU2_matrix, su2_coefficients, angle, axis, quaternion_vector, quaternion_matrix, str_type, str_Hermann_Mauguin, str_Schoenflies, flag_inversion, is_pgroup, \
      is_pgroup_xtal, is_pgroupk_Patterson, is_pgroupk, is_pgroupk_xtal, ctau, ftau, basis_atoms_map, basis_types_map, basis_map_calculated, is_fgroup, ctrasl, ftrasl, is_sgroup, site, is_agroup
};

// DX201801107 - add _kpoint class - START
//  --------------------------------------------------------------------------
class _kpoint {
public:
  // constructor destructor
  _kpoint(); // default, just allocate
  ~_kpoint(); // default, just allocate
  // content
  char iomode; // store format (not used yet)
  aurostd::xmatrix<double> klattice; // reciprocal lattice
  aurostd::xvector<double> fpos; // fractional position of kpoint
  aurostd::xvector<double> cpos; // Cartesian position of kpoint (not used yet)
  std::string label; // kpoint label (i.e., high-symmetry point labels)
  bool is_transformed; // indicates if kpoint is transformed from AFLOW standard
  const _kpoint& operator=(const _kpoint& b); // assignment operator
  // operators/functions                               // operator/functions
  [[nodiscard]] std::string str() const; // prints "fpos ! label" (e.g., 0.0000 0.0000 0.0000 ! \\Gamma)
  void TransformKpoint(const aurostd::xmatrix<double>& P);
  // transforms kpoint via P matrix (k'=k*P) and klattice via Q matrix (L_recip'=Q*L_recip) (see ITC-A pg. 79)
  friend std::ostream& operator<<(std::ostream&, const _kpoint&); // ostream operator
private:
  void free();
};

// DX201801107 - add _kpoint class - END

// --------------------------------------------------------------------------
class wyckoffsite_ITC : public JsonSerializable<wyckoffsite_ITC> {
  // Also for wyckoff sites
public:
  wyckoffsite_ITC();
  wyckoffsite_ITC(const wyckoffsite_ITC& b);
  ~wyckoffsite_ITC();
  // OPERATORS                                                  // --------------------------------------
  const wyckoffsite_ITC& operator=(const wyckoffsite_ITC& b); // some operators
  bool operator<(const wyckoffsite_ITC& b) const;
  // < operator //DX20190130 - so we can sort by Wyckoff letter, then by species
  friend std::ostream& operator<<(std::ostream&, const wyckoffsite_ITC&); // ostream
  // CONTENT
  aurostd::xvector<double> coord;
  uint index; // index //DX20200427
  std::string type; // chemical label etc //DX20200427
  std::string wyckoffSymbol;
  std::string letter; // DX20190128 - add Wyckoff letter
  std::string site_symmetry; // DX20190128 - add Wyckoff site symmetry
  uint multiplicity; // DX20190128 - add Wyckoff multiplicity
  double site_occupation; // DX20190128 - add Wyckoff site occupation
  std::vector<std::vector<std::string>> equations; // DX20190128 - add Wyckoff equations
  uint parameter_index; // DX20200513 - for ANRL parameter
  // initializers
  void getWyckoffFromLetter(uint space_group_number, // DX20200501
                            const std::string& Wyckoff_letter,
                            int setting = SG_SETTING_1);
  void getWyckoffFromLetter(const std::string& space_group_string, // DX20200501
                            const std::string& Wyckoff_letter);

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  wyckoffsite_ITC deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "wyckoffsite_ITC"; }

private: // ---------------------------------------
  void free(); // to free everything

  // SERIALIZATION MEMBERS
#define JSON_wyckoffsite_MEMBERS coord, index, type, wyckoffSymbol, letter, site_symmetry, multiplicity, site_occupation, equations, parameter_index
};

extern std::vector<std::string> vatom_symbol;             // store starting from ONE
extern std::vector<std::string> vatom_name;               // store starting from ONE
extern std::vector<double> vatom_mass;               // store starting from ONE
extern std::vector<double> vatom_volume;             // store starting from ONE
extern std::vector<int> vatom_valence_iupac;         // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
extern std::vector<int> vatom_valence_std;           // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
extern std::vector<double> vatom_miedema_phi_star;       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28 10.1016/0378-4363(80)90054-6
extern std::vector<double> vatom_miedema_nws;            // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28 10.1016/0378-4363(80)90054-6
extern std::vector<double> vatom_miedema_Vm;             // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28 10.1016/0378-4363(80)90054-6
extern std::vector<double> vatom_miedema_gamma_s;        // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28 10.1016/0378-4363(80)90054-6
extern std::vector<double> vatom_miedema_BVm;            // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28 10.1016/0378-4363(80)90054-6
extern std::vector<double> vatom_radius;             // store starting from ONE - Saxena
extern std::vector<double> vatom_radius_covalent;    // store starting from ONE - Codero, Covalent radii revisited, DOI: 10.1039/b801115j //DX+CO20170904
extern std::vector<double> vatom_electronegativity;  // store starting from ONE - Saxena
extern std::vector<std::string> vatom_crystal;            // store starting from ONE - Ashcroft Mermin
extern std::vector<double> vatom_xray_scatt;              // store starting from ONE
extern std::vector<double> vatom_pettifor_scale;              // store starting from ONE - Chemical Scale Pettifor Solid State Communications 51 31-34 1984
extern std::vector<double> vatom_pearson_coefficient;         // ME20181020

#endif // AFLOW_XATOM_H
