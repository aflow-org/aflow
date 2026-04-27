#ifndef AFLOW_XSTRUCTURE_H
#define AFLOW_XSTRUCTURE_H

#include <deque>
#include <iosfwd>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xplotter.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "structure/aflow_xatom.h"

class xstructure : public JsonSerializable<xstructure> {
public:
  // constructors/destructors                                   // --------------------------------------
  xstructure(const std::string& title = ""); // constructor default
  xstructure(const xstructure& b); // constructor copy
  xstructure(std::istream& input, int = IOVASP_POSCAR); // constructor from istream
  xstructure(std::ifstream& input, int = IOVASP_POSCAR); // constructor from ifstream
  xstructure(const std::stringstream& input, int = IOVASP_POSCAR);
  // constructor from stringstream //DX20210129 - added const
  xstructure(const std::string& input, int); // constructor from file
  xstructure(const std::string& url, const std::string& file, int = IOVASP_POSCAR); // constructor from URL
  ~xstructure(); // destructor
  // I/O, mutators                                              // --------------------------------------
  void initialize(const std::string& structure_title = "");
  // initialize xstructure based on input (avoids copying xstructure); //CO20211122
  void initialize(std::istream& input, int = IOVASP_POSCAR);
  // initialize xstructure based on input (avoids copying xstructure); //DX20210129
  void initialize(std::ifstream& input, int = IOVASP_POSCAR);
  // initialize xstructure based on input (avoids copying xstructure); //DX20210129
  void initialize(const std::stringstream& input, int = IOVASP_POSCAR);
  // initialize xstructure based on input (avoids copying xstructure); //DX20210129
  void initialize(const std::string& input, int);
  // initialize xstructure based on input (avoids copying xstructure); //CO20211122
  void initialize(const std::string& url, const std::string& file, int = IOVASP_POSCAR);
  // initialize xstructure based on input (avoids copying xstructure); //CO20211122
  bool GetStoich(); // get stoich_each_type - CO20170724
  bool sortAtomsEquivalent(); // sort by equivalent atoms - CO20190116
  bool FixLattices(); // Reciprocal/f2c/c2f
  void SetCoordinates(int mode); // change coordinates
  void MakeBasis(); // make basis for atoms (basis and number)
  void MakeTypes(); // refresh types based on num_each_type  //CO20180420
  void AddAtom(const _atom& atom,
               bool add_species,
               // adding an atom //CO20230319 - adding add_species which assumes the species information of the incoming atoms is correct, otherwise add by type
               bool check_present = true);
  void AddAtom(const std::deque<_atom>& atom,
               bool add_species,
               // adding a deque<_atom> //CO20210129 //DX20210201  //CO20230319 - adding add_species which assumes the species information of the incoming atoms is correct, otherwise add by type
               bool check_present = true);
  void AddAtom_POCC(const _atom& atom); // adding an atom FOR POCC ONLY
  void RemoveAtom(const uint& iat); // deleting an atom (index)
  void RemoveAtom(std::vector<uint>& v_atoms_to_remove); // deleting many atoms (indices)
  void RemoveAtom(); // removes all atoms //DX20210129
  void ReplaceAtoms(const std::deque<_atom>& new_atoms,
                    bool check_atom_overlap = true,
                    // replace all atoms SAFELY/CLEANLY //DX20210129 - added option to check atom overlap
                    bool sort_species = true); // CO20230319 - sort species
  void RemoveCopies(double = 1.0e-3); // deleting atoms too close F/C
  void RemoveFractionalCopies(double = 1.0e-3); // deleting atoms too close F
  void RemoveCartesianCopies(double = 1.0e-3); // deleting atoms too close C
  void AddCorners(); // for picturing purpose
  void clear(); // clear everything //DX20191220 - uppercase to lowercase clear
  void clean(); // performs stringstream clean //DX20191220 - uppercase to lowercase clean
  void ClearSpecies(); // Clear all the symmetry
  void CleanStructure(); // Fix up structure - ME20211004
  void ShiftOriginToAtom(const int& iat); // Shift the origin to atom(iat)
  void IdenticalAtoms(); // Make identical atoms
  void SwapCoordinates(const uint& i, const uint& j); // Permute Coordinates i with j
  std::string SpeciesLabel(const uint& A); // Returns the Label of the specie A (if available)
  void SpeciesSwap(const uint& A, const uint& B); // Permute Species A with B (safe for species C).
  bool SpeciesGetAlphabetic(); // Check is species are in alphabetic order
  bool SpeciesPutAlphabetic(); // Put Species in alphabetic
  std::string SpeciesString(); // Gives a string with the list of all the species
  uint SetSpecies(const std::deque<std::string>& vspecies, // Set the species
                  bool sort_species = true); // CO20230319 - sort species
  void UpdateSpecies(const _atom& atom, bool add_species);
  void GetLatticeType(double sym_eps = AUROSTD_MAX_DOUBLE,
                      bool no_scan = false); // Get all lattices
  void GetLatticeType(xstructure& sp, xstructure& sc, double sym_eps = AUROSTD_MAX_DOUBLE,
                      bool no_scan = false); // Get all lattices
  void GetExtendedCrystallographicData(double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = SG_SETTING_1);
  void GetExtendedCrystallographicData(xstructure& sp, xstructure& sc, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = SG_SETTING_1);
  void GetRealLatticeType(xstructure& sp, xstructure& sc,
                          double sym_eps = AUROSTD_MAX_DOUBLE); // Get real lattice type //DX2021011
  void GetRealLatticeType(double sym_eps = AUROSTD_MAX_DOUBLE); // Get real lattice type //DX20210211
  void GetReciprocalLatticeType(xstructure& sp, xstructure& sc,
                                double sym_eps = AUROSTD_MAX_DOUBLE); // Get reciprocal lattice type //DX20210209
  void GetReciprocalLatticeType(double sym_eps = AUROSTD_MAX_DOUBLE); // Get reciprocal lattice type //DX20210209
  void GetSuperlatticeType(xstructure& sp, xstructure& sc,
                           double sym_eps = AUROSTD_MAX_DOUBLE); // Get superlattice type //DX20210209
  void GetSuperlatticeType(double sym_eps = AUROSTD_MAX_DOUBLE); // Get superlattice type //DX20210209
  void Standard_Primitive_UnitCellForm(); // Reduce the Unit Cell to Standard Primitive Form
  void GetStandardPrimitive(); // stub for void Standard_Primitive_UnitCellForm(void);
  void Standard_Conventional_UnitCellForm(); // Reduce the Unit Cell to Standard Conventional Form
  void GetStandardConventional(); // stub for void Standard_Conventional_UnitCellForm(void);
  void NiggliUnitCellForm(); // Reduce the Unit Cell to Niggli Form
  void GetNiggliStructures(std::vector<xstructure>& structures, uint start_index = 0,
                           uint end_index = AUROSTD_MAX_UINT); // DX20201006
  void MinkowskiBasisReduction(); // Reduce the Basis to the max orthogonality (Minkowski)
  void GetMinkowskiStructures(std::vector<xstructure>& structures, uint start_index = 0,
                              uint end_index = AUROSTD_MAX_UINT); // DX20201006
  void LatticeReduction(); // Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
  void BringInCell(double tolerance = _ZERO_TOL_, double upper_bound = 1.0, double lower_bound = 0.0); // DX20190904
  void BringInCompact(); // Bring all the atoms near the origin
  void BringInWignerSeitz(); // Bring all the atoms in the Wigner Seitz Cell
  void GetPrimitive_20210322(double eps = AUROSTD_MAX_DOUBLE); // Make it primitive, if possible //DX20210323
  void GetPrimitive(); // Make it primitive, if possible
  void GetPrimitive(double tol); // Make it primitive, if possible
  void GetPrimitive1(); // Make it primitive, if possible
  void GetPrimitive2(); // Make it primitive, if possible
  void GetPrimitive3(); // Make it primitive, if possible
  void GetPrimitiveStructures(std::vector<xstructure>& structures, uint start_index = 0,
                              uint end_index = AUROSTD_MAX_UINT); // DX20201006
  uint GetPrimitiveCell(); // Make it primitive, if possible. Returns 1 if routine fails (RHT)   //RHT
  double MinDist(); // get minimum interatomic distance -- CO20171024
  void ReScale(const double& in_scale); // Change scale but keep volume fixed
  void SetScale(const double& in_scale); // Change scale
  void UpdateCartesianCoordinates(); // AS20200514
  void ChangeBasis(const aurostd::xmatrix<double>& transformation_matrix); // DX20201215
  void Rotate(const aurostd::xmatrix<double>& rm); // DX20201215 - added modify-in-place variant
  void TransformStructure(const aurostd::xmatrix<double>& transformation_matrix, const aurostd::xmatrix<double>& rotation);
  void TransformStructure(const aurostd::xmatrix<double>& transformation_matrix, const aurostd::xmatrix<double>& rotation, const aurostd::xvector<double>& origin_shift, bool is_shift_frac = false);
  void ShiftPos(const aurostd::xvector<double>& shift, bool is_frac);
  // Shift origin by vector (Cartesian/fractional boolean) //DX20201215 - added modify-in-place variant
  void ShiftCPos(const aurostd::xvector<double>& shift);
  // Shift origin by Cartesian vector //DX20201215 - added modify-in-place variant
  void ShiftFPos(const aurostd::xvector<double>& shift);
  // Shift origin by fractional vector //DX20201215 - added modify-in-place variant
  void SetVolume(const double& in_volume); // Change volume
  void SetAutoVolume(bool use_AFLOW_defaults_in = false); // Change volume to sum of atoms  //CO20191010
  void InflateLattice(const double& coefficient); // Inflate lattice
  void InflateVolume(const double& coefficient); // Inflate volume
  void foldAtomsInCell( // fold atoms into new cell representation //DX20210113
      const aurostd::xmatrix<double>& lattice_new,
      bool skew,
      double tol,
      bool check_min_dists = true);
  std::string platon2print(bool, bool, double, double, double, double); // Create Platon input file >=51108
  void DecorateWithElements(); // Decorate with elements (alphabetic order) - useful for platon
  void DecorateWithFakeElements(); // Decorate with fake elements - useful for prototypes //DX20200727
  [[nodiscard]] std::vector<std::string> GetElements(bool clean_name = false,
                                                     bool fake_names = false) const; // DX20200724 //SD20220222 - made function const
  [[nodiscard]] std::vector<std::string> GetElementsFromAtomNames(bool clean_name = true) const; // DX20200724 //SD20220222 - made function const
  std::vector<uint> GetReducedComposition(bool numerical_sort = false); // DX20200724
  std::string platon2sg(bool P_EQUAL = DEFAULT_PLATON_P_EQUAL, bool P_EXACT = DEFAULT_PLATON_P_EXACT, double P_ang = DEFAULT_PLATON_P_ANG, double P_d1 = DEFAULT_PLATON_P_D1, double P_d2 = DEFAULT_PLATON_P_D2, double P_d3 = DEFAULT_PLATON_P_D3);
  std::string findsym2sg(double tolerance = DEFAULT_FINDSYM_TOL);
  std::string findsym2execute(double tolerance = DEFAULT_FINDSYM_TOL);
  std::string findsym2print(double tolerance = DEFAULT_FINDSYM_TOL);
  //  string platon2sg(void);
  [[nodiscard]] double GetVolume() const; // Return volume  //CO20200201
  [[nodiscard]] double Volume() const; // Return volume  //CO20200201
  double GetZVAL(const std::vector<double>& vZVAL); // Given the ZVAL of each species, it returns total ZVAL of cell
  double GetPOMASS(const std::vector<double>& vPOMASS);
  // Given the POMASS of each species, it returns total POMASS of cell
  void ClearSymmetry(); // Clear all the symmetry
  bool CalculateSymmetry(bool, double); // Calculate the symmetry
  bool CalculateSymmetry(); // Calculate the symmetry
  void CalculateSymmetryPointGroup(bool); // Calculate the symmetry
  void CalculateSymmetryPointGroup(); // Calculate the symmetry
  void CalculateSymmetryFactorGroup(bool); // Calculate the symmetry
  void CalculateSymmetryFactorGroup(); // Calculate the symmetry
  void CalculateSymmetryPointGroupCrystal(bool); // Calculate the symmetry
  void CalculateSymmetryPointGroupCrystal(); // Calculate the symmetry
  void CalculateSymmetryPointGroupKLattice(bool); // Calculate the symmetry
  void CalculateSymmetryPointGroupKLattice(); // Calculate the symmetry
  void CalculateSymmetryPointGroupKCrystal(bool); // Calculate the symmetry  //ME20200114
  void CalculateSymmetryPointGroupKCrystal(); // Calculate the symmetry  //ME20200114
  void CalculateSymmetryPointGroupKPatterson(bool); // Calculate the symmetry  //ME20200129
  void CalculateSymmetryPointGroupKPatterson(); // Calculate the symmetry  //ME20200129
  int GenerateGridAtoms(double); // generate grid of atoms
  int GenerateGridAtoms(int); // generate grid of atoms
  int GenerateGridAtoms(int, int, int); // generate grid of atoms
  int GenerateGridAtoms(const aurostd::xvector<int>& dims); // generate grid of atoms
  int GenerateGridAtoms(int, int, int, int, int, int); // generate grid of atoms
  int GenerateLIJK(double); // generate lijk look up table
  // QUANTUM ESPRESSO AND ABINIT AND AIMS                       // --------------------------------------
  void fixEmptyAtomNames(bool force_fix = false); // CO20200829
  void buildGenericTitle(bool vasp_input = false, bool force_fix = false); // build a nice title with atoms
  void xstructure2qe(); // some wrap up IOs to convert format to QE
  void xstructure2vasp(); // some wrap up IOs to convert format to VASP
  void xstructure2itc(); // some wrap up IOs to convert format to ITC  //CO20220613
  void xstructure2abinit(); // some wrap up IOs to convert format to ABINIT
  void xstructure2aims(); // some wrap up IOs to convert format to AIMS
  void xstructure2cif(); // some wrap up IOs to convert format to CIF //DX20190123
  void xstructure2abccar(); // some wrap up IOs to convert format to ABCCAR //DX20190123
  void xstructure2elk(); // some wrap up IOs to convert format to ELK //DX20200313
  void xstructure2atat(); // some wrap up IOs to convert format to ATAT //SD20220123
  void xstructure2lmp(); // some wrap up IOs to convert format to LAMMPS //SD20240111

  [[nodiscard]] aurostd::x3DWriter render() const;
  //[CO20180420 - moved outside of xstructure]bool sortAtomsTypes(const _atom& a1,const _atom& a2);		// sort atoms by types
  //[CO20180420 - moved outside of xstructure]bool sortAtomsNames(const _atom& a1,const _atom& a2);		// sort atoms by names
  // OPERATORS                                                  // --------------------------------------
  const xstructure& operator=(const xstructure& b); // some operators
  friend std::istream& operator>>(std::istream&, xstructure&); // istream
  friend std::ostream& operator<<(std::ostream&, const xstructure&); // ostream
  // CONTENT                                                    // --------------------------------------
  std::string title; // Title of the structure
  std::string directory; // Directory where xstructure came from //DX
  std::string prototype; // Prototype of the structure
  std::string info; // Info of the structure
  int iomode; // IOVASP_POSCAR/IOXXXX
  // int num_types=num_each_type.size();                        // old useless stuff
  // int num_atoms=atoms.size();                                // old useless stuff
  bool neg_scale; // flag for negative scale (for printing)
  double scale; // scale (always linear A)
  bool neg_scale_second; // POCC (hnf vs. tol) //CO20180409
  double scale_second; // POCC hnf/stoich tol/site tol //CO20180409
  aurostd::xoption scale_third;
  // if there is a third scale number provided, we use isentry and content_double //CO20170803 - site tol
  char coord_type[2]; // type of coordinates
  bool coord_flag; // _COORDS_FRACTIONAL_ (0) fractional, _COORDS_CARTESIAN_ (1) cartesian.
  bool isd; // true=Selective dynamics; false=no selective dynamics.
  aurostd::xmatrix<double> lattice;
  // LATTICE in REAL SPACE (meters)            // vector per RAW (must trasp per algebra)
  double a, b, c, alpha, beta, gamma; // LATTICE in a,b,c,alpha,beta,gamma
  aurostd::xmatrix<double> klattice;
  // LATTICE in MOMENTUM SPACE (1/meters)      // vevror per RAW (must trasp per algebra)
  aurostd::xvector<double> origin; // origin
  aurostd::xmatrix<double> f2c; // transformation matrix for F vector per COLUM f2c=trasp(lattice)
  aurostd::xmatrix<double> c2f; // transformation matrix for C vector per ROW   c2f=inv(trasp(lattice))
  double equiv_fpos_epsilon; // when they are the same DEFAULT _EQUIV_FPOS_EPS_
  std::deque<int> num_each_type; // WARNING: we use starting from 0
  std::deque<double> comp_each_type; // WARNING: we use starting from 0
  std::deque<double> stoich_each_type; // WARNING: we use starting from 0 - 20170724
  std::deque<_atom> atoms; // WARNING: we use starting from 0
  std::deque<std::string> species, species_pp, species_pp_type, species_pp_version;
  // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<double> species_pp_ZVAL; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<std::deque<double>> species_pp_vLDAU; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<double> species_volume; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  std::deque<double> species_mass; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
  //  ----------------------------------------------------------------------------------------
  // SYMBOLIC MATH stuff
  bool symbolic_math_representation_only; // print symbolic math representation only //DX20180618
  bool constrained_symmetry_calculation;
  // append symbolic math representation for constrained symmetry calculation //DX20180618
  std::vector<std::vector<std::string>> symbolic_math_lattice; // symbolic math representation of lattice //DX20180618
  uint num_parameters; // number of parameters ANRL 20180618
  uint num_lattice_parameters; // number of lattice parameters ANRL 20180618
  std::vector<std::string> prototype_parameter_list; // prototype parameter list ANRL 20180618
  std::vector<double> prototype_parameter_values; // prototype parameters values ANRL 20180618
  //  ----------------------------------------------------------------------------------------
  bool is_vasp4_poscar_format; // flags for VASP4*
  bool is_vasp5_poscar_format; // flags for VASP5*
  bool primitive_calculated; // flags for calculation //DX20201007
  bool Niggli_calculated; // flags for calculation
  bool Niggli_avoid; // flags for avoiding the calculation
  bool Minkowski_calculated; // flags for calculation
  bool Minkowski_avoid; // flags for avoiding the calculation
  bool LatticeReduction_calculated; // flags for calculation
  bool LatticeReduction_avoid; // flags for avoiding the calculation
  //  ----------------------------------------------------------------------------------------
  // PRINTING stuff
  std::string PrintSymbolicMathRepresentation(); // Print symbolic math representation of structure //DX20180618
  std::string PrintUNCLE(); // Print in UNCLE format
  //  ----------------------------------------------------------------------------------------
  // LATTICE stuff
  bool Standard_Lattice_calculated; // flags for calculation
  bool Standard_Lattice_avoid; // flags for avoiding the calculation
  bool Standard_Lattice_primitive; // flags for calculation
  bool Standard_Lattice_conventional; // flags for calculation
  bool Standard_Lattice_has_failed; // flags for Lattice has failed ?
  std::string bravais_lattice_type; // lattice type as a string  (14)
  std::string bravais_lattice_variation_type; // lattice type as a string WSETYAWAN mod  (with the mods of WSETYAWAN)
  std::string bravais_lattice_system; // lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
  std::string bravais_lattice_lattice_type; // lattice_lattice type as a string  (14)
  std::string bravais_lattice_lattice_variation_type;
  // lattice_lattice type as a string WSETYAWAN mod  (with the mods of WSETYAWAN)
  std::string bravais_lattice_lattice_system; // lattice_lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
  std::string pearson_symbol; // pearson symbol as a string
  std::string reciprocal_lattice_type; // reciprocal lattice type as a string
  std::string reciprocal_lattice_variation_type; // reciprocal lattice type as a string WSETYAWAN mod
  // string reciprocal_conventional_lattice_type;                // reciprocal lattice type as a string
  aurostd::xmatrix<double> bravais_superlattice_lattice; // superlattice lattice (identical atoms) //DX20210209
  std::string bravais_superlattice_type; // super lattice type as a string (identical atoms)
  std::string bravais_superlattice_variation_type; // super lattice type as a string (identical atoms) WSETYAWAN mod
  std::string bravais_superlattice_system; // lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
  std::string pearson_symbol_superlattice; // pearson symbol of the superlattice (identical atoms)
  bool volume_changed_original2new;
  // flag for volume has changed between original and new (i.e., transformation won't work) //DX20181105
  aurostd::xmatrix<double> transform_coordinates_original2new;
  // transform coordinate system from original to new; (Q in ITC notation) //DX20181105
  aurostd::xmatrix<double> transform_coordinates_new2original;
  // transform coordinate system from new to original; (Q^-1 in ITC notation) //DX20181105
  aurostd::xmatrix<double> rotate_lattice_original2new;
  // rotate from original to new lattice; (P in ITC notation) //DX20181105
  aurostd::xmatrix<double> rotate_lattice_new2original;
  // rotate from new to original lattice; (P^-1 in ITC notation) //DX20181105
  //  ----------------------------------------------------------------------------------------
  // GENERAL PURPOSE LABELS                                     // general purpose label
  uint label_uint; // general purpose label_uint
  int label_int; // general purpose label_int
  double label_double; // general purpose label_double
  // ----------------------------------------------------------------------------------------
  // ORDER PARAMETER                                            // order parameter for xstructure
  bool order_parameter_structure; // order parameter for xstructure
  std::vector<uint> order_parameter_atoms; // indices of atoms to be shuffled
  uint order_parameter_orbit; // number of equivalent configurations with the factor group
  int order_parameter_sum; // sum of all the order parameters
  // ----------------------------------------------------------------------------------------
  // PARTIAL OCCUPATION                                         // partial occupation for xstructure
  bool partial_occupation_flag; // flags for partial occupation true/false
  double partial_occupation_site_tol; // tolerance for partial occupation site >=0.0 <=1.0   //CO20180409
  double partial_occupation_stoich_tol; // tolerance for partial occupation stoich >=0.0 <=1.0 //CO20180409
  int partial_occupation_HNF; // volume HNF size
  std::vector<int> partial_occupation_sublattice;
  // contains the information about the specie# of the sublattice in the partial occupation otherwise _pocc_no_sublattice_
  // ----------------------------------------------------------------------------------------
  // GEOMETRY ENERGETICS after the QM calculations              // --------------------------------------
  void qm_clear(); // QM create/clean all the vectors
  void qm_recycle(); // QM shift data from QM to GEOM
  void qm_load(const std::string& directory, const std::string& suffix = "", int = IOVASP_POSCAR);
  // QM results load from an ab-initio calculation
  bool qm_calculated; // QM calculation
  double qm_scale; // QM scale (always linear A)
  aurostd::xmatrix<double> qm_lattice; // QM LATTICE in REAL SPACE (meters)
  aurostd::xmatrix<double> qm_klattice; // QM LATTICE in MOMENTUM SPACE (1/meters)     // SAVED TRASP
  aurostd::xvector<double> qm_origin; // QM origin
  aurostd::xmatrix<double> qm_f2c; // QM transformation matrix for F vector per COLUM f2c=trasp(lattice)
  aurostd::xmatrix<double> qm_c2f; // QM transformation matrix for C vector per ROW   c2f=inv(trasp(lattice))
  std::deque<_atom> qm_atoms; // QM WARNING: we use starting from 0
  std::vector<aurostd::xvector<double>> qm_forces; // QM FORCES calculation
  bool qm_forces_write; // QM FORCES calculation
  std::vector<aurostd::xvector<double>> qm_positions; // QM POSITIONS calculation
  bool qm_positions_write; // QM POSITIONS calculation
  double qm_E_cell, qm_dE_cell, qm_H_cell, qm_PV_cell, qm_P, qm_mag_cell; // QM energetics PER CELL
  double qm_E_atom, qm_dE_atom, qm_H_atom, qm_PV_atom, qm_mag_atom; // QM energetics ATOMIC
  // ----------------------------------------------------------------------------------------
  // KPOINTS                                                    // --------------------------------------
  int kpoints_mode; // mode of kpoints
  int kpoints_k1, kpoints_k2, kpoints_k3; // parameters that are plug during
  double kpoints_s1, kpoints_s2, kpoints_s3; // parameters that are plug during
  int kpoints_kmax, kpoints_kppra; // load/unload and calculations
  std::string kpoints_kscheme; // of ab-initio
  // ---------------------- SYMMETRY --------------------------------------------------------
  // A=U*B but in A and B we plug vectors as columns watch lattice is per row Uc=A*inv(B)
  // A is the lattice (vectors per colum), B is the test lattice (epr column)
  // Uc is the point_group operation which operates AFTER the vector (row)
  // as: new_vector_row=old_vector_row*Uc  and point_group is the list of all the Uc !!!
  // DX+CO START
  double dist_nn_min;
  // SYMMETRY TOLERANCE ----------------------------
  bool sym_eps_calculated; // was it calculated automatically per symmetry operations (aflowSYM)?
  double sym_eps; // universal tolerance for symmetry (dictates resolution and mapping tolerances)
  uint sym_eps_change_count; // universal tolerance count for symmetry //DX20180223 - added count to xstructure
  bool sym_eps_no_scan; // do not use tolerance scan (forced by user or because the scan terminated) //DX20210331
  // DX+CO END
  //  POINT GROUP                                                // POINT GROUP LATTICE
  std::vector<_sym_op> pgroup; // rotations/inversions operations
  bool pgroup_calculated; // WARNING: we use starting from 0
  // POINT GROUP CRYSTAL                                        // POINT GROUP CRYSTAL
  std::vector<_sym_op> pgroup_xtal; // rotations/inversions operations
  bool pgroup_xtal_calculated; // WARNING: we use starting from 0
  std::string crystal_family; // crystal and point group properties
  std::string crystal_system; // crystal and point group properties
  std::string point_group_crystal_class; // crystal and point group properties
  std::string point_group_Shoenflies; // crystal and point group properties
  std::string point_group_Hermann_Mauguin; // crystal and point group properties
  std::string point_group_orbifold; // crystal and point group properties
  std::string point_group_type; // crystal and point group properties
  std::string point_group_order; // crystal and point group properties
  std::string point_group_structure; // crystal and point group properties
  // POINT GROUP PATTERSON                                      // POINT GROUP PATTERSON //DX20200129
  std::vector<_sym_op> pgroupk_Patterson; // rotations/inversions operations
  bool pgroupk_Patterson_calculated; // WARNING: we use starting from 0
  // POINT GROUP KLATTICE                                       // POINT GROUP
  std::vector<_sym_op> pgroupk; // rotations/inversions operations
  bool pgroupk_calculated; // WARNING: we use starting from 0
  // POINT GROUP KCRYSTAL                                       // POINT GROUP
  std::vector<_sym_op> pgroupk_xtal; // rotations/inversions operations
  bool pgroupk_xtal_calculated; // WARNING: we use starting from 0
  // FACTOR GROUP                                               // FACTOR GROUP
  std::vector<_sym_op> fgroup; // rotations/inversions + incell_translations operations
  bool fgroup_calculated; // WARNING: we use starting from 0
  // SPACE GROUP                                                // SPACE GROUP
  std::vector<_sym_op> sgroup; // rotations/inversions + incell        //outcell_translations operations
  bool sgroup_calculated; // WARNING: we use starting from 0
  double sgroup_radius; // radius of application (all ops connecting objects inside the sphere)
  aurostd::xvector<int> sgroup_radius_dims; // dimension of the radius (in +- integers)
  // SITE POINT GROUP                                           // SITE POINT GROUP
  bool agroup_calculated; // WARNING: we use starting from 0
  std::vector<std::vector<_sym_op>> agroup; // rotations/inversions operations on sites, has one for each atom (all)
  // INEQUIVALENTE ATOMS                                        // --------------------------------------
  bool iatoms_calculated; // given the symmetry, the atoms are mapped in inequivalent
  std::vector<std::vector<int>> iatoms; // equivalent/inequivalent atoms lookup table
  // SPACE GROUP CALCULATION WITH PLATON/FINDSYM                // with platon >= 51108
  std::string spacegroup; // space group symbol
  std::string spacegrouplabel; // the number with #
  std::string spacegroupoption; // origin, axes and so on
  int spacegroupnumber; // the number
  int spacegroupnumberoption; // the option as number
  bool is_spacegroup_platon, is_spacegroup_findsym, is_spacegroup_aflow; // got spacegroup
  // SPACE GROUP ITC
  // DX+CO START
  std::string crystal_system_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  std::string point_group_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  char bravais_label_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  char lattice_label_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  uint space_group_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  std::string wyckoff_library_entry_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  int setting_ITC; // aflow_symmetry_spacegroup.cpp (RHT) //DX20170830 - SGDATA
  aurostd::xvector<double> origin_ITC; // aflow_symmetry_spacegroup.cpp (RHT) //DX20170830 - SGDATA
  std::vector<std::string> general_position_ITC; // aflow_symmetry_spacegroup.cpp (RHT) //DX20170830 - SGDATA
  // DX+CO END
  //  double volume;  USE double GetVolume(const xstructure& a);  and SetVolume //
  //   int number_of_atoms; USE (int) atoms.size();  looks you do not need.
  // DX+CO START
  std::vector<std::string> wyccar_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  aurostd::xmatrix<double> standard_lattice_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  std::deque<_atom> standard_basis_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  std::vector<wyckoffsite_ITC> wyckoff_sites_ITC; // aflow_symmetry_spacegroup.cpp (RHT) //(x,y,z) XX
  std::vector<std::string> wyckoff_symbols_ITC; // aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(); // aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol); // aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol, bool& no_scan); // aflow_symmetry_spacegroup.cpp (RHT)
  // uint SpaceGroup_ITC(double& use_tol, const int& manual_it);// aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol, const int& setting); // aflow_symmetry_spacegroup.cpp (RHT) //DX20180806
  // uint SpaceGroup_ITC(double& use_tol,const int& manual_it,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT)
  uint SpaceGroup_ITC(double& use_tol, const int& setting, bool& no_scan);
  // aflow_symmetry_spacegroup.cpp (RHT) //DX20180806
  uint SpaceGroup_ITC(double& use_tol, const int& manual_it, const int& setting, bool& no_scan);
  // aflow_symmetry_spacegroup.cpp (RHT) //DX20180806
  std::string aflow2sg(); // aflow_symmetry_spacegroup.cpp (DX)
  std::string aflow2sg(double& use_tol); // aflow_symmetry_spacegroup.cpp (DX)
  std::string aflow2sg(double& use_tol, const int& manual_it); // aflow_symmetry_spacegroup.cpp (DX)
  // DX+CO END
  //  ---------------------- FROZSL ----------------------
  std::ostringstream FROZSL_output(std::vector<std::string> Kvectors);
  // ---------------------- PHONONS ---------------------------------------------------------
  // based on the Maradudin harmonic and deformation analysis
  // LIJK OBEJCTS                                               // LIJK OBEJCTS   WORKING
  bool lijk_calculated; // is calculated ?
  std::vector<aurostd::xvector<int>> lijk_table; // bravais lattice look up table l <-> i,j,k
  std::vector<aurostd::xvector<double>> lijk_fpos; // bravais lattice look up table fpos same as lijk_table
  std::vector<aurostd::xvector<double>> lijk_cpos; // bravais lattice look up table cpos
  aurostd::xvector<int> lijk_dims; // dimension
  // ----------------------------------------------------------------------------------------
  // GRID ATOMS                                                 // --------------------------------------
  bool grid_atoms_calculated; // GRID ATOMS from dimsL to dimsH has been calculated ?
  aurostd::xvector<int> grid_atoms_dimsL; // dims low of i,j,k (minus is NOT included)
  aurostd::xvector<int> grid_atoms_dimsH; // dims high of i,j,k (plus is NOT included)
  std::deque<_atom> grid_atoms; // WARNING: we use starting from 0
  int grid_atoms_number; // how many...
  std::vector<int> grid_atoms_sc2pcMap; // CO20170804 - mapping between grid_atoms (sc) and atoms (pc)
  std::vector<int> grid_atoms_pc2scMap; // CO20170804 - mapping between grid_atoms (sc) and atoms (pc)
  // ----------------------------------------------------------------------------------------
  // NEIGHBORS OBEJCTS EXPERIMENTAL/UNFINISHED                 // NEIGHBORS OBEJCTS   WORKING EXPERIMENTAL
  // aurostd::xvector<int> ndims;                                        // dimension of the radius (in +- integers)
  // std::deque<_atom> ashell;                                 // all the atoms in the shell
  // std::deque<deque<_atom> > natoms;                        // vector of vectors
  // std::vector<vector<double> > rshell;                       // vector of shells
  // std::vector<vector<int> > nshell;                          // vector of density in shells
  // int nbins;                                                 // number of bins
  void checkStructure(); // RF20200831; rescale structure to 1 and check whether e.g. species and atoms are present
  //
  void GetNeighbors(std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, double rmin = 0.0, bool prim = true, bool unique_only = true); // CO20200912
  void GetNeighbors(std::deque<_atom>& atoms_cell,
                    std::deque<std::deque<uint>>& i_neighbors,
                    std::deque<std::deque<double>>& distances,
                    double rmin = 0.0,
                    bool prim = true,
                    bool unique_only = true); // CO20200912
  void GetNeighbors(std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, double rmax, double rmin = 0.0, bool prim = true, bool unique_only = true); // CO20200912
  void GetNeighbors(std::deque<_atom>& atoms_cell,
                    std::deque<std::deque<uint>>& i_neighbors,
                    std::deque<std::deque<double>>& distances,
                    double rmax,
                    double rmin = 0.0,
                    bool prim = true,
                    bool unique_only = true); // CO20200912
  //
  void GetCoordinations(std::deque<std::deque<uint>>& coordinations, double rmin = 0.0, double tol = 0.5, bool prim = true, bool unique_only = true); // CO20200912
  void GetCoordinations(std::deque<_atom>& atoms_cell, std::deque<std::deque<uint>>& coordinations, double rmin = 0.0, double tol = 0.5, bool prim = true, bool unique_only = true); // CO20200912
  void GetCoordinations(std::deque<std::deque<uint>>& coordinations, double rmax, double rmin = 0.0, double tol = 0.5, bool prim = true, bool unique_only = true); // CO20200912
  void GetCoordinations(std::deque<_atom>& atoms_cell, std::deque<std::deque<uint>>& coordinations, double rmax, double rmin = 0.0, double tol = 0.5, bool prim = true, bool unique_only = true); // CO20200912
  //
  void GetCoordinations(std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, std::deque<std::deque<uint>>& coordinations, double rmin = 0.0, double tol = 0.5, bool prim = true, bool unique_only = true); // CO20200912
  void GetCoordinations(std::deque<_atom>& atoms_cell,
                        std::deque<std::deque<uint>>& i_neighbors,
                        std::deque<std::deque<double>>& distances,
                        std::deque<std::deque<uint>>& coordinations,
                        double rmin = 0.0,
                        double tol = 0.5,
                        bool prim = true,
                        bool unique_only = true); // CO20200912
  void GetCoordinations(
      std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, std::deque<std::deque<uint>>& coordinations, double rmax, double rmin = 0.0, double tol = 0.5, bool prim = true, bool unique_only = true); // CO20200912
  void GetCoordinations(std::deque<_atom>& atoms_cell,
                        std::deque<std::deque<uint>>& i_neighbors,
                        std::deque<std::deque<double>>& distances,
                        std::deque<std::deque<uint>>& coordinations,
                        double rmax,
                        double rmin = 0.0,
                        double tol = 0.5,
                        bool prim = true,
                        bool unique_only = true); // CO20200912
  //
  // NEIGHBORS OBEJCTS OLD-ACONVASP BUT WORKS                  // NEIGHBORS OBEJCTS
  // GetNeighData collects all the neighbor data between rmin and rmax and stores it for each atom in a vector of atom objects in order of increasing distance.
  void GetNeighData(const double rmax, std::deque<std::deque<_atom>>& neigh_mat, const double rmin = _ZERO_TOL_) const;
  // CO20220623 - use this function
  void GetNeighData(std::deque<_atom>& atoms_cell, const double rmax, std::deque<std::deque<_atom>>& neigh_mat,
                    const double rmin = _ZERO_TOL_) const; // CO20220623 - use this function
  void GetNeighData_20220101(const std::deque<_atom>& in_atom_vec, const double rmin, const double rmax, std::deque<std::deque<_atom>>& neigh_mat) const;
  // CO20220623 - AVOID this function, use one above

  // ----------------------------------------------------------------------------------------
  // OUTPUT/ERROR FLAGS                                         // --------------------------------------
  bool Niggli_has_failed; // Niggli has failed ?
  bool Minkowski_has_failed; // Minkowski has failed ?
  bool LatticeReduction_has_failed; // LatticeReduction has failed ?
  bool write_lattice_flag; // flag for OUTPUT printing
  bool write_klattice_flag; // flag for OUTPUT printing
  bool write_inequivalent_flag; // flag for OUTPUT printing
  bool write_DEBUG_flag; // flag for OUTPUT printing
  bool error_flag; // flag true is error
  std::string error_string; // contains type of error
  // END OF CONTENT                                             //

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  xstructure deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "xstructure"; }

private: // ---------------------------------------
  void free(); // to free everything
  void copy(const xstructure& b); // the flag is necessary because sometimes you need to allocate the space.

  // SERIALIZATION MEMBERS ignoring following vars: scale_third, coord_type, \*/
#define JSON_xstructure_MEMBERS                                                                                                                                                                                      \
  title, directory, prototype, info, iomode, neg_scale, scale, neg_scale_second, scale_second, coord_flag, isd, lattice, a, b, c, alpha, beta, gamma, klattice, origin, f2c, c2f, equiv_fpos_epsilon, num_each_type, \
      comp_each_type, stoich_each_type, atoms, species, species_pp, species_pp_type, species_pp_version, species_pp_ZVAL, species_pp_vLDAU, species_volume, species_mass, symbolic_math_representation_only,         \
      constrained_symmetry_calculation, symbolic_math_lattice, num_parameters, num_lattice_parameters, prototype_parameter_list, prototype_parameter_values, is_vasp4_poscar_format, is_vasp5_poscar_format,         \
      primitive_calculated, Niggli_calculated, Niggli_avoid, Minkowski_calculated, Minkowski_avoid, LatticeReduction_calculated, LatticeReduction_avoid, Standard_Lattice_calculated, Standard_Lattice_avoid,        \
      Standard_Lattice_primitive, Standard_Lattice_conventional, Standard_Lattice_has_failed, bravais_lattice_type, bravais_lattice_variation_type, bravais_lattice_system, bravais_lattice_lattice_type,            \
      bravais_lattice_lattice_variation_type, bravais_lattice_lattice_system, pearson_symbol, reciprocal_lattice_type, reciprocal_lattice_variation_type, bravais_superlattice_lattice,                              \
      bravais_superlattice_type, bravais_superlattice_variation_type, bravais_superlattice_system, pearson_symbol_superlattice, volume_changed_original2new, transform_coordinates_original2new,                     \
      transform_coordinates_new2original, rotate_lattice_original2new, rotate_lattice_new2original, label_uint, label_int, label_double, order_parameter_structure, order_parameter_atoms, order_parameter_orbit,    \
      order_parameter_sum, partial_occupation_flag, partial_occupation_site_tol, partial_occupation_stoich_tol, partial_occupation_HNF, partial_occupation_sublattice, qm_calculated, qm_scale, qm_lattice,          \
      qm_klattice, qm_origin, qm_f2c, qm_c2f, qm_atoms, qm_forces, qm_forces_write, qm_positions, qm_positions_write, qm_E_cell, qm_dE_cell, qm_H_cell, qm_PV_cell, qm_P, qm_mag_cell, qm_E_atom, qm_dE_atom,        \
      qm_H_atom, qm_PV_atom, qm_mag_atom, kpoints_mode, kpoints_k1, kpoints_k2, kpoints_k3, kpoints_s1, kpoints_s2, kpoints_s3, kpoints_kmax, kpoints_kppra, kpoints_kscheme, dist_nn_min, sym_eps_calculated,       \
      sym_eps, sym_eps_change_count, sym_eps_no_scan, pgroup, pgroup_calculated, pgroup_xtal, pgroup_xtal_calculated, crystal_family, crystal_system, point_group_crystal_class, point_group_Shoenflies,             \
      point_group_Hermann_Mauguin, point_group_orbifold, point_group_type, point_group_order, point_group_structure, pgroupk_Patterson, pgroupk_Patterson_calculated, pgroupk, pgroupk_calculated, pgroupk_xtal,     \
      pgroupk_xtal_calculated, fgroup, fgroup_calculated, sgroup, sgroup_calculated, sgroup_radius, sgroup_radius_dims, agroup_calculated, agroup, iatoms_calculated, iatoms, spacegroup, spacegrouplabel,           \
      spacegroupoption, spacegroupnumber, spacegroupnumberoption, is_spacegroup_platon, is_spacegroup_findsym, is_spacegroup_aflow, crystal_system_ITC, point_group_ITC, bravais_label_ITC, lattice_label_ITC,       \
      space_group_ITC, wyckoff_library_entry_ITC, setting_ITC, origin_ITC, general_position_ITC, wyccar_ITC, standard_lattice_ITC, standard_basis_ITC, wyckoff_sites_ITC, wyckoff_symbols_ITC,                       \
      Niggli_has_failed, Minkowski_has_failed, LatticeReduction_has_failed, write_lattice_flag, write_klattice_flag, write_inequivalent_flag, write_DEBUG_flag, error_flag, error_string
};

void GetNeighbors(const xstructure& xstr_in, std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, double rmin, bool prim, bool unique_only);
void GetNeighbors(const xstructure& xstr_in, std::deque<_atom>& atoms_cell, std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, double rmin, bool prim, bool unique_only);
void GetNeighbors(const xstructure& xstr_in, std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, double rmax, double rmin, bool prim, bool unique_only);
void GetNeighbors(const xstructure& xstr, std::deque<_atom>& atoms_cell, std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, double rmax, double rmin, bool prim, bool unique_only);
//
void GetCoordinations(const xstructure& xstr_in, std::deque<std::deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in, std::deque<_atom>& atoms_cell, std::deque<std::deque<uint>>& coordinations, double rmin, double tol, bool prim,
                      bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in, std::deque<std::deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in, std::deque<_atom>& atoms_cell, std::deque<std::deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim,
                      bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in, std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, std::deque<std::deque<uint>>& coordinations, double rmin, double tol, bool prim, bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in,
                      std::deque<_atom>& atoms_cell,
                      std::deque<std::deque<uint>>& i_neighbors,
                      std::deque<std::deque<double>>& distances,
                      std::deque<std::deque<uint>>& coordinations,
                      double rmin,
                      double tol,
                      bool prim,
                      bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in, std::deque<std::deque<uint>>& i_neighbors, std::deque<std::deque<double>>& distances, std::deque<std::deque<uint>>& coordinations, double rmax, double rmin, double tol, bool prim, bool unique_only); // CO2020914
void GetCoordinations(const xstructure& xstr_in,
                      std::deque<_atom>& atoms_cell,
                      std::deque<std::deque<uint>>& i_neighbors,
                      std::deque<std::deque<double>>& distances,
                      std::deque<std::deque<uint>>& coordinations,
                      double rmax,
                      double rmin,
                      double tol,
                      bool prim,
                      bool unique_only); // CO2020914

void LightCopy(const xstructure&, xstructure&); // ME20200220

aurostd::JSON::object xstructure2json(const xstructure& xstr); // DX20170831 - xstructure2json

// --------------------------------------------------------------------------
// AtomEnvironment Class //DX20191120
class AtomEnvironment {
public:
  AtomEnvironment(); // constructor operator
  ~AtomEnvironment(); // destructor operator
  friend std::ostream& operator<<(std::ostream& oss, const AtomEnvironment& AtomEnvironment);
  // stringstream operator (printing)
  const AtomEnvironment& operator=(const AtomEnvironment& b); // assignment operator
  AtomEnvironment(const AtomEnvironment& b); // copy constructor
  std::string element_center;
  // species/element at center of environment
  uint num_neighbors;
  uint type_center; // type (uint) at center of environment
  uint mode; // AE mode
  uint num_types;
  std::vector<std::string> elements_neighbor; // species/element of atoms neighboring center atom
  std::vector<uint> types_neighbor; // types (uint) of atoms neighboring center atom
  std::vector<double> distances_neighbor;
  // distances to atoms neighboring atoms (typically put in a bin with small tolerance threshold)
  std::vector<uint> coordinations_neighbor;
  // coordination of neighboring distance
  std::vector<std::vector<aurostd::xvector<double>>> coordinates_neighbor;
  // coordinates of atoms neighboring atoms (center is assumed to be zero,i.e. coord=neighbor-origin)
  std::vector<std::vector<uint>> facets;
  // list of facet vertices in order (coordinates_neighbor_flat indexes)  //HE20210408
  std::vector<double> facet_area; // area of each facet //HE20210408
  std::vector<uint> facet_order;
  // count of facet with (index+3) vertices - index 7 counts facets with more than 9 vertices //HE20210408
  bool has_hull;
  double area; // surface area of each environment //HE20210408
  double volume; // volume of each environment //HE20210408
  // functions
  void getAtomEnvironment(const xstructure& xstr, uint center_index, uint ae_mode = ATOM_ENVIRONMENT_MODE_1);
  // get environment around atom index
  void getAtomEnvironment(const xstructure& xstr, uint center_index, const std::vector<std::string>& neighbor_elements, uint ae_mode = ATOM_ENVIRONMENT_MODE_1);
  // get restricted environment (via specified elements) around atom index
  void constructAtomEnvironmentHull(); // construct hull around an environment //HE20210408
  aurostd::xvector<double> index2Point(uint index); // flat view on coordinates_neighbor //HE20210408
  [[nodiscard]] aurostd::JSON::object toJSON(bool full = true) const; // HE20210408 //DX20210624 - added mode input
private:
  void free(); // free operator
  void copy(const AtomEnvironment& b); // copy constructor
};

// --------------------------------------------------------------------------

bool sortAtomsTypes(const _atom& a1, const _atom& a2); // sort atoms by types
bool sortAtomsNames(const _atom& a1, const _atom& a2); // sort atoms by names
bool sortAtomsDist(const _atom& a1, const _atom& a2); // sort atoms by dist  //CO20180420
bool sortAtomsEquiv(const _atom& a1, const _atom& a2); // cluster by equivalent atoms //CO20190116
// sort Wyckoff positions //DX20200515
bool sortWyckoffByLetter(const wyckoffsite_ITC& a, const wyckoffsite_ITC& b);
// sort Wyckoff positions by Wyckoff letter
bool sortWyckoffByType(const wyckoffsite_ITC& a, const wyckoffsite_ITC& b); // sort Wyckoff positions by atom type

// LATTICE/BASIS TRANSFORMATIONS
aurostd::xmatrix<double> GetBasisTransformation(const aurostd::xmatrix<double>& lattice_original,
                                                const aurostd::xmatrix<double>& lattice_new); // DX20201015
std::vector<aurostd::xvector<double>> GetBasisTransformationInternalTranslations(const aurostd::xmatrix<double>& basis_transformation); // DX20201124
aurostd::xmatrix<double> GetRotation(const aurostd::xmatrix<double>& lattice_original,
                                     const aurostd::xmatrix<double>& lattice_new); // DX20201015
xstructure ChangeBasis(const xstructure& xstr, const aurostd::xmatrix<double>& transformation_matrix); // DX20201015
xstructure TransformStructure(const xstructure& xstr, const aurostd::xmatrix<double>& transformation_matrix,
                              const aurostd::xmatrix<double>& rotation); // DX20201125
xstructure TransformStructure(const xstructure& xstr,
                              const aurostd::xmatrix<double>& transformation_matrix,
                              const aurostd::xmatrix<double>& rotation,
                              const aurostd::xvector<double>& origin_shift,
                              bool is_shift_frac = false); // DX20201125
aurostd::matrix<double> GetRotationMatrix(const std::vector<double>& angles);
// CO20200404 pflow::matrix()->aurostd::matrix()  //DX20210127 - moved from pflow
void RotateStrVec(std::vector<xstructure>& str_vec, const std::vector<double>& rot); // DX20210127 - moved from pflow

#endif // AFLOW_XSTRUCTURE_H
