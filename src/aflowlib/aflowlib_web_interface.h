
#ifndef AFLOWLIB_WEB_INTERFACE_H
#define AFLOWLIB_WEB_INTERFACE_H

#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"
#include "flow/aflow_xclasses.h"

namespace aflowlib {
  class _aflowlib_entry {
  public:
      // constructor destructor                                 // constructor/destructor
    _aflowlib_entry();                                        // default, just allocate
    ~_aflowlib_entry();                                       // kill everything
    _aflowlib_entry(const _aflowlib_entry& b);                // constructor copy
    _aflowlib_entry(const std::string& file);                      // constructor from file
    const _aflowlib_entry& operator=(const _aflowlib_entry& b); // copy
      // CONTROL
    std::string entry;
    std::vector<std::string> ventry;                       // ventry split by "|"
    std::string auid;
    std::deque<std::string> vauid;                          // AFLOW UNIQUE IDENTIFIER // AFLOW UNIQUE IDENTIFIER SPLIT
    std::string aurl;
    std::deque<std::string> vaurl;                          // AFLOW RESEARCH LOCATOR and TOKENS
    std::string system_name;                                       // ME20190125 - system_name of the calculation
    std::string keywords;
    std::deque<std::string> vkeywords;                  // keywords inside
    std::string aflowlib_date;
    std::vector<std::string> vaflowlib_date;       // CONTAINS 2 DATES: [0]=lock_date, [1]=lib2raw_date
    std::string aflowlib_version;                                  // version
    std::string aflowlib_entries;
    std::vector<std::string> vaflowlib_entries; // this contains the subdirectories that can be associated
    int aflowlib_entries_number;                              // their number
    std::string aflow_version;                                     // version
    std::string catalog;                                           // ICSD,LIB2, etc.
    std::string data_api, data_source, data_language;                // version/source/language
    std::string error_status;                                            // ERROR ??
    std::string author;
    std::vector<std::string> vauthor;
    int calculation_cores;
    double calculation_memory, calculation_time;
    std::string corresponding;
    std::vector<std::string> vcorresponding;
    std::string loop;
    std::vector<std::string> vloop;                         // postprocessing
    int node_CPU_Cores;                                       // computer
    double node_CPU_MHz;                                      // computer
    std::string node_CPU_Model;                                    // computer
    double node_RAM_GB;                                       // computer
      // materials
    std::string Bravais_lattice_orig, Bravais_lattice_relax;        // structure
    std::string code;
    std::string composition;
    std::vector<double> vcomposition;
    bool pocc_parent; //Indicates whether entry is a POCC parent which should be ignored by entry loader. This occurs when vcomposition.empty() and when entry.prototype contains POCC substring but not ARUN substring
    std::string compound;
    double density;
    double density_orig; // DX20190124 - add original crystal info
    std::string dft_type;
    std::vector<std::string> vdft_type;
    double eentropy_cell, eentropy_atom;
    double Egap, Egap_fit;
    std::string Egap_type;
    double energy_cell, energy_atom, energy_atom_relax1;
    double energy_cutoff;
    double delta_electronic_energy_convergence;
    double delta_electronic_energy_threshold;
    uint nkpoints, nkpoints_irreducible, kppra;
    std::string kpoints;
    aurostd::xvector<int> kpoints_nnn_relax, kpoints_nnn_static;
    std::vector<std::string> kpoints_pairs;
    double kpoints_bands_path_grid;
    double enthalpy_cell, enthalpy_atom;
    double enthalpy_formation_cell, enthalpy_formation_atom;
    double enthalpy_formation_cce_300K_cell, enthalpy_formation_cce_300K_atom; // CO20200624
    double enthalpy_formation_cce_0K_cell, enthalpy_formation_cce_0K_atom; // CO20200624
    double entropy_forming_ability; // CO20200624
    double entropic_temperature;
    std::string files;
    std::vector<std::string> vfiles;
    std::string files_LIB;
    std::vector<std::string> vfiles_LIB;
    std::string files_RAW;
    std::vector<std::string> vfiles_RAW;
    std::string files_WEB;
    std::vector<std::string> vfiles_WEB;
    std::string forces;
    std::vector<aurostd::xvector<double>> vforces;
    std::string geometry;
    std::vector<double> vgeometry; // a,b,c and unit_cell_angles (b,c) (a,c) (a,b)
    std::string geometry_orig;
    std::vector<double> vgeometry_orig; // a,b,c and unit_cell_angles (b,c) (a,c) (a,b) //DX20190124 - add original crystal info
    std::string lattice_system_orig, lattice_variation_orig, lattice_system_relax, lattice_variation_relax;
    std::string ldau_TLUJ;
    std::vector<std::vector<double>> vLDAU; // ME20190124
    uint natoms;
    uint natoms_orig; // DX20190124 - add original crystal info
    std::string nbondxx;
    std::vector<double> vnbondxx;
    uint nspecies;
    std::string Pearson_symbol_orig, Pearson_symbol_relax;
    std::string positions_cartesian;
    std::vector<aurostd::xvector<double>> vpositions_cartesian;
    std::string positions_fractional;
    std::vector<aurostd::xvector<double>> vpositions_fractional;
    double pressure; // the true applied pressure (PSTRESS)
    std::string stress_tensor;
    std::vector<double> vstress_tensor; // (1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)
    double pressure_residual; // the leftover pressure due to convergence
    double Pulay_stress; // the leftover pressure due to incomplete basis set
    std::string prototype;
    double PV_cell, PV_atom;
    double scintillation_attenuation_length;
    std::string sg, sg2;
    std::vector<std::string> vsg, vsg2; // CO20180101
    uint spacegroup_orig, spacegroup_relax;
    std::string species;
    std::vector<std::string> vspecies;
    std::string species_pp;
    std::vector<std::string> vspecies_pp;
    std::string species_pp_version;
    std::vector<std::string> vspecies_pp_version;
    std::string species_pp_ZVAL;
    std::vector<double> vspecies_pp_ZVAL;
    std::string species_pp_AUID;
    std::vector<std::string> vspecies_pp_AUID;
    std::string METAGGA; // empty if none, potential type is in xPOTCAR/xOUTCAR
    double spin_cell, spin_atom;
    std::string spinD;
    std::vector<double> vspinD;
    std::string spinD_magmom_orig;
    std::vector<double> vspinD_magmom_orig;
    double spinF;
    std::string sponsor;
    std::vector<std::string> vsponsor;
    std::string stoichiometry;
    std::vector<double> vstoichiometry;
    double valence_cell_std, valence_cell_iupac;
    double volume_cell, volume_atom;
    double volume_cell_orig, volume_atom_orig; // DX20190124 - add original crystal info
    // DX20190124 - added original symmetry info - START
    //  SYMMETRY
    std::string crystal_family_orig;
    std::string crystal_system_orig;
    std::string crystal_class_orig;
    std::string point_group_Hermann_Mauguin_orig;
    std::string point_group_Schoenflies_orig;
    std::string point_group_orbifold_orig;
    std::string point_group_type_orig;
    uint point_group_order_orig;
    std::string point_group_structure_orig;
    std::string Bravais_lattice_lattice_type_orig;
    std::string Bravais_lattice_lattice_variation_type_orig;
    std::string Bravais_lattice_lattice_system_orig;
    std::string Bravais_superlattice_lattice_type_orig;
    std::string Bravais_superlattice_lattice_variation_type_orig;
    std::string Bravais_superlattice_lattice_system_orig;
    std::string Pearson_symbol_superlattice_orig;
    std::string reciprocal_geometry_orig;
    std::vector<double> vreciprocal_geometry_orig;
    double reciprocal_volume_cell_orig;
    std::string reciprocal_lattice_type_orig;
    std::string reciprocal_lattice_variation_type_orig;
    std::string Wyckoff_letters_orig;
    std::string Wyckoff_multiplicities_orig;
    std::string Wyckoff_site_symmetries_orig;
    // DX20190124 - added original symmetry info - END
    // DX20180823 - added more symmetry info - START
    //  SYMMETRY
    std::string crystal_family;
    std::string crystal_system;
    std::string crystal_class;
    std::string point_group_Hermann_Mauguin;
    std::string point_group_Schoenflies;
    std::string point_group_orbifold;
    std::string point_group_type;
    uint point_group_order;
    std::string point_group_structure;
    std::string Bravais_lattice_lattice_type;
    std::string Bravais_lattice_lattice_variation_type;
    std::string Bravais_lattice_lattice_system;
    std::string Bravais_superlattice_lattice_type;
    std::string Bravais_superlattice_lattice_variation_type;
    std::string Bravais_superlattice_lattice_system;
    std::string Pearson_symbol_superlattice;
    std::string reciprocal_geometry_relax;
    std::vector<double> vreciprocal_geometry_relax; // CO20220719 _relax
    double reciprocal_volume_cell;
    std::string reciprocal_lattice_type;
    std::string reciprocal_lattice_variation_type;
    std::string Wyckoff_letters;
    std::string Wyckoff_multiplicities;
    std::string Wyckoff_site_symmetries;
    // DX20180823 - added more symmetry info - END
    // DX20190208 - added anrl info - START
    std::string aflow_prototype_label_orig; // DX20201001 - renamed
    std::string aflow_prototype_params_list_orig; // DX20201001 - renamed
    std::string aflow_prototype_params_values_orig; // DX20201001 - renamed
    std::string aflow_prototype_label_relax; // DX20201001 - renamed
    std::string aflow_prototype_params_list_relax; // DX20201001 - renamed
    std::string aflow_prototype_params_values_relax; // DX20201001 - renamed
    // DX20190208 - added anrl info - END
    std::string pocc_parameters; // CO20200731
    // AGL/AEL
    double agl_thermal_conductivity_300K; //  (W/m*K)
    double agl_debye; //  (K)
    double agl_acoustic_debye; //  (K)
    double agl_gruneisen; //
    double agl_heat_capacity_Cv_300K; //  (kB/cell)
    double agl_heat_capacity_Cp_300K; //  (kB/cell)
    double agl_thermal_expansion_300K; //  (1/K)
    double agl_bulk_modulus_static_300K; //  (GPa)
    double agl_bulk_modulus_isothermal_300K; //  (GPa)
    std::string agl_poisson_ratio_source; //
    double agl_vibrational_free_energy_300K_cell; // (meV/cell) //CT20181212
    double agl_vibrational_free_energy_300K_atom; // (meV/atom) //CT20181212
    double agl_vibrational_entropy_300K_cell; // (meV/cell*K) //CT20181212
    double agl_vibrational_entropy_300K_atom; // (meV/atom*K) //CT20181212
    double ael_poisson_ratio; //
    double ael_bulk_modulus_voigt; //  (GPa)
    double ael_bulk_modulus_reuss; //  (GPa)
    double ael_shear_modulus_voigt; //  (GPa)
    double ael_shear_modulus_reuss; //  (GPa)
    double ael_bulk_modulus_vrh; //  (GPa)
    double ael_shear_modulus_vrh; //  (GPa)
    double ael_elastic_anisotropy; // //CO20181128
    double ael_youngs_modulus_vrh; //  (GPa) //CT20181212
    double ael_speed_sound_transverse; // (m/s) //CT20181212
    double ael_speed_sound_longitudinal; // (m/s) //CT20181212
    double ael_speed_sound_average; // (m/s) //CT20181212
    double ael_pughs_modulus_ratio; // //CT20181212
    double ael_debye_temperature; // (K) //CT20181212
    double ael_applied_pressure; // (GPa) //CT20181212
    double ael_average_external_pressure; // (GPa) //CT20181212
    aurostd::xmatrix<double> ael_stiffness_tensor; // ME20191105
    aurostd::xmatrix<double> ael_compliance_tensor; // ME20191105
    // APL //ME20210927
    double energy_free_vibrational_cell_apl_300K;
    double energy_free_vibrational_atom_apl_300K;
    double entropy_vibrational_cell_apl_300K;
    double entropy_vibrational_atom_apl_300K;
    double energy_internal_vibrational_cell_apl_300K;
    double energy_internal_vibrational_atom_apl_300K;
    double energy_zero_point_cell_apl;
    double energy_zero_point_atom_apl;
    double heat_capacity_Cv_cell_apl_300K;
    double heat_capacity_Cv_atom_apl_300K;
    // QHA  //AS20200831
    double gruneisen_qha; // AS20200831
    double gruneisen_qha_300K; // AS20200903
    double thermal_expansion_qha_300K; // AS20200831
    double modulus_bulk_qha_300K; // AS20200831
    double modulus_bulk_derivative_pressure_qha_300K; // AS20201008
    double heat_capacity_Cv_atom_qha_300K; // AS20201008
    double heat_capacity_Cv_cell_qha_300K; // AS20201207
    double heat_capacity_Cp_atom_qha_300K; // AS20201008
    double heat_capacity_Cp_cell_qha_300K; // AS20201207
    double volume_atom_qha_300K; // AS20201008
    double energy_free_atom_qha_300K; // AS20201008
    double energy_free_cell_qha_300K; // AS20201207
    // BADER
    std::string bader_net_charges;
    std::vector<double> vbader_net_charges; // electrons
    std::string bader_atomic_volumes;
    std::vector<double> vbader_atomic_volumes; // Angst^3
    // legacy
    std::string server;
    std::vector<std::string> vserver;
    std::vector<std::vector<std::string>> vserverdir;
    std::string icsd;
    std::string stoich;
    std::vector<double> vstoich;
    std::string structure_name;
    std::string structure_description;
    double distance_gnd; // distance_gnd
    double distance_tie; // distance_tie
    bool pureA, pureB; // pureA,pureB
    bool fcc, bcc, hcp; // options for lattices
    double stoich_a, stoich_b; // stoich_b,stoich_b
    double bond_aa, bond_ab, bond_bb; // bond_xx // BOND_XX [norm V_ATOM^0.33]
    std::vector<uint> vNsgroup; // vNsgroups
    std::vector<std::string> vsgroup; // vsgroups
    std::vector<xstructure> vstr; // vstructures
    // details from EntryLoader //HE20220913
    std::string el_source_type;
    std::string el_source;
    // dielectric
    aurostd::xmatrix<double> freq_plasma;
    aurostd::xmatrix<double> dielectric_static;
    // functions
    bool FixDescription(); // fix description names
    void clear(); // free space
    uint Load(const std::stringstream& stream, std::ostream& oss); // load from std::stringstream it std is std::cout
    uint Load(const std::string& entry, std::ostream& oss); // load from std::string it std is std::cout
    uint Load(const std::vector<uint64_t>& key_hash, const std::vector<std::string>& content); // Load variant for EntryLoader //HE20220404
    void Set(const std::string& keyword, const std::string& content); // set a class member to content // HE20220404
    void SetByHash(const uint64_t key_hash, const std::string& content); // set a class member to content based on crc64 keyword hash // HE20220404
    uint file2aflowlib(const std::string& file, std::ostream& oss = std::cout); // load from file
    uint url2aflowlib(const std::string& url, std::ostream& oss, bool = true); // load from the web (VERBOSE)
    std::string aflowlib2string(std::string = "out", bool = false); // ME20210408 - added PRINT_NULL
    std::string aflowlib2file(std::string file, std::string = "out"); //
    std::string POCCdirectory2MetadataAUIDjsonfile(const std::string& directory, uint salt = 0); // CO20200624 - get contents of auid_metadata.json
    bool directory2auid(const std::string& directory); // from directory and AURL gives AUID and VAUID
    [[nodiscard]] double enthalpyFormationCell(int T = 300) const; // CO20200624 - CCE correction added

    [[nodiscard]] double enthalpyFormationAtom(bool& cce_used, int T) const; // CO20200624 - CCE correction added
    [[nodiscard]] bool ignoreBadDatabase() const; // CO20171202 - apennsy fixes
    bool ignoreBadDatabase(std::string& reason) const; // CO20171202 - apennsy fixes
    std::string getPathAURL(std::ostream& oss = std::cout, bool load_from_common = false) const; // converts entry.aurl to url/path (common)
    std::string getPathAURL(std::ofstream& FileMESSAGE, std::ostream& oss, bool load_from_common = false) const; // converts entry.aurl to url/path (common)
    [[nodiscard]] std::vector<std::string> getSpecies() const; // CO20221110 - extracts species from restapi
    std::vector<std::string> getSpeciesAURL(std::ostream& oss) const; // CO20210201 - extracts species from aurl
    std::vector<std::string> getSpeciesAURL(std::ofstream& FileMESSAGE, std::ostream& oss) const; // CO20210201 - extracts species from aurl
    // ML stoich features
    void getStoichFeatures(std::vector<std::string>& vheaders, const std::string& e_props = _AFLOW_XELEMENT_PROPERTIES_ALL_);
    void getStoichFeatures(std::vector<std::string>& vheaders, std::vector<double>& vfeatures, bool vheaders_only = false, const std::string& e_props = _AFLOW_XELEMENT_PROPERTIES_ALL_);

  private: //
    void free(); // free space
    void copy(const _aflowlib_entry& b); //
    void LoadCleanup();
  };

  ////////////////////////////////////////////////////////////////////////////////
  // merge vector entries lists
  // 3-vec
  bool mergeEntries(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries_new);
  bool mergeEntries(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, const std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries_new, bool entries_new_same_type = false);
  bool mergeEntries(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, const std::vector<aflowlib::_aflowlib_entry>& entries_new, bool entries_new_same_type = false);
  bool mergeEntries(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, const aflowlib::_aflowlib_entry& entries_new);
  bool mergeEntries(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, const aflowlib::_aflowlib_entry& entries_new, int& index_layer1, int& index_layer2);
  // 2-vec
  bool mergeEntries(std::vector<std::vector<aflowlib::_aflowlib_entry>>& naries, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries_new, bool sort_by_species = true);
  bool mergeEntries(std::vector<std::vector<aflowlib::_aflowlib_entry>>& naries, const std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries_new, bool entries_new_same_type = false, bool sort_by_species = true);
  bool mergeEntries(std::vector<std::vector<aflowlib::_aflowlib_entry>>& naries, const std::vector<aflowlib::_aflowlib_entry>& entries_new, bool entries_new_same_type = false, bool sort_by_species = true);
  bool mergeEntries(std::vector<std::vector<aflowlib::_aflowlib_entry>>& naries, const aflowlib::_aflowlib_entry& entry_new, bool sort_by_species = true);
  bool mergeEntries(std::vector<std::vector<aflowlib::_aflowlib_entry>>& naries, const aflowlib::_aflowlib_entry& entry_new, int& index, bool sort_by_species = true);
  // 1-vec
  bool mergeEntries(std::vector<aflowlib::_aflowlib_entry>& naries, const std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries_new);
  bool mergeEntries(std::vector<aflowlib::_aflowlib_entry>& naries, const std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries_new);
  bool mergeEntries(std::vector<aflowlib::_aflowlib_entry>& naries, const std::vector<aflowlib::_aflowlib_entry>& entries_new);
  bool mergeEntries(std::vector<aflowlib::_aflowlib_entry>& naries, const aflowlib::_aflowlib_entry& entry_new);

  // ***************************************************************************
  bool json2aflowlib(const std::string& json, std::string key, std::string& value); // ME I do not understand why you have out this function as a private member... it would have been simpler to be able to use everywhere
  uint auid2present(std::string auid, std::string& aurl, int mode = 1); // returns json.size() if found...
  std::map<std::string, std::string> AflowlibLocator(const std::vector<std::string>&, const std::string&); // HE20240324 change to use entry loader
  bool AflowlibLocator(const std::string&, std::string&, const std::string&);
  std::string AflowlibLocator(const std::string& options, const std::string& mode);
  std::string AFLUXCall(const aurostd::xoption& vpflow); // DX20190206 - add AFLUX functionality for command line   //CO20200520
  std::string AFLUXCall(const std::vector<std::string>& matchbook); // DX20190206 - add AFLUX functionality //CO20200520
  std::string AFLUXCall(const std::string& summons); // DX20190206 - add AFLUX functionality   //CO20200520
  std::vector<std::vector<std::pair<std::string, std::string>>> getPropertiesFromAFLUXResponse(const std::string& response); // DX20190206 - get properties from AFLUX response  //CO20200520
  std::string getSpaceGroupAFLUXSummons(const std::vector<uint>& space_groups, uint relaxation_step); // DX20200929
  std::string getSpaceGroupAFLUXSummons(uint space_group_number, uint relaxation_step, bool only_one_sg = true); // DX20200929
  std::string AFLUXCall(const aurostd::xoption& vpflow); // DX20190206 - add AFLUX functionality for command line
  std::string AFLUXCall(const std::vector<std::string>& matchbook); // DX20190206 - add AFLUX functionality
  std::string AFLUXCall(const std::string& summons); // DX20190206 - add AFLUX functionality
  std::vector<std::vector<std::pair<std::string, std::string>>> getPropertiesFromAFLUXResponse(const std::string& response); // DX20190206 - get properties from AFLUX response
  std::string getSpaceGroupMatchbook(const std::vector<uint>& space_groups, uint relaxation_step); // DX20200929
  std::string getSpaceGroupMatchbook(uint space_group_number, uint relaxation_step, bool only_one_sg = true); // DX20200929
  void WEB_Aflowlib_Entry(const std::string& option, std::ostream& oss);

  std::string VASPdirectory2auid(const std::string& directory, const std::string& aurl); // CO20200624 - moving from inside _aflowlib_entry
  uint auid2vauid(const std::string auid, std::deque<std::string>& vauid); // splits the auid into vauid
  std::string auid2directory(const std::string auid); // gives AUID directory from existence of vauid

} // namespace aflowlib

#endif // AFLOWLIB_WEB_INTERFACE_H
