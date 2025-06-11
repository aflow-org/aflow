// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Rico Friedrich, Corey Oses, and Marco Esters
// rico.friedrich@duke.edu

#ifndef _AFLOW_CCE_H_
#define _AFLOW_CCE_H_

#include <deque>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"

#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "structure/aflow_xstructure.h"

// when adding a new functional also introduce new 'offset' in get_offset function needed for reading corrections from lookup table
static const std::vector<std::string> CCE_vallowed_functionals{"PBE", "LDA", "SCAN", "PBE+U:ICSD", "exp"};
// corrections are given for these functionals if only a structure is given as input for the command line and web tools (i.e. --functionals= is not set)
static const std::vector<std::string> CCE_vdefault_output_functionals{"PBE", "LDA", "SCAN", "exp"};
// needs to be extended when adding new corrections for other temperatures
static const std::vector<std::string> CCE_vtemperatures{"298.15", "0"};
static const double _CCE_NN_DIST_TOL_ = 0.5; // 0.5 Ang tolerance between shortest and longest bonds for each cation-anion pair; works best up to now; in future maybe bonding could be explicitly determined via Bader analysis
//[RF20200912 - MOVED TO AFLOW.RC]static const double _CCE_NN_DIST_TOL_MULTI_ANION_ = 0.4; // 0.4 Ang tolerance between shortest and longest bonds for each bond when testing for multi-anion compound; it was
// found that the standard 0.5 Ang tol. is too large such that different anions appear to be bonded, which would prevent anions to be detected as such [RF20200912 - MOVED TO AFLOW.RC]static const double
//_CCE_OX_TOL_ = 0.001; // choose small finite value since sum of oxidation states might not be exactly zero due to numerics [RF20200912 - MOVED TO AFLOW.RC]static const double _CCE_SELF_DIST_TOL_ = 0.001; //
// distance tolerance in Ang for neighbor screening to savely exclude the cation itself having distance zero to itself [RF20200912 - MOVED TO AFLOW.RC]static const double _CCE_perox_cutoff_=1.6; // O-O bonds in
// peroxides for the studied examples are all shorter than 1.6 Ang [RF20200912 - MOVED TO AFLOW.RC]static const double _CCE_superox_cutoff_=1.4; // O-O bonds in superoxides for the studied examples are all
// shorter than 1.4 Ang [RF20200912 - MOVED TO AFLOW.RC]static const double _CCE_O2_molecule_cutoff_=1.2; // O-O bonds in the O2 molecule is about 1.21 Ang.
static const uint CCE_num_functionals_Bader = 4; // Currently, Bader charges used to determine oxidation states are available for 4 functionals: PBE, LDA, SCAN, and PBE+U:ICSD and ONLY for oxides, see get_Bader_templates function

namespace cce {
  struct CCE_Variables {
    std::vector<double> enthalpies_dft;
    std::vector<std::string> vfunctionals; // should be needed as long as output for corrected dft formation enthalpies is based on vfunctionals
    std::vector<int> offset; // needed for reading corrections from lookup table for different functionals
    std::vector<std::string> vtemperatures;
    double standard_anion_charge;
    std::vector<double> electronegativities;
    std::vector<uint> multi_anion_atoms; // vector in which elements will be 1 for multi_anion atoms and 0 otherwise
    std::vector<double> oxidation_states;
    std::string anion_species;
    std::vector<double> cutoffs;
    xstructure xstr_neighbors;
    std::deque<std::deque<uint>> i_neighbors; // CO20200914
    std::deque<std::deque<double>> distances; // CO20200914
    std::vector<std::string> multi_anion_species; // vector storing all the multi anion species
    uint num_perox_bonds;
    uint num_superox_bonds;
    std::vector<uint> perox_indices; // vector in which elements will be 1 for peroxide O atoms and 0 otherwise; needed for correct setting of oxidation numbers below
    std::vector<uint> superox_indices; // vector in which elements will be 1 for superoxide O atoms and 0 otherwise; needed for correct setting of oxidation numbers below
    std::vector<uint> num_neighbors;
    std::vector<std::string> species_electronegativity_sorted;
    std::vector<int> num_pref_ox_states_electronegativity_sorted;
    std::vector<std::vector<double>> pref_ox_states_electronegativity_sorted;
    std::vector<int> num_all_ox_states_electronegativity_sorted;
    std::vector<std::vector<double>> all_ox_states_electronegativity_sorted;
    std::vector<std::vector<int>> cations_map;
    double oxidation_sum; // double because for superoxides O ox. number is -0.5
    std::vector<double> Bader_charges;
    std::vector<std::vector<double>> corrections_atom; // 1st dim. is number of functionals*2 (vfunctionals.size()*2) i.e. 298.15 & 0K corrections for each functional; 2nd dimension for corrections for each atom
    std::vector<std::vector<std::vector<double>>> multi_anion_corrections_atom; // 1st dim. for multi_anion_species, 2nd dim. for functionals and temperatures as above, 3rd dim. for corrections for each atom
    std::vector<double> perox_correction; // peroxide correction per cell for functionals and temperatures as above
    std::vector<double> superox_correction; // superoxide correction per cell for functionals and temperatures as above
    std::vector<double> cce_correction; // total correction per cell for functionals and temperatures as above
    std::vector<double> enthalpy_formation_cell_cce; // CCE formation enthalpy per cell for functionals and temperatures as above
  };

  // main CCE functions
  // for command line use,
  // use inside AFLOW providing directory path or xstructure & functional string or flags and istream for web tool,
  // and CCE core function called by all other main CCE functions
  void run(aurostd::xoption& flags, std::ostream& oss = std::cout); // CO20201105
  void run(aurostd::xoption& flags, std::istream& ist, std::ostream& oss = std::cout); // CO20201105
  void print_corrections(aurostd::xoption& flags, std::ostream& oss = std::cout);
  void print_corrections(xstructure& structure, aurostd::xoption& flags, std::ostream& oss = std::cout); // CO20201105
  void print_corrections(xstructure& structure, aurostd::xoption& flags, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void print_corrections(aurostd::xoption& flags, std::istream& ist, std::ostream& oss = std::cout); // ME20200213 //CO20201105
  void print_cation_coordination_numbers(aurostd::xoption& flags, std::istream& ist, std::ostream& oss = std::cout);
  void print_oxidation_numbers(aurostd::xoption& flags, std::istream& ist, std::ostream& oss = std::cout);
  std::vector<double> calculate_corrections(const std::string& directory_path);
  std::vector<double> calculate_corrections(const xstructure& structure, std::string functional, const std::string& directory_path = aurostd::getPWD(), std::ostream& oss = std::cout);
  std::vector<double> calculate_corrections(const xstructure& structure, std::string functional, std::ofstream& FileMESSAGE, const std::string& directory_path = aurostd::getPWD(), std::ostream& oss = std::cout);
  // void CCE_core(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, const string& directory_path=aurostd::getPWD());
  void CCE_core(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::vector<std::vector<uint>>& multi_anion_num_neighbors, const std::string& directory_path = aurostd::getPWD());
  // read user input (from command line or directory path)
  xstructure read_structure(const std::string& structure_file, int = IOAFLOW_AUTO); // set xstructure mode argument only here and it is automoatically recognized in the main CCE cpp file
  xstructure read_structure(std::istream& ist);
  void get_dft_form_enthalpies_functionals(const std::string& enthalpies_dft_input_str, const std::string& functionals_input_str, CCE_Variables& cce_vars);
  int get_offset(const std::string& functional);
  std::vector<double> get_oxidation_states(const std::string& oxidation_numbers_input_str, const xstructure& structure, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  std::string get_functional_from_aflow_in_outcar(const xstructure& structure, std::string& aflowin_file, std::string& outcar_file);
  // initialise flags and variables
  aurostd::xoption init_flags(); // ME20200213
  CCE_Variables init_variables(const xstructure&); // ME20200213
  // structural analysis
  std::string determine_anion_species(const xstructure& structure, CCE_Variables& cce_vars);
  std::vector<uint> check_for_multi_anion_system(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, double tolerance = DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION);
  std::vector<uint> get_num_neighbors(const xstructure& structure, double tolerance = _CCE_NN_DIST_TOL_);
  std::vector<uint> get_num_neighbors(const xstructure& structure, const std::string& anion_species, double tolerance = _CCE_NN_DIST_TOL_);
  std::vector<uint> get_num_neighbors(const xstructure& structure, const std::string& anion_species, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, double tolerance = _CCE_NN_DIST_TOL_);
  std::vector<double> get_dist_cutoffs(const xstructure& structure);
  void check_per_super_oxides(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars);
  // determine oxidation numbers from electronegativities
  std::vector<double> get_oxidation_states_from_electronegativities(xstructure& structure);
  std::vector<double> get_oxidation_states_from_electronegativities(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void set_anion_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars);
  void sort_species_by_electronegativity(const xstructure& structure, CCE_Variables& cce_vars);
  void load_ox_states_templates_each_species(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars);
  void try_preferred_oxidation_states(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars);
  void treat_SbO2_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void treat_Pb3O4_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void treat_Ti_O_Magneli_phase_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void treat_Fe3O4_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void treat_X3O4_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, const std::string& cation_species, std::ostream& oss = std::cout);
  void treat_alkali_sesquioxide_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  // following special cases only needed when determining oxidation states from Bader charges
  void treat_MnMoO4_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void treat_FeTiO3_LDA_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void check_ox_nums_special_case(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  void try_all_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars);
  void determine_cation_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars, const std::vector<std::vector<double>>& possible_ox_states); // ME20191101
  double get_oxidation_states_sum(CCE_Variables& cce_vars);
  // determine oxidation numbers from Bader charges (outdated)
  std::vector<double> get_oxidation_states_from_Bader(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, const std::string& directory_path = aurostd::getPWD(), std::ostream& oss = std::cout);
  void get_system_name_functional_from_aflow_in(
      const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::string& system_name, std::string& functional, const std::string& directory_path = aurostd::getPWD());
  std::vector<double> get_Bader_charges_from_Bader_file(const xstructure& structure, CCE_Variables& cce_vars, const std::string& Bader_file);
  std::vector<double> Bader_charges_to_oxidation_states(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::string& functional, std::ostream& oss = std::cout);
  void general_attempt_fixing_oxidation_states(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars);
  // assign corrections
  void get_corrections(const xstructure& structure,
                       aurostd::xoption& cce_flags,
                       CCE_Variables& cce_vars,
                       const std::string& considered_anion_species,
                       const std::vector<uint>& num_neighbors,
                       std::vector<std::vector<double>>& corrections_atom,
                       std::ostream& oss = std::cout);
  void load_cation_corrections(const xstructure& structure, CCE_Variables& cce_vars, const std::string& corrections_line, std::vector<std::vector<double>>& corrections_atom, uint i);
  void set_anion_corrections(const xstructure& structure, CCE_Variables& cce_vars, std::vector<std::vector<double>>& corrections_atom, uint i);
  void check_get_per_super_ox_corrections(CCE_Variables& cce_vars);
  // apply corrections and get corrected formation enthalpies
  void check_apply_per_super_ox_corrections(CCE_Variables& cce_vars);
  void apply_pbe_u_icsd_shifts(const xstructure& structure, CCE_Variables& cce_vars, std::ostream& oss = std::cout);
  // print output and citation
  std::string print_output(const xstructure& structure, aurostd::xoption& flags, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::vector<std::vector<uint>>& multi_anion_num_neighbors);
  std::string print_JSON_cation_coordination_numbers(const xstructure& structure, aurostd::xoption& cce_flags, const CCE_Variables& cce_vars, std::vector<std::vector<uint>>& multi_anion_num_neighbors);
  std::string print_JSON_ox_nums(const xstructure& structure, const CCE_Variables& cce_vars);
  std::string print_JSON_corrections(const xstructure& structure, const CCE_Variables& cce_vars); // ME20200213
  std::string print_output_cation_coordination_numbers(const xstructure& structure, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, std::vector<std::vector<uint>>& multi_anion_num_neighbors);
  std::string print_output_oxidation_numbers(const xstructure& structure, CCE_Variables& cce_vars);
  std::string print_output_corrections(const xstructure& structure, CCE_Variables& cce_vars, const std::vector<double>& enthalpy_formation_cell_cce);
  std::string print_test_output(CCE_Variables& cce_vars, const std::vector<double>& enthalpy_formation_cell_cce);
  std::string print_citation();
  // print user instructions
  std::string print_usage();
  // look up tables to load corrections and other data
  std::string get_corrections_line(const std::string& considered_anion_species, const std::string& cor_identifier);
  std::string get_corrections_line_O(const std::string& cor_identifier);
  std::string get_corrections_line_N(const std::string& cor_identifier);
  std::string get_Bader_templates(const std::string& element);
  double get_ref_enthalpy_shift_pbe_u_icsd(const std::string& element);
} // namespace cce

#endif // _AFLOW_CCE_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2021              *
// *                                                                         *
// ***************************************************************************
