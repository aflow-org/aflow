// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalFinder (compare crystal structures) - Functions
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo

#ifndef __AFLOW_COMPARE_STRUCTURE_H_
#define __AFLOW_COMPARE_STRUCTURE_H_

#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"
#include "modules/PROTOTYPES/aflow_anrl.h"
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

// Create the version of GCC, we will uset it for multithread parts of code,
// to check if the current compiling version of gcc is able to compile the
// std::thead features
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

#define JSON_MODE 0
#define _CALCULATE_MAGNETIC_MISFIT_ false
#define _SPIN_TOL_ 0.1

#define _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ 0
#define _COMPARE_DATABASE_GEOMETRY_RELAX1_ 1
#define _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_ 2

#define _DEBUG_COMPARE_ false  // DX20201223

// TODO: for this file, split classes into own files, structs can stay

static std::string GENERAL_XTALFINDER_OPTIONS_LIST =
    "general_options: [--usage] [--misfit_match=0.1] [--misfit_match=0.2] [--np=16|--num_proc=16] [--optimize_match] [--no_scale_volume] [--ignore_symmetry] [--ignore_Wyckoff] [--ignore_environment] "
    "[--keep_unmatched!] [--match_to_aflow_prototypes!] [--magmom=<m1,m2,...|INCAR|OUTCAR>:...] [--add_aflow_prototype_designation] [--remove_duplicate_compounds] [--ICSD] [--print_mapping|--print] "
    "[--print=TEXT|JSON] [--quiet|--q] [--screen_only] [--primitivize] [--minkowski] [--niggli]";

//// ===== GroupedWyckoffPosition Class ===== //
// DX20191120 [MOVED TO aflow.h]

// ===== AtomEnvironment Class ===== //
// DX20191120 [MOVED TO aflow.h]

// ***************************************************************************
// matched_structure_type_xtalfinder (enum) //DX20210112
// ***************************************************************************
// Enum to easily toggle between duplicate and same family structure matches
// Added "xf" to the end of the variable names to help avoid clashing in
// global namespace
enum matched_structure_type_xtalfinder {
  duplicate_structure_xf,                // duplicate structures (misfit < misfit_match)
  family_structure_xf                    // same family structures (misfit_match < misfit < misfit_family)
};

// ***************************************************************************
// output_file_xtalfinder (enum) //DX20210112
// ***************************************************************************
// Determines the file prefix for writing the results
// Added "xf" to the end of the variable names to help avoid clashing in
// global namespace
enum output_file_xtalfinder {
  compare_input_xf,                // DEFAULT_XTALFINDER_FILE_MATERIAL prefix
  duplicate_compounds_xf,          // DEFAULT_XTALFINDER_FILE_DUPLICATE prefix
  compare2database_xf,             // DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE prefix
  compare_database_entries_xf      // DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE prefix
};

// ***************************************************************************
// Struct: structure_mapping_info  //DX20191212 //DX20201218 - updated
// ***************************************************************************
// Stores the structural similarity between two structures
// (i.e., the representative structure and a matched structure)
// The basis transformation, rotation, origin shift, and atom mapping
// information is also stored (new as of 20201218).
struct structure_mapping_info {
  // quantitative similarity/cost functions
  double misfit;                             // structural misfit (=1-(1-lattice_deviation)(1-coordinate_displacement)(1-failure)) (see Burzlaff)
  double lattice_deviation;                  // lattice deviation, captures differences between lattices
  double coordinate_displacement;            // coordinate displacement; captures differences between atom positions (relatively close together)
  double failure;                            // figure of failure; captures differences between atom positions (significantly far apart)
  bool is_magnetic_misfit;                   // boolean indicating if a magnetic system and using magnetic as misfit
  double magnetic_misfit;                    // magnetic misfit (DX inspired by Burzlaff's misfit; =1-(1-magnetic_displacement)(1-magnetic_failure))
  double magnetic_displacement;              // magnetic displacement; captures differences between magnetic moment magnitude (and angle for non-collinear)
  double magnetic_failure;                   // magnetic failure; captures spin flip differences
  // transformation info
  double rescale_factor;                     // rescaling factor for duplicate structure
  aurostd::xmatrix<double> rotation;                  // rotation for duplicate structure (rotate into representative reference frame)
  aurostd::xmatrix<double> basis_transformation;      // basis transformation for duplicate structure (convert to lattice similar to reference)
  aurostd::xvector<double> origin_shift;              // shift origin for duplicate structure (shift to common origin with representative structure)
  // mapping info
  std::vector<uint> atom_map;                     // atom index map: position in vector corresponds to reference index, value corresponds to duplicate structure index
  std::vector<uint> basis_map;                    // atom type map: position in vector corresponds to reference index, value corresponds to duplicate structure type
  std::vector<double> distances_mapped;           // distances between mapped atoms between reference and duplciate structures
  std::vector<aurostd::xvector<double>> vectors_mapped;   // mapping vectors between mapped atoms between reference and duplicate structures
};

// ***************************************************************************
// structure_mapping_info functions
// ***************************************************************************
namespace compare {
  structure_mapping_info initialize_misfit_struct(bool magnetic = false);
  void resetMisfitInfo(structure_mapping_info& str_mis, bool magnetic = false); // DX20220406
  std::string printAtomMappings(const structure_mapping_info& misfit_info);
  std::string printUnmatchedAtoms(const structure_mapping_info& misfit_info, const xstructure& xstr1, const xstructure& xstr2);
  void resetMappingInfo(structure_mapping_info& misfit_info); // DX20220406
  void resizeMappingInfo(structure_mapping_info& str_mis, uint size); // DX20220406
} // namespace compare

// DX20200225 - temp struct; working on more robust scheme
struct matching_structure {
  std::string name;
  double misfit;
};

// ***************************************************************************
// Struct: structure_container //DX20201218
// ***************************************************************************
struct structure_container {
  // intialization info
  std::string name;                                                  // name of structure
  std::string link;
  std::string compound;                                              // compound name of structure
  std::vector<uint> stoichiometry;                                   // reduced stoichiometry of structure
  uint natoms;                                                  // number of atoms in unit cell
  uint ntypes;                                                  // number of atom types
  std::vector<std::string> elements;                                      // element names in structure
  std::string source;                                                // indicating where structure came from, i.e., input, auid, aflow protos (for structure regeneration if necessary)
  bool is_structure_generated;                                  // boolean indicating if the structure has been generated
  uint relaxation_step;                                         // specifies the relaxation step of the structure (0=original, 1=one relaxation, 2=most relaxed)
  // structure info
  xstructure structure;                                         // xstructure
  std::vector<double> nearest_neighbor_distances;                    // nearest neighbor distances for atoms, vector position corresponds to atom index
  std::vector<AtomEnvironment> environments_LFA;                     // LFA atom environments (for near isoconfigurational comparison)
  // symmetry info
  std::string Pearson;                                               // Pearson symbol
  uint space_group;                                             // space group
  std::vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;     // Wyckoff positions grouped by species/types
  uint number_compounds_matching_structure;                     // number of compounds that match to this structure
  // property info
  std::vector<std::string> properties_names;                              // API keyword for calculated properties (database comparisons only)
  std::vector<std::string> properties_units;                              // units for calculated properties (database comparisons only)
  std::vector<std::string> properties_types;                              // datatype for calculated properties (database comparison only)
  std::vector<std::string> properties;                                    // value(s) of calculated properties (database comparison only)
};

// ---------------------------------------------------------------------------
namespace compare {
  structure_container initializeStructureContainer();
  structure_container initializeStructureContainer(const xstructure& structure, bool same_species);
} // namespace compare

// DX 20201119
//  ===== CompapareStructureContainers Class ===== //
//  This class is necessary to pass an argument to the sorting function;
//  the other solution is to use std::bind, but this is not available for early G++ versions
//  see https://stackoverflow.com/questions/26444216/is-it-possible-to-use-stdsort-with-a-sort-function-that-takes-extra-arguments
class CompareStructureContainers {
public:
  CompareStructureContainers(const std::vector<std::string>& sorting_attributes) : sorting_attributes(sorting_attributes) {}
  bool operator()(const structure_container* a, const structure_container* b);

private:
  std::vector<std::string> sorting_attributes;
};

// ***************************************************************************
// Class: StructurePrototype
// ***************************************************************************
// This class provides the unique prototype information, including the
// symmetry, environments, and structures that exhibit this prototype
// ---------------------------------------------------------------------------
class StructurePrototype {
public:
  StructurePrototype() = default;                                                                   // constructor operator
  void clear();                                                                           // clear
  friend std::ostream& operator<<(std::ostream& oss, const StructurePrototype& StructurePrototype); // stringstream operator (printing)
  int iomode = JSON_MODE;                                                                             // mode for printing

    // ---------------------------------------------------------------------------
    // structure prototype information
  uint natoms = 0;                                                                            // number of atoms in the prototype (from the representative structure; not necessarily reduced)
  int ntypes = 0;                                                                             // number of types in prototype
  std::vector<std::string> elements;                                                                // list of elements exhibiting in this protoype (from representative and duplicate structures)
  std::vector<uint> stoichiometry;                                                             // reduced stoichiometry of prototype
  std::vector<std::vector<std::string>> atom_decorations_equivalent;                                    // equivalent atom decorations (permutations) of prototype
  std::string Pearson;                                                                         // Pearson symbol of prototype
  uint space_group = 0;                                                                       // space group number of prototype
  std::vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;                               // Wyckoff positions grouped by site type
  std::vector<std::string> wyckoff_site_symmetry;                                                   // vector of Wyckoff site symmetries of prototype
  std::vector<int> wyckoff_multiplicity;                                                       // vector of Wyckoff multiplicities of prototype
  std::vector<int> wyckoff_letter;                                                             // vector of Wyckoff letters of prototype
  std::vector<AtomEnvironment> environments_LFA;                                               // vector of LFA atom environments
  std::string aflow_label;                                                                     // AFLOW label designation
  std::vector<std::string> aflow_parameter_list;                                                    // vector of strings corresponding to AFLOW parameter variables
  std::vector<double> aflow_parameter_values;                                                  // vector of doubles corresponding to AFLOW parameter values
  std::vector<std::string> matching_aflow_prototypes;                                               // vector of strings indicating matching AFLOW prototype labels

    // ---------------------------------------------------------------------------
    // structures
  structure_container* structure_representative;                                          // structural information for representative structure
  std::vector<structure_container*> structures_duplicate;                                      // structural information for the duplicate structures
  std::vector<structure_container*> structures_family;                                         // structural information for the same family structures
  std::vector<structure_mapping_info> mapping_info_duplicate;                                  // structural misfit and transformation information between the duplicate and representative structures
  std::vector<structure_mapping_info> mapping_info_family;                                     // structural misfit and transformation information between the family and representative structures

    // ---------------------------------------------------------------------------
    // properties stored in structure containers
  std::vector<std::string> property_names;                                                          // vector of property names (if using AFLUX)
  std::vector<std::string> property_units;                                                          // vector of property units (if using AFLUX)
  std::vector<std::string> property_types;                                                          // vector of property types (if using AFLUX) //DX20201230

    // ---------------------------------------------------------------------------
    // functions
  [[nodiscard]] uint numberOfDuplicates() const; // DX20190506                                           // return the number of duplicate structures for this prototype (i.e., checks misfit value)
  [[nodiscard]] std::string printRepresentativeStructure() const; // DX20201028                               // return the representative structure in a JSON format
  [[nodiscard]] std::string printMatchedStructures(matched_structure_type_xtalfinder mode) const;            // return the matched structures in a JSON format //DX20201028
  std::string printPropertiesOfStructure(structure_container* str_pointer) const;              // return properties of structure in a JSON format
  [[nodiscard]] std::string printStructureTransformationInformation(const structure_mapping_info& misfit_info) const; // return structure transformation information between structure and representative in JSON format
  [[nodiscard]] uint numberOfComparisons() const; // DX20181221                                                // return the number of comparisons for this prototype
  [[nodiscard]] bool isSymmetryCalculated() const; // DX20190228
  [[nodiscard]] bool isLFAEnvironmentCalculated() const; // DX20191105
  bool calculateSymmetry(); // DX20190118                                                  // calculate space group symmetry and populate symmetry attributes for StructurePrototype
  void putDuplicateAsFamily(uint index, bool keep_generated = false); // make duplicate structure a same family structure in the same StructurePrototypeObject //DX20190814
  void copyPrototypeInformation(const StructurePrototype& b); // copy prototype information from one StructurePrototype object to another
  void copyDuplicate(const StructurePrototype& b, uint index, bool copy_misfit = false); // copy duplicate info to another StructurePrototype object
  void removeNonDuplicate(uint index); // remove non-duplicate structures
};

// ***************************************************************************
// Class: XtalFinderCalculator //DX20201201
// ***************************************************************************
// This is the main calculator class for XtalFinder
// Carries "universal" attributes to the comparison functions:
//   misfit_match         : misfit threshold for matching structures
//   misfit_family        : misfit threshold for structures of the same family
//   num_proc             : number of processors (parallel)
//   structure_containers : relevant crystal structures to analyze to
//                          efficiently pass structures so we do not
//                          constantly copy or move them around.
// The StructurePrototype class points to the structure_containers, and
// aggregates the pointers based on structural simlarity.
// The XtalFinderCalculator class "owns" the address of the structure,
// to prevent memory leaks (e.g., if StructurePrototype "owned" it, we would
// have a memory leak once the instance went out of scope).
// Inheritance with xStream
// ---------------------------------------------------------------------------
class XtalFinderCalculator : public xStream {
public:
  // ---------------------------------------------------------------------------
  // constructors
  XtalFinderCalculator(uint num_proc_input = 1, std::ostream& oss = std::cout) : xStream(oss), num_proc(num_proc_input) {}
  XtalFinderCalculator(std::ofstream& FileMESSAGE, uint num_proc_input = 1, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss), num_proc(num_proc_input) {}
  XtalFinderCalculator(double misfit_match_input, double misfit_family_input, uint num_proc_input = 1, std::ostream& oss = std::cout) :
      xStream(oss), num_proc(num_proc_input), misfit_match(misfit_match_input), misfit_family(misfit_family_input) {}
  XtalFinderCalculator(double misfit_match_input, double misfit_family_input, std::ofstream& FileMESSAGE, uint num_proc_input = 1, std::ostream& oss = std::cout) :
      xStream(FileMESSAGE, oss), num_proc(num_proc_input), misfit_match(misfit_match_input), misfit_family(misfit_family_input) {}
  // ---------------------------------------------------------------------------
  // clear
  void clear();
  // ---------------------------------------------------------------------------
  // operator<<
  friend std::ostream& operator<<(std::ostream& oss, const XtalFinderCalculator& XtalFinderCalculator);
  // ---------------------------------------------------------------------------
  // attributes
  double misfit_match = DEFAULT_XTALFINDER_MISFIT_MATCH;
  double misfit_family = DEFAULT_XTALFINDER_MISFIT_FAMILY;
  uint num_proc = 1;
  std::vector<structure_container> structure_containers; // stores structures in a container (pointer for easy manipulation and mobility)
  anrl::ProtoData pd = anrl::ProtoData::get();

  // ---------------------------------------------------------------------------
  // compare methods
  void compareStructures(structure_container& str_rep, structure_container& str_matched, structure_mapping_info& match_info, bool same_species, bool scale_volume, bool optimize_match);
  std::vector<StructurePrototype> compareStructuresFromStructureList(const std::vector<std::string>& filenames, std::vector<std::string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); // DX20200103 - condensed bools to xoptions
  std::vector<StructurePrototype> compareStructuresFromDirectory(const std::string& directory, std::vector<std::string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); // DX20200103 - condensed bools to xoptions
  std::vector<StructurePrototype> compareStructuresFromString(const std::string& structures_string, std::vector<std::string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); // ME202010206
  std::vector<StructurePrototype> compareStructuresFromFile(const std::string& filename, std::vector<std::string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); // DX20200103 - condensed bools to xoptions

  // ---------------------------------------------------------------------------
  // compare multiple structures
  std::vector<StructurePrototype> compareMultipleStructures(uint num_proc, bool same_species, const std::string& directory);
  std::vector<StructurePrototype> compareMultipleStructures(uint num_proc, bool same_species, const std::string& directory, const aurostd::xoption& comparison_options);

  // ---------------------------------------------------------------------------
  // compare2prototypes
  std::string printMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow);
  std::vector<StructurePrototype> compare2prototypes(const xstructure& xstrIN, const aurostd::xoption& vpflow);
  void calculateMatchingAFLOWPrototypes(std::vector<StructurePrototype>& prototypes, uint num_proc);
  void getMatchingAFLOWPrototypes(uint i, std::vector<StructurePrototype>& prototypes, const aurostd::xoption& vpflow_protos);

  // ---------------------------------------------------------------------------
  // compare2database
  std::vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow);

  // ---------------------------------------------------------------------------
  // compare permuations
  std::vector<std::string> getUniquePermutations(xstructure& xstr, uint num_proc = 1, bool optimize_match = false);
  std::vector<std::string> getUniquePermutations(xstructure& xstr, aurostd::xoption& comparison_options, std::stringstream& results_ss, uint num_proc = 1, bool print_misfit = false);
  void compareAtomDecorations(StructurePrototype& structure, uint num_proc, bool optimize_match);
  void compareAtomDecorations(StructurePrototype& structure, uint num_proc, aurostd::xoption& permutation_options);
  void compareAtomDecorations(StructurePrototype& structure, std::string& misfit_results, uint num_proc, aurostd::xoption& permutation_options);
  void generateAtomPermutedStructures(structure_container& structure);
  std::vector<std::string> getSpeciesPermutedStrings(const std::deque<uint>& stoichiometry);
  std::vector<std::string> getSpeciesPermutedStrings(const std::vector<uint>& stoichiometry);

  // ---------------------------------------------------------------------------
  // get command line options
  void getOptions(const aurostd::xoption& vpflow, aurostd::xoption& comparison_options);

  // ---------------------------------------------------------------------------
  // get command line options
  std::string getSpaceGroupMatchbookFromOptions(const aurostd::xoption& vpflow, uint relaxation_step); // DX20210615 - uint not bool

  // ---------------------------------------------------------------------------
  // add structures to container
  void addStructure2container(const xstructure& xstr, const std::string& structure_name, const std::string& source, uint relaxation_step, bool same_species);

  void addAFLOWPrototypes2container(const std::set<std::string>& vuid);

  void addDatabaseEntry2container(aflowlib::_aflowlib_entry& entry, const std::vector<std::string>& species, uint relaxation_step, bool same_species);

  // ---------------------------------------------------------------------------
  // remove methods
  void removeStructureFromContainerByName(const std::string& structure_name);

  // ---------------------------------------------------------------------------
  // set as structure representative
  void setStructureAsRepresentative(StructurePrototype& structure_tmp, uint container_index);
  void setStructureAsRepresentative(StructurePrototype& structure_tmp, structure_container* str_pointer);
  // ---------------------------------------------------------------------------
  // set as duplicate structure
  void addStructure2duplicatesList(StructurePrototype& structure_tmp, uint container_index);
  void addStructure2duplicatesList(StructurePrototype& structure_tmp, structure_container* str_pointer);
  // ---------------------------------------------------------------------------
  // set as same family structure
  void addStructure2sameFamilyList(StructurePrototype& structure_tmp, uint container_index);
  void addStructure2sameFamilyList(StructurePrototype& structure_tmp, structure_container* str_pointer);

  // ---------------------------------------------------------------------------
  // load structure methods
  void loadStructuresFromStructureList(const std::vector<std::string>& filenames,
                                       const std::vector<std::string>& magmoms_for_systems,
                                       bool same_species); // DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
  void loadStructuresFromDirectory(const std::string& directory, const std::vector<std::string>& magmoms_for_systems, bool same_species); // DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
  void loadStructuresFromFile(const std::string& filename, const std::vector<std::string>& magmoms_for_systems, bool same_species); // DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
  void loadStructuresFromStringstream(std::stringstream& input_stream, const std::vector<std::string>& magmoms_for_systems, bool same_species); // ME20210206
  void loadStructuresFromAflowlibEntries(const std::vector<aflowlib::_aflowlib_entry>& entries, const std::vector<std::string>& magmoms_for_systems, bool same_species); // DX20201201

  // ---------------------------------------------------------------------------
  // transform structures
  void performStructureConversions(uint i, // ME20220207 - new xThread scheme
                                   const std::vector<bool>& calculate_primitive_vec,
                                   const std::vector<bool>& calculate_Minkowski_vec,
                                   const std::vector<bool>& calculate_Niggli_vec); // DX20210113
  void convertStructures(const aurostd::xoption& comparison_options, uint num_proc); // DX20201005
  void getPrimitiveStructures(uint start_index = 0, uint end_index = AUROSTD_MAX_UINT); // DX20201005
  void getMinkowskiStructures(uint start_index = 0, uint end_index = AUROSTD_MAX_UINT); // DX20201005
  void getNiggliStructures(uint start_index = 0, uint end_index = AUROSTD_MAX_UINT); // DX20201005

  // ---------------------------------------------------------------------------
  // analyze symmetry
  bool isSymmetryCalculated(structure_container& structure);
  void calculateSymmetry(structure_container& str_rep);
  void calculateSymmetries(uint num_proc);
  void calculateSpaceGroups(uint start_index = 0, uint end_index = AUROSTD_MAX_UINT, uint setting = 0); // DX20191230 added setting
  void setSymmetryPlaceholders();
  void getSymmetryInfoFromXstructure(structure_container& str_rep); // DX20210104
  // ---------------------------------------------------------------------------
  // analyze environment
  bool isLFAEnvironmentCalculated(structure_container& structure);
  void computeLFAEnvironment(structure_container& str_rep, bool unique_only = true);
  void calculateLFAEnvironments(uint num_proc);
  void computeLFAEnvironments(uint start_index = 0, uint end_index = AUROSTD_MAX_UINT);
  // ---------------------------------------------------------------------------
  // analyze neighbors
  bool areNearestNeighborsCalculated(structure_container& structure);
  void getNearestNeighbors(uint num_proc);
  void calculateNearestNeighbors(uint start_index = 0, uint end_index = AUROSTD_MAX_UINT);

  // ---------------------------------------------------------------------------
  // group structures
  std::vector<StructurePrototype> groupStructurePrototypes(bool same_species,
                                                           bool ignore_symmetry,
                                                           bool ignore_Wyckoff,
                                                           bool ignore_environment,
                                                           bool ignore_environment_angles,
                                                           bool duplicates_removed,
                                                           bool quiet = false); // DX20190731 - remove const and & //DX20190830 - added duplicates_removed //DX20200320 - added environment angles

  // ---------------------------------------------------------------------------
  // reorder structures
  void representativePrototypeForICSDRunsNEW(std::vector<StructurePrototype>& comparison_schemes);
  void makeRepresentativeEvenPermutation(std::vector<StructurePrototype>& comparison_schemes, const std::vector<std::string>& name_order);

  // ---------------------------------------------------------------------------
  // find duplicate compounds //DX20210114
  void findDuplicateCompounds(uint num_proc, bool remove_duplicates, const std::string& directory, const aurostd::xoption& comparison_options);

  // ---------------------------------------------------------------------------
  // check for better structure matches
  std::vector<StructurePrototype> checkForBetterMatches(std::vector<StructurePrototype>& prototype_schemes, uint num_proc, bool check_for_better_matches, bool same_species, const aurostd::xoption& comparison_options, bool quiet);
  void combinePrototypesOfDifferentSymmetry(std::vector<StructurePrototype>& prototypes_final, bool same_species, uint num_proc);

  // ---------------------------------------------------------------------------
  // append unmatched structures into new groups
  void appendStructurePrototypes(std::vector<StructurePrototype>& comparison_schemes,
                                 std::vector<StructurePrototype>& prototypes_final,
                                 bool clean_unmatched, // DX20190506
                                 bool quiet);

  // ---------------------------------------------------------------------------
  // count matched/unmatched
  [[nodiscard]] uint numberOfMismatches(const std::vector<StructurePrototype>& comparison_schemes) const;
  [[nodiscard]] uint numberOfDuplicates(const StructurePrototype& prototype) const;

  // ---------------------------------------------------------------------------
  // split comparisons into threads (2D-array splitting)
  bool splitComparisonIntoThreads(std::vector<StructurePrototype>& comparison_schemes, uint num_proc, std::vector<std::pair<uint, uint>>& start_indices, std::vector<std::pair<uint, uint>>& end_indices);

  // ---------------------------------------------------------------------------
  // run multiple structures
  std::vector<StructurePrototype> runComparisonScheme(std::vector<StructurePrototype>& comparison_schemes,
                                                      bool same_species,
                                                      uint num_proc,
                                                      const aurostd::xoption& comparison_options,
                                                      bool quiet = false); // DX20200103 - condensed bools to xoptions
  void runComparisons(std::vector<StructurePrototype>& comparison_schemes, bool same_species, bool scale_volume, bool optimize_match);
  void runComparisonThreads(uint index,
                            std::vector<StructurePrototype>& comparison_schemes,
                            const std::vector<std::pair<uint, uint>>& vstart_indices,
                            const std::vector<std::pair<uint, uint>>& vend_indices,
                            bool same_species,
                            bool scale_volume,
                            bool optimize_match); // DX20190822 - added comparison log bool

  // ---------------------------------------------------------------------------
  // get aflow label
  void calculatePrototypeDesignations(std::vector<StructurePrototype>& prototypes, uint num_proc);
  // ME20220207 - replaced with iterator for xThread
  void getPrototypeDesignations(std::vector<StructurePrototype>::iterator& prototypes);
  // void getPrototypeDesignations(vector<StructurePrototype>& prototypes);
  void getPrototypeDesignationsPreDistributed(std::vector<StructurePrototype>& prototypes, uint start_index = 0,
                                              uint end_index = AUROSTD_MAX_UINT); // DX20191122

  // ---------------------------------------------------------------------------
  // print
  [[nodiscard]] std::string printResults(const std::vector<StructurePrototype>& prototypes_final, bool same_species, filetype format) const;
  std::string printStructureMappingResults(const structure_mapping_info& misfit_info, const xstructure& xstr_reference, const xstructure& xstr_mapped, const std::string& mode = "TEXT");
  std::string printAtomMappings(const structure_mapping_info& misfit_info, const xstructure& xstr1, const xstructure& xstr2);
  std::string printUnmatchedAtoms(const structure_mapping_info& misfit_info, const xstructure& xstr1, const xstructure& xstr2);

  // ---------------------------------------------------------------------------
  // lattice search
  void latticeSearch(structure_container& xstr_rep,
                     structure_container& xstr_match,
                     structure_mapping_info& match_info,
                     bool same_species,
                     bool optimize_match,
                     bool scale_volume, // DX20200422
                     uint num_proc); // DX20201123

  // ---------------------------------------------------------------------------
  // find translation vectors
  void findSimilarTranslationVectors(const aurostd::xmatrix<double>& q1, const xstructure& xstr_LFA_supercell, const xstructure& xstr, std::vector<aurostd::xvector<double>>& lattice_vecs);

  // ---------------------------------------------------------------------------
  // similar lattices
  bool buildSimilarLattices(const std::vector<aurostd::xvector<double>>& translation_vectors,
                            const aurostd::xmatrix<double>& q1,
                            std::vector<aurostd::xmatrix<double>>& lattices,
                            std::vector<double>& latt_devs,
                            bool optimize_match,
                            bool scale_volume) const;

  // ---------------------------------------------------------------------------
  // search for atom mappings
  bool searchAtomMappings(uint start_index,
                          uint end_index,
                          const xstructure& xstr1,
                          const std::vector<double>& all_nn1,
                          const xstructure& xstr2,
                          const std::string& lfa,
                          std::vector<aurostd::xmatrix<double>>& lattices,
                          std::vector<structure_mapping_info>& vstrs_matched,
                          bool same_species,
                          bool optimize_match);

  // ---------------------------------------------------------------------------
  // find matches (atoms)
  bool findMatch(const xstructure& xstr1,
                 const xstructure& xstr2,
                 const std::vector<uint>& atom_indices_xstr1,
                 const std::vector<uint>& atom_indices_xstr2,
                 double minimum_interatomic_distance, // DX20200622
                 structure_mapping_info& mapping_info,
                 bool same_species);

  // ---------------------------------------------------------------------------
  // helper functions for external use
  std::vector<std::vector<uint>> groupSimilarXstructures(const std::vector<xstructure>& vxstrs, aurostd::xoption& comparison_options, bool same_species = true);
  std::vector<std::vector<uint>> groupSimilarXstructures(const std::vector<xstructure>& vxstrs, bool same_species = true);

  // ---------------------------------------------------------------------------
  // write output to file
  void writeComparisonOutputFile(const std::string& output, const std::string& directory, filetype format, output_file_xtalfinder comparison_mode, bool same_species);
};

namespace compare {
  // ---------------------------------------------------------------------------
  // pair-wise comparisons (for use by other AFLOW processes)
  std::string compareInputStructures(aurostd::xoption& vpflow, std::istream& cin, std::ostream& logstream = std::cout); // ME20210206
  std::string compareInputStructures(const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX //DX20190425
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optimize_match, double& final_misfit, uint num_proc = 1); // Main function //DX20191108 - remove const & from bools
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optimize_match, double& final_misfit, structure_mapping_info& final_misfit_info, uint num_proc = 1); // Main function //DX20191108 - remove const & from bools
  // bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, uint num_proc=1); //Overco, returns true (match), false (no match) //DX20191108 - remove const & from bools
  // bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optmize_match, uint num_proc=1);  //DX20191108 - remove const & from bools
  bool structuresMatch(const xstructure& xstr1, const xstructure& xstr2, bool same_species, uint num_proc = 1); // Overco, returns true (match), false (no match) //DX20191108 - remove const & from bools
  bool structuresMatch(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optmize_match, uint num_proc = 1); // DX20191108 - remove const & from bools
  // double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match, uint num_proc=1); //Overloaded, returns misfit value //DX20191108 - remove const & from bools
  double getMisfitBetweenStructures(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match, uint num_proc = 1); // Overloaded, returns misfit value //DX20191108 - remove const & from bools
  structure_mapping_info getTransformationBetweenStructures(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match, uint num_proc = 1); // Overloaded, returns misfit value //DX20191108 - remove const & from bools

  // ---------------------------------------------------------------------------
  // comparisons to AFLOW database
  std::vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX20200225
  std::vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // DX20200225
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX20200225
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // DX20200225
  std::vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX20200225
  std::vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // DX20200225
  std::string printCompare2Database(std::istream& input, const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX20200225
  std::string printCompare2Database(std::istream& input, const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // CO20200225
  std::string printCompare2Database(const xstructure& xstrIN, const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // DX20200225

  // ---------------------------------------------------------------------------
  // comparisons between entries in AFLOW database
  std::string compareDatabaseEntries(const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX20191125
  std::string compareDatabaseEntries(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // DX20191125

  // ---------------------------------------------------------------------------
  // comparisons aflowlib entries
  std::vector<aflowlib::_aflowlib_entry> getUniqueEntries(std::vector<aflowlib::_aflowlib_entry>& entries, uint num_proc, bool same_species, bool scale_volume, bool optimize_match); // DX20201111
  std::vector<aflowlib::_aflowlib_entry> getUniqueEntries(std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss, std::ofstream& FileMESSAGE, uint num_proc, bool same_species, bool scale_volume, bool optimize_match); // DX20201111

  // ---------------------------------------------------------------------------
  // comparisons to AFLOW prototype library
  std::vector<StructurePrototype> compare2prototypes(std::istream& input, const aurostd::xoption& vpflow, std::ostream& logstream = std::cout); // DX20181004 //DX20190314 - changed return value
  std::vector<StructurePrototype> compare2prototypes(std::istream& input, const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout); // DX20181004 //DX20190314 - changed return value
  std::string printMatchingPrototypes(std::istream& cin, const aurostd::xoption& vpflow); // DX20190314
  std::vector<std::string> getMatchingPrototypes(xstructure& xstr, const std::string& catalog); // DX20190314
  std::vector<std::string> getMatchingPrototypes(xstructure& xstr, const std::string& catalog); // DX20190314
  std::vector<std::string> getMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow, const std::string& catalog); // DX20210421

  // ---------------------------------------------------------------------------
  // permutaion comparisons
  std::string compareAtomDecorations(std::istream& input, const aurostd::xoption& vpflow); // DX20181004

  // ---------------------------------------------------------------------------
  // isopointal AFLOW prototype functions
  std::string isopointalPrototypes(std::istream& input, const aurostd::xoption& vpflow); // DX20200131
  std::set<std::string> getIsopointalPrototypes(xstructure& xstr, const std::string& catalog); // DX20200131

  // ---------------------------------------------------------------------------
  // comparison options
  aurostd::xoption loadDefaultComparisonOptions(const std::string& mode = ""); // DX20200103

  // ---------------------------------------------------------------------------
  // geneate structures
  void generateStructures(std::vector<StructurePrototype>& structures, std::ostream& oss = std::cout, uint start_index = 0, uint end_index = AUROSTD_MAX_UINT); // DX20191122
  bool generateStructure(const std::string& structure_name, const std::string& structure_source, uint relaxation_step, xstructure& structure, std::ostream& oss); // DX20200429 - added relaxation_step
  void removeNonGeneratedStructures(std::vector<StructurePrototype>& structures); // DX20191105

  // ---------------------------------------------------------------------------
  // functions for determining isopointal structures (same/compatible symmetry)
  void groupWyckoffPositions(const xstructure& xstr, std::vector<GroupedWyckoffPosition>& grouped_positions);
  void groupWyckoffPositions(const std::vector<wyckoffsite_ITC>& wyckoff_sites_ITC, std::vector<GroupedWyckoffPosition>& grouped_positions); // DX20200512
  void groupWyckoffPositionsFromGroupedString(uint space_group_number, uint setting, const std::vector<std::vector<std::string>>& grouped_Wyckoff_string, std::vector<GroupedWyckoffPosition>& grouped_positions); // DX20200622 - removed pointer to uints
  std::string printWyckoffString(const std::vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize = false);
  std::vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(const std::vector<GroupedWyckoffPosition>& grouped_Wyckoffs); // DX20190219
  bool matchableWyckoffPositions(const std::vector<GroupedWyckoffPosition>& grouped_Wyckoffs_str1, const std::vector<GroupedWyckoffPosition>& grouped_Wyckoffs_str2, bool same_species);
  bool matchableWyckoffPositionSet(const std::vector<std::vector<std::vector<std::string>>>& grouped_possible_Wyckoff_letters, const std::vector<std::vector<std::string>>& grouped_Wyckoff_letters);
  std::vector<std::vector<std::string>> convertANRLWyckoffString2GroupedPositions(const std::string& label);
  std::vector<std::vector<std::string>> convertWyckoffString2GroupedPositions(const std::string& Wyckoff_letter_string);
  bool sameStoichiometry(const std::vector<uint>& stoich1, const std::vector<uint>& stoich2);
  bool matchableSpaceGroups(uint space_group_1, uint space_group_2);
  bool matchableEnantiomorphicSpaceGroups(uint space_group_1, uint space_group_2);
  bool filterPrototypes(uint& species_count,
                        std::string& reduced_stoichiometry,
                        uint& space_group,
                        std::vector<std::vector<std::vector<std::string>>>& grouped_possible_Wyckoff_letters,
                        std::vector<std::string>& prototype_labels,
                        std::vector<uint>& species_counts,
                        std::vector<uint>& space_groups);
  bool structuresCompatible(const structure_container& structure1, const structure_container& structure2, bool same_species, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool ignore_environment_angles); // DX20201207

  // ---------------------------------------------------------------------------
  // comparing permutations/atom decorations
  std::vector<std::pair<uint, uint>> calculateDivisors(uint number);
  bool checkNumberOfGroupings(const std::vector<StructurePrototype>& comparison_schemes, uint number);
  void generatePermutationString(const std::deque<uint>& stoichiometry, std::vector<std::string>& permutation); // DX20190508 //DX20191125 - changed from std::vector to std::deque
  void generatePermutationString(const std::vector<uint>& stoichiometry, std::vector<std::string>& permutation); // DX20190508
  bool generatePermutations(uint& num_elements,
                            std::vector<uint>& indices,
                            std::vector<std::string>& names,
                            std::vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions,
                            std::vector<std::vector<uint>>& permutations,
                            std::vector<std::vector<std::string>>& name_order,
                            std::vector<std::vector<GroupedWyckoffPosition>>& permutation_grouped_Wyckoff_positions);
  bool arePermutationsComparableViaComposition(const xstructure& xstr); // DX20190624
  bool arePermutationsComparableViaComposition(const std::vector<uint>& composition, bool reduce_composition = false); // DX20190624
  bool arePermutationsComparableViaSymmetry(const std::vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions); // DX20190624

  // ---------------------------------------------------------------------------
  // ICSD comparisons
  std::string findICSDName(const std::string& name);
  std::string findMinimumICSDEntry(const std::vector<std::string>& ICSD_entries);

  // ---------------------------------------------------------------------------
  // matchable species/types
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species);
  bool sameSpecies(const xstructure& x1, const xstructure& x2, bool display);

  // ---------------------------------------------------------------------------
  // structure scaling
  void rescaleStructure(xstructure& x1, xstructure& x2);
  void atomicNumberDensity(const xstructure& xstr1, xstructure& xstr2);
  void atomicNumberDensity(const xstructure& xstr1, xstructure& xstr2, double& rescale_factor); // DX20201215
  void printParameters(const xstructure& xstr, std::ostream& oss);

  // ---------------------------------------------------------------------------
  // least-frequently occurring atom (LFA) functions
  bool sortBySecondPair(const std::pair<std::string, uint>& a, const std::pair<std::string, uint>& b);
  std::vector<std::string> sortSpeciesByFrequency(const xstructure& xstr);
  std::vector<uint> atomIndicesSortedByFrequency(const xstructure& xstr);
  bool similarLatticeParameters(const aurostd::xvector<double> d1, const aurostd::xmatrix<double> d2);

  // ---------------------------------------------------------------------------
  // atom mapping functions
  bool consistentAtomMappingType(const _atom& atom1, const _atom& atom2, uint index_x1, uint index_x2, bool same_species, bool is_collinear,
                                 bool is_non_collinear); // DX20201209
  bool consistentAtomMappingIndex(uint index1, uint index2, const std::vector<uint>& index1_list,
                                  const std::vector<uint>& index2_list); // DX20201209
  bool consistentAtomSetMappings(const std::string& atom1_name, const std::string& atom2_name, const std::vector<std::string>& vatoms1_name,
                                 const std::vector<std::string>& vatoms2_name); // DX20201209
  std::vector<aurostd::xvector<double>> minimizeMappingDistances(const std::vector<aurostd::xvector<double>>& distance_vectors); // DX20200909
  std::vector<aurostd::xvector<double>> minimizeMappingDistances(const std::vector<aurostd::xvector<double>>& distance_vectors, aurostd::xvector<double>& origin_shift); // DX20200909

  // ---------------------------------------------------------------------------
  // lattice similarity
  void cellDiagonal(const xstructure& xstr, std::vector<double>& diag_sum, std::vector<double>& diag_diff, const double& scale);
  void cellDiagonal(const aurostd::xmatrix<double>& lattice, std::vector<double>& diag_sum, std::vector<double>& diag_diff, const double& scale);

  // ---------------------------------------------------------------------------
  // environment analysis (near isoconfigurational analysis)
  void computeLFAEnvironments(std::vector<StructurePrototype>& structures, uint start_index = 0, uint end_index = AUROSTD_MAX_UINT); // DX20191122
  std::vector<AtomEnvironment> computeLFAEnvironment(const xstructure& xstr, bool unique_only = true); // DX20190711
  bool compatibleEnvironmentSets(const std::vector<AtomEnvironment>& env_set1, const std::vector<AtomEnvironment>& env_set2, bool same_species, bool ignore_environment_angles, bool exact_match); // DX20190711 //DX20200320 - added environment angles
  bool compatibleEnvironments(const AtomEnvironment& env_1,
                              const AtomEnvironment& env_2,
                              bool same_species,
                              bool ignore_environment_angles,
                              bool exact_match); // DX20190711 //DX20200320 - added environment angles
  bool compatibleEnvironments(const AtomEnvironment& env_1,
                              const AtomEnvironment& env_2,
                              std::vector<std::vector<std::string>>& matched_species,
                              bool same_species,
                              bool ignore_environment_angles,
                              bool exact_match); // DX20190711 //DX20200320 - added environment angles
  std::vector<std::vector<double>> getAnglesBetweenMixedSpeciesEnvironments(const std::vector<std::vector<aurostd::xvector<double>>>& neighbor_coordinates); // DX20190715
  bool compatibleNearestNeighborTypesEnvironments(const std::vector<std::vector<double>>& nn_lfa_with_types_1, const std::vector<std::vector<double>>& nn_lfa_with_types_2, int type_match);

  // ---------------------------------------------------------------------------
  // cost functions
  double latticeDeviation(const std::vector<double>& diag_sum1, const std::vector<double>& diag_sum2, const std::vector<double>& diag_diff1, const std::vector<double>& diag_diff2);
  void coordinateDeviation(structure_mapping_info& mapping_info, const std::vector<double>& nn_xstr1,
                           const std::vector<double>& nn_xstr2); // DX20201210
  void magneticDeviation(const xstructure& xstr1, const xstructure& xstr2,
                         structure_mapping_info& mapping_info); // DX20201210
  double computeMisfit(const structure_mapping_info& mapping_info); // DX20201210
  double computeMisfit(double dev, double dis, double fail);
  double computeMisfitMagnetic(const structure_mapping_info& mapping_info); // DX20201210
  double computeMisfitMagnetic(double dev, double dis, double fail, double mag_dis, double mag_fail); // DX20190801
  double checkLatticeDeviation(const double& xstr1_vol, const aurostd::xmatrix<double>& q2, const std::vector<double>& D1, const std::vector<double>& F1);

  // ---------------------------------------------------------------------------
  // structure manipulation
  xstructure GetLFASupercell(const xstructure& xstr, const aurostd::xvector<int>& dims, const std::string& lfa_name); // DX20190530
  void getLatticeTransformations(const aurostd::xmatrix<double>& lattice_original,
                                 const aurostd::xmatrix<double>& lattice_ideal,
                                 const std::vector<aurostd::xmatrix<double>>& candidate_lattices,
                                 std::vector<aurostd::xmatrix<double>>& basis_transformations,
                                 std::vector<aurostd::xmatrix<double>>& rotations); // DX20201015
  void getLatticeTransformation(const aurostd::xmatrix<double>& lattice_original,
                                const aurostd::xmatrix<double>& lattice_ideal,
                                const aurostd::xmatrix<double>& candidate_lattice,
                                aurostd::xmatrix<double>& basis_transformation,
                                aurostd::xmatrix<double>& rotation); // DX20201015
  std::vector<xstructure> getTransformedStructures(const xstructure& xstr,
                                                   const std::vector<aurostd::xmatrix<double>>& basis_transformations,
                                                   const std::vector<aurostd::xmatrix<double>>& rotations); // DX20201119

  xstructure supercell2newRepresentation(const xstructure& xstr_supercell, const aurostd::xmatrix<double>& lattice);

  void writePythonScript(std::ostream& oss); // DX20201228

} // namespace compare

#endif

// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
