// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalFinder (compare crystal structures) - Functions
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo

#include "aflow.h"
#include "aflow_pflow.h"

#ifndef __AFLOW_COMPARE_STRUCTURE_H_
#define __AFLOW_COMPARE_STRUCTURE_H_

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

#define _DEBUG_COMPARE_ false  //DX20201223

static string GENERAL_XTALFINDER_OPTIONS_LIST = "general_options: [--usage] [--misfit_match=0.1] [--misfit_match=0.2] [--np=16|--num_proc=16] [--optimize_match] [--no_scale_volume] [--ignore_symmetry] [--ignore_Wyckoff] [--ignore_environment] [--keep_unmatched!] [--match_to_aflow_prototypes!] [--magmom=<m1,m2,...|INCAR|OUTCAR>:...] [--add_aflow_prototype_designation] [--remove_duplicate_compounds] [--ICSD] [--print_mapping|--print] [--print=TEXT|JSON] [--quiet|--q] [--screen_only] [--primitivize] [--minkowski] [--niggli]";

//// ===== GroupedWyckoffPosition Class ===== //
//DX20191120 [MOVED TO aflow.h]

// ===== AtomEnvironment Class ===== //
//DX20191120 [MOVED TO aflow.h]

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
  xmatrix<double> rotation;                  // rotation for duplicate structure (rotate into representative reference frame)
  xmatrix<double> basis_transformation;      // basis transformation for duplicate structure (convert to lattice similar to reference)
  xvector<double> origin_shift;              // shift origin for duplicate structure (shift to common origin with representative structure)
  // mapping info                          
  vector<uint> atom_map;                     // atom index map: position in vector corresponds to reference index, value corresponds to duplicate structure index
  vector<uint> basis_map;                    // atom type map: position in vector corresponds to reference index, value corresponds to duplicate structure type 
  vector<double> distances_mapped;           // distances between mapped atoms between reference and duplciate structures
  vector<xvector<double> > vectors_mapped;   // mapping vectors between mapped atoms between reference and duplicate structures
};

// ***************************************************************************
// structure_mapping_info functions
// ***************************************************************************
namespace compare{
  structure_mapping_info initialize_misfit_struct(bool magnetic=false);
  void resetMisfitInfo(structure_mapping_info& str_mis, bool magnetic=false); //DX20220406
  string printAtomMappings(const structure_mapping_info& misfit_info);
  string printUnmatchedAtoms(const structure_mapping_info& misfit_info,const xstructure& xstr1,const xstructure& xstr2);
  void resetMappingInfo(structure_mapping_info& misfit_info); //DX20220406
  void resizeMappingInfo(structure_mapping_info& str_mis, uint size); //DX20220406
}

//DX20200225 - temp struct; working on more robust scheme
struct matching_structure {
  string name;
  double misfit;
};

// ***************************************************************************
// Struct: structure_container //DX20201218
// ***************************************************************************
struct structure_container {
  // intialization info
  string name;                                                  // name of structure
  string compound;                                              // compound name of structure
  vector<uint> stoichiometry;                                   // reduced stoichiometry of structure 
  uint natoms;                                                  // number of atoms in unit cell
  uint ntypes;                                                  // number of atom types
  vector<string> elements;                                      // element names in structure
  string source;                                                // indicating where structure came from, i.e., input, auid, aflow protos (for structure regeneration if necessary)
  bool is_structure_generated;                                  // boolean indicating if the structure has been generated
  uint relaxation_step;                                         // specifies the relaxation step of the structure (0=original, 1=one relaxation, 2=most relaxed)  
  // structure info
  xstructure structure;                                         // xstructure
  vector<double> nearest_neighbor_distances;                    // nearest neighbor distances for atoms, vector position corresponds to atom index 
  vector<AtomEnvironment> environments_LFA;                     // LFA atom environments (for near isoconfigurational comparison) 
  // symmetry info
  string Pearson;                                               // Pearson symbol
  uint space_group;                                             // space group
  vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;     // Wyckoff positions grouped by species/types
  uint number_compounds_matching_structure;                     // number of compounds that match to this structure
  // property info
  vector<string> properties_names;                              // API keyword for calculated properties (database comparisons only)
  vector<string> properties_units;                              // units for calculated properties (database comparisons only)
  vector<string> properties_types;                              // datatype for calculated properties (database comparison only)
  vector<string> properties;                                    // value(s) of calculated properties (database comparison only)
};

// ---------------------------------------------------------------------------
namespace compare{
  structure_container initializeStructureContainer();
  structure_container initializeStructureContainer(const xstructure& structure, bool same_species);
}

//DX 20201119
// ===== CompapareStructureContainers Class ===== //
// This class is necessary to pass an argument to the sorting function;
// the other solution is to use std::bind, but this is not available for early G++ versions
// see https://stackoverflow.com/questions/26444216/is-it-possible-to-use-stdsort-with-a-sort-function-that-takes-extra-arguments
class CompareStructureContainers{
  public:
    CompareStructureContainers(const vector<string>& sorting_attributes) : sorting_attributes(sorting_attributes){}
    bool operator()(const structure_container *a, const structure_container *b);
  private:
    vector<string> sorting_attributes;
};

// ***************************************************************************
// Class: StructurePrototype
// ***************************************************************************
// This class provides the unique prototype information, including the 
// symmetry, environments, and structures that exhibit this prototype 
// ---------------------------------------------------------------------------
class StructurePrototype{
  public:
    // ---------------------------------------------------------------------------
    // constructor/destructor/copy/assignment/operator<<
    StructurePrototype();                                                                   // constructor operator
    ~StructurePrototype();                                                                  // destructor operator
    void clear();                                                                           // clear
    friend ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype); // stringstream operator (printing)
    const StructurePrototype& operator=(const StructurePrototype& b);                       // assignment operator
    StructurePrototype(const StructurePrototype& b);                                        // copy constructor
    int iomode;                                                                             // mode for printing

    // ---------------------------------------------------------------------------
    // structure prototype information
    uint natoms;                                                                            // number of atoms in the prototype (from the representative structure; not necessarily reduced)
    int ntypes;                                                                             // number of types in prototype
    vector<string> elements;                                                                // list of elements exhibiting in this protoype (from representative and duplicate structures)
    vector<uint> stoichiometry;                                                             // reduced stoichiometry of prototype
    vector<vector<string> > atom_decorations_equivalent;                                    // equivalent atom decorations (permutations) of prototype
    string Pearson;                                                                         // Pearson symbol of prototype
    uint space_group;                                                                       // space group number of prototype
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;                               // Wyckoff positions grouped by site type
    vector<string> wyckoff_site_symmetry;                                                   // vector of Wyckoff site symmetries of prototype
    vector<int> wyckoff_multiplicity;                                                       // vector of Wyckoff multiplicities of prototype
    vector<int> wyckoff_letter;                                                             // vector of Wyckoff letters of prototype
    vector<AtomEnvironment> environments_LFA;                                               // vector of LFA atom environments
    string aflow_label;                                                                     // AFLOW label designation  
    vector<string> aflow_parameter_list;                                                    // vector of strings corresponding to AFLOW parameter variables 
    vector<double> aflow_parameter_values;                                                  // vector of doubles corresponding to AFLOW parameter values
    vector<string> matching_aflow_prototypes;                                               // vector of strings indicating matching AFLOW prototype labels

    // ---------------------------------------------------------------------------
    // structures
    structure_container *structure_representative;                                          // structural information for representative structure 
    vector<structure_container*> structures_duplicate;                                      // structural information for the duplicate structures    
    vector<structure_container*> structures_family;                                         // structural information for the same family structures
    vector<structure_mapping_info> mapping_info_duplicate;                                  // structural misfit and transformation information between the duplicate and representative structures
    vector<structure_mapping_info> mapping_info_family;                                     // structural misfit and transformation information between the family and representative structures

    // ---------------------------------------------------------------------------
    // properties stored in structure containers
    vector<string> property_names;                                                          // vector of property names (if using AFLUX)
    vector<string> property_units;                                                          // vector of property units (if using AFLUX)
    vector<string> property_types;                                                          // vector of property types (if using AFLUX) //DX20201230

    // ---------------------------------------------------------------------------
    // functions
    uint numberOfDuplicates() const; //DX20190506                                           // return the number of duplicate structures for this prototype (i.e., checks misfit value)
    string printRepresentativeStructure() const; //DX20201028                               // return the representative structure in a JSON format
    string printMatchedStructures(matched_structure_type_xtalfinder mode) const;            // return the matched structures in a JSON format //DX20201028
    string printPropertiesOfStructure(structure_container* str_pointer) const;              // return properties of structure in a JSON format
    string printStructureTransformationInformation(const structure_mapping_info& misfit_info) const; // return structure transformation information between structure and representative in JSON format
    uint numberOfComparisons(); //DX20181221                                                // return the number of comparisons for this prototype
    bool isSymmetryCalculated(); //DX20190228
    bool isLFAEnvironmentCalculated(); //DX20191105
    bool calculateSymmetry(); //DX20190118                                                  // calculate space group symmetry and populate symmetry attributes for StructurePrototype
    void putDuplicateAsFamily(uint index, bool keep_generated=false);                       // make duplicate structure a same family structure in the same StructurePrototypeObject //DX20190814 
    void copyPrototypeInformation(const StructurePrototype& b);                             // copy prototype information from one StructurePrototype object to another
    void copyDuplicate(const StructurePrototype& b, uint index, bool copy_misfit=false);    // copy duplicate info to another StructurePrototype object
    void removeNonDuplicate(uint index);                                                    // remove non-duplicate structures 
  private:
    void free();                                                                            // free operator
    void copy(const StructurePrototype& b);                                                 // copy constructor
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
    XtalFinderCalculator(uint num_proc_input=1, ostream& oss=cout);
    XtalFinderCalculator(ofstream& FileMESSAGE, uint num_proc_input=1, ostream& oss=cout);
    XtalFinderCalculator(double misfit_match_input,double misfit_family_input,uint num_proc_input=1,ostream& oss=cout);
    XtalFinderCalculator(double misfit_match_input,double misfit_family_input,ofstream& FileMESSAGE,uint num_proc_input=1,ostream& oss=cout);
    // ---------------------------------------------------------------------------
    // destructors
    ~XtalFinderCalculator();
    // ---------------------------------------------------------------------------
    // clear
    void clear();
    // ---------------------------------------------------------------------------
    // operator<<
    friend ostream& operator<<(ostream& oss, const XtalFinderCalculator& XtalFinderCalculator);
    // ---------------------------------------------------------------------------
    // copy/assignment
    const XtalFinderCalculator& operator=(const XtalFinderCalculator& b);
    XtalFinderCalculator(const XtalFinderCalculator& b);

    // ---------------------------------------------------------------------------
    // attributes
    double misfit_match;
    double misfit_family;
    uint num_proc;
    vector<structure_container> structure_containers;  // stores structures in a container (pointer for easy manipulation and mobility)

    // ---------------------------------------------------------------------------
    // compare methods
    void compareStructures(structure_container& str_rep,
        structure_container& str_matched, 
        structure_mapping_info& match_info, 
        bool same_species,
        bool scale_volume,
        bool optimize_match);
    vector<StructurePrototype> compareStructuresFromStructureList(const vector<string>& filenames, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); //DX20200103 - condensed bools to xoptions
    vector<StructurePrototype> compareStructuresFromDirectory(const string& directory, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); //DX20200103 - condensed bools to xoptions
    vector<StructurePrototype> compareStructuresFromString(const string& structures_string, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); //ME202010206
    vector<StructurePrototype> compareStructuresFromFile(const string& filename, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options); //DX20200103 - condensed bools to xoptions

    // ---------------------------------------------------------------------------
    // compare multiple structures
    vector<StructurePrototype> compareMultipleStructures(uint num_proc, bool same_species, const string& directory);
    vector<StructurePrototype> compareMultipleStructures(uint num_proc, bool same_species, const string& directory, const aurostd::xoption& comparison_options);

    // ---------------------------------------------------------------------------
    // compare2prototypes
    string printMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow);
    vector<StructurePrototype> compare2prototypes(const xstructure& xstrIN, const aurostd::xoption& vpflow); 
    void calculateMatchingAFLOWPrototypes(
        vector<StructurePrototype>& prototypes,
        uint num_proc);
    void getMatchingAFLOWPrototypes(
        uint i,
        vector<StructurePrototype>& prototypes,
        const aurostd::xoption& vpflow_protos);

    // ---------------------------------------------------------------------------
    // compare2database
    vector<StructurePrototype> compare2database(
        const xstructure& xstrIN, const aurostd::xoption& vpflow);

    // ---------------------------------------------------------------------------
    // compare permuations
    vector<string> getUniquePermutations(xstructure& xstr, uint num_proc=1, bool optimize_match=false);
    vector<string> getUniquePermutations(xstructure& xstr, aurostd::xoption& comparison_options, stringstream& results_ss, uint num_proc=1, bool print_misfit=false);
    void compareAtomDecorations(StructurePrototype& structure, uint num_proc, bool optimize_match);
    void compareAtomDecorations(StructurePrototype& structure, uint num_proc, aurostd::xoption& permutation_options);
    void compareAtomDecorations(StructurePrototype& structure, string& misfit_results, uint num_proc, aurostd::xoption& permutation_options);
    void generateAtomPermutedStructures(structure_container& structure);
    vector<string> getSpeciesPermutedStrings(const deque<uint>& stoichiometry);
    vector<string> getSpeciesPermutedStrings(const vector<uint>& stoichiometry);

    // ---------------------------------------------------------------------------
    // get command line options
    void getOptions(const aurostd::xoption& vpflow, aurostd::xoption& comparison_options);

    // ---------------------------------------------------------------------------
    // get command line options
    string getSpaceGroupMatchbookFromOptions(const aurostd::xoption& vpflow, uint relaxation_step); //DX20210615 - uint not bool

    // ---------------------------------------------------------------------------
    // add structures to container
    void addStructure2container(const xstructure& xstr,
        const string& structure_name,
        const string& source,
        uint relaxation_step,
        bool same_species);

    void addAFLOWPrototypes2container(const vector<string>& vlabel);

    void addDatabaseEntry2container(
        aflowlib::_aflowlib_entry& entry,
        const vector<string>& species,
        uint relaxation_step,
        bool same_species);

    // ---------------------------------------------------------------------------
    // remove methods
    void removeStructureFromContainerByName(const string& structure_name);

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
    void loadStructuresFromStructureList(const vector<string>& filenames, const vector<string>& magmoms_for_systems, bool same_species); //DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
    void loadStructuresFromDirectory(const string& directory, const vector<string>& magmoms_for_systems, bool same_species); //DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
    void loadStructuresFromFile(const string& filename, const vector<string>& magmoms_for_systems, bool same_species); //DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
    void loadStructuresFromStringstream(stringstream& input_stream, const vector<string>& magmoms_for_systems, bool same_species); //ME20210206
    void loadStructuresFromAflowlibEntries(const vector<aflowlib::_aflowlib_entry>& entries, const vector<string>& magmoms_for_systems, bool same_species); //DX20201201

    // ---------------------------------------------------------------------------
    // transform structures
    void performStructureConversions(
        uint i,  //ME20220207 - new xThread scheme
        const vector<bool>& calculate_primitive_vec,
        const vector<bool>& calculate_Minkowski_vec,
        const vector<bool>& calculate_Niggli_vec); //DX20210113
    void convertStructures(const aurostd::xoption& comparison_options, uint num_proc); //DX20201005
    void getPrimitiveStructures(uint start_index=0, uint end_index=AUROSTD_MAX_UINT); //DX20201005
    void getMinkowskiStructures(uint start_index=0, uint end_index=AUROSTD_MAX_UINT); //DX20201005
    void getNiggliStructures(uint start_index=0, uint end_index=AUROSTD_MAX_UINT); //DX20201005

    // ---------------------------------------------------------------------------
    // analyze symmetry
    bool isSymmetryCalculated(structure_container& structure);
    void calculateSymmetry(structure_container& str_rep);
    void calculateSymmetries(uint num_proc);
    void calculateSpaceGroups(uint start_index=0, uint end_index=AUROSTD_MAX_UINT, uint setting=0); //DX20191230 added setting
    void setSymmetryPlaceholders();
    void getSymmetryInfoFromXstructure(structure_container& str_rep); //DX20210104
    // ---------------------------------------------------------------------------
    // analyze environment
    bool isLFAEnvironmentCalculated(structure_container& structure);
    void computeLFAEnvironment(structure_container& str_rep, bool unique_only=true);
    void calculateLFAEnvironments(uint num_proc);
    void computeLFAEnvironments(uint start_index=0, uint end_index=AUROSTD_MAX_UINT);
    // ---------------------------------------------------------------------------
    // analyze neighbors 
    bool areNearestNeighborsCalculated(structure_container& structure);
    void getNearestNeighbors(uint num_proc);
    void calculateNearestNeighbors(uint start_index=0, uint end_index=AUROSTD_MAX_UINT);

    // ---------------------------------------------------------------------------
    // group structures
    vector<StructurePrototype> groupStructurePrototypes(
        bool same_species,
        bool ignore_symmetry,
        bool ignore_Wyckoff,
        bool ignore_environment,
        bool ignore_environment_angles,
        bool duplicates_removed,
        bool quiet=false); //DX20190731 - remove const and & //DX20190830 - added duplicates_removed //DX20200320 - added environment angles

    // ---------------------------------------------------------------------------
    // reorder structures
    void representativePrototypeForICSDRunsNEW(vector<StructurePrototype>& comparison_schemes);
    void makeRepresentativeEvenPermutation(vector<StructurePrototype>& comparison_schemes, const vector<string>& name_order);

    // ---------------------------------------------------------------------------
    // find duplicate compounds //DX20210114
    void findDuplicateCompounds(
        uint num_proc,
        bool remove_duplicates,
        const string& directory,
        const aurostd::xoption& comparison_options);

    // ---------------------------------------------------------------------------
    // check for better structure matches 
    vector<StructurePrototype> checkForBetterMatches(vector<StructurePrototype>& prototype_schemes, 
        uint num_proc, 
        bool check_for_better_matches, 
        bool same_species,
        const aurostd::xoption& comparison_options, 
        bool quiet);
    void combinePrototypesOfDifferentSymmetry(vector<StructurePrototype>& prototypes_final, bool same_species, uint num_proc);

    // ---------------------------------------------------------------------------
    // append unmatched structures into new groups
    void appendStructurePrototypes(
        vector<StructurePrototype>& comparison_schemes,
        vector<StructurePrototype>& prototypes_final,
        bool clean_unmatched, //DX20190506
        bool quiet);

    // ---------------------------------------------------------------------------
    // count matched/unmatched
    uint numberOfMismatches(const vector<StructurePrototype>& comparison_schemes);
    uint numberOfDuplicates(const StructurePrototype& prototype);

    // ---------------------------------------------------------------------------
    // split comparisons into threads (2D-array splitting)
    bool splitComparisonIntoThreads(
        vector<StructurePrototype>& comparison_schemes,
        uint num_proc,
        vector<std::pair<uint,uint> >& start_indices,
        vector<std::pair<uint,uint> >& end_indices);

    // ---------------------------------------------------------------------------
    // run multiple structures
    vector<StructurePrototype> runComparisonScheme(vector<StructurePrototype>& comparison_schemes,
        bool same_species, uint num_proc, const aurostd::xoption& comparison_options, 
        bool quiet=false); //DX20200103 - condensed bools to xoptions
    void runComparisons(
        vector<StructurePrototype>& comparison_schemes, 
        bool same_species, 
        bool scale_volume,
        bool optimize_match); 
    void runComparisonThreads(uint index, vector<StructurePrototype>& comparison_schemes,
        const vector<std::pair<uint,uint> >& vstart_indices,
        const vector<std::pair<uint,uint> >& vend_indices,
        bool same_species, 
        bool scale_volume, bool optimize_match); //DX20190822 - added comparison log bool


    // ---------------------------------------------------------------------------
    // get aflow label
    void calculatePrototypeDesignations(vector<StructurePrototype>& prototypes,uint num_proc);
    //ME20220207 - replaced with iterator for xThread
    void getPrototypeDesignations(vector<StructurePrototype>::iterator& prototypes);
    //void getPrototypeDesignations(vector<StructurePrototype>& prototypes);
    void getPrototypeDesignationsPreDistributed(
        vector<StructurePrototype>& prototypes,
        uint start_index=0,
        uint end_index=AUROSTD_MAX_UINT); //DX20191122

    // ---------------------------------------------------------------------------
    // print
    string printResults(
        const vector<StructurePrototype>& prototypes_final,
        bool same_species,
        filetype format);
    string printStructureMappingResults(
        const structure_mapping_info& misfit_info,
        const xstructure& xstr_reference,
        const xstructure& xstr_mapped,
        const string& mode="TEXT");
    string printAtomMappings(const structure_mapping_info& misfit_info,const xstructure& xstr1,const xstructure& xstr2);
    string printUnmatchedAtoms(const structure_mapping_info& misfit_info,const xstructure& xstr1,const xstructure& xstr2);

    // ---------------------------------------------------------------------------
    // lattice search
    void latticeSearch(
        structure_container& xstr_rep,
        structure_container& xstr_match,
        structure_mapping_info& match_info,
        bool same_species,
        bool optimize_match,
        bool scale_volume, //DX20200422
        uint num_proc); //DX20201123

    // ---------------------------------------------------------------------------
    // find translation vectors
    void findSimilarTranslationVectors(
        const xmatrix<double>& q1,
        const xstructure& xstr_LFA_supercell, 
        const xstructure& xstr,
        vector<xvector<double> >& lattice_vecs);

    // ---------------------------------------------------------------------------
    // similar lattices
    bool buildSimilarLattices(
        const vector<xvector<double> >& translation_vectors,
        const xmatrix<double>& q1,
        vector<xmatrix<double> >& lattices,
        vector<double>& latt_devs,
        bool optimize_match,
        bool scale_volume);

    // ---------------------------------------------------------------------------
    // search for atom mappings 
    bool searchAtomMappings(
        uint start_index, uint end_index,
        const xstructure& xstr1,
        const vector<double>& all_nn1,
        const xstructure& xstr2,
        const string& lfa,
        vector<xmatrix<double> >& lattices,
        vector<structure_mapping_info>& vstrs_matched,
        bool same_species,
        bool optimize_match);

    // ---------------------------------------------------------------------------
    // find matches (atoms) 
    bool findMatch(
        const xstructure& xstr1,
        const xstructure& xstr2,
        const vector<uint>& atom_indices_xstr1,
        const vector<uint>& atom_indices_xstr2,
        double minimum_interatomic_distance, //DX20200622
        structure_mapping_info& mapping_info,
        bool same_species);

    // ---------------------------------------------------------------------------
    // helper functions for external use
    vector<vector<uint> > groupSimilarXstructures(const vector<xstructure>& vxstrs, aurostd::xoption& comparison_options, bool same_species=true);
    vector<vector<uint> > groupSimilarXstructures(const vector<xstructure>& vxstrs, bool same_species=true);

    // ---------------------------------------------------------------------------
    // write output to file 
    void writeComparisonOutputFile(const string& output,
        const string& directory,
        filetype format,
        output_file_xtalfinder comparison_mode,
        bool same_species);


  private:
    void free();                              // free operator
    void copy(const XtalFinderCalculator& b); // copy constructor
};

namespace compare{
  // ---------------------------------------------------------------------------
  // pair-wise comparisons (for use by other AFLOW processes) 
  string compareInputStructures(aurostd::xoption& vpflow, std::istream& cin, ostream& logstream=cout); //ME20210206
  string compareInputStructures(const aurostd::xoption& vpflow, ostream& logstream=cout); //DX //DX20190425
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optimize_match, double& final_misfit, uint num_proc=1); //Main function //DX20191108 - remove const & from bools
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optimize_match, double& final_misfit, structure_mapping_info& final_misfit_info, uint num_proc=1); //Main function //DX20191108 - remove const & from bools 
  //bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, uint num_proc=1); //Overco, returns true (match), false (no match) //DX20191108 - remove const & from bools
  //bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optmize_match, uint num_proc=1);  //DX20191108 - remove const & from bools
  bool structuresMatch(const xstructure& xstr1, const xstructure& xstr2, bool same_species, uint num_proc=1); //Overco, returns true (match), false (no match) //DX20191108 - remove const & from bools
  bool structuresMatch(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optmize_match, uint num_proc=1);  //DX20191108 - remove const & from bools
  //double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match, uint num_proc=1); //Overloaded, returns misfit value //DX20191108 - remove const & from bools
  double getMisfitBetweenStructures(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match, uint num_proc=1); //Overloaded, returns misfit value //DX20191108 - remove const & from bools
  structure_mapping_info getTransformationBetweenStructures(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match, uint num_proc=1); //Overloaded, returns misfit value //DX20191108 - remove const & from bools

  // ---------------------------------------------------------------------------
  // comparisons to AFLOW database 
  vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout);  //CO20200225
  string printCompare2Database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225

  // ---------------------------------------------------------------------------
  // comparisons between entries in AFLOW database 
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20191125
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20191125

  // ---------------------------------------------------------------------------
  // comparisons aflowlib entries
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries, uint num_proc, bool same_species, bool scale_volume, bool optimize_match); //DX20201111
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries, ostream& oss, ofstream& FileMESSAGE, uint num_proc, bool same_species, bool scale_volume, bool optimize_match); //DX20201111

  // ---------------------------------------------------------------------------
  // comparisons to AFLOW prototype library 
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20181004 //DX20190314 - changed return value
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20181004 //DX20190314 - changed return value
  string printMatchingPrototypes(istream& cin, const aurostd::xoption& vpflow); //DX20190314
  vector<string> getMatchingPrototypes(xstructure& xstr, const string& catalog); //DX20190314
  vector<string> getMatchingPrototypes(xstructure& xstr, const string& catalog); //DX20190314
  vector<string> getMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow, const string& catalog); //DX20210421

  // ---------------------------------------------------------------------------
  // permutaion comparisons
  string compareAtomDecorations(istream& input, const aurostd::xoption& vpflow); //DX20181004

  // ---------------------------------------------------------------------------
  // isopointal AFLOW prototype functions 
  string isopointalPrototypes(istream& input, const aurostd::xoption& vpflow); //DX20200131
  vector<string> getIsopointalPrototypes(xstructure& xstr, const string& catalog); //DX20200131

  // ---------------------------------------------------------------------------
  // comparison options
  aurostd::xoption loadDefaultComparisonOptions(const string& mode=""); //DX20200103

  // ---------------------------------------------------------------------------
  // geneate structures 
  void generateStructures(vector<StructurePrototype>& structures, ostream& oss=cout, uint start_index=0, uint end_index=AUROSTD_MAX_UINT); //DX20191122
  bool generateStructure(const string& structure_name, const string& structure_source, uint relaxation_step, xstructure& structure, ostream& oss); //DX20200429 - added relaxation_step
  void removeNonGeneratedStructures(vector<StructurePrototype>& structures); //DX20191105

  // ---------------------------------------------------------------------------
  // functions for determining isopointal structures (same/compatible symmetry)
  void groupWyckoffPositions(const xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions);
  void groupWyckoffPositions(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, vector<GroupedWyckoffPosition>& grouped_positions); //DX20200512
  void groupWyckoffPositionsFromGroupedString(uint space_group_number, uint setting, const vector<vector<string> >& grouped_Wyckoff_string, vector<GroupedWyckoffPosition>& grouped_positions); //DX20200622 - removed pointer to uints
  string printWyckoffString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize=false);
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(const vector<GroupedWyckoffPosition>& grouped_Wyckoffs); //DX20190219
  bool matchableWyckoffPositions(
      const vector<GroupedWyckoffPosition>& grouped_Wyckoffs_str1,
      const vector<GroupedWyckoffPosition>& grouped_Wyckoffs_str2,
      bool same_species);
  bool matchableWyckoffPositionSet(const vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters,
      const vector<vector<string> >& grouped_Wyckoff_letters);
  vector<vector<string> > convertANRLWyckoffString2GroupedPositions(const string& label);
  vector<vector<string> > convertWyckoffString2GroupedPositions(const string& Wyckoff_letter_string);
  bool sameStoichiometry(const vector<uint>& stoich1, const vector<uint>& stoich2);
  bool matchableSpaceGroups(uint space_group_1, uint space_group_2);
  bool matchableEnantiomorphicSpaceGroups(uint space_group_1, uint space_group_2);
  bool filterPrototypes(uint& species_count, string& reduced_stoichiometry, uint& space_group,
      vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters,
      vector<string>& prototype_labels, vector<uint>& species_counts,
      vector<uint>& space_groups);
  bool structuresCompatible(const structure_container& structure1,
      const structure_container& structure2, bool same_species, 
      bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool ignore_environment_angles); //DX20201207

  // ---------------------------------------------------------------------------
  // comparing permutations/atom decorations 
  vector<std::pair<uint,uint> > calculateDivisors(uint number);
  bool checkNumberOfGroupings(const vector<StructurePrototype>& comparison_schemes, uint number);
  void generatePermutationString(const deque<uint>& stoichiometry, vector<string>& permutation); //DX20190508 //DX20191125 - changed from vector to deque
  void generatePermutationString(const vector<uint>& stoichiometry, vector<string>& permutation); //DX20190508
  bool generatePermutations(uint& num_elements, vector<uint>& indices, vector<string>& names, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, vector<vector<uint> >& permutations, vector<vector<string> >&name_order, vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions);
  bool arePermutationsComparableViaComposition(const xstructure& xstr); //DX20190624 
  bool arePermutationsComparableViaComposition(const vector<uint>& composition, bool reduce_composition=false); //DX20190624
  bool arePermutationsComparableViaSymmetry(const vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions); //DX20190624

  // ---------------------------------------------------------------------------
  // ICSD comparisons 
  string findICSDName(const string& name);
  string findMinimumICSDEntry(const vector<string>& ICSD_entries);

  // ---------------------------------------------------------------------------
  // matchable species/types 
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species);
  bool sameSpecies(const xstructure& x1, const xstructure& x2, bool display); 

  // ---------------------------------------------------------------------------
  // structure scaling
  void rescaleStructure(xstructure& x1, xstructure& x2);
  void atomicNumberDensity(const xstructure& xstr1, xstructure& xstr2);
  void atomicNumberDensity(const xstructure& xstr1, xstructure& xstr2, double& rescale_factor); //DX20201215
  void printParameters(const xstructure& xstr, ostream& oss);

  // ---------------------------------------------------------------------------
  // least-frequently occurring atom (LFA) functions
  bool sortBySecondPair(const std::pair<string,uint>& a, const std::pair<string,uint>& b);
  vector<string> sortSpeciesByFrequency(const xstructure& xstr);
  vector<uint> atomIndicesSortedByFrequency(const xstructure& xstr);
  bool similarLatticeParameters(const xvector<double> d1, const xmatrix<double> d2);

  // ---------------------------------------------------------------------------
  // atom mapping functions 
  bool consistentAtomMappingType(
      const _atom& atom1,
      const _atom& atom2,
      uint index_x1,
      uint index_x2,
      bool same_species,
      bool is_collinear,
      bool is_non_collinear); //DX20201209
  bool consistentAtomMappingIndex(
      uint index1,
      uint index2,
      const vector<uint>& index1_list,
      const vector<uint>& index2_list); //DX20201209
  bool consistentAtomSetMappings(
      const string& atom1_name,
      const string& atom2_name,
      const vector<string>& vatoms1_name,
      const vector<string>& vatoms2_name); //DX20201209
  vector<xvector<double> > minimizeMappingDistances(const vector<xvector<double> >& distance_vectors); //DX20200909
  vector<xvector<double> > minimizeMappingDistances(const vector<xvector<double> >& distance_vectors, xvector<double>& origin_shift); //DX20200909

  // ---------------------------------------------------------------------------
  // lattice similarity 
  void cellDiagonal(const xstructure& xstr, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  void cellDiagonal(const xmatrix<double>& lattice, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);

  // ---------------------------------------------------------------------------
  // environment analysis (near isoconfigurational analysis) 
  void computeLFAEnvironments(vector<StructurePrototype>& structures, uint start_index=0, uint end_index=AUROSTD_MAX_UINT);  //DX20191122
  vector<AtomEnvironment> computeLFAEnvironment(const xstructure& xstr, bool unique_only=true); //DX20190711
  bool compatibleEnvironmentSets(const vector<AtomEnvironment>& env_set1, 
      const vector<AtomEnvironment>& env_set2, bool same_species, bool ignore_environment_angles, bool exact_match); //DX20190711 //DX20200320 - added environment angles
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2,
      bool same_species,
      bool ignore_environment_angles,
      bool exact_match); //DX20190711 //DX20200320 - added environment angles
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2,
      vector<vector<string> > & matched_species, 
      bool same_species,
      bool ignore_environment_angles,
      bool exact_match); //DX20190711 //DX20200320 - added environment angles
  vector<vector<double> > getAnglesBetweenMixedSpeciesEnvironments(const vector<vector<xvector<double> > >& neighbor_coordinates); //DX20190715
  bool compatibleNearestNeighborTypesEnvironments(const vector<vector<double> >& nn_lfa_with_types_1,
      const vector<vector<double> >& nn_lfa_with_types_2,
      int type_match);

  // ---------------------------------------------------------------------------
  // cost functions
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
      const vector<double>& diag_diff1,const vector<double>& diag_diff2); 
  void coordinateDeviation(
      structure_mapping_info& mapping_info,
      const vector<double>& nn_xstr1,
      const vector<double>& nn_xstr2); //DX20201210
  void magneticDeviation(
      const xstructure& xstr1, const xstructure& xstr2, 
      structure_mapping_info& mapping_info); //DX20201210
  double computeMisfit(const structure_mapping_info& mapping_info); //DX20201210
  double computeMisfit(double dev, double dis, double fail);
  double computeMisfitMagnetic(const structure_mapping_info& mapping_info); //DX20201210
  double computeMisfitMagnetic(double dev,double dis,double fail,double mag_dis,double mag_fail); //DX20190801
  double checkLatticeDeviation(const double& xstr1_vol, const xmatrix<double>& q2, const vector<double>& D1, const vector<double>& F1);

  // ---------------------------------------------------------------------------
  // structure manipulation 
  xstructure GetLFASupercell(const xstructure& xstr, const xvector<int>& dims, const string& lfa_name); //DX20190530
  void getLatticeTransformations(const xmatrix<double>& lattice_original, 
      const xmatrix<double>& lattice_ideal,
      const vector<xmatrix<double> >& candidate_lattices,
      vector<xmatrix<double> >& basis_transformations,
      vector<xmatrix<double> >& rotations); //DX20201015
  void getLatticeTransformation(const xmatrix<double>& lattice_original, 
      const xmatrix<double>& lattice_ideal,
      const xmatrix<double>& candidate_lattice,
      xmatrix<double>& basis_transformation,
      xmatrix<double>& rotation); //DX20201015
  vector<xstructure> getTransformedStructures(
      const xstructure& xstr,
      const vector<xmatrix<double> >& basis_transformations,
      const vector<xmatrix<double> >& rotations); //DX20201119


  xstructure supercell2newRepresentation(const xstructure& xstr_supercell, const xmatrix<double>& lattice);

  void writePythonScript(ostream& oss); //DX20201228

} // end namespace

#endif 

// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
