// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalFinder (identify prototypes and compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
#include "aflow_symmetry_spacegroup.h"

//[OBSOLETE - ME20220128] #undef AFLOW_MULTITHREADS_ENABLE

//[OBSOLETE - ME20220128] #if GCC_VERSION >= 40400   // added two zeros
//[OBSOLETE - ME20220128] #define AFLOW_MULTITHREADS_ENABLE 1
//[OBSOLETE - ME20220128] #include <thread>
//[OBSOLETE - ME20220128] #include <mutex>
//[OBSOLETE - ME20220128] static std::mutex _mutex_;
//[OBSOLETE - ME20220128] #else
//[OBSOLETE - ME20220128] #warning "The multithread parts of AFLOW-XtalFinder will be not included, since they need gcc 4.4 and higher (C++0x support)."
//[OBSOLETE - ME20220128] #endif

//[OBSOLETE - ME20220128] // for multi-threads on-the-fly scheme (explanation in AAPL/aflow_aapl_tcond.cpp, developed by M. Esters (ME))
//[OBSOLETE - ME20220128] static int task_counter = 0;

//ME20220207 - Changed all functions to use xThread

using namespace std::placeholders;  //ME20220207

// ***************************************************************************
// XtalFinderCalculator::compareStructures() - MAIN COMPARISON FUNCTION
// ***************************************************************************
void XtalFinderCalculator::compareStructures(
    structure_container& str_rep,
    structure_container& str_matched,
    structure_mapping_info& match_info,
    bool same_species,
    bool scale_volume,
    bool optimize_match) {

  // This is the main comparison function that analyzes similarity between
  // two crystal structures

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  // ---------------------------------------------------------------------------
  // determine minimum interatomic distances of structures (resolution of atoms) //DX20200623
  // //DX20200715 - may need to rescale this if the structures are being scaled later....
  if(str_rep.structure.dist_nn_min==AUROSTD_NAN){
    if(str_rep.nearest_neighbor_distances.size() > 0){ str_rep.structure.dist_nn_min=aurostd::min(str_rep.nearest_neighbor_distances); }
    else{ str_rep.structure.dist_nn_min=SYM::minimumDistance(str_rep.structure); }
  }
  if(str_matched.structure.dist_nn_min==AUROSTD_NAN){
    if(str_matched.nearest_neighbor_distances.size() > 0){ str_matched.structure.dist_nn_min=aurostd::min(str_matched.nearest_neighbor_distances); }
    else{ str_matched.structure.dist_nn_min=SYM::minimumDistance(str_matched.structure); }
  }

  if(compare::matchableSpecies(str_rep.structure,str_matched.structure,same_species)){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Searching for new representation of test structure ..."<<endl;}
    latticeSearch(str_rep,str_matched,match_info,same_species,optimize_match,scale_volume,num_proc); //DX20190530 //DX20200422 - scale_volume added
  }
}

// ***************************************************************************
// XtalFinderCalculator::compareMultipleStructures() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareMultipleStructures(
    uint num_proc, bool same_species,
    const string& directory){

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
  if(!same_species){
    comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
  }

  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

vector<StructurePrototype> XtalFinderCalculator::compareMultipleStructures(
    uint num_proc,
    bool same_species,
    const string& directory,
    const aurostd::xoption& comparison_options){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;
  bool quiet = false;

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}  //CO20200508

  message << "Total number of structures to compare: " << structure_containers.size();
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  if(structure_containers.size() == 0){
    message << "No structures to compare.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  // ---------------------------------------------------------------------------
  // convert to structures to certain representations //DX20201006
  // conversion type(s) is indicated in the comparison_options flag
  if(comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE") ||
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI") ||
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI")){
    convertStructures(comparison_options,num_proc);
  }

  // ---------------------------------------------------------------------------
  // calculate symmetries of structures
  // if already calculated, do not recalculate
  bool all_symmetries_calculated = true;
  for(uint i=0;i<structure_containers.size();i++){ all_symmetries_calculated = (all_symmetries_calculated&&isSymmetryCalculated(structure_containers[i])); } //DX20200810

  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && !all_symmetries_calculated){
    calculateSymmetries(num_proc);
  }
  else if(!all_symmetries_calculated){
    setSymmetryPlaceholders();
  }

  // ---------------------------------------------------------------------------
  // calculate LFA environments of database entries
  // if already calculated, do not recalculate
  bool all_environments_calculated = true;
  for(uint i=0;i<structure_containers.size();i++){ all_environments_calculated = (all_environments_calculated&&isLFAEnvironmentCalculated(structure_containers[i])); } //DX20200810
  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS") && !all_environments_calculated){
    calculateLFAEnvironments(num_proc);
  }

  // ---------------------------------------------------------------------------
  // get nearest neighbor info, can perhaps only calculate this if necessary
  // (i.e., if we know we will perform a comparison)
  // if already calculated, do not recalculate //DX20201201
  bool all_neighbors_calculated = true;
  for(uint i=0;i<structure_containers.size();i++){ all_neighbors_calculated = (all_neighbors_calculated&&(structure_containers[i].nearest_neighbor_distances.size()!=0)); } //DX20200925 - gcc-10 warnings
  if(!all_neighbors_calculated){
    getNearestNeighbors(num_proc);
  }

  // ---------------------------------------------------------------------------
  // remove duplicate compounds first; recursion of this function
  if(comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS")){
    findDuplicateCompounds(num_proc, true, directory, comparison_options); //true: remove duplicates
  }

  // ---------------------------------------------------------------------------
  // group structures based on stoichiometry and symmetry (unless ignoring symmetry/Wyckoff)
  vector<StructurePrototype> comparison_schemes = groupStructurePrototypes(
      same_species,
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
      false,
      quiet); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // if ICSD comparison, make structure with minimum ICSD number the representative structure
  if(comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON")){
    representativePrototypeForICSDRunsNEW(comparison_schemes);
  }

  // ---------------------------------------------------------------------------
  // compare structures
  vector<StructurePrototype> prototypes_final = runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, quiet); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // return if there are no similar structures
  if(prototypes_final.size()==0){ return prototypes_final; }
  comparison_schemes.clear();

  // ---------------------------------------------------------------------------
  // combine prototypes regardless of having different space groups
  // BETA TESTING (perhaps we shouldn't do this) - compare::checkPrototypes(num_proc,same_species,prototypes_final);

  message << "Number of unique prototypes: " << prototypes_final.size() << " (out of " << structure_containers.size() << " structures).";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);

  // ---------------------------------------------------------------------------
  // get unique atom decorations prototype (representative) structures
  if(!same_species && comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS")){

    message << "Determining the unique atom decorations for each prototype.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    // find unique atom decorations of prototype
    for(uint i=0;i<prototypes_final.size();i++){
      XtalFinderCalculator xtal_finder_permutations;
      xtal_finder_permutations.compareAtomDecorations(prototypes_final[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
    }
  }

  if(comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION")){
    calculatePrototypeDesignations(prototypes_final,num_proc);
  }

  // ---------------------------------------------------------------------------
  // for testing/development; in case the subsequent analyses fail, checkpoint file
  bool store_checkpoint=false;
  if(store_checkpoint){
    string results_json = printResults(prototypes_final, same_species, json_ft);
    writeComparisonOutputFile(results_json, directory, json_ft, compare_input_xf, same_species);
    string results_txt = printResults(prototypes_final, same_species, txt_ft);
    writeComparisonOutputFile(results_txt, directory, txt_ft, compare_input_xf, same_species);
  }

  if(comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS")){
    calculateMatchingAFLOWPrototypes(prototypes_final,num_proc);
  }

  return prototypes_final;
}



// ***************************************************************************
// XtalFinderCalculator::getOptions() //DX20201201
// ***************************************************************************
void XtalFinderCalculator::getOptions(
    const aurostd::xoption& vpflow,
    aurostd::xoption& comparison_options){

  // Get options from vpflow (i.e., command-line)
  // Contains the options common to nearly all comparison functions

  stringstream message;

  // ---------------------------------------------------------------------------
  // FLAG: misfit threshold //DX20201119
  if(vpflow.flag("COMPARE::MISFIT_MATCH")) {
    misfit_match = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE::MISFIT_MATCH"));
    comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE::MISFIT_MATCH"));
  }
  if(vpflow.flag("COMPARE::MISFIT_FAMILY")) {
    misfit_family = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE::MISFIT_FAMILY"));
    comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE::MISFIT_FAMILY"));
  }
  // match threshold must be less than family threshold
  if(misfit_match>misfit_family){
    message << "Matching misfit threshold must be less than the same family threshold:"
      << " misfit match threshold: " << misfit_match
      << " misfit family threshold: " << misfit_family;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
  }
  message << "Misfit threshold for matched structures: " << misfit_match << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  message << "Misfit threshold for structures in the same family: " << misfit_family << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  // ---------------------------------------------------------------------------
  // FLAG: number of processors (multithreading)
  if(vpflow.flag("COMPARE::NP")) {
    num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE::NP"));
  }

  // ---------------------------------------------------------------------------
  // FLAG: optimize match
  if(vpflow.flag("COMPARE::OPTIMIZE_MATCH")) {
    comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",TRUE);
    message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons).";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: no volume scaling
  if(vpflow.flag("COMPARE::NO_SCALE_VOLUME")) {
    comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
    message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: do not remove unmatched structures from the StructurePrototype Object
  // keeps results of each comparison
  if(vpflow.flag("COMPARE::KEEP_UNMATCHED")) {
    comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
    message << "OPTIONS: Keep unmatched structure comparisons (reveals mappings despite being above the misfit_match threshold).";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: ignore Wyckoff positions
  if(vpflow.flag("COMPARE::IGNORE_WYCKOFF")) {
    comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
    message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons).";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: ignore symmetry
  if(vpflow.flag("COMPARE::IGNORE_SYMMETRY")) {
    comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
    comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
    message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons).";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: ignore LFA environment analysis
  if(vpflow.flag("COMPARE::IGNORE_ENVIRONMENT_ANALYSIS")) {
    comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
    message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: remove duplicate compounds (useful for non-biased statistics)
  if(vpflow.flag("COMPARE::REMOVE_DUPLICATE_COMPOUNDS")) {
    comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
    message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: match unique structures to the AFLOW prototypes
  if(vpflow.flag("COMPARE::MATCH_TO_AFLOW_PROTOS")) {
    comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",TRUE);
    message << "OPTIONS: Compare unique structures to the AFLOW prototypes.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: get AFLOW ANRL designation for unique structures
  if(vpflow.flag("COMPARE::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
    comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",TRUE);
    message << "OPTIONS: Cast unique structures into AFLOW standard designation.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: UNDECORATED
  if(vpflow.flag("COMPARE::UNDECORATED_COMPARISON")) {
    comparison_options.flag("COMPARISON_OPTIONS::UNDECORATED_COMPARISON",TRUE);
    message << "OPTIONS: Undecorated comparison; compare structures as if they had one atom-type.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: primitivize structures //DX20201005
  if(vpflow.flag("COMPARE::PRIMITIVIZE")) {
    comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",TRUE);
    message << "OPTIONS: Converting all structures to a primitive representation.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: Minkowski reduction //DX20201005
  if(vpflow.flag("COMPARE::MINKOWSKI")) {
    comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",TRUE);
    message << "OPTIONS: Performing Minkowski lattice reduction on all structures.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: Niggli reduction //DX20201005
  if(vpflow.flag("COMPARE::NIGGLI")) {
    comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",TRUE);
    message << "OPTIONS: Performing Niggli lattice reduction on all structures.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: ICSD comparison - structure with minimum ICSD number as representative prototype
  // in general: smaller ICSD number = older = more reliable
  if(vpflow.flag("COMPARE::ICSD_COMPARISON")) {
    comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",TRUE);
    message << "OPTIONS: Running on ICSD structures; use oldest ICSD number as representative prototype.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // FLAG: do not calculate unique atom decorations
  if(vpflow.flag("COMPARE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
    comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
    message << "OPTIONS: Do not calculate unique atom decorations.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

}

// ***************************************************************************
// XtalFinderCalculator::getSpaceGroupMatchbookFromOptions() //DX20210112
// ***************************************************************************
string XtalFinderCalculator::getSpaceGroupMatchbookFromOptions(
    const aurostd::xoption& vpflow,
    uint relaxation_step){ //DX20210615 - uint not bool

  stringstream message;

  // ---------------------------------------------------------------------------
  // FLAG: specify space group
  vector<uint> vspace_groups_uint;
  if(vpflow.flag("COMPARE_DATABASE_ENTRIES::SPACE_GROUP")){

    string space_group_input = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::SPACE_GROUP");
    message << "OPTIONS: Requesting the following space groups: " << space_group_input;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<string> vspace_group_strings;
    aurostd::string2tokens(space_group_input,vspace_group_strings,",");

    uint sg_tmp = 0;
    for(uint i=0;i<vspace_group_strings.size();i++){
      sg_tmp = aurostd::string2utype<uint>(vspace_group_strings[i]);
      if(sg_tmp<1 || sg_tmp>230){
        message << "Invalid space group requested: " << vspace_group_strings[i] << ". Please check input.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
      }
      vspace_groups_uint.push_back(sg_tmp);
    }
  }

  // if nothing, return empty matchbook
  if(vspace_groups_uint.size()==0){ return ""; }

  // ---------------------------------------------------------------------------
  // AFLUX matchbook preparations: get space group(s)
  return aflowlib::getSpaceGroupMatchbook(vspace_groups_uint, relaxation_step);
}

// ***************************************************************************
// initializeStructureContainer() //DX20201115
// ***************************************************************************
namespace compare {
  structure_container initializeStructureContainer(
      const xstructure& structure,
      bool same_species){

    structure_container str_rep;
    str_rep.structure = structure;

    // ---------------------------------------------------------------------------
    // pre-condition structure
    str_rep.structure.FixLattices();
    str_rep.structure.ReScale(1.0);
    str_rep.structure.BringInCell();

    // ---------------------------------------------------------------------------
    // set element/stoichiometry attributes
    str_rep.stoichiometry = str_rep.structure.GetReducedComposition(!same_species);
    if(str_rep.structure.atoms.size() && str_rep.structure.atoms[0].name.empty()){str_rep.structure.DecorateWithFakeElements();} // GetElements does not change xstructure //SD20220221
    str_rep.elements = str_rep.structure.GetElements(true,true); // true: clean names and assign fake names
    str_rep.compound = pflow::prettyPrintCompound(str_rep.elements,str_rep.stoichiometry,no_vrt,false,txt_ft); //remove ones is true  //DX20190311 //DX20190313 - use xstr //eventually redundant
    // update xstructure species
    if(str_rep.structure.species.size()==0){
      deque<string> deque_species; for(uint j=0;j<str_rep.elements.size();j++){deque_species.push_back(str_rep.elements[j]);}
      str_rep.structure.SetSpecies(deque_species);
      str_rep.structure.SpeciesPutAlphabetic();
    }
    // clean species
    else{
      for(uint s=0;s<str_rep.structure.species.size();s++){str_rep.structure.species[s]=KBIN::VASP_PseudoPotential_CleanName(str_rep.structure.species[s]); } //DX20190711
      str_rep.structure.SetSpecies(str_rep.structure.species);
    }

    // ---------------------------------------------------------------------------
    // initialize attributes
    str_rep.name = "";
    str_rep.natoms = str_rep.structure.atoms.size();
    str_rep.ntypes = str_rep.structure.num_each_type.size();
    str_rep.number_compounds_matching_structure = 0;
    str_rep.space_group = 0; //DX20210126

    vector<double> zero_vec_double;
    str_rep.nearest_neighbor_distances = zero_vec_double;
    return str_rep;
  }
}

// ***************************************************************************
// StructurePrototype Class
// ***************************************************************************

// ***************************************************************************
// StructurePrototype::StructurePrototype() (constructor)
// ***************************************************************************
StructurePrototype::StructurePrototype(){
  free();
}

// ***************************************************************************
// StructurePrototype::free()
// ***************************************************************************
void StructurePrototype::free(){
  iomode=JSON_MODE;
  natoms=0;
  ntypes=0;
  elements.clear();
  stoichiometry.clear();
  atom_decorations_equivalent.clear();
  Pearson="";
  space_group=0;
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  aflow_label=""; //DX20190724
  aflow_parameter_list.clear(); //DX20190724
  aflow_parameter_values.clear(); //DX20190724
  matching_aflow_prototypes.clear(); //DX20190724
  environments_LFA.clear(); //DX20190711
  structure_representative = NULL; //DX20201204
  structures_duplicate.clear(); //DX20201207
  structures_family.clear(); //DX20201207
  mapping_info_duplicate.clear(); //DX20191217
  mapping_info_family.clear(); //DX20191217
  property_names.clear();
  property_units.clear();
}

// ***************************************************************************
// StructurePrototype::~StructurePrototype() (destructor)
// ***************************************************************************
StructurePrototype::~StructurePrototype(){
  free();
}

// ***************************************************************************
// StructurePrototype::clear() (clear)
// ***************************************************************************
void StructurePrototype::clear(){
  free();
}

// ***************************************************************************
// StructurePrototype::StructurePrototype(...) (copy constructor)
// ***************************************************************************
StructurePrototype::StructurePrototype(const StructurePrototype& b){
  if(this!=&b){
    copy(b);
  }
}

// ***************************************************************************
// StructurePrototype::copy()
// ***************************************************************************
void StructurePrototype::copy(const StructurePrototype& b) {
  iomode=b.iomode;
  ntypes=b.ntypes;
  elements=b.elements;
  stoichiometry=b.stoichiometry;
  natoms=b.natoms;
  atom_decorations_equivalent=b.atom_decorations_equivalent;
  Pearson=b.Pearson;
  space_group=b.space_group;
  grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
  wyckoff_site_symmetry=b.wyckoff_site_symmetry;
  wyckoff_multiplicity=b.wyckoff_multiplicity;
  wyckoff_letter=b.wyckoff_letter;
  aflow_label=b.aflow_label; //DX20190724
  aflow_parameter_list=b.aflow_parameter_list; //DX20190724
  aflow_parameter_values=b.aflow_parameter_values; //DX20190724
  matching_aflow_prototypes=b.matching_aflow_prototypes; //DX20190724
  environments_LFA=b.environments_LFA; //DX20190711
  structure_representative=b.structure_representative; //DX20201204
  structures_duplicate=b.structures_duplicate; //DX20201207
  structures_family=b.structures_family; //DX20201207
  mapping_info_duplicate=b.mapping_info_duplicate; //DX20191217
  mapping_info_family=b.mapping_info_family; //DX20191217
  property_names=b.property_names;
  property_units=b.property_units;
}

// ***************************************************************************
// StructurePrototype::operator= (assignment)
// ***************************************************************************
const StructurePrototype& StructurePrototype::operator=(const StructurePrototype& b){
  if(this!=&b){
    free();
    copy(b);
  }
  return *this;
}

// ***************************************************************************
// StructurePrototype::operator<< (output/print)
// ***************************************************************************
ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype){

  // operator<< for the StruturePrototype object (looks like a JSON)

  //if(StructurePrototype.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  StructurePrototype.iomode=JSON_MODE;
  //}

  if(StructurePrototype.iomode==JSON_MODE){
    string eendl="";
    bool roff=true; //round off
    stringstream sscontent_json;
    vector<string> vcontent_json, tmp_vstring;

    // structure_representative
    sscontent_json << "\"structure_representative\":" << StructurePrototype.printRepresentativeStructure() << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // ntypes
    sscontent_json << "\"ntypes\":" << StructurePrototype.ntypes << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // elements
    sscontent_json << "\"elements\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.elements,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // stoichiometry
    sscontent_json << "\"stoichiometry\":[" << aurostd::joinWDelimiter(StructurePrototype.stoichiometry,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // natoms
    sscontent_json << "\"natoms\":" << StructurePrototype.natoms << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // atom_decorations_equivalent
    if(StructurePrototype.atom_decorations_equivalent.size()!=0){ //DX20190425 - only print if calculated
      sscontent_json << "\"atom_decorations_equivalent\":[";
      tmp_vstring.clear();
      for(uint i=0;i<StructurePrototype.atom_decorations_equivalent.size();i++){
        tmp_vstring.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.atom_decorations_equivalent[i],"\""),",")+"]");
      }
      sscontent_json << aurostd::joinWDelimiter(tmp_vstring,",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    } //DX20190425

    //DX20190425 [OBSOLETE] // Pearson
    //DX20190425 [OBSOLETE] sscontent_json << "\"Pearson\":\"" << StructurePrototype.Pearson << "\"" << eendl;
    //DX20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // space_group
    sscontent_json << "\"space_group\":" << StructurePrototype.space_group << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // grouped_Wyckoff_positions
    sscontent_json << "\"grouped_Wyckoff_positions\":[";
    tmp_vstring.clear();
    for(uint i=0;i<StructurePrototype.grouped_Wyckoff_positions.size();i++){
      stringstream ss_tmp;
      ss_tmp << StructurePrototype.grouped_Wyckoff_positions[i];
      tmp_vstring.push_back(ss_tmp.str());
    }
    sscontent_json << aurostd::joinWDelimiter(tmp_vstring,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    //DX20190425 [OBSOLETE] // Wyckoff_multiplicities
    //DX20190425 [OBSOLETE] sscontent_json << "\"Wyckoff_multiplicities\":[" << aurostd::joinWDelimiter(StructurePrototype.wyckoff_multiplicity,",") << "]" << eendl;
    //DX20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    //DX20190425 [OBSOLETE] // Wyckoff_site_symmetries
    //DX20190425 [OBSOLETE] sscontent_json << "\"Wyckoff_site_symmetries\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.wyckoff_site_symmetry,"\""),",") << "]" << eendl;
    //DX20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    //DX20190425 [OBSOLETE] // Wyckoff_letters
    //DX20190425 [OBSOLETE] sscontent_json << "\"Wyckoff_letters\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.wyckoff_site_symmetry,"\""),",") << "]" << eendl;
    //DX20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    if(StructurePrototype.aflow_label.size()!=0){
      // aflow_label
      sscontent_json << "\"aflow_prototype_label\":\"" << StructurePrototype.aflow_label << "\"" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // aflow_parameter_list
      sscontent_json << "\"aflow_prototype_params_list\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.aflow_parameter_list,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // aflow_parameter_values
      sscontent_json << "\"aflow_prototype_params_values\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(StructurePrototype.aflow_parameter_values,8,roff),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }

    // matching_aflow_prototypes
    if(StructurePrototype.matching_aflow_prototypes.size()!=0){
      sscontent_json << "\"matching_aflow_prototypes\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.matching_aflow_prototypes,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }

    // environments_LFA
    sscontent_json << "\"environments_LFA\":[";
    tmp_vstring.clear();
    for(uint i=0;i<StructurePrototype.environments_LFA.size();i++){
      stringstream ss_tmp;
      //DX20210624 [OBSOLETE] ss_tmp << StructurePrototype.environments_LFA[i];
      ss_tmp << StructurePrototype.environments_LFA[i].toJSON().toString(); //DX20210624
      tmp_vstring.push_back(ss_tmp.str());
    }
    sscontent_json << aurostd::joinWDelimiter(tmp_vstring,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // duplicate info //DX20201218 - put into separate function
    sscontent_json << "\"structures_duplicate\":[" << StructurePrototype.printMatchedStructures(duplicate_structure_xf) << "]";
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // family info  //DX20201218 - put into separate function
    sscontent_json << "\"structures_family\":[" << StructurePrototype.printMatchedStructures(family_structure_xf) << "]";
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    if(StructurePrototype.property_names.size()!=0){ //DX20190425 - only print if requested
      // property_names
      sscontent_json << "\"property_names\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.property_names,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // property_units
      sscontent_json << "\"property_units\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.property_units,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    } //DX20190425

    // Put into json StructurePrototype object
    oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
    vcontent_json.clear();
  }
  return oss;
}

// ***************************************************************************
// StructurePrototype::numberOfDuplicates() //DX20201115
// ***************************************************************************
uint XtalFinderCalculator::numberOfDuplicates(const StructurePrototype& prototype){

  // Count the number of duplicate structures

  uint number_of_duplicates = 0;
  for(uint i=0;i<prototype.mapping_info_duplicate.size();i++){
    if(prototype.mapping_info_duplicate[i].misfit<=misfit_match &&
        (prototype.mapping_info_duplicate[i].misfit+1.0)>1e-3 ){
      number_of_duplicates+=1;
    }
  }

  return number_of_duplicates;
}

// ***************************************************************************
// StructurePrototype::printPropertiesOfStructure() //DX20201230
// ***************************************************************************
string StructurePrototype::printPropertiesOfStructure(
    structure_container* str_pointer) const {

  // Print the properties of the structure (JSON format)

  vector<string> tokens;
  aurostd::JSONwriter json;

  if(str_pointer->properties.size()>0){
    for(uint i=0;i<str_pointer->properties.size();i++){
      if(str_pointer->properties_types[i] == "number"){ json.addNumber(str_pointer->properties_names[i], str_pointer->properties[i]); } 
      else if(str_pointer->properties_types[i] == "numbers" || str_pointer->properties_types[i] == "strings"){ 
        aurostd::string2tokens(str_pointer->properties[i],tokens,",");
        json.addVector(str_pointer->properties_names[i], tokens); 
      } 
      else { json.addString(str_pointer->properties_names[i], str_pointer->properties[i]); }  // when all else fails: treat as string
    }
  }

  return json.toString(false);

}

// ***************************************************************************
// StructurePrototype::printStructureTransformationInformation() //DX20210115
// ***************************************************************************
string StructurePrototype::printStructureTransformationInformation(
    const structure_mapping_info& misfit_info) const {

  // Print the structure transformation information
  // between this structure and the representative structure (JSON format)

  vector<string> tokens;
  aurostd::JSONwriter json;

  bool roff = true;
  // basis transformation
  json.addMatrix("basis_transformation", misfit_info.basis_transformation, 5, roff);

  // rotation
  json.addMatrix("rotation", misfit_info.rotation, 5, roff);

  // origin_shift
  json.addVector("origin_shift", misfit_info.origin_shift, 5, roff);

  // atom map
  json.addVector("atom_map", misfit_info.atom_map);

  // basis map
  json.addVector("basis_map", misfit_info.basis_map);

  return json.toString(false);

}

// ***************************************************************************
// StructurePrototype::printRepresentativeStructure() //DX20201115
// ***************************************************************************
string StructurePrototype::printRepresentativeStructure() const {

  // Print the representative structure information (JSON format)

  string eendl="";
  stringstream sscontent_json;
  vector<string> vcontent_json;

  sscontent_json << "\"name\":" << "\"" << structure_representative->name << "\"" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  sscontent_json << "\"number_compounds_matching_structure\":" << structure_representative->number_compounds_matching_structure << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // representative structure may not have properties
  if(structure_representative->properties.size()>0){
    vcontent_json.push_back(printPropertiesOfStructure(structure_representative));
  }

  return "{" + aurostd::joinWDelimiter(vcontent_json,",") + "}";

}

// ***************************************************************************
// StructurePrototype::printMatchedStructure() //DX20201115
// ***************************************************************************
string StructurePrototype::printMatchedStructures(matched_structure_type_xtalfinder mode) const {

  // Print matched structure information (JSON format)
  // Generalized for duplicate or same family structure

  bool roff=true;
  string eendl="";
  stringstream sscontent_json;
  vector<string> vcontent_json, vstructures;

  bool print_transformation = true;

  // print duplicate structure information
  if(mode == duplicate_structure_xf){
    for(uint j=0;j<structures_duplicate.size();j++){
      sscontent_json << "\"name\":" << "\"" << structures_duplicate[j]->name << "\"" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"misfit\":" << aurostd::utype2string<double>(mapping_info_duplicate[j].misfit,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"lattice_deviation\":" << aurostd::utype2string<double>(mapping_info_duplicate[j].lattice_deviation,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"coordinate_displacement\":" << aurostd::utype2string<double>(mapping_info_duplicate[j].coordinate_displacement,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"failure\":" << aurostd::utype2string<double>(mapping_info_duplicate[j].failure,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      if(mapping_info_duplicate[j].is_magnetic_misfit){
        // [BETA] sscontent_json << "\"misfit_magnetic\":" << mapping_info_duplicate[j].magnetic_misfit << eendl;
        // [BETA] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        sscontent_json << "\"displacement_magnetic\":" << aurostd::utype2string<double>(mapping_info_duplicate[j].magnetic_displacement,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        sscontent_json << "\"failure_magnetic\":" << aurostd::utype2string<double>(mapping_info_duplicate[j].magnetic_failure,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      }
      sscontent_json << "\"number_compounds_matching_structure\":" << structures_duplicate[j]->number_compounds_matching_structure << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      if(structures_duplicate[j]->properties.size()>0){
        vcontent_json.push_back(printPropertiesOfStructure(structures_duplicate[j]));
      }
      // print structure transformation info
      if(print_transformation){
        vcontent_json.push_back(printStructureTransformationInformation(mapping_info_duplicate[j])); //DX20210115
      }
      vstructures.push_back("{" + aurostd::joinWDelimiter(vcontent_json,",") + "}");
      vcontent_json.clear();
    }
  }
  // print same family structure information
  else if(mode == family_structure_xf){
    for(uint j=0;j<structures_family.size();j++){
      sscontent_json << "\"name\":" << "\"" << structures_family[j]->name << "\"" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"misfit\":" << aurostd::utype2string<double>(mapping_info_family[j].misfit,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"lattice_deviation\":" << aurostd::utype2string<double>(mapping_info_family[j].lattice_deviation,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"coordinate_displacement\":" << aurostd::utype2string<double>(mapping_info_family[j].coordinate_displacement,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      sscontent_json << "\"failure\":" << aurostd::utype2string<double>(mapping_info_family[j].failure,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      if(mapping_info_family[j].is_magnetic_misfit){
        // [BETA] sscontent_json << "\"misfit_magnetic\":" << mapping_info_family[j].magnetic_misfit << eendl;
        // [BETA] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        sscontent_json << "\"displacement_magnetic\":" << aurostd::utype2string<double>(mapping_info_family[j].magnetic_displacement,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
        sscontent_json << "\"failure_magnetic\":" << aurostd::utype2string<double>(mapping_info_family[j].magnetic_failure,(double)AUROSTD_ROUNDOFF_TOL,roff) << eendl;
        vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      }
      sscontent_json << "\"number_compounds_matching_entry\":" << structures_family[j]->number_compounds_matching_structure << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      if(structures_family[j]->properties.size()>0){
        vcontent_json.push_back(printPropertiesOfStructure(structures_family[j]));
      }
      // print structure transformation info
      if(print_transformation){
        vcontent_json.push_back(printStructureTransformationInformation(mapping_info_family[j])); //DX20210115
      }
      vstructures.push_back("{" + aurostd::joinWDelimiter(vcontent_json,",") + "}");
      vcontent_json.clear();
    }
  }

  return aurostd::joinWDelimiter(vstructures,",");
}

// ***************************************************************************
// StructurePrototype::numberOfComparisons()
// ***************************************************************************
uint StructurePrototype::numberOfComparisons(){

  // Count the number of comparisons

  return structures_duplicate.size();
}

// ***************************************************************************
// StructurePrototype::isSymmetryCalculated()
// ***************************************************************************
bool StructurePrototype::isSymmetryCalculated(){

  // Check if the space group symmetry is calculated for the prototype
  // If the space group is calculated, then the Wyckoff positions are also
  // calculated (by construction)

  return (space_group>0 && space_group<231);
}

// ***************************************************************************
// XtalFinderCalculator::isSymmetryCalculated() //DX20201115
// ***************************************************************************
bool XtalFinderCalculator::isSymmetryCalculated(structure_container& structure){

  // Check if the space group symmetry is calculated for a particular
  // structure in the container
  // If the space group is calculated, then the Wyckoff positions are also
  // calculated (by construction)

  return (structure.space_group>0 && structure.space_group<231);
}

// ***************************************************************************
// StructurePrototype::isLFAEnvironmentCalculated()
// ***************************************************************************
bool StructurePrototype::isLFAEnvironmentCalculated(){

  // Check if the LFA environment is calculated for the prototype

  return (environments_LFA.size()>0);
}

// ***************************************************************************
// XtalFinderCalculator::isLFAEnvironmentCalculated()
// ***************************************************************************
bool XtalFinderCalculator::isLFAEnvironmentCalculated(
    structure_container& structure){

  // Check if the LFA environment is calculated for a particular
  // structure in the container

  return (structure.environments_LFA.size()>0);
}

// ***************************************************************************
// XtalFinderCalculator::areNearestNeighborsCalculated()
// ***************************************************************************
bool XtalFinderCalculator::areNearestNeighborsCalculated(
    structure_container& structure){

  // Check if the nearest neighbor distances have been calculated for
  // a particular structure in the container

  return (structure.nearest_neighbor_distances.size()>0);
}

// ***************************************************************************
// StructurePrototype::calculateSymmetry()
// ***************************************************************************
bool StructurePrototype::calculateSymmetry(){

  // Calculate the symmetry of the representative structure
  // (i.e., space group and Wyckoff positions)
  // The Pearson symbol calculation takes more time, so I do not calculate it

  Pearson = "";
  bool no_scan = false; //DX20191230
  // use sym_eps if given, otherwise use default
  if(!structure_representative->structure.sym_eps_calculated || structure_representative->structure.sym_eps==AUROSTD_NAN) { structure_representative->structure.sym_eps = SYM::defaultTolerance(structure_representative->structure); } //DX20210421
  double use_tol = structure_representative->structure.sym_eps; //DX20210421
  space_group = structure_representative->structure.SpaceGroup_ITC(use_tol,-1,SG_SETTING_ANRL,no_scan); //DX20191230
  vector<GroupedWyckoffPosition> tmp_grouped_Wyckoff_positions;
  compare::groupWyckoffPositions(structure_representative->structure, tmp_grouped_Wyckoff_positions);
  grouped_Wyckoff_positions = tmp_grouped_Wyckoff_positions;
  return true;
}

// ***************************************************************************
// XtalFinderCalculator::calculateSymmetry()
// ***************************************************************************
void XtalFinderCalculator::calculateSymmetry(structure_container& str_rep){

  // Calculate the symmetry of the representative structure
  // (i.e., space group and Wyckoff positions)
  // The Pearson symbol calculation takes more time, so I do not calculate it

  str_rep.Pearson = "";
  bool no_scan = false; //DX20191230
  // use sym_eps if given, otherwise use default
  if(!str_rep.structure.sym_eps_calculated || str_rep.structure.sym_eps==AUROSTD_NAN) { str_rep.structure.sym_eps = SYM::defaultTolerance(str_rep.structure); } //DX20210421
  double use_tol = str_rep.structure.sym_eps; //DX20210421
  str_rep.space_group = str_rep.structure.SpaceGroup_ITC(use_tol,-1,SG_SETTING_ANRL,no_scan); //DX20191230
  vector<GroupedWyckoffPosition> tmp_grouped_Wyckoff_positions;
  compare::groupWyckoffPositions(str_rep.structure, tmp_grouped_Wyckoff_positions);
  str_rep.grouped_Wyckoff_positions = tmp_grouped_Wyckoff_positions;
}

// ***************************************************************************
// XtalFinderCalculator::getSymmetryInfoFromXstructure() //DX20210104
// ***************************************************************************
void XtalFinderCalculator::getSymmetryInfoFromXstructure(structure_container& str_rep){

  // Get the symmetry info from xstructure

  str_rep.Pearson = "";
  str_rep.space_group = str_rep.structure.space_group_ITC;
  vector<GroupedWyckoffPosition> tmp_grouped_Wyckoff_positions;
  compare::groupWyckoffPositions(str_rep.structure, tmp_grouped_Wyckoff_positions);
  str_rep.grouped_Wyckoff_positions = tmp_grouped_Wyckoff_positions;
}

// ***************************************************************************
// XtalFinderCalculator::addStructure2container()
// ***************************************************************************
void XtalFinderCalculator::addStructure2container(
    const xstructure& xstr,
    const string& structure_name,
    const string& source,
    uint relaxation_step,
    bool same_species){

  // Add crystal structure to XtalFinderCalculator container to
  // easily pass between functions.

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  structure_container str_rep_tmp;

  // ---------------------------------------------------------------------------
  // adds information to container and pre-conditions xstructure
  str_rep_tmp = compare::initializeStructureContainer(xstr,same_species);
  str_rep_tmp.name = structure_name;
  str_rep_tmp.is_structure_generated = true;
  str_rep_tmp.source = source;
  str_rep_tmp.relaxation_step = relaxation_step;

  // ---------------------------------------------------------------------------
  // check if fake names for same species comparison
  if(same_species && !pflow::hasRealElements(str_rep_tmp.structure)){
    message << "Atomic species are not real/physical " << str_rep_tmp.name << " cannot perform material comparison; skipping structure.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
    return; // not storing structure
  }

  if(LDEBUG){ cerr << __AFLOW_FUNC__ << " loaded " << structure_name << " structure." << endl; }

  structure_containers.push_back(str_rep_tmp);

}

// ***************************************************************************
// XtalFinderCalculator::removeStructureFromContainerByName()
// ***************************************************************************
void XtalFinderCalculator::removeStructureFromContainerByName(
    const string& structure_name){

  // Remove crystal structure from XtalFinderCalculator container.
  // CAREFUL: This will rewrite the structure container vector and
  // any pointers to this vector will be invalid.
  // Removing structures by name instead of pointer for this same reason
  // (i.e., once we remove by one element, the address/index are no longer
  // the same).

  //bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  //stringstream message;

  for(uint i=0;i<structure_containers.size();i++){
    if(structure_containers[i].name == structure_name){
      structure_containers.erase(structure_containers.begin()+i);
      break;
    }
  }

}

// ***************************************************************************
// XtalFinderCalculator::setStructureAsRepresentative() //DX20201204
// ***************************************************************************
void XtalFinderCalculator::setStructureAsRepresentative(
    StructurePrototype& structure_tmp,
    uint container_index){

  setStructureAsRepresentative(structure_tmp, &structure_containers[container_index]);
}

void XtalFinderCalculator::setStructureAsRepresentative(StructurePrototype& structure_tmp,
    structure_container* str_pointer){

  // Point structure representative to a particular structure in the
  // XtalFinderCalculator container.

  //bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  // check pointer is not null
  if(str_pointer == NULL){
    stringstream message;
    message << "Input is a null pointer.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }

  // point to the input
  structure_tmp.structure_representative = str_pointer;

  // update atom/type counts and stoichiometry
  structure_tmp.stoichiometry = structure_tmp.structure_representative->stoichiometry;
  structure_tmp.natoms = structure_tmp.structure_representative->structure.atoms.size();
  structure_tmp.ntypes = structure_tmp.structure_representative->structure.num_each_type.size();
  structure_tmp.elements = structure_tmp.structure_representative->elements;

  // update symmetry and environment info
  structure_tmp.Pearson = structure_tmp.structure_representative->Pearson;
  structure_tmp.space_group = structure_tmp.structure_representative->space_group;
  structure_tmp.grouped_Wyckoff_positions = structure_tmp.structure_representative->grouped_Wyckoff_positions;
  structure_tmp.environments_LFA = structure_tmp.structure_representative->environments_LFA;

  // properties
  structure_tmp.property_names = structure_tmp.structure_representative->properties_names;
  structure_tmp.property_units = structure_tmp.structure_representative->properties_units;
  structure_tmp.property_types = structure_tmp.structure_representative->properties_types;

}

// ***************************************************************************
// XtalFinderCalculator::addStructure2duplicatesList() //DX20201204
// ***************************************************************************
// Add structure to the duplicate list for the StructurePrototype object
// index input
void XtalFinderCalculator::addStructure2duplicatesList(StructurePrototype& structure_tmp,
    uint container_index){

  addStructure2duplicatesList(structure_tmp,&structure_containers[container_index]);
}

// structure_container pointer input
void XtalFinderCalculator::addStructure2duplicatesList(StructurePrototype& structure_tmp,
    structure_container* str_pointer){

  // check pointer is not null
  if(str_pointer == NULL){
    stringstream message;
    message << "Input is a null pointer.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }

  structure_tmp.structures_duplicate.push_back(str_pointer);

  // initialize
  structure_mapping_info temp_misfit_info = compare::initialize_misfit_struct();
  structure_tmp.mapping_info_duplicate.push_back(temp_misfit_info);

  // update properties
  if(structure_tmp.structures_duplicate.back()->properties.size()>0){
    for(uint i=0;i<structure_tmp.structures_duplicate.back()->properties_names.size();i++){
      if(!aurostd::WithinList(structure_tmp.property_names, structure_tmp.structures_duplicate.back()->properties_names[i])){
        structure_tmp.property_names.push_back(structure_tmp.structures_duplicate.back()->properties_names[i]);
        structure_tmp.property_units.push_back(structure_tmp.structures_duplicate.back()->properties_units[i]);
      }
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::addStructure2sameFamilyList() //DX20201204
// ***************************************************************************
// Add structure to the same family list for the StructurePrototype object
void XtalFinderCalculator::addStructure2sameFamilyList(StructurePrototype& structure_tmp,
    uint index){

  // Index corresponds to element in StructurePrototype object.

  // check pointer is not null
  if(structure_tmp.structures_duplicate[index] == NULL){
    stringstream message;
    message << "Pointer at index " << index << " is a null pointer.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }

  structure_tmp.structures_family.push_back(structure_tmp.structures_duplicate[index]);

  // add misfit
  structure_tmp.mapping_info_family.push_back(structure_tmp.mapping_info_duplicate[index]);

}

void XtalFinderCalculator::addStructure2sameFamilyList(StructurePrototype& structure_tmp,
    structure_container* str_pointer){

  // Similar to variant above, but this resets the misfit object
  // (so it cannot be combined with function above).

  // check pointer is not null
  if(str_pointer == NULL){
    stringstream message;
    message << "Input is a null pointer.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }

  structure_tmp.structures_family.push_back(str_pointer);

  // initialize
  structure_mapping_info temp_misfit_info = compare::initialize_misfit_struct();
  structure_tmp.mapping_info_family.push_back(temp_misfit_info);

}

// ***************************************************************************
// StructurePrototype::copyPrototypeInformation()
// ***************************************************************************
void StructurePrototype::copyPrototypeInformation(const StructurePrototype& b){

  // Copy prototype information to new StructurePrototype object
  // (i.e., stoichiometry, number of types, Pearson symbol, space group,
  // and Wyckoff positions). When a structure does not match with its
  // original comparison group, it is moved to its own StructurePrototype
  // object and we need to update the symmetry and environment info.
  // This is different than the copy constructor, which would copy over all
  // attributes (including duplicate structures, mapping info, etc.).

  //elements=comparison_schemes[i].elements;
  stoichiometry=b.stoichiometry;
  ntypes=b.stoichiometry.size();
  Pearson=b.Pearson;
  space_group=b.space_group;
  grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
  environments_LFA=b.environments_LFA; //DX20190711
}

// ***************************************************************************
// StructurePrototype::copyDuplicate()
// ***************************************************************************
void StructurePrototype::copyDuplicate(const StructurePrototype& b,
    uint index,
    bool copy_misfit){ //DX20190730 - added copy_misfit option

  // Copy structure at index as a potential duplicate in this object

  // check pointer is not null
  if(b.structures_duplicate[index] == NULL){
    stringstream message;
    message << "Pointer at index " << index << " is a null pointer.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
  }

  // point to the structure in the container at correponding index
  structures_duplicate.push_back(b.structures_duplicate[index]);

  // signifies comparison has not been performed yet
  if(!copy_misfit){
    structure_mapping_info temp_misfit_info = compare::initialize_misfit_struct(); //DX20191217
    mapping_info_duplicate.push_back(temp_misfit_info); //DX20191217
  }
  else{ //DX20190730
    mapping_info_duplicate.push_back(b.mapping_info_duplicate[index]); //DX20191217
  }

}

// ***************************************************************************
// StructurePrototype::removeNonDuplicate()
// ***************************************************************************
void StructurePrototype::removeNonDuplicate(uint index){

  // If structure (structures_duplicate[index]) does not match with the
  // original representative structure (i.e., misfit>misfit_match),
  // remove from scheme

  structures_duplicate.erase(structures_duplicate.begin()+index);
  mapping_info_duplicate.erase(mapping_info_duplicate.begin()+index); //DX20191217

}
// ***************************************************************************
// END:: StructurePrototype Class
// ***************************************************************************

// ***************************************************************************
// XtalFinderCalculator::XtalFinderCalculator() (Constructors)
// ***************************************************************************
XtalFinderCalculator::XtalFinderCalculator(uint num_proc_input, ostream& oss): xStream(oss) {
  free();
  num_proc=num_proc_input;
}

XtalFinderCalculator::XtalFinderCalculator(ofstream& FileMESSAGE, uint num_proc_input, ostream& oss): xStream(FileMESSAGE,oss) {
  free();
  num_proc=num_proc_input;
}

XtalFinderCalculator::XtalFinderCalculator(double misfit_match_input, double misfit_family_input, uint num_proc_input, ostream& oss) : xStream(oss) {
  free();
  num_proc=num_proc_input;
  misfit_match=misfit_match_input;
  misfit_family=misfit_family_input;
}

XtalFinderCalculator::XtalFinderCalculator(double misfit_match_input, double misfit_family_input, ofstream& FileMESSAGE, uint num_proc_input, ostream& oss) : xStream(FileMESSAGE,oss) {
  free();
  num_proc=num_proc_input;
  misfit_match=misfit_match_input;
  misfit_family=misfit_family_input;
}

// ***************************************************************************
// XtalFinderCalculator::free()
// ***************************************************************************
void XtalFinderCalculator::free(){
  misfit_match = DEFAULT_XTALFINDER_MISFIT_MATCH;
  misfit_family = DEFAULT_XTALFINDER_MISFIT_FAMILY;
  num_proc=1;
  structure_containers.clear();
}

// ***************************************************************************
// XtalFinderCalculator::~XtalFinderCalculator() (Destructor)
// ***************************************************************************
XtalFinderCalculator::~XtalFinderCalculator(){
  xStream::free();
  free();
}

// ***************************************************************************
// XtalFinderCalculator::clear() (clear)
// ***************************************************************************
void XtalFinderCalculator::clear(){
  free();
}

// ***************************************************************************
// XtalFinderCalculator::XtalFinderCalculator(...) (copy)
// ***************************************************************************
XtalFinderCalculator::XtalFinderCalculator(const XtalFinderCalculator& b): xStream(*b.getOFStream(),*b.getOSS()) {
  if(this!=&b){
    copy(b);
  }
}

// ***************************************************************************
// XtalFinderCalculator::copy()
// ***************************************************************************
void XtalFinderCalculator::copy(const XtalFinderCalculator& b) {
  if(this==&b){ return; }
  xStream::copy(b);
  misfit_match=b.misfit_match;
  misfit_family=b.misfit_family;
  num_proc=b.num_proc;
  structure_containers=b.structure_containers;
}

// ***************************************************************************
// XtalFinderCalculator::operator= (assignment)
// ***************************************************************************
const XtalFinderCalculator& XtalFinderCalculator::operator=(const XtalFinderCalculator& b){
  if(this!=&b){
    copy(b);
  }
  return *this;
}

// ***************************************************************************
// XtalFinderCalculator::operator<< (output/print)
// ***************************************************************************
ostream& operator<<(ostream& oss, const XtalFinderCalculator& xtalfinder){

  oss << "XtalFinderCalculator:" << endl;
  oss << " misfit match: " << xtalfinder.misfit_match << endl;
  oss << " misfit family: " << xtalfinder.misfit_family << endl;
  oss << " number of processors: " << xtalfinder.num_proc << endl;
  oss << " structures in container (" << xtalfinder.structure_containers.size() << "):" << endl;
  for(uint i=0;i<xtalfinder.structure_containers.size();i++){
    oss << "  " << xtalfinder.structure_containers[i].name
      << " (compound=" << xtalfinder.structure_containers[i].compound
      << ", SG=" << xtalfinder.structure_containers[i].space_group << ")" << endl;
  }
  return oss;
}

// ***************************************************************************
// START:: GroupedWyckoffPosition Class
// ***************************************************************************

// ***************************************************************************
// GroupedWyckoffPosition::GroupedWyckoffPosition() (constructor)
// ***************************************************************************
GroupedWyckoffPosition::GroupedWyckoffPosition(){
  free();
}

// ***************************************************************************
// GroupedWyckoffPosition::free()
// ***************************************************************************
void GroupedWyckoffPosition::free(){
  type = 0;
  element = "";
  site_symmetries.clear();
  multiplicities.clear();
  letters.clear();
}

// ***************************************************************************
// GroupedWyckoffPosition::~GroupedWyckoffPosition() (Destructor)
// ***************************************************************************
GroupedWyckoffPosition::~GroupedWyckoffPosition(){
  free();
}

// ***************************************************************************
// GroupedWyckoffPosition::GroupedWyckoffPosition(...) (copy)
// ***************************************************************************
GroupedWyckoffPosition::GroupedWyckoffPosition(const GroupedWyckoffPosition& b){
  if(this!=&b){
    copy(b);
  }
}

// ***************************************************************************
// GroupedWyckoffPosition::copy()
// ***************************************************************************
void GroupedWyckoffPosition::copy(const GroupedWyckoffPosition& b) {
  if(this != &b){
    type=b.type;
    element=b.element;
    site_symmetries=b.site_symmetries;
    multiplicities=b.multiplicities;
    letters=b.letters;
  }
}

// ***************************************************************************
// GroupedWyckoffPosition::operator= (assignment)
// ***************************************************************************
const GroupedWyckoffPosition& GroupedWyckoffPosition::operator=(const GroupedWyckoffPosition& b){
  if(this!=&b){
    free();
    copy(b);
  }
  return *this;
}

// ***************************************************************************
// GroupedWyckoffPosition::operator<< (output/print)
// ***************************************************************************
ostream& operator<<(ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition){

  // operator<< for the GroupedWyckoffPosition object (looks like a JSON)

  //if(StructurePrototype.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  StructurePrototype.iomode=JSON_MODE;
  //}

  string eendl="";
  stringstream sscontent_json;
  vector<string> vcontent_json;

  // type
  sscontent_json << "\"type\":" << GroupedWyckoffPosition.type << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // element
  sscontent_json << "\"element\":\"" << GroupedWyckoffPosition.element << "\"" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // site_symmetries
  sscontent_json << "\"site_symmetries\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(GroupedWyckoffPosition.site_symmetries,"\""),",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // multiplicities
  sscontent_json << "\"multiplicities\":[" << aurostd::joinWDelimiter(GroupedWyckoffPosition.multiplicities,",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // letters
  sscontent_json << "\"letters\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(GroupedWyckoffPosition.letters,"\""),",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // Put into json GroupedWyckoffPosition object
  oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
  vcontent_json.clear();

  return oss;
}

// ***************************************************************************
// GroupedWyckoffPosition::operator< (less than operator, for sorting)
// ***************************************************************************
bool GroupedWyckoffPosition::operator<(const GroupedWyckoffPosition& b) const{

  // allows sorting vector by element name

  if(element<b.element){
    return true;
  }
  return false;
}

// ***************************************************************************
// END:: GroupedWyckoffPosition Class
// ***************************************************************************

// ***************************************************************************
// AtomEnvironment Class
// ***************************************************************************
// MOVED CLASS TO aflow_xatom.cpp (DX20191120)

// ***************************************************************************
// initialize_misfit_struct() //DX20201218
// ***************************************************************************
namespace compare{
  structure_mapping_info initialize_misfit_struct(bool magnetic){

    // Initialize misfit structure object
    // DX20200317 - set attributes explicitly
    // { .attribute=<>, .attribute=<>, ...} doesn't work for old GCC versions

    structure_mapping_info misfit_info;
    resetMisfitInfo(misfit_info, magnetic); //DX20220406 - consolidated
    return misfit_info;
  }
}

// ***************************************************************************
// compare::resetMisfitInfo() //DX20220406
// ***************************************************************************
namespace compare{
  void resetMisfitInfo(structure_mapping_info& misfit_info, bool magnetic){

    // Reset misfit_info struct to default values
    // DX20200317 - set attributes explicitly
    // { .attribute=<>, .attribute=<>, ...} doesn't work for old GCC versions

    misfit_info.is_magnetic_misfit=magnetic;
    misfit_info.misfit=AUROSTD_MAX_DOUBLE;
    misfit_info.lattice_deviation=AUROSTD_MAX_DOUBLE;
    misfit_info.coordinate_displacement=AUROSTD_MAX_DOUBLE;
    misfit_info.failure=AUROSTD_MAX_DOUBLE;
    misfit_info.magnetic_misfit=AUROSTD_MAX_DOUBLE;
    misfit_info.magnetic_displacement=AUROSTD_MAX_DOUBLE;
    misfit_info.magnetic_failure=AUROSTD_MAX_DOUBLE;

    uint resize=0;
    resizeMappingInfo(misfit_info,resize);
  }
}

// ***************************************************************************
// compare::resetMappingInfo() //DX20220406
// ***************************************************************************
namespace compare{
  void resetMappingInfo(structure_mapping_info& misfit_info){

    // Reset anything related to the mapping info
    // including the misfit and coordinate devation
    // NOT the lattice displacement

    misfit_info.misfit=AUROSTD_MAX_DOUBLE;
    misfit_info.coordinate_displacement=AUROSTD_MAX_DOUBLE;
    misfit_info.failure=AUROSTD_MAX_DOUBLE;
    misfit_info.magnetic_misfit=AUROSTD_MAX_DOUBLE;
    misfit_info.magnetic_displacement=AUROSTD_MAX_DOUBLE;
    misfit_info.magnetic_failure=AUROSTD_MAX_DOUBLE;

    uint resize=0;
    resizeMappingInfo(misfit_info,resize);
  }
}

// ***************************************************************************
// compare::resizeMappingInfo() //DX20201218
// ***************************************************************************
namespace compare{
  void resizeMappingInfo(structure_mapping_info& str_mis, uint size){

    // Set the size of the atom mapping information

    str_mis.atom_map.clear();
    str_mis.basis_map.clear();
    str_mis.distances_mapped.clear();
    str_mis.vectors_mapped.clear();

    str_mis.atom_map.resize(size);
    str_mis.basis_map.resize(size);
    str_mis.distances_mapped.resize(size);
    str_mis.vectors_mapped.resize(size);
  }
}


// ***************************************************************************
// *                                                                         *
// *                             FUNCTIONS                                   *
// *                                                                         *
// ***************************************************************************

// ***************************************************************************
// loadDefaultComparisonOptions()
// ***************************************************************************
namespace compare {
  aurostd::xoption loadDefaultComparisonOptions(const string& mode){ //DX20200103

    aurostd::xoption comparison_options;

    if(mode=="permutation"){
      comparison_options.flag("COMPARISON_OPTIONS::MISFIT_MATCH",DEFAULT_XTALFINDER_MISFIT_MATCH); //DX20201119
      comparison_options.flag("COMPARISON_OPTIONS::MISFIT_FAMILY",DEFAULT_XTALFINDER_MISFIT_FAMILY); //DX20201119
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",TRUE); // permutations are generated from the same structure, so they will have the same volume anyway
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",FALSE); // duplicate permutations should have the same symmetry
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",FALSE); // duplicate permutations should have the same symmetry
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",FALSE); // duplicate permutations should have the same environment
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES",FALSE); //DX20200320
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",TRUE); // remove unmatched structures from object
      comparison_options.flag("COMPARISON_OPTIONS::PRINT_MATCHES_TO_INPUT_ONLY",FALSE); // only print matches to the input structure (database comparisons) //DX20201230
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::CHECK_OTHER_GROUPINGS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::STORE_COMPARISON_LOGS",FALSE); // do not store comparison logs
      comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",FALSE); // compare all permutations until matched or exhausted all comparisons
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",FALSE); //DX20201006
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",FALSE); //DX20201006
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",FALSE); //DX20201006
    }
    else{
      comparison_options.flag("COMPARISON_OPTIONS::MISFIT_MATCH",DEFAULT_XTALFINDER_MISFIT_MATCH); //DX20201119
      comparison_options.flag("COMPARISON_OPTIONS::MISFIT_FAMILY",DEFAULT_XTALFINDER_MISFIT_FAMILY); //DX20201119
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES",FALSE); //DX20200320
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::PRINT_MATCHES_TO_INPUT_ONLY",FALSE); // only print matches to the input structure (database comparisons) //DX20201230
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::CHECK_OTHER_GROUPINGS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::STORE_COMPARISON_LOGS",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",FALSE); //DX20201006
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",FALSE); //DX20201006
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",FALSE); //DX20201006
    }

    return comparison_options;
  }
}

// ***************************************************************************
// loadStructuresFromDirectory()
// ***************************************************************************
void XtalFinderCalculator::loadStructuresFromDirectory(
    const string& directory,
    const vector<string>& magmoms_for_systems,
    bool same_species) {

  // Load all structures from a directory into
  // XtalFinderCalculator.structure_containers

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  string source = "file"; // all structures from file in the directory
  uint relaxation_step = 0; // assume unrelaxed; impossible to know from this input-type

  string structure_name = "";
  vector<string> vfiles;
  aurostd::DirectoryLS(directory, vfiles);
  std::sort(vfiles.begin(),vfiles.end()); //CO20180830

  message << "Loading " << vfiles.size() << " files in directory ... ";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  for(uint i=0; i<vfiles.size(); i++){
    if(LDEBUG) {cerr << "compare:: " << i << "/" << vfiles.size() << " " << vfiles[i] << endl;}
    if(vfiles[i].find("comparison_output.json") != std::string::npos ||
        vfiles[i].find("comparison_output.out") != std::string::npos ||
        vfiles[i].find("comparison_no_duplicate_compounds_output.json") != std::string::npos ||
        vfiles[i].find("comparison_no_duplicate_compounds_output.out") != std::string::npos ||
        vfiles[i].find("duplicate_compounds_output.json") != std::string::npos ||
        vfiles[i].find("duplicate_compounds_output.out") != std::string::npos ||
        vfiles[i].find("nohup.out") != std::string::npos){
      message << "Ignoring file=" << vfiles[i] << endl;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      vfiles.erase(vfiles.begin()+i);
      i--;
    }
    else {
      structure_name = directory+"/"+vfiles[i];
      xstructure xstr;
      if(!compare::generateStructure(structure_name, source, 0, xstr, *p_oss)){ continue; }

      // add magnetic
      if(magmoms_for_systems.size()==vfiles.size()){
        try { pflow::ProcessAndAddSpinToXstructure(xstr, magmoms_for_systems[i]); } //DX20190801
        catch(aurostd::xerror& excpt) { message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "."; throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); } //DX20190801
      }

      // adds the structure to a container in XtalFinderCalculator to be passed easily by reference
      addStructure2container(xstr, structure_name, source, relaxation_step, same_species);
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::loadStructuresFromAflowlibEntries() //DX20201222
// ***************************************************************************
void XtalFinderCalculator::loadStructuresFromAflowlibEntries(
    const vector<aflowlib::_aflowlib_entry>& entries,
    const vector<string>& magmoms_for_systems,
    bool same_species){

  // Load all structures from aflowlib entries into a vector of
  // XtalFinderCalculator.structure_containers

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  string structure_name = "", source = "aflowlib";
  uint relaxation_step = 2; // input assumed to be most relaxed

  // ---------------------------------------------------------------------------
  // loop through aflowlib entries
  for(uint i=0;i<entries.size();i++){
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Loading entry " << i << ": " << entries[i].auid << endl;}

    // load structure (if possible with try-catch)
    xstructure xstr_tmp;
    try { xstr_tmp = entries[i].vstr.back(); } //back() is the most relaxed structure
    catch(aurostd::xerror& excpt) {
      message << "Could not load entry " << i << ": " << entries[i].auid << "...skipping entry";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      continue;
    }
    if(magmoms_for_systems.size()==entries.size()){
      try { pflow::ProcessAndAddSpinToXstructure(xstr_tmp, magmoms_for_systems[i]); }
      catch(aurostd::xerror& excpt) {
        message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "...skipping structure";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }
    }
    structure_name = entries[i].auid;

    // adds the structure to a container in XtalFinderCalculator to be passed easily by reference
    addStructure2container(xstr_tmp, structure_name, source, relaxation_step, same_species);

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Successfully loaded entry " << i << ": " << entries[i].auid << endl;}
  }
}

// ***************************************************************************
// XtalFinderCalculator::loadStructuresFromFile()
// ***************************************************************************
void XtalFinderCalculator::loadStructuresFromFile(
    const string& filename,
    const vector<string>& magmoms_for_systems,
    bool same_species){

  // Load all structures from a file into
  // XtalFinderCalculator.structure_containers.
  // Useful for reading in aflow.in relaxation steps or pocc structures

  // ---------------------------------------------------------------------------
  // file to stringstream
  stringstream input_file;
  aurostd::efile2stringstream(filename, input_file);
  loadStructuresFromStringstream(input_file, magmoms_for_systems, same_species);
}

// ME20210206 - Added stringstream variant
void XtalFinderCalculator::loadStructuresFromStringstream(
    stringstream& input_stream,
    const vector<string>& magmoms_for_systems,
    bool same_species){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  // ---------------------------------------------------------------------------
  // tokenize stringstream by newline
  vector<string> lines;
  aurostd::string2tokens(input_stream.str(),lines,"\n");

  // ---------------------------------------------------------------------------
  // structure delimiters
  string START=_VASP_POSCAR_MODE_EXPLICIT_START_;
  string STOP=_VASP_POSCAR_MODE_EXPLICIT_STOP_;

  // ---------------------------------------------------------------------------
  // used to find the total number of structures
  vector<string> start_string;
  aurostd::substring2strings(input_stream.str(),start_string,string(START),string(STOP)); //CO20230502

  message << "Loading " << start_string.size() << " structures in file ... ";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  // parse all lines
  vector<string> structure_lines;
  bool structure_lines_found = false;
  uint structure_count = 0;
  stringstream geometry;geometry.clear();geometry.str(std::string());
  for(uint i=0;i<lines.size();i++){
    if(lines[i].find(START) != std::string::npos){
      geometry.clear();geometry.str("");
      structure_lines_found = true;
      structure_count+=1;
    }
    else if(structure_lines_found && lines[i].find(STOP) == std::string::npos){
      geometry << lines[i] << endl;
    }
    else if(lines[i].find(STOP) != std::string::npos){
      structure_lines_found = false;
      structure_lines.push_back(geometry.str());
    }
  }
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " structure_count=" << structure_count << endl;}

  uint relaxation_step = 0;
  for(uint i=0;i<structure_lines.size();i++){
    if(LDEBUG) {cerr << "compare:: loading " << i << "/" << structure_lines.size() << endl;}

    stringstream designation; designation << "file structure # " << i+1 << "/" << start_string.size();  //CO20211118 - otherwise last structure is 81/82
    xstructure xstr;
    compare::generateStructure("input geometry", structure_lines[i], relaxation_step, xstr, *p_oss);

    // add magnetic
    if(magmoms_for_systems.size()==structure_lines.size()){
      try { pflow::ProcessAndAddSpinToXstructure(xstr, magmoms_for_systems[i]); } //DX20190801
      catch(aurostd::xerror& excpt) { message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "."; throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); } //DX20190801
    }

    // adds the structure to a container in XtalFinderCalculator to be passed easily by reference
    addStructure2container(xstr, designation.str(), structure_lines[i], relaxation_step, same_species);
  }
}

//DX20190424 START
// ***************************************************************************
// XtalFinderCalculator::loadStructuresFromStructureList()
// ***************************************************************************
void XtalFinderCalculator::loadStructuresFromStructureList(
    const vector<string>& filenames,
    const vector<string>& magmoms_for_systems,
    bool same_species){

  // Load all structures from a vector of filenames into
  // XtalFinderCalculator.structure_containers

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  string source = "file";
  uint relaxation_step = 0; // assume unrelaxed; impossible to know from this input-type

  // ---------------------------------------------------------------------------
  // read in files
  for(uint i=0;i<filenames.size();i++){
    if(!aurostd::FileExist(filenames[i])){
      message << filenames[i] << " file not found.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_FILE_NOT_FOUND_);
    }

    // generate structure
    xstructure xstr;
    compare::generateStructure(filenames[i], source, 0, xstr, *p_oss);

    // add magnetic
    if(magmoms_for_systems.size()==filenames.size()){
      try { pflow::ProcessAndAddSpinToXstructure(xstr, magmoms_for_systems[i]); } //DX20190801
      catch(aurostd::xerror& excpt) { message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "."; throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_); } //DX20190801
    }

    // adds the structure to a container in XtalFinderCalculator to be passed easily by reference
    addStructure2container(xstr, filenames[i], source, relaxation_step, same_species);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << ": loaded structure " << i << endl;
    }
  }
}

//DX20191105 - load multiple structures (useful for multithreaded structure loading) - START
// ***************************************************************************
// compare::generateStructures()
// ***************************************************************************
namespace compare {
  void generateStructures(vector<StructurePrototype>& structures, ostream& oss, uint start_index, uint end_index){

    // if end index is greater than structures.size(), then generate over entire range
    if(end_index > structures.size()){ end_index=structures.size(); }

    for(uint i=start_index;i<end_index;i++){
      if(!structures[i].structure_representative->is_structure_generated){
        structures[i].structure_representative->is_structure_generated = generateStructure(
            structures[i].structure_representative->name,
            structures[i].structure_representative->source,
            structures[i].structure_representative->relaxation_step, //DX20200429
            structures[i].structure_representative->structure,
            oss);
      }
    }
  }
}
//DX20191105 - load multiple structures (useful for multithreaded structure loading) - END

// ***************************************************************************
// generateStructure
// ***************************************************************************
namespace compare {
  bool generateStructure(const string& structure_name,
      const string& structure_source,
      uint relaxation_step,
      xstructure& structure,
      ostream& oss){ //DX20200429 - added relaxation_step

    // generate the xstructure object, having this separate function allows us to load structures in
    // a threaded environment

    // there are multiple modes of generating xstructures
    // 1) AFLOW prototypes
    // 2) AURL
    // 3) input (cin)
    // 4) from files (directory, list of files, single file)

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
    ofstream FileMESSAGE;
    vector<string> tokens;

    if(relaxation_step){} //CO20200508 - keep it busy

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " generating structure: " << structure_name << " from " << structure_source << endl;
    }

    // ---------------------------------------------------------------------------
    // load from AFLOW prototypes
    if(structure_source=="aflow_prototypes"){
      // htqc or anrl
      structure = aflowlib::PrototypeLibraries(oss,structure_name,"",2); //DX20200103 - cout to oss
    }
    // ---------------------------------------------------------------------------
    // load from AURL
    else if(structure_source.find("aurl") != std::string::npos){
      //DX20200225 - check if relaxation_step is appended
      uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_; //default
      bool load_most_relaxed_structure_only = true;
      if(relaxation_step != _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){ load_most_relaxed_structure_only = false; }
      aflowlib::_aflowlib_entry entry; entry.aurl = structure_name;
      vector<string> structure_files;
      //DX20190326 - need to put url path, i.e., structure name, [OBSOLETE] if(!pflow::loadXstructures(entry,FileMESSAGE,oss,true,structure_name,true)){ cerr << __AFLOW_FUNC__ << "WARNING::Could not load structure via aurl..." << endl; return false;}
      //DX ORIG B4 20191105 - if(!pflow::loadXstructures(entry,FileMESSAGE,oss,true,structure_name,true)){ cerr << __AFLOW_FUNC__ << "WARNING::Could not load structure via aurl..." << endl; return false;} //DX20190326
      if(!pflow::loadXstructures(entry,structure_files,FileMESSAGE,oss,load_most_relaxed_structure_only)){
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Could not load structure (aurl="+entry.aurl+") ... skipping...", FileMESSAGE, oss, _LOGGER_WARNING_);
        return false;
      }
      //DX20200225 - added compare to particular geometry files - START
      bool found_structure = false;
      uint structure_index = 0;
      if(load_most_relaxed_structure_only && entry.vstr.size()==1){
        found_structure = true;
        structure_index = 0;
      }
      else if(!load_most_relaxed_structure_only){
        if(entry.vstr.size()==3 && structure_files.size()==3){
          if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ &&
              (structure_files[0] == "POSCAR.orig" ||
               structure_files[0] == "POSCAR.relax1")){
            structure_index = 0;
            found_structure = true;
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded original structure: " << structure_files[0] << endl;}
          }
          else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ &&
              (structure_files[1] == "POSCAR.relax2" ||
               structure_files[1] == "CONTCAR.relax1")){
            structure_index = 1;
            found_structure = true;
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded relax1 structure: " << structure_files[1] << endl;}
          }
        }
      }
      //DX20200225 - added compare to particular geometry files - END
      if(found_structure){
        structure = entry.vstr[structure_index];
      }
      else {
        cerr << __AFLOW_FUNC__ << "::WARNING: More structures loaded than anticipated." << endl;
        return false;
      }
    }
    // ---------------------------------------------------------------------------
    // load from file
    else if(structure_source=="file"){
      stringstream sss;
      aurostd::efile2stringstream(structure_name,sss);
      try{
        //DX20210129 [OBSOLETE] xstructure xstr(sss);
        //DX20210129 [OBSOLETE] structure = xstr;
        structure.clear(); //DX20210129
        structure.initialize(sss); //DX20210129
      }
      catch(aurostd::xerror& excpt) { cerr << "Could not load structure " << structure_name << "...skipping structure"; return false; } //DX20190718
    }
    // ---------------------------------------------------------------------------
    // load from file containing list of structures (like aflow.in or pocc)
    else if(structure_name.find("file structure ") != std::string::npos){
      aurostd::string2tokens(structure_name,tokens,"#");
      string structure_designation = aurostd::RemoveWhiteSpaces(tokens[1]);
      uint structure_number=0; uint number_of_structures=0;
      tokens.clear(); aurostd::string2tokens(structure_designation,tokens,"/");
      structure_number=aurostd::string2utype<uint>(tokens[0]);
      number_of_structures=aurostd::string2utype<uint>(tokens[1]);
      if(number_of_structures){} //CO20200508 - keep it busy

      // ---------------------------------------------------------------------------
      // tokenize stringstream by newline
      vector<string> lines;
      aurostd::string2tokens(structure_source,lines,"\n");

      // ---------------------------------------------------------------------------
      // structure delimiters
      bool structure_lines = false;
      uint structure_count = 0;
      stringstream geometry;geometry.clear();geometry.str(std::string());
      for(uint i=0;i<lines.size();i++){
        if(lines[i].find(_VASP_POSCAR_MODE_EXPLICIT_START_) != std::string::npos){
          stringstream geometry;geometry.clear();geometry.str(std::string());
          structure_lines = true;
          structure_count+=1;
        }
        else if(structure_lines && structure_count==structure_number && lines[i].find(_VASP_POSCAR_MODE_EXPLICIT_STOP_) == std::string::npos){
          geometry << lines[i] << endl;
        }
        else if(lines[i].find(_VASP_POSCAR_MODE_EXPLICIT_STOP_) != std::string::npos && structure_lines){
          structure_lines = false;
          xstructure xstr(geometry);
          structure = xstr;
          break;
        }
      }
    }
    // ---------------------------------------------------------------------------
    // load from input (istream, e.g., from 'cat' or redirect '<')
    else if(structure_name=="input geometry"){
      stringstream sss; sss << structure_source;
      //DX20210129 [OBSOLETE] xstructure xstr(sss);
      //DX20210129 [OBSOLETE] structure.clear(); //DX20210127
      //DX20210129 [OBSOLETE] structure = xstr;
      structure.clear(); //DX20210129
      structure.initialize(sss); //DX20210129
    }
    // ---------------------------------------------------------------------------
    // load permutation
    else if(structure_source.find("permutation of: ") != std::string::npos){
      //cerr << "permutation generator: " << structure_source << endl;
      //cerr << "permutation of: " << structure_name << endl;
      string tmp_source = structure_source;
      stringstream sss; sss << aurostd::StringSubst(tmp_source, "permutation of: ", "");
      //DX20210129 [OBSOLETE] xstructure xstr(sss);
      structure.clear(); //DX20210129
      structure.initialize(sss); //DX20210129
      deque<string> species;
      for(uint j=0;j<structure_name.size();j++){stringstream ss_site; ss_site << structure_name[j]; species.push_back(ss_site.str());}
      structure.SetSpecies(species);
      structure.SpeciesPutAlphabetic();
      //DX20210202 [OBSOLETE] structure.SetNumEachType();
      //DX20210129 [OBSOLETE] structure = xstr;
    }
    // ---------------------------------------------------------------------------
    // load-type not accounted for
    else {
      stringstream message;
      message << "Structure location (from=" << structure_source << ") is not specified correctly for " << structure_name << " (i.e., input, aflow_prototype, aurl, etc.).";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, _LOGGER_WARNING_); //DX20200103 - cerr to logger
      return false;
    }

    // ---------------------------------------------------------------------------
    // pre-condition structures
    structure.ReScale(1.0); //DX20200707
    structure.BringInCell(); //DX20200707

    return true;
  }
}

// ***************************************************************************
// Find ICSD name - Find ICSD name
// ***************************************************************************
namespace compare{
  string findICSDName(const string& name){

    // Find ICSD substring within path name
    // In order for this to work, the following ICSD name format must be
    // present in the string (e.g., Ag1_ICSD_#####)

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    string ICSD_substring = "";
    bool ICSD_substring_found = false;
    if(name.find("/") != std::string::npos){ //DX20200709 - aurostd::substring2bool to find
      vector<string> tokens;
      aurostd::string2tokens(name,tokens,"/");
      for(uint i=0;i<tokens.size();i++){
        if(tokens[i].find("_ICSD_") != std::string::npos){ //DX20200709 - aurostd::substring2bool to find
          ICSD_substring = tokens[i];
          ICSD_substring_found = true;
        }
      }
    }
    else {
      if(name.find("_ICSD_") != std::string::npos){ //DX20200709 - aurostd::substring2bool to find
        ICSD_substring = name;
        ICSD_substring_found = true;
      }
    }
    if(!ICSD_substring_found){
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " WARNING: Could not find ICSD substring in name. Representative prototype will not necessarily be the minimum ICSD number." << endl;
        cerr << "string: " << name << endl;
      }
    }
    return ICSD_substring;
  }
}

// ***************************************************************************
// Find Minimum ICSD Entry - Find the minimum ICSD number in a set of names
// ***************************************************************************
namespace compare{
  string findMinimumICSDEntry(const vector<string>& ICSD_entries){

    // Identify the structure with the minimum ICSD number (to use
    // as the representative structure in the comparison)

    string min_ICSD = "";
    int min_num = 1e9;
    for(uint i=0;i<ICSD_entries.size();i++){
      if(ICSD_entries[i].empty()){ continue; } //DX20191108 - if not an ICSD, skip
      vector<string> tokens;
      aurostd::string2tokens(ICSD_entries[i],tokens,"_");
      //DX20200706 [find ICSD number, more robust] - START
      string num_string = "";
      for(uint j=0;j<tokens.size();j++){
        if(tokens[j].find("ICSD") != std::string::npos){ //DX20200709 - aurostd::substring2bool to find
          if(j+1<tokens.size()){
            num_string = tokens[j+1];
            break;
          }
          else{ break; }
        }
      }
      int num=AUROSTD_MAX_INT;
      if(!num_string.empty()){
        if(num_string.find(".") != std::string::npos){ //DX20200709 - aurostd::substring2bool to find
          vector<string> sub_tokens;
          aurostd::string2tokens(num_string,sub_tokens,".");
          num_string = sub_tokens[0];
        }
        num = aurostd::string2utype<int>(num_string);
      }
      if(num < min_num){
        min_ICSD = ICSD_entries[i];
        min_num = num;
      }
    }
    return min_ICSD;
  }
}

// ***************************************************************************
// XtalFinderCalculator::compareAtomDecorations() //DX20201215
// ***************************************************************************
void XtalFinderCalculator::compareAtomDecorations(
    StructurePrototype& structure,
    uint num_proc,
    bool optimize_match){

  aurostd::xoption permutation_options = compare::loadDefaultComparisonOptions("permutation");
  permutation_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",optimize_match);

  return compareAtomDecorations(structure,num_proc,permutation_options);
}

void XtalFinderCalculator::compareAtomDecorations(
    StructurePrototype& structure,
    uint num_proc,
    aurostd::xoption& permutation_options){

  string misfit_result = "";
  return compareAtomDecorations(structure,misfit_result,num_proc,permutation_options);
}

void XtalFinderCalculator::compareAtomDecorations(
    StructurePrototype& structure,
    string& misfit_results,
    uint num_proc,
    aurostd::xoption& permutation_options){

  // Compare the atom decorations for a given structure.
  // Calculates the symmetry, if not already calculated (hence pass by
  // non-const reference).

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  bool VERBOSE=false;
  stringstream message;

  // ---------------------------------------------------------------------------
  // options for atom decoration (permutation) comparisons
  bool same_species = true; // permutation comparisons must compare the same species
  bool quiet = true; //true

  vector<StructurePrototype> final_permutations;

  // ---------------------------------------------------------------------------
  // get stoichiometry
  vector<uint> stoichiometry = structure.structure_representative->structure.GetReducedComposition();

  // ---------------------------------------------------------------------------
  // calculate symmetry (if not already calculated)
  bool ignore_symmetry = permutation_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY");
  if(!ignore_symmetry && !isSymmetryCalculated(*structure.structure_representative)){
    calculateSymmetry(*structure.structure_representative);
  }

  // ---------------------------------------------------------------------------
  // check if permutations of the structure are compatible via
  // stoichiometry and/or symmetry, if not, return early
  if(!compare::arePermutationsComparableViaComposition(stoichiometry) ||
      (!ignore_symmetry && !compare::arePermutationsComparableViaSymmetry(structure.structure_representative->grouped_Wyckoff_positions))){

    vector<string> unique_permutations = getSpeciesPermutedStrings(stoichiometry); //DX20191125
    // store permutation results in main StructurePrototype object
    for(uint j=0;j<unique_permutations.size();j++){
      vector<string> permutations_tmp; permutations_tmp.push_back(unique_permutations[j]);
      structure.atom_decorations_equivalent.push_back(permutations_tmp);
    }
    return;
  }

  // ---------------------------------------------------------------------------
  // get LFA environments
  if(!isLFAEnvironmentCalculated(*structure.structure_representative)){
    computeLFAEnvironment(*structure.structure_representative);
  }

  // ---------------------------------------------------------------------------
  // generate all permutation structures
  generateAtomPermutedStructures(*structure.structure_representative);

  // ---------------------------------------------------------------------------
  // store decoration names
  vector<string> decoration_names;
  for(uint i=0;i<structure_containers.size();i++){
    decoration_names.push_back(structure_containers[i].name);
  }

  // ---------------------------------------------------------------------------
  // loop over grouping modes
  // atom permutations must obey Lagrange's theorem (divisor theory)
  // initial atom groupings by LFA environment are fast, but may result in
  // groupings that violate divisor theory; if this is the case we need to
  // regroup without considering the LFA environment
  //
  // mode=0: use LFA environment to filter
  // mode=1: if incommensurate groupings, then ignore LFA environment in grouping
  for(uint mode=0;mode<2;mode++){
    if(mode==0){
      if(!quiet || LDEBUG){
        message << "Considering environment analysis in grouping permutations (mode=0)." << endl;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
    }
    if(mode==1){
      //ignore_environment=true;
      permutation_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      permutation_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES",TRUE); //DX20200320
      if(!quiet || LDEBUG){
        message << "Could not find commensurate pemutations when grouping via environment. Ignoring environment analysis in grouping permutations (mode=1)." << endl;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
    }
    final_permutations.clear();

    // ---------------------------------------------------------------------------
    // group comparable permutations
    vector<StructurePrototype> permutation_comparisons;
    permutation_comparisons = groupStructurePrototypes(same_species,
        permutation_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"),
        permutation_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"),
        permutation_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"),
        permutation_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320
        false,
        quiet); //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // [COME BACK DX] ensure the representative stucture is an even permutation
    makeRepresentativeEvenPermutation(permutation_comparisons, decoration_names);

    if(VERBOSE){ for(uint i=0;i<permutation_comparisons.size();i++){ cerr << "Initial permutation groupings: " << permutation_comparisons[i] << endl; } }

    // ---------------------------------------------------------------------------
    // compare permutations
    final_permutations = runComparisonScheme(
        permutation_comparisons,
        same_species,
        num_proc,
        permutation_options,
        quiet); //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // check if matched permutations are physically possible
    if(!compare::checkNumberOfGroupings(final_permutations, decoration_names.size())){
      if(!quiet || LDEBUG){
        message << "Compared groupings of permutations do not follow number theory (# unique=" << final_permutations.size() << " vs # total=" << decoration_names.size() << ")" << endl;
        // comprehensive output
        if(LDEBUG){
          for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i] << endl; }
        }
        // minimal output
        else{
          for(uint i=0;i<final_permutations.size();i++){ //DX20191218 - new misfit struct version
            message << final_permutations[i].structure_representative->name << " = ";
            vector<string> duplicate_info;
            for(uint j=0;j<final_permutations[i].structures_duplicate.size();j++){
              stringstream ss_tmp; ss_tmp << final_permutations[i].structures_duplicate[j]->name << " (misfit_duplicate=" << final_permutations[i].mapping_info_duplicate[j].misfit << ")";
              duplicate_info.push_back(ss_tmp.str());
            }
            message << aurostd::joinWDelimiter(duplicate_info,",") << endl;
          }
        }
        message << "Trying to check if duplicates match better with other representative structures ... " << endl;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }

      // ---------------------------------------------------------------------------
      // check if better matchings; perhaps matched structures would have smaller misfits if matched to different representatives
      final_permutations = checkForBetterMatches(final_permutations,
          num_proc,
          true, same_species,
          permutation_options,
          quiet); //DX20200103 - condensed booleans to xoptions

      // ---------------------------------------------------------------------------
      // check if NEW matched permutations are physically possible
      if(!compare::checkNumberOfGroupings(final_permutations, decoration_names.size())){
        message << "Compared groupings of permutations do not follow number theory (# unique=" << final_permutations.size() << " vs # total=" << decoration_names.size() << ")" << endl;
        // comprehensive output
        if(LDEBUG){
          for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i] << endl; }
        }
        // minimal output
        else{
          //DX20191218 [ORIG] for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i].structure_representative_name << " = " << aurostd::joinWDelimiter(final_permutations[i].structures_duplicate_names,",") << " (misfits_duplicate=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(final_permutations[i].misfits_duplicate,8,true),",") << ")" << endl; }
          for(uint i=0;i<final_permutations.size();i++){ //DX20191218 - new misfit struct version
            message << final_permutations[i].structure_representative->name << " = ";
            vector<string> duplicate_info;
            for(uint j=0;j<final_permutations[i].structures_duplicate.size();j++){
              stringstream ss_tmp; ss_tmp << final_permutations[i].structures_duplicate[j]->name << " (misfit_duplicate=" << final_permutations[i].mapping_info_duplicate[j].misfit << ")";
              duplicate_info.push_back(ss_tmp.str());
            }
            message << aurostd::joinWDelimiter(duplicate_info,",") << endl;
          }
        }
        if(mode==1){  // exhausted checks
          message << "Please email aflow@groups.io and provide the corresponding example." << endl;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
      }
      else{ break; }
    }
    else{ break; }
  }

  if(VERBOSE){ for(uint i=0;i<final_permutations.size();i++){ cerr << "Final permutation groupings: " << final_permutations[i] << endl; } }

  // ---------------------------------------------------------------------------
  // store results
  for(uint j=0;j<final_permutations.size();j++){
    vector<string> permutations_tmp;
    permutations_tmp.push_back(final_permutations[j].structure_representative->name); //push back representative permutation
    for(uint d=0;d<final_permutations[j].structures_duplicate.size();d++){ permutations_tmp.push_back(final_permutations[j].structures_duplicate[d]->name); } //push back equivalent permutations
    structure.atom_decorations_equivalent.push_back(permutations_tmp);
  }

  // ---------------------------------------------------------------------------
  // print format
  filetype format = txt_ft;
  if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
    format = txt_ft;
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
    format = json_ft;
  }
  misfit_results = printResults(final_permutations, true, format);

  //DX20210115 [OBSOLETE] return final_permutations;
}

// ***************************************************************************
// XtalFinderCalculator::generateAtomPermutedStructures()
// ***************************************************************************
void XtalFinderCalculator::generateAtomPermutedStructures(
    structure_container& structure){

  // Generate all atom permuted variants of a given structure.
  // Moved Heap's algorithm into aurostd::xcombos.

  deque<string> species;
  bool is_symmetry_calculated = isSymmetryCalculated(structure);

  // ---------------------------------------------------------------------------
  // get permuted species strings via Heap's algorithm:
  // swap lowest position index first (left-most)
  vector<string> species_permuted = getSpeciesPermutedStrings(structure.stoichiometry);

  vector<string> names = pflow::getFakeElements(structure.stoichiometry.size()); //DX20200728 - now in pflow
  vector<uint> indices = structure.stoichiometry;
  uint num_elements = structure.stoichiometry.size();

  vector<vector<int> > all_indices;
  vector<int> _indices;
  for(uint i=0;i<num_elements;i++){_indices.push_back(i);}

  // ---------------------------------------------------------------------------
  // use Heap's algorithm: swap lowest position index first (left-most)
  // this is the preferred order for the representative atom decorations
  aurostd::xcombos indices_combos(_indices, true, 'P', heap_alg_xcombos);
  while (indices_combos.increment()) all_indices.push_back(indices_combos.getCombo());

  // ---------------------------------------------------------------------------
  // create permuted structure
  for(uint i=0;i<all_indices.size();i++){
    xstructure xstr_tmp = structure.structure;
    species.clear();
    for(uint j=0;j<all_indices[i].size();j++){
      species.push_back(names[all_indices[i][j]]);
    }
    xstr_tmp.SetSpecies(species);
    std::stable_sort(xstr_tmp.atoms.begin(),xstr_tmp.atoms.end(),sortAtomsNames); //DX20230227 - with updated SetSpecies using ReplaceAtoms, need to resort atoms based on names
    xstr_tmp.ReScale(1.0); //DX20190715
    xstr_tmp.BringInCell(); //DX20200707

    string system_name = aurostd::joinWDelimiter(species,"");
    stringstream ss_str; ss_str << "permutation of: " << xstr_tmp; //DX20190730
    addStructure2container(xstr_tmp,system_name,ss_str.str(),0,false);

    // for now compute, perhaps we can use index swap
    computeLFAEnvironment(structure_containers.back());

    if(is_symmetry_calculated){
      vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
      grouped_Wyckoff_positions = structure.grouped_Wyckoff_positions;
      for(uint j=0;j<grouped_Wyckoff_positions.size();j++){
        grouped_Wyckoff_positions[j].element=names[all_indices[i][j]];
      }
      structure_containers.back().grouped_Wyckoff_positions = grouped_Wyckoff_positions;
      structure_containers.back().space_group = structure.space_group;
    }

    // copy over neighbor distances (it will be the same as the parent structure)
    // CANNOT DO OUTRIGHT (update neighbor species). Need to think about this more carefully ...
    // structure_containers.back().nearest_neighbor_distances = structure.nearest_neighbor_distances;
  }
}

// ***************************************************************************
// XtalFinderCalculator::getSpeciesPermutedStrings() //DX20201222
// ***************************************************************************
// deque input
vector<string> XtalFinderCalculator::getSpeciesPermutedStrings(
    const deque<uint>& stoichiometry){
  vector<uint> stoichiometry_vstring = aurostd::deque2vector(stoichiometry);
  return getSpeciesPermutedStrings(stoichiometry_vstring);
}

// vector input
vector<string> XtalFinderCalculator::getSpeciesPermutedStrings(
    const vector<uint>& stoichiometry){

  vector<string> species, species_permuted;

  vector<string> names = pflow::getFakeElements(stoichiometry.size()); //DX20200728 - now in pflow
  vector<uint> indices = stoichiometry;
  uint num_elements = stoichiometry.size();

  vector<vector<int> > all_indices;
  vector<int> _indices;
  for(uint i=0;i<num_elements;i++){_indices.push_back(i);}

  // ---------------------------------------------------------------------------
  // use Heap's algorithm: swap lowest position index first (left-most)
  // this is the preferred order for the representative atom decorations
  aurostd::xcombos indices_combos(_indices, true, 'P', heap_alg_xcombos);
  while (indices_combos.increment()) all_indices.push_back(indices_combos.getCombo());

  // ---------------------------------------------------------------------------
  // create permuted species strings
  for(uint i=0;i<all_indices.size();i++){
    species.clear();
    for(uint j=0;j<all_indices[i].size();j++){
      species.push_back(names[all_indices[i][j]]);
    }
    species_permuted.push_back(aurostd::joinWDelimiter(species,""));
  }

  return species_permuted;
}

// ***************************************************************************
// arePermutationsComparableViaComposition()
// ***************************************************************************
// Check if permutations of a structure are possible via stoichiometry
// i.e., check if stoichiometric ratio value occurs more than once

// xstructure version
namespace compare{
  bool arePermutationsComparableViaComposition(const xstructure& xstr){

    bool matchable_sites=false;
    for(uint i=0;i<xstr.num_each_type.size();i++){
      for(uint j=i+1;j<xstr.num_each_type.size();j++){
        if(xstr.num_each_type[i]==xstr.num_each_type[j]){ matchable_sites=true; break; }
      }
      if(matchable_sites){ break; }
    }
    return matchable_sites;
  }
}

// vector<uint> stoichiometry version
namespace compare{
  bool arePermutationsComparableViaComposition(const vector<uint>& composition, bool reduce_composition){

    // ---------------------------------------------------------------------------
    // reduce stoichiometry first if necessary
    vector<uint> tmp_composition=composition;
    if(reduce_composition){ aurostd::reduceByGCD(composition, tmp_composition); } //DX20191125
    else{ tmp_composition = composition; } // assuming it is already reduced

    // ---------------------------------------------------------------------------
    // check if stoichiometries occur more than once
    bool matchable_sites=false;
    for(uint i=0;i<tmp_composition.size();i++){
      for(uint j=i+1;j<tmp_composition.size();j++){
        if(tmp_composition[i]==tmp_composition[j]){ matchable_sites=true; break; }
      }
      if(matchable_sites){ break; }
    }
    return matchable_sites;
  }
}

// ***************************************************************************
// arePermutationsComparableViaSymmetry()
// ***************************************************************************
namespace compare{
  bool arePermutationsComparableViaSymmetry(const vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions){

    // Check if permutations of a structure are possible via symmetry
    // (Wyckoff positions) i.e., check if there are matchable Wyckoff
    // positions in a given structure

    // ---------------------------------------------------------------------------
    // re-write site symmetries to account for cell choice differences
    // note: the resulting site symmetries may NOT be physical
    // (this is just a speed up; fast way to compare)
    vector<GroupedWyckoffPosition> sorted_temp_grouped_Wyckoffs = compare::sortSiteSymmetryOfGroupedWyckoffPositions(grouped_Wyckoff_positions);

    // ---------------------------------------------------------------------------
    // check if Wyckoff site symmetries and multiplicites occur more than once
    for(uint i=0;i<sorted_temp_grouped_Wyckoffs.size();i++){
      string set_1_site_symmetry = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[i].site_symmetries,",");
      string set_1_multiplicity = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[i].multiplicities,",");
      for(uint j=i+1;j<sorted_temp_grouped_Wyckoffs.size();j++){
        string set_2_site_symmetry = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[j].site_symmetries,",");
        string set_2_multiplicity = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[j].multiplicities,",");
        if(set_1_site_symmetry == set_2_site_symmetry && set_1_multiplicity == set_2_multiplicity){
          return true;
        }
      }
    }
    return false;
  }
}


// ***************************************************************************
// XtalFinderCalculator::addAFLOWPrototypes2container()
// ***************************************************************************
void XtalFinderCalculator::addAFLOWPrototypes2container(
    const vector<string>& vlabel){

  // Add the AFLOW prototypes to the vector of structure_container objects
  // Note: The AFLOW labels should already be filtered to the relevant
  // strutures for comparison.
  // This does NOT generate the AFLOW prototype structure; it only stores
  // the structural information. It will be created just prior to the
  // comparison routine (more efficient, since we may not end up comparing all
  // of the prototype structures).

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  for(uint i=0;i<vlabel.size();i++){
    if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Storing AFLOW prototype information for " << vlabel[i] << endl; }

    // anrl prototypes
    vector<string> tokens;
    if(aurostd::string2tokens(vlabel[i],tokens,"_")>=4) {
      vector<string> vparameter_values = anrl::getANRLParameters(vlabel[i],"all");
      for(uint j=0;j<vparameter_values.size();j++){

        structure_container structure_tmp;
        //if one degree of freedom, then no number scheme required
        if(aurostd::string2tokens(vparameter_values[j],tokens,",")==1){
          structure_tmp.name = vlabel[i];
        }
        //if multiple degrees of freedom, then number scheme is required
        else {
          stringstream tmp; tmp << vlabel[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
          structure_tmp.name = tmp.str();
        }
        structure_tmp.stoichiometry = structure_containers[0].stoichiometry;
        structure_tmp.space_group = structure_containers[0].space_group; // same as representative structure (either will be the same, or we are forcing it to be for the ignore_symmetry/ignore_Wyckoff run)
        structure_tmp.grouped_Wyckoff_positions = structure_containers[0].grouped_Wyckoff_positions;
        structure_tmp.elements = pflow::getFakeElements(structure_containers[0].stoichiometry.size()); //DX20200728 - now in pflow
        structure_tmp.compound = pflow::prettyPrintCompound(structure_tmp.elements,structure_tmp.stoichiometry,no_vrt,true,txt_ft); //remove ones is true
        structure_tmp.is_structure_generated = false; // will be generated on-the-fly (faster)
        structure_tmp.source = "aflow_prototypes";
        structure_tmp.relaxation_step = 0; //DX20200429 - prototypes are unrelaxed
        structure_tmp.number_compounds_matching_structure = 0;
        structure_containers.push_back(structure_tmp);
      }
    }
    // htqc prototypes
    else {
      structure_container structure_tmp;
      structure_tmp.name = vlabel[i];
      structure_tmp.stoichiometry = structure_containers[0].stoichiometry;
      structure_tmp.space_group = structure_containers[0].space_group; // same as representative structure (either will be the same, or we are forcing it to be for the ignore_symmetry/ignore_Wyckoff run)
      structure_tmp.grouped_Wyckoff_positions = structure_containers[0].grouped_Wyckoff_positions;
      structure_tmp.elements = pflow::getFakeElements(structure_containers[0].stoichiometry.size()); //DX20200728 - now in pflow
      structure_tmp.compound = pflow::prettyPrintCompound(structure_tmp.elements,structure_tmp.stoichiometry,no_vrt,true,txt_ft); //remove ones is true
      structure_tmp.is_structure_generated = false; // will be generated on-the-fly
      structure_tmp.source = "aflow_prototypes";
      structure_tmp.relaxation_step = 0; //DX20200429 - prototypes are unrelaxed
      structure_tmp.number_compounds_matching_structure = 0;
      structure_containers.push_back(structure_tmp);
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::addDatabaseEntry2container
// ***************************************************************************
void XtalFinderCalculator::addDatabaseEntry2container(
    aflowlib::_aflowlib_entry& entry,
    const vector<string>& species,
    uint relaxation_step,
    bool same_species){

  // Add the AFLOW entry to the vector of structure containers objects
  // Note: this may change with new AFLUX integration.

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;  

  // ---------------------------------------------------------------------------
  // load xstructures from pflow::loadXstructures()
  bool load_most_relaxed_structure_only = false;
  if(relaxation_step==_COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){ load_most_relaxed_structure_only = true; }

  vector<string> structure_files;
  if(!pflow::loadXstructures(entry,structure_files,*p_FileMESSAGE,*p_oss,load_most_relaxed_structure_only)){
    message << "Could not load structure (auid=" << entry.auid << ") ... skipping...";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  // ---------------------------------------------------------------------------
  // ensure the correct relaxation type was extracted
  if(load_most_relaxed_structure_only && entry.vstr.size()!=1){
    message << "Expected only one structure to be loaded (the most relaxed structure). " << entry.vstr.size() << " were returned. Unexpected code behavior.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  uint structure_index = 0;
  if(!load_most_relaxed_structure_only){
    if(entry.vstr.size()==3 && structure_files.size()==3){
      if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ &&
          (structure_files[0] == "POSCAR.orig" ||
           structure_files[0] == "POSCAR.relax1")){
        structure_index = 0;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded original structure: " << structure_files[0] << endl;}
      }
      else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ &&
          (structure_files[1] == "POSCAR.relax2" ||
           structure_files[1] == "CONTCAR.relax1")){
        structure_index = 1;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded relax1 structure: " << structure_files[1] << endl;}
      }
      else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_ &&
          (structure_files[2] == "POSCAR.relax2" ||
           structure_files[2] == "CONTCAR.relax1")){
        structure_index = 2;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded most relaxed structure: " << structure_files[2] << endl;}
      }
      else{
        message << "Unexpected file names: " << aurostd::joinWDelimiter(structure_files, ",") << ". Cannot verify relaxation based on filename.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
      }
    }
    else{
      message << "Expected three structures to be loaded (original, relax1, and the most relaxed). " << entry.vstr.size() << " were returned. Unexpected code behavior.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
    }
  }

  // ---------------------------------------------------------------------------
  // store entry from database into the structure_container object
  deque<string> deque_species = aurostd::vector2deque(species);
  entry.vstr[structure_index].SetSpecies(deque_species);
  string str_path = entry.getPathAURL(*p_FileMESSAGE,*p_oss,false);
  addStructure2container(entry.vstr[structure_index], str_path, "aurl", relaxation_step, same_species);
}

// ***************************************************************************
//  Stoichiometry, Element, and GCD functions
// ***************************************************************************
// DX20200728 [OBSOLETE - MOVED OR CONSOLIDATED INTO AUROSTD FUNCTIONS]

// ***************************************************************************
// XtalFinderCalculator::splitComparisonIntoThreads()
// ***************************************************************************
bool XtalFinderCalculator::splitComparisonIntoThreads(
    vector<StructurePrototype>& comparison_schemes,
    uint num_proc,
    vector<std::pair<uint,uint> >& start_indices,
    vector<std::pair<uint,uint> >& end_indices){

  // Partitioning comparisons onto threads.
  // This is splitting on a 2D-array (different than getThreadDistribution).

  // ---------------------------------------------------------------------------
  // split comparisons into threads via indices
  stringstream message;
  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  bool safety_check=false; // safety check if split incorrectly

  uint number_of_comparisons = 0;
  for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
  //cerr << "# of comparisons: " << number_of_comparisons << endl;

  if(number_of_comparisons==0){
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Number of comparisons is zero, no need to split into threads." << endl;
    }
    return true;
  }

  uint num_per_thread = number_of_comparisons/num_proc;
  uint residual = number_of_comparisons%num_proc;
  bool accounted_for_residual=false;
  if(residual!=0){num_per_thread+=1;}

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Number of comparisons per thread: " << num_per_thread << endl;
  }

  //[CO20220625 - not used]uint tmp =0;

  uint count = 0;
  uint thread_count = 0;
  std::pair<uint,uint> tmp_start, tmp_end;
  std::pair<uint,uint> indices;
  for(uint i=0;i<comparison_schemes.size();i++){
    //DEBUG cerr << "splitting comparison indices: i: " << i << "/" << comparison_schemes.size() << endl;
    for(uint j=0;j<comparison_schemes[i].structures_duplicate.size();j++){
      indices.first=i, indices.second=j;
      count+=1;
      //[CO20220625 - not used]tmp+=1;
      if(count == num_per_thread && thread_count<num_proc-1){
        thread_count+=1;
        start_indices.push_back(tmp_start);
        //update tmp_start
        if(j+1>=comparison_schemes[i].structures_duplicate.size()-1 && i+1<comparison_schemes.size()-1){
          tmp_start.first=i+1; tmp_start.second=0;
          tmp_end.first=i+1; tmp_end.second=0;
        }
        else {
          tmp_start.first=i; tmp_start.second=j+1;
          tmp_end.first=i; tmp_end.second=j+1;
        }
        end_indices.push_back(tmp_end);
        count = 0;
      }
      else if(thread_count==num_proc-1 && i==comparison_schemes.size()-1 && j==comparison_schemes[i].structures_duplicate.size()-1){
        thread_count+=1;
        start_indices.push_back(tmp_start);
        //update tmp_start
        if(j+1>=comparison_schemes[i].structures_duplicate.size()-1 && i+1<comparison_schemes.size()-1){
          tmp_start.first=i+1; tmp_start.second=0;
          tmp_end.first=i+1; tmp_end.second=0;
        }
        else {
          tmp_start.first=i; tmp_start.second=j+1;
          tmp_end.first=i; tmp_end.second=j+1;
        }
        end_indices.push_back(tmp_end);
        count = 0;
      }
      if(!accounted_for_residual && residual!=0 && thread_count==residual){
        accounted_for_residual=true;
        num_per_thread=num_per_thread-1;
      }
    }
  }
  // need to add last if it was not already added, using "indices" variable to populate last indices
  if(count){
    start_indices.push_back(tmp_start);
    //tmp_end.first=indices.first; tmp_end.second=indices.second; //+1
    tmp_end.first=comparison_schemes.size()-1; tmp_end.second=comparison_schemes[tmp_end.first].structures_duplicate.size();
    end_indices.push_back(tmp_end);
  }

  // ---------------------------------------------------------------------------
  // check if split correctly
  // in case the number of threads is greater than the number of xstructures to test
  // put inside an if-statement (on 20190715) to save time; the function has been well-tested
  if(safety_check){
    uint recovered=0;
    uint num_of_threads=0;
    if(start_indices.size()>=num_proc){
      num_of_threads=num_proc;
    }
    else if(start_indices.size()<num_proc){
      num_of_threads=start_indices.size();
    }

    for(uint n=0; n<num_of_threads; n++){
      uint i_min=start_indices[n].first; uint i_max=end_indices[n].first;
      uint j_min=0; uint j_max=0;
      for(uint i=0;i<comparison_schemes.size();i++){
        // to loop properly
        if(i==i_min){
          j_min=start_indices[n].second;
          if(i==i_max){j_max=end_indices[n].second;}
          else {j_max=comparison_schemes[i].structures_duplicate.size();} //-1 since in loop: j<=j_max
        }
        else if(i==i_max){j_min=0; j_max=end_indices[n].second;}
        else {j_min=0; j_max=comparison_schemes[i].structures_duplicate.size();} //-1 since in loop: j<=j_max
        for(uint j=0;j<comparison_schemes[i].structures_duplicate.size();j++){
          if(i>=i_min && j>=j_min &&
              i<=i_max && j<j_max){
            recovered+=1;
          }
        }
      }
    }
    if(recovered != number_of_comparisons){
      message << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << number_of_comparisons;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }
  return true;
}

// ***************************************************************************
// XtalFinderCalculator::performStructureConversions() //DX20210113
// ***************************************************************************
void XtalFinderCalculator::performStructureConversions(
    uint i,
    const vector<bool>& calculate_primitive_vec,
    const vector<bool>& calculate_Minkowski_vec,
    const vector<bool>& calculate_Niggli_vec){

  // Perform primitivizations, Minkowski reductions, or Niggli reductions
  // on the structures indicated by the corresponding vector of booleans

  // ---------------------------------------------------------------------------
  // primitivize
  if(calculate_primitive_vec[i]){
    structure_containers[i].structure.GetPrimitive();
    structure_containers[i].natoms = structure_containers[i].structure.atoms.size(); //DX20210316 - updated number of atoms
  }
  // ---------------------------------------------------------------------------
  // Minkowski
  if(calculate_Minkowski_vec[i]){ structure_containers[i].structure.MinkowskiBasisReduction(); }
  // ---------------------------------------------------------------------------
  // Niggli
  if(calculate_Niggli_vec[i]){ structure_containers[i].structure.NiggliUnitCellForm(); }

}

// ***************************************************************************
// XtalFinderCalculator::convertStructures()
// ***************************************************************************
void XtalFinderCalculator::convertStructures(
    const aurostd::xoption& comparison_options,
    uint num_proc){

  // Converts all structures to certain unit cell representations
  // generally offering a speed increase for comparisons.
  // The unit cell conversion routines are in aflow_xatom; this function
  // is a wrapper for the XtalFinderCalculator class.

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  message << "Converting structures standard representation (primitive, Minkowski, and/or Niggli).";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  uint number_of_structures = structure_containers.size();

  // ---------------------------------------------------------------------------
  // primitivize (do this first)
  vector<bool> calculate_primitive_vec;
  calculate_primitive_vec.resize(number_of_structures,false);
  if(comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE")){
    for(uint i=0;i<number_of_structures;i++){
      calculate_primitive_vec[i] = !structure_containers[i].structure.primitive_calculated;
    }
  }

  // ---------------------------------------------------------------------------
  // Minkowski (second)
  vector<bool> calculate_Minkowski_vec;
  calculate_Minkowski_vec.resize(number_of_structures,false);
  if(comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI")){
    for(uint i=0;i<number_of_structures;i++){
      calculate_Minkowski_vec[i] = !structure_containers[i].structure.Minkowski_calculated;
    }
  }

  // ---------------------------------------------------------------------------
  // Niggli (last)
  vector<bool> calculate_Niggli_vec;
  calculate_Niggli_vec.resize(number_of_structures,false);
  if(comparison_options.flag("COMPARISON_OPTIONS::NIGGLI")){
    for(uint i=0;i<number_of_structures;i++){
      calculate_Niggli_vec[i] = !structure_containers[i].structure.Niggli_calculated;
    }
  }

#ifdef AFLOW_MULTITHREADS_ENABLE
  // THREADED VERSION - START
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of threads=" << num_proc << endl;}

  xthread::xThread xt(num_proc);
  std::function<void(uint, const vector<bool>&, const vector<bool>&, const vector<bool>&)> fn =
    std::bind(&XtalFinderCalculator::performStructureConversions, this, _1, _2, _3, _4);
  xt.run(number_of_structures, fn, calculate_primitive_vec, calculate_Minkowski_vec, calculate_Niggli_vec);

  // THREADED VERSION - END
#else
  // NONTHREADS - START
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Non-threaded version. Number of threads=" << num_proc << endl;}

  // ---------------------------------------------------------------------------
  // perform the relevant structure conversions
  for (uint i = 0; i < number_of_structures; i++) {
    performStructureConversions(
        i,
        calculate_primitive_vec,
        calculate_Minkowski_vec,
        calculate_Niggli_vec);
  }
#endif
  message << "All structures converted.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
}

// ***************************************************************************
// XtalFinderCalculator::getPrimitiveStructures() //DX20201006
// ***************************************************************************
void XtalFinderCalculator::getPrimitiveStructures(uint start_index, uint end_index){

  // If end index is greater than structure_containers.size(), then compute
  // Primitive cell for all structures
  // Used for pre-distributed threads
  if(end_index > structure_containers.size()){ end_index=structure_containers.size(); }

  for(uint i=start_index;i<end_index;i++){
    structure_containers[i].structure.GetPrimitive();
  }
}

// ***************************************************************************
// XtalFinderCalculator::getMinkowskiStructures() //DX20201006
// ***************************************************************************
void XtalFinderCalculator::getMinkowskiStructures(uint start_index, uint end_index){

  // If end index is greater than structure_containers.size(), then compute
  // Minkowski cell for all structures
  // Used for pre-distributed threads
  if(end_index > structure_containers.size()){ end_index=structure_containers.size(); }

  for(uint i=start_index;i<end_index;i++){
    structure_containers[i].structure.MinkowskiBasisReduction();
  }
}

// ***************************************************************************
// XtalFinderCalculator::getNiggliStructures() //DX20201006
// ***************************************************************************
void XtalFinderCalculator::getNiggliStructures(uint start_index, uint end_index){

  // If end index is greater than structure_containers.size(), then compute
  // Niggli cell for all structures
  // Used for pre-distributed threads
  if(end_index > structure_containers.size()){ end_index=structure_containers.size(); }

  for(uint i=start_index;i<end_index;i++){
    structure_containers[i].structure.NiggliUnitCellForm();
  }
}


// ***************************************************************************
// compare::splitTaskIntoThreads()
// ***************************************************************************
//DX20191108 [OBSOLETE - switching to getThreadDistribution]

// ***************************************************************************
// XtalFinderCalculator::calculateSpaceGroups()
// ***************************************************************************
void XtalFinderCalculator::calculateSpaceGroups(uint start_index, uint end_index, uint setting){ //DX20191230 - added setting

  // Calculates the space group and Wyckoff positions for the structure in
  // the XtalFinderCalculator structure container.
  // Wraps around SYM::calculateSpaceGroups(), and it is specific for
  // XtalFinderCalculator objects (populates relevant attributes)

  bool no_scan = false; //DX20191230

  // If end index is greater than structure_containers.size(), then compute
  // symmetry for all structures in the containers
  if(end_index > structure_containers.size()){ end_index=structure_containers.size(); }

  for(uint i=start_index;i<end_index;i++){ //DX20191107 - switching convention <= vs <
    // use sym_eps if given, otherwise use default
    if(!structure_containers[i].structure.sym_eps_calculated || structure_containers[i].structure.sym_eps==AUROSTD_NAN) { structure_containers[i].structure.sym_eps = SYM::defaultTolerance(structure_containers[i].structure); } //DX20210421
    double use_tol = structure_containers[i].structure.sym_eps; //DX20210421
    structure_containers[i].space_group = structure_containers[i].structure.SpaceGroup_ITC(use_tol, -1, setting, no_scan); //DX20191230 - added arguments
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    compare::groupWyckoffPositions(structure_containers[i].structure, grouped_Wyckoff_positions);
    structure_containers[i].grouped_Wyckoff_positions=grouped_Wyckoff_positions;
  }
}

// ***************************************************************************
// calculateSymmetries - Calculate Symmetries
// ***************************************************************************
void XtalFinderCalculator::calculateSymmetries(uint num_proc){

  // Calculates the symmetry (space group and Wyckoff positions) of each
  // structure

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  message << "Calculating the symmetries of the structure(s).";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

#ifdef AFLOW_MULTITHREADS_ENABLE
  // THREADED VERISON - START
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of threads=" << num_proc << endl;}

  // Distribute threads via indices
  uint number_of_structures = structure_containers.size();
  xthread::xThread xt(num_proc);
  std::function<void(uint, uint, uint)> fn = std::bind(&XtalFinderCalculator::calculateSpaceGroups, this, _1, _2, _3);
  uint setting = SG_SETTING_ANRL;
  xt.runPredistributed(number_of_structures, fn, setting);
  // THREADED VERISON - END

#else

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Non-threaded version: " << num_proc << endl;}

  // NON-THREADED VERSION - START
  for(uint i=0; i<structure_containers.size(); i++){
    calculateSymmetry(structure_containers[i]);
  }
  // NON-THREADED VERSION - END

#endif

  message << "Symmetries calculated.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);

}

// ***************************************************************************
// XtalFinderCalculator::setSymmetryPlaceholders()
// ***************************************************************************
void XtalFinderCalculator::setSymmetryPlaceholders(){

  // If the symmetry analysis is suppressed, then set the symmetries
  // attributes to their default values

  for(uint i=0;i<structure_containers.size();i++){
    structure_containers[i].Pearson = "xX";
    structure_containers[i].space_group = 0;
    vector<GroupedWyckoffPosition> vGWyckoffPos_tmp;
    structure_containers[i].grouped_Wyckoff_positions = vGWyckoffPos_tmp;
  }

}

// ***************************************************************************
// XtalFinderCalculator::calculateLFAEnvironments()
// ***************************************************************************
void XtalFinderCalculator::computeLFAEnvironments(uint start_index, uint end_index){

  // Calculates the LFA environments for a structure and
  // stores it in the structure container

  // If end index is greater than structure_containers.size(), then compute
  // the LFA environment analysis for all structures
  if(end_index > structure_containers.size()){ end_index=structure_containers.size(); }

  for(uint i=start_index;i<end_index;i++){ //DX20191107 switching end index convention <= vs <
    computeLFAEnvironment(structure_containers[i]);
  }

}

// ***************************************************************************
// XtalFinderCalculator::calculateLFAEnvironments()
// ***************************************************************************
void XtalFinderCalculator::calculateLFAEnvironments(uint num_proc){

  // Calculates the LFA environments for a structure and
  // stores it in the structure container

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  message << "Calculating the environments of the structure(s).";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

#ifdef AFLOW_MULTITHREADS_ENABLE
  // THREADED VERISON - START

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of threads=" << num_proc << endl;}

  // Distribute threads via indices
  uint number_of_structures = structure_containers.size();
  xthread::xThread xt(num_proc);
  std::function<void(uint, uint)> fn = std::bind(&XtalFinderCalculator::computeLFAEnvironments, this, _1, _2);
  xt.runPredistributed(number_of_structures, fn);
  // THREADED VERISON - END

#else
  // NON-THREADED VERSION
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Non-threaded version: " << num_proc << endl;}
  computeLFAEnvironments(); //DX20191122 - for all structures

#endif
  message << "Environments calculated.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);

}

// ***************************************************************************
// XtalFinderCalculator::getNearestNeighbors()
// ***************************************************************************
void XtalFinderCalculator::getNearestNeighbors(uint num_proc){

  // Calculates the LFA environments for a structure and
  // stores it in the structure container

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  message << "Calculating the nearest neighbors of all the structure(s).";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

#ifdef AFLOW_MULTITHREADS_ENABLE
  // THREADED VERISON - START

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of threads=" << num_proc << endl;}
  // Distribute threads via indices
  uint number_of_structures = structure_containers.size();
  xthread::xThread xt(num_proc);
  std::function<void(uint, uint)> fn = std::bind(&XtalFinderCalculator::calculateNearestNeighbors, this, _1, _2);
  xt.runPredistributed(number_of_structures, fn);
  // THREADED VERISON - END

#else
  // NON-THREADED VERSION - START
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Non-threaded version: " << num_proc << endl;}
  for(uint i=0;i<structure_containers.size();i++){
    structure_containers[i].nearest_neighbor_distances = NearestNeighbors(structure_containers[i].structure); // nearest neighbor distances (invariant of origin shifts)
  }
  // NON-THREADED VERSION - END

#endif
  message << "Nearest neighbors information calculated.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);

}

// ***************************************************************************
// XtalFinderCalculator::calculateNearestNeighbors()
// ***************************************************************************
void XtalFinderCalculator::calculateNearestNeighbors(uint start_index, uint end_index){

  // Calculates the nearest neighbor information for a structure and
  // stores it in the structure container

  // If end index is greater than structure_containers.size(), then compute
  // nearest neighbor analysis for all structures
  if(end_index > structure_containers.size()){ end_index=structure_containers.size(); }

  for(uint i=start_index;i<end_index;i++){ //DX20191107 switching end index convention <= vs <
    structure_containers[i].nearest_neighbor_distances = NearestNeighbors(structure_containers[i].structure); // nearest neighbor distances (invariant of origin shifts)
  }
}

// ***************************************************************************
// compare::groupWyckoffPositions()
// ***************************************************************************
// Groups the Wyckoff positions via species/types
// Obtains information from xstructure
// Assumes xstr.SpaceGroup_ITC() has been called, otherwise, this will fail
namespace compare{
  void groupWyckoffPositions(const xstructure& xstr,
      vector<GroupedWyckoffPosition>& grouped_positions){

    groupWyckoffPositions(xstr.wyckoff_sites_ITC, grouped_positions);
  }
}

namespace compare{
  void groupWyckoffPositions(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC,
      vector<GroupedWyckoffPosition>& grouped_positions){

    uint type_count = 0; //DX20190425 - add type indicator

    for(uint i=0;i<wyckoff_sites_ITC.size();i++){
      uint multiplicity = wyckoff_sites_ITC[i].multiplicity; //DX20191031
      string letter = wyckoff_sites_ITC[i].letter; //DX20191031
      string site_symmetry = wyckoff_sites_ITC[i].site_symmetry; //DX20191031

      bool element_found = false;
      uint element_index = 0;
      for(uint j=0;j<grouped_positions.size();j++){
        if(KBIN::VASP_PseudoPotential_CleanName(wyckoff_sites_ITC[i].type) == KBIN::VASP_PseudoPotential_CleanName(grouped_positions[j].element)){ //DX20190329 - remove pseudopotential info
          element_found = true;
          element_index = j;
          break;
        }
      }
      if(element_found == false){
        GroupedWyckoffPosition tmp;
        tmp.type = type_count; //DX20190425 - added type
        tmp.element = KBIN::VASP_PseudoPotential_CleanName(wyckoff_sites_ITC[i].type); //DX20190329 - remove pseudopotential info
        tmp.site_symmetries.push_back(site_symmetry);
        tmp.multiplicities.push_back(multiplicity);
        tmp.letters.push_back(letter); //DX20190208 - add Wyckoff letters
        grouped_positions.push_back(tmp);
        type_count++; //DX20190425
      }
      else {
        grouped_positions[element_index].site_symmetries.push_back(site_symmetry);
        grouped_positions[element_index].multiplicities.push_back(multiplicity);
        grouped_positions[element_index].letters.push_back(letter); //DX20190208 - add Wyckoff letters
      }
    }

  }
}

// ***************************************************************************
// compare::groupWyckoffPositionsFromGroupedString()
// ***************************************************************************
namespace compare{
  void groupWyckoffPositionsFromGroupedString(
      uint space_group_number,
      uint setting,
      const vector<vector<string> >& grouped_Wyckoff_string,
      vector<GroupedWyckoffPosition>& grouped_positions){

    // Groups the Wyckoff positions via species
    // Obtains information from the string of the following form: a,f,g;b,c;a,a
    // (i.e., ";" separates by species, and "," splits the Wyckoff letters for that species)

    stringstream axis_cell;
    axis_cell.str(std::string());
    axis_cell << setting; //DX20180806 - use setting
    SymmetryInformationITC ITC_sym_info; //DX20190215
    ITC_sym_info.initsgs(axis_cell.str()); //DX20190215
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1]; //DX20190215
    for(uint i=0;i<grouped_Wyckoff_string.size();i++){
      GroupedWyckoffPosition tmp;
      tmp.type = i; //DX20200625
      tmp.element = aurostd::utype2string<uint>(i); //DX20200625 - this is just a placeholder since atoms were not fed in
      for(uint j=0;j<grouped_Wyckoff_string[i].size();j++){
        string Wyckoff_letter = grouped_Wyckoff_string[i][j];
        uint Wyckoff_multiplicity = 0;
        string Wyckoff_site_symmetry = "";
        vector<string> positions;
        SYM::get_Wyckoff_from_letter(spacegroupstring, Wyckoff_letter, Wyckoff_multiplicity, Wyckoff_site_symmetry, positions);
        tmp.site_symmetries.push_back(Wyckoff_site_symmetry);
        tmp.multiplicities.push_back(Wyckoff_multiplicity);
        tmp.letters.push_back(Wyckoff_letter);
      }
      grouped_positions.push_back(tmp);
    }
  }
}

// ***************************************************************************
// compare::printWyckoffString()
// ***************************************************************************
namespace compare{
  string printWyckoffString(
      const vector<GroupedWyckoffPosition>& grouped_positions,
      bool alphabetize){

    // Prints the grouped Wyckoff positions as a string of the following form: a,f,g;b,c;a,a
    // (i.e., ";" separates by species, and "," splits the Wyckoff letters for that species)

    vector<string> all_Wyckoff_sets;
    for(uint i=0;i<grouped_positions.size();i++){
      vector<string> Wyckoff_set;
      for(uint j=0;j<grouped_positions[i].letters.size();j++){
        Wyckoff_set.push_back(grouped_positions[i].letters[j]);
      }
      if(alphabetize){
        std::sort(Wyckoff_set.begin(),Wyckoff_set.end());
      }
      all_Wyckoff_sets.push_back(aurostd::joinWDelimiter(Wyckoff_set,","));
    }
    return aurostd::joinWDelimiter(all_Wyckoff_sets,";");
  }
}

// ***************************************************************************
// compare::sortSiteSymmetryOfGroupedWyckoffPositions()
// ***************************************************************************
namespace compare{
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(
      const vector<GroupedWyckoffPosition>& grouped_Wyckoffs){

    // Sort they Wyckoff site symmetry for each Wyckoff position.
    // It is possible to match Wyckoff positions with different ordering of
    // the site symmetry depending on the cell choice.
    // To account for this, we split the site symmetry into a vector
    // containing the symmetries along the primary, secondary, and tertiary
    // directions.  Then, we sort the vector (alphabetic will suffice,
    // just need a standard comparison method).
    // NOTE: Sorting is only for comparing site symmetry strings to negate
    // choice of directions. Sorting may yield "site symmetries" that are
    // not physical.

    vector<GroupedWyckoffPosition> sorted_site_symmetry_Wyckoff_positions;

    for(uint i=0;i<grouped_Wyckoffs.size();i++){
      sorted_site_symmetry_Wyckoff_positions.push_back(grouped_Wyckoffs[i]);
      for(uint j=0;j<grouped_Wyckoffs[i].site_symmetries.size();j++){
        vector<string> split_site_symmetry = SYM::splitSiteSymmetry(grouped_Wyckoffs[i].site_symmetries[j]);
        std::sort(split_site_symmetry.begin(),split_site_symmetry.end());
        sorted_site_symmetry_Wyckoff_positions[i].site_symmetries[j] = aurostd::joinWDelimiter(split_site_symmetry,"");
      }
    }

    return sorted_site_symmetry_Wyckoff_positions;
  }
}

// ***************************************************************************
// compare::matchableWyckoffPositions()
// ***************************************************************************
namespace compare{
  bool matchableWyckoffPositions(
      const vector<GroupedWyckoffPosition>& grouped_Wyckoffs_str1,
      const vector<GroupedWyckoffPosition>& grouped_Wyckoffs_str2,
      bool same_species){

    // Determines if two sets of grouped Wyckoff positions are commensurate
    // with one another, i.e., checks if the Wyckoff multiplicities and site
    // symmetries are the same.  Comparing the same species is optional.

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    // ---------------------------------------------------------------------------
    // quick check: are number of Wyckoff positions the same? cannot match otherwise
    if(grouped_Wyckoffs_str1.size() != grouped_Wyckoffs_str2.size()){
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " # of Wyckoff positions does not match ("
          << grouped_Wyckoffs_str1.size() << " vs " << grouped_Wyckoffs_str2.size() << endl;
      }
      return false;
    }

    // ---------------------------------------------------------------------------
    // sort site symmetries to account for different cell choices
    vector<GroupedWyckoffPosition> grouped_Wyckoffs_str1_sorted = compare::sortSiteSymmetryOfGroupedWyckoffPositions(grouped_Wyckoffs_str1);
    vector<GroupedWyckoffPosition> grouped_Wyckoffs_str2_sorted = compare::sortSiteSymmetryOfGroupedWyckoffPositions(grouped_Wyckoffs_str2);

    // ---------------------------------------------------------------------------
    // DEBUG: print sorted Wyckoff positions
    if(LDEBUG) {
      for(uint i=0;i<grouped_Wyckoffs_str1_sorted.size();i++){ cerr << "str1: " << grouped_Wyckoffs_str1_sorted[i] << endl; }
      for(uint i=0;i<grouped_Wyckoffs_str2_sorted.size();i++){ cerr << "str2: " << grouped_Wyckoffs_str2_sorted[i] << endl; }
    }

    // ---------------------------------------------------------------------------
    // check if multiplicities and site symmetries are commensurate
    vector<bool> matched_species_str2;
    matched_species_str2.resize(grouped_Wyckoffs_str2_sorted.size(), false); // set all to false
    vector<bool> Wyckoff_subset_matches;

    bool match_set_str1 = false, match_single_str1 = false;
    uint match_count = 0;
    for(uint i=0;i<grouped_Wyckoffs_str1_sorted.size();i++){
      match_set_str1 = false;
      for(uint j=0;j<grouped_Wyckoffs_str2_sorted.size();j++){
        // ---------------------------------------------------------------------------
        // check if species in str2 have already been matched
        if(matched_species_str2[j]){ continue; }
        // ---------------------------------------------------------------------------
        // check if multiplicities are the same size
        // if considering species, check the element names
        if(((same_species && grouped_Wyckoffs_str1_sorted[i].element == grouped_Wyckoffs_str2_sorted[j].element) || !same_species) && // check species requirement
            grouped_Wyckoffs_str1_sorted[i].multiplicities.size() == grouped_Wyckoffs_str2_sorted[j].multiplicities.size()){ // check multiplicities
          // ---------------------------------------------------------------------------
          // loop through multiplicies and ensure each Wyckoff position is matched to
          // only once (via vector<bool> Wyckoff_subset_matches)
          match_count = 0;
          Wyckoff_subset_matches.clear(); Wyckoff_subset_matches.resize(grouped_Wyckoffs_str2_sorted[j].multiplicities.size(), false);
          for(uint m=0;m<grouped_Wyckoffs_str1_sorted[i].multiplicities.size();m++){
            match_single_str1 = false;
            for(uint n=0;n<grouped_Wyckoffs_str2_sorted[j].multiplicities.size();n++){
              if(grouped_Wyckoffs_str1_sorted[i].multiplicities[m] == grouped_Wyckoffs_str2_sorted[j].multiplicities[n] &&
                  grouped_Wyckoffs_str1_sorted[i].site_symmetries[m] == grouped_Wyckoffs_str2_sorted[j].site_symmetries[n] &&
                  !Wyckoff_subset_matches[n]){
                match_count++;
                match_single_str1 = true;
                Wyckoff_subset_matches[n] = true;
                break; // need to break to prevent a many-to-one mapping of Wyckoff sites
              }
            }
            if(!match_single_str1){ break; } // if no matches for any Wyckoff positions in str1, then we cannot match to this Wyckoff set
          }
          // ---------------------------------------------------------------------------
          // if any match, all need to match; otherwise the Wyckoff positions are not matchable
          // Note: last two if-statements determine if all of the Wyckoff subsets in str2
          // have been matched
          if(match_count == grouped_Wyckoffs_str1_sorted[i].multiplicities.size() && std::equal(Wyckoff_subset_matches.begin(),Wyckoff_subset_matches.end(), Wyckoff_subset_matches.begin()) && Wyckoff_subset_matches[0]){
            match_set_str1 = true;
            matched_species_str2[j] = true;
          }
        }
        // ---------------------------------------------------------------------------
        // if we've matched the ith Wyckoff position in str1, we do not need to
        // continue checking against the representative Wyckoff positions
        if(match_set_str1) { break; }
      }
      // ---------------------------------------------------------------------------
      // if we cannot match the ith Wyckoff position in str1, then the Wyckoff
      // sequences cannot be matched
      if(!match_set_str1) { return false; }
    }

    // ---------------------------------------------------------------------------
    // if all species in str2 have not been mapped, sequences do not match
    for(uint i=0;i<matched_species_str2.size();i++){
      if(!matched_species_str2[i]){ return false; }
    }

    return true;
  }
}

// ***************************************************************************
// compare::convertANRLWyckoffString2GroupedPositions()
// ***************************************************************************
namespace compare{
  vector<vector<string> > convertANRLWyckoffString2GroupedPositions(const string& label){

    // Converts the ANRL Wyckoff string to a grouped Wyckoff position object
    // ANRL Wyckoff string example: a2b_4bcd_e2f3g

    vector<vector<string> > anrl_grouped_Wyckoff_letters;
    vector<string> anrl_Wyckoff_set, tmp;
    aurostd::string2tokens(label, tmp, "_");
    for(uint i=0;i<tmp.size();i++){ if(i>2){ anrl_Wyckoff_set.push_back(tmp[i]); }}
    // DEBUGGING - ::print(anrl_Wyckoff_set);
    for(uint i=0;i<anrl_Wyckoff_set.size();i++){
      vector<string> Wyckoff_letters_for_species;
      uint Wyckoff_count = 1;
      bool is_previous_digit = false;
      for(uint j=0;j<anrl_Wyckoff_set[i].size();j++){
        if(isdigit(anrl_Wyckoff_set[i][j])){
          stringstream ss_tmp;
          if(is_previous_digit){
            ss_tmp << Wyckoff_count << anrl_Wyckoff_set[i][j];
          }
          else {
            ss_tmp << anrl_Wyckoff_set[i][j];
          }
          Wyckoff_count = aurostd::string2utype<uint>(ss_tmp.str());
          is_previous_digit=true;
          continue;
        }
        else {
          for(uint c=0;c<Wyckoff_count;c++){
            stringstream ss_tmp; ss_tmp << anrl_Wyckoff_set[i][j];
            Wyckoff_letters_for_species.push_back(ss_tmp.str());
          }
          Wyckoff_count = 1; //reset
          is_previous_digit=false;
        }
      }
      anrl_grouped_Wyckoff_letters.push_back(Wyckoff_letters_for_species);
    }
    //DEBUG for(uint i=0;i<anrl_grouped_Wyckoff_letters.size();i++){
    //DEBUG  ::print(anrl_grouped_Wyckoff_letters[i]);
    //DEBUG }
    return anrl_grouped_Wyckoff_letters;
  }
}

// ***************************************************************************
// compare::convertWyckoffString2GroupedPositions()
// ***************************************************************************
namespace compare{
  vector<vector<string> > convertWyckoffString2GroupedPositions(
      const string& Wyckoff_letter_string){

    // Groups the Wyckoff positions via species
    // Obtains information from the string of the following form: a,f,g;b,c;a,a
    // (i.e., ";" separates by species, and "," splits the Wyckoff letters for that species)

    vector<vector<string> > grouped_Wyckoff_letters;
    vector<string> Wyckoff_set, token;
    aurostd::string2tokens(Wyckoff_letter_string, token, ";");
    for(uint i=0;i<token.size();i++){ Wyckoff_set.push_back(token[i]); }
    for(uint i=0;i<Wyckoff_set.size();i++){
      aurostd::string2tokens(Wyckoff_set[i], token, ",");
      vector<string> Wyckoff_letters_for_species;
      for(uint j=0;j<token.size();j++){ Wyckoff_letters_for_species.push_back(token[j]); }
      grouped_Wyckoff_letters.push_back(Wyckoff_letters_for_species);
    }
    return grouped_Wyckoff_letters;
  }
}

// ***************************************************************************
// compare::matchableWyckoffPositionSet()
// ***************************************************************************
namespace compare{
  bool matchableWyckoffPositionSet(
      const vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters,
      const vector<vector<string> >& grouped_Wyckoff_letters){

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    bool all_Wyckoffs_matched = false;

    // quick check
    if(grouped_possible_Wyckoff_letters.size() != grouped_Wyckoff_letters.size()){
      //cerr << "sizes do not match" << endl;
      return false;
    }

    //check if each set has the same number of Wyckoff letters for a given species
    vector<int> number_of_letters_1, number_of_letters_2;
    for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
      number_of_letters_1.push_back((int)grouped_possible_Wyckoff_letters[i].size());
      number_of_letters_2.push_back((int)grouped_Wyckoff_letters[i].size());
    }
    aurostd::sort(number_of_letters_1);
    aurostd::sort(number_of_letters_2);

    for(uint i=0;i<number_of_letters_1.size();i++){
      if(number_of_letters_1[i]!=number_of_letters_2[i]){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of Wyckoff letters does not match between structures." << endl;}
        return false;
      }
    }


    vector<std::pair<uint,uint> > matchable_indices;
    uint number_of_matched_sets = 0;
    for(uint m=0;m<grouped_Wyckoff_letters.size();m++){
      bool matched_set = false;
      for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
        if(grouped_possible_Wyckoff_letters[i].size()==grouped_Wyckoff_letters[m].size()){
          uint matched_Wyckoffs = 0;
          vector<uint> index_matched;
          for(uint j=0;j<grouped_possible_Wyckoff_letters[i].size();j++){
            bool Wyckoff_matched = false;
            for(uint k=0;k<grouped_possible_Wyckoff_letters[i][j].size();k++){
              for(uint n=0;n<grouped_Wyckoff_letters[m].size();n++){
                bool index_used = false;
                for(uint t=0;t<index_matched.size();t++){if(n==index_matched[t]){index_used=true;}}
                if(index_used){break;}
                if(grouped_possible_Wyckoff_letters[i][j][k] == grouped_Wyckoff_letters[m][n]){
                  index_matched.push_back(n);
                  Wyckoff_matched = true;
                  matched_Wyckoffs++;
                  break;
                }
              }
              if(Wyckoff_matched){break;}
            }
          }
          if(matched_Wyckoffs==grouped_possible_Wyckoff_letters[i].size() && matched_Wyckoffs==grouped_Wyckoff_letters[m].size()){
            matched_set=true;
          }
        }
        if(matched_set){
          number_of_matched_sets++;
          matched_set=false;
        }
      }
    }
    if(number_of_matched_sets==grouped_possible_Wyckoff_letters.size()){
      //cerr << "possible match!" << endl;
      return true;
    }
    cerr << number_of_matched_sets << " == " << grouped_possible_Wyckoff_letters.size() << endl;

    //cerr << "all_Wyckoffs_matched: " << all_Wyckoffs_matched << endl;
    return all_Wyckoffs_matched;
  }
}

// ***************************************************************************
// compare::matchableSpaceGroups()
// ***************************************************************************
namespace compare{
  bool matchableSpaceGroups(uint space_group_1, uint space_group_2){

    // Determines if two space groups are commensurate.
    // The space groups must be the same or enantiomorphic pairs

    if(space_group_1 == space_group_2){return true;}
    else {
      return matchableEnantiomorphicSpaceGroups(space_group_1, space_group_2);
    }
    return false;
  }
}

// ***************************************************************************
// compare::matchableEnantiomorphicSpaceGroups()
// ***************************************************************************
namespace compare{
  bool matchableEnantiomorphicSpaceGroups(
      uint space_group_1,
      uint space_group_2){

    // Check if the space group has an enantimorphic pair.
    // If it does have an enantiomorphic pair, and it is equal to the
    // second space group, then return true.
    // If it does not, the function below returns the input space group and
    // compares it to the second space group.

    return (SYM::getEnantiomorphSpaceGroupNumber(space_group_1)==space_group_2);
  }
}

// ***************************************************************************
// compare::structuresCompatible()
// ***************************************************************************
namespace compare{
  bool structuresCompatible(
      const structure_container& structure1,
      const structure_container& structure2,
      bool same_species,
      bool ignore_symmetry,
      bool ignore_Wyckoff,
      bool ignore_environment,
      bool ignore_environment_angles, //DX20200320
      bool duplicates_removed){ //DX20190829 - added duplicates_removed

    // Check compatiblity of structures: check species, stoichiometry,
    // number of types, local environments, space group, and Wyckoff
    // positions. Possibly check volumes (i.e. integer multiples, but this
    // is not well tested yet).

    // ---------------------------------------------------------------------------
    // check if species/stoichiometries are compatible
    //DX20190430 - this may take longer, use compound if(same_species==true && matchableSpecies(structures[i].structure_representative,comparison_schemes[j].structure_representative,same_species)==true)
    if(same_species && structure1.compound!=structure2.compound) //DX20190430 - quicker //DX20190702 - changed to "!=" and "false" for speed increase
    { //CO20200106 - patching for auto-indenting
      return false;
    }
    else if(!same_species && structure1.stoichiometry!=structure2.stoichiometry){ //DX20190702 - changed to "!=" and "false" for speed increase
      return false;
    }
    // if already removed duplicate compounds, then structures were already compared, so don't compare again
    else if(!same_species && duplicates_removed && structure1.compound==structure2.compound){
      return false;
    }
    // check if number of atoms are integer multiples of one another (elements only; compounds verified with stoich) //DX20200421
    //DX20210316 [TOO STRICT] if(structure1.ntypes == 1 && structure1.is_structure_generated && structure2.is_structure_generated){
    //DX20210316 [TOO STRICT]   if(structure1.natoms%structure2.natoms!=0 && structure2.natoms%structure1.natoms!=0){
    //DX20210316 [TOO STRICT]     return false;
    //DX20210316 [TOO STRICT]   }
    //DX20210316 [TOO STRICT] }
    // ---------------------------------------------------------------------------
    // check if LFA environments are compatible - DX20190711
    if(!ignore_environment && !compatibleEnvironmentSets(structure1.environments_LFA,structure2.environments_LFA,same_species,ignore_environment_angles,false)){ //DX20200320
      return false;
    }
    // ---------------------------------------------------------------------------
    // check symmetry (if applicable)
    //DX20190702 - checking stoich is redundant for compound checking - if(same_material_stoich==true && structures[i].stoichiometry==comparison_schemes[j].stoichiometry &&
    //DX20190702 [OBSOLETE] if(same_material_stoich==true &&  //DX20190702 - moved stoichiometry up
    if(((ignore_symmetry && ignore_Wyckoff) ||
          (!ignore_symmetry && ignore_Wyckoff &&
           structure1.Pearson == structure2.Pearson &&
           matchableSpaceGroups(structure1.space_group,structure2.space_group)) ||
          (!ignore_symmetry && !ignore_Wyckoff &&
           structure1.Pearson == structure2.Pearson &&
           matchableSpaceGroups(structure1.space_group,structure2.space_group) &&
           matchableWyckoffPositions(structure1.grouped_Wyckoff_positions,structure2.grouped_Wyckoff_positions,same_species)))){
      return true;
    }
    return false;
  }
}

// ***************************************************************************
// XtalFinderCalculator::createStructurePrototypes()
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::groupStructurePrototypes(
    bool same_species,
    bool ignore_symmetry,
    bool ignore_Wyckoff,
    bool ignore_environment,
    bool ignore_environment_angles, //DX20200320
    bool duplicates_removed,
    bool quiet){ //DX20190829 - added duplicates_removed

  // Populates the structure information into the StructurePrototype object.
  // It groups structure based on their species, stoichiometry, number of
  // types, local environments, space group, and Wyckoff positions.
  // A "representative" structure is chosen and will be compared to the
  // possible "duplicates". The misfit values are set to AUROSTD_MAX_DOUBLE
  // until compared.

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  if(!quiet){
    message << "Grouping sets of comparisons.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // variable to store structure sets to compare
  vector<StructurePrototype> comparison_schemes;

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of structures to group: " << structure_containers.size() << endl;}

  // Loop over structures.
  // Group structures that have comparable by stoichiometry and symmetry
  // Optional booleans control certain grouping requirements:
  //   same_species    : groups structures comprised of the same species
  //   ignore_symmetry : ignores space group when grouping (i.e., can group
  //                     structures with different space groups)
  //   ignore_Wyckoff : ignores Wyckoff positions when grouping (i.e., can group
  //                     structures with the same space group number, but different Wyckoff positions)
  //   ignore_environment : ignores LFA environment analysis

  for(uint i=0;i<structure_containers.size(); i++){
    bool scheme_created=false;
    for(uint j=0; j<comparison_schemes.size(); j++){
      //bool same_material_stoich=false;

      if(compare::structuresCompatible(structure_containers[i],
            *comparison_schemes[j].structure_representative,
            same_species,
            ignore_symmetry,
            ignore_Wyckoff,
            ignore_environment,
            ignore_environment_angles,
            duplicates_removed)){ //DX20190829 - added duplicates_removed //DX20200320 - added environment angles

        if(!same_species){
          for(uint e=0;e<structure_containers[i].elements.size();e++){
            if(!aurostd::WithinList(comparison_schemes[j].elements,structure_containers[i].elements[e])){
              comparison_schemes[j].elements.push_back(structure_containers[i].elements[e]);
            }
          }
        }
        addStructure2duplicatesList(comparison_schemes[j], i);
        scheme_created=true;
        break;
      }
    }
    if(!scheme_created){
      StructurePrototype str_proto_tmp;
      setStructureAsRepresentative(str_proto_tmp,i);
      comparison_schemes.push_back(str_proto_tmp);
    }
  }
  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Prepared comparison sets: " << endl;
    cerr << printResults(comparison_schemes, same_species, txt_ft) << endl;
  }

  if(!quiet){
    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  return comparison_schemes;
}

// ***************************************************************************
// XtalFinderCalculator::findDuplicateCompounds()
// ***************************************************************************
void XtalFinderCalculator::findDuplicateCompounds(
    uint num_proc,
    bool remove_duplicates,
    const string& directory,
    const aurostd::xoption& comparison_options){

  // This function finds the duplicate compounds (material type comparison)
  // and writes the results to a file and stores the duplicate count information
  // for the structures.
  // Includes an option to remove the duplicates to avoid biased statistics

  //bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  // ---------------------------------------------------------------------------
  // prepare options 
  bool tmp_same_species = true;
  aurostd::xoption remove_duplicates_options = comparison_options;
  remove_duplicates_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
  remove_duplicates_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",FALSE);
  remove_duplicates_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
  remove_duplicates_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",FALSE);
  remove_duplicates_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",FALSE);

  // ---------------------------------------------------------------------------
  // run multiple comparison function 
  message << "Comparing to remove duplicate materials.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  vector<StructurePrototype> unique_compounds = compareMultipleStructures(
      num_proc,
      tmp_same_species,
      directory,
      remove_duplicates_options);

  message << "Duplicate materials removed.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  // ---------------------------------------------------------------------------
  // include duplicate compounds count in object
  for(uint i=0;i<unique_compounds.size();i++){
    unique_compounds[i].structure_representative->number_compounds_matching_structure = unique_compounds[i].numberOfComparisons();
  }

  // ---------------------------------------------------------------------------
  // write results to files //DX20201229 - consolidated into functions
  if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
    string results_json = printResults(unique_compounds, tmp_same_species, json_ft);
    writeComparisonOutputFile(results_json, directory, json_ft, duplicate_compounds_xf, true);
  }
  else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
    string results_txt = printResults(unique_compounds, tmp_same_species, txt_ft);
    writeComparisonOutputFile(results_txt, directory, txt_ft, duplicate_compounds_xf, true);
  }
  else{
    string results_json = printResults(unique_compounds, tmp_same_species, json_ft);
    writeComparisonOutputFile(results_json, directory, json_ft, duplicate_compounds_xf, true);
    string results_txt = printResults(unique_compounds, tmp_same_species, txt_ft);
    writeComparisonOutputFile(results_txt, directory, txt_ft, duplicate_compounds_xf, true);
  }

  // ---------------------------------------------------------------------------
  // remove duplicates in structure container
  if(remove_duplicates){
    vector<string> duplicates_name;
    for(uint i=0;i<unique_compounds.size();i++){
      for(uint j=0;j<unique_compounds[i].structures_duplicate.size();j++){
        duplicates_name.push_back(unique_compounds[i].structures_duplicate[j]->name);
      }
    }
    for(uint i=0;i<duplicates_name.size();i++){
      removeStructureFromContainerByName(duplicates_name[i]);
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::checkForBetterMatches()
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::checkForBetterMatches(
    vector<StructurePrototype>& prototype_schemes,
    uint num_proc,
    bool check_for_better_matches,
    bool same_species,
    const aurostd::xoption& comparison_options,
    bool quiet){

  // Checks if compounds/structures match better with another prototype group
  // i.e., checks if misfit is smaller when grouped to another prototype

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  if(!quiet || LDEBUG){
    message << "Check if initial structure groupings match better (lower similarity metric) with other groups." << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption check_better_matches_options = comparison_options; //DX20200103
  check_better_matches_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE); // always true for this function
  check_better_matches_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE); //always false for this function
  check_better_matches_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES",TRUE); //always true for this function //DX20200320

  double misfit_min = 0.01; // used for check_for_better_matches : this is quite strict; if too expensive, make more loose
  double misfit_max = misfit_match; // used for check_for_better_matches : otherwise we will compare same family structures which have already been moved (if using !clean_unmatched)

  // ---------------------------------------------------------------------------
  // check which structures could potentially match to others based on misfit
  vector<StructurePrototype> comparison_groups;

  for(uint i=0;i<prototype_schemes.size();i++){
    for(uint j=0;j<prototype_schemes[i].structures_duplicate.size();j++){
      if((check_for_better_matches && prototype_schemes[i].mapping_info_duplicate[j].misfit > misfit_min &&
            prototype_schemes[i].mapping_info_duplicate[j].misfit < misfit_max) || // find better match
          (!check_for_better_matches && (prototype_schemes[i].mapping_info_duplicate[j].misfit > misfit_max || aurostd::isequal(prototype_schemes[i].mapping_info_duplicate[j].misfit,1.0,1e-6) || aurostd::isequal(prototype_schemes[i].mapping_info_duplicate[j].misfit,AUROSTD_MAX_DOUBLE,1e-6)))){   // find a match //DX20191218
        StructurePrototype str_proto_tmp;
        bool found_new_match=false;
        // ---------------------------------------------------------------------------
        // check for other compatible representative structures
        // start_index=i+1 : (only need to search forward for better matches, due to appendStructurePrototypes() scheme)
        uint start_index = 0;
        if(check_for_better_matches){ start_index = i+1; }
        for(uint k=start_index;k<prototype_schemes.size();k++){
          if(i!=k){ // don't perform the same comparison again //DX20200414
            if(compare::structuresCompatible(*prototype_schemes[i].structure_representative, //compares StructurePrototype objects
                  *prototype_schemes[k].structure_representative,
                  same_species,
                  check_better_matches_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"),
                  check_better_matches_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"),
                  check_better_matches_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"),
                  check_better_matches_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320
                  false)){ // can check based on representatives; duplicate info matches its representative info //DX20200103 - condensed booleans to xoptions
              if(!quiet || LDEBUG){
                message << "Found potential match for " << prototype_schemes[i].structures_duplicate[j]->name << ": " << prototype_schemes[k].structure_representative->name;
                pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
              }

              // ---------------------------------------------------------------------------
              // reverse the convention, single structurePrototype object to find best match
              // i.e., duplicate -> representative and representatives -> duplicates
              if(!found_new_match){
                str_proto_tmp.copyPrototypeInformation(prototype_schemes[i]);
                setStructureAsRepresentative(str_proto_tmp, prototype_schemes[i].structures_duplicate[j]);
                // store the current match so we can check fast if it matches to any other
                addStructure2duplicatesList(str_proto_tmp, prototype_schemes[i].structure_representative); // store the current match structure
                found_new_match=true;
              }
              if(k!=j){
                addStructure2duplicatesList(str_proto_tmp,prototype_schemes[k].structure_representative); // store the current match structure
              }
            }
          }
        }
        if(found_new_match){
          comparison_groups.push_back(str_proto_tmp);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // if no alternative matches found, return immediately //DX20210111
  if(comparison_groups.size() == 0){ return prototype_schemes; }

  // ---------------------------------------------------------------------------
  // compare structures
  vector<StructurePrototype> other_matches_schemes = runComparisonScheme(
      comparison_groups,
      same_species,
      num_proc,
      check_better_matches_options,
      quiet); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // check if there are any better matches and reorganize if necessary
  // the original match is stored in the first position
  for(uint i=0;i<other_matches_schemes.size();i++){
    double min_misfit = aurostd::abs(other_matches_schemes[i].mapping_info_duplicate[0].misfit); // put first as min //abs to turn -1 into 1 for comparison //DX20191218
    uint min_index = 0;
    for(uint j=1;j<other_matches_schemes[i].mapping_info_duplicate.size();j++){ //DX20191218
      if(other_matches_schemes[i].mapping_info_duplicate[j].misfit<min_misfit && aurostd::isdifferent(other_matches_schemes[i].mapping_info_duplicate[j].misfit,AUROSTD_MAX_DOUBLE,1e-6)){
        min_misfit=other_matches_schemes[i].mapping_info_duplicate[j].misfit;
        min_index=j;
      }
    }
    // ---------------------------------------------------------------------------
    // move if found better match
    if(min_index!=0){
      for(uint j=0;j<prototype_schemes.size();j++){
        // add structure to its better matching representative
        if(prototype_schemes[j].structure_representative->name == other_matches_schemes[i].structures_duplicate[min_index]->name){
          if(!quiet || LDEBUG){
            message << other_matches_schemes[i].structure_representative->name << " matches better with " << prototype_schemes[j].structure_representative->name;
            pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
          }
          addStructure2duplicatesList(prototype_schemes[j],other_matches_schemes[i].structure_representative); // store the current match structure
          prototype_schemes[j].mapping_info_duplicate.back()=other_matches_schemes[i].mapping_info_duplicate[min_index]; //DX20191218
        }
        // remove from old representative
        //DX20201223 [OBSOLETE] if(check_better_matches_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED") && prototype_schemes[j].structure_representative->name == other_matches_schemes[i].structures_duplicate[0]->name)
        if(prototype_schemes[j].structure_representative->name == other_matches_schemes[i].structures_duplicate[0]->name){ //DX20201223 - always remove the old match
          for(uint k=0;k<prototype_schemes[j].structures_duplicate.size();k++){
            if(prototype_schemes[j].structures_duplicate[k]->name == other_matches_schemes[i].structure_representative->name){
              if(!quiet || LDEBUG){
                message << "removing " << other_matches_schemes[i].structure_representative->name << " from " << prototype_schemes[j].structure_representative->name << " set";
                pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
              }
              prototype_schemes[j].removeNonDuplicate(k);
              break;
            }
          }
        }
      }
    }
    else{
      if(!quiet){
        message << other_matches_schemes[i].structure_representative->name << " matches better with original set " << other_matches_schemes[i].structures_duplicate[0]->name;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
    }
  }

  if(LDEBUG){
    for(uint i=0;i<prototype_schemes.size();i++){
      cerr << __AFLOW_FUNC__ << " prototype_schemes[i]: " << prototype_schemes[i] << endl;
    }
  }

  return prototype_schemes;
}

// ***************************************************************************
// XtalFinderCalculator::representativePrototypeForICSDRuns()
// ***************************************************************************
void XtalFinderCalculator::representativePrototypeForICSDRunsNEW(
    vector<StructurePrototype>& comparison_schemes){

  // For ICSD comparisons; make the structure with the smallest ICSD number the
  // representative structure (normally the oldest).

  for(uint i=0;i<comparison_schemes.size();i++){
    if(comparison_schemes[i].structures_duplicate.size()>0){
      vector<string> ICSD_entries;
      ICSD_entries.push_back(compare::findICSDName(comparison_schemes[i].structure_representative->name));
      for(uint j=0;j<comparison_schemes[i].structures_duplicate.size();j++){
        ICSD_entries.push_back(compare::findICSDName(comparison_schemes[i].structures_duplicate[j]->name));
      }
      string min_ICSD_entry = compare::findMinimumICSDEntry(ICSD_entries);
      if((comparison_schemes[i].structure_representative->name.find(min_ICSD_entry) == std::string::npos) && !min_ICSD_entry.empty()){ //DX20191108 - add not empty case //DX20200709 - aurostd::substring2bool to find
        for(uint j=0;j<comparison_schemes[i].structures_duplicate.size();j++){
          if(comparison_schemes[i].structures_duplicate[j]->name.find(min_ICSD_entry) != std::string::npos){ //DX20200709 - aurostd::substring2bool to find

            structure_container *str_container_tmp;
            str_container_tmp = comparison_schemes[i].structure_representative;
            setStructureAsRepresentative(comparison_schemes[i],comparison_schemes[i].structures_duplicate[j]);
            comparison_schemes[i].structures_duplicate[j] = str_container_tmp;
            break;
          }
        }
      }
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::runComparisonThreads()
// ***************************************************************************
void XtalFinderCalculator::runComparisonThreads(
    uint index,
    vector<StructurePrototype>& comparison_schemes,
    const vector<std::pair<uint,uint> >& vstart_indices,
    const vector<std::pair<uint,uint> >& vend_indices,
    bool same_species,
    bool scale_volume,
    bool optimize_match){

  // Run comparisons on a paricular thread
  // If the xstructure is not generated, it will generate a local copy for
  // the single comparison only (prevents threads from accessing on the same
  // memory block)

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  const std::pair<uint,uint>& start_indices = vstart_indices[index];
  const std::pair<uint,uint>& end_indices = vend_indices[index];
  uint i_min=start_indices.first; uint i_max=end_indices.first;
  uint j_min=0; uint j_max=0;

  // to prevent nested multi-processes
  uint num_proc_orig = num_proc;
  num_proc=1;

  for(uint i=i_min; i<=i_max; i++){
    // ---------------------------------------------------------------------------
    // copy representative structure
    structure_container structure_rep_tmp = *comparison_schemes[i].structure_representative;

    // ---------------------------------------------------------------------------
    // to loop properly
    if(i==i_min){
      j_min=start_indices.second;
      if(i==i_max){j_max=end_indices.second;}
      else {j_max=comparison_schemes[i].structures_duplicate.size();} //-1 since in loop: j<=j_max
    }
    else if(i==i_max){j_min=0; j_max=end_indices.second;}
    else {j_min=0; j_max=comparison_schemes[i].structures_duplicate.size();} //-1 since in loop: j<=j_max

    // ---------------------------------------------------------------------------
    // begin comparison loop
    for(uint j=j_min; j<j_max; j++){
      if(!structure_rep_tmp.is_structure_generated){
        if(!compare::generateStructure(structure_rep_tmp.name,structure_rep_tmp.source,structure_rep_tmp.relaxation_step,structure_rep_tmp.structure,*p_oss)){ //DX20200429
          message << "Could not generate representative structure (" << structure_rep_tmp.name << ").";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
        }
      }
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Loaded representative structure = " << structure_rep_tmp.name << endl; }

      // ---------------------------------------------------------------------------
      // copy duplicate structure
      structure_container structure_dup_tmp = *comparison_schemes[i].structures_duplicate[j];
      if(!structure_dup_tmp.is_structure_generated){
        if(!compare::generateStructure(structure_dup_tmp.name,structure_dup_tmp.source,structure_dup_tmp.relaxation_step,structure_dup_tmp.structure,*p_oss)){ //DX20200429
          message << "Could not generate duplicate structure (" << structure_dup_tmp.name << ").";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
        }
      }
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Loaded duplicate structure = " << structure_dup_tmp.name << endl; }

      // ---------------------------------------------------------------------------
      // call the main comparison function
      structure_mapping_info final_misfit_info = compare::initialize_misfit_struct(); //DX20191218
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Comparing " << structure_rep_tmp.name << " and " << structure_dup_tmp.name <<  endl; }
      compareStructures(structure_rep_tmp,
          structure_dup_tmp,
          final_misfit_info,
          same_species,
          scale_volume,
          optimize_match); //DX20200103 - condensed booleans to xoptions

      // ---------------------------------------------------------------------------
      // store the figure of misfit
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Finished comparing " << structure_rep_tmp.name << " and " << structure_dup_tmp.name << endl;
        cerr << __AFLOW_FUNC__ << " Comparison complete, misfit = " << final_misfit_info.misfit << "." << endl;
      }
      comparison_schemes[i].mapping_info_duplicate[j]=final_misfit_info; //DX20191218
    }
  }

  // reset number of processes
  num_proc=num_proc_orig;
}

// ***************************************************************************
// XtalFinderCalculator::runComparisons [NON-THREADED VERSION]
// ***************************************************************************
void XtalFinderCalculator::runComparisons(
    vector<StructurePrototype>& comparison_schemes,
    bool same_species,
    bool scale_volume,
    bool optimize_match){

  // Run comparisons [NON-THREADED VERSION]
  // If the xstructure is not generated, it will generate and store in
  // the structure_container object (so we only generate once).
  // This is only possible for a non-threaded process, otherwise we may
  // run into thread overwriting problems.

  stringstream message;

  for(uint i=0;i<comparison_schemes.size(); i++){
    for(uint j=0;j<comparison_schemes[i].structures_duplicate.size();j++){
      if(j==0){
        // ---------------------------------------------------------------------------
        // generate the representative structure if necessary
        if(!comparison_schemes[i].structure_representative->is_structure_generated){
          if(!compare::generateStructure(comparison_schemes[i].structure_representative->name,comparison_schemes[i].structure_representative->source,comparison_schemes[i].structure_representative->relaxation_step,comparison_schemes[i].structure_representative->structure,*p_oss)){
            message << "Could not generate representative structure (" << comparison_schemes[i].structure_representative->name << ").";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
          }
        }
      }
      structure_mapping_info final_misfit_info = compare::initialize_misfit_struct();
      // ---------------------------------------------------------------------------
      // generate the duplicate structure if necessary
      if(!comparison_schemes[i].structures_duplicate[j]->is_structure_generated){
        if(!compare::generateStructure(comparison_schemes[i].structures_duplicate[j]->name,comparison_schemes[i].structures_duplicate[j]->source,comparison_schemes[i].structures_duplicate[j]->relaxation_step,comparison_schemes[i].structures_duplicate[j]->structure,*p_oss)){ //DX20200429
          message << "Could not generate duplicate structure (" << comparison_schemes[i].structures_duplicate[j]->name << ").";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ERROR_);
        }
      }

      // ---------------------------------------------------------------------------
      // call the main comparison function
      compareStructures(*comparison_schemes[i].structure_representative,
          *comparison_schemes[i].structures_duplicate[j],
          final_misfit_info,
          same_species,
          scale_volume,
          optimize_match); //DX20200103 - condensed booleans to xoptions

      // ---------------------------------------------------------------------------
      // store the figure of misfit
      comparison_schemes[i].mapping_info_duplicate[j]=final_misfit_info;
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::runComparisonScheme():
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::runComparisonScheme(
    vector<StructurePrototype>& comparison_schemes,
    bool same_species,
    uint num_proc,
    const aurostd::xoption& comparison_options,
    bool quiet){

  // Runs all comparisons automatically.
  // Compares and regroups similar structures until all structures are
  // matched or the comparisons are exhausted.

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  // create new object for comparisons
  vector<StructurePrototype> prototypes_final;

  if(!quiet){
    message << "Running comparisons ...";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // print initial grouped sets of comparisons
  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Number of comparison sets: " << comparison_schemes.size() << endl;
    cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cerr << printResults(comparison_schemes, same_species, txt_ft) << endl;
  }

  uint number_of_comparisons = 0;
  vector<std::pair<uint,uint> > start_indices, end_indices;

  bool scale_volume = comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME");
  bool optimize_match = comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH");
#ifdef AFLOW_MULTITHREADS_ENABLE
  uint num_comparison_threads = 1;

  // ---------------------------------------------------------------------------
  // THREADED VERSION - START
  number_of_comparisons = 0;
  for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
  num_comparison_threads = aurostd::min(num_proc,number_of_comparisons);
  splitComparisonIntoThreads(comparison_schemes, num_comparison_threads, start_indices, end_indices);

  xthread::xThread xt(num_comparison_threads);
  std::function<void(uint, vector<StructurePrototype>&,
      const vector<std::pair<uint, uint> >&, const vector<std::pair<uint, uint> >&,
      bool, bool, bool)> fn = std::bind(&XtalFinderCalculator::runComparisonThreads, this, _1, _2, _3, _4, _5, _6, _7);
  xt.run(num_comparison_threads, fn,
      comparison_schemes,
      start_indices,
      end_indices,
      same_species,
      scale_volume,
      optimize_match);

#else

  // ---------------------------------------------------------------------------
  // NON-THREADED VERISON - START
  if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Non-threaded version." << endl; }
  runComparisons(comparison_schemes,
      same_species,
      scale_volume,
      optimize_match);

#endif

  // ---------------------------------------------------------------------------
  // count the number of mismatches (i.e. mis > 0.1)
  uint num_mismatches_orig=numberOfMismatches(comparison_schemes);
  uint num_mismatches=num_mismatches_orig;

  // ---------------------------------------------------------------------------
  //DX20190504 - added clean unmatched option - START
  if(!comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED") &&
      comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND")){
    return comparison_schemes;
  }
  //DX20190504 - added clean unmatched option - END

  // ---------------------------------------------------------------------------
  // if print only matches to the input structure //DX20201230
  if(comparison_options.flag("COMPARISON_OPTIONS::PRINT_MATCHES_TO_INPUT_ONLY") && comparison_schemes.size()>0){
    appendStructurePrototypes(comparison_schemes,
        prototypes_final,
        comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED"),
        quiet);
    vector<StructurePrototype> input_prototype_only;
    input_prototype_only.push_back(prototypes_final[0]);
    return input_prototype_only;
  }

  if(num_mismatches > 0 && !comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND") && !quiet){
    message << "Number of unmatched structures: " << num_mismatches << ". Continuing comparisons ...";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // For this first iteration: LFA may have incorrectly grouped, so we need to check if the structures belong to other groups
  // OR after removing duplicate compounds these compounds remain separate, so for !same_species comparisons we need to check
  // if they should match with other groups
  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS") ||
      comparison_options.flag("COMPARISON_OPTIONS::CHECK_OTHER_GROUPINGS")){
    aurostd::xoption check_better_matches_options = comparison_options;
    check_better_matches_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",FALSE); //DX20200320 - changed from true to false
    comparison_schemes = checkForBetterMatches(comparison_schemes, num_proc,
        false,
        same_species,
        check_better_matches_options,
        quiet); //DX20200103 - condensed booleans to xoptions
  }

  // ---------------------------------------------------------------------------
  // regroup comparisons based on misfit value
  if(num_mismatches==0 && !comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND")){
    appendStructurePrototypes(comparison_schemes,
        prototypes_final,
        comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED"),
        quiet); //DX20200103
  }

  // ---------------------------------------------------------------------------
  // Loop: continue comparison until all strucutures are matched or all comparisons schemes exhaused
  while(num_mismatches!=0){
    // regroup comparisons based on misfit value
    appendStructurePrototypes(
        comparison_schemes,
        prototypes_final,
        comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED"),
        quiet); //DX20200103

    // return if only one round of comparison is requested
    if(comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND")){return prototypes_final;}

    // reorder structures so minimum ICSD is the representative structure
    if(comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON")){ representativePrototypeForICSDRunsNEW(comparison_schemes); }

    // split into threads
    number_of_comparisons=0;
    for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " number_of_comparisons: " << number_of_comparisons << endl; }

    if(number_of_comparisons>0){
      if(!quiet){
        message << "Continuing comparisons to match " << num_mismatches << " structures ...";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      }
      if(LDEBUG) { cerr << __AFLOW_FUNC__ << ": Number of comparisons is not zero... " << number_of_comparisons << endl; }
#ifdef AFLOW_MULTITHREADS_ENABLE
      num_comparison_threads = aurostd::min(num_proc,number_of_comparisons);

      // THREADED VERISON - START
      // split into threads
      start_indices.clear(); end_indices.clear();
      splitComparisonIntoThreads(comparison_schemes, num_comparison_threads, start_indices, end_indices);
      xt.run(num_comparison_threads, fn,
          comparison_schemes,
          start_indices,
          end_indices,
          same_species,
          scale_volume,
          optimize_match);
      // THREADED VERISON - END
#else
      //SINGLE THREAD - START
      start_indices.clear(); end_indices.clear();
      uint single_thread=1;
      splitComparisonIntoThreads(comparison_schemes, single_thread, start_indices, end_indices);
      runComparisonThreads(0, comparison_schemes,
          start_indices,
          end_indices,
          same_species,
          scale_volume,
          optimize_match);
      //SINGLE THREAD - END
#endif
    }

    // update number of mismatches
    num_mismatches_orig=num_mismatches;
    num_mismatches=numberOfMismatches(comparison_schemes);

    if(num_mismatches > 0 && !quiet){
      message << "Number of unmatched structures: " << num_mismatches << ". Continuing comparisons ...";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ensure while loop is controlled
    if(num_mismatches>=num_mismatches_orig){
      *p_oss << "The number of mismatches is increasing ... impossible (bug in comparision framework). "
        << "Email aflow@groups.io." << endl;
      prototypes_final.clear();
      return prototypes_final;
    }
  }
  // end of while loop

  // append new prototype groupings
  appendStructurePrototypes(
      comparison_schemes,
      prototypes_final,
      comparison_options.flag("COMPARE_STRUCTURE::CLEAN_UNMATCHED"),
      quiet); //DX20200103

  if(!quiet){
    message << "Comparisons complete!";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
  }

  return prototypes_final;
}

// ***************************************************************************
// compare::calculateDivisors()
// ***************************************************************************
namespace compare{
  vector<std::pair<uint,uint> > calculateDivisors(uint number){

    // Determine the divisors for a number
    // Used to check if the number of permutations follows number theory

    vector<std::pair<uint,uint> > divisor_pairs;
    std::pair<uint,uint> tmp_pair;
    for(uint i=1; i<=number; i++){
      if(number%i==0){
        tmp_pair.first = i; tmp_pair.second = number/i;
        divisor_pairs.push_back(tmp_pair);
      }
    }
    return divisor_pairs;
  }
}

// ***************************************************************************
// compare::checkNumberOfGroupings()
// ***************************************************************************
namespace compare{
  bool checkNumberOfGroupings(
      const vector<StructurePrototype>& comparison_schemes,
      uint number){

    // Check if the number of comparisons for the atom decoration comparisons
    // are commensurate. The number of unqiue atom decorations must be
    // divisors of the total number of permutations,
    // e.g., 6 permutations:
    //   1 unique in a group of 6;
    //   2 unique in groups of 3;
    //   3 unique in groups of 2; or
    //   6 unique in groups of 1.
    // If this is not the case, then there is something wrong with the
    // comparison result (the atom decoration may match better with another
    // decoration; see checkForBetterMatches()).

    vector<std::pair<uint,uint> > divisor_pairs = calculateDivisors(number);
    uint divisor_index = 0;
    for(uint i=0;i<divisor_pairs.size();i++){
      if(comparison_schemes.size() == divisor_pairs[i].first){
        divisor_index = i;
        break;
      }
    }
    uint num_consistent = 0;
    for(uint j=0; j<comparison_schemes.size(); j++){
      if(comparison_schemes[j].structures_duplicate.size()+1 == divisor_pairs[divisor_index].second){ //+1 to include representative
        num_consistent++;
      }
    }
    if(num_consistent != comparison_schemes.size()){
      return false;
    }
    return true;
  }
}

// ***************************************************************************
// XtalFinderCalculator::makeRepresentativeEvenPermutation()
// ***************************************************************************
void XtalFinderCalculator::makeRepresentativeEvenPermutation(
    vector<StructurePrototype>& comparison_schemes,
    const vector<string>& name_order){

  // Make sure the even permutation is the representative.
  // If there are multiple even permutations in a given set of comparisons,
  // default to the mininum even permutation. (DX, may want to change default)

  uint representative_permutation_num = 0, min_even_duplicate_permutation_num = AUROSTD_MAX_UINT;
  uint duplicate_permutation_num = 0, min_duplicate_index = 0;

  for(uint i=0; i<comparison_schemes.size(); i++){	
    //Find representative permutation number
    representative_permutation_num = 0;
    for(uint j=0; j<name_order.size(); j++){			
      if(comparison_schemes[i].structure_representative->name == name_order[j]){
        representative_permutation_num = j;
        break;
      }
    }
    // ---------------------------------------------------------------------------
    // find duplicate permutation number
    min_even_duplicate_permutation_num = AUROSTD_MAX_UINT;
    duplicate_permutation_num = 0;
    min_duplicate_index = 0;
    for(uint p=0; p<comparison_schemes[i].structures_duplicate.size(); p++){
      for(uint j=0; j<name_order.size(); j++){			
        if(comparison_schemes[i].structures_duplicate[p]->name == name_order[j]){
          duplicate_permutation_num = j;
          // Check if proto is minimum/even permutation
          if(duplicate_permutation_num < min_even_duplicate_permutation_num && duplicate_permutation_num%2 == 0){
            min_duplicate_index = p;
            min_even_duplicate_permutation_num = duplicate_permutation_num;
          }
        }
      }
    }
    // ---------------------------------------------------------------------------
    // if representative permutation is already even, only swap representative
    // and duplicate proto if proto permutation is less and even
    // else replace automatically if proto permutation is not AUROSTD_MAX_DOUBLE
    if((representative_permutation_num%2 == 0 && min_even_duplicate_permutation_num < representative_permutation_num) ||
        (representative_permutation_num%2 != 0 && min_even_duplicate_permutation_num != AUROSTD_MAX_UINT)){
      structure_container *str_container_tmp;
      str_container_tmp = comparison_schemes[i].structure_representative;
      setStructureAsRepresentative(comparison_schemes[i],comparison_schemes[i].structures_duplicate[min_duplicate_index]);
      comparison_schemes[i].structures_duplicate[min_duplicate_index] = str_container_tmp;
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::numberOfMismatches()
// ***************************************************************************
uint XtalFinderCalculator::numberOfMismatches(
    const vector<StructurePrototype>& comparison_schemes){

  // Count the number of comparisons that have a misfit greater than
  // misfit_match (default: 0.1) i.e., not a match

  int num_mismatches=0;
  for(uint i=0; i<comparison_schemes.size(); i++){
    for(uint j=0; j<comparison_schemes[i].mapping_info_duplicate.size(); j++){
      if(comparison_schemes[i].mapping_info_duplicate[j].misfit > misfit_match || aurostd::isequal(comparison_schemes[i].mapping_info_duplicate[j].misfit,AUROSTD_MAX_DOUBLE,1e-6)){
        num_mismatches+=1;
      }
    }
  }
  return num_mismatches;
}

// ***************************************************************************
// XtalFinderCalculator::appendStructurePrototypes()
// ***************************************************************************
void XtalFinderCalculator::appendStructurePrototypes(
    vector<StructurePrototype>& comparison_schemes,
    vector<StructurePrototype>& prototypes_final,
    bool clean_unmatched, //DX20190506
    bool quiet){

  // This "cleans" the StrucuturePrototype objects by removing all the
  // structures that do not match with the representative structure.
  // Then, it takes the mismatched structures and puts them into new
  // StructurePrototype objects to be compared.

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  uint DUPLICATE_NAME_WIDTH = 100, MISFIT_WIDTH=15;
  string dash_line_separator = std::string(DUPLICATE_NAME_WIDTH+MISFIT_WIDTH, '-');

  if(LDEBUG){
    cerr << printResults(comparison_schemes, true, txt_ft) << endl;
  }

  uint number_of_comparisons = comparison_schemes.size();

  vector<StructurePrototype> tmp_list;
  for(uint i=0; i<number_of_comparisons; i++){
    bool first_mismatch=true;
    for(uint j=0; j<comparison_schemes[i].mapping_info_duplicate.size(); j++){
      if(comparison_schemes[i].mapping_info_duplicate[j].misfit > misfit_match){
        // First, store any family prototype information
        if(comparison_schemes[i].mapping_info_duplicate[j].misfit <= misfit_family){
          addStructure2sameFamilyList(comparison_schemes[i],j); //DX20190814 - consolidated below into single function
        }
        // Take first mismatch and make as the representative structure in the new object
        if(first_mismatch){
          StructurePrototype str_proto_tmp;
          str_proto_tmp.copyPrototypeInformation(comparison_schemes[i]);
          setStructureAsRepresentative(str_proto_tmp, comparison_schemes[i].structures_duplicate[j]);
          tmp_list.push_back(str_proto_tmp);
          if(clean_unmatched){ comparison_schemes[i].removeNonDuplicate(j); j--; } //DX20190504 - put in if-statement
          first_mismatch=false;
        }
        // If not the first mismatch, add as a proto structure in the new object
        else if(!first_mismatch){
          tmp_list.back().copyDuplicate(comparison_schemes[i],j);
          if(clean_unmatched){ comparison_schemes[i].removeNonDuplicate(j); j--; } //DX20190504 - put in if-statement
        }
      }
    }

    // if not quiet, print the comparison results to the screen
    // (useful for long comparison times or if the program terminates early)
    if(!quiet){
      message << "Identified unique prototype: " << endl;
      message << "   prototype=" << comparison_schemes[i].structure_representative->name << endl;
      if(comparison_schemes[i].structures_duplicate.size()==0){
        message << "   No duplicates. " << endl;
      }
      else {
        message << "   " << setw(DUPLICATE_NAME_WIDTH) << std::left << "List of duplicates"
          << setw(MISFIT_WIDTH) << std::left << "misfit value" << endl;
        message << "   " << setw(DUPLICATE_NAME_WIDTH) << std::left
          << dash_line_separator << endl;
        for(uint d=0;d<comparison_schemes[i].structures_duplicate.size();d++){
          message << "   " << setw(DUPLICATE_NAME_WIDTH) << std::left << comparison_schemes[i].structures_duplicate[d]->name
            << setw(MISFIT_WIDTH) << std::left << comparison_schemes[i].mapping_info_duplicate[d].misfit << endl;
        }
      }
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_); //DX+CO20201119
    }

    // Store finished (already compared) schemes in prototypes_final
    prototypes_final.push_back(comparison_schemes[i]);
  }

  // Store newly generated schemes (not compared yet) into comparison_schemes
  comparison_schemes=tmp_list;

  if(LDEBUG){
    cerr << printResults(comparison_schemes, true, txt_ft) << endl;
  }
}

// ***************************************************************************
// XtalFinderCalculator::combinePrototypesOfDifferentSymmetry()
// ***************************************************************************
void XtalFinderCalculator::combinePrototypesOfDifferentSymmetry(
    vector<StructurePrototype>& prototypes_final,
    bool same_species,
    uint num_proc) {

  // Checks to see if prototypes of different space groups are similar.
  // If they are, combine the StructurePrototype objects into one.
  // When combining, we keep the "representative" prototype as the one with a higher
  // symmetry (i.e. higher space group).
  // This is an optional function; we may not want to do this when
  // comparing prototypes or comparing material properties

  uint nprototypes=prototypes_final.size();

  for(uint i=0;i<nprototypes;i++){
    int min_index=-1;
    double min_misfit=AUROSTD_MAX_DOUBLE;
    structure_mapping_info min_misfit_info = compare::initialize_misfit_struct(); //DX20191218
    for(uint j=i;j<prototypes_final.size();j++){
      if(
          // If same_species==true
          (same_species &&
           compare::matchableSpecies(prototypes_final[i].structure_representative->structure,prototypes_final[j].structure_representative->structure,same_species) &&
           prototypes_final[i].stoichiometry==prototypes_final[j].stoichiometry &&
           prototypes_final[i].Pearson==prototypes_final[j].Pearson &&
           !compare::matchableSpaceGroups(prototypes_final[i].space_group,prototypes_final[j].space_group)) ||
          // If same_species==false
          (!same_species &&
           prototypes_final[i].stoichiometry==prototypes_final[j].stoichiometry &&
           prototypes_final[i].Pearson==prototypes_final[j].Pearson &&
           !compare::matchableSpaceGroups(prototypes_final[i].space_group,prototypes_final[j].space_group))
        ){
        double final_misfit=AUROSTD_MAX_DOUBLE;
        structure_mapping_info final_misfit_info = compare::initialize_misfit_struct(); //DX20191218
        bool scale_volume=true; //default is true
        bool optimize_match=false; //default is false
        compare::aflowCompareStructure(
            prototypes_final[i].structure_representative->structure,
            prototypes_final[j].structure_representative->structure,
            same_species,
            scale_volume,
            optimize_match,
            final_misfit,
            final_misfit_info,
            num_proc); //DX20191122 - move ostream to end  //DX20191218 - added misfit_info
        if(final_misfit < min_misfit){
          min_misfit_info=final_misfit_info; //DX20191218
          min_misfit=final_misfit;
          min_index=j;
        }
      }
    }
    // If one prototype is similar to another, add to one with higher space group
    if(min_misfit!=AUROSTD_MAX_DOUBLE){
      int sg_ind=-1;
      int other_ind=-1;
      if(prototypes_final[i].space_group > prototypes_final[min_index].space_group){
        sg_ind=i;
        other_ind=min_index;
      }
      else {
        sg_ind=min_index;
        other_ind=i;
      }
      // Transfer info to prototype with higher space group
      addStructure2duplicatesList(prototypes_final[sg_ind], prototypes_final[other_ind].structure_representative);
      prototypes_final[sg_ind].mapping_info_duplicate.push_back(min_misfit_info);
      for(uint j=0;j<prototypes_final[other_ind].structures_duplicate.size();j++){
        addStructure2duplicatesList(prototypes_final[sg_ind], prototypes_final[other_ind].structures_duplicate[j]);
        prototypes_final[sg_ind].mapping_info_duplicate.push_back(prototypes_final[other_ind].mapping_info_duplicate[j]); //this is the approximate misfit to the new representative ...
      }
      // Delete the prototype with the lower space group
      prototypes_final.erase(prototypes_final.begin()+other_ind);
      // If the index deleted was less than the initial loop (i), then need to reduce iterator
      if(other_ind<=(int)i){
        i--;
      }
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::printResults()
// ***************************************************************************
string XtalFinderCalculator::printResults(
    const vector<StructurePrototype>& prototypes_final,
    bool same_species,
    filetype format){

  // Print the comparison results in a JSON (json_ft) or TXT (txt_ft) format
  // In general, the structure along with the misfit value is printed
  // If material properties are provided, then the properties will be displayed
  // next to the misfit value

  stringstream ss_out;
  bool roff=true; //round off

  // ---------------------------------------------------------------------------
  // json format
  if(format==json_ft){
    ss_out << "[" << endl;
    for(uint j=0; j<prototypes_final.size(); j++){
      ss_out << prototypes_final[j];
      if(j!=prototypes_final.size()-1){
        ss_out << "," << endl;
      }
    }
    ss_out << endl << "]" << endl;
  }

  // ---------------------------------------------------------------------------
  // text format
  else if(format==txt_ft){
    int indent_spacing = 2;
    int structure_spacing = 100; // structure name spacing
    int misfit_spacing = 15;
    int property_spacing = 35;
    int num_properties = 0;
    for(uint j=0; j<prototypes_final.size(); j++){ num_properties = aurostd::max((int)prototypes_final[j].property_names.size(),num_properties); }

    // ---------------------------------------------------------------------------
    // set length of separators based on table width
    string equal_line_separator = std::string(indent_spacing+structure_spacing+misfit_spacing+(num_properties*property_spacing), '=');
    string dash_line_separator = std::string(indent_spacing+structure_spacing+misfit_spacing+(num_properties*property_spacing), '-');

    for(uint j=0; j<prototypes_final.size(); j++){

      // ---------------------------------------------------------------------------
      // start of new prototype
      ss_out << equal_line_separator << endl;
      ss_out << "# ";

      // ---------------------------------------------------------------------------
      // print compound or stoichometry of prototype
      if(same_species){ ss_out << prototypes_final[j].structure_representative->compound; }
      else { ss_out << aurostd::joinWDelimiter(prototypes_final[j].stoichiometry,":"); }

      // ---------------------------------------------------------------------------
      // print space group and Wyckoff positions of prototype
      ss_out << "  SG=#" << prototypes_final[j].space_group;
      ss_out << "  Wyckoffs=" << compare::printWyckoffString(prototypes_final[j].grouped_Wyckoff_positions,true); //DX20190228 - remove count

      // ---------------------------------------------------------------------------
      // print number of duplicate compounds/stuctures (material-/structure-type)
      if(same_species){
        uint number_of_duplicates = numberOfDuplicates(prototypes_final[j]); //DX20190506 - made function
        ss_out << "  compounds_duplicate=" << number_of_duplicates << endl; //DX20190228 - add count
      }
      else {
        uint number_of_duplicates = numberOfDuplicates(prototypes_final[j]); //DX20190506 - made function
        ss_out << "  structures_duplicate=" << number_of_duplicates;
        uint number_duplicate_compounds = 0;
        for(uint k=0;k<prototypes_final[j].structures_duplicate.size();k++){
          number_duplicate_compounds+=prototypes_final[j].structures_duplicate[k]->number_compounds_matching_structure;
        }
        number_duplicate_compounds+= number_of_duplicates+prototypes_final[j].structure_representative->number_compounds_matching_structure; //DX20190321 - need to update variable, otherwise may not enter if statement
        if(number_duplicate_compounds!=0){
          ss_out << "  duplicate_compounds=" << number_duplicate_compounds; //DX20190228 - add count
        }
        ss_out << endl;
      }

      // ---------------------------------------------------------------------------
      // print prototype designation info
      if(prototypes_final[j].aflow_label.size()!=0){
        ss_out << "  aflow_label=" << prototypes_final[j].aflow_label << endl;
        ss_out << "  aflow_parameter_list=" << aurostd::joinWDelimiter(prototypes_final[j].aflow_parameter_list,",") << endl;
        ss_out << "  aflow_parameter_values=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(prototypes_final[j].aflow_parameter_values,8,roff),",") << endl;
      }

      // ---------------------------------------------------------------------------
      // print matching AFLOW prototype labels
      if(prototypes_final[j].matching_aflow_prototypes.size()!=0){
        ss_out << "  matching_aflow_prototypes=" << aurostd::joinWDelimiter(prototypes_final[j].matching_aflow_prototypes,",") << endl;
      }

      // ---------------------------------------------------------------------------
      // print unique atom decorations
      // perhaps add which permutations are duplicates
      if(prototypes_final[j].atom_decorations_equivalent.size()!=0){
        vector<string> unique_decorations;
        for(uint d=0;d<prototypes_final[j].atom_decorations_equivalent.size();d++){ unique_decorations.push_back(prototypes_final[j].atom_decorations_equivalent[d][0]); }
        ss_out << "  " << setw(structure_spacing) << std::left << "atom_decorations_unique="+aurostd::joinWDelimiter(unique_decorations,",") << endl;
      }

      // ---------------------------------------------------------------------------
      // print properties of representative structure (database comparisons only)
      if(prototypes_final[j].structure_representative->properties.size()!=0){
        ss_out << "  " << setw(structure_spacing) << std::left << "structure";
        ss_out << setw(misfit_spacing) << std::right << "misfit";
        for(uint l=0;l<prototypes_final[j].structure_representative->properties_names.size();l++){
          if(prototypes_final[j].structure_representative->properties_units[l].size()!=0){
            ss_out << setw(property_spacing) << std::right
              << prototypes_final[j].structure_representative->properties_names[l]+"("+prototypes_final[j].structure_representative->properties_units[l]+")";
          }
          else {ss_out << setw(property_spacing) << std::right << prototypes_final[j].structure_representative->properties_names[l];}
        }
        ss_out << endl;
        ss_out << dash_line_separator << endl;
      }

      // ---------------------------------------------------------------------------
      // print representative structure
      ss_out << "  " << setw(structure_spacing) << std::left << "prototype="+prototypes_final[j].structure_representative->name; 
      if(prototypes_final[j].structure_representative->properties.size()!=0){
        ss_out << setw(misfit_spacing) << std::right << "-";
        for(uint l=0;l<prototypes_final[j].structure_representative->properties.size();l++){
          ss_out << setw(property_spacing) << std::right << prototypes_final[j].structure_representative->properties[l];
        }
      }
      ss_out << endl;
      ss_out << dash_line_separator << endl;

      // ---------------------------------------------------------------------------
      // print duplicate structures [header]
      if(prototypes_final[j].structures_duplicate.size()!=0 && prototypes_final[j].structure_representative->properties.size()==0){
        ss_out << "  " << setw(structure_spacing) << std::left << "list of duplicates";
        ss_out << setw(misfit_spacing) << std::right << "misfit";
      }
      else if(prototypes_final[j].structures_duplicate.size()==0){
        ss_out << "  " << setw(structure_spacing) << std::left << "no duplicates";
      }
      if(prototypes_final[j].structure_representative->properties.size()==0){
        for(uint l=0;l<prototypes_final[j].property_names.size();l++){
          if(prototypes_final[j].property_units[l].size()!=0){
            ss_out << setw(property_spacing) << std::right
              << prototypes_final[j].property_names[l]+"("+prototypes_final[j].property_units[l]+")";
          }
          else {ss_out << setw(property_spacing) << std::right << prototypes_final[j].property_names[l];}
        }
      }
      ss_out << endl;
      ss_out << dash_line_separator << endl;

      // ---------------------------------------------------------------------------
      // print duplicate structures [content]
      if(prototypes_final[j].structures_duplicate.size()>0){
        for(uint k=0;k<prototypes_final[j].structures_duplicate.size();k++){
          ss_out << "  " << setw(structure_spacing) << std::left << prototypes_final[j].structures_duplicate[k]->name
            << setw(misfit_spacing) << std::right << prototypes_final[j].mapping_info_duplicate[k].misfit;
          if(prototypes_final[j].property_names.size()!=0){
            for(uint l=0;l<prototypes_final[j].structures_duplicate[k]->properties.size();l++){
              ss_out << setw(property_spacing) << std::right << prototypes_final[j].structures_duplicate[k]->properties[l];
            }
          }
          ss_out << endl;
        }
      }
      else {
        ss_out << endl;
      }
    }
  }
  return ss_out.str();
}

// ***************************************************************************
// XtalFinderCalculator::printStructureMappingResults()
// ***************************************************************************
string XtalFinderCalculator::printStructureMappingResults(
    const structure_mapping_info& misfit_info,
    const xstructure& xstr_reference,
    const xstructure& xstr_mapped,
    const string& mode){

  // Prints comprehensive mapping information: misfit, lattice deviation,
  // coordinate displacement, failure, basis transformation, rotation,
  // origin shift, and atom mapping info (matched indices/types, mapping
  // distances, mapping vectors)

  stringstream output;

  if(aurostd::toupper(mode) == "TEXT" || aurostd::toupper(mode) == "TXT"){
    output << endl <<"**************************** MAPPING RESULTS ****************************"<<endl;
    // structural misfit information
    if(misfit_info.misfit<=misfit_match){
      output << endl << "MISFIT:			 " << aurostd::roundoff(misfit_info.misfit,AUROSTD_ROUNDOFF_TOL);
      output << "  STRUCTURES ARE COMPATIBLE (0<misfit<=" << aurostd::roundoff(misfit_match,AUROSTD_ROUNDOFF_TOL) << ")" << endl;
    }
    else if(misfit_info.misfit<=misfit_family){
      output << endl << "MISFIT:       " << aurostd::roundoff(misfit_info.misfit,AUROSTD_ROUNDOFF_TOL);
      output << "  STRUCTURES ARE IN THE SAME FAMILY (" << aurostd::roundoff(misfit_match,AUROSTD_ROUNDOFF_TOL) << "<misfit<=" << aurostd::roundoff(misfit_family,AUROSTD_ROUNDOFF_TOL) << ")" << endl;
    }
    else {
      output << endl <<"MISFIT:			 " << aurostd::roundoff(misfit_info.misfit,AUROSTD_ROUNDOFF_TOL);
      output << "  STRUCTURES ARE INCOMPATIBLE (misfit>" << aurostd::roundoff(misfit_family,AUROSTD_ROUNDOFF_TOL) << ", or no match found)" << endl;
    }
    output << "-------------------------------------------------------------------------"<<endl;
    output << "Lattice Deviation:       " << aurostd::roundoff(misfit_info.lattice_deviation,AUROSTD_ROUNDOFF_TOL) << endl;
    output << "Coordinate Displacement: " << aurostd::roundoff(misfit_info.coordinate_displacement,AUROSTD_ROUNDOFF_TOL) << endl;
    output << "Figure of Failure:       " << aurostd::roundoff(misfit_info.failure,AUROSTD_ROUNDOFF_TOL) << endl;
    // magnetic misfit information
    if(misfit_info.is_magnetic_misfit){
      output << "Figure of Magnetic Displacement:	" << aurostd::roundoff(misfit_info.magnetic_displacement,AUROSTD_ROUNDOFF_TOL) << endl;
      output << "Figure of Magnetic Failure:	    " << aurostd::roundoff(misfit_info.magnetic_failure,AUROSTD_ROUNDOFF_TOL) << endl;
    }
    // transformation information (basis transformation, rotation, and origin shift)
    // (note: transformation info is not available with the supercell method)
    if(!DEFAULT_XTALFINDER_SUPERCELL_METHOD){
      output << "-------------------------------------------------------------------------"<<endl;
      output << "STRUCTURE TRANSFORMATION (test structure -> reference structure)" << endl;
      output << "Volume scaling factor:" << endl;
      output << aurostd::roundoff(misfit_info.rescale_factor,AUROSTD_ROUNDOFF_TOL) << endl;
      output << "Basis Transformation:" << endl;
      output << aurostd::roundoff(misfit_info.basis_transformation,AUROSTD_ROUNDOFF_TOL) << endl;
      output << "Rotation:" << endl;
      output << aurostd::roundoff(misfit_info.rotation,AUROSTD_ROUNDOFF_TOL) << endl;
      output << "Origin Shift: " << ((xstr_mapped.coord_flag==_COORDS_CARTESIAN_) ? "(Cart.)" : "(Frac.)") << endl;
      output << aurostd::roundoff(misfit_info.origin_shift,AUROSTD_ROUNDOFF_TOL) << endl;

      xstructure xstr_transformed = xstr_mapped;
      // rescale transformed structure
      xstr_transformed.InflateVolume(misfit_info.rescale_factor);
      // apply transformations to xstructure
      xstr_transformed.TransformStructure(
          misfit_info.basis_transformation,
          misfit_info.rotation,
          misfit_info.origin_shift,
          (xstr_mapped.coord_flag==_COORDS_FRACTIONAL_)); //DX20210222

      // mapping information
      output << "-------------------------------------------------------------------------"<<endl;
      output << printAtomMappings(misfit_info,xstr_reference,xstr_transformed);
      output << printUnmatchedAtoms(misfit_info,xstr_reference,xstr_transformed);
      // closest matching representation of structures
      output << "-------------------------------------------------------------------------"<<endl;
      output << "FINAL - REFERENCE STRUCTURE: " << endl;	
      output << xstr_reference << endl;
      output << "-------------------------------------------------------------------------"<<endl;
      output << "FINAL - MAPPED STRUCTURE: " << endl;
      output << xstr_transformed;
    }
  }

  return output.str();
}

// ***************************************************************************
// compare::sameStoichiometry()
// ***************************************************************************
namespace compare{
  bool sameStoichiometry(const vector<uint>& stoich1, const vector<uint>& stoich2){

    // Determine if two stoichiometries are equivalent.
    // Stoichiometries must be in the same order to match.

    // quick check
    if(stoich1.size()!=stoich2.size()){return false;}

    for(uint i=0;i<stoich1.size();i++){
      if(stoich1[i]!=stoich2[i]){
        return false;
      }
    }
    return true;
  }
}

// ***************************************************************************
// compare::matchableSpecies()
// ***************************************************************************
namespace compare{
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2,
      const bool& same_species){

    // Determine if it is possible to match species based on the number of
    // atom types (i.e, reduced stoichiometries are equal)

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    deque<int> stoich1; //DX20191125
    deque<int> stoich2; //DX20191125
    if(xstr1.species.size()==xstr2.species.size()){
      if(xstr1.species.size()==1){
        stoich1.push_back(1); stoich2.push_back(1);
      }
      else {
        aurostd::reduceByGCD(xstr1.num_each_type, stoich1); //DX20191125
        aurostd::reduceByGCD(xstr2.num_each_type, stoich2); //DX20191125
      }
      // ---------------------------------------------------------------------------
      // same species: check if atoms and reduced stoichiometries are the same
      if(same_species){
        bool commensurate=false;
        for(uint i=0; i<stoich1.size(); i++){
          for(uint j=0; j<stoich2.size(); j++){
            if(stoich1[i]==stoich2[j] &&
                KBIN::VASP_PseudoPotential_CleanName(xstr1.species[i])==KBIN::VASP_PseudoPotential_CleanName(xstr2.species[j])){ //DX20190329 - remove pseudopotential information
              //cerr << "matching: " << stoich1[i] << "==" << stoich2[j] << " && " << xstr1.species[i] << "==" << xstr2.species[j] << endl;
              commensurate=true;
              break;
            }
          }
          if(!commensurate){ return false; }
        }
        return true;
      }
      // ---------------------------------------------------------------------------
      // ignore species: check if sorted reduced stoichiometries are the same
      else {
        for(uint i=0; i<stoich1.size(); i++){
          std::sort(stoich1.begin(),stoich1.end());
          std::sort(stoich2.begin(),stoich2.end());
        }
        for(uint i=0; i<stoich1.size(); i++){
          if(stoich1[i]!=stoich2[i]){ return false; }
        }
        return true;
      }
    }
    else {
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " NUMBER OF TYPES OF ATOMIC SPECIES IS NOT THE SAME." << endl;
        cerr << " xstr1: " << xstr1.num_each_type.size() << " " << xstr1.title << endl << xstr1 << endl;
        cerr << " xstr2: " << xstr2.num_each_type.size() << " " << xstr2.title << endl << xstr2 << endl;
      }
      return false;
    }
  }
}

// ***************************************************************************
// Same Species 	
// ***************************************************************************
namespace compare{
  bool sameSpecies(const xstructure& xstr1,
      const xstructure& xstr2,
      bool display){

    // Determine if the structures have the same types and counts of species

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    bool VERBOSE = (display && LDEBUG); //DX20191125

    // ---------------------------------------------------------------------------
    // Check number of types
    if(xstr1.num_each_type.size() != xstr2.num_each_type.size()){
      // Display counts
      if(VERBOSE) { //DX20190702 - condense if-statements
        cerr << __AFLOW_FUNC__ << ": Number of element types are not the same."
          << " xstr 1: " << xstr1.num_each_type.size()
          << " and xstr2: " << xstr2.num_each_type.size() << endl;
      }
      return false; //DX20190702 - bug fix, should not be in if-statement
    }

    // ---------------------------------------------------------------------------
    // Check counts and species
    // sort, then check (to check if matchable) = fast
    deque<int> xstr1_num_each_type = xstr1.num_each_type;
    deque<int> xstr2_num_each_type = xstr2.num_each_type;

    std::sort(xstr1_num_each_type.begin(), xstr1_num_each_type.end());
    std::sort(xstr2_num_each_type.begin(), xstr2_num_each_type.end());

    if(xstr1_num_each_type!=xstr2_num_each_type){
      if(VERBOSE) {
        cerr << __AFLOW_FUNC__ << " Number of each type of element are incompatible." << endl;
      }
      return false;
    }
    if(VERBOSE) {
      cerr << __AFLOW_FUNC__ << " Number of each type of element are compatible; proceeding." << endl;
    }

    return true;
  }
}

// ***************************************************************************
// compare::rescaleStructures()
// ***************************************************************************
namespace compare{
  void rescaleStructure(xstructure& xstr1, xstructure& xstr2){

    // If the scaling factors are different, the two structures are rescaled
    // to 1.00

    if(abs(xstr1.scale-xstr2.scale)>0.001){
      xstr1.ReScale(1.0);
      xstr2.ReScale(1.0);
      xstr1.FixLattices();
      xstr2.FixLattices();
    }
  }
}

// ***************************************************************************
// compare::atomicNumberDensity()
// ***************************************************************************
// To compare structure with different volumes
// we rescale the second cell so that the volume divided by
// the number of atoms is the as the first structure.
namespace compare{
  void atomicNumberDensity(const xstructure& xstr1, xstructure& xstr2) {
    double rescale_factor = 1.0;
    return atomicNumberDensity(xstr1, xstr2, rescale_factor);
  }
}

namespace compare{
  void atomicNumberDensity(const xstructure& xstr1,
      xstructure& xstr2,
      double& rescale_factor) {

    //cerr << xstr1.Volume()/xstr1.atoms.size() << " vs " << xstr2.Volume()/xstr2.atoms.size() << endl;
    rescale_factor=(xstr1.Volume()/xstr1.atoms.size())/(xstr2.Volume()/xstr2.atoms.size()); //DX20201215 - store rescale factor
    xstr2.InflateVolume(rescale_factor); //already updates cartesian coordinates
    xstr2.dist_nn_min*=std::pow(rescale_factor,(double) 1/3); //DX20210505 - need to update minimum distance with scaling factor
    // update Cartesian coordinates
    //DX20201210 [OBSOLETE - INFLATE VOLUME ACCOUNTS FOR THIS NOW] for(uint i=0; i<xstr2.atoms.size(); i++){
    //DX20201210 [OBSOLETE - INFLATE VOLUME ACCOUNTS FOR THIS NOW]  xstr2.atoms[i].cpos=F2C(xstr2.lattice,xstr2.atoms[i].fpos);
    //DX20201210 [OBSOLETE - INFLATE VOLUME ACCOUNTS FOR THIS NOW] }
  }
}

// ***************************************************************************
// Fake Atoms Name
// ***************************************************************************
//DX20200728 [OBSOLETE - moved to pflow]

// ***************************************************************************
// compare::printParameters()
// ***************************************************************************
namespace compare{
  void printParameters(const xstructure& xstr, ostream& oss) {

    // Print lattice parameters and volume for the xstructure

    oss << "========================================================" << endl;
    oss << xstr.title << endl;
    oss << "Parameters:" << endl;
    oss << "a:     " << xstr.a << endl;
    oss << "b:     " << xstr.b << endl;
    oss << "c:     " << xstr.c << endl;
    oss << "alpha: " << xstr.alpha << endl;
    oss << "beta:  " << xstr.beta << endl;
    oss << "gamma: " << xstr.gamma << endl;
    oss << "volume:" << xstr.GetVolume() << endl;
  }
}

// ***************************************************************************
// DX20201230 - moved the following functions to aflow_xatom.cpp
// compare::getLeastFrequentAtomSpecies()
// ***************************************************************************

// ***************************************************************************
// compare::sortBySecondPair()
// ***************************************************************************
namespace compare{
  bool sortBySecondPair(const std::pair<string,uint>& a,
      const std::pair<string,uint>& b) {
    return (a.second<b.second);
  }
}

// ***************************************************************************
// compare::sortSpeciesByFrequency()
// ***************************************************************************
namespace compare{
  vector<string> sortSpeciesByFrequency(const xstructure& xstr) {

    vector<std::pair<string,uint> > species_and_counts;

    for(uint i=0;i<xstr.num_each_type.size();i++){
      std::pair<string,uint> tmp_pair;
      tmp_pair.first=xstr.species[i]; tmp_pair.second=xstr.num_each_type[i];
      species_and_counts.push_back(tmp_pair);
    }

    std::stable_sort(species_and_counts.begin(),species_and_counts.end(),compare::sortBySecondPair);

    vector<string> reordered_species;
    for(uint i=0;i<species_and_counts.size();i++){ reordered_species.push_back(species_and_counts[i].first); }

    return reordered_species;
  }
}

// ***************************************************************************
// compare::atomIndicesSortedByFrequency()
// ***************************************************************************
namespace compare{
  vector<uint> atomIndicesSortedByFrequency(const xstructure& xstr) {

    vector<uint> atom_index_sorted;

    vector<string> species_str=sortSpeciesByFrequency(xstr);

    for(uint i=0;i<species_str.size();i++){
      for(uint j=0;j<xstr.atoms.size();j++){
        if(species_str[i]==xstr.atoms[j].name){
          atom_index_sorted.push_back(j);
        }
      }
    }

    return atom_index_sorted;
  }
}

// ***************************************************************************
// compare::similarLatticeParameters()
// ***************************************************************************
namespace compare{
  bool similarLatticeParameters(const xvector<double> d1,
      const xvector<double> d2){

    // Look for 2 corresponding reference frames, check that
    // the length of the 3 vectors and the angles between them are within
    // a given tolerance (in this case 30% of the reference value has been set).

    // ---------------------------------------------------------------------------
    // switch between relative and absolute tolerance
    bool relative = false;

    // ---------------------------------------------------------------------------
    // relative
    if(relative){
      double tol_length=0.3, tol_angle=0.3;

      if( abs(d1(1)-d2(1)) < tol_length*abs(d1(1)) &&
          abs(d1(2)-d2(2)) < tol_length*abs(d1(2)) &&
          abs(d1(3)-d2(3)) < tol_length*abs(d1(3)) &&
          abs(d1(4)-d2(4)) < tol_angle*abs(d1(4)) &&
          abs(d1(5)-d2(5)) < tol_angle*abs(d1(5)) &&
          abs(d1(6)-d2(6)) < tol_angle*abs(d1(6))
        ){
        return true;
      }
      else {	
        return false;
      }
    }
    // ---------------------------------------------------------------------------
    // absolute
    else{
      double tol_length=1.0, tol_angle=5.0; // 1 Angstrom; 5 degrees

      if( abs(d1(1)-d2(1)) < tol_length &&
          abs(d1(2)-d2(2)) < tol_length &&
          abs(d1(3)-d2(3)) < tol_length &&
          abs(d1(4)-d2(4)) < tol_angle &&
          abs(d1(5)-d2(5)) < tol_angle &&
          abs(d1(6)-d2(6)) < tol_angle
        ){
        return true;
      }
      else {	
        return false;
      }
    }
  }
}

// ***************************************************************************
// DX20191122 [MOVED THE FOLLOWING FUNCTIONS TO XATOM]
// resetLatticeDimensions(),
// minimumCoordinationShellLatticeOnly()
// minimumCoordinationShell()
// ***************************************************************************

// ***************************************************************************
// DX20200728 [Moved centroid (NON-PBC and PBC) into XATOM and extended AUROSTD getCentroid()]
// ***************************************************************************

// ***************************************************************************
// XtalFinderCalculator::findMatch()
// ***************************************************************************
bool XtalFinderCalculator::findMatch(
    const xstructure& xstr1,
    const xstructure& xstr2,
    const vector<uint>& atom_indices_xstr1,
    const vector<uint>& atom_indices_xstr2,
    double minimum_interatomic_distance, //DX20200622
    structure_mapping_info& mapping_info,
    bool same_species){ //DX20200910 - added origin_shift

  // To find the best atom matches, the routine computes
  // the distance between an atom in structure1 to all the others,
  // in structure2. The minimum distance is chosen as the match.
  // There are consistency checks to ensure
  //  1) mapped atoms are of similar types,
  //  2) one-to-one atom mappings, and
  //  3) no cross-matching of types.
  // Once the set of mappings is identified, the best origin choice
  // is found by removing the residual from the geometric center.

  // A | 1    A1, A2     I can check which is the best matching
  // B | 2 -> B1, B2  -> for the atom 1 and 2 in the structure
  // C |      C1, C2     with A,B,C,D
  // D |      D1, D2

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  bool VERBOSE=false;


  // ---------------------------------------------------------------------------
  // Determines cutoff distance in which atoms map onto one another and are
  // unlikely to match with other atoms based on the resolution of atoms
  // (i.e., minimum interatomic distance). If mapping distances are below this
  // value, then we do not need to check other possible atom mappings (offering
  // a speed increase).
  // DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING = 4.0 (default)
  // Increase scaling to have a more stringent cutoff (slower, but more robust)
  double _SAFE_MATCH_CUTOFF_ = minimum_interatomic_distance/DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING; //DX20200623

  // sets the size of the mapping info
  uint natoms1 = xstr1.atoms.size();
  uint natoms2 = xstr2.atoms.size();
  compare::resizeMappingInfo(mapping_info, natoms1);


  // magnetic info
  bool is_non_collinear = xstr1.atoms[0].noncoll_spin_is_given; //DX20191213
  bool is_collinear = xstr1.atoms[0].spin_is_given; //DX20191213

  // to store mapping info to check
  vector<uint> mapped_indices_1, mapped_indices_2;
  vector<string> mapped_names_1, mapped_names_2;

  xmatrix<double> lattice=xstr2.lattice;

  vector<xvector<double> > l1, l2, l3;
  vector<int> a_index, b_index, c_index;
  xvector<int> dims(3); //DX20190701 - use robust method
  uint l1_size=0, l2_size=0, l3_size=0;

  // ---------------------------------------------------------------------------
  // declare variables outside of loop (efficiency) //DX20200401
  xvector<double> min_xvec, incell_dist, tmp_xvec, a_component, ab_component;
  uint j=0, k=0, x1=0, x2=0;
  int i2=0;       // indices of atoms (index after sorting)
  double tmp=AUROSTD_MAX_DOUBLE, dist=AUROSTD_MAX_DOUBLE, match_dist=AUROSTD_MAX_DOUBLE, incell_mod=AUROSTD_MAX_DOUBLE;
  double prev_match_dist=0; //to avoid recalculating dims if nothing has changed

  // ---------------------------------------------------------------------------
  // loop
  for(j=0;j<natoms1;j++){
    x1=atom_indices_xstr1[j];

    match_dist=AUROSTD_MAX_DOUBLE;
    prev_match_dist=0; //to avoid recalculating dims if nothing has changed
    dims[1]=dims[2]=dims[3]=0; //reset
    for(k=0;k<natoms2;k++){
      x2 = atom_indices_xstr2[k];
      if(match_dist<prev_match_dist){
        if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
          resetLatticeDimensions(lattice,match_dist,dims,l1,l2,l3,a_index,b_index,c_index);
          prev_match_dist=match_dist;
          l1_size=l1.size(); l2_size=l2.size(); l3_size=l3.size();
        }
      }
      dist=AUROSTD_MAX_DOUBLE;
      incell_dist = xstr1.atoms[x1].cpos-xstr2.atoms[x2].cpos;
      incell_mod = aurostd::modulus(incell_dist);

      if(incell_mod < match_dist){
        if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
          resetLatticeDimensions(lattice,incell_mod,dims,l1,l2,l3,a_index,b_index,c_index);
          prev_match_dist=incell_mod;
          l1_size=l1.size(); l2_size=l2.size(); l3_size=l3.size();
        }
      }
      // ---------------------------------------------------------------------------
      // Find the min distance; thus check distance between neighboring cells to find true minimum.
      // DX - running vector in each loop saves computations; fewer duplicate operations
      if(incell_mod>_SAFE_MATCH_CUTOFF_){
        for(uint m=0;m<l1_size;m++){
          a_component = incell_dist + l1[m];    //DX : coord1-coord2+a*lattice(1)
          for(uint n=0;n<l2_size;n++){
            ab_component = a_component + l2[n]; //DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
            for(uint p=0;p<l3_size;p++){
              tmp_xvec = ab_component + l3[p];  //DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
              tmp=aurostd::modulus(tmp_xvec);
              if(tmp < dist){
                i2 = x2;
                dist = tmp;
                min_xvec = tmp_xvec;
                if(dist<_SAFE_MATCH_CUTOFF_){ break; } // only check if we lower, not every time
              }
              //if(dist<_SAFE_MATCH_CUTOFF_){ break; }
            }
            if(dist<_SAFE_MATCH_CUTOFF_){ break; }
          }
          if(dist<_SAFE_MATCH_CUTOFF_){ break; }
        }
      }
      else{
        if(incell_mod < dist){
          i2 = x2;
          dist = incell_mod;
          min_xvec = incell_dist;
        }
      }
      if(dist<match_dist){
        mapping_info.atom_map[x1]=i2; // perhaps this could just be [x1]?
        mapping_info.basis_map[x1]=xstr2.atoms[i2].type;
        mapping_info.distances_mapped[x1]=dist;
        mapping_info.vectors_mapped[x1]=min_xvec;
        match_dist = dist;
      }
      if(dist<_SAFE_MATCH_CUTOFF_){
        break;
      }
    }

    if(VERBOSE){
      cerr << __AFLOW_FUNC__ << " Mappings information (iterative): " << endl;
      cerr << compare::printAtomMappings(mapping_info) << endl;
    }

    // ---------------------------------------------------------------------------
    // check if mapped atoms are of the same species/type and spin (if applicable)
    if(!compare::consistentAtomMappingType(
          xstr1.atoms[x1],
          xstr2.atoms[mapping_info.atom_map[x1]],
          x1,
          mapping_info.atom_map[x1],
          same_species,
          is_collinear,
          is_non_collinear)){
      return false;
    }

    // ---------------------------------------------------------------------------
    // check if indices used more than once (i.e., many-to-one mapping)
    if(!compare::consistentAtomMappingIndex(
          x1,
          mapping_info.atom_map[x1],
          mapped_indices_1,
          mapped_indices_2)){
      return false;
    }

    // ---------------------------------------------------------------------------
    // check if sets of mapped atoms are consistent (no cross-matching of types)
    if(!compare::consistentAtomSetMappings(
          xstr1.atoms[x1].name,
          xstr2.atoms[mapping_info.atom_map[x1]].name,
          mapped_names_1,
          mapped_names_2)){
      //cerr << "x1: " << x1 << endl;
      //cerr << "mapping_info.atom_map[x1]: " << mapping_info.atom_map[x1] << endl;
      return false;
    }

    // ---------------------------------------------------------------------------
    // store to check consistent index mappings as we go
    mapped_indices_1.push_back(x1);
    mapped_indices_2.push_back(mapping_info.atom_map[x1]);
    // store to check consistent atom set mappings as we go
    mapped_names_1.push_back(xstr1.atoms[x1].name);
    mapped_names_2.push_back(xstr2.atoms[mapping_info.atom_map[x1]].name);
  }

  // ---------------------------------------------------------------------------
  // try to minimize mapping distances //DX20200910
  xvector<double> origin_shift_test;
  vector<xvector<double> > new_mapping_vectors = compare::minimizeMappingDistances(
      mapping_info.vectors_mapped,
      origin_shift_test);

  uint num_distances = new_mapping_vectors.size(); //DX20200922
  vector<double> new_mapping_distances(num_distances); //DX20200922 - set size; fixed
  for(uint i=0;i<num_distances;i++){ new_mapping_distances[i] = aurostd::modulus(new_mapping_vectors[i]); } //DX20200922 - assign with [i], no longer dynamic

  if(LDEBUG){
    for(uint i=0;i<num_distances;i++){ //DX20200922 - use uint instead of vector.size(); efficiency
      cerr << __AFLOW_FUNC__ << " minimum distance: " << mapping_info.distances_mapped[i] << " (before) --> " << new_mapping_distances[i] << " (after)" << endl;
    }
    cerr << __AFLOW_FUNC__ << " sum(distances_orig)=" << aurostd::sum(mapping_info.distances_mapped) << " sum(distances_new)=" << aurostd::sum(new_mapping_distances) << endl;
  }

  // ---------------------------------------------------------------------------
  // check if the distances were minimized or stayed the same
  // if they increased, then something went wrong and we default to the original
  // mapping distances //DX20200910
  if((aurostd::sum(new_mapping_distances)-aurostd::sum(mapping_info.distances_mapped)) < _ZERO_TOL_){
    mapping_info.distances_mapped = new_mapping_distances;
    mapping_info.vectors_mapped = new_mapping_vectors;
    mapping_info.origin_shift -= xstr2.c2f*origin_shift_test;
    mapping_info.origin_shift = BringInCell(mapping_info.origin_shift);
  }
  else{
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " the minimization method did not reduce the mapping distances; use the original mapping distances." << endl;
    }
  }
  return true;
}

// ***************************************************************************
// compare::mimimizeMatchingDistance() //DX20200910
// ***************************************************************************
namespace compare {
  vector<xvector<double> > minimizeMappingDistances(
      const vector<xvector<double> >& distance_vectors){

    xvector<double> origin_shift;
    return minimizeMappingDistances(distance_vectors,origin_shift);
  }
}

namespace compare {
  vector<xvector<double> > minimizeMappingDistances(
      const vector<xvector<double> >& distance_vectors,
      xvector<double>& origin_shift){

    // After identifying the mapped positions (distance vectors for mapped
    // positions), try to minimize the mapping vectors.
    // Do this by seeing if there is a drift/residual between the mapping
    // vectors and subtract the drift (rigid shift) from all positions

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    // ---------------------------------------------------------------------------
    // calculate the drift/residual of the mapping distances (and normalize)
    // via centroid method
    origin_shift = aurostd::getCentroid(distance_vectors);

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " origin_shift " << origin_shift << " mod=" << aurostd::modulus(origin_shift) << endl; }

    // ---------------------------------------------------------------------------
    // subtract off the residuals (rigid shift)
    uint num_distances = distance_vectors.size(); //DX20200922
    vector<xvector<double> > new_distance_vectors(num_distances); //DX20200922 - set size; fixed
    for(uint i=0;i<num_distances;i++){ //DX20200922 - use uint instead of vector.size(); efficiency
      new_distance_vectors[i] = distance_vectors[i]-origin_shift; //DX20200922 - assign with [i], no longer dynamic
    }

    return new_distance_vectors;
  }
}

// ***************************************************************************
// compare::cellDiagonal()
// ***************************************************************************
// Vectorial sums and differences of the lattice basis vectors.
namespace compare{
  void cellDiagonal(const xstructure& xstr,
      vector<double>& diag_sum,
      vector<double>& diag_diff,
      const double& scale) {

    cellDiagonal(xstr.lattice, diag_sum, diag_diff, scale);
  }
}
namespace compare{
  void cellDiagonal(const xmatrix<double>& lattice,
      vector<double>& diag_sum,
      vector<double>& diag_diff,
      const double& scale) {


    xvector<double> origin;

    diag_sum.push_back(abs(distance(origin,(lattice(1)+lattice(2))*scale)));
    diag_sum.push_back(abs(distance(origin,(lattice(1)+lattice(3))*scale)));
    diag_sum.push_back(abs(distance(origin,(lattice(2)+lattice(3))*scale)));
    diag_diff.push_back(abs(distance(lattice(1),lattice(2))*scale));
    diag_diff.push_back(abs(distance(lattice(1),lattice(3))*scale));
    diag_diff.push_back(abs(distance(lattice(2),lattice(3))*scale));
  }
}

// ***************************************************************************
// compare::latticeDeviation()
// ***************************************************************************
namespace compare{
  double latticeDeviation(const vector<double>& diag_sum1,
      const vector<double>& diag_sum2,
      const vector<double>& diag_diff1,
      const vector<double>& diag_diff2) {

    // The lattice deviation is the difference between each face-diagonal,
    // normalized by the face-diagonals of the reference structure.

    uint i=0;
    vector<double> dev;

    for(i=0;i<3;i++)
      //DX20191112 [ORIG] dev.push_back((abs(diag_sum1[i]-diag_sum2[i])+abs(diag_diff1[i]-diag_diff2[i]))/(diag_sum2[i]+diag_diff2[i]));
      dev.push_back((abs(diag_sum2[i]-diag_sum1[i])+abs(diag_diff2[i]-diag_diff1[i]))/(diag_sum1[i]+diag_diff1[i])); //DX20191112 - should be compared to the reference, (1) not (2)?

    double d=1.0;
    for(i=0;i<dev.size();i++)
      d=d*(1.0-dev[i]);
    d=1.0-d;

    return d;
  }
}

// ***************************************************************************
// compare::computeLFAEnvironment()
// ***************************************************************************
// Computes an abridged atomic environment around all LFA atom centers
// only determines the closest distance fore each atom type (i.e., one for each species)
// this is a quick way to look at the atom environment before trying to compare
// atom environment is invariant of unit cell representation/origin choice

// xstructure version
namespace compare{
  vector<AtomEnvironment> computeLFAEnvironment(const xstructure& xstr, bool unique_only){

    // computes an abridged atomic environment around all LFA atom centers
    // only determines the closest distance fore each atom type (i.e., one for each species)
    // this is a quick way to look at the atom environment before trying to compare
    // atom environment is invariant of unit cell representation/origin choice

    // ---------------------------------------------------------------------------
    // determine all LFA atoms in the structure (could be more than one)
    vector<string> LFAs=getLeastFrequentAtomSpecies(xstr);

    // ---------------------------------------------------------------------------
    // compute all LFA environments, looping through each LFA type
    vector<AtomEnvironment> all_environments_LFA;

    for(uint i=0;i<LFAs.size();i++){
      //DX20191122 [OBSOLETE, moved functionality to XATOM] vector<AtomEnvironment> environments_LFA = getUniqueTypesAtomEnvironmentForLFA(xstr, LFAs[i], LFAs);
      vector<AtomEnvironment> environments_LFA = getLFAAtomEnvironments(xstr, LFAs[i], LFAs, ATOM_ENVIRONMENT_MODE_1); //DX20191122
      // ---------------------------------------------------------------------------
      // may have non-primitive cell, but we only want unique/smallest set of information (fast)
      if(unique_only){
        for(uint j=0;j<environments_LFA.size();j++){
          bool duplicate = false;
          for(uint k=0;k<all_environments_LFA.size();k++){
            if(compatibleEnvironments(environments_LFA[j],all_environments_LFA[k],true,false,true)){ //DX20200401 - ignore_environment_angles=false
              duplicate = true;
              break;
            }
          }
          if(!duplicate){
            all_environments_LFA.push_back(environments_LFA[j]);
          }
        }
      }
      // ---------------------------------------------------------------------------
      // primitivized cell, push back all
      else{
        all_environments_LFA.insert(all_environments_LFA.end(),environments_LFA.begin(),environments_LFA.end());
      }
    }

    return all_environments_LFA;
  }
}

// XtalFinderCalculator structure container version
void XtalFinderCalculator::computeLFAEnvironment(structure_container& str_rep, bool unique_only){

  // ---------------------------------------------------------------------------
  // determine all LFA atoms in the structure (could be more than one)
  vector<string> LFAs=getLeastFrequentAtomSpecies(str_rep.structure);

  // ---------------------------------------------------------------------------
  // compute all LFA environments, looping through each LFA type
  vector<AtomEnvironment> all_environments_LFA;

  for(uint i=0;i<LFAs.size();i++){
    //DX20191122 [OBSOLETE, moved functionality to XATOM] vector<AtomEnvironment> environments_LFA = getUniqueTypesAtomEnvironmentForLFA(xstr, LFAs[i], LFAs);
    vector<AtomEnvironment> environments_LFA = getLFAAtomEnvironments(str_rep.structure, LFAs[i], LFAs, ATOM_ENVIRONMENT_MODE_1); //DX20191122
    // ---------------------------------------------------------------------------
    // may have non-primitive cell, but we only want unique/smallest set of information (fast)
    if(unique_only){
      for(uint j=0;j<environments_LFA.size();j++){
        bool duplicate = false;
        for(uint k=0;k<all_environments_LFA.size();k++){
          if(compare::compatibleEnvironments(environments_LFA[j],all_environments_LFA[k],true,false,true)){ //DX20200401 - ignore_environment_angles=false
            duplicate = true;
            break;
          }
        }
        if(!duplicate){
          all_environments_LFA.push_back(environments_LFA[j]);
        }
      }
    }
    // ---------------------------------------------------------------------------
    // primitivized cell, push back all
    else{
      all_environments_LFA.insert(all_environments_LFA.end(),environments_LFA.begin(),environments_LFA.end());
    }
  }

  str_rep.environments_LFA=all_environments_LFA;
}

// ***************************************************************************
// compare::compatibleEnvironmentSets()
// ***************************************************************************
namespace compare{
  bool compatibleEnvironmentSets(const vector<AtomEnvironment>& env_set1,
      const vector<AtomEnvironment>& env_set2,
      bool same_species,
      bool ignore_environment_angles,
      bool exact_match){ //DX20200320 - added environment angles

    // Determines if sets of LFA environments are similar
    // same_species : requires the same atom decorations/types
    // exact_match  : signals if an exact match (remove duplicates) or conversely
    //                if it possible to match later via structure comparison
    // If this used elsewhere, it may need to be modified since it only considers
    // the special case of comparing LFA environments

    // ---------------------------------------------------------------------------
    // if one is not calculated then we need to assume they may be compatible
    if(env_set1.size()==0 || env_set2.size()==0){
      return true;
    }

    // ---------------------------------------------------------------------------
    // check if atom environment set has only one atom
    // (signal more comprehensive environment comparison)
    //DX20200320 [OBSOLETE] bool only_single_LFA_atoms = true;
    //DX20200320 [OBSOLETE] string first_lfa_element = env_set1[0].element_center; //safe to acess element since I checked earlier
    //DX20200320 [OBSOLETE] for(uint i=1;i<env_set1.size();i++){ // start after first 1
    //DX20200320 [OBSOLETE]   if(first_lfa_element ==env_set1[i].element_center){
    //DX20200320 [OBSOLETE]     only_single_LFA_atoms = false;
    //DX20200320 [OBSOLETE]     break;
    //DX20200320 [OBSOLETE]   }
    //DX20200320 [OBSOLETE] }
    //DX20200320 [OBSOLETE] bool compare_frequency = only_single_LFA_atoms;
    // ---------------------------------------------------------------------------
    // if same species
    if(same_species){
      for(uint i=0;i<env_set1.size();i++){
        bool matched_set = false;
        vector<vector<string> > matched_species;
        for(uint j=0;j<env_set2.size();j++){
          matched_set = compatibleEnvironments(env_set1[i],env_set2[j],matched_species,same_species,ignore_environment_angles,exact_match); //DX20200320 - changed compare_frequency to ignore_environment_angles
          if(matched_set) { break; }
        }
        if(!matched_set){ return false; }
      }
    }
    // ---------------------------------------------------------------------------
    // any species
    else{
      // this is a bit more complicated
      // need to think about this; come back later
      // for now return true so we do a robust match
      return true;
    }
    return true;
  }
}

// ***************************************************************************
// compare::compatibleEnvironments()
// ***************************************************************************
namespace compare{
  bool compatibleEnvironments(const AtomEnvironment& env_1,
      const AtomEnvironment& env_2, bool same_species, bool ignore_environment_angles, //DX20200320 - added environment angles
      bool exact_match){

    vector<vector<string> > matched_species;

    return compatibleEnvironments(env_1,env_2,matched_species,same_species,ignore_environment_angles,exact_match); //DX20200320 - added environment angles
  }
}

namespace compare{
  bool compatibleEnvironments(const AtomEnvironment& env_1,
      const AtomEnvironment& env_2,
      vector<vector<string> > & matched_species,
      bool same_species,
      bool ignore_environment_angles,
      bool exact_match){ //DX20200320 - added environment angles

    // Determines if sets of LFA environments are similar
    // same_species     : requires the same atom decorations/types
    // exact_match      : signals if an exact match (remove duplicates) or conversely
    //                    if it possible to match later via structure comparison

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
    bool VERBOSE=false;

    double _TOL_EXACT_MATCH_ = 0.01; // hundredth of an Angstrom, perhaps put in header?
    double _TOL_LOOSE_MATCH_ = 0.2; // ten percent, perhaps put in header? //DX20190724 - changed from 0.25 to 0.1 //DX20200421 - changed to 0.2 with new relative matching
    double max_distance_env1 = aurostd::max(env_1.distances_neighbor); // normalize distances for relative comparisons (needed for volume scaling) //DX20200421
    double max_distance_env2 = aurostd::max(env_2.distances_neighbor); // normalize distances for relative comparisons (needed for volume scaling) //DX20200421

    // ---------------------------------------------------------------------------
    // check for element for center first (fast)
    if(same_species){
      if(env_1.element_center!=env_2.element_center){ return false; }
    }

    // ---------------------------------------------------------------------------
    // check neighboring sites
    for(uint i=0;i<env_1.elements_neighbor.size();i++){
      bool match_found = false;
      vector<string> species;
      for(uint j=0;j<env_2.elements_neighbor.size();j++){
        // ---------------------------------------------------------------------------
        // check same species, if applicable
        if(same_species && env_1.elements_neighbor[i]!=env_2.elements_neighbor[j]){ continue; }

        // ---------------------------------------------------------------------------
        // check frequency of distance
        // this is sensitive to tolerance of cutoff; use with caution
        //DX - THIS IS TOO SENSITIVE - if(compare_frequency && env_1.coordinations_neighbor[i]!=env_2.coordinations_neighbor[j]){ continue; }
        if(exact_match && env_1.coordinations_neighbor[i]!=env_2.coordinations_neighbor[j]){ continue; } //DX20210514 - if exact match, then this is needed

        // ---------------------------------------------------------------------------
        // exact match
        if(exact_match && aurostd::abs(env_1.distances_neighbor[i]-env_2.distances_neighbor[j])<_TOL_EXACT_MATCH_){
          match_found = true; species.push_back(env_2.elements_neighbor[j]);
        }

        // ---------------------------------------------------------------------------
        // relative match
        else if(!exact_match &&
            //aurostd::abs(env_1.distances_neighbor[i]-env_2.distances_neighbor[j])/(env_1.distances_neighbor[i]+env_2.distances_neighbor[j])<_TOL_RELATIVE_MATCH_) //DX20190730 - too strict
            //DX20200416 [OBSOLETE]  TEST aurostd::abs(env_1.distances_neighbor[i]-env_2.distances_neighbor[j])<_TOL_LOOSE_MATCH_) //DX20190730
          aurostd::abs((env_1.distances_neighbor[i]/max_distance_env1)-(env_2.distances_neighbor[j]/max_distance_env2))<_TOL_LOOSE_MATCH_) //DX20200421
          { //CO20200106 - patching for auto-indenting
            match_found = true; species.push_back(env_2.elements_neighbor[j]);
          }

      }
      if(!match_found){ return false; }
      matched_species.push_back(species);
    }

    // ---------------------------------------------------------------------------
    // check relationship between LFA environments
    // only check for on one atom center (cheaper)
    // check this if all checks have passed thus far
    // this is sensitive to tolerance of cutoff; use with caution
    if(!ignore_environment_angles && same_species){ //same species for now
      vector<vector<double> > angles_sets_1 = getAnglesBetweenMixedSpeciesEnvironments(env_1.coordinates_neighbor);
      vector<vector<double> > angles_sets_2 = getAnglesBetweenMixedSpeciesEnvironments(env_2.coordinates_neighbor);

      for(uint i=0;i<angles_sets_1.size();i++){
        // ---------------------------------------------------------------------------
        // use soft-cutoff; if number of coordinates is not equal, then check that
        // the smallest set matches; a better alternative then the hard-cutoff freqency match
        // case 1) if set 1 < set 2
        if(angles_sets_1[i].size()<=angles_sets_2[i].size()){
          for(uint j=0;j<angles_sets_1[i].size();j++){
            bool matched=false;
            for(uint k=0;k<angles_sets_2[i].size();k++){
              if(aurostd::isequal(angles_sets_1[i][j],angles_sets_2[i][k],10.0)) //equal within 10 degrees (//DX20200320 - used to be 20, but 10 is sufficient since we have checkForBetterMatches())
              { //CO20200106 - patching for auto-indenting
                matched=true;
                break;
              }
              if(VERBOSE){cerr << __AFLOW_FUNC__ << " angles dont match: " << angles_sets_1[i][j] << " vs " << angles_sets_2[i][k] << " | diff=" << angles_sets_1[i][j]-angles_sets_2[i][k] << endl;}
            }
            if(!matched){ matched_species.clear(); return false; }
          }
        }
        // ---------------------------------------------------------------------------
        // case 2) if set 1 > set 2
        else if(angles_sets_1[i].size()>angles_sets_2[i].size()){
          for(uint j=0;j<angles_sets_2[i].size();j++){
            bool matched=false;
            for(uint k=0;k<angles_sets_1[i].size();k++){
              if(aurostd::isequal(angles_sets_2[i][j],angles_sets_1[i][k],10.0)){ //equal within 10 degrees
                matched=true;
                break;
              }
            }
            if(!matched){ matched_species.clear(); return false; }
          }
        }
      }
    }

    if(LDEBUG){ cerr << __AFLOW_FUNC__ << " environments are compatible." << endl; }

    return true;
  }
}

// ***************************************************************************
// compare::getAnglesBetweenMixedSpeciesEnvironments()
// ***************************************************************************
namespace compare{
  vector<vector<double> > getAnglesBetweenMixedSpeciesEnvironments(
      const vector<vector<xvector<double> > >& coordinates_neighbor){

    // Gets the angles between mixed species environments

    vector<vector<double> > angles_sets;
    for(uint j=1;j<coordinates_neighbor.size();j++){
      vector<double> angles;
      for(uint i=0;i<coordinates_neighbor[0].size();i++){
        for(uint k=0;k<coordinates_neighbor[j].size();k++){
          angles.push_back(aurostd::angle(coordinates_neighbor[0][i],coordinates_neighbor[j][k])*rad2deg);
        }
      }
      angles_sets.push_back(angles);
    }
    return angles_sets;
  }
}

// ***************************************************************************
// compatible::compatibleNearestNeighborTypesEnvironments()
// ***************************************************************************
namespace compare{
  bool compatibleNearestNeighborTypesEnvironments(
      const vector<vector<double> >& nn_lfa_with_types_1,
      const vector<vector<double> >& nn_lfa_with_types_2,
      int type_match){

    // Determine if LFA nearest neighbor distances are compatible
    // i.e., quick filter for determining if we can match atoms
    // hinges on alphabetic, perhaps make more robust

    bool VERBOSE=false;

    if(VERBOSE){
      for(uint i=0;i<nn_lfa_with_types_1.size();i++){
        cerr <<"1 " << i << ": " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nn_lfa_with_types_1[i]),",") << endl;
      }
      for(uint i=0;i<nn_lfa_with_types_2.size();i++){
        cerr <<"2 " << i << ": " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nn_lfa_with_types_2[i]),",") << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // if same species comparison
    if(type_match==2){
      if(nn_lfa_with_types_1.size()==nn_lfa_with_types_2.size()){
        for(uint i=0;i<nn_lfa_with_types_1.size();i++){
          bool matched_set=false;
          for(uint j=0;j<nn_lfa_with_types_2.size();j++){
            bool matched = true;
            for(uint ii=0;ii<nn_lfa_with_types_1[i].size();ii++){
              if(aurostd::abs(nn_lfa_with_types_1[i][ii]-nn_lfa_with_types_2[j][ii])/(nn_lfa_with_types_1[i][ii]+nn_lfa_with_types_2[j][ii])>0.1){ // less than 10% relative error
                matched=false;
                break;
              }
            }
            if(matched){
              matched_set=true;
            }
          }
          if(!matched_set){
            return false;
          }
        }
      }
      else{
        // come back here to deal with comparing supercells
        // for now, putting true and forcing comparison
        return true;
      }
    }
    // ---------------------------------------------------------------------------
    // if not same species
    else{
      // come back here later (this is a bit more complicated)
      return true;
    }

    return true;
  }
}

// ***************************************************************************
// DX20191122 [MOVED THE FOLLOWING FUNCTIONS INTO XATOM]:
// getUniqueTypesAtomEnvironmentForLFA()
// shortestDistanceRestrictType()
// ***************************************************************************


// ***************************************************************************
// DX20191122 [MOVED THE FOLLOWING FUNCTIONS INTO XPROTO]:
// compare::shortestDistance() -> double NearestNeighborToAtom()
// ***************************************************************************

// ***************************************************************************
// Coordinates Deviation
// ***************************************************************************
namespace compare{
  void coordinateDeviation(
      structure_mapping_info& mapping_info,
      const vector<double>& nn_xstr1,
      const vector<double>& nn_xstr2){

    // Compute the coordinates deviation by looking at each pair of atoms from
    // the reference and mapped structure
    // updates coordinate displacement and fail

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    uint j=0;
    double num=0.0, den=0.0, nfail=0.0;
    double dd=0.0, nn1=0.0, nn2=0.0; //dd=delta distance, nn=nearest neighbor
    int fail1=0, fail2=0;
    for(j=0; j<mapping_info.atom_map.size(); j++){
      nn1 = nn_xstr1[j]; // placement is x1
      nn2 = nn_xstr2[mapping_info.atom_map[j]]; // value at placement is x2
      dd = mapping_info.distances_mapped[j];

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " xstr1 index: " << j << endl;
        cerr << __AFLOW_FUNC__ << " xstr2 index: " << mapping_info.atom_map[j] << endl;
        cerr << __AFLOW_FUNC__ << " nn1: " << nn1 << endl;
        cerr << __AFLOW_FUNC__ << " nn2: " << nn2 << endl;
        cerr << __AFLOW_FUNC__ << " dd: " << dd << endl;
      }
      if(dd<=0.5*nn1){ fail1=0; }
      else { fail1=1; }
      if(dd<=0.5*nn2) fail2=0;
      else { fail2=1; }

      if(fail1==0){
        num=num+dd;
        den=den+nn1;
      }
      if(fail2==0){
        num=num+dd;
        den=den+nn2;
      }
      if(fail1==1) nfail++;
      if(fail2==1) nfail++;
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " cumulative num: " << num << endl;
      cerr << __AFLOW_FUNC__ << " cumulative den: " << den << endl;
    }

    if(den==0){ mapping_info.coordinate_displacement=1; }
    else{ mapping_info.coordinate_displacement=num/den; }

    // consider unmatched atoms
    // DX20201221 - speed up, just check size differences...
    nfail+=nn_xstr1.size()-mapping_info.atom_map.size();
    nfail+=nn_xstr2.size()-mapping_info.atom_map.size();

    mapping_info.failure=(nfail/(nn_xstr1.size()+nn_xstr2.size()));
  }
}

// ***************************************************************************
// compare::magneticDeviation() (BETA FUNCTIONALITY)
// ***************************************************************************
namespace compare{
  void magneticDeviation(
      const xstructure& xstr1,
      const xstructure& xstr2,
      structure_mapping_info& mapping_info){

    // BETA functionality: determine magnetic deviation between magnetic
    // structures

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    double _NON_COLLINEAR_ANGLE_DEGREE_TOL_ = 10.0;
    double magmom_num = 0.0;
    double magmom_den = 0.0;
    uint mag_fail_1 = 0, mag_fail_2 = 0;
    bool is_non_collinear = xstr1.atoms[0].noncoll_spin_is_given;

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " is_non_collinear = " << is_non_collinear << endl;
    }

    for(uint j=0; j<mapping_info.atom_map.size(); j++){

      // ---------------------------------------------------------------------------
      // collinear
      if(!is_non_collinear){
        double magmom_diff = aurostd::abs(xstr1.atoms[j].spin-xstr2.atoms[mapping_info.atom_map[j]].spin);
        if(magmom_diff>_ZERO_TOL_){ //to account for -0.0 vs 0.0
          if(std::signbit(xstr1.atoms[j].spin) == std::signbit(xstr2.atoms[mapping_info.atom_map[j]].spin)){
            magmom_num += magmom_diff;
            magmom_den += aurostd::abs(xstr1.atoms[j].spin)+aurostd::abs(xstr2.atoms[mapping_info.atom_map[j]].spin);
          }
          else{
            mag_fail_1++;
            mag_fail_2++;
          }
        }
        else{
          magmom_num += magmom_diff;
          magmom_den += aurostd::abs(xstr1.atoms[j].spin)+aurostd::abs(xstr2.atoms[mapping_info.atom_map[j]].spin);
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " matching xstr1 (" << j << ") mag=" << xstr1.atoms[j].spin << " to xstr2 (" << mapping_info.atom_map[j] << ") mag=" << xstr2.atoms[mapping_info.atom_map[j]].spin << endl;
          cerr << __AFLOW_FUNC__ << " magmom_diff: " << magmom_diff << endl;
          cerr << __AFLOW_FUNC__ << " spin structure 1 (signbit): " << std::signbit(xstr1.atoms[j].spin) << endl;
          cerr << __AFLOW_FUNC__ << " spin structure 2 (signbit): " << std::signbit(xstr2.atoms[mapping_info.atom_map[j]].spin) << endl;
        }

      }

      // ---------------------------------------------------------------------------
      // non-collinear
      else if(is_non_collinear){
        double magmom_diff = aurostd::modulus(xstr1.atoms[j].noncoll_spin-xstr2.atoms[mapping_info.atom_map[j]].noncoll_spin);
        double angle_between_non_collinear_spins = rad2deg*aurostd::angle(xstr1.atoms[j].noncoll_spin,xstr2.atoms[mapping_info.atom_map[j]].noncoll_spin);
        if(magmom_diff>_ZERO_TOL_){ //to account for -0.0 vs 0.0
          if(angle_between_non_collinear_spins<=_NON_COLLINEAR_ANGLE_DEGREE_TOL_){
            magmom_num += magmom_diff;
            magmom_den += aurostd::modulus(xstr1.atoms[j].noncoll_spin)+aurostd::modulus(xstr2.atoms[mapping_info.atom_map[j]].noncoll_spin);
          }
          else{
            mag_fail_1++;
            mag_fail_2++;
          }
        }
        else{
          magmom_num += magmom_diff;
          magmom_den += aurostd::modulus(xstr1.atoms[j].noncoll_spin)+aurostd::modulus(xstr2.atoms[mapping_info.atom_map[j]].noncoll_spin);
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " matching xstr1 (" << j << ") mag=" << xstr1.atoms[mapping_info.atom_map[j]].noncoll_spin << " to xstr2 (" << mapping_info.atom_map[j] << ") mag=" << xstr2.atoms[mapping_info.atom_map[j]].noncoll_spin << endl;
          cerr << __AFLOW_FUNC__ << " magmom_diff: " << magmom_diff << endl;
          cerr << __AFLOW_FUNC__ << " angle between two structures " << angle_between_non_collinear_spins << endl;
        }
      }
    }

    mapping_info.magnetic_displacement = magmom_num/magmom_den;
    mapping_info.magnetic_failure = (double)(mag_fail_1+mag_fail_2)/(double)(xstr1.atoms.size()+xstr2.atoms.size());
  }
}

// ***************************************************************************
// compare::computeMisfit()
// ***************************************************************************
// Combines differences between all aspects of crystal structure into a
// figure of misfit. (See Burzlaff)
namespace compare{
  double computeMisfit(const structure_mapping_info& mapping_info){
    return computeMisfit(mapping_info.lattice_deviation,mapping_info.coordinate_displacement,mapping_info.failure);
  }
}

namespace compare{
  double computeMisfit(double dev, double dis, double fail) {

    return 1-((1-dev)*(1-dis)*(1-fail));
  }
}

// ***************************************************************************
// compare::computeMisfitMagnetic() (BETA functionality)
// ***************************************************************************
// Combines differences between all aspects of crystal structure into a
// figure of misfit. (See Burzlaff)
// It also adds in a new magnetic contribution (created by us), in the
// same spirit as the Burzlaff criteria
// BETA FUNCTIONALITY
namespace compare{
  double computeMisfitMagnetic(const structure_mapping_info& mapping_info){

    return computeMisfitMagnetic(
        mapping_info.lattice_deviation,
        mapping_info.coordinate_displacement,
        mapping_info.failure,
        mapping_info.magnetic_displacement,
        mapping_info.magnetic_failure);
  }
}

namespace compare{
  double computeMisfitMagnetic(
      double dev, double dis, double fail, double mag_dis, double mag_fail) {

    return 1-((1-dev)*(1-dis)*(1-fail)*(1-mag_dis)*(1-mag_fail));
  }
}

// ***************************************************************************
// XtalFinderCalculator::printAtomMappings() //DX20201214
// ***************************************************************************
string XtalFinderCalculator::printAtomMappings(
    const structure_mapping_info& misfit_info,
    const xstructure& xstr1,
    const xstructure& xstr2){

  // Print the atom mapping information: mapping indices, types, distances
  // and vectors

  stringstream output, index_map, element_map;

  output << std::setiosflags(std::ios::fixed | std::ios::left);

  // header line
  output << std::setw(15) << "Indices"
    << std::setw(15) << "Types"
    << std::setw(25) << "Distances (Angst.)"
    << "Mapping vector (Cart.)" << endl;

  // content
  for(uint i=0;i<misfit_info.atom_map.size();i++){
    index_map.str(""); element_map.str("");

    index_map << i << "-" << misfit_info.atom_map[i];
    element_map << xstr1.atoms[i].name << "-" << xstr2.atoms[misfit_info.atom_map[i]].name;

    output << std::setw(15) << index_map.str()
      << std::setw(15) << element_map.str()
      << std::setw(25) << aurostd::roundoff(misfit_info.distances_mapped[i],AUROSTD_ROUNDOFF_TOL)
      << aurostd::roundoff(misfit_info.vectors_mapped[i],AUROSTD_ROUNDOFF_TOL) << endl;
  }
  return output.str();
}

// ***************************************************************************
// compare::printAtomMappings() //DX20201214
// ***************************************************************************
namespace compare{
  string printAtomMappings(const structure_mapping_info& misfit_info){

    stringstream output, index_map;

    output << std::setiosflags(std::ios::fixed | std::ios::left);

    // header line
    output << std::setw(15) << "Indices"
      << std::setw(25) << "Distances (Angst.)"
      << "Mapping vector (Cart.)" << endl;

    // content
    for(uint i=0;i<misfit_info.atom_map.size();i++){
      index_map.str("");

      index_map << i << "-" << misfit_info.atom_map[i];

      output << std::setw(15) << index_map.str()
        << std::setw(25) << misfit_info.distances_mapped[i]
        << misfit_info.vectors_mapped[i] << endl;
    }
    return output.str();
  }
}

// ***************************************************************************
// compare::printUnmatchedAtoms() //DX20201214
// ***************************************************************************
string XtalFinderCalculator::printUnmatchedAtoms(
    const structure_mapping_info& misfit_info,
    const xstructure& xstr1,
    const xstructure& xstr2){

  stringstream output;

  uint i=0, j=0;
  bool is_mapped = false;

  output << "-------------------------------------------------------------------------"<<endl;
  output << "Missing Atoms in Reference structure:"<< endl;
  for(i=0;i<xstr1.atoms.size();i++){
    is_mapped=false;
    for(j=0;j<misfit_info.atom_map.size();j++){
      if(i==j){ is_mapped=true; }
    }
    if(!is_mapped){
      output << "# "<< i << "   " << xstr1.atoms[i].cpos << "   " << xstr1.atoms[i].name << endl;
    }
  }

  output << "Missing Atoms in Mapped structure:" << endl;
  for(i=0;i<xstr2.atoms.size();i++){
    is_mapped=false;
    for(j=0;j<misfit_info.atom_map.size();j++){
      if(misfit_info.atom_map[j]==i){ is_mapped=true; }
    }
    if(!is_mapped){
      output << "# "<< i << "   " << xstr2.atoms[i].cpos << "   " << xstr2.atoms[i].name << endl;
    }
  }
  return output.str();
}

// ***************************************************************************
// DX20201218 [OBSOLETE - FOLLOWING FUNCTIONS ARE GIVEN BY XATOM VARIANTS]
// bringCoordinateInCell()
// atomInCell()
// ***************************************************************************

// ***************************************************************************
// compare::vectorPeriodic() -> isTranslationVector() in XATOM //DX20210316
// ***************************************************************************

// ***************************************************************************
// compare::GetLFASupercell()
// ***************************************************************************
namespace compare{
  xstructure GetLFASupercell(const xstructure& xstr,
      const xvector<int>& dims,
      const string& lfa_name){

    // Build a supercell comprised only of LFA atoms
    // to speed up translation vector search

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    // ---------------------------------------------------------------------------
    // remove all atoms that are not of the LFA type
    xstructure xstr_LFA_only=xstr;
    xstr_LFA_only.ClearSymmetry(); //DX20181022
    uint num_atoms=xstr_LFA_only.atoms.size();
    for(uint i=0;i<num_atoms;i++){
      if(xstr_LFA_only.atoms[i].name!=lfa_name){
        xstr_LFA_only.RemoveAtom(i);
        num_atoms--;
        i--;
      }
    }

    // ---------------------------------------------------------------------------
    // shift an LFA atom to origin
    for(uint a=0;a<xstr_LFA_only.atoms.size();a++){
      if(xstr_LFA_only.atoms[a].name == lfa_name){
        xstr_LFA_only.ShiftOriginToAtom(a);
        break;
      }
    }

    // ---------------------------------------------------------------------------
    // create supercell (fast)
    //DX [OBSOLETE] vector<int> sc2pcMap, pc2scMap;
    //DX [OBSOLETE] bool get_symmetry=false;
    //DX [OBSOLETE] bool get_full_basis=false;
    //DX [OBSOLETE] bool force_supercell_matrix=true;
    //DX [OBSOLETE] xstructure xstr_LFA_supercell=GetSuperCell(xstr_LFA_only,3,0,0,0,3,0,0,0,3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix); //DX20190319 - use supercell matrix in expansion

    // ---------------------------------------------------------------------------
    // create supercell (fast/robust)
    xstructure xstr_LFA_supercell=xstr_LFA_only;
    GenerateGridAtoms(xstr_LFA_supercell,dims(1),dims(2),dims(3));

    // ---------------------------------------------------------------------------
    // update atoms
    xstr_LFA_supercell.atoms = xstr_LFA_supercell.grid_atoms;
    xstr_LFA_supercell.grid_atoms.clear();

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " Number of LFAs in supercell: " << xstr_LFA_supercell.atoms.size() << endl;}

    return xstr_LFA_supercell;
  }
}

// ***************************************************************************
// XtalFinderCalculator::latticeSearch()
// ***************************************************************************
void XtalFinderCalculator::latticeSearch(
    structure_container& xstr_rep,
    structure_container& xstr_match,
    structure_mapping_info& match_info,
    bool same_species,
    bool optimize_match,
    bool scale_volume, //DX20200422
    uint num_proc){ //DX20201123

  // Performs lattice and origin search

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  bool test_one_lfa_only = false; //DX20190318
  bool test_one_origin_only = false; //DX20200715
  //DX - SPEED UP BUT NOT ROBUST - if(type_match==2){ test_one_lfa_only=true;} //DX20190318

  xstructure xstr1 = xstr_rep.structure;
  xstructure xstr2 = xstr_match.structure;
  //vector<xstructure> vprotos; //tmp - will get rid of

  // ---------------------------------------------------------------------------
  // assign fake atom names
  if(!same_species){
    xstr1.DecorateWithFakeElements();
    xstr2.DecorateWithFakeElements();
  }

  // ---------------------------------------------------------------------------
  // scale volumes of structures
  if(scale_volume){compare::atomicNumberDensity(xstr1, xstr2, match_info.rescale_factor);}

  // ensure cleared
  //xstr_match.misfit_info = compare::initialize_misfit_struct();
  bool magnetic_analysis = (xstr1.atoms[0].spin_is_given || xstr1.atoms[0].noncoll_spin_is_given);
  match_info.is_magnetic_misfit=(magnetic_analysis && _CALCULATE_MAGNETIC_MISFIT_); //DX20191218

  // ---------------------------------------------------------------------------
  // determine least-frequently occuring atom type (LFA) for each structure
  // (there may be more than one)
  // perhaps put in _structure_rep object
  vector<string> LFA_str1=getLeastFrequentAtomSpecies(xstr1);
  vector<string> LFA_str2=getLeastFrequentAtomSpecies(xstr2);
  string lfa_str1=LFA_str1[0]; //initialize
  string lfa_str2=LFA_str2[0]; //initialize

  xvector<double> shift_xstr1; //DX20201215

  xvector<double> abc_angles_q1=Getabc_angles(xstr1.lattice,DEGREES); // lattice parameters

  // ---------------------------------------------------------------------------
  // determine supercell size via search radius/dims
  double search_radius = aurostd::max(abc_angles_q1(1),abc_angles_q1(2),abc_angles_q1(3));
  xvector<int> dims = LatticeDimensionSphere(xstr2.lattice,search_radius);

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " lattice search radius: " << search_radius << endl;}
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " lattice dims : " << dims << endl;}

  // ---------------------------------------------------------------------------
  // peform supercell expansion on LFA atoms in structure2
  xstructure xstr_LFA_supercell = compare::GetLFASupercell(xstr2, dims, lfa_str2);

  // ---------------------------------------------------------------------------
  // find possible translation vectors
  vector<xvector<double> > translation_vectors;
  findSimilarTranslationVectors(
      xstr1.lattice,
      xstr_LFA_supercell,
      xstr2,
      translation_vectors);

  // ---------------------------------------------------------------------------
  // build possible lattices
  vector<xmatrix<double> > lattices;
  //vector<xmatrix<double> > clattices;
  vector<double> latt_devs;

  buildSimilarLattices(translation_vectors,
      xstr1.lattice,
      lattices,
      latt_devs,
      optimize_match,
      scale_volume); //DX20200422

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " Number of lattices to compare: " << lattices.size() << endl;}

  if(lattices.size()>0){

    // ---------------------------------------------------------------------------
    // calculate attributes of structure 1 (volume, lattice parameters, nearest neighbor distances, etc.)
    vector<double> all_nn1;
    if(xstr_rep.nearest_neighbor_distances.size()==0){ // use xstr_rep so we only calculate once
      all_nn1 = NearestNeighbors(xstr_rep.structure); // nearest neighbor distances (invariant of origin shifts)
      xstr_rep.nearest_neighbor_distances = all_nn1;
    }
    else{
      all_nn1 = xstr_rep.nearest_neighbor_distances;
    }
    // ---------------------------------------------------------------------------
    // CALCULATED LATER calculate attributes of structure 2 (volume, lattice parameters, nearest neighbor distances, etc.)
    //vector<double> all_nn2;
    //if(xstr_match.nearest_neighbor_distances.size()==0){
    //  all_nn2 = NearestNeighbors(xstr_match.structure); // nearest neighbor distances (invariant of origin shifts)
    //  xstr_match.nearest_neighbor_distances = all_nn2;
    //}
    //else{
    //  all_nn2 = xstr_match.nearest_neighbor_distances;
    //}

    if(DEFAULT_XTALFINDER_SUPERCELL_METHOD){ // default = false
      // ---------------------------------------------------------------------------
      // calculate attributes of structure 1 (volume, lattice parameters, nearest neighbor distances, etc.)
      vector<double> D1,F1;
      compare::cellDiagonal(xstr1,D1,F1,1); // cell diagonals
      // convert to clattice representation
      xstr1.lattice=GetClat(xstr1.a,xstr1.b,xstr1.c,xstr1.alpha,xstr1.beta,xstr1.gamma);
      xmatrix<double> f2c = xstr1.scale*trasp(xstr1.lattice); //DX+ME20210113 - calculate outside, fast

      for(uint iat=0; iat<xstr1.atoms.size(); iat++){
        xstr1.atoms[iat].cpos=f2c*xstr1.atoms[iat].fpos;
      }

      // makes xstr2 a supercell
      GenerateGridAtoms(xstr2,dims(1),dims(2),dims(3));
      // ---------------------------------------------------------------------------
      // update atoms
      xstr2.atoms = xstr2.grid_atoms;
      xstr2.grid_atoms.clear();
    }

    // ---------------------------------------------------------------------------
    // create structure misfit object for each lattice and add lattice deviation
    vector<structure_mapping_info> vstrs_matched;
    for(uint i=0;i<lattices.size();i++){
      structure_mapping_info str_misfit_tmp = compare::initialize_misfit_struct();
      str_misfit_tmp.rescale_factor = match_info.rescale_factor; // all have been rescaled to this factor
      str_misfit_tmp.lattice_deviation = latt_devs[i];
      vstrs_matched.push_back(str_misfit_tmp);
    }

#ifdef AFLOW_MULTITHREADS_ENABLE
    // ---------------------------------------------------------------------------
    // split task into threads
    uint number_of_structures = vstrs_matched.size();
    xthread::xThread xt(num_proc);
    std::function<bool(const uint, const uint, const xstructure&,
        const vector<double>&, const xstructure&, const string&,
        vector<xmatrix<double> >&, vector<structure_mapping_info>&, bool, bool)> search_atom_mappings
      = std::bind(&XtalFinderCalculator::searchAtomMappings, this, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10);
#else
    if (num_proc) {} // Suppress compiler warnings
#endif

    // ---------------------------------------------------------------------------
    // test origin shifts
    for(uint y=0;y<LFA_str2.size();y++){
      for(uint x=0;x<LFA_str1.size();x++){
        lfa_str1=LFA_str1[x];
        lfa_str2=LFA_str2[y];
        if(same_species && lfa_str1 != lfa_str2){ continue;}

        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " LFA (structure 1): " << lfa_str1 << endl;
          cerr << __AFLOW_FUNC__ << " LFA (structure 2): " << lfa_str2 << endl;
        }

        // ---------------------------------------------------------------------------
        // shift representative structure to LFA
        // NEED TO SHIFT origin of xstr1_tmp to one of the LFA (this was missing before and caused ICSD_102428.BCA, and CBA to not match, but they should
        // //DX20200715 - now explore all shifts, cannot just test one
        for(uint i=0;i<xstr1.atoms.size();i++){
          if(xstr1.atoms[i].name==lfa_str1){
            //CART - shift_xstr1 += xstr1.atoms[i].cpos; //need to store all shifts
            //CART - shift_xstr1 = xstr1.f2c*BringInCell(xstr1.c2f*shift_xstr1); // convert to fractional, bring in cell, convert to cartesian
            shift_xstr1 += xstr1.atoms[i].fpos; //need to store all shifts
            shift_xstr1 = BringInCell(shift_xstr1);
            xstr1.ShiftOriginToAtom(i);
            xstr1.BringInCell(1e-10);

            // ---------------------------------------------------------------------------
            // ME20200207 - No need for vector if it's all the same and xstr1 ist const
            // create vector of variables for each thread
            //vector<xstructure> xstr1_for_thread;
            //for(uint n=0; n<num_proc; n++){
            //  xstr1_for_thread.push_back(xstr1);
            //}

#ifdef AFLOW_MULTITHREADS_ENABLE
            // ---------------------------------------------------------------------------
            // threaded (DX20191107 thread pointer)
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " Searching for possible matching structures [THREADED VERSION]" << endl;}
            xt.runPredistributed(number_of_structures, search_atom_mappings,
                xstr1,
                all_nn1,
                xstr2,
                lfa_str2,
                lattices,
                vstrs_matched,
                same_species,
                optimize_match
                );
#else
            // ---------------------------------------------------------------------------
            // non-threaded
            //uint n=0;
            uint start_index=0;
            uint end_index=vstrs_matched.size();  //DX20191107 switching end point convention
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " Searching for possible matching structures [NON-THREADED VERSION]" << endl;}
            searchAtomMappings(
                start_index,
                end_index,
                xstr1,
                all_nn1,
                xstr2, //DX20201208
                lfa_str2,
                lattices,
                vstrs_matched,
                same_species,
                optimize_match);
#endif

            // ---------------------------------------------------------------------------
            // collect misfits and matching structure representations
            for(uint p=0;p<vstrs_matched.size();p++){
              if(!aurostd::isequal(vstrs_matched[p].misfit,AUROSTD_MAX_DOUBLE) && vstrs_matched[p].misfit<=match_info.misfit){ //DX20220406 - check if AUROSTD_MAX_DOUBLE
                match_info = vstrs_matched[p];
                match_info.origin_shift = BringInCell(match_info.origin_shift+shift_xstr1);
                // if xstr2 was given in Cartesian coordinates, convert shift //DX20210116
                if(xstr2.coord_flag==_COORDS_CARTESIAN_){
                  match_info.origin_shift = F2C(trasp(match_info.rotation*trasp((match_info.basis_transformation*xstr2.lattice))),match_info.origin_shift); // convert lattice to new basis, rotate, then do F2C
                }
              }
            }

            // ---------------------------------------------------------------------------
            // quick return if found a match
            if(match_info.misfit<misfit_match && !optimize_match){ //DX20220406 - 0.1 to misfit_match (tunable)
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " Found match (misfit = " << match_info.misfit << ")! Terminating search early." << endl;}
              return;
            }

            // ---------------------------------------------------------------------------
            // quick return if testing only one origin //DX20200715
            //if(!optimize_match && aurostd::isequal(match_info.misfit,AUROSTD_MAX_DOUBLE) && same_species){ test_one_origin_only=true;} //DX20201102 - need type_match==2 (otherwise we don't check different types)
            if(!optimize_match && aurostd::isequal(match_info.misfit,AUROSTD_MAX_DOUBLE)){
              if(same_species){ test_one_origin_only=true;} //DX20190809 - need type match here; otherwise we may miss structure-type matches
              else { break; } //DX20201217 - move on to next LFA
            }
            if(test_one_origin_only){
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " No mapping found. Searched only one origin. Terminating search early." << endl;}
              return;
            }
          }
        }

        // ---------------------------------------------------------------------------
        // quick return if testing only one LFA set
        //DX20190702 - can i do this: if(!optimize_match && minMis==1){ test_one_lfa_only=true;}
        //if(!optimize_match && aurostd::isequal(match_info.misfit,AUROSTD_MAX_DOUBLE) && same_species){ test_one_lfa_only=true;} //DX20190809 - need type match here; otherwise we may miss structure-type matches
        if(!optimize_match && aurostd::isequal(match_info.misfit,AUROSTD_MAX_DOUBLE)){
          if(same_species){ test_one_lfa_only=true;} //DX20190809 - need type match here; otherwise we may miss structure-type matches
          else { break; } //DX20201217 - move on to next LFA
        }
        if(test_one_lfa_only){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " No match found. Searched only one LFA set. Terminating search early." << endl;}
          return;
        }
      }
    }
  }
  //if(aurostd::isdifferent(match_info.misfit,AUROSTD_MAX_DOUBLE)){
  //  //NEW printStructureMappingResults(xstr_rep,xstr_match);
  //}
}

// ***************************************************************************
// compare::supercell2newRepresentation() //DX20201231
// ***************************************************************************
namespace compare {
  xstructure supercell2newRepresentation(const xstructure& xstr_supercell, const xmatrix<double>& lattice){

    // Converts a supercell to a new representation
    // Wrapper for foldAtomsInCell

    // ---------------------------------------------------------------------------
    // make smaller lattice the new lattice in the supercell structure
    // note: lattices[p] are oriented wrt to supercell (it has to be), otherwise could break periodicity
    xstructure proto=xstr_supercell; //DX20190530 - added "_supercell"; more descriptive
    proto.lattice=lattice;

    xmatrix<double> f2c, c2f;
    // ---------------------------------------------------------------------------
    // C2F - (i.e., will provide fractional coordinates wrt to new lattice)
    // AND remove all atoms outside unit cell based on fractional coordinates
    // speed increase: ensure this is in cell before computing F2C
    // (don't calculate unnecessary matrix-vector multiplication)
    // Note: C2F (done later) changes lattice to one that is aligned with Cartesian directions (a along +X, etc.)
    //       this is like rotating the global coordinates, therefore, fpos does not change

    deque<_atom> new_basis_2;
    // ---------------------------------------------------------------------------
    // supercell method : orig, slow
    c2f=inverse(proto.scale*trasp(proto.lattice)); //DX+CO20200429 - calculate outside loop [speed]
    for(uint iat=0;iat<proto.atoms.size();iat++){
      proto.atoms[iat].fpos=c2f*proto.atoms[iat].cpos; //DX+CO20200429 - C2F (matrix inverse + matrix multiplication) -> c2f (matrix multiplication)
      if(atomInCell(proto.atoms[iat],_ZERO_TOL_,1.05,-0.05)){ //DX20210203 - modified to new function inputs; soft cutoff between 1.05 to -0.05, using robust MapAtom later on resulting subset
        new_basis_2.push_back(proto.atoms[iat]);
      }
    }
    xstructure proto_new;
    proto_new.title=proto.title;

    xvector<double> abc_angles=Getabc_angles(proto.lattice,DEGREES);
    xmatrix<double> clattice = GetClat(abc_angles);
    proto_new.lattice=clattice;

    //DX NEW - START =======================
    f2c = trasp(proto_new.lattice); //DX20190717
    c2f = aurostd::inverse(trasp(proto_new.lattice)); //DX20190717
    //DX20190717 [OBSOLETE] xmatrix<double> f2c = trasp(proto.lattice);
    //DX20190717 [OBSOLETE] xmatrix<double> c2f = aurostd::inverse(trasp(proto.lattice));
    bool skew = false;
    double tol=0.01;

    //DX2021029 [OBSOLETE - doesn't work because of pre-filtering above] deque<_atom> new_basis = foldAtomsInCell(new_basis_2, proto.lattice, proto_new.lattice, skew, tol, false);
    deque<_atom> new_basis;
    xvector<double> tmp; //DX20200330 - declare outside loop
    for(uint j=0;j<new_basis_2.size();j++){
      bool duplicate_lattice_point=false;
      for(uint a=0; a<new_basis.size(); a++){
        tmp = BringInCell(new_basis_2[j].fpos,1e-10);
        if(SYM::MapAtom(new_basis[a].fpos,tmp,proto_new.lattice,f2c,skew,tol)){
          duplicate_lattice_point=true;
          break;
        }
      }
      if(duplicate_lattice_point==false){
        new_basis_2[j].fpos = BringInCell(new_basis_2[j].fpos,1e-10);
        new_basis_2[j].cpos = f2c*new_basis_2[j].fpos;
        new_basis.push_back(new_basis_2[j]);
      }
    }
    std::stable_sort(new_basis.begin(),new_basis.end(),sortAtomsNames); //DX20190709 - need to sort now
    for(uint i=0;i<new_basis.size();i++){
      proto_new.AddAtom(new_basis[i],true); //CO20230319 - add by species
    }
    proto = proto_new;

    return proto;
  }
}

// ***************************************************************************
// XtalFinderCalculator::searchAtomMappings()
// ***************************************************************************
bool XtalFinderCalculator::searchAtomMappings(
    uint start_index, uint end_index,
    const xstructure& xstr1,
    const vector<double>& all_nn1,
    const xstructure& xstr2,
    const string& lfa,
    vector<xmatrix<double> >& lattices,
    vector<structure_mapping_info>& vstrs_matched,
    bool same_species,
    bool optimize_match){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  bool VERBOSE=FALSE;
  stringstream message;

  // ---------------------------------------------------------------------------
  // check if magnetic comparison
  bool calculate_magnetic_misfit =(_CALCULATE_MAGNETIC_MISFIT_&&
      ((xstr1.atoms[0].spin_is_given && xstr2.atoms[0].spin_is_given) ||
       (xstr1.atoms[0].noncoll_spin_is_given && xstr2.atoms[0].noncoll_spin_is_given)));

  // ---------------------------------------------------------------------------
  // sort atom index by species frequency (speed increase)
  vector<uint> atom_index_xstr1 = compare::atomIndicesSortedByFrequency(xstr1);

  xstructure xstr2_tmp;
  vector<double> all_nn2_test;
  // ---------------------------------------------------------------------------
  // loop through possible lattices/structures
  for(uint p=start_index;p<end_index;p++){

    xstr2_tmp = xstr2;
    all_nn2_test.clear();

    // ---------------------------------------------------------------------------
    // two possible comparison methods
    // 1) supercell method: find new unit cell representation via supercell
    //    expansion (slow, robust, well-tested)
    // 2) transformation method: perform basis transformation and rotation
    //    to a new unit cell representation (fast, robust, new)
    // 1 = supercell method
    if(DEFAULT_XTALFINDER_SUPERCELL_METHOD){
      xstr2_tmp = compare::supercell2newRepresentation(xstr2, lattices[p]);
      //xstr2_tmp = foldAtomsInCellXstructure(xstr2, lattices[p], false, 0.01, false);
      xstr2_tmp.dist_nn_min=SYM::minimumDistance(xstr2_tmp);
    }
    // 2 = transformation method (default)
    else{
      // ---------------------------------------------------------------------------
      // identify basis transformation from xstr2 to its new lattice
      // then determine the rotation between xstr2's new lattice and xstr1
      compare::getLatticeTransformation(xstr2_tmp.lattice,
          xstr1.lattice,
          lattices[p],
          vstrs_matched[p].basis_transformation,
          vstrs_matched[p].rotation);

      // ---------------------------------------------------------------------------
      // now transform xstr2 into its new representation
      try{
        xstr2_tmp.TransformStructure(
            vstrs_matched[p].basis_transformation,
            vstrs_matched[p].rotation);
      }
      catch(aurostd::xerror& re){
        if(LDEBUG){
          message << "The basis transformation does not preserve crystal periodicity (different number of atoms). Skipping transformation.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_WARNING_);
          continue;
        }
      }
    }

    if(LDEBUG){
      cerr << "compare::structureSearch: Trying structure " << p << endl;
      cerr << "structure=" << xstr2_tmp << endl;
    }

    double minimum_interatomic_distance = aurostd::min(xstr1.dist_nn_min,xstr2_tmp.dist_nn_min); //DX20200622

    if(compare::sameSpecies(xstr2_tmp,xstr1,false)){

      // ---------------------------------------------------------------------------
      // sort atom index by species frequency (speed increase)
      // this will remain the same regardless of origin shifts below
      vector<uint> atom_index_xstr2 = compare::atomIndicesSortedByFrequency(xstr2_tmp);

      bool all_nn_calculated = false;
      xvector<double> shift_xstr2;
      for(uint iat=0; iat<xstr2_tmp.atoms.size();iat++){
        if(xstr2_tmp.atoms[iat].name==lfa){
          shift_xstr2 += -xstr2_tmp.atoms[iat].fpos; //DX20201215
          shift_xstr2 = BringInCell(shift_xstr2); // convert to fractional, bring in cell, convert to cartesian
          xstr2_tmp.ShiftOriginToAtom(iat);
          vstrs_matched[p].origin_shift = shift_xstr2; //DX20201215
          // need to get shift from here
          xstr2_tmp.BringInCell(1e-10);
          if(VERBOSE){
            cerr << "compare::structureSearch: orig structure " << xstr1 << endl;
            cerr << "compare::structureSearch: test structure shifted " << xstr2_tmp << endl;
          }

          // ---------------------------------------------------------------------------
          // find matched atoms
          if(findMatch(
                xstr1,
                xstr2_tmp,
                atom_index_xstr1,
                atom_index_xstr2,
                minimum_interatomic_distance,
                vstrs_matched[p],
                same_species)){

            if(VERBOSE){
              for(uint m=0;m<vstrs_matched[p].atom_map.size();m++){
                cerr << "compare::structureSearch: " << m << " == " << vstrs_matched[p].atom_map[m] << " : dist=" << vstrs_matched[p].distances_mapped[m] << endl;
              }
            }

            // ---------------------------------------------------------------------------
            // Only calculate the NN for the proto if we found suitable matches.
            // Only calculate once, nothing changes between shifts to origin (affine)
            // MAY NOT NEED THIS, TRANSFORMATION DOESN'T (SHOULDN'T) CHANGE NN INFO
            // (Can save up to a few seconds if omitted ...)
            // BUT, we would need to rescale since atomicNumberDensity can change the
            // nn distances
            if(!all_nn_calculated){
              all_nn2_test = NearestNeighbors(xstr2_tmp);
              all_nn_calculated = true;
              //cerr << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(all_nn2_test),",") << endl;
            }

            // ---------------------------------------------------------------------------
            compare::coordinateDeviation(vstrs_matched[p],all_nn1,all_nn2_test);

            // ---------------------------------------------------------------------------
            // calculate misfit
            double mis = AUROSTD_MAX_DOUBLE;
            if(calculate_magnetic_misfit){
              compare::magneticDeviation(xstr1,xstr2_tmp,vstrs_matched[p]);
              mis=vstrs_matched[p].magnetic_misfit=compare::computeMisfitMagnetic(vstrs_matched[p]);
              if(LDEBUG){
                cerr << "with spin: mis,latt_dev,cd,f,mag_dis,mag_fail: "
                  << vstrs_matched[p].magnetic_misfit << ", "
                  << vstrs_matched[p].lattice_deviation << ", "
                  << vstrs_matched[p].coordinate_displacement << ", "
                  << vstrs_matched[p].failure << ", "
                  << vstrs_matched[p].magnetic_displacement << ", "
                  << vstrs_matched[p].magnetic_failure <<  endl;
              }
            }
            else{
              mis=vstrs_matched[p].misfit=compare::computeMisfit(vstrs_matched[p]);
              if(LDEBUG){
                cerr << "mis,latt_dev,cd,f: "
                  << vstrs_matched[p].misfit << ", "
                  << vstrs_matched[p].lattice_deviation << ", "
                  << vstrs_matched[p].coordinate_displacement << ", "
                  << vstrs_matched[p].failure << endl;
              }
            }

            // If we want to simply find a match and not find the best match, return early
            if(mis<misfit_match && !optimize_match) {
              return true;
            }
          }
          // ---------------------------------------------------------------------------
          // could not map atoms with this origin choice
          else{
            if(LDEBUG){ cerr << __AFLOW_FUNC__ << " Could not match atom positions. Try new origin choice." << endl; }
            compare::resetMappingInfo(vstrs_matched[p]); //DX20220406 - need to reset if no maps were found
          }
        }
      }
    }

    // ---------------------------------------------------------------------------
    // transformation does not yield commensurate atom counts
    else{
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " Structure transformation does not yield commensurate atom counts:"
          << " ref=" << aurostd::joinWDelimiter(xstr1.num_each_type,",")
          << " vs test=" << aurostd::joinWDelimiter(xstr2_tmp.num_each_type,",") << endl;
      }
    }
  }
  return true;
}

// ***************************************************************************
// compare::consistentAtomMapping()
// ***************************************************************************
namespace compare{
  bool consistentAtomMappingType(
      const _atom& atom1,
      const _atom& atom2,
      uint index_x1,
      uint index_x2,
      bool same_species,
      bool is_collinear,
      bool is_non_collinear){

    bool VERBOSE=FALSE;

    // ---------------------------------------------------------------------------
    // check species mapping
    if(same_species && (atom1.name != atom2.name)){
      if(VERBOSE){
        cerr << "xstr1 atom " << index_x1 << ": " << atom1.name << "; xstr2 atom " << index_x2 << " " << atom2.name << endl;
        cerr << "xstr1 basis " << index_x1 << ": " << atom1.type << "; xstr2 basis " << index_x2 << " " << atom2.type << endl;
        cerr << __AFLOW_FUNC__ << " WARNING: Matching species are not the same type, throwing out match (same species comparison)" << endl;
      }
      return false;
    }

    // ---------------------------------------------------------------------------
    // check spin mapping
    if(!_CALCULATE_MAGNETIC_MISFIT_){
      // check collinear spin
      if(is_collinear && (aurostd::abs(atom1.spin-atom2.spin)>_SPIN_TOL_)){
        if(VERBOSE){
          cerr << __AFLOW_FUNC__ << " WARNING: Matching atoms do not have the same collinear spin, throwing out match" << endl;
        }
        return false;
      }
      // check non_collinear spin
      else if(is_non_collinear &&
          (aurostd::abs(atom1.noncoll_spin(1)-atom2.noncoll_spin(1))>_SPIN_TOL_ ||
           aurostd::abs(atom1.noncoll_spin(2)-atom2.noncoll_spin(2))>_SPIN_TOL_ ||
           aurostd::abs(atom1.noncoll_spin(3)-atom2.noncoll_spin(3))>_SPIN_TOL_)){
        if(VERBOSE){
          cerr << __AFLOW_FUNC__ << " WARNING: Matching atoms do not have the same non-collinear spin, throwing out match" << endl;
        }
        return false;
      }
    }
    return true;
  }
}

// ***************************************************************************
// compare::consistentAtomMappingIndex()
// ***************************************************************************
namespace compare{
  bool consistentAtomMappingIndex(
      uint index1,
      uint index2,
      const vector<uint>& index1_list,
      const vector<uint>& index2_list){

    bool VERBOSE=FALSE;

    // ---------------------------------------------------------------------------
    // check if index in structure 1 has been mapped to already
    if(aurostd::WithinList(index1_list, index1)){
      if(VERBOSE){
        cerr << "WARNING: STRUCTURE 1: index " << index1
          << " has already been mapped (stored indices: " << aurostd::joinWDelimiter(index1_list,",") << endl;
      }
      return false;
    }
    // ---------------------------------------------------------------------------
    // check if index in structure 2 has been mapped to already
    else if(aurostd::WithinList(index2_list, index2)){
      if(VERBOSE){
        cerr << __AFLOW_FUNC__ << " WARNING: STRUCTURE 2: index " << index2
          << " has already been mapped (stored indices: " << aurostd::joinWDelimiter(index2_list,",") << endl;
      }
      return false;
    }
    return true;
  }
}

// ***************************************************************************
// compare::consistentAtomSetMappings()
// ***************************************************************************
namespace compare{
  bool consistentAtomSetMappings(
      const string& atom1_name,
      const string& atom2_name,
      const vector<string>& vatoms1_name,
      const vector<string>& vatoms2_name){

    bool VERBOSE=FALSE;

    for(uint i=0;i<vatoms1_name.size();i++){
      if(atom1_name == vatoms1_name[i] && atom2_name != vatoms2_name[i]){
        if(VERBOSE){
          cerr << __AFLOW_FUNC__ << " WARNING: Matching one type of atom to more than one type: "
            << atom1_name << " == " << atom2_name << " | "
            << vatoms1_name[i] << " == " << vatoms2_name[i] << endl;
        }
        return false;
      }
    }
    return true;
  }
}

// ***************************************************************************
// XtalFinderCalculator::findSimilarTranslationVectors()
// ***************************************************************************
void XtalFinderCalculator::findSimilarTranslationVectors(
    const xmatrix<double>& q1,
    const xstructure& xstr_LFA_supercell,
    const xstructure& xstr,
    vector<xvector<double> >& lattice_vecs){

  // This function scans through the possible lattice points to find a
  // lattice that is commensurate with the reference structure (xstr1).

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  bool relative_tolerance=true;

  double min_q1_a = 0.0; double max_q1_a = 0.0;
  double min_q1_b = 0.0; double max_q1_b = 0.0;
  double min_q1_c = 0.0; double max_q1_c = 0.0;

  if(relative_tolerance){
    min_q1_a = aurostd::modulus(q1(1))-aurostd::modulus(q1(1))*0.1;
    max_q1_a = aurostd::modulus(q1(1))+aurostd::modulus(q1(1))*0.1;
    min_q1_b = aurostd::modulus(q1(2))-aurostd::modulus(q1(2))*0.1;
    max_q1_b = aurostd::modulus(q1(2))+aurostd::modulus(q1(2))*0.1;
    min_q1_c = aurostd::modulus(q1(3))-aurostd::modulus(q1(3))*0.1;
    max_q1_c = aurostd::modulus(q1(3))+aurostd::modulus(q1(3))*0.1;
  }
  else{
    min_q1_a = aurostd::modulus(q1(1))-1.0;
    max_q1_a = aurostd::modulus(q1(1))+1.0;
    min_q1_b = aurostd::modulus(q1(2))-1.0;
    max_q1_b = aurostd::modulus(q1(2))+1.0;
    min_q1_c = aurostd::modulus(q1(3))-1.0;
    max_q1_c = aurostd::modulus(q1(3))+1.0;
  }

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Lattice parameters: " << aurostd::modulus(q1(1)) << ", " << aurostd::modulus(q1(2)) << ", " << aurostd::modulus(q1(3)) << endl;
    cerr << __AFLOW_FUNC__ << " Modulus search range for lattice vector a: " << min_q1_a << " - " << max_q1_a << endl;
    cerr << __AFLOW_FUNC__ << " Modulus search range for lattice vector b: " << min_q1_b << " - " << max_q1_b << endl;
    cerr << __AFLOW_FUNC__ << " Modulus search range for lattice vector c: " << min_q1_c << " - " << max_q1_c << endl;
  }

  xvector<double> tmp_vec;
  double tmp_mod = 0.0;

  // Search all possible vectors with modulus comparable to one of the lattice vectors
  uint natoms_lfa_supercell = xstr_LFA_supercell.atoms.size(); //DX20201221 - reduce cost in loop
  for(uint i=0;i<natoms_lfa_supercell;i++){
    tmp_vec = xstr_LFA_supercell.atoms[i].cpos;
    tmp_mod = aurostd::modulus(tmp_vec);
    if((tmp_mod <= max_q1_a && tmp_mod >= min_q1_a) ||
        (tmp_mod <= max_q1_b && tmp_mod >= min_q1_b) ||
        (tmp_mod <= max_q1_c && tmp_mod >= min_q1_c)){
      lattice_vecs.push_back(tmp_vec);
    }
  }

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Number of potential lattice vectors: " << lattice_vecs.size() << endl;
  }

  // ---------------------------------------------------------------------------
  // check if vectors preserve crystal periodicity
  // if only one LFA atom in unit cell -> lattice point -> vectors=lattice vectors (by definition)
  // DX20201210 - the vector periodic check can be expensive for large systems
  // (consider speed increase)
  if(aurostd::min(xstr.num_each_type)!=1){
    // Removing non-periodic lattice vectors
    vector<xvector<double> > lattice_vecs_periodic;
    for(uint i=0;i<lattice_vecs.size();i++){
      if(isTranslationVector(xstr,lattice_vecs[i],0.5,false)){
        lattice_vecs_periodic.push_back(lattice_vecs[i]);
      }
    }
    lattice_vecs = lattice_vecs_periodic; //DX20190320
  }
  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Number of lattice vectors (preserves periodicity): " << lattice_vecs.size() << endl;
    for(uint i=0;i<lattice_vecs.size();i++){
      cerr << __AFLOW_FUNC__ << " lattice vector " << i << ": " << lattice_vecs[i] << " (" << aurostd::modulus(lattice_vecs[i]) << ")" << endl;
    }
  }
}

// ***************************************************************************
// XtalFinderCalculator::buildSimilarLattices()
// ***************************************************************************
bool XtalFinderCalculator::buildSimilarLattices(
    const vector<xvector<double> >& translation_vectors,
    const xmatrix<double>& q1,
    vector<xmatrix<double> >& lattices,
    vector<double>& latt_devs,
    bool optimize_match,
    bool scale_volume){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  bool VERBOSE=false;

  // ---------------------------------------------------------------------------
  // sort via smallest misfit for speed up
  bool sort_via_lattice_deviation = true; //DX20190626 - speed increase
  bool relative_tolerance = false; //DX20190703

  vector<double> D1,F1;
  compare::cellDiagonal(q1,D1,F1,1);

  // ---------------------------------------------------------------------------
  // volume of unit cell tolerance
  //DX20200422 - tol_vol used to be 0.1; now if we allow for volume scaling
  // then we make it larger to find matches in the same family misfit range (0.1<=misfit<=0.2)
  double tol_vol=0.1;
  double abs_det_q1=abs(det(q1)); // volume //DX20201130
  if(scale_volume){ tol_vol=(1.0/3.0); }
  double det_tol=tol_vol*abs_det_q1;

  // ---------------------------------------------------------------------------
  // tolerance for lattice vectors: relative or absolute
  double tol_a=1e9; double tol_b=1e9; double tol_c=1e9;
  xvector<double> abc_angles_q1=Getabc_angles(q1,DEGREES); // lattice parameters //DX20201130
  if(relative_tolerance){
    tol_a=abc_angles_q1[1]*0.3;
    tol_b=abc_angles_q1[2]*0.3;
    tol_c=abc_angles_q1[3]*0.3;
  }
  else{
    tol_a=1.0; tol_b=1.0; tol_c=1.0; // 1 Angstrom
  }

  // ---------------------------------------------------------------------------
  // LDEBUG: print lattice tolerances
  if(LDEBUG) {
    if(relative_tolerance){
      cerr << __AFLOW_FUNC__ << " Tolerance for a (Angstroms): " << tol_a << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for b (Angstroms): " << tol_b << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for c (Angstroms): " << tol_c << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for alpha (degrees): " << abc_angles_q1[4]*0.3 << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for beta (degrees): " << abc_angles_q1[5]*0.3 << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for gamma (degrees): " << abc_angles_q1[6]*0.3 << endl;
    }
    else{
      cerr << __AFLOW_FUNC__ << " Tolerance for a (Angstroms): " << tol_a << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for b (Angstroms): " << tol_b << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for c (Angstroms): " << tol_c << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for alpha (degrees): " << 5 << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for beta (degrees): " << 5 << endl;
      cerr << __AFLOW_FUNC__ << " Tolerance for gamma (degrees): " << 5 << endl;
    }
  }

  xmatrix<double> tmp_lattice(3,3);
  //DX20201130 xmatrix<double> tmp_clatt(3,3);

  uint n_translations = translation_vectors.size();

  // ---------------------------------------------------------------------------
  // compute lenths of possible lattices before-hand (faster than on-the-fly)
  int store=0;
  vector<double> translations_mod;
  for(uint i=0;i<n_translations;i++){
    translations_mod.push_back(aurostd::modulus(translation_vectors[i]));
  }
  if(LDEBUG) { cerr << __AFLOW_FUNC__ << " Number of lattice vectors: " << n_translations << endl; }

  // ---------------------------------------------------------------------------
  // build all possible unit cells with combinations of lattice vectors
  // (order matters, hence not upper triangular)
  for(uint i=0;i<n_translations;i++){
    // ---------------------------------------------------------------------------
    // check lattice vector length: a
    if(abs(translations_mod[i]-abc_angles_q1[1])<tol_a){ //check a
      for(uint j=0;j<n_translations;j++){
        if(j!=i){
          // ---------------------------------------------------------------------------
          // check lattice vector length: b
          if(abs(translations_mod[j]-abc_angles_q1[2])<tol_b){ // check b
            for(uint k=0;k<n_translations;k++){
              if(k!=i && k!=j){
                // ---------------------------------------------------------------------------
                // check lattice vector length: c
                if(abs(translations_mod[k]-abc_angles_q1[3])<tol_c){ //check c
                  tmp_lattice = SYM::xvec2xmat(translation_vectors[i],translation_vectors[j],translation_vectors[k]);
                  // ---------------------------------------------------------------------------
                  // check determinant
                  if(abs(abs_det_q1-abs(det(tmp_lattice))) < det_tol){ //check determinant/volume
                    xvector<double> abc_angles_q2=Getabc_angles(tmp_lattice,DEGREES);
                    // ---------------------------------------------------------------------------
                    // check angles
                    //if(compare::checkTolerance(abc_angles_q1,abc_angles_q2)==false)
                    if(compare::similarLatticeParameters(abc_angles_q1,abc_angles_q2)){
                      double tmp_latt_dev = compare::checkLatticeDeviation(abs_det_q1,tmp_lattice,D1,F1);
                      // ---------------------------------------------------------------------------
                      // check lattice deviation (time-saver)
                      // case 1: NOT optimize match: keep lattices with deviation smaller than Burzlaff's matching requirement)
                      //         otherwise, there is no possible way that it could match with anything
                      // case 2: optimize match: keep lattices with deviation smaller than Burzlaff's same-family requirement)
                      // otherwise, there is no possible way that it could match with anything or be in the same-family
                      if((!optimize_match && tmp_latt_dev <= misfit_match) || (optimize_match && tmp_latt_dev <= misfit_family)) { //fast match doesn't care about finding same family information //DX20190318 - removed unique since it doesn't exist yet
                        // ---------------------------------------------------------------------------
                        // now check uniqueness (this is more expensive than checking lattice deviation, hence why it is further in nesting)
                        bool unique = true;
                        uint placement_index = lattices.size(); //DX20190626 //default to the end
                        for(uint t=0;t<lattices.size();t++){
                          if(identical(tmp_lattice,lattices[t],1e-10)){
                            unique=false;
                            break;
                          }
                          // if lattices are equal (no transformation/rotation) between the structures, make this first //DX20210203
                          if(identical(tmp_lattice,q1)){
                            placement_index = 0;
                          }
                          // store placement/index based on lattice deviation
                          else if(sort_via_lattice_deviation && tmp_latt_dev < latt_devs[t] && placement_index == lattices.size()){ // third condition necessary; otherwise it could move down in order
                            placement_index = t;
                          }
                        }
                        if(unique){
                          // ---------------------------------------------------------------------------
                          // order/re-order by minimum lattice deviation (speed increase: more likely to find matches faster)
                          // append to the end
                          if(placement_index==lattices.size()){ // push_back: either it is the first or it goes to the end
                            lattices.push_back(tmp_lattice); // stores original original orientation
                            //DX20201130 tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
                            //DX20201130 clattices.push_back(tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
                            latt_devs.push_back(tmp_latt_dev);
                          }
                          // insert to certain location via index
                          else{
                            lattices.insert(lattices.begin()+placement_index, tmp_lattice); // stores original original orientation
                            //DX20201130 tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
                            //DX20201130 clattices.insert(clattices.begin()+placement_index, tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
                            latt_devs.insert(latt_devs.begin()+placement_index, tmp_latt_dev);
                          }
                          store++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // if any are identity, put it first

  // ---------------------------------------------------------------------------
  // print lattices and corresponding lattice deviation
  if(VERBOSE){
    cerr << "q1: " << q1 << endl;
    cerr << "det(q1): " << det(q1) << endl;
    cerr << "abc angles q1: " << abc_angles_q1 << endl;
    for(uint i=0;i<lattices.size();i++){
      //DX20201130 cerr << __AFLOW_FUNC__ << endl << " lattice: " << endl << lattices[i] << endl << " clattice: " << clattices[i] << endl << " volume: " << det(lattices[i]) << endl << " lattice deviation: " << endl << latt_devs[i] << endl;
      cerr << __AFLOW_FUNC__ << endl << " lattice: " << endl << lattices[i] << endl << " volume: " << det(lattices[i]) << endl << " lattice deviation: " << endl << latt_devs[i] << endl;
      xvector<double> abc_angles_q2=Getabc_angles(lattices[i],DEGREES);
      cerr << "abc angles: " << abc_angles_q2 << endl;
    }
  }

  return true;
}

// ***************************************************************************
// compare::getLatticeTransformations()
// ***************************************************************************
namespace compare{
  void getLatticeTransformations(const xmatrix<double>& lattice_original,
      const xmatrix<double>& lattice_ideal,
      const vector<xmatrix<double> >& candidate_lattices,
      vector<xmatrix<double> >& basis_transformations,
      vector<xmatrix<double> >& rotations){

    //bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    xmatrix<double> basis_transformation, rotation;

    // ---------------------------------------------------------------------------
    // cycle through candidate lattices and determine transformations
    for(uint i=0;i<candidate_lattices.size();i++){

      getLatticeTransformation(lattice_original,
          lattice_ideal,
          candidate_lattices[i],
          basis_transformation,
          rotation);

      basis_transformations.push_back(basis_transformation);
      rotations.push_back(rotation);

    }
  }
}

// ***************************************************************************
// compare::getLatticeTransformation()
// ***************************************************************************
namespace compare{
  void getLatticeTransformation(const xmatrix<double>& lattice_original,
      const xmatrix<double>& lattice_ideal,
      const xmatrix<double>& candidate_lattice,
      xmatrix<double>& basis_transformation,
      xmatrix<double>& rotation){

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    xmatrix<double> deformation;

    // ---------------------------------------------------------------------------
    // calculate the metric tensor of the original lattice
    xmatrix<double> metric_tensor_original = MetricTensor(lattice_original);
    xmatrix<double> metric_tensor_ideal = MetricTensor(lattice_ideal);

    // ---------------------------------------------------------------------------
    // checks to determine if the lattice basis changed (as opposed to rotated)
    // 1) check if the metric tensors are different
    //    DX20201203 - metric tensors do not identify reflections
    // 2) check for negative determinant of basis transformation, indicates reflection
    //    DX2021015 - could have two negatives in determinant
    // 3) check if R*R^T!=I, property of SO(3)
    // 4) check if rotation angle is pi
    xmatrix<double> metric_tensor_candidate = MetricTensor(candidate_lattice);
    xmatrix<double> basis_transformation_tmp = GetBasisTransformation(lattice_original,candidate_lattice);
    double theta=acos((aurostd::trace(basis_transformation_tmp)-1.0)/2.0);

    if(!aurostd::identical(metric_tensor_original,metric_tensor_candidate) ||
        aurostd::det(basis_transformation_tmp) < 0.0 ||
        !aurostd::isidentity(basis_transformation_tmp*trasp(basis_transformation_tmp)) ||
        aurostd::isequal(theta,PI,1e-2)){

      // ---------------------------------------------------------------------------
      // if the volume change is not an integer, the basis transformation may include a deformation
      // component which must be removed
      if(!aurostd::isinteger(aurostd::det(basis_transformation_tmp)) && aurostd::det(basis_transformation_tmp) > 1.0){
        aurostd::polarDecomposition(basis_transformation_tmp, basis_transformation, deformation);
      }
      else{
        basis_transformation = basis_transformation_tmp;
      }

      // ---------------------------------------------------------------------------
      // convert to new lattice
      xmatrix<double> lattice_new = basis_transformation*lattice_original;

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " lattice new: " << lattice_new << endl;
        cerr << __AFLOW_FUNC__ << " lattice ideal: " << lattice_ideal << endl;
        cerr << __AFLOW_FUNC__ << " basis transformation: " << basis_transformation << endl;
        cerr << __AFLOW_FUNC__ << " det(basis transformation): " << aurostd::det(basis_transformation) << endl;

        // ---------------------------------------------------------------------------
        // calculate the metric tensors, they should be equal after the basis
        // transformation
        xmatrix<double> metric_tensor_new = MetricTensor(lattice_new);
        xmatrix<double> metric_tensor_ideal = MetricTensor(lattice_ideal);
        cerr << __AFLOW_FUNC__ << " metric_tensor_new: " << metric_tensor_new << endl;
        cerr << __AFLOW_FUNC__ << " metric_tensor_ideal: " << metric_tensor_ideal << endl;
      }

      // ---------------------------------------------------------------------------
      // then rotate to the ideal lattice
      xmatrix<double> rotation_tmp = trasp(GetRotation(lattice_new,lattice_ideal)); // use trasp for AFLOW convention

      // ---------------------------------------------------------------------------
      // since we are rotating to the ideal lattice, the GetRotation() function
      // may incorporate a "deformation" component in the matrix
      // we can differentiate this with a polar decomposition T=R*U
      // T: original matrix, R: pure rotation, U: deformation matrix
      aurostd::polarDecomposition(rotation_tmp, rotation, deformation);

    }

    // ---------------------------------------------------------------------------
    // if the metric tensors ARE equal: simple rotation between lattices
    else{

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " rotation only (no basis transformation)!" << endl;
      }
      basis_transformation = aurostd::eye<double>();

      // ---------------------------------------------------------------------------
      // since we are rotating to the ideal lattice, the GetRotation() function
      // may incorporate a "deformation" component in the matrix
      // we can differentiate this with a polar decomposition T=R*U
      // T: original matrix, R: pure rotation, U: deformation matrix
      xmatrix<double> rotation_tmp = trasp(GetRotation(lattice_original,lattice_ideal)); // use trasp for AFLOW convention
      aurostd::polarDecomposition(rotation_tmp, rotation, deformation);

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " rotation: " << rotation << endl;
      }

    }
  }
}

// ***************************************************************************
// compare::getTransformedStructures()
// ***************************************************************************
namespace compare{
  vector<xstructure> getTransformedStructures(
      const xstructure& xstr,
      const vector<xmatrix<double> >& basis_transformations,
      const vector<xmatrix<double> >& rotations){

    //bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    vector<xstructure> vxstrs_transformed;
    xstructure xstr_transformed_tmp;

    for(uint i=0;i<basis_transformations.size();i++){
      vxstrs_transformed.push_back(TransformStructure(xstr,basis_transformations[i],rotations[i]));
    }

    return vxstrs_transformed;
  }
}

// ***************************************************************************
// checkLatticeDeviation
// ***************************************************************************
namespace compare{
  double checkLatticeDeviation(
      const double& xstr1_vol,
      const xmatrix<double>& q2,
      const vector<double>& D1,
      const vector<double>& F1){

    double scale=xstr1_vol/(aurostd::abs(aurostd::det(q2)));
    scale=pow(scale,0.3333);
    vector<double> D2,F2;
    cellDiagonal(q2,D2,F2,scale);
    return latticeDeviation(D1,D2,F1,F2);
  }
}

// ***************************************************************************
// XtalFinderCalculator::calculatePrototypeDesignations()
// ***************************************************************************
void XtalFinderCalculator::calculatePrototypeDesignations(
    vector<StructurePrototype>& prototypes,
    uint num_proc){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  message << "Determining the AFLOW standard designation.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

#ifdef AFLOW_MULTITHREADS_ENABLE
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Number of threads=" << num_proc << endl;}
  // ---------------------------------------------------------------------------
  // [THREADED] determine AFLOW standard designation
  xthread::xThread xt(num_proc);
  std::function<void(vector<StructurePrototype>::iterator&)> fn =
    std::bind(&XtalFinderCalculator::getPrototypeDesignations, this, _1);
  xt.run(prototypes, fn);
  // ---------------------------------------------------------------------------
#else
  // NON-THREADED
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " Non-threaded version. Number of threads=" << num_proc << endl;}
  for (auto it = prototypes.begin(); it != prototypes.end(); ++it) getPrototypeDesignations(it);
#endif

}

// ***************************************************************************
// XtalFinderCalculator::getPrototypeDesignations()
// ***************************************************************************
// OBSOLETE ME20220207 - Replace with iterator for xThread
//void XtalFinderCalculator::getPrototypeDesignations(
//    vector<StructurePrototype>::iterator& prototypes){
//
//  // NOTE: This on-the-fly threaded scheme follows the procedure
//  // discussed in AAPL/aflow_aapl_tcond.cpp, developed by M. Esters (ME).
//
//  int i = AUROSTD_MAX_INT;
//  int nstructures = prototypes.size();
//
//  if(task_counter < nstructures){
//#ifdef AFLOW_MULTITHREADS_ENABLE
//    std::unique_lock<std::mutex> lock(_mutex_);
//#endif
//    i = task_counter++;
//  }
//  else {
//    return;
//  }
//
//  while (i < nstructures){
//    anrl::structure2anrl(prototypes[i].structure_representative->structure,false); //DX20190829 - false for do not recalulate symmetry, save time
//
//    prototypes[i].aflow_label = prototypes[i].structure_representative->structure.prototype;
//    prototypes[i].aflow_parameter_list = prototypes[i].structure_representative->structure.prototype_parameter_list;
//    prototypes[i].aflow_parameter_values = prototypes[i].structure_representative->structure.prototype_parameter_values;
//#ifdef AFLOW_MULTITHREADS_ENABLE
//    std::unique_lock<std::mutex> lock(_mutex_);
//#endif
//    i = task_counter++;
//  }
//}

void XtalFinderCalculator::getPrototypeDesignations(
    vector<StructurePrototype>::iterator& proto){
  anrl::structure2anrl((*proto).structure_representative->structure,false); //DX20190829 - false for do not recalulate symmetry, save time
  (*proto).aflow_label = (*proto).structure_representative->structure.prototype;
  (*proto).aflow_parameter_list = (*proto).structure_representative->structure.prototype_parameter_list;
  (*proto).aflow_parameter_values = (*proto).structure_representative->structure.prototype_parameter_values;
}

// ***************************************************************************
// XtalFinderCalculator::getPrototypeDesignationsPreDistributed()
// ***************************************************************************
void XtalFinderCalculator::getPrototypeDesignationsPreDistributed(
    vector<StructurePrototype>& prototypes,
    uint start_index,
    uint end_index){

  // If end index is greater than prototypes.size(), then
  // determine prototoype designation for all structures
  if(end_index > prototypes.size()){ end_index=prototypes.size(); }

  for(uint i=start_index;i<end_index;i++){ //DX20191107 switching end index convention <= vs <
    anrl::structure2anrl(prototypes[i].structure_representative->structure,false); //DX20190829 - false for do not recalulate symmetry, save time
    prototypes[i].aflow_label = prototypes[i].structure_representative->structure.prototype;
    prototypes[i].aflow_parameter_list = prototypes[i].structure_representative->structure.prototype_parameter_list;
    prototypes[i].aflow_parameter_values = prototypes[i].structure_representative->structure.prototype_parameter_values;
  }
}

// ***************************************************************************
// XtalFinderCalculator::calculateMatchingAFLOWPrototypes()
// ***************************************************************************
void XtalFinderCalculator::calculateMatchingAFLOWPrototypes(
    vector<StructurePrototype>& prototypes,
    uint num_proc){

  //bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;

  message << "Determining if representative structures map to any of the AFLOW prototypes.";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  aurostd::xoption vpflow_protos;
  vpflow_protos.flag("COMPARE2PROTOTYPES",TRUE);

  // ---------------------------------------------------------------------------
  // specify catalog
  vpflow_protos.flag("COMPARE2PROTOTYPES::CATALOG",TRUE);
  vpflow_protos.push_attached("COMPARE2PROTOTYPES::CATALOG","all");

  // ---------------------------------------------------------------------------
  // specify number of processors
  vpflow_protos.flag("COMPARISON_OPTIONS::NP",TRUE);
  vpflow_protos.push_attached("COMPARISON_OPTIONS::NP",aurostd::utype2string<uint>(num_proc));

  // ---------------------------------------------------------------------------
  // do not calculate unique atom decorations since this was already done
  vpflow_protos.flag("COMPARE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);

  bool quiet_orig = XHOST.QUIET;
  XHOST.QUIET=true;

  uint number_of_structures = prototypes.size();
#ifdef AFLOW_MULTITHREADS_ENABLE
  // split task into threads
  xthread::xThread xt(num_proc);
  std::function<void(uint, vector<StructurePrototype>&, const aurostd::xoption&)> fn =
    std::bind(&XtalFinderCalculator::getMatchingAFLOWPrototypes, this, _1, _2, _3);
  xt.run(number_of_structures, fn, prototypes, vpflow_protos);
#else
  // NON-THREADED
  for (uint i = 0; i < number_of_structures; i++) getMatchingAFLOWPrototypes(i,prototypes,vpflow_protos);
#endif
  XHOST.QUIET=quiet_orig;
}

// ***************************************************************************
// XtalFinderCalculator::getMatchingAFLOWPrototypes()
// ***************************************************************************
void XtalFinderCalculator::getMatchingAFLOWPrototypes(
    uint i,
    vector<StructurePrototype>& prototypes,
    const aurostd::xoption& vpflow_protos){

  XtalFinderCalculator xtal_finder_protos(misfit_match,misfit_family,*p_FileMESSAGE,num_proc,*p_oss);
  vector<StructurePrototype> matching_protos = xtal_finder_protos.compare2prototypes(prototypes[i].structure_representative->structure, vpflow_protos);
  if(matching_protos.size()>0){
    for(uint j=0;j<matching_protos[0].structures_duplicate.size();j++){
      prototypes[i].matching_aflow_prototypes.push_back(matching_protos[0].structures_duplicate[j]->name);
    }
  }
}

// ******************************************************************************
// XtalFinderCalculator::writeComparisonOutputFile //DX20201229
// ******************************************************************************
void XtalFinderCalculator::writeComparisonOutputFile(const string& output,
    const string& directory,
    filetype format,
    output_file_xtalfinder comparison_mode,
    bool same_species){

  // Writes comparison results to an output file

  stringstream message;

  string file_prefix = "", contents_info = "";

  // ---------------------------------------------------------------------------
  // determine file prefix based on the comparison type (same species or not)
  // and the mode (structures being compared)
  if(same_species){
    contents_info = "materials";
    switch(comparison_mode){
      case(compare_input_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_MATERIAL;
        break;
      case(duplicate_compounds_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_DUPLICATE;
        break;
      case(compare2database_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE;
        contents_info += " in the database";
        break;
      case(compare_database_entries_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE;
        contents_info += " in the database";
        break;
      default:
        // if not specified, write somewhere as opposed to throwing an error
        file_prefix = DEFAULT_XTALFINDER_FILE_MATERIAL;
        break;
    }
  }
  else {
    contents_info = "structures";
    switch(comparison_mode){
      case(compare_input_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_STRUCTURE;
        break;
      case(duplicate_compounds_xf):
        break;
      case(compare2database_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE;
        contents_info += " in the database";
        break;
      case(compare_database_entries_xf):
        file_prefix = DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE;
        contents_info += " in the database";
        break;
      default:
        // if not specified, write somewhere as opposed to throwing an error
        file_prefix = DEFAULT_XTALFINDER_FILE_STRUCTURE;
        break;
    }
  }

  // ---------------------------------------------------------------------------
  // write JSON
  if(format==json_ft){
    aurostd::string2file(output,directory+"/"+file_prefix+".json");
    message << "RESULTS: See " << directory << "/"+file_prefix+".json" << " for list of unique/duplicate " << contents_info << ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
  }
  // ---------------------------------------------------------------------------
  // write TEXT file (human-readable)
  else if(format==txt_ft){
    aurostd::string2file(output,directory+"/"+file_prefix+".out");
    message << "RESULTS: See " << directory << "/"+file_prefix+".out" << " for list of unique/duplicate " << contents_info << ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
  }
  // ---------------------------------------------------------------------------
  // unexpected file specifications, write to logger rather than lose the
  // information
  else{
    message << "Unexpected file specifications. Printing results to the log rather than writing to a file.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
    message << output;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_);
  }
}

// ******************************************************************************
// compare::writePythonScript() //DX20201228
// ******************************************************************************
namespace compare {
  void writePythonScript(ostream& oss){

    // Writes AFLOW-XtalFinder Python script in a subdirectory


    string directory = aurostd::getPWD();
    string xtalfinder_python_subdir = "AFLOW_XTALFINDER_PYTHON";
    string python_dir = directory + "/" + xtalfinder_python_subdir;

    aurostd::DirectoryMake(python_dir);

    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Writing out python script to: "+python_dir, oss, _LOGGER_NOTICE_);

    aurostd::EmbData::save_to_file("aflow_xtalfinder_python.py", python_dir+"/"+"aflow_xtalfinder_python.py");
  }
}

// AFLOW-XtalFinder (identify prototypes and compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
