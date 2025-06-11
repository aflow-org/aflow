// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalFinder (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
#include "aflow_symmetry_spacegroup.h"

//[OBSOLETE ME20220129 - not used] #undef AFLOW_COMPARE_MULTITHREADS_ENABLE

//[OBSOLETE ME20220129 - not used] #if GCC_VERSION >= 40400   // added two zeros
//[OBSOLETE ME20220129 - not used] #define AFLOW_COMPARE_MULTITHREADS_ENABLE 1
//[OBSOLETE ME20220129 - not used] #include <thread>
//[OBSOLETE ME20220129 - not used] #else
//[OBSOLETE ME20220129 - not used] #warning "The multithread parts of AFLOW-XtalFinder will be not included, since they need gcc 4.4 and higher (C++0x support)."
//[OBSOLETE ME20220129 - not used] #endif

  static uint DEFAULT_COMPARE_AFLUX_PAGE_SIZE = 75000;  // Appears to be small enough to prevent timeouts while being large enough to limit the number of requests

// ***************************************************************************
// SEE README_AFLOW_COMPARE.TXT FOR THE FULL LIST OF AFLOW COMMANDS
// ***************************************************************************

// ***************************************************************************
// compare::compareAtomDecorations()
// ***************************************************************************
namespace compare {
  string compareAtomDecorations(istream& input, const aurostd::xoption& vpflow){
    ostringstream oss;
    ofstream FileMESSAGE; //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare_atom_decorations|--unique_atom_decorations [GENERAL_COMPARISON_OPTIONS] < file";
      string options_function_string = "unique_atom_decorations_options: [--print_misfit]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_XTALFINDER_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::compareAtomDecorations()",options);
    }

    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator
    uint ncpus_default = 1; // initialize with one cpu; default
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        ncpus_default,
        oss);

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions("permutation"); //DX20200103
    xtal_finder.getOptions(vpflow, comparison_options);

    // ---------------------------------------------------------------------------
    // FLAG: print misfit values for duplicate permutations
    bool print_misfit=false; //defalut=false
    if(vpflow.flag("COMPARE_PERMUTATION::PRINT")) { print_misfit=true; }

    // ---------------------------------------------------------------------------
    // load structure
    xstructure xstr(input,IOAFLOW_AUTO);

    // ---------------------------------------------------------------------------
    // calculate unique/duplicate permutations
    stringstream results_ss;
    vector<string> unique_permutations = xtal_finder.getUniquePermutations(xstr, comparison_options, results_ss, xtal_finder.num_proc, print_misfit);

    return results_ss.str();
  }
}


// ***************************************************************************
// XtalFinderCalculator::getUniquePermutations() //DX20201201
// ***************************************************************************
vector<string> XtalFinderCalculator::getUniquePermutations(xstructure& xstr, uint num_proc, bool optimize_match){
  bool print_misfit=false;
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions("permutation");
  comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",optimize_match);
  stringstream results_ss;
  return getUniquePermutations(xstr, comparison_options, results_ss, num_proc, print_misfit);
}

vector<string> XtalFinderCalculator::getUniquePermutations(
    xstructure& xstr,
    aurostd::xoption& comparison_options,
    stringstream& results_ss,
    uint num_proc,
    bool print_misfit){

  vector<string> unique_permutations;
  stringstream ss_output; //DX20190506

  //DX20190506 START
  // ---------------------------------------------------------------------------
  // print format
  filetype format = txt_ft;
  if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
    format = txt_ft;
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
    format = json_ft;
  }
  //DX20190506 END

  // ---------------------------------------------------------------------------
  // permutation comparisons must compare the same species
  //bool same_species=true;

  // ---------------------------------------------------------------------------
  // quick check: check if any sites have the same number of atoms; if not, then no need to try comparing
  if(!print_misfit){
    if(!compare::arePermutationsComparableViaComposition(xstr)){ //DX20190624 - put into function
      deque<int> reduced_stoichiometry; aurostd::reduceByGCD(xstr.num_each_type, reduced_stoichiometry); //DX20191125
      deque<uint> reduced_stoichiometry_uint; for(uint i=0;i<reduced_stoichiometry.size(); i++){ reduced_stoichiometry_uint.push_back((uint)reduced_stoichiometry[i]); } //DX20191125
      unique_permutations = getSpeciesPermutedStrings(reduced_stoichiometry_uint); //DX20191125
      if(format==txt_ft){ //DX20190506
        ss_output << "Unique atom decorations (" << unique_permutations.size() << "): " << endl;
        ss_output << " " << aurostd::joinWDelimiter(unique_permutations,"\n ") << endl;
      }
      if(format==json_ft){ //DX20190506
        vector<string> vcontent;
        for(uint j=0;j<unique_permutations.size();j++){ vcontent.push_back("[\""+unique_permutations[j]+"\"]"); } //DX20210517 - fixed printing error
        ss_output << "{\"atom_decorations_equivalent\":[";
        ss_output << aurostd::joinWDelimiter(vcontent,","); //DX20191125 - Vec to Dec //DX20210517 - fixed printing error
        ss_output << "]}" << endl;
      }
      results_ss << ss_output.str();
      return unique_permutations;
    }
  }

  // ---------------------------------------------------------------------------
  // load input structure
  stringstream xstr_ss; xstr_ss << xstr;
  addStructure2container(xstr, "input geometry", xstr_ss.str(), 0, false); // false: so we can find decorations on systems without atoms
  StructurePrototype structure;
  uint container_index = 0;
  setStructureAsRepresentative(structure,container_index);

  // ---------------------------------------------------------------------------
  // get the unique atom decorations for the structure
  XtalFinderCalculator xtal_finder_permutations;
  xtal_finder_permutations.misfit_match = misfit_match; //copy misfit_match
  xtal_finder_permutations.misfit_family = misfit_family; //copy misfit_family
  xtal_finder_permutations.num_proc = num_proc; //copy num_proc

  string misfit_results = "";
  //vector<StructurePrototype> final_permutations = xtal_finder_permutations.compareAtomDecorations(structure,num_proc,comparison_options);
  xtal_finder_permutations.compareAtomDecorations(structure,misfit_results,num_proc,comparison_options);

  // ---------------------------------------------------------------------------
  // print results
  if(format==txt_ft){ //DX20190506
    ss_output << "Unique atom decorations (" << structure.atom_decorations_equivalent.size() << "): " << endl;
    for(uint j=0;j<structure.atom_decorations_equivalent.size();j++){
      ss_output << " " << aurostd::joinWDelimiter(structure.atom_decorations_equivalent[j], " = ") << endl;
    }
  }
  //DX20190506 START
  else if(format==json_ft){
    stringstream sscontent_json;
    vector<string> vcontent_json;
    sscontent_json << "\"atom_decorations_equivalent\":[";
    for(uint j=0;j<structure.atom_decorations_equivalent.size();j++){
      stringstream sstmp;
      sstmp << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(structure.atom_decorations_equivalent[j],"\""),",") << "]";
      vcontent_json.push_back(sstmp.str()); sstmp.str("");
    }
    sscontent_json << aurostd::joinWDelimiter(vcontent_json,",");
    sscontent_json << "]";
    ss_output << "{" << sscontent_json.str() << "}" << endl;
  }
  //DX20190506 END

  // ---------------------------------------------------------------------------
  // print misfit results
  if(print_misfit){
    if(format==txt_ft){ //DX20190506
      ss_output << "Misfit values: " << endl;
      ss_output << misfit_results << endl;
    }
    else if(format==json_ft){ //DX20190506
      ss_output.str(""); // need to clear content abbreviated content from above
      ss_output << misfit_results << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // store unique atom decorations in vector
  for(uint j=0;j<structure.atom_decorations_equivalent.size();j++){
    unique_permutations.push_back(structure.atom_decorations_equivalent[j][0]);
  }

  // update results_ss
  results_ss << ss_output.str();
  return unique_permutations;
}

// ***************************************************************************
// compare::compareInputStructures()
// ***************************************************************************
namespace compare {
  string compareInputStructures(aurostd::xoption& vpflow, istream& input_stream, ostream& logstream){ //ME20210206 - std::cin variant
    string input_string;
    aurostd::stream2string(input_stream, input_string);
    vpflow.push_attached("COMPARE_STRUCTURE::STRUCTURE_STRING", input_string);
    vpflow.pop_attached("COMPARE_STRUCTURE::STRUCTURE_LIST");
    vpflow.pop_attached("COMPARE_STRUCTURE::DIRECTORY");
    vpflow.pop_attached("COMPARE_STRUCTURE::FILE");
    return compareInputStructures(vpflow, logstream);
  }

  string compareInputStructures(const aurostd::xoption& vpflow, ostream& logstream){ //DX20190425 - changed name, more general

    // This function compares multiple structures (i.e., more than two).

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
    ostringstream oss;
    stringstream message;
    ofstream FileMESSAGE;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}  //CO20200508

    uint ncpus_default = 1; // initialize with one cpu; default
    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        ncpus_default,
        logstream);

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
    // get options from command line and pass to comparison options //DX20201218
    xtal_finder.getOptions(vpflow, comparison_options);

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage = "";
      // material-type comparisons
      if(vpflow.flag("COMPARE_MATERIAL")){
        usage="aflow --compare_materials=str1,str2,str3,... | aflow --compare_materials -D <dir_path> | aflow --compare_materials -F=<filename>";
      }
      // structure-type comparisons
      else if(vpflow.flag("COMPARE_STRUCTURE")){
        usage="aflow --compare_structures=str1,str2,str3,... | aflow --compare_structures -D <dir_path> | aflow --compare_structures -F=<filename>";
      }
      vector<string> options, options_general;
      aurostd::string2tokens(GENERAL_XTALFINDER_OPTIONS_LIST,options_general," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      init::ErrorOption("--usage","compare::compareInputStructures()",options);
    }

    // ---------------------------------------------------------------------------
    // distinguish structures coming from directory or file
    string structures_source = ""; // "structure_list", "directory" or "file"

    // from list appended to command
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST").empty()){
      structures_source = "structure_list";
    }
    // from directory
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY").empty()){
      structures_source = "directory";
    }
    // from file
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE").empty()){
      structures_source = "file";
    }
    // ME20210206 - from stream
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_STRING").empty()){
      structures_source = "string";
    }

    // ---------------------------------------------------------------------------
    // FLAG: directory of structures to compare
    vector<string> file_list; //DX20190424
    string directory=".";
    string filename="";
    string structures_string="";  // ME20210206
    //DX20190424 START
    if(structures_source=="structure_list") {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST"),file_list,",");
      message << "List of files to compare: " << aurostd::joinWDelimiter(file_list,",");
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    //DX20190424 END
    else if(structures_source=="directory") {
      directory=vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY");
      if(!aurostd::FileExist(directory)) {
        message << "Unable to locate directory: " << directory << "." << endl;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Comparison directory: " << directory;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    // ---------------------------------------------------------------------------
    // FLAG: file of structures to compare
    else if(structures_source=="file") {
      filename=vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE");
      if(!aurostd::FileExist(filename)) {
        message << "Unable to locate file: " << filename << "." << endl;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Comparison file: " << filename;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    // ---------------------------------------------------------------------------
    // FLAG: structures stored in a string (ME20210206)
    else if(structures_source=="string") {
      structures_string = vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_STRING");
      if (structures_string.empty()) {
        message << "Input string empty.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Structures taken from string.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else {
      message << "Need to specify location of structures to compare: -D <directory> or -F=<filename>." << endl;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species=false;
    if(vpflow.flag("COMPARE_MATERIAL")){  // material comparisons (find duplicates) //DX20190429 - generalized
      same_species=true;
      message << "OPTIONS: Performing material type comparisons (comparing alike atomic species).";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else if (vpflow.flag("COMPARE_STRUCTURE")){ // structure comparisons (find structure prototypes) //DX20190425 - generalized
      same_species=false;
      message << "OPTIONS: Performing structure type comparisons (any atomic species).";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison
    bool magnetic_comparison= vpflow.flag("COMPARE::MAGNETIC");
    vector<string> magmoms_for_systems;
    if(magnetic_comparison){
      string magnetic_info=vpflow.getattachedscheme("COMPARE::MAGNETIC");
      aurostd::string2tokens(magnetic_info,magmoms_for_systems,":");
      message << "OPTIONS: Including magnetic moment information in comparisons. Magnetic input detected for " << magmoms_for_systems.size() << " systems.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    //DX20190425 - added print and screen only flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format (default is to write both), set the relevant bools
    // to false when another is specified
    bool write_txt = false;
    bool write_json = false;
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){ write_txt = true; }
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){ write_json = true; }
    // if not specified, write both by default
    if(!write_txt && !write_json){ write_txt = true; write_json = true; }

    // ---------------------------------------------------------------------------
    // FLAG: print2screen
    bool screen_only = false;
    if(vpflow.flag("COMPARE::SCREEN_ONLY")) {
      screen_only=true;
    }

    // ---------------------------------------------------------------------------
    // FLAG: print mapping details (for two-structure comparison only)
    bool print=false;
    if(vpflow.flag("COMPARE_STRUCTURE::PRINT")) {
      print=true;
    }

    // ---------------------------------------------------------------------------
    // check if two-structure comparison
    if(file_list.size()==2){
      comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::STORE_COMPARISON_LOGS",TRUE);  //DX20200113 - fixed typo
    }

    // ---------------------------------------------------------------------------
    // load structures
    vector<StructurePrototype> prototypes_final;
    if(structures_source=="structure_list") {
      prototypes_final = xtal_finder.compareStructuresFromStructureList(file_list, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }
    else if(structures_source=="directory") {
      prototypes_final = xtal_finder.compareStructuresFromDirectory(directory, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }
    else if(structures_source=="string") {  // ME20210205
      prototypes_final = xtal_finder.compareStructuresFromString(structures_string, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options);
    }
    if(structures_source=="file") {
      prototypes_final = xtal_finder.compareStructuresFromFile(filename, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }

    // ---------------------------------------------------------------------------
    // prepare both JSON and TEXT outputs (we may end up printing both)
    string results_json = xtal_finder.printResults(prototypes_final, same_species, json_ft);
    string results_txt = xtal_finder.printResults(prototypes_final, same_species, txt_ft);

    //DX20190429 - added format options - START
    // ---------------------------------------------------------------------------
    // if only two comparisons and text only, print mismatch information
    if(file_list.size()==2 && prototypes_final.size() > 0){
      // return abbreviated results (i.e., misfit value along with match, same family, or no match text
      if(prototypes_final[0].mapping_info_duplicate.size()==1){
        double misfit_final = prototypes_final[0].mapping_info_duplicate[0].misfit;
        if(misfit_final <= xtal_finder.misfit_match && (misfit_final+1.0)> 1e-3){
          message << misfit_final << " : " << "MATCH" << endl;
        }
        else if(misfit_final > xtal_finder.misfit_match && misfit_final <= xtal_finder.misfit_family){
          message << misfit_final << " : " << "SAME FAMILY" << endl;
        }
        else if(misfit_final > xtal_finder.misfit_family && misfit_final <= 1.0){
          message << misfit_final << " : " << "NOT A MATCH" << endl;
        }
        else{
          message << "UNMATCHABLE" << endl;
        }
        if(XHOST.QUIET){
          oss << message.str();
        }
        else {
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
        }
        if(print){
          oss << xtal_finder.printStructureMappingResults(prototypes_final[0].mapping_info_duplicate[0],
              prototypes_final[0].structure_representative->structure,
              prototypes_final[0].structures_duplicate[0]->structure);
        }
      }
      else{
        message << "UNMATCHABLE" << endl;
        if(XHOST.QUIET){ oss << message.str(); }
        else{
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
        }
      }
      if(screen_only) {
        if(write_json){ return results_json; }
        if(write_txt){ return results_txt; }
      }
      return oss.str();
    }

    //DX20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only) {
      if(write_json){ return results_json; }
      if(write_txt){ return results_txt; }
    }
    // ---------------------------------------------------------------------------
    // write results to files //DX20201229 - consolidated into functions
    else{
      if(write_json){ xtal_finder.writeComparisonOutputFile(results_json, directory, json_ft, compare_input_xf, same_species); }
      if(write_txt){ xtal_finder.writeComparisonOutputFile(results_txt, directory, txt_ft, compare_input_xf, same_species); }
    }

    return oss.str();
  }
}

// ***************************************************************************
// compare::getIsopointalPrototypes - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  string isopointalPrototypes(istream& input, const aurostd::xoption& vpflow){


    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("ISOPOINTAL_PROTOTYPES::USAGE")) {
      string usage="aflow --get_isopointal_prototypes|--isopointal_prototypes|--get_same_symmetry_prototypes < file";
      string options_function_string = "options: [--catalog=aflow|htqc|all]";
      vector<string> options, options_function;
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::isopointalPrototypes()",options);
    }

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    // ---------------------------------------------------------------------------
    // FLAG: catalog (htqc, anrl, or all)
    string catalog="all";
    if(vpflow.flag("ISOPOINTAL_PROTOTYPES::CATALOG")) {
      catalog=aurostd::tolower(vpflow.getattachedscheme("ISOPOINTAL_PROTOTYPES::CATALOG"));
      if(catalog!="htqc" && catalog!="anrl" && catalog!="all"){
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Catalog/library can only be htqc, anrl, or all.",_INPUT_ILLEGAL_);
      }
    }

    // ---------------------------------------------------------------------------
    // get isopointal structures
    // (calculates symmetry of input structure and grabs symmetrically similar prototypes)
    vector<string> isopointal_prototypes = compare::getIsopointalPrototypes(xstr, catalog);

    // ---------------------------------------------------------------------------
    // print format //DX20210208
    bool write_txt = XHOST.vflag_control.flag("PRINT_MODE::TXT");
    bool write_json = XHOST.vflag_control.flag("PRINT_MODE::JSON");
    // if not specified, write text by default
    if(!write_txt && !write_json){ write_txt = true; }

    stringstream ss_output;
    if(write_txt){
      if(isopointal_prototypes.size()==0){ ss_output << "no isopointal prototypes in AFLOW"; }
      else{ ss_output << aurostd::joinWDelimiter(isopointal_prototypes,","); }
    }
    if(write_json){
      if(!ss_output.str().empty()){ ss_output << endl; } // if printing both text and json, add a newline
      ss_output << "{\"prototypes_isopointal\":";
      ss_output << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(isopointal_prototypes,"\""),",") << "]";
      ss_output << "}";
    }

    return ss_output.str();
  }
}

// ***************************************************************************
// compare::getIsopointalPrototypes - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  vector<string> getIsopointalPrototypes(xstructure& xstr, const string& catalog){


    // ---------------------------------------------------------------------------
    // stoichiometry
    vector<uint> stoichiometry = xstr.GetReducedComposition(true);

    // ---------------------------------------------------------------------------
    // symmetry
    if(xstr.space_group_ITC<1 || xstr.space_group_ITC>230){ // don't recalculate symmetry if already calculated
      double use_tol = SYM::defaultTolerance(xstr); //DX20200821
      xstr.SpaceGroup_ITC(use_tol, SG_SETTING_ANRL); //DX20200821 - added ANRL setting
    }
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions);

    // ---------------------------------------------------------------------------
    // extract isopointal prototypes from AFLOW
    vector<string> vlabel;
    vector<uint> prototype_space_groups;
    vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, xstr.space_group_ITC, grouped_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);

    return vlabel;

  }
}

//DX20190314 - added new function - START
// ***************************************************************************
// compare::getMatchingPrototype - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  vector<string> getMatchingPrototypes(xstructure& xstr, const string& catalog){
    aurostd::xoption vpflow_input;
    return getMatchingPrototypes(xstr,vpflow_input,catalog);
  }
}
namespace compare {
  vector<string> getMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow_input, const string& catalog){ //DX20210421 - added vpflow variant

    // Returns the matching prototype label, if any exists

    aurostd::xoption vpflow = vpflow_input; //DX20210421 - transfer input options to vpflow
    vpflow.flag("COMPARE2PROTOTYPES",TRUE);

    // ---------------------------------------------------------------------------
    // specify catalog
    vpflow.flag("COMPARE2PROTOTYPES::CATALOG",TRUE); //DX20190329 - need to make scheme before attaching, otherwise it doesn't work
    vpflow.push_attached("COMPARE2PROTOTYPES::CATALOG",catalog);

    // ---------------------------------------------------------------------------
    // do not calculate unique atom decorations
    vpflow.flag("COMPARE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);

    // ---------------------------------------------------------------------------
    // quiet output
    bool original_quiet = XHOST.QUIET;
    XHOST.QUIET=TRUE;

    // ---------------------------------------------------------------------------
    // compare structure to AFLOW prototypes
    XtalFinderCalculator xtal_finder; //DX20201218
    vector<StructurePrototype> prototypes = xtal_finder.compare2prototypes(xstr,vpflow);

    // ---------------------------------------------------------------------------
    // global quiet back to default
    XHOST.QUIET=original_quiet;

    vector<string> matching_prototypes;
    if(prototypes.size()==0){ return matching_prototypes; } //DX20210421 - protect against no matching entries

    // ---------------------------------------------------------------------------
    // sort, put best matches first
    vector<double> misfit_matched;
    uint placement_index = prototypes[0].structures_duplicate.size();
    for(uint i=0;i<prototypes[0].structures_duplicate.size();i++){
      for(uint j=0;j<misfit_matched.size();j++){
        if(prototypes[0].mapping_info_duplicate[i].misfit<misfit_matched[j]){ placement_index = 0; }
      }
      if(placement_index == prototypes[0].structures_duplicate.size()){
        matching_prototypes.push_back(prototypes[0].structures_duplicate[i]->name);
        misfit_matched.push_back(prototypes[0].mapping_info_duplicate[i].misfit);
      }
      else{
        matching_prototypes.insert(matching_prototypes.begin()+placement_index,prototypes[0].structures_duplicate[i]->name);
        misfit_matched.insert(misfit_matched.begin()+placement_index,prototypes[0].mapping_info_duplicate[i].misfit);
      }
    }
    return matching_prototypes; // duplicates names are prototype labels
  }
}
//DX20190314 - added new function - START

//DX20190314 - added overloads for compare2prototypes - START
// ***************************************************************************
// compare::compare2prototypes - identifies corresponding protos
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ostream& logstream){
    ofstream FileMESSAGE;
    return compare2prototypes(input, vpflow, FileMESSAGE, logstream);
  }
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    uint ncpus_default = 1; // initialize with one cpu; default
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        ncpus_default,
        logstream);
    return xtal_finder.compare2prototypes(xstr,vpflow);
  }
}

// ***************************************************************************
// compare::printMatchingPrototypes - returns list of matching structures
// ***************************************************************************
namespace compare {
  string printMatchingPrototypes(istream& input, const aurostd::xoption& vpflow){

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare2prototypes|--compare2protos [GENERAL_COMPARISON_OPTIONS] [COMPARE2PROTOTYPES_OPTIONS] < file";
      string options_function_string = "compare2protos_options: [--catalog=aflow|htqc|all]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_XTALFINDER_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::printMatchingPrototypes()",options);
    }

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    XtalFinderCalculator xtal_finder;
    return xtal_finder.printMatchingPrototypes(xstr,vpflow);

  }
}

// ***************************************************************************
// XtalFinderCalculator::printMatchingPrototypes - returns list of matching structures
// ***************************************************************************
string XtalFinderCalculator::printMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow){

  vector<StructurePrototype> prototypes = compare2prototypes(xstr,vpflow);

  // ---------------------------------------------------------------------------
  // print results (use XHOST flags since only checking once)
  bool same_species = false; //default for prototypes
  if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){ return printResults(prototypes, same_species, json_ft); }
  else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){ return printResults(prototypes, same_species, txt_ft); }

  // return txt if no format was given
  return printResults(prototypes, same_species, txt_ft);
}

// ***************************************************************************
// XtalFinderCalculator::compare2prototypes - identifies corresponding protos
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compare2prototypes(
    const xstructure& xstrIN,
    const aurostd::xoption& vpflow){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;
  bool quiet = false;

  xstructure xstr = xstrIN; //DX20200226 - copy

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
  getOptions(vpflow, comparison_options);

  // ---------------------------------------------------------------------------
  // FLAG: catalog (htqc, anrl, or all)
  string catalog="all";
  if(vpflow.flag("COMPARE2PROTOTYPES::CATALOG")) {
    catalog=aurostd::tolower(vpflow.getattachedscheme("COMPARE2PROTOTYPES::CATALOG"));
    if(catalog != "aflow" && catalog!="htqc" && catalog!="anrl" && catalog!="all"){
      message << "Catalog/library can only be htqc, anrl, or all.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    message << "OPTIONS: Catalog/library (htqc, aflow, or all): " << catalog << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // same spaces
  bool same_species=false; //compare2prototype: by definition, want to compare backbone structure, i.e., ignore species

  // ---------------------------------------------------------------------------
  // single round of comparisons
  comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);
  comparison_options.flag("COMPARISON_OPTIONS::PRINT_MATCHES_TO_INPUT_ONLY",TRUE);

  // ---------------------------------------------------------------------------
  // add structure to container
  stringstream ss_input; ss_input << xstr;
  addStructure2container(xstr, "input geometry", ss_input.str(), 0, false);

  // check if structure is loaded;
  if(structure_containers.size() != 1){
    message << "Input structure was not loaded. Nothing to compare.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  // ---------------------------------------------------------------------------
  // convert the input structure based on the comparison_options
  // if a conversion is requested, this should be done before the symmetry;
  // otherwise the symmetry analysis could be very slow
  convertStructures(comparison_options, num_proc);

  vector<StructurePrototype> prototypes_final;

  // ---------------------------------------------------------------------------
  // symmetry
  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && (xstr.space_group_ITC<1 || xstr.space_group_ITC>230)){ //DX20190829 - don't recalculate symmetry if already calculated //DX20191220 - put range instead of ==0
    calculateSymmetries(1);  //1: one structure -> one processor
  }
  else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && xstr.space_group_ITC>=1 && xstr.space_group_ITC<=230){ //DX20191220 - put range instead of !=0
    getSymmetryInfoFromXstructure(structure_containers[0]);
  }
  else{
    setSymmetryPlaceholders();
  }

  if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " Wyckoff positions of input structure:" << endl;
    for(uint i=0;i<structure_containers[0].grouped_Wyckoff_positions.size();i++){
      cerr << structure_containers[0].grouped_Wyckoff_positions[i] << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // variables to pass to GetPrototypes functions
  vector<uint> stoichiometry = structure_containers[0].stoichiometry;
  std::sort(stoichiometry.begin(),stoichiometry.end()); // order stoichiometry, so we can match to AFLOW prototypes more easily
  uint space_group_num = structure_containers[0].space_group;
  vector<GroupedWyckoffPosition> grouped_Wyckoff_positions = structure_containers[0].grouped_Wyckoff_positions;

  vector<string> vlabel;
  vector<uint> prototype_space_groups;

  // ---------------------------------------------------------------------------
  // find prototypes based on stoichiometry only
  if(comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
    message << "Load prototypes with the same stoichiometry as the input.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vlabel = aflowlib::GetPrototypesByStoichiometry(stoichiometry, catalog);
  }
  // find prototypes based on stoichiometry and space group only
  else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF")){
    message << "Load prototypes with the same stoichiometry and space group as the input.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<GroupedWyckoffPosition> empty_Wyckoff_positions;
    vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, empty_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);
  }
  // find prototypes based on stoichiometry, space group, and Wyckoff positions only (recommended/default)
  else {
    message << "Load prototypes with the same stoichiometry (" << aurostd::joinWDelimiter(stoichiometry,":") << ")"
      << ", space group (" << space_group_num << "), and Wyckoff positions as the input.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, grouped_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);
  }

  if(vlabel.size()>0){
    message << "Potential compatible prototypes: " << vlabel.size() << " (" << aurostd::joinWDelimiter(vlabel,",") << ").";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }
  else{
    message << "No compatible prototypes found.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    return prototypes_final;
  }

  //DX20190830 - to avoid multiple threads being spun-up (here and in aflow_xproto.cpp), turn of aflow_pthreads
  uint nthreads_original=AFLOW_PTHREADS::MAX_PTHREADS;
  AFLOW_PTHREADS::MAX_PTHREADS=1;

  //compare::addAFLOWPrototypes2StructurePrototypeVector(all_structures, vlabel);
  addAFLOWPrototypes2container(vlabel);

  // ---------------------------------------------------------------------------
  // group into objects based on stoichiometry and symmetry (Pearson and space group)
  vector<StructurePrototype> comparison_schemes = groupStructurePrototypes(
      same_species,
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
      false,
      quiet); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // compare structures
  prototypes_final = runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, quiet); //DX20200103 - condensed booleans to xoptions

  AFLOW_PTHREADS::MAX_PTHREADS = nthreads_original; //DX20190830 - set back to original setting

  // ---------------------------------------------------------------------------
  // return if there are no similar structures
  if(prototypes_final.size()==0){ return prototypes_final; } //DX20190314 originally : return oss.str()

  comparison_schemes.clear();

  // ---------------------------------------------------------------------------
  // get unique atom decorations
  if(comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS")){
    message << "Identifying unique atom decorations for representative structures ...";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    XtalFinderCalculator xtal_finder_permutations;

    for(uint i=0;i<prototypes_final.size();i++){
      // check if xstructure is generated; if not, make it
      if(!prototypes_final[i].structure_representative->is_structure_generated){
        if(!compare::generateStructure(prototypes_final[i].structure_representative->name,prototypes_final[i].structure_representative->source,prototypes_final[i].structure_representative->relaxation_step,prototypes_final[i].structure_representative->structure,*p_oss)){ //DX20200429
          message << "Could not generate structure (" << prototypes_final[i].structure_representative->name << ").";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        }
      }
      if(LDEBUG){ //DX20190601 - added LDEBUG
        cerr << "Finding unique atom decorations for " << prototypes_final[i].structure_representative->name << ".";
      }
      xtal_finder_permutations.clear();
      xtal_finder_permutations.compareAtomDecorations(prototypes_final[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
      //vector<StructurePrototype> final_permutations = xtal_finder_permutations.compareAtomDecorations(prototypes_final[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
      //for(uint j=0;j<final_permutations.size();j++){
      //  vector<string> tmp_permutations;
      //  tmp_permutations.push_back(final_permutations[j].structure_representative->name); //push back representative permutation
      //  for(uint d=0;d<final_permutations[j].structures_duplicate.size();d++){ tmp_permutations.push_back(final_permutations[j].structures_duplicate[d]->name); } //push back equivalent permutations
      //  prototypes_final[i].atom_decorations_equivalent.push_back(tmp_permutations);
      //}
    }
    message << "Unique atom decorations found.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
  }

  return prototypes_final; //DX20190314 - new return type
}

// ***************************************************************************
// compare::isMatchingStructureInDatabase - boolean if match found in database
// ***************************************************************************
namespace compare {
  // load input structure
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream){ ofstream FileMESSAGE; return isMatchingStructureInDatabase(xstrIN, vpflow, FileMESSAGE, logstream); }
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    // ---------------------------------------------------------------------------
    // main compare2database() function
    uint ncpus_default = 1; // initialize with one cpu; default
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        ncpus_default,
        logstream);
    vector<StructurePrototype> prototypes_final = xtal_finder.compare2database(xstrIN, vpflow);

    // ---------------------------------------------------------------------------
    // safety against bad input geometry files
    if(prototypes_final.empty()){
      stringstream message;
      message << "The input geometry file is invalid (could not be read, corrupt, etc.); it could not be compared to the database.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_FILE_ERROR_);
    }

    // ---------------------------------------------------------------------------
    // return if the database contains an equivalent structure to the input
    return (prototypes_final[0].structures_duplicate.size() != 0);

  }
}

// ***************************************************************************
// compare::matchingStructuresInDatabase - returns matching structures in database
// ***************************************************************************
namespace compare {
  // load input structure
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream){ ofstream FileMESSAGE; return matchingStructuresInDatabase(xstrIN, vpflow, FileMESSAGE, logstream); }
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    // ---------------------------------------------------------------------------
    // instantiate class
    uint ncpus_default = 1; // initialize with one cpu; default
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        ncpus_default,
        logstream);

    // ---------------------------------------------------------------------------
    // main compare2database() function
    vector<StructurePrototype> prototypes_final = xtal_finder.compare2database(xstrIN, vpflow);

    vector<matching_structure> matched_database_structures;

    // ---------------------------------------------------------------------------
    // database DOESN'T contain equivalent structure to input
    if(prototypes_final[0].structures_duplicate.size() == 0){ //DX20210122 - "==" not "!="
      return matched_database_structures;
    }

    // ---------------------------------------------------------------------------
    // return equivalent structures to input
    for(uint i=0;i<prototypes_final[0].structures_duplicate.size();i++){
      matching_structure database_entry;
      database_entry.name = prototypes_final[0].structures_duplicate[i]->name;
      database_entry.misfit = prototypes_final[0].mapping_info_duplicate[i].misfit;
      matched_database_structures.push_back(database_entry);
    }

    return matched_database_structures;
  }
}

// ***************************************************************************
// XtalFinderCalculator::compare2database - compares database
// ***************************************************************************
// load input structure
vector<StructurePrototype> XtalFinderCalculator::compare2database(
    const xstructure& xstrIN,
    const aurostd::xoption& vpflow){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  stringstream message;
  string directory = aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
  // get options from vpflow/command line
  getOptions(vpflow,comparison_options); //DX20200103

  // ---------------------------------------------------------------------------
  // single round of comparisons, only want to match to the input structure
  comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

  // ---------------------------------------------------------------------------
  // only compare entries to the input representation, the rest are extraneous comparisons
  comparison_options.flag("COMPARISON_OPTIONS::PRINT_MATCHES_TO_INPUT_ONLY",TRUE);

  // ---------------------------------------------------------------------------
  // FLAG: do not remove unmatched structures from the StructurePrototype Object
  // keeps results of each comparison
  if(vpflow.flag("COMPARE2DATABASE::KEEP_UNMATCHED")) {
    comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
  }

  // ---------------------------------------------------------------------------
  // FLAG: type of comparison (material-type or structure-type)
  bool same_species = true;
  if(vpflow.flag("COMPARE2DATABASE::STRUCTURE")) {
    same_species = false;
    message << "OPTIONS: Structure-type comparison, i.e., ignore atomic species.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // add input structure to container //DX20210219 - put after same_species flag
  stringstream ss_input; ss_input << xstrIN;
  addStructure2container(xstrIN, "input geometry", ss_input.str(), 0, same_species); //DX20210219 - use same_species flag

  // check if structure is loaded;
  if(structure_containers.size() != 1){
    message << "Input structure was not loaded. Nothing to compare.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
  }

  // ---------------------------------------------------------------------------
  // convert the input structure based on the comparison_options
  // if a conversion is requested, this should be done before the symmetry;
  // otherwise the symmetry analysis could be very slow
  convertStructures(comparison_options, num_proc);

  // ---------------------------------------------------------------------------
  // start building AFLUX query
  vector<string> vmatchbook;

  // ---------------------------------------------------------------------------
  // FLAG: property list to extract from database (using AFLUX)
  vector<string> property_list;
  if(vpflow.flag("COMPARE2DATABASE::PROPERTY_LIST")) {
    aurostd::string2tokens(vpflow.getattachedscheme("COMPARE2DATABASE::PROPERTY_LIST"),property_list,",");
    message << "OPTIONS: Extracting the following properties: " << aurostd::joinWDelimiter(property_list,", ");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vmatchbook.insert(vmatchbook.end(), property_list.begin(), property_list.end());
  }

  // ---------------------------------------------------------------------------
  // FLAG: get relaxation step
  uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_;
  if(vpflow.flag("COMPARE2DATABASE::RELAXATION_STEP")) {
    if(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")).find("orig") != std::string::npos ||
        vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP") == "0"){
      relaxation_step = _COMPARE_DATABASE_GEOMETRY_ORIGINAL_;
    }
    else if(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")).find("relax1") != std::string::npos ||
        aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")).find("middle_relax") != std::string::npos ||
        vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP") == "1"){
      relaxation_step = _COMPARE_DATABASE_GEOMETRY_RELAX1_;
    }
    message << "OPTIONS: Relaxation step (0=original, 1=relax1, 2=most_relaxed): " << relaxation_step << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // should not grab properties and compare structures other than the most relaxed structure
  if(property_list.size() && relaxation_step != _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){
    string relaxation_name = "";
    if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_){ relaxation_name = "original"; }
    else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_){ relaxation_name = "relax1"; }
    else { throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unexpected relaxation step input: "+aurostd::utype2string<uint>(relaxation_step)+".", _INPUT_ERROR_); }
    message << "The " << relaxation_name << " structures will be extracted; the properties may not correspond to these structures. Proceed with caution.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
  }

  // ---------------------------------------------------------------------------
  // get schema from xoptions, i.e., metadata (for the units and types)
  // DX20201230 - added type
  vector<string> property_units, property_types;
  string schema_unit = "", schema_type = "";
  for(uint p=0;p<property_list.size();p++){
    schema_unit = "SCHEMA::UNIT:" + aurostd::toupper(property_list[p]);
    schema_type = "SCHEMA::TYPE:" + aurostd::toupper(property_list[p]);
    property_units.push_back(XHOST.vschema.getattachedscheme(schema_unit));
    property_types.push_back(XHOST.vschema.getattachedscheme(schema_type));
  }

  // ---------------------------------------------------------------------------
  // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
  if(vpflow.flag("COMPARE2DATABASE::CATALOG")) {
    string catalog = aurostd::toupper(vpflow.getattachedscheme("COMPARE2DATABASE::CATALOG")); //DX20210615 - older versions of aflux require uppercase; safety
    if(catalog != "ALL"){ vmatchbook.push_back("catalog(\'" + catalog + "\')"); }
    // ME20220419 - Add MUTE operator to reduce response size
    if (!aurostd::WithinList(property_list, "catalog")) {
      vmatchbook.back() = "$" + vmatchbook.back();
    }
    message << "OPTIONS: Specify catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << " (default=all)" << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
  }

  // ---------------------------------------------------------------------------
  // AFLUX matchbook preparations: add number of species and species
  vmatchbook.push_back("nspecies(" + aurostd::utype2string<uint>(structure_containers[0].ntypes) + ")");
  // ME20220419 - Add MUTE operator to reduce response size
  if (!aurostd::WithinList(property_list, "nspecies")) {
    vmatchbook.back() = "$" + vmatchbook.back();
  }
  if(same_species){
    vmatchbook.push_back("species(" + aurostd::joinWDelimiter(structure_containers[0].elements,",") + ")");
    // ME20220419 - Add MUTE operator to reduce response size
    if (!aurostd::WithinList(property_list, "species")) {
      vmatchbook.back() = "$" + vmatchbook.back();
    }
  }

  // ---------------------------------------------------------------------------
  // AFLUX matchbook preparations: add space group symmetry
  bool ignore_symmetry = comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY");
  string sg_matchbook = "";
  if(!ignore_symmetry){
    if(!isSymmetryCalculated(structure_containers[0])){ calculateSymmetries(1); }  //1: one structure -> one processor

    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<structure_containers[0].grouped_Wyckoff_positions.size();i++){
        cerr << structure_containers[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // add space group query to AFLUX matchbook: get entries with compatible space groups, i.e., same or enantiomorph
    sg_matchbook = aflowlib::getSpaceGroupMatchbook(structure_containers[0].space_group,relaxation_step);
    // ME20220419 - Add MUTE operator to reduce response size
    vmatchbook.push_back(sg_matchbook);
    if (!aurostd::WithinList(property_list, sg_matchbook.substr(0, sg_matchbook.find("(")))) {
      vmatchbook.back() = "$" + vmatchbook.back();
    }
  }
  else { setSymmetryPlaceholders(); }

  // ---------------------------------------------------------------------------
  // AFLUX matchbook preparations: add aurl for entry
  //vmatchbook.push_back("aurl");  // ME20220419 - Part of the response already; slows down AFLUX
  // ME20220419 - Suppress unused default properties
  if ((sg_matchbook.find("spacegroup_relax") == string::npos) && !aurostd::WithinList(property_list, "spacegroup_relax")) {
    vmatchbook.push_back("$spacegroup_relax");
  }
  if (!aurostd::WithinList(property_list, "Pearson_symbol_relax")) {
    vmatchbook.push_back("$Pearson_symbol_relax");
  }

  // ---------------------------------------------------------------------------
  // construct aflux summons, i.e., combine matchbook
  string Summons = aurostd::joinWDelimiter(vmatchbook,",");
  message << "AFLUX matchbook request: " << Summons;
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  // ---------------------------------------------------------------------------
  // AFLUX matchbook preparations: format AFLUX output
  vmatchbook.push_back("format(aflow)");
  // ME20220420 - Use paged requests to avoid timeouts
  vmatchbook.push_back("");
  uint page = 1;
  string page_size_str = vpflow.getattachedscheme("COMPARE::PAGE_SIZE");
  uint page_size = page_size_str.empty()?DEFAULT_COMPARE_AFLUX_PAGE_SIZE:aurostd::string2utype<uint>(page_size_str);
  string response = "", response_page = "";
  do {
    vmatchbook.back() = "paging(" + aurostd::utype2string<uint>(page) + "," + aurostd::utype2string<uint>(page_size) + ")";
    Summons = aurostd::joinWDelimiter(vmatchbook, ",");

    // ---------------------------------------------------------------------------
    // call AFLUX
    response_page = aflowlib::AFLUXCall(Summons);
    response += response_page;
    page++;

  } while (!response_page.empty());

  vector<string> response_tokens;
  message << "Number of entries returned: " << aurostd::string2tokens(response,response_tokens,"\n");
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " AFLUX response:" << endl << response << endl;}

  // ---------------------------------------------------------------------------
  // extract properties from AFLUX response
  // (will switch over to CO's AFLUX+AFLOW code when available)
  vector<vector<std::pair<string,string> > > properties_response = aflowlib::getPropertiesFromAFLUXResponse(response);
  if(LDEBUG) {
    for(uint i=0;i<properties_response.size();i++){
      for(uint j=0;j<properties_response[i].size();j++){
        cerr << properties_response[i][j].first << " = " << properties_response[i][j].second << ", ";
      }
      cerr << endl;
    }
  }

  // ---------------------------------------------------------------------------
  // extract aurl, auid, and compound type from properties variable
  vector<string> auids, aurls, compounds;
  for(uint i=0;i<properties_response.size();i++){
    for(uint j=0;j<properties_response[i].size();j++){
      if(properties_response[i][j].first=="aurl"){
        aurls.push_back(properties_response[i][j].second);
      }
      if(properties_response[i][j].first=="auid"){
        auids.push_back(properties_response[i][j].second);
      }
      if(properties_response[i][j].first=="compound"){
        compounds.push_back(properties_response[i][j].second);
      }
    }
  }
  //cerr << "==============================" << endl;
  //::print(auids);
  //::print(aurls);
  //::print(compounds);

  message << "Loading structures ... ";
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  // ---------------------------------------------------------------------------
  // load and store entries from the database
  for(uint i=0; i<auids.size(); i++){

    // ---------------------------------------------------------------------------
    // first, get stoichiometry from entry
    vector<double> vcomposition;
    vector<string> species = aurostd::getElements(compounds[i], vcomposition);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " species=" << aurostd::joinWDelimiter(species,",") << endl;}
    vector<uint> tmp_stoich;
    // ME20220420 - Skip non-integer stoichiometry instead of throwing
    uint j = 0;
    for(;j<vcomposition.size();j++){
      if(aurostd::isinteger(vcomposition[j])){ tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j])); }
      else {
        //message << "Expected natoms in " << auids[i] << " to be an integer.";
        //throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
        break;
      }
    }
    if (j != vcomposition.size()) continue;
    vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
    if(!same_species){ std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end()); }

    // ---------------------------------------------------------------------------
    // second, check if stoichiometries are compatible
    // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
    // instead: filter on stoichiometry after recieving AFLUX response
    if(compare::sameStoichiometry(structure_containers[0].stoichiometry,tmp_reduced_stoich)){

      // ---------------------------------------------------------------------------
      // load structures from aflowlib_entry
      try {
        aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i];
        addDatabaseEntry2container(entry, species, relaxation_step, same_species);
      }
      // if cannot load one, keep going
      catch(aurostd::xerror& re){
        continue;
      }

      // store any properties
      for(uint l=0;l<properties_response[i].size();l++){
        bool property_requested = false;
        for(uint m=0;m<property_list.size();m++){
          if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
        }
        if(property_requested){
          structure_containers.back().properties.push_back(properties_response[i][l].second);
          structure_containers.back().properties_names = property_list;
          structure_containers.back().properties_units = property_units;
          structure_containers.back().properties_types = property_types;
        }
      }
    }
  }
  message << "Total number of candidate structures loaded: " << structure_containers.size(); //DX20190403
  pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_); //DX20190403

  return compareMultipleStructures(num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions
}

// ***************************************************************************
// compare::compare2database - compares database
// ***************************************************************************
namespace compare {
  // load input structure
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ostream& logstream){

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare2database [GENERAL_COMPARISON_OPTIONS] [COMPARE2DATABASE_OPTIONS] < file";
      string options_function_string = "compare2database_options: [--catalog=lib1|lib2|lib3|lib4|lib6|lib7|icsd] [--properties=enthalpy_atom,natoms,...] [--relaxation_step=original|relax1|most_relaxed]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_XTALFINDER_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::printCompare2Database()",options);
    }

    xstructure xstr(input,IOAFLOW_AUTO);
    ofstream FileMESSAGE;
    return printCompare2Database(xstr,vpflow,FileMESSAGE,logstream);
  }  //DX20200225
}

namespace compare {
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){
    xstructure xstr(input,IOAFLOW_AUTO);
    return printCompare2Database(xstr,vpflow,FileMESSAGE,logstream);
  }  //CO20200225
}

namespace compare {
  string printCompare2Database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
    stringstream message;
    ostringstream oss;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}  //CO20200508

    // ---------------------------------------------------------------------------
    // main compare2database() function
    XtalFinderCalculator xtal_finder_database(DEFAULT_XTALFINDER_MISFIT_MATCH,DEFAULT_XTALFINDER_MISFIT_FAMILY,FileMESSAGE,1,logstream);
    vector<StructurePrototype> prototypes_final = xtal_finder_database.compare2database(xstrIN, vpflow);

    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(prototypes_final.size()==0){
      return oss.str();
    }

    // ---------------------------------------------------------------------------
    // FLAG: print format (default is to write both), set the relevant bools
    // to false when another is specified
    bool write_txt = false;
    bool write_json = false;
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){ write_txt = true; }
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){ write_json = true; }
    // if not specified, write both by default
    if(!write_txt && !write_json){ write_txt = true; write_json = true; }

    // ---------------------------------------------------------------------------
    // FLAG: print2screen
    bool screen_only = false;
    if(vpflow.flag("COMPARE::SCREEN_ONLY")) {
      screen_only=true;
    }

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species = true;
    if(vpflow.flag("COMPARE2DATABASE::STRUCTURE")) {
      same_species = false;
    }

    // DEBUG oss << ss_out.str();
    message << "Number of structures in database matching with the input structure: " << prototypes_final[0].structures_duplicate.size() << "." << endl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // print results
    string results_txt = xtal_finder_database.printResults(prototypes_final, same_species, txt_ft);
    string results_json = xtal_finder_database.printResults(prototypes_final, same_species, json_ft);

    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(write_json){ return results_json; }
      if(write_txt){ return results_txt; }
    }
    // ---------------------------------------------------------------------------
    // write results to files //DX20201229 - consolidated into functions
    else{
      string directory = aurostd::getPWD();
      if(write_json){ xtal_finder_database.writeComparisonOutputFile(results_json, directory, json_ft, compare2database_xf, same_species); }
      if(write_txt){ xtal_finder_database.writeComparisonOutputFile(results_txt, directory, txt_ft, compare2database_xf, same_species); }
    }

    return oss.str();

  }
}

//DX - COMPARE DATABASE ENTRIES - START
// ***************************************************************************
// compare::compareDatabaseEntries - compares database entries
// ***************************************************************************
namespace compare {
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ostream& logstream){ //DX20191125 - added ofstream overload and added ostream as input
    ofstream FileMESSAGE;
    return compareDatabaseEntries(vpflow, FileMESSAGE, logstream);
  }

  string compareDatabaseEntries(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){ //DX20191125 - added ofstream and ostream
    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    string directory = aurostd::getPWD(); //DX20201229
    stringstream message;
    stringstream oss;

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare_database_entries [GENERAL_COMPARISON_OPTIONS] [COMPARE_DATABASE_ENTRIES_OPTIONS] < file";
      string options_function_string = "compare_database_entries_options: [--alloy=AgAlMn...] [--nspecies=3] [--catalog=lib1|lib2|lib3|lib4|lib6|lib7|icsd] [--properties=enthalpy_atom,natoms,...] [--relaxation_step=original|relax1|most_relaxed] [--space_group=225,186,227,...] [--stoichiometry=1,2,3,...]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_XTALFINDER_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::compareDatabaseEntries()",options);
    }

    // ---------------------------------------------------------------------------
    // FLAG: print format (default is to write both), set the relevant bools
    // to false when another is specified
    bool write_txt = false;
    bool write_json = false;
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){ write_txt = true; }
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){ write_json = true; }
    // if not specified, write both by default
    if(!write_txt && !write_json){ write_txt = true; write_json = true; }

    // ---------------------------------------------------------------------------
    // FLAG: print2screen
    bool screen_only = false;
    if(vpflow.flag("COMPARE::SCREEN_ONLY")) {
      screen_only=true;
    }

    vector<string> vmatchbook; //aflux - filter/get properties

    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator
    uint ncpus_default = 1; // initialize with one cpu; default
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        ncpus_default,
        logstream);

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
    xtal_finder.getOptions(vpflow, comparison_options);

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species = true;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::STRUCTURE")) {
      same_species = false;
      message << "OPTIONS: Structure-type comparison, i.e., ignore atomic species.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison //DX20191212
    bool magnetic_comparison = vpflow.flag("COMPARE::MAGNETIC");
    vector<string> magmoms_for_systems;
    if(magnetic_comparison){
      message << "OPTIONS: Including magnetic moment information in comparisons.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: property list to extract from database (using AFLUX)
    vector<string> property_list;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST")) {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST"),property_list,",");
      message << "OPTIONS: Extracting the following properties: " << aurostd::joinWDelimiter(property_list,", ");
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      vmatchbook.insert(vmatchbook.end(), property_list.begin(), property_list.end());
    }

    // ---------------------------------------------------------------------------
    // get schema from xoptions, i.e., metadata (for the units and types)
    // DX20201230 - added type
    vector<string> property_units, property_types;
    string schema_unit = "", schema_type = "";
    for(uint p=0;p<property_list.size();p++){
      schema_unit = "SCHEMA::UNIT:" + aurostd::toupper(property_list[p]);
      schema_type = "SCHEMA::TYPE:" + aurostd::toupper(property_list[p]);
      property_units.push_back(XHOST.vschema.getattachedscheme(schema_unit));
      property_types.push_back(XHOST.vschema.getattachedscheme(schema_type));
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the relaxation step to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")) {
      if(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")).find("orig") != std::string::npos ||
          vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP") == "0" ){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_ORIGINAL_;
      }
      else if(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")).find("relax1") != std::string::npos ||
          aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")).find("middle_relax") != std::string::npos ||
          vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP") == "1"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_RELAX1_;
      }
      message << "OPTIONS: Relaxation step (0=original, 1=relax1, 2=most_relaxed): " << relaxation_step << endl;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::CATALOG")) {
      catalog = aurostd::toupper(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::CATALOG")); //DX20190718 //DX20210615 - older versions of aflux require uppercase; safety
      if(catalog != "ALL"){
        vmatchbook.push_back("catalog(\'" + catalog + "\')");
        // ME20220419 - Add MUTE operator to reduce response size
        if (!aurostd::WithinList(property_list, "catalog")) {
          vmatchbook.back() = "$" + vmatchbook.back();
        }
      }
      message << "OPTIONS: Specify catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << " (default=all)" << endl;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(catalog=="" || catalog=="icsd" || catalog=="all"){ //DX20191108 - needs to be outside of loop
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",TRUE);
    }
    // ---------------------------------------------------------------------------
    // FLAG: specify number of species
    uint arity=0; //Defalut=0 : all
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ARITY")) {
      arity=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::ARITY"));
      message << "OPTIONS: Getting entries with nspecies=" << arity;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      if(arity!=0){ vmatchbook.push_back("nspecies(" + aurostd::utype2string<uint>(arity) + ")"); }
      // ME20220419 - Add MUTE operator to reduce response size
      if (!aurostd::WithinList(property_list, "nspecies")) {
        vmatchbook.back() = "$" + vmatchbook.back();
      }
    }

    // ---------------------------------------------------------------------------
    // FLAG : stoichiometry
    // Note: not added to vmatchbook, checked when entries are returned (more robust)
    vector<uint> stoichiometry_reduced, stoichiometry_input;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::STOICHIOMETRY")){
      vector<string> stoichiometry_vstring;
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::STOICHIOMETRY"), stoichiometry_vstring,",");
      for(uint i=0;i<stoichiometry_vstring.size();i++){ stoichiometry_input.push_back(aurostd::string2utype<uint>(stoichiometry_vstring[i])); }
      // check input for consistency with arity (if given)
      if(arity!=0 && arity!=stoichiometry_input.size()){
        message << "arity=" << arity << " and stoichiometry=" << aurostd::joinWDelimiter(stoichiometry_input,",") << " do not match.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ERROR_);
      }
      aurostd::reduceByGCD(stoichiometry_input, stoichiometry_reduced);
      // if structure comparison, sort stoichiometry, otherwise perserve order
      if(!same_species){ std::sort(stoichiometry_reduced.begin(),stoichiometry_reduced.end()); }
      message << "OPTIONS: Getting entries with stoichiometry=" << aurostd::joinWDelimiter(stoichiometry_reduced,",");
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify alloy systems
    vector<string> species;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ALLOY")){
      string alloy_string = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::ALLOY");
      // split by comma
      if(alloy_string.find(",") != std::string::npos){
        aurostd::string2tokens(alloy_string,species,",");
      }
      // split by colon
      else if(alloy_string.find(":") != std::string::npos){
        aurostd::string2tokens(alloy_string,species,":");
      }
      // split by alloy species (no delimiter)
      else{
        species = aurostd::getElements(alloy_string); //DX20191106
      }
      if(species.size()!=0){
        vmatchbook.push_back("species(" + aurostd::joinWDelimiter(species,",") + ")");
        // ME20220419 - Add MUTE operator to reduce response size
        if (!aurostd::WithinList(property_list, "species")) {
          vmatchbook.back() = "$" + vmatchbook.back();
        }
      }
    }

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get space group(s) //DX20200929
    string sg_matchbook = xtal_finder.getSpaceGroupMatchbookFromOptions(vpflow, relaxation_step);
    if (!sg_matchbook.empty()) {
      vmatchbook.push_back(sg_matchbook);
      // ME20220419 - Add MUTE operator to reduce response size
      if (!aurostd::WithinList(property_list, sg_matchbook.substr(0, sg_matchbook.find("(")))) {
        vmatchbook.back() = "$" + vmatchbook.back();
      }
    }

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get aurl for entry
    //vmatchbook.push_back("aurl");  // ME20220419 - Part of the response already; slows down AFLUX
    // ME20220419 - Suppress unused default properties
    if ((sg_matchbook.find("spacegroup_relax") == string::npos) && !aurostd::WithinList(property_list, "spacegroup_relax")) {
      vmatchbook.push_back("$spacegroup_relax");
    }
    if (!aurostd::WithinList(property_list, "Pearson_symbol_relax")) {
      vmatchbook.push_back("$Pearson_symbol_relax");
    }

    // ---------------------------------------------------------------------------
    // construct aflux summons, i.e., combine matchbook
    string Summons = aurostd::joinWDelimiter(vmatchbook,",");
    message << "AFLUX matchbook request: " << Summons;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: format AFLUX output
    vmatchbook.push_back("format(aflow)");
    // ME20220420 - Use paged requests to avoid timeouts
    vmatchbook.push_back("");
    uint page = 1;
    string page_size_str = vpflow.getattachedscheme("COMPARE::PAGE_SIZE");
    uint page_size = page_size_str.empty()?DEFAULT_COMPARE_AFLUX_PAGE_SIZE:aurostd::string2utype<uint>(page_size_str);
    string response = "", response_page = "";
    do {
      vmatchbook.back() = "paging(" + aurostd::utype2string<uint>(page) + "," + aurostd::utype2string<uint>(page_size) + ")";
      Summons = aurostd::joinWDelimiter(vmatchbook, ",");

      // ---------------------------------------------------------------------------
      // call AFLUX
      response_page = aflowlib::AFLUXCall(Summons);
      response += response_page;
      page++;

    } while (!response_page.empty());

    if(LDEBUG){cerr << __AFLOW_FUNC__ << "::AFLUX response:" << endl << response << endl;}

    vector<string> response_tokens;
    message << "Number of entries returned: " << aurostd::string2tokens(response,response_tokens,"\n");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // extract properties from AFLUX response
    vector<vector<std::pair<string,string> > > properties_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    if(LDEBUG){
      for(uint i=0;i<properties_response.size();i++){
        for(uint j=0;j<properties_response[i].size();j++){
          cerr << properties_response[i][j].first << " = " << properties_response[i][j].second << ", ";
        }
        cerr << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // extract aurl, auid, and compound type from properties variable
    vector<string> auids, aurls, compounds;
    for(uint i=0;i<properties_response.size();i++){
      for(uint j=0;j<properties_response[i].size();j++){
        if(properties_response[i][j].first=="aurl"){
          aurls.push_back(properties_response[i][j].second);
        }
        if(properties_response[i][j].first=="auid"){
          auids.push_back(properties_response[i][j].second);
        }
        if(properties_response[i][j].first=="compound"){
          compounds.push_back(properties_response[i][j].second);
        }
      }
    }
    //cerr << "==============================" << endl;
    //::print(auids);
    //::print(aurls);
    //::print(compounds);

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ..." << auids.size();
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // load and store entries from the database
    for(uint i=0; i<auids.size(); i++){

      // ---------------------------------------------------------------------------
      // first, get stoichiometry from entry
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " species=" << aurostd::joinWDelimiter(species,",") << endl;}

      vector<uint> tmp_stoich;
      // ME20220420 - Skip non-integer stoichiometry instead of throwing
      uint j = 0;
      for(;j<vcomposition.size();j++){
        if(aurostd::isinteger(vcomposition[j])){ tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j])); }
        else {
          //message << "Expected natoms in " << auids[i] << " to be an integer.";
          //throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__,message,_RUNTIME_ERROR_);
          break;
        }
      }
      if (j != vcomposition.size()) continue;
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      if(!same_species){ std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end()); }

      // ---------------------------------------------------------------------------
      // second, check if stoichiometries are compatible (if given)
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      if(stoichiometry_reduced.size()==0 || compare::sameStoichiometry(stoichiometry_reduced,tmp_reduced_stoich)){

        // ---------------------------------------------------------------------------
        // load structures from aflowlib_entry
        try {
          aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i];
          xtal_finder.addDatabaseEntry2container(entry, species, relaxation_step, same_species);
        }
        // if structure cannot be loaded, keep going
        catch(aurostd::xerror& re){
          continue;
        }

        // store any properties
        for(uint l=0;l<properties_response[i].size();l++){
          bool property_requested = false;
          for(uint m=0;m<property_list.size();m++){
            if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
          }
          if(property_requested){
            xtal_finder.structure_containers.back().properties.push_back(properties_response[i][l].second);
            xtal_finder.structure_containers.back().properties_names = property_list;
            xtal_finder.structure_containers.back().properties_units = property_units;
            xtal_finder.structure_containers.back().properties_types = property_types;
          }
        }
      }
    }


    message << "Total number of candidate structures loaded: " << xtal_finder.structure_containers.size(); //DX20190403
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_); //DX20190403

    vector<StructurePrototype> prototypes_final = xtal_finder.compareMultipleStructures(xtal_finder.num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // print results
    string results_txt = xtal_finder.printResults(prototypes_final, same_species, txt_ft);
    string results_json = xtal_finder.printResults(prototypes_final, same_species, json_ft);

    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only) {
      if(write_json){ return results_json; }
      if(write_txt){ return results_txt; }
    }
    // ---------------------------------------------------------------------------
    // write results to files //DX20201229 - consolidated into functions
    else{
      if(write_txt){ xtal_finder.writeComparisonOutputFile(results_txt, directory, txt_ft, compare_database_entries_xf, same_species); }
      if(write_json){ xtal_finder.writeComparisonOutputFile(results_json, directory, json_ft, compare_database_entries_xf, same_species); }
    }

    return oss.str();

  }
}
//DX - COMPARE DATABASE ENTRIES - END

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromStructureList() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromStructureList(
    const vector<string>& filenames,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // directory to write results
  string directory = aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // load structures appended to command
  loadStructuresFromStructureList(
      filenames,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromDirectory() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromDirectory(
    const string& directory,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // load structures in directory
  loadStructuresFromDirectory(
      directory,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromFile() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromFile(
    const string& filename,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // directory to write results
  string directory = aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // load structures in file
  loadStructuresFromFile(filename,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromString() //ME20201206
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromString(
    const string& structures_string,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // directory to write results
  string directory = aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // load structures in string
  stringstream structures;
  structures << structures_string;
  loadStructuresFromStringstream(structures,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info

  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// compare::getUniqueEntries() //DX20201111
// ***************************************************************************
// Get structurally unique aflowlib entries. Helper function to GFA-code.
// If needed, this function can be extended to read in some of the
// structural properties from the database (i.e., symmetry) and suppress the
// on-the-fly analysis (potential speed-up).
namespace compare {
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries, uint num_proc, bool same_species, bool scale_volume, bool optimize_match){
    ostringstream oss;
    ofstream FileMESSAGE;
    return getUniqueEntries(entries, oss, FileMESSAGE, num_proc, same_species, scale_volume, optimize_match);
  }
}

namespace compare {
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries,
      ostream& oss,
      ofstream& FileMESSAGE,
      uint num_proc,
      bool same_species,
      bool scale_volume,
      bool optimize_match){

    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        num_proc,
        oss);

    // ---------------------------------------------------------------------------
    // load structures from aflowlib entries
    vector<string> magmoms_for_systems; //DX20201111 - not included for now
    xtal_finder.loadStructuresFromAflowlibEntries(
        entries,
        magmoms_for_systems,
        same_species);

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = aurostd::getPWD();

    // ---------------------------------------------------------------------------
    // default comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
    comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",scale_volume);
    comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",optimize_match);

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    vector<StructurePrototype> grouped_structures = xtal_finder.compareMultipleStructures(
        xtal_finder.num_proc,
        same_species,
        directory,
        comparison_options);

    // ---------------------------------------------------------------------------
    // filter out duplicate aflowlib entries
    vector<aflowlib::_aflowlib_entry> entries_unique;
    for(uint i=0;i<grouped_structures.size();i++){
      for(uint j=0;j<entries.size();j++){
        if(grouped_structures[i].structure_representative->name==entries[j].auid){
          entries_unique.push_back(entries[j]);
          break;
        }
      }
    }

    return entries_unique;
  }
}

// ***************************************************************************
// XtalFinderCalculator::groupSimilarXstructures() //DX20201229
// ***************************************************************************
/// @brief Compares a vector of xstructures and groups the indices of the
/// xstructures by structural similarity
///
/// @param vxstrs vector of xstructures
/// @param same_species require the same atom decorations/types
/// @param comparison_options options to be used for the comparison
///
/// @return vector of a vector of indicies of the grouped structures
///
/// @authors
/// @mod{DX,20201229,created function}
/// @mod{SD,20221117,require explicit comparison options input}
vector<vector<uint> > XtalFinderCalculator::groupSimilarXstructures(
    const vector<xstructure>& vxstrs,
    aurostd::xoption& comparison_options,
    bool same_species) {

  // Compares a vector of xstructures and groups the indices of the
  // xstructures by structural similarity

  uint nstrs= vxstrs.size();

  stringstream structure_name;
  string source = "input";
  uint relaxation_step = 0; // unknowable from input
  string directory = aurostd::getPWD();

  // add structures to containers
  for(uint i=0;i<nstrs;i++){
    structure_name.str("");
    structure_name << i; // name by index
    addStructure2container(vxstrs[i], structure_name.str(), source, relaxation_step, same_species);
  }

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  vector<StructurePrototype> grouped_structures = compareMultipleStructures(
      num_proc,
      same_species,
      directory,
      comparison_options); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // group indices based on structural similarity
  vector<vector<uint> > grouped_indices;
  vector<uint> structure_indices_equivalent;
  for(uint i=0;i<grouped_structures.size();i++){
    structure_indices_equivalent.clear();
    structure_indices_equivalent.push_back(aurostd::string2utype<uint>(grouped_structures[i].structure_representative->name));
    for(uint j=0;j<grouped_structures[i].structures_duplicate.size();j++){
      structure_indices_equivalent.push_back(aurostd::string2utype<uint>(grouped_structures[i].structures_duplicate[j]->name));
    }
    grouped_indices.push_back(structure_indices_equivalent);
  }

  return grouped_indices;
}

// ***************************************************************************
// XtalFinderCalculator::groupSimilarXstructures()
// ***************************************************************************
/// @brief Compares a vector of xstructures and groups the indices of the
/// xstructures by structural similarity
///
/// @param vxstrs vector of xstructures
/// @param same_species require the same atom decorations/types
///
/// @return vector of a vector of indicies of the grouped structures
///
/// @authors
/// @mod{SD,20221117,created function}
vector<vector<uint> > XtalFinderCalculator::groupSimilarXstructures(
    const vector<xstructure>& vxstrs,
    bool same_species) {
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
  return groupSimilarXstructures(vxstrs, comparison_options, same_species);
}

// ***************************************************************************
// compare::structuresMatch() //DX20201228 - renamed
// ***************************************************************************
namespace compare {
  bool structuresMatch(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      uint num_proc) {

    double misfit_final = AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    return aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, misfit_final, num_proc); //DX20191122 - move ostream to end
  }
}

namespace compare {
  bool structuresMatch(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      bool scale_volume,
      bool optimize_match,
      uint num_proc) {

    double misfit_final = AUROSTD_MAX_DOUBLE;
    return aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, misfit_final, num_proc); //DX20191122 - move ostream to end and add default
  }
}


// ***************************************************************************
// compare::getMisfitBetweenStructures() //DX20201228 - renamed
// ***************************************************************************
namespace compare {
  double getMisfitBetweenStructures(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      uint num_proc) {

    double misfit_final=AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, misfit_final, num_proc); //DX20191122 - move ostream to end
    return misfit_final;
  }
}

// ***************************************************************************
// compare::getTransformationBetweenStructures() //DX20201229
// ***************************************************************************
namespace compare {
  structure_mapping_info getTransformationBetweenStructures(const xstructure& xstr1, const xstructure& xstr2, bool same_species, uint num_proc) { //DX20191108 - remove const & from bools
    double misfit_final=AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    structure_mapping_info misfit_info_final = compare::initialize_misfit_struct();
    aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, misfit_final, misfit_info_final, num_proc); //DX20191122 - move ostream to end
    return misfit_info_final;
  }
}

// ***************************************************************************
// compare::aflowCompareStructure
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      bool scale_volume,
      bool optimize_match,
      double& misfit_final,
      uint num_proc) {

    structure_mapping_info misfit_info_final = compare::initialize_misfit_struct(); //DX20191218

    return aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, misfit_final, misfit_info_final, num_proc); //DX20191122 - move ostream to end
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      bool scale_volume,
      bool optimize_match,
      double& misfit_final,
      structure_mapping_info& misfit_info_final,
      uint num_proc) { //DX20191108 - remove const & from bools //DX20191122 - move ostream to end and add default

    structure_container str_rep = compare::initializeStructureContainer(xstr1, same_species);
    structure_container str_matched = compare::initializeStructureContainer(xstr2, same_species);

    XtalFinderCalculator xtal_finder(num_proc);
    xtal_finder.compareStructures(str_rep,str_matched,misfit_info_final,same_species,scale_volume,optimize_match);

    misfit_final = misfit_info_final.misfit;
    return(misfit_final<=xtal_finder.misfit_match);

  }
}

// AFLOW-XtalFinder (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
