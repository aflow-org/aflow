// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2021              *
// *                                                                         *
// ***************************************************************************
// Written by Rico Friedrich, Corey Oses, and Marco Esters
// rico.friedrich@duke.edu

#ifndef _AFLOW_CCE_CPP_
#define _AFLOW_CCE_CPP_

#include "aflow.h"
#include "aflowlib.h"
#include "aflow_pflow.h"
#include "aflow_cce.h"

using std::cout;
using std::cerr;
using std::endl;

#define CCE_DEBUG false
#define OX_NUMS_FROM_EN_ALLEN 1
#define OX_NUMS_FROM_BADER_CHARGES 2

namespace cce {

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           MAIN CCE FUNCTIONS                            //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // for command line use, 
  // use inside AFLOW providing directory path or xstructure & functional string or flags and istream for web tool, 
  // and CCE core function called by all other main CCE functions

  //CO20201105 START
  void run(aurostd::xoption& flags, ostream& oss) { //CO20201105
    if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "PYTHON") {
      string directory=aurostd::getPWD(); //CO20201126 //"."; //can change later with a flag input
      string aflow_cce_python_subdir = "AFLOW_CCE_PYTHON";
      string python_directory=directory + '/' + aflow_cce_python_subdir;
      aurostd::DirectoryMake(python_directory);
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Writing out python script to: "+python_directory, oss, _LOGGER_NOTICE_);
      stringstream output;
      aurostd::EmbData::save_to_file("aflow_cce_python.py", python_directory+'/'+"aflow_cce_python.py");
      // print CCE citation
      oss << print_citation();
    } else {
      if(flags.flag("CCE_CORRECTION")){cce::print_corrections(flags,oss);}
    }
  }
  void run(aurostd::xoption& flags, std::istream& ist, ostream& oss) {  //CO20201105
    if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "PYTHON") {
      string directory=aurostd::getPWD(); //CO20201126 //"."; //can change later with a flag input
      string aflow_cce_python_subdir = "AFLOW_CCE_PYTHON";
      string python_directory=directory + '/' + aflow_cce_python_subdir;
      aurostd::DirectoryMake(python_directory);
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Writing out python script to: "+python_directory, oss, _LOGGER_NOTICE_);
      aurostd::EmbData::save_to_file("aflow_cce_python.py", python_directory+'/'+"aflow_cce_python.py");
      // print CCE citation
      oss << print_citation();
    } else {
      if(flags.flag("CCE_CORRECTION::POSCAR2CCE")) {cce::print_corrections(flags, ist, oss);}
      if(flags.flag("CCE_CORRECTION::GET_CCE_CORRECTION")) {cce::print_corrections(flags, ist, oss);}
      if(flags.flag("CCE_CORRECTION::GET_OXIDATION_NUMBERS")) {cce::print_oxidation_numbers(flags, ist, oss);}
      if(flags.flag("CCE_CORRECTION::GET_CATION_COORDINATION_NUMBERS")) {cce::print_cation_coordination_numbers(flags, ist, oss);}
    }
  }
  //CO20201105 END

  //print_corrections////////////////////////////////////////////////////////
  // main CCE function for command line use 
  // for reading input, analyzing structure, determining oxidation numbers, assigning corrections, 
  // calculating total corrections and corrected formation enthalpies, and writing output
  void print_corrections(aurostd::xoption& flags, ostream& oss) {

    /************************************/
    // Print user instructions
    /************************************/
    // option to print user instructions and exit upon completion
    if(flags.flag("CCE_CORRECTION::USAGE")){
      oss << print_usage();
      return;
    }

    /************************************/
    // Read structure
    /************************************/
    //read structural data from structure file provided on command line
    xstructure structure=read_structure(flags.getattachedscheme("CCE_CORRECTION::POSCAR_PATH"));

    print_corrections(structure, flags);
  }

  void print_corrections(xstructure& structure, aurostd::xoption& flags, ostream& oss) {  //CO20201105
    aurostd::xoption cce_flags = init_flags();
    if (aurostd::toupper(flags.flag("CCE_CORRECTION::UNIT_TEST"))) {
      cce_flags.flag("UNIT_TEST",TRUE);
    }
    CCE_Variables cce_vars = init_variables(structure);

    /********************************************************/
    // Read DFT formation enthalpies and functionals if provided
    /********************************************************/
    get_dft_form_enthalpies_functionals(flags.getattachedscheme("CCE_CORRECTION::ENTHALPIES_FORMATION_DFT"), flags.getattachedscheme("CCE_CORRECTION::FUNCTIONALS"), cce_vars); //provide precalc. DFT formation enthalpies & corresponding functionals

    /********************************************************/
    // Read oxidation numbers if provided
    /********************************************************/
    if(flags.flag("CCE_CORRECTION::OXIDATION_NUMBERS")){
      cce_flags.flag("OX_NUMS_PROVIDED",TRUE);
      cce_vars.oxidation_states = get_oxidation_states(flags.getattachedscheme("CCE_CORRECTION::OXIDATION_NUMBERS"), structure, cce_vars);
    }

    print_corrections(structure, flags, cce_flags, cce_vars, oss);  //CO20201105
  }

  void print_corrections(xstructure& structure, aurostd::xoption& flags, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    /********************************************************/
    // Get CCE corrections
    /********************************************************/
    // obtain total CCE corrections per cell from CCE core function only based on
    // structure, oxidation numbers and functional information
    vector<vector<uint> > multi_anion_num_neighbors;
    CCE_core(structure, cce_flags, cce_vars, multi_anion_num_neighbors);

    // CALCULATE CORRECTED FORMATION ENTHALPIES AT 298.15 AND 0K ###################################
    if (cce_flags.flag("CORRECTABLE")){
      //calculate CCE corrected DFT formation enthalpies if precalculated DFT formation enthalpies are provided
      uint num_funcs=cce_vars.vfunctionals.size();
      uint num_temps=cce_vars.vtemperatures.size();
      for(uint k=0,ksize=num_funcs;k<ksize;k++){ // looping over and checking of vfunctionals is necessary to ensure correct correspondence between given formation enthalpies [k] and corrections with respect to the functional
        if (cce_vars.vfunctionals[k] != "exp") {
          if(cce_vars.enthalpies_dft.size()!=0){ // only when precalculated DFT formation enthalpies are provided, calculate corrected values
            for (uint l = 0; l < num_temps; l++) {
              cce_vars.enthalpy_formation_cell_cce[num_temps*k+l] = cce_vars.enthalpies_dft[k] - cce_vars.cce_correction[num_temps*k+l] ; // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
            }
          }
        } else { // for CCE@exp the formation enthalpy (only at 298.15 K) is directly calculated from the exp. formation enthalpies per bond times the number of cation-anion bonds and is not subtracted from a given value
          cce_vars.enthalpy_formation_cell_cce[num_temps*k] = cce_vars.cce_correction[num_temps*k] ;
        }
      }

      //// print output
      oss << print_output(structure, flags, cce_flags, cce_vars, multi_anion_num_neighbors);
    }
  } // main CCE function for command line use



  //print_corrections///////////////////////////////////////////////////////////////////////
  //ME20200213
  // For poscar2cce
  void print_corrections(aurostd::xoption& flags, std::istream& ist, ostream& oss) {  //CO20201105
    // read structure
    xstructure structure=read_structure(ist);

    print_corrections(structure, flags, oss);
  }



  //print_cation_coordination_numbers///////////////////////////////////////////////////////////////////////
  // For determining the the cation coordination numbers, i.e. the number of anion neighbors for each cation
  void print_cation_coordination_numbers(aurostd::xoption& flags, std::istream& ist, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    // read structure
    xstructure structure=read_structure(ist);

    // initialise flags and variables
    aurostd::xoption cce_flags = init_flags();
    CCE_Variables cce_vars = init_variables(structure);

    // determine which species is the anion (for single anion systems)
    cce_vars.anion_species=determine_anion_species(structure, cce_vars);

    // check whether it is a multi-anion system and set variables accordingly
    cce_vars.multi_anion_atoms=check_for_multi_anion_system(structure, cce_flags, cce_vars);

    // tolerance for counting anion neighbors should be adjusted if requested by user
    double tolerance=_CCE_NN_DIST_TOL_;
    if(!flags.getattachedscheme("CCE_CORRECTION::DIST_TOL").empty()){
      tolerance=aurostd::string2utype<double>(flags.getattachedscheme("CCE_CORRECTION::DIST_TOL"));
    }

    // apply species selective cutoffs to determine only nearest neighbors within respective cutoff
    cce_vars.num_neighbors=get_num_neighbors(structure, cce_vars.anion_species, cce_flags, cce_vars, tolerance);

    // determine anion nearest neighbors for cations bound to multi anion atoms if needed
    vector<vector<uint> > multi_anion_num_neighbors(cce_vars.multi_anion_species.size(), vector<uint>(structure.atoms.size()));
    if (cce_flags.flag("MULTI_ANION_SYSTEM")){
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " getting neighbors for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
        }
        multi_anion_num_neighbors[k]=get_num_neighbors(structure, cce_vars.multi_anion_species[k], cce_flags, cce_vars, tolerance);
      }
    }

    // check for per- and superoxide ions based on O-O bond length and set number of (su-)peroxide bonds 
    // and indices accordingly if ions are detected
    // added here to throw warning when (su-)peroxide is detected
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      check_per_super_oxides(structure, cce_flags, cce_vars);
    }

    if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") {
      oss << print_JSON_cation_coordination_numbers(structure, cce_flags, cce_vars, multi_anion_num_neighbors) << std::endl;
    } else {
      // print cation coordination numbers
      oss << print_output_cation_coordination_numbers(structure, cce_flags, cce_vars, multi_anion_num_neighbors);
      // print CCE citation
      oss << print_citation();
    }
  }



  //print_oxidation_numbers///////////////////////////////////////////////////////////////////////
  // For determining only oxidation numbers of the system (from Allen electronegativities) without corrections
  void print_oxidation_numbers(aurostd::xoption& flags, std::istream& ist, ostream& oss) {
    // read structure
    xstructure structure=read_structure(ist);

    // initialise flags and variables
    aurostd::xoption cce_flags = init_flags();
    CCE_Variables cce_vars = init_variables(structure);

    // determine oxidation numbers
    cce_vars.oxidation_states=get_oxidation_states_from_electronegativities(structure);

    if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") {
      oss << print_JSON_ox_nums(structure, cce_vars) << std::endl;
    } else {
      // print oxidation numbers
      oss << print_output_oxidation_numbers(structure, cce_vars);
      // print CCE citation
      oss << print_citation();
    }
  }



  //calculate_corrections////////////////////////////////////////////////////////
  // main CCE function for calling inside AFLOW providing only directory path where data neded for correction (POSCAR.static & aflow.in) are located
  // for setting parameters, analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, converting correction vector, and returning corrections
  vector<double> calculate_corrections(const string& directory_path) {
    stringstream message;
    // get structure
    string poscar="";
    try { // try most relaxed CONTCAR first
      poscar=aflowlib::vaspfile2stringstream(directory_path,"CONTCAR");
    } catch (aurostd::xerror& e) { // if that doesn't work try POSCAR
      poscar=aflowlib::vaspfile2stringstream(directory_path,"POSCAR");
    }
    xstructure structure=read_structure(poscar); // AFLOW seems to automatically unzip and rezip zipped files so that only the file name without zipping extension needs to be given
    // get functional from aflow.in in directory
    string aflowin_file= directory_path + "/" + _AFLOWIN_;
    string outcar_file= directory_path + "/" + "OUTCAR.relax1";
    string functional=get_functional_from_aflow_in_outcar(structure, aflowin_file, outcar_file);
    if (functional.empty()) {
      message << " Functional cannot be determined from aflow.in. Corrections are available for PBE, LDA, SCAN, or PBE+U:ICSD.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    return calculate_corrections(structure, functional, directory_path); // directory path must be propagated if ox. states are determined from Bader charges
  } // main CCE function for calling inside AFLOW with directory path



  //calculate_corrections////////////////////////////////////////////////////////
  // main CCE function for calling inside AFLOW providing only structure (and functional) assuming that oxidation numbers will be obtained from Bader charges or electronegativities
  // for setting parameters, analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, converting correction vector, and returning corrections
  //vector<double> CCE(xstructure& structure) // OLD: functional will be automatically determined during Bader charge analysis for the current implementation, later when using e.g. electronegativities, it might be needed as input
  vector<double> calculate_corrections(const xstructure& structure, string functional, const string& directory_path, ostream& oss) {ofstream FileMESSAGE;return calculate_corrections(structure, functional, FileMESSAGE, directory_path, oss);} // functional needed as input when determining oxidation numbers from electronegativities
  vector<double> calculate_corrections(const xstructure& structure, string functional, ofstream& FileMESSAGE, const string& directory_path, ostream& oss) { // functional needed as input when determining oxidation numbers from electronegativities
    stringstream message;
    // copy structure to structure_to_use since checkStructure includes rescaling to 1
    xstructure structure_to_use=structure;
    // check structure
    structure_to_use.checkStructure(); //includes rescaling the structure to 1
    // init variables and flags
    CCE_Variables cce_vars = init_variables(structure_to_use);
    aurostd::xoption cce_flags = init_flags();
    // determine functional
    vector<string> CCE_vallowed_functionals;
    aurostd::string2tokens(CCE_allowed_functionals, CCE_vallowed_functionals, ",");
    if(functional!="exp"){
      functional=aurostd::toupper(functional);
    }
    if (!aurostd::WithinList(CCE_vallowed_functionals, functional) || get_offset(functional) == -1) {
      message << " Unknown functional " << functional << ". Please choose PBE, LDA, or SCAN.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    cce_vars.vfunctionals.push_back(functional);
    cce_vars.offset.push_back(get_offset(functional));
    // run main CCE function to determine correction
    vector<vector<uint> > multi_anion_num_neighbors;
    CCE_core(structure_to_use, cce_flags, cce_vars, multi_anion_num_neighbors, directory_path);
    // print oxidation_numbers
    message << print_output_oxidation_numbers(structure_to_use, cce_vars);
    _aflags aflags;aflags.Directory=aurostd::getPWD();
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    // cce_vars.cce_corrections can be returned directly since there is always only one functional for this CCE function
    return cce_vars.cce_correction;
  } // main CCE function for calling inside AFLOW



  //CCE_core////////////////////////////////////////////////////////
  // main CCE function core called by all other main CCE functions
  // analyzing structure, determining oxidation numbers, assigning corrections,
  // and calculating total corrections
  void CCE_core(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, vector<vector<uint> >& multi_anion_num_neighbors, const string& directory_path) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    // if there is only one species, it must be an elemental phase and is hence not correctable
    if (structure.species.size() == 1){
      message << " BAD NEWS: Only one species found. Enthalpies of elemental systems cannot be corrected with the CCE methodology.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    // set flag that full correction scheme is run and not just determination of ox. nums. or num. anion neighbors
    cce_flags.flag("RUN_FULL_CCE",TRUE);
    // determine temperatures
    aurostd::string2tokens(CCE_temperatures, cce_vars.vtemperatures, ",");
    uint num_temps=cce_vars.vtemperatures.size();
    // corrections_atom, (su-)perox_correction, cce_correction and enthalpy_formation_cell_cce can only be resized after vfunctionals.size() is known from previous main CCE functions calling CCE_core
    // vfunctionals.size()*num_temps for corrections for different temperatures (currently 298.15 and 0K)
    cce_vars.corrections_atom.resize(cce_vars.vfunctionals.size()*num_temps, vector<double>(structure.atoms.size()));
    cce_vars.perox_correction.resize(cce_vars.vfunctionals.size()*num_temps);
    cce_vars.superox_correction.resize(cce_vars.vfunctionals.size()*num_temps);
    cce_vars.cce_correction.resize(cce_vars.vfunctionals.size()*num_temps);
    cce_vars.enthalpy_formation_cell_cce.resize(cce_vars.vfunctionals.size()*num_temps);

    // DETERMINE NUMBER OF NEAREST ANION NEIGHBORS FOR EACH CATION: STRUCTURAL PART OF CORRECTION ##
    /********************************************************/
    // Determine anion species
    /********************************************************/
    // determine which species is the anion (for single anion systems)
    cce_vars.anion_species=determine_anion_species(structure, cce_vars);

    /********************************************************/
    // Check for multi anion system
    /********************************************************/
    // check whether it is a multi-anion system and set variables accordingly
    cce_flags.flag("MULTI_ANION_SYSTEM",FALSE);
    cce_flags.flag("O_MULTI_ANION_SPECIES",FALSE); // whether one of the multi anion species is O for which per/superoxide tests need to be made
    cce_vars.multi_anion_atoms=check_for_multi_anion_system(structure, cce_flags, cce_vars);
    // multi anion corrections can only be resized after number of multi anion species is known from check for multi anion system
    cce_vars.multi_anion_corrections_atom.resize(cce_vars.multi_anion_species.size(), vector<vector<double> >(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size())));

    /********************************************************/
    // Determine anion nearest neighbors for each cation
    /********************************************************/
    // apply species selective cutoffs to determine only nearest neighbors within respective cutoff
    cce_vars.num_neighbors=get_num_neighbors(structure, cce_vars.anion_species, cce_flags, cce_vars);

    // determine anion nearest neighbors for cations bound to multi anion atoms if needed
    //vector<vector<uint> > multi_anion_num_neighbors(cce_vars.multi_anion_species.size(), vector<uint>(structure.atoms.size()));
    multi_anion_num_neighbors.resize(cce_vars.multi_anion_species.size(), vector<uint>(structure.atoms.size()));
    if (cce_flags.flag("MULTI_ANION_SYSTEM")){
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " getting neighbors for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
        }
        multi_anion_num_neighbors[k]=get_num_neighbors(structure, cce_vars.multi_anion_species[k], cce_flags, cce_vars);
      }
    }

    /********************************************************/
    // Check for per- and superoxide ions
    /********************************************************/
    // check for per- and superoxide ions based on O-O bond length and set number of (su-)peroxide bonds 
    // and indices accordingly if ions are detected
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      check_per_super_oxides(structure, cce_flags, cce_vars);
    }


    // DETERMINE CORRECTION FOR EACH CATION: OXIDATION STATE DEPENDENT PART OF COORECTION ##########
    /********************************************************/
    // Assign oxidation numbers
    /********************************************************/
    // determine oxidation numbers automatically from structure and Allen electronegativities if not provided on command line
    if(!cce_flags.flag("OX_NUMS_PROVIDED")) {
      if(DEFAULT_CCE_OX_METHOD == OX_NUMS_FROM_EN_ALLEN) { // determining oxidation numbers from Allen electronegativities is the default
        cce_vars.oxidation_states=get_oxidation_states_from_electronegativities(structure, cce_flags, cce_vars);
      } else if(DEFAULT_CCE_OX_METHOD == OX_NUMS_FROM_BADER_CHARGES) { // obtaining oxidation states from Bader charges is outdated but the functionality is kept mainly for test purposes
        cce_vars.oxidation_states=get_oxidation_states_from_Bader(structure, cce_flags, cce_vars, directory_path);
      }
    }

    /********************************************************/
    // Assign corrections
    /********************************************************/
    // assigning corrections for each atom after oxidation numbers are determined if the system is correctable
    if (cce_flags.flag("CORRECTABLE")){
      get_corrections(structure, cce_flags, cce_vars, cce_vars.anion_species, cce_vars.num_neighbors, cce_vars.corrections_atom);
      // determine corrections for cations bound to multi anion atoms if needed
      if (cce_flags.flag("MULTI_ANION_SYSTEM")){
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
          if(LDEBUG){
            cerr << endl;
            cerr << __AFLOW_FUNC__ << " getting corrections for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
          }
          get_corrections(structure, cce_flags, cce_vars, cce_vars.multi_anion_species[k], multi_anion_num_neighbors[k], cce_vars.multi_anion_corrections_atom[k]);
        }
      }
      // load per- & superox. corrections if needed;
      if (cce_vars.num_perox_bonds > 0 || cce_vars.num_superox_bonds > 0){
        check_get_per_super_ox_corrections(cce_vars);
      }
    }


    // CALCULATE CCE CORRECTIONS AT 298.15 AND 0K ##################################################
    // calculate CCE correction if the system is correctable
    if (cce_flags.flag("CORRECTABLE")){
      // calculate total correction per cell only using cations although 
      // corrections for anions should be set to zero
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){
        if ((structure.atoms[i].cleanname != cce_vars.anion_species) && (cce_vars.multi_anion_atoms[i] != 1)){ // exclude main anion species and multi anion atoms detected previously
          if (cce_vars.num_neighbors[i] > 0){ // are there actually bonds between the cation and the (main) anion species
            uint num_funcs=cce_vars.vfunctionals.size();
            for (uint k = 0; k < num_funcs; k++) {
              if (cce_vars.vfunctionals[k] != "exp") {
                for (uint l = 0; l < num_temps; l++) {
                  cce_vars.cce_correction[num_temps*k+l] += (cce_vars.num_neighbors[i]*cce_vars.corrections_atom[num_temps*k+l][i]) ; // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
                }
              } else { // for exp this calculates already the CCE@exp formation enthalpy (only at 298.15 K) from the exp. formation enthalpies per bond
                cce_vars.cce_correction[num_temps*k] += (cce_vars.num_neighbors[i]*cce_vars.corrections_atom[num_temps*k][i]) ;
              }
            }
          }
        }
      }
      // add correction for cations bound to multi anion atoms if needed
      if (cce_flags.flag("MULTI_ANION_SYSTEM")){
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
          if(LDEBUG){
            cerr << endl;
            cerr << __AFLOW_FUNC__ << " adding corrections for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
          }
          for(uint i=0,isize=structure.atoms.size();i<isize;i++){
            if ((structure.atoms[i].cleanname != cce_vars.anion_species) && (cce_vars.multi_anion_atoms[i] != 1)){ // exclude main anion species and multi anion atoms detected previously
              if (multi_anion_num_neighbors[k][i] > 0){ // are there actually bonds between the cation and the anion species under consideration
                uint num_funcs=cce_vars.vfunctionals.size();
                for (uint l = 0; l < num_funcs; l++) {
                  if (cce_vars.vfunctionals[l] != "exp") {
                    for (uint m = 0; m < num_temps; m++) {
                      cce_vars.cce_correction[num_temps*l+m] += (multi_anion_num_neighbors[k][i]*cce_vars.multi_anion_corrections_atom[k][num_temps*l+m][i]) ;
                    }
                  } else { // for exp this calculates already the CCE@exp formation enthalpy (only at 298.15 K) from the exp. formation enthalpies per bond
                    cce_vars.cce_correction[num_temps*l] += (multi_anion_num_neighbors[k][i]*cce_vars.multi_anion_corrections_atom[k][num_temps*l][i]) ;
                  }
                }
              }
            }
          }
        }
      }
      // add per- and superoxide corrections if needed
      if (cce_vars.num_perox_bonds > 0 || cce_vars.num_superox_bonds > 0){
        check_apply_per_super_ox_corrections(cce_vars);
      }
      // add ref. enthalpy shifts for PBE+U:ICSD if needed
      apply_pbe_u_icsd_shifts(structure, cce_vars);
    }
  } // main CCE function core



  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                   READ USER INPUT (FROM COMMAND LINE)                   //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //read_structure////////////////////////////////////////////////////////
  // read structural data from structure file provided on command line
  xstructure read_structure(const string& structure_file, int mode){ // first argument can be directly structure_file and not structure_file_path
    stringstream message;
    //string structure_file = aurostd::file2string(structure_file_path); // first argument of read_structure_function does not need to be converted to string since it contains already the file content and not only the file name
    xstructure structure(structure_file, mode);
    structure.checkStructure();
    return structure;
  }

  //read_structure////////////////////////////////////////////////////////
  // read structural data from istream
  xstructure read_structure(std::istream& ist){
    xstructure structure(ist);
    structure.checkStructure();
    return structure;
  }

  //get_dft_form_enthalpies_functionals////////////////////////////////////////////////////////
  // set the DFT formation enthalpies and functionals according to the input and check consistency of functionals and formation enthalpies
  void get_dft_form_enthalpies_functionals(const string& enthalpies_dft_input_str, const string& functionals_input_str, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    // vectorize precalculated DFT formation enthalpies if provided
    if(!enthalpies_dft_input_str.empty()){ //if the DFT enthalpies input string is not empty, convert it to vector of doubles
      aurostd::string2tokens<double>(enthalpies_dft_input_str,cce_vars.enthalpies_dft,",");
    } 
    // vectorize corresponding functionals
    if(!functionals_input_str.empty()){ // if functionals input string is not empty, convert it to vector of strings
      aurostd::string2tokens(functionals_input_str,cce_vars.vfunctionals,",");
    }
    //use PBE as default if only 1 DFT formation enthalpy is provided
    ostream& oss = cout;
    ofstream FileMESSAGE;
    _aflags aflags;aflags.Directory=aurostd::getPWD();
    if(functionals_input_str.empty() && cce_vars.enthalpies_dft.size() == 1){
      message << " Setting functionals=PBE since only 1 DFT formation enthalpy is provided and PBE is the default functional!";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      cce_vars.vfunctionals.push_back("PBE");
      // otherwise if sizes of provided DFT formation enthalpies and functionals do not match, throw error
      // if only functional argument is set corrections should only be returned for desired functional
    } else if(cce_vars.enthalpies_dft.size()!=cce_vars.vfunctionals.size() && !enthalpies_dft_input_str.empty() ){ 
      message << " BAD NEWS: The number of provided precalculated DFT formation enthalpies and functionals must match.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    //let the program spit out what it thinks
    if(LDEBUG){
      bool roff = false;
      cerr << __AFLOW_FUNC__ << " INPUT DFT FORMATION ENTHALPIES & FUNCTIONALS:" << endl;
      cerr << __AFLOW_FUNC__ << "  input dft formation enthalpies=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(cce_vars.enthalpies_dft,6,roff),",") << " (assumed to be in eV/cell)" << endl;
      cerr << __AFLOW_FUNC__ << "  input functionals=" << functionals_input_str << endl;
    }
    // determine whether it is a functional for which corrections are available
    vector<string> CCE_vallowed_functionals;
    aurostd::string2tokens(CCE_allowed_functionals, CCE_vallowed_functionals, ",");
    uint num_funcs=cce_vars.vfunctionals.size();
    for(uint k=0,ksize=num_funcs;k<ksize;k++){
      if(cce_vars.vfunctionals[k]!="exp"){
        cce_vars.vfunctionals[k]=aurostd::toupper(cce_vars.vfunctionals[k]);
      }
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " cce_vars.vfunctionals[" << k << "]: " << cce_vars.vfunctionals[k] << endl;
      }
      if (!aurostd::WithinList(CCE_vallowed_functionals, cce_vars.vfunctionals[k]) || get_offset(cce_vars.vfunctionals[k]) == -1) {
        message << " Unknown functional " << cce_vars.vfunctionals[k] << ". Please choose PBE, LDA, or SCAN.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
      cce_vars.offset.push_back(get_offset(cce_vars.vfunctionals[k]));
    }
    // if only structure (and oxidation numbers) are provided, corrections should be given for PBE, LDA, and SCAN 
    // and the CCE@exp formation enthalpy from the exp. formation enthalpies per bond
    // ICSD correction should only be returned when explicitly asked for
    if(functionals_input_str.empty() && enthalpies_dft_input_str.empty()){
      vector<string> CCE_vdefault_output_functionals;
      aurostd::string2tokens(CCE_default_output_functionals, CCE_vdefault_output_functionals, ",");
      for(uint k=0,ksize=CCE_vdefault_output_functionals.size();k<ksize;k++){
        if (get_offset(CCE_vdefault_output_functionals[k]) == -1) {
          message << " Unknown functional " << cce_vars.vfunctionals[k] << ". Please choose PBE, LDA, or SCAN.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        }
        cce_vars.vfunctionals.push_back(CCE_vdefault_output_functionals[k]); cce_vars.offset.push_back(get_offset(CCE_vdefault_output_functionals[k]));
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " PBE: " << aurostd::WithinList(cce_vars.vfunctionals, "PBE") << endl;
      cerr << __AFLOW_FUNC__ << " LDA: " << aurostd::WithinList(cce_vars.vfunctionals, "LDA") << endl;
      cerr << __AFLOW_FUNC__ << " SCAN: " << aurostd::WithinList(cce_vars.vfunctionals, "SCAN") << endl;
      cerr << __AFLOW_FUNC__ << " PBE+U:ICSD: " << aurostd::WithinList(cce_vars.vfunctionals, "PBE+U:ICSD") << endl;
      cerr << __AFLOW_FUNC__ << " exp: " << aurostd::WithinList(cce_vars.vfunctionals, "exp") << endl;
      cerr << endl;
    }
  }

  //get_offset////////////////////////////////////////////////////////
  // get offset needed for reading corrections from lookup table for different functionals
  int get_offset(const string& functional) {
    if (functional=="PBE")          {return 0;}
    else if (functional=="LDA")          {return 2;}
    else if (functional=="SCAN")         {return 4;}
    else if (functional=="PBE+U:ICSD")   {return 6;}
    else if (functional=="exp")          {return 8;}
    else {return -1;}
  }

  //get_oxidation_states////////////////////////////////////////////////////////
  // Retrieves the oxidation states of the material.
  vector<double> get_oxidation_states(const string& oxidation_numbers_input_str, const xstructure& structure, CCE_Variables& cce_vars, ostream& oss) {
    stringstream message;
    if(!oxidation_numbers_input_str.empty()){
      aurostd::string2tokens<double>(oxidation_numbers_input_str,cce_vars.oxidation_states,","); //if the oxidation numbers input string is not empty, convert it to vector of doubles
      //sizes of oxidation numbers and atoms must match
      if(cce_vars.oxidation_states.size()!=structure.atoms.size()){ 
        message << " BAD NEWS: The number of provided oxidation numbers does not match the number of atoms in the structure! Please correct and rerun.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
      // get sum of oxidation numbers and validate (system should not be regarded correctable if sum over oxidation states is not zero)
      cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars); // double because for superoxides O ox. number is -0.5
      if (std::abs(cce_vars.oxidation_sum) > DEFAULT_CCE_OX_TOL) {
        oss << print_output_oxidation_numbers(structure, cce_vars);
        message << " BAD NEWS: The formation enthalpy of this system is not correctable!" << endl;
        message << " The oxidation numbers that you provided do not add up to zero!" << endl;
        message << " Sum over all oxidation numbers is: " << cce_vars.oxidation_sum << endl;
        message << " Please correct and rerun." << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
    } else {
      message << " It seems you forgot to provide the oxidation numbers after \"--oxidation_numbers=\". Please add them or omit the option.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    return cce_vars.oxidation_states;
  }

  //get_functional_from_aflow_in_outcar////////////////////////////////////////////////////////
  // determine the functional from the aflow.in or the outcar file if the PP information cannot be determined from the aflow.in
  string get_functional_from_aflow_in_outcar(const xstructure& structure, string& aflowin_file, string& outcar_file) {
    stringstream message;
    string functional = "";
    string aflowIn = aurostd::RemoveComments(aurostd::file2string(aflowin_file));
    vector<string> vlines = aurostd::string2vectorstring(aflowIn);
    string line_a = "";
    bool pbe = false;
    bool ldau = false;
    bool pbe_u_icsd = false;
    bool ldau2 = false;
    bool lda = false;
    bool scan = false;
    _kflags kflags;
    _aflags aflags;
    _vflags vflags=KBIN::VASP_Get_Vflags_from_AflowIN(aflowIn,aflags,kflags);
    //_xvasp xvasp;xvasp.clear();
    //KBIN::readModulesFromAflowIn(aflowIn, kflags, xvasp);  //ME20181027
    //cout << "xvasp.str.species_pp: " << xvasp.str.species_pp << endl;
    //cout << "xvasp.str: " << xvasp.str << endl;
    // try to get PP information from aflow.in
    for (uint i = 0; i < vlines.size(); i++) {
      line_a = aurostd::RemoveSpaces(vlines[i]);
      // check whether it is a PBE calculation; PBE+U will be checked later
      if (line_a.find("=potpaw_PBE") != string::npos || line_a.find("potpaw_PBE/") != string::npos){
        pbe = true;
      }
      // check whether it is an LDA calculation
      if (line_a.find("=potpaw_LDA") != string::npos || line_a.find("potpaw_LDA/") != string::npos){
        lda = true;
      } // if both pbe and lda remain false and also no other functional can be determined, an error will be thrown at the end of this function
    }
    // if both functionals: pbe and lda are still false, the PP information was maybe not obtainable from the aflow.in and the OUTCAR should be checked, if it exists
    if (!pbe && !lda) {
      if (aurostd::FileExist(outcar_file) || aurostd::EFileExist(outcar_file)) {
        xOUTCAR outcar;
        outcar.GetPropertiesFile(outcar_file);
        uint pbe_count=0;
        uint lda_count=0;
        for(uint i=0;i<outcar.species_pp_type.size();i++) {
          if (outcar.species_pp_type.at(i).find("PAW_PBE") != string::npos){
            pbe_count+=1;
          } else if (outcar.species_pp_type.at(i).find("PAW_LDA") != string::npos){
            lda_count+=1;
          }
        }
        if (outcar.species_pp_type.size() == pbe_count){
          pbe = true;
        } else if (outcar.species_pp_type.size() == lda_count){
          lda = true;
        } // if both pbe and lda remain false and also no other functional can be determined, an error will be thrown at the end of this function
      }
    }
    //check for SCAN
    if (vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme == "SCAN"){
      scan = true;
    }
    // check whether it is a DFT+U calculation with parameters as for the ICSD (PBE+U:ICSD calculation)
    // new implementation checking Us explicitly
    if (!vflags.KBIN_VASP_LDAU_PARAMETERS.empty() && vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry && !vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry && !vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry){
      ldau2=true;
      vector<string> ldau_params_vector;
      aurostd::string2tokens(vflags.KBIN_VASP_LDAU_PARAMETERS, ldau_params_vector, ";");
      if (ldau_params_vector.size() != 4){
        message << " BAD NEWS: The LDAU parameters are not of the right size. Please adapt and rerun.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
      // get species
      string species_part = ldau_params_vector[0];
      vector<string> species_vector;
      aurostd::string2tokens(species_part, species_vector, ",");
      if (species_vector.size() != structure.species.size()){
        message << " BAD NEWS: The number of species in the DFT+U settings differs from the total number of species for this structure. Please adapt and rerun.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
      // get Us
      string Us_part = ldau_params_vector[2];
      vector<string> Us_string_vector;
      vector<double> Us_vector;
      aurostd::string2tokens(Us_part, Us_vector,",");
      if (species_vector.size() != Us_vector.size()){
        message << " BAD NEWS: The number of species in the DFT+U settings differs from the number of provided U values. Please adapt and rerun.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      }
      // get standard Us for PBE+U:ICSD
      vector<double> standard_ICSD_Us_vector;
      bool LDAU=true;
      vector<string> vLDAUspecies;
      vector<uint> vLDAUtype;
      vector<int> vLDAUL;
      vector<double> vLDAUU;
      vector<double> vLDAUJ;
      for (uint k = 0; k < species_vector.size(); k++) {
        AVASP_Get_LDAU_Parameters(species_vector[k],LDAU,vLDAUspecies,vLDAUtype,vLDAUL,vLDAUU,vLDAUJ);
        if (vLDAUtype[k] != 2 && vLDAUU[k] !=0){
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          message << " BAD NEWS: It seems the standard DFT+U method for " << species_vector[k] << " was changed from Dudarev's approach (LDAU2=ON, vLDAUtype[" << k << "]=2)" << std::endl;
          message << "as used for the AFLOW ICSD database when obtaining the corrections to vLDAUtype[" << k << "]=" << vLDAUtype[k] << " now." << std::endl;
          message << "If the standard U values have also been changed, then the corrections might have been constructed for other U values than used in this calculation" << std::endl;
          message << "and should not be applied. Please check this carefully!";
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
        standard_ICSD_Us_vector.push_back(vLDAUU[k]);
      }
      // compare read and standard Us
      for (uint k = 0; k < species_vector.size(); k++) {
        if (Us_vector[k] != standard_ICSD_Us_vector[k]){
          message << " BAD NEWS: For this DFT+U calculation with Dudarev's method the provided U value of " << Us_vector[k] << " eV for " << species_vector[k] << " does not match the standard value of " << standard_ICSD_Us_vector[k] << " eV. There are no corrections for this case.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        }
      }
      pbe_u_icsd = true;
    }
    // check whether it is a DFT+U calculation with different parameters than for PBE+U:ICSD
    if ((vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) && !pbe_u_icsd){
      message << " BAD NEWS: It seems you are providing an aflow.in for a DFT+U calculation with different parameters than for the AFLOW ICSD database (Dudarev's approach, LDAU2=ON). There are no corrections for this case.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }

    if (pbe){
      if (pbe_u_icsd){
        functional = "PBE+U:ICSD";
      } else if (scan){ //SCAN calculations are usually done with PBE PPs
        functional = "SCAN";
      } else if (!ldau && !ldau2){
        functional = "PBE";
      }
    } else if (lda){
      if (scan){
        functional = "SCAN"; //SCAN calcs. with LDA PPs should be okay
      } else if (!ldau && !ldau2){
        functional = "LDA";
      }
    } else if (scan){
      if (!ldau && !ldau2){
        functional = "SCAN";
      }
    }
    // if functional is still empty, i.e. cannot be determined from aflow.in, throw error
    if (functional.empty()) {
      message << " Functional cannot be determined from aflow.in. Corrections are available for PBE, LDA, SCAN, or PBE+U:ICSD.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    return functional;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                      INITIALISE FLAGS AND VARIABLES                     //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //init_flags////////////////////////////////////////////////////////////
  //ME20200213
  // Initializes the CCE flags to their default values.
  aurostd::xoption init_flags() {
    aurostd::xoption flags;
    flags.flag("UNIT_TEST", false);
    flags.flag("RUN_FULL_CCE", false);
    flags.flag("CORRECTABLE", true); // first assuming that formation enthalpy of system IS correctable; will be set to not correctable if, for any atom, no correction can be identified
    flags.flag("OX_NUMS_PROVIDED", false);
    flags.flag("MULTI_ANION_SYSTEM", false);
    flags.flag("O_MULTI_ANION_SPECIES", false); // whether one of the multi anion species is O for which per/superoxide tests need to be made
    return flags;
  }

  //init_variables////////////////////////////////////////////////////////
  //ME20200213
  // Initializes the variables struct.
  CCE_Variables init_variables(const xstructure& structure) {
    CCE_Variables cce_vars;
    cce_vars.enthalpies_dft.clear();
    cce_vars.vfunctionals.clear();
    cce_vars.offset.clear();
    cce_vars.num_perox_bonds=0;
    cce_vars.num_superox_bonds=0;
    // clear and resize all vectors
    cce_vars.electronegativities.clear();
    cce_vars.electronegativities.resize(structure.species.size());
    cce_vars.multi_anion_atoms.clear();
    cce_vars.multi_anion_atoms.resize(structure.atoms.size());
    cce_vars.oxidation_states.clear();
    cce_vars.oxidation_states.resize(structure.atoms.size());
    cce_vars.multi_anion_species.clear();
    cce_vars.perox_indices.clear();
    cce_vars.perox_indices.resize(structure.atoms.size());
    cce_vars.superox_indices.clear();
    cce_vars.superox_indices.resize(structure.atoms.size());
    cce_vars.num_neighbors.clear();
    cce_vars.num_neighbors.resize(structure.atoms.size());
    cce_vars.species_electronegativity_sorted.clear();
    cce_vars.species_electronegativity_sorted.resize(structure.species.size());
    cce_vars.num_pref_ox_states_electronegativity_sorted.clear();
    cce_vars.num_pref_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.num_all_ox_states_electronegativity_sorted.clear();
    cce_vars.num_all_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.pref_ox_states_electronegativity_sorted.clear();
    cce_vars.pref_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.all_ox_states_electronegativity_sorted.clear();
    cce_vars.all_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.cations_map.clear();
    cce_vars.Bader_charges.clear();
    cce_vars.corrections_atom.clear();
    cce_vars.multi_anion_corrections_atom.clear();
    cce_vars.perox_correction.clear();
    cce_vars.superox_correction.clear();
    cce_vars.cce_correction.clear();
    cce_vars.enthalpy_formation_cell_cce.clear();
    return cce_vars;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           STRUCTURAL ANALYSIS                           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // counting the number of anion neighbors (multiple anion species possible) for each cation

  //determine_anion_species////////////////////////////////////////////////////////
  // determine anion species
  string determine_anion_species(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ANION SPECIES FROM ALLEN ELECTRONEGATIVITIES:" << endl;
    }
    _atom atom;
    uint z = 0;

    double anion_electronegativity = 0;
    for(uint k=0,ksize=structure.species.size();k<ksize;k++){
      z = GetAtomNumber(KBIN::VASP_PseudoPotential_CleanName(structure.species[k]));
      xelement::xelement element(z);
      if (element.electronegativity_Allen == NNN) {
        message << " VERY BAD NEWS: There is no known electronegativity value for " << KBIN::VASP_PseudoPotential_CleanName(structure.species[k]) << ".";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
      } else{
        cce_vars.electronegativities[k] = element.electronegativity_Allen;
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " electronegativity of species " << k << " (" << KBIN::VASP_PseudoPotential_CleanName(structure.species[k]) << "): " << cce_vars.electronegativities[k] << endl;
        }
        if (cce_vars.electronegativities[k] > anion_electronegativity) {
          anion_electronegativity = cce_vars.electronegativities[k];
          cce_vars.anion_species = KBIN::VASP_PseudoPotential_CleanName(structure.species[k]);
        }
      }
    }
    // set anion charge and check whether it is negative
    z = GetAtomNumber(cce_vars.anion_species);
    xelement::xelement element(z);
    cce_vars.standard_anion_charge = element.oxidation_states[element.oxidation_states.size()-1];
    if (cce_vars.standard_anion_charge > 0) {
      message << " VERY BAD NEWS: There is no known negative oxidation number for " << cce_vars.anion_species << " detected as anion species.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " anion electronegativity: " << anion_electronegativity << endl;
      cerr << __AFLOW_FUNC__ << " anion species: " << cce_vars.anion_species << endl;
      cerr << __AFLOW_FUNC__ << " anion charge: " << cce_vars.standard_anion_charge << endl;
      cerr << endl;
    }
    return cce_vars.anion_species;
  }

  //check_for_multi_anion_system////////////////////////////////////////////////////////
  // check whether it is a multi-anion system, i. e. whether atoms of another species than the anion_species are only bound to atoms of lower electronegativity or of the same type
  vector<uint> check_for_multi_anion_system(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, double tolerance) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " CHECKING FOR MULTI-ANION SYSTEM:" << endl;
    }
    for ( uint i = 0; i < structure.atoms.size(); i++ ) { // initialize elements of vector to 0 
      cce_vars.multi_anion_atoms[i] = 0; 
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " structure=" << endl;cerr << structure << endl;}

    //CO20200914 - START
    if(cce_vars.xstr_neighbors.atoms.size()==0 || cce_vars.i_neighbors.size()!=structure.atoms.size() || cce_vars.i_neighbors.size()!=cce_vars.distances.size()){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEIGHBORS ANALYSIS START" << endl;}
      cce_vars.xstr_neighbors=structure;
      cce_vars.xstr_neighbors.GetNeighbors(cce_vars.i_neighbors,cce_vars.distances,0.0,false,false);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEIGHBORS ANALYSIS STOP" << endl;}
    }
    //CO20200914 - END

    uint i=0,isize=0,j=0,jsize=0,neighbors_count=0,multi_anion_count=0;
    double cutoff=0.0,electronegativity_atom=0.0,electronegativity_neighbor=0.0;
    xelement::xelement atom_element,neigh_element;
    for(i=0,isize=cce_vars.i_neighbors.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      neighbors_count=0;
      multi_anion_count=0;
      cutoff=cce_vars.distances[i][0]+tolerance;  //add tolerance
      for(j=0,jsize=cce_vars.i_neighbors[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by cutoffs_max
        if(cce_vars.distances[i][j]>cutoff){break;}  //we only need first nearest neighbors
        const _atom& atom=cce_vars.xstr_neighbors.grid_atoms[cce_vars.i_neighbors[i][j]];
        if (cce_vars.distances[i][j] <= cutoff) { // distance must be larger than DEFAULT_CCE_SELF_DIST_TOL since GetStrNeighData includes also structure.atoms[i] itself as neighbor having distance zero to itself
          if (atom.cleanname == cce_vars.anion_species){
            neighbors_count+=1;
          } else if ((atom.cleanname != cce_vars.anion_species) && (structure.atoms[i].cleanname != cce_vars.anion_species)){ // second condition set since the anion_species cannot be set as a multi-anion species again
            neighbors_count+=1;
            atom_element.populate(structure.atoms[i].cleanname);
            neigh_element.populate(atom.cleanname);
            electronegativity_atom = atom_element.electronegativity_Allen;
            electronegativity_neighbor = neigh_element.electronegativity_Allen;
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " electronegativity of atom " << i << ": " << electronegativity_atom << endl;
              cerr << __AFLOW_FUNC__ << " electronegativity of neighbor " << j << ": " << electronegativity_neighbor << endl;
            }
            if((electronegativity_neighbor < electronegativity_atom) || (structure.atoms[i].cleanname == atom.cleanname)){ // could be multi-anion atom if bound to only atoms of lower electroneg. or of same species
              multi_anion_count+=1;
            }
          }
        }
      }
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " multi_anion_count for atom " << i << " (" << structure.atoms[i].cleanname << "): " << multi_anion_count << endl;
        cerr << __AFLOW_FUNC__ << " neighbors_count for atom " << i << " (" << structure.atoms[i].cleanname << "): " << neighbors_count << endl;
      }
      if ((multi_anion_count == neighbors_count) && (structure.atoms[i].cleanname != cce_vars.anion_species)){ // anion_species should not be detected again as multi_anion_species
        if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
          cce_flags.flag("MULTI_ANION_SYSTEM",TRUE);
        }
        // set multi anion atoms
        cce_vars.multi_anion_atoms[i]=1;
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " Atom " << i << " (" << structure.atoms[i].cleanname << ") has been detected as a multi-anion atom." << endl;
        }
        // set multi anion species
        uint multi_anion_species_count=0;
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ // have to make sure that no anion species is added twice
          if(cce_vars.multi_anion_species[k] != structure.atoms[i].cleanname){
            multi_anion_species_count+=1;
          }
        }
        if(multi_anion_species_count == cce_vars.multi_anion_species.size()){
          cce_vars.multi_anion_species.push_back(structure.atoms[i].cleanname);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " New multi anion species: " << structure.atoms[i].cleanname << endl;
          }
          if(structure.atoms[i].cleanname == "O"){ // check whether one of the multi anion species is O for which per/superoxide test needs to be made
            cce_flags.flag("O_MULTI_ANION_SPECIES",TRUE);
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " Oxygen was detected as multi anion species, i.e. system needs to be tested for (su-)peroxide ions." << endl;
            }
          }
        }
        // set multi anion oxidation numbers and check whether it is negative
        atom_element.populate(structure.atoms[i].cleanname);
        cce_vars.oxidation_states[i] = atom_element.oxidation_states[atom_element.oxidation_states.size()-1];
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " Oxidation state for atom " << i << " (" << structure.atoms[i].cleanname << ") has been set to: " << cce_vars.oxidation_states[i] << endl;
        }
        if (cce_vars.oxidation_states[i] > 0) {
          message << " VERY BAD NEWS: There is no known negative oxidation number for " << structure.atoms[i].cleanname << " detected as multi anion species.";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        }
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of multi anion species (1-total_number_anion_species since the main anion species with largest electronegativity is NOT counted as multi anion species): " << cce_vars.multi_anion_species.size() << endl;
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        cerr << __AFLOW_FUNC__ << " multi anion species " << k << ": " << cce_vars.multi_anion_species[k] << endl;
      }
      cerr << endl;
    }
    return cce_vars.multi_anion_atoms;
  }

  //get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors for each atom within respective species selective cutoff
  vector<uint> get_num_neighbors(const xstructure& structure, double tolerance) {
    CCE_Variables cce_vars = init_variables(structure);
    string anion_species="";
    aurostd::xoption cce_flags = init_flags();
    return get_num_neighbors(structure, anion_species, cce_flags, cce_vars, tolerance);
  }

  //get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors of a special anion_species for each atom within respective species selective cutoff
  vector<uint> get_num_neighbors(const xstructure& structure, const string& anion_species, double tolerance) {
    CCE_Variables cce_vars = init_variables(structure);
    aurostd::xoption cce_flags = init_flags();
    return get_num_neighbors(structure, anion_species, cce_flags, cce_vars, tolerance);
  }

  //get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors for each atom within respective species selective cutoff
  vector<uint> get_num_neighbors(const xstructure& structure, const string& anion_species, xoption& cce_flags, CCE_Variables& cce_vars, double tolerance) { // anion_species here cannot be taken from cce_vars since function is also used to determine multi anion num_neighbors for which anion_species is the respective multi_anion_species
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " STRUCTURAL ANALYSIS:" << endl;
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " structure=" << endl;cerr << structure << endl;}

    vector<uint> num_neighbors(structure.atoms.size());

    //CO20200914 - START
    if(cce_vars.xstr_neighbors.atoms.size()==0 || cce_vars.i_neighbors.size()!=structure.atoms.size() || cce_vars.i_neighbors.size()!=cce_vars.distances.size()){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEIGHBORS ANALYSIS START" << endl;}
      cce_vars.xstr_neighbors=structure;
      cce_vars.xstr_neighbors.GetNeighbors(cce_vars.i_neighbors,cce_vars.distances,0.0,false,false);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEIGHBORS ANALYSIS STOP" << endl;}
    }
    //CO20200914 - END

    uint i=0,isize=0,j=0,jsize=0,neighbors_count=0;
    bool warning=false;
    stringstream message,other_neighbors;
    double cutoff=0.0;
    for(i=0,isize=cce_vars.i_neighbors.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      neighbors_count=0;
      warning = false;
      aurostd::StringstreamClean(other_neighbors);
      cutoff=cce_vars.distances[i][0]+tolerance;  //add tolerance
      for(j=0,jsize=cce_vars.i_neighbors[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by cutoffs_max
        if(cce_vars.distances[i][j]>cutoff){break;}  //we only need first nearest neighbors
        const _atom& atom=cce_vars.xstr_neighbors.grid_atoms[cce_vars.i_neighbors[i][j]];
        if (cce_vars.distances[i][j] <= cutoff) { // distance must be larger than DEFAULT_CCE_SELF_DIST_TOL since GetStrNeighData includes also structure.atoms[i] itself as neighbor having distance zero to itself
          if (!anion_species.empty()){ // variable called anion type since function was developed for CCE for polar materials but it can be used to check for any atom type and only include those as neighbors
            // implement check whether each nearest neighbor is of the anion_species, otherwise throw warning; 
            if (atom.cleanname == anion_species){
              neighbors_count+=1;
            } else if ((atom.cleanname != anion_species) && (structure.atoms[i].cleanname != anion_species)){ // second condition set since it is expected that the anion has predominantly other neighbors than its own type 
              if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
                if (!warning) { warning = true; }
                other_neighbors << "Not all nearest neighbors of " << structure.atoms[i].cleanname << " (ATOM[" << i << "]) within the distance tolerance of " << tolerance << " Ang. are " << anion_species << ", there is also " << atom.cleanname << endl;
              }
            }
          } else {
            neighbors_count+=1;
          }
        }
      }
      if (!cce_flags.flag("MULTI_ANION_SYSTEM") && warning){
        message << " Not all nearest neighbors of " << structure.atoms[i].cleanname << " (ATOM[" << i << "]) within the distance tolerance are " << anion_species << "!" << endl;
        message << other_neighbors.str();
        if (!cce_flags.flag("UNIT_TEST")){
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
      }
      num_neighbors[i]=neighbors_count; // zero-based counting as for cutoffs array above
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " number of " << anion_species << " nearest neighbors within " << tolerance << " Ang. tolerance of " << structure.atoms[i].cleanname << " (ATOM[" << i << "]): " << num_neighbors[i] << endl;
        cerr << endl;
      }
    }
    return num_neighbors;
  }

  //get_dist_cutoffs////////////////////////////////////////////////////////
  // determine species selective nearest neighbor distances and then cutoffs accordingly
  vector<double> get_dist_cutoffs(const xstructure& structure) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    vector<double> cutoffs(structure.species.size());
    xmatrix<double> distsij=GetDistMatrix(structure); // gives matrix with nearest neighbor distances between all species pairs with species running over rows and columns
    vector<double> near_neigh_dist(structure.species.size());
    for(int i=1,isize=distsij.rows+1;i<isize;i++){
      for(int j=1,jsize=distsij.cols+1;j<jsize;j++){
        if (j==1){
          near_neigh_dist[i-1]=distsij(i,j);
        }
        if (distsij(i,j) < near_neigh_dist[i-1]){
          near_neigh_dist[i-1]=distsij(i,j);
        }
      }
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " nearest neighbor distance for species " << i << " is " << near_neigh_dist[i-1] << " Ang." << endl;
      }
      cutoffs[i-1]=near_neigh_dist[i-1];  //CO20200914 - NO TOLERANCE ADDED HERE, this is local to the function
    }
    return cutoffs;
  }

  //check_per_super_oxides////////////////////////////////////////////////////////
  // check whether the system contains per- or superoxide ions based on the O-O bond length and set variables accordingly
  void check_per_super_oxides(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " CHECKING FOR (SU-)PEROXIDES:" << endl;
    }
    uint perox_count=0;
    uint superox_count=0;
    for ( uint i = 0; i < structure.atoms.size(); i++ ) { // initialize elements of vectors to 0 
      cce_vars.perox_indices[i] = 0; 
      cce_vars.superox_indices[i] = 0; 
    }

    //CO20200914 - START
    if(cce_vars.xstr_neighbors.atoms.size()==0 || cce_vars.i_neighbors.size()!=structure.atoms.size() || cce_vars.i_neighbors.size()!=cce_vars.distances.size()){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEIGHBORS ANALYSIS START" << endl;}
      cce_vars.xstr_neighbors=structure;
      cce_vars.xstr_neighbors.GetNeighbors(cce_vars.i_neighbors,cce_vars.distances,0.0,false,false);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " NEIGHBORS ANALYSIS STOP" << endl;}
    }
    //CO20200914 - END

    uint i=0,isize=0,j=0,jsize=0;
    for(i=0,isize=cce_vars.i_neighbors.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      if (structure.atoms[i].cleanname == "O"){ // identify per- and superoxides by O-O bond length
        for(j=0,jsize=cce_vars.i_neighbors[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by the cutoffs_max
          if(cce_vars.distances[i][j]>cce_vars.distances[i][0]){break;}  //we only need first nearest neighbors
          const _atom& atom=cce_vars.xstr_neighbors.grid_atoms[cce_vars.i_neighbors[i][j]];
          if (atom.cleanname == "O"){
            if (cce_vars.distances[i][j] <= DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF) { // distance must be larger than DEFAULT_CCE_SELF_DIST_TOL to savely exclude the anion itself having distance zero to itself; if O-O bond is shorter than in O2 molecule (approx. 1.21 Ang) the result of the structural relaxation is most likely wrong
              message << " THE DETERMINED OXYGEN-OXYGEN BOND LENGTH IS SHORTER THAN IN THE O2 MOLECULE; CHECK YOUR STRUCTURE! THE O-O BOND LENGTH IS: " << cce_vars.distances[i][j] << " Ang.";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
            } else if ((DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF < cce_vars.distances[i][j]) && (cce_vars.distances[i][j] < DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF) ){
              message << " THE DETERMINED OXYGEN-OXYGEN BOND LENGTH IS ABOUT THE SAME AS IN THE O2 MOLECULE, I.E. THE STRUCTURE SEEMS TO INCLUDE MOLECULAR OXYGEN FOR WHICH NO CCE CORRECTION IS AVAILABLE! THE O-O BOND LENGTH IS: " << cce_vars.distances[i][j] << " Ang.";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
            } else if ((DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF <= cce_vars.distances[i][j]) && (cce_vars.distances[i][j] <= DEFAULT_CCE_SUPEROX_CUTOFF) ){
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " WARNING: This should be a superoxide; the O-O bond length is: " << cce_vars.distances[i][j] << " Ang." << endl;
              }
              superox_count+=1;
              cce_vars.superox_indices[i]=1;
            } else if ((DEFAULT_CCE_SUPEROX_CUTOFF < cce_vars.distances[i][j]) && (cce_vars.distances[i][j] <= DEFAULT_CCE_PEROX_CUTOFF) ){
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " WARNING: This should be a peroxide; the O-O bond length is: " << cce_vars.distances[i][j] << " Ang." << endl;
              }
              perox_count+=1;
              cce_vars.perox_indices[i]=1;
            }
          }
        }
      }
    }
    cce_vars.num_perox_bonds=perox_count/2; // needs to be divided by two due to double counting of O-O bonds when looping over all atoms of the structure
    cce_vars.num_superox_bonds=superox_count/2;
    // print out the number of per- & superoxide O-O bonds;
    if (cce_vars.num_perox_bonds > 0){
      stringstream message;
      message << " This should be a peroxide!" << endl;
      message << "Number of peroxide O-O bonds in cell: " << cce_vars.num_perox_bonds << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
    } else if (cce_vars.num_superox_bonds > 0) {
      stringstream message_so;
      message_so << " This should be a superoxide!" << endl;
      message_so << "Number of superoxide O-O bonds in cell: " << cce_vars.num_superox_bonds << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message_so, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //           DETERMINE OXIDATION STATES FROM ELECTRONEGATIVITIES           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // determine oxidation numbers based on Allen electronegativities and treat mixed valence binary systems specially

  //get_oxidation_states_from_electronegativities////////////////////////////////////////////////////////
  // function overloading for below function to be able to use oxidation number determination independently of CCE
  vector<double> get_oxidation_states_from_electronegativities(xstructure& structure) {
    // check structure
    structure.checkStructure(); //includes rescaling the structure to 1
    CCE_Variables cce_vars = init_variables(structure);
    cce_vars.anion_species=determine_anion_species(structure, cce_vars);
    aurostd::xoption cce_flags = init_flags();
    cce_vars.multi_anion_atoms=check_for_multi_anion_system(structure, cce_flags, cce_vars);
    // multi anion corrections can only be resized after number of multi anion species is known from check for multi anion system
    cce_vars.multi_anion_corrections_atom.resize(cce_vars.multi_anion_species.size(), vector<vector<double> >(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size())));
    //cce_vars.num_neighbors=get_num_neighbors(structure, cce_vars.anion_species, cce_flags, cce_vars);
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      check_per_super_oxides(structure, cce_flags, cce_vars);
    }
    cce_flags.flag("CORRECTABLE",TRUE);
    return get_oxidation_states_from_electronegativities(structure, cce_flags, cce_vars);
  }

  //get_oxidation_states_from_electronegativities////////////////////////////////////////////////////////
  // determine the oxidation numbers of the ions using preferred/all known oxidation numbers, electronegativities and structural information
  vector<double> get_oxidation_states_from_electronegativities(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " DETERMINATION OF OXIDATION NUMBERS FROM PREFERRED/ALL KNOWN OXIDATION STATES, ALLEN ELECTRONEGATIVITIES, AND STRUCTURE:" << endl;
    }
    // deal with anion charge in cell first
    set_anion_oxidation_states(structure, cce_vars);
    // treat cations
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " CATION PART:" << endl;
    }
    // determine number of cation species (currently just all species minus one anion species, multi-anion atoms are dealt with separately)
    uint num_cation_species = structure.species.size()-1;
    // sort species ascending by electronegativity
    // preferred and all known oxidation states will be automatically sorted by electronegativity during subsequent loading
    // the anion species should always be the last one with the highest electronegativity
    sort_species_by_electronegativity(structure, cce_vars);
    // load needed info about preferred and all known oxidation states for all species (sorted by electronegativity)
    cce_flags.flag("NO_PREF_OX_STATES",FALSE);
    cce_flags.flag("NO_OX_STATES",FALSE);
    cce_flags.flag("OX_STATES_DETERMINED",FALSE);
    load_ox_states_templates_each_species(structure, cce_flags, cce_vars);
    //ME20191101 for getting cations_map: a vector of vectors that lists for each cation species the atom numbers of the structure that are of this species (for Fe2ZnO4 there might be two Fe atoms at positions 0 and 1 in the structure)
    uint natoms = structure.atoms.size();
    cce_vars.cations_map.resize(num_cation_species);
    uint i = 0;
    for (uint at = 0; at < natoms; at++) {
      for (i = 0; i < num_cation_species; i++) {
        if (structure.atoms[at].cleanname == cce_vars.species_electronegativity_sorted[i]) break; // if it finds that the atom belongs to the ith species sorted by electronegativities, break to insert it at the proper place for the ith species into cation_map
      }
      if (i < num_cation_species) cce_vars.cations_map[i].push_back(at);
    }
    // try to find proper oxidation states (making oxidation sum zero) by using preferred oxidation states
    if(!cce_flags.flag("NO_PREF_OX_STATES")){ // for He, Ne, Ar, and Xe there are no preferred ox. states
      try_preferred_oxidation_states(structure, cce_flags, cce_vars);
    }
    // if preferred oxidation states do not work, it could be a mixed valence (binary) system
    // treat them here as special cases
    if(!cce_flags.flag("OX_STATES_DETERMINED")){
      // for SbO2 the oxidation states are not identified properly, Bader analysis finds them automatically
      treat_SbO2_special_case(structure, cce_flags, cce_vars);
      // for Pb3O4 the oxidation states are not identified properly, Bader analysis finds them automatically
      treat_Pb3O4_special_case(structure, cce_flags, cce_vars);
      // Ti-O Magneli phases need to be treated specially since oxidation numbers are not recognized appropriately by above approach
      treat_Ti_O_Magneli_phase_special_case(structure, cce_flags, cce_vars);
      // for Fe3O4 in inverse spinel structure, the oxidation states are not identified properly
      treat_Fe3O4_special_case(structure, cce_flags, cce_vars);
      // for Mn3O4 in spinel structure, the oxidation states are not identified properly
      treat_X3O4_special_case(structure, cce_flags, cce_vars, "Mn");
      // for Co3O4 in spinel structure, the oxidation states are not identified properly
      treat_X3O4_special_case(structure, cce_flags, cce_vars, "Co");
      // alkali metal sesquioxides need to be treated specially since oxidation numbers and number of per- and superoxide bonds are not recognized appropriately
      treat_alkali_sesquioxide_special_case(structure, cce_flags, cce_vars);
    }
    // if preferred oxidation states approach and mixed valence doesn't work, try to find proper oxidation states (making oxidation sum zero) by using all known oxidation states
    if(!cce_flags.flag("OX_STATES_DETERMINED")){
      if(!cce_flags.flag("NO_OX_STATES")){ // for He, Ne, and Ar there are no known oxidation states
        try_all_oxidation_states(structure, cce_vars);
        // calculate sum of oxidation numbers
        cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > DEFAULT_CCE_OX_TOL) {
          cce_flags.flag("CORRECTABLE",FALSE);
          oss << print_output_oxidation_numbers(structure, cce_vars); // print previously gathered output e.g. from determination of oxidation numbers
          message << "BAD NEWS: The determined oxidation numbers do not add up to zero!"  << endl;
          message << "Sum over all oxidation numbers is: " << cce_vars.oxidation_sum << endl;
          if(cce_flags.flag("RUN_FULL_CCE")){
            message << "The formation enthalpy of this system is hence not correctable!"  << endl;
            message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
          }
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
        }
      } else{ // error message that oxidation numbers cannot be determined since for at least one species there are no known oxidation numbers is already included in load_ox_states_templates_each_species function
        cce_flags.flag("CORRECTABLE",FALSE);
        if(cce_flags.flag("RUN_FULL_CCE")){
          message << "The formation enthalpy of this system is hence not correctable!"  << endl;
          message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
        }
      }
    }
    return cce_vars.oxidation_states;
  }

  //set_anion_oxidation_states////////////////////////////////////////////////////////
  // determine the oxidation numbers of the anions
  void set_anion_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << __AFLOW_FUNC__ << " ANION PART:" << endl;
    }
    double total_anion_charge=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if ((structure.atoms[i].cleanname == cce_vars.anion_species) || (cce_vars.multi_anion_atoms[i] == 1)){
        if (structure.atoms[i].cleanname == cce_vars.anion_species){ // for multi-anion atoms oxidation states have been assigned previously
          cce_vars.oxidation_states[i] = cce_vars.standard_anion_charge; // oxidation numbers for O are assumed to be -2 (other anions accordingly) and are corrected below if it is a per- or superoxide O atom as identified from the structural analysis in other functions
        }
        if (cce_vars.num_perox_bonds > 0){
          if (cce_vars.perox_indices[i]==1){
            cce_vars.oxidation_states[i]=-1;
          }
        }
        if (cce_vars.num_superox_bonds > 0){
          if (cce_vars.superox_indices[i]==1){
            cce_vars.oxidation_states[i]=-0.5;
          }
        }
        if(LDEBUG){
          if (cce_vars.multi_anion_atoms[i] == 1){ // for multi anion atoms oxidation states have been assigned previously
            cerr << __AFLOW_FUNC__ << " anion oxidation number for multi-anion ATOM[" << i << "] (" << structure.atoms[i].cleanname << ") has been assigned previously to: " << cce_vars.oxidation_states[i] << endl;
          } else {
            cerr << __AFLOW_FUNC__ << " anion oxidation number for ATOM[" << i << "] (" << structure.atoms[i].cleanname << "): " << cce_vars.oxidation_states[i] << endl;
          }
        }
        total_anion_charge += cce_vars.oxidation_states[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " Total anion charge in cell: " << total_anion_charge << endl;
    }
  }

  //sort_species_by_electronegativity////////////////////////////////////////////////////////
  // sort species ascending by electronegativity
  void sort_species_by_electronegativity(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    // using aurostd sort functions
    vector<double> electronegativities_sorted = cce_vars.electronegativities;
    for(uint j=0,jsize=structure.species.size();j<jsize;j++){ //loop over all species
      cce_vars.species_electronegativity_sorted[j] = KBIN::VASP_PseudoPotential_CleanName(structure.species[j]);
    }
    // for the oxidation state algorithm the species must be sorted by electronegativity and the preferred 
    // and all known oxidation states will be automatically sorted by electronegativity in the subsequent loading
    // sort species by electronegativity
    aurostd::sort(electronegativities_sorted,cce_vars.species_electronegativity_sorted);
    for(uint j=0,jsize=structure.species.size();j<jsize;j++){ //loop over all species
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " species_electronegativity_sorted[" << j << "]: " << cce_vars.species_electronegativity_sorted[j] << endl;
        cerr << __AFLOW_FUNC__ << " electronegativities_sorted[" << j << "]: " << electronegativities_sorted[j] << endl;
      }
    }
  }

  //load_ox_states_templates_each_species////////////////////////////////////////////////////////
  // load templates for preferred and other oxidation states
  void load_ox_states_templates_each_species(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    _atom atom;
    uint z = 0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){ 
      z = GetAtomNumber(KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]));
      xelement::xelement element(z);
      // load preferred oxidation states for each species
      cce_vars.num_pref_ox_states_electronegativity_sorted[i] = element.oxidation_states_preferred.size();
      if(element.oxidation_states_preferred[0] != NNN){
        for (uint k=0,ksize=cce_vars.num_pref_ox_states_electronegativity_sorted[i];k<ksize;k++) {
          cce_vars.pref_ox_states_electronegativity_sorted[i].push_back(element.oxidation_states_preferred[k]);
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " num_pref_ox_states_electronegativity_sorted[" << i << "]: " << cce_vars.num_pref_ox_states_electronegativity_sorted[i] << endl;
          for (uint k=0,ksize=cce_vars.num_pref_ox_states_electronegativity_sorted[i];k<ksize;k++) {
            cerr << __AFLOW_FUNC__ << " preferred oxidation state " << k << " of species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << "): " <<  cce_vars.pref_ox_states_electronegativity_sorted[i][k] << endl;
          }
        }
      } else{
        cce_vars.num_pref_ox_states_electronegativity_sorted[i] = 0;
        cce_flags.flag("NO_PREF_OX_STATES",TRUE);
        if(LDEBUG){
          cerr << endl;
          cerr << __AFLOW_FUNC__ << " BAD NEWS: There are no preferred oxidation states for species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << ")."  << endl;
          cerr << __AFLOW_FUNC__ << " Therefore the oxidation states cannot be determined on this basis." << endl;
        }
      }
      // load all oxidation states for each species
      cce_vars.num_all_ox_states_electronegativity_sorted[i] = element.oxidation_states.size();
      if(element.oxidation_states[0] != NNN){
        for (uint k=0,ksize=cce_vars.num_all_ox_states_electronegativity_sorted[i];k<ksize;k++) {
          cce_vars.all_ox_states_electronegativity_sorted[i].push_back(element.oxidation_states[k]);
        }
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " num_all_ox_states_electronegativity_sorted[" << i << "]: " << cce_vars.num_all_ox_states_electronegativity_sorted[i] << endl;
          for (uint k=0,ksize=cce_vars.num_all_ox_states_electronegativity_sorted[i];k<ksize;k++) {
            cerr << __AFLOW_FUNC__ << " all oxidation state " << k << " of species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << "): " << cce_vars.all_ox_states_electronegativity_sorted[i][k] << endl;
          }
          cerr << endl;
        }
      } else{
        cce_vars.num_all_ox_states_electronegativity_sorted[i] = 0;
        cce_flags.flag("NO_OX_STATES",TRUE);
        message << "BAD NEWS: There are no known oxidation states for species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << ")."  << endl;
        message << "Therefore the oxidation states cannot be determined on this basis." << endl;
        if(cce_flags.flag("RUN_FULL_CCE")){
          message << "The formation enthalpy of this system is hence not correctable!"  << endl;
          message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
    }
  }

  //try_preferred_oxidation_states////////////////////////////////////////////////////////
  // try to determine the oxidation numbers using the preferred oxidation states for all cation species
  void try_preferred_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " Trying preferred oxidation numbers:" << endl;
    }
    // use ME's implementation of my algorithm to determine oxidation_numbers
    determine_cation_oxidation_states(structure, cce_vars, cce_vars.pref_ox_states_electronegativity_sorted);
    // print oxidation numbers and calculate sum
    cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars);
    // if sum of oxidation numbers is essentially zero, oxidation states should be regarded as determined correctly
    if (std::abs(cce_vars.oxidation_sum) <= DEFAULT_CCE_OX_TOL) {
      cce_flags.flag("OX_STATES_DETERMINED",TRUE);
    }
  }

  //treat_SbO2_special_case////////////////////////////////////////////////////////
  // for SbO2 the oxidation states are not identified properly
  // with the actual formula Sb2O4 it is a mixed valence oxide with one Sb+3 with 4 Sb-O bonds and one Sb+5 with 6 Sb-O bonds
  void treat_SbO2_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Sb") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Sb" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Sb = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Sb" ) {
      num_O_before_Sb=4;
    } else {
      num_O_before_Sb=0;
    }
    double amount_Sb=0;
    double amount_O=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "Sb"){
        amount_Sb=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of Sb ions= " << amount_Sb << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
    }
    double Sb_O_ratio = 0.0;
    if ( amount_O != 0 ){
      Sb_O_ratio=amount_Sb/amount_O;
    } else {
      message << " SbO2 special case. Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Sb_O_ratio,0.5) ){
      stringstream message;
      message << " This system is identified as a mixed valence compound." << endl;
      message << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
      message << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
      message << "Sb2O4 with ratio of Sb/O= " << Sb_O_ratio << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if ((structure.atoms.size() % 6) == 0) {
        uint num_formula_units_in_cell=structure.atoms.size()/6; // 6 for Sb2O4
        if(LDEBUG){
          cerr << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Sb"){
            // taking the first Sb ions as +3 (without knowing whether they are actually +3); works only because
            // the number of bonds for them will be adjusted to 4 as needed for Sb2O4 disregarding the actual
            // number of Sb-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Sb)*num_formula_units_in_cell ){  // (1 Sb3+ ions per formula unit + 4*O listed before in alphabetic order) * number of formula units
              cce_vars.oxidation_states[i]=3;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Sb+3 " << endl;
              }
              cce_vars.num_neighbors[i]=4; // for Sb2O4 the Sb3+ ions are 4-fold coordinated by oxygen
              if(LDEBUG){
                cerr << "Modified number of neighbors of Sb (atom[" << i << "]) taken as Sb+3: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=5;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Sb+5 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Sb2O4 the Sb5+ ions are 6-fold coordinated by oxygen
              if(LDEBUG){
                cerr << "Modified number of neighbors of Sb (atom[" << i << "]) taken as Sb+5: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
      } else {
        message << " The total number of atoms is not divisible by 6 as needed for Sb2O4. Please check your structure.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
    }
  }

  //treat_Pb3O4_special_case////////////////////////////////////////////////////////
  // for Pb3O4 the oxidation states are not identified properly
  // https://en.wikipedia.org/wiki/Lead(II,IV)_oxide
  void treat_Pb3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    stringstream message;
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Pb") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Pb" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Pb = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Pb" ) {
      num_O_before_Pb=4;
    } else {
      num_O_before_Pb=0;
    }
    double amount_Pb=0;
    double amount_O=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "Pb"){
        amount_Pb=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of Pb ions= " << amount_Pb << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
    }
    double Pb_O_ratio = 0.0;
    if ( amount_O != 0 ){
      Pb_O_ratio=amount_Pb/amount_O;
    } else {
      message << " Pb3O4 special case. Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Pb_O_ratio,0.75) ){
      stringstream message;
      message << " This system is identified as a mixed valence compound." << endl;
      message << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
      message << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
      message << "Pb3O4 with ratio of Pb/O= " << Pb_O_ratio << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if ((structure.atoms.size() % 7) == 0) {
        uint num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Pb3O4
        if(LDEBUG){
          cerr << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Pb"){
            // taking the first Pb ions as +4 (without knowing whether they are actually +4); works only because
            // the number of bonds for them will be adjusted to 6 as needed for Pb3O4 disregarding the actual
            // number of Pb-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Pb)*num_formula_units_in_cell ){  // (1 Pb4+ ions per formula unit + 4*O listed before in alphabetic order) * number of formula units
              cce_vars.oxidation_states[i]=4;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Pb+4 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Pb3O4 the Pb4+ ions are 6-fold coordinated by oxygen
              if(LDEBUG){
                cerr << "Modified number of neighbors of Pb (atom[" << i << "]) taken as Pb+4: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=2;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Pb+2 " << endl;
              }
              cce_vars.num_neighbors[i]=3; // for Pb3O4 the Pb2+ ions are 3-fold coordinated by oxygen
              if(LDEBUG){
                cerr << "Modified number of neighbors of Pb (atom[" << i << "]) taken as Pb+2: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
      } else {
        message << " The total number of atoms is not divisible by 7 as needed for Pb3O4. Please check your structure.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
    }
  }

  //treat_Ti_O_Magneli_phase_special_case////////////////////////////////////////////////////////
  // treat Ti-O Magneli phases; there are always 2xTi+3 per formula unit and the rest is Ti+4; 
  // fortunately, both Ti+3 and Ti+4 have 6 Ti-O bonds, hence one only needs to know how many ions 
  // of the respective oxidation state there are, not which one is which
  void treat_Ti_O_Magneli_phase_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Ti") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Ti" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Ti = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Ti" ) {
      num_O_before_Ti=5;
    } else {
      num_O_before_Ti=0;
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " Ti-O system, Magneli for Ti_nO_(2n-1), i.e. Ti-O ratio= " << 3.0/5 << ", " << 4.0/7 << ", " << 5.0/9 << ", " << 6.0/11 << ", " << 7.0/13 << ", " << 8.0/15 << ", " << 9.0/17 << ", " << 10.0/19 << "..." << endl;
    }
    double amount_O=0;
    double amount_Ti=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      } else if (structure.species[i] == "Ti"){
        amount_Ti=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
      cerr << __AFLOW_FUNC__ << " number of Ti ions= " << amount_Ti << endl;
    }
    double Ti_O_ratio = 0.0;
    if ( amount_O != 0 ){
      Ti_O_ratio=amount_Ti/amount_O;
    } else {
      message << " Ti-O Magneli phases special case. Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ratio of Ti/O= " << Ti_O_ratio << endl;
    }
    uint num_formula_units_in_cell = 0;
    // check for Magneli composition Ti_(n)O_(2n-1)
    double n = 0.0;
    bool magneli = false;
    for(n=3;n<101;n++){ // looping up to 101 should be enough
      if(LDEBUG){
        cerr << "n/(2*n-1)= " << n/(2*n-1) << endl;
      }
      if ( aurostd::isequal(Ti_O_ratio,n/(2*n-1)) ){
        stringstream message;
        message << " This system is identified as a mixed valence compound." << endl;
        message << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
        message << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
        message << "n= " << n << " Magneli composition Ti_nO_(2n-1)" << endl;
        if (!cce_flags.flag("UNIT_TEST")) {
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
        magneli = true;
        num_formula_units_in_cell=amount_Ti/n;
        if(LDEBUG){
          cerr << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
      }
    }
    if ( magneli == false){
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " Not a Magneli composition." << endl;
      }
    }
    if (magneli){
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Ti"){
          // taking the first Ti as +3 (without knowing whether they are actually +3) works only because
          // for the Magneli phases both Ti+3 and Ti+4 have both always 6 Ti-O bonds
          if ( i < (2+num_O_before_Ti)*num_formula_units_in_cell ){
            cce_vars.oxidation_states[i]=3;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+3 " << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=4;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+4 " << endl;
            }
          }
        }
      }
      check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
    }
  }

  //treat_Fe3O4_special_case////////////////////////////////////////////////////////
  // for Fe3O4 in inverse spinel structure the oxidation states are not identified properly
  // According to Wikipedia the Fe2+ ions are octahedrally coordinated and the Fe3+ ions are evenly distributed 
  // between the octahedral and tetrahedral sites and one hence needs per formula unit 6x the Fe2+ correction 
  // and 4+6=10x the Fe3+ correction, see:
  // https://en.wikipedia.org/wiki/Iron(II,III)_oxide
  // for Co3O4 and Mn3O4 this is different; in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different but at the present stage 
  // there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void treat_Fe3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Fe") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Fe" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Fe = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Fe" ) {
      num_O_before_Fe=4;
    } else {
      num_O_before_Fe=0;
    }
    double amount_Fe=0;
    double amount_O=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "Fe"){
        amount_Fe=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of Fe ions= " << amount_Fe << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
    }
    double Fe_O_ratio = 0.0;
    if ( amount_O != 0 ){
      Fe_O_ratio=amount_Fe/amount_O;
    } else {
      message << " Fe3O4 special case. Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Fe_O_ratio,0.75) ){
      stringstream message;
      message << " This system is identified as a mixed valence compound." << endl;
      message << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
      message << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
      message << "Fe3O4 with ratio of Fe/O= " << Fe_O_ratio << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if ((structure.atoms.size() % 7) == 0) {
        uint num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Fe3O4
        if(LDEBUG){
          cerr << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Fe"){
            // taking the first Fe ions as +2 (without knowing whether they are actually +2); works only because
            // the number of bonds for them will be adjusted to 6 (octahedral) as needed for Fe3O4 disregarding the actual
            // number of Fe-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Fe)*num_formula_units_in_cell ){  // 1 Fe2+ ions per formula unit * number of formula units
              cce_vars.oxidation_states[i]=2;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+2 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Fe3O4 the Fe2+ ions are 6-fold coordinated by oxygen according to Wikipedia
              if(LDEBUG){
                cerr << "Modified number of neighbors of Fe (atom[" << i << "]) taken as Fe+2: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=3;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
              }
              cce_vars.num_neighbors[i]=5; // for Fe3O4 the Fe3+ ions are evenly 6- and 4-fold, so on average 5-fold (set here as a hack) coordinated by oxygen according to Wikipedia
              if(LDEBUG){
                cerr << "Modified number of neighbors of Fe (atom[" << i << "]) taken as Fe+3 (average between even 6- and 4-fold coordination): " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
      } else {
        message << " The total number of atoms is not divisible by 7 as needed for Fe3O4. Please check your structure.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
    }
  }

  //treat_X3O4_special_case////////////////////////////////////////////////////////
  // for Co3O4 and Mn3O4 in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different than for Fe3O4
  // but at the present stage there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void treat_X3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, const string& cation_species, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == cation_species) || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == cation_species && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_cation_species = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == cation_species ) {
      num_O_before_cation_species=4;
    } else {
      num_O_before_cation_species=0;
    }
    double amount_cation_species=0;
    double amount_O=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == cation_species){
        amount_cation_species=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of " << cation_species << " ions = " << amount_cation_species << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
    }
    double cation_species_O_ratio = 0.0;
    if ( amount_O != 0 ){
      cation_species_O_ratio=amount_cation_species/amount_O;
    } else {
      message << " " << cation_species << "3O4 special case. Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(cation_species_O_ratio,0.75) ){
      stringstream message;
      message << " This system is identified as a mixed valence compound." << endl;
      message << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
      message << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
      message << cation_species << "3O4 with ratio of " << cation_species << "/O= " << cation_species_O_ratio << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if ((structure.atoms.size() % 7) == 0) {
        uint num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Mn3O4 and Co3O4
        if(LDEBUG){
          cerr << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == cation_species){
            // taking the first Mn/Co ions as +2 (without knowing whether they are actually +2); works only because
            // the number of bonds for them will be adjusted to 4 (tetrahedral) as needed for Mn3O4/Co3O4 disregarding the actual
            // number of Mn/Co-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_cation_species)*num_formula_units_in_cell ){  // 1 Mn/Co2+ ions per formula unit * number of formula units
              cce_vars.oxidation_states[i]=2;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to " << cation_species << "+2 " << endl;
              }
              cce_vars.num_neighbors[i]=4; // for Mn3O4/Co3O4 the Mn/Co2+ ions are 4-fold coordinated by oxygen according to Wikipedia
              if(LDEBUG){
                cerr << "Modified number of neighbors of " << cation_species << " (atom[" << i << "]) taken as " << cation_species << "+2: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=3;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to " << cation_species << "+3 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Mn3O4/Co3O4 the Mn/Co3+ ions are 6-fold coordinated by oxygen according to Wikipedia
              if(LDEBUG){
                cerr << "Modified number of neighbors of " << cation_species << " (atom[" << i << "]) taken as " << cation_species << "+3: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
      } else {
        message << " The total number of atoms is not divisible by 7 as needed for " << cation_species << "3O4. Please check your structure.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
    }
  }

  //treat_alkali_sesquioxide_special_case////////////////////////////////////////////////////////
  // for alkali metal sesquioxides (X2O3) the oxidation states and corrections might not be found correctly
  // since these systems are formulated to contain both per- and superoxide ions
  // see: https://en.wikipedia.org/wiki/Sesquioxide
  // the only known example to date for which also a structure was found on Springer Materials and the AFLOW ICSD
  // is Rb2O3 formulated as [(Rb+)4(O2−2)(O2-1)2] (for Cs4O6 there is only a single entry in a scanned pdf on Springer Materials)
  // this function will treat it as an exceptional case for all possible alkali metal sesquioxides
  // the oxidation numbers and per- as well as superoxide corrections will be adjusted accordingly
  void treat_alkali_sesquioxide_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    string alkali_metals = "Li,Na,K,Rb,Cs,Fr";
    vector<string> valkali_metals;
    aurostd::string2tokens(alkali_metals, valkali_metals, ",");
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && aurostd::WithinList(valkali_metals, KBIN::VASP_PseudoPotential_CleanName(structure.species[1]))) || (aurostd::WithinList(valkali_metals, KBIN::VASP_PseudoPotential_CleanName(structure.species[0])) && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;} // check whether it is a binary alkali metal oxide
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " This is a binary alkali metal oxide, checking whether it is an alkali metal sesquioxide..." << endl;
    }
    uint num_alkali_before_O = 0; // num cations before O not O before cations since setting oxidation states of anions below, not for cations as in other cases
    string alkali_metal = "";
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && aurostd::WithinList(valkali_metals, KBIN::VASP_PseudoPotential_CleanName(structure.species[1])) ) {
      num_alkali_before_O=0;
      alkali_metal=KBIN::VASP_PseudoPotential_CleanName(structure.species[1]);
    } else {
      num_alkali_before_O=2;
      alkali_metal=KBIN::VASP_PseudoPotential_CleanName(structure.species[0]);
    }
    double amount_O=0;
    double amount_alkali=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      } else {
        amount_alkali=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
      cerr << __AFLOW_FUNC__ << " number of alkali (" << alkali_metal << ") ions= " << amount_alkali << endl;
    }
    double O_alkali_ratio = 0.0;
    if ( amount_alkali != 0 ){
      O_alkali_ratio=amount_O/amount_alkali;
    } else {
      message << " Alkali metal sesquioxide special case. Amount of alkali atoms determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ratio of O/" << alkali_metal << "= " << O_alkali_ratio << endl;
    }
    // check for sesqui-composition alkali_metal2O3
    if ( aurostd::isequal(O_alkali_ratio,1.5) ){
      stringstream message;
      message << " This system is identified as an alkali metal sesquioxide (formally) containing both per- and superoxide ions." << endl;
      message << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
      message << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
      message << alkali_metal << "2O3 with ratio of O/" << alkali_metal << "= " << O_alkali_ratio << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
      if ((structure.atoms.size() % 5) == 0) {
        uint num_formula_units_in_cell=structure.atoms.size()/5; // 5 for alkali_metal2O3
        if(LDEBUG){
          cerr << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "O"){
            // taking the first O ions as -0.5 (superox. Os, without knowing whether they are actually -0.5); works only because
            // the number of superoxide bonds will be adjusted to 1 per formula unit as needed for alkali metal sesquioxides
            // taking last O ions as -1 (perox. Os, without knowing whether they are actually -1); works only because
            // the number of peroxide bonds will be adjusted to 0.5 per formula unit as needed for alkali metal sesquioxides
            // (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (2+num_alkali_before_O)*num_formula_units_in_cell ){  // 2 superoxide O-atoms per formula unit * number of formula units
              cce_vars.oxidation_states[i]=-0.5;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to -0.5 " << endl;
              }
              cce_vars.superox_indices[i]=1;
            } else {
              cce_vars.oxidation_states[i]=-1;
              if(LDEBUG){
                cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to -1 " << endl;
              }
              cce_vars.perox_indices[i]=1;
            }
          }
        }
        cce_vars.num_superox_bonds=1*num_formula_units_in_cell; // 1 superox. bonds per formula unit * num. formula units
        if(LDEBUG){
          cerr << "Number of superoxide bonds set to: " << cce_vars.num_superox_bonds << endl;
        }
        cce_vars.num_perox_bonds=0.5*num_formula_units_in_cell; // 0.5 perox. bonds per formula unit * num. formula units; this should neve yield a non-integer since there should be at least one complete peroxide ion per cell
        if(LDEBUG){
          cerr << "Number of peroxide bonds set to: " << cce_vars.num_perox_bonds << endl;
        }
      } else {
        message << " The total number of atoms is not divisible by 5 as needed for alkali metal sesquioxide. Please check your structure.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
    }
  }

  // following special cases only needed when determining oxidation states from Bader charges

  //treat_MnMoO4_special_case////////////////////////////////////////////////////////
  // for MnMoO4 the oxidation numbers are determined to be +4 for both Mn and Mo from the Bader charges; 
  // it should be Mn+2 & Mo+6; however since the sum of the falsely determined oxidation numbers 
  // is accidentally 0 it needs to be corrected individually
  void treat_MnMoO4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 3 && KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Mn" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Mo" && KBIN::VASP_PseudoPotential_CleanName(structure.species[2]) == "O" )) {return;}
    double amount_Mn=0;
    double amount_Mo=0;
    double amount_O=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "Mn"){
        amount_Mn=structure.num_each_type[i];
      } else if (structure.species[i] == "Mo"){
        amount_Mo=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of Mn ions= " << amount_Mn << endl;
      cerr << __AFLOW_FUNC__ << " number of Mo ions= " << amount_Mo << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
    }
    if (amount_Mo != 0 && amount_O != 0) {
      if (aurostd::isequal(amount_Mn/amount_Mo,1.0) && aurostd::isequal(amount_Mn/amount_O,0.25) && aurostd::isequal(amount_Mo/amount_O,0.25)) {
        stringstream message;
        message << " MnMoO4 special treatment since sum over oxidation states is zero but individual oxidation numbers are wrong!!!" << endl;
        if (!cce_flags.flag("UNIT_TEST")) {
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Mn"){
            cce_vars.oxidation_states[i]=2;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Mn+2 " << endl;
            }
          }
          if (structure.atoms[i].cleanname == "Mo"){
            cce_vars.oxidation_states[i]=6;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Mo+6 " << endl;
            }
          }
        }
        check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
      }
    } else {
      message << " MnMoO4 special case. Amount of Mo or O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
  }

  //treat_Ca2Fe2O5_CaFe2O4_LDA_special_case////////////////////////////////////////////////////////
  // for Ca2Fe2O5 and CaFe2O4 for LDA the oxidation numbers of Fe are not correctly determined 
  // to be Fe+3 but are partly Fe+2 which will be corrected here
  void treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 3 && KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Ca" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Fe" && KBIN::VASP_PseudoPotential_CleanName(structure.species[2]) == "O" )) {return;}
    double amount_Ca=0;
    double amount_Fe=0;
    double amount_O=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "Ca"){
        amount_Ca=structure.num_each_type[i];
      } else if (structure.species[i] == "Fe"){
        amount_Fe=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of Ca ions= " << amount_Ca << endl;
      cerr << __AFLOW_FUNC__ << " number of Fe ions= " << amount_Fe << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
    }
    //making sure it is Ca2Fe2O5
    if (amount_Fe != 0 && amount_O != 0) {
      if (aurostd::isequal(amount_Ca/amount_Fe,1.0) && aurostd::isequal(amount_Ca/amount_O,0.4) && aurostd::isequal(amount_Fe/amount_O,0.4)) {
        stringstream message;
        message << " Ca2Fe2O5 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
        if (!cce_flags.flag("UNIT_TEST")){
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Fe"){
            cce_vars.oxidation_states[i]=3;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
            }
          }
        }
        check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
        //making sure it is CaFe2O4
      } else if (amount_Ca/amount_Fe == 0.5 && amount_Ca/amount_O == 0.25 && amount_Fe/amount_O == 0.5) {
        stringstream message;
        message << " CaFe2O4 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
        if (!cce_flags.flag("UNIT_TEST")){
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Fe"){
            cce_vars.oxidation_states[i]=3;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
            }
          }
        }
        check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
      }
    } else {
      message << " Ca2Fe2O5/CaFe2O4 special case. Amount of Fe or O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
  }

  //treat_FeTiO3_LDA_special_case////////////////////////////////////////////////////////
  // for FeTiO3 for LDA the oxidation numbers of Ti are not correctly determined to be Ti+4 and using 
  // the general fixes would modify both the Ti AND the Fe oxidation numbers resulting again 
  // in non-zero oxidation number sum, which is fixed here
  void treat_FeTiO3_LDA_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if (! ( structure.species.size() == 3 && KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Fe" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[2]) == "Ti" )) {return;}
    double amount_Fe=0;
    double amount_O=0;
    double amount_Ti=0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){
      if (structure.species[i] == "Fe"){
        amount_Fe=structure.num_each_type[i];
      } else if (structure.species[i] == "O"){
        amount_O=structure.num_each_type[i];
      } else if (structure.species[i] == "Ti"){
        amount_Ti=structure.num_each_type[i];
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " number of Fe ions= " << amount_Fe << endl;
      cerr << __AFLOW_FUNC__ << " number of O ions= " << amount_O << endl;
      cerr << __AFLOW_FUNC__ << " number of Ti ions= " << amount_Ti << endl;
    }
    //making sure it is FeTiO3
    if (amount_Ti != 0 && amount_O != 0) {
      if (aurostd::isequal(amount_Fe/amount_Ti,1.0) && aurostd::isequal(amount_Fe/amount_O,1.0/3) && aurostd::isequal(amount_Ti/amount_O,1.0/3)) {
        stringstream message;
        message << " FeTiO3 special treatment for LDA since oxidation numbers for Ti, which should be Ti+4, are not correctly determined from Bader charges and using other fixing would also change Fe oxidation numbers!!!" << endl;
        if (!cce_flags.flag("UNIT_TEST")){
          ostream& oss = cout;
          ofstream FileMESSAGE;
          _aflags aflags;aflags.Directory=aurostd::getPWD();
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Ti"){
            cce_vars.oxidation_states[i]=4;
            if(LDEBUG){
              cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+4 " << endl;
            }
          }
        }
        check_ox_nums_special_case(structure, cce_flags, cce_vars, oss);
      }
    } else {
      message << " FeTiO3 special case. Amount of Ti or O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
  }

  //check_ox_nums_special_case////////////////////////////////////////////////////////
  // part checking oxidation numbers for special case treatment
  void check_ox_nums_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    stringstream message;
    // calculate sum of oxidation numbers
    cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars);
    // system should not be regarded correctable if sum over oxidation states is not zero
    if (std::abs(cce_vars.oxidation_sum) > DEFAULT_CCE_OX_TOL) { // this case should never occur since set oxidation states should always work
      cce_flags.flag("CORRECTABLE",FALSE);
      oss << print_output_oxidation_numbers(structure, cce_vars); // print previously gathered output e.g. from determination of oxidation numbers
      message << "BAD NEWS: The formation enthalpy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
      message << "Sum over all oxidation numbers is: " << cce_vars.oxidation_sum << endl;
      message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    } else {
      cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
    }
  }

  //try_all_oxidation_states////////////////////////////////////////////////////////
  // try to determine the oxidation numbers using all known oxidation states for all cation species
  void try_all_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " Trying all known oxidation numbers:" << endl;
    }
    // use ME's implementation of my algorithm to determine oxidation_numbers
    determine_cation_oxidation_states(structure, cce_vars, cce_vars.all_ox_states_electronegativity_sorted);
  }

  //determine_cation_oxidation_states////////////////////////////////////////////////////////
  //ME20191101
  // for avoiding recursion algorithm to determine cation oxidation_numbers
  // determine the cation oxidation numbers by using possible oxidation states for each species 
  // which can be either the preferred or all known oxidation numbers
  void determine_cation_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars, const vector<vector<double> >& possible_ox_states) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    uint num_cation_species = cce_vars.cations_map.size(); // the number of cation species
    uint natoms = cce_vars.oxidation_states.size();

    // Initialize
    for (uint i = 0; i < num_cation_species; i++) { // loop over all cations
      for (uint j = 0; j < cce_vars.cations_map[i].size(); j++) { // loop over atoms of the ith cation type that is given by the second index of cation_map
        if (cce_vars.multi_anion_atoms[cce_vars.cations_map[i][j]] != 1){ // exclude atoms that have been identified as multi anion atoms previously
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " Updating oxidation number for atom " << cce_vars.cations_map[i][j] << " (" << structure.atoms[cce_vars.cations_map[i][j]].cleanname << ") to " << possible_ox_states[i][0] << endl;
          }
          cce_vars.oxidation_states[cce_vars.cations_map[i][j]] = possible_ox_states[i][0]; // possible_ox_states (either preferred or all) for all species should be electronegativity sorted
        }
      }
    }
    if(LDEBUG){
      for(uint n=0,nsize=structure.atoms.size();n<nsize;n++){ 
        cerr << __AFLOW_FUNC__ << " chosen oxidation state for atom " << n << " (" << structure.atoms[n].cleanname << "): " << cce_vars.oxidation_states[n] << endl;
      }
    }
    // check
    double total_ox = 0.0;
    for (uint at = 0; at < natoms; at++) total_ox += cce_vars.oxidation_states[at];
    bool iszero = aurostd::isequal(total_ox, 0.0, DEFAULT_CCE_OX_TOL);
    // indicate whether these oxidation states were successful or not
    if (iszero) {
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " Oxidation numbers successfully determined since oxidation sum is " << total_ox << endl;
      }
    } else{
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " No successful determination of oxidation numbers since oxidation sum is " << total_ox << endl;
      }
    }

    // Cycle through oxidation states
    if (!iszero) { // maybe the oxidation numbers from the initialization work already, then nothing else needs to be done
      uint maxsize = 0; // maximum number of preferred or all oxidation states for a cation species in the system
      for (uint j = 0; j < num_cation_species; j++) {
        if (possible_ox_states[j].size() > maxsize) maxsize = possible_ox_states[j].size();
      }

      uint num_ox_states = 0;
      uint i = 1; // i: number of tested ox state; can start at 1 since 0th pref./all ox states has been used during initialization
      while (!iszero && (i < maxsize)) {
        for (uint j = num_cation_species; j > 0; j--) { // have to loop downward since species are sorted ascending by electronegativity and oxidation number should be changed for less electronegative species first
          num_ox_states = possible_ox_states[j-1].size(); // j-1 since species loop goes downward and j_min should be comapred to be > 0 for uint, how many pref./all ox. states there are for the jth cation species
          if (i < num_ox_states) { // the index i looping over all possible ox states for species j must be smaller than the maximum number of pref./all ox states for this species 
            for (uint k = 0; k < cce_vars.cations_map[j-1].size(); k++) {
              if (cce_vars.multi_anion_atoms[cce_vars.cations_map[j-1][k]] != 1){ // exclude atoms that have been identified as multi anion atoms previously
                if(LDEBUG){
                  cerr << __AFLOW_FUNC__ << " Updating oxidation number for atom " << cce_vars.cations_map[j-1][k] << "(" << structure.atoms[cce_vars.cations_map[j-1][k]].cleanname << ") to " << possible_ox_states[j-1][i] << endl;
                }
                cce_vars.oxidation_states[cce_vars.cations_map[j-1][k]] = possible_ox_states[j-1][i]; // i goes over all preferred/all oxidation states
              }
            }
            if(LDEBUG){
              for(uint n=0,nsize=structure.atoms.size();n<nsize;n++){ 
                cerr << __AFLOW_FUNC__ << " chosen oxidation state for atom " << n << " (" << structure.atoms[n].cleanname << "): " << cce_vars.oxidation_states[n] << endl;
              }
            }
            // check
            total_ox = 0.0;
            for (uint at = 0; at < natoms; at++) total_ox += cce_vars.oxidation_states[at];
            iszero = aurostd::isequal(total_ox, 0.0, DEFAULT_CCE_OX_TOL);
            // indicate whether these oxidation states were successful or not
            if (iszero) {
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " Oxidation numbers successfully determined since oxidation sum is " << total_ox << endl;
              }
              break; // BREAK for loop over possible oxidation states k for species j-1 upon success (sum over oxidation numbers is zero)
            } else{
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " No successful determination of oxidation numbers since oxidation sum is " << total_ox << endl;
              }
            }
          }
        }
        i++;
      }
    }
  }

  //get_oxidation_states_sum////////////////////////////////////////////////////////
  // Calculate and return the sum of all oxidation states.
  double get_oxidation_states_sum(CCE_Variables& cce_vars) {
    cce_vars.oxidation_sum=0;
    for (uint k=0,ksize=cce_vars.oxidation_states.size();k<ksize;k++){
      cce_vars.oxidation_sum += cce_vars.oxidation_states[k];
    }
    return cce_vars.oxidation_sum;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //              DETERMINE OXIDATION STATES FROM BADER CHARGES              //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // determine oxidation numbers based on Bader charges (outdated)
  // functions for treating known special cases not included for the electronegativity based implementation are included here

  //get_oxidation_states_from_Bader////////////////////////////////////////////////////////
  // determine the oxidation numbers of the ions by an analysis of the Bader charges of the system and handle known exceptional cases explicitly
  vector<double> get_oxidation_states_from_Bader(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, const string& directory_path, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " DETERMINATION OF OXIDATION NUMBERS FROM BADER CHARGES:" << endl;
    }
    if (cce_vars.anion_species == "O") {
      // in this part functional dependent settings will be primarily evaluated with "if, else if" conditions, 
      // since the calculation (for getting the Bader charges) could have only been done with one functional and not more than one
      // oxidation numbers should not be functional dependent since results are for a specific calculation, i.e. PBE or LDA or SCAN or PBE+U:ICSD and not more than one
      // for a given structure there should be only one correct assignment of oxidation numbers
      string system_name = "";
      string functional = "";
      // check whether aflow.in exists
      if (aurostd::FileExist(directory_path + "/" + _AFLOWIN_)) {
        // get system name to identify Bader charges file and functional to distinguish corrections to be loaded from aflow.in
        get_system_name_functional_from_aflow_in(structure, cce_flags, cce_vars, system_name, functional, directory_path);
        if (cce_flags.flag("CORRECTABLE")){
          // check whether Bader file exists
          string Bader_file = directory_path + "/" + system_name + "_abader.out";
          if (aurostd::FileExist(Bader_file) || aurostd::EFileExist(Bader_file)) {
            // determine Bader charges from Bader file and store in array/vector Bader_charges; O Bader charges will be included although they might not be used
            cce_vars.Bader_charges = get_Bader_charges_from_Bader_file(structure, cce_vars, Bader_file);
            // determine oxidation numbers from Bader charges
            cce_vars.oxidation_states = Bader_charges_to_oxidation_states(structure, cce_flags, cce_vars, functional);
            // check whether sum of oxidation numbers is equal zero and try to correct for known special cases (for some cases even if it is zero)
            if (cce_flags.flag("CORRECTABLE")){
              // calculate sum of oxidation numbers
              cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars);
              // SECOND point to deal with special cases known from (binary+ternary) oxides
              // Ti-O Magneli phases need to be treated specially since oxidation numbers are not recognized appropriately for all functionals
              cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
              treat_Ti_O_Magneli_phase_special_case(structure, cce_flags, cce_vars);
              // for Fe3O4 in inverse spinel structure, the oxidation states are not identified properly via the Bader charges.
              treat_Fe3O4_special_case(structure, cce_flags, cce_vars);
              // for Mn3O4 in spinel structure, the oxidation states are not identified properly
              treat_X3O4_special_case(structure, cce_flags, cce_vars, "Mn");
              // for Co3O4 in spinel structure, the oxidation states are not identified properly
              treat_X3O4_special_case(structure, cce_flags, cce_vars, "Co");
              // MnMoO4 needs to be treated specially since oxidation numbers are not recognized appropriately
              treat_MnMoO4_special_case(structure, cce_flags, cce_vars);
              if (functional == "LDA") {
                // Ca2Fe2O5 & CaFe2O4 need to be treated specially for LDA since oxidation numbers are not recognized appropriately
                treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(structure, cce_flags, cce_vars);
                // FeTiO3 needs to be treated specially for LDA since oxidation numbers are not recognized appropriately
                treat_FeTiO3_LDA_special_case(structure, cce_flags, cce_vars);
              }
              if (std::abs(cce_vars.oxidation_sum) > DEFAULT_CCE_OX_TOL) {
                // general scheme to repair wrong oxidation numbers based on changing oxidation state of several ions known to be problematic
                general_attempt_fixing_oxidation_states(structure, cce_flags, cce_vars);
                // calculate sum of oxidation numbers
                cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars);
                // system should not be regarded correctable if sum over oxidation states is not zero
                if (std::abs(cce_vars.oxidation_sum) > DEFAULT_CCE_OX_TOL) {
                  cce_flags.flag("CORRECTABLE",FALSE);
                  oss << print_output_oxidation_numbers(structure, cce_vars); // print previously gathered output e.g. from determination of oxidation numbers
                  message << "BAD NEWS: The formation enthalpy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
                  message << "Sum over all oxidation numbers is: " << cce_vars.oxidation_sum << endl;
                  message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
                  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
                }
              }
            }
          } else {
            cce_flags.flag("CORRECTABLE",FALSE);
            message << "Bader file " << Bader_file << " (or xz, bz2, gz version) not found. A Bader file is required to determine the oxidation numbers." << endl;
            message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
          }
        }
      } else {
        cce_flags.flag("CORRECTABLE",FALSE);
        message << "aflow.in file not found. An aflow.in file is required to identify the functional and to find the Bader charges file to determine the oxidation numbers." << endl;
        message << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
    } else {
      message << " BAD NEWS: The determination of oxidation numbers from Bader charges works currently only for oxides. Here the anion is determined to be " << cce_vars.anion_species << ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    return cce_vars.oxidation_states;
  }

  //get_system_name_functional_from_aflow_in////////////////////////////////////////////////////////
  // determine the system name and functional from the aflow_in
  void get_system_name_functional_from_aflow_in(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, string& system_name, string& functional, const string& directory_path) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    string aflowin_file= directory_path + "/" + _AFLOWIN_;
    string outcar_file= directory_path + "/" + "OUTCAR.relax1";
    functional=get_functional_from_aflow_in_outcar(structure, aflowin_file, outcar_file);
    if (functional.empty()) {
      message << " Functional cannot be determined from aflow.in. Corrections are available for PBE, LDA, SCAN, or PBE+U:ICSD.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    // Bader implementation works currently only with Bader file and aflow.in in directory where AFLOW is run
    system_name=KBIN::ExtractSystemName(directory_path);  //CO20200928
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " functional determined from aflow.in: " << functional << endl;
      cerr << __AFLOW_FUNC__ << " PBE: " << aurostd::WithinList(cce_vars.vfunctionals, "PBE") << endl;
      cerr << __AFLOW_FUNC__ << " LDA: " << aurostd::WithinList(cce_vars.vfunctionals, "LDA") << endl;
      cerr << __AFLOW_FUNC__ << " SCAN: " << aurostd::WithinList(cce_vars.vfunctionals, "SCAN") << endl;
      cerr << __AFLOW_FUNC__ << " PBE+U:ICSD: " << aurostd::WithinList(cce_vars.vfunctionals, "PBE+U:ICSD") << endl;
      cerr << __AFLOW_FUNC__ << " exp: " << aurostd::WithinList(cce_vars.vfunctionals, "exp") << endl;
    }
    // if functional determined from aflow.in is different from the ones given by the input options, 
    // throw warning that oxidation numbers are only determined on the basis of a specific functional
    // also when only the formation enthalpy from the exprimental values per bond shall be obtained,
    // a warning should be given from which functional the data are used to determine the ox nums.
    if(!(cce_vars.vfunctionals.size() == 1 && cce_vars.vfunctionals[0] == functional)){
      stringstream message;
      message << " The oxidation numbers are only determined on the basis of a" << (functional == "LDA"?"n ":" ") << functional << " calculation." << endl;
      if (!cce_flags.flag("UNIT_TEST")) {
        ostream& oss = cout;
        ofstream FileMESSAGE;
        _aflags aflags;aflags.Directory=aurostd::getPWD();
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
    }
  }

  //get_Bader_charges_from_Bader_file////////////////////////////////////////////////////////
  // determine the Bader charges from the AFLOW Bader file (_abader.out)
  vector<double> get_Bader_charges_from_Bader_file(const xstructure& structure, CCE_Variables& cce_vars, const string& Bader_file) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    vector<string> vlines; 
    aurostd::efile2vectorstring(Bader_file,vlines);
    vector<string> tokens;
    for (uint i = 1; i < vlines.size(); i++) {
      aurostd::string2tokens(vlines[i], tokens, " ");
      cce_vars.Bader_charges.push_back(aurostd::string2utype<double>(tokens[2]));
    }
    for (uint k=0,ksize=structure.atoms.size();k<ksize;k++){ // there should always be as many Bader charges as atoms in the structure
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " Bader_charges[" << k << "]: " << cce_vars.Bader_charges[k] << endl;
      }
    }
    return cce_vars.Bader_charges;
  }

  //Bader_charges_to_oxidation_states////////////////////////////////////////////////////////
  // compare Bader charges to template values from fitting set and set oxidation numbers accordingly
  vector<double> Bader_charges_to_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, string& functional, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    bool print_empty_line=TRUE;
    bool error1=FALSE;
    bool error2=FALSE;
    vector<string> species_missing_corrections;
    vector<string> undetermined_ox_states;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname != cce_vars.anion_species){
        string Bader_templ_line = "";
        if (get_Bader_templates(structure.atoms[i].cleanname) == "") {
          cce_flags.flag("CORRECTABLE",FALSE);
          error1=TRUE;
          if (print_empty_line) { // print only one empty line before error blog
            oss << endl;
            print_empty_line=FALSE; // previously gathered output should only be printed once
          }
          oss << "VERY BAD NEWS: There is no correction for " << structure.atoms[i].cleanname << " (ATOM[" << i << "])" << " coordinated by " << cce_vars.anion_species << " since this species was not included in the set for deducing corrections!"  << endl;
          oss << endl;
          if (species_missing_corrections.size() == 0) {
            species_missing_corrections.push_back(structure.atoms[i].cleanname);
          } else if (species_missing_corrections[species_missing_corrections.size()-1] != structure.atoms[i].cleanname) {
            species_missing_corrections.push_back(structure.atoms[i].cleanname);
          }
        } else {
          Bader_templ_line=get_Bader_templates(structure.atoms[i].cleanname);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " Bader templates: " << Bader_templ_line << endl;
          }
          vector<string> Bader_tokens;
          aurostd::string2tokens(Bader_templ_line, Bader_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. Bader_tokens[1] is not a space but the oxidation number
          uint num_ox_states = aurostd::string2utype<uint>(Bader_tokens[0]);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " number of oxidation states for which Bader charges are available: " << num_ox_states << endl;
          }
          double Bader_template=0; // Bader cation charge from the binary oxide from which the correction was deduced
          double Bader_tolerance=0.26; // tolerance with which the Bader charge of each atom is allowed to deviate from any of the template values to still assign the correction (and oxidation number)
          double Bader_deviation=0.5; // deviation of the Bader charge of the cation from the template value; initial value safely larger than Bader tolerance for identifying oxidation number so that if value is not changed, no correction is made
          double Bader_deviation_0=Bader_deviation; // initial Bader deviation for later check whether any corrections were found, i. e. Bader deviation was changed
          double Bader_deviation_min=17; // initialize with high value; used later to output smallest Bader deviation if no correction according to the Bader tolerance was identified
          for(uint n=0;n<num_ox_states;n++){ //loop over all oxidation number for which Bader charges are available and read the Bader charges from the respective positions
            // here the functional determined from the aflow.in is decisive; exp should not occur
            vector<string> vfunctional_aflow_in;
            vfunctional_aflow_in.push_back(functional);
            if (get_offset(functional) == -1) {
              message << " Unknown functional " << functional << ". Please choose PBE, LDA, or SCAN.";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
            }
            int offset=get_offset(functional); // define also local offset variable since offset must correspond to functional detected from aflow.in
            for (uint j = 0; j < vfunctional_aflow_in.size(); j++) {
              Bader_template = aurostd::string2utype<double>(Bader_tokens[offset/2+2+(CCE_num_functionals_Bader+1)*n]);
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " Bader_template " << vfunctional_aflow_in[j] << ": " << Bader_template << endl;
              }
            }
            if ( std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_tolerance && std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_deviation ){ // must be compared to Bader_deviation to load correction for Bader charge closest to template and not to be overloaded by other corrections for later tested Bader charges (oxidation states)
              Bader_deviation= std::abs(Bader_template-cce_vars.Bader_charges[i]);
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " Bader_deviation: " << Bader_deviation << endl;
              }
              cce_vars.oxidation_states[i]= aurostd::string2utype<double>(Bader_tokens[1+(CCE_num_functionals_Bader+1)*n]); // (CCE_num_functionals_Bader+1)*n because oxidation numbers in Bader_templ_line are separated by corrections for four different functionals (PBE, LDA, SCAN, and PBE+U:ICSD) and there is 1 formal oxidation number in front of the Bader charges
            } else { // only if element is found but no corrections can be assigned because above conditions are not met, update Bader_deviation_min
              if ( std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_deviation_min ) {
                Bader_deviation_min= std::abs(Bader_template-cce_vars.Bader_charges[i]);
              }
            }
          }
          if ( Bader_deviation == Bader_deviation_0 ){ // only if no oxidation state was found, i. e. Bader deviation was not changed but there are corrections (Bader charges) for this species (get_Bader_templates returns non-empty string), display the following error message after handling the special cases
            // FIRST point to deal with special cases known from (ternary) oxides; other cases will be dealt with when looking at the oxidation numbers and checking the sum
            // since it can happen that the Bader charge for W in the compound (especially for W+6) is too far away 
            // from the Bader charge template, here W+6 is set and later the sum over the oxidation states will still be checked
            if (structure.atoms[i].cleanname == "W") {
              cce_vars.oxidation_states[i]=6;
              // since it is possible that the Bader charge for Pb in the compound (especially for Pb+2) is too far away 
              // from the Bader charge template, here Pb+2 is set and later the sum over the oxidation states will still be checked
            } else if (structure.atoms[i].cleanname == "Pb") {
              cce_vars.oxidation_states[i]=2;
              // since it can happen that the Bader charge for Ag in the compound (especially for Ag+1) is too far away 
              // from the Bader charge template, here Ag+1 is set and later the sum over the oxidation states will still be checked
            } else if (structure.atoms[i].cleanname == "Ag") {
              cce_vars.oxidation_states[i]=1;
            } else {
              cce_flags.flag("CORRECTABLE",FALSE);
              error2=TRUE;
              if (print_empty_line) { // print only one empty line before error blog
                oss << endl;
                print_empty_line=FALSE; // previously gathered output should only be printed once
              }
              undetermined_ox_states.push_back(structure.atoms[i].cleanname + " (ATOM[" + aurostd::utype2string<uint>(i) + "])");
              oss << "BAD NEWS: The oxidation number (and hence the correction) for " << structure.atoms[i].cleanname << " (ATOM[" << i << "])" << " cannot be identified from the Bader charges!"  << endl;
              oss << "The deviation of the Bader charge from the closest tested template value is: " << Bader_deviation_min << " electrons. This is larger than the tolerance: " << Bader_tolerance  << " electrons." << endl;
              // list all oxidation states of the element for which corrections are available
              string ox_nums_avail="";
              string separator=", ";
              vector<string> ox_nums_avail_vec; //RF20200826
              for(uint n=0;n<num_ox_states;n++){ 
                ox_nums_avail_vec.push_back(Bader_tokens[1+(CCE_num_functionals_Bader+1)*n]); // (CCE_num_functionals_Bader+1)*n because oxidation numbers in Bader_templ_line are separated by corrections for four different functionals (PBE, LDA, SCAN, and PBE+U:ICSD) and there is 1 formal oxidation number in front of the Bader charges  //RF20200826
                //[RF20200826 - OBSOLETE]if (n<num_ox_states-1){
                //[RF20200826 - OBSOLETE]  ox_nums_avail+= Bader_tokens[1+(CCE_num_functionals_Bader+1)*n] + separator ; // (CCE_num_functionals_Bader+1)*n because oxidation numbers in Bader_templ_line are separated by corrections for four different functionals (PBE, LDA, SCAN, and PBE+U:ICSD) and there is 1 formal oxidation number in front of the Bader charges
                //[RF20200826 - OBSOLETE]} else if (n==num_ox_states-1){
                //[RF20200826 - OBSOLETE]  ox_nums_avail+= Bader_tokens[1+(CCE_num_functionals_Bader+1)*n]; // (CCE_num_functionals_Bader+1)*n because oxidation numbers in Bader_templ_line are separated by corrections for four different functionals (PBE, LDA, SCAN, and PBE+U:ICSD) and there is 1 formal oxidation number in front of the Bader charges
                //[RF20200826 - OBSOLETE]}
              }
              ox_nums_avail = aurostd::joinWDelimiter(ox_nums_avail_vec, separator);
              oss << "Corrections for " << structure.atoms[i].cleanname << " coordinated by " << cce_vars.anion_species << " are available for oxidation states: " << ox_nums_avail << endl;
              oss << "If the desired oxidation state is listed but it is just not correctly determined from the Bader charges," << endl;
              oss << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
              oss << endl;
            }
          }
        }
      } else if (structure.atoms[i].cleanname == cce_vars.anion_species) {
        cce_vars.oxidation_states[i] = cce_vars.standard_anion_charge; // oxidation numbers for O are assumed to be -2 (other anions accordingly) and are corrected below if it is a per- or superoxide O atom as identified from the structural analysis in other functions
        if (cce_vars.num_perox_bonds > 0){
          for (uint j=0,jsize=structure.atoms.size();j<jsize;j++){
            if (cce_vars.perox_indices[j]==1 && j == i){
              cce_vars.oxidation_states[i]=-1;
            }
          }
        }
        if (cce_vars.num_superox_bonds > 0){
          for (uint j=0,jsize=structure.atoms.size();j<jsize;j++){
            if (cce_vars.superox_indices[j]==1 && j == i){
              cce_vars.oxidation_states[i]=-0.5;
            }
          }
        }
      }
    }
    if (error1) { // errors can only be thrown after loop over atoms is complete since the output should indicate all species for which corrections might be missing/cannot be identified
      message << " VERY BAD NEWS: The formation enthalpy of this system is not correctable since there are no corrections for " << aurostd::joinWDelimiter(species_missing_corrections, ", ") << " coordinated by " << cce_vars.anion_species << "!" << endl;
      message << " See also the output for details.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    } else if (error2) {
      message << " BAD NEWS: The oxidation numbers (and hence the corrections) of " << aurostd::joinWDelimiter(undetermined_ox_states, ", ") << " cannot be identified from the Bader charges!"  << endl;
      message << " See the output for details.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
    return cce_vars.oxidation_states;
  }

  //general_attempt_fixing_oxidation_states////////////////////////////////////////////////////////
  // for some ions known to be problematic for the Bader analysis (e.g. V+5, Fe+2, Fe+3, Ti+4) 
  // a brute force ansatz can be implemented to just change their oxidation state and later check 
  // whether it fixes the oxidation number sum rule
  void general_attempt_fixing_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    message << " The sum over all oxidation numbers (as determined from the Bader charges) is NOT zero, trying to repair that based on known problematic cases (Ti, V, Fe). This may or may not work." << endl;
    if (!cce_flags.flag("UNIT_TEST")) {
      ostream& oss = cout;
      ofstream FileMESSAGE;
      _aflags aflags;aflags.Directory=aurostd::getPWD();
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
    }
    // repairing only by considering non mixed valence oxides
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ 
      if (structure.atoms[i].cleanname == "Ti"){
        if ( cce_vars.oxidation_states[i] == 3){
          cce_vars.oxidation_states[i]=4;
          if(LDEBUG){
            cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+4 " << endl;
          }
        }
      }
      if (structure.atoms[i].cleanname == "Fe"){
        if ( cce_vars.oxidation_states[i] == 2){
          cce_vars.oxidation_states[i]=3;
          if(LDEBUG){
            cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
          }
        } else if ( cce_vars.oxidation_states[i] == 3){
          cce_vars.oxidation_states[i]=2;
          if(LDEBUG){
            cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+2 " << endl;
          }
        }
      }
      if (structure.atoms[i].cleanname == "V"){
        if ( cce_vars.oxidation_states[i] == 4){
          cce_vars.oxidation_states[i]=5;
          if(LDEBUG){
            cerr << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to V+5 " << endl;
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                            ASSIGN CORRECTIONS                           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //get_corrections////////////////////////////////////////////////////////
  // determine the corrections after oxidation numbers are known (from input, EN, or Bader)
  void get_corrections(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, const string& considered_anion_species, const vector<uint>& num_neighbors, vector<vector<double> >& corrections_atom, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    bool error=FALSE;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << "  ASSIGNMENT OF CORRECTIONS:" << endl;
    }
    bool print_previous_output=TRUE;
    vector<string> missing_corrections;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if ((structure.atoms[i].cleanname != cce_vars.anion_species) && (cce_vars.multi_anion_atoms[i] != 1)){ // exclude main anion species and multi anion atoms detected previously
        string corrections_line = "";
        string cor_identifier = structure.atoms[i].cleanname + "_+" + aurostd::utype2string<double>(cce_vars.oxidation_states[i]) + "_" + considered_anion_species;
        if (num_neighbors[i] > 0){ // are there actually bonds between the cation and the anion species under consideration (to load only needed corrections in multi-anion systems)
          if ( get_corrections_line(considered_anion_species, cor_identifier) == "") { // the considered anion species can be the main anion species or a multi anion species
            cce_flags.flag("CORRECTABLE",FALSE);
            error=TRUE;
            if (print_previous_output && cce_vars.anion_species == considered_anion_species) { //second condition should make sure that for multi-anion systems oxidationnnumbers are only printed once
              oss << print_output_oxidation_numbers(structure, cce_vars);
              print_previous_output=FALSE; // previously gathered output should only be printed once
            }
            string info_missing_corrections=structure.atoms[i].cleanname + " in oxidation state " + "+" + aurostd::utype2string<double>(cce_vars.oxidation_states[i]) + " coordinated by " + considered_anion_species;
            oss << "BAD NEWS: No correction available for " << structure.atoms[i].cleanname << " (atom " << i+1 << ")" << " in oxidation state +" << cce_vars.oxidation_states[i] << " when coordinated by " << considered_anion_species << "." << endl; // i+1 to convert to human based counting from 1
            if (missing_corrections.size() == 0) {
              missing_corrections.push_back(info_missing_corrections);
            } else if (missing_corrections[missing_corrections.size()-1] != info_missing_corrections) {
              missing_corrections.push_back(info_missing_corrections);
            }
            //checking for which oxidation states corrections are available and throw out errors accordingly
            uint ox_nums_count=0;
            vector<uint> ox_nums_avail_vec;
            for(uint k=0;k<12;k++){ // larger than ox. num. +12 should not occur
              if ( get_corrections_line(considered_anion_species, structure.atoms[i].cleanname + "_+" + aurostd::utype2string<uint>(k) + "_" + considered_anion_species) != "") {
                ox_nums_count+=1;
                ox_nums_avail_vec.push_back(k);
              }
            }
            if ( ox_nums_count == 0) {
              oss << "Currently no corrections available for " << structure.atoms[i].cleanname << " when coordinated by "<< considered_anion_species << "." << endl;
            } else {
              // list all oxidation states of the element for which corrections are available
              string ox_nums_avail="";
              string separator=", +";
              ox_nums_avail = aurostd::joinWDelimiter(ox_nums_avail_vec, separator);
              oss << "Corrections for " << structure.atoms[i].cleanname << " coordinated by " << considered_anion_species << " are available for oxidation states: " << "+" << ox_nums_avail << endl;
              oss << "If the desired oxidation state is listed but it is just not correctly determined," << endl;
              oss << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
            }
            oss << endl;
          } else {
            corrections_line=get_corrections_line(considered_anion_species, cor_identifier);
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " corrections line: " << corrections_line << endl;
            }
            if ( corrections_line.find("*") != string::npos ) {
              message << " The correction for " << structure.atoms[i].cleanname << " (atom " << i+1 << ")" << " in oxidation state +" << cce_vars.oxidation_states[i] << " when coordinated by " << considered_anion_species << " might be less accurate since it was obtained from less well validated experimental data as the other corrections!";
              pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, oss, _LOGGER_WARNING_,XHOST.vflag_control.flag("PRINT_MODE::JSON")); //CO20210623 - silent with json output
            }
            if ( corrections_line.find("^") != string::npos ) {
              message << " The correction for " << structure.atoms[i].cleanname << " (atom " << i+1 << ")" << " in oxidation state +" << cce_vars.oxidation_states[i] << " when coordinated by " << considered_anion_species << " is only approximate since in this case the explicit oxidation state dependence was lifted!";
              pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__, message, oss, _LOGGER_WARNING_,XHOST.vflag_control.flag("PRINT_MODE::JSON"));  //CO20210623 - silent with json output
            }
            // load cation corrections
            load_cation_corrections(structure, cce_vars, corrections_line, corrections_atom, i);
          }
        }
      } else if ((structure.atoms[i].cleanname == cce_vars.anion_species) || (cce_vars.multi_anion_atoms[i] == 1)) {
        // set anion corrections (to zero)
        set_anion_corrections(structure, cce_vars, corrections_atom, i);
      }
    }
    if(error){
      // errors can only be thrown after loop over atoms is complete since the output should indicate all species for which corrections might be missing/cannot be identified
      message << " BAD NEWS: No correction available for " << std::showpos << aurostd::joinWDelimiter(missing_corrections, ", ") << "." << endl;
      message << " See also the output for details.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
  }

  //load_cation_corrections////////////////////////////////////////////////////////
  // load corrections for the cations from the corrections line according to the functional
  // here for every functional there will be a separate "if" evaluation, since when only a structure (+ oxidation numbers) 
  // are provided as input, one might want to get the corrections for more than one functional
  void load_cation_corrections(const xstructure& structure, CCE_Variables& cce_vars, const string& corrections_line, vector<vector<double> >& corrections_atom, uint i) { // also provide increment for loop of get_corrections function
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    vector<string> corrections_tokens;
    aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        for (uint l = 0; l < num_temps; l++) {
          corrections_atom[num_temps*k+l][i]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+l+2]); // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for " << cce_vars.vtemperatures[l] << "K: " << corrections_atom[num_temps*k+l][i] << " eV/bond" << endl;
          }
        }
      } else {
        corrections_atom[num_temps*k][i]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for 298.15K: " << corrections_atom[num_temps*k][i] << " eV/bond" << endl;
        }
      }
    }
  }

  //set_anion_corrections////////////////////////////////////////////////////////
  // set corrections for anion atoms to zero since only for cations number of bonds with anions are needed; 
  // per- & superoxides will be dealt with in other function
  void set_anion_corrections(const xstructure& structure, CCE_Variables& cce_vars, vector<vector<double> >& corrections_atom, uint i) { // also provide increment for loop of get_corrections function
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        for (uint l = 0; l < num_temps; l++) {
          corrections_atom[num_temps*k+l][i]=0; // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for " << cce_vars.vtemperatures[l] << "K: " << corrections_atom[num_temps*k+l][i] << " eV/bond" << endl;
          }
        }
      } else {
        corrections_atom[num_temps*k][i]=0;
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for 298.15K: " << corrections_atom[num_temps*k][i] << " eV/bond" << endl;
        }
      }
    }
  }

  //check_get_per_super_ox_corrections////////////////////////////////////////////////////////
  // check whether corrections for per- or superoxide ions are needed and assign accordingly
  void check_get_per_super_ox_corrections(CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    // check for per- and superoxides should always be kept separate since a compound could have both per- and superoxide ions (alkali-metal sesquioxides)
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    if (cce_vars.num_perox_bonds > 0){
      string corrections_line="";
      corrections_line=get_corrections_line_O("O2_-2_O");
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " corrections line peroxides: " << corrections_line << endl;
      }
      vector<string> corrections_tokens;
      aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          for (uint l = 0; l < num_temps; l++) {
            cce_vars.perox_correction[num_temps*k+l]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+l+2]); // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " peroxide correction for " << cce_vars.vtemperatures[l] << "K: " << cce_vars.perox_correction[num_temps*k+l] << " eV/bond" << endl;
            }
          }
        } else {
          cce_vars.perox_correction[num_temps*k]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " peroxide correction for 298.15K: " << cce_vars.perox_correction[num_temps*k] << " eV/bond" << endl;
          }
        }
      }
    } 
    if (cce_vars.num_superox_bonds > 0) {
      string corrections_line="";
      corrections_line=get_corrections_line_O("O2_-1_O");
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " corrections line superoxides: " << corrections_line << endl;
      }
      vector<string> corrections_tokens;
      aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          for (uint l = 0; l < num_temps; l++) {
            cce_vars.superox_correction[num_temps*k+l]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+l+2]);
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " superoxide correction for " << cce_vars.vtemperatures[l] << "K: " << cce_vars.superox_correction[num_temps*k+l] << " eV/bond" << endl;
            }
          }
        } else {
          cce_vars.superox_correction[num_temps*k]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " superoxide correction for 298.15K: " << cce_vars.superox_correction[num_temps*k] << " eV/bond" << endl;
          }
        }
      }
    }
    if(LDEBUG){
      cerr << endl;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //         APPLY CORRECTIONS AND GET CORRECTED FORMATION ENTHALPIES        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //check_apply_per_super_ox_corrections////////////////////////////////////////////////////////
  // check whether corrections for per- or superoxide ions are needed and apply accordingly
  void check_apply_per_super_ox_corrections(CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    if (cce_vars.num_perox_bonds > 0){
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          for (uint l = 0; l < num_temps; l++) {
            cce_vars.cce_correction[num_temps*k+l] += (cce_vars.num_perox_bonds * cce_vars.perox_correction[num_temps*k+l]) ; // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " peroxide correction for " << cce_vars.vtemperatures[l] << "K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_perox_bonds * cce_vars.perox_correction[num_temps*k+l]) << " eV" << endl;
              cerr << endl;
            }
          }
        } else {
          cce_vars.cce_correction[num_temps*k] += (cce_vars.num_perox_bonds * cce_vars.perox_correction[num_temps*k]) ;
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " peroxide correction for 298.15K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_perox_bonds * cce_vars.perox_correction[num_temps*k]) << " eV" << endl;
          }
        }
      }
    }
    if (cce_vars.num_superox_bonds > 0){
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          for (uint l = 0; l < num_temps; l++) {
            cce_vars.cce_correction[num_temps*k+l] += (cce_vars.num_superox_bonds * cce_vars.superox_correction[num_temps*k+l]) ; // num_temps*k since to each functional belong num_temps corrections (currently for 298.15 and 0K)
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " superoxide correction for " << cce_vars.vtemperatures[l] << "K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_superox_bonds * cce_vars.superox_correction[num_temps*k+l]) << " eV" << endl;
              cerr << endl;
            }
          }
        } else {
          cce_vars.cce_correction[num_temps*k] += (cce_vars.num_superox_bonds * cce_vars.superox_correction[num_temps*k]) ;
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " superoxide correction for 298.15K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_superox_bonds * cce_vars.superox_correction[num_temps*k]) << " eV" << endl;
          }
        }
      }
    }
  }

  //apply_pbe_u_icsd_shifts////////////////////////////////////////////////////////
  // apply the shifts for the ref. enthalpies for PBE+U:ICSD if needed
  void apply_pbe_u_icsd_shifts(const xstructure& structure, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    stringstream message;
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] == "PBE+U:ICSD") { // for PBE+U:ICSD additional ref. enthalpy shifts must be applied since corrections have been fitted to ref. enthalpies calculated with U but internal AFLOW ref. enthalpies are from plain PBE
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << cce_vars.vfunctionals[k] << " applying ref. enthalpy shifts for PBE+U:ICSD." << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){
          for (uint l = 0; l < num_temps; l++) {
            double ref_enthalpy_shift=get_ref_enthalpy_shift_pbe_u_icsd(structure.atoms[i].cleanname);
            if (ref_enthalpy_shift == AUROSTD_NAN) { // for some species needing shifts, there is no reference (ground state) energy yet
              oss << print_output_oxidation_numbers(structure, cce_vars);
              message << " No ref. enthalpy shift for " << structure.atoms[i].cleanname << " for PBE+U:ICSD yet.";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
            } else if (ref_enthalpy_shift > 0) { // consider only species that need a ref. enthalpy shift, i.e. for which a U is used
              cce_vars.cce_correction[num_temps*k+l] += ref_enthalpy_shift;
            }
          }
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                        PRINT OUTPUT AND CITATION                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //print_output/////////////////////////////////////////////////////////////
  // Decides which output to print when running corrections
  string print_output(const xstructure& structure, aurostd::xoption& flags, xoption& cce_flags, CCE_Variables& cce_vars, vector<vector<uint> >& multi_anion_num_neighbors) {
    stringstream output;
    if ((flags.flag("CCE_CORRECTION::POSCAR_PATH") || flags.flag("CCE_CORRECTION::POSCAR2CCE")) && !flags.flag("CCE_CORRECTION::UNIT_TEST")) {
      if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") {
        output << "{";
        output << "\"cation_coordination_numbers\":";
        output << print_JSON_cation_coordination_numbers(structure, cce_flags, cce_vars, multi_anion_num_neighbors);
        output << ",\"oxidation_states\":";
        output << print_JSON_ox_nums(structure, cce_vars) << ",";
      } else {
        // print cation coordination numbers
        output << "CATION COORDINATION NUMBERS:" << std::endl;
        output << print_output_cation_coordination_numbers(structure, cce_flags, cce_vars, multi_anion_num_neighbors);
        // print oxidation numbers
        output << "OXIDATION NUMBERS:" << std::endl;
        output << print_output_oxidation_numbers(structure, cce_vars);
      }
    }
    if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") {
      if (!((flags.flag("CCE_CORRECTION::POSCAR_PATH") || flags.flag("CCE_CORRECTION::POSCAR2CCE")) && !flags.flag("CCE_CORRECTION::UNIT_TEST"))) {
        output << "{";
      }
      output << print_JSON_corrections(structure, cce_vars);
      output << "}" << std::endl;
    } else if (flags.flag("CCE_CORRECTION::UNIT_TEST")) {
      output << print_test_output(cce_vars, cce_vars.enthalpy_formation_cell_cce) << std::endl;
    } else {
      // print CCE corrections & corrected formation enthalpies per cell and atom
      output << print_output_corrections(structure, cce_vars, cce_vars.enthalpy_formation_cell_cce);
      // print CCE citation
      output << print_citation();
    }
    return output.str();
  }

  //print_JSON_cation_coordination_numbers/////////////////////////////////////////////////////////////
  // Returns cation coordination numbers, i.e. number of anions coordinating each cation in JSON format
  string print_JSON_cation_coordination_numbers(const xstructure& structure, xoption& cce_flags, const CCE_Variables& cce_vars, vector<vector<uint> >& multi_anion_num_neighbors) {
    stringstream json;
    vector<string> cations_neighbors_vector;
    vector<string> atoms_neighbors_vector;
    vector<string> anion_neighbors_vector;
    uint nspecies = structure.species.size();

    json << "{";
    uint l=0;
    for (uint i=0;i<nspecies;i++){
      atoms_neighbors_vector.clear();
      string cation_info="";
      for (int k = 0; k < structure.num_each_type[i]; k++) {
        if ((structure.atoms[l].cleanname != cce_vars.anion_species) && (cce_vars.multi_anion_atoms[l] != 1)){ // exclude main anion species and multi anion atoms detected previously
          anion_neighbors_vector.clear();
          cation_info="";
          string atom_info="";
          if (cce_vars.num_neighbors[l] > 0){ // are there actually bonds between the cation and the (main) anion species
            cation_info="\"" + structure.species[i] + "\":";
            anion_neighbors_vector.push_back("\"" + aurostd::utype2string<uint>(l) + "\":{" + "\"" + cce_vars.anion_species + "\":" + aurostd::utype2string<uint>(cce_vars.num_neighbors[l]));
          }
          if (cce_flags.flag("MULTI_ANION_SYSTEM")){
            for(uint j=0,jsize=cce_vars.multi_anion_species.size();j<jsize;j++){ 
              if (multi_anion_num_neighbors[j][l] > 0){ // are there actually bonds between the cation and the multi anion species
                if (cation_info == ""){ // if no neighbors from main anion species, append basic cation info first
                  cation_info="\"" + structure.species[i] + "\":";
                }
                if (cce_vars.num_neighbors[l] == 0){
                  atom_info="\"" + aurostd::utype2string<uint>(l) + "\":{";
                }
                anion_neighbors_vector.push_back("\"" + cce_vars.multi_anion_species[j] + "\":" + aurostd::utype2string<uint>(multi_anion_num_neighbors[j][l]));
              }
            }
          }
          atoms_neighbors_vector.push_back(atom_info + aurostd::joinWDelimiter(anion_neighbors_vector,",") + "}");
        }
        l++;
      }
      if (atoms_neighbors_vector.size() != 0) {
        cations_neighbors_vector.push_back(cation_info + "{" + aurostd::joinWDelimiter(atoms_neighbors_vector,",") + "}");
      }
    }
    json << aurostd::joinWDelimiter(cations_neighbors_vector,",");
    json << "}";
    return json.str();
  }

  //print_JSON_ox_nums/////////////////////////////////////////////////////////////
  // Returns oxidation numbers in JSON format
  string print_JSON_ox_nums(const xstructure& structure, const CCE_Variables& cce_vars) {
    stringstream json;
    uint nspecies = structure.species.size();

    json << "{";
    uint l=0;
    for (uint i = 0; i < nspecies; i++) {
      json << "\"" << structure.species[i] << "\":{";
      for (int k = 0; k < structure.num_each_type[i]; k++) {
        json << "\"" << aurostd::utype2string<uint>(l) << "\":" << cce_vars.oxidation_states[l];
        if (k < structure.num_each_type[i] - 1) json << ",";
        l++;
      }
      json << "}";
      if (i < nspecies - 1) json << ",";
    }
    json << "}";
    return json.str();
  }

  //print_JSON_corrections/////////////////////////////////////////////////////////////
  //ME20200213
  // Returns CCE results in JSON format
  string print_JSON_corrections(const xstructure& structure, const CCE_Variables& cce_vars) {
    stringstream json;
    bool print_Hf = (cce_vars.enthalpies_dft.size() > 0);
    uint nfuncs = cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    uint natoms = structure.atoms.size();

    //json << "{";
    //json << "\"oxidation_states\":";
    //json << "[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(cce_vars.oxidation_states),",") << "],";
    json << "\"CCE\":{";
    for (uint i = 0; i < nfuncs; i++) {
      json << "\"" << cce_vars.vfunctionals[i] << "\":{";
      if (cce_vars.vfunctionals[i] != "exp") {
        for (uint l = 0; l < num_temps; l++) {
          json << (l!=0?",":"") << "\"" << cce_vars.vtemperatures[l] << "K\":{";
          json << "\"cce_correction_cell\":" << cce_vars.cce_correction[num_temps*i+l] << ",";
          json << "\"cce_correction_atom\":" << (cce_vars.cce_correction[num_temps*i+l]/natoms);
          if (print_Hf) {
            json << ",";
            json << "\"formation_enthalpy_cell\":" << cce_vars.enthalpy_formation_cell_cce[num_temps*i+l] << ",";
            json << "\"formation_enthalpy_atom\":" << (cce_vars.enthalpy_formation_cell_cce[num_temps*i+l]/natoms);
          }
          json << "}";
        }
      } else {
        json << "\"298.15K\":{";
        json << "\"formation_enthalpy_cell\":" << cce_vars.enthalpy_formation_cell_cce[num_temps*i] << ",";
        json << "\"formation_enthalpy_atom\":" << (cce_vars.enthalpy_formation_cell_cce[num_temps*i]/natoms);
        json << "}";
      }
      json << "}";
      if (i < nfuncs - 1) json << ",";
    }
    json << "},";
    json << "\"publication\":";
    json << "\"Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019). ";
    json << "https://doi.org/10.1038/s41524-019-0192-1\"";
    //json << "}";
    return json.str();
  }

  //print_output_cation_coordination_numbers////////////////////////////////////////////////////////
  // print cation coordination numbers, i.e. number of anions coordinating each cation
  string print_output_cation_coordination_numbers(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, vector<vector<uint> >& multi_anion_num_neighbors) {
    stringstream output;
    output << std::setw(4) << std::right << "atom" << std::setw(13) << std::left << "   species" << std::setw(8) << "anion" << std::setw(13) << std::right << "coord. number" << endl;
    for (uint i=0,isize=structure.atoms.size();i<isize;i++){
      if ((structure.atoms[i].cleanname != cce_vars.anion_species) && (cce_vars.multi_anion_atoms[i] != 1)){ // exclude main anion species and multi anion atoms detected previously
        if (cce_vars.num_neighbors[i] > 0){ // are there actually bonds between the cation and the (main) anion species
          output << std::setw(4) << std::right << i+1 << std::setw(13) << std::left << "   " + structure.atoms[i].cleanname << std::setw(8) << cce_vars.anion_species << std::setw(13) << std::right << cce_vars.num_neighbors[i] << endl; // i+1: convert to 1-based counting
        }
        if (cce_flags.flag("MULTI_ANION_SYSTEM")){
          for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
            if (multi_anion_num_neighbors[k][i] > 0){ // are there actually bonds between the cation and the multi anion species
              output << std::setw(4) << std::right << i+1 << std::setw(13) << std::left << "   " + structure.atoms[i].cleanname << std::setw(8) << cce_vars.multi_anion_species[k] << std::setw(13) << std::right << multi_anion_num_neighbors[k][i] << endl; // i+1: convert to 1-based counting
            }
          }
        }
      }
    }
    output << endl;
    return output.str();
  }

  //print_output_oxidation_numbers////////////////////////////////////////////////////////
  // print oxidation numbers
  string print_output_oxidation_numbers(const xstructure& structure, CCE_Variables& cce_vars) {
    stringstream output;
    // print oxidation numbers
    output << std::setw(4) << std::right << "atom" << std::setw(13) << std::left << "   species" << std::setw(15) << std::right << "oxidation state" << endl;
    for (uint k=0,ksize=cce_vars.oxidation_states.size();k<ksize;k++){
      output << std::showpos << std::setw(4) << std::right << k+1 << std::setw(13) << std::left << "   " + structure.atoms[k].cleanname << std::setw(15) << std::right << cce_vars.oxidation_states[k] << endl; // k+1: convert to 1-based counting
    }
    output << endl;
    return output.str();
  }

  //print_output_corrections////////////////////////////////////////////////////////
  // print CCE corrections and corrected formation enthalpies if precalculated DFT values are provided
  string print_output_corrections(const xstructure& structure, CCE_Variables& cce_vars, const vector<double>& enthalpy_formation_cell_cce) {
    stringstream output;
    // print out CCE corrections per cell and atom for functionals selected
    if (!(cce_vars.vfunctionals.size() == 1 && cce_vars.vfunctionals[0] == "exp")){ // if only exp is set as functional CCE CORRECTIONS: should not be written
      output << "CCE CORRECTIONS (to be subtracted from precalculated DFT formation enthalpies):" << endl;
      output << std::setw(10) << "functional" << std::setw(14) << "temperature" << std::setw(13) << "correction" << std::setw(13) << "correction" << endl;
      output << std::setw(24) << "(K)" << std::setw(13) << "(eV/cell)" << std::setw(13) << "(eV/atom)" << endl;
    }
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        for (uint l = 0; l < num_temps; l++) {
          output << std::showpos << std::setprecision(3) << std::fixed << std::setw(10) << std::left << cce_vars.vfunctionals[k] << std::setw(14) << std::right << cce_vars.vtemperatures[l] << std::setw(13) << cce_vars.cce_correction[num_temps*k+l] << std::setw(13) << cce_vars.cce_correction[num_temps*k+l]/structure.atoms.size() << endl;
        }
      }
    }
    output << endl;
    // exp result should always be written at the end, hence print only after writing output for other functionals
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] == "exp" && cce_vars.enthalpies_dft.size()==0) { // second condition for that if precalc. form. enthalpies are given and asking for exp., exp. result is not written twice
        output << "CCE FORMATION ENTHALPIES:" << endl;
        output << std::setw(10) << "functional" << std::setw(14) << "temperature" << std::setw(17) << "form. enthalpy" << std::setw(17) << "form. enthalpy" << endl;
        output << std::setw(24) << "(K)" << std::setw(17) << "(eV/cell)" << std::setw(17) << "(eV/atom)" << endl;
        output << std::showpos << std::setprecision(3) << std::fixed << std::setw(10) << std::left << "CCE@exp" << std::setw(14) << std::right << "298.15" << std::setw(17) << enthalpy_formation_cell_cce[num_temps*k] << std::setw(17) << enthalpy_formation_cell_cce[num_temps*k]/structure.atoms.size() << endl;
        output << "Note that CCE@exp provides A ROUGH GUESS with an estimated average accuracy of only about 250 meV/atom (from test for ternary oxides)!" << endl;
      }
    }
    // print CCE formation enthalpies per cell and atom for functionals selected 
    // if precalculated DFT values are provided
    if(cce_vars.enthalpies_dft.size()!=0){ 
      output << "CCE FORMATION ENTHALPIES:" << endl;
      output << std::setw(10) <<  "functional" << std::setw(14) << "temperature" << std::setw(17) << "form. enthalpy" << std::setw(17) << "form. enthalpy" << endl;
      output << std::setw(24) << "(K)" << std::setw(17) << "(eV/cell)" << std::setw(17) << "(eV/atom)" << endl;
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          for (uint l = 0; l < num_temps; l++) {
            output << std::showpos << std::setprecision(3) << std::fixed << std::setw(10) << std::left << cce_vars.vfunctionals[k] << std::setw(14) << std::right << cce_vars.vtemperatures[l] << std::setw(17) << enthalpy_formation_cell_cce[num_temps*k+l] << std::setw(17) << enthalpy_formation_cell_cce[num_temps*k+l]/structure.atoms.size() << endl;
          }
        } else if (cce_vars.vfunctionals[k] == "exp") {
          output << endl;
          output << std::showpos << std::setprecision(3) << std::fixed << std::setw(10) << std::left << "CCE@exp" << std::setw(14) << std::right << "298.15" << std::setw(17) << enthalpy_formation_cell_cce[num_temps*k] << std::setw(17) << enthalpy_formation_cell_cce[num_temps*k]/structure.atoms.size() << endl;
          output << "Note that CCE@exp provides A ROUGH GUESS with an estimated average accuracy of only about 250 meV/atom (from test for ternary oxides)!" << endl;
        }
      }
    }
    output << endl;
    return output.str();
  }

  //print_test_output////////////////////////////////////////////////////////
  // print CCE corrections and corrected formation enthalpies for testing
  string print_test_output(CCE_Variables& cce_vars, const vector<double>& enthalpy_formation_cell_cce) {
    stringstream output;
    uint num_funcs=cce_vars.vfunctionals.size();
    uint num_temps=cce_vars.vtemperatures.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        for (uint l = 0; l < num_temps; l++) {
          output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[num_temps*k+l] << (k<num_funcs-1 || l==0 || cce_vars.enthalpies_dft.size()!=0?",":"");
        }
      }
    }
    // exp part
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] == "exp" && cce_vars.enthalpies_dft.size()==0) { // second condition for that if precalc. form. enthalpies are given and asking for exp., exp. result is not written twice
        output << std::setprecision(3) << std::fixed << enthalpy_formation_cell_cce[num_temps*k];
      }
    }
    // corrected DFT formation enthalpies
    if(cce_vars.enthalpies_dft.size()!=0){
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          for (uint l = 0; l < num_temps; l++) {
            output << std::setprecision(3) << std::fixed << enthalpy_formation_cell_cce[num_temps*k+l] << (k<num_funcs-1 || l==0?",":"");
          }
        } else if (cce_vars.vfunctionals[k] == "exp") {
          output << std::setprecision(3) << std::fixed << enthalpy_formation_cell_cce[num_temps*k] << (k<num_funcs-1?",":"");
        }
      }
    }
    return output.str();
  }

  //print_citation////////////////////////////////////////////////////////
  // print citation information at end of output
  string print_citation() {
    stringstream oss;
    oss << "#############################################################################################" << endl;
    oss << "When you use results from CCE and/or this implementation, please cite the following articles:" << endl;
    oss << "https://doi.org/10.1038/s41524-019-0192-1; https://doi.org/10.1103/PhysRevMaterials.5.043803" << endl;
    oss << "#############################################################################################" << endl;
    return oss.str();
  }

  //print_usage////////////////////////////////////////////////////////
  // For printing user instructions if no additional input is provided.
  string print_usage() {
    stringstream oss;
    oss << endl;
    oss << "Written by Rico Friedrich, Corey Oses, and Marco Esters, 2018-2021" << endl;
    oss << endl;
    oss << "USER INSTRUCTIONS:" << endl;
    oss << endl;
    oss << "(i) GENERAL INFORMATION:" << endl;
    oss << "Implementation to obtain corrected DFT formation enthalpies based on the coordination corrected" << endl; 
    oss << "enthalpies (CCE) methodology described in:" << endl; 
    oss << "Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019);" << endl;
    oss << "https://doi.org/10.1038/s41524-019-0192-1" << endl;
    oss << "and" << endl;
    oss << "Friedrich et al., Automated coordination corrected enthalpies with AFLOW-CCE, Phys. Rev. Mater. 5, 043803 (2021);" << endl;
    oss << "https://doi.org/10.1103/PhysRevMaterials.5.043803" << endl;
    oss << "Please cite these articles when using this method and/or this implementation." << endl;
    oss << "The corrections depend on the number of cation-anion bonds and on the cation oxidation state." << endl;
    oss << "More general information and a list of elements and oxidation states for which corrections are available" << endl; 
    oss << "can be found in README_AFLOW_CCE.TXT." << endl;
    oss << endl;
    oss << "(ii) AVAILABLE OPTIONS:" << endl;
    oss << "--cce                            Prints these user instructions." << endl;
    oss << endl;
    oss << "--cce=STRUCTURE_FILE_PATH        Provide the path to the structure file. It can be in any structure" << endl;
    oss << "                                 format that AFLOW supports, e.g. VASP POSCAR, QE, AIMS, ABINIT, ELK, and CIF." << endl;
    oss << "                                 For VASP, a VASP5 POSCAR is required or, if a VASP4 POSCAR is used, the species" << endl;
    oss << "                                 must be written on the right side next to the coordinates for each atom" << endl;
    oss << "                                 just as for the EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO below." << endl;
    oss << endl;
    oss << "--enthalpies_formation_dft=|--dfte=" << endl;
    oss << "                                 Provide a comma separated list of precalculated DFT formation enthalpies," << endl; 
    oss << "                                 they are assumed to be: (i) negative for compounds lower in enthalpy" << endl; 
    oss << "                                 than the elements, (ii) in eV/cell. Currently, corrections are available" << endl; 
    oss << "                                 for PBE, LDA and SCAN." << endl;
    oss << endl;
    oss << "--functionals=|--func=|--functional=" << endl;
    oss << "                                 Provide a comma separated list of functionals for which corrections" << endl;
    oss << "                                 should be returned. If used together with --dft_formation enthalpies," << endl;
    oss << "                                 the functionals must be in the same sequence as the DFT formation" << endl;
    oss << "                                 enthalpies they correspond to. Available functionals are:" << endl;
    oss << "                                 (i) PBE, (ii) LDA or (iii) SCAN. Default: PBE (if only one DFT formation" << endl; 
    oss << "                                 enthalpy is provided)." << endl;
    oss << endl;
    oss << "--oxidation_numbers=|--ox_nums=|--oxidation_number=" << endl;
    oss << "                                 Provide as a comma separated list the oxidation numbers. It is" << endl;
    oss << "                                 assumed that: (i) one is provided for each atom of the structure and" << endl; 
    oss << "                                 (ii) they are in the same sequence as the corresponding atoms in the" << endl;
    oss << "                                 provided structure file." << endl;
    oss << endl;
    oss << "--get_cce_correction|--get_cce_cor|--poscar2cce < STRUCTURE_FILE_PATH" << endl;
    oss << "                                 Determines the CCE corrections for the structure in STRUCTURE_FILE_PATH." << endl;
    oss << "                                 It can be in any structure format that AFLOW supports, e.g. VASP POSCAR," << endl;
    oss << "                                 QE, AIMS, ABINIT, ELK, and CIF. For VASP, a VASP5 POSCAR is required or, if a" << endl;
    oss << "                                 VASP4 POSCAR is used, the species must be written on the right side next to" << endl;
    oss << "                                 the coordinates for each atom just as for the EXAMPLE INPUT STRUCTURE FOR" << endl;
    oss << "                                 ROCKSALT MgO below." << endl;
    oss << endl;
    oss << "--get_oxidation_numbers|--get_ox_nums|--poscar2ox_nums < STRUCTURE_FILE_PATH" << endl;
    oss << "                                 Determines the oxidation numbers for the structure in STRUCTURE_FILE_PATH." << endl;
    oss << "                                 It can be in any structure format that AFLOW supports, e.g. VASP POSCAR," << endl;
    oss << "                                 QE, AIMS, ABINIT, ELK, and CIF. For VASP, a VASP5 POSCAR is required or if a" << endl;
    oss << "                                 VASP4 POSCAR is used, the species must be written on the right side next to" << endl;
    oss << "                                 the coordinates for each atom just as for the EXAMPLE INPUT STRUCTURE FOR" << endl;
    oss << "                                 ROCKSALT MgO below." << endl;
    oss << endl;
    oss << "--get_cation_coordination_numbers|--get_cation_coord_nums|--poscar2cation_coord_nums < STRUCTURE_FILE_PATH" << endl;
    oss << "                                 Determines the number of anion neighbors for each cation for the structure" << endl;
    oss << "                                 in STRUCTURE_FILE_PATH. It can be in any structure format that AFLOW" << endl;
    oss << "                                 supports, e.g. VASP POSCAR, QE, AIMS, ABINIT, ELK, and CIF. For VASP, a" << endl;
    oss << "                                 VASP5 POSCAR is required or if a VASP4 POSCAR is used, the species must" << endl;
    oss << "                                 be written on the right side next to the coordinates for each atom" << endl;
    oss << "                                 just as for the EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO below." << endl;
    oss << endl;
    oss << "--print=                         Obtain output in standard format (--print=out) or as json (--print=json)." << endl;
    oss << "                                 Default: out." << endl;
    oss << endl;
    oss << "(iii) EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO:" << endl;
    oss << "Mg1O1   [FCC,FCC,cF8] (STD_PRIM doi:10.1  [FCC,FCC,cF8] (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)" << endl;
    oss << "1.224745" << endl;
    oss << "   0.00000000000000   1.73568248770103   1.73568248770103" << endl;
    oss << "   1.73568248770103   0.00000000000000   1.73568248770103" << endl;
    oss << "   1.73568248770103   1.73568248770103   0.00000000000000" << endl;
    oss << "Mg O" << endl;
    oss << "1 1" << endl;
    oss << "Direct(2) [A1B1]" << endl;
    oss << "   0.00000000000000   0.00000000000000   0.00000000000000  Mg" << endl;
    oss << "   0.50000000000000   0.50000000000000   0.50000000000000  O" << endl;
    oss << endl;
    oss << "(iv) EXAMPLE COMMANDS:" << endl;
    oss << "Assuming that AFLOW is in your PATH and you saved the above example structure file for MgO" << endl; 
    oss << "in the current directory as POSCAR, the following commands can be executed:" << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --enthalpies_formation_dft=-5.434,-6.220,-6.249 --functionals=PBE,LDA,SCAN" << endl;
    oss << "This will give you the CCE corrections and CCE formation enthalpies for PBE, LDA, and SCAN for MgO." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --enthalpies_formation_dft=-6.220 --functionals=LDA" << endl;
    oss << "This gives you only the CCE corrections and CCE formation enthalpies for LDA." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --enthalpies_formation_dft=-5.434" << endl;
    oss << "This gives you the CCE corrections and CCE formation enthalpies for PBE with a warning that" << endl; 
    oss << "PBE is assumed as functional." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR" << endl;
    oss << "This gives you the CCE corrections for PBE, LDA, and SCAN and a rough guess of the formation" << endl;
    oss << "enthalpy based on experimental formation enthalpies per bond." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --oxidation_numbers=2,-2" << endl;
    oss << "Oxidation numbers for each atom can also be provided as input." << endl;
    oss << endl;
    return oss.str();
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //       LOOK UP TABLES FOR THE ACTUAL CORRECTIONS AND BADER CHARGES       //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //get_corrections_line////////////////////////////////////////////////////////
  // function to get corrections
  //CCE corrections per bond for DFT formation enthalpies for polar materials from binary data using AFLOW (PBE, LDA, SCAN, and PBE+U:ICSD) 
  //with PAW data sets for VASP 5.4.4 and measured experimental values according to the CCE paper 
  //Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019). 
  //The extensions after the binary formula indicate the source from which the experimental value was taken: 
  //NJ=NIST_JANAF 1998; Ba=Barin 1995; if no extension is given, the value from Kubaschewski et al. 1993 was used. 
  //For DFT+U corrections will only be included for the PBE(+U) calculations used to generate the AFLOW ICSD data
  //DFT+U calculations differ not only by the U value but also by other specifications of the implementation, see:
  //Kick, Reuter, and Oberhofer; Intricacies of DFT+U, Not Only in a Numeric Atom Centered Orbital Framework, JCTC (2019); DOI: 10.1021/acs.jctc.8b01211
  //The formation enthalpies per bond from the experimental data (Exp.) need to be multiplied by the number of M-anion bonds 
  //and added (sign change compared to usual CCE corrections to be substracted from DFT formation enthalpies) 
  //for all cations to get an estimate for the formation enthalpy only based on experimental data; 
  //positive values for the peroxide and superoxide O-O bond formation enthalpy indicate that taking the cation-O2- 
  //correction for the cation-O bonds is not a very good approximation for per- and superoxides, 
  //i.e. the oxidation state dependence with respect to the oxidation state of O is relevant!)
  string get_corrections_line(const string& considered_anion_species, const string& cor_identifier) {
    if (considered_anion_species=="O") {return get_corrections_line_O(cor_identifier);}
    if (considered_anion_species=="N") {return get_corrections_line_N(cor_identifier);}
    else {return "";}
  }

  //get_corrections_line_O////////////////////////////////////////////////////////
  string get_corrections_line_O(const string& cor_identifier) {
    //The number of M-anion bonds per formula unit in the last column are for the compounds used to get the corrections (at the start of the return string).
    //They can be used to recalculate the used experimental values by multiplying with the exp. corrections per bond in the previous column.
    //For per- and superoxides the number given is the number of O-O bonds per formula unit of the compound used to get the correction.
    //                                            cation species          ox. state  PBE (for 298.15K)  PBE (for 0K)  LDA (for 298.15K)  LDA (for 0K)  SCAN (for 298.15K)  SCAN (for 0K)  PBE(+U)_ICSD (for 298.15K)  PBE(+U)_ICSD (for 0K)  Exp. (for 298.15K)  M-X bonds        Protos                              n cations   n anions        H_f^exp  (for 298.15K)
    //                                                                               (eV/bond)          (eV/bond)     (eV/bond)          (eV/bond)     (eV/bond)           (eV/bond)      (eV/bond)                   (eV/bond)              (eV/bond)           per f. u.                                            per f. u.   per f. u.       (eV/bond) TABLE IV in supp.
    //OXIDES                                                                                                                                                                                                                                                                                                                                 
    if (cor_identifier=="Ag_+1_O")       {return "Ag_+1_Ag2O             +1         -0.00825           -0.00700      -0.05400           -0.05375      -0.06525            -0.06450       -0.06150                    -0.06000               -0.08050            4                 AgO/A2B_cP6_224_b_a.AB              2           1                -0.107   ";}
    else if (cor_identifier=="Al_+3_O")  {return "Al_+3_Al2O3            +3          0.18692            0.17783      -0.00733           -0.01683      -0.01242            -0.02217        0.18667                     0.17758               -1.44725            12                AlO/A2B3_hR10_167_c_e-001.AB        2           3                -3.473   ";}
    else if (cor_identifier=="As_+5_O")  {return "As_+5_As2O5            +5          0.20220            0.19190      -0.06360           -0.07520       0.02120             0.00920        0.19770                     0.18760               -0.95360            10                AsO/A2B5_oP28_19_2a_5a-001.AB       2           5                -1.362   ";}
    else if (cor_identifier=="B_+3_O")   {return "B_+3_B2O3              +3          0.19517            0.18250      -0.06933           -0.08350      -0.04717            -0.06117        0.19317                     0.18067               -2.19983            6                 B_hO/A2B3_hP15_144_2a_3a-001.AB     2           3                -2.64    ";}
    else if (cor_identifier=="Ba_+2_O")  {return "Ba_+2_BaO              +2          0.11833            0.11650       0.00983            0.00750       0.00417             0.00200        0.11817                     0.11633               -0.94683            6                 Ba_svO/AB_cF8_225_a_b.AB            1           1                -2.841   ";}
    else if (cor_identifier=="Be_+2_O")  {return "Be_+2_BeO              +2          0.19525            0.18750       0.00825            0.00000       0.00600            -0.00225        0.19475                     0.18700               -1.57900            4                 Be_svO/AB_hP4_186_b_b-001.AB        1           1                -3.158   ";}
    else if (cor_identifier=="Bi_+3_O")  {return "Bi_+3_Bi2O3            +3         -0.02580           -0.02760      -0.17520           -0.17780      -0.03560            -0.03790       -0.02830                    -0.03010               -0.59150            10                Bi_dO/A2B3_mP20_14_2e_3e-001.AB     2           3                -1.183   ";}
    else if (cor_identifier=="Ca_+2_O")  {return "Ca_+2_CaO              +2          0.10567            0.10017      -0.03317           -0.03950      -0.02217            -0.02800        0.10567                     0.10017               -1.09667            6                 Ca_svO/AB_cF8_225_a_b.AB            1           1                -3.29    ";}
    else if (cor_identifier=="Cd_+2_O")  {return "Cd_+2_CdO              +2          0.10417            0.10133       0.01617            0.01300       0.00533             0.00233        0.07417                     0.07117               -0.44633            6                 CdO/AB_cF8_225_a_b.AB               1           1                -1.339   ";}
    else if (cor_identifier=="Co_+2_O")  {return "Co_+2_CoO              +2          0.23967            0.23733       0.16617            0.16450       0.12300             0.11933        0.00300                    -0.00017               -0.41067            6                 CoO/AB_cF8_225_a_b.AB               1           1                -1.232   ";}
    else if (cor_identifier=="Cr_+3_O")  {return "Cr_+3_Cr2O3            +3          0.15283            0.14733       0.04542            0.03908      -0.01892            -0.02467       -0.18417                    -0.18958               -0.98000            12                Cr_pvO/A2B3_hR10_167_c_e-001.AB     2           3                -2.352   ";}
    else if (cor_identifier=="Cr_+6_O")  {return "Cr_+6_CrO3             +6         -0.13225           -0.14425      -0.30525           -0.32100      -0.28125            -0.29675       -0.06775                    -0.07950               -1.52100            4                 Cr_pvO/AB3_oC16_40_b_a2b-001.AB     1           3                -1.521   ";}
    else if (cor_identifier=="Cs_+1_O")  {return "Cs_+1_Cs2O             +1          0.10083            0.10083      -0.05667           -0.05833      -0.00500            -0.00600        0.10117                     0.10117               -0.59767            6                 Cs_svO/A2B_hR3_166_c_a-001.AB       2           1                -1.195   ";}
    else if (cor_identifier=="Cu_+1_O")  {return "Cu_+1_Cu2O_NJ          +1          0.13100            0.12950       0.03400            0.03175       0.06350             0.06175        0.06800                     0.06625               -0.44225            4                 Cu_pvO/A2B_cP6_224_b_a.AB           2           1                -0.59    ";}
    else if (cor_identifier=="Cu_+2_O")  {return "Cu_+2_CuO_NJ           +2          0.10150            0.09725      -0.02575           -0.03075       0.00675             0.00200        0.06525                     0.06075               -0.40425            4                 Cu_pvO/AB_mC8_15_a_e-001.AB         1           1                -0.808   ";}
    else if (cor_identifier=="Fe_+2_O")  {return "Fe_+2_FeO_NJ           +2          0.17633            0.17283       0.13100            0.12867       0.01883             0.01433       -0.00900                    -0.01300               -0.47000            6                 Fe_pvO/AB_cF8_225_a_b.AB            1           1                -1.41    ";}
    else if (cor_identifier=="Fe_+3_O")  {return "Fe_+3_Fe2O3            +3          0.16325            0.15858       0.01300            0.00550      -0.06550            -0.07175       -0.03417                    -0.04517               -0.71117            12                Fe_pvO/A2B3_hR10_167_c_e-001.AB     2           3                -1.707   ";}
    else if (cor_identifier=="Ga_+3_O")  {return "Ga_+3_Ga2O3            +3          0.20090            0.19250       0.00680           -0.00220       0.03860             0.02890        0.16220                     0.15370               -1.12880            10                Ga_hO/A2B3_mC20_12_2i_3i-001.AB     2           3                -2.258   ";}
    else if (cor_identifier=="Ge_+4_O")  {return "Ge_+4_GeO2             +4          0.19917            0.18950      -0.03567           -0.04617       0.03967             0.02900        0.19883                     0.18917               -1.00183            6                 Ge_hO/A2B_tP6_136_f_a-001.BA        1           2                -2.004   ";}
    else if (cor_identifier=="Hf_+4_O")  {return "Hf_+4_HfO2_Ba          +4          0.16171            0.15657      -0.02957           -0.03529      -0.02686            -0.03257        0.16286                     0.15786               -1.69486            7                 Hf_pvO/A2B_mP12_14_2e_e-001.BA      1           2                -3.955   ";}
    else if (cor_identifier=="Hg_+2_O")  {return "Hg_+2_HgO              +2          0.17000            0.15250      -0.06500           -0.08650       0.05800             0.03800        0.16050                     0.14300               -0.47050            2                 HgO/AB_oP8_62_c_c-006.AB            1           1                -0.47    ";}
    else if (cor_identifier=="In_+3_O")  {return "In_+3_In2O3            +3          0.13525            0.13025      -0.01667           -0.02217      -0.01300            -0.01858        0.09750                     0.09258               -0.79967            12                In_dO/A2B3_cI80_199_a2b_2c-001.AB   2           3                -1.919   ";}
    else if (cor_identifier=="Ir_+4_O")  {return "Ir_+4_IrO2_Ba          +4          0.02017            0.01567      -0.18133           -0.18683       0.01800             0.01350       -0.04283                    -0.04717               -0.41917            6                 IrO/A2B_tP6_136_f_a-001.BA          1           2                -0.838   ";}
    else if (cor_identifier=="K_+1_O")   {return "K_+1_K2O               +1          0.08300            0.08025      -0.02625           -0.03013      -0.00350            -0.00712        0.08262                     0.07987               -0.47050            8                 K_svO/AB2_cF12_225_a_c.BA           2           1                -1.255   ";}
    else if (cor_identifier=="Li_+1_O")  {return "Li_+1_Li2O             +1          0.07663            0.07038      -0.01538           -0.02225      -0.01175            -0.01863        0.08400                     0.07775               -0.77463            8                 Li_svO/AB2_cF12_225_a_c.BA          2           1                -2.066   ";}
    else if (cor_identifier=="Mg_+2_O")  {return "Mg_+2_MgO              +2          0.13350            0.12717       0.00250           -0.00417      -0.00233            -0.00900        0.13183                     0.12550               -1.03917            6                 Mg_pvO/AB_cF8_225_a_b.AB            1           1                -3.118   ";}
    else if (cor_identifier=="Mn_+2_O")  {return "Mn_+2_MnO              +2          0.25133            0.25133       0.27000            0.26933      -0.03250            -0.03400       -0.01650                    -0.02033               -0.66483            6                 Mn_pvO/AB_cF8_225_a_b.AB            1           1                -1.994   ";}
    else if (cor_identifier=="Mn_+4_O")  {return "Mn_+4_MnO2             +4          0.06000            0.05233      -0.09433           -0.10300      -0.15817            -0.16667        0.11017                     0.10233               -0.89983            6                 Mn_pvO/A2B_tP6_136_f_a-001.BA       1           2                -1.8     ";}
    else if (cor_identifier=="Mo_+4_O")  {return "Mo_+4_MoO2             +4          0.02917            0.02150      -0.18433           -0.19267      -0.10533            -0.11367        0.06400                     0.05650               -1.01550            6                 Mo_pvO/A2B_mP12_14_2e_e-002.BA      1           2                -2.031   ";}
    else if (cor_identifier=="Mo_+6_O")  {return "Mo_+6_MoO3             +6         -0.04700           -0.06025      -0.34100           -0.35750      -0.25575            -0.27175        0.07800                     0.06375               -1.93075            4                 Mo_pvO/AB3_oP16_62_c_3c-001.AB      1           3                -1.931   ";}
    else if (cor_identifier=="Na_+1_O")  {return "Na_+1_Na2O_NJ          +1          0.08225            0.07762      -0.00425           -0.00962      -0.01125            -0.01675        0.08200                     0.07737               -0.54150            8                 Na_svO/AB2_cF12_225_a_c.BA          2           1                -1.444   ";}
    else if (cor_identifier=="Nb_+2_O")  {return "Nb_+2_NbO              +2          0.05925            0.05325      -0.12625           -0.13275      -0.08450            -0.09100        0.03975                     0.03375               -1.08750            4                 Nb_svO/AB_cP6_221_c_d.AB            1           1                -2.175   ";}
    else if (cor_identifier=="Ni_+2_O")  {return "Ni_+2_NiO              +2          0.25583            0.25367       0.15517            0.15117       0.19800             0.19550       -0.07683                    -0.08017               -0.41400            6                 Ni_pvO/AB_cF8_225_a_b.AB            1           1                -1.242   ";}
    else if (cor_identifier=="Os_+4_O")  {return "Os_+4_OsO2             +4          0.06167            0.05700      -0.14433           -0.14983       0.00117            -0.00400       -0.01267                    -0.01717               -0.50883            6                 OOs_pv/A2B_tP6_136_f_a-001.AB       1           2                -1.018   ";}
    else if (cor_identifier=="Os_+8_O")  {return "Os_+8_OsO4             +8         -0.22250           -0.22950      -0.37925           -0.39200      -0.27725            -0.28800       -0.22450                    -0.23175               -1.02000            4                 OOs_pv/A4B_mC20_15_2f_e-001.AB      1           4                -0.816   ";}
    else if (cor_identifier=="Pb_+2_O")  {return "Pb_+2_PbO              +2          0.00275            0.00325      -0.10875           -0.10925      -0.05475            -0.05500       -0.00075                    -0.00025               -0.56850            4                 OPb_d/AB_tP4_129_a_c-001.AB         1           1                -1.137   ";}
    else if (cor_identifier=="Pb_+4_O")  {return "Pb_+4_PbO2             +4          0.05683            0.05450      -0.12500           -0.12850      -0.02500            -0.02833        0.04683                     0.04467               -0.47417            6                 OPb_d/A2B_tP6_136_f_a-001.AB        1           2                -0.948   ";}
    else if (cor_identifier=="Pd_+2_O")  {return "Pd_+2_PdO              +2          0.05775            0.05475      -0.07925           -0.08300      -0.02450            -0.02800       -0.02650                    -0.02950               -0.29925            4                 OPd_pv/AB_tP4_131_c_e-002.BA        1           1                -0.599   ";}
    else if (cor_identifier=="Rb_+1_O")  {return "Rb_+1_Rb2O             +1          0.09500            0.09400      -0.02150           -0.02350       0.00488             0.00313        0.09450                     0.09350               -0.43900            8                 ORb_sv/AB2_cF12_225_a_c.AB          2           1                -1.171   ";}
    else if (cor_identifier=="Re_+4_O")  {return "Re_+4_ReO2_Ba          +4          0.08967            0.08450      -0.12417           -0.13017       0.02117             0.01550        0.05067                     0.04650               -0.77550            6                 ORe_pv/A2B_mP12_14_2e_e-002.AB      1           2                -1.551   ";}
    else if (cor_identifier=="Re_+6_O")  {return "Re_+6_ReO3_Ba          +6         -0.07267           -0.08033      -0.30400           -0.31250      -0.15967            -0.16817       -0.08233                    -0.08967               -1.01767            6                 ORe_pv/A3B_cP4_221_d_a.AB           1           3                -1.526   ";}
    else if (cor_identifier=="Rh_+3_O")  {return "Rh_+3_Rh2O3            +3          0.01408            0.00650      -0.13633           -0.14150      -0.05833            -0.06308       -0.11975                    -0.12783               -0.30717            12                ORh_pv/A2B3_hR10_167_c_e-001.BA     2           3                -0.737   ";}
    else if (cor_identifier=="Ru_+4_O")  {return "Ru_+4_RuO2             +4         -0.00483           -0.01150      -0.20567           -0.21333      -0.11183            -0.11917       -0.05283                    -0.05800               -0.52683            6                 ORu_pv/A2B_tP6_136_f_a-001.AB       1           2                -1.054   ";}
    else if (cor_identifier=="Sb_+3_O")  {return "Sb_+3_Sb2O3            +3          0.12067            0.11533      -0.11350           -0.12117      -0.01983            -0.02667        0.11733                     0.11217               -1.23700            6                 OSb/A3B2_cF80_227_f_e-001.AB        2           3                -1.484   ";}
    else if (cor_identifier=="Sb_+5_O")  {return "Sb_+5_Sb2O5_Ba         +5          0.10517            0.09700      -0.13233           -0.14175      -0.05725            -0.06692        0.09825                     0.09008               -0.83942            12                OSb/A2B5_mC28_15_f_e2f-002.BA       2           5                -1.439   ";}
    else if (cor_identifier=="Sc_+3_O")  {return "Sc_+3_Sc2O3            +3          0.16175            0.15408      -0.00833           -0.01658      -0.02567            -0.03383        0.16608                     0.15867               -1.64817            12                OSc_sv/A2B3_cI80_206_ad_e-001.BA    2           3                -3.956   ";}
    else if (cor_identifier=="Se_+4_O")  {return "Se_+4_SeO2             +4          0.07500            0.06567      -0.22767           -0.23967      -0.09633            -0.10833        0.07300                     0.06367               -0.77767            3                 OSe/A2B_tP24_135_gh_h-001.AB        1           2                -0.778   ";}
    else if (cor_identifier=="Si_+4_O")  {return "Si_+4_SiO2(al-quartz)  +4          0.25300            0.23800      -0.02325           -0.03900      -0.02075            -0.03675        0.25275                     0.23775               -2.36025            4                 OSi/A2B_hP9_152_c_a-001.AB          1           2                -3.147   ";}
    else if (cor_identifier=="Sn_+2_O")  {return "Sn_+2_SnO_Ba           +2          0.06700            0.06500      -0.06375           -0.06650      -0.01300            -0.01575        0.12025                     0.11825               -0.74050            4                 OSn/AB_tP4_129_a_c-001.AB           1           1                -1.481   ";}
    else if (cor_identifier=="Sn_+4_O")  {return "Sn_+4_SnO2_Ba          +4          0.15050            0.14333      -0.04600           -0.05400      -0.01533            -0.02367        0.28650                     0.28000               -1.00333            6                 OSn/A2B_tP6_136_f_a-001.AB          1           2                -2.007   ";}
    else if (cor_identifier=="Sr_+2_O")  {return "Sr_+2_SrO              +2          0.10817            0.10483      -0.01917           -0.02317      -0.01550            -0.01933        0.10767                     0.10433               -1.02267            6                 OSr_sv/AB_cF8_225_a_b.AB            1           1                -3.068   ";}
    else if (cor_identifier=="Te_+4_O")  {return "Te_+4_TeO2_Ba          +4          0.06300            0.05575      -0.20325           -0.21225      -0.08850            -0.09700        0.06275                     0.05575               -0.83800            4                 OTe/A2B_tP12_92_b_a-002.AB          1           2                -1.117   ";}
    else if (cor_identifier=="Ti_+2_O")  {return "Ti_+2_TiO              +2          0.11313            0.10667      -0.06667           -0.07375      -0.02646            -0.03312        0.10479                     0.09958               -1.17188            4.8               OTi_sv/AB_mC20_12_b2i_c2i-001.AB    1           1                -2.812   ";}
    else if (cor_identifier=="Ti_+3_O")  {return "Ti_+3_Ti2O3            +3          0.09800            0.09050      -0.08550           -0.09358      -0.06875            -0.07667        0.07933                     0.07233               -1.31358            12                OTi_sv/A2B3_hR10_167_c_e-001.BA     2           3                -3.153   ";}
    else if (cor_identifier=="Ti_+4_O")  {return "Ti_+4_TiO2(rutile)     +4          0.10717            0.09717      -0.09650           -0.10717      -0.12367            -0.13450        0.05217                     0.04233               -1.63067            6                 OTi_sv/A2B_tP6_136_f_a-001.AB       1           2                -3.261   ";}
    else if (cor_identifier=="Tl_+1_O")  {return "Tl_+1_Tl2O             +1         -0.00650           -0.00533      -0.06650           -0.06600      -0.06117            -0.06050       -0.00783                    -0.00667               -0.28917            6                 OTl_d/AB2_hR6_166_c_2c-001.AB       2           1                -0.578   ";}
    else if (cor_identifier=="Tl_+3_O")  {return "Tl_+3_Tl2O3            +3          0.05358            0.05175      -0.09358           -0.09617      -0.01400            -0.01658        0.04492                     0.04300               -0.33717            12                OTl_d/A2B3_cI80_206_ad_e-001.BA     2           3                -0.809   ";}
    else if (cor_identifier=="V_+2_O")   {return "V_+2_VO                +2          0.26367            0.26200       0.12033            0.11517       0.15683             0.15467        0.24633                     0.24167               -0.74583            6                 OV_sv/AB_cF8_225_a_b.AB             1           1                -2.237   ";}
    else if (cor_identifier=="V_+3_O")   {return "V_+3_V2O3              +3          0.09825            0.09183      -0.06608           -0.07342      -0.05350            -0.06000       -0.04833                    -0.05458               -1.05267            12                OV_sv/A2B3_hR10_167_c_e-002.BA      2           3                -2.526   ";}
    else if (cor_identifier=="V_+4_O")   {return "V_+4_VO2               +4          0.04667            0.03750      -0.14967           -0.15983      -0.15367            -0.16367       -0.02633                    -0.03600               -1.23300            6                 OV_sv/A2B_mP12_14_2e_e-003.AB       1           2                -2.466   ";}
    else if (cor_identifier=="V_+5_O")   {return "V_+5_V2O5              +5         -0.00820           -0.01890      -0.21180           -0.22480      -0.21790            -0.23070        0.00730                    -0.00370               -1.60670            10                OV_sv/A5B2_oP14_59_a2e_e-001.AB     2           5                -2.295   ";}
    else if (cor_identifier=="W_+4_O")   {return "W_+4_WO2               +4          0.05667            0.05117      -0.15867           -0.16483      -0.04833            -0.05433        0.05417                     0.04883               -1.01833            6                 OW_pv/A2B_mP12_14_2e_e-002.AB       1           2                -2.037   ";}
    else if (cor_identifier=="W_+6_O")   {return "W_+6_WO3               +6          0.00500           -0.00250      -0.20783           -0.21650      -0.13617            -0.14483        0.00867                     0.00117               -1.45567            6                 OW_pv/A3B_mP16_14_3e_e-001.AB       1           3                -2.183   ";}
    else if (cor_identifier=="Y_+3_O")   {return "Y_+3_Y2O3              +3          0.13583            0.13042      -0.02192           -0.02800      -0.04825            -0.05425        0.15800                     0.15250               -1.64533            12                OY_sv/A2B3_cI80_206_ad_e-001.BA     2           3                -3.949   ";}
    else if (cor_identifier=="Zn_+2_O")  {return "Zn_+2_ZnO              +2          0.18525            0.18025       0.04225            0.03675       0.04550             0.03975        0.00300                    -0.00250               -0.90825            4                 OZn/AB_hP4_186_b_b-001.AB           1           1                -1.817   ";}
    else if (cor_identifier=="Zr_+4_O")  {return "Zr_+4_ZrO2_NJ          +4          0.13929            0.13200      -0.04514           -0.05300      -0.06314            -0.07086        0.18086                     0.17371               -1.62486            7                 OZr_sv/A2B_mP12_14_2e_e-001.AB      1           2                -3.791   ";}
    // corrections for per- and superoxides                                                                                                                                                                                                                                                                                                 
    else if (cor_identifier=="O2_-2_O")  {return "O2_-2_Li2O2            -1         -0.09300           -0.08556      -0.11600           -0.11100       0.24000             0.24756       -0.12700                    -0.12000                2.72556            1                 -                                   2           2                0.0      ";}
    else if (cor_identifier=="O2_-1_O")  {return "O2_-1_KO2              -0.5       -0.54500           -0.54350      -0.26900           -0.26970      -0.04900            -0.04680       -0.54320                    -0.54170                1.75600            1                 -                                   1           2                0.0      ";}
    else {return "";}                                                                                                      
  }

  //get_corrections_line_N////////////////////////////////////////////////////////
  // function to get corrections
  string get_corrections_line_N(const string& cor_identifier) {
    //For GaN Lany, PRB 78, 245207 2008 mentions "that the long-time tabulated value Hf =−1.10 eV (Refs. 10–12) for GaN has recently been corrected to −1.62 eV (Ref. 21)". 
    //However when I compare the values with my calculated LDA, PBE and SCAN data, the -1.10 eV/formula unit fits better to the overall fact that LDA and SCAN overestimate the formation enthalpy (too negative value) and PBE does the opposite. 
    //With the newer value of -1.62 eV/formula unit this would not be the case. Therefore I took the old tabulated value for deducing the corrections (January 07 2020).
    //                                          cation species            ox. state  PBE (for 298.15K)  PBE (for 0K)  LDA (for 298.15K)  LDA (for 0K)  SCAN (for 298.15K)  SCAN (for 0K)  PBE(+U)_ICSD (for 298.15K)  PBE(+U)_ICSD (for 0K)  Exp. (for 298.15K)  M-X bonds       Protos                               n cations   n anions        H_f^exp (for 298.15K)
    //                                                                               (eV/bond)          (eV/bond)     (eV/bond)          (eV/bond)     (eV/bond)           (eV/bond)      (eV/bond)                   (eV/bond)              (eV/bond)           per f. u.                                            per f. u.   per f. u.       (eV/bond) TABLE IV in supp.
    //NITRIDES                                                                                                                                                                                                                                                                                                                               
    if (cor_identifier=="Al_+3_N")       {return "Al_+3_AlN              +3          0.11875            0.10850      -0.04575           -0.05575      -0.03850            -0.04875        0.11800                     0.10825               -0.82500            4                AlN/AB_hP4_186_b_b-001.AB            1            1                0.0    ";}
    else if (cor_identifier=="B_+3_N")   {return "B_+3_BN_NJ             +3          0.00100           -0.00867      -0.15933           -0.17000      -0.11000            -0.12067        0.00167                    -0.00767               -0.86700            3                B_hN/AB_hP4_194_c_d-001.AB           1            1                0.0    ";}
    else if (cor_identifier=="Ca_+2_N")  {return "Ca_+2_Ca3N2_Ba         +2          0.03967            0.03400      -0.09225           -0.09792      -0.04392            -0.04942        0.03958                     0.03442               -0.37225            12               Ca_svN/A2B3_cI80_206_ad_e-001.BA     3            2                0.0    ";}
    else if (cor_identifier=="Li_+1_N")  {return "Li_+1_Li3N             +1          0.01663            0.00912      -0.06850           -0.07562      -0.03588            -0.04338        0.02725                     0.02013               -0.21350            8                Li_svN/A3B_hP4_191_bc_a-001.AB       3            1                0.0    ";}
    else if (cor_identifier=="Zn_+2_N")  {return "Zn_+2_Zn3N2            +2          0.06733            0.06225      -0.03492           -0.03908       0.00417            -0.00025       -0.06792                    -0.07242               -0.01950            12               NZn/A2B3_cI80_206_ad_e-001.AB        3            2                0.0    ";}
    else if (cor_identifier=="Be_+2_N")  {return "Be_+2_Be3N2            +2          0.05908            0.05292      -0.06500           -0.07133      -0.03433            -0.04075        0.05875                     0.05258               -0.50917            12               Be_svN/A2B3_cI80_206_ad_e-001.BA     3            2                0.0    ";}
    else if (cor_identifier=="Cr_+3_N")  {return "Cr_+3_CrN              +3          0.08417            0.08167      -0.02100           -0.02583      -0.01867            -0.02233       -0.22883                    -0.23283               -0.20250            6                Cr_pvN/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    else if (cor_identifier=="Ga_+3_N")  {return "Ga_+3_GaN_!!!          +3          0.04725            0.03950      -0.11725           -0.12550      -0.04200            -0.05075        0.01550                     0.00750               -0.28400            4                Ga_hN/AB_hP4_186_b_b-001.AB          1            1                0.0    ";}
    else if (cor_identifier=="Hf_+3_N")  {return "Hf_+3_HfN              +3          0.05983            0.05633      -0.09283           -0.09683      -0.01433            -0.01817        0.05933                     0.05583               -0.64533            6                Hf_pvN/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    else if (cor_identifier=="In_+3_N")  {return "In_+3_InN_Ba           +3          0.08825            0.08350      -0.05700           -0.06225      -0.01125            -0.01650        0.04800                     0.04325               -0.04450            4                In_dN/AB_hP4_186_b_b-001.AB          1            1                0.0    ";}
    else if (cor_identifier=="La_+3_N")  {return "La_+3_LaN_Ba           +3          0.08300            0.07983      -0.00100           -0.00433      -0.01600            -0.01950        0.13183                     0.12883               -0.52400            6                LaN/AB_cF8_225_a_b.AB                1            1                0.0    ";}
    else if (cor_identifier=="Mg_+2_N")  {return "Mg_+2_Mg3N2            +2          0.08300            0.07642      -0.03475           -0.04167      -0.01242            -0.01950        0.08083                     0.07425               -0.39858            12               Mg_pvN/A2B3_cI80_206_ad_e-001.BA     3            2                0.0    ";}
    else if (cor_identifier=="Nb_+3_N")  {return "Nb_+3_NbN              +3          0.07467            0.06967      -0.06550           -0.07083       0.02667             0.02150        0.07400                     0.06900               -0.40833            6                NNb_sv/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    else if (cor_identifier=="Sc_+3_N")  {return "Sc_+3_ScN              +3         -0.09850           -0.10500      -0.23150           -0.23833      -0.20667            -0.21350       -0.07883                    -0.08517               -0.54200            6                NSc_sv/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    else if (cor_identifier=="Si_+4_N")  {return "Si_+4_Si3N4            +4          0.00000           -0.01300      -0.22408           -0.23733      -0.15183            -0.16508       -0.00058                    -0.01358               -0.64325            12               NSi/A4B3_hP14_173_bc_c-001.AB        3            4                0.0    ";}
    else if (cor_identifier=="Ta_+3_N")  {return "Ta_+3_TaN              +3          0.15283            0.14933       0.00617            0.00233       0.13100             0.12733        0.15650                     0.15317               -0.43583            6                NTa_pv/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    else if (cor_identifier=="Ti_+3_N")  {return "Ti_+3_TiN              +3          0.00367           -0.00283      -0.14550           -0.15233      -0.08733            -0.09400        0.03083                     0.02467               -0.58400            6                NTi_sv/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    else if (cor_identifier=="V_+3_N")   {return "V_+3_VN                +3          0.04500            0.03917      -0.09483           -0.10100      -0.01683            -0.02250        0.03283                     0.03100               -0.37650            6                NV_sv/AB_cF8_225_a_b.AB              1            1                0.0    ";}
    else if (cor_identifier=="Y_+3_N")   {return "Y_+3_YN                +3         -0.06283           -0.06750      -0.18267           -0.18783      -0.16950            -0.17450       -0.04933                    -0.05400               -0.51683            6                NY_sv/AB_cF8_225_a_b.AB              1            1                0.0    ";}
    else if (cor_identifier=="Zr_+3_N")  {return "Zr_+3_ZrN_NJ           +3          0.04250            0.03733      -0.09933           -0.10500      -0.04067            -0.04617        0.06783                     0.06267               -0.63100            6                NZr_sv/AB_cF8_225_a_b.AB             1            1                0.0    ";}
    // following corrections (marked with *) might be less accurate since they were fitted to less well validated experimental data; Elder93: https://doi.org/10.1021/cm00034a027; Vajenine07: https://doi.org/10.1021/ic700406q
    else if (cor_identifier=="Na_+1_N")  {return "Na_+1_Na3N_Vajenine07* +1          0.03850            0.03400      -0.04367           -0.04883      -0.01150            -0.01683        0.03850                     0.03400                0.11050            6                -                                    3            1                0.663  ";}
    else if (cor_identifier=="W_+6_N")   {return "W_+6_Na3WN3_Elder93*   +6         -0.02325           -0.03050      -0.08401           -0.10626      -0.02850            -0.03251       -0.00025                    -0.00750               -1.26100            4                -                                    1            3                -3.718 ";}
    else if (cor_identifier=="Mo_+5_N")  {return "Mo_+5_LiMoN2_Elder93*  +5          0.28370            0.28038       0.15583            0.15129       0.24888             0.24488        0.28625                     0.28304               -0.45333            6                -                                    1            2                -4.001 ";}
    // following corrections (marked with ^) are only approximate as they were obtained from previous ones by lifting the explicit oxidation state dependence of the corrections
    else if (cor_identifier=="W_+4_N")   {return "W_+4_Na3WN3_Elder93*^  +4         -0.02325           -0.03050      -0.08401           -0.10626      -0.02850            -0.03251       -0.00025                    -0.00750               -1.26100            4                -                                    1            3                -3.718 ";}
    else if (cor_identifier=="Mo_+4_N")  {return "Mo_+4_LiMoN2_Elder93*^ +4          0.28370            0.28038       0.15583            0.15129       0.24888             0.24488        0.28625                     0.28304               -0.45333            6                -                                    1            2                -4.001 ";}
    else {return "";}                                                                                                      
  }

  //get_Bader_templates////////////////////////////////////////////////////////
  // function to get Bader charges from the binaries used to determine the corrections
  string get_Bader_templates(const string& element) {
    //Bader charges of the elements obtained from binary oxides
    //for DFT+U corrections will only be included for the PBE(+U) calculations used to generate the AFLOW ICSD data
    //DFT+U calculations differ not only by the U value but also by other specifications of the implementation, see:
    //Kick, Reuter, and Oberhofer; Intricacies of DFT+U, Not Only in a Numeric Atom Centered Orbital Framework, JCTC (2019); DOI: 10.1021/acs.jctc.8b01211
    //cation species           #ox.  ox       Bader charge (e charges)
    //                         stat  nr     PBE      LDA     SCAN    PBE(+U)_ICSD 
    if (element=="Ag")  {return "1   +1   0.4794   0.4659   0.4992   0.4513";}
    if (element=="Al")  {return "1   +3   2.4800   2.4533   2.5188   2.4706";}
    if (element=="As")  {return "1   +5   2.5884   2.6167   2.7515   2.5858";}
    if (element=="B")   {return "1   +3   2.3417   2.3087   2.4076   2.3261";}
    if (element=="Ba")  {return "1   +2   1.4375   1.4110   1.4636   1.4534";}
    if (element=="Be")  {return "1   +2   1.7043   1.6823   1.7276   1.7159";}
    if (element=="Bi")  {return "1   +3   1.7693   1.7505   1.7962   1.7707";}
    if (element=="Ca")  {return "1   +2   1.4772   1.4152   1.5117   1.4828";}
    if (element=="Cd")  {return "1   +2   1.1479   1.0778   1.2298   1.1791";}
    if (element=="Co")  {return "1   +2   1.1974   1.1351   1.2440   1.2665";}
    if (element=="Cr")  {return "2   +3   1.7009   1.6240   1.7560   1.7671   +6   2.0198   1.9856   2.1249   2.0498";}
    if (element=="Cs")  {return "1   +1   0.6818   0.6747   0.7132   0.6845";}
    if (element=="Cu")  {return "2   +1   0.5287   0.5236   0.5464   0.5105   +2   0.9651   0.9614   1.0033   1.0437";}
    if (element=="Fe")  {return "2   +2   1.2937   1.2688   1.3440   1.3784   +3   1.6233   1.3824   1.7882   1.8411";}
    if (element=="Ga")  {return "1   +3   1.8610   1.8334   1.9507   1.8600";}
    if (element=="Ge")  {return "1   +4   2.3971   2.3777   2.5471   2.4253";}
    if (element=="Hf")  {return "1   +4   2.5659   2.4994   2.6054   2.5587";}
    if (element=="Hg")  {return "1   +2   0.8813   0.8887   0.9201   0.8884";}
    if (element=="In")  {return "1   +3   1.8191   1.7843   1.9155   1.8147";}
    if (element=="Ir")  {return "1   +4   1.6523   1.6448   1.6469   1.6083";}
    if (element=="K")   {return "1   +1   0.7458   0.7101   0.7660   0.7265";}
    if (element=="Li")  {return "1   +1   0.8162   0.7974   0.8223   0.8372";}
    if (element=="Mg")  {return "1   +2   1.6973   1.6729   1.7147   1.6821";}
    if (element=="Mn")  {return "2   +2   1.3286   1.2956   1.4272   1.3775   +4   1.8460   1.7701   1.9347   1.8942";}
    if (element=="Mo")  {return "2   +4   2.0898   2.0799   2.1917   2.0983   +6   2.6270   2.6499   2.7532   2.5918";}
    if (element=="Na")  {return "1   +1   0.7868   0.7676   0.8073   0.7888";}
    if (element=="Nb")  {return "1   +2   1.3151   1.3781   1.3329   1.2806";}
    if (element=="Ni")  {return "1   +2   1.1075   1.0731   1.1578   1.2148";}
    if (element=="Os")  {return "2   +4   1.8518   1.8489   1.9032   1.8615   +8   2.6008   2.5884   2.6683   2.5591";}
    if (element=="Pb")  {return "2   +2   1.1513   1.1301   1.1698   1.1493   +4   2.0173   2.0324   2.1245   2.0469";}
    if (element=="Pd")  {return "1   +2   0.8312   0.8326   0.8630   0.8242";}
    if (element=="Rb")  {return "1   +1   0.7401   0.7124   0.7634   0.7192";}
    if (element=="Re")  {return "2   +4   2.0264   2.0190   2.0771   2.0627   +6   2.9344   2.9106   3.0224   2.8782";}
    if (element=="Rh")  {return "1   +3   1.2622   1.2465   1.3048   1.2867";}
    if (element=="Ru")  {return "1   +4   1.7022   1.6972   1.7991   1.8402";}
    if (element=="Sb")  {return "2   +3   1.7754   1.7743   1.8537   1.7742   +5   2.7273   2.7153   2.8873   2.7311";}
    if (element=="Sc")  {return "1   +3   2.0251   1.9722   2.0868   2.1215";}
    if (element=="Se")  {return "1   +4   1.9132   1.9535   2.0425   1.8766";}
    if (element=="Si")  {return "1   +4   3.2043   3.1655   3.2758   3.1716";}
    if (element=="Sn")  {return "2   +2   1.1884   1.1704   1.2280   1.1716   +4   2.2780   2.2608   2.4336   2.2097";}
    if (element=="Sr")  {return "1   +2   1.4895   1.4274   1.5335   1.4767";}
    if (element=="Te")  {return "1   +4   2.2551   2.2661   2.3516   2.2452";}
    if (element=="Ti")  {return "3   +2   1.2335   1.2183   1.2637   1.2658   +3   1.9228   1.8736   1.9955   2.0175   +4   2.2327   2.1844   2.3356   2.3956";}
    if (element=="Tl")  {return "2   +1   0.5628   0.5488   0.5779   0.5640   +3   1.4739   1.4718   1.5568   1.4781";}
    if (element=="V")   {return "4   +2   1.4271   1.4437   1.4306   1.5312   +3   1.7739   1.7286   1.8327   1.8783   +4   2.0728   2.0077   2.1737   2.1641   +5   2.1742   2.1430   2.2781   2.1982";}
    if (element=="W")   {return "2   +4   2.2270   2.1974   2.2933   2.2501   +6   2.9727   2.9284   3.0812   2.9838";}
    if (element=="Y")   {return "1   +3   2.1451   2.1027   2.2107   2.1473";}
    if (element=="Zn")  {return "1   +2   1.2134   1.2091   1.2764   1.1739";}
    if (element=="Zr")  {return "1   +4   2.5731   2.5281   2.6697   2.5711";}
    // corrections for per- (Li2O2) and superoxides (KO2)
    if (element=="O")   {return "2   -1   -0.8505  -0.8275  -0.8534  -0.8765  -0.5 -0.4334  -0.4157  -0.4338  -0.4313";}
    else {return "";}                                                                                                      
  }

  //get_ref_enthalpy_shift_pbe_u_icsd////////////////////////////////////////////////////////
  // Since the corrections for the PBE+U:ICSD calculations in the AFLOW ICSD were fit to ref. enthalpies
  // calculated with the same Us as for the species in the compound, a shift needs to be applied to 
  // the PBE ref. enthalpies that are used for formation enthalpy calculations to bring them to the 
  // PBE+U:ICSD ref. enthalpy values. These shifts will be summed for all atoms in the structure as 
  // needed and included in the CCE correction
  double get_ref_enthalpy_shift_pbe_u_icsd(const string& element) {
    if (element=="Sc")          {return 1.1610235;} // AUID="4b48722b5ae8f9e9"
    if (element=="Ti")          {return 3.12868285;} // AUID="9cdd346558c528c6"
    if (element=="V")           {return 2.7806973;} // AUID="94da8dda3e65c762"
    if (element=="Cr")          {return 4.0654442;} // AUID="f020a26862c0a532"
    if (element=="Mn")          {return 2.5458603448275854;} // AUID="99be850476e2dfb3"
    if (element=="Fe")          {return 3.4640268;} // AUID="15e65f5e50047695"
    if (element=="Co")          {return 4.0265057;} // AUID="9c4653e42c9f870b"
    if (element=="Ni")          {return 3.7335372;} // AUID="ead949d4e6fa93d1"
    if (element=="Cu")          {return 1.887473;} // AUID="62263cf2c885649c"
    if (element=="Zn")          {return 0.80815757;} // AUID="d4ad14ef9da00329"
    if (element=="Ga")          {return 0.486032725;} // AUID="567d70c3d85e2c7f"
    if (element=="Nb")          {return 1.8050886;} // AUID="141a7bea5919e0b3"
    if (element=="Mo")          {return 2.7869557;} // AUID="444a5f723c4c0576"
    if (element=="Tc")          {return 3.3680085;} // AUID="c33964edb4f74ebf"
    if (element=="Ru")          {return 3.613056;} // AUID="a3d039d8d88a4a86"
    if (element=="Rh")          {return 3.436254;} // AUID="b1547f3221c0c660"
    if (element=="Pd")          {return 2.7531925;} // AUID="4af3f66c14a34cee"
    if (element=="Ag")          {return 2.40546596;} // AUID="eca65d7d992efb48"
    if (element=="Cd")          {return 0.47503727;} // AUID="3388a058d6fd344d"
    if (element=="In")          {return 0.1350879;} // AUID="5ce1eee07a5df3a7"
    if (element=="Sn")          {return 0.1359635;} // AUID="37c5d78c5699ee73"
    if (element=="Ta")          {return 1.772861;} // AUID="926ce9924f7a5fd2"
    if (element=="W")           {return 2.388796;} // AUID="cad81c4b63d2d96e"
    if (element=="Re")          {return 2.957064;} // AUID="d72276b27490b853"
    if (element=="Os")          {return 3.1945335;} // AUID="216ecda62b9b32bf"
    if (element=="Ir")          {return 2.8837526;} // AUID="9c41e654d31db16a"
    if (element=="Pt")          {return 2.4902537;} // AUID="59081d68fabd053c"
    if (element=="La")          {return 0.139627;} // AUID="8353806ee117d7fd"
    if (element=="Ce")          {return AUROSTD_NAN;} // AUID="6e846dc10a52b11b"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Pr")          {return AUROSTD_NAN;} // AUID="82dbd5ee06192a4c"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Nd")          {return AUROSTD_NAN;} // AUID="81811d7c33a85a8b"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Sm")          {return AUROSTD_NAN;} // AUID="e892f33edd3901b8"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Eu")          {return 0.3926978;} // AUID="9b801fd8fd52e116"
    if (element=="Gd")          {return AUROSTD_NAN;} // AUID="c0f01c59eb971d30"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Dy")          {return 0.00023385;} // AUID="0373199190a8abb7"; calc. with and without U gives basically same energy
    if (element=="Tm")          {return AUROSTD_NAN;} // AUID="766b44cce10b0315"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Yb")          {return 0.5301314;} // AUID="de1c22a384a8d945"
    if (element=="Lu")          {return AUROSTD_NAN;} // AUID="98b3a5d589b252b5"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="Th")          {return AUROSTD_NAN;} // AUID="78f7995fc479457d"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    if (element=="U")           {return AUROSTD_NAN;} // AUID="2f0fba94372190f4"; no groundstate_energy in aflow_xpseudopotentials_data.cpp yet
    else {return -1;}
  }
} // namespace cce

#endif // _AFLOW_CCE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2021              *
// *                                                                         *
// ***************************************************************************
