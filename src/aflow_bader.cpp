// ***************************************************************************
// *                                                                         *
// *              Aflow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
// aflow_bader.cpp
// functions written by KEVIN RASCH
// UPDATED by COREY OSES
// 2013: kevin.rasch@duke.edu
// 2015: corey.oses@duke.edu

#include "aflow.h"
#include "aflow_bader.h"

// ***************************************************************************//
// bader_functions::BaderCalc
// ***************************************************************************//
namespace bader_functions {
  string BaderCalc(aurostd::xoption vpflow) {  //CO //keep as non-const, non-reference
    ostringstream oss;
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    string usage_usage = "aflow --bader [bader_options] -D DIRECTORY";
    vector<string> usage_options;
    usage_options.push_back(usage_usage);
    usage_options.push_back("");
    usage_options.push_back("options:");
    usage_options.push_back("--usage");
    usage_options.push_back("--critical_points|--cp");
    usage_options.push_back("--calculate=|--calc=bader|voronoi");
    usage_options.push_back("--nocalculate=|nocalc=bader|voronoi");
    usage_options.push_back("--partition=|--part=neargrid|ongrid");
    usage_options.push_back("--refine_edge_method=|--rem=|--r=-1,-2,-3");
    usage_options.push_back("--reference=|--ref=CHGCAR FILE");
    usage_options.push_back("--vacuum=|--vac=off|auto|VACUUM DENSITY");
    usage_options.push_back("--terminate|--term=known|max");
    usage_options.push_back("--print_all=atom|bader");
    usage_options.push_back("--print_index=|--print_ind=atom|bader");
    usage_options.push_back("--print_select_atom=|--print_sel_atom=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
    usage_options.push_back("--print_select_bader=|--print_sel_bader=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
    usage_options.push_back("--print_sum_atom=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
    usage_options.push_back("--print_sum_bader=1,2,3...|1-5 [volume list (separated by commas) or range (separated by hyphen)]");
    usage_options.push_back("--quiet|--q");
    // output usage
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK USAGE" << endl;
    if(vpflow.flag("BADER::USAGE")) {
      init::ErrorOption( "--usage", "pflow::BADER()", usage_options);
      return oss.str();
    }
    //FIRST STOP FOR CALLS MADE FROM COMMAND LINE
    string bader_options;  //to store bader options
    //directory
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    FixDirectory(directory);
    if(!aurostd::FileExist(directory)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to locate " << directory << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    if(!Flags2BaderCommands(vpflow, bader_options, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to gather flags for bader command." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    if(!BaderCalc(vpflow, bader_options, directory, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to perform bader calculation as specified." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    if(vpflow.flag("BADER::QUIET")) { oss.str(""); }
    return oss.str();
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::BaderCalc
// ***************************************************************************//
namespace bader_functions {
  bool BaderCalc(aurostd::xoption& vpflow,
      const string& bader_options,
      const string& _directory,
      ostream& oss) {
    //SECOND STOP FOR CALLS MADE FROM COMMAND_LINE
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    string empty_string(__AFLOW_FUNC__.size(), ' ');
    string directory=_directory;  //CO20200624
    FixDirectory(directory);  //CO20200624

    // get species names and valence charges from OUTCAR
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOOKING FOR SUITABLE OUTCAR" << endl;
    string outcar_file;
    if(!BaderExtensionFound("OUTCAR", outcar_file, directory)) {  //replace in place
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Cannot find suitable OUTCAR (OUTCAR.static, OUTCAR, or compressed variants)." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << __AFLOW_FUNC__ << " using " << outcar_file << " to extract system name, species, and valence charges." << endl;
    xOUTCAR outcar(outcar_file);  //can handle compressed files, automatically pulls in system, species (strings), and vZVAL (vector<doubles>)
    //[OBSOLETE as it is automatic now]outcar.GetProperties();                       //we need SYSTEM

    // get number of each species and valence charges from POSCAR
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOOKING FOR SUITABLE POSCAR" << endl;
    string poscar_file;
    if(!BaderExtensionFound("POSCAR", poscar_file, directory)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Cannot find suitable POSCAR (POSCAR.static, POSCAR, or compressed variants)." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << __AFLOW_FUNC__ << " using " << poscar_file << " to extract atom count for each species." << endl;
    xstructure xstr_bader(poscar_file, IOVASP_POSCAR);  //automatically pulls in num_each_type (ints); can handle compressed files

    // get cutoffs and downsample_ratios from FLAGS
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK FOR JVXL COMMAND" << endl;
    vector<double> cutoffs;
    vector<int> downsample_ratios;
    vpflow.flag("BADER::JVXL_CYCLIC", FALSE);
    if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::JVXL_ALL_SPECIES\")=" << vpflow.flag("BADER::JVXL_ALL_SPECIES") << endl;
      vector<string> sets, tokens;
      string input = vpflow.getattachedscheme("BADER::JVXL_ALL_SPECIES");

      // case 1: cyclic variable
      // 2a:  CUTOFF1,CUTOFF2…[::DOWNSAMPLE]
      // 2b:  CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...]
      // 2c:  some combination of above
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK INPUT TYPE" << endl;
      if(aurostd::substring2bool(input, "::")) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CYCLIC FOUND" << endl;
        vpflow.flag("BADER::JVXL_CYCLIC",TRUE);
        aurostd::string2tokens(input, tokens, "::");
        if(tokens.size() < 1 || tokens.size() > 2) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input - number of tokens (::)." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Needs to be either: " << endl;
          oss << empty_string << "CUTOFF1,CUTOFF2…[::DOWNSAMPLE], or " << endl;
          oss << empty_string << "CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...], or" << endl;
          oss << empty_string << "some combination." << endl;
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        //make sure there's no mixing
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK THAT THERE'S NO MIXING BETWEEN CYCLIC AND SETS" << endl;
        for (uint i = 0; i < tokens.size(); i++) {
          if(aurostd::substring2bool(tokens.at(i), ":")) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input, cannot specify sets (:) and cyclic (::) parameters." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
            return FALSE;
          }
        }
        if(aurostd::substring2bool(tokens.at(0), ",")) {
          aurostd::string2tokens<double>(tokens.at(0), cutoffs, ",");
        } else {
          cutoffs.push_back(aurostd::string2utype<double>(tokens.at(0)));
        }
        if(tokens.size() == 2) {
          if(aurostd::substring2bool(tokens.at(1), ",")) {
            aurostd::string2tokens<int>(tokens.at(1), downsample_ratios, ",");
          } else {
            downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(1)));
          }
        }
        // case 2, sets, CUTOFF1[,DOWNSAMPLE1]:CUTOFF2[,DOWNSAMPLE2]:...
      } else if(aurostd::substring2bool(input, ":")) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " SETS FOUND" << endl;
        aurostd::string2tokens(input, sets, ":");
        for (uint i = 0; i < sets.size(); i++) {
          aurostd::string2tokens(sets.at(i), tokens, ",");
          if(tokens.size() < 1 || tokens.size() > 2) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input " << i + 1 << "." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Needs to be CUTOFF1[,DOWNSAMPLE1]:CUTOFF2[,DOWNSAMPLE2]:..." << endl;
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
            return FALSE;
          }
          cutoffs.push_back(aurostd::string2utype<double>(tokens.at(0)));
          if(tokens.size() == 2) { downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(1))); }
        }
        // case 3: single set
      } else {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " SINGLE SET FOUND" << endl;
        if(aurostd::substring2bool(input, ",")) {
          aurostd::string2tokens(input, tokens, ",");
          if(tokens.size() != 2) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input - number of tokens." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Needs to be CUTOFF[,DOWNSAMPLE]." << endl;
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
            return FALSE;
          }
          cutoffs.push_back(aurostd::string2utype<double>(tokens.at(0)));
          downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(1)));
        } else {
          cutoffs.push_back(aurostd::string2utype<double>(input));
        }
      }
    }
    //check if extracted successfully
    if(outcar.SYSTEM.empty()) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to extract system name from " << outcar_file << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << outcar_file << " may be corrupted." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE; //CO20180220
    }
    if(outcar.species.empty()) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to extract species types from " << outcar_file << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << outcar_file << " may be corrupted." << endl;
      oss << __AFLOW_FUNC__ << " OUTCAR may be corrupted." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE; //CO20180220
    }
    if(outcar.vZVAL.empty()) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to extract valence charges from " << outcar_file << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << outcar_file << " may be corrupted." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE; //CO20180220
    }
    if(xstr_bader.num_each_type.empty()) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to extract atom count for each species from " << outcar_file << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << poscar_file << " may be corrupted." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE; //CO20180220
    }
    //get vector for vZVAL, we need it for VVdiff later
    //it needs to be as big as the number of atoms in cell
    vector<double> vZVAL;
    for (uint i = 0; i < outcar.vZVAL.size(); i++) {
      for (uint j = 0; j < (uint)xstr_bader.num_each_type.at(i); j++) {
        vZVAL.push_back(outcar.vZVAL.at(i));
      }
    }
    string system_name=KBIN::ExtractSystemName(directory);
    //return BaderCalc(vpflow, bader_options, outcar.SYSTEM, system_name, outcar.species, xstr_bader.num_each_type, vZVAL, cutoffs, downsample_ratios, directory, oss);
    return BaderCalc(vpflow, bader_options, system_name, outcar.species, xstr_bader.num_each_type, vZVAL, cutoffs, downsample_ratios, directory, oss);
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::BaderCalc
// ***************************************************************************//
namespace bader_functions {
  bool BaderCalc(aurostd::xoption& vpflow,
      const string& bader_options,
      const string& prototype,
      const vector<string>& vspecies,
      const deque<int>& num_each_type,
      const vector<double>& vZVAL,
      const vector<double>& cutoffs,
      const vector<int>& downsample_ratios,
      const string& _directory,
      ostream& oss) {
    //MAIN FUNCTION
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    //CO20180220 moved up from below
    string execution_path = aurostd::getPWD(); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd");
    string directory=_directory;  //CO20200624
    FixDirectory(directory);
    oss << __AFLOW_FUNC__ << " working within " << directory << "." << endl;

    //debug
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO fast checks START - 20170613
    if(vspecies.size()!=num_each_type.size()){
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input incorrect, vspecies.size()!=num_each_type.size() (" << vspecies.size() << "!=" << num_each_type.size() << ")" << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    uint natoms=0;
    for(uint i=0;i<vspecies.size();i++){
      for(uint j=0;j<(uint)num_each_type[i];j++){
        natoms++;
      }
    }
    if(vZVAL.size()!=natoms){
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input incorrect, vZVAL.size()!=natoms (" << vZVAL.size() << "!=" << natoms << ")" << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    //CO fast checks END - 20170613

    vector<string> compressed_files, remove_files, move_files;
    string compressed_file;
    //compressed_files=files to compress later
    //move_files=files to move later
    //remove_files=files to remove, include directory

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //  Start Directory
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // translate "./" into a full path, might be done twice
    //[CO20180220 - moved up]string execution_path = aurostd::getPWD(); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd");
    //[CO20180220 - moved up]FixDirectory(directory);

    //NO LONGER AN ISSUE, CO
    //FORTRAN has a problem with strings/paths longer than 128 characters, check
    //if(directory.length()>=128) {
    //  oss << endl;
    //  oss << __AFLOW_FUNC__ << " ERROR: Directory string length longer than 128. This will cause issues in Fortran." << " "; //<< endl //CO20180502;
    //  oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
    //  oss << __AFLOW_FUNC__ << " Exiting." << endl;
    //  oss << endl;
    //  return FALSE;
    //} else if(directory.length()>100) {
    //  oss << __AFLOW_FUNC__ << " WARNING: Directory string length longer than 100. This may cause issues with Fortran." << endl;
    //}

    // announce the directory to work in
    // unfortunate location for this message if you are coming from command line, but this is central to both command line/library
    // so I need it here
    //[CO20180220 - moved up]oss << __AFLOW_FUNC__ << " working within " << directory << "." << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //  End Directory
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // store expected name of the VASP data files required for Bader analysis
    // in a string

    //populate required_files
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " GATHERING REQUIRED FILES" << endl;
    vector<string> required_files;  //set up vector with TEMPORARY_DECOMPRESSED_FILES that we move later, and a files to move container
    oss << __AFLOW_FUNC__ << " adding CHGCAR to required files." << endl;
    required_files.push_back("CHGCAR");  //might need to decompress
    //DEFAULT, take AECCAR0+AECCAR2=REFERENCE_FILE (aflow.CHGCAR_sum), otherwise get from FLAGS
    //get REFERENCE_FILE if specified, otherwise DEFAULT
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK REFERENCE FILE" << endl;
    string ref_file="", tmp_file="";
    bool tempfile_created=false;
    if(vpflow.flag("BADER::REFERENCE")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::REFERENCE\")=" << vpflow.flag("BADER::REFERENCE") << endl;
      ref_file = execution_path + "/" + vpflow.getattachedscheme("BADER::REFERENCE");  //might need to decompress
      oss << __AFLOW_FUNC__ << " looking to use " << ref_file << " as reference file." << endl;
      if((!aurostd::FileExist(ref_file)) || aurostd::FileEmpty(ref_file)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: The reference charge file " << ref_file << " isn't present/is empty." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
        return FALSE;
      }
      if(aurostd::IsCompressed(ref_file)) {
        oss << __AFLOW_FUNC__ << " decompressing " << ref_file << "." << endl;
        if(!aurostd::efile2tempfile(ref_file, ref_file, tempfile_created)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to decompress " << ref_file << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
          return FALSE;
        }
        oss << __AFLOW_FUNC__ << " created temporary file " << ref_file << "." << endl;
        if(tempfile_created){remove_files.push_back(ref_file);}
      }
      oss << __AFLOW_FUNC__ << " referencing " << ref_file << "." << endl;
    } else {
      required_files.push_back("AECCAR0");
      required_files.push_back("AECCAR2");
      oss << __AFLOW_FUNC__ << " looking to use AECCAR0 and AECCAR2 to generate reference file." << endl;
    }

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " SEARCHING FOR / DECOMPRESSING REQUIRED FILES" << endl;
    for (uint i = 0; i < required_files.size(); i++) {
      if(!BaderExtensionFound(required_files.at(i), required_files.at(i), directory)) {  //replace in place
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Cannot find required file " << required_files.at(i) << " (";
        oss << required_files.at(i) << ".static, " << required_files.at(i) << ", or compressed variants)." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
        return FALSE;
      }
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " FOUND " << required_files.at(i) << endl;
      if(aurostd::IsCompressed(required_files.at(i))) {
        oss << __AFLOW_FUNC__ << " decompressing " << required_files.at(i) << "." << endl;
        if(!aurostd::efile2tempfile(required_files.at(i), tmp_file, tempfile_created)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to decompress " << required_files.at(i) << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
          return FALSE;
        }
        oss << __AFLOW_FUNC__ << " created temporary file " << tmp_file << "." << endl;
        if(tempfile_created){remove_files.push_back(tmp_file);}
        if(aurostd::FileEmpty(tmp_file)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: " << required_files.at(i) << " is empty." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
          return FALSE;
        }
        required_files.at(i) = tmp_file;
      }
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " FOUND ALL REQUIRED FILES" << endl;

    //// Make sub-directory to hold output Bader files
    //string bader_directory=directory + "/BADER";
    //oss << __AFLOW_FUNC__ << " creating BADER sub-directory." << endl;
    //if(!aurostd::DirectoryMake(bader_directory)) {
    //oss << endl;
    //oss << __AFLOW_FUNC__ << " ERROR: Could not open sub-directory in " << directory << "." << " "; //<< endl //CO20180502;
    //oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
    //oss << __AFLOW_FUNC__ << " Exiting." << endl;
    //oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
    //aurostd::RemoveFile(remove_files);   //CLEAN UP BEFORE RETURNING FALSE
#endif
    //return FALSE;
    //}
    //else {
    //oss << __AFLOW_FUNC__ << " Saving Bader analysis data in " << bader_directory << "." << endl;
    //}

    //first get species_header
    string species_header;
    for (uint i = 0; i < vspecies.size(); i++) {
      species_header.append(vspecies.at(i));
      if(i < vspecies.size() - 1) { species_header.append(" "); }
    }
    //species_header=species_header.substr(0,species_header.size()-1);  //FIX ME

    if(!vpflow.flag("BADER::REFERENCE")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " BUILDING REFERENCE FILE" << endl;
      // build CHGCAR_sum to use as a reference when running Bader
      //ref_file=bader_directory + "/aflow.CHGCAR_sum";  // path and name for the charge density file
      ref_file = directory + "/aflow.CHGCAR_sum";  // path and name for the charge density file
      oss << __AFLOW_FUNC__ << " summing core and valence charge densities for reference file." << endl;
      if(!pflow::CHGSUM(species_header, required_files.at(1), required_files.at(2), ref_file, oss)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to sum core and valence charge densities." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
        return FALSE;
      }
      // check files exist
      if((!aurostd::FileExist(ref_file)) || aurostd::FileEmpty(ref_file)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: The reference charge file " << ref_file << " isn't present/is empty." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
        return FALSE;
      }
      compressed_files.push_back("aflow.CHGCAR_sum");
      oss << __AFLOW_FUNC__ << " referencing " << ref_file << "." << endl;
    }

    // run the bader code
    if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
      oss << __AFLOW_FUNC__ << " calling BADER code." << endl;
    } else {
      oss << __AFLOW_FUNC__ << " calling \"GRID BASED BADER ANALYSIS\" with" << endl;
      oss << "                             " << "bader " << bader_options << endl;
      oss << "                             " << required_files.at(0) << endl;  //CHGCAR file (decompressed)
      oss << "                             " << "-ref " << ref_file << "." << endl;
    }

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK BADER COMMAND" << endl;
    if(!aurostd::IsCommandAvailable("bader")) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: bader command NOT found." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
      return FALSE;
    }
    string BADER_RUN_DIRECTORY = aurostd::TmpDirectoryCreate("BADER");
    oss << __AFLOW_FUNC__ << " changing directories to " << BADER_RUN_DIRECTORY << " to execute bader command." << endl;
    string work_dir=aurostd::getPWD(); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd");
    chdir(BADER_RUN_DIRECTORY.c_str());
    //FORTRAN has a problem with strings/paths longer than 128 characters, so make everything local
    aurostd::CopyFile(required_files.at(0), "CHGCAR");
    aurostd::CopyFile(ref_file, "REFERENCE");
    stringstream BADER_CODE_OUTPUT;
    string bader_command="";
    //do critical points first, so normal ACF, AVF, etc. get overwritten by the usual bader command (with bader CHGCAR -ref SUM)
    if(vpflow.flag("BADER::CRITICAL_POINTS")) {
      bader_command=XHOST.command("bader") + " -vac auto -cp REFERENCE";
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " bader_command=\"" << bader_command << "\"" << endl;}
      BADER_CODE_OUTPUT << endl << aurostd::execute2string(bader_command) << endl;  //local variant, must be SUM_CHARGE file (sum of valence+core charges)
      //not sure if we expect a particular type of error
      //I'm guessing since the didn't exit for maxima in edge refinement (original bader command), we're okay
      if(aurostd::substring2bool(BADER_CODE_OUTPUT.str(), "ERROR: should be no new maxima in edge refinement")) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: bader command issued error" << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " BADER CODE FAILURE: should be no new maxima in edge refinement" << endl;
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
        oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
        chdir(work_dir.c_str());//directory.c_str());
        oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
          oss << __AFLOW_FUNC__ << " Please check." << endl;
          oss << endl;
        }
#endif
        //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
        for (uint i = 0; i < compressed_files.size(); i++) {
          //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
          compressed_file = directory + "/" + compressed_files.at(i);
          if(aurostd::FileExist(compressed_file)) {
            oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
            aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
          }
        }
        return FALSE;
      }
    }
    // we don't want output in library runs
    if(!vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
      oss << BADER_CODE_OUTPUT.str();
    }
    BADER_CODE_OUTPUT.str("");

    //BADER_CODE_OUTPUT << endl << aurostd::execute2string(XHOST.command("bader")+" "+bader_options+" "+required_files.at(0)+" -ref "+ref_file) << endl;
    bader_command=XHOST.command("bader") + " " + bader_options + " CHGCAR -ref REFERENCE";
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " bader_command=\"" << bader_command << "\"" << endl;}
    BADER_CODE_OUTPUT << endl << aurostd::execute2string(bader_command) << endl;  //local variant
    //BADER_CODE_OUTPUT << endl << aurostd::execute2string("~/bin/bader "+bader_options+" "+required_files.at(0)+" -ref "+ref_file) << endl;
    if(aurostd::substring2bool(BADER_CODE_OUTPUT.str(), "ERROR: should be no new maxima in edge refinement")) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: bader command issued error" << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " BADER CODE FAILURE: should be no new maxima in edge refinement" << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
      oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
      chdir(work_dir.c_str());//directory.c_str());
      oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
        oss << __AFLOW_FUNC__ << " Please check." << endl;
        oss << endl;
      }
#endif
      //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
      for (uint i = 0; i < compressed_files.size(); i++) {
        //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
        compressed_file = directory + "/" + compressed_files.at(i);
        if(aurostd::FileExist(compressed_file)) {
          oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
          aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
        }
      }
      return FALSE;
    }
    // we don't want output in library runs
    if(!vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
      oss << BADER_CODE_OUTPUT.str();
    }
    BADER_CODE_OUTPUT.str("");

#ifndef _AFLOW_TEMP_PRESERVE_
    //REMOVE FILES, round 1
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " removing files" << endl;
    for (uint i = 0; i < remove_files.size(); i++) {
      oss << __AFLOW_FUNC__ << " removing " << remove_files.at(i) << "." << endl;
      aurostd::RemoveFile(remove_files.at(i));
    }
    remove_files.clear();
#endif

    //no bader files expected, dumb case but it's there
    if(vpflow.flag("BADER::NOCALCULATE") && vpflow.getattachedscheme("BADER::NOCALCULATE") == "bader") {
      oss << __AFLOW_FUNC__ << " no bader files expected." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
      oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
      chdir(work_dir.c_str());//directory.c_str());
      oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Please check." << endl;
        oss << endl;
      }
#endif
      //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
      for (uint i = 0; i < compressed_files.size(); i++) {
        //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
        compressed_file = directory + "/" + compressed_files.at(i);
        if(aurostd::FileExist(compressed_file)) {
          oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
          aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
        }
      }
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      //oss << __AFLOW_FUNC__ << " Exiting." << endl;
      return TRUE;
    }

    //GET LIST OF FILES CREATED TO MOVE AND BZIP
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " GET LIST OF FILES CREATED TO MOVE AND BZIP" << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK STANDARD BADER FILES" << endl;
    vector<string> standard_bader_files;
    standard_bader_files.push_back("ACF.dat");
    standard_bader_files.push_back("AVF.dat");
    standard_bader_files.push_back("BCF.dat");
    for (uint i = 0; i < standard_bader_files.size(); i++) {
      //created in current directory
      if(aurostd::FileExist(standard_bader_files.at(i)) && aurostd::FileNotEmpty(standard_bader_files.at(i))) {
        compressed_files.push_back(standard_bader_files.at(i));
        move_files.push_back(standard_bader_files.at(i));
        oss << __AFLOW_FUNC__ << " file created: " << standard_bader_files.at(i) << "." << endl;
      } else {
        //echo warning, but don't exist
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: " << standard_bader_files.at(i) << " was not created/is empty." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
    }
    //critical points file
    if(vpflow.flag("BADER::CRITICAL_POINTS")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK CRITICAL POINTS FILE" << endl;
      if(aurostd::FileExist("CPF.dat") && aurostd::FileNotEmpty("CPF.dat")) {
        compressed_files.push_back("CPF.dat");
        move_files.push_back("CPF.dat");
        oss << __AFLOW_FUNC__ << " file created: CPF.dat." << endl;
      } else {
        //echo warning, but don't exist
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: CPF.dat was not created/is empty." << endl;
        oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
    }
    vector<int> range;
    int atomCOUNT = 1;
    string bader_file;
    //we don't want duplicate files, and print_all_atom and print_select_atom will produce duplicates
    //we also don't want to REedit files
    if(vpflow.flag("BADER::PRINT_ALL") || vpflow.flag("BADER::PRINT_SELECT_ATOM") || vpflow.flag("BADER::PRINT_SELECT_BADER")) {
      //BvAtxxxx.dat
      //know this from vpspecies and num_each_type
      if(vpflow.getattachedscheme("BADER::PRINT_ALL") == "atom" || vpflow.getattachedscheme("BADER::PRINT_ALL") == "both") {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT ALL ATOM FILES" << endl;
        atomCOUNT = 1;
        for (uint i = 0; i < vspecies.size(); i++) {
          for (uint speciesCOUNT = 0; speciesCOUNT < (uint)num_each_type.at(i); speciesCOUNT++) {
            bader_file = "BvAt" + aurostd::PaddedNumString(atomCOUNT, 4) + ".dat";
            if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
              compressed_files.push_back(bader_file);
              move_files.push_back(bader_file);
              oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
              //if we keep these files, then edit their header
              if(!vpflow.flag("BADER::REMOVE_BADER_ATOMS")) {
                if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
                  oss << endl;
                  oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
                  oss << __AFLOW_FUNC__ << " Please check." << endl;
                  oss << endl;
                }
              }
            } else {
              //echo warning, do not exit
              oss << endl;
              oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
              oss << __AFLOW_FUNC__ << " This is VERY unusual." << endl;
              oss << endl;
            }
            atomCOUNT++;
          }
        }
        //LIST/RANGE
      } else if(vpflow.flag("BADER::PRINT_SELECT_ATOM")) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING SELECT ATOM FILES" << endl;
        if(!listORrange2vec(vpflow.getattachedscheme("BADER::PRINT_SELECT_ATOM"), range, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to extract list/range command for print select atom." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
          oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
          chdir(work_dir.c_str());//directory.c_str());
          oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
            oss << __AFLOW_FUNC__ << " Please check." << endl;
            oss << endl;
          }
#endif
          //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
          for (uint i = 0; i < compressed_files.size(); i++) {
            //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
            compressed_file = directory + "/" + compressed_files.at(i);
            if(aurostd::FileExist(compressed_file)) {
              oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
              aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
            }
          }
          return FALSE;
        }
        for (uint i = 0; i < range.size(); i++) {
          bader_file = "BvAt" + aurostd::PaddedNumString(range.at(i), 4) + ".dat";
          if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
            compressed_files.push_back(bader_file);
            move_files.push_back(bader_file);
            oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
            if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
              oss << endl;
              oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
              oss << __AFLOW_FUNC__ << " Please check." << endl;
              oss << endl;
            }
          } else {
            //echo warning, do not exit
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
            oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
            oss << endl;
          }
        }
      }
      //Bvolxxxx.dat
      //no way to know unless we check what bader spits out
      if(vpflow.getattachedscheme("BADER::PRINT_ALL") == "bader" || vpflow.getattachedscheme("BADER::PRINT_ALL") == "both") {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT ALL BADER FILES" << endl;
        atomCOUNT = 1;
        //check that first exists, other throw an error
        bader_file = "Bvol" + aurostd::PaddedNumString(atomCOUNT, 4) + ".dat";
        if((!aurostd::FileExist(bader_file)) || aurostd::FileEmpty(bader_file)) {
          //echo warning, but do not exist
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: First expected Bvolxxxx.dat file NOT FOUND/IS EMPTY." << endl;
          oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
          oss << endl;
        }
        while (aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
          compressed_files.push_back(bader_file);
          move_files.push_back(bader_file);
          oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
          if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
            oss << __AFLOW_FUNC__ << " Please check." << endl;
            oss << endl;
          }
          bader_file = "Bvol" + aurostd::PaddedNumString(++atomCOUNT, 4) + ".dat";  //CO20200404 - add BEFORE
        }
        //B_wexxxx.dat
        //these have different file names than Bvolxxx.dat even though they are the same files (artifact of bader code?)
        //LIST/RANGE
      } else if(vpflow.flag("BADER::PRINT_SELECT_BADER")) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT SELECT BADER FILES" << endl;
        if(!listORrange2vec(vpflow.getattachedscheme("BADER::PRINT_SELECT_BADER"), range, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to extract list/range command for print select bader." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
          oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
          chdir(work_dir.c_str());//directory.c_str());
          oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
            oss << __AFLOW_FUNC__ << " Please check." << endl;
            oss << endl;
          }
#endif
          //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
          for (uint i = 0; i < compressed_files.size(); i++) {
            //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
            compressed_file = directory + "/" + compressed_files.at(i);
            if(aurostd::FileExist(compressed_file)) {
              oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
              aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
            }
          }
          return FALSE;
        }
        for (uint i = 0; i < range.size(); i++) {
          bader_file = "B_we" + aurostd::PaddedNumString(range.at(i), 4) + ".dat";
          if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
            compressed_files.push_back(bader_file);
            move_files.push_back(bader_file);
            oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
            if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
              oss << endl;
              oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
              oss << __AFLOW_FUNC__ << " Please check." << endl;
              oss << endl;
            }
          } else {
            //echo warning, do not exit
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
            oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
            oss << endl;
          }
        }
      }
    }
    //BvAt_summed.dat
    if(vpflow.flag("BADER::PRINT_SUM_ATOM")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT SUM ATOM FILE" << endl;
      bader_file = "BvAt_summed.dat";
      if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
        compressed_files.push_back(bader_file);
        move_files.push_back(bader_file);
        oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
        if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
          oss << __AFLOW_FUNC__ << " Please check." << endl;
          oss << endl;
        }
      } else {
        //echo warning, do not exit
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
        oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
    }
    //Bvol_summed.dat
    if(vpflow.flag("BADER::PRINT_SUM_BADER")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT SUM BADER FILE" << endl;
      bader_file = "Bvol_summed.dat";
      if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
        compressed_files.push_back(bader_file);
        move_files.push_back(bader_file);
        oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
        if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
          oss << __AFLOW_FUNC__ << " Please check." << endl;
          oss << endl;
        }
      } else {
        //echo warning, do not exit
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
        oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
        oss << endl;
      }
    }
    if(vpflow.flag("BADER::PRINT_INDEX")) {
      //AtIndex.dat
      if(vpflow.getattachedscheme("BADER::PRINT_INDEX") == "atom" || vpflow.getattachedscheme("BADER::PRINT_INDEX") == "both") {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT ATOM INDEX FILE" << endl;
        bader_file = "AtIndex.dat";
        if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
          compressed_files.push_back(bader_file);
          move_files.push_back(bader_file);
          oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
          if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
            oss << __AFLOW_FUNC__ << " Please check." << endl;
            oss << endl;
          }
        } else {
          //echo warning, do not exit
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
          oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
          oss << endl;
        }
      }
      //BvIndex.dat
      if(vpflow.getattachedscheme("BADER::PRINT_INDEX") == "bader" || vpflow.getattachedscheme("BADER::PRINT_INDEX") == "both") {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECKING PRINT BADER INDEX FILE" << endl;
        bader_file = "BvIndex.dat";
        if(aurostd::FileExist(bader_file) && aurostd::FileNotEmpty(bader_file)) {
          compressed_files.push_back(bader_file);
          move_files.push_back(bader_file);
          oss << __AFLOW_FUNC__ << " file created: " << bader_file << "." << endl;
          if(!prepare_CHGCAR_4_Jmol(bader_file, species_header, FALSE, oss)) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " WARNING: Unable to fix header for " << bader_file << "." << endl;
            oss << __AFLOW_FUNC__ << " Please check." << endl;
            oss << endl;
          }
        } else {
          //echo warning, do not exit
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: " << bader_file << " was not created/is empty." << endl;
          oss << __AFLOW_FUNC__ << " This is VERY unusual and likely an issue with the Bader analysis." << endl;
          oss << endl;
        }
      }
    }
    //move created files into correct directory, we can now work exclusively in this directory
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " MOVING FILES CREATED BY BADER TO DIRECTORY" << endl;
    for (uint i = 0; i < move_files.size(); i++) {
      oss << __AFLOW_FUNC__ << " moving " << move_files.at(i) << " to " << directory << "." << endl;
      if(!aurostd::file2directory(move_files.at(i), directory)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to move " << move_files.at(i) << "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
        oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
        chdir(work_dir.c_str());//directory.c_str());
        oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
        if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {  //CLEAN UP BEFORE RETURNING FALSE
          oss << endl;
          oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
          oss << __AFLOW_FUNC__ << " Please check." << endl;
          oss << endl;
        }
#endif
        //BZIP aflow.CHGCAR_sum (only file created), and it's created in the right directory, so no need to move
        for (uint i = 0; i < compressed_files.size(); i++) {
          //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
          compressed_file = directory + "/" + compressed_files.at(i);
          if(aurostd::FileExist(compressed_file)) {
            oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
            aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
          }
        }
        return FALSE;
      }
    }
    oss << __AFLOW_FUNC__ << " changing directories back to original directory " << directory << "." << endl;
    chdir(work_dir.c_str());//directory.c_str());
    oss << __AFLOW_FUNC__ << " removing temporary directory " << BADER_RUN_DIRECTORY << "." << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
    if(!aurostd::RemoveDirectory(BADER_RUN_DIRECTORY)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " WARNING: issues deleting directory " << BADER_RUN_DIRECTORY << "." << endl;
      oss << __AFLOW_FUNC__ << " Please check." << endl;
      oss << endl;
    }
#endif

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " READING ACF.DAT FILE FOR BADER VOLUME AND CHARGE" << endl;
    // containers for AFLOW net charge data
    vector<string> lines, tokens;
    vector<double> bader_charge, volume, net_charge;

    // get bader charge and atomic volume results for each atom from ACF.dat
    //if(aurostd::FileExist(bader_directory+"/ACF.dat"))
    if(aurostd::FileExist(directory + "/ACF.dat") && aurostd::FileNotEmpty(directory + "/ACF.dat"))
    { //CO20200106 - patching for auto-indenting
      //aurostd::efile2vectorstring(bader_directory+"/ACF.dat",lines);
      aurostd::efile2vectorstring(directory + "/ACF.dat", lines);
      for (uint i = 0; i < lines.size(); i++) {
        aurostd::string2tokens(lines.at(i), tokens, " ");
        if(tokens.size() == 7 && tokens.at(4) != "CHARGE") {
          bader_charge.push_back(aurostd::string2utype<double>(tokens.at(4)));
          volume.push_back(aurostd::string2utype<double>(tokens.at(6)));
        }
      }
    } else {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to locate ACF.dat (or its empty)." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveFile(remove_files);
#endif
      //BZIP those that remain
      for (uint i = 0; i < compressed_files.size(); i++) {
        //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
        compressed_file = directory + "/" + compressed_files.at(i);
        if(aurostd::FileExist(compressed_file)) {
          oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
          aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
        }
      }
      return FALSE;
    }

    if(bader_charge.size() == 0 || volume.size() == 0 || bader_charge.size() != volume.size()) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to extract bader_charge/volume correctly." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Perhaps the file is corrupt?" << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      aurostd::RemoveFile(remove_files);
#endif
      //BZIP those that remain
      for (uint i = 0; i < compressed_files.size(); i++) {
        //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
        compressed_file = directory + "/" + compressed_files.at(i);
        if(aurostd::FileExist(compressed_file)) {
          oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
          aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
        }
      }
      return FALSE;
    }

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " TAKING DIFFERENCE BETWEEN VALENCE CHARGE AND BADER CHARGE" << endl;
    // calculate net charge on each atom
    // vZVAL is read from POTCAR/OUTCAR and is the number of valence electrons used by that particular pseudopotential (core ignored)
    // bader_charge is the charge (population) found in the volume
    // (vZVAL - bader_charge) yields residual (net) charge, positive if charge has been removed, negative if charge has been added
    // we only care about valence here, but we use the total charge (core + valence) as a reference for location of maxima (located at atomic centers)
    // bader still provides population wrt CHGCAR (valence ONLY)
    net_charge = pflow::VVdiff(vZVAL, bader_charge);

    // make sure it's zero! otherwise spit warning
    double net_sum = 0;
    for (uint i = 0; i < net_charge.size(); i++) {
      net_sum += net_charge.at(i);
    }
    if(abs(net_sum) > 1e-3) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " WARNING: The sum of the partial charges is non-zero (sum=" << net_sum << ",tol=1e-3)." << endl;
      oss << endl;
    }

    // write net charge on each atom to file abader.out
    string abader_out = prototype + "_abader.out";
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " WRITING RESULTS TO " << abader_out << endl;
    string net_charges_string = "bader_net_charges=";
    string volume_string = "bader_atomic_volumes=";
    string abader_spacer = "[AFLOW] **************************************************************************************************************************";
    //abader_out works for both command line users (non-library use) and library use
    stringstream abader_ss;

    oss << __AFLOW_FUNC__ << " writing net charges to " << abader_out << "." << endl;
    string column1 = "Atom #";
    string column2 = "Atom Type";
    string column3 = "Net Charge (electrons)";
    string column4 = "Volume (Angst^3)";
    string spaces = "    ";  //4 spaces
    if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
      abader_ss << abader_spacer << endl;
      abader_ss << "[BADER_RESULTS]START" << endl;
    } else {
      abader_ss << column1 << spaces << column2 << spaces << column3 << spaces << column4 << endl;  // header of labels
    }

    stringstream num_prec;
    atomCOUNT = 0;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vspecies().size()=" << vspecies.size() << endl;
    for (uint i = 0; i < vspecies.size(); i++) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " (uint)num_each_type.at(i=" << i << ")=" << (uint)num_each_type.at(i) << endl;
      for (uint speciesCOUNT = 0; speciesCOUNT < (uint)num_each_type.at(i); speciesCOUNT++) {
        if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
          num_prec.str("");
          num_prec << std::fixed << setprecision(4) << net_charge.at(atomCOUNT);
          net_charges_string.append(num_prec.str());
          num_prec.str("");
          num_prec << std::fixed << setprecision(4) << volume.at(atomCOUNT);
          volume_string.append(num_prec.str());
          if(!((i == vspecies.size() - 1) && (speciesCOUNT == (uint)num_each_type.at(i) - 1))) {
            net_charges_string.append(",");
            volume_string.append(",");
          }
        } else {
          abader_ss << aurostd::PaddedPRE(aurostd::utype2string(atomCOUNT + 1), column1.length(), " ") << spaces;  //right align
          abader_ss << aurostd::PaddedPOST(vspecies.at(i), column2.length(), " ") << spaces;                       //left align
          num_prec.str("");
          num_prec << std::fixed << setprecision(4) << net_charge.at(atomCOUNT);
          abader_ss << aurostd::PaddedPRE(num_prec.str(), column3.length(), " ") << spaces;  //right align
          num_prec.str("");
          num_prec << std::fixed << setprecision(4) << volume.at(atomCOUNT);
          abader_ss << aurostd::PaddedPRE(num_prec.str(), column4.length(), " ") << endl;  //right align
        }
        atomCOUNT++;
      }
    }

    //FOR AFLOWLIB LABELS
    //vpflow.flag("BADER::AFLOWLIB_LIBRARY",TRUE);
    if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
      //units
      net_charges_string.append(" (electrons)");
      volume_string.append(" (Angst^3)");
      abader_ss << net_charges_string << endl;
      abader_ss << volume_string << endl;
      abader_ss << "[BADER_RESULTS]STOP" << endl;
      abader_ss << abader_spacer << endl;
      aurostd::stringstream2file(abader_ss, directory + "/" + abader_out);
    } else {
      //aurostd::stringstream2file(abader_ss,bader_directory+"/"+abader_out);
      aurostd::stringstream2file(abader_ss, directory + "/" + abader_out);
      compressed_files.push_back(abader_out);  //DO NOT ZIP FOR LIBRARY
    }

    //consolidate, if specified by FLAG
    string species_file;
    if(vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) {
      //use sum_files
      vector<string> species_files, bader_atom_files, sum_files;
      atomCOUNT = 1;
      oss << __AFLOW_FUNC__ << " attempting to add individual atom contributions together for each species." << endl;
      //cycle through to CONSOLIDATE_ATOMS3SPECIES
      for (uint i = 0; i < vspecies.size(); i++) {
        species_file = "BvAt_" + vspecies.at(i) + ".dat";
        compressed_files.push_back(species_file);
        species_files.push_back(directory + "/" + species_file);
        //species_files.push_back(bader_directory+"/"+species_file);
        for (uint speciesCOUNT = 0; speciesCOUNT < (uint)num_each_type.at(i); speciesCOUNT++) {
          bader_file = "BvAt" + aurostd::PaddedNumString(atomCOUNT, 4) + ".dat";
          //bader_atom_files.push_back(bader_directory+"/"+bader_file);
          bader_atom_files.push_back(directory + "/" + bader_file);
          sum_files.push_back(bader_atom_files.back());
          if((!aurostd::FileExist(bader_atom_files.back())) || aurostd::FileEmpty(bader_atom_files.back())) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: " << bader_atom_files.back() << " not found/is empty." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
            aurostd::RemoveFile(remove_files);
#endif
            //BZIP those that remain
            for (uint i = 0; i < compressed_files.size(); i++) {
              //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
              compressed_file = directory + "/" + compressed_files.at(i);
              if(aurostd::FileExist(compressed_file)) {
                oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
                aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
              }
            }
            return FALSE;
          }
          atomCOUNT++;
        }
        if(!pflow::CHGSUM(species_header, sum_files, species_files.back(), oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to sum atom files to create " << species_files.back() << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
          aurostd::RemoveFile(remove_files);
#endif
          //BZIP those that remain
          for (uint i = 0; i < compressed_files.size(); i++) {
            //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
            compressed_file = directory + "/" + compressed_files.at(i);
            if(aurostd::FileExist(compressed_file)) {
              oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
              aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
            }
          }
          return FALSE;
        }
        oss << __AFLOW_FUNC__ << " file created: " << species_files.back() << "." << endl;
        sum_files.clear();
      }
      //if REMOVE_BADER_ATOMS
      if(vpflow.flag("BADER::REMOVE_BADER_ATOMS")) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " REMOVE BADER ATOMS FLAG SPECIFIED" << endl;
        remove_files.insert(remove_files.end(), bader_atom_files.begin(), bader_atom_files.end());
      }
      if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
        oss << __AFLOW_FUNC__ << " attempting to create .jvxl files for visualization in JMOL." << endl;
        string output_file, cutoff_string;
        bool cyclic = FALSE;
        if(vpflow.flag("BADER::JVXL_CYCLIC")) { cyclic = TRUE; }
        //if(downsample_ratios.size()==0) cyclic=FALSE; //OBVIOUS
        //check that they are same size if not cyclic
        if(!cyclic) {
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " PERFORMING SETS ROUTINE WITH CUTOFFS AND DOWNSAMPLE RATIOS" << endl;
          //cutoffs and downsample_ratios share same index
          if(downsample_ratios.size() != 0) {
            if(LDEBUG) cerr << __AFLOW_FUNC__ << " NO DOWNSAMPLE RATIOS SPECIFIED" << endl;
            if(cutoffs.size() != downsample_ratios.size()) {
              oss << endl;
              oss << __AFLOW_FUNC__ << " ERROR: Number of cutoffs != number of downsample ratios." << " "; //<< endl //CO20180502;
              oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
              oss << __AFLOW_FUNC__ << " This is necessary for sets." << endl;
              oss << __AFLOW_FUNC__ << " Exiting." << endl;
              oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
              aurostd::RemoveFile(remove_files);  //CLEAN UP BEFORE RETURNING FALSE
#endif
              //BZIP those that remain
              for (uint i = 0; i < compressed_files.size(); i++) {
                //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
                compressed_file = directory + "/" + compressed_files.at(i);
                if(aurostd::FileExist(compressed_file)) {
                  oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
                  aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
                }
              }
              return FALSE;
            }
          }
          for (uint i = 0; i < species_files.size(); i++) {
            for (uint j = 0; j < cutoffs.size(); j++) {
              //get cutoff string, 25 instead of 0.25
              num_prec.str("");
              num_prec << std::fixed << setprecision(2) << cutoffs.at(j);  //avoid issues with setprecision(20) in string2utype<double>
              cutoff_string = num_prec.str();
              //cutoff_string=aurostd::utype2string(cutoffs.at(j));
              if(aurostd::substring2bool(cutoff_string, "0.")) { cutoff_string = aurostd::StringSubst(cutoff_string, "0.", ""); }
              if(aurostd::substring2bool(cutoff_string, ".")) { cutoff_string = aurostd::StringSubst(cutoff_string, ".", ""); }
              //output_file=bader_directory+"/"+prototype+"_Bader_"+cutoff_string+"_";
              output_file = directory + "/" + prototype + "_Bader_" + cutoff_string + "_";
              if(downsample_ratios.size() != 0 && !(vpflow.flag("BADER::AFLOWLIB_LIBRARY"))) {
                output_file.append(aurostd::utype2string(downsample_ratios.at(j)) + "_");
              }
              output_file.append(vspecies.at(i) + ".jvxl");
              oss << __AFLOW_FUNC__ << " output file will be " << output_file << "." << endl;
              if(downsample_ratios.size() != 0) {
                if(!pflow::CHGCAR2JVXL(species_files.at(i), cutoffs.at(j), downsample_ratios.at(j), output_file, oss)) {
                  oss << endl;
                  oss << __AFLOW_FUNC__ << " ERROR: Unable to create " << output_file << "." << " "; //<< endl //CO20180502;
                  oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
                  oss << __AFLOW_FUNC__ << " Exiting." << endl;
                  oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
                  aurostd::RemoveFile(remove_files);
#endif
                  //BZIP those that remain
                  for (uint i = 0; i < compressed_files.size(); i++) {
                    //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
                    compressed_file = directory + "/" + compressed_files.at(i);
                    if(aurostd::FileExist(compressed_file)) {
                      oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
                      aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
                    }
                  }
                  return FALSE;
                }
              } else {
                if(!pflow::CHGCAR2JVXL(species_files.at(i), cutoffs.at(j), output_file, oss)) {
                  oss << endl;
                  oss << __AFLOW_FUNC__ << " ERROR: Unable to create " << output_file << "." << " "; //<< endl //CO20180502;
                  oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
                  oss << __AFLOW_FUNC__ << " Exiting." << endl;
                  oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
                  aurostd::RemoveFile(remove_files);
#endif
                  //BZIP those that remain
                  for (uint i = 0; i < compressed_files.size(); i++) {
                    //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
                    compressed_file = directory + "/" + compressed_files.at(i);
                    if(aurostd::FileExist(compressed_file)) {
                      oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
                      aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
                    }
                  }
                  return FALSE;
                }
              }
            }
          }
        } else {
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " PERFORMING CYCLIC ROUTINE WITH CUTOFFS AND DOWNSAMPLE RATIOS" << endl;
          //if cyclic specified, it must have downsample_ratios, what else is it going to cycle over?
          for (uint i = 0; i < species_files.size(); i++) {
            for (uint j = 0; j < cutoffs.size(); j++) {
              for (uint k = 0; k < downsample_ratios.size(); k++) {
                //get cutoff string, 25 instead of 0.25
                num_prec.str("");
                num_prec << std::fixed << setprecision(2) << cutoffs.at(j);  //avoid issues with setprecision(20) in string2utype<double>
                cutoff_string = num_prec.str();
                //cutoff_string=aurostd::utype2string(cutoffs.at(j));
                if(aurostd::substring2bool(cutoff_string, "0.")) { cutoff_string = aurostd::StringSubst(cutoff_string, "0.", ""); }
                if(aurostd::substring2bool(cutoff_string, ".")) { cutoff_string = aurostd::StringSubst(cutoff_string, ".", ""); }
                //output_file=bader_directory+"/"+prototype+"_Bader_"+cutoff_string+"_"
                output_file = directory + "/" + prototype + "_Bader_" + cutoff_string + "_";
                if(!(vpflow.flag("BADER::AFLOWLIB_LIBRARY"))) {
                  output_file.append(aurostd::utype2string(downsample_ratios.at(k)) + "_");
                }
                output_file.append(vspecies.at(i) + ".jvxl");
                oss << __AFLOW_FUNC__ << " output file will be " << output_file << "." << endl;
                if(!pflow::CHGCAR2JVXL(species_files.at(i), cutoffs.at(j), downsample_ratios.at(k), output_file, oss)) {
                  oss << endl;
                  oss << __AFLOW_FUNC__ << " ERROR: Unable to create " << output_file << "." << " "; //<< endl //CO20180502;
                  oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
                  oss << __AFLOW_FUNC__ << " Exiting." << endl;
                  oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
                  aurostd::RemoveFile(remove_files);
#endif
                  //BZIP those that remain
                  for (uint i = 0; i < compressed_files.size(); i++) {
                    //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
                    compressed_file = directory + "/" + compressed_files.at(i);
                    if(aurostd::FileExist(compressed_file)) {
                      oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
                      aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
                    }
                  }
                  return FALSE;
                }
              }
            }
          }
        }
        //if KEEP::JVXL_ONLY
        if(vpflow.flag("BADER::KEEP::JVXL_ONLY")) {
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " KEEP JVXL FILES ONLY FLAG SPECIFIED" << endl;
          remove_files.insert(remove_files.end(), species_files.begin(), species_files.end());
        }
      }
    }
#ifndef _AFLOW_TEMP_PRESERVE_
    //REMOVE FILES, round 2
    for (uint i = 0; i < remove_files.size(); i++) {
      oss << __AFLOW_FUNC__ << " removing " << remove_files.at(i) << "." << endl;
      aurostd::RemoveFile(remove_files.at(i));
    }
    remove_files.clear();
#endif

    //BZIP those that remain
    for (uint i = 0; i < compressed_files.size(); i++) {
      //if(aurostd::FileExist(bader_directory+"/"+compressed_files.at(i))) {  //[CO20200106 - close bracket for indenting]}
      compressed_file = directory + "/" + compressed_files.at(i);
      if(aurostd::FileExist(compressed_file)) {
        oss << __AFLOW_FUNC__ << " compressing " << compressed_files.at(i) << "." << endl;
        aurostd::CompressFile(compressed_file,DEFAULT_KZIP_BIN);
      }
    }
    return TRUE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::BaderCalc
// ***************************************************************************//
namespace bader_functions {
  bool BaderCalc(const string& bader_options,
      const string& prototype,
      const vector<string>& vspecies,
      const deque<int>& num_each_type,
      const vector<double>& vZVAL,
      const vector<double>& cutoffs,
      const vector<int>& downsample_ratios,
      const string& directory,
      ostream& oss) {
    //NO FLAGS GIVEN
    aurostd::xoption bader_flags;
    return BaderCalc(bader_flags, bader_options, prototype, vspecies, num_each_type, vZVAL, cutoffs, downsample_ratios, directory, oss);
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::FixDirectory
// ***************************************************************************//
namespace bader_functions {
  void FixDirectory(string& directory) {
    string execution_path = aurostd::getPWD(); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd");
    // translate "./" into a full path
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK DIRECTORY" << endl;
    if(directory.at(directory.length() - 1) == '/') { directory = directory.substr(0, directory.size() - 1); }
    if(directory == "./" || directory == ".") {
      directory = execution_path;
    } else if(directory.at(0) != '/') {
      directory = execution_path + "/" + directory;
    }
  }
}

// ***************************************************************************//
// bader_functions::Flags2BaderCommands
// ***************************************************************************//
namespace bader_functions {
  bool Flags2BaderCommands(aurostd::xoption& vpflow, string& bader_options, ostream& oss) {  //CO
    // Perform Bader analysis by means of the code from Henkelman Group at UT, Austin
    // Results for the Bader volume and net charge in the volume are stored in
    // the vectors passed into the function, volume and charge.

    // DEBUG
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    FixDirectory(directory);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //  Start FLAG automation
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //this is for ease of specifying
    if(vpflow.flag("BADER::AFLOWLIB_LIBRARY")) {
      //vpflow.flag("BADER::JVXL_ALL_SPECIES") has to be turned on previously, need input CUTOFF and DOWNSAMPLES
      //vpflow.flag("BADER::JVXL_CYCLIC")      has to be turned on previously, need input CUTOFF and DOWNSAMPLES
      vpflow.flag("BADER::KEEP::JVXL_ONLY",TRUE);  //only matters if JVXL_ALL_SPECIES turned on already
      //vpflow.flag("BADER::CRITICAL_POINTS",TRUE); //not ready yet
      vpflow.flag("BADER::QUIET",TRUE);
    }
    //if we JVXL_ALL_SPECIES, it should automatically CONSOLIDATE_ATOMS2SPECIES
    //if we KEEP::JVXL_ONLY, it should REMOVE_BADER_ATOMS
    if(vpflow.flag("BADER::JVXL_ALL_SPECIES")) {
      if(!vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) { vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES",TRUE); }  //set consolidate_atoms2species if not already set
      if(vpflow.flag("BADER::KEEP::JVXL_ONLY")) {
        if(!vpflow.flag("BADER::REMOVE_BADER_ATOMS")) { vpflow.flag("BADER::REMOVE_BADER_ATOMS",TRUE); }  //set consolidate_atoms4species if not already set
      }
    }
    //same here
    //we can JUST specify CONSOLIDATE_ATOMS2SPECIES,and it automatically knows to print_all=atoms
    if(vpflow.flag("BADER::CONSOLIDATE_ATOMS2SPECIES")) {
      if(!vpflow.flag("BADER::PRINT_ALL")) {
        vpflow.flag("BADER::PRINT_ALL",TRUE);  //set print_all if not already set
        vpflow.push_attached("BADER::PRINT_ALL","atom");
      } else {
        if(vpflow.getattachedscheme("BADER::PRINT_ALL") == "bader") {
          vpflow.pop_attached("BADER::PRINT_ALL");
          vpflow.push_attached("BADER::PRINT_ALL","both");
        }
      }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //  End FLAG automation
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Bader Options -- build a string of optional flags accepted by the bader code
    // discard any option not understood
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string misc_option;

    // calculate/nocalculate option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK CALCULATE | NOCALCULATE" << endl;
    if(vpflow.flag("BADER::CALCULATE")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::CALCULATE\")=" << vpflow.flag("BADER::CALCULATE") << endl;
      misc_option = vpflow.getattachedscheme("BADER::CALCULATE");
      if(misc_option == "bader") {
        oss << __AFLOW_FUNC__ << " calculating Bader atoms in molecules." << endl;
        bader_options.append("-c bader ");
      } else if(misc_option == "voronoi") {
        oss << __AFLOW_FUNC__ << " calculating population analysis based on distance." << endl;
        bader_options.append("-c voronoi ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown calculate option." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    } else if(vpflow.flag("BADER::NOCALCULATE")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::NOCALCULATE\")=" << vpflow.flag("BADER::NOCALCULATE") << endl;
      misc_option = vpflow.getattachedscheme("BADER::NOCALCULATE");
      if(misc_option == "bader") {
        oss << __AFLOW_FUNC__ << " NOT calculating Bader atoms in molecules." << endl;
        bader_options.append("-n bader ");
      } else if(misc_option == "voronoi") {
        oss << __AFLOW_FUNC__ << " NOT calculating population analysis based on distance." << endl;
        bader_options.append("-n voronoi ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown nocalculate option." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    } else {
      oss << __AFLOW_FUNC__ << " calculating default Bader atoms in molecules." << endl;
    }

    // partition option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK PARTITION" << endl;
    if(vpflow.flag("BADER::PARTITION")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::NOCALCULATE\")=" << vpflow.flag("BADER::PARTITION") << endl;
      misc_option = vpflow.getattachedscheme("BADER::PARTITION");
      if(misc_option == "neargrid") {
        oss << __AFLOW_FUNC__ << " using the default near-grid bader partitioning algorithm." << endl;
        bader_options.append("-b neargrid ");
      } else if(misc_option == "ongrid") {
        oss << __AFLOW_FUNC__ << " using the on-grid bader partitioning algorithm." << endl;
        bader_options.append("-b ongrid ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown partition option." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    } else {
      oss << __AFLOW_FUNC__ << " using the default near-grid bader partitioning algorithm." << endl;
    }

    // refine_edge_method option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK REFINE_EDGE_METHOD" << endl;
    if(vpflow.flag("BADER::REFINE_EDGE_METHOD")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::REFINE_EDGE_METHOD\")=" << vpflow.flag("BADER::REFINE_EDGE_METHOD") << endl;
      misc_option = vpflow.getattachedscheme("BADER::REFINE_EDGE_METHOD");
      if(misc_option == "-1" || misc_option == "-2" || misc_option == "-3") {
        oss << __AFLOW_FUNC__ << " requesting refinement method " << misc_option << "." << endl;
        bader_options.append("-r " + misc_option + " ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown refine_edge_method option " + misc_option + "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    } else {
      oss << __AFLOW_FUNC__ << " using default refinement method -1." << endl;
    }

    // reference option handled internally in BaderCalc function

    // vacuum option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK VACUUM" << endl;
    if(vpflow.flag("BADER::VACUUM")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::VACUUM\")=" << vpflow.flag("BADER::VACUUM") << endl;
      misc_option = vpflow.getattachedscheme("BADER::VACUUM");
      if(misc_option == "off") {
        oss << __AFLOW_FUNC__ << " NOT assigning low density points to vacuum." << endl;
        bader_options.append("-vac off ");
      } else if(misc_option == "auto") {
        oss << __AFLOW_FUNC__ << " automatically assigning density points below 1E-3 e/Ang^3 to vacuum." << endl;
        bader_options.append("-vac auto ");
      } else if(misc_option == "") {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: No density threshold given." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      } else {
        double vac_test = aurostd::string2utype<double>(misc_option);
        if(vac_test == 0) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Incorrect input for vacuum density. Must be a number (double)." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        oss << __AFLOW_FUNC__ << " assigning density points lower than " << misc_option << " to vacuum." << endl;
        bader_options.append("-vac " + misc_option + " ");
      }
    } else {
      oss << __AFLOW_FUNC__ << " automatically assigning density points below 1E-3 e/Ang^3 to vacuum." << endl;
      bader_options.append("-vac auto ");  //not actually default in bader code
    }

    // terminate option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK TERMINATE" << endl;
    if(vpflow.flag("BADER::TERMINATE")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::TERMINATE\")=" << vpflow.flag("BADER::TERMINATE") << endl;
      misc_option = vpflow.getattachedscheme("BADER::TERMINATE");
      if(misc_option == "known") {
        oss << __AFLOW_FUNC__ << " terminating trajectories when a point is surrounded by known points." << endl;
        bader_options.append("-m known ");
      } else if(misc_option == "max") {
        oss << __AFLOW_FUNC__ << " terminating trajectories when charge density maximum is reached." << endl;
        bader_options.append("-m max ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown termination option " + misc_option + "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }

    // critical points option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK CRITICAL POINTS" << endl;
    if(vpflow.flag("BADER::CRITICAL_POINTS")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.flag(\"BADER::CRITICAL_POINTS\")" << endl;
      oss << __AFLOW_FUNC__ << " finding the critical points." << endl;  // not shown anyway (only in library)
      //bader_options.append("-cp ");   //need to run as a separate command
    }

    // print_all option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK PRINT_ALL" << endl;
    if(vpflow.flag("BADER::PRINT_ALL")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::PRINT_ALL\")=" << vpflow.flag("BADER::PRINT_ALL") << endl;
      misc_option = vpflow.getattachedscheme("BADER::PRINT_ALL");
      //BvAtxxxx.dat
      if(misc_option == "atom") {
        oss << __AFLOW_FUNC__ << " combining all volumes associated with an atom and writing them to BvAtxxxx.dat." << endl;
        bader_options.append("-p all_atom ");
        //Bvolxxxx.dat
      } else if(misc_option == "bader") {
        oss << __AFLOW_FUNC__ << " writing each Bader volume to Bvolxxxx.dat." << endl;
        bader_options.append("-p all_bader ");
      } else if(misc_option == "both") {
        oss << __AFLOW_FUNC__ << " combining all volumes associated with an atom and writing them to BvAtxxxx.dat AND write each Bader volume to Bvolxxxx.dat." << endl;
        bader_options.append("-p all_atom -p all_bader ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown print_all option " + misc_option + "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }

    // print_index option
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK PRINT_INDEX" << endl;
    if(vpflow.flag("BADER::PRINT_INDEX")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::PRINT_INDEX\")=" << vpflow.flag("BADER::PRINT_INDEX") << endl;
      misc_option = vpflow.getattachedscheme("BADER::PRINT_INDEX");
      //AtIndex.dat
      if(misc_option == "atom") {
        oss << __AFLOW_FUNC__ << " printing atomic volume indices and writing them to AtIndex.dat." << endl;
        bader_options.append("-p atom_index ");
        //BvIndex.dat
      } else if(misc_option == "bader") {
        oss << __AFLOW_FUNC__ << " printing bader volume indices and writing them to BvIndex.dat." << endl;
        bader_options.append("-p bader_index ");
      } else if(misc_option == "both") {
        oss << __AFLOW_FUNC__ << " printing both atomic and bader volume indices and writing them to AtIndex.dat and BvIndex.dat, respectively." << endl;
        bader_options.append("-p atom_index -p bader_index ");
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unknown print_index option " + misc_option + "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }

    // print_select_atom/print_select_bader option
    vector<string> tokens;
    string push_command;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK PRINT_SELECT_ATOM | PRINT_SELECT_BADER" << endl;
    //BvAtxxxx.dat
    if(vpflow.flag("BADER::PRINT_SELECT_ATOM")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::PRINT_SELECT_ATOM\")=" << vpflow.flag("BADER::PRINT_SELECT_ATOM") << endl;
      if(vpflow.flag("BADER::PRINT_ALL") && vpflow.getattachedscheme("BADER::PRINT_ALL") == "atom") {
        oss << __AFLOW_FUNC__ << " ignoring print select atom command--already specified print ALL atoms command." << endl;
      } else {
        if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SELECT_ATOM"), push_command, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to gather input for print select atom command." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        oss << __AFLOW_FUNC__ << " combining all volumes associated with atom(s) " + push_command + "and writing them to BvAtxxxx.dat." << endl;  //push command comes with POST-space
        bader_options.append("-p sel_atom " + push_command);
      }
    }
    //B_wexxxx.dat
    if(vpflow.flag("BADER::PRINT_SELECT_BADER")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::PRINT_SELECT_BADER\")=" << vpflow.flag("BADER::PRINT_SELECT_BADER") << endl;
      if(vpflow.flag("BADER::PRINT_ALL") && vpflow.getattachedscheme("BADER::PRINT_ALL") == "bader") {
        oss << __AFLOW_FUNC__ << " ignoring print select bader volume command--already specified print ALL bader volumes command." << endl;
      } else {
        if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SELECT_BADER"), push_command, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to gather input for print select bader command." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        oss << __AFLOW_FUNC__ << " writing bader volume(s) " + push_command + "to Bvolxxxx.dat." << endl;  //push command comes with POST-space
        bader_options.append("-p sel_bader " + push_command);
      }
    }

    // print_sum_atom/print_sum_bader option
    //BvAt_summed.dat
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK PRINT_SUM_ATOM | PRINT_SUM_BADER" << endl;
    if(vpflow.flag("BADER::PRINT_SUM_ATOM")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::PRINT_SUM_ATOM\")=" << vpflow.flag("BADER::PRINT_SUM_ATOM") << endl;
      if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SUM_ATOM"), push_command, oss)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to gather input for print sum atom command." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      oss << __AFLOW_FUNC__ << " combining all volumes associated with atoms " + push_command + "and writing them to BvAt_summed.dat." << endl;  //push command comes with POST-space
      bader_options.append("-p sum_atom " + push_command);
      //Bvol_summed.dat
    }
    if(vpflow.flag("BADER::PRINT_SUM_BADER")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.getattachedscheme(\"BADER::PRINT_SUM_BADER\")=" << vpflow.flag("BADER::PRINT_SUM_BADER") << endl;
      if(!getPushCommand(vpflow.getattachedscheme("BADER::PRINT_SUM_BADER"), push_command, oss)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to gather input for print sum bader command." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      oss << __AFLOW_FUNC__ << " combining bader volumes " + push_command + "and writing them to Bvol_summed.dat." << endl;  //push command comes with POST-space
      bader_options.append("-p sum_bader " + push_command);
    }

    // quiet/verbose option, different from AFLOW's silent (which would silent EVERYTHING, including other functions in AFLOW)
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK QUIET" << endl;
    if(vpflow.flag("BADER::QUIET")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.flag(\"BADER::QUIET\")" << endl;
      oss << __AFLOW_FUNC__ << " silencing output." << endl;  // not shown anyway (only in library)
    } else {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vpflow.flag(\"BADER::VERBOSE\")" << endl;
      oss << __AFLOW_FUNC__ << " displaying output." << endl;
      bader_options.append("-v ");
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // END Bader Options
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return TRUE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::getPushCommand
// ***************************************************************************//
namespace bader_functions {
  bool getPushCommand(const string& misc_option, string& push_command, ostream& oss) {
    //this converts aflow lists/ranges to bader lists/ranges
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    FixDirectory(directory);

    vector<int> tokens;
    int test_value;
    push_command = "";
    //check that there is no mixing between list and range
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK THAT THERE'S NO MIXING BETWEEN LIST AND RANGE" << endl;
    if(aurostd::substring2bool(misc_option, ",") && aurostd::substring2bool(misc_option, "-")) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Do not mix notation for LIST and RANGE." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " You may only have '::' OR ':'." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    //list
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK FOR LIST" << endl;
    if(aurostd::substring2bool(misc_option, ",")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " LIST FOUND" << endl;
      aurostd::string2tokens<int>(misc_option, tokens, ",");
      for (uint i = 0; i < tokens.size(); i++) {
        if(tokens.at(i) == 0) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Input for list must be numerical and above 0: " << tokens.at(i) << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
        push_command.append(aurostd::utype2string(tokens.at(i)) + " ");
      }
      return TRUE;
    }
    //range or single item, no need to change syntax
    //do some checks
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK FOR RANGE" << endl;
    if(aurostd::substring2bool(misc_option, "-")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " RANGE FOUND" << endl;
      //check that range only has two items
      aurostd::string2tokens<int>(misc_option, tokens, "-");
      if(tokens.size() != 2) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Cannot specify more than one range." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      for (uint i = 0; i < tokens.size(); i++) {
        if(tokens.at(i) == 0) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Input for range must be numerical and above 0: " << tokens.at(i) << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
      //check that first entry !> second entry
      if(tokens.at(0) > tokens.at(1)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: First range element > second range element." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      push_command.append(aurostd::utype2string(tokens.at(0)) + "-" + aurostd::utype2string(tokens.at(1)) + " ");
      return TRUE;
    }
    //must be singular entry
    test_value = aurostd::string2utype<int>(misc_option);
    if(test_value == 0) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input for range must be numerical and above 0: " << test_value << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    push_command.append(aurostd::utype2string(test_value) + " ");
    return TRUE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::listORrange2vec
// ***************************************************************************//
namespace bader_functions {
  bool listORrange2vec(const string& misc_option, vector<int>& vout, ostream& oss) {
    int test_value;
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    FixDirectory(directory);

    vout.clear();
    //check that there is no mixing between list and range
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK THAT THERE'S NO MIXING BETWEEN LIST AND RANGE" << endl;
    if(aurostd::substring2bool(misc_option, ",") && aurostd::substring2bool(misc_option, "-")) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Do not mix notation for LIST and RANGE." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " You may only have '::' OR ':'." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    //list
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK FOR LIST" << endl;
    if(aurostd::substring2bool(misc_option, ",")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " LIST FOUND" << endl;
      aurostd::string2tokens<int>(misc_option, vout, ",");
      for (uint i = 0; i < vout.size(); i++) {
        if(vout.at(i) == 0) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Input for list must be numerical and above 0: " << vout.at(i) << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
      return TRUE;
    }
    //range
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK FOR RANGE" << endl;
    if(aurostd::substring2bool(misc_option, "-")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " RANGE FOUND" << endl;
      //check that range only has two items
      vector<uint> tokens;
      aurostd::string2tokens<uint>(misc_option, tokens, "-");
      if(tokens.size() != 2) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Cannot specify more than one range." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      for (uint i = 0; i < tokens.size(); i++) {
        if(tokens.at(i) == 0) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Input for range must be numerical and above 0: " << tokens.at(i) << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
      //check that first entry !> second entry
      if(tokens.at(0) > tokens.at(1)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: First range element > second range element." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      for (uint i = tokens.at(0); i < tokens.at(1) + 1; i++) {
        vout.push_back(i);
      }
      return TRUE;
    }
    //or single item
    test_value = aurostd::string2utype<int>(misc_option);
    if(test_value == 0) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input must be numerical and above 0: " << misc_option << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    vout.push_back(test_value);
    return TRUE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::BaderExtensionFound
// ***************************************************************************//
namespace bader_functions {
  bool BaderExtensionFound(const string& FileNameIN, string& FileNameOUT, const string& directory) {
    //input filename and directory, output full path of compressed/unccompressed file
    //we only care about .static runs (aflow), so we will only look for files with .static extension
    //or no extension
    string full_FileNameIN="";
    full_FileNameIN=directory + "/" + FileNameIN + ".static";
    if(aurostd::EFileExist(full_FileNameIN, FileNameOUT) && aurostd::FileNotEmpty(FileNameOUT)) { return TRUE; }
    if(aurostd::FileExist(full_FileNameIN) && aurostd::FileNotEmpty(full_FileNameIN)) { FileNameOUT=full_FileNameIN; return TRUE; } //CO20180524
    full_FileNameIN=directory + "/" + FileNameIN;
    if(aurostd::EFileExist(full_FileNameIN, FileNameOUT) && aurostd::FileNotEmpty(FileNameOUT)) { return TRUE; }
    if(aurostd::FileExist(full_FileNameIN) && aurostd::FileNotEmpty(full_FileNameIN)) { FileNameOUT=full_FileNameIN; return TRUE; } //CO20180524
    return FALSE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::BaderExtensionFound
// ***************************************************************************//
namespace bader_functions {
  bool BaderExtensionFound(const string& FileNameIN, const string& directory) {
    string FileNameOUT;
    return BaderExtensionFound(FileNameIN, FileNameOUT, directory);
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::adjust_header
// ***************************************************************************//
namespace bader_functions {
  void adjust_header(string& new_header, stringstream& FileIN_ss) {
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //adds new_header at top of stringstream
    string content;
    vector<string> vcontent;
    content = FileIN_ss.str();
    //cerr << "LOOOOK " << endl;
    //cerr << content << endl;
    //cerr << new_header << endl;
    aurostd::string2vectorstring(content, vcontent);
    //get right number of spaces for formatting
    uint white_space_counter1 = 0, white_space_counter2 = 0;
    while (isspace(vcontent.at(0).at(white_space_counter1))) { white_space_counter1++; }
    while (isspace(new_header.at(white_space_counter2))) { white_space_counter2++; }
    if(white_space_counter1 != white_space_counter2) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " ADDING " << white_space_counter1 << " TO HEADER" << endl;
      vcontent.at(0) = new_header.insert(0, white_space_counter1, ' ');
    } else {
      vcontent.at(0) = new_header;
    }
    //rebuild
    FileIN_ss.str("");
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " WRITING OUT STRINGSTREAM WITH ADJUSTED HEADER" << endl;
    for (uint i = 0; i < vcontent.size(); i++) { FileIN_ss << vcontent.at(i) << endl; }
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
// ***************************************************************************//
namespace bader_functions {
  string prepare_CHGCAR_4_Jmol(aurostd::xoption vpflow) {
    ostringstream oss;
    //adjust header of CHGCAR for reading in Jmol
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    string usage_usage = "aflow --prep4jmol [options]";
    string usage_options = aurostd::liststring2string("options = --usage",
        "          --outcar=OUTCAR1,OUTCAR2,...",
        "          --zip");
    // output usage
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHECK USAGE" << endl;
    if(vpflow.flag("PREPARE_CHGCAR_4_JMOL::USAGE")) {
      init::ErrorOption( vpflow.getattachedscheme("PREPARE_CHGCAR_4_JMOL"), "bader_functions::prepare_CHGCAR_4_Jmol()", aurostd::liststring2string(usage_usage, usage_options));
      return oss.str();
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " GATHER CHGCAR_FILES" << endl;
    vector<string> chgcar_files;
    string directory = aurostd::getPWD(); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd");
    string input = vpflow.getattachedscheme("PREPARE_CHGCAR_4_JMOL");
    if(aurostd::substring2bool(input, ",")) {
      aurostd::string2tokens(input, chgcar_files, ",");
    } else {
      chgcar_files.push_back(input);
    }

    string outcar_file;
    if(vpflow.flag("PREPARE_CHGCAR_4_JMOL::OUTCAR")) {
      outcar_file = vpflow.getattachedscheme("PREPARE_CHGCAR_4_JMOL::OUTCAR");
      if(outcar_file.at(0) != '/') {  //not root dir, so we can append directory
        outcar_file = directory + "/" + outcar_file;
      }
      oss << __AFLOW_FUNC__ << " OUTCAR file specied: " << outcar_file << "." << endl;
    }

    bool zip_file = false;
    if(vpflow.flag("PREPARE_CHGCAR_4_JMOL::ZIP")) { zip_file = true; }

    string species_header, path, file, chgcar_file;
    vector<string> _path;
    //TO RUN QUICKLY, SPECIFY OUTCAR
    for (uint i = 0; i < chgcar_files.size(); i++) {
      path = "";
      oss << __AFLOW_FUNC__ << " working on " << chgcar_files.at(i) << "." << endl;
      //remove file from path so we can search path
      chgcar_file = chgcar_files.at(i);
      if(chgcar_file.at(0) != '/') {  //not root dir, so we can append directory
        chgcar_file = directory + "/" + chgcar_file;
      }
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CHGCAR " << chgcar_file << endl;
      //only look every time if outcar_file not specified, otherwise, only look once
      if(outcar_file.empty() || species_header.empty()) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " LOOKING FOR OUTCAR" << endl;
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " GETTING PATH FROM CHGCAR" << endl;
        if(aurostd::substring2bool(chgcar_file, "/")) {
          aurostd::string2tokens(chgcar_file, _path, "/");
          if(chgcar_file.at(0) == '/') {  //if root, which it will be
            path.append("/");
          }
          for (uint i = 0; i < _path.size() - 1; i++) {
            path.append(_path.at(i) + "/");  //last character will be '/', chgcar_file==_path.back();
          }
          path = path.substr(0, path.size() - 1);  //pop off last '/'
          file = _path.back();
        } else {
          path = ".";
        }  //chgcar_file==chgcar_file;}
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " PATH " << path << endl;
      if(!get_species_string(outcar_file, species_header, path, file, oss)) {  //look in same path as CHGCAR
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to gather new header for file." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      if(!vpflow.flag("PREPARE_CHGCAR_4_JMOL::OUTCAR")) {
        outcar_file = "";  //reset so we can look again later (once per CHGCAR)
      }
    }
    if(!prepare_CHGCAR_4_Jmol(chgcar_file, species_header, zip_file, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to rewrite file with appropriate header." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    oss << __AFLOW_FUNC__ << " finished rewriting " << chgcar_file << "." << endl;
  }
  return oss.str();
}
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::get_species_string
// ***************************************************************************//
namespace bader_functions {
  bool get_species_string(string& outcar_file, string& species_string, const string& dir_to_look, const string& file, ostream& oss) {
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    FixDirectory(directory);

    species_string = "";  //clear
    string here = ".";
    //CHECK IF OUTCAR CANNOT BE FOUND
    if(!outcar_file.empty()) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " OUTCAR SPECIFIED, NO NEED TO LOOK" << endl;
      if((!aurostd::FileExist(outcar_file)) || aurostd::FileEmpty(outcar_file)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to locate specified OUTCAR (or its empty): " << outcar_file << "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      //FIRST CHECK PATH, OTHERWISE CHECK CURRENT DIRECTORY
      //C++ SHOULD SHORT CIRCUIT HERE to yield correct path
    } else {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " OUTCAR NOT SPECIFIED, LOOKING FOR IT" << endl;
      //first check in CHGCAR file path
      //otherwise check in current path
      if(!BaderExtensionFound("OUTCAR", outcar_file, dir_to_look) && !BaderExtensionFound("OUTCAR", outcar_file, ".")) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to locate OUTCAR/OUTCAR.static (compressed or otherwise) in path(" << file << ") or pwd." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
    }
    oss << __AFLOW_FUNC__ << " using " << outcar_file << " to gather species." << endl;
    xOUTCAR outcar(outcar_file);
    for (uint i = 0; i < outcar.species.size(); i++) {
      species_string.append(outcar.species.at(i));
      if(i < outcar.species.size() - 1) { species_string.append(" "); }
    }
    //species_string=species_string.substr(0,species_string.size()-1);  //FIX ME
    oss << __AFLOW_FUNC__ << " new header: " << species_string << "." << endl;
    return TRUE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
// ***************************************************************************//
namespace bader_functions {
  bool prepare_CHGCAR_4_Jmol(const string& _chgcar_file, string& species_header, bool zip_file, ostream& oss) {
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    FixDirectory(directory);

    //overwrite file so header has species
    string chgcar_file = aurostd::CleanFileName(_chgcar_file), chgcar_file_compressed, chgcar_file_uncompressed;
    bool any_file_exists=false;
    any_file_exists=(any_file_exists || (aurostd::EFileExist(chgcar_file,chgcar_file_compressed) && aurostd::FileNotEmpty(chgcar_file_compressed)));
    any_file_exists=(any_file_exists || (aurostd::FileExist(chgcar_file) && aurostd::FileNotEmpty(chgcar_file)));
    if(!any_file_exists){
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to find " << chgcar_file << " (compressed or otherwise, or its empty)." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << __AFLOW_FUNC__ << " editing header of " << chgcar_file << " for Jmol visualization." << endl;
    aurostd::IsCompressed(chgcar_file, chgcar_file_uncompressed);  //get uncompressed variant
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " UNCOMPRESSED VARIANT OF CHGCAR_FILE IS " << chgcar_file_uncompressed << endl;
    stringstream chgcar_ss;
    //get stringstream
    if(!chgcar_file_compressed.empty()){aurostd::efile2stringstream(chgcar_file, chgcar_ss);}  //already checked if EFileExist()
    else {aurostd::file2stringstream(chgcar_file, chgcar_ss);}
    adjust_header(species_header, chgcar_ss);        //fix it
    aurostd::stringstream2file(chgcar_ss, chgcar_file_uncompressed);   //write back out
    if(zip_file && !aurostd::MatchCompressed(chgcar_file, chgcar_file_uncompressed)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to zip " << chgcar_file_uncompressed << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << __AFLOW_FUNC__ << " finish editing " << chgcar_file_uncompressed << "." << endl;
    oss << __AFLOW_FUNC__ << chgcar_file_uncompressed << " should be ready for Jmol." << endl;
    return TRUE;
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
// ***************************************************************************//
namespace bader_functions {
  bool prepare_CHGCAR_4_Jmol(string& _chgcar_file, string& species_header, ostream& oss) {
    bool zip_file = false;
    return prepare_CHGCAR_4_Jmol(_chgcar_file, species_header, zip_file, oss);
  }
}  // namespace bader_functions

// ***************************************************************************//
// bader_functions::prepare_CHGCAR_4_Jmol
// ***************************************************************************//
namespace bader_functions {
  bool prepare_CHGCAR_4_Jmol(string& _chgcar_file, string& species_header) {
    ostringstream oss;
    return prepare_CHGCAR_4_Jmol(_chgcar_file, species_header, oss);
  }
}  // namespace bader_functions

//[CO20180220 - moved to aurostd]***************************************************************************//
//[CO20180220 - moved to aurostd] bader_functions::efile2tempfile
//[CO20180220 - moved to aurostd]***************************************************************************//
//[CO20180220 - moved to aurostd]namespace bader_functions {
//[CO20180220 - moved to aurostd]}  // namespace bader_functions

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  string CHGCAR2JVXL(aurostd::xoption& vpflow) {  //CO
    ostringstream oss;
    // handles flags for CHGCAR2JVXL
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    bader_functions::FixDirectory(directory);

    vector<string> tokens;
    string input = vpflow.getattachedscheme("CHGCAR2JVXL");

    // case 1: cyclic variable
    // 2a:  FILE1,FILE2,FILE3...::CUTOFF[::DOWNSAMPLE]
    // 2b:  FILE::CUTOFF1,CUTOFF2…[::DOWNSAMPLE]
    // 2c:  FILE::CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...]
    // 2d:  some combination of above
    vector<string> chgcar_files, sets;
    vector<double> cutoffs;
    vector<int> downsample_ratios;
    bool cyclic = FALSE;
    if(aurostd::substring2bool(input, "::")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CYCLIC ROUTINE DETECTED" << endl;
      cyclic = TRUE;
      aurostd::string2tokens(input, tokens, "::");
      if(tokens.size() < 2 || tokens.size() > 3) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input - number of tokens (::)." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Needs to be either: " << endl;
        oss << "                             " << "FILE1,FILE2,FILE3...::CUTOFF[::DOWNSAMPLE], " << endl;
        oss << "                             " << "FILE::CUTOFF1,CUTOFF2…[::DOWNSAMPLE], or " << endl;
        oss << "                             " << "FILE::CUTOFF[::DOWNSAMPLE1,DOWNSAMPLE2,...]." << endl;
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      //make sure there's no mixing
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CYCLIC ROUTINE DETECTED" << endl;
      for (uint i = 0; i < tokens.size(); i++) {
        if(aurostd::substring2bool(tokens.at(i), ":")) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input, cannot specify sets (:) and cyclic (::) parameters." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return oss.str();
        }
      }
      if(aurostd::substring2bool(tokens.at(0), ",")) {
        aurostd::string2tokens(tokens.at(0), chgcar_files, ",");
      } else {
        chgcar_files.push_back(tokens.at(0));
      }
      if(aurostd::substring2bool(tokens.at(1), ",")) {
        aurostd::string2tokens<double>(tokens.at(1), cutoffs, ",");
      } else {
        cutoffs.push_back(aurostd::string2utype<double>(tokens.at(1)));
      }
      if(tokens.size() == 3) {
        if(aurostd::substring2bool(tokens.at(2), ",")) {
          aurostd::string2tokens<int>(tokens.at(2), downsample_ratios, ",");
        } else {
          downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(2)));
        }
      }
      if(downsample_ratios.empty()) {
        if(chgcar_files.size() == 1 && chgcar_files.size() == cutoffs.size()) {
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " CYCLIC ROUTINE IGNORED, ONLY ONE SET" << endl;
          cyclic = FALSE;
        }
      } else {
        if(chgcar_files.size() == 1 && chgcar_files.size() == cutoffs.size() && chgcar_files.size() == downsample_ratios.size()) {
          if(LDEBUG) cerr << __AFLOW_FUNC__ << " CYCLIC ROUTINE IGNORED, ONLY ONE SET" << endl;
          cyclic = FALSE;
        }
      }
      //case 2: sets, FILE1,CUTOFF1[,DOWNSAMPLE1]:FILE2,CUTOFF2[,DOWNSAMPLE2]:...
    } else if(aurostd::substring2bool(input, ":")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " SETS ROUTINE DETECTED" << endl;
      aurostd::string2tokens(input, sets, ":");
      for (uint i = 0; i < sets.size(); i++) {
        aurostd::string2tokens(sets.at(i), tokens, ",");
        if(tokens.size() < 2 || tokens.size() > 3) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input " << i + 1 << "." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Needs to be FILE1,CUTOFF1[,DOWNSAMPLE1]:FILE2,CUTOFF2[,DOWNSAMPLE2]:..." << endl;
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return oss.str();
        }
        chgcar_files.push_back(tokens.at(0));
        cutoffs.push_back(aurostd::string2utype<double>(tokens.at(1)));
        if(tokens.size() == 3) { downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(2))); }
      }
      // case 3: single set
    } else if(aurostd::substring2bool(input, ",")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " SINGLE SET ROUTINE DETECTED" << endl;
      aurostd::string2tokens(input, tokens, ",");
      if(tokens.size() < 2 || tokens.size() > 3) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Incorrect format for input - number of tokens (,)." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Needs to be: FILE,CUTOFF[,DOWNSAMPLE]." << endl;
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return oss.str();
      }
      chgcar_files.push_back(tokens.at(0));
      cutoffs.push_back(aurostd::string2utype<double>(tokens.at(1)));
      if(tokens.size() == 3) { downsample_ratios.push_back(aurostd::string2utype<int>(tokens.at(2))); }
    } else {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Incorrect input, need at least FILE1,CUTOFF1." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return oss.str();
    }

    if(vpflow.flag("CHGCAR2JVXL::OUTPUT")) {
      if(!cyclic) {
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " OUTPUT FILES SPECIFIED" << endl;
        //they must be sets, deal with this case first
        vector<string> output_files;
        string output = vpflow.getattachedscheme("CHGCAR2JVXL::OUTPUT");
        if(aurostd::substring2bool(output, "::")) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Incorrect input delimiter for output_files." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Please remove any \"::\"." << endl;
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return oss.str();
        }
        if(aurostd::substring2bool(output, ",")) {
          aurostd::string2tokens(output, output_files, ",");
        } else {
          output_files.push_back(output);
        }
        if(!CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to convert CHGCAR files to .jvxl's." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return oss.str();
        }
        return oss.str();
      } else {
        oss << endl;
        oss << __AFLOW_FUNC__ << " WARNING: An output file name can only be specified for sets." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " You have specified CYCLIC mode (::)." << endl;
        oss << __AFLOW_FUNC__ << " Ignoring output file name input." << endl;
        oss << endl;
      }
    }
    if(!CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, cyclic, oss)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to convert CHGCAR files to .jvxl's." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return oss.str();
    }
    return oss.str();
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios, const bool& cyclic, ostream& oss) {  //CO
    //NO OUTPUT_FILES SPECIFIED
    vector<string> output_files;
    return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, cyclic, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files,
      const vector<double>& cutoffs,
      const vector<int>& downsample_ratios,
      vector<string>& output_files,
      const bool& cyclic,
      ostream& oss) {                 //CO
    //debug
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    bader_functions::FixDirectory(directory);

    //NO OUTPUT_FILES SPECIFIED, so have to specify whether to cycle through variables or not
    if(cyclic) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " CYCLIC ROUTINE DETECTED" << endl;

      if(downsample_ratios.size() != 0) {
        for (uint i = 0; i < chgcar_files.size(); i++) {
          for (uint j = 0; j < cutoffs.size(); j++) {
            for (uint k = 0; k < downsample_ratios.size(); k++) {
              if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(j), downsample_ratios.at(k), oss)) {
                oss << endl;
                oss << __AFLOW_FUNC__ << " ERROR: Unable to convert CHGCAR files to .jvxl's." << " "; //<< endl //CO20180502;
                oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
                oss << __AFLOW_FUNC__ << " Exiting." << endl;
                oss << endl;
                return FALSE;
              }
            }
          }
        }
      } else {
        for (uint i = 0; i < chgcar_files.size(); i++) {
          for (uint j = 0; j < cutoffs.size(); j++) {
            if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(j), oss)) {
              oss << endl;
              oss << __AFLOW_FUNC__ << " ERROR: Unable to convert CHGCAR files to .jvxl's." << " "; //<< endl //CO20180502;
              oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
              oss << __AFLOW_FUNC__ << " Exiting." << endl;
              oss << endl;
              return FALSE;
            }
          }
        }
      }
    } else {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " SETS ROUTINE DETECTED" << endl;

      //first check if we have the right number of inputs
      if(downsample_ratios.size() != 0) {
        if(chgcar_files.size() != cutoffs.size() || chgcar_files.size() != downsample_ratios.size()) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unequal set parameters found." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Must have equal number of CHGCAR files, cutoffs, and downsample ratios." << endl;
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      } else {
        if(chgcar_files.size() != cutoffs.size()) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unequal set parameters found." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Must have equal number of CHGCAR files and cutoffs." << endl;
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }

      //figure out if we need output_filenames, otherwise check if we have the right size
      if(output_files.empty()) {
        for (uint i = 0; i < chgcar_files.size(); i++) {
          if(downsample_ratios.size() != 0) {
            output_files.push_back(CHGCAR2JVXL_get_output_filename(chgcar_files.at(i), cutoffs.at(i), downsample_ratios.at(i)));
          } else {
            output_files.push_back(CHGCAR2JVXL_get_output_filename(chgcar_files.at(i), cutoffs.at(i)));
          }
        }
      } else {
        if(downsample_ratios.size() != 0) {
          if(chgcar_files.size() != cutoffs.size() || chgcar_files.size() != downsample_ratios.size() || chgcar_files.size() != output_files.size()) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: Unequal set parameters found." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Number of OUTPUT files must match number of CHGCAR files, cutoffs, and downsample ratios." << endl;
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
            return FALSE;
          }
        } else {
          //downsample_ratios.size()==0
          if(chgcar_files.size() != cutoffs.size() || chgcar_files.size() != output_files.size()) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: Unequal set parameters found." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Number of OUTPUT files must match number of CHGCAR files and cutoffs." << endl;
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
            return FALSE;
          }
        }
      }

      //call function
      if(downsample_ratios.size() != 0) {
        for (uint i = 0; i < chgcar_files.size(); i++) {
          if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(i), downsample_ratios.at(i), output_files.at(i), oss)) {
            oss << endl;
            oss << __AFLOW_FUNC__ << " ERROR: Unable to convert CHGCAR files to .jvxl's." << " "; //<< endl //CO20180502;
            oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
            oss << __AFLOW_FUNC__ << " Exiting." << endl;
            oss << endl;
            return FALSE;
          }
        }
        return TRUE;
      }
      //downsample_ratios.size()==0
      for (uint i = 0; i < chgcar_files.size(); i++) {
        if(!CHGCAR2JVXL(chgcar_files.at(i), cutoffs.at(i), output_files.at(i), oss)) {
          oss << endl;
          oss << __AFLOW_FUNC__ << " ERROR: Unable to convert CHGCAR files to .jvxl's." << " "; //<< endl //CO20180502;
          oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
          oss << __AFLOW_FUNC__ << " Exiting." << endl;
          oss << endl;
          return FALSE;
        }
      }
      return TRUE;
    }
    return TRUE;
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const vector<int>& downsample_ratios, ostream& oss) {  //CO
    //NO CYCLIC VARIABLE, so assume sets
    bool cyclic = FALSE;
    return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, cyclic, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, const bool& cyclic, ostream& oss) {  //CO
    //NO DOWNSAMPLE
    vector<int> downsample_ratios;
    return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, cyclic, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files, const vector<double>& cutoffs, ostream& oss) {  //CO
    //NO DOWNSAMPLE OR CYCLIC
    bool cyclic = FALSE;
    return CHGCAR2JVXL(chgcar_files, cutoffs, cyclic, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files,
      const vector<double>& cutoffs,
      const vector<int>& downsample_ratios,
      vector<string>& output_files,
      ostringstream& oss) {  //CO
    //OUTPUT_FILES SPECIFIED, so must be sets
    bool cyclic = FALSE;
    return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, cyclic, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(vector<string>& chgcar_files,
      const vector<double>& cutoffs,
      vector<string>& output_files,
      ostringstream& oss) {  //CO
    //NO DOWNSAMPLE SPECIFIED
    vector<int> downsample_ratios;
    return CHGCAR2JVXL(chgcar_files, cutoffs, downsample_ratios, output_files, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(string chgcar_file, const double cutoff, const int downsample_ratio, string output_file, ostream& oss) {
    //ACTUALLY HANDLES JMOL INTERACTION
    bool LDEBUG = (FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;

    //CO20180220 - directory stuff for logging
    string directory = ".";  //default
    if(XHOST.vflag_control.flag("DIRECTORY")) {
      directory = XHOST.vflag_control.getattachedscheme("DIRECTORY");  //XHOST.vflag_control.flag("DIRECTORY");//aflags.Directory;//vpflow.getattachedscheme("BADER::DIRECTORY");
    }
    bader_functions::FixDirectory(directory);

    //string command="java -jar Jmol.jar -ionx";
    string command;
    //Look for JmolData.jar first
    if((XHOST.hostname == "nietzsche.mems.duke.edu" || XHOST.hostname == "aflowlib.duke.edu") && aurostd::FileExist("/usr/local/bin/JmolData.jar")) {
      command = "java -jar /usr/local/bin/JmolData.jar";
      //command="/usr/local/bin/jmol";
    } else if(aurostd::IsCommandAvailable("JmolData.jar")) {
      command = "java -jar " + XHOST.command("JmolData.jar");
      //Look for jmol instead
    } else if(aurostd::IsCommandAvailable("jmol")) {
      command = XHOST.command("jmol");
    } else if(aurostd::IsCommandAvailable("jmol.sh")) {
      command = XHOST.command("jmol.sh");
    } else if(aurostd::IsCommandAvailable("Jmol.jar")) {
      command = "java -jar " + XHOST.command("Jmol.jar");
    } else {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Unable to locate either JmolData (preferred), jmol, jmol.sh, or Jmol.jar in path." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Please install Jmol." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    oss << __AFLOW_FUNC__ << " using " + command + " to run Jmol." << endl;
    command.append(" -ionx");

    //check if files exist/no negatives are present, throw error
    if(chgcar_file.empty() || (!aurostd::FileExist(chgcar_file)) || aurostd::FileEmpty(chgcar_file)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input CHGCAR file does not exist (or its empty): " << chgcar_file << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(cutoff < 0) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Cutoff cannot be negative: " << cutoff << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(cutoff == 0) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input for cutoff be numerical and above 0: " << cutoff << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(downsample_ratio < 0) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Downsample ratio cannot be negative: " << downsample_ratio << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(downsample_ratio == 0) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Input for downsample ratio be numerical and above 0: " << downsample_ratio << "." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(aurostd::FileExist(output_file)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " WARNING: Output file " << output_file << " EXISTS." << endl;
      oss << __AFLOW_FUNC__ << " WARNING: Overwriting output file " << output_file << "." << endl;
      oss << endl;
    }
    if(output_file.substr(output_file.length() - 5, 5) != ".jvxl") {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: OUTPUT file must have .jvxl extension." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
      return FALSE;
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " INPUTS SEEM TO BE ENTERED CORRECTLY" << endl;

    bool clean_up = false;
    //check if INPUT is compressed, and decompress
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " FIND / DECOMPRESS INPUT" << endl;
    bool tempfile_created=false;
    if(aurostd::IsCompressed(chgcar_file)) {
      clean_up = true;
      oss << __AFLOW_FUNC__ << " decompressing " << chgcar_file << "." << endl;
      if(!aurostd::efile2tempfile(chgcar_file, chgcar_file, tempfile_created)) {
        oss << endl;
        oss << __AFLOW_FUNC__ << " ERROR: Unable to decompress " << chgcar_file << "." << " "; //<< endl //CO20180502;
        oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
        oss << __AFLOW_FUNC__ << " Exiting." << endl;
        oss << endl;
        return FALSE;
      }
      oss << __AFLOW_FUNC__ << " created temporary file " << chgcar_file << "." << endl;
    }
    //NO LONGER NEEDED IF WORKING IN NO DISPLAY VERSION
    //if(!jmoldata) {
    //  command.append(" \""+chgcar_file+"\"");
    //}
    command.append(" -j 'isosurface");
    if(downsample_ratio != 1) {  //already took care of 0
      command.append(" downsample " + aurostd::utype2string(downsample_ratio));
    }
    //command.append(" cutoff "+aurostd::utype2string(cutoff)+" \""+chgcar_file+"\"; write ISOSURFACE");
    stringstream num_prec;
    num_prec << cutoff;
    command.append(" cutoff " + num_prec.str() + " \"" + chgcar_file + "\"; write");
    num_prec.str("");
    command.append(" \"" + output_file + "\";'");
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Jmol command: " << command << endl;
    string jmol_output;
    jmol_output = aurostd::execute2string(command);
    //jmol_output=aurostd::RemoveWhiteSpaces(aurostd::execute2string(command)); //no longer needed
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " Jmol output:  " << endl;
    if(LDEBUG) cerr << jmol_output;
    //I WOULD LIKE TO DO THIS, BUT JMOL'S OUTPUT MAY CHANGE WITH DIFFERENT VERSIONS, SO JUST CHECK FOR JVXL FILE
    //check output and that jvxl file exists
    //if(!(aurostd::substring2bool(jmol_output,"isosurface1created")&&aurostd::substring2bool(jmol_output,"OKXJVXL")&&
    //  aurostd::FileExist(output_file))) {
    //  oss << endl;
    //  oss << __AFLOW_FUNC__ << " ERROR: Jmol unable to create .jvxl file." << " "; //<< endl //CO20180502;
    //  oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
    //  oss << __AFLOW_FUNC__ << " Please correctly set your display variable and check your version of JmolData/Jmol." << endl;
    //  oss << __AFLOW_FUNC__ << " Exiting." << endl;
    //  oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
    //  if(tempfile_created){aurostd::RemoveFile(chgcar_file);}
#endif
    //  return FALSE;
    //}
    //isosurface1 created
    //NO LONGER APPLICABLE
    //if(jmol_output!="-5") {    //expected command from Jmol, if it doesn't come out, then likely display variable off
    //  oss << endl;
    //  oss << __AFLOW_FUNC__ << " ERROR: Jmol unable to create .jvxl file." << " "; //<< endl //CO20180502;
    //  oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
    //  oss << __AFLOW_FUNC__ << " Please correctly set your display variable." << endl;
    //  oss << __AFLOW_FUNC__ << " Exiting." << endl;
    //  oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
    //  if(tempfile_created){aurostd::RemoveFile(chgcar_file);}
#endif
    //  return FALSE;
    //}
    if((!aurostd::FileExist(output_file)) || aurostd::FileEmpty(output_file)) {
      oss << endl;
      oss << __AFLOW_FUNC__ << " ERROR: Jmol unable to create .jvxl file." << " "; //<< endl //CO20180502;
      oss << __AFLOW_FUNC__ << " [dir=" << directory << "]" << endl; //CO20180220
      oss << __AFLOW_FUNC__ << " Please correctly set your display variable and check your version of JmolData/Jmol." << endl;
      oss << __AFLOW_FUNC__ << " Exiting." << endl;
      oss << endl;
#ifndef _AFLOW_TEMP_PRESERVE_
      if(tempfile_created){aurostd::RemoveFile(chgcar_file);}
#endif
      return FALSE;
    }
    oss << __AFLOW_FUNC__ << output_file << " successfully created." << endl;
    //REMOVE TEMP CHGCAR_FILE
#ifndef _AFLOW_TEMP_PRESERVE_
    if(clean_up) {
      oss << __AFLOW_FUNC__ << " removing temporary file " << chgcar_file << "." << endl;
      if(tempfile_created){aurostd::RemoveFile(chgcar_file);}
    }
#endif
    return TRUE;
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, const int& downsample_ratio, ostream& oss) {
    //MISSING OUTPUT
    string output_file = CHGCAR2JVXL_get_output_filename(chgcar_file, cutoff, downsample_ratio);
    return CHGCAR2JVXL(chgcar_file, cutoff, downsample_ratio, output_file, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, string& output_file, ostream& oss) {
    //MISSING DOWNSAMPLE
    int downsample_ratio = 1;
    return CHGCAR2JVXL(chgcar_file, cutoff, downsample_ratio, output_file, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL
// ***************************************************************************//
namespace pflow {
  bool CHGCAR2JVXL(string chgcar_file, const double& cutoff, ostream& oss) {
    //MISSING DOWNSAMPLE AND OUTPUT
    int downsample_ratio = 1;
    return CHGCAR2JVXL(chgcar_file, cutoff, downsample_ratio, oss);
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL_get_output_filename
// ***************************************************************************//
namespace pflow {
  string CHGCAR2JVXL_get_output_filename(string chgcar_file, const double& cutoff, const int& downsample_ratio) {
    //for handling compressed files, e.g. we don't want .EXT to be part of output name
    string FileNameIN_decompressed;
    aurostd::IsCompressed(chgcar_file, FileNameIN_decompressed);
    stringstream num_prec;
    num_prec << std::fixed << setprecision(2) << cutoff;
    string cutoff_string = num_prec.str();  //avoid issues with setprecision(20) in string2utype<double>
    //string cutoff_string=aurostd::utype2string(cutoff);
    if(aurostd::substring2bool(cutoff_string, "0.")) { cutoff_string = aurostd::StringSubst(cutoff_string, "0.", ""); }
    if(aurostd::substring2bool(cutoff_string, ".")) { cutoff_string = aurostd::StringSubst(cutoff_string, ".", ""); }
    string output_file = FileNameIN_decompressed + "_Bader_" + cutoff_string;
    if(downsample_ratio != 0) {
      output_file.append("_" + aurostd::utype2string(downsample_ratio));
    }
    output_file.append(".jvxl");
    return output_file;
  }
}  // namespace pflow

// ***************************************************************************//
// pflow::CHGCAR2JVXL_get_output_filename
// ***************************************************************************//
namespace pflow {
  string CHGCAR2JVXL_get_output_filename(string chgcar_file, const double& cutoff) {
    //NO DOWNSAMPLE
    int downsample_ratio = 0;
    return CHGCAR2JVXL_get_output_filename(chgcar_file, cutoff, downsample_ratio);
  }
}  // namespace pflow

// ***************************************************************************
// *                                                                         *
// *              AFlow KEVIN RASCH - Duke University 2013                   *
// *                                                                         *
// ***************************************************************************
