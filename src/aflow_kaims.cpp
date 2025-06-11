// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to run AIMS in KBIN
// Stefano Curtarolo - 2007-2019 Duke
// Corey Oses - 2015-2019 Duke

#ifndef _AFLOW_KBIN_CPP
#define _AFLOW_KBIN_CPP

#include "aflow.h"

using aurostd::RemoveWhiteSpaces;
using aurostd::RemoveWhiteSpacesFromTheBack;

#define _STROPT_ string("[AIMS_FORCE_OPTION]")

namespace KBIN{
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(string& AflowIn,_aflags& aflags,_kflags& kflags) {
    ofstream FileMESSAGE("/dev/null");
    return KBIN::AIMS_Get_AIMSflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags);
  }

  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(string& AflowIn,ofstream& FileMESSAGE,_aflags& aflags,_kflags& kflags) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    _aimsflags aimsflags;
    AflowIn=aurostd::RemoveComments(AflowIn); // for safety
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
    string BflowIn=AflowIn;

    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (START)" << endl;

    aimsflags.KBIN_AIMS_RUN.clear();
    if((aurostd::substring2bool(AflowIn,"[AIMS_RUN_GENERATE]") || aurostd::substring2bool(AflowIn,"[AIMS_RUN]GENERATE")) || aflags.KBIN_GEN_AIMS_FROM_AFLOWIN) 
      aimsflags.KBIN_AIMS_RUN.push("GENERATE");

    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_FILE]")) aimsflags.KBIN_AIMS_CONTROL_FILE.push("KEYWORD");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_FILE]SYSTEM_AUTO",TRUE)) aimsflags.KBIN_AIMS_CONTROL_FILE.push("SYSTEM_AUTO");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_FILE]FILE=",TRUE)) aimsflags.KBIN_AIMS_CONTROL_FILE.push("FILE");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_FILE]COMMAND=",TRUE)) aimsflags.KBIN_AIMS_CONTROL_FILE.push("COMMAND");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_EXPLICIT]")) aimsflags.KBIN_AIMS_CONTROL_MODE.push("EXPLICIT");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_EXPLICIT]STOP")) aimsflags.KBIN_AIMS_CONTROL_MODE.push("EXPLICIT_START_STOP");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_IMPLICIT]")) aimsflags.KBIN_AIMS_CONTROL_MODE.push("IMPLICIT");
    if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_EXTERNAL]")) aimsflags.KBIN_AIMS_CONTROL_MODE.push("EXTERNAL");

    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]"))  aimsflags.KBIN_AIMS_GEOM_FILE.push("KEYWORD");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]")) aimsflags.KBIN_AIMS_GEOM_MODE.push("EXPLICIT");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]START") &&  aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]STOP")) aimsflags.KBIN_AIMS_GEOM_MODE.push("EXPLICIT_START_STOP");
    if((aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]START.") && aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]STOP."))) aimsflags.KBIN_AIMS_GEOM_MODE.push("EXPLICIT_START_STOP_POINT");

    if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT") && !kflags.KBIN_FROZSL) {  // NO FROZSL
      if(LDEBUG) cerr << "DEBUG: aimsflags.KBIN_AIMS_GEOM_MODE.flag(\"EXPLICIT_START_STOP_POINT\")=" << aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT") << endl;
      if(LDEBUG) cerr << "DEBUG: kflags.KBIN_PHONONS_CALCULATION_FROZSL=" << kflags.KBIN_PHONONS_CALCULATION_FROZSL << endl;
      if(LDEBUG) cerr << "DEBUG: kflags.KBIN_FROZSL_DOWNLOAD=" << kflags.KBIN_FROZSL_DOWNLOAD << endl;
      if(LDEBUG) cerr << "DEBUG: kflags.KBIN_FROZSL_FILE=" << kflags.KBIN_FROZSL_FILE << endl;
      stringstream input_file;
      // loading
      input_file.clear();
      if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
        FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags);
        FROZSL::Extract_INPUT(AflowIn,FileMESSAGE,input_file,aflags,kflags);
      }
      if(kflags.KBIN_FROZSL_DOWNLOAD)    FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags);
      if(kflags.KBIN_FROZSL_FILE)        FROZSL::File_INPUT(AflowIn,FileMESSAGE,input_file,aflags,kflags);

      if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) input_file.str(AflowIn);

      aimsflags.KBIN_AIMS_GEOM_MODE.push("EXPLICIT_START_STOP_POINT");
      // done loading now load structures up
      aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP",FALSE); // some default
      aurostd::substring2strings(input_file.str(),aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING,"[AIMS_GEOM_MODE_EXPLICIT]START.");
      // some verbose
      for(uint i=0;i<aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size();i++)
        if(LDEBUG) cerr << "DEBUG= " << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.at(i) << endl;
      // load up the structures
      for(uint i=0;i<aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size();i++) {
        string START="[AIMS_GEOM_MODE_EXPLICIT]START";
        string STOP="[AIMS_GEOM_MODE_EXPLICIT]STOP";
        START="[AIMS_GEOM_MODE_EXPLICIT]START."+aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.at(i);
        STOP="[AIMS_GEOM_MODE_EXPLICIT]STOP."+aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.at(i);
        stringstream GEOM;
        GEOM.str(aurostd::substring2string(input_file.str(),START,STOP,-1));
        //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(input_file.str(),START) && aurostd::substring2bool(input_file.str(),STOP))
          //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(input_file.str(),GEOM,START,STOP);
        if(!GEOM.str().empty()) aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.push_back(xstructure(GEOM,IOAIMS_AUTO));
      }
      if(LDEBUG) cerr << "DEBUG " << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size() << endl;
      if(LDEBUG) cerr << "DEBUG " << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.size() << endl;
      if(aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size() != aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.size()) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size() != aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.size()",_FILE_CORRUPT_);  //CO20200624
      }
      for(uint i=0;i<aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size();i++)
        if(LDEBUG) cerr << "DEBUG= " << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.at(i) << endl;
    } else {
      aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.clear();
      aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.clear();
    }

    // the rest for GEOM
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_IMPLICIT]")) aimsflags.KBIN_AIMS_GEOM_MODE.push("IMPLICIT");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]PROTOTYPE=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE.push("PROTOTYPE");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXTERNAL]")) aimsflags.KBIN_AIMS_GEOM_MODE.push("EXTERNAL");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]FILE=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE.push("FILE");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]COMMAND=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE.push("COMMAND");

    // VOLUMES
    aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.clear();
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]VOLUME=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.push("EQUAL_EQUAL");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]VOLUME+=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.push("PLUS_EQUAL");
    // [OBSOLETE] if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]VOLUME-=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.push("MINUS_EQUAL");
    if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]VOLUME*=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.push("MULTIPLY_EQUAL");
    // [OBSOLETE] if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_FILE]VOLUME/=",TRUE)) aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.push("DIVIDE_EQUAL");
    if(aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.xscheme!="") aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.isentry=TRUE;

    // CONVERT_UNIT_CELL stuff
    BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"CONVERT_UNIT_CELL_","CONVERT_UNIT_CELL="); // bypass for getting all "_"
    aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.options2entry(BflowIn,string(_STROPT_+"CONVERT_UNIT_CELL="),aurostd_xoptionMULTI,""); // stack them all
    // // PRIORITIES
    if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") || aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
      aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI",FALSE);
      aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL",FALSE);
      aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT",FALSE);
      aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ",FALSE);
    } // some PRIORITIES
    if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
      aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL",FALSE);
    }

    // DEBUG
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.content_string=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.content_string  << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"MINKOWSKI\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" <<  aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;

    // VOLUMES
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME_EQUAL_EQUAL      =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME=",TRUE);
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME_PLUS_EQUAL       =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME+=",TRUE);
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME_MINUS_EQUAL      =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME-=",TRUE);
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME_MULTIPLY_EQUAL   =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME*=",TRUE);
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME_DIVIDE_EQUAL     =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME/=",TRUE);
    aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.clear();

    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME=",TRUE)) aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.push("EQUAL_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME+=",TRUE)) aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.push("PLUS_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME-=",TRUE)) aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.push("MINUS_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME*=",TRUE)) aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.push("MULTIPLY_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME/=",TRUE)) aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.push("DIVIDE_EQUAL");

    aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"EQUAL_EQUAL",_STROPT_+"VOLUME=","");
    aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"PLUS_EQUAL",_STROPT_+"VOLUME+=","");
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"MINUS_EQUAL",_STROPT_+"VOLUME-=","");
    aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"MULTIPLY_EQUAL",_STROPT_+"VOLUME*=","");
    // [OBSOLETE] aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"DIVIDE_EQUAL",_STROPT_+"VOLUME/=","");
    if(LDEBUG) cerr << "CORMAC STUFF  " << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedscheme(\"EQUAL_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedscheme("EQUAL_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag(\"PLUS_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedscheme(\"PLUS_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedscheme("PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedscheme(\"MULTIPLY_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedscheme("MULTIPLY_EQUAL") << endl;

    if(aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.xscheme!="") aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.isentry=TRUE;

    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (START)" << endl;

    return aimsflags;
  }

} // namespace KBIN

namespace KBIN{
  bool AIMS_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << "DEBUG: KBIN::AIMS_Directory (BEGIN)" << endl;
    string subS,subS1,subS2;
    ostringstream aus;
    string::iterator pos;
    bool Krun=TRUE;

    ifstream FileAFLOWIN;
    string FileNameAFLOWIN;
    string AflowIn;
    FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
    FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
    FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
    //DX20190125 [OBSOLETE] - need to remove null bytes : AflowIn="";char c; while (FileAFLOWIN.get(c)) AflowIn+=c;  // READ _AFLOWIN_ and put into AflowIn
    AflowIn="";char c; while (FileAFLOWIN.get(c)) if(c!='\0'){ AflowIn+=c; }               // READ AFLOW.IN and put into AflowIn //DX20190125 - remove null bytes
    FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
    AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
      aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      return FALSE;
    }
    aflags.QUIET=FALSE;
    _aimsflags aimsflags;
    aimsflags=KBIN::AIMS_Get_AIMSflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags);

    // *********************************************************************************************************************
    // OPERATIONS related to PARTICULAR MACHINES ***************************************************************************

    if(LDEBUG) cerr << "[DEBUG] aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL.getattachedscheme("NAME") << endl; //HE20220309 use machine name
    if(LDEBUG) cerr << "[DEBUG] aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << endl; //HE20220309 use machine name

    // ***************************************************************************
    // Get the KBIN_BIN name
    aurostd::StringstreamClean(aus);
    aus << "00000  MESSAGE KBIN::AIMS_Directory Running KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // ***************************************************************************
    // Some verbose
    if(kflags.KBIN_PHONONS_CALCULATION_APL) aus << "00000  MESSAGE KBIN::AIMS_Directory Running PHONONS_CALCULATION_APL" << Message(__AFLOW_FILE__,aflags) << endl;
    if(kflags.KBIN_PHONONS_CALCULATION_QHA) aus << "00000  MESSAGE KBIN::AIMS_Directory Running PHONONS_CALCULATION_QHA" << Message(__AFLOW_FILE__,aflags) << endl;   //CO20170601
    if(kflags.KBIN_PHONONS_CALCULATION_AAPL) aus << "00000  MESSAGE KBIN::AIMS_Directory Running PHONONS_CALCULATION_AAPL" << Message(__AFLOW_FILE__,aflags) << endl; //CO20170601
    if(kflags.KBIN_PHONONS_CALCULATION_AGL) aus << "00000  MESSAGE KBIN::AIMS_Directory Running PHONONS_CALCULATION_AGL (Debye Model)" << Message(__AFLOW_FILE__,aflags) << endl;
    if(kflags.KBIN_PHONONS_CALCULATION_AEL) aus << "00000  MESSAGE KBIN::AIMS_Directory Running PHONONS_CALCULATION_AEL (Elastic constants)" << Message(__AFLOW_FILE__,aflags) << endl;
    if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) aus << "00000  MESSAGE KBIN::AIMS_Directory Running PHONONS_CALCULATION_FROZSL" << Message(__AFLOW_FILE__,aflags) << endl;
    if(aimsflags.KBIN_AIMS_RUN.flag("GENERATE")) aus << "00000  MESSAGE KBIN::AIMS_Directory Running RUN_GENERATE" << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // ***************************************************************************
    uint ntasks=0;
    ntasks=1; // default
    if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
      ntasks=aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size();
      aus << "00000  MESSAGE Loaded ntasks = " << ntasks << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      for(uint i=0;i<aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size();i++) {
        aus << "00000  MESSAGE task " << i << "/" << ntasks << " in subdirectory " << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.at(i) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // ***************************************************************************
    // start the loop !
    _aflags aflags_backup;aflags_backup=aflags;
    _kflags kflags_backup;kflags_backup=kflags;


    for(uint ixinput=0;ixinput<ntasks;ixinput++) {  // LOOP ixinput
      // declarations
      _xaims xaims;xaims.clear();
      xaims.GEOM_index=ixinput;
      aflags=aflags_backup;kflags=kflags_backup; // load it up
      // some verbose
      if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus << "00000  MESSAGE START loop " << xaims.GEOM_index << "/" << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size() << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(LDEBUG) cerr << XPID << "KBIN::AIMS_Directory: [1]" << xaims.str << endl; 
      // ------------------------------------------
      // now start for each xaims
      Krun=TRUE;  // guess everything is intelligent !
      xaims.Directory=aflags.Directory;
      if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        xaims.Directory=aflags.Directory+"/"+KBIN_SUBDIRECTORIES+aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.at(xaims.GEOM_index);
        aus << "00000  MESSAGE Taking loop directory = " << xaims.Directory << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // check for directory KY CHECK THIS (if Krun=FALSE, everything stops).
      if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        if(aurostd::FileExist(xaims.Directory)) {
          Krun=FALSE; // avoid rerunning
          aus << "00000  MESSAGE Skipping loop directory = " << xaims.Directory << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          // before making it, check it again... NFS problem... check LOCK again
          if(Krun && aurostd::FileExist(xaims.Directory+"/"+_AFLOWLOCK_)) Krun=FALSE;    // to fight against NFS cache
          if(Krun && aurostd::EFileExist(xaims.Directory+"/"+_AFLOWLOCK_)) Krun=FALSE;     // to fight against NFS cache
          if(Krun && aurostd::FileExist(xaims.Directory+"/LLOCK")) Krun=FALSE;     // to fight against NFS cache
          if(Krun && aurostd::EFileExist(xaims.Directory+"/LLOCK")) Krun=FALSE;     // to fight against NFS cache
          if(Krun) {
            aurostd::DirectoryMake(xaims.Directory);
            aus << "00000  MESSAGE Creating loop directory = " << xaims.Directory << " - " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            aurostd::execute("echo \"NNNNN  KBIN LLOCK ASAP for NFS concurrent jobs (aflow"+string(AFLOW_VERSION)+")\" >> "+xaims.Directory+"/LLOCK");
          }
        }
      }


      if(Krun) {
        aflags.Directory=xaims.Directory; // so we are set ! since there are plenty of routines with aflags.Directory inside
        aus << "00000  MESSAGE Performing loop directory = " << xaims.Directory << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // ------------------------------------------
      // do the flags
      if(LDEBUG) cerr << XPID << "KBIN::AIMS_Directory: [2]" << xaims.str << endl;
      aimsflags.KBIN_AIMS_CONTROL_VERBOSE=TRUE; // ALWAYS

      // produce BEFORE NOMIX
      if(!kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
        if(Krun) Krun=(Krun && KBIN::AIMS_Produce_INPUT(xaims,AflowIn,FileMESSAGE,aflags,kflags,aimsflags));
        if(Krun) Krun=(Krun && KBIN::AIMS_Modify_INPUT(xaims,FileMESSAGE,aflags,kflags,aimsflags));
        if(Krun && kflags.KBIN_QSUB) Krun=(Krun && KBIN::QSUB_Extract(xaims.xqsub,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags));
        if(Krun && kflags.KBIN_QSUB_MODE1) Krun=(Krun && KBIN::QSUB_Extract_Mode1(xaims.xqsub,FileMESSAGE,aflags,kflags));
        if(Krun && kflags.KBIN_QSUB_MODE2) Krun=(Krun && KBIN::QSUB_Extract_Mode2(xaims.xqsub,FileMESSAGE,aflags,kflags));
        if(Krun && kflags.KBIN_QSUB_MODE3) Krun=(Krun && KBIN::QSUB_Extract_Mode3(xaims.xqsub,FileMESSAGE,aflags,kflags));
      }

      // ***************************************************************************
      // READY TO RUN
      if(Krun) {
        if(LDEBUG) cerr << XPID << "KBIN::AIMS_Directory: [3]" << endl;
        if(LDEBUG) cerr << xaims.str << endl;
        bool Krun=true;
        ostringstream aus;
        // ***************************************************************************
        // directory check
        ifstream DirectoryStream;
        DirectoryStream.open(xaims.Directory.c_str(),std::ios::in);
        if(!DirectoryStream) {
          //   aus << "EEEEE  DIRECTORY_NOT_FOUND = " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "XXXXX  MAKING DIRECTORY = " << xaims.Directory << "  " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
          string str="mkdir "+xaims.Directory;
          system(str.c_str());
        }
        // some check
        if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
          //    aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__,aflags) << endl;
          //    aurostd::PrintMessageStream(aus,XHOST.QUIET);
          //    return FALSE;
        }
        // ***************************************************************************
        // DO THE SYMMETRY NEIGHBORS CALCULATION
        //if(!kflags.KBIN_PHONONS_CALCULATION_FROZSL)
        //DX
        if(!(kflags.KBIN_PHONONS_CALCULATION_FROZSL || 
              kflags.KBIN_PHONONS_CALCULATION_APL ||
              kflags.KBIN_PHONONS_CALCULATION_QHA||     //CO20170601
              kflags.KBIN_PHONONS_CALCULATION_AAPL) ||  //CO20170601
            aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN ) {  //CO, do internally
          //DX
          if(Krun) Krun=KBIN_StepSymmetryPerform(xaims.str,AflowIn,FileMESSAGE,aflags,kflags,TRUE,cout); // DO THE SYMMETRY CALCULATION
          //DX20210122 [OBSOLETE - function doesn't calculate anything, removed] if(Krun) Krun=StepNeighborsPerform(xaims.str,AflowIn,FileMESSAGE,aflags,kflags); // DO THE NEIGHBORS CALCULATION
          //DX
          //cerr << "KBIN GEN SYMMETRY OF AFLOWIN: " << aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN << endl;
          if(aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN){
            return Krun;
          }
          //DX
        }
        // AIMS AIMS WRITE
        //   if(Krun) Krun=(Krun && KBIN::AIMS_Write_INPUT(xaims,aimsflags));
        // ***************************************************************************
        // AIMS INPUT FILES ARE DONE, NOW WE CAN USE OR MODYFYING THEM
        if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_NOTUNE.isentry) {
          aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]NOTUNE, no tuning xCARs - ";
          aus << xaims.Directory << " - K=[" << xaims.str.kpoints_k1 << " " << xaims.str.kpoints_k2 << " " << xaims.str.kpoints_k3 << " " << xaims.str.kpoints_kmax << "] - ";
          aus << XHOST.hostname << " - " << aflow_get_time_string() << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        // ***************************************************************************
        // AIMS HOW TO RUN ??
        // GENERATE ONLY -------------------------------------------------------------
        if(aimsflags.KBIN_AIMS_RUN.flag("GENERATE")) {
          KBIN::AIMS_Write_INPUT(xaims,aimsflags); // AIMS AIMS WRITE
          aus << "00000  MESSAGE AIMS generation files ONLY " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
        } else {
          // RUN SOMETHING
          if(kflags.KBIN_PHONONS_CALCULATION_APL) {  // RUN PHONONS APL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_APL " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          //CO START 20170601
          if(kflags.KBIN_PHONONS_CALCULATION_QHA) {  // RUN PHONONS QHA ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_QHA " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if(kflags.KBIN_PHONONS_CALCULATION_AAPL) {  // RUN PHONONS AAPL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AAPL " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          //CO END 20170601
          if(kflags.KBIN_PHONONS_CALCULATION_AGL) {  // RUN PHONONS AGL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AGL (Debye Model) " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if(kflags.KBIN_PHONONS_CALCULATION_AEL) {  // RUN PHONONS AEL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AEL (Elastic constants) " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {  // RUN PHONONS FROZSL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_FROZSL " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }

          // ***************************************************************************
          // READY TO RUN
          if(Krun) {   // survived all troubles
            // ***************************************************************************
            // START
            if(LDEBUG) cerr << XPID << "KBIN::AIMS_Directory: [4]" << xaims.str << endl;
            // ***************************************************************************
            // PRESCRIPT
            if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
              KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_PRESCRIPT_COMMAND,DEFAULT_AFLOW_PRESCRIPT_OUT);
            // ***************************************************************************
            // PHONONIC PHONONIC PHONONIC
            if(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL) { //CO20170601
              _xinput xinput(xaims);
              readModulesFromAflowIn(AflowIn, kflags, xinput);  //ME20181027
              _xflags xflags(aimsflags);
              //ME20200107 - Wrap in a try statement so that faulty APL runs don't kill other post-processing
              try {
                KBIN::RunPhonons_APL(xinput,AflowIn,aflags,kflags,xflags,FileMESSAGE);  //now it's general
              } catch (aurostd::xerror e) {
                pflow::logger(e.whereFileName(), e.whereFunction(), e.buildMessageString(), aflags.Directory, FileMESSAGE, std::cout, _LOGGER_ERROR_);
              }
              //KBIN::RunPhonons_APL(xaims,AflowIn,aflags,kflags,aimsflags,FileMESSAGE);
            }
            //[MAKE XINPUT] here CORMAC WILL POOP an IF=>AGL
            //[MAKE XINPUT]if(kflags.KBIN_PHONONS_CALCULATION_AGL==TRUE) {
            //[MAKE XINPUT]  KBIN::VASP_RunPhonons_AGL(xaims,AflowIn,aflags,kflags,aimsflags,FileMESSAGE);
            //[MAKE XINPUT]}
            //[MAKE XINPUT]if(kflags.KBIN_PHONONS_CALCULATION_AEL==TRUE) {
            //[MAKE XINPUT]  KBIN::VASP_RunPhonons_AEL(xaims,AflowIn,aflags,kflags,aimsflags,FileMESSAGE);
            //[MAKE XINPUT]}
            //[MAKE XINPUT]if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
            //[MAKE XINPUT]  KBIN::VASP_RunPhonons_FROZSL(xaims,AflowIn,aflags,kflags,aimsflags,FileMESSAGE);
            //[MAKE XINPUT]  //  return Krun;
            //[MAKE XINPUT]}
            if(LDEBUG) cerr << XPID << "KBIN::AIMS_Directory: [5] xaims.str.species.size()=" << xaims.str.species.size() << endl;
            if(LDEBUG) for(uint i=0;i<xaims.str.species.size();i++) cerr << XPID << "KBIN::AIMS_Directory: [5] xaims.str.species.at(i)=[" << xaims.str.species.at(i) << "]" << endl;
            if(LDEBUG) cerr << XPID << "KBIN::AIMS_Directory: [6]" << xaims.str << endl;
            // --------------------------------------------------------------------------------------------------------------------
            // ***************************************************************************
            // POSTSCRIPT
            if(!aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT"))
              if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
                KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
            // ***************************************************************************
          }
        }
      }
      // **********
      // some verbose
      if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus << "00000  MESSAGE END loop " << xaims.GEOM_index << "/" << aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size()
          << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE END loop in directory =" << xaims.Directory << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // compress the subdirectories
        if(Krun && kflags.KZIP_COMPRESS) KBIN::CompressDirectory(aflags,kflags);
      }
      aflags=aflags_backup;kflags=kflags_backup; // RESTORE
    } // LOOP ixinput
    // ***************************************************************************
    aflags=aflags_backup;kflags=kflags_backup; // RESTORE
    // POSTSCRIPT
    if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT"))
      if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
        KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
    // ***************************************************************************
    FileAFLOWIN.clear();FileAFLOWIN.close();
    return Krun;
  }
} // namespace KBIN

#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
