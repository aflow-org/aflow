// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to run VASP in KBIN
// Stefano Curtarolo - 2007-2011 Duke

#ifndef _AFLOW_KALIEN_CPP
#define _AFLOW_KALIEN_CPP

#include "aflow.h"

namespace ALIEN {
  _alienflags Get_Alienflags_from_AflowIN(string& AflowIn) {
    _alienflags alienflags;
    // what to run
    alienflags.KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG=
      aurostd::substring2bool(AflowIn,"[ALIEN_COMMAND]START",TRUE) &&
      aurostd::substring2bool(AflowIn,"[ALIEN_COMMAND]STOP",TRUE);
    if(alienflags.KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG==FALSE) {
      alienflags.KBIN_ALIEN_COMMAND_BINARY_FLAG =
        aurostd::substring2bool(AflowIn,"[ALIEN_COMMAND]",TRUE);
      if(alienflags.KBIN_ALIEN_COMMAND_BINARY_FLAG)
        alienflags.KBIN_ALIEN_COMMAND_BINARY_VALUE =
          aurostd::substring2string(AflowIn,"[ALIEN_COMMAND]",1,FALSE);
    } else {
      alienflags.KBIN_ALIEN_COMMAND_BINARY_FLAG = FALSE;
      alienflags.KBIN_ALIEN_COMMAND_BINARY_VALUE=aurostd::substring2string(AflowIn,"[ALIEN_COMMAND]START","[ALIEN_COMMAND]STOP",0);
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringEXPLICIT(AflowIn,alienflags.KBIN_ALIEN_COMMAND_BINARY_VALUE,"[ALIEN_COMMAND]START","[ALIEN_COMMAND]STOP");
    }
    // HOW TO RUN
    // just run it
    // FORCE OPTIONS
    alienflags.KBIN_ALIEN_FORCE_OPTION_NOTUNE   =
      aurostd::substring2bool(AflowIn,"[ALIEN_FORCE_OPTION_NOTUNE]")   ||
      aurostd::substring2bool(AflowIn,"[ALIEN_FORCE_OPTION]NOTUNE",TRUE);
    // FORCE OPTIONS SOMETHING
    alienflags.KBIN_ALIEN_FORCE_OPTION_SOMETHING   =
      aurostd::substring2bool(AflowIn,"[ALIEN_FORCE_OPTION_SOMETHING]")   ||
      aurostd::substring2bool(AflowIn,"[ALIEN_FORCE_OPTION]SOMETHING",TRUE);
    // INPUT FILES
    alienflags.KBIN_ALIEN_INPUT_FILE                      =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE]");          
    alienflags.KBIN_ALIEN_INPUT_MODE_EXPLICIT             =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_EXPLICIT]");
    alienflags.KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP   =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_EXPLICIT]START") &&
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_EXPLICIT]STOP");
    alienflags.KBIN_ALIEN_INPUT_MODE_IMPLICIT             =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_IMPLICIT]");
    alienflags.KBIN_ALIEN_INPUT_MODE_EXTERNAL             =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_EXTERNAL]");
    alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG       =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE]FILE=",TRUE);
    if(alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG)
      alienflags.KBIN_ALIEN_INPUT_FILE_FILE_VALUE =
        aurostd::substring2string(AflowIn,"[ALIEN_INPUT_FILE]FILE=",1,TRUE);
    alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG    =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE]COMMAND=",TRUE);
    if(alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG)
      alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE =
        aurostd::substring2string(AflowIn,"[ALIEN_INPUT_FILE]COMMAND=",1,FALSE);
    alienflags.KBIN_ALIEN_INPUT_MODE_INPUT_FLAG    =
      aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_NAME]INPUT=",TRUE);
    if(alienflags.KBIN_ALIEN_INPUT_MODE_INPUT_FLAG)
      alienflags.KBIN_ALIEN_INPUT_MODE_INPUT_VALUE =
        aurostd::substring2string(AflowIn,"[ALIEN_INPUT_FILE_NAME]INPUT=",1,TRUE);
    alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG    =
      aurostd::substring2bool(AflowIn,"[ALIEN_OUTPUT_FILE_NAME]OUTPUT=",TRUE);
    if(alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG)
      alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE =
        aurostd::substring2string(AflowIn,"[ALIEN_OUTPUT_FILE_NAME]OUTPUT=",1,TRUE);

    alienflags.KBIN_ALIEN_INPUT_MODE_EXTERNAL=FALSE;
    if(alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG) alienflags.KBIN_ALIEN_INPUT_MODE_EXTERNAL=TRUE;
    if(alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG) alienflags.KBIN_ALIEN_INPUT_MODE_EXTERNAL=TRUE;

    return alienflags;
  }
} // namespace ALIEN


// ***************************************************************************
namespace ALIEN {
  bool Run_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    // bool KBIN_MPI_LOCAL;KBIN_MPI_LOCAL=MPI;
    string subS,subS1,subS2;
    ostringstream aus;
    string::iterator pos;
    bool Krun=TRUE;
    _alienflags alienflags;
    _xalien xalien;
    xalien.clear();

    ifstream FileAFLOWIN;
    string FileNameAFLOWIN;
    string AflowIn;
    FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
    FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
    FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
    //DX20190125 [OBSOLETE] - need to remove null bytes : AflowIn="";char c; while (FileAFLOWIN.get(c)) AflowIn+=c;               // READ AFLOW.IN and put into AflowIn
    AflowIn="";char c; while (FileAFLOWIN.get(c)) if(c!='\0'){ AflowIn+=c; }               // READ AFLOW.IN and put into AflowIn //DX20190125 - remove null bytes
    FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
    if(!FileAFLOWIN) {                                                                                      // ******* Aflow.in does not exist
      aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      return FALSE;
    }
    aflags.QUIET=FALSE;
    alienflags=ALIEN::Get_Alienflags_from_AflowIN(AflowIn);
    if(alienflags.KBIN_ALIEN_COMMAND_BINARY_FLAG) kflags.KBIN_BIN=alienflags.KBIN_ALIEN_COMMAND_BINARY_VALUE;

    // ***************************************************************************
    // Get the KBIN_BIN name
    aurostd::StringstreamClean(aus);
    aus << "00000  MESSAGE ALIEN::Directory Running KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // ***************************************************************************
    Krun=TRUE;  // guess everything is intelligent !
    xalien.Directory=aflags.Directory;
    if(Krun) Krun=(Krun && ALIEN::Produce_INPUT(xalien,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,alienflags));
    if(Krun) Krun=(Krun && ALIEN::Modify_INPUT(xalien,FileMESSAGE,aflags,alienflags));
    if(Krun && kflags.KBIN_QSUB) Krun=(Krun && KBIN::QSUB_Extract(xalien.xqsub,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags));
    if(Krun && kflags.KBIN_QSUB_MODE1) Krun=(Krun && KBIN::QSUB_Extract_Mode1(xalien.xqsub,FileMESSAGE,aflags,kflags));
    if(Krun && kflags.KBIN_QSUB_MODE2) Krun=(Krun && KBIN::QSUB_Extract_Mode2(xalien.xqsub,FileMESSAGE,aflags,kflags));
    if(Krun && kflags.KBIN_QSUB_MODE3) Krun=(Krun && KBIN::QSUB_Extract_Mode3(xalien.xqsub,FileMESSAGE,aflags,kflags));
    // ***************************************************************************
    // READY TO RUN
    if(Krun) {
      bool Krun=true;
      stringstream aus_exec;
      // ***************************************************************************
      // directory check
      ifstream DirectoryStream;
      DirectoryStream.open(xalien.Directory.c_str(),std::ios::in);
      if(!DirectoryStream) {
        aus << "XXXXX  MAKING DIRECTORY = " << xalien.Directory << "  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
        string str="mkdir "+xalien.Directory;
        system(str.c_str());
      }
      // ***************************************************************************
      // READY TO RUN
      if(Krun==TRUE) {   // survived all troubles
        // ***************************************************************************
        // START
        // ***************************************************************************
        // PRESCRIPT
        if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
          KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_PRESCRIPT_COMMAND,DEFAULT_AFLOW_PRESCRIPT_OUT);
        // ***************************************************************************
        // RUN
        ALIEN::Write_INPUT(xalien); // ALIEN WRITE INPUT
        aus << "AAAAA  ALIEN RUN - " <<  xalien.Directory << " - \"" << kflags.KBIN_BIN << "\" - " << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // single line command
        if(alienflags.KBIN_ALIEN_COMMAND_BINARY_FLAG==TRUE) {
          aus_exec.clear();
          aus_exec << "cd " << xalien.Directory << endl;
          aus_exec << "rm -f " << alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE << endl;
          aus_exec << kflags.KBIN_BIN << " >> " << alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE << endl;
          //   cout << aus_exec.str()  << endl;
          aurostd::execute(aus_exec);
        }
        // start stop command
        if(alienflags.KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG==TRUE) {
          aus_exec.clear();
          aus_exec << "cd " << xalien.Directory << endl;
          aus_exec << "rm -f " << alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE << endl;
          aurostd::string2file(alienflags.KBIN_ALIEN_COMMAND_BINARY_VALUE,xalien.Directory+"/alien.bin");
          aus_exec << "sh ./alien.bin >> " << alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE << endl;
          aus_exec << "rm -f ./alien.bin" << endl;
          //   cout << aus_exec.str()  << endl;
          aurostd::execute(aus_exec);
        }
        //   ALIEN::Run(xalien,aflags,kflags,alienflags,FileMESSAGE);
        // ***************************************************************************
        // POSTSCRIPT
        if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
          KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
        // ***************************************************************************
      }
    }
    FileAFLOWIN.clear();FileAFLOWIN.close();
    return Krun;
  }
} // namespace ALIEN


#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
