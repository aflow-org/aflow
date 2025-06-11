// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to prepare ALIEN input files
// Stefano Curtarolo - 2008 Duke

#ifndef _AFLOW_IALIEN_CPP
#define _AFLOW_IALIEN_CPP

#include "aflow.h"

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INPUT
namespace ALIEN {
  bool Produce_INPUT(_xalien& xalien,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_alienflags &alienflags) {
    bool Krun=TRUE;
    if(Krun) Krun=(Krun && ALIEN::Produce_INPUT_FILE(xalien,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,alienflags));     // produce INPUT
    return Krun;
  }
} // namespace ALIEN

namespace ALIEN {
  bool Modify_INPUT(_xalien& xalien,ofstream &FileMESSAGE,_aflags &aflags,_alienflags &alienflags) {
    bool Krun=TRUE;
    if(Krun) Krun=(Krun && ALIEN::Modify_INPUT_FILE(xalien,FileMESSAGE,aflags,alienflags));
    return Krun;
  }
} // namespace ALIEN

namespace ALIEN {
  bool Write_INPUT(_xalien& xalien) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool Krun=TRUE;
    ifstream DirectoryStream;
    DirectoryStream.open(xalien.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      ostringstream aus;
      //   aus << "EEEEE  DIRECTORY_NOT_FOUND = " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "XXXXX  MAKING DIRECTORY = " << xalien.Directory << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+xalien.Directory;
      system(str.c_str());
    }
    DirectoryStream.close();
    // ALIEN ALIEN WRITE
    if(Krun && xalien.INPUT_generated) Krun=(Krun && aurostd::stringstream2file(xalien.INPUT,string(xalien.Directory+"/"+xalien.INPUT_FILE_NAME)));
    // ALIEN BACKUP ALIEN WRITE
    if(Krun && xalien.INPUT_generated && xalien.INPUT_changed)  Krun=(Krun && aurostd::stringstream2file(xalien.INPUT_orig,string(xalien.Directory+"/"+xalien.INPUT_FILE_NAME+".orig")));
    return Krun;
  }
} // namespace ALIEN

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INPUT
namespace ALIEN {
  bool Produce_INPUT_FILE(_xalien& xalien,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_alienflags &alienflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()) {;} // phony, just to keep AflowIn busy
    if(!kflags.AFLOW_MODE_ALIEN) {cerr << "ALIEN::Produce_INPUT: should kflags.AFLOW_MODE_ALIEN be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    xalien.INPUT.str(std::string());xalien.INPUT.clear();
    xalien.INPUT_orig.str(std::string());xalien.INPUT_orig.clear();
    xalien.INPUT_generated=FALSE;
    xalien.INPUT_changed=FALSE;
    // IMPLICIT or EXPLICIT or EXTERNAL for INPUT
    Krun=(Krun && (alienflags.KBIN_ALIEN_INPUT_MODE_IMPLICIT ||
          alienflags.KBIN_ALIEN_INPUT_MODE_EXPLICIT ||
          alienflags.KBIN_ALIEN_INPUT_MODE_EXTERNAL));
    if(!Krun)  {
      aurostd::StringstreamClean(aus);
      //    aus << "EEEEE  [ALIEN_INPUT_MODE_IMPLICIT] or [ALIEN_INPUT_MODE_EXPLICIT] or [ALIEN_INPUT_MODE_EXPLICIT] must be specified "  << Message(__AFLOW_FILE__,aflags) << endl;
      //    aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
      aus << "AAAAA  ALIEN, no input file specified "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      xalien.INPUT_generated=FALSE;            // no need of generating input
      Krun=TRUE;                               // and run anyway
      return Krun;
    }
    // EXPLICIT **************************************************
    if(Krun && alienflags.KBIN_ALIEN_INPUT_MODE_EXPLICIT) {  // [ALIEN_INPUT_MODE_EXPLICIT] construction
      if(alienflags.KBIN_ALIEN_INPUT_FILE && !alienflags.KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP) {
        aus << "00000  MESSAGE INPUT  generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
        xalien.INPUT.str(aurostd::substring2string(FileAFLOWIN,"[ALIEN_INPUT_FILE]",0));
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xalien.INPUT,"[ALIEN_INPUT_FILE]");
        // DO SOME LOADING
        // DEBUG xalien.str=xstructure(xalien.INPUT,IOALIEN);  // load structure
      } else if(!alienflags.KBIN_ALIEN_INPUT_FILE && alienflags.KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP) {
        aus << "00000  MESSAGE INPUT  generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
        xalien.INPUT.str(aurostd::substring2string(FileAFLOWIN,"[ALIEN_INPUT_FILE_EXPLICIT]START","[ALIEN_INPUT_FILE_EXPLICIT]STOP",0));
        //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[ALIEN_INPUT_FILE_EXPLICIT]STOP"))
          //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xalien.INPUT,"[ALIEN_INPUT_FILE_EXPLICIT]START","[ALIEN_INPUT_FILE_EXPLICIT]STOP");
        // DO SOME LOADING
        // DEBUG  xalien.str=xstructure(xalien.INPUT,IOALIEN);   // load structure
      } else {
        aus << "EEEEE  [ALIEN_INPUT_FILE_EXPLICIT] do not confuse aflow !!"  << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [ALIEN_INPUT_FILE_EXPLICIT] Possible modes "  << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] INPUT EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[ALIEN_INPUT_FILE_EXPLICIT]" << endl;
        aus << "[ALIEN_INPUT_FILE]INPUT of the structure example" << endl;
        aus << "[ALIEN_INPUT_FILE]-98.5397" << endl;
        aus << "[ALIEN_INPUT_FILE]   4.18890 0.00000 0.00000" << endl;
        aus << "[ALIEN_INPUT_FILE]  -2.09445 3.62769 0.00000" << endl;
        aus << "[ALIEN_INPUT_FILE]   0.00000 0.00000 5.12300" << endl;
        aus << "[ALIEN_INPUT_FILE]2 4" << endl;
        aus << "[ALIEN_INPUT_FILE]Direct" << endl;
        aus << "[ALIEN_INPUT_FILE]0.33333 0.66666 0.25000 Au" << endl;
        aus << "[ALIEN_INPUT_FILE]0.66666 0.33333 0.75000 Au" << endl;
        aus << "[ALIEN_INPUT_FILE]0.00000 0.00000 0.00000 Ti" << endl;
        aus << "[ALIEN_INPUT_FILE]0.00000 0.00000 0.50000 Ti" << endl;
        aus << "[ALIEN_INPUT_FILE]0.33333 0.66666 0.75000 Ti" << endl;
        aus << "[ALIEN_INPUT_FILE]0.66666 0.33333 0.25000 Ti" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] INPUT EXPLICIT MODE with START/STOP" << endl;
        aus << "[ALIEN_INPUT_FILE_EXPLICIT]" << endl;
        aus << "[ALIEN_INPUT_FILE_EXPLICIT]START" << endl;
        aus << "INPUT of the structure example with START/STOP" << endl;
        aus << "-98.5397" << endl;
        aus << "   4.18890 0.00000 0.00000" << endl;
        aus << "  -2.09445 3.62769 0.00000" << endl;
        aus << "   0.00000 0.00000 5.12300" << endl;
        aus << "2 4" << endl;
        aus << "Direct" << endl;
        aus << "0.33333 0.66666 0.25000 Au" << endl;
        aus << "0.66666 0.33333 0.75000 Au" << endl;
        aus << "0.00000 0.00000 0.00000 Ti" << endl;
        aus << "0.00000 0.00000 0.50000 Ti" << endl;
        aus << "0.33333 0.66666 0.75000 Ti" << endl;
        aus << "0.66666 0.33333 0.25000 Ti" << endl;
        aus << "[ALIEN_INPUT_FILE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [ALIEN_INPUT_FILE_EXPLICIT] Note "  << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [ALIEN_INPUT_FILE_EXPLICIT]START must be present and no [ALIEN_INPUT_FILE]"  << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [ALIEN_INPUT_FILE_EXPLICIT]STOP  must be present and no [ALIEN_INPUT_FILE]"  << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  or [ALIEN_INPUT_FILE] present and NO START/STOP"  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
    }
    // IMPLICIT **************************************************
    if(Krun && alienflags.KBIN_ALIEN_INPUT_MODE_IMPLICIT) {
      aus << "EEEEE  [ALIEN_INPUT_FILE_IMPLICIT] is not supporded      "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [ALIEN_INPUT_FILE_EXPLICIT] is supported mode     "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [ALIEN_INPUT_FILE_EXTERNAL] is supported mode     "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
      Krun=FALSE;  // ALIEN_INPUT_MODE_IMPLICIT
      return Krun;
    }
    // EXTERNAL **************************************************
    if(Krun && alienflags.KBIN_ALIEN_INPUT_MODE_EXTERNAL) {  // [ALIEN_INPUT_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE INPUT  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
      if(alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG && alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG) {
        aus << "EEEEE   [ALIEN_INPUT_MODE]FILE=  and  [ALIEN_INPUT_MODE]COMMAND=  can not be used together "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
      if(!alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG && (alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG || !alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG)) {
        if(alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG) {
          file=alienflags.KBIN_ALIEN_INPUT_FILE_FILE_VALUE;
          aus << "00000  MESSAGE INPUT  generation EXTERNAL from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
        } else {
          file=ALIEN_EXTERNAL_INPUT_DEFAULT;
          aus << "00000  MESSAGE INPUT  generation EXTERNAL from DEFAULT file=" << ALIEN_EXTERNAL_INPUT_DEFAULT << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR INPUT file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR INPUT file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
          Krun=FALSE;
          return Krun;
        }
        xalien.INPUT << aurostd::file2string(file);
        // DO SOME LOADING
        // xalien.str=xstructure(xalien.INPUT,IOALIEN);  // load structure
      }
      if(alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG && !alienflags.KBIN_ALIEN_INPUT_FILE_FILE_FLAG) {
        file=alienflags.KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE;
        aus << "00000  MESSAGE INPUT  generation EXTERNAL from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
        file=file+" > ./_aflow_INPUT."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp  //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_INPUT."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name  //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR INPUT file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR INPUT file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
          Krun=FALSE;
          return Krun;
        }
        xalien.INPUT << aurostd::file2string(file);       // load INPUT
        // DO SOME LOADING
        // xalien.str=xstructure(xalien.INPUT,IOALIEN);              // load structure
        file="rm -f ./_aflow_INPUT."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";     // remove temp //CO20200502 - threadID
        aurostd::execute(file);                          // remove temp
      }
    }
    // INPUT DONE **************************************************
    xalien.INPUT_orig << xalien.INPUT.str();
    xalien.INPUT_generated=TRUE;
    xalien.INPUT_changed=FALSE;

    if(!alienflags.KBIN_ALIEN_INPUT_MODE_INPUT_FLAG) {
      xalien.INPUT_FILE_NAME=alienflags.KBIN_ALIEN_INPUT_MODE_INPUT_VALUE;
      aus << "00000  MESSAGE INPUT  generation to default filename " << xalien.INPUT_FILE_NAME << "   " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
    } else {
      xalien.INPUT_FILE_NAME=alienflags.KBIN_ALIEN_INPUT_MODE_INPUT_VALUE;
      aus << "00000  MESSAGE INPUT  generation to filename " << xalien.INPUT_FILE_NAME << "   " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
    }
    if(!alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG) {
      xalien.OUTPUT_FILE_NAME=alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE;
      aus << "00000  MESSAGE OUTPUT  generation to default filename " << xalien.OUTPUT_FILE_NAME << "   " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
    } else {
      xalien.OUTPUT_FILE_NAME=alienflags.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE;
      aus << "00000  MESSAGE OUTPUT  generation to filename " << xalien.OUTPUT_FILE_NAME << "   " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
    }
    // done produced and modified
    return Krun;
  };  // ALIEN::Produce_INPUT
} // namespace ALIEN


namespace ALIEN {
  bool Modify_INPUT_FILE(_xalien& xalien,ofstream &FileMESSAGE,_aflags &aflags,_alienflags &alienflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(xalien.INPUT_generated==FALSE) {
      // it is not necessary to have generated the input for ALIEN
      return Krun;
    }
    //   if(xalien.INPUT_generated==FALSE) {
    //     aus << "EEEEE  ALIEN::Modify_INPUT: can`t modify INPUT if it does not exist ! "  << Message(__AFLOW_FILE__,aflags) << endl;
    //     aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
    //     Krun=FALSE;
    //     return Krun;
    //   }

    xalien.INPUT_generated=FALSE;

    xalien.INPUT_orig.str(std::string());xalien.INPUT_orig.clear();
    xalien.INPUT_orig << xalien.INPUT.str();

    if(Krun && alienflags.KBIN_ALIEN_FORCE_OPTION_SOMETHING) {  // [ALIEN_FORCE_OPTION]FORCE_OPTION_SOMETHING construction
      aus << "00000  MESSAGE INPUT   FORCE_OPTION_SOMETHING " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
      xalien.INPUT.str(std::string());xalien.INPUT.clear();
      // some modification
      xalien.INPUT_generated=TRUE;
      xalien.INPUT_changed=TRUE;
    }
    if(0) {
      cout <<"DEBUG artificially modified to see if it works" << endl;
      xalien.INPUT << "modified" << endl;
      xalien.INPUT_generated=TRUE;
      xalien.INPUT_changed=TRUE;
    }
    // INPUT done
    if(Krun && alienflags.KBIN_ALIEN_FORCE_OPTION_NOTUNE==FALSE) {
      if(0) {
        aus << "00000  MESSAGE Option XXXXX"  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
    }
    // ------------------------------------
    // end
    xalien.INPUT_generated=TRUE;
    return Krun;
  }; // ALIEN::Modify_INPUT
} // namespace ALIEN


// ---------------------------------------------------------------------------------------------------------------------------------------------------------

#endif  // _AFLOW_IALIEN_CPP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

