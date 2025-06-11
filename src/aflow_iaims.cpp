// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to prepare AIMS input files
// Stefano Curtarolo - 2007 Duke
// Corey Oses - 2017 Duke

#ifndef _AFLOW_IAIMS_CPP
#define _AFLOW_IAIMS_CPP

#include "aflow.h"
#define _controlinpad_ 48
#define _IAIMS_DOUBLE2STRING_PRECISION_ 7

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INPUT
namespace KBIN {
  bool AIMS_Produce_INPUT(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags) {
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::AIMS_Produce_INPUT():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    bool Krun=TRUE;
    if(Krun) Krun=(Krun && KBIN::AIMS_Produce_GEOM(xaims,AflowIn,FileMESSAGE,aflags,kflags,aimsflags));     // produce GEOM before KPOINTS
    if(Krun) Krun=(Krun && KBIN::AIMS_Produce_CONTROL(xaims,AflowIn,FileMESSAGE,aflags,kflags,aimsflags));
    return Krun;
  }
}

namespace KBIN {
  bool AIMS_Modify_INPUT(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags) {
    bool Krun=TRUE;
    if(Krun) Krun=(Krun && KBIN::AIMS_Modify_CONTROL(xaims,FileMESSAGE,aflags,kflags,aimsflags));
    return Krun;
  }
}

namespace KBIN {
  bool AIMS_Write_INPUT(_xaims& xaims,_aimsflags &aimsflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ifstream DirectoryStream;
    DirectoryStream.open(xaims.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      ostringstream aus;
      aus << "XXXXX  MAKING DIRECTORY = " << xaims.Directory << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+xaims.Directory;
      system(str.c_str());
    }
    DirectoryStream.close();
    bool Krun=TRUE;
    // AIMS AIMS WRITE
    if(Krun) Krun=(Krun && KBIN::AIMS_Write_CONTROL(xaims,aimsflags));
    if(Krun) Krun=(Krun && KBIN::AIMS_Write_GEOM(xaims,aimsflags));
    return Krun;
  }
}

namespace KBIN {
  bool AIMS_Write_CONTROL(_xaims& xaims,_aimsflags &aimsflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ifstream DirectoryStream;
    DirectoryStream.open(xaims.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      ostringstream aus;
      aus << "XXXXX  MAKING DIRECTORY = " << xaims.Directory << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+xaims.Directory;
      system(str.c_str());
    }
    DirectoryStream.close();
    bool Krun=TRUE;
    // AIMS AIMS WRITE
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xaims.CONTROL,string(xaims.Directory+"/control.in")));
    // AIMS BACKUP AIMS WRITE
    if(Krun && xaims.aopts.flag("FLAG::XAIMS_CONTROL_changed"))   Krun=(Krun && aurostd::stringstream2file(xaims.CONTROL_orig,string(xaims.Directory+"/control.in.orig")));
    if(aimsflags.KBIN_AIMS_CONTROL_VERBOSE) {;} // DUMMY
    return Krun;
  }
}

namespace KBIN {
  bool AIMS_Write_GEOM(_xaims& xaims,_aimsflags &aimsflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ifstream DirectoryStream;
    DirectoryStream.open(xaims.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      ostringstream aus;
      aus << "XXXXX  MAKING DIRECTORY = " << xaims.Directory << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+xaims.Directory;
      system(str.c_str());
    }
    DirectoryStream.close();
    bool Krun=TRUE;
    // AIMS AIMS WRITE
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xaims.GEOM,string(xaims.Directory+"/geometry.in")));
    // AIMS BACKUP AIMS WRITE
    if(Krun && xaims.aopts.flag("FLAG::XAIMS_GEOM_changed"))  Krun=(Krun && aurostd::stringstream2file(xaims.GEOM_orig,string(xaims.Directory+"/geometry.in.orig")));
    if(aimsflags.KBIN_AIMS_CONTROL_VERBOSE) {;} // DUMMY
    return Krun;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// CONTROL
namespace KBIN {
  bool AIMS_Produce_CONTROL(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::AIMS_Produce_CONTROL():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_AIMS) {cerr << XPID << "KBIN::AIMS_Produce_CONTROL: should kflags.AFLOW_MODE_AIMS be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    xaims.CONTROL.str(std::string());
    xaims.aopts.flag("FLAG::XAIMS_CONTROL_generated",FALSE);
    xaims.aopts.flag("FLAG::XAIMS_CONTROL_changed",FALSE);

    aus << "00000  MESSAGE CONTROL   generation in " << xaims.Directory << "  " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    bool KBIN_AIMS_CONTROL_MODE_EMPTY=!aimsflags.KBIN_AIMS_CONTROL_MODE.flag("IMPLICIT") && !aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXPLICIT") && !aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXTERNAL");

    // IMPLICIT or EXPLICIT or EXTERNAL for CONTROL
    Krun=(Krun && (aimsflags.KBIN_AIMS_CONTROL_MODE.flag("IMPLICIT") ||
          aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXPLICIT") ||
          aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXTERNAL") || KBIN_AIMS_CONTROL_MODE_EMPTY));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [AIMS_CONTROL_MODE_IMPLICIT] or [AIMS_CONTROL_MODE_EXPLICIT] or [AIMS_CONTROL_MODE_EXPLICIT] must be specified " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // EMPTY ************************************************** CONTROL
    if(Krun && KBIN_AIMS_CONTROL_MODE_EMPTY) {  // [AIMS_CONTROL_MODE_EMPTY] construction
      aus << "00000  MESSAGE CONTROL   generation EMPTY file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.CONTROL << "#AFLOW CONTROL automatically generated" << endl;
    }
    // IMPLICIT ************************************************** CONTROL
    if(Krun && aimsflags.KBIN_AIMS_CONTROL_MODE.flag("IMPLICIT")) {  // [AIMS_CONTROL_MODE_IMPLICIT] construction
      aus << "00000  MESSAGE CONTROL   generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(aimsflags.KBIN_AIMS_CONTROL_FILE.flag("SYSTEM_AUTO")) {
        xaims.CONTROL << "#AFLOW CONTROL automatically generated" << endl;
      }
    }
    // EXPLICIT ************************************************** CONTROL
    if(Krun && aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXPLICIT")) {  // [AIMS_CONTROL_MODE_EXPLICIT] construction
      if(aimsflags.KBIN_AIMS_CONTROL_FILE.flag("KEYWORD") && !aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE CONTROL   generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xaims.CONTROL.str(aurostd::substring2string(AflowIn,"[AIMS_CONTROL_FILE]",0));
        // [OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xaims.CONTROL,"[AIMS_CONTROL_FILE]");
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xaims.CONTROL,"[AIMS_CONTROL_FILE]");
      } else if(!aimsflags.KBIN_AIMS_CONTROL_FILE.flag("KEYWORD") && aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE CONTROL   generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xaims.CONTROL.str(aurostd::substring2string(AflowIn,"[AIMS_CONTROL_MODE_EXPLICIT]START","[AIMS_CONTROL_MODE_EXPLICIT]STOP",0));
        //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[AIMS_CONTROL_MODE_EXPLICIT]STOP"))
          // [OBSOLETE]	aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xaims.CONTROL,"[AIMS_CONTROL_MODE_EXPLICIT]START","[AIMS_CONTROL_MODE_EXPLICIT]STOP");
          //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xaims.CONTROL,"[AIMS_CONTROL_MODE_EXPLICIT]START","[AIMS_CONTROL_MODE_EXPLICIT]STOP");
      } else {
        aus << "EEEEE  [AIMS_CONTROL_MODE_EXPLICIT] do not confuse aflow !!" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [AIMS_CONTROL_MODE_EXPLICIT] Possible modes " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] CONTROL EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[AIMS_CONTROL_MODE_EXPLICIT]" << endl;
        aus << "[AIMS_CONTROL_FILE]xc               pw−lda                #Physical settings" << endl;
        aus << "[AIMS_CONTROL_FILE]spin             none                  #Physical settings" << endl;
        aus << "[AIMS_CONTROL_FILE]relativistic     atomic_zora scalar    #Physical settings" << endl;
        aus << "[AIMS_CONTROL_FILE]sc_accuracy_rho  1E-4                  #SCF settings" << endl;
        aus << "[AIMS_CONTROL_FILE]sc_accuracy_eev  1E-2                  #SCF settings" << endl;
        aus << "[AIMS_CONTROL_FILE]sc_accuracy_etot 1E-5                  #SCF settings" << endl;
        aus << "[AIMS_CONTROL_FILE]k_grid           15 15 15              #k-grid settings" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] CONTROL EXPLICIT MODE with START/STOP" << endl;
        aus << "[AIMS_CONTROL_MODE_EXPLICIT]" << endl;
        aus << "[AIMS_CONTROL_MODE_EXPLICIT]START" << endl;
        aus << "#Physical settings" << endl;
        aus << "xc               pw−lda" << endl;
        aus << "spin             none" << endl;
        aus << "relativistic     atomic_zora scalar" << endl;
        aus << "#SCF settings" << endl;
        aus << "sc_accuracy_rho  1E-4" << endl;
        aus << "sc_accuracy_eev  1E-2" << endl;
        aus << "sc_accuracy_etot 1E-5" << endl;
        aus << "#k-grid settings" << endl;
        aus << "k_grid           15 15 15" << endl;
        aus << "[AIMS_CONTROL_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [AIMS_CONTROL_MODE_EXPLICIT] Note " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [AIMS_CONTROL_MODE_EXPLICIT]START must be present and no [AIMS_CONTROL_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [AIMS_CONTROL_MODE_EXPLICIT]STOP  must be present and no [AIMS_CONTROL_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  or [AIMS_CONTROL_FILE] present and NO START/STOP" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && aimsflags.KBIN_AIMS_CONTROL_MODE.flag("EXTERNAL")) {  // [AIMS_CONTROL_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE CONTROL   generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(aimsflags.KBIN_AIMS_CONTROL_FILE.flag("COMMAND") && aimsflags.KBIN_AIMS_CONTROL_FILE.flag("FILE")) {
        aus << "EEEEE   [AIMS_CONTROL_MODE]FILE=  and  [AIMS_CONTROL_MODE]COMMAND=  can not be used together " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!aimsflags.KBIN_AIMS_CONTROL_FILE.flag("COMMAND") && (aimsflags.KBIN_AIMS_CONTROL_FILE.flag("FILE") || !aimsflags.KBIN_AIMS_CONTROL_FILE.flag("FILE"))) {
        if(aimsflags.KBIN_AIMS_CONTROL_FILE.flag("FILE")) {
          file=aurostd::substring2string(AflowIn,"[AIMS_CONTROL_FILE]FILE=",1,TRUE);
          aus << "00000  MESSAGE CONTROL   generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_AIMS_EXTERNAL_CONTROL;
          aus << "00000  MESSAGE CONTROL   generation from DEFAULT file=" << DEFAULT_AIMS_EXTERNAL_CONTROL << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR CONTROL file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR CONTROL file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xaims.CONTROL << aurostd::file2string(file);
      }
      if(aimsflags.KBIN_AIMS_CONTROL_FILE.flag("COMMAND") && !aimsflags.KBIN_AIMS_CONTROL_FILE.flag("FILE")) {
        file=aurostd::substring2string(AflowIn,"[AIMS_CONTROL_FILE]COMMAND=",1,FALSE);
        aus << "00000  MESSAGE CONTROL   generation from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_CONTROL."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp  //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_CONTROL."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name  //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR CONTROL file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR CONTROL file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xaims.CONTROL << aurostd::file2string(file);       // load CONTROL
        file="rm -f ./_aflow_CONTROL."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";     // remove temp //CO20200502 - threadID
        aurostd::execute(file);                          // remove temp
      }
    }
    // CONTROL DONE **************************************************
    xaims.CONTROL << "#control.in" << endl;
    xaims.CONTROL_orig << xaims.CONTROL.str();
    xaims.aopts.flag("FLAG::XAIMS_CONTROL_generated",TRUE);
    return Krun;
  };  // KBIN::AIMS_Produce_CONTROL
}

namespace KBIN {
  bool AIMS_Modify_CONTROL(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;

    if(Krun && kflags.KBIN_MPI) {
      xaims.NCPUS=kflags.KBIN_MPI_NCPUS;
      if(kflags.KBIN_MPI_AUTOTUNE) {
        aus << "00000  MESSAGE CONTROL-MPI: found AUTOTUNE option " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE CONTROL-MPI: input files WILL be auto-tuned for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        //KBIN::AIMS_MPI_Autotune(xaims,vflags.KBIN_AIMS_CONTROL_VERBOSE);
        //xaims.aopts.flag("FLAG::xaims_CONTROL_changed",TRUE);
      } else {
        aus << "00000  MESSAGE CONTROL-MPI: AUTOTUNE option NOT found! (aflow_iaims.cpp) " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE CONTROL-MPI: input files MUST be appropriate for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // CONVERT_UNIT_CELL
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** CONTROL **************/
      if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.isentry) {                                                              /*************** CONTROL **************/
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=STANDARD_PRIMITIVE - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=STANDARD_CONVENTIONAL - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=NIGGLI - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=MINKOWSKI - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=INCELL - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=COMPACT - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=WIGNERSEITZ - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=CARTESIAN - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=FRACTIONAL - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE"))  aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=PRESERVE - "<< Message(__AFLOW_FILE__,aflags) << endl; //CO
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // ------------------------------------
    // end
    xaims.aopts.flag("FLAG::XAIMS_CONTROL_generated",TRUE);
    return Krun;
  };  // KBIN::AIMS_Produce_CONTROL
}


namespace KBIN {
  bool AIMS_Reread_CONTROL(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(!aurostd::FileExist(xaims.Directory+"/CONTROL")) {
      aus << "EEEEE  KBIN::AIMS_Reread_CONTROL: CONTROL not present in directory: " << xaims.Directory << " - "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    xaims.CONTROL_orig.str(std::string()); xaims.CONTROL_orig << xaims.CONTROL.str();
    xaims.CONTROL.str(std::string()); xaims.CONTROL << aurostd::file2string(xaims.Directory+"/CONTROL"); // DID REREAD
    // xaims.aopts.flag("FLAG::XAIMS_CONTROL_changed",TRUE);
    return Krun;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// GEOM
namespace KBIN {
  bool AIMS_Produce_GEOM(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::AIMS_Produce_GEOM():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_AIMS) {cerr << XPID << "KBIN::AIMS_Produce_GEOM: should kflags.AFLOW_MODE_AIMS be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    xaims.GEOM.str(std::string());xaims.GEOM.clear();
    xaims.GEOM_orig.str(std::string());xaims.GEOM_orig.clear();
    xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",FALSE);
    xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",FALSE);
    //
    // IMPLICIT or EXPLICIT or EXTERNAL for GEOM
    Krun=(Krun && (aimsflags.KBIN_AIMS_GEOM_MODE.flag("IMPLICIT") ||
          aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT") ||
          aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [AIMS_GEOM_MODE_IMPLICIT] or [AIMS_GEOM_MODE_EXPLICIT] or [AIMS_GEOM_MODE_EXPLICIT] must be specified " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // IMPLICIT **************************************************
    if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("IMPLICIT")) {  // [AIMS_GEOM_MODE_IMPLICIT] construction
      aus << "00000  MESSAGE GEOM  generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(!aimsflags.KBIN_AIMS_GEOM_FILE.flag("PROTOTYPE")) {
        aus << "EEEEE  [AIMS_GEOM_FILE] In GEOM_MODE_IMPLICIT you must specify PROTOTYPE " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      } else {
        std::vector<string> tokens,tokens2,atomABC;
        std::string structure,label,parameters="";  // FIX NRL
        vector<double> volumeABC;
        structure=aurostd::substring2string(AflowIn,"[AIMS_GEOM_FILE]PROTOTYPE=",1,TRUE);
        aurostd::string2tokens(structure,tokens,";");
        label=tokens[0];
        for(uint i=1;i<tokens.size();i++) {
          // find SPECIES
          if(aurostd::substring2bool(tokens[i],"SPECIES=",TRUE) || aurostd::substring2bool(tokens[i],"SPECIE=",TRUE)) {
            aurostd::string2tokens(aurostd::substring2string(tokens[i],"=",1,TRUE),tokens2,",");
            for(uint j=0;j<tokens2.size();j++)
              atomABC.push_back(tokens2[j]);
          }
          // find VOLUMES
          if(aurostd::substring2bool(tokens[i],"VOLUMES=",TRUE) || aurostd::substring2bool(tokens[i],"VOLUME=",TRUE)) {
            aurostd::string2tokens(aurostd::substring2string(tokens[i],"=",1,TRUE),tokens2,",");
            for(uint j=0;j<tokens2.size();j++)
              volumeABC.push_back(aurostd::string2utype<double>(tokens2[j]));
          }
        }
        // for(uint j=0;j<atomABC.size();j++) cerr << atomABC.at(j) << endl;
        // for(uint j=0;j<volumeABC.size();j++) cerr << volumeABC.at(j) << endl;
        bool done=FALSE;
        if(atomABC.size()==2 && volumeABC.size()==0) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(GetAtomVolume(atomABC[isp]));} //KBIN::AIMS_PseudoPotential_CleanName(atomABC[isp]))); //CO20200106 - patching for auto-indenting
          xaims.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);
        }
        if(atomABC.size()==2 && volumeABC.size()==1) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(0.0);}
          // xaims.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],0.0,atomABC[1],0.0,volumeABC[0]); // OLD WAY
          xaims.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,volumeABC[0],LIBRARY_MODE_HTQC);
        }
        if(atomABC.size()==2 && volumeABC.size()==2) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(volumeABC[isp]);};
          // xaims.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],volumeABC[0],atomABC[1],volumeABC[1],-1.0); // OLD WAY
          xaims.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);
        }
        if(done==FALSE) {
          aus << "EEEEE  GEOM_MODE_IMPLICIT error in the PROTOTYPE definition" << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "EEEEE  [AIMS_GEOM_FILE]PROTOTYPE=" << structure << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        // done
        xaims.GEOM << xaims.str;
        // cerr << xaims.GEOM.str() << endl;
      }
    }
    // EXPLICIT **************************************************
    if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT")) {  // [AIMS_GEOM_MODE_EXPLICIT] construction
      if(aimsflags.KBIN_AIMS_GEOM_FILE.flag("KEYWORD") && !aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP") && !aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus << "00000  MESSAGE GEOM  generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xaims.GEOM.str(aurostd::substring2string(AflowIn,"[AIMS_GEOM_FILE]",0));
        // [OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xaims.GEOM,"[AIMS_GEOM_FILE]");
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xaims.GEOM,"[AIMS_GEOM_FILE]");
        xaims.str=xstructure(xaims.GEOM,IOAIMS_AUTO);  // load structure
        xaims.GEOM.str(std::string());xaims.GEOM.clear();
        xaims.str.iomode=IOAIMS_GEOM;
        xaims.GEOM << xaims.str;
      } else if(!aimsflags.KBIN_AIMS_GEOM_FILE.flag("KEYWORD") && (aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP") || aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT"))) {
        aus << "00000  MESSAGE GEOM  generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // normal get ONE of ONE
        if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP")) {
          //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]START") &&
          //[SD20220520 - OBSOLETE]    aurostd::substring2bool(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]STOP"))
          xaims.GEOM.str(aurostd::substring2string(AflowIn,"[AIMS_GEOM_MODE_EXPLICIT]START","[AIMS_GEOM_MODE_EXPLICIT]STOP",-1));
            // [OBSOLETE]	  aurostd::ExtractLastToStringstreamEXPLICIT(FileAFLOWIN,xaims.GEOM,"[AIMS_GEOM_MODE_EXPLICIT]START","[AIMS_GEOM_MODE_EXPLICIT]STOP");
            //[SD20220520 - OBSOLETE]aurostd::ExtractLastToStringstreamEXPLICIT(AflowIn,xaims.GEOM,"[AIMS_GEOM_MODE_EXPLICIT]START","[AIMS_GEOM_MODE_EXPLICIT]STOP");
          if(!xaims.GEOM.str().empty()) {xaims.str=xstructure(xaims.GEOM,IOAIMS_AUTO);}   // load structure
        }
        // get ONE of MANY
        if(aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXPLICIT_START_STOP_POINT")) {
          xaims.str=aimsflags.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.at(xaims.GEOM_index);
        }
        // GOT IT
        xaims.GEOM.str(std::string());xaims.GEOM.clear();
        xaims.str.iomode=IOAIMS_GEOM;
        xaims.GEOM << xaims.str;
      } else {
        aus << "EEEEE  [AIMS_GEOM_MODE_EXPLICIT] do not confuse aflow !!" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [AIMS_GEOM_MODE_EXPLICIT] Possible modes " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] GEOM EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[AIMS_GEOM_MODE_EXPLICIT]" << endl;
        aus << "[AIMS_GEOM_FILE]GEOM of the structure example" << endl;
        aus << "[AIMS_GEOM_FILE]-98.5397" << endl;
        aus << "[AIMS_GEOM_FILE]   4.18890 0.00000 0.00000" << endl;
        aus << "[AIMS_GEOM_FILE]  -2.09445 3.62769 0.00000" << endl;
        aus << "[AIMS_GEOM_FILE]   0.00000 0.00000 5.12300" << endl;
        aus << "[AIMS_GEOM_FILE]2 4" << endl;
        aus << "[AIMS_GEOM_FILE]Direct" << endl;
        aus << "[AIMS_GEOM_FILE]0.33333 0.66666 0.25000 Au" << endl;
        aus << "[AIMS_GEOM_FILE]0.66666 0.33333 0.75000 Au" << endl;
        aus << "[AIMS_GEOM_FILE]0.00000 0.00000 0.00000 Ti" << endl;
        aus << "[AIMS_GEOM_FILE]0.00000 0.00000 0.50000 Ti" << endl;
        aus << "[AIMS_GEOM_FILE]0.33333 0.66666 0.75000 Ti" << endl;
        aus << "[AIMS_GEOM_FILE]0.66666 0.33333 0.25000 Ti" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] GEOM EXPLICIT MODE with START/STOP" << endl;
        aus << "[AIMS_GEOM_MODE_EXPLICIT]" << endl;
        aus << "[AIMS_GEOM_MODE_EXPLICIT]START" << endl;
        aus << "GEOM of the structure example with START/STOP" << endl;
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
        aus << "[AIMS_GEOM_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [AIMS_GEOM_MODE_EXPLICIT] Note " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [AIMS_GEOM_MODE_EXPLICIT]START must be present and no [AIMS_GEOM_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [AIMS_GEOM_MODE_EXPLICIT]STOP  must be present and no [AIMS_GEOM_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  or [AIMS_GEOM_FILE] present and NO START/STOP" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("EXTERNAL")) {  // [AIMS_GEOM_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE GEOM  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(aimsflags.KBIN_AIMS_GEOM_FILE.flag("COMMAND") && aimsflags.KBIN_AIMS_GEOM_FILE.flag("FILE")) {
        aus << "EEEEE   [AIMS_GEOM_MODE]FILE=  and  [AIMS_GEOM_MODE]COMMAND=  can not be used together " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!aimsflags.KBIN_AIMS_GEOM_FILE.flag("COMMAND") && (aimsflags.KBIN_AIMS_GEOM_FILE.flag("FILE") || !aimsflags.KBIN_AIMS_GEOM_FILE.flag("FILE"))) {
        if(aimsflags.KBIN_AIMS_GEOM_FILE.flag("FILE")) {
          file=aurostd::substring2string(AflowIn,"[AIMS_GEOM_FILE]FILE=",1,TRUE);
          aus << "00000  MESSAGE GEOM  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_AIMS_EXTERNAL_GEOM;
          aus << "00000  MESSAGE GEOM  generation from DEFAULT file=" << DEFAULT_AIMS_EXTERNAL_GEOM << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR GEOM file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR GEOM file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xaims.GEOM << aurostd::file2string(file);
        xaims.str=xstructure(xaims.GEOM,IOAIMS_GEOM);  // load structure
      }
      if(aimsflags.KBIN_AIMS_GEOM_FILE.flag("COMMAND") && !aimsflags.KBIN_AIMS_GEOM_FILE.flag("FILE")) {
        file=aurostd::substring2string(AflowIn,"[AIMS_GEOM_FILE]COMMAND=",1,FALSE);
        aus << "00000  MESSAGE GEOM  generation from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_GEOM."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_GEOM."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR GEOM file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR GEOM file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xaims.GEOM << aurostd::file2string(file);       // load GEOM
        xaims.str=xstructure(xaims.GEOM,IOAIMS_GEOM);              // load structure
        file="rm -f ./_aflow_GEOM."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";     // remove temp  //CO20200502 - threadID
        aurostd::execute(file);                          // remove temp
      }
    }
    // GEOM DONE **************************************************
    xaims.GEOM_orig << xaims.GEOM.str();
    xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
    xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",FALSE);

    // must modify GEOM before calculating everything else
    if(Krun) Krun=(Krun && KBIN::AIMS_Modify_GEOM(xaims,AflowIn,FileMESSAGE,aflags,aimsflags));

    // CHECK for negative determinant
    if(det(xaims.str.scale*xaims.str.lattice)<0.0) {
      aus << "EEEEE  GEOM ERROR: the triple product of the basis vectors is negative                      " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  GEOM ERROR: exchange two basis vectors and adjust the atomic positions accordingly   " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }

    // some useful LDEBUG
    if(Krun && LDEBUG) {
      // STRUCTURE IS GENERATED
      FileMESSAGE <<  endl;
      FileMESSAGE <<  "******** STRUCTURE IN CARTESIAN *****************************" << endl;
      xaims.str.SetCoordinates(_COORDS_CARTESIAN_);
      FileMESSAGE <<  xaims.str << endl;
      FileMESSAGE <<  "******** STRUCTURE IN FRACTIONAL ****************************" << endl;
      xaims.str.SetCoordinates(_COORDS_FRACTIONAL_);
      FileMESSAGE <<  xaims.str << endl;
      FileMESSAGE <<  "*************************************************************" << endl;
      // xaims.str.write_klattice_flag=TRUE;
      FileMESSAGE <<  "SCALE" << endl;
      FileMESSAGE <<  xaims.str.scale << endl;
      FileMESSAGE <<  "DIRECT LATTICE (with scale)" << endl;
      FileMESSAGE <<  xaims.str.scale*(xaims.str.lattice) << endl;
      FileMESSAGE <<  "RECIPROCAL LATTICE" << endl;
      FileMESSAGE <<  (xaims.str.klattice) << endl;
      FileMESSAGE <<  "ORTOGONALITY (a*b')/2pi=I" << endl;
      FileMESSAGE <<  (xaims.str.scale*(xaims.str.lattice))*trasp(xaims.str.klattice)/(2.0*pi) << endl;
      FileMESSAGE <<  "*************************************************************" << endl;
    }
    // done produced and modified
    return Krun;
  };  // KBIN::AIMS_Produce_GEOM
}

namespace KBIN {
  bool AIMS_Produce_GEOM(_xaims& xaims) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool Krun=TRUE;
    xaims.GEOM.str(std::string());xaims.GEOM.clear();
    xaims.GEOM_orig.str(std::string());xaims.GEOM_orig.clear();
    xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",FALSE);
    xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",FALSE);
    xaims.GEOM << xaims.str;
    // GEOM done
    xaims.GEOM_orig << xaims.GEOM.str();
    xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
    return Krun;
  };  // KBIN::AIMS_Produce_GEOM
}

namespace KBIN {
  bool AIMS_Modify_GEOM(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_aimsflags &aimsflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    bool Krun=TRUE;
    if(xaims.aopts.flag("FLAG::XAIMS_GEOM_generated")==FALSE) {
      aus << "EEEEE  KBIN::AIMS_Modify_GEOM: can`t modify GEOM if it does not exist ! " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // return Krun;
    // xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",FALSE);
    //LDEBUG=TRUE;


    // CONVERT_UNIT_CELL STUFF
    // GEOM must be modified before doing the KPOINTS
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.isentry) {        /*************** GEOM **************/
      //    aus << "00000  MESSAGE GEOM  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "00000  MESSAGE-OPTION  [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;
    // GEOM must be modified before doing the KPOINTS
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_PRESERVE construction /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  PRESERVE Unit Cell " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
    // GEOM must be modified before doing the KPOINTS
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_STANDARD_PRIMITIVE construction /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  STANDARD_PRIMITIVE Unit Cell " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.Standard_Primitive_UnitCellForm();
      //CO START
      //CO, fix issue with iatoms tag becoming atom names
      bool write_inequivalent_flag=xaims.str.write_inequivalent_flag;
      //CO END
      string bravais_lattice_type=xaims.str.bravais_lattice_type,bravais_lattice_variation_type=xaims.str.bravais_lattice_variation_type,pearson_symbol=xaims.str.pearson_symbol;
      xaims.str.title=xaims.str.title+" [Standard_Primitive Unit Cell Form]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      //CO START
      //CO, fix if write_inequivalent_flag is present
      xaims.str.write_inequivalent_flag=FALSE;
      //CO, fix if write_inequivalent_flag is present
      xaims.GEOM << xaims.str;
      xaims.str.clear();xaims.GEOM >> xaims.str;  //CO, this is important, clear all symmetry stuff as the whole lattice has changed //DX20191220 - uppercase to lowercase clear
      //CO add these flags to prevent recalculation and wasted effort
      xaims.str.Standard_Lattice_calculated=TRUE;
      xaims.str.Standard_Lattice_primitive=TRUE;
      //CO add these flags to prevent recalculation and wasted effort
      //CO, fix issue with iatoms tag becoming atom names
      xaims.str.write_inequivalent_flag=write_inequivalent_flag;
      //CO END
      xaims.str.bravais_lattice_type=bravais_lattice_type;xaims.str.bravais_lattice_variation_type=bravais_lattice_variation_type;xaims.str.pearson_symbol=pearson_symbol;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
      aus << "00000  MESSAGE GEOM  STANDARD_PRIMITIVE Unit Cell Lattice = ["+xaims.str.bravais_lattice_type << "," << xaims.str.bravais_lattice_variation_type << "," << xaims.str.pearson_symbol << "]" << "  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
    // GEOM must be modified before doing the KPOINTS
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL construction /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  STANDARD_CONVENTIONAL Unit Cell " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.Standard_Conventional_UnitCellForm();
      //CO START
      //CO, fix issue with iatoms tag becoming atom names
      bool write_inequivalent_flag=xaims.str.write_inequivalent_flag;
      //CO END
      string bravais_lattice_type=xaims.str.bravais_lattice_type,bravais_lattice_variation_type=xaims.str.bravais_lattice_variation_type,pearson_symbol=xaims.str.pearson_symbol;
      xaims.str.title=xaims.str.title+" [Standard_Conventional Unit Cell Form]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      //CO START
      //CO, fix if write_inequivalent_flag is present
      xaims.str.write_inequivalent_flag=FALSE;
      //CO, fix if write_inequivalent_flag is present
      xaims.GEOM << xaims.str;
      xaims.str.clear();xaims.GEOM >> xaims.str;  //CO, this is important, clear all symmetry stuff as the whole lattice has change //DX20191220 - uppercase to lowercase clear
      //CO add these flags to prevent recalculation and wasted effort
      xaims.str.Standard_Lattice_calculated=TRUE;
      xaims.str.Standard_Lattice_conventional=TRUE;
      //CO add these flags to prevent recalculation and wasted effort
      //CO, fix issue with iatoms tag becoming atom names
      xaims.str.write_inequivalent_flag=write_inequivalent_flag;
      //CO END
      xaims.str.bravais_lattice_type=bravais_lattice_type;xaims.str.bravais_lattice_variation_type=bravais_lattice_variation_type;xaims.str.pearson_symbol=pearson_symbol;
      // xaims.str.clear();xaims.GEOM >> xaims.str; //DX20191220 - uppercase to lowercase clear
      // cout << xaims.str << endl;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
      aus << "00000  MESSAGE GEOM  STANDARD_CONVENTIONAL Unit Cell Lattice = ["+xaims.str.bravais_lattice_type << "," << xaims.str.bravais_lattice_variation_type << "," << xaims.str.pearson_symbol << "]" << "  " << Message(__AFLOW_FILE__,aflags) << endl; //CO
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //CO
      // cerr << det(xaims.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
    // GEOM must be modified before doing the KPOINTS
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_NIGGLI construction                     /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  NIGGLI Unit Cell Reduction " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.NiggliUnitCellForm();
      xaims.str.title=xaims.str.title+" [Niggli Unit Cell Form]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
      // cerr << det(xaims.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(MINKOWSKI)=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_MINKOWSKI construction               /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  MINKOWSKI Basis Reduction " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.MinkowskiBasisReduction();
      xaims.str.title=xaims.str.title+" [Minkowski Basis Reduction]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
      // cerr << det(xaims.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_INCELL construction                     /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  INCELL Unit Cell Basis " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.BringInCell();
      xaims.str.title=xaims.str.title+" [Bring In Cell Basis]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_COMPACT construction                   /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  COMPACT Unit Cell Basis " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.BringInCompact();
      xaims.str.title=xaims.str.title+" [Bring In Compact Basis]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_WIGNERSEITZ construction           /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  WIGNERSEITZ Unit Cell Basis " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.BringInWignerSeitz();
      xaims.str.title=xaims.str.title+" [WignerSeitz Basis]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_CARTESIAN construction               /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  CARTESIAN Basis Coordinates" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetCoordinates(_COORDS_CARTESIAN_);
      xaims.str.title=xaims.str.title+" [WignerSeitz Basis]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_FRACTIONAL construction            /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  FRACTIONAL Basis Coordinate" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetCoordinates(_COORDS_FRACTIONAL_);
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT")) {  // [AIMS_FORCE_OPTION]CONVERT_UNIT_CELL_DIRECT construction                    /*************** GEOM **************/
      aus << "00000  MESSAGE GEOM  DIRECT Basis Coordinate" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetCoordinates(_COORDS_FRACTIONAL_);
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag(\"EQUAL_EQUAL\")=" << aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("IMPLICIT") && aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag("EQUAL_EQUAL")) {  // [AIMS_GEOM_FILE]VOLUME=                       /*************** GEOM **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[AIMS_GEOM_FILE]VOLUME=",FALSE);
      aus << "00000  MESSAGE GEOM  IMPLICIT Volume = " << factor << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::AIMS_Modify_GEOM Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE GEOM  IMPLICIT Old Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetVolume(factor);
      aus << "00000  MESSAGE GEOM  IMPLICIT New Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.title=xaims.str.title+" [Forced Volume =]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL")) {  // [AIMS_FORCE_OPTION]VOLUME=                                                    /*************** GEOM **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[AIMS_FORCE_OPTION]VOLUME=",FALSE);
      double factor=aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedutype<double>("EQUAL_EQUAL");
      aus << "00000  MESSAGE GEOM  FORCE Volume = " << factor << " (factor1=" << factor1 << ")  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::AIMS_Modify_GEOM Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE GEOM  FORCE Old Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetVolume(factor);
      aus << "00000  MESSAGE GEOM  FORCE New Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.title=xaims.str.title+" [Forced Volume =]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("IMPLICIT") && aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag("MULTIPLY_EQUAL")) {  // [AIMS_GEOM_FILE]VOLUME*=             /*************** GEOM **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[AIMS_GEOM_FILE]VOLUME*=",FALSE);
      //     double factor=aurostd::string2utype<double>(aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.getattachedscheme("MULTIPLY_EQUAL"));
      //      cerr << "CORMAC MULTIPLY_EQUAL=" << factor << endl;
      aus << "00000  MESSAGE GEOM  IMPLICIT Volume *= " << factor << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::AIMS_Modify_GEOM Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE GEOM  IMPLICIT Old Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetVolume(xaims.str.Volume()*factor);
      aus << "00000  MESSAGE GEOM  IMPLICIT New Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.title=xaims.str.title+" [Forced Volume *=]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL")) {  // [AIMS_FORCE_OPTION]VOLUME*=                                                    /*************** GEOM **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[AIMS_FORCE_OPTION]VOLUME*=",FALSE);
      double factor=aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedutype<double>("MULTIPLY_EQUAL");
      aus << "00000  MESSAGE GEOM  FORCE Volume *= " << factor << " (factor1=" << factor1 << ")  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::AIMS_Modify_GEOM Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE GEOM  FORCE Old Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetVolume(xaims.str.Volume()*factor);
      aus << "00000  MESSAGE GEOM  FORCE New Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.title=xaims.str.title+" [Forced Volume *=]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag(PLUS_EQUAL)=" << aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag("PLUS_EQUAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_GEOM_MODE.flag("IMPLICIT") && aimsflags.KBIN_AIMS_GEOM_FILE_VOLUME.flag("PLUS_EQUAL")) {  // [AIMS_GEOM_FILE]VOLUME+=               /*************** GEOM **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[AIMS_GEOM_FILE]VOLUME+=",FALSE);
      aus << "00000  MESSAGE GEOM  IMPLICIT Volume += " << factor << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE GEOM  IMPLICIT Old Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetVolume(xaims.str.Volume()+factor);
      aus << "00000  MESSAGE GEOM  IMPLICIT New Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.title=xaims.str.title+" [Forced Volume +=]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    if(LDEBUG) cerr << "aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag(PLUS_EQUAL)=" << aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL")) {  // [AIMS_FORCE_OPTION]VOLUME+=                                                    /*************** GEOM **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[AIMS_FORCE_OPTION]VOLUME+=",FALSE);
      double factor=aimsflags.KBIN_AIMS_FORCE_OPTION_VOLUME.getattachedutype<double>("PLUS_EQUAL");
      aus << "00000  MESSAGE GEOM  FORCE Volume += " << factor << " (factor1=" << factor1 << ")  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE GEOM  FORCE Old Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.SetVolume(xaims.str.Volume()+factor);
      aus << "00000  MESSAGE GEOM  FORCE New Volume= " << xaims.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xaims.str.title=xaims.str.title+" [Forced Volume +=]";
      xaims.GEOM.str(std::string());xaims.GEOM.clear();
      xaims.GEOM << xaims.str;
      xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
      xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    }

    // GEOM done
    if(Krun && aimsflags.KBIN_AIMS_FORCE_OPTION_NOTUNE.isentry==FALSE) {
      if(0) {
        aus << "00000  MESSAGE-OPTION  XXXXX" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // ------------------------------------
    // end
    // xaims.aopts.flag("FLAG::XAIMS_GEOM_generated",TRUE);
    return Krun;
  }; // KBIN::AIMS_Produce_GEOM
}

namespace KBIN {
  bool AIMS_Reread_GEOM(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(!aurostd::FileExist(xaims.Directory+"/GEOM")) {
      aus << "EEEEE  KBIN::AIMS_Reread_GEOM: GEOM not present in directory: " << xaims.Directory << " - "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    xaims.GEOM_orig.str(std::string()); xaims.GEOM_orig << xaims.GEOM.str();
    xaims.GEOM.str(std::string()); xaims.GEOM << aurostd::file2string(xaims.Directory+"/GEOM"); // DID REREAD
    xaims.aopts.flag("FLAG::XAIMS_GEOM_changed",TRUE);
    return Krun;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// CONTROL MODIFICATIONS

// ***************************************************************************
// KBIN::XAIMS_CONTROL_PSTRESS
namespace KBIN {
  bool XAIMS_CONTROL_PREPARE_GENERIC(string command,_xaims& xaims,_aimsflags& aimsflags,string svalue,int ivalue,double dvalue,bool OPTION) {
    bool DONE=FALSE;
    bool LDEBUG=FALSE;
    string FileContent,strline;
    FileContent=xaims.CONTROL.str();
    xaims.CONTROL.str(std::string());
    xaims.aopts.flag("FLAG::XAIMS_CONTROL_changed",TRUE);

    if(aimsflags.KBIN_AIMS_CONTROL_VERBOSE) xaims.CONTROL << "# Preparing generic control.in" << endl;

    // ***************************************************************************
    // GENERIC GENERIC GENERIC GENERIC GENERIC GENERIC
    if(command=="GENERIC") {
      DONE=TRUE;
    }

    // ***************************************************************************
    // GENERIC GENERIC GENERIC GENERIC GENERIC GENERIC
    if(command=="GENERIC") {
      DONE=TRUE;
    }

    // ***************************************************************************
    // GENERIC GENERIC GENERIC GENERIC GENERIC GENERIC
    if(command=="GENERIC") {
      DONE=TRUE;
    }

    // ***************************************************************************

    if(LDEBUG) cout << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: RETURNING CORRECTLY" << endl;
    if(LDEBUG) cout << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: command=" << command << endl;
    if(LDEBUG) cout << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: svalue=" << svalue << endl;
    if(LDEBUG) cout << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: ivalue=" << ivalue << endl;
    if(LDEBUG) cout << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: dvalue=" << dvalue << endl;
    if(LDEBUG) cout << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: OPTION=" << OPTION << endl;

    if(DONE) return TRUE;

    cerr << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: ERROR" << endl;
    cerr << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: command=" << command << " not found " << endl;
    cerr << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: svalue=" << svalue << endl;
    cerr << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: ivalue=" << ivalue << endl;
    cerr << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: dvalue=" << dvalue << endl;
    cerr << XPID << "KBIN::XAIMS_CONTROL_PREPARE_GENERIC: OPTION=" << OPTION << endl;

    return FALSE; //previously an exit
  }
} 


// ***************************************************************************
// KBIN::XAIMS_CONTROL_REMOVE_ENTRY
namespace KBIN {
  void XAIMS_CONTROL_REMOVE_ENTRY(_xaims& xaims,string ENTRY,string COMMENT,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xaims.CONTROL.str();
    xaims.CONTROL.str(std::string());
    xaims.aopts.flag("FLAG::XAIMS_CONTROL_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,ENTRY,TRUE) || aurostd::substring2bool(strline,"#"+ENTRY,TRUE)) {
        if(VERBOSE) xaims.CONTROL << "# " << strline << " # AFLOW REMOVED (KBIN::XAIMS_CONTROL_REMOVE_ENTRY) " << COMMENT << endl;
      } else {
        if(!VERBOSE && strline.length()) xaims.CONTROL << strline << endl;
        if(VERBOSE) xaims.CONTROL << strline << endl;
      }
    }
  }
}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
