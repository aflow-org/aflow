// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// contains routines to extract MPI information from _AFLOWIN_


#ifndef _AFLOW_MPI_CPP
#define _AFLOW_MPI_CPP
#include "aflow.h"

#define _KQSUB_CHECK_SLEEP_  120
#define QSUB_walltime_days   16

// ***************************************************************************
// KBIN::QSUB_Extract
// ***************************************************************************
// This function extracts from _AFLOWIN_ the parameters for QSUB run
namespace KBIN {
  bool QSUB_Extract(_xqsub& xqsub,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()) {;} // phony, just to keep AflowIn busy
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_QSUB: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    xqsub.QSUB.str(std::string());xqsub.QSUB.clear();
    xqsub.QSUB_orig.str(std::string());xqsub.QSUB_orig.clear();
    xqsub.QSUB_generated=FALSE;
    xqsub.QSUB_changed=FALSE;
    //
    kflags.KBIN_QSUB=aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]");                         // search for QSUB string
    if(kflags.KBIN_QSUB) {
      aus      << "00000  [AFLOW_MODE_QSUB] found in " << _AFLOWIN_ << " "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=kflags.KBIN_QSUB;

      // INPUT FILES
      kflags.KBIN_QSUB_FILE                      =
        aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_FILE]");          
      kflags.KBIN_QSUB_MODE_EXPLICIT             =
        aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE_EXPLICIT]");
      kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP   =
        aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE_EXPLICIT]START") &&
        aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE_EXPLICIT]STOP");
      kflags.KBIN_QSUB_MODE_IMPLICIT             =
        aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE_IMPLICIT]");

      if(0) {
        cerr << " kflags.KBIN_QSUB                          =  " << kflags.KBIN_QSUB << endl;
        cerr << " kflags.KBIN_QSUB_FILE                     =  " << kflags.KBIN_QSUB_FILE << endl;
        cerr << " kflags.KBIN_QSUB_MODE_EXPLICIT            =  " << kflags.KBIN_QSUB_MODE_EXPLICIT << endl;
        cerr << " kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP =  " << kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP << endl;
        cerr << " kflags.KBIN_QSUB_MODE_IMPLICIT            =  " << kflags.KBIN_QSUB_MODE_IMPLICIT << endl;
      }

      // get (string) kflags.KBIN_QSUB_COMMAND
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE]COMMAND=",TRUE)) {
        kflags.KBIN_QSUB_COMMAND=QSUB_COMMAND_DEFAULT;
        aus << "00000  MESSAGE QSUB: COMMAND string is missing, taking COMMAND=\"" << kflags.KBIN_QSUB_COMMAND << "\"  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_QSUB_COMMAND=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_QSUB_MODE]COMMAND=",1,TRUE),'"');
        aus << "00000  MESSAGE QSUB: found COMMAND=\"" << kflags.KBIN_QSUB_COMMAND << "\"  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_QSUB_PARAMS
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE]PARAMS=",TRUE)) {
        kflags.KBIN_QSUB_PARAMS=QSUB_PARAMS_DEFAULT;
        aus << "00000  MESSAGE QSUB: PARAMS string is missing, taking PARAMS=\"" << kflags.KBIN_QSUB_PARAMS << "\"  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_QSUB_PARAMS=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_QSUB_MODE]PARAMS=",1,TRUE),'"');
        aus << "00000  MESSAGE QSUB: found PARAMS=\"" << kflags.KBIN_QSUB_PARAMS << "\"  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }

      Krun=(Krun && kflags.KBIN_QSUB_MODE_EXPLICIT);
      if(!Krun) {
        aus << "EEEEE  [AFLOW_QSUB_MODE_IMPLICIT] is the only supported mode "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
      if(Krun && kflags.KBIN_QSUB_MODE_EXPLICIT) {  // [AFLOW_QSUB_MODE_EXPLICIT] construction
        if(kflags.KBIN_QSUB_FILE && !kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP) {
          aus << "00000  MESSAGE QSUB   generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
          xqsub.QSUB.str(aurostd::substring2string(FileAFLOWIN,"[AFLOW_QSUB_FILE]"));
          //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xqsub.QSUB,"[AFLOW_QSUB_FILE]");
        } else if(!kflags.KBIN_QSUB_FILE && kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP) {
          aus << "00000  MESSAGE QSUB   generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if(aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[AFLOW_QSUB_MODE_EXPLICIT]STOP"))
          xqsub.QSUB.str(aurostd::substring2string(FileAFLOWIN,"[AFLOW_QSUB_MODE_EXPLICIT]START","[AFLOW_QSUB_MODE_EXPLICIT]STOP"));
            //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xqsub.QSUB,"[AFLOW_QSUB_MODE_EXPLICIT]START","[AFLOW_QSUB_MODE_EXPLICIT]STOP");
        } else {
          aus << "EEEEE  [AFLOW_QSUB_MODE_EXPLICIT] do not confuse aflow !!"  << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "EEEEE  [AFLOW_QSUB_MODE_EXPLICIT] Possible modes "  << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "----------------------------------------------------------------------------------------------------" << endl;
          aus << "[AFLOW] QSUB EXPLICIT MODE without START/STOP (default)" << endl;
          aus << "[AFLOW_QSUB_MODE_EXPLICIT] " << endl;
          aus << "[AFLOW_QSUB_FILE]#!/usr/local/bin/bash " << endl;
          aus << "[AFLOW_QSUB_FILE]# automatic qsub by AFLOW " << endl;
          aus << "[AFLOW_QSUB_FILE]# For Marylou4 " << endl;
          aus << "[AFLOW_QSUB_FILE]#PBS -l walltime=16:00:00:00 " << endl;
          aus << "[AFLOW_QSUB_FILE]#PBS -N pdpt.1 " << endl;
          // aus << "[AFLOW_QSUB_FILE]#PBS -m abe " << endl;
          aus << "[AFLOW_QSUB_FILE]#PBS -V " << endl;
          aus << "[AFLOW_QSUB_FILE]PROG=~/bin/vasp46s " << endl;
          aus << "[AFLOW_QSUB_FILE]WDIR=$PBS_O_WORKDIR " << endl;
          aus << "[AFLOW_QSUB_FILE]SCRATCH=~/compute/$PBS_JOBID " << endl;
          aus << "[AFLOW_QSUB_FILE]mkdir -p $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]if [ $? -ne 0 ]; then " << endl;
          aus << "[AFLOW_QSUB_FILE]    exit 1 " << endl;
          aus << "[AFLOW_QSUB_FILE]fi " << endl;
          aus << "[AFLOW_QSUB_FILE]cp -p $WDIR/POTCAR  $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]cp -p $WDIR/INCAR   $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]cp -p $WDIR/POSCAR  $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]cp -p $WDIR/KPOINTS $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]cd $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]if [ $? -ne 0 ]; then " << endl;
          aus << "[AFLOW_QSUB_FILE]    exit 1 " << endl;
          aus << "[AFLOW_QSUB_FILE]fi " << endl;
          aus << "[AFLOW_QSUB_FILE]# run " << endl;
          aus << "[AFLOW_QSUB_FILE]rm -f " << DEFAULT_VASP_OUT << " " << endl;
          aus << "[AFLOW_QSUB_FILE]$PROG >> " << DEFAULT_VASP_OUT << " " << endl;
          aus << "[AFLOW_QSUB_FILE]rm -f WAVECAR core " << endl;
          aus << "[AFLOW_QSUB_FILE]mv * $WDIR/ " << endl;
          aus << "[AFLOW_QSUB_FILE]cd $WDIR  " << endl;
          aus << "[AFLOW_QSUB_FILE]rm -r $SCRATCH " << endl;
          aus << "[AFLOW_QSUB_FILE]echo \"DONE\" > aflow.qsub.done " << endl;
          aus << "[AFLOW_QSUB_FILE]exit 0  " << endl;
          aus << "[AFLOW_QSUB_FILE]# automatic qsub by AFLOW " << endl;
          aus << "[AFLOW]" << endl;
          aus << "----------------------------------------------------------------------------------------------------" << endl;
          aus << "[AFLOW] QSUB EXPLICIT MODE with START/STOP" << endl;
          aus << "[AFLOW_QSUB_MODE_EXPLICIT]" << endl;
          aus << "[AFLOW_QSUB_MODE_EXPLICIT]START" << endl;
          aus << "#!/usr/local/bin/bash " << endl;
          aus << "# automatic qsub by AFLOW " << endl;
          aus << "# For Marylou4 " << endl;
          aus << "#PBS -l walltime=16:00:00:00 " << endl;
          aus << "#PBS -N pdpt.1 " << endl;
          // aus << "#PBS -m abe " << endl;
          aus << "#PBS -V " << endl;
          aus << "PROG=~/bin/vasp46s " << endl;
          aus << "WDIR=$PBS_O_WORKDIR " << endl;
          aus << "SCRATCH=~/compute/$PBS_JOBID " << endl;
          aus << "mkdir -p $SCRATCH " << endl;
          aus << "if [ $? -ne 0 ]; then " << endl;
          aus << "    exit 1 " << endl;
          aus << "fi " << endl;
          aus << "cp -p $WDIR/POTCAR  $SCRATCH " << endl;
          aus << "cp -p $WDIR/INCAR   $SCRATCH " << endl;
          aus << "cp -p $WDIR/POSCAR  $SCRATCH " << endl;
          aus << "cp -p $WDIR/KPOINTS $SCRATCH " << endl;
          aus << "cd $SCRATCH " << endl;
          aus << "if [ $? -ne 0 ]; then " << endl;
          aus << "    exit 1 " << endl;
          aus << "fi " << endl;
          aus << "# run " << endl;
          aus << "rm -f " << DEFAULT_VASP_OUT << " " << endl;
          aus << "$PROG >> " << DEFAULT_VASP_OUT << " " << endl;
          aus << "rm -f WAVECAR core " << endl;
          aus << "mv * $WDIR/ " << endl;
          aus << "cd $WDIR  " << endl;
          aus << "rm -r $SCRATCH " << endl;
          aus << "echo \"DONE\" > aflow.qsub.done " << endl;
          aus << "exit 0  " << endl;
          aus << "# automatic qsub by AFLOW " << endl;
          aus << "[AFLOW_QSUB_MODE_EXPLICIT]STOP" << endl;
          aus << "[AFLOW]" << endl;
          aus << "----------------------------------------------------------------------------------------------------" << endl;
          aus << "EEEEE  [AFLOW_QSUB_MODE_EXPLICIT] Note "  << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "EEEEE  [AFLOW_QSUB_MODE_EXPLICIT]START must be present and no [AFLOW_QSUB_FILE]"  << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "EEEEE  [AFLOW_QSUB_MODE_EXPLICIT]STOP  must be present and no [AFLOW_QSUB_FILE]"  << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "EEEEE  or [AFLOW_QSUB_FILE] present and NO START/STOP"  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
          Krun=FALSE;
          return Krun;
        }      
      }
    }
    // extra
    xqsub.QSUB << "echo \"DONE\" > aflow.qsub.done " << endl;
    xqsub.QSUB << "exit 0  " << endl;
    xqsub.QSUB << "# automatic qsub by AFLOW " << endl;
    // done
    xqsub.QSUB_orig << xqsub.QSUB.str();
    xqsub.QSUB_generated=TRUE;
    return Krun;
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::QSUB_RunFinished
// ***************************************************************************
// This function looks for QSUB completion (to fix)
namespace KBIN {
  bool QSUB_RunFinished(_aflags &aflags,ofstream &FileMESSAGE,bool verbose) {
    ostringstream aus_exec,aus;
    aurostd::StringstreamClean(aus_exec);
    aurostd::StringstreamClean(aus);
    // if(verbose) aus << "00000  MESSAGE RUN CHECK FINISHED : " << Message(__AFLOW_FILE__,aflags) << endl;
    // if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // aflow.qsub.done DOES NOT EXIST
    if(!aurostd::FileExist(aflags.Directory+"/aflow.qsub.done")) {
      if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (aflow.qsub.done does not exist) : " << Message(__AFLOW_FILE__,aflags) << endl;
      if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return FALSE;
    }
    // OTUCAR ESISTS BUT EMPTY
    if(aurostd::FileEmpty(aflags.Directory+"/aflow.qsub.done")) {
      if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (aflow.qsub.done is empty) : " << Message(__AFLOW_FILE__,aflags) << endl;
      if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return FALSE;
    }
    // aflow.qsub.done EXISTS
    if(aurostd::substring_present_file(aflags.Directory+"/aflow.qsub.done","DONE",TRUE) ||
        aurostd::substring_present_file(aflags.Directory+"/aflow.qsub.done","done",TRUE)) {
      if(verbose) aus << "00000  MESSAGE RUN FINISHED (aflow.qsub.done is complete) : " << Message(__AFLOW_FILE__,aflags) << endl;
      if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aurostd::execute(aus_exec);    
      return TRUE;
    }
    aurostd::execute(aus_exec);
    if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (aflow.qsub.done is incomplete) : " << Message(__AFLOW_FILE__,aflags) << endl;
    if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    return FALSE;
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::QSUB_WaitFinished
// ***************************************************************************
// This function looks for QSUB completion (to fix)
namespace KBIN {
  void QSUB_WaitFinished(_aflags &aflags,ofstream &FileMESSAGE,bool verbose) {
    while(!KBIN::QSUB_RunFinished(aflags,FileMESSAGE,verbose)) {
      aurostd::Sleep(_KQSUB_CHECK_SLEEP_);
    }
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::QSUB_Extract_Mode1
// ***************************************************************************
namespace KBIN {
  bool QSUB_Extract_Mode1(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {        // AFLOW_FUNCTION_IMPLEMENTATION  
    ostringstream aus;
    bool Krun=TRUE;
    xqsub.QSUB.str(std::string());xqsub.QSUB.clear();
    xqsub.QSUB_orig.str(std::string());xqsub.QSUB_orig.clear();
    xqsub.QSUB_generated=FALSE;
    xqsub.QSUB_changed=FALSE;
    //
    kflags.KBIN_QSUB=TRUE;                         // search for QSUB string
    aus      << "00000  [AFLOW_MODE_QSUB] QSUB Mode1 "  << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    Krun=kflags.KBIN_QSUB;

    // INPUT FILES
    kflags.KBIN_QSUB_FILE = TRUE;
    kflags.KBIN_QSUB_MODE_EXPLICIT = TRUE;
    kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP = TRUE;
    kflags.KBIN_QSUB_MODE_IMPLICIT = FALSE;

    kflags.KBIN_QSUB_COMMAND="qsub";
    kflags.KBIN_QSUB_PARAMS=" ";

    xqsub.QSUB << "#!/usr/local/bin/bash " << endl;
    xqsub.QSUB << "# automatic qsub by AFLOW " << endl;
    xqsub.QSUB << "# For Marylou4 " << endl;
    xqsub.QSUB << "# QSUB MODE 1 " << endl;  
    xqsub.QSUB << "#PBS -l walltime=" << QSUB_walltime_days << ":00:00:00 " << endl;
    xqsub.QSUB << "#PBS -N pdpt.1 " << endl;
    // xqsub.QSUB << "#PBS -m abe " << endl;
    xqsub.QSUB << "#PBS -V " << endl;
    xqsub.QSUB << "PROG=~/bin/vasp46s " << endl;
    xqsub.QSUB << "WDIR=$PBS_O_WORKDIR " << endl;
    xqsub.QSUB << "SCRATCH=~/compute/$PBS_JOBID " << endl;
    xqsub.QSUB << "mkdir -p $SCRATCH " << endl;
    xqsub.QSUB << "if [ $? -ne 0 ]; then " << endl;
    xqsub.QSUB << "    exit 1 " << endl;
    xqsub.QSUB << "fi " << endl;
    xqsub.QSUB << "cp -p $WDIR/POTCAR  $SCRATCH " << endl;
    xqsub.QSUB << "cp -p $WDIR/INCAR   $SCRATCH " << endl;
    xqsub.QSUB << "cp -p $WDIR/POSCAR  $SCRATCH " << endl;
    xqsub.QSUB << "cp -p $WDIR/KPOINTS $SCRATCH " << endl;
    xqsub.QSUB << "cd $SCRATCH " << endl;
    xqsub.QSUB << "if [ $? -ne 0 ]; then " << endl;
    xqsub.QSUB << "    exit 1 " << endl;
    xqsub.QSUB << "fi " << endl;
    xqsub.QSUB << "# run " << endl;
    xqsub.QSUB << "rm -f " << DEFAULT_VASP_OUT << " " << endl;
    xqsub.QSUB << "$PROG >> " << DEFAULT_VASP_OUT << " " << endl;
    xqsub.QSUB << "rm -f WAVECAR core " << endl;
    xqsub.QSUB << "mv * $WDIR/ " << endl;
    xqsub.QSUB << "cd $WDIR  " << endl;
    xqsub.QSUB << "rm -r $SCRATCH " << endl;
    xqsub.QSUB << "echo \"DONE\" > aflow.qsub.done " << endl;
    xqsub.QSUB << "exit 0  " << endl;
    xqsub.QSUB << "# automatic qsub by AFLOW " << endl;

    // done
    xqsub.QSUB_orig << xqsub.QSUB.str();
    xqsub.QSUB_generated=TRUE;
    return Krun;
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::QSUB_Extract_Mode2
// ***************************************************************************
namespace KBIN {
  bool QSUB_Extract_Mode2(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {        // AFLOW_FUNCTION_IMPLEMENTATION  
    ostringstream aus;
    bool Krun=TRUE;
    xqsub.QSUB.str(std::string());xqsub.QSUB.clear();
    xqsub.QSUB_orig.str(std::string());xqsub.QSUB_orig.clear();
    xqsub.QSUB_generated=FALSE;
    xqsub.QSUB_changed=FALSE;
    //
    kflags.KBIN_QSUB=TRUE;                         // search for QSUB string
    aus      << "00000  [AFLOW_MODE_QSUB] QSUB Mode2 "  << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    Krun=kflags.KBIN_QSUB;

    // INPUT FILES
    kflags.KBIN_QSUB_FILE = TRUE;
    kflags.KBIN_QSUB_MODE_EXPLICIT = TRUE;
    kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP  = TRUE;
    kflags.KBIN_QSUB_MODE_IMPLICIT = FALSE;

    kflags.KBIN_QSUB_COMMAND="qsub";
    kflags.KBIN_QSUB_PARAMS=" ";

    xqsub.QSUB << "#!/usr/local/bin/bash " << endl;
    xqsub.QSUB << "# automatic qsub by AFLOW " << endl;
    xqsub.QSUB << "# For Marylou4 " << endl;
    xqsub.QSUB << "# QSUB MODE 2 " << endl;  
    xqsub.QSUB << "#PBS -l walltime=" << QSUB_walltime_days << ":00:00:00 " << endl;
    xqsub.QSUB << "#PBS -N pdpt.1 " << endl;
    // xqsub.QSUB << "#PBS -m abe " << endl;
    xqsub.QSUB << "#PBS -V " << endl;
    xqsub.QSUB << "PROG=~/bin/vasp46s " << endl;
    xqsub.QSUB << "$PROG > " << DEFAULT_VASP_OUT << " " << endl;
    xqsub.QSUB << "echo \"DONE\" > aflow.qsub.done " << endl;
    xqsub.QSUB << "exit 0  " << endl;
    xqsub.QSUB << "# automatic qsub by AFLOW " << endl;

    // done
    xqsub.QSUB_orig << xqsub.QSUB.str();
    xqsub.QSUB_generated=TRUE;
    return Krun;
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::QSUB_Extract_Mode3
// ***************************************************************************
namespace KBIN {
  bool QSUB_Extract_Mode3(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {        // AFLOW_FUNCTION_IMPLEMENTATION  
    ostringstream aus;
    bool Krun=TRUE;
    xqsub.QSUB.str(std::string());xqsub.QSUB.clear();
    xqsub.QSUB_orig.str(std::string());xqsub.QSUB_orig.clear();
    xqsub.QSUB_generated=FALSE;
    xqsub.QSUB_changed=FALSE;
    //
    kflags.KBIN_QSUB=TRUE;                         // search for QSUB string
    aus      << "00000  [AFLOW_MODE_QSUB] QSUB Mode3 "  << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    Krun=kflags.KBIN_QSUB;

    // INPUT FILES
    kflags.KBIN_QSUB_FILE = TRUE;
    kflags.KBIN_QSUB_MODE_EXPLICIT = TRUE;
    kflags.KBIN_QSUB_MODE_EXPLICIT_START_STOP  = TRUE;
    kflags.KBIN_QSUB_MODE_IMPLICIT = FALSE;

    kflags.KBIN_QSUB_COMMAND="qsub";
    kflags.KBIN_QSUB_PARAMS=" ";

    xqsub.QSUB << "#!/usr/local/bin/bash " << endl;
    xqsub.QSUB << "# QSUB MODE 3 " << endl;  
    xqsub.QSUB << "# automatic qsub by AFLOW " << endl;

    // done
    xqsub.QSUB_orig << xqsub.QSUB.str();
    xqsub.QSUB_generated=TRUE;
    return Krun;
  }
}  // namespace KBIN

#endif  // _QSUB_IMPLEMENTATIONS_


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************

//   bool   KBIN_QSUB;
//   string KBIN_QSUB_COMMAND;
//   string KBIN_QSUB_PARAMS;
//   bool   KBIN_QSUB_MODE_EXPLICIT;
//   bool   KBIN_QSUB_MODE_EXPLICIT_START;
//   bool   KBIN_QSUB_MODE_EXPLICIT_STOP;
//   bool   KBIN_QSUB_MODE_IMPLICIT;
