// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// contains routines to extract MATLAB information from _AFLOWIN_

#ifndef _AFLOW_MATLAB_CPP
#define _AFLOW_MATLAB_CPP

#include "aflow.h"

bool KBIN_MATLAB_Extract(string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
  if(AflowIn.length()) {;} // phony, just to keep AflowIn busy
  ostringstream aus;
  stringstream aflowm;
  bool Krun=TRUE;
  Krun=(Krun && (kflags.AFLOW_MATLAB_MODE_EXPLICIT || kflags.AFLOW_MATLAB_MODE_EXTERNAL));
  if(!Krun) {
    aus << "EEEEE  [AFLOW_MATLAB_MODE_EXPLICIT] and [AFLOW_MATLAB_MODE_EXTERNAL] are the only supported modes "  << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
    Krun=FALSE;
    return Krun;
  }
  // EXPLICIT
  if(Krun && kflags.AFLOW_MATLAB_MODE_EXPLICIT) {  // [MATLAB_MODE_EXPLICIT] construction
    if(kflags.AFLOW_MATLAB_FILE && !kflags.AFLOW_MATLAB_MODE_EXPLICIT_START_STOP) {
      aus << "00000  MESSAGE MATLAB  generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      aurostd::string2file(aurostd::substring2string(FileAFLOWIN,"[AFLOW_MATLAB_FILE]",0),string(aflags.Directory+"/aflow.m"));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToFileEXPLICIT(FileAFLOWIN,string(aflags.Directory+"/aflow.m"),"[AFLOW_MATLAB_FILE]");
    } else if(!kflags.AFLOW_MATLAB_FILE && kflags.AFLOW_MATLAB_MODE_EXPLICIT_START_STOP) {
      aus << "00000  MESSAGE MATLAB  generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      aurostd::string2file(aurostd::substring2string(FileAFLOWIN,"[MATLAB_MODE_EXPLICIT]START","[MATLAB_MODE_EXPLICIT]STOP",0),string(aflags.Directory+"/aflow.m"));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToFileEXPLICIT(FileAFLOWIN,string(aflags.Directory+"/aflow.m"),"[MATLAB_MODE_EXPLICIT]START","[MATLAB_MODE_EXPLICIT]STOP");
    } else {
      aus << "EEEEE  [MATLAB_MODE_EXPLICIT] do not confuse aflow !!"  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [MATLAB_MODE_EXPLICIT] Possible modes "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "----------------------------------------------------------------------------------------------------" << endl;
      aus << "[AFLOW] MATLAB EXPLICIT MODE without START/STOP (default)" << endl;
      aus << "[MATLAB_MODE_EXPLICIT]" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]&& MATLAB of the structure example with [AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]&& alpha-Boron SG 166 hR32 with two atoms in 18h" << endl;
      aus << "[AFLOW_MATLAB_FILE]&& and symmetry from the international crystallography book" << endl;
      aus << "[AFLOW_MATLAB_FILE]clear" << endl;
      aus << "[AFLOW_MATLAB_FILE]a=4.927" << endl;
      aus << "[AFLOW_MATLAB_FILE]ca=2.550;" << endl;
      aus << "[AFLOW_MATLAB_FILE]L=[a             0 0;" << endl;
      aus << "[AFLOW_MATLAB_FILE]  -a/2 a*sqrt(3)/2 0;" << endl;
      aus << "[AFLOW_MATLAB_FILE]   0             0 a*ca];" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]L=L/a;L" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]%Space group 166 R-3m" << endl;
      aus << "[AFLOW_MATLAB_FILE]ATOMS_DIRECT=[];" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]% B1" << endl;
      aus << "[AFLOW_MATLAB_FILE]x= 0.1187  % from pauling" << endl;
      aus << "[AFLOW_MATLAB_FILE]y=-0.1187  % from pauling" << endl;
      aus << "[AFLOW_MATLAB_FILE]z=-0.1088  % from pauling" << endl;
      aus << "[AFLOW_MATLAB_FILE]% position 18h" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]shift=[0 0 0];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]shift=[2/3 1/3 1/3];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]shift=[1/3 2/3 2/3];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]% B2" << endl;
      aus << "[AFLOW_MATLAB_FILE]x= 0.8035  % from pauling" << endl;
      aus << "[AFLOW_MATLAB_FILE]y= 0.1965  % from pauling" << endl;
      aus << "[AFLOW_MATLAB_FILE]z= 0.9757  % from pauling" << endl;
      aus << "[AFLOW_MATLAB_FILE]% position 18h" << endl;
      aus << "[AFLOW_MATLAB_FILE] " << endl;
      aus << "[AFLOW_MATLAB_FILE]shift=[0 0 0];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]shift=[2/3 1/3 1/3];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]shift=[1/3 2/3 2/3];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW_MATLAB_FILE]ATOMS_DIRECT" << endl;
      aus << "[AFLOW_MATLAB_FILE]" << endl;
      aus << "[AFLOW]" << endl;
      aus << "----------------------------------------------------------------------------------------------------" << endl;
      aus << "[AFLOW] MATLAB EXPLICIT MODE with START/STOP" << endl;
      aus << "[MATLAB_MODE_EXPLICIT]" << endl;
      aus << "[MATLAB_MODE_EXPLICIT]START" << endl;
      aus << "&& MATLAB of the structure example with START/STOP" << endl;
      aus << "&& alpha-Boron SG 166 hR32 with two atoms in 18h" << endl;
      aus << "&& and symmetry from the international crystallography book" << endl;
      aus << "clear" << endl;
      aus << "a=4.927" << endl;
      aus << "ca=2.550;" << endl;
      aus << "L=[a             0 0;" << endl;
      aus << "  -a/2 a*sqrt(3)/2 0;" << endl;
      aus << "   0             0 a*ca];" << endl;
      aus << "" << endl;
      aus << "L=L/a;L" << endl;
      aus << "" << endl;
      aus << "%Space group 166 R-3m" << endl;
      aus << "ATOMS_DIRECT=[];" << endl;
      aus << "" << endl;
      aus << "% B1" << endl;
      aus << "x= 0.1187  % from pauling" << endl;
      aus << "y=-0.1187  % from pauling" << endl;
      aus << "z=-0.1088  % from pauling" << endl;
      aus << "% position 18h" << endl;
      aus << "" << endl;
      aus << "shift=[0 0 0];" << endl;
      aus << "atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "shift=[2/3 1/3 1/3];" << endl;
      aus << "atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "shift=[1/3 2/3 2/3];" << endl;
      aus << "atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "% B2" << endl;
      aus << "x= 0.8035  % from pauling" << endl;
      aus << "y= 0.1965  % from pauling" << endl;
      aus << "z= 0.9757  % from pauling" << endl;
      aus << "% position 18h" << endl;
      aus << " " << endl;
      aus << "shift=[0 0 0];" << endl;
      aus << "atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "shift=[2/3 1/3 1/3];" << endl;
      aus << "atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "shift=[1/3 2/3 2/3];" << endl;
      aus << "atom=+[   x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[   x  2*x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[-2*x   -x  z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[  -x -2*x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "atom=+[ 2*x    x -z];ATOMS_DIRECT=[ATOMS_DIRECT; atom+shift];" << endl;
      aus << "" << endl;
      aus << "ATOMS_DIRECT" << endl;
      aus << "[MATLAB_MODE_EXPLICIT]STOP" << endl;
      aus << "[AFLOW]" << endl;
      aus << "----------------------------------------------------------------------------------------------------" << endl;
      aus << "EEEEE  [MATLAB_MODE_EXPLICIT] Note "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [MATLAB_MODE_EXPLICIT]START must be present and no [AFLOW_MATLAB_FILE]"  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [MATLAB_MODE_EXPLICIT]STOP  must be present and no [AFLOW_MATLAB_FILE]"  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  or [AFLOW_MATLAB_FILE] present and NO START/STOP"  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
      Krun=FALSE;
      return Krun;
    }      
  }
  // aflags.Directory+"/aflow.m"
  // EXTERNAL **************************************************
  if(Krun && kflags.AFLOW_MATLAB_MODE_EXTERNAL) {  // [AFLOW_MATLAB_MODE_EXTERNAL] construction
    string file;
    aus << "00000  MESSAGE aflow.m   generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
    if(kflags.AFLOW_MATLAB_FILE_COMMAND && kflags.AFLOW_MATLAB_FILE_FILE) {
      aus << "EEEEE   [AFLOW_MATLAB_MODE]FILE=  and  [AFLOW_MATLAB_MODE]COMMAND=  can not be used together "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
      Krun=FALSE;
      return Krun;
    }
    if(!kflags.AFLOW_MATLAB_FILE_COMMAND && (kflags.AFLOW_MATLAB_FILE_FILE || !kflags.AFLOW_MATLAB_FILE_FILE)) {
      if(kflags.AFLOW_MATLAB_FILE_FILE) {
        file=aurostd::substring2string(AflowIn,"[AFLOW_MATLAB_FILE]FILE=",1,TRUE);
        aus << "00000  MESSAGE aflow.m   generation from file=" << file << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
      } else {
        file=DEFAULT_VASP_EXTERNAL_INCAR;
        aus << "00000  MESSAGE aflow.m   generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_INCAR << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
      }
      if(!aurostd::FileExist(file)) {
        aus << "EEEEE  ERROR aflow.m file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
      if(aurostd::FileEmpty(file)) {
        aus << "EEEEE  ERROR aflow.m file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
      aflowm.clear();aflowm.str(std::string());
      aflowm << aurostd::file2string(file);
      if(Krun) Krun=(Krun && aurostd::stringstream2file(aflowm,string(aflags.Directory+"/aflow.m")));
    }
    if(kflags.AFLOW_MATLAB_FILE_COMMAND && !kflags.AFLOW_MATLAB_FILE_FILE) {
      file=aurostd::substring2string(AflowIn,"[AFLOW_MATLAB_FILE]COMMAND=",1,FALSE);
      aus << "00000  MESSAGE aflow.m   generation from command= '" << file << "' " << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
      file=file+" > ./_aflow_AFLOWM."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp //CO20200502 - threadID
      aurostd::execute(file);                           // create temp
      file="./_aflow_AFLOWM."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name //CO20200502 - threadID
      if(!aurostd::FileExist(file)) {  // could not write (directory protected)
        aus << "EEEEE  ERROR aflow.m file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
      if(aurostd::FileEmpty(file)) {  // contains nothing good
        aus << "EEEEE  ERROR aflow.m file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
        Krun=FALSE;
        return Krun;
      }
      aflowm.clear();aflowm.str(std::string());
      aflowm << aurostd::file2string(file);
      if(Krun) Krun=(Krun && aurostd::stringstream2file(aflowm,string(aflags.Directory+"/aflow.m")));
      file="rm -f ./_aflow_AFLOWM."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";     // remove temp  //CO20200502 - threadID
      aurostd::execute(file);                          // remove temp
    }
  }
  // DONE
  return Krun;
}

// ***************************************************************************
bool KBIN_MATLAB_Run(_aflags &aflags) {
  ostringstream aus;
  bool Krun=TRUE;
  if(Krun) {
    aurostd::StringstreamClean(aus);
    aus << "cd " << string(aflags.Directory) << endl;
    aus << endl;
    aus << "echo \"exit\" >> aflow.m"  << endl;
    aus << endl;
    aus << DEFAULT_KBIN_MATLAB_BIN << " -r aflow > aflow.out" << endl;
    aurostd::execute(aus);
  }
  return Krun;
  //  return FALSE;
};

// ***************************************************************************
_kflags KBIN_MATLAB_Get_Matlabflags_from_AflowIN(string &AflowIn) {
  _kflags kflags;
  // HOW TO RUN
  // INPUT FILES
  kflags.AFLOW_MATLAB_FILE                      =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_FILE]");          
  kflags.AFLOW_MATLAB_MODE_EXPLICIT             =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_MODE_EXPLICIT]");
  kflags.AFLOW_MATLAB_MODE_EXPLICIT_START_STOP   =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_MODE_EXPLICIT]START") &&
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_MODE_EXPLICIT]STOP");
  kflags.AFLOW_MATLAB_MODE_IMPLICIT             =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_MODE_IMPLICIT]");
  kflags.AFLOW_MATLAB_MODE_EXTERNAL             =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_MODE_EXTERNAL]");
  kflags.AFLOW_MATLAB_FILE_FILE       =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_FILE]FILE=",TRUE);
  kflags.AFLOW_MATLAB_FILE_COMMAND       =
    aurostd::substring2bool(AflowIn,"[AFLOW_MATLAB_FILE]COMMAND=",TRUE);

  return kflags;
}

// ***************************************************************************
bool KBIN_MATLAB_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) { // AFLOW_FUNCTION_IMPLEMENTATION
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
  //DX20190125 [OBSOLETE] - need to remove null bytes : AflowIn="";char c; while (FileAFLOWIN.get(c)) AflowIn+=c;               // READ _AFLOWIN_ and put into AflowIn
  AflowIn="";char c; while (FileAFLOWIN.get(c)) if(c!='\0'){ AflowIn+=c; }              // READ _AFLOWIN_ and put into AflowIn //DX20190125 remove null bytes
  FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
  if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
    aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(aus,XHOST.QUIET);
    return FALSE;
  }
  aflags.QUIET=FALSE;
  kflags=KBIN_MATLAB_Get_Matlabflags_from_AflowIN(AflowIn);

  // ***************************************************************************
  // Get the KBIN_BIN name
  aurostd::StringstreamClean(aus);
  aus << "00000  MESSAGE KBIN_MATLAB_Directory Running KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__,aflags) << endl;
  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
  // ***************************************************************************
  Krun=TRUE;  // guess everything is intelligent !
  // ***************************************************************************
  // READY TO RUN
  if(Krun) {
    bool Krun=true;
    ostringstream aus;
    // ***************************************************************************
    // directory check
    ifstream DirectoryStream;
    DirectoryStream.open(aflags.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      aus << "XXXXX  MAKING DIRECTORY = " << aflags.Directory << "  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+aflags.Directory;
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
      aus << "AAAAA  MATLAB RUN - " <<  aflags.Directory << " - " << kflags.KBIN_BIN << " - " << Message(__AFLOW_FILE__) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=(Krun && KBIN_MATLAB_Run(aflags));
      // ***************************************************************************
      // POSTSCRIPT
      if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
        KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
      // ***************************************************************************
    }
  }
  FileAFLOWIN.clear();FileAFLOWIN.close();
  return Krun;
}





#endif  // _AFLOW_MATLAB_MODE_CPP


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
