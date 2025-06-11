// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AEL_GET_STRESS_CPP
#define _AFLOW_AEL_GET_STRESS_CPP
#include "aflow.h"
#include "aflow_ael_elasticity.h"

// ###############################################################################
//                  AFLOW Automatic Elasticity Library (AEL) (2014-2021)
// ###############################################################################
//
// Uses strain-stress calculations to obtain elastic constants of materials
// Based on original Python program written by M. de Jong et al.
// See Scientific Data 2, 150009 (2015) for details of original program
// See Phys. Rev. Materials 1, 015401 (2017) for details of this implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AEL
//

// *****************************************************************************************************************
// The following functions are for setting up AEL inputs for postprocessing runs called from other parts of AFLOW
// *****************************************************************************************************************
namespace AEL_functions {
  uint AEL_xvasp_flags_populate(_xvasp& xvasp, string& AflowIn, string& AflowInName, string& FileLockName, const string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) {
    ifstream FileAFLOWIN, FileAFLOWINcheck;
    string FileNameAFLOWIN = "", FileNameAFLOWINcheck = "", AflowInCheck = "";
    string FileNameMessage = "";
    ostringstream aus;
    vector<string> vAflowInCheck;
    bool ael_aflowin_found = false;
    bool Krun = true;
    bool load_POSCAR_from_xvasp = false;
    aurostd::xoption USER_AEL_POISSON_RATIO;
    USER_AEL_POISSON_RATIO.option = false;
    // Set aflags
    aflags.Directory=directory_LIB;
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
    aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
    aus << _AELSTR_MESSAGE_ << "AflowInName = " << AflowInName << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
    if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;
    aflags.KBIN_RUN_AFLOWIN=TRUE;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN=FALSE;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP=FALSE;
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN=FALSE;
    aflags.KBIN_DELETE_AFLOWIN=FALSE;
    //CT20200721 The following code has been moved to the function AEL_Get_AflowInName
    // Set FileMESSAGE name
    // [OBSOLETE] if(FileLockName.length() > 0) { //CT20200624 Moved down: only rename LOCK files if it is confirmed that this is an AEL main directory and not an ARUN.AEL directory
    // [OBSOLETE]   if (aurostd::FileExist(directory_LIB+"/"+FileLockName)) {
    // [OBSOLETE]     aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+" ")+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+".run"));   
    // [OBSOLETE]   }
    // [OBSOLETE]   string FileNameMessage=directory_LIB+"/"+FileLockName;
    // [OBSOLETE]   FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
    // [OBSOLETE] } else {
    // [OBSOLETE]   if (aurostd::FileExist(directory_LIB+"/ael.LOCK")) {
    // [OBSOLETE]     aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/ael.LOCK ")+aurostd::CleanFileName(directory_LIB+"/ael.LOCK.run"));
    // [OBSOLETE]     string FileNameMessage=directory_LIB+"/ael.LOCK";
    // [OBSOLETE]     FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
    // [OBSOLETE]   } else if (aurostd::FileExist(directory_LIB+"/agl.LOCK")) {
    // [OBSOLETE]     aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/agl.LOCK ")+aurostd::CleanFileName(directory_LIB+"/agl.LOCK.run"));
    // [OBSOLETE]     string FileNameMessage=directory_LIB+"/agl.LOCK";
    // [OBSOLETE]     FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
    // [OBSOLETE]   } else {
    // [OBSOLETE]     string FileNameMessage=directory_LIB+"/ael.LOCK";
    // [OBSOLETE]     FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
    // [OBSOLETE]   }
    // [OBSOLETE] }
    //CO20200502 START - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future
    //CT20200720 Moved finding aflow.in file to separate function
    // [OBSOLETE] vector<string> vaflowins;
    // [OBSOLETE] if(AflowInName.length()>0){vaflowins.push_back(AflowInName);}
    // [OBSOLETE] if(_AFLOWIN_.length()>0){vaflowins.push_back(_AFLOWIN_);}
    // [OBSOLETE] vaflowins.push_back(_AFLOWIN_AEL_DEFAULT_);
    // [OBSOLETE] vaflowins.push_back(_AFLOWIN_AGL_DEFAULT_);
    // [OBSOLETE] for(uint iaf=0;iaf<vaflowins.size()&&!ael_aflowin_found;iaf++){
    // [OBSOLETE]   const string& aflowinname = vaflowins.at(iaf);
    // [OBSOLETE]   if((!ael_aflowin_found) && (aurostd::FileExist(directory_LIB+"/"+aflowinname))) {
    // [OBSOLETE]     FileNameAFLOWINcheck = directory_LIB+"/"+aflowinname;
    // [OBSOLETE]     FileAFLOWINcheck.open(FileNameAFLOWINcheck.c_str(),std::ios::in);
    // [OBSOLETE]     FileAFLOWINcheck.clear();
    // [OBSOLETE]     FileAFLOWINcheck.seekg(0);
    // [OBSOLETE]     AflowInCheck="";
    // [OBSOLETE]     char c;
    // [OBSOLETE]     // READ aflowinname and put into AflowInCheck
    // [OBSOLETE]     while (FileAFLOWINcheck.get(c)) {
    // [OBSOLETE]       AflowInCheck+=c;
    // [OBSOLETE]     }
    // [OBSOLETE]     FileAFLOWINcheck.clear();
    // [OBSOLETE]     FileAFLOWINcheck.seekg(0);
    // [OBSOLETE]     AflowInCheck=aurostd::RemoveComments(AflowInCheck); // NOW Clean AFLOWIN
    // [OBSOLETE]     vAflowInCheck.clear();
    // [OBSOLETE]     aurostd::string2vectorstring(AflowInCheck,vAflowInCheck); 
    // [OBSOLETE]     // Check if aflowinname contains command to run AEL
    // [OBSOLETE]     for(uint i=0;i<vAflowInCheck.size()&&!ael_aflowin_found;i++){
    // [OBSOLETE]       if((aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AEL]CALC",TRUE) || aurostd::substring2bool(AflowInCheck,"[VASP_AEL]CALC",TRUE)) &&
    // [OBSOLETE]           !(aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AEL]CALC_",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AEL]CALC_",TRUE) ||
    // [OBSOLETE]             aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AEL]CALCS",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AEL]CALCS",TRUE) || FALSE)){
    // [OBSOLETE]         FileNameAFLOWIN = FileNameAFLOWINcheck;
    // [OBSOLETE]         ael_aflowin_found = true;
    // [OBSOLETE]       } else if((aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowInCheck,"[VASP_AGL]CALC",TRUE)) &&
    // [OBSOLETE]           !(aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC_",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALC_",TRUE) ||
    // [OBSOLETE]             aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALCS",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALCS",TRUE) || FALSE)){
    // [OBSOLETE]         if (aurostd::substring2bool(AflowInCheck,"[AFLOW_AGL]AEL_POISSON_RATIO=",TRUE) ) {
    // [OBSOLETE]           USER_AEL_POISSON_RATIO.options2entry(AflowInCheck,"[AFLOW_AGL]AEL_POISSON_RATIO=",USER_AEL_POISSON_RATIO.option);
    // [OBSOLETE]         } else if( aurostd::substring2bool(AflowInCheck,"[AFLOW_AGL]AELPOISSONRATIO=",TRUE) ) {
    // [OBSOLETE]           USER_AEL_POISSON_RATIO.options2entry(AflowInCheck,"[AFLOW_AGL]AELPOISSONRATIO=",USER_AEL_POISSON_RATIO.option);
    // [OBSOLETE]         }        	
    // [OBSOLETE]         if (USER_AEL_POISSON_RATIO.option) {
    // [OBSOLETE]           FileNameAFLOWIN = FileNameAFLOWINcheck;
    // [OBSOLETE]           ael_aflowin_found = true;
    // [OBSOLETE]         }
    // [OBSOLETE]       }
    // [OBSOLETE]     }
    // [OBSOLETE]   }
    // [OBSOLETE] }

    //CT20200721 Calls function AEL_Get_AflowInName to find correct aflow.in filename
    AEL_functions::AEL_Get_AflowInName(AflowInName, directory_LIB, ael_aflowin_found);
    FileNameAFLOWIN = directory_LIB+"/"+AflowInName;

    //CO20200502 STOP - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future

    if (ael_aflowin_found) {
      // Set FileMESSAGE name
      //CT20200624 Moved down so LOCK file is only renamed if this is actually an AEL main directory and not an ARUN.AEL directory
      if(FileLockName.length() > 0) {
        if (aurostd::FileExist(directory_LIB+"/"+FileLockName)) {
          aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+" ")+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+".run"));   
        }
        string FileNameMessage=directory_LIB+"/"+FileLockName;
        FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
      } else {
        if (aurostd::FileExist(directory_LIB+"/ael.LOCK")) {
          aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/ael.LOCK ")+aurostd::CleanFileName(directory_LIB+"/ael.LOCK.run"));
          string FileNameMessage=directory_LIB+"/ael.LOCK";
          FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
        } else if (aurostd::FileExist(directory_LIB+"/agl.LOCK")) {
          aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/agl.LOCK ")+aurostd::CleanFileName(directory_LIB+"/agl.LOCK.run"));
          string FileNameMessage=directory_LIB+"/agl.LOCK";
          FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
        } else {
          string FileNameMessage=directory_LIB+"/ael.LOCK";
          FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
        }
      }
      // Search for AEL aflow.in filename      
      FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
      FileAFLOWIN.clear();
      FileAFLOWIN.seekg(0);
      AflowIn="";
      char c;
      // READ _AFLOWIN_ and put into AflowInCheck
      while (FileAFLOWIN.get(c)) {
        AflowIn+=c;
      }
      FileAFLOWIN.clear();
      FileAFLOWIN.seekg(0);
      AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
      vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn); 
      // Set kflags
      kflags.KBIN_MPI=aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI]");
      kflags.AFLOW_MODE_VASP=aurostd::substring2bool(AflowIn,"[AFLOW_MODE=VASP]") || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_VASP]") || aurostd::substring2bool(AflowIn,"[AFLOW_MODE]VASP");                 // check VASP
      if(kflags.AFLOW_MODE_VASP && !aflags.KBIN_GEN_VASP_FROM_AFLOWIN){aflags.KBIN_GEN_VASP_FROM_AFLOWIN=true;} //do vasp last, default
      kflags.KBIN_SYMMETRY_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]CALC",TRUE);
      kflags.KBIN_SYMMETRY_NO_SCAN  = aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]NO_SCAN",TRUE);
      if(aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]SYM_EPS=",TRUE)){
        kflags.KBIN_SYMMETRY_EPS = aurostd::substring2utype<double>(AflowIn,"[AFLOW_SYMMETRY]SYM_EPS=",TRUE);
      }
      // parameters for zip/compression
      kflags.KZIP_COMPRESS=TRUE;
      aurostd::StringstreamClean(aus);
      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=none]") ||
          aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=NONE]") ||
          !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP")) {
        kflags.KZIP_COMPRESS=FALSE;
        for(int i=0;i<1;i++) {
          aus << "WWWWW  Warning no compression of output files... " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintWarningStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      } else {
        if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP")) { // "[AFLOW_MODE_ZIP=" not found
          kflags.KZIP_BIN=DEFAULT_KZIP_BIN;  // take default
          aus << "00000  MESSAGE Taking DEFAULT KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP]")) { // "[AFLOW_MODE_ZIP]" not found
          kflags.KZIP_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_ZIP]");
          aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=")) { // "[AFLOW_MODE_ZIP=" found
          kflags.KZIP_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_ZIP="),']');
          aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      // parameters for AAPL
      aurostd::xoption KBIN_PHONONS_CALCULATION_AAPL;
      KBIN_PHONONS_CALCULATION_AAPL.option=false;
      KBIN_PHONONS_CALCULATION_AAPL.options2entry(AflowIn, string("[AFLOW_AAPL]KAPPA=|[AFLOW_PHONONS]KAPPA="), KBIN_PHONONS_CALCULATION_AAPL.option, KBIN_PHONONS_CALCULATION_AAPL.xscheme); KBIN_PHONONS_CALCULATION_AAPL.option |= aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AAPL]CALC",TRUE);  //legacy
      kflags.KBIN_PHONONS_CALCULATION_AAPL  = KBIN_PHONONS_CALCULATION_AAPL.option;
      // Parameters for QHA-APL
      if (!kflags.KBIN_PHONONS_CALCULATION_AAPL) {
        kflags.KBIN_PHONONS_CALCULATION_QHA  = aurostd::substring2bool(AflowIn,"[AFLOW_QHA]CALC",TRUE) || aurostd::substring2bool(AflowIn,"VASP_QHA]CALC",TRUE);
      }
      // parameters for APL
      if(!(kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA)){ //mutually exclusive
        kflags.KBIN_PHONONS_CALCULATION_APL  = aurostd::substring2bool(AflowIn,"[AFLOW_APL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_PHONONS]CALC",TRUE);
      }
      // parameters for AGL (Debye Model)
      for(uint i=0;i<vAflowIn.size()&&!kflags.KBIN_PHONONS_CALCULATION_AGL;i++){
        if((aurostd::substring2bool(vAflowIn[i],"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AGL]CALC",TRUE)) 
            && !(aurostd::substring2bool(vAflowIn[i],"[AFLOW_AGL]CALC_",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AGL]CALC_",TRUE) ||
              aurostd::substring2bool(vAflowIn[i],"[AFLOW_AGL]CALCS",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AGL]CALCS",TRUE) || FALSE)){
          kflags.KBIN_PHONONS_CALCULATION_AGL=true;
        }
      }
      // parameters for AEL (Elastic constants)
      for(uint i=0;i<vAflowIn.size()&&!kflags.KBIN_PHONONS_CALCULATION_AEL;i++){
        if((aurostd::substring2bool(vAflowIn[i],"[AFLOW_AEL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AEL]CALC",TRUE)) 
            && !(aurostd::substring2bool(vAflowIn[i],"[AFLOW_AEL]CALC_",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AEL]CALC_",TRUE) ||
              aurostd::substring2bool(vAflowIn[i],"[AFLOW_AEL]CALCS",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AEL]CALCS",TRUE) || FALSE)){
          kflags.KBIN_PHONONS_CALCULATION_AEL=true;
        }
      }
      // parameters for POCC CALCULATIONS
      kflags.KBIN_POCC=FALSE;
      kflags.KBIN_POCC_CALCULATION=aurostd::substring2bool(AflowIn,"[AFLOW_POCC]CALC",TRUE) && (aurostd::substring2bool(AflowIn,"[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE",TRUE) && aurostd::substring2bool(AflowIn,"[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE",TRUE)); //CO20180419
      if(kflags.KBIN_POCC_CALCULATION) {
        kflags.KBIN_POCC=TRUE;
      } 
      // parameters for FROZSL
      kflags.KBIN_FROZSL=FALSE;
      kflags.KBIN_PHONONS_CALCULATION_FROZSL = aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]CALC",TRUE);
      kflags.KBIN_FROZSL_DOWNLOAD=(aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]DOWN",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]DOWNLOAD",TRUE));
      kflags.KBIN_FROZSL_FILE = aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]FILE",TRUE); 
      if(kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_FROZSL_DOWNLOAD|| kflags.KBIN_FROZSL_FILE) kflags.KBIN_FROZSL=TRUE;
      // Set KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_VASP_BIN;
      KBIN::MPI_Extract(AflowIn, FileMESSAGE, aflags, kflags);
      // Set vflags from AflowIN
      vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, FileMESSAGE, aflags, kflags);

      // Set-up xvasp
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      xvasp.clear();
      uint ixvasp=0;
      xvasp.POSCAR_index=ixvasp;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      KBIN::readModulesFromAflowIn(AflowIn, kflags, xvasp);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      xvasp.Directory=aflags.Directory;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags,load_POSCAR_from_xvasp));
      if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);      
      // Fix blank species
      if(xvasp.str.species.size()>0) {
        if(xvasp.str.species.at(0)=="") {
          pflow::fixEmptyAtomNames(xvasp.str);  
        }
      }
      if (Krun) {
        return 0;
      } else {
        return 1;
      }
    } else {
      // [OBSOLETE] aurostd::StringstreamClean(aus);
      cerr << _AELSTR_MESSAGE_ << "AEL input file not found!" << endl; //CT20200624 If no AEL file present, then this is not an AEL main directory: write to cerr and return 2 to inform calling function
      cerr << _AELSTR_MESSAGE_ << "Not an AEL main directory" << endl;      
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    }     
  }
}

// *****************************************************************************************************************
// Finds aflow.in file for AEL calculation (if it exists) //CT20200715
// *****************************************************************************************************************
namespace AEL_functions {
  uint AEL_Get_AflowInName(string& AflowInName, const string& directory_LIB, bool& ael_aflowin_found) {
    ifstream FileAFLOWIN, FileAFLOWINcheck;
    string FileNameAFLOWIN = "", FileNameAFLOWINcheck = "", AflowInCheck = "";
    string FileNameMessage = "";
    string aflowinname = "";
    string stmp="";
    uint aelerror = 0, filelength = 0;
    vector<string> vAflowInCheck;
    ael_aflowin_found = false;
    aurostd::xoption USER_AEL_POISSON_RATIO;
    USER_AEL_POISSON_RATIO.option = false;
    vector<string> vaflowins;
    if(AflowInName.length()>0){vaflowins.push_back(AflowInName);} // Check if AflowInName exists
    if(_AFLOWIN_.length()>0){vaflowins.push_back(_AFLOWIN_);} // Otherwise, check if _AFLOWIN_ file is AEL input file
    vaflowins.push_back(_AFLOWIN_AEL_DEFAULT_); // Otherwise, check for other commonly used names for AEL aflow.in file
    vaflowins.push_back(_AFLOWIN_AGL_DEFAULT_);
    for(uint iaf=0;iaf<vaflowins.size()&&!ael_aflowin_found;iaf++){
      aflowinname = vaflowins.at(iaf);
      if(aurostd::EFileExist(directory_LIB+"/"+aflowinname,stmp)&&aurostd::IsCompressed(stmp)){aurostd::UncompressFile(stmp);}  //CO20210204 - fix aflow.in.xz
      if((!ael_aflowin_found) && (aurostd::FileExist(directory_LIB+"/"+aflowinname))) {
        FileNameAFLOWINcheck = directory_LIB+"/"+aflowinname;
        // [OBSOLETE] FileAFLOWINcheck.open(FileNameAFLOWINcheck.c_str(),std::ios::in);
        // [OBSOLETE] FileAFLOWINcheck.clear();
        // [OBSOLETE] FileAFLOWINcheck.seekg(0);
        AflowInCheck="";
        // [OBSOLETE] char c;
        // READ aflowinname and put into AflowInCheck
        // [OBSOLETE] while (FileAFLOWINcheck.get(c)) {
        // [OBSOLETE]   AflowInCheck+=c;
        // [OBSOLETE] }
        // [OBSOLETE] FileAFLOWINcheck.clear();
        // [OBSOLETE] FileAFLOWINcheck.seekg(0);
        filelength=aurostd::file2string(FileNameAFLOWINcheck, AflowInCheck);
        if (filelength > 0) {
          aelerror = 0;
        } else {
          aelerror = 1;
        }	
        AflowInCheck=aurostd::RemoveComments(AflowInCheck); // NOW Clean AFLOWIN
        vAflowInCheck.clear();
        aurostd::string2vectorstring(AflowInCheck,vAflowInCheck); 
        // Check if aflowinname contains command to run AEL
        for(uint i=0;i<vAflowInCheck.size()&&!ael_aflowin_found;i++){
          if((aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AEL]CALC",TRUE) || aurostd::substring2bool(AflowInCheck,"[VASP_AEL]CALC",TRUE)) &&
              !(aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AEL]CALC_",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AEL]CALC_",TRUE) ||
                aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AEL]CALCS",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AEL]CALCS",TRUE) || FALSE)){
            FileNameAFLOWIN = FileNameAFLOWINcheck;
            ael_aflowin_found = true;
            AflowInName = aflowinname;
          } else if((aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowInCheck,"[VASP_AGL]CALC",TRUE)) &&
              !(aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC_",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALC_",TRUE) ||
                aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALCS",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALCS",TRUE) || FALSE)){
            if (aurostd::substring2bool(AflowInCheck,"[AFLOW_AGL]AEL_POISSON_RATIO=",TRUE) ) {
              USER_AEL_POISSON_RATIO.options2entry(AflowInCheck,"[AFLOW_AGL]AEL_POISSON_RATIO=",USER_AEL_POISSON_RATIO.option);
            } else if( aurostd::substring2bool(AflowInCheck,"[AFLOW_AGL]AELPOISSONRATIO=",TRUE) ) {
              USER_AEL_POISSON_RATIO.options2entry(AflowInCheck,"[AFLOW_AGL]AELPOISSONRATIO=",USER_AEL_POISSON_RATIO.option);
            }        	
            if (USER_AEL_POISSON_RATIO.option) {
              FileNameAFLOWIN = FileNameAFLOWINcheck;
              ael_aflowin_found = true;
              AflowInName = aflowinname;
            }
          }
        }
        FileAFLOWINcheck.close();
      }
    }
    return aelerror;
  }
}

// *******************************************************************************
// The following functions are for generating _AFLOWIN_ files
// *******************************************************************************

// ***************************************************************************
// AEL_functions::aelvaspflags
// ***************************************************************************
namespace AEL_functions {
  //
  // Function to assign values for VASP input flags from aflow.in file to vaspRun _xvasp class
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint aelvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, string& dirrunname, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
    ostringstream aus;
    vector<string> vfile;
    string vfilename;
    bool vfileexist = false;
    if(AEL_data.relax_static || AEL_data.static_only) {
      aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"),vfile,",");
      for(uint ij=0;ij<vfile.size();ij++) {
        if(aurostd::FileExist(dirrunname+"/"+vfile.at(ij))) {
          vfilename = vfile.at(ij);
          vfileexist = true;
        }    
      }  
    } else {
      aurostd::string2tokens(string("OUTCAR.relax2.bz2,OUTCAR.relax2.gz,OUTCAR.relax2.xz,OUTCAR.relax2"),vfile,",");
      for(uint ij=0;ij<vfile.size();ij++) {
        if(aurostd::FileExist(dirrunname+"/"+vfile.at(ij))) {
          vfilename = vfile.at(ij);
          vfileexist = true;
        }    
      }  
    }
    // SOME WARNINGS: check existence of LOCK and OUTCAR.relax2 files
    // [OBSOLETE] if( !aurostd::FileExist( vaspRun.Directory + string("/LOCK") ) &&
    if( !(aurostd::FileExist( dirrunname + "/" + _AFLOWLOCK_ ) ||
          ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/ael.LOCK") || aurostd::FileExist(dirrunname + "/LOCK")))) &&
        ( vfileexist ) ) {
      //[OBSOLETE]	aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.relax2") ) )
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AELSTR_WARNING_ + "found OUTCAR.static but no LOCK in " <<  vaspRun.Directory << endl;
      aus << _AELSTR_WARNING_ + "found " << vfilename << " but no " << _AFLOWLOCK_ << " in " <<  dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    // [OBSOLETE] if( aurostd::FileExist( vaspRun.Directory + string("/LOCK") ) &&
    if( (aurostd::FileExist( dirrunname + "/" + _AFLOWLOCK_ ) ||
          ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/ael.LOCK") || aurostd::FileExist(dirrunname + "/LOCK")))) && !(vfileexist) ) {
      //[OBSOLETE]	!aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.relax2") ) )
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AELSTR_WARNING_ + "found LOCK but no OUTCAR.static in " <<  vaspRun.Directory << endl;
      aus << _AELSTR_WARNING_ + "found " << _AFLOWLOCK_ << " but no OUTCAR in " <<  dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    // Switch off autotune
    kbinFlags.KBIN_MPI_AUTOTUNE = true;

    if(AEL_data.relax_static) {
      vaspRun.AVASP_flag_RUN_RELAX_STATIC = true;
      vaspRun.AVASP_flag_GENERATE = false;
      vaspRun.AVASP_flag_RUN_RELAX = false;
      vaspRun.AVASP_flag_RUN_STATIC = false;
      vaspRun.AVASP_flag_RUN_STATIC_BANDS = false;
      vaspRun.AVASP_flag_RUN_RELAX_STATIC_BANDS = false;
      vaspFlags.KBIN_VASP_RUN.flag("RELAX_STATIC", true);
      vaspFlags.KBIN_VASP_RUN.flag("GENERATE", false);
      vaspFlags.KBIN_VASP_RUN.flag("RELAX",false);
      vaspFlags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",false);
      vaspFlags.KBIN_VASP_RUN.flag("STATIC",false);
      vaspFlags.KBIN_VASP_RUN.flag("STATIC_BANDS",false);
      vaspFlags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry = true;      
      vaspFlags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.content_string = "IONS";
      vaspRun.aopts.flag("FLAG::VOLUME_PRESERVED",TRUE);
    } else {
      vaspRun.AVASP_flag_RUN_STATIC = true;
      vaspRun.AVASP_flag_RUN_RELAX_STATIC = false;
      vaspRun.AVASP_flag_GENERATE = false;
      vaspRun.AVASP_flag_RUN_RELAX = false;
      vaspRun.AVASP_flag_RUN_STATIC_BANDS = false;
      vaspRun.AVASP_flag_RUN_RELAX_STATIC_BANDS = false;
      vaspFlags.KBIN_VASP_RUN.flag("STATIC", true);
      vaspFlags.KBIN_VASP_RUN.flag("RELAX_STATIC", false);
      vaspFlags.KBIN_VASP_RUN.flag("GENERATE", false);
      vaspFlags.KBIN_VASP_RUN.flag("RELAX", false);
      vaspFlags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS", false);
      vaspFlags.KBIN_VASP_RUN.flag("STATIC_BANDS", false);      
    }

    // Set unit cell conversion to "PRESERVE"
    setPreserveUnitCell(vaspRun);
    vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry = true;
    vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.pop("STANDARD_PRIMITIVE");    
    vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.push("PRESERVE");

    if(AEL_data.precaccalgonorm) {
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.isentry = true;      
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.content_string = "ACCURATE";
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.isentry = true;      
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.content_string = "NORMAL";
    }

    // Switch off VASP symmetry - this can help when applied strains break the symmetry of the primitive cell
    if(!AEL_data.vasp_symmetry) {
      vaspFlags.KBIN_VASP_FORCE_OPTION_SYM.option = false;
    }

    // Common KPOINTS settings and OVERRIDES
    // [OBSOLETE] vaspRun.AVASP_KSCHEME = _vaspFlags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
    // [OBSOLETE] vaspRun.AVASP_value_KPPRA = _vaspFlags.KBIN_VASP_KPOINTS_KPPRA.content_int;
    // [OBSOLETE] vaspRun.AVASP_STATIC_KSCHEME = _vaspFlags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;
    // [OBSOLETE] vaspRun.AVASP_value_KPPRA_STATIC = _vaspFlags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_int;


    // Clear old INCAR and set it as we want...
    // Want to relax ions only keeping cell size and shape fixed
    // Might want to create new relaxation type instead of creating INCAR by hand
    // [OBSOLETE] vaspRun.INCAR.str(std::string());
    // [OBSOLETE] string system;
    // [OBSOLETE] for(uint j=0; j < vaspRun.str.species.size(); j++)
    // [OBSOLETE]   system = system + vaspRun.str.species_pp.at(j);
    // [OBSOLETE] system = system + "@" + runname;
    // [OBSOLETE] vaspRun.INCAR << "SYSTEM=" << system << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "# Added by [AFLOW_AEL] begin" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "ISIF=2           # To calculate stress tensor including ion relaxation only" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "IBRION=2         # Ion relaxation using conjugate gradient" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "NSW=51           # Relax ions for long" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "# Added by [AFLOW_AEL] end" << std::endl;

    // Change format of POSCAR
    if( ( !kbinFlags.KBIN_MPI && ( kbinFlags.KBIN_BIN.find("46") != string::npos ) ) ||
        (  kbinFlags.KBIN_MPI && ( kbinFlags.KBIN_MPI_BIN.find("46") != string::npos ) ) ) {
      vaspRun.str.is_vasp5_poscar_format = false; 
    }

    return 0;
  }
} // namespace AEL_functions

//[CT20200502 - OBSOLETE]  // ***************************************************************************
//[CT20200502 - OBSOLETE]  // AEL_functions::createAFLOWIN
//[CT20200502 - OBSOLETE]  // ***************************************************************************
//[CT20200502 - OBSOLETE]  namespace AEL_functions {
//[CT20200502 - OBSOLETE]  //
//[CT20200502 - OBSOLETE]  // Create aflow.in file: makes new directory and writes aflow.in for strained structure file inside it 
//[CT20200502 - OBSOLETE]  // Adapted from that in AFLOW APL function PhononCalculator::createAFLOWIN()
//[CT20200502 - OBSOLETE]  //
//[CT20200502 - OBSOLETE]  uint createAFLOWIN(_xvasp& vaspRun, _xvasp& xvasp, _kflags& _kbinFlags, _vflags& _vaspFlags, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
//[CT20200502 - OBSOLETE]  bool AFLOWIN_QE_FLAG=FALSE;
//[CT20200502 - OBSOLETE]  bool SPACES=FALSE;
//[CT20200502 - OBSOLETE]  ostringstream aus;
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if( !aurostd::FileExist( vaspRun.Directory) ) {
//[CT20200502 - OBSOLETE]  aurostd::DirectoryMake( vaspRun.Directory );
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  // CHMOD Directory 777: change directory permissions to read+write+execute for all users
//[CT20200502 - OBSOLETE]  aurostd::DirectoryChmod("777", vaspRun.Directory);
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Create file
//[CT20200502 - OBSOLETE]  // [OBSOLETE] string filename =  vaspRun.Directory + string("/aflow.in");
//[CT20200502 - OBSOLETE]  string filename =  vaspRun.Directory + "/" + _AFLOWIN_;
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Check if aflow.in file exists in the directory     
//[CT20200502 - OBSOLETE]  // If the aflow.in file does exist and the overwrite has not been enabled, exit createAFLOWIN for this structure
//[CT20200502 - OBSOLETE]  if( aurostd::FileExist( filename) && (!AEL_data.aflowin_overwrite) ) {
//[CT20200502 - OBSOLETE]  aurostd::StringstreamClean(aus);
//[CT20200502 - OBSOLETE]  aus << _AELSTR_MESSAGE_ + "Not overwriting existing file " << _AFLOWIN_ << endl;  
//[CT20200502 - OBSOLETE]  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]  return 0;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  ofstream outfile(filename.c_str(),ios_base::out);
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Check aflow.in file is open
//[CT20200502 - OBSOLETE]  if( !outfile.is_open() ) {
//[CT20200502 - OBSOLETE]  aurostd::StringstreamClean(aus);
//[CT20200502 - OBSOLETE]  // [OBSOLETE] aus << _AELSTR_WARNING_ + "Cannot create [aflow.in] file" << endl;
//[CT20200502 - OBSOLETE]  aus << _AELSTR_WARNING_ + "Cannot create [" << _AFLOWIN_ << "] file" << endl;
//[CT20200502 - OBSOLETE]  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]  return 1;
//[CT20200502 - OBSOLETE]  } 
//[CT20200502 - OBSOLETE]  // CHMOD a+rw aflow.in: change permissions on aflow.in file
//[CT20200502 - OBSOLETE]  aurostd::ChmodFile("a+rw",filename);
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write to aflow.in file
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW]  _  ___ _" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] / \\|   | \\ |" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] | o |-- | " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] |_n_|__ |___| automatic generated file" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW]" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE=VASP]" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  //CO20180130 START
//[CT20200502 - OBSOLETE]  //adding aflow.rc stuff
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE_BINARY=";
//[CT20200502 - OBSOLETE]  if(!_kbinFlags.KBIN_BIN.empty()){outfile << _kbinFlags.KBIN_BIN;}
//[CT20200502 - OBSOLETE]  else {outfile << DEFAULT_VASP_BIN;}
//[CT20200502 - OBSOLETE]  outfile << "]" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(!(_kbinFlags.KBIN_MPI || XHOST.MPI)){outfile << "#";}
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE_MPI]" << std::endl;
//[CT20200502 - OBSOLETE]  //be super cautious and avoid empty tags here
//[CT20200502 - OBSOLETE]  string NCPUS_VAL="MAX";
//[CT20200502 - OBSOLETE]  if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")){NCPUS_VAL=XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS");}
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=" << NCPUS_VAL << " " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]COMMAND =\"" << MPI_COMMAND_DEFAULT << "\" " << std::endl;
//[CT20200502 - OBSOLETE]  if( _kbinFlags.KBIN_MPI_AUTOTUNE ) {outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE " << std::endl;}
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]BINARY=\"";
//[CT20200502 - OBSOLETE]  if(!_kbinFlags.KBIN_MPI_BIN.empty()){outfile << _kbinFlags.KBIN_MPI_BIN;}
//[CT20200502 - OBSOLETE]  else {outfile << DEFAULT_VASP_MPI_BIN;}
//[CT20200502 - OBSOLETE]  outfile << "\"" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  //CO20180130 STOP
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  //CO20180130 - making obsolete with lines above
//[CT20200502 - OBSOLETE]  //[OBSOLETE]if( _kbinFlags.KBIN_MPI ) {
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI]" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  if( _kbinFlags.KBIN_MPI_AUTOTUNE ) {
//[CT20200502 - OBSOLETE]  //[OBSOLETE]outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  } else {
//[CT20200502 - OBSOLETE]  //[OBSOLETE]outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=MAX" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  }
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]BINARY=" << _kbinFlags.KBIN_MPI_BIN << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]} else {
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_BINARY=" << _kbinFlags.KBIN_BIN << "]" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]BINARY=\"mpi" << _kbinFlags.KBIN_BIN << "\"" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]NCPUS=MAX" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]COMMAND=\"mpirun -np\" " << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]  outfile << "[AFLOW_MODE_MPI_MODE]AUTOTUNE" << std::endl;
//[CT20200502 - OBSOLETE]  //[OBSOLETE]}
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write INCAR lines to aflow.in file 
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(AEL_data.relax_static) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_RUN]RELAX_STATIC=2" << std::endl;
//[CT20200502 - OBSOLETE]  } else if(AEL_data.static_only) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_RUN]STATIC" << std::endl;      
//[CT20200502 - OBSOLETE]  } else if(AEL_data.relax_only) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_RUN]RELAX=2" << std::endl;
//[CT20200502 - OBSOLETE]  } else {
//[CT20200502 - OBSOLETE]    aurostd::StringstreamClean(aus);
//[CT20200502 - OBSOLETE]    aus << _AELSTR_ERROR_ + "No run type selected" << endl;  
//[CT20200502 - OBSOLETE]    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]    return 1;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU1=ON"  << std::endl;
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " << _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU2=ON " << std::endl;
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " <<  _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]RELAX_IONS" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]WAVECAR=OFF" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl;
//[CT20200502 - OBSOLETE]  if(AEL_data.precaccalgonorm) {      
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]PREC=ACCURATE" << std::endl;
//[CT20200502 - OBSOLETE]    // [OBSOLETE] outfile << "[VASP_FORCE_OPTION]PREC=HIGH" << std::endl;
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]ALGO=NORMAL" << std::endl;
//[CT20200502 - OBSOLETE]    // [OBSOLETE] outfile << "[VASP_FORCE_OPTION]ALGO=FAST" << std::endl;
//[CT20200502 - OBSOLETE]  } else {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]PREC=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.content_string << std::endl;
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]ALGO=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.content_string << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  //outfile << "[VASP_FORCE_OPTION]RELAX_MODE=ENERGY" << std::endl;
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Switch off VASP symmetry - this can help when applied strains break the symmetry of the primitive cell
//[CT20200502 - OBSOLETE]  if(!AEL_data.vasp_symmetry) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]SYM=OFF" << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // [OBSOLETE] if( _vaspFlags.KBIN_VASP_FORCE_OPTION_SYM.option ) outfile << "[VASP_FORCE_OPTION]SYM=ON" << std::endl;
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry ) { outfile << "[VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << std::endl; }
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry ) { outfile << "[VASP_FORCE_OPTION]ABMIX=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme << std::endl; }
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_TYPE.isentry ) { outfile << "[VASP_FORCE_OPTION]TYPE=" << _vaspFlags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << std::endl; }
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry ) { outfile << "[VASP_FORCE_OPTION]AUTO_MAGMOM=" << (_vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option?"ON":"OFF") << std::endl; }
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.isentry &&_vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.option) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]BADER=ON" << std::endl; 
//[CT20200502 - OBSOLETE]  } else { 
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]BADER=OFF" << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_SPIN.isentry ) {
//[CT20200502 - OBSOLETE]    if(_vaspFlags.KBIN_VASP_FORCE_OPTION_SPIN.option) { 
//[CT20200502 - OBSOLETE]      outfile << "[VASP_FORCE_OPTION]SPIN=ON" << std::endl;
//[CT20200502 - OBSOLETE]    } else {
//[CT20200502 - OBSOLETE]      outfile << "[VASP_FORCE_OPTION]SPIN=OFF" << std::endl;
//[CT20200502 - OBSOLETE]    }
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  else { outfile << "[VASP_FORCE_OPTION]IGNORE_AFIX=NPARC" << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_INCAR_MODE_EXPLICIT]START" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << vaspRun.INCAR.str();
//[CT20200502 - OBSOLETE]  outfile << "[VASP_INCAR_MODE_EXPLICIT]STOP" << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write KPOINTS related lines to aflow.in file
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_MODE_IMPLICIT] " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]KSCHEME=" << vaspRun.AVASP_KSCHEME << " " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]KPPRA=" << vaspRun.AVASP_value_KPPRA << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << vaspRun.AVASP_STATIC_KSCHEME << " " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << vaspRun.AVASP_value_KPPRA_STATIC << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write POTCAR related lines to aflow.in file
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_POTCAR_MODE_IMPLICIT] " << std::endl;
//[CT20200502 - OBSOLETE]  string pp;
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
//[CT20200502 - OBSOLETE]    if(!_vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) { outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(j) << std::endl; }
//[CT20200502 - OBSOLETE]    if(_vaspFlags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry)  { outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species.at(j) << std::endl; }
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   aurostd::StringstreamClean(aus);
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
//[CT20200502 - OBSOLETE]  // [OBSOLETE]     aus << _AELSTR_MESSAGE_ + "Species_pp " << j << " = " << xvasp.str.species_pp.at(j) << endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE]     aus << _AELSTR_MESSAGE_ + "Species " << j << " = " << xvasp.str.species.at(j) << endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   }
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
//[CT20200502 - OBSOLETE]  // [OBSOLETE]  outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(j) << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE]}
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write POSCAR lines to aflow.in file
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
//[CT20200502 - OBSOLETE]  vaspRun.str.is_vasp4_poscar_format=TRUE;
//[CT20200502 - OBSOLETE]  vaspRun.str.is_vasp5_poscar_format=FALSE;
//[CT20200502 - OBSOLETE]  outfile << vaspRun.str;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  //
//[CT20200502 - OBSOLETE]  // [OBSOLETE] outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE] if(SPACES) outfile << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE] outfile << "[VASP_RUN]STATIC" << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE] if(SPACES) outfile << std::endl;
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if(AFLOWIN_QE_FLAG) {
//[CT20200502 - OBSOLETE]    outfile << AFLOWIN_SEPARATION_LINE << std::endl; 
//[CT20200502 - OBSOLETE]    outfile << "[QE_GEOM_MODE_EXPLICIT]START " << std::endl;
//[CT20200502 - OBSOLETE]    // [OBSOLETE] xstructure qestr(vaspRun.str);qestr.vasp2qe();
//[CT20200502 - OBSOLETE]    xstructure qestr(vaspRun.str);qestr.xstructure2qe();
//[CT20200502 - OBSOLETE]    outfile << qestr;
//[CT20200502 - OBSOLETE]    outfile << "[QE_GEOM_MODE_EXPLICIT]STOP " << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile.close();
//[CT20200502 - OBSOLETE]  outfile.clear();
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  return 0;
//[CT20200502 - OBSOLETE]}
//[CT20200502 - OBSOLETE]}  // namespace AEL_functions

// ************************************************************************************************
// This set of functions extract stress tensor data from VASP runs
// ************************************************************************************************

// ***************************************************************************
// AEL_functions::extractstress
// ***************************************************************************
namespace AEL_functions {
  //
  // extractstress: Extract stress tensors from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractstress(vector<_xvasp>& vaspRuns, _AEL_data& AEL_data, vector<string>& dirrunname, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile, dfile;
    xOUTCAR outcar;
    xVASPRUNXML vasprunxml;
    string vfilename, dfilename, ffilename;
    bool skipdir = false;
    aurostd::xmatrix<double> stress_tensor(3, 3);
    double pressure_val = 0.0, energy_cell_val = 0.0;
    AEL_data.normal_stress.clear();
    AEL_data.shear_stress.clear();
    AEL_data.normal_stress.resize(AEL_data.normal_strain.size());
    AEL_data.shear_stress.resize(AEL_data.shear_strain.size());
    AEL_data.normal_deformations_complete.resize(AEL_data.normal_strain.size());
    AEL_data.shear_deformations_complete.resize(AEL_data.shear_strain.size());
    AEL_data.energycalculated.clear();
    AEL_data.pressurecalculated.clear();
    AEL_data.stresscalculated.clear();
    AEL_data.structurecalculated.clear();    
    if(AEL_data.relax_static || AEL_data.static_only) {
      if(AEL_data.vasprunxmlstress) {
        aurostd::string2tokens(string("vasprun.xml.static.bz2,vasprun.xml.static.gz,vasprun.xml.static.xz,vasprun.xml.static"),vfile,",");
      } else {
        aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"),vfile,",");
      }
    } else if(AEL_data.relax_only) {
      if(AEL_data.vasprunxmlstress) {
        aurostd::string2tokens(string("vasprun.xml.relax2.bz2,vasprun.xml.relax2.gz,vasprun.xml.relax2.xz,vasprun.xml.relax2"),vfile,",");
      } else {
        aurostd::string2tokens(string("OUTCAR.relax2.bz2,OUTCAR.relax2.gz,OUTCAR.relax2.xz,OUTCAR.relax2"),vfile,",");
      }
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "No run type selected" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 1;
    }

    // Loop over normal strains and read stress tensors
    uint idVaspRun = 0;
    for (uint i = 1; i <= AEL_data.normal_strain.size(); i++) {
      for (uint j = 0; j < AEL_data.normal_deformations.size(); j++) {   
        skipdir = false;
        if(idVaspRun > vaspRuns.size()) {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        double strainfactor = 1.0 + AEL_data.normal_deformations.at(j);
        aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory << endl;
        aus << _AELSTR_MESSAGE_ + "Normal stress: i = " << i << ", j = " << j << ", strain factor = " << strainfactor << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
        dfilename = dfile.at(dfile.size()-1);
        aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
        aus << _AELSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // Check if structure is on list of failed runs to be skipped
        // If so, then skip reading and continue to next structure
        for (uint ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
          ffilename = AEL_data.failed_arun_list.at(ij);
          if(LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
            aus << _AELSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          // [OBSOLETE] if(dfilename == AEL_data.failed_arun_list.at(ij)) continue;
          // [OBSOLETE] if(strncmp(dfilename.c_str(), ffilename.c_str(), 12) == 0) continue;
          if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            skipdir = true;
          }
        }
        if(skipdir) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          idVaspRun++;
          skipdir = false;
          continue;
        }
        // if(skipdir) continue;
        // continue;

        // If tarred and compressed directory exists...
        if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.bz2") ) {
          aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.bz2" );
        } else if ( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2") ) { 
          aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2" );
        } // Extract all...
        if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.gz") ) {
          aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.gz" );
        } else if ( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.gz") ) {
          aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.gz" );
        } // Extract all...
        if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.xz") ) {
          aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.xz" );
        } else if ( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.xz") ) {
          aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.xz" );
        } // Extract all...

        // If the LOCK file is missing, then it is probably a corrupted run
        // Do not accept it and wait for the new run
        // [OBSOLETE] if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/LOCK") ) )
        aurostd::StringstreamClean(aus);
        aus <<  _AELSTR_MESSAGE_ + "LOCK file path = " << dirrunname.at(idVaspRun) + "/" + _AFLOWLOCK_ << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) && !aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) &&
            !(((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
                ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) || (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
                 (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/ael.LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
                 (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/ael.LOCK")))))) {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if(AEL_data.autoskipfailedaruns) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            idVaspRun++;
            continue;
          } else {
            throw AELStageBreak();
          }
        }
        if(AEL_data.vasprunxmlstress) {
          for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
            if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            } else if (aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname.at(idVaspRun)+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            }	  
          }
          if(vasprunxml.content=="") {
            aurostd::StringstreamClean(aus);
            aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -vasprunxml.stress;
          vasprunxml.clear();
        } else {
          for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
            if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname.at(idVaspRun)+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              outcar.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            }
          }
          if(outcar.content=="") {
            aurostd::StringstreamClean(aus);
            aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -outcar.stress;
          energy_cell_val = outcar.energy_cell;
          pressure_val = outcar.pressure_residual;
          outcar.clear();
        }
        AEL_data.normal_stress.at(i-1).push_back(stress_tensor);
        AEL_data.normal_deformations_complete.at(i-1).push_back(AEL_data.normal_deformations.at(j));
        AEL_data.energycalculated.push_back(energy_cell_val);
        AEL_data.pressurecalculated.push_back(pressure_val);
        AEL_data.stresscalculated.push_back(stress_tensor);
        AEL_data.structurecalculated.push_back(idVaspRun);
        idVaspRun++;
        // Print out stress tensor
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Normal stress tensor = " << AEL_data.normal_stress.at(i-1).at(AEL_data.normal_stress.at(i-1).size()-1) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    // Loop over shear strains and read stress tensors
    for (uint i = 1; i <= AEL_data.shear_strain.size(); i++) {
      for (uint j = 0; j < AEL_data.shear_deformations.size(); j++) {   
        skipdir = false;
        if(idVaspRun > vaspRuns.size()) {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          return 1;
        }
        double strainfactor = 1.0 + AEL_data.shear_deformations.at(j);
        aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory << endl;
        aus << _AELSTR_MESSAGE_ + "Shear stress: i = " << i << ", j = " << j << ", strain factor = " << strainfactor << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
        dfilename = dfile.at(dfile.size()-1);
        aurostd::StringstreamClean(aus);
        // [OBSOLETE] aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory.at(idVaspRun) << endl;
        aus << _AELSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // Check if structure is on list of failed runs to be skipped
        // If so, then skip reading and continue to next structure
        for (uint ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
          ffilename = AEL_data.failed_arun_list.at(ij);
          if(LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
            aus << _AELSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          }
          // [OBSOLETE] if(dfilename == AEL_data.failed_arun_list.at(ij)) continue;
          // [OBSOLETE] if(strncmp(dfilename.c_str(), ffilename.c_str(), 12) == 0) continue;
          if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            skipdir = true;
          }
        }
        if(skipdir) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = false;
          idVaspRun++;
          continue;
        }

        // If tarred and compressed directory exists...
        if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.bz2") ) {
          aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.bz2" );
        } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2") ) {
          aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2" );
        } // Extract all...
        if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.gz") ) {
          aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.gz" );
        } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.gz") ) {
          aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.gz" );
        } // Extract all...
        if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.xz") ) {
          aurostd::execute( string("tar -xf") + vaspRuns.at(idVaspRun).Directory + ".tar.xz" );
        } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.xz") ) {
          aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.xz" );
        } // Extract all...

        // If the LOCK file is missing, then it is probably a corrupted run
        // Do not accept it and wait for the new run
        // [OBSOLETE] if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + string("/LOCK") ) )
        if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ )  && !aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) &&
            !(((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
                ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) || (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
                 (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/ael.LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
                 (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/ael.LOCK")))))) {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          if(AEL_data.autoskipfailedaruns) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            idVaspRun++;
            continue;
          } else {
            throw AELStageBreak();
          }
        }
        if(AEL_data.vasprunxmlstress) {
          for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
            if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname.at(idVaspRun)+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            }
          }
          if(vasprunxml.content=="") {
            aurostd::StringstreamClean(aus);
            aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -vasprunxml.stress;
          vasprunxml.clear();
        } else {
          for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
            if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(ij))) {
              if(LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname.at(idVaspRun)+"/"+vfile.at(ij) << endl;	    
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              outcar.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(ij));
              vfilename = vfile.at(ij);
            }
          }
          if(outcar.content=="") {
            aurostd::StringstreamClean(aus);
            aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            if(AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -outcar.stress;
          energy_cell_val = outcar.energy_cell;
          pressure_val = outcar.pressure_residual;
          outcar.clear();
        }
        AEL_data.shear_stress.at(i-1).push_back(stress_tensor);
        AEL_data.shear_deformations_complete.at(i-1).push_back(AEL_data.shear_deformations.at(j));
        AEL_data.energycalculated.push_back(energy_cell_val);
        AEL_data.pressurecalculated.push_back(pressure_val);
        AEL_data.stresscalculated.push_back(stress_tensor);
        AEL_data.structurecalculated.push_back(idVaspRun);
        // Print out stress tensor
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Shear stress tensor = " << AEL_data.shear_stress.at(i-1).at(AEL_data.shear_stress.at(i-1).size()-1)  << endl;
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", energy = " << outcar.energy_cell << endl;
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", pressure = " << outcar.pressure_residual << endl;		
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        idVaspRun++;
      }
    }
    if(AEL_data.calcstrainorigin) {
      if(idVaspRun > vaspRuns.size()) {
        aurostd::StringstreamClean(aus);
        aus <<  _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns.at(idVaspRun).Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // If tarred and compressed directory exists...
      // [OBSOLETE]    string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.bz2") ) {
        aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.bz2" );
      } else if ( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2") ) { 
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2" );
      } // Extract all...
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.gz") ) {
        aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.gz" );
      } else if ( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2") ) { 
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2" );
      } // Extract all...
      if( aurostd::FileExist(vaspRuns.at(idVaspRun).Directory + ".tar.xz") ) {
        aurostd::execute( string("tar -xf ") + vaspRuns.at(idVaspRun).Directory + ".tar.xz" );
      }  else if ( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2") ) { 
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2" );
      } // Extract all...

      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/" + _AFLOWLOCK_ ) && !aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ )  &&
          !(((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) &&
              ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) || (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) ||
               (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/ael.LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
               (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/ael.LOCK")))))) {
        aurostd::StringstreamClean(aus);
        aus <<  _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        throw AELStageBreak();
      }
      if(AEL_data.vasprunxmlstress) {
        for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
          if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            vasprunxml.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
            vfilename = vfile.at(ij);
          } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(ij))) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << dirrunname.at(idVaspRun)+"/"+vfile.at(ij) << endl;	    
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            vasprunxml.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(ij));
            vfilename = vfile.at(ij);
          }
        }
        if(vasprunxml.content=="") {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -vasprunxml.stress;
        vasprunxml.clear();
      } else {
        for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
          if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij))) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij) << endl;	    
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(ij));
            vfilename = vfile.at(ij);
          } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(ij))) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << dirrunname.at(idVaspRun)+"/"+vfile.at(ij) << endl;	    
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            outcar.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(ij));
            vfilename = vfile.at(ij);
          }
        }
        if(outcar.content=="") {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -outcar.stress;
        energy_cell_val = outcar.energy_cell;
        pressure_val = outcar.pressure_residual;
        outcar.clear();
      }
      AEL_data.origin_stress.push_back(stress_tensor);
      AEL_data.energycalculated.push_back(energy_cell_val);
      AEL_data.pressurecalculated.push_back(pressure_val);
      AEL_data.stresscalculated.push_back(stress_tensor);
      AEL_data.structurecalculated.push_back(idVaspRun); 
      // Print out stress tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Origin stress tensor = " << AEL_data.origin_stress.at(0)  << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      idVaspRun++;
    } else if(AEL_data.fitrelaxedstruct) {
      string relaxedstructdirname = AEL_data.dirpathname + "/";
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Initial relaxed structure directory = " << relaxedstructdirname << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

      // If tarred and compressed directory exists...
      // [OBSOLETE] string tarfilename = relaxedstructdirname + ".tar.bz2";
      if( aurostd::FileExist(relaxedstructdirname + ".tar.bz2") ) { aurostd::execute( string("tar -xf ") + relaxedstructdirname + ".tar.bz2" ); } // Extract all...
      if( aurostd::FileExist(relaxedstructdirname + ".tar.gz") ) { aurostd::execute( string("tar -xf ") + relaxedstructdirname + ".tar.gz" ); } // Extract all...
      if( aurostd::FileExist(relaxedstructdirname + ".tar.xz") ) { aurostd::execute( string("tar -xf ") + relaxedstructdirname + ".tar.xz" ); } // Extract all...

      if(AEL_data.vasprunxmlstress) {
        for(uint ij=0;ij<vfile.size()&&(vasprunxml.content=="");ij++) {
          if(aurostd::FileExist(relaxedstructdirname+"/"+vfile.at(ij))) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << relaxedstructdirname+"/"+vfile.at(ij) << endl;	    
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            vasprunxml.GetPropertiesFile(relaxedstructdirname+"/"+vfile.at(ij));
            vfilename = vfile.at(ij);
          }
        }
        if(vasprunxml.content=="") {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << relaxedstructdirname << " directory is missing." << endl;
          aus <<  _AELSTR_ERROR_ + "The flag [AFLOW_AEL]FITRELAXEDSTRUCT=ON is set in the input file." << endl;
          aus <<  _AELSTR_ERROR_ + "This requires the results of a relaxed calculation to be present in the directory " << relaxedstructdirname << endl;
          aus <<  _AELSTR_ERROR_ + "Either set the flag to OFF or run AEL in the appropriate directory." << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -vasprunxml.stress;
        vasprunxml.clear();
      } else {
        for(uint ij=0;ij<vfile.size()&&(outcar.content=="");ij++) {
          if(aurostd::FileExist(relaxedstructdirname+"/"+vfile.at(ij))) {
            if(LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile.at(ij) << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << relaxedstructdirname+"/"+vfile.at(ij) << endl;	    
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            outcar.GetPropertiesFile(relaxedstructdirname+"/"+vfile.at(ij));
            vfilename = vfile.at(ij);
          }
        }
        if(outcar.content=="") {
          aurostd::StringstreamClean(aus);
          aus <<  _AELSTR_WARNING_ + "The " << vfilename << " file in " << relaxedstructdirname << " directory is missing." << endl;
          aus <<  _AELSTR_ERROR_ + "The flag [AFLOW_AEL]FITRELAXEDSTRUCT=ON is set in the input file." << endl;
          aus <<  _AELSTR_ERROR_ + "This requires the results of a relaxed calculation to be present in the directory " << relaxedstructdirname << endl;
          aus <<  _AELSTR_ERROR_ + "Either set the flag to OFF or run AEL in the appropriate directory." << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -outcar.stress;
        energy_cell_val = outcar.energy_cell;
        pressure_val = outcar.pressure_residual;
        outcar.clear();
      }
      AEL_data.origin_stress.push_back(stress_tensor);
      AEL_data.energycalculated.push_back(energy_cell_val);
      AEL_data.pressurecalculated.push_back(pressure_val);
      AEL_data.stresscalculated.push_back(stress_tensor);
      AEL_data.structurecalculated.push_back(idVaspRun);
      // Print out stress tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Origin stress tensor = " << AEL_data.origin_stress.at(0)  << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      idVaspRun++;
    }        
    return 0;
  }
} // namespace AEL_functions

// **************************************************************************
//  End of AFLOW AEL set-up and extract stress-strain data
// **************************************************************************

#endif  // _AFLOW_AEL_GET_STRESS_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
