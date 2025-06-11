// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_GET_EV_CPP
#define _AFLOW_AGL_GET_EV_CPP
#include "aflow.h"
#include "aflow_agl_debye.h"


// ###############################################################################
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2021)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//

// *****************************************************************************************************************
// The following functions are for setting up AGL inputs for postprocessing runs called from other parts of AFLOW
// *****************************************************************************************************************
namespace AGL_functions {
  uint AGL_xvasp_flags_populate(_xvasp& xvasp, string& AflowIn, string& AflowInName, string& FileLockName, const string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) {
    ifstream FileAFLOWIN, FileAFLOWINcheck;
    string FileNameAFLOWIN = "", FileNameAFLOWINcheck = "", AflowInCheck = "";
    string FileNameMessage = "";
    ostringstream aus;
    vector<string> vAflowInCheck;
    bool agl_aflowin_found = false;
    bool Krun = true;
    bool load_POSCAR_from_xvasp = false;
    // Set aflags
    aflags.Directory=directory_LIB;
    if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;
    aflags.KBIN_RUN_AFLOWIN=TRUE;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN=FALSE;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP=FALSE;
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN=FALSE;
    aflags.KBIN_DELETE_AFLOWIN=FALSE;
    //CT20200721 The following code has been moved to the function AGL_Get_AflowInName
    // Set FileMESSAGE name
    // [OBSOLETE] if(FileLockName.length() > 0) { //CT20200624 Moved down: only rename LOCK files if it is confirmed that this is an AGL main directory and not an ARUN.AGL directory
    // [OBSOLETE]   if (aurostd::FileExist(directory_LIB+"/"+FileLockName)) {
    // [OBSOLETE]     aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+" ")+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+".run"));   
    // [OBSOLETE]   }
    // [OBSOLETE]   string FileNameMessage=directory_LIB+"/"+FileLockName;
    // [OBSOLETE]   FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
    // [OBSOLETE] } else {
    // [OBSOLETE]   if (aurostd::FileExist(directory_LIB+"/agl.LOCK")) {
    // [OBSOLETE]     aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/agl.LOCK ")+aurostd::CleanFileName(directory_LIB+"/agl.LOCK.run"));
    // [OBSOLETE]   }
    // [OBSOLETE]   string FileNameMessage=directory_LIB+"/agl.LOCK";
    // [OBSOLETE]   FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
    // [OBSOLETE] }
    //CO20200502 START - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future
    // [OBSOLETE] vector<string> vaflowins;
    // [OBSOLETE] if(AflowInName.length()>0){vaflowins.push_back(AflowInName);} // Check if AflowInName exists
    // [OBSOLETE] if(_AFLOWIN_.length()>0){vaflowins.push_back(_AFLOWIN_);} // Otherwise, check if _AFLOWIN_ file is AGL input file
    // [OBSOLETE] vaflowins.push_back(_AFLOWIN_AGL_DEFAULT_);  // Otherwise, check for other commonly used names for AGL aflow.in file:
    // [OBSOLETE] for(uint iaf=0;iaf<vaflowins.size()&&!agl_aflowin_found;iaf++){
    // [OBSOLETE]   const string& aflowinname = vaflowins.at(iaf);
    // [OBSOLETE]   if((!agl_aflowin_found) && (aurostd::FileExist(directory_LIB+"/"+aflowinname))) {
    // [OBSOLETE]     FileNameAFLOWINcheck = directory_LIB+"/"+aflowinname;
    // [OBSOLETE]     cerr << "FileNameAFLOWINcheck = " << FileNameAFLOWINcheck << endl;
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
    // [OBSOLETE]     // Check if aflowinname contains command to run AGL
    // [OBSOLETE]     for(uint i=0;i<vAflowInCheck.size()&&!agl_aflowin_found;i++){
    // [OBSOLETE]       if((aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowInCheck,"[VASP_AGL]CALC",TRUE)) &&
    // [OBSOLETE]           !(aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC_",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALC_",TRUE) ||
    // [OBSOLETE]             aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALCS",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALCS",TRUE) || FALSE)){
    // [OBSOLETE]         FileNameAFLOWIN = FileNameAFLOWINcheck;
    // [OBSOLETE]         agl_aflowin_found = true;
    // [OBSOLETE]       }
    // [OBSOLETE]     }
    // [OBSOLETE]     FileAFLOWINcheck.close();
    // [OBSOLETE]   }
    // [OBSOLETE] }

    //CT20200721 Calls function AGL_Get_AflowInName to find correct aflow.in filename
    AGL_functions::AGL_Get_AflowInName(AflowInName, directory_LIB, agl_aflowin_found);
    FileNameAFLOWIN = directory_LIB+"/"+AflowInName;

    //CO20200502 STOP - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future
    if (agl_aflowin_found) {
      // Set FileMESSAGE name
      //CT20200624 Moved down so LOCK file is only renamed if this is actually an AGL main directory and not an ARUN.AGL directory
      // [OBSOLETE] cerr << _AGLSTR_MESSAGE_ << "FileLockName = " << FileLockName << endl;
      if(FileLockName.length() > 0) {
        if (aurostd::FileExist(directory_LIB+"/"+FileLockName)) {
          aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+" ")+aurostd::CleanFileName(directory_LIB+"/"+FileLockName+".run"));   
        }
        string FileNameMessage=directory_LIB+"/"+FileLockName;
        FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
      } else {
        if (aurostd::FileExist(directory_LIB+"/agl.LOCK")) {
          aurostd::execute("mv "+aurostd::CleanFileName(directory_LIB+"/agl.LOCK ")+aurostd::CleanFileName(directory_LIB+"/agl.LOCK.run"));
        }
        string FileNameMessage=directory_LIB+"/agl.LOCK";
        FileMESSAGE.open(FileNameMessage.c_str(),std::ios::app);
      }
      // Search for AGL aflow.in filename
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AFLOW Input file name = " << FileNameAFLOWIN << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
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
      xvasp.clear();
      uint ixvasp=0;
      xvasp.POSCAR_index=ixvasp;
      KBIN::readModulesFromAflowIn(AflowIn, kflags, xvasp);
      xvasp.Directory=aflags.Directory;
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags,load_POSCAR_from_xvasp));
      if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
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
      cerr << _AGLSTR_MESSAGE_ << "AGL input file not found!" << endl; //CT20200624 If no AGL file present, then this is not an AGL main directory: write to cerr and return 2 to inform calling function
      cerr << _AGLSTR_MESSAGE_ << "Not an AGL main directory" << endl;
      // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 2;
    }     
  }
}


// *****************************************************************************************************************
// Finds aflow.in file for AGL calculation (if it exists) //CT20200713
// *****************************************************************************************************************
namespace AGL_functions {
  uint AGL_Get_AflowInName(string& AflowInName, const string& directory_LIB, bool& agl_aflowin_found) {
    ifstream FileAFLOWIN, FileAFLOWINcheck;
    string FileNameAFLOWIN = "", FileNameAFLOWINcheck = "", AflowInCheck = "";
    string FileNameMessage = "";
    string aflowinname = "";
    string stmp="";
    uint aglerror = 0, filelength = 0;
    vector<string> vAflowInCheck;
    agl_aflowin_found = false;
    vector<string> vaflowins;
    if(AflowInName.length()>0){vaflowins.push_back(AflowInName);} // Check if AflowInName exists
    if(_AFLOWIN_.length()>0){vaflowins.push_back(_AFLOWIN_);} // Otherwise, check if _AFLOWIN_ file is AGL input file
    vaflowins.push_back(_AFLOWIN_AGL_DEFAULT_);  // Otherwise, check for other commonly used names for AGL aflow.in file
    for(uint iaf=0;iaf<vaflowins.size()&&!agl_aflowin_found;iaf++){
      aflowinname = vaflowins.at(iaf);
      if(aurostd::EFileExist(directory_LIB+"/"+aflowinname,stmp)&&aurostd::IsCompressed(stmp)){aurostd::UncompressFile(stmp);}  //CO20210204 - fix aflow.in.xz
      if((!agl_aflowin_found) && (aurostd::FileExist(directory_LIB+"/"+aflowinname))) {
        FileNameAFLOWINcheck = directory_LIB+"/"+aflowinname;
        // [OBSOLETE] FileAFLOWINcheck.open(FileNameAFLOWINcheck.c_str(),std::ios::in);
        // [OBSOLETE] FileAFLOWINcheck.clear();
        // [OBSOLETE] FileAFLOWINcheck.seekg(0);
        AflowInCheck="";
        // [OBSOLETE] char c;
        // READ aflowinname and put into AflowInCheck
        // [OBSOLETE] while (FileAFLOWINcheck.get(c)) {
        // [OBSOLETE]  AflowInCheck+=c;
        // [OBSOLETE] }
        // [OBSOLETE] FileAFLOWINcheck.clear();
        // [OBSOLETE] FileAFLOWINcheck.seekg(0);
        filelength=aurostd::file2string(FileNameAFLOWINcheck, AflowInCheck);
        if (filelength > 0) {
          aglerror = 0;
        } else {
          aglerror = 1;
        }
        AflowInCheck=aurostd::RemoveComments(AflowInCheck); // NOW Clean AFLOWIN
        vAflowInCheck.clear();
        aurostd::string2vectorstring(AflowInCheck,vAflowInCheck); 
        // Check if aflowinname contains command to run AGL
        for(uint i=0;i<vAflowInCheck.size()&&!agl_aflowin_found;i++){
          if((aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowInCheck,"[VASP_AGL]CALC",TRUE)) &&
              !(aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALC_",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALC_",TRUE) ||
                aurostd::substring2bool(vAflowInCheck[i],"[AFLOW_AGL]CALCS",TRUE) || aurostd::substring2bool(vAflowInCheck[i],"[VASP_AGL]CALCS",TRUE) || FALSE)){
            FileNameAFLOWIN = FileNameAFLOWINcheck;
            agl_aflowin_found = true;
            AflowInName = aflowinname;
          }
        }
        FileAFLOWINcheck.close();
      }
    }
    return aglerror;
  }
}

// *******************************************************************************
// The following functions are for generating _AFLOWIN_ files
// *******************************************************************************

// ***************************************************************************
// AGL_functions::aglvaspflags
// ***************************************************************************
namespace AGL_functions {
  //
  // Function to assign values for VASP input flags from _AFLOWIN_ file to vaspRun _xvasp class
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  // [OBSOLETE] uint aglvaspflags(_xvasp& vaspRun, _vflags& _vaspFlags, _kflags& _kbinFlags, string& runname, ofstream& FileMESSAGE)
  uint aglvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, string& dirrunname, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile, dfile;
    string vfilename = "", dfilename = "", ffilename = "";
    bool vfileexist = false;
    bool skipdir = false;
    aurostd::string2tokens(dirrunname, dfile, "/");
    dfilename = dfile.at(dfile.size()-1);    
    if(AGL_data.relax_static || AGL_data.static_only) {
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
    // SOME WARNINGS: check existence of LOCK and OUTCAR.static files
    // [OBSOLETE] if( !aurostd::FileExist( vaspRun.Directory + "/"+_AFLOWLOCK_ ) && aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.static") ) )
    if( !(aurostd::FileExist( dirrunname + "/" + _AFLOWLOCK_ ) ||
          ((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/LOCK")))) &&
        ( vfileexist ) ) {
      aurostd::StringstreamClean(aus);
      // [OBOLSETE] aus << _AGLSTR_WARNING_ + "found OUTCAR.static but no LOCK in " <<  vaspRun.Directory << endl;
      aus << _AGLSTR_WARNING_ + "found " << vfilename << " but no " << _AFLOWLOCK_ << " in " <<  dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list.at(ij);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = true;
        }
      }
      if(skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Directory: " << dfilename << " will be skipped." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        skipdir = false;
      } else {
        return 1;
      }
    }

    // [OBSOLETE] if( aurostd::FileExist( vaspRun.Directory + "/"+_AFLOWLOCK_ ) && !aurostd::FileExist( vaspRun.Directory + string("/OUTCAR.static") ) )
    if( (aurostd::FileExist( dirrunname + "/" + _AFLOWLOCK_ ) ||
          ((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/LOCK")))) &&
        !(vfileexist) ) {
      aurostd::StringstreamClean(aus);
      // [OBSOLETE] aus << _AGLSTR_WARNING_ + "found LOCK but no OUTCAR.static in " <<  vaspRun.Directory << endl;
      aus << _AGLSTR_WARNING_ + "found " << _AFLOWLOCK_ << " but no OUTCAR in " << dirrunname << endl;      
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list.at(ij);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = true;
        }
      }
      if(skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Directory: " << dfilename << " will be skipped." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        skipdir = false;
      } else {      
        return 1;
      }
    }

    // Switch off autotune
    kbinFlags.KBIN_MPI_AUTOTUNE = true;

    if(AGL_data.relax_static) {
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
      // Increase DOSCAR density for STATIC-only runs
      vaspRun.aopts.flag("FLAG::EXTRA_INCAR", true);
      vaspRun.AVASP_EXTRA_INCAR << "# Added by [AFLOW_AGL] begin" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "EMIN= -30.0    # For finer DOS grid" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "EMAX=  45.0    # For finer DOS grid" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "NEDOS= 5000    # For finer DOS grid" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "# Added by [AFLOW_AGL] end" << std::endl;    
    }

    if(AGL_data.precaccalgonorm) {
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.isentry = true;      
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.content_string = "ACCURATE";
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.isentry = true;      
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.content_string = "NORMAL";
    }

    // Common KPOINTS settings and OVERRIDES
    // [OBSOLETE] vaspRun.AVASP_KSCHEME = _vaspFlags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
    // [OBSOLETE] vaspRun.AVASP_value_KPPRA = _vaspFlags.KBIN_VASP_KPOINTS_KPPRA.content_int;
    // [OBSOLETE] vaspRun.AVASP_STATIC_KSCHEME = _vaspFlags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;
    // [OBSOLETE] vaspRun.AVASP_value_KPPRA_STATIC = _vaspFlags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_int;

    // Clear old INCAR and set it as we want...
    // [OBSOLETE] vaspRun.INCAR.str(std::string());
    // [OBSOLETE] string system;
    // [OBSOLETE] for(uint j=0; j < vaspRun.str.species.size(); j++)
    // [OBSOLETE]   system = system + vaspRun.str.species_pp.at(j);
    // [OBSOLETE] system = system + "@" + runname;
    // [OBSOLETE] vaspRun.INCAR << "SYSTEM=" << system << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "# Added by [AFLOW_AGL] begin" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
    // Added July 2018 to increase DOSCAR grid density
    // Temporary measure until aflow.in generation with AVASP is working with APL, AEL and AGL
    // [OBSOLETE] vaspRun.INCAR << "EMIN= -30.0    # For finer DOS grid" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "EMAX=  45.0    # For finer DOS grid" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "NEDOS= 5000    # For finer DOS grid" << std::endl;
    // [OBSOLETE] vaspRun.INCAR << "# Added by [AFLOW_AGL] end" << std::endl;    

    // Change format of POSCAR
    if( ( !kbinFlags.KBIN_MPI && ( kbinFlags.KBIN_BIN.find("46") != string::npos ) ) ||
        (  kbinFlags.KBIN_MPI && ( kbinFlags.KBIN_MPI_BIN.find("46") != string::npos ) ) ) {
      vaspRun.str.is_vasp5_poscar_format = false; 
    }
    return 0;
  }
} // namespace AGL_functions

//[CT20200502 - OBSOLETE]  // ***************************************************************************
//[CT20200502 - OBSOLETE]  // AGL_functions::createAFLOWIN
//[CT20200502 - OBSOLETE]  // ***************************************************************************
//[CT20200502 - OBSOLETE]  namespace AGL_functions {
//[CT20200502 - OBSOLETE]  //
//[CT20200502 - OBSOLETE]  // Create _AFLOWIN_ file: makes new directory and writes _AFLOWIN_ for strained structure file inside it 
//[CT20200502 - OBSOLETE]  // Adapted from that in AFLOW APL function PhononCalculator::createAFLOWIN()
//[CT20200502 - OBSOLETE]  //
//[CT20200502 - OBSOLETE]  uint createAFLOWIN(_xvasp& vaspRun, _xvasp& xvasp, _kflags& _kbinFlags, _vflags& _vaspFlags, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
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
//[CT20200502 - OBSOLETE]  string filename =  vaspRun.Directory + "/"+_AFLOWIN_;
//[CT20200502 - OBSOLETE]  ofstream outfile(filename.c_str(),ios_base::out);
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Check _AFLOWIN_ file is open
//[CT20200502 - OBSOLETE]  if( !outfile.is_open() ) {
//[CT20200502 - OBSOLETE]  aurostd::StringstreamClean(aus);
//[CT20200502 - OBSOLETE]  aus << _AGLSTR_WARNING_ + "Cannot create [" << _AFLOWIN_ << "] file" << endl;
//[CT20200502 - OBSOLETE]  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]  return 1;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  //CHMOD a+rw _AFLOWIN_: change permissions on _AFLOWIN_ file
//[CT20200502 - OBSOLETE]  aurostd::ChmodFile("a+rw",filename);
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write to _AFLOWIN_ file
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] _ ___ _" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] / \\|  || \\ |" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] | o | _ | " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[AFLOW] |_n_|__||___| automatic generated file" << std::endl;
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
//[CT20200502 - OBSOLETE]  // Write INCAR lines to _AFLOWIN_ file 
//[CT20200502 - OBSOLETE]  //
//[CT20200502 - OBSOLETE]  // [OBSOLETE] outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE] if(SPACES) outfile << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE] outfile << "[VASP_RUN]STATIC" << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE] if(SPACES) outfile << std::endl;
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(AGL_data.relax_static) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_RUN]RELAX_STATIC=2" << std::endl;
//[CT20200502 - OBSOLETE]  } else if(AGL_data.static_only) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_RUN]STATIC" << std::endl;      
//[CT20200502 - OBSOLETE]  } else {
//[CT20200502 - OBSOLETE]    aurostd::StringstreamClean(aus);
//[CT20200502 - OBSOLETE]    aus << _AGLSTR_ERROR_ + "No run type selected" << endl;  
//[CT20200502 - OBSOLETE]    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]    return 1;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU1= ON"  << std::endl;
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " << _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  if(_vaspFlags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU2=ON " << std::endl;
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]LDAU_PARAMETERS= " <<  _vaspFlags.KBIN_VASP_LDAU_PARAMETERS  << std::endl;
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]RELAX_IONS" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]WAVECAR=OFF" << std::endl;
//[CT20200502 - OBSOLETE]  if( _vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.isentry &&_vaspFlags.KBIN_VASP_FORCE_OPTION_BADER.option) { 
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]CHGCAR=ON" << std::endl; 
//[CT20200502 - OBSOLETE]  } else {
//[CT20200502 - OBSOLETE]    outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl; 
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  // [OBSOLETE] outfile << "[VASP_FORCE_OPTION]CHGCAR=OFF" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]PREC=ACCURATE" << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_FORCE_OPTION]ALGO=NORMAL" << std::endl;
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
//[CT20200502 - OBSOLETE]  // Write KPOINTS related lines to _AFLOWIN_ file
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_MODE_IMPLICIT] " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]KSCHEME=" << vaspRun.AVASP_KSCHEME << " " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]KPPRA=" << vaspRun.AVASP_value_KPPRA << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]STATIC_KSCHEME=" << vaspRun.AVASP_STATIC_KSCHEME << " " << std::endl;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_KPOINTS_FILE]STATIC_KPPRA=" << vaspRun.AVASP_value_KPPRA_STATIC << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write POTCAR related lines to _AFLOWIN_ file
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
//[CT20200502 - OBSOLETE]  // [OBSOLETE]     aus << _AGLSTR_MESSAGE_ << "Species_pp " << j << " = " << xvasp.str.species_pp.at(j) << endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE]     aus << _AGLSTR_MESSAGE_ << "Species " << j << " = " << xvasp.str.species.at(j) << endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   }
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CT20200502 - OBSOLETE]  // [OBSOLETE]   for(uint j=0; j < xvasp.str.species_pp.size(); j++) {
//[CT20200502 - OBSOLETE]  // [OBSOLETE]  outfile << "[VASP_POTCAR_FILE]" << xvasp.str.species_pp.at(j) << std::endl;
//[CT20200502 - OBSOLETE]  // [OBSOLETE]}
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]
//[CT20200502 - OBSOLETE]  // Write POSCAR lines to _AFLOWIN_ file
//[CT20200502 - OBSOLETE]  outfile << AFLOWIN_SEPARATION_LINE << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
//[CT20200502 - OBSOLETE]  outfile << "[VASP_POSCAR_MODE_EXPLICIT]START " << std::endl;
//[CT20200502 - OBSOLETE]  vaspRun.str.is_vasp4_poscar_format=TRUE;
//[CT20200502 - OBSOLETE]  vaspRun.str.is_vasp5_poscar_format=FALSE;
//[CT20200502 - OBSOLETE]  outfile << vaspRun.str;
//[CT20200502 - OBSOLETE]  outfile << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << std::endl;
//[CT20200502 - OBSOLETE]  if(SPACES) { outfile << std::endl; }
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
//[CT20200502 - OBSOLETE]  }
//[CT20200502 - OBSOLETE]  } // namespace AGL_functions

// ************************************************************************************************
// This set of functions extract, sort and check (E, V) data from VASP runs
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::extractenerg
// ***************************************************************************
namespace AGL_functions {
  //
  // extractenerg: Extract final energies from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractenerg(vector<_xvasp>& vaspRuns, _AGL_data& AGL_data, vector<string>& dirrunname, ofstream& FileMESSAGE) {
    bool LVERBOSE=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile, dfile;
    aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"),vfile,",");
    xOUTCAR outcar;
    string vfilename, dfilename, ffilename;
    bool skipdir = false;
    aurostd::xmatrix<double> stress_tensor(3, 3);
    AGL_data.volumeinput.clear();
    AGL_data.energyinput.clear();
    AGL_data.pressurecalculated.clear();
    AGL_data.stresscalculated.clear();
    for(uint idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
      skipdir = false;
      aurostd::StringstreamClean(aus);
      // Print out total energy
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aurostd::string2tokens(vaspRuns.at(idVaspRun).Directory, dfile, "/");
      dfilename = dfile.at(dfile.size()-1);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (uint ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list.at(ij);
        if(LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(aurostd::substring2bool(dfilename,ffilename,TRUE)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          skipdir = true;
        }
      }
      if(skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        skipdir = false;
        continue;
      }

      // If tarred and compressed directory exists...
      // [OBSOLETE]      string tarfilename = vaspRuns.at(idVaspRun).Directory + ".tar.bz2";
      // [OBSOLETE]    if( aurostd::FileExist(tarfilename) ) {
      // [OBSOLETE]   	// Extract all...
      // [OBSOLETE]   	aurostd::execute( string("tar -xf ") + tarfilename );
      // [OBSOLETE]     }
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
      } else if( aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.xz") ) {
        aurostd::execute( string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.xz" );
      } // Extract all...

      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if( !aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/"+_AFLOWLOCK_ ) && !aurostd::FileExist( dirrunname.at(idVaspRun) + "/"+_AFLOWLOCK_ ) && !((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) && ((aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/agl.LOCK")) || (aurostd::FileExist( vaspRuns.at(idVaspRun).Directory + "/LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/agl.LOCK")) || (aurostd::FileExist( dirrunname.at(idVaspRun) + "/LOCK"))))) { //CT20200625 Modify to check for other LOCK file names if postprocessing run
        aurostd::StringstreamClean(aus);
        aus <<  _AGLSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }

      //for(uint i=0;i<vfile.size()&&(outcar.outcar=="");i++)
      for(uint i=0;i<vfile.size()&&(outcar.content=="");i++) {
        if(aurostd::FileExist(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(i))) {
          outcar.GetPropertiesFile(vaspRuns.at(idVaspRun).Directory+"/"+vfile.at(i));
        } else if(aurostd::FileExist(dirrunname.at(idVaspRun)+"/"+vfile.at(i))) {
          outcar.GetPropertiesFile(dirrunname.at(idVaspRun)+"/"+vfile.at(i));
        }
      }
      //if(outcar.outcar=="")
      if(outcar.content=="") {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "The OUTCAR.static file in " << vaspRuns.at(idVaspRun).Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }

      AGL_data.energyinput.push_back(outcar.energy_cell);
      AGL_data.volumeinput.push_back(vaspRuns.at(idVaspRun).str.Volume());
      AGL_data.pressurecalculated.push_back(outcar.pressure_residual);
      stress_tensor = -outcar.stress;
      AGL_data.stresscalculated.push_back(stress_tensor);
      AGL_data.structurecalculated.push_back(idVaspRun);

      aurostd::StringstreamClean(aus);
      // Print out total energy, volume and calculated residual pressure
      // [OBSOLETE] aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(idVaspRun) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(AGL_data.energyinput.size()-1) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Volume (Ang^3) = " << AGL_data.volumeinput.at(AGL_data.volumeinput.size()-1) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Pressure (kB) = " << AGL_data.pressurecalculated.at(AGL_data.pressurecalculated.size()-1) << endl;      
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      outcar.clear();
    }
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::checkmin
// ***************************************************************************
namespace AGL_functions {
  //
  // checkmin: checks position of lowest energy in (E, V)
  // If lowest energy corresponds to smallest volume, sets cmerr = 1
  // If lowest energy corresponds to largest volume, sets cmerr = 2
  // Otherwise, sets cmerr = 0
  //
  uint checkmin(_AGL_data& AGL_data, int& cmerr, ofstream& FileMESSAGE) {
    ostringstream aus;
    vector<double> energy(AGL_data.energyinput.size()), volume(AGL_data.volumeinput.size());
    cmerr = 0;
    uint aglerror = 0;

    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      energy.at(i) = AGL_data.energyinput.at(i);
      volume.at(i) = AGL_data.volumeinput.at(i);
    }
    // Sorts (E, V) data into order of increasing volume
    //AGL_functions::qcksort (vol, idx, 0, is-1, FileMESSAGE);
    aglerror = AGL_functions::qcksortev (volume, energy, FileMESSAGE);
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V)" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aglerror;
    }
    // Find global minimum of energy data points
    double etref = energy.at(0);
    uint itmin = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      if(energy.at(i) < etref) {
        etref = energy.at(i);
        itmin = i;
      }
    }
    // Check that the minimum energy does not correspond to the largest or smallest volume
    // If it does, then this suggests that a larger or smaller volume may have a lower energy
    if(itmin == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Minimum energy is for smallest volume" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      cmerr = 1;
      return aglerror;
    } else if(itmin == (AGL_data.energyinput.size() - 1) ) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Minimum energy is for largest volume" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      cmerr = 2;
      return aglerror;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Energy minimum is contained in (E, V) data" << endl;  
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      cmerr = 0;
      return aglerror;
    }
    energy.clear();
    volume.clear();
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::checkconcav
// ***************************************************************************
namespace AGL_functions {
  //
  // checkconcav: checks concavity of (E, V) data
  // If data around global minimum of supplied (E, V) data is not concave, sets ccerr = 1
  // Otherwise, sets ccerr = 0
  //
  uint checkconcav(_AGL_data& AGL_data, int& ccerr, ofstream& FileMESSAGE) {
    ostringstream aus;
    uint j = 1;
    vector<double> energy(AGL_data.energyinput.size()), volume(AGL_data.volumeinput.size());
    ccerr = 0;
    uint aglerror = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      energy.at(i) = AGL_data.energyinput.at(i);
      volume.at(i) = AGL_data.volumeinput.at(i);
    }
    // Sort (E, V) data in order of increasing volume
    aglerror = AGL_functions::qcksortev (volume, energy, FileMESSAGE);
    if(aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V)" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return aglerror;
    }
    // Find global minimum of energy data points
    double etref = energy.at(0);
    uint itmin = 0;
    for (uint i = 0; i < AGL_data.energyinput.size(); i++) {
      if(energy.at(i) < etref) {
        etref = energy.at(i);
        itmin = i;
      }
    }
    // Finds first acceptable point (first point where E-V curve is concave)
    while ( ( (energy.at(j) - energy.at(j-1))/(volume.at(j) - volume.at(j-1)) >= (energy.at(j) - energy.at(j+1))/(volume.at(j) - volume.at(j+1)) ) && j < AGL_data.energyinput.size() - 2) {
      j = j + 1;
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Concavity check: j = " << j << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // If only the last three points are concave, then the (E, V) data cannot be used for GIBBS
    // A warning is given and the signal is given to rerun the VASP calculations with more k-points
    if(j == AGL_data.energyinput.size() - 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "All points show convex patterns" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    // If the first accepted point is already past the global minimum, the (E, V) could cause problems for GIBBS
    // Gives a warning and sends the signal to rerun the VASP calculations with more k-points
    if(j >= itmin) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "First concave point is already passed the global minimum" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    j = j + 1;
    uint jtmax = 2;

    // Point j marks the last accepted point, i the new trial point
    for (uint i=j+1; i <= AGL_data.energyinput.size()-1; i++) {
      if( (energy.at(j) - energy.at(j-1))/(volume.at(j) - volume.at(j-1)) < (energy.at(j) - energy.at(i))/(volume.at(j) - volume.at(i)) ) {
        j = j + 1;
        jtmax = i;
      }
    }

    // If the global minimum lies outside of the range of the accepted points, then there will be problems using the (E, V) data for GIBBS
    // Gives a warning and then gives the signal to rerun the VASP calculations with more k-points
    // Problems with noise in the (E, V) data are often caused by insufficient K-points or basis set
    if(jtmax <= itmin) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Global minimum of (E, V) data lies outside of initial range of points accepted for concavity" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    // If the data is concave and includes the global minimum, then data should be good for GIBBS
    // Gives success message and sends signal to continue to GIBBS method
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "(E, V) data is concave" << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    ccerr = 0;
    return aglerror;
  }
} // namespace AGL_functions

// ************************************************************************************************
//  This set of functions sorts the (E, V) data, fits it by a polynomial, and finds the minimum
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::qcksortev
// ***************************************************************************
namespace AGL_functions {
  //
  // Quick-sort algorithm: sorts the elements of array in ascending order
  // Calls aurostd::quicksort to actually sort data
  // Returns vol and energ sorted in order of increasing vol
  //
  uint qcksortev(vector<double>& vol, vector<double>& energ, ofstream& FileMESSAGE) {
    int icheck = 0;
    ostringstream aus;
    // Check if data is already in correct order
    for (uint i = 0; i < (vol.size()-1); i++) {
      if(vol.at(i+1) < vol.at(i)) {
        icheck = 1;
      }
    }
    // If data is already in correct order, exits function without sorting
    if(icheck == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "qcksort: Data is already arranged in increasing order" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    }
    // If data is not in the correct order, calls quicksort2 from aurostd library to sort it
    // Sorts both vol and energ so that they are in the order of increasing vol
    aurostd::sort (vol, energ);
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::qcksortevt
// ***************************************************************************
namespace AGL_functions {
  //
  // Quick-sort algorithm: sorts the elements of array in ascending order
  // Calls aurostd::quicksort to actually sort data
  // Returns vol, energ and tdebye sorted in order of increasing vol
  //
  uint qcksortevt(vector<double>& vol, vector<double>& energ, vector<double>& tdebye, ofstream& FileMESSAGE) {
    int icheck = 0;
    ostringstream aus;
    // Check if data is already in correct order
    for (uint i = 0; i < (vol.size()-1); i++) {
      if(vol.at(i+1) < vol.at(i)) {
        icheck = 1;
      }
    }
    // If data is already in correct order, exits function without sorting
    if(icheck == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "qcksort: Data is already arranged in increasing order" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return 0;
    }
    // If data is not in the correct order, calls quicksort2 from aurostd library to sort it
    // Sorts both vol and energ so that they are in the order of increasing vol
    aurostd::sort (vol, energ, tdebye);
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL set-up, extract, sort and check (E, V) data
// **************************************************************************

#endif  // _AFLOW_AGL_GET_EV_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
