// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to prepare VASP input files
// Stefano Curtarolo - 2007 Duke
// fixed for xz - 2008 (SC)

#ifndef _AFLOW_IVASP_CPP
#define _AFLOW_IVASP_CPP

#include "aflow.h"
#define _incarpad_ 48
#define _IVASP_DOUBLE2STRING_PRECISION_ 7
#define DIELECTRIC_DK 0.1
#define DEFAULT_EFIELD_PEAD 0.001

#define _DEBUG_IVASP_ false  //CO20190116

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INPUT
namespace KBIN {
  bool VASP_Produce_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::VASP_Produce_INPUT():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    bool Krun=TRUE;
    if(load_POSCAR_from_xvasp){
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp));     // produce POSCAR before KPOINTS  //CO20180420 - good for POCC
    } else {
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));     // produce POSCAR before KPOINTS
    }
    if(Krun) Krun=(Krun && KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
    if(LDEBUG) cerr << "VASP_Produce_INPUT: Checking if generate aflow.in only before calling VASP_Produce_POTCAR" << endl;  //CT20180719
    if (!XHOST.GENERATE_AFLOWIN_ONLY) { //CT20180719
      if(LDEBUG) cerr << "VASP_Produce_INPUT: Calling VASP_Produce_POTCAR" << endl;  //CT20180719
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_POTCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
      if(LDEBUG) {  //CO20181226 check species_pp patching (may cause problems with auto aflow.in generation)
        cerr << "VASP_Produce_INPUT: Checking species_pp patch" << endl;
        for(uint i=0;i<xvasp.str.species_pp.size();i++){
          cerr << "xstr.species[" << i << "]=" << xvasp.str.species[i] << endl;
          cerr << "xstr.species_pp[" << i << "]=" << xvasp.str.species_pp[i] << endl;
        }
      }
      convertPOSCARFormat(xvasp, kflags); //ME20190220
    } //CT20180719
    if(Krun) Krun=(Krun && KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
    return Krun;
  }
}

namespace KBIN {
  bool VASP_Modify_INPUT(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    bool Krun=TRUE;
    // if(Krun) Krun=(Krun && KBIN::VASP_Modify_POSCAR(xvasp,FileMESSAGE,aflags,vflags));  // moved up
    if(Krun) Krun=(Krun && KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags));
    if(LDEBUG) cerr << "VASP_Modify_INPUT: Checking if generate aflow.in only before calling VASP_Produce_POTCAR" << endl; //CT20180719
    if (!XHOST.GENERATE_AFLOWIN_ONLY) { //CT20180719
      if(LDEBUG) cerr << "VASP_Modify_INPUT: Calling VASP_Modify_POTCAR" << endl;  //CT20180719 
      if(Krun) Krun=(Krun && KBIN::VASP_Modify_POTCAR(xvasp,FileMESSAGE,aflags,vflags));
    } //CT20180719
    if(Krun) Krun=(Krun && KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags));
    return Krun;
  }
}

namespace KBIN {
  bool VASP_Produce_and_Modify_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp) {  //CO20180418
    bool Krun=TRUE;
    if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags,load_POSCAR_from_xvasp));
    if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
    return Krun;
  }
}


namespace KBIN {
  bool VASP_Write_ppAUID_FILE(const string& directory,const vector<string>& vppAUIDs,const vector<string>& species) {
    vector<string> vWRITE;
    vWRITE.push_back("AFLOW: AUID of pseudopotentials: NAUID="+aurostd::utype2string<uint>(vppAUIDs.size()));
    if(vppAUIDs.size()!=species.size()) {
      cout << "WARNING: KBIN::VASP_Write_ppAUID_FILE: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      cerr << "WARNING: KBIN::VASP_Write_ppAUID_FILE: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      return FALSE;
    }
    for(uint i=0;i<vppAUIDs.size();i++)
      vWRITE.push_back(vppAUIDs.at(i)+" "+aurostd::PaddedPOST(KBIN::VASP_PseudoPotential_CleanName(species.at(i)),3));
    aurostd::vectorstring2file(vWRITE,directory+"/"+DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT);
    return TRUE;
  }
  bool VASP_Write_ppAUID_FILE(const string& directory,const deque<string>& vppAUIDs,const deque<string>& species) {
    return VASP_Write_ppAUID_FILE(directory,aurostd::deque2vector(vppAUIDs),aurostd::deque2vector(species));
  }
  bool VASP_Write_ppAUID_AFLOWIN(const string& directory,const vector<string>& vppAUIDs,const vector<string>& species) {
    if(vppAUIDs.size()!=species.size()) {
      cout << "WARNING: KBIN::VASP_Write_ppAUID_AFLOWIN: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      cerr << "WARNING: KBIN::VASP_Write_ppAUID_AFLOWIN: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      return FALSE;
    }
    struct stat fileInfo;
    stat(string(directory+"/"+_AFLOWIN_).c_str(), &fileInfo);
    string date=std::ctime(&fileInfo.st_mtime);
    if (!date.empty() && date[date.length()-1] == '\n') date.erase(date.length()-1); // remove last newline

    string WRITE="[VASP_POTCAR_AUID]"+aurostd::joinWDelimiter(vppAUIDs,",");
    if(1||!aurostd::substring2bool(aurostd::file2string(directory+"/"+_AFLOWIN_),WRITE)){  //CO20210315 - better to write many times in case the file is run on another computer
      KBIN::AFLOWIN_ADD(directory+"/"+_AFLOWIN_,WRITE,"");
      //[CO20210315 - OBSOLETE]for(uint i=0;i<vppAUIDs.size();i++) {
      //[CO20210315 - OBSOLETE]  WRITE+=vppAUIDs.at(i);
      //[CO20210315 - OBSOLETE]  if(i<vppAUIDs.size()-1) WRITE+=",";
      //[CO20210315 - OBSOLETE]}
      //[CO20210315 - OBSOLETE]aurostd::execute("echo \""+WRITE+"\" >> "+directory+"/"+_AFLOWIN_);
      aurostd::execute("touch -m --date=\""+date+"\" "+directory+"/"+_AFLOWIN_);
    }
    return TRUE;
  }
  bool VASP_Write_ppAUID_AFLOWIN(const string& directory,const deque<string>& vppAUIDs,const deque<string>& species) {
    return VASP_Write_ppAUID_AFLOWIN(directory,aurostd::deque2vector(vppAUIDs),aurostd::deque2vector(species));
  }
}

namespace KBIN {
  bool VASP_Write_INPUT(_xvasp& xvasp,_vflags &vflags,const string& ext_module) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ifstream DirectoryStream;
    DirectoryStream.open(xvasp.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      ostringstream aus;
      aus << "XXXXX  MAKING DIRECTORY = " << xvasp.Directory << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+xvasp.Directory;
      system(str.c_str());
    }
    DirectoryStream.close();
    if(XHOST.AVOID_RUNNING_VASP){  //CO20200624
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"VASP should NOT be running",_INPUT_ILLEGAL_);  //better to throw to avoid VASP_Backup(), etc.
      //return false;
    }
    bool Krun=TRUE;
    // VASP VASP WRITE
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR")));
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR")));
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS")));
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.POTCAR,string(xvasp.Directory+"/POTCAR")));
    if(Krun) Krun=(Krun && VASP_Write_ppAUID_FILE(xvasp.Directory,xvasp.POTCAR_AUID,xvasp.str.species));
    if(Krun) Krun=(Krun && VASP_Write_ppAUID_AFLOWIN(xvasp.Directory,xvasp.POTCAR_AUID,xvasp.str.species));

    // VASP BACKUP VASP WRITE
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed"))  Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR_orig,string(xvasp.Directory+"/POSCAR.orig"+ext_module)));
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_INCAR_changed"))   Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR_orig,string(xvasp.Directory+"/INCAR.orig"+ext_module)));
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed")) Krun=(Krun && aurostd::stringstream2file(xvasp.KPOINTS_orig,string(xvasp.Directory+"/KPOINTS.orig"+ext_module)));
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_POTCAR_changed"))  Krun=(Krun && aurostd::stringstream2file(xvasp.POTCAR_orig,string(xvasp.Directory+"/POTCAR.orig"+ext_module)));
    // AS20210302 adding ext_module to be able to distinguish between different orig files (for example, we do not want APL to overwrite an existing POSCAR.orig file)

    if(vflags.KBIN_VASP_INCAR_VERBOSE) {;} // DUMMY

    return Krun;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// PseudoPotential_CleanName
// gets rid of all junk in the name
namespace KBIN {
  string VASP_PseudoPotential_CleanName(const string& speciesIN) {return VASP_PseudoPotential_CleanName_20190712(speciesIN);} //CO20190712
  string VASP_PseudoPotential_CleanName_20190712(const string& speciesIN) { //CO20190712
    //the old function assumed only a single species input
    //now, the species input can be a full species string: Mn_pvPt, etc.
    //need to remove ALL instances of pp info
    //no longer need for loops
    //also, this new function leverages InPlace substitution, instead of creating new strings every time (FASTER)

    string species=speciesIN;
    aurostd::VASP_PseudoPotential_CleanName_InPlace(species);
    return species;
  }
  string VASP_PseudoPotential_CleanName_20190101(const string& speciesIN) {
    string species=speciesIN;
    uint i=1,imax=2;
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_GW");  //CO20190712
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_new");  //CO20190712
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_AE");  //CO20190712
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_sv");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_pv");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_d2");  //CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.05May2010/As_d2_GW
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_d");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_soft");  //CO20190712 - BEFORE _s
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_s");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_200eV");
    //[CO20190712 - moved up before _s]for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_soft");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_2_n");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_h");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_1");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_2");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_3");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.75"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.75
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.66
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H1.33
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.25"); //CO20190712 - before all other decimal numbers
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.5"); //CO20190712 - potpaw_PBE/potpaw_PBE.06May2010/H1.5
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".75");  //CO20190712 - before 0.5
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".25");  //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.25
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.66
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.33
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".42"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.42
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".58"); //CO20190712 - before 0.5 //potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.58
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".5");
    //[CO20190712 - moved up 0.5]for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".75");
    //[CO20190712 - moved up before all other decimal numbers]for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.25");

    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+1");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+3");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+5");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+7");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-1");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-3");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-5");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-7");
    //  for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1");
    //CO START
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POT_LDA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POT_GGA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POT_PBE+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN+"/");
    //CO END
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"__"); //CO20190712 - BEFORE _ - potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW__
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_");  //CO20190712  //potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW_
    return species;
  }
}

namespace KBIN {
  bool VASP_PseudoPotential_CleanName_TEST(void){ //CO20190712
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    stringstream message;

    ostream& oss=cout;
    ofstream FileMESSAGE;
    _aflags aflags; aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");

    string compound_to_test="";
    vector<string> velements;
    vector<double> vcomposition;

    compound_to_test="W_sv_GWCs_sv_GWCr_sv_GWPt_ZORALi_AE_GW2";
    bool keep_pp=false;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="CrCsLiPtW"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="W_sv_GWCs_sv_GWCr_sv_GWPt_ZORALi_AE_GW2";
    keep_pp=true;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="Cr_sv_GWCs_sv_GWLi_AE_GW2Pt_ZORAW_sv_GW"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="H1.25";
    keep_pp=false;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="H"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="H1.25";
    keep_pp=true;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="H1.25"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="Mn.5W5H1.25";
    keep_pp=false;
    velements=aurostd::getElements(compound_to_test,vcomposition,composition_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << ", vcomposition=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vcomposition,5),",") << endl;}
    if(!(aurostd::joinWDelimiter(velements,"")=="HMnW" && vcomposition.size()==3 && aurostd::isequal(vcomposition[0],1.25) && aurostd::isequal(vcomposition[1],0.5) && aurostd::isequal(vcomposition[2],5.0) )){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    //return true;

    //full test of pp available
    vector<string> directories_2_search;
    vector<string> potcar_paths;
    vector<string> path_parts;
    vector<string> element_parts,_element_parts;
    string element_raw="",element_clean="";
    string element="",suffix="";
    aurostd::string2tokens(DEFAULT_VASP_POTCAR_DIRECTORIES,directories_2_search,",");
    for(uint i=0;i<directories_2_search.size();i++){
      if(aurostd::IsDirectory(directories_2_search[i])){
        const string& dir=directories_2_search[i];
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " searching dir=" << dir << endl;}
        aurostd::string2vectorstring(aurostd::execute2string("find "+dir+" -name POTCAR*"),potcar_paths);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " potcar_paths=" << aurostd::joinWDelimiter(potcar_paths," ") << endl;}
        for(uint j=0;j<potcar_paths.size();j++){
          if(!aurostd::RemoveWhiteSpacesFromTheBack(potcar_paths[j]).empty()){
            const string& potcar_path=potcar_paths[j];
            aurostd::string2tokens(potcar_path,path_parts,"/");
            if(aurostd::WithinList(path_parts,"TESTS")){continue;}  //skip VASP/TESTS
            if(aurostd::WithinList(path_parts,"BENCHS")){continue;}  //skip VASP/BENCHS
            if(path_parts.size()<2){continue;}
            element_raw=path_parts[path_parts.size()-2];
            if(aurostd::substring2bool(element_raw,"runelements")){continue;} //garbage pp
            if(aurostd::substring2bool(element_raw,"Free")){continue;} //garbage pp
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " element_raw=" << element_raw << endl;}
            aurostd::string2tokens(element_raw,element_parts,"_");
            if(element_parts.size()==0){continue;}
            element=suffix="";
            if(element_parts.size()>1){
              _element_parts.clear();
              for(uint k=1;k<element_parts.size();k++){_element_parts.push_back(element_parts[k]);} //skip element
              suffix="_"+aurostd::joinWDelimiter(_element_parts,"_");
            }
            element=element_parts[0];
            element=KBIN::VASP_PseudoPotential_CleanName(element);  //clean H1.25
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " element=\"" << element << "\", suffix=\"" << suffix << "\"" << endl;}
            element_clean=KBIN::VASP_PseudoPotential_CleanName(element_raw);
            if(element!=element_clean){ //idempotent
              message << "\"" << element << "\" != \"" << element_clean << "\" == KBIN::VASP_PseudoPotential_CleanName(\"" << element_raw << "\")";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
              return false;
            }
            if(GetAtomNumber(element)==0){
              message << "(GetAtomNumber(\"" << element << "\")==0 [element_raw=\"" << element_raw << "\",element_clean=\"" << element_clean << "\"" << "]";
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
              return false;
            }
          }
        }
      }
    }
    return true;
  }
}

namespace KBIN {
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX) {
    return XATOM_SplitAlloySpecies(alloy_in,speciesX);
  }
}

namespace KBIN {
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX) {
    return XATOM_SplitAlloySpecies(alloy_in,speciesX,natomsX);
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(string alloy_in, string &speciesA, string &speciesB) {
    string alloy=KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(alloy_in));
    if(0) {
      int i=0;
      speciesA=alloy[i++];
      if(alloy[i]>='a' && alloy[i]<='z') speciesA=speciesA+alloy[i++];
      speciesB=alloy[i++];
      if(alloy[i]>='a' && alloy[i]<='z') speciesB=speciesB+alloy[i++];
    }
    if(1) {
      speciesA="";speciesB="";
      int speciesN=0;
      for(uint i=0;i<alloy.length();i++) {
        if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
        if(speciesN==1) speciesA+=alloy[i];
        if(speciesN==2) speciesB+=alloy[i];
      }
    }
    speciesA=aurostd::CleanStringASCII(speciesA);
    speciesB=aurostd::CleanStringASCII(speciesB);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(string alloy_in, string &speciesA, string &speciesB, string &speciesC) {
    string alloy=KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(alloy_in));
    speciesA="";speciesB="";speciesC="";
    int speciesN=0;
    for(uint i=0;i<alloy.length();i++) {
      if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
      if(speciesN==1) speciesA+=alloy[i];
      if(speciesN==2) speciesB+=alloy[i];
      if(speciesN==3) speciesC+=alloy[i];
    }
    speciesA=aurostd::CleanStringASCII(speciesA);
    speciesB=aurostd::CleanStringASCII(speciesB);
    speciesC=aurostd::CleanStringASCII(speciesC);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB) {
    for (uint k=0;k<alloy.size();k++) {
      speciesA.push_back("");speciesB.push_back("");
      KBIN::VASP_SplitAlloySpecies(alloy.at(k),speciesA.at(k),speciesB.at(k));
    }
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB, vector<string> &speciesC) {
    for (uint k=0;k<alloy.size();k++) {
      speciesA.push_back("");speciesB.push_back("");
      KBIN::VASP_SplitAlloySpecies(alloy.at(k),speciesA.at(k),speciesB.at(k),speciesC.at(k));
    }
    return TRUE;
  }
}

namespace KBIN {
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX) {
    return XATOM_SplitAlloyPseudoPotentials(alloy_in,species_ppX);
  }
}

namespace KBIN {
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX) {
    return XATOM_SplitAlloyPseudoPotentials(alloy_in,species_ppX,natomsX);
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB) {
    if(0) {
      uint i=0;
      species_ppA=alloy[i++];
      while((alloy[i]>='Z' || alloy[i]<='A') && i<=alloy.length()) species_ppA=species_ppA+alloy[i++];
      species_ppB=alloy[i++];
      while((alloy[i]>='Z' || alloy[i]<='A') && i<=alloy.length()) species_ppB=species_ppB+alloy[i++];
    }
    if(1) {
      species_ppA="";species_ppB="";
      int speciesN=0;
      for(uint i=0;i<alloy.length();i++) {
        if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
        if(speciesN==1) species_ppA+=alloy[i];
        if(speciesN==2) species_ppB+=alloy[i];
      }
    }
    species_ppA=aurostd::CleanStringASCII(species_ppA);
    species_ppB=aurostd::CleanStringASCII(species_ppB);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB, string &species_ppC) {
    species_ppA="";species_ppB="";species_ppC="";
    int speciesN=0;
    for(uint i=0;i<alloy.length();i++) {
      if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
      if(speciesN==1) species_ppA+=alloy[i];
      if(speciesN==2) species_ppB+=alloy[i];
      if(speciesN==3) species_ppC+=alloy[i];
    }
    species_ppA=aurostd::CleanStringASCII(species_ppA);
    species_ppB=aurostd::CleanStringASCII(species_ppB);
    species_ppC=aurostd::CleanStringASCII(species_ppC);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &pseudosA, vector<string> &pseudosB) {
    for (uint k=0;k<alloy.size();k++) {
      pseudosA.push_back("");pseudosB.push_back("");
      KBIN::VASP_SplitAlloyPseudoPotentials(alloy.at(k),pseudosA.at(k),pseudosB.at(k));
    }
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &pseudosA, vector<string> &pseudosB, vector<string> &pseudosC) {
    for (uint k=0;k<alloy.size();k++) {
      pseudosA.push_back("");pseudosB.push_back("");
      KBIN::VASP_SplitAlloyPseudoPotentials(alloy.at(k),pseudosA.at(k),pseudosB.at(k),pseudosC.at(k));
    }
    return TRUE;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INCAR
namespace KBIN {
  bool VASP_Produce_INCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"empty AflowIn",_FILE_NOT_FOUND_);}
    //[CO20200624 - OBSOLETE]if(!kflags.AFLOW_MODE_VASP) {cerr << __AFLOW_FUNC__ << " should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    if(!kflags.AFLOW_MODE_VASP) {pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"kflags.AFLOW_MODE_VASP NOT set",std::cout,_LOGGER_WARNING_);}
    ostringstream aus;
    bool Krun=TRUE;
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.NCPUS=0;
    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",FALSE);

    aus << "00000  MESSAGE INCAR   generation in " << xvasp.Directory << "  " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    bool KBIN_VASP_INCAR_MODE_EMPTY=!vflags.KBIN_VASP_INCAR_MODE.flag("IMPLICIT") && !vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT") && !vflags.KBIN_VASP_INCAR_MODE.flag("EXTERNAL");

    // IMPLICIT or EXPLICIT or EXTERNAL for INCAR
    Krun=(Krun && (vflags.KBIN_VASP_INCAR_MODE.flag("IMPLICIT") ||
          vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT") ||
          vflags.KBIN_VASP_INCAR_MODE.flag("EXTERNAL") || KBIN_VASP_INCAR_MODE_EMPTY));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_INCAR_MODE_IMPLICIT] or [VASP_INCAR_MODE_EXPLICIT] or [VASP_INCAR_MODE_EXTERNAL] must be specified " << Message(__AFLOW_FILE__,aflags) << endl; //ME20181113
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // EMPTY ************************************************** INCAR
    if(Krun && KBIN_VASP_INCAR_MODE_EMPTY) {  // [VASP_INCAR_MODE_EMPTY] construction
      aus << "00000  MESSAGE INCAR   generation EMPTY file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.INCAR << "#AFLOW INCAR automatically generated (INCAR_MODE_EMPTY)" << endl; //CO20200624 - FIRST LINE!
    }
    // IMPLICIT ************************************************** INCAR
    if(Krun && vflags.KBIN_VASP_INCAR_MODE.flag("IMPLICIT")) {  // [VASP_INCAR_MODE_IMPLICIT] construction
      aus << "00000  MESSAGE INCAR   generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.INCAR << "#AFLOW INCAR automatically generated (INCAR_MODE_IMPLICIT)" << endl; //CO20200624 - FIRST LINE!
      if(vflags.KBIN_VASP_INCAR_FILE.flag("SYSTEM_AUTO")) {
        xvasp.INCAR << "SYSTEM=" << xvasp.str.title << endl;
        xvasp.INCAR << "#PROTOTYPE=" << xvasp.str.prototype << endl;
        xvasp.INCAR << "#INFO=" << xvasp.str.info << endl;
        // if(LDEBUG) cerr << xvasp.INCAR.str() << endl;
      } else if (vflags.AFLOW_SYSTEM.isentry) {  //ME20181121
        xvasp.INCAR << "SYSTEM=" << vflags.AFLOW_SYSTEM.content_string << std::endl;
      }
    }
    // EXPLICIT ************************************************** INCAR
    if(Krun && vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT")) {  // [VASP_INCAR_MODE_EXPLICIT] construction
      if(vflags.KBIN_VASP_INCAR_FILE.flag("KEYWORD") && !vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE INCAR   generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.INCAR << "#AFLOW INCAR automatically generated (INCAR_MODE_EXPLICIT)" << endl; //CO20200624 - FIRST LINE!
        xvasp.INCAR << vflags.KBIN_VASP_INCAR_EXPLICIT.str(); //ME20181113
        // [OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.INCAR,"[VASP_INCAR_FILE]");
        // [ME20181226 OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.INCAR,"[VASP_INCAR_FILE]");
      } else if(!vflags.KBIN_VASP_INCAR_FILE.flag("KEYWORD") && vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE INCAR   generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.INCAR << "#AFLOW INCAR automatically generated (INCAR_MODE_EXPLICIT_START_STOP)" << endl; //CO20200624 - FIRST LINE!
        xvasp.INCAR << vflags.KBIN_VASP_INCAR_EXPLICIT_START_STOP.str();
        // [ME20181226 OBSOLETE] if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]STOP"))
        // [OBSOLETE]	aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.INCAR,"[VASP_INCAR_MODE_EXPLICIT]START","[VASP_INCAR_MODE_EXPLICIT]STOP");
        // [ME20181226 OBSOLETE] aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.INCAR,"[VASP_INCAR_MODE_EXPLICIT]START","[VASP_INCAR_MODE_EXPLICIT]STOP");
      } else {
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT] do not confuse aflow !!" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT] Possible modes " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] INCAR EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_INCAR_FILE]SYSTEM=AuTi.274_LDA" << endl;
        aus << "[VASP_INCAR_FILE]#Prototype Ni2In" << endl;
        aus << "[VASP_INCAR_FILE]PREC=med" << endl;
        aus << "[VASP_INCAR_FILE]ISMEAR=1" << endl;
        aus << "[VASP_INCAR_FILE]SIGMA=0.2" << endl;
        aus << "[VASP_INCAR_FILE]IBRION=2" << endl;
        aus << "[VASP_INCAR_FILE]NSW=51" << endl;
        aus << "[VASP_INCAR_FILE]ISIF=3" << endl;
        aus << "[VASP_INCAR_FILE]ENMAX=333.666" << endl;
        aus << "[VASP_INCAR_FILE]NBANDS=47" << endl;
        aus << "[VASP_INCAR_FILE]MAGMOM=  5 5 5 5 5 5" << endl;
        aus << "[VASP_INCAR_FILE]ISPIND=2" << endl;
        aus << "[VASP_INCAR_FILE]ISPIN=2" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] INCAR EXPLICIT MODE with START/STOP" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]START" << endl;
        aus << "SYSTEM=AuTi.274_LDA" << endl;
        aus << "#Prototype Ni2In" << endl;
        aus << "PREC=med" << endl;
        aus << "ISMEAR=1" << endl;
        aus << "SIGMA=0.2" << endl;
        aus << "IBRION=2" << endl;
        aus << "NSW=51" << endl;
        aus << "ISIF=3" << endl;
        aus << "ENMAX=333.666" << endl;
        aus << "NBANDS=47" << endl;
        aus << "MAGMOM=  5 5 5 5 5 5" << endl;
        aus << "ISPIND=2" << endl;
        aus << "ISPIN=2" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT] Note " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT]START must be present and no [VASP_INCAR_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT]STOP  must be present and no [VASP_INCAR_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  or [VASP_INCAR_FILE] present and NO START/STOP" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_INCAR_MODE.flag("EXTERNAL")) {  // [VASP_INCAR_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE INCAR   generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.INCAR << "#AFLOW INCAR automatically generated (INCAR_MODE_EXTERNAL)" << endl; //CO20200624 - FIRST LINE!
      if(vflags.KBIN_VASP_INCAR_FILE.flag("COMMAND") && vflags.KBIN_VASP_INCAR_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_INCAR_MODE]FILE=  and  [VASP_INCAR_MODE]COMMAND=  can not be used together " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_INCAR_FILE.flag("COMMAND") && (vflags.KBIN_VASP_INCAR_FILE.flag("FILE") || !vflags.KBIN_VASP_INCAR_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_INCAR_FILE.flag("FILE")) {
          // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]FILE=",1,TRUE);
          file = vflags.KBIN_VASP_INCAR_FILE.getattachedscheme("FILE"); //ME20181113
          aus << "00000  MESSAGE INCAR   generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_INCAR;
          aus << "00000  MESSAGE INCAR   generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_INCAR << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR INCAR file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR INCAR file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.INCAR << aurostd::file2string(file);
      }
      if(vflags.KBIN_VASP_INCAR_FILE.flag("COMMAND") && !vflags.KBIN_VASP_INCAR_FILE.flag("FILE")) {
        // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]COMMAND=",1,FALSE);
        file = vflags.KBIN_VASP_INCAR_FILE.getattachedscheme("COMMAND"); //ME20181113
        aus << "00000  MESSAGE INCAR   generation from command= '" << file << "'" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_INCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp  //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_INCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name  //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR INCAR file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR INCAR file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.INCAR << aurostd::file2string(file);       // load INCAR
        aurostd::RemoveFile("./_aflow_INCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp //CO20200502 - threadID
      }
    }
    // INCAR DONE **************************************************
    xvasp.INCAR << "#incar" << endl;
    VASP_CleanUp_INCAR(xvasp);  //CO20200624
    xvasp.INCAR_orig << xvasp.INCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);
    return Krun;
  };  // KBIN::VASP_Produce_INCAR
}

namespace KBIN {
  bool VASP_Modify_INCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;

    // bool vflags.KBIN_VASP_INCAR_VERBOSE=TRUE;
    if(Krun && kflags.KBIN_MPI) {
      xvasp.NCPUS=kflags.KBIN_MPI_NCPUS;
      if(kflags.KBIN_MPI_AUTOTUNE) {
        aus << "00000  MESSAGE INCAR-MPI: found AUTOTUNE option " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE INCAR-MPI: input files WILL be auto-tuned for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_MPI_Autotune(xvasp,aflags,vflags.KBIN_VASP_INCAR_VERBOSE);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      } else {
        aus << "00000  MESSAGE INCAR-MPI: AUTOTUNE option NOT found! (aflow_ivasp.cpp) " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE INCAR-MPI: input files MUST be appropriate for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO.isentry) {                                                        /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SYSTEM_AUTO - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_System_Auto(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    } else if (vflags.AFLOW_SYSTEM.isentry) {  //ME20181121
      xvasp.INCAR << "SYSTEM=" << vflags.AFLOW_SYSTEM.content_string << std::endl;
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // if(Krun && !vflags.KBIN_VASP_RUN.flag("STATIC") && !vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC") && ! vflags.KBIN_VASP_RUN_RELAX_STATIC_PATCH_STATIC) //o change for relax if STATIC !!!
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry && 
        (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL")               || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS") ||
         vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")        || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME") ||
         vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))) {      /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL"))                 aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_ALL - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS"))                aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_IONS - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE"))          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_CELL_SHAPE - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME"))         aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_CELL_VOLUME - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))    aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_IONS_CELL_VOLUME - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.isentry) {
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="ENERGY")         aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=ENERGY (default) - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES")         aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=FORCES - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="ENERGY_FORCES")  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=ENERGY_FORCES - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES_ENERGY")  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=FORCES_ENERGY - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT RELAX_MODE=" << "ENERGY" << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Relax_ON(xvasp,vflags,1); // relax number (start)
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]EDIFFG=" << vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.content_double << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("EDIFFG",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry) {             /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NBANDS - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("NBANDS",xvasp,vflags,"",KBIN::XVASP_INCAR_GetNBANDS(xvasp,aflags,vflags.KBIN_VASP_FORCE_OPTION_SPIN.option),0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry && vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int>0) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NBANDS=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("NBANDS",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]PSTRESS=" << vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.content_double << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("PSTRESS",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry && vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double>0) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]POTIM=" << vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("POTIM",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // SPIN
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SPIN=" << (vflags.KBIN_VASP_FORCE_OPTION_SPIN.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        //      aus << "00000  MESSAGE-DEFAULT SPIN=" << (DEFAULT_VASP_FORCE_OPTION_SPIN?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE-DEFAULT SPIN=" << "NEGLECT" << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) {
        KBIN::XVASP_INCAR_PREPARE_GENERIC("SPIN",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_SPIN.option);  // CHANGE ONLY IF SPIN IS MENTIONED
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SPIN=REMOVE_RELAX_1 - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SPIN=REMOVE_RELAX_2 - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    // LSCOUPLING must be after spin
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LSCOUPLING=" << (vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT LSCOUPLING=" << (DEFAULT_VASP_FORCE_OPTION_LSCOUPLING?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.isentry || DEFAULT_VASP_FORCE_OPTION_LSCOUPLING) {
        KBIN::XVASP_INCAR_PREPARE_GENERIC("LS_COUPLING",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // AUTO_MAGMOM must be after LSCOUPLING, AFTER SPIN
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry&&vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {   //CO, MAGMOM should only be on if SPIN is SPECIFIED and ON (as default is to neglect)
        if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry) {
          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]AUTO_MAGMOM=" << (vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        } else {
          aus << "00000  MESSAGE-DEFAULT AUTO_MAGMOM=" << (DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        }
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if((vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry && vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option) || DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM) {
          KBIN::XVASP_INCAR_PREPARE_GENERIC("AUTO_MAGMOM",xvasp,vflags,"",0,0.0,FALSE);
          xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
        }
      }
    }

    // BADER 
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]BADER=" << (vflags.KBIN_VASP_FORCE_OPTION_BADER.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT BADER=" << (DEFAULT_VASP_FORCE_OPTION_BADER?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry || DEFAULT_VASP_FORCE_OPTION_BADER) KBIN::XVASP_INCAR_BADER(xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_BADER.option); WILL BE DONE WHEN NEEDED
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // ELF 
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ELF=" << (vflags.KBIN_VASP_FORCE_OPTION_ELF.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT ELF=" << (DEFAULT_VASP_FORCE_OPTION_ELF?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry || DEFAULT_VASP_FORCE_OPTION_ELF) KBIN::XVASP_INCAR_ELF(xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_ELF.option); WILL BE DONE WHEN NEEDED
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20200624 - START
    // NELM_EQUAL
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.isentry && xvasp.NRELAXING<=xvasp.NRELAX) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NELM=" << vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("NELM",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.content_int,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // NELM_STATIC_EQUAL
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.isentry && xvasp.NRELAXING>xvasp.NRELAX) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NELM_STATIC=" << vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("NELM_STATIC",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.content_int,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }
    //CO20200624 - END

    // NSW_EQUAL
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL) {           /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE>0) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NSW=" << vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("NSW",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE,0.0,FALSE);	
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // SYM
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_SYM.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SYM=" << (vflags.KBIN_VASP_FORCE_OPTION_SYM.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT SYM=" << (DEFAULT_VASP_FORCE_OPTION_SYM?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SYM",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_SYM.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // WAVECAR
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]WAVECAR=" << (vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT WAVECAR=" << (DEFAULT_VASP_FORCE_OPTION_WAVECAR?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("WAVECAR",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // CHGCAR
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CHGCAR=" << (vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT CHGCAR=" << (DEFAULT_VASP_FORCE_OPTION_CHGCAR?"ON":"OFF") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("LCHARG",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //ME20191028
    if (Krun && vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry) {
      string chgcar = vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.content_string;
      if (chgcar[0] != '/') chgcar = aurostd::CleanFileName(aflags.Directory + "/" + chgcar);  // relative path
      if (aurostd::FileExist(chgcar)) {
        if (aurostd::IsCompressed(chgcar)) {
          string ext = aurostd::GetCompressionExtension(chgcar);
          aurostd::CopyFile(chgcar, aflags.Directory + "/CHGCAR" + ext);  // relative path
          aurostd::UncompressFile(aflags.Directory + "/CHGCAR" + ext);
        } else {
          aurostd::CopyFile(chgcar, aflags.Directory + "/CHGCAR");
        }
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("ICHARG", xvasp, vflags, vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.content_string, 1, 0.0, true);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      } else {
        string message = "Cannot use CHGCAR file " + chgcar + ". File not found.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }
    }

    // LDAU0
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU0.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU=OFF - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_LDAU_OFF(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // LDAU1
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU1=ON - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_LDAU_SPECIES!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_SPECIES=\"" << vflags.KBIN_VASP_LDAU_SPECIES << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_LDAU_PARAMETERS!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_PARAMETERS=\"" << vflags.KBIN_VASP_LDAU_PARAMETERS << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_AFLOW_AUTO_flag=\"" << vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_LDAU_ON(xvasp,vflags,1);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // LDAU2
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU2=ON - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_LDAU_SPECIES!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_SPECIES=\"" << vflags.KBIN_VASP_LDAU_SPECIES << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_LDAU_PARAMETERS!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_PARAMETERS=\"" << vflags.KBIN_VASP_LDAU_PARAMETERS << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_AFLOW_AUTO_flag=\"" << vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_LDAU_ON(xvasp,vflags,2);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // LDAU_ADIABATIC
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU=ADIABATIC (steps=" << vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int << ") - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  KBIN::XVASP_INCAR_LDAU_ADIABATIC_OFF(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE); // NOTHING TO CHANGE, IT IS ON THE FLY
      //  xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);                                                  // NOTHING TO CHANGE, IT IS ON THE FLY
    }

    // LDAU_CUTOFF
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry) {          /*************** INCAR **************/
      //    cerr << "vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU=CUTOFF - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  KBIN::XVASP_INCAR_LDAU_CUTOFF_OFF(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE); // NOTHING TO CHANGE, IT IS ON THE FLY
      //  xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);                                                  // NOTHING TO CHANGE, IT IS ON THE FLY
    }

    // PREC
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry==TRUE) {                                                         /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]PREC=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved)   aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]PREC_preserved - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT PREC=" << DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Precision(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.isentry
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.isentry) {             /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ENMAX_MULTIPLY="<< vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.content_double << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ENMAX_MULTIPLY",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // ALGO
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      //ME20181108 - PREC=PHONONS overrides ALGO
      if (vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme == "PHONONS") {
        aus << "00000  MESSAGE-OPTION ALGO = NORMAL (override by PREC=PHONONS) - " << Message(__AFLOW_FILE__,aflags) << endl;
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme, "NORMAL");
      } else if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ALGO=" << vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved)   aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ALGO_PRESERVED - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT ALGO=" << DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ALGO",xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme,0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // METAGGA
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]METAGGA=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << " - " << Message(__AFLOW_FILE__,aflags) << endl;;
      } else {
        aus << "00000  MESSAGE-DEFAULT METAGGA=" << DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Metagga(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // IVDW
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]IVDW=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << " - " << Message(__AFLOW_FILE__,aflags) << endl;;
      } else {
        aus << "00000  MESSAGE-DEFAULT IVDW=" << DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Ivdw(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // ABMIX
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ABMIX=" << vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        // the rest is neglected... no AMIX BMIX AMIX_MAG BMIX_MAG
      } else {
        // neglected     aus << "00000  MESSAGE-DEFAULT ABMIX=" << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_ABMIX(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // TYPE
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry) {                                                              /*************** INCAR **************/
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='D') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='M') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='S') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='I') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(__AFLOW_FILE__,aflags) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT TYPE=" << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("TYPE",xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme,0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.isentry && xvasp.NRELAXING<=xvasp.NRELAX) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ISMEAR=" << vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ISMEAR",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.content_int
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.isentry && xvasp.NRELAXING==xvasp.NRELAX+1) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ISMEAR_STATIC=" << vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ISMEAR_STATIC",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.content_int,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.content_int
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.isentry && xvasp.NRELAXING>xvasp.NRELAX+1) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ISMEAR_BANDS=" << vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ISMEAR_BANDS",xvasp,vflags,"",vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.content_int,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.isentry && xvasp.NRELAXING<=xvasp.NRELAX) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SIGMA=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_) << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SIGMA",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.isentry && xvasp.NRELAXING==xvasp.NRELAX+1) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SIGMA_STATIC=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_) << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SIGMA_STATIC",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.isentry && xvasp.NRELAXING>xvasp.NRELAX+1) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SIGMA_BANDS=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_) << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SIGMA_BANDS",xvasp,vflags,"",0,vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.content_double,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // CONVERT_UNIT_CELL
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry) {                                                              /*************** INCAR **************/
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=STANDARD_PRIMITIVE - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=STANDARD_CONVENTIONAL - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=NIGGLI - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=MINKOWSKI - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=INCELL - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=COMPACT - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=WIGNERSEITZ - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=CARTESIAN - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=FRACTIONAL - "<< Message(__AFLOW_FILE__,aflags) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=PRESERVE - "<< Message(__AFLOW_FILE__,aflags) << endl; //CO
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // print the AFIX
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.isentry) {                                                    /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]IGNORE_AFIX=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.content_string << "  - " << Message(__AFLOW_FILE__,aflags) << endl;
      for(uint i=0;i<vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme.size();i++) 
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]IGNORE_AFIX=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme.at(i) << "  - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // DO THE PAW CORRECTIONS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && 0) {                                                    /*************** INCAR **************/
      if(xvasp.POTCAR_PAW==TRUE) {
        aus << "00000  MESSAGE-DEFAULT PAW_CORRECTIONS" << " - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("LASPH",xvasp,vflags,"",0,0.0,ON);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // STATIC at the end after all the relax stuff has been added
    if(Krun && vflags.KBIN_VASP_RUN.flag("STATIC")) {                                                                      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_RUN_STATIC] - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Static_ON(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC")) {                                                             /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]STATIC - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Static_ON(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // these functions are forced in aflow_kvasp.cpp inside the RELAX_STATIC_BANDS run
    // if(Krun && vflags.KBIN_VASP_RUN_RELAX_STATIC_PATCH_STATIC) {                                                   /*************** INCAR **************/
    // aus << "00000  MESSAGE-OPTION  [VASP_RUN_RELAX_STATIC] - " << Message(__AFLOW_FILE__,aflags) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // KBIN::XVASP_INCAR_Relax_Static_ON(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);
    // xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    // }
    // if(Krun && vflags.KBIN_VASP_RUN_RELAX_STATIC_BANDS_PATCH_STATIC_BANDS) {                                       /*************** INCAR **************/
    // aus << "00000  MESSAGE-OPTION  [VASP_RUN_RELAX_STATIC_BANDS] - " << Message(__AFLOW_FILE__,aflags) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // KBIN::XVASP_INCAR_Bands_ON(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);     // FIX
    // xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    // }


    // INTERCEPT ERRORS AND WRAP UP after all the incar has been prepared
    if(Krun && aurostd::substring2bool(xvasp.INCAR,"LEPSILON=.TRUE.",true)) {  /*************** INCAR **************/  //CO20210315
      //ME20191205 - also remove NCORE
      aus << "00000  MESSAGE REMOVE ENTRIES NPAR and NCORE because of LEPSILON - " << Message(__AFLOW_FILE__,aflags) << endl;  //ME20191205
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR,NCORE","LEPSILON=.TRUE.",vflags.KBIN_VASP_INCAR_VERBOSE);  //ME20191205 //CO20210313
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //repetita iuvant
    }
    if(Krun && aurostd::substring2bool(xvasp.INCAR,"LCALCEPS=.TRUE.",true)) {  /*************** INCAR **************/  //CO20210315
      //ME20191205 - also remove NCORE
      aus << "00000  MESSAGE REMOVE ENTRIES NPAR and NCORE because of LCALCEPS - " << Message(__AFLOW_FILE__,aflags) << endl;  //ME20191205
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR,NCORE","LCALCEPS=.TRUE.",vflags.KBIN_VASP_INCAR_VERBOSE);  //ME20191205
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //repetita iuvant
    } // KEVIN
    if(Krun && aurostd::substring2bool(xvasp.INCAR,"IBRION=",true)) {  /*************** INCAR **************/  //CO20210315
      uint IBRION=aurostd::kvpair2utype<uint>(xvasp.INCAR,"IBRION","=");  //CO20210315
      if(IBRION==8) {
        aus << "00000  MESSAGE REMOVE ENTRIES NPAR and NCORE because of IBRION=" << IBRION << "  - " << Message(__AFLOW_FILE__,aflags) << endl;  //ME20191205
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR,NCORE","IBRION=8",vflags.KBIN_VASP_INCAR_VERBOSE);  //ME20191205 //CO20210313
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //repetita iuvant
      }
    }

    // ------------------------------------
    // end
    VASP_CleanUp_INCAR(xvasp);  //CO20200624
    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);
    return Krun;
  };  // KBIN::VASP_Modify_INCAR
}

namespace KBIN {
  void VASP_CleanUp_INCAR(_xvasp& xvasp){ //CO20200624
    //remove duplicate lines - round 1 START
    vector<string> vlines_orig,vlines;
    aurostd::stream2vectorstring(xvasp.INCAR,vlines_orig);
    bool found=false;
    for(uint iline1=0;iline1<vlines_orig.size();iline1++){
      found=false;
      for(uint iline2=0;iline2<vlines.size()&&found==false;iline2++){
        if(vlines_orig[iline1]==vlines[iline2]){found=true;} //find STRICT duplicate lines first
      }
      if(found){continue;}
      vlines.push_back(vlines_orig[iline1]);
    }
    //remove duplicate lines - round 1 END
    //ADD MORE CLEANING HERE

    //write out clean INCAR
    aurostd::StringstreamClean(xvasp.INCAR);
    //[need additional endl]xvasp.INCAR << aurostd::joinWDelimiter(vlines,"\n");
    for(uint iline=0;iline<vlines.size();iline++){
      if(iline==vlines.size()-1 && aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vlines[iline]).empty()){continue;} //remove last extra line, without this line, clean will add newlines everytime (growing INCAR unnecessarily)
      xvasp.INCAR << vlines[iline] << endl;
    }
  }
}

namespace KBIN {
  bool VASP_Reread_INCAR(_xvasp& xvasp) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    ofstream FileMESSAGE;
    _aflags aflags;aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");
    return VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags);
  }
  bool VASP_Reread_INCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    if(!aurostd::FileExist(xvasp.Directory+"/INCAR")) {
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"INCAR not present in directory: "+xvasp.Directory,aflags,FileMESSAGE,std::cout,_LOGGER_ERROR_);
      return false;
    }
    aurostd::StringstreamClean(xvasp.INCAR_orig); xvasp.INCAR_orig << xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR"); // DID REREAD
    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    return true;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// POSCAR
namespace KBIN {
  bool VASP_Produce_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::VASP_Produce_POSCAR():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_POSCAR: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    aurostd::StringstreamClean(xvasp.POSCAR);
    aurostd::StringstreamClean(xvasp.POSCAR_orig);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);
    //
    // IMPLICIT or EXPLICIT or EXTERNAL for POSCAR
    Krun=(Krun && (vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") ||
          vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT") ||
          vflags.KBIN_VASP_POSCAR_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_POSCAR_MODE_IMPLICIT] or [VASP_POSCAR_MODE_EXPLICIT] or [VASP_POSCAR_MODE_EXPLICIT] must be specified " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }

    // IMPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT")) {  // [VASP_POSCAR_MODE_IMPLICIT] construction
      aus << "00000  MESSAGE POSCAR  generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(!vflags.KBIN_VASP_POSCAR_FILE.flag("PROTOTYPE")) {
        aus << "EEEEE  [VASP_POSCAR_FILE] In POSCAR_MODE_IMPLICIT you must specify PROTOTYPE " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      } else {
        std::vector<string> tokens,tokens2,atomABC;
        std::string structure,label,parameters="";  // FIX NRL
        vector<double> volumeABC;
        structure=aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]PROTOTYPE=",1,TRUE);
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
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(GetAtomVolume(atomABC[isp]));} //CO20181129
          //[CO20181106]for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomABC[isp])));}
          xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);
          // xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomABC[0])),atomABC[1],GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomABC[1])),-1.0); // OLD WAY
        }
        if(atomABC.size()==2 && volumeABC.size()==1) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(0.0);}
          // xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],0.0,atomABC[1],0.0,volumeABC[0]); // OLD WAY
          xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,volumeABC[0],LIBRARY_MODE_HTQC);
        }
        if(atomABC.size()==2 && volumeABC.size()==2) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(volumeABC[isp]);};
          // xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],volumeABC[0],atomABC[1],volumeABC[1],-1.0); // OLD WAY
          xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);
        }
        if(done==FALSE) {
          aus << "EEEEE  POSCAR_MODE_IMPLICIT error in the PROTOTYPE definition" << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "EEEEE  [VASP_POSCAR_FILE]PROTOTYPE=" << structure << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        // done

        xvasp.POSCAR << xvasp.str;
        // cerr << xvasp.POSCAR.str() << endl;
      }
    }
    // EXPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT")) {  // [VASP_POSCAR_MODE_EXPLICIT] construction
      if(vflags.KBIN_VASP_POSCAR_FILE.flag("KEYWORD") && !vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP") && !vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus << "00000  MESSAGE POSCAR  generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.POSCAR.str(aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]",0));
        // [OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.POSCAR,"[VASP_POSCAR_FILE]");
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.POSCAR,"[VASP_POSCAR_FILE]");
        xvasp.str=xstructure(xvasp.POSCAR,IOVASP_AUTO);  // load structure
        aurostd::StringstreamClean(xvasp.POSCAR);
        xvasp.str.iomode=IOVASP_POSCAR;
        xvasp.POSCAR << xvasp.str;
      } else if(!vflags.KBIN_VASP_POSCAR_FILE.flag("KEYWORD") && (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP") || vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT"))) {
        aus << "00000  MESSAGE POSCAR  generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // normal get ONE of ONE
        xvasp.POSCAR.str(aurostd::substring2string(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_START_,_VASP_POSCAR_MODE_EXPLICIT_STOP_,-1));
        //[SD20220520 - OBSOLETE]if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP")) {
        //[SD20220520 - OBSOLETE]  if(aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_START_) &&
        //[SD20220520 - OBSOLETE]      aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_STOP_))
        // [OBSOLETE]	  aurostd::ExtractLastToStringstreamEXPLICIT(FileAFLOWIN,xvasp.POSCAR,_VASP_POSCAR_MODE_EXPLICIT_START_,_VASP_POSCAR_MODE_EXPLICIT_STOP_);
        //[SD20220520 - OBSOLETE]    aurostd::ExtractLastToStringstreamEXPLICIT(AflowIn,xvasp.POSCAR,_VASP_POSCAR_MODE_EXPLICIT_START_,_VASP_POSCAR_MODE_EXPLICIT_STOP_);
        //[SD20220520 - OBSOLETE]}
        if(!xvasp.POSCAR.str().empty()) {xvasp.str=xstructure(xvasp.POSCAR,IOVASP_AUTO);}   // load structure
        // get ONE of MANY
        if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
          xvasp.str=vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(xvasp.POSCAR_index);
        }
        // GOT IT
        aurostd::StringstreamClean(xvasp.POSCAR);
        xvasp.str.iomode=IOVASP_POSCAR;
        xvasp.POSCAR << xvasp.str;
      } else {
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT] do not confuse aflow !!" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT] Possible modes " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] POSCAR EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_POSCAR_FILE]POSCAR of the structure example" << endl;
        aus << "[VASP_POSCAR_FILE]-98.5397" << endl;
        aus << "[VASP_POSCAR_FILE]   4.18890 0.00000 0.00000" << endl;
        aus << "[VASP_POSCAR_FILE]  -2.09445 3.62769 0.00000" << endl;
        aus << "[VASP_POSCAR_FILE]   0.00000 0.00000 5.12300" << endl;
        aus << "[VASP_POSCAR_FILE]2 4" << endl;
        aus << "[VASP_POSCAR_FILE]Direct" << endl;
        aus << "[VASP_POSCAR_FILE]0.33333 0.66666 0.25000 Au" << endl;
        aus << "[VASP_POSCAR_FILE]0.66666 0.33333 0.75000 Au" << endl;
        aus << "[VASP_POSCAR_FILE]0.00000 0.00000 0.00000 Ti" << endl;
        aus << "[VASP_POSCAR_FILE]0.00000 0.00000 0.50000 Ti" << endl;
        aus << "[VASP_POSCAR_FILE]0.33333 0.66666 0.75000 Ti" << endl;
        aus << "[VASP_POSCAR_FILE]0.66666 0.33333 0.25000 Ti" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] POSCAR EXPLICIT MODE with START/STOP" << endl;
        aus << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
        aus << _VASP_POSCAR_MODE_EXPLICIT_START_ << endl;
        aus << "POSCAR of the structure example with START/STOP" << endl;
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
        aus << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT] Note " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT]START must be present and no [VASP_POSCAR_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT]STOP  must be present and no [VASP_POSCAR_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  or [VASP_POSCAR_FILE] present and NO START/STOP" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("EXTERNAL")) {  // [VASP_POSCAR_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE POSCAR  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_POSCAR_FILE.flag("COMMAND") && vflags.KBIN_VASP_POSCAR_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_POSCAR_MODE]FILE=  and  [VASP_POSCAR_MODE]COMMAND=  can not be used together " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_POSCAR_FILE.flag("COMMAND") && (vflags.KBIN_VASP_POSCAR_FILE.flag("FILE") || !vflags.KBIN_VASP_POSCAR_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_POSCAR_FILE.flag("FILE")) {
          file=aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]FILE=",1,TRUE);
          aus << "00000  MESSAGE POSCAR  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_POSCAR;
          aus << "00000  MESSAGE POSCAR  generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_POSCAR << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR POSCAR file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR POSCAR file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POSCAR << aurostd::file2string(file);
        xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);  // load structure
      }
      if(vflags.KBIN_VASP_POSCAR_FILE.flag("COMMAND") && !vflags.KBIN_VASP_POSCAR_FILE.flag("FILE")) {
        file=aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]COMMAND=",1,FALSE);
        aus << "00000  MESSAGE POSCAR  generation from command= '" << file << "'" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_POSCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_POSCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR POSCAR file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR POSCAR file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POSCAR << aurostd::file2string(file);       // load POSCAR
        xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);              // load structure
        aurostd::RemoveFile("./_aflow_POSCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp  //CO20200502 - threadID
      }
    }
    // POSCAR DONE **************************************************
    xvasp.POSCAR_orig << xvasp.POSCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);

    // must modify POSCAR before calculating everything else
    if(!kflags.KBIN_POCC){  //CO20180420 - VERY SLOW AND UNNECESSARY - let's do these mods when we read-in/run in the ARUNS
      if(Krun) Krun=(Krun && KBIN::VASP_Modify_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,vflags));
    }

    // CHECK for negative determinant
    if(det(xvasp.str.scale*xvasp.str.lattice)<0.0) {
      aus << "EEEEE  POSCAR ERROR: the triple product of the basis vectors is negative                      " << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  POSCAR ERROR: exchange two basis vectors and adjust the atomic positions accordingly   " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }

    // some useful LDEBUG
    if(Krun && LDEBUG) {
      // STRUCTURE IS GENERATED
      FileMESSAGE <<  endl;
      FileMESSAGE <<  "******** STRUCTURE IN CARTESIAN *****************************" << endl;
      xvasp.str.SetCoordinates(_COORDS_CARTESIAN_);
      FileMESSAGE <<  xvasp.str << endl;
      FileMESSAGE <<  "******** STRUCTURE IN FRACTIONAL ****************************" << endl;
      xvasp.str.SetCoordinates(_COORDS_FRACTIONAL_);
      FileMESSAGE <<  xvasp.str << endl;
      FileMESSAGE <<  "*************************************************************" << endl;
      // xvasp.str.write_klattice_flag=TRUE;
      FileMESSAGE <<  "SCALE" << endl;
      FileMESSAGE <<  xvasp.str.scale << endl;
      FileMESSAGE <<  "DIRECT LATTICE (with scale)" << endl;
      FileMESSAGE <<  xvasp.str.scale*(xvasp.str.lattice) << endl;
      FileMESSAGE <<  "RECIPROCAL LATTICE" << endl;
      FileMESSAGE <<  (xvasp.str.klattice) << endl;
      FileMESSAGE <<  "ORTHOGONALITY (a*b')/2pi=I" << endl;
      FileMESSAGE <<  (xvasp.str.scale*(xvasp.str.lattice))*trasp(xvasp.str.klattice)/(2.0*pi) << endl;
      FileMESSAGE <<  "*************************************************************" << endl;
    }
    // done produced and modified
    return Krun;
  };  // KBIN::VASP_Produce_POSCAR
}

namespace KBIN {
  bool VASP_Produce_POSCAR(_xvasp& xvasp) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool Krun=TRUE;
    aurostd::StringstreamClean(xvasp.POSCAR);
    aurostd::StringstreamClean(xvasp.POSCAR_orig);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);
    xvasp.POSCAR << xvasp.str;
    // POSCAR done
    xvasp.POSCAR_orig << xvasp.POSCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
    return Krun;
  };  // KBIN::VASP_Produce_POSCAR
}

namespace KBIN {
  bool VASP_Modify_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    ostringstream aus;
    bool Krun=TRUE;
    if(xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated")==FALSE) {
      aus << "EEEEE  KBIN::VASP_Modify_POSCAR: can`t modify POSCAR if it does not exist ! " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // return Krun;
    // xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
    //LDEBUG=TRUE;

    // POSCAR must be modified before doing the KPOINTS
    //ME20181121 - Outsourced to be used by other functions
    Krun = VASP_Convert_Unit_Cell(xvasp, vflags, aflags, FileMESSAGE, aus);

    if(LDEBUG) cerr << "vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("EQUAL_EQUAL")) {  // [VASP_POSCAR_FILE]VOLUME=                       /*************** POSCAR **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[VASP_POSCAR_FILE]VOLUME=",FALSE);
      aus << "00000  MESSAGE POSCAR  IMPLICIT Volume = " << factor << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  IMPLITIC Old Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(factor);
      aus << "00000  MESSAGE POSCAR  IMPLITIC New Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume =]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL")) {  // [VASP_FORCE_OPTION]VOLUME=                                                    /*************** POSCAR **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[VASP_FORCE_OPTION]VOLUME=",FALSE);
      double factor=vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedutype<double>("EQUAL_EQUAL");
      aus << "00000  MESSAGE POSCAR  FORCE Volume = " << factor << " (factor1=" << factor1 << ")  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  FORCE Old Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(factor);
      aus << "00000  MESSAGE POSCAR  FORCE New Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume =]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("MULTIPLY_EQUAL")) {  // [VASP_POSCAR_FILE]VOLUME*=             /*************** POSCAR **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[VASP_POSCAR_FILE]VOLUME*=",FALSE);
      //     double factor=aurostd::string2utype<double>(vflags.KBIN_VASP_POSCAR_FILE_VOLUME.getattachedscheme("MULTIPLY_EQUAL"));
      //      cerr << "CORMAC MULTIPLY_EQUAL=" << factor << endl;
      aus << "00000  MESSAGE POSCAR  IMPLICIT Volume *= " << factor << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  IMPLITIC Old Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()*factor);
      aus << "00000  MESSAGE POSCAR  IMPLITIC New Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume *=]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL")) {  // [VASP_FORCE_OPTION]VOLUME*=                                                    /*************** POSCAR **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[VASP_FORCE_OPTION]VOLUME*=",FALSE);
      double factor=vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedutype<double>("MULTIPLY_EQUAL");
      aus << "00000  MESSAGE POSCAR  FORCE Volume *= " << factor << " (factor1=" << factor1 << ")  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  FORCE Old Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()*factor);
      aus << "00000  MESSAGE POSCAR  FORCE New Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume *=]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag(PLUS_EQUAL)=" << vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("PLUS_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("PLUS_EQUAL")) {  // [VASP_POSCAR_FILE]VOLUME+=               /*************** POSCAR **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[VASP_POSCAR_FILE]VOLUME+=",FALSE);
      aus << "00000  MESSAGE POSCAR  IMPLICIT Volume += " << factor << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE POSCAR  IMPLITIC Old Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()+factor);
      aus << "00000  MESSAGE POSCAR  IMPLITIC New Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume +=]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(PLUS_EQUAL)=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL")) {  // [VASP_FORCE_OPTION]VOLUME+=                                                    /*************** POSCAR **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[VASP_FORCE_OPTION]VOLUME+=",FALSE);
      double factor=vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedutype<double>("PLUS_EQUAL");
      aus << "00000  MESSAGE POSCAR  FORCE Volume += " << factor << " (factor1=" << factor1 << ")  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE POSCAR  FORCE Old Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()+factor);
      aus << "00000  MESSAGE POSCAR  FORCE New Volume= " << xvasp.str.Volume() << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume +=]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    // POSCAR done
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {
      if(0) {
        aus << "00000  MESSAGE-OPTION  XXXXX" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // ------------------------------------
    // end
    // xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
    return Krun;
  }; // KBIN::VASP_Produce_POSCAR
}

namespace KBIN {
  /// @brief convert POSCAR into the format that is consistent with the binary version
  ///
  /// @param xvasp xvasp object
  /// @param kflags KBIN flags
  ///
  /// @note Versions 6.3+ of VASP can only be compiled in parallel. Attemping to run VASP in serial
  /// will cause VASP, as well as AFLOW, to hang (not crash), wasting precious CPU hours.
  /// Therefore, either one would need to always call VASP using MPI commands or read the version
  /// from the binary, which would need to be updated periodically. We have chosen to go with the
  /// latter, as the worst case is an unknown VASP version.
  ///
  /// @authors
  /// @mod{SD,20220923,updated for vasp6.3}
  /// @mod{CO,20210713,added aflags}
  /// @mod{ME,20190220,created function}
  void convertPOSCARFormat(_xvasp& xvasp, const _kflags& kflags) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string vasp_path_full=kflags.KBIN_BIN;
    double vaspVersion=KBIN::getVASPVersionDouble(vasp_path_full); //SD20220331
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vasp_path_full=\"" << vasp_path_full << "\"" << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vaspVersion=" << vaspVersion << endl;}
    if (aurostd::isequal(vaspVersion,0.0)) {  //CO20210713
      stringstream message;
      message << "Could not determine VASP version " << kflags.KBIN_BIN << "." << " POSCAR may have the wrong format.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cout, _LOGGER_WARNING_);
      return;
    }
    //SD20220923 - if version is not found, default to VASP5 format
    if ((vaspVersion>=5.0 || vaspVersion==0.0) && !xvasp.str.is_vasp5_poscar_format) { //CO20210713
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " converting POSCAR to vasp5 format" << endl;} 
      xvasp.str.is_vasp4_poscar_format = false;
      xvasp.str.is_vasp5_poscar_format = true;
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.POSCAR=" << endl << xvasp.POSCAR.str() << endl;}
    }
    else if ((vaspVersion<5.0) && !xvasp.str.is_vasp4_poscar_format) { //CO20210713
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " converting POSCAR to vasp4 format" << endl;}
      xvasp.str.is_vasp4_poscar_format = true;
      xvasp.str.is_vasp5_poscar_format = false;
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.POSCAR=" << endl << xvasp.POSCAR.str() << endl;}
    } 
  }

  // CONVERT_UNIT_CELL STUFF
  bool VASP_Convert_Unit_Cell(_xvasp& xvasp, _vflags& vflags, _aflags& aflags, ofstream& FileMESSAGE, ostringstream& aus) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    bool Krun = true;

    if(LDEBUG) aus << "00000  MESSAGE-OPTION  KBIN::VASP_Convert_Unit_Cell: BEGIN" << endl;  //SC20200410
    if(LDEBUG) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);                      //SC20200410

    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry) {        /*************** POSCAR **************/
      //    aus << "00000  MESSAGE POSCAR  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_PRESERVE construction /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  PRESERVE Unit Cell " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_STANDARD_PRIMITIVE construction /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  STANDARD_PRIMITIVE Unit Cell " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.Standard_Primitive_UnitCellForm();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      //CO START
      //CO, fix issue with iatoms tag becoming atom names
      bool write_inequivalent_flag=xvasp.str.write_inequivalent_flag;
      //CO END
      string bravais_lattice_type=xvasp.str.bravais_lattice_type,bravais_lattice_variation_type=xvasp.str.bravais_lattice_variation_type,pearson_symbol=xvasp.str.pearson_symbol;
      xvasp.str.title=xvasp.str.title+" [Standard_Primitive Unit Cell Form]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      //CO START
      //CO, fix if write_inequivalent_flag is present
      xvasp.str.write_inequivalent_flag=FALSE;
      //CO, fix if write_inequivalent_flag is present
      xvasp.POSCAR << xvasp.str;
      xvasp.str.clear();xvasp.POSCAR >> xvasp.str;  //CO, this is important, clear all symmetry stuff as the whole lattice has changed //DX20191220 - uppercase to lowercase clear
      //CO add these flags to prevent recalculation and wasted effort
      xvasp.str.Standard_Lattice_calculated=TRUE;
      xvasp.str.Standard_Lattice_primitive=TRUE;
      //CO add these flags to prevent recalculation and wasted effort
      //CO, fix issue with iatoms tag becoming atom names
      xvasp.str.write_inequivalent_flag=write_inequivalent_flag;
      //CO END
      xvasp.str.bravais_lattice_type=bravais_lattice_type;xvasp.str.bravais_lattice_variation_type=bravais_lattice_variation_type;xvasp.str.pearson_symbol=pearson_symbol;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      aus << "00000  MESSAGE POSCAR  STANDARD_PRIMITIVE Unit Cell Lattice = ["+xvasp.str.bravais_lattice_type << "," << xvasp.str.bravais_lattice_variation_type << "," << xvasp.str.pearson_symbol << "]" << "  " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL construction /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  STANDARD_CONVENTIONAL Unit Cell " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.Standard_Conventional_UnitCellForm();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      //CO START
      //CO, fix issue with iatoms tag becoming atom names
      bool write_inequivalent_flag=xvasp.str.write_inequivalent_flag;
      //CO END
      string bravais_lattice_type=xvasp.str.bravais_lattice_type,bravais_lattice_variation_type=xvasp.str.bravais_lattice_variation_type,pearson_symbol=xvasp.str.pearson_symbol;
      xvasp.str.title=xvasp.str.title+" [Standard_Conventional Unit Cell Form]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      //CO START
      //CO, fix if write_inequivalent_flag is present
      xvasp.str.write_inequivalent_flag=FALSE;
      //CO, fix if write_inequivalent_flag is present
      xvasp.POSCAR << xvasp.str;
      xvasp.str.clear();xvasp.POSCAR >> xvasp.str;  //CO, this is important, clear all symmetry stuff as the whole lattice has change //DX20191220 - uppercase to lowercase clear
      //CO add these flags to prevent recalculation and wasted effort
      xvasp.str.Standard_Lattice_calculated=TRUE;
      xvasp.str.Standard_Lattice_conventional=TRUE;
      //CO add these flags to prevent recalculation and wasted effort
      //CO, fix issue with iatoms tag becoming atom names
      xvasp.str.write_inequivalent_flag=write_inequivalent_flag;
      //CO END
      xvasp.str.bravais_lattice_type=bravais_lattice_type;xvasp.str.bravais_lattice_variation_type=bravais_lattice_variation_type;xvasp.str.pearson_symbol=pearson_symbol;
      // xvasp.str.clear();xvasp.POSCAR >> xvasp.str; //DX20191220 - uppercase to lowercase clear
      // cout << xvasp.str << endl;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      aus << "00000  MESSAGE POSCAR  STANDARD_CONVENTIONAL Unit Cell Lattice = ["+xvasp.str.bravais_lattice_type << "," << xvasp.str.bravais_lattice_variation_type << "," << xvasp.str.pearson_symbol << "]" << "  " << Message(__AFLOW_FILE__,aflags) << endl; //CO
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //CO
      // cerr << det(xvasp.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_NIGGLI construction                     /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  NIGGLI Unit Cell Reduction " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.NiggliUnitCellForm();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      xvasp.str.title=xvasp.str.title+" [Niggli Unit Cell Form]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      // cerr << det(xvasp.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(MINKOWSKI)=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_MINKOWSKI construction               /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  MINKOWSKI Basis Reduction " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.MinkowskiBasisReduction();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      xvasp.str.title=xvasp.str.title+" [Minkowski Basis Reduction]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      // cerr << det(xvasp.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_INCELL construction                     /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  INCELL Unit Cell Basis " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.BringInCell();
      xvasp.str.title=xvasp.str.title+" [Bring In Cell Basis]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_COMPACT construction                   /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  COMPACT Unit Cell Basis " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.BringInCompact();
      xvasp.str.title=xvasp.str.title+" [Bring In Compact Basis]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_WIGNERSEITZ construction           /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  WIGNERSEITZ Unit Cell Basis " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.BringInWignerSeitz();
      xvasp.str.title=xvasp.str.title+" [WignerSeitz Basis]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_CARTESIAN construction               /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  CARTESIAN Basis Coordinates" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetCoordinates(_COORDS_CARTESIAN_);
      //[CO20190216 - seems mistake]xvasp.str.title=xvasp.str.title+" [WignerSeitz Basis]";
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_FRACTIONAL construction            /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  FRACTIONAL Basis Coordinate" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetCoordinates(_COORDS_FRACTIONAL_);
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_DIRECT construction                    /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  DIRECT Basis Coordinate" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetCoordinates(_COORDS_FRACTIONAL_);
      aurostd::StringstreamClean(xvasp.POSCAR);
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) aus << "00000  MESSAGE-OPTION  KBIN::VASP_Convert_Unit_Cell: END" << endl;  //SC20200410
    if(LDEBUG) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);                    //SC20200410

    return Krun;
  }
}

namespace KBIN {
  bool VASP_Reread_POSCAR(_xvasp& xvasp) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    ofstream FileMESSAGE;
    _aflags aflags;aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");
    return VASP_Reread_POSCAR(xvasp,FileMESSAGE,aflags);
  }
  bool VASP_Reread_POSCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " reading POSCAR from: " << xvasp.Directory << endl;}
    if(!aurostd::FileExist(xvasp.Directory+"/POSCAR")) {
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR not present in directory: "+xvasp.Directory,aflags,FileMESSAGE,std::cout,_LOGGER_ERROR_);
      return false;
    }
    if(aurostd::FileEmpty(xvasp.Directory+"/POSCAR")) {
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR empty in directory: "+xvasp.Directory,aflags,FileMESSAGE,std::cout,_LOGGER_ERROR_);
      return false;
    }
    aurostd::StringstreamClean(xvasp.POSCAR_orig); xvasp.POSCAR_orig << xvasp.POSCAR.str();
    aurostd::StringstreamClean(xvasp.POSCAR); xvasp.POSCAR << aurostd::file2string(xvasp.Directory+"/POSCAR"); // DID REREAD
    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.POSCAR=" << endl << xvasp.POSCAR.str() << endl;}
    if(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xvasp.POSCAR.str()).empty()){
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR file contents empty in directory: "+xvasp.Directory,aflags,FileMESSAGE,std::cout,_LOGGER_ERROR_);
      return false;
    }
    return true;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// KPOINTS
namespace KBIN {
  bool VASP_Produce_KPOINTS(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::VASP_Produce_KPOINTS():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_KPOINTS: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    bool IMPLICIT=FALSE;
    aurostd::StringstreamClean(xvasp.KPOINTS);
    xvasp.str.kpoints_k1=0;xvasp.str.kpoints_k2=0;xvasp.str.kpoints_k3=0;        // RESET KPOINTS
    xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;  // RESET SHIFTS

    xvasp.str.kpoints_kmax=0;
    xvasp.str.kpoints_kscheme.clear();
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",FALSE);

    // IMPLICIT or EXPLICIT or EXTERNAL for KPOINTS
    Krun=(Krun && (vflags.KBIN_VASP_KPOINTS_MODE.flag("IMPLICIT") ||
          vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT") ||
          vflags.KBIN_VASP_KPOINTS_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_KPOINTS_MODE_IMPLICIT] or [VASP_KPOINTS_MODE_EXPLICIT] or [VASP_KPOINTS_MODE_EXTERNAL] must be specified " << Message(__AFLOW_FILE__,aflags) << endl; //ME20181113
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated")==FALSE) {
      aus      << "EEEEE  INCAR   generation: POSCAR must be called before " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // EXPLICIT **************************************************
    // kpoints explicit through string xvasp.KPOINTS
    if(Krun && vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT")) {  // [VASP_KPOINTS_MODE_EXPLICIT] construction
      IMPLICIT=FALSE;
      if(vflags.KBIN_VASP_KPOINTS_FILE.flag("KEYWORD") && !vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE KPOINTS generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.KPOINTS << vflags.KBIN_VASP_KPOINTS_EXPLICIT.str(); //ME20181113
        // [OBSOLETE]      aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.KPOINTS,"[VASP_KPOINTS_FILE]");
        // [ME20181226 OBSOLETE]      aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.KPOINTS,"[VASP_KPOINTS_FILE]");
        KBIN::XVASP_KPOINTS_string2numbers(xvasp);
      } else if(!vflags.KBIN_VASP_KPOINTS_FILE.flag("KEYWORD") && vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")) {
        string START="[VASP_KPOINTS_MODE_EXPLICIT]START";
        string STOP="[VASP_KPOINTS_MODE_EXPLICIT]STOP";
        aus << "00000  MESSAGE KPOINTS generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.KPOINTS << vflags.KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(); //ME20181113
        // [ME20181226 OBSOLETE] if(aurostd::substring2bool(AflowIn,START) && aurostd::substring2bool(AflowIn,STOP))
        // [OBSOLETE]	aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.KPOINTS,START,STOP);
        // [ME20181226 OBSOLETE] aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.KPOINTS,START,STOP);
        KBIN::XVASP_KPOINTS_string2numbers(xvasp);
      } else {
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT] do not confuse aflow !!" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT] Possible modes " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] KPOINTS EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
        aus << "[VASP_KPOINTS_FILE]KPOINTS of the structure example with START/STOP" << endl;
        aus << "[VASP_KPOINTS_FILE]0" << endl;
        aus << "[VASP_KPOINTS_FILE]Monkhorst-Pack" << endl;
        aus << "[VASP_KPOINTS_FILE]7 7 5" << endl;
        aus << "[VASP_KPOINTS_FILE]0 0 0" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] KPOINTS EXPLICIT MODE with START/STOP" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]START" << endl;
        aus << "KPOINTS of the structure example with START/STOP" << endl;
        aus << "0" << endl;
        aus << "Monkhorst-Pack" << endl;
        aus << "7 7 5" << endl;
        aus << "0 0 0" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT] Note " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT]START must be present and no [VASP_KPOINTS_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT]STOP  must be present and no [VASP_KPOINTS_FILE]" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "EEEEE  or [VASP_KPOINTS_FILE] present and NO START/STOP" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // IMPLICIT **************************************************
    // kpoints implicit through string xvasp.KPOINTS
    if(Krun && vflags.KBIN_VASP_KPOINTS_MODE.flag("IMPLICIT")) // IT MIGHT NOT CONTAIN OPTION SO IT IS ALL DEFAULT && vflags.KBIN_VASP_KPOINTS_FILE)  // [VASP_KPOINTS_MODE_IMPLICIT] construction
    { //CO20200106 - patching for auto-indenting
      IMPLICIT=TRUE;
      aus << "00000  MESSAGE KPOINTS generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      string stringKPPRA;
      int NK=1;
      xvasp.str.kpoints_k1=1;xvasp.str.kpoints_k2=1;xvasp.str.kpoints_k3=1;
      xvasp.str.kpoints_mode=0;
      xvasp.str.kpoints_kscheme.clear();
      bool LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry=vflags.KBIN_VASP_KPOINTS_KMODE.isentry;  // DEFAULT
      bool LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry=vflags.KBIN_VASP_KPOINTS_KPPRA.isentry;  // DEFAULT
      bool LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry;  // DEFAULT
      bool LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry=vflags.KBIN_VASP_KPOINTS_KSHIFT.isentry;  // DEFAULT
      int LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int=vflags.KBIN_VASP_KPOINTS_KMODE.content_int;  // DEFAULT
      int LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int=vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;  // DEFAULT
      string LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;  // DEFAULT
      string LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string=vflags.KBIN_VASP_KPOINTS_KSHIFT.content_string;  // DEFAULT

      string STRING_KPOINTS_TO_SHOW="KPOINTS";

      if(vflags.KBIN_VASP_RUN.flag("STATIC")) { // do the switching
        STRING_KPOINTS_TO_SHOW="KPOINTS_STATIC";
        LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry;  // overrides kscheme
        LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;  // overrides kscheme
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int=vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.content_int;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KMODE overrides KPOINTS_KMODE generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KMODE: LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KMODE: LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.content_int << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int=vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_uint;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KPPRA overrides KPOINTS_KPPRA generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KPPRA: LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KPPRA: LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_uint << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME overrides KPOINTS_KSCHEME generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME: LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME: LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string=vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.content_string;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSHIFT overrides KPOINTS_KSHIFT generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSHIFT: LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSHIFT: LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.content_string << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }

      if(LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int==0 || LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry==FALSE) {
        // KSCHEME ******************************
        if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry==FALSE) {
          xvasp.str.kpoints_kscheme="Monkhorst-Pack";  // DEFAULT FIX
          xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;  // DEFAULT FIX
          aus << "00000  MESSAGE-DEFAULT " << STRING_KPOINTS_TO_SHOW << "_KSCHEME=" << xvasp.str.kpoints_kscheme << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          xvasp.str.kpoints_kscheme="Monkhorst-Pack";  // DEFAULT FIX
          xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;  // DEFAULT FIX
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='M' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='m') xvasp.str.kpoints_kscheme="Monkhorst-Pack";
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='G' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='g') xvasp.str.kpoints_kscheme="Gamma";
          //	if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='A' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='a') xvasp.str.kpoints_kscheme="Auto";
          //	if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='L' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='l') xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='A' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='a') xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='C' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='c') xvasp.str.kpoints_kscheme="Cartesian";
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='K' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='k') xvasp.str.kpoints_kscheme="Cartesian";
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME=" << xvasp.str.kpoints_kscheme << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(xvasp.str.kpoints_kscheme==DEFAULT_KSCHEME) {
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " Calculating structure lattice" << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          xvasp.str.GetLatticeType();
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " Found Lattice=" << xvasp.str.bravais_lattice_variation_type << " - " << Message(__AFLOW_FILE__,aflags) << endl;
          xvasp.str.kpoints_kscheme="Monkhorst-Pack";
          //          if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC") xvasp.str.kpoints_kscheme="Gamma";  // also add RHL
          if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC" || xvasp.str.bravais_lattice_variation_type=="RHL") xvasp.str.kpoints_kscheme="Gamma";  //CO+DX20200404 - 3-fold symmetries REQUIRE Gamma-centered
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME=\"" << xvasp.str.kpoints_kscheme << "\" " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // 
        // KPPRA ******************************
        if(LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry==FALSE) {
          NK=1; // DEFAULT FIX
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KPPRA=NNNN is missing, taking NNNN=" << NK << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          NK=LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KPPRA=" << NK << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        stringKPPRA=KPPRA(xvasp.str,NK);
        // VERBOSE on LOCK /SCREEN
        if(TRUE) FileMESSAGE << stringKPPRA;
        if(TRUE) {
          aus.precision(5);
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KPPRA routine ["
            << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << "]=" << xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3 << "=["
            << modulus(xvasp.str.klattice(1)/((double) xvasp.str.kpoints_k1)) << ","
            << modulus(xvasp.str.klattice(2)/((double) xvasp.str.kpoints_k2)) << ","
            << modulus(xvasp.str.klattice(3)/((double) xvasp.str.kpoints_k3)) << "]" << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // KSHIFT ******************************
        if(LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry==FALSE) {
          xvasp.str.kpoints_s1=xvasp.str.kpoints_s2=xvasp.str.kpoints_s3=0.0;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KSHIFT= X X X  is missing, taking X X X = 0.0 0.0 0.0 " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          xvasp.str.kpoints_s1=xvasp.str.kpoints_s2=xvasp.str.kpoints_s3=0.0;
          stringstream oaus;
          oaus.str(LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string);
          oaus >> xvasp.str.kpoints_s1 >> xvasp.str.kpoints_s2 >> xvasp.str.kpoints_s3;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KSHIFT= " << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << " " << " " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // create KPOINTS ******************************
        aurostd::StringstreamClean(xvasp.KPOINTS);
        aurostd::StringstreamClean(xvasp.KPOINTS_orig);
        xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA - aflow with Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK << "]" << endl;
        xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;  // MODE AUTO 0
        xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
        xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
        xvasp.KPOINTS.precision(3);
        xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;

        aurostd::StringstreamClean(aus);
        aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK  << "]    with " << xvasp.str.kpoints_kscheme << "   " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
      } else {
        aurostd::StringstreamClean(aus);
        aus << "EEEEE             Only [VASP_KPOINTS_FILE]KMODE=0 is supported ! " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_KPOINTS_MODE.flag("EXTERNAL")) {  // [VASP_KPOINTS_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE KPOINTS  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_KPOINTS_FILE.flag("COMMAND") && vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_KPOINTS_MODE]FILE=  and  [VASP_KPOINTS_MODE]COMMAND=  can not be used together " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_KPOINTS_FILE.flag("COMMAND") && (vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE") || !vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE")) {
          // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]FILE=",1,TRUE);
          file = vflags.KBIN_VASP_KPOINTS_FILE.getattachedscheme("FILE"); //ME20181113
          aus << "00000  MESSAGE KPOINTS  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_KPOINTS;
          aus << "00000  MESSAGE KPOINTS  generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_KPOINTS << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR KPOINTS file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR KPOINTS file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.KPOINTS << aurostd::file2string(file);
      }
      if(vflags.KBIN_VASP_KPOINTS_FILE.flag("COMMAND") && !vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE")) {
        // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=",1,FALSE);
        file = vflags.KBIN_VASP_KPOINTS_FILE.getattachedscheme("COMMAND"); //ME20181113
        aus << "00000  MESSAGE KPOINTS  generation from command= '" << file << "'" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_KPOINTS."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp  //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_KPOINTS."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name  //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR KPOINTS file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR KPOINTS file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.KPOINTS << aurostd::file2string(file);       // load KPOINTS
        aurostd::RemoveFile("./_aflow_KPOINTS."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp //CO20200502 - threadID
      }
      KBIN::XVASP_KPOINTS_string2numbers(xvasp);                          // create KPOINTS numbers
    }
    // KPOINTS DONE **************************************************
    if(IMPLICIT==FALSE) KBIN::XVASP_KPOINTS_string2numbers(xvasp);
    xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
    xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
    // KPOINTS done
    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);
    // cerr << xvasp.str.kpoints_kscheme << endl;
    return Krun;
  };  // KBIN::VASP_Produce_KPOINTS
}

namespace KBIN {
  bool VASP_Modify_KPOINTS(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;

    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KEEPK")==FALSE && 0) {             /*************** KPOINTS **************/
      vflags.KBIN_VASP_WRITE_KPOINTS=TRUE;
      if(_isodd(xvasp.str.kpoints_k1)) xvasp.str.kpoints_k1++;
      if(_isodd(xvasp.str.kpoints_k2)) xvasp.str.kpoints_k2++;
      if(_isodd(xvasp.str.kpoints_k3)) xvasp.str.kpoints_k3++;
      xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
      xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
      KBIN::XVASP_KPOINTS_numbers2string(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      aus << "00000  MESSAGE KPOINTS Option Tune - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.isentry) {                        /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.content_string << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KEEPK")) {                        /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KEEPK - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  KBIN::XVASP_KPOINTS_OPERATION(xvasp,""); //    KBIN::XVASP_KPOINTS_KEEPK(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("EVEN")) {                        /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=EVEN - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xeven,Yeven,Zeven"); //    KBIN::XVASP_KPOINTS_EVEN(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("ODD")) {                         /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=ODD - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_EVEN")) {           /*************** KPOINTS **************/ //CO20210315 - the 'GAMMA' here is not appropriate unless 'Gamma' is set below, shifting the even by 0.5 puts it on the same 'grid' as Gamma, see https://www.vasp.at/wiki/index.php/KPOINTS#Monkhorst-Pack, not changing for legacy
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSHIFT_GAMMA_EVEN - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xevenshift,Yevenshift,Zevenshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp); //CO20210315 - 'Gamma' should be set here, see https://www.vasp.at/wiki/index.php/KPOINTS#Monkhorst-Pack, not changing for legacy
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_ODD")) {            /*************** KPOINTS **************/ //CO20210315 - the 'GAMMA' here is not appropriate unless 'Gamma' is set below, shifting the odd by 0.5 offsets it from the Gamma 'grid', see https://www.vasp.at/wiki/index.php/KPOINTS#Monkhorst-Pack, not changing for legacy
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSHIFT_GAMMA_ODD - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xoddshift,Yoddshift,Zoddshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_ODD(xvasp); //CO20210315 - 'Gamma' should be set here, see https://www.vasp.at/wiki/index.php/KPOINTS#Monkhorst-Pack, not changing for legacy
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_MONKHORST_PACK")) {      /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_MONKHORST_PACK - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Monkhorst-Pack");
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_GAMMA")) {               /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_GAMMA - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Gamma");
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && (vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_AUTO"))) {              /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_AUTO  Calculating structure lattice" << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.GetLatticeType();
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_AUTO  Found Lattice=" << xvasp.str.bravais_lattice_variation_type << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      xvasp.str.kpoints_kscheme="Monkhorst-Pack";
      //     if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC") xvasp.str.kpoints_kscheme="Gamma"; // add RHL
      if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC" || xvasp.str.bravais_lattice_variation_type=="RHL") xvasp.str.kpoints_kscheme="Gamma";  //CO+DX20200404 - 3-fold symmetries REQUIRE Gamma-centered
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_" << xvasp.str.kpoints_kscheme << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,xvasp.str.kpoints_kscheme);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("GAMMA")) {
      vflags.KBIN_VASP_WRITE_KPOINTS=TRUE;
      xvasp.str.kpoints_k1=1;xvasp.str.kpoints_k2=1;xvasp.str.kpoints_k3=1;
      xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
      xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
      KBIN::XVASP_KPOINTS_numbers2string(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=GAMMA - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT")) {
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=IBZKPT - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // ------------------------------------
    // end
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);
    return Krun;
  }; // KBIN::VASP_Modify_KPOINTS
}

namespace KBIN {
  bool VASP_Reread_KPOINTS(_xvasp& xvasp) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    ofstream FileMESSAGE;
    _aflags aflags;aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");
    return VASP_Reread_KPOINTS(xvasp,FileMESSAGE,aflags);
  }
  bool VASP_Reread_KPOINTS(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    if(!aurostd::FileExist(xvasp.Directory+"/KPOINTS")) {
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"KPOINTS not present in directory: "+xvasp.Directory,aflags,FileMESSAGE,std::cout,_LOGGER_ERROR_);
      return false;
    }
    aurostd::StringstreamClean(xvasp.KPOINTS_orig); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    aurostd::StringstreamClean(xvasp.KPOINTS); xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS"); // DID REREAD
    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    return true;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// POTCAR
namespace KBIN {
  bool VASP_Find_DATA_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar,string &AUIDPotcar) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string pseudopotential="nothing";
    FilePotcar="";DataPotcar="";
    if(aurostd::substring2bool(species_pp,"pot_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_LDA;
    if(aurostd::substring2bool(species_pp,"pot_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_GGA;
    if(aurostd::substring2bool(species_pp,"pot_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA;
    if(aurostd::substring2bool(species_pp,"potpaw_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN;

    if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: species_pp=" << species_pp << endl; 
    if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: pseudopotential=" << pseudopotential << endl;

    bool found=FALSE;

    string string2load_list="AFLOW_PSEUDOPOTENTIALS_LIST_TXT:"+species_pp+"/POTCAR";
    string string2load_data="AFLOW_PSEUDOPOTENTIALS_TXT:"+species_pp+"/POTCAR";
    if(aurostd::substring2bool(init::InitLoadString(string2load_list,FALSE),"POTCAR")) {
      found=TRUE;
      FilePotcar=init::InitLoadString(string2load_list,FALSE);
      DataPotcar=init::InitLoadString(string2load_data,FALSE);
      string file_tmp=aurostd::TmpFileCreate("DataPotcar");
      AUIDPotcar=aurostd::file2auid(file_tmp);
      if(aurostd::FileExist(file_tmp)) aurostd::RemoveFile(file_tmp);

      if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: species_pp=" << species_pp << endl; 
      if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: pseudopotential=" << pseudopotential << endl;
      if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: AUIDPotcar=" << AUIDPotcar << endl; 
    }
    return found;
  }
}

namespace KBIN {
  bool VASP_Find_FILE_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar,string &AUIDPotcar) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string pseudopotential="nothing";
    FilePotcar="";DataPotcar="";
    if(aurostd::substring2bool(species_pp,"pot_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_LDA;
    if(aurostd::substring2bool(species_pp,"pot_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_GGA;
    if(aurostd::substring2bool(species_pp,"pot_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA;
    if(aurostd::substring2bool(species_pp,"potpaw_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN;

    if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: species_pp=" << species_pp << endl; 
    if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: pseudopotential=" << pseudopotential << endl;

    bool found=FALSE;
    FilePotcar="";
    for(uint j=0;j<vVASP_POTCAR_DIRECTORIES.size()&&!found;j++) {   // cycle through possible directories
      for(uint k=0;k<=9&&!found;k++) {                     // cycle through fixing current before/after
        string FilePotcark;
        if(k==0) FilePotcark=species_pp;      // k==0 dont touch
        if(k>0) FilePotcark=vVASP_POTCAR_DIRECTORIES.at(j)+"/"+species_pp;
        //if(k==1) FilePotcark=FilePotcark; // nothing done //CO20190329 - compiler doesn't like FilePotcark=FilePotcark
        if(k==2) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/"+DEFAULT_VASP_POTCAR_DATE);  // current POST
        if(k==3) aurostd::StringSubst(FilePotcark,pseudopotential,DEFAULT_VASP_POTCAR_DATE+"/"+pseudopotential);  // current PRE
        if(k==4) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/dateP1");   // dateP1 POST
        if(k==5) aurostd::StringSubst(FilePotcark,pseudopotential,"dateP1/"+pseudopotential);   // dateP1 PRE
        if(k==6) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/dateP2");   // dateP2 POST
        if(k==7) aurostd::StringSubst(FilePotcark,pseudopotential,"dateP2/"+pseudopotential);   // dateP2 PRE
        if(k==8) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/dateP3");   // dateP3 POST
        if(k==9) aurostd::StringSubst(FilePotcark,pseudopotential,"dateP3/"+pseudopotential);   // dateP3 PRE
        for(uint l=0;l<=1&&!found;l++) {                  // cycle through POTCAR as postfix
          FilePotcar=FilePotcark;                      // default
          if(l==0) FilePotcar=FilePotcark;             // current NO
          if(l==1) FilePotcar=FilePotcark+"/POTCAR";   // POTCAR after
          // clean up to avoid comments //
          // TESTING POTCARS
          if(!aurostd::substring2bool(FilePotcar,"POTCAR")) FilePotcar=FilePotcar+"/POTCAR"; // some fix, _AFLOWIN_ might have some personality issues
          // DEBUG=TRUE;
          FilePotcar=aurostd::CleanFileName(FilePotcar);
          if(aurostd::FileExist(FilePotcar) && !aurostd::FileEmpty(FilePotcar)) found=TRUE;
          if(LDEBUG) cout << "DDDDD  POTCAR  (" << species_pp << ",j=" << j << ",k=" << k << ",l=" << l << ")  [FilePotcar=" <<  FilePotcar << "] found=" << found << " aurostd::FileExist(FilePotcar)=" << aurostd::FileExist(FilePotcar) << endl;
        } // l-cycle through POTCAR as postfix
      } // k-cycle through fixing current before/after
    } // j-cycle through possible directories

    if(found) {
      stringstream aus;
      ifstream FileINPUT;
      FileINPUT.clear();
      FileINPUT.open(FilePotcar.c_str(),std::ios::in);
      char c;
      while (FileINPUT.get(c)) aus.put(c);
      FileINPUT.clear();FileINPUT.close();
      DataPotcar=aus.str();

      AUIDPotcar=aurostd::file2auid(FilePotcar);
      if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: species_pp=" << species_pp << endl; 
      if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: pseudopotential=" << pseudopotential << endl;
      if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: AUIDPotcar=" << AUIDPotcar << endl; 
    }

    return found;
  }
}

namespace KBIN {
  bool VASP_Produce_POTCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    if(AflowIn.length()==0) {throw aurostd::xerror(__AFLOW_FILE__,"KBIN::VASP_Produce_POTCAR():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_POTCAR: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    string::size_type sub_size1,sub_size2;
    string subS,subS1,subS2,subSDIR="$POTCARDIR";
    ostringstream aus;
    bool Krun=TRUE;
    aurostd::StringstreamClean(xvasp.POTCAR);
    xvasp.POTCAR_AUID.clear();
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_changed",FALSE);
    aurostd::StringstreamClean(xvasp.POTCAR_POTENTIALS);

    // IMPLICIT or EXPLICIT or EXTERNAL for POTCAR
    Krun=(Krun && (vflags.KBIN_VASP_POTCAR_MODE.flag("IMPLICIT") || vflags.KBIN_VASP_POTCAR_MODE.flag("EXPLICIT") || vflags.KBIN_VASP_POTCAR_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_POTCAR_MODE_IMPLICIT] or [VASP_POTCAR_MODE_EXPLICIT] or [VASP_POTCAR_MODE_EXPLICIT] must be specified " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // IMPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POTCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POTCAR_FILE.flag("KEYWORD")) { // [VASP_POTCAR_MODE_IMPLICIT] construction
      // Prepare POTCAR
      aus << "00000  MESSAGE POTCAR  generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) {         /*************** POTCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=\"" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << "\" - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // Prepare POTCAR
      ifstream FileINPUT;
      vector<string> FilePotcars;
      string FilePotcar,DataPotcar,AUIDPotcar;
      if(!vflags.KBIN_VASP_POTCAR_FILE.flag("SYSTEM_AUTO")) {
        subS2="[VASP_POTCAR_FILE]"; // STRING TO SEARCH
        subS1=AflowIn.substr(AflowIn.find(subS2));
        while (aurostd::substring2bool(subS1,subS2)) {
          subS=subS1;
          sub_size1=subS.find(subS2);
          sub_size2=(subS.substr(sub_size1)).find("\n");
          subS=subS.substr(sub_size1+subS2.length(),sub_size2-subS2.length());
          subS1=subS1.substr(sub_size1+1);
          if(aurostd::substring2bool(subS,subSDIR)) {
            string subSpre,subSpost;
            sub_size1=subS.find(subSDIR)+subSDIR.length();
            sub_size2=(subS.substr(sub_size1)).find("\n");
            subSpost=subS.substr(sub_size1,sub_size2);
            sub_size1=0;
            sub_size2=(subS.substr(sub_size1)).find(subSDIR);
            subSpre=subS.substr(sub_size1,sub_size2);
            subS.clear();
            subS=subSpre+vVASP_POTCAR_DIRECTORIES.at(0)+subSpost;
          }
          FilePotcars.push_back(subS);
        }
      }
      if(vflags.KBIN_VASP_POTCAR_FILE.flag("SYSTEM_AUTO")) {
        for(uint i=0;i<xvasp.str.species.size();i++) {
          if(xvasp.str.species.at(i)!="") {
            subS="";
            if(!vflags.KBIN_VASP_POTCAR_FILE.flag("PREFIX")) {
              subS+=vVASP_POTCAR_DIRECTORIES.at(0);
              subS+=DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE+"/"+DEFAULT_VASP_POTCAR_DATE+"/";
            }
            if(vflags.KBIN_VASP_POTCAR_FILE.flag("PREFIX")) {
              subS=vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("PREFIX");  //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",1,TRUE);  //CO20181226
              if(LDEBUG) cerr << "DEBUG subS=" << subS << endl;
              if(aurostd::substring2bool(subS,subSDIR)) {
                subS=vVASP_POTCAR_DIRECTORIES.at(0)+aurostd::substring2string(subS,subSDIR,1,TRUE);
              }
              if(LDEBUG) cerr << "DEBUG subS=" << subS << endl;
            }
            subS+="/"+xvasp.str.species.at(i);
            if(!vflags.KBIN_VASP_POTCAR_FILE.flag("SUFFIX"))
              subS+=DEFAULT_VASP_POTCAR_SUFFIX;
            if(vflags.KBIN_VASP_POTCAR_FILE.flag("SUFFIX"))
              subS+="/"+vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("SUFFIX"); //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",1,TRUE);  //CO20181226
            FilePotcars.push_back(subS);
          }
        }
      }
      // FIX PSEUDOPOTENTIALS
      //CO20181226 - originally this was applied regardless of AUTO_PSEUDOPOTENTIALS
      //this causes problems, as AUTO_PSEUDOPOTENTIALS strips PP information ([VASP_POTCAR_FILE]) in automatic aflow generation (see avasp)
      //therefore, Mn_pv becomes Mn (WRONG)
      //I believe the patch is to do this ONLY if !AUTO_PSEUDOPOTENTIALS
      //if this guard is removed, please patch fixEmptyAtomNames() to be force_fix=false (APL/aflow_apl_kphonons)
      //ME20190211 - This patch causes species_pp to be unoccupied if AUTO_PSEUDOPOTENTIALS (fixed below)
      //ME20190313 - This breaks when [VASP_POTCAR_FILE] is a path - need to break the path apart
      if(!vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) {  //CO20181226
        vector<string> tokens;  //ME20190313
        for(uint i=0;i<FilePotcars.size() && i<xvasp.str.species_pp.size();i++) {
          aurostd::string2tokens(FilePotcars.at(i), tokens, "/");  //ME20190313
          xvasp.str.species_pp.at(i) = tokens.back();  //ME20190313
          //xvasp.str.species_pp.at(i)=FilePotcars.at(i);
        }
      }
      // IF AUTO_PSEUDOPOTENTIALS
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry << endl;}
      if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) {
        for(uint i=0;i<FilePotcars.size();i++) {
          string cleanname=KBIN::VASP_PseudoPotential_CleanName(FilePotcars.at(i));
          aurostd::StringSubst(cleanname,"/POTCAR","");
          vector<string> tokens;aurostd::string2tokens(cleanname,tokens,"/");
          if(tokens.size()>0) cleanname=tokens.at(tokens.size()-1);
          if(LDEBUG) cerr << "[VASP_Produce_POTCAR]: vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_PBE_KIN")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN+"/"+AVASP_Get_PseudoPotential_PAW_PBE_KIN(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_LDA_KIN")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN+"/"+AVASP_Get_PseudoPotential_PAW_LDA_KIN(cleanname);

          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_PBE")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+"/"+AVASP_Get_PseudoPotential_PAW_PBE(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_GGA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA+"/"+AVASP_Get_PseudoPotential_PAW_GGA(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_LDA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA+"/"+AVASP_Get_PseudoPotential_PAW_LDA(cleanname);

          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="pot_PBE")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POT_PBE+"/"+AVASP_Get_PseudoPotential_PBE(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="pot_GGA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POT_GGA+"/"+AVASP_Get_PseudoPotential_GGA(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="pot_LDA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POT_LDA+"/"+AVASP_Get_PseudoPotential_LDA(cleanname);

          //ME20190211 - Occupy species_pp
          aurostd::string2tokens(FilePotcars[i], tokens, "/");
          if (tokens.size() > 0) xvasp.str.species_pp[i] = tokens.back();

          if(LDEBUG) cerr << "[VASP_Produce_POTCAR]: cleanname=" << cleanname << endl;
          if(LDEBUG) cerr << "[VASP_Produce_POTCAR]: FilePotcars.at(" << i << ")=" << FilePotcars.at(i) << endl;
        }
      }

      // PRINT POTCARS
      for(uint i=0;i<FilePotcars.size();i++)
        aus << "00000  MESSAGE POTCAR  Loading potcar file = [" << FilePotcars.at(i) << "]" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // LOAD POTCARS (test some fixes because of different machines)

      string FilePotcarNaked;
      vector<string> tokens;

      bool AFLOW_PSEUDOPOTENTIALS=FALSE;
      if(aurostd::substring2bool(init::InitGlobalObject("AFLOW_PSEUDOPOTENTIALS"),"1")) {
        AFLOW_PSEUDOPOTENTIALS=TRUE;
        aus << "00000  MESSAGE POTCAR  AFLOW_PSEUDOPOTENTIALS=TRUE - " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      //	AFLOW_PSEUDOPOTENTIALS=FALSE;

      for(uint i=0;i<FilePotcars.size();i++) {                        // cycle through all POTCARS
        aurostd::StringSubst(FilePotcars.at(i),"//","/");
        // strip everything
        tokens.clear();
        FilePotcarNaked=FilePotcars.at(i);
        for(uint iext=0;iext<XHOST.vext.size();iext++) { 
          FilePotcarNaked=aurostd::RemoveSubStringFirst(FilePotcarNaked,"/POTCAR"+XHOST.vext.at(iext));
        }
        aurostd::string2tokens(FilePotcarNaked,tokens,"/");
        xvasp.POTCAR_POTENTIALS << tokens.at(tokens.size()-1);
        // now build it

        bool found_DATA=FALSE;
        bool found_FILE=FALSE;
        if(AFLOW_PSEUDOPOTENTIALS) {
          found_DATA=KBIN::VASP_Find_DATA_POTCAR(FilePotcars.at(i),FilePotcar,DataPotcar,AUIDPotcar);
          if(found_DATA) {
            aus << "00000  MESSAGE POTCAR  DATA: Found potcar FilePotcar=" << FilePotcar << " - DataPotcar.size()=" << DataPotcar.size() << " - AUIDPotcar=" << AUIDPotcar << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.POTCAR << DataPotcar << endl;
            xvasp.POTCAR_AUID.push_back(AUIDPotcar);
            // [OBSOLETE] aus << "00000  MESSAGE POTCAR  DATA: Found potcar xvasp.POTCAR.str().size()=" << xvasp.POTCAR.str().size() << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	      
          }
        } else {
          found_FILE=KBIN::VASP_Find_FILE_POTCAR(FilePotcars.at(i),FilePotcar,DataPotcar,AUIDPotcar);
          if(found_FILE) {
            // FileINPUT.clear();FileINPUT.open(FilePotcar.c_str(),std::ios::in);
            // char c;while (FileINPUT.get(c)) xvasp.POTCAR.put(c);
            // FileINPUT.clear();FileINPUT.close();
            aus << "00000  MESSAGE POTCAR  FILE: Found potcar FilePotcar=" << FilePotcar << " - DataPotcar.size()=" << DataPotcar.size() << " - AUIDPotcar=" << AUIDPotcar << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.POTCAR << DataPotcar;// << endl; file has already new line
            xvasp.POTCAR_AUID.push_back(AUIDPotcar);
            // [OBSOLETE] aus << "00000  MESSAGE POTCAR  FILE: Found potcar xvasp.POTCAR.str().size()=" << xvasp.POTCAR.str().size() << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	      
          } else {
            aus << "EEEEE  POTCAR [" << FilePotcars.at(i).c_str() << "] not found! (aflow_ivasp.cpp) " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
            Krun=FALSE;                                         // DO NOT RUN ANYMORE
            return Krun;
          }
        }
      }
    }
    // EXPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POTCAR_MODE.flag("EXPLICIT")) {   // [VASP_POTCAR_MODE_EXPLICIT] construction
      // Prepare POTCAR
      aus << "00000  MESSAGE POTCAR  generation EXPLICIT " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      {
        ifstream FileINPUT;
        string FileNameINPUT=xvasp.Directory+"/"+_AFLOWIN_;
        string line;bool lflag=FALSE;
        FileINPUT.open(FileNameINPUT.c_str(),std::ios::in);
        while (getline(FileINPUT,line)) {
          if(lflag==TRUE) xvasp.POTCAR << line << endl;
          //	if(lflag==TRUE)  xvasp.POTCAR+=line+"\n";
          if(aurostd::substring2bool(line,"[VASP_POTCAR_MODE_EXPLICIT]")) lflag=TRUE;
        }
        FileINPUT.clear();FileINPUT.close();
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_POTCAR_MODE.flag("EXTERNAL")) {  // [VASP_POTCAR_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE POTCAR  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_POTCAR_FILE.flag("COMMAND") && vflags.KBIN_VASP_POTCAR_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_POTCAR_MODE]FILE=  and  [VASP_POTCAR_MODE]COMMAND=  can not be used together " << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_POTCAR_FILE.flag("COMMAND") && (vflags.KBIN_VASP_POTCAR_FILE.flag("FILE") || !vflags.KBIN_VASP_POTCAR_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_POTCAR_FILE.flag("FILE")) {
          file=vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("FILE"); //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]FILE=",1,TRUE); //CO20181226
          aus << "00000  MESSAGE POTCAR  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_POTCAR;
          aus << "00000  MESSAGE POTCAR  generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_POTCAR << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR POTCAR file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR POTCAR file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POTCAR << aurostd::file2string(file);
      }
      if(vflags.KBIN_VASP_POTCAR_FILE.flag("COMMAND") && !vflags.KBIN_VASP_POTCAR_FILE.flag("FILE")) {
        file=vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("COMMAND"); //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",1,FALSE);  //CO20181226
        aus << "00000  MESSAGE POTCAR  generation from command= '" << file << "'" << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_POTCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_POTCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR POTCAR file=" << file << " does not exist! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR POTCAR file=" << file << " is empty! " << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POTCAR << aurostd::file2string(file);       // load POTCAR
        aurostd::RemoveFile("./_aflow_POTCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp  //CO20200502 - threadID
      }
    }

    // get xvasp.POTCAR_ENMAX **************************************************
    xPOTCAR potcar;potcar.GetProperties(xvasp.POTCAR);
    xvasp.POTCAR_ENMAX=potcar.ENMAX;
    xvasp.POTCAR_ENMIN=potcar.ENMIN;
    aus << "00000  MESSAGE POTCAR  ENMAX max   = " << xvasp.POTCAR_ENMAX << " - " << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "00000  MESSAGE POTCAR  ENMIN min   = " << xvasp.POTCAR_ENMIN << " - " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // get xvasp.POTCAR_PAW **************************************************
    xvasp.POTCAR_PAW=potcar.POTCAR_PAW;
    xvasp.POTCAR_TYPE=potcar.POTCAR_TYPE;
    aus << "00000  MESSAGE POTCAR  PAW = " << (xvasp.POTCAR_PAW ? "TRUE":"FALSE") << " - " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aus << "00000  MESSAGE POTCAR  POTCAR_TYPE=" <<  xvasp.POTCAR_TYPE << " - " << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    xvasp.POTCAR_orig << xvasp.POTCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_generated",TRUE);

    // done now return
    return Krun;
  };  // KBIN::VASP_Produce_POTCAR
}



namespace KBIN {
  bool VASP_Modify_POTCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    // ------------------------------------
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {
      if(0) {
        aus << "00000  MESSAGE-OPTION  XXXXX" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // end
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_generated",TRUE);
    return Krun;
  }; // KBIN::VASP_Produce_POTCAR
}


namespace KBIN {
  bool VASP_Reread_POTCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    if(!aurostd::FileExist(xvasp.Directory+"/POTCAR")) {
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"POTCAR not present in directory: "+xvasp.Directory,aflags,FileMESSAGE,std::cout,_LOGGER_ERROR_);
      return false;
    }
    aurostd::StringstreamClean(xvasp.POTCAR_orig); xvasp.POTCAR_orig << xvasp.POTCAR.str();
    aurostd::StringstreamClean(xvasp.POTCAR); xvasp.POTCAR << aurostd::file2string(xvasp.Directory+"/POTCAR"); // DID REREAD
    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_POTCAR_changed",TRUE);
    return true;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INCAR MODIFICATIONS

// ***************************************************************************
// KBIN::XVASP_Get_NPAR_NCORE
// ***************************************************************************
namespace KBIN {
  void XVASP_Get_NPAR_NCORE(const _xvasp& xvasp,const _aflags& aflags,int& NPAR,int& NCORE){  //CO20210315
    NPAR=0,NCORE=1; //reset

    if(xvasp.NCPUS >   0 && xvasp.NCPUS <=  32) {NPAR=2;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS >  32 && xvasp.NCPUS <=  64) {NPAR=4;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS >  64 && xvasp.NCPUS <= 128) {NPAR=8;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS > 128 && xvasp.NCPUS <= 256) {NPAR=16;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS > 256 && xvasp.NCPUS <= 512) {NPAR=32;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS > 512 && xvasp.NCPUS <=1024) {NPAR=64;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS >1024) {NPAR=128;NCORE=xvasp.NCPUS/NPAR;}
    if(NPAR==0) {NPAR=4;NCORE=1;}

    //DX COME BACK: should these be set by machine? looks like they are being overridden below
    if(xvasp.NCPUS==32)  {NPAR=4;NCORE=32;} // test for DX
    if(xvasp.NCPUS==48)  {NPAR=4;NCORE=48;} // test for DX
    if(xvasp.NCPUS==44)  {NPAR=4;NCORE=44;} // test for DX

    // marylou fulton super computer center
    if(xvasp.NCPUS==4)  { // best 4 times
      //  NCPUS04 NPAR01 NCORE04 18576.588 20.6407
      //  NCPUS04 NPAR01 NCORE03 18289.541 20.3217
      //  NCPUS04 NPAR01 NCORE01 18223.492 20.2483
      //  NCPUS04 NPAR01 NCORE02 18176.841 20.1965
      NPAR=1;NCORE=4;  // trying to make NPAR*NCORE=NCPUS 
    }
    if(xvasp.NCPUS==6)  { // best 4 times
      //  NCPUS06 NPAR03 NCORE02 12758.774 21.2646
      //  NCPUS06 NPAR03 NCORE04 12724.680 21.2078
      //  NCPUS06 NPAR03 NCORE05 12678.754 21.1313
      //  NCPUS06 NPAR03 NCORE01 12669.103 21.1152
      NPAR=3;NCORE=2;  // trying to make NPAR*NCORE=NCPUS 
    }
    if(xvasp.NCPUS==8)  { // best 4 times   
      //  NCPUS08 NPAR02 NCORE06 10757.264 23.905
      //  NCPUS08 NPAR02 NCORE04 10731.053 23.8468
      //  NCPUS08 NPAR02 NCORE07 10728.917 23.842
      //  NCPUS08 NPAR02 NCORE01 10707.747 23.795
      NPAR=2;NCORE=4;  // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==12)  { // best times
      //  NCPUS12 NPAR03 NCORE07 7256.892 24.1896
      //  NCPUS12 NPAR03 NCORE11 7253.081 24.1769
      //  NCPUS12 NPAR03 NCORE04 7251.042 24.1701
      //  NCPUS12 NPAR03 NCORE09 7248.894 24.163
      //  NCPUS12 NPAR03 NCORE05 7245.983 24.1533
      //  NCPUS12 NPAR03 NCORE06 7242.303 24.141
      //  NCPUS12 NPAR03 NCORE02 7239.272 24.1309
      //  NCPUS12 NPAR03 NCORE12 7239.034 24.1301
      //  NCPUS12 NPAR03 NCORE03 7234.666 24.1156
      //  NCPUS12 NPAR03 NCORE10 7230.139 24.1005
      //  NCPUS12 NPAR03 NCORE08 7175.204 23.9173
      //  NCPUS12 NPAR03 NCORE01 7169.198 23.8973
      NPAR=3;NCORE=2; // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==14)  { // best times
      // NCPUS14 NPAR02 NCORE13 6647.652 25.852
      // NCPUS14 NPAR02 NCORE04 6625.906 25.7674
      // NCPUS14 NPAR02 NCORE01 6622.861 25.7556
      // NCPUS14 NPAR02 NCORE12 6606.536 25.6921
      // NCPUS14 NPAR02 NCORE14 6605.233 25.687
      // NCPUS14 NPAR02 NCORE06 6602.645 25.677
      // NCPUS14 NPAR02 NCORE02 6583.709 25.6033
      // NCPUS14 NPAR02 NCORE08 6581.525 25.5948
      // NCPUS14 NPAR02 NCORE05 6577.889 25.5807
      // NCPUS14 NPAR02 NCORE03 6573.488 25.5636
      // NCPUS14 NPAR02 NCORE07 6562.956 25.5226
      // NCPUS14 NPAR02 NCORE10 6536.529 25.4198
      // NCPUS14 NPAR02 NCORE11 6518.928 25.3514
      NPAR=2;NCORE=7; // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==16)  { // best times
      //  NCPUS16 NPAR04 NCORE04 6423.776 28.5501
      //  NCPUS16 NPAR04 NCORE02 6418.438 28.5264
      //  NCPUS16 NPAR04 NCORE13 6416.928 28.5197
      //  NCPUS16 NPAR04 NCORE05 6416.292 28.5169
      //  NCPUS16 NPAR04 NCORE14 6415.565 28.5136
      //  NCPUS16 NPAR04 NCORE10 6407.344 28.4771
      //  NCPUS16 NPAR04 NCORE01 6405.943 28.4709
      //  NCPUS16 NPAR04 NCORE16 6403.422 28.4597
      //  NCPUS16 NPAR04 NCORE06 6323.368 28.1039
      //  NCPUS16 NPAR04 NCORE08 6323.029 28.1024
      NPAR=4;NCORE=8; // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==24)  { // best times
      //  NCPUS24 NPAR03 NCORE02 5293.432 35.2895
      //  NCPUS24 NPAR03 NCORE05 5286.059 35.2404
      //  NCPUS24 NPAR03 NCORE12 5277.060 35.1804
      //  NCPUS24 NPAR03 NCORE09 5274.808 35.1654
      //  NCPUS24 NPAR03 NCORE17 5267.672 35.1178
      //  NCPUS24 NPAR03 NCORE15 5265.254 35.1017
      //  NCPUS24 NPAR03 NCORE01 5262.930 35.0862
      //  NCPUS24 NPAR03 NCORE23 5254.012 35.0267
      //  NCPUS24 NPAR03 NCORE06 5247.212 34.9814
      //  NCPUS24 NPAR03 NCORE24 5219.655 34.7977
      //  NCPUS24 NPAR03 NCORE07 5185.569 34.5705
      NPAR=3;NCORE=12;  // trying to make NPAR*NCORE=NCPUS if possible
    }

    if(xvasp.NCPUS==28)  { // best times  (delay %)
      //  NCPUS28 NPAR04 NCORE12 4817.701 37.471  1.97962
      //  NCPUS28 NPAR04 NCORE07 4809.654 37.4084 1.80929
      //  NCPUS28 NPAR04 NCORE13 4799.799 37.3318 1.60068
      //  NCPUS28 NPAR04 NCORE05 4799.391 37.3286 1.59204
      //  NCPUS28 NPAR04 NCORE20 4797.978 37.3176 1.56213
      //  NCPUS28 NPAR04 NCORE22 4794.220 37.2884 1.48259
      //  NCPUS28 NPAR04 NCORE01 4785.839 37.2232 1.30518
      //  NCPUS28 NPAR04 NCORE18 4769.633 37.0971 0.962135
      //  NCPUS28 NPAR04 NCORE23 4767.426 37.08   0.915418
      //  NCPUS28 NPAR04 NCORE06 4724.180 36.7436
      NPAR=4;NCORE=6;  // trying to make NPAR*NCORE=NCPUS if possible
    }

    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS") && xvasp.NCPUS==32)  {
      NPAR=4;NCORE=32;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO") && xvasp.NCPUS==32)  {
      NPAR=4;NCORE=32;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA") && xvasp.NCPUS==32)  {
      NPAR=4;NCORE=32;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA") && xvasp.NCPUS==40)  {
      NPAR=4;NCORE=40;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA") && xvasp.NCPUS==32)  {
      NPAR=4;NCORE=32;
    }

    //  cerr << "xvasp.NCPUS=" << xvasp.NCPUS << endl;;
    // cerr << "NPAR=" << NPAR << endl;
  }
}

// ***************************************************************************
// KBIN::XVASP_MPI_Autotune
// ***************************************************************************
namespace KBIN {
  void XVASP_MPI_Autotune(_xvasp& xvasp,_aflags &aflags,bool VERBOSE) {
    string FileContent="",strline="";
    int imax=0;

    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"IALGO",TRUE)  || aurostd::substring2bool(strline,"#IALGO",TRUE) ||
          aurostd::substring2bool(strline,"LPLANE",TRUE) || aurostd::substring2bool(strline,"#LPLANE",TRUE) ||
          aurostd::substring2bool(strline,"NPAR",TRUE)   || aurostd::substring2bool(strline,"#NPAR",TRUE)   ||
          aurostd::substring2bool(strline,"LSCALU",TRUE) || aurostd::substring2bool(strline,"#LSCALU",TRUE) ||
          aurostd::substring2bool(strline,"NSIM",TRUE)   || aurostd::substring2bool(strline,"#NSIM",TRUE)) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::VASP_MPI_Autotune)" << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing MPI-AUTOTUNE [AFLOW] begin" << endl;
    if(VERBOSE) xvasp.INCAR << "# MPI tuning for Linux cluster with CPUS=" << xvasp.NCPUS << endl;
    if(VERBOSE) xvasp.INCAR << "# look at the VASP manual" << endl;
    //[CO20210315 - OBSOLETE link]if(VERBOSE) xvasp.INCAR << "# http://www.teragridforum.org/mediawiki/images/7/7e/Performance_on_Ranger_of_the_medium_size_problem.pdf" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("IALGO=48",_incarpad_) << " # MPI algorithm" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LPLANE=.TRUE.",_incarpad_) << " # parallelization" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(xvasp.NCPUS),_incarpad_) << " # number of CPUs" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(sqrt(xvasp.NCPUS)),_incarpad_) << " # number of CPUs" << endl;
    int NPAR=0,NCORE=1;

    KBIN::XVASP_Get_NPAR_NCORE(xvasp,aflags,NPAR,NCORE); //CO20210315

    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_EOS and xvasp.NCPUS==32",_incarpad_) << endl; 
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_DRACO and xvasp.NCPUS==32",_incarpad_) << endl; 
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_COBRA and xvasp.NCPUS==32",_incarpad_) << endl; 
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA") && xvasp.NCPUS==40)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_COBRA and xvasp.NCPUS==40",_incarpad_) << endl; 
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_HYDRA and xvasp.NCPUS==32",_incarpad_) << endl; 
    }

    xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(NPAR),_incarpad_) << " # lookup table from NCPUS=" << xvasp.NCPUS << endl;    // 6 times faster than NPAR=cpus.
    xvasp.INCAR << aurostd::PaddedPOST("NCORE="+aurostd::utype2string(NCORE),_incarpad_) << " # lookup table from NCPUS=" << xvasp.NCPUS << endl;    
    //  xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(NPAR),_incarpad_) << " # seems to give good results with all MPI" << endl;    // 6 times faster than NPAR=cpus.
    //  xvasp.INCAR << aurostd::PaddedPOST("NPAR=4",_incarpad_) << " # seems to give good results with all MPI" << endl;    // 6 times faster than NPAR=cpus..
    xvasp.INCAR << aurostd::PaddedPOST("NSIM=4",_incarpad_) << " # tune size vector/matrix product" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LSCALU=.FALSE.",_incarpad_) << " # neglect scalapack" << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing MPI-AUTOTUNE [AFLOW] end" << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_System_Auto
// ***************************************************************************
namespace KBIN {
  void XVASP_INCAR_System_Auto(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"SYSTEM",TRUE)    || aurostd::substring2bool(strline,"#SYSTEM",TRUE) ||
          aurostd::substring2bool(strline,"INFO",TRUE)      || aurostd::substring2bool(strline,"#INFO",TRUE) ||
          aurostd::substring2bool(strline,"PROTOTYPE",TRUE) || aurostd::substring2bool(strline,"#PROTOTYPE",TRUE)) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_System_Auto)" << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.INCAR << "#AFLOW INCAR automatically generated" << endl;
    xvasp.INCAR << "SYSTEM=" << xvasp.str.title << endl;
    xvasp.INCAR << "#PROTOTYPE=" << xvasp.str.prototype << endl;
    xvasp.INCAR << "#INFO=" << xvasp.str.info << endl;
    xvasp.INCAR << FileContent;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Relax_ON
// ***************************************************************************
namespace KBIN {
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,_vflags& vflags,int number) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    int imax,isif=3;
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")==FALSE) {  // whatever is the number
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==FALSE) isif=2; // 0,1,2 FIX
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==TRUE)  isif=3;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==FALSE) isif=4;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==FALSE) isif=5;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==TRUE)  isif=6;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==TRUE)  isif=7;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL")) isif=3;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_SHAPE")) isif=4;//AS20201208
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")) {  // whatever is the number
      if(aurostd::_isodd(number)) isif=7;  // VOLUME
      if(aurostd::_iseven(number)) isif=2; // IONS
    }

    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);

    //CO20210315 - nelm is a special case
    //it was not previously written to INCAR for relaxations, allowing WSETYAWAN to write it in [VASP_INCAR_MODE_EXPLICIT] section of aflow.in for ICSD runs (aflow_avasp.cpp), see xvasp.AVASP_EXTRA_INCAR
    //however, now we write it in as it is an aflowrc parameter, as well as a VASP_FORCE_OPTION
    //only for this relaxation function (applies for AUTOTUNE), check if NELM is specified in the INCAR already
    //if so, do not override
    //the VASP_FORCE_OPTION will overwrite later, so don't worry about here
    //this hack is required to be backwards-compatible with all ICSD aflow.in's in database
    //this is the ONLY exception for not overwriting a property, use the VASP_FORCE_OPTIONs otherwise
    //the only properties that should go in the [VASP_INCAR_MODE_EXPLICIT] section are those NOT overwritten by AUTOTUNE
    bool nelm_patch=(aurostd::kvpair2bool(FileContent,"NELM","=")==false);
    //[CO20210315 - always patch, you might have NELM in INCAR section of aflow.in](vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.content_int!=AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM); //CO20200624 - default for VASP is 60, don't add the line if unnecessary

    // RELAX_MODE=ENERGY mode
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="ENERGY") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
            aurostd::substring2bool(strline,"NSW",TRUE)    || aurostd::substring2bool(strline,"#NSW",TRUE)    ||
            aurostd::substring2bool(strline,"LORBIT",TRUE)    || aurostd::substring2bool(strline,"#LORBIT",TRUE)    || //CO20180130 get spinD
            aurostd::substring2bool(strline,"ISIF",TRUE)   || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
            ( nelm_patch && (aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
            FALSE) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_ON)" << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX=ENERGY [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("IBRION=2",_incarpad_) << " # relax with Conjugate Gradient" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NSW=51",_incarpad_) << " # relax for long" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ISIF="+aurostd::utype2string(isif),_incarpad_) << " # relax appropriately" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << " # get spin decomposition" << endl; //CO20180130 get spinD
      if(nelm_patch){ //avoid writing NELM twice
        xvasp.INCAR << aurostd::PaddedPOST("NELM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.content_int),_incarpad_) << " # set max SC electronic steps" << endl;  //CO20200624
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX=ENERGY [AFLOW] end" << endl;
    }
    // RELAX=MODE=FORCES mode
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
            aurostd::substring2bool(strline,"NELMIN",TRUE)    || aurostd::substring2bool(strline,"#NELMIN",TRUE)    ||
            aurostd::substring2bool(strline,"ADDGRID",TRUE)    || aurostd::substring2bool(strline,"#ADDGRID",TRUE)    ||
            aurostd::substring2bool(strline,"EDIFFG",TRUE)    || aurostd::substring2bool(strline,"#EDIFFG",TRUE)    ||
            aurostd::substring2bool(strline,"LORBIT",TRUE)    || aurostd::substring2bool(strline,"#LORBIT",TRUE)    || //CO20180130 get spinD
            aurostd::substring2bool(strline,"NSW",TRUE)    || aurostd::substring2bool(strline,"#NSW",TRUE)    ||
            aurostd::substring2bool(strline,"ISIF",TRUE)   || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
            ( nelm_patch && (aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
            FALSE) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_ON)" << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX_MODE=FORCES [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELMIN=4",_incarpad_) << " # The forces have to be well converged" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ADDGRID=.TRUE.",_incarpad_) << " # To support finer forces calculation" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFFG="+aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG,12),_incarpad_) << " # The final structure has to have zero forces!" << endl;
      //    xvasp.INCAR << aurostd::PaddedPOST("EDIFFG=-1E-5",_incarpad_) << " # The final structure has to have zero forces!" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("IBRION=1",_incarpad_) << " # More stable algorithm" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NSW=100",_incarpad_) << " # relax for very long" << endl;
      // xvasp.INCAR << "ISIF=3         # relax everything" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ISIF="+aurostd::utype2string(isif),_incarpad_)     << " # relax appropriately" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << " # get spin decomposition" << endl; //CO20180130 get spinD
      if(nelm_patch){ //avoid writing NELM twice
        xvasp.INCAR << aurostd::PaddedPOST("NELM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.content_int),_incarpad_) << " # set max SC electronic steps" << endl;  //CO20200624
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX_MODE=FORCES [AFLOW] end" << endl;
    }
    // done now write if necessary
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME") && number>1) {  // whatever is the number
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",FALSE);
      aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
    }
  }
}

// ***************************************************************************
namespace KBIN {
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    _vflags vflags;
    vflags.KBIN_VASP_INCAR_VERBOSE=VERBOSE;
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("ALL");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_SHAPE");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_VOLUME");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME",FALSE);
    KBIN::XVASP_INCAR_Relax_ON(xvasp,vflags,1);
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Static_ON
namespace KBIN {
  void XVASP_INCAR_Static_ON(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i); // ibrion = -1,5,6,7,8 dont need to be removed
      if(aurostd::substring2bool(strline,"IBRION=0",TRUE) || aurostd::substring2bool(strline,"#IBRION=0",TRUE) || // std md
          aurostd::substring2bool(strline,"IBRION=1",TRUE) || aurostd::substring2bool(strline,"#IBRION=1",TRUE) || // quasi-newton algo to relax ions
          aurostd::substring2bool(strline,"IBRION=2",TRUE) || aurostd::substring2bool(strline,"#IBRION=2",TRUE) || // conj-grad algo to relax
          aurostd::substring2bool(strline,"IBRION=3",TRUE) || aurostd::substring2bool(strline,"#IBRION=3",TRUE) || // damping-factor algo to relax
          aurostd::substring2bool(strline,"NSW",TRUE)    || aurostd::substring2bool(strline,"#NSW",TRUE)    ||
          aurostd::substring2bool(strline,"LORBIT",TRUE)    || aurostd::substring2bool(strline,"#LORBIT",TRUE)    || //CO20180130 get spinD
          (aurostd::substring2bool(strline,"LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) || (aurostd::substring2bool(strline,"#LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry)||
          (aurostd::substring2bool(strline,"LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) || (aurostd::substring2bool(strline,"#LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry)||
          aurostd::substring2bool(strline,"IALGO",TRUE)    || aurostd::substring2bool(strline,"#IALGO",TRUE) || //SC from ICSD
          (aurostd::substring2bool(strline,"ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) ||  //SC from ICSD
          (aurostd::substring2bool(strline,"#ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) || //SC from ICSD
          aurostd::substring2bool(strline,"ISIF",TRUE)   || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
          //[CO20200624 - OBSOLETE]aurostd::substring2bool(strline,"NELM",TRUE)   || aurostd::substring2bool(strline,"#NELM",TRUE) ||  //CO20200624
          ((aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
          FALSE) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Static_ON)" << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing STATIC [AFLOW]" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.TRUE.",_incarpad_) << " # Performing STATIC  (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && !vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.FALSE.",_incarpad_) << " # Performing STATIC  (Bader OFF)" << endl;
    // if(vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("ADDGRID=.TRUE.",_incarpad_) << " # Performing STATIC  (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.TRUE.",_incarpad_) << " # Performing STATIC  (Elf ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && !vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.FALSE.",_incarpad_) << " # Performing STATIC  (Elf OFF)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("ALGO=Normal",_incarpad_) << " # Performing STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << " # Performing STATIC" << endl;  //CO20180130 get spinD
    xvasp.INCAR << aurostd::PaddedPOST("NELM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.content_int),_incarpad_) << " # Performing STATIC" << endl;  //CO20200624
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "#removing IBRION, NSW, ISIF" << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_Relax_Static_ON
namespace KBIN {
  void XVASP_INCAR_Relax_Static_ON(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(
          ((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || //CO20171003 - don't confuse PREC and SYMPREC
          (aurostd::substring2bool(strline,"#PREC",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || // Aleksey
          aurostd::substring2bool(strline,"EMIN",TRUE) || aurostd::substring2bool(strline,"#EMIN",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"EMAX",TRUE) || aurostd::substring2bool(strline,"#EMAX",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"IALGO",TRUE)    || aurostd::substring2bool(strline,"#IALGO",TRUE) || //SC from ICSD
          (aurostd::substring2bool(strline,"ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) ||  //SC from ICSD
          (aurostd::substring2bool(strline,"#ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) || //SC from ICSD
          aurostd::substring2bool(strline,"ISMEAR",TRUE) || aurostd::substring2bool(strline,"#ISMEAR",TRUE) ||
          aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE) ||
          aurostd::substring2bool(strline,"ISIF",TRUE) || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
          aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
          aurostd::substring2bool(strline,"NSW",TRUE) || aurostd::substring2bool(strline,"#NSW",TRUE) ||
          //[CO20200624 - OBSOLETE]aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) ||
          ((aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
          aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) ||
          aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT",TRUE) ||
          aurostd::substring2bool(strline,"RWIGS",TRUE) || aurostd::substring2bool(strline,"#RWIGS",TRUE) ||
          aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE) ||
          (aurostd::substring2bool(strline,"LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) || (aurostd::substring2bool(strline,"#LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry)||
          aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE) ||
          (aurostd::substring2bool(strline,"LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) || (aurostd::substring2bool(strline,"#LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry)||
          aurostd::substring2bool(strline,"LPARD",TRUE) || aurostd::substring2bool(strline,"#LPARD",TRUE) ||
          aurostd::substring2bool(strline,"IBAND",TRUE) || aurostd::substring2bool(strline,"#IBAND",TRUE) ||
          aurostd::substring2bool(strline,"LSEPB",TRUE) || aurostd::substring2bool(strline,"#LSEPB",TRUE)) {
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_Static_ON)" << endl;
          } else {
            if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
          }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX_STATIC [AFLOW]" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=0",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISMEAR="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.content_int),_incarpad_) << " # Performing RELAX_STATIC (tetrahedron method with Blochl corrections : wahyu/stefano 2009/10/09)" << endl; //CO20200624
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("",_incarpad_) << " # with -5 forces are not correct in metals, so it is good" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("",_incarpad_) << " # for semiconductors/insulator or for NO relaxations" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("SIGMA="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # Performing RELAX_STATIC (so the DOS will not spill too much from band edges)" << endl; //CO20200624
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("ALGO=Normal",_incarpad_) << " # Performing RELAX_STATIC (stefano from ICSD)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISIF=2",_incarpad_) << " # Performing RELAX_STATIC (relax ions, but single step so no relax)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("IBRION=2",_incarpad_) << " # Performing RELAX_STATIC (relax conj. grad, but single step so no relax)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NSW=0",_incarpad_) << " # Performing RELAX_STATIC (zero ionic steps, so no relax, just static)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NELM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.content_int),_incarpad_) << " # Performing RELAX_STATIC" << endl;  //CO20200624
    xvasp.INCAR << aurostd::PaddedPOST("NELMIN=2",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LCHARG=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC  (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && !vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC  (Bader OFF)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC  (Elf ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && !vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC  (Elf OFF)" << endl;
    // if(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC) {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC (forced TRUE by DIELECTRIC)" << endl;
    // } else {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    // }
    xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("#LPARD=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("#IBAND=8 9",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("#LSEPB=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("PREC=Accurate",_incarpad_) << " # Performing RELAX_STATIC (aleksey)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("EMIN= -14.0",_incarpad_) << " # Performing RELAX_STATIC (aleksey)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("EMAX=  12.0",_incarpad_) << " # Performing RELAX_STATIC (aleksey)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("EMIN= -35.0",_incarpad_) << " # Performing RELAX_STATIC (aleksey)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("EMIN= -30.0",_incarpad_) << " # Performing RELAX_STATIC (aleksey) force search for EMIN" << endl; // was 35
    xvasp.INCAR << aurostd::PaddedPOST("EMAX=  45.0",_incarpad_) << " # Performing RELAX_STATIC (aleksey)" << endl;  // was 35
    xvasp.INCAR << aurostd::PaddedPOST("NEDOS= 5000",_incarpad_) << " # Performing RELAX_STATIC (aleksey)" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "#adjusting ISMEAR,SIGMA,ISIF,IBRION,NSW,NELM,NELMIN,LORBIT,LCHARG,LAECHG,LWAVE,LELF,LPARD,IBAND,LSEPB" << endl;

    // check for LDA/GGA
    if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA") {
      xvasp.INCAR << "#Fixing RWIGS/LORBIG for LDA/GGA" << endl;
      ofstream oss;KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,oss,ON);}
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_Bands_ON
namespace KBIN {
  void XVASP_INCAR_Bands_ON(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(
          ((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || //CO20171003 - don't confuse PREC and SYMPREC
          (aurostd::substring2bool(strline,"#PREC",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || // Aleksey
          aurostd::substring2bool(strline,"ISMEAR",TRUE) || aurostd::substring2bool(strline,"#ISMEAR",TRUE) ||
          aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE) ||
          aurostd::substring2bool(strline,"IALGO",TRUE)    || aurostd::substring2bool(strline,"#IALGO",TRUE) || //SC from ICSD
          (aurostd::substring2bool(strline,"ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) ||  //SC from ICSD
          (aurostd::substring2bool(strline,"#ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) || //SC from ICSD
          aurostd::substring2bool(strline,"ISIF",TRUE) || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
          aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
          aurostd::substring2bool(strline,"NSW",TRUE) || aurostd::substring2bool(strline,"#NSW",TRUE) ||
          //[CO20200624 - OBSOLETE]aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) ||
          ((aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
          aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) ||
          aurostd::substring2bool(strline,"EMIN",TRUE) || aurostd::substring2bool(strline,"#EMIN",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"EMAX",TRUE) || aurostd::substring2bool(strline,"#EMAX",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT (Wahyu REMOVED, CORRECT)",TRUE) ||
          aurostd::substring2bool(strline,"RWIGS",TRUE) || aurostd::substring2bool(strline,"#RWIGS",TRUE) ||
          aurostd::substring2bool(strline,"LELF",TRUE) || aurostd::substring2bool(strline,"#LELF (Wahyu REMOVED, CORRECT)",TRUE) ||
          aurostd::substring2bool(strline,"ICHARG",TRUE) || aurostd::substring2bool(strline,"#ICHARG",TRUE) ||
          aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE) ||
          (aurostd::substring2bool(strline,"LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) || (aurostd::substring2bool(strline,"#LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry)||
          (aurostd::substring2bool(strline,"LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) || (aurostd::substring2bool(strline,"#LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry)||
          aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE)) {
            if(vflags.KBIN_VASP_INCAR_VERBOSE){xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Bands_ON)" << endl;}
          } else {
            if(vflags.KBIN_VASP_INCAR_VERBOSE||strline.length()){xvasp.INCAR << strline << endl;} //CO20210315
          }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX_STATIC_BANDS [AFLOW]" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("PREC=Accurate",_incarpad_) << " # Performing RELAX_STATIC_BANDS (aleksey)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISMEAR="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.content_int),_incarpad_) << " # Performing RELAX_STATIC_BANDS (put back Gauss, use 1 if metals)" << endl; //CO20200624
    // xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=-5",_incarpad_) << " # Performing RELAX_STATIC (tetrahedron method with Blochl corrections : WSETYAWAN+SC20091009)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("SIGMA="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # Performing RELAX_STATIC_BANDS (so the DOS will not spill too much from band edges)" << endl; //CO20200624
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("ALGO=Normal",_incarpad_) << " # Performing RELAX_STATIC_BANDS (stefano from ICSD)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISIF=2",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("IBRION=2",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NSW=0 ",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NELM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.content_int),_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;  //CO20200624
    xvasp.INCAR << aurostd::PaddedPOST("NELMIN=2",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;  //CO20180130 get spinD
    xvasp.INCAR << aurostd::PaddedPOST("ICHARG=11",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LCHARG=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && !vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS (Bader OFF)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS (Elf ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && !vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS (Elf OFF)" << endl;

    // if(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC) {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS (forced TRUE by DIELECTRIC)" << endl;
    // } else {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    // }
    xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << " # Performing RELAX_STATIC_BANDS" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "#adjusting ISMEAR,SIGMA,ISIF,IBRION,NSW,NELM,NELMIN,ICHARG,LCHARG,LAECHG,LWAVE" << endl;
    // cerr << "FIX KBIN::XVASP_INCAR_Relax_Static_Bands_Bands_ON(_xvasp& xvasp,bool vflags.KBIN_VASP_INCAR_VERBOSE)" << endl;

    // check for LDA/GGA
    // check for LDA/GGA
    if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA") {
      xvasp.INCAR << "#Fixing RWIGS/LORBIG for LDA/GGA" << endl;
      ofstream oss;KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,oss,ON);
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_RWIGS_Static
namespace KBIN {
  void XVASP_INCAR_RWIGS_Static(_xvasp& xvasp,_vflags& vflags,ofstream &FileMESSAGE,bool OPERATION) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // cerr << "DEBUG RWIGS" << endl;
    ostringstream aus;
    string FileContent,strline;
    int imax;
    if(OPERATION==ON)  aus << "00000  MESSAGE FORCE RWIGS_STATIC  ON : [VASP_FORCE_OPTION]RWIGS_STATIC " << " - " << Message(__AFLOW_FILE__) << endl;
    if(OPERATION==OFF) aus << "00000  MESSAGE FORCE RWIGS_STATIC  OFF: [VASP_FORCE_OPTION]RWIGS_STATIC " << " - " << Message(__AFLOW_FILE__) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"RWIGS",TRUE) || aurostd::substring2bool(strline,"#RWIGS",TRUE) ||
          aurostd::substring2bool(strline,"LORBIT",TRUE)   || aurostd::substring2bool(strline,"#LORBIT",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==ON ) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_RWIGS_Static_ON)" << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_RWIGS_Static_OFF)" << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==ON ) xvasp.INCAR << "# Performing RWIGS_STATIC=ON [AFLOW]" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==OFF) xvasp.INCAR << "# Performing RWIGS_STATIC=OFF [AFLOW]" << endl;	
    if(OPERATION==ON) {
      xvasp.INCAR << "RWIGS=";
      vector<string> tokens,tokens2;
      aurostd::string2tokens(xvasp.POTCAR.str(),tokens,"\n");
      for(uint i=0;i<tokens.size();i++)
        if(aurostd::substring2bool(tokens.at(i),"RWIGS")) {
          aurostd::string2tokens(tokens.at(i),tokens2);
          xvasp.INCAR << tokens2.at(5) << " ";
        }
      xvasp.INCAR << "# Performing Density of states" << endl;
      if(xvasp.POTCAR_TYPE!="LDA" && xvasp.POTCAR_TYPE!="GGA")
        xvasp.INCAR << aurostd::PaddedPOST("LORBIT=0",_incarpad_) << " # Performing Density of states (PAW)" << endl;
      if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA")
        xvasp.INCAR << aurostd::PaddedPOST("LORBIT=1",_incarpad_) << " # Performing Density of states (LDA/GGA)" << endl;
    }
    if(OPERATION==OFF) {;} // dummy nothng to do
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==ON ) xvasp.INCAR << "# adjusting RWIGS/LORBIT" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==OFF) xvasp.INCAR << "# removing RWIGS/LORBIT" << endl;

    if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {;} // dummy
  }
}


// double round(double x,double eps) {
// double out=x;
// out=out/eps;
// out=round(out)*eps;
// return out;
// }

// ***************************************************************************
// KBIN::XVASP_INCAR_Precision
namespace KBIN {
  void XVASP_INCAR_Precision(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    int imax;
    //  cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="LOW" || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="MEDIUM" || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="NORMAL") {
        if((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) || aurostd::substring2bool(strline,"#PREC",TRUE) || //CO20171003 - don't confuse PREC and SYMPREC
            aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="LOW") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=LOW)" << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="MEDIUM") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=MEDIUM)" << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="NORMAL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=NORMAL)" << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=HIGH)" << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=ACCURATE)" << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE" || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") {
        if(((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) || aurostd::substring2bool(strline,"#PREC",TRUE) ||  //CO20171003 - don't confuse PREC and SYMPREC
              aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE) ||
              aurostd::substring2bool(strline,"LREAL",TRUE) || aurostd::substring2bool(strline,"#LREAL",TRUE) ||
              aurostd::substring2bool(strline,"EDIFF",TRUE) || aurostd::substring2bool(strline,"#EDIFF",TRUE) ||
              aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE)) &&
            (!aurostd::substring2bool(strline,"IALGO",TRUE)  && !aurostd::substring2bool(strline,"#IALGO",TRUE)) &&
            (!aurostd::substring2bool(strline,"EDIFFG",TRUE)  && !aurostd::substring2bool(strline,"#EDIFFG",TRUE))) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=HIGH)" << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=ACCURATE)" << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      //BEGIN JJPR
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="PHONONS") {
        if((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) || aurostd::substring2bool(strline,"#PREC",TRUE) ||  //CO20171003 - don't confuse PREC and SYMPREC
            aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE) ||
            aurostd::substring2bool(strline,"LREAL",TRUE) || aurostd::substring2bool(strline,"#LREAL",TRUE) ||
            aurostd::substring2bool(strline,"EDIFF",TRUE) || aurostd::substring2bool(strline,"#EDIFF",TRUE) ||
            aurostd::substring2bool(strline,"EDIFFG",TRUE) || aurostd::substring2bool(strline,"#EDIFFG",TRUE) ||  //CO20210315
            ((aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
            aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) ||  //CO20210315
            aurostd::substring2bool(strline,"ADDGRID",TRUE) || aurostd::substring2bool(strline,"#ADDGRID",TRUE) ||  //CO20210315
            //[CO20210315 - no ALGO replaced below]  aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE)) &&
            //[CO20210315 - no ALGO replaced below](!aurostd::substring2bool(strline,"IALGO",TRUE)  && !aurostd::substring2bool(strline,"#IALGO",TRUE)) &&
            //[CO20210315 - specified below](!aurostd::substring2bool(strline,"EDIFFG",TRUE)  && !aurostd::substring2bool(strline,"#EDIFFG",TRUE)))
          FALSE) {
            if(vflags.KBIN_VASP_INCAR_VERBOSE) {
              if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="PHONONS") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=PHONONS)" << endl;
            }
          } else {
            if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
          }
      }
      //END JJPR 
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Precision PREC=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " [AFLOW] begin" << endl;

    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="LOW") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Low",_incarpad_) << " # low and fast" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_LOW,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # " << DEFAULT_VASP_PREC_ENMAX_LOW << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="MEDIUM") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Medium",_incarpad_) << " # medium reasonable" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_MEDIUM,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # " << DEFAULT_VASP_PREC_ENMAX_MEDIUM << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="NORMAL") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Normal",_incarpad_) << " # avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_NORMAL,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # " << DEFAULT_VASP_PREC_ENMAX_NORMAL << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl; 
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=High",_incarpad_) << " # avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_HIGH,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # " << DEFAULT_VASP_PREC_ENMAX_HIGH << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE") // || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH")
    { //CO20200106 - patching for auto-indenting
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Accurate",_incarpad_) << " # avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_ACCURATE,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # " << DEFAULT_VASP_PREC_ENMAX_ACCURATE << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LREAL=.FALSE.",_incarpad_) << " # reciprocal space projection technique" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFF=1E-6",_incarpad_) << " # high accuracy required" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ALGO=Fast",_incarpad_) << " # fast determination of ground state" << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="PHONONS") { 
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Accurate",_incarpad_) << " # avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_ACCURATE,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << " # " << DEFAULT_VASP_PREC_ENMAX_ACCURATE << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LREAL=.FALSE.",_incarpad_) << " # reciprocal space projection technique" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFF=1E-8",_incarpad_) << " # high accuracy required" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFFG=-0.001",_incarpad_) << " # high accuracy required" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.content_int),_incarpad_) << " # Many electronic steps" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELMIN=4",_incarpad_) << " # The forces have to be well converged" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ADDGRID=.TRUE.",_incarpad_) << " # For finer forces" << endl;
    };
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Precision PREC=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " [AFLOW] end" << endl;
    //  cerr << xvasp.INCAR.str() << endl;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Metagga
namespace KBIN {
  void XVASP_INCAR_Metagga(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION   TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE
    string FileContent,strline;
    int imax;
    // cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"METAGGA",TRUE) || aurostd::substring2bool(strline,"#METAGGA",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) {
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="TPSS") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=TPSS)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="RTPSS") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=RTPSS)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="M06L") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=M06L)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MBJL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MBJL)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="SCAN") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=SCAN)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS0") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MS0)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS1") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MS1)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS2") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MS2)" << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="NONE") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=NONE)" << endl;
        }
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme!="NONE") {
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Metagga METAGGA=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << " [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="TPSS") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=TPSS",_incarpad_) << " # METAGGA = TPSS" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="RTPSS") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=RTPSS",_incarpad_) << " # METAGGA = RTPSS" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="M06L") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=M06L",_incarpad_) << " # METAGGA = M06L" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MBJL") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MBJL",_incarpad_) << " # METAGGA = MBJL" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="SCAN") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=SCAN",_incarpad_) << " # METAGGA = SCAN" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS0") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MS0",_incarpad_) << " # METAGGA = MS0" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS1") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MS1",_incarpad_) << " # METAGGA = MS1" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS2") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MS2",_incarpad_) << " # METAGGA = MS2" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="NONE") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=NONE",_incarpad_) << " # METAGGA = NONE" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Metagga METAGGA=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << " [AFLOW] end" << endl;
    }
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Ivdw
namespace KBIN {
  void XVASP_INCAR_Ivdw(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION   number_for_VASP_see_manual_for_IVDW | 0
    string FileContent,strline;
    int imax;
    // cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"IVDW",TRUE) || aurostd::substring2bool(strline,"#IVDW",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) {
          xvasp.INCAR << "# " << strline << " # AFLOW REMOVED" << endl;
        }
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme!="0") {
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Ivdw IVDW=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << " [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="0" || vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="OFF" ||
          vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="NONE" || vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="") {
        xvasp.INCAR << aurostd::PaddedPOST("#IVDW=0",_incarpad_) << " # IVDW = 0" << endl;
      } else {	
        xvasp.INCAR << aurostd::PaddedPOST("IVDW="+vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme,_incarpad_) << " # IVDW = " << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << endl;
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Ivdw IVDW=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << " [AFLOW] end" << endl;
    }
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_ABMIX
namespace KBIN {
  void XVASP_INCAR_ABMIX(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    if(!vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry) return;
    string FileContent,strline;
    int imax;
    // cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"AMIX",TRUE) || aurostd::substring2bool(strline,"#AMIX",TRUE) ||
          aurostd::substring2bool(strline,"BMIX",TRUE) || aurostd::substring2bool(strline,"#BMIX",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_ABMIX)" << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }

    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing KBIN::XVASP_INCAR_ABMIX [AFLOW]" << endl;
    uint isUS=FALSE,isPAW=FALSE;
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme.at(0)=='A') { // AUTO
      // needs to check PPs
      //   cerr << "xvasp.POTCAR_TYPE=" << xvasp.POTCAR_TYPE << endl;
      if(xvasp.POTCAR_TYPE=="LDA") isUS=TRUE;
      if(xvasp.POTCAR_TYPE=="GGA") isUS=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_LDA") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_GGA") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_PBE") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_LDA_KIN") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_PBE_KIN") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_RPBE") isPAW=TRUE;
      if(!isUS && !isPAW) {
        isPAW=TRUE; // some default
        xvasp.INCAR << "# ABMIX=AUTO => POTCAR type not found, switching to PAW (aflow_ivasp.cpp)." << endl;   
      }      
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme=="US" || isUS) { // US
      if(isUS) xvasp.INCAR << "# ABMIX=AUTO => Found US pseudopotentials" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("AMIX=0.5",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX=1.0",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("AMIX_MAG=2.0",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX_MAG=1.0",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      // needs to check PPs
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme=="PAW" || isPAW) { // PAW
      if(isPAW) xvasp.INCAR << "# ABMIX=AUTO => Found PAW pseudopotentials" << endl;   
      xvasp.INCAR << aurostd::PaddedPOST("AMIX=0.2",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX=0.00001",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("AMIX_MAG=0.8",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX_MAG=0.00001",_incarpad_) << " # Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      // needs to check PPs
    }
    // check for LDA/GGA
    if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA") {
      //    xvasp.INCAR << "#Fixing RWIGS/LORBIG for LDA/GGA" << endl;
      ofstream oss;KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,oss,ON);
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_GetNBANDS
namespace KBIN {
  int XVASP_INCAR_GetNBANDS(const _xvasp& xvasp,const _aflags& aflags,bool ispin) {
    vector<double> vZVAL;
    double ZVAL=GetZVAL(xvasp.POTCAR,vZVAL);
    // // cout << XPID << "00000  MESSAGE POTCAR ZVAL max   = " << ZVAL << endl;
    int nbands=0;
    ispin=true; //CO20210315 - override for safety
    int NPAR=0,NCORE=1;
    KBIN::XVASP_Get_NPAR_NCORE(xvasp,aflags,NPAR,NCORE); //CO20210315
    nbands=GetNBANDS((int) ZVAL,(int) xvasp.str.atoms.size(),5,ispin,NPAR);
    // nbands=nbands+20+nbands/5; // MORE SAFETY
    // cout << XPID << "00000  MESSAGE POTCAR NBANDS = " << nbands << endl;
    if(!XHOST.QUIET) cout << XPID << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NBANDS  = " << nbands << endl;
    return nbands;
  }
}

// ***************************************************************************
// KBIN::XVASP_IALGO2ALGO
namespace KBIN {
  string INCAR_IALGO2ALGO(int ialgo){ //CO20210315
    //https://www.vasp.at/wiki/index.php/ALGO
    if(ialgo==38){return "Normal";}
    else if(ialgo==48){return "Veryfast";}
    else if(ialgo==58){return "All";}
    else if(ialgo==53){return "Damped";}
    else if(ialgo==90){return "Exact";}
    else if(ialgo==4){return "Subrot";}
    else if(ialgo==3){return "Eigenval";}
    else if(ialgo==2){return "None";}
    return "";
  }
} // namespace KBIN

namespace KBIN {
  bool XVASP_INCAR_Read_MAGMOM(_xvasp& xvasp){  //CO20210315
    //needs support for LNONCOLLINEAR=.TRUE., add later
    //order_parameter_value should probably become double, not necessary for now
    if(aurostd::kvpair2bool(xvasp.INCAR,"LNONCOLLINEAR","=") && aurostd::kvpair2string(xvasp.INCAR,"LNONCOLLINEAR","=")==".TRUE."){return false;}
    if(!aurostd::kvpair2bool(xvasp.INCAR,"MAGMOM","=")){return false;}
    string value=aurostd::kvpair2string(xvasp.INCAR,"MAGMOM","=");
    vector<string> vtokens,vtokens2;
    vector<int> vmags;
    aurostd::string2tokens(value,vtokens," ");
    uint i=0,j=0;
    int natoms=0;
    for(i=0;i<vtokens.size();i++){
      if(vtokens[i].find("*")==string::npos){
        if(!aurostd::isfloat(vtokens[i])){return false;}
        vmags.push_back(aurostd::string2utype<int>(vtokens[i]));
      }else{
        aurostd::string2tokens(vtokens[i],vtokens2,"*");
        if(vtokens2.size()!=2){return false;}
        if((!aurostd::isfloat(vtokens2[0]))||(!aurostd::isfloat(vtokens2[1]))){return false;}
        //assume natoms is first, not last
        natoms=aurostd::string2utype<int>(vtokens2[0]);
        if(natoms<1){return false;}
        for(j=0;j<(uint)natoms;j++){vmags.push_back(aurostd::string2utype<int>(vtokens2[1]));}
      }
    }
    if(vmags.size()!=xvasp.str.atoms.size()){return false;} //check before you apply
    for(i=0;i<vmags.size();i++){
      xvasp.str.atoms[i].order_parameter_atom=true;
      xvasp.str.atoms[i].order_parameter_value=vmags[i];
    }
    return true;
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::XVASP_INCAR_PREPARE_GENERIC
namespace KBIN {
  bool XVASP_INCAR_PREPARE_GENERIC(const string& command,_xvasp& xvasp,const _vflags& vflags,const string& svalue,int ivalue,double dvalue,bool OPTION){
    //CO20210315 - extensive rewrite
    //the schemes below check if the modification needs to be made (has it already been made?)
    //maintain this feedback system to ensure aflow doesn't keep spinning its wheels on the same fixes
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::XVASP_INCAR_PREPARE_GENERIC";
    string operation=function+"("+command+")";
    string operation_svalue=function+"("+command+"="+svalue+")";
    string operation_ivalue=function+"("+command+"="+aurostd::utype2string(ivalue)+")";
    string operation_dvalue=function+"("+command+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_)+")";
    string operation_option=function+"("+command+"="+(OPTION==ON?string("ON"):string("OFF"))+")";

    //this function only changes xvasp.INCAR
    //no pre-load of INCAR, i.e., assume current INCAR is inside xvasp.INCAR
    //no rewrite of INCAR

    //return true/false indicates whether INCAR was modified

    bool Krun=true; //indicates whether INCAR was modified
    bool VERBOSE=(FALSE || XHOST.vflag_control.flag("MONITOR_VASP")==false || LDEBUG);
    string keyword="",incar_input="",incar_input_old="",incar_comment="";

    // ***************************************************************************
    // ALGO ALGO ALGO ALGO ALGO ALGO
    if(command=="ALGO") {
      keyword=command;

      if(Krun){
        if(svalue.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"svalue empty (command="+command+")",_INPUT_ILLEGAL_);}
        if(svalue=="NORMAL")        {incar_input=keyword+"=Normal";   incar_comment="blocked Davidson block iteration scheme (IALGO=38)";}
        else if(svalue=="VERYFAST") {incar_input=keyword+"=Veryfast"; incar_comment="RMM-DIIS (IALGO=48)";}
        else if(svalue=="FAST")     {incar_input=keyword+"=Fast";     incar_comment="IALGO=38 is used for the initial phase, then VASP switches to IALGO=48";}
        else if(svalue=="ALL")      {incar_input=keyword+"=All";      incar_comment="all band simultaneous update of orbitals (IALGO=58)";}
        else if(svalue=="DAMPED")   {incar_input=keyword+"=Damped";   incar_comment="a damped velocity friction algorithm (IALGO=53)";}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"svalue="+svalue+" unknown (command="+command+")",_INPUT_ILLEGAL_);}

        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword+",IALGO",operation_svalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_svalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_svalue << " " << incar_comment << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_svalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // AMIN AMIN AMIN AMIN AMIN AMIN
    else if(command=="AMIN") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(std::signbit(dvalue)==false){  //use dvalue<0 to remove the key only
          if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
          xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
          if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
        }
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // AMIX AMIX AMIX AMIX AMIX AMIX
    else if(command=="AMIX") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(std::signbit(dvalue)==false){  //use dvalue<0 to remove the key only
          if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
          xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
          if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
        }
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM
    else if(command=="AUTO_MAGMOM") {
      keyword="MAGMOM";

      if(Krun){
        incar_input=keyword+"=";
        for(uint i=0;i<xvasp.str.atoms.size();i++) {
          incar_input+=(i>0?" ":"");
          if(xvasp.str.atoms[i].order_parameter_atom==true){ 
            incar_input+=aurostd::utype2string(xvasp.str.atoms[i].order_parameter_value);
          }else{  //xvasp.str.atoms[i].order_parameter_atom==false
            incar_input+="0";
          }
        }
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"="); //add space 
          if(incar_input==incar_input_old){Krun=false;}
        }
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation << " MAGMOM " << xvasp.str.atoms.size() << " atoms" << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // BMIX BMIX BMIX BMIX BMIX BMIX
    else if(command=="BMIX") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(std::signbit(dvalue)==false){  //use dvalue<0 to remove the key only
          if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
          xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
          if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
        }
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // EDIFF EDIFF EDIFF EDIFF EDIFF EDIFF
    else if(command=="EDIFF") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // EDIFFG EDIFFG EDIFFG EDIFFG EDIFFG EDIFFG
    else if(command=="EDIFFG") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // ENMAX ENMAX ENMAX ENMAX ENMAX ENMAX
    else if(command=="ENMAX") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY
    else if(command=="ENMAX_MULTIPLY") {
      keyword="ENMAX";

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue*xvasp.POTCAR_ENMAX,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << " " << dvalue << "*"+keyword+" (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials" << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // IALGO IALGO IALGO IALGO IALGO IALGO
    else if(command=="IALGO") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword+",ALGO",operation,VERBOSE);  //CO20200624 //CO20210315 - why IMIX too?
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation << endl;  //IALGOX1?
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // ICHARG ICHARG ICHARG ICHARG ICHARG ICHARG
    else if(command=="ICHARG") {  //ME20191028
      keyword=command;

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_svalue,VERBOSE);  //CO20200624

        // Use negative values to just remove ICHARG
        if(ivalue>=0) {
          incar_input=keyword+"="+aurostd::utype2string(ivalue);
          if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
            incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
            if(incar_input==incar_input_old){Krun=false;}
          }
          //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces

          if(Krun){
            const string& chgcar=svalue;
            if(VERBOSE) xvasp.INCAR << "# Performing " << operation_svalue << " CHGCAR_FILE=" << chgcar << " [AFLOW] begin" << std::endl;
            xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_svalue << std::endl;
            if(VERBOSE) xvasp.INCAR << "# Performing " << operation_svalue << " CHGCAR_FILE=" << chgcar << " [AFLOW] end" << std::endl;
          }
        }
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // IMIX IMIX IMIX IMIX IMIX IMIX
    else if(command=="IMIX") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // ISMEAR ISMEAR ISMEAR ISMEAR ISMEAR ISMEAR
    else if(command=="ISMEAR" || command=="ISMEAR_STATIC" || command=="ISMEAR_BANDS") { //CO20181129
      //https://www.vasp.at/wiki/index.php/ISMEAR
      keyword="ISMEAR";

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        if(ivalue==0)      {incar_comment="Gaussian smearing";}
        else if(ivalue==1) {incar_comment="Methfessel-Paxton smearing order 1";}
        else if(ivalue==2) {incar_comment="Methfessel-Paxton smearing order 2";}
        else if(ivalue==-1){incar_comment="Fermi smearing";}
        else if(ivalue==-4){incar_comment="tetrahedron smearing";}
        else if(ivalue==-5){incar_comment="tetrahedron smearing with Bloechl corrections";}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ivalue="+aurostd::utype2string(ivalue)+" unknown (command="+command+")",_INPUT_ILLEGAL_);}

        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_ivalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_ivalue << " " << incar_comment << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // LASPH LASPH LASPH LASPH LASPH LASPH
    else if(command=="LASPH") {
      //PAW_CORRECTIONS
      //http://cms.mpi.univie.ac.at/vasp/vasp/LASPH_tag.html#incar-lasph
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+(OPTION==ON?string(".TRUE."):string(".FALSE."));
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,"LASPH",operation_option,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE)  xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option <<  endl;
        if(VERBOSE)  xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // LCHARG LCHARG LCHARG LCHARG LCHARG LCHARG
    else if(command=="LCHARG") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+(OPTION==ON?string(".TRUE."):string(".FALSE."));
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_option,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // LREAL LREAL LREAL LREAL LREAL LREAL
    else if(command=="LREAL") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+(OPTION==ON?string(".TRUE."):string(".FALSE."));
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_option,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // LS_COUPLING LS_COUPLING LS_COUPLING LS_COUPLING LS_COUPLING LS_COUPLING
    else if(command=="LS_COUPLING") {
      keyword="LSORBIT";

      //do a check that we're not already doing the ls-coupling-type calculation, LSORBIT is a good enough check
      if(Krun){
        incar_input=keyword+"="+(OPTION==ON?".TRUE.":".FALSE.");
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
      }

      if(Krun){
        //REMOVE LINES
        // strip magmom from incar and make the magmom from the order parameter
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword+",LNONCOLLINEAR,MAGMOM",operation_option,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << endl;
        if(OPTION==ON){
          xvasp.INCAR << aurostd::PaddedPOST("LNONCOLLINEAR=.TRUE.",_incarpad_) << " # " << operation_option << endl;
          incar_input="MAGMOM=";
          for(uint i=0;i<xvasp.str.atoms.size();i++){
            incar_input+=(i>0?" ":"");
            incar_input+="0 0 5";
          }
          xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << " MAGMOM 3*" << xvasp.str.atoms.size() << " atoms" << endl;
        }else{ //OPTION==OFF
          xvasp.INCAR << aurostd::PaddedPOST("LNONCOLLINEAR=.FALSE.",_incarpad_) << " # " << operation_option << endl;
        }
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // NBANDS NBANDS NBANDS NBANDS NBANDS NBANDS
    else if(command=="NBANDS") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_ivalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_ivalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] end" << endl;
      }
    }

    // ***************************************************************************
    // NELM NELM NELM NELM NELM NELM
    else if(command=="NELM" || command=="NELM_STATIC") { //CO20200624
      keyword="NELM";

      if(Krun){
        bool nelm_fix=true; //[CO20210315 - always patch, you might have NELM in INCAR section of aflow.in](ivalue!=AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM); //CO20200624 - default for VASP is 60, don't add the line if unnecessary
        if(nelm_fix==false){Krun=false;}
      }

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_ivalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_ivalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // NPAR NPAR NPAR NPAR NPAR NPAR
    else if(command=="NPAR") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_ivalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_ivalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // NSW NSW NSW NSW NSW NSW
    else if(command=="NSW") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(ivalue);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_ivalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_ivalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_ivalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // POTIM POTIM POTIM POTIM POTIM POTIM
    else if(command=="POTIM") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // PSTRESS PSTRESS PSTRESS PSTRESS PSTRESS PSTRESS
    else if(command=="PSTRESS") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // SIGMA SIGMA SIGMA SIGMA SIGMA SIGMA
    else if(command=="SIGMA" || command=="SIGMA_STATIC" || command=="SIGMA_BANDS") { //CO20181129
      keyword="SIGMA";

      if(Krun){
        incar_input=keyword+"="+aurostd::utype2string(dvalue,_IVASP_DOUBLE2STRING_PRECISION_);
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_dvalue,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_dvalue << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_dvalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // SPIN SPIN SPIN SPIN SPIN SPIN
    else if(command=="SPIN") {
      keyword="ISPIN";

      //do a check that we're not already doing the spin-type calculation, ISPIN is a good enough check
      if(Krun){
        incar_input=keyword+"="+(OPTION==ON?"2":"1");
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword+",ISPIND",operation_option,VERBOSE);  //CO20200624
        bool remove_magmom=false;
        if(OPTION==OFF){remove_magmom=true;}
        else{ //OPTION==ON
          //if AUTO_MAGMOM, remove MAGMOM from INCAR
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry && vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option){remove_magmom=true;}
        }
        if(remove_magmom){XVASP_INCAR_REMOVE_ENTRY(xvasp,"MAGMOM",operation_option,VERBOSE);}  //CO20200624
        //otherwise, check if has already been specified (e.g., in INCAR section of _AFLOWIN_)
        bool magmom_already_specified=aurostd::kvpair2bool(xvasp.INCAR,"MAGMOM","=");  //CO20210315 - must use kvpair2bool instead of substring2bool() since magmom values have spaces
        //ADD LINES
        if(VERBOSE)  xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << endl;
        if(OPTION==ON){xvasp.INCAR << aurostd::PaddedPOST("ISPIND=2",_incarpad_) << " # " << operation_option << endl;}
        if(OPTION==ON && magmom_already_specified==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option==FALSE) { //CO
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry && vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option){  // with spin, make the mag mom +5
            xvasp.INCAR << aurostd::PaddedPOST("MAGMOM="+aurostd::utype2string(xvasp.str.atoms.size())+"*5",_incarpad_) << " # " << operation_option << " " << xvasp.str.atoms.size() << " atom" << (xvasp.str.atoms.size()>1?"s":"") << endl;
          }else{  // otherwise, make the magmom +1
            xvasp.INCAR << aurostd::PaddedPOST("MAGMOM="+aurostd::utype2string(xvasp.str.atoms.size())+"*1",_incarpad_) << " # " << operation_option << " " << xvasp.str.atoms.size() << " atom" << (xvasp.str.atoms.size()>1?"s":"") << " (default)" << endl;
          }
        }
        if(VERBOSE)  xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // SYM SYM SYM SYM SYM SYM
    else if(command=="SYM") {
      //https://www.vasp.at/wiki/index.php/ISYM
      keyword="ISYM";

      if(Krun){
        if(OPTION==ON){incar_input=keyword+"=2";incar_comment="SYMMETRY=ON";}
        else{ //OPTION==OFF
          //expand this list as needed in the future: 
          if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry && vflags.KBIN_VASP_FORCE_OPTION_SPIN.option){incar_input=keyword+"=-1";incar_comment="SYMMETRY=OFF (MAGNETIC)";} //SPIN==ON
          else{incar_input=keyword+"=0";incar_comment="SYMMETRY=OFF";}  //SPIN==OFF
        }
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_option,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << " " << incar_comment << endl; //CO20181113
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // SYMPREC SYMPREC SYMPREC SYMPREC SYMPREC SYMPREC
    else if(command=="SYMPREC") {
      keyword=command;

      if(Krun){
        incar_input=keyword+"=1E-7";  //no other value worth trying. hacking SYMPREC is a cheap way to avoid the error, it's better to get to more robust solutions
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //CO20210315 - remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation << endl;  //no real point trying other values, better to try another fix
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // TYPE TYPE TYPE TYPE TYPE TYPE
    else if(command=="TYPE") {
      keyword="ISMEAR";
      string keyword2="SIGMA",incar_input2="";

      if(Krun){
        if(svalue=="DEFAULT") {incar_input=keyword+"=1";incar_input2=keyword2+"=0.1";incar_comment="for default (as metal)";}
        else if(svalue=="METAL") {incar_input=keyword+"=1";incar_input2=keyword2+"=0.1";incar_comment="as metal";}
        else if(svalue=="SEMICONDUCTOR" || svalue=="INSULATOR") {incar_input=keyword+"=0";incar_input2=keyword2+"=0.05";incar_comment="for insulators/semiconductors";}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"svalue="+svalue+" unknown (command="+command+")",_INPUT_ILLEGAL_);}
        bool Krun1=true,Krun2=true;
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun1=false;}
        }
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword2,"=")){
          incar_input_old=keyword2+"="+aurostd::kvpair2string(xvasp.INCAR,keyword2,"=");
          if(incar_input2==incar_input_old){Krun2=false;}
        }
        if(Krun1==false&&Krun2==false){Krun=false;}
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword+","+keyword2,operation,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_svalue << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation << " " << incar_comment << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input2,_incarpad_) << " # " << operation << " " << incar_comment << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_svalue << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // WAVECAR WAVECAR WAVECAR WAVECAR WAVECAR WAVECAR
    else if(command=="WAVECAR") {
      keyword="LWAVE";

      if(Krun){
        incar_input=keyword+"="+(OPTION==ON?string(".TRUE."):string(".FALSE."));
        if(aurostd::kvpair2bool(xvasp.INCAR,keyword,"=")){
          incar_input_old=keyword+"="+aurostd::kvpair2string(xvasp.INCAR,keyword,"=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //CO20210315 - remove whitespaces
      }

      if(Krun){
        //REMOVE LINES
        XVASP_INCAR_REMOVE_ENTRY(xvasp,keyword,operation_option,VERBOSE);  //CO20200624
        //ADD LINES
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] begin" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation_option << endl;
        if(VERBOSE) xvasp.INCAR << "# Performing " << operation_option << " [AFLOW] end" << endl;
      }
    }
    // ***************************************************************************

    // ***************************************************************************
    // UNKNOWN UNKNOWN UNKNOWN UNKNOWN UNKNOWN UNKNOWN
    else {
      stringstream message;
      message << "command=" << command << " not found" << endl;
      message << "svalue=" << svalue << endl;
      message << "ivalue=" << ivalue << endl;
      message << "dvalue=" << dvalue << endl;
      message << "OPTION=" << OPTION << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
    }
    // ***************************************************************************

    // ***************************************************************************

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " RETURNING CORRECTLY" << endl;
      cerr << __AFLOW_FUNC__ << " command=" << command << endl;
      cerr << __AFLOW_FUNC__ << " svalue=" << svalue << endl;
      cerr << __AFLOW_FUNC__ << " ivalue=" << ivalue << endl;
      cerr << __AFLOW_FUNC__ << " dvalue=" << dvalue << endl;
      cerr << __AFLOW_FUNC__ << " OPTION=" << OPTION << endl;
    }

    if(Krun==false){return false;}

    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    return true;
  }
} 

// ***************************************************************************
// KBIN::XVASP_INCAR_ADJUST_ICHARG
namespace KBIN {
  //ME20191028
  // When a CHGCAR file is specified in the aflow.in file, it is really only
  // useful in the first relaxation calculations. For all other calculations,
  // it should not read from that CHGCAR file, so set ICHARG = 2 (default, but do not write) and comment
  // out the original CHGCAR file. Exception: if a CHGCAR file is output
  // after the relaxation, assume that the user wants to reuse it.
  void XVASP_INCAR_ADJUST_ICHARG(_xvasp& xvasp, _vflags& vflags, _aflags& aflags, int step, bool write_incar, ofstream& FileMESSAGE) {  //CO20210315 - write_input
    string function="KBIN::XVASP_INCAR_ADJUST_ICHARG";
    if ((step == 1) && vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry) {
      // Do not set ICHARG when a CHGCAR file is output
      if (!vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option) {
        ostringstream aus;
        aus << "00000  MESSAGE ICHARG: Removing ICHARG - " << Message(__AFLOW_FILE__,aflags) << std::endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("ICHARG", xvasp, vflags, vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.content_string, -1, 0.0, false); //remove any ICHARG from INCAR, leaving default ==2
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed", true);
      }
      if (xvasp.aopts.flag("FLAG::XVASP_INCAR_changed")) {
        xvasp.aopts.flag("FLAG::XVASP_INCAR_generated", true);
        aurostd::StringstreamClean(xvasp.INCAR_orig);
        xvasp.INCAR_orig << xvasp.INCAR.str();
        if(write_incar){  //[CO20210315 - we don't know if there's a static run afterwards]if (step < xvasp.NRELAX)  // Do not write when at the last step or there will be an extra INCAR
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory+"/INCAR"));
        }
      }

      // Comment out CHGCAR file in aflow.in
      if (vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry) {
        KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]CHGCAR_FILE=",function);
        //CO20210315 - ME is this right? mimicking what we do in XVASP_INCAR_SPIN_REMOVE_RELAX()
        vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option=false;
        vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry=false;
        //[CO20210315 - OBSOLETE]stringstream aflowin_fixed;
        //[CO20210315 - OBSOLETE]string filename = aurostd::CleanFileName(xvasp.Directory + "/" + _AFLOWIN_);
        //[CO20210315 - OBSOLETE]string filecontent = aurostd::file2string(aurostd::CleanFileName(filename));
        //[CO20210315 - OBSOLETE]int nlines = aurostd::GetNLinesString(filecontent);
        //[CO20210315 - OBSOLETE]string line = "";
        //[CO20210315 - OBSOLETE]for (int l = 1; l <= nlines; l++) {
        //[CO20210315 - OBSOLETE]  line = aurostd::GetLineString(filecontent, l);
        //[CO20210315 - OBSOLETE]  if (aurostd::substring2bool("[VASP_FORCE_OPTION]CHGCAR_FILE=", line)) {
        //[CO20210315 - OBSOLETE]    string line_fixed = aurostd::RemoveWhiteSpacesFromTheFront(line);
        //[CO20210315 - OBSOLETE]    if (line_fixed[0] != '#') line = "#" + line_fixed;
        //[CO20210315 - OBSOLETE]  }
        //[CO20210315 - OBSOLETE]  aflowin_fixed << line << std::endl;
        //[CO20210315 - OBSOLETE]}
        //[CO20210315 - OBSOLETE]aurostd::stringstream2file(aflowin_fixed, filename);
      }
    }
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX
namespace KBIN {
  void XVASP_INCAR_SPIN_REMOVE_RELAX(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,bool write_incar,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315 - write_input
    string function="KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX";
    ostringstream aus;
    bool fix_aflowin=FALSE;
    _kflags kflags;
    if(step==1 && vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 && abs(xvasp.str.qm_mag_atom)<DEFAULT_VASP_SPIN_REMOVE_CUTOFF) {
      aus << "00000  MESSAGE FORCE SPIN_OFF: [VASP_FORCE_OPTION]SPIN_REMOVE_RELAX_1  SPIN= " << xvasp.str.qm_mag_atom << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SPIN",xvasp,vflags,"",0,0.0,OFF);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      fix_aflowin=TRUE;
    }
    if(step==2 && vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 && abs(xvasp.str.qm_mag_atom)<DEFAULT_VASP_SPIN_REMOVE_CUTOFF) {
      aus << "00000  MESSAGE FORCE SPIN_OFF: [VASP_FORCE_OPTION]SPIN_REMOVE_RELAX_2  SPIN= " << xvasp.str.qm_mag_atom << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SPIN",xvasp,vflags,"",0,0.0,OFF);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      fix_aflowin=TRUE;
    }
    if(xvasp.aopts.flag("FLAG::XVASP_INCAR_changed")) {
      xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);
      aurostd::StringstreamClean(xvasp.INCAR_orig); xvasp.INCAR_orig << xvasp.INCAR.str();
      if(write_incar){ //[CO20210315 - we don't know if there's a static run afterwards]if (step < xvasp.NRELAX)  //ME20200107 - do not write when at the last step or there will be an extra INCAR
        aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
      }
      // xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR"); // DID REREAD
    }
    if(fix_aflowin) {
      // fix aflow.in
      KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]SPIN=",function); //CO20210315
      KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]SPIN=OFF",function);  //CO20210315
      //[CO20210315 - OBSOLETE]stringstream aus_exec;
      //[CO20210315 - OBSOLETE]aus_exec << "cd " << xvasp.Directory << endl;
      //[CO20210315 - OBSOLETE]aus_exec << "cat " << _AFLOWIN_ << " | sed \"s/\\[VASP_FORCE_OPTION\\]SPIN=/#\\[VASP_FORCE_OPTION\\]SPIN=/g\" | sed \"s/##\\[/#\\[/g\" > aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << endl;
      //[CO20210315 - OBSOLETE]aus_exec << "echo \"[VASP_FORCE_OPTION]SPIN=OFF" << "      // Self Correction\"" << " >> " << _AFLOWIN_ << endl;
      //[CO20210315 - OBSOLETE]aurostd::execute(aus_exec);
      //ME20181117 - Also change vflags
      vflags.KBIN_VASP_FORCE_OPTION_SPIN.option = false;
      vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = false;
      vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = false;
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_KPOINTS_IBZKPT_UPDATE
namespace KBIN {
  void XVASP_KPOINTS_IBZKPT_UPDATE(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,bool write_incar,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210315 - write_input
    ostringstream aus;
    // aus << "00000  MESSAGE IBZKPT.relax1 - step=" << step << " " << Message(__AFLOW_FILE__,aflags) << endl;
    // aus << "00000  MESSAGE IBZKPT.relax1 - vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT") << " " << Message(__AFLOW_FILE__,aflags) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    if(step>=1 && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT")) {
      aus << "00000  MESSAGE FORCE KPOINTS from IBZKPT.relax1 - " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aurostd::StringstreamClean(xvasp.KPOINTS);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      aurostd::file2stringstream(string(xvasp.Directory+"/IBZKPT.relax1"),xvasp.KPOINTS);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed")) {
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);
      aurostd::StringstreamClean(xvasp.KPOINTS_orig); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
      if(write_incar){  //[CO20210315 - we don't know if there's a static run afterwards]if (step < xvasp.NRELAX)  //ME20200107 - do not write when at the last step or there will be an extra INCAR
        aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      }
      // xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS"); // DID REREAD
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_OFF
namespace KBIN {
  void XVASP_INCAR_LDAU_OFF(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(!aurostd::substring2bool(strline,"SYSTEM", TRUE) &&
          (aurostd::substring2bool(strline,"LDAU",TRUE) || aurostd::substring2bool(strline,"#LDAU",TRUE) ||
           aurostd::substring2bool(strline,"LMAXMIX",TRUE) || aurostd::substring2bool(strline,"#LMAXMIX",TRUE))) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_LDAU_OFF)" << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    } // xvasp.INCAR << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=OFF [AFLOW] begin" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LDAU=.FALSE.",_incarpad_) << " # LDAU=OFF" << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=OFF [AFLOW] end" << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_ON
namespace KBIN {
  void XVASP_INCAR_LDAU_ON(_xvasp& xvasp,_vflags& vflags,uint type) {        // AFLOW_FUNCTION_IMPLEMENTATION
    stringstream message;
    /// is LDAU necessary ?
    if(type!=1 && type!=2) return; // NO LDAU
    // YES
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(!aurostd::substring2bool(strline,"SYSTEM", TRUE) &&
          !aurostd::substring2bool(strline,"LDAUL",TRUE) && !aurostd::substring2bool(strline,"#LDAUL",TRUE) &&
          !aurostd::substring2bool(strline,"LDAUU",TRUE) && !aurostd::substring2bool(strline,"#LDAUU",TRUE) &&
          !aurostd::substring2bool(strline,"LDAUJ",TRUE) && !aurostd::substring2bool(strline,"#LDAUJ",TRUE) &&
          (aurostd::substring2bool(strline,"LDAU",TRUE) || aurostd::substring2bool(strline,"#LDAU",TRUE) ||
           aurostd::substring2bool(strline,"LMAXMIX",TRUE) || aurostd::substring2bool(strline,"#LMAXMIX",TRUE))) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_LDAU" << type << "_ON)" << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    } // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing LDAU" << type << "=ON [AFLOW] begin" << endl;
    int LMAXMIX=6;
    xvasp.INCAR << aurostd::PaddedPOST("LDAU=.TRUE.",_incarpad_) << " # AFLOW LSDA+U" << endl;
    // LDAU generation

    vector<string> vLDAUspecies;vector<uint> vLDAUtype;vector<int> vLDAUL;vector<double> vLDAUU,vLDAUJ;
    bool LDAU=FALSE;vLDAUtype.clear();vLDAUL.clear();vLDAUU.clear();vLDAUJ.clear();
    uint nspecies=0;
    vector<string> tokens_group,tokens;

    if(vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag==TRUE) { // get them from the names
      // cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag==TRUE => get them from the names " << endl;
      if(vflags.KBIN_VASP_LDAU_SPECIES.length()>1) {
        aurostd::string2tokens(vflags.KBIN_VASP_LDAU_SPECIES,tokens," ");
        nspecies=tokens.size();
        if(xvasp.str.species.size()!=nspecies) {
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"xvasp.str.species.size() != nspecies=tokens.size()",_INDEX_MISMATCH_); //CO20200624
        }
        for(uint i=0;i<xvasp.str.species.size();i++)
          AVASP_Get_LDAU_Parameters(aurostd::CleanStringASCII(KBIN::VASP_PseudoPotential_CleanName(tokens.at(i))),LDAU,vLDAUspecies,vLDAUtype,vLDAUL,vLDAUU,vLDAUJ); // parameters for LDAU2
      }
    } else { // get them from PARAMETERS
      // cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag==FALSE => get them from the parameters:= " << vflags.KBIN_VASP_LDAU_PARAMETERS << endl;
      LDAU=TRUE;
      aurostd::string2tokens(vflags.KBIN_VASP_LDAU_PARAMETERS,tokens_group,";");
      if(tokens_group.size()<3) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"must specify 4 entries in "+vflags.KBIN_VASP_LDAU_PARAMETERS,_INPUT_NUMBER_); //CO20200624
      }
      aurostd::string2tokens(tokens_group.at(0),tokens,",");
      nspecies=tokens.size();
      for(uint j=0;j<nspecies;j++) { vLDAUspecies.push_back("");vLDAUtype.push_back(2);vLDAUL.push_back(-1);vLDAUU.push_back(0.0);vLDAUJ.push_back(0.0);} // make space
      for(uint i=0;i<tokens_group.size();i++) {
        // cerr << i << " " << tokens_group.at(i) << endl;
        aurostd::string2tokens(tokens_group.at(i),tokens,",");
        if(nspecies!=tokens.size()) {
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"must specify "+aurostd::utype2string(nspecies)+" entries in "+tokens_group.at(i),_INPUT_NUMBER_); //CO20200624
        }
        //	for(uint j=0;j<nspecies;j++) cerr << tokens.at(j) << " "; cerr << endl;
        if(i==0) for(uint j=0;j<nspecies;j++) vLDAUspecies.at(j)=tokens.at(j);
        if(i==1) for(uint j=0;j<nspecies;j++) vLDAUL.at(j)=aurostd::string2utype<int>(tokens.at(j));
        if(i==2) for(uint j=0;j<nspecies;j++) vLDAUU.at(j)=aurostd::string2utype<double>(tokens.at(j));
        if(i==3) for(uint j=0;j<nspecies;j++) vLDAUJ.at(j)=aurostd::string2utype<double>(tokens.at(j));
      }
    }
    // get type
    if(xvasp.str.species.size()!=nspecies) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"xvasp.str.species.size() != nspecies=tokens.size()",_INPUT_NUMBER_); //CO20200624
    }

    bool type2=FALSE,type1=FALSE;
    for(uint i=0;i<xvasp.str.species.size();i++) {
      if(vLDAUtype.at(i)==0) {;} // do nothing
      if(vLDAUtype.at(i)==1) type1=TRUE;
      if(vLDAUtype.at(i)==2) type2=TRUE;
    }
    if(type1==FALSE && type2==FALSE)             // no LDAU
      cerr << __AFLOW_FUNC__ << " no need for LDAU, all species do not need or are not in the table : " << xvasp.str.title << endl;
    if(type1==TRUE && type2==FALSE && type!=1)   // LDAU type 1
      cerr << __AFLOW_FUNC__ << " your type=" << type << " is incompatible with type1 of table : " << xvasp.str.title << endl;
    if(type1==FALSE && type2==TRUE && type!=2)   // LDAU type 2
      cerr << __AFLOW_FUNC__ << " your type=" << type << " is incompatible with type2 of table : " << xvasp.str.title << endl;
    if(type1==TRUE && type2==TRUE) {
      message << "Your type=" << type << " is incompatible with type1 and type2 of table : " << xvasp.str.title << endl;
      message << "  You can NOT HAVE type1 and type2 simultaneously from table.... My suggestion is that you remove LDAU_SPECIES from " << _AFLOWIN_ << "" << endl;
      message << "  and specify LDAUL, LDAUU, LDAUJ in the INCAR part of aflow,in manually" << endl;
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_); //CO20200624
    }

    LMAXMIX=0;
    for(uint i=0;i<vLDAUL.size();i++) if(2*vLDAUL.at(i)>LMAXMIX) LMAXMIX=2*vLDAUL.at(i);        // get the MAX MIX for the orbitals
    // this is to max all the orbitals... if necessary..
    if(0) { //WSETYAWAN DISAGREES // do I need to max the orbitals ?
      for(uint i=0;i<vLDAUL.size();i++) if(vLDAUL.at(i)>vLDAUL.at(0)) vLDAUL.at(0)=vLDAUL.at(i);  // MAX the orbitals
      for(uint i=0;i<vLDAUL.size();i++) vLDAUL.at(i)=vLDAUL.at(0);                                // MAX the orbitals
    }
    // FIX NEGATIVES
    for(uint i=0;i<vLDAUL.size();i++) if(vLDAUL.at(i)<0) vLDAUL.at(i)=0;

    if(1) {
      stringstream straus;
      straus.clear();straus.str(std::string());
      straus << "#LDAU_SPECIES=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUspecies.at(i) << " ";
      xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # LDAU species" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int==0) {
        straus.clear();straus.str(std::string());
        straus << "LDAUL=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUL.at(i) << " ";
        xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # l-quantum number for which the on site interaction is added (Default 2) automatic LDAUL table" << endl;
        straus.clear();straus.str(std::string());
        straus << "LDAUU=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUU.at(i) << " ";
        xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # UEFF parameter. Automatic LDAUU table" << endl;
        straus.clear();straus.str(std::string());
        straus << "LDAUJ=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUJ.at(i) << " ";
        xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # J parameter (if used). Automatic LDAUJ table" << endl;
      } else {
        xvasp.INCAR << "#ADIABATIC_LDAU table = {nsteps,nspecies,{{LUJ}_nspecies}_steps}" << endl;
        straus.clear();straus.str(std::string());
        straus << "#ADIABATIC_LDAU="<< vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int << "," << nspecies << ",";
        for(int j=1;j<=vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int;j++) {
          for(uint i=0;i<nspecies;i++) {
            double cff=0.0;
            if(vLDAUL.at(i)==0) cff=0.0;                                                        // s orbital no U
            if(vLDAUL.at(i)==1) cff=0.0;                                                        // p orbital no U
            if(vLDAUL.at(i)==2) cff=(double) ((j+1.5)/(double) (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int));  // d orbital U
            if(vLDAUL.at(i)==3) cff=(double) ((j+1.5)/(double) (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int));  // f orbital U
            if(cff<0.0) cff=0.0;
            if(cff>1.0) cff=1.0;
            straus << vLDAUL.at(i) << "," << vLDAUU.at(i)*cff << "," <<  vLDAUJ.at(i)*cff;
            if(j<vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int || i<nspecies-1) straus << ",";
          }
        }
        xvasp.INCAR << straus.str() << endl; //  << " # ADIABATIC LDAU table = nsteps,nspecies,{{LUJ}_nspecies}_steps" << endl;
      }
    }

    if(0) {
      xvasp.INCAR << "#LDAU_SPECIES=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUspecies.at(i) << " ";xvasp.INCAR << "# LDAU species" << endl;
      xvasp.INCAR << "LDAUL=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUL.at(i) << " ";xvasp.INCAR << "# l-quantum number for which the on site interaction is added (Default 2) automatic LDAUL table" << endl;
      xvasp.INCAR << "LDAUU=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUU.at(i) << " ";xvasp.INCAR << "# UEFF parameter. Automatic LDAUU table" << endl;
      xvasp.INCAR << "LDAUJ=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUJ.at(i) << " ";xvasp.INCAR << "# J parameter (if used). Automatic LDAUJ table" << endl;
    }

    xvasp.INCAR << aurostd::PaddedPOST("LDAUTYPE="+aurostd::utype2string(type),_incarpad_) << " # Type of LDA+U." << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LMAXMIX="+aurostd::utype2string(LMAXMIX),_incarpad_) << " # Controls up to which l-quantum number the onsite PAW charge densities are passed through the charge density mixer." << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("LDAUPRINT=1",_incarpad_) << " # Controls verbosity of the L(S)DA+U module. (Default 0) # AFLOW LSDA+U" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LDAUPRINT=0",_incarpad_) << " # Controls verbosity of the L(S)DA+U module. (Default 0) # AFLOW LSDA+U" << endl;

    // cerr << xvasp.INCAR.str() << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing LDAU" << type << "=ON [AFLOW] end" << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_ADIABATIC
namespace KBIN {
  void XVASP_INCAR_LDAU_ADIABATIC(_xvasp& xvasp,int step) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    // reload incar
    xvasp.INCAR_orig.clear(); xvasp.INCAR_orig.str(xvasp.INCAR.str());
    xvasp.INCAR.clear(); xvasp.INCAR.str(aurostd::file2string(xvasp.Directory+"/INCAR"));
    //  xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    // operate
    string FileContent=xvasp.INCAR.str(),strline;
    int imax=aurostd::GetNLinesString(FileContent);
    stringstream straus;
    straus.clear();straus.str(std::string());

    for(int iline=1;iline<=imax;iline++) {
      strline=aurostd::GetLineString(FileContent,iline);
      if(!aurostd::substring2bool(strline,"LDAUL=",TRUE) && !aurostd::substring2bool(strline,"LDAUU=",TRUE) && !aurostd::substring2bool(strline,"LDAUJ=",TRUE)) {
        xvasp.INCAR << strline << endl;
        if(strline.find("#ADIABATIC_LDAU=")!=string::npos) {
          // FOUND write LDAUL LDAUU LDAUJ
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_ADIABATIC" << endl;
          vector<string> tokens1,tokens;
          aurostd::string2tokens(strline,tokens,"=");
          aurostd::string2tokens(tokens.at(1),tokens1,"#");
          aurostd::string2tokens(tokens1.at(0),tokens,",");
          uint t=0; // start of tokens. j
          uint nsteps=aurostd::string2utype<int>(tokens.at(t++));
          uint nspecies=aurostd::string2utype<int>(tokens.at(t++));
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_STEP=" << step << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_ADIABATIC_NSTEPS=" << nsteps << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_ADIABATIC_NSPECIES=" << nspecies << endl;
          if(tokens.size() != (uint) nsteps*nspecies*3+2) {
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"wrong number of ADIABATIC_LDAU entries ("+aurostd::utype2string(tokens.size())+")",_INPUT_NUMBER_); //CO20200624
          }
          vector<int> vLDAUL;vector<double> vLDAUU,vLDAUJ;
          for(uint k=1;k<=nsteps;k++) {
            for(uint j=0;j<nspecies;j++) {
              if((int) k==step) {
                vLDAUL.push_back(aurostd::string2utype<int>(tokens.at(t++)));
                vLDAUU.push_back(aurostd::string2utype<double>(tokens.at(t++)));
                vLDAUJ.push_back(aurostd::string2utype<double>(tokens.at(t++)));
              } else {
                t+=3;
              }
            }
          }
          straus.clear();straus.str(std::string());
          straus << "LDAUL=";for(uint j=0;j<nspecies;j++) straus <<  vLDAUL.at(j) << " ";
          xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # l-quantum number for which the on site interaction is added (Default 2) automatic LDAUL table" << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] " << aurostd::PaddedPOST(straus.str(),_incarpad_) << endl;
          straus.clear();straus.str(std::string());
          straus << "LDAUU=";for(uint j=0;j<nspecies;j++) straus <<  vLDAUU.at(j) << " ";
          xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # UEFF parameter. Automatic LDAUU table" << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] " << aurostd::PaddedPOST(straus.str(),_incarpad_) << endl;
          straus.clear();straus.str(std::string());
          straus << "LDAUJ=";for(uint j=0;j<nspecies;j++) straus <<  vLDAUJ.at(j) << " ";
          xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << " # J parameter (if used). Automatic LDAUJ table" << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] " << aurostd::PaddedPOST(straus.str(),_incarpad_) << endl;
        }
      }
    }
    // xvasp.INCAR << endl;
    // rewrite incar
    aurostd::RemoveFile(string(xvasp.Directory+"/INCAR"));
    aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_CUTOFF
namespace KBIN {
  void XVASP_INCAR_LDAU_CUTOFF(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    aurostd::StringstreamClean(xvasp.INCAR);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(!aurostd::substring2bool(strline,"SYSTEM", TRUE) &&
          (aurostd::substring2bool(strline,"LDAU",TRUE) || aurostd::substring2bool(strline,"#LDAU",TRUE))) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_LDAU_CUTOFF)" << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=CUTOFF [AFLOW] begin" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LDAU=.FALSE.",_incarpad_) << " # LDAU=CUTOFF" << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=CUTOFF [AFLOW] end" << endl;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_KPOINTS_Dielectric  DIELECTRIC
namespace KBIN {
  void XVASP_INCAR_KPOINTS_Dielectric_SET(_xvasp& xvasp,_kflags &kflags,_vflags& vflags,string svalue) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    // STATIC
    if(svalue=="STATIC" || svalue=="static") {
      //.  Retain the following static run entries and their values: ALGO, LREAL, NSIM, ISYM, IBRION, NSW, NELM, NELMIN, ENMAX, ISPIN, ISMEAR, SIGMA, and everything LDA+U related.
      // d.  Eliminate PSTRESS, EMIN, EMAX, LORBIT, ISIF, NEDOS.
      // NELM = 0

      string FileContent,strline;
      FileContent=xvasp.INCAR.str();
      aurostd::StringstreamClean(xvasp.INCAR);xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      int imax=aurostd::GetNLinesString(FileContent);
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) // Camilo
            || aurostd::substring2bool(strline,"ISIF",TRUE) || aurostd::substring2bool(strline,"#ISIF",TRUE) // Camilo
            || aurostd::substring2bool(strline,"EMIN",TRUE) || aurostd::substring2bool(strline,"#EMIN",TRUE) // Camilo
            || aurostd::substring2bool(strline,"EMAX",TRUE) || aurostd::substring2bool(strline,"#EMAX",TRUE) // Camilo
            // ||  aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) // Camilo
            || aurostd::substring2bool(strline,"ICHARG",TRUE) || aurostd::substring2bool(strline,"#ICHARG",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LEPSILON",TRUE) || aurostd::substring2bool(strline,"#LEPSILON",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LCALCEPS",TRUE) || aurostd::substring2bool(strline,"#LCALCEPS",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LRPA",TRUE) || aurostd::substring2bool(strline,"#LRPA",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT",TRUE)  // Camilo
            || aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE)  // Camilo
            || aurostd::substring2bool(strline,"NPAR",TRUE) || aurostd::substring2bool(strline,"#NPAR",TRUE)   // NPAR UNCOMPATIBLE WITH Dielectric
            // hope...	 || aurostd::substring2bool(strline,"NBANDS",TRUE) || aurostd::substring2bool(strline,"#NBANDS",TRUE) || aurostd::substring2bool(strline,"#fixed nbands",TRUE)   // NPAR UNCOMPATIBLE WITH Dielectric
          ) { // Camilo
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET)" << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_STATIC [AFLOW] begin" << endl;
      // xvasp.INCAR << aurostd::PaddedPOST("NELM=0",_incarpad_) << " # Performing DIELECTRIC_STATIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LEPSILON=.TRUE.",_incarpad_) << " # Performing DIELECTRIC_STATIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LRPA=.TRUE.",_incarpad_) << " # Performing DIELECTRIC_STATIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << " # Performing DIELECTRIC_STATIC (Camilo)" << endl;
      // FIX ISMEAR=0 and SIGMA=0.01 ?
      if(0) { // NO NPAR
        if(kflags.KBIN_MPI)
          //      xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(kflags.KBIN_MPI_NCPUS),_incarpad_) << " # Performing DIELECTRIC_STATIC (Camilo)" << endl;
          xvasp.INCAR << aurostd::PaddedPOST("NPAR=1",_incarpad_) << " # Performing DIELECTRIC_STATIC (Camilo)" << endl;
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_STATIC [AFLOW] end" << endl;

      // now to the KPOINTS
      FileContent=xvasp.KPOINTS.str();
      imax=aurostd::GetNLinesString(FileContent);
      aurostd::StringstreamClean(xvasp.KPOINTS);xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA_DELTA - aflow with DK=" << DIELECTRIC_DK << endl;
      xvasp.KPOINTS << aurostd::GetLineString(FileContent,2) << endl;
      KPPRA_DELTA(xvasp.str,DIELECTRIC_DK);  // COUT if you want...
      if(aurostd::FileExist(xvasp.Directory+string("/KPOINTS.bands"))) {  // with some LUCK I`ve already nailed the FCC/HEX
        xvasp.str.kpoints_kscheme=aurostd::GetLineString(FileContent,3);
        string KPOINTS_bands=aurostd::file2string(xvasp.Directory+"/KPOINTS.bands"); // cerr << KPOINTS_bands << endl;
        bool isFCC=aurostd::substring2bool(KPOINTS_bands,"FCC","fcc"); // cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET: isFCC=" << isFCC << endl;
        bool isHEX=aurostd::substring2bool(KPOINTS_bands,"HEX","hex"); // cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET: isHEX=" << isHEX << endl;
        bool isRHL=aurostd::substring2bool(KPOINTS_bands,"RHL","rhl"); // cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET: isRHL=" << isRHL << endl;
        if(abs(xvasp.str.kpoints_s1)+abs(xvasp.str.kpoints_s2)+abs(xvasp.str.kpoints_s3)>0.1) {isFCC=FALSE;isHEX=FALSE;isRHL=FALSE;} // no shift if already shifted
        //  if(isFCC || isHEX) {cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET switching to GAMMA" << endl;}
        if(isFCC || isHEX || isRHL) {
          xvasp.str.kpoints_kscheme="Gamma";
          if(_iseven(xvasp.str.kpoints_k1)) {xvasp.str.kpoints_k1++;xvasp.str.kpoints_s1=0.0;if(LDEBUG) cout << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}
          if(_iseven(xvasp.str.kpoints_k2)) {xvasp.str.kpoints_k2++;xvasp.str.kpoints_s2=0.0;if(LDEBUG) cout << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}
          if(_iseven(xvasp.str.kpoints_k3)) {xvasp.str.kpoints_k3++;xvasp.str.kpoints_s3=0.0;if(LDEBUG) cout << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}
        }  // else dont touch it
      }
      xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
      xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
      xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;
      //    xvasp.KPOINTS << aurostd::GetLineString(FileContent,5) << endl;
    }
    // DYNAMIC
    if(svalue=="DYNAMIC" || svalue=="dynamic") {
      cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET DYNAMIC" << endl;
      // a.  Reuse the STEP 01 WAVECAR
      // b.  ALGO=EXACT NELM=1 LOPTICS=.TRUE.   CSHIFT=0.15  OMEGAMAX=25   NEDOS=12500
      // i.  Remove LEPSILON and LRPA
      string FileContent,strline;
      FileContent=xvasp.INCAR.str();
      aurostd::StringstreamClean(xvasp.INCAR);xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      int imax=aurostd::GetNLinesString(FileContent);
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE) || // Camilo
            //[CO20200624 - OBSOLETE]aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) || // Camilo
            ((aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE)) && (!( aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) || aurostd::substring2bool(strline,"NELMDL",TRUE) || aurostd::substring2bool(strline,"#NELMDL",TRUE) )) ) ||  //CO20200624
            aurostd::substring2bool(strline,"LOPTICS",TRUE) || aurostd::substring2bool(strline,"#LOPTICS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"CSHIFT",TRUE) || aurostd::substring2bool(strline,"#CSHIFT",TRUE) || // Camilo
            aurostd::substring2bool(strline,"OMEGAMAX",TRUE) || aurostd::substring2bool(strline,"#OMEGAMAX",TRUE) || // Camilo
            aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"NBANDS",TRUE) || aurostd::substring2bool(strline,"#NBANDS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LEPSILON",TRUE) || aurostd::substring2bool(strline,"#LEPSILON",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LCALCEPS",TRUE) || aurostd::substring2bool(strline,"#LCALCEPS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LRPA",TRUE) || aurostd::substring2bool(strline,"#LRPA",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT",TRUE)) { // Camilo
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET)" << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // get NBANDS from OUTCAR.dielectric_static
      int NBANDS_OUTCAR=0;
      xOUTCAR OUTCAR_NBANDS(xvasp.Directory+"/OUTCAR.dielectric_static",true);  //quiet, there might be issues with halfway-written OUTCARs
      NBANDS_OUTCAR=OUTCAR_NBANDS.NBANDS;

      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_DYNAMIC [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ALGO=EXACT",_incarpad_) << " # Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELM=1",_incarpad_) << " # Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LOPTICS=.TRUE.",_incarpad_) << " # Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("CSHIFT=0.15",_incarpad_) << " # Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("OMEGAMAX=25",_incarpad_) << " # Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NEDOS=12500",_incarpad_) << " # Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NBANDS="+aurostd::utype2string((double) 4*NBANDS_OUTCAR),_incarpad_) << " # Performing DIELECTRIC_DYNAMIC NBANDS_OUTCAR="  << NBANDS_OUTCAR << " (Camilo)" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_DYNAMIC [AFLOW] end" << endl;
    }
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_REMOVE_ENTRY
namespace KBIN {
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp,const string& keyword,const string& comment,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION  
    vector<string> vkeywords;aurostd::string2tokens(keyword,vkeywords,",");
    return XVASP_INCAR_REMOVE_ENTRY(xvasp,vkeywords,comment,VERBOSE);
  }
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp,const vector<string>& vkeywords,const string& comment,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION 
    //CO20210315 - INCAR always uses keyword=value: https://www.vasp.at/wiki/index.php/INCAR
    //CO20210315 - modifying this function, not using substring2bool() which might match IALGO with ALGO or NELM with NELMIN, keyword is specific key to remove
    //do NOT preload/rewrite the INCAR file, only modify xvasp.INCAR
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::XVASP_INCAR_REMOVE_ENTRY";

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    string patch="#[AFLOW REMOVED"; //("+function+") //CO20210315 - this new comment style avoids figuring out how much to pad previously-padded line
    if(!comment.empty()){patch+="("+comment+")";}
    patch+="] ";

    vector<string> vlines;
    aurostd::string2vectorstring(xvasp.INCAR.str(),vlines);
    aurostd::StringstreamClean(xvasp.INCAR);
    uint iline=0,i=0;
    bool found=false;
    for(iline=0;iline<vlines.size();iline++) {
      const string& strline=vlines[iline];
      found=false;
      for(i=0;i<vkeywords.size()&&!found;i++){
        if(aurostd::kvpair2bool(strline,vkeywords[i],"=")){found=true;}
      }
      if(found){
        if(VERBOSE){xvasp.INCAR << patch << strline << endl;}  //aurostd::PaddedPOST("# "+aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline),_incarpad_) << " # AFLOW REMOVED (" << function << ") " << comment << endl;
      }else{
        if(VERBOSE||strline.length()){xvasp.INCAR << strline << endl;}
      }
    }
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::AFLOWIN_REMOVE and AFLOWIN_ADD
namespace KBIN {
  bool AFLOWIN_REMOVE(const string& aflowin_file,const string& keyword,const string& comment){  //CO20210313
    vector<string> vkeywords;
    aurostd::string2tokens(keyword,vkeywords,",");
    return AFLOWIN_REMOVE(aflowin_file,vkeywords,comment);
  }
  bool AFLOWIN_REMOVE(const string& aflowin_file,const vector<string>& vkeywords,const string& comment){  //CO20210313
    vector<string> vkeywords2ignore;
    return AFLOWIN_REMOVE(aflowin_file,vkeywords,vkeywords2ignore,comment);
  }
  bool AFLOWIN_REMOVE(const string& aflowin_file,const string& keyword,const string& keyword2avoid,const string& comment){  //CO20210313
    vector<string> vkeywords,vkeywords2ignore;
    aurostd::string2tokens(keyword,vkeywords,",");aurostd::string2tokens(keyword2avoid,vkeywords2ignore,",");
    return AFLOWIN_REMOVE(aflowin_file,vkeywords,vkeywords2ignore,comment);
  }
  bool AFLOWIN_REMOVE(const string& aflowin_file,const vector<string>& vkeywords,const vector<string>& vkeywords2ignore,const string& comment){  //CO20210313
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::AFLOWIN_REMOVE";

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vkeywords=" << aurostd::joinWDelimiter(vkeywords,",") << endl;
      cerr << __AFLOW_FUNC__ << " vkeywords2ignore=" << aurostd::joinWDelimiter(vkeywords2ignore,",") << endl;
    }

    if(!aurostd::FileExist(aflowin_file)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,aflowin_file+" not found",_FILE_NOT_FOUND_);}

    vector<string> vlines,vlines_new;
    aurostd::file2vectorstring(aflowin_file,vlines,true,true); //consecutive==true,trim_edges==true (PRESERVE FILE AS MUCH AS POSSIBLE without avoid new newlines at the end)

    uint iline=0,i=0,j=0;

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " FILE(orig)=" << aflowin_file << endl;
      for(iline=0;iline<vlines.size();iline++) {
        cerr << "[iline=" << iline << "]=\"" << vlines[iline] << "\"" << endl;
      }
    }

    string patch="#[AFLOW REMOVED"; //("+function+")
    if(!comment.empty()){patch+=" ("+comment+")";}
    patch+="] ";

    //vkeywords2ignore avoids matching ALGO and IAGO, AFLOW has no particular delimiter (kvpair2bool)

    bool mod_made=false;
    bool found=false,found_ignore=false;
    string line="";
    for(iline=0;iline<vlines.size();iline++) {
      line=aurostd::RemoveComments(vlines[iline]);
      found=false;
      for(i=0;i<vkeywords.size()&&!found;i++){
        if(line.find(vkeywords[i])!=string::npos){
          found_ignore=false;
          for(j=0;j<vkeywords2ignore.size()&&!found_ignore;j++){
            if(line.find(vkeywords2ignore[j])!=string::npos){found_ignore=true;}
          }
          if(found_ignore){continue;}
          else{found=true;}
        }
      }
      if(found==false){vlines_new.push_back(vlines[iline]);}
      else{ //found==true
        vlines_new.push_back( patch+vlines[iline] );
        mod_made=true;
      }
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " FILE(new)=" << aflowin_file << endl;
      for(iline=0;iline<vlines.size();iline++) {
        cerr << "[iline=" << iline << "]=\"" << vlines[iline] << "\"" << endl;
      }
    }

    if(!aurostd::vectorstring2file(vlines_new,aflowin_file)){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"issues writing new "+aflowin_file,_FILE_CORRUPT_);
    }

    return mod_made;
  }
  void AFLOWIN_ADD(const string& aflowin_file,const stringstream& streamin,const string& comment){  //CO20210313
    string lines=streamin.str();
    return AFLOWIN_ADD(aflowin_file,lines,comment);
  }
  void AFLOWIN_ADD(const string& aflowin_file,const ostringstream& streamin,const string& comment){  //CO20210313
    string lines=streamin.str();
    return AFLOWIN_ADD(aflowin_file,lines,comment);
  }
  void AFLOWIN_ADD(const string& aflowin_file,const string& line,const string& comment){  //CO20210313
    vector<string> vlines2add;
    aurostd::string2vectorstring(line,vlines2add,true,true);  //keep lines in the middle, trim edges
    return AFLOWIN_ADD(aflowin_file,vlines2add,comment);
  }
  void AFLOWIN_ADD(const string& aflowin_file,const vector<string>& vlines2add,const string& comment){  //CO20210313
    //always use this function to add data at the end of the _AFLOWIN_
    //echo "" >> _AFLOWIN_ does NOT work in general, might attach data to the end of the previous line (depends on whether there was a final endl)
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    uint iline=0;

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vlines2add=" << endl;
      for(iline=0;iline<vlines2add.size();iline++){cerr << "\"" << vlines2add[iline] << "\"" << endl;}
    }

    if(!aurostd::FileExist(aflowin_file)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,aflowin_file+" not found",_FILE_NOT_FOUND_);}

    vector<string> vlines;
    aurostd::file2vectorstring(aflowin_file,vlines,true,true); //consecutive==true,trim_edges==true (PRESERVE FILE AS MUCH AS POSSIBLE without avoid new newlines at the end)

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " FILE(orig)=" << aflowin_file << endl;
      for(iline=0;iline<vlines.size();iline++) {
        cerr << "[iline=" << iline << "]=\"" << vlines[iline] << "\"" << endl;
      }
    }

    if(!comment.empty()){
      for(iline=0;iline<vlines2add.size();iline++) {
        vlines.push_back( aurostd::PaddedPOST(vlines2add[iline],_AFLOWINPAD_)+" // Self Correction ("+comment+")");
      }
    }
    else{vlines.insert(vlines.end(),vlines2add.begin(),vlines2add.end());}

    if(!aurostd::vectorstring2file(vlines,aflowin_file)){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"issues writing new "+aflowin_file,_FILE_CORRUPT_);
    }
  }
} // namespace KBIN


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// KPOINTS MODIFICATIONS
namespace KBIN {
  void XVASP_KPOINTS_KPOINTS(_xvasp &xvasp,ofstream &FileMESSAGE,bool VERBOSE) {
    //[CO20200624 - OBSOLETE]aurostd::StringstreamClean(xvasp.KPOINTS);
    //[CO20200624 - OBSOLETE]aurostd::StringstreamClean(xvasp.KPOINTS_orig);
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA - aflow with Kpoints=" << xvasp.str.kpoints_kppra << endl;
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;  // MODE AUTO 0
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS.precision(3);
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;
    XVASP_KPOINTS_KPOINTS(xvasp); //CO20200624

    if(VERBOSE) {
      ostringstream aus;
      aurostd::StringstreamClean(aus);
      aus << "00000  MESSAGE KPOINTS K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " ";
      aus << " Kshift=[" << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << "]" << " ";
      aus << " Kpoints=" << xvasp.str.kpoints_kppra  << " with " << xvasp.str.kpoints_kscheme << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
    }
    // KPOINTS done
    //[CO20200624 - OBSOLETE]xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
  }
}

namespace KBIN {
  void XVASP_KPOINTS_KPOINTS(_xvasp &xvasp) {
    aurostd::StringstreamClean(xvasp.KPOINTS);
    aurostd::StringstreamClean(xvasp.KPOINTS_orig);
    //[CO20200624]xvasp.KPOINTS << "KPOINTS File, automatically generated by KBIN::XVASP_KPOINTS_KPOINTS" << endl;
    xvasp.KPOINTS << "KPOINTS File, automatically generated by KBIN::XVASP_KPOINTS_KPOINTS - KPPRA=" << xvasp.str.kpoints_kppra << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
    xvasp.KPOINTS.precision(3);
    xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;
    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
  }
}

namespace KBIN {
  bool XVASP_KPOINTS_IncludesGamma(const _xvasp& xvasp){  //CO20210315
    //https://www.vasp.at/wiki/index.php/KPOINTS#Monkhorst-Pack
    if((!xvasp.str.kpoints_kscheme.empty() && aurostd::toupper(xvasp.str.kpoints_kscheme[0])=='G') &&
        (aurostd::isequal(xvasp.str.kpoints_s1,0.0) && aurostd::isequal(xvasp.str.kpoints_s2,0.0) && aurostd::isequal(xvasp.str.kpoints_s3,0.0))){
      return true; //Gamma with no shift ensures Gamma
    }
    //otherwise, check that we have an odd scheme with no shift, this is also Gamma-centered
    if(!xvasp.str.kpoints_kscheme.empty() && (aurostd::toupper(xvasp.str.kpoints_kscheme[0])=='G' || aurostd::toupper(xvasp.str.kpoints_kscheme[0])=='M') &&
        (_isodd(xvasp.str.kpoints_k1) && _isodd(xvasp.str.kpoints_k2) && _isodd(xvasp.str.kpoints_k3)) &&
        (aurostd::isequal(xvasp.str.kpoints_s1,0.0) && aurostd::isequal(xvasp.str.kpoints_s2,0.0) && aurostd::isequal(xvasp.str.kpoints_s3,0.0))){
      return true;
    }
    return false;
  }
}

namespace KBIN {
  bool XVASP_KPOINTS_OPERATION(_xvasp& xvasp,const string& _operation) {
    //CO20210315 - extensive rewrite
    //the schemes below check if the modification needs to be made (has it already been made?)
    //maintain this feedback system to ensure aflow doesn't keep spinning its wheels on the same fixes
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG); // TRUE;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " operation=" << _operation << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_kscheme=(" << xvasp.str.kpoints_kscheme << ")" << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_k*=(" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ")" << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_s*=(" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << ")" << endl;
    }
    if(_operation.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"empty operation",_INPUT_ILLEGAL_);}  //CO20210315
    string operation=aurostd::toupper(_operation);
    aurostd::StringSubst(operation,"+=1","++"); //backwards compatibility
    aurostd::StringSubst(operation,"-=1","--"); //backwards compatibility
    aurostd::StringSubst(operation,"KMAX","XMAX,YMAX,ZMAX");  //expand
    aurostd::StringSubst(operation,"GAMMA_SHIFT","GAMMA,XZEROSHIFT,YZEROSHIFT,ZZEROSHIFT");  //expand  //repetita iuvant  //can also do XODD,YODD,ZODD
    uint nchanges_made=0; //CO20210315 - prevent aflow from spinning its wheels on an identical input
    //i--: decrease ki by 1 (if possible)
    if(operation.find("X--")!=string::npos){if(xvasp.str.kpoints_k1>=2){xvasp.str.kpoints_k1--;nchanges_made++;if(LDEBUG) cerr << "X-- k1=" << xvasp.str.kpoints_k1 << endl;}}
    if(operation.find("Y--")!=string::npos){if(xvasp.str.kpoints_k2>=2){xvasp.str.kpoints_k2--;nchanges_made++;if(LDEBUG) cerr << "Y-- k2=" << xvasp.str.kpoints_k2 << endl;}}
    if(operation.find("Z--")!=string::npos){if(xvasp.str.kpoints_k3>=2){xvasp.str.kpoints_k3--;nchanges_made++;if(LDEBUG) cerr << "Z-- k3=" << xvasp.str.kpoints_k3 << endl;}}
    //i++: increase ki by 1 (always possible)
    if(operation.find("X++")!=string::npos){xvasp.str.kpoints_k1++;nchanges_made++;if(LDEBUG) cerr << "X++ k1=" << xvasp.str.kpoints_k1 << endl;}
    if(operation.find("Y++")!=string::npos){xvasp.str.kpoints_k2++;nchanges_made++;if(LDEBUG) cerr << "Y++ k2=" << xvasp.str.kpoints_k2 << endl;}
    if(operation.find("Z++")!=string::npos){xvasp.str.kpoints_k3++;nchanges_made++;if(LDEBUG) cerr << "Z++ k3=" << xvasp.str.kpoints_k3 << endl;}
    //i-=2: decrease ki by 2 (if possible)
    //CO20210315 - -=2 and +=2 are important if we want to keep even/odd parity
    //if k1==2, running X-- twice will stop at 1, changing parity
    //these routines fix this
    if(operation.find("X-=2")!=string::npos){if(xvasp.str.kpoints_k1>=3){xvasp.str.kpoints_k1-=2;nchanges_made++;if(LDEBUG) cerr << "X-=2 k1=" << xvasp.str.kpoints_k1 << endl;}} //CO20210315
    if(operation.find("Y-=2")!=string::npos){if(xvasp.str.kpoints_k2>=3){xvasp.str.kpoints_k2-=2;nchanges_made++;if(LDEBUG) cerr << "Y-=2 k2=" << xvasp.str.kpoints_k2 << endl;}} //CO20210315
    if(operation.find("Z-=2")!=string::npos){if(xvasp.str.kpoints_k3>=3){xvasp.str.kpoints_k3-=2;nchanges_made++;if(LDEBUG) cerr << "Z-=2 k3=" << xvasp.str.kpoints_k3 << endl;}} //CO20210315
    //i+=2: increase ki by 2 (always possible)
    if(operation.find("X+=2")!=string::npos){xvasp.str.kpoints_k1+=2;nchanges_made++;if(LDEBUG) cerr << "X+=2 k1=" << xvasp.str.kpoints_k1 << endl;}  //CO20210315
    if(operation.find("Y+=2")!=string::npos){xvasp.str.kpoints_k2+=2;nchanges_made++;if(LDEBUG) cerr << "Y+=2 k2=" << xvasp.str.kpoints_k2 << endl;}  //CO20210315
    if(operation.find("Z+=2")!=string::npos){xvasp.str.kpoints_k3+=2;nchanges_made++;if(LDEBUG) cerr << "Z+=2 k3=" << xvasp.str.kpoints_k3 << endl;}  //CO20210315
    //ieven: make ki even (if not already) AND si=0 (if not already)
    if(operation.find("XEVEN")!=string::npos){
      if(_isodd(xvasp.str.kpoints_k1)){xvasp.str.kpoints_k1++;nchanges_made++;if(LDEBUG) cerr << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}
      if(!aurostd::isequal(xvasp.str.kpoints_s1,0.0)){xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xeven s1=" << xvasp.str.kpoints_s1 << endl;}
    }
    if(operation.find("YEVEN")!=string::npos){
      if(_isodd(xvasp.str.kpoints_k2)){xvasp.str.kpoints_k2++;nchanges_made++;if(LDEBUG) cerr << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}
      if(!aurostd::isequal(xvasp.str.kpoints_s2,0.0)){xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yeven s2=" << xvasp.str.kpoints_s2 << endl;}
    }
    if(operation.find("ZEVEN")!=string::npos){
      if(_isodd(xvasp.str.kpoints_k3)){xvasp.str.kpoints_k3++;nchanges_made++;if(LDEBUG) cerr << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}
      if(!aurostd::isequal(xvasp.str.kpoints_s3,0.0)){xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zeven s3=" << xvasp.str.kpoints_s3 << endl;}
    }
    //iodd: make ki odd (if not already) AND si=0 (if not already)
    if(operation.find("XODD")!=string::npos){
      if(_iseven(xvasp.str.kpoints_k1)){xvasp.str.kpoints_k1++;nchanges_made++;if(LDEBUG) cerr << "Xodd k1=" << xvasp.str.kpoints_k1 << endl;}
      if(!aurostd::isequal(xvasp.str.kpoints_s1,0.0)){xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xodd s1=" << xvasp.str.kpoints_s1 << endl;}
    }
    if(operation.find("YODD")!=string::npos){
      if(_iseven(xvasp.str.kpoints_k2)){xvasp.str.kpoints_k2++;nchanges_made++;if(LDEBUG) cerr << "Yodd k2=" << xvasp.str.kpoints_k2 << endl;}
      if(!aurostd::isequal(xvasp.str.kpoints_s2,0.0)){xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yodd s2=" << xvasp.str.kpoints_s2 << endl;}
    }
    if(operation.find("ZODD")!=string::npos){
      if(_iseven(xvasp.str.kpoints_k3)){xvasp.str.kpoints_k3++;nchanges_made++;if(LDEBUG) cerr << "Zodd k3=" << xvasp.str.kpoints_k3 << endl;}
      if(!aurostd::isequal(xvasp.str.kpoints_s3,0.0)){xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zodd s3=" << xvasp.str.kpoints_s3 << endl;}
    }
    //
    //CO20210315 - below here, do not change _k* (only for MAX)
    xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
    //imax: set ki to kmax (if not already)
    if(operation.find("XMAX")!=string::npos){if(xvasp.str.kpoints_k1!=xvasp.str.kpoints_kmax){xvasp.str.kpoints_k1=xvasp.str.kpoints_kmax;nchanges_made++;if(LDEBUG) cerr << "Xmax k1=" << xvasp.str.kpoints_k1 << endl;}}
    if(operation.find("YMAX")!=string::npos){if(xvasp.str.kpoints_k2!=xvasp.str.kpoints_kmax){xvasp.str.kpoints_k2=xvasp.str.kpoints_kmax;nchanges_made++;if(LDEBUG) cerr << "Ymax k2=" << xvasp.str.kpoints_k2 << endl;}}
    if(operation.find("ZMAX")!=string::npos){if(xvasp.str.kpoints_k3!=xvasp.str.kpoints_kmax){xvasp.str.kpoints_k3=xvasp.str.kpoints_kmax;nchanges_made++;if(LDEBUG) cerr << "Zmax k3=" << xvasp.str.kpoints_k3 << endl;}}
    //ievenshift: even ki get si=0.5 and odd k1 si=0 (if not already)
    if(operation.find("XEVENSHIFT")!=string::npos){
      if(_iseven(xvasp.str.kpoints_k1)){
        if(!aurostd::isequal(xvasp.str.kpoints_s1,0.5)){xvasp.str.kpoints_s1=0.5;nchanges_made++;if(LDEBUG) cerr << "Xevenshift s1=" << xvasp.str.kpoints_s1 << endl;}
      } 
      else{ //_isodd()
        if(!aurostd::isequal(xvasp.str.kpoints_s1,0.0)){xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xevenshift s1=" << xvasp.str.kpoints_s1 << endl;}
      }
    }
    if(operation.find("YEVENSHIFT")!=string::npos){
      if(_iseven(xvasp.str.kpoints_k2)){
        if(!aurostd::isequal(xvasp.str.kpoints_s2,0.5)){xvasp.str.kpoints_s2=0.5;nchanges_made++;if(LDEBUG) cerr << "Yevenshift s2=" << xvasp.str.kpoints_s2 << endl;}
      }
      else{ //_isodd()
        if(!aurostd::isequal(xvasp.str.kpoints_s2,0.0)){xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yevenshift s2=" << xvasp.str.kpoints_s2 << endl;}
      }
    }
    if(operation.find("ZEVENSHIFT")!=string::npos){
      if(_iseven(xvasp.str.kpoints_k3)){
        if(!aurostd::isequal(xvasp.str.kpoints_s3,0.5)){xvasp.str.kpoints_s3=0.5;nchanges_made++;if(LDEBUG) cerr << "Zevenshift s3=" << xvasp.str.kpoints_s3 << endl;}
      } 
      else{ //_isodd()
        if(!aurostd::isequal(xvasp.str.kpoints_s3,0.0)){xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zevenshift s3=" << xvasp.str.kpoints_s3 << endl;}
      }
    }
    //ioddshift: odd ki get si=0.5 and even k1 si=0 (if not already)
    if(operation.find("XODDSHIFT")!=string::npos){
      if(_isodd(xvasp.str.kpoints_k1)){
        if(!aurostd::isequal(xvasp.str.kpoints_s1,0.5)){xvasp.str.kpoints_s1=0.5;nchanges_made++;if(LDEBUG) cerr << "Xoddshift s1=" << xvasp.str.kpoints_s1 << endl;}
      } 
      else{ //_iseven()
        if(!aurostd::isequal(xvasp.str.kpoints_s1,0.0)){xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xoddshift s1=" << xvasp.str.kpoints_s1 << endl;}
      }
    }
    if(operation.find("YODDSHIFT")!=string::npos){
      if(_isodd(xvasp.str.kpoints_k2)){
        if(!aurostd::isequal(xvasp.str.kpoints_s2,0.5)){xvasp.str.kpoints_s2=0.5;nchanges_made++;if(LDEBUG) cerr << "Yoddshift s2=" << xvasp.str.kpoints_s2 << endl;}
      }
      else{ //_iseven()
        if(!aurostd::isequal(xvasp.str.kpoints_s2,0.0)){xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yoddshift s2=" << xvasp.str.kpoints_s2 << endl;}
      }
    }
    if(operation.find("ZODDSHIFT")!=string::npos){
      if(_isodd(xvasp.str.kpoints_k3)){
        if(!aurostd::isequal(xvasp.str.kpoints_s3,0.5)){xvasp.str.kpoints_s3=0.5;nchanges_made++;if(LDEBUG) cerr << "Zoddshift s3=" << xvasp.str.kpoints_s3 << endl;}
      } 
      else{ //_iseven()
        if(!aurostd::isequal(xvasp.str.kpoints_s3,0.0)){xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zoddshift s3=" << xvasp.str.kpoints_s3 << endl;}
      }
    }
    //izeroshift: si=0 (if not already)
    if(operation.find("XZEROSHIFT")!=string::npos){if(!aurostd::isequal(xvasp.str.kpoints_s1,0.0)){xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xzeroshift s1=" << xvasp.str.kpoints_s1 << endl;}}
    if(operation.find("YZEROSHIFT")!=string::npos){if(!aurostd::isequal(xvasp.str.kpoints_s2,0.0)){xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yzeroshift s2=" << xvasp.str.kpoints_s2 << endl;}}
    if(operation.find("ZZEROSHIFT")!=string::npos){if(!aurostd::isequal(xvasp.str.kpoints_s3,0.0)){xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zzeroshift s3=" << xvasp.str.kpoints_s3 << endl;}}
    //M vs. G: set kscheme accordingly (if not already)
    if(operation.at(0)=='M' || operation.find("MONKHORST-PACK")!=string::npos){if(xvasp.str.kpoints_kscheme!="Monkhorst-Pack"){xvasp.str.kpoints_kscheme="Monkhorst-Pack";nchanges_made++;if(LDEBUG) cerr << "Monkhorst-Pack" << endl;}} // "Monkhorst-Pack";
    if(operation.at(0)=='G' || operation.find("GAMMA")!=string::npos){if(xvasp.str.kpoints_kscheme!="Gamma"){xvasp.str.kpoints_kscheme="Gamma";nchanges_made++;if(LDEBUG) cerr << "Gamma" << endl;}} // "Gamma";  //CO20210315 - do not check if !KBIN::XVASP_KPOINTS_IncludesGamma(xvasp) here, changing to 'G' if 'G' is not set is an entirely new scheme: https://www.vasp.at/wiki/index.php/KPOINTS#Monkhorst-Pack
    //
    //[CO20210315 - better not to do something not requested, just return false]//CO20210315 - below here, we attempt to make a change if a change was not effected earlier
    //[CO20210315 - better not to do something not requested, just return false]if((nchanges_made==0)&&(operation.find("XEVEN")!=string::npos||operation.find("YEVEN")!=string::npos||operation.find("ZEVEN")!=string::npos)){  //CO20210201 - was already even before, increment by 2
    //[CO20210315 - better not to do something not requested, just return false]  if(operation.find("XEVEN")!=string::npos) {xvasp.str.kpoints_k1+=2;xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}
    //[CO20210315 - better not to do something not requested, just return false]  if(operation.find("YEVEN")!=string::npos) {xvasp.str.kpoints_k2+=2;xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}
    //[CO20210315 - better not to do something not requested, just return false]  if(operation.find("ZEVEN")!=string::npos) {xvasp.str.kpoints_k3+=2;xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}
    //[CO20210315 - better not to do something not requested, just return false]}
    //[CO20210315 - better not to do something not requested, just return false]if((nchanges_made==0)&&(operation.find("XZERO")!=string::npos||operation.find("YZERO")!=string::npos||operation.find("ZZERO")!=string::npos)){  //CO20210201 - was already odd before, increment by 2
    //[CO20210315 - better not to do something not requested, just return false]  if(operation.find("XZERO")!=string::npos) {xvasp.str.kpoints_k1+=2;xvasp.str.kpoints_s1=0.0;nchanges_made++;if(LDEBUG) cerr << "Xodd k1=" << xvasp.str.kpoints_k1 << endl;}
    //[CO20210315 - better not to do something not requested, just return false]  if(operation.find("YZERO")!=string::npos) {xvasp.str.kpoints_k2+=2;xvasp.str.kpoints_s2=0.0;nchanges_made++;if(LDEBUG) cerr << "Yodd k2=" << xvasp.str.kpoints_k2 << endl;}
    //[CO20210315 - better not to do something not requested, just return false]  if(operation.find("ZZERO")!=string::npos) {xvasp.str.kpoints_k3+=2;xvasp.str.kpoints_s3=0.0;nchanges_made++;if(LDEBUG) cerr << "Zodd k3=" << xvasp.str.kpoints_k3 << endl;}
    //[CO20210315 - better not to do something not requested, just return false]}
    //above here, modify change_made
    if(xvasp.str.kpoints_k1<1) xvasp.str.kpoints_k1=1;
    if(xvasp.str.kpoints_k2<1) xvasp.str.kpoints_k2=1;
    if(xvasp.str.kpoints_k3<1) xvasp.str.kpoints_k3=1;
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " nchanges_made=" << nchanges_made << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_kscheme=(" << xvasp.str.kpoints_kscheme << ")" << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_k*=(" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ")" << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_s*=(" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << ")" << endl;
    }
    xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
    xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
    KBIN::XVASP_KPOINTS_KPOINTS(xvasp);
    return (nchanges_made>0);
  }
}

// bool XVASP_KPOINTS_Kscheme(_xvasp& xvasp,string kscheme) {
//   xvasp.str.kpoints_kscheme=kscheme;
//   if(kscheme.at(0)=='M' || kscheme.at(0)=='m') xvasp.str.kpoints_kscheme="Monkhorst-Pack";// "Monkhorst-Pack";
//   if(kscheme.at(0)=='G' || kscheme.at(0)=='g') xvasp.str.kpoints_kscheme="Gamma";// "Gamma";
//   KBIN::XVASP_KPOINTS_KPOINTS(xvasp);
//   cerr << XPID << "KBIN::XVASP_KPOINTS_Kscheme" << endl << xvasp.KPOINTS.str() << endl;
//   return TRUE;
// }


//[CO20210315 - OBSOLETE]namespace KBIN {
//[CO20210315 - OBSOLETE]  bool XVASP_KPOINTS_Fix_KPPRA(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE) {
//[CO20210315 - OBSOLETE]    ostringstream aus;
//[CO20210315 - OBSOLETE]    string stringKPPRA;
//[CO20210315 - OBSOLETE]    stringKPPRA=KPPRA(xvasp.str,NK);
//[CO20210315 - OBSOLETE]    if(VERBOSE) FileMESSAGE << stringKPPRA;
//[CO20210315 - OBSOLETE]    // VERBOSE on the SCREEN
//[CO20210315 - OBSOLETE]    if(VERBOSE) {
//[CO20210315 - OBSOLETE]      aus.precision(5);
//[CO20210315 - OBSOLETE]      aus << "00000  MESSAGE KPOINTS KPPRA routine ["
//[CO20210315 - OBSOLETE]        <<  xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << "]=" << xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3 << "=["
//[CO20210315 - OBSOLETE]        <<  modulus(xvasp.str.klattice(1)/((double) xvasp.str.kpoints_k1)) << "," << modulus(xvasp.str.klattice(2)/((double) xvasp.str.kpoints_k2)) << ","
//[CO20210315 - OBSOLETE]        <<  modulus(xvasp.str.klattice(3)/((double) xvasp.str.kpoints_k3)) << "]" << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    return TRUE;
//[CO20210315 - OBSOLETE]  }
//[CO20210315 - OBSOLETE]}

//[CO20210315 - OBSOLETE]namespace KBIN {
//[CO20210315 - OBSOLETE]  bool XVASP_KPOINTS_Fix_KSHIFT(_xvasp &xvasp,_xvasp &rxvasp,bool KAUTOSHIFT,bool VERBOSE) {
//[CO20210315 - OBSOLETE]    stringstream aus;
//[CO20210315 - OBSOLETE]    if(KAUTOSHIFT) {
//[CO20210315 - OBSOLETE]      if(_isodd(xvasp.str.kpoints_k1+rxvasp.str.kpoints_k1)) { xvasp.str.kpoints_s1=0.5; } else { xvasp.str.kpoints_s1=0.0; } 
//[CO20210315 - OBSOLETE]      if(_isodd(xvasp.str.kpoints_k2+rxvasp.str.kpoints_k2)) { xvasp.str.kpoints_s2=0.5; } else { xvasp.str.kpoints_s2=0.0; } 
//[CO20210315 - OBSOLETE]      if(_isodd(xvasp.str.kpoints_k3+rxvasp.str.kpoints_k3)) { xvasp.str.kpoints_s3=0.5; } else { xvasp.str.kpoints_s3=0.0; } 
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    if(VERBOSE) {
//[CO20210315 - OBSOLETE]      aus << "00000  MESSAGE KPOINTS KSHIFT=[" << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << "]" << endl;
//[CO20210315 - OBSOLETE]      // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
//[CO20210315 - OBSOLETE]      cout << aus.str();
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    return TRUE;
//[CO20210315 - OBSOLETE]  }
//[CO20210315 - OBSOLETE]}

//[CO20210315 - OBSOLETE]namespace KBIN {
//[CO20210315 - OBSOLETE]  bool XVASP_KPOINTS_Fix_KPOINTS(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE) {
//[CO20210315 - OBSOLETE]    stringstream aus;
//[CO20210315 - OBSOLETE]    // create KPOINTS stringstream from KPPRA
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(xvasp.KPOINTS);
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(xvasp.KPOINTS_orig);
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA - aflow with Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK << "]" << endl;
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;  // MODE AUTO 0
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS.precision(3);
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;
//[CO20210315 - OBSOLETE]
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(aus);
//[CO20210315 - OBSOLETE]    if(VERBOSE) {
//[CO20210315 - OBSOLETE]      aus << "00000  MESSAGE KPOINTS K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << endl;
//[CO20210315 - OBSOLETE]      aus << "00000  MESSAGE KPOINTS Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK  << "]    with " << xvasp.str.kpoints_kscheme << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // KPOINTS done
//[CO20210315 - OBSOLETE]    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
//[CO20210315 - OBSOLETE]    return TRUE;
//[CO20210315 - OBSOLETE]  }
//[CO20210315 - OBSOLETE]}

namespace KBIN {
  bool XVASP_KPOINTS_isAutoMesh(const _xvasp& xvasp){ //CO20210315
    //returns if KPOINTS is automatic generation scheme: https://www.vasp.at/wiki/index.php/KPOINTS
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_IVASP_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.KPOINTS=" << endl;cerr << xvasp.KPOINTS.str() << endl;}

    vector<string> vlines;
    aurostd::string2vectorstring(aurostd::RemoveComments(xvasp.KPOINTS.str()),vlines);
    if(vlines.size()<2){return false;}
    string nkpoints=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vlines[1]);
    if(aurostd::isfloat(nkpoints) && aurostd::string2utype<uint>(nkpoints)==0){return true;} //is an auto-mesh

    return false;
  }
} // namespace KBIN

namespace KBIN {
  bool XVASP_KPOINTS_string2numbers(_xvasp& xvasp) {  //CO20210315 - cleaned up
    //CO20210315 - can only read auto-meshes of MP or G
    //do not patch for generic KPOINTS files, this will break Afix which might change KPOINTS when it needs to keep it fixed for BANDS runs
    //https://www.vasp.at/wiki/index.php/KPOINTS
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.KPOINTS=" << endl;cerr << xvasp.KPOINTS.str() << endl;}

    vector<string> vlines;
    aurostd::string2vectorstring(aurostd::RemoveComments(xvasp.KPOINTS.str()),vlines);
    if(vlines.size()<3){return false;}
    string nkpoints=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vlines[1]);
    if(!(aurostd::isfloat(nkpoints) && aurostd::string2utype<uint>(nkpoints)==0)){return false;} //not an auto-mesh
    string scheme=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vlines[2]);
    if(scheme.empty()||(aurostd::toupper(scheme[0])!='M' && aurostd::toupper(scheme[0])!='G')){return false;}

    xvasp.str.kpoints_kscheme=scheme;

    vector<string> vtokens;
    if(vlines.size()>3){  //kgrid
      aurostd::string2tokens(vlines[3],vtokens," ");
      if(vtokens.size()==3){
        xvasp.str.kpoints_k1=aurostd::string2utype<int>(vtokens[0]);
        xvasp.str.kpoints_k2=aurostd::string2utype<int>(vtokens[1]);
        xvasp.str.kpoints_k3=aurostd::string2utype<int>(vtokens[2]);
      }
    }
    if(vlines.size()>4){  //kshift
      aurostd::string2tokens(vlines[4],vtokens," ");
      if(vtokens.size()==3){
        xvasp.str.kpoints_s1=aurostd::string2utype<int>(vtokens[0]);
        xvasp.str.kpoints_s2=aurostd::string2utype<int>(vtokens[1]);
        xvasp.str.kpoints_s3=aurostd::string2utype<int>(vtokens[2]);
      }
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_kscheme=" << xvasp.str.kpoints_kscheme << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_k*=(" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ")" << endl;
      cerr << __AFLOW_FUNC__ << " xvasp.str.kpoints_s*=(" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << ")" << endl;
    }

    return true;
  }
}
namespace KBIN {
  bool XVASP_KPOINTS_numbers2string(_xvasp& xvasp) {  //CO20210315 - cleaned up
    //CO20210315 - can only read auto-meshes of MP or G
    //https://www.vasp.at/wiki/index.php/KPOINTS
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.KPOINTS(pre)=" << endl;cerr << xvasp.KPOINTS.str() << endl;}

    if(xvasp.str.kpoints_kscheme.empty()||(aurostd::toupper(xvasp.str.kpoints_kscheme[0])!='M' && aurostd::toupper(xvasp.str.kpoints_kscheme[0])!='G')){return false;}

    vector<string> vlines;
    aurostd::string2vectorstring(xvasp.KPOINTS.str(),vlines);
    if(vlines.size()==0){return false;}

    aurostd::StringstreamClean(xvasp.KPOINTS);

    xvasp.KPOINTS << vlines[0] << endl;
    xvasp.KPOINTS << 0 << endl; //auto generation (M or G)  //vlines[1]
    xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl; //CO20210315 - adding shift

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xvasp.KPOINTS(post)=" << endl;cerr << xvasp.KPOINTS.str() << endl;}
    return true;
  }
}
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// FIX MODIFICATIONS

namespace KBIN {
  void XVASP_Afix_Clean(const _xvasp& xvasp,const string& preserve_name) {
    aurostd::file2file(xvasp.Directory+"/"+DEFAULT_VASP_OUT,xvasp.Directory+"/"+preserve_name);
    deque<string> vrem;aurostd::string2tokens("aflow.tmp,CHG,CONTCAR,DOSCAR,EIGENVAL,IBZKPT,OUTCAR,OSZICAR,PCDAT,WAVECAR,aflow.qsub*,XDATCAR,vasprun.xml,"+DEFAULT_VASP_OUT+",core*",vrem,",");
    for (uint i=0;i<vrem.size();i++)
      aurostd::RemoveFile(xvasp.Directory+"/"+vrem.at(i));
    if(!xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED")) aurostd::RemoveFile(xvasp.Directory+"/CHGCAR"); 
  }
}
//[CO20210315 - OBSOLETE]namespace KBIN {
//[CO20210315 - OBSOLETE]  bool XVASP_Afix_ROTMAT(_xvasp& xvasp,int mode,_aflags& aflags,bool VERBOSE,ofstream &FileMESSAGE) {  //CO20200624
//[CO20210315 - OBSOLETE]    _kflags kflags;_vflags vflags; // temporary
//[CO20210315 - OBSOLETE]    return XVASP_Afix_ROTMAT(xvasp,mode,kflags,vflags,aflags,VERBOSE,FileMESSAGE);
//[CO20210315 - OBSOLETE]  }
//[CO20210315 - OBSOLETE]  bool XVASP_Afix_ROTMAT(_xvasp& xvasp,int mode,_kflags kflags,_vflags vflags,_aflags& aflags,bool VERBOSE,ofstream &FileMESSAGE) {
//[CO20210315 - OBSOLETE]    // https://www.vasp.at/wiki/index.php/KPOINTS
//[CO20210315 - OBSOLETE]    // cerr << XPID << "KBIN::XVASP_Afix_ROTMAT mode=" << mode << endl;
//[CO20210315 - OBSOLETE]    // mode=0 put gamma and even
//[CO20210315 - OBSOLETE]    // mode=1 put gamma and odd
//[CO20210315 - OBSOLETE]    // mode=2 put kmax
//[CO20210315 - OBSOLETE]    // mode=3 put odd
//[CO20210315 - OBSOLETE]    // mode=4 put even
//[CO20210315 - OBSOLETE]    // mode=5 SYMPREC
//[CO20210315 - OBSOLETE]    // mode=6 ISYM=0
//[CO20210315 - OBSOLETE]    // AFIX *************************************
//[CO20210315 - OBSOLETE]    if(!kflags.AFLOW_MODE_VASP){;}  //keep busy
//[CO20210315 - OBSOLETE]    int kmax=xvasp.str.kpoints_kmax;
//[CO20210315 - OBSOLETE]    // int k1=xvasp.str.kpoints_k1,k2=xvasp.str.kpoints_k2,k3=xvasp.str.kpoints_k3,kppra=xvasp.str.kpoints_kppra;
//[CO20210315 - OBSOLETE]    // double s1=xvasp.str.kpoints_s1,s2=xvasp.str.kpoints_s2,s3=xvasp.str.kpoints_s3;
//[CO20210315 - OBSOLETE]    bool LQUIET=!VERBOSE;
//[CO20210315 - OBSOLETE]    ostringstream aus;
//[CO20210315 - OBSOLETE]    bool DEBUG_WRITE_OUT_KPOINTS_FILE=FALSE;
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]if(mode<=0) {	
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]  KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (XXXX) ");
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]  aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (XXXX) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]  aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]}
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==0) { // put gamma and even
//[CO20210315 - OBSOLETE]      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return false; // don't touch kpoints
//[CO20210315 - OBSOLETE]      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xevenshift,Yevenshift,Zevenshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp);
//[CO20210315 - OBSOLETE]      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
//[CO20210315 - OBSOLETE]      if(DEBUG_WRITE_OUT_KPOINTS_FILE) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_Kshift_Gamma_EVEN"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_KSHIFT_EVEN) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_KSHIFT_EVEN) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==1) { // put gamma and odd
//[CO20210315 - OBSOLETE]      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return false; // don't touch kpoints
//[CO20210315 - OBSOLETE]      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xoddshift,Yoddshift,Zoddshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_ODD(xvasp);
//[CO20210315 - OBSOLETE]      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
//[CO20210315 - OBSOLETE]      if(DEBUG_WRITE_OUT_KPOINTS_FILE) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_Kshift_Gamma_ODD"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_KSHIFT_ODD) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_KSHIFT_ODD) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==2) { // put kmax
//[CO20210315 - OBSOLETE]      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return false; // don't touch kpoints
//[CO20210315 - OBSOLETE]      if(xvasp.str.atoms.size()>50) return false;  //CO20200624 - skipping kmax for big cells
//[CO20210315 - OBSOLETE]      xvasp.str.kpoints_k1=kmax;xvasp.str.kpoints_k2=kmax;xvasp.str.kpoints_k3=kmax;
//[CO20210315 - OBSOLETE]      xvasp.str.kpoints_kppra=kmax*kmax*kmax*xvasp.str.atoms.size();
//[CO20210315 - OBSOLETE]      xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;
//[CO20210315 - OBSOLETE]      // if(_isodd(xvasp.str.kpoints_k1))  {xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;}
//[CO20210315 - OBSOLETE]      // if(_iseven(xvasp.str.kpoints_k1)) {xvasp.str.kpoints_s1=0.5;xvasp.str.kpoints_s2=0.5;xvasp.str.kpoints_s3=0.5;}
//[CO20210315 - OBSOLETE]      KBIN::XVASP_KPOINTS_KPOINTS(xvasp);
//[CO20210315 - OBSOLETE]      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
//[CO20210315 - OBSOLETE]      if(DEBUG_WRITE_OUT_KPOINTS_FILE) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_KPOINTS_MAX"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_MAX) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_MAX) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==3) { // put odd
//[CO20210315 - OBSOLETE]      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return false; // don't touch kpoints
//[CO20210315 - OBSOLETE]      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);
//[CO20210315 - OBSOLETE]      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
//[CO20210315 - OBSOLETE]      if(DEBUG_WRITE_OUT_KPOINTS_FILE) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_ODD"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_ODD) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_ODD) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==4) { // put even
//[CO20210315 - OBSOLETE]      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return false; // don't touch kpoints
//[CO20210315 - OBSOLETE]      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xeven,Yeven,Zeven"); //    KBIN::XVASP_KPOINTS_EVEN(xvasp);
//[CO20210315 - OBSOLETE]      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
//[CO20210315 - OBSOLETE]      if(DEBUG_WRITE_OUT_KPOINTS_FILE) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_EVEN"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_EVEN) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_EVEN) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==5) { // SYMPREC
//[CO20210315 - OBSOLETE]      KBIN::XVASP_INCAR_PREPARE_GENERIC("SYMPREC",xvasp,vflags,"",0,0.0,OFF);
//[CO20210315 - OBSOLETE]      // [OBSOLETE]      KBIN::XVASP_INCAR_SYM(xvasp,OFF,TRUE);
//[CO20210315 - OBSOLETE]      //[CO20200624 - OBSOLETE]aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (SYM=OFF) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (SYMPREC) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // mode *************************************
//[CO20210315 - OBSOLETE]    if(mode==6) { // SYM is the desperate case
//[CO20210315 - OBSOLETE]      KBIN::XVASP_INCAR_PREPARE_GENERIC("SYM",xvasp,vflags,"",0,0.0,OFF);
//[CO20210315 - OBSOLETE]      // [OBSOLETE]      KBIN::XVASP_INCAR_SYM(xvasp,OFF,TRUE);
//[CO20210315 - OBSOLETE]      //[CO20200624 - OBSOLETE]aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
//[CO20210315 - OBSOLETE]      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (SYM=OFF) ");
//[CO20210315 - OBSOLETE]      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (SYM=OFF) - " << Message(__AFLOW_FILE__,aflags) << endl;
//[CO20210315 - OBSOLETE]      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
//[CO20210315 - OBSOLETE]    }
//[CO20210315 - OBSOLETE]    // reload to restart ---------------------------------
//[CO20210315 - OBSOLETE]
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(xvasp.KPOINTS_orig); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(xvasp.KPOINTS); xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS");
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(xvasp.INCAR_orig); xvasp.INCAR_orig << xvasp.INCAR.str();
//[CO20210315 - OBSOLETE]    aurostd::StringstreamClean(xvasp.INCAR); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR");
//[CO20210315 - OBSOLETE]    // clean to restart ----------------------------------    
//[CO20210315 - OBSOLETE]    XVASP_Afix_Clean(xvasp,"aflow.error.rotmat"); //CO20200624
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]aurostd::execute("cat "+xvasp.Directory+"/"+DEFAULT_VASP_OUT+" >> "+xvasp.Directory+"/aflow.error.rotmat");
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]deque<string> vrem;aurostd::string2tokens("aflow.tmp,CHG,CONTCAR,DOSCAR,EIGENVAL,IBZKPT,OUTCAR,OSZICAR,PCDAT,WAVECAR,aflow.qsub*,XDATCAR,vasprun.xml,"+DEFAULT_VASP_OUT+",core*",vrem,",");
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]for (uint i=0;i<vrem.size();i++)
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]  aurostd::RemoveFile(xvasp.Directory+"/"+vrem.at(i));
//[CO20210315 - OBSOLETE]    //[CO20200624 - OBSOLETE]if(!xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED")) aurostd::RemoveFile(xvasp.Directory+"/CHGCAR");
//[CO20210315 - OBSOLETE]    return true;
//[CO20210315 - OBSOLETE]  }
//[CO20210315 - OBSOLETE]}

namespace KBIN {
  bool XVASP_Afix_NBANDS(_xvasp& xvasp,const _aflags& aflags,const _vflags& vflags,bool increase) {
    int nbands=0;
    return XVASP_Afix_NBANDS(xvasp,aflags,vflags,nbands,increase);
  }
  bool XVASP_Afix_NBANDS(_xvasp& xvasp,const _aflags& aflags,const _vflags& vflags,int& nbands,bool increase) {
    //note, this function is part of the XVASP_Afix_* series
    //it will assume xvasp.INCAR has been pre-loaded and will NOT rewrite the INCAR
    //this is all done inside the main XVASP_Afix() function
    //BE CAREFUL not to overwrite xvasp.INCAR
    //for now the default is to increase NBANDS, we might decrease later as a fix to MEMORY
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::XVASP_Afix_NBANDS";

    // get NBANDS from OUTCAR
    aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    nbands=0;
    int nelectrons=0;
    if(aurostd::FileExist(xvasp.Directory+"/OUTCAR")&&aurostd::FileNotEmpty(xvasp.Directory+"/OUTCAR")){
      xOUTCAR OUTCAR_NBANDS(xvasp.Directory+"/OUTCAR",true); //quiet, there might be issues with halfway-written OUTCARs
      nbands=OUTCAR_NBANDS.NBANDS;
      nelectrons=OUTCAR_NBANDS.nelectrons;
    }
    if(nbands==0 && aurostd::kvpair2bool(xvasp.INCAR,"NBANDS","=")) {  // GET BANDS FROM INCAR //CO20210315 - the kvpair2bool() is more robust than substring2bool
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " GET NBANDS FROM INCAR" << endl;
      nbands=aurostd::kvpair2utype<int>(xvasp.INCAR,"NBANDS","=");  //CO20210315
      //no need to spit error, if it doesn't find NBANDS in INCAR, then use defaults (below)
    }
    //CO20210315 - get NPAR, as NBANDS will be set internally to something divisible by NPAR
    int npar=1;  //default in vasp in ncore, which is 1
    if(aurostd::kvpair2bool(xvasp.INCAR,"NPAR","=")) {
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " GET NPAR FROM INCAR" << endl;
      npar=aurostd::kvpair2utype<int>(xvasp.INCAR,"NPAR","=");  //CO20210315
    }
    else if(aurostd::kvpair2bool(xvasp.INCAR,"NCORE","=")) {  //default is npar=ncore
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " GET NCORE FROM INCAR" << endl;
      npar=aurostd::kvpair2utype<int>(xvasp.INCAR,"NCORE","=");  //CO20210315
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " nbands=" << nbands << endl;
      cerr << __AFLOW_FUNC__ << " nelectrons=" << nelectrons << endl;
      cerr << __AFLOW_FUNC__ << " npar=" << npar << endl;
    }
    if(nbands==0){nbands=KBIN::XVASP_INCAR_GetNBANDS(xvasp,aflags,TRUE);}
    else{
      //increment choices are arbitrary, the increase is historical
      //decrease reduces by 20% every time
      //it is good to reduce more slowly than we increase, decreasing NBANDS can lead to issues with the calculation
      //increasing NBANDS can only lead to memory issues
      if(increase){nbands+=(int)(20+(double)nbands*0.2);}
      else{nbands=(int)((double)nbands*0.8);}  //CO20210315 //CO20210713 - increasing rate of decrease for NBANDS to 20%, 10% is too slow, issues with walltime
    }
    //make sure nbands is divisible by npar (default =1)
    int increment=(increase?+1:-1);
    while((nbands%npar)!=0){nbands+=increment;}

    if(LDEBUG) cerr << __AFLOW_FUNC__ << " nbands=" << nbands << endl;

    if(nelectrons!=0 && nbands<(nelectrons+1)/2){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " nbands < " << (nelectrons+1)/2 << " = (nelectrons+1)/2 (too small)" << endl;}
      return false;
    }

    if(!KBIN::XVASP_INCAR_PREPARE_GENERIC("NBANDS",xvasp,vflags,"",nbands,0.0,FALSE)){return false;}

    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //CO20200624

    // fix aflow.in
    if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
      KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]NBANDS",function); //CO20210315 - no equals here, remove ALL NBANDS
      KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]NBANDS="+aurostd::utype2string(nbands,_IVASP_DOUBLE2STRING_PRECISION_),function);
    }

    return true;
  }
}

namespace KBIN {
  bool XVASP_Afix_POTIM(_xvasp& xvasp,const _vflags& vflags) {
    double potim=0.0;
    return XVASP_Afix_POTIM(xvasp,vflags,potim);
  }
  bool XVASP_Afix_POTIM(_xvasp& xvasp,const _vflags& vflags,double& potim) {
    //note, this function is part of the XVASP_Afix_* series
    //it will assume xvasp.INCAR has been pre-loaded and will NOT rewrite the INCAR
    //this is all done inside the main XVASP_Afix() function
    //BE CAREFUL not to overwrite xvasp.INCAR
    //there is no need for increasing POTIM, as increasing POTIM simply speeds up the calculation
    //you can ignore "BRIONS problems: POTIM should be increased", not a problem
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::XVASP_Afix_POTIM";

    // get POTIM from OUTCAR
    aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    potim=0.0;
    if(aurostd::isequal(potim,0.0) && aurostd::FileExist(xvasp.Directory+"/OUTCAR")&&aurostd::FileNotEmpty(xvasp.Directory+"/OUTCAR")){
      xOUTCAR OUTCAR_POTIM(xvasp.Directory+"/OUTCAR",true);  //quiet, there might be issues with halfway-written OUTCARs
      potim=OUTCAR_POTIM.POTIM;
    }
    if(aurostd::isequal(potim,0.0) && aurostd::kvpair2bool(xvasp.INCAR,"POTIM","=")) {  // GET POTIM FROM INCAR  //CO20210315 - the kvpair2bool() is more robust than substring2bool
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " GET POTIM FROM INCAR" << endl;
      potim=aurostd::kvpair2utype<double>(xvasp.INCAR,"POTIM","="); //CO20210315
    }
    if(aurostd::isequal(potim,0.0)){potim=0.5;} //default - https://www.vasp.at/wiki/index.php/POTIM
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " potim=" << potim << endl;
    potim/=2.0;
    if(potim<0.01){return false;} //too low, not worth it
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " potim=" << potim << endl;

    if(!KBIN::XVASP_INCAR_PREPARE_GENERIC("POTIM",xvasp,vflags,"",0,potim,FALSE)){return false;}

    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //CO20200624

    // fix aflow.in
    if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
      KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]POTIM=",function);
      KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]POTIM="+aurostd::utype2string(potim,_IVASP_DOUBLE2STRING_PRECISION_),function);
    }

    return true;
  }
}

namespace KBIN {
  bool XVASP_Afix_NELM(_xvasp& xvasp,const _vflags& vflags) {  //CO20200624
    int nelm=0;
    return XVASP_Afix_NELM(xvasp,vflags,nelm);
  }
  bool XVASP_Afix_NELM(_xvasp& xvasp,const _vflags& vflags,int& nelm) {  //CO20200624
    //note, this function is part of the XVASP_Afix_* series
    //it will assume xvasp.INCAR has been pre-loaded and will NOT rewrite the INCAR
    //this is all done inside the main XVASP_Afix() function
    //BE CAREFUL not to overwrite xvasp.INCAR
    //there is no need to decrease NELM, as it will stop when converged no matter what
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::XVASP_Afix_NELM";

    // get NELM from OUTCAR
    aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    nelm=0;
    if(aurostd::FileExist(xvasp.Directory+"/OUTCAR")&&aurostd::FileNotEmpty(xvasp.Directory+"/OUTCAR")){
      xOUTCAR OUTCAR_NELM(xvasp.Directory+"/OUTCAR",true); //quiet, there might be issues with halfway-written OUTCARs
      nelm=OUTCAR_NELM.NELM;
    }
    if(nelm==0 && aurostd::kvpair2bool(xvasp.INCAR,"NELM","=")) {  // GET NELM FROM INCAR  //CO20210315 - the kvpair2bool() is more robust than substring2bool
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " GET NELM FROM INCAR" << endl;
      nelm=aurostd::kvpair2utype<int>(xvasp.INCAR,"NELM","=");  //CO20210315
    }
    if(nelm==0){nelm=60;} //default - https://www.vasp.at/wiki/index.php/NELM
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " nelm=" << nelm << endl;
    //[CO20200624 - no point in setting over and over again, just do in one shot]nelm=nelm*2.0;
    if(nelm>=MAX_VASP_NELM){return false;}  //too high as set by user
    nelm=MAX_VASP_NELM;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " nelm=" << nelm << endl;

    if(!KBIN::XVASP_INCAR_PREPARE_GENERIC("NELM",xvasp,vflags,"",nelm,0.0,FALSE)){return false;}

    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //CO20200624
    //[CO20200624 - do NOT fix aflow.in, keep fix LOCAL]// fix aflow.in

    return true;
  }
}

// ***************************************************************************
// KBIN::XVASP_Afix_EFIELD_PEAD
namespace KBIN {
  bool XVASP_Afix_EFIELD_PEAD(_xvasp& xvasp,const _vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    xvector<double> E;
    return XVASP_Afix_EFIELD_PEAD(xvasp,vflags,E);
  }
  bool XVASP_Afix_EFIELD_PEAD(_xvasp& xvasp,const _vflags& vflags,xvector<double>& E) {        // AFLOW_FUNCTION_IMPLEMENTATION
    //note, this function is part of the XVASP_Afix_* series
    //it will assume xvasp.INCAR has been pre-loaded and will NOT rewrite the INCAR
    //this is all done inside the main XVASP_Afix() function
    //BE CAREFUL not to overwrite xvasp.INCAR
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(aurostd::FileExist(xvasp.Directory+"/OUTCAR")&&aurostd::FileNotEmpty(xvasp.Directory+"/OUTCAR")){
      xOUTCAR OUTCAR_EP(xvasp.Directory+"/OUTCAR",true); //quiet, there might be issues with halfway-written OUTCARs
      E=OUTCAR_EP.efield_pead;
    }
    if(E.lrows!=1||E.urows!=3){E.clear();E.resize(3);}  //safety
    for(int i=E.lrows;i<=E.urows;i++){if(E[i]<=1.0E-6){E[i]=DEFAULT_EFIELD_PEAD;}}
    E*=0.2;  //dvalue used to be an input, not used right now, maybe make an aflowrc value in the future

    //[CO20210315 - OBSOLETE]E(1)=DEFAULT_EFIELD_PEAD;
    //[CO20210315 - OBSOLETE]E(2)=DEFAULT_EFIELD_PEAD;
    //[CO20210315 - OBSOLETE]E(3)=DEFAULT_EFIELD_PEAD;
    //[CO20210315 - OBSOLETE]vector<string> tokens;
    //[CO20210315 - OBSOLETE]// from vasp output file
    //[CO20210315 - OBSOLETE]if(0) for(uint i=1;i<=3;i++) {
    //[CO20210315 - OBSOLETE]  aurostd::string2tokens(aurostd::execute2string("cat "+xvasp.Directory+"/"+DEFAULT_VASP_OUT+" | grep E_g/N_"+aurostd::utype2string(i)),tokens,">");
    //[CO20210315 - OBSOLETE]  if(tokens.size()>0) {
    //[CO20210315 - OBSOLETE]    aurostd::string2tokens(string(tokens.at(1)),tokens,"=");
    //[CO20210315 - OBSOLETE]    if(tokens.size()>0) E(i)=aurostd::string2utype<double>(tokens.at(1));
    //[CO20210315 - OBSOLETE]    if(E(i)<=1.0E-6) E(i)=DEFAULT_EFIELD_PEAD;
    //[CO20210315 - OBSOLETE]  }
    //[CO20210315 - OBSOLETE]  //    cerr << "E(" << i << ")=" << E(i) << endl;
    //[CO20210315 - OBSOLETE]}
    //[CO20210315 - OBSOLETE]// from OUTCAR
    //[CO20210315 - OBSOLETE]aurostd::string2tokens(aurostd::execute2string("cat "+xvasp.Directory+"/OUTCAR | grep EFIELD_PEAD"),tokens,"=");
    //[CO20210315 - OBSOLETE]if(tokens.size()>0) {
    //[CO20210315 - OBSOLETE]  //   cerr << tokens.at(1) << endl;
    //[CO20210315 - OBSOLETE]  if(tokens.size()>0) aurostd::string2tokens(string(tokens.at(1)),tokens);
    //[CO20210315 - OBSOLETE]  if(tokens.size()>2) 
    //[CO20210315 - OBSOLETE]    for(uint i=1;i<=3;i++) {
    //[CO20210315 - OBSOLETE]      E(i)=aurostd::string2utype<double>(tokens.at(i-1));
    //[CO20210315 - OBSOLETE]      if(E(i)<=1.0E-6) E(i)=DEFAULT_EFIELD_PEAD;
    //[CO20210315 - OBSOLETE]      //	cerr << "E(" << i << ")=" << E(i) << endl;
    //[CO20210315 - OBSOLETE]    }
    //[CO20210315 - OBSOLETE]}

    string incar_input="EFIELD_PEAD="+aurostd::utype2string(E[1],6)+" "+aurostd::utype2string(E[2],6)+" "+aurostd::utype2string(E[3],6);
    if(aurostd::kvpair2bool(xvasp.INCAR,"EFIELD_PEAD","=")){
      string incar_input_old="EFIELD_PEAD="+aurostd::kvpair2string(xvasp.INCAR,"EFIELD_PEAD","=");
      if(incar_input==incar_input_old){return false;}
    }
    //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){return false;}  //remove whitespaces

    KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"EFIELD_PEAD",__AFLOW_FUNC__,vflags.KBIN_VASP_INCAR_VERBOSE);  //CO20200624
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing " << __AFLOW_FUNC__ << " [AFLOW] start" << endl;
    xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << __AFLOW_FUNC__ << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing " << __AFLOW_FUNC__ << " [AFLOW] end" << endl;

    //[CO20210315 - avoid writing orig]xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE); //CO20200624

    return true;
  }
}

namespace KBIN {
  bool XVASP_Afix_ApplyFix(const string& _fix,aurostd::xoption& xfixed,_xvasp& xvasp,_kflags& kflags,_vflags& vflags,_aflags &aflags,ofstream& FileMESSAGE) {  //CO20200624
    //CO20210315 - extensive rewrite
    //the schemes below check if the modification needs to be made (has it already been made?)
    //maintain this feedback system to ensure aflow doesn't keep spinning its wheels on the same fixes
    string function="KBIN::XVASP_Afix_ApplyFix";
    string fix=aurostd::toupper(_fix);
    string fix_str="("+fix+")";
    string operation=function+fix_str;
    string incar_input="";
    stringstream aus;
    bool VERBOSE=(FALSE || VERBOSE_MONITOR_VASP || XHOST.vflag_control.flag("MONITOR_VASP")==false);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //checks and quick return false here

    //check whether we need to apply once or can apply many times
    bool apply_once=true;
    if(fix=="EFIELD_PEAD"){apply_once=false;}
    else if(fix=="ENMAX"){apply_once=false;}
    else if(fix.find("KPOINTS+")!=string::npos||fix.find("KPOINTS-")!=string::npos){apply_once=false;}
    else if(fix.find("NBANDS+")!=string::npos||fix.find("NBANDS-")!=string::npos){apply_once=false;}
    else if(fix=="NBANDS"){apply_once=false;}
    else if(fix=="POTIM"){apply_once=false;}
    else if(fix=="RECYCLE_CONTCAR"){apply_once=false;}
    //add others here
    if(apply_once==true && xfixed.flag(fix)){
      if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": already applied" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      return false;
    }

    //check whether fix can be applied given current state
    bool Krun=true;
    //[CO20210315 - OBSOLETE]if(xvasp.NRELAXING<=xvasp.NRELAX){  //RELAX
    //[CO20210315 - OBSOLETE]  //add here any fixes that would not apply for relaxations
    //[CO20210315 - OBSOLETE]}else{ //STATIC, BANDS, etc.
    //[CO20210315 - OBSOLETE]  if(fix.find("ISMEAR=")!=string::npos){Krun=false;}  //specific smearing selected is very important
    //[CO20210315 - OBSOLETE]  //[CO20210315 - better to look for IBRION and NSW (linear-response calcs)]else if(fix=="POTIM"){Krun=false;}  //only applies to MD or ionic relaxations
    //[CO20210315 - OBSOLETE]}
    //[CO20210315 - OBSOLETE]if(Krun && (fix.find("RELAX_MODE")!=string::npos || fix=="POTIM")){ //check for IBRION/NSW
    //[CO20210315 - OBSOLETE]  VASP_Reread_INCAR(xvasp);
    //[CO20210315 - OBSOLETE]  if(!aurostd::kvpair2bool(xvasp.INCAR,"IBRION","=")){Krun=false;}
    //[CO20210315 - OBSOLETE]  else{ //STATIC also has IBRION, but NSW=0 so POTIM has no effect
    //[CO20210315 - OBSOLETE]    if(!aurostd::kvpair2bool(xvasp.INCAR,"NSW","=")){Krun=false;} //NSW=0 by default
    //[CO20210315 - OBSOLETE]    if(aurostd::kvpair2bool(xvasp.INCAR,"NSW","=") && aurostd::kvpair2utype<int>(xvasp.INCAR,"NSW","=")==0){Krun=false;}
    //[CO20210315 - OBSOLETE]  }
    //[CO20210315 - OBSOLETE]}
    //do not rely on xvasp.NRELAXING, etc., access to this information is limited for vasp monitor
    //instead, read INCAR and look for IBRION
    //for linear-response calcs, you also have IBRION, we must also look for NSW and make sure it is bigger than 0
    bool relaxing=false;
    VASP_Reread_INCAR(xvasp);
    if(aurostd::kvpair2bool(xvasp.INCAR,"IBRION","=")){
      if(aurostd::kvpair2bool(xvasp.INCAR,"NSW","=") && aurostd::kvpair2utype<int>(xvasp.INCAR,"NSW","=")>0){relaxing=true;}
    }
    if(relaxing){
      //add here any fixes that would not apply for relaxations
    }else{  //static, bands, etc.
      if(fix.find("ISMEAR=")!=string::npos){ //specific smearing selected is very important for static/bands
        if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static/bands calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
        Krun=false;
      }
      //[CO20210315 - can change KPOINTS for static]else if(fix.find("KPOINTS")!=string::npos){  //do not change the KPOINTS, only during relaxation
      //[CO20210315 - can change KPOINTS for static]  if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static/bands calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //[CO20210315 - can change KPOINTS for static]  Krun=false;
      //[CO20210315 - can change KPOINTS for static]}
      else if(fix.find("POSCAR")!=string::npos){  //do not change the POSCAR, only during relaxation
        if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static/bands calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
        Krun=false;
      }
      else if(fix=="POTIM"){  //has no effect for these runs
        if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static/bands calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
        Krun=false;
      }
      else if(fix.find("RECYCLE_CONTCAR")!=string::npos){  //do not change the POSCAR, only during relaxation
        if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static/bands calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
        Krun=false;
      }
      else if(fix.find("RELAX_MODE")!=string::npos){ //does not apply
        if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static/bands calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
        Krun=false;
      }
    }
    //some safety
    VASP_Reread_KPOINTS(xvasp);XVASP_KPOINTS_string2numbers(xvasp); //preload kpoints, load into xvasp.str
    if(fix.find("KPOINTS")!=string::npos && XVASP_KPOINTS_isAutoMesh(xvasp)==false){  //cannot change non-auto-mesh (bands)
      if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": cannot change KPOINTS for non-auto-mesh" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      Krun=false;
    }
    //[CO20210315 - OBSOLETE]//ismear=-5 is very important for static calcs
    //[CO20210315 - OBSOLETE]bool not_relaxing_or_bands=false; //looking for static-type calc, might improve this in the future
    //[CO20210315 - OBSOLETE]VASP_Reread_KPOINTS(xvasp);
    //[CO20210315 - OBSOLETE]if(relaxing==false && XVASP_KPOINTS_isAutoMesh(xvasp)){not_relaxing_or_bands=true;}  //bands is NOT auto-mesh
    //[CO20210315 - OBSOLETE]if(not_relaxing_or_bands){
    //[CO20210315 - OBSOLETE]  if(fix.find("ISMEAR=")!=string::npos){ //specific smearing selected is very important for static
    //[CO20210315 - OBSOLETE]    if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": does not apply for static calculations" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    //[CO20210315 - OBSOLETE]    Krun=false;
    //[CO20210315 - OBSOLETE]  }
    //[CO20210315 - OBSOLETE]}
    //add others here
    if(Krun==false){return false;}

    bool found_up_down=false;

    //check KPOINTS: avoid up and down (more than once)
    found_up_down=false; //if true, adds "KPOINTS_EXHAUSTED" to xfixed. we want to allow one up and down which might fix two different errors
    for(uint i=0;i<xfixed.vxscheme.size()&&Krun;i++){
      const string& scheme=xfixed.vxscheme[i];
      if(scheme.find("KPOINTS+")!=string::npos && fix.find("KPOINTS-")!=string::npos){
        if(xfixed.flag("KPOINTS_EXHAUSTED")){Krun=false;}
        else{found_up_down=true;}
      }
      if(scheme.find("KPOINTS-")!=string::npos && fix.find("KPOINTS+")!=string::npos){
        if(xfixed.flag("KPOINTS_EXHAUSTED")){Krun=false;}
        else{found_up_down=true;}
      }
    }
    if(found_up_down){xfixed.flag("KPOINTS_EXHAUSTED",true);}
    if(Krun==false){
      if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": avoiding KPOINTS++ and KPOINTS--" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      return false;
    }

    //check NBANDS: avoid up and down (more than once)
    found_up_down=false; //if true, adds "NBANDS_EXHAUSTED" to xfixed. we want to allow one up and down which might fix two different errors
    for(uint i=0;i<xfixed.vxscheme.size()&&Krun;i++){
      const string& scheme=xfixed.vxscheme[i];
      if(scheme.find("NBANDS+")!=string::npos && fix.find("NBANDS-")!=string::npos){
        if(xfixed.flag("NBANDS_EXHAUSTED")){Krun=false;}
        else{found_up_down=true;}
      }
      if(scheme.find("NBANDS-")!=string::npos && fix.find("NBANDS+")!=string::npos){
        if(xfixed.flag("NBANDS_EXHAUSTED")){Krun=false;}
        else{found_up_down=true;}
      }
    }
    if(found_up_down){xfixed.flag("NBANDS_EXHAUSTED",true);}
    if(Krun==false){
      if(VERBOSE){aus << "MMMMM  MESSAGE ignoring FIX=\"" << fix << "\"" << ": avoiding NBANDS++ and NBANDS--" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      return false;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////

    _xvasp xvasp_orig(xvasp);     //keep original, restore later if necessary
    _kflags kflags_orig(kflags);  //keep original, restore later if necessary
    _vflags vflags_orig(vflags);  //keep original, restore later if necessary

    vflags.KBIN_VASP_INCAR_VERBOSE=TRUE;  //VERBOSE, will restore original later
    string::size_type loc=0;

    string param_string="";
    double param_double=0.0;
    int param_int=0;
    xvector<double> param_xvd(3);

    //NB: xvasp.aopts.flag("FLAG::AFIX_DRYRUN") is very important
    //make sure that any changes made to files (INCAR, KPOINTS, aflow.in) are only done if xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false
    //otherwise, the monitor might make changes before they are needed

    if(fix.find("ALGO=")!=string::npos && fix.find("IALGO=")==string::npos) {
      if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved){Krun=false;} // don't touch kpoints
      if(Krun){
        loc=fix.find('=');
        if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ALGO mode not found",_INPUT_ILLEGAL_);}
        param_string=fix.substr(loc+1);  //get everything after the '='
      }

      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("ALGO",xvasp,vflags,param_string,0,0.0,FALSE));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - add fix to _AFLOWIN_
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        if(Krun && relaxing==true){ //the ALGO command in the aflow.in only affects relaxations
          KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]ALGO=",operation);
          KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]ALGO="+param_string,operation);
        }
      }
      //END - add fix to _AFLOWIN_ 
    }
    else if(fix.find("AMIN=")!=string::npos) {
      loc=fix.find('=');
      if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AMIN mode not found",_INPUT_ILLEGAL_);}
      param_double=aurostd::string2utype<double>(fix.substr(loc+1)); //get everything after the '='

      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("AMIN",xvasp,vflags,"",0,param_double,ON));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("AMIX=")!=string::npos) {
      loc=fix.find('=');
      if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AMIX mode not found",_INPUT_ILLEGAL_);}
      param_double=aurostd::string2utype<double>(fix.substr(loc+1)); //get everything after the '='

      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("AMIX",xvasp,vflags,"",0,param_double,ON));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("BMIX=")!=string::npos) {
      loc=fix.find('=');
      if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"BMIX mode not found",_INPUT_ILLEGAL_);}
      param_double=aurostd::string2utype<double>(fix.substr(loc+1)); //get everything after the '='

      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("BMIX",xvasp,vflags,"",0,param_double,ON));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="EDIFF") { //CO20210315
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      if(Krun){
        param_double=AUROSTD_MAX_DOUBLE;
        //try OUTCAR first
        if(param_double==AUROSTD_MAX_DOUBLE && aurostd::FileExist(xvasp.Directory+"/OUTCAR")&&aurostd::FileNotEmpty(xvasp.Directory+"/OUTCAR")){
          xOUTCAR OUTCAR(xvasp.Directory+"/OUTCAR",true); //quiet, there might be issues with halfway-written OUTCARs
          param_double=OUTCAR.EDIFF;
        }
        //try INCAR next
        if(param_double==AUROSTD_MAX_DOUBLE && aurostd::kvpair2bool(xvasp.INCAR,"EDIFF","=")){
          param_double=aurostd::kvpair2utype<double>(xvasp.INCAR,"EDIFF","=");
        }
        if(param_double==AUROSTD_MAX_DOUBLE){param_double=1e-4;}  //default for vasp
        param_double*=10; //increase by an order of magnitdue
      }
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("EDIFF",xvasp,vflags,"",0,param_double,FALSE));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"=" << param_double << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="EFIELD_PEAD") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_Afix_EFIELD_PEAD(xvasp,vflags,param_xvd));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"=[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(param_xvd,6)," ") << "]" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="ENMAX") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      if(Krun){
        param_double=AUROSTD_MAX_DOUBLE;
        if(aurostd::kvpair2bool(xvasp.INCAR,"ENMAX","=")){param_double=aurostd::kvpair2utype<int>(xvasp.INCAR,"ENMAX","=");}  //CO20210315 - the kvpair2bool check makes sure we don't get enmax=0 for lack of input
        else{param_double=DEFAULT_VASP_PREC_ENMAX_ACCURATE*xvasp.POTCAR_ENMAX;}
        if(param_double==AUROSTD_MAX_DOUBLE || aurostd::isequal(param_double,0.0)){param_double=500.0;}  //CO20210315
        param_double*=0.97; // reduce 3%.... if enough // param_double*=0.99; // reduce 1%.... if enough
        if(aurostd::isequal(xvasp.POTCAR_ENMAX,0.0)==false){
          if(param_double<DEFAULT_VASP_ENMAX_MINIMUM*xvasp.POTCAR_ENMAX){Krun=false;} //so it doesn't run forever
        }
        KBIN::XVASP_INCAR_PREPARE_GENERIC("LREAL",xvasp,vflags,"",0,0.0,ON);  //no Krun here, if it's already there, we're set
      }
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("ENMAX",xvasp,vflags,"",0,param_double,FALSE));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"=" << param_double << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("ISMEAR=")!=string::npos) {
      loc=fix.find('=');
      if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ISMEAR mode not found",_INPUT_ILLEGAL_);}
      param_int=aurostd::string2utype<int>(fix.substr(loc+1)); //get everything after the '='

      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("ISMEAR",xvasp,vflags,"",param_int,0.0,FALSE));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("KPOINTS")!=string::npos) {
      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")){Krun=false;} // don't touch kpoints
      if(Krun){
        string key="KPOINTS";
        loc=fix.find(key);
        if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"KPOINTS mode not found",_INPUT_ILLEGAL_);}
        string value="";
        value=fix.substr(loc+key.length());  //get everything after 'KPOINTS', including '='
        if(value=="=GAMMA"){param_string="Gamma_Shift";}
        else if(value=="=GAMMA_ODD"){param_string="Gamma,Xodd,Yodd,Zodd";}
        else if(value=="=GAMMA_EVEN"){param_string="Gamma,Xeven,Yeven,Zeven";}
        else if(value=="=KMAX"){param_string="Kmax";}
        else if(value=="++"){param_string="X++,Y++,Z++";}
        else if(value=="+=2"){param_string="X+=2,Y+=2,Z+=2";}
        else if(value=="--"){param_string="X--,Y--,Z--";}
        else if(value=="-=2"){param_string="X-=2,Y-=2,Z-=2";}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"KPOINTS unknown mode: \""+value+"\"",_INPUT_ILLEGAL_);}

        Krun=(Krun && VASP_Reread_KPOINTS(xvasp) && XVASP_KPOINTS_string2numbers(xvasp)); //preload kpoints, load into xvasp.str
        if(Krun && value=="=GAMMA" && KBIN::XVASP_KPOINTS_IncludesGamma(xvasp)){Krun=false;} //already done
      }

      if(Krun && VERBOSE){
        aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aus << "MMMMM  MESSAGE KPOINTS(pre )=[" << xvasp.str.kpoints_kscheme << ";" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ";" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << "]" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      //START - load KPOINTS into xvasp, modify xvasp.str.kpoints*, and write out new KPOINTS
      //[CO20210315 - move to before XVASP_KPOINTS_IncludesGamma()]Krun=(Krun && VASP_Reread_KPOINTS(xvasp) && XVASP_KPOINTS_string2numbers(xvasp)); //preload kpoints, load into xvasp.str
      Krun=(Krun && KBIN::XVASP_KPOINTS_OPERATION(xvasp,param_string));  //CO20210315
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS")));
      }
      //END - load KPOINTS into xvasp, modify xvasp.str.kpoints*, and write out new KPOINTS
      if(Krun && VERBOSE){
        aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aus << "MMMMM  MESSAGE KPOINTS(post)=[" << xvasp.str.kpoints_kscheme << ";" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ";" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << "]" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    else if(fix=="LREAL") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("LREAL",xvasp,vflags,"",0,0.0,ON));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("NBANDS")!=string::npos) {
      bool increase=(fix.find("--")==string::npos); //increase if "--" NOT found, so NBANDS or NBANDS++ works
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_Afix_NBANDS(xvasp,aflags,vflags,param_int,increase));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"=" << param_int << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="NCPUS") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\" kflags.KBIN_MPI_NCPUS(pre)=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //do not check for KBIN_MPI, since --monitor_vasp is running separate binary
      //we are coming from an MPI error, so we know this is running with NCPUS>1
      if(Krun){
        param_int=kflags.KBIN_MPI_NCPUS/2;
        if(param_int<1){Krun=false;}
      }
      if(Krun){
        kflags.KBIN_MPI_NCPUS_ORIG=kflags.KBIN_MPI_NCPUS;
        kflags.KBIN_MPI_NCPUS=param_int;
      }
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\" kflags.KBIN_MPI_NCPUS(post)=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="NCPUS_RESTORE") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\" kflags.KBIN_MPI_NCPUS(pre)=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //do not check for KBIN_MPI, since --monitor_vasp is running separate binary
      //we are coming from an MPI error, so we know this is running with NCPUS>1
      if(Krun){
        param_int=kflags.KBIN_MPI_NCPUS_ORIG;
        if(param_int<1){Krun=false;}
      }
      if(Krun){
        //[original is only set in "NCPUS"]kflags.KBIN_MPI_NCPUS_ORIG=kflags.KBIN_MPI_NCPUS;
        kflags.KBIN_MPI_NCPUS=param_int;
      }
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\" kflags.KBIN_MPI_NCPUS(post)=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="NELM") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_Afix_NELM(xvasp,vflags,param_int));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=" << fix << "=" << param_int << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="NPAR_REMOVE") { //before 'NPAR='
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      if(Krun){
        KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR",operation,vflags.KBIN_VASP_INCAR_VERBOSE); //CO20200624
      }
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("NPAR=")!=string::npos) {
      loc=fix.find('=');
      if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"NPAR mode not found",_INPUT_ILLEGAL_);}
      param_int=aurostd::string2utype<int>(fix.substr(loc+1)); //get everything after the '='

      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("NPAR",xvasp,vflags,"",param_int,0.0,FALSE));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="POTIM") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_Afix_POTIM(xvasp,vflags,param_double));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"=" << param_double << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("POSCAR_SCALE")!=string::npos) {
      if(xvasp.aopts.flag("FLAG::POSCAR_PRESERVED")){Krun=false;} // don't touch poscar
      if(Krun){
        string key="POSCAR_SCALE";
        loc=fix.find(key);
        if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR_SCALE mode not found",_INPUT_ILLEGAL_);}
        string value="";
        value=fix.substr(loc+key.length());  //get everything after 'POSCAR_SCALE', including '='
        if(value.find("*=")!=string::npos){param_string=value;aurostd::StringSubst(param_string,"*=","");param_double=aurostd::string2utype<double>(param_string);}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR_SCALE unknown mode: \""+value+"\"",_INPUT_ILLEGAL_);}
      }
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load POSCAR, modify xvasp.str, and write out new POSCAR
      Krun=(Krun && KBIN::VASP_Reread_POSCAR(xvasp)); //preload POSCAR
      if(Krun){xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);} //plug inside xvasp.str (WARNING: this overwrites xvasp.str_kpoints*)
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR.orig")));  //CO20210315 - POSCAR.orig here is NOT the original structure, but a saved state, come back later
      }
      if(Krun){
        if(VERBOSE){
          aus << "00000  MESSAGE WARNING changing scale" << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE BEFORE: volume=" << xvasp.str.GetVolume() << endl;
          aus << "00000  MESSAGE BEFORE: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      if(Krun){xvasp.str.scale*=param_double;}
      if(Krun){
        if(VERBOSE){
          aus << "00000  MESSAGE AFTER: volume=" << xvasp.str.GetVolume() << endl;
          aus << "00000  MESSAGE AFTER: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp)); //creates xvasp.POSCAR
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR")));
      }
      //END - load POSCAR, modify xvasp.str, and write out new POSCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix.find("POSCAR_VOLUME")!=string::npos) {
      if(xvasp.aopts.flag("FLAG::POSCAR_PRESERVED")){Krun=false;} // don't touch poscar
      if(Krun){
        string key="POSCAR_VOLUME";
        loc=fix.find(key);
        if(loc==fix.size()-1||loc==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR_VOLUME mode not found",_INPUT_ILLEGAL_);}
        string value="";
        value=fix.substr(loc+key.length());  //get everything after 'POSCAR_VOLUME', including '='
        if(value.find("*=")!=string::npos){param_string=value;aurostd::StringSubst(param_string,"*=","");param_double=aurostd::string2utype<double>(param_string);}
        else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR_VOLUME unknown mode: \""+value+"\"",_INPUT_ILLEGAL_);}
      }
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load POSCAR, modify xvasp.str, and write out new POSCAR
      Krun=(Krun && KBIN::VASP_Reread_POSCAR(xvasp)); //preload POSCAR
      if(Krun){xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);} //plug inside xvasp.str (WARNING: this overwrites xvasp.str_kpoints*)
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR.orig")));  //CO20210315 - POSCAR.orig here is NOT the original structure, but a saved state, come back later
      }
      if(Krun){
        if(VERBOSE){
          aus << "00000  MESSAGE WARNING changing volume" << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE BEFORE: volume=" << xvasp.str.GetVolume() << endl;
          aus << "00000  MESSAGE BEFORE: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      if(Krun){xvasp.str.scale*=std::pow(param_double,(1.0/3.0));}
      if(Krun){
        if(VERBOSE){
          aus << "00000  MESSAGE AFTER: volume=" << xvasp.str.GetVolume() << endl;
          aus << "00000  MESSAGE AFTER: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp)); //creates xvasp.POSCAR
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR")));
      }
      //END - load POSCAR, modify xvasp.str, and write out new POSCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="POSCAR=STANDARD_CONVENTIONAL") {
      if(xvasp.aopts.flag("FLAG::POSCAR_PRESERVED")){Krun=false;} // don't touch poscar
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load POSCAR, modify xvasp.str, and write out new POSCAR
      Krun=(Krun && KBIN::VASP_Reread_POSCAR(xvasp)); //preload POSCAR
      if(Krun){xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);} //plug inside xvasp.str (WARNING: this overwrites xvasp.str_kpoints*)
      if(Krun && xvasp.str.Standard_Lattice_conventional){Krun=false;}  //DX is this right?  //CO20210315 - OBSOLETE as we read in fresh POSCAR above, so this setting will not be set, we need a function here
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR.orig")));  //CO20210315 - POSCAR.orig here is NOT the original structure, but a saved state, come back later
      }
      if(Krun){
        //CO20210315 - need to make sure fix magmom too (if it's there)
        //load into xvasp.str.atoms[i].order_parameter_atom/value, convert to sconv, then print out
        VASP_Reread_INCAR(xvasp);  //preload incar - MAGMOM
        bool remove_magmom=aurostd::kvpair2bool(xvasp.INCAR,"MAGMOM","=");
        bool write_magmom=XVASP_INCAR_Read_MAGMOM(xvasp);
        if(VERBOSE){
          aus << "00000  MESSAGE WARNING changing to the conventional cell" << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE BEFORE: a,b,c,alpha,beta,gamma=" << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
          aus << "00000  MESSAGE BEFORE: volume=" << xvasp.str.GetVolume() << endl;
          if(write_magmom){aus << "00000  MESSAGE BEFORE: MAGMOM=" << aurostd::kvpair2string(xvasp.INCAR,"MAGMOM","=") << endl;}
          aus << "00000  MESSAGE BEFORE: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        //
        xvasp.str.Standard_Conventional_UnitCellForm(); //CO20210315 - previously missing before
        if(remove_magmom){KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"MAGMOM",operation,vflags.KBIN_VASP_INCAR_VERBOSE);} //CO20200624
        if(write_magmom){KBIN::XVASP_INCAR_PREPARE_GENERIC("AUTO_MAGMOM",xvasp,vflags,"",0,0.0,FALSE);}
        if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
          if(remove_magmom || write_magmom){aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));}  //write out incar
        }
        //
        if(VERBOSE){
          aus << "00000  MESSAGE AFTER: a,b,c,alpha,beta,gamma=" << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
          aus << "00000  MESSAGE AFTER: volume=" << xvasp.str.GetVolume() << endl;
          if(write_magmom){aus << "00000  MESSAGE AFTER: MAGMOM=" << aurostd::kvpair2string(xvasp.INCAR,"MAGMOM","=") << endl;}
          aus << "00000  MESSAGE AFTER: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp)); //CO20210315 - previously missing before  //creates xvasp.POSCAR
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR")));
      }
      //END - load POSCAR, modify xvasp.str, and write out new POSCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - add fix to _AFLOWIN_ 
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        if(Krun){
          KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=",operation);
          KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SCONV",operation);
        }
      }
      //END - add fix to _AFLOWIN_ 
    }
    else if(fix=="RECYCLE_CONTCAR") {
      if(xvasp.aopts.flag("FLAG::POSCAR_PRESERVED")){Krun=false;} // don't touch poscar
      param_string=xvasp.Directory+string("/CONTCAR");
      Krun=(Krun && aurostd::FileExist(param_string) && !aurostd::FileEmpty(param_string));
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load POSCAR, modify xvasp.str, and write out new POSCAR
      Krun=(Krun && KBIN::VASP_Reread_POSCAR(xvasp)); //preload POSCAR
      if(Krun){xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);} //plug inside xvasp.str (WARNING: this overwrites xvasp.str_kpoints*)
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR.orig")));  //CO20210315 - POSCAR.orig here is NOT the original structure, but a saved state, come back later
      }
      if(Krun){
        if(VERBOSE){
          aus << "00000  MESSAGE WARNING recycling CONTCAR" << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "00000  MESSAGE BEFORE: a,b,c,alpha,beta,gamma=" << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
          aus << "00000  MESSAGE BEFORE: volume=" << xvasp.str.GetVolume() << endl;
          aus << "00000  MESSAGE BEFORE: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      if(Krun){xvasp.str=xstructure(param_string,IOVASP_POSCAR);}
      if(Krun){
        if(VERBOSE){
          aus << "00000  MESSAGE AFTER: a,b,c,alpha,beta,gamma=" << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
          aus << "00000  MESSAGE AFTER: volume=" << xvasp.str.GetVolume() << endl;
          aus << "00000  MESSAGE AFTER: structure: " << endl;
          aus << xvasp.str;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
      Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp)); //CO20210315 - previously missing before  //creates xvasp.POSCAR
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR")));
      }
      //END - load POSCAR, modify xvasp.str, and write out new POSCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - add fix to _AFLOWIN_ 
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        if(Krun){KBIN::VASP_CONTCAR_Save(xvasp);}
      }
      //END - add fix to _AFLOWIN_ 
    }
    else if(fix=="RELAX_MODE=FORCES") {
      //do not move to XVASP_INCAR_PREPARE_GENERIC, this fix will affect many parameters
      //aflow relies on specialized function to initialize INCAR with these kind of settings (XVASP_INCAR_Relax_ON())
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      if(Krun){
        //refer to vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES" above
        //we're probably coming from RELAX=ENERGY, so only change non-overlapping settings
        //first do a check that we're not already doing RELAX_MODE=FORCES, IBRION is a good enough check
        incar_input="IBRION=1";
        if(aurostd::kvpair2bool(xvasp.INCAR,"IBRION","=")){
          string incar_input_old="IBRION="+aurostd::kvpair2string(xvasp.INCAR,"IBRION","=");
          if(incar_input==incar_input_old){Krun=false;}
        }
        //[CO20210315 - substring2bool() can match ENMAX=2 with ENMAX=20]if(aurostd::substring2bool(xvasp.INCAR,incar_input,true)){Krun=false;}  //remove whitespace
      }
      if(Krun){
        //REMOVE LINES
        KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"IBRION,NELMIN,ADDGRID,EDIFFG,NSW",operation,vflags.KBIN_VASP_INCAR_VERBOSE); //CO20200624
        //ADD LINES
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] start" << endl;
        xvasp.INCAR << aurostd::PaddedPOST(incar_input,_incarpad_) << " # " << operation << endl; //CO20200624
        xvasp.INCAR << aurostd::PaddedPOST("NELMIN=4",_incarpad_) << " # " << operation << endl; //CO20200624
        xvasp.INCAR << aurostd::PaddedPOST("ADDGRID=.TRUE.",_incarpad_) << " # " << operation << endl; //CO20200624
        xvasp.INCAR << aurostd::PaddedPOST("EDIFFG="+aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG,12),_incarpad_) << " # " << operation << endl; //CO20200624
        xvasp.INCAR << aurostd::PaddedPOST("NSW=100",_incarpad_) << " # " << operation << endl; //CO20200624
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing " << operation << " [AFLOW] end" << endl;
      }
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - add fix to _AFLOWIN_
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        if(Krun){
          KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]RELAX_MODE=",operation);
          KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]RELAX_MODE=FORCES",operation);
        }
      }
      //END - add fix to _AFLOWIN_ 
    }
    else if(fix=="RESTART_CALC") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="SKIP_RUN") { //CO20210315
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - one-off solution for memory
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        if(Krun){
          stringstream command("");
          command << "cat " << xvasp.Directory << "/" << DEFAULT_VASP_OUT << " | grep AFLOW > " << xvasp.Directory << "/" << DEFAULT_AFLOW_MEMORY_OUT << endl;
          command << "cat " << xvasp.Directory << "/" << DEFAULT_VASP_OUT << " | grep AFLOW > " << xvasp.Directory << "/SKIP" << endl;	
          command << "cat " << xvasp.Directory << "/" << DEFAULT_VASP_OUT << " | grep AFLOW >> " << xvasp.Directory << DEFAULT_AFLOW_ERVASP_OUT << endl;	
          aurostd::execute(command);
        }
      }
      //END - one-off solution for memory
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="SYMPREC") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("SYMPREC",xvasp,vflags,"",0,0.0,OFF));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else if(fix=="SYM=OFF") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - load INCAR into xvasp, modify, then write out new INCAR
      Krun=(Krun && VASP_Reread_INCAR(xvasp));  //preload incar
      Krun=(Krun && KBIN::XVASP_INCAR_PREPARE_GENERIC("SYM",xvasp,vflags,"",0,0.0,OFF));
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"))); //write out incar
      }
      //END - load INCAR into xvasp, modify, then write out new INCAR
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      //START - add fix to _AFLOWIN_
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        if(Krun){
          KBIN::AFLOWIN_REMOVE(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]SYM=",operation);
          KBIN::AFLOWIN_ADD(xvasp.Directory+"/"+_AFLOWIN_,"[VASP_FORCE_OPTION]SYM=OFF",operation);
        }
      }
      //END - add fix to _AFLOWIN_ 
    }
    else if(fix=="ULIMIT") {
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE attempting FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      if(Krun){kflags.KBIN_MPI_OPTIONS=string("ulimit -s unlimited");} //+string(" && ")+kflags.KBIN_MPI_OPTIONS;
      if(Krun && VERBOSE){aus << "MMMMM  MESSAGE applied FIX=\"" << fix << "\"" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown fix: \""+fix+"\"",_INPUT_ILLEGAL_);  //CO20210315
    }

    if(Krun==false){
      xvasp=xvasp_orig;     //restore original
      kflags=kflags_orig;   //restore original
      vflags=vflags_orig;   //restore original
      return false;
    }

    //if the following fixes were applied, reset xfixed for ALGO, some fixes that did not work previously might work now
    bool reset_xfixes=false;
    if(fix.find("KPOINTS")!=string::npos){reset_xfixes=true;}
    else if(fix.find("POSCAR_SCALE")!=string::npos){reset_xfixes=true;}
    else if(fix.find("POSCAR_VOLUME")!=string::npos){reset_xfixes=true;}
    else if(fix.find("RECYCLE_CONTCAR")!=string::npos){reset_xfixes=true;}
    //add others here
    if(reset_xfixes){
      for(uint i=xfixed.vxscheme.size()-1;i<xfixed.vxscheme.size();i--){ //go backwards
        const string& scheme=xfixed.vxscheme[i];
        if(scheme.find("ALGO=")!=string::npos){
          if(VERBOSE){aus << "MMMMM  MESSAGE clearing \"" << scheme << "\" from xfixed (FIX=\"" << fix << "\")" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
          xfixed.vxscheme.erase(xfixed.vxscheme.begin()+i); //erase is really inefficient, but it is okay here because xfixed.vxscheme is not large
        }
      }
    }

    vflags.KBIN_VASP_INCAR_VERBOSE=vflags_orig.KBIN_VASP_INCAR_VERBOSE; //restore original

    xfixed.flag(fix,true);  //always set, whether we check depends on the fix
    return true;  //CO20200624
  }
}

namespace KBIN {
  bool XVASP_Afix_IgnoreFix(const string& _fix,const _vflags& vflags){
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_IVASP_ || XHOST.DEBUG);

    if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ALL")){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ALL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+_fix) << endl;}
      return true;
    }

    if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+_fix)){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:"+_fix+"\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+_fix) << endl;}
      return true;
    }

    //try removing all math first
    string fix=_fix;
    string::size_type loc=0;
    loc=fix.find("=");if(loc!=string::npos){fix=fix.substr(0,loc);}
    loc=fix.find("-");if(loc!=string::npos){fix=fix.substr(0,loc);}
    loc=fix.find("+");if(loc!=string::npos){fix=fix.substr(0,loc);}
    loc=fix.find("*");if(loc!=string::npos){fix=fix.substr(0,loc);}
    loc=fix.find("/");if(loc!=string::npos){fix=fix.substr(0,loc);}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " generic_fix_string[1]=" << fix << endl;}

    if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+fix)){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:"+fix+"\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+fix) << endl;}
      return true;
    }

    //also try removing _, not good for RELAX_MODE or SKIP_RUN, but good for POSCAR_VOLUME
    loc=fix.find("_");if(loc!=string::npos){fix=fix.substr(0,loc);}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " generic_fix_string[2]=" << fix << endl;}

    if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+fix)){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:"+fix+"\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+fix) << endl;}
      return true;
    }

    //nparc and nparn
    aurostd::StringSubst(fix,"NPARC","NPAR");
    aurostd::StringSubst(fix,"NPARN","NPAR");

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " generic_fix_string[3]=" << fix << endl;}

    if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+fix)){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:"+fix+"\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:"+fix) << endl;}
      return true;
    }

    return false;
  }
} // namespace KBIN

namespace KBIN {
  bool XVASP_Afix(const string& mode,int& submode,bool try_last_ditch_efforts,aurostd::xoption& xfixed,_xvasp& xvasp,_kflags& kflags,_vflags& vflags,_aflags &aflags,ofstream& FileMESSAGE) {  //CO20200624 - adding submode
    //CO20210315 - extensive rewrite
    //the schemes below check if the modification needs to be made (has it already been made?)
    //maintain this feedback system to ensure aflow doesn't keep spinning its wheels on the same fixes
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_IVASP_ || XHOST.DEBUG);
    string function="KBIN::XVASP_Afix";
    string file_error="",incar_input="";
    stringstream aus;

    int submode_orig=submode;     //keep original, restore later if necessary
    _xvasp xvasp_orig(xvasp);     //keep original, restore later if necessary
    _kflags kflags_orig(kflags);  //keep original, restore later if necessary
    _vflags vflags_orig(vflags);  //keep original, restore later if necessary

    vflags.KBIN_VASP_INCAR_VERBOSE=TRUE;  //VERBOSE, will restore original later
    bool Krun=true;
    int submode_increment=1;

    string fix="";
    file_error="aflow.error."+aurostd::tolower(mode); //if mode is empty we throw later

    //some intuition about the schemes below
    //changing KPOINTS will remove any 'ALGO' from xfixed, because they might work with a different set of KPOINTS
    //therefore, the submodes are set up so ALGO changes go first and KPOINTS changes occur later
    //an exception is KPOINTS=GAMMA, always try this first
    //ending with KPOINTS++ is similar to resetting submode to 0
    //try GAMMA_EVEN before GAMMA_ODD, as later solutions would benefit from and odd scheme
    //both GAMMA_EVEN and GAMMA_ODD are expected to increase at least some of the KPOINTS

    if(mode=="BRMIX"){
      //https://www.vasp.at/forum/viewtopic.php?f=3&t=1417 (BRMIX)
      //https://www.vasp.at/forum/viewtopic.php?t=7949 (BRMIX)
      //CO20210315 - from several tests, "ALGO=VERYFAST" is the only solution that has worked
      //keep this as the only solution and then try last-ditch efforts
      fix="ALGO=VERYFAST";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="CSLOSHING") {
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //try ALGO=NORMAL //CO20200624
        fix="ALGO=NORMAL";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==1){ //try ALGO=FAST //CO20200624
        fix="ALGO=FAST";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=2){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="DAV" || mode=="EDDDAV") {
      //https://www.vasp.at/forum/viewtopic.php?t=1192&p=2321 (DAV)
      //https://www.vasp.at/forum/viewtopic.php?f=2&t=10409&p=10434#:~:text=7%3A18%20am-,Re%3A%20on%20solving%20%22Error%20EDDDAV%3A%20Call%20to%20ZHEGV%20failed,Returncode%20%3D%20xx%22&text=This%20is%20an%20error%20of,of%20the%20very%20last%20one. (EDDDAV)
      //problems with DAV, avoid NORMAL
      //solution: try flipping around different ALGOs (EDDDAV)
      //previously it was IALGO=48, better to use ALGO=VERYFAST (same thing: https://www.vasp.at/wiki/index.php/ALGO) (EDDDAV)
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //try ALGO=VERYFAST //CO20200624
        fix="ALGO=VERYFAST";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==1){ //try ALGO=FAST //CO20200624
        fix="ALGO=FAST";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=2){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="DENTET") {
      //https://www.vasp.at/forum/viewtopic.php?f=3&t=416
      fix="ISMEAR=2";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="EFIELD_PEAD") {
      fix="EFIELD_PEAD";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="EXCCOR") {
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //try ALGO=VERYFAST //CO20200624
        fix="POSCAR_SCALE*=1.2"; //volume*=1.2^3
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==1){  //CO20210315 - try lowering NBANDS, worked for LIB2/LIB/Ca_svCr_pv/64
        fix="NBANDS--";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=2){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="GAMMA_SHIFT" || mode=="IBZKPT") {  //CO20210315 - does a simple shift work for IBZKPT? come back
      fix="KPOINTS=GAMMA";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="KKSYM") {
      //CO20210315 - tried SYMPREC and KPOINTS=GAMMA first in several tests, they never worked. KPOINTS=KMAX works every time
      //ROTMAT is trigged for KKSYM, it will try these other procedures anyway in the off change KPOINTS=KMAX doesn't work
      fix="KPOINTS=KMAX";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="LREAL") {
      fix="LREAL";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="MEMORY") { //CO20210315
      //for memory issues, always try to save the CONTCAR, it should be a good starting point for the next calc
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //lower NCPUS
        bool Krun1=true;fix="ULIMIT"; //always apply ULIMIT for memory stuff
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="NCPUS";
        bool ignorefix2=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix2){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun3=true;fix="RECYCLE_CONTCAR";  //recycle contcar if possible
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun3=false;}
        Krun3=((Krun1||Krun2) && Krun3 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));  //only run if (Krun1||Krun2)

        Krun=(Krun1||Krun2||Krun3);
        if(!Krun){  //remove fixes before going to next submode
          if(!ignorefix2){fix="NCPUS_RESTORE";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
        }
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==1){ //lower NBANDS, try before lowering k-points
        bool Krun1=true;fix="ULIMIT"; //always apply ULIMIT for memory stuff
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="NBANDS--";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun3=true;fix="RECYCLE_CONTCAR";  //recycle contcar if possible
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun3=false;}
        Krun3=((Krun1||Krun2) && Krun3 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));  //only run if (Krun1||Krun2)

        Krun=(Krun1||Krun2||Krun3);
        //no need to remove these settings if it fails
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==2){ //lower k-points
        bool Krun1=true;fix="ULIMIT"; //always apply ULIMIT for memory stuff
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="KPOINTS--";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun3=true;fix="RECYCLE_CONTCAR";  //recycle contcar if possible
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun3=false;}
        Krun3=((Krun1||Krun2) && Krun3 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));  //only run if (Krun1||Krun2)

        Krun=(Krun1||Krun2||Krun3);
        //no need to remove these settings if it fails
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      //we might add lowering NGX, need to test...
      if(submode==3){ //skip run
        fix="SKIP_RUN";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=4){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="MPICH11") {
      fix="ULIMIT"; //always apply ULIMIT for memory stuff
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="MPICH139") {
      bool Krun1=true;fix="ULIMIT"; //always apply ULIMIT for memory stuff
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun1=false;}
      Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

      bool Krun2=true;fix="KPOINTS-=2";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun2=false;}
      Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

      Krun=(Krun1||Krun2);
    }
    else if(mode=="MPICH174") { //CO20210315
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //lower NCPUS
        bool Krun1=true;fix="ULIMIT"; //always apply ULIMIT for memory stuff
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="NCPUS";
        bool ignorefix2=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix2){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        Krun=(Krun1||Krun2);
        if(!Krun){  //remove fixes before going to next submode
          if(!ignorefix2){fix="NCPUS_RESTORE";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
        }
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      //this might be a memory issue, might not be
      //so let's not waste time lowering NBANDS and KPOINTS, just fail
      if(submode>=1){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="NATOMS") {
      fix="POSCAR_VOLUME*=2";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="NBANDS") {
      fix="NBANDS++";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      if(Krun){xfixed.flag("NBANDS_EXHAUSTED",false);}  //if the NBANDS error comes up, the ONLY solution is to increase NBANDS //let it continue resetting "NBANDS_EXHAUSTED", no problem... //remember, there are other reasons for fix=="NBANDS++"
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="NELM") { //CSLOSHING solutions should be tried first
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //AMIX/BMIX fixes, helps here: ICSD/LIB/CUB/Pd1Tm1_ICSD_649071
        bool Krun1=true;fix="AMIX=0.1";
        bool ignorefix1=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix1){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="BMIX=0.01";
        bool ignorefix2=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix2){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        Krun=(Krun1||Krun2);
        if(!Krun){  //remove fixes before going to next submode
          if(!ignorefix1){fix="AMIX=-1";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
          if(!ignorefix2){fix="BMIX=-1";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
        }
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==1){ //BMIX/AMIN fixes
        bool Krun1=true;fix="BMIX=3"; //3.0
        bool ignorefix1=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix1){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="AMIN=0.01";
        bool ignorefix2=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix2){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        Krun=(Krun1||Krun2);
        if(!Krun){  //remove fixes before going to next submode
          if(!ignorefix1){fix="BMIX=-1";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
          if(!ignorefix2){fix="AMIN=-1";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
        }
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==2){ //desperate attempt, increase NELM  //CO20211017 - do BEFORE you increase KPOINTS (KPOINTS=GAMMA_ODD), as the longer run will be more expensive with more KPOINTS
        fix="NELM";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==3){ //Gamma-odd (might increase ki) //CO20211017 - worked for Cr_pvHf_pvMo_pvNV_svZr_sv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xD_P1-0.2xA-0.2xB-0.2xC-0.2xE-0.2xF/ARUN.POCC_05_H0C4
        fix="KPOINTS=GAMMA_ODD";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==4){ //desperate attempt 2, also increase EDIFF, sometimes the threshold is just a bit too high
        bool Krun1=true;fix="NELM";
        bool ignorefix1=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix1){Krun1=false;}
        Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        bool Krun2=true;fix="EDIFF";
        bool ignorefix2=XVASP_Afix_IgnoreFix(fix,vflags);
        if(ignorefix2){Krun2=false;}
        Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

        Krun=(Krun1||Krun2);
        //no need to remove these settings if it fails
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=5){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="NKXYZ_IKPTD") {
      //https://www.vasp.at/forum/viewtopic.php?t=1228
      fix="KPOINTS--";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="NPAR") {
      fix="NPAR=1";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="NPARC" || mode=="NPARN") {
      fix="NPAR=4";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="NPAR_REMOVE") {
      fix="NPAR_REMOVE";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="PSMAXN") {
      fix="ENMAX";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="READ_KPOINTS_RD_SYM") {
      fix="SYM=OFF";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="RECYCLE_CONTCAR") {
      fix="RECYCLE_CONTCAR";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="RESTART_CALC") {
      fix="RESTART_CALC";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="RMM_DIIS") {
      //problems with RMM-DIIS, avoid VERYFAST
      //https://www.vasp.at/forum/viewtopic.php?t=214 (EDDRMM)
      //https://www.vasp.at/forum/viewtopic.php?f=3&t=18028 (NUM_PROB)
      //https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //CO20200624 - try ALGO=FAST (fast is faster than normal, try first)
        //https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS
        fix="ALGO=FAST";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      //else if(submode==1) //DO NOT USE else if, we may change submode inside
      if(submode==1){ //CO20200624 - try ALGO=NORMAL
        fix="ALGO=NORMAL";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==2){ //Gamma-odd (might increase ki) //CO20200624
        fix="KPOINTS=GAMMA_ODD";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==3){ //CO20200624 - try RELAX_MODE=FORCES
        fix="RELAX_MODE=FORCES";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==4){ //CO20200624 - try fixing POTIM
        fix="POTIM";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==5){ //CO20200624 - try Fermi smearing
        fix="ISMEAR=-1";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=6){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="ROTMAT") {
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){ //try SYMPREC first (easy)
        //SYMPREC has been shown to work: LIB2/LIB/CeO/AB_tI4_139_b_a-001.AB
        fix="SYMPREC";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      //else if(submode==1) //DO NOT USE else if, we may change submode inside
      if(submode==1){ //Gamma without increasing KPOINTS  //CO20210315 - does a simple shift fix ROTMAT? come back
        fix="KPOINTS=GAMMA";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==2){ //Gamma-even  //CO20210315 - does this ever work? or is it just increase KPOINTS helps? come back
        fix="KPOINTS=GAMMA_EVEN";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==3){ //Gamma-odd
        fix="KPOINTS=GAMMA_ODD";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==4){ //Kmax
        fix="KPOINTS=KMAX";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==5){ //ISYM=0
        fix="SYM=OFF";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=6){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else if(mode=="SYMPREC") {
      fix="SYMPREC";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="THREADS") {
      bool Krun1=true;fix="ULIMIT"; //always apply ULIMIT for memory stuff
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun1=false;}
      Krun1=(Krun1 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

      bool Krun2=true;fix="NCPUS";
      bool ignorefix2=XVASP_Afix_IgnoreFix(fix,vflags);
      if(ignorefix2){Krun2=false;}
      Krun2=(Krun2 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));

      bool Krun3=true;fix="RECYCLE_CONTCAR";  //recycle contcar if possible
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun3=false;}
      Krun3=((Krun1||Krun2) && Krun3 && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));  //only run if (Krun1||Krun2)

      Krun=(Krun1||Krun2||Krun3);
      if(!Krun){  //remove fixes before going to next submode
        if(!ignorefix2){fix="NCPUS_RESTORE";XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE);xfixed.flag(fix,false);}
      }
    }
    else if(mode=="ZBRENT") { //other RMM-DIIS patches will follow
      //https://www.vasp.at/forum/viewtopic.php?t=1856
      fix="RELAX_MODE=FORCES";
      if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }
    else if(mode=="ZPOTRF") {
      if(submode<0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no submode set: \""+mode+"\"",_INPUT_ILLEGAL_);}  //CO20210315
      if(submode==0){  //CO20210315 - POTIM usually works
        fix="POTIM";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode==1){  //CO20210315 - increasing volume might get you out of rut
        fix="POSCAR_VOLUME*=1.05";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      //previously we tried "POSCAR=STANDARD_CONVENTIONAL", but it was turned off and has not been shown to work
      //[CO20210315 - doesn't work]//we keep here as a last-ditch effort
      //[CO20210315 - doesn't work]if(submode==2){  //CO20210315 - try converting to standard conventional (previously not applied)
      //[CO20210315 - doesn't work]  fix="POSCAR=STANDARD_CONVENTIONAL";
      //[CO20210315 - doesn't work]  if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
      //[CO20210315 - doesn't work]  Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
      //[CO20210315 - doesn't work]  if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      //[CO20210315 - doesn't work]}
      if(submode==2){  //CO20210315 - try lowering NBANDS
        fix="NBANDS--";
        if(XVASP_Afix_IgnoreFix(fix,vflags)){Krun=false;}
        Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
        if(!Krun){Krun=true;submode++;} //reset and go to the next solution
      }
      if(submode>=3){Krun=false;}
      submode+=submode_increment;submode_increment=1;  //increment and reset
    }
    else{
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"unknown mode: \""+mode+"\"",_INPUT_ILLEGAL_);  //CO20210315
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //last-ditch efforts

    //try_last_ditch_efforts controls whether we try ANY last_ditch_efforts
    //the local variable try_last_ditch_effort targets each individual effort
    //think XHOST.DBEUG vs. LDEBUG
    bool try_last_ditch_effort=try_last_ditch_efforts;

    if(LDEBUG){aus << __AFLOW_FUNC__ << " Krun=" << Krun << " [1]" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    try_last_ditch_effort=try_last_ditch_efforts;
    if(mode=="MEMORY") {try_last_ditch_effort=false;} //changing POSCAR doesn't help
    else if(mode=="MPICH11") {try_last_ditch_effort=false;} //changing POSCAR doesn't help
    else if(mode=="MPICH139") {try_last_ditch_effort=false;} //changing POSCAR doesn't help
    else if(mode=="MPICH174") {try_last_ditch_effort=false;} //changing POSCAR doesn't help
    else if(mode=="NBANDS") {try_last_ditch_effort=false;} //changing POSCAR doesn't help

    if(Krun==false && try_last_ditch_effort){
      //last-ditch effort, increase volume
      fix="POSCAR_VOLUME*=1.05";  //this is a great "general" solution
      Krun=(xfixed.flag(fix)==false); //only try if it has not already been tried before
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }

    if(LDEBUG){aus << __AFLOW_FUNC__ << " Krun=" << Krun << " [2]" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    try_last_ditch_effort=try_last_ditch_efforts;
    if(mode=="EXCCOR") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MEMORY") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MPICH11") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MPICH139") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MPICH174") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="NATOMS") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="NBANDS") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help

    if(Krun==false && try_last_ditch_effort){
      //last-ditch effort, increase KPOINTS
      fix="KPOINTS=GAMMA_ODD";  //this is a great "general" solution
      Krun=(xfixed.flag(fix)==false); //only try if it has not already been tried before
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }

    if(LDEBUG){aus << __AFLOW_FUNC__ << " Krun=" << Krun << " [3]" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    try_last_ditch_effort=try_last_ditch_efforts;
    if(mode=="EXCCOR") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MEMORY") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MPICH11") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MPICH139") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="MPICH174") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="NATOMS") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help
    else if(mode=="NBANDS") {try_last_ditch_effort=false;} //changing KPOINTS doesn't help

    if(Krun==false && try_last_ditch_effort){
      //last-ditch effort, increase KPOINTS
      fix="KPOINTS++";
      Krun=(xfixed.flag(fix)==false); //only try if it has not already been tried before
      Krun=(Krun && XVASP_Afix_ApplyFix(fix,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE));
    }

    if(LDEBUG){aus << __AFLOW_FUNC__ << " Krun=" << Krun << " [4]" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    ///////////////////////////////////////////////////////////////////////////////////////////////

    vflags.KBIN_VASP_INCAR_VERBOSE=vflags_orig.KBIN_VASP_INCAR_VERBOSE; //restore original (always)

    if(Krun==false||xvasp.aopts.flag("FLAG::AFIX_DRYRUN")){ //restore original for dry-run too
      xvasp=xvasp_orig;     //restore original
      kflags=kflags_orig;   //restore original
      vflags=vflags_orig;   //restore original
    }

    if(Krun==false){
      submode=submode_orig; //restore original
      return false;
    }

    if(file_error.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No file_error specified",_INPUT_ILLEGAL_);} //CO20200624
    if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
      KBIN::XVASP_Afix_Clean(xvasp,file_error);
    }
    return true;  //CO20200624
  }
}
// ---------------------------------------------------------------------------------------------------------------------------------------------------------

// ***************************************************************************//
// KBIN::GetMostRelaxedStructure
// ***************************************************************************//
namespace KBIN {
  xstructure GetMostRelaxedStructure(string directory) {
    string POSCARfile;

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "XHOST.vext.size()!=XHOST.vcat.size(), aborting", _RUNTIME_ERROR_); //CO20200404
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //READ POSCAR.orig
    string file_poscar_tmp=aurostd::TmpFileCreate("POSCAR.tmp");

    //GET its lattice vectors and reciprocal lattice vectors
    bool found_POSCAR=FALSE;
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.bands";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.static")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.static";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax2")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.relax2";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax1")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.relax1";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR")) { 
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR";
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.bands"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.static"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.static"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax1"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.relax1"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax2"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.relax2"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    if(!found_POSCAR)  {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No POSCAR[.bands|.static|.relax2|.relax1][.EXT] found in the directory, aborting", _INPUT_MISSING_); //CO20200404
    }

    xstructure xstr_name(POSCARfile, IOVASP_POSCAR); //import xstructure from filename
    if(aurostd::FileExist(file_poscar_tmp)) aurostd::RemoveFile(file_poscar_tmp);
    return xstr_name;
  }
}

// ***************************************************************************//
// KBIN::ExtractAtomicSpecies
// ***************************************************************************//
namespace KBIN {
  vector<string> ExtractAtomicSpecies(const string& directory,ostream& oss) {ofstream FileMESSAGE;return ExtractAtomicSpecies(directory,FileMESSAGE,oss);} //CO20200404 - added ofstream
  vector<string> ExtractAtomicSpecies(const string& directory,ofstream& FileMESSAGE,ostream& oss) { //CO20200404 - added ofstream
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    xstructure xstr_name=GetMostRelaxedStructure(directory); //CO20180626

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " most relaxed structure:" << endl << xstr_name << endl;}

    //WARNING!!!!
    if(0){  //CO20210623 - this warning is useless, there are good reasons the two will not be same
      if(aurostd::FileExist(directory+"/POSCAR.orig")) {
        xstructure xstr_poscar_orig(directory+"/POSCAR.orig", IOVASP_POSCAR);
        if(xstr_poscar_orig.atoms.size()!=xstr_name.atoms.size()) {
          //[CO20200404 - OBSOLETE]cerr << "!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!" << endl;
          //[CO20200404 - OBSOLETE]cerr << "I tell you--'POSCAR.orig' has different atoms number from 'POSCAR.bands', though I can still produce PEDOS plots for you!" << endl;
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "POSCAR.orig has different atoms number from POSCAR.bands", directory, FileMESSAGE, oss, _LOGGER_WARNING_);
        }
      }
    }

    vector<string> AtomNames;
    uint i=0;
    int j=0;

    // ********************************************************************************  
    //READ XSTRUCTURE
    if(AtomNames.empty()){
      for(i=0;i<xstr_name.species.size();i++) {
        if(!aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr_name.species[i]).empty()) {
          aurostd::VASP_PseudoPotential_CleanName_InPlace(xstr_name.species[i]);
          for(j=0;j<xstr_name.num_each_type[i];j++){AtomNames.push_back(xstr_name.species[i]);}
        }else{AtomNames.clear();break;} //go on to next extraction method
      }
    }
    // ********************************************************************************  
    //READ OUTCAR
    if(AtomNames.empty()){
      if(XHOST.vext.size()!=XHOST.vcat.size()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "XHOST.vext.size()!=XHOST.vcat.size(), aborting", _RUNTIME_ERROR_); //CO20200404
      }
      if(aurostd::FileExist(directory+"./OUTCARfile.tmp")) aurostd::RemoveFile(directory+"./OUTCARfile.tmp"); //keep for legacy
      string OUTCARfile="";
      vector<string> vstages;
      aurostd::string2tokens(".bands,.static,.relax,.relax2,.relax1,.orig",vstages,",");
      vstages.insert(vstages.begin(),"");  //look for just OUTCAR
      uint istage=0,iext=0;
      for(istage=0;istage<vstages.size()&&OUTCARfile.empty();istage++){
        if(aurostd::FileExist(directory+"/OUTCAR"+vstages[istage])){OUTCARfile=directory+"/OUTCAR"+vstages[istage];}
      }
      string file_outcar_tmp="";
      bool delete_OUTCARfile=false;
      for(istage=0;istage<vstages.size()&&OUTCARfile.empty();istage++){
        for(iext=0;iext<XHOST.vext.size()&&OUTCARfile.empty();iext++){
          if(aurostd::FileExist(directory+"/OUTCAR"+vstages[istage]+XHOST.vext[iext])){
            aurostd::efile2tempfile(directory+"/OUTCAR"+vstages[istage]+XHOST.vext[iext],OUTCARfile,delete_OUTCARfile);
          }
        }
      }
      if(OUTCARfile.empty()){
        if(delete_OUTCARfile){aurostd::RemoveFile(OUTCARfile);}
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No OUTCAR[.static][.EXT] found in the directory, aborting", _INPUT_MISSING_);  //CO20200404
      }

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " OUTCARfile=" << OUTCARfile << endl;}

      vector<string> vpotcars;
      if(0){  //avoid making system call if we can
        string str_grep=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(aurostd::execute2string("grep TITEL "+OUTCARfile+" 2>/dev/null"));
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " grep output [1]=\"" << str_grep << "\"" << endl;}
        if(str_grep.empty()){
          str_grep=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(aurostd::execute2string("grep --text TITEL "+OUTCARfile+" 2>/dev/null"));  //CO20210623 - sometimes grep finds binary characters inside, --text will force grep to read file as text and return results
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " grep output [2]=\"" << str_grep << "\"" << endl;}
        }
        aurostd::string2tokens(str_grep,vpotcars,"\n");
      }

      vpotcars=aurostd::GrepFile(OUTCARfile,"TITEL");
      if(delete_OUTCARfile){aurostd::RemoveFile(OUTCARfile);}

      vector<string> AtomSpeciesName,vpp;
      for(i=0;i<vpotcars.size();i++){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found \"TITEL\" line=\"" << vpotcars[i] << "\"" << endl;}
        aurostd::string2tokens(vpotcars[i],vpp," ");
        if(vpp.size()<4){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,OUTCARfile+" has weird \"TITEL\" line: \""+vpotcars[i]+"\"",_FILE_CORRUPT_);  //CO20200404
        }
        AtomSpeciesName.push_back(vpp[3]);
        aurostd::VASP_PseudoPotential_CleanName_InPlace(AtomSpeciesName.back());
        for(j=0;j<xstr_name.num_each_type[i];j++){AtomNames.push_back(AtomSpeciesName.back());}
      }
    }
    // ********************************************************************************  
    //ERROR
    if(AtomNames.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AtomNames cannot be extracted",_RUNTIME_ERROR_);}  //CO20200404

    for(i=0;i<AtomNames.size();i++){aurostd::VASP_PseudoPotential_CleanName_InPlace(AtomNames[i]);} //repetita iuvant
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " AtomNames=" << aurostd::joinWDelimiter(AtomNames,",") << endl;}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}

    return AtomNames;
  }
} // namespace KBIN

// ***************************************************************************//
// KBIN::ExtractEfermiOUTCAR
// ***************************************************************************//
namespace KBIN {
  double ExtractEfermiOUTCAR(string directory) {
    double Efermi=0;
    string OUTCARfile, stmp, line;
    stringstream  strline;

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "XHOST.vext.size()!=XHOST.vcat.size(), aborting", _RUNTIME_ERROR_); //CO20200404
    }

    bool found_OUTCAR=FALSE;

    //READ OUTCAR, get the number of IONS, Fermi Level	
    string file_outcar_tmp=aurostd::TmpFileCreate("OUTCAR.tmp");
    if(aurostd::FileExist(directory+"./OUTCARfile.tmp")) aurostd::RemoveFile(directory+"./OUTCARfile.tmp");
    if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR")) {
      found_OUTCAR=TRUE;
      OUTCARfile=directory+"/OUTCAR";
    }
    if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static")) {
      found_OUTCAR=TRUE;
      OUTCARfile=directory+"/OUTCAR.static";
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static"+XHOST.vext.at(iext))) {
        found_OUTCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/OUTCAR.static"+XHOST.vext.at(iext)+" > "+file_outcar_tmp);
        OUTCARfile=file_outcar_tmp;
      }
    }
    if(!found_OUTCAR) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No OUTCAR[.static][.EXT] found in the directory, aborting", _INPUT_MISSING_);  //CO20200404
    }

    ////GET Number of IONS and Fermi
    string anchor_word_Efermi="E-fermi";
    vector<string> vlines;
    aurostd::file2vectorstring(OUTCARfile, vlines);
    for(unsigned int i=0; i<vlines.size(); i++) {
      if(vlines.at(i).find(anchor_word_Efermi) !=string::npos) {
        strline.clear();
        strline.str(vlines.at(i));
        strline >> stmp >> stmp >> Efermi;
      }
    }
    if(aurostd::FileExist(file_outcar_tmp)) aurostd::RemoveFile(file_outcar_tmp);
    //cout << Efermi << "fermi" << endl;
    return Efermi;
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::ExtractSystemName(string directory)
// ***************************************************************************
namespace KBIN {
  //ME20200217 - the old ExtractSystemName function (now ExtractSystemNameVASP)
  // breaks for long system names since VASP cuts them off. The system in the
  // aflow.in file should be the canonical system name since all VASP output
  // ultimately depends on it.
  string ExtractSystemName(const string& directory) {
    string system_name = ExtractSystemNameFromAFLOWIN(directory);
    if (system_name.empty()) {
      string message = "Could not extract system.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    aurostd::StringSubst(system_name,":.","."); //CO20200731 - patching bug in system name generation
    aurostd::StringSubst(system_name,"@ARUN.A",":ARUN.A");  //CO20200731 - patching for AEL/AGL ARUNs - @ in filenames are NOT good
    return system_name;
  }

  string ExtractSystemNameFromAFLOWIN(const string& _directory) {
    //CO20200731 - adjusting for old aflow.ins that only contain 'SYSTEM=' in INCAR section, or
    //[AFLOW] SYSTEM = Mo_pvRu_pv.A30B0
    string directory=_directory;  //CO20200624
    if (directory.empty()) directory = ".";
    string system_name = "";
    string aflowin_path = directory + "/" + _AFLOWIN_;
    if(!aurostd::FileExist(aflowin_path)) {
      string message = "Could not find file " + aflowin_path + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    // ME20200525 - Need to take potential white space between [AFLOW],
    // SYSTEM, and = into account.
    string AflowIn=aurostd::RemoveComments(aurostd::file2string(aflowin_path)); // for safety //CO20180502
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
    uint nlines = vAflowIn.size();
    for (uint iline = 0; iline < nlines && system_name.empty(); iline++) {
      if (vAflowIn[iline].find("SYSTEM") != string::npos) {
        const string& line = vAflowIn[iline];
        string::size_type t = line.find_first_of("=");
        system_name = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line.substr(t + 1));
        if(!system_name.empty()){return system_name;}
      }
    }

    return system_name;
  }

  string ExtractSystemNameFromVASP(const string& _directory) {
    bool LDEBUG=(FALSE || _DEBUG_IVASP_ || XHOST.DEBUG);
    string directory="", SystemName="", stmp="", DOSCARfile="";
    stringstream strline;

    directory=_directory; //Get directory
    if(directory=="") directory="."; // here

    if(LDEBUG) {cerr << XPID << "KBIN::ExtractSystemName: directory=" << directory << endl;}

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"XHOST.vext.size()!=XHOST.vcat.size(), aborting",_RUNTIME_ERROR_); //CO20200404
    }

    bool found=FALSE;

    //READ DOSCAR.static; Firstly, we try to get SystemName from DOSCAR
    if(LDEBUG) {cerr << XPID << "KBIN::ExtractSystemName: [1]" << endl;}
    string doscarfile_tmp=aurostd::TmpFileCreate("DOSCAR.tmp");
    if(aurostd::FileExist(directory+"/DOSCARfile.tmp")) aurostd::RemoveFile(directory+"/DOSCARfile.tmp");
    if(!found && aurostd::FileExist(directory+"/DOSCAR")) {
      found=TRUE;
      DOSCARfile=directory+"/DOSCAR";
    }
    if(!found && aurostd::FileExist(directory+"/DOSCAR.static")) {
      found=TRUE;
      DOSCARfile=directory+"/DOSCAR.static";
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/DOSCAR.static"+XHOST.vext.at(iext))) {
        found=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" \""+directory+"/DOSCAR.static"+XHOST.vext.at(iext)+"\""+" > "+doscarfile_tmp);
        DOSCARfile=doscarfile_tmp;found=TRUE;
      } 
    }    
    if(found) {
      //OPEN DOSCAR file  
      //
      vector<string> vlines;
      aurostd::file2vectorstring(DOSCARfile, vlines);
      strline.clear();
      strline.str(vlines.at(4));
      strline >> SystemName;
      //cout <<SystemName << " Is this a systemname?" << endl;
      if(aurostd::FileExist(doscarfile_tmp)) aurostd::RemoveFile(doscarfile_tmp);

      //Get the working directory
      //[CO20191112 - OBSOLETE]const int PATH_LENGTH_MAX=1024;
      //[CO20191112 - OBSOLETE]char work_dir[PATH_LENGTH_MAX];
      //[CO20191112 - OBSOLETE]getcwd(work_dir, PATH_LENGTH_MAX);  
      string work_dir=aurostd::getPWD();  //CO20191112

      //Get the data directory
      chdir(directory.c_str());  //Jump into the data directory
      //[CO20191112 - OBSOLETE]char data_dir[PATH_LENGTH_MAX];
      //[CO20191112 - OBSOLETE]getcwd(data_dir, PATH_LENGTH_MAX);  
      string data_dir=aurostd::getPWD();  //CO20191112
      chdir(work_dir.c_str()); //Jump back into the work directory  //CO20191112

      //Check whether it is MAGNETIC DATABASE
      bool FLAG_MAGNETIC=FALSE;
      string abs_dir(data_dir);  //convert the char type of data_dir into string type
      if(abs_dir.find("LIB3")!=string::npos) FLAG_MAGNETIC=TRUE;

      if(FLAG_MAGNETIC) {
        return SystemName;
      } else {
        vector<string> data;
        aurostd::string2tokens(data_dir.c_str(), data, "/");  //CO20191112
        int datasize=data.size();
        bool FLAG_DIRECTORY_ICSD=FALSE;
        for (uint i=0; i<data.size();i++) {
          if((data.at(i).find("ICSD")!=string::npos)) {
            FLAG_DIRECTORY_ICSD =TRUE;
          }
        }

        //Check whether the last entry is "BANDS"
        string Last_Entry=data.back();
        bool FLAG_LAST_ENTRY_BANDS=FALSE;
        if(Last_Entry.compare("BANDS")==0) FLAG_LAST_ENTRY_BANDS=TRUE;

        //Check whether the SystemName is "unknown", to fix the bug of aflow
        bool FLAG_SYSTEMNAME_UNKNOWN=FALSE;
        if(SystemName.compare("unknown")==0) FLAG_SYSTEMNAME_UNKNOWN=TRUE;
        //if ICSD is not found in system name
        if((SystemName.find("ICSD")==string::npos)) {
          if(datasize >=2) {
            if(FLAG_LAST_ENTRY_BANDS) {
              if(FLAG_SYSTEMNAME_UNKNOWN && FLAG_DIRECTORY_ICSD) {
                SystemName=data.at(datasize-2);
              } else SystemName=data.at(datasize-3)+"."+data.at(datasize-2);
            } else {
              if(FLAG_SYSTEMNAME_UNKNOWN && FLAG_DIRECTORY_ICSD) {
                SystemName=Last_Entry;
              } else SystemName=data.at(datasize-2)+"."+data.back();
            }
          } else {
            SystemName=data.back();
          }
        }
        return SystemName;
      }
    }

    // TRY OUTCAR.relax2
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/OUTCAR.relax2"+XHOST.vext.at(iext))) {
        found=TRUE;
        xOUTCAR outcar(directory+"/OUTCAR.relax2"+XHOST.vext.at(iext));
        SystemName=outcar.SYSTEM;
        return SystemName;
      }
    }
    // TRY OUTCAR.relax1
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/OUTCAR.relax1"+XHOST.vext.at(iext))) {
        found=TRUE;
        xOUTCAR outcar(directory+"/OUTCAR.relax1"+XHOST.vext.at(iext));
        SystemName=outcar.SYSTEM;
        return SystemName;
      }
    }
    // NOTHING FOUND !!!
    if(!found) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No DOSCAR[.static][.EXT], OUTCAR[.relax1|.relax2]][.EXT] found in the directory, aborting",_FILE_CORRUPT_); //CO20200624
    }
    return "";
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::ExtractPOSCARFromAFLOWIN
// ***************************************************************************
namespace KBIN {
  // SD20220228 - Extract the nth POSCAR from the AFLOWIN, negative values go backwards
  bool ExtractPOSCARStringStreamFromDirectory(const string& directory, stringstream& poscar, const int index) {
    string AflowIn = "";
    if(!aurostd::file2string(directory + "/" + _AFLOWIN_, AflowIn)) {return false;}
    return ExtractPOSCARStringStreamFromAFLOWIN(AflowIn, poscar, index);
  }

  xstructure GetPOSCARXStructureFromDirectory(const string& directory, const int iomode, const int index) {
    stringstream poscar;
    if(!ExtractPOSCARStringStreamFromDirectory(directory, poscar, index)) {
      string message = "Could not find valid " + _AFLOWIN_ + " in " + directory;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    return xstructure(poscar, iomode);
  }

  bool ExtractPOSCARStringStreamFromAFLOWIN(const string& AflowIn, stringstream& poscar, const int index) {
    if(index==0) {
      string message = "Index cannot be 0";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    poscar.str(aurostd::substring2string(AflowIn, _VASP_POSCAR_MODE_EXPLICIT_START_, _VASP_POSCAR_MODE_EXPLICIT_STOP_, index));
    if(poscar.str().empty()) {return false;}
    return true;
  }

  xstructure GetPOSCARXStructureFromAFLOWIN(const string& AflowIn, const int iomode, const int index) {
    stringstream poscar;
    if(!ExtractPOSCARStringStreamFromAFLOWIN(AflowIn, poscar, index)) {
      string message = "Invalid " + _AFLOWIN_;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
    return xstructure(poscar, iomode);
  }

} // namespace KBIN

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
