// ***************************************************************************
// *                                                                         *
// *              AFlow STEFANO CURTAROLO  Duke University 2003-2023         *
// *              AFlow COREY OSES  Duke University 2013-2021                *
// *                                                                         *
// ***************************************************************************
// aflow_pocc.cpp
// completely revised approach to KESONG YANG's original implementation
// no recursion needed --- too memory intensive
// issues with UFF (structure comparison), multiply occupied sites, and 
// vacancies are all addressed robustly here
// 2013-2019: corey.oses@duke.edu
// 2010-2011: kesong.yang@gmail.com (LEGACY)

// This file contains the routines to prepare partial occupation input files.

#ifndef _AFLOW_POCC_CPP_
#define _AFLOW_POCC_CPP_

#include "aflow.h"
#include "aflow_pocc.h"
#include "aflow_compare_structure.h"

#define _DEBUG_POCC_ false  //CO20190116
#define _DEBUG_POCC_CLUSTER_ANALYSIS_ false && _DEBUG_POCC_  //CO20190116

const string POSCAR_START_tag=_VASP_POSCAR_MODE_EXPLICIT_START_; //no-period is important
const string POSCAR_STOP_tag=_VASP_POSCAR_MODE_EXPLICIT_STOP_; //no-period is important
const string POSCAR_series_START_tag=_VASP_POSCAR_MODE_EXPLICIT_START_P_; //period is important
const string POSCAR_series_STOP_tag=_VASP_POSCAR_MODE_EXPLICIT_STOP_P_; //period is important
const string POSCAR_POCC_series_START_tag=POSCAR_series_START_tag+"POCC_";
const string POSCAR_POCC_series_STOP_tag=POSCAR_series_STOP_tag+"POCC_";

const string POCC_AFLOWIN_tag="[AFLOW_POCC]";

//make defaults in AFLOW_RC
const double ENERGY_RADIUS = 10; //angstroms  //keep, so we can still compare with KY

//some constants
const int A_START = 1, C_START = 1, F_START = 1;
const int B_START = 0, D_START = 0, E_START = 0;

bool PRIMITIVIZE=false;

//testing
bool COMPARE_WITH_KESONG=false;
bool ENUMERATE_ALL_HNF=false;

namespace pocc {
  //temporary function to replace KY's code, but should never be called
  //does NOT handle aflow.in generation well at all
  bool poccInput() {
    //aflow_pocc.log
    ofstream FileMESSAGE; FileMESSAGE.open("aflow_pocc.log");
    //directory
    string directory="";
    if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");} //CO20190402
    if(directory.empty() || directory=="./" || directory=="."){directory=aurostd::getPWD()+"/";} //".";  //CO20180220 //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
    KBIN::VASP_RunPOCC(directory,FileMESSAGE);
    return true;
  }
} // namespace pocc

namespace pocc {
  void poccOld2New(ostream& oss){ofstream FileMESSAGE;return poccOld2New(FileMESSAGE,oss);}
  void poccOld2New(ofstream& FileMESSAGE,ostream& oss){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    //aflags
    _aflags aflags;
    if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");} //CO20190402
    if(aflags.Directory.empty() || aflags.Directory=="./" || aflags.Directory=="."){aflags.Directory=aurostd::getPWD()+"/";} //".";  //CO20180220 //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))

    //aflow.in
    string AflowIn_file="",AflowIn="";
    try{KBIN::getAflowInFromAFlags(aflags,AflowIn_file,AflowIn,FileMESSAGE,oss);}
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return;
    }

    //other flags
    _kflags kflags=KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,oss);
    _vflags vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags,oss);

    if(!kflags.KBIN_POCC){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Not a POCC run",_INPUT_MISSING_);}

    if(!vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No POCC ordered structures found in "+_AFLOWIN_,_INPUT_MISSING_);}

    if(LDEBUG){
      for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++){
        cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i=" << i << "]=" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i] << endl;
        cerr << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE[i] << endl;
      }
    }

    message << "Writing " << POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    stringstream unique_structures_ss;
    for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++){
      unique_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
      unique_structures_ss << _VASP_POSCAR_MODE_EXPLICIT_START_P_ << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i] << endl;
      unique_structures_ss << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE[i];
      unique_structures_ss << _VASP_POSCAR_MODE_EXPLICIT_STOP_P_ << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i] << endl;
      unique_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
    }
    aurostd::stringstream2file(unique_structures_ss,aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE);

    vector<string> vlines_orig,vlines;
    aurostd::file2vectorstring(AflowIn_file,vlines_orig); //keep comments
    bool reading_explicit=false;
    uint iline=0;
    for(iline=0;iline<vlines_orig.size();iline++){
      const string& line=vlines_orig[iline];
      if(line.find(_VASP_POSCAR_MODE_EXPLICIT_STOP_P_)!=string::npos){reading_explicit=false;}
      else if(line.find(_VASP_POSCAR_MODE_EXPLICIT_START_P_)!=string::npos){reading_explicit=true;}
      else{
        if(reading_explicit==false){vlines.push_back(line);}
      }
    }
    string search="[AFLOW] **********";  //in case AFLOWIN_SEPARATION_LINE changes in length
    for(iline=vlines.size()-1;iline<vlines.size();iline--){
      if(vlines[iline].find(search)==string::npos){break;} //separation line not found
      if((iline-1)<vlines.size() && vlines[iline-1].find(search)!=string::npos){vlines.pop_back();} //double separation line found, remove last one
    }
    string AflowIn_file_ORIG=AflowIn_file;aurostd::StringSubst(AflowIn_file_ORIG,".in",".ORIG.in");
    aurostd::file2file(AflowIn_file,AflowIn_file_ORIG);
    message << "Saving " << AflowIn_file << " as " << AflowIn_file_ORIG;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    message << "Rewriting " << AflowIn_file;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
    stringstream aflowin_ss;for(uint i=0;i<vlines.size();i++){aflowin_ss << vlines[i] << endl;}
    aurostd::stringstream2file(aflowin_ss,AflowIn_file);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}

    return;
  }
} // namespace pocc

namespace pocc {
  void POCC_Convolution(const aurostd::xoption& vpflow){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    //CO+ME20220118
    //this code needs to be revisited
    //checking for XHOST.DIRECTORY would make command line syntax "aflow ... -D DIR"
    //can this take multiple inputs?
    //would a --conv_dir=... flag be better?
    //address for lib2raw
    _aflags aflags;aflags.Directory=".";  //make current directory for now, we can change in the future

    ofstream FileMESSAGE;
    string FileNameLOCK=aflags.Directory+"/"+_AFLOWLOCK_;
    FileMESSAGE.open(FileNameLOCK.c_str(),std::ios::out);
    POccCalculator pcalc(vpflow,FileMESSAGE);
    pcalc.m_aflags.Directory=aflags.Directory;
    pcalc.CleanPostProcessing();
    pcalc.convolution();
    FileMESSAGE.flush();FileMESSAGE.clear();FileMESSAGE.close();
    stringstream command;
    KBIN::CompressDirectory(aflags);
    aurostd::execute(command);
  }
} // namespace pocc

namespace pocc {
  string addDefaultPOCCTOL2string(const string& input){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    string output=input;
    if(output.find(TAG_TITLE_POCC_TOL)==string::npos){
      int prec=3;
      prec=(int)ceil(log10(1.0/DEFAULT_POCC_SITE_TOL));
      output+=TAG_TITLE_POCC_TOL+aurostd::utype2string(DEFAULT_POCC_SITE_TOL,prec);
      prec=(int)ceil(log10(1.0/DEFAULT_POCC_STOICH_TOL));
      output+="_"+aurostd::utype2string(DEFAULT_POCC_STOICH_TOL,prec);
    }
    return output;
  }
}

namespace pocc {
  double poccDOSCAR2temperature(const string& doscar_path){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    //get temperature from title
    //DOSCAR.pocc_T0000K
    //ME20210927 - Added PHDOSCAR support
    string pocc_doscar_start = "";
    if(doscar_path.find(POCC_DOSCAR_PREFIX) != string::npos){
      pocc_doscar_start = POCC_DOSCAR_PREFIX;
    } else if (doscar_path.find(POCC_PHDOSCAR_PREFIX) != string::npos) {
      pocc_doscar_start = POCC_PHDOSCAR_PREFIX;
    } else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"odd DOSCAR filename, format should be "+POCC_DOSCAR_PREFIX+"0000K",_FILE_CORRUPT_);
    }
    vector<string> vtokens;
    aurostd::string2tokens(doscar_path,vtokens,"/");
    string pocc_doscar_end="K";
    string::size_type loc_start=vtokens.back().find(pocc_doscar_start);
    string::size_type loc_end=vtokens.back().find(pocc_doscar_end);
    if(loc_start==string::npos||loc_end==string::npos){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"cannot get temperature_str",_FILE_CORRUPT_);
    }
    string temperature_str=vtokens.back().substr(loc_start+(pocc_doscar_start.size()),loc_end-(loc_start+(pocc_doscar_start.size())));
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " temperature_str=" << temperature_str << endl;}
    if(!aurostd::isfloat(temperature_str)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"temperature_str is not a float",_FILE_CORRUPT_);}
    //[CO+ME20200921 - just look at vtokens.back()]uint i=0,j=0;
    //[CO+ME20200921 - just look at vtokens.back()]for(i=0;i<vtokens.size();i++){
    //[CO+ME20200921 - just look at vtokens.back()]  if(vtokens[i].find("pocc_T")!=string::npos && vtokens[i].find("K")!=string::npos){
    //[CO+ME20200921 - just look at vtokens.back()]    temperature_str=vtokens[i];
    //[CO+ME20200921 - just look at vtokens.back()]    aurostd::StringSubst(temperature_str,POCC_DOSCAR_PREFIX,"");
    //[CO+ME20200921 - just look at vtokens.back()]    aurostd::StringSubst(temperature_str,"K","");
    //[CO+ME20200921 - just look at vtokens.back()]    for(j=0;j<XHOST.vext.size();j++){aurostd::StringSubst(temperature_str,XHOST.vext[j],"");} //remove compression extension
    //[CO+ME20200921 - just look at vtokens.back()]    if(LDEBUG){cerr << __AFLOW_FUNC__ << " temperature_str=" << temperature_str << endl;}
    //[CO+ME20200921 - just look at vtokens.back()]    if(aurostd::isfloat(temperature_str)){break;}
    //[CO+ME20200921 - just look at vtokens.back()]  }
    //[CO+ME20200921 - just look at vtokens.back()]}
    //[CO+ME20200921 - just look at vtokens.back()]if(doscar_path.find(POCC_DOSCAR_PREFIX)==string::npos){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"cannot get temperature_str",_FILE_CORRUPT_);}
    double temperature=aurostd::string2utype<double>(temperature_str);

    return temperature;
  }
}

namespace KBIN {
  void VASP_RunPOCC(const string& directory,ostream& oss) {ofstream FileMESSAGE;return VASP_RunPOCC(directory,FileMESSAGE,oss);}
  void VASP_RunPOCC(const string& directory,ofstream& FileMESSAGE,ostream& oss) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    _xvasp xvasp;
    string AflowIn_file="",AflowIn="";
    _aflags aflags;aflags.Directory=directory;

    //get AflowIn
    try{KBIN::getAflowInFromAFlags(aflags,AflowIn_file,AflowIn,FileMESSAGE,oss);}
    catch(aurostd::xerror& err){
      pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      FileMESSAGE.close();
      return;
    }

    _kflags kflags=KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn,aflags);
    _vflags vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,aflags,kflags);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " kflags.KBIN_PHONONS_CALCULATION_AEL=" << kflags.KBIN_PHONONS_CALCULATION_AEL << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " kflags.KBIN_PHONONS_CALCULATION_AGL=" << kflags.KBIN_PHONONS_CALCULATION_AGL << endl;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " kflags.KBIN_PHONONS_CALCULATION_APL=" << kflags.KBIN_PHONONS_CALCULATION_APL << endl;}

    return VASP_RunPOCC(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE,oss);
  }
  void VASP_RunPOCC(const _xvasp& xvasp,const string& AflowIn,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    stringstream message;

    try{
      xstructure xstr_pocc=pocc::extractPARTCAR(AflowIn); //prefer to pull from AflowIn input vs. aflags
      pocc::POccCalculator pcalc(xstr_pocc,aflags,kflags,vflags,FileMESSAGE,oss);

      if(!pocc::structuresGenerated(aflags.Directory)){
        try{pcalc.generateStructures(xvasp);}
        catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);}
        //ME2021104
        if (kflags.KBIN_PHONONS_CALCULATION_AEL
            || kflags.KBIN_PHONONS_CALCULATION_AGL
            || kflags.KBIN_PHONONS_CALCULATION_APL) {
          string message = "Cannot run AFLOW modules without fully completed POCC calculations."
            " Please see the README for more information.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, aflags, FileMESSAGE, oss, _LOGGER_NOTICE_);
        }
        return;
      }

      //ME20211004
      if (kflags.KBIN_PHONONS_CALCULATION_APL) {
        // Simple check: if APL input files are not present anywhere,
        // assume that we need to create aflow.in files
        if (!pcalc.inputFilesFoundAnywhereAPL()) {
          pcalc.createModuleAflowIns(xvasp, "APL");
          return;
        }
      }

      //post-processing
      pcalc.CleanPostProcessing();
      //pocc::patchStructuresFile(aflags,FileMESSAGE,oss);  //patch if needed
      pcalc.postProcessing();
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);return;}

    return;
  }
} // namespace pocc

namespace pocc {
  string POCC_MINIMUM_CONFIGURATION(const aurostd::xoption& vpflow){return POCC_MINIMUM_CONFIGURATION(vpflow.getattachedscheme("POCC_MINIMUM_CONFIGURATION"));}
  string POCC_MINIMUM_CONFIGURATION(const string& directory){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(!aurostd::IsDirectory(directory)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Directory="+directory+" does not exist",_FILE_NOT_FOUND_);}
    vector<string> ls_contents;
    aurostd::DirectoryLS(directory,ls_contents);
    vector<xQMVASP> v_qmvasp;
    for(uint i=0;i<ls_contents.size();i++){
      if(aurostd::IsDirectory(ls_contents[i]) && aurostd::substring2bool(ls_contents[i],"ARUN.POCC")){
        if(!(aurostd::EFileExist(ls_contents[i]+"/"+DEFAULT_AFLOW_QMVASP_OUT))){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,ls_contents[i]+" does not contain "+DEFAULT_AFLOW_QMVASP_OUT,_FILE_NOT_FOUND_);}
        v_qmvasp.push_back(xQMVASP(ls_contents[i]+"/"+DEFAULT_AFLOW_QMVASP_OUT));
        if(v_qmvasp.back().H_atom_relax==AUROSTD_NAN){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,ls_contents[i]+"/"+DEFAULT_AFLOW_QMVASP_OUT+" does not show a relaxation run",_FILE_ERROR_);}
      }
    }
    uint v_qmvasp_size=v_qmvasp.size(); //NO MORE PUSH_BACK BELOW! Will save time for AAPL calculations
    if(v_qmvasp_size<1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No ARUN.POCC found in Directory="+directory,_FILE_NOT_FOUND_);}
    //sort by H_atom_relax
    xQMVASP qmvasp_tmp;
    for(uint i=0;i<v_qmvasp_size-1;i++){
      for(uint j=i;j<v_qmvasp_size;j++){
        if(v_qmvasp[i].H_atom_relax>v_qmvasp[j].H_atom_relax){
          qmvasp_tmp=v_qmvasp[i];
          v_qmvasp[i]=v_qmvasp[j];
          v_qmvasp[j]=qmvasp_tmp;
        }
      }
    }
    string qmvasp_filename="";
    if(LDEBUG){
      for(uint i=0;i<v_qmvasp_size-1;i++){
        qmvasp_filename=v_qmvasp[i].filename;
        aurostd::StringSubst(qmvasp_filename,"/"+DEFAULT_AFLOW_QMVASP_OUT,"");
        cerr << __AFLOW_FUNC__ << " H_atom_relax(" << qmvasp_filename << ")=" << v_qmvasp[i].H_atom_relax << endl;
      }
    }
    qmvasp_filename=v_qmvasp[0].filename;
    aurostd::StringSubst(qmvasp_filename,"/"+DEFAULT_AFLOW_QMVASP_OUT,"");
    return qmvasp_filename;
  }
} // namespace pocc

namespace pocc {
  bool structuresGenerated(const string& directory){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    string file="";
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " aurostd::EFileExist(" << directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE << ")=" << aurostd::EFileExist(directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,file) << endl;}
    if(!aurostd::EFileExist(directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,file)){return false;} //CO20200606 - necessary because efile2tempfile is verbose
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " aurostd::EFileNotEmpty(" << file << ")=" << aurostd::EFileNotEmpty(file) << endl;}
    if(aurostd::EFileNotEmpty(file)){return true;}
    return false;
  }
  xstructure extractPARTCAR(const string& AflowIn){
    stringstream ss_pocc_structure;
    ss_pocc_structure.str(aurostd::substring2string(AflowIn,"[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE","[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE",-1));
    //[SD20220520 - OBSOLETE]aurostd::ExtractLastToStringstreamEXPLICIT(AflowIn,ss_pocc_structure, "[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE", "[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE");
    if(ss_pocc_structure.str().empty()){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No PARTCAR found (looking for [POCC_MODE_EXPLICIT]START.POCC_STRUCTURE / [POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE)");
    }
    xstructure xstr_pocc(ss_pocc_structure, IOVASP_POSCAR);
    return xstr_pocc;
  }
} // namespace pocc

namespace pocc {
  void parsePOccHashFromXStructureTitle(const string& title,string& pocc_hash){
    string hnf_index_str="",site_config_index_str="";
    return parsePOccHashFromXStructureTitle(title,pocc_hash,hnf_index_str,site_config_index_str);
  }
  void parsePOccHashFromXStructureTitle(const string& title,string& pocc_hash,string& hnf_index_str,string& site_config_index_str){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " title=" << title << endl;}
    vector<string> vtokens,vtokens2;
    string tmp_str;
    aurostd::string2tokens(title,vtokens," ");
    pocc_hash="",hnf_index_str="",site_config_index_str="";
    for(uint t=0;t<vtokens.size();t++){ //find pocc hash directly
      if(aurostd::substring2bool(vtokens[t],"H") && aurostd::substring2bool(vtokens[t],"C")){
        pocc_hash=vtokens[t];
        aurostd::string2tokens(pocc_hash,vtokens2,"C");
        if(vtokens2.size()!=2){pocc_hash="";continue;}  //did not find pocc_hash, keep going
        //check hnf_index
        tmp_str=vtokens2[0]; //should contain an H
        if(!aurostd::substring2bool(tmp_str,"H")){pocc_hash="";continue;}  //did not find pocc_hash, keep going
        aurostd::StringSubst(tmp_str,"H",""); //should now be hnf_index
        hnf_index_str=tmp_str;
        if(!aurostd::isfloat(hnf_index_str)){pocc_hash="";continue;}  //did not find pocc_hash, keep going
        //now check site_config_index
        site_config_index_str=vtokens2[1];
        if(!aurostd::isfloat(site_config_index_str)){pocc_hash="";continue;}  //did not find pocc_hash, keep going
        break;
      }
    }
    if(pocc_hash.empty()){  //otherwise build pocc hash from components in the title
      vector<string> vtokens3;
      for(uint t=0;t<vtokens.size();t++){
        if(aurostd::substring2bool(vtokens[t],"HNF(n=")){ //HNF(n=2,#2/3)
          aurostd::string2tokens(vtokens[t],vtokens2,",");
          if(vtokens2.size()!=2){continue;}    //did not find pocc_hash, keep going
          aurostd::string2tokens(vtokens2[1],vtokens3,"/");
          if(vtokens3.size()!=2){continue;}    //did not find pocc_hash, keep going
          aurostd::StringSubst(vtokens3[0],"#","");
          if(!aurostd::isfloat(vtokens3[0])){continue;}   //did not find pocc_hash, keep going
          hnf_index_str=aurostd::utype2string( aurostd::string2utype<int>(vtokens3[0])-1 );
          continue;
        }
        if(aurostd::substring2bool(vtokens[t],"site_config(")){ //site_config(#1/2)
          aurostd::string2tokens(vtokens[t],vtokens2,"/");
          if(vtokens2.size()!=2){continue;}    //did not find pocc_hash, keep going
          aurostd::StringSubst(vtokens2[0],"site_config(","");
          aurostd::StringSubst(vtokens2[0],"#","");
          if(!aurostd::isfloat(vtokens2[0])){continue;}   //did not find pocc_hash, keep going
          site_config_index_str=aurostd::utype2string( aurostd::string2utype<int>(vtokens2[0])-1 );
        }
      }
      pocc_hash="H"+hnf_index_str+"C"+site_config_index_str;
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " hnf_index_str=" << hnf_index_str << endl;
      cerr << __AFLOW_FUNC__ << " site_config_index_str=" << site_config_index_str << endl;
      cerr << __AFLOW_FUNC__ << " pocc_hash=" << pocc_hash << endl;
    }
  }
  unsigned long long int getDGFromXStructureTitle(const string& title){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    if(!aurostd::substring2bool(title,"DG=")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"\"DG=\" not found in title [title=\""+title+"\"]",_FILE_CORRUPT_);}
    unsigned long long int dg=1;
    vector<string> vtokens,vtokens2;
    aurostd::string2tokens(title,vtokens," ");
    for(uint i=0;i<vtokens.size();i++){
      if(aurostd::substring2bool(vtokens[i],"DG=")){
        aurostd::string2tokens(vtokens[i],vtokens2,"=");
        if(vtokens2.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot parse at DG= (vtokens2.size()=="+aurostd::utype2string(vtokens2.size())+")",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(vtokens2[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"degeneracy is not an integer",_VALUE_ILLEGAL_);}
        dg*=aurostd::string2utype<unsigned long long int>(vtokens2[1]); //backwards compatible with old pocc scheme, but this should be obsolete
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " DG=" << dg << endl;}
    return dg;
  }
  void parsePropertyByTag(const string& line,const string& tag,double& prop){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    vector<string> vtokens;
    aurostd::string2tokens(line,vtokens,"=");
    if(vtokens.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown "+tag+" line format",_FILE_CORRUPT_);}
    if(!aurostd::isfloat(vtokens[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,tag+" is not a double",_VALUE_ERROR_);}
    prop=aurostd::string2utype<double>(vtokens[1]);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << tag << "=" << prop << endl;}
  }
  bool patchStructuresAllFile(const _aflags& aflags,string& structures_file,stringstream& structures_file_ss,ofstream& FileMESSAGE,ostream& oss){
    bool found_file=false;
    if(!found_file && aurostd::EFileExist(aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE,structures_file)){aurostd::efile2stringstream(structures_file,structures_file_ss);found_file=true;}
    if(!found_file && aurostd::FileExist(aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE,structures_file)){aurostd::file2stringstream(structures_file,structures_file_ss);found_file=true;}
    if(!found_file){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE+" does not exist",_FILE_NOT_FOUND_);}
    if(aurostd::substring2bool(structures_file_ss,POSCAR_series_START_tag)){return false;}

    stringstream message;
    message << "Patching " << POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    return patchStructuresAllFile(structures_file_ss);
  }
  bool patchStructuresAllFile(stringstream& structures_file_ss){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    vector<string> vlines,vtokens,vtokens2;
    string pocc_hash="";
    stringstream pocc_addendum_ss;

    unsigned long long int l_supercell_sets_size=0;
    vector<unsigned long long int> pscs_psc_set_sizes;
    unsigned long long int i_l_supercell_sets_size=0,i_pscs_psc_set=0;

    string POccStructure_title_tag="HNF(n="; //in case we cannot find POSCAR_START_tag, SHACHAR OLD FORMAT
    bool add_POSCAR_START_tag=false;  //SHACHAR OLD FORMAT
    bool add_l_supercell_sets_size=false; //SHACHAR OLD FORMAT
    string structures_file_str=structures_file_ss.str();
    bool found_tag_sg=false;
    string structure_group_tag="";
    //STRUCTURES_GROUP
    if(found_tag_sg==false){
      structure_group_tag=POCC_AFLOWIN_tag+"STRUCTURES_GROUP";
      if(aurostd::substring2bool(structures_file_str,structure_group_tag)){found_tag_sg=true;}
    }
    //STRUCTURE_GROUP
    if(found_tag_sg==false){
      structure_group_tag=POCC_AFLOWIN_tag+"STRUCTURE_GROUP";
      if(aurostd::substring2bool(structures_file_str,structure_group_tag)){found_tag_sg=true;}
    }
    //STRUCTURE GROUP
    if(found_tag_sg==false){
      structure_group_tag=POCC_AFLOWIN_tag+"STRUCTURE GROUP";
      if(aurostd::substring2bool(structures_file_str,structure_group_tag)){found_tag_sg=true;}
    }
    //not found
    if(found_tag_sg==false){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find structure group tag",_FILE_CORRUPT_);}

    if(!aurostd::substring2bool(structures_file_str,POSCAR_series_START_tag)){
      if(!aurostd::substring2bool(structures_file_str,POSCAR_START_tag)){
        add_POSCAR_START_tag=true;
        if(!aurostd::substring2bool(structures_file_str,POccStructure_title_tag)){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown format for "+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE+" file: cannot find POSCAR title tag",_FILE_ERROR_);
        }
      }
      aurostd::string2vectorstring(structures_file_str,vlines); //aurostd::efile2vectorstring(structures_file,vlines);
      for(uint iline=0;iline<vlines.size();iline++){  //see if you can get l_supercell_sets_size from first structure_group_tag line
        if(aurostd::substring2bool(vlines[iline],structure_group_tag)){
          if(l_supercell_sets_size==0 && aurostd::substring2bool(vlines[iline],"/")){
            aurostd::string2tokens(vlines[iline],vtokens," ");
            aurostd::string2tokens(vtokens.back(),vtokens2,"/");
            if(vtokens2.size()!=2){break;}  //didn't work //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown format for line containing l_supercell_sets_size",_FILE_CORRUPT_);
            if(!aurostd::isfloat(vtokens2.back())){break;}  //didn't work //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,vtokens2.back()+" is not a number",_FILE_CORRUPT_);
            l_supercell_sets_size=aurostd::string2utype<int>(vtokens2.back());
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " l_supercell_sets_size=" << l_supercell_sets_size << endl;}
            break;
          }
        }
      }
      string poscar_start_tag=(add_POSCAR_START_tag ? POccStructure_title_tag : POSCAR_START_tag );
      string poscar_stop_tag=(add_POSCAR_START_tag ? "[AFLOW] **" : POSCAR_STOP_tag );
      for(uint iline=0;iline<vlines.size();iline++){
        if(aurostd::substring2bool(vlines[iline],structure_group_tag)){pscs_psc_set_sizes.push_back(0);}
        if(aurostd::substring2bool(vlines[iline],poscar_start_tag)){pscs_psc_set_sizes.back()++;}
      }
      if(l_supercell_sets_size==0){add_l_supercell_sets_size=true;l_supercell_sets_size=pscs_psc_set_sizes.size();} //could not find from structure_group_tag
      if(l_supercell_sets_size==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find l_supercell_sets_size",_RUNTIME_INIT_);}
      if(pscs_psc_set_sizes.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find pscs_psc_set_sizes",_RUNTIME_INIT_);}
      if((uint)l_supercell_sets_size!=pscs_psc_set_sizes.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"l_supercell_sets_size[="+aurostd::utype2string((uint)l_supercell_sets_size)+"]!=pscs_psc_set_sizes.size()[="+aurostd::utype2string(pscs_psc_set_sizes.size())+"]",_INDEX_MISMATCH_);}
      for(uint i=0;i<pscs_psc_set_sizes.size();i++){
        if(pscs_psc_set_sizes[i]==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"pscs_psc_set_sizes[i="+aurostd::utype2string(i)+"]==0",_RUNTIME_INIT_);}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " pscs_psc_set_sizes[i=" << i << "]=" << pscs_psc_set_sizes[i] << endl;}
      }
      i_l_supercell_sets_size=-1;i_pscs_psc_set=0;
      bool found_POccStructure=false;
      for(uint iline=0;iline<vlines.size();iline++){
        if(aurostd::substring2bool(vlines[iline],structure_group_tag)){ //reset
          i_l_supercell_sets_size++;
          i_pscs_psc_set=0;
          if(add_l_supercell_sets_size){
            aurostd::StringSubst(vlines[iline],POCC_AFLOWIN_tag+"STRUCTURE_GROUP",POCC_AFLOWIN_tag+"STRUCTURES_GROUP");
            aurostd::StringSubst(vlines[iline],POCC_AFLOWIN_tag+"STRUCTURE GROUP",POCC_AFLOWIN_tag+"STRUCTURES_GROUP");
            aurostd::StringSubst(vlines[iline],aurostd::RemoveWhiteSpacesFromTheBack(vlines[iline]),aurostd::RemoveWhiteSpacesFromTheBack(vlines[iline])+"/"+aurostd::utype2string(l_supercell_sets_size));
          }
        }
        if(found_POccStructure==true && aurostd::substring2bool(vlines[iline],poscar_stop_tag)){
          if(pocc_addendum_ss.str().empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"STOP tag found before START",_FILE_CORRUPT_);}
          if(aurostd::substring2bool(vlines[iline],POSCAR_STOP_tag)){aurostd::StringSubst(vlines[iline],POSCAR_STOP_tag,POSCAR_series_STOP_tag+pocc_addendum_ss.str());}
          else{vlines.insert(vlines.begin()+iline++,POSCAR_series_STOP_tag+pocc_addendum_ss.str());}
          pocc_addendum_ss.str("");
          found_POccStructure=false;
        }
        if(aurostd::substring2bool(vlines[iline],poscar_start_tag)){
          found_POccStructure=true;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " vlines[iline=" << iline << "]=\"" << vlines[iline] << "\"" << endl;}
          if(i_l_supercell_sets_size>=l_supercell_sets_size){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found more unique supercells than expected",_INDEX_MISMATCH_);}
          if(aurostd::substring2bool(vlines[iline],POSCAR_START_tag)){
            if(iline+1>vlines.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot grab title at iline="+aurostd::utype2string(iline),_FILE_CORRUPT_);} //see if title is accessible
            parsePOccHashFromXStructureTitle(vlines[iline+1],pocc_hash);
          }else{parsePOccHashFromXStructureTitle(vlines[iline],pocc_hash);} //SHACHAR OLD FORMAT
          if(pocc_hash.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not find pocc_hash at iline="+aurostd::utype2string(iline+1),_FILE_CORRUPT_);}
          pocc_addendum_ss << "POCC_";
          pocc_addendum_ss << std::setfill('0') << std::setw(aurostd::getZeroPadding(l_supercell_sets_size)) << i_l_supercell_sets_size+1 << "_";
          if((uint)i_l_supercell_sets_size>pscs_psc_set_sizes.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"i_l_supercell_sets_size>pscs_psc_set_sizes.size()-1",_INDEX_MISMATCH_);}
          pocc_addendum_ss << std::setfill('0') << std::setw(aurostd::getZeroPadding(pscs_psc_set_sizes[(uint)i_l_supercell_sets_size])) << (i_pscs_psc_set++)+1 << "_";
          pocc_addendum_ss << pocc_hash;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_addendum=" << pocc_addendum_ss.str() << endl;}
          //everything checks out, patch the line
          if(aurostd::substring2bool(vlines[iline],POSCAR_START_tag)){aurostd::StringSubst(vlines[iline],POSCAR_START_tag,POSCAR_series_START_tag+pocc_addendum_ss.str());}
          else{vlines.insert(vlines.begin()+iline++,POSCAR_series_START_tag+pocc_addendum_ss.str());} //++ so we don't find the same line
          if(!aurostd::substring2bool(vlines[iline],pocc_hash)){aurostd::StringSubst(vlines[iline],POccStructure_title_tag,pocc_hash+" "+POccStructure_title_tag);} //SHACHAR OLD FORMAT
          if(!aurostd::substring2bool(vlines[iline],DOI_POCC)){vlines[iline]=aurostd::RemoveWhiteSpacesFromTheBack(vlines[iline])+DOI_POCC;}  //SHACHAR OLD FORMAT
        }
        aurostd::StringSubst(vlines[iline],"#",""); //NO continue
      }
      //stringstream structures_file_ss_new;
      structures_file_ss.str("");
      for(uint iline=0;iline<vlines.size();iline++){structures_file_ss << vlines[iline] << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " patched file=" << endl;cerr << structures_file_ss.str() << endl;}
      return true;
    }
    return false;
  }

  bool patchStructuresUniqueFile(const _aflags& aflags,string& structures_file,stringstream& structures_file_ss,ofstream& FileMESSAGE,ostream& oss){
    bool found_file=false;
    if(!found_file && aurostd::EFileExist(aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,structures_file)){aurostd::efile2stringstream(structures_file,structures_file_ss);found_file=true;}
    if(!found_file && aurostd::FileExist(aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,structures_file)){aurostd::file2stringstream(structures_file,structures_file_ss);found_file=true;}
    if(!found_file){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE+" does not exist",_FILE_NOT_FOUND_);}
    if(aurostd::substring2bool(structures_file_ss,POSCAR_series_START_tag)){return false;}

    stringstream message;
    message << "Patching " << POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    return patchStructuresUniqueFile(structures_file_ss);
  }
  bool patchStructuresUniqueFile(stringstream& structures_file_ss){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    vector<string> vlines;
    string pocc_hash="";
    stringstream pocc_addendum_ss;

    unsigned long long int l_supercell_sets_size=0;
    unsigned long long int i_l_supercell_sets_size=0;

    string POccStructure_title_tag="HNF(n="; //in case we cannot find POSCAR_START_tag, SHACHAR OLD FORMAT
    string structures_file_str=structures_file_ss.str();
    if(!aurostd::substring2bool(structures_file_str,POSCAR_series_START_tag)){
      if(!aurostd::substring2bool(structures_file_str,POSCAR_START_tag)){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown format for "+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE+" file",_FILE_WRONG_FORMAT_);
      }
      aurostd::string2vectorstring(structures_file_str,vlines); //aurostd::efile2vectorstring(structures_file,vlines);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for l_supercell_sets_size" << endl;}
      for(uint iline=0;iline<vlines.size();iline++){
        if(aurostd::substring2bool(vlines[iline],POSCAR_START_tag)){l_supercell_sets_size++;}
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " l_supercell_sets_size=" << l_supercell_sets_size << endl;}
      if(l_supercell_sets_size==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find l_supercell_sets_size",_RUNTIME_INIT_);}
      i_l_supercell_sets_size=0;
      for(uint iline=0;iline<vlines.size();iline++){
        if(aurostd::substring2bool(vlines[iline],POSCAR_STOP_tag)){
          if(pocc_addendum_ss.str().empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"STOP tag found before START",_FILE_CORRUPT_);}
          aurostd::StringSubst(vlines[iline],POSCAR_STOP_tag,POSCAR_series_STOP_tag+pocc_addendum_ss.str());
          pocc_addendum_ss.str("");
        }
        if(aurostd::substring2bool(vlines[iline],POSCAR_START_tag)){
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " vlines[iline=" << iline << "]=\"" << vlines[iline] << "\"" << endl;}
          if(i_l_supercell_sets_size>=l_supercell_sets_size){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found more unique supercells than expected",_INDEX_MISMATCH_);}
          if(iline+1>vlines.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot grab title at iline="+aurostd::utype2string(iline),_FILE_CORRUPT_);} //see if title is accessible
          parsePOccHashFromXStructureTitle(vlines[iline+1],pocc_hash);
          if(pocc_hash.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not find pocc_hash at iline="+aurostd::utype2string(iline+1),_FILE_CORRUPT_);}
          pocc_addendum_ss << "POCC_";
          pocc_addendum_ss << std::setfill('0') << std::setw(aurostd::getZeroPadding(l_supercell_sets_size)) << (i_l_supercell_sets_size++)+1 << "_";
          pocc_addendum_ss << pocc_hash;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_addendum=" << pocc_addendum_ss.str() << endl;}
          //everything checks out, patch the line
          aurostd::StringSubst(vlines[iline],POSCAR_START_tag,POSCAR_series_START_tag+pocc_addendum_ss.str());
          if(!aurostd::substring2bool(vlines[iline+1],pocc_hash)){aurostd::StringSubst(vlines[iline+1],POccStructure_title_tag,pocc_hash+" "+POccStructure_title_tag);} //SHACHAR OLD FORMAT  //iline+1 for title line
          if(!aurostd::substring2bool(vlines[iline+1],DOI_POCC)){vlines[iline+1]=aurostd::RemoveWhiteSpacesFromTheBack(vlines[iline+1])+DOI_POCC;}  //SHACHAR OLD FORMAT  //iline+1 for title line
        }
        aurostd::StringSubst(vlines[iline],"#",""); //NO continue
      }
      //stringstream structures_file_ss_new;
      structures_file_ss.str("");
      for(uint iline=0;iline<vlines.size();iline++){structures_file_ss << vlines[iline] << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " patched file=" << endl;cerr << structures_file_ss.str() << endl;}
      return true;
    }
    return false;
  }

  void patchStructuresFile(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    //POCC_ALL_SUPERCELLS_FILE
    string structures_file="";
    stringstream structures_file_ss;
    if(patchStructuresAllFile(aflags,structures_file,structures_file_ss,FileMESSAGE,oss)){
      message << "Writing out patched " << POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      string structures_file_new=structures_file;
      aurostd::StringSubst(structures_file_new,POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE,POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE+".OLD");
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " structures_file=" << structures_file << endl;
        cerr << __AFLOW_FUNC__ << " structures_file_new=" << structures_file_new << endl;
      }
      if(!aurostd::file2file(structures_file,structures_file_new)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not move structures_file to structures_file_new="+structures_file_new,_RUNTIME_ERROR_);}
      aurostd::stringstream2file(structures_file_ss,aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE);structures_file_ss.str("");
    }

    //POCC_UNIQUE_SUPERCELLS_FILE
    if(patchStructuresUniqueFile(aflags,structures_file,structures_file_ss,FileMESSAGE,oss)){
      message << "Writing out patched " << POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      string structures_file_new=structures_file;
      aurostd::StringSubst(structures_file_new,POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE+".OLD");
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " structures_file=" << structures_file << endl;
        cerr << __AFLOW_FUNC__ << " structures_file_new=" << structures_file_new << endl;
      }
      if(!aurostd::file2file(structures_file,structures_file_new)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Could not move structures_file to structures_file_new="+structures_file_new,_RUNTIME_ERROR_);}
      aurostd::stringstream2file(structures_file_ss,aflags.Directory+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE);structures_file_ss.str("");
    }

    if(!(aurostd::EFileExist(aflags.Directory+"/PARTCAR") || aurostd::FileExist(aflags.Directory+"/PARTCAR"))){   //patch: write out PARTCAR, should really only happen during generateStructures()
      message << "Writing out PARTCAR";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      string AflowIn_file,AflowIn;
      KBIN::getAflowInFromAFlags(aflags,AflowIn_file,AflowIn,FileMESSAGE,oss);
      xstructure xstr_pocc=pocc::extractPARTCAR(AflowIn); //prefer to pull from AflowIn input vs. aflags
      stringstream xstr_pocc_ss;
      xstr_pocc_ss << xstr_pocc;
      aurostd::stringstream2file(xstr_pocc_ss,aflags.Directory+"/PARTCAR");
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }
} // namespace pocc

namespace pocc {
  vector<string> getElementsList() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    uint skip_every=2;
    vector<string> elements;
    vector<string> VSTD_ELEMENTS_LIST;
    aurostd::string2tokens(STD_ELEMENTS_LIST,VSTD_ELEMENTS_LIST," ");
    for(uint i=0;i<skip_every+1;i++){
      for(uint j=i;j<VSTD_ELEMENTS_LIST.size();j+=skip_every+1){
        elements.push_back(VSTD_ELEMENTS_LIST[j]);
      }
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " elements_list=" << aurostd::joinWDelimiter(elements,",") << endl;}
    return elements;
  }
} // namespace pocc

namespace pocc {
  bool sortPSCSetsUFFEnergy(const POccSuperCellSet& a, const POccSuperCellSet& b){return a.getUFFEnergy()<b.getUFFEnergy();}
} // namespace pocc

namespace pocc {
  double getHmix(const xvector<double>& xv_energies,const xvector<double>& xv_dgs){double dg_total;return getHmix(xv_energies,xv_dgs,dg_total);}
  double getHmix(const xvector<double>& xv_energies,const xvector<double>& xv_dgs,double& dg_total){return aurostd::meanWeighted(xv_energies,xv_dgs,dg_total);}
  double getEFA(const xvector<double>& xv_energies,const xvector<double>& xv_dgs){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    double dg_total=0.0;
    double Hmix=getHmix(xv_energies,xv_dgs,dg_total);
    double sigma=0.0;
    for(int i=xv_dgs.lrows;i<=xv_dgs.urows;i++){sigma+=xv_dgs[i]*pow(xv_energies[i]-Hmix,2.0);}
    sigma/=(dg_total-1);
    sigma=sqrt(sigma);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " sigma=" << sigma << endl;}
    double efa=1.0/sigma;
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " efa=" << efa << endl;}
    return efa;
  }

  void POccCalculator::StructuresAllFile2SupercellSets(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    vector<string> vlines,vtokens,vtokens2;

    string structures_file="";
    stringstream structures_file_ss;
    patchStructuresAllFile(m_aflags,structures_file,structures_file_ss,*p_FileMESSAGE,*p_oss);
    aurostd::string2vectorstring(structures_file_ss.str(),vlines);
    string arun_directory="",pocc_directory_abs="",qmvasp_path="";
    xQMVASP qmvasp;
    m_ARUN_directories.clear();
    for(uint iline=0;iline<vlines.size();iline++){
      if(aurostd::substring2bool(vlines[iline],POCC_AFLOWIN_tag+"STRUCTURES_GROUP")){
        l_supercell_sets.push_back(POccSuperCellSet());
        continue;
      }
      if(aurostd::substring2bool(vlines[iline],POCC_AFLOWIN_tag+"UFF_ENERGY=")){
        aurostd::string2tokens(vlines[iline],vtokens,"=");
        if(vtokens.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown UFF_ENERGY line format",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(vtokens[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"UFF_ENERGY is not a double",_VALUE_ERROR_);}
        l_supercell_sets.back().m_psc_set.push_back(POccSuperCell()); //this comes before structure
        l_supercell_sets.back().m_psc_set.back().m_energy_uff=aurostd::string2utype<double>(vtokens[1]);
        continue;
      }
      if(aurostd::substring2bool(vlines[iline],POSCAR_POCC_series_START_tag)){
        //hnf_matrix
        aurostd::string2tokens(vlines[iline],vtokens,"H");
        if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (H1) line",_FILE_CORRUPT_);}
        aurostd::string2tokens(vtokens.back(),vtokens2,"C");
        if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (H2) line",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(vtokens2[0])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"hnf_index is not an integer",_VALUE_ILLEGAL_);}
        l_supercell_sets.back().m_psc_set.back().m_hnf_index=aurostd::string2utype<unsigned long long int>(vtokens2[0]);
        //site_configuration
        aurostd::string2tokens(vlines[iline],vtokens,"C");
        if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (C1) line",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(aurostd::RemoveWhiteSpaces(vtokens.back()))){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"site_config_index is not an integer",_VALUE_ERROR_);}
        l_supercell_sets.back().m_psc_set.back().m_site_config_index=aurostd::string2utype<unsigned long long int>(aurostd::RemoveWhiteSpaces(vtokens.back()));
        //degeneracy
        l_supercell_sets.back().m_psc_set.back().m_degeneracy=1;
        if(iline+1>vlines.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot grab title at iline="+aurostd::utype2string(iline),_FILE_CORRUPT_);} //see if title is accessible
        if(!aurostd::substring2bool(vlines[iline+1],"DG=")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"DG= is not present in POSCAR title line",_FILE_CORRUPT_);}
        aurostd::string2tokens(vlines[iline+1],vtokens," ");
        for(uint i=0;i<vtokens.size();i++){
          if(aurostd::substring2bool(vtokens[i],"DG=")){
            aurostd::string2tokens(vtokens[i],vtokens2,"=");
            if(vtokens2.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot parse at DG= (vtokens2.size()=="+aurostd::utype2string(vtokens2.size())+")",_FILE_CORRUPT_);}
            if(!aurostd::isfloat(vtokens2[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"degeneracy is not an integer",_VALUE_ILLEGAL_);}
            l_supercell_sets.back().m_psc_set.back().m_degeneracy*=aurostd::string2utype<unsigned long long int>(vtokens2[1]); //backwards compatible with old pocc scheme, but this should be obsolete
          }
        }
        //directory
        arun_directory=vlines[iline];
        aurostd::StringSubst(arun_directory,POSCAR_series_START_tag,"");
        aurostd::string2tokens(arun_directory,vtokens,"_");  //POCC_01_01_H0C0
        if(vtokens.size()!=4){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown POcc hash format (_)",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(vtokens[2])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine POcc structure group index",_VALUE_ILLEGAL_);}
        if(aurostd::string2utype<int>(vtokens[2])==1){
          arun_directory="ARUN."+vtokens[0]+"_"+vtokens[1]+"_"+vtokens[3];
          m_ARUN_directories.push_back(arun_directory);
          pocc_directory_abs=m_aflags.Directory+"/"+arun_directory;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_directory_abs=" << pocc_directory_abs << endl;}
          if(!aurostd::IsDirectory(pocc_directory_abs)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POcc directory [dir="+pocc_directory_abs+"] not found",_FILE_NOT_FOUND_);}
          //m_ARUN_directories.push_back(pocc_directory_abs); //do not save full path
        }
        //debug
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " found hnf_index=" << l_supercell_sets.back().m_psc_set.back().m_hnf_index;
          cerr << " site_config_index=" << l_supercell_sets.back().m_psc_set.back().m_site_config_index;
          cerr << " degeneracy=" << l_supercell_sets.back().m_psc_set.back().m_degeneracy;
          cerr << endl;
        }
      }
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " l_supercell_sets.size()=" << l_supercell_sets.size() << endl;
      unsigned long long int isupercell=0;
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        isupercell=std::distance(l_supercell_sets.begin(),it);
        cerr << __AFLOW_FUNC__ << " l_supercell_sets[i=" << isupercell << "].m_psc_set.size()=" << (*it).m_psc_set.size() << endl;
      }
    }
  }

  void POccCalculator::StructuresUniqueFile2SupercellSets(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    vector<string> vlines,vtokens,vtokens2;

    string structures_file="";
    stringstream structures_file_ss;
    patchStructuresUniqueFile(m_aflags,structures_file,structures_file_ss,*p_FileMESSAGE,*p_oss);
    aurostd::string2vectorstring(structures_file_ss.str(),vlines);
    string arun_directory="",pocc_directory_abs="",qmvasp_path="";
    xQMVASP qmvasp;
    m_ARUN_directories.clear();
    for(uint iline=0;iline<vlines.size();iline++){
      if(aurostd::substring2bool(vlines[iline],POSCAR_POCC_series_START_tag)){
        l_supercell_sets.push_back(POccSuperCellSet());
        l_supercell_sets.back().m_psc_set.push_back(POccSuperCell());
        //hnf_matrix
        aurostd::string2tokens(vlines[iline],vtokens,"H");
        if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (H1) line",_FILE_CORRUPT_);}
        aurostd::string2tokens(vtokens.back(),vtokens2,"C");
        if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (H2) line",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(vtokens2[0])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"hnf_index is not an integer",_VALUE_ILLEGAL_);}
        l_supercell_sets.back().m_psc_set.back().m_hnf_index=aurostd::string2utype<unsigned long long int>(vtokens2[0]);
        //site_configuration
        aurostd::string2tokens(vlines[iline],vtokens,"C");
        if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (C1) line",_FILE_CORRUPT_);}
        if(!aurostd::isfloat(aurostd::RemoveWhiteSpaces(vtokens.back()))){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"site_config_index is not an integer",_VALUE_ILLEGAL_);}
        l_supercell_sets.back().m_psc_set.back().m_site_config_index=aurostd::string2utype<unsigned long long int>(aurostd::RemoveWhiteSpaces(vtokens.back()));
        //degeneracy
        l_supercell_sets.back().m_psc_set.back().m_degeneracy=1;
        if(iline+1>vlines.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot grab title at iline="+aurostd::utype2string(iline),_FILE_CORRUPT_);} //see if title is accessible
        if(!aurostd::substring2bool(vlines[iline+1],"DG=")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"DG= is not present in POSCAR title line",_FILE_CORRUPT_);}
        aurostd::string2tokens(vlines[iline+1],vtokens," ");
        for(uint i=0;i<vtokens.size();i++){
          if(aurostd::substring2bool(vtokens[i],"DG=")){
            aurostd::string2tokens(vtokens[i],vtokens2,"=");
            if(vtokens2.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot parse at DG= (vtokens2.size()=="+aurostd::utype2string(vtokens2.size())+")",_FILE_CORRUPT_);}
            if(!aurostd::isfloat(vtokens2[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"degeneracy is not an integer",_VALUE_ILLEGAL_);}
            l_supercell_sets.back().m_psc_set.back().m_degeneracy*=aurostd::string2utype<unsigned long long int>(vtokens2[1]); //backwards compatible with old pocc scheme, but this should be obsolete
          }
        }
        //directory
        arun_directory=vlines[iline];
        aurostd::StringSubst(arun_directory,POSCAR_series_START_tag,"");  //POCC_01_H0C0
        arun_directory="ARUN."+arun_directory;
        m_ARUN_directories.push_back(arun_directory);
        pocc_directory_abs=m_aflags.Directory+"/"+arun_directory;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_directory_abs=" << pocc_directory_abs << endl;}
        if(!aurostd::IsDirectory(pocc_directory_abs)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POcc directory [dir="+pocc_directory_abs+"] not found",_FILE_CORRUPT_);}
        //m_ARUN_directories.push_back(pocc_directory_abs); //do not save full path
        //debug
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " found hnf_index=" << l_supercell_sets.back().m_psc_set.back().m_hnf_index;
          cerr << " site_config_index=" << l_supercell_sets.back().m_psc_set.back().m_site_config_index;
          cerr << " degeneracy=" << l_supercell_sets.back().m_psc_set.back().m_degeneracy;
          cerr << endl;
        }
      }
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " l_supercell_sets.size()=" << l_supercell_sets.size() << endl;
      unsigned long long int isupercell=0;
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        isupercell=std::distance(l_supercell_sets.begin(),it);
        cerr << __AFLOW_FUNC__ << " l_supercell_sets[i=" << isupercell << "].m_psc_set.size()=" << (*it).m_psc_set.size() << endl;
      }
    }
  }

  void POccCalculator::setDFTEnergies() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}
    if(m_ARUN_directories.size()!=l_supercell_sets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()!=l_supercell_sets.size()",_RUNTIME_ERROR_);}

    xQMVASP qmvasp;
    unsigned long long int isupercell=0;
    string qmvasp_filename="";
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      const string& pocc_directory_abs=getARUNDirectoryPath(isupercell);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " look at pocc_directory_abs=" << pocc_directory_abs << endl;}
      //energy_dft
      qmvasp_filename=pocc_directory_abs+"/"+DEFAULT_AFLOW_QMVASP_OUT;
      if(!aurostd::EFileExist(qmvasp_filename,qmvasp_filename)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No qmvasp file found [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);}
      qmvasp.GetPropertiesFile(qmvasp_filename);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " qmvasp.H_atom_relax=" << qmvasp.H_atom_relax << endl;
        cerr << __AFLOW_FUNC__ << " qmvasp.H_atom_static=" << qmvasp.H_atom_static << endl;
      }
      //CO20200921 - we pulled H_atom_static previously, but this does not match what is being done in aflowlib_libraries
      //we pull relax and not static because relaxation runs come with stresses, pressures, etc. from which we can derive PV (for enthalpies)
      //use H_atom_relax instead
      (*it).m_energy_dft=qmvasp.H_atom_relax; //this will be the LAST relax, relax matches with aflowlib_libraries (NOT static)
      if((*it).m_energy_dft==AUROSTD_NAN){(*it).m_energy_dft=qmvasp.H_atom_relax;}
      if((*it).m_energy_dft==AUROSTD_NAN){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No H_atom found in qmvasp [dir="+pocc_directory_abs+"]",_FILE_CORRUPT_);}
      //energy_dft_ground
      if((*it).m_energy_dft<m_energy_dft_ground){
        m_energy_dft_ground=(*it).m_energy_dft;
        m_ARUN_directory_ground=isupercell;
      }
      //m_energy_dft_ground=std::min(m_energy_dft_ground,(*it).m_energy_dft);
      if(LDEBUG){cerr << " H_atom[isupercell=" << isupercell << "]=" << (*it).m_energy_dft << endl;}
    }
    if(LDEBUG){cerr << " H_atom_ground=" << m_energy_dft_ground << endl;}

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }

  void POccCalculator::setPOccStructureProbabilities(double temperature) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " temperature=" << temperature << endl;}
    if(l_supercell_sets.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"l_supercell_sets.size()==0",_RUNTIME_ERROR_);}
    if(m_energy_dft_ground==AUROSTD_MAX_DOUBLE){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_energy_dft_ground not set",_RUNTIME_ERROR_);}

    if(std::signbit(temperature)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Negative temperature found",_VALUE_ILLEGAL_);}
    if(aurostd::isequal(temperature,0.0)){
      temperature=_ZERO_TOL_; //1e-6;
      message << "Converting T=0 -> T=" << temperature << " for Boltzmann probabilities";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); //WARNING
    }

    double denom=0.0;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      denom+=(*it).getDegeneracy()*exp( -( (*it).m_energy_dft-m_energy_dft_ground ) / (KBOLTZEV*temperature) );
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " denom=" << denom << endl;}

    double prob_total=0.0;
    unsigned long long int isupercell=0;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " degeneracy[isupercell=" << isupercell << "]=" << (*it).getDegeneracy() << endl;}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_energy_dft[isupercell=" << isupercell << "]=" << (*it).m_energy_dft << endl;}
      (*it).m_probability=(*it).getDegeneracy()*exp( -( (*it).m_energy_dft-m_energy_dft_ground ) / (KBOLTZEV*temperature) ) / denom;
      prob_total+=(*it).m_probability;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " prob[isupercell=" << isupercell << "]=" << (*it).m_probability << endl;}
    }
    if(!aurostd::isequal(prob_total,1.0)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"prob_total != 1.0",_RUNTIME_ERROR_);}
  }

  void POccCalculator::setEFA(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(m_convolution){
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"skipping cEFA calculation, not coded yet",m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      return;
    }
    if(l_supercell_sets.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"l_supercell_sets.size()==0",_RUNTIME_ERROR_);}
    xvector<double> xv_dgs(l_supercell_sets.size()),xv_energies(l_supercell_sets.size());
    unsigned long long int isupercell=0;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      xv_dgs[isupercell+xv_dgs.lrows]=(*it).getDegeneracy();
      xv_energies[isupercell+xv_energies.lrows]=(*it).m_energy_dft;
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " xv_dgs[i=" << isupercell+xv_dgs.lrows << "]=" << xv_dgs[isupercell+xv_dgs.lrows];
        cerr << " xv_energies[" << isupercell+xv_energies.lrows << "]=" << xv_energies[isupercell+xv_energies.lrows] << endl;
      }
    }
    m_Hmix=pocc::getHmix(xv_energies,xv_dgs);
    m_efa=pocc::getEFA(xv_energies,xv_dgs);
  }

  void getTemperatureStringParameters(int& temperature_precision,bool& temperatures_int,int& zero_padding_temperature){
    vector<double> v_temperatures=getVTemperatures(DEFAULT_POCC_TEMPERATURE_STRING);
    return getTemperatureStringParameters(v_temperatures,temperature_precision,temperatures_int,zero_padding_temperature);
  }
  void getTemperatureStringParameters(const vector<double>& v_temperatures,int& temperature_precision,bool& temperatures_int,int& zero_padding_temperature){
    temperature_precision=TEMPERATURE_PRECISION;
    temperatures_int=true;
    zero_padding_temperature=4;
    if(!v_temperatures.empty()){
      temperatures_int=true;
      for(uint itemp=0;itemp<v_temperatures.size();itemp++){if(!aurostd::isinteger(v_temperatures[itemp])){temperatures_int=false;break;}}  //found a non-int temperature
      zero_padding_temperature=aurostd::getZeroPadding(max(v_temperatures))+(temperatures_int ? 0 : TEMPERATURE_PRECISION+1); //+1 for decimal place
    }
  }

  string POccCalculator::getTemperatureString(double temperature) const {
    return pocc::getTemperatureString(temperature,m_temperature_precision,m_temperatures_int,m_zero_padding_temperature);
  }
  string getTemperatureString(double temperature,int precision,bool temperatures_int,int zero_padding_temperature) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream t_ss;
    if(temperatures_int==false){t_ss.setf(std::ios::fixed,std::ios::floatfield);t_ss.precision(precision);}
    t_ss.width(zero_padding_temperature);t_ss.fill('0');
    t_ss << temperature;

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " temperature=" << temperature << endl;
      cerr << __AFLOW_FUNC__ << " temperatures_int=" << temperatures_int << endl;
      cerr << __AFLOW_FUNC__ << " zero_padding=" << zero_padding_temperature << endl;
      cerr << __AFLOW_FUNC__ << " temperature_precision=" << precision << endl;
      cerr << __AFLOW_FUNC__ << " temperature_string=" << t_ss.str() << endl;
    }

    return t_ss.str();
  }

  void POccCalculator::setAvgRDF(double temperature){
    bool LDEBUG=(true || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}
    if(m_ARUN_directories.size()!=l_supercell_sets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()!=l_supercell_sets.size()",_RUNTIME_ERROR_);}

    setPOccStructureProbabilities(temperature); //done in calculateRELAXProperties() - repetita iuvant

    //check if pdos can be averaged
    //depends not only the number of atoms (after relaxation), but also whether species/num_each_type are identical
    bool perform_rdf=true;
    vector<string> vspecies,_vspecies;
    deque<int> num_each_type;
    xstructure xstr;
    unsigned long long int isupercell=0;
    string POSCAR_file="";
    bool found_POSCAR_file=false;
    uint i=0;
    bool raw_counts=false;  //make options for these later
    double sigma=0; //make options for these later
    int window_gaussian=0;  //make options for these later
    aurostd::xmatrix<double> rdf_all;
    m_rdf_all.clear();
    m_rdf_rmax=DEFAULT_POCC_RDF_RMAX;
    m_rdf_nbins=DEFAULT_POCC_RDF_NBINS;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end()&&perform_rdf==true;++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      if(aflowlib::GetSpeciesDirectory(getARUNDirectoryPath(isupercell),_vspecies)==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot extract vspecies from "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
      found_POSCAR_file=false;
      if(found_POSCAR_file==false){found_POSCAR_file=aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/CONTCAR.relax",POSCAR_file);}
      if(found_POSCAR_file==false){found_POSCAR_file=aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/CONTCAR.relax"+aurostd::utype2string(m_relaxation_max),POSCAR_file);}
      if(found_POSCAR_file==false){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"CONTCAR.relax"+aurostd::utype2string(m_relaxation_max)+" not found in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);
      }
      message << "Processing CONTCAR.relax of " << m_ARUN_directories[isupercell] << " for the RDF (takes some time)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      xstr.clear();
      xstr.initialize(POSCAR_file,IOVASP_POSCAR);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " xstr[isupercell=" << isupercell << "]=" << endl << xstr;
        cerr << __AFLOW_FUNC__ << " xstr.atoms.size()=" << xstr.atoms.size() << endl;
        cerr << __AFLOW_FUNC__ << " xstr.num_each_type=" << aurostd::joinWDelimiter(xstr.num_each_type,",") << endl;
        cerr << __AFLOW_FUNC__ << " vspecies=" << aurostd::joinWDelimiter(_vspecies,",") << endl;
      }
      if(xstr.species.empty()||(xstr.species.size()>0 && xstr.species[0].empty())){xstr.SetSpecies(aurostd::vector2deque(_vspecies));}
      if(isupercell==0){
        vspecies.clear();for(i=0;i<_vspecies.size();i++){vspecies.push_back(_vspecies[i]);}
        num_each_type.clear();for(i=0;i<xstr.num_each_type.size();i++){num_each_type.push_back(xstr.num_each_type[i]);}
        //
        pflow::GetRDF(xstr,rdf_all,m_rdf_rmax,m_rdf_nbins,raw_counts,sigma,window_gaussian);
        m_rdf_all=(*it).m_probability * rdf_all;
      }else{
        if(perform_rdf && (vspecies!=_vspecies)){
          message << "Found a mismatch in vspecies in " << m_ARUN_directories[isupercell] << ", not performing pdos averaging";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
          perform_rdf=false;
        }
        if(perform_rdf && (num_each_type!=xstr.num_each_type)){
          message << "Found a mismatch in num_each_type in " << m_ARUN_directories[isupercell] << ", not performing pdos averaging";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
          perform_rdf=false;
        }
        //
        pflow::GetRDF(xstr,rdf_all,m_rdf_rmax,m_rdf_nbins,raw_counts,sigma,window_gaussian);
        m_rdf_all+=(*it).m_probability * rdf_all;
      }
      PrintRDF(xstr,m_rdf_rmax,m_rdf_nbins,rdf_all,*p_oss);
    }
    if(perform_rdf==false){m_rdf_all.clear();return;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_rdf_all=" << endl << m_rdf_all << endl;}
  }

  void POccCalculator::setAvgDOSCAR(double temperature){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}
    if(m_ARUN_directories.size()!=l_supercell_sets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()!=l_supercell_sets.size()",_RUNTIME_ERROR_);}

    setPOccStructureProbabilities(temperature); //done in calculateRELAXProperties() - repetita iuvant

    //check if pdos can be averaged
    //depends not only the number of atoms (after relaxation), but also whether species/num_each_type are identical
    bool perform_pdos_avg=true;
    vector<string> vspecies,_vspecies;
    deque<int> num_each_type;
    vector<uint> vnatoms;
    xstructure xstr;
    unsigned long long int isupercell=0;
    string POSCAR_file="";
    uint i=0;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      if(aflowlib::GetSpeciesDirectory(getARUNDirectoryPath(isupercell),_vspecies)==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot extract vspecies from "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
      if(!aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/POSCAR.static",POSCAR_file)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POSCAR.static not found in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
      xstr.clear();
      xstr.initialize(POSCAR_file,IOVASP_POSCAR);
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " xstr[isupercell=" << isupercell << "]=" << endl << xstr;
        cerr << __AFLOW_FUNC__ << " xstr.atoms.size()=" << xstr.atoms.size() << endl;
        cerr << __AFLOW_FUNC__ << " xstr.num_each_type=" << aurostd::joinWDelimiter(xstr.num_each_type,",") << endl;
        cerr << __AFLOW_FUNC__ << " vspecies=" << aurostd::joinWDelimiter(_vspecies,",") << endl;
      }
      vnatoms.push_back(xstr.atoms.size());
      if(isupercell==0){
        vspecies.clear();for(i=0;i<_vspecies.size();i++){vspecies.push_back(_vspecies[i]);}
        num_each_type.clear();for(i=0;i<xstr.num_each_type.size();i++){num_each_type.push_back(xstr.num_each_type[i]);}
      }else{
        if(perform_pdos_avg && (vspecies!=_vspecies)){
          message << "Found a mismatch in vspecies in " << m_ARUN_directories[isupercell] << ", not performing pdos averaging";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
          perform_pdos_avg=false;
        }
        if(perform_pdos_avg && (num_each_type!=xstr.num_each_type)){
          message << "Found a mismatch in num_each_type in " << m_ARUN_directories[isupercell] << ", not performing pdos averaging";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
          perform_pdos_avg=false;
        }
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vnatoms=" << aurostd::joinWDelimiter(vnatoms,",") << endl;}

    xDOSCAR xdoscar;
    string DOSCAR_file="";
    m_Egap_DOS_net=0.0;
    bool metal_found=false,insulator_found=false;
    uint ispin=0,iatom=0,iorbital=0,vDOS_avg_size=0;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      if(!aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/DOSCAR.static",DOSCAR_file)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"DOSCAR.static not found in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
      message << "Processing xDOSCAR of " << m_ARUN_directories[isupercell];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      xdoscar.GetPropertiesFile(DOSCAR_file);
      message << "xDOSCAR read";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      xdoscar.GetBandGap();
      for(ispin=0;ispin<xdoscar.Egap.size();ispin++){
        message << "ISPIN=" << ispin << ", Egap_DOS=" << xdoscar.Egap[ispin];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      }
      if(isupercell==0){ //set all to 0
        m_xdoscar=xdoscar;  //import all properties
        //set all values to 0 for ensemble averages
        //[NOT necessary, should be the same across all ARUNS]std::fill(m_xdoscar.venergy.begin(),m_xdoscar.venergy.end(),0.0);
        std::fill(m_xdoscar.venergyEf.begin(),m_xdoscar.venergyEf.end(),0.0); //NOT the same across all ARUNS
        for(ispin=0;ispin<m_xdoscar.viDOS.size();ispin++){std::fill(m_xdoscar.viDOS[ispin].begin(),m_xdoscar.viDOS[ispin].end(),0.0);}

        if(m_xdoscar.vDOS.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.vDOS.size()==0",_INDEX_MISMATCH_);}
        if(perform_pdos_avg){vDOS_avg_size=m_xdoscar.vDOS.size();}
        else{
          m_xdoscar.vDOS.resize(1); //resize to length 1 (total index)
          m_xdoscar.number_atoms=0;
          vDOS_avg_size=m_xdoscar.vDOS.size();  //only avg the totals, the extra pdos channel is a hack to print orbital sums
          m_xdoscar.addAtomChannel(); //add single pdos for orbital sums
        }
        m_xdoscar.resetVDOS();
        m_xdoscar.Efermi=0.0;
        m_Egap_DOS.assign(xdoscar.Egap.size(),0.0);
      }else{ //check that all dimensions match
        //these are VERY important checks, make sure venergyEf can be ensemble averaged
        //otherwise none of the properties can be averaged (comparing values at two different energies)
        if(m_xdoscar.energy_max!=xdoscar.energy_max){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.energy_max!=xdoscar.energy_max",_INDEX_MISMATCH_);}
        if(m_xdoscar.energy_min!=xdoscar.energy_min){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.energy_min!=xdoscar.energy_min",_INDEX_MISMATCH_);}
        if(m_xdoscar.number_energies!=xdoscar.number_energies){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.number_energies!=xdoscar.number_energies",_INDEX_MISMATCH_);}

        if(m_xdoscar.venergy.size()!=xdoscar.venergy.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.venergy.size()!=xdoscar.venergy.size()",_INDEX_MISMATCH_);}
        if(m_xdoscar.venergyEf.size()!=xdoscar.venergyEf.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.venergyEf.size()!=xdoscar.venergyEf.size()",_INDEX_MISMATCH_);}
        if(m_xdoscar.viDOS.size()!=xdoscar.viDOS.size()){
          message << "Mismatch SPIN-ON/SPIN-OFF settings, attempting to rectify";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
          if(m_xdoscar.viDOS.size()==1){ //make duplicate for spin-off
            m_xdoscar.convertSpinOFF2ON();
            m_Egap_DOS.push_back(m_Egap_DOS.back());  //need to extend m_Egap_DOS too
          }
          else if(xdoscar.viDOS.size()==1){xdoscar.convertSpinOFF2ON();}
          if(m_xdoscar.viDOS.size()!=xdoscar.viDOS.size()){
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.viDOS.size()!=xdoscar.viDOS.size() ["+aurostd::utype2string(m_xdoscar.viDOS.size())+"!="+aurostd::utype2string(xdoscar.viDOS.size())+"]",_INDEX_MISMATCH_);
          }
        }
        for(ispin=0;ispin<m_xdoscar.viDOS.size();ispin++){
          if(m_xdoscar.viDOS[ispin].size()!=xdoscar.viDOS[ispin].size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.viDOS[ispin].size()!=xdoscar.viDOS[ispin].size()",_INDEX_MISMATCH_);}
        }
        if(perform_pdos_avg){
          if(m_xdoscar.vDOS.size()!=xdoscar.vDOS.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.vDOS.size()!=xdoscar.vDOS.size()",_INDEX_MISMATCH_);}
        }else{
          if(xdoscar.vDOS.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"xdoscar.vDOS.size()==0",_INDEX_MISMATCH_);}
        }
        for(iatom=0;iatom<vDOS_avg_size;iatom++){
          if(m_xdoscar.vDOS[iatom].size()!=xdoscar.vDOS[iatom].size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.vDOS[iatom].size()!=xdoscar.vDOS[iatom].size()",_INDEX_MISMATCH_);}
          for(iorbital=0;iorbital<m_xdoscar.vDOS[iatom].size();iorbital++){
            if(m_xdoscar.vDOS[iatom][iorbital].size()!=xdoscar.vDOS[iatom][iorbital].size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.vDOS[iatom][iorbital].size()!=xdoscar.vDOS[iatom][iorbital].size()",_INDEX_MISMATCH_);}
            for(ispin=0;ispin<m_xdoscar.vDOS[iatom][iorbital].size();ispin++){
              if(m_xdoscar.vDOS[iatom][iorbital][ispin].size()!=xdoscar.vDOS[iatom][iorbital][ispin].size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_xdoscar.vDOS[iatom][iorbital][ispin].size()!=xdoscar.vDOS[iatom][iorbital][ispin].size()",_INDEX_MISMATCH_);}
            }
          }
        }
        if(m_Egap_DOS.size()!=xdoscar.Egap.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_Egap_DOS.size()!=xdoscar.Egap.size()",_INDEX_MISMATCH_);}
      }

      //these are all VECTORS and not xvectors, if we change in the future, this can be reduced by ienergy for loop
      //[NOT necessary, should be the same across all ARUNS]for(ienergy=0;ienergy<m_xdoscar.venergy.size();ienergy++){m_xdoscar.venergy[ienergy]+=( (*it).m_probability*xdoscar.venergy[ienergy] );}
      uint ienergy=0;
      for(ienergy=0;ienergy<m_xdoscar.venergyEf.size();ienergy++){m_xdoscar.venergyEf[ienergy]+=( (*it).m_probability*xdoscar.venergyEf[ienergy] );} //NOT the same across all ARUNS
      for(ispin=0;ispin<m_xdoscar.viDOS.size();ispin++){
        for(ienergy=0;ienergy<m_xdoscar.viDOS[ispin].size();ienergy++){
          m_xdoscar.viDOS[ispin][ienergy]+=( (*it).m_probability * xdoscar.viDOS[ispin][ienergy] / vnatoms[isupercell] );
        }
      }
      //[sum of pdos != total: need to consider interstitials]if(0){
      for(iatom=0;iatom<vDOS_avg_size;iatom++){  //iatom==0 is total, we can average this too
        for(iorbital=0;iorbital<m_xdoscar.vDOS[iatom].size();iorbital++){  //iorbital==0 is total, we can average this too
          for(ispin=0;ispin<m_xdoscar.vDOS[iatom][iorbital].size();ispin++){
            for(ienergy=0;ienergy<m_xdoscar.vDOS[iatom][iorbital][ispin].size();ienergy++){
              m_xdoscar.vDOS[iatom][iorbital][ispin][ienergy]+=( (*it).m_probability * xdoscar.vDOS[iatom][iorbital][ispin][ienergy] / vnatoms[isupercell] );
            }
          }
        }
      }
      if(!perform_pdos_avg){
        for(iatom=1;iatom<xdoscar.vDOS.size();iatom++){  //iatom==1 starts first atom, sum from here onward
          while(xdoscar.vDOS[iatom].size()>m_xdoscar.vDOS[iatom].size()){m_xdoscar.addOrbitalChannel();}
          for(iorbital=0;iorbital<xdoscar.vDOS[iatom].size();iorbital++){  //iorbital==0 is total, we can average this too
            for(ispin=0;ispin<xdoscar.vDOS[iatom][iorbital].size();ispin++){
              for(ienergy=0;ienergy<xdoscar.vDOS[iatom][iorbital][ispin].size();ienergy++){
                m_xdoscar.vDOS[1][iorbital][ispin][ienergy]+=( (*it).m_probability * xdoscar.vDOS[iatom][iorbital][ispin][ienergy] / vnatoms[isupercell] );
              }
            }
          }
        }
      }
      //[sum of pdos != total: need to consider interstitials}else{
      //[sum of pdos != total: need to consider interstitials  double dos;
      //[sum of pdos != total: need to consider interstitials  for(iatom=1;iatom<vDOS_avg_size;iatom++){  //iatom==0 is total, sum at the end
      //[sum of pdos != total: need to consider interstitials    for(iorbital=1;iorbital<m_xdoscar.vDOS[iatom].size();iorbital++){  //iorbital==0 is total, sum at the end
      //[sum of pdos != total: need to consider interstitials      for(ispin=0;ispin<m_xdoscar.vDOS[iatom][iorbital].size();ispin++){
      //[sum of pdos != total: need to consider interstitials        for(ienergy=0;ienergy<m_xdoscar.vDOS[iatom][iorbital][ispin].size();ienergy++){
      //[sum of pdos != total: need to consider interstitials          dos=( (*it).m_probability*xdoscar.vDOS[iatom][iorbital][ispin][ienergy] );
      //[sum of pdos != total: need to consider interstitials          m_xdoscar.vDOS[iatom][iorbital][ispin][ienergy]+=dos;
      //[sum of pdos != total: need to consider interstitials          m_xdoscar.vDOS[0][iorbital][ispin][ienergy]+=dos;
      //[sum of pdos != total: need to consider interstitials          m_xdoscar.vDOS[0][0][ispin][ienergy]+=dos;
      //[sum of pdos != total: need to consider interstitials        }
      //[sum of pdos != total: need to consider interstitials      }
      //[sum of pdos != total: need to consider interstitials    }
      //[sum of pdos != total: need to consider interstitials  }
      //[sum of pdos != total: need to consider interstitials}

      m_xdoscar.Efermi+=( (*it).m_probability*xdoscar.Efermi ); //venergyEf ensemble average should give the SAME energies as subtracting averaged Ef from venergy

      for(ispin=0;ispin<m_Egap_DOS.size();ispin++){
        //might consider adding vbm/cbm later
        if(xdoscar.Egap[ispin]==_METALGAP_){metal_found=true;} //just be safe, check spin-polarized
        else{insulator_found=true;}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " xdoscar.Egap[ispin=" << ispin << "]=" << xdoscar.Egap[ispin] << endl;}
        m_Egap_DOS[ispin]+=( (*it).m_probability*xdoscar.Egap[ispin]);
      }
      if(xdoscar.Egap_net==_METALGAP_){metal_found=true;}  //just be safe, check _net
      else{insulator_found=true;}
      m_Egap_DOS_net+=( (*it).m_probability*xdoscar.Egap_net);
    }
    if(metal_found && insulator_found){
      message << "Mixed metal and insulator states found, averaging to a metallic gap";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
    }
    if(metal_found){
      for(ispin=0;ispin<m_Egap_DOS.size();ispin++){m_Egap_DOS[ispin]=_METALGAP_;}
      m_Egap_DOS_net=_METALGAP_;
    }

    //set some attributes
    if(m_vflags.AFLOW_SYSTEM.isentry && !m_vflags.AFLOW_SYSTEM.content_string.empty()){m_xdoscar.title=m_vflags.AFLOW_SYSTEM.content_string;}
    m_xdoscar.carstring="POCC"; //important for plotting

    m_xdoscar.GetBandGap(); //gap of averaged vDOS, we don't use this definition anymore

    //from m_xdoscar
    message << "Egap of average DOS (OBSOLETE)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    for(ispin=0;ispin<m_xdoscar.Egap.size();ispin++){
      message << "ISPIN=" << ispin << ", Egap_DOS=" << m_xdoscar.Egap[ispin];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    }
    message << "Egap_DOS_net=" << xdoscar.Egap_net;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    //from POccCalculator
    message << "Average of Egap_DOSs";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    for(ispin=0;ispin<xdoscar.Egap.size();ispin++){
      message << "ISPIN=" << ispin << ", Egap_DOS=" << m_Egap_DOS[ispin];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    }
    message << "Egap_DOS_net=" << m_Egap_DOS_net;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    //remaining properties
    //should really print RELAXED values, but these are still being worked out for pocc
    //default to initial values
    m_xdoscar.Vol=xstr_pocc.GetVolume();
    m_xdoscar.lattice=Getabc_angles(xstr_pocc.scale*xstr_pocc.lattice,DEGREES); //will be size 6, but only the first 3 are printed/read
    for(int i=1;i<=3;i++){m_xdoscar.lattice[i]*=1e-10;}  //convert Angstroms to meters

    //XDOSCAR.OUT
    //string xdoscar_filename=POCC_DOSCAR_FILE+"_T"+aurostd::utype2string(temperature,TEMPERATURE_PRECISION)+"K";
    string xdoscar_filename=POCC_DOSCAR_FILE+"_T"+(*this).getTemperatureString(temperature)+"K";
    message << "Writing out " << xdoscar_filename;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    stringstream pocc_xdoscar_ss;
    pocc_xdoscar_ss << m_xdoscar;
    aurostd::stringstream2file(pocc_xdoscar_ss,getOutputPath()+"/"+xdoscar_filename);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }

  void POccCalculator::setAvgPlasmonicData(double temperature){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}
    if(m_ARUN_directories.size()!=l_supercell_sets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()!=l_supercell_sets.size()",_RUNTIME_ERROR_);}

    setPOccStructureProbabilities(temperature); //done in calculateRELAXProperties() - repetita iuvant

    m_vxplasm.clear();
    xPLASMONICS xplasm;
    string PLASM_file="";
    xcomplex<double> xcomp_tmp(0.0,0.0);
    vector<string> vlines,vtokens;
    uint ieps=0,i=0;
    unsigned long long int isupercell=0;
    for(ieps=0;ieps<m_veps_plasm.size();ieps++){
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        isupercell=std::distance(l_supercell_sets.begin(),it);
        if(!aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/"+m_vfilenames_plasm[ieps],PLASM_file)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,m_vfilenames_plasm[ieps]+" not found in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
        message << "Processing "+m_vfilenames_plasm[ieps]+" (eps=" << m_veps_plasm[ieps] << ") of " << m_ARUN_directories[isupercell];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
        xplasm.GetPropertiesFile(PLASM_file);
        message << "xPLASM read (eps=" << xplasm.eps << ")";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
        if(isupercell==0){  //set all to 0
          m_vxplasm[ieps].venergy.clear();
          m_vxplasm[ieps].veels.clear();
          m_vxplasm[ieps].vdielectric.clear();
          for(i=0;i<xplasm.venergy.size();i++){
            m_vxplasm[ieps].venergy.push_back(xplasm.venergy[i]);
            m_vxplasm[ieps].veels.push_back(0.0);
            xcomp_tmp.re=0.0;xcomp_tmp.im=0.0;
            m_vxplasm[ieps].vdielectric.push_back(xcomp_tmp);
          }
        }else{ //check that all dimensions match
          if(m_vxplasm[ieps].venergy.size()!=xplasm.venergy.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_vxplasm[ieps].venergy.size()!=xplasm.venergy.size() for isupercell="+aurostd::utype2string(isupercell)+" in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
          if(m_vxplasm[ieps].veels.size()!=xplasm.veels.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_vxplasm[ieps].veels.size()!=xplasm.veels.size() for isupercell="+aurostd::utype2string(isupercell)+" in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
          if(m_vxplasm[ieps].vdielectric.size()!=xplasm.vdielectric.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_vxplasm[ieps].vdielectric.size()!=xplasm.vdielectric.size() for isupercell="+aurostd::utype2string(isupercell)+" in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);}
          for(i=0;i<m_vxplasm[ieps].venergy.size();i++){
            if(!aurostd::isequal(m_vxplasm[ieps].venergy[i],xplasm.venergy[i])){
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " venergy: " << m_vxplasm[ieps].venergy[i] << " != " << xplasm.venergy[i] << endl;}
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_vxplasm[ieps].venergy[i="+aurostd::utype2string(i)+"]!=xplasm.venergy[i] for isupercell="+aurostd::utype2string(isupercell)+" in "+m_ARUN_directories[isupercell],_FILE_NOT_FOUND_);
            }
          }
        }
        //do averaging
        for(i=0;i<m_vxplasm[ieps].venergy.size();i++){
          m_vxplasm[ieps].veels[i]+=( (*it).m_probability*xplasm.veels[i] );
          m_vxplasm[ieps].vdielectric[i].re+=( (*it).m_probability*xplasm.vdielectric[i].re );
          m_vxplasm[ieps].vdielectric[i].im+=( (*it).m_probability*xplasm.vdielectric[i].im );
        }
      }
    }
  }

  void POccCalculator::calculateRELAXProperties(double temperature){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}

    //get most relaxed outcar
    uint i=0,max=10;
    m_relaxation_max=AUROSTD_MAX_UINT;
    for(i=1;i<=max;i++){  //i=max;i<=max;i--
      if(aurostd::EFileExist(getARUNDirectoryPath(0)+"/OUTCAR.relax"+aurostd::utype2string(i))){m_relaxation_max=i;}
      else{break;}
    }
    if(m_relaxation_max==AUROSTD_MAX_UINT){return;} //no runs completed

    string OUTCAR_relax="OUTCAR.relax"+aurostd::utype2string(m_relaxation_max);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " OUTCAR_relax=" << OUTCAR_relax << endl;}

    bool found_all_OUTCARs=true;
    for(unsigned long long int isupercell=0;isupercell<m_ARUN_directories.size()&&found_all_OUTCARs;isupercell++){
      if(!aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/"+OUTCAR_relax)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << OUTCAR_relax << " not found in "+m_ARUN_directories[isupercell] << endl;}
        found_all_OUTCARs=false;
      }
    }
    if(!found_all_OUTCARs){return;}

    setPOccStructureProbabilities(temperature);
    setAvgRDF(temperature);
  }

  void POccCalculator::calculateSTATICProperties(double temperature){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}

    bool found_all_OUTCARs=true;
    for(unsigned long long int isupercell=0;isupercell<m_ARUN_directories.size()&&found_all_OUTCARs==true;isupercell++){
      if(!aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/OUTCAR.static")){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " OUTCAR.static not found in "+m_ARUN_directories[isupercell] << endl;}
        found_all_OUTCARs=false;
      }
    }
    if(!found_all_OUTCARs){return;}

    setAvgDOSCAR(temperature);
    if(0){plotAvgDOSCAR(temperature);}  //do not plot as part of LIB2LIB, leave for LIB2RAW
  }

  void POccCalculator::calculatePlasmonicProperties(double temperature){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}

    vector<string> vfiles;

    bool found_all_eps_dat_files=true;
    uint i=0;
    string tmp_str="";
    m_vfilenames_plasm.clear();
    for(unsigned long long int isupercell=0;isupercell<m_ARUN_directories.size()&&found_all_eps_dat_files==true;isupercell++){
      aurostd::DirectoryLS(getARUNDirectoryPath(isupercell),vfiles);
      if(isupercell==0){
        for(i=0;i<vfiles.size();i++){
          if(vfiles[i].find(DEFAULT_AFLOW_PLASMONICS_FILE)!=string::npos){m_vfilenames_plasm.push_back(vfiles[i]);}
        }
      }else{
        for(i=0;i<m_vfilenames_plasm.size();i++){
          if(!aurostd::EFileExist(getARUNDirectoryPath(isupercell)+"/"+m_vfilenames_plasm[i])){
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " "+m_vfilenames_plasm[i]+" not found in "+m_ARUN_directories[isupercell] << endl;}
            found_all_eps_dat_files=false;
          }
        }
      }
    }
    if(!found_all_eps_dat_files){return;}

    //check eps can be extracted
    xPLASMONICS xplasm;
    for(i=0;i<m_vfilenames_plasm.size();i++){
      xplasm.filename=m_vfilenames_plasm[i];xplasm.getEPS();
      if(!aurostd::isfloat(xplasm.eps)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " eps cannot be extracted from "+m_vfilenames_plasm[i] << endl;}
        found_all_eps_dat_files=false;
      }
    }
    if(!found_all_eps_dat_files){return;}

    //sort m_veps_plasm, do manually as string sorting is different than double sorting
    uint j=0;
    double eps1=0.0,eps2=0.0;
    for(i=0;i<m_vfilenames_plasm.size();i++){
      for(j=i;j<m_vfilenames_plasm.size();j++){
        xplasm.filename=m_vfilenames_plasm[i];xplasm.getEPS();
        eps1=aurostd::string2utype<double>(xplasm.eps);
        xplasm.filename=m_vfilenames_plasm[j];xplasm.getEPS();
        eps2=aurostd::string2utype<double>(xplasm.eps);
        if(eps2<eps1){
          tmp_str=m_vfilenames_plasm[i];
          m_vfilenames_plasm[i]=m_vfilenames_plasm[j];
          m_vfilenames_plasm[j]=tmp_str;
        }
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_vfilenames_plasm=" << aurostd::joinWDelimiter(m_vfilenames_plasm,",") << endl;}

    m_veps_plasm.clear();
    m_vxplasm.clear();
    for(i=0;i<m_vfilenames_plasm.size();i++){
      xplasm.filename=m_vfilenames_plasm[i];xplasm.getEPS();
      m_veps_plasm.push_back(xplasm.eps);
      m_vxplasm.push_back(xPLASMONICS());
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_veps_plasm=" << aurostd::joinWDelimiter(m_veps_plasm,",") << endl;}
    setAvgPlasmonicData(temperature);
  }

  void POccCalculator::plotAvgDOSCAR(double temperature) const {return plotAvgDOSCAR(m_xdoscar,temperature,m_aflags.Directory);}
  void POccCalculator::plotAvgDOSCAR(const string& doscar_path,const string& directory) const {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    xDOSCAR xdos(doscar_path,*p_FileMESSAGE,true,*p_oss);
    double temperature=poccDOSCAR2temperature(doscar_path);
    return plotAvgDOSCAR(xdos,temperature,directory);
  }
  void POccCalculator::plotAvgDOSCAR(const xDOSCAR& xdos,double temperature,const string& directory) const {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}
    aurostd::xoption cmdline_opts, plot_opts;
    cmdline_opts.push_attached("PLOTTER::PRINT", "png");
    //ME20210927 - use carstring to distingish between types of DOSCARs
    if(xdos.carstring == "POCC") {  //turn off for ME - POCC+APL
      if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}
      cmdline_opts.push_attached("PLOT_DOS", directory);
      plot_opts = plotter::getPlotOptionsEStructure(cmdline_opts, "PLOT_DOS");
      plot_opts.push_attached("NORMALIZATION","ATOM");  //CO20211124
      plot_opts.push_attached("PROJECTION","ORBITALS");
      //plot_opts.push_attached("EXTENSION","dos_orbitals_T"+aurostd::utype2string(temperature,TEMPERATURE_PRECISION)+"K");
      plot_opts.push_attached("EXTENSION","dos_orbitals_T"+(*this).getTemperatureString(temperature)+"K");
      plotter::PLOT_DOS(plot_opts,xdos,*p_FileMESSAGE,*p_oss);

      //try species-projected
      plot_opts.push_attached("ARUN_DIRECTORY",m_ARUN_directories[0]);
      plot_opts.pop_attached("PROJECTION");
      plot_opts.push_attached("PROJECTION","SPECIES");
      plot_opts.pop_attached("EXTENSION");
      //plot_opts.push_attached("EXTENSION","dos_species_T"+aurostd::utype2string(temperature,TEMPERATURE_PRECISION)+"K");
      plot_opts.push_attached("EXTENSION","dos_species_T"+(*this).getTemperatureString(temperature)+"K");
      xstructure xstr=plotter::getStructureWithNames(plot_opts,*p_FileMESSAGE,xdos.carstring,*p_oss);
      if(xdos.vDOS.size()==xstr.atoms.size()){plotter::PLOT_DOS(plot_opts,xdos,*p_FileMESSAGE,*p_oss);}
      else{
        message << "Skipping species-projected plots, there are varying num_each_type among POCC.ARUN's" << endl;
        pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      }

      if(0){  //turn on for SG - plot orbitals for each species near fermi energy (dos_species_T0K_Cs_1)
        plot_opts.pop_attached("PROJECTION");
        plot_opts.push_attached("PROJECTION","ORBITALS");
        plot_opts.push_attached("DATATYPE","SPECIES");
        plot_opts.pop_attached("XMIN"); plot_opts.push_attached("XMIN","-1.5");
        plot_opts.pop_attached("XMAX"); plot_opts.push_attached("XMAX","1.5");
        for(uint ispecies=0;ispecies<xstr_pocc.species.size();ispecies++){
          plot_opts.pop_attached("DATASET"); plot_opts.push_attached("DATASET",aurostd::utype2string(ispecies+1));
          plot_opts.pop_attached("DATALABEL"); plot_opts.push_attached("DATALABEL",KBIN::VASP_PseudoPotential_CleanName(xstr_pocc.species[0]));
          plotter::PLOT_PDOS(plot_opts,xdos,*p_FileMESSAGE,*p_oss);
        }
      }
    } else if (xdos.carstring == "PHON") {
      cmdline_opts.push_attached("PLOT_PHDOS", directory);
      plot_opts = plotter::getPlotOptionsPhonons(cmdline_opts, "PLOT_PHDOS");
      plot_opts.push_attached("PROJECTION", "ATOMS");
      plot_opts.push_attached("EXTENSION", "phdos_T"+getTemperatureString(xdos.temperature)+"K");
      plot_opts.flag("PLOT_ALL_ATOMS", true);
      plot_opts.pop_attached("XMIN");plot_opts.pop_attached("XMAX");
      plotter::PLOT_PHDOS(plot_opts, xdos);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Invalid xDOSCAR format for POCC.", _RUNTIME_ERROR_);
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }

#define pocc_precision _DOUBLE_WRITE_PRECISION_ //12
#define pocc_roundoff_tol 5.0*pow(10,-((int)pocc_precision)-1)
  void POccCalculator::writeResults() const { //TEMPERATURE-INDEPENDENT
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    string enthalpy_tag="enthalpy"; //H
    unsigned long long int isupercell=0;

    //int pocc_precision=12;
    //double pocc_roundoff_tol=5.0*pow(10,-((int)pocc_precision)-1);

    //POCC.OUT
    message << "Writing out " << POCC_FILE_PREFIX+POCC_OUT_FILE << " (temperature-independent properties)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    stringstream pocc_out_ss;
    pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;
    pocc_out_ss << POCC_AFLOWIN_tag << "START_TEMPERATURE=ALL" << endl;  //"  (K)"
    //[CO20200502 - removed unnecessary separation line]pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;

    //supercell degeneracy
    isupercell=0;
    for(std::list<POccSuperCellSet>::const_iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      pocc_out_ss << "degeneracy_supercell_" << std::setfill('0') << std::setw(aurostd::getZeroPadding(l_supercell_sets.size())) << isupercell+1  << "=" << aurostd::utype2string((*it).getDegeneracy()) << "  [" << m_ARUN_directories[isupercell] << "]" << endl; //+1 so we start at 1, not 0 (count)
    }
    //supercell enthalpy_atom
    isupercell=0;
    for(std::list<POccSuperCellSet>::const_iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      pocc_out_ss << enthalpy_tag << "_atom_supercell_" << std::setfill('0') << std::setw(aurostd::getZeroPadding(l_supercell_sets.size())) << isupercell+1  << "=" << aurostd::utype2string((*it).m_energy_dft,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV/at)  " << "[" << m_ARUN_directories[isupercell] << "]";
      if(isupercell==m_ARUN_directory_ground){pocc_out_ss << "  [ground]";}
      pocc_out_ss << endl; //+1 so we start at 1, not 0 (count)
    }
    //if(m_energy_dft_ground!=AUROSTD_MAX_DOUBLE) pocc_out_ss << enthalpy_tag << "_atom_ground=" << aurostd::utype2string(m_energy_dft_ground,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV/at)  " << "[" << m_ARUN_directories[m_ARUN_directory_ground] << "]" << endl;
    if(m_Hmix!=AUROSTD_MAX_DOUBLE) pocc_out_ss << enthalpy_tag << "_mix_atom=" << aurostd::utype2string(m_Hmix,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV/at)" << endl;
    if(m_efa!=AUROSTD_MAX_DOUBLE) pocc_out_ss << "entropy_forming_ability=" << aurostd::utype2string(m_efa,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV/at)^{-1}" << endl;

    //[CO20200502 - removed unnecessary separation line]pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;
    pocc_out_ss << POCC_AFLOWIN_tag << "STOP_TEMPERATURE=ALL" << endl;  //"  (K)"
    pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;

    aurostd::stringstream2file(pocc_out_ss,getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_OUT_FILE);  //NO APPEND, write clean

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }

  void POccCalculator::writeResults(double temperature) const { //TEMPERATURE-DEPENDENT
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    //int pocc_precision=12;
    //double pocc_roundoff_tol=5.0*pow(10,-((int)pocc_precision)-1);

    //POCC.OUT
    message << "Writing out " << POCC_FILE_PREFIX+POCC_OUT_FILE << " (T=" << temperature << "K properties)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); //CO20200502 - no getTemperatureString(T) needed here
    stringstream pocc_out_ss;
    pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;
    pocc_out_ss << POCC_AFLOWIN_tag << "START_TEMPERATURE=" << (*this).getTemperatureString(temperature) << "_K" << endl;  //"  (K)"
    //[CO20200502 - removed unnecessary separation line]pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;

    //supercell probabilities
    unsigned long long int isupercell=0;
    for(std::list<POccSuperCellSet>::const_iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      //pocc_out_ss << "probability[" << m_ARUN_directories[isupercell] << "]=" << aurostd::utype2string((*it).m_probability,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << endl;
      pocc_out_ss << "probability_supercell_" << std::setfill('0') << std::setw(aurostd::getZeroPadding(l_supercell_sets.size())) << isupercell+1  << "=" << aurostd::utype2string((*it).m_probability,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  [" << m_ARUN_directories[isupercell] << "]" << endl; //+1 so we start at 1, not 0 (count)
    }
    
    if(m_rdf_all.rows!=0){
      pocc_out_ss << POCC_AFLOWIN_tag << "START_RDF" << endl;
      pocc_out_ss << RDF2string(xstr_pocc,m_rdf_rmax,m_rdf_nbins,m_rdf_all);
      pocc_out_ss << POCC_AFLOWIN_tag << "STOP_RDF" << endl;
    }

    //fix m_Egap and m_Egap_net for metals (Egap==0!)
    vector<double> Egap;Egap.assign(m_Egap.size(),0.0);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " m_Egap.size()=" << m_Egap.size() << endl;
      cerr << __AFLOW_FUNC__ << " Egap.size()=" << Egap.size() << endl;
    }
    double Egap_net=0.0;
    for(uint ispin=0;ispin<m_Egap.size();ispin++){
      if(m_Egap[ispin]!=_METALGAP_){Egap[ispin]=m_Egap[ispin];}
    }
    if(m_Egap_net!=_METALGAP_){Egap_net=m_Egap_net;}

    if(Egap.size()==2){
      pocc_out_ss << "Egap_spin_up=" << aurostd::utype2string(Egap[0],pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV)" << endl;
      pocc_out_ss << "Egap_spin_dn=" << aurostd::utype2string(Egap[1],pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV)" << endl;
    }
    if(Egap_net!=AUROSTD_MAX_DOUBLE) pocc_out_ss << "Egap_net=" << aurostd::utype2string(Egap_net,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM) << "  (eV)" << endl;
    //[CO20200502 - removed unnecessary separation line]pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;
    //PLASMONICS
    for(uint ieps=0;ieps<m_veps_plasm.size();ieps++){
      if(m_vxplasm[ieps].venergy.size()!=0){
        uint i=0,padding=25;
        pocc_out_ss << POCC_AFLOWIN_tag << "START_PLASMONICS_EPS_" << m_veps_plasm[ieps] << endl;
        //header
        pocc_out_ss << "#";
        pocc_out_ss << aurostd::PaddedPRE("energy (eV)",padding) << " ";
        pocc_out_ss << aurostd::PaddedPRE("EELS",padding) << " ";
        pocc_out_ss << aurostd::PaddedPRE("Re(dielectric)",padding) << " ";
        pocc_out_ss << aurostd::PaddedPRE("Im(dielectric)",padding) << " ";
        pocc_out_ss << endl;
        for(i=0;i<m_vxplasm[ieps].venergy.size();i++){
          pocc_out_ss << " "; //spacing
          pocc_out_ss << aurostd::PaddedPRE(aurostd::utype2string(m_vxplasm[ieps].venergy[i],pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM),padding) << " ";
          pocc_out_ss << aurostd::PaddedPRE(aurostd::utype2string(m_vxplasm[ieps].veels[i],pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM),padding) << " ";
          pocc_out_ss << aurostd::PaddedPRE(aurostd::utype2string(m_vxplasm[ieps].vdielectric[i].re,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM),padding) << " ";
          pocc_out_ss << aurostd::PaddedPRE(aurostd::utype2string(m_vxplasm[ieps].vdielectric[i].im,pocc_precision,true,pocc_roundoff_tol,SCIENTIFIC_STREAM),padding) << " ";
          pocc_out_ss << endl;
        }
        pocc_out_ss << POCC_AFLOWIN_tag << "STOP_PLASMONICS_EPS_" << m_veps_plasm[ieps] << endl;
      }
    }
    pocc_out_ss << POCC_AFLOWIN_tag << "STOP_TEMPERATURE=" << (*this).getTemperatureString(temperature) << "_K" << endl;  //"  (K)"
    pocc_out_ss << AFLOWIN_SEPARATION_LINE << endl;

    aurostd::stringstream2file(pocc_out_ss,getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_OUT_FILE,"APPEND");

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }

  bool POccCalculator::QMVASPsFound() const{
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()==0",_RUNTIME_ERROR_);}
    if(m_ARUN_directories.size()!=l_supercell_sets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()!=l_supercell_sets.size()",_RUNTIME_ERROR_);}

    bool found_all_QMVASPs=true;
    string qmvasp_file="";
    unsigned long long int isupercell=0;
    for(isupercell=0;isupercell<m_ARUN_directories.size()&&found_all_QMVASPs==true;isupercell++){
      qmvasp_file=getARUNDirectoryPath(isupercell)+"/"+DEFAULT_AFLOW_QMVASP_OUT;
      if(!aurostd::EFileExist(qmvasp_file)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << DEFAULT_AFLOW_QMVASP_OUT << " not found in "+m_ARUN_directories[isupercell] << endl;}
        found_all_QMVASPs=false;
      }
    }
    if(!found_all_QMVASPs){
      message << "Waiting for complete VASP calculations (" << qmvasp_file << " is missing)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
      return false;
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
    return true;
  }

  vector<double> getVTemperatures(const string& temp_string){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " temp_string=" << temp_string << endl;}
    if(temp_string.empty()){
      vector<double> tmp;return tmp;  //return something, it will try again
      //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"temp_string.empty()",_INPUT_ILLEGAL_);
    }

    vector<string> tokens;
    if(temp_string.find(":")!=string::npos){aurostd::string2tokens(temp_string,tokens,":");}
    else if(temp_string.find("-")!=string::npos){aurostd::string2tokens(temp_string,tokens,"-");}
    else if(temp_string.find(",")!=string::npos){aurostd::string2tokens(temp_string,tokens,",");}
    else{tokens.push_back(temp_string);}

    if(tokens.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown temp_string format",_INPUT_ILLEGAL_);}

    //check that they are all doubles
    vector<double> dtokens;
    for(uint i=0;i<tokens.size();i++){
      if(!aurostd::isfloat(tokens[i])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"temp_string token[i="+tokens[i]+"] is not a float",_INPUT_ILLEGAL_);} //check that all doubles
      dtokens.push_back(aurostd::string2utype<double>(tokens[i]));
      if(std::signbit(dtokens.back())){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Negative temperature found in temp_string: "+tokens[i],_INPUT_ILLEGAL_);}
    }

    //simple comma-separated list
    if(dtokens.size()==1 || temp_string.find(",")!=string::npos){return dtokens;}

    //treat ':' and '-' as delimiter for ranges
    if(tokens.size()<1 || tokens.size()>3){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown temp_string format",_INPUT_ILLEGAL_);}

    double interval=100.0;
    if(dtokens.size()==3){interval=dtokens.back();}
    vector<double> vtemperatures;
    for(double temp=dtokens[0];temp<=dtokens[1];temp+=interval){vtemperatures.push_back(temp);}
    if(vtemperatures.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No temperatures extracted from temp_string",_INPUT_ILLEGAL_);}

    if(LDEBUG){
      for(uint i=0;i<vtemperatures.size();i++){
        cerr << __AFLOW_FUNC__ << " vtemperatures[i=" << i << "]=" << vtemperatures[i] << endl;
      }
    }

    return vtemperatures;
  }

  void POccCalculator::CleanPostProcessing(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(LDEBUG){;} //keep busy

    vector<string> vfiles;
    aurostd::DirectoryLS(m_aflags.Directory,vfiles);

    for(uint i=0;i<vfiles.size();i++){
      //ME20211019 - Added APL options
      if (m_kflags.KBIN_PHONONS_CALCULATION_APL) {
        if(vfiles[i].find(POCC_PHDOSCAR_FILE)!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //PHDOSCAR.pocc_T0000K.xz
        if(vfiles[i].find(DEFAULT_APL_PHPOSCAR_FILE)!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //PHDOSCAR.xz
        if(vfiles[i].find(POCC_FILE_PREFIX+POCC_APL_OUT_FILE)!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //aflow.pocc.apl.out
      } else {
        if(vfiles[i].find(POCC_FILE_PREFIX+POCC_OUT_FILE)!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //aflow.pocc.out
        if(vfiles[i].find(POCC_DOSCAR_FILE)!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //DOSCAR.pocc_T0000K.xz
        if(vfiles[i].find("dos_orbitals_")!=string::npos && vfiles[i].find(".png")!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //SYSTEM_dos_orbitals_T2400K.png.xz
        if(vfiles[i].find("dos_species_")!=string::npos && vfiles[i].find(".png")!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //SYSTEM_dos_species_T2400K.png.xz
        if(vfiles[i].find("dos_atoms_")!=string::npos && vfiles[i].find(".png")!=string::npos){message << "Removing old postprocessing file: " << vfiles[i];pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);aurostd::RemoveFile(aurostd::CleanFileName(getOutputPath()+"/"+vfiles[i]));continue;}  //SYSTEM_dos_atoms_T2400K.png.xz
      }
    }
  }

  void POccCalculator::convolution(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    string inputs=m_p_flags.getattachedscheme("POCC::CONVOLUTION");
    vector<string> vinputs;
    aurostd::string2tokens(inputs,vinputs,",");

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
      cerr << __AFLOW_FUNC__ << " vinputs=" << aurostd::joinWDelimiter(vinputs,",") << endl;
    }

    if(vinputs.size()<2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"vinputs.size()<2",_INPUT_NUMBER_);}

    uint idir=0,i=0;
    POccCalculator pcalc;
    _aflags aflags;
    l_supercell_sets.clear();
    m_ARUN_directories.clear();
    vector<string> vsystems;
    for(idir=0;idir<vinputs.size();idir++){
      aflags.Directory=vinputs[idir];
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " analyzing directory=" << aflags.Directory << endl;}
      if(!pocc::structuresGenerated(aflags.Directory)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,vinputs[idir]+" is not a POCC directory",_INPUT_ILLEGAL_);}
      pcalc.initialize(aflags,*p_FileMESSAGE,*p_oss);
      if(!pcalc.m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccCalculator failed to initialize");}
      vsystems.push_back(pcalc.m_vflags.AFLOW_SYSTEM.content_string);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " system_name=" << vsystems.back() << endl;}
      //CO20211130 - come back, need to define xstr_pocc for xDOSCAR
      //requires knowing the convolved concentration
      //will solve later, for now, just load first pcalc.xstr_pocc
      if(idir==0){xstr_pocc=pcalc.xstr_pocc;}
      pcalc.loadDataIntoCalculator();
      if(!pcalc.QMVASPsFound()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"QMVASP files not found");}
      //add to pcalc 
      for(std::list<POccSuperCellSet>::const_iterator it=pcalc.l_supercell_sets.begin();it!=pcalc.l_supercell_sets.end();++it){l_supercell_sets.push_back(*it);}
      for(i=0;i<pcalc.m_ARUN_directories.size();i++){
        m_ARUN_directories.push_back(pcalc.m_aflags.Directory+"/"+pcalc.m_ARUN_directories[i]);
        message << "adding m_ARUN_directories[" << m_ARUN_directories.size()-1 << "]=" << m_ARUN_directories.back();pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      }
    }
    if(m_ARUN_directories.size()!=l_supercell_sets.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_ARUN_directories.size()!=l_supercell_sets.size()",_RUNTIME_ERROR_);}
    m_initialized=true;

    bool found_pocc=false;
    aurostd::xoption pocc_settings;
    string pps="",proto="",pocc_params="",pocc_tol="",pocc_arun="",module_arun="";
    string pps0="",proto0="",pocc_params0="";
    uint natoms0=pcalc.xstr_pocc.atoms.size();
    bool same_pps_proto=true;
    xstructure xstr;
    vector<xstructure> vxstr;
    for(i=0;i<vsystems.size();i++){
      found_pocc=pflow::POccInputs2Xstr(vsystems[i],pocc_settings,xstr,*p_FileMESSAGE,*p_oss);
      if(!found_pocc){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to process system name: "+vsystems[i]);}
      vxstr.push_back(xstr);
      pps=pocc_settings.getattachedscheme("PPS");
      proto=pocc_settings.getattachedscheme("PROTO");
      pocc_params=pocc_settings.getattachedscheme("POCC_PARAMS");
      pocc_tol=pocc_settings.getattachedscheme("POCC_TOL");
      pocc_arun=pocc_settings.getattachedscheme("POCC_ARUN");
      module_arun=pocc_settings.getattachedscheme("MODULE_ARUN");
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " xstr=" << endl << xstr;
        cerr << __AFLOW_FUNC__ << " pps=" << pps << endl;
        cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;
        cerr << __AFLOW_FUNC__ << " pocc_params=" << pocc_params << endl;
        cerr << __AFLOW_FUNC__ << " pocc_tol=" << pocc_tol << endl;
        cerr << __AFLOW_FUNC__ << " pocc_arun=" << pocc_arun << endl;
        cerr << __AFLOW_FUNC__ << " module_arun=" << module_arun << endl;
      }
      if(i==0){
        pps0=pps;
        proto0=proto;
        pocc_params0=pocc_params;
      }
      if(pps!=pps0||proto0!=proto||natoms0!=xstr.atoms.size()){same_pps_proto=false;}
    }
    if(same_pps_proto==false){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"found mismatch in pps or proto or natoms");}

    //same pps and proto mean they have the same number of atoms, checked anyway above
    uint iatom=0;
    deque<_atom> datoms;
    for(iatom=0;iatom<natoms0;iatom++){datoms.push_back(xstr_pocc.atoms[iatom]);datoms.back().partial_occupation_value=0.0;}  //add and reset
    for(iatom=0;iatom<natoms0;iatom++){
      for(i=0;i<vxstr.size();i++){datoms[iatom].partial_occupation_value+=(1.0/vxstr.size())*vxstr[i].atoms[iatom].partial_occupation_value;}
    }
    xstr_pocc.ReplaceAtoms(datoms,false); //do NOT check overlap (POCC)
    xstr_pocc.SetAutoVolume();
    xstr_pocc.ReScale(std::pow(xstr_pocc.GetVolume()/GetVol(vxstr.front().lattice),1.0/3.0));
    xstr_pocc.neg_scale=true;

    vector<POccUnit> pocc_sites=getPOccSites(xstr_pocc,m_aflags,*p_FileMESSAGE,*p_oss);
    string pocc_params_new="";
    stringstream ss_tmp;
    uint isite=0,iocc=0;
    for(isite=0;isite<pocc_sites.size();isite++){
      if(isite>0){pocc_params_new+="_";}
      pocc_params_new+="P"+aurostd::utype2string(isite);
      for(iocc=0;iocc<pocc_sites[isite].v_occupants.size();iocc++){
        pocc_params_new+="-"+aurostd::utype2string(xstr_pocc.atoms[pocc_sites[isite].v_occupants[iocc]].partial_occupation_value);
        ss_tmp << char('A'+pocc_sites[isite].v_types[iocc]);
        pocc_params_new+="x"+ss_tmp.str();
        aurostd::StringstreamClean(ss_tmp);
      }
    }
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " pocc_params0=" << pocc_params0 << endl;
      cerr << __AFLOW_FUNC__ << " pocc_params_new(pre)=" << pocc_params_new << endl;
    }

    xstr_nopocc=createNonPOccStructure(xstr_pocc,pocc_sites);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr_nopocc=" << endl << xstr_nopocc;}
    pflow::FIX_POCC_PARAMS(xstr_nopocc,pocc_params_new);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_params_new(post)=" << pocc_params_new << endl;}
    aurostd::StringSubst(xstr_pocc.title,pocc_params0,pocc_params_new);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr_pocc=" << endl << xstr_pocc;}

    m_vflags.AFLOW_SYSTEM.isentry=true;
    m_vflags.AFLOW_SYSTEM.content_string=vsystems.front();
    aurostd::StringSubst(m_vflags.AFLOW_SYSTEM.content_string,pocc_params0,pocc_params_new);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_vflags.AFLOW_SYSTEM.content_string=" << m_vflags.AFLOW_SYSTEM.content_string << endl;}

    setPOccStructure(xstr_pocc);
    m_convolution=true;
    writePARTCAR();
    postProcessing();
  }

  void POccCalculator::loadDataIntoCalculator(){
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccCalculator failed to initialize");}

    POccStructuresFile psf(*p_FileMESSAGE,*p_oss);
    psf.initialize(POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,m_aflags);
    if(!psf.loadDataIntoCalculator((*this),false)){  //do not try DirectoryLS() yet
      //if we get here, then postprocessing file is TOO old, try rewriting after we read in ALL_STRUCTURES file
      psf.initialize(POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE,m_aflags);
      if(!psf.loadDataIntoCalculator((*this),true)){  //try DirectoryLS() here
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot load data from structures files",_FILE_NOT_FOUND_);
      }
      //move old file and write new one
    }
    //StructuresAllFile2SupercellSets();  //this also works, just slower
    //StructuresUniqueFile2SupercellSets();


    //skip files if you need to
    vector<string> v_aruns2skip;
    if(!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()){aurostd::string2tokens(m_kflags.KBIN_POCC_ARUNS2SKIP_STRING,v_aruns2skip,",");}
    if(XHOST.vflag_control.flag("ARUNS2SKIP")){aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("ARUNS2SKIP"),v_aruns2skip,",");}
    if(v_aruns2skip.size()){
      string index_str="";
      uint index=0;
      vector<uint> v_indices;
      uint i=0;
      for(i=0;i<v_aruns2skip.size();i++){
        string index_str=v_aruns2skip[i];
        if(!aurostd::isfloat(index_str)){  //input looks likes ARUN.POCC_..., convert to 1,2,3...
          aurostd::StringSubst(index_str,"ARUN.POCC_","");
          vector<string> tokens;
          aurostd::string2tokens(index_str,tokens,"_");
          index_str=tokens[0];
          if(!aurostd::isfloat(index_str)){  //non-recognizable input
            pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown format for ARUN2SKIP[i="+aurostd::utype2string(i)+"]="+v_aruns2skip[i],m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_); 
            continue;
          }
        }
        //now index_str should be good - assume index_str = index_str+1
        index=aurostd::string2utype<uint>(index_str);
        if((index-1)>=m_ARUN_directories.size()){
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"ARUN2SKIP[i="+aurostd::utype2string(i)+"]="+v_aruns2skip[i]+" is out of range",m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_); 
          continue;
        }
        v_indices.push_back(index-1);
      }
      if(v_indices.size()){
        std::sort(v_indices.rbegin(),v_indices.rend());  //sort backwards so we erase
        std::list<POccSuperCellSet>::iterator it;
        for(i=0;i<v_indices.size();i++){
          pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,"Ignoring "+m_ARUN_directories[v_indices[i]],m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_); 
          m_ARUN_directories.erase(m_ARUN_directories.begin()+v_indices[i]);

          //erase
          it=l_supercell_sets.begin();std::advance(it,v_indices[i]);
          l_supercell_sets.erase(it);
        }
      }
    }
  }

  //CT20200319 - added AEL/AGL option
  void POccCalculator::setTemperatureStringParameters(){vector<double> v_temperatures;return setTemperatureStringParameters(v_temperatures);}
  void POccCalculator::setTemperatureStringParameters(vector<double>& v_temperatures){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccCalculator failed to initialized");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    //START: TEMPERATURE DEPENDENT PROPERTIES
    v_temperatures.clear();
    //v_temperatures.push_back(300);  //1000
    v_temperatures=getVTemperatures(m_kflags.KBIN_POCC_TEMPERATURE_STRING);
    if(XHOST.vflag_control.flag("CALCULATION_TEMPERATURE")){v_temperatures.clear();v_temperatures=getVTemperatures(XHOST.vflag_control.getattachedscheme("CALCULATION_TEMPERATURE"));}  //command line input
    if(v_temperatures.empty()){v_temperatures.clear();v_temperatures=getVTemperatures(DEFAULT_POCC_TEMPERATURE_STRING);}  //user aflow.rc
    if(v_temperatures.empty()){v_temperatures.clear();v_temperatures=getVTemperatures(AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING);}  //internal aflow.rc, will always work
    if(v_temperatures.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"v_temperatures.empty()",_INPUT_ILLEGAL_);}

    //get zero-padding, since temperature cannot be negative, just get max
    getTemperatureStringParameters(v_temperatures,m_temperature_precision,m_temperatures_int,m_zero_padding_temperature);
  }

  void POccCalculator::postProcessing(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccCalculator failed to initialize");}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    if(m_ARUN_directories.size() == 0) loadDataIntoCalculator();
    if(!QMVASPsFound()){return;}
    setDFTEnergies();
    setEFA();

    vector<double> v_temperatures;
    setTemperatureStringParameters(v_temperatures);

    message << "Performing POCC post-processing for these temperatures: " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(v_temperatures,5),",");pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    if(1 && !m_kflags.KBIN_PHONONS_CALCULATION_APL){  //turn off slow pocc-postprocessing  //ME20211019 - do not run regular POCC postprocessing for APL or LIB2LIB will run this twice
      aurostd::RemoveFile(getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_OUT_FILE); //clear file
      writeResults(); //write temperature-independent properties first
      for(uint itemp=0;itemp<v_temperatures.size();itemp++){
        calculateRELAXProperties(v_temperatures[itemp]);
        calculateSTATICProperties(v_temperatures[itemp]);
        calculatePlasmonicProperties(v_temperatures[itemp]);
        writeResults(v_temperatures[itemp]);  //write temperature-dependent properties next
      }
    }
    if (m_kflags.KBIN_PHONONS_CALCULATION_AEL) { //CT20200319
      if(LDEBUG){cerr << __AFLOW_FUNC__ << "Running AEL postprocessing" << endl;}
      calculateElasticProperties(v_temperatures);
    }
    if (m_kflags.KBIN_PHONONS_CALCULATION_AGL) { //CT20200323
      if(LDEBUG){cerr << "Running AGL postprocessing" << endl;}
      calculateDebyeThermalProperties(v_temperatures);
    }
    if (m_kflags.KBIN_PHONONS_CALCULATION_APL) { //ME20210927
      if(LDEBUG){cerr << "Running APL postprocessing" << endl;}
      calculatePhononPropertiesAPL(v_temperatures);
    }
    if (m_kflags.KBIN_PHONONS_CALCULATION_QHA) { //AS20210402
      if(LDEBUG){cerr << "Running QHA postprocessing" << endl;}
      calculateQHAProperties();
      calculateQHAPropertiesAVG(v_temperatures);
    }

    //END: TEMPERATURE DEPENDENT PROPERTIES

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
  }
} // namespace pocc

namespace pocc {
  ////--------------------------------------------------------------------------------
  //// class POccGroup
  ////--------------------------------------------------------------------------------
  ////--------------------------------------------------------------------------------
  //// constructor
  ////--------------------------------------------------------------------------------
  //POccGroup::POccGroup(ostream& oss) : xStream(),m_initialized(false) {initialize(oss);}
  //POccGroup::POccGroup(ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(FileMESSAGE,oss);}
  //POccGroup::POccGroup(const _aflags& aflags,ostream& oss) : xStream(),m_initialized(false) {initialize(aflags,oss);}
  //POccGroup::POccGroup(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(aflags,FileMESSAGE,oss);}
  //POccGroup::POccGroup(const POccGroup& b) {copy(b);} // copy PUBLIC

  //POccGroup::~POccGroup() {xStream::free();free();}

  //const POccGroup& POccGroup::operator=(const POccGroup& b) {  // operator= PUBLIC
  //if(this!=&b) {copy(b);}
  //return *this;
  //}
  //bool POccGroup::operator<(const POccGroup& other) const {return partial_occupation_value<other.partial_occupation_value;}

  //void POccGroup::clear() {POccGroup a;copy(a);}
  //void POccGroup::free() {
  //m_initialized=false;
  //m_aflags.clear();
  //site=AUROSTD_MAX_UINT;
  //partial_occupation_value=AUROSTD_MAX_DOUBLE;
  //v_occupants.clear();
  //v_types.clear();
  //}

  //void POccGroup::copy(const POccGroup& b) { // copy PRIVATE
  //xStream::copy(b);
  //m_initialized=b.m_initialized;
  //m_aflags=b.m_aflags;
  //site=b.site;
  //partial_occupation_value=b.partial_occupation_value;
  //v_occupants.clear();for(uint i=0;i<b.v_occupants.size();i++){v_occupants.push_back(b.v_occupants[i]);}
  //v_types.clear();for(uint i=0;i<b.v_types.size();i++){v_types.push_back(b.v_types[i]);}
  //}

  //void POccGroup::initialize(ostream& oss) {
  //xStream::free();
  //ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
  //initialize(*_p_FileMESSAGE,oss);
  //f_new_ofstream=true;  //override
  //}
  //void POccGroup::initialize(ofstream& FileMESSAGE,ostream& oss) {
  //free();
  //try{
  //setOFStream(FileMESSAGE); f_new_ofstream=false;
  //setOSS(oss);
  //m_initialized=false;  //no point
  //}
  //catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  //}

  //void POccGroup::initialize(const _aflags& aflags,ostream& oss) {
  //xStream::free();
  //ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
  //initialize(aflags,*_p_FileMESSAGE,oss);
  //f_new_ofstream=true;  //override
  //}
  //void POccGroup::initialize(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
  //free();
  //try{
  //setOFStream(FileMESSAGE); f_new_ofstream=false;
  //setOSS(oss);
  //setAFlags(aflags);
  //m_initialized=false;  //no point
  //}
  //catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  //}

  //void POccGroup::setAFlags(const _aflags& aflags) {m_aflags=aflags;}
} // namespace pocc

namespace pocc {
  bool SiteVacancyCount::operator<(const SiteVacancyCount& other) const {return vacancy_count<other.vacancy_count;}
} // namespace pocc

namespace pocc {
  POccSuperCell::POccSuperCell() {free();}
  POccSuperCell::POccSuperCell(const POccSuperCell& b) {copy(b);}
  POccSuperCell::~POccSuperCell() {free();}
  void POccSuperCell::free() {
    m_hnf_index=AUROSTD_MAX_ULLINT;
    m_site_config_index=AUROSTD_MAX_ULLINT;
    m_degeneracy=0; //VERY IMPORTANT, 0 m_degeneracy will kill division
    m_energy_uff=AUROSTD_MAX_DOUBLE;
  }
  void POccSuperCell::copy(const POccSuperCell& b){
    m_hnf_index=b.m_hnf_index;
    m_site_config_index=b.m_site_config_index;
    m_degeneracy=b.m_degeneracy;
    m_energy_uff=b.m_energy_uff;
  }
  void POccSuperCell::clear() {POccSuperCell b;copy(b);}
  const POccSuperCell& POccSuperCell::operator=(const POccSuperCell& b){
    if(this!=&b){copy(b);}
    return *this;
  }
  bool POccSuperCell::operator<(const POccSuperCell& other) const {
    if(m_hnf_index!=other.m_hnf_index){return m_hnf_index<other.m_hnf_index;}
    return m_site_config_index<other.m_site_config_index;
    //if(m_site_config_index!=other.m_site_config_index){return m_site_config_index<other.m_site_config_index;} //if hnf + site_config is same, it's the same structure!
    //return energy<other.energy; //not smart to sort by double as last criteria (consider precision)
  }
  bool sortPSCsUFFEnergy(const POccSuperCell& a, const POccSuperCell& b){return a.m_energy_uff<b.m_energy_uff;}
} // namespace pocc

namespace pocc {
  POccSuperCellSet::POccSuperCellSet() {free();}
  POccSuperCellSet::POccSuperCellSet(const POccSuperCellSet& b) {copy(b);}
  POccSuperCellSet::~POccSuperCellSet() {free();}
  void POccSuperCellSet::free() {
    m_psc_set.clear();
    m_energy_dft=AUROSTD_MAX_DOUBLE;
    m_probability=AUROSTD_MAX_DOUBLE;
  }
  void POccSuperCellSet::copy(const POccSuperCellSet& b){
    m_psc_set.clear();for(uint i=0;i<b.m_psc_set.size();i++){m_psc_set.push_back(b.m_psc_set[i]);}
    m_energy_dft=b.m_energy_dft;
    m_probability=b.m_probability;
  }
  void POccSuperCellSet::clear() {free();}
  const POccSuperCellSet& POccSuperCellSet::operator=(const POccSuperCellSet& b){
    if(this!=&b){copy(b);}
    return *this;
  }
  bool POccSuperCellSet::operator<(const POccSuperCellSet& other) const {
    if(getHNFIndex()!=other.getHNFIndex()){return getHNFIndex()<other.getHNFIndex();}
    return getSiteConfigIndex()<other.getSiteConfigIndex();
    //if(getSiteConfigIndex()!=other.getSiteConfigIndex()){return getSiteConfigIndex()<other.getSiteConfigIndex();} //if hnf + site_config is same, it's the same structure!
    //return getUFFEnergy()<other.getUFFEnergy(); //not smart to sort by double as last criteria (consider precision)
  }
  unsigned long long int POccSuperCellSet::getDegeneracy() const {
    unsigned long long int degeneracy=0;
    for(uint i=0;i<m_psc_set.size();i++){degeneracy+=m_psc_set[i].m_degeneracy;}
    if(degeneracy==0){throw aurostd::xerror(__AFLOW_FILE__,"POccSuperCellSet::getDegeneracy():","degeneracy==0",_RUNTIME_ERROR_);}
    return degeneracy;
    //return m_psc_set.size();
  }
  const POccSuperCell& POccSuperCellSet::getSuperCell() const {
    if(m_psc_set.size()==0){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_psc_set.size()==0",_INDEX_ERROR_);
    }
    return m_psc_set[0];
  }
  double POccSuperCellSet::getHNFIndex() const {return getSuperCell().m_hnf_index;} //fix me eventually to do robust matrix comparison
  double POccSuperCellSet::getSiteConfigIndex() const {return getSuperCell().m_site_config_index;} //fix me eventually to do robust matrix comparison
  double POccSuperCellSet::getUFFEnergy() const {return getSuperCell().m_energy_uff;}
} // namespace pocc

namespace pocc {
  UFFParamAtom::UFFParamAtom() {free();}
  UFFParamAtom::UFFParamAtom(const UFFParamAtom& b) {copy(b);}
  UFFParamAtom::~UFFParamAtom() {free();}
  void UFFParamAtom::free() {
    symbol.clear();
    r1=AUROSTD_MAX_DOUBLE;
    theta0=AUROSTD_MAX_DOUBLE;
    x1=AUROSTD_MAX_DOUBLE;
    D1=AUROSTD_MAX_DOUBLE;
    zeta=AUROSTD_MAX_DOUBLE;
    Z1=AUROSTD_MAX_DOUBLE;
    Vi=AUROSTD_MAX_DOUBLE;
    Uj=AUROSTD_MAX_DOUBLE;
    ChiI=AUROSTD_MAX_DOUBLE;
    hard=AUROSTD_MAX_DOUBLE;
    radius=AUROSTD_MAX_DOUBLE;
  }
  void UFFParamAtom::copy(const UFFParamAtom& b){
    symbol=b.symbol;
    r1=b.r1;
    theta0=b.theta0;
    x1=b.x1;
    D1=b.D1;
    zeta=b.zeta;
    Z1=b.Z1;
    Vi=b.Vi;
    Uj=b.Uj;
    ChiI=b.ChiI;
    hard=b.hard;
    radius=b.radius;
  }
  void UFFParamAtom::clear() {free();}
  const UFFParamAtom& UFFParamAtom::operator=(const UFFParamAtom& b){
    if(this!=&b){copy(b);}
    return *this;
  }
} // namespace pocc

namespace pocc {
  const UFFParamBond& UFFParamBond::operator=(const UFFParamBond& b){
    if(this!=&b){
      ren=b.ren;
      R0=b.R0;
      Kij=b.Kij;
      Xij=b.Xij;
      Dij=b.Dij;
      delta=b.delta;
      X6=b.X6;
      X12=b.X12;
    }
    return *this;
  }

  void UFFParamBond::calculate(UFFParamAtom& uffp1,UFFParamAtom& uffp2,double distij){
    //uffp1=_uffp1; uffp2=_uffp2; distij=_distij;
    //cerr << "SYMBOL " << _uffp1.symbol << " " << uffp1.r1 << "  " << _uffp2.symbol << " " << uffp2.r1 << endl;
    //SKIP equation 3 - zero for single bond (rbo)
    //equation 4
    ren = uffp1.r1 * uffp2.r1 * (pow( sqrt(uffp1.ChiI) - sqrt(uffp2.ChiI),2.0)) / (uffp1.ChiI*uffp1.r1 + uffp2.ChiI*uffp2.r1);
    //equation 2
    //NOTE: See http://towhee.sourceforge.net/forcefields/uff.html: there is a typo in the published paper
    R0 = uffp1.r1 + uffp2.r1 - ren;
    //equation 6
    Kij = 664.12 * uffp1.Z1 * uffp2.Z1 / (R0 * R0 * R0);
    Xij = sqrt(uffp1.x1 * uffp2.x1);
    Dij = sqrt(uffp1.D1 * uffp2.D1);
    delta = distij - R0;
    //cerr << "x1 " << uffp1.x1 << endl;
    //cerr << "x2 " << uffp2.x1 << endl;
    //cerr << "Xij " << Xij << endl;
    //cerr << "distij " << distij << endl;
    //cerr << "Xij/distij " << Xij/distij << endl;
    X6 = pow(Xij/distij,6);
    X12 = pow(Xij/distij,12);
    //cerr << "Dij " << Dij << endl;
    //cerr << "X12 " << X12<< endl;
    //cerr << "X6 " << X6<< endl;
  }

} // namespace pocc

namespace pocc {
  //--------------------------------------------------------------------------------
  // class POccUnit
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  POccUnit::POccUnit(ostream& oss) : xStream(oss),m_initialized(false) {initialize();}
  POccUnit::POccUnit(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
  POccUnit::POccUnit(const _aflags& aflags,ostream& oss) : xStream(oss),m_initialized(false) {initialize(aflags);}
  POccUnit::POccUnit(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags);}
  POccUnit::POccUnit(const POccUnit& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  POccUnit::~POccUnit() {xStream::free();free();}

  const POccUnit& POccUnit::operator=(const POccUnit& b) {  // operator= PUBLIC
    if(this!=&b) {copy(b);}
    return *this;
  }

  //   bool POccUnit::operator<(const POccUnit& other) const {
  //   if(v_occupants.size()!=other.v_occupants.size()){
  //   if(equivalent!=other.equivalent){
  //   return v_occupants.size()>other.v_occupants.size(); //sort in ascending order
  //   }
  //   //return equivalent<other.equivalent;
  //   }
  //   if(equivalent!=other.equivalent){return equivalent<other.equivalent;}
  //   return site<other.site;
  //   }

  void POccUnit::clear() {POccUnit a;copy(a);}
  void POccUnit::free() {
    m_initialized=false;
    m_aflags.clear();
    site=AUROSTD_MAX_UINT;
    partial_occupation_flag=false;
    partial_occupation_value=AUROSTD_MAX_DOUBLE;
    v_occupants.clear();
    v_types.clear();
    m_pocc_groups.clear();
    is_inequivalent=false;
    equivalent=AUROSTD_MAX_UINT;
  }

  void POccUnit::copy(const POccUnit& b) { // copy PRIVATE
    xStream::copy(b);
    m_initialized=b.m_initialized;
    m_aflags=b.m_aflags;
    site=b.site;
    partial_occupation_flag=b.partial_occupation_flag;
    partial_occupation_value=b.partial_occupation_value;
    v_occupants.clear();for(uint i=0;i<b.v_occupants.size();i++){v_occupants.push_back(b.v_occupants[i]);}
    v_types.clear();for(uint i=0;i<b.v_types.size();i++){v_types.push_back(b.v_types[i]);}
    m_pocc_groups.clear();for(uint i=0;i<b.m_pocc_groups.size();i++){m_pocc_groups.push_back(b.m_pocc_groups[i]);}
    is_inequivalent=b.is_inequivalent;
    equivalent=b.equivalent;
  }

  bool POccUnit::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }

  bool POccUnit::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }

  bool POccUnit::initialize() {
    free();
    m_initialized=false;  //no point
    return m_initialized;
  }

  bool POccUnit::initialize(const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(aflags);
  }

  bool POccUnit::initialize(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags);
  }

  bool POccUnit::initialize(const _aflags& aflags) {
    free();
    setAFlags(aflags);
    m_initialized=false;  //no point
    return m_initialized;
  }

  void POccUnit::setAFlags(const _aflags& aflags) {m_aflags=aflags;}

} // namespace pocc

namespace pocc {
  //--------------------------------------------------------------------------------
  // class POccCalculator
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  POccCalculator::POccCalculator(ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize();}
  POccCalculator::POccCalculator(const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(aflags);}
  POccCalculator::POccCalculator(const _aflags& aflags,const _kflags& kflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(aflags,kflags);}
  POccCalculator::POccCalculator(const _aflags& aflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(aflags,vflags);}
  POccCalculator::POccCalculator(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,aflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,aflags,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,aflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,aflags,kflags,vflags);}
  POccCalculator::POccCalculator(ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
  POccCalculator::POccCalculator(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags);}
  POccCalculator::POccCalculator(const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags,kflags);}
  POccCalculator::POccCalculator(const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags,vflags);}
  POccCalculator::POccCalculator(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,aflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,aflags,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,aflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags,aflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags,aflags,kflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags,aflags,vflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags,aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags,aflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags,aflags,kflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags,aflags,vflags);}
  POccCalculator::POccCalculator(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags,aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,kflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,kflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,vflags);}
  POccCalculator::POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,kflags,vflags);}
  POccCalculator::POccCalculator(const POccCalculator& b) : POccCalculatorTemplate(),xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  POccCalculator::~POccCalculator() {xStream::free();free();}

  const POccCalculator& POccCalculator::operator=(const POccCalculator& b) {  // operator= PUBLIC
    if(this!=&b) {copy(b);}
    return *this;
  }

  void POccCalculator::clear() {POccCalculator a;copy(a);}
  void POccCalculator::free() {
    POccCalculatorTemplate::free();
    m_initialized=false;
    m_p_flags.clear();
    m_aflags.clear();
    m_kflags.clear();
    m_vflags.clear();
    //p_str.clear();
    energy_analyzer.clear();
    resetHNFMatrices();
    resetSiteConfigurations();
    hnf_count=0;
    types_config_permutations_count=0;
    total_permutations_count=0;
    //for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){(*it).clear();}
    l_supercell_sets.clear();
    m_convolution=false;
    m_count_unique_fast=false;
    m_ARUN_directories.clear();
    m_Hmix=AUROSTD_MAX_DOUBLE;
    m_efa=AUROSTD_MAX_DOUBLE;
    m_temperature_precision=TEMPERATURE_PRECISION;
    m_zero_padding_temperature=0;
    m_temperatures_int=false;
    m_relaxation_max=AUROSTD_MAX_UINT;
    m_energy_dft_ground=AUROSTD_MAX_DOUBLE;
    m_ARUN_directory_ground=AUROSTD_MAX_UINT;
    m_rdf_all.clear();
    m_rdf_rmax=AUROSTD_MAX_DOUBLE;
    m_rdf_nbins=AUROSTD_MAX_INT;
    m_xdoscar.clear();
    m_Egap_DOS.clear();
    m_Egap.clear();
    m_Egap_DOS_net=AUROSTD_MAX_DOUBLE;
    m_Egap_net=AUROSTD_MAX_DOUBLE;
    m_vfilenames_plasm.clear();
    m_veps_plasm.clear();
    m_vxplasm.clear();
    enumerator_mode.clear();

    m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
  }

  void POccCalculator::copy(const POccCalculator& b) { // copy PRIVATE
    POccCalculatorTemplate::copy(b);
    xStream::copy(b);
    m_initialized=b.m_initialized;
    m_p_flags=b.m_p_flags;
    m_aflags=b.m_aflags;
    m_kflags=b.m_kflags;
    m_vflags=b.m_vflags;
    //p_str=b.p_str;
    energy_analyzer=b.energy_analyzer;
    hnf_mat=b.hnf_mat;
    a_start=b.a_start;b_start=b.b_start;c_start=b.c_start;
    d_start=b.d_start;e_start=b.e_start;f_start=b.f_start;
    v_unique_superlattices=b.v_unique_superlattices;
    v_types_config=b.v_types_config;
    config_iterator=b.config_iterator;
    v_config_order=b.v_config_order;
    m_energy_uff_tolerance=b.m_energy_uff_tolerance;
    hnf_count=b.hnf_count;
    types_config_permutations_count=b.types_config_permutations_count;
    total_permutations_count=b.total_permutations_count;
    //for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){(*it).clear();}
    l_supercell_sets.clear();for(std::list<POccSuperCellSet>::const_iterator it=b.l_supercell_sets.begin();it!=b.l_supercell_sets.end();++it){l_supercell_sets.push_back(*it);}
    m_convolution=b.m_convolution;
    m_count_unique_fast=b.m_count_unique_fast;
    m_Hmix=b.m_Hmix;
    m_efa=b.m_efa;
    m_temperature_precision=b.m_temperature_precision;
    m_zero_padding_temperature=b.m_zero_padding_temperature;
    m_temperatures_int=b.m_temperatures_int;
    m_relaxation_max=b.m_relaxation_max;
    m_energy_dft_ground=b.m_energy_dft_ground;
    m_convolution=b.m_convolution;
    m_ARUN_directory_ground=b.m_ARUN_directory_ground;
    m_ARUN_directories.clear();for(uint i=0;i<b.m_ARUN_directories.size();i++){m_ARUN_directories.push_back(b.m_ARUN_directories[i]);}
    m_rdf_all=b.m_rdf_all;
    m_rdf_rmax=b.m_rdf_rmax;
    m_rdf_nbins=b.m_rdf_nbins;
    m_xdoscar=b.m_xdoscar;
    m_Egap_DOS.clear();for(uint ispin=0;ispin<b.m_Egap_DOS.size();ispin++){m_Egap_DOS.push_back(b.m_Egap_DOS[ispin]);}
    m_Egap.clear();for(uint ispin=0;ispin<b.m_Egap.size();ispin++){m_Egap.push_back(b.m_Egap[ispin]);}
    m_Egap_DOS_net=b.m_Egap_DOS_net;
    m_Egap_net=b.m_Egap_net;
    m_vfilenames_plasm.clear();for(uint i=0;i<b.m_vfilenames_plasm.size();i++){m_vfilenames_plasm.push_back(b.m_vfilenames_plasm[i]);}
    m_veps_plasm.clear();for(uint i=0;i<b.m_veps_plasm.size();i++){m_veps_plasm.push_back(b.m_veps_plasm[i]);}
    m_vxplasm.clear();for(uint i=0;i<b.m_vxplasm.size();i++){m_vxplasm.push_back(b.m_vxplasm[i]);}
    enumerator_mode=b.enumerator_mode;
  }

  bool POccCalculator::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }
  bool POccCalculator::initialize(const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(aflags);
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _kflags& kflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(aflags,kflags);
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(aflags,vflags);
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,aflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,aflags,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,aflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }
  bool POccCalculator::initialize(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags);
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags,kflags);
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags,vflags);
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,aflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,aflags,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,aflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,aflags,kflags,vflags);
  }
  bool POccCalculator::initialize() {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const _aflags& aflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _kflags& kflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setVFlags(vflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      setVFlags(vflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _kflags& kflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setKFlags(kflags);
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(pocc_flags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(pocc_flags,aflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(pocc_flags,aflags,kflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(pocc_flags,aflags,vflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(pocc_flags,aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,aflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,aflags,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,aflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(xstr_pocc,pocc_flags,aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags,aflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags,aflags,kflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags,aflags,vflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags,aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,aflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,kflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,aflags,kflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,aflags,vflags);
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,pocc_flags,aflags,kflags,vflags);
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setPOccFlags(pocc_flags);
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setVFlags(vflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      setVFlags(vflags);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",true);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setPOccStructure(xstr_pocc);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setKFlags(kflags);
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",true);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",true);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccCalculator::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags) {
    free();
    try{
      m_energy_uff_tolerance=DEFAULT_UFF_ENERGY_TOLERANCE;
      setAFlags(aflags);
      setKFlags(kflags);
      setVFlags(vflags);
      setPOccStructure(xstr_pocc);
      aurostd::xoption loader;loader.flag("LOAD::KFLAGS",false);loader.flag("LOAD::VFLAGS",false);loader.flag("LOAD::PARTCAR",false);
      loadFromAFlags(loader);
      setPOccFlags(pocc_flags);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  void POccCalculator::setPOccFlags(const aurostd::xoption& pocc_flags) {m_p_flags=pocc_flags;}

  void POccCalculator::loadFromAFlags() { //grabs from m_aflags
    aurostd::xoption loader;
    loader.flag("LOAD::KFLAGS",true);
    loader.flag("LOAD::VFLAGS",true);
    loader.flag("LOAD::PARTCAR",true);
    return loadFromAFlags(loader);
  }
  void POccCalculator::loadFromAFlags(const aurostd::xoption& loader) { //grabs from m_aflags
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    string AflowIn_file="",AflowIn="";
    KBIN::getAflowInFromAFlags(m_aflags,AflowIn_file,AflowIn,*p_FileMESSAGE,*p_oss);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded aflow.in" << endl;}
    if(loader.flag("LOAD::KFLAGS")){
      m_kflags=KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn,*p_FileMESSAGE,m_aflags,*p_oss);  //set them here if we can, they will get overwritten with input kflags
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUP_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_PGROUP_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUPK_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_PGROUPK_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_FGROUP_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_FGROUP_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_IATOMS_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_IATOMS_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_AGROUP_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_AGROUP_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_SGROUP_WRITE(pre)=" << m_kflags.KBIN_SYMMETRY_SGROUP_WRITE << endl;
      }
      pflow::defaultKFlags4SymWrite(m_kflags,false);  //if loading from aflags, we don't want to write new files in that directory
      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUP_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_PGROUP_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUPK_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_PGROUPK_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_FGROUP_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_FGROUP_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_PGROUP_XTAL_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_IATOMS_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_IATOMS_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_AGROUP_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_AGROUP_WRITE << endl;
        cerr << __AFLOW_FUNC__ << " m_kflags.KBIN_SYMMETRY_SGROUP_WRITE(post)=" << m_kflags.KBIN_SYMMETRY_SGROUP_WRITE << endl;
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded kflags" << endl;}
    }
    if(loader.flag("LOAD::VFLAGS")){
      m_vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,*p_FileMESSAGE,m_aflags,m_kflags,*p_oss); //set them here if we can, they will get overwritten with input vflags
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded vflags" << endl;}
    }
    if(loader.flag("LOAD::PARTCAR")){
      setPOccStructure(pocc::extractPARTCAR(AflowIn));
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " loaded PARTCAR" << endl;}
    }
  }

  void POccCalculator::setPOccStructure(const xstructure& xstr_pocc) {
    POccCalculatorTemplate::setPOccStructure(xstr_pocc);
    initializePOccStructure();
  }

  void POccCalculator::writePARTCAR() const {
    stringstream message;

    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccCalculator failed to initialize");}

    message << "Writing out PARTCAR";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    stringstream partcar_ss;
    partcar_ss << xstr_pocc;
    aurostd::stringstream2file(partcar_ss,m_aflags.Directory+"/PARTCAR");
  }

  void POccCalculator::setAFlags(const _aflags& aflags) {
    POccCalculatorTemplate::setAFlags(aflags);
  }
  void POccCalculator::setKFlags(const _kflags& kflags) {m_kflags.clear();m_kflags=kflags;}
  void POccCalculator::setVFlags(const _vflags& vflags) {m_vflags.clear();m_vflags=vflags;}

  string POccCalculator::getARUNDirectoryPath(uint isupercell) const {
    if(isupercell>m_ARUN_directories.size()-1){throw aurostd::xerror(__AFLOW_FILE__,XPID+"getARUNDirectoryPath():","isupercell>m_ARUN_directories.size()-1",_INDEX_BOUNDS_);}
    if(m_convolution){return m_ARUN_directories[isupercell];}
    return m_aflags.Directory+"/"+m_ARUN_directories[isupercell];
  }
  string POccCalculator::getOutputPath() const {
    if(m_convolution){return ".";} //we can change later
    return m_aflags.Directory;
  }

  bool sortPOccSites(const POccUnit& p1,const POccUnit& p2){
    if(p1.v_occupants.size()!=p2.v_occupants.size()){
      if(p1.equivalent!=p2.equivalent){
        return p1.v_occupants.size()>p2.v_occupants.size(); //sort in ascending order
      }
      //return equivalent<p2.equivalent;
    }
    if(p1.equivalent!=p2.equivalent){return p1.equivalent<p2.equivalent;}
    return p1.site<p2.site;
  }

  bool sortPOccGroups(const POccUnit& p1,const POccUnit& p2){return p1.partial_occupation_value<p2.partial_occupation_value;}

  //bool POccCalculator::updatePOccValues() {
  //  //This function is only used to update partial occupation value i.e., comp_each_type
  //  xstr_pocc.comp_each_type.clear();
  //  for(uint i=0;i<xstr_pocc.num_each_type.size();i++) {xstr_pocc.comp_each_type.push_back(0.0);}
  //
  //  for(uint i=0;i<xstr_pocc.atoms.size();i++) {
  //    xstr_pocc.atoms[i].partial_occupation_flag = !aurostd::isequal(xstr_pocc.atoms[i].partial_occupation_value,1.0,_AFLOW_POCC_ZERO_TOL_);
  //    xstr_pocc.atoms[i].partial_occupation_value = ( xstr_pocc.atoms[i].partial_occupation_flag ? xstr_pocc.atoms[i].partial_occupation_value : 1.0 );
  //    xstr_pocc.comp_each_type[xstr_pocc.atoms[i].type] += xstr_pocc.atoms[i].partial_occupation_value;
  //  }
  //
  //  //CO add check for if all sites are fully occupied
  //  return true;
  //}

  //void POccCalculator::RemoveAtom(xstructure& xstr,vector<uint>& v_atoms_to_remove){
  //  std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  //  for(uint atom=0;atom<v_atoms_to_remove.size();atom++){xstr.RemoveAtom(v_atoms_to_remove[atom]);}
  //}

  //void POccCalculator::resetMaxSLRadius() {
  //  max_superlattice_radius=0.0;
  //}

  void POccCalculator::resetHNFMatrices() {
    hnf_mat.clear();
    //a_start=1;c_start=1;f_start=1;
    //b_start=0;d_start=0;e_start=0;
    a_start=A_START;c_start=C_START;f_start=F_START;
    b_start=B_START;d_start=D_START;e_start=E_START;
    v_unique_superlattices.clear();
  }

  bool POccCalculator::iterateHNFMatrix(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || ENUMERATE_ALL_HNF || XHOST.DEBUG);

    xmatrix<double> _hnf_mat(3,3),duplicate_mat(3,3);
    xmatrix<double> superlattice(3,3), rotated_superlattice(3,3);
    bool duplicate_H;
    //int n_hnf=p_str.n_hnf;
    const xmatrix<double>& lattice=getLattice();
    const vector<_sym_op>& pgroup=getPGroup();

    bool eliminate_by_pgroup=true;

    //this allows us to reset completely everytime (SAFE but inefficient)
    //int a_start=1;int c_start=1;int f_start=1;
    //int b_start=0;int d_start=0;int e_start=0;

    // PHYSICAL REVIEW B 77, 224115 2008
    // Algorithm for generating derivative structures
    // Gus L. W. Hart and Rodney W. Forcade
    // IMPORTANT DISTINCTION FROM PAPER, our lattice is in row matrix, not columns
    // let's work in row space to avoid so many transpose operations
    // also, slight optimization for both speed and memory:
    // we leverage a scheme that remembers where it was last, and picks up here in the nested for-loop
    // it relies on saving a to a_start, b to b_start, etc. when it finds a unique hnf_mat
    // then we start the loops at these values only for the initial loop through
    // starting_config ensures that we reset the loops back to A_START, B_START, etc. after the initial loop
    // that way, we explore all possibilities efficiently without having to save the hnf's (only corresponding superlattices)
    // this iterative scheme helps save some memory, which we need elsewhere in the code
    // starting_config ensures we don't ALSO lose calculation iterations in recalculating hnf possibilities we already explored
    // REMEMBER when comparing with KY's hnf's, his are TRANSPOSED!
    bool starting_config;
    for(int a=a_start;a<=n_hnf;a++){
      starting_config=(a==a_start);
      for(int c=( starting_config ? c_start : C_START );c<=n_hnf/a;c++) //for(int c=( a==a_start ? c_start : C_START );c<=n_hnf/a;c++) //for(int c=c_start;c<=n_hnf/a;c++)
      {   //CO20200106 - patching for auto-indenting
        starting_config=(starting_config && (c==c_start));
        for(int f=( starting_config ? f_start : F_START );f<=n_hnf/a/c;f++) //for(int f=( (a==a_start && c==c_start) ? f_start : F_START );f<=n_hnf/a/c;f++) //for(int f=f_start;f<=n_hnf/a/c;f++)
        { //CO20200106 - patching for auto-indenting
          if(a*c*f==n_hnf){
            //found new viable diagonal set, now enumerate based on off-diagonals
            starting_config=(starting_config && (f==f_start));
            for(int b=( starting_config ? b_start : B_START );b<c;b++) //for(int b=( (a==a_start && c==c_start && f==f_start) ? b_start : B_START );b<c;b++) //for(int b=b_start;b<c;b++)
            { //CO20200106 - patching for auto-indenting
              starting_config=(starting_config && (b==b_start));
              for(int d=( starting_config ? d_start : D_START );d<f;d++) //for(int d=( (a==a_start && c==c_start && f==f_start && b==b_start) ? d_start : D_START );d<f;d++) //for(int d=d_start;d<f;d++)
              { //CO20200106 - patching for auto-indenting
                starting_config=(starting_config && (d==d_start));
                for(int e=( starting_config ? e_start : E_START );e<f;e++) //for(int e=( (a==a_start && c==c_start && f==f_start && b==b_start && d==d_start) ? e_start : E_START );e<f;e++) //for(int e=e_start;e<f;e++)
                { //CO20200106 - patching for auto-indenting
                  if(LDEBUG) {
                    cerr << " a=" << a << "(max=n_hnf=" << n_hnf << ")";
                    cerr << " c=" << c << "(max=n_hnf/a=" << n_hnf/a << ")";
                    cerr << " f=" << f << "(max=n_hnf/a/c=" << n_hnf/a/c << ")";
                    cerr << " b=" << b << "(max=c=" << c-1 << ")";
                    cerr << " d=" << d << "(max=f=" << f-1 << ")";
                    cerr << " e=" << e << "(max=f=" << f-1 << ")";
                    cerr << endl;
                  }

                  _hnf_mat(1,1)=a;_hnf_mat(1,2)=b;_hnf_mat(1,3)=d;
                  _hnf_mat(2,1)=0;_hnf_mat(2,2)=c;_hnf_mat(2,3)=e;
                  _hnf_mat(3,1)=0;_hnf_mat(3,2)=0;_hnf_mat(3,3)=f;

                  if(LDEBUG) {
                    cerr << "HNF MAT " << endl;
                    cerr << _hnf_mat << endl;
                  }

                  superlattice=_hnf_mat*lattice;

                  if(LDEBUG) {
                    cerr << "SUPERLATTICE " << endl;
                    cerr << superlattice << endl;
                  }

                  duplicate_H=false;
                  for(uint i=0;i<v_unique_superlattices.size() && !duplicate_H;i++){
                    if(LDEBUG) {
                      cerr << "SUPERLATTICES OLD " << endl;
                      cerr << v_unique_superlattices[i] << endl;
                    }
                    for(uint pg=0;pg<(eliminate_by_pgroup==true ? pgroup.size() : 1) && !duplicate_H;pg++){ //identity is ALWAYS the first pgroup by construction (see TRICK)
                      if(LDEBUG) {
                        cerr << "POINT GROUP " << endl;
                        cerr << pgroup[pg].Uc << endl;
                      }
                      //original matrix in column space: Bi = R * Bj * H, Gus derives H = Bj^-1 * R^-1 * Bi
                      //in row space:  
                      //(Bi)^T = (R * Bj * H)^T
                      //Bi^T = H^T * Bj^T * R^T
                      //Bi^T * (R^T)^-1 * (Bj^T)^-1 = H^T
                      duplicate_mat = (superlattice * aurostd::inverse(pgroup[pg].Uc) * aurostd::inverse(v_unique_superlattices[i]));
                      duplicate_H = aurostd::isinteger(duplicate_mat);
                      if(LDEBUG) {
                        cerr << "DUPLICATE MATRIX" << endl;
                        cerr << duplicate_mat << endl;
                      }
                      //duplicate_H=false;
                    }
                  }
                  //cerr << "n_hnf " << n_hnf << endl;
                  //cerr << "lattice " << endl;
                  //cerr << lattice << endl;
                  //cerr << "hnf_mat" << endl;
                  //cerr << _hnf_mat << endl;
                  //cerr << "superlattice " << endl;
                  //cerr << superlattice << endl;
                  ////duplicate_H=true;
                  if(ENUMERATE_ALL_HNF){
                    if(!duplicate_H){ //still get FOUND output
                      cerr << "----------------------" << endl;
                      cerr << "FOUND " << endl;
                      cerr << _hnf_mat << endl;
                      cerr << "----------------------" << endl;
                      duplicate_H=true; //so we go through all iterations
                    }
                  }
                  if(!duplicate_H){
                    if(LDEBUG) {
                      cerr << "----------------------" << endl;
                      cerr << "FOUND " << endl;
                      cerr << _hnf_mat << endl;
                      cerr << "----------------------" << endl;
                    }
                    //found new matrix, set new starting conditions and return
                    a_start=a;b_start=b;c_start=c;
                    d_start=d;e_start=e;f_start=f;
                    hnf_mat=_hnf_mat;
                    v_unique_superlattices.push_back(superlattice);
                    return true;

                  }
                }
              }
            }
          }
        }
      }
    }

    return false;
  }

  void POccCalculator::resetSiteConfigurations() {
    v_types_config.clear();
    config_iterator=0;
    //v_config_iterators.clear();
    v_config_order.clear();
  }

  void POccCalculator::setConfigOrder(){
    vector<uint> getConfigOrder(vector<vector<int> >& v_types_config);
    v_config_order=getConfigOrder(v_types_config);
  }

  //CO20190205 - come back!
  //we do config_order based on count of vacancies
  //not sure this is optimal
  //we want to do it based on which vacancy configurations are most similar (so as to not recalculate the cluster bonding)
  //so consider creating a new vector which is 0 everywhere EXCEPT where there is a vacancy
  //and sort this vector lexicographically
  vector<uint> getConfigOrder(vector<vector<int> >& v_types_config){
    vector<SiteVacancyCount> v_svc;
    for(uint site=0;site<v_types_config.size();site++){
      v_svc.push_back(SiteVacancyCount());
      v_svc.back().site=site;
      v_svc.back().vacancy_count=0;
      for(uint pos=0;pos<v_types_config[site].size();pos++){
        if(v_types_config[site][pos]==-1){v_svc.back().vacancy_count++;}
      }
    }
    std::sort(v_svc.begin(),v_svc.end()); //most vacancies should go last
    vector<uint> v_config_order;
    for(uint i=0;i<v_svc.size();i++){v_config_order.push_back(v_svc[i].site);}
    return v_config_order;
  }

  bool POccCalculator::getNextSiteConfiguration() {
    //initialize everything
    if(v_types_config.size()==0){
      config_iterator=0;
      //const vector<StructureConfiguration>& v_str_configs=p_str.v_str_configs;
      for(uint site=0;site<v_str_configs[config_iterator].site_configs.size();site++){
        v_types_config.push_back(v_str_configs[config_iterator].site_configs[site].getStartingTypesConfiguration());
      }
      setConfigOrder();
      return true;
    }
    bool getNextSiteConfiguration(vector<uint>& v_config_order,vector<vector<int> >& v_types_config);
    if(!getNextSiteConfiguration(v_config_order,v_types_config)){
      //if(config_iterator>=p_str.v_str_configs.size()-1){return false;}
      if(config_iterator>=v_str_configs.size()-1){return false;}
      config_iterator++;
      //const vector<StructureConfiguration>& v_str_configs=p_str.v_str_configs;
      for(uint site=0;site<v_types_config.size();site++){v_types_config[site].clear();} v_types_config.clear();
      for(uint site=0;site<v_str_configs[config_iterator].site_configs.size();site++){
        //v_config_iterators.push_back(0);
        v_types_config.push_back(v_str_configs[config_iterator].site_configs[site].getStartingTypesConfiguration());
      }
      setConfigOrder();
      return true;
    }
    return true;
  }

  bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,
      vector<int>& v_config_iterators, 
      vector<vector<int> >& v_types_config) {
    if(v_types_config.size()==0){
      v_config_iterators.clear();
      //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
      //cerr << "SEE " << vv_count_configs.size() << endl;
      for(uint site=0;site<vv_count_configs.size();site++){
        v_config_iterators.push_back(0);
        v_types_config.push_back(vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration());
      }
      //v_config_order=getConfigOrder(v_types_config);
      //cerr << "WOW2 " << v_types_config.size() << endl;
      return true;
    }
    for(uint site=0;site<v_config_iterators.size();site++){
      if(v_config_iterators[site]>=(int)vv_count_configs[site].size()-1){  //greater than not necessary but safe
        if(site==v_config_iterators.size()-1){break;} //stop condition
        v_config_iterators[site]=0; 
        v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();
        continue;
      }
      v_config_iterators[site]++;
      v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();
      return true;
    }
    return false;
  }

  //bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,
  //    vector<uint>& v_config_order,
  //    vector<int>& v_config_iterators, 
  //    vector<vector<int> >& v_types_config) {
  //  //starting condition START - v_types_config is empty!
  //  if(v_types_config.size()==0){
  //    v_config_iterators.clear();
  //    //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
  //    //cerr << "SEE " << vv_count_configs.size() << endl;
  //    for(uint site=0;site<vv_count_configs.size();site++){
  //      v_config_iterators.push_back(0);
  //      v_types_config.push_back(vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration());
  //    }
  //    v_config_order=getConfigOrder(v_types_config);
  //    //cerr << "WOW2 " << v_types_config.size() << endl;
  //    return true;
  //  }
  //  //starting condition STOP - v_types_config is empty!
  //  //otherwise, similar next bitstring generator
  //  if(!getNextSiteConfiguration(v_config_order,v_types_config)){
  //    //const vector<vector<POccSiteConfiguration> >& vv_count_configs=p_str.v_str_configs;
  //    //if(1){
  //    //for(uint site=0;site<vv_count_configs.size();site++){
  //    //  if(v_config_iterators[site]>=(int)vv_count_configs[site].size()-1){
  //    //    if(site==vv_count_configs.size()-1){break;} //stop condition
  //    //    v_config_iterators[site]=0;
  //    //    v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
  //    //    continue;
  //    //  }
  //    //  v_config_iterators[site]++;
  //    //  v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
  //    //  v_config_order=getConfigOrder(v_types_config);
  //    //  return true;
  //    //}
  //    //}
  //    //if(0){
  //    uint site;
  //    for(uint index=0;index<v_config_order.size();index++){
  //      site=v_config_order[index];
  //      //cerr << "index = " << index << " site=" << site << " v_config_iterators[site]=" << v_config_iterators[site] << endl;
  //      if(v_config_iterators[site]>=(int)vv_count_configs[site].size()-1){
  //        //if(site==vv_count_configs.size()-1){break;} //stop condition
  //        if(index==vv_count_configs.size()-1){break;} //stop condition
  //        v_config_iterators[site]=0;
  //        v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
  //        continue;
  //      }
  //      //cerr << "HERE" << endl;
  //      v_config_iterators[site]++;
  //      v_types_config[site]=vv_count_configs[site][v_config_iterators[site]].getStartingTypesConfiguration();  //reset configuration
  //      v_config_order=getConfigOrder(v_types_config);
  //      return true;
  //    }
  //    //}
  //    return false; //stop condition SAFE
  //  }
  //  return true;
  //}

  bool POccCalculator::getNextSiteConfiguration(vector<vector<int> >& v_types_config) {
    bool getNextSiteConfiguration(vector<uint>& v_config_order,vector<vector<int> >& v_types_config);
    return getNextSiteConfiguration(v_config_order,v_types_config);
  }

  bool getNextSiteConfiguration(vector<uint>& v_config_order,
      vector<vector<int> >& v_types_config) {
    //bool getNextSiteConfiguration(vector<int>& types_config);
    //if(1){
    //for(uint site=0;site<v_types_config.size();site++){
    //  if(!getNextSiteConfiguration(v_types_config[site])){
    //    //reset site first
    //    //we know that the last permutation is the reverse sorted of original vector
    //    //we can either a) completely sort or b) simply reverse (slightly faster)
    //    //std::sort(v_types_config[site].begin(),v_types_config[site].end());
    //    std::reverse(v_types_config[site].begin(),v_types_config[site].end());
    //    if(site==v_types_config.size()-1){return false;} //stop condition
    //    continue;
    //  }
    //  return true;
    //}
    //return false; //stop condition SAFE
    //}
    //here we optimized ordering to avoid fewest number of bonding recalculations
    //if(0){
    uint site;
    for(uint index=0;index<v_config_order.size();index++){
      site=v_config_order[index];
      if(!getNextSiteConfiguration(v_types_config[site])){
        //reset site first
        //we know that the last permutation is the reverse sorted of original vector
        //we can either a) completely sort or b) simply reverse (slightly faster)
        //std::sort(v_types_config[site].begin(),v_types_config[site].end());
        std::reverse(v_types_config[site].begin(),v_types_config[site].end());
        //if(site==v_types_config.size()-1){return false;} //stop condition
        if(index==v_types_config.size()-1){return false;} //stop condition
        continue;
      }
      return true;
    }
    //}
    return false; //stop condition SAFE
  }

  bool getNextSiteConfiguration(vector<int>& types_config) {
    //Shen, MK. BIT (1962) 2: 228. doi:10.1007/BF01940170
    int _i=-1;
    int _j=-1;
    for(int i=1;i<(int)types_config.size();i++){if(types_config[i-1]<types_config[i]&&(i>_i)){_i=i;}}
    if(_i==-1){return false;} //stop condition
    for(int j=0;j<(int)types_config.size();j++){if(types_config[_i-1]<types_config[j]&&(j>_j)){_j=j;}}
    //cerr << "i=" << _i << "  j=" << _j << " ";
    std::swap(types_config[_i-1],types_config[_j]);
    for(int i=0;i<((int)types_config.size()-_i+1)/2;i++){std::swap(types_config[_i+i],types_config[types_config.size()-i-1]);}
    return true;
  }

  ////simultaneously define and calculate
  //double POccCalculator::getUFFEnergy() {energy_analyzer.setBonds(v_types_config);return energy_analyzer.getUFFEnergy();}

  //void POccCalculator::calculateHNF(){p_str.calculateHNF();}               //get n_hnf

  void POccCalculator::getTotalPermutationsCount(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    message << "Getting total number of supercell decoration permutations";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    //unsigned long long int hnf_count=0;
    //unsigned long long int types_config_permutations_count=0;   //per hnf matrix

    //starting criteria for HNF matrices
    //int a_start=1;int c_start=1;int f_start=1;
    //int b_start=0;int d_start=0;int e_start=0;
    xmatrix<double> hnf_mat;                            //really xmatrix of int's, but we rule for int * double in xmatrix, no big deal
    //vector<xmatrix<double> > v_unique_superlattices;    //only store locally, as we need to make comparisons

    //START KY
    //cerr << "START KY " << endl;
    //CalculateHNF(p_str.xstr_pocc,p_str.n_hnf);

    //cerr << endl;
    //cerr << "START CO " << endl;
    //resetMaxSLRadius();
    hnf_count=0;
    resetHNFMatrices();
    while(iterateHNFMatrix()){hnf_count++;}
    //cerr << "max_radius " << max_superlattice_radius << endl;
    //for(uint i=0;i<v_unique_superlattices.size();i++){
    //  //cerr << LatticeDimensionSphere(v_unique_superlattices[i],max_superlattice_radius/24.0) << endl;
    //  xvector<int> dims = LatticeDimensionSphere(v_unique_superlattices[i],10.0);
    //  cerr << "dims = " << dims << endl;
    //  xmatrix<double> scell(3,3);
    //  scell(1,1)=dims(1);scell(2,2)=dims(2);scell(3,3)=dims(3);
    //  xmatrix<double> newlattice = scell * v_unique_superlattices[i];
    //  cerr << "radius = " << RadiusSphereLattice(newlattice) << endl;
    //}
    message << "Total count of unique HNF matrices = " << hnf_count;
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    //starting criteria for site combinations
    //vector<vector<int> > v_types_config;
    //vector<int> v_config_iterators; //, v_site_iterators;

    //while(getNextSiteConfiguration(v_types_config,v_config_iterators,v_site_iterators)){  //[CO20200106 - close bracket for indenting]}

    //CO MAKE THIS A TRUE CALCULATOR, simply go through each site, add between configs, 
    //multiple across sites
    types_config_permutations_count=0;   //per hnf matrix
    unsigned long long int str_config_permutations_count;
    unsigned long long int config_permutations_count=0;         //per hnf matrix
    total_permutations_count=0;  //per hnf matrix
    vector<int> current_config;
    //const vector<StructureConfiguration>& v_str_configs=p_str.v_str_configs;
    bool getNextSiteConfiguration(vector<int>& types_config);
    for(uint str_config=0;str_config<v_str_configs.size();str_config++){
      str_config_permutations_count=1;
      for(uint site=0;site<v_str_configs[str_config].site_configs.size();site++){
        //for(uint config=0;config<v_str_configs[site].size();config++){
        config_permutations_count=0;
        current_config=v_str_configs[str_config].site_configs[site].getStartingTypesConfiguration();
        //if(config>0){config_permutations_count++;}
        config_permutations_count++;	//for starting config
        while(getNextSiteConfiguration(current_config)){config_permutations_count++;}
        str_config_permutations_count*=config_permutations_count;
        //}
        //str_config_permutations_count*=config_permutations_count;
      }
      types_config_permutations_count+=str_config_permutations_count;
    }
    //cerr << types_config_permutations_count << endl;
    message << "Total count of unique types-configuration permutations = " << types_config_permutations_count;
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Checking unique types-configuration permutation count" << endl;
      types_config_permutations_count=0;
      resetSiteConfigurations();
      while(POccCalculator::getNextSiteConfiguration()){
        types_config_permutations_count++; 
        //cerr << "WOW1 " << v_types_config.size() << endl;
        //for(uint i=0;i<v_types_config.size();i++){
        //  cerr << "LOOK i=" << i << " " << v_types_config[i].size() << endl;
        //}
        cerr << __AFLOW_FUNC__ << " Permutation " << types_config_permutations_count << "  "; 
        for(uint i=0;i<v_types_config.size();i++){
          for(uint j=0;j<v_types_config[i].size();j++){
            cerr << v_types_config[i][j] << " ";
          } 
          cerr << "    ";
        }
        cerr << endl;
      }
      cerr << __AFLOW_FUNC__ << "Total count of unique types-configuration permutations = " << types_config_permutations_count << endl;
    }

    total_permutations_count=hnf_count*types_config_permutations_count;

    message << "Total number of supercell decoration permutations = " <<total_permutations_count;
    char LOGGER_TYPE=_LOGGER_MESSAGE_;
    if(m_p_flags.flag("POCC_COUNT_TOTAL")){LOGGER_TYPE=_LOGGER_COMPLETE_;}
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,LOGGER_TYPE);
  }


  //void POccCalculator::getUFFParamBond(uint type1,uint type2,UFFParamBond& uffb){
  //}

  //double POccCalculator::getUFFBondEnergy(xstructure& xstr,vector<uint>& v_vacancies,
  //xmatrix<double>& distance_matrix,vector<vector<uint> >& v_bonded_atom_indices,
  //uint MODE){

  //double energy=0.0;
  //bool found_vacancy;
  //uint atom1,atom2,type1,type2;
  //UFFParamBond uffb;

  //for(uint i=0;i<v_bonded_atom_indices.size();i++){
  //found_vacancy=false;
  //for(uint vi=0;vi<v_vacancies.size() && !found_vacancy;vi++){
  //found_vacancy=(v_bonded_atom_indices[i][0]==v_vacancies[vi] ||
  //v_bonded_atom_indices[i][1]==v_vacancies[vi]);
  //}
  //if(found_vacancy){cerr<<"VACANCYYYY"<<endl;continue;}
  //atom1=v_bonded_atom_indices[i][0]; atom2=v_bonded_atom_indices[i][1];
  //type1=xstr.atoms[atom1].type; type2=xstr.atoms[atom2].type;
  ////uffb.uffp1=&types2uffparams_map[type1]; uffb.uffp2=&types2uffparams_map[type2];
  ////uffb.distij=distance_matrix(atom1,atom2);
  ////simultaneously initialize and calculate uffb
  //uffb.calculate(types2uffparams_map[type1],types2uffparams_map[type2],distance_matrix(atom1,atom2));
  //if(MODE==BOND_MODE){energy += 0.5 * uffb.Kij * uffb.delta * uffb.delta;}
  //else if(MODE==NONBOND_MODE){energy += uffb.Dij * (uffb.X12 - 2.0 * uffb.X6);}
  //}
  //energy*=kcal2eV;
  //return energy;
  //}

  //xstructure POccUFFEnergyAnalyzer::getXStructure(vector<vector<int> >& _v_types_config){
  //v_types_config=&_v_types_config;

  //const xstructure& xstr_nopocc = p_str.xstr_nopocc;
  //const vector<Site>& pocc_sites = p_str.m_pocc_sites;
  //vector<int> sc2pc_map, pc2sc_map;
  //xstructure supercell=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false);

  //uint starting_supercell_atom_index,pocc_atom_index,supercell_atom_index;
  //vector<uint> v_atoms_to_remove;
  //for(uint site=0;site<v_types_config.size();site++){
  //if(pocc_sites[site].partial_occupation_flag){  //slight optimization
  ////test of stupidity
  //if(pocc_sites[site].v_occupants.size()<2){throw aurostd::xerror(__AFLOW_FILE__,"pocc::POccUFFEnergyAnalyzer::getXStructure():","pocc_sites[site].v_occupants.size()<2",_RUNTIME_ERROR_);}
  ////find index of atom first in v_occupants list, that's the one that remained
  //starting_supercell_atom_index=pc2sc_map[pocc_sites[site].v_occupants[0]];
  //for(uint i=0;i<v_types_config[site].size();i++){
  ////switch these atoms
  ////supercell.atoms[starting_supercell_atom_index+i]
  ////xstr_pocc.atoms[types2pc_map[v_types_config[site][i]]]
  ////pocc_atom_index=types2pc_map[v_types_config[site][i]];
  //pocc_atom_index=p_str.types2pc_map[v_types_config[site][i]];
  //supercell_atom_index=starting_supercell_atom_index+i;
  //const _atom& pocc_atom = p_str.xstr_pocc.atoms[pocc_atom_index];
  //if(v_types_config[site][i]>=0){copyAtomAttributes(pocc_atom,supercell.atoms[supercell_atom_index]);}
  //else {v_atoms_to_remove.push_back(supercell_atom_index);}
  //}
  //}
  //}
  ////std::sort(v_atoms_to_remove.rbegin(),v_atoms_to_remove.rend()); //NOTE the r, reverse sort, that way when we remove, it doesn't affect other indices
  ////for(uint atom=0;atom<v_atoms_to_remove.size();atom++){supercell.RemoveAtom(v_atoms_to_remove[atom]);}

  ////calculate UFF energy here before you rearrange the structure
  //double energy=0.0;
  ////energy += getUFFBondEnergy(supercell,v_atoms_to_remove,distance_matrix,v_bonded_atom_indices,BOND_MODE);
  ////energy += getUFFBondEnergy(supercell,v_atoms_to_remove,distance_matrix,v_nonbonded_atom_indices,NONBOND_MODE);

  //cerr << energy << endl;



  ////rearrange the structure below here

  ////RemoveAtom(supercell,v_atoms_to_remove);
  //supercell.RemoveAtom(v_atoms_to_remove);
  ////supercell.SpeciesPutAlphabetic();

  ////cerr << pocc::CalculateUFFEnergy(supercell) << endl;
  ////cerr << supercell << endl;
  //return supercell;
  //}

  string getUFFParameterString(const string& atom) {
    if (atom=="H")  {return "H   0.354   180.000 2.886   0.044   12.000  0.712   0.000   0.000   4.528   6.945   0.371   ";}   
    if (atom=="He") {return "He  0.849   90.000  2.362   0.056   15.240  0.098   0.000   0.000   9.660   14.920  1.300   ";}   
    if (atom=="Li") {return "Li  1.336   180.000 2.451   0.025   12.000  1.026   0.000   2.000   3.006   2.386   1.557   ";}   
    if (atom=="Be") {return "Be  1.074   109.470 2.745   0.085   12.000  1.565   0.000   2.000   4.877   4.443   1.240   ";}  
    if (atom=="B")  {return "B   0.838   109.470 4.083   0.180   12.052  1.755   0.000   2.000   5.110   4.750   0.822   ";}  
    if (atom=="C")  {return "C   0.706   180.000 3.851   0.105   12.730  1.912   0.000   2.000   5.343   5.063   0.759   ";}  
    if (atom=="N")  {return "N   0.656   180.000 3.660   0.069   13.407  2.544   0.000   2.000   6.899   5.880   0.715   ";}  
    if (atom=="O")  {return "O   0.639   180.000 3.500   0.060   14.085  2.300   0.000   2.000   8.741   6.682   0.669   ";}  
    if (atom=="F")  {return "F   0.668   180.000 3.364   0.050   14.762  1.735   0.000   2.000   10.874  7.474   0.706   ";}  
    if (atom=="Ne") {return "Ne  0.920   90.000  3.243   0.042   15.440  0.194   0.000   2.000   11.040  10.550  1.768   ";}  
    if (atom=="Na") {return "Na  1.539   180.000 2.983   0.030   12.000  1.081   0.000   1.250   2.843   2.296   2.085   ";}  
    if (atom=="Mg") {return "Mg  1.421   109.470 3.021   0.111   12.000  1.787   0.000   1.250   3.951   3.693   1.500   ";}  
    if (atom=="Al") {return "Al  1.244   109.470 4.499   0.505   11.278  1.792   0.000   1.250   4.060   3.590   1.201   ";}  
    if (atom=="Si") {return "Si  1.117   109.470 4.295   0.402   12.175  2.323   1.225   1.250   4.168   3.487   1.176   ";}  
    if (atom=="P")  {return "P   1.101   93.800  4.147   0.305   13.072  2.863   2.400   1.250   5.463   4.000   1.102   ";}  
    if (atom=="S")  {return "S   1.064   92.1000 4.035   0.274   13.969  2.703   0.000   1.250   6.928   4.486   1.047   ";}  
    if (atom=="Cl") {return "Cl  1.044   180.000 3.947   0.227   14.866  2.348   0.000   1.250   8.564   4.946   0.994   ";}  
    if (atom=="Ar") {return "Ar  1.032   90.000  3.868   0.185   15.763  0.300   0.000   1.250   9.465   6.355   2.108   ";}  
    if (atom=="K")  {return "K   1.953   180.000 3.812   0.035   12.000  1.165   0.000   0.700   2.421   1.920   2.586   ";}  
    if (atom=="Ca") {return "Ca  1.761   90.000  3.399   0.238   12.000  2.141   0.000   0.700   3.231   2.880   2.000   ";}  
    if (atom=="Sc") {return "Sc  1.513   109.470 3.295   0.019   12.000  2.592   0.000   0.700   3.395   3.080   1.750   ";}  
    if (atom=="Ti") {return "Ti  1.412   109.470 3.175   0.017   12.000  2.659   0.000   0.700   3.470   3.380   1.607   ";}  
    if (atom=="V")  {return "V   1.402   109.470 3.144   0.016   12.000  2.679   0.000   0.700   3.650   3.410   1.470   ";}  
    if (atom=="Cr") {return "Cr  1.345   90.000  3.023   0.015   12.000  2.463   0.000   0.700   3.415   3.865   1.402   ";}  
    if (atom=="Mn") {return "Mn  1.382   90.000  2.961   0.013   12.000  2.430   0.000   0.700   3.325   4.105   1.533   ";}  
    if (atom=="Fe") {return "Fe  1.270   109.470 2.912   0.013   12.000  2.430   0.000   0.700   3.760   4.140   1.393   ";}  
    if (atom=="Co") {return "Co  1.241   90.000  2.872   0.014   12.000  2.430   0.000   0.700   4.105   4.175   1.406   ";}  
    if (atom=="Ni") {return "Ni  1.164   90.000  2.834   0.015   12.000  2.430   0.000   0.700   4.465   4.205   1.398   ";}  
    if (atom=="Cu") {return "Cu  1.302   109.470 3.495   0.005   12.000  1.756   0.000   0.700   4.200   4.220   1.434   ";}  
    if (atom=="Zn") {return "Zn  1.193   109.470 2.763   0.124   12.000  1.308   0.000   0.700   5.106   4.285   1.400   ";}  
    if (atom=="Ga") {return "Ga  1.260   109.470 4.383   0.415   11.000  1.821   0.000   0.700   3.641   3.160   1.211   ";}  
    if (atom=="Ge") {return "Ge  1.197   109.470 4.280   0.379   12.000  2.789   0.701   0.700   4.051   3.438   1.189   ";}  
    if (atom=="As") {return "As  1.211   92.100  4.230   0.309   13.000  2.864   1.500   0.700   5.188   3.809   1.204   ";}  
    if (atom=="Se") {return "Se  1.190   90.600  4.205   0.291   14.000  2.764   0.335   0.700   6.428   4.131   1.224   ";}  
    if (atom=="Br") {return "Br  1.192   180.000 4.189   0.251   15.000  2.519   0.000   0.700   7.790   4.425   1.141   ";}  
    if (atom=="Kr") {return "Kr  1.147   90.000  4.141   0.220   16.000  0.452   0.000   0.700   8.505   5.715   2.270   ";}  
    if (atom=="Rb") {return "Rb  2.260   180.000 4.114   0.040   12.000  1.592   0.000   0.200   2.331   1.846   2.770   ";}  
    if (atom=="Sr") {return "Sr  2.052   90.000  3.641   0.235   12.000  2.449   0.000   0.200   3.024   2.440   2.415   ";}  
    if (atom=="Y")  {return "Y   1.698   109.470 3.345   0.072   12.000  3.257   0.000   0.200   3.830   2.810   1.998   ";}  
    if (atom=="Zr") {return "Zr  1.564   109.470 3.124   0.069   12.000  3.667   0.000   0.200   3.400   3.550   1.758   ";}  
    if (atom=="Nb") {return "Nb  1.473   109.470 3.165   0.059   12.000  3.618   0.000   0.200   3.550   3.380   1.603   ";}  
    if (atom=="Mo") {return "Mo  1.467   90.000  3.052   0.056   12.000  3.400   0.000   0.200   3.465   3.755   1.530   ";}  
    if (atom=="Tc") {return "Tc  1.322   90.000  2.998   0.048   12.000  3.400   0.000   0.200   3.290   3.990   1.500   ";}  
    if (atom=="Ru") {return "Ru  1.478   90.000  2.963   0.056   12.000  3.400   0.000   0.200   3.575   4.015   1.500   ";}  
    if (atom=="Rh") {return "Rh  1.332   90.000  2.929   0.053   12.000  3.500   0.000   0.200   3.975   4.005   1.509   ";}  
    if (atom=="Pd") {return "Pd  1.338   90.000  2.899   0.048   12.000  3.210   0.000   0.200   4.320   4.000   1.544   ";}  
    if (atom=="Ag") {return "Ag  1.386   180.000 3.148   0.036   12.000  1.956   0.000   0.200   4.436   3.134   1.622   ";}  
    if (atom=="Cd") {return "Cd  1.403   109.470 2.848   0.228   12.000  1.650   0.000   0.200   5.034   3.957   1.600   ";}  
    if (atom=="In") {return "In  1.459   109.470 4.463   0.599   11.000  2.070   0.000   0.200   3.506   2.896   1.404   ";}  
    if (atom=="Sn") {return "Sn  1.398   109.470 4.392   0.567   12.000  2.961   0.199   0.200   3.987   3.124   1.354   ";}  
    if (atom=="Sb") {return "Sb  1.407   91.600  4.420   0.449   13.000  2.704   1.100   0.200   4.899   3.342   1.404   ";}  
    if (atom=="Te") {return "Te  1.386   90.250  4.470   0.398   14.000  2.882   0.300   0.200   5.816   3.526   1.380   ";}  
    if (atom=="I")  {return "I   1.382   180.000 4.500   0.339   15.000  2.650   0.000   0.200   6.822   3.762   1.333   ";}  
    if (atom=="Xe") {return "Xe  1.267   90.000  4.404   0.332   12.000  0.556   0.000   0.200   7.595   4.975   2.459   ";}  
    if (atom=="Cs") {return "Cs  2.570   180.000 4.517   0.045   12.000  1.573   0.000   0.100   2.183   1.711   2.984   ";}  
    if (atom=="Ba") {return "Ba  2.277   90.000  3.703   0.364   12.000  2.727   0.000   0.100   2.814   2.396   2.442   ";}  
    if (atom=="La") {return "La  1.943   109.470 3.522   0.017   12.000  3.300   0.000   0.100   2.836   2.742   2.071   ";}  
    if (atom=="Ce") {return "Ce  1.841   90.000  3.556   0.013   12.000  3.300   0.000   0.100   2.774   2.692   1.925   ";}  
    if (atom=="Pr") {return "Pr  1.823   90.000  3.606   0.010   12.000  3.300   0.000   0.100   2.858   2.564   2.007   ";}  
    if (atom=="Nd") {return "Nd  1.816   90.000  3.575   0.010   12.000  3.300   0.000   0.100   2.869   2.621   2.007   ";}  
    if (atom=="Pm") {return "Pm  1.801   90.000  3.547   0.009   12.000  3.300   0.000   0.100   2.881   2.673   2.000   ";}  
    if (atom=="Sm") {return "Sm  1.780   90.000  3.520   0.008   12.000  3.300   0.000   0.100   2.912   2.720   1.978   ";}  
    if (atom=="Eu") {return "Eu  1.771   90.000  3.493   0.008   12.000  3.300   0.000   0.100   2.879   2.788   2.227   ";}  
    if (atom=="Gd") {return "Gd  1.735   90.000  3.368   0.009   12.000  3.300   0.000   0.100   3.167   2.975   1.968   ";}  
    if (atom=="Tb") {return "Tb  1.732   90.000  3.451   0.007   12.000  3.300   0.000   0.100   3.018   2.834   1.954   ";}  
    if (atom=="Dy") {return "Dy  1.710   90.000  3.428   0.007   12.000  3.300   0.000   0.100   3.056   2.872   1.934   ";}  
    if (atom=="Ho") {return "Ho  1.696   90.000  3.409   0.007   12.000  3.416   0.000   0.100   3.127   2.891   1.925   ";}  
    if (atom=="Er") {return "Er  1.673   90.000  3.391   0.007   12.000  3.300   0.000   0.100   3.187   2.915   1.915   ";}  
    if (atom=="Tm") {return "Tm  1.660   90.000  3.374   0.006   12.000  3.300   0.000   0.100   3.251   2.933   2.000   ";}  
    if (atom=="Yb") {return "Yb  1.637   90.000  3.355   0.228   12.000  2.618   0.000   0.100   3.289   2.965   2.158   ";}  
    if (atom=="Lu") {return "Lu  1.671   90.000  3.640   0.041   12.000  3.271   0.000   0.100   2.963   2.463   1.896   ";}  
    if (atom=="Hf") {return "Hf  1.611   109.470 3.141   0.072   12.000  3.921   0.000   0.100   3.700   3.400   1.759   ";}  
    if (atom=="Ta") {return "Ta  1.511   109.470 3.170   0.081   12.000  4.075   0.000   0.100   5.100   2.850   1.605   ";}  
    if (atom=="W")  {return "W   1.392   90.000  3.069   0.067   12.000  3.700   0.000   0.100   4.630   3.310   1.538   ";}  
    if (atom=="Re") {return "Re  1.372   90.000  2.954   0.066   12.000  3.700   0.000   0.100   3.960   3.920   1.600   ";}  
    if (atom=="Os") {return "Os  1.372   90.000  3.120   0.037   12.000  3.700   0.000   0.100   5.140   3.630   1.700   ";}  
    if (atom=="Ir") {return "Ir  1.371   90.000  2.840   0.073   12.000  3.731   0.000   0.100   5.000   4.000   1.866   ";}  
    if (atom=="Pt") {return "Pt  1.364   90.000  2.754   0.080   12.000  3.382   0.000   0.100   4.790   4.430   1.557   ";}  
    if (atom=="Au") {return "Au  1.262   90.000  3.293   0.039   12.000  2.625   0.000   0.100   4.894   2.586   1.618   ";}  
    if (atom=="Hg") {return "Hg  1.340   180.000 2.705   0.385   12.000  1.750   0.000   0.100   6.270   4.160   1.600   ";}  
    if (atom=="Tl") {return "Tl  1.518   120.000 4.347   0.680   11.000  2.068   0.000   0.100   3.200   2.900   1.530   ";}  
    if (atom=="Pb") {return "Pb  1.459   109.470 4.297   0.663   12.000  2.846   0.100   0.100   3.900   3.530   1.444   ";}  
    if (atom=="Bi") {return "Bi  1.512   90.000  4.370   0.518   13.000  2.470   1.000   0.100   4.690   3.740   1.514   ";}  
    if (atom=="Po") {return "Po  1.500   90.000  4.709   0.325   14.000  2.330   0.300   0.100   4.210   4.210   1.480   ";}  
    if (atom=="At") {return "At  1.545   180.000 4.750   0.284   15.000  2.240   0.000   0.100   4.750   4.750   1.470   ";}  
    if (atom=="Rn") {return "Rn  1.420   90.000  4.765   0.248   16.000  0.583   0.000   0.100   5.370   5.370   2.200   ";}  
    if (atom=="Fr") {return "Fr  2.880   180.000 4.900   0.050   12.000  1.847   0.000   0.000   2.000   2.000   2.300   ";}  
    if (atom=="Ra") {return "Ra  2.512   90.000  3.677   0.404   12.000  2.920   0.000   0.000   2.843   2.434   2.200   ";}  
    if (atom=="Ac") {return "Ac  1.983   90.000  3.478   0.033   12.000  3.900   0.000   0.000   2.835   2.835   2.108   ";}  
    if (atom=="Th") {return "Th  1.721   90.000  3.396   0.026   12.000  4.202   0.000   0.000   3.175   2.905   2.018   ";}  
    if (atom=="Pa") {return "Pa  1.711   90.000  3.424   0.022   12.000  3.900   0.000   0.000   2.985   2.905   1.800   ";}  
    if (atom=="U")  {return "U   1.684   90.000  3.395   0.022   12.000  3.900   0.000   0.000   3.341   2.853   1.713   ";}  
    if (atom=="Np") {return "Np  1.666   90.000  3.424   0.019   12.000  3.900   0.000   0.000   3.549   2.717   1.800   ";}  
    if (atom=="Pu") {return "Pu  1.657   90.000  3.424   0.016   12.000  3.900   0.000   0.000   3.243   2.819   1.840   ";}  
    if (atom=="Am") {return "Am  1.660   90.000  3.381   0.014   12.000  3.900   0.000   0.000   2.990   3.004   1.942   ";}  
    if (atom=="Cm") {return "Cm  1.801   90.000  3.326   0.013   12.000  3.900   0.000   0.000   2.832   3.190   1.900   ";}  
    if (atom=="Bk") {return "Bk  1.761   90.000  3.339   0.013   12.000  3.900   0.000   0.000   3.194   3.036   1.900   ";}  
    if (atom=="Cf") {return "Cf  1.750   90.000  3.313   0.013   12.000  3.900   0.000   0.000   3.197   3.101   1.900   ";}  
    if (atom=="Es") {return "Es  1.724   90.000  3.299   0.012   12.000  3.900   0.000   0.000   3.333   3.089   1.900   ";}  
    if (atom=="Fm") {return "Fm  1.712   90.000  3.286   0.012   12.000  3.900   0.000   0.000   3.400   3.100   1.900   ";}  
    if (atom=="Md") {return "Md  1.689   90.000  3.274   0.011   12.000  3.900   0.000   0.000   3.470   3.110   1.900   ";}  
    if (atom=="No") {return "No  1.679   90.000  3.248   0.011   12.000  3.900   0.000   0.000   3.475   3.175   1.900   ";}  
    if (atom=="Lw") {return "Lw  1.698   90.000  3.236   0.011   12.000  3.900   0.000   0.000   3.500   3.200   1.900   ";}  
    else {return "";}                                                                                                      
  }

  vector<UFFParamAtom> getTypes2UFFParamsMap(const vector<string>& elements){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    vector<UFFParamAtom> types2uffparams_map;

    string uff_parameters_string="";
    UFFParamAtom uffp;
    vector<string> tokens;
    for(uint i=0;i<elements.size();i++){
      const string& element=elements[i];

      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " element[i=" << i << "]=" << element << endl;}
      uff_parameters_string=getUFFParameterString(element);
      if(uff_parameters_string.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to fetch UFF parameters (requested "+element+")");}
      aurostd::string2tokens(uff_parameters_string,tokens," ");
      if(tokens.size()!=12){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected UFF Parameters size",_VALUE_ILLEGAL_);}

      types2uffparams_map.push_back(uffp);
      types2uffparams_map.back().symbol=tokens[0];
      types2uffparams_map.back().r1=aurostd::string2utype<double>(tokens[1]);
      types2uffparams_map.back().theta0=aurostd::string2utype<double>(tokens[2]);
      types2uffparams_map.back().x1=aurostd::string2utype<double>(tokens[3]);
      types2uffparams_map.back().D1=aurostd::string2utype<double>(tokens[4]);
      types2uffparams_map.back().zeta=aurostd::string2utype<double>(tokens[5]);
      types2uffparams_map.back().Z1=aurostd::string2utype<double>(tokens[6]);
      types2uffparams_map.back().Vi=aurostd::string2utype<double>(tokens[7]);
      types2uffparams_map.back().Uj=aurostd::string2utype<double>(tokens[8]);
      types2uffparams_map.back().ChiI=aurostd::string2utype<double>(tokens[9]);
      types2uffparams_map.back().hard=aurostd::string2utype<double>(tokens[10]);
      types2uffparams_map.back().radius=aurostd::string2utype<double>(tokens[11]);
    }

    return types2uffparams_map;
  }

  //void POccCalculator::getBonding(xmatrix<double>& hnf_mat,xmatrix<double>& _distance_matrix,vector<vector<uint> >& v_bonded_atom_indices,
  //    vector<vector<uint> >& v_nonbonded_atom_indices) {
  //vector<int> sc2pc_map, pc2sc_map;
  //xstructure supercell=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false);

  //[OBSOLETE]xvector<double> min_vec; xvector<int> ijk;  //dummy
  //[OBSOLETE]xmatrix<double> distance_matrix(supercell.atoms.size()-1,supercell.atoms.size()-1,0,0);
  //[OBSOLETE]vector<double> v_nn_dists;
  //[OBSOLETE]for(uint atom1=0;atom1<supercell.atoms.size();atom1++){v_nn_dists.push_back(std::numeric_limits<double>::max());} //initialize with big double
  //[OBSOLETE]for(uint atom1=0;atom1<supercell.atoms.size();atom1++){
  //[OBSOLETE]  for(uint atom2=atom1+1;atom2<supercell.atoms.size();atom2++){
  //[OBSOLETE]    distance_matrix(atom1,atom2)=distance_matrix(atom2,atom1)=SYM::minimumCartesianDistance(supercell.atoms[atom1].cpos,supercell.atoms[atom2].cpos,supercell.lattice,min_vec,ijk);
  //[OBSOLETE]    if(distance_matrix(atom1,atom2)<v_nn_dists[atom1]){v_nn_dists[atom1]=distance_matrix(atom1,atom2);}
  //[OBSOLETE]  }
  //[OBSOLETE]}
  //[OBSOLETE]_distance_matrix=distance_matrix;

  //[OBSOLETE]for(uint atom1=0;atom1<supercell.atoms.size();atom1++){
  //[OBSOLETE]  for(uint atom2=atom1+1;atom2<supercell.atoms.size();atom2++){
  //[OBSOLETE]    if(abs(distance_matrix(atom1,atom2)-v_nn_dists[atom1])<0.5){  //KY standard for bonding, keep for now
  //[OBSOLETE]      v_bonded_atom_indices.push_back(vector<uint>(0));
  //[OBSOLETE]      v_bonded_atom_indices.back().push_back(atom1);
  //[OBSOLETE]      v_bonded_atom_indices.back().push_back(atom2);
  //[OBSOLETE]      //cerr << "BOND   " << atom1 << " " << atom2 << endl;
  //[OBSOLETE]    } else {
  //[OBSOLETE]      v_nonbonded_atom_indices.push_back(vector<uint>(0));
  //[OBSOLETE]      v_nonbonded_atom_indices.back().push_back(atom1);
  //[OBSOLETE]      v_nonbonded_atom_indices.back().push_back(atom2);
  //[OBSOLETE]      //cerr << "NOBOND " << atom1 << " " << atom2 << endl;
  //[OBSOLETE]    }
  //[OBSOLETE]  }
  //[OBSOLETE]}

  //[OBSOLETE]//cerr << distance_matrix << endl;


  //}

  bool POccCalculator::areEquivalentStructuresByUFF(std::list<POccSuperCellSet>::iterator it, const POccSuperCell& psc) const {
    return ((*it).getHNFIndex()==psc.m_hnf_index) && aurostd::isequal((*it).getUFFEnergy(),psc.m_energy_uff,m_energy_uff_tolerance);
  }

  //NEW
  void POccCalculator::add2DerivativeStructuresList(const POccSuperCell& psc,
      std::list<POccSuperCellSet>::iterator i_start,
      std::list<POccSuperCellSet>::iterator i_end){
    //  cerr << "FULL LIST START" << endl;
    //  for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
    //    cerr << std::fixed << std::setprecision(9) << (*it).getUFFEnergy() << " ";
    //  }
    //  cerr << endl;
    //  cerr << "FULL LIST END" << endl;
    //  cerr << "CURRENT " << psc.energy << endl;
    //  cerr << "LOOK start " << std::distance(l_supercell_sets.begin(),i_start) << " " << (*i_start).getUFFEnergy() << endl;
    //  cerr << "LOOK end " << std::distance(l_supercell_sets.begin(),i_end) << " " << (*i_end).getUFFEnergy() << endl;
    if(i_start==i_end) //std::distance(i_start,i_end)==0)
    { //CO20200106 - patching for auto-indenting
      //    cerr << "start==stop" << endl;
      //
      if(i_start==l_supercell_sets.end()){--i_start;}

      if(areEquivalentStructuresByUFF(i_start,psc)) //aurostd::isequal(psc.energy,(*i_start).getUFFEnergy(),m_energy_uff_tolerance))
      { //CO20200106 - patching for auto-indenting
        //      cerr << "found degeneracy with " << std::distance(l_supercell_sets.begin(),i_start) << endl;
        (*i_start).m_psc_set.push_back(psc);
        //(*i_start).m_degeneracy++;
        return;
      }

      if(psc.m_energy_uff>(*i_start).getUFFEnergy()){++i_start;}
      //  if(i_start!=l_supercell_sets.end()){++i_start;}
      //    cerr << "Adding to " << std::distance(l_supercell_sets.begin(),i_start) << endl;
      POccSuperCellSet pscs; pscs.m_psc_set.push_back(psc);
      l_supercell_sets.insert(i_start,pscs);
      //l_supercell_sets.insert(i_start,psc);
      return;
    }
    //  cerr << "start!=stop" << endl;
    std::list<POccSuperCellSet>::iterator i_middle=i_start;
    //std::advance(i_middle,std::distance(l_supercell_sets.begin(),i_start));
    std::advance(i_middle,std::distance(i_start,i_end)/2);
    //  cerr << "Looking at middle=" << std::distance(l_supercell_sets.begin(),i_middle) << " " << (*i_middle).getUFFEnergy() << endl;
    //cerr << "SEE HERE " << i_middle << endl;
    if(areEquivalentStructuresByUFF(i_middle,psc)) //aurostd::isequal(psc.energy,(*i_middle).getUFFEnergy(),m_energy_uff_tolerance))
    { //CO20200106 - patching for auto-indenting
      //    cerr << "found degeneracy with " << std::distance(l_supercell_sets.begin(),i_middle) << endl;
      (*i_middle).m_psc_set.push_back(psc);
      //(*i_middle).m_degeneracy++;
      return;
    }
    if(psc.m_energy_uff<(*i_middle).getUFFEnergy()){
      //    cerr << "less than middle" << endl;
      return add2DerivativeStructuresList(psc,i_start,i_middle);
    }
    //  cerr << "greater than middle" << endl;
    return add2DerivativeStructuresList(psc,++i_middle,i_end);
  }

  //ORIGINAL
  //void POccCalculator::add2DerivativeStructuresList(ABCDStructureSuper& pds,
  //    std::list<ABCDStructureSuper>::iterator i_start,
  //    std::list<ABCDStructureSuper>::iterator i_end){
  ////  cerr << "FULL LIST START" << endl;
  ////  for(std::list<ABCDStructureSuper>::iterator it=l_derivative_structures.begin();it!=l_derivative_structures.end();++it){
  ////    cerr << std::fixed << std::setprecision(9) << (*it).energy << " ";
  ////  }
  ////  cerr << endl;
  ////  cerr << "FULL LIST END" << endl;
  ////  cerr << "CURRENT " << pds.energy << endl;
  ////  cerr << "LOOK start " << std::distance(l_derivative_structures.begin(),i_start) << " " << (*i_start).energy << endl;
  ////  cerr << "LOOK end " << std::distance(l_derivative_structures.begin(),i_end) << " " << (*i_end).energy << endl;
  //  if(i_start==i_end){ //std::distance(i_start,i_end)==0){
  ////    cerr << "start==stop" << endl;
  //    if(aurostd::isequal(pds.energy,(*i_start).energy,m_energy_uff_tolerance)){
  ////      cerr << "found degeneracy with " << std::distance(l_derivative_structures.begin(),i_start) << endl;
  //      (*i_start).m_degeneracy++;
  //      return;
  //    }
  //    if(pds.energy>(*i_start).energy){
  //      if(i_start!=l_derivative_structures.end()){++i_start;}
  //    }
  ////    cerr << "Adding to " << std::distance(l_derivative_structures.begin(),i_start) << endl;
  //    l_derivative_structures.insert(i_start,pds);
  //    return;
  //  }
  ////  cerr << "start!=stop" << endl;
  //  std::list<ABCDStructureSuper>::iterator i_middle=i_start;
  //  //std::advance(i_middle,std::distance(l_derivative_structures.begin(),i_start));
  //  std::advance(i_middle,std::distance(i_start,i_end)/2);
  ////  cerr << "Looking at middle=" << std::distance(l_derivative_structures.begin(),i_middle) << " " << (*i_middle).energy << endl;
  //  //cerr << "SEE HERE " << i_middle << endl;
  //  if(aurostd::isequal(pds.energy,(*i_middle).energy,m_energy_uff_tolerance)){
  ////    cerr << "found degeneracy with " << std::distance(l_derivative_structures.begin(),i_middle) << endl;
  //    (*i_middle).m_degeneracy++;
  //    return;
  //  }
  //  if(pds.energy<(*i_middle).energy){
  ////    cerr << "less than middle" << endl;
  //    return add2DerivativeStructuresList(pds,i_start,i_middle);
  //  }
  ////  cerr << "greater than middle" << endl;
  //  return add2DerivativeStructuresList(pds,++i_middle,i_end);
  //}
  //}

  void POccCalculator::add2DerivativeStructuresList(const POccSuperCell& psc){
    //  cerr << "SIZE " << l_supercell_sets.size() << endl;
    if(l_supercell_sets.size()==0){
      POccSuperCellSet pscs; pscs.m_psc_set.push_back(psc);
      l_supercell_sets.push_back(pscs);
      return;
    }
    return add2DerivativeStructuresList(psc,l_supercell_sets.begin(),l_supercell_sets.end());
  }

  void POccCalculator::getHNFMatSiteConfig(const POccSuperCell& psc,
      xmatrix<double>& _hnf_mat,
      vector<vector<int> >& _v_types_config){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " psc.m_hnf_index=" << psc.m_hnf_index << endl;
      cerr << __AFLOW_FUNC__ << " psc.m_site_config_index=" << psc.m_site_config_index << endl;
      cerr << __AFLOW_FUNC__ << " psc.m_energy_uff=" << psc.m_energy_uff << endl;
    }
    unsigned long long int hnf_index=0;
    unsigned long long int site_config_index=0;
    //add an initialization flag for energy_analyzer eventually, and check here
    resetHNFMatrices();
    resetSiteConfigurations();
    while(iterateHNFMatrix()){  //sets hnf_mat
      if(hnf_index!=psc.m_hnf_index){hnf_index++;continue;}
      //hnf_mat set here
      _hnf_mat=hnf_mat;
      while(getNextSiteConfiguration()){  //sets v_types_config
        if(site_config_index!=psc.m_site_config_index){site_config_index++;continue;}
        //v_types_config set here
        _v_types_config=v_types_config;
        return;
      }
    }
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid hnf site configuration indices");
  }

  bool POccCalculator::areEquivalentStructures(const POccSuperCell& psc_a,const POccSuperCell& psc_b) {
    //ABCDStructureSuper pds_a,pds_b;
    //pds_a=getParticularABCDStructureSuper(psc_a);
    //pds_b=getParticularABCDStructureSuper(psc_b);

    xstructure a=createXStructure(psc_a,n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
    xstructure b=createXStructure(psc_b,n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
    return areEquivalentStructures(a,b);
  }

  bool POccCalculator::areEquivalentStructures(const xstructure& a, const xstructure& b) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " comparing structure:" << endl;
      cerr << a;
      cerr << "========================================== vs. ==========================================" << endl;
      cerr << b;
    }
    //cerr << aa << endl;
    //cerr << bb << endl;
    bool are_equivalent=compare::structuresMatch(a,b,true,false,false); //match species and use fast match, but not scale volume, two structures with different volumes (pressures) are different! //DX20180123 - added fast_match = true //DX20190318 - not fast_match but optimized_match=false
    //cerr << are_equivalent << endl;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structures are " << (are_equivalent?"":"NOT ") << "equivalent" << endl;}
    return are_equivalent;
  }

  unsigned long long int POccCalculator::runRobustStructureComparison(std::list<POccSuperCellSet>::iterator it_pscs){
    //NOTE, this function needs UFF comparisons first
    //it assumes distinct structures by UFF are indeed distinct, and so no comparisons are made cross-bin, only within bins
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    const vector<POccSuperCell>& vpsc=(*it_pscs).m_psc_set;
    vector<vector<uint> > unique_structure_bins;
    unique_structure_bins.push_back(vector<uint>(0));
    unique_structure_bins.back().push_back(0);  //initialize with first structure

    //when you parallelize this, make a new structure that holds bool, equivalent or not
    //and then run through remaining set of structures

    unsigned long long int i_pscs=std::distance(l_supercell_sets.begin(),it_pscs);
    message << "Performing robust structure comparison for structure group[" << i_pscs+1 << "/" << std::distance(l_supercell_sets.begin(),l_supercell_sets.end()) << "]";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); 

    bool are_equivalent=false,found_equivalent=false;
    bool test_iterator_insertion=false; //short circuit
    xstructure a,b;
    uint starting_index=1;
    pflow::updateProgressBar(0,vpsc.size()-starting_index,*p_oss);
    for(uint i=starting_index;i<vpsc.size();i++){
      const POccSuperCell& psc_b=vpsc[i];
      b=createXStructure(psc_b,n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
      found_equivalent=false;
      for(uint j=0;j<unique_structure_bins.size()&&!found_equivalent;j++){
        if(unique_structure_bins[j].size()==0){
          message << "No structure groups found";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
        }
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " comparing structure[" << i << "] with structure[bin=" << j << ",i=" << 0 << "] (max=" << vpsc.size() << ")" << endl;}
        const POccSuperCell& psc_a=vpsc[unique_structure_bins[j][0]];
        a=createXStructure(psc_a,n_hnf,hnf_count,types_config_permutations_count,true,false);  //PRIMITIVIZE==false, in general it is faster to find whether two structures are equivalent than it is to find primitive cell
        are_equivalent=(!test_iterator_insertion && areEquivalentStructures(a,b));//vpsc[unique_structure_bins[j][0]],vpsc[i]);
        if(are_equivalent){
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " found equivalent structure" << endl;}
          unique_structure_bins[j].push_back(i);
          found_equivalent=true;
          break;
        }
        //if(!are_equivalent){cerr << "AHHHHHHHHHHHHHHHHHH" << endl;}
        //if(j==unique_structure_bins.size()-1){  //we've reached the end and no match, need to make a new bin
        //  unique_structure_bins.push_back(vector<uint>(0));
        //  unique_structure_bins.back().push_back(i);
        //}
      }
      if(!found_equivalent){
        unique_structure_bins.push_back(vector<uint>(0));
        unique_structure_bins.back().push_back(i);
      }
      pflow::updateProgressBar(i,vpsc.size(),*p_oss);
    }

    //test of stupidity
    uint count_structures=0;
    for(uint i=0;i<unique_structure_bins.size();i++){
      count_structures+=unique_structure_bins[i].size();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " unique_structure_bins[" << i << "]=" << unique_structure_bins[i].size() << endl;}
    }
    if(count_structures!=vpsc.size()){
      message << "count_structures == " << count_structures << " != " << vpsc.size() << " == vpsc.size()";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
    }

    if(unique_structure_bins.size()==1){return unique_structure_bins.size();}  //they are all equivalent as dictated by UFF
    //otherwise, we need to rearrange structures

    message << "Splitting structure group[" << std::distance(l_supercell_sets.begin(),it_pscs)+1;
    message << "] into " << unique_structure_bins.size() << " groups as determined by the robust structure comparison";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_); 

    vector<uint> vi_psc_to_remove;
    std::list<POccSuperCellSet>::iterator it=it_pscs;
    it++; //insert continuously at this iterator, do NOT increment inside loop
    for(uint i=1;i<unique_structure_bins.size();i++){
      POccSuperCellSet pscs;
      for(uint j=0;j<unique_structure_bins[i].size();j++){
        pscs.m_psc_set.push_back(vpsc[unique_structure_bins[i][j]]);
        vi_psc_to_remove.push_back(unique_structure_bins[i][j]);
      }
      if(pscs.m_psc_set.size()){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " pre-insertion it=" << std::distance(l_supercell_sets.begin(),it) << endl;}
        l_supercell_sets.insert(it,pscs);
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " post-insertion it=" << std::distance(l_supercell_sets.begin(),it) << endl;}
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " it_pscs=" << std::distance(l_supercell_sets.begin(),it_pscs) << endl;
          cerr << __AFLOW_FUNC__ << " NEW l_supercell_sets.size()=" << l_supercell_sets.size() << endl;
        }
      }
    }

    std::sort(vi_psc_to_remove.rbegin(),vi_psc_to_remove.rend()); //sort in reverse for removal

    for(uint i=0;i<vi_psc_to_remove.size();i++){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " removing structure " << vi_psc_to_remove[i] << " from bin[" << std::distance(l_supercell_sets.begin(),it_pscs)+1 << "]" << endl;}
      (*it_pscs).m_psc_set.erase((*it_pscs).m_psc_set.begin()+vi_psc_to_remove[i]);  //remove!
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " new bin size=" << (*it_pscs).m_psc_set.size() << endl;}

    return unique_structure_bins.size();
  }

  /// @brief thread function which calculates the UFF energy for a vector of POccSuperCell
  /// @param thread_id dummy variable
  /// @param vpsc vector of POccSuperCell
  /// @param v_energy_analyzer vector of POccUFFEnergyAnalyzer
  /// @param vv_types_config vector of v_types_config
  /// @param npsc_queue number of POccSuperCell left to calculate
  /// @param m_save save mutex
  /// @param m_job job mutex
  /// @note two mutexes are needed, m_job locks when determining the index of the POccSuperCell, m_save locks when writing to the vector of POccSuperCell
  ///
  /// @authors
  /// @mod{SD+HE,20230604,created}
  void POccCalculator::calculatePOccSuperCellUFF(int thread_id, vector<POccSuperCell>& vpsc, const vector<POccUFFEnergyAnalyzer>& v_energy_analyzer, const vector<vector<vector<int>>>& vv_types_config, size_t& npsc_queue, std::mutex& m_save, std::mutex& m_job){
    (void)thread_id;
    size_t ipsc = 0;
    std::map<uint, POccUFFEnergyAnalyzer> local_ea;
    const size_t step_size = std::max((size_t)1, vpsc.size() / 100);
    while(true){
      {
        std::lock_guard<std::mutex> lg_job(m_job);
        if (npsc_queue <= 0){break;}
        npsc_queue--;
        ipsc = npsc_queue;
        if(npsc_queue % step_size == 0){pflow::updateProgressBar(vpsc.size() - npsc_queue, vpsc.size(), *p_oss);} // use a step size to avoid updating the progress bar too often
      }
      if(local_ea.count(vpsc[ipsc].m_hnf_index) == 0){
        local_ea[vpsc[ipsc].m_hnf_index] = v_energy_analyzer[vpsc[ipsc].m_hnf_index];
      }
      local_ea[vpsc[ipsc].m_hnf_index].setBonds(vv_types_config[ipsc]);
      double uff_energy = aurostd::round(local_ea[vpsc[ipsc].m_hnf_index].getUFFEnergy(), 12);
      std::lock_guard<std::mutex> lg_save(m_save);
      vpsc[ipsc].m_energy_uff = uff_energy;
    }
  }

  /// @brief thread function which counts the number of unique POCC tiles using the UFF energy and an unordered map
  /// @param thread_id dummy variable
  /// @param map_unique unordered map of the UFF energies
  /// @param vpsc vector of POccSuperCell
  /// @param v_energy_analyzer vector of POccUFFEnergyAnalyzer
  /// @param vv_types_config vector of v_types_config
  /// @param npsc_queue number of POccSuperCell left to calculate
  /// @param m_save save mutex
  /// @param m_job job mutex
  /// @note two mutexes are needed, m_job locks when determining the index of the POccSuperCell, m_save locks when writing to the vector of POccSuperCell
  ///
  /// @authors
  /// @mod{SD+HE,20230609,created}
  void POccCalculator::countUniquePOccSuperCellUFF(int thread_id, std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>& map_unique, const vector<POccSuperCell>& vpsc, const vector<POccUFFEnergyAnalyzer>& v_energy_analyzer, const vector<vector<vector<int>>>& vv_types_config, size_t& npsc_queue, std::mutex& m_save, std::mutex& m_job){
    (void)thread_id;
    const double uff_tol = DEFAULT_UFF_ENERGY_TOLERANCE;
    size_t ipsc = 0;
    std::map<uint, POccUFFEnergyAnalyzer> local_ea;
    std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>> local_mu;
    const size_t step_size = std::max((size_t)1, vpsc.size() / 100);
    while(true){
      { 
        std::lock_guard<std::mutex> lg_job(m_job);
        if (npsc_queue <= 0){break;}
        npsc_queue--;
        ipsc = npsc_queue;
        if (npsc_queue % step_size == 0){pflow::updateProgressBar(vpsc.size() - npsc_queue, vpsc.size(), *p_oss);} // use a step size to avoid updating the progress bar too often
      }
      if (local_ea.count(vpsc[ipsc].m_hnf_index) == 0){
        local_ea[vpsc[ipsc].m_hnf_index] = v_energy_analyzer[vpsc[ipsc].m_hnf_index];
      }
      local_ea[vpsc[ipsc].m_hnf_index].setBonds(vv_types_config[ipsc]);
      unsigned long int uff_energy = std::ceil(local_ea[vpsc[ipsc].m_hnf_index].getUFFEnergy() / uff_tol);
      local_mu[vpsc[ipsc].m_hnf_index][uff_energy]++;
    }
    std::lock_guard<std::mutex> lg_save(m_save);
    for(std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>::iterator it_out = local_mu.begin(); it_out != local_mu.end(); it_out++){
      for(std::unordered_map<unsigned long int, unsigned long int>::iterator it = it_out->second.begin(); it != it_out->second.end(); it++){
        map_unique[it_out->first][it->first] += it->second;
      }
    }
  }

  void POccCalculator::calculate(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    //starting criteria for HNF matrices
    //xmatrix<double> hnf_mat;                            //really xmatrix of int's, but we rule for int * double in xmatrix, no big deal
    //vector<xmatrix<double> > v_unique_superlattices;    //only store locally, as we need to make comparisons

    //starting criteria for site combinations
    //vector<vector<int> > v_types_config;
    //vector<int> v_config_iterators, v_site_iterators;

    //debug
    //for(uint i=0;i<vv_count_configs.size();i++){vv_count_configs[i].clear();}
    //vv_count_configs.clear();
    //vv_count_configs.push_back(vector<POccSiteConfiguration>(0));
    //vv_count_configs.push_back(vector<POccSiteConfiguration>(0));
    //POccSiteConfiguration test;
    //test.types_configuration_debug.push_back(-1);
    //test.types_configuration_debug.push_back(-1);
    //test.types_configuration_debug.push_back(0);
    //test.types_configuration_debug.push_back(0);
    //vv_count_configs[0].push_back(test);
    //test.types_configuration_debug.clear();
    //test.types_configuration_debug.push_back(1);
    //test.types_configuration_debug.push_back(1);
    //test.types_configuration_debug.push_back(1);
    //test.types_configuration_debug.push_back(2);
    //vv_count_configs[0].push_back(test);
    //test.types_configuration_debug.clear();
    //test.types_configuration_debug.push_back(3);
    //test.types_configuration_debug.push_back(3);
    //test.types_configuration_debug.push_back(4);
    //test.types_configuration_debug.push_back(4);
    //vv_count_configs[1].push_back(test);
    //test.types_configuration_debug.clear();
    //test.types_configuration_debug.push_back(5);
    //test.types_configuration_debug.push_back(6);
    //test.types_configuration_debug.push_back(7);
    //test.types_configuration_debug.push_back(8);
    //vv_count_configs[1].push_back(test);
    //for(uint i=0;i<vv_count_configs.size();i++){
    //for(uint j=0;j<vv_count_configs[i].size();j++){
    //cerr << "site " << i << " config " << j <<  endl;
    //for(uint k=0;k<vv_count_configs[i][j].types_configuration_debug.size();k++){
    //cerr << vv_count_configs[i][j].types_configuration_debug[k] << " ";
    //}
    //cerr << endl;
    //}
    //}

    //getTotalPermutationsCount();  //done separately

    //POccUFFEnergyAnalyzer energy_analyzer;
    //const xstructure& xstr_nopocc = p_str.xstr_nopocc;

    //xmatrix<double> distance_matrix;
    //vector<vector<uint> > v_bonded_atom_indices, v_nonbonded_atom_indices; 

    //BIG REVELATION: nearest neighbor distances change with vacancies! 
    //this is a problem, because we need to recalculate bonding for every new site configuration
    //if no vacancies, only perform bonding analysis once
    //if no vacancies, we can extrapolate nn distances from primitive cell
    //need to set two booleans, one for if vacancies, and if bonding set
    //if vacancies, bonding set won't matter, as we redo every time
    //remember, wait to do minkowski reduction until the end
    //it does reduce the size of the cluster, but it also changes the number of atoms (BAD for v_types_config)
    //the cluster only needs to be calculated once, from there we can determine each new bonding
    //make functions for getting nn distances and then determining bonding, which simply returns if no vacancies and already set

    message << "Calculating unique supercells. Please be patient";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    ////[CO20181226 - need to experiment]double energy_radius=RadiusSphereLattice(getLattice())*1.5; //figure this out CO, only works if energy radius = 10
    //double energy_radius=DEFAULT_UFF_CLUSTER_RADIUS;
    ////if(SET_KESONG_STANDARD_DIST){energy_radius=ENERGY_RADIUS;}
    //if(COMPARE_WITH_KESONG){energy_radius=ENERGY_RADIUS;}

    ////energy_radius=12;

    //energy_analyzer.initialize(xstr_pocc,xstr_nopocc,m_species_redecoration,m_p_flags,m_aflags,*p_FileMESSAGE,*p_oss); //p_str,energy_radius);
    //energy_analyzer.m_species_redecoration=m_species_redecoration;  //clean me up, fix how we integrate these two objects
    //energy_analyzer.types2uffparams_map=getTypes2UFFParamsMap(m_species_redecoration); //clean me up, fix how we integrate these two objects

    unsigned long long int hnf_index=0;
    unsigned long long int site_config_index;

    //get algorithm settings
    string struct_gen_algo=DEFAULT_POCC_STRUCTURE_GENERATION_ALGO;

    if(struct_gen_algo=="UFF"){
      POccSuperCell psc;
      vector<POccSuperCell> vpsc;
      vector<POccUFFEnergyAnalyzer> v_energy_analyzer;
      vector<vector<vector<int>>> vv_types_config;
      energy_analyzer.initialize(xstr_pocc,xstr_nopocc,m_species_redecoration,m_p_flags,m_aflags,*p_FileMESSAGE,*p_oss); //p_str,energy_radius);
      resetHNFMatrices();
      while(iterateHNFMatrix()){
        psc.m_hnf_index = hnf_index;
        psc.m_degeneracy = 1;
        resetSiteConfigurations();
        site_config_index = 0;
        POccUFFEnergyAnalyzer ea(energy_analyzer);
        ea.getCluster(hnf_mat);
        v_energy_analyzer.push_back(ea);
        while(getNextSiteConfiguration()){
          psc.m_site_config_index=site_config_index;
          vpsc.push_back(psc);
          vv_types_config.push_back(v_types_config);
          site_config_index++;
        }
        hnf_index++;
      }
      std::mutex m_save, m_job;
      xthread::xThread xt(KBIN::get_NCPUS(m_kflags),1);
      size_t npsc_queue = vpsc.size();
      if(m_count_unique_fast){
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Count POCC unique FAST", m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        const double uff_tol = DEFAULT_UFF_ENERGY_TOLERANCE;
        std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>> map_unique;
        std::function<void(int,
            std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>&,
            const vector<POccSuperCell>&,
            const vector<POccUFFEnergyAnalyzer>&,
            const vector<vector<vector<int>>>&,
            size_t&,
            std::mutex&,
            std::mutex&)> fn=std::bind(&POccCalculator::countUniquePOccSuperCellUFF,this,
          std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5,std::placeholders::_6,std::placeholders::_7,std::placeholders::_8);
        pflow::updateProgressBar(0, vpsc.size(), *p_oss);
        xt.run(KBIN::get_NCPUS(m_kflags), fn, map_unique, vpsc, v_energy_analyzer, vv_types_config, npsc_queue, m_save, m_job);
        unsigned long long int count_tot = 0, count_bin = 0;
        for(std::map<unsigned long long int, std::unordered_map<unsigned long int, unsigned long int>>::iterator it_out = map_unique.begin(); it_out != map_unique.end(); it_out++){
          cout << "\tHFN index " << it_out->first << endl;
          for(std::unordered_map<unsigned long int, unsigned long int>::iterator it = it_out->second.begin(); it != it_out->second.end(); it++){
            count_bin++;
            cout << "\t\tStructure bin " << count_bin << ": energy=" << (it->first) * uff_tol << ", " << "degeneracy=" << it->second << endl;
            count_tot += it->second;
          }
        }
        cout << "Total structures: " << count_tot << endl;
        cout << "Total unique structures: " << count_bin << endl;
        return;
      }
      else{
        std::function<void(int,
            vector<POccSuperCell>&,
            const vector<POccUFFEnergyAnalyzer>&,
            const vector<vector<vector<int>>>&,
            size_t&,
            std::mutex&,
            std::mutex&)> fn=std::bind(&POccCalculator::calculatePOccSuperCellUFF,this,
          std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5,std::placeholders::_6,std::placeholders::_7);
        // the first number sets the number of threads that are created overall (in this case one thread per CPU)
        pflow::updateProgressBar(0, vpsc.size(), *p_oss);
        xt.run(KBIN::get_NCPUS(m_kflags), fn, vpsc, v_energy_analyzer, vv_types_config, npsc_queue, m_save, m_job);
        for(size_t i = 0; i < vpsc.size(); i++){add2DerivativeStructuresList(vpsc[i]);}
      }

    }else{  //group theory approach

      uint i=0,j=0,k=0,site=0,occ=0,permut=0;
      int gi=0;

      vector<xmatrix<double> > vhnf_mats;
      vector<vector<vector<int> > > vv_types_config;
      resetHNFMatrices();
      while(iterateHNFMatrix()){vhnf_mats.push_back(hnf_mat);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vhnf_mats.size()=" << vhnf_mats.size() << endl;}

      xmatrix<long long int> H,U,V,S;
      xmatrix<double> V_double;
      vector<xmatrix<long long int> > vS_uniq;
      vector<uint> viH_uniqSmith;
      bool match=false;
      for(i=0;i<vhnf_mats.size();i++){
        H=aurostd::xmatrixdouble2utype<long long int>(vhnf_mats[i]);
        getSmithNormalForm(H,U,V,S);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " H=" << endl;cerr << H << endl;
          cerr << __AFLOW_FUNC__ << " S=" << endl;cerr << S << endl;
        }
        match=false;
        for(j=0;j<vS_uniq.size() && match==false;j++){
          if(aurostd::isequal(S,vS_uniq[j])){match=true;}
        }
        if(match==false){
          vS_uniq.push_back(S);
          viH_uniqSmith.push_back(i);
        }
      }

      if(LDEBUG){
        cerr << __AFLOW_FUNC__ << " vS_uniq.size()=" << vS_uniq.size() << endl;
        cerr << __AFLOW_FUNC__ << " viH_uniqSmith.size()=" << viH_uniqSmith.size() << endl;
        for(i=0;i<vS_uniq.size();i++){
          cerr << __AFLOW_FUNC__ << " viH_uniqSmith[i=" << i <<"]=" << endl;cerr << viH_uniqSmith[i] << endl;
          cerr << __AFLOW_FUNC__ << " vS_uniq[i=" << i <<"]=" << endl;cerr << vS_uniq[i] << endl;
        }
      }

      xstructure xstr_ss;
      vector<int> sc2pc_map,pc2sc_map;
      xvector<double> fpos,d,x;
      vector<GroupMember> vg;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr_nopocc.lattice=" << endl;cerr << xstr_nopocc.lattice << endl;}
      for(i=0;i<vhnf_mats.size();i++){
        H=aurostd::xmatrixdouble2utype<long long int>(vhnf_mats[i]);
        V_double=aurostd::xmatrixutype2double(V);
        getSmithNormalForm(H,U,V,S);
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " hnf_mat=" << endl;cerr << vhnf_mats[i] << endl;
          cerr << __AFLOW_FUNC__ << " H=" << endl;cerr << H << endl;
          cerr << __AFLOW_FUNC__ << " U=" << endl;cerr << U << endl;
          cerr << __AFLOW_FUNC__ << " V=" << endl;cerr << V << endl;
        }
        xstr_ss=GetSuperCell(xstr_nopocc,vhnf_mats[i],sc2pc_map,pc2sc_map,false,false,false,false);  //force_strict_pc2scMap==false is very important here, we want the FIRST equivalent atom
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " xstr_ss.lattice=" << endl;cerr << xstr_ss.lattice << endl;
          cerr << __AFLOW_FUNC__ << " sc2pc_map=" << aurostd::joinWDelimiter(sc2pc_map,",") << endl;
          cerr << __AFLOW_FUNC__ << " pc2sc_map=" << aurostd::joinWDelimiter(pc2sc_map,",") << endl;
        }
        vg.clear();
        for(j=0;j<sc2pc_map.size();j++){
          vg.push_back(GroupMember());
          //create map h from parent lattice onto the group (g are elements of group)
          //h(x)=[UA^{-1}x]_S:
          //c2f is A^{-1} (done for fpos)
          vg.back().basis=xstr_nopocc.atoms[sc2pc_map[j]].basis;
          fpos=xstr_nopocc.c2f*xstr_ss.atoms[j].cpos; //coordinates of sc_atom[j] in basis of primitive lattice
          //fpos=x+d
          //where d is fraction part and x is integer
          //x is the lattice
          //d is the basis
          x=aurostd::floor(fpos); //int only
          d=aurostd::mod(fpos,1.0); //between 0 and 1
          if(LDEBUG){
            //cerr << __AFLOW_FUNC__ << " atom[pc=" << sc2pc_map[j] << ",sc=" << j << "].ijk=" << xstr_ss.atoms[j].ijk << endl;
            //cerr << __AFLOW_FUNC__ << " atom[pc=" << sc2pc_map[j] << "].cpos=" << xstr_nopocc.atoms[sc2pc_map[j]].cpos << endl;
            //cerr << __AFLOW_FUNC__ << " atom[sc=" << j << "].cpos=" << xstr_ss.atoms[j].cpos << endl;
            cerr << __AFLOW_FUNC__ << " fpos=" << fpos << endl;
            cerr << __AFLOW_FUNC__ << " d=" << d << endl;
            cerr << __AFLOW_FUNC__ << " x=" << x << endl;
          }
          //g is the group element
          vg.back().cpos=x*V_double;
          //mod: translation symmetry (next cell)
          for(gi=vg.back().cpos.lrows;gi<=vg.back().cpos.urows;gi++){
            vg.back().cpos[gi]=aurostd::mod(vg.back().cpos[gi],(double)S[gi][gi]);
          }
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " g.basis=" << vg.back().basis << endl;
            cerr << __AFLOW_FUNC__ << " g.cpos=" << vg.back().cpos << endl;
          }
        }
      }

      resetSiteConfigurations();
      while(getNextSiteConfiguration()){vv_types_config.push_back(v_types_config);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vv_types_config.size()=" << vv_types_config.size() << endl;}

      //eliminate translation duplicates
      vector<int> types_config0,types_config1,types_config2;
      uint natoms_super=(uint)n_hnf*xstr_nopocc.atoms.size();
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " natoms_super=" << natoms_super << endl;}
      if(natoms_super==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"natoms_super==0",_RUNTIME_ERROR_);}
      types_config0.resize(natoms_super);types_config1.resize(natoms_super);types_config2.resize(natoms_super);

      //aurostd::xcombos xc;  //this is NOT an xcombos permutations, these are 1->2, 2->3, strictly
      bool found=false;
      for(i=0;i<vv_types_config.size();i++){  //do NOT use a constant vv_types_config.size(), it will change within loop
        //build types_config0
        k=0;
        for(site=0;site<vv_types_config[i].size();site++){
          for(occ=0;occ<vv_types_config[i][site].size();occ++){
            types_config0[k++]=vv_types_config[i][site][occ];
          }
        }
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " types_config0=" << aurostd::joinWDelimiter(types_config0,",") << endl;}
        for(j=vv_types_config.size()-1;j>i;j--){  //do NOT use a constant vv_types_config.size(), it will change within loop
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " comparing i=" << i << " with j=" << j << endl;}
          //build types_config2
          k=0;
          for(site=0;site<vv_types_config[j].size();site++){
            for(occ=0;occ<vv_types_config[j][site].size();occ++){
              types_config2[k++]=vv_types_config[j][site][occ];
            }
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " types_config2=" << aurostd::joinWDelimiter(types_config2,",") << endl;}
          found=false;
          //build types_config1
          for(k=0;k<types_config0.size();k++){types_config1[k]=types_config0[k];}
          for(permut=0;permut<natoms_super;permut++){
            if(permut>0){ //actually permute
              std::rotate(types_config1.begin(),types_config1.begin()+1,types_config1.end());
            }
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " comparing:" << endl;
              cerr << __AFLOW_FUNC__ << " types_config1=" << aurostd::joinWDelimiter(types_config1,",") << endl;
              cerr << __AFLOW_FUNC__ << " types_config2=" << aurostd::joinWDelimiter(types_config2,",") << endl;
            }
            if(types_config1==types_config2){found=true;break;} //vector<int> can be compared exactly
          }
          if(found==true){
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " erasing j=" << j << endl;}
            vv_types_config.erase(vv_types_config.begin()+j);
            continue;
          }
        }
      }

      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vv_types_config.size()=" << vv_types_config.size() << endl;}

      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"approach not complete yet",_RUNTIME_ERROR_);
    }

    //ORIGINAL
    //ABCDStructureSuper pds;
    //pds.reset();
    //resetHNFMatrices();
    ////cerr << "HERE2" << endl;
    //while(iterateHNFMatrix()){
    //  //cerr << "HERE3" << endl;
    //  energy_analyzer.getCluster(hnf_mat);
    //  pds.m_hnf_index=hnf_index;
    //  //cerr << "HERE4" << endl;
    //  pds.hnf_mat=hnf_mat;
    //  resetSiteConfigurations();
    //  site_config_index=0;
    //  //cerr << "HERE5" << endl;
    //  while(getNextSiteConfiguration()){
    //    //cerr << "HERE6" << endl;
    //    pds.m_site_config_index=site_config_index;
    //    //cerr << pds.site_config_index << endl;
    //    pds.v_types_config=v_types_config;
    //    //cerr << "HERE7" << endl;
    //    getBonds();
    //    //cerr << "HERE8" << endl;
    //    pds.energy=getUFFEnergy();
    //    pds.m_degeneracy=1;
    //    //cerr << "HERE9" << endl;
    //    add2DerivativeStructuresList(pds);
    //    //cerr << "DONE ADDING " << endl;
    //    site_config_index++;
    //		pflow::updateProgressBar(++current_iteration,total_permutations_count,*p_oss);
    //    //cerr << "HERE10" << endl;
    //  }
    //  hnf_index++;
    //}

    //tests of stupidity - START
    if(!l_supercell_sets.size()){
      message << "No supercells detected by UFF. Please run calculate() to determine unique supercells";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
      //pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
      //cerr << "probably broken pointers" << endl;
    }


    if(DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON){
      bool DEBUG_RSC=false;
      message << "Calculated " << l_supercell_sets.size() << " unique supercell groups via UFF, now performing robust structure comparison";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);//_LOGGER_COMPLETE_);

      unsigned long long int count_new_groups;
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        if(DEBUG_RSC||LDEBUG){cerr << __AFLOW_FUNC__ << " current it=" << std::distance(l_supercell_sets.begin(),it)+1 << "/" << std::distance(l_supercell_sets.begin(),l_supercell_sets.end()) << " [1]" << endl;}
        POccSuperCellSet& pscs=(*it);
        std::sort(pscs.m_psc_set.begin(),pscs.m_psc_set.end(),sortPSCsUFFEnergy); //sort by energies (should all be about the same, ensures splits via structure comparison stay sorted by energy
        count_new_groups=runRobustStructureComparison(it);  //l_supercell_sets may change size, get count of the increase
        if(DEBUG_RSC||LDEBUG){cerr << __AFLOW_FUNC__ << " advancing by " << count_new_groups-1 << endl;}
        std::advance(it,count_new_groups-1);  //adjust iterator appropriately
        if(DEBUG_RSC||LDEBUG){cerr << __AFLOW_FUNC__ << " current it=" << std::distance(l_supercell_sets.begin(),it)+1 << "/" << std::distance(l_supercell_sets.begin(),l_supercell_sets.end()) << " [2]" << endl;}
        if(DEBUG_RSC||LDEBUG){cerr << __AFLOW_FUNC__ << " it=" << std::distance(l_supercell_sets.begin(),it)  << " vs. l_supercell_sets.end()=" << std::distance(l_supercell_sets.begin(),l_supercell_sets.end()) << endl;}
      }
    }

    message << "Calculated " << l_supercell_sets.size() << " unique supercells";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_COMPLETE_);

    unsigned long long int total_degeneracy=0;

    stringstream message_prec;
    int prec=(int)ceil(log10(1.0/m_energy_uff_tolerance));
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      //cerr << "UNIQUE " << std::distance(l_supercell_sets.begin(),it) << endl;
      message_prec << std::fixed << std::setprecision(prec) << "Structure bin "  << std::distance(l_supercell_sets.begin(),it)+1 << ": energy=" << (*it).getUFFEnergy() << ", " << "degeneracy=" << (*it).getDegeneracy();
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message_prec,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      //cerr << endl;
      total_degeneracy+=(*it).getDegeneracy();
    }

    if(total_permutations_count!=total_degeneracy){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected degeneracy count (does not match expected total permutations count)");
    }

    if(m_p_flags.flag("POCC_COUNT_UNIQUE")){return;}  //so we don't bombard user with too much verbosity

    message << "Resorting unique supercell order by HNF matrix/site configuration indices";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    //go through sets and sort by hnf/site_config index
    //this ensures that even if we mess with underlying UFF/structure comparison algorithm, the same structure is chosen everytime
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      POccSuperCellSet& pscs=(*it);
      std::sort(pscs.m_psc_set.begin(),pscs.m_psc_set.end()); //sort by indices
    }
    l_supercell_sets.sort();  //finally sort on indices so the list is always fixed

    if(true && !m_aflags.Directory.empty() && !m_p_flags.flag("POCC_SKIP_WRITING_FILES")){ //move me so it runs whenever pocc runs (even post-processing)
      message << "Writing out " << POCC_ALL_HNF_MATRICES_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      stringstream all_hnf_mat_ss;
      resetHNFMatrices();
      hnf_index=0;
      while(iterateHNFMatrix()){
        all_hnf_mat_ss << AFLOWIN_SEPARATION_LINE << endl;
        all_hnf_mat_ss << AFLOW_POCC_TAG << "HNF_MATRIX " << hnf_index+1 << "/" << hnf_count << endl;
        all_hnf_mat_ss << hnf_mat << endl;
        all_hnf_mat_ss << AFLOWIN_SEPARATION_LINE << endl;
        hnf_index++;
      }
      aurostd::stringstream2file(all_hnf_mat_ss,getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_ALL_HNF_MATRICES_FILE);
    }

    if(true && !m_aflags.Directory.empty() && !m_p_flags.flag("POCC_SKIP_WRITING_FILES")){ //move me so it runs whenever pocc runs (even post-processing)
      message << "Writing out " << POCC_ALL_SITE_CONFIGURATIONS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      stringstream all_site_configs_ss;
      resetSiteConfigurations();
      site_config_index=0;
      while(getNextSiteConfiguration()){
        all_site_configs_ss << AFLOWIN_SEPARATION_LINE << endl;
        all_site_configs_ss << AFLOW_POCC_TAG << "SITE_CONFIGURATION " << site_config_index+1 << "/" << types_config_permutations_count << endl;
        for(uint site=0;site<v_types_config.size();site++){
          for(uint atom=0;atom<v_types_config[site].size();atom++){
            all_site_configs_ss << (site==0&&atom==0?"":" ") << aurostd::utype2string(v_types_config[site][atom]);
          }
          all_site_configs_ss << ((site!=v_types_config.size()-1)?";":"");
        }
        all_site_configs_ss << endl;
        const deque<string>& species=xstr_pocc.species;
        for(uint site=0;site<v_types_config.size();site++){
          for(uint atom=0;atom<v_types_config[site].size();atom++){
            if(v_types_config[site][atom]>=0){
              if(v_types_config[site][atom]<(int)species.size()){all_site_configs_ss << (site==0&&atom==0?"":" ") << species[v_types_config[site][atom]];}
              else {all_site_configs_ss << (site==0&&atom==0?"":" ") << "XX";}
            }
            else {all_site_configs_ss << (site==0&&atom==0?"":" ") << "VAC";}
          }
          all_site_configs_ss << ((site!=v_types_config.size()-1)?";":"");
        }
        all_site_configs_ss << endl;
        all_site_configs_ss << AFLOWIN_SEPARATION_LINE << endl;
        site_config_index++;
      }
      aurostd::stringstream2file(all_site_configs_ss,getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_ALL_SITE_CONFIGURATIONS_FILE);
    }

    if(DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS && !m_aflags.Directory.empty() && !m_p_flags.flag("POCC_SKIP_WRITING_FILES")){
      message << "Writing out " << POCC_ALL_SUPERCELLS_FILE << ". Please be patient.";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      stringstream all_supercells_ss;
      //string POSCAR_START=POSCAR_POCC_series_START_tag;
      //string POSCAR_STOP=POSCAR_POCC_series_STOP_tag;
      //stringstream POSCAR_strtag_ss;
      string POSCAR_strtag="";
      POccSuperCell psc;
      unsigned long long int isupercell=0;
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        isupercell=std::distance(l_supercell_sets.begin(),it);
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " isupercell=" << isupercell << endl;}
        const POccSuperCellSet& pscs=(*it);
        //all_supercells_ss << "---------------------------------------------------------------------------------------------" << endl;
        all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
        all_supercells_ss << AFLOW_POCC_TAG << "STRUCTURES_GROUP " << isupercell+1 << "/" << l_supercell_sets.size() << endl;
        all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
        for(uint i=0;i<pscs.m_psc_set.size();i++){
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " i=" << i << endl;}
          psc=pscs.m_psc_set[i];
          psc.m_degeneracy=1;
          POSCAR_strtag=pocc::getARUNString(isupercell,l_supercell_sets.size(),i,pscs.m_psc_set.size(),pscs.m_psc_set[i].m_hnf_index,pscs.m_psc_set[i].m_site_config_index,true);
          all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
          all_supercells_ss << AFLOW_POCC_TAG << "STRUCTURE " << i+1 << "/" << pscs.m_psc_set.size() << endl;
          all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
          //[OBSOLETE]POSCAR_strtag_ss.str("");
          //[OBSOLETE]POSCAR_strtag_ss << std::setfill('0') << std::setw(aurostd::getZeroPadding(l_supercell_sets.size())) << isupercell+1 << "_"; //+1 so we start at 1, not 0 (count)
          //[OBSOLETE]POSCAR_strtag_ss << std::setfill('0') << std::setw(aurostd::getZeroPadding(pscs.m_psc_set.size())) << (i+1) << "_"; //+1 so we start at 1, not 0 (count)
          //[OBSOLETE]POSCAR_strtag_ss << "H";
          //[OBSOLETE]POSCAR_strtag_ss << pscs.m_psc_set[i].m_hnf_index;
          //[OBSOLETE]POSCAR_strtag_ss << "C";
          //[OBSOLETE]POSCAR_strtag_ss << pscs.m_psc_set[i].m_site_config_index;
          all_supercells_ss << AFLOW_POCC_TAG << "UFF_ENERGY=" << std::fixed << std::setprecision(15) << psc.m_energy_uff << endl; all_supercells_ss.unsetf(std::ios_base::floatfield);
          all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
          //[OBSOLETE]all_supercells_ss << POSCAR_POCC_series_START_tag << POSCAR_strtag_ss.str() << endl;
          all_supercells_ss << POSCAR_POCC_series_START_tag << POSCAR_strtag << endl;
          all_supercells_ss << createXStructure(psc,n_hnf,hnf_count,types_config_permutations_count,true,PRIMITIVIZE); // << endl;
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structure created" << endl;}
          //[OBSOLETE]all_supercells_ss << POSCAR_POCC_series_STOP_tag << POSCAR_strtag_ss.str() << endl;
          all_supercells_ss << POSCAR_POCC_series_STOP_tag << POSCAR_strtag << endl;
          all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
        }
        //all_supercells_ss << AFLOWIN_SEPARATION_LINE << endl;
      }
      aurostd::stringstream2file(all_supercells_ss,getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_ALL_SUPERCELLS_FILE);
    }

    //tests of stupidity - END
  }

  unsigned long long int POccCalculator::getUniqueSuperCellsCount() const {return l_supercell_sets.size();}

  //int getZeroPadding(unsigned long long int num){return int(log10(num))+1;} ME20190108 - moved to aurostd

  string POccCalculator::getARUNString(unsigned long long int i) {return pocc::getARUNString(l_supercell_sets,i);}

  string getARUNString(const std::list<POccSuperCellSet>& l_supercell_sets,unsigned long long int i) {
    std::list<POccSuperCellSet>::const_iterator it=l_supercell_sets.begin();
    std::advance(it,i);
    return pocc::getARUNString(i,l_supercell_sets.size(),0,(*it).m_psc_set.size(),(*it).getSuperCell().m_hnf_index,(*it).getSuperCell().m_site_config_index,false);
  }

  string getARUNString(
      unsigned long long int index_structure_group,
      unsigned long long int vstructure_groups_size,
      unsigned long long int index_structure,
      unsigned long long int vstructures_size,
      unsigned long long int index_hnf,
      unsigned long long int index_site_config,
      bool include_strgrp) {
    stringstream ARUN;

    //[CO20181226 - done automatically in AVASP_populateXVASP()]ARUN << "ARUN.POCC_";
    ARUN << std::setfill('0') << std::setw(aurostd::getZeroPadding(vstructure_groups_size)) << index_structure_group+1 << "_";  //+1 so we start at 1, not 0 (count)
    if(include_strgrp){ARUN << std::setfill('0') << std::setw(aurostd::getZeroPadding(vstructures_size)) << index_structure+1 << "_";}
    ARUN << "H";
    ARUN << index_hnf; //-1 because this is treated as a HASH and it starts at 0
    ARUN << "C";
    ARUN << index_site_config; //-1 because this is treated as a HASH and it starts at 0

    //no padding the hash, unnecessary
    //ARUN << "H";
    //ARUN << std::setfill('0') << std::setw(aurostd::getZeroPadding(hnf_count-1)) << (*it).getSuperCell().m_hnf_index; //-1 because this is treated as a HASH and it starts at 0
    //ARUN << "C";
    //ARUN << std::setfill('0') << std::setw(aurostd::getZeroPadding(types_config_permutations_count-1)) << (*it).getSuperCell().m_site_config_index; //-1 because this is treated as a HASH and it starts at 0

    //string ARUN;
    //ARUN="ARUN.POCC"+aurostd::utype2string(i);
    //ARUN+="H"+aurostd::utype2string((*it).getSuperCell().m_hnf_index);//+1);
    //ARUN+="C"+aurostd::utype2string((*it).getSuperCell().m_site_config_index);//+1);
    //return ARUN;

    return ARUN.str();
  }

  xstructure POccCalculator::getUniqueSuperCell(unsigned long long int i) {
    stringstream message;

    if(i>l_supercell_sets.size()-1){
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Invalid supercell index");
    }

    std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();
    std::advance(it,i);

    POccSuperCell psc=(*it).getSuperCell();
    psc.m_degeneracy=(*it).getDegeneracy(); //BEWARE OF DEGENERACY of this special POccSuperCell, representative of all supercells in that set

    return createXStructure(psc,n_hnf,hnf_count,types_config_permutations_count,true,PRIMITIVIZE);
  }

  vector<xstructure> POccCalculator::getUniqueDerivativeStructures() {
    stringstream message;

    vector<xstructure> v_xstr;
    if(!l_supercell_sets.size()){
      message << "No supercells detected. Please run calculate() to determine unique supercells";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
      return v_xstr;
    }

    unsigned long long int current_iteration=0;
    pflow::updateProgressBar(current_iteration,l_supercell_sets.size()-1,*p_oss);
    POccSuperCell psc;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      psc=(*it).getSuperCell();
      psc.m_degeneracy=(*it).getDegeneracy(); //BEWARE OF DEGENERACY of this special POccSuperCell, representative of all supercells in that set
      v_xstr.push_back(createXStructure(psc,n_hnf,hnf_count,types_config_permutations_count,true,PRIMITIVIZE));
      pflow::updateProgressBar(++current_iteration,l_supercell_sets.size()-1,*p_oss);
      //cout << AFLOWIN_SEPARATION_LINE << endl;
      //cout << createXStructure((*it),true);
      //cout << AFLOWIN_SEPARATION_LINE << endl;
    }

    message << "Full list generated";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    return v_xstr;

    //  resetHNFMatrices();
    //  while(iterateHNFMatrix()){
    ////unsigned long long int total_permutations_count=0;
    //  //cerr << "NEW HNF" << endl;
    //  //cerr << hnf_mat << endl;
    //  
    //  v_types_config.clear();
    //  

    //getBonding(hnf_mat,distance_matrix,v_bonded_atom_indices,v_nonbonded_atom_indices);
    //energy_analyzer.getBonding(p_str,hnf_mat);
    //    getBonds();

    //    resetSiteConfigurations();
    //    while(getNextSiteConfiguration()){
    //getSpecificBonding();
    //      getUFFEnergy();
    //getXStructure();
    //    //total_permutations_count++;
    //    //cerr << "Possibility " << total_permutations_count << "    ";
    //    //for(uint i=0;i<v_types_config.size();i++){
    //    //  for(uint j=0;j<v_types_config[i].size();j++){
    //    //    cerr << v_types_config[i][j] << " ";
    //    //  }
    //    //  cerr << "    ";
    //    //}
    //    //cerr << endl;
    //    }
    //
    //  //cerr << "end hnf" << endl << endl;
    //  }



    //ostream& oss=cout;
    //message << "Total number of derivative structure possibilities = " << total_permutations_count;
    //pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);


    //return true;
  }

  // ME20211008 - Map an input structure to the parent PARTCAR
  // It builds this map by creating a map from the input structure
  // to the unrelaxed derivative structure first. Then, it creates
  // a map from the unrelaxed derivative structure to the PARTCAR.
  // The maps are combined in the end.
  vector<uint> POccCalculator::getMapToPARTCAR(unsigned long long int isupercell, const xstructure& xstr_in) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG); 
    string message = "";
    if (isupercell > l_supercell_sets.size() - 1) {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, "Invalid supercell index.");
    }

    std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();
    std::advance(it, isupercell);
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Input structure:" << std::endl;
      std::cerr << xstr_in << std::endl;
    }

    // Get the derivatice structure from the unique structure file
    vector<string> vlines;
    string structures_file = aurostd::CleanFileName(m_aflags.Directory + "/" + POCC_FILE_PREFIX + POCC_UNIQUE_SUPERCELLS_FILE);
    if (!aurostd::EFileExist(structures_file)) {
      message = "Could not load derivative structure."
        " Unique structures file " + structures_file + " not found.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    aurostd::efile2vectorstring(structures_file, vlines);
    stringstream poscar;
    bool found_POSCAR_START_tag = false, found_POSCAR_STOP_tag = false;
    unsigned long long int structure_counter = 0;
    for (uint i = 0; i < vlines.size(); i++) {
      if (!found_POSCAR_START_tag) {
        found_POSCAR_START_tag = aurostd::substring2bool(vlines[i], POSCAR_START_tag);
        if (found_POSCAR_START_tag) continue;
      }
      found_POSCAR_STOP_tag = aurostd::substring2bool(vlines[i], POSCAR_STOP_tag);
      if (found_POSCAR_STOP_tag) {
        if (!found_POSCAR_START_tag) {
          message = "Could not load derivative structure."
            "Found POSCAR stop tag before start tag while reading " + structures_file + ".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
        }
        found_POSCAR_START_tag = false;
        if (structure_counter == isupercell) break;
        else structure_counter++;
      }
      if (found_POSCAR_START_tag && (structure_counter == isupercell)) {
        poscar << vlines[i] << std::endl;
      }
    }
    xstructure xstr_derivative(poscar);
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Derivative structure:" << std::endl;
      std::cerr << xstr_derivative << std::endl;
    }

    // Match structures
    bool same_species = true;
    bool scale_volume = true;
    bool optimize_match = false;
    double misfit = 0.0;
    structure_mapping_info mapping_info = compare::initialize_misfit_struct();
    bool matched = false;
    uint natoms_in = xstr_in.atoms.size();
    uint natoms_derivative = xstr_derivative.atoms.size();
    if (natoms_in < natoms_derivative) {
      matched = compare::aflowCompareStructure(xstr_derivative, xstr_in, same_species, scale_volume, optimize_match, misfit, mapping_info);
    } else {
      matched = compare::aflowCompareStructure(xstr_in, xstr_derivative, same_species, scale_volume, optimize_match, misfit, mapping_info);
    }
    if (!matched) {
      message = "Could not map structure to derivative structure.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << "atom_map = " << aurostd::joinWDelimiter(mapping_info.atom_map, " ") << std::endl;

    // Create map from the relaxed structrue to the derivative structure
    // The number of atoms in the relaxed structure may be bigger larger
    // than in the derivative structure (sconv vs. sprim). In that case,
    // structure comparison only returns the maps to the supercell of the
    // derivative structure with the same number of atoms.
    // Luckily, structure comparison creates supercells in a predictable
    // way (0 0 1 1 2 2 instead of,  e.g. 0 1 2 0 1 2) so we can use the
    // ratio of atoms to map to the actual derivative structure.
    uint atom_ratio = 0;
    vector<uint> map_to_derivative(natoms_in);
    if (natoms_in < natoms_derivative) {
      // NOT TESTED!
      // If the input structure has fewer atoms, then we need to do a
      // one-to-many map. It should not matter which equivalent atoms
      // in the derivative structure it maps to.
      atom_ratio = natoms_derivative/natoms_in;
      int index = -1;
      for (uint i = 0; i < natoms_in; i++) {
        aurostd::WithinList(mapping_info.atom_map, i * atom_ratio, index);
        map_to_derivative[i] = index;
      }
    } else {
      atom_ratio = natoms_in/natoms_derivative;
      for (uint i = 0; i < mapping_info.atom_map.size(); i++) {
        map_to_derivative[i] = mapping_info.atom_map[i]/atom_ratio;
      }
    }

    if (LDEBUG) std::cerr << "map_to_derivative = " << aurostd::joinWDelimiter(map_to_derivative, " ") << std::endl;

    // To get the map from the derivative structure to the PARTCAR, get the
    // full site configuration vector from the derivative structure title
    vector<string> tokens, tokens2;
    string config_str = "";
    aurostd::string2tokens(xstr_derivative.title, tokens, " ");
    for (uint i = 0; tokens.size(); i++) {
      if (tokens[i].find("site_config(") != string::npos) {
        aurostd::string2tokens(tokens[i], tokens2, "/");
        if (tokens2.size() != 2) continue;  // not the right site_config string
        aurostd::StringSubst(tokens2[0], "site_config(", "");
        aurostd::StringSubst(tokens2[0], "#", "");
        if (!aurostd::isfloat(tokens2[0])) continue;  // not the right site_config string
        aurostd::string2tokens(tokens[i], tokens2, "=");
        config_str += tokens2[1];
        i++;
        do {
          config_str += " " + tokens[i];
          if (i == tokens.size() - 1) break;
        } while(tokens[i++].find("]") == string::npos);
        break;
      }
    }
    if (config_str.empty()) {
      message = "Could not determine full site configuration.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::string2tokens(config_str.substr(1, config_str.length() - 2), tokens, ";");
    vector<vector<int> > config(tokens.size());
    for (uint i = 0; i < tokens.size(); i++) {
      aurostd::string2tokens(tokens[i], tokens2, " ");
      for (uint j = 0; j < tokens2.size(); j++) {
        config[i].push_back(aurostd::string2utype<int>(tokens2[j]));
      }
    }
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Site configuration:" << std::endl;
      for (uint i = 0; i < config.size(); i++) {
        std::cerr << aurostd::joinWDelimiter(config[i], " ") << std::endl;
      }
    }

    // Site configurations are stored by type, so we need to create a
    // type to PARTCAR map for each site
    uint nsites = m_pocc_sites.size();
    uint ntypes = xstr_pocc.species.size();
    vector<vector<int> > type2occMap(nsites, vector<int>(ntypes, -1)); // -1 means no atoms of this type on that site
    int type = 0, occ = 0;
    for (uint site = 0; site < m_pocc_sites.size(); site++) {
      for (uint t = 0; t < m_pocc_sites[site].v_types.size(); t++) {
        type = (int) m_pocc_sites[site].v_types[t];
        occ = (int) m_pocc_sites[site].v_occupants[t];
        type2occMap[site][type] = occ;
      }
    }
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " type2occMap:" << std::endl;
      for (uint i = 0; i < type2occMap.size(); i++) {
        std::cerr << aurostd::joinWDelimiter(type2occMap[i], " ") << std::endl;
      }
    }

    // Create a flat map from the derivative structure to PARTCAR
    vector<uint> derivative_partcar_map;
    for (uint site = 0; site < config.size(); site++) {
      std::sort(config[site].begin(), config[site].end());  // Site configurations are not sorted, but structures always are
      for (uint at = 0; at < config[site].size(); at++) {
        type = config[site][at];
        if (type < 0) continue;  // type can be -1 when sites are not used
        occ = type2occMap[site][type];
        if (occ < 0) {
          message = "Mismatch between derivative structure and PARTCAR.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        } else {
          derivative_partcar_map.push_back((uint) occ);
        }
      }
    }

    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " derivative_partcar_map = " << aurostd::joinWDelimiter(derivative_partcar_map, " ") << std::endl;

    // Finally, combine the maps to map the input structure to PARTCAR
    vector<uint> map_to_partcar;
    for (uint i = 0; i < map_to_derivative.size(); i++) {
      map_to_partcar.push_back(derivative_partcar_map[map_to_derivative[i]]);
    }
    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " map_to_partcar = " << aurostd::joinWDelimiter(map_to_partcar, " ") << std::endl;

    return map_to_partcar;
  }

} // namespace pocc

namespace pocc {
  //--------------------------------------------------------------------------------
  // class POccSiteConfiguration (nested in POccCalculator)
  //--------------------------------------------------------------------------------
  POccSiteConfiguration::POccSiteConfiguration() {
    free();
  }

  POccSiteConfiguration::POccSiteConfiguration(int _site,int _i_hnf,vector<POccUnit>& pocc_groups) {
    free();
    site=_site;
    i_hnf=_i_hnf;
    m_pocc_groups=pocc_groups;
  }

  POccSiteConfiguration::POccSiteConfiguration(const POccSiteConfiguration& b){copy(b);}

  POccSiteConfiguration::~POccSiteConfiguration() {free();}

  void POccSiteConfiguration::free() {
    site=0;
    i_hnf=0;
    m_pocc_groups.clear();
    //cerr << "WOO1" << endl;
    //delete(v_pocc_groups);
    //cerr << "WOO2" << endl;
    //vector<POccUnit> tmp; v_pocc_groups=tmp;
    //delete v_pocc_groups;
    partial_occupation_flag=false;
    xv_occupation_count_input.clear();
    xv_occupation_multiple.clear();
    xv_occupation_count_supercell.clear();
    xv_partial_occupation_value.clear();
    xv_site_error.clear();
    occupation_count_total=0;
    vacancy_count=0;
    max_site_error=0.0;
    //error_total=0.0;
    //types_configuration.clear();
  }

  void POccSiteConfiguration::copy(const POccSiteConfiguration& b) {
    site=b.site;
    i_hnf=b.i_hnf;
    m_pocc_groups=b.m_pocc_groups;
    partial_occupation_flag=b.partial_occupation_flag;
    xv_occupation_count_input=b.xv_occupation_count_input;
    xv_occupation_multiple=b.xv_occupation_multiple;
    xv_occupation_count_supercell=b.xv_occupation_count_supercell;
    xv_partial_occupation_value=b.xv_partial_occupation_value;
    xv_site_error=b.xv_site_error;
    occupation_count_total=b.occupation_count_total;
    vacancy_count=b.vacancy_count;
    //error_total=b.error_total;
    max_site_error=b.max_site_error;
    //types_configuration_debug=b.types_configuration_debug;
  }

  const POccSiteConfiguration& POccSiteConfiguration::operator=(const POccSiteConfiguration& b) { // operator= PUBLIC
    if(this!=&b) {free();copy(b);}
    return *this;
  }

  void POccSiteConfiguration::clear() {POccSiteConfiguration a;copy(a);} // clear PRIVATE

  void POccSiteConfiguration::prepareNoPOccConfig() {
    //////////////////////////////////////////////////////////////////////////////
    //clear all fields, but save that which is declared in constructor
    int _site=site, _i_hnf=i_hnf;
    vector<POccUnit> pocc_groups=m_pocc_groups;
    clear();
    site=_site; i_hnf=_i_hnf;
    m_pocc_groups=pocc_groups;
    //////////////////////////////////////////////////////////////////////////////
    partial_occupation_flag = false;
    xvector<int> xvi_dummy(0,0); 
    xvector<double> xvd_dummy(0,0); 
    xvi_dummy(0)=1; xv_occupation_count_input=xvi_dummy;
    xvi_dummy(0)=i_hnf; xv_occupation_multiple=xvi_dummy; xv_occupation_count_supercell=xvi_dummy;
    occupation_count_total=i_hnf;
    xv_site_error=xvd_dummy;
    xv_partial_occupation_value=xvd_dummy;
  }

  void POccSiteConfiguration::preparePOccConfig() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    //////////////////////////////////////////////////////////////////////////////
    //clear all fields, but save that which is declared in constructor
    int _site=site, _i_hnf=i_hnf;
    vector<POccUnit> pocc_groups=m_pocc_groups;
    clear();
    site=_site; i_hnf=_i_hnf;
    m_pocc_groups=pocc_groups;
    //////////////////////////////////////////////////////////////////////////////
    partial_occupation_flag = true;
    int pocc_groups_count=m_pocc_groups.size();
    xvector<int> xvi_dummy(pocc_groups_count-1,0);
    xvector<double> xvd_dummy(pocc_groups_count-1,0);
    //cerr << pocc_groups_count << endl;
    //cerr << xvd_dummy.rows << endl;
    xv_occupation_count_input=xvi_dummy;
    xv_occupation_multiple=xvi_dummy;
    xv_occupation_count_supercell=xvi_dummy;
    xv_partial_occupation_value=xvd_dummy;
    xv_site_error=xvd_dummy;
    //create xv of counts of _pocc_group occupants
    //for(int i=0;i<xv_occupation_count_input.rows;i++){xv_occupation_count_input(i)=(*m_pocc_groups)[i].v_occupants.size();cerr << "xv_occupation_count_input(" << i << ")=" << xv_occupation_count_input(i) << endl;}
    for(int i=0;i<xv_occupation_count_input.rows;i++){xv_occupation_count_input[i]=m_pocc_groups[i].v_occupants.size();}
    //debug
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " site=" << site << ", i_hnf=" << i_hnf << endl;
      for(uint i=0;i<(uint)pocc_groups_count;i++){
        cerr << __AFLOW_FUNC__ << " pocc_group=" << i+1 << "/" << m_pocc_groups.size() << ": ";
        cerr << "v_occupants.size()=" << m_pocc_groups[i].v_occupants.size() << " ?= ";
        cerr << "xv_occupation_count_input=" << xv_occupation_count_input[i] << " ";
        cerr << endl;
      }
    }
    //debug
  }

  int POccSiteConfiguration::getNextOccupationMultiple(int i_hnf, xvector<int>& xv_next_occupation_multiple) {
    //this is a "dummy" bit-string constructor, except canonical bit-string constructors rely on recursion to enumerate ALL possibilities
    //in general, this is not memory-safe
    //so we enumerator ONE-BY-ONE, with a stop condition in place
    //also, most bit-string constructors start with the last digit, and iterate to the first digit (this is how binary numbers are read)
    //we do the opposite, iterate first digit first
    //it conveniently returns index (starts at 0) of _pocc_group impacted
    //i_hnf=_i_hnf;
    //cerr << "HERE  " << xv_next_occupation_multiple << endl;
    for(int i=0;i<xv_next_occupation_multiple.rows;i++){
      if(xv_next_occupation_multiple(i)>=i_hnf){  //greater than not necessary, but safe
        if(i==xv_next_occupation_multiple.rows-1){break;}  //STOP condition, safe
        xv_next_occupation_multiple(i)=0; continue;
      }
      xv_next_occupation_multiple(i)++;
      //cerr << xv_next_occupation_multiple << endl;
      return i;
    }
    return -1;  //STOP condition, safe
  }

  int POccSiteConfiguration::calculateOccupationCountTotal(xvector<int>& xv_next_occupation_multiple){
    return aurostd::scalar_product(xv_next_occupation_multiple,xv_occupation_count_input);
  }

  void POccSiteConfiguration::updateOccupationCounts(int _i_hnf, xvector<int> & xv_next_occupation_multiple){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    i_hnf=_i_hnf;
    xv_occupation_multiple=xv_next_occupation_multiple;
    occupation_count_total=calculateOccupationCountTotal(xv_occupation_multiple); //aurostd::scalar_product(xv_occupation_multiple,xv_occupation_count_input);
    vacancy_count=i_hnf-occupation_count_total;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xv_occupation_multiple=" << xv_occupation_multiple << endl;
      cerr << __AFLOW_FUNC__ << " xv_occupation_count_input" << xv_occupation_count_input << endl;
      cerr << __AFLOW_FUNC__ << " i_hnf=" << i_hnf << endl;
      cerr << __AFLOW_FUNC__ << " occupation_count_total=" << occupation_count_total << endl;
      cerr << __AFLOW_FUNC__ << " vacancy_count=" << vacancy_count << endl;
    }
  }

  void POccSiteConfiguration::calculateError() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    //now check max_site_error
    max_site_error=0.0;
    for(uint pocc_group=0;pocc_group<m_pocc_groups.size();pocc_group++){
      xv_occupation_count_supercell[pocc_group]=xv_occupation_multiple[pocc_group]*xv_occupation_count_input[pocc_group];
      xv_partial_occupation_value[pocc_group] = (xv_occupation_multiple[pocc_group]/(double)i_hnf); //0.0;
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " m_pocc_groups[" << pocc_group << "].partial_occupation_value=" << m_pocc_groups[pocc_group].partial_occupation_value << endl;
        cerr << __AFLOW_FUNC__ << " xv_occupation_count_supercell[" << pocc_group << "]=" << xv_occupation_count_supercell[pocc_group] << endl;
        cerr << __AFLOW_FUNC__ << " xv_partial_occupation_value[" << pocc_group << "]=" << xv_partial_occupation_value[pocc_group] << endl;
        cerr << __AFLOW_FUNC__ << " occupation_count_total=" << occupation_count_total << endl;
        cerr << __AFLOW_FUNC__ << " xv_occupation_count_input[" << pocc_group << "]=" << xv_occupation_count_input[pocc_group] << endl;
      }
      //err = (double)xv_occupation_count_input[pocc_group] * m_pocc_sites[site].m_pocc_groups[pocc_group].partial_occupation_value;
      xv_site_error[pocc_group] = m_pocc_groups[pocc_group].partial_occupation_value;
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " xv_site_error[" << pocc_group << "] = m_pocc_groups[" << pocc_group << "].partial_occupation_value = " << m_pocc_groups[pocc_group].partial_occupation_value << endl;
      }
      //if(occupation_count_total > 0){
      //  xv_partial_occupation_value[pocc_group] = ((double)xv_occupation_count_supercell[pocc_group] / (double)i_hnf); //occupation_count_total);
      //  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xv_partial_occupation_value[pocc_group]  = ((double)xv_occupation_count_supercell[pocc_group] / (double)i_hnf) = " << xv_partial_occupation_value[pocc_group] << endl;}
      //  xv_partial_occupation_value[pocc_group] /= (double)xv_occupation_count_input[pocc_group];  //NORMALIZE to 1.0
      //  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xv_partial_occupation_value[pocc_group] /= (double)xv_occupation_count_input[pocc_group] = " << xv_partial_occupation_value[pocc_group] << endl;}
      xv_site_error[pocc_group] -= xv_partial_occupation_value[pocc_group];
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xv_site_error[" << pocc_group << "] -= xv_partial_occupation_value[" << pocc_group << "] = " << xv_site_error[pocc_group] << endl;}
      //}
      //err/=(double)xv_occupation_count_input[pocc_group];  //NORMALIZE to 1.0
      //xv_site_error[pocc_group] = abs(xv_site_error[pocc_group]); //stability
      //error_total += xv_site_error[pocc_group];
      //cerr << "GROUP " << i << endl;
      ////cerr << "POCC " << m_pocc_sites[site].m_pocc_groups[pocc_group].partial_occupation_value << endl;
      //cerr << "MULTIPLE " << xv_occupation_multiple[pocc_group] << endl;
      //cerr << "COUNT " << xv_occupation_count_input[pocc_group] << endl;
      //cerr << "TOTAL " << xv_occupation_count_supercell[pocc_group] << endl;
      //cerr << "ERR1 " << xv_site_error[pocc_group] << endl;
      ////cerr << "ERR " << err << endl;
      //cerr << endl;
    }
    xv_site_error=aurostd::abs(xv_site_error);  //could be negative
    max_site_error = aurostd::max(xv_site_error);
    //cerr << "ERROR " << max_site_error << endl;
    //cerr << endl << endl;
  }

  //double POccSiteConfiguration::getErrorTotal() const {
  //  return error_total;
  //}

  bool POccSiteConfiguration::isPartiallyOccupied() const {
    return partial_occupation_flag;
  }

  vector<int> POccSiteConfiguration::getStartingTypesConfiguration() const {
    //return types_configuration_debug; //debug
    vector<int> starting_config(occupation_count_total+vacancy_count);
    //cerr << "WOOOOT " << starting_config.size() << endl;
    //cerr << "occupation_count_total=" << occupation_count_total << endl;
    //cerr << "vacancy_count=" << vacancy_count << endl;
    int site_count=0;
    for(int i=0;i<vacancy_count;i++){starting_config[i]=-1;} site_count+=vacancy_count;
    //for(int i=0;i<vacancy_count;i++){starting_config[site_count]=-1; site_count++;}
    for(uint group=0;group<m_pocc_groups.size();group++){
      for(uint occ=0;occ<m_pocc_groups[group].v_types.size();occ++){
        for(int i=0;i<xv_occupation_multiple[group];i++){
          starting_config[site_count]=m_pocc_groups[group].v_types[occ]; site_count++;
        }
      }
    }
    std::sort(starting_config.begin(),starting_config.end());
    //cerr << "site " << site << ": "; 
    //for(uint i=0;i<starting_config.size();i++){
    //  cerr << starting_config[i] << " ";
    //}
    //cerr << endl;
    //types_configuration=starting_configuration;
    return starting_config;
  }

} // namespace pocc

namespace pocc {
  vector<uint> getVacanciesSuperCell(const vector<int>& pc2sc_map,const vector<vector<int> >& v_types_config) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    vector<uint> v_vacancies;

    int type;
    uint starting_supercell_atom_index,supercell_atom_index;
    for(uint site=0;site<v_types_config.size();site++){
      starting_supercell_atom_index=pc2sc_map[site];
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " starting_supercell_atom_index=" << starting_supercell_atom_index << endl;}
      for(uint i=0;i<v_types_config[site].size();i++){
        type=v_types_config[site][i];
        if(type<0){
          supercell_atom_index=starting_supercell_atom_index+i;
          if(LDEBUG) {cerr << __AFLOW_FUNC__ << " supercell_atom_index=" << supercell_atom_index << endl;}
          //get index of corresponding site in supercell
          //then simply iterate to the particular site
          //remember, v_types_config[site][type]
          v_vacancies.push_back(supercell_atom_index);
        }
      }
    }

    std::sort(v_vacancies.rbegin(),v_vacancies.rend()); //sort greater first
    return v_vacancies;
  }

  void replaceRandomSitesSuperCell(
      const xstructure& xstr_pocc,
      const vector<uint>& types2pc_map,
      const vector<int>& pc2sc_map,
      const vector<vector<int> >& v_types_config,
      xstructure& supercell){
#if _DEBUG_POCC_CLUSTER_ANALYSIS_
    stringstream message;
#endif
    uint starting_supercell_atom_index,pocc_atom_index,supercell_atom_index;
    int type;
    xvector<double> fpos,cpos;

#if _DEBUG_POCC_CLUSTER_ANALYSIS_
    if(types2pc_map.size()==0){
      message << "types2pc_map.size()==0";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }

    if(pc2sc_map.size()==0){
      message << "pc2sc_map.size()==0";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }
#endif

    for(uint site=0;site<v_types_config.size();site++){
      starting_supercell_atom_index=pc2sc_map[site];
#if _DEBUG_POCC_CLUSTER_ANALYSIS_
      cerr << __AFLOW_FUNC__ << " starting_supercell_atom_index=" << starting_supercell_atom_index << endl;
#endif
      for(uint i=0;i<v_types_config[site].size();i++){
        type=v_types_config[site][i];
        //switch these atoms
        //supercell.atoms[starting_supercell_atom_index+i]
        //xstr_pocc.atoms[types2pc_map[(*v_types_config)[site][i]]]
        //pocc_atom_index=types2pc_map[(*v_types_config)[site][i]];
#if _DEBUG_POCC_CLUSTER_ANALYSIS_
        cerr << __AFLOW_FUNC__ << " type=" << type << endl;
#endif
        if(type>=0){  //not a vacancy
          //if(replace_types){
          supercell_atom_index=starting_supercell_atom_index+i;
          pocc_atom_index=types2pc_map[type]; //p_str->types2pcMap(type);
          const _atom& pocc_atom = xstr_pocc.atoms[pocc_atom_index]; //xstr_pocc.atoms[pocc_atom_index]; //p_str->xstr_pocc.atoms[pocc_atom_index];
          _atom& supercell_atom = supercell.atoms[supercell_atom_index];
          //cerr << pocc_atom << endl;
          //cerr << supercell_atom_index << endl;
          //cerr << supercell.atoms.size() << endl;

          //updates everything of pocc_atom except positions (fpos,cpos)
          fpos=supercell_atom.fpos; cpos=supercell_atom.cpos;
          supercell_atom=pocc_atom;
          supercell_atom.fpos=fpos; supercell_atom.cpos=cpos;

          //update right type
          supercell_atom.type=type;
          supercell_atom.partial_occupation_flag=false;
          supercell_atom.partial_occupation_value=1.0;
          //supercell.num_each_type[type]++;
          //supercell.comp_each_type[type]+=supercell_atom.partial_occupation_value;
          //cerr << "supercell.atoms[" << supercell_atom_index << "].type=" << supercell_atom.type << endl;
          //}
        }
      }
    }
  }

  void rebuildSuperCell(const xstructure& xstr_pocc,const vector<uint>& v_vacancies,xstructure& supercell) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    //remove vacancies from atoms first
    //removeVacancies(supercell,v_vacancies);
    //we have the right atoms, we just need to rearrange and get right counts
    //std::sort(supercell.atoms.begin(),supercell.atoms.end(),sortAtomsDist); //pocc::sortAtoms); //safe because we do Addatom() below
    //if(LDEBUG) {cerr << __AFLOW_FUNC__ << " new atom count=" << supercell.atoms.size() << endl;}

    deque<_atom> atoms=supercell.atoms;

    //remove vacancies
    if(v_vacancies.size()){
      //std::sort(v_vacancies.rbegin(),v_vacancies.rend()); //USE getVacanciesSuperCell(), does automatically!
      for(uint i=0;i<v_vacancies.size();i++){atoms.erase(atoms.begin()+v_vacancies[i]);}
    }

    std::stable_sort(atoms.begin(),atoms.end(),sortAtomsDist); //pocc::sortAtoms);  //safe because we do Addatom() below
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " new atom count=" << atoms.size() << endl;}

    //I have two options here, I can:
    //1) make a full copy of atoms and use AddAtom(), followed by copying over all pp info, or 
    //2) I can simply copy all the stuff over from xstr_pocc

    //get vector of applicable types
    vector<uint> new_types; //vacancies are -1, and they have been removed
    vector<uint> new_types_atoms;
    vector<uint> new_species;
    bool found;
    for(uint atom=0;atom<atoms.size();atom++){
      found=false;
      for(uint type=0;type<new_types.size()&&!found;type++){
        if(atoms[atom].type==(int)new_types[type]){found=true;}
      }
      if(!found){
        new_types.push_back(atoms[atom].type);
        new_types_atoms.push_back(atom);
      }
    }

    //make sure they are in order
    for(uint type=1;type<new_types.size();type++){
      if(new_types[type-1]>new_types[type]){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Error in atom types reassignment, atoms are out of order");}
    }

    //load in to new clean structure
    xstructure supercell_clean;
    supercell_clean.lattice=supercell.lattice;
    supercell_clean.scale=supercell.scale;
    supercell_clean.FixLattices();

    for(uint atom=0;atom<atoms.size();atom++){supercell_clean.AddAtom(atoms[atom],false);}  //CO20230319 - add by type

    //we need to copy the follow species info over, which is not copied with AddAtom(), space is created for them though
    //species_pp_type
    //species_pp_version
    //species_pp_ZVAL
    //species_pp_vLDAU

    uint atom,type;
    for(uint it=0;it<new_types.size();it++){
      type=new_types[it];
      atom=new_types_atoms[it];
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " type[" << it << "]=" << type << endl;
        cerr << __AFLOW_FUNC__ << " atom[" << it << "]=" << atom << endl;
      }
      //the if statement ensures we don't seg fault
      if((uint)atoms[atom].type<xstr_pocc.species_pp_type.size()){supercell_clean.species_pp_type[it]=xstr_pocc.species_pp_type[type];}  //species_pp_type
      if((uint)atoms[atom].type<xstr_pocc.species_pp_version.size()){supercell_clean.species_pp_version[it]=xstr_pocc.species_pp_version[type];}  //species_pp_version
      if((uint)atoms[atom].type<xstr_pocc.species_pp_ZVAL.size()){supercell_clean.species_pp_ZVAL[it]=xstr_pocc.species_pp_ZVAL[type];}  //species_pp_ZVAL
      if((uint)atoms[atom].type<xstr_pocc.species_pp_vLDAU.size()){supercell_clean.species_pp_vLDAU[it]=xstr_pocc.species_pp_vLDAU[type];}  //species_pp_vLDAU
    }

    //if(0){  //this works, but it's not stable, what if we add more species properties to xstructure? better to start fresh and load in what you need!
    //supercell.ClearSpecies();
    //for(uint type=0;type<new_types.size();type++){
    //supercell.num_each_type.push_back(0);
    //supercell.comp_each_type.push_back(0.0);
    //}
    //
    //uint last_type=0;
    //for(uint atom=0;atom<supercell.atoms.size();atom++){
    //found=false;
    //for(uint type=last_type;type<new_types.size()&&!found;type++){
    //if(supercell.atoms[atom].type==(int)new_types[type]){
    //supercell.num_each_type[type]++;
    //supercell.comp_each_type[type]+=supercell.atoms[atom].partial_occupation_value;
    //last_type=type;
    //found=true;
    //}
    //}
    //}
    //
    ////clear what you can with ClearSpecies()
    ////add here what's needed for species specification in POSCAR
    //for(uint type=0;type<new_types.size();type++){
    //if(type>xstr_pocc.species.size()-1){
    //message << "Bad index with supercell.species[i=" << type << ">" << xstr_pocc.species.size()-1 << "=species.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species.push_back(xstr_pocc.species[type]);
    //
    //if(type>xstr_pocc.species_pp.size()-1){
    //message << "Bad index with supercell.species_pp[i=" << type << ">" << xstr_pocc.species_pp.size()-1 << "=species_pp.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_pp.push_back(xstr_pocc.species_pp[type]);
    //
    //if(type>xstr_pocc.species_pp_type.size()-1){
    //message << "Bad index with supercell.species_pp_type[i=" << type << ">" << xstr_pocc.species_pp_type.size()-1 << "=species_pp_type.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_pp_type.push_back(xstr_pocc.species_pp_type[type]);
    //
    //if(type>xstr_pocc.species_pp_version.size()-1){
    //message << "Bad index with supercell.species_pp_version[i=" << type << ">" << xstr_pocc.species_pp_version.size()-1 << "=species_pp_version.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_pp_version.push_back(xstr_pocc.species_pp_version[type]);
    //
    //if(type>xstr_pocc.species_pp_ZVAL.size()-1){
    //message << "Bad index with supercell.species_pp_ZVAL[i=" << type << ">" << xstr_pocc.species_pp_ZVAL.size()-1 << "=species_pp_ZVAL.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_pp_ZVAL.push_back(xstr_pocc.species_pp_ZVAL[type]);
    //
    //if(type>xstr_pocc.species_pp_vLDAU.size()-1){
    //message << "Bad index with supercell.species_pp_vLDAU[i=" << type << ">" << xstr_pocc.species_pp_vLDAU.size()-1 << "=species_pp_vLDAU.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_pp_vLDAU.push_back(xstr_pocc.species_pp_vLDAU[type]);
    //
    //if(type>xstr_pocc.species_volume.size()-1){
    //message << "Bad index with supercell.species_volume[i=" << type << ">" << xstr_pocc.species_volume.size()-1 << "=species_volume.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_volume.push_back(xstr_pocc.species_volume[type]);
    //
    //if(type>xstr_pocc.species_mass.size()-1){
    //message << "Bad index with supercell.species_mass[i=" << type << ">" << xstr_pocc.species_mass.size()-1 << "=species_mass.size()-1], likely not properly copied over";
    //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    //}
    //supercell.species_mass.push_back(xstr_pocc.species_mass[type]);
    //}
    //
    ////fix volume and mass if they are missing
    //if(xstr_pocc.species.size()){
    //  //volume
    //  double volume_sum=0.0;
    //  for(uint species=0;species<xstr_pocc.species_volume.size();species++){volume_sum+=xstr_pocc.species_volume[species];}
    //  if(volume_sum<_ZERO_TOL_){  //hard zero
    //    xstr_pocc.species_volume.clear();
    //    for(uint species=0;species<xstr_pocc.species.size();species++){
    //      xstr_pocc.species_volume.push_back(GetAtomVolume(xstr_pocc.species[species]));
    //      if(xstr_pocc.species_volume.back()==NNN){
    //        message << "No volume found for " << xstr_pocc.species[species];
    //        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    //      }
    //    }
    //  }
    //  //mass
    //  double mass_sum=0.0;
    //  for(uint species=0;species<xstr_pocc.species_mass.size();species++){mass_sum+=xstr_pocc.species_mass[species];}
    //  if(mass_sum<_ZERO_TOL_){  //hard zero
    //    xstr_pocc.species_mass.clear();
    //    for(uint species=0;species<xstr_pocc.species.size();species++){
    //      xstr_pocc.species_mass.push_back(GetAtomMass(xstr_pocc.species[species]));
    //      if(xstr_pocc.species_mass.back()==NNN){
    //        message << "No mass found for " << xstr_pocc.species[species];
    //        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    //      }
    //    }
    //  }
    //}
    ////leave for end - fix types
    //last_type=0;
    //for(uint atom=0;atom<supercell.atoms.size();atom++){
    //  found=false;
    //  for(uint type=last_type;type<new_types.size()&&!found;type++){
    //    if(supercell.atoms[atom].type==(int)new_types[type]){
    //      supercell.atoms[atom].type=type;
    //      found=true;
    //    }
    //  }
    //}
    //}

    //if(0){  //this assumes that ALL atom types are created... not true, depends on tol/hnf/vacancies
    //  for(uint i=0;i<types2pc_map.size();i++){
    //    supercell.num_each_type.push_back(0);
    //    supercell.comp_each_type.push_back(0.0);
    //  }
    //
    //  for(uint atom=0;atom<supercell.atoms.size();atom++){
    //    supercell.num_each_type[supercell.atoms[atom].type]++;
    //    supercell.comp_each_type[supercell.atoms[atom].type]+=supercell.atoms[atom].partial_occupation_value;
    //  }
    //
    //  //assume we can just copy from xstr_pocc, which should be true unless you screwed something up royally with input
    //  supercell.species=xstr_pocc.species;
    //  supercell.species_pp=xstr_pocc.species_pp;
    //  supercell.species_pp_type=xstr_pocc.species_pp_type;
    //  supercell.species_pp_version=xstr_pocc.species_pp_version;
    //  supercell.species_pp_ZVAL=xstr_pocc.species_pp_ZVAL;
    //  supercell.species_pp_vLDAU=xstr_pocc.species_pp_vLDAU;
    //  supercell.species_volume=xstr_pocc.species_volume;
    //  supercell.species_mass=xstr_pocc.species_mass;
    //  supercell.order_parameter_atoms=xstr_pocc.order_parameter_atoms;
    //}

    supercell_clean.MakeBasis();
    supercell_clean.MakeTypes();

    //change title
    supercell_clean.title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(supercell.title);
    if(supercell_clean.title.empty()){supercell_clean.buildGenericTitle(true,false);}

    //fix volume
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " fixing volume" << endl;
      cerr << __AFLOW_FUNC__ << " OLD VOLUME=" << supercell_clean.Volume() << endl;
      cerr << __AFLOW_FUNC__ << " supercell_clean.scale=" << supercell_clean.scale << endl;
      cerr << __AFLOW_FUNC__ << " supercell_clean.lattice=" << endl;cerr << supercell_clean.lattice << endl;
    }
    supercell_clean.SetAutoVolume();  //defaults to autovolume if no volume in 
    //supercell_clean.neg_scale=TRUE; //no need, scale is already taken care of (no need to DISPLAY volume)
    supercell_clean.ReScale(1.0); //return back to 1.0 to adjust lattice
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " NEW VOLUME=" << supercell_clean.Volume() << endl;
      cerr << __AFLOW_FUNC__ << " supercell_clean.scale=" << supercell_clean.scale << endl;
      cerr << __AFLOW_FUNC__ << " supercell_clean.lattice=" << endl;cerr << supercell_clean.lattice << endl;
    }

    //full copy over
    supercell=supercell_clean;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " final supercell=" << endl;cerr << supercell << endl;}
  }

} // namespace pocc

namespace pocc {
  //--------------------------------------------------------------------------------
  // class POccCalculatorTemplate
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  POccCalculatorTemplate::POccCalculatorTemplate() {;}
  POccCalculatorTemplate::POccCalculatorTemplate(const POccCalculatorTemplate& b) {copy(b);} // copy PUBLIC
  POccCalculatorTemplate::~POccCalculatorTemplate() {free();}

  void POccCalculatorTemplate::free() {
    xstr_pocc.clear(); //DX20191220 - uppercase to lowercase clear
    m_p_flags.clear();
    stoich_each_type.clear();
    xstr_nopocc.clear(); //DX20191220 - uppercase to lowercase clear
    types2pc_map.clear();
    m_species_redecoration.clear();
    //types2uffparams_map.clear();
  }

  void POccCalculatorTemplate::copy(const POccCalculatorTemplate& b){ // copy PRIVATE
    xstr_pocc=b.xstr_pocc;
    m_p_flags=b.m_p_flags;
    stoich_each_type=b.stoich_each_type;
    xstr_nopocc=b.xstr_nopocc;
    types2pc_map.clear(); for(uint i=0;i<b.types2pc_map.size();i++){types2pc_map.push_back(b.types2pc_map[i]);}
    m_species_redecoration.clear(); for(uint i=0;i<b.m_species_redecoration.size();i++){m_species_redecoration.push_back(b.m_species_redecoration[i]);}
    //types2uffparams_map.clear(); for(uint i=0;i<b.types2uffparams_map.size();i++){types2uffparams_map.push_back(b.types2uffparams_map[i]);}
  }

  void POccCalculatorTemplate::setPOccFlags(const aurostd::xoption& pocc_flags){m_p_flags.clear();m_p_flags=pocc_flags;}
  void POccCalculatorTemplate::setAFlags(const _aflags& aflags) {m_aflags.clear();m_aflags=aflags;}

  void POccCalculatorTemplate::setPOccStructure(const xstructure& _xstr_pocc){
    xstr_pocc=_xstr_pocc;
    if(!pflow::checkAnionSublattice(xstr_pocc)){  //CO20210201
      if(!XHOST.vflag_control.flag("FORCE_POCC")){
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found non-anion in anion sublattice. Please check (and run with --force_pocc).",_VALUE_ILLEGAL_);
      }
    }
    types2pc_map=getTypes2PCMap(xstr_pocc);     //get types2pc_map
    xvector<double> _stoich_each_type(xstr_pocc.stoich_each_type.size()-1,0);stoich_each_type=_stoich_each_type;
    for(uint i=0;i<xstr_pocc.stoich_each_type.size();i++){stoich_each_type[i]=xstr_pocc.stoich_each_type[i];}
  }

  void POccCalculatorTemplate::setNonPOccStructure(const xstructure& _xstr_nopocc){xstr_nopocc=_xstr_nopocc;}

  void POccCalculatorTemplate::setSpeciesRedecoration(const vector<string>& species_redecoration){m_species_redecoration=species_redecoration;}
} // namespace pocc

namespace pocc {
  /*
  //--------------------------------------------------------------------------------
  // class POccStructure
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  POccStructure::POccStructure(ostream& oss) : xStream(),m_initialized(false) {initialize(oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,aflags,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,kflags,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,aflags,kflags,oss);}
  POccStructure::POccStructure(ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,aflags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,kflags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,aflags,kflags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const aurostd::xoption& pocc_flags,ostream& oss) : xStream(),m_initialized(false) {initialize(pocc_flags,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,kflags,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,kflags,oss);}
  POccStructure::POccStructure(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(pocc_flags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,kflags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(xstr_pocc,pocc_flags,aflags,kflags,FileMESSAGE,oss);}
  POccStructure::POccStructure(const POccStructure& b) {copy(b);} // copy PUBLIC

  POccStructure::~POccStructure() {xStream::free();free();}

  const POccStructure& POccStructure::operator=(const POccStructure& b) {  // operator= PUBLIC
  if(this!=&b) {copy(b);}
  return *this;
  }

  void POccStructure::clear() {POccStructure a;copy(a);}  // clear PRIVATE
  void POccStructure::free() {
  POccCalculatorTemplate::free();
  m_initialized=false;
  m_aflags.clear();
  m_kflags.clear();
  xstr_sym.clear(); //DX20191220 - uppercase to lowercase clear
  n_hnf=0;
  for(uint site=0;site<m_pocc_sites.size();site++){m_pocc_sites[site].clear();} m_pocc_sites.clear();
  pocc_atoms_total=0;
  v_str_configs.clear();
  }

  void POccStructure::copy(const POccStructure& b) { // copy PRIVATE
  POccCalculatorTemplate::copy(b);
  xStream::copy(b);
  m_initialized=b.m_initialized;
  m_aflags=b.m_aflags;
  m_kflags=b.m_kflags;
  xstr_sym=b.xstr_sym;
  n_hnf=b.n_hnf;
  m_pocc_sites=b.m_pocc_sites;
  pocc_atoms_total=b.pocc_atoms_total;
  //[OBSOLETE - inside m_pocc_sites ]vv_pocc_groups.clear(); for(uint i=0;i<vv_pocc_groups.size();i++){vv_pocc_groups[i].clear();}
  //[OBSOLETE - inside m_pocc_sites ]for(uint i=0;i<b.vv_pocc_groups.size();i++){vv_pocc_groups.push_back(vector<POccUnit>(0));for(uint j=0;j<b.vv_pocc_groups[i].size();j++){vv_pocc_groups[i].push_back(b.vv_pocc_groups[i][j]);}}
  v_str_configs.clear(); for(uint i=0;i<b.v_str_configs.size();i++){v_str_configs.push_back(b.v_str_configs[i]);}

  //hnf_table_general_precision=b.hnf_table_general_precision;
  //hnf_table_iteration_padding=b.hnf_table_iteration_padding;
  //hnf_table_error_padding=b.hnf_table_error_padding;
  //hnf_table_column_padding=b.hnf_table_column_padding;
  //header_max_stoich_error=b.header_max_stoich_error;
  //header_max_site_error=b.header_max_site_error;
  }

  void POccStructure::initialize(ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,aflags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,kflags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,aflags,kflags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setAFlags(aflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setAFlags(aflags);
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }

  void POccStructure::initialize(const aurostd::xoption& pocc_flags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(pocc_flags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,pocc_flags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,pocc_flags,aflags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,pocc_flags,kflags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(xstr_pocc,pocc_flags,aflags,kflags,*_p_FileMESSAGE,oss);
    f_new_ofstream=true;  //override
  }
  void POccStructure::initialize(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setPOccFlags(pocc_flags);
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setPOccFlags(pocc_flags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setPOccFlags(pocc_flags);
      setAFlags(aflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setPOccFlags(pocc_flags);
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }
  void POccStructure::initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    try{
      setOFStream(FileMESSAGE); f_new_ofstream=false;
      setOSS(oss);
      setPOccFlags(pocc_flags);
      setAFlags(aflags);
      setKFlags(kflags);
      setPOccStructure(xstr_pocc);
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
  }

  void POccStructure::setAFlags(const _aflags& aflags) {m_aflags=aflags;}
  void POccStructure::setKFlags(const _kflags& Kflags) {m_kflags=Kflags;}
  */

    void POccCalculator::generateStructures(const _xvasp& in_xvasp) {
      bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
      stringstream message;
      //[NO NEED, AVASP_populateXVASP() is now const]_aflags aflags=m_aflags; //make copies
      _xvasp xvasp=in_xvasp; //[NO NEED, AVASP_populateXVASP() is now const]_kflags kflags=in_kflags;_vflags vflags=in_vflags; //make copies

      if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccCalculator failed to initialize");}
      writePARTCAR();
      calculateHNF();
      getTotalPermutationsCount();
      calculate();

      message << "Writing out " << POCC_UNIQUE_SUPERCELLS_FILE;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

      stringstream unique_derivative_structures_ss;
      stringstream new_AflowIn_ss;
      bool original_quiet_flag=XHOST.QUIET;
      unsigned long long int isupercell=0;

      //for(unsigned long long int i=0;i<getUniqueSuperCellsCount();i++)
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it)
      { //CO20200106 - patching for auto-indenting
        isupercell=std::distance(l_supercell_sets.begin(),it);
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " starting structure[" << isupercell << "]" << endl;}

        //do first so we only need to run getUniqueSuperCell() once
        xvasp.str.clear(); //DX20191220 - uppercase to lowercase clear
        xvasp.str=getUniqueSuperCell(isupercell);

        //populate POCC_UNIQUE_SUPERCELLS_FILE
        unique_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
        unique_derivative_structures_ss << AFLOW_POCC_TAG << "STRUCTURE " << isupercell+1 << "/" << l_supercell_sets.size() << endl;
        unique_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
        unique_derivative_structures_ss << AFLOW_POCC_TAG << "UFF_ENERGY=" << std::fixed << std::setprecision(15) << (*it).getUFFEnergy() << endl; unique_derivative_structures_ss.unsetf(std::ios_base::floatfield);
        unique_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;
        unique_derivative_structures_ss << POSCAR_POCC_series_START_tag << getARUNString(isupercell) << endl; // ." << ss_pocc_count.str() << endl;
        unique_derivative_structures_ss << xvasp.str; //getUniqueSuperCell(isupercell);
        unique_derivative_structures_ss << POSCAR_POCC_series_STOP_tag << getARUNString(isupercell) << endl; // ." << ss_pocc_count.str() << endl;
        unique_derivative_structures_ss << AFLOWIN_SEPARATION_LINE << endl;

        xvasp.Directory=m_aflags.Directory; //arun_directory;
        xvasp.AVASP_arun=true;
        xvasp.AVASP_arun_mode = "POCC"; //CO20181226
        xvasp.AVASP_arun_runname = getARUNString(isupercell);
        xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX",FALSE); //[CO20191110]high-entropy offers increased solubility over binaries

        XHOST.QUIET=true;

        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " printing xstructure" << endl;
          cerr << xvasp.str << endl;
          cerr << __AFLOW_FUNC__ << " printing POSCAR" << endl;
          cerr << xvasp.POSCAR.str() << endl;
          cerr << __AFLOW_FUNC__ << " xvasp.str.species.size()=" << xvasp.str.species.size() << endl;
        }

        AVASP_populateXVASP(m_aflags,m_kflags,m_vflags,xvasp);
        AVASP_MakeSingleAFLOWIN(xvasp,new_AflowIn_ss,true); //false,-1,false);  //don't write/print and hence don't pthread

        XHOST.QUIET=original_quiet_flag;
      }
      aurostd::stringstream2file(unique_derivative_structures_ss,getOutputPath()+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE);
    }

  // ME20211004
  void POccCalculator::createModuleAflowIns(const _xvasp& xvasp_in, const string& MODULE) {
    string message = "";

    if (m_ARUN_directories.size() == 0) loadDataIntoCalculator();

    // We need relaxed structures to propagate into the modules.
    // If not all directories have a qmvasp file, it is safe to
    // assume that not all have relaxed structures.
    if(!QMVASPsFound()) {
      message = "Not all POCC calculations have finished."
        " Please run AFLOW in all subdirectories before running AFLOW modules."
        " See the README for more information.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    _xvasp xvasp;
    for (uint i = 0; i < m_ARUN_directories.size(); i++) {
      xvasp = xvasp_in;
      try {
        xvasp.str = KBIN::GetMostRelaxedStructure(m_ARUN_directories[i]);
      } catch (aurostd::xerror& e) {
        message = "Could not find valid structure file inside directory " + m_ARUN_directories[i] + ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }
      xvasp.str.CleanStructure();
      xvasp.aopts.push_attached("AFLOWIN_FLAG::MODULE", MODULE);
      xvasp.AVASP_arun=true;
      xvasp.AVASP_arun_mode = "POCC";
      xvasp.AVASP_arun_runname = getARUNString(i);
      xvasp.aopts.flag("FLAG::AVASP_SKIP_NOMIX", false);
      AVASP_populateXVASP(m_aflags, m_kflags, m_vflags, xvasp);
      stringstream aflowin;
      AVASP_MakeSingleAFLOWIN(xvasp, aflowin, true);
    }
  }

  void POccCalculator::cleanPOccStructure() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    xstr_pocc.CleanStructure();  // ME20211004 - Moved to xstructure

    //[OBSOLETE ME20211004]//fix up structure
    //[OBSOLETE ME20211004]xstr_pocc.neg_scale=false;  //NO negative scale... doesn't really matter, scale is one variable
    //[OBSOLETE ME20211004]xstr_pocc.ReScale(1.0);
    //[OBSOLETE ME20211004]xstr_pocc.ShiftOriginToAtom(0);
    //[OBSOLETE ME20211004]xstr_pocc.BringInCell();
    //[OBSOLETE ME20211004]xstr_pocc.clean(); //DX20191220 - uppercase to lowercase clean

    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cleaned PARTCAR" << endl;
      cerr << xstr_pocc << endl;
    }
  }

  vector<uint> getTypes2PCMap(const xstructure& xstr) {
    stringstream message;
    vector<uint> types2pc_map;
    uint atom_index=0;
    for(uint i=0;i<xstr.num_each_type.size();i++){
      types2pc_map.push_back(atom_index);
      atom_index+=xstr.num_each_type[i];
    }

    //quick test of stupidity - START
    //see if types are all unique
    for(uint i=0;i<types2pc_map.size();i++){
      for(uint j=i+1;j<types2pc_map.size();j++){
        if(xstr.atoms[types2pc_map[i]].name==xstr.atoms[types2pc_map[j]].name){
          if(xstr.atoms[types2pc_map[i]].name==""){
            message << "Bad input - please assign names to all inputs. Results vary with different atom types.";
          } else {
            message << "Bad input - atoms " << types2pc_map[i]+1 << "/" << xstr.atoms.size() << " ";
            message << "and " << types2pc_map[j]+1 << "/" << xstr.atoms.size() << " ";
            message << "show to have the same name (" << xstr.atoms[types2pc_map[i]].name << "). ";
            message << "Please reorganize your input by type.";
          }
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
        }
      }
    }
    return types2pc_map;
    //quick test of stupidity - END
  }

  void POccCalculator::preparePOccStructure() {
    cleanPOccStructure(); //cosmetics
    partitionPOccSites(); //get m_pocc_sites
    types2pc_map=getTypes2PCMap(xstr_pocc);     //get types2pc_map
    //START CREATING xstr_nopocc
    //this structure is VERY important, as it is the basis for the symmetry analysis + cluster and bonding analysis
    //we need to STANDARDIZE this structure so UFF energies are chemically-agnostic (as much as possible)
    //this includes relabeling the atoms, volume, and sorting atoms
    //sorting atoms explained below
    xstr_nopocc=createNonPOccStructure(xstr_pocc,m_pocc_sites);  //get parent structure, needs xstr_pocc and m_pocc_sites
    calculateSymNonPOccStructure(false);  //calculate symmetry of parent structure, no verbose
  }

  void POccCalculator::initializePOccStructure() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG); 
    stringstream message;

    if(!xstr_pocc.partial_occupation_flag){
      message << "No partially occupied sites detected in this structure";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    }

    message << "Input PARTCAR" << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    message << xstr_pocc;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);

    //no need to assign fake names to xstr_pocc, we already do in xstr_nopocc

    //clean atom/species names
    pflow::fixEmptyAtomNames(xstr_pocc, false); //true);  //aflow likes pp info //CO come back - use fakenames() too if needed
    //add here fakenames()

    //quick tests of stupidity
    if(xstr_pocc.species.size()==0){
      message << "xstr_pocc.species.size()==0";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INDEX_MISMATCH_);
    }
    if(xstr_pocc.species.size()!=xstr_pocc.num_each_type.size()){
      message << "xstr_pocc.species.size() != xstr_pocc.num_each_type.size()";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INDEX_MISMATCH_);
    }
    if(aurostd::StringsAlphabetic(xstr_pocc.species,false)==false){
      message << "xstr_pocc.species() not alphabetized!";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INDEX_MISMATCH_);
    }

    preparePOccStructure();  //first pass to prepare pocc structure, need symmetry to reshuffle atoms by equivalent sites

    //resort atoms by equivalent atoms
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARTCAR pre-sorted by equivalent atoms" << endl;
      cerr << xstr_pocc << endl;
    }
    propagateEquivalentAtoms2POccStructure();
    xstr_pocc.sortAtomsEquivalent();  //basis safe (AddAtom() inside)!
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " PARTCAR post-sorted by equivalent atoms" << endl;
      cerr << xstr_pocc << endl;
    }
    preparePOccStructure(); //need to redo this chunk after resorting atoms

    message << "Clean PARTCAR" << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    message << xstr_pocc;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);

    message << "Parent structure" << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    message << xstr_nopocc;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);

    redecorateXStructures();  //redecorate to create standard parent structure, eventually we will also standardize xstr_pocc by permutation symmetry

    message << "Standardized (re-decorated) parent structure" << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    message << xstr_nopocc;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);

    ///////////////////////////////////////////////////////////////////
    //
    //xstr_pocc and xstr_nopocc are SET above here, do not modify below
    //
    ///////////////////////////////////////////////////////////////////

    calculateSymNonPOccStructure(true); //calculate standard parent structure symmetry, report results
  }

  const xmatrix<double>& POccCalculator::getLattice() const {return xstr_sym.lattice;}
  const vector<_sym_op>& POccCalculator::getPGroup() const {return xstr_sym.pgroup;}

  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,ostream& oss) {
    _aflags aflags;
    if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");} //CO20190402
    if(aflags.Directory.empty() || aflags.Directory=="./" || aflags.Directory=="."){aflags.Directory=aurostd::getPWD()+"/";} //".";  //CO20180220 //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
    return getPOccSites(xstr_pocc,aflags,oss);
  }
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss) {ofstream FileMESSAGE;return getPOccSites(xstr_pocc,aflags,FileMESSAGE,oss);}
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss) {
    _aflags aflags;
    if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");} //CO20190402
    if(aflags.Directory.empty() || aflags.Directory=="./" || aflags.Directory=="."){aflags.Directory=aurostd::getPWD()+"/";} //".";  //CO20180220 //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
    return getPOccSites(xstr_pocc,aflags,FileMESSAGE,oss);
  }
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    //occupy pocc_sites
    vector<POccUnit> pocc_sites;
    bool new_site=true;
    for(uint atom=0;atom<xstr_pocc.atoms.size();atom++){
      new_site=true;
      const _atom& pocc_atom=xstr_pocc.atoms[atom]; 
      for(uint site=0;site<pocc_sites.size() && new_site;site++){
        POccUnit& pocc_site=pocc_sites[site];
        if(pocc_site.v_occupants.size() && aurostd::modulus(pocc_atom.cpos-xstr_pocc.atoms[pocc_site.v_occupants[0]].cpos)<_AFLOW_POCC_ZERO_TOL_){
          pocc_site.v_occupants.push_back(atom);
          pocc_site.v_types.push_back(pocc_atom.type);
          new_site=false;
        }
      }
      if(new_site){
        pocc_sites.push_back(POccUnit(aflags,FileMESSAGE,oss));
        pocc_sites.back().site=pocc_sites.size()-1;
        pocc_sites.back().v_occupants.push_back(atom);
        pocc_sites.back().v_types.push_back(pocc_atom.type);
        pocc_sites.back().partial_occupation_flag=!aurostd::isequal(pocc_atom.partial_occupation_value,1.0,_AFLOW_POCC_ZERO_TOL_); //if the only other occupant is a vacancy
        pocc_sites.back().m_initialized=true;
      }
    }
    for(uint site=0;site<pocc_sites.size();site++){
      pocc_sites[site].partial_occupation_flag=bool(pocc_sites[site].partial_occupation_flag||pocc_sites[site].v_occupants.size()>1); //if the only other occupant is a vacancy
    } //is it a pocc site?

    if(LDEBUG) { //make some ldebug here
      for(uint site=0;site<pocc_sites.size();site++){
        const POccUnit& pocc_site=pocc_sites[site];
        cerr << __AFLOW_FUNC__ << " pocc_site[site=" << site << "].site=" << pocc_site.site << endl;
        cerr << __AFLOW_FUNC__ << " pocc_site[site=" << site << "].partial_occupation_flag=" << pocc_site.partial_occupation_flag << endl;
        cerr << __AFLOW_FUNC__ << " pocc_site[site=" << site << "].v_occupants=" << aurostd::joinWDelimiter(pocc_site.v_occupants,",") << endl;
        cerr << __AFLOW_FUNC__ << " pocc_site[site=" << site << "].v_types=" << aurostd::joinWDelimiter(pocc_site.v_types,",") << endl;
      }
    }

    //get pocc_groups
    bool new_group=true;
    for(uint site=0;site<pocc_sites.size();site++){
      POccUnit& pocc_site=pocc_sites[site];
      for(uint occ=0;occ<pocc_site.v_occupants.size();occ++){
        new_group=true;
        const _atom& pocc_atom=xstr_pocc.atoms[pocc_site.v_occupants[occ]];
        for(uint group=0;group<pocc_site.m_pocc_groups.size() && new_group;group++){
          POccUnit& pocc_group=pocc_site.m_pocc_groups[group];
          if(aurostd::isequal(pocc_atom.partial_occupation_value,pocc_group.partial_occupation_value,_AFLOW_POCC_ZERO_TOL_)){
            //cerr << pocc_atom.partial_occupation_value << " vs. " << pocc_group.partial_occupation_value << endl;
            pocc_group.v_occupants.push_back(pocc_site.v_occupants[occ]);
            pocc_group.v_types.push_back(pocc_site.v_types[occ]);
            new_group=false;
          }
        }
        if(new_group){
          pocc_site.m_pocc_groups.push_back(POccUnit(aflags,FileMESSAGE,oss));
          pocc_site.m_pocc_groups.back().site=site;
          pocc_site.m_pocc_groups.back().partial_occupation_value=pocc_atom.partial_occupation_value;
          pocc_site.m_pocc_groups.back().v_occupants.push_back(pocc_site.v_occupants[occ]);
          pocc_site.m_pocc_groups.back().v_types.push_back(pocc_site.v_types[occ]);
          pocc_site.m_pocc_groups.back().m_initialized=true;
        }
      }
      //std::sort(pocc_site[site].m_pocc_groups.rbegin(),pocc_site[site].m_pocc_groups.rend(),sortPOccGroup); //NOTE the r, reverse sort, bigger poccs first
    }

    if(LDEBUG) {
      cerr.unsetf(std::ios_base::floatfield); //no fixed
      for(uint site=0;site<pocc_sites.size();site++){
        POccUnit& pocc_site=pocc_sites[site];
        for(uint group=0;group<pocc_site.m_pocc_groups.size();group++){
          POccUnit& pocc_group=pocc_site.m_pocc_groups[group];
          cerr << __AFLOW_FUNC__ << " pocc_group[site=" << site << ",group=" << group << "].";
          cerr << "partial_occupation_value=" << pocc_group.partial_occupation_value << endl;
          cerr << __AFLOW_FUNC__ << " pocc_group[site=" << site << ",group=" << group << "].";
          cerr << "v_occupants=" << aurostd::joinWDelimiter(pocc_group.v_occupants,",") << endl;
          cerr << __AFLOW_FUNC__ << " pocc_group[site=" << site << ",group=" << group << "].";
          cerr << "v_types=" << aurostd::joinWDelimiter(pocc_group.v_types,",") << endl;
        }
      }
    }
    return pocc_sites;
  }

  void POccCalculator::partitionPOccSites() {
    //this function gets m_pocc_sites
    m_pocc_sites=getPOccSites(xstr_pocc);

    //get pocc_atoms_total
    pocc_atoms_total=0;
    for(uint site=0;site<m_pocc_sites.size();site++){
      if(m_pocc_sites[site].partial_occupation_flag){pocc_atoms_total+=m_pocc_sites[site].v_occupants.size();}
    }
  }

  void POccCalculator::getOptimizedSiteCountConfigurations(
      int site,
      int i_hnf,
      vector<POccSiteConfiguration>& v_site_configs) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    //reset
    v_site_configs.clear();

    //test of stupidity
    if(m_pocc_sites[site].m_pocc_groups.size()<1){
      message << "Invalid count of POCC site groups (m_pocc_sites[" << site << "].m_pocc_groups.size()=" << m_pocc_sites[site].m_pocc_groups.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    }

    //trivial site_config, no pocc
    POccSiteConfiguration site_config(site,i_hnf,m_pocc_sites[site].m_pocc_groups); //sc_dummy;

    //if it's truly non-pocc, then there's only one m_pocc_sites[site].m_pocc_groups; 
    //otherwise, there's more than 1 and it doesn't matter
    site_config.prepareNoPOccConfig();

    if(!m_pocc_sites[site].partial_occupation_flag){v_site_configs.push_back(site_config);return;}  //no pocc, trivial

    //otherwise, pocc!
    site_config.preparePOccConfig();
    //int pocc_groups_count=m_pocc_sites[site].v_pocc_groups.size();
    //already initialized to 0, xvector<uint> is ill-defined (in AUROSTD), so just use int, it's okay
    //aurostd::xvector<int> xv_occupation_count_input(pocc_groups_count);
    //aurostd::xvector<int> xv_occupation_multiple(pocc_groups_count), xv_occupation_multiple_last(pocc_groups_count);


    //this is the stop condition
    //for(int i=0;i<site_config.xv_occupation_multiple_last.rows;i++){site_config.xv_occupation_multiple_last(i)=(int)i_hnf;}

    //cerr << "START, n=" << pocc_groups_count << ",k=" << i_hnf << endl;
    int pocc_group_iterator=0;    //current group impacted
    int pocc_group_to_skip=-1;    //slight optimizer, if occupation_number hits i_hnf early, skip remaining possibilities
    bool skipping = false;        //needed for resetting pocc_group_to_skip

    double eps=1.0;
    xvector<int> xv_next_occupation_multiple=site_config.xv_occupation_multiple;
    //int minimum_vacancy_count=i_hnf;
    while(pocc_group_iterator!=-1){
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " looking at xv_next_occupation_multiple=" << xv_next_occupation_multiple << endl;
        cerr << __AFLOW_FUNC__ << " pocc_group_to_skip=" << pocc_group_to_skip << " ?= pocc_group_iterator=" << pocc_group_iterator << endl;
      }
      if(pocc_group_to_skip!=pocc_group_iterator){  //slight optimizer
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " exploring xv_next_occupation_multiple=" << xv_next_occupation_multiple << endl;
        }
        //    cerr << "calculateOccupationCountTotal(xv_next_occupation_multiple)=" << site_config.calculateOccupationCountTotal(xv_next_occupation_multiple) << endl;

        if(site_config.calculateOccupationCountTotal(xv_next_occupation_multiple)>i_hnf){ //start ignoring
          if(LDEBUG) {
            cerr << __AFLOW_FUNC__ << " will skip multiples from now on" << endl;
            cerr << __AFLOW_FUNC__ << " site_config.calculateOccupationCountTotal(xv_next_occupation_multiple)=" << site_config.calculateOccupationCountTotal(xv_next_occupation_multiple) << " > " << i_hnf << "=i_hnf" << endl;
          }
          pocc_group_to_skip=pocc_group_iterator; skipping=true; continue;
        }

        site_config.updateOccupationCounts(i_hnf,xv_next_occupation_multiple);

        //found viable configuration
        site_config.calculateError();

        //    cerr << "eps=" << eps << endl;
        //    cerr << "site_config.max_site_error=" << site_config.max_site_error << endl;
        if(aurostd::isequal(site_config.max_site_error,eps,_AFLOW_POCC_ZERO_TOL_) || site_config.max_site_error<eps){
          if(!aurostd::isequal(site_config.max_site_error,eps,_AFLOW_POCC_ZERO_TOL_)){
            eps=site_config.max_site_error;
            v_site_configs.clear();
          }
          //      cerr << "adding one!!!!!!1" << endl;
          v_site_configs.push_back(site_config);
          //minimum_vacancy_count=min(minimum_vacancy_count,site_config.vacancy_count);
        }
      }
      //getNextOccupationMultiple() is basically a bitstring iterator:  000, 001, 010, ...

      //it will increment xv_occupation_multiple and return the index of the digit that was modified
      //a return value of -1 means it's exhausted all possibilities 
      pocc_group_iterator = site_config.getNextOccupationMultiple(i_hnf,xv_next_occupation_multiple);
      if(skipping && pocc_group_iterator!=pocc_group_to_skip){pocc_group_to_skip=-1; skipping=false;}
      //  cerr << "i_hnf=" << i_hnf << endl;
      //  cerr << "xv_next_occupation_multiple=" << xv_next_occupation_multiple << endl;
      //  cerr << "pocc_group_iterator=" << pocc_group_iterator << endl;
      //  cerr << endl;
    }
    //  for(uint i=0;i<v_site_configs.size();i++){
    //    cerr << "site config " << i << endl;
    //    cerr << "v_site_configs[i].xv_occupation_count_input = " << v_site_configs[i].xv_occupation_count_input << endl;
    //    cerr << "v_site_configs[i].xv_occupation_multiple = " << v_site_configs[i].xv_occupation_multiple << endl;
    //    cerr << "v_site_configs[i].xv_occupation_count_supercell = " << v_site_configs[i].xv_occupation_count_supercell << endl;
    //    cerr << "v_site_configs[i].occupation_count_total = " << v_site_configs[i].occupation_count_total << endl;
    //    cerr << "v_site_configs[i].vacancy_count = " << v_site_configs[i].vacancy_count << endl;
    //    cerr << "v_site_configs[i].max_site_error = " << v_site_configs[i].max_site_error << endl;
    //  }
    //  cerr << "END" << endl;
    //
    //if(REDUCE_VACANCIES){
    //  if(v_site_configs.size()){
    //    for(uint config=v_site_configs.size()-1;config-- > 0; ){  //will stop at 0
    //      if(v_site_configs[config].vacancy_count>minimum_vacancy_count){
    //        v_site_configs.erase(v_site_configs.begin()+config);
    //      }
    //    }
    //  }
    //}

  }

  xvector<double> POccCalculator::calculateStoichEachType(vector<vector<int> >& v_types_config) {
    xvector<double> stoich_each_type(types2pc_map.size()-1,0);
    //for(uint type=0;type<types2pc_map.size();type++){stoich_each_type.push_back(0.0);}
    double total_atoms=0.0;
    for(uint site=0;site<v_types_config.size();site++){
      for(uint atom=0;atom<v_types_config[site].size();atom++){
        if(v_types_config[site][atom]!=-1){stoich_each_type[v_types_config[site][atom]]+=1.0;total_atoms+=1.0;}
      }
    }
    if(total_atoms>0.0){
      for(int type=0;type<stoich_each_type.rows;type++){stoich_each_type[type]/=total_atoms;}
    }
    return stoich_each_type;
  }

  //double POccStructure::calculateStoichDiff(deque<double>& s1,deque<double>& s2){
  //  if(s1.size()!=s2.size()){throw aurostd::xerror(__AFLOW_FILE__,"pocc::POccStructure::calculateStoichDiff():","s1.size()!=s2.size(),_INPUT_ILLEGAL_);}
  //  double error=0.0;
  //  for(uint i=0;i<s1.size();i++){
  //    error+=abs(s1[i]-s2[i]);
  //  }
  //  return error;
  //}

  //void POccStructure::getSiteCountConfigurations(int i_hnf, double& min_stoich_error)
  void POccCalculator::getSiteCountConfigurations(int i_hnf)
  { //CO20200106 - patching for auto-indenting
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    //return these values, do not necessarily store as "global" values

    //reset
    //for(uint str_config=0;str_config<v_str_configs.size();str_config++){
    //  v_str_configs[str_config].clear();
    //}
    v_str_configs.clear();
    //v_site_occupancy.clear();
    //v_site_effective_pocc.clear();
    //v_site_error.clear();
    //v_vacancy_count.clear();

    //only useful here, but super important!
    //look at case i_hnf = 5, and a site with 8 occupants, evenly occupied (easy)
    //the original algorithm (look at each pocc atom individually), 
    //would want to fill this site with one atom (+8 atoms total, not even space!)
    //so we group by occupation as well, and allowed occupations are when all equally occupants can fill together

    //create space
    //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {v_site_error.push_back(1.0);}
    //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {v_site_effective_pocc.push_back(0.0);}
    //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {v_site_occupancy.push_back(0);}
    //for(uint i=0;i<m_pocc_sites.size();i++){v_vacancy_count.push_back(i_hnf);}

    //perform calculation
    //cerr << "HERE4" << endl;
    vector<vector<POccSiteConfiguration> > vv_count_configs;  //first layer site, second layer config (opposite from v_str_configs)
    vector<POccSiteConfiguration> v_site_configs;
    for(uint site=0;site<m_pocc_sites.size();site++){
      //cerr << "site=" << site <<endl;
      //cerr << "i_hnf=" << i_hnf <<endl;
      //cerr << "v_site_configs.size()=" << v_site_configs.size() <<endl;
      getOptimizedSiteCountConfigurations(site,i_hnf,v_site_configs);
      //  cerr << "site = " << site << " count: " << v_site_configs.size() << endl;
      vv_count_configs.push_back(v_site_configs);
      //cerr << "HERE8" << endl;
      //cerr << "site = " << site << endl;
    }
    if(LDEBUG) {
      for(uint site=0;site<vv_count_configs.size();site++){
        cerr << __AFLOW_FUNC__ << " site " << site << endl;
        for(uint config=0;config<vv_count_configs[site].size();config++){
          cerr << __AFLOW_FUNC__ << " config " << config << endl;    //
          cerr << __AFLOW_FUNC__ << " vv_count_configs[" << site << "][" << config << "].xv_site_error=" << vv_count_configs[site][config].xv_site_error << endl;
          cerr << __AFLOW_FUNC__ << " vv_count_configs[" << site << "][" << config << "].max_site_error=" << vv_count_configs[site][config].max_site_error << endl;
          cerr << __AFLOW_FUNC__ << " vv_count_configs[" << site << "][" << config << "].xv_occupation_count_input=" << vv_count_configs[site][config].xv_occupation_count_input << endl;
          cerr << __AFLOW_FUNC__ << " vv_count_configs[" << site << "][" << config << "].xv_occupation_multiple=" << vv_count_configs[site][config].xv_occupation_multiple << endl;
          cerr << __AFLOW_FUNC__ << " vv_count_configs[" << site << "][" << config << "].xv_partial_occupation_value=" << vv_count_configs[site][config].xv_partial_occupation_value << endl;
        }
        //cerr << endl << endl;
      }
    }
    //we need to fundamentally change scheme here to satisfy STRICT_STOICH_EQUIVALENCE
    //we cannot simply permute through each config on each site
    //only certain configs work together, we must define these
    //vector<uint> v_config_order;    //dummy, don't really care about order just yet
    vector<int> v_config_iterators;
    vector<vector<int> > v_types_config;
    xvector<double> stoich_each_type_new;
    double stoich_error=0.0;
    double eps=1.0; //(double)m_pocc_sites.size();   //number of sites
    //cerr << "CO START" << endl;
    //std::vector<double>::iterator it_max_stoich;
    //while(getNextSiteConfiguration(vv_count_configs,v_config_order,v_config_iterators,v_types_config))
    while(pocc::getNextSiteConfiguration(vv_count_configs,v_config_iterators,v_types_config))
    { //CO20200106 - patching for auto-indenting
      stoich_each_type_new=calculateStoichEachType(v_types_config);
      stoich_error=aurostd::max(aurostd::abs(stoich_each_type-stoich_each_type_new));
      //stoich_error=calculateStoichDiff(stoich_each_type,xstr_pocc.stoich_each_type);
      //cerr << "TYPES" << endl;
      //for(uint i=0;i<v_types_config.size();i++){
      //  for(uint j=0;j<v_types_config[i].size();j++){
      //    cerr << v_types_config[i][j] << " ";
      //  }
      //  cerr << endl;
      //}
      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " stoichiometry=" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(stoich_each_type_new,_AFLOW_POCC_PRECISION_),",") << endl;
        //for(int i=0;i<stoich_each_type_new.rows; i++){
        //  cerr << stoich_each_type_new[i] << endl;
        //}
      }
      //cerr << "STOICH_SET" << endl;
      //for(int i=0;i<stoich_each_type.rows; i++){
      //  cerr << stoich_each_type[i] << endl;
      //}
      //cerr << "ERROR " << stoich_error << " eps" << eps  << endl;

      if(aurostd::isequal(stoich_error,eps,_AFLOW_POCC_ZERO_TOL_) || stoich_error<eps){
        //cerr << "FOUND 1" << endl;
        if(!aurostd::isequal(stoich_error,eps,_AFLOW_POCC_ZERO_TOL_)){
          eps=stoich_error;
          //for(uint str_config=0;str_config<v_str_configs.size();str_config++){v_str_configs[str_config].clear();} 
          v_str_configs.clear();
        }
        v_str_configs.push_back(StructureConfiguration());
        v_str_configs.back().max_site_error=0.0;
        for(uint site=0;site<v_config_iterators.size();site++){
          v_str_configs.back().site_configs.push_back(vv_count_configs[site][v_config_iterators[site]]);
          v_str_configs.back().max_site_error=max(v_str_configs.back().max_site_error,v_str_configs.back().site_configs.back().max_site_error);
        }
        //cerr << "SEE " << v_str_configs.back().max_site_error << endl;
        v_str_configs.back().max_stoich_error=stoich_error;
      }
    }
    //cerr << "CO STOP" << endl;

    //vector<int> starting_config;
    //for(uint str_config=0;str_config<v_str_configs.size();str_config++){
    //  //cerr << "config " << str_config << endl;
    //  for(uint site=0;site<v_str_configs[str_config].size();site++){
    //    //cerr << "site " << site << ":  ";
    //    starting_config=v_str_configs[str_config][site].getStartingTypesConfiguration();
    //    //for(uint i=0;i<starting_config.size();i++){
    //    //  cerr << starting_config[i] << " ";
    //    //}
    //    //cerr << endl;
    //  }
    //}
    //cerr << "HERE" << endl;

    //perform calculation
    //double eps;
    //for(uint iatom=0;iatom<xstr_pocc.atoms.size();iatom++) {
    //  v_site_error[iatom]=1.0;
    //  v_site_occupancy[iatom]=0;
    //  for(uint j=0;j<=i_hnf;j++) {
    //    eps=abs((double) j/i_hnf-xstr_pocc.atoms[iatom].partial_occupation_value);
    //    if(eps<v_site_error[iatom]) {
    //      v_site_error[iatom]=eps;
    //      v_site_occupancy[iatom]=j;
    //      v_site_effective_pocc[iatom]=(double) v_site_occupancy[iatom]/i_hnf;
    //    }
    //  }
    //}

    //min_stoich_error=eps;
  }

  xstructure createNonPOccStructure(const xstructure& xstr_pocc,const vector<POccUnit>& pocc_sites,bool use_automatic_volumes_in){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    xstructure xstr_nopocc(xstr_pocc);
    //clear structure pocc properties, set defaults
    xstr_nopocc.scale_second=DEFAULT_POCC_SITE_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
    xstr_nopocc.neg_scale_second=FALSE;
    xstr_nopocc.partial_occupation_flag=FALSE;
    xstr_nopocc.partial_occupation_site_tol=DEFAULT_POCC_SITE_TOL; //DEFAULT_PARTIAL_OCCUPATION_TOLERANCE;
    xstr_nopocc.partial_occupation_HNF=0;
    xstr_nopocc.partial_occupation_sublattice.clear();

    //make pocc sites non-pocc sites by removing all but one occupant
    //gather indices of atoms to remove, then remove backwards
    vector<uint> v_atoms_to_remove;
    for(uint site=0;site<pocc_sites.size();site++){
      if(pocc_sites[site].partial_occupation_flag){
        //test of stupidity
        //[CO COME BACK, add another reasonable check that includes vacancies]if(pocc_sites[site].v_occupants.size()<2){return false;}
        //set the 0th index to non-pocc status
        xstr_nopocc.atoms[pocc_sites[site].v_occupants[0]].partial_occupation_value=1.0;
        xstr_nopocc.atoms[pocc_sites[site].v_occupants[0]].partial_occupation_flag=false;
        //add others to list to be removed
        for(uint atom=1;atom<pocc_sites[site].v_occupants.size();atom++){
          v_atoms_to_remove.push_back(pocc_sites[site].v_occupants[atom]);
        }
      }
    }
    xstr_nopocc.RemoveAtom(v_atoms_to_remove);

    //recalculate volume, assume automatic volumes
    //BIG REVELATION - since POCC now relies on automatic volumes, we need to standardize the volume
    //of cluster for calculating the energy as we did the species so UFF energy is species-agnostic
    //CHf_pvNb_svTa_pvTi_svZr_sv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF
    //this configuration produces 47 configurations, while all other specie combinations (carbides paper) make 49
    //CHANGE SPECIES LATER
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " fixing volume" << endl;
      cerr << __AFLOW_FUNC__ << " OLD VOLUME=" << xstr_nopocc.Volume() << endl;
      cerr << __AFLOW_FUNC__ << " xstr_nopocc.scale=" << xstr_nopocc.scale << endl;
      cerr << __AFLOW_FUNC__ << " xstr_nopocc.lattice=" << endl;cerr << xstr_nopocc.lattice << endl;
      cerr << __AFLOW_FUNC__ << " xstr_nopocc.atoms.size()=" << xstr_nopocc.atoms.size() << endl;
    }
    xstr_nopocc.SetAutoVolume(use_automatic_volumes_in);  //always use automatic_volumes_in
    //xstr_nopocc.neg_scale=TRUE; //no need, scale is already taken care of (no need to DISPLAY volume)
    xstr_nopocc.ReScale(1.0); //return back to 1.0 to adjust lattice
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " NEW VOLUME=" << xstr_nopocc.Volume() << endl;
      cerr << __AFLOW_FUNC__ << " xstr_nopocc.scale=" << xstr_nopocc.scale << endl;
      cerr << __AFLOW_FUNC__ << " xstr_nopocc.lattice=" << endl;cerr << xstr_nopocc.lattice << endl;
    }

    /*
    //CHANGE SPECIES LATER
    //recalculate volume, assume automatic volumes
    //BIG REVELATION - since POCC now relies on automatic volumes, we need to standardize the volume
    //of cluster for calculating the energy as we did the species so UFF energy is species-agnostic
    //CHf_pvNb_svTa_pvTi_svZr_sv:PAW_PBE/AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF
    //this configuration produces 47 configurations, while all other specie combinations (carbides paper) make 49
    if(!(COMPARE_WITH_KESONG||m_use_mixed_set_elements==false)){  //need to fix volume of cell
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " need to fix volume" << endl;}
    vector<string> elements=getElementsList();
    double volume=0,voli=0;
    for(uint i=0;i<xstr_nopocc.atoms.size();i++){
    for(uint j=0;j<xstr_nopocc.num_each_type.size();j++){
    if(xstr_nopocc.atoms[i].name==xstr_nopocc.species[j]){
    //[ORIGINAL]volume+=xstr_nopocc.species_volume[j]*xstr_nopocc.atoms[i].partial_occupation_value;
    voli=GetAtomVolume(elements[j]);
    if(voli==NNN){
    message << "No volume found for " << elements[j];
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
    }
    volume+=voli*xstr_nopocc.atoms[i].partial_occupation_value;
    }
    }
    }
    //xstr_nopocc.scale=std::pow((double) (abs(volume)/det(xstr_nopocc.lattice)),(double) 1.0/3.0);
    xstr_nopocc.SetVolume(volume);  //more robust
    //xstr_nopocc.neg_scale=TRUE; //no need, scale is already taken care of (no need to DISPLAY volume)
    xstr_nopocc.ReScale(1.0); //return back to 1.0 to adjust lattice
    if(LDEBUG) {
    cerr << __AFLOW_FUNC__ << " adjusting xstr_nopocc.volume=" << volume << endl;
    cerr << __AFLOW_FUNC__ << " adjusting xstr_nopocc.scale=" << xstr_nopocc.scale << endl;
    cerr << __AFLOW_FUNC__ << " xstr_nopocc.lattice=" << endl;
    cerr << xstr_nopocc.lattice << endl;
    }
    }
    */

    xstr_nopocc.title="Parent structure";
    xstr_nopocc.clean();    //real true clean //DX20191220 - uppercase to lowercase clean
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " final non-pocc structure" << endl;
      cerr << xstr_nopocc << endl;
    }

    return xstr_nopocc;
  }

  void POccCalculator::redecorateXStructures(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;

    //check that code can handle arbitrary reassignment of random elements
    vector<string> elements_std=getElementsList();
    bool use_standard_elements_set=true;
    if(xstr_pocc.num_each_type.size()>elements_std.size()){
      use_standard_elements_set=false;
      message << "Unexpectedly large count of species in structure (nspecies=" << xstr_pocc.num_each_type.size();
      message <<  ">" << elements_std.size()<< "), resorting to exact UFF energy determination";
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_);
    }
    if(COMPARE_WITH_KESONG){use_standard_elements_set=false;}
    vector<string>& species=m_species_redecoration;

    for(uint i=0;i<types2pc_map.size();i++){species.push_back(xstr_pocc.atoms[types2pc_map[i]].cleanname);}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species             =" << aurostd::joinWDelimiter(species,",") << endl;}

    if(use_standard_elements_set==true){
      species.clear();
      for(uint i=0;i<types2pc_map.size();i++){species.push_back(elements_std[i]);}
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_std         =" << aurostd::joinWDelimiter(species,",") << endl;}

      //RENAMING SPECIES TO ELIMINATE CHEMISTRY-SPECIFIC UFF ENERGIES
      //this code is NOT a general fix: true fix is structure comparison
      //the biggest problem: the parent structure defines the volume for the cluster used to find the UFF energies
      //the volume defines the bond distances and thus which pairs are bonded/non-bonded
      //we need to make sure the parent structure is re-decorated as consistently as possible given different input decorations
      //useful test:
      // Ba_svHf_pvNb_svOSnTi_svZr_sv:PAW_PBE/AB3C_cP5_221_a_c_b.BAC:POCC_S0-1xD_P3-1xA_P4-0.2xB-0.2xC-0.2xE-0.2xF-0.2xG
      // vs.
      // Ba_svHf_pvOSnTi_svY_svZr_sv:PAW_PBE/AB3C_cP5_221_a_c_b.BAC:POCC_S0-1xC_P3-1xA_P4-0.2xB-0.2xD-0.2xE-0.2xF-0.2xG
      bool resort_standard_elements_set=true;  //the approach here is WRONG, we need structure comparison
      if(resort_standard_elements_set==true){
        uint site=0,itype=0,occ=0;
        vector<string> species_orig;  //keep orig species
        for(uint i=0;i<species.size();i++){species_orig.push_back(species[i]);}
        vector<uint> occupant2species_map;  //maps occupants (m_pocc_sites) to species (types) in order of m_pocc_sites, NOT NEEDED
        species.clear();species.resize(species_orig.size());  //make room for species before you fill, need occupant2species_map to dictate where new species go
        vector<POccUnit> m_pocc_sites_resorted;   //make copy of m_pocc_sites for reshuffling
        for(site=0;site<m_pocc_sites.size();site++){
          for(occ=0;occ<m_pocc_sites[site].v_occupants.size();occ++){occupant2species_map.push_back(m_pocc_sites[site].v_types[occ]);}
          m_pocc_sites_resorted.push_back(m_pocc_sites[site]);
        }
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " occupant2species_map=" << aurostd::joinWDelimiter(occupant2species_map,",") << endl;
          cerr << __AFLOW_FUNC__ << " original m_pocc_sites" << endl;
          for(site=0;site<m_pocc_sites_resorted.size();site++){
            cerr << __AFLOW_FUNC__ << " site=" << site << " occupants: ";
            for(occ=0;occ<m_pocc_sites_resorted[site].v_occupants.size();occ++){
              cerr << xstr_pocc.atoms[m_pocc_sites_resorted[site].v_occupants[occ]].name << " ";
            }
            cerr << endl;
          }
        }

        std::stable_sort(m_pocc_sites_resorted.rbegin(),m_pocc_sites_resorted.rend(),sortPOccSites);  //ascending order, fewer occupant-sites get populated first (sublattice anions)

        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " sorted m_pocc_sites" << endl;
          for(site=0;site<m_pocc_sites_resorted.size();site++){
            cerr << __AFLOW_FUNC__ << " site=" << site << " occupants: ";
            for(occ=0;occ<m_pocc_sites_resorted[site].v_occupants.size();occ++){
              cerr << xstr_pocc.atoms[m_pocc_sites_resorted[site].v_occupants[occ]].name << " ";
            }
            cerr << endl;
          }
        }

        uint iocc2spec=0;
        vector<uint> types_added;
        uint type=0;

        //0. prioritize anion sublattice, must find a sublattice full of anions
        vector<string> vanions;
        aurostd::string2tokens(POCC_ANIONS_LIST,vanions,",");
        bool sublattice_found=false;
        for(site=0;site<m_pocc_sites_resorted.size();site++){
          if(m_pocc_sites_resorted[site].v_types.size()!=m_pocc_sites_resorted[site].v_occupants.size()){
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_pocc_sites_resorted[site].v_types.size()!=m_pocc_sites_resorted[site].v_occupants.size()",_RUNTIME_ERROR_);
          }
          sublattice_found=(m_pocc_sites_resorted[site].v_types.size()>0);  //set to true first
          for(itype=0;itype<m_pocc_sites_resorted[site].v_types.size()&&sublattice_found;itype++){
            type=m_pocc_sites_resorted[site].v_types[itype];
            const string& occupant=xstr_pocc.atoms[m_pocc_sites_resorted[site].v_occupants[itype]].name;
            if(!aurostd::WithinList(vanions,occupant)){sublattice_found=false;}
          }
          if(sublattice_found){
            itype=0;  //prioritize first species of sublattice for parent structure
            type=m_pocc_sites_resorted[site].v_types[itype];
            if(!aurostd::WithinList(types_added,type)){
              const string& occupant=xstr_pocc.atoms[m_pocc_sites_resorted[site].v_occupants[itype]].name;
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping " << occupant << " to " << species_orig[iocc2spec] << " [sublattice]" << endl;}
              species[type]=species_orig[iocc2spec++];
              types_added.push_back(type);
            }
          }
        }

        //1. prioritize the first occupant of every site, these occupants will define the parent lattice volume (and thus bond lengths)
        for(site=0;site<m_pocc_sites_resorted.size();site++){
          itype=0;
          type=m_pocc_sites_resorted[site].v_types[itype];
          if(!aurostd::WithinList(types_added,type)){
            const string& occupant=xstr_pocc.atoms[m_pocc_sites_resorted[site].v_occupants[itype]].name;
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping " << occupant << " to " << species_orig[iocc2spec] << " [first-occupant]" << endl;}
            species[type]=species_orig[iocc2spec++];
            types_added.push_back(type);
          }
        }

        //LAST. get remaining species
        for(site=0;site<m_pocc_sites_resorted.size();site++){
          for(itype=0;itype<m_pocc_sites_resorted[site].v_types.size();itype++){
            type=m_pocc_sites_resorted[site].v_types[itype];
            if(!aurostd::WithinList(types_added,type)){
              const string& occupant=xstr_pocc.atoms[m_pocc_sites_resorted[site].v_occupants[itype]].name;
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " mapping " << occupant << " to " << species_orig[iocc2spec] << " [remaining]" << endl;}
              //species[occupant2species_map[iocc2spec++]]=species_orig[type];
              //species[type]=species_orig[occupant2species_map[iocc2spec++]];
              species[type]=species_orig[iocc2spec++];
              types_added.push_back(type);
            }
          }
        }
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species_std_resorted=" << aurostd::joinWDelimiter(species,",") << endl;}
      }

      //redecorate xstr_pocc (tmp) to create xstr_nopocc
      xstructure xstr_pocc_tmp(xstr_pocc);
      xstr_pocc_tmp.SetSpecies(aurostd::vector2deque(species),false); //do NOT change order of species assignment

      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]deque<_atom> atoms;
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]uint iatom=0;
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]for(uint ispec=0;ispec<xstr_pocc_tmp.num_each_type.size();ispec++){ //re-label all atoms
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]  for(uint i=0;i<(uint)xstr_pocc_tmp.num_each_type[ispec];i++){
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " converting " << xstr_pocc_tmp.atoms[iatom].name << " to " << species[ispec] << endl;}
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]    atoms.push_back(xstr_pocc_tmp.atoms[iatom++]);
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]    atoms.back().name=species[ispec];
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]    atoms.back().CleanName();
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]    atoms.back().name_is_given=(!atoms.back().name.empty());
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]  }
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]}
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]std::stable_sort(atoms.begin(),atoms.end(),sortAtomsNames); //safe because we do AddAtom() below
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]//assign types
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]uint itype=0;
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]atoms[0].type=itype;
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]for(uint i=1;i<atoms.size();i++){
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]	if(atoms[i].name!=atoms[i-1].name){itype++;}
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]	atoms[i].type=itype;
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]}
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]for(uint i=xstr_pocc_tmp.atoms.size()-1;i<xstr_pocc_tmp.atoms.size();i--){  //remove atoms
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]  if(LDEBUG) cerr << __AFLOW_FUNC__ << " removing atom[" << i << "]" << endl;
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]  xstr_pocc_tmp.RemoveAtom(i);
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]}
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]for(uint i=0;i<atoms.size();i++){xstr_pocc_tmp.AddAtom(atoms[i],true);}  //CO20230319 - add by species
      //[CO20190317 - do NOT sort atoms (AddAtom()), this will confuse the order that pseudonames are assigned]atoms.clear();

      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " xstr_pocc_tmp" << endl << xstr_pocc_tmp << endl;
        cerr << __AFLOW_FUNC__ << " xstr_pocc_tmp.species=" << aurostd::joinWDelimiter(xstr_pocc_tmp.species,",") << endl;
      }

      vector<POccUnit> pocc_sites_tmp=getPOccSites(xstr_pocc_tmp);
      xstr_nopocc=createNonPOccStructure(xstr_pocc_tmp,pocc_sites_tmp); //replace xstr_nopocc (NOT tmp)
      xstr_nopocc.title="Standardized (re-decorated) parent structure";  //just for printing

      if(LDEBUG) {
        cerr << __AFLOW_FUNC__ << " xstr_pocc pre-decoration" << endl;
        cerr << xstr_pocc << endl;
        cerr << __AFLOW_FUNC__ << " xstr_pocc post-decoration" << endl;
        cerr << xstr_pocc_tmp << endl;
        cerr << __AFLOW_FUNC__ << " xstr_nopocc post-decoration" << endl;
        cerr << xstr_nopocc << endl;
      }
    }

    if(species.size()!=xstr_pocc.species.size()){
      message << "species.size()!=xstr_pocc.species.size()";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_RANGE_);
    }

    ////DO UFF STUFF IF IT'S THE RIGHT MODE
    ////MAKE DIFFERENT MODES HERE
    //string uff_parameters_string;
    //vector<string> uff_params;
    //UFFParamAtom uffp;
    //for(uint i=0;i<species.size();i++){
    //const string& element=species[i];

    //if(LDEBUG) {cerr << __AFLOW_FUNC__ << " element[i=" << i << "]=" << element << endl;}
    //uff_parameters_string=getUFFParameterString(element);
    //if(uff_parameters_string.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unable to fetch UFF parameters (requested "+element+")");}
    //aurostd::string2tokens(uff_parameters_string,uff_params," ");
    //if(uff_params.size()!=12){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unexpected UFF Parameters size",_VALUE_ILLEGAL_);}

    ////where to start tomorrow
    ////we get types2uffparams_map here (poccstructure) but use in uffenergyanalyzer
    ////move this to uffenergyanalyzer
    ////fix types2pcmap getter and setter
    ////continue decoupling objects


    //types2uffparams_map.push_back(uffp);
    //types2uffparams_map.back().symbol=uff_params[0];
    ////cerr << "i=" << i << " " << uff_params[0] << endl;
    //types2uffparams_map.back().r1=aurostd::string2utype<double>(uff_params[1]);
    //types2uffparams_map.back().theta0=aurostd::string2utype<double>(uff_params[2]);
    //types2uffparams_map.back().x1=aurostd::string2utype<double>(uff_params[3]);
    //types2uffparams_map.back().D1=aurostd::string2utype<double>(uff_params[4]);
    //types2uffparams_map.back().zeta=aurostd::string2utype<double>(uff_params[5]);
    //types2uffparams_map.back().Z1=aurostd::string2utype<double>(uff_params[6]);
    //types2uffparams_map.back().Vi=aurostd::string2utype<double>(uff_params[7]);
    //types2uffparams_map.back().Uj=aurostd::string2utype<double>(uff_params[8]);
    //types2uffparams_map.back().ChiI=aurostd::string2utype<double>(uff_params[9]);
    //types2uffparams_map.back().hard=aurostd::string2utype<double>(uff_params[10]);
    //types2uffparams_map.back().radius=aurostd::string2utype<double>(uff_params[11]);
    //}
  }

  void POccCalculator::calculateSymNonPOccStructure(bool verbose) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    //CO - default should always be to turn everything on, modify input
    //get ROBUST symmetry determination
    xstr_sym=xstr_nopocc; //make copy to avoid carrying extra sym_ops (heavy) if not needed
    _kflags kflags_sym=m_kflags;  //make a copy so as to not overwrite m_kflags
    pflow::defaultKFlags4SymCalc(kflags_sym,true);
    //[this needs to be specified by m_kflags (do NOT override), otherwise LIB2RAW writes in LIB]pflow::defaultKFlags4SymWrite(kflags_sym,true);

    //streams
    ostream* p_oss_sym=p_oss;
    bool badbit_set=(*p_oss_sym).bad(); //((*p_oss_sym).rdstate()==std::ios_base::badbit);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " badbit_set=" << badbit_set << endl;}
    ofstream* p_FileMESSAGE_sym=p_FileMESSAGE;
    ofstream devNull("/dev/null");  //must be outside if statement (scope)
    //verbose=true;
    bool osswrite=verbose;
    if(!verbose){
      pflow::defaultKFlags4SymWrite(kflags_sym,false);  //setting to false is always ok
      if(!badbit_set){(*p_oss_sym).setstate(std::ios_base::badbit);}  //like NULL
      p_FileMESSAGE_sym=&devNull;
    }
    message << "Calculating symmetry of parent structure" << endl;
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE_sym,*p_oss_sym,_LOGGER_MESSAGE_);
    bool sym_calculated=pflow::PerformFullSymmetry(xstr_sym,*p_FileMESSAGE_sym,m_aflags,kflags_sym,osswrite,*p_oss_sym);
    if(!badbit_set){(*p_oss_sym).clear();}  //clear NULL
    if(sym_calculated==false){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"AFLOW-SYM failed to calculate the symmetry of the clean (non-pocc'd) structure");}

    if(LDEBUG) {
      xstructure xstr_sym2=xstr_sym;
      xstr_sym2.sortAtomsEquivalent();
      cerr << __AFLOW_FUNC__ << " non-pocc structure sorted" << endl;
      cerr << xstr_sym2 << endl;
    }

    //load iatoms info to m_pocc_sites
    for(uint site=0;site<m_pocc_sites.size();site++){
      m_pocc_sites[site].is_inequivalent=xstr_sym.atoms[site].is_inequivalent;
      m_pocc_sites[site].equivalent=xstr_sym.atoms[site].equivalent;
    }

  }

  void POccCalculator::propagateEquivalentAtoms2POccStructure() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    //CO - default should always be to turn everything on, modify input
    //assign iatoms info to pocc structure as well, allows sorting later
    //remember, the 0th occupant of xstr_pocc remains in xstr_sym (xstr_nopocc)
    for(uint site=0;site<m_pocc_sites.size();site++){
      for(uint atom=0;atom<m_pocc_sites[site].v_occupants.size();atom++){
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " mapping iatoms info of xstr_sym[" << site << "] to xstr_pocc[" << m_pocc_sites[site].v_occupants[atom] << "]" << endl;}
        xstr_pocc.atoms[m_pocc_sites[site].v_occupants[atom]].is_inequivalent=xstr_sym.atoms[site].is_inequivalent;  //not really important (or necessarily correct), but we propagate anyway
        xstr_pocc.atoms[m_pocc_sites[site].v_occupants[atom]].equivalent=xstr_sym.atoms[site].equivalent;  //allows sorting via sortAtomsEquiv, use stable_sort()
      }
    }
  }

  uint POccCalculator::getHNFTabelPOCCPrecision() const {return (uint)ceil(log10(1.0/min(xstr_pocc.partial_occupation_site_tol,xstr_pocc.partial_occupation_stoich_tol)));}
  uint POccCalculator::getHNFTableGeneralPrecision() const {return getHNFTabelPOCCPrecision();}
  uint POccCalculator::getHNFTableIterationPadding() const {return 4;}
  uint POccCalculator::getHNFTableErrorPadding() const {return getHNFTabelPOCCPrecision()+2;}
  uint POccCalculator::getHNFTableColumnPadding() const {return 2*getHNFTableGeneralPrecision()+2*getHNFTableIterationPadding()+15;}
  string POccCalculator::getHeaderMaxSiteError() const {return "site_error";}
  string POccCalculator::getHeaderMaxStoichError() const {return "stoich_error";}
  string POccCalculator::getHeaderStoichiometry() const {return "stoichiometry";}

  string POccCalculator::hnfTableHeader(){

    stringstream ss_header;

    //HEADER - START
    ss_header << aurostd::PaddedPRE("n",getHNFTableIterationPadding()) << "  " << string(getHNFTableErrorPadding(),'-') << " ";
    uint pocc_atom=1;   //start from 1
    string atom_name;
    for(uint site=0;site<m_pocc_sites.size();site++){
      if(m_pocc_sites[site].partial_occupation_flag){
        for(uint occ=0;occ<m_pocc_sites[site].v_occupants.size();occ++){
          atom_name=xstr_pocc.atoms[m_pocc_sites[site].v_occupants[occ]].name;
          if(atom_name.empty()){
            atom_name="type="+aurostd::utype2string(xstr_pocc.atoms[m_pocc_sites[site].v_occupants[occ]].type);
          }
          //output - START
          ss_header << "| ";
          ss_header << aurostd::PaddedCENTER("pocc_atom = " + aurostd::utype2string(pocc_atom) + "/" + aurostd::utype2string(pocc_atoms_total) + (!atom_name.empty()?" ["+atom_name+"]":""),getHNFTableColumnPadding()+2);  //need to +2, bug in PaddedCENTER()
          ss_header << " " ;  //CO20170629
          //output - END
          pocc_atom++;
        }
      }
    }
    ss_header << " | " << getHeaderMaxSiteError(); //header_max_site_error;
    ss_header << " | " << getHeaderMaxStoichError(); //header_max_stoich_error;
    //ss_header << " | " << getHeaderStoichiometry(); //CO20190131 - fix later
    ss_header << " " << endl;

    return ss_header.str();
  }

  string POccCalculator::hnfTableLineOutput(int i_hnf,int str_config){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    stringstream line_output;
    stringstream tmp_ss;tmp_ss.precision(getHNFTableGeneralPrecision());tmp_ss.setf(std::ios::fixed,std::ios::floatfield);

    double stoich_error,site_error=1.0/i_hnf;
    POccSiteConfiguration site_config;

    line_output << aurostd::PaddedPRE(aurostd::utype2string(i_hnf),getHNFTableIterationPadding());
    line_output << "  ";
    tmp_ss.str(""); tmp_ss << site_error;
    line_output << aurostd::PaddedPOST(tmp_ss.str(),getHNFTableErrorPadding());
    line_output << " ";

    stoich_error=v_str_configs[str_config].max_stoich_error;
    site_error=v_str_configs[str_config].max_site_error;
    if(LDEBUG) {
      if(v_str_configs.size()>1){
        //cerr << "CO CHECK IT OUT" << endl;
        for(uint str_config=0;str_config<v_str_configs.size();str_config++){
          cerr << "config " << str_config << endl;
          for(uint site=0;site<v_str_configs[str_config].site_configs.size();site++){
            cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_site_error=" << v_str_configs[str_config].site_configs[site].xv_site_error << endl;
            cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].max_site_error=" << v_str_configs[str_config].site_configs[site].max_site_error << endl;
            cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_occupation_count_input=" << v_str_configs[str_config].site_configs[site].xv_occupation_count_input << endl;
            cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_occupation_multiple=" << v_str_configs[str_config].site_configs[site].xv_occupation_multiple << endl;
            cerr << "v_str_configs[" << str_config << "].site_configs[" << site << "].xv_partial_occupation_value=" << v_str_configs[str_config].site_configs[site].xv_partial_occupation_value << endl;
            cerr << endl << endl;
          }
        }
      }  //i am interested in finding an example!
    }

    //print out results
    for(uint site=0;site<v_str_configs[str_config].site_configs.size();site++){
      site_config = v_str_configs[str_config].site_configs[site];
      if(site_config.isPartiallyOccupied()){
        for(int pocc_group=0;pocc_group<site_config.xv_occupation_count_input.rows;pocc_group++){
          for(int pocc_atom=0;pocc_atom<site_config.xv_occupation_count_input[pocc_group];pocc_atom++){
            //IMPORTANT CHECK
            if(!aurostd::isequal(site_config.xv_occupation_multiple[pocc_group]/((double)i_hnf),site_config.xv_partial_occupation_value[pocc_group],_ZERO_TOL_)){ //hard zero
              message << "site_config.xv_occupation_multiple[pocc_group]/((double)i_hnf)=" << site_config.xv_occupation_multiple[pocc_group]/((double)i_hnf);
              message << " != ";
              message << site_config.xv_partial_occupation_value[pocc_group] << "=site_config.xv_partial_occupation_value[pocc_group]";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
            }
            //output - START
            tmp_ss.str("");
            tmp_ss << aurostd::PaddedPRE(aurostd::utype2string(site_config.xv_occupation_multiple[pocc_group]),getHNFTableIterationPadding());
            tmp_ss << "/";
            tmp_ss << aurostd::PaddedPOST(aurostd::utype2string(i_hnf),getHNFTableIterationPadding());
            tmp_ss << " ";
            tmp_ss << site_config.xv_partial_occupation_value[pocc_group];
            tmp_ss << " ";
            tmp_ss << "(err = " << site_config.xv_site_error[pocc_group] << ")";
            line_output << "| " << aurostd::PaddedPOST(tmp_ss.str(),getHNFTableColumnPadding()) << " ";
            //output - END
          }
        }
      }
    }
    tmp_ss.str(""); tmp_ss << site_error;
    line_output << " | " << aurostd::PaddedPOST(tmp_ss.str(),getHeaderMaxSiteError().size()); //header_max_site_error.size()//getHNFTableErrorPadding());
    tmp_ss.str(""); tmp_ss << stoich_error;
    line_output << " | " << aurostd::PaddedPOST(tmp_ss.str(),getHeaderMaxStoichError().size()); //header_max_stoich_error.size()//getHNFTableErrorPadding());
    tmp_ss.str(""); tmp_ss << stoich_error;
    //line_output << " | " << aurostd::PaddedPOST(tmp_ss.str(),getHeaderStoichiometry().size()); //CO20190131 - fix later
    //tmp_ss.str(""); tmp_ss << stoich_error; //CO20190131 - fix later
    line_output << " " << endl;

    return line_output.str();

  }

  //void POccStructure::setHNFTablePadding(int _AFLOW_POCC_PRECISION_){
  //  //set some nice printing precisions and paddings
  //  hnf_table_general_precision=_AFLOW_POCC_PRECISION_;
  //  hnf_table_iteration_padding = 4;
  //  hnf_table_error_padding=_AFLOW_POCC_PRECISION_+2;
  //  hnf_table_column_padding=2*getHNFTableGeneralPrecision()+2*getHNFTableIterationPadding()+15;
  //  //header_max_stoich_error="stoich_error";
  //  //header_max_site_error="site_error";
  //}

  void POccCalculator::writeHNFTableOutput(int i_hnf,double& stoich_error,double& site_error) {
    stringstream message;
    stringstream multiple_configs_ss;

    if(v_str_configs.size()==0){
      message << "No structural configurations found";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }

    stoich_error=v_str_configs[0].max_stoich_error; //same for all
    site_error=v_str_configs[0].max_site_error;     //same for all

    int centering_padding=AFLOWIN_SEPARATION_LINE.size();
    for(uint str_config=0;str_config<v_str_configs.size();str_config++){
      if(!v_str_configs[str_config].site_configs.size()){
        message << "No site configurations found";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      centering_padding=hnfTableLineOutput(i_hnf,str_config).size();
      break;  //get one, they are all the same
    }

    //if(v_str_configs.size()>1){message_raw << multiple_configs_partiion << endl;}
    if(v_str_configs.size()>1){
      multiple_configs_ss.str("");
      multiple_configs_ss << "*** MULTIPLE CONFIGURATIONS POSSIBLE FOR HNF=" << i_hnf << " - START ***" << endl;
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,aurostd::PaddedCENTER(multiple_configs_ss.str(),centering_padding+2),m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    }

    for(uint str_config=0;str_config<v_str_configs.size();str_config++){
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,hnfTableLineOutput(i_hnf,str_config),m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    }

    if(v_str_configs.size()>1){
      multiple_configs_ss.str("");
      multiple_configs_ss << "*** MULTIPLE CONFIGURATIONS POSSIBLE FOR HNF=" << i_hnf << " - END ***" << endl;
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,aurostd::PaddedCENTER(multiple_configs_ss.str(),centering_padding+2),m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    }
  }

  //tol -> HNF optimizer
  void POccCalculator::calculateHNF(){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    stringstream message,message_raw;
    message_raw.precision(getHNFTableGeneralPrecision());
    message << "Fetching HNF value";
    pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    const double& stoich_tolerance=xstr_pocc.partial_occupation_stoich_tol;
    const double& site_tolerance=xstr_pocc.partial_occupation_site_tol;

    int i_hnf=0;
    double stoich_error,site_error;
    stoich_error=site_error=1.0;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " POcc structure:" << endl;
      cerr << xstr_pocc;
    }
    //[OLD, SHOW EVERYTHING NOW]n_hnf=xstr_pocc.partial_occupation_HNF;
    //[OLD, SHOW EVERYTHING NOW]message << "Using input HNF value = " << n_hnf;
    //[OLD, SHOW EVERYTHING NOW]pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    //[OLD, SHOW EVERYTHING NOW]//getSiteCountConfigurations(n_hnf,stoich_error);
    //[OLD, SHOW EVERYTHING NOW]getSiteCountConfigurations(n_hnf);
    //[OLD, SHOW EVERYTHING NOW]message_raw << AFLOWIN_SEPARATION_LINE << endl; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message_raw,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    //[OLD, SHOW EVERYTHING NOW]pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__,hnfTableHeader(),m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    //[OLD, SHOW EVERYTHING NOW]if(!writeHNFTableOutput(n_hnf,stoich_error,site_error)){return false;}
    //[OLD, SHOW EVERYTHING NOW]message_raw << AFLOWIN_SEPARATION_LINE << endl; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message_raw,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);

    bool hnf_already_provided=(xstr_pocc.partial_occupation_HNF!=0);

    if (hnf_already_provided) { //if HNF exists, no need to optimize pocc values
      message << "Using input HNF value = " << xstr_pocc.partial_occupation_HNF; 
      char LOGGER_TYPE=_LOGGER_MESSAGE_;
      if(m_p_flags.flag("HNF")){LOGGER_TYPE=_LOGGER_COMPLETE_;}
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,LOGGER_TYPE);
    } else {
      message << "Stoichiometry tolerance = " << stoich_tolerance; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      message << "Site tolerance = " << site_tolerance; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      message << "Optimizing HNF value"; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
    }

    i_hnf=0;
    message_raw << AFLOWIN_SEPARATION_LINE  << endl; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message_raw,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__,hnfTableHeader(),m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    while( hnf_already_provided ? i_hnf<xstr_pocc.partial_occupation_HNF : stoich_error>stoich_tolerance || site_error>site_tolerance ){
      i_hnf++;
      getSiteCountConfigurations(i_hnf);
      writeHNFTableOutput(i_hnf,stoich_error,site_error);
    }
    message_raw << AFLOWIN_SEPARATION_LINE << endl; pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message_raw,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    n_hnf=i_hnf;

    if(hnf_already_provided==false) { //if HNF exists, no need to optimize pocc values
      message << "Optimized HNF value = " << n_hnf;
      char LOGGER_TYPE=_LOGGER_MESSAGE_;
      if(m_p_flags.flag("HNF")){LOGGER_TYPE=_LOGGER_COMPLETE_;}
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,LOGGER_TYPE);
    }
  }

} //namespace pocc

namespace pocc {
  //--------------------------------------------------------------------------------
  // class POccUFFEnergyAnalyzer
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize();}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(pocc_flags,aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(pocc_flags,aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration,aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss) : POccCalculatorTemplate(),xStream(oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags,aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration,aflags);}
  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : POccCalculatorTemplate(),xStream(FileMESSAGE,oss),m_initialized(false) {initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags,aflags);}

  POccUFFEnergyAnalyzer::POccUFFEnergyAnalyzer(const POccUFFEnergyAnalyzer& b) : POccCalculatorTemplate(),xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  POccUFFEnergyAnalyzer::~POccUFFEnergyAnalyzer() {free();}

  void POccUFFEnergyAnalyzer::free() {
    POccCalculatorTemplate::free();
    m_initialized=false;
    hnf_mat.clear();
    for(uint i=0;i<m_types_config.size();i++){m_types_config[i].clear();} m_types_config.clear();
    types2uffparams_map.clear();
    m_vacancies.clear();
    m_exploration_radius=AUROSTD_MAX_DOUBLE;
    distance_matrix.clear();
    v_dist_nn.clear();
    xstr_ss.clear(); //DX20191220 - uppercase to lowercase clear
    sc2pc_map.clear();
    pc2sc_map.clear();
    xstr_cluster.clear(); //DX20191220 - uppercase to lowercase clear
    for(uint i=0;i<v_bonded_atom_indices.size();i++){v_bonded_atom_indices[i].clear();} v_bonded_atom_indices.clear();
    for(uint i=0;i<v_nonbonded_atom_indices.size();i++){v_nonbonded_atom_indices[i].clear();} v_nonbonded_atom_indices.clear();
    has_vacancies=false;
    bonding_set=false;
    m_energy_uff=AUROSTD_MAX_DOUBLE;
  }

  void POccUFFEnergyAnalyzer::copy(const POccUFFEnergyAnalyzer& b) { // copy PRIVATE
    POccCalculatorTemplate::copy(b);
    xStream::copy(b);
    m_initialized=b.m_initialized;
    hnf_mat=b.hnf_mat;
    for(uint i=0;i<m_types_config.size();i++){m_types_config[i].clear();} m_types_config.clear();for(uint i=0;i<b.m_types_config.size();i++){m_types_config.push_back(vector<int>(0));for(uint j=0;j<b.m_types_config[i].size();j++){m_types_config.back().push_back(b.m_types_config[i][j]);}}
    types2uffparams_map.clear();for(uint i=0;i<b.types2uffparams_map.size();i++){types2uffparams_map.push_back(b.types2uffparams_map[i]);}
    m_vacancies.clear();for(uint i=0;i<b.m_vacancies.size();i++){m_vacancies.push_back(b.m_vacancies[i]);}
    m_exploration_radius=b.m_exploration_radius;
    distance_matrix=b.distance_matrix;
    v_dist_nn.clear();for(uint i=0;i<b.v_dist_nn.size();i++){v_dist_nn.push_back(b.v_dist_nn[i]);}
    xstr_ss=b.xstr_ss;
    sc2pc_map.clear();for(uint i=0;i<b.sc2pc_map.size();i++){sc2pc_map.push_back(b.sc2pc_map[i]);}
    pc2sc_map.clear();for(uint i=0;i<b.pc2sc_map.size();i++){pc2sc_map.push_back(b.pc2sc_map[i]);}
    xstr_cluster=b.xstr_cluster;
    for(uint i=0;i<v_bonded_atom_indices.size();i++){v_bonded_atom_indices[i].clear();} v_bonded_atom_indices.clear();for(uint i=0;i<b.v_bonded_atom_indices.size();i++){v_bonded_atom_indices.push_back(vector<uint>(0));for(uint j=0;j<b.v_bonded_atom_indices[i].size();j++){v_bonded_atom_indices.back().push_back(b.v_bonded_atom_indices[i][j]);}}
    for(uint i=0;i<v_nonbonded_atom_indices.size();i++){v_nonbonded_atom_indices[i].clear();} v_nonbonded_atom_indices.clear();for(uint i=0;i<b.v_nonbonded_atom_indices.size();i++){v_nonbonded_atom_indices.push_back(vector<uint>(0));for(uint j=0;j<b.v_nonbonded_atom_indices[i].size();j++){v_nonbonded_atom_indices.back().push_back(b.v_nonbonded_atom_indices[i][j]);}}
    has_vacancies=b.has_vacancies;
    bonding_set=b.bonding_set;
    m_energy_uff=b.m_energy_uff;
  }

  const POccUFFEnergyAnalyzer& POccUFFEnergyAnalyzer::operator=(const POccUFFEnergyAnalyzer& b) {  // operator= PUBLIC
    if(this!=&b) {copy(b);}
    return *this;
  }

  void POccUFFEnergyAnalyzer::clear() {POccUFFEnergyAnalyzer a;copy(a);}  // clear PRIVATE

  bool POccUFFEnergyAnalyzer::initialize(ostream& oss){
    xStream::initialize(oss);
    return initialize();
  }
  bool POccUFFEnergyAnalyzer::initialize(const aurostd::xoption& pocc_flags,ostream& oss){
    xStream::initialize(oss);
    return initialize(pocc_flags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const _aflags& aflags,ostream& oss){
    xStream::initialize(oss);
    return initialize(aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss){
    xStream::initialize(oss);
    return initialize(pocc_flags,aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }
  bool POccUFFEnergyAnalyzer::initialize(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(pocc_flags,aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(){
    free();
    try{
      setExplorationRadius();
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const aurostd::xoption& pocc_flags){
    free();
    try{
      setPOccFlags(pocc_flags);
      setExplorationRadius();
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const _aflags& aflags){
    free();
    try{
      setAFlags(aflags);
      setExplorationRadius();
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags){
    free();
    try{
      setPOccFlags(pocc_flags);
      setAFlags(aflags);
      setExplorationRadius();
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ostream& oss){
    xStream::initialize(oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ostream& oss){
    xStream::initialize(oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ostream& oss){
    xStream::initialize(oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration,aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss){
    xStream::initialize(oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags,aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration,aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss){
    xStream::initialize(FileMESSAGE,oss);
    return initialize(xstr_pocc,xstr_nopocc,species_redecoration,pocc_flags,aflags);
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration){
    free();
    try{
      setPOccStructure(xstr_pocc);
      setNonPOccStructure(xstr_nopocc);
      setSpeciesRedecoration(species_redecoration);
      setExplorationRadius();
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags){
    free();
    try{
      setPOccFlags(pocc_flags);
      setPOccStructure(xstr_pocc);
      setNonPOccStructure(xstr_nopocc);
      setSpeciesRedecoration(species_redecoration);
      setExplorationRadius();
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags){
    free();
    try{
      setAFlags(aflags);
      setPOccStructure(xstr_pocc);
      setNonPOccStructure(xstr_nopocc);
      setSpeciesRedecoration(species_redecoration);
      setExplorationRadius();
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccUFFEnergyAnalyzer::initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags){
    free();
    try{
      setPOccFlags(pocc_flags);
      setAFlags(aflags);
      setPOccStructure(xstr_pocc);
      setNonPOccStructure(xstr_nopocc);
      setSpeciesRedecoration(species_redecoration);
      setExplorationRadius();
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  void POccUFFEnergyAnalyzer::calculateNNDistances(xstructure& xstr,vector<uint>& v_vacancies){
    xmatrix<double> _distance_matrix(xstr.atoms.size()-1,xstr.atoms.size()-1,0,0); distance_matrix=_distance_matrix;

    //get distance matrix first, then find nearest-neighbor distances
    //DX20190619 [OBSOLETE] xvector<double> min_vec; xvector<int> ijk;  //dummy

    //get distance matrix first
    for(uint atom1=0;atom1<xstr.atoms.size();atom1++){
      if(isVacancy(v_vacancies,atom1)){continue;}
      for(uint atom2=atom1+1;atom2<xstr.atoms.size();atom2++){
        if(isVacancy(v_vacancies,atom2)){continue;}
        distance_matrix(atom1,atom2)=distance_matrix(atom2,atom1)=aurostd::modulus(SYM::minimizeDistanceCartesianMethod(xstr.atoms[atom1].cpos,xstr.atoms[atom2].cpos,xstr.lattice)); //DX20190619 - updated function name and take the modulus of minimum cpos distance
      }
    }

    //now get n-n distances
    v_dist_nn.clear();
    for(uint i=0;i<xstr.atoms.size();i++){v_dist_nn.push_back(std::numeric_limits<double>::max());} //initialize with big number
    for(uint atom1=0;atom1<xstr.atoms.size();atom1++){
      if(isVacancy(v_vacancies,atom1)){continue;}
      for(uint atom2=0;atom2<xstr.atoms.size();atom2++){
        if(isVacancy(v_vacancies,atom2)){continue;}
        if(atom1!=atom2 && distance_matrix(atom1,atom2)<v_dist_nn[atom1]){v_dist_nn[atom1]=distance_matrix(atom1,atom2);}
      }
    }
  }

  //   void POccUFFEnergyAnalyzer::initialize(const POccCalculator& p_calc){
  //   setPOccFlags(p_calc.m_p_flags);
  //   setPOccStructure(p_calc.xstr_pocc);
  //   setNonPOccStructure(p_calc.xstr_nopocc);
  //   types2pc_map=getTypes2PCMap(xstr_pocc);     //get types2pc_map
  //   types2uffparams_map=getTypes2UFFParamsMap(m_species_redecoration);
  //   setExplorationRadius();
  //   }

  void POccUFFEnergyAnalyzer::setSpeciesRedecoration(const vector<string>& species_redecoration){
    POccCalculatorTemplate::setSpeciesRedecoration(species_redecoration);
    types2uffparams_map=getTypes2UFFParamsMap(m_species_redecoration);
  }

  void POccUFFEnergyAnalyzer::setExplorationRadius(){
    //[CO20181226 - need to experiment]double m_exploration_radius=RadiusSphereLattice(getLattice())*1.5; //figure this out CO, only works if energy radius = 10
    m_exploration_radius=DEFAULT_UFF_CLUSTER_RADIUS;
    //if(SET_KESONG_STANDARD_DIST){m_exploration_radius=ENERGY_RADIUS;}
    if(COMPARE_WITH_KESONG){m_exploration_radius=ENERGY_RADIUS;}
    //m_exploration_radius=12;
  }

  void POccUFFEnergyAnalyzer::getCluster(xmatrix<double>& _hnf_mat){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    hnf_mat=_hnf_mat;
    sc2pc_map.clear(); pc2sc_map.clear();
    xstr_ss=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false,false,false);  //force_strict_pc2scMap==false is very important here, we want the FIRST equivalent atom
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " xstr_ss=" << endl;cerr << xstr_ss << endl;}


    //double radius=10;
    //xvector<int> dims(3);
    //double radius = ENERGY_RADIUS/2.0;  //LatticeDimensionSphere radius fits within 2*dim
    //double radius = m_exploration_radius/2.0;  //LatticeDimensionSphere radius fits within 2*dim
    //cerr << "exploration radius = " << m_exploration_radius << endl;
    //cerr << RadiusSphereLattice(xstr_ss.lattice) << endl;
    //safety, no need, gives 1x1x1 as long as radius >= 0.0
    //if(RadiusSphereLattice(xstr_ss.lattice)<radius){dims(1)=1;dims(2)=1;dims(3)=3;}
    //else {dims=LatticeDimensionSphere(xstr_ss.lattice,radius);}
    //cerr << xstr_ss.lattice << endl;
    //dims=LatticeDimensionSphere(xstr_ss.lattice,m_exploration_radius);
    //cerr << "dims = " << dims << endl;

    //xmatrix<double> slattice(3,3);
    //slattice(1,1)=dims(1);slattice(2,2)=dims(2);slattice(3,3)=dims(3);
    //cerr << "START CLUSTER" << endl;
    //cluster=GetCluster(xstr_ss,radius,ssc2sc_map,sc2ssc_map); //GetSuperCell(xstr_ss,slattice,ssc2sc_map,sc2ssc_map,false,false);
    //xstr_cluster=GetSuperCell(xstr_ss,slattice,ssc2sc_map,sc2ssc_map,false,false);

    //[OBSOLETE]ssc2sc_map.clear(); sc2ssc_map.clear();
    //[OBSOLETE]cluster=GetCluster(xstr_ss,2.0*m_exploration_radius,ssc2sc_map,sc2ssc_map);

    xstr_cluster=xstr_ss;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_cluster=" << endl;cerr << xstr_cluster << endl;
      cerr << __AFLOW_FUNC__ << " m_exploration_radius=" << m_exploration_radius << endl;
    }

    GenerateGridAtoms(xstr_cluster,m_exploration_radius);
    if(LDEBUG) {
      for(uint i=0;i<xstr_cluster.grid_atoms.size();i++){
        cerr << __AFLOW_FUNC__ << " xstr_cluster.grid_atoms[i]=(" << xstr_cluster.grid_atoms[i] << ")" << endl;
        cerr << __AFLOW_FUNC__ << " xstr_cluster.grid_atoms[i].cpos=" << xstr_cluster.grid_atoms[i].cpos << endl;
        cerr << __AFLOW_FUNC__ << " xstr_cluster.grid_atoms[i].fpos=" << xstr_cluster.grid_atoms[i].fpos << endl;
      }
    }

    bonding_set=false;

    //cerr << xstr_cluster.grid_atoms_number << endl;
    //cerr << "STOP CLUSTER" << endl;
    //cerr << "CLUSTER " << endl;

    //xstr_bonding=xstr_nopocc; //will be updated if vacancies present
  }

  //nice function pointers so we don't need to rewrite so much code, and for loops are fast
  uint POccUFFEnergyAnalyzer::NNDistancesMapPC(uint atom){return sc2pc_map[xstr_cluster.grid_atoms_sc2pcMap[atom]];}
  uint POccUFFEnergyAnalyzer::NNDistancesMapSC(uint atom){return xstr_cluster.grid_atoms_sc2pcMap[atom];}

  //shortcut for recalculating bonds - only if bonds have been calculated already
  //take in current v_types_config
  //we will generate a list of vacancies
  //then we compare with previously determined vacancies
  //if they match, do not recalculate
  void POccUFFEnergyAnalyzer::setBonds(const vector<vector<int> >& v_types_config){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    //cout << "m_types_config OLD ";
    //for(uint i=0;i<m_types_config.size();i++){for(uint j=0;j<m_types_config[i].size();j++){cout << m_types_config[i][j] << " ";} }
    //cout << endl;
    //cout << "v_types_config         NEW ";
    //for(uint i=0;i<v_types_config.size();i++){for(uint j=0;j<v_types_config[i].size();j++){cout << v_types_config[i][j] << " ";} }
    //cout << endl;

    m_types_config=v_types_config;
    //vector<uint> v_vacancies=getVacanciesSuperCell(xstr_ss, v_types_config, v_vacancies, false);  //do not modify xstr_ss (no replace_types, just get vacancies)
    vector<uint> v_vacancies=getVacanciesSuperCell(pc2sc_map,v_types_config);  //get vacancies //do not modify xstr_ss (no replace_types, just get vacancies)
    has_vacancies=(v_vacancies.size()>0);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " has_vacancies=" << has_vacancies << endl;}

    //cout << "m_vacancies OLD ";
    //for(uint i=0;i<m_vacancies.size();i++){cout << m_vacancies[i] << " ";}
    //cout << endl;
    //cout << "v_vacancies         NEW ";
    //for(uint i=0;i<v_vacancies.size();i++){cout << v_vacancies[i] << " ";}
    //cout << endl;

    xstructure* ptr_xstr_bonding=&xstr_nopocc;  //slight optimization, if no vacancies, use primitive cell n-n distances
    uint (POccUFFEnergyAnalyzer::*NNDistancesMap)(uint)=NULL; //nice function pointers so we don't need to rewrite so much code, and for loops are fast
    NNDistancesMap=&POccUFFEnergyAnalyzer::NNDistancesMapPC;

    if(has_vacancies){
      if(bonding_set){bonding_set=(v_vacancies==m_vacancies);}  //bonding_set==false when cluster is generated, check to recalculate only if has been calculated before
      ptr_xstr_bonding=&xstr_ss;
      NNDistancesMap=&POccUFFEnergyAnalyzer::NNDistancesMapSC;
      //xstr_bonding=xstr_ss; //copies it all over
      //removeVacancies(xstr_bonding,v_vacancies);
      //} else {
  }

  xstructure& xstr_bonding = *ptr_xstr_bonding;

  if(!bonding_set){
    //set m_types_config and m_vacancies!!!!!
    m_vacancies=v_vacancies;
    bonding_set=true;
    calculateNNDistances(xstr_bonding,m_vacancies); //m_vacancies is empty if xstr_bonding==xstr_nopocc, so this is okay

    //now get bonding indices for super-superstructure
    double distij;
    for(uint i=0;i<v_bonded_atom_indices.size();i++){v_bonded_atom_indices[i].clear();} v_bonded_atom_indices.clear();
    for(uint i=0;i<v_nonbonded_atom_indices.size();i++){v_nonbonded_atom_indices[i].clear();} v_nonbonded_atom_indices.clear();
    //double uff_bonding_distance=DEFAULT_UFF_BONDING_DISTANCE;
    //if(0){
    //for(uint atom1=0;atom1<xstr_cluster.atoms.size();atom1++){
    //  for(uint atom2=atom1+1;atom2<xstr_cluster.atoms.size();atom2++){
    //    distij=AtomDist(xstr_cluster.atoms[atom1],xstr_cluster.atoms[atom2]);
    //    //if(distij<=radius){
    //    if(abs(distij-v_dist_nn[sc2pc_map[xstr_cluster.grid_atoms_sc2pcMap[atom1]]])<uff_bonding_distance){
    //      v_bonded_atom_indices.push_back(vector<uint>(0));
    //      v_bonded_atom_indices.back().push_back(atom1);
    //      v_bonded_atom_indices.back().push_back(atom2);
    //      //cerr << "BONDING    " << atom1 << " " << atom2 << endl;
    //    } else {
    //      v_nonbonded_atom_indices.push_back(vector<uint>(0));
    //      v_nonbonded_atom_indices.back().push_back(atom1);
    //      v_nonbonded_atom_indices.back().push_back(atom2);
    //      //cerr << "NONBONDING " << atom1 << " " << atom2 << endl;
    //    }
    //    //}
    //  }
    //}
    //}
    uint atom1;
    double uff_bonding_distance=DEFAULT_UFF_BONDING_DISTANCE; //CO20190208 - remember, this is not a raw distance, but a distance to within the nn distance
    for(uint sc_atom1=0;sc_atom1<xstr_ss.atoms.size();sc_atom1++){
      if(isVacancy(m_vacancies,sc_atom1)){continue;}
      atom1=xstr_cluster.grid_atoms_pc2scMap[sc_atom1];
      for(uint atom2=0;atom2<(uint)xstr_cluster.grid_atoms_number;atom2++){
        if(atom1!=atom2){
          //distij=AtomDist(xstr_cluster.atoms[atom1],xstr_cluster.atoms[atom2]);
          distij=AtomDist(xstr_cluster.grid_atoms[atom1],xstr_cluster.grid_atoms[atom2]);
          if(distij<=m_exploration_radius){ //KY (via SC) used <= as well (GetNeighData)
            //  cerr << "OK cluster atom 1 " << atom1 << " " << xstr_cluster.grid_atoms[atom1].cpos << endl;
            //  cerr << "OK cluster atom 2 " << atom2 << " " << xstr_cluster.grid_atoms[atom2].cpos << endl;
            if(abs(distij-v_dist_nn[(this->*NNDistancesMap)(atom1)])<=uff_bonding_distance){  //KY only looked at <, not <= (ExtractBonds)
              v_bonded_atom_indices.push_back(vector<uint>(0));
              v_bonded_atom_indices.back().push_back(atom1);
              v_bonded_atom_indices.back().push_back(atom2);
              //cerr << "BOND " << atom1 << " " << atom2 << endl;
              //count++;
              //cerr << "BONDING    " << atom1 << " " << atom2 << endl;
            } else {
              v_nonbonded_atom_indices.push_back(vector<uint>(0));
              v_nonbonded_atom_indices.back().push_back(atom1);
              v_nonbonded_atom_indices.back().push_back(atom2);
              //cerr << "NONBONDING " << atom1 << " " << atom2 << endl;
            }
          }
        }
      }
      //cerr << "_atom1=" << _atom1 << " " << count << endl;
    }

    //for(uint i=0;i<v_dist_nn.size();i++){cerr << "i=" << i << " " << v_dist_nn[i] << endl;}
    //cerr << distance_matrix << endl;

    //cerr << *xstr_nopocc << endl;
    //cerr << xstr_ss << endl;
    //cerr << dims << endl;
    //cerr << xstr_cluster.atoms << endl;
    //cerr << "DONE" << endl;
  }
  if(LDEBUG) {cerr << __AFLOW_FUNC__ << " bonding set" << endl;}
  }

  bool POccUFFEnergyAnalyzer::isVacancy(vector<uint>& v_vacancies, uint atom){ //this function must remain fast
    for(uint vi=0,v_vacancies_size=v_vacancies.size();vi<v_vacancies_size;vi++){if(atom==v_vacancies[vi]){return true;}}
    return false;
  }

  double POccUFFEnergyAnalyzer::bondEnergyBond(const UFFParamBond& uffb){return 0.5 * uffb.Kij * uffb.delta * uffb.delta;}
  double POccUFFEnergyAnalyzer::bondEnergyNoBond(const UFFParamBond& uffb){return uffb.Dij * (uffb.X12 - 2.0 * uffb.X6);}

  double POccUFFEnergyAnalyzer::getUFFBondEnergy(xstructure& xstr,vector<vector<uint> >& v_bonded_atom_indices,uint MODE) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    double energy=0.0;
    double distij;
    uint atom1,atom2,type1,type2;
    UFFParamAtom uffai, uffaj;
    UFFParamBond uffb;
    uint count=0;

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " bonding mode = " << (MODE==BOND_MODE ? "bonding" : "non-bonding") << endl;}

    //pointer's for bond/non-bonded energy contributions
    double (POccUFFEnergyAnalyzer::*BondEnergy)(const UFFParamBond& uffb)=NULL;
    if(MODE==BOND_MODE){BondEnergy=&POccUFFEnergyAnalyzer::bondEnergyBond;}  //spring potential
    else if(MODE==NONBOND_MODE){BondEnergy=&POccUFFEnergyAnalyzer::bondEnergyNoBond;}  //lennard-jones potential
    else {
      message << "Invalid energy mode (MODE=" << MODE << ")";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " v_bonded_atom_indices.size()=" << v_bonded_atom_indices.size() << endl;}

    if(types2uffparams_map.size()==0){
      message << "types2uffparams_map.size()==0";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
    }

    for(uint i=0;i<v_bonded_atom_indices.size();i++){

      atom1=xstr_cluster.grid_atoms_sc2pcMap[v_bonded_atom_indices[i][0]];
      if(isVacancy(m_vacancies,atom1)){continue;}

      atom2=xstr_cluster.grid_atoms_sc2pcMap[v_bonded_atom_indices[i][1]];
      if(isVacancy(m_vacancies,atom2)){continue;}

      //cerr << "HERE2" << endl;
      //found_vacancy=false;
      //skip anything with a vacancy - START
      //for(uint vi=0;vi<v_vacancies.size() && !found_vacancy;vi++){if(atom1==v_vacancies[vi] || atom2==v_vacancies[vi]){found_vacancy=true;}}
      //if(found_vacancy){
      //cerr << endl;
      //cerr << "SKIPPING!" << endl;
      //cerr << "atom1 " << v_bonded_atom_indices[i][0] << "(" << atom1 << ")  " << "(cluster " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]].cpos << ") (supercell " << xstr.atoms[atom1].cpos << ")" << endl;
      //cerr << "atom3 " << v_bonded_atom_indices[i][1] << "(" << atom2 << ")  " << "(cluster " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]].cpos << ") (supercell " << xstr.atoms[atom2].cpos << ")" << endl;
      //cerr << "SKIPPING!" << endl;
      //cerr << endl;
      //  continue;
      //}
      count++;
      //skip anything with a vacancy - END
      //cerr << "atom1=" << atom1 << endl;
      //cerr << "atom2=" << atom2 << endl;
      //cerr << xstr.atoms.size() << endl;
      type1=xstr.atoms[atom1].type;
      type2=xstr.atoms[atom2].type;
      //cerr << "v_bonded_atom_indices[i][0]=" << v_bonded_atom_indices[i][0] << endl;
      //cerr << "v_bonded_atom_indices[i][1]=" << v_bonded_atom_indices[i][1] << endl;
      //cerr << xstr_cluster.grid_atoms_number << endl;
      //distij=AtomDist(xstr_cluster.atoms[v_bonded_atom_indices[i][0]],xstr_cluster.atoms[v_bonded_atom_indices[i][1]]);
      distij=AtomDist(xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]],xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]]);
      if(distij>m_exploration_radius){
        message << "Attempting to explore atoms outside of the search radius ";
        message << "(";
        message << "cluster atom 1 " << v_bonded_atom_indices[i][0] << " " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]].cpos;
        message << " vs. ";
        message << "cluster atom 2 " << v_bonded_atom_indices[i][1] << " " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]].cpos;
        message << ")";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message);
      }
      //cerr << "OK " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][0]].cpos << " " << xstr_cluster.grid_atoms[v_bonded_atom_indices[i][1]]  << " " ;
      //cerr << "DIST " << distij <<endl;
      //cerr << "HERE4" << endl;
      //cerr << "sss atom1=" << v_bonded_atom_indices[i][0] << " type=" << xstr_cluster.atoms[v_bonded_atom_indices[i][0]].type << " ";
      //cerr << "sss atom2=" << v_bonded_atom_indices[i][1] << " type=" << xstr_cluster.atoms[v_bonded_atom_indices[i][1]].type << endl;
      //cerr << "ss  atom1=" << atom1 << " type=" << type1 << " ";
      //cerr << "ss  atom2=" << atom2 << " type=" << type2 << endl;
      //cerr << "distij=" << distij << " " << (MODE==BOND_MODE ? "BOND" : "NONBOND") << endl;
      uffai=types2uffparams_map[type1];
      uffaj=types2uffparams_map[type2];

#if _DEBUG_POCC_CLUSTER_ANALYSIS_
      cerr << __AFLOW_FUNC__ << " uffai.symbol=" << uffai.symbol << endl;
      cerr << __AFLOW_FUNC__ << " uffaj.symbol=" << uffaj.symbol << endl;
#endif

      uffb.calculate(uffai,uffaj,distij);
      energy+=(this->*BondEnergy)(uffb);

#if _DEBUG_POCC_CLUSTER_ANALYSIS_
      //spring potential
      cerr << __AFLOW_FUNC__ << " spring potential debug" << endl;
      cerr << __AFLOW_FUNC__ << " atomi=" << xstr.atoms[atom1] << endl;
      cerr << __AFLOW_FUNC__ << " atomj=" << xstr.atoms[atom2] << endl;
      cerr << __AFLOW_FUNC__ << " distij=" << distij << endl;
      cerr << __AFLOW_FUNC__ << " uffb.Kij=" << uffb.Kij << endl;
      cerr << __AFLOW_FUNC__ << " uffb.delta=" << uffb.delta << endl;
      cerr << __AFLOW_FUNC__ << " uffb.ren=" << uffb.ren << endl;
      cerr << __AFLOW_FUNC__ << " uffb.R0=" << uffb.R0 << endl;
      cerr << __AFLOW_FUNC__ << " energy=" << kcal2eV * 0.5 * uffb.Kij * uffb.delta * uffb.delta << endl;

      //lennard-jones potential
      cerr << __AFLOW_FUNC__ << " lennard-jones potential debug" << endl;
      cerr << __AFLOW_FUNC__ << " type1=" << type1 << endl;
      cerr << __AFLOW_FUNC__ << " type2=" << type2 << endl;
      //cerr << __AFLOW_FUNC__ << " types2uffparams_map[type1]=" << types2uffparams_map[type1] << endl;
      //cerr << __AFLOW_FUNC__ << " types2uffparams_map[type2]=" << types2uffparams_map[type2] << endl;
      cerr << __AFLOW_FUNC__ << " atomi=" << xstr.atoms[atom1] << endl;
      cerr << __AFLOW_FUNC__ << " atomj=" << xstr.atoms[atom2] << endl;
      cerr << __AFLOW_FUNC__ << " atomi.type=" << xstr.atoms[atom2].type << endl;
      cerr << __AFLOW_FUNC__ << " atomj.type=" << xstr.atoms[atom2].type << endl;
      cerr << __AFLOW_FUNC__ << " distij=" << distij << endl;
      cerr << __AFLOW_FUNC__ << " energy=" << kcal2eV * uffb.Dij * (uffb.X12 - 2.0 * uffb.X6) << endl;
#endif

    }
    //if(MODE==BOND_MODE) cerr << "BONDING LENGTH C " << count << endl;
    energy*=kcal2eV;
    return energy;
  }

  xstructure POccCalculator::createXStructure(const POccSuperCell& psc,int n_hnf,
      unsigned long long int hnf_count,
      unsigned long long int types_config_permutations_count,
      bool clean_structure,bool primitivize){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    vector<int> sc2pc_map,pc2sc_map;

    xmatrix<double> hnf_mat;
    vector<vector<int> > v_types_config;
    vector<string> vstmp;
    getHNFMatSiteConfig(psc,hnf_mat,v_types_config);
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " hnf_mat=" << endl;
      cerr << hnf_mat << endl;
      vstmp.clear();
      for(uint site=0;site<v_types_config.size();site++){vstmp.push_back(aurostd::joinWDelimiter(v_types_config[site],","));}
      cerr << __AFLOW_FUNC__ << " v_types_config=" << aurostd::joinWDelimiter(vstmp,";") << endl;
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " creating superstructure" << endl;}
    xstructure supercell=GetSuperCell(xstr_nopocc,hnf_mat,sc2pc_map,pc2sc_map,false,false,false,false); //force_strict_pc2scMap==false is very important here, we want the FIRST equivalent atom
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " sc2pc_map.size()=" << sc2pc_map.size() << endl;
      cerr << __AFLOW_FUNC__ << " pc2sc_map.size()=" << pc2sc_map.size() << endl;
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " resolving vacancies (if any)" << endl;}
    vector<uint> v_vacancies=getVacanciesSuperCell(pc2sc_map,v_types_config);  //get vacancies
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " vacancies=" << aurostd::joinWDelimiter(v_vacancies,",") << endl;}
    replaceRandomSitesSuperCell(xstr_pocc,types2pc_map,pc2sc_map,v_types_config,supercell); //replace random sites (ignore vacancies for now, deal in rebuildSuperCell ONLY)

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " rebuilding structure" << endl;}
    rebuildSuperCell(xstr_pocc,v_vacancies,supercell);    //remove vacancies and clean all xstructure properties
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structure built!" << endl;}

    if(clean_structure){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structure pre-cleaning" << endl;cerr << supercell << endl;}
      if(primitivize){
        supercell.GetStandardPrimitive();
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structure post-primitivization" << endl;cerr << supercell << endl;}
      }
      supercell.ReScale(1.0);
      supercell.ShiftOriginToAtom(0);
      supercell.BringInCell();
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structure post-cleaning" << endl;cerr << supercell << endl;}
    }
    //[HNF(4)=7/7= 1 0 0; 1 2 0; 1 0 2]
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " rewriting title" << endl;}
    supercell.title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr_pocc.title); //undo "Standardized (re-decorated) parent structure"
    supercell.title+=" H"+aurostd::utype2string(psc.m_hnf_index)+"C"+aurostd::utype2string(psc.m_site_config_index);
    supercell.title+=" HNF(n="+aurostd::utype2string(n_hnf)+","+aurostd::utype2string(psc.m_hnf_index+1)+"/"+aurostd::utype2string(hnf_count)+")";
    supercell.title+="=["+aurostd::utype2string(hnf_mat(1,1))+" "+aurostd::utype2string(hnf_mat(1,2))+" "+aurostd::utype2string(hnf_mat(1,3))+";";
    supercell.title+= " "+aurostd::utype2string(hnf_mat(2,1))+" "+aurostd::utype2string(hnf_mat(2,2))+" "+aurostd::utype2string(hnf_mat(2,3))+";";
    supercell.title+= " "+aurostd::utype2string(hnf_mat(3,1))+" "+aurostd::utype2string(hnf_mat(3,2))+" "+aurostd::utype2string(hnf_mat(3,3))+"]";
    supercell.title+=" site_config("+aurostd::utype2string(psc.m_site_config_index+1)+"/"+aurostd::utype2string(types_config_permutations_count)+")";
    supercell.title+="=[";
    for(uint site=0;site<v_types_config.size();site++){
      for(uint atom=0;atom<v_types_config[site].size();atom++){
        supercell.title+=(site==0&&atom==0?"":" ")+aurostd::utype2string(v_types_config[site][atom]);
      }
      supercell.title+=((site!=v_types_config.size()-1)?";":"");
    }
    supercell.title+="]";
    supercell.title+=" DG="+aurostd::utype2string(psc.m_degeneracy);
    supercell.title+=string(DOI_POCC);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " final structure" << endl;cerr << supercell << endl;}
    return supercell;
  }

  double POccUFFEnergyAnalyzer::getUFFEnergy() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    stringstream message;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " starting" << endl;}

    xstructure supercell(xstr_ss);
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " supercell(start)=" << endl;cerr << supercell << endl;}
    replaceRandomSitesSuperCell(xstr_pocc,types2pc_map,pc2sc_map,m_types_config,supercell); //replace random sites (ignore vacancies for now, deal in rebuildSuperCell ONLY)
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " supercell(random-sites-only)=" << endl;cerr << supercell << endl;}

    m_energy_uff=0.0;
    m_energy_uff+=getUFFBondEnergy(supercell,v_bonded_atom_indices,BOND_MODE);
    //cerr << "BONDING ENERGY C " << energy << endl;
    m_energy_uff+=getUFFBondEnergy(supercell,v_nonbonded_atom_indices,NONBOND_MODE);

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " structure energy=" << m_energy_uff << endl;}

    if(COMPARE_WITH_KESONG){
      cout << "COREY  " << std::fixed << std::setprecision(15) << m_energy_uff << endl;
      //if(SET_KESONG_STANDARD_DIST){
      rebuildSuperCell(xstr_pocc,m_vacancies,supercell); 
      double uff_energy_kesong=pocc::CalculateUFFEnergy(supercell); 
      cout << "KESONG " << std::fixed << std::setprecision(15) << uff_energy_kesong << endl; 
      double diff=abs(m_energy_uff-uff_energy_kesong);
      cout << "DIFF " << std::fixed << std::setprecision(15) << diff << endl;
      if(diff>=DEFAULT_UFF_ENERGY_TOLERANCE){
        message << "Difference in energy between pocc_old and pocc_new schemes is TOO large, please check.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_VALUE_ILLEGAL_);
      }
      cout << supercell << endl;
      //}
    }

    return m_energy_uff;
  }

} // namespace pocc

namespace pocc {
  //--------------------------------------------------------------------------------
  // class POccStructuresFile
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  POccStructuresFile::POccStructuresFile(ostream& oss): xStream(oss),m_initialized(false) {initialize();}
  POccStructuresFile::POccStructuresFile(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize();}
  POccStructuresFile::POccStructuresFile(const string& fileIN,ostream& oss) : xStream(oss),m_initialized(false) {initialize(fileIN);}
  POccStructuresFile::POccStructuresFile(const string& fileIN,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(fileIN);}
  POccStructuresFile::POccStructuresFile(const _aflags& aflags,ostream& oss): xStream(oss),m_initialized(false) {initialize(aflags);}
  POccStructuresFile::POccStructuresFile(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(aflags);}
  POccStructuresFile::POccStructuresFile(const string& fileIN,const _aflags& aflags,ostream& oss) : xStream(oss),m_initialized(false) {initialize(fileIN,aflags);}
  POccStructuresFile::POccStructuresFile(const string& fileIN,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {initialize(fileIN,aflags);}
  POccStructuresFile::POccStructuresFile(const POccStructuresFile& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  POccStructuresFile::~POccStructuresFile() {xStream::free();free();}

  const POccStructuresFile& POccStructuresFile::operator=(const POccStructuresFile& b) {  // operator= PUBLIC
    if(this!=&b) {copy(b);}
    return *this;
  }

  void POccStructuresFile::clear() {POccStructuresFile a;copy(a);}

  void POccStructuresFile::free() {
    m_initialized=false;
    m_content="";
    m_vcontent.clear();
    m_filename="";
    m_aflags.clear();
    //m_ARUN_directories.clear();
    l_supercell_sets.clear();
    m_fileoptions.clear();
    m_vPOSCAR_lines.clear();
  }

  void POccStructuresFile::copy(const POccStructuresFile& b) { // copy PRIVATE
    xStream::copy(b);
    uint i=0,j=0,k=0;
    m_initialized=b.m_initialized;
    m_content=b.m_content;
    m_vcontent.clear();for(i=0;i<b.m_vcontent.size();i++){m_vcontent.push_back(b.m_vcontent[i]);}
    m_filename=b.m_filename;
    m_aflags=b.m_aflags;
    //m_ARUN_directories.clear();for(i=0;i<b.m_ARUN_directories.size();i++){m_ARUN_directories.push_back(b.m_ARUN_directories[i]);}
    l_supercell_sets.clear();for(std::list<POccSuperCellSet>::const_iterator it=b.l_supercell_sets.begin();it!=b.l_supercell_sets.end();++it){l_supercell_sets.push_back(*it);}
    m_fileoptions=b.m_fileoptions;
    m_vPOSCAR_lines.clear();
    for(i=0;i<b.m_vPOSCAR_lines.size();i++){
      m_vPOSCAR_lines.push_back(vector<vector<uint> >(0));
      for(j=0;j<b.m_vPOSCAR_lines[i].size();j++){
        m_vPOSCAR_lines[i].push_back(vector<uint>(0));
        for(k=0;k<b.m_vPOSCAR_lines[i][j].size();k++){
          m_vPOSCAR_lines[i][j].push_back(b.m_vPOSCAR_lines[i][j][k]);
        }
      }
    }
  }

  bool POccStructuresFile::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }
  bool POccStructuresFile::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }
  bool POccStructuresFile::initialize() {
    free();
    m_initialized=false;  //no point
    return m_initialized;
  }
  bool POccStructuresFile::initialize(const string& fileIN,ostream& oss) {
    xStream::initialize(oss);
    return initialize(fileIN);
  }
  bool POccStructuresFile::initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(fileIN);
  }
  bool POccStructuresFile::initialize(const string& fileIN) {
    free();
    try{
      readFile(fileIN);
      processFile();
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccStructuresFile::initialize(const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(aflags);
  }
  bool POccStructuresFile::initialize(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(aflags);
  }
  bool POccStructuresFile::initialize(const _aflags& aflags) {
    free();
    try{
      setAFlags(aflags);
      m_initialized=false;  //no point
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }
  bool POccStructuresFile::initialize(const string& fileIN,const _aflags& aflags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(fileIN,aflags);
  }
  bool POccStructuresFile::initialize(const string& fileIN,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(fileIN,aflags);
  }
  bool POccStructuresFile::initialize(const string& fileIN,const _aflags& aflags) {
    free();
    try{
      setAFlags(aflags);
      readFile(fileIN);
      processFile();
      m_initialized=true;
    }
    catch(aurostd::xerror& err){pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), m_aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);}
    return m_initialized;
  }

  void POccStructuresFile::setAFlags(const _aflags& aflags) {m_aflags=aflags;}

  void POccStructuresFile::readFile(const string& fileIN) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    m_filename="";
    stringstream structures_file_ss;

    bool found_file=false;
    if(!found_file && aurostd::EFileExist(m_aflags.Directory+"/"+fileIN,m_filename)){aurostd::efile2stringstream(m_filename,structures_file_ss);found_file=true;}
    if(!found_file && aurostd::FileExist(m_aflags.Directory+"/"+fileIN,m_filename)){aurostd::file2stringstream(m_filename,structures_file_ss);found_file=true;}
    if(!found_file){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,m_aflags.Directory+"/"+fileIN+" does not exist",_FILE_NOT_FOUND_);}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_filename=" << m_filename << endl;}

    m_content=structures_file_ss.str();
    aurostd::string2vectorstring(m_content,m_vcontent);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " m_vcontent.size()=" << m_vcontent.size() << endl;}
  }

  void POccStructuresFile::processFile() {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    //get poscar_start_tag and poscar_stop_tag
    string POccStructure_title_tag="HNF(n="; //in case we cannot find POSCAR_START_tag, SHACHAR OLD FORMAT
    string poscar_start_tag=POSCAR_START_tag;
    bool found_POSCAR_START_tag=aurostd::substring2bool(m_content,poscar_start_tag);
    if(!found_POSCAR_START_tag){poscar_start_tag=POccStructure_title_tag;}
    if(!aurostd::substring2bool(m_content,poscar_start_tag)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find POSCAR start tag",_FILE_CORRUPT_);}
    string poscar_stop_tag=POSCAR_STOP_tag;
    bool found_POSCAR_STOP_tag=aurostd::substring2bool(m_content,poscar_stop_tag);
    if(!found_POSCAR_STOP_tag){poscar_stop_tag="[AFLOW] **";}
    if(!aurostd::substring2bool(m_content,poscar_stop_tag)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot find POSCAR stop tag",_FILE_CORRUPT_);}
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " poscar_start_tag=\"" << poscar_start_tag << "\"" << endl;
      cerr << __AFLOW_FUNC__ << " poscar_stop_tag=\"" << poscar_stop_tag << "\"" << endl;
    }

    //get structure_group_tag and m_fileoptions.flag("TYPE::XXXX")
    bool found_tag_sg=false;
    string structure_group_tag="";

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " PRE" << endl;
      cerr << __AFLOW_FUNC__ << " structure_group_tag=" << structure_group_tag << endl;
      cerr << __AFLOW_FUNC__ << " m_fileoptions.flag(\"TYPE::ALL_SUPERCELLS\")=" << m_fileoptions.flag("TYPE::ALL_SUPERCELLS") << endl;
      cerr << __AFLOW_FUNC__ << " m_fileoptions.flag(\"TYPE::UNIQUE_SUPERCELLS\")=" << m_fileoptions.flag("TYPE::UNIQUE_SUPERCELLS") << endl;
    }

    //STRUCTURES_GROUP
    if(found_tag_sg==false){
      structure_group_tag=POCC_AFLOWIN_tag+"STRUCTURES_GROUP";
      if(aurostd::substring2bool(m_content,structure_group_tag)){
        m_fileoptions.flag("TYPE::ALL_SUPERCELLS",TRUE);
        found_tag_sg=true;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " \"STRUCTURES_GROUP\" tag found" << endl;}
      }
    }
    //STRUCTURE_GROUP
    if(found_tag_sg==false){
      structure_group_tag=POCC_AFLOWIN_tag+"STRUCTURE_GROUP";
      if(aurostd::substring2bool(m_content,structure_group_tag)){
        m_fileoptions.flag("TYPE::ALL_SUPERCELLS",TRUE);
        found_tag_sg=true;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " \"STRUCTURE_GROUP\" tag found" << endl;}
      }
    }
    //STRUCTURE GROUP
    if(found_tag_sg==false){
      structure_group_tag=POCC_AFLOWIN_tag+"STRUCTURE GROUP";
      if(aurostd::substring2bool(m_content,structure_group_tag)){
        m_fileoptions.flag("TYPE::ALL_SUPERCELLS",TRUE);
        found_tag_sg=true;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " \"STRUCTURE GROUP\" tag found" << endl;}
      }
    }
    //TYPE::UNIQUE_SUPERCELLS
    if(found_tag_sg==false){
      structure_group_tag=poscar_start_tag;
      m_fileoptions.flag("TYPE::UNIQUE_SUPERCELLS",TRUE);
    }

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " POST" << endl;
      cerr << __AFLOW_FUNC__ << " structure_group_tag=" << structure_group_tag << endl;
      cerr << __AFLOW_FUNC__ << " m_fileoptions.flag(\"TYPE::ALL_SUPERCELLS\")=" << m_fileoptions.flag("TYPE::ALL_SUPERCELLS") << endl;
      cerr << __AFLOW_FUNC__ << " m_fileoptions.flag(\"TYPE::UNIQUE_SUPERCELLS\")=" << m_fileoptions.flag("TYPE::UNIQUE_SUPERCELLS") << endl;
    }

    if(!(m_fileoptions.flag("TYPE::ALL_SUPERCELLS") || m_fileoptions.flag("TYPE::UNIQUE_SUPERCELLS"))){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown structures file type (none set)",_FILE_CORRUPT_);}
    if(m_fileoptions.flag("TYPE::ALL_SUPERCELLS") && m_fileoptions.flag("TYPE::UNIQUE_SUPERCELLS")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown structures file type (all set)",_FILE_CORRUPT_);}

    unsigned long long int iline=0;
    unsigned long long int ititle=0;
    bool reading_POSCAR=false;
    string pocc_hash="",hnf_index_str="",site_config_index_str="";
    for(iline=0;iline<m_vcontent.size();iline++){
      if(aurostd::substring2bool(m_vcontent[iline],structure_group_tag)){
        l_supercell_sets.push_back(POccSuperCellSet());
        m_vPOSCAR_lines.push_back(vector<vector<uint> >(0));
      }
      if(aurostd::substring2bool(m_vcontent[iline],poscar_start_tag)){  //do NOT use "else if" as structure_group_tag might equal poscar_start_tag
        reading_POSCAR=true;
        if(l_supercell_sets.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found POSCAR before structures tag",_FILE_CORRUPT_);}
        l_supercell_sets.back().m_psc_set.push_back(POccSuperCell());
        if(m_vPOSCAR_lines.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found POSCAR stop tag before start tag",_FILE_CORRUPT_);}
        m_vPOSCAR_lines.back().push_back(vector<uint>(0));
        ititle=(found_POSCAR_START_tag ? iline+1 : iline);
        if(ititle>=m_vcontent.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot grab title at iline="+aurostd::utype2string(ititle),_FILE_CORRUPT_);}
        m_vPOSCAR_lines.back().back().push_back( ititle );
        //process title
        parsePOccHashFromXStructureTitle(m_vcontent[ititle],pocc_hash,hnf_index_str,site_config_index_str);
        l_supercell_sets.back().m_psc_set.back().m_hnf_index=aurostd::string2utype<unsigned long long int>(hnf_index_str);
        l_supercell_sets.back().m_psc_set.back().m_site_config_index=aurostd::string2utype<unsigned long long int>(site_config_index_str);
        l_supercell_sets.back().m_psc_set.back().m_degeneracy=getDGFromXStructureTitle(m_vcontent[ititle]);
      }
      if(reading_POSCAR && aurostd::substring2bool(m_vcontent[iline],poscar_stop_tag)){
        reading_POSCAR=false;
        if(m_vPOSCAR_lines.size()==0 || m_vPOSCAR_lines.back().size()==0 || m_vPOSCAR_lines.back().back().size()!=1){
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Found POSCAR stop tag before start tag",_FILE_CORRUPT_);
        }
        m_vPOSCAR_lines.back().back().push_back( found_POSCAR_STOP_tag ? iline : iline-1 );
      }
      //do NOT read UFF_ENERGY yet, as it might be above the POSCAR, we need to go through m_vcontent again after populating l_supercell_sets sizes
    }

    reading_POSCAR=false;
    bool found_first=false;
    std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();
    unsigned long long int index_str=0;  //index of structures within structures group
    for(iline=0;iline<m_vcontent.size();iline++){
      if(aurostd::substring2bool(m_vcontent[iline],structure_group_tag)){if(found_first==true){++it;}found_first=true;index_str=0;}
      if(aurostd::substring2bool(m_vcontent[iline],POCC_AFLOWIN_tag+"UFF_ENERGY=")){
        if(index_str>=(*it).m_psc_set.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"index_str>=(*it).m_psc_set.size()",_INDEX_BOUNDS_);}
        parsePropertyByTag(m_vcontent[iline],POCC_AFLOWIN_tag+"UFF_ENERGY=",(*it).m_psc_set[index_str].m_energy_uff);
        index_str++;  //while UFF_ENERGY is the only property read in, just iterate to the next structure once it is read in
      }
      ////add quick checks to see that ALL properties read in before going to the next structure
      //if(
      //    (*it).m_psc_set[index_str].m_energy_uff!=AUROSTD_MAX_DOUBLE &&
      //    TRUE
      //  ){index_str++;}
    }

    //vector<string> vtokens,vtokens2;
    //string arun_directory="",pocc_directory_abs="",qmvasp_path="";
    //m_ARUN_directories.clear();
    //
    //for(uint iline=0;iline<m_vcontent.size();iline++){
    //  if(aurostd::substring2bool(m_vcontent[iline],POCC_AFLOWIN_tag+"UFF_ENERGY=")){
    //    aurostd::string2tokens(m_vcontent[iline],vtokens,"=");
    //    if(vtokens.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown UFF_ENERGY line format",_FILE_CORRUPT_);}
    //    if(!aurostd::isfloat(vtokens[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"UFF_ENERGY is not a double",_VALUE_ERROR_);}
    //    l_supercell_sets.back().m_psc_set.push_back(POccSuperCell()); //this comes before structure
    //    l_supercell_sets.back().m_psc_set.back().m_energy_uff=aurostd::string2utype<double>(vtokens[1]);
    //    continue;
    //  }
    //  if(aurostd::substring2bool(m_vcontent[iline],POSCAR_POCC_series_START_tag)){
    //    //hnf_matrix
    //    aurostd::string2tokens(m_vcontent[iline],vtokens,"H");
    //    if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (H1) line",_FILE_CORRUPT_);}
    //    aurostd::string2tokens(vtokens.back(),vtokens2,"C");
    //    if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (H2) line",_FILE_CORRUPT_);}
    //    if(!aurostd::isfloat(vtokens2[0])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"hnf_index is not an integer",_VALUE_ILLEGAL_);}
    //    l_supercell_sets.back().m_psc_set.back().m_hnf_index=aurostd::string2utype<unsigned long long int>(vtokens2[0]);
    //    //site_configuration
    //    aurostd::string2tokens(m_vcontent[iline],vtokens,"C");
    //    if(vtokens.size()<=1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown START.POCC_ (C1) line",_FILE_CORRUPT_);}
    //    if(!aurostd::isfloat(aurostd::RemoveWhiteSpaces(vtokens.back()))){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"site_config_index is not an integer",_VALUE_ERROR_);}
    //    l_supercell_sets.back().m_psc_set.back().m_site_config_index=aurostd::string2utype<unsigned long long int>(aurostd::RemoveWhiteSpaces(vtokens.back()));
    //    //degeneracy
    //    l_supercell_sets.back().m_psc_set.back().m_degeneracy=1;
    //    if(iline+1>m_vcontent.size()-1){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot grab title at iline="+aurostd::utype2string(iline),_FILE_CORRUPT_);} //see if title is accessible
    //    if(!aurostd::substring2bool(m_vcontent[iline+1],"DG=")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"DG= is not present in POSCAR title line",_FILE_CORRUPT_);}
    //    aurostd::string2tokens(m_vcontent[iline+1],vtokens," ");
    //    for(uint i=0;i<vtokens.size();i++){
    //      if(aurostd::substring2bool(vtokens[i],"DG=")){
    //        aurostd::string2tokens(vtokens[i],vtokens2,"=");
    //        if(vtokens2.size()!=2){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot parse at DG= (vtokens2.size()=="+aurostd::utype2string(vtokens2.size())+")",_FILE_CORRUPT_);}
    //        if(!aurostd::isfloat(vtokens2[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"degeneracy is not an integer",_VALUE_ILLEGAL_);}
    //        l_supercell_sets.back().m_psc_set.back().m_degeneracy*=aurostd::string2utype<unsigned long long int>(vtokens2[1]); //backwards compatible with old pocc scheme, but this should be obsolete
    //      }
    //    }
    //    //directory
    //    pocc_directory_abs=m_vcontent[iline];
    //    aurostd::StringSubst(pocc_directory_abs,POSCAR_series_START_tag,"");
    //    aurostd::string2tokens(pocc_directory_abs,vtokens,"_");  //POCC_01_01_H0C0
    //    if(vtokens.size()!=4){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown POcc hash format (_)",_FILE_CORRUPT_);}
    //    if(!aurostd::isfloat(vtokens[2])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine POcc degeneracy index",_VALUE_ILLEGAL_);}
    //    if(aurostd::string2utype<int>(vtokens[2])==1){
    //      arun_directory="ARUN."+vtokens[0]+"_"+vtokens[1]+"_"+vtokens[3];
    //      m_ARUN_directories.push_back(arun_directory);
    //      pocc_directory_abs=m_aflags.Directory+"/"+arun_directory;
    //      if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_directory_abs=" << pocc_directory_abs << endl;}
    //      if(!aurostd::IsDirectory(pocc_directory_abs)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POcc directory [dir="+pocc_directory_abs+"] not found",_FILE_NOT_FOUND_);}
    //      //m_ARUN_directories.push_back(pocc_directory_abs); //do not save full path
    //    }
    //    //debug
    //    if(LDEBUG){
    //      cerr << __AFLOW_FUNC__ << " found hnf_index=" << l_supercell_sets.back().m_psc_set.back().m_hnf_index;
    //      cerr << " site_config_index=" << l_supercell_sets.back().m_psc_set.back().m_site_config_index;
    //      cerr << " degeneracy=" << l_supercell_sets.back().m_psc_set.back().m_degeneracy;
    //      cerr << endl;
    //    }
    //  }
    //}

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " l_supercell_sets.size()=" << l_supercell_sets.size() << endl;
      unsigned long long int isupercell=0;
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        isupercell=std::distance(l_supercell_sets.begin(),it);
        cerr << __AFLOW_FUNC__ << " l_supercell_sets[i=" << isupercell << "].m_psc_set.size()=" << (*it).m_psc_set.size() << endl;
      }
    }
  }

  bool POccStructuresFile::getARUNDirectories(vector<string>& ARUN_directories,bool tryDirectoryLS){
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);
    if(!m_initialized){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POccStructuresFile not initialized",_RUNTIME_ERROR_);}

    ARUN_directories.clear();
    uint i=0,iline=0;
    bool missing_tag=false; //very important, if we cannot find the tag for the structures_unique, we must read in structures_all to get total number of structures
    string arun_directory="",pocc_directory_abs="";
    vector<string> vtokens;
    for(i=0;i<m_vPOSCAR_lines.size()&&missing_tag==false;i++){
      if(m_vPOSCAR_lines[i].size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_vPOSCAR_lines["+aurostd::utype2string(i)+"].size()==0",_INDEX_ILLEGAL_);}
      if(m_vPOSCAR_lines[i][0].size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"m_vPOSCAR_lines["+aurostd::utype2string(i)+"][0].size()==0",_INDEX_ILLEGAL_);}
      iline=m_vPOSCAR_lines[i][0][0];
      iline-=1; //the line we are looking for is just above the POSCAR title
      if(iline>=m_vcontent.size()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"iline>=m_vcontent.size()",_INDEX_ILLEGAL_);}
      if(!aurostd::substring2bool(m_vcontent[iline],POSCAR_POCC_series_START_tag)){missing_tag=true;break;}
      //TAG FOUND: grab it
      arun_directory=m_vcontent[iline];
      aurostd::StringSubst(arun_directory,POSCAR_series_START_tag,"");
      aurostd::string2tokens(arun_directory,vtokens,"_");  //POCC_01_01_H0C0
      if(vtokens.size()==2){  //POCC_01
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " POCC_NN setting found (old POCC type)" << endl;}
        if(!aurostd::isfloat(vtokens[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine POcc structure group index",_VALUE_ILLEGAL_);}
        arun_directory="ARUN."+vtokens[0]+"_"+vtokens[1];
      }else if(vtokens.size()==3){  //POCC_01_H0C0
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " POCC_NN_HASH setting found (unique_supercells filetype)" << endl;}
        if(m_fileoptions.flag("TYPE::ALL_SUPERCELLS")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"ALL_STRUCTURES setting but detected UNIQUE_STRUCTURES format",_FILE_WRONG_FORMAT_);}
        if(!aurostd::isfloat(vtokens[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine POcc structure group index",_VALUE_ILLEGAL_);}
        arun_directory="ARUN."+vtokens[0]+"_"+vtokens[1]+"_"+vtokens[2];
      }else if(vtokens.size()==4){  //POCC_01_01_H0C0
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " POCC_NN_NN_HASH setting found (all_supercells filetype)" << endl;}
        if(m_fileoptions.flag("TYPE::UNIQUE_SUPERCELLS")){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"UNIQUE_STRUCTURES setting but detected ALL_STRUCTURES format",_FILE_WRONG_FORMAT_);}
        if(!aurostd::isfloat(vtokens[1])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine POcc structure group index",_VALUE_ILLEGAL_);}
        if(!aurostd::isfloat(vtokens[2])){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Cannot determine POcc structure index",_VALUE_ILLEGAL_);}
        arun_directory="ARUN."+vtokens[0]+"_"+vtokens[1]+"_"+vtokens[3];
      }else{
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Unknown POcc hash format (_)",_FILE_CORRUPT_);
      }
      ARUN_directories.push_back(arun_directory);
      pocc_directory_abs=m_aflags.Directory+"/"+arun_directory;
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_directory_abs=" << pocc_directory_abs << endl;}
      if(!aurostd::IsDirectory(pocc_directory_abs)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POcc directory [dir="+pocc_directory_abs+"] not found",_FILE_NOT_FOUND_);}
    }
    if(missing_tag==false && ARUN_directories.size()==l_supercell_sets.size()){return true;}

    if(m_fileoptions.flag("TYPE::ALL_SUPERCELLS")){  //build titles, but we need vstructures_size, which requires TYPE::ALL_SUPERCELLS
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " trying to build directories from information in ALL_STRUCTURES file" << endl;}
      unsigned long long int isupercell=0;
      for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
        isupercell=std::distance(l_supercell_sets.begin(),it);
        arun_directory="ARUN.POCC_"+pocc::getARUNString(l_supercell_sets,isupercell);
        ARUN_directories.push_back(arun_directory);
        pocc_directory_abs=m_aflags.Directory+"/"+arun_directory;
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " pocc_directory_abs=" << pocc_directory_abs << endl;}
        if(!aurostd::IsDirectory(pocc_directory_abs)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"POcc directory [dir="+pocc_directory_abs+"] not found",_FILE_NOT_FOUND_);}
      }
    }
    if(ARUN_directories.size()==l_supercell_sets.size()){return true;}

    //last attempt - DirectoryLS()
    if(tryDirectoryLS==false){return false;}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " trying DirectoryLS() (last attempt)" << endl;}
    aurostd::DirectoryLS(m_aflags.Directory,vtokens);
    for(uint i=0;i<vtokens.size();i++){
      if(aurostd::substring2bool(vtokens[i],"ARUN.POCC_") && aurostd::IsDirectory(m_aflags.Directory+"/"+vtokens[i])){ARUN_directories.push_back(vtokens[i]);}
    }
    std::sort(ARUN_directories.begin(),ARUN_directories.end()); //sort as DirectoryLS() may not return in order of 01 ... XX
    if(ARUN_directories.size()==l_supercell_sets.size()){return true;}

    return false;
  }

  bool POccStructuresFile::loadDataIntoCalculator(POccCalculator& pcalc,bool tryDirectoryLS) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_ || XHOST.DEBUG);

    //copy over l_supercell_sets
    pcalc.l_supercell_sets.clear();for(std::list<POccSuperCellSet>::const_iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){pcalc.l_supercell_sets.push_back(*it);}
    if(!getARUNDirectories(pcalc.m_ARUN_directories,tryDirectoryLS)){return false;}
    if(LDEBUG){
      for(uint i=0;i<pcalc.m_ARUN_directories.size();i++){cerr << __AFLOW_FUNC__ << " m_ARUN_directories[" << i << "]=" << pcalc.m_ARUN_directories[i] << endl;}
    }
    return true;
  }

} // namespace pocc

ostream& operator<<(ostream& oss, const pocc::POccStructuresFile& psf) {
  if(!psf.m_initialized){}  //keep busy
  //will do for the next round
  return oss;
}

#endif  // _AFLOW_POCC_CPP_

// ***************************************************************************
// *                                                                         *
// *              AFlow STEFANO CURTAROLO  Duke University 2003-2023         *
// *              AFlow COREY OSES  Duke University 2013-2021                *
// *                                                                         *
// ***************************************************************************
