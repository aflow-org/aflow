//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow Cormac Toher - Duke University 2013-2021           *
// *                                                                         *
//****************************************************************************
// Written by Cormac Toher, 2019-2021.
//
// This file provides a framework to calculate thermo-mechanical properties for
// disordered materials modeled using the POCC algorithm.

#ifndef _AFLOW_POCC_AEL_AGL_CPP_ //CO20200404
#define _AFLOW_POCC_AEL_AGL_CPP_ //CO20200404

#include "aflow.h"
#include "aflow_pocc.h"
#include "aflow_ael_elasticity.h"
#include "aflow_agl_debye.h"

#define _DEBUG_POCC_AEL_AGL_ false  //CO20190116

//CO20200502 START - may move to AUROSTD/aurostd_xscalar.h if these constants are needed elsewhere
//adding hardness value constants from models of Chen et al., Teter, and Tian et al.
//Chen et al. //doi:10.1016/J.INTERMET.2011.03.026
#define HARDNESS_CHEN_EXP 0.585
//Teter //doi:10.1557/S0883769400031420
#define HARDNESS_TETER_COEF 0.151
//Tian et al. //doi:10.1016/J.IJRMHM.2012.02.021
#define HARDNESS_TIAN_COEF 0.92
#define HARDNESS_TIAN_EXP_PUGHSMOD 1.137
#define HARDNESS_TIAN_EXP_SHEARMOD 0.708
//CO20200502 END - may move to AUROSTD/aurostd_xscalar.h if these constants are needed elsewhere


//CT20200319 - POCC+AEL functions
// ***************************************************************************
//  void POccCalculator::calculateElasticProperties(const vector<double>& v_temperatures)
//  Runs POCC+AEL
// ***************************************************************************
namespace pocc {
  void POccCalculator::calculateElasticProperties(const vector<double>& v_temperatures) {
    bool ael_run_postprocess = false, ael_write_full_results = false;
    vector<double> Bvoigt, Breuss, Bvrh, Gvoigt, Greuss, Gvrh, Poisson_ratio;
    vector<vector<vector<double> > > elastic_tensor_list, compliance_tensor_list;
    // for (uint i = 0; i < v_temperatures.size(); i++) {
    //  cerr << v_temperatures.at(i) << endl;
    //}
    setAELOptions(ael_run_postprocess, ael_write_full_results);
    if (ael_run_postprocess) {
      generateElasticProperties(Bvoigt, Breuss, Bvrh, Gvoigt, Greuss, Gvrh, Poisson_ratio, elastic_tensor_list, compliance_tensor_list);
    } else {
      getElasticProperties(Bvoigt, Breuss, Bvrh, Gvoigt, Greuss, Gvrh, Poisson_ratio, elastic_tensor_list, compliance_tensor_list);
    }
    getAverageElasticProperties(v_temperatures, ael_write_full_results, Bvoigt, Breuss, Bvrh, Gvoigt, Greuss, Gvrh, Poisson_ratio, elastic_tensor_list, compliance_tensor_list);
  } 
} // namespace pocc

// ***************************************************************************
// void POccCalculator::setAELOptions(bool& ael_run_postprocess, bool& ael_write_full_results)
// Sets options for running POCC+AEL from command line
// ael_run_postprocess: runs AEL postprocessing; otherwise reads results from aflow.ael.out
// ael_write_full_results: writes out results for all POCC configurations
// ***************************************************************************
namespace pocc {
  void POccCalculator::setAELOptions(bool& ael_run_postprocess, bool& ael_write_full_results) {
    // bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    ael_run_postprocess = false;
    ael_write_full_results = false;

    //run postprocessing
    ael_run_postprocess=(ael_run_postprocess || XHOST.vflag_control.flag("AEL_RUN_POSTPROCESSING") ); //CO20201026
    ael_run_postprocess=(ael_run_postprocess || XHOST.POSTPROCESS ); //CO20201026

    //write extra output (debugging mostly)
    ael_write_full_results=(ael_write_full_results || XHOST.vflag_control.flag("AEL_WRITE_FULL_RESULTS") ); //CO20201026
  }
} // namespace pocc

// ***************************************************************************
//  void POccCalculator::generateElasticProperties(vector<double>& Bvoigt, vector<double>& Breuss, vector<double>& Bvrh, vector<double>& Gvoigt, vector<double>& Greuss, vector<double>& Gvrh,vector<double>& Poisson_ratio, vector<vector<vector<double> > >& elastic_tensor_list,  vector<vector<vector<double> > >& compliance_tensor_list)
//  Runs AEL postprocessing to generate elastic properties 
// ***************************************************************************
namespace pocc {
  void POccCalculator::generateElasticProperties(vector<double>& Bvoigt, vector<double>& Breuss, vector<double>& Bvrh, vector<double>& Gvoigt, vector<double>& Greuss, vector<double>& Gvrh,vector<double>& Poisson_ratio, vector<vector<vector<double> > >& elastic_tensor_list,  vector<vector<vector<double> > >& compliance_tensor_list) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    uint aelerror = 0;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " directory=[" << m_aflags.Directory << "]" << endl;
    }
    double ael_bulk_modulus_voigt = 0.0, ael_bulk_modulus_reuss = 0.0, ael_bulk_modulus_vrh = 0.0;
    double ael_shear_modulus_voigt = 0.0, ael_shear_modulus_reuss = 0.0, ael_shear_modulus_vrh = 0.0;
    double ael_poisson_ratio = 0.0;
    double elastic_const = 0.0;
    vector<double> elastic_const_row;
    vector<vector<double> > elastic_tensor;
    double compliance_const=0.0;
    vector<double> compliance_const_row;
    vector<vector<double> > compliance_tensor;
    unsigned long long int isupercell = 0;
    string pocc_directory_abs="";
    string ael_filename="";
    stringstream aflow_ael_out;
    stringstream aflow_elastic_tensor;
    stringstream aflow_compliance_tensor;
    vector<string> vline, tokens, vfiles;
    bool run_directory=false;
    bool ael_aflowin_found = false;
    string AflowInName = _AFLOWIN_;
    string FileLockName = _AFLOWLOCK_;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      pocc_directory_abs=m_aflags.Directory+"/"+m_ARUN_directories[isupercell];
      aurostd::StringstreamClean(aflow_ael_out);
      aurostd::StringstreamClean(aflow_elastic_tensor);
      aurostd::StringstreamClean(aflow_compliance_tensor);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " look at pocc_directory_abs=" << pocc_directory_abs << endl;}
      //energy_dft
      ael_filename=pocc_directory_abs+"/aflow.ael.out";
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ael_filename=[" << ael_filename << "]" << endl;}
      // if(!aurostd::EFileExist(ael_filename,ael_filename)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No aflow.ael.out file found [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);}
      // [OBSOLETE] aelerror = AEL_functions::Get_ElasticProperties_AEL_postprocess(pocc_directory_abs, ael_bulk_modulus_voigt, ael_bulk_modulus_reuss, ael_bulk_modulus_vrh, ael_shear_modulus_voigt, ael_shear_modulus_reuss, ael_shear_modulus_vrh, ael_poisson_ratio, elastic_tensor, compliance_tensor);
      // [OBSOLETE] if (aelerror > 0) {
      // [OBSOLETE]   throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem with AEL calculation: [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);
      // [OBSOLETE] }
      aelerror = AEL_functions::AEL_Get_AflowInName(AflowInName, pocc_directory_abs, ael_aflowin_found);
      if (aelerror != 0) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem finding AEL aflow.in filename [dir="+pocc_directory_abs+"]",_FILE_ERROR_);
      }
      if (ael_aflowin_found) {
        if(aurostd::FileExist(pocc_directory_abs+"/ael.LOCK")) {
          FileLockName = "ael.LOCK";
        } else if(aurostd::FileExist(pocc_directory_abs+"/agl.LOCK")) {
          FileLockName = "agl.LOCK";
        }
        run_directory=true;
      }
      if(run_directory){
        _aflags aflags;
        aflags.Directory=pocc_directory_abs;

        //save originals
        string _AFLOWIN_orig=_AFLOWIN_;
        string _AFLOWLOCK_orig=_AFLOWLOCK_;

        //set env for RUN_Directory()
        _AFLOWIN_=AflowInName;
        _AFLOWLOCK_=FileLockName;
        if(aurostd::FileExist(pocc_directory_abs+"/"+_AFLOWLOCK_)){aurostd::file2file(pocc_directory_abs+"/"+_AFLOWLOCK_,pocc_directory_abs+"/"+_AFLOWLOCK_+".run");} //keep original LOCK
        KBIN::RUN_Directory(aflags);

        //return to original
        _AFLOWIN_=_AFLOWIN_orig;
        _AFLOWLOCK_=_AFLOWLOCK_orig;
      }
      if(aurostd::EFileExist(pocc_directory_abs+"/aflow.ael.out")) {
        aflow_ael_out.str(aurostd::substring2string(aurostd::efile2string(pocc_directory_abs+"/"+"aflow.ael.out"),"[AEL_RESULTS]START","[AEL_RESULTS]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::StringstreamClean(aflow_ael_out);
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(pocc_directory_abs+"/"+"aflow.ael.out"),aflow_ael_out,"[AEL_RESULTS]START","[AEL_RESULTS]STOP");
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " aflow.ael.out:" << endl << aflow_ael_out.str() << endl;
        }
        aurostd::stream2vectorstring(aflow_ael_out,vline);
        for (uint i=0;i<vline.size();i++) {
          aurostd::StringSubst(vline.at(i),"="," ");
          aurostd::string2tokens(vline.at(i),tokens," ");
          if(tokens.size()>=2) {
            if(tokens.at(0)=="ael_poisson_ratio") ael_poisson_ratio=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="ael_bulk_modulus_voigt") ael_bulk_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="ael_bulk_modulus_reuss") ael_bulk_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="ael_shear_modulus_voigt") ael_shear_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="ael_shear_modulus_reuss") ael_shear_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="ael_bulk_modulus_vrh") ael_bulk_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="ael_shear_modulus_vrh") ael_shear_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
          }
        }
        aflow_elastic_tensor.str(aurostd::substring2string(aurostd::efile2string(ael_filename),"[AEL_STIFFNESS_TENSOR]START","[AEL_STIFFNESS_TENSOR]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(ael_filename),aflow_elastic_tensor,"[AEL_STIFFNESS_TENSOR]START","[AEL_STIFFNESS_TENSOR]STOP");
        if(LDEBUG) {
          //aurostd::StringstreamClean(aus);
          cerr << __AFLOW_FUNC__ << " aflow_elastic_tensor = " << endl << aflow_elastic_tensor.str() << endl;
          //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        aurostd::stream2vectorstring(aflow_elastic_tensor,vline);
        for (uint i=0;i<vline.size();i++) {
          aurostd::string2tokens(vline.at(i),tokens," ");
          if(tokens.size()>=6) {
            elastic_const_row.clear();
            for(uint j=0; j<tokens.size(); j++) {
              elastic_const=aurostd::string2utype<double>(tokens.at(j));
              elastic_const_row.push_back(elastic_const);
            }
            elastic_tensor.push_back(elastic_const_row);
          }
        }
        aflow_compliance_tensor.str(aurostd::substring2string(aurostd::efile2string(ael_filename),"[AEL_COMPLIANCE_TENSOR]START","[AEL_COMPLIANCE_TENSOR]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(ael_filename),aflow_compliance_tensor,"[AEL_COMPLIANCE_TENSOR]START","[AEL_COMPLIANCE_TENSOR]STOP");
        if(LDEBUG) {
          //aurostd::StringstreamClean(aus);
          cerr << __AFLOW_FUNC__ << " aflow_compliance_tensor = " << endl << aflow_compliance_tensor.str() << endl;
          //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        aurostd::stream2vectorstring(aflow_compliance_tensor,vline);
        for (uint i=0;i<vline.size();i++) {
          aurostd::string2tokens(vline.at(i),tokens," ");
          if(tokens.size()>=6) {
            compliance_const_row.clear();
            for(uint j=0; j<tokens.size(); j++) {
              compliance_const=aurostd::string2utype<double>(tokens.at(j));
              compliance_const_row.push_back(compliance_const);
            }
            compliance_tensor.push_back(compliance_const_row);
          }
        }      
        Bvoigt.push_back(ael_bulk_modulus_voigt);
        Breuss.push_back(ael_bulk_modulus_reuss);
        Bvrh.push_back(ael_bulk_modulus_vrh);
        Gvoigt.push_back(ael_shear_modulus_voigt);
        Greuss.push_back(ael_shear_modulus_reuss);
        Gvrh.push_back(ael_shear_modulus_vrh);
        Poisson_ratio.push_back(ael_poisson_ratio);
        elastic_tensor_list.push_back(elastic_tensor);
        compliance_tensor_list.push_back(compliance_tensor);
        elastic_tensor.clear();
        compliance_tensor.clear();
      } else {
        // [OBSOLETE] cerr << __AFLOW_FUNC__ << pocc_directory_abs << "/aflow.ael.out: File not found" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem with AEL calculation: [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);
      }
    }
  }
} // namespace pocc

// ***************************************************************************
//  void POccCalculator::getElasticProperties(vector<double>& Bvoigt, vector<double>& Breuss, vector<double>& Bvrh, vector<double>& Gvoigt, vector<double>& Greuss, vector<double>& Gvrh,vector<double>& Poisson_ratio, vector<vector<vector<double> > >& elastic_tensor_list,  vector<vector<vector<double> > >& compliance_tensor_list)
//  Retrieves elastic properties from aflow.ael.out files 
// ***************************************************************************
namespace pocc {
  void POccCalculator::getElasticProperties(vector<double>& Bvoigt, vector<double>& Breuss, vector<double>& Bvrh, vector<double>& Gvoigt, vector<double>& Greuss, vector<double>& Gvrh,vector<double>& Poisson_ratio, vector<vector<vector<double> > >& elastic_tensor_list,  vector<vector<vector<double> > >& compliance_tensor_list) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    stringstream message;
    //uint aelerror = 0;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " directory=[" << m_aflags.Directory << "]" << endl;
    }
    double ael_bulk_modulus_voigt = 0.0, ael_bulk_modulus_reuss = 0.0, ael_bulk_modulus_vrh = 0.0;
    double ael_shear_modulus_voigt = 0.0, ael_shear_modulus_reuss = 0.0, ael_shear_modulus_vrh = 0.0;
    double ael_poisson_ratio = 0.0;
    double elastic_const = 0.0;
    vector<double> elastic_const_row;
    vector<vector<double> > elastic_tensor;
    double compliance_const = 0.0;
    vector<double> compliance_const_row;
    vector<vector<double> > compliance_tensor;    
    unsigned long long int isupercell = 0;
    string pocc_directory_abs="";
    string ael_filename="";
    stringstream aflow_ael_out;
    stringstream aflow_elastic_tensor;
    stringstream aflow_compliance_tensor;
    vector<string> vline, tokens, vfiles;
    bool aelfilezipped = false;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      pocc_directory_abs=m_aflags.Directory+"/"+m_ARUN_directories[isupercell];
      aurostd::StringstreamClean(aflow_ael_out);
      aurostd::StringstreamClean(aflow_elastic_tensor);
      aurostd::StringstreamClean(aflow_compliance_tensor);
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " look at pocc_directory_abs=" << pocc_directory_abs << endl;}
      //energy_dft
      ael_filename=pocc_directory_abs+"/aflow.ael.out";
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " ael_filename=[" << ael_filename << "]" << endl;}
      // if(!aurostd::EFileExist(ael_filename,ael_filename)){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No aflow.ael.out file found [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);}
      if(aurostd::EFileExist(ael_filename,ael_filename)){
        message << "Reading " << m_ARUN_directories[isupercell]+"/aflow.ael.out";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); //CO20200404
        aflow_ael_out.str(aurostd::substring2string(aurostd::efile2string(ael_filename),"[AEL_RESULTS]START","[AEL_RESULTS]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(ael_filename),aflow_ael_out,"[AEL_RESULTS]START","[AEL_RESULTS]STOP");
        aelfilezipped = true;
        //[CO20200404 - OBSOLETE]} else if (aurostd::FileExist(ael_filename))
        //[CO20200404 - OBSOLETE]  aurostd::ExtractToStringstreamEXPLICIT(aurostd::file2string(ael_filename),aflow_ael_out,"[AEL_RESULTS]START","[AEL_RESULTS]STOP");
        //[CO20200404 - OBSOLETE]  aelfilezipped = false;
    } else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No aflow.ael.out file found [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);
    }
    if(LDEBUG) {
      //aurostd::StringstreamClean(aus);
      cerr << __AFLOW_FUNC__ << " aflow_ael_out:" << endl << aflow_ael_out.str() << endl;
      //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    aurostd::stream2vectorstring(aflow_ael_out,vline);
    for (uint i=0;i<vline.size();i++) {
      aurostd::StringSubst(vline.at(i),"="," ");
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(tokens.size()>=2) {
        if(tokens.at(0)=="ael_poisson_ratio") ael_poisson_ratio=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="ael_bulk_modulus_voigt") ael_bulk_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="ael_bulk_modulus_reuss") ael_bulk_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="ael_shear_modulus_voigt") ael_shear_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="ael_shear_modulus_reuss") ael_shear_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="ael_bulk_modulus_vrh") ael_bulk_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="ael_shear_modulus_vrh") ael_shear_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
      }
    }
    if (aelfilezipped) {
      aflow_elastic_tensor.str(aurostd::substring2string(aurostd::efile2string(ael_filename),"[AEL_STIFFNESS_TENSOR]START","[AEL_STIFFNESS_TENSOR]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(ael_filename),aflow_elastic_tensor,"[AEL_STIFFNESS_TENSOR]START","[AEL_STIFFNESS_TENSOR]STOP");
    } else {
      aflow_elastic_tensor.str(aurostd::substring2string(aurostd::file2string(ael_filename),"[AEL_STIFFNESS_TENSOR]START","[AEL_STIFFNESS_TENSOR]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::file2string(ael_filename),aflow_elastic_tensor,"[AEL_STIFFNESS_TENSOR]START","[AEL_STIFFNESS_TENSOR]STOP");
    }
    if(LDEBUG) {
      //aurostd::StringstreamClean(aus);
      cerr << __AFLOW_FUNC__ << " aflow_elastic_tensor = " << endl << aflow_elastic_tensor.str() << endl;
      //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    aurostd::stream2vectorstring(aflow_elastic_tensor,vline);
    for (uint i=0;i<vline.size();i++) {
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(tokens.size()>=6) {
        elastic_const_row.clear();
        for(uint j=0; j<tokens.size(); j++) {
          elastic_const=aurostd::string2utype<double>(tokens.at(j));
          elastic_const_row.push_back(elastic_const);
        }
        elastic_tensor.push_back(elastic_const_row);
      }
    }
    if (aelfilezipped) {
      aflow_compliance_tensor.str(aurostd::substring2string(aurostd::efile2string(ael_filename),"[AEL_COMPLIANCE_TENSOR]START","[AEL_COMPLIANCE_TENSOR]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(ael_filename),aflow_compliance_tensor,"[AEL_COMPLIANCE_TENSOR]START","[AEL_COMPLIANCE_TENSOR]STOP");
    } else {
      aflow_compliance_tensor.str(aurostd::substring2string(aurostd::file2string(ael_filename),"[AEL_COMPLIANCE_TENSOR]START","[AEL_COMPLIANCE_TENSOR]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::file2string(ael_filename),aflow_compliance_tensor,"[AEL_COMPLIANCE_TENSOR]START","[AEL_COMPLIANCE_TENSOR]STOP");
    }
    if(LDEBUG) {
      //aurostd::StringstreamClean(aus);
      cerr << __AFLOW_FUNC__ << " aflow_compliance_tensor = " << endl << aflow_compliance_tensor.str() << endl;
      //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    aurostd::stream2vectorstring(aflow_compliance_tensor,vline);
    for (uint i=0;i<vline.size();i++) {
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(tokens.size()>=6) {
        compliance_const_row.clear();
        for(uint j=0; j<tokens.size(); j++) {
          compliance_const=aurostd::string2utype<double>(tokens.at(j));
          compliance_const_row.push_back(compliance_const);
        }
        compliance_tensor.push_back(compliance_const_row);
      }
    }
    Bvoigt.push_back(ael_bulk_modulus_voigt);
    Breuss.push_back(ael_bulk_modulus_reuss);
    Bvrh.push_back(ael_bulk_modulus_vrh);
    Gvoigt.push_back(ael_shear_modulus_voigt);
    Greuss.push_back(ael_shear_modulus_reuss);
    Gvrh.push_back(ael_shear_modulus_vrh);
    Poisson_ratio.push_back(ael_poisson_ratio);
    elastic_tensor_list.push_back(elastic_tensor);
    compliance_tensor_list.push_back(compliance_tensor);
    aflow_ael_out.clear();
    aflow_elastic_tensor.clear();
    aflow_compliance_tensor.clear();
    elastic_tensor.clear();
    compliance_tensor.clear();
  }
}
} // namespace pocc

// ***************************************************************************
//   void POccCalculator::getAverageElasticProperties(const vector<double>& v_temperatures, bool& ael_write_full_results, vector<double>& Bvoigt, vector<double>& Breuss, vector<double>& Bvrh, vector<double>& Gvoigt, vector<double>& Greuss, vector<double>& Gvrh,vector<double>& Poisson_ratio, vector<vector<vector<double> > >& elastic_tensor_list,  vector<vector<vector<double> > >& compliance_tensor_list) 
// Calculates thermal average of elastic properties for all POCC configurations
// Writes results to aflow.pocc_ael.out
// ***************************************************************************
namespace pocc {
  void POccCalculator::getAverageElasticProperties(const vector<double>& v_temperatures, bool ael_write_full_results, vector<double>& Bvoigt, vector<double>& Breuss, vector<double>& Bvrh, vector<double>& Gvoigt, vector<double>& Greuss, vector<double>& Gvrh,vector<double>& Poisson_ratio, vector<vector<vector<double> > >& elastic_tensor_list,  vector<vector<vector<double> > >& compliance_tensor_list) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    stringstream message;
    ostringstream aus;
    ofstream FileMESSAGE;
    double T = 0.0;
    vector<double> Yvoigt, Yreuss, Yvrh;
    vector<double> Hv_Teter_voigt, Hv_Teter_reuss, Hv_Teter_vrh;
    vector<double> Hv_Chen_voigt, Hv_Chen_reuss, Hv_Chen_vrh;
    vector<double> Hv_Tian_voigt, Hv_Tian_reuss, Hv_Tian_vrh;
    vector<double> Pugh_ratio, Cauchy_pressure;
    vector<vector<double> > elastic_tensor_ave, compliance_tensor_ave;
    vector<double> elastic_tensor_row, compliance_tensor_row;
    double Bvoigt_ave=0.0, Breuss_ave=0.0, Bvrh_ave=0.0, Gvoigt_ave=0.0, Greuss_ave=0.0, Gvrh_ave=0.0;
    double Hv_Teter_voigt_val=0.0, Hv_Teter_reuss_val=0.0, Hv_Teter_vrh_val=0.0;
    double Hv_Chen_voigt_val=0.0, Hv_Chen_reuss_val=0.0, Hv_Chen_vrh_val=0.0;
    double Hv_Tian_voigt_val=0.0, Hv_Tian_reuss_val=0.0, Hv_Tian_vrh_val=0.0;
    double Hv_Teter_voigt_ave=0.0, Hv_Teter_reuss_ave=0.0, Hv_Teter_vrh_ave=0.0;
    double Hv_Chen_voigt_ave=0.0, Hv_Chen_reuss_ave=0.0, Hv_Chen_vrh_ave=0.0;
    double Hv_Tian_voigt_ave=0.0, Hv_Tian_reuss_ave=0.0, Hv_Tian_vrh_ave=0.0;
    double Yvoigt_ave=0.0, Yreuss_ave=0.0, Yvrh_ave=0.0, Yvoigt_val=0.0, Yreuss_val=0.0, Yvrh_val=0.0;
    double Bvoigt_max=0.0, Breuss_max=0.0, Bvrh_max=0.0, Gvoigt_max=0.0, Greuss_max=0.0, Gvrh_max=0.0, Yvoigt_max=0.0, Yreuss_max=0.0, Yvrh_max=0.0;
    double Bvoigt_min=0.0, Breuss_min=0.0, Bvrh_min=0.0, Gvoigt_min=0.0, Greuss_min=0.0, Gvrh_min=0.0, Yvoigt_min=0.0, Yreuss_min=0.0, Yvrh_min=0.0;
    double Hv_Teter_voigt_max=0.0, Hv_Teter_reuss_max=0.0, Hv_Teter_vrh_max=0.0, Hv_Chen_voigt_max=0.0, Hv_Chen_reuss_max=0.0, Hv_Chen_vrh_max=0.0, Hv_Tian_voigt_max=0.0, Hv_Tian_reuss_max=0.0, Hv_Tian_vrh_max=0.0;
    double Hv_Teter_voigt_min=0.0, Hv_Teter_reuss_min=0.0, Hv_Teter_vrh_min=0.0, Hv_Chen_voigt_min=0.0, Hv_Chen_reuss_min=0.0, Hv_Chen_vrh_min=0.0, Hv_Tian_voigt_min=0.0, Hv_Tian_reuss_min=0.0, Hv_Tian_vrh_min=0.0;
    double Poisson_ratio_ave=0.0, Poisson_ratio_min=0.0, Poisson_ratio_max=0.0;
    double Pugh_val=0.0, Pugh_ratio_ave=0.0, Pugh_ratio_min=0.0, Pugh_ratio_max=0.0;
    double Cauchy_ave=0.0, Cauchy_val=0.0, Cauchy_min=0.0, Cauchy_max=0.0;
    xmatrix<double> elastic_const(6, 6), elastic_const_inv(6, 6);

    double kval = 0.0, kvalsq = 0.0;
    stringstream ossmain, ossfull;

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " v_temperatures=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(v_temperatures,5),",") << endl;}  //CO20200404

    // Calculate elastic modulus, hardness and ductility 
    unsigned long long int isupercell = 0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);
      // Young's modulus
      Yvoigt_val = (9.0 * Bvoigt.at(isupercell) * Gvoigt.at(isupercell)) / (3.0 * Bvoigt.at(isupercell) + Gvoigt.at(isupercell));
      Yreuss_val = (9.0 * Breuss.at(isupercell) * Greuss.at(isupercell)) / (3.0 * Breuss.at(isupercell) + Greuss.at(isupercell));
      Yvrh_val = (9.0 * Bvrh.at(isupercell) * Gvrh.at(isupercell)) / (3.0 * Bvrh.at(isupercell) + Gvrh.at(isupercell));
      Yvoigt.push_back(Yvoigt_val);
      Yreuss.push_back(Yreuss_val);
      Yvrh.push_back(Yvrh_val);
      // Pugh's modulus ratio
      Pugh_val = Gvrh.at(isupercell) / Bvrh.at(isupercell);
      Pugh_ratio.push_back(Pugh_val);
      // Cauchy pressure
      Cauchy_val = ((elastic_tensor_list.at(isupercell).at(0).at(1) + elastic_tensor_list.at(isupercell).at(0).at(2) + elastic_tensor_list.at(isupercell).at(1).at(2)) / 3.0) - ((elastic_tensor_list.at(isupercell).at(3).at(3) + elastic_tensor_list.at(isupercell).at(4).at(4) + elastic_tensor_list.at(isupercell).at(5).at(5)) / 3.0);
      Cauchy_pressure.push_back(Cauchy_val);
      // Hardness: Teter model
      Hv_Teter_voigt_val = HARDNESS_TETER_COEF * Gvoigt.at(isupercell);
      Hv_Teter_reuss_val = HARDNESS_TETER_COEF * Greuss.at(isupercell);      
      Hv_Teter_vrh_val = HARDNESS_TETER_COEF * Gvrh.at(isupercell);
      Hv_Teter_voigt.push_back(Hv_Teter_voigt_val);
      Hv_Teter_reuss.push_back(Hv_Teter_reuss_val);
      Hv_Teter_vrh.push_back(Hv_Teter_vrh_val);
      // Hardness: Chen and Voigt models
      kval = Gvoigt.at(isupercell) / Bvoigt.at(isupercell);
      kvalsq = kval * kval;
      Hv_Chen_voigt_val = (2.0 * pow(kvalsq * Gvoigt.at(isupercell), HARDNESS_CHEN_EXP)) - 3.0;
      Hv_Tian_voigt_val = HARDNESS_TIAN_COEF * pow(kval, HARDNESS_TIAN_EXP_PUGHSMOD) * pow(Gvoigt.at(isupercell), HARDNESS_TIAN_EXP_SHEARMOD);
      kval = Greuss.at(isupercell) / Breuss.at(isupercell);
      kvalsq = kval * kval;      
      Hv_Chen_reuss_val = (2.0 * pow(kvalsq * Greuss.at(isupercell), HARDNESS_CHEN_EXP)) - 3.0;
      Hv_Tian_reuss_val = HARDNESS_TIAN_COEF * pow(kval, HARDNESS_TIAN_EXP_PUGHSMOD) * pow(Greuss.at(isupercell), HARDNESS_TIAN_EXP_SHEARMOD);      
      kval = Gvrh.at(isupercell) / Bvrh.at(isupercell);
      kvalsq = kval * kval;      
      Hv_Chen_vrh_val = (2.0 * pow(kvalsq * Gvrh.at(isupercell), HARDNESS_CHEN_EXP)) - 3.0;
      Hv_Tian_vrh_val = HARDNESS_TIAN_COEF * pow(kval, HARDNESS_TIAN_EXP_PUGHSMOD) * pow(Gvrh.at(isupercell), HARDNESS_TIAN_EXP_SHEARMOD);
      Hv_Chen_voigt.push_back(Hv_Chen_voigt_val);
      Hv_Chen_reuss.push_back(Hv_Chen_reuss_val);
      Hv_Chen_vrh.push_back(Hv_Chen_vrh_val);
      Hv_Tian_voigt.push_back(Hv_Tian_voigt_val);
      Hv_Tian_reuss.push_back(Hv_Tian_reuss_val);
      Hv_Tian_vrh.push_back(Hv_Tian_vrh_val);
    }


    // Find maximum and minimum values of elastic properties
    Bvoigt_min = aurostd::min(Bvoigt);
    Bvoigt_max = aurostd::max(Bvoigt); 
    Breuss_min = aurostd::min(Breuss);
    Breuss_max = aurostd::max(Breuss);
    Bvrh_min = aurostd::min(Bvrh);
    Bvrh_max = aurostd::max(Bvrh); 
    Gvoigt_min = aurostd::min(Gvoigt);
    Gvoigt_max = aurostd::max(Gvoigt); 
    Greuss_min = aurostd::min(Greuss);
    Greuss_max = aurostd::max(Greuss);
    Gvrh_min = aurostd::min(Gvrh);
    Gvrh_max = aurostd::max(Gvrh); 
    Yvoigt_min = aurostd::min(Yvoigt);
    Yvoigt_max = aurostd::max(Yvoigt); 
    Yreuss_min = aurostd::min(Yreuss);
    Yreuss_max = aurostd::max(Yreuss);
    Yvrh_min = aurostd::min(Yvrh);
    Yvrh_max = aurostd::max(Yvrh); 
    Poisson_ratio_min = aurostd::min(Poisson_ratio);
    Poisson_ratio_max = aurostd::max(Poisson_ratio);    
    Hv_Teter_voigt_min = aurostd::min(Hv_Teter_voigt);
    Hv_Teter_voigt_max = aurostd::max(Hv_Teter_voigt);
    Hv_Teter_reuss_min = aurostd::min(Hv_Teter_reuss);
    Hv_Teter_reuss_max = aurostd::max(Hv_Teter_reuss);
    Hv_Teter_vrh_min = aurostd::min(Hv_Teter_vrh);
    Hv_Teter_vrh_max = aurostd::max(Hv_Teter_vrh);
    Hv_Chen_voigt_min = aurostd::min(Hv_Chen_voigt);
    Hv_Chen_voigt_max = aurostd::max(Hv_Chen_voigt);
    Hv_Chen_reuss_min = aurostd::min(Hv_Chen_reuss);
    Hv_Chen_reuss_max = aurostd::max(Hv_Chen_reuss);
    Hv_Chen_vrh_min = aurostd::min(Hv_Chen_vrh);
    Hv_Chen_vrh_max = aurostd::max(Hv_Chen_vrh);
    Hv_Tian_voigt_min = aurostd::min(Hv_Tian_voigt);
    Hv_Tian_voigt_max = aurostd::max(Hv_Tian_voigt);
    Hv_Tian_reuss_min = aurostd::min(Hv_Tian_reuss);
    Hv_Tian_reuss_max = aurostd::max(Hv_Tian_reuss);
    Hv_Tian_vrh_min = aurostd::min(Hv_Tian_vrh);
    Hv_Tian_vrh_max = aurostd::max(Hv_Tian_vrh);
    Cauchy_min = aurostd::min(Cauchy_pressure);
    Cauchy_max = aurostd::max(Cauchy_pressure);
    Pugh_ratio_min = aurostd::min(Pugh_ratio);
    Pugh_ratio_max = aurostd::max(Pugh_ratio);    

    // Boltzmann average properties
    ossmain.clear();
    ossmain.str(std::string());
    ossfull.clear();
    ossfull.str(std::string());
    for (uint j=0; j<elastic_tensor_list.at(0).size(); j++) {
      elastic_tensor_row.clear();
      for (uint k=0; k<elastic_tensor_list.at(0).at(j).size(); k++) {
        elastic_tensor_row.push_back(0.0);
      }
      elastic_tensor_ave.push_back(elastic_tensor_row);
    }

    for (uint j=0; j<compliance_tensor_list.at(0).size(); j++) {
      compliance_tensor_row.clear();
      for (uint k=0; k<compliance_tensor_list.at(0).at(j).size(); k++) {
        compliance_tensor_row.push_back(0.0);
      }
      compliance_tensor_ave.push_back(compliance_tensor_row);
    }

    for (uint jt = 0; jt < v_temperatures.size(); jt++) {
      T = v_temperatures.at(jt);
      message << "Averaging elastic properties at T=" << T << "K";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); //CO20200404 //CO20200502 - no getTemperatureString(T) needed here
      setPOccStructureProbabilities(T);
      Bvoigt_ave = 0.0;
      Breuss_ave = 0.0;
      Bvrh_ave = 0.0;
      Gvoigt_ave = 0.0;
      Greuss_ave = 0.0;
      Gvrh_ave = 0.0;
      Yvoigt_ave = 0.0;
      Yreuss_ave = 0.0;
      Yvrh_ave = 0.0;
      Poisson_ratio_ave = 0.0;
      Pugh_ratio_ave = 0.0;
      Cauchy_ave = 0.0;
      Hv_Teter_voigt_ave = 0.0;
      Hv_Teter_reuss_ave = 0.0;
      Hv_Teter_vrh_ave = 0.0;
      Hv_Chen_voigt_ave = 0.0;
      Hv_Chen_reuss_ave = 0.0;
      Hv_Chen_vrh_ave = 0.0;
      Hv_Tian_voigt_ave = 0.0;
      Hv_Tian_reuss_ave = 0.0;
      Hv_Tian_vrh_ave = 0.0;
      for (uint j=0; j<elastic_tensor_ave.size(); j++) {
        for (uint k=0; k<elastic_tensor_ave.at(j).size(); k++) {
          elastic_tensor_ave.at(j).at(k) = 0.0;
        }
      }

      for (uint j=0; j<compliance_tensor_ave.size(); j++) {
        for (uint k=0; k<compliance_tensor_ave.at(j).size(); k++) {
          compliance_tensor_ave.at(j).at(k) = 0.0;
        }
      }

      isupercell = 0;
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
        isupercell = std::distance(l_supercell_sets.begin(), it);
        // Boltzmann weighted average for bulk modulus
        Bvoigt_ave = Bvoigt_ave + (*it).m_probability * Bvoigt.at(isupercell);
        Breuss_ave = Breuss_ave + (*it).m_probability * Breuss.at(isupercell);
        Bvrh_ave = Bvrh_ave + (*it).m_probability * Bvrh.at(isupercell);
        // Boltzmann weighted average for shear modulus      
        Gvoigt_ave = Gvoigt_ave + (*it).m_probability * Gvoigt.at(isupercell);
        Greuss_ave = Greuss_ave + (*it).m_probability * Greuss.at(isupercell);
        Gvrh_ave = Gvrh_ave + (*it).m_probability * Gvrh.at(isupercell);
        // Boltzmann weighted average for Young's modulus      
        Yvoigt_ave = Yvoigt_ave + (*it).m_probability * Yvoigt.at(isupercell);
        Yreuss_ave = Yreuss_ave + (*it).m_probability * Yreuss.at(isupercell);
        Yvrh_ave = Yvrh_ave + (*it).m_probability * Yvrh.at(isupercell);
        // Boltzmann weighted average for Poisson ratio
        Poisson_ratio_ave = Poisson_ratio_ave + (*it).m_probability * Poisson_ratio.at(isupercell);
        // Boltzmann weighted average for Pugh's modulus ratio
        Pugh_ratio_ave = Pugh_ratio_ave + (*it).m_probability * Pugh_ratio.at(isupercell);
        // Boltzmann weighted average for Cauchy pressure
        Cauchy_ave = Cauchy_ave + (*it).m_probability * Cauchy_pressure.at(isupercell);
        // Boltzmann weighted average for hardness    
        Hv_Teter_voigt_ave = Hv_Teter_voigt_ave + (*it).m_probability * Hv_Teter_voigt_val;
        Hv_Teter_reuss_ave = Hv_Teter_reuss_ave + (*it).m_probability * Hv_Teter_reuss_val;
        Hv_Teter_vrh_ave = Hv_Teter_vrh_ave + (*it).m_probability * Hv_Teter_vrh_val;
        Hv_Chen_voigt_ave = Hv_Chen_voigt_ave + (*it).m_probability * Hv_Chen_voigt_val;
        Hv_Chen_reuss_ave = Hv_Chen_reuss_ave + (*it).m_probability * Hv_Chen_reuss_val;
        Hv_Chen_vrh_ave = Hv_Chen_vrh_ave + (*it).m_probability * Hv_Chen_vrh_val;
        Hv_Tian_voigt_ave = Hv_Tian_voigt_ave + (*it).m_probability * Hv_Tian_voigt_val;
        Hv_Tian_reuss_ave = Hv_Tian_reuss_ave + (*it).m_probability * Hv_Tian_reuss_val;
        Hv_Tian_vrh_ave = Hv_Tian_vrh_ave + (*it).m_probability * Hv_Tian_vrh_val;
        for (uint j=0; j<elastic_tensor_list.at(0).size(); j++) {
          for (uint k=0; k<elastic_tensor_list.at(0).at(j).size(); k++) {
            elastic_tensor_ave.at(j).at(k) = elastic_tensor_ave.at(j).at(k) + ((*it).m_probability * elastic_tensor_list.at(isupercell).at(j).at(k));
          }
        }
        for (uint j=0; j<compliance_tensor_list.at(0).size(); j++) {
          for (uint k=0; k<compliance_tensor_list.at(0).at(j).size(); k++) {
            compliance_tensor_ave.at(j).at(k) = compliance_tensor_ave.at(j).at(k) + ((*it).m_probability * compliance_tensor_list.at(isupercell).at(j).at(k));
          }
        }
      }
      // Write main results to a file
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AEL_RESULTS]START" << endl;
      ossmain << "[POCC_AEL_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << "pocc_ael_temperature=" << T << "  (K)" << endl; //CO20200502 - no (*this).getTemperatureString(T) needed here
      ossmain << "pocc_ael_poisson_ratio=" << Poisson_ratio_ave << endl;    
      ossmain << "pocc_ael_bulk_modulus_voigt=" << Bvoigt_ave << "  (GPa)" << endl;
      ossmain << "pocc_ael_bulk_modulus_reuss=" << Breuss_ave  << "  (GPa)" << endl;
      ossmain << "pocc_ael_shear_modulus_voigt=" << Gvoigt_ave << "  (GPa)" << endl;
      ossmain << "pocc_ael_shear_modulus_reuss=" << Greuss_ave << "  (GPa)" << endl;
      ossmain << "pocc_ael_bulk_modulus_vrh=" << Bvrh_ave << "  (GPa)" << endl;    
      ossmain << "pocc_ael_shear_modulus_vrh=" << Gvrh_ave << "  (GPa)" << endl;
      ossmain << "pocc_ael_youngs_modulus_vrh=" << Yvrh_ave << "  (GPa)" << endl;
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AEL_RESULTS]STOP" << endl;
      ossmain << "[POCC_AEL_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AEL_STIFFNESS_TENSOR]START" << endl;
      ossmain << "[POCC_AEL_STIFFNESS_TENSOR]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      for (uint j = 0; j < elastic_tensor_ave.size(); j++) {
        for (uint k = 0; k < elastic_tensor_ave.at(j).size(); k++) {
          ossmain << setw(15) << setprecision(6) << (elastic_tensor_ave.at(j).at(k));
        }
        ossmain << endl;
      }
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AEL_STIFFNESS_TENSOR]STOP" << endl;
      ossmain << "[POCC_AEL_STIFFNESS_TENSOR]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AEL_COMPLIANCE_TENSOR]START" << endl;
      ossmain << "[POCC_AEL_COMPLIANCE_TENSOR]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      for (uint j = 0; j < compliance_tensor_ave.size(); j++) {
        for (uint k = 0; k < compliance_tensor_ave.at(j).size(); k++) {
          ossmain << setw(15) << setprecision(6) << (compliance_tensor_ave.at(j).at(k));
        }
        ossmain << endl;
      }
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AEL_COMPLIANCE_TENSOR]STOP" << endl;
      ossmain << "[POCC_AEL_COMPLIANCE_TENSOR]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      // Writes all results to a file
      if (ael_write_full_results) {
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_RESULTS]START" << endl;
        ossfull << "[POCC_AEL_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << "pocc_ael_temperature=" << T << "  (K)" << endl; //CO20200502 - no (*this).getTemperatureString(T) needed here
        ossfull << "pocc_ael_bulk_modulus_voigt=" << Bvoigt_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_bulk_modulus_reuss=" << Breuss_ave  << "  (GPa)" << endl;
        ossfull << "pocc_ael_bulk_modulus_vrh=" << Bvrh_ave << "  (GPa)" << endl;    
        ossfull << "pocc_ael_shear_modulus_voigt=" << Gvoigt_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_shear_modulus_reuss=" << Greuss_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_shear_modulus_vrh=" << Gvrh_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_youngs_modulus_voigt=" << Yvoigt_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_youngs_modulus_reuss=" << Yreuss_ave  << "  (GPa)" << endl;
        ossfull << "pocc_ael_youngs_modulus_vrh=" << Yvrh_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_poisson_ratio=" << Poisson_ratio_ave << endl;    
        ossfull << "pocc_ael_hardness_teter_voigt=" << Hv_Teter_voigt_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_teter_reuss=" << Hv_Teter_reuss_ave  << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_teter_vrh=" << Hv_Teter_vrh_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_chen_voigt=" << Hv_Chen_voigt_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_chen_reuss=" << Hv_Chen_reuss_ave  << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_chen_vrh=" << Hv_Chen_vrh_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_tian_voigt=" << Hv_Tian_voigt_ave << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_tian_reuss=" << Hv_Tian_reuss_ave  << "  (GPa)" << endl;
        ossfull << "pocc_ael_hardness_tian_vrh=" << Hv_Tian_vrh_ave << "  (GPa)" << endl;
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_RESULTS]STOP" << endl;
        ossfull << "[POCC_AEL_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_STIFFNESS_TENSOR]START" << endl;
        ossfull << "[POCC_AEL_STIFFNESS_TENSOR]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        for (uint j = 0; j < elastic_tensor_ave.size(); j++) {
          for (uint k = 0; k < elastic_tensor_ave.at(j).size(); k++) {
            ossfull << setw(15) << setprecision(6) << (elastic_tensor_ave.at(j).at(k));
          }
          ossfull << endl;
        }
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_STIFFNESS_TENSOR]STOP" << endl;
        ossfull << "[POCC_AEL_STIFFNESS_TENSOR]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_COMPLIANCE_TENSOR]START" << endl;
        ossfull << "[POCC_AEL_COMPLIANCE_TENSOR]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        for (uint j = 0; j < compliance_tensor_ave.size(); j++) {
          for (uint k = 0; k < compliance_tensor_ave.at(j).size(); k++) {
            ossfull << setw(15) << setprecision(6) << (compliance_tensor_ave.at(j).at(k));
          }
          ossfull << endl;
        }
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_COMPLIANCE_TENSOR]STOP" << endl;
        ossfull << "[POCC_AEL_COMPLIANCE_TENSOR]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
      }
    }
    string ofileafpoccaelname = m_aflags.Directory + "/aflow.pocc_ael.out";
    aurostd::stringstream2file(ossmain, ofileafpoccaelname, "WRITE");
    if (ael_write_full_results) {
      ossfull << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_RESULTS]START" << endl;
      ossfull << "[POCC_AEL_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossfull << "pocc_ael_bulk_modulus_voigt_min=" << Bvoigt_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_bulk_modulus_voigt_max=" << Bvoigt_max << "  (GPa)" << endl;    
      ossfull << "pocc_ael_bulk_modulus_reuss_min=" << Breuss_min  << "  (GPa)" << endl;
      ossfull << "pocc_ael_bulk_modulus_reuss_max=" << Breuss_max  << "  (GPa)" << endl;    
      ossfull << "pocc_ael_bulk_modulus_vrh_min=" << Bvrh_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_bulk_modulus_vrh_max=" << Bvrh_max << "  (GPa)" << endl;
      ossfull << "pocc_ael_shear_modulus_voigt_min=" << Gvoigt_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_shear_modulus_voigt_max=" << Gvoigt_max << "  (GPa)" << endl;    
      ossfull << "pocc_ael_shear_modulus_reuss_min=" << Greuss_min  << "  (GPa)" << endl;
      ossfull << "pocc_ael_shear_modulus_reuss_max=" << Greuss_max  << "  (GPa)" << endl;    
      ossfull << "pocc_ael_shear_modulus_vrh_min=" << Gvrh_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_shear_modulus_vrh_max=" << Gvrh_max << "  (GPa)" << endl;
      ossfull << "pocc_ael_youngs_modulus_voigt_min=" << Yvoigt_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_youngs_modulus_voigt_max=" << Yvoigt_max << "  (GPa)" << endl;    
      ossfull << "pocc_ael_youngs_modulus_reuss_min=" << Yreuss_min  << "  (GPa)" << endl;
      ossfull << "pocc_ael_youngs_modulus_reuss_max=" << Yreuss_max  << "  (GPa)" << endl;    
      ossfull << "pocc_ael_youngs_modulus_vrh_min=" << Yvrh_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_youngs_modulus_vrh_max=" << Yvrh_max << "  (GPa)" << endl;
      ossfull << "pocc_ael_poisson_ratio_min=" << Poisson_ratio_min << endl;
      ossfull << "pocc_ael_poisson_ratio_max=" << Poisson_ratio_max << endl;    
      ossfull << "pocc_ael_hardness_teter_voigt_min=" << Hv_Teter_voigt_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_teter_voigt_max=" << Hv_Teter_voigt_max << "  (GPa)" << endl;    
      ossfull << "pocc_ael_hardness_teter_reuss_min=" << Hv_Teter_reuss_min  << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_teter_reuss_max=" << Hv_Teter_reuss_max  << "  (GPa)" << endl;    
      ossfull << "pocc_ael_hardness_teter_vrh_min=" << Hv_Teter_vrh_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_teter_vrh_max=" << Hv_Teter_vrh_max << "  (GPa)" << endl;       
      ossfull << "pocc_ael_hardness_chen_voigt_min=" << Hv_Chen_voigt_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_chen_voigt_max=" << Hv_Chen_voigt_max << "  (GPa)" << endl;    
      ossfull << "pocc_ael_hardness_chen_reuss_min=" << Hv_Chen_reuss_min  << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_chen_reuss_max=" << Hv_Chen_reuss_max  << "  (GPa)" << endl;    
      ossfull << "pocc_ael_hardness_chen_vrh_min=" << Hv_Chen_vrh_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_chen_vrh_max=" << Hv_Chen_vrh_max << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_Tian_voigt_min=" << Hv_Tian_voigt_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_Tian_voigt_max=" << Hv_Tian_voigt_max << "  (GPa)" << endl;    
      ossfull << "pocc_ael_hardness_Tian_reuss_min=" << Hv_Tian_reuss_min  << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_Tian_reuss_max=" << Hv_Tian_reuss_max  << "  (GPa)" << endl;    
      ossfull << "pocc_ael_hardness_Tian_vrh_min=" << Hv_Tian_vrh_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_hardness_Tian_vrh_max=" << Hv_Tian_vrh_max << "  (GPa)" << endl;
      ossfull << "pocc_ael_pughs_modulus_ratio=" << Pugh_ratio_ave << endl;
      ossfull << "pocc_ael_pughs_modulus_ratio_min=" << Pugh_ratio_min << endl;
      ossfull << "pocc_ael_pughs_modulus_ratio_max=" << Pugh_ratio_max << endl;    
      ossfull << "pocc_ael_Cauchy_pressure=" << Cauchy_ave << endl;
      ossfull << "pocc_ael_cauchy_pressure_min=" << Cauchy_min << "  (GPa)" << endl;
      ossfull << "pocc_ael_cauchy_pressure_max=" << Cauchy_max << "  (GPa)" << endl;
      //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_RESULTS]STOP" << endl;
      ossfull << "[POCC_AEL_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossfull << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_CONFIG_RESULTS]START" << endl;
      ossfull << "[POCC_AEL_CONFIG_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossfull << "# Config" << "  " << aurostd::PaddedPRE("Bvoigt(GPa)",2," ") << "  " << aurostd::PaddedPRE("Breuss(GPa)",2," ") << "  " << aurostd::PaddedPRE("Bvrh(GPa)",2," ") << "  " << aurostd::PaddedPRE("Gvoigt(GPa)",2," ") << "  " << aurostd::PaddedPRE("Greuss(GPa)",2," ") << "  " << aurostd::PaddedPRE("Gvrh(GPa)",2," ") << "  " << aurostd::PaddedPRE("Evoigt(GPa)",2," ") << "  " << aurostd::PaddedPRE("Ereuss(GPa)",2," ") << "  " << aurostd::PaddedPRE("Evrh(GPa)",2," ") << "  " << aurostd::PaddedPRE("Poisson_ratio",2," ") << "  " << aurostd::PaddedPRE("Hv_Teter(GPa)",2," ") << "  " << aurostd::PaddedPRE("Hv_Chen(GPa)",2," ") << "  " << aurostd::PaddedPRE("Hv_Tian(GPa)",2," ")  << "  " << aurostd::PaddedPRE("Pugh_modulus_ratio",2," ")  << "  " << aurostd::PaddedPRE("Cauchy_pressure(GPa)",2," ")<< endl;
      for (uint i = 0; i < Bvoigt.size(); i++) {
        ossfull << setw(8) << fixed << i << "  " << setw(11) << setprecision(3) << Bvoigt.at(i) << "  " << setw(11) << setprecision(3) << Breuss.at(i) << "  " << setw(9) << setprecision(3) << Bvrh.at(i)  << "  " << setw(11) << setprecision(3) << Gvoigt.at(i) << "  " << setw(11) << setprecision(3) << Greuss.at(i) << "  " << setw(9) << setprecision(3) << Gvrh.at(i) << "  " << setw(11) << setprecision(3) << Yvoigt.at(i) << "  " << setw(11) << setprecision(3) << Yreuss.at(i) << "  " << setw(9) << setprecision(3) << Yvrh.at(i)  << "  " << setw(13) << setprecision(6) << Poisson_ratio.at(i) << "  " << setw(13) << setprecision(3) << Hv_Teter_vrh.at(i) << "  " << setw(12) << setprecision(3) << Hv_Chen_vrh.at(i) << "  " << setw(12) << setprecision(3) << Hv_Tian_vrh.at(i) << "  " << setw(18) << setprecision(6) << Pugh_ratio.at(i)  << "  " << setw(20) << setprecision(3) << Cauchy_pressure.at(i) << endl;
      }
      //[CO20200731 - OBSOLETE]ossfull << "[POCC_AEL_CONFIG_RESULTS]STOP" << endl;
      ossfull << "[POCC_AEL_CONFIG_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossfull << AFLOWIN_SEPARATION_LINE << endl;
      string ofileafpoccaelfullname = m_aflags.Directory + "/aflow.pocc_ael_full.out";
      aurostd::stringstream2file(ossfull, ofileafpoccaelfullname, "WRITE");
      ossfull.clear();
      ossfull.str(std::string());
    }
  }
} // namespace pocc


//CT20200323 - POCC+AGL functions
// ***************************************************************************
//  void POccCalculator::calculateDebyeThermalProperties(const vector<double>& v_temperatures)
//  Runs POCC+AGL
// ***************************************************************************
namespace pocc {
  void POccCalculator::calculateDebyeThermalProperties(const vector<double>& v_temperatures) {
    bool agl_run_postprocess, agl_write_full_results;
    // [OBSOLETE] uint ntemperature = 0, npressure = 0;
    // [OBSOLETE] double stemperature = 0.0, spressure = 0.0;
    vector<double> Debye_temperature, Debye_acoustic, Gruneisen, Cv300K, Cp300K, Fvib300K_atom, Fvib300K_cell, Svib300K_atom,  Svib300K_cell, kappa300K;
    vector<vector<double> > agl_temperatures, agl_gibbs_energies_atom, agl_vibrational_energies_atom;
    // [OBSOLETE] setAGLOptions(agl_run_postprocess, agl_write_full_results, ntemperature, stemperature, npressure, spressure);
    setAGLOptions(agl_run_postprocess, agl_write_full_results);
    if (agl_run_postprocess) {
      // [OBSOLETE] generateDebyeThermalProperties(ntemperature, stemperature, npressure, spressure, Debye_temperature, Debye_acoustic, Gruneisen, Cv300K, Cp300K, Fvib300K_atom, Fvib300K_cell, Svib300K_atom, Svib300K_cell, kappa300K, agl_temperatures, agl_gibbs_energies_atom, agl_vibrational_energies_atom);
      generateDebyeThermalProperties(Debye_temperature, Debye_acoustic, Gruneisen, Cv300K, Cp300K, Fvib300K_atom, Fvib300K_cell, Svib300K_atom, Svib300K_cell, kappa300K, agl_temperatures, agl_gibbs_energies_atom, agl_vibrational_energies_atom);

    } else {
      getDebyeThermalProperties(Debye_temperature, Debye_acoustic, Gruneisen, Cv300K, Cp300K, Fvib300K_atom, Fvib300K_cell, Svib300K_atom, Svib300K_cell, kappa300K, agl_temperatures, agl_gibbs_energies_atom, agl_vibrational_energies_atom);
    }
    getAverageDebyeThermalProperties(v_temperatures, agl_write_full_results, Debye_temperature, Debye_acoustic, Gruneisen, Cv300K, Cp300K, Fvib300K_atom, Fvib300K_cell, Svib300K_atom, Svib300K_cell, kappa300K, agl_temperatures, agl_gibbs_energies_atom, agl_vibrational_energies_atom);
  }
} // namespace pocc

// ***************************************************************************
// void POccCalculator::setAGLOptions(bool& agl_run_postprocess, bool& agl_write_full_results, uint& ntemperature, double& stemperature, uint& npressure, double& spressure)
// Sets options for running POCC+AGL from command line
// agl_run_postprocess: runs AGL postprocessing; otherwise reads results from aflow.agl.out
// agl_write_full_results: writes out results for all POCC configurations
// ntemperature: number of temperature steps for AGL postprocessing
// stemperature: size of temperature steps for AGL postprocessing
// npressure: number of pressure steps for AGL postprocessing
// spressure: size of pressure steps for AGL postprocessing
// ***************************************************************************
namespace pocc {
  // [OBSOLETE] void POccCalculator::setAGLOptions(bool& agl_run_postprocess, bool& agl_write_full_results, uint& ntemperature, double& stemperature, uint& npressure, double& spressure)
  void POccCalculator::setAGLOptions(bool& agl_run_postprocess, bool& agl_write_full_results) {
    // bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    agl_run_postprocess = false;
    agl_write_full_results = false;
    // [OBSOLETE] ntemperature = 0;
    // [OBSOLETE] stemperature = 0.0;
    // [OBSOLETE] npressure = 0;
    // [OBSOLETE] spressure = 0.0;

    //run postprocessing
    agl_run_postprocess=( agl_run_postprocess || XHOST.vflag_control.flag("AGL_RUN_POSTPROCESSING") );  //CO20201026
    agl_run_postprocess=( agl_run_postprocess || XHOST.POSTPROCESS );  //CO20201026

    //write extra output (debugging mostly)
    agl_write_full_results=( agl_write_full_results || XHOST.vflag_control.flag("AGL_WRITE_FULL_RESULTS") );

    // [OBSOLETE] if(XHOST.vflag_control.flag("AGL_NTEMPERATURE")){
    // [OBSOLETE]   ntemperature=aurostd::string2utype<uint>(XHOST.vflag_control.getattachedscheme("AGL_NTEMPERATURE"));
    // [OBSOLETE] }  //command line input

    // [OBSOLETE] if(XHOST.vflag_control.flag("AGL_STEMPERATURE")){
    // [OBSOLETE]   stemperature=aurostd::string2utype<double>(XHOST.vflag_control.getattachedscheme("AGL_STEMPERATURE"));
    // [OBSOLETE] }  //command line input

    // [OBSOLETE] if(XHOST.vflag_control.flag("AGL_NPRESSURE")){
    // [OBSOLETE]   npressure=aurostd::string2utype<uint>(XHOST.vflag_control.getattachedscheme("AGL_NPRESSURE"));
    // [OBSOLETE] }  //command line input

    // [OBSOLETE] if(XHOST.vflag_control.flag("AGL_SPRESSURE")){
    // [OBSOLETE]   spressure=aurostd::string2utype<double>(XHOST.vflag_control.getattachedscheme("AGL_SPRESSURE"));
    // [OBSOLETE] }  //command line input
  }
} // namespace pocc

// ***************************************************************************
//  void POccCalculator::generateDebyeThermalProperties(const uint& ntemperature, const uint& stemperature, const uint& npressure, const double& spressure, vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom)
//  Runs AGL postprocessing to generate thermal properties 
// ***************************************************************************
namespace pocc {
  // [OBSOLETE] void POccCalculator::generateDebyeThermalProperties(uint ntemperature, double stemperature, uint npressure, double spressure, vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom)
  void POccCalculator::generateDebyeThermalProperties(vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    uint aglerror = 0;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " directory=[" << m_aflags.Directory << "]" << endl;
    }
    // string aflowin;
    // aflowin=string(directory +"/"+_AFLOWIN_);
    // if(!aurostd::FileExist(aflowin)) {throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"pocc::POccCalculator::generateDebyeThermalProperties():",_RUNTIME_ERROR_);}
    // string str_AflowIn; aurostd::file2string(aflowin, str_AflowIn);
    double agl_debye_temperature = 0.0, agl_debye_acoustic = 0.0, agl_gruneisen = 0.0;
    double agl_heat_capacity_Cv_300K = 0.0, agl_heat_capacity_Cp_300K = 0.0;
    double agl_vibrational_free_energy_300K_cell = 0.0, agl_vibrational_free_energy_300K_atom = 0.0;
    double agl_vibrational_entropy_300K_cell = 0.0, agl_vibrational_entropy_300K_atom = 0.0;	
    double agl_thermal_conductivity_300K = 0.0;
    vector<double> agl_temperature;
    vector<double> agl_gibbs_energy_atom;
    vector<double> agl_vibrational_energy_atom;
    vector<string> vline, tokens;
    stringstream aflow_agl_out;
    stringstream aflow_agl_energies_out;
    unsigned long long int isupercell = 0;
    string pocc_directory_abs="";
    bool run_directory=false;
    bool agl_aflowin_found = false;
    string AflowInName = _AFLOWIN_;
    string FileLockName = _AFLOWLOCK_;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      pocc_directory_abs=m_aflags.Directory+"/"+m_ARUN_directories[isupercell];
      run_directory=false;
      agl_aflowin_found = false;
      AflowInName = _AFLOWIN_;
      FileLockName = _AFLOWLOCK_;
      agl_temperature.clear(); //CT20200729
      agl_gibbs_energy_atom.clear(); //CT20200729
      agl_vibrational_energy_atom.clear(); //CT20200729
      // [OBSOLETE] aglerror = AGL_functions::Get_ThermalProperties_AGL_postprocess(pocc_directory_abs, ntemperature, stemperature, npressure, spressure, agl_temperature, agl_gibbs_energy_atom, agl_vibrational_energy_atom);
      // [OBSOLETE] if (aglerror > 0) {
      // [OBSOLETE]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem with AGL calculation: [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);
      // [OBSOLETE] }
      //CT20200722 Run AGL postprocessing through KBIN
      aglerror = AGL_functions::AGL_Get_AflowInName(AflowInName, pocc_directory_abs, agl_aflowin_found);
      if (aglerror != 0) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem finding AGL aflow.in filename [dir="+pocc_directory_abs+"]",_FILE_ERROR_);
      }
      if (agl_aflowin_found) {
        if(aurostd::FileExist(pocc_directory_abs+"/agl.LOCK")) {
          FileLockName = "agl.LOCK";
        }
        run_directory=true;
      }
      if(run_directory){
        _aflags aflags;
        aflags.Directory=pocc_directory_abs;

        //save originals
        string _AFLOWIN_orig=_AFLOWIN_;
        string _AFLOWLOCK_orig=_AFLOWLOCK_;

        //set env for RUN_Directory()
        _AFLOWIN_=AflowInName;
        _AFLOWLOCK_=FileLockName;
        if(aurostd::FileExist(pocc_directory_abs+"/"+_AFLOWLOCK_)){aurostd::file2file(pocc_directory_abs+"/"+_AFLOWLOCK_,pocc_directory_abs+"/"+_AFLOWLOCK_+".run");} //keep original LOCK
        KBIN::RUN_Directory(aflags);

        //return to original
        _AFLOWIN_=_AFLOWIN_orig;
        _AFLOWLOCK_=_AFLOWLOCK_orig;
      }
      if(aurostd::EFileExist(pocc_directory_abs+"/aflow.agl.out")) {
        aflow_agl_out.str(aurostd::substring2string(aurostd::efile2string(pocc_directory_abs+"/"+"aflow.agl.out"),"[AGL_RESULTS]START","[AGL_RESULTS]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::StringstreamClean(aflow_agl_out);
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(pocc_directory_abs+"/"+"aflow.agl.out"),aflow_agl_out,"[AGL_RESULTS]START","[AGL_RESULTS]STOP");
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " aflow.agl.out:" << endl << aflow_agl_out.str() << endl;
        }
        aurostd::stream2vectorstring(aflow_agl_out,vline);
        for (uint i=0;i<vline.size();i++) {
          aurostd::StringSubst(vline.at(i),"="," ");
          aurostd::string2tokens(vline.at(i),tokens," ");
          // cout << "Tokens = " << tokens << endl;
          if(tokens.size()>=2) {
            if(tokens.at(0)=="agl_debye") agl_debye_temperature=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_acoustic_debye") agl_debye_acoustic=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_gruneisen") agl_gruneisen=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_heat_capacity_Cv_300K") agl_heat_capacity_Cv_300K=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_heat_capacity_Cp_300K") agl_heat_capacity_Cp_300K=aurostd::string2utype<double>(tokens.at(1));	      
            if(tokens.at(0)=="agl_vibrational_free_energy_300K_cell") agl_vibrational_free_energy_300K_cell=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_vibrational_free_energy_300K_atom") agl_vibrational_free_energy_300K_atom=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_vibrational_entropy_300K_cell") agl_vibrational_entropy_300K_cell=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_vibrational_entropy_300K_atom") agl_vibrational_entropy_300K_atom=aurostd::string2utype<double>(tokens.at(1));
            if(tokens.at(0)=="agl_thermal_conductivity_300K") agl_thermal_conductivity_300K=aurostd::string2utype<double>(tokens.at(1));
          }
        }
        aflow_agl_energies_out.str(aurostd::substring2string(aurostd::efile2string(pocc_directory_abs+"/"+"aflow.agl.out"),"[AGL_ENERGIES_TEMPERATURE]START","[AGL_ENERGIES_TEMPERATURE]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::StringstreamClean(aflow_agl_energies_out);
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(pocc_directory_abs+"/"+"aflow.agl.out"),aflow_agl_energies_out,"[AGL_ENERGIES_TEMPERATURE]START","[AGL_ENERGIES_TEMPERATURE]STOP");
        if(LDEBUG) {
          cerr << __AFLOW_FUNC__ << " aflow_agl_energies:" << endl << aflow_agl_energies_out.str() << endl;
        }
        aurostd::stream2vectorstring(aflow_agl_energies_out,vline);
        for (uint i=0;i<vline.size();i++) {
          aurostd::StringSubst(vline.at(i),"="," ");
          aurostd::string2tokens(vline.at(i),tokens," ");
          if(tokens.size()>=9) {
            if (tokens.at(0) == "#") {
              continue;
            } else {
              agl_temperature.push_back(aurostd::string2utype<double>(tokens.at(0)));
              agl_gibbs_energy_atom.push_back(aurostd::string2utype<double>(tokens.at(5)));
              agl_vibrational_energy_atom.push_back(aurostd::string2utype<double>(tokens.at(6)));
            }
          }
        }
        Debye_temperature.push_back(agl_debye_temperature);
        Debye_acoustic.push_back(agl_debye_acoustic);
        Gruneisen.push_back(agl_gruneisen);
        Cv300K.push_back(agl_heat_capacity_Cv_300K);
        Cp300K.push_back(agl_heat_capacity_Cp_300K);
        Fvib300K_atom.push_back(agl_vibrational_free_energy_300K_atom);
        Fvib300K_cell.push_back(agl_vibrational_free_energy_300K_cell);
        Svib300K_atom.push_back(agl_vibrational_entropy_300K_atom);
        Svib300K_cell.push_back(agl_vibrational_entropy_300K_cell);	
        kappa300K.push_back(agl_thermal_conductivity_300K);
        agl_temperatures.push_back(agl_temperature);
        agl_gibbs_energies_atom.push_back(agl_gibbs_energy_atom);
        agl_vibrational_energies_atom.push_back(agl_vibrational_energy_atom);
      } else {
        // [OBSOLETE] cerr << __AFLOW_FUNC__ << pocc_directory_abs << "/aflow.agl.out: File not found" << endl;
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Problem with AGL calculation: [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);
      }
    }
  }
} // namespace pocc

// ***************************************************************************
//  void POccCalculator::getDebyeThermalProperties(vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom)
//  Retrieves thermal properties from aflow.agl.out files 
// ***************************************************************************
namespace pocc {
  void POccCalculator::getDebyeThermalProperties(vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom) {
    bool LDEBUG=(FALSE || _DEBUG_POCC_AEL_AGL_ || XHOST.DEBUG);
    stringstream message;
    // [OBSOLETE] uint aglerror = 0;
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " directory=[" << m_aflags.Directory << "]" << endl;
    }
    double agl_debye_temperature = 0.0, agl_debye_acoustic = 0.0, agl_gruneisen = 0.0;
    double agl_heat_capacity_Cv_300K = 0.0, agl_heat_capacity_Cp_300K = 0.0;
    double agl_vibrational_free_energy_300K_cell = 0.0, agl_vibrational_free_energy_300K_atom = 0.0;
    double agl_vibrational_entropy_300K_cell = 0.0, agl_vibrational_entropy_300K_atom = 0.0;	
    double agl_thermal_conductivity_300K = 0.0;
    vector<double> agl_temperature;
    vector<double> agl_gibbs_energy_atom;
    vector<double> agl_vibrational_energy_atom;
    vector<string> vline, tokens;
    stringstream aflow_agl_out;
    stringstream aflow_agl_energies_out;
    unsigned long long int isupercell = 0;
    string pocc_directory_abs="";
    string agl_filename="";
    bool aglfilezipped;
    for(std::list<POccSuperCellSet>::iterator it=l_supercell_sets.begin();it!=l_supercell_sets.end();++it){
      isupercell=std::distance(l_supercell_sets.begin(),it);
      pocc_directory_abs=m_aflags.Directory+"/"+m_ARUN_directories[isupercell];
      aurostd::StringstreamClean(aflow_agl_out);
      agl_temperature.clear(); //CT20200729
      agl_gibbs_energy_atom.clear(); //CT20200729
      agl_vibrational_energy_atom.clear(); //CT20200729
      agl_filename=pocc_directory_abs+"/aflow.agl.out";
      if(aurostd::EFileExist(agl_filename,agl_filename)){
        message << "Reading " << m_ARUN_directories[isupercell]+"/aflow.agl.out";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); //CO20200404
        aflow_agl_out.str(aurostd::substring2string(aurostd::efile2string(agl_filename),"[AGL_RESULTS]START","[AGL_RESULTS]STOP",0));
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(agl_filename),aflow_agl_out,"[AGL_RESULTS]START","[AGL_RESULTS]STOP");
        aglfilezipped = true;
        //[CO20200404 - OBSOLETE]} else if (aurostd::FileExist(agl_filename))
        //[CO20200404 - OBSOLETE]  aurostd::ExtractToStringstreamEXPLICIT(aurostd::file2string(agl_filename),aflow_agl_out,"[AGL_RESULTS]START","[AGL_RESULTS]STOP");
        //[CO20200404 - OBSOLETE]  aglfilezipped = false;
    } else {
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No aflow.agl.out file found [dir="+pocc_directory_abs+"]",_FILE_NOT_FOUND_);
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " aflow.agl.out:" << endl << aflow_agl_out.str() << endl;
    }
    aurostd::stream2vectorstring(aflow_agl_out,vline);
    for (uint i=0;i<vline.size();i++) {
      aurostd::StringSubst(vline.at(i),"="," ");
      aurostd::string2tokens(vline.at(i),tokens," ");
      // cout << "Tokens = " << tokens << endl;
      if(tokens.size()>=2) {
        if(tokens.at(0)=="agl_debye") agl_debye_temperature=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_acoustic_debye") agl_debye_acoustic=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_gruneisen") agl_gruneisen=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_heat_capacity_Cv_300K") agl_heat_capacity_Cv_300K=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_heat_capacity_Cp_300K") agl_heat_capacity_Cp_300K=aurostd::string2utype<double>(tokens.at(1));	      
        if(tokens.at(0)=="agl_vibrational_free_energy_300K_cell") agl_vibrational_free_energy_300K_cell=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_vibrational_free_energy_300K_atom") agl_vibrational_free_energy_300K_atom=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_vibrational_entropy_300K_cell") agl_vibrational_entropy_300K_cell=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_vibrational_entropy_300K_atom") agl_vibrational_entropy_300K_atom=aurostd::string2utype<double>(tokens.at(1));
        if(tokens.at(0)=="agl_thermal_conductivity_300K") agl_thermal_conductivity_300K=aurostd::string2utype<double>(tokens.at(1));
      }
    }
    aurostd::StringstreamClean(aflow_agl_energies_out);
    if (aglfilezipped) {
      aflow_agl_energies_out.str(aurostd::substring2string(aurostd::efile2string(agl_filename),"[AGL_ENERGIES_TEMPERATURE]START","[AGL_ENERGIES_TEMPERATURE]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(agl_filename),aflow_agl_energies_out,"[AGL_ENERGIES_TEMPERATURE]START","[AGL_ENERGIES_TEMPERATURE]STOP");
    } else {
      aflow_agl_energies_out.str(aurostd::substring2string(aurostd::file2string(agl_filename),"[AGL_ENERGIES_TEMPERATURE]START","[AGL_ENERGIES_TEMPERATURE]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(aurostd::file2string(agl_filename),aflow_agl_energies_out,"[AGL_ENERGIES_TEMPERATURE]START","[AGL_ENERGIES_TEMPERATURE]STOP");
    }
    if(LDEBUG) {
      cerr << __AFLOW_FUNC__ << " aflow_agl_energies:" << endl << aflow_agl_energies_out.str() << endl;
    }
    aurostd::stream2vectorstring(aflow_agl_energies_out,vline);
    for (uint i=0;i<vline.size();i++) {
      aurostd::StringSubst(vline.at(i),"="," ");
      aurostd::string2tokens(vline.at(i),tokens," ");
      if(tokens.size()>=9) {
        if (tokens.at(0) == "#") {
          continue;
        } else {
          agl_temperature.push_back(aurostd::string2utype<double>(tokens.at(0)));
          agl_gibbs_energy_atom.push_back(aurostd::string2utype<double>(tokens.at(5)));
          agl_vibrational_energy_atom.push_back(aurostd::string2utype<double>(tokens.at(6)));
        }
      }
    }
    Debye_temperature.push_back(agl_debye_temperature);
    Debye_acoustic.push_back(agl_debye_acoustic);
    Gruneisen.push_back(agl_gruneisen);
    Cv300K.push_back(agl_heat_capacity_Cv_300K);
    Cp300K.push_back(agl_heat_capacity_Cp_300K);
    Fvib300K_atom.push_back(agl_vibrational_free_energy_300K_atom);
    Fvib300K_cell.push_back(agl_vibrational_free_energy_300K_cell);
    Svib300K_atom.push_back(agl_vibrational_entropy_300K_atom);
    Svib300K_cell.push_back(agl_vibrational_entropy_300K_cell);	
    kappa300K.push_back(agl_thermal_conductivity_300K);
    agl_temperatures.push_back(agl_temperature);
    agl_gibbs_energies_atom.push_back(agl_gibbs_energy_atom);
    agl_vibrational_energies_atom.push_back(agl_vibrational_energy_atom);
  }
}
} // namespace pocc

// ***************************************************************************
//   void POccCalculator::getAverageDebyeThermalProperties(const vector<double>& v_temperatures, bool& agl_write_full_results, vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom)
// Calculates thermal average of thermal properties for all POCC configurations
// Writes results to aflow.pocc_agl.out
// ***************************************************************************
namespace pocc {
  void POccCalculator::getAverageDebyeThermalProperties(const vector<double>& v_temperatures, bool agl_write_full_results, vector<double>& Debye_temperature, vector<double>& Debye_acoustic, vector<double>& Gruneisen, vector<double>& Cv300K, vector<double>& Cp300K, vector<double>& Fvib300K_atom, vector<double>& Fvib300K_cell, vector<double>& Svib300K_atom, vector<double>& Svib300K_cell, vector<double>& kappa300K, vector<vector<double> >& agl_temperatures, vector<vector<double> >& agl_gibbs_energies_atom, vector<vector<double> >& agl_vibrational_energies_atom) {
    stringstream message;
    ostringstream aus;
    ofstream FileMESSAGE;
    double T=0.0;
    stringstream ossmain, ossfull;
    vector<double> agl_gibbs_energies_atom_ave, agl_vibrational_energies_atom_ave, agl_temperatures_ave;
    double Debye_temperature_ave=0.0;
    double Debye_acoustic_ave=0.0, Gruneisen_ave=0.0;
    double Cv300K_ave=0.0, Cp300K_ave=0.0;
    double Fvib300K_atom_ave=0.0, Fvib300K_cell_ave=0.0, Svib300K_atom_ave=0.0, Svib300K_cell_ave=0.0;
    double kappa300K_ave=0.0;
    double Debyemin=0.0, Debyemax=0.0;
    double gibbs_energy_atom_ave=0.0, vibrational_energy_atom_ave=0.0, gibbs_energy_atom_ave_temp=0.0, vibrational_energy_atom_ave_temp=0.0;
    double tdiff=0.0, tdiffmin=0.0, tmaxlow=0.0;
    uint jtp = 0, jtmaxlow = 0, imaxlow = 0;
    string aelaflowin = "agl_aflow.in";

    // First, find lowest maximum AGL temperature reached for the set of configurations
    jtmaxlow = agl_temperatures.at(0).size();
    imaxlow = 0;
    if (jtmaxlow > 0) {
      tmaxlow = agl_temperatures.at(0).at(jtmaxlow-1);
    } else {
      tmaxlow = 0.0;
    }
    for (uint isupercell=1; isupercell<agl_temperatures.size(); isupercell++) {
      if (agl_temperatures.at(isupercell).size() < jtmaxlow) {
        jtmaxlow = agl_temperatures.at(isupercell).size();
        imaxlow = isupercell;
        if (jtmaxlow > 0) {
          tmaxlow = agl_temperatures.at(isupercell).at(jtmaxlow-1);
        } else {
          tmaxlow = 0.0;
        }
      }
    }
    if (tmaxlow < 1.0) {
      message << "Maximum temperature = " << tmaxlow;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_WARNING_); //CO20200404
    }
    // Next, set temperatures for averaged thermal properties equal to temperatures of smallest range
    for (uint j = 0; j < jtmaxlow; j++) {
      agl_temperatures_ave.push_back(agl_temperatures.at(imaxlow).at(j));
    }

    ossmain.clear();
    ossmain.str(std::string());
    ossfull.clear();
    ossfull.str(std::string());    
    for (uint jt = 0; jt < v_temperatures.size(); jt++) {
      T = v_temperatures.at(jt);
      message << "Averaging Debye thermal properties at T=" << T << "K";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_); //CO20200404 //CO20200502 - no (*this).getTemperatureString(T) needed here
      setPOccStructureProbabilities(T);
      Debye_temperature_ave = 0.0;
      Debye_acoustic_ave = 0.0;
      Gruneisen_ave = 0.0;
      Cv300K_ave = 0.0;
      Cp300K_ave = 0.0;
      Fvib300K_atom_ave = 0.0;
      Fvib300K_cell_ave = 0.0;
      Svib300K_atom_ave = 0.0;
      Svib300K_cell_ave = 0.0;
      kappa300K_ave = 0.0;
      gibbs_energy_atom_ave = 0.0;
      vibrational_energy_atom_ave = 0.0;
      agl_gibbs_energies_atom_ave.clear(); //CT20200729
      agl_vibrational_energies_atom_ave.clear(); //CT20200729
      unsigned long long int isupercell = 0;
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
        isupercell = std::distance(l_supercell_sets.begin(), it);
        // Boltzmann weighted averages
        Debye_temperature_ave = Debye_temperature_ave + (*it).m_probability * Debye_temperature.at(isupercell);
        Debye_acoustic_ave = Debye_acoustic_ave + (*it).m_probability * Debye_acoustic.at(isupercell);
        Gruneisen_ave = Gruneisen_ave + (*it).m_probability * Gruneisen.at(isupercell);
        Cv300K_ave = Cv300K_ave + (*it).m_probability * Cv300K.at(isupercell);
        Cp300K_ave = Cp300K_ave + (*it).m_probability * Cp300K.at(isupercell);
        Fvib300K_atom_ave = Fvib300K_atom_ave + (*it).m_probability * Fvib300K_atom.at(isupercell);
        Fvib300K_cell_ave = Fvib300K_cell_ave + (*it).m_probability * Fvib300K_cell.at(isupercell);
        Svib300K_atom_ave = Svib300K_atom_ave + (*it).m_probability * Svib300K_atom.at(isupercell);
        Svib300K_cell_ave = Svib300K_cell_ave + (*it).m_probability * Svib300K_cell.at(isupercell);
        kappa300K_ave = kappa300K_ave + (*it).m_probability * kappa300K.at(isupercell);
        // Find closest vibrational temperature to POCC temperature
        tdiffmin = fabs(agl_temperatures.at(isupercell).at(0) - T);
        jtp = 0;
        for (uint j = 0; j < agl_temperatures.at(isupercell).size(); j++) {
          tdiff = fabs(agl_temperatures.at(isupercell).at(j) - T);
          if (tdiff < tdiffmin) {
            jtp = j;
            tdiffmin = tdiff;
          }
        }
        gibbs_energy_atom_ave = gibbs_energy_atom_ave + (*it).m_probability * agl_gibbs_energies_atom.at(isupercell).at(jtp);
        vibrational_energy_atom_ave = vibrational_energy_atom_ave + (*it).m_probability * agl_vibrational_energies_atom.at(isupercell).at(jtp);
      }
      // Average Gibbs and vibrational free energies at each temperature
      // Use smallest AGL temperature range
      // Next, for each temperature in this range, average over energies of all configurations
      for (uint j = 0; j < jtmaxlow; j++) {
        gibbs_energy_atom_ave_temp = 0.0;
        vibrational_energy_atom_ave_temp = 0.0;
        isupercell = 0;
        for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
          isupercell = std::distance(l_supercell_sets.begin(), it);      
          gibbs_energy_atom_ave_temp = gibbs_energy_atom_ave_temp + ((*it).m_probability * agl_gibbs_energies_atom.at(isupercell).at(j));
          vibrational_energy_atom_ave_temp = vibrational_energy_atom_ave_temp + ((*it).m_probability * agl_vibrational_energies_atom.at(isupercell).at(j));
        }
        agl_gibbs_energies_atom_ave.push_back(gibbs_energy_atom_ave_temp);
        agl_vibrational_energies_atom_ave.push_back(vibrational_energy_atom_ave_temp);
      }
      // Find maximum and minimum values of Debye temperature
      Debyemin = Debye_temperature.at(0);
      Debyemax = Debye_temperature.at(0);
      isupercell = 0;
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
        isupercell = std::distance(l_supercell_sets.begin(), it);
        if (Debye_temperature.at(isupercell) < Debyemin) {
          Debyemin = Debye_temperature.at(isupercell);
        }
        if (Debye_temperature.at(isupercell) > Debyemax) {
          Debyemax = Debye_temperature.at(isupercell);
        }      
      }

      // Writes all results to a file
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AGL_RESULTS]START" << endl;
      ossmain << "[POCC_AGL_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << "pocc_temperature=" << T << "  (K)" << endl;
      ossmain << "pocc_agl_debye=" << Debye_temperature_ave << "  (K)" << endl; //CO20200731 - debye_temperature -> debye (matching old AGL designation)
      ossmain << "pocc_agl_acoustic_debye=" << Debye_acoustic_ave << "  (K)" << endl; //CO20200731 - debye_acoustic -> acoustic_debye (matching old AGL designation)
      ossmain << "pocc_agl_gruneisen=" << Gruneisen_ave << endl;
      ossmain << "pocc_agl_heat_capacity_Cv_300K=" << Cv300K_ave << "  (kB/cell)" << endl;
      ossmain << "pocc_agl_heat_capacity_Cp_300K=" << Cp300K_ave << "  (kB/cell)" << endl;
      ossmain << "pocc_agl_vibrational_free_energy_300K_cell=" << Fvib300K_cell_ave << "  (meV/cell)" << endl;
      ossmain << "pocc_agl_vibrational_free_energy_300K_atom=" << Fvib300K_atom_ave << "  (meV/atom)" << endl;
      ossmain << "pocc_agl_vibrational_entropy_300K_cell=" << Svib300K_cell_ave << "  (meV/cell)" << endl;
      ossmain << "pocc_agl_vibrational_entropy_300K_atom=" << Svib300K_atom_ave << "  (meV/atom)" << endl;
      ossmain << "pocc_agl_thermal_conductivity_300K=" << kappa300K_ave << "  (W/m*K)" << endl;
      ossmain << "pocc_agl_gibbs_energy_atom_ave=" << gibbs_energy_atom_ave << "  (eV/atom)" << endl;
      ossmain << "pocc_agl_vibrational_energy_atom_ave=" << vibrational_energy_atom_ave << "  (meV/atom)" << endl;    
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AGL_RESULTS]STOP" << endl;
      ossmain << "[POCC_AGL_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AGL_ENERGIES_TEMPERATURE]START" << endl;
      ossmain << "[POCC_AGL_ENERGIES_TEMPERATURE]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << "#   T(K)" << "  " << aurostd::PaddedPRE("G(eV/atom)",2," ") << "  " << aurostd::PaddedPRE("Fvib(meV/atom)",2," ") << endl;
      for (uint j = 0; j < agl_temperatures_ave.size(); j++) {
        ossmain << setw(8) << setprecision(2) << fixed << agl_temperatures_ave.at(j) << "  " << setw(10) << setprecision(4) << agl_gibbs_energies_atom_ave.at(j) << "  " << setw(14) << setprecision(4) << agl_vibrational_energies_atom_ave.at(j) << endl;
      }
      //[CO20200731 - OBSOLETE]ossmain << "[POCC_AGL_ENERGIES_TEMPERATURE]STOP" << endl;
      ossmain << "[POCC_AGL_ENERGIES_TEMPERATURE]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossmain << AFLOWIN_SEPARATION_LINE << endl;
      if (agl_write_full_results) {
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AGL_RESULTS]START" << endl;
        ossfull << "[POCC_AGL_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << "pocc_temperature=" << T << "  (K)" << endl;
        ossfull << "pocc_agl_debye_temperature=" << Debye_temperature_ave << "  (K)" << endl;
        ossfull << "pocc_agl_debye_acoustic=" << Debye_acoustic_ave << "  (K)" << endl;
        ossfull << "pocc_agl_gruneisen=" << Gruneisen_ave << endl;
        ossfull << "pocc_agl_heat_capacity_Cv_300K=" << Cv300K_ave << "  (kB/cell)" << endl;
        ossfull << "pocc_agl_heat_capacity_Cp_300K=" << Cp300K_ave << "  (kB/cell)" << endl;
        ossfull << "pocc_agl_vibrational_free_energy_300K_cell=" << Fvib300K_cell_ave << "  (meV/cell)" << endl;
        ossfull << "pocc_agl_vibrational_free_energy_300K_atom=" << Fvib300K_atom_ave << "  (meV/atom)" << endl;
        ossfull << "pocc_agl_vibrational_entropy_300K_cell=" << Svib300K_cell_ave << "  (meV/cell)" << endl;
        ossfull << "pocc_agl_vibrational_entropy_300K_atom=" << Svib300K_atom_ave << "  (meV/atom)" << endl;
        ossfull << "pocc_agl_thermal_conductivity_300K=" << kappa300K_ave << "  (W/m*K)" << endl;
        ossfull << "pocc_agl_gibbs_energy_atom_ave=" << gibbs_energy_atom_ave << "  (eV/atom)" << endl;
        ossfull << "pocc_agl_vibrational_energy_atom_ave=" << vibrational_energy_atom_ave << "  (meV/atom)" << endl;    
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AGL_RESULTS]STOP" << endl;
        ossfull << "[POCC_AGL_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AGL_ENERGIES_TEMPERATURE]START" << endl;
        ossfull << "[POCC_AGL_ENERGIES_TEMPERATURE]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        // Writes Gibbs and vibrational free energies to a file
        ossfull << "#   T(K)" << "  " << aurostd::PaddedPRE("G(eV/atom)",2," ") << "  " << aurostd::PaddedPRE("Fvib(meV/atom)",2," ") << endl;
        for (uint j = 0; j < agl_temperatures_ave.size(); j++) {
          ossfull << setw(8) << setprecision(2) << fixed << agl_temperatures_ave.at(j) << "  " << setw(10) << setprecision(4) << agl_gibbs_energies_atom_ave.at(j) << "  " << setw(14) << setprecision(4) << agl_vibrational_energies_atom_ave.at(j) << endl;
        }
        //[CO20200731 - OBSOLETE]ossfull << "[POCC_AGL_ENERGIES_TEMPERATURE]STOP" << endl;
        ossfull << "[POCC_AGL_ENERGIES_TEMPERATURE]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
        ossfull << AFLOWIN_SEPARATION_LINE << endl;
      }
    }
    string ofileafpoccaglname = m_aflags.Directory + "/aflow.pocc_agl.out";
    aurostd::stringstream2file(ossmain, ofileafpoccaglname, "WRITE");
    ossmain.clear();
    ossmain.str(std::string());
    if (agl_write_full_results) {
      ossfull << AFLOWIN_SEPARATION_LINE << endl;
      //[CO20200731 - OBSOLETE]ossfull << "[POCC_AGL_CONFIG_RESULTS]START" << endl;
      ossfull << "[POCC_AGL_CONFIG_RESULTS]START_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossfull << "# Config" << "  " << aurostd::PaddedPRE("Debye_temperature(K)",2," ") << "  " << aurostd::PaddedPRE("Debye_acoustic(K)",2," ") << "  " << aurostd::PaddedPRE("Gruneisen_parameter",2," ") << "  " << aurostd::PaddedPRE("Cv(kB/cell)",2," ") << "  " << aurostd::PaddedPRE("Cp(kB/cell)",2," ") << "  " << aurostd::PaddedPRE("Fvib/cell(meV/cell)",2," ") << "  " << aurostd::PaddedPRE("Fvib/atom(meV/atom)",2," ") << "  " << aurostd::PaddedPRE("Svib/cell(meV/cell*K)",2," ") << "  " << aurostd::PaddedPRE("Svib/atom(meV/atom*K)",2," ") << "  " << aurostd::PaddedPRE("Thermal_conductivity(W/(m*K))",2," ") << endl;
      for (uint i = 0; i < Debye_temperature.size(); i++) {
        ossfull << setw(8) << fixed << i << "  " << setw(20) << setprecision(3) << Debye_temperature.at(i) << "  " << setw(17) << setprecision(3) << Debye_acoustic.at(i) << "  " << setw(19) << setprecision(5) << Gruneisen.at(i)  << "  " << setw(11) << setprecision(4) << Cv300K.at(i) << "  " << setw(11) << setprecision(4) << Cp300K.at(i) << "  " << setw(19) << setprecision(4) << Fvib300K_cell.at(i) << "  " << setw(19) << setprecision(4) << Fvib300K_atom.at(i) << "  " << setw(21) << setprecision(5) << Svib300K_cell.at(i) << "  " << setw(21) << setprecision(6) << Svib300K_atom.at(i)  << "  " << setw(29) << setprecision(6) << kappa300K.at(i) << endl;
      }
      //[CO20200731 - OBSOLETE]ossfull << "[POCC_AGL_CONFIG_RESULTS]STOP" << endl;
      ossfull << "[POCC_AGL_CONFIG_RESULTS]STOP_TEMPERATURE=" << (*this).getTemperatureString(T) << "_K" << endl;
      ossfull << AFLOWIN_SEPARATION_LINE << endl;
      string ofileafpoccaglenergname = m_aflags.Directory + "/aflow.pocc_agl_full.out";
      aurostd::stringstream2file(ossfull, ofileafpoccaglenergname, "WRITE");
      ossfull.clear();
      ossfull.str(std::string());
    }
  }
} // namespace pocc    

#endif  // _AFLOW_POCC_AEL_AGL_CPP_  //CO20200404

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                Aflow Cormac Toher - Duke University 2013-2021           *
// *                                                                         *
//****************************************************************************
