// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2020
#ifndef _AFLOWLIB_LIBRARIES_SCRUBBER_CPP
#define _AFLOWLIB_LIBRARIES_SCRUBBER_CPP
#include "aflow.h"

#define acerr std::cerr
#define acout std::cout

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

namespace aflowlib {
  bool aflowlib2stream(const aflowlib::_aflowlib_entry& data,const string& file,stringstream& stream,bool VERBOSE) {
    return aurostd::url2stringstream(data.aurl+"/"+file,stream,VERBOSE);
  }
}

namespace aflowlib {
  bool aflowlib2file(const aflowlib::_aflowlib_entry& data,const string& file,bool VERBOSE) {
    stringstream stream;
    bool out=aurostd::url2stringstream(data.aurl+"/"+file,stream,VERBOSE);
    aurostd::stringstream2file(stream,file);
    return out;
  }
}

// will be moved near LI2RAW
namespace aflowlib {
  uint LIB2SCRUB(string library,bool VERBOSE) {
    if(VERBOSE) acerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    vector<string> vlib;
    uint fixes=0;
    if(library=="ALL" || library=="all") {
      aurostd::string2tokens("LIB0,LIB1,LIB2,LIB3,LIB4,LIB5,LIB6,LIB7,LIB8,LIB9,ICSD",vlib,","); // not default vlib.push_back("AUID");
    }
    if(library=="LIB0" || library=="lib0") { vlib.push_back("LIB0");}
    if(library=="LIB1" || library=="lib1") { vlib.push_back("LIB1");}
    if(library=="LIB2" || library=="lib2") { vlib.push_back("LIB2");}
    if(library=="LIB3" || library=="lib3") { vlib.push_back("LIB3");}
    if(library=="LIB4" || library=="lib4") { vlib.push_back("LIB4");}
    if(library=="LIB5" || library=="lib5") { vlib.push_back("LIB5");}
    if(library=="LIB6" || library=="lib6") { vlib.push_back("LIB6");}
    if(library=="LIB7" || library=="lib7") { vlib.push_back("LIB7");}
    if(library=="LIB8" || library=="lib8") { vlib.push_back("LIB8");}
    if(library=="LIB9" || library=="lib9") { vlib.push_back("LIB9");}
    if(library=="ICSD" || library=="icsd") { vlib.push_back("ICSD");}
    if(library=="AUID" || library=="auid") { vlib.push_back("AUID");}

    for(uint i=0;i<vlib.size();i++) {
      if(VERBOSE) acerr << __AFLOW_FUNC__ << " **********************************" << endl;
      if(VERBOSE) acerr << __AFLOW_FUNC__ << " TESTING vlib.at(i)=" << vlib.at(i) << endl;
      vector<string> list2found;
      vector<string> listLIB2RAW,listRM,listANRL,listANRL_subst,listENTHALPY,listppAUID,listINCOMPLETE,listAGL2FIX,listTOUCH,listLIB2AUID,listREMOVE_MARYLOU,listICSD2LINK,listZIPNOMIX,listBROKEN,listAUIDerrors;
      stringstream ossLIB2RAW,ossRM,ossANRL,ossANRL_subst,ossENTHALPY,ossppAUID,ossINCOMPLETE,ossAGL2FIX,ossTOUCH,ossLIB2AUID,ossREMOVE_MARYLOU,ossICSD2LINK,ossZIPNOMIX,ossBROKEN,ossAUIDerrors;
      vector<string> tokens;

      if(vlib.at(i)=="AUID") {
        uint AUID_found=0;
        uint AUID_errors=0;
        _aflags aflags;aflags.Directory=".";
        aurostd::string2tokens("0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f",tokens,",");
        for(uint j=0;j<tokens.size();j++) {
          for(uint k=0;k<tokens.size();k++) {
            for(uint l=0;l<tokens.size();l++) {
              for(uint m=0;m<tokens.size();m++) {
                string digits_jklm="aflow:"+tokens.at(j)+tokens.at(k)+"/"+tokens.at(l)+tokens.at(m)+"/";
                // 	acout << __AFLOW_FUNC__ << " TESTING /common/" << vlib.at(i) << "/RAW/"  << digits_jklm << "*.json links" << endl;
                //	acout << __AFLOW_FUNC__ << " TIME" << Message(__AFLOW_FILE__,aflags,"user,host,time") << endl;
                aurostd::string2vectorstring(aurostd::execute2string("find /common/AUID/"+digits_jklm+"* -name RAW"),list2found);
                //	acout << __AFLOW_FUNC__ << " TIME" << Message(__AFLOW_FILE__,aflags,"user,host,time") << endl;
                //	acout << __AFLOW_FUNC__ << " list2found.size()=" << list2found.size() << endl;
                //	acout << __AFLOW_FUNC__ << " ordering " << endl;
                sort(list2found.begin(),list2found.end());
                //	acout << __AFLOW_FUNC__ << " list2found.size()=" << list2found.size() << endl;
                //  acout << __AFLOW_FUNC__ << " TIME" << Message(__AFLOW_FILE__,aflags,"user,host,time") << endl;
                for(uint n=0;n<list2found.size();n++) {
                  AUID_found++;
                  string filename_auid_json=list2found.at(n),filename_json;
                  aurostd::StringSubst(filename_auid_json,"/common/AUID/","");
                  aurostd::StringSubst(filename_auid_json,"RAW","");
                  aurostd::StringSubst(filename_auid_json,"/","");
                  filename_auid_json=list2found.at(n)+"/"+filename_auid_json+".json";
                  filename_json=list2found.at(n)+"/aflowlib.json";	  
                  if(!aurostd::FileExist(filename_auid_json)) {
                    AUID_errors++;
                    acout << __AFLOW_FUNC__ << " MISSING " << filename_auid_json << "   AUID_found=" << AUID_found << "   AUID_errors=" << AUID_errors << "  filename_json=" << filename_json << endl;
                    listAUIDerrors.push_back(list2found.at(n));
                    ossAUIDerrors << "fix " << list2found.at(n) << endl;
                  }
                }
              }
            }
          }
        }
        acout << __AFLOW_FUNC__ << " AUID_found=" << AUID_found << endl;
        acout << __AFLOW_FUNC__ << " AUID_errors=" << AUID_errors << endl;
      }      
      if(vlib.at(i)!="AUID") {
        //   if(level==0) // just LIB/aflow.in RAW/aflowlib.out
        // testing all libraries
        acerr << __AFLOW_FUNC__ << " TESTING /common/" << vlib.at(i) << "/LIB/*/" << _AFLOWIN_ << " <=> /common/" << vlib.at(i) << "/RAW/*/aflowlib.out" << endl;
        aurostd::string2vectorstring(aurostd::execute2string("find /common/"+vlib.at(i)+"/LIB -name "+_AFLOWIN_),list2found);
        acerr << __AFLOW_FUNC__ << " list2found.size()=" << list2found.size() << endl;
        //      acerr << __AFLOW_FUNC__ << " ordering " << endl;
        sort(list2found.begin(),list2found.end());
        //     acerr << __AFLOW_FUNC__ << " list2found.size()=" << list2found.size() << endl;

        vector<string> vremoveALL;
        aurostd::string2tokens("aflow.in~,agl_aflow.in~,WAVECAR.xz,REPORT.xz,POTCAR.relax1.xz,POTCAR.relax2.xz,POTCAR.relax3.xz,POTCAR.static.xz,POTCAR.bands.xz,AECCAR0.xz,AECCAR1.xz,AECCAR2.xz,AECCAR1.static.xz,AECCAR0.bands.xz,AECCAR1.bands.xz,AECCAR2.bands.xz,AECCAR0.relax1.xz,AECCAR1.relax1.xz,AECCAR2.relax1.xz,AECCAR0.relax2.xz,AECCAR1.relax2.xz,AECCAR2.relax2.xz",vremoveALL,",");

        vector<string> vremoveLIB6_LIB7;
        aurostd::string2tokens("CHGCAR,CHG,EIGENVAL,PROCAR",vremoveLIB6_LIB7,",");

        bool LIB2AUID=FALSE;//TRUE;
        bool REMOVE_MARYLOU=FALSE;//TRUE;
        bool ICSD2LINK=TRUE;//TRUE;
        bool ZIPNOMIX=TRUE;
        bool BROKEN=FALSE;//TRUE; // VERY VERY SLOW
        bool TOUCH=FALSE;

        deque<string> vext; aurostd::string2tokens(".xz",vext,",");
        deque<string> vcmd; aurostd::string2tokens("xzcat",vcmd,",");
        //  deque<string> vrelax; aurostd::string2tokens(".relax1,.relax2,.relax3,.static,.bands",vrelax,",");
        deque<string> vrelax; aurostd::string2tokens(".relax1",vrelax,",");
        //     deque<string> vbroken;aurostd::string2tokens("OUTCAR,CHG,CHGCAR,PROCAR,EIGENVAL,vasprun.xml",vbroken,",");
        deque<string> vbroken;aurostd::string2tokens("OUTCAR",vbroken,",");

        for(uint j=0;j<list2found.size();j++) {
          string directory_LIB=list2found.at(j);
          aurostd::StringSubst(directory_LIB,"/"+_AFLOWIN_,"");
          string directory_RAW=directory_LIB;
          aurostd::StringSubst(directory_RAW,"/LIB/","/RAW/");
          string directory_WEB=directory_LIB;
          aurostd::StringSubst(directory_WEB,"/LIB/","/WEB/");

          bool FileExist_directory_LIB_AFLOW_IN=aurostd::FileExist(directory_LIB+"/"+_AFLOWIN_);  // 6 times
          bool FileExist_directory_RAW_AFLOW_IN=aurostd::FileExist(directory_RAW+"/"+_AFLOWIN_);  // 5 times
          bool FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT=aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT); // 8 times
          bool FileExist_directory_RAW_AFLOWLIB_ENTRY_JSON=aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON); // 1 times

          vector<string> files2found;
          //	aurostd::execute2string("find "+directory_RAW);
          // aurostd::string2vectorstring(aurostd::execute2string("find "+directory_RAW),files2found);

          // emergency check LIB2AUID from old to new

          if(ICSD2LINK && vlib.at(i)=="ICSD") {
            aurostd::string2tokens(directory_LIB,tokens,"_");
            if(tokens.size()>2) {
              if(tokens.at(tokens.size()-2)=="ICSD") {
                if(FileExist_directory_LIB_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT && aurostd::FileExist(directory_WEB+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT))  // check aflowlib.out
                { //CO20200106 - patching for auto-indenting
                  // acerr << directory_LIB << " " << directory_RAW << " " << directory_WEB << endl;
                  string directory_ICSD2LINK=init::AFLOW_Projects_Directories("AUID")+"/icsd:/"+tokens.at(tokens.size()-1);
                  if(!aurostd::FileExist(directory_ICSD2LINK+"/LIB")) {
                    //	    acerr << directory_ICSD2LINK << endl;
                    listICSD2LINK.push_back(directory_ICSD2LINK);
                    ossICSD2LINK << "mkdir -pv " << directory_ICSD2LINK << endl;
                    ossICSD2LINK << "rm -fv " << directory_ICSD2LINK << "/LIB" << endl;
                    ossICSD2LINK << "ln -sfv " << directory_LIB << " " << directory_ICSD2LINK << "/LIB" << endl;
                    ossICSD2LINK << "rm -fv " << directory_ICSD2LINK << "/RAW" << endl;
                    ossICSD2LINK << "ln -sfv " << directory_RAW << " " << directory_ICSD2LINK << "/RAW" << endl;
                    ossICSD2LINK << "rm -fv " << directory_ICSD2LINK << "/WEB" << endl;
                    ossICSD2LINK << "ln -sfv " << directory_WEB << " " << directory_ICSD2LINK << "/WEB" << endl;
                    fixes++; 
                  }
                }
              }
            }
          }

          // check REMOVE MARYLOU	
          if(REMOVE_MARYLOU) {
            if(FileExist_directory_LIB_AFLOW_IN) {
              string directory_MARYLOU="~/LIBS/"+directory_LIB;
              aurostd::StringSubst(directory_MARYLOU,"common","");
              aurostd::StringSubst(directory_MARYLOU,"//","");
              aurostd::StringSubst(directory_MARYLOU,"//","");

              //    acerr << __AFLOW_FUNC__ << " fixing " << directory_LIB << endl;
              listREMOVE_MARYLOU.push_back(directory_MARYLOU);
              ossREMOVE_MARYLOU << "rm -rfv \"" << directory_MARYLOU << "\"" << endl;
              //	    fixes++;
            }
          }

          // check aflow.immiscibility.out   //SC20200318
          if(ZIPNOMIX) { //SC20200318
            if(FileExist_directory_LIB_AFLOW_IN) {
              if(aurostd::FileExist(directory_LIB+"/"+"aflow.immiscibility.out")) {
                //    acerr << __AFLOW_FUNC__ << " fixing " << directory_LIB << endl;
                listZIPNOMIX.push_back(directory_LIB);
                ossZIPNOMIX << "zip -9rmv /common/" << vlib.at(i) << "/nomix.zip \"" << directory_LIB << "\"" << endl;
                fixes++;
              }
            }
          }

          // check BROKEN   //SC20200319
          if(BROKEN) {    //SC20200319
            if(FileExist_directory_LIB_AFLOW_IN) {
              bool failed=FALSE;
              for(uint ibroken=0;ibroken<vbroken.size()&&!failed;ibroken++) {
                for(uint irelax=0;irelax<vrelax.size()&&!failed;irelax++) {
                  for(uint iext=0;iext<vext.size()&&!failed;iext++) {
                    if(!failed) {
                      if(aurostd::FileExist(directory_LIB+"/"+vbroken.at(ibroken)+vrelax.at(irelax)+vext.at(iext))) {
                        int answer=aurostd::execute2utype<int>(vcmd.at(iext)+" \""+directory_LIB+"/"+vbroken.at(ibroken)+vrelax.at(irelax)+vext.at(iext)+"\" 2>&1 | grep -c \"Unexpected end of input\" ");
                        if(answer!=0) {
                          failed=TRUE; // so I step out quicker
                          listBROKEN.push_back(directory_LIB);
                          ossBROKEN << "zip -9rmv /common/" << vlib.at(i) << "/nomix.zip \"" << directory_LIB << "\"" << endl;
                          fixes++;
                          acerr << "BROKEN=" << directory_LIB << endl;
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          // check LIB2AUID MISSING	
          if(LIB2AUID) {
            if(FileExist_directory_LIB_AFLOW_IN || FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT) {
              if(aflowlib::LIB2AUID(directory_LIB,TRUE,FALSE)) {	    
                //    acerr << __AFLOW_FUNC__ << " fixing " << directory_LIB << endl;
                listLIB2AUID.push_back(directory_LIB);
                if(_AFLOWIN_=="aflow.in") {
                  ossLIB2AUID << "aflow "<< "--use_aflow.in=" << _AFLOWIN_ << " --lib2auid=\"" << directory_LIB << "\"" << endl;
                  fixes++;
                }
                //	  if(_AFLOWIN_==_AFLOWIN_AGL_DEFAULT_ && !aurostd::FileExist(directory_LIB+"/LOCK") && aurostd::FileExist(directory_LIB+"/agl.LOCK"))
                if(_AFLOWIN_==_AFLOWIN_AGL_DEFAULT_ && aurostd::FileExist(directory_LIB+"/agl.LOCK")) 
                { //CO20200106 - patching for auto-indenting
                  ossLIB2AUID << "aflow "<< "--use_aflow.in=agl_aflow.in --use_LOCK=agl.LOCK " << " --lib2auid=\"" << directory_LIB << "\"" << endl;
                  fixes++;
                }
              }
            }
          }

          // check LIB2RAW MISSING	
          if(!FileExist_directory_RAW_AFLOW_IN || !FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT || !aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON)) {
            if(!aurostd::FileExist(directory_LIB+"/"+"aflow.immiscibility.out")) {
              //    acerr << __AFLOW_FUNC__ << " fixing " << directory_LIB << endl;
              listLIB2RAW.push_back(directory_LIB);
              if(_AFLOWIN_=="aflow.in") {
                ossLIB2RAW << "aflow "<< "--use_aflow.in=" << _AFLOWIN_ << " --beep --force --showPID --lib2raw=\"" << directory_LIB << "\"" << endl;
                fixes++;
              }
              //	  if(_AFLOWIN_==_AFLOWIN_AGL_DEFAULT_ && !aurostd::FileExist(directory_LIB+"/LOCK") && aurostd::FileExist(directory_LIB+"/agl.LOCK"))
              if(_AFLOWIN_==_AFLOWIN_AGL_DEFAULT_ && aurostd::FileExist(directory_LIB+"/agl.LOCK")) 
              { //CO20200106 - patching for auto-indenting
                ossLIB2RAW << "aflow "<< "--use_aflow.in=agl_aflow.in --use_LOCK=agl.LOCK " << " --beep --force --showPID --lib2raw=\"" << directory_LIB << "\"" << endl;
                fixes++;
              }
            }
          }
          // check LIB2RAW EXISTENT BUT MESSED UP
          if(FileExist_directory_RAW_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT) {
            //    acerr << __AFLOW_FUNC__ << " fixing " << directory_LIB << endl;
            if(aurostd::FileExist(directory_RAW+"/aflow.fgroup.orig.json") ||
                aurostd::FileExist(directory_RAW+"/aflow.pgroupk_xtal.relax.json") || // force
                aurostd::FileExist(directory_RAW+"/aflow.pgroup_xtal.relax.out")) {
              //	    acerr << __AFLOW_FUNC__ << " FOUND OVERWRITTEN = " << directory_RAW << " " << endl;
              listINCOMPLETE.push_back(directory_LIB);
              ossINCOMPLETE << "aflow "<< "--use_aflow.in=" << _AFLOWIN_ << " --beep --force --showPID --lib2raw=\"" << directory_LIB << "\"" << endl;
              fixes++;
            }	
          }

          // check REMOVE normal
          for(uint k=0;k<vremoveALL.size();k++) {
            if(aurostd::FileExist(directory_LIB+"/"+vremoveALL.at(k))) {
              //   acerr << __AFLOW_FUNC__ << " removing " << directory_LIB << "/" << vremoveALL.at(k) << endl;
              listRM.push_back(directory_LIB+"/"+vremoveALL.at(k));
              ossRM << "rm -fv \"" << directory_LIB << "/" << vremoveALL.at(k) << "\"" << endl;
              fixes++;	    
              // [SAFETY]	    aurostd::RemoveFile(directory_LIB+"/"+vremoveALL.at(k));
            }
          }

          // check REMOVE LIB6_LIB7
          //	if(vlib.at(i)=="LIB5" || vlib.at(i)=="LIB6" || vlib.at(i)=="LIB7")
          if(vlib.at(i)=="LIB6" || vlib.at(i)=="LIB7") {
            //	  acerr << "LIB6_LIB7" << endl;
            for(uint k=0;k<vremoveLIB6_LIB7.size();k++) {
              string FILE_relax1,FILE_relax2,FILE_static;
              FILE_relax1=directory_LIB+"/"+vremoveLIB6_LIB7.at(k)+".relax1.xz";
              FILE_relax2=directory_LIB+"/"+vremoveLIB6_LIB7.at(k)+".relax2.xz";
              FILE_static=directory_LIB+"/"+vremoveLIB6_LIB7.at(k)+".static.xz";
              if(aurostd::FileExist(FILE_static)) {
                if(aurostd::FileExist(FILE_relax1)) {
                  //   acerr << __AFLOW_FUNC__ << " removing " << FILE_relax1 << endl;
                  listRM.push_back(FILE_relax1);
                  ossRM << "rm -fv \"" << FILE_relax1 << "\"" << endl;
                  fixes++; }
                if(aurostd::FileExist(FILE_relax2)) {
                  //   acerr << __AFLOW_FUNC__ << " removing " << FILE_relax2 << endl;
                  listRM.push_back(FILE_relax2);
                  ossRM << "rm -fv \"" << FILE_relax2 << "\"" << endl;
                  fixes++; }
              }
            }
          }

          // check AGL_FIX
          if(aurostd::FileExist(directory_LIB+"/"+_AFLOWIN_AGL_DEFAULT_) && aurostd::FileExist(directory_LIB+"/LOCK") && !aurostd::FileExist(directory_LIB+"/agl.LOCK")) {
            listAGL2FIX.push_back(directory_LIB+"/"+_AFLOWIN_AGL_DEFAULT_);
            ossAGL2FIX << "cp \"" << directory_LIB << "/" << "LOCK\"" << " \"" << directory_LIB << "/" << "agl.LOCK\"" << endl;
            fixes++;	    
          }

          // check LIB2RAW - ANRL	
          if(FileExist_directory_RAW_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT) {
            if(!aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"anrl_label") && !aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"aflow_prototype_label")) {
              listANRL.push_back(directory_LIB);
              ossANRL << "aflow "<< "--use_aflow.in=" << _AFLOWIN_ << " --beep --force --showPID --lib2raw=\"" << directory_LIB << "\"" << endl;
              fixes++;
            }	
          }
          // check LIB2RAW - ANRL_subst
          if(FileExist_directory_RAW_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT) {
            if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"anrl")) { // to speed up
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"anrl_label")) { // anrl_label
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst anrl_label aflow_prototype_label " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"anrl_parameter_list")) { // anrl_parameter_list
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst anrl_parameter_list aflow_prototype_params_list " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"anrl_parameter_values")) { // anrl_parameter_values
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst anrl_parameter_values aflow_prototype_params_values " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"aflow_prototype_parameter_list")) { // aflow_prototype_parameter_list
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst aflow_prototype_parameter_list aflow_prototype_params_list " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"aflow_prototype_parameter_values")) { // aflow_prototype_parameter_values
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst aflow_prototype_parameter_values aflow_prototype_params_values " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "~" << endl;
                fixes++;
              }	    
            }
          }
          if(FileExist_directory_RAW_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_JSON) {
            if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON),"anrl")) { // to speed up
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON),"anrl_label")) { // anrl_label
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst anrl_label aflow_prototype_label " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON),"anrl_parameter_list")) { // anrl_parameter_list
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst anrl_parameter_list aflow_prototype_params_list " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON),"anrl_parameter_values")) { // anrl_parameter_values
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst anrl_parameter_values aflow_prototype_params_values " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON),"aflow_prototype_parameter_list")) { // aflow_prototype_parameter_list
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst aflow_prototype_parameter_list aflow_prototype_params_list " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "~" << endl;
                fixes++;
              }
              if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON),"aflow_prototype_parameter_values")) { // aflow_prototype_parameter_values
                listANRL_subst.push_back(directory_LIB);
                ossANRL_subst << "subst aflow_prototype_parameter_values aflow_prototype_params_values " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << " && " << " rm -fv " << directory_RAW << "/" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "~" << endl;
                fixes++;
              }	    
            }
          }

          // check LIB2RAW - ENTHALPY	
          if(FileExist_directory_RAW_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT) {
            if(!aurostd::substring2bool(directory_LIB,"LDAU2")) {
              if(!aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"enthalpy_formation_atom")) {
                listENTHALPY.push_back(directory_LIB);
                ossENTHALPY << "aflow "<< "--use_aflow.in=" << _AFLOWIN_ << " --beep --force --showPID --lib2raw=\"" << directory_LIB << "\"" << endl;
                fixes++;
              }	
            }
          }

          // check LIB2RAW - ppAUID	
          if(FileExist_directory_RAW_AFLOW_IN && FileExist_directory_RAW_AFLOWLIB_ENTRY_OUT) {
            if(!aurostd::substring2bool(directory_LIB,"LDAU2")) {
              bool issueFOUND=FALSE;
              if(!issueFOUND) issueFOUND=!aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"species_pp_AUID");
              //  [OBSOLETE]              if(!issueFOUND) { if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"6bf17162620b7ce3")) { issueFOUND=TRUE; }} // Cl:PAW_PBE:17Jan2003 changed to diatom
              //  [OBSOLETE] if(!issueFOUND) { if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"85ded82734544fa9")) { issueFOUND=TRUE; }} // O:PAW_PBE:08Apr2002 changed to diatom
              //  [OBSOLETE] if(!issueFOUND) { if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"5c530b0eaed90a00")) { issueFOUND=TRUE; }} // O:PAW_PBE_KIN:08Apr2002 changed to diatom
              //  [OBSOLETE] if(!issueFOUND) { if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"4ce4637079c90d89")) { issueFOUND=TRUE; }} // N:PAW_PBE_KIN:08Apr2002 changed to diatom
              //  [OBSOLETE] if(!issueFOUND) { if(aurostd::substring2bool(aurostd::file2string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT),"3cd4003f27559a49")) { issueFOUND=TRUE; }} // N:PAW_PBE:08Apr2002 changed to diatom
              if(issueFOUND) {
                listppAUID.push_back(directory_LIB);
                ossppAUID << "aflow "<< "--use_aflow.in=" << _AFLOWIN_ << " --beep --force --showPID --lib2raw=\"" << directory_LIB << "\"" << endl;
                fixes++;
              }
            }
          }

          // RETOUCHING DATE OF AFLOW.IN TO REPRESENT AFLOW.END.OUT
          if(TOUCH) { // TOO MANY
            if(FileExist_directory_LIB_AFLOW_IN && aurostd::FileExist(directory_LIB+"/aflow.end.out")) {
              struct stat fileInfo_IN,fileInfo_OUT;
              stat(string(directory_LIB+"/"+_AFLOWIN_).c_str(), &fileInfo_IN);
              stat(string(directory_LIB+"/aflow.end.out").c_str(), &fileInfo_OUT);
              if(fileInfo_IN.st_mtime!=fileInfo_OUT.st_mtime) {
                // acout << "FILE=" << string(directory_LIB+"/aflow.in") << ":" << fileInfo_IN.st_mtime << endl;
                // acout << "FILE=" << string(directory_LIB+"/aflow.end.out") << ":" << fileInfo_OUT.st_mtime << endl;
                listTOUCH.push_back(directory_LIB);
                string date=std::ctime(&fileInfo_OUT.st_mtime);
                if (!date.empty() && date[date.length()-1] == '\n') date.erase(date.length()-1); // remove last newline
                ossTOUCH << "echo " << directory_LIB << " && " << "touch -m --date=\"" << date << "\" " << string(directory_LIB+"/aflow.in") << " " << string(directory_LIB+"/aflow.end.out") << " " << string(directory_LIB+"/LOCK*") << " " << string(directory_LIB+"/*.xz") << endl;
                fixes++;
                stringstream sss;
                sss << "touch -m --date=\"" << date << "\" " << string(directory_LIB+"/aflow.in") << " " << string(directory_LIB+"/aflow.end.out") << " " << string(directory_LIB+"/LOCK*") << " " << string(directory_LIB+"/*.xz");
                //	    aurostd::execute(sss);
                //  acout << "FIXED " << directory_LIB << endl;
              }
            }
          }

          // some step debug
          aurostd::ProgressBar(acerr,"aflowlib::LIB2SCRUB ",j,list2found.size(),1,1,1);
        }
      }
      acerr << __AFLOW_FUNC__ << " listLIB2RAW.size()=" << listLIB2RAW.size() << endl;
      if(listLIB2RAW.size()) {
        aurostd::stringstream2file(ossLIB2RAW,XHOST.tmpfs+"/xscrubber."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listICSD2LINK.size()=" << listICSD2LINK.size() << endl;
      if(listICSD2LINK.size()) {
        aurostd::stringstream2file(ossICSD2LINK,XHOST.tmpfs+"/xscrubber_ICSD2LINK."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_ICSD2LINK."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listZIPNOMIX.size()=" << listZIPNOMIX.size() << endl;
      if(listZIPNOMIX.size()) {
        aurostd::stringstream2file(ossZIPNOMIX,XHOST.tmpfs+"/xscrubber_ZIPNOMIX."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_ZIPNOMIX."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listBROKEN.size()=" << listBROKEN.size() << endl;
      if(listBROKEN.size()) {
        aurostd::stringstream2file(ossBROKEN,XHOST.tmpfs+"/xscrubber_BROKEN."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_BROKEN."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listLIB2AUID.size()=" << listLIB2AUID.size() << endl;
      if(listLIB2AUID.size()) {
        aurostd::stringstream2file(ossLIB2AUID,XHOST.tmpfs+"/xscrubber_LIB2AUID."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_LIB2AUID."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listREMOVE_MARYLOU.size()=" << listREMOVE_MARYLOU.size() << endl;
      if(listREMOVE_MARYLOU.size()) {
        aurostd::stringstream2file(ossREMOVE_MARYLOU,XHOST.tmpfs+"/xscrubber_REMOVE_MARYLOU."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_REMOVE_MARYLOU."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listRM.size()=" << listRM.size() << endl;
      if(listRM.size()) {
        aurostd::stringstream2file(ossRM,XHOST.tmpfs+"/xscrubber_RM."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_RM."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listANRL.size()=" << listANRL.size() << endl;
      if(listANRL.size()) {
        aurostd::stringstream2file(ossANRL,XHOST.tmpfs+"/xscrubber_ANRL."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_ANRL."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listANRL_subst.size()=" << listANRL_subst.size() << endl;
      if(listANRL_subst.size()) {
        aurostd::stringstream2file(ossANRL_subst,XHOST.tmpfs+"/xscrubber_ANRL_subst."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_ANRL_subst."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listENTHALPY.size()=" << listENTHALPY.size() << endl;
      if(listENTHALPY.size()) {
        aurostd::stringstream2file(ossENTHALPY,XHOST.tmpfs+"/xscrubber_ENTHALPY."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_ENTHALPY."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listppAUID.size()=" << listppAUID.size() << endl;
      if(listppAUID.size()) {
        aurostd::stringstream2file(ossppAUID,XHOST.tmpfs+"/xscrubber_ppAUID."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_ppAUID."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listINCOMPLETE.size()=" << listINCOMPLETE.size() << endl;
      if(listINCOMPLETE.size()) {
        aurostd::stringstream2file(ossINCOMPLETE,XHOST.tmpfs+"/xscrubber_INCOMPLETE."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_INCOMPLETE."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listAGL2FIX.size()=" << listAGL2FIX.size() << endl;
      if(listAGL2FIX.size()) {
        aurostd::stringstream2file(ossAGL2FIX,XHOST.tmpfs+"/xscrubber_AGL2FIX."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_AGL2FIX."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listTOUCH.size()=" << listTOUCH.size() << endl;
      if(listTOUCH.size()) {
        aurostd::stringstream2file(ossTOUCH,XHOST.tmpfs+"/xscrubber_TOUCH."+vlib.at(i));
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_TOUCH."+vlib.at(i));
      }
      acerr << __AFLOW_FUNC__ << " listAUIDerrors.size()=" << listAUIDerrors.size() << endl;
      if(listAUIDerrors.size()) {
        aurostd::stringstream2file(ossAUIDerrors,XHOST.tmpfs+"/xscrubber_AUIDerrors");
        aurostd::ChmodFile("755",XHOST.tmpfs+"/xscrubber_AUIDerrors");
      }
    }
    if(VERBOSE) acerr << __AFLOW_FUNC__ << " fixes=" << fixes << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ << " END" << endl;
    return fixes;
  }
}

// will be moved near LIB2AUID
namespace aflowlib {
  bool LIB2AUID(string entry,bool TEST,bool _VERBOSE) {
    if(_VERBOSE){;} //CO20190906 - keep _VERBOSE busy
    bool VERBOSE=FALSE;// _VERBOSE;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " BEGIN" << endl;
    string _entry=entry,directory_LIB,directory_RAW,directory_WEB;
    aurostd::StringSubst(_entry,"/aflow.in","");
    aurostd::StringSubst(_entry,"/"+_AFLOWIN_AEL_DEFAULT_,"");
    aurostd::StringSubst(_entry,"/"+_AFLOWIN_AGL_DEFAULT_,"");
    aurostd::StringSubst(_entry,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"");
    aurostd::StringSubst(_entry,"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON,"");
    aurostd::StringSubst(_entry,"RAW/","LIB/");
    aurostd::StringSubst(_entry,"WEB/","LIB/");
    directory_LIB=_entry;
    directory_RAW=_entry;aurostd::StringSubst(directory_RAW,"LIB/","RAW/");
    directory_WEB=_entry;aurostd::StringSubst(directory_WEB,"LIB/","WEB/");
    // acout << __AFLOW_FUNC__ + " entry=" << entry << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_LIB=" << directory_LIB << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_RAW=" << directory_RAW << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_WEB=" << directory_WEB << endl;

    // if(aurostd::FileExist(directory_LIB)) {
    //   // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_LIB << endl;
    // } else {
    //   // acout << __AFLOW_FUNC__ + " MISSING = " << directory_LIB << endl;
    // }
    // if(aurostd::FileExist(directory_RAW)) {
    //   // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_RAW << endl;
    // } else {
    //   // acout << __AFLOW_FUNC__ + " MISSING = " << directory_RAW << endl;
    // }
    // if(aurostd::FileExist(directory_WEB)) {
    //   // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_WEB << endl;
    // } else {
    //   // acout << __AFLOW_FUNC__ + " MISSING = " << directory_WEB << endl;
    // }

    string directory_old_LIB_AUID,directory_old_RAW_AUID,directory_old_WEB_AUID;
    string directory_new_LIB_AUID,directory_new_RAW_AUID,directory_new_WEB_AUID,directory_new_AUID;
    if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
      _aflowlib_entry entry_tmp(string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT));
      string auid=entry_tmp.auid;
      if(auid.size()!=22) {
        string message = "error on size of auid=";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      directory_old_LIB_AUID=init::AFLOW_Projects_Directories("AUID")+"/LIB";
      directory_old_RAW_AUID=init::AFLOW_Projects_Directories("AUID")+"/RAW";
      directory_old_WEB_AUID=init::AFLOW_Projects_Directories("AUID")+"/WEB";
      directory_new_LIB_AUID=init::AFLOW_Projects_Directories("AUID");
      directory_new_RAW_AUID=init::AFLOW_Projects_Directories("AUID");
      directory_new_WEB_AUID=init::AFLOW_Projects_Directories("AUID");
      directory_new_AUID=init::AFLOW_Projects_Directories("AUID");
      for(uint i=0;i<entry_tmp.vauid.size();i++) {
        directory_old_LIB_AUID+="/"+entry_tmp.vauid.at(i);
        directory_old_RAW_AUID+="/"+entry_tmp.vauid.at(i);
        directory_old_WEB_AUID+="/"+entry_tmp.vauid.at(i);
        directory_new_LIB_AUID+="/"+entry_tmp.vauid.at(i);
        directory_new_RAW_AUID+="/"+entry_tmp.vauid.at(i);
        directory_new_WEB_AUID+="/"+entry_tmp.vauid.at(i);
        directory_new_AUID+="/"+entry_tmp.vauid.at(i);
      }
      directory_new_LIB_AUID+="/LIB";
      directory_new_RAW_AUID+="/RAW";
      directory_new_WEB_AUID+="/WEB";

      if(!TEST) {
        // [OBSOLETE] if(aurostd::FileExist(directory_old_LIB_AUID)) {
        // [OBSOLETE]   // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_old_LIB_AUID << endl;
        // [OBSOLETE] } else {
        // [OBSOLETE]   acout << __AFLOW_FUNC__ + " MISSING = " << directory_old_LIB_AUID << endl;
        // [OBSOLETE] }
        // [OBSOLETE] if(aurostd::FileExist(directory_old_RAW_AUID)) {
        // [OBSOLETE]   // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_old_RAW_AUID << endl;
        // [OBSOLETE] } else {
        // [OBSOLETE]   acout << __AFLOW_FUNC__ + " MISSING = " << directory_old_RAW_AUID << endl;
        // [OBSOLETE] }
        // [OBSOLETE] if(aurostd::FileExist(directory_old_WEB_AUID)) {
        // [OBSOLETE]   // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_old_WEB_AUID << endl;
        // [OBSOLETE] } else {
        // [OBSOLETE]   acout << __AFLOW_FUNC__ + " MISSING = " << directory_old_WEB_AUID << endl;
        // [OBSOLETE] }
        if(aurostd::FileExist(directory_new_LIB_AUID)) {
          // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_new_LIB_AUID << endl;
        } else {
          acout << __AFLOW_FUNC__ + " MISSING = " << directory_new_LIB_AUID << endl;
        }
        if(aurostd::FileExist(directory_new_RAW_AUID)) {
          // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_new_RAW_AUID << endl;
        } else {
          acout << __AFLOW_FUNC__ + " MISSING = " << directory_new_RAW_AUID << endl;
        }
        if(aurostd::FileExist(directory_new_WEB_AUID)) {
          // acout << __AFLOW_FUNC__ + " EXIST   = " << directory_new_WEB_AUID << endl;
        } else {
          acout << __AFLOW_FUNC__ + " MISSING = " << directory_new_WEB_AUID << endl;
        }
      }

      if(aurostd::FileExist(directory_RAW) && !aurostd::FileExist(directory_new_RAW_AUID)) {
        //	acout << __AFLOW_FUNC__ + ": directory_AUID_RAW=" << directory_new_RAW_AUID << " -> " << directory_RAW << endl;
        if(TEST) {
          return true;
        } else {
          acout << __AFLOW_FUNC__ + ": linking file AUID_LIB->LIB: " << directory_new_LIB_AUID << " -> " << directory_LIB << endl; acout.flush();
          aurostd::DirectoryMake(directory_new_AUID);
          aurostd::LinkFile(directory_LIB,directory_new_LIB_AUID);         // LINK
          //	aurostd::execute(XHOST.command("beep")+" -f 2500 -l 1");
        }
      }
      if(aurostd::FileExist(directory_RAW) && !aurostd::FileExist(directory_new_RAW_AUID)) {
        //	acout << __AFLOW_FUNC__ + ": directory_AUID_RAW=" << directory_new_RAW_AUID << " -> " << directory_RAW << endl;
        if(TEST) {
          return true;
        } else {
          acout << __AFLOW_FUNC__ + ": linking file AUID_RAW->RAW: " << directory_new_RAW_AUID << " -> " << directory_RAW << endl; acout.flush();
          aurostd::DirectoryMake(directory_new_AUID);
          aurostd::LinkFile(directory_RAW,directory_new_RAW_AUID);         // LINK
          //	aurostd::execute(XHOST.command("beep")+" -f 2600 -l 1");
        }
      }
      if(aurostd::FileExist(directory_WEB) && !aurostd::FileExist(directory_new_WEB_AUID)) {
        //	acout << __AFLOW_FUNC__ + ": directory_AUID_WEB=" << directory_new_WEB_AUID << " -> " << directory_WEB << endl;
        if(TEST) {
          return true;
        } else {
          acout << __AFLOW_FUNC__ + ": linking file AUID_WEB->WEB: " << directory_new_WEB_AUID << " -> " << directory_WEB << endl; acout.flush();
          aurostd::DirectoryMake(directory_new_AUID);
          aurostd::LinkFile(directory_WEB,directory_new_WEB_AUID);         // LINK
          //	aurostd::execute(XHOST.command("beep")+" -f 2500 -l 1");
        }
      }
      if(!aurostd::FileExist(directory_WEB) && aurostd::FileExist(directory_RAW) && !aurostd::FileExist(directory_new_WEB_AUID)) { // no WEB so point to RAW
        //	acout << __AFLOW_FUNC__ + ": directory_AUID_WEB=" << directory_new_WEB_AUID << " -> " << directory_RAW << endl;
        if(TEST) {
          return true;
        } else {
          acout << __AFLOW_FUNC__ + ": linking file AUID_WEB->RAW: " << directory_new_WEB_AUID << " -> " << directory_RAW << endl; acout.flush();
          aurostd::DirectoryMake(directory_new_AUID);
          aurostd::LinkFile(directory_RAW,directory_new_WEB_AUID);         // LINK
          //	aurostd::execute(XHOST.command("beep")+" -f 2500 -l 1");
        }
      }

      //      directory_old_LIB_AUID=init::AFLOW_Projects_Directories("AUID")+"/LIB/"+auid.substr(0,8); for(uint i=8;i<=20;i+=2) directory_old_LIB_AUID+="/"+auid.substr(i,2);  // splitting aflow:ab/cd..
    }
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_old_LIB_AUID=" << directory_old_LIB_AUID << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_old_RAW_AUID=" << directory_old_RAW_AUID << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_old_WEB_AUID=" << directory_old_WEB_AUID << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_new_LIB_AUID=" << directory_new_LIB_AUID << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_new_RAW_AUID=" << directory_new_RAW_AUID << endl;
    if(VERBOSE) acerr << __AFLOW_FUNC__ + " directory_new_WEB_AUID=" << directory_new_WEB_AUID << endl;

    if(VERBOSE) acerr << __AFLOW_FUNC__ + " END" << endl;
    return false;
  }
}

#endif // _AFLOWLIB_LIBRARIES_SCRUBBER_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
