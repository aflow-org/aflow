// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AUROSTD_ARGV_CPP_
#define _AUROSTD_ARGV_CPP_
//#include "aurostd.h"

#define VERBOSE_ARGV false //DX20200907

namespace aurostd {  // namespace aurostd
  string attach(const string& s1) { return s1;}
  string attach(const string& s1,const string& s2) { return s1+"|"+s2;}
  string attach(const string& s1,const string& s2,const string& s3) { return s1+"|"+s2+"|"+s3;}
  string attach(const string& s1,const string& s2,const string& s3,const string& s4) { return s1+"|"+s2+"|"+s3+"|"+s4;}
  string attach(const string& s1,const string& s2,const string& s3,const string& s4,const string& s5) { return s1+"|"+s2+"|"+s3+"|"+s4+"|"+s5;}
}

// ***************************************************************************
// Function get_arguments_from_input
// ***************************************************************************
// this function creates the argv vector as vector<string> easier to handle than the char **argv
namespace aurostd {  // namespace aurostd
  vector<string> get_arguments_from_input(int _argc,char **_argv) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    vector<string> out_argv;
    // [OBSOLETE]    for(int i=0;i<_argc;i++) out_argv.push_back(_argv[i]);
    for(int i=0;i<_argc;i++) {
      string argi(_argv[i]);
      if(argi=="-np" || argi=="-npmax") argi=string("-")+argi;
      if(argi=="-f" || argi=="-F") argi=string("-")+argi;
      if(argi=="-d" || argi=="-D") argi=string("-")+argi;
      if(argi.size()>=2) if(argi.at(0)=='-' && argi.at(1)=='n') argi=string("-")+argi;  // for -np
      //   if(argi.size()>=2) if(argi.at(0)=='-' && argi.at(1)=='m') argi=string("-")+argi;  // for -machine
      if(argi.size()>=2) if(argi.at(0)=='-' && argi.at(1)=='f') argi=string("-")+argi;  // for -f
      if(argi.size()>=2) if(argi.at(0)=='-' && argi.at(1)=='F') argi=string("-")+argi;  // for -F
      if(argi.size()>=2) if(argi.at(0)=='-' && argi.at(1)=='d') argi=string("-")+argi;  // for -d
      if(argi.size()>=2) if(argi.at(0)=='-' && argi.at(1)=='D') argi=string("-")+argi;  // for -D

      if(argi=="--machine") argi+=string("=");                  // forcing "=" after machine !
      if(argi=="--machine_name") argi+=string("=");             // forcing "=" after machine_name ! //HE20220309
      if(argi=="--aflowlib") argi+=string("=");                 // forcing "=" after aflowlib !
      if(argi=="--np") argi+=string("=");                       // forcing "=" after np !
      if(argi.at(argi.size()-1)=='=' && i<_argc-1) {argi+=string(_argv[i+1]);i++;}  // fixing space after "= "
      out_argv.push_back(argi);
    }
    if(VERBOSE) for(uint i=0;i<out_argv.size();i++) cerr << "out_argv.at(" << i << ")=" << out_argv.at(i) << endl;
    return out_argv;
  }
}

// ***************************************************************************
// Function args2flag without/with  commands
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  bool args2flag(const vector<string>& argv,const string& s0) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint j=0;j<tokens.size();j++)
      for(uint i=0;i<argv.size();i++)
        if(argv.at(i)==tokens.at(j)) return TRUE;
    return FALSE;
  }

  bool args2flag(const vector<string>& argv,std::vector<string>& cmds,const string& s0) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint j=0;j<tokens.size();j++)
      cmds.push_back(tokens.at(j));
    for(uint i=0;i<argv.size();i++)
      for(uint j=0;j<tokens.size();j++)     
        if(argv.at(i)==tokens.at(j)) return TRUE;
    return FALSE;
  }
}

// ***************************************************************************
// args2utype of utype type
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template<class utype> utype args2utype(const vector<string>& argv,const string& s0,utype def_out) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    utype out=def_out;
    for(uint i=1;i<argv.size()-1;i++)
      for(uint j=0;j<tokens.size();j++)     
        if(argv.at(i)==tokens.at(j)) { if(_isfloat(out)) out=(utype) atof(argv.at(i+1).c_str()); else out=(utype) atoi(argv.at(i+1).c_str());} // OLD
    return (utype) out;
  }
}

// ***************************************************************************
// Function get xvector from input
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype> args2xvectorutype(const vector<string>& argv,const string& s0, const xvector<utype>& def_out) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    xvector<utype> out(def_out.lrows,def_out.urows);
    out=def_out;
    for(uint i=0;i<argv.size();i++) {
      if(i+out.rows<argv.size()) {
        for(uint j=0;j<tokens.size();j++) {
          if(argv.at(i)==tokens.at(j)) {
            for(int k=0;k<out.rows;k++) {
              if(_isfloat(out(1))) {
                out(k+out.lrows)=(utype) atof(argv.at(i+k+1).c_str());
              } else {
                out(k+out.lrows)=(utype) atoi(argv.at(i+k+1).c_str());
              }
            }
          }
        }
      }
    }
    return out;
  }
  template<class utype> xvector<utype> args2xvectorutype(const vector<string>& argv,const string& s0,int dim) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    xvector<utype> out(1,dim);
    for(uint i=0;i<argv.size();i++) {
      if(i+out.rows<argv.size()) {
        for(uint j=0;j<tokens.size();j++) {
          if(argv.at(i)==tokens.at(j)) {
            for(int k=0;k<out.rows;k++) {
              if(_isfloat(out(1))) {
                out(k+out.lrows)=(utype) atof(argv.at(i+k+1).c_str());
              } else {
                out(k+out.lrows)=(utype) atoi(argv.at(i+k+1).c_str());
              }
            }
          }
        }
      }
    }
    return out; // something phony to keep t used
  }
}
// ***************************************************************************
// Function get vector/deque from input
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template<class utype> vector<utype>
    args2vectorutype(const vector<string>& argv,const string& s0) {
      bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
      string s=aurostd::RemoveWhiteSpaces(s0);
      vector<string> tokens;
      aurostd::string2tokens(s,tokens,"|"); 
      if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
      vector<utype> out;
      for(uint i=0;i<argv.size();i++) {
        for(uint j=0;j<tokens.size();j++) {
          if(argv.at(i)==tokens.at(j)) {
            for(uint k=0;k<argv.size();k++) {
              if(_isfloat(out.at(0))) {
                out.push_back((utype) atof(argv.at(i+k+1).c_str()));
              } else {
                out.push_back((utype) atoi(argv.at(i+k+1).c_str()));
              }
            }
          }
        }
      }  
      return out; 
    }

  template<class utype> deque<utype>
    args2dequeutype(const deque<string>& argv,const string& s0) {
      bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
      string s=aurostd::RemoveWhiteSpaces(s0);
      deque<string> tokens;
      aurostd::string2tokens(s,tokens,"|"); 
      if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
      deque<utype> out;
      for(uint i=0;i<argv.size();i++) {
        for(uint j=0;j<tokens.size();j++) {
          if(argv.at(i)==tokens.at(j)) {
            for(uint k=0;k<argv.size();k++) {
              if(_isfloat(out.at(0))) {
                out.push_back((utype) atof(argv.at(i+k+1).c_str()));
              } else {
                out.push_back((utype) atoi(argv.at(i+k+1).c_str()));
              }
            }
          }
        }
      }
      return out;
    }
}

// ***************************************************************************
// Functions args2string functions
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  string args2string(const vector<string>& argv,const string& s0,const string& s_def) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint i=0;i<argv.size()-1;i++)
      for(uint j=0;j<tokens.size();j++) 
        if(argv.at(i)==tokens.at(j)) return argv.at(i+1);
    return s_def;
  }

  string args2string(const vector<string>& argv,vector<string>& cmds,const string& s0,const string& s_def) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint j=0;j<tokens.size();j++) 
      cmds.push_back(tokens.at(j));
    for(uint i=0;i<argv.size()-1;i++)
      for(uint j=0;j<tokens.size();j++) 
        if(argv.at(i)==tokens.at(j)) return argv.at(i+1);
    return s_def;
  }
}

// ***************************************************************************
// Functions args2vectorstring
// ***************************************************************************
namespace aurostd {
  vector<string> string2vstring(const string& str_in) {
    vector<string> str_out;
    str_out.push_back(str_in);
    return str_out;
  }

  vector<string> args2vectorstring(const vector<string>& argv,const string& s0,const string& s_def) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint i=0;i<argv.size()-1;i++)
      for(uint j=0;j<tokens.size();j++) 
        if(argv.at(i)==tokens.at(j))
          return vector<string>(argv.begin()+(i+1),argv.end());
    return string2vstring(s_def);
  }
}

// ***************************************************************************
// Functions get_itemized_vector_string stuff
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) {// =":")
    if(aurostd::substring2bool(s0,"|")) {
      string message = "not ported to \"|\"";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    uint icount=0;
    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
    if(aurostd::args2attachedflag(argv,s0equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
    if(tokens.size()==0) return FALSE;
    if(icount==0) return FALSE;
    return TRUE;
  }
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) {// =":")
    uint icount=0;
    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
    string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
    tokens.clear();
    if(!tokens.size() && aurostd::args2attachedflag(argv,s0equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
    if(!tokens.size() && aurostd::args2attachedflag(argv,s1equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
    if(!tokens.size() && aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
    if(!tokens.size() && aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
    if(tokens.size()==0) return FALSE;
    if(icount==0) return FALSE;
    return TRUE;
  }
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) {// =":")
    uint icount=0;
    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
    string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
    string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
    if(aurostd::args2attachedflag(argv,s0equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2attachedflag(argv,s1equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2attachedflag(argv,s2equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
    if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
    if(aurostd::args2flag(argv,s2neq)) {aurostd::args2string(argv,s2neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s2neq,EMPTY_WORDING);}
    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
    if(tokens.size()==0) return FALSE;
    if(icount==0) return FALSE;
    return TRUE;
  }
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) {// =":")
    uint icount=0;
    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
    string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
    string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
    string s3neq=s3,s3equ;aurostd::StringSubst(s3neq,"=","");s3equ=s3neq+"=";
    if(aurostd::args2attachedflag(argv,s0equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2attachedflag(argv,s1equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2attachedflag(argv,s2equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2attachedflag(argv,s3equ)) {icount+=aurostd::string2tokens(aurostd::args2attachedstring(argv,s3equ,EMPTY_WORDING),tokens,delimiter);}
    if(aurostd::args2flag(argv,s0neq)) {aurostd::args2string(argv,s0neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s0neq,EMPTY_WORDING);}
    if(aurostd::args2flag(argv,s1neq)) {aurostd::args2string(argv,s1neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s1neq,EMPTY_WORDING);}
    if(aurostd::args2flag(argv,s2neq)) {aurostd::args2string(argv,s2neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s2neq,EMPTY_WORDING);}
    if(aurostd::args2flag(argv,s3neq)) {aurostd::args2string(argv,s3neq,EMPTY_WORDING);tokens=aurostd::args2vectorstring(argv,s3neq,EMPTY_WORDING);}
    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
    if(tokens.size()==0) return FALSE;
    if(icount==0) return FALSE;
    return TRUE;
  }
}

//[CO20200624 - OBSOLETE]// ***************************************************************************
//[CO20200624 - OBSOLETE]// Function getproto_itemized_vector_string_from_input
//[CO20200624 - OBSOLETE]// ***************************************************************************
//[CO20200624 - OBSOLETE]namespace aurostd {
//[CO20200624 - OBSOLETE]  uint LOCAL_XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX) {
//[CO20200624 - OBSOLETE]    string alloy=alloy_in;
//[CO20200624 - OBSOLETE]    alloy=aurostd::RemoveNumbers(alloy);              // remove composition
//[CO20200624 - OBSOLETE]    species_ppX.clear();
//[CO20200624 - OBSOLETE]    for(uint i=0;i<alloy.length();i++) {
//[CO20200624 - OBSOLETE]      if(alloy[i]>='A' && alloy[i]<='Z') species_ppX.push_back("");
//[CO20200624 - OBSOLETE]      species_ppX.at(species_ppX.size()-1)+=alloy[i];
//[CO20200624 - OBSOLETE]    }
//[CO20200624 - OBSOLETE]    for(uint i=0;i<species_ppX.size();i++)
//[CO20200624 - OBSOLETE]      species_ppX.at(i)=aurostd::CleanStringASCII(species_ppX.at(i));
//[CO20200624 - OBSOLETE]    return species_ppX.size();
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]  
//[CO20200624 - OBSOLETE]  bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) {// =":") {
//[CO20200624 - OBSOLETE]    string ss;tokens.clear();vector<string> stokens;
//[CO20200624 - OBSOLETE]    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//[CO20200624 - OBSOLETE]    if(ss!="") {
//[CO20200624 - OBSOLETE]      if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,tokens,delimiter);
//[CO20200624 - OBSOLETE]      } else {
//[CO20200624 - OBSOLETE]	aurostd::string2tokens(ss,stokens,"/");
//[CO20200624 - OBSOLETE]	tokens.push_back(stokens.at(1));
//[CO20200624 - OBSOLETE]	aurostd::LOCAL_XATOM_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
//[CO20200624 - OBSOLETE]	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//[CO20200624 - OBSOLETE]      }
//[CO20200624 - OBSOLETE]    }
//[CO20200624 - OBSOLETE]    uint icount=0;
//[CO20200624 - OBSOLETE]    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
//[CO20200624 - OBSOLETE]    if(tokens.size()==0) return FALSE;
//[CO20200624 - OBSOLETE]    return TRUE;
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]  bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) {// =":") {
//[CO20200624 - OBSOLETE]     string ss;tokens.clear();vector<string> stokens;
//[CO20200624 - OBSOLETE]    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//[CO20200624 - OBSOLETE]    string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//[CO20200624 - OBSOLETE]    if(ss!="") {
//[CO20200624 - OBSOLETE]      if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,tokens,delimiter);
//[CO20200624 - OBSOLETE]      } else {
//[CO20200624 - OBSOLETE]	aurostd::string2tokens(ss,stokens,"/");
//[CO20200624 - OBSOLETE]	tokens.push_back(stokens.at(1));
//[CO20200624 - OBSOLETE]	aurostd::LOCAL_XATOM_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
//[CO20200624 - OBSOLETE]	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//[CO20200624 - OBSOLETE]      }
//[CO20200624 - OBSOLETE]    }
//[CO20200624 - OBSOLETE]    uint icount=0;
//[CO20200624 - OBSOLETE]    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
//[CO20200624 - OBSOLETE]    if(tokens.size()==0) return FALSE;
//[CO20200624 - OBSOLETE]    return TRUE;
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]  bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) {// =":") {
//[CO20200624 - OBSOLETE]    string ss;tokens.clear();vector<string> stokens;
//[CO20200624 - OBSOLETE]    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//[CO20200624 - OBSOLETE]    string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//[CO20200624 - OBSOLETE]    string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s2equ)) ss=aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s2neq)) ss=aurostd::args2string(argv,s2neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//[CO20200624 - OBSOLETE]    if(ss!="") {
//[CO20200624 - OBSOLETE]      if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,s2,tokens,delimiter);
//[CO20200624 - OBSOLETE]      } else {
//[CO20200624 - OBSOLETE]	aurostd::string2tokens(ss,stokens,"/");
//[CO20200624 - OBSOLETE]	tokens.push_back(stokens.at(1));
//[CO20200624 - OBSOLETE]	aurostd::LOCAL_XATOM_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
//[CO20200624 - OBSOLETE]	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//[CO20200624 - OBSOLETE]      }
//[CO20200624 - OBSOLETE]    }
//[CO20200624 - OBSOLETE]    uint icount=0;
//[CO20200624 - OBSOLETE]    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
//[CO20200624 - OBSOLETE]    if(tokens.size()==0) return FALSE;
//[CO20200624 - OBSOLETE]    return TRUE;
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]  bool getproto_itemized_vector_string_from_input(vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) {// =":") {
//[CO20200624 - OBSOLETE]    string ss;tokens.clear();vector<string> stokens;
//[CO20200624 - OBSOLETE]    string s0neq=s0,s0equ;aurostd::StringSubst(s0neq,"=","");s0equ=s0neq+"=";
//[CO20200624 - OBSOLETE]    string s1neq=s1,s1equ;aurostd::StringSubst(s1neq,"=","");s1equ=s1neq+"=";
//[CO20200624 - OBSOLETE]    string s2neq=s2,s2equ;aurostd::StringSubst(s2neq,"=","");s2equ=s2neq+"=";
//[CO20200624 - OBSOLETE]    string s3neq=s3,s3equ;aurostd::StringSubst(s3neq,"=","");s3equ=s3neq+"=";
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s0equ)) ss=aurostd::args2attachedstring(argv,s0equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s0neq)) ss=aurostd::args2string(argv,s0neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s1equ)) ss=aurostd::args2attachedstring(argv,s1equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s1neq)) ss=aurostd::args2string(argv,s1neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s2equ)) ss=aurostd::args2attachedstring(argv,s2equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s2neq)) ss=aurostd::args2string(argv,s2neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2attachedflag(argv,s3equ)) ss=aurostd::args2attachedstring(argv,s3equ,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::args2flag(argv,s3neq)) ss=aurostd::args2string(argv,s3neq,EMPTY_WORDING);
//[CO20200624 - OBSOLETE]    if(aurostd::substring2bool(ss,"./")) aurostd::StringSubst(ss,"./","");
//[CO20200624 - OBSOLETE]    if(ss!="") {
//[CO20200624 - OBSOLETE]      if(!aurostd::substring2bool(ss,"/")) { return get_itemized_vector_string_from_input(argv,s0,s1,s2,s3,tokens,delimiter);
//[CO20200624 - OBSOLETE]      } else {
//[CO20200624 - OBSOLETE]	aurostd::string2tokens(ss,stokens,"/");
//[CO20200624 - OBSOLETE]	tokens.push_back(stokens.at(1));
//[CO20200624 - OBSOLETE]	aurostd::LOCAL_XATOM_SplitAlloyPseudoPotentials(stokens.at(0),stokens);
//[CO20200624 - OBSOLETE]	for(uint i=0;i<stokens.size();i++) tokens.push_back(stokens.at(i));
//[CO20200624 - OBSOLETE]      }
//[CO20200624 - OBSOLETE]    }
//[CO20200624 - OBSOLETE]    uint icount=0;
//[CO20200624 - OBSOLETE]    if(tokens.size()==1 && aurostd::substring2bool(tokens.at(0),delimiter)) {s0equ=tokens.at(0);icount+=aurostd::string2tokens(s0equ,tokens,delimiter);}
//[CO20200624 - OBSOLETE]    if(tokens.size()==0) return FALSE;
//[CO20200624 - OBSOLETE]    return TRUE;
//[CO20200624 - OBSOLETE]  }
//[CO20200624 - OBSOLETE]}

// ***************************************************************************
// args2attachedflag without/with commands
// ***************************************************************************
namespace aurostd {  // namespace aurostd

  bool args2attachedflag(const vector<string>& argv,const string& s0) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint i=0;i<argv.size();i++)
      for(uint j=0;j<tokens.size();j++)
        if(aurostd::substring2bool(argv.at(i),tokens.at(j))) return TRUE;
    return FALSE;
  }

  bool args2attachedflag(const vector<string>& argv,std::vector<string>& cmds,const string& s0) {
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;
    for(uint j=0;j<tokens.size();j++)
      cmds.push_back(tokens.at(j));
    for(uint i=0;i<argv.size();i++)
      for(uint j=0;j<tokens.size();j++)
        if(aurostd::substring2bool(argv.at(i),tokens.at(j))) return TRUE;
    return FALSE;
  }
}

// ***************************************************************************
// args2attachedstring(argv,"-xxx=",abc) returns the content of xxx= or abc
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  string args2attachedstring(const vector<string>& argv,const string& s0,string s_def) { // string=""
    bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s=aurostd::RemoveWhiteSpaces(s0),output="";
    vector<string> tokens;
    aurostd::string2tokens(s,tokens,"|"); 
    if(VERBOSE) {cerr << "argv.size()=" << argv.size() << endl;for(uint j=0;j<argv.size();j++) cerr << "[" << argv.at(j) << "]" << endl;}
    if(VERBOSE) {cerr << "tokens.size()=" << tokens.size() << endl;for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;}
    for(uint i=1;i<argv.size();i++)
      for(uint j=0;j<tokens.size();j++)
        if(argv.at(i).find(tokens.at(j))!=string::npos) {
          output=argv.at(i).substr(argv.at(i).find(tokens.at(j))+tokens.at(j).length()); 
          if(output.length()>0) return output;
        }
    return s_def;
  }
  // [OBSOLETE]  string args2attachedstring(const vector<string>& argv,const string& s0,const string& s_def) {
  // [OBSOLETE]    string s;
  // [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
  // [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s0)+s0.length()); if(s.length()>0) return s;}
  // [OBSOLETE]    }
  // [OBSOLETE]    return s_def;
  // [OBSOLETE]  }
  // [OBSOLETE]  string args2attachedstring(const vector<string>& argv,const string& s0,const string& s1,const string& s_def) {
  // [OBSOLETE]    string s;
  // [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
  // [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s0)+s0.length()); if(s.length()>0) return s;}
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s1)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s1.length());
  // [OBSOLETE]      if(argv.at(i).find(s1)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s1)+s1.length()); if(s.length()>0) return s;}
  // [OBSOLETE]    }
  // [OBSOLETE]    return s_def;
  // [OBSOLETE]  }
  // [OBSOLETE]  string args2attachedstring(const vector<string>& argv,const string& s0,const string& s1,const string& s2,const string& s_def) {
  // [OBSOLETE]    string s;
  // [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
  // [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s0)+s0.length()); if(s.length()>0) return s;}
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s1)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s1.length());
  // [OBSOLETE]      if(argv.at(i).find(s1)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s1)+s1.length()); if(s.length()>0) return s;}
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s2)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s2.length());
  // [OBSOLETE]      if(argv.at(i).find(s2)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s2)+s2.length()); if(s.length()>0) return s;}
  // [OBSOLETE]    }
  // [OBSOLETE]    return s_def;
  // [OBSOLETE]  }
  // [OBSOLETE]  string args2attachedstring(const vector<string>& argv,const string& s0,const string& s1,const string& s2,const string& s3,const string& s_def) {
  // [OBSOLETE]    string s;
  // [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
  // [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s0)+s0.length()); if(s.length()>0) return s;}
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s1)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s1.length());
  // [OBSOLETE]      if(argv.at(i).find(s1)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s1)+s1.length()); if(s.length()>0) return s;}
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s2)!=string::npos) return argv.at(i).substr(argv.at(i).find(s2)+s2.length());
  // [OBSOLETE]      if(argv.at(i).find(s2)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s2)+s2.length()); if(s.length()>0) return s;}
  // [OBSOLETE]      // [OLDSTYLE] if(argv.at(i).find(s3)!=string::npos) return argv.at(i).substr(argv.at(i).find(s3)+s3.length());
  // [OBSOLETE]      if(argv.at(i).find(s3)!=string::npos) {s=argv.at(i).substr(argv.at(i).find(s3)+s3.length()); if(s.length()>0) return s;}
  // [OBSOLETE]    }
  // [OBSOLETE]    return s_def;
  // [OBSOLETE]  }

  template<typename string>
    string args2attachedutype(const vector<string>& argv,const string& s0,const string& s_def) {
      bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
      string s=aurostd::RemoveWhiteSpaces(s0),output="";
      vector<string> tokens;
      aurostd::string2tokens(s,tokens,"|"); 
      if(VERBOSE) {cerr << "argv.size()=" << argv.size() << endl;for(uint j=0;j<argv.size();j++) cerr << "[" << argv.at(j) << "]" << endl;}
      if(VERBOSE) {cerr << "tokens.size()=" << tokens.size() << endl;for(uint j=0;j<tokens.size();j++) cerr << "[" << tokens.at(j) << "]" << endl;}
      for(uint i=1;i<argv.size();i++)
        for(uint j=0;j<tokens.size();j++)
          if(argv.at(i).find(tokens.at(j))!=string::npos) return argv.at(i).substr(argv.at(i).find(tokens.at(j))+tokens.at(j).length());
      return s_def;
    }
}
// [OBSOLETE]    template<typename string>
// [OBSOLETE]    string args2attachedutype(const vector<string>& argv,const string& s0,const string& s_def) {
// [OBSOLETE]      vector<string> tokens;
// [OBSOLETE]      for(uint i=0;i<argv.size();i++)
// [OBSOLETE]        if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
// [OBSOLETE]      return s_def;
// [OBSOLETE]    }
// [OBSOLETE]  template<typename string>
// [OBSOLETE]  string args2attachedutype(const vector<string>& argv,const string& s0,const string& s1,const string& s_def) {
// [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
// [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
// [OBSOLETE]      if(argv.at(i).find(s1)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s1.length());
// [OBSOLETE]    }
// [OBSOLETE]    return s_def;
// [OBSOLETE]  }
// [OBSOLETE]  template<typename string>
// [OBSOLETE]  string args2attachedutype(const vector<string>& argv,const string& s0,const string& s1,const string& s2,const string& s_def) {
// [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
// [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
// [OBSOLETE]      if(argv.at(i).find(s1)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s1.length());
// [OBSOLETE]      if(argv.at(i).find(s2)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s2.length());
// [OBSOLETE]    }
// [OBSOLETE]    return s_def;
// [OBSOLETE]  }
// [OBSOLETE]  template<typename string>
// [OBSOLETE]  string args2attachedutype(const vector<string>& argv,const string& s0,const string& s1,const string& s2,const string& s3,const string& s_def) {
// [OBSOLETE]    for(uint i=0;i<argv.size();i++) {
// [OBSOLETE]      if(argv.at(i).find(s0)!=string::npos) return argv.at(i).substr(argv.at(i).find(s0)+s0.length());
// [OBSOLETE]      if(argv.at(i).find(s1)!=string::npos) return argv.at(i).substr(argv.at(i).find(s1)+s1.length());
// [OBSOLETE]      if(argv.at(i).find(s2)!=string::npos) return argv.at(i).substr(argv.at(i).find(s2)+s2.length());
// [OBSOLETE]      if(argv.at(i).find(s3)!=string::npos) return argv.at(i).substr(argv.at(i).find(s3)+s3.length());
// [OBSOLETE]    }
// [OBSOLETE]    return s_def;
// [OBSOLETE]  }
// [OBSOLETE]}

// ***************************************************************************
// args2attachedutype
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  template<typename utype>
    utype args2attachedutype(const vector<string>& argv,const string& str1,const utype& utype_default) {
      bool VERBOSE=(FALSE || VERBOSE_ARGV); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
      vector<string> tokens1;
      aurostd::string2tokens(aurostd::RemoveWhiteSpaces(str1),tokens1,"|"); 
      if(VERBOSE) {cerr << "argv.size()=" << argv.size() << endl;for(uint j=0;j<argv.size();j++) cerr << "[" << argv.at(j) << "]" << endl;}
      if(VERBOSE) {cerr << "tokens1.size()=" << tokens1.size() << endl;for(uint j=0;j<tokens1.size();j++) cerr << "[" << tokens1.at(j) << "]" << endl;}
      utype out=utype_default;
      for(uint j=0;j<tokens1.size();j++) {
        string s1=tokens1.at(j),s1eq,s1neq;     //   s1=aurostd::RemoveSubString(s1,"-");s1=aurostd::RemoveSubString(s1,"-");
        s1=aurostd::RemoveSubString(s1,"=");
        s1eq=s1+"=";s1neq=s1;//cerr << s1eq << " " << s1neq << endl;
        if(aurostd::args2flag(argv,s1neq) || aurostd::args2attachedflag(argv,s1eq)) {
          if(aurostd::args2flag(argv,s1neq)) out=aurostd::args2utype(argv,s1neq,(utype) out);
          if(aurostd::args2attachedflag(argv,s1eq)) {
            vector<string> tokens2;
            get_itemized_vector_string_from_input(argv,s1neq,tokens2,",");
            if(tokens2.size()>0) out=aurostd::string2utype<utype>(tokens2.at(0));
          }
        }
      }
      return (utype) out;
    }
}

// [OBSOLETE]  template<typename utype>
// [OBSOLETE]  utype args2attachedutype(const vector<string>& argv,const string& str1,const string& str2,const utype& utype_default) {
// [OBSOLETE]    string s1=str1,s1eq,s1neq;     //   s1=aurostd::RemoveSubString(s1,"-");s1=aurostd::RemoveSubString(s1,"-");
// [OBSOLETE]    s1=aurostd::RemoveSubString(s1,"=");
// [OBSOLETE]    s1eq=s1+"=";s1neq=s1;//cerr << s1eq << " " << s1neq << endl;
// [OBSOLETE]    string s2=str2,s2eq,s2neq;     //   s2=aurostd::RemoveSubString(s2,"-");s2=aurostd::RemoveSubString(s2,"-");
// [OBSOLETE]    s2=aurostd::RemoveSubString(s2,"=");
// [OBSOLETE]    s2eq=s2+"=";s2neq=s2;//cerr << s2eq << " " << s2neq << endl;
// [OBSOLETE]    utype out=utype_default;
// [OBSOLETE]    if(aurostd::args2flag(argv,aurostd::attach(s1neq,s2neq)) || aurostd::args2attachedflag(argv,aurostd::attach(s1eq,s2eq))) {
// [OBSOLETE]      if(aurostd::args2flag(argv,aurostd::attach(s1neq,s2neq))) out=aurostd::args2utype(argv,aurostd::attach(s1neq,s2neq),(utype) out);
// [OBSOLETE]      if(aurostd::args2attachedflag(argv,aurostd::attach(s1eq,s2eq))) {
// [OBSOLETE]	vector<string> tokens;
// [OBSOLETE]	get_itemized_vector_string_from_input(argv,s1neq,s2neq,tokens,",");
// [OBSOLETE]	if(tokens.size()>0) out=aurostd::string2utype<utype>(tokens.at(0));
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]    return (utype) out;
// [OBSOLETE]  }
// [OBSOLETE]  
// [OBSOLETE]  template<typename utype>
// [OBSOLETE]  utype args2attachedutype(const vector<string>& argv,const string& str1,const string& str2,const string& str3,const utype& utype_default) {
// [OBSOLETE]    string s1=str1,s1eq,s1neq;     //   s1=aurostd::RemoveSubString(s1,"-");s1=aurostd::RemoveSubString(s1,"-");
// [OBSOLETE]    s1=aurostd::RemoveSubString(s1,"=");
// [OBSOLETE]    s1eq=s1+"=";s1neq=s1;//cerr << s1eq << " " << s1neq << endl;
// [OBSOLETE]    string s2=str2,s2eq,s2neq;     //   s2=aurostd::RemoveSubString(s2,"-");s2=aurostd::RemoveSubString(s2,"-");
// [OBSOLETE]    s2=aurostd::RemoveSubString(s2,"=");
// [OBSOLETE]    s2eq=s2+"=";s2neq=s2;//cerr << s2eq << " " << s2neq << endl;
// [OBSOLETE]    string s3=str3,s3eq,s3neq;     //   s3=aurostd::RemoveSubString(s3,"-");s3=aurostd::RemoveSubString(s3,"-");
// [OBSOLETE]    s3=aurostd::RemoveSubString(s3,"=");
// [OBSOLETE]    s3eq=s3+"=";s3neq=s3;//cerr << s3eq << " " << s3neq << endl;
// [OBSOLETE]    utype out=utype_default;
// [OBSOLETE]    if(aurostd::args2flag(argv,aurostd::attach(s1neq,s2neq,s3neq)) || aurostd::args2attachedflag(argv,aurostd::attach(s1eq,s2eq,s3eq))) {
// [OBSOLETE]      if(aurostd::args2flag(argv,aurostd::attach(s1neq,s2neq,s3neq))) out=aurostd::args2utype(argv,aurostd::attach(s1neq,s2neq,s3neq),(utype) out);
// [OBSOLETE]      if(aurostd::args2attachedflag(argv,aurostd::attach(s1eq,s2eq,s3eq))) {
// [OBSOLETE]	vector<string> tokens;
// [OBSOLETE]	get_itemized_vector_string_from_input(argv,s1neq,s2neq,s3neq,tokens,",");
// [OBSOLETE]	if(tokens.size()>0) out=aurostd::string2utype<utype>(tokens.at(0));
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]    return (utype) out;
// [OBSOLETE]  }
// [OBSOLETE]  
// [OBSOLETE]  template<typename utype>
// [OBSOLETE]  utype args2attachedutype(const vector<string>& argv,const string& str1,const string& str2,const string& str3,const string& str4,const utype& utype_default) {
// [OBSOLETE]    string s1=str1,s1eq,s1neq;     //   s1=aurostd::RemoveSubString(s1,"-");s1=aurostd::RemoveSubString(s1,"-");
// [OBSOLETE]    s1=aurostd::RemoveSubString(s1,"=");
// [OBSOLETE]    s1eq=s1+"=";s1neq=s1;//cerr << s1eq << " " << s1neq << endl;
// [OBSOLETE]    string s2=str2,s2eq,s2neq;     //   s2=aurostd::RemoveSubString(s2,"-");s2=aurostd::RemoveSubString(s2,"-");
// [OBSOLETE]    s2=aurostd::RemoveSubString(s2,"=");
// [OBSOLETE]    s2eq=s2+"=";s2neq=s2;//cerr << s2eq << " " << s2neq << endl;
// [OBSOLETE]    string s3=str3,s3eq,s3neq;     //   s3=aurostd::RemoveSubString(s3,"-");s3=aurostd::RemoveSubString(s3,"-");
// [OBSOLETE]    s3=aurostd::RemoveSubString(s3,"=");
// [OBSOLETE]    s3eq=s3+"=";s3neq=s3;//cerr << s3eq << " " << s3neq << endl;
// [OBSOLETE]    string s4=str4,s4eq,s4neq;     //   s4=aurostd::RemoveSubString(s4,"-");s4=aurostd::RemoveSubString(s4,"-");
// [OBSOLETE]    s4=aurostd::RemoveSubString(s4,"=");
// [OBSOLETE]    s4eq=s4+"=";s4neq=s4;//cerr << s4eq << " " << s4neq << endl;
// [OBSOLETE]    utype out=utype_default;
// [OBSOLETE]    if(aurostd::args2flag(argv,aurostd::attach(s1neq,s2neq,s3neq,s4neq)) || aurostd::args2attachedflag(argv,aurostd::attach(s1eq,s2eq,s3eq,s4eq))) {
// [OBSOLETE]      if(aurostd::args2flag(argv,aurostd::attach(s1neq,s2neq,s3neq,s4neq))) out=aurostd::args2utype(argv,aurostd::attach(s1neq,s2neq,s3neq,s4neq),(utype) out);
// [OBSOLETE]      if(aurostd::args2attachedflag(argv,aurostd::attach(s1eq,s2eq,s3eq,s4eq))) {
// [OBSOLETE]	vector<string> tokens;
// [OBSOLETE]	get_itemized_vector_string_from_input(argv,s1neq,s2neq,s3neq,s4neq,tokens,",");
// [OBSOLETE]	if(tokens.size()>0) out=aurostd::string2utype<utype>(tokens.at(0));
// [OBSOLETE]      }
// [OBSOLETE]    }
// [OBSOLETE]    return (utype) out;
// [OBSOLETE]  }
// [OBSOLETE]}

// *******************************************************************************************
// *******************************************************************************************

#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

