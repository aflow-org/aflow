// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017
// streamline schemes SC 2017

#ifndef _AUROSTD_XOPTION_CPP_
#define _AUROSTD_XOPTION_CPP_
//#include "aflow.h"
using std::ostream;

#define VERBOSE_XOPTION false //DX20200907

namespace aurostd {
  // ***************************************************************************

  // constructure
  xoption::xoption() {free();}  //CO20200624 - moved to free()

  // destructor
  xoption::~xoption() {free();} //CO20200624 - moved to free()

  // free
  void xoption::free() {
    keyword=""; //DX20180824 - missing from constructor
    isentry=FALSE; 
    content_string="";
    content_double=0.0;
    content_int=0;
    content_uint=0;
    option=FALSE;
    option_default=FALSE;
    xscheme="";
    vxscheme.clear();
    vxsghost.clear();
    preserved=FALSE;
  }

  // copy fuction
  void xoption::copy(const xoption& b) {
    keyword=b.keyword; //DX20180824 - missing from copy constructor
    isentry=b.isentry;
    content_string=b.content_string;
    content_double=b.content_double;
    content_int=b.content_int;
    content_uint=b.content_uint;
    option=b.option;
    option_default=b.option_default;
    xscheme=b.xscheme;
    vxscheme.clear();for(uint i=0;i<b.vxscheme.size();i++) vxscheme.push_back(b.vxscheme.at(i));
    vxsghost.clear();for(uint i=0;i<b.vxsghost.size();i++) vxsghost.push_back(b.vxsghost.at(i));
    preserved=b.preserved;
  }

  // copy conctructor
  xoption::xoption(const xoption& b) {
    //  free();
    // *this=b;
    copy(b);
  }

  // copy operator b=a
  const xoption& xoption::operator=(const xoption& b) {  // operator=
    if(this!=&b) {
      free();
      copy(b);
    }
    return *this;
  }

  std::ostream& operator<<(std::ostream& oss,const xoption& a) {
    for(uint i=0;i<a.vxscheme.size();i++)
      oss << a.vxscheme.at(i) << (i<a.vxscheme.size()-1?",":"");
    return oss;
  }

  void xoption::clear() {
    //[CO20200624 - creating objects is SLOW]xoption aflow_option_temp;
    //[CO20200624 - creating objects is SLOW]copy(aflow_option_temp);
    free();
  }

  // **************************************************************************
  //void xoption::options2entry(const string& options_FILE,const string& input_keyword,int _option_DEFAULT,const string& xscheme_DEFAULT) //CO20210805 - const&
  void xoption::options2entry(const string& options_FILE_IN,const string& input_keyword_IN,int option_DEFAULT_IN,const string& xscheme_DEFAULT_IN) {
    bool VERBOSE=(FALSE || VERBOSE_XOPTION); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;

    //CO20210909 - BIG BUG HERE
    //the following clear() will reset all of the internal xoption variables
    //if we pass one of these variables into options2entry(), it is reset as well
    //see for example, this construction:
    //vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.options2entry(AflowIn,_STROPT_+"NELM=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.xscheme);
    //xscheme gets cleared before it's set
    //in order to preserve this construction and prevent headaches, make copies of the inputs BEFORE the clear
    string options_FILE=options_FILE_IN;
    string input_keyword=input_keyword_IN;
    int _option_DEFAULT=option_DEFAULT_IN;
    string xscheme_DEFAULT=xscheme_DEFAULT_IN;

    clear();  //CO20210909 - DANGEROUS! see note above

    bool option_DEFAULT=FALSE; 
    if(_option_DEFAULT==0) option_DEFAULT=FALSE; // it is a int.. it might be -1
    if(_option_DEFAULT==1) option_DEFAULT=TRUE; // it is a int.. it might be -1
    isentry=option_DEFAULT;option=option_DEFAULT;content_string=xscheme_DEFAULT;xscheme=xscheme_DEFAULT;preserved=FALSE;   // DEFAULT
    if(VERBOSE){
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " BEGIN " << endl;
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " input_keyword=\"" << input_keyword << "\"" << endl;
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " option_DEFAULT=" << (option_DEFAULT?"TRUE":"FALSE") << endl;
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme_DEFAULT=\"" << xscheme_DEFAULT << "\"" << endl;
    }
    // start the scan
    //string keyword; //CO20180404 - now a member of the object
    vector<string> vkeyword;
    // tokenize the option
    aurostd::string2tokens(input_keyword,vkeyword,"|"); 
    if(VERBOSE) for(uint i=0;i<vkeyword.size();i++) cerr << "\"" << vkeyword.at(i) << "\"" << endl;
    // loop through the scan
    if(vkeyword.size()>0) {
      // some default
      keyword=vkeyword.at(0);
      for(uint i=0;i<vkeyword.size();i++) 
        if(aurostd::substring2bool(options_FILE,vkeyword.at(i),TRUE))
          keyword=vkeyword.at(i);
      // found one keyword
      if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " keyword=\"" << keyword << "\"" << endl;
      // LOOK FOR EXIST/!EXIST ENTRY
      if(_option_DEFAULT==aurostd_xoptionONOFF) {
        isentry=aurostd::substring2bool(options_FILE,keyword,TRUE);
        if(isentry) {option=TRUE;content_string="ON";}
        if(!isentry) {option=FALSE;content_string="OFF";}
      } // aurostd_xoptionONOFF exit/~exit
      // LOOK FOR ON/OFF MODE WITH strings/schemes.
      if(_option_DEFAULT==0 || _option_DEFAULT==1) {
        if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " LOOK FOR ON/OFF MODE WITH strings/schemes" << endl;
        // start the scan
        isentry=aurostd::substring2bool(options_FILE,keyword,TRUE);
        if(isentry && xscheme_DEFAULT.empty()) {
          content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,1,FALSE));
          if(content_string.empty()){content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,1,TRUE));}  //CO20200731 - "[AFLOW]SYSTEM=" vs. "[AFLOW] SYSTEM = "
          string saus=content_string;content_string="";
          if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " found saus=" << saus << endl;
          vector<string> tokens;aurostd::string2tokens(saus,tokens,",");
          for(uint i=0;i<tokens.size();i++) { //      c<< tokens.at(i) << endl;
            if(tokens.at(i)=="ON" || tokens.at(i)[0]=='T' || tokens.at(i)[0]=='t' || tokens.at(i)[0]=='1' || tokens.at(i)[0]=='Y' || tokens.at(i)[0]=='y') {
              option=TRUE;content_string=saus;} // modify option and value
            if(tokens.at(i)=="OFF"|| tokens.at(i)[0]=='F' || tokens.at(i)[0]=='f' || tokens.at(i)[0]=='0' || tokens.at(i)[0]=='N' || tokens.at(i)[0]=='n') {
              option=FALSE;content_string=saus;} // modify option and value
            // if(tokens.at(i)=="REMOVE_RELAX_1") {option=TRUE;content_string=saus;} // modify option and value // compatibility with SPIN
            // if(tokens.at(i)=="REMOVE_RELAX_2") {option=TRUE;content_string=saus;} // modify option and value // compatibility with SPIN
            if(tokens.at(i)=="REMOVE_RELAX_1") {content_string=saus;} // modify option and value // compatibility with SPIN but dont touch ON/OFF
            if(tokens.at(i)=="REMOVE_RELAX_2") {content_string=saus;} // modify option and value // compatibility with SPIN but dont touch ON/OFF
          }
        }
        // SCHEME MODE
        if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme_DEFAULT=\"" << xscheme_DEFAULT << "\"" << endl;
        if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme_DEFAULT.empty()=" << xscheme_DEFAULT.empty() << endl;
        if(isentry && !xscheme_DEFAULT.empty()) {
          if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " SCHEME MODE" << endl;
          content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,1,FALSE));
          if(content_string.empty()){content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,1,TRUE));}  //CO20200731 - "[AFLOW]SYSTEM=" vs. "[AFLOW] SYSTEM = "
          //ME20181030 - Special case: if the scheme is a Boolean keyword, unset option
          //ME20190107 - Cannot use N or F because it's ambiguous (nitrogen, fluorine)
          string content = aurostd::toupper(content_string);
          if ((content == "OFF") || (content == "FALSE") || (content == "NO")) {
            option = false;
          } else {
            option=isentry;
          }
        }
        if(isentry && (xscheme_DEFAULT.empty() && content_string.empty())) {
          if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " SCHEME MODE EMPTY DEFAULT STILL EMPTY CONTENT" << endl;
          content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,1,FALSE));
          if(content_string.empty()){content_string=aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE,keyword,1,TRUE));}  //CO20200731 - "[AFLOW]SYSTEM=" vs. "[AFLOW] SYSTEM = "
          option=isentry;
        }
        if(!isentry && option_DEFAULT) {
          option=TRUE;content_string="ON";
        }
      } // 0/1 on/off mode
      // LOOK FOR EXIST/!EXIST ENTRY
      if(_option_DEFAULT==aurostd_xoptionMULTI) {
        vector<string> voptions_FILE,vcontent;
        aurostd::string2vectorstring(options_FILE,voptions_FILE);
        isentry=FALSE;content_string="";
        for(uint i=0;i<voptions_FILE.size();i++) {
          if(aurostd::substring2bool(voptions_FILE.at(i),keyword,TRUE)) {
            vector<string> vstrcheck;
            string strcheck=aurostd::toupper(aurostd::RemoveWhiteSpaces(aurostd::substring2string(voptions_FILE.at(i),keyword,1,FALSE)));
            aurostd::StringSubst(strcheck,";",",");
            aurostd::string2tokens(strcheck,vstrcheck,","); 
            for(uint j=0;j<vstrcheck.size();j++) {
              if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " BEFORE keyword=" << keyword << "   " << "vstrcheck.at(j)=" << vstrcheck.at(j) << endl;
              if(aurostd::substring2bool(keyword,"KPOINTS")) {
                if(vstrcheck.at(j)=="A") vstrcheck.at(j)="AUTO";
                if(vstrcheck.at(j)=="G") vstrcheck.at(j)="GAMMA";
                if(vstrcheck.at(j)=="M") vstrcheck.at(j)="MONKHORST_PACK";
              }	 
              if(aurostd::substring2bool(keyword,"IGNORE_AFIX")) {
                ;
              }; // dummy load
              if(aurostd::substring2bool(keyword,"CONVERT_UNIT_CELL")) {
                if(vstrcheck.at(j)=="SPRIM") vstrcheck.at(j)="STANDARD_PRIMITIVE";
                if(vstrcheck.at(j)=="STD_PRIM") vstrcheck.at(j)="STANDARD_PRIMITIVE";
                if(vstrcheck.at(j)=="STANDARD_PRIMITIVE") vstrcheck.at(j)="STANDARD_PRIMITIVE";
                if(vstrcheck.at(j)=="SCONV") vstrcheck.at(j)="STANDARD_CONVENTIONAL";
                if(vstrcheck.at(j)=="STD_CONV") vstrcheck.at(j)="STANDARD_CONVENTIONAL";
                if(vstrcheck.at(j)=="STANDARD_CONVENTIONAL") vstrcheck.at(j)="STANDARD_CONVENTIONAL";
                if(vstrcheck.at(j)=="NIGGLI") vstrcheck.at(j)="NIGGLI";
                if(vstrcheck.at(j)=="MINK") vstrcheck.at(j)="MINKOWSKI";
                if(vstrcheck.at(j)=="MINKOWSKI") vstrcheck.at(j)="MINKOWSKI";
                if(vstrcheck.at(j)=="INCELL") vstrcheck.at(j)="INCELL";
                if(vstrcheck.at(j)=="COMPACT") vstrcheck.at(j)="COMPACT";
                if(vstrcheck.at(j)=="INCOMPACT") vstrcheck.at(j)="COMPACT";
                if(vstrcheck.at(j)=="INWIGNERSEITZ") vstrcheck.at(j)="WIGNERSEITZ";
                if(vstrcheck.at(j)=="WS") vstrcheck.at(j)="WIGNERSEITZ";
                if(vstrcheck.at(j)=="WIGNERSEITZ") vstrcheck.at(j)="WIGNERSEITZ";
                if(vstrcheck.at(j)=="C") vstrcheck.at(j)="CARTESIAN";
                if(vstrcheck.at(j)=="CART") vstrcheck.at(j)="CARTESIAN";
                if(vstrcheck.at(j)=="CARTESIAN") vstrcheck.at(j)="CARTESIAN";
                if(vstrcheck.at(j)=="F") vstrcheck.at(j)="FRACTIONAL";
                if(vstrcheck.at(j)=="FRAC") vstrcheck.at(j)="FRACTIONAL";
                if(vstrcheck.at(j)=="FRACTIONAL") vstrcheck.at(j)="FRACTIONAL";
                if(vstrcheck.at(j)=="D") vstrcheck.at(j)="DIRECT";
                if(vstrcheck.at(j)=="DIR") vstrcheck.at(j)="DIRECT";
                if(vstrcheck.at(j)=="DIRECT") vstrcheck.at(j)="DIRECT";
                if(vstrcheck.at(j)=="PRE") vstrcheck.at(j)="PRESERVE";
                if(vstrcheck.at(j)=="PRES") vstrcheck.at(j)="PRESERVE";
                if(vstrcheck.at(j)=="PRESERVE") vstrcheck.at(j)="PRESERVE";
              }	
              if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " AFTER keyword=" << keyword << "   " << "vstrcheck.at(j)=" << vstrcheck.at(j) << endl;
              vcontent.push_back(vstrcheck.at(j));
            }
          }
        }
        for(uint i=0;i<vcontent.size();i++)
          content_string+=vcontent.at(i)+(i<vcontent.size()-1?",":"");
        aurostd::StringSubst(content_string,"=","_");aurostd::StringSubst(content_string,";",",");
        if(vcontent.size()) isentry=TRUE;
      } // aurostd_xoptionMULTI list
    }
    content_double=aurostd::string2utype<double>(content_string);
    content_int=aurostd::string2utype<int>(content_string);
    content_uint=aurostd::string2utype<uint>(content_string);
    xscheme=content_string;
    aurostd::string2tokens(xscheme,vxscheme,","); 
    if(VERBOSE) if(_option_DEFAULT==aurostd_xoptionMULTI) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - " << __AFLOW_FUNC__ << " vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;

    preserved=FALSE;
    for(uint i=0;i<vxscheme.size()&&!preserved;i++) preserved=(vxscheme.at(i)=="PRESERVED");
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " isentry=" << (isentry?"TRUE":"FALSE") << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_string=\"" << content_string << "\"" << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_double=\"" << content_double << "\"" << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_int=\"" << content_int << "\"" << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_uint=\"" << content_uint << "\"" << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " option=" << (option?"TRUE":"FALSE") << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " preserved=" << (preserved?"TRUE":"FALSE") << endl;
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme=\"" << xscheme << "\"" << endl;
    if(isentry && content_string.empty()) {
      stringstream message;
      message << "Content string empty. content_string=" <<  content_string
        << ", content_double=" <<  content_double
        << ", content_int=" <<  content_int
        << ", content_uint=" << content_uint
        << ", keyword=" << keyword
        << ", isentry=" << isentry;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if(VERBOSE) cerr << "DEBUG - " << __AFLOW_FUNC__ << " END" << endl;
    // return isentry;
  }

  void xoption::scheme2scheme(char c,const string& s) { //CO20210805 - const&
    for(uint i=0;i<vxscheme.size();i++) {
      if(vxscheme.at(i).at(0)==c ||
          vxscheme.at(i).at(0)==aurostd::tolower(c) ||
          vxscheme.at(i).at(0)==aurostd::toupper(c)) {
        xscheme=s;
      }
    }
  }

  void xoption::scheme2scheme(const string& s1,const string& s2) {  //CO20210805 - const&
    for(uint i=0;i<vxscheme.size();i++) {
      if(vxscheme.at(i)==s1 ||
          vxscheme.at(i)==aurostd::tolower(s1) ||
          vxscheme.at(i)==aurostd::toupper(s1)) {
        xscheme=s2;
        //  for(uint i=0;i<vxscheme.size();i++) if(vxscheme.at(i)==s1) scheme=s2;
      }
    }
  }

  bool xoption::isscheme(const string& check) const {                     //CO20180101  //CO20210805 - const&
    // ISSCHEME and FLAG checks only vxscheme... does not manage the ghost, so for example    //SC20200114
    // CONVERT_UNIT_CELL (as flag) will not be confused with CONVERT_UNIT_CELL=STANDARD as method.   //SC20200114
    // Thanks to Marco Esters for getting this bug.   //SC20200114
    string a,b;
    // check schemes list going through vxscheme 1 by 1
    for(uint i=0;i<vxscheme.size();i++) {
      a=aurostd::toupper(vxscheme.at(i));                         // shortcuts
      b=aurostd::toupper(check);                                  // shortcuts
      // cerr << "xoption::isscheme for scheme i=" << i << " " << a << " " << b << endl;  
      aurostd::StringSubst(a,"GAMMA","G");aurostd::StringSubst(b,"GAMMA","G");                      // shortcuts
      aurostd::StringSubst(a,"MONKHORST_PACK","M");aurostd::StringSubst(b,"MONKHORST_PACK","M");    // shortcuts
      aurostd::StringSubst(a,"MP","M");aurostd::StringSubst(b,"MP","M");                            // shortcuts
      aurostd::StringSubst(a,"AUTO","A");aurostd::StringSubst(b,"AUTO","A");                        // shortcuts
      if(a==b) {
        //	cerr << "xoption::isscheme BINGO FOUND SCHEME " << a << " " << b << endl;  
        return TRUE;
      }
    }
    //SC20200310 // THIS IS INCORRECT 
    //SC20200310 // check attached schemes list going through vxsghost 2 by 2  //SC20191227
    //SC20200310 for(uint i=0;i<vxsghost.size();i+=2) {
    //SC20200310   //    cerr << "xoption::isscheme for attached scheme i=" << i << " " << a << " " << b << endl;  
    //SC20200310   a=aurostd::toupper(vxsghost.at(i));                         // shortcuts
    //SC20200310   b=aurostd::toupper(check);                                  // shortcuts
    //SC20200310   if(a==b) {
    //SC20200310     //	cerr << "xoption::isscheme BINGO FOUND ATTACHED SCHEME" << a << " " << b << endl;  
    //SC20200310    return TRUE;
    //SC20200310   }
    //SC20200310 }
    //SC20200310 // nor in scheme nor in attached scheme... exit
    return FALSE;
  }

  bool xoption::refresh() {
    content_string="";
    for(uint i=0;i<vxscheme.size();i++)
      content_string+=vxscheme.at(i)+(i<vxscheme.size()-1?",":"");
    aurostd::StringSubst(content_string,"=","_");aurostd::StringSubst(content_string,";",",");
    content_double=aurostd::string2utype<double>(content_string);
    content_int=aurostd::string2utype<int>(content_string);
    content_uint=aurostd::string2utype<uint>(content_string);
    xscheme=content_string;
    return TRUE;
  }

  // [OBSOLETE] uint xoption::addscheme(string _xscheme)   { return opscheme(_xscheme,TRUE); }
  // [OBSOLETE] uint xoption::purgescheme(string _xscheme) { return opscheme(_xscheme,FALSE); }
  uint xoption::push(const string& _xscheme)        { return opscheme(_xscheme,TRUE); } //CO20210805 - const&
  uint xoption::pop(const string& _xscheme)         { return opscheme(_xscheme,FALSE); }  //CO20210805 - const&

  uint xoption::opscheme(const string& _xscheme,bool operation) { //CO20210805 - const&
    bool VERBOSE=(FALSE || VERBOSE_XOPTION); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    if(operation==TRUE) {
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::opscheme: ADD=" << aurostd::toupper(_xscheme) << endl;
      if(VERBOSE) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: ADD_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
      //CO20181226 START - check that it doesn't already exist, multiples don't affect isscheme, but affects how we iterate through aplopts
      for(uint i=0;i<vxscheme.size();i++){
        if(aurostd::toupper(vxscheme.at(i))==aurostd::toupper(_xscheme)){opscheme(_xscheme,FALSE);} //recursion is GNU's pleasure
      }
      //CO20181226 STOP
      vxscheme.push_back(aurostd::toupper(_xscheme));
      if(VERBOSE) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: ADD_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
    } else {
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::opscheme: PURGE=" << aurostd::toupper(_xscheme) << endl;
      if(VERBOSE) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: PURGE_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
      vector<string> _vxscheme(vxscheme);
      vxscheme.clear();
      for(uint i=0;i<_vxscheme.size();i++) {
        if(aurostd::toupper(_vxscheme.at(i))!=aurostd::toupper(_xscheme)) vxscheme.push_back(_vxscheme.at(i));
      }
      if(VERBOSE) for(uint i=0;i<vxscheme.size();i++) cerr << "DEBUG - aurostd::xoption::opscheme: PURGE_AFTER vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
    }
    refresh();
    return vxscheme.size();
  }

  bool xoption::flag(const string& _xscheme,bool operation) { //CO20210805 - const&
    // [OBSOLETE]    if(operation) push(_xscheme);
    // [OBSOLETE]    if(!operation) pop(_xscheme);
    if(operation) opscheme(_xscheme,TRUE);  // push
    if(!operation) opscheme(_xscheme,FALSE);  // pop
    return operation;
  }

  bool xoption::flag(const string& xscheme) const { return isscheme(xscheme);  } // same as ischeme //CO20210805 - const&

  bool xoption::flag(void) const {  // same as ischeme
    if(vxscheme.size()>0) return TRUE;
    // NO NEED ANYMORE SC20200114    if(vxsghost.size()>0) return TRUE;  //SC20191227
    return FALSE;
  }

  // now for the attached ones.

  bool xoption::isdefined(const string& check) const {                        //SC20200114  //CO20210805 - const&
    // checks only scheme (vxscheme) it does not go through the attached schemes (vxghost).   //SC20200114
    string a,b;   //SC20200114
    // check schemes list going through vxscheme 1 by 1   //SC20200114
    // check attached schemes list going through vxsghost 2 by 2  //SC20191227    //SC20200114
    b=aurostd::toupper(check);                                  // shortcuts   //SC20200114 //CO20220630 - no need to redo over and over again
    for(uint i=0;i<vxsghost.size();i+=2) {   //SC20200114
      //    cerr << "xoption::isscheme for attached scheme i=" << i << " " << a << " " << b << endl;     //SC20200114
      a=aurostd::toupper(vxsghost.at(i));                         // shortcuts   //SC20200114
      if(a==b) {   //SC20200114
        //	cerr << "xoption::isscheme BINGO FOUND ATTACHED SCHEME" << a << " " << b << endl;     //SC20200114
        return TRUE;   //SC20200114
      }   //SC20200114
    }   //SC20200114
    return FALSE;   //SC20200114
  }   //SC20200114

  string xoption::getattachedscheme(const string& xscheme) const {
    bool VERBOSE=(FALSE || VERBOSE_XOPTION); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    if(vxsghost.size()==0) return "";
    for(uint i=0;i<vxsghost.size()-1;i+=2) {
      if(VERBOSE) cerr << i << " --- [" << aurostd::toupper(xscheme) << "] --- [" << aurostd::toupper(vxsghost.at(i)) << "] --- [" << aurostd::toupper(vxsghost.at(i+1)) << "]" << endl;
      if(aurostd::toupper(xscheme)==aurostd::toupper(vxsghost.at(i))) 
        return vxsghost.at(i+1);
    }
    return "";
  }
  template<class utype> utype xoption::getattachedutype(const string& xscheme) const { //CO20200731
    return aurostd::string2utype<utype>(getattachedscheme(xscheme));
  }

  uint xoption::opattachedscheme(const string& _xscheme,const string& attached,bool operation) {  //CO20210805 - const&
    bool VERBOSE=(FALSE || VERBOSE_XOPTION); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    if(operation==TRUE) {
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::opattachedscheme: ADD=" << aurostd::toupper(_xscheme) << endl;
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::opattachedscheme: GHOST=" << attached << endl;
      //CO20181226 START - check that it doesn't already exist, multiples affect getattachedscheme
      for(uint i=0;i<vxsghost.size();i+=2) {
        if(aurostd::toupper(vxsghost.at(i))==aurostd::toupper(_xscheme)) {opattachedscheme(_xscheme,attached,FALSE);} //recursion is GNU's pleasure
      }
      //CO20181226 STOP
      vxsghost.push_back(aurostd::toupper(_xscheme));
      vxsghost.push_back(attached);
    }  else {
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::opattachedscheme: PURGE=" << aurostd::toupper(_xscheme) << endl;
      vector<string> _vxsghost(vxsghost);
      vxsghost.clear();
      for(uint i=0;i<_vxsghost.size();i+=2) {
        if(aurostd::toupper(_vxsghost.at(i))!=aurostd::toupper(_xscheme)) {
          vxsghost.push_back(_vxsghost.at(i));vxsghost.push_back(_vxsghost.at(i+1));
        }
      }
      if(VERBOSE) for(uint i=0;i<vxsghost.size();i++) cerr << "PURGEATTACHED_AFTER vxsghost.at(" << i << ")=" << vxsghost.at(i) << endl;
    } 
    refresh();
    return vxsghost.size();
  }  

  uint xoption::addattachedscheme(const string& _xscheme,const string& attached,bool operation) { //CO20210805 - const&
    if(operation) return opattachedscheme(_xscheme,attached,TRUE);
    return vxsghost.size();
  }

  // [OBSOLETE] uint xoption::purgeattachedscheme(string _xscheme) {
  // [OBSOLETE]   return opattachedscheme(_xscheme,"",FALSE);
  // [OBSOLETE] }

  uint xoption::push_attached(const string& _xscheme,const string& attached) {  //CO20210805 - const&
    return opattachedscheme(_xscheme,attached,TRUE);
  }

  uint xoption::pop_attached(const string& _xscheme) {  //CO20210805 - const&
    return opattachedscheme(_xscheme,"",FALSE);
  }

  bool xoption::args2addattachedscheme(vector<string> &argv,const string _xscheme,const string& _s_search,string string_default) {
    vector<string> cmds;
    return args2addattachedscheme(argv,cmds,_xscheme,_s_search,string_default);
  }

  bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string xscheme,const string& _s_search,string string_default) {
    bool VERBOSE=(FALSE || VERBOSE_XOPTION); //DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    string s_search(_s_search);
    if(aurostd::args2attachedflag(argv,cmds,s_search)) {
      flag(xscheme,TRUE);
      addattachedscheme(xscheme,aurostd::args2attachedstring(argv,s_search,string_default),TRUE);
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " attached=" << aurostd::args2attachedstring(argv,s_search,string_default) << endl;
      return TRUE;
    } 
    aurostd::StringSubst(s_search,"=","");
    if(aurostd::args2flag(argv,cmds,s_search)) {
      //    cerr << aurostd::args2string(argv,s_search,string_default) << endl;
      flag(xscheme,TRUE);
      // [OBSOLETE]      addattachedscheme(xscheme,string_default,TRUE);
      // [OBSOLETE] if(VERBOSE) cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " taking=" << string_default << endl;
      addattachedscheme(xscheme,aurostd::args2string(argv,s_search,string_default),TRUE);
      if(VERBOSE) cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " taking=" << aurostd::args2string(argv,s_search,string_default) << endl;
      return TRUE;
    }
    return FALSE;
  }

  bool xoption::args2addattachedscheme(vector<string> &argv,const string xscheme,const string& _s_search,char const* string_default) {
    return args2addattachedscheme(argv,xscheme,_s_search,string(string_default));
  }

  bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string xscheme,const string& _s_search,char const* string_default) {
    return args2addattachedscheme(argv,cmds,xscheme,_s_search,string(string_default));
  }

  template<class utype> bool xoption::args2addattachedscheme(vector<string> &argv,const string xscheme,const string& _s_search,utype utype_default) {
    return args2addattachedscheme(argv,xscheme,_s_search,aurostd::utype2string(utype_default));
  }

  template<class utype> bool xoption::args2addattachedscheme(vector<string> &argv,vector<string> &cmds,const string xscheme,const string& _s_search,utype utype_default) {
    return args2addattachedscheme(argv,cmds,xscheme,_s_search,aurostd::utype2string(utype_default));
  }

}

#endif  // _AUROSTD_XOPTION_CPP_


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
