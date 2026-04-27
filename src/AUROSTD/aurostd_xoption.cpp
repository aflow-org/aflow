// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 2013-2014
// added template<class utype> bool xoption::args2addattachedscheme SC 2017
// streamline schemes SC 2017

#include "aurostd_xoption.h"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "aurostd.h"
#include "aurostd_argv.h"
#include "aurostd_automatic_template.h"
#include "aurostd_xerror.h"

using std::cerr;
using std::endl;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

#define VERBOSE_XOPTION false // DX20200907

namespace aurostd {
  // ***************************************************************************

  // constructure
  xoption::xoption() {
    free();
  }  // CO20200624 - moved to free()

  // destructor
  xoption::~xoption() {
    free();
  } // CO20200624 - moved to free()

  // free
  void xoption::free() {
    keyword = ""; // DX20180824 - missing from constructor
    isentry = false;
    content_string = "";
    content_double = 0.0;
    content_int = 0;
    content_uint = 0;
    option = false;
    option_default = false;
    xscheme = "";
    vxscheme.clear();
    vxsghost.clear();
    preserved = false;
  }

  // copy fuction
  void xoption::copy(const xoption& b) {
    keyword = b.keyword; // DX20180824 - missing from copy constructor
    isentry = b.isentry;
    content_string = b.content_string;
    content_double = b.content_double;
    content_int = b.content_int;
    content_uint = b.content_uint;
    option = b.option;
    option_default = b.option_default;
    xscheme = b.xscheme;
    vxscheme.clear();
    for (size_t i = 0; i < b.vxscheme.size(); i++) {
      vxscheme.push_back(b.vxscheme[i]);
    }
    vxsghost.clear();
    for (size_t i = 0; i < b.vxsghost.size(); i++) {
      vxsghost.push_back(b.vxsghost[i]);
    }
    preserved = b.preserved;
  }

  // copy conctructor
  xoption::xoption(const xoption& b) {
    //  free();
    // *this=b;
    copy(b);
  }

  // copy operator b=a
  const xoption& xoption::operator=(const xoption& b) {  // operator=
    if (this != &b) {
      free();
      copy(b);
    }
    return *this;
  }

  std::ostream& operator<<(std::ostream& oss, const xoption& a) {
    for (size_t i = 0; i < a.vxscheme.size(); i++) {
      oss << a.vxscheme[i] << (i < a.vxscheme.size() - 1 ? "," : "");
    }
    return oss;
  }

  void xoption::clear() {
    //[CO20200624 - creating objects is SLOW]xoption aflow_option_temp;
    //[CO20200624 - creating objects is SLOW]copy(aflow_option_temp);
    free();
  }

  // **************************************************************************
  // void xoption::options2entry(const string& options_FILE,const string& input_keyword,int _option_DEFAULT,const string& xscheme_DEFAULT) //CO20210805 - const&
  void xoption::options2entry(const string& options_FILE_IN, const string& input_keyword_IN, int option_DEFAULT_IN, const string& xscheme_DEFAULT_IN) {
    const bool VERBOSE = (false || VERBOSE_XOPTION); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;

    // CO20210909 - BIG BUG HERE
    // the following clear() will reset all of the internal xoption variables
    // if we pass one of these variables into options2entry(), it is reset as well
    // see for example, this construction:
    // vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.options2entry(AflowIn,_STROPT_+"NELM=",false,vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.xscheme);
    // xscheme gets cleared before it's set
    // in order to preserve this construction and prevent headaches, make copies of the inputs BEFORE the clear
    const string options_FILE = options_FILE_IN;
    const string input_keyword = input_keyword_IN;
    const int _option_DEFAULT = option_DEFAULT_IN;
    const string xscheme_DEFAULT = xscheme_DEFAULT_IN;

    clear();  // CO20210909 - DANGEROUS! see note above

    bool option_DEFAULT = false;
    if (_option_DEFAULT == 0) {
      option_DEFAULT = false; // it is a int.. it might be -1
    }
    if (_option_DEFAULT == 1) {
      option_DEFAULT = true; // it is a int.. it might be -1
    }
    isentry = option_DEFAULT;
    option = option_DEFAULT;
    content_string = xscheme_DEFAULT;
    xscheme = xscheme_DEFAULT;
    preserved = false;   // DEFAULT
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " BEGIN " << endl;
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " input_keyword=\"" << input_keyword << "\"" << endl;
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " option_DEFAULT=" << (option_DEFAULT ? "true" : "false") << endl;
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme_DEFAULT=\"" << xscheme_DEFAULT << "\"" << endl;
    }
    // start the scan
    // string keyword; //CO20180404 - now a member of the object
    vector<string> vkeyword;
    // tokenize the option
    aurostd::string2tokens(input_keyword, vkeyword, "|");
    if (VERBOSE) {
      for (size_t i = 0; i < vkeyword.size(); i++) {
        cerr << "\"" << vkeyword[i] << "\"" << endl;
      }
    }
    // loop through the scan
    if (!vkeyword.empty()) {
      // some default
      keyword = vkeyword[0];
      for (size_t i = 0; i < vkeyword.size(); i++) {
        if (aurostd::substring2bool(options_FILE, vkeyword[i], true)) {
          keyword = vkeyword[i];
        }
      }
      // found one keyword
      if (VERBOSE) {
        cerr << "DEBUG - " << __AFLOW_FUNC__ << " keyword=\"" << keyword << "\"" << endl;
      }
      // LOOK FOR EXIST/!EXIST ENTRY
      if (_option_DEFAULT == aurostd_xoptionONOFF) {
        isentry = aurostd::substring2bool(options_FILE, keyword, true);
        if (isentry) {
          option = true;
          content_string = "ON";
        }
        if (!isentry) {
          option = false;
          content_string = "OFF";
        }
      } // aurostd_xoptionONOFF exit/~exit
      // LOOK FOR ON/OFF MODE WITH strings/schemes.
      if (_option_DEFAULT == 0 || _option_DEFAULT == 1) {
        if (VERBOSE) {
          cerr << "DEBUG - " << __AFLOW_FUNC__ << " LOOK FOR ON/OFF MODE WITH strings/schemes" << endl;
        }
        // start the scan
        isentry = aurostd::substring2bool(options_FILE, keyword, true);
        if (isentry && xscheme_DEFAULT.empty()) {
          content_string = aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE, keyword, 1, false));
          if (content_string.empty()) {
            content_string = aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE, keyword, 1, true));
          }  // CO20200731 - "[AFLOW]SYSTEM=" vs. "[AFLOW] SYSTEM = "
          const string saus = content_string;
          content_string = "";
          if (VERBOSE) {
            cerr << "DEBUG - " << __AFLOW_FUNC__ << " found saus=" << saus << endl;
          }
          vector<string> tokens;
          aurostd::string2tokens(saus, tokens, ",");
          for (size_t i = 0; i < tokens.size(); i++) { //      c<< tokens[i] << endl;
            if (tokens[i] == "ON" || tokens[i][0] == 'T' || tokens[i][0] == 't' || tokens[i][0] == '1' || tokens[i][0] == 'Y' || tokens[i][0] == 'y') {
              option = true;
              content_string = saus;
            } // modify option and value
            if (tokens[i] == "OFF" || tokens[i][0] == 'F' || tokens[i][0] == 'f' || tokens[i][0] == '0' || tokens[i][0] == 'N' || tokens[i][0] == 'n') {
              option = false;
              content_string = saus;
            } // modify option and value
            // if(tokens[i]=="REMOVE_RELAX_1") {option=true;content_string=saus;} // modify option and value // compatibility with SPIN
            // if(tokens[i]=="REMOVE_RELAX_2") {option=true;content_string=saus;} // modify option and value // compatibility with SPIN
            if (tokens[i] == "REMOVE_RELAX_1") {
              content_string = saus;
            } // modify option and value // compatibility with SPIN but dont touch ON/OFF
            if (tokens[i] == "REMOVE_RELAX_2") {
              content_string = saus;
            } // modify option and value // compatibility with SPIN but dont touch ON/OFF
          }
        }
        // SCHEME MODE
        if (VERBOSE) {
          cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme_DEFAULT=\"" << xscheme_DEFAULT << "\"" << endl;
        }
        if (VERBOSE) {
          cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme_DEFAULT.empty()=" << xscheme_DEFAULT.empty() << endl;
        }
        if (isentry && !xscheme_DEFAULT.empty()) {
          if (VERBOSE) {
            cerr << "DEBUG - " << __AFLOW_FUNC__ << " SCHEME MODE" << endl;
          }
          content_string = aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE, keyword, 1, false));
          if (content_string.empty()) {
            content_string = aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE, keyword, 1, true));
          }  // CO20200731 - "[AFLOW]SYSTEM=" vs. "[AFLOW] SYSTEM = "
          // ME20181030 - Special case: if the scheme is a Boolean keyword, unset option
          // ME20190107 - Cannot use N or F because it's ambiguous (nitrogen, fluorine)
          const string content = aurostd::toupper(content_string);
          if ((content == "OFF") || (content == "false") || (content == "NO")) {
            option = false;
          } else {
            option = isentry;
          }
        }
        if (isentry && (xscheme_DEFAULT.empty() && content_string.empty())) {
          if (VERBOSE) {
            cerr << "DEBUG - " << __AFLOW_FUNC__ << " SCHEME MODE EMPTY DEFAULT STILL EMPTY CONTENT" << endl;
          }
          content_string = aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE, keyword, 1, false));
          if (content_string.empty()) {
            content_string = aurostd::RemoveWhiteSpaces(aurostd::substring2string(options_FILE, keyword, 1, true));
          }  // CO20200731 - "[AFLOW]SYSTEM=" vs. "[AFLOW] SYSTEM = "
          option = isentry;
        }
        if (!isentry && option_DEFAULT) {
          option = true;
          content_string = "ON";
        }
      } // 0/1 on/off mode
      // LOOK FOR EXIST/!EXIST ENTRY
      if (_option_DEFAULT == aurostd_xoptionMULTI) {
        vector<string> voptions_FILE;
        vector<string> vcontent;
        aurostd::string2vectorstring(options_FILE, voptions_FILE);
        isentry = false;
        content_string = "";
        for (size_t i = 0; i < voptions_FILE.size(); i++) {
          if (aurostd::substring2bool(voptions_FILE[i], keyword, true)) {
            vector<string> vstrcheck;
            string strcheck = aurostd::toupper(aurostd::RemoveWhiteSpaces(aurostd::substring2string(voptions_FILE[i], keyword, 1, false)));
            aurostd::StringSubstInPlace(strcheck, ";", ",");
            aurostd::string2tokens(strcheck, vstrcheck, ",");
            for (size_t j = 0; j < vstrcheck.size(); j++) {
              if (VERBOSE) {
                cerr << "DEBUG - " << __AFLOW_FUNC__ << " BEFORE keyword=" << keyword << "   " << "vstrcheck[j]=" << vstrcheck[j] << endl;
              }
              if (aurostd::substring2bool(keyword, "KPOINTS")) {
                if (vstrcheck[j] == "A") {
                  vstrcheck[j] = "AUTO";
                }
                if (vstrcheck[j] == "G") {
                  vstrcheck[j] = "GAMMA";
                }
                if (vstrcheck[j] == "M") {
                  vstrcheck[j] = "MONKHORST_PACK";
                }
              }
              if (aurostd::substring2bool(keyword, "IGNORE_AFIX")) {
                ;
              }; // dummy load
              if (aurostd::substring2bool(keyword, "CONVERT_UNIT_CELL")) {
                if (vstrcheck[j] == "SPRIM") {
                  vstrcheck[j] = "STANDARD_PRIMITIVE";
                }
                if (vstrcheck[j] == "STD_PRIM") {
                  vstrcheck[j] = "STANDARD_PRIMITIVE";
                }
                if (vstrcheck[j] == "STANDARD_PRIMITIVE") {
                  vstrcheck[j] = "STANDARD_PRIMITIVE";
                }
                if (vstrcheck[j] == "SCONV") {
                  vstrcheck[j] = "STANDARD_CONVENTIONAL";
                }
                if (vstrcheck[j] == "STD_CONV") {
                  vstrcheck[j] = "STANDARD_CONVENTIONAL";
                }
                if (vstrcheck[j] == "STANDARD_CONVENTIONAL") {
                  vstrcheck[j] = "STANDARD_CONVENTIONAL";
                }
                if (vstrcheck[j] == "NIGGLI") {
                  vstrcheck[j] = "NIGGLI";
                }
                if (vstrcheck[j] == "MINK") {
                  vstrcheck[j] = "MINKOWSKI";
                }
                if (vstrcheck[j] == "MINKOWSKI") {
                  vstrcheck[j] = "MINKOWSKI";
                }
                if (vstrcheck[j] == "INCELL") {
                  vstrcheck[j] = "INCELL";
                }
                if (vstrcheck[j] == "COMPACT") {
                  vstrcheck[j] = "COMPACT";
                }
                if (vstrcheck[j] == "INCOMPACT") {
                  vstrcheck[j] = "COMPACT";
                }
                if (vstrcheck[j] == "INWIGNERSEITZ") {
                  vstrcheck[j] = "WIGNERSEITZ";
                }
                if (vstrcheck[j] == "WS") {
                  vstrcheck[j] = "WIGNERSEITZ";
                }
                if (vstrcheck[j] == "WIGNERSEITZ") {
                  vstrcheck[j] = "WIGNERSEITZ";
                }
                if (vstrcheck[j] == "C") {
                  vstrcheck[j] = "CARTESIAN";
                }
                if (vstrcheck[j] == "CART") {
                  vstrcheck[j] = "CARTESIAN";
                }
                if (vstrcheck[j] == "CARTESIAN") {
                  vstrcheck[j] = "CARTESIAN";
                }
                if (vstrcheck[j] == "F") {
                  vstrcheck[j] = "FRACTIONAL";
                }
                if (vstrcheck[j] == "FRAC") {
                  vstrcheck[j] = "FRACTIONAL";
                }
                if (vstrcheck[j] == "FRACTIONAL") {
                  vstrcheck[j] = "FRACTIONAL";
                }
                if (vstrcheck[j] == "D") {
                  vstrcheck[j] = "DIRECT";
                }
                if (vstrcheck[j] == "DIR") {
                  vstrcheck[j] = "DIRECT";
                }
                if (vstrcheck[j] == "DIRECT") {
                  vstrcheck[j] = "DIRECT";
                }
                if (vstrcheck[j] == "PRE") {
                  vstrcheck[j] = "PRESERVE";
                }
                if (vstrcheck[j] == "PRES") {
                  vstrcheck[j] = "PRESERVE";
                }
                if (vstrcheck[j] == "PRESERVE") {
                  vstrcheck[j] = "PRESERVE";
                }
              }
              if (VERBOSE) {
                cerr << "DEBUG - " << __AFLOW_FUNC__ << " AFTER keyword=" << keyword << "   " << "vstrcheck[j]=" << vstrcheck[j] << endl;
              }
              vcontent.push_back(vstrcheck[j]);
            }
          }
        }
        for (size_t i = 0; i < vcontent.size(); i++) {
          content_string += vcontent[i] + (i < vcontent.size() - 1 ? "," : "");
        }
        aurostd::StringSubstInPlace(content_string, "=", "_");
        aurostd::StringSubstInPlace(content_string, ";", ",");
        if (!vcontent.empty()) {
          isentry = true;
        }
      } // aurostd_xoptionMULTI list
    }
    content_double = aurostd::string2utype<double>(content_string);
    content_int = aurostd::string2utype<int>(content_string);
    content_uint = aurostd::string2utype<uint>(content_string);
    xscheme = content_string;
    aurostd::string2tokens(xscheme, vxscheme, ",");
    if (VERBOSE) {
      if (_option_DEFAULT == aurostd_xoptionMULTI) {
        for (size_t i = 0; i < vxscheme.size(); i++) {
          cerr << "DEBUG - " << __AFLOW_FUNC__ << " vxscheme.at(" << i << ")=" << vxscheme[i] << endl;
        }
      }
    }

    preserved = false;
    for (size_t i = 0; i < vxscheme.size() && !preserved; i++) {
      preserved = (vxscheme[i] == "PRESERVED");
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " isentry=" << (isentry ? "true" : "false") << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_string=\"" << content_string << "\"" << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_double=\"" << content_double << "\"" << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_int=\"" << content_int << "\"" << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " content_uint=\"" << content_uint << "\"" << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " option=" << (option ? "true" : "false") << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " preserved=" << (preserved ? "true" : "false") << endl;
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " xscheme=\"" << xscheme << "\"" << endl;
    }
    if (isentry && content_string.empty()) {
      stringstream message;
      message << "Content string empty. content_string=" << content_string << ", content_double=" << content_double << ", content_int=" << content_int << ", content_uint=" << content_uint
              << ", keyword=" << keyword << ", isentry=" << isentry;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if (VERBOSE) {
      cerr << "DEBUG - " << __AFLOW_FUNC__ << " END" << endl;
    }
    // return isentry;
  }

  void xoption::scheme2scheme(char c, const string& s) { // CO20210805 - const&
    for (size_t i = 0; i < vxscheme.size(); i++) {
      if (vxscheme[i].at(0) == c || vxscheme[i].at(0) == aurostd::tolower(c) || vxscheme[i].at(0) == aurostd::toupper(c)) {
        xscheme = s;
      }
    }
  }

  void xoption::scheme2scheme(const string& s1, const string& s2) {  // CO20210805 - const&
    for (size_t i = 0; i < vxscheme.size(); i++) {
      if (vxscheme[i] == s1 || vxscheme[i] == aurostd::tolower(s1) || vxscheme[i] == aurostd::toupper(s1)) {
        xscheme = s2;
        //  for(size_t i=0;i<vxscheme.size();i++) if(vxscheme[i]==s1) scheme=s2;
      }
    }
  }

  bool xoption::isscheme(const string& check) const {                     // CO20180101  //CO20210805 - const&
    // ISSCHEME and FLAG checks only vxscheme... does not manage the ghost, so for example    //SC20200114
    // CONVERT_UNIT_CELL (as flag) will not be confused with CONVERT_UNIT_CELL=STANDARD as method.   //SC20200114
    // Thanks to Marco Esters for getting this bug.   //SC20200114
    string a;
    string b;
    // check schemes list going through vxscheme 1 by 1
    for (size_t i = 0; i < vxscheme.size(); i++) {
      a = aurostd::toupper(vxscheme[i]);                         // shortcuts
      b = aurostd::toupper(check);                                  // shortcuts
      // cerr << "xoption::isscheme for scheme i=" << i << " " << a << " " << b << endl;
      aurostd::StringSubstInPlace(a, "GAMMA", "G");
      aurostd::StringSubstInPlace(b, "GAMMA", "G");                      // shortcuts
      aurostd::StringSubstInPlace(a, "MONKHORST_PACK", "M");
      aurostd::StringSubstInPlace(b, "MONKHORST_PACK", "M");    // shortcuts
      aurostd::StringSubstInPlace(a, "MP", "M");
      aurostd::StringSubstInPlace(b, "MP", "M");                            // shortcuts
      aurostd::StringSubstInPlace(a, "AUTO", "A");
      aurostd::StringSubstInPlace(b, "AUTO", "A");                        // shortcuts
      if (a == b) {
        //	cerr << "xoption::isscheme BINGO FOUND SCHEME " << a << " " << b << endl;
        return true;
      }
    }
    // SC20200310 // THIS IS INCORRECT
    // SC20200310 // check attached schemes list going through vxsghost 2 by 2  //SC20191227
    // SC20200310 for(size_t i=0;i<vxsghost.size();i+=2) {
    // SC20200310   //    cerr << "xoption::isscheme for attached scheme i=" << i << " " << a << " " << b << endl;
    // SC20200310   a=aurostd::toupper(vxsghost.at(i));                         // shortcuts
    // SC20200310   b=aurostd::toupper(check);                                  // shortcuts
    // SC20200310   if(a==b) {
    // SC20200310     //	cerr << "xoption::isscheme BINGO FOUND ATTACHED SCHEME" << a << " " << b << endl;
    // SC20200310    return true;
    // SC20200310   }
    // SC20200310 }
    // SC20200310 // nor in scheme nor in attached scheme... exit
    return false;
  }

  bool xoption::refresh() {
    content_string = "";
    for (size_t i = 0; i < vxscheme.size(); i++) {
      content_string += vxscheme[i] + (i < vxscheme.size() - 1 ? "," : "");
    }
    aurostd::StringSubstInPlace(content_string, "=", "_");
    aurostd::StringSubstInPlace(content_string, ";", ",");
    content_double = aurostd::string2utype<double>(content_string);
    content_int = aurostd::string2utype<int>(content_string);
    content_uint = aurostd::string2utype<uint>(content_string);
    xscheme = content_string;
    return true;
  }

  uint xoption::push(const string& _xscheme) {
    return opscheme(_xscheme, true);
  } // CO20210805 - const&
  uint xoption::pop(const string& _xscheme) {
    return opscheme(_xscheme, false);
  }  // CO20210805 - const&

  uint xoption::opscheme(const string& _xscheme, bool operation) { // CO20210805 - const&
    const bool VERBOSE = (false || VERBOSE_XOPTION); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    if (operation == true) {
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::opscheme: ADD=" << aurostd::toupper(_xscheme) << endl;
      }
      if (VERBOSE) {
        for (size_t i = 0; i < vxscheme.size(); i++) {
          cerr << "DEBUG - aurostd::xoption::opscheme: ADD_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
        }
      }
      // CO20181226 START - check that it doesn't already exist, multiples don't affect isscheme, but affects how we iterate through aplopts
      for (size_t i = 0; i < vxscheme.size(); i++) {
        if (aurostd::toupper(vxscheme[i]) == aurostd::toupper(_xscheme)) {
          opscheme(_xscheme, false);
        } // recursion is GNU's pleasure
      }
      // CO20181226 STOP
      vxscheme.push_back(aurostd::toupper(_xscheme));
      if (VERBOSE) {
        for (size_t i = 0; i < vxscheme.size(); i++) {
          cerr << "DEBUG - aurostd::xoption::opscheme: ADD_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
        }
      }
    } else {
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::opscheme: PURGE=" << aurostd::toupper(_xscheme) << endl;
      }
      if (VERBOSE) {
        for (size_t i = 0; i < vxscheme.size(); i++) {
          cerr << "DEBUG - aurostd::xoption::opscheme: PURGE_BEFORE vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
        }
      }
      vector<string> _vxscheme(vxscheme);
      vxscheme.clear();
      for (size_t i = 0; i < _vxscheme.size(); i++) {
        if (aurostd::toupper(_vxscheme[i]) != aurostd::toupper(_xscheme)) {
          vxscheme.push_back(_vxscheme[i]);
        }
      }
      if (VERBOSE) {
        for (size_t i = 0; i < vxscheme.size(); i++) {
          cerr << "DEBUG - aurostd::xoption::opscheme: PURGE_AFTER vxscheme.at(" << i << ")=" << vxscheme.at(i) << endl;
        }
      }
    }
    refresh();
    return vxscheme.size();
  }

  bool xoption::flag(const string& _xscheme, bool operation) { // CO20210805 - const&
    if (operation) {
      opscheme(_xscheme, true);  // push
    }
    if (!operation) {
      opscheme(_xscheme, false);  // pop
    }
    return operation;
  }

  bool xoption::flag(const string& xscheme) const {
    return isscheme(xscheme);
  } // same as ischeme //CO20210805 - const&

  bool xoption::flag() const {  // same as ischeme
    if (!vxscheme.empty()) {
      return true;
    }
    // NO NEED ANYMORE SC20200114    if(vxsghost.size()>0) return true;  //SC20191227
    return false;
  }

  // now for the attached ones.

  bool xoption::isdefined(const string& check) const {                        // SC20200114  //CO20210805 - const&
    // checks only scheme (vxscheme) it does not go through the attached schemes (vxghost).   //SC20200114
    string a;
    string b;   // SC20200114
    // check schemes list going through vxscheme 1 by 1   //SC20200114
    // check attached schemes list going through vxsghost 2 by 2  //SC20191227    //SC20200114
    b = aurostd::toupper(check);                                  // shortcuts   //SC20200114 //CO20220630 - no need to redo over and over again
    for (size_t i = 0; i < vxsghost.size(); i += 2) {   // SC20200114
      //    cerr << "xoption::isscheme for attached scheme i=" << i << " " << a << " " << b << endl;     //SC20200114
      a = aurostd::toupper(vxsghost[i]);                         // shortcuts   //SC20200114
      if (a == b) {   // SC20200114
        //	cerr << "xoption::isscheme BINGO FOUND ATTACHED SCHEME" << a << " " << b << endl;     //SC20200114
        return true;   // SC20200114
      }   // SC20200114
    }   // SC20200114
    return false;   // SC20200114
  }   // SC20200114

  string xoption::getattachedscheme(const string& xscheme) const {
    const bool VERBOSE = (false || VERBOSE_XOPTION); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    if (vxsghost.empty()) {
      return "";
    }
    for (size_t i = 0; i < vxsghost.size() - 1; i += 2) {
      if (VERBOSE) {
        cerr << i << " --- [" << aurostd::toupper(xscheme) << "] --- [" << aurostd::toupper(vxsghost.at(i)) << "] --- [" << aurostd::toupper(vxsghost.at(i + 1)) << "]" << endl;
      }
      if (aurostd::toupper(xscheme) == aurostd::toupper(vxsghost[i])) {
        return vxsghost.at(i + 1);
      }
    }
    return "";
  }
  template <class utype> utype xoption::getattachedutype(const string& xscheme) const { // CO20200731
    return aurostd::string2utype<utype>(getattachedscheme(xscheme));
  }
#define AST_TEMPLATE(utype) template utype xoption::getattachedutype(const string&) const;
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  uint xoption::opattachedscheme(const string& _xscheme, const string& attached, bool operation) {  // CO20210805 - const&
    const bool VERBOSE = (false || VERBOSE_XOPTION); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    if (operation == true) {
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::opattachedscheme: ADD=" << aurostd::toupper(_xscheme) << endl;
      }
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::opattachedscheme: GHOST=" << attached << endl;
      }
      // CO20181226 START - check that it doesn't already exist, multiples affect getattachedscheme
      for (size_t i = 0; i < vxsghost.size(); i += 2) {
        if (aurostd::toupper(vxsghost[i]) == aurostd::toupper(_xscheme)) {
          opattachedscheme(_xscheme, attached, false);
        } // recursion is GNU's pleasure
      }
      // CO20181226 STOP
      vxsghost.push_back(aurostd::toupper(_xscheme));
      vxsghost.push_back(attached);
    } else {
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::opattachedscheme: PURGE=" << aurostd::toupper(_xscheme) << endl;
      }
      vector<string> _vxsghost(vxsghost);
      vxsghost.clear();
      for (size_t i = 0; i < _vxsghost.size(); i += 2) {
        if (aurostd::toupper(_vxsghost[i]) != aurostd::toupper(_xscheme)) {
          vxsghost.push_back(_vxsghost[i]);
          vxsghost.push_back(_vxsghost.at(i + 1));
        }
      }
      if (VERBOSE) {
        for (size_t i = 0; i < vxsghost.size(); i++) {
          cerr << "PURGEATTACHED_AFTER vxsghost.at(" << i << ")=" << vxsghost.at(i) << endl;
        }
      }
    }
    refresh();
    return vxsghost.size();
  }

  uint xoption::addattachedscheme(const string& _xscheme, const string& attached, bool operation) { // CO20210805 - const&
    if (operation) {
      return opattachedscheme(_xscheme, attached, true);
    }
    return vxsghost.size();
  }

  uint xoption::push_attached(const string& _xscheme, const string& attached) {  // CO20210805 - const&
    return opattachedscheme(_xscheme, attached, true);
  }

  uint xoption::pop_attached(const string& _xscheme) {  // CO20210805 - const&
    return opattachedscheme(_xscheme, "", false);
  }

  bool xoption::args2addattachedscheme(vector<string>& argv, const string _xscheme, const string& _s_search, string string_default) {
    vector<string> cmds;
    return args2addattachedscheme(argv, cmds, _xscheme, _s_search, string_default);
  }

  bool xoption::args2addattachedscheme(vector<string>& argv, vector<string>& cmds, const string xscheme, const string& _s_search, string string_default) {
    const bool VERBOSE = (false || VERBOSE_XOPTION); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG;
    string s_search(_s_search);
    if (aurostd::args2attachedflag(argv, cmds, s_search)) {
      flag(xscheme, true);
      addattachedscheme(xscheme, aurostd::args2attachedstring(argv, s_search, string_default), true);
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " attached=" << aurostd::args2attachedstring(argv, s_search, string_default) << endl;
      }
      return true;
    }
    aurostd::StringSubstInPlace(s_search, "=", "");
    if (aurostd::args2flag(argv, cmds, s_search)) {
      //    cerr << aurostd::args2string(argv,s_search,string_default) << endl;
      flag(xscheme, true);
      addattachedscheme(xscheme, aurostd::args2string(argv, s_search, string_default), true);
      if (VERBOSE) {
        cerr << "DEBUG - aurostd::xoption::args2addscheme: xscheme=" << xscheme << " s_search=" << s_search << " taking=" << aurostd::args2string(argv, s_search, string_default) << endl;
      }
      return true;
    }
    return false;
  }

  bool xoption::args2addattachedscheme(vector<string>& argv, const string xscheme, const string& _s_search, const char* string_default) {
    return args2addattachedscheme(argv, xscheme, _s_search, string(string_default));
  }

  bool xoption::args2addattachedscheme(vector<string>& argv, vector<string>& cmds, const string xscheme, const string& _s_search, const char* string_default) {
    return args2addattachedscheme(argv, cmds, xscheme, _s_search, string(string_default));
  }

  template <class utype> bool xoption::args2addattachedscheme(vector<string>& argv, const string xscheme, const string& _s_search, utype utype_default) {
    return args2addattachedscheme(argv, xscheme, _s_search, aurostd::utype2string(utype_default));
  }

  template <class utype> bool xoption::args2addattachedscheme(vector<string>& argv, vector<string>& cmds, const string xscheme, const string& _s_search, utype utype_default) {
    return args2addattachedscheme(argv, cmds, xscheme, _s_search, aurostd::utype2string(utype_default));
  }

} // namespace aurostd

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
