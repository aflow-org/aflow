// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AUROSTD_ARGV_CPP_
#define _AUROSTD_ARGV_CPP_

#include "aurostd_argv.h"

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "aurostd.h"
#include "aurostd_automatic_template.h"
#include "aurostd_defs.h"
#include "aurostd_xerror.h"
#include "aurostd_xscalar.h"
#include "aurostd_xvector.h"

using std::cerr;
using std::endl;
using std::string;
using std::vector;

#define VERBOSE_ARGV false // DX20200907

namespace aurostd {  // namespace aurostd
  string attach(const string& s1) {
    return s1;
  }
  string attach(const string& s1, const string& s2) {
    return s1 + "|" + s2;
  }
  string attach(const string& s1, const string& s2, const string& s3) {
    return s1 + "|" + s2 + "|" + s3;
  }
  string attach(const string& s1, const string& s2, const string& s3, const string& s4) {
    return s1 + "|" + s2 + "|" + s3 + "|" + s4;
  }
  string attach(const string& s1, const string& s2, const string& s3, const string& s4, const string& s5) {
    return s1 + "|" + s2 + "|" + s3 + "|" + s4 + "|" + s5;
  }
} // namespace aurostd

// ***************************************************************************
// Function get_arguments_from_input
// ***************************************************************************
// this function creates the argv vector as vector<string> easier to handle than the char **argv
namespace aurostd {  // namespace aurostd
  vector<string> get_arguments_from_input(int _argc, char** _argv) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    vector<string> out_argv;
    for (int i = 0; i < _argc; i++) {
      string argi(_argv[i]);
      if (argi == "-np" || argi == "-npmax") {
        argi = string("-") + argi;
      }
      if (argi == "-f" || argi == "-F") {
        argi = string("-") + argi;
      }
      if (argi == "-d" || argi == "-D") {
        argi = string("-") + argi;
      }
      if (argi.size() >= 2) {
        if (argi[0] == '-' && argi[1] == 'n') {
          argi = string("-") + argi;  // for -np
        }
      }
      //   if(argi.size()>=2) if(argi[0]=='-' && argi[1]=='m') argi=string("-")+argi;  // for -machine
      if (argi.size() >= 2) {
        if (argi[0] == '-' && argi[1] == 'f') {
          argi = string("-") + argi;  // for -f
        }
      }
      if (argi.size() >= 2) {
        if (argi[0] == '-' && argi[1] == 'F') {
          argi = string("-") + argi;  // for -F
        }
      }
      if (argi.size() >= 2) {
        if (argi[0] == '-' && argi[1] == 'd') {
          argi = string("-") + argi;  // for -d
        }
      }
      if (argi.size() >= 2) {
        if (argi[0] == '-' && argi[1] == 'D') {
          argi = string("-") + argi;  // for -D
        }
      }

      if (argi == "--machine") {
        argi += string("=");                  // forcing "=" after machine !
      }
      if (argi == "--machine_name") {
        argi += string("=");             // forcing "=" after machine_name ! //HE20220309
      }
      if (argi == "--aflowlib") {
        argi += string("=");                 // forcing "=" after aflowlib !
      }
      if (argi == "--np") {
        argi += string("=");                       // forcing "=" after np !
      }
      if (argi.at(argi.size() - 1) == '=' && i < _argc - 1) {
        argi += string(_argv[i + 1]);
        i++;
      }  // fixing space after "= "
      out_argv.push_back(argi);
    }
    if (VERBOSE) {
      for (size_t i = 0; i < out_argv.size(); i++) {
        cerr << "out_argv.at(" << i << ")=" << out_argv[i] << endl;
      }
    }
    return out_argv;
  }
} // namespace aurostd

// ***************************************************************************
// Function args2flag without/with  commands
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  bool args2flag(const vector<string>& argv, const string& s0) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t j = 0; j < tokens.size(); j++) {
      for (size_t i = 0; i < argv.size(); i++) {
        if (argv[i] == tokens[j]) {
          return true;
        }
      }
    }
    return false;
  }

  bool args2flag(const vector<string>& argv, std::vector<string>& cmds, const string& s0) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t j = 0; j < tokens.size(); j++) {
      cmds.push_back(tokens[j]);
    }
    for (size_t i = 0; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          return true;
        }
      }
    }
    return false;
  }
} // namespace aurostd

// ***************************************************************************
// args2utype of utype type
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template <class utype> utype args2utype(const vector<string>& argv, const string& s0, utype def_out) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    utype out = def_out;
    for (size_t i = 1; i < argv.size() - 1; i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          if (_isfloat(out)) {
            out = (utype) atof(argv.at(i + 1).c_str());
          } else {
            out = (utype) atoi(argv.at(i + 1).c_str());
          }
        } // OLD
      }
    }
    return out;
  }
#define AST_TEMPLATE(atype) template atype args2utype(const vector<string>&, const string&, atype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
} // namespace aurostd

// ***************************************************************************
// Function get xvector from input
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  template <class utype> xvector<utype> args2xvectorutype(const vector<string>& argv, const string& s0, const xvector<utype>& def_out) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    xvector<utype> out(def_out.lrows, def_out.urows);
    out = def_out;
    for (size_t i = 0; i < argv.size(); i++) {
      if (i + out.rows < argv.size()) {
        for (size_t j = 0; j < tokens.size(); j++) {
          if (argv[i] == tokens[j]) {
            for (int k = 0; k < out.rows; k++) {
              if (_isfloat(out(1))) {
                out(k + out.lrows) = (utype) atof(argv.at(i + k + 1).c_str());
              } else {
                out(k + out.lrows) = (utype) atoi(argv.at(i + k + 1).c_str());
              }
            }
          }
        }
      }
    }
    return out;
  }
#define AST_TEMPLATE(atype) template xvector<atype> args2xvectorutype(const vector<string>&, const string&, const xvector<atype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> xvector<utype> args2xvectorutype(const vector<string>& argv, const string& s0, int dim) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    const xvector<utype> out(1, dim);
    for (size_t i = 0; i < argv.size(); i++) {
      if (i + out.rows < argv.size()) {
        for (size_t j = 0; j < tokens.size(); j++) {
          if (argv[i] == tokens[j]) {
            for (int k = 0; k < out.rows; k++) {
              if (_isfloat(out(1))) {
                out(k + out.lrows) = (utype) atof(argv.at(i + k + 1).c_str());
              } else {
                out(k + out.lrows) = (utype) atoi(argv.at(i + k + 1).c_str());
              }
            }
          }
        }
      }
    }
    return out; // something phony to keep t used
  }
#define AST_TEMPLATE(atype) template xvector<atype> args2xvectorutype(const vector<string>&, const string&, int);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
} // namespace aurostd
// ***************************************************************************
// Function get vector/deque from input
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template <class utype> vector<utype> args2vectorutype(const vector<string>& argv, const string& s0) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    vector<utype> out;
    for (size_t i = 0; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          for (size_t k = 0; k < argv.size(); k++) {
            if (_isfloat(out.at(0))) {
              out.push_back((utype) atof(argv.at(i + k + 1).c_str()));
            } else {
              out.push_back((utype) atoi(argv.at(i + k + 1).c_str()));
            }
          }
        }
      }
    }
    return out;
  }

  template <class utype> deque<utype> args2dequeutype(const deque<string>& argv, const string& s0) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    deque<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    deque<utype> out;
    for (size_t i = 0; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          for (size_t k = 0; k < argv.size(); k++) {
            if (_isfloat(out.at(0))) {
              out.push_back((utype) atof(argv.at(i + k + 1).c_str()));
            } else {
              out.push_back((utype) atoi(argv.at(i + k + 1).c_str()));
            }
          }
        }
      }
    }
    return out;
  }
} // namespace aurostd

// ***************************************************************************
// Functions args2string functions
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  string args2string(const vector<string>& argv, const string& s0, const string& s_def) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t i = 0; i < argv.size() - 1; i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          return argv.at(i + 1);
        }
      }
    }
    return s_def;
  }

  string args2string(const vector<string>& argv, vector<string>& cmds, const string& s0, const string& s_def) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t j = 0; j < tokens.size(); j++) {
      cmds.push_back(tokens[j]);
    }
    for (size_t i = 0; i < argv.size() - 1; i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          return argv.at(i + 1);
        }
      }
    }
    return s_def;
  }
} // namespace aurostd

// ***************************************************************************
// Functions args2vectorstring
// ***************************************************************************
namespace aurostd {
  vector<string> string2vstring(const string& str_in) {
    vector<string> str_out;
    str_out.push_back(str_in);
    return str_out;
  }

  vector<string> args2vectorstring(const vector<string>& argv, const string& s0, const string& s_def) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t i = 0; i < argv.size() - 1; i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i] == tokens[j]) {
          return vector<string>(argv.begin() + (i + 1), argv.end());
        }
      }
    }
    return string2vstring(s_def);
  }
} // namespace aurostd

// ***************************************************************************
// Functions get_itemized_vector_string stuff
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  bool get_itemized_vector_string_from_input(const vector<string>& argv, const string& s0, vector<string>& tokens, const string& delimiter) {// =":")
    if (aurostd::substring2bool(s0, "|")) {
      const string message = "not ported to \"|\"";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    uint icount = 0;
    string s0neq = s0;
    string s0equ;
    aurostd::StringSubstInPlace(s0neq, "=", "");
    s0equ = s0neq + "=";
    if (aurostd::args2attachedflag(argv, s0equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s0equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2flag(argv, s0neq)) {
      aurostd::args2string(argv, s0neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s0neq, EMPTY_WORDING);
    }
    if (tokens.size() == 1 && aurostd::substring2bool(tokens[0], delimiter)) {
      s0equ = tokens[0];
      icount += aurostd::string2tokens(s0equ, tokens, delimiter);
    }
    if (tokens.empty()) {
      return false;
    }
    if (icount == 0) {
      return false;
    }
    return true;
  }
  bool get_itemized_vector_string_from_input(const vector<string>& argv, const string& s0, const string& s1, vector<string>& tokens, const string& delimiter) {// =":")
    uint icount = 0;
    string s0neq = s0;
    string s0equ;
    aurostd::StringSubstInPlace(s0neq, "=", "");
    s0equ = s0neq + "=";
    string s1neq = s1;
    string s1equ;
    aurostd::StringSubstInPlace(s1neq, "=", "");
    s1equ = s1neq + "=";
    tokens.clear();
    if (tokens.empty() && aurostd::args2attachedflag(argv, s0equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s0equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (tokens.empty() && aurostd::args2attachedflag(argv, s1equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s1equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (tokens.empty() && aurostd::args2flag(argv, s0neq)) {
      aurostd::args2string(argv, s0neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s0neq, EMPTY_WORDING);
    }
    if (tokens.empty() && aurostd::args2flag(argv, s1neq)) {
      aurostd::args2string(argv, s1neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s1neq, EMPTY_WORDING);
    }
    if (tokens.size() == 1 && aurostd::substring2bool(tokens[0], delimiter)) {
      s0equ = tokens[0];
      icount += aurostd::string2tokens(s0equ, tokens, delimiter);
    }
    if (tokens.empty()) {
      return false;
    }
    if (icount == 0) {
      return false;
    }
    return true;
  }
  bool get_itemized_vector_string_from_input(const vector<string>& argv, const string& s0, const string& s1, const string& s2, vector<string>& tokens, const string& delimiter) {// =":")
    uint icount = 0;
    string s0neq = s0;
    string s0equ;
    aurostd::StringSubstInPlace(s0neq, "=", "");
    s0equ = s0neq + "=";
    string s1neq = s1;
    string s1equ;
    aurostd::StringSubstInPlace(s1neq, "=", "");
    s1equ = s1neq + "=";
    string s2neq = s2;
    string s2equ;
    aurostd::StringSubstInPlace(s2neq, "=", "");
    s2equ = s2neq + "=";
    if (aurostd::args2attachedflag(argv, s0equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s0equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2attachedflag(argv, s1equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s1equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2attachedflag(argv, s2equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s2equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2flag(argv, s0neq)) {
      aurostd::args2string(argv, s0neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s0neq, EMPTY_WORDING);
    }
    if (aurostd::args2flag(argv, s1neq)) {
      aurostd::args2string(argv, s1neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s1neq, EMPTY_WORDING);
    }
    if (aurostd::args2flag(argv, s2neq)) {
      aurostd::args2string(argv, s2neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s2neq, EMPTY_WORDING);
    }
    if (tokens.size() == 1 && aurostd::substring2bool(tokens[0], delimiter)) {
      s0equ = tokens[0];
      icount += aurostd::string2tokens(s0equ, tokens, delimiter);
    }
    if (tokens.empty()) {
      return false;
    }
    if (icount == 0) {
      return false;
    }
    return true;
  }
  bool get_itemized_vector_string_from_input(const vector<string>& argv, const string& s0, const string& s1, const string& s2, const string& s3, vector<string>& tokens, const string& delimiter) {// =":")
    uint icount = 0;
    string s0neq = s0;
    string s0equ;
    aurostd::StringSubstInPlace(s0neq, "=", "");
    s0equ = s0neq + "=";
    string s1neq = s1;
    string s1equ;
    aurostd::StringSubstInPlace(s1neq, "=", "");
    s1equ = s1neq + "=";
    string s2neq = s2;
    string s2equ;
    aurostd::StringSubstInPlace(s2neq, "=", "");
    s2equ = s2neq + "=";
    string s3neq = s3;
    string s3equ;
    aurostd::StringSubstInPlace(s3neq, "=", "");
    s3equ = s3neq + "=";
    if (aurostd::args2attachedflag(argv, s0equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s0equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2attachedflag(argv, s1equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s1equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2attachedflag(argv, s2equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s2equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2attachedflag(argv, s3equ)) {
      icount += aurostd::string2tokens(aurostd::args2attachedstring(argv, s3equ, EMPTY_WORDING), tokens, delimiter);
    }
    if (aurostd::args2flag(argv, s0neq)) {
      aurostd::args2string(argv, s0neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s0neq, EMPTY_WORDING);
    }
    if (aurostd::args2flag(argv, s1neq)) {
      aurostd::args2string(argv, s1neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s1neq, EMPTY_WORDING);
    }
    if (aurostd::args2flag(argv, s2neq)) {
      aurostd::args2string(argv, s2neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s2neq, EMPTY_WORDING);
    }
    if (aurostd::args2flag(argv, s3neq)) {
      aurostd::args2string(argv, s3neq, EMPTY_WORDING);
      tokens = aurostd::args2vectorstring(argv, s3neq, EMPTY_WORDING);
    }
    if (tokens.size() == 1 && aurostd::substring2bool(tokens[0], delimiter)) {
      s0equ = tokens[0];
      icount += aurostd::string2tokens(s0equ, tokens, delimiter);
    }
    if (tokens.empty()) {
      return false;
    }
    if (icount == 0) {
      return false;
    }
    return true;
  }
} // namespace aurostd

// ***************************************************************************
// args2attachedflag without/with commands
// ***************************************************************************
namespace aurostd {  // namespace aurostd

  bool args2attachedflag(const vector<string>& argv, const string& s0) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t i = 0; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (aurostd::substring2bool(argv[i], tokens[j])) {
          return true;
        }
      }
    }
    return false;
  }

  bool args2attachedflag(const vector<string>& argv, std::vector<string>& cmds, const string& s0) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveSpaces(s0);
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t j = 0; j < tokens.size(); j++) {
      cmds.push_back(tokens[j]);
    }
    for (size_t i = 0; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (aurostd::substring2bool(argv[i], tokens[j])) {
          return true;
        }
      }
    }
    return false;
  }
} // namespace aurostd

// ***************************************************************************
// args2attachedstring(argv,"-xxx=",abc) returns the content of xxx= or abc
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  string args2attachedstring(const vector<string>& argv, const string& s0, string s_def) { // string=""
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    const string s = aurostd::RemoveWhiteSpaces(s0);
    string output;
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      cerr << "argv.size()=" << argv.size() << endl;
      for (size_t j = 0; j < argv.size(); j++) {
        cerr << "[" << argv[j] << "]" << endl;
      }
    }
    if (VERBOSE) {
      cerr << "tokens.size()=" << tokens.size() << endl;
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t i = 1; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i].find(tokens[j]) != string::npos) {
          output = argv[i].substr(argv[i].find(tokens[j]) + tokens[j].length());
          if (!output.empty()) {
            return output;
          }
        }
      }
    }
    return s_def;
  }

  template <typename string> string args2attachedutype(const vector<string>& argv, const string& s0, const string& s_def) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    string s = aurostd::RemoveWhiteSpaces(s0);
    string output = "";
    vector<string> tokens;
    aurostd::string2tokens(s, tokens, "|");
    if (VERBOSE) {
      cerr << "argv.size()=" << argv.size() << endl;
      for (size_t j = 0; j < argv.size(); j++) {
        cerr << "[" << argv[j] << "]" << endl;
      }
    }
    if (VERBOSE) {
      cerr << "tokens.size()=" << tokens.size() << endl;
      for (size_t j = 0; j < tokens.size(); j++) {
        cerr << "[" << tokens[j] << "]" << endl;
      }
    }
    for (size_t i = 1; i < argv.size(); i++) {
      for (size_t j = 0; j < tokens.size(); j++) {
        if (argv[i].find(tokens[j]) != string::npos) {
          return argv[i].substr(argv[i].find(tokens[j]) + tokens[j].length());
        }
      }
    }
    return s_def;
  }
} // namespace aurostd

// ***************************************************************************
// args2attachedutype
// ***************************************************************************
namespace aurostd {  // namespace aurostd
  template <typename utype> utype args2attachedutype(const vector<string>& argv, const string& str1, const utype& utype_default) {
    const bool VERBOSE = (false || VERBOSE_ARGV); // DX20200907 - LDEBUG to VERBOSE; decouple from XHOST.DEBUG
    vector<string> tokens1;
    aurostd::string2tokens(aurostd::RemoveWhiteSpaces(str1), tokens1, "|");
    if (VERBOSE) {
      cerr << "argv.size()=" << argv.size() << endl;
      for (size_t j = 0; j < argv.size(); j++) {
        cerr << "[" << argv[j] << "]" << endl;
      }
    }
    if (VERBOSE) {
      cerr << "tokens1.size()=" << tokens1.size() << endl;
      for (size_t j = 0; j < tokens1.size(); j++) {
        cerr << "[" << tokens1[j] << "]" << endl;
      }
    }
    utype out = utype_default;
    for (size_t j = 0; j < tokens1.size(); j++) {
      string s1 = tokens1[j];
      string s1eq;
      string s1neq;     //   s1=aurostd::RemoveSubString(s1,"-");s1=aurostd::RemoveSubString(s1,"-");
      s1 = aurostd::RemoveSubString(s1, "=");
      s1eq = s1 + "=";
      s1neq = s1;// cerr << s1eq << " " << s1neq << endl;
      if (aurostd::args2flag(argv, s1neq) || aurostd::args2attachedflag(argv, s1eq)) {
        if (aurostd::args2flag(argv, s1neq)) {
          out = aurostd::args2utype(argv, s1neq, out);
        }
        if (aurostd::args2attachedflag(argv, s1eq)) {
          vector<string> tokens2;
          get_itemized_vector_string_from_input(argv, s1neq, tokens2, ",");
          if (!tokens2.empty()) {
            out = aurostd::string2utype<utype>(tokens2[0]);
          }
        }
      }
    }
    return out;
  }
#define AST_TEMPLATE(atype) template atype args2attachedutype(const vector<string>&, const string&, const atype&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
} // namespace aurostd

// *******************************************************************************************
// *******************************************************************************************

#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
