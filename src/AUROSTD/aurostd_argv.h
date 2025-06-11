// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_ARGV_H_
#define _AUROSTD_ARGV_H_

#include "aurostd.h"
//#include "aflow.h"

// ***************************************************************************
// GET WORLD
namespace aurostd {
  // ATTACH
  string attach(const string&) __xprototype;
  string attach(const string&,const string&) __xprototype;
  string attach(const string&,const string&,const string&) __xprototype;
  string attach(const string&,const string&,const string&,const string&) __xprototype;
  string attach(const string&,const string&,const string&,const string&,const string&) __xprototype;
  // get_arguments_from_input
  vector<string> get_arguments_from_input(int argc,char **argv) __xprototype;
  // get_flag without/with command list
  bool args2flag(const vector<string>& argv,const string&) __xprototype;
  bool args2flag(const vector<string>& argv,vector<string>&,const string&) __xprototype;
  // args2utype of utype type
  template<class utype> utype args2utype(const vector<string>& argv,const string&,utype) __xprototype;
  // get_xvector get_vector/deque
  template<class utype> xvector<utype> args2xvectorutype(const vector<string>& argv,const string&, const xvector<utype>&) __xprototype;
  template<class utype> xvector<utype> args2xvectorutype(const vector<string>& argv,const string&,int) __xprototype;
  template<class utype> vector<utype> args2vectorutype(const vector<string>& argv,const string&) __xprototype;
  template<class utype> deque<utype> args2dequeutype(const deque<string>& argv,const string&) __xprototype;
  // args2vectorstring
  string args2string(const vector<string>& argv,const string&,const string&) __xprototype;
  string args2string(const vector<string>& argv,vector<string>&,const string&,const string&) __xprototype;
  // args2vectorstring
  vector<string> args2vectorstring(const vector<string>& argv,const string&,const string&) __xprototype;

  //__get_itemized_vector_string
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter) __xprototype;
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter) __xprototype;
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter) __xprototype;
  bool get_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter) __xprototype;
  //[CO20200624 - OBSOLETE]// getproto_itemized_vector_string_from_input
  //[CO20200624 - OBSOLETE]bool getproto_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,vector<string>& tokens,const string& delimiter=":") __xprototype;
  //[CO20200624 - OBSOLETE]bool getproto_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,vector<string>& tokens,const string& delimiter=":") __xprototype;
  //[CO20200624 - OBSOLETE]bool getproto_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,const string& s2,vector<string>& tokens,const string& delimiter=":") __xprototype;
  //[CO20200624 - OBSOLETE]bool getproto_itemized_vector_string_from_input(const vector<string> &argv,const string& s0,const string& s1,const string& s2,const string& s3,vector<string>& tokens,const string& delimiter=":") __xprototype;

  // args2attachedflag without/with commands
  bool args2attachedflag(const vector<string>&,const string&) __xprototype;
  // args2attachedflag with commands
  bool args2attachedflag(const vector<string>&,vector<string>&,const string&) __xprototype;
  // args2attachedstring
  // [OBSOLETE] string args2attachedstring(vector<string>,const string&) __xprototype;
  // [OBSOLETE] string args2attachedstring(vector<string>,const string&,const string&) __xprototype;
  string args2attachedstring(const vector<string>&,const string&,string="") __xprototype;
  // [OBSOLETE] string args2attachedstring(vector<string>,const string&,const string&,const string&) __xprototype;
  // [OBSOLETE] string args2attachedstring(vector<string>,const string&,const string&,const string&,const string&) __xprototype;
  // [OBSOLETE] string args2attachedstring(vector<string>,const string&,const string&,const string&,const string&,const string&) __xprototype;
  // args2attachedint
  // [OBSOLETE] int args2attachedint(const vector<string>& argv,const string&,const int&) __xprototype;
  // [OBSOLETE] int args2attachedint(const vector<string>& argv,const string&,const string&,const int&) __xprototype;
  // args2attacheddouble
  // [OBSOLETE] double args2attacheddouble(const vector<string>& argv,const string&,const double&) __xprototype;
  // [OBSOLETE] double args2attacheddouble(const vector<string>& argv,const string&,const string&,const double&) __xprototype;
  // args2attachedutype
  template<typename utype> utype args2attachedutype(const vector<string>& argv,const string&,const utype&) __xprototype;
  // [OBSOLETE] template<typename utype> utype args2attachedutype(const vector<string>& argv,const string&,const string&,const utype&) __xprototype;
  // [OBSOLETE] template<typename utype> utype args2attachedutype(const vector<string>& argv,const string&,const string&,const string&,const utype&) __xprototype;
  // [OBSOLETE] template<typename utype> utype args2attachedutype(const vector<string>& argv,const string&,const string&,const string&,const string&,const utype&) __xprototype;
  string args2attachedutype(const vector<string>& argv,const string&,const string&) __xprototype;
  // [OBSOLETE] string args2attachedutype(const vector<string>& argv,const string&,const string&,const string&) __xprototype;
  // [OBSOLETE] string args2attachedutype(const vector<string>& argv,const string&,const string&,const string&,const string&) __xprototype;
  // [OBSOLETE] string args2attachedutype(const vector<string>& argv,const string&,const string&,const string&,const string&,const string&) __xprototype;
}

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
