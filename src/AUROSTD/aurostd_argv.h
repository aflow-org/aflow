// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_ARGV_H_
#define _AUROSTD_ARGV_H_

#include <deque>
#include <string>
#include <vector>

#include "aurostd_defs.h"
#include "aurostd_xvector.h"

// ***************************************************************************
// GET WORLD
namespace aurostd {
  // ATTACH
  std::string attach(const std::string&) __xprototype;
  std::string attach(const std::string&, const std::string&) __xprototype;
  std::string attach(const std::string&, const std::string&, const std::string&) __xprototype;
  std::string attach(const std::string&, const std::string&, const std::string&, const std::string&) __xprototype;
  std::string attach(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&) __xprototype;
  // get_arguments_from_input
  std::vector<std::string> get_arguments_from_input(int argc, char** argv) __xprototype;
  // get_flag without/with command list
  bool args2flag(const std::vector<std::string>& argv, const std::string&) __xprototype;
  bool args2flag(const std::vector<std::string>& argv, std::vector<std::string>&, const std::string&) __xprototype;
  // args2utype of utype type
  template <class utype> utype args2utype(const std::vector<std::string>& argv, const std::string&, utype) __xprototype;
  // get_xvector get_vector/deque
  template <class utype> xvector<utype> args2xvectorutype(const std::vector<std::string>& argv, const std::string&, const xvector<utype>&) __xprototype;
  template <class utype> xvector<utype> args2xvectorutype(const std::vector<std::string>& argv, const std::string&, int) __xprototype;
  template <class utype> std::vector<utype> args2vectorutype(const std::vector<std::string>& argv, const std::string&) __xprototype;
  template <class utype> std::deque<utype> args2dequeutype(const std::deque<std::string>& argv, const std::string&) __xprototype;
  // args2vectorstring
  std::string args2string(const std::vector<std::string>& argv, const std::string&, const std::string&) __xprototype;
  std::string args2string(const std::vector<std::string>& argv, std::vector<std::string>&, const std::string&, const std::string&) __xprototype;
  // args2vectorstring
  std::vector<std::string> args2vectorstring(const std::vector<std::string>& argv, const std::string&, const std::string&) __xprototype;

  //__get_itemized_vector_string
  bool get_itemized_vector_string_from_input(const std::vector<std::string>& argv, const std::string& s0, std::vector<std::string>& tokens, const std::string& delimiter) __xprototype;
  bool get_itemized_vector_string_from_input(const std::vector<std::string>& argv, const std::string& s0, const std::string& s1, std::vector<std::string>& tokens, const std::string& delimiter) __xprototype;
  bool get_itemized_vector_string_from_input(const std::vector<std::string>& argv, const std::string& s0, const std::string& s1, const std::string& s2, std::vector<std::string>& tokens, const std::string& delimiter) __xprototype;
  bool get_itemized_vector_string_from_input(const std::vector<std::string>& argv, const std::string& s0, const std::string& s1, const std::string& s2, const std::string& s3, std::vector<std::string>& tokens, const std::string& delimiter)
      __xprototype;

  // args2attachedflag without/with commands
  bool args2attachedflag(const std::vector<std::string>&, const std::string&) __xprototype;
  // args2attachedflag with commands
  bool args2attachedflag(const std::vector<std::string>&, std::vector<std::string>&, const std::string&) __xprototype;
  // args2attachedstring
  std::string args2attachedstring(const std::vector<std::string>&, const std::string&, std::string = "") __xprototype;
  // args2attachedint
  // args2attacheddouble
  // args2attachedutype
  template <typename utype> utype args2attachedutype(const std::vector<std::string>& argv, const std::string&, const utype&) __xprototype;
  std::string args2attachedutype(const std::vector<std::string>& argv, const std::string&, const std::string&) __xprototype;
} // namespace aurostd

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
