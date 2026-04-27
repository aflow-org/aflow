// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_MAIN_H_
#define _AUROSTD_MAIN_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include "aurostd_defs.h"
#include "aurostd_xcomplex.h"
#include "aurostd_xerror.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xoption.h"
#include "aurostd_xvector.h"

namespace fs = std::filesystem;

#ifndef uint
typedef unsigned uint;
#endif

// todo ST: this should be moved somewhere dedicated for util or logging
#define __AFLOW_FUNC__ aflowFunc(__PRETTY_FUNCTION__, __func__)
#define __AFLOW_FILE__ aflowFile(__FILE__, __LINE__)

// ME+HE20220321
// Based on https://stackoverflow.com/questions/1666802/is-there-a-class-macro-in-c
// Get full formatted name of function
inline std::string aflowFunc(const std::string& pretty_func, const std::string& func) {
  const size_t end = pretty_func.find(func);
  // Everything between the function name and the last space character
  // are namespace and class name, if present
  const size_t begin = pretty_func.rfind(" ", end) + 1;
  return aurostd::xerror_PID + pretty_func.substr(begin, end - begin) + func + "():";
}

// SD+HE20220914
// Get filename and line number
inline std::string aflowFile(const std::string& file_name, const size_t line_number) {
  return file_name + ":" + std::to_string(line_number);
}

// ----------------------------------------------------------------------------

// TODO an implementation for these is in aflow_xproto.cpp, should be fixed
template <class utype> std::ostream& operator<<(std::ostream&, const std::vector<utype>&);// __xprototype;
template <class utype> std::ostream& operator<<(std::ostream&, const std::deque<utype>&);// __xprototype;

// ----------------------------------------------------------------------------
// threadID stuff
namespace aurostd {
  unsigned long long int getTID(); // CO20200502 - threadID
}
// ----------------------------------------------------------------------------

namespace aurostd {
  void sizes() __xprototype;
  // aflow_aurostd.cpp
  // int Nint(double x);
  // int Sign(const double& x);
  // int SignNoZero(const double& x);
  template <class utype> void aswap(utype& a, utype& b);                     // some dumb algebra
  template <class utype> utype max(const std::vector<utype> vec);           // some dumb algebra
  template <class utype> utype max(const std::deque<utype> vec);            // some dumb algebra
  template <class utype> utype max(const std::vector<std::vector<utype>> mat);   // some dumb algebra
  template <class utype> utype min(const std::vector<utype> vec);           // some dumb algebra
  template <class utype> utype min(const std::deque<utype> vec);            // some dumb algebra
  template <class utype> utype min(const std::vector<std::vector<utype>> mat);   // some dumb algebra
  template <class utype> utype sum(const std::vector<utype> vec);           // some dumb algebra
  template <class utype> utype sum(const std::deque<utype> vec);            // some dumb algebra
  template <class utype> utype sum(const std::vector<std::vector<utype>> mat);   // some dumb algebra
  template <class utype> utype mean(const std::vector<utype> vec);          // some dumb algebra
  template <class utype> utype mean(const std::deque<utype> vec);           // some dumb algebra
  template <class utype> utype mean(const std::vector<std::vector<utype>> mat);  // some dumb algebra
  template <class utype> std::vector<utype> reset(std::vector<utype>& v);        // some dumb algebra
  template <class utype> std::deque<utype> reset(std::deque<utype>& v);          // some dumb algebra
  template <class utype> std::vector<std::vector<utype>> reset(std::vector<std::vector<utype>> m);  // some dumb algebra
  template <class utype> std::vector<utype> clear(std::vector<utype>& v);        // some dumb algebra
  template <class utype> std::deque<utype> clear(std::deque<utype>& v);          // some dumb algebra
  template <class utype> std::vector<std::vector<utype>> clear(std::vector<std::vector<utype>> m);  // some dumb algebra
  template <class utype> void random_shuffle(std::vector<utype>& v);        // some dumb algebra
  template <class utype> void random_shuffle(std::deque<utype>& v);         // some dumb algebra
  template <class utype> std::vector<utype> getEveryNth(const std::vector<utype>& vec, uint n); // SD20230216
  std::vector<double> getEveryNth(const std::vector<double>& vec, uint n);
  template <class utype> bool identical(std::vector<utype> v1, std::vector<utype> v2, utype epsilon);
  template <class utype> bool identical(std::deque<utype> v1, std::deque<utype> v2, utype epsilon);
  bool identical(std::vector<int> v1, std::vector<int> v2, int epsilon);
  bool identical(std::deque<int> v1, std::deque<int> v2, int epsilon);
  template <class utype> bool identical(const std::vector<std::vector<utype>>& m1, const std::vector<std::vector<utype>>& m2, utype epsilon);
  template <class utype> bool identical(const std::vector<utype>& vec, utype eps = (utype) AUROSTD_IDENTITY_TOL); // DX20210422 - checks if all values are the same
  template <class utype> bool identical(const std::deque<utype>& vec, utype eps = (utype) AUROSTD_IDENTITY_TOL); // DX20210422 - checks if all values are the same
  std::string toupper(const std::string& in) __xprototype;
  std::string tolower(const std::string& in) __xprototype;
  char toupper(const char& in) __xprototype;
  char tolower(const char& in) __xprototype;
  int GetNumFields(const std::string& s);
  std::string GetNextVal(const std::string& s, int& id);
  std::string PaddedNumString(const int num, const int ndigits);
  int getZeroPadding(double num);  // CO20191217
  int getZeroPadding(int num);  // CO20191217
  int getZeroPadding(uint num); // CO20191217
  int getZeroPadding(long int num); // CO20191217
  int getZeroPadding(unsigned long int num);  // CO20191217
  int getZeroPadding(long long int num);  // CO20191217
  int getZeroPadding(unsigned long long int num);  // ME20190108
  template <class utype> std::string PaddedPRE(utype, int, std::string = " ");
  std::string PaddedPRE(std::string, int, std::string = " ");
  template <class utype> std::string PaddedPOST(utype, int, std::string = " ");
  std::string PaddedPOST(std::string, int, std::string = " ");
  template <class utype> std::string PaddedCENTER(utype, int, std::string = " ");
  // write progresses
  std::string PaddedCENTER(std::string, int, std::string = " ");
  uint ProgressBar(std::ostream& oss, std::string prelim, uint j, uint jmax, bool VERBOSE_PERCENTAGE, bool VERBOSE_ROLLER, bool VERBOSE_CURSOR);
  uint ProgressBar(std::ostream& oss, std::string prelim, uint j, uint jmax);
  uint ProgressBar(std::ostream& oss, std::string prelim, double j, bool VERBOSE_PERCENTAGE, bool VERBOSE_ROLLER, bool VERBOSE_CURSOR);
  uint ProgressBar(std::ostream& oss, std::string prelim, double j);
  // about cleaning up strings
  bool RemoveControlCodeCharactersFromString(const std::string& in, std::string& out); // DX20190516  //CO20190620
  bool RemoveControlCodeCharactersFromStringstream(std::stringstream& ss_in, std::stringstream& ss_out); // DX20190516
  bool RemoveControlCodeCharactersFromFile(const std::string& directory, const std::string& filename, bool keep_orig_file = true); // DX20190516
  bool isNullByte(char c); // DX20190131
  std::string removeNullBytes(std::string in); // DX20190131
  bool RemoveBinaryCharactersFromFile(const std::string& directory, const std::string& filename); // DX20190211 //CO20210315
  std::string PercentEncodeASCII(const char c) __xprototype; // DX20210706
  std::string CleanStringASCII(const std::string& s) __xprototype;
  void CleanStringASCII_InPlace(std::string& s) __xprototype;  // CO20190712
  std::string RemoveTrailingCharacter(const std::string& s, char c); // CO+ME20200825
  void RemoveTrailingCharacter_InPlace(std::string& s, char c); // CO+ME20200825
  std::string CGI_StringClean(const std::string& stringIN) __xprototype;
  std::string RemoveWhiteSpaces(const std::string& s) __xprototype;
  std::string RemoveWhiteSpaces(const std::string& s, const char toogle) __xprototype;
  std::string RemoveWhiteSpacesFromTheBack(const std::string& s) __xprototype;
  std::string RemoveWhiteSpacesFromTheFront(const std::string& s) __xprototype;
  std::string RemoveWhiteSpacesFromTheFrontAndBack(const std::string& s) __xprototype;
  std::string RemoveSpaces(const std::string& s) __xprototype;
  std::string RemoveSpaces(const std::string& s, const char toogle) __xprototype;
  std::string RemoveSpacesFromTheBack(const std::string& s) __xprototype;
  std::string RemoveTabs(const std::string& s) __xprototype;
  std::string RemoveTabs(const std::string& s, const char toogle) __xprototype;
  std::string RemoveTabsFromTheBack(const std::string& s) __xprototype;
  std::string RemoveComments(const std::string& s) __xprototype;
  std::vector<std::string> RemoveComments(const std::vector<std::string>&) __xprototype;  // ME20190614
  std::deque<std::string> RemoveComments(const std::deque<std::string>&) __xprototype;  // ME20190614
  std::string RemoveCharacter(const std::string& s, const char character) __xprototype;
  void RemoveCharacterInPlace(std::string& s, const char character) __xprototype;  // CO20190712
  std::string RemoveCharacterFromTheBack(const std::string& s, const char character);
  __xprototype; // DX20190708
  std::string RemoveCharacterFromTheFront(const std::string& s, const char character);
  __xprototype; // DX20190708
  std::string RemoveCharacterFromTheFrontAndBack(const std::string& s, const char character);
  __xprototype; // DX20190708
  std::string RemoveNumbers(const std::string& s) __xprototype; // CO20190712
  void RemoveNumbersInPlace(std::string& s) __xprototype;  // CO20190712
  std::string RemoveRounding(const std::string& s) __xprototype;
  // std::string RemoveCharacter(const std::string& s, const char character, const char toogle) __xprototype;
  std::string RemoveSubStringFirst(const std::string& str_orig, const std::string& str_rm) __xprototype;
  void RemoveSubStringFirstInPlace(std::string& str_orig, const std::string& str_rm) __xprototype;  // CO20190712
  std::string RemoveSubString(const std::string& str_orig, const std::string& str_rm) __xprototype;
  void RemoveSubStringInPlace(std::string& str_orig, const std::string& str_rm) __xprototype; // CO20190712
  double VersionString2Double(const std::string& version_str); // SD20220331
  std::vector<std::string> ProcessPIDs(const std::string& process, bool user_specific = true); // CO20210315
  std::vector<std::string> ProcessPIDs(const std::string& process, std::string& output_syscall, bool user_specific = true); // CO20210315
  std::vector<std::string> ProcessPIDs(const std::string& process, const std::string& pgid, std::string& output_syscall, bool user_specific = true); // SD20220329
  bool ProcessRunning(const std::string& process, bool user_specific = true); // CO20210315
  bool ProcessRunning(const std::string& process, const std::string& pgid, bool user_specific = true); // SD20220329
  void ProcessKill(const std::string& process, bool user_specific = true, uint signal = 9); // CO20210315 //SD20220627 - 9 = SIGKILL, 15 = SIGTERM
  void ProcessKill(const std::string& process, const std::string& pgid, bool user_specific = true, uint signal = 9); // SD20220329 - 9 = SIGKILL, 15 = SIGTERM
  bool ProcessRenice(const std::string& process, int nvalue, bool user_specific = true, const std::string& pgid = ""); // CO20210315
  bool ReniceAvailable(); // CO20221029

  // about printing
  void PrintANSIEscapeSequence(const aurostd::xoption& color, FILE* fstr);
  void PrintMessageStream(std::ostringstream& stream, bool quiet, std::ostream& oss = std::cout); // CO20200624
  void PrintMessageStream(std::ofstream& FileMESSAGE, std::ostringstream& stream, bool quiet, std::ostream& oss = std::cout); // CO20200624
  void PrintMessageStream(std::ofstream& FileMESSAGE, std::ostringstream& stream, bool quiet, bool osswrite, std::ostream& oss = std::cout);
  void PrintErrorStream(std::ostringstream& stream, bool quiet); // CO20200624
  void PrintErrorStream(std::ofstream& FileERROR, std::ostringstream& stream, bool quiet); // CO20200624
  void PrintErrorStream(std::ofstream& FileERROR, std::ostringstream& stream, bool quiet, bool osswrite);
  void PrintWarningStream(std::ostringstream& stream, bool quiet); // CO20200624
  void PrintWarningStream(std::ofstream& FileWARNING, std::ostringstream& stream, bool quiet); // CO20200624
  void PrintWarningStream(std::ofstream& FileWARNING, std::ostringstream& stream, bool quiet, bool osswrite);

  void PrintMessageStream(std::stringstream& stream, bool quiet, std::ostream& oss = std::cout);  // CO20200624
  void PrintMessageStream(std::ofstream& FileMESSAGE, std::stringstream& stream, bool quiet, std::ostream& oss = std::cout);  // CO20200624
  void PrintMessageStream(std::ofstream& FileMESSAGE, std::stringstream& stream, bool quiet, bool osswrite, std::ostream& oss = std::cout);
  void PrintErrorStream(std::stringstream& stream, bool quiet);  // CO20200624
  void PrintErrorStream(std::ofstream& FileERROR, std::stringstream& stream, bool quiet);  // CO20200624
  void PrintErrorStream(std::ofstream& FileERROR, std::stringstream& stream, bool quiet, bool osswrite);
  void PrintWarningStream(std::stringstream& stream, bool quiet);  // CO20200624
  void PrintWarningStream(std::ofstream& FileWARNING, std::stringstream& stream, bool quiet);  // CO20200624
  void PrintWarningStream(std::ofstream& FileWARNING, std::stringstream& stream, bool quiet, bool osswrite);

  // about executing

  bool IsCommandAvailable(const std::string& command);
  bool IsCommandAvailable(const std::string& command, std::string& position);
  bool IsCommandAvailableModify(std::string& command);
  bool CommandRequired(const std::string& command);
  bool CommandRequired(const std::string& command, std::string& position);
  bool IsExecutableAvailable(const std::string& executable);
  bool IsExecutableAvailable(const std::string& executable, std::string& position);
  bool ExecutableRequired(const std::string& executable);
  bool ExecutableRequired(const std::string& executable, std::string& position);

  bool execute(std::ostringstream& command);
  bool execute(std::stringstream& command);
  bool execute(const std::string& command);
  bool execute(const std::vector<std::string>& vcommand);
  bool execute(const std::deque<std::string>& dcommand);
  bool execute_thread_safe(const std::string& command); // HE20240220

  // about cleaning, higiene is important
  void StringstreamClean(std::ostringstream& stream);
  void StringstreamClean(std::stringstream& stream);
  int FindIfStringInStream(const std::string& key, std::istream& instream);

  bool GetMemoryUsagePercentage(double& usage_percentage_ram, double& usage_percentage_swap);  // CO20210601
  bool GetMemory(unsigned long long int& free_ram, unsigned long long int& total_ram, unsigned long long int& free_swap, unsigned long long int& total_swap); // CO20210315

#ifdef _stringcharstar_
  bool execute(char* command);
#endif
  // Execute and report
  std::pair<std::string, std::string> execute2OutErrPair(const std::string& command); // HE20240220
  std::string execute2string(std::ostringstream& command);
  std::string execute2string(std::stringstream& command);
  std::string execute2string(const std::string& command);
  std::vector<std::string> execute2string(const std::vector<std::string>& vcommand);
  std::deque<std::string> execute2string(const std::deque<std::string>& dcommand);
#ifdef _stringcharstar_
  std::string execute2string(char* command);
#endif
  std::string CleanCommand4Execute(const std::string& command); // CO20200624
  template <class utype> utype execute2utype(std::ostringstream& command);
  template <class utype> utype execute2utype(std::stringstream& command);
  template <class utype> utype execute2utype(std::string command);
  template <class utype> std::vector<utype> execute2utype(std::vector<std::string> vcommand);
  template <class utype> std::deque<utype> execute2utype(std::deque<std::string> dcommand);
#ifdef _stringcharstar_
  template <class utype> utype execute2utype(char* command);
#endif
  // about sleeping
  unsigned int Sleep(unsigned int seconds);
  // about extracting from to files
  std::vector<std::string> GrepFile(const std::string& filename, const std::string& keyword, bool RemoveWS = false, bool RemoveComments = true); // CO20210623
  // take just after
  bool ExtractJustAfterToStringstreamEXPLICIT(std::ifstream& FileIN, std::stringstream& StringstreamOUTPUT, const std::string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(std::stringstream& StringStreamIN, std::stringstream& StringstreamOUTPUT, const std::string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(const std::string& StringIN, std::stringstream& StringstreamOUTPUT, const std::string& Keyword_start);
  bool ExtractJustAfterToFileEXPLICIT(std::ifstream& FileIN, const std::string& FileNameOUTPUT, const std::string& Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(std::ifstream& FileIN, std::string& StringOUTPUT, const std::string& Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(const std::string& StringIN, std::string& StringOUTPUT, const std::string& Keyword_start);
  // about taking in istreams and std::stringstream and strings
  size_t stream2vectorstring(std::istream& istreamIN, std::vector<std::string>& vstringout);
  size_t stream2vectorstring(std::ifstream& ifstreamIN, std::vector<std::string>& vstringout);
  size_t stream2vectorstring(std::stringstream& stringstreamIN, std::vector<std::string>& vstringout);
  size_t string2vectorstring(const std::string& stringIN, std::vector<std::string>& vstringout, bool consecutive = false, bool trim_edges = false); // CO20170613, defaults to usual string2tokens() behavior
  std::vector<std::string> stream2vectorstring(std::istream& istreamIN);
  std::vector<std::string> stream2vectorstring(std::ifstream& ifstreamIN);
  std::vector<std::string> stream2vectorstring(std::stringstream& stringstreamIN);
  std::vector<std::string> string2vectorstring(const std::string& stringIN, bool consecutive = false, bool trim_edges = false);  // CO20170613, defaults to usual string2tokens() behavior
  std::string liststring2string(std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "",
                                std::string = "");
  uint stream2dequestring(std::istream& istreamIN, std::deque<std::string>& vstringout);
  uint stream2dequestring(std::ifstream& ifstreamIN, std::deque<std::string>& vstringout);
  uint stream2dequestring(std::stringstream& stringstreamIN, std::deque<std::string>& vstringout);
  uint string2dequestring(const std::string& stringIN, std::deque<std::string>& vstringout);
  std::deque<std::string> stream2dequestring(std::istream& istreamIN);
  std::deque<std::string> stream2dequestring(std::ifstream& ifstreamIN);
  std::deque<std::string> stream2dequestring(std::stringstream& stringstreamIN);
  std::deque<std::string> string2dequestring(const std::string& stringIN);

  // about std::istream/std::ostream
  std::string ostream2string(std::ostream& oss);
  uint stream2string(std::istream& istreamIN, std::string& vstringout);
  uint stream2string(std::ifstream& ifstreamIN, std::string& vstringout);
  uint stream2string(std::stringstream& stringstreamIN, std::string& vstringout);
  // about environments
  std::string getenv2string(const std::string& str);
  int getenv2int(const std::string& str);
  uint getenv2uint(const std::string& str);
  double getenv2double(const std::string& str);
  bool Chmod(const uint chmod, const std::string& path, const std::filesystem::perm_options perm_opt = std::filesystem::perm_options::replace);
  bool ChmodRecursive(const uint chmod_dir,
                      const uint chmod_file,
                      const std::string& directory,
                      const std::filesystem::perm_options perm_opt_dir = std::filesystem::perm_options::replace,
                      const std::filesystem::perm_options perm_opt_file = std::filesystem::perm_options::replace);

  // about getting info from strings
  uint string2tokens(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ", bool consecutive = false) __xprototype; // CO20170613, defaults to usual string2tokens() behavior
  uint string2tokens(const std::string& str, std::deque<std::string>& tokens, const std::string& delimiters = " ", bool consecutive = false) __xprototype; // CO20170613, defaults to usual string2tokens() behavior
  template <class utype>
  uint string2tokens(const std::string& str, std::vector<utype>& tokens, const std::string& delimiters = " ", bool consecutive = false) __xprototype; // CO20170613, defaults to usual string2tokens() behavior
  template <class utype>
  uint string2tokens(const std::string& str, std::deque<utype>& tokens, const std::string& delimiters = " ", bool consecutive = false) __xprototype; // CO20170613, defaults to usual string2tokens() behavior
  uint string2tokensAdd(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ") __xprototype;
  uint string2tokensAdd(const std::string& str, std::deque<std::string>& tokens, const std::string& delimiters = " ") __xprototype;
  template <class utype> uint string2tokensAdd(const std::string& str, std::vector<utype>& tokens, const std::string& delimiters = " ") __xprototype;
  template <class utype> uint string2tokensAdd(const std::string& str, std::deque<utype>& tokens, const std::string& delimiters = " ") __xprototype;
  uint string2tokensByDelimiter(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiter); // SD20220504
  uint string2tokensByDelimiter(const std::string& str, std::deque<std::string>& tokens, const std::string& delimiter); // SD20220504

  //[CO20210315 - not defined]template<typename typeTo, typename typeFrom> typeTo NumberStreamConvert(const typeFrom& from);  //CO20210315 - cleaned up

  template <typename utype> std::vector<utype> vectorstring2vectorutype(const std::vector<std::string>& from); // SD20220520
  std::vector<double> vectorstring2vectordouble(const std::vector<std::string>& from); // CO20210315 - cleaned up
  std::string string2string(const std::string& from) __xprototype;
  template <typename utype> utype string2utype(const std::string& from, uint base = 10); // CO20210315 - cleaned up //HE20220324 add base option
  std::vector<int> vectorstring2vectorint(const std::vector<std::string>& from); // CO20210315 - cleaned up
  std::vector<uint> vectorstring2vectoruint(const std::vector<std::string>& from); // CO20210315 - cleaned up

  std::vector<float> vectorstring2vectorfloat(const std::vector<std::string>& from); // CO20210315 - cleaned up
  std::string vectorstring2string(const std::vector<std::string>& vstrings);
  std::string vectorstring2string(const std::deque<std::string>& vstrings);

  template <typename utype>
  std::string utype2string(const utype& from, int precision = AUROSTD_DEFAULT_PRECISION, char FORMAT = DEFAULT_STREAM) __xprototype; // DX20201028 - this declaration was missing //DX20210128 - add defaults
  std::string utype2string(double from, bool roff);
  std::string utype2string(double from, int precision, bool roff);
  std::string utype2string(double from, bool roff, double tol);
  std::string utype2string(double from, int precision, bool roff, double tol);
  std::string utype2string(double from, bool roff, char FORMAT);
  std::string utype2string(double from, int precision, char FORMAT, bool roff = false); // CO20200624
  std::string utype2string(double from, int precision, bool roff, char FORMAT);
  std::string utype2string(double from, bool roff, double tol, char FORMAT);
  std::string utype2string(double from, int precision, bool roff, double tol, char FORMAT);
  std::string bool2string(bool from);

  template <class utype> std::deque<utype> utypes2deque(utype u1) __xprototype;
  template <class utype> std::deque<utype> utypes2deque(utype u1, utype u2) __xprototype;
  template <class utype> std::deque<utype> utypes2deque(utype u1, utype u2, utype u3) __xprototype;
  template <class utype> std::deque<utype> utypes2deque(utype u1, utype u2, utype u3, utype u4) __xprototype;

  void StringCommasColumsVectorInt(std::string vstring, std::vector<int>& vint) __xprototype;
  void StringCommasColumsVectorUnsignedInt(std::string vstring, std::vector<uint>& vuint) __xprototype;
  void StringCommasColumsVectorFloat(std::string vstring, std::vector<float>& vfloat) __xprototype;
  void StringCommasColumsVectorDouble(std::string vstring, std::vector<double>& vdouble) __xprototype;
  size_t GetNLinesString(const std::string& str) __xprototype;
  size_t GetNLinesString(const std::stringstream& strstream) __xprototype;
  size_t GetNLinesFile(const std::string& file_name) __xprototype;
  std::string GetLineString(const std::string& strstream, int line);
  std::string GetLineString(const std::stringstream& strstream, int line);
  // substitute strings in strings and stringstreams
  bool StringsAlphabetic(const std::string& A, const std::string& B, bool allow_identical = true); // CO20180801
  bool StringsAlphabetic(const std::vector<std::string>& input, bool allow_identical = true); // CO20180801
  bool StringsAlphabetic(const std::deque<std::string>& input, bool allow_identical = true); // CO20180801
  void StringSubstInPlace(std::string& work_string, const std::string& old_string, const std::string& new_string); // HE20240908
  void StringSubstInPlace(std::string& work_string, char old_char, char new_char); // HE20240908
  std::string StringSubst(const std::string& work_string, const std::string& old_string, const std::string& new_string); // HE20220321
  std::string StringSubst(const std::string& work_string, char old_char, char new_char); // HE20220321
  void StringStreamSubst(std::stringstream& strstring, const std::string& strfind, const std::string& strreplace); // ME20190128 - fixed type declaration
  // about present substrings
  bool substring2bool(const std::string& strstream, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up
  bool substring2bool(const std::vector<std::string>& vstrstream, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up
  bool substring2bool(const std::deque<std::string>& vstrstream, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up
  bool substring2bool(const std::stringstream& strstream, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up
  bool substringlist2bool(const std::string& strin, const std::vector<std::string>& substrings, bool match_all = true); // ME20220505
  bool substringlist2bool(const std::string& strin, const std::deque<std::string>& substrings, bool match_all = true); // ME20220505
  bool substring_present_file(const std::string& FileName, const std::boyer_moore_searcher<std::string::const_iterator>& boyer_moore, const size_t strsub_size, size_t& offset);
  bool substring_present_file(const std::string& FileName, const std::string& strsub, size_t& offset);
  bool substring_present_file(const std::string& FileName, const std::string& strsub);
  std::vector<bool> substrings_present_file(const std::string& FileName, const std::vector<std::boyer_moore_searcher<std::string::const_iterator>>& vboyer_moore, const size_t strsub_size, size_t& offset);
  std::vector<bool> substrings_present_file(const std::string& FileName, const std::vector<std::string>& vstrsub, size_t& offset);
  std::vector<bool> substrings_present_file(const std::string& FileName, const std::vector<std::string>& vstrsub);
  void substrings_map_present_file(const std::string& FileName, std::map<std::string, std::pair<std::boyer_moore_searcher<std::string::const_iterator>, bool>>& map_boyer_moore, const size_t strsub_size, size_t& offset);
  void substrings_map_present_file(const std::string& FileName, std::unordered_map<std::string, bool>& map_strsub, size_t& offset);
  void substrings_map_present_file(const std::string& FileName, std::unordered_map<std::string, bool>& map_strsub);
  template <class utype> bool WithinList(const std::vector<utype>& list, const utype& input, size_t& index, bool sorted = false); // SD20220705
  bool WithinList(const std::vector<std::string>& list, const std::string& input, size_t& index, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const std::deque<utype>& list, const utype& input, size_t& index, bool sorted = false); // SD20220705
  bool WithinList(const std::deque<std::string>& list, const std::string& input, size_t& index, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const std::vector<utype>& list, const utype& input, bool sorted = false); // SD20220705
  bool WithinList(const std::vector<std::string>& list, const std::string& input, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const std::deque<utype>& list, const utype& input, bool sorted = false); // SD20220705
  bool WithinList(const std::deque<std::string>& list, const std::string& input, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const std::vector<utype>& list, const utype& input, std::vector<size_t>& index, bool sorted = false); // SD20220705
  bool WithinList(const std::vector<std::string>& list, const std::string& input, std::vector<size_t>& index, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const std::deque<utype>& list, const utype& input, std::vector<size_t>& index, bool sorted = false); // SD20220705
  bool WithinList(const std::deque<std::string>& list, const std::string& input, std::vector<size_t>& index, bool sorted = false); // SD20220705
  bool SubstringWithinList(const std::deque<std::string>& list, const std::string& input); // ME20220503
  bool SubstringWithinList(const std::deque<std::string>& list, const std::string& input, int& index); // ME20220503
  bool SubstringWithinList(const std::vector<std::string>& list, const std::string& input); // ME20220503
  bool SubstringWithinList(const std::vector<std::string>& list, const std::string& input, int& index); // ME20220503
  // about present substrings and taking off the value
  std::string substring2string(std::ifstream& input, const std::string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  std::string substring2string(const std::string& input, const std::string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  std::string substring2string(const std::stringstream& input, const std::string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  std::string substring2string(std::ifstream& input, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  std::string substring2string(const std::string& input, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SG20240402
  std::string substring2string(const std::stringstream& input, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype> utype substring2utype(std::ifstream& input, const std::string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype substring2utype(const std::string& input, const std::string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  template <typename utype>
  utype substring2utype(const std::stringstream& input, const std::string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  template <typename utype>
  utype substring2utype(std::ifstream& input, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype substring2utype(const std::string& input, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype substring2utype(const std::stringstream& input, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520

  bool kvpair2bool(std::ifstream& input, const std::string& keyword, const std::string& delim, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  bool kvpair2bool(const std::string& input, const std::string& keyword, const std::string& delim, bool RemoveWS = false, bool RemoveComments = true); // CO20210315
  bool kvpair2bool(const std::stringstream& input, const std::string& keyword, const std::string& delim, bool RemoveWS = false, bool RemoveComments = true); // CO20210315
  std::string kvpair2string(std::ifstream& input, const std::string& keyword, const std::string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  std::string kvpair2string(const std::string& input, const std::string& keyword, const std::string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 //SD20220520 - rewritten
  std::string kvpair2string(const std::stringstream& input, const std::string& keyword, const std::string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 //SD20220520 - rewritten
  template <typename utype>
  utype kvpair2utype(std::ifstream& input, const std::string& keyword, const std::string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype kvpair2utype(const std::string& input, const std::string& keyword, const std::string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  template <typename utype>
  utype kvpair2utype(const std::stringstream& input, const std::string& keyword, const std::string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten

  uint substring2strings(std::ifstream& input, std::vector<std::string>& vstringout, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520 - rewritten
  uint substring2strings(const std::string& input, std::vector<std::string>& vstringout, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  uint substring2strings(const std::stringstream& input, std::vector<std::string>& vstringout, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  uint substring2strings(std::ifstream& input, std::vector<std::string>& vstringout, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true);
  uint substring2strings(const std::string& input, std::vector<std::string>& vstringout, const std::string& strsub1, const std::string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SG20240402
  uint substring2strings(const std::stringstream& input,
                         std::vector<std::string>& vstringout,
                         const std::string& strsub1,
                         const std::string& strsub2,
                         const int instance = 1,
                         bool RemoveWS = false,
                         bool RemoveComments = true); // SD20220520  //CO20230502 - give instance here to speed up aurostd::substring2string(), default instance==0 //CO20230502 - trim_edges is the default behavior of substring2string()
  template <typename utype> uint substring2utypes(std::ifstream& input, std::vector<utype>& vutypeout, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // SD202205020
  template <typename utype>
  uint substring2utypes(const std::string& input, std::vector<utype>& vutypeout, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD202205020 - added utype
  template <typename utype>
  uint substring2utypes(const std::stringstream& input, std::vector<utype>& vutypeout, const std::string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD202205020 - added utype
  template <typename utype>
  uint substring2utypes(std::ifstream& input, std::vector<utype>& vutypeout, const std::string& strsub1, const std::string& strsub2, bool RemoveWS = false, bool RemoveComments = true); // SD202205020 - added utype
  template <typename utype>
  uint substring2utypes(const std::string& input, std::vector<utype>& vutypeout, const std::string& strsub1, const std::string& strsub2, bool RemoveWS = false, bool RemoveComments = true); // SD202205020
  template <typename utype>
  uint substring2utypes(const std::stringstream& input, std::vector<utype>& vutypeout, const std::string& strsub1, const std::string& strsub2, bool RemoveWS = false, bool RemoveComments = true); // SD202205020
} // namespace aurostd

// ***************************************************************************

namespace aurostd {
  std::string text2html(const std::string& str) __xprototype; // ME20200921
  std::string html2latex(const std::string& str) __xprototype;
  std::string html2txt(const std::string& str) __xprototype;
  std::string string2latex(const std::string& str) __xprototype;
  std::string latex2html(const std::string& str) __xprototype;
  std::string latex2txt(const std::string& str) __xprototype;
  std::string fixStringLatex(const std::string& input, bool double_back_slash = false, bool symmetry_string = false); // CO20190419
} // namespace aurostd

// ***************************************************************************
// SORT WORLD
// double
// ----------------------------------------------------------------------------
// sort for vectors

namespace aurostd {
  template <class utype1> void sort(std::vector<utype1>& arr);
  template <class utype1> void sort_remove_duplicates(std::vector<utype1>& arr);
  template <class utype1, class utype2> void sort(std::vector<utype1>& arr, std::vector<utype2>& brr);
  template <class utype1, class utype2> void sort(std::deque<utype1>& arr, std::deque<utype2>& brr); // CO20200915
  template <class utype1, class utype2, class utype3> void sort(std::vector<utype1>& arr, std::vector<utype2>& brr, std::vector<utype3>& crr);
  template <class utype1, class utype2, class utype3, class utype4> void sort(std::vector<utype1>& arr, std::vector<utype2>& brr, std::vector<utype3>& crr, std::vector<utype4>& drr);
} // namespace aurostd

namespace aurostd { // DOUBLE
  class _sort_double_value0 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] < v2[0]); }
  };
  class _isort_double_value0 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] > v2[0]); }
  };
  class _sort_double_value1 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[1] < v2[1]); }
  };
  class _isort_double_value1 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[1] > v2[1]); }
  };
  class _sort_double_value2 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[2] < v2[2]); }
  };
  class _isort_double_value2 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[2] > v2[2]); }
  };
  class _sort_double_value3 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[3] < v2[3]); }
  };
  class _isort_double_value3 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[3] > v2[3]); }
  };
  class _sort_double_value4 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[4] < v2[4]); }
  };
  class _isort_double_value4 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[4] > v2[4]); }
  };
  class _sort_double_value5 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[5] < v2[5]); }
  };
  class _isort_double_value5 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[5] > v2[5]); }
  };
  class _sort_double_value6 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[6] < v2[6]); }
  };
  class _isort_double_value6 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[6] > v2[6]); }
  };
  class _sort_double_value7 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[7] < v2[7]); }
  };
  class _isort_double_value7 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[7] > v2[7]); }
  };
  class _sort_double_value8 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[8] < v2[8]); }
  };
  class _isort_double_value8 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[8] > v2[8]); }
  };
  class _sort_double_value9 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[9] < v2[9]); }
  };
  class _isort_double_value9 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[9] > v2[9]); }
  };
  class _sort_double_value01 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] < v2[0] + v2[1]); }
  };
  class _isort_double_value01 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] > v2[0] + v2[1]); }
  };
  class _sort_double_value012 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] < v2[0] + v2[1] + v2[2]); }
  };
  class _isort_double_value012 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] > v2[0] + v2[1] + v2[2]); }
  };
  class _sort_double_value0123 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] < v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _isort_double_value0123 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] > v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _sort_double_value01234 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] < v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
  class _isort_double_value01234 { // sorting through reference
  public:
    bool operator()(const std::vector<double>& v1, const std::vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] > v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
} // namespace aurostd
// int
namespace aurostd { // INT
  class _sort_int_value0 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] < v2[0]); }
  };
  class _isort_int_value0 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] > v2[0]); }
  };
  class _sort_int_value1 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[1] < v2[1]); }
  };
  class _isort_int_value1 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[1] > v2[1]); }
  };
  class _sort_int_value2 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[2] < v2[2]); }
  };
  class _isort_int_value2 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[2] > v2[2]); }
  };
  class _sort_int_value3 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[3] < v2[3]); }
  };
  class _isort_int_value3 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[3] > v2[3]); }
  };
  class _sort_int_value4 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[4] < v2[4]); }
  };
  class _isort_int_value4 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[4] > v2[4]); }
  };
  class _sort_int_value5 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[5] < v2[5]); }
  };
  class _isort_int_value5 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[5] > v2[5]); }
  };
  class _sort_int_value6 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[6] < v2[6]); }
  };
  class _isort_int_value6 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[6] > v2[6]); }
  };
  class _sort_int_value7 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[7] < v2[7]); }
  };
  class _isort_int_value7 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[7] > v2[7]); }
  };
  class _sort_int_value8 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[8] < v2[8]); }
  };
  class _isort_int_value8 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[8] > v2[8]); }
  };
  class _sort_int_value9 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[9] < v2[9]); }
  };
  class _isort_int_value9 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[9] > v2[9]); }
  };
  class _sort_int_value01 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] < v2[0] + v2[1]); }
  };
  class _isort_int_value01 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] > v2[0] + v2[1]); }
  };
  class _sort_int_value012 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] < v2[0] + v2[1] + v2[2]); }
  };
  class _isort_int_value012 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] > v2[0] + v2[1] + v2[2]); }
  };
  class _sort_int_value0123 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] < v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _isort_int_value0123 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] > v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _sort_int_value01234 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] < v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
  class _isort_int_value01234 { // sorting through reference
  public:
    bool operator()(const std::vector<int>& v1, const std::vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] > v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
  // STRING
  void sort(std::vector<std::string>& arg);
  void sort(std::deque<std::string>& arg);
  void sort_remove_duplicates(std::vector<std::string>& arg);
  void sort_remove_duplicates(std::deque<std::string>& arg);
  class _sort_string_ { // sorting through reference
  public:
    bool operator()(const std::string& str1, const std::string& str2) const { return (str1 < str2); }
  };
  void rsort(std::vector<std::string>& arg);
  void rsort(std::deque<std::string>& arg);
  void rsort_remove_duplicates(std::vector<std::string>& arg);
  void rsort_remove_duplicates(std::deque<std::string>& arg);

  // _STRING_INT_
  void sort(std::vector<std::string>& varg1, std::vector<int>& varg2);
  void sort(std::deque<std::string>& varg1, std::deque<int>& varg2);
  struct _string_int_ {
    std::string arg1;
    int arg2;
  };
  class _sort_string_int_ { // sorting through reference
  public:
    bool operator()(const _string_int_& x1, const _string_int_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_DOUBLE_
  void sort(std::vector<std::string>& varg1, std::vector<double>& varg2);
  void sort(std::deque<std::string>& varg1, std::deque<double>& varg2);
  struct _string_double_ {
    std::string arg1;
    double arg2;
  };
  class _sort_string_double_ { // sorting through reference
  public:
    bool operator()(const _string_double_& x1, const _string_double_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_
  void sort(std::vector<std::string>& varg1, std::vector<std::string>& varg2);
  void sort(std::deque<std::string>& varg1, std::deque<std::string>& varg2);
  struct _string_string_ {
    std::string arg1;
    std::string arg2;
  };
  class _sort_string_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_& x1, const _string_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _DOUBLE_INT_
  void sort(std::vector<double>& varg1, std::vector<int>& varg2);
  void sort(std::deque<double>& varg1, std::deque<int>& varg2);
  struct _double_int_ {
    double arg1;
    int arg2;
  };
  class _sort_double_int_ { // sorting through reference
  public:
    bool operator()(const _double_int_& x1, const _double_int_& x2) const { return (bool) (x1.arg1 < x2.arg1); }
  };
  // _DOUBLE_DOUBLE_
  void sort(std::vector<double>& varg1, std::vector<double>& varg2);
  void sort(std::deque<double>& varg1, std::deque<double>& varg2);
  struct _double_double_ {
    double arg1;
    double arg2;
  };
  class _sort_double_double_ { // sorting through reference
  public:
    bool operator()(const _double_double_& x1, const _double_double_& x2) const { return (bool) (x1.arg1 < x2.arg1); }
  };
  // _DOUBLE_STRING_
  void sort(std::vector<double>& varg1, std::vector<std::string>& varg2);
  void sort(std::deque<double>& varg1, std::deque<std::string>& varg2);
  struct _double_string_ {
    double arg1;
    std::string arg2;
  };
  class _sort_double_string_ { // sorting through reference
  public:
    bool operator()(const _double_string_& x1, const _double_string_& x2) const { return (bool) (x1.arg1 < x2.arg1); }
  };
  // _STRING_INT_STRING
  void sort(std::vector<std::string>& varg1, std::vector<int>& varg2, std::vector<std::string>& varg3);
  void sort(std::deque<std::string>& varg1, std::deque<int>& varg2, std::deque<std::string>& varg3);
  struct _string_int_string_ {
    std::string arg1;
    int arg2;
    std::string arg3;
  };
  class _sort_string_int_string_ { // sorting through reference
  public:
    bool operator()(const _string_int_string_& x1, const _string_int_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_DOUBLE_STRING
  void sort(std::vector<std::string>& varg1, std::vector<double>& varg2, std::vector<std::string>& varg3);
  void sort(std::deque<std::string>& varg1, std::deque<double>& varg2, std::deque<std::string>& varg3);
  struct _string_double_string_ {
    std::string arg1;
    double arg2;
    std::string arg3;
  };
  class _sort_string_double_string_ { // sorting through reference
  public:
    bool operator()(const _string_double_string_& x1, const _string_double_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_STRING
  void sort(std::vector<std::string>& varg1, std::vector<std::string>& varg2, std::vector<std::string>& varg3);
  void sort(std::deque<std::string>& varg1, std::deque<std::string>& varg2, std::deque<std::string>& varg3);
  struct _string_string_string_ {
    std::string arg1;
    std::string arg2;
    std::string arg3;
  };
  class _sort_string_string_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_string_& x1, const _string_string_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_DOUBLE_STRING
  void sort(std::vector<std::string>& varg1, std::vector<std::string>& varg2, std::vector<double>& varg3, std::vector<std::string>& varg4);
  void sort(std::deque<std::string>& varg1, std::deque<std::string>& varg2, std::deque<double>& varg3, std::deque<std::string>& varg4);
  struct _string_string_double_string_ {
    std::string arg1;
    std::string arg2;
    double arg3;
    std::string arg4;
  };
  class _sort_string_string_double_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_double_string_& x1, const _string_string_double_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_DOUBLE_DOUBLE_STRING
  void sort(std::vector<std::string>& varg1, std::vector<std::string>& varg2, std::vector<double>& varg3, std::vector<double>& varg4, std::vector<std::string>& varg5);
  void sort(std::deque<std::string>& varg1, std::deque<std::string>& varg2, std::deque<double>& varg3, std::deque<double>& varg4, std::deque<std::string>& varg5);
  struct _string_string_double_double_string_ {
    std::string arg1;
    std::string arg2;
    double arg3;
    double arg4;
    std::string arg5;
  };
  class _sort_string_string_double_double_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_double_double_string_& x1, const _string_string_double_double_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
} // namespace aurostd

// ***************************************************************************
// reorder //CO20221111
namespace aurostd {
  template <class utype> void reorder(std::vector<utype>& vec, std::vector<uint>& vorder, uint mode = 1);
}
// ***************************************************************************
// some statistical stuff
namespace aurostd {
  template <class utype> utype combinations(utype n, utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination
  template <class utype> utype Cnk(utype n, utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination
} // namespace aurostd

// ***************************************************************************
namespace aurostd {
  // template <typename utype> utype sum(std::vector<utype>& a) {
  // utype result = 0;
  // for (unsigned int i=0; i<a.size();i++) result += a.at(i);
  // return result;
  // }
  std::vector<std::vector<double>> ShiftFirstColumn(const std::vector<std::vector<double>>& vva, const double& value);
  std::vector<std::vector<double>> ShrinkValuesExceptFirstColumn(const std::vector<std::vector<double>>& vva, const double& Fi);
  std::vector<std::vector<double>> NormalizeAndSum3DVector(const std::vector<std::vector<std::vector<double>>>& vvva, const std::vector<double>& vFi);
  std::vector<std::vector<double>> Sum3DVectorAndReduce2D(const std::vector<std::vector<std::vector<double>>>& vvva);
  std::vector<std::vector<double>> Sum2DVectorExceptFirstColumn(const std::vector<std::vector<double>>& vva, const std::vector<std::vector<double>>& vvb);
  std::string vector2string(const std::vector<std::vector<double>>& vva);
  template <typename utype> std::deque<utype> vector2deque(const std::vector<utype>& vin); // CO20181226
  template <typename utype> std::vector<utype> deque2vector(const std::deque<utype>& din); // CO20181226
  std::vector<std::vector<double>> ReduceVector(const std::vector<std::vector<double>>& vva, const int& n);
  double CalculateIntegrate(const std::vector<std::vector<double>>& vva, const int& n);
  double CalculateIntegrate(const std::vector<std::vector<double>>& vva, const int& n, const double& Emin, const double& Emax);
  double CalculateIntegrate(const std::vector<std::vector<double>>& vva);
  double CalculateIntegrate(const std::vector<std::vector<double>>& vva, const double& Emin, const double& Emax);
  double FindMaxIn2DvectorExcept1stColumn(const std::vector<std::vector<double>>& vva);
  double FindMaxIn2DvectorExcept1stColumn(const std::vector<std::vector<double>>& vva, const double& min, const double& max);
  double FindMaxInTDOS(const std::vector<std::vector<double>>& vva, const double& min, const double& max);
} // namespace aurostd

// ***************************************************************************

// ***************************************************************************
bool initialize_templates_never_call_this_procedure(bool flag);

////////////////////////////////////////////////////////////////////////////////
namespace aurostd {
  // joins int/string type of objects together by a delimiter
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const char delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const std::string& delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const std::string& delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const std::stringstream& delimiter);
  template <class utype> std::string joinWDelimiter(const xvector<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const char delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::string& delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::string& delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::stringstream& delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::vector<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);

  std::string joinWDelimiter(const std::vector<std::string>& sentries, const char delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const char delimiter, const char l_delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const std::string& delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const std::string& delimiter, const std::string& l_delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const std::stringstream& delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  std::string joinWDelimiter(const std::vector<std::string>& sentries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);

  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const char delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const std::string& delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const std::string& delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const std::stringstream& delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::deque<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);

  std::string joinWDelimiter(const std::deque<std::string>& sentries, const char delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const char delimiter, const char l_delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const std::string& delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const std::string& delimiter, const std::string& l_delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const std::stringstream& delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  std::string joinWDelimiter(const std::deque<std::string>& sentries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);

  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const char delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const std::string& delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const std::string& delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const std::stringstream& delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  template <class utype> std::string joinWDelimiter(const std::set<utype>& ientries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);

  std::string joinWDelimiter(const std::set<std::string>& sentries, const char delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const char delimiter, const char l_delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const std::string& delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const std::string& delimiter, const std::string& l_delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const std::string& delimiter, const std::string& m_delimiter, const std::string& l_delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const std::stringstream& delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const std::stringstream& delimiter, const std::stringstream& l_delimiter);
  std::string joinWDelimiter(const std::set<std::string>& sentries, const std::stringstream& delimiter, const std::stringstream& m_delimiter, const std::stringstream& l_delimiter);
} // namespace aurostd
////////////////////////////////////////////////////////////////////////////////

// DX20180118 - Add xcomplex to json
namespace aurostd {
  template <typename utype> std::string _xcomplex2json(xcomplex<utype>& number) __xprototype;
  std::string xcomplex2json(xcomplex<double>& number);
} // namespace aurostd

// DX20170803 - Add Matrix print out
namespace aurostd {
  std::string xmatDouble2String(const xmatrix<double>& xmat_in, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  // ME20220324
  template <typename utype> std::string xmat2String(const xmatrix<utype>& mat_in);
} // namespace aurostd

// CO20171215 - more json functionality
namespace aurostd {
  std::string wrapString(const std::string& input, const std::string& wrapper);
  std::string wrapString(const std::string& input, const std::string& wrapper_start, const std::string& wrapper_end);
} // namespace aurostd

namespace aurostd {
  std::vector<std::string> vecDouble2vecString(const std::vector<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  std::string vecDouble2String(const std::vector<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  std::vector<std::string> xvecDouble2vecString(const xvector<double>& vin,
                                                int precision = AUROSTD_DEFAULT_PRECISION,
                                                bool roff = false,
                                                double tol = (double) AUROSTD_ROUNDOFF_TOL,
                                                char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  std::string xvecDouble2String(const xvector<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  std::deque<std::string> vecDouble2vecString(const std::deque<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  std::string vecDouble2String(const std::deque<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  // std::deque<std::string> deqDouble2deqString(const std::deque<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=(double)AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM); //CO20230313 - utype clang warnings
} // namespace aurostd

namespace aurostd {
  std::vector<std::string> wrapVecEntries(const std::vector<std::string>& vin, const std::string& wrap);
  std::vector<std::string> wrapVecEntries(const std::vector<std::string>& vin, const std::string& wrap_start, const std::string& wrap_end);
  std::deque<std::string> wrapVecEntries(const std::deque<std::string>& vin, const std::string& wrap); // SC20200329 nice overload to deal with ME
  std::deque<std::string> wrapVecEntries(const std::deque<std::string>& vin, const std::string& wrap_start, const std::string& wrap_end); // SC20200329 nice overload to deal with ME
} // namespace aurostd

// base64 stuff
// CO START
namespace aurostd {
  // http://www.adp-gmbh.ch/cpp/common/base64.html for base64 encoding/decoding
  // http://stackoverflow.com/questions/535444/custom-manipulator-for-c-iostream for fancy stream manipulators
  // static inline bool isBase64(unsigned char c); //determines if char is base64
  inline bool isBase64(unsigned char c); // determines if char is base64
  std::string base64Encoder(const unsigned char* bytes_to_encode, unsigned int in_len); // encodes bytes to base64
  std::string base64Decoder(const std::string& encoded_string); // decodes base64 to bytes
  bool bin2base64(const std::string& b_file, std::string& b64String); // converts binary file to base64 string
  bool base642bin(const std::string& b64String, const std::string& b_file); // converts base64 string to binary file

  // structure that allows "cout << b64_encoder << ifstream"
  struct b64_encoder_proxy {
    explicit b64_encoder_proxy(std::ostream& os) : os(os) {}

    template <typename Rhs> friend std::ostream& operator<<(const b64_encoder_proxy& b, const Rhs& rhs) { return b.os << rhs; }

    friend std::ostream& operator<<(const b64_encoder_proxy& b, std::ifstream& file) {
      // Stop eating new lines in binary mode!!!
      file.unsetf(std::ios::skipws);

      // get its size:
      std::streampos fileSize;

      file.seekg(0, std::ios::end);
      fileSize = file.tellg();
      file.seekg(0, std::ios::beg);

      // read the data:
      std::vector<unsigned char> vec((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

      return b.os << base64Encoder(reinterpret_cast<const unsigned char*>(&vec[0]), vec.size());
    }

  private:
    std::ostream& os;
  };

  // structure that allows "ofstream << b64_decoder << string"
  // WARNING, SPITS BINARY CHARACTERS OUT IF COUT
  struct b64_decoder_proxy {
    explicit b64_decoder_proxy(std::ostream& os) : os(os) {}

    template <typename Rhs> friend std::ostream& operator<<(const b64_decoder_proxy& b, const Rhs& rhs) { return b.os << rhs; }

    friend std::ostream& operator<<(const b64_decoder_proxy& b, const std::string& rhs) { return b.os << base64Decoder(rhs); }

  private:
    std::ostream& os;
  };

  // structure that allows "cout << b64_encoder << ifstream"
  struct b64_encoder_creator {};
  extern b64_encoder_creator b64_encoder;
  b64_encoder_proxy operator<<(std::ostream& os, b64_encoder_creator);

  // structure that allows "ofstream << b64_decoder << string"
  // WARNING, SPITS BINARY CHARACTERS OUT IF COUT
  struct b64_decoder_creator {};
  extern b64_decoder_creator b64_decoder;
  b64_decoder_proxy operator<<(std::ostream& os, b64_decoder_creator);
} // namespace aurostd

namespace aurostd {
  int CountWordsinString(std::string& input); // put aurostd
  int CountWordsinString_web(std::string input); // put aurostd
} // namespace aurostd

// binary to base64 conversion
const std::string base64_chars =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";
// CO END

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
