// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_MAIN_H_
#define _AUROSTD_MAIN_H_

#include <cassert>
#include <cstdio>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>  //ME20220127 - for unit tests and multithreading
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>       //HE20210609 - for pairs in chull (C++98 changes, already included in SYMBOLICCPLUSPLUS)
#include <vector>

#include <sys/types.h>

#include "aurostd_defs.h"
#include "aurostd_xcomplex.h"
#include "aurostd_xerror.h"
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
  using std::cout;
  using std::deque;
  using std::ifstream;
  using std::ofstream;
  using std::ostream;
  using std::ostringstream;
  using std::string;
  using std::stringstream;
  using std::vector;
  void sizes() __xprototype;
  // aflow_aurostd.cpp
  // int Nint(double x);
  // int Sign(const double& x);
  // int SignNoZero(const double& x);
  template <class utype> void aswap(utype& a, utype& b);                     // some dumb algebra
  template <class utype> utype max(const vector<utype> vec);           // some dumb algebra
  template <class utype> utype max(const deque<utype> vec);            // some dumb algebra
  template <class utype> utype max(const vector<vector<utype>> mat);   // some dumb algebra
  template <class utype> utype min(const vector<utype> vec);           // some dumb algebra
  template <class utype> utype min(const deque<utype> vec);            // some dumb algebra
  template <class utype> utype min(const vector<vector<utype>> mat);   // some dumb algebra
  template <class utype> utype sum(const vector<utype> vec);           // some dumb algebra
  template <class utype> utype sum(const deque<utype> vec);            // some dumb algebra
  template <class utype> utype sum(const vector<vector<utype>> mat);   // some dumb algebra
  template <class utype> utype mean(const vector<utype> vec);          // some dumb algebra
  template <class utype> utype mean(const deque<utype> vec);           // some dumb algebra
  template <class utype> utype mean(const vector<vector<utype>> mat);  // some dumb algebra
  template <class utype> vector<utype> reset(vector<utype>& v);        // some dumb algebra
  template <class utype> deque<utype> reset(deque<utype>& v);          // some dumb algebra
  template <class utype> vector<vector<utype>> reset(vector<vector<utype>> m);  // some dumb algebra
  template <class utype> vector<utype> clear(vector<utype>& v);        // some dumb algebra
  template <class utype> deque<utype> clear(deque<utype>& v);          // some dumb algebra
  template <class utype> vector<vector<utype>> clear(vector<vector<utype>> m);  // some dumb algebra
  template <class utype> void random_shuffle(vector<utype>& v);        // some dumb algebra
  template <class utype> void random_shuffle(deque<utype>& v);         // some dumb algebra
  template <class utype> bool identical(vector<utype> v1, vector<utype> v2, utype epsilon);
  template <class utype> bool identical(deque<utype> v1, deque<utype> v2, utype epsilon);
  bool identical(vector<int> v1, vector<int> v2, int epsilon);
  bool identical(deque<int> v1, deque<int> v2, int epsilon);
  template <class utype> bool identical(const vector<vector<utype>>& m1, const vector<vector<utype>>& m2, utype epsilon);
  template <class utype> bool identical(const vector<utype>& vec, utype eps = (utype) AUROSTD_IDENTITY_TOL); // DX20210422 - checks if all values are the same
  template <class utype> bool identical(const deque<utype>& vec, utype eps = (utype) AUROSTD_IDENTITY_TOL); // DX20210422 - checks if all values are the same
  string toupper(const string& in) __xprototype;
  string tolower(const string& in) __xprototype;
  char toupper(const char& in) __xprototype;
  char tolower(const char& in) __xprototype;
  string getPWD();  // CO20191112
  int GetNumFields(const string& s);
  string GetNextVal(const string& s, int& id);
  string PaddedNumString(const int num, const int ndigits);
  int getZeroPadding(double num);  // CO20191217
  int getZeroPadding(int num);  // CO20191217
  int getZeroPadding(uint num); // CO20191217
  int getZeroPadding(long int num); // CO20191217
  int getZeroPadding(unsigned long int num);  // CO20191217
  int getZeroPadding(long long int num);  // CO20191217
  int getZeroPadding(unsigned long long int num);  // ME20190108
  template <class utype> string PaddedPRE(utype, int, string = " ");
  string PaddedPRE(string, int, string = " ");
  template <class utype> string PaddedPOST(utype, int, string = " ");
  string PaddedPOST(string, int, string = " ");
  template <class utype> string PaddedCENTER(utype, int, string = " ");
  // write progresses
  string PaddedCENTER(string, int, string = " ");
  uint ProgressBar(ostream& oss, string prelim, uint j, uint jmax, bool VERBOSE_PERCENTAGE, bool VERBOSE_ROLLER, bool VERBOSE_CURSOR);
  uint ProgressBar(ostream& oss, string prelim, uint j, uint jmax);
  uint ProgressBar(ostream& oss, string prelim, double j, bool VERBOSE_PERCENTAGE, bool VERBOSE_ROLLER, bool VERBOSE_CURSOR);
  uint ProgressBar(ostream& oss, string prelim, double j);
  // about cleaning up strings
  bool RemoveControlCodeCharactersFromString(const string& in, string& out); // DX20190516  //CO20190620
  bool RemoveControlCodeCharactersFromStringstream(stringstream& ss_in, stringstream& ss_out); // DX20190516
  bool RemoveControlCodeCharactersFromFile(const string& directory, const string& filename, bool keep_orig_file = true); // DX20190516
  bool isNullByte(char c); // DX20190131
  string removeNullBytes(string in); // DX20190131
  bool RemoveBinaryCharactersFromFile(const string& directory, const string& filename); // DX20190211 //CO20210315
  string PercentEncodeASCII(const char c) __xprototype; // DX20210706
  string CleanStringASCII(const string& s) __xprototype;
  void CleanStringASCII_InPlace(string& s) __xprototype;  // CO20190712
  string RemoveTrailingCharacter(const string& s, char c); // CO+ME20200825
  void RemoveTrailingCharacter_InPlace(string& s, char c); // CO+ME20200825
  string CGI_StringClean(const string& stringIN) __xprototype;
  string RemoveWhiteSpaces(const string& s) __xprototype;
  string RemoveWhiteSpaces(const string& s, const char toogle) __xprototype;
  string RemoveWhiteSpacesFromTheBack(const string& s) __xprototype;
  string RemoveWhiteSpacesFromTheFront(const string& s) __xprototype;
  string RemoveWhiteSpacesFromTheFrontAndBack(const string& s) __xprototype;
  string RemoveSpaces(const string& s) __xprototype;
  string RemoveSpaces(const string& s, const char toogle) __xprototype;
  string RemoveSpacesFromTheBack(const string& s) __xprototype;
  string RemoveTabs(const string& s) __xprototype;
  string RemoveTabs(const string& s, const char toogle) __xprototype;
  string RemoveTabsFromTheBack(const string& s) __xprototype;
  string RemoveComments(const string& s) __xprototype;
  vector<string> RemoveComments(const vector<string>&) __xprototype;  // ME20190614
  deque<string> RemoveComments(const deque<string>&) __xprototype;  // ME20190614
  string RemoveCharacter(const string& s, const char character) __xprototype;
  void RemoveCharacterInPlace(string& s, const char character) __xprototype;  // CO20190712
  string RemoveCharacterFromTheBack(const string& s, const char character);
  __xprototype; // DX20190708
  string RemoveCharacterFromTheFront(const string& s, const char character);
  __xprototype; // DX20190708
  string RemoveCharacterFromTheFrontAndBack(const string& s, const char character);
  __xprototype; // DX20190708
  string RemoveNumbers(const string& s) __xprototype; // CO20190712
  void RemoveNumbersInPlace(string& s) __xprototype;  // CO20190712
  string RemoveRounding(const string& s) __xprototype;
  // string RemoveCharacter(const string& s, const char character, const char toogle) __xprototype;
  string RemoveSubStringFirst(const string& str_orig, const string& str_rm) __xprototype;
  void RemoveSubStringFirstInPlace(string& str_orig, const string& str_rm) __xprototype;  // CO20190712
  string RemoveSubString(const string& str_orig, const string& str_rm) __xprototype;
  void RemoveSubStringInPlace(string& str_orig, const string& str_rm) __xprototype; // CO20190712
  double VersionString2Double(const string& version_str); // SD20220331
  vector<string> ProcessPIDs(const string& process, bool user_specific = true); // CO20210315
  vector<string> ProcessPIDs(const string& process, string& output_syscall, bool user_specific = true); // CO20210315
  vector<string> ProcessPIDs(const string& process, const string& pgid, string& output_syscall, bool user_specific = true); // SD20220329
  bool ProcessRunning(const string& process, bool user_specific = true); // CO20210315
  bool ProcessRunning(const string& process, const string& pgid, bool user_specific = true); // SD20220329
  void ProcessKill(const string& process, bool user_specific = true, uint signal = 9); // CO20210315 //SD20220627 - 9 = SIGKILL, 15 = SIGTERM
  void ProcessKill(const string& process, const string& pgid, bool user_specific = true, uint signal = 9); // SD20220329 - 9 = SIGKILL, 15 = SIGTERM
  bool ProcessRenice(const string& process, int nvalue, bool user_specific = true, const string& pgid = ""); // CO20210315
  bool ReniceAvailable(); // CO20221029

  // about printing
  void PrintANSIEscapeSequence(const aurostd::xoption& color, FILE* fstr);
  void PrintMessageStream(ostringstream& stream, bool quiet, std::ostream& oss = cout); // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet, std::ostream& oss = cout); // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet, bool osswrite, std::ostream& oss = cout);
  void PrintErrorStream(ostringstream& stream, bool quiet); // CO20200624
  void PrintErrorStream(ofstream& FileERROR, ostringstream& stream, bool quiet); // CO20200624
  void PrintErrorStream(ofstream& FileERROR, ostringstream& stream, bool quiet, bool osswrite);
  void PrintWarningStream(ostringstream& stream, bool quiet); // CO20200624
  void PrintWarningStream(ofstream& FileWARNING, ostringstream& stream, bool quiet); // CO20200624
  void PrintWarningStream(ofstream& FileWARNING, ostringstream& stream, bool quiet, bool osswrite);

  void PrintMessageStream(stringstream& stream, bool quiet, std::ostream& oss = cout);  // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet, std::ostream& oss = cout);  // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet, bool osswrite, std::ostream& oss = cout);
  void PrintErrorStream(stringstream& stream, bool quiet);  // CO20200624
  void PrintErrorStream(ofstream& FileERROR, stringstream& stream, bool quiet);  // CO20200624
  void PrintErrorStream(ofstream& FileERROR, stringstream& stream, bool quiet, bool osswrite);
  void PrintWarningStream(stringstream& stream, bool quiet);  // CO20200624
  void PrintWarningStream(ofstream& FileWARNING, stringstream& stream, bool quiet);  // CO20200624
  void PrintWarningStream(ofstream& FileWARNING, stringstream& stream, bool quiet, bool osswrite);

  // about executing

  bool IsCommandAvailable(const string& command);
  bool IsCommandAvailable(const string& command, string& position);
  bool IsCommandAvailableModify(string& command);
  bool CommandRequired(const string& command);
  bool CommandRequired(const string& command, string& position);
  bool IsExecutableAvailable(const string& executable);
  bool IsExecutableAvailable(const string& executable, string& position);
  bool ExecutableRequired(const string& executable);
  bool ExecutableRequired(const string& executable, string& position);

  bool execute(ostringstream& command);
  bool execute(stringstream& command);
  bool execute(const string& command);
  bool execute(const vector<string>& vcommand);
  bool execute(const deque<string>& dcommand);
  bool execute_thread_safe(const string& command); // HE20240220

  // about cleaning, higiene is important
  void StringstreamClean(ostringstream& stream);
  void StringstreamClean(stringstream& stream);
  int FindIfStringInStream(const string& key, std::istream& instream);

  bool GetMemoryUsagePercentage(double& usage_percentage_ram, double& usage_percentage_swap);  // CO20210601
  bool GetMemory(unsigned long long int& free_ram, unsigned long long int& total_ram, unsigned long long int& free_swap, unsigned long long int& total_swap); // CO20210315

#ifdef _stringcharstar_
  bool execute(char* command);
#endif
  // Execute and report
  std::pair<string, string> execute2OutErrPair(const string& command); // HE20240220
  string execute2string(ostringstream& command);
  string execute2string(stringstream& command);
  string execute2string(const string& command);
  vector<string> execute2string(const vector<string>& vcommand);
  deque<string> execute2string(const deque<string>& dcommand);
#ifdef _stringcharstar_
  string execute2string(char* command);
#endif
  string CleanCommand4Execute(const string& command); // CO20200624
  template <class utype> utype execute2utype(ostringstream& command);
  template <class utype> utype execute2utype(stringstream& command);
  template <class utype> utype execute2utype(string command);
  template <class utype> vector<utype> execute2utype(vector<string> vcommand);
  template <class utype> deque<utype> execute2utype(deque<string> dcommand);
#ifdef _stringcharstar_
  template <class utype> utype execute2utype(char* command);
#endif
  // about sleeping
  unsigned int Sleep(unsigned int seconds);
  // about extracting from to files
  vector<string> GrepFile(const string& filename, const string& keyword, bool RemoveWS = false, bool RemoveComments = true); // CO20210623
  // take just after
  bool ExtractJustAfterToStringstreamEXPLICIT(ifstream& FileIN, stringstream& StringstreamOUTPUT, const string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(stringstream& StringStreamIN, stringstream& StringstreamOUTPUT, const string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(const string& StringIN, stringstream& StringstreamOUTPUT, const string& Keyword_start);
  bool ExtractJustAfterToFileEXPLICIT(ifstream& FileIN, const string& FileNameOUTPUT, const string& Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(ifstream& FileIN, string& StringOUTPUT, const string& Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(const string& StringIN, string& StringOUTPUT, const string& Keyword_start);
  // about taking in istreams and stringstream and strings
  size_t stream2vectorstring(std::istream& istreamIN, vector<string>& vstringout);
  size_t stream2vectorstring(std::ifstream& ifstreamIN, vector<string>& vstringout);
  size_t stream2vectorstring(stringstream& stringstreamIN, vector<string>& vstringout);
  size_t string2vectorstring(const string& stringIN, vector<string>& vstringout, bool consecutive = false, bool trim_edges = false); // CO20170613, defaults to usual string2tokens() behavior
  vector<string> stream2vectorstring(std::istream& istreamIN);
  vector<string> stream2vectorstring(std::ifstream& ifstreamIN);
  vector<string> stream2vectorstring(stringstream& stringstreamIN);
  vector<string> string2vectorstring(const string& stringIN, bool consecutive = false, bool trim_edges = false);  // CO20170613, defaults to usual string2tokens() behavior
  string liststring2string(string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "",
                           string = "");
  uint stream2dequestring(std::istream& istreamIN, deque<string>& vstringout);
  uint stream2dequestring(std::ifstream& ifstreamIN, deque<string>& vstringout);
  uint stream2dequestring(stringstream& stringstreamIN, deque<string>& vstringout);
  uint string2dequestring(const string& stringIN, deque<string>& vstringout);
  deque<string> stream2dequestring(std::istream& istreamIN);
  deque<string> stream2dequestring(std::ifstream& ifstreamIN);
  deque<string> stream2dequestring(stringstream& stringstreamIN);
  deque<string> string2dequestring(const string& stringIN);

  // about istream/ostream
  string ostream2string(std::ostream& oss);
  uint stream2string(std::istream& istreamIN, string& vstringout);
  uint stream2string(std::ifstream& ifstreamIN, string& vstringout);
  uint stream2string(stringstream& stringstreamIN, string& vstringout);
  // about environments
  string getenv2string(const string& str);
  int getenv2int(const string& str);
  uint getenv2uint(const string& str);
  double getenv2double(const string& str);
  bool Chmod(const uint chmod, const string& path, const std::filesystem::perm_options perm_opt = std::filesystem::perm_options::replace);
  bool ChmodRecursive(const uint chmod_dir,
                      const uint chmod_file,
                      const string& directory,
                      const std::filesystem::perm_options perm_opt_dir = std::filesystem::perm_options::replace,
                      const std::filesystem::perm_options perm_opt_file = std::filesystem::perm_options::replace);

  // about getting info from strings
  uint string2tokens(const string& str, vector<string>& tokens, const string& delimiters = " ", bool consecutive = false) __xprototype;  // CO20170613, defaults to usual string2tokens() behavior
  uint string2tokens(const string& str, deque<string>& tokens, const string& delimiters = " ", bool consecutive = false) __xprototype; // CO20170613, defaults to usual string2tokens() behavior
  template <class utype>
  uint string2tokens(const string& str, std::vector<utype>& tokens, const string& delimiters = " ", bool consecutive = false) __xprototype;  // CO20170613, defaults to usual string2tokens() behavior
  template <class utype>
  uint string2tokens(const string& str, std::deque<utype>& tokens, const string& delimiters = " ", bool consecutive = false) __xprototype; // CO20170613, defaults to usual string2tokens() behavior
  uint string2tokensAdd(const string& str, vector<string>& tokens, const string& delimiters = " ") __xprototype;
  uint string2tokensAdd(const string& str, deque<string>& tokens, const string& delimiters = " ") __xprototype;
  template <class utype> uint string2tokensAdd(const string& str, std::vector<utype>& tokens, const string& delimiters = " ") __xprototype;
  template <class utype> uint string2tokensAdd(const string& str, std::deque<utype>& tokens, const string& delimiters = " ") __xprototype;
  uint string2tokensByDelimiter(const string& str, vector<string>& tokens, const string& delimiter); // SD20220504
  uint string2tokensByDelimiter(const string& str, deque<string>& tokens, const string& delimiter); // SD20220504

  //[CO20210315 - not defined]template<typename typeTo, typename typeFrom> typeTo NumberStreamConvert(const typeFrom& from);  //CO20210315 - cleaned up

  template <typename utype> vector<utype> vectorstring2vectorutype(const vector<string>& from); // SD20220520
  vector<double> vectorstring2vectordouble(const vector<string>& from); // CO20210315 - cleaned up
  string string2string(const string& from) __xprototype;
  template <typename utype> utype string2utype(const string& from, const uint base = 10);  // CO20210315 - cleaned up //HE20220324 add base option
  vector<int> vectorstring2vectorint(const vector<string>& from); // CO20210315 - cleaned up
  vector<uint> vectorstring2vectoruint(const vector<string>& from); // CO20210315 - cleaned up

  vector<float> vectorstring2vectorfloat(const vector<string>& from);  // CO20210315 - cleaned up
  string vectorstring2string(const vector<string>& vstrings);
  string vectorstring2string(const deque<string>& vstrings);

  template <typename utype>
  string utype2string(const utype& from, int precision = AUROSTD_DEFAULT_PRECISION, char FORMAT = DEFAULT_STREAM) __xprototype; // DX20201028 - this declaration was missing //DX20210128 - add defaults
  string utype2string(double from, bool roff);
  string utype2string(double from, int precision, bool roff);
  string utype2string(double from, bool roff, double tol);
  string utype2string(double from, int precision, bool roff, double tol);
  string utype2string(double from, bool roff, char FORMAT);
  string utype2string(double from, int precision, char FORMAT, bool roff = false); // CO20200624
  string utype2string(double from, int precision, bool roff, char FORMAT);
  string utype2string(double from, bool roff, double tol, char FORMAT);
  string utype2string(double from, int precision, bool roff, double tol, char FORMAT);
  string bool2string(bool from);

  template <class utype> deque<utype> utypes2deque(utype u1) __xprototype;
  template <class utype> deque<utype> utypes2deque(utype u1, utype u2) __xprototype;
  template <class utype> deque<utype> utypes2deque(utype u1, utype u2, utype u3) __xprototype;
  template <class utype> deque<utype> utypes2deque(utype u1, utype u2, utype u3, utype u4) __xprototype;

  void StringCommasColumsVectorInt(string vstring, vector<int>& vint) __xprototype;
  void StringCommasColumsVectorUnsignedInt(string vstring, vector<uint>& vuint) __xprototype;
  void StringCommasColumsVectorFloat(string vstring, vector<float>& vfloat) __xprototype;
  void StringCommasColumsVectorDouble(string vstring, vector<double>& vdouble) __xprototype;
  size_t GetNLinesString(const string& str) __xprototype;
  size_t GetNLinesString(const stringstream& strstream) __xprototype;
  size_t GetNLinesFile(const string& file_name) __xprototype;
  string GetLineString(const string& strstream, int line);
  string GetLineString(const stringstream& strstream, int line);
  // substitute strings in strings and stringstreams
  bool StringsAlphabetic(const string& A, const string& B, bool allow_identical = true); // CO20180801
  bool StringsAlphabetic(const vector<string>& input, bool allow_identical = true);  // CO20180801
  bool StringsAlphabetic(const deque<string>& input, bool allow_identical = true);  // CO20180801
  void StringSubstInPlace(string& work_string, const string& old_string, const string& new_string); // HE20240908
  void StringSubstInPlace(string& work_string, char old_char, char new_char); // HE20240908
  string StringSubst(const string& work_string, const string& old_string, const string& new_string); // HE20220321
  string StringSubst(const string& work_string, char old_char, char new_char); // HE20220321
  void StringStreamSubst(stringstream& strstring, const string& strfind, const string& strreplace);  // ME20190128 - fixed type declaration
  // about present substrings
  bool substring2bool(const string& strstream, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315 - cleaned up
  bool substring2bool(const vector<string>& vstrstream, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up
  bool substring2bool(const deque<string>& vstrstream, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315 - cleaned up
  bool substring2bool(const stringstream& strstream, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315 - cleaned up
  bool substringlist2bool(const string& strin, const vector<string>& substrings, bool match_all = true); // ME20220505
  bool substringlist2bool(const string& strin, const deque<string>& substrings, bool match_all = true);  // ME20220505
  bool substring_present_file(const string& FileName, const std::boyer_moore_searcher<std::string::const_iterator>& boyer_moore, const size_t strsub_size, size_t& offset);
  bool substring_present_file(const string& FileName, const string& strsub, size_t& offset);
  bool substring_present_file(const string& FileName, const string& strsub);
  vector<bool> substrings_present_file(const string& FileName, const vector<std::boyer_moore_searcher<std::string::const_iterator>>& vboyer_moore, const size_t strsub_size, size_t& offset);
  vector<bool> substrings_present_file(const string& FileName, const vector<string>& vstrsub, size_t& offset);
  vector<bool> substrings_present_file(const string& FileName, const vector<string>& vstrsub);
  void substrings_map_present_file(const string& FileName, std::map<string, std::pair<std::boyer_moore_searcher<std::string::const_iterator>, bool>>& map_boyer_moore, const size_t strsub_size, size_t& offset);
  void substrings_map_present_file(const string& FileName, std::unordered_map<string, bool>& map_strsub, size_t& offset);
  void substrings_map_present_file(const string& FileName, std::unordered_map<string, bool>& map_strsub);
  template <class utype> bool WithinList(const vector<utype>& list, const utype& input, size_t& index, bool sorted = false); // SD20220705
  bool WithinList(const vector<string>& list, const string& input, size_t& index, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const deque<utype>& list, const utype& input, size_t& index, bool sorted = false); // SD20220705
  bool WithinList(const deque<string>& list, const string& input, size_t& index, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const vector<utype>& list, const utype& input, bool sorted = false); // SD20220705
  bool WithinList(const vector<string>& list, const string& input, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const deque<utype>& list, const utype& input, bool sorted = false); // SD20220705
  bool WithinList(const deque<string>& list, const string& input, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const vector<utype>& list, const utype& input, vector<size_t>& index, bool sorted = false); // SD20220705
  bool WithinList(const vector<string>& list, const string& input, vector<size_t>& index, bool sorted = false); // SD20220705
  template <class utype> bool WithinList(const deque<utype>& list, const utype& input, vector<size_t>& index, bool sorted = false); // SD20220705
  bool WithinList(const deque<string>& list, const string& input, vector<size_t>& index, bool sorted = false); // SD20220705
  bool SubstringWithinList(const deque<string>& list, const string& input);  // ME20220503
  bool SubstringWithinList(const deque<string>& list, const string& input, int& index);  // ME20220503
  bool SubstringWithinList(const vector<string>& list, const string& input);  // ME20220503
  bool SubstringWithinList(const vector<string>& list, const string& input, int& index);  // ME20220503
  // about present substrings and taking off the value
  string substring2string(ifstream& input, const string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  string substring2string(const string& input, const string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315 - cleaned up //SD20220520 - rewritten
  string substring2string(const stringstream& input, const string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315 - cleaned up //SD20220520 - rewritten
  string substring2string(ifstream& input, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  string substring2string(const string& input, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SG20240402
  string substring2string(const stringstream& input, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype> utype substring2utype(ifstream& input, const string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype substring2utype(const string& input, const string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  template <typename utype>
  utype substring2utype(const stringstream& input, const string& strsub1, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  template <typename utype> utype substring2utype(ifstream& input, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype> utype substring2utype(const string& input, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype substring2utype(const stringstream& input, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520

  bool kvpair2bool(ifstream& input, const string& keyword, const string& delim, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  bool kvpair2bool(const string& input, const string& keyword, const string& delim, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315
  bool kvpair2bool(const stringstream& input, const string& keyword, const string& delim, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315
  string kvpair2string(ifstream& input, const string& keyword, const string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  string kvpair2string(const string& input, const string& keyword, const string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true);  // CO20210315 //SD20220520 - rewritten
  string kvpair2string(const stringstream& input, const string& keyword, const string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 //SD20220520 - rewritten
  template <typename utype> utype kvpair2utype(ifstream& input, const string& keyword, const string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520
  template <typename utype>
  utype kvpair2utype(const string& input, const string& keyword, const string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  template <typename utype>
  utype kvpair2utype(const stringstream& input, const string& keyword, const string& delim, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten

  uint substring2strings(ifstream& input, vector<string>& vstringout, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // SD20220520 - rewritten
  uint substring2strings(const string& input, vector<string>& vstringout, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  uint substring2strings(const stringstream& input, vector<string>& vstringout, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD20220520 - rewritten
  uint substring2strings(ifstream& input, vector<string>& vstringout, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true);
  uint substring2strings(const string& input, vector<string>& vstringout, const string& strsub1, const string& strsub2, const int instance = 1, bool RemoveWS = false, bool RemoveComments = true); // SG20240402
  uint substring2strings(const stringstream& input,
                         vector<string>& vstringout,
                         const string& strsub1,
                         const string& strsub2,
                         const int instance = 1,
                         bool RemoveWS = false,
                         bool RemoveComments = true); // SD20220520  //CO20230502 - give instance here to speed up aurostd::substring2string(), default instance==0 //CO20230502 - trim_edges is the default behavior of substring2string()
  template <typename utype> uint substring2utypes(ifstream& input, vector<utype>& vutypeout, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // SD202205020
  template <typename utype>
  uint substring2utypes(const string& input, vector<utype>& vutypeout, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD202205020 - added utype
  template <typename utype>
  uint substring2utypes(const stringstream& input, vector<utype>& vutypeout, const string& strsub1, bool RemoveWS = false, bool RemoveComments = true); // CO20210315 - cleaned up //SD202205020 - added utype
  template <typename utype>
  uint substring2utypes(ifstream& input, vector<utype>& vutypeout, const string& strsub1, const string& strsub2, bool RemoveWS = false, bool RemoveComments = true); // SD202205020 - added utype
  template <typename utype> uint substring2utypes(const string& input, vector<utype>& vutypeout, const string& strsub1, const string& strsub2, bool RemoveWS = false, bool RemoveComments = true); // SD202205020
  template <typename utype>
  uint substring2utypes(const stringstream& input, vector<utype>& vutypeout, const string& strsub1, const string& strsub2, bool RemoveWS = false, bool RemoveComments = true); // SD202205020
} // namespace aurostd

// ***************************************************************************

namespace aurostd {
  string text2html(const string& str) __xprototype; // ME20200921
  string html2latex(const string& str) __xprototype;
  string html2txt(const string& str) __xprototype;
  string string2latex(const string& str) __xprototype;
  string latex2html(const string& str) __xprototype;
  string latex2txt(const string& str) __xprototype;
  string fixStringLatex(const string& input, bool double_back_slash = false, bool symmetry_string = false); // CO20190419
} // namespace aurostd

// ***************************************************************************
// SORT WORLD
// double
// ----------------------------------------------------------------------------
// sort for vectors

namespace aurostd {
  template <class utype1> void sort(vector<utype1>& arr);
  template <class utype1> void sort_remove_duplicates(vector<utype1>& arr);
  template <class utype1, class utype2> void sort(vector<utype1>& arr, vector<utype2>& brr);
  template <class utype1, class utype2> void sort(deque<utype1>& arr, deque<utype2>& brr); // CO20200915
  template <class utype1, class utype2, class utype3> void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr);
  template <class utype1, class utype2, class utype3, class utype4> void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr, vector<utype4>& drr);
} // namespace aurostd

namespace aurostd { // DOUBLE
  class _sort_double_value0 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] < v2[0]); }
  };
  class _isort_double_value0 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] > v2[0]); }
  };
  class _sort_double_value1 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[1] < v2[1]); }
  };
  class _isort_double_value1 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[1] > v2[1]); }
  };
  class _sort_double_value2 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[2] < v2[2]); }
  };
  class _isort_double_value2 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[2] > v2[2]); }
  };
  class _sort_double_value3 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[3] < v2[3]); }
  };
  class _isort_double_value3 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[3] > v2[3]); }
  };
  class _sort_double_value4 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[4] < v2[4]); }
  };
  class _isort_double_value4 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[4] > v2[4]); }
  };
  class _sort_double_value5 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[5] < v2[5]); }
  };
  class _isort_double_value5 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[5] > v2[5]); }
  };
  class _sort_double_value6 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[6] < v2[6]); }
  };
  class _isort_double_value6 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[6] > v2[6]); }
  };
  class _sort_double_value7 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[7] < v2[7]); }
  };
  class _isort_double_value7 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[7] > v2[7]); }
  };
  class _sort_double_value8 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[8] < v2[8]); }
  };
  class _isort_double_value8 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[8] > v2[8]); }
  };
  class _sort_double_value9 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[9] < v2[9]); }
  };
  class _isort_double_value9 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[9] > v2[9]); }
  };
  class _sort_double_value01 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] < v2[0] + v2[1]); }
  };
  class _isort_double_value01 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] > v2[0] + v2[1]); }
  };
  class _sort_double_value012 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] < v2[0] + v2[1] + v2[2]); }
  };
  class _isort_double_value012 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] > v2[0] + v2[1] + v2[2]); }
  };
  class _sort_double_value0123 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] < v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _isort_double_value0123 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] > v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _sort_double_value01234 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] < v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
  class _isort_double_value01234 { // sorting through reference
  public:
    bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] > v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
} // namespace aurostd
// int
namespace aurostd { // INT
  class _sort_int_value0 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] < v2[0]); }
  };
  class _isort_int_value0 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] > v2[0]); }
  };
  class _sort_int_value1 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[1] < v2[1]); }
  };
  class _isort_int_value1 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[1] > v2[1]); }
  };
  class _sort_int_value2 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[2] < v2[2]); }
  };
  class _isort_int_value2 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[2] > v2[2]); }
  };
  class _sort_int_value3 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[3] < v2[3]); }
  };
  class _isort_int_value3 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[3] > v2[3]); }
  };
  class _sort_int_value4 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[4] < v2[4]); }
  };
  class _isort_int_value4 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[4] > v2[4]); }
  };
  class _sort_int_value5 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[5] < v2[5]); }
  };
  class _isort_int_value5 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[5] > v2[5]); }
  };
  class _sort_int_value6 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[6] < v2[6]); }
  };
  class _isort_int_value6 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[6] > v2[6]); }
  };
  class _sort_int_value7 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[7] < v2[7]); }
  };
  class _isort_int_value7 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[7] > v2[7]); }
  };
  class _sort_int_value8 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[8] < v2[8]); }
  };
  class _isort_int_value8 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[8] > v2[8]); }
  };
  class _sort_int_value9 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[9] < v2[9]); }
  };
  class _isort_int_value9 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[9] > v2[9]); }
  };
  class _sort_int_value01 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] < v2[0] + v2[1]); }
  };
  class _isort_int_value01 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] > v2[0] + v2[1]); }
  };
  class _sort_int_value012 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] < v2[0] + v2[1] + v2[2]); }
  };
  class _isort_int_value012 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] > v2[0] + v2[1] + v2[2]); }
  };
  class _sort_int_value0123 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] < v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _isort_int_value0123 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] > v2[0] + v2[1] + v2[2] + v2[3]); }
  };
  class _sort_int_value01234 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] < v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
  class _isort_int_value01234 { // sorting through reference
  public:
    bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0] + v1[1] + v1[2] + v1[3] + v1[4] > v2[0] + v2[1] + v2[2] + v2[3] + v2[4]); }
  };
  // STRING
  void sort(vector<string>& arg);
  void sort(deque<string>& arg);
  void sort_remove_duplicates(vector<string>& arg);
  void sort_remove_duplicates(deque<string>& arg);
  class _sort_string_ { // sorting through reference
  public:
    bool operator()(const string& str1, const string& str2) const { return (str1 < str2); }
  };
  void rsort(vector<string>& arg);
  void rsort(deque<string>& arg);
  void rsort_remove_duplicates(vector<string>& arg);
  void rsort_remove_duplicates(deque<string>& arg);

  // _STRING_INT_
  void sort(vector<string>& varg1, vector<int>& varg2);
  void sort(deque<string>& varg1, deque<int>& varg2);
  struct _string_int_ {
    string arg1;
    int arg2;
  };
  class _sort_string_int_ { // sorting through reference
  public:
    bool operator()(const _string_int_& x1, const _string_int_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_DOUBLE_
  void sort(vector<string>& varg1, vector<double>& varg2);
  void sort(deque<string>& varg1, deque<double>& varg2);
  struct _string_double_ {
    string arg1;
    double arg2;
  };
  class _sort_string_double_ { // sorting through reference
  public:
    bool operator()(const _string_double_& x1, const _string_double_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_
  void sort(vector<string>& varg1, vector<string>& varg2);
  void sort(deque<string>& varg1, deque<string>& varg2);
  struct _string_string_ {
    string arg1;
    string arg2;
  };
  class _sort_string_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_& x1, const _string_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _DOUBLE_INT_
  void sort(vector<double>& varg1, vector<int>& varg2);
  void sort(deque<double>& varg1, deque<int>& varg2);
  struct _double_int_ {
    double arg1;
    int arg2;
  };
  class _sort_double_int_ { // sorting through reference
  public:
    bool operator()(const _double_int_& x1, const _double_int_& x2) const { return (bool) (x1.arg1 < x2.arg1); }
  };
  // _DOUBLE_DOUBLE_
  void sort(vector<double>& varg1, vector<double>& varg2);
  void sort(deque<double>& varg1, deque<double>& varg2);
  struct _double_double_ {
    double arg1;
    double arg2;
  };
  class _sort_double_double_ { // sorting through reference
  public:
    bool operator()(const _double_double_& x1, const _double_double_& x2) const { return (bool) (x1.arg1 < x2.arg1); }
  };
  // _DOUBLE_STRING_
  void sort(vector<double>& varg1, vector<string>& varg2);
  void sort(deque<double>& varg1, deque<string>& varg2);
  struct _double_string_ {
    double arg1;
    string arg2;
  };
  class _sort_double_string_ { // sorting through reference
  public:
    bool operator()(const _double_string_& x1, const _double_string_& x2) const { return (bool) (x1.arg1 < x2.arg1); }
  };
  // _STRING_INT_STRING
  void sort(vector<string>& varg1, vector<int>& varg2, vector<string>& varg3);
  void sort(deque<string>& varg1, deque<int>& varg2, deque<string>& varg3);
  struct _string_int_string_ {
    string arg1;
    int arg2;
    string arg3;
  };
  class _sort_string_int_string_ { // sorting through reference
  public:
    bool operator()(const _string_int_string_& x1, const _string_int_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_DOUBLE_STRING
  void sort(vector<string>& varg1, vector<double>& varg2, vector<string>& varg3);
  void sort(deque<string>& varg1, deque<double>& varg2, deque<string>& varg3);
  struct _string_double_string_ {
    string arg1;
    double arg2;
    string arg3;
  };
  class _sort_string_double_string_ { // sorting through reference
  public:
    bool operator()(const _string_double_string_& x1, const _string_double_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_STRING
  void sort(vector<string>& varg1, vector<string>& varg2, vector<string>& varg3);
  void sort(deque<string>& varg1, deque<string>& varg2, deque<string>& varg3);
  struct _string_string_string_ {
    string arg1;
    string arg2;
    string arg3;
  };
  class _sort_string_string_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_string_& x1, const _string_string_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_DOUBLE_STRING
  void sort(vector<string>& varg1, vector<string>& varg2, vector<double>& varg3, vector<string>& varg4);
  void sort(deque<string>& varg1, deque<string>& varg2, deque<double>& varg3, deque<string>& varg4);
  struct _string_string_double_string_ {
    string arg1;
    string arg2;
    double arg3;
    string arg4;
  };
  class _sort_string_string_double_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_double_string_& x1, const _string_string_double_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
  // _STRING_STRING_DOUBLE_DOUBLE_STRING
  void sort(vector<string>& varg1, vector<string>& varg2, vector<double>& varg3, vector<double>& varg4, vector<string>& varg5);
  void sort(deque<string>& varg1, deque<string>& varg2, deque<double>& varg3, deque<double>& varg4, deque<string>& varg5);
  struct _string_string_double_double_string_ {
    string arg1;
    string arg2;
    double arg3;
    double arg4;
    string arg5;
  };
  class _sort_string_string_double_double_string_ { // sorting through reference
  public:
    bool operator()(const _string_string_double_double_string_& x1, const _string_string_double_double_string_& x2) const { return (x1.arg1 < x2.arg1); }
  };
} // namespace aurostd

// ***************************************************************************
// reorder //CO20221111
namespace aurostd {
  template <class utype> void reorder(vector<utype>& vec, vector<uint>& vorder, uint mode = 1);
}
// ***************************************************************************
// some statistical stuff
namespace aurostd {
  template <class utype> utype combinations(utype n, utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination
  template <class utype> utype Cnk(utype n, utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination
} // namespace aurostd

// ***************************************************************************
namespace aurostd {
  // template <typename utype> utype sum(vector<utype>& a) {
  // utype result = 0;
  // for (unsigned int i=0; i<a.size();i++) result += a.at(i);
  // return result;
  // }
  vector<vector<double>> ShiftFirstColumn(const vector<vector<double>>& vva, const double& value);
  vector<vector<double>> ShrinkValuesExceptFirstColumn(const vector<vector<double>>& vva, const double& Fi);
  vector<vector<double>> NormalizeAndSum3DVector(const vector<vector<vector<double>>>& vvva, const vector<double>& vFi);
  vector<vector<double>> Sum3DVectorAndReduce2D(const vector<vector<vector<double>>>& vvva);
  vector<vector<double>> Sum2DVectorExceptFirstColumn(const vector<vector<double>>& vva, const vector<vector<double>>& vvb);
  string vector2string(const vector<vector<double>>& vva);
  template <typename utype> deque<utype> vector2deque(const vector<utype>& vin); // CO20181226
  template <typename utype> vector<utype> deque2vector(const deque<utype>& din); // CO20181226
  vector<vector<double>> ReduceVector(const vector<vector<double>>& vva, const int& n);
  double CalculateIntegrate(const vector<vector<double>>& vva, const int& n);
  double CalculateIntegrate(const vector<vector<double>>& vva, const int& n, const double& Emin, const double& Emax);
  double CalculateIntegrate(const vector<vector<double>>& vva);
  double CalculateIntegrate(const vector<vector<double>>& vva, const double& Emin, const double& Emax);
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>>& vva);
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>>& vva, const double& min, const double& max);
  double FindMaxInTDOS(const vector<vector<double>>& vva, const double& min, const double& max);
} // namespace aurostd

// ***************************************************************************

// ***************************************************************************
bool initialize_templates_never_call_this_procedure(bool flag);

////////////////////////////////////////////////////////////////////////////////
namespace aurostd {
  // joins int/string type of objects together by a delimiter
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const char delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const string& delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const string& delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const stringstream& delimiter);
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const char delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const string& delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const string& delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter);
  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);

  string joinWDelimiter(const vector<string>& sentries, const char delimiter);
  string joinWDelimiter(const vector<string>& sentries, const char delimiter, const char l_delimiter);
  string joinWDelimiter(const vector<string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter);
  string joinWDelimiter(const vector<string>& sentries, const string& delimiter);
  string joinWDelimiter(const vector<string>& sentries, const string& delimiter, const string& l_delimiter);
  string joinWDelimiter(const vector<string>& sentries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  string joinWDelimiter(const vector<string>& sentries, const stringstream& delimiter);
  string joinWDelimiter(const vector<string>& sentries, const stringstream& delimiter, const stringstream& l_delimiter);
  string joinWDelimiter(const vector<string>& sentries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const char delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const string& delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const string& delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const stringstream& delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter);
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);

  string joinWDelimiter(const deque<string>& sentries, const char delimiter);
  string joinWDelimiter(const deque<string>& sentries, const char delimiter, const char l_delimiter);
  string joinWDelimiter(const deque<string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter);
  string joinWDelimiter(const deque<string>& sentries, const string& delimiter);
  string joinWDelimiter(const deque<string>& sentries, const string& delimiter, const string& l_delimiter);
  string joinWDelimiter(const deque<string>& sentries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  string joinWDelimiter(const deque<string>& sentries, const stringstream& delimiter);
  string joinWDelimiter(const deque<string>& sentries, const stringstream& delimiter, const stringstream& l_delimiter);
  string joinWDelimiter(const deque<string>& sentries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const char delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const char delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const string& delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const string& delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const stringstream& delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter);
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);

  string joinWDelimiter(const std::set<string>& sentries, const char delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const char delimiter, const char l_delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const string& delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const string& delimiter, const string& l_delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const string& delimiter, const string& m_delimiter, const string& l_delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const stringstream& delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const stringstream& delimiter, const stringstream& l_delimiter);
  string joinWDelimiter(const std::set<string>& sentries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter);
} // namespace aurostd
////////////////////////////////////////////////////////////////////////////////

// DX20180118 - Add xcomplex to json
namespace aurostd {
  template <typename utype> string _xcomplex2json(xcomplex<utype>& number) __xprototype;
  string xcomplex2json(xcomplex<double>& number);
} // namespace aurostd

// DX20170803 - Add Matrix print out
namespace aurostd {
  string xmatDouble2String(const xmatrix<double>& xmat_in, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  // ME20220324
  template <typename utype> string xmat2String(const xmatrix<utype>& mat_in);
} // namespace aurostd

// CO20171215 - more json functionality
namespace aurostd {
  string wrapString(const string& input, const string& wrapper);
  string wrapString(const string& input, const string& wrapper_start, const string& wrapper_end);
} // namespace aurostd

namespace aurostd {
  vector<string> vecDouble2vecString(const vector<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  string vecDouble2String(const vector<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  vector<string> xvecDouble2vecString(const xvector<double>& vin,
                                      int precision = AUROSTD_DEFAULT_PRECISION,
                                      bool roff = false,
                                      double tol = (double) AUROSTD_ROUNDOFF_TOL,
                                      char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  string xvecDouble2String(const xvector<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  deque<string> vecDouble2vecString(const deque<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  string vecDouble2String(const deque<double>& vin, int precision = AUROSTD_DEFAULT_PRECISION, bool roff = false, double tol = (double) AUROSTD_ROUNDOFF_TOL, char FORMAT = DEFAULT_STREAM); // CO20230313 - utype clang warnings
  // deque<string> deqDouble2deqString(const deque<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=(double)AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM); //CO20230313 - utype clang warnings
} // namespace aurostd

namespace aurostd {
  vector<string> wrapVecEntries(const vector<string>& vin, const string& wrap);
  vector<string> wrapVecEntries(const vector<string>& vin, const string& wrap_start, const string& wrap_end);
  deque<string> wrapVecEntries(const deque<string>& vin, const string& wrap); // SC20200329 nice overload to deal with ME
  deque<string> wrapVecEntries(const deque<string>& vin, const string& wrap_start, const string& wrap_end); // SC20200329 nice overload to deal with ME
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

    friend std::ostream& operator<<(const b64_encoder_proxy& b, ifstream& file) {
      // Stop eating new lines in binary mode!!!
      file.unsetf(std::ios::skipws);

      // get its size:
      std::streampos fileSize;

      file.seekg(0, std::ios::end);
      fileSize = file.tellg();
      file.seekg(0, std::ios::beg);

      // read the data:
      vector<unsigned char> vec((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

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
