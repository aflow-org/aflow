// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <spawn.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include "aurostd.h"
#include "aurostd_automatic_template.h"
#include "aurostd_defs.h"
#include "aurostd_time.h"
#include "aurostd_xcomplex.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"
#include "aurostd_xhttp.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xoption.h"
#include "aurostd_xscalar.h"
#include "aurostd_xvector.h"

#include "aflow_xhost.h" // todo required for XPID use and XHOST.DEBUG use
#include "aflowlib/aflowlib_web_interface.h" // NOLINT // todo required for a template declaration, should be in header only

#define _CIN_LINE_BUFFER_LENGTH_ 16384
#ifndef CHMOD_BIN
#define CHMOD_BIN string("chmod")
#endif

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xmatrix;
using aurostd::xvector;

extern char** environ; // HE20240220 for posix_spawn()

#define COMMENT_NEGLECT_1 string("#")
// #define COMMENT_NEGLECT_2 string("// ")
#define COMMENT_NEGLECT_2 string("//")
#define COMMENT_NEGLECT_3 string("!")

#ifdef AFLOW_MULTITHREADS_ENABLE  // CO+HE20221116
#include <mutex>
// Global mutex that prevents two xThread instances from executing a system call.
// system calls are not generally thread-safe: https://stackoverflow.com/questions/12025640/how-can-i-know-whether-a-linux-syscall-is-thread-safe
static std::mutex xthread_execute;
#endif

// ***************************************************************************
// get threadID
namespace aurostd {
  unsigned long long int getTID() { // CO20200502 - threadID
    // for mac these numbers can be QUITE large, so better to be safe and return unsigned long long int
    // see here: http://elliotth.blogspot.com/2012/04/gettid-on-mac-os.html
    // also for macs: pid!=tid
#ifdef _MACOSX_
#if MAC_OS_X_VERSION_MAX_ALLOWED >= MAC_OS_X_VERSION_10_12
    uint64_t tid64;
    pthread_threadid_np(NULL, &tid64);
    pid_t tid = (pid_t) tid64;
    return (unsigned long long int) tid;
#else
    //////////////////////////////////////////////////////////
#ifdef __GLIBC__
#include <sys/syscall.h>  //CO20200502 - need for gettid()
    pid_t tid = syscall(__NR_gettid);
    return (unsigned long long int) tid;
#else // ONLY if _MACOSX_ AND not __GLIBC__
    return (unsigned long long int) getpid();
#endif
    //////////////////////////////////////////////////////////
#endif  // END _MACOSX_
#else // if NOT _MACOSX_
    //////////////////////////////////////////////////////////
#ifdef __GLIBC__
    return (unsigned long long int) gettid();
#else // for example CYGWIN
    return (unsigned long long int) getpid();
#endif
    //////////////////////////////////////////////////////////
#endif
  }
} // namespace aurostd

// ***************************************************************************
// Function extra operator << for vector
template <class utype>                            // operator <<  vector<>
std::ostream& operator<<(std::ostream& buf, const std::vector<utype>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    buf << x[i] << " ";
  }
  return buf;
}
// ***************************************************************************
// Function extra operator << for deque
template <class utype>                            // operator <<  deque<>
std::ostream& operator<<(std::ostream& buf, const std::deque<utype>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    buf << x[i] << " ";
  }
  return buf;
}
// ***************************************************************************
std::ostream& operator<<(std::ostream& b, const vector<uint>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<uint>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const vector<char>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<char>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const vector<int>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<int>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const vector<long>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<long>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const vector<double>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<double>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const vector<long double>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<long double>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const vector<string>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}
std::ostream& operator<<(std::ostream& b, const deque<string>& x) {
  for (size_t i = 0; i < x.size(); i++) {
    b << x[i] << " ";
  }
  return b;
}

namespace aurostd {
  // ***************************************************************************
  // Function aswap
  // ***************************************************************************
  // namespace aurostd
  template <class utype> void aswap(utype& a, utype& b) {
    utype temp = a;
    a = b;
    b = temp;
  }

  // ***************************************************************************
  // max functions
  // ***************************************************************************
  // SC
  template <class utype> utype max(const vector<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = vec.at(0);
    for (size_t i = 0; i < vec.size(); i++) {
      if (vec[i] >= out) {
        out = vec[i];
      }
    }
    return out;
  }

  template <class utype> utype max(const deque<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = vec.at(0);
    for (size_t i = 0; i < vec.size(); i++) {
      if (vec[i] >= out) {
        out = vec[i];
      }
    }
    return out;
  }

  template <class utype> utype max(const vector<vector<utype>> mat) {
    if (mat.size() == 0) {
      return (utype) 0;
    }
    if (mat.at(0).size() == 0) {
      return (utype) 0;
    }
    utype out = mat.at(0).at(0);
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        if (mat[i][j] >= out) {
          out = mat[i][j];
        }
      }
    }
    return out;
  }

#define AST_TEMPLATE(utype) template utype max(const vector<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype max(const deque<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype max(const vector<vector<utype>>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ***************************************************************************
  // min functions
  // ***************************************************************************
  // SC
  template <class utype> utype min(const vector<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = vec.at(0);
    for (size_t i = 0; i < vec.size(); i++) {
      if (vec[i] <= out) {
        out = vec[i];
      }
    }
    return out;
  }

  template <class utype> utype min(const deque<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = vec.at(0);
    for (size_t i = 0; i < vec.size(); i++) {
      if (vec[i] <= out) {
        out = vec[i];
      }
    }
    return out;
  }

  template <class utype> utype min(const vector<vector<utype>> mat) {
    if (mat.size() == 0) {
      return (utype) 0;
    }
    if (mat.at(0).size() == 0) {
      return (utype) 0;
    }
    utype out = mat.at(0).at(0);
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        if (mat[i][j] <= out) {
          out = mat[i][j];
        }
      }
    }
    return out;
  }

#define AST_TEMPLATE(utype) template utype min(const vector<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype min(const deque<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype min(const vector<vector<utype>>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ***************************************************************************
  // sum functions
  // ***************************************************************************
  template <class utype> utype sum(const vector<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = 0;
    for (size_t i = 0; i < vec.size(); i++) {
      out += vec[i];
    }
    return out;
  }

  template <class utype> utype sum(const deque<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = 0;
    for (size_t i = 0; i < vec.size(); i++) {
      out += vec[i];
    }
    return out;
  }

  template <class utype> utype sum(const vector<vector<utype>> mat) {
    if (mat.size() == 0) {
      return (utype) 0;
    }
    if (mat.at(0).size() == 0) {
      return (utype) 0;
    }
    utype out = 0;
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        out += mat[i][j];
      }
    }
    return out;
  }

#define AST_TEMPLATE(utype) template utype sum(const vector<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype sum(const deque<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype sum(const vector<vector<utype>>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ***************************************************************************
  // mean functions
  // ***************************************************************************
  template <class utype> utype mean(const vector<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = 0;
    for (size_t i = 0; i < vec.size(); i++) {
      out += vec[i];
    }
    return out / ((utype) vec.size());
  }

  template <class utype> utype mean(const deque<utype> vec) {
    if (vec.size() == 0) {
      return (utype) 0;
    }
    utype out = 0;
    for (size_t i = 0; i < vec.size(); i++) {
      out += vec[i];
    }
    return out / ((utype) vec.size());
  }

  template <class utype> utype mean(const vector<vector<utype>> mat) {
    if (mat.size() == 0) {
      return (utype) 0;
    }
    if (mat.at(0).size() == 0) {
      return (utype) 0;
    }
    utype out = 0;
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        out += mat[i][j];
      }
    }
    return out / ((utype) mat.size() * mat.at(0).size());
  }

#define AST_TEMPLATE(utype) template utype mean(const vector<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype mean(const deque<utype>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template utype mean(const vector<vector<utype>>);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ***************************************************************************
  // Function reset of a vector/deque
  // ***************************************************************************
  template <class utype> vector<utype> reset(vector<utype>& vec) {
    for (size_t i = 0; i < vec.size(); i++) {
      vec[i] = (utype) 0;
    }
    return vec;
  }

  template <class utype> deque<utype> reset(deque<utype>& vec) {
    for (size_t i = 0; i < vec.size(); i++) {
      vec[i] = (utype) 0;
    }
    return vec;
  }

  // ***************************************************************************
  // Function reset of a vector<vector<>>
  // ***************************************************************************
  template <class utype> vector<vector<utype>> reset(vector<vector<utype>> mat) {
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        mat[i][j] = (utype) 0;
      }
    }
    return mat;
  }

  // ***************************************************************************
  // Function clear of a vector/deque
  // ***************************************************************************
  template <class utype> vector<utype> clear(vector<utype>& vec) {
    for (size_t i = 0; i < vec.size(); i++) {
      vec[i] = (utype) 0;
    }
    return vec;
  }

  template <class utype> deque<utype> clear(deque<utype>& vec) {
    for (size_t i = 0; i < vec.size(); i++) {
      vec[i] = (utype) 0;
    }
    return vec;
  }

  // ***************************************************************************
  // Function clear of a vector<vector<>>
  // ***************************************************************************
  template <class utype> vector<vector<utype>> clear(vector<vector<utype>> mat) {
    for (size_t i = 0; i < mat.size(); i++) {
      for (size_t j = 0; j < mat[i].size(); j++) {
        mat[i][j] = (utype) 0;
      }
    }
    return mat;
  }

  // ***************************************************************************
  // Function random_shuffle of a vector/deque
  // ***************************************************************************
  template <class utype> void random_shuffle(vector<utype>& vec) {
    // switch std::random_shuffle to not deprecated std::shuffle //HE20230620
    // https://en.cppreference.com/w/cpp/algorithm/random_shuffle
    std::random_device rd;
    std::shuffle(vec.begin(), vec.end(), rd);
  }

  template <class utype> void random_shuffle(deque<utype>& vec) {
    // switch std::random_shuffle to not deprecated std::shuffle //HE20230620
    // https://en.cppreference.com/w/cpp/algorithm/random_shuffle
    std::random_device rd;
    std::shuffle(vec.begin(), vec.end(), rd);
  }

#define AST_TEMPLATE(utype) template void random_shuffle(vector<utype>& vec);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template void random_shuffle(deque<utype>& vec);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE

  // ***************************************************************************
  // Function isequal of vector vector
  // ***************************************************************************
  template <class utype> bool identical(vector<utype> vec1, vector<utype> vec2, utype epsilon) {
    if (vec1.size() != vec2.size()) {
      return false;
    }
    for (size_t i = 0; i < vec1.size(); i++) {
      if (aurostd::abs(vec1[i] - vec2[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool identical(vector<utype>, vector<utype>, utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> bool identical(deque<utype> vec1, deque<utype> vec2, utype epsilon) {
    if (vec1.size() != vec2.size()) {
      return false;
    }
    for (size_t i = 0; i < vec1.size(); i++) {
      if (aurostd::abs(vec1[i] - vec2[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }
#define AST_TEMPLATE(utype) template bool identical(deque<utype>, deque<utype>, utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  bool identical(vector<int> vec1, vector<int> vec2, int epsilon) {
    if (vec1.size() != vec2.size()) {
      return false;
    }
    for (size_t i = 0; i < vec1.size(); i++) {
      if (aurostd::abs(vec1[i] - vec2[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  bool identical(deque<int> vec1, deque<int> vec2, int epsilon) {
    if (vec1.size() != vec2.size()) {
      return false;
    }
    for (size_t i = 0; i < vec1.size(); i++) {
      if (aurostd::abs(vec1[i] - vec2[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  // ***************************************************************************
  // Function isequal of vector<vector<>> vector<vector<>>
  // ***************************************************************************
  template <class utype> bool identical(const vector<vector<utype>>& mat1, const vector<vector<utype>>& mat2, utype epsilon) {
    if (mat1.size() != mat2.size()) {
      return false;
    }
    for (size_t i = 0; i < mat1.size(); i++) {
      if (mat1[i].size() != mat2[i].size()) {
        return false;
      }
      for (size_t j = 0; j < mat1[i].size(); j++) {
        if (aurostd::abs(mat1[i][j] - mat2[i][j]) > epsilon) {
          return false;
        }
      }
    }
    return true;
  }

  // ***************************************************************************
  // Function isequal of vector<vector<vector<>>> vector<vector<vector<>>>
  // ***************************************************************************
  template <class utype> bool identical(const vector<vector<vector<utype>>>& t1, const vector<vector<vector<utype>>>& t2, utype epsilon) {
    if (t1.size() != t2.size()) {
      return false;
    }
    for (size_t i = 0; i < t1.size(); i++) {
      if (t1[i].size() != t2[i].size()) {
        return false;
      }
      for (size_t j = 0; j < t1[i].size(); j++) {
        if (t1[i][j].size() != t2[i][j].size()) {
          return false;
        }
        for (size_t k = 0; k < t1[i][i].size(); k++) {
          if (aurostd::abs(t1[i][j][k] - t2[i][j][k]) > epsilon) {
            return false;
          }
        }
      }
    }
    return true;
  }

  // ***************************************************************************
  // Function all identical vector (one argument, checks if all entries are equal) //DX20210422
  // ***************************************************************************
  template <class utype> bool identical(const vector<utype>& vec, utype eps) {
    for (size_t i = 0; i < vec.size(); i++) {
      if (isdifferent(vec[i], vec[0], eps)) {
        return false;
      }
    }
    return true; // includes case when vec is empty
  }
#define AST_TEMPLATE(utype) template bool identical(const vector<utype>&, utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ***************************************************************************
  // Function all identical deque (one argument, checks if all entries are equal) //DX20210422
  // ***************************************************************************
  template <class utype> bool identical(const deque<utype>& deq, utype eps) {
    for (size_t i = 0; i < deq.size(); i++) {
      if (isdifferent(deq[i], deq[0], eps)) {
        return false;
      }
    }
    return true; // includes case when deq is empty
  }
#define AST_TEMPLATE(utype) template bool identical(const deque<utype>&, utype);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ***************************************************************************
  // Function toupper/tolower
  // ***************************************************************************
  string toupper(const string& in) {
    string out(in);
    std::transform(out.begin(), out.end(), out.begin(), ::toupper);
    return out;
  }

  string tolower(const string& in) {
    string out(in);
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
  }

  char toupper(const char& in) {
    return std::toupper(in);
  }

  char tolower(const char& in) {
    return std::tolower(in);
  }

  // ***************************************************************************
  // Function getEveryNth //SD20230216
  // ***************************************************************************
  template <class utype> vector<utype> getEveryNth(const vector<utype>& vec, uint n) {
    if (n == 0) {
      string message = "Cannot increment by 0";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    vector<utype> vec_n;
    for (size_t i = 0; i < vec.size(); i += n) {
      vec_n.push_back(vec[i]);
    }
    return vec_n;
  }

  vector<double> getEveryNth(const vector<double>& vec, uint n) {
    if (n == 0) {
      string message = "Cannot increment by 0";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    vector<double> vec_n;
    for (size_t i = 0; i < vec.size(); i += n) {
      vec_n.push_back(vec[i]);
    }
    return vec_n;
  }

//#define AST_TEMPLATE(utype) template<class utype> vector<utype> getEveryNth(const vector<utype>& vec, uint n);
//  AST_GEN_1(AST_UTYPE_NUM)
//  AST_GEN_1(AST_UTYPE_STRING)
//#undef AST_TEMPLATE

  // ***************************************************************************
  // Function GetNumFields
  // ***************************************************************************
  // Dane Morgan
  int GetNumFields(const string& s) {
    int nf = 0;
    int in_a_field = 0;
    for (size_t i = 0; i < s.size(); i++) {
      if (!(in_a_field) && s[i] != ' ') {
        in_a_field = 1;
        nf++;
      }
      if (in_a_field && s[i] == ' ') {
        in_a_field = 0;
      }
    }
    return nf;
  }

  // ***************************************************************************
  // Function GetNextVal
  // ***************************************************************************
  // Dane Morgan - Stefano Curtarolo
  string GetNextVal(const string& s, int& id) {
    string ss;
    int i = id;
    while (s[i] == ' ') {
      i++;
    } // ignore leading spaces.
    while (i < (int) s.size() && s[i] != ' ') { // pull out all text until next space.
      ss += s[i];
      i++;
    }
    id = i;
    return ss;
  }

  // ***************************************************************************
  // Function PaddedNumString
  // ***************************************************************************
  // Dane Morgan - Stefano Curtarolo
  string PaddedNumString(const int num, const int ndigits) {
    ostringstream oss;
    oss << std::setw(ndigits) << std::setfill('0') << num;// << ends;
    return oss.str();
  }

  // ***************************************************************************
  // Function getZeroPadding
  // ***************************************************************************
  // Corey Oses
  int getZeroPadding(double d) {
    return int(log10(d)) + 1;
  }
  int getZeroPadding(int num) {
    return getZeroPadding((double) num);
  }
  int getZeroPadding(uint num) {
    return getZeroPadding((double) num);
  }
  int getZeroPadding(long int num) {
    return getZeroPadding((double) num);
  }
  int getZeroPadding(unsigned long int num) {
    return getZeroPadding((double) num);
  }
  int getZeroPadding(long long int num) {
    return getZeroPadding((double) num);
  }
  int getZeroPadding(unsigned long long int num) {
    return getZeroPadding((double) num);
  }

  // ***************************************************************************
  // Function PaddedPRE
  // ***************************************************************************
  // Add PRE characters to pad
  string PaddedPRE(string input, int depth, string ch) {
    stringstream aus("");
    aus << input;
    string strout;
    for (int i = 0; i < depth - (int) aus.str().length(); i++) {
      strout += ch;
    }
    strout += aus.str();
    return strout;
  }
  template <class utype> string PaddedPRE(utype input, int depth, string ch) {
    stringstream sss;
    sss << input;
    return PaddedPRE(sss.str(), depth, ch);
  }
#define AST_TEMPLATE(utype) template string PaddedPRE(utype input, int depth, string ch);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_TEXT)
#undef AST_TEMPLATE

  // ***************************************************************************
  // Function PaddedPOST
  // ***************************************************************************
  // Add POST characters to pad
  string PaddedPOST(string input, int depth, string ch) {
    stringstream aus("");
    aus << input;
    string strout = aus.str();
    for (int i = 0; i < depth - (int) aus.str().length(); i++) {
      strout += ch;
    }
    return strout;
  }
  template <class utype> string PaddedPOST(utype input, int depth, string ch) {
    stringstream sss;
    sss << input;
    return PaddedPOST(sss.str(), depth, ch);
  }

#define AST_TEMPLATE(utype) template string PaddedPOST(utype, int, string);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_TEXT)
#undef AST_TEMPLATE

  // ***************************************************************************
  // Function PaddedCENTER
  // ***************************************************************************
  // Add PRE AND POST characters to pad so that string is in the center
  string PaddedCENTER(string input, int depth, string ch) {
    stringstream aus("");
    const int pre = (depth - (int) input.length()) / 2;
    const int post = depth - pre - (int) input.length();
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: input.length()=" << input.length() << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: pre=" << pre << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: post=" << post << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: depth=" << depth << endl;
    // if(DEBUG) cerr << "aurostd::PaddedCENTER: pre+post+input.length()=" << pre+post+input.length() << endl;
    for (int i = 1; i < pre; i++) {
      aus << ch;
    }
    aus << input;
    for (int i = 1; i < post; i++) {
      aus << ch;
    }
    return aus.str();
  }
  template <class utype> string PaddedCENTER(utype input, int depth, string ch) {
    stringstream sss;
    sss << input;
    return PaddedCENTER(sss.str(), depth, ch);
  }

  // ***************************************************************************
  // Function ProgressBar
  // ***************************************************************************
  uint ProgressBar(std::ostream& oss, string prelim, uint j, uint jmax, bool VERBOSE_PERCENTAGE, bool VERBOSE_ROLLER, bool VERBOSE_CURSOR) {
    uint position = 0;
    const double percentage = double(j) / double(jmax);
    if (j == 0) {
      oss << prelim; // position+=prelim.size();
    }
    // VERBOSE PERCENTAGE
    if (!mod<uint>(j, 50) || j == jmax - 1 || j == jmax) {
      if (VERBOSE_PERCENTAGE) {
        if (j == jmax - 1 || j == jmax) {
          oss << "[100.0%]";
          position += 8;
        } else {
          // 99.99%
          oss << "[" << (percentage < 0.1 ? " " : "") << mod<uint>(uint(percentage * 100), 100) << "." << mod<uint>(uint(percentage * 1000), 10) << mod<uint>(uint(percentage * 10000), 10) << "%]";
          position += 8;
        }
        oss << " ";
        position++;
      }
    }
    // VERBOSE_ROLLER
    if (!mod<uint>(j, 50) || j == jmax - 1 || j == jmax) {
      if (VERBOSE_ROLLER) {
        if (j == jmax - 1 || j == jmax) {
          oss << "[=]";
        } else {
          if (mod<uint>(j / 513, 4) == 0) {
            oss << "[\\]";
          }
          if (mod<uint>(j / 513, 4) == 1) {
            oss << "[|]";
          }
          if (mod<uint>(j / 513, 4) == 2) {
            oss << "[/]";
          }
          if (mod<uint>(j / 513, 4) == 3) {
            oss << "[-]";
          }
        }
        position += 3;
        oss << " ";
        position++;
      }
    }
    // VERBOSE CURSOR
    if (j == 0 || !mod<uint>(j, 50) || j == jmax - 1 || j == jmax) {
      if (VERBOSE_CURSOR) {
        if (j == jmax - 1 || j == jmax) {
          oss << "[======================================================================================================]";
          position += 102;
        } else {
          oss << "[";
          position++;
          for (double k = 0; k < percentage * 100; k += 1.0) {
            oss << "=";
            position++;
          }
          // if(mod<uint>(j,500)==0)
          {
            if (mod<uint>(j / 478, 4) == 0) {
              oss << "\\";
              position++;
            }
            if (mod<uint>(j / 478, 4) == 1) {
              oss << "|";
              position++;
            }
            if (mod<uint>(j / 478, 4) == 2) {
              oss << "/";
              position++;
            }
            if (mod<uint>(j / 478, 4) == 3) {
              oss << "-";
              position++;
            }
          }
          if (j == 0) {
            for (double k = 0; k < (1.0 - percentage) * 100.0 + 0.01; k += 1.0) {
              oss << " ";
              position++;
            }
            oss << "]";
            position++;
          }
        }
        oss << " ";
        position++;
      }
      // NOW GO BACK
      for (uint k = 0; k < position; k++) {
        oss << "\b";
      }
      if (j == jmax - 1 || j == jmax) {
        oss << endl;
      }
    }
    return position;
  }

  uint ProgressBar(std::ostream& oss, string prelim, uint j, uint jmax) {
    return ProgressBar(oss, prelim, j, jmax, true, true, true);
  }

  uint ProgressBar(std::ostream& oss, string prelim, double j, bool VERBOSE_PERCENTAGE, bool VERBOSE_ROLLER, bool VERBOSE_CURSOR) {
    return ProgressBar(oss, prelim, uint(double(j * 100)), 100, VERBOSE_PERCENTAGE, VERBOSE_ROLLER, VERBOSE_CURSOR);
  }

  uint ProgressBar(std::ostream& oss, string prelim, double j) {
    return ProgressBar(oss, prelim, uint(double(j * 100)), 100, true, true, true);
  }

  // ***************************************************************************
  // Function PercentEncodeASCII //DX20210706
  // ***************************************************************************
  // Converts a single ASCII character into its percent-encoded form
  string PercentEncodeASCII(const char c) {
    stringstream char_percent_encoded;
    char_percent_encoded << "%" << std::hex << (int) c;
    return char_percent_encoded.str();
  }

  // ***************************************************************************
  // Function CleanStringASCII
  // ***************************************************************************
  // Clean a string from ASCII junk
  // Stefano Curtarolo
  string CleanStringASCII(const string& s) { // CO20190712
    string ss = s;
    CleanStringASCII_InPlace(ss);
    return ss;
  }

  // ***************************************************************************
  // Function CleanStringASCII_InPlace
  // ***************************************************************************
  // Similar to CleanStringASCII, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void CleanStringASCII_InPlace(string& s) {
    //[CO20190712 - slight optimization if we go backwards]for(uint i=0;i<s.length();i++)
    for (uint i = s.length() - 1; i < s.length(); i--) { // CO20200106 - patching for auto-indenting
      //[CO20200624 - not inclusive enough]if(!(
      //[CO20200624 - not inclusive enough]      (s[i]>='A' && s[i]<='Z') || //LETTERS
      //[CO20200624 - not inclusive enough]      (s[i]>='a' && s[i]<='z') || //letters
      //[CO20200624 - not inclusive enough]      (s[i]>='0' && s[i]<='9') || //numbers
      //[CO20200624 - not inclusive enough]      (s[i]=='.' || s[i]=='+' || s[i]=='-' || s[i]=='*' || s[i]=='/') ||  //operations
      //[CO20200624 - not inclusive enough]      (s[i]=='_' || s[i]=='#' || s[i]=='&' || s[i]==':' || s[i]==',' || s[i]=='@' || s[i]=='$') ||  //punctuation1
      //[CO20200624 - not inclusive enough]      (s[i]=='=' || s[i]=='|' || s[i]=='\'' || s[i]=='\"' || s[i]==' ') ||  //punctuation2
      //[CO20200624 - not inclusive enough]      false)
      //[CO20200624 - not inclusive enough]  ){RemoveCharacterInPlace(s,s[i]);}
      // https://stackoverflow.com/questions/48212992/how-to-find-out-if-there-is-any-non-ascii-character-in-a-string-with-a-file-path
      // cerr << s[i] << " " << static_cast<unsigned int>(s[i]) << endl;
      if (static_cast<unsigned int>(s[i]) > 127) {
        RemoveCharacterInPlace(s, s[i]);
      }
    }
  }

  // ***************************************************************************
  // Function RemoveTrailingCharacter
  // ***************************************************************************
  // Removes trailing character
  // CO+ME20200825
  string RemoveTrailingCharacter(const string& s, char c) {
    string ss = s;
    RemoveTrailingCharacter_InPlace(ss, c);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTrailingCharacter_InPlace
  // ***************************************************************************
  // Similar to RemoveTrailingCharacter, but does NOT create a new string (costly if done MANY times)
  // CO+ME20200825
  void RemoveTrailingCharacter_InPlace(string& s, char c) {
    while (!s.empty() && s.at(s.size() - 1) == c) {
      s = s.substr(0, s.size() - 1);
    }
  }

  // DX20190516 - remove control code characters - START
  //  ***************************************************************************
  //  Function removeControlCodeCharactersFromString
  //  ***************************************************************************
  bool RemoveControlCodeCharactersFromString(const string& in, string& out) {  // CO20190620

    // removes control code and backspace characters (e.g., NUL, DEL, etc.)
    // only keep printable characters (i.e., digits, letters, punctuation, and spaces)
    // and white space characters (i.e., space, newline, tabs, and carrage returns)
    // a boolean indicates if the stringstream contained a control code character
    // string input version

    stringstream ss_in;
    stringstream ss_out;
    ss_in << in;
    const bool detected_control_char = RemoveControlCodeCharactersFromStringstream(ss_in, ss_out);
    out = ss_out.str();
    return detected_control_char;
  }

  // ***************************************************************************
  // Function removeControlCodeCharactersFromStringStream
  // ***************************************************************************
  bool RemoveControlCodeCharactersFromStringstream(std::stringstream& ss_in, std::stringstream& ss_out) {
    // removes control code and backspace characters (e.g., NUL, DEL, etc.)
    // only keep printable characters (i.e., digits, letters, punctuation, and spaces)
    // and white space characters (i.e., space, newline, tabs)
    // a boolean indicates if the stringstream contained a control code character
    // stringstream input version
    //
    // ME20190614: We don't want carriage returns either because they mess up string additions.
    // Since they point to the beginning of the string, adding to a string with a carriage
    // return would overwrite instead of append

    bool detected_control_char = false;
    char c;
    // char c1;
    // char c2;

    // stringstream tmp; tmp << ss_in.str();
    while (ss_in.get(c)) {
      //[CO20190620 - still doesn't work]if(isprint(c) || isspace(c) || (c != '\r')) {  //ME20190614 //[CO20200106 - close bracket for indenting]}
      if ((isprint(c) || isspace(c) || false) && ((c != '\r') || false)) {
        ss_out << c;
      }  // CO20190620 - add more cases before false
      else {
        detected_control_char = true;
      }
    }

    return detected_control_char;
  }
  // DX20190516 - remove control code characters - END

  // DX20190211 - remove control code characters from file - START
  //  ***************************************************************************
  //  Function RemoveControlCodeCharactersFromFile
  //  ***************************************************************************
  bool RemoveControlCodeCharactersFromFile(const string& directory, const string& filename, bool keep_orig_file) {  // CO20210315

    // removes control code and backspace characters (e.g., NUL, DEL, etc.)
    // overwrites file if control characters are detected, otherwise the file is untouched (preserve original timestamp)
    // matches original compression

    if (aurostd::FileExist(directory + "/" + filename)) {
      stringstream ss_in;
      stringstream ss_out;
      aurostd::compressfile2stringstream(directory + "/" + filename, ss_in);
      // if file contains control code characters, then overwrite file
      if (RemoveControlCodeCharactersFromStringstream(ss_in, ss_out)) {
        string uncompressed_filename;
        // compressed files
        if (IsCompressed(filename, uncompressed_filename)) {
          stringstream2file(ss_out, directory + "/" + uncompressed_filename + "_tmp");
          const compression_type ct = GetCompressionType(filename);
          CompressFile(directory + "/" + uncompressed_filename + "_tmp", ct);
          if (keep_orig_file) {
            file2file(directory + "/" + filename, directory + "/" + uncompressed_filename + "_old" + compression_suffix[ct]);
          } // move original
          file2file(directory + "/" + uncompressed_filename + "_tmp" + compression_suffix[ct], directory + "/" + filename); // overwrite
        }
        // uncompressed files
        else {
          stringstream2file(ss_out, directory + "/" + filename + ".tmp");
          if (keep_orig_file) {
            file2file(directory + "/" + filename, directory + "/" + filename + "_old");
          } // move original
          file2file(directory + "/" + filename + ".tmp", directory + "/" + filename); // overwrite
        }
        return true;
      }
      // file is ok, do not update
      else {
        return false;
      } // signals file is unchanged
    } else {
      const string message = "File does not exist: " + directory + "/" + filename;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    return false;
  }
  // DX20190211 - remove control code characters from file - END

  // DX20190125 - remove null bytes - START
  //  ***************************************************************************
  //  Function isNullbyte
  //  ***************************************************************************
  //  Deterine if char is a null byte (e.g., ^@)
  bool isNullByte(char c) {
    return (c == '\0');
  }

  // ***************************************************************************
  // Function removeNullBytes
  // ***************************************************************************
  // Remove all null bytes in string (e.g., ^@)
  string removeNullBytes(string in) {
    string out = in;
    out.erase(remove_if(out.begin(), out.end(), isNullByte), out.end());
    return out;
  }
  // DX20190125 - remove null bytes - END

  // DX20190211 - remove null characters from file - START
  //  ***************************************************************************
  //  Function RemoveBinaryCharactersFromFile()
  //  ***************************************************************************
  //  Remove all null bytes in file
  bool RemoveBinaryCharactersFromFile(const string& directory, const string& filename) {  // CO20210315
    stringstream aus_exec;
    const vector<string> vext{"", ".bz2", ".xz", ".gz"};
    const vector<string> vcmd{"cat", "bzcat", "xzcat", "gzcat"};
    const vector<string> vzip{"", "bzip2", "xz", "gzip"};
    if (vext.size() != vcmd.size()) {
      const string message = "vext.size()!=vcmd.size()";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    for (size_t iext = 0; iext < vext.size(); iext++) { // check filename.EXT
      if (aurostd::FileExist(directory + "/" + filename + vext[iext])) {
        aus_exec << "cd \"" << directory << "\"" << endl;
        aus_exec << vcmd[iext] << " " << filename << vext[iext] << R"( | sed "s/[^[:print:]\r\t]//g" > )" << filename << ".tmp && mv " << filename << ".tmp " << filename << endl;
        if (!vext[iext].empty()) {
          aus_exec << vzip[iext] << " " << filename << endl;
        }
        aurostd::execute(aus_exec);
      }
    }
    return true;
  }
  // DX20190211 - remove binary characters from file - END

  // ***************************************************************************
  // Function CGI_StringClean
  // ***************************************************************************
  // Clean a string from CGI junk
  string CGI_StringClean(const string& stringIN) {
    string stringOUT = stringIN;
    aurostd::StringSubstInPlace(stringOUT, "%0D%0A", "\n");
    aurostd::StringSubstInPlace(stringOUT, "%0d%0a", "\n");   // newlines
    aurostd::StringSubstInPlace(stringOUT, "+", " ");    // spaces
    aurostd::StringSubstInPlace(stringOUT, "%28", "(");
    aurostd::StringSubstInPlace(stringOUT, "%29", ")");   // ()
    aurostd::StringSubstInPlace(stringOUT, "%5B", "[");
    aurostd::StringSubstInPlace(stringOUT, "%5D", "]");   // []
    aurostd::StringSubstInPlace(stringOUT, "%7B", "{");
    aurostd::StringSubstInPlace(stringOUT, "%7D", "}");   // brackets (do not write, it screws up indent)
    aurostd::StringSubstInPlace(stringOUT, "%2B", "+");
    aurostd::StringSubstInPlace(stringOUT, "%2F", "/");   //  operations
    aurostd::StringSubstInPlace(stringOUT, "%23", "#");
    aurostd::StringSubstInPlace(stringOUT, "%21", "!");
    aurostd::StringSubstInPlace(stringOUT, "%3F", "?");
    aurostd::StringSubstInPlace(stringOUT, "%2C", ",");
    aurostd::StringSubstInPlace(stringOUT, "%3A", ":");
    aurostd::StringSubstInPlace(stringOUT, "%3B", ";");
    aurostd::StringSubstInPlace(stringOUT, "%27", "'");
    aurostd::StringSubstInPlace(stringOUT, "%22", "\"");
    aurostd::StringSubstInPlace(stringOUT, "%60", "`");
    aurostd::StringSubstInPlace(stringOUT, "%40", "@");
    aurostd::StringSubstInPlace(stringOUT, "%24", "$");
    aurostd::StringSubstInPlace(stringOUT, "%25", "%");
    aurostd::StringSubstInPlace(stringOUT, "%5E", "^");
    aurostd::StringSubstInPlace(stringOUT, "%26", "&");
    aurostd::StringSubstInPlace(stringOUT, "%3D", "=");
    aurostd::StringSubstInPlace(stringOUT, "%7E", "~");
    aurostd::StringSubstInPlace(stringOUT, "%5C", "\\");
    aurostd::StringSubstInPlace(stringOUT, "%7C", "|");
    aurostd::StringSubstInPlace(stringOUT, "%3C", "<");
    aurostd::StringSubstInPlace(stringOUT, "%3E", ">");  // <>
    aurostd::StringSubstInPlace(stringOUT, "\n\n", "\n");
    return stringOUT;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpaces
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpaces(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss;
    for (size_t i = 0; i < s.size(); i++) {
      if (s[i] != ' ' && s[i] != '\t') {
        ss += s[i];
      }
    }
    return ss;
  }
  string RemoveWhiteSpaces(const string& s, const char toggle) {  // CO20190710
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss;
    bool copy = true;
    for (size_t i = 0; i < s.size(); i++) {
      if (s[i] == toggle) {
        copy = !copy;
      }  // CO20190710
      if (copy) {
        if (s[i] != ' ' && s[i] != '\t') {
          ss += s[i];
        }
      } else {
        ss += s[i];
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveWhiteSpacesFromTheBack(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    while (ss[ss.size() - 1] == ' ' || ss[ss.size() - 1] == '\t') {
      ss.erase(ss.size() - 1, 1);
      if (ss.empty()) {
        return ss;  // nothing to do
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheFront
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Oses
  string RemoveWhiteSpacesFromTheFront(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    while (ss[0] == ' ' || ss[ss.size() - 1] == '\t') {
      ss.erase(0, 1);
      if (ss.empty()) {
        return ss;  // nothing to do
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveWhiteSpacesFromTheFrontAndBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Oses
  string RemoveWhiteSpacesFromTheFrontAndBack(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    ss = RemoveWhiteSpacesFromTheBack(ss);
    ss = RemoveWhiteSpacesFromTheFront(ss);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSpaces
  // ***************************************************************************
  // Removes all spaces from a string. Morgan / Curtarolo
  string RemoveSpaces(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss;
    for (size_t i = 0; i < s.size(); i++) {
      if (s[i] != ' ') {
        ss += s[i];
      }
    }
    return ss;
  }
  string RemoveSpaces(const string& s, const char toggle) { // CO20190710
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss;
    bool copy = true;
    for (size_t i = 0; i < s.size(); i++) {
      if (s[i] == toggle) {
        copy = !copy;
      }  // CO20190710
      if (copy) {
        if (s[i] != ' ') {
          ss += s[i];
        }
      } else {
        ss += s[i];
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSpacesFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string. Morgan / Curtarolo
  string RemoveSpacesFromTheBack(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    while (ss[ss.size() - 1] == ' ') {
      ss.erase(ss.size() - 1, 1);
      if (ss.empty()) {
        return ss;  // nothing to do
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTabs
  // ***************************************************************************
  // Removes all tabs from a string. Stefano Curtarolo
  string RemoveTabs(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss;
    for (size_t i = 0; i < s.size(); i++) {
      if (s[i] != '\t') {
        ss += s[i];
      }
    }
    return ss;
  }
  string RemoveTabs(const string& s, const char toggle) { // CO20190710
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss;
    bool copy = true;
    for (size_t i = 0; i < s.size(); i++) {
      if (s[i] == toggle) {
        copy = !copy;
      }  // CO20190710
      if (copy) {
        if (s[i] != '\t') {
          ss += s[i];
        }
      } else {
        ss += s[i];
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveTabsFromTheBack
  // ***************************************************************************
  // Removes all white spaces (spaces, tabs) from a string.
  // Dane Morgan / Stefano Curtarolo
  string RemoveTabsFromTheBack(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    while (ss[ss.size() - 1] == '\t') {
      ss.erase(ss.size() - 1, 1);
      if (ss.empty()) {
        return ss;  // nothing to do
      }
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveComments
  // ***************************************************************************
  // Removes all comments from a string.
  // Stefano Curtarolo
  //   string RemoveComments(const string& s) {
  //     if(s.size()==0) return s;  // nothing to do
  //     string ss;
  //     bool copy=true;
  //     for (size_t i=0;i<s.size();i++) {
  //       if(s[i]=='#')  copy=false;
  //       if(s[i]=='\n') copy=true;
  //       if(copy) ss+=s[i];
  //     }
  //     return ss;
  //   }

  // ME20190614 - added vector<string> version of RemoveComments
  vector<string> RemoveComments(const vector<string>& vstrin) { // CO20210315 - cleaned up
    vector<string> vstrout;
    string::size_type loc;
    string line;
    for (size_t i = 0; i < vstrin.size(); i++) {
      line = vstrin[i];
      // COMMENT_NEGLECT_1
      loc = line.find(COMMENT_NEGLECT_1);
      while (loc != string::npos) {
        // Do not remove #[1-9] since it is not a comment (spacegroup)
        if (!((loc > 0) && (loc < line.size()) && (isdigit(line[loc + 1])))) {
          line = line.substr(0, loc);
          break;
        }
        loc = line.find(COMMENT_NEGLECT_1, loc + 1);
      }
      // COMMENT_NEGLECT_2
      loc = line.find(COMMENT_NEGLECT_2);
      while (loc != string::npos) {
        // Do not remove :// since it is not a comment (web address)
        if (!((loc > 0) && (loc < line.size()) && (line[loc - 1] == ':'))) {
          line = line.substr(0, loc);
          break;
        }
        loc = line.find(COMMENT_NEGLECT_2, loc + 1);
      }
      // COMMENT_NEGLECT_3
      loc = line.find(COMMENT_NEGLECT_3);
      line = line.substr(0, loc);
      if (!line.empty()) {
        vstrout.push_back(line);
      }
    }
    return vstrout;
  }
  deque<string> RemoveComments(const deque<string>& vstrin) {
    return aurostd::vector2deque(RemoveComments(aurostd::deque2vector(vstrin)));
  }

  string RemoveComments(const string& strin) {  // CO20210315 - cleaned up
    vector<string> vlines;
    aurostd::string2vectorstring(strin, vlines);
    vlines = RemoveComments(vlines);
    if (vlines.empty()) {
      return "";
    }
    if (vlines.size() == 1) {
      return vlines[0];
    }
    string strout;
    for (size_t i = 0; i < vlines.size(); i++) {
      strout += vlines[i] + '\n';
    }
    return strout;
  }

  // ***************************************************************************
  // Function RemoveCharacter
  // ***************************************************************************
  // Removes charecters from string
  // Stefano Curtarolo
  string RemoveCharacter(const string& s, const char character) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    RemoveCharacterInPlace(ss, character);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveCharacterInPlace
  // ***************************************************************************
  // Similar to RemoveCharacter, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveCharacterInPlace(string& t, const char character) {
    t.erase(std::remove(t.begin(), t.end(), character), t.end());
  }

  // ***************************************************************************
  // Function RemoveCharacterFromTheBack
  // ***************************************************************************
  // Remove character from the back of a string. DX (Hicks)
  string RemoveCharacterFromTheBack(const string& s, const char character) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    if (ss[ss.size() - 1] == character) {
      ss.erase(ss.size() - 1, 1);
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveCharacterFromTheFront
  // ***************************************************************************
  // Removes character from the front of a string. DX (Hicks)
  string RemoveCharacterFromTheFront(const string& s, const char character) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    if (ss[0] == character) {
      ss.erase(0, 1);
    }
    return ss;
  }

  // ***************************************************************************
  // Function RemoveCharacterFromTheFrontAndBack
  // ***************************************************************************
  // Removes character from the front and back of a string. DX (Hicks)
  string RemoveCharacterFromTheFrontAndBack(const string& s, const char character) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    ss = RemoveCharacterFromTheBack(ss, character);
    if (s.empty()) {
      return s; // cannot remove anything else
    }
    ss = RemoveCharacterFromTheFront(ss, character);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveNumbers
  // ***************************************************************************
  // Removes numbers from string
  // Stefano Curtarolo
  string RemoveNumbers(const string& s) {  // CO20190712 - avoids creating many copies of string
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    RemoveNumbersInPlace(ss);
    return ss;
  }

  // ***************************************************************************
  // Function RemoveNumbersInPlace
  // ***************************************************************************
  // Similar to RemoveNumbers, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveNumbersInPlace(string& s) { // CO20190712
    RemoveCharacterInPlace(s, '0');
    RemoveCharacterInPlace(s, '1');
    RemoveCharacterInPlace(s, '2');
    RemoveCharacterInPlace(s, '3');
    RemoveCharacterInPlace(s, '4');
    RemoveCharacterInPlace(s, '5');
    RemoveCharacterInPlace(s, '6');
    RemoveCharacterInPlace(s, '7');
    RemoveCharacterInPlace(s, '8');
    RemoveCharacterInPlace(s, '9');
    RemoveCharacterInPlace(s, '.');
  }

  // ***************************************************************************
  // Function RemoveRounding
  // ***************************************************************************
  // Removes rounding from string
  // Stefano Curtarolo
  string RemoveRounding(const string& s) {
    if (s.empty()) {
      return s;  // nothing to do
    }
    string ss = s;
    ss = RemoveSubString(ss, "(0)");
    ss = RemoveSubString(ss, "(1)");
    ss = RemoveSubString(ss, "(2)");
    ss = RemoveSubString(ss, "(3)");
    ss = RemoveSubString(ss, "(4)");
    ss = RemoveSubString(ss, "(5)");
    ss = RemoveSubString(ss, "(6)");
    ss = RemoveSubString(ss, "(7)");
    ss = RemoveSubString(ss, "(8)");
    ss = RemoveSubString(ss, "(9)");
    return ss;
  }

  // ***************************************************************************
  // Function RemoveSubStringFirst
  // ***************************************************************************
  // Removes the first substring from string
  // Stefano Curtarolo
  string RemoveSubStringFirst(const string& str_orig, const string& str_rm) {
    string t = str_orig;
    RemoveSubStringFirstInPlace(t, str_rm);  // CO20190712
    return t;
  }

  // ***************************************************************************
  // Function RemoveSubStringFirstInPlace
  // ***************************************************************************
  // Similar to RemoveSubStringFirst, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveSubStringFirstInPlace(string& t, const string& str_rm) {
    const std::string::size_type i = t.find(str_rm);
    if (i != std::string::npos) {
      t.erase(i, str_rm.length());
    }
  }

  // ***************************************************************************
  // Function RemoveSubString
  // ***************************************************************************
  // Removes all instances of substring from string
  // Stefano Curtarolo
  string RemoveSubString(const string& str_orig, const string& str_rm) {
    string t = str_orig;
    RemoveSubStringInPlace(t, str_rm); // CO20190712
    return t;
  }

  // ***************************************************************************
  // Function RemoveSubStringInPlace
  // ***************************************************************************
  // Similar to RemoveSubString, but does NOT create a new string (costly if done MANY times)
  // Corey Oses 20190712
  void RemoveSubStringInPlace(string& t, const string& str_rm) {
    string::size_type i = t.find(str_rm); // CO20190712 - fewer operations
    while (i != string::npos) {  // CO20190712 - fewer operations
      t.erase(i, str_rm.length()); // CO20190712 - fewer operations
      i = t.find(str_rm); // CO20190712 - fewer operations
    }
  }

  // ***************************************************************************
  // Function VersionString2Double
  // ***************************************************************************
  // 5.1.311 -> 5.0013311
  // 4.2.34 -> 4.002034
  double VersionString2Double(const string& version_str) { // SD20220331
    vector<string> tokens;
    aurostd::string2tokens(version_str, tokens, ".");
    double version = 0.0;
    for (size_t i = 0; i < tokens.size(); i++) {
      version += aurostd::string2utype<double>(tokens[i]) * std::pow(10.0, -3.0 * i);
    }
    return version;
  }

  // ***************************************************************************
  // Function ProcessPIDs
  // ***************************************************************************
  // CO20210315
  vector<string> ProcessPIDs(const string& process, bool user_specific) { // CO20210315
    string output_syscall;
    return ProcessPIDs(process, output_syscall, user_specific);
  }
  vector<string> ProcessPIDs(const string& process, string& output_syscall, bool user_specific) { // CO20210315
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " looking for process=" << process << endl;
    }

    string command;
    vector<string> vlines;
    vector<string> vtokens;
    vector<string> vpids;
    uint i = 0;
    uint j = 0;
    if (aurostd::IsCommandAvailable("pgrep")) {
      string command_pgrep = "pgrep -a";  // needed over -l on fossies installs: https://fossies.org/linux/procps-ng/pgrep.1
      if (!aurostd::execute2OutErrPair(command_pgrep + " test").second.empty()) {
        command_pgrep = "pgrep -l";
      }  // should work on all linux  //the "test" is a dummy to see if the pgrep command works
      if (aurostd::execute2OutErrPair(command_pgrep + " test").second.empty()) { // the "test" is a dummy to see if the pgrep command works
        command = command_pgrep; // the -a/-l is important, we will need to neglect the subshell call below
        if (user_specific && !XHOST.user.empty()) {
          command += " -u " + XHOST.user;
        }
        command += " -f " + process + " 2> /dev/null";  // the -f is important, will match mpivasp46s in /usr/bin/mpivasp46s
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;
        }
        const string output = output_syscall = aurostd::execute2string(command);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " pgrep output:" << endl << "\"" << output << "\"" << endl;
        }
        aurostd::string2vectorstring(output, vlines);
        for (i = 0; i < vlines.size(); i++) {
          aurostd::string2tokens(vlines[i], vtokens, " ");
          if (vtokens.size() < 2) {
            continue;
          }
          const string& pid = vtokens[0];
          string proc = vtokens[1]; // since we split on " ", we need to join columns 11-onward
          for (j = 2; j < vtokens.size(); j++) {
            proc += " " + vtokens[j];
          }
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " proc[i=" << i << "]=\"" << proc << "\"" << endl;
          }
          if (proc.find(process) == string::npos) {
            continue;
          }
          if (proc.find(command) != string::npos) {
            continue;
          } // ps aux | grep ... always returns itself, neglect  //do a find() instead of == here
          vpids.push_back(pid);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " vpids=" << aurostd::joinWDelimiter(vpids, ",") << endl;
          cerr << __AFLOW_FUNC__ << " vpids.empty()=" << vpids.empty() << endl;
        }
        return vpids;
      }
    }

    if (aurostd::IsCommandAvailable("ps") && aurostd::IsCommandAvailable("grep")) {
      // FR recommends ps aux vs. ps -e
      // tested on linux and mac, PIDs are in second column, process is the last column
      const string command_grep = "grep " + process;
      command = "ps";
      if (user_specific) {
        command += " ux";
      } else {
        command += " aux";
      }
      command += " 2>/dev/null | " + command_grep + " 2> /dev/null";
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;
      }
      const string output = output_syscall = aurostd::execute2string(command);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ps/grep output:" << endl << output << endl;
      }
      aurostd::string2vectorstring(output, vlines);
      for (i = 0; i < vlines.size(); i++) {
        aurostd::string2tokens(vlines[i], vtokens, " ");
        if (vtokens.size() < 11) {
          continue;
        }
        const string& pid = vtokens[1];
        string proc = vtokens[10]; // since we split on " ", we need to join columns 11-onward
        for (j = 11; j < vtokens.size(); j++) {
          proc += " " + vtokens[j];
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " proc[i=" << i << "]=\"" << proc << "\"" << endl;
        }
        if (proc.find(process) == string::npos) {
          continue;
        }
        if (proc.find(command) != string::npos) {
          continue;
        } // ps aux | grep ... always returns itself, neglect  //do a find() instead of == here
        if (proc == command_grep) {
          continue;
        } // ps aux | grep ... always returns itself, neglect  //do a == instead a find() here
        vpids.push_back(pid);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vpids=" << aurostd::joinWDelimiter(vpids, ",") << endl;
        cerr << __AFLOW_FUNC__ << " vpids.empty()=" << vpids.empty() << endl;
      }
      return vpids;
    }
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "\"pgrep\"-type command not found", _INPUT_ILLEGAL_);
    return vpids;
  }

  // SD20220329 - overload to allow for only getting the PIDs with a specific PGID
  vector<string> ProcessPIDs(const string& process, const string& pgid, string& output_syscall, bool user_specific) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (pgid.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "PGID is empty", _INPUT_ILLEGAL_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " looking for pgid=" << pgid << endl;
      cerr << __AFLOW_FUNC__ << " looking for process=" << process << endl;
    }
    if (!aurostd::IsCommandAvailable("ps") || !aurostd::IsCommandAvailable("grep")) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "\"pgrep\"-type command not found", _INPUT_ILLEGAL_);
    }
    const string ps_opts = " uid,pgid,pid,etime,pcpu,pmem,args"; // user-defined options, since just "u" or "j" might not be good enough
    string command = "ps";
    vector<string> vlines;
    vector<string> vtokens;
    vector<string> vpids;
    uint i = 0;
    uint j = 0;
    aurostd::string2tokens(ps_opts, vtokens, ",");
    const uint nopts = vtokens.size();
    const string command_grep = "grep " + process;
    if (user_specific) {
      command += " xo";
    } else {
      command += " axo";
    }
    command += ps_opts;
    if (!aurostd::execute2OutErrPair(command + " > /dev/null").second.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown options in \"ps\"", _INPUT_ILLEGAL_);
    }
    command += " 2>/dev/null | " + command_grep + " 2> /dev/null";
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;
    }
    const string output = output_syscall = aurostd::execute2string(command);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " ps/grep output:" << endl << output << endl;
    }
    aurostd::string2vectorstring(output, vlines);
    for (i = 0; i < vlines.size(); i++) {
      aurostd::string2tokens(vlines[i], vtokens, " ");
      if (vtokens.size() < nopts) {
        continue;
      } // set by ps_opts
      const string& pid = vtokens[2]; // set by ps_opts
      if (vtokens[1] == pgid) { // set by ps_opts
        string proc = vtokens[nopts - 1]; // set by ps_opts
        for (j = nopts; j < vtokens.size(); j++) {
          proc += " " + vtokens[j];
        } // set by ps_opts
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " proc[i=" << i << "]=\"" << proc << "\"" << endl;
        }
        if (proc.find(process) == string::npos) {
          continue;
        }
        if (proc.find(command) != string::npos) {
          continue;
        }
        if (proc == command_grep) {
          continue;
        }
        vpids.push_back(pid);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vpids=" << aurostd::joinWDelimiter(vpids, ",") << endl;
      cerr << __AFLOW_FUNC__ << " vpids.empty()=" << vpids.empty() << endl;
    }
    return vpids;
  }

  // ***************************************************************************
  // Function ProcessRunning
  // ***************************************************************************
  // CO20210315
  bool ProcessRunning(const string& process, bool user_specific) {
    return !aurostd::ProcessPIDs(process, user_specific).empty();
  } // CO20210315

  // SD20220329 - overload to allow for only getting the PIDs with a specific PGID
  bool ProcessRunning(const string& process, const string& pgid, bool user_specific) {
    string output_syscall;
    return !aurostd::ProcessPIDs(process, pgid, output_syscall, user_specific).empty();
  }

  // ***************************************************************************
  // Function ProcessKill
  // ***************************************************************************
  // CO20210315
  // SD20220627 - Typicall signals that we use to kill processes are 9 (SIGKILL), 15 (SIGTERM) and 5 (SIGTRAP)
  // Do not throw an error since only SIGKILL and SIGTERM actually kill the process
  void ProcessKill(const string& process, bool user_specific, uint signal) { // CO20210315
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (signal < 1 || signal > 64) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "invalid signal specification", _VALUE_ILLEGAL_);
    }
    const vector<string> vpids = aurostd::ProcessPIDs(process, user_specific);
    if (vpids.empty()) {
      return;
    }
    string command = "kill";
    //[CO20210315 - does not work, user-specific comes from PID search]if(user_specific && !XHOST.user.empty()){command+=" -u "+XHOST.user;}
    command += " -" + aurostd::utype2string(signal) + " " + aurostd::joinWDelimiter(vpids, " ") + " 2>/dev/null";
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;
    }
    aurostd::execute(command);
  }

  // SD20220329 - overload to allow for only killing the PIDs with a specific PGID
  void ProcessKill(const string& process, const string& pgid, bool user_specific, uint signal) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (signal < 1 || signal > 64) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "invalid signal specification", _VALUE_ILLEGAL_);
    }
    if (pgid.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "PGID is empty", _INPUT_ILLEGAL_);
    }
    string output_syscall;
    const vector<string> vpids = aurostd::ProcessPIDs(process, pgid, output_syscall, user_specific);
    if (vpids.empty()) {
      return;
    }
    string command = "kill";
    command += " -" + aurostd::utype2string(signal) + " " + aurostd::joinWDelimiter(vpids, " ") + " 2>/dev/null";
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " running command=\"" << command << "\"" << endl;
    }
    aurostd::execute(command);
  }

  // ***************************************************************************
  // Function ProcessRenice
  // ***************************************************************************
  // CO20210315
  bool ProcessRenice(const string& process, int nvalue, bool user_specific, const string& pgid) { // CO20210315 //CO20221029 - void->bool to catch if renice worked and added pgid to renice executable-specific processes
    vector<string> vpids;
    if (pgid.empty()) {
      vpids = ProcessPIDs(process, user_specific);
    } // CO20221028
    else {
      string output_syscall;
      vpids = ProcessPIDs(process, pgid, output_syscall, user_specific);
    } // CO20221028
    if (vpids.empty()) {
      return false;
    }
    const string command = "renice " + aurostd::utype2string(nvalue) + " " + aurostd::joinWDelimiter(vpids, " ");
    const string err = aurostd::execute2OutErrPair(command).second;
    return err.empty();
  }

  // ***************************************************************************
  // Function ReniceAvailable
  // ***************************************************************************
  // CO20221029 - checks if renice up/down is allowed on the system
  bool ReniceAvailable() {
    aurostd::execute("sleep 10 &"); // to background
    if (!ProcessRenice("sleep", 19, true, aurostd::utype2string(getpgrp()))) {
      return false;
    } // try making more nice, this should almost always work
    if (!ProcessRenice("sleep", 0, true, aurostd::utype2string(getpgrp()))) {
      return false;
    } // try making default nice, this may not always work: https://superuser.com/questions/88542/why-cant-unix-users-renice-downwards
    return true;
  }

  bool GetMemoryUsagePercentage(double& usage_percentage_ram, double& usage_percentage_swap) { // CO20210601
    const bool LDEBUG = (false || XHOST.DEBUG);

    unsigned long long int free_ram = 0;
    unsigned long long int total_ram = 0;
    unsigned long long int free_swap = 0;
    unsigned long long int total_swap = 0;
    usage_percentage_ram = 0.0;
    usage_percentage_swap = 0.0;
    const bool memory_read = aurostd::GetMemory(free_ram, total_ram, free_swap, total_swap);
    if (memory_read) {
      if (total_ram > 0) {
        usage_percentage_ram = 100.0 * (((double) (total_ram - free_ram)) / ((double) (total_ram)));
      }
      if (total_swap > 0) {
        usage_percentage_swap = 100.0 * (((double) (total_swap - free_swap)) / ((double) (total_swap)));
      } // some qrats nodes have no swap
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [date=" << aflow_get_time_string() << "]" << endl; // helps debugging
        cerr << __AFLOW_FUNC__ << " free_ram=" << free_ram << endl;
        cerr << __AFLOW_FUNC__ << " used_ram=" << total_ram - free_ram << endl;
        cerr << __AFLOW_FUNC__ << " total_ram=" << total_ram << endl;
        cerr << __AFLOW_FUNC__ << " usage_percentage_ram=" << usage_percentage_ram << endl;
        cerr << __AFLOW_FUNC__ << " free_swap=" << free_swap << endl;
        cerr << __AFLOW_FUNC__ << " used_swap=" << total_swap - free_swap << endl;
        cerr << __AFLOW_FUNC__ << " total_swap=" << total_swap << endl;
        cerr << __AFLOW_FUNC__ << " usage_percentage_swap=" << usage_percentage_swap << endl;
        cerr << endl; // helps debugging
      }
    } else {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " unable to query memory on the node" << endl;
      }
    }
    return memory_read;
  }

  bool GetMemory(unsigned long long int& free_ram, unsigned long long int& total_ram, unsigned long long int& free_swap, unsigned long long int& total_swap) { // CO20210315 - only works for linux: needs `free` command
    // https://www.howtogeek.com/456943/how-to-use-the-free-command-on-linux/
    // will grab the total and the free
    // the free is the memory unused by anything
    // used column includes buff/cache, some of which the kernel can sacrifice for other applications if necessary
    // available column is an "estimate" of what could become available if needed
    // it's best to make decisions based on the free column
    // https://unix.stackexchange.com/questions/14102/real-memory-usage
    // free will follow real memory (physical RAM), using the available column will follow the actual memory (what could become available if necessary)
    const bool LDEBUG = (false || XHOST.DEBUG);

    if (!aurostd::IsCommandAvailable("free")) {
      return false;
    }
    const string output = aurostd::execute2string("free");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " free output:" << endl << output << endl;
    }
    // on most linux machines:
    //                             total        used        free      shared  buff/cache   available
    //               Mem:      395654628    41363940    30948848     4106252   323341840   349143640
    //               Swap:       2097148           0     2097148
    // on qrats:
    //                           total       used       free     shared    buffers     cached
    //              Mem:     264523076  255134588    9388488         36    1745836  232705576
    //              -/+ buffers/cache:   20683176  243839900
    //              Swap:      4194300      71436    4122864
    vector<string> vlines;
    aurostd::string2vectorstring(output, vlines);
    if (vlines.size() < 3) {
      return false;
    }
    vector<string> vtokens;
    uint iline = 0;
    // ram
    iline = 1;
    if (vlines[iline].find("Mem:") == string::npos) {
      return false;
    }
    aurostd::string2tokens(vlines[iline], vtokens, " ");
    if (!aurostd::isfloat(vtokens[1])) {
      return false;
    }
    total_ram = aurostd::string2utype<unsigned long long int>(vtokens[1]);
    if (!aurostd::isfloat(vtokens[3])) {
      return false;
    }
    free_ram = aurostd::string2utype<unsigned long long int>(vtokens[3]);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " free_ram=" << free_ram << " total_ram=" << total_ram << endl;
    }
    // swap
    iline = 2;
    if (vlines[iline].find("Swap:") == string::npos) {
      iline++;
    } // try next line
    if (vlines[iline].find("Swap:") == string::npos) {
      return false;
    }
    aurostd::string2tokens(vlines[iline], vtokens, " ");
    if (!aurostd::isfloat(vtokens[1])) {
      return false;
    }
    total_swap = aurostd::string2utype<unsigned long long int>(vtokens[1]);
    if (!aurostd::isfloat(vtokens[3])) {
      return false;
    }
    free_swap = aurostd::string2utype<unsigned long long int>(vtokens[3]);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " free_swap=" << free_swap << " total_swap=" << total_swap << endl;
    }
    //
    return true;
  }

  // ***************************************************************************
  // Function IsCommandAvailable
  // ***************************************************************************
  // tells you if the command is available
  bool IsCommandAvailable(const string& command, string& position) {
    // position=aurostd::execute2string("which "+command+" 2>&1 2> /dev/null");
    position = aurostd::execute2string("bash -c \"which " + command + " 2> /dev/null\""); // CO20210315 - put stderr to /dev/null //2>&1
    // cerr << position.length() << endl;
    aurostd::StringSubstInPlace(position, "\n", "");
    position = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(position); // CO20210315 - remove white spaces
    if (!position.empty()) {
      return true;
    }
    if (aurostd::FileExist("./" + command)) {
      position = "./" + command;
      return true;
    }
    if (aurostd::FileExist("/bin/" + command)) {
      position = "/bin/" + command;
      return true;
    }
    if (aurostd::FileExist("/sbin/" + command)) {
      position = "/sbin/" + command;
      return true;
    } // go around path
    if (aurostd::FileExist("/usr/bin/" + command)) {
      position = "/usr/bin/" + command;
      return true;
    } // go around path
    if (aurostd::FileExist("/usr/sbin/" + command)) {
      position = "/usr/sbin/" + command;
      return true;
    } // go around path
    if (aurostd::FileExist("/usr/local/bin/" + command)) {
      position = "/usr/local/bin/" + command;
      return true;
    } // go around path
    if (aurostd::FileExist("/usr/local/sbin/" + command)) {
      position = "/usr/local/sbin/" + command;
      return true;
    } // go around path
    if (aurostd::FileExist("/usr/local/maui/bin/" + command)) {
      position = "/usr/local/maui/bin/" + command;
      return true;
    } // go around path  //CO20200526
    position = "";
    return false;
  }

  bool IsCommandAvailable(const string& command) {
    string position;
    return aurostd::IsCommandAvailable(command, position);
  }

  // CO20180706 - fixed this function, previously command/position trampled all over each other
  bool IsCommandAvailableModify(string& command) {
    string position;
    if (!aurostd::IsCommandAvailable(command, position)) {
      return false;
    }
    command = position;
    return true;
  }

  // ***************************************************************************
  // Function CommandRequired
  // ***************************************************************************
  // tells you if the command is available
  bool CommandRequired(const string& command, string& position) {
    position = aurostd::execute2string("which " + command);
    aurostd::StringSubstInPlace(position, "\n", "");
    if (!position.empty()) {
      return true;
    }
    const string message = "\"" + command + "\" is not available";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    return false;
  }

  bool CommandRequired(const string& command) {
    string position;
    return CommandRequired(command, position);
  }

  // ***************************************************************************
  // Function IsExecutableAvailable
  // ***************************************************************************
  // tells you if the executable is available
  bool IsExecutableAvailable(const string& executable, string& position) {
    return IsCommandAvailable(executable, position);
  }

  bool IsExecutableAvailable(const string& executable) {
    string position;
    return IsCommandAvailable(executable, position);
  }

  // ***************************************************************************
  // Function ExecutableRequired
  // ***************************************************************************
  // tells you if the executable is available
  bool ExecutableRequired(const string& executable, string& position) {
    return CommandRequired(executable, position);
  }

  bool ExecutableRequired(const string& executable) {
    string position;
    return CommandRequired(executable, position);
  }

  // ***************************************************************************
  // DeleteOstringStreams
  // ***************************************************************************
  void StringstreamClean(ostringstream& aus) {
    aus.str(std::string()); // CO20200624
    aus.clear(); // CO20200624
  }
  void StringstreamClean(stringstream& aus) {
    aus.str(std::string()); // CO20200624
    aus.clear(); // CO20200624
  }

  // ***************************************************************************
  // Function FindIfStringInStream
  // ***************************************************************************
  //  This function returns true if string is in stream
  //  (on one line), or otherwise false. The search starts
  //  at the present file pointer location.  Note that this
  //  does alter the input string, resetting the file pointer
  //  to the input value at the end.
  // Dane Morgan style

  int FindIfStringInStream(const string& key, std::istream& instream) {
    // Get file pointer location at entry
    const int loc = instream.tellg();
    int found_match = 0;
    string s;
    getline(instream, s);
    int cont = 0;
    if (getline(instream, s)) {
      cont = 1;
    }
    while (cont) {
      const int id = s.find(key);
      if (id != (int) s.npos) { // Found key
        cont = 0;
        found_match = 1;
      }
      if (!getline(instream, s)) {
        cont = 0;
      }
    }
    // Clear any fail bits associated with searching.
    instream.clear();
    // Set file pointer location to entry value
    instream.seekg(loc);
    return found_match;
  }

  // ***************************************************************************
  // Print Messages Errors and Warnings on and off streams.
  // ***************************************************************************
#define ErrorBarString "EEEEE  ---------------------------------------------------------------------------------------------------------------------------- "
#define WarningBarString "WWWWW  ---------------------------------------------------------------------------------------------------------------------------- "

  void PrintANSIEscapeSequence(const aurostd::xoption& color, FILE* fstr) {
    if (color.option == false) {
      return;
    }
    if (color.flag("COLOR==GREEN")) {
      cursor_fore_green(fstr);
      return;
    }
    if (color.flag("COLOR==CYAN")) {
      cursor_fore_cyan(fstr);
      return;
    }
    if (color.flag("COLOR==YELLOW")) {
      cursor_fore_yellow(fstr);
      return;
    }
    if (color.flag("COLOR==RED")) {
      cursor_fore_red(fstr);
      return;
    }
  }

  void PrintMessageStream(ostringstream& stream, bool quiet, std::ostream& oss) {
    ofstream FileMESSAGE;
    return PrintMessageStream(FileMESSAGE, stream, quiet, oss);
  } // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet, std::ostream& oss) {
    const bool osswrite = true;
    return PrintMessageStream(FileMESSAGE, stream, quiet, osswrite, oss);
  } // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet, bool osswrite, std::ostream& oss) {
    // CO20181226 - split by newlines and print separately
    vector<string> message_parts;
    vector<string> _message_parts;
    const string stream_str = stream.str();
    aurostd::StringstreamClean(stream);
    aurostd::string2vectorstring(stream_str, _message_parts);
    for (size_t i = 0; i < _message_parts.size(); i++) {
      if (!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()) {
        message_parts.push_back(_message_parts[i]);
      }
    }
    if (message_parts.empty()) {
      return;
    }

    // CO20220129 - input quiet is always XHOST.QUIET, this is legacy: aflow.h (XHOST) was previously not made available in aurostd
    // keep redundant for now, just in case in the future we change our mind about including aflow.h
    // CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    // but it would require a full overhaul of many critical aflow printing functions
    // better not to touch and leave the overloads
    // bool verbose=(!XHOST.QUIET && !quiet && osswrite);  //ME20220503 - XHOST.QUIET should be part of quiet to allow for whitelisting
    bool verbose = (!quiet && osswrite);
    if (XHOST.QUIET_GLOBAL) {
      verbose = false;
    } // CO20220630
    if ((&oss == &cout) && XHOST.QUIET_COUT) {
      verbose = false;
    }
    if ((&oss == &cerr) && XHOST.QUIET_CERR) {
      verbose = false;
    }
    bool fancy_print = (!XHOST.vflag_control.flag("WWW") && !XHOST.vflag_control.flag("NO_FANCY_PRINT")); // CO20200404 - new web flag

    FILE* fstr = stdout;
    if (&oss == &std::cerr) {
      fstr = stderr;
    }

    // COLOR CANNOT BE A STRING, this construction will cause errors for the compiler
    // string color="\033[32m";
    // printf(color.c_str());
    // the compiler needs to verify that you are not printf'ing junk
    // so it needs to be a direct injection of code that the compiler can check
    // reference aurostd.h: CO20200624 START - adding from Jahnatek
    for (size_t i = 0; i < message_parts.size(); i++) {
      FileMESSAGE << message_parts[i] << endl;
    } // flush included in endl
    if (verbose) {
      string::size_type loc;
      string str2search; // replicate old behavior, look for ERROR coming from logger() which has two pre spaces
      aurostd::xoption color;
      color.clear(); // use xoption: .option is global color flag (do we have color?), and .vxscheme tells me which color
      if (fancy_print) {
        const string message = stream_str;
        // COMPLETE - START
        if (color.option == false) {
          str2search = "  COMPLETE "; // replicate old behavior, look for ERROR coming from logger() which has two pre spaces
          if (message.find(str2search) != string::npos) {
            color.option = true;
            color.flag("COLOR==GREEN", true);
          } // green
        }
        // COMPLETE - END
        // NOTICE - START
        if (color.option == false) {
          str2search = "  NOTICE "; // replicate old behavior, look for ERROR coming from logger() which has two pre spaces
          if (message.find(str2search) != string::npos) {
            color.option = true;
            color.flag("COLOR==CYAN", true);
          } // cyan
        }
        // NOTICE - END
      }

      // cursor_fore_green(fstr)
      // cursor_fore_cyan(fstr)

      if (color.option == false) {
        fancy_print = false;
      } // add others as needed
      if (fancy_print) {
        PrintANSIEscapeSequence(color, fstr);
      }
      for (size_t i = 0; i < message_parts.size(); i++) {
        loc = (!str2search.empty() ? message_parts[i].find(str2search) : string::npos);
        oss << message_parts[i].substr(0, loc);
        if (loc != string::npos) {
          // colors see here: https://en.m.wikipedia.org/wiki/ANSI_escape_code
          if (fancy_print) {
            cursor_attr_none(fstr); // turn off all cursor attributes
          }
          if (fancy_print) {
            PrintANSIEscapeSequence(color, fstr); // color
          }
          if (fancy_print) {
            cursor_attr_blink(fstr);
            cursor_attr_bold(fstr);
          } // bold+blink
          oss << str2search;
          if (fancy_print) {
            cursor_attr_none(fstr); // turn off all cursor attributes
          }
          if (fancy_print) {
            PrintANSIEscapeSequence(color, fstr); // color
          }
          oss << message_parts[i].substr(loc + str2search.size(), string::npos);
        }
        oss << endl; // flush included in endl
      }
      if (fancy_print) {
        cursor_attr_none(fstr); // turn off all cursor attributes
      }
    }
  }

  // CO20200624 - no std::ostream& oss input: THIS MUST GO TO CERR
  void PrintErrorStream(ostringstream& stream, bool quiet) {
    ofstream FileMESSAGE;
    return PrintErrorStream(FileMESSAGE, stream, quiet);
  } // CO20200624
  void PrintErrorStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet) {
    const bool osswrite = true;
    return PrintErrorStream(FileMESSAGE, stream, quiet, osswrite);
  } // CO20200624
  void PrintErrorStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet, bool osswrite) {
    // CO20181226 - split by newlines and print separately
    vector<string> message_parts;
    vector<string> _message_parts;
    const string stream_str = stream.str();
    aurostd::StringstreamClean(stream);
    aurostd::string2vectorstring(stream_str, _message_parts);
    for (size_t i = 0; i < _message_parts.size(); i++) {
      if (!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()) {
        message_parts.push_back(_message_parts[i]);
      }
    }
    if (message_parts.empty()) {
      return;
    }

    // CO20220129 - input quiet is always XHOST.QUIET, this is legacy: aflow.h (XHOST) was previously not made available in aurostd
    // keep redundant for now, just in case in the future we change our mind about including aflow.h
    // CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    // but it would require a full overhaul of many critical aflow printing functions
    // better not to touch and leave the overloads
    // bool verbose=(!XHOST.QUIET && !quiet && osswrite);  //[CO2010315 - not always, removing for OUTCARs read during vasp runs]verbose=true; //ALWAYS! //ME20220503 - XHOST.QUIET should be part of quiet to allow for whitelisting
    bool verbose = (!quiet && osswrite);
    if (XHOST.QUIET_GLOBAL) {
      verbose = false;
    } // CO20220630
    if (XHOST.QUIET_CERR) {
      verbose = false;
    }
    const bool fancy_print = (!XHOST.vflag_control.flag("WWW") && !XHOST.vflag_control.flag("NO_FANCY_PRINT")); // CO20200404 - new web flag

    FILE* fstr = stderr;

    FileMESSAGE << ErrorBarString << endl;
    for (size_t i = 0; i < message_parts.size(); i++) {
      FileMESSAGE << message_parts[i] << endl;
    } // flush included in endl
    FileMESSAGE << ErrorBarString << endl;
    if (verbose) {
      string::size_type loc;
      const string str2search = "  ERROR "; // replicate old behavior, look for ERROR coming from logger() which has two pre spaces
      std::ostream& oss = std::cerr;
      if (fancy_print) {
        cursor_fore_red(fstr); // red
      }
      oss << ErrorBarString << endl; // flush included in endl
      for (size_t i = 0; i < message_parts.size(); i++) {
        loc = message_parts[i].find(str2search);
        oss << message_parts[i].substr(0, loc);
        if (loc != string::npos) {
          if (fancy_print) {
            cursor_attr_none(fstr); // turn off all cursor attributes
          }
          if (fancy_print) {
            cursor_fore_red(fstr); // red
          }
          if (fancy_print) {
            cursor_attr_blink(fstr);
            cursor_attr_bold(fstr);
          } // bold+blink
          oss << str2search;
          if (fancy_print) {
            cursor_attr_none(fstr); // turn off all cursor attributes
          }
          if (fancy_print) {
            cursor_fore_red(fstr); // red
          }
          oss << message_parts[i].substr(loc + str2search.size(), string::npos);
        }
        oss << endl; // flush included in endl
      }
      oss << ErrorBarString << endl; // flush included in endl
      if (fancy_print) {
        cursor_attr_none(fstr); // turn off all cursor attributes
      }
    }
  }

  // CO20200624 - no std::ostream& oss input: THIS MUST GO TO CERR
  void PrintWarningStream(ostringstream& stream, bool quiet) {
    ofstream FileMESSAGE;
    return PrintWarningStream(FileMESSAGE, stream, quiet);
  } // CO20200624
  void PrintWarningStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet) {
    const bool osswrite = true;
    return PrintWarningStream(FileMESSAGE, stream, quiet, osswrite);
  } // CO20200624
  void PrintWarningStream(ofstream& FileMESSAGE, ostringstream& stream, bool quiet, bool osswrite) {
    // CO20181226 - split by newlines and print separately
    vector<string> message_parts;
    vector<string> _message_parts;
    const string stream_str = stream.str();
    aurostd::StringstreamClean(stream);
    aurostd::string2vectorstring(stream_str, _message_parts);
    for (size_t i = 0; i < _message_parts.size(); i++) {
      if (!aurostd::RemoveWhiteSpacesFromTheBack(_message_parts[i]).empty()) {
        message_parts.push_back(_message_parts[i]);
      }
    }
    if (message_parts.empty()) {
      return;
    }

    // CO20220129 - input quiet is always XHOST.QUIET, this is legacy: aflow.h (XHOST) was previously not made available in aurostd
    // keep redundant for now, just in case in the future we change our mind about including aflow.h
    // CO20220630 - note about osswrite, it is redundant with quiet, so it would be nice to get rid of it in the future
    // but it would require a full overhaul of many critical aflow printing functions
    // better not to touch and leave the overloads
    // bool verbose=(!XHOST.QUIET && !quiet && osswrite);  //[CO2010315 - not always, removing for OUTCARs read during vasp runs]verbose=true; //ALWAYS! //ME20220503 - XHOST.QUIET should be part of quiet to allow for whitelisting
    bool verbose = (!quiet && osswrite);
    if (XHOST.QUIET_GLOBAL) {
      verbose = false;
    } // CO20220630
    if (XHOST.QUIET_CERR) {
      verbose = false;
    }
    const bool fancy_print = (!XHOST.vflag_control.flag("WWW") && !XHOST.vflag_control.flag("NO_FANCY_PRINT")); // CO20200404 - new web flag

    FILE* fstr = stderr;

    FileMESSAGE << WarningBarString << endl;
    for (size_t i = 0; i < message_parts.size(); i++) {
      FileMESSAGE << message_parts[i] << endl;
    } // flush included in endl
    FileMESSAGE << WarningBarString << endl;
    if (verbose) {
      string::size_type loc;
      const string str2search = "  WARNING "; // replicate old behavior, look for WARNING coming from logger() which has two pre spaces
      std::ostream& oss = std::cerr;
      if (fancy_print) {
        cursor_fore_yellow(fstr); // yellow
      }
      oss << WarningBarString << endl; // flush included in endl
      for (size_t i = 0; i < message_parts.size(); i++) {
        loc = message_parts[i].find(str2search);
        oss << message_parts[i].substr(0, loc);
        if (loc != string::npos) {
          if (fancy_print) {
            cursor_attr_none(fstr); // turn off all cursor attributes
          }
          if (fancy_print) {
            cursor_fore_yellow(fstr); // yellow
          }
          if (fancy_print) {
            cursor_attr_blink(fstr);
            cursor_attr_bold(fstr);
          } // bold+blink
          oss << str2search;
          if (fancy_print) {
            cursor_attr_none(fstr); // turn off all cursor attributes
          }
          if (fancy_print) {
            cursor_fore_yellow(fstr); // yellow
          }
          oss << message_parts[i].substr(loc + str2search.size(), string::npos);
        }
        oss << endl; // flush included in endl
      }
      oss << WarningBarString << endl; // flush included in endl
      if (fancy_print) {
        cursor_attr_none(fstr); // turn off all cursor attributes
      }
    }
  }

  void PrintMessageStream(stringstream& stream, bool quiet, std::ostream& oss) {
    ofstream FileMESSAGE;
    return PrintMessageStream(FileMESSAGE, stream, quiet, oss);
  } // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet, std::ostream& oss) {
    const bool osswrite = true;
    return PrintMessageStream(FileMESSAGE, stream, quiet, osswrite, oss);
  } // CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet, bool osswrite, std::ostream& oss) {
    ostringstream omess;
    omess << stream.str();
    aurostd::StringstreamClean(stream);
    return PrintMessageStream(FileMESSAGE, omess, quiet, osswrite, oss);
  }

  void PrintErrorStream(stringstream& stream, bool quiet) {
    ofstream FileMESSAGE;
    return PrintErrorStream(FileMESSAGE, stream, quiet);
  } // CO20200624
  void PrintErrorStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet) {
    const bool osswrite = true;
    return PrintErrorStream(FileMESSAGE, stream, quiet, osswrite);
  } // CO20200624
  void PrintErrorStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet, bool osswrite) {
    ostringstream omess;
    omess << stream.str();
    aurostd::StringstreamClean(stream);
    return PrintErrorStream(FileMESSAGE, omess, quiet, osswrite);
  }

  void PrintWarningStream(stringstream& stream, bool quiet) {
    ofstream FileMESSAGE;
    return PrintWarningStream(FileMESSAGE, stream, quiet);
  } // CO20200624
  void PrintWarningStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet) {
    const bool osswrite = true;
    return PrintWarningStream(FileMESSAGE, stream, quiet, osswrite);
  } // CO20200624
  void PrintWarningStream(ofstream& FileMESSAGE, stringstream& stream, bool quiet, bool osswrite) {
    ostringstream omess;
    omess << stream.str();
    aurostd::StringstreamClean(stream);
    return PrintWarningStream(FileMESSAGE, omess, quiet, osswrite);
  }

  // ***************************************************************************
  // Execute Streams/Strings/C_strings
  // ***************************************************************************
  bool execute(ostringstream& command) {
    // cerr << "COMMAND " <<  command.str().c_str() << endl;
    execute(command.str()); // CO20200624
    aurostd::StringstreamClean(command);
    return true;
  }

  bool execute(stringstream& command) {
    // cerr << "COMMAND " <<  command.str().c_str() << endl;
    execute(command.str()); // CO20200624
    aurostd::StringstreamClean(command);
    return true;
  }

  bool execute(const string& command_raw) {
#ifdef AFLOW_MULTITHREADS_ENABLE // CO+HE20221116
    std::lock_guard<std::mutex> const lk(xthread_execute); // prevents race conditions likely caused by system calls
#endif
    bool const LDEBUG = (false || XHOST.DEBUG);
    // cerr << "COMMAND " <<  command.c_str() << endl;
    const string command = aurostd::CleanCommand4Execute(command_raw); // CO20200624
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " command.c_str()=\"" << command.c_str() << "\"" << endl;
    }
    const int rvalue = system(command.c_str());
    return rvalue == 0;
  }

  /// @brief execute a command that is itself thread safe and does not need a lock guard
  /// @param command_raw command to run
  /// @return true if command ran successfully else false
  /// @note this code replaces std::system as it has in modern implementation (UNIX03) a mutex lock
  /// @note the command is executed in a `/bin/sh` environment
  ///
  /// @xlink{FreeBSD/Apple LIBC, https://github.com/apple-oss-distributions/Libc/blob/main/stdlib/FreeBSD/system.c}
  /// @xlink{GNU LIBC, https://github.com/bminor/glibc/blob/master/sysdeps/posix/system.c}
  /// @authors
  /// @mod{HE,20240219,created}
  bool execute_thread_safe(const string& command_raw) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    const string command = aurostd::CleanCommand4Execute(command_raw);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " command.c_str()=\"" << command.c_str() << "\"" << endl;
    }
    pid_t pid;
    int status;

    const char* argv[] = {"/bin/sh", "-c", command.data(), nullptr};
    if (posix_spawn(&pid, argv[0], nullptr, nullptr, const_cast<char**>(argv), environ) != 0) {
      return false;
    }
    if (waitpid(pid, &status, 0) != pid) {
      return false;
    }
    if (status != 0) {
      return false;
    }
    return true;
  }

#ifdef _stringcharstar_
  bool execute(char* command_raw) {
    // cerr << "COMMAND " <<  command << endl;
    string command = std::string(command_raw); // CO20200624
    execute(command);
    return true;
  }
#endif

  // ***************************************************************************
  // Execute vectors/deque of Strings
  // ***************************************************************************
  bool execute(const deque<string>& vcommand) {
    for (size_t i = 0; i < vcommand.size(); i++) {
      execute(vcommand[i]);
    }
    return true;
  }
  bool execute(const vector<string>& vcommand) {
    for (size_t i = 0; i < vcommand.size(); i++) {
      execute(vcommand[i]);
    }
    return true;
  }

  /// @brief execute a command and returns the results as strings
  /// @param command_raw command to run
  /// @return std_out and std_err as string pair
  /// @note this code replaces the old approach using temporary files and `std::system`
  /// @note the command is executed in a `/bin/sh` environment
  ///
  /// @xlink{Pipe capacity, https://man7.org/linux/man-pages/man7/pipe.7.html}
  /// @authors
  /// @mod{HE,20240221,created}
  std::pair<std::string, std::string> execute2OutErrPair(const string& command_raw) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    int exit_code;
    pid_t pid;
    int cout_pipe[2];
    int cerr_pipe[2];
    posix_spawn_file_actions_t action;

    // Create two pipes to collect data
    if (pipe(cout_pipe) || pipe(cerr_pipe)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Could not create pipes for redirect", _RUNTIME_ERROR_);
    }

    // Connect the pipes
    posix_spawn_file_actions_init(&action);
    posix_spawn_file_actions_addclose(&action, cout_pipe[0]);
    posix_spawn_file_actions_addclose(&action, cerr_pipe[0]);
    posix_spawn_file_actions_adddup2(&action, cout_pipe[1], 1);
    posix_spawn_file_actions_adddup2(&action, cerr_pipe[1], 2);
    posix_spawn_file_actions_addclose(&action, cout_pipe[1]);
    posix_spawn_file_actions_addclose(&action, cerr_pipe[1]);

    const string command = aurostd::CleanCommand4Execute(command_raw);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cmdstream=\"" << command << "\"" << endl;
    }

    // The command is executed by `/bin/sh -c` to mirror the behavior of std::system()
    const char* argv[] = {"/bin/sh", "-c", command.data(), nullptr};

    // Run the command and check if the spawn process was successful
    if (posix_spawnp(&pid, argv[0], &action, nullptr, const_cast<char**>(argv), environ) != 0) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " posix_spawnp failed with error: " << std::strerror(errno) << endl;
      }
    }

    // The child-side of both pipes can now be closed
    close(cout_pipe[1]);
    close(cerr_pipe[1]);

    // Buffer and target to read pipe data
    string buffer(64 * 1024, ' '); // read in 64 kbytes chunks (should be the pipe capacity on modern systems - 16 pages)
    string std_out_content;
    string std_err_content;

    // Both pipes need to be read concurrently to avoid blocking the called program due to a full pipe
    // Set the pipes to non-blocking while reading
    fcntl(cout_pipe[0], F_SETFL, fcntl(cout_pipe[0], F_GETFL) | O_NONBLOCK);
    fcntl(cerr_pipe[0], F_SETFL, fcntl(cout_pipe[0], F_GETFL) | O_NONBLOCK);

    // Read all data until both pipes are closed by the called program
    ssize_t read_count;
    bool cout_pipe_active = true;
    bool cerr_pipe_active = true;
    while (cout_pipe_active | cerr_pipe_active) {
      if (cout_pipe_active) {
        // read count: -1 pipe has no data, 0 pipe closed
        read_count = read(cout_pipe[0], &buffer[0], buffer.length());
        if (read_count > 0) {
          std_out_content += buffer.substr(0, read_count);
        }
        if (read_count == 0) {
          cout_pipe_active = false;
        }
      }
      if (cerr_pipe_active) {
        read_count = read(cerr_pipe[0], &buffer[0], buffer.length());
        if (read_count > 0) {
          std_err_content += buffer.substr(0, read_count);
        }
        if (read_count == 0) {
          cerr_pipe_active = false;
        }
      }
    }

    // Wait until process is completely closed
    waitpid(pid, &exit_code, 0);

    if (exit_code != 0) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " command exited with code " << exit_code << endl;
      }
    }

    // Close all remaining connections
    close(cout_pipe[0]);
    close(cerr_pipe[0]);
    posix_spawn_file_actions_destroy(&action);

    return {std_out_content, std_err_content};
  }

  // ***************************************************************************
  // Execute & Report Streams/Strings/C_strings
  // ***************************************************************************
  string execute2string(const string& command) { // HE20240220 switched to execute2OutErrPair
    return execute2OutErrPair(command).first;
  }

  string execute2string(ostringstream& command) { // CO20200624 - added file system IO mode
    const string command_str = command.str();
    aurostd::StringstreamClean(command);
    return execute2string(command_str); // CO20200624
  }

  string execute2string(stringstream& command) { // CO20200624 - added file system IO mode
    const string command_str = command.str();
    aurostd::StringstreamClean(command);
    return execute2string(command_str); // CO20200624
  }

  vector<string> execute2string(const vector<string>& vcommand) { // CO20200624 - added file system IO mode
    vector<string> out;
    for (size_t i = 0; i < vcommand.size(); i++) {
      out.push_back(execute2string(vcommand[i])); // CO20200624
    }
    return out;
  }

  deque<string> execute2string(const deque<string>& vcommand) { // CO20200624 - added file system IO mode
    deque<string> out;
    for (size_t i = 0; i < vcommand.size(); i++) {
      out.push_back(execute2string(vcommand[i])); // CO20200624
    }
    return out;
  }

#ifdef _stringcharstar_
  string execute2string(char* command) { // HE20240220 switched to execute2OutErrPair
    string command_str = string(command);
    return execute2OutErrPair(command_str).first; // CO20200624
  }
#endif

  string CleanCommand4Execute(const string& command_raw) { // CO20200624
    const bool LDEBUG = (false || XHOST.DEBUG);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " command(pre )=\"" << command_raw << "\"" << endl;
    }
    // CO20200624 START - some command cleanup
    vector<string> vtokens;
    vector<string> vtokens_new;
    aurostd::string2vectorstring(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(command_raw), vtokens);
    string tmp;
    uint i = 0;
    for (i = 0; i < vtokens.size(); i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vtokens[i=" << i << "](pre )=\"" << vtokens[i] << "\"" << endl;
      }
      tmp = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vtokens[i]);
      aurostd::CleanStringASCII_InPlace(tmp);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vtokens[i=" << i << "](post)=\"" << tmp << "\"" << endl;
      }
      if (!tmp.empty()) {
        vtokens_new.push_back(tmp);
      }
    }
    if (vtokens_new.empty()) {
      return "";
    }
    //[CO20210312 - must be smarter, could end with && or ;]string command=aurostd::joinWDelimiter(vtokens_new,"; ");
    string command;
    uint len = 0;
    bool add_semicolon = false;
    for (i = 0; i < vtokens_new.size(); i++) { // vtokens_new has no empty entries, so we don't need to check again
      const string& cmd = vtokens_new[i];
      command += cmd;
      if (i < vtokens_new.size() - 1) {
        add_semicolon = false;
        len = cmd.size();
        if (!(cmd[len - 1] == '&' || cmd[len - 1] == '|' || cmd[len - 1] == ';')) {
          add_semicolon = true;
        }
        if (add_semicolon) {
          command += "; ";
        } else {
          command += " ";
        } // add a space, looks good for ' this && that '
      }
    }
    // CO20200624 END - some command cleanup
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " command(post)=\"" << command << "\"" << endl;
    }
    return command;
  }

  // ***************************************************************************
  // Execute & Report Int Streams/Strings/C_strings
  // ***************************************************************************
  template <class utype> utype execute2utype(ostringstream& command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template <class utype> utype execute2utype(stringstream& command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template <class utype> utype execute2utype(string command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }

  template <class utype> vector<utype> execute2utype(vector<string> vcommand) {
    vector<utype> out;
    for (size_t i = 0; i < vcommand.size(); i++) {
      out.push_back((utype) execute2utype<utype>(vcommand[i]));
    }
    return out;
  }

  template <class utype> deque<utype> execute2utype(deque<string> vcommand) {
    deque<utype> out;
    for (size_t i = 0; i < vcommand.size(); i++) {
      out.push_back((utype) execute2utype<utype>(vcommand[i]));
    }
    return out;
  }
#define AST_TEMPLATE(atype) template atype execute2utype(ostringstream& command);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype execute2utype(stringstream& command);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype execute2utype(string command);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template vector<atype> execute2utype(vector<string> vcommand);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template deque<atype> execute2utype(deque<string> vcommand);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

#ifdef _stringcharstar_
  template <class utype> utype execute2utype(char* command) {
    return (utype) aurostd::string2utype<utype>(execute2string(command));
  }
#endif

  // ***************************************************************************
  // Sleep
  // ***************************************************************************
  unsigned int Sleep(unsigned int seconds) {
    //  ostringstream aus;
    // aus << "sleep " << (int) seconds << " " << endl;
    // aurostd::execute(aus);
    return sleep(seconds);
  }

  // *******************************************************************************************
  // *******************************************************************************************
  vector<string> GrepFile(const string& filename, const string& keyword, bool RemoveWS, bool RemoveComments) { // CO20210623 - update after integrating new substring2bool (RemoveWS,RemoveComments)
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    vector<string> vout;
    if (aurostd::FileExist(filename) == false) {
      return vout;
    }

    if (RemoveWS && RemoveComments) {
      ;
    } // keep busy until we update substring2bool

    string strline;
    ifstream FileStream;
    FileStream.open(filename.c_str(), std::ios::in);
    while (getline(FileStream, strline)) {
      if (aurostd::substring2bool(strline, keyword, RemoveWS)) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " found matching line: \"" << strline << "\"" << endl;
        }
        vout.push_back(strline);
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return vout;
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(ifstream& FileIN, stringstream& StringstreamOUTPUT, const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    aurostd::StringstreamClean(StringstreamOUTPUT);
    string strline;
    FileIN.clear();
    FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    bool status = false;
    while (getline(FileIN, strline)) {
      if (status) {
        StringstreamOUTPUT << strline << endl;
      }
      if (aurostd::substring2bool(strline, Keyword_start)) {
        status = true;
      }
    }
    FileIN.clear();
    FileIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return status; // return false if something got messed up
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(stringstream& StringStreamIN, stringstream& StringstreamOUTPUT, const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    const string StringIN = StringStreamIN.str();
    return ExtractJustAfterToStringstreamEXPLICIT(StringIN, StringstreamOUTPUT, Keyword_start);
  }

  bool ExtractJustAfterToStringstreamEXPLICIT(const string& StringIN, stringstream& StringstreamOUTPUT, const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    aurostd::StringstreamClean(StringstreamOUTPUT);
    bool status = false;
    vector<string> tokens;
    aurostd::string2tokens(StringIN, tokens, "\n");
    for (size_t i = 0; i < tokens.size(); i++) {
      if (status) {
        StringstreamOUTPUT << tokens[i] << endl;
      }
      if (aurostd::substring2bool(tokens[i], Keyword_start)) {
        status = true;
      }
    }
    return status; // return false if something got messed up
  }

  bool ExtractJustAfterToFileEXPLICIT(ifstream& FileIN, const string& FileNameOUTPUT, const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    stringstream StringstreamOUTPUT;
    if (ExtractJustAfterToStringstreamEXPLICIT(FileIN, StringstreamOUTPUT, Keyword_start)) {
      return aurostd::stringstream2file(StringstreamOUTPUT, FileNameOUTPUT);
    }
    return false;
  }

  bool ExtractJustAfterToStringEXPLICIT(ifstream& FileIN, string& StringOUTPUT, const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    stringstream StringstreamOUTPUT;
    if (ExtractJustAfterToStringstreamEXPLICIT(FileIN, StringstreamOUTPUT, Keyword_start)) {
      StringOUTPUT = StringstreamOUTPUT.str();
      return true;
    }
    return false;
  }

  bool ExtractJustAfterToStringEXPLICIT(const string& StringIN, string& StringOUTPUT, const string& Keyword_start) { // AFLOW_FUNCTION_IMPLEMENTATION
    stringstream StringstreamOUTPUT;
    if (ExtractJustAfterToStringstreamEXPLICIT(StringIN, StringstreamOUTPUT, Keyword_start)) {
      StringOUTPUT = StringstreamOUTPUT.str();
      return true;
    }
    return false;
  }

  // ***************************************************************************
  // Function stream2vectorstring return UINT
  // ***************************************************************************
  // take istream into a vector strings - Stefano Curtarolo

  ///
  size_t stream2vectorstring(std::istream& istreamIN, vector<string>& vstringout) {
    vstringout.clear();
    for (std::string line; std::getline(istreamIN, line);) {
      vstringout.push_back(line);
    }
    return vstringout.size();
  }
  // take ifstream into a vector strings - Stefano Curtarolo
  size_t stream2vectorstring(std::ifstream& ifstreamIN, vector<string>& vstringout) {
    vstringout.clear();
    for (std::string line; std::getline(ifstreamIN, line);) {
      vstringout.push_back(line);
    }
    return vstringout.size();
  }

  size_t stream2vectorstring(std::stringstream& stringstreamIN, vector<string>& vstringout) {
    vstringout.clear();
    for (std::string line; std::getline(stringstreamIN, line);) {
      vstringout.push_back(line);
    }
    return vstringout.size();
  }
  uint trimEmptyEdges(vector<string>& vstringout) { // CO20230502
    uint count = vstringout.size();
    // start with front
    while (!vstringout.empty()) {
      if (vstringout.front().find_first_not_of("\t\n ") != string::npos) { // CO20230502 - fast way to check for empty string //!aurostd::RemoveWhiteSpaces(vstringout.front()).empty())
        break;
      }
      vstringout.erase(vstringout.begin());
      count--;
    }
    // now back
    while (!vstringout.empty()) {
      if (vstringout.back().find_first_not_of("\t\n ") != string::npos) { // CO20230502 - fast way to check for empty string //aurostd::RemoveWhiteSpaces(vstringout.back()).empty())
        break;
      }
      vstringout.pop_back();
      count--;
    }
    return count;
  }
  size_t string2vectorstring(const string& stringIN, vector<string>& vstringout, bool consecutive, bool trim_edges) { // CO20170613
    // CO mods 20170613
    // we are adding functionality here, because string2tokens will treat "\n\n" same as "\n", but not "\n \n"
    // consecutive will do the following: "sssss" -> <"s","s",...>
    // trim_edges will remove delimiters from beginning and end, similar to consecutive=false behavior
    // return aurostd::string2tokens(stringIN,vstringout,"\n",true);
    aurostd::string2tokens(stringIN, vstringout, "\n", consecutive);
    if (trim_edges) {
      trimEmptyEdges(vstringout);
    }
    return vstringout.size();
  }

  // ***************************************************************************
  // Function string2vectorstring return VECTOR
  // ***************************************************************************
  // take sitring into a vector strings - Stefano Curtarolo
  vector<string> stream2vectorstring(std::istream& istreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(istreamIN, vstringout);
    return vstringout;
  }
  vector<string> stream2vectorstring(std::ifstream& iftreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(iftreamIN, vstringout);
    return vstringout;
  }
  vector<string> stream2vectorstring(std::stringstream& stringstreamIN) {
    vector<string> vstringout;
    aurostd::stream2vectorstring(stringstreamIN, vstringout);
    return vstringout;
  }
  vector<string> string2vectorstring(const string& stringIN, bool consecutive, bool trim_edges) { // CO20170613
    vector<string> vstringout;
    aurostd::string2vectorstring(stringIN, vstringout, consecutive, trim_edges); // CO20170613
    return vstringout;
  }

  // ***************************************************************************
  // Function liststring2string return string
  // ***************************************************************************
  string liststring2string(string s00,
                           string s01,
                           string s02,
                           string s03,
                           string s04,
                           string s05,
                           string s06,
                           string s07,
                           string s08,
                           string s09,
                           string s0A,
                           string s0B,
                           string s0C,
                           string s0D,
                           string s0E,
                           string s0F,
                           string s10,
                           string s11,
                           string s12,
                           string s13,
                           string s14,
                           string s15,
                           string s16,
                           string s17,
                           string s18,
                           string s19,
                           string s1A,
                           string s1B,
                           string s1C,
                           string s1D,
                           string s1E,
                           string s1F,
                           string s20,
                           string s21,
                           string s22,
                           string s23,
                           string s24,
                           string s25,
                           string s26,
                           string s27,
                           string s28,
                           string s29,
                           string s2A,
                           string s2B,
                           string s2C,
                           string s2D,
                           string s2E,
                           string s2F) {
    string out;
    if (!s00.empty()) {
      out += s00 + "\n";
    }
    if (!s01.empty()) {
      out += s01 + "\n";
    }
    if (!s02.empty()) {
      out += s02 + "\n";
    }
    if (!s03.empty()) {
      out += s03 + "\n";
    }
    if (!s04.empty()) {
      out += s04 + "\n";
    }
    if (!s05.empty()) {
      out += s05 + "\n";
    }
    if (!s06.empty()) {
      out += s06 + "\n";
    }
    if (!s07.empty()) {
      out += s07 + "\n";
    }
    if (!s08.empty()) {
      out += s08 + "\n";
    }
    if (!s09.empty()) {
      out += s09 + "\n";
    }
    if (!s0A.empty()) {
      out += s0A + "\n";
    }
    if (!s0B.empty()) {
      out += s0B + "\n";
    }
    if (!s0C.empty()) {
      out += s0C + "\n";
    }
    if (!s0D.empty()) {
      out += s0D + "\n";
    }
    if (!s0E.empty()) {
      out += s0E + "\n";
    }
    if (!s0F.empty()) {
      out += s0F + "\n";
    }
    if (!s10.empty()) {
      out += s10 + "\n";
    }
    if (!s11.empty()) {
      out += s11 + "\n";
    }
    if (!s12.empty()) {
      out += s12 + "\n";
    }
    if (!s13.empty()) {
      out += s13 + "\n";
    }
    if (!s14.empty()) {
      out += s14 + "\n";
    }
    if (!s15.empty()) {
      out += s15 + "\n";
    }
    if (!s16.empty()) {
      out += s16 + "\n";
    }
    if (!s17.empty()) {
      out += s17 + "\n";
    }
    if (!s18.empty()) {
      out += s18 + "\n";
    }
    if (!s19.empty()) {
      out += s19 + "\n";
    }
    if (!s1A.empty()) {
      out += s1A + "\n";
    }
    if (!s1B.empty()) {
      out += s1B + "\n";
    }
    if (!s1C.empty()) {
      out += s1C + "\n";
    }
    if (!s1D.empty()) {
      out += s1D + "\n";
    }
    if (!s1E.empty()) {
      out += s1E + "\n";
    }
    if (!s1F.empty()) {
      out += s1F + "\n";
    }
    if (!s20.empty()) {
      out += s20 + "\n";
    }
    if (!s21.empty()) {
      out += s21 + "\n";
    }
    if (!s22.empty()) {
      out += s22 + "\n";
    }
    if (!s23.empty()) {
      out += s23 + "\n";
    }
    if (!s24.empty()) {
      out += s24 + "\n";
    }
    if (!s25.empty()) {
      out += s25 + "\n";
    }
    if (!s26.empty()) {
      out += s26 + "\n";
    }
    if (!s27.empty()) {
      out += s27 + "\n";
    }
    if (!s28.empty()) {
      out += s28 + "\n";
    }
    if (!s29.empty()) {
      out += s29 + "\n";
    }
    if (!s2A.empty()) {
      out += s2A + "\n";
    }
    if (!s2B.empty()) {
      out += s2B + "\n";
    }
    if (!s2C.empty()) {
      out += s2C + "\n";
    }
    if (!s2D.empty()) {
      out += s2D + "\n";
    }
    if (!s2E.empty()) {
      out += s2E + "\n";
    }
    if (!s2F.empty()) {
      out += s2F + "\n";
    }
    return out;
  }

  // ***************************************************************************
  // Function stream2dequestring return UINT
  // ***************************************************************************
  // take istream into a deque strings - Stefano Curtarolo
  uint stream2dequestring(std::istream& istreamIN, deque<string>& vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while (!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp, _CIN_LINE_BUFFER_LENGTH_ - 1);
      vstringout.emplace_back(tmp);
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size(); // return false if something got messed up
  }
  // take ifstream into a deque strings - Stefano Curtarolo
  uint stream2dequestring(std::ifstream& ifstreamIN, deque<string>& vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while (!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp, _CIN_LINE_BUFFER_LENGTH_ - 1);
      vstringout.emplace_back(tmp);
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size(); // return false if something got messed up
  }
  uint stream2dequestring(std::stringstream& stringstreamIN, deque<string>& vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while (!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp, _CIN_LINE_BUFFER_LENGTH_ - 1);
      vstringout.emplace_back(tmp);
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.size(); // return false if something got messed up
  }
  uint string2dequestring(const string& stringIN, deque<string>& vstringout) {
    return aurostd::string2tokens(stringIN, vstringout, "\n");
  }

  // ***************************************************************************
  // Function string2dequestring return DEQUE
  // ***************************************************************************
  // take sitring into a deque strings - Stefano Curtarolo
  deque<string> stream2dequestring(std::istream& istreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(istreamIN, vstringout);
    return vstringout;
  }
  deque<string> stream2dequestring(std::ifstream& iftreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(iftreamIN, vstringout);
    return vstringout;
  }
  deque<string> stream2dequestring(std::stringstream& stringstreamIN) {
    deque<string> vstringout;
    aurostd::stream2dequestring(stringstreamIN, vstringout);
    return vstringout;
  }
  deque<string> string2dequestring(const string& stringIN) {
    deque<string> vstringout;
    aurostd::string2dequestring(stringIN, vstringout);
    return vstringout;
  }

  // ***************************************************************************
  // Function ostream2string  istream2string
  // ***************************************************************************
  // convert ostream/istream to string - Stefano Curtarolo
  std::string ostream2string(std::ostream& oss) {
    std::stringstream soss;
    soss << oss.rdbuf();
    return soss.str();
  }

  // ***************************************************************************
  // Function stream2string
  // ***************************************************************************
  // take istream into a  strings - Stefano Curtarolo
  uint stream2string(std::istream& istreamIN, string& vstringout) {
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while (!istreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      istreamIN.getline(tmp, _CIN_LINE_BUFFER_LENGTH_ - 1);
      vstringout += string(tmp) + "\n"; // ME20210206 - fixed line break
    }
    // istreamIN.clear(); // istreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length(); // return false if something got messed up
  }

  // take ifstream into a  strings - Stefano Curtarolo
  uint stream2string(std::ifstream& ifstreamIN, string& vstringout) {
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while (!ifstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      ifstreamIN.getline(tmp, _CIN_LINE_BUFFER_LENGTH_ - 1);
      vstringout += string(tmp) + "\n"; // ME20210206 - fixed line break
    }
    // ifstreamIN.clear(); // ifstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length(); // return false if something got messed up
  }

  uint stream2string(std::stringstream& stringstreamIN, string& vstringout) {
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    vstringout.clear();
    while (!stringstreamIN.eof()) {
      char tmp[_CIN_LINE_BUFFER_LENGTH_];
      stringstreamIN.getline(tmp, _CIN_LINE_BUFFER_LENGTH_ - 1);
      vstringout += string(tmp) + "\n"; // ME20210206 - fixed line break
    }
    // stringstreamIN.clear();stringstreamIN.seekg(0); // ******* INPUT FILE goes at the beginning
    return vstringout.length(); // return false if something got messed up
  }

  // ***************************************************************************
  // Function getenv2string getenv2int getenv2uint getenv2double
  // ***************************************************************************
  // convert environments to string;
  string getenv2string(const string& str) {
    if (getenv(str.c_str()) == nullptr) {
      return string("");
    }
    return string(getenv(str.c_str()));
  }
  int getenv2int(const string& str) {
    if (getenv(str.c_str()) == nullptr) {
      return 0;
    }
    return aurostd::string2utype<int>(getenv(str.c_str()));
  }
  uint getenv2uint(const string& str) {
    if (getenv(str.c_str()) == nullptr) {
      return uint(0);
    }
    return aurostd::string2utype<uint>(getenv(str.c_str()));
  }
  double getenv2double(const string& str) {
    if (getenv(str.c_str()) == nullptr) {
      return double(0);
    }
    return aurostd::string2utype<double>(getenv(str.c_str()));
  }

  // ***************************************************************************
  // Function Chmod
  // ***************************************************************************
  /// @brief change permissions of the path
  /// @param chmod path chmod
  /// @param path_raw path to object
  /// @param perm_opt_file path permission options
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @mod{SG,20240531,added check if path exists}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool Chmod(const uint chmod, const string& path_raw, const std::filesystem::perm_options perm_opt) {
    const std::filesystem::path path = CleanFileName(path_raw);
    if (path.empty() || !std::filesystem::exists(path)) {
      return false;
    }
    std::filesystem::permissions(path, static_cast<std::filesystem::perms>(chmod), perm_opt);
    return true;
  }

  // ***************************************************************************
  // Function ChmodRecursive
  // ***************************************************************************
  /// @brief change permissions of the directory's content recursively
  /// @param chmod_dir directory chmod
  /// @param chmod_file file chmod
  /// @param _directory path to directory
  /// @param perm_opt_dir directory permission options
  /// @param perm_opt_file file permission options
  /// @authors
  /// @mod{SD,20240312,created function}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool ChmodRecursive(const uint chmod_dir, const uint chmod_file, const string& _directory, const std::filesystem::perm_options perm_opt_dir, const std::filesystem::perm_options perm_opt_file) {
    const string directory = CleanFileName(_directory);
    if (directory.empty()) {
      return false;
    }
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::recursive_directory_iterator(directory)) {
      if (dir_entry.is_directory()) {
        std::filesystem::permissions(dir_entry.path(), static_cast<std::filesystem::perms>(chmod_dir), perm_opt_dir);
      } else if (dir_entry.is_regular_file()) {
        std::filesystem::permissions(dir_entry.path(), static_cast<std::filesystem::perms>(chmod_file), perm_opt_file);
      } else {
        continue;
      }
    }
    return true;
  }

  // #define DEBUG_STRING2TOKENS
  //  ***************************************************************************
  //  Function string2tokens string2tokens<utype>
  //  ***************************************************************************
  //  Finds string2tokens to split strings in tokens
  //  Stefano Curtarolo
  //  void string2tokens(const string& str,vector<string>& tokens,const string& delimiters=" ") {  //[CO20200106 - close bracket for indenting]}
  uint string2tokens(const string& str, std::vector<string>& tokens, const string& delimiters, bool consecutive) { // CO20170613
    // CO mods 20170613
    // we are adding functionality here, because string2tokens will treat "\n\n" same as "\n", but not "\n \n"
    // consecutive will do the following: "sssss" -> <"s","s",...>
    // consecutive ALSO starts at 0, not first_not_of
    // return aurostd::string2tokens(stringIN,vstringout,"\n",true);
    tokens.clear(); // clear in the case there was something already in!
    string::size_type lastPos = (consecutive ? 0 : str.find_first_not_of(delimiters, 0)); // Skip delimiters at beginning.
    string::size_type pos = str.find_first_of(delimiters, lastPos); // Find first "non-delimiter".
#ifdef DEBUG_STRING2TOKENS
    cerr << __AFLOW_FUNC__ << " delimiters=" << delimiters << endl;
    cerr << __AFLOW_FUNC__ << " consecutive=" << consecutive << endl;
    cerr << __AFLOW_FUNC__ << " lastPos=" << lastPos << endl;
    cerr << __AFLOW_FUNC__ << " pos=" << pos << endl;
#endif
    while (pos != string::npos || lastPos != string::npos) {
      tokens.push_back(str.substr(lastPos, pos - lastPos)); // Found a token, add it to the vector.
#ifdef DEBUG_STRING2TOKENS
      cerr << __AFLOW_FUNC__ << " tokens.back()=\"" << tokens.back() << "\"" << endl;
#endif
      if (consecutive) {
        lastPos = (pos != string::npos ? pos + 1 : string::npos);
      } else {
        lastPos = str.find_first_not_of(delimiters, pos);
      } // Skip delimiters.  Note the "not_of"
      pos = str.find_first_of(delimiters, lastPos); // Find next "non-delimiter"
#ifdef DEBUG_STRING2TOKENS
      cerr << __AFLOW_FUNC__ << " lastPos=" << lastPos << endl;
      cerr << __AFLOW_FUNC__ << " pos=" << pos << endl;
#endif
    }
    return tokens.size();
  }

  // void string2tokens(const string& str,deque<string>& tokens,const string& delimiters=" ") { //[CO20200106 - close bracket for indenting]}
  uint string2tokens(const string& str, std::deque<string>& tokens, const string& delimiters, bool consecutive) { // CO20170613
    vector<string> vtokens;
    uint i = aurostd::string2tokens(str, vtokens, delimiters, consecutive); // CO20170613
    tokens.clear();
    for (i = 0; i < vtokens.size(); i++) {
      tokens.push_back(vtokens[i]);
    }
    return tokens.size();
  }
  template <class utype> uint string2tokens(const string& str, std::vector<utype>& tokens, const string& delimiters, bool consecutive) { // CO20170613
    vector<string> stokens;
    const uint out = aurostd::string2tokens(str, stokens, delimiters, consecutive); // CO20170613
    tokens.clear();
    for (size_t i = 0; i < stokens.size(); i++) {
      tokens.push_back(aurostd::string2utype<utype>(stokens[i]));
    }
    return out;
  }
#define AST_TEMPLATE(utype) template uint string2tokens(const string&, std::vector<utype>&, const string&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> uint string2tokens(const string& str, std::deque<utype>& tokens, const string& delimiters, bool consecutive) { // CO20170613
    deque<string> stokens;
    const uint out = aurostd::string2tokens(str, stokens, delimiters, consecutive); // CO20170613
    tokens.clear();
    for (size_t i = 0; i < stokens.size(); i++) {
      tokens.push_back(aurostd::string2utype<utype>(stokens[i]));
    }
    return out;
  }
#define AST_TEMPLATE(utype) template uint string2tokens(const string&, std::deque<utype>&, const string&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  // SD20220504 - string2tokens, but using a single delimiter that can have more than one char. Makes use of "find" rather than "find_first_of"
  uint string2tokensByDelimiter(const string& str, vector<string>& tokens, const string& delimiter) { // SD20220504
    tokens.clear();
    const uint dlen = delimiter.length();
    if (dlen == 1) {
      return aurostd::string2tokens(str, tokens, delimiter);
    }
    string::size_type lpos = 0;
    string::size_type cpos = str.find(delimiter, lpos);
    while (cpos != string::npos) {
      tokens.push_back(str.substr(lpos, cpos - lpos));
      lpos = cpos + dlen;
      cpos = str.find(delimiter, lpos);
    }
    tokens.push_back(str.substr(lpos, str.length() - lpos));
    return tokens.size();
  }
  uint string2tokensByDelimiter(const string& str, deque<string>& tokens, const string& delimiter) { // SD20220504
    tokens.clear();
    vector<string> vtokens;
    aurostd::string2tokensByDelimiter(str, vtokens, delimiter);
    tokens = aurostd::vector2deque(vtokens);
    return tokens.size();
  }

  // ***************************************************************************
  // Function string2tokensAdd string2tokensAdd<utype>
  // ***************************************************************************

  uint string2tokensAdd(const string& str, std::vector<string>& tokens, const string& delimiters) {
    vector<string> vtokens;
    uint i = aurostd::string2tokens(str, vtokens, delimiters);
    for (i = 0; i < vtokens.size(); i++) {
      tokens.push_back(vtokens[i]);
    }
    return tokens.size();
  }
  uint string2tokensAdd(const string& str, std::deque<string>& tokens, const string& delimiters) {
    vector<string> vtokens;
    uint i = aurostd::string2tokens(str, vtokens, delimiters);
    for (i = 0; i < vtokens.size(); i++) {
      tokens.push_back(vtokens[i]);
    }
    return tokens.size();
  }
  template <class utype> uint string2tokensAdd(const string& str, std::vector<utype>& tokens, const string& delimiters) {
    vector<string> vtokens;
    uint i = aurostd::string2tokens(str, vtokens, delimiters);
    for (i = 0; i < vtokens.size(); i++) {
      tokens.push_back(aurostd::string2utype<utype>(vtokens[i]));
    }
    return tokens.size();
  }
  template <class utype> uint string2tokensAdd(const string& str, std::deque<utype>& tokens, const string& delimiters) {
    deque<string> vtokens;
    uint i = aurostd::string2tokens(str, vtokens, delimiters);
    for (i = 0; i < vtokens.size(); i++) {
      tokens.push_back(aurostd::string2utype<utype>(vtokens[i]));
    }
    return tokens.size();
  }

  // ***************************************************************************
  // Function stream2stream
  // ***************************************************************************
  // convert whatever into a string !
  template <typename typeTo, typename typeFrom> typeTo stream2stream(const typeFrom& from, int precision, char FORMAT) { // CO20210315 - cleaned up
    std::stringstream temp;
    if (FORMAT == DEFAULT_STREAM) {
      ;
    } // default
    if (FORMAT == FIXED_STREAM) {
      temp << std::fixed;
    }
    if (FORMAT == SCIENTIFIC_STREAM) {
      temp << std::scientific;
    }
    temp.precision(precision);
    temp << from;
    typeTo to = typeTo();
    temp >> to;
    return to;
  }
  template <typename typeTo, typename typeFrom> typeTo stream2stream(const typeFrom& from, int precision) { // CO20210315 - cleaned up
    return (typeTo) stream2stream<typeTo>(from, precision, DEFAULT_STREAM);
  }
  template <typename typeTo, typename typeFrom> typeTo stream2stream(const typeFrom& from) { // CO20210315 - cleaned up
    return (typeTo) stream2stream<typeTo>(from, AUROSTD_DEFAULT_PRECISION, DEFAULT_STREAM);
  }

  /// @brief convert a string to an utype (number)
  /// @param from string to convert
  /// @param base number base - default 10
  /// @return parsed number
  /// @authors
  /// @mod{HE,20260218,complete rewrite using C++99/C++11}
  /// @note failed conversions or empty strings returns a 0 \n
  /// @note only with base 10 decimals are possible \n
  /// @note all characters after a valid number are ignored\n
  /// @note prefixes (0x and 0X) are filtered when using base 16 \n
  template <typename utype> utype string2utype(const std::string& from, const uint base) {
    if (from.empty()) {
      return 0;
    }

    // static bool sets
    static const std::set<std::string, std::less<>> true_set = {"TRUE", "True", "true", "T", "t", ".TRUE.", ".True.", ".true.", "1"};
    static const std::set<std::string, std::less<>> false_set = {"FALSE", "False", "false", "F", "f", ".FALSE.", ".False.", ".false.", "0"};

    if (const auto it = true_set.find(from); it != true_set.end()) {
      return 1;
    }
    if (const auto it = false_set.find(from); it != false_set.end()) {
      return 0;
    }

    char* p_end{};
    // if the user wants integer (not bool) values we can use strtol directly
    if constexpr (std::is_integral_v<utype> and !std::is_same_v<utype, bool>) {
      return std::strtol(from.data(), &p_end, base);
    }
    // for float double and bool it depends on the base
    if (base == 10) {
      // decimals are well-defined for base 10
      return std::strtod(from.data(), &p_end);
    }
    // for all others decimal points will be ignored
    return std::strtol(from.data(), &p_end, base);
  }

#define AST_TEMPLATE(atype) template atype string2utype(const std::string& from, const uint base);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  string string2string(const string& from) {
    return from;
  }

  template <typename utype> vector<utype> vectorstring2vectorutype(const vector<string>& from) { // SD20220520
    vector<utype> vout;
    for (size_t i = 0; i < from.size(); i++) {
      vout.push_back(aurostd::string2utype<utype>(from[i]));
    }
    return vout;
  }
#define AST_TEMPLATE(utype) template vector<utype> vectorstring2vectorutype(const vector<string>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <typename utype> vector<utype> vectorstring2vectorutype(const deque<string>& from) { // SD20220520
    vector<utype> vout;
    for (size_t i = 0; i < from.size(); i++) {
      vout.push_back(aurostd::string2utype<utype>(from[i]));
    }
    return vout;
  }
#define AST_TEMPLATE(utype) template vector<utype> vectorstring2vectorutype(const deque<string>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  string vectorstring2string(const vector<string>& vstrings) {
    string out;
    for (size_t istr = 0; istr < vstrings.size(); istr++) {
      out += vstrings[istr];
    }
    return out;
  }
  string vectorstring2string(const deque<string>& vstrings) {
    string out;
    for (size_t istr = 0; istr < vstrings.size(); istr++) {
      out += vstrings[istr];
    }
    return out;
  }

  // ***************************************************************************
  // Function utype2string
  // ***************************************************************************

  //  template<typename string> string utype2string(const string& from) {
  //   return (string) stream2stream<string>(from);
  // }

  template <typename utype> string utype2string(const utype& from, int precision, char FORMAT) { // see DEFAULT_STREAM, FIXED_STREAM, SCIENTIFIC_STREAM
    return (string) stream2stream<string>(from, precision, FORMAT);
  }

#define AST_TEMPLATE(utype) template string utype2string(const utype& from, int precision, char FORMAT);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
  AST_GEN_1(AST_UTYPE_CHAR)
#undef AST_TEMPLATE

  //  string utype2string(const string& from) {
  //    return (string) from;
  //  }
  //  string utype2string(const std::basic_string<char, std::char_traits<char>, std::allocator<char> >& from) {    return (string) from;  }
  //  string utype2string(std::basic_string<char, std::char_traits<char>, std::allocator<char> > from) {    return (string) from;  }

  // cannot template this the same as others, char's don't make sense with roff
  string utype2string(double from, bool roff) {
    return utype2string(from, AUROSTD_DEFAULT_PRECISION, roff, DEFAULT_STREAM);
  }
  string utype2string(double from, int precision, bool roff) {
    return utype2string(from, precision, roff, AUROSTD_ROUNDOFF_TOL, DEFAULT_STREAM);
  }
  string utype2string(double from, bool roff, double tol) {
    return utype2string(from, AUROSTD_DEFAULT_PRECISION, roff, tol, DEFAULT_STREAM);
  }
  string utype2string(double from, int precision, bool roff, double tol) {
    return utype2string(from, precision, roff, tol, DEFAULT_STREAM);
  }
  string utype2string(double from, bool roff, char FORMAT) {
    return utype2string(from, AUROSTD_DEFAULT_PRECISION, roff, FORMAT);
  }
  string utype2string(double from, int precision, char FORMAT, bool roff) {
    return utype2string(from, precision, roff, FORMAT);
  } // CO20200624
  string utype2string(double from, int precision, bool roff, char FORMAT) {
    return utype2string(from, precision, roff, AUROSTD_ROUNDOFF_TOL, FORMAT);
  }
  string utype2string(double from, bool roff, double tol, char FORMAT) {
    return utype2string(from, AUROSTD_DEFAULT_PRECISION, roff, tol, FORMAT);
  }
  string utype2string(double from, int precision, bool roff, double tol, char FORMAT) {
    double tmp = from;
    if (roff) {
      tmp = roundoff(from, tol);
    }
    return (string) stream2stream<string>(tmp, precision, FORMAT);
  }
  string bool2string(bool from) {
    if (from) {
      return "true";
    }
    return "false";
  }

  // ***************************************************************************
  // Function utypes2deque
  // ***************************************************************************
  template <class utype> deque<utype> utypes2deque(utype u1) {
    deque<utype> out;
    out.push_back(u1);
    return out;
  }
  template <class utype> deque<utype> utypes2deque(utype u1, utype u2) {
    deque<utype> out;
    out.push_back(u1);
    out.push_back(u2);
    return out;
  }
  template <class utype> deque<utype> utypes2deque(utype u1, utype u2, utype u3) {
    deque<utype> out;
    out.push_back(u1);
    out.push_back(u2);
    out.push_back(u3);
    return out;
  }
  template <class utype> deque<utype> utypes2deque(utype u1, utype u2, utype u3, utype u4) {
    deque<utype> out;
    out.push_back(u1);
    out.push_back(u2);
    out.push_back(u3);
    out.push_back(u4);
    return out;
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorInt
  // ***************************************************************************
  void StringCommasColumsVectorInt(string vstring, vector<int>& vint) {
    vector<string> tokens_commas;
    vector<string> tokens_colums;
    vint.clear();
    string2tokens(vstring, tokens_commas, ",");
    for (size_t i = 0; i < tokens_commas.size(); i++) {
      tokens_colums.clear();
      if (aurostd::substring2bool(tokens_commas[i], ":")) {
        string2tokens(tokens_commas[i], tokens_colums, ":");
        for (int j = aurostd::string2utype<int>(tokens_colums[0]); j <= aurostd::string2utype<int>(tokens_colums.at(tokens_colums.size() - 1)); j++) {
          vint.push_back(j);
        }
      } else {
        vint.push_back(aurostd::string2utype<int>(tokens_commas[i]));
      }
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorUnsignedInt
  // ***************************************************************************
  void StringCommasColumsVectorUnsignedInt(string vstring, vector<uint>& vuint) {
    vector<string> tokens_commas;
    vector<string> tokens_colums;
    vuint.clear();
    string2tokens(vstring, tokens_commas, ",");
    for (size_t i = 0; i < tokens_commas.size(); i++) {
      tokens_colums.clear();
      if (aurostd::substring2bool(tokens_commas[i], ":")) {
        string2tokens(tokens_commas[i], tokens_colums, ":");
        for (uint j = aurostd::string2utype<uint>(tokens_colums[0]); j <= (uint) aurostd::string2utype<uint>(tokens_colums.at(tokens_colums.size() - 1)); j++) {
          vuint.push_back(j);
        }
      } else {
        vuint.push_back(aurostd::string2utype<uint>(tokens_commas[i]));
      }
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorFloat
  // ***************************************************************************
  void StringCommasColumsVectorFloat(string vstring, vector<float>& vfloat) {
    vector<string> tokens_commas;
    vector<string> tokens_colums;
    vfloat.clear();
    string2tokens(vstring, tokens_commas, ",");
    for (size_t i = 0; i < tokens_commas.size(); i++) {
      tokens_colums.clear();
      if (aurostd::substring2bool(tokens_commas[i], ":")) {
        string2tokens(tokens_commas[i], tokens_colums, ":");
        for (float j = aurostd::string2utype<float>(tokens_colums[0]); j <= aurostd::string2utype<float>(tokens_colums.at(tokens_colums.size() - 1)); j++) {
          vfloat.push_back(j);
        }
      } else {
        vfloat.push_back(aurostd::string2utype<float>(tokens_commas[i]));
      }
    }
  }

  // ***************************************************************************
  // Function StringCommasColumsVectorDouble
  // ***************************************************************************
  void StringCommasColumsVectorDouble(string vstring, vector<double>& vdouble) {
    vector<string> tokens_commas;
    vector<string> tokens_colums;
    vdouble.clear();
    string2tokens(vstring, tokens_commas, ",");
    for (size_t i = 0; i < tokens_commas.size(); i++) {
      tokens_colums.clear();
      if (aurostd::substring2bool(tokens_commas[i], ":")) {
        string2tokens(tokens_commas[i], tokens_colums, ":");
        for (double j = aurostd::string2utype<double>(tokens_colums[0]); j <= aurostd::string2utype<double>(tokens_colums.at(tokens_colums.size() - 1)); j++) {
          vdouble.push_back(j);
        }
      } else {
        vdouble.push_back(aurostd::string2utype<double>(tokens_commas[i]));
      }
    }
  }

  // ***************************************************************************
  // Function StringsAlphabetic
  // ***************************************************************************
  // says if two strings are in alphabetical order
  bool StringsAlphabetic(const string& A, const string& B, bool allow_identical) { // CO20181019
    if (A < B) {
      return true; // cerr << "A<B" << "  " << A << "<" << B << endl;
    }
    if (A > B) {
      return false; // cerr << "A>B" << "  " << A << ">" << B << endl;
    }
    if (A == B) {
      return allow_identical; // cerr << "A==B" << "  " << A << "==" << B << endl; //CO20181019
    }
    return true;
  }

  // ***************************************************************************
  // Function StringsAlphabetic
  // ***************************************************************************
  // says if two strings are in alphabetical order  //CO20180801
  bool StringsAlphabetic(const vector<string>& input, bool allow_identical) {
    for (size_t i = 1; i < input.size(); i++) { // CO20190218
      if (!StringsAlphabetic(input[i - 1], input[i], allow_identical)) {
        return false;
      }
    }
    return true;
  }
  bool StringsAlphabetic(const deque<string>& input, bool allow_identical) { // CO20190218
    for (size_t i = 1; i < input.size(); i++) {
      if (!StringsAlphabetic(input[i - 1], input[i], allow_identical)) {
        return false;
      }
    }
    return true;
  }

  /// @brief substitute parts of a strings in place
  /// @param work_string string to change in place
  /// @param old_string the string to search for
  /// @param new_string the string to replace the `old_string` with
  /// @authors
  /// @mod{SC,,created}
  /// @mod{HE,20240908,inPlace variant}
  void StringSubstInPlace(string& work_string, const string& old_string, const string& new_string) {
    string::size_type pos = 0;
    while ((pos = work_string.find(old_string, pos)) != std::string::npos) {
      work_string.replace(pos, old_string.length(), new_string);
      pos += new_string.length(); //  avoids replacing parts of `new_string` again
    }
  }

  /// @brief create a string with parts substituted
  /// @param work_string string to alter
  /// @param old_string the string to search for
  /// @param new_string the string to replace the `old_string` with
  /// @return altered string
  /// @authors
  /// @mod{HE,20220321,created}
  /// @mod{HE,20240908,updated naming}
  string StringSubst(const string& work_string, const string& old_string, const string& new_string) {
    std::string work_copy = work_string;
    StringSubstInPlace(work_copy, old_string, new_string);
    return work_copy;
  }

  /// @brief change a string in place by replacing a single character by a new one
  /// @param work_string string to change in place
  /// @param old_char the single char to search for
  /// @param new_char the char to replace `find_char` with
  /// @authors
  /// @mod{SC,,created}
  /// @mod{HE,20240908,simplefied logic}
  void StringSubstInPlace(string& work_string, const char old_char, const char new_char) {
    std::replace(work_string.begin(), work_string.end(), old_char, new_char);
  }

  /// @brief create a string with a single character replaced by a new one
  /// @param work_string string to change in place
  /// @param old_char the single char to search for
  /// @param new_char the char to replace `find_char` with
  /// @authors
  /// @mod{HE,20220321,created}
  /// @mod{HE,20240908,updated naming}
  string StringSubst(const string& work_string, const char old_char, const char new_char) {
    std::string work_copy = work_string;
    StringSubstInPlace(work_copy, old_char, new_char);
    return work_copy;
  }

  void StringStreamSubst(stringstream& strstringstream, const string& strfind, const string& strreplace) {
    string strstring = strstringstream.str();
    StringSubstInPlace(strstring, strfind, strreplace);
    aurostd::StringstreamClean(strstringstream);
    strstringstream << strstring;
  }

  // ***************************************************************************
  // Function SubStrings
  // ***************************************************************************

  /// @brief count number of newlines in a string
  /// @param str to count
  /// @authors
  /// @mod{SC,20190409,created}
  /// @mod{SD,20250331,updated using standard function}
  size_t GetNLinesString(const string& str) {
    const std::string::difference_type n = std::count(str.begin(), str.end(), '\n');
    return static_cast<size_t>(n);
  }
  size_t GetNLinesString(const stringstream& strstream) {
    return aurostd::GetNLinesString(strstream.str());
  }

  size_t GetNLinesFile(const string& file_name) {
    string stringFILE;
    aurostd::file2string(file_name, stringFILE);
    return GetNLinesString(stringFILE);
  }

  string GetLineString(const string& strstream, int line) {
    string _strstream(strstream);
    string _strline;
    //  if(line>aurostd::GetNLinesString(_strstream)) return (string) "";   // TOO SLOW IF THE STRING IS LONG !
    for (int i = 0; i < line; i++) {
      _strline = _strstream.substr(0, _strstream.find("\n"));
      _strstream = _strstream.substr(_strstream.find("\n") + 1);
    }
    return _strline;
  }
  string GetLineString(const stringstream& strstream, int line) {
    return aurostd::GetLineString(strstream.str(), line);
  }

  // ***************************************************************************
  // Function SubStringsPresent
  // ***************************************************************************
  bool substring2bool(const string& strstream, const string& strsub1, bool RemoveWS, bool RemoveComments) { // CO20210315 - cleaned up
    // substring2bool and kvpair2bool are similar but distinct
    // substring2bool will match any strsub1 and return true
    // kvpair2bool will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return true
    // matching KEY exactly is useful, e.g.:
    //_FILE_START_
    // IALGO=48
    //_FILE_END_
    // strsub1="ALGO": substring2bool will return true
    // keyword="ALGO,delim="=": kvpair2bool will return false
    // kvpair2bool must match KEY exactly! skips the rest
    // substring2bool is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    const bool LDEBUG = false; // true;
    if (LDEBUG) {
      cerr << XPID << "aurostd::substring2bool(): BEGIN [substring=\"" << strsub1 << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    string _strstream(strstream);
    if (RemoveWS == true) {
      _strstream = aurostd::RemoveWhiteSpaces(_strstream, '"');
    }
    if (LDEBUG) {
      cerr << XPID << "aurostd::substring2bool(): [input=\"" << strstream << "\"], [substring=\"" << strsub1 << "\"]" << endl;
    }
    if (_strstream.find(strsub1) == string::npos) {
      return false;
    }
    if (RemoveComments) { // SD20220403 - substring exists, but now check if it exists outside of comments
      vector<string> tokens;
      aurostd::string2tokens(_strstream, tokens, "\n");
      string strline;
      for (size_t i = 0; i < tokens.size(); i++) {
        strline = aurostd::RemoveComments(tokens[i]); // CO20210315
        if (strline.find(strsub1) != string::npos) {
          if (LDEBUG) {
            cerr << XPID << "aurostd::substring2bool(): END [substring=\"" << strsub1 << "\" found] [RemoveWS=" << RemoveWS << "]" << endl;
          }
          return true;
        }
      }
      if (LDEBUG) {
        cerr << XPID << "aurostd::substring2bool(): END [substring=" << strsub1 << " NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return false;
    }
    return true; // SD20220403 - since substring exists, return true
  }

  bool substring2bool(const vector<string>& vstrstream, const string& strsub1, bool RemoveWS, bool RemoveComments) {
    for (size_t i = 0; i < vstrstream.size(); i++) {
      if (aurostd::substring2bool(vstrstream[i], strsub1, RemoveWS, RemoveComments)) {
        return true;
      }
    }
    return false;
  }
  bool substring2bool(const deque<string>& vstrstream, const string& strsub1, bool RemoveWS, bool RemoveComments) {
    for (size_t i = 0; i < vstrstream.size(); i++) {
      if (aurostd::substring2bool(vstrstream[i], strsub1, RemoveWS, RemoveComments)) {
        return true;
      }
    }
    return false;
  }

  // ME20220505
  // Matches a list of substrings to a string
  // match_all: only substring must be inside the string
  bool substringlist2bool(const string& strin, const vector<string>& substrings, bool match_all) {
    for (size_t i = 0; i < substrings.size(); i++) {
      if (strin.find(substrings[i]) == string::npos) {
        // Didn't find substring, but need all
        if (match_all) {
          return false;
        }
      } else if (!match_all) {
        // Found something and only need one
        return true;
      }
    }
    // Code only gets here when all substrings are
    // in the string (match_all) or none have (!match_all)
    return match_all;
  }

  bool substring2bool(const stringstream& strstream, const string& strsub1, bool RemoveWS, bool RemoveComments) {
    return aurostd::substring2bool(strstream.str(), strsub1, RemoveWS, RemoveComments);
  }

  bool substringlist2bool(const string& strin, const deque<string>& substrings, bool match_all) {
    for (size_t i = 0; i < substrings.size(); i++) {
      if (strin.find(substrings[i]) == string::npos) {
        // Didn't find substring, but need all
        if (match_all) {
          return false;
        }
      } else if (!match_all) {
        // Found something and only need one
        return true;
      }
    }
    // Code only gets here when all substrings are
    // in the string (match_all) or none have (!match_all)
    return match_all;
  }

  /// @brief Determines whether a specific value exists inside a list
  ///
  /// @param list List of values
  /// @param input Specific value to match the values within the list
  /// @param index Index of the element that matched the input
  /// @param sorted Whether the list is sorted
  /// @return Boolean if the input was found in the list or not
  ///
  /// @note index is assigned the value of 0 if the input is not found in the list
  template <class utype> bool WithinList(const vector<utype>& list, const utype& input, size_t& index, bool sorted) { // SD20220705
    for (size_t i = 0; i < list.size(); i++) {
      if (sorted && list[i] > input) {
        break;
      }
      if (aurostd::isequal(list[i], input)) {
        index = i;
        return true;
      }
    }
    index = 0;
    return false;
  }
  bool WithinList(const vector<string>& list, const string& input, size_t& index, bool sorted) { // SD20220705
    for (size_t i = 0; i < list.size(); i++) {
      if (sorted && list[i] > input) {
        break;
      }
      if (aurostd::isequal(list[i], input)) {
        index = i;
        return true;
      }
    }
    index = 0;
    return false;
  }
  template <class utype> bool WithinList(const deque<utype>& list, const utype& input, size_t& index, bool sorted) { // SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
  bool WithinList(const deque<string>& list, const string& input, size_t& index, bool sorted) { // SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
  template <class utype> bool WithinList(const vector<utype>& list, const utype& input, bool sorted) { // SD20220705
    size_t index = 0;
    return WithinList(list, input, index, sorted);
  }
  bool WithinList(const vector<string>& list, const string& input, bool sorted) { // SD20220705
    size_t index = 0;
    return WithinList(list, input, index, sorted);
  }
  template <class utype> bool WithinList(const deque<utype>& list, const utype& input, bool sorted) { // SD20220705
    return WithinList(aurostd::deque2vector(list), input, sorted);
  }
  bool WithinList(const deque<string>& list, const string& input, bool sorted) { // SD20220705
    return WithinList(aurostd::deque2vector(list), input, sorted);
  }
  template <class utype> bool WithinList(const vector<utype>& list, const utype& input, vector<size_t>& index, bool sorted) { // SD20220705
    index.clear();
    for (size_t i = 0; i < list.size(); i++) {
      if (sorted && list[i] > input) {
        break;
      }
      if (aurostd::isequal(list[i], input)) {
        index.push_back(i);
      }
    }
    return !index.empty();
  }
  bool WithinList(const vector<string>& list, const string& input, vector<size_t>& index, bool sorted) { // SD20220705
    index.clear();
    for (size_t i = 0; i < list.size(); i++) {
      if (sorted && list[i] > input) {
        break;
      }
      if (aurostd::isequal(list[i], input)) {
        index.push_back(i);
      }
    }
    return !index.empty();
  }
  template <class utype> bool WithinList(const deque<utype>& list, const utype& input, vector<size_t>& index, bool sorted) { // SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
  bool WithinList(const deque<string>& list, const string& input, vector<size_t>& index, bool sorted) { // SD20220705
    return WithinList(aurostd::deque2vector(list), input, index, sorted);
  }
#define AST_TEMPLATE(utype) template bool WithinList(const deque<utype>&, const utype&, size_t&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template bool WithinList(const vector<utype>&, const utype&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template bool WithinList(const deque<utype>&, const utype&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template bool WithinList(const vector<utype>&, const utype&, vector<size_t>&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template bool WithinList(const deque<utype>&, const utype&, vector<size_t>&, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // ME20220503
  bool SubstringWithinList(const deque<string>& list, const string& input) {
    int index = -1;
    return SubstringWithinList(list, input, index);
  }

  bool SubstringWithinList(const deque<string>& list, const string& input, int& index) {
    for (deque<string>::const_iterator it = list.begin(); it != list.end(); ++it) {
      if ((*it).find(input) != string::npos) {
        index = std::distance(list.begin(), it);
        return true;
      }
    }
    index = -1;
    return false;
  }

  bool SubstringWithinList(const vector<string>& list, const string& input) {
    int index = -1;
    return SubstringWithinList(list, input, index);
  }

  bool SubstringWithinList(const vector<string>& list, const string& input, int& index) {
    for (vector<string>::const_iterator it = list.begin(); it != list.end(); ++it) {
      if ((*it).find(input) != string::npos) {
        index = std::distance(list.begin(), it);
        return true;
      }
    }
    index = -1;
    return false;
  }

  // ***************************************************************************
  // Function substring_present_file
  // ***************************************************************************
  /// @brief determine if a substring is present in the file
  /// @param FileNameRaw file to inspect
  /// @param boyer_moore boyer_moore_searcher of the substring
  /// @param strsub_size substring size
  /// @param offset where to start reading the file
  /// @return true if the substring was found
  /// @authors
  /// @mod{CO,20210315,created function}
  /// @mod{SD,20240514,rewritten using ifstream and boyer_moore_searcher}
  bool substring_present_file(const string& FileNameRaw, const std::boyer_moore_searcher<std::string::const_iterator>& boyer_moore, const size_t strsub_size, size_t& offset) {
    const string FileName(aurostd::CleanFileName(FileNameRaw));
    static char buff[10485760]; // 1MB
    const size_t count = sizeof(buff);
    if (!aurostd::FileExist(FileName)) {
      const string message = "file input not found =" + FileName;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    ifstream is(FileName, std::ios::binary);
    bool found = false;
    while (!found && is) {
      is.seekg(offset);
      is.read(buff, sizeof(buff));
      if (const auto it = std::search(std::begin(buff), std::begin(buff) + is.gcount(), boyer_moore); it != std::begin(buff) + is.gcount()) {
        found = true;
      }
      offset = static_cast<std::size_t>(is.tellg()) - strsub_size;
    }
    is.close();
    return found;
  }

  /// @brief determine if a substring is present in the file
  /// @param FileName file to inspect
  /// @param strsub substring to find
  /// @param offset where to start reading the file
  /// @return true if the substring was found
  /// @authors
  /// @mod{SD,20240514,created function}
  bool substring_present_file(const string& FileName, const string& strsub, size_t& offset) {
    const std::boyer_moore_searcher boyer_moore(strsub.begin(), strsub.end());
    return substring_present_file(FileName, boyer_moore, strsub.size(), offset);
  }

  /// @brief determine if a substring is present in the file
  /// @param FileName file to inspect
  /// @param strsub substring to find
  /// @return true if the substring was found
  /// @authors
  /// @mod{SD,20240514,created function}
  bool substring_present_file(const string& FileName, const string& strsub) {
    size_t offset = 0;
    return substring_present_file(FileName, strsub, offset);
  }

  // ***************************************************************************
  // Function substrings_present_file
  // ***************************************************************************
  /// @brief determine if substrings are present in the file
  /// @param FileNameRaw file to inspect
  /// @param vboyer_moore vector of boyer_moore_searchers of the substrings
  /// @param strsub_size maximum substring size
  /// @param offset where to start reading the file
  /// @return vector of booleans, true if a substring was found
  /// @authors
  /// @mod{SD,20240514,created function}
  vector<bool> substrings_present_file(const string& FileNameRaw, const vector<std::boyer_moore_searcher<std::string::const_iterator>>& vboyer_moore, const size_t strsub_size, size_t& offset) {
    const string FileName(aurostd::CleanFileName(FileNameRaw));
    static char buff[10485760]; // 1MB
    if (!aurostd::FileExist(FileName)) {
      const string message = "file input not found =" + FileName;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    ifstream is(FileName, std::ios::binary);
    vector<bool> found(vboyer_moore.size(), false);
    while (std::all_of(found.begin(), found.end(), [](bool x) { return !x; }) && is) {
      is.seekg(offset);
      is.read(buff, sizeof(buff));
      for (size_t i = 0; i < vboyer_moore.size(); i++) {
        if (const auto it = std::search(std::begin(buff), std::begin(buff) + is.gcount(), vboyer_moore[i]); it != std::begin(buff) + is.gcount()) {
          found[i] = true;
        }
      }
      offset = static_cast<std::size_t>(is.tellg()) - strsub_size;
    }
    is.close();
    return found;
  }

  /// @brief determine if substrings are present in the file
  /// @param FileName file to inspect
  /// @param vstrsub substrings to find
  /// @param offset where to start reading the file
  /// @return vector of booleans, true if a substring was found
  /// @authors
  /// @mod{SD,20240514,created function}
  vector<bool> substrings_present_file(const string& FileName, const vector<string>& vstrsub, size_t& offset) {
    vector<std::boyer_moore_searcher<std::string::const_iterator>> vboyer_moore;
    for (const string& strsub : vstrsub) {
      vboyer_moore.emplace_back(strsub.begin(), strsub.end());
    }
    const size_t strsub_size = (*std::max_element(vstrsub.begin(), vstrsub.end(), [](const string& x1, const string& x2) { return x1.size() < x2.size(); })).size();
    return substrings_present_file(FileName, vboyer_moore, strsub_size, offset);
  }

  /// @brief determine if substrings are present in the file
  /// @param FileName file to inspect
  /// @param vstrsub substrings to find
  /// @return vector of booleans, true if a substring was found
  /// @authors
  /// @mod{SD,20240514,created function}
  vector<bool> substrings_present_file(const string& FileName, const vector<string>& vstrsub) {
    size_t offset = 0;
    return substrings_present_file(FileName, vstrsub, offset);
  }

  // ***************************************************************************
  // Function substrings_map_present_file
  // ***************************************************************************
  /// @brief determine if substrings are present in the file
  /// @param FileNameRaw file to inspect
  /// @param map_boyer_moore map of <substring, pair of <boyer_moore_searcher of substring, substring present>>
  /// @param strsub_size maximum substring size
  /// @param offset where to start reading the file
  /// @authors
  /// @mod{SD,20240514,created function}
  void substrings_map_present_file(const string& FileNameRaw, std::map<string, std::pair<std::boyer_moore_searcher<std::string::const_iterator>, bool>>& map_boyer_moore, const size_t strsub_size, size_t& offset) {
    const string FileName(aurostd::CleanFileName(FileNameRaw));
    static char buff[10485760]; // 1MB
    if (!aurostd::FileExist(FileName)) {
      const string message = "file input not found =" + FileName;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    ifstream is(FileName, std::ios::binary);
    while (std::all_of(map_boyer_moore.begin(), map_boyer_moore.end(), [](const auto& x) { return !x.second.second; }) && is) {
      is.seekg(offset);
      is.read(buff, sizeof(buff));
      for (auto& entry : map_boyer_moore) {
        if (entry.second.second) {
          continue;
        }
        if (const auto it = std::search(std::begin(buff), std::begin(buff) + is.gcount(), entry.second.first); it != std::begin(buff) + is.gcount()) {
          entry.second.second = true;
        }
      }
      offset = static_cast<std::size_t>(is.tellg()) - strsub_size;
    }
    is.close();
  }

  /// @brief determine if substrings are present in the file
  /// @param FileName file to inspect
  /// @param map_strsub unordered map of <substring, substring present>
  /// @param offset where to start reading the file
  /// @authors
  /// @mod{SD,20240514,created function}
  void substrings_map_present_file(const string& FileName, std::unordered_map<string, bool>& map_strsub, size_t& offset) {
    std::map<string, std::pair<std::boyer_moore_searcher<std::string::const_iterator>, bool>> map_boyer_moore;
    for (const auto& entry : map_strsub) {
      map_boyer_moore.emplace(entry.first, std::pair<std::boyer_moore_searcher<std::string::const_iterator>, bool>(std::boyer_moore_searcher(entry.first.begin(), entry.first.end()), entry.second));
    }
    const size_t strsub_size = (*std::max_element(map_strsub.begin(), map_strsub.end(), [](const auto& x1, const auto& x2) { return x1.first.size() < x2.first.size(); })).first.size();
    substrings_map_present_file(FileName, map_boyer_moore, strsub_size, offset);
    for (const auto& entry : map_boyer_moore) {
      map_strsub[entry.first] = entry.second.second;
    }
  }

  /// @brief determine if substrings are present in the file
  /// @param FileName file to inspect
  /// @param map_strsub unordered map of <substring, substring present>
  /// @authors
  /// @mod{SD,20240514,created function}
  void substrings_map_present_file(const string& FileName, std::unordered_map<string, bool>& map_strsub) {
    size_t offset = 0;
    substrings_map_present_file(FileName, map_strsub, offset);
  }

  // ***************************************************************************
  // Function SubStringsPresent and EXTRACT
  // ***************************************************************************
  // SD20220520
  // Rewritten substring2string to be more understandable, accept more input types, and incorporate extracting Nth entries
  // n==0 returns all entries, starts counting from 1, negative numbers go backwards
  // Original substring2string always returned just the first entry
  // When two substrings are present strsub1 is the start keyword and strsub2 is the stop keyword
  // CO20230502 - another rewrite for strsub1 and strsub2
  // substring2string() with strsub1 grabs everything on the line
  // substring2string() with strsub1 and strsub2 grabs everything in between
  // substring2strings() with strsub1 and strsub2 now handles newlines within entries of the return vector, as well as start/stop tags within the same line
  string substring2string(ifstream& input, const string& strsub1, const int instance, bool RemoveWS, bool RemoveComments) {
    // substring2string and kvpair2string are similar but distinct
    // substring2string will match any strsub1 and return everything AFTER strsub1
    // kvpair2string will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return VALUE
    // matching KEY exactly is useful, e.g.:
    //_FILE_START_
    // IALGO==48
    // ALGO==FAST
    //_FILE_END_
    // strsub1="ALGO",instance=1: substring2string will return "==48"
    // keyword="ALGO,delim="==": kvpair2string will return "FAST"
    // kvpair2string must match KEY exactly! skips the rest
    // substring2string is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [substring=\"" << strsub1 << "\"]" << endl;
    }
    string strline;
    input.clear();
    input.seekg(0);
    vector<string> tokens;
    int iter = 0;
    while (getline(input, strline) && (instance == 0 || iter != instance)) {
      if (RemoveWS) {
        strline = aurostd::RemoveWhiteSpaces(strline, '"');
      }
      if (RemoveComments) {
        strline = aurostd::RemoveComments(strline);
      }
      if (strline.find(strsub1) != string::npos) {
        tokens.push_back(strline.substr(strline.find(strsub1) + strsub1.length()));
        iter++;
      }
    }
    input.clear();
    input.seekg(0);
    if (tokens.empty() || (uint) aurostd::abs(instance) > tokens.size()) {
      return "";
    }
    stringstream strstream;
    if (instance == 0) {
      for (size_t i = 0; i < tokens.size(); i++) {
        strstream << tokens[i] << endl;
      }
    } else if (instance > 0) {
      strstream << tokens[tokens.size() - 1];
    } else if (instance < 0) {
      const uint i = (uint) aurostd::boundary_conditions_periodic(0, tokens.size() - 1, instance);
      strstream << tokens[i];
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    return strstream.str();
  }

  string substring2string(const string& _input, const string& strsub1, const int instance, bool RemoveWS, bool RemoveComments) {
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    string input = _input;
    if (RemoveWS) {
      input = aurostd::RemoveWhiteSpaces(_input, '"');
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [input=\"" << input << "\"] [substring=\"" << strsub1 << "\"]" << endl;
    }
    if (input.find(strsub1) == string::npos) {
      return "";
    }
    stringstream strstream;
    vector<string> tokens;
    aurostd::string2vectorstring(input, tokens);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [tokens.size()=" << tokens.size() << "]" << endl;
    }
    int iter = 0;
    if (instance > 0) {
      for (size_t i = 0; i < tokens.size(); i++) {
        if (RemoveComments) {
          tokens[i] = aurostd::RemoveComments(tokens[i]);
        }
        if (tokens[i].find(strsub1) != string::npos) {
          iter++;
          if (instance == iter) {
            strstream << tokens[i].substr(tokens[i].find(strsub1) + strsub1.length());
            break;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return strstream.str();
    } else if (instance < 0) {
      for (int i = tokens.size() - 1; i >= 0; i--) {
        if (RemoveComments) {
          tokens[i] = aurostd::RemoveComments(tokens[i]);
        }
        if (tokens[i].find(strsub1) != string::npos) {
          iter--;
          if (instance == iter) {
            strstream << tokens[i].substr(tokens[i].find(strsub1) + strsub1.length());
            break;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return strstream.str();
    } else { // instance==0
      for (size_t i = 0; i < tokens.size(); i++) {
        if (RemoveComments) {
          tokens[i] = aurostd::RemoveComments(tokens[i]);
        }
        if (tokens[i].find(strsub1) != string::npos) {
          strstream << tokens[i].substr(tokens[i].find(strsub1) + strsub1.length()) << endl;
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " END [substring=\"" << strsub1 << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return strstream.str();
    }
    return "";
  }

  string substring2string(const stringstream& input, const string& strsub1, const int instance, bool RemoveWS, bool RemoveComments) {
    return substring2string(input.str(), strsub1, instance, RemoveWS, RemoveComments);
  }

  /// @brief find string between two substrings of a text
  /// @param input input file
  /// @param strsub1 first substring
  /// @param strsub2 second substring
  /// @param instance which instance (nonzero; if negative, search from the back)
  /// @param RemoveWS whether to remove white spaces
  /// @param RemoveComments whether to remove comments
  /// @return string that was found (empty if substrings not found)
  ///
  /// @authors
  /// @mod{SG,20240517,created}
  string substring2string(ifstream& input, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;
    }
    if (instance == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "instance must be nonzero integer", _INPUT_ILLEGAL_);
    }
    if (aurostd::substring2bool(strsub1, "\n") || aurostd::substring2bool(strsub2, "\n")) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "substring2string with ifstream input does not support newline characters inside of substrings", _INPUT_ILLEGAL_);
    }
    size_t start_pos;
    size_t pos1 = 0;
    size_t pos2 = 0;
    string substr;
    string strline;
    bool one_or_two = true; // if true, start searching the line for strsub1; if false, strsub2
    // bool read = false;  // not sure if it's important whether the ifstream lines have been "read"
    input.clear();
    int current_instance = 0;
    if (instance > 0) {
      input.seekg(0);
      while (getline(input, strline)) { // loop over lines in file
        start_pos = 0;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        if (RemoveComments) {
          strline = aurostd::RemoveComments(strline);
        }
        if (RemoveWS) {
          strline = aurostd::RemoveWhiteSpaces(strline, '"');
        }
        while (true) { // loop over instances of strsub1 and strsub2 in line. If break is called, then we move on to the next line (unless we're done).
          if (one_or_two) {
            // start searching the line for strsub1
            pos1 = strline.find(strsub1, start_pos);
            if (pos1 == string::npos) {
              // if strsub1 not found, move on to the next line and search for strsub1 again
              break;
            }
            pos2 = strline.find(strsub2, pos1 + strsub1.size());
            if (pos2 == string::npos) {
              // if strsub1 found but not strsub2, move on to the next line and search for strsub2
              one_or_two = false;
              substr = strline.substr(pos1 + strsub1.size(), string::npos);
              break;
            }
            // if both strsub1 and strsub2 found in line, increase the instance and if we're not done, continue searching the same line
            substr = strline.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
            current_instance++;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos2 + strsub2.size();
          } else {
            // start searching the line for strsub2
            pos1 = strline.find(strsub2, 0);
            if (pos1 == string::npos) {
              // if strsub2 not found in line, move on to the next line and search for strsub2 again
              substr += strline;
              break;
            }
            // if strsub2 found, increase the instance and if we're not done, search the same line for strsub1
            substr += strline.substr(0, pos1);
            current_instance++;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos1 + strsub2.size();
            one_or_two = true;
          }
          if (current_instance == instance) {
            break;
          }
        }
        if (current_instance == instance) {
          break;
        }
        if (!substr.empty()) {
          substr += "\n";
        } // preserve the newlines, because getline(ifstream) removes them
      }
      if (current_instance < instance) {
        substr = "";
      } // this line probably not needed
    } else if (instance < 0) {
      // we have to read the file backwards!
      string thischar;
      input.seekg(-1, std::ios::end);
      const int numchars = input.tellg();
      for (int i = 0; i < numchars; i++) { // loop over instances of strsub1 and strsub2 in line. If break is called, then we move on to the previous line (unless we're done).
        thischar = input.get(); // note that causes current line (std::ios::cur) to move forward by 1
        if (thischar != "\n") {
          // we have to read char by char, and assemble the line string until a newline is encountered
          strline = thischar + strline;
          input.seekg(-2, std::ios::cur);
          continue;
        }
        start_pos = strline.size();
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        if (RemoveComments) {
          strline = aurostd::RemoveComments(strline);
        }
        if (RemoveWS) {
          strline = aurostd::RemoveWhiteSpaces(strline, '"');
        }
        while (true) { // loop over instances of strsub1 and strsub2 in line. If break is called, then we move on to the previous line (unless we're done).
          if (one_or_two) {
            // start searching the line (backwards) for strsub2
            pos2 = strline.rfind(strsub2, start_pos);
            if (pos2 == string::npos) {
              // if strsub2 not found, move on to the previous line and search for strsub2 again
              break;
            }
            pos1 = strline.rfind(strsub1, pos2);
            if (pos1 == string::npos) {
              // if strsub2 found but not strsub1, move on to the previous line and search for strsub1
              one_or_two = false;
              substr = strline.substr(0, pos2);
              break;
            }
            // if both strsub2 and strsub1 found in line, decrease the instance and if we're not done, continue searching the same line
            substr = strline.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
            current_instance--;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos1;
          } else {
            // start searching the line (backwards) for strsub1
            pos1 = strline.rfind(strsub1, start_pos);
            if (pos1 == string::npos) {
              // if strsub1 not found in line, move on to the previous line and search for strsub1 again
              substr = strline + substr;
              break;
            }
            // if strsub1 found, decrease the instance and if we're not done, search the same line for strsub2
            substr = strline.substr(0, strline.size() - pos1) + substr;
            current_instance--;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos1;
            one_or_two = true;
          }
          if (current_instance == instance) {
            break;
          }
        }
        if (current_instance == instance) {
          break;
        }
        if (!substr.empty()) {
          substr = "\n" + substr;
        } // preserve the newlines
        input.seekg(-2, std::ios::cur); // move backwards 2 lines
        strline = "";
      }
      if (current_instance > instance) {
        substr = "";
      } // this line probably not needed
      if (!substr.empty()) {
        // In the reverse case, the substr will include strsub1 and strsub2, so we need to remove them
        if (substr.substr(0, strsub1.size()) == strsub1) {
          substr.erase(0, strsub1.size());
        }
        if (substr.substr(substr.size() - strsub2.size(), strsub2.size()) == strsub2) {
          substr.erase(substr.size() - strsub2.size(), strsub2.size());
        }
      }
    }
    substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
    substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [output=\"" << substr << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    return substr;
  }

  /// @brief find string between two substrings of a text
  /// @param input input text
  /// @param strsub1 first substring
  /// @param strsub2 second substring
  /// @param instance which instance (nonzero; if negative, search from the back)
  /// @param RemoveWS whether to remove white spaces
  /// @param RemoveComments whether to remove comments
  /// @return string that was found (empty if substrings not found)
  ///
  /// @authors
  /// @mod{SG,20240517,created}
  string substring2string(const std::string& input, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;
    }
    if (instance == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "instance must be nonzero integer", _INPUT_ILLEGAL_);
    }
    std::string input_mod = input;
    if (RemoveComments) {
      input_mod = aurostd::RemoveComments(input_mod);
    }
    if (RemoveWS) {
      input_mod = aurostd::RemoveWhiteSpaces(input_mod, '"');
    }
    size_t start_pos;
    size_t pos1 = 0;
    size_t pos2 = 0;
    std::string substr;
    if (instance > 0) {
      start_pos = 0;
      for (size_t i = 1; i <= instance; i++) {
        pos1 = input_mod.find(strsub1, start_pos);
        pos2 = input_mod.find(strsub2, pos1 + strsub1.size());
        if (pos1 == string::npos or pos2 == string::npos) {
          substr = "";
          break;
        }
        substr = input_mod.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
        start_pos = pos2 + strsub2.size();
      }
    } else if (instance < 0) {
      start_pos = input_mod.size();
      for (size_t i = 1; i <= (-1) * instance; i++) {
        pos2 = input_mod.rfind(strsub2, start_pos);
        pos1 = input_mod.rfind(strsub1, pos2);
        if (pos1 == string::npos or pos2 == string::npos) {
          substr = "";
          break;
        }
        substr = input_mod.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
        start_pos = pos1;
      }
    }
    substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
    substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [output=\"" << substr << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    return substr;
  }

  /// @brief find string between two substrings of a text
  /// @param input input text
  /// @param strsub1 first substring
  /// @param strsub2 second substring
  /// @param instance which instance (nonzero; if negative, search from the back)
  /// @param RemoveWS whether to remove white spaces
  /// @param RemoveComments whether to remove comments
  /// @return string that was found (empty if substrings not found)
  ///
  /// @authors
  /// @mod{SD,20220530,created}
  /// @mod{SG,20240517,added documentation and removed trim_edges}
  string substring2string(const stringstream& input, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    return substring2string(input.str(), strsub1, strsub2, instance, RemoveWS, RemoveComments);
  }

  template <typename utype> utype substring2utype(ifstream& input, const string& strsub1, const int instance, bool RemoveWS, bool RemoveComments) {
    return string2utype<utype>(substring2string(input, strsub1, instance, RemoveWS, RemoveComments));
  }
  template <typename utype> utype substring2utype(const string& input, const string& strsub1, const int instance, bool RemoveWS, bool RemoveComments) {
    return string2utype<utype>(substring2string(input, strsub1, instance, RemoveWS, RemoveComments));
  }
  template <typename utype> utype substring2utype(const stringstream& input, const string& strsub1, const int instance, bool RemoveWS, bool RemoveComments) {
    return substring2utype<utype>(input.str(), strsub1, instance, RemoveWS, RemoveComments);
  }
  template <typename utype> utype substring2utype(ifstream& input, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    return string2utype<utype>(substring2string(input, strsub1, strsub2, instance, RemoveWS, RemoveComments));
  }
  template <typename utype> utype substring2utype(const string& input, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    return string2utype<utype>(substring2string(input, strsub1, strsub2, instance, RemoveWS, RemoveComments));
  }
  template <typename utype> utype substring2utype(const stringstream& input, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    return substring2utype<utype>(input.str(), strsub1, strsub2, instance, RemoveWS, RemoveComments);
  }

#define AST_TEMPLATE(atype) template atype substring2utype(ifstream&, const string&, const int, bool, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype substring2utype(const string&, const string&, const int, bool, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype substring2utype(const stringstream&, const string&, const int, bool, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype substring2utype(ifstream&, const string&, const string&, const int, bool, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype substring2utype(const string&, const string&, const string&, const int, bool, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype substring2utype(const stringstream&, const string&, const string&, const int, bool, bool);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  bool kvpair2bool(ifstream& input, const string& keyword, const string& delim, bool RemoveWS, bool RemoveComments) { // SD20220520
    // substring2bool and kvpair2bool are similar but distinct
    // substring2bool will match any strsub1 and return true
    // kvpair2bool will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return true
    // matching KEY exactly is useful, e.g.:
    //_FILE_START_
    // IALGO==48
    //_FILE_END_
    // strsub1="ALGO": substring2bool will return true
    // keyword="ALGO,delim="==": kvpair2bool will return false
    // kvpair2bool must match KEY exactly! skips the rest
    // substring2bool is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;
    }
    string strline;
    string _keyword;
    string::size_type idx = 0;
    while (getline(input, strline)) {
      if (RemoveWS) {
        strline = aurostd::RemoveWhiteSpaces(strline, '"');
      }
      if (RemoveComments) {
        strline = aurostd::RemoveComments(strline);
      }
      idx = strline.find(delim);
      if (idx != string::npos) {
        _keyword = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0, idx));
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << "_keyword=\"" << _keyword << "\"" << endl;
        }
        if (_keyword == keyword) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" found] [RemoveWS=" << RemoveWS << "]" << endl;
          }
          return true;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    return false;
  }

  bool kvpair2bool(const string& input, const string& keyword, const string& delim, bool RemoveWS, bool RemoveComments) { // CO20210315
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    string _input(input);
    if (RemoveWS == true) {
      _input = aurostd::RemoveWhiteSpaces(_input, '"');
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"], [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;
    }
    if (_input.find(keyword) == string::npos) {
      return false;
    }

    if (RemoveComments) { // SD20220403 - substring exists, but now check if it exists outside of comments
      vector<string> tokens;
      aurostd::string2tokens(_input, tokens, "\n");
      string strline;
      string _keyword;
      string::size_type idx = 0;
      for (size_t i = 0; i < tokens.size(); i++) {
        strline = aurostd::RemoveComments(tokens[i]);
        idx = strline.find(delim);
        if (idx != string::npos) {
          _keyword = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0, idx));
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << "_keyword=\"" << _keyword << "\"" << endl;
          }
          if (_keyword == keyword) {
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" found] [RemoveWS=" << RemoveWS << "]" << endl;
            }
            return true;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << "END [keyword=\"" << keyword << "\" NOT found] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return false;
    }
    return true; // SD20220403 - since substring exists, return true
  }

  bool kvpair2bool(const stringstream& input, const string& keyword, const string& delim, bool RemoveWS, bool RemoveComments) { // CO20210315 - cleaned up
    return kvpair2bool(input.str(), keyword, delim, RemoveWS, RemoveComments);
  }

  // SD20220520
  // Rewritten kvpair2string to be more understandable, accept more input types, allow for multi-char delimiters and incorporate extracting Nth entries
  // n==0 returns all entries, starts counting from 1, negative numbers go backwards
  // Original kvpair2string always returned just the first entry
  string kvpair2string(ifstream& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments) {
    // substring2string and kvpair2string are similar but distinct
    // substring2string will match any strsub1 and return everything AFTER strsub1
    // kvpair2string will assume the line is written KEY+DELIMITER+VALUE, if it matches KEY and DELIMITER exactly, it will return VALUE
    // matching KEY exactly is useful, e.g.:
    //_FILE_START_
    // IALGO==48
    // ALGO==FAST
    //_FILE_END_
    // strsub1="ALGO",instance=1: substring2string will return "==48"
    // keyword="ALGO,delim="==": kvpair2string will return "FAST"
    // kvpair2string must match KEY exactly! skips the rest
    // substring2string is good for aflow.in's which has no set delimiter style: [AFLOW_BIN_XZ] vs. [AFLOW_BIN=XZ] vs. [AFLOW_BIN]XZ vs. [AFLOW]BIN=XZ
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;
    }
    string strline;
    string _keyword;
    input.clear();
    input.seekg(0);
    vector<string> tokens;
    string::size_type idx = 0;
    int iter = 0;
    while (getline(input, strline) && (instance == 0 || iter != instance)) {
      if (RemoveWS) {
        strline = aurostd::RemoveWhiteSpaces(strline, '"');
      }
      if (RemoveComments) {
        strline = aurostd::RemoveComments(strline);
      }
      idx = strline.find(delim);
      if (idx != string::npos) {
        _keyword = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(0, idx));
        if (_keyword == keyword) {
          tokens.push_back(aurostd::RemoveWhiteSpacesFromTheFrontAndBack(strline.substr(idx + delim.length())));
          iter++;
        }
      }
    }
    input.clear();
    input.seekg(0);
    if (tokens.empty() || (uint) aurostd::abs(instance) > tokens.size()) {
      return "";
    }
    stringstream strstream;
    if (instance == 0) {
      for (size_t i = 0; i < tokens.size(); i++) {
        strstream << tokens[i] << endl;
      }
    } else if (instance > 0) {
      strstream << tokens[tokens.size() - 1];
    } else if (instance < 0) {
      const uint i = (uint) aurostd::boundary_conditions_periodic(0, tokens.size() - 1, instance);
      strstream << tokens[i];
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    return strstream.str();
  }

  string kvpair2string(const string& _input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments) {
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    string input = _input;
    if (RemoveWS) {
      input = aurostd::RemoveWhiteSpaces(_input, '"');
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"]" << endl;
    }
    if (input.find(delim) == string::npos || input.find(keyword) == string::npos) {
      return "";
    }
    stringstream strstream;
    vector<string> tokens;
    aurostd::string2vectorstring(input, tokens);
    string _keyword;
    string::size_type idx = 0;
    int iter = 0;
    if (instance > 0) {
      for (size_t i = 0; i < tokens.size(); i++) {
        if (RemoveComments) {
          tokens[i] = aurostd::RemoveComments(tokens[i]);
        }
        idx = tokens[i].find(delim);
        if (idx != string::npos) {
          _keyword = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(0, idx));
          if (_keyword == keyword) {
            iter++;
            if (instance == iter) {
              strstream << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(idx + delim.length()));
              break;
            }
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return strstream.str();
    } else if (instance < 0) {
      for (int i = tokens.size() - 1; i >= 0; i--) {
        if (RemoveComments) {
          tokens[i] = aurostd::RemoveComments(tokens[i]);
        }
        idx = tokens[i].find(delim);
        if (idx != string::npos) {
          _keyword = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(0, idx));
          if (_keyword == keyword) {
            iter--;
            if (instance == iter) {
              strstream << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(idx + delim.length()));
              break;
            }
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return strstream.str();
    } else { // instance==0
      for (size_t i = 0; i < tokens.size(); i++) {
        if (RemoveComments) {
          tokens[i] = aurostd::RemoveComments(tokens[i]);
        }
        idx = tokens[i].find(delim);
        if (idx != string::npos) {
          _keyword = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(0, idx));
          if (_keyword == keyword) {
            strstream << aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tokens[i].substr(idx + delim.length())) << endl;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << "BEGIN [keyword=\"" << keyword << "\"] [delimiter=\"" << delim << "\"] [instance=" << instance << "] [output=\"" << strstream.str() << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
      }
      return strstream.str();
    }
    return "";
  }

  string kvpair2string(const stringstream& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments) {
    return kvpair2string(input.str(), keyword, delim, instance, RemoveWS, RemoveComments);
  }

  template <typename utype> utype kvpair2utype(ifstream& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments) {
    return string2utype<utype>(kvpair2string(input, keyword, delim, instance, RemoveWS, RemoveComments));
  }
  template <typename utype> utype kvpair2utype(const string& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments) {
    return string2utype<utype>(kvpair2string(input, keyword, delim, instance, RemoveWS, RemoveComments));
  }
  template <typename utype> utype kvpair2utype(const stringstream& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments) {
    return kvpair2utype<utype>(input.str(), keyword, delim, instance, RemoveWS, RemoveComments);
  }
#define AST_TEMPLATE(atype) template atype kvpair2utype(ifstream& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype kvpair2utype(const string& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(atype) template atype kvpair2utype(const stringstream& input, const string& keyword, const string& delim, const int instance, bool RemoveWS, bool RemoveComments);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  // ***************************************************************************
  // Function SubStringsPresentExtractString and other
  // ***************************************************************************
  uint substring2strings(ifstream& input, vector<string>& vstringout, const string& strsub1, bool RemoveWS, bool RemoveComments) { // SD20220520
    vstringout = aurostd::string2vectorstring(aurostd::substring2string(input, strsub1, 0, RemoveWS, RemoveComments));
    return vstringout.size();
  }
  uint substring2strings(const string& input, vector<string>& vstringout, const string& strsub1, bool RemoveWS, bool RemoveComments) {
    vstringout = aurostd::string2vectorstring(aurostd::substring2string(input, strsub1, 0, RemoveWS, RemoveComments));
    return vstringout.size();
  }
  uint substring2strings(const stringstream& input, vector<string>& vstringout, const string& strsub1, bool RemoveWS, bool RemoveComments) {
    vstringout = aurostd::string2vectorstring(aurostd::substring2string(input, strsub1, 0, RemoveWS, RemoveComments));
    return vstringout.size();
  }

  /// @brief find strings between consecutive instances of two substrings of a text
  /// @param input input file
  /// @param vstringout container to store the strings that were found
  /// @param strsub1 first substring
  /// @param strsub2 second substring
  /// @param instance (=n) if positive, find the first n instances; if negative, find the last n instances; cannot be zero
  /// @param RemoveWS whether to remove white spaces
  /// @param RemoveComments whether to remove comments
  /// @return number of strings that were found
  ///
  /// @authors
  /// @mod{SG,20240517,created}
  uint substring2strings(ifstream& input, vector<string>& vstringout, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input.rdbuf() << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;
    }
    if (instance == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "instance must be nonzero integer", _INPUT_ILLEGAL_);
    }
    if (aurostd::substring2bool(strsub1, "\n") || aurostd::substring2bool(strsub2, "\n")) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "substring2strings with ifstream input does not support newline characters inside of substrings", _INPUT_ILLEGAL_);
    }
    vstringout.clear();
    size_t start_pos;
    size_t pos1 = 0;
    size_t pos2 = 0;
    string substr;
    string strline;
    bool one_or_two = true; // if true, start searching the line for strsub1; if false, strsub2
    // bool read = false;  // not sure if it's important whether the ifstream lines have been "read"
    input.clear();
    int current_instance = 0;
    if (instance > 0) {
      input.seekg(0);
      while (getline(input, strline)) { // loop over lines in file
        start_pos = 0;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        if (RemoveComments) {
          strline = aurostd::RemoveComments(strline);
        }
        if (RemoveWS) {
          strline = aurostd::RemoveWhiteSpaces(strline, '"');
        }
        while (true) {
          if (one_or_two) {
            // start searching the line for strsub1
            pos1 = strline.find(strsub1, start_pos);
            if (pos1 == string::npos) {
              // if strsub1 not found, move on to the next line and search for strsub1 again
              if (!substr.empty()) {
                substr += "\n";
              }
              break;
            }
            pos2 = strline.find(strsub2, pos1 + strsub1.size());
            if (pos2 == string::npos) {
              // if strsub1 found but not strsub2, move on to the next line and search for strsub2
              one_or_two = false;
              substr = strline.substr(pos1 + strsub1.size(), string::npos) + "\n";
              break;
            }
            // if both strsub1 and strsub2 found in line, increase the instance and if we're not done, continue searching the same line
            substr = strline.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
            substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
            substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
            vstringout.push_back(substr);
            current_instance++;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos2 + strsub2.size();
          } else {
            // start searching the line for strsub2
            pos1 = strline.find(strsub2, 0);
            if (pos1 == string::npos) {
              // if strsub2 not found in line, move on to the next line and search for strsub2 again
              substr += strline + "\n";
              break;
            }
            // if strsub2 found, increase the instance and if we're not done, search the same line for strsub1
            substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
            substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
            vstringout.push_back(substr);
            current_instance++;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos1 + strsub2.size();
            one_or_two = true;
          }
          if (current_instance == instance) {
            break;
          }
        }
        if (current_instance == instance) {
          break;
        }
      }
    } else if (instance < 0) {
      // we have to read the file backwards!
      string thischar;
      input.seekg(-1, std::ios::end);
      const int numchars = input.tellg();
      for (int i = 0; i < numchars; i++) {
        thischar = input.get(); // note that causes current line (std::ios::cur) to move forward by 1
        if (thischar != "\n") {
          // we have to read char by char, and assemble the line string until a newline is encountered
          strline = thischar + strline;
          input.seekg(-2, std::ios::cur);
          continue;
        }
        start_pos = strline.size();
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " strline=" << strline << endl;
        }
        if (RemoveComments) {
          strline = aurostd::RemoveComments(strline);
        }
        if (RemoveWS) {
          strline = aurostd::RemoveWhiteSpaces(strline, '"');
        }
        while (true) {
          if (one_or_two) {
            // start searching the line (backwards) for strsub2
            pos2 = strline.rfind(strsub2, start_pos);
            if (pos2 == string::npos) {
              // if strsub2 not found, move on to the previous line and search for strsub2 again
              if (!substr.empty()) {
                substr = "\n" + substr;
              }
              break;
            }
            pos1 = strline.rfind(strsub1, pos2);
            if (pos1 == string::npos) {
              // if strsub2 found but not strsub1, move on to the previous line and search for strsub1
              one_or_two = false;
              substr = "\n" + strline.substr(0, pos2);
              break;
            }
            // if both strsub2 and strsub1 found in line, decrease the instance and if we're not done, continue searching the same line
            substr = strline.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
            if (!substr.empty()) {
              // In the reverse case, the substr will include strsub1 and strsub2, so we need to remove them
              if (substr.substr(0, strsub1.size()) == strsub1) {
                substr.erase(0, strsub1.size());
              }
              if (substr.substr(substr.size() - strsub2.size(), strsub2.size()) == strsub2) {
                substr.erase(substr.size() - strsub2.size(), strsub2.size());
              }
            }
            substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
            substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
            vstringout.push_back(substr);
            current_instance--;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos1;
          } else {
            // start searching the line (backwards) for strsub1
            pos1 = strline.rfind(strsub1, start_pos);
            if (pos1 == string::npos) {
              // if strsub1 not found in line, move on to the previous line and search for strsub1 again
              substr = "\n" + strline + substr;
              break;
            }
            // if strsub1 found, decrease the instance and if we're not done, search the same line for strsub2
            substr = strline.substr(0, strline.size() - pos1) + substr;
            if (!substr.empty()) {
              // In the reverse case, the substr will include strsub1 and strsub2, so we need to remove them
              if (substr.substr(0, strsub1.size()) == strsub1) {
                substr.erase(0, strsub1.size());
              }
              if (substr.substr(substr.size() - strsub2.size(), strsub2.size()) == strsub2) {
                substr.erase(substr.size() - strsub2.size(), strsub2.size());
              }
            }
            substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
            substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
            vstringout.push_back(substr);
            current_instance--;
            if (current_instance == instance) {
              break;
            }
            substr = "";
            start_pos = pos1;
            one_or_two = true;
          }
          if (current_instance == instance) {
            break;
          }
        }
        if (current_instance == instance) {
          break;
        }
        input.seekg(-2, std::ios::cur); // move backwards 2 lines
        strline = "";
      }
      std::reverse(vstringout.begin(), vstringout.end());
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    const uint count = vstringout.size();
    return count;
  }

  /// @brief find strings between consecutive instances of two substrings of a text
  /// @param input input text
  /// @param vstringout container to store the strings that were found
  /// @param strsub1 first substring
  /// @param strsub2 second substring
  /// @param instance (=n) if positive, find the first n instances; if negative, find the last n instances; cannot be zero
  /// @param RemoveWS whether to remove white spaces
  /// @param RemoveComments whether to remove comments
  /// @return number of strings that were found
  ///
  /// @authors
  /// @mod{SG,20240517,created}
  uint substring2strings(const std::string& input, vector<std::string>& vstringout, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) {
    const bool LDEBUG = false;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "BEGIN [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "[input=\"" << input << "\"] [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"]" << endl;
    }
    if (instance == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "instance must be nonzero integer", _INPUT_ILLEGAL_);
    }
    vstringout.clear();
    std::string input_mod = input;
    if (RemoveComments) {
      input_mod = aurostd::RemoveComments(input_mod);
    }
    if (RemoveWS) {
      input_mod = aurostd::RemoveWhiteSpaces(input_mod, '"');
    }
    size_t start_pos = 0;
    size_t pos1 = 0;
    size_t pos2 = 0;
    string substr;
    int this_instance;
    if (instance >= 0) {
      this_instance = 1;
      start_pos = 0;
      while (true) {
        pos1 = input_mod.find(strsub1, start_pos);
        pos2 = input_mod.find(strsub2, pos1 + strsub1.size());
        if (pos1 == string::npos or pos2 == string::npos) {
          break;
        }
        substr = input_mod.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
        substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
        substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
        vstringout.push_back(substr);
        if (this_instance == instance) {
          break;
        }
        this_instance++;
        start_pos = pos2 + strsub2.size();
      }
    } else if (instance < 0) {
      this_instance = -1;
      start_pos = input_mod.size();
      while (true) {
        pos2 = input_mod.rfind(strsub2, start_pos);
        pos1 = input_mod.rfind(strsub1, pos2);
        if (pos1 == string::npos or pos2 == string::npos) {
          break;
        }
        substr = input_mod.substr(pos1 + strsub1.size(), pos2 - pos1 - strsub1.size());
        substr = RemoveWhiteSpacesFromTheFrontAndBack(substr);
        substr = RemoveCharacterFromTheFrontAndBack(substr, '\n');
        vstringout.push_back(substr);
        if (this_instance == instance) {
          break;
        }
        this_instance--;
        start_pos = pos1;
      }
      std::reverse(vstringout.begin(), vstringout.end());
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "END [substring=\"" << strsub1 << "\"] [substring=\"" << strsub2 << "\"] [instance=" << instance << "] [RemoveWS=" << RemoveWS << "]" << endl;
    }
    const uint count = vstringout.size();
    return count;
  }

  /// @brief find strings between consecutive instances of two substrings of a text
  /// @param input input text
  /// @param vstringout container to store the strings that were found
  /// @param strsub1 first substring
  /// @param strsub2 second substring
  /// @param instance (=n) if positive, first n instances; if negative, last n instances; cannot be zero
  /// @param RemoveWS whether to remove white spaces
  /// @param RemoveComments whether to remove comments
  /// @return number of strings that were found
  ///
  /// @authors
  /// @mod{CO,20190617,created}
  /// @mod{SG,20240517,added documentation and removed trim_edges}
  uint substring2strings(const stringstream& input, vector<string>& vstringout, const string& strsub1, const string& strsub2, const int instance, bool RemoveWS, bool RemoveComments) { // CO20230502
    return aurostd::substring2strings(input.str(), vstringout, strsub1, strsub2, instance, RemoveWS, RemoveComments); // CO20230502
  }

  template <typename utype> uint substring2utypes(ifstream& input, vector<utype>& vutypeout, const string& strsub1, bool RemoveWS, bool RemoveComments) { // SD202205020
    vector<string> vstringout;
    aurostd::substring2strings(input, vstringout, strsub1, RemoveWS, RemoveComments);
    vutypeout = aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template <typename utype>
  uint substring2utypes(const string& input, vector<utype>& vutypeout, const string& strsub1, bool RemoveWS, bool RemoveComments) { // CO20210315 - cleaned up //SD202205020 - added utype
    vector<string> vstringout;
    aurostd::substring2strings(input, vstringout, strsub1, RemoveWS, RemoveComments);
    vutypeout = aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template <typename utype>
  uint substring2utypes(const stringstream& input, vector<utype>& vutypeout, const string& strsub1, bool RemoveWS, bool RemoveComments) { // CO20210315 - cleaned up //SD202205020 - added utype
    return substring2utypes<utype>(input.str(), vutypeout, strsub1, RemoveWS, RemoveComments);
  }

  template <typename utype> uint substring2utypes(ifstream& input, vector<utype>& vutypeout, const string& strsub1, const string& strsub2, bool RemoveWS, bool RemoveComments) { // SD202205020
    vector<string> vstringout;
    aurostd::substring2strings(input, vstringout, strsub1, strsub2, RemoveWS, RemoveComments);
    vutypeout = aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template <typename utype> uint substring2utypes(const string& input, vector<utype>& vutypeout, const string& strsub1, const string& strsub2, bool RemoveWS, bool RemoveComments) { // SD202205020 - added utype
    vector<string> vstringout;
    aurostd::substring2strings(input, vstringout, strsub1, strsub2, RemoveWS, RemoveComments);
    vutypeout = aurostd::vectorstring2vectorutype<utype>(vstringout);
    return vutypeout.size();
  }
  template <typename utype> uint substring2utypes(const stringstream& input, vector<utype>& vutypeout, const string& strsub1, const string& strsub2, bool RemoveWS, bool RemoveComments) { // SD202205020
    return substring2utypes<utype>(input.str(), vutypeout, strsub1, strsub2, RemoveWS, RemoveComments);
  }
} // namespace aurostd

// ***************************************************************************
// FUNCTION HTML LATEX TXT

namespace // anonymous-namespace
{
  /// @brief convert between latex and HTML accents
  /// @authors
  /// @mod{HE,20260327,created}
  /// @see
  /// @xlink{http://en.wikibooks.org/wiki/LaTeX/Accents}
  /// @xlink{https://en.wikipedia.org/wiki/List_of_XML_and_HTML_character_entity_references#List_of_character_entity_references_in_HTML}
  void accent_convert(std::string& work_str, const std::string& direction) {
    const std::vector<std::pair<std::string, std::string>> accent_list({
        { "\\`", "grave"}, // ò grave accent
        { "\\'", "acute"}, // ó acute accent
        { "\\^",  "circ"}, // ô circumflex
        {"\\\"",   "uml"}, // ö umlaut
        { "\\H", "dblac"}, // ő	long Hungarian umlaut (double acute)
        { "\\~", "tilde"}, // õ tilde
        { "\\c", "cedil"}, // ç cedilla
        { "\\k",  "ogon"}, // ą	ogonek
        { "\\=",  "macr"}, // ō	macron accent (a bar over the letter)
        { "\\.",   "dot"}, // ȯ dot over the letter
        { "\\r",  "ring"}, // å ring over the letter
        { "\\u", "breve"}, // ŏ breve over the letter
        { "\\v", "caron"}, // š caron/háček ("v") over the letter
    });
    const std::vector<std::string> letters({"a", "A", "c", "C", "e", "E", "i", "I", "o", "O", "u", "U", "n", "N", "z", "Z"});

    for (const auto& [latex, html] : accent_list) {
      for (const string& letter : letters) {
        if (direction == "latex2html") {
          aurostd::StringSubstInPlace(work_str, latex + "{" + letter + "}", "&" + letter + html + ";");
        } else if (direction == "html2latex") {
          aurostd::StringSubstInPlace(work_str, "&" + letter + html + ";", latex + "{" + letter + "}"); // umlaut
        }
      }
    }
  }
} // namespace

namespace aurostd {

  // http://www.w3schools.com/html/html_entities.asp
  // http://en.wikibooks.org/wiki/LaTeX/Accents

  // ME20200921 - Replaces HTML special characters with the correct entity name
  string text2html(const string& str) {
    string out = str;
    // Ampersand must come first since it is in the entity name!
    aurostd::StringSubstInPlace(out, "&", "&amp;");
    aurostd::StringSubstInPlace(out, "<", "&lt;");
    aurostd::StringSubstInPlace(out, ">", "&gt;");
    aurostd::StringSubstInPlace(out, "\"", "&quot;");
    aurostd::StringSubstInPlace(out, "'", "&apos;");
    return out;
  }

  /// @brief convert HTML characters to their latex variant
  /// @authors
  /// @mod{HE,20241031,compacting and bug fixes}
  /// @mod{HE,20260327,use unified accent function}
  string html2latex(const string& str) {
    string out = str;
    const std::vector<std::pair<std::string, std::string>> simple_replacements = {
        {       "_",        "\\_"},
        {   "<sub>",        "$_{"},
        {  "</sub>",         "}$"},
        {   "<sup>",        "$^{"},
        {  "</sup>",         "}$"},
        {     "<i>",  "\\textit{"},
        {    "</i>",          "}"},
        {     "<b>",  "\\textbf{"},
        {    "</b>",          "}"},
        { "<blink>",  "\\textbf{"},
        {"</blink>",          "}"},
        {    "MgB2",    "MgB$_2$"},
        {  "Csanyi",  "Cs\\'anyi"},
        {   "Pólya", "P\\'{o}lya"}  // Legacy special cases
    };
    for (const auto& [old, replacements] : simple_replacements) {
      aurostd::StringSubstInPlace(out, old, replacements);
    }

    if (!aurostd::substring2bool(out, "Rosenbrock")) {
      aurostd::StringSubstInPlace(out, "Rosen", "Ros\\'en");
    }

    // Accents
    accent_convert(out, "html2latex");

    aurostd::StringSubstInPlace(out, "&lstrok;", "\\l{}");
    aurostd::StringSubstInPlace(out, "&Lstrok;", "\\L{}"); // strok
    aurostd::StringSubstInPlace(out, "&oslash;", "{\\o{}}");
    aurostd::StringSubstInPlace(out, "&Oslash;", "{\\O{}}"); // slash

    // Greek symbols
    for (string letter :
         {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega"}) {
      aurostd::StringSubstInPlace(out, "&" + letter + ";", "$\\" + letter + "$");
      letter[0] = std::toupper(letter[0]);
      aurostd::StringSubstInPlace(out, "&" + letter + ";", "$\\" + letter + "$");
    }
    aurostd::StringSubstInPlace(out, "&thetasym", "$\\vartheta$");

    // FINAL
    aurostd::StringSubstInPlace(out, "&", "\\&");

    return out;
  }

  string html2txt(const string& str) {
    string out = str;
    aurostd::StringSubstInPlace(out, "<sub>", "");
    aurostd::StringSubstInPlace(out, "</sub>", "");
    aurostd::StringSubstInPlace(out, "<i>", "");
    aurostd::StringSubstInPlace(out, "</i>", "");
    aurostd::StringSubstInPlace(out, "<b>", "");
    aurostd::StringSubstInPlace(out, "</b>", "");
    aurostd::StringSubstInPlace(out, "MgB2", "MgB2");
    aurostd::StringSubstInPlace(out, "&", "&");
    aurostd::StringSubstInPlace(out, "_", "");
    aurostd::StringSubstInPlace(out, "\\", "");
    return out;
  }

  // ***************************************************************************
  // Function aurostd::string2latex
  // ***************************************************************************
  string string2latex(const string& str) {
    string out = str;
    aurostd::StringSubstInPlace(out, "_pv", "_{pv}");
    aurostd::StringSubstInPlace(out, "_sv", "_{sv}");
    aurostd::StringSubstInPlace(out, "_h", "_{h}");
    aurostd::StringSubstInPlace(out, "_d", "_{d}");
    aurostd::StringSubstInPlace(out, "_s", "_{s}");
    aurostd::StringSubstInPlace(out, "_1", "_{1}");
    aurostd::StringSubstInPlace(out, "_2", "_{2}");
    aurostd::StringSubstInPlace(out, "_3", "_{3}");
    return out;
  }

  /// @brief convert latex characters to their HTML varient
  /// @authors
  /// @mod{HE,20241031,compacting and bug fixes}
  /// @mod{HE,20260327,use unified accent function}
  std::string latex2html(const std::string& str) {
    string out = str;

    accent_convert(out, "latex2html");

    aurostd::StringSubstInPlace(out, "\\l{}", "&lstrok;");
    aurostd::StringSubstInPlace(out, "\\L{}", "&Lstrok;"); // strok
    aurostd::StringSubstInPlace(out, "{\\o{}}", "&oslash;");
    aurostd::StringSubstInPlace(out, "{\\O{}}", "&Oslash;"); // slash

    // Greek symbols
    for (string letter :
         {"alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega"}) {
      aurostd::StringSubstInPlace(out, "$\\" + letter + "$", "&" + letter + ";");
      letter[0] = std::toupper(letter[0]);
      aurostd::StringSubstInPlace(out, "$\\" + letter + "$", "&" + letter + ";");
    }
    aurostd::StringSubstInPlace(out, "$\\vartheta$", "&thetasym");

    // subscript replace with regex
    const std::regex sub_regex(R"(\$_\{(.*?)\}\$)");
    out = std::regex_replace(out, sub_regex, "<sub>$1</sub>");
    const std::regex sub_short_regex(R"(\$_(.)\$)");
    out = std::regex_replace(out, sub_short_regex, "<sub>$1</sub>");

    // superscipt replace with regex
    const std::regex sup_regex(R"(\$\^\{(.*?)\}\$)");
    out = std::regex_replace(out, sub_regex, "<sup>$1</sup>");
    const std::regex sup_short_regex(R"(\$\^(.)\$)");
    out = std::regex_replace(out, sup_short_regex, "<sup>$1</sup>");

    // italics replace with regex
    const std::regex italics_regex(R"(\textit{(.*?)\})");
    out = std::regex_replace(out, sub_regex, "<i>$1</i>");

    // bold replace with regex
    const std::regex bold_regex(R"(\textbf{(.*?)\})");
    out = std::regex_replace(out, sub_regex, "<b>$1</b>");

    // Final
    aurostd::StringSubstInPlace(out, "\\&", "&");

    return out;
  }

  string latex2txt(const string& str) {
    string out = str;
    aurostd::StringSubstInPlace(out, "\\&", "&");
    aurostd::StringSubstInPlace(out, "MgB$_2$", "MgB2");
    aurostd::StringSubstInPlace(out, "<sub>", "");
    aurostd::StringSubstInPlace(out, "</sub>", "");
    aurostd::StringSubstInPlace(out, "<i>", "");
    aurostd::StringSubstInPlace(out, "</i>", "");
    aurostd::StringSubstInPlace(out, "<b>", "");
    aurostd::StringSubstInPlace(out, "</b>", "");
    return out;
  }

  // CO20190419
  string fixStringLatex(const string& input, bool double_back_slash, bool symmetry_string) {
    // deals with special characters for LaTeX, like some characters in prototype
    // see http://tex.stackexchange.com/questions/34580/escape-character-in-latex
    // double_back_slash was needed SOMETIMES for gnuplot output, as one backslash
    string output;
    vector<char> problem_characters;
    problem_characters.push_back('&');
    problem_characters.push_back('%');
    problem_characters.push_back('$');
    problem_characters.push_back('#');
    if (!symmetry_string) {
      problem_characters.push_back('_');
      problem_characters.push_back('{');
      problem_characters.push_back('}');
    }
    problem_characters.push_back('~'); // different fix
    problem_characters.push_back('^'); // different fix
    string solution_string;
    solution_string = "\\\\"; // has to be string, \\ char does not work
    bool found_escaped_char;
    bool found_hyphen_symmetry = false;
    bool solved_hyphen_symmetry = false;
    for (uint i = 0; i < input.length(); i++) {
      // we first enter this loop because symmetry_string and input[i]=='-'
      // second enter loop because symmetry_string and found_hyphen_symmetry
      if (symmetry_string && (input[i] == '-' || found_hyphen_symmetry)) {
        if (!found_hyphen_symmetry) {
          // first enter loop, come here
          found_hyphen_symmetry = true;
          output.append((double_back_slash ? string("\\") : string("")) + string("\\overline{"));
          // very important, we don't want to add hyphen, just replace
          // with overline, so continue
          continue;
        } else {
          // second enter loop, do nothing but turn this flag on
          // allow us to add input[i]
          found_hyphen_symmetry = false;
          solved_hyphen_symmetry = true;
        }
      } else {
        if (symmetry_string && solved_hyphen_symmetry) {
          // last step of symmetry_string fix, but we have to do this in part of
          // the loop to allow for next character to be identified as problem
          // character as well
          output.append(1, '}');
          solved_hyphen_symmetry = false;
        }
        // go through all problem characters
        for (size_t j = 0, fl_size_j = problem_characters.size(); j < fl_size_j; j++) {
          if (input[i] == problem_characters[j]) {
            if (double_back_slash) {
              // if we find one, but it has double backslash, leave alone
              // doesn't matter what it is, if it has double backslash it's good
              // if we find one, but it only has single backslash, add one
              if (i && i - 1 && input[i - 1] == '\\' && input[i - 2] == '\\') {
                break;
              } else if (i && input[i - 1] == '\\') {
                output.append(1, '\\'); // just add one
                break;
              }
              // if we find one, give two backslashes
              output.append("\\\\");
              break;
            } else {
              // if we find one, but it has single backslash, leave alone
              // doesn't matter what it is, if it has single backslash it's good
              // if we find one, give single backslash
              if (i && input[i - 1] == '\\') {
                break;
              }
              output.append(1, '\\');
              break;
            }
          }
        }
        // we also have to add {} for these characters
        if (input[i] == '~' || input[i] == '^') {
          output.append("{}");
        }
        found_escaped_char = false;
        if (input[i] == '\\') {
          for (size_t j = 0, fl_size_j = problem_characters.size(); j < fl_size_j; j++) {
            // the only way this works if it's serving as an escape for a character
            // don't worry about double backslash here, we get to that when we find
            // the actual character
            if (i != (input.length() - 1) && input[i + 1] == problem_characters[j]) {
              found_escaped_char = true;
              break; // doesn't matter what it is, if it has backslash it's good
            }
          }
          // this is a problem, no way around it--we cannot output single backslash
          if (!found_escaped_char) {
            stringstream message;
            message << "Extraneous backslash found in \"" << input << "\" which may cause problems for LaTeX/gnuplot";
            cerr << __AFLOW_FUNC__ << " ERROR - " << message.str() << endl;
            return input;
          }
        }
      }
      // add in character from input
      output.append(1, input[i]);
    }
    return output;
  }
} // namespace aurostd

// ***************************************************************************
// SORT WORLD
// ----------------------------------------------------------------------------
// sort for vector (starting from xvector)

namespace aurostd {
  template <class utype1> // function quicksort
  void sort(vector<utype1>& arr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    aurostd::sort(xarr);
    arr = aurostd::xvector2vector(xarr);
  }
#define AST_TEMPLATE(utype) template void sort(vector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype1, class utype2> // function quicksort
  void sort(vector<utype1>& arr, vector<utype2>& brr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    xvector<utype2> xbrr = aurostd::vector2xvector(brr);
    aurostd::sort2(xarr.rows, xarr, xbrr);
    arr = aurostd::xvector2vector(xarr);
    brr = aurostd::xvector2vector(xbrr);
  }

#define AST_TEMPLATE(utype1, utype2) template void sort(vector<utype1>&, vector<utype2>&);
  AST_GEN_2(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype1, class utype2> // function quicksort //CO20200915
  void sort(deque<utype1>& arr, deque<utype2>& brr) {
    xvector<utype1> xarr(arr.size());
    xvector<utype2> xbrr(brr.size());
    for (size_t i = 0; i < arr.size(); i++) {
      xarr[i + 1] = arr[i];
    }
    for (size_t i = 0; i < brr.size(); i++) {
      xbrr[i + 1] = brr[i];
    }
    aurostd::sort2(xarr.rows, xarr, xbrr);
    // aurostd::sort2(xarr,xbrr);
    arr.clear();
    brr.clear();
    for (int i = 0; i < xarr.rows; i++) {
      arr.push_back(xarr[i + 1]);
      brr.push_back(xbrr[i + 1]);
    }
  }
#define AST_TEMPLATE(utype1, utype2) template void sort(deque<utype1>&, deque<utype2>&);
  AST_GEN_2(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype1, class utype2, class utype3> // function quicksort
  void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    xvector<utype2> xbrr = aurostd::vector2xvector(brr);
    xvector<utype3> xcrr = aurostd::vector2xvector(crr);
    aurostd::sort3(xarr.rows, xarr, xbrr, xcrr);
    arr = aurostd::xvector2vector(xarr);
    brr = aurostd::xvector2vector(xbrr);
    crr = aurostd::xvector2vector(xcrr);
  }

#define AST_TEMPLATE(utype) template void sort(vector<utype>&, vector<utype>&, vector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype1, class utype2, class utype3, class utype4> // function quicksort
  void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr, vector<utype4>& drr) {
    xvector<utype1> xarr = aurostd::vector2xvector(arr);
    xvector<utype2> xbrr = aurostd::vector2xvector(brr);
    xvector<utype3> xcrr = aurostd::vector2xvector(crr);
    xvector<utype4> xdrr = aurostd::vector2xvector(drr);
    aurostd::sort4(xarr.rows, xarr, xbrr, xcrr, xdrr);
    arr = aurostd::xvector2vector(xarr);
    brr = aurostd::xvector2vector(xbrr);
    crr = aurostd::xvector2vector(xcrr);
    drr = aurostd::xvector2vector(xdrr);
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// sort for vector of strings
namespace aurostd {
  void sort(vector<string>& arg) {
    sort(arg.begin(), arg.end(), aurostd::_sort_string_());
  }
  void sort(deque<string>& arg) {
    std::sort(arg.begin(), arg.end(), aurostd::_sort_string_());
  }
  void rsort(vector<string>& arg) {
    std::reverse(arg.begin(), arg.end()); //,aurostd::_sort_string_());
  }
  void rsort(deque<string>& arg) {
    std::reverse(arg.begin(), arg.end()); //,aurostd::_sort_string_());
  }
} // namespace aurostd

// sort_remove_duplicates for vector of strings
namespace aurostd {
  void sort_remove_duplicates(vector<string>& arg) {
    sort(arg.begin(), arg.end(), aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(), arg.end()), arg.end());
  }
  void sort_remove_duplicates(deque<string>& arg) {
    std::sort(arg.begin(), arg.end(), aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(), arg.end()), arg.end());
  }
  void rsort_remove_duplicates(vector<string>& arg) {
    std::reverse(arg.begin(), arg.end()); //,aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(), arg.end()), arg.end());
  }
  void rsort_remove_duplicates(deque<string>& arg) {
    std::reverse(arg.begin(), arg.end()); //,aurostd::_sort_string_());
    arg.erase(std::unique(arg.begin(), arg.end()), arg.end());
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_int
namespace aurostd {
  void sort(vector<string>& varg1, vector<int>& varg2) {
    vector<aurostd::_string_int_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_int_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
  void sort(deque<string>& varg1, deque<int>& varg2) {
    deque<aurostd::_string_int_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_int_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_double
namespace aurostd {
  void sort(vector<string>& varg1, vector<double>& varg2) {
    vector<aurostd::_string_double_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_double_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
  void sort(deque<string>& varg1, deque<double>& varg2) {
    deque<aurostd::_string_double_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_double_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string
namespace aurostd {
  void sort(vector<string>& varg1, vector<string>& varg2) {
    vector<aurostd::_string_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
  void sort(deque<string>& varg1, deque<string>& varg2) {
    deque<aurostd::_string_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of double_int
// HERE THEY ARE

namespace aurostd {
  void sort(vector<double>& varg1, vector<int>& varg2) {
    vector<aurostd::_double_int_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_double_int_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
  void sort(deque<double>& varg1, deque<int>& varg2) {
    deque<aurostd::_double_int_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_double_int_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of double_double
namespace aurostd {
  void sort(vector<double>& varg1, vector<double>& varg2) {
    vector<aurostd::_double_double_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_double_double_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
  void sort(deque<double>& varg1, deque<double>& varg2) {
    deque<aurostd::_double_double_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_double_double_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of double_string
namespace aurostd {
  void sort(vector<double>& varg1, vector<string>& varg2) {
    vector<aurostd::_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
  void sort(deque<double>& varg1, deque<string>& varg2) {
    deque<aurostd::_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
    }
    sort(vv.begin(), vv.end(), _sort_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_int_string
namespace aurostd {
  void sort(vector<string>& varg1, vector<int>& varg2, vector<string>& varg3) {
    vector<aurostd::_string_int_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_int_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
    }
  }
  void sort(deque<string>& varg1, deque<int>& varg2, deque<string>& varg3) {
    deque<aurostd::_string_int_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_int_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_double_string
namespace aurostd {
  void sort(vector<string>& varg1, vector<double>& varg2, vector<string>& varg3) {
    vector<aurostd::_string_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
    }
  }
  void sort(deque<string>& varg1, deque<double>& varg2, deque<string>& varg3) {
    deque<aurostd::_string_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_string
namespace aurostd {
  void sort(vector<string>& varg1, vector<string>& varg2, vector<string>& varg3) {
    vector<aurostd::_string_string_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
    }
  }
  void sort(deque<string>& varg1, deque<string>& varg2, deque<string>& varg3) {
    deque<aurostd::_string_string_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_double_string
namespace aurostd {
  void sort(vector<string>& varg1, vector<string>& varg2, vector<double>& varg3, vector<string>& varg4) {
    vector<aurostd::_string_string_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
      vv[i].arg4 = varg4[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
      varg4[i] = vv[i].arg4;
    }
  }
  void sort(deque<string>& varg1, deque<string>& varg2, deque<double>& varg3, deque<string>& varg4) {
    deque<aurostd::_string_string_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
      vv[i].arg4 = varg4[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
      varg4[i] = vv[i].arg4;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// sort for vector/deque of string_string_double_double_string
namespace aurostd {
  void sort(vector<string>& varg1, vector<string>& varg2, vector<double>& varg3, vector<double>& varg4, vector<string>& varg5) {
    vector<aurostd::_string_string_double_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
      vv[i].arg4 = varg4[i];
      vv[i].arg5 = varg5[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_double_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
      varg4[i] = vv[i].arg4;
      varg5[i] = vv[i].arg5;
    }
  }
  void sort(deque<string>& varg1, deque<string>& varg2, deque<double>& varg3, deque<double>& varg4, deque<string>& varg5) {
    deque<aurostd::_string_string_double_double_string_> vv(varg1.size());
    for (size_t i = 0; i < varg1.size(); i++) {
      vv[i].arg1 = varg1[i];
      vv[i].arg2 = varg2[i];
      vv[i].arg3 = varg3[i];
      vv[i].arg4 = varg4[i];
      vv[i].arg5 = varg5[i];
    }
    sort(vv.begin(), vv.end(), _sort_string_string_double_double_string_());
    for (size_t i = 0; i < varg1.size(); i++) {
      varg1[i] = vv[i].arg1;
      varg2[i] = vv[i].arg2;
      varg3[i] = vv[i].arg3;
      varg4[i] = vv[i].arg4;
      varg5[i] = vv[i].arg5;
    }
  }
} // namespace aurostd

// ----------------------------------------------------------------------------
// reorder vector //CO20221111
namespace aurostd {
  template <class utype> // function quicksort
  void reorder(vector<utype>& vec, vector<uint>& vorder, uint mode) { // CO20221111
    // algorithms and discussion from here: https://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices (very good!)
    // solution by chmike
    // reorder a vector given input indices
    // there are two ways this can be done depending on what is inside vorder
    // input: vec={7,5,9,6}; vorder={1,3,0,2}
    //
    // mode 1: ``draw the elements of vector from the position of the indices''
    // result: {5,6,7,9}
    // NOTE: this is the default mode
    //
    // mode 2: ``move elements of vector to the position of the indices''
    // result: {9,7,6,5}
    // NOTE: this can also be accomplished with aurostd::sort(vorder,vec) but it requires vec to be C++ type or string
    // this function seems to run faster than aurostd::sort() as well
    uint i = 0;
    uint j = 0;
    if (mode == 1) {
      for (i = 0; i < vec.size() - 1; i++) {
        if (vorder[i] == i) {
          continue;
        }
        for (j = i + 1; j < vorder.size(); j++) {
          if (vorder[j] == i) {
            break;
          }
        }
        std::iter_swap(vec.begin() + i, vec.begin() + vorder[i]);
        std::iter_swap(vorder.begin() + i, vorder.begin() + j);
      }
      return;
    } else if (mode == 2) {
      uint alt = 0;
      // for all elements to put in place
      for (i = 0; i < vec.size() - 1; ++i) {
        // while the element i is not yet in place
        while (i != vorder[i]) {
          // swap it with the element at its final place
          alt = vorder[i];
          std::iter_swap(vec.begin() + i, vec.begin() + alt);
          std::iter_swap(vorder.begin() + i, vorder.begin() + alt);
        }
      }
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown mode", _INPUT_ILLEGAL_);
    }
  }
#define AST_TEMPLATE(utype) template void reorder(vector<utype>&, vector<uint>&, uint);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
  template void reorder(vector<aflowlib::_aflowlib_entry>&, vector<uint>&, uint);
} // namespace aurostd

// ***************************************************************************
// Function some statistical stuff
// combinations
// ***************************************************************************
namespace aurostd {
  template <class utype> utype combinations(utype n, utype k) { // http://en.wikipedia.org/wiki/Combination // C^n_k=n!/k!(n-k)!   hard to calculate
    double cnk = 1.0;
    for (utype i = 0; i <= k - 1; i++) {
      cnk = cnk * (n - i) / (k - i);
    }
    return (utype) cnk;
  }
#define AST_TEMPLATE(utype) template utype combinations(utype n, utype k);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> utype Cnk(utype n, utype k) {
    return combinations(n, k);
  } // http://en.wikipedia.org/wiki/Combination
} // namespace aurostd

// ***************************************************************************
// aurostd::ShiftFirstColumn(const vector<vector<double> >& a, const double& value)
// ***************************************************************************
namespace aurostd {
  vector<vector<double>> ShiftFirstColumn(const vector<vector<double>>& vva, const double& value) {
    // change value in the first column (usually menas energy in DOS)
    vector<vector<double>> vvb = vva;
    for (size_t i = 0; i < vvb.size(); i++) {
      vvb[i].at(0) = vvb[i].at(0) - value;
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& value)
// ***************************************************************************
namespace aurostd {
  vector<vector<double>> ShrinkValuesExceptFirstColumn(const vector<vector<double>>& vva, const double& Fi) {
    // shrink Fis (usually means DOS Fis in DOS); Fi means probability
    vector<vector<double>> vvb = vva;
    for (size_t i = 0; i < vvb.size(); i++) {
      for (size_t j = 1; j < vvb[i].size(); j++) {
        vvb[i][j] *= Fi;
      }
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// vector<vector<double> > aurostd::NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<vector<double> >& vFi)
// ***************************************************************************
namespace aurostd {
  vector<vector<double>> NormalizeAndSum3DVector(const vector<vector<vector<double>>>& vvva, const vector<double>& vFi) {
    // normalize DOS and sum
    if (vvva.size() != vFi.size()) {
      const string message = "Vector sizes are not equal.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    vector<vector<double>> vvb;
    vector<vector<double>> vv_tmp;
    vector<vector<double>> vv_tmp_shrinked;
    vector<vector<vector<double>>> vvvc;
    double Fi;
    for (size_t i = 0; i < vvva.size(); i++) {
      vv_tmp = vvva[i];
      Fi = vFi[i];
      vv_tmp_shrinked = aurostd::ShrinkValuesExceptFirstColumn(vv_tmp, Fi);
      vvvc.push_back(vv_tmp_shrinked);
    }
    vvb = aurostd::Sum3DVectorAndReduce2D(vvvc);
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva)
// ***************************************************************************
namespace aurostd {
  vector<vector<double>> Sum3DVectorAndReduce2D(const vector<vector<vector<double>>>& vvva) {
    // The first column will not change! (For example, PDOS into TOTALPDOS)
    vector<vector<double>> vvtmp;
    vector<vector<double>> vv_sum;
    vv_sum = vvva.at(0);
    for (size_t i = 1; i < vvva.size(); i++) {
      vvtmp = vvva[i];
      vv_sum = aurostd::Sum2DVectorExceptFirstColumn(vv_sum, vvtmp);
    }
    return vv_sum;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva)
// ***************************************************************************
namespace aurostd {
  vector<vector<double>> Sum2DVectorExceptFirstColumn(const vector<vector<double>>& vva, const vector<vector<double>>& vvb) {
    if ((vva.size() != vvb.size()) && (vva.at(0).size() != vvb.at(0).size())) {
      const string message = "Vector sizes are not equal.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    vector<vector<double>> vv_sum;
    vv_sum.resize(vva.size());
    for (size_t i = 0; i < vva.size(); i++) {
      const int N = vva[i].size();
      vv_sum[i].resize(N);
    }

    for (size_t i = 0; i < vva.size(); i++) {
      vv_sum[i][0] = vva[i].at(0);
      for (size_t j = 1; j < vva[i].size(); j++) {
        vv_sum[i][j] = vva[i][j] + vvb[i][j];
      }
    }
    return vv_sum;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::ReduceVector(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  vector<vector<double>> ReduceVector(const vector<vector<double>>& vva, const int& n) {
    // Pick up the first (begin from 0) and the nth column of 2D vector
    vector<vector<double>> vvb;
    vvb.clear();
    vector<double> vtmp;
    for (size_t i = 0; i < vva.size(); i++) {
      vtmp.clear();
      vtmp.push_back(vva[i].at(0));
      vtmp.push_back(vva[i].at(n));
      vvb.push_back(vtmp);
    }
    return vvb;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  double CalculateIntegrate(const vector<vector<double>>& vva, const int& n) {
    // Calculate integration of vva, the 0st column is x0, x1..., the n column is y1, y2 ...
    // begin from 0
    const vector<vector<double>> vvb = aurostd::ReduceVector(vva, n);
    return aurostd::CalculateIntegrate(vvb);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  double CalculateIntegrate(const vector<vector<double>>& vva, const int& n, const double& Emin, const double& Emax) {
    // Calculate integration of vva, the 0st column is x0, x1..., the n column is y1, y2 ...
    // begin from 0
    const vector<vector<double>> vvb = aurostd::ReduceVector(vva, n);
    return aurostd::CalculateIntegrate(vvb, Emin, Emax);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  double CalculateIntegrate(const vector<vector<double>>& vva) {
    const double Emin = -100;
    const double Emax = 0.0; // default setting
    return aurostd::CalculateIntegrate(vva, Emin, Emax);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CalculateIntegrate(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  double CalculateIntegrate(const vector<vector<double>>& vva, const double& Emin, const double& Emax) {
    // Integral function
    // format of vva: x0, y0; x1, y1; x2, y2
    double integral_result = 0.0;
    double area_tmp = 0.0;
    double xbeg;
    double xend;
    double ybeg;
    double yend;
    for (size_t i = 0; i < vva.size() - 1; i++) {
      xbeg = vva[i].at(0);
      xend = vva.at(i + 1).at(0);
      ybeg = vva[i].at(1);
      yend = vva.at(i + 1).at(1);
      if (xbeg >= Emin && xend <= Emax) {
        area_tmp = 0.5 * (ybeg + yend) * (xend - xbeg);
        integral_result += area_tmp;
      }
    }
    return integral_result;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2string(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  string vector2string(const vector<vector<double>>& vva) {
    stringstream ss_vva;
    aurostd::StringstreamClean(ss_vva);
    ss_vva << std::scientific;
    for (size_t i = 0; i < vva.size(); i++) {
      for (size_t j = 0; j < vva[i].size(); j++) {
        ss_vva << vva[i][j] << "   ";
      }
      ss_vva << endl;
    }
    return ss_vva.str();
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2deque(const vector<utype>& vin)
// ***************************************************************************
// CO20181226
namespace aurostd {
  template <class utype> deque<utype> vector2deque(const vector<utype>& vin) {
    deque<utype> dout;
    for (size_t i = 0; i < vin.size(); i++) {
      dout.push_back(vin[i]);
    }
    return dout;
  }
#define AST_TEMPLATE(utype) template deque<utype> vector2deque(const vector<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE
} // namespace aurostd

// ***************************************************************************
// aurostd::vector2deque(const vector<utype>& vin)
// ***************************************************************************
// CO20181226
namespace aurostd {
  template <class utype> vector<utype> deque2vector(const deque<utype>& din) {
    vector<utype> vout;
    for (size_t i = 0; i < din.size(); i++) {
      vout.push_back(din[i]);
    }
    return vout;
  }
#define AST_TEMPLATE(utype) template vector<utype> deque2vector(const deque<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE

} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva)
// ***************************************************************************
namespace aurostd {
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>>& vva) {
    const double min = -10; // default
    const double max = 10;
    return aurostd::FindMaxIn2DvectorExcept1stColumn(vva, min, max);
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>& vva, const double& min, const double& max)
// ***************************************************************************
namespace aurostd {
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double>>& vva, const double& min, const double& max) {
    double max_value = 0.0;
    for (size_t i = 0; i < vva.size(); i++) {
      const double E_tmp = vva[i].at(0);
      if (E_tmp >= min && E_tmp <= max) {
        for (size_t j = 1; j < vva[i].size(); j++) {
          const double db_tmp = vva[i][j];
          if (abs(db_tmp) > max_value) {
            max_value = abs(db_tmp);
          }
        }
      }
    }
    return max_value;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max)
// ***************************************************************************
namespace aurostd {
  double FindMaxInTDOS(const vector<vector<double>>& vva, const double& min, const double& max) {
    double max_value = 0.0;
    for (size_t i = 0; i < vva.size(); i++) {
      const double E_tmp = vva[i].at(0);
      if (E_tmp >= min && E_tmp <= max) {
        int column_max = 0; // some default
        if (vva.at(0).size() == 3) {
          column_max = 2; // get rid of the sum of TDOS
        }
        if (vva.at(0).size() == 5) {
          column_max = 3;
        }
        for (int j = 1; j < column_max; j++) {
          const double db_tmp = vva[i][j];
          if (abs(db_tmp) > max_value) {
            max_value = db_tmp;
          }
        }
      }
    }
    return max_value;
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(vector<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const char delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const char);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const char, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const string& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const string&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const xvector<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();

    if (ientries.rows > 2) {
      for (int i = ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) { // CO20180216 - added -1
          output << lDelim;
        } else if (i != ientries.urows) {
          output << delim;
        }
      }
    } else {
      for (int i = ientries.lrows; i <= ientries.urows; i++) {
        output << ientries[i];
        if (i == ientries.urows - 1) { // CO20180216 - added -1
          output << mDelim;
        } else if (i != ientries.urows) {
          output << delim;
        }
      }
    }
    return output.str();
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const xvector<utype>&, const stringstream&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const char delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const char);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const char, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const string& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const string& delimiter, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const vector<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();

    if (ientries.size() > 2) {
      for (size_t i = 0; i < ientries.size(); i++) {
        if constexpr (std::is_same_v<utype, bool>) {
          if (ientries[i]) {
            output << "true";
          } else {
            output << "false";
          }
        } else {
          output << ientries[i];
        }
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (size_t i = 0; i < ientries.size(); i++) {
        if constexpr (std::is_same_v<utype, bool>) {
          if (ientries[i]) {
            output << "true";
          } else {
            output << "false";
          }
        } else {
          output << ientries[i];
        }
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const vector<utype>&, const stringstream&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_BOOL)
#undef AST_TEMPLATE

} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(vector<string>& sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const vector<string>& sentries, const char delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<string>& sentries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(sentries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
  string joinWDelimiter(const vector<string>& sentries, const string& delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<string>& sentries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& sentries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(sentries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
  string joinWDelimiter(const vector<string>& sentries, const stringstream& delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const vector<string>& sentries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const vector<string>& sentries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries_processed;
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (size_t i = 0; i < sentries.size(); i++) {
      if (!sentries[i].empty()) {
        sentries_processed.push_back(sentries[i]);
      }
    }
    if (sentries_processed.size() > 2) {
      for (size_t i = 0; i < sentries_processed.size(); i++) {
        output << sentries_processed[i];
        if (i == sentries_processed.size() - 2) {
          output << lDelim;
        } else if (i != sentries_processed.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (size_t i = 0; i < sentries_processed.size(); i++) {
        output << sentries_processed[i];
        if (i == sentries_processed.size() - 2) {
          output << mDelim;
        } else if (i != sentries_processed.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(deque<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const char delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const char);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const char, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const string& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const string&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const deque<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();
    if (ientries.size() > 2) {
      for (size_t i = 0; i < ientries.size(); i++) {
        output << ientries[i];
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (size_t i = 0; i < ientries.size(); i++) {
        output << ientries[i];
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const deque<utype>&, const stringstream&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(std::set<uint>& uientries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters int/uint type of objects together by a delimiter
  // no point for double objects, faster to just do it on the spot with
  // setprecision,fixed, etc.
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const char delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const char);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const char, const char, const char);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const string& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(ientries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const string&, const string&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const stringstream& delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(ientries, delimiter, delimiter, l_delimiter);
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  template <class utype> string joinWDelimiter(const std::set<utype>& ientries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();
    size_t i;
    if (ientries.size() > 2) {
      i = 0;
      for (auto& entry : ientries) {
        output << entry;
        if (i == ientries.size() - 2) {
          output << lDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
        i++;
      }
    } else {
      i = 0;
      for (auto& entry : ientries) {
        output << entry;
        if (i == ientries.size() - 2) {
          output << mDelim;
        } else if (i != ientries.size() - 1) {
          output << delim;
        }
        i++;
      }
    }
    return output.str();
  }
#define AST_TEMPLATE(utype) template string joinWDelimiter(const std::set<utype>&, const stringstream&, const stringstream&, const stringstream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  //***************************************************************************//
  // aurostd::joinWDelimiter(set<string>& sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const std::set<string>& sentries, const char delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const std::set<string>& sentries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const std::set<string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(sentries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
  string joinWDelimiter(const std::set<string>& sentries, const string& delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const std::set<string>& sentries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const std::set<string>& sentries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(sentries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
  string joinWDelimiter(const std::set<string>& sentries, const stringstream& delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const std::set<string>& sentries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const std::set<string>& sentries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries_processed;
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (const string& entry : sentries) {
      if (!entry.empty()) {
        sentries_processed.push_back(entry);
      }
    }
    if (sentries_processed.size() > 2) {
      for (size_t i = 0; i < sentries_processed.size(); i++) {
        output << sentries_processed[i];
        if (i == sentries_processed.size() - 2) {
          output << lDelim;
        } else if (i != sentries_processed.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (size_t i = 0; i < sentries_processed.size(); i++) {
        output << sentries_processed[i];
        if (i == sentries_processed.size() - 2) {
          output << mDelim;
        } else if (i != sentries_processed.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::joinWDelimiter(deque<string>& sentries,const stringstream&
  // delimiter,const stringstream& m_delimiter,const stringstream& l_delimiter)
  //***************************************************************************//
  // joinWDelimiters string type of objects together by a delimiter
  // m_delimiter is used if input is exactly length 2
  // l_delimiter otherwise
  string joinWDelimiter(const deque<string>& sentries, const char delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<string>& sentries, const char delimiter, const char l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& sentries, const char delimiter, const char m_delimiter, const char l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(sentries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
  string joinWDelimiter(const deque<string>& sentries, const string& delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<string>& sentries, const string& delimiter, const string& l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& sentries, const string& delimiter, const string& m_delimiter, const string& l_delimiter) {
    stringstream delimiter_processed;
    stringstream m_delimiter_processed;
    stringstream l_delimiter_processed;
    delimiter_processed << delimiter;
    m_delimiter_processed << m_delimiter;
    l_delimiter_processed << l_delimiter;
    return joinWDelimiter(sentries, delimiter_processed, m_delimiter_processed, l_delimiter_processed);
  }
  string joinWDelimiter(const deque<string>& sentries, const stringstream& delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, delimiter);
  }
  string joinWDelimiter(const deque<string>& sentries, const stringstream& delimiter, const stringstream& l_delimiter) {
    return joinWDelimiter(sentries, delimiter, delimiter, l_delimiter);
  }
  string joinWDelimiter(const deque<string>& sentries, const stringstream& delimiter, const stringstream& m_delimiter, const stringstream& l_delimiter) {
    stringstream output;
    vector<string> sentries_processed; // no point working with deque
    const string delim = delimiter.str();
    const string mDelim = m_delimiter.str();
    const string lDelim = l_delimiter.str();
    // go through once to eliminate empty strings
    for (size_t i = 0; i < sentries.size(); i++) { // CO20200106 - patching for auto-indenting
      if (!sentries[i].empty()) {
        sentries_processed.push_back(sentries[i]);
      }
    }
    if (sentries_processed.size() > 2) {
      for (size_t i = 0; i < sentries_processed.size(); i++) {
        output << sentries_processed[i];
        if (i == sentries_processed.size() - 2) {
          output << lDelim;
        } else if (i != sentries_processed.size() - 1) {
          output << delim;
        }
      }
    } else {
      for (size_t i = 0; i < sentries_processed.size(); i++) {
        output << sentries_processed[i];
        if (i == sentries_processed.size() - 2) {
          output << mDelim;
        } else if (i != sentries_processed.size() - 1) {
          output << delim;
        }
      }
    }
    return output.str();
  }
} // namespace aurostd

namespace aurostd {
  string wrapString(const string& input, const string& wrapper) {
    return wrapString(input, wrapper, wrapper);
  }
  string wrapString(const string& input, const string& wrapper_start, const string& wrapper_end) {
    if (input.empty()) {
      return input;
    }
    return wrapper_start + input + wrapper_end;
  }
} // namespace aurostd

// DX20180118 START: XCOMPLEX TO JSON
namespace aurostd {
  //***************************************************************************//
  // aurostd::xcomplex2json
  //***************************************************************************//
  template <typename utype> string _xcomplex2json(xcomplex<utype>& number) {
    const string eendl;
    const bool roff = true; // round off
    stringstream sss;
    stringstream sscontent_json;
    vector<string> vcontent_json;
    // real
    sscontent_json << R"("real":")" << aurostd::utype2string(number.re, 5, roff) << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str());
    aurostd::StringstreamClean(sscontent_json);
    // imaginary
    sscontent_json << R"("imag":")" << aurostd::utype2string(number.im, 5, roff) << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str());
    aurostd::StringstreamClean(sscontent_json);

    sss << "{" << aurostd::joinWDelimiter(vcontent_json, ",") << "}" << eendl;
    return sss.str();
  }
} // namespace aurostd

// Need to initalize
namespace aurostd {
  string xcomplex2json(xcomplex<double>& number) {
    return _xcomplex2json(number);
  }
} // namespace aurostd

// DX20180118 END: XCOMPLEX TO JSON

// DX20170803 START: Matrix to JSON
namespace aurostd {
  //***************************************************************************//
  // aurostd::xmatDouble2String(xmatrix<double>& xmat_in)
  //***************************************************************************//
  // converts xmatrix<double> to json string
  string xmatDouble2String(const xmatrix<double>& xmat_in, int precision, bool roff, double tol, char FORMAT) {
    stringstream output;
    vector<string> rows;
    for (int i = xmat_in.lrows; i <= xmat_in.urows; i++) { // DX20180323 - fixed typo for initial index "int i=1" not "int i=xmat_in.urows" //ME20220324 - changed to lrows
      stringstream row;
      const xvector<double> xvec = xmat_in(i); // DX20170822 - added roundoff
      // if(roff){ xvec = roundoff(xvec,tol);} //DX20170822 - added roundoff
      row << "[" << joinWDelimiter(xvecDouble2vecString(xvec, precision, roff, tol, FORMAT), ",") << "]";
      rows.push_back(row.str());
      // cerr << i << "row.str(): " << row.str() << endl;
    }
    output << joinWDelimiter(rows, ",");
    return output.str();
  }

  // ME20220324
  template <typename utype> string xmat2String(const xmatrix<utype>& xmat_in) {
    vector<string> rows;
    for (int i = xmat_in.lrows; i <= xmat_in.urows; i++) {
      rows.push_back("[" + joinWDelimiter(xmat_in(i), ",") + "]");
    }
    return joinWDelimiter(rows, ",");
  }
#define AST_TEMPLATE(utype) template string xmat2String(const xmatrix<utype>&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
} // namespace aurostd
// DX20170803 START: Matrix to END

namespace aurostd {
  //***************************************************************************//
  // aurostd::vecDouble2vecString(vector<double>& vin,int precision)
  //***************************************************************************//
  // converts vector<double> to vector<string> with precision
  // also works for xvectors and deques
  vector<string> vecDouble2vecString(const vector<double>& vin, int precision, bool roff, double tol, char FORMAT) {
    vector<string> vout;
    for (size_t i = 0; i < vin.size(); i++) {
      // double tmp = vin[i]; //DX20170822 - add roundoff
      // if(roff){ tmp=roundoff(tmp,tol); } //DX20170822 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i], precision, roff, tol, FORMAT)); // DX20170822 - add roundoff
    }
    return vout;
  }
  string vecDouble2String(const vector<double>& vin, int precision, bool roff, double tol, char FORMAT) {
    return aurostd::joinWDelimiter(vecDouble2vecString(vin, precision, roff, tol, FORMAT), ",");
  }
  vector<string> xvecDouble2vecString(const xvector<double>& vin, int precision, bool roff, double tol, char FORMAT) {
    vector<string> vout;
    if (vin.rows == 0) {
      return vout;
    }
    for (int i = vin.lrows; i <= vin.urows; i++) {
      // double tmp = vin(i); //DX20170822 - add roundoff
      // if(roff){ tmp=roundoff(tmp,tol); } //DX20170822 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i], precision, roff, tol, FORMAT)); // DX20170822 - add roundoff
    }
    return vout;
  }
  string xvecDouble2String(const xvector<double>& vin, int precision, bool roff, double tol, char FORMAT) {
    return aurostd::joinWDelimiter(xvecDouble2vecString(vin, precision, roff, tol, FORMAT), ",");
  }
  deque<string> vecDouble2vecString(const deque<double>& vin, int precision, bool roff, double tol, char FORMAT) { // SC20200330
    deque<string> vout;
    for (size_t i = 0; i < vin.size(); i++) {
      // double tmp = vin[i]; //DX20170822 - add roundoff
      // if(roff){ tmp=roundoff(tmp,tol); } //DX20170822 - add roundoff
      vout.push_back(aurostd::utype2string(vin[i], precision, roff, tol, FORMAT)); // DX20170822 - add roundoff
    }
    return vout;
  }
  string vecDouble2String(const deque<double>& vin, int precision, bool roff, double tol, char FORMAT) {
    return aurostd::joinWDelimiter(vecDouble2vecString(vin, precision, roff, tol, FORMAT), ",");
  }
} // namespace aurostd

namespace aurostd {
  //***************************************************************************//
  // aurostd::wrapVecEntries(vector<string>& vin,string wrap)
  //***************************************************************************//
  // individually wraps entries of vector with specified string
  // converts <a,b,c> to <'a','b','c'>
  // also works for deques
  vector<string> wrapVecEntries(const vector<string>& vin, const string& wrap) {
    return wrapVecEntries(vin, wrap, wrap);
  }
  vector<string> wrapVecEntries(const vector<string>& vin, const string& wrap_start, const string& wrap_end) {
    vector<string> vout;
    for (size_t i = 0; i < vin.size(); i++) {
      if (!vin[i].empty()) {
        vout.push_back(wrap_start + vin[i] + wrap_end);
      }
    }
    return vout;
  }
  deque<string> wrapVecEntries(const deque<string>& vin, const string& wrap) {
    return wrapVecEntries(vin, wrap, wrap);
  }
  deque<string> wrapVecEntries(const deque<string>& vin, const string& wrap_start, const string& wrap_end) {
    deque<string> vout;
    for (size_t i = 0; i < vin.size(); i++) {
      if (!vin[i].empty()) {
        vout.push_back(wrap_start + vin[i] + wrap_end);
      }
    }
    return vout;
  }
} // namespace aurostd

// base64 stuff
// CO START
namespace aurostd {
  // ***************************************************************************
  // aurostd::isBase64(unsigned char c)
  // ***************************************************************************
  // determines if char is base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  // static inline bool isBase64(unsigned char c)
  inline bool isBase64(unsigned char c) { // CO20200106 - patching for auto-indenting
    return (isalnum(c) || (c == '+') || (c == '/'));
  }

  // ***************************************************************************
  // aurostd::base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len)
  // ***************************************************************************
  // encodes bytes to base64
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  std::string base64Encoder(const unsigned char* bytes_to_encode, unsigned int in_len) {
    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
      char_array_3[i++] = *(bytes_to_encode++);
      if (i == 3) {
        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (i = 0; (i < 4); i++) {
          ret += base64_chars[char_array_4[i]];
        }
        i = 0;
      }
    }

    if (i) {
      for (j = i; j < 3; j++) {
        char_array_3[j] = '\0';
      }

      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;

      for (j = 0; (j < i + 1); j++) {
        ret += base64_chars[char_array_4[j]];
      }

      while ((i++ < 3)) {
        ret += '=';
      }
    }

    return ret;
  }

  // ***************************************************************************
  // aurostd::base64Decoder(std::string const& encoded_string)
  // ***************************************************************************
  // decodes base64 to bytes
  // http://www.adp-gmbh.ch/cpp/common/base64.html
  std::string base64Decoder(const std::string& encoded_string) {
    int in_len = encoded_string.size();
    int i = 0;
    int j = 0;
    int in_ = 0;
    unsigned char char_array_4[4];
    unsigned char char_array_3[3];
    std::string ret;

    while (in_len-- && (encoded_string[in_] != '=') && isBase64(encoded_string[in_])) {
      char_array_4[i++] = encoded_string[in_];
      in_++;
      if (i == 4) {
        for (i = 0; i < 4; i++) {
          char_array_4[i] = base64_chars.find(char_array_4[i]);
        }

        char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (i = 0; (i < 3); i++) {
          ret += char_array_3[i];
        }
        i = 0;
      }
    }

    if (i) {
      for (j = i; j < 4; j++) {
        char_array_4[j] = 0;
      }

      for (j = 0; j < 4; j++) {
        char_array_4[j] = base64_chars.find(char_array_4[j]);
      }

      char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

      for (j = 0; (j < i - 1); j++) {
        ret += char_array_3[j];
      }
    }

    return ret;
  }

  // ***************************************************************************
  // aurostd::bin2base64(const std::string& b_file, std::string& b64String)
  // ***************************************************************************
  // converts binary file to base64 string
  bool bin2base64(const std::string& b_file, std::string& b64String) {
    stringstream output;
    if (!aurostd::FileExist(b_file)) {
      cerr << "ERROR - aurostd::bin2base64: Binary file " << b_file << " does not exist!";
      return false;
    }
    ifstream file(b_file.c_str(), std::ios::in | std::ios::binary);
    output << b64_encoder << file;
    b64String = output.str();
    return true;
  }

  // ***************************************************************************
  // aurostd::base642bin(const std::string& b64String, const std::string& b_file)
  // ***************************************************************************
  // converts base64 string to binary file
  bool base642bin(const std::string& b64String, const std::string& b_file) {
    ofstream output;
    output.open(b_file.c_str(), std::ios::out | std::ios::binary);
    output << b64_decoder << b64String;
    output.flush();
    output.clear();
    output.close();
    return true;
  }

  b64_encoder_proxy operator<<(std::ostream& os, b64_encoder_creator) {
    return b64_encoder_proxy(os);
  }

  b64_decoder_proxy operator<<(std::ostream& os, b64_decoder_creator) {
    return b64_decoder_proxy(os);
  }

} // namespace aurostd
// CO END

// ***************************************************************************
// aurostd::CountofWords
// ***************************************************************************
namespace aurostd {
  int CountWordsinString(string& input) {
    // return the number of words in a string;
    int number = 0;
    string::iterator str_it;
    for (str_it = input.begin(); str_it < input.end(); str_it++) {
      if (!(*str_it == ' ' || *str_it == '\t' || *str_it == '\n')) {
        number++;
        while (str_it != input.end() && !(*str_it == ' ' || *str_it == '\t' || *str_it == '\n')) {
          str_it++;
        }
      }
    }
    return number;
  }
} // namespace aurostd

// ***************************************************************************
// aurostd::CountofWords
// ***************************************************************************
namespace aurostd {
  int CountWordsinString_web(string input) {
    // return the number of words in a string;
    vector<string> vstr;
    aurostd::string2tokens(input, vstr, " ");
    return vstr.size();
  }
} // namespace aurostd

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
