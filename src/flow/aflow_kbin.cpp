// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#include <cstddef>
#include <deque>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <regex>
#include <set>
#include <sstream>
#include <vector>

#include <pthread.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_argv.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_time.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xrandom.h"
#include "AUROSTD/aurostd_xscalar.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "flow/aflow_kvasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "interfaces/aflow_pthreads.h"

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
using std::stringstream;
using std::vector;

#define __XLIBS_LINK
#define cdebug cerr

// bool nocase_compare(char c1,char c2) {return toupper(c1)==toupper(c2);}

// #define MaxAflowInSize 65535
// string AflowIn; //[MaxAflowInSize];

#define VRUNS_MAX_CUTOFF 32768
#define DEFAULT_KILL_MEM_CUTOFF 1.50

namespace aurostd {
  // ***************************************************************************
  // Function DirectoryAlreadyInDabatase
  // ***************************************************************************
  bool DirectoryAlreadyInDatabase(string directory, bool FORCE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (FORCE) {
      return false; // neglect already in the database
    }

    // already scanned
    if (aurostd::FileExist(directory + "/ALREADY_IN_DATABASE") || aurostd::CompressFileExist(directory + "/ALREADY_IN_DATABASE")) {
      return true;
    }

    // no, then scan
    if (LDEBUG) {
      cerr << "SEARCHING" << endl;
    }
    bool already_in_database = false;
    string tmp_directory = directory;
    aurostd::StringSubstInPlace(tmp_directory, _AFLOWIN_, " ");
    aurostd::StringSubstInPlace(tmp_directory, "./", " ");
    aurostd::StringSubstInPlace(tmp_directory, "/", " ");
    vector<string> tokens;
    aurostd::string2tokens(tmp_directory, tokens);
    if (tokens.size() >= 2) {
      tmp_directory = tokens[tokens.size() - 2] + "/" + tokens[tokens.size() - 1];
    }
    if (tokens.size() == 1) {
      tmp_directory = tokens[tokens.size() - 1];
    }

    uint library = LIBRARY_NOTHING;
    // XHOST_LIBRARY_LIB0
    if (aurostd::substring2bool(directory, "LIB0")) {
      library = XHOST_LIBRARY_LIB0;
    }
    // XHOST_LIBRARY_LIB1
    if (aurostd::substring2bool(directory, "LIB1")) {
      library = XHOST_LIBRARY_LIB1;
    }
    // XHOST_LIBRARY_AURO
    if (aurostd::substring2bool(directory, "AURO")) {
      library = XHOST_LIBRARY_LIB1;    // [HISTORIC]
    }
    // XHOST_LIBRARY_LIB2
    if (aurostd::substring2bool(directory, "LIBRARYU")) {
      library = XHOST_LIBRARY_LIB2;
    }
    if (aurostd::substring2bool(directory, "LIBRARYX")) {
      library = XHOST_LIBRARY_LIB2;    // [HISTORIC]
    }
    if (aurostd::substring2bool(directory, "LIB2")) {
      library = XHOST_LIBRARY_LIB2;
    }
    // XHOST_LIBRARY_ICDS
    if (aurostd::substring2bool(directory, "ICSD")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "SCINT")) {
      library = XHOST_LIBRARY_ICSD;    // [HISTORIC]
    }
    if (aurostd::substring2bool(directory, "BCC")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "BCT")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "CUB")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "FCC")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "HEX")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "MCL")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "MCLC")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "ORC")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "ORCC")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "ORCF")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "ORCI")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "RHL")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "TET")) {
      library = XHOST_LIBRARY_ICSD;
    }
    if (aurostd::substring2bool(directory, "TRI")) {
      library = XHOST_LIBRARY_ICSD;
    }
    // XHOST_LIBRARY_ICSD
    if (aurostd::substring2bool(directory, "MAGNETIC")) {
      library = XHOST_LIBRARY_LIB3;    // [HISTORIC]
    }
    if (aurostd::substring2bool(directory, "LIB3")) {
      library = XHOST_LIBRARY_LIB3;
    }
    if (aurostd::substring2bool(directory, "T0001")) {
      library = XHOST_LIBRARY_LIB3;
    }
    if (aurostd::substring2bool(directory, "T0002")) {
      library = XHOST_LIBRARY_LIB3;
    }
    if (aurostd::substring2bool(directory, "T0001")) {
      library = XHOST_LIBRARY_LIB3;
    }
    if (aurostd::substring2bool(directory, "T0004")) {
      library = XHOST_LIBRARY_LIB3;
    }
    if (aurostd::substring2bool(directory, "T0005")) {
      library = XHOST_LIBRARY_LIB3;
    }
    if (aurostd::substring2bool(directory, "T0006")) {
      library = XHOST_LIBRARY_LIB3;
    }
    // XHOST_LIBRARY_LIB4
    if (aurostd::substring2bool(directory, "LIB4")) {
      library = XHOST_LIBRARY_LIB4;
    }
    if (aurostd::substring2bool(directory, "Q0001")) {
      library = XHOST_LIBRARY_LIB4;
    }
    // XHOST_LIBRARY_LIB5
    if (aurostd::substring2bool(directory, "LIB5")) {
      library = XHOST_LIBRARY_LIB5;
    }
    if (aurostd::substring2bool(directory, "P0001")) {
      library = XHOST_LIBRARY_LIB5;
    }
    // XHOST_LIBRARY_LIB6
    if (aurostd::substring2bool(directory, "LIB6")) {
      library = XHOST_LIBRARY_LIB6;
    }
    if (aurostd::substring2bool(directory, "H0001")) {
      library = XHOST_LIBRARY_LIB6;
    }
    // XHOST_LIBRARY_LIB7
    if (aurostd::substring2bool(directory, "LIB7")) {
      library = XHOST_LIBRARY_LIB7;
    }
    // XHOST_LIBRARY_LIB8
    if (aurostd::substring2bool(directory, "LIB8")) {
      library = XHOST_LIBRARY_LIB8;
    }
    // XHOST_LIBRARY_LIB9
    if (aurostd::substring2bool(directory, "LIB9")) {
      library = XHOST_LIBRARY_LIB9;
    }

    // found something
    if (library != LIBRARY_NOTHING) {
      // TODO replace the check with lookup in local SQLITE DB
      tokens.clear();
      string tmp;
      vector<string> vLibrary;
      if (library == XHOST_LIBRARY_LIB0) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB0" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB0_RAW;
      }
      if (library == XHOST_LIBRARY_LIB1) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB1" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB1_RAW;
      }
      if (library == XHOST_LIBRARY_LIB2) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB2" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB2_RAW;
      }
      if (library == XHOST_LIBRARY_LIB3) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB3" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB3_RAW;
      }
      if (library == XHOST_LIBRARY_LIB4) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB4" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB4_RAW;
      }
      if (library == XHOST_LIBRARY_LIB5) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB5" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB5_RAW;
      }
      if (library == XHOST_LIBRARY_LIB6) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB6" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB6_RAW;
      }
      if (library == XHOST_LIBRARY_LIB7) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB7" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB7_RAW;
      }
      if (library == XHOST_LIBRARY_LIB8) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB8" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB8_RAW;
      }
      if (library == XHOST_LIBRARY_LIB9) {
        if (LDEBUG) {
          cerr << "library==XHOST_LIBRARY_LIB9" << endl;
        }
        tokens = XHOST_Library_CALCULATED_LIB9_RAW;
      }
      for (size_t i = 0; i < tokens.size(); i++) {
        if (aurostd::substring2bool(tokens[i], "/")) {
          tmp = tokens[i];
          aurostd::StringSubstInPlace(tmp, " ", "");
          vLibrary.push_back(tmp);
        }
      }
      for (int i = vLibrary.size() - 1; i >= 0; i--) {
        if (tmp_directory == vLibrary[i]) {
          already_in_database = true;
          if (LDEBUG) {
            cerr << vLibrary[i] << " FOUND .." << endl; // NEW
          }
        }
      }
    }

    if (LDEBUG) {
      cerr << "DONE..." << endl;
    }
    if (LDEBUG) {
      cerr << tmp_directory << endl;
    }

    if (already_in_database) {
      //    cerr << directory << " already in database" << endl;
      aurostd::CopyFile(directory + "/" + _AFLOWIN_, directory + "/ALREADY_IN_DATABASE");
      aurostd::execute(DEFAULT_KZIP_BIN + " -f " + directory + "/" + _AFLOWIN_);
    }

    return already_in_database;
  }
} // namespace aurostd

using aurostd::DirectoryAlreadyInDatabase;

// int KBIN_MODE;

// #define KBIN_VOID_MODE 0             // just a shift
// #define KBIN_VASP_MODE 2             // for ab-initio VASP mode
// #define KBIN_XXXX_MODE 3             // for XXXX program
// #define KBIN_GRND_MODE 4             // for classical monte carlo

// GND MODE
#define KBIN_VASP_N_VPARS 32
#define _KBIN_VASP_SLEEP_ 2
#define _KBIN_LOOP_SLEEP_ 300
// PRIORITY
#define PRIORITY_PROBABILITY 0.2000
// #define PRIORITY_GREP_STRING string("grep -vi xxxxx ")
#define PRIORITY_GREP_STRING string("grep -vi PRIORITY ")

// ***************************************************************************
// KBIN::Legitimate_krun
// ***************************************************************************
namespace KBIN {
  bool Legitimate_krun(const _aflags& aflags, const bool osswrite, ostringstream& oss) {
    string aflowindir = aflags.Directory;
    string filename;
    aurostd::StringSubstInPlace(aflowindir, "//", "/");
    if (!Legitimate_aflowdir(aflowindir, aflags, osswrite, oss)) { // aflowdir not legitimate
      return false;
    }
    if (aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
      filename = aflowindir + "/INCAR";
    } else if (aflags.KBIN_RUN_AFLOWIN || aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AIMS_FROM_AFLOWIN) {
      filename = aflowindir + "/" + _AFLOWIN_;
    } else { // SD20220319 - catch-all
      return false;
    }
    if (!aurostd::FileExist(filename)) { // file does not exist
      if (osswrite) {
        oss << "MMMMM  File does not exist = " << filename << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    return true;
  }
} // namespace KBIN

namespace KBIN {
  bool Legitimate_krun(const _aflags& aflags) {
    ostringstream aus;
    return KBIN::Legitimate_krun(aflags, false, aus);
  };
} // namespace KBIN

// ***************************************************************************
// KBIN::Legitimate_aflowin
// ***************************************************************************
namespace KBIN {
  bool Legitimate_aflowin(const string& _aflowindir, const bool osswrite, ostringstream& oss) {
    string aflowindir = _aflowindir; // SD20220321 - aflowindir = aflowdir/_AFLOWIN_
    aurostd::StringSubstInPlace(aflowindir, "//", "/");
    if (!aurostd::FileExist(aflowindir)) { // file does not exist
      if (osswrite) {
        oss << "MMMMM  File does not exist = " << aflowindir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    if (aurostd::FileEmpty(aflowindir)) { // file is empty
      if (osswrite) {
        oss << "MMMMM  File is empty = " << aflowindir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    if (!aurostd::substring2bool(aflowindir, _AFLOWIN_)) { // file name does not contain _AFLOWIN_
      if (osswrite) {
        oss << "MMMMM  File name does not contain _AFLOWIN_ " << _AFLOWIN_ << " = " << aflowindir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    aurostd::StringSubstInPlace(aflowindir, _AFLOWIN_, "");
    if (aurostd::DirectoryLocked(aflowindir, _AFLOWLOCK_)) { // LOCK file exists
      if (osswrite) {
        oss << "MMMMM  Directory locked = " << aflowindir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    return true;
  }
} // namespace KBIN

namespace KBIN {
  bool Legitimate_aflowin(const string aflowindir) {
    ostringstream aus;
    return KBIN::Legitimate_aflowin(aflowindir, false, aus);
  };
} // namespace KBIN

// ***************************************************************************
// KBIN::Legitimate_aflowdir
// ***************************************************************************
namespace KBIN {
  bool Legitimate_aflowdir(const string& aflowdir, const _aflags& aflags, const bool osswrite, ostringstream& oss) {
    if (!aurostd::IsDirectory(aflowdir)) { // directory does not exist
      if (osswrite) {
        oss << "MMMMM  Directory does not exist = " << aflowdir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    if (!aurostd::DirectoryWritable(aflowdir)) { // directory unwritable
      if (osswrite) {
        oss << "MMMMM  Directory unwritable = " << aflowdir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    if (aurostd::DirectoryLocked(aflowdir, _AFLOWLOCK_)) { // directory locked
      if (osswrite) {
        oss << "MMMMM  Directory locked = " << aflowdir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    if (aurostd::DirectorySkipped(aflowdir)) { // directory skipped
      if (osswrite) {
        oss << "MMMMM  Directory skipped = " << aflowdir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    if (DirectoryAlreadyInDatabase(aflowdir, aflags.AFLOW_FORCE_RUN)) { // directory already in database
      if (osswrite) {
        oss << "MMMMM  Directory already in database = " << aflowdir << Message(__AFLOW_FILE__) << endl;
        aurostd::PrintMessageStream(oss, XHOST.QUIET);
      };
      return false;
    }
    return true;
  }
} // namespace KBIN

namespace KBIN {
  bool Legitimate_aflowdir(const string& aflowdir, const _aflags& aflags) {
    ostringstream aus;
    return KBIN::Legitimate_aflowdir(aflowdir, aflags, false, aus);
  };
} // namespace KBIN

namespace KBIN {
  void getAflowInFromAFlags(const _aflags& aflags, string& AflowIn_file, string& AflowIn, ostream& oss) {
    ofstream FileMESSAGE;
    return getAflowInFromAFlags(aflags, AflowIn_file, AflowIn, FileMESSAGE, oss);
  }  // CO20191110
  void getAflowInFromAFlags(const _aflags& aflags, string& AflowIn_file, string& AflowIn, ofstream& FileMESSAGE, ostream& oss) { // CO20191110
    return getAflowInFromDirectory(aflags.Directory, AflowIn_file, AflowIn, FileMESSAGE, oss);
  }
  void getAflowInFromDirectory(const string& directory, string& AflowIn_file, string& AflowIn, ostream& oss) {
    ofstream FileMESSAGE;
    return getAflowInFromDirectory(directory, AflowIn_file, AflowIn, FileMESSAGE, oss);
  }  // CO20191110
  void getAflowInFromDirectory(const string& directory, string& AflowIn_file, string& AflowIn, ofstream& FileMESSAGE, ostream& oss) { // CO20191110
    AflowIn_file = aurostd::CleanFileName(directory + "/" + _AFLOWIN_); // CO20200624
    if (!aurostd::FileExist(AflowIn_file)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Input file does not exist: " + AflowIn_file, _INPUT_ERROR_);
    }
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Using input file: " + AflowIn_file, FileMESSAGE, oss, _LOGGER_MESSAGE_);
    aurostd::file2string(AflowIn_file, AflowIn);
    AflowIn = aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "AflowIn.size()=" + aurostd::utype2string(AflowIn.size()), FileMESSAGE, oss, _LOGGER_MESSAGE_); // CO20200624 - check size!=0
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::Main
// ***************************************************************************
namespace KBIN {
  int KBIN_Main(vector<string> argv) {        // AFLOW_FUNCTION_IMPLEMENTATION
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (argv.empty()) {
      return 0;
    } // SD2022024 - check argv!=0
    //  string Directory;
    int i;
    ostringstream aus;
    const ifstream FileAUS;
    _aflags aflags;
    aflags.Directory = XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");  // CO20190629 - a good default until we get down into vDirectory
    aurostd::StringstreamClean(aus);
    const bool _VERBOSE_ = false;
    const int XHOST_AFLOW_RUNXnumber_multiplier = 3;

    const std::deque<_aflags> qaflags;

    AFLOW_PTHREADS::Clean_Threads();                                    // clean threads
    // _aflags taflags[MAX_ALLOCATABLE_PTHREADS];
    // _threaded_KBIN_params params[MAX_ALLOCATABLE_PTHREADS];

    // cerr << "GMODE" << endl;
    // check BlackList **************************************************
    if (AFLOW_BlackList(XHOST.hostname)) {
      aus << "MMMMM  HOSTNAME BLACKLISTED = " << XHOST.hostname << " - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      return 0;
    }
    // get KBIN **************************************************
    //  aus << "MMMMM  AFLOW: running " << _AFLOWIN_ << ": " << " "  << " - " << Message(__AFLOW_FILE__,aflags) << endl;   // too much verbosity is annoying
    //  aurostd::PrintMessageStream(aus,XHOST.QUIET);                                                  // too much verbosity is annoying

    // do some running

    //   cerr << __AFLOW_FUNC__ << " XHOST.AFLOW_RUNXflag=" << XHOST.AFLOW_RUNXflag << endl; //

    if (XHOST.AFLOW_RUNXflag) {
      vector<string> tokens;
      if (aurostd::args2attachedflag(argv, "--run=")) {
        aurostd::string2tokens(aurostd::args2attachedstring(argv, "--run=", "1"), tokens, "=");
        XHOST.AFLOW_RUNXnumber = aurostd::string2utype<uint>(tokens[tokens.size() - 1]);
      }
      if (aurostd::args2attachedflag(argv, "-run=")) {
        aurostd::string2tokens(aurostd::args2attachedstring(argv, "-run=", "1"), tokens, "=");
        XHOST.AFLOW_RUNXnumber = aurostd::string2utype<uint>(tokens[tokens.size() - 1]);
      }
      if (XHOST.AFLOW_RUNXnumber < 1) {
        XHOST.AFLOW_RUNXnumber = 1;
      }
    }

    if (aurostd::args2flag(argv, "-runone|--runone|-run_one|--run_one|-run1|--run1|--run=1|-run=1")) { // RUNONE COMPATIBILITY
      XHOST.AFLOW_RUNXflag = true;
      XHOST.AFLOW_RUNXnumber = 1;
    }

    if (XHOST.AFLOW_RUNXflag) {
      XHOST.AFLOW_RUNDIRflag = false;
      XHOST.AFLOW_MULTIflag = false;
    }
    if (XHOST.AFLOW_MULTIflag) {
      XHOST.AFLOW_RUNDIRflag = false;
      XHOST.AFLOW_RUNXflag = false;
    }

    //    cerr << __AFLOW_FUNC__ << " XHOST.AFLOW_RUNXflag=" << XHOST.AFLOW_RUNXflag << endl; //

    if (XHOST.vflag_aflow.flag("LOOP")) {
      aus << "MMMMM  KBIN option XHOST.vflag_aflow.flag(\"LOOP\") - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (XHOST.AFLOW_RUNDIRflag) {
      aus << "MMMMM  KBIN option [--run] (XHOST.AFLOW_RUNDIRflag=true) - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (XHOST.AFLOW_MULTIflag) {
      aus << "MMMMM  KBIN option [--run=multi] (XHOST.AFLOW_MULTIflag=true) - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (XHOST.AFLOW_RUNXflag && XHOST.AFLOW_RUNXnumber == 1) {
      aus << "MMMMM  KBIN option [--run=1] (XHOST.AFLOW_RUNXflag=true, XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (XHOST.AFLOW_RUNXflag && XHOST.AFLOW_RUNXnumber > 1) {
      aus << "MMMMM  KBIN option [--run=N] (XHOST.AFLOW_RUNXflag=true, XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }

    aflags.KBIN_RUN_AFLOWIN = true;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN = false;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP = false;
    // DX
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN = false;
    // DX
    aflags.KBIN_DELETE_AFLOWIN = false;

    //--generate could be for VASP/AIMS/etc., need to read aflow.in later
    aflags.KBIN_GEN_GENERAL = aurostd::args2flag(argv, "--generate");                             // CO20180402 - we will use this to modify other GEN flags later, when reading AFLOWIN
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN = aurostd::args2flag(argv, "--generate_vasp_from_aflowin|--generate");// CO20180402 - --generate assumes vasp generation FOR NOW
    if (aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
      aflags.KBIN_RUN_AFLOWIN = false;
    }
    aflags.KBIN_GEN_AIMS_FROM_AFLOWIN = aurostd::args2flag(argv, "--generate_aims_from_aflowin");// CO20180402
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP = aurostd::args2flag(argv, "--generate_aflowin_from_vasp");
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN = aurostd::args2flag(argv, "--generate_symmetry|--generate_sym"); // DX

    // DX
    // DX if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AFLOWIN_FROM_VASP)
    if (aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AFLOWIN_FROM_VASP || aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN || aflags.KBIN_GEN_AIMS_FROM_AFLOWIN) // CO20180409
    { // CO20200106 - patching for auto-indenting
      // DX
      XHOST.AFLOW_RUNXflag = true;
      XHOST.AFLOW_RUNXnumber = 1;
      XHOST.AFLOW_RUNDIRflag = true;
      if (aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
        aflags.KBIN_RUN_AFLOWIN = false;
      }
    }

    aflags.KBIN_DELETE_AFLOWIN = aurostd::args2flag(argv, "--delete_aflowin");

    aflags.AFLOW_MODE_QSUB_MODE1 = aurostd::args2flag(argv, "--qsub1|-qsub1");
    aflags.AFLOW_MODE_QSUB_MODE2 = aurostd::args2flag(argv, "--qsub2|-qsub2");
    aflags.AFLOW_MODE_QSUB_MODE3 = aurostd::args2flag(argv, "--qsub3|-qsub3");

    aflags.AFLOW_FORCE_MPI = XHOST.MPI; //[CO20210315]aurostd::args2flag(argv,"--MPI|--mpi");
    aflags.AFLOW_FORCE_SERIAL = aurostd::args2flag(argv, "--nompi|-nompi|--serial|-serial");
    aflags.AFLOW_GLOBAL_NCPUS = aurostd::args2attachedutype<int>(argv, "--np=", 0);
    // if(aflags.AFLOW_GLOBAL_NCPUS && !MPI && !aflags.AFLOW_FORCE_MPI) {}
    if (XHOST.MPI || aflags.AFLOW_FORCE_MPI) {
      AFLOW_PTHREADS::No_Threads();
    }

    aflags.AFLOW_MACHINE_GLOBAL.clear();
    // "MACHINE::DUKE_BETA_MPICH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH", aurostd::args2flag(argv, "--machine=beta|--machine=duke_beta|--beta|--duke_beta|--machine=beta_mpich|--machine=duke_beta_mpich"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::DUKE_BETA_OPENMPI"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI", aurostd::args2flag(argv, "--machine=beta_openmpi|--machine=duke_beta_openmpi"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::DUKE_QRATS_MPICH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH", aurostd::args2flag(argv, "--machine=qrats|--machine=duke_qrats|--machine=qrats_mpich|--machine=duke_qrats_mpich"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::DUKE_QFLOW_OPENMPI"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QFLOW_OPENMPI", aurostd::args2flag(argv,
                                                                                       "--machine=qflow|--machine=duke_qflow|--machine=qflow_openmpi|--machine=duke_qflow_openmpi|--machine=quser|--machine=duke_"
                                                                                       "quser|--machine=quser_openmpi|--machine=duke_quser_openmpi")); // backwards compatible //CO20180409
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // CO20201220 X START
    //  "MACHINE::DUKE_X_X"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_X", aurostd::args2flag(argv, "--machine=x_x|--machine=duke_x_x")); // backwards compatible //CO20180409
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_X")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::DUKE_X_CRAY"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_CRAY", aurostd::args2flag(argv, "--machine=x_cray|--machine=duke_x_cray")); // SD20221006
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_CRAY")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::DUKE_X_OLDCRAY"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_OLDCRAY", aurostd::args2flag(argv, "--machine=x_oldcray|--machine=duke_x_oldcray")); // SD20221006
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_OLDCRAY")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::DUKE_X_SMB"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_SMB", aurostd::args2flag(argv, "--machine=x_smb|--machine=duke_x_smb")); // SD20221006
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_SMB")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // CO20201220 X STOP
    // CO20220818 JHU_ROCKFISH START
    //  "MACHINE::JHU_ROCKFISH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::JHU_ROCKFISH", aurostd::args2flag(argv, "--machine=rockfish|--machine=jhu_rockfish")); // backwards compatible //CO20180409
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::JHU_ROCKFISH")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // CO20220818 JHU_ROCKFISH STOP
    //  "MACHINE::MPCDF_EOS"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS", aurostd::args2flag(argv, "--machine=eos|--machine=mpcdf_eos|--machine=eos_mpiifort|--machine=mpcdf_eos_mpiifort"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::MPCDF_DRACO"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO", aurostd::args2flag(argv, "--machine=draco|--machine=mpcdf_draco|--machine=draco_mpiifort|--machine=mpcdf_draco_mpiifort"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::MPCDF_COBRA"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_COBRA", aurostd::args2flag(argv, "--machine=cobra|--machine=mpcdf_cobra|--machine=cobra_mpiifort|--machine=mpcdf_cobra_mpiifort"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_COBRA")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // "MACHINE::MPCDF_HYDRA"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA", aurostd::args2flag(argv, "--machine=hydra|--machine=mpcdf_hydra|--machine=hydra_mpiifort|--machine=mpcdf_hydra_mpiifort"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // DX20190509 - MACHINE001 - START
    //  "MACHINE::MACHINE001"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001", aurostd::args2flag(argv, "--machine=machine001"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // DX20190509 - MACHINE001 - END
    // DX20190509 - MACHINE002 - START
    //  "MACHINE::MACHINE002"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002", aurostd::args2flag(argv, "--machine=machine002"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // DX20190509 - MACHINE002 - END
    // DX20201005 - MACHINE003 - START
    //  "MACHINE::MACHINE003"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE003", aurostd::args2flag(argv, "--machine=machine003"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE003")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // DX20201005 - MACHINE003 - END
    // DX20211011 - MACHINE004 - START
    //  "MACHINE::MACHINE004"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE004", aurostd::args2flag(argv, "--machine=machine004"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE004")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // DX20211011 - MACHINE004 - END
    //  DUKE_MATERIALS
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_MATERIALS", aurostd::args2flag(argv, "--machine=materials|--machine=duke_materials"));
    // DUKE_AFLOWLIB
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_AFLOWLIB", aurostd::args2flag(argv, "--machine=aflowlib|--machine=duke_aflowlib"));
    // DUKE_HABANA
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_HABANA", aurostd::args2flag(argv, "--machine=habana|--machine=duke_habana"));
    // FULTON_MARYLOU
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::FULTON_MARYLOU", aurostd::args2flag(argv, "--machine=marylou|--machine=fulton_marylou"));
    // MACHINE2
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::OHAD", aurostd::args2flag(argv, "--machine=ohad|--machine=machine2")); // CO20181113
    // MACHINE1
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::HOST1", aurostd::args2flag(argv, "--machine=host1|--machine=machine1")); // CO20181113
    // DX20190107 - CMU EULER - START
    //  "MACHINE::CMU_EULER"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::CMU_EULER", aurostd::args2flag(argv, "--machine=euler|--machine=cmu_euler"));
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::CMU_EULER")) {
      XHOST.maxmem = DEFAULT_KILL_MEM_CUTOFF;
    }
    // DX20190107 - CMU EULER - END

    // HE20220309 - Machine Name - START
    //  allow the use of generic machine templates while using descriptive names
    stringstream machine_name;
    const string argv_machine_name = aurostd::args2attachedstring(argv, "--machine_name=", "");
    if (!argv_machine_name.empty()) {
      machine_name << argv_machine_name << "(" << aflags.AFLOW_MACHINE_GLOBAL << ")";
    } else {
      machine_name << aflags.AFLOW_MACHINE_GLOBAL;
    }
    aflags.AFLOW_MACHINE_GLOBAL.addattachedscheme("NAME", machine_name.str(), true);
    // HE20220309 - Machine Name - END

    // turn off multi pthreads on specific machines
    if (aflags.AFLOW_MACHINE_GLOBAL.flag()) {
      AFLOW_PTHREADS::No_Threads();
      AFLOW_PTHREADS::FLAG = false; // safety...
      AFLOW_PTHREADS::MAX_PTHREADS = 1; // safety...
      //    kflags.KBIN_MPI=true; // overrides the MPI for machines
      XHOST.MPI = true;
    }

    vector<string> vruns;

    aflags.AFLOW_PERFORM_CLEAN = XHOST.vflag_aflow.flag("CLEAN");// || XHOST.vflag_aflow.flag("XCLEAN"));
    aflags.AFLOW_PERFORM_DIRECTORY = XHOST.vflag_control.flag("VDIR");
    aflags.AFLOW_PERFORM_FILE = XHOST.vflag_control.flag("FILE");
    aflags.AFLOW_PERFORM_ORDER_SORT = aurostd::args2flag(argv, "--sort|-sort");                    // Sorts the _AFLOWIN_ in the list
    aflags.AFLOW_PERFORM_ORDER_REVERSE = aurostd::args2flag(argv, "--reverse|--rsort|-reverse|-rsort"); // Reverse the _AFLOWIN_ in the list
    aflags.AFLOW_PERFORM_ORDER_RANDOM = aurostd::args2flag(argv, "--random|--rnd|-random|-rnd"); // Randomize the _AFLOWIN_ in the list
    aflags.AFLOW_FORCE_RUN = aurostd::args2flag(argv, "--force|-force");

    if (aflags.AFLOW_PERFORM_DIRECTORY) {
      aus << "MMMMM  KBIN option PERFORM_DIRECTORY - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (aflags.AFLOW_PERFORM_FILE) {
      aus << "MMMMM  KBIN option PERFORM_FILE - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (aflags.AFLOW_PERFORM_ORDER_SORT) {
      aus << "MMMMM  KBIN option ORDER_SORT - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (aflags.AFLOW_PERFORM_ORDER_REVERSE) {
      aus << "MMMMM  KBIN option ORDER_REVERSE - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (aflags.AFLOW_PERFORM_ORDER_RANDOM) {
      aus << "MMMMM  KBIN option ORDER_RANDOM - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    if (aflags.AFLOW_FORCE_RUN) {
      aus << "MMMMM  KBIN option FORCE_RUN - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }

    uint maxcheck = VRUNS_MAX_CUTOFF;
    if (XHOST.AFLOW_MULTIflag) {
      maxcheck = VRUNS_MAX_CUTOFF;
    }
    if (XHOST.AFLOW_RUNXflag) {
      maxcheck = XHOST.AFLOW_RUNXnumber;
    }
    if (maxcheck == 0) {
      maxcheck = 1; // safety
    }

    // simple commands
    if (XHOST.vflag_aflow.flag("XCLEAN")) {
      KBIN::XClean(XHOST.vflag_aflow.getattachedscheme("XCLEAN"));
      return 1;
    }

    // for directory mode load them all
    if (aflags.AFLOW_PERFORM_DIRECTORY) {
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VDIR"), vruns, ",");
      if (LDEBUG) {
        for (size_t i = 0; i < vruns.size(); i++) {
          cerr << XPID << "KBIN::Main: vruns[i]=" << vruns[i] << endl;
        }
      }
    } else {
      if (!aflags.AFLOW_PERFORM_FILE) {
        vruns.push_back(aurostd::getPWD()); // CO20191112
      }
      // vruns.push_back(aurostd::execute2string(XHOST.command("pwd"))+" ./");
    }
    // if file found
    if (aflags.AFLOW_PERFORM_FILE) {
      const string file_name = XHOST.vflag_control.getattachedscheme("FILE");
      if (!aurostd::FileExist(file_name)) {
        aus << "EEEEE  FILE_NOT_FOUND = " << file_name << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        return 1;  // CO20200624 - previously exit
      }
      if (aurostd::FileEmpty(file_name)) {
        aus << "EEEEE  FILE_EMPTY = " << file_name << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        return 1;  // CO20200624 - previously exit
      }
      aus << "MMMMM  Loading File = " << file_name << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      vector<string> vlines;
      vlines.clear();
      aurostd::file2vectorstring(file_name, vlines);

      aus << "MMMMM  " << aurostd::PaddedPOST("Legitimate VLINES = " + aurostd::utype2string(vlines.size()), 40) << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      aus << "MMMMM  " << aurostd::PaddedPOST("         maxcheck = " + aurostd::utype2string(maxcheck), 40) << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);

      if (aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
        aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        aurostd::sort(vlines);
      }
      if (aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
        aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        aurostd::rsort(vlines);
      }
      if (aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
        aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] start (XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        aurostd::random_shuffle(vlines);  // uses the std library but the seed is initialized in xrandom too
        aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] stop (XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
      }

      for (size_t i = 0; (i < vlines.size() && vruns.size() < XHOST_AFLOW_RUNXnumber_multiplier * maxcheck && vruns.size() < VRUNS_MAX_CUTOFF); i++) { // XHOST_AFLOW_RUNXnumber_multiplier times more... for safety
        if (Legitimate_aflowin(vlines[i], false, aus)) {
          vruns.push_back(vlines[i]);
        } // true puts too much verbosity
      }

      aus << "MMMMM  " << aurostd::PaddedPOST("Legitimate VRUNS = " + aurostd::utype2string(vruns.size()), 40) << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      // aus << "MMMMM  Legitimate VRUNS = " << vruns.size() << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // cerr << vlines.size() << endl;
    }

    if (aflags.AFLOW_PERFORM_CLEAN && !XHOST.AFLOW_RUNDIRflag && !XHOST.AFLOW_MULTIflag && !XHOST.AFLOW_RUNXflag) {
      XHOST.AFLOW_RUNDIRflag = true; // give something to clean
    }

    if (aflags.AFLOW_PERFORM_CLEAN && !aflags.AFLOW_PERFORM_DIRECTORY) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "to use --clean, you must specify one or more directories", _INPUT_MISSING_); // CO20200624
    }
    if (aflags.AFLOW_PERFORM_CLEAN && aflags.AFLOW_PERFORM_DIRECTORY) {
      //    cout << "DEBUG CLEAN = " << vruns.size() << endl;
    }

    const bool STOP_DEBUG = aurostd::args2flag(argv, "--STOP|--stop");

    // cdebug << "aflags.KBIN_RUN_AFLOWIN=" << aflags.KBIN_RUN_AFLOWIN << endl;
    // cdebug << "aflags.KBIN_GEN_VASP_FROM_AFLOWIN=" << aflags.KBIN_GEN_VASP_FROM_AFLOWIN << endl;
    // cdebug << "aflags.KBIN_GEN_AFLOWIN_FROM_VASP=" << aflags.KBIN_GEN_AFLOWIN_FROM_VASP << endl;

    // ------------------------------------------------------------------------------------------------------------------------------------
    // nothing to run
    if (!XHOST.AFLOW_RUNDIRflag && !XHOST.AFLOW_MULTIflag && !XHOST.AFLOW_RUNXflag && !aflags.AFLOW_PERFORM_CLEAN && !aflags.AFLOW_PERFORM_DIRECTORY) {
      aus << "MMMMM  KBIN option nothing to run [--run , --run=multi, --run=1, --run=N] - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      //    return 0;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // nothing specified : CHECK IF DAEMON
    if (XHOST.AFLOW_RUNDIRflag) { // check if daemaon aflowd
      const string progname = argv[0];
      if (aurostd::substring2bool(progname, "aflowd")) {
        XHOST.AFLOW_MULTIflag = true;
        XHOST.AFLOW_RUNXflag = false;
        XHOST.vflag_aflow.flag("LOOP", true); // add automatically
        aus << "MMMMM  AFLOW: running as DAEMON (aflow --kmode --multi --loop): " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
      }
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // run specified directory: XHOST.AFLOW_RUNDIRflag
    if (XHOST.AFLOW_RUNDIRflag) {
      //    bool krun=true;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " STEP0b" << endl;
      }
      vector<string> vDirectory = vruns;
      // fix the RUNS
      for (size_t ii = 0; ii < vDirectory.size(); ii++) {
        if (aurostd::substring2bool(vDirectory[ii], _AFLOWIN_)) {
          aurostd::StringSubstInPlace(vDirectory[ii], _AFLOWIN_, "");
        }
      }

      // if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
      //   aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::sort(vDirectory);}
      // if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
      //   aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::rsort(vDirectory);}
      // if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
      //   aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] - " << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::random_shuffle(vDirectory);}  // uses the std library but the seed is initialized in xrandom too

      // aurostd::random_shuffle(vDirectory);
      // std::random_shuffle(vDirectory.begin(),vDirectory.end());

      aurostd::xoption opts_clean; // CO20210716
      opts_clean.flag("SAVE_CONTCAR", aurostd::args2flag(argv, "--contcar_save|--save_contcar")); // CO20210716 - saves contcar no matter what
      opts_clean.flag("SAVE_CONTCAR_OUTCAR_COMPLETE", aurostd::args2flag(argv, "--contcar_save_outcar_complete|--save_contcar_outcar_complete")); // CO20210716 - saves contcar only if outcar is complete

      for (size_t idir = 0; idir < vDirectory.size(); idir++) {
        bool krun = true;
        aflags.Directory = vDirectory[idir];
        aus << "MMMMM  AFLOW: running " << _AFLOWIN_ << ", directory" << "=\"" << aflags.Directory << "\" - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        // If necessary PERFORM CLEAN
        if (krun && aflags.AFLOW_PERFORM_CLEAN) {
          aflags.Directory = vDirectory[idir];
          //  cerr << aflags.Directory << endl;
          KBIN::Clean(aflags, opts_clean);
          krun = false;
        }
        // RUN
        // cerr << "STEP0b" << endl;
        if (krun) {
          if (Legitimate_krun(aflags)) {
            if (aflags.KBIN_RUN_AFLOWIN) {
              KBIN::RUN_Directory(aflags);
            }
            if (aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
              KBIN::GenerateAflowinFromVASPDirectory(aflags);
            }
            aus << "MMMMM  AFLOW: Done " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(aus, XHOST.QUIET);
          }
        }
      } // idir
      // aus << "MMMMM  AFLOW: Done " << " - " << Message(__AFLOW_FILE__,aflags) << endl;
      // aurostd::PrintMessageStream(aus,XHOST.QUIET);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " STEP2" << endl;
      }
    }

    // ERRORS ------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------------------------------
    // run MULTI and XHOST.AFLOW_RUNXflag (in XHOST.AFLOW_RUNXflag, runs only XHOST.AFLOW_RUNXnumber and then dies)
    // MULTI with SINGLE AND MULTI THREAD VERSION -----------------------------------------------------------------------------------------
    if (XHOST.AFLOW_MULTIflag || XHOST.AFLOW_RUNXflag) {
      uint RUN_times = 0;
      const bool MULTI_DEBUG = false;
      if (MULTI_DEBUG) {
        aus << "MMMMM  AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
      }
      if (MULTI_DEBUG) {
        aus << "MMMMM  AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
      }
      bool EMPTY = false;
      bool FOUND = false;
      bool free_thread;
      int ithread = 0;
      if (AFLOW_PTHREADS::FLAG && AFLOW_PTHREADS::MAX_PTHREADS <= 1) {
        aus << "EEEEE  ERROR PTHREADS" << endl;
        aus << "MMMMM  AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << endl;
        aus << "MMMMM  AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        return 1; // CO20200624 - previously exit
      }
      if (AFLOW_PTHREADS::MAX_PTHREADS > 1) {
        aus << "MMMMM  AFLOW: MULTI THREAD START: phread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
      }
      aus << "MMMMM  AFLOW: searching subdirectories [d1] " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      while (!EMPTY) {
        vector<string> vaflowin;
        //   cerr << "aflags.AFLOW_PERFORM_DIRECTORY=" << aflags.AFLOW_PERFORM_DIRECTORY << endl;
        // cerr << "aflags.AFLOW_PERFORM_FILE=" << aflags.AFLOW_PERFORM_FILE << endl;

        if (aflags.AFLOW_PERFORM_FILE == false) { // NO FILE SPECIFIED = standard  RUN DIRECTORY
          aus << "MMMMM  AFLOW: aflags.AFLOW_PERFORM_FILE==false" << " - " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);

          const string FileNameSUBDIR = aurostd::TmpFileCreate("RUN");
          aurostd::RemoveFile(FileNameSUBDIR);

          bool isPRIORITY = false;
          // NEW, the sorting is done internally (speed and reliability)
          for (int ifind = 0; ifind < (int) vruns.size(); ifind++) {
            isPRIORITY = (aurostd::substring2bool(vruns[ifind], "PRIORITY") || aurostd::substring2bool(vruns[ifind], "priority"));
            if ((isPRIORITY && aurostd::uniform(1.0) <= PRIORITY_PROBABILITY) || !isPRIORITY) {
              aus << "find " << vruns[ifind] << " " << XHOST.Find_Parameters;
              if (aflags.KBIN_RUN_AFLOWIN) {
                aus << " -name \"" << _AFLOWIN_ << "\" ";
              }
              // if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
              // if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
              // if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)  aus << " -name \"" << _AFLOWIN_ << "\" ";   // IS THIS CORRECT ?? CHECK !!!
              if (aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
                aus << " -name \"INCAR\" ";
              }
              aus << " | " << PRIORITY_GREP_STRING << " >> ";
              aus << FileNameSUBDIR << endl;
            }
          }
          // perform the command
          aurostd::execute(aus); // RESET  // RESET

          if (_VERBOSE_) {
            aus << "MMMMM  AFLOW: searching subdirectories [d2] " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
          }
          if (_VERBOSE_) {
            aurostd::PrintMessageStream(aus, XHOST.QUIET);
          }
          aurostd::string2tokens(aurostd::file2string(FileNameSUBDIR), vaflowin, "\n");
          for (i = 0; i < (int) vaflowin.size(); i++) {
            aurostd::StringSubstInPlace(vaflowin[i], _AFLOWIN_, "");
          }
          for (i = 0; i < (int) vaflowin.size(); i++) {
            aurostd::StringSubstInPlace(vaflowin[i], "INCAR", "");
          }
          if (STOP_DEBUG) {
            for (i = 0; i < (int) vaflowin.size(); i++) {
              cout << vaflowin[i] << endl;
            }
          }
          // RANDOMIZING priority
          if (vaflowin.size() > 1) { // only if I can poll
            // if(aurostd::substring2bool(vaflowin[0],"PRIORITY") || aurostd::substring2bool(vaflowin[0],"priority")) aurostd::random_shuffle(vaflowin);
          }
          // loaded up
          aurostd::RemoveFile(FileNameSUBDIR);
        }

        // FILE SPECIFIED
        if (aflags.AFLOW_PERFORM_FILE) {
          aus << "MMMMM  AFLOW: aflags.AFLOW_PERFORM_FILE==true" << " - " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
          vaflowin.clear();
          for (size_t i = 0; (i < vruns.size() && vaflowin.size() < maxcheck); i++) {
            if (Legitimate_aflowin(vruns[i], false, aus)) {
              vaflowin.push_back(vruns[i]);
            } // true puts too much verbosity
            // vaflowin.push_back(vruns[i]); // just load them up... they were checked before //OLD MUST RECHECH THEM as things change on the fly
          }
        }
        // NOW TIME OF SORTING/RANDOMIZING
        if (aflags.AFLOW_PERFORM_ORDER_SORT) { // SORT do something
          aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
          aurostd::sort(vaflowin);
        }
        if (aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
          aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
          aurostd::sort(vaflowin); // PN modification
          aurostd::rsort(vaflowin);
        }
        if (aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
          aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] - " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
          aurostd::random_shuffle(vaflowin);
        } // uses the std library but the seed is initialized in xrandom too
        aus << "MMMMM  " << aurostd::PaddedPOST("Legitimate VAFLOWIN = " + aurostd::utype2string(vaflowin.size()), 40) << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        //      aus << "MMMMM  Legitimate VAFLOWIN = " << vaflowin.size() << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);

        //	for(size_t i=0;i<vaflowin.size();i++) cerr << i << " " << vaflowin[i] << endl;

        // clean AFLOWIN // SAFETY
        for (size_t i = 0; i < vaflowin.size(); i++) {
          aurostd::StringSubstInPlace(vaflowin[i], _AFLOWIN_, "");
        }

        if (MULTI_DEBUG) {
          aus << "MMMMM  SIZE vaflowin=" << vaflowin.size() << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
        }
        if (MULTI_DEBUG) {
          for (size_t i = 0; i < vaflowin.size(); i++) {
            aus << "MMMMM  vaflowin[i]=" << vaflowin[i] << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(aus, XHOST.QUIET);
          }
        }
        // cerr << vaflowin.at(10) << " " << vaflowin.at(20) << " " << vaflowin.at(30) << endl;
        // if(aurostd::substring2bool(vaflowin[0],"PRIORITY") || aurostd::substring2bool(vaflowin[0],"priority")) aurostd::random_shuffle(vaflowin);

        FOUND = false;
        for (size_t i = 0; i < vaflowin.size() && !FOUND; i++) {
          aflags.Directory = "NULL";
          if (Legitimate_aflowdir(vaflowin[i], aflags)) {
            aflags.Directory = vaflowin[i];
            FOUND = true;
          }
        }
        // exiting if STOP_DEBUG
        if (STOP_DEBUG) {
          cout << "aflow_kbin.cpp: STOP_DEBUG" << endl;
          return 1;
        } // CO20200624 - previously exit

        // FOUND SOMETHING
        if (FOUND == false) {
          EMPTY = true;
          if (AFLOW_PTHREADS::FLAG) {
            aus << "MMMMM  AFLOW: MULTI-THREADED: FLUSHING PTHREADS - " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(aus, XHOST.QUIET);
            for (ithread = 0; ithread < AFLOW_PTHREADS::MAX_PTHREADS; ithread++) {
              if (AFLOW_PTHREADS::vpthread_busy[ithread]) {
                aus << "MMMMM  AFLOW: MULTI-THREADED: Flushing   pthread=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
                aurostd::PrintMessageStream(aus, XHOST.QUIET);
                pthread_join(AFLOW_PTHREADS::vpthread[ithread], nullptr);
              }
            }
          }
        }
        // again another check for LOCK, because NFS (network file system might be slow in concurrent seaches
        //     if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) && FOUND) {cerr << "AFLOW EXCEPTION on concurrent LOCK: " << aflags.Directory << endl; FOUND=false;}
        if (aurostd::DirectoryLocked(aflags.Directory, _AFLOWLOCK_) && FOUND) {
          aus << "AFLOW EXCEPTION on concurrent LOCK: " << aflags.Directory << endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
          FOUND = false;
        }

        //  ---------------------------------------------------------------------------- START RUNNING
        if (FOUND) {
          // again another check for LOCK, because NFS (network file system might be slow in concurrent seaches
          EMPTY = false;
          //  ---------------------------------------------------------------------------- KIN_RUN_AFLOWIN
          if (aflags.KBIN_RUN_AFLOWIN) {
            //	  bool PHONONS;
            //  -------------------------------------------------------------------------- KIN_RUN_AFLOWIN multithreaded
            // 	  bool found;
            // 	  for(size_t ii=0;ii<qaflags.size()&&!found;ii++)
            // 	    found=(qaflags[ii].Directory==aflags.Directory);       // look in all the list of operations
            // 	  if(found==false) {                                 // new operation, generate and save it
            // 	    qaflags.push_back(aflags);
            // 	  }
            //	  cerr << qaflags.size() << endl;
            if (AFLOW_PTHREADS::FLAG) {
              // there is something to run in aflags.
              // wait and put in ithread there is the number of the thread
              free_thread = AFLOW_PTHREADS::Wait_Available_Free_Threads(ithread, _VERBOSE_); // WAIT A WHILE !!
              if (free_thread) {
                aus << "MMMMM  AFLOW: Found subdirectory to run " << aflags.Directory << " - " << Message(__AFLOW_FILE__, aflags) << endl;
                aus << "MMMMM  AFLOW: MULTI-THREADED: Starting    pthread_free=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
                aurostd::PrintMessageStream(aus, XHOST.QUIET);
                aflags.AFLOW_PTHREADS_NUMBER = ithread;
                KBIN::RUN_Directory_PTHREADS(aflags);
                RUN_times++;
                if (XHOST.AFLOW_RUNXflag) {
                  aus << "MMMMM  AFLOW: RUNFX finished running " << RUN_times << "/" << XHOST.AFLOW_RUNXnumber << " " << aflags.Directory << " - " << Message(__AFLOW_FILE__, aflags) << endl;
                  aurostd::PrintMessageStream(aus, XHOST.QUIET);
                }
                if (XHOST.AFLOW_RUNXflag && RUN_times == XHOST.AFLOW_RUNXnumber) {
                  EMPTY = true; // force to end if RUXN reached
                }
              }
            }
            //  -------------------------------------------------------------------------- KIN_RUN_AFLOWIN normal
            if (!AFLOW_PTHREADS::FLAG) {
              KBIN::RUN_Directory(aflags);
              RUN_times++;
              if (XHOST.AFLOW_RUNXflag) {
                aus << "MMMMM  AFLOW: RUNFX finished running " << RUN_times << "/" << XHOST.AFLOW_RUNXnumber << " " << aflags.Directory << " - " << Message(__AFLOW_FILE__, aflags) << endl;
                aurostd::PrintMessageStream(aus, XHOST.QUIET);
              }
              if (XHOST.AFLOW_RUNXflag && RUN_times == XHOST.AFLOW_RUNXnumber) {
                EMPTY = true; // force to end if RUXN reached
              }
            }
          }
          //  ---------------------------------------------------------------------------- KBIN_GEN_VASP_FROM_AFLOWIN normal
          //	if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     KBIN::RUN_Directory(argv,aflags);  // IS THIS CORRECT ?? CHECK !!!
          //  ---------------------------------------------------------------------------- KBIN_GEN_AFLOWIN_FROM_VASP normal
          if (aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
            KBIN::GenerateAflowinFromVASPDirectory(aflags);
          }
        }
        if (XHOST.vflag_aflow.flag("LOOP") && EMPTY && XHOST.AFLOW_RUNXflag == false) {
          EMPTY = false;
          aus << "MMMMM  AFLOW: waiting for new subdirectories: " << (int) _KBIN_LOOP_SLEEP_ / 60 << "mins ";
          aus << " - " << XHOST.hostname << " - " << aflow_get_time_string() << endl; // endl;
          aurostd::PrintMessageStream(aus, XHOST.QUIET);
          aurostd::Sleep(_KBIN_LOOP_SLEEP_);
        }
      }
      aus << "MMMMM  AFLOW: no more subdirectories to run " << " - " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }

    // ------------------------------------------------------------------------------------------------------------------------------------

    return 1;
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::MPI_Extract
// ***************************************************************************
// This function extracts from _AFLOWIN_ the parameters for MPI run
namespace KBIN {
  void MPI_Extract(string AflowIn, ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags) {
    ostringstream aus;
    bool Kmpi = true;
    kflags.KBIN_MPI_NCPUS = 0;
    aus << "00000  [AFLOW_MODE_MPI] found in " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__, aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    // get (integer) kflags.KBIN_MPI_NCPUS

    if (aflags.AFLOW_GLOBAL_NCPUS < 1) {
      if (Kmpi && !aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=", true) && !aflags.AFLOW_MACHINE_LOCAL.flag()) { // DEFAULT NO CPU SPECIFIED
        kflags.KBIN_MPI_NCPUS = MPI_NCPUS_DEFAULT;
        aus << "00000  MESSAGE MPI: NCPUS=NNNN is missing, taking NCPUS=" << kflags.KBIN_MPI_NCPUS << "  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        Kmpi = false;
      }
      if (Kmpi && (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=MAX", true) || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=AUTO", true) || aflags.AFLOW_MACHINE_LOCAL.flag())) { // DEFAULT NCPUS=MAX
        kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores;
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=MAX", true)) {
          kflags.KBIN_MPI_NCPUS_STRING = "MAX";
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=AUTO", true)) {
          kflags.KBIN_MPI_NCPUS_STRING = "AUTO";
        }

        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // CO
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // with DUKE_BETA force NCPUS from QUEUE
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // CO
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
          kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores; // with DUKE_MATERIALS force NCPUS
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
          kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores; // with DUKE_AFLOWLIB force NCPUS
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
          kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores; // with DUKE_HABANA force NCPUS
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // with MACHINE001; DX added 20190509
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // with MACHINE002; DX added 20190509
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE003")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // with MACHINE003; DX added 20201005
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE004")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN; // with MACHINE004; DX added 20211011
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) {
          kflags.KBIN_MPI_NCPUS = XHOST.PBS_NUM_PPN;
        }; // DX20190107 - CMU EULER // with CMU_EULER force NCPUS //DX20181113
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) {
          kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores; // MACHINE2 has only NCPUS //CO20181113
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) {
          kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores; // MACHINE1 has only NCPUS //CO20181113
        }
        if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
          kflags.KBIN_MPI_NCPUS = XHOST.SLURM_NTASKS; // with FULTON_MARYLOU force NCPUS
        }
        if (kflags.KBIN_MPI_NCPUS < 1) {
          kflags.KBIN_MPI_NCPUS = XHOST.CPU_Cores; // SAFE
        }

        aus << "00000  MESSAGE MPI: found NCPUS=MAX  NCPUS=" << kflags.KBIN_MPI_NCPUS << " " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        Kmpi = false;
      }
      if (Kmpi && aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=", true)) { // DEFAULT NCPUS=XXX
        kflags.KBIN_MPI_NCPUS = aurostd::substring2utype<int>(AflowIn, "[AFLOW_MODE_MPI_MODE]NCPUS=", true);
        if (kflags.KBIN_MPI_NCPUS > 0) {
          kflags.KBIN_MPI_NCPUS_STRING = aurostd::utype2string<int>(kflags.KBIN_MPI_NCPUS); // ME20181113
          aus << "00000  MESSAGE MPI: found NCPUS=" << kflags.KBIN_MPI_NCPUS << " " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        Kmpi = false;
      }
    } else {
      kflags.KBIN_MPI_NCPUS = aflags.AFLOW_GLOBAL_NCPUS;
      aus << "00000  MESSAGE MPI: NCPUS is overriden, taking NCPUS=" << kflags.KBIN_MPI_NCPUS << "  " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    if (kflags.KBIN_MPI_NCPUS < 1) {
      kflags.KBIN_MPI_NCPUS = 1; // DEFAULT NCPUS=troubles
    }

    if (kflags.KBIN_MPI) {
      // get (string) kflags.KBIN_MPI_START
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]START=", true)) {
        kflags.KBIN_MPI_START = MPI_START_DEFAULT;
        aus << "00000  MESSAGE MPI: START string is missing, taking START=\"" << kflags.KBIN_MPI_START << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } else {
        kflags.KBIN_MPI_START = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_MPI_MODE]START=", 1, true), '"');
        aus << "00000  MESSAGE MPI: found START=\"" << kflags.KBIN_MPI_START << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      // get (string) kflags.KBIN_MPI_STOP
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]STOP=", true)) {
        kflags.KBIN_MPI_STOP = MPI_STOP_DEFAULT;
        aus << "00000  MESSAGE MPI: STOP string is missing, taking STOP=\"" << kflags.KBIN_MPI_STOP << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } else {
        kflags.KBIN_MPI_STOP = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_MPI_MODE]STOP=", 1, true), '"');
        aus << "00000  MESSAGE MPI: found STOP=\"" << kflags.KBIN_MPI_STOP << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      // get (string) kflags.KBIN_MPI_COMMAND
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]COMMAND=", true)) {
        kflags.KBIN_MPI_COMMAND = MPI_COMMAND_DEFAULT;
        aus << "00000  MESSAGE MPI: COMMAND string is missing, taking COMMAND=\"" << kflags.KBIN_MPI_COMMAND << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } else {
        kflags.KBIN_MPI_COMMAND = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_MPI_MODE]COMMAND=", 1, true), '"');
        aus << "00000  MESSAGE MPI: found COMMAND=\"" << kflags.KBIN_MPI_COMMAND << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      kflags.KBIN_MPI_AUTOTUNE = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]AUTOTUNE", true);
      if (kflags.KBIN_MPI_AUTOTUNE) {
        aus << "00000  MESSAGE MPI: found AUTOTUNE option " << Message(__AFLOW_FILE__, aflags) << endl;
        aus << "00000  MESSAGE MPI: input files WILL be auto-tuned for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } else {
        aus << "00000  MESSAGE MPI: AUTOTUNE option NOT found " << Message(__AFLOW_FILE__, aflags) << endl;
        aus << "00000  MESSAGE MPI: input files MUST be appropriate for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      // get (string) kflags.KBIN_MPI_BIN
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]BINARY=", true)) {
        kflags.KBIN_MPI_BIN = DEFAULT_VASP_MPI_BIN;
        aus << "00000  MESSAGE MPI: BINARY string is missing, taking BIN=\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } else {
        kflags.KBIN_MPI_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_MPI_MODE]BINARY=", 1, true), '"');
        aus << "00000  MESSAGE MPI: found BINARY=\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      // ME20190107 - Grab the serial binary to propagate into child aflow.in files
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY]")) {
        kflags.KBIN_SERIAL_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY]");
      } else if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY=")) {
        kflags.KBIN_SERIAL_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY="), ']');
      }
      aus << "00000  MESSAGE MPI: Overriding BINARY=\"" << kflags.KBIN_BIN << "\" to BINARY =\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      kflags.KBIN_BIN = kflags.KBIN_MPI_BIN;

      // get (string) kflags.KBIN_MPI_OPTIONS
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI_MODE]OPTIONS=", true)) {
        kflags.KBIN_MPI_OPTIONS = VASP_OPTIONS_MPI_DEFAULT;
        aus << "00000  MESSAGE MPI: OPTIONS string is missing, taking OPTIONS=\"" << kflags.KBIN_MPI_OPTIONS << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } else {
        kflags.KBIN_MPI_OPTIONS = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_MPI_MODE]OPTIONS=", 1, true), '"');
        aus << "00000  MESSAGE MPI: found OPTIONS=\"" << kflags.KBIN_MPI_OPTIONS << "\"  " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
    }
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::StartStopCheck
// ***************************************************************************
namespace KBIN {
  void StartStopCheck(const string& AflowIn, string str1, string str2, bool& flag, bool& flagS) {
    flag = aurostd::substring2bool(AflowIn, str1) || aurostd::substring2bool(AflowIn, str2);
    flagS = (aurostd::substring2bool(AflowIn, str1 + "START") && aurostd::substring2bool(AflowIn, str1 + "STOP")) || (aurostd::substring2bool(AflowIn, str2 + "_START") && aurostd::substring2bool(AflowIn, str2 + "_STOP"));
    if (flagS) {
      flag = false;
    }
  }
} // namespace KBIN

namespace KBIN {
  void StartStopCheck(const string& AflowIn, string str1, bool& flag, bool& flagS) {
    flag = aurostd::substring2bool(AflowIn, str1);
    flagS = aurostd::substring2bool(AflowIn, str1 + "START") && aurostd::substring2bool(AflowIn, str1 + "STOP");
    if (flagS) {
      flag = false;
    }
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::MoveRun2NewDirectory() //DX20210901
// ***************************************************************************
// Moves an aflow run (i.e., aflow.in) to a new directory and adds a LOCK
// to the original directory to prevent aflow from re-running the system.
// Some machines have scrubbers on certain filesystems/workspaces that remove
// old/untouched files, i.e., removes/deletes aflow.ins that are waiting to be
// run. To avoid this, we store the aflow.in in a filesystem/workspace "safe"
// from the scrubber, and move the aflow.in to the "run" filesystem/workspace
// once it has been selected from the aflow daemon.
// This is needed for machine001, machine002, machine003
// SD20220319 - changed function to return a bool
namespace KBIN {
  bool MoveRun2NewDirectory(_aflags& aflags, const string& subdirectory_orig, const string& subdirectory_new) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    ostringstream message;

    // ---------------------------------------------------------------------------
    // create lock immediately // CO20210901
    if (!aurostd::LinkFile(aflags.Directory + "/" + _AFLOWIN_, aflags.Directory + "/" + _AFLOWLOCK_ + _LOCK_LINK_SUFFIX_, false)) {
      return false;
    } // create LOCK link //SD20220224
    message << "MMMMM  created link \"" << aflags.Directory + "/" + _AFLOWLOCK_ + _LOCK_LINK_SUFFIX_ << "\"" << Message(__AFLOW_FILE__, aflags, "user,host,time") << endl;
    aurostd::PrintMessageStream(message, XHOST.QUIET);
    message.clear();
    message.str(std::string());

    // ---------------------------------------------------------------------------
    // Changing the run directory from the "original" to a "new" directory
    const string directory_orig = aflags.Directory;
    aurostd::StringSubstInPlace(aflags.Directory, subdirectory_orig, subdirectory_new);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " original full directory " << directory_orig << endl;
      cerr << __AFLOW_FUNC__ << " changing subdirectory " << subdirectory_orig << " to " << subdirectory_new << endl;
      cerr << __AFLOW_FUNC__ << " new full directory " << aflags.Directory << endl;
    }

    // ---------------------------------------------------------------------------
    // make new directory
    message << "MMMMM create directory \"" << aflags.Directory << "\"" << Message(__AFLOW_FILE__, aflags, "user,host,time") << endl;
    aurostd::PrintMessageStream(message, XHOST.QUIET);
    message.clear();
    message.str(std::string());
    aurostd::DirectoryMake(aflags.Directory);

    // ---------------------------------------------------------------------------
    // copy aflow.in to new directory
    message << "MMMMM copying aflowin from \"" << directory_orig << "\"" << Message(__AFLOW_FILE__, aflags, "user,host,time") << endl;
    aurostd::PrintMessageStream(message, XHOST.QUIET);
    message.clear();
    message.str(std::string());
    aurostd::CopyFile(directory_orig + "/" + _AFLOWIN_, aflags.Directory + "/" + _AFLOWIN_);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " new full directory " << aflags.Directory << endl;
    }
    return true;
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::RUN_Directory
// ***************************************************************************
namespace KBIN {
  void RUN_Directory(_aflags& aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    const bool LDEBUG = (false || XHOST.DEBUG);
    ostringstream aus;

    // ---------------------------------------------------------------------------
    // Move aflow run (i.e., aflow.in) to a new directory and add a LOCK
    // to the original directory to prevent machine scrubbers from removing
    // DX20210901
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001") || aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002") || aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE003") ||
        aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE004")) {
      const std::string subdirectory_orig = aurostd::getenv2string("HOME"); // $HOME    : environment variable pointing to "home" filesystem (specific to machine001/002/003)
      const std::string subdirectory_new = aurostd::getenv2string("WORKDIR"); // $WORKDIR : environment variable pointing to "work" filesystem (specific to machine001/002/003)
      if (subdirectory_new.empty() || subdirectory_orig.empty()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Required environment variables ($HOME, $WORKDIR) are not completely defined", _RUNTIME_INIT_);
      }
      if (!KBIN::MoveRun2NewDirectory(aflags, subdirectory_orig, subdirectory_new)) {
        return;
      }
    }

    // string::size_type sub_size1,sub_size2;
    string AflowIn;
    string agllock;
    const string::iterator pos;
    bool Krun = true;
    //  int i;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (aflags.Directory.empty() || (aflags.Directory[0] != '/' && aflags.Directory[0] != '.' && aflags.Directory[0] != ' ')) {
      aflags.Directory = "./" + aflags.Directory;
    }

    if (!Legitimate_aflowdir(aflags.Directory, aflags) || !aurostd::LinkFile(aflags.Directory + "/" + _AFLOWIN_, aflags.Directory + "/" + _AFLOWLOCK_ + _LOCK_LINK_SUFFIX_, false)) {
      return;
    } // create LOCK link //SD20220224
    // ***************************************************************************
    // RESET LOCK
    ofstream FileLOCK;
    const string FileNameLOCK = aflags.Directory + "/" + _AFLOWLOCK_;
    //	FileLOCK.open(FileNameLOCK.c_str(),std::ios::out);
    FileLOCK.open(FileNameLOCK.c_str(), std::ios::app);
    // ***************************************************************************
    // START DIRECTORY
    aus << "XXXXX  KBIN DIRECTORY BEGIN (aflow" << string(AFLOW_VERSION) << ")  " << Message(__AFLOW_FILE__, aflags) << endl;
    //	aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
    aus << "XXXXX  KBIN XHOST.CPU_Model : " << XHOST.CPU_Model << "" << endl; // << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "XXXXX  KBIN XHOST.CPU_Cores : " << XHOST.CPU_Cores << "" << endl; // << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "XXXXX  KBIN XHOST.CPU_MHz   : " << XHOST.CPU_MHz << "" << endl; // << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "XXXXX  KBIN XHOST.RAM_GB    : " << XHOST.RAM_GB << "" << endl; // << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileLOCK, aus, XHOST.QUIET);
    // ***************************************************************************
    // FLUSH & REOPEN to avoid double writing
    FileLOCK.flush();
    FileLOCK.close();
    FileLOCK.open(FileNameLOCK.c_str(), std::ios::app);
    aurostd::RemoveFile(aflags.Directory + "/" + _AFLOWLOCK_ + _LOCK_LINK_SUFFIX_); // Remove LOCK link //SD20220207
    // ***************************************************************************
    // NOW Digest AFLOWIN
    ifstream FileAFLOWIN;
    const string FileNameAFLOWIN = aflags.Directory + "/" + _AFLOWIN_; // SD20220224
    FileAFLOWIN.open(FileNameAFLOWIN.c_str(), std::ios::in); // SD20220224
    AflowIn.clear();
    char c;
    while (FileAFLOWIN.get(c)) {
      if (c != '\0') {
        AflowIn += c;
      } // READ _AFLOWIN_ and put into AflowIn //DX20190125 - remove null bytes from AflowIn
    }
    FileAFLOWIN.clear();
    FileAFLOWIN.seekg(0);
    AflowIn = aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    vector<string> vAflowIn;
    aurostd::string2vectorstring(AflowIn, vAflowIn); // CO20181226

    _kflags kflags = KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn, FileLOCK, aflags); // CO20200624 - made separate function for getting kflags

    if (kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) { // [AFLOW_MODE_PRESCRIPT] construction
      aurostd::stringstream2file(kflags.AFLOW_MODE_PRESCRIPT, string(aflags.Directory + "/" + DEFAULT_AFLOW_PRESCRIPT_COMMAND));
      aurostd::Chmod(0755, string(aflags.Directory + "/" + DEFAULT_AFLOW_PRESCRIPT_COMMAND));
    }
    if (kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) { // [AFLOW_MODE_POSTSCRIPT] construction
      aurostd::stringstream2file(kflags.AFLOW_MODE_POSTSCRIPT, string(aflags.Directory + "/" + DEFAULT_AFLOW_POSTSCRIPT_COMMAND));
      aurostd::Chmod(0755, string(aflags.Directory + "/" + DEFAULT_AFLOW_POSTSCRIPT_COMMAND));
    }
    // ************************************************************************************************************************************
    // ALIEN MODE
    if (Krun && kflags.AFLOW_MODE_ALIEN) {
      // ***************************************************************************
      // ALIEN MODE  // must contain EMAIL perform
      if (Krun) {
        Krun = (Krun && ALIEN::Run_Directory(FileLOCK, aflags, kflags));
      }
      // ***************************************************************************
      // COMPRESS
      if (Krun && kflags.KZIP_COMPRESS) {
        Krun = (Krun && KBIN::CompressDirectory(aflags));
      }
    }
    // ************************************************************************************************************************************
    // MATLAB MODE
    if (Krun && kflags.AFLOW_MODE_MATLAB) {
      if (Krun) {
        aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
        Krun = (Krun && KBIN_MATLAB_Directory(FileLOCK, aflags, kflags));
      }
      // ***************************************************************************
      // COMPRESS
      if (Krun && kflags.KZIP_COMPRESS) {
        Krun = (Krun && KBIN::CompressDirectory(aflags));
      }
      Krun = false;
    }
    // ************************************************************************************************************************************
    // AIMS MODE
    if (Krun && kflags.AFLOW_MODE_AIMS) {
      // ***************************************************************************
      // AIMS MODE  // must contain EMAIL perform
      if (Krun) {
        Krun = (Krun && KBIN::AIMS_Directory(FileLOCK, aflags, kflags));
      }
      // ***************************************************************************
      // COMPRESS
      if (Krun && kflags.KZIP_COMPRESS) {
        Krun = (Krun && KBIN::CompressDirectory(aflags));
      }
    }
    // ************************************************************************************************************************************
    // ************************************************************************************************************************************
    // VASP MODE
    if (Krun && kflags.AFLOW_MODE_VASP) {
      // ***************************************************************************
      // VASP MODE  // must contain EMAIL perform
      if (Krun) {
        Krun = (Krun && KBIN::VASP_Directory(FileLOCK, aflags, kflags));
      }
      // ***************************************************************************
      // COMPRESS
      if (Krun && kflags.KZIP_COMPRESS) {
        // cerr << aurostd::execute2string("ls -las "+aflags.Directory) << endl;
        Krun = (Krun && KBIN::CompressDirectory(aflags));
        // cerr << aurostd::execute2string("ls -las "+aflags.Directory) << endl;
      }
    }
    // ************************************************************************************************************************************
    // MATLAB MODE
    if (Krun && kflags.KBIN_PHONONS_CALCULATION_FROZSL && !kflags.AFLOW_MODE_VASP) {
      // PRESCRIPT
      if (Krun && (kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)) {
        KBIN::RUN_DirectoryScript(aflags, DEFAULT_AFLOW_PRESCRIPT_COMMAND, DEFAULT_AFLOW_PRESCRIPT_OUT);
      }
      // POSTSCRIPT
      if (Krun && (kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)) {
        KBIN::RUN_DirectoryScript(aflags, DEFAULT_AFLOW_POSTSCRIPT_COMMAND, DEFAULT_AFLOW_POSTSCRIPT_OUT);
      }
    }
    // ***************************************************************************
    // FINALIZE AFLOWIN
    // DONE, turn OFF the flag
    kflags.AFLOW_MODE_VASP = false;
    // ***************************************************************************
    // NFS cache-cleaning HACK, opendir() and closedir() to invalidate cache
    // https://stackoverflow.com/questions/8311710/nfs-cache-cleaning-command
    // CO20171106
    vector<string> vfiles_dummyls;
    aurostd::DirectoryLS(aflags.Directory, vfiles_dummyls);
    // ***************************************************************************
    // WRITE END
    aurostd::string2file("AFLOW calculation complete" + string(Message(__AFLOW_FILE__, aflags) + "\n"), string(aflags.Directory + "/" + DEFAULT_AFLOW_END_OUT));
    // ***************************************************************************
    // MAKE READEABLE
    agllock = aflags.Directory + "/" + _AFLOWLOCK_;
    aurostd::CompressFileExist(agllock, agllock);
    aurostd::Chmod(0664, string(agllock));
    aurostd::CopyFile(agllock, agllock + "." + string(AFLOW_VERSION));
    aurostd::ChmodRecursive(0777, 0664, aflags.Directory);

    FileAFLOWIN.clear();
    FileAFLOWIN.close();
  };
} // namespace KBIN

// *******************************************************************************************
namespace KBIN {
  void AFLOW_RUN_Directory(const _aflags& aflags) {
    aurostd::execute(XHOST.command("aflow") + " --run=1 --DIRECTORY=" + aflags.Directory); // run it OUTSIDE
    // this is cool as if the particula program explodes, there is still aflow running
  }
} // namespace KBIN

// *******************************************************************************************
namespace KBIN {
  void RUN_DirectoryScript(const _aflags& aflags, const string& script, const string& output) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << "cd " << aflags.Directory << endl;
    aus << "chmod 755 " << script << endl;
    aus << "./" << script << " >> " << output << endl;
    aurostd::execute(aus);
  }
} // namespace KBIN

// *******************************************************************************************
// KBIN::CompressDirectory
// *******************************************************************************************
namespace KBIN {
  /// @brief compress AFLOW run directory
  /// @param directory to compress
  /// @authors
  /// @mod{ME,20210927,created function}
  /// @mod{SD,20240503,rewritten using filesystem}
  /// @mod{SG,20240531,patched filename bug}
  bool CompressDirectory(const string& directory) {
    const std::set<string> skip_files = {_AFLOWIN_, _AFLOWLOCK_, DEFAULT_AFLOW_END_OUT, KBIN_SUBDIRECTORIES, "aflow.in", "aflow.end.out", "LOCK", "SKIP"};
    const std::filesystem::path path(aurostd::CleanFileName(directory));
    if (!std::filesystem::exists(path)) {
      return false;
    }

    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(path)) {
      if (dir_entry.is_regular_file()) {
        if (aurostd::IsCompressed(dir_entry.path().string())) {
          continue;
        }
        if (skip_files.count(dir_entry.path().filename().string())) {
          continue;
        }
        if (dir_entry.path().filename().string()[0] == '.') {
          continue;
        }
        if (std::regex_match(dir_entry.path().filename().string(), std::regex("aflow_.*.in"))) {
          continue;
        }
        aurostd::CompressFile(dir_entry.path().string());
      }
    }

    return true;
  }

  bool CompressDirectory(const _aflags& aflags) {
    return CompressDirectory(aflags.Directory);
  }
} // namespace KBIN

// *******************************************************************************************
// KBIN::DecompressDirectory
// *******************************************************************************************
namespace KBIN {
  /// @brief decompress AFLOW run directory
  /// @param directory to decompress
  /// @authors
  /// @mod{SD,20250318,created function}
  bool DecompressDirectory(const string& directory) {
    const std::filesystem::path path(aurostd::CleanFileName(directory));
    if (!std::filesystem::exists(path)) {
      return false;
    }

    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(path)) {
      if (dir_entry.is_regular_file()) {
        if (aurostd::IsCompressed(dir_entry.path().string())) {
          aurostd::DecompressFile(dir_entry.path().string());
        }
      }
    }

    return true;
  }

  bool DecompressDirectory(const _aflags& aflags) {
    return DecompressDirectory(aflags.Directory);
  }
} // namespace KBIN

// *******************************************************************************************
// KBIN::Clean
// *******************************************************************************************
namespace KBIN {
  void Clean(const _aflags& aflags) {
    const aurostd::xoption opts_clean;
    return KBIN::Clean(aflags, opts_clean);
  } // AFLOW_FUNCTION_IMPLEMENTATION
  void Clean(const _aflags& aflags, const aurostd::xoption& opts_clean) { // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210901
    return KBIN::Clean(aflags.Directory, opts_clean);
  }
} // namespace KBIN

namespace KBIN {
  /// @brief clean AFLOW run directory
  /// @param directory to compress
  /// @param opts_clean command line options
  /// @authors
  /// @mod{CO,20210901,created function}
  /// @mod{SD,20240417,rewritten using filesystem}
  void Clean(const string& directory, const aurostd::xoption& opts_clean) {
    const std::set<string> skip_files = {"NOCLEAN", "noclean", "NOCLEAR", "noclear"};
    const std::set<string> input_files = {_AFLOWIN_, _AFLOWIN_AGL_DEFAULT_, _AFLOWIN_AEL_DEFAULT_, "aflow.in"};
    const std::filesystem::path path(aurostd::CleanFileName(directory));
    if (!std::filesystem::exists(path)) {
      return;
    }

    vector<std::filesystem::directory_entry> vfile;
    std::set<string> sfile;

    // Collect the files as a set and vector
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(path)) {
      if (dir_entry.is_regular_file() || dir_entry.is_directory()) {
        sfile.insert(dir_entry.path().filename().string());
        vfile.push_back(dir_entry);
      }
    }

    // Check if we skip
    for (const string& skip_file : skip_files) {
      if (sfile.count(skip_file)) {
        return;
      }
    }

    // Check if input exists
    bool valid = false;
    for (const string& input_file : input_files) {
      if (sfile.count(input_file)) {
        valid = true;
        break;
      }
    }
    if (!valid) {
      return;
    }

    // Save CONTCAR
    if (opts_clean.flag("SAVE_CONTCAR")) {
      KBIN::VASP_CONTCAR_Save(directory);
    } else if (opts_clean.flag("SAVE_CONTCAR_OUTCAR_COMPLETE")) {
      _xvasp xvasp;
      xvasp.Directory = path.parent_path().string();
      _aflags aflags;
      ofstream FileMESSAGE;
      if (KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, false)) {
        KBIN::VASP_CONTCAR_Save(directory);
      }
    }

    // Clean directory
    for (const std::filesystem::directory_entry& file : vfile) {
      if (file.is_regular_file() && input_files.find(file.path().filename().string()) == input_files.end()) {
        std::filesystem::remove(file.path());
      } else if (file.is_directory()) {
        std::filesystem::remove_all(file.path());
      }
    }
  }

  void Clean(const string& _directory) {
    const aurostd::xoption opts_clean;
    return KBIN::Clean(_directory, opts_clean);
  }
} // namespace KBIN

// *******************************************************************************************
// KBIN::XClean
// *******************************************************************************************
namespace KBIN {
  void XClean(string options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");
    if (!tokens.empty()) {
      init::ErrorOption(options, __AFLOW_FUNC__, "aflow --xclean");
    }

    vector<string> vcheck1;
    aurostd::string2tokens(string("OUTCAR.static,OUTCAR.relax2,OUTCAR.relax1"), vcheck1, ",");
    vector<string> vcheck2;
    aurostd::string2tokens(string("OUTCAR,OUTCAR,OUTCAR"), vcheck2, ",");

    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      vcheck1.push_back("OUTCAR.relax1" + XHOST.vext[iext]);
    }

    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      vcheck2.push_back("OUTCAR.relax2" + XHOST.vext[iext]);
    }

    vector<string> vfile;
    const bool test = false;

    cout << __AFLOW_FUNC__ << " checking missing " << "OUTCAR*" << " with " << _AFLOWLOCK_ << endl; // check OUTCAR.static
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find") + " ./ -name " + _AFLOWLOCK_), vfile);
    for (size_t j = 0; j < vfile.size(); j++) {
      aurostd::StringSubstInPlace(vfile[j], _AFLOWLOCK_, "");
      if (!aurostd::FileExist(vfile[j] + "OUTCAR") && !aurostd::CompressFileExist(vfile[j] + "OUTCAR.relax1")) {
        cout << __AFLOW_FUNC__ << " cleaning=" << vfile[j] << endl;
        if (!test) {
          KBIN::Clean(vfile[j]);
        }
      }
    }
    for (size_t i = 0; i < vcheck1.size(); i++) {
      cout << __AFLOW_FUNC__ << " checking missing " << vcheck2[i] << " with " << vcheck1[i] << endl; // check OUTCAR.static
      aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find") + " ./ -name " + vcheck1[i]), vfile);
      for (size_t j = 0; j < vfile.size(); j++) {
        aurostd::StringSubstInPlace(vfile[j], vcheck1[i], "");
        if (!aurostd::FileExist(vfile[j] + vcheck2[i])) {
          cout << __AFLOW_FUNC__ << " cleaning=" << vfile[j] << endl;
          if (!test) {
            KBIN::Clean(vfile[j]);
          }
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
  }
} // namespace KBIN

// *******************************************************************************************
// KBIN::XClean
// *******************************************************************************************

// ME20200219 - based on CO's code in old PhononCalculator
namespace KBIN {
  int get_NCPUS() {
    const _kflags kflags;
    return get_NCPUS(kflags);
  }

  int get_NCPUS(const _kflags& kflags) { // CO20220630 - rewritten
    int ncpus = 1;
    if (kflags.KBIN_MPI_NCPUS > 0) {
      ncpus = kflags.KBIN_MPI_NCPUS;
    } // aflow.in is lowest priority
    if (XHOST.vflag_control.flag("NUM_THREADS")) { // command line argument overrides
      const string ncpus_str = XHOST.vflag_control.getattachedscheme("NUM_THREADS");
      if (aurostd::isfloat(ncpus_str)) {
        ncpus = aurostd::string2utype<int>(ncpus_str);
      }
      if (ncpus_str == "MAX") {
        ncpus = MPI_NCPUS_MAX;
      }
    }
    if (ncpus < 1) {
      ncpus = 1;
    }
    return ncpus;
  }
} // namespace KBIN

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
