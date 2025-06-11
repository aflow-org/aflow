// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AFLOWLIB_LIBRARIES_CPP_
#define _AFLOWLIB_LIBRARIES_CPP_

#include "aflowlib/aflowlib_libraries.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <iterator>
#include <ostream>
#include <regex>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_hash.h"
#include "AUROSTD/aurostd_time.h"
#include "AUROSTD/aurostd_xcombos.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "aflowlib/aflowlib.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_bader.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_kvasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"
#include "interfaces/aflow_pthreads.h"
#include "modules/AEL/aflow_ael_elasticity.h" //CT20200713
#include "modules/AGL/aflow_agl_debye.h" //CT20200713
#include "modules/APL/aflow_apl.h"
#include "modules/CCE/aflow_cce.h" //CO20200624
#include "modules/POCC/aflow_pocc.h" //CO20200624
#include "modules/PROTOTYPES/aflow_anrl.h"
#include "modules/SYM/aflow_symmetry.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using std::cerr;
using std::cout;
using std::deque;
using std::endl;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::stringstream;
using std::vector;

vector<string> vLibrary_ALL;
vector<vector<string>> vLibrary_ALL_tokens;
vector<string> vLibrary_ICSD;
vector<vector<string>> vLibrary_ICSD_tokens;
vector<string> vLibrary_LIB0;
vector<vector<string>> vLibrary_LIB0_tokens;
vector<string> vLibrary_LIB1;
vector<vector<string>> vLibrary_LIB1_tokens;
vector<string> vLibrary_LIB2;
vector<vector<string>> vLibrary_LIB2_tokens;
vector<string> vLibrary_LIB3;
vector<vector<string>> vLibrary_LIB3_tokens;
vector<string> vLibrary_LIB4;
vector<vector<string>> vLibrary_LIB4_tokens;
vector<string> vLibrary_LIB5;
vector<vector<string>> vLibrary_LIB5_tokens;
vector<string> vLibrary_LIB6;
vector<vector<string>> vLibrary_LIB6_tokens;
vector<string> vLibrary_LIB7;
vector<vector<string>> vLibrary_LIB7_tokens;
vector<string> vLibrary_LIB8;
vector<vector<string>> vLibrary_LIB8_tokens;
vector<string> vLibrary_LIB9;
vector<vector<string>> vLibrary_LIB9_tokens;

bool AFLOWLIB_VERBOSE = true; // false;
#define _EPSILON_COMPOSITION_ DEFAULT_POCC_STOICH_TOL // 0.001

#define USE_AFLOW_SG
// #define USE_PLATON_SG

#define AFLOW_CORE_TEMPERATURE_LIB2RAW 46.0

#define RELAX_MAX 10  // CO20200829

// ******************************************************************************************************************************************************

// ***************************************************************************
// XPLUG/XUPDATE STUFF
namespace aflowlib {  // CO20210302
  //--------------------------------------------------------------------------------
  // constructor
  //--------------------------------------------------------------------------------
  ARun::ARun(const string& dir, ostream& oss) : xStream(oss), m_initialized(false) {
    initialize(dir);
  }
  ARun::ARun(const string& dir, ofstream& FileMESSAGE, ostream& oss) : xStream(FileMESSAGE, oss), m_initialized(false) {
    initialize(dir);
  }
  ARun::ARun(const ARun& b) : xStream(*b.getOFStream(), *b.getOSS()) {
    copy(b);
  } // copy PUBLIC

  ARun::~ARun() {
    xStream::free();
    free();
  }

  const ARun& ARun::operator=(const ARun& other) {
    if (this != &other) {
      copy(other);
    }
    return *this;
  }

  void ARun::clear() {
    free();
  }  // clear PUBLIC

  void ARun::free() {
    m_initialized = false;
    m_dir.clear();
    m_aflowin.clear();
    m_LOCK.clear();
    m_completed = false;
  }

  void ARun::copy(const ARun& b) {  // copy PRIVATE
    xStream::copy(b);
    m_initialized = b.m_initialized;
    m_dir = b.m_dir;
    m_aflowin = b.m_aflowin;
    m_LOCK = b.m_LOCK;
    m_completed = b.m_completed;
  }

  bool ARun::initialize(const string& dir) {
    m_dir = dir;
    m_initialized = (!m_dir.empty());

    const vector<string> vLOCKs;
    const vector<string> vaflowins = {"aflow.in", _AFLOWIN_AGL_DEFAULT_, ""};
    return true;
  }

} // namespace aflowlib

// ***************************************************************************

// ***************************************************************************
// aflowlib::TokenPresentAFLOWLIB
// ***************************************************************************
namespace aflowlib {
  /// @brief tests existence of substring query in line
  /// @param line string to search in
  /// @param query string to search for
  /// @return true if query is present in line
  /// @authors
  /// @mod{ST,20241021,created doxy}
  bool TokenPresentAFLOWLIB(const string& line, const string& query) {
    return aurostd::substring2bool(line, query, true);
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::XLIB2RAW_CheckProjectFromDirectory
// ***************************************************************************
namespace aflowlib {
  /// @brief find the appropriate aflow project library from the given directory
  /// @param directory the directory to match against the aflow project libraries
  /// @return the string of the matching project library
  /// @authors
  /// @mod{ST,20241022,created doxy\, use func macro\, [] on vector loops}
  string LIB2RAW_CheckProjectFromDirectory(const string& directory) {
    const bool LDEBUG = XHOST.DEBUG;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << endl;
    }
    CheckMaterialServer(__AFLOW_FUNC__); // must be in AFLOW_MATERIALS_SERVER
    // find from PWD
    string PROJECT_LIBRARY = "NOTHING";
    string directory_pwd = aurostd::getPWD();
    CleanDirectoryLIB(directory_pwd);

    // if not found by PWD switch to directory
    for (size_t i = 0; i < vAFLOW_PROJECTS_DIRECTORIES.size(); i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " looking for " << vAFLOW_PROJECTS_DIRECTORIES[i] << " in " << directory << endl;
      }  // CO20200624
      if (aurostd::substring2bool(directory, vAFLOW_PROJECTS_DIRECTORIES[i])) {
        PROJECT_LIBRARY = vAFLOW_PROJECTS_DIRECTORIES[i];
      }
    }
    if (PROJECT_LIBRARY != "NOTHING") {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FOUND from directory: " << PROJECT_LIBRARY << endl;
      }
      return PROJECT_LIBRARY;
    }

    for (size_t i = 0; i < vAFLOW_PROJECTS_DIRECTORIES.size(); i++) {
      if (LDEBUG) {
        cerr << XPID << vAFLOW_PROJECTS_DIRECTORIES[i] << endl;
      }
      if (aurostd::substring2bool(directory_pwd, vAFLOW_PROJECTS_DIRECTORIES[i])) {
        PROJECT_LIBRARY = vAFLOW_PROJECTS_DIRECTORIES[i];
      }
    }
    if (PROJECT_LIBRARY != "NOTHING") {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FOUND from pwd: " << PROJECT_LIBRARY << endl;
      }
      return PROJECT_LIBRARY;
    }

    if (PROJECT_LIBRARY == "NOTHING") {
      const string message = "Nothing found from pwd or directory [directory=" + directory + "]   [pwd=" + directory_pwd + "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return XHOST.tmpfs;
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::XFIX_LIBRARY_ALL
// ***************************************************************************
namespace aflowlib {
  /// @brief copies given system in the library to a fix folder, lib to fix, raw to tmpfs, removes raw
  /// @param LIBRARY_IN library to apply fix for
  /// @param argv system/structure to apply fix for
  /// @authors
  /// @mod{ST,20241022,created doxy\, use aurostd filesys funcs}
  void XFIX_LIBRARY_ALL(const string& LIBRARY_IN, vector<string> argv) {
    CheckMaterialServer(__AFLOW_FUNC__); // must be in AFLOW_MATERIALS_SERVER
    const string& system = argv.at(2);
    const string& structure = argv.at(3);  // if 3 inputs need to fix with fewer inputs
    const string systemstructure = system + "/" + structure;
    const string lib_fix = string(LIBRARY_IN).append("/FIX/").append(system);
    const string lib_lib = string(LIBRARY_IN).append("/LIB/").append(systemstructure);
    const string lib_raw = string(LIBRARY_IN).append("/RAW/").append(systemstructure);
    cout << "Fixing " << lib_lib << endl;
    if (!aurostd::FileExist(lib_fix)) {
      aurostd::DirectoryMake(lib_fix);
    }
    aurostd::file2file(lib_lib, lib_fix);
    aurostd::file2file(lib_raw, string(XHOST.tmpfs.append("/").append(system)));
    aurostd::RemoveDirectory(lib_raw);
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_ALL
// ***************************************************************************
namespace aflowlib {
  /// @brief run LIB2RAW for everything with an _AFLOWIN_
  /// @param options options in the form of dir or all[,dir]
  /// @param flag_FORCE whether to force if target directory already exists
  /// @return true on succesful completion, otherwise false
  /// @authors
  /// @mod{ST,20241022,created doxy\, optimize\, cleanup}
  bool LIB2RAW_ALL(const string& options, bool flag_FORCE) {
    const bool LDEBUG = XHOST.DEBUG;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " options=" << options << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ","); // QUICK FIX

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
    }
    if (tokens.size() > 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, aurostd::liststring2string("aflow --lib2raw=directory", "aflow --lib2raw=all[,dir]"));
    }
    if (!tokens.empty()) {
      if (tokens[0] != "all") {
        init::ErrorOption(options, __AFLOW_FUNC__, aurostd::liststring2string("aflow --lib2raw=directory", "aflow --lib2raw=all[,dir]"));
      }
    }

    CheckMaterialServer(__AFLOW_FUNC__); // must be in AFLOW_MATERIALS_SERVER
    string PROJECT_LIBRARY;
    if (tokens.size() == 2) {
      PROJECT_LIBRARY = aflowlib::LIB2RAW_CheckProjectFromDirectory(tokens[1]);
    } else {
      PROJECT_LIBRARY = aflowlib::LIB2RAW_CheckProjectFromDirectory(aurostd::getPWD());
    }
    cerr << __AFLOW_FUNC__ << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;

    int multi_sh_value = XHOST.CPU_Cores;
    // ME20181109 - Handle NCPUS=MAX
    if (XHOST.vflag_control.flag("NUM_THREADS") && !(XHOST.vflag_control.flag("NUM_THREADS_MAX"))) {
      multi_sh_value = aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("NUM_THREADS"));
    }

    cerr << XPID << "multi_sh_value=" << multi_sh_value << endl;
    deque<string> dcmds;

    if (PROJECT_LIBRARY == init::AFLOW_Projects_Directories("ICSD") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB0") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB1") ||
        PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB2") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB3") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB4") ||
        PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB5") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB6") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB7") ||
        PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB8") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB9")) {
      vector<string> file_list;
      for (const auto& dir_entry : std::filesystem::recursive_directory_iterator(PROJECT_LIBRARY + "/LIB/")) {
        if (aurostd::substring2bool(dir_entry.path().filename().string(), _AFLOWIN_)) {
          file_list.emplace_back(dir_entry.path().string());
        }
      }
      std::sort(file_list.begin(), file_list.end());
      vector<string> tokens = file_list;
      uint tokens_i_size = 0;
      for (size_t i = 0; i < tokens.size(); i++) {
        aurostd::StringSubstInPlace(tokens[i], "/" + _AFLOWIN_, "");
        aurostd::StringSubstInPlace(tokens[i], "/" + _AFLOWLOCK_, "");
        tokens_i_size = tokens[i].size();
        if (tokens_i_size >= 5 && tokens[i][tokens_i_size - 5] == '/' && tokens[i][tokens_i_size - 4] == 'c' && tokens[i][tokens_i_size - 3] == 'o' && tokens[i][tokens_i_size - 2] == 'r' &&
            tokens[i][tokens_i_size - 1] == 'e') {
          tokens[i] = tokens[i].substr(0, tokens_i_size - 5);
        } // aurostd::StringSubst(tokens[i],"/core",""); //CO20200624 - prevent /home/corey -> /homey
        string cmd = "aflow";
        if (XHOST.vflag_control.flag("BEEP")) {
          cmd += " --beep";
        }
        if (flag_FORCE) {
          cmd += " --force";
        }
        cmd += " --lib2raw=" + tokens[i];
        dcmds.push_back(cmd);
      };
      if (multi_sh_value == 0 || multi_sh_value == 1) {
        long double reference_seconds = aurostd::get_seconds();
        vector<long double> vdelta_seconds;
        for (size_t i = 0; i < dcmds.size(); i++) {
          aflowlib::LIB2RAW(tokens[i], flag_FORCE);
          const long double delta_seconds = aurostd::get_delta_seconds(reference_seconds);
          if (delta_seconds > 1.0) {
            vdelta_seconds.push_back(delta_seconds);
            const long double eta_seconds = aurostd::mean(vdelta_seconds) * (dcmds.size() - vdelta_seconds.size());
            cout << __AFLOW_FUNC__ << " [STEP]"
                 << "  DONE= " << vdelta_seconds.size() << " / " << dcmds.size() - vdelta_seconds.size() << " (" << 100 * vdelta_seconds.size() / dcmds.size() << ")"
                 << "  iSEC=" << vdelta_seconds.at(vdelta_seconds.size() - 1) << "  aSEC=" << aurostd::mean(vdelta_seconds) << "  TO_DO="
                 << dcmds.size() - vdelta_seconds.size()
              //       << "  ETA(secs)=" << eta_seconds
              //       << "  ETA(mins)=" << eta_seconds/(60.0)
              //       << "  ETA(hours)=" << eta_seconds/(60*60)
                 << "  ETA(days)="
                 << eta_seconds / (60 * 60 * 24)
              //       << "  ETA(weeks)=" << eta_seconds/(60*60*24*7)
              //       << "  ETA(years)=" << eta_seconds/(60*60*24*365)
                 << endl;
          }
        }
      } else {
        // DO MULTI
        aurostd::multithread_execute(dcmds, multi_sh_value, false);
      }
      return true;
    }
    const string message = "Project Not Found";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    return false;
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::GetSpeciesDirectory
// ***************************************************************************
namespace aflowlib {
  /// @brief get species in the direccttory
  /// @param directory[in] directory to look in
  /// @param vspecies[out] species in the directory
  /// @return number of species found
  /// @authors
  /// @mod{ST,20241022,created doxy\, optimize\, cleanup}
  uint GetSpeciesDirectory(const string& directory, vector<string>& vspecies) {
    const bool LDEBUG = XHOST.DEBUG;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vspecies.clear();
    vector<string> vs;
    vector<string> tokens;
    stringstream oss;

    if (XHOST.vext.size() != XHOST.vcat.size()) {
      const string message = "XHOST.vext.size()!=XHOST.vcat.size()";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vspecies.size()=" << vspecies.size() << " [1]" << endl;
    }
    for (const auto& iext : XHOST.vext) { // check _AFLOWIN_.EXT
      const string path = string(directory).append("/").append(_AFLOWIN_).append(iext);
      if (vspecies.empty() && aurostd::FileExist(path)) {
        vector<string> vlines;
        aurostd::compressfile2vectorstring(path, vlines);
        auto predicate = [](const string& line) { return !std::regex_search(line, std::regex("VASP_POTCAR_FILE")); };
        vlines.erase(std::remove_if(vlines.begin(), vlines.end(), predicate), vlines.end());
        vspecies = vlines;
        for (auto& vspecie : vspecies) {
          aurostd::StringSubstInPlace(vspecie, "[VASP_POTCAR_FILE]", "");
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " vspecies=" << aurostd::joinWDelimiter(vspecies, ",") << endl;  // CO20210315
        }
        for (auto& vspecie : vspecies) {
          vspecie = KBIN::VASP_PseudoPotential_CleanName(vspecie); // CO20210315
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " vspecies=" << aurostd::joinWDelimiter(vspecies, ",") << endl;  // CO20210315
        }
        break;
      }
    }
    vector<string> poscars{"POSCAR.orig", "POSCAR.relax1", "POSCAR.bands"};
    for (size_t j = 0; j < poscars.size(); j++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vspecies.size()=" << vspecies.size() << " [" << j + 2 << "]" << endl;
      }
      for (const auto& iext : XHOST.vext) {
        const string path = string(directory).append("/").append(poscars[j]).append(iext);
        if (vspecies.empty() && aurostd::FileExist(path)) {
          aurostd::compressfile2stringstream(path, oss);
          xstructure xstr(oss, IOVASP_POSCAR);
          if (!xstr.species.empty() && !xstr.species[0].empty()) {
            vspecies.assign(xstr.species.begin(), xstr.species.end());
            break;
          }
        }
      }
    }
    vector<string> outcars{"OUTCAR.relax1", "OUTCAR.static"};
    for (size_t j = 0; j < outcars.size(); j++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vspecies.size()=" << vspecies.size() << " [" << j + 2 + poscars.size() << "]" << endl;
      }
      for (const auto& iext : XHOST.vext) {
        const string path = string(directory).append("/").append(outcars[j]).append(iext);
        if (vspecies.empty() && aurostd::FileExist(path)) {
          vector<string> vlines;
          aurostd::compressfile2vectorstring(path, vlines);
          auto predicate = [](const string& line) { return !std::regex_search(line, std::regex("TITEL")); };
          vlines.erase(std::remove_if(vlines.begin(), vlines.end(), predicate), vlines.end());
          vs = vlines;
          vspecies.reserve(vs.size());
          for (const auto& v : vs) {
            aurostd::string2tokens(v, tokens, " ");
            vspecies.emplace_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(3)));
          }
          break;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return vspecies.size();
  }
} // namespace aflowlib

namespace aflowlib {
  /// @brief cleans the directory of certain paths and whitespace
  /// @param directory[in/out] to be cleaned
  /// @authors
  /// @mod{CO,20200624,created function}
  /// @mod{ST,20241022,created doxy\, optimize\, cleanup\, use loops}
  void CleanDirectoryLIB(string& directory) {
    if (directory.empty()) {
      return;
    }

    const vector<std::pair<string, string>> substitution_pairs{
        {     "/" + _AFLOWIN_,             ""},
        {   "/" + _AFLOWLOCK_,             ""},
        {  "/common/GNDSTATE", "/common/LIB2"},
        {      "common/SCINT",  "common/ICSD"},
        {             "SCINT",         "ICSD"},
        {"common/ELPASOLITES",  "common/AURO"},
        {       "ELPASOLITES",         "AURO"},
        {      "LIBRARYX/RAW", "LIBRARYX/LIB"},
        {          "LIB2/RAW",     "LIB2/LIB"},
        {             "/RAW/",        "/LIB/"}
    };

    for (const auto& [find, replace] : substitution_pairs) {
      aurostd::StringSubstInPlace(directory, find, replace);
    }

    aurostd::RemoveTrailingCharacter_InPlace(directory, '/');
    directory = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(directory);

    if (directory.empty()) {
      return;
    }
    if (directory[0] != '/') {
      directory = aurostd::getPWD() + "/" + directory;
      aurostd::StringSubstInPlace(directory, "//", "/");
    }
  }
} // namespace aflowlib

namespace aflowlib {
  /// @brief sets the aurl based on the given directory
  /// @param aflowlib_data data that gets its aurl set
  /// @param directory_LIB direcctory to use for setting aurl
  /// @param LOCAL whether to perform locally
  /// @authors
  /// @mod{CO,20220124,created function}
  /// @mod{ST,20241022,created doxy\, optimize\, cleanup\, use loops}
  void setAURL(aflowlib::_aflowlib_entry& aflowlib_data, const string& directory_LIB, bool LOCAL) {
    const bool LDEBUG = XHOST.DEBUG;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }

    if (LOCAL) {
      aflowlib_data.aurl = directory_LIB; // dummy
    } else {
      // build aflowlib_data.aurl
      aflowlib_data.aurl = AFLOWLIB_SERVER_DEFAULT + ":" + directory_LIB;

      const vector<string> substitutions{"/common", "/scratch/common", "/SCRATCH/common", "/work/common", "/WORK/common", "/archive/common", "/ARCHIVE/common", "/home/" + XHOST.user + "/common"};
      for (const auto& substitution : substitutions) {
        aurostd::StringSubstInPlace(aflowlib_data.aurl, XHOST.home + substitution, "common");
      }

      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " aurl(PRE )=" << aflowlib_data.aurl << endl;
      }

      if (aurostd::substring2bool(aflowlib_data.aurl, "LIBRARYX")) {
        aurostd::StringSubstInPlace(aflowlib_data.aurl, "common/GNDSTATE/LIBRARYX/LIB", "AFLOWDATA/LIB2_RAW");
        aflowlib_data.catalog = "LIBRARYX";
      } // [HISTORIC]
      const vector<string> aurls = {"ICSD", "LIB0", "LIB1", "LIB2", "LIB3", "LIB4", "LIB5", "LIB6", "LIB7", "LIB8", "LIB9", "AURO"};
      for (auto& aurl : aurls) {
        if (aurostd::substring2bool(aflowlib_data.aurl, aurl)) {
          aurostd::StringSubstInPlace(aflowlib_data.aurl, "common/" + aurl + "/LIB", "AFLOWDATA/" + aurl + "_RAW");
          aflowlib_data.catalog = aurl;
        }
      }
      aurostd::StringSubstInPlace(aflowlib_data.aurl, ":/AFLOWDATA", ":AFLOWDATA");

      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " aurl(POST)=" << aflowlib_data.aurl << endl;
      }
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW
// ***************************************************************************
namespace aflowlib {
  /// @brief adds new calcs to aflow database
  /// @param options options in the form of dir or all[,dir]
  /// @param flag_FORCE whether to force if target directory already exists
  /// @param LOCAL perform database operations locally only, useful for testing
  /// @authors
  /// @mod{ST,20241022,created doxy\, changed signature\, optimize\, cleanup}
  void LIB2RAW(const string& options, bool flag_FORCE, bool LOCAL) {
    const bool LDEBUG = XHOST.DEBUG;
    const long double seconds_begin = aurostd::get_seconds();
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " options=" << options << endl;
    }

    // Call to LIB2LIB function to run postprocessing - added by CT20181212
    const bool perform_LIB2LIB = true; // CO20200624 - good for debugging the rest of lib2raw, lib2lib can be slow
    if (perform_LIB2LIB && !LIB2LIB(options, flag_FORCE, LOCAL)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "aflowlib::LIB2LIB() failed", _RUNTIME_ERROR_);
    }

    vector<string> tokens;
    if (!aurostd::substring2bool(options, "POCC")) {
      aurostd::string2tokens(options, tokens, ","); // QUICK FIX
    } else {
      tokens.clear();
      tokens.push_back(options);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
    }
    if (tokens.empty() || tokens.size() > 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, aurostd::liststring2string("aflow --lib2raw=directory", "aflow --lib2raw=all[,dir]"));
    }

    if (!tokens.empty()) {
      if (tokens[0] == "all") {
        XHOST.sensors_allowed = false;
        XHOST.vflag_pflow.flag("MULTI=SH", false);
        LIB2RAW_ALL(options, flag_FORCE);
        XHOST.sensors_allowed = true;
        return;
      }
    }

    string directory = options;
    CleanDirectoryLIB(directory);

    string directory_LIB;
    string directory_RAW;
    string directory_WEB;
    bool flag_WEB = false;
    const bool flag_files_LIB = false;
    bool flag_files_RAW = false;
    bool flag_files_WEB = false;
    string PROJECT_LIBRARY;

    if (LOCAL) {
      flag_FORCE = true;
      directory_LIB = directory;
      directory_RAW = directory + "/RAW";
      directory_WEB = directory + "/WEB";
      PROJECT_LIBRARY = directory_LIB;
    } else {
      CheckMaterialServer(__AFLOW_FUNC__);
      PROJECT_LIBRARY = LIB2RAW_CheckProjectFromDirectory(directory);
      cout << __AFLOW_FUNC__ << " directory=" << directory << endl;

      //  bool flag_LIB=false;
      const vector<std::tuple<string, uint, bool>> vmap_libs{
        // lib  lib uint  web flag
          {"LIB0", XHOST_LIBRARY_LIB0, false},
          {"LIB1", XHOST_LIBRARY_LIB1, false},
          {"LIB2", XHOST_LIBRARY_LIB2, false},
          {"LIB3", XHOST_LIBRARY_LIB3,  true},
          {"LIB4", XHOST_LIBRARY_LIB4,  true},
          {"LIB5", XHOST_LIBRARY_LIB5,  true},
          {"LIB6", XHOST_LIBRARY_LIB6,  true},
          {"LIB7", XHOST_LIBRARY_LIB7,  true},
          {"LIB8", XHOST_LIBRARY_LIB8,  true},
          {"LIB9", XHOST_LIBRARY_LIB9,  true}
      };
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " scan libraries BEGIN [1]" << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " PROJECT_LIBRARY=" << PROJECT_LIBRARY << endl;
      }
      for (const auto& [lib, libint, web_flag] : vmap_libs) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " XHOST_LIBRARY_" << lib << "=" << libint << endl;
        }
        if (PROJECT_LIBRARY == init::AFLOW_Projects_Directories(lib)) {
          flag_WEB = web_flag;
          flag_files_RAW = true;
        }
      }
      // different for ICSD
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " XHOST_LIBRARY_ICSD=" << XHOST_LIBRARY_ICSD << endl;
      }
      if (PROJECT_LIBRARY == init::AFLOW_Projects_Directories("ICSD")) {
        flag_WEB = true;
        flag_files_WEB = true;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " scan libraries END [2]" << endl;
      }

      if (PROJECT_LIBRARY == init::AFLOW_Projects_Directories("ICSD")) {
        vector<string> path_parts;
        aurostd::string2tokens(directory, path_parts, "/");
        auto predicate_icsd = [](const string& part) { return aurostd::substring2bool(part, "_ICSD_"); };
        const bool find_icsd = std::any_of(path_parts.begin() + 1, path_parts.end(), predicate_icsd);
        if (!find_icsd) {
          cerr << __AFLOW_FUNC__ << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
          cerr << __AFLOW_FUNC__ << " you must specify the directory including the whole lattice type" << endl;
          cerr << " such as  aflow --lib2raw=FCC/La1Se1_ICSD_27104  " << endl;
          return;
        }
      }

      string dir_base = directory;
      if (aurostd::substring2bool(dir_base, "LIB/")) {
        dir_base = aurostd::substring2string(dir_base, "LIB/", 1, false);
      }
      directory_LIB = aurostd::CleanFileName(PROJECT_LIBRARY + "/LIB/" + dir_base);
      directory_RAW = aurostd::CleanFileName(PROJECT_LIBRARY + "/RAW/" + dir_base);
      directory_WEB = aurostd::CleanFileName(PROJECT_LIBRARY + "/WEB/" + dir_base);

      directory_LIB.erase(directory_LIB.find_last_not_of('/') + 1, string::npos);

      aurostd::StringSubstInPlace(directory_RAW, "RAW/LIB", "RAW");
      aurostd::StringSubstInPlace(directory_WEB, "WEB/LIB", "WEB");

      // CO20200624 - BELOW HERE DIRECTORY_LIB IS FIXED!!!!!!!!!

      if (flag_WEB == false) {
        directory_WEB = aurostd::CleanFileName(PROJECT_LIBRARY + "/WEB/");
      }

      if (PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB0") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB1") || PROJECT_LIBRARY == init::AFLOW_Projects_Directories("LIB2")) {
        //      directory_WEB=directory_RAW;
      }
      cout << __AFLOW_FUNC__ << " PROJECT_LIBRARY=" << PROJECT_LIBRARY << endl;
    }
    cout << __AFLOW_FUNC__ << " directory_LIB=" << directory_LIB << endl;
    cout << __AFLOW_FUNC__ << " directory_RAW=" << directory_RAW << endl;
    cout << __AFLOW_FUNC__ << " directory_WEB=" << directory_WEB << endl;

    if (!directory_LIB.empty()) {
      XHOST.vflag_control.flag("DIRECTORY_CLEAN", true);
      XHOST.vflag_control.push_attached("DIRECTORY_CLEAN", directory_LIB);
    } // CO20200624 - fix error messages to point to directory_LIB

    // NOLINTBEGIN(*-const-correctness)
    bool perform_LOCK = false;  // CO20200624 - turning off in general, check below
    bool perform_STATIC = false;
    bool perform_BANDS = false;
    bool perform_DIELECTRIC = false;
    bool perform_BADER = false;
    bool perform_THERMODYNAMICS = false;
    bool perform_AGL = false;
    bool perform_AEL = false;
    bool perform_APL = false; // ME20210901
    bool perform_QHA = false; // AS20200831
    bool perform_POCC = false;  // CO20200624
    bool perform_PATCH = false; // to inject updates while LIB2RAW  //CO20200624 - turning off in general, check below
    // NOLINTEND(*-const-correctness)

    const vector<std::pair<bool*, vector<string>>> vmap_task2files{
        {          &perform_LOCK,                                                        {_AFLOWLOCK_}},
        {        &perform_STATIC,                                                    {"OUTCAR.static"}},
        {         &perform_BANDS,                                                     {"OUTCAR.bands"}},
        {    &perform_DIELECTRIC,                                                {"OUTCAR.dielectric"}},
        {         &perform_BADER,                                 {"AECCAR0.static", "AECCAR2.static"}},
        {&perform_THERMODYNAMICS, {"OUTCAR.relax1", "OUTCAR.relax2", "OUTCAR.relax3", "OUTCAR.static"}},
        {           &perform_AGL,                                                    {"aflow.agl.out"}},
        {           &perform_AEL,                                                    {"aflow.ael.out"}},
        {           &perform_APL,                     {DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_OUT_FILE}},
        {           &perform_QHA,                                    {DEFAULT_QHA_FILE_PREFIX + "out"}},
        {          &perform_POCC,                     {POCC_FILE_PREFIX + POCC_UNIQUE_SUPERCELLS_FILE}}
    }; // ST20240912 - assign to bools based on existence of any files in respective lists
    // CO20200624 - do NOT use CompressFileExist(), we want to EXPLICITLY catch those that are zipped
    for (const auto& [perform_task, files] : vmap_task2files) {
      for (const auto& file : files) {
        for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
          const string path = string(directory_LIB).append("/").append(file).append(XHOST.vext[iext]);
          *perform_task |= aurostd::FileExist(path);
        }
      }
    }
    perform_PATCH = perform_THERMODYNAMICS || perform_BANDS || perform_STATIC;

    if (!aurostd::FileExist(directory_LIB + "/" + _AFLOWIN_)) {
      cout << __AFLOW_FUNC__ << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
      cout << __AFLOW_FUNC__ << " directory does not exist: " << directory_LIB << endl;
      return;
    }
    if (flag_FORCE == false) { // only goes in this loop if no --force, this happens BEFORE directory_RAW is deleted
      // CO20200624 - checking files in RAW from previous lib2raw run
      if (aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {  // directory_RAW+"/"+_AFLOWIN_)
        const _aflags aflags;
        aurostd::compressfiles2compressfiles(directory_LIB, aurostd::compression_type::XZ);

        cout << __AFLOW_FUNC__ << " ALREADY CALCULATED = " << directory_RAW << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << Message(__AFLOW_FILE__, aflags, "time") << endl;
        return;
      }
      if (perform_BANDS) {
        if (aurostd::CompressFileExist(directory_RAW + "/OUTCAR.bands")) { // CO20200624 - EIGENVAL.bands->OUTCAR.bands
          cout << __AFLOW_FUNC__ << " directory is skipped because of BANDS: " << directory_RAW << endl;
          return;
        }
      }
    }
    // directory_LIB exists and directory_RAW does not exist can move on:
    if (!aurostd::DirectoryMake(directory_RAW)) {
      cout << __AFLOW_FUNC__ << " directory is skipped because directory_RAW cannot be created: " << directory_RAW << endl;
      return;
    }
    aurostd::RemoveFiles(directory_RAW);
    if (flag_WEB) {
      if (!aurostd::DirectoryMake(directory_WEB)) {
        cout << __AFLOW_FUNC__ << " directory is skipped because directory_WEB cannot be created: " << directory_WEB << endl;
        return;
      }
      aurostd::RemoveFiles(directory_RAW);
    }

    cout << __AFLOW_FUNC__ << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;

    const _aflags aflags;
    cout << __AFLOW_FUNC__ << " dir=" << directory_LIB << " BEGIN_DATE = " << Message(__AFLOW_FILE__, aflags) << endl;
    aurostd::compressfiles2compressfiles(directory_LIB, aurostd::compression_type::XZ);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    // DX20190516 [ADD THIS LINE ONCE TESTED WITH REAL-WORLD CASE] if(aurostd::RemoveControlCodeCharactersFromFile(directory_LIB,_AFLOWIN_)) { //DX20190516 - Ensure no control characters in aflow.in; only
    // modifies if control characters are detected DX20190516 [ADD THIS LINE ONCE TESTED WITH REAL-WORLD CASE] }

    vector<string> vfile;   // the needed files
    aflowlib::_aflowlib_entry aflowlib_data;

    // CO20200731 - species may get overwritten later
    GetSpeciesDirectory(directory_LIB, aflowlib_data.vspecies);
    aflowlib_data.species = aurostd::joinWDelimiter(aflowlib_data.vspecies, ",");
    aflowlib_data.nspecies = aflowlib_data.vspecies.size();

    cout << __AFLOW_FUNC__ << " nspecies=" << aflowlib_data.nspecies << endl;
    cout << __AFLOW_FUNC__ << " species=" << aflowlib_data.species << endl;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [2]" << endl;
    }

    aflowlib_data.system_name = KBIN::ExtractSystemName(directory_LIB); // ME20200207 - the system name is the canonical title

    if (perform_THERMODYNAMICS || perform_POCC) {
      aurostd::string2tokens(directory_LIB, tokens, "/");
      aflowlib_data.prototype = tokens.at(tokens.size() - 1);
      aurostd::StringSubstInPlace(aflowlib_data.prototype, ":LDAU2", "");
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [3]" << endl;
    }

    aflowlib::setAURL(aflowlib_data, directory_LIB, LOCAL); // build aflowlib_data.aurl

    // build aflowlib_data.auid

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [AUID=0] directory_LIB=" << directory_LIB << endl;
    }
    if (perform_POCC) { // CO20200624
      const string metadata_auid_json = aflowlib_data.POCCdirectory2MetadataAUIDjsonfile(directory_LIB);
      aurostd::string2file(metadata_auid_json, aurostd::CleanFileName(directory_RAW + "/" + "metadata_auid.json"));
    } else {
      aflowlib_data.directory2auid(directory_LIB); // NEW STYLE
      if (LOCAL) {
        aurostd::StringSubstInPlace(aflowlib_data.auid, "aflow", "local");
      } // SD20240112 - new AUID style for local lib2raw
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [AUID=0b] directory_LIB=" << directory_LIB << endl;
    }
    cout << __AFLOW_FUNC__ << " AURL  = " << aurostd::PaddedPOST(aflowlib_data.aurl, 60) << endl;
    cout << __AFLOW_FUNC__ << " AUID  = " << aurostd::PaddedPOST(aflowlib_data.auid, 60) << endl;
    cout << __AFLOW_FUNC__ << " VAUID = " << aflowlib::auid2directory(aflowlib_data.auid) << endl;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [AUID=2]" << endl;
    }
    cout << __AFLOW_FUNC__ << " CATALOG = " << aurostd::PaddedPOST(aflowlib_data.catalog, 60) << endl;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [AUID=3]" << endl;
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the THERMODYNAMICS
    if (perform_THERMODYNAMICS) {
      cout << __AFLOW_FUNC__ << " THERMODYNAMIC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Thermodynamics(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (thermodynamics):", LOCAL); // identifier inside
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [4]" << endl;
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the STATIC
    if (perform_STATIC) {
      cout << __AFLOW_FUNC__ << " STATIC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Static(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (static):");
      // MOVE/LINK PICS data
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the BANDS
    if (perform_BANDS) {
      cout << __AFLOW_FUNC__ << " BANDS LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Bands(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (bands):");
      // MOVE/LINK PICS data
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the STATIC
    if ((perform_STATIC || perform_BANDS)) { // CO20200731 - MAGNETIC->STATIC //JX
      cout << __AFLOW_FUNC__ << " MAGNETIC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Magnetic(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (magnetic):");
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the DIELECTRIC
    if (perform_DIELECTRIC) {
      cout << __AFLOW_FUNC__ << " DIELECTRIC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Dielectric(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (dielectric):");
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the BADER
    //
    if (perform_BADER) {
      cout << __AFLOW_FUNC__ << " BADER LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Bader(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (bader):");
      // MOVE/LINK PICS data
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the AGL
    if (perform_AGL) {
      cout << __AFLOW_FUNC__ << " AGL LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_AGL(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (agl):");
      if (flag_WEB) {
        const vector<string> agl_file_list{"aflow.agl.out",
                                           "AGL.out",
                                           "AGL_energies_temperature.out",
                                           "AGL_thermal_properties_temperature.out",
                                           "AGL_edos_gap_pressure.out",
                                           "AGL_edos_gap_pressure.json",
                                           "AGL_energy.json",
                                           "AGL_energy_structures.json",
                                           "AGL_energy_volume.out",
                                           "AGL_gibbs_energy_pT.out",
                                           "AGL_Hugoniot.out"};
        for (const auto& file : agl_file_list) {
          const string path = string(directory_RAW).append("/").append(file);
          if (aurostd::FileExist(path)) {  // LINK  //CO20200624 - adding FileExist() check
            aurostd::LinkFile(path, directory_WEB);
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the AEL
    if (perform_AEL) {
      cout << __AFLOW_FUNC__ << " AEL LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_AEL(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (ael):");
      if (flag_WEB) {
        const vector<string> ael_file_list{"aflow.ael.out", "AEL_Elastic_constants.out", "AEL_Compliance_tensor.out", "AEL_elastic_tensor.json", "AEL_energy_structures.json"};
        for (const auto& file : ael_file_list) {
          const string path = string(directory_RAW).append("/").append(file);
          if (aurostd::FileExist(path)) {  // LINK  //CO20200624 - adding FileExist() check
            aurostd::LinkFile(path, directory_WEB);
          }
        }
      }
    }
    // BEGIN ME20210901
    // do the APL loop
    if (perform_APL) {
      cout << __AFLOW_FUNC__ << " APL LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_APL(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (apl):");
      if (flag_WEB) {
        const vector<string> apl_file_list{DEFAULT_APL_PHDOSCAR_FILE,
                                           DEFAULT_APL_PHKPOINTS_FILE,
                                           DEFAULT_APL_PHEIGENVAL_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_JSON,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_MSQRDISP_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_AAPL_GVEL_FILE,
                                           aflowlib_data.system_name + "_phdosdata.json"};
        for (const auto& file : apl_file_list) {
          const string path = string(directory_RAW).append("/").append(file);
          if (aurostd::FileExist(path)) {
            aurostd::LinkFile(path, directory_WEB);
          }
        }
        vector<string> files;
        aurostd::DirectoryLS(directory_LIB, files);
        for (size_t i = 0; i < files.size(); i++) {
          if ((files[i].find("phdispdos.png") != string::npos) || (files[i].find("phdisp.png") != string::npos) || (files[i].find("phdos.png") != string::npos)) {
            aurostd::LinkFile(directory_RAW + "/" + files[i], directory_WEB);
          }
        }
      }
    }
    // END ME20210901
    // BEGIN AS20200831
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the QHA
    if (perform_QHA) {
      cout << __AFLOW_FUNC__ << " QHA LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_QHA(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (qha):");
      if (flag_WEB) {
        const vector<string> qha_file_list{DEFAULT_QHA_FILE_PREFIX + "out", DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_THERMO_FILE, DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_FVT_FILE,
                                           DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_PDIS_FILE + ".T300K.out", DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_PDIS_FILE + ".T300K.json"};
        for (const auto& file : qha_file_list) {
          const string path = string(directory_RAW).append("/").append(file);
          if (aurostd::FileExist(path)) {
            aurostd::LinkFile(path, directory_WEB);
          }
        }
        // link all QHA plots
        vector<string> files;
        aurostd::DirectoryLS(directory_LIB, files);
        for (size_t i = 0; i < files.size(); i++) {
          if (files[i].find("qha") != string::npos && files[i].find(".png") != string::npos) {
            aurostd::LinkFile(directory_RAW + "/" + files[i], directory_WEB);
          }
        }
      }
    }
    // END AS20200831
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the POCC
    if (perform_POCC) {
      cout << __AFLOW_FUNC__ << " POCC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_POCC(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (POCC):");
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the LOCK
    if (perform_LOCK) {
      cout << __AFLOW_FUNC__ << " LOCK LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_LOCK(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (LOCK):");
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the PATCH
    if (perform_PATCH) {
      cout << __AFLOW_FUNC__ << " PATCH LOOP --------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_PATCH(directory_LIB, directory_RAW, vfile, aflowlib_data, __AFLOW_FUNC__ + " (PATCH):");
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // write DOS + BANDS JSON //CO20171025
    if ((aurostd::FileExist(directory_LIB + "/DOSCAR.static") || aurostd::CompressFileExist(directory_LIB + "/DOSCAR.static")) &&
        (aurostd::FileExist(directory_LIB + "/POSCAR.static") || aurostd::CompressFileExist(directory_LIB + "/POSCAR.static"))) {
      cout << __AFLOW_FUNC__ << " DOS + BANDS JSON  ---------------------------------------------------------------------------------" << endl;
      stringstream _json_file;
      stringstream json_file;
      aurostd::xoption vpflow;  // dummy
      if (estructure::DOSDATA_JSON(vpflow, directory_LIB, _json_file, true)) {
        json_file << _json_file.str() << endl;
        _json_file.str("");
        cout << __AFLOW_FUNC__ << " compressing file: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata.json") << endl;
        cout.flush();
        aurostd::stringstream2compressfile(json_file, directory_RAW + "/" + aflowlib_data.system_name + "_dosdata.json");
        json_file.str("");
        if (aurostd::CompressFileExist(directory_LIB + "/EIGENVAL.bands") && aurostd::CompressFileExist(directory_LIB + "/KPOINTS.bands")) {
          if (estructure::BANDSDATA_JSON(vpflow, directory_LIB, _json_file, true)) {
            json_file << _json_file.str() << endl;
            _json_file.str("");
            cout << __AFLOW_FUNC__ << " compressing file: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata.json") << endl;
            cout.flush();
            aurostd::stringstream2compressfile(json_file, directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata.json");
            json_file.str("");
          }
        }
      }
    }

    // ---------------------------------------------------------------------------------------------------------------------------------
    // DO THE COMPRESSING OF VASP FILES
    //      cout << "COMPRESSING" << endl;
    cout << __AFLOW_FUNC__ << " COMPRESSING REMOVING LINKING --------------------------------------------------------------" << endl;

    // generic for compressing and linking
    deque<string> vtype = {".orig", ".relax", ".relax1", ".relax2", ".relax3", ".relax4", ".static", ".bands", ".dielectric"};
    deque<string> vout = {".out", ".json"};

    // CO20200624 - removed COMPRESSED eps
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) {
      aurostd::RemoveFile(directory_RAW, std::regex(".*.eps" + XHOST.vext[iext]));
    }

    // FILES to remove
    deque<string> vfile2remove = {"KPOINTS.bands.old", "EIGENVAL.bands.old", "OUTCAR.relax1"}; // CO20200404 - removed OUTCAR.bands from this list, we need for Egap. EIGENVAL.bands does not have occupancies
    vfile2remove.emplace_back("aflow.pgroupk_xtal.out"); // comes from nowhere DX - we have variants for each VASP run (relax, static, bands), so this is NOW redundant
    for (size_t iremove = 0; iremove < vfile2remove.size(); iremove++) {
      if (aurostd::FileExist(directory_RAW + "/" + vfile2remove[iremove])) { // need to be present
        cout << __AFLOW_FUNC__ << " removing file: " << string(directory_RAW + "/" + vfile2remove[iremove]) << endl;
        aurostd::RemoveFile(directory_RAW + "/" + vfile2remove[iremove]);
      } // FileExist
    } // iremove
    // FILES to compress if not compressed already (linked), in such case they will be deleted.
    deque<string> vfile2compress0 = {"OUTCAR.relax"};
    for (size_t icompress = 0; icompress < vfile2compress0.size(); icompress++) {
      if (aurostd::FileExist(directory_RAW + "/" + vfile2compress0[icompress] + "." + DEFAULT_KZIP_BIN)) { // test if there is one already compressed (possibly linked)
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " found compressed file: " << string(directory_RAW + "/" + vfile2compress0[icompress] + "." + DEFAULT_KZIP_BIN) << endl;
          cout.flush();
        }
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " removing file: " << string(directory_RAW + "/" + vfile2compress0[icompress]) << endl;
          cout.flush();
        }
        aurostd::RemoveFile(directory_RAW + "/" + vfile2compress0[icompress]);
      }
      if (aurostd::FileExist(directory_RAW + "/" + vfile2compress0[icompress])) { // need to be present
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " compressing file: " << string(directory_RAW + "/" + vfile2compress0[icompress]) << endl;
          cout.flush();
        }
        aurostd::CompressFile(directory_RAW + "/" + vfile2compress0[icompress]);
      } // FILES to compress
    } // icompress
    // FILES to compress
    deque<string> vfile2compress1 = {"aflow.pgroup", "aflow.pgroup_xtal", "aflow.pgroupk", "aflow.pgroupk_xtal", "aflow.pgroupk_Patterson",
                                     "aflow.fgroup", "aflow.iatoms",      "aflow.agroup"}; // DX20200520 - added Patterson
    for (size_t ilink = 0; ilink < vfile2compress1.size(); ilink++) {
      for (size_t iout = 0; iout < vout.size(); iout++) {
        for (size_t itype = 0; itype < vtype.size(); itype++) {
          if (aurostd::FileExist(directory_RAW + "/" + vfile2compress1[ilink] + vtype[itype] + vout[iout])) {
            if (LDEBUG) {
              cout << __AFLOW_FUNC__ << " compressing file (" << ilink << "," << iout << "," << itype << "): " << string(directory_RAW + "/" + vfile2compress1[ilink] + vtype[itype] + vout[iout]) << endl;
              cout.flush();
            }
            aurostd::CompressFile(directory_RAW + "/" + vfile2compress1[ilink] + vtype[itype] + vout[iout]);
          } // FileExist
        } // itype
      } // iout
    } // ilink

    // FILES to link LIB to RAW
    for (size_t ifile = 0; ifile < vfile.size(); ifile++) {
      if (aurostd::FileExist(directory_RAW + "/" + vfile[ifile])) { // need to be present
        deque<string> vfile2link0 = {"DOSCAR.static", "OUTCAR.static", "CHGCAR.static", "AECCAR0.static", "AECCAR2.static", "EIGENVAL.bands", "OSZICAR.bands", "OUTCAR.delectric"};
        for (size_t ilink = 0; ilink < vfile2link0.size(); ilink++) {
          if (vfile[ifile] == vfile2link0[ilink]) {
            if (aurostd::FileExist(directory_LIB + "/" + vfile[ifile]) || aurostd::CompressFileExist(directory_LIB + "/" + vfile[ifile])) { // need to be present in LIB also
              aurostd::RemoveFile(directory_RAW + "/" + vfile[ifile]); // remove RAW original
              if (aurostd::FileExist(directory_LIB + "/" + vfile[ifile])) {
                cout << __AFLOW_FUNC__ << " linking file RAW->LIB: " << string(directory_LIB + "/" + vfile[ifile]) << endl;
                cout.flush();
                aurostd::LinkFile(directory_LIB + "/" + vfile[ifile], directory_RAW); // link LIB to RAW (save space)
              }
              for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
                if (aurostd::FileExist(directory_LIB + "/" + vfile[ifile] + XHOST.vext[iext])) {
                  cout << __AFLOW_FUNC__ << " linking file RAW->LIB: " << string(directory_LIB + "/" + vfile[ifile] + XHOST.vext[iext]) << endl;
                  cout.flush();
                  aurostd::LinkFile(directory_LIB + "/" + vfile[ifile] + XHOST.vext[iext], directory_RAW); // link LIB to RAW (save space)
                }
              }
            }
          }
        } // files to link LIB to RAW
      } // File Exist
    } // ifile

    // DO THE FINISH LINK/COPY FOR WEB
    cout << __AFLOW_FUNC__ << " linking stuff flag_WEB=" << flag_WEB << endl;
    // MOVE/LINK PICS data

    if (flag_WEB) {
      cout << __AFLOW_FUNC__ << " linking SYSTEM=" << aflowlib_data.system_name << endl;
      if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + ".png") || aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_banddos.png") || // ME20190621 - new file name convention
          false) {
        aurostd::LinkFile(directory_RAW + "/*png", directory_WEB); // LINK
      }
      if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + ".cif")) {
        aurostd::LinkFile(directory_RAW + "/*cif", directory_WEB); // LINK //CO20200624 - adding FileExist() check
      }

      if (aurostd::FileExist(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php")) {
        aurostd::LinkFile(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php", directory_RAW + "/index.php"); // LINK //CO20200624 - adding FileExist() check
      }
      if (aurostd::FileExist(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php")) {
        aurostd::LinkFile(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php", directory_WEB + "/index.php"); // LINK //CO20200624 - adding FileExist() check
      }
      // CO20200624 - these are PRE-LINKS: these files have NOT been written yet
      // if we check if they exist first, then the links will not be created
      // we want to do this NOW so that it can be captured by vfiles
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT))
      aurostd::LinkFile(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT, directory_WEB); // LINK  //CO20200624 - adding FileExist() check
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON))
      aurostd::LinkFile(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_JSON, directory_WEB); // LINK //CO20200624 - adding FileExist() check
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+aflowlib_data.auid+".out"))
      aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.auid + ".out", directory_WEB); // LINK for AUID  //CO20200624 - adding FileExist() check
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+aflowlib_data.auid+".json"))
      aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.auid + ".json", directory_WEB); // LINK for AUID  //CO20200624 - adding FileExist() check

      deque<string> vfile2link1 = {"aflow.pgroup", "aflow.pgroup_xtal", "aflow.pgroupk", "aflow.pgroupk_xtal", "aflow.fgroup", "aflow.iatoms", "aflow.agroup", "edata", "data"};
      for (size_t ilink = 0; ilink < vfile2link1.size(); ilink++) {
        for (size_t iout = 0; iout < vout.size(); iout++) {
          for (size_t itype = 0; itype < vtype.size(); itype++) {
            if (aurostd::FileExist(directory_RAW + "/" + vfile2link1[ilink] + vtype[itype] + vout[iout])) { // no compression
              if (LDEBUG) {
                cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + vfile2link1[ilink] + vtype[itype] + vout[iout]) << endl;
                cout.flush();
              }
              aurostd::LinkFile(directory_RAW + "/" + vfile2link1[ilink] + vtype[itype] + vout[iout], directory_WEB);
            } // FileExist
            for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
              if (aurostd::FileExist(directory_RAW + "/" + vfile2link1[ilink] + vtype[itype] + vout[iout] + XHOST.vext[iext])) { // with compression
                if (LDEBUG) {
                  cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + vfile2link1[ilink] + vtype[itype] + vout[iout] + XHOST.vext[iext]) << endl;
                  cout.flush();
                }
                aurostd::LinkFile(directory_RAW + "/" + vfile2link1[ilink] + vtype[itype] + vout[iout] + XHOST.vext[iext], directory_WEB);
              } // FileExist
            } // iext
          } // itype
        } // iout
      } // ilink
      for (size_t iout = 0; iout < vout.size(); iout++) {
        if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_structure_relax" + vout[iout])) {
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_structure_relax" + vout[iout]) << endl;
            cout.flush();
          }
          aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_structure_relax" + vout[iout], directory_WEB); // CO20171024
        } // FileExist
      } // iout
      for (size_t iout = 0; iout < vout.size(); iout++) {
        if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_structure_relax1" + vout[iout])) {
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_structure_relax1" + vout[iout]) << endl;
            cout.flush();
          }
          aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_structure_relax1" + vout[iout], directory_WEB); // CO20171024
        } // FileExist
      } // iout
      for (size_t iout = 0; iout < vout.size(); iout++) {
        if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata" + vout[iout])) { // NO EXTENSION
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata" + vout[iout]) << endl;
            cout.flush();
          }
          aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata" + vout[iout], directory_WEB); // CO20171024
        }
        for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
          if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata" + vout[iout] + XHOST.vext[iext])) {
            if (LDEBUG) {
              cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata" + vout[iout] + XHOST.vext[iext]) << endl;
              cout.flush();
            }
            aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_dosdata" + vout[iout] + XHOST.vext[iext], directory_WEB); // CO20171024
          } // FileExist
        } // iext
      } // iout
      for (size_t iout = 0; iout < vout.size(); iout++) {
        if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata" + vout[iout])) { // NO EXTENSION
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata" + vout[iout]) << endl;
            cout.flush();
          }
          aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata" + vout[iout], directory_WEB); // CO20171024
        } // FileExist
        for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
          if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata" + vout[iout] + XHOST.vext[iext])) {
            if (LDEBUG) {
              cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata" + vout[iout] + XHOST.vext[iext]) << endl;
              cout.flush();
            }
            aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_bandsdata" + vout[iout] + XHOST.vext[iext], directory_WEB); // CO20171024
          } // FileExist
        } // iext
      } // iout
      deque<string> vfile2link2 = {"EIGENVAL.bands",   "DOSCAR.static",        "OUTCAR.static",      "CONTCAR.relax",      "CONTCAR.relax1",   "POSCAR.bands",     "CONTCAR.relax.vasp",
                                   "CONTCAR.relax.qe", "CONTCAR.relax.abinit", "CONTCAR.relax.aims", "KPOINTS.relax",      "KPOINTS.static",   "KPOINTS.bands",    "INCAR.bands",
                                   "AECCAR0.static",   "AECCAR2.static",       "CHGCAR.static",      "KPOINTS.dielectric", "INCAR.dielectric", "OUTCAR.dielectric"};
      for (size_t ilink = 0; ilink < vfile2link2.size(); ilink++) {
        if (aurostd::FileExist(directory_RAW + "/" + vfile2link2[ilink])) { // NO EXTENSION
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + vfile2link2[ilink]) << endl;
            cout.flush();
          }
          aurostd::LinkFile(directory_RAW + "/" + vfile2link2[ilink], directory_WEB); // LINK
        }
        for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
          if (aurostd::FileExist(directory_RAW + "/" + vfile2link2[ilink] + XHOST.vext[iext])) {
            if (LDEBUG) {
              cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + vfile2link2[ilink] + XHOST.vext[iext]) << endl;
              cout.flush();
            }
            aurostd::LinkFile(directory_RAW + "/" + vfile2link2[ilink] + XHOST.vext[iext], directory_WEB); // LINK
          }
        } // iext
      } // ilink

      if (perform_BADER) {
        // CO20200624 - we don't need to worry about zipped variants, this analysis is done in RAW only (no LIB compression)
        // H15Th4_ICSD_638495_Bader_50_H.jvxl - 50 will always exist
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + "*jvxl") << endl;
          cout.flush();
        }
        if (!aflowlib_data.vspecies.empty() && aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_Bader_50_" + aflowlib_data.vspecies[0] + ".jvxl")) { // CO20200624 - adding FileExist() check
          aurostd::LinkFile(directory_RAW + "/" + "*jvxl*", directory_WEB); // LINK
        }
        // H15Th4_ICSD_638495_abader.out
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " linking file WEB->RAW: " << string(directory_RAW + "/" + aflowlib_data.system_name + "_abader.out") << endl;
          cout.flush();
        } // CO20200624
        if (aurostd::FileExist(directory_RAW + "/" + aflowlib_data.system_name + "_abader.out")) {
          aurostd::LinkFile(directory_RAW + "/" + aflowlib_data.system_name + "_abader.out", directory_WEB); // LINK //CO20200624  //CO20200624 - adding FileExist() check
        }
      }
    } // flag_WEB

    // DONE
    // write files if necessary
    vector<string> vdirectory;
    // do the directories
    if (flag_files_LIB) {
      aurostd::DirectoryLS(directory_LIB, vdirectory);
      for (size_t i = 0; i < vdirectory.size(); i++) {
        aflowlib_data.vfiles.emplace_back(vdirectory[i]);
      }
    }
    if (flag_files_RAW) {
      aurostd::DirectoryLS(directory_RAW, vdirectory);
      for (size_t i = 0; i < vdirectory.size(); i++) {
        aflowlib_data.vfiles.emplace_back(vdirectory[i]);
      }
    }
    if (flag_files_WEB) {
      aurostd::DirectoryLS(directory_WEB, vdirectory);
      for (size_t i = 0; i < vdirectory.size(); i++) {
        aflowlib_data.vfiles.emplace_back(vdirectory[i]);
      }
    }
    // DO THE FINAL WRITING

    if (aflowlib_data.vaflowlib_date.size() != 2) { // CO20200624 - this means we didn't get the LOCK dates, spit warning
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "LOCK dates NOT found", _LOGGER_WARNING_);
    }
    aflowlib_data.vaflowlib_date.emplace_back(aurostd::get_datetime(true)); // CO20200624 - adding LOCK date  //CO20210624 - adding UTC offset

    if (!LOCAL) { // CO20171025
      if (aurostd::FileExist(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php")) {
        aurostd::LinkFile(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php", directory_RAW + "/" + _XENTRY_); // CO20200624 - adding FileExist() check
      }
    }

    // write aflowlib.out
    cout << __AFLOW_FUNC__ << " writing " << string(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT) << endl;
    cout.flush();
    cout << XPID << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT, "out");
    cout << __AFLOW_FUNC__ << " linking file RAW->RAW: " << string(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT) << " -> " << string(directory_RAW + "/" + aflowlib_data.auid + ".out") << endl;
    cout.flush();
    if (aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {
      aurostd::LinkFile(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_OUT, directory_RAW + "/" + aflowlib_data.auid + ".out"); // LINK //CO20200624 - adding FileExist() check
    }

    // write aflowlib.json
    cout << __AFLOW_FUNC__ << " writing " << string(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_JSON) << endl;
    cout.flush();
    cout << XPID << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ": " << aflowlib_data.aflowlib2file(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_JSON, "json");
    cout << __AFLOW_FUNC__ << " linking file RAW->RAW: " << string(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_JSON) << " -> " << string(directory_RAW + "/" + aflowlib_data.auid + ".json") << endl;
    cout.flush();
    if (aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_JSON)) {
      aurostd::LinkFile(directory_RAW + "/" + DEFAULT_FILE_AFLOWLIB_ENTRY_JSON, directory_RAW + "/" + aflowlib_data.auid + ".json"); // LINK  //CO20200624 - adding FileExist() check
    }

    if (flag_WEB) {
      if (!LOCAL) { // CO20171025
        if (aurostd::FileExist(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php")) {
          aurostd::LinkFile(DEFAULT_AFLOWDATA_WEB_DIRECTORY + "/api_index.php", directory_WEB + "/" + _XENTRY_); // CO20200624 - adding FileExist() check
        }
      }
    }

    // CHECK FOR AUID-WEB-LINKS
    if (!LOCAL) {
      if (LDEBUG) {
        cout << __AFLOW_FUNC__ << " flag_WEB=" << flag_WEB << endl;
      }
      if (LDEBUG) {
        cout << __AFLOW_FUNC__ << " directory_LIB=" << directory_LIB << endl;
      }
      if (LDEBUG) {
        cout << __AFLOW_FUNC__ << " directory_RAW=" << directory_RAW << endl;
      }
      if (LDEBUG) {
        cout << __AFLOW_FUNC__ << " directory_WEB=" << directory_WEB << endl;
      }
      // if(LDEBUG)
      cout << __AFLOW_FUNC__ << " aflowlib_data.auid=" << aflowlib_data.auid << endl;
      // **----------------------------------------------------------------------------
      if (XHOST.hostname == "nietzsche.mems.duke.edu" && (XHOST.user == "auro" || XHOST.user == "common" || XHOST.user == "stefano")) { // CO20200624 - cannot create these links otherwise
        // **----------------------------------------------------------------------------
        // NEW BUT STILL DOING
        string directory_AUID;
        string directory_AUID_LIB;
        string directory_AUID_RAW;
        string directory_AUID_WEB;
        directory_AUID = init::AFLOW_Projects_Directories("AUID") + "/" + aflowlib::auid2directory(aflowlib_data.auid);
        aurostd::DirectoryMake(directory_AUID);
        directory_AUID_LIB = directory_AUID + "/LIB";
        directory_AUID_RAW = directory_AUID + "/RAW";
        directory_AUID_WEB = directory_AUID + "/WEB";
        aurostd::RemoveFile(directory_AUID_LIB); // to avoid auto-linking SC20181205
        aurostd::RemoveFile(directory_AUID_RAW); // to avoid auto-linking SC20181205
        aurostd::RemoveFile(directory_AUID_WEB); // to avoid auto-linking SC20181205

        // directory_AUID_LIB
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " (AUID_NEW) directory_AUID_LIB=" << directory_AUID_LIB << " -> " << directory_LIB << endl;
        }
        cout << __AFLOW_FUNC__ << " (AUID_NEW) linking file AUID_LIB->LIB: " << directory_AUID_LIB << " -> " << directory_LIB << endl;
        cout.flush();
        if (aurostd::IsDirectory(directory_LIB)) {
          aurostd::LinkFile(directory_LIB, directory_AUID_LIB); // LINK  //CO20200624 - adding IsDirectory() check
        }
        // directory_AUID_RAW
        if (LDEBUG) {
          cout << __AFLOW_FUNC__ << " (AUID_NEW) directory_AUID_RAW=" << directory_AUID_RAW << " -> " << directory_RAW << endl;
        }
        cout << __AFLOW_FUNC__ << " (AUID_NEW) linking file AUID_RAW->LIB: " << directory_AUID_RAW << " -> " << directory_RAW << endl;
        cout.flush();
        if (aurostd::IsDirectory(directory_RAW)) {
          aurostd::LinkFile(directory_RAW, directory_AUID_RAW); // LINK  //CO20200624 - adding IsDirectory() check
        }

        if (flag_WEB) {
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " (AUID_NEW) directory_AUID_WEB=" << directory_AUID_WEB << " -> " << directory_WEB << endl;
          }
          cout << __AFLOW_FUNC__ << " (AUID_NEW) linking file AUID_WEB->LIB: " << directory_AUID_WEB << " -> " << directory_WEB << endl;
          cout.flush();
          if (aurostd::IsDirectory(directory_WEB)) {
            aurostd::LinkFile(directory_WEB, directory_AUID_WEB); // LINK  //CO20200624 - adding IsDirectory() check
          }
        } else {
          if (LDEBUG) {
            cout << __AFLOW_FUNC__ << " (AUID_NEW) directory_AUID_WEB=" << directory_AUID_WEB << " -> " << directory_RAW << endl;
          }
          cout << __AFLOW_FUNC__ << " (AUID_NEW) linking file AUID_WEB->LIB: " << directory_AUID_WEB << " -> " << directory_RAW << endl;
          cout.flush();
          if (aurostd::IsDirectory(directory_RAW)) {
            aurostd::LinkFile(directory_RAW, directory_AUID_WEB); // LINK  //CO20200624 - adding IsDirectory() check
          }
        }

        // ICSD2LINK
        if (aflowlib_data.catalog == "ICSD") {
          aurostd::string2tokens(directory_LIB, tokens, "_");
          if (tokens.size() > 2) {
            if (tokens[tokens.size() - 2] == "ICSD") {
              const string directory_ICSD2LINK = init::AFLOW_Projects_Directories("AUID") + "/icsd:/" + tokens[tokens.size() - 1];
              aurostd::DirectoryMake(directory_ICSD2LINK);
              cout << __AFLOW_FUNC__ << " (ICSD2LINK) making ICSD2LINK: " << directory_ICSD2LINK << endl;
              cout.flush();
              // LIB
              aurostd::RemoveFile(directory_ICSD2LINK + "/LIB"); // to avoid auto-linking SC20190830
              if (aurostd::IsDirectory(directory_LIB)) {
                aurostd::LinkFile(directory_LIB, directory_ICSD2LINK + "/LIB"); // LINK  //CO20200624 - adding IsDirectory() check
              }
              cout << __AFLOW_FUNC__ << " (ICSD2LINK) linking file LIB->ICSD2LINK/LIB: " << directory_LIB << " -> " << directory_ICSD2LINK << "/LIB" << endl;
              cout.flush();
              // RAW
              aurostd::RemoveFile(directory_ICSD2LINK + "/RAW"); // to avoid auto-linking SC20190830
              if (aurostd::IsDirectory(directory_RAW)) {
                aurostd::LinkFile(directory_RAW, directory_ICSD2LINK + "/RAW"); // LINK  //CO20200624 - adding IsDirectory() check
              }
              cout << __AFLOW_FUNC__ << " (ICSD2LINK) linking file RAW->ICSD2LINK/RAW: " << directory_RAW << " -> " << directory_ICSD2LINK << "/RAW" << endl;
              cout.flush();
              // WEB
              aurostd::RemoveFile(directory_ICSD2LINK + "/WEB"); // to avoid auto-linking SC20190830
              if (aurostd::IsDirectory(directory_WEB)) {
                aurostd::LinkFile(directory_WEB, directory_ICSD2LINK + "/WEB"); // LINK  //CO20200624 - adding IsDirectory() check
              }
              cout << __AFLOW_FUNC__ << " (ICSD2LINK) linking file WEB->ICSD2LINK/WEB: " << directory_WEB << " -> " << directory_ICSD2LINK << "/WEB" << endl;
              cout.flush();
            }
          }
        }
      }

      // DONE
      cout << __AFLOW_FUNC__ << " dir=" << directory_LIB << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << Message(__AFLOW_FILE__, aflags, "time")
           << " [time=" << aurostd::utype2string(aurostd::get_seconds(seconds_begin), 2, FIXED_STREAM) << "] " << endl; // CO20200624 - added FIXED_STREAM
      if (XHOST.vflag_control.flag("BEEP")) {
        aurostd::beep(aurostd::min(6000, aurostd::abs(int(1 * aflowlib_data.aflowlib2string().length() - 2000))), 50);
      }
    }

    // CHANGE PERMISSIONS
    // changing order of permission editing, if LOCAL, then order matters, otherwise NOT REALLY
    // files first, since we do /* (just in case directory_RAW is inside directory_LIB)
    // then directories
    // FILES

    const vector<string> Chmod_Files;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " pwd=" << aurostd::getPWD() << endl;
    }

    const std::filesystem::perms dir_perms =
        (std::filesystem::perms::owner_all | std::filesystem::perms::group_read | std::filesystem::perms::group_exec | std::filesystem::perms::others_read | std::filesystem::perms::others_exec);
    const std::filesystem::perms file_perms = (std::filesystem::perms::owner_read | std::filesystem::perms::owner_write | std::filesystem::perms::group_read | std::filesystem::perms::others_read);
    for (const auto& dir_entry : std::filesystem::directory_iterator(aurostd::CleanFileName(directory_RAW))) {
      if (dir_entry.is_directory()) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << "chmod " << static_cast<unsigned>(dir_perms) << " " << dir_entry.path().filename().string() << endl;
        }
        std::filesystem::permissions(dir_entry.path().string(), dir_perms);
      } else if (dir_entry.is_regular_file()) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << "chmod " << static_cast<unsigned>(file_perms) << " " << dir_entry.path().filename().string() << endl;
        }
        std::filesystem::permissions(dir_entry.path().string(), file_perms);
      }
    }

    if (CHMODWEB) {
      for (const auto& dir_entry : std::filesystem::directory_iterator(aurostd::CleanFileName(directory_WEB))) {
        if (dir_entry.is_directory()) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << "chmod " << static_cast<unsigned>(dir_perms) << " " << dir_entry.path().filename().string() << endl;
          }
          std::filesystem::permissions(dir_entry.path().string(), dir_perms);
        } else if (dir_entry.is_regular_file()) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << "chmod " << static_cast<unsigned>(file_perms) << " " << dir_entry.path().filename().string() << endl;
          }
          std::filesystem::permissions(dir_entry.path().string(), file_perms);
        }
      }
    }

    for (const auto& dir_entry : std::filesystem::directory_iterator(aurostd::CleanFileName(directory_LIB))) {
      if (dir_entry.is_directory()) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << "chmod " << static_cast<unsigned>(dir_perms) << " " << dir_entry.path().filename().string() << endl;
        }
        std::filesystem::permissions(dir_entry.path().string(), dir_perms);
      } else if (dir_entry.is_regular_file()) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << "chmod " << static_cast<unsigned>(file_perms) << " " << dir_entry.path().filename().string() << endl;
        }
        std::filesystem::permissions(dir_entry.path().string(), file_perms);
      }
    }
    // done
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_FileNeeded
// ***************************************************************************
namespace aflowlib {
  /// @brief check for files that need to exist in the library, if they exist then copy to raw
  /// @param directory_LIB lib directory to check
  /// @param fileLIB lib file to check
  /// @param directory_RAW raw directory to copy to
  /// @param fileRAW raw file to copy to
  /// @param vfile[out] add the copied file to this list
  /// @param MESSAGE message to use for errors when the lib file doesn't exist
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature}
  void LIB2RAW_FileNeeded(const string& directory_LIB, const string& fileLIB, const string& directory_RAW, const string& fileRAW, vector<string>& vfile, const string& MESSAGE) {
    string error_message;

    if (XHOST.vext.size() != XHOST.vzip.size()) {
      error_message = "XHOST.vext.size()!=XHOST.vzip.size()";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_message, _INDEX_MISMATCH_);
    }

    const string file_LIB = directory_LIB + "/" + fileLIB;
    const string file_RAW = directory_RAW + "/" + fileRAW;
    string file_LIB_nocompress = directory_LIB + "/" + fileLIB;
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) {
      aurostd::StringSubstInPlace(file_LIB_nocompress, XHOST.vext[iext], ""); // SKIP uncompressed
    }
    string file_RAW_nocompress = directory_RAW + "/" + fileRAW;
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) {
      aurostd::StringSubstInPlace(file_RAW_nocompress, XHOST.vext[iext], ""); // SKIP uncompressed
    }
    if (aurostd::FileExist(file_RAW)) {
      return; // already there
    }
    if (aurostd::FileExist(file_RAW_nocompress)) {
      return; // already there
    }

    if (!aurostd::FileExist(file_LIB) && !aurostd::FileExist(file_LIB_nocompress) && !aurostd::CompressFileExist(file_LIB_nocompress)) {
      if (!aurostd::FileExist(file_LIB)) {
        error_message = MESSAGE + " file not found " + file_LIB;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_message, _FILE_NOT_FOUND_);
      }
      if (!aurostd::FileExist(file_LIB_nocompress)) {
        error_message = MESSAGE + " file not found " + file_LIB_nocompress;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_message, _FILE_NOT_FOUND_);
      }
      if (!aurostd::CompressFileExist(file_LIB_nocompress)) {
        error_message = MESSAGE + " file not found " + file_LIB_nocompress + ".EXT";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, error_message, _FILE_NOT_FOUND_);
      }
    }
    if (aurostd::FileExist(file_LIB)) {
      aurostd::CopyFile(file_LIB, file_RAW);
    }
    if (aurostd::FileExist(file_LIB_nocompress)) {
      aurostd::CopyFile(file_LIB_nocompress, file_RAW_nocompress);
    }
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      if (aurostd::FileExist(file_LIB_nocompress + XHOST.vext[iext])) {
        aurostd::CopyFile(file_LIB_nocompress + XHOST.vext[iext], file_RAW_nocompress + XHOST.vext[iext]);
      }
    }
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      if (aurostd::FileExist(file_RAW) && aurostd::substring2bool(file_RAW, XHOST.vext[iext])) {
        aurostd::DecompressFile(file_RAW);
      }
      if (aurostd::FileExist(file_RAW_nocompress + XHOST.vext[iext])) {
        aurostd::DecompressFile(file_RAW_nocompress + XHOST.vext[iext]);
      }
    }
    string file2add = fileRAW;
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      aurostd::StringSubstInPlace(file2add, XHOST.vext[iext], "");
    }
    vfile.push_back(file2add);
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Static
// ***************************************************************************
namespace aflowlib {
  /// @brief runs the LIB2RAW subroutine for static files
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_Static(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = XHOST.DEBUG;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }

    cout << MESSAGE << " " << __AFLOW_FUNC__ << " begin " << directory_LIB << endl;
    cout << MESSAGE << " " << __AFLOW_FUNC__ << " species = " << data.vspecies.size() << endl;
    data.vloop.emplace_back("static");

    vector<string> files{_AFLOWIN_, "DOSCAR.static", "OUTCAR.static", "OSZICAR.static", "KPOINTS.static"};
    // ONLY GET ONE POSCAR, bands is better as we need it later
    if (aurostd::CompressFileExist(directory_LIB + "/POSCAR.bands")) {
      files.emplace_back("POSCAR.bands");
    } else {
      files.emplace_back("POSCAR.static");
    }
    for (const auto& file : files) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    }

    const bool flag_use_GNUPLOT = true;

    if (flag_use_GNUPLOT) {
      // ME20190614 BEGIN
      //  This has to come first because FIXBANDS messes up the EIGENVAL files
      aurostd::xoption opts;
      aurostd::xoption plotoptions;
      const string dosscale = aurostd::utype2string<double>(DEFAULT_DOS_SCALE);
      opts.push_attached("PLOT_DOS", directory_RAW + ",,," + dosscale);
      opts.push_attached("PLOT_PDOS", directory_RAW + ",-1,,," + dosscale);
      opts.push_attached("PLOTTER::PRINT", "png");
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_DOS");
      plotter::PLOT_DOS(plotoptions);
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_PDOS", true);
      plotter::PLOT_PDOS(plotoptions);
      // ME20190614 END
    }

    xKPOINTS kpoints_static;
    kpoints_static.GetPropertiesFile(directory_RAW + "/KPOINTS.static");
    data.kpoints_nnn_static = kpoints_static.nnn_kpoints;
    if (!data.kpoints.empty()) {
      data.kpoints += ";";
    } // CO20220706
    data.kpoints += aurostd::utype2string(kpoints_static.nnn_kpoints[1]) + "," + aurostd::utype2string(kpoints_static.nnn_kpoints[2]) + "," + aurostd::utype2string(kpoints_static.nnn_kpoints[3]);

    cout << MESSAGE << " " << __AFLOW_FUNC__ << " end " << directory_LIB << endl;
  }
} // namespace aflowlib
// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bands
// ***************************************************************************
namespace aflowlib {
  /// @brief runs the LIB2RAW subroutine for bands files
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_Bands(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = XHOST.DEBUG;
    // LDEBUG=true;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }

    cout << MESSAGE << " " << __AFLOW_FUNC__ << " begin " << directory_LIB << endl;
    cout << MESSAGE << " " << __AFLOW_FUNC__ << " species = " << data.vspecies.size() << endl;
    data.vloop.emplace_back("bands");

    stringstream command;
    command.clear();
    command.str(std::string());
    // directories must exist already
    const bool flag_DATA_BANDS_ = false;
    const bool flag_use_MATLAB = false;
    const bool flag_use_GNUPLOT = !flag_use_MATLAB; // KY

    // copy _AFLOWIN_ LOCK DOSCAR.static.EXT EIGENVAL.bands.EXT KPOINTS.bands.EXT POSCAR.bands.EXT
    const vector<string> files{_AFLOWIN_, _AFLOWLOCK_, "DOSCAR.static", "OUTCAR.static", "OSZICAR.static", "OSZICAR.bands", "EIGENVAL.bands", "KPOINTS.static", "KPOINTS.bands",
                               // "INCAR.static",  // not needed but good for show off SC 0914
                               "INCAR.bands",
                               // "POSCAR.relax1",
                               "POSCAR.bands"};
    for (const auto& file : files) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    }

    xKPOINTS kpoints_bands;
    kpoints_bands.GetPropertiesFile(directory_RAW + "/KPOINTS.bands");
    // get pairs
    data.kpoints_pairs.clear();
    if (kpoints_bands.vpath.size() % 2 == 0) { // if even
      for (size_t i = 0; i < kpoints_bands.vpath.size(); i += 2) {
        data.kpoints_pairs.emplace_back(kpoints_bands.vpath[i] + "-" + kpoints_bands.vpath[i + 1]);
      }
    }
    data.kpoints_bands_path_grid = kpoints_bands.path_grid;
    if (!data.kpoints.empty()) {
      data.kpoints += ";";
    } // CO20220706
    data.kpoints += kpoints_bands.path + ";" + aurostd::utype2string(kpoints_bands.path_grid);
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " KPOINTS = " << data.kpoints << endl;
    }

    if (flag_use_MATLAB) { // MATLAB STUFF  OLD WSETYAWAN+SC
      // PERFORM THE MATLAB STEP
      cout << MESSAGE << " MATLAB start: " << directory_RAW << endl;
      // WRITE plotbz.sh
      stringstream gnuplot_plotbz;
      gnuplot_plotbz.clear();
      gnuplot_plotbz.str(std::string());
      gnuplot_plotbz << aurostd::EmbData::get_content("GNUPLOT_FUNCS_plotbz.sh", "SCRIPTS") << endl;
      aurostd::stringstream2file(gnuplot_plotbz, string(directory_RAW + "/plotbz.sh"));

      // WRITE PARAM.M
      stringstream matlab_param;
      matlab_param.clear();
      matlab_param.str(std::string());
      matlab_param << aurostd::EmbData::get_content("MATLAB_FUNCS_param.mat", "SCRIPTS") << endl;
      // WRITE PLOTBAND.M
      stringstream matlab_plotband;
      matlab_plotband.clear();
      matlab_plotband.str(std::string());
      matlab_plotband << aurostd::EmbData::get_content("MATLAB_FUNCS_plotband.mat", "SCRIPTS") << endl; // normal OR log
      matlab_plotband << "exit;" << endl;
      aurostd::stringstream2file(matlab_plotband, string(directory_RAW + "/plotband.m"));

      // NEW STUFF
      aurostd::file2file(directory_RAW + "/KPOINTS.bands", directory_RAW + "/KPOINTS.bands.old");
      vfiles.emplace_back("KPOINTS.bands.old");
      aurostd::file2file(directory_RAW + "/EIGENVAL.bands", directory_RAW + "/EIGENVAL.bands.old");
      vfiles.emplace_back("EIGENVAL.bands.old");
      _aflags aflags;
      aflags.Directory = directory_RAW;

      if (pflow::FIXBANDS(aflags, "POSCAR.bands,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands,EIGENVAL.bands") == false) {
        cout << "ERROR_RERUN " << directory_LIB << endl;
        return;
      }

      // EXECUTE PLOTBZ.SH using ksh
      command.clear();
      command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);

      // EXECUTE MATLAB
      command.clear();
      command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "export DISPLAY=:0.0" << endl;
      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
      command << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("plotband") << endl;
      aurostd::execute(command);
    }

    if (flag_use_GNUPLOT) { // GNUPLOT STUFF NEW KY+SC
      // ME20190614 BEGIN
      //  This has to come first because FIXBANDS messes up the EIGENVAL files
      aurostd::xoption opts;
      aurostd::xoption plotoptions;
      const string dosscale = aurostd::utype2string<double>(DEFAULT_DOS_SCALE);
      opts.push_attached("PLOT_DOS", directory_RAW + ",,," + dosscale);
      opts.push_attached("PLOT_BANDDOS", directory_RAW + ",,," + dosscale);
      opts.push_attached("PLOT_PDOS", directory_RAW + ",-1,,," + dosscale);
      opts.push_attached("PLOTTER::PRINT", "png");
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_DOS");
      plotter::PLOT_DOS(plotoptions);
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_BANDDOS");
      plotter::PLOT_BANDDOS(plotoptions);
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_PDOS", true);
      plotter::PLOT_PDOS(plotoptions);
      // ME20190614 END

      // KY WRITE THE CODE HERE
      cout << MESSAGE << " GNUPLOT start: " << directory_RAW << endl;
      // WRITE plotbz.sh
      stringstream gnuplot_plotbz;
      gnuplot_plotbz.clear();
      gnuplot_plotbz.str(std::string());
      gnuplot_plotbz << aurostd::EmbData::get_content("GNUPLOT_FUNCS_plotbz.sh", "SCRIPTS") << endl;
      aurostd::stringstream2file(gnuplot_plotbz, string(directory_RAW + "/plotbz.sh"));

      // NEW STUFF
      aurostd::file2file(directory_RAW + "/KPOINTS.bands", directory_RAW + "/KPOINTS.bands.old");
      vfiles.emplace_back("KPOINTS.bands.old");
      aurostd::file2file(directory_RAW + "/EIGENVAL.bands", directory_RAW + "/EIGENVAL.bands.old");
      vfiles.emplace_back("EIGENVAL.bands.old");
      _aflags aflags;
      aflags.Directory = directory_RAW;
      if (pflow::FIXBANDS(aflags, "POSCAR.bands,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands,EIGENVAL.bands") == false) {
        cout << "ERROR_RERUN " << directory_LIB << endl;
        return;
      }

      // EXECUTE PLOTBZ.SH using ksh
      command.clear();
      command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);

      // EXECUTE PLOTBAND
    }

    // DONE
    cout << MESSAGE << " " << __AFLOW_FUNC__ << " end " << directory_LIB << endl;
  }
} // namespace aflowlib
// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bands
// ***************************************************************************
namespace aflowlib {
  /// @brief runs the LIB2RAW subroutine for bands files
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{SD,20250409,created function}
  void LIB2RAW_Loop_Dielectric(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = XHOST.DEBUG;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }

    cout << MESSAGE << " " << __AFLOW_FUNC__ << " begin " << directory_LIB << endl;
    cout << MESSAGE << " " << __AFLOW_FUNC__ << " species = " << data.vspecies.size() << endl;
    data.vloop.emplace_back("dielectric");

    const vector<string> files{_AFLOWIN_, "OUTCAR.dielectric"};
    for (const auto& file : files) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    }

    xOUTCAR xOUT;
    xOUT.GetPropertiesFile(directory_RAW + "/OUTCAR.dielectric");
    data.freq_plasma = xOUT.freq_plasma;
    data.dielectric_static = xOUT.dielectric_static;

    cout << MESSAGE << " " << __AFLOW_FUNC__ << " end " << directory_LIB << endl;
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Thermodynamics
// ***************************************************************************

namespace aflowlib {
  /// @brief appends a suffix to the base of a filename while keeping the file extension at the end
  /// @param _file the file to append to
  /// @param addendum the suffix to add the base file name
  /// @param out_file the resulting file name
  /// @return true on success, false on failure
  /// @authors
  /// @mod{CO,20171025,created function}
  /// @mod{ST,20241022,created doxy\, optimize\, cleanup}
  bool AddFileNameBeforeExtension(const string& _file, const string& addendum, string& out_file) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    std::filesystem::path path = aurostd::CleanFileName(_file);
    out_file = path.string();
    if (is_directory(path)) {
      return false;
    } // not doing directories
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " path=" << path << endl;
      cerr << __AFLOW_FUNC__ << " file=" << path.filename() << endl;
    }
    out_file = aurostd::CleanFileName(path.replace_extension(addendum + path.extension().string()));
    return true;
  }
} // namespace aflowlib

namespace aflowlib {
  /// @brief the LIB2RAW subroutine for caluclating formation enthalpies
  /// @param data data structure for the library entry
  /// @param xstr the structure for the entry
  /// @param MESSAGE message to use for logging errors
  /// @return true if formation enthalpy was calculated
  /// @authors
  /// @mod{CO,20200731,created function}
  /// @mod{ST,20241022,created doxy\, cleanup}
  bool LIB2RAW_Calculate_FormationEnthalpy(aflowlib::_aflowlib_entry& data, const xstructure& xstr, const string& MESSAGE) {
    //   make 2 flags
    //   formation_calc_cce which allows or not calculaiton of cce having the right functional-cce
    //   formation_calc_U which allows or not calculation of Href based on the fact that we have - or not - U
    //   then a combination of AND/OR will do the if(combination) so that  lines XXX1 and XXX2 are always printed

    const bool LDEBUG = (false || XHOST.DEBUG);
    // reference
    bool FORMATION_CALC = true;
    bool formation_calc_U = false;
    if (!xstr.species_pp_vLDAU.empty() && !xstr.species_pp_vLDAU[0].empty() && !aurostd::isequal(xstr.species_pp_vLDAU[0][0], 0.0)) {
      formation_calc_U = true; // CO20210712 - new style, vLDAU will be populated even if no +U, check type==0
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " formation_calc_U=" << formation_calc_U << endl;
    }

    // CO20200624 START - CCE
    bool formation_calc_cce = false;
    if (aurostd::WithinList(data.vspecies, "O")) {
      formation_calc_cce = true;
    }
    if (aurostd::WithinList(data.vspecies, "N")) {
      formation_calc_cce = true;
    }
    string functional_cce;
    if (formation_calc_cce) {
      // CO20200624 - don't use else if's, we might prefer that which comes later
      // CCE cannot correct 'GGA'
      if (data.dft_type.find("PBE") != string::npos) {
        functional_cce = "PBE";
      }
      if (data.dft_type.find("LDA") != string::npos) {
        functional_cce = "LDA";
      }
      if (functional_cce == "PBE") {
        if (data.dft_type.find("SCAN") != string::npos) {
          functional_cce = "SCAN";
        } // PAW_PBE_KIN:SCAN
        if (data.catalog == "ICSD" && formation_calc_U == true) {
          functional_cce = "PBE+U:ICSD";
        }
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " FUNCTIONAL_CCE=" << functional_cce << endl;
    }
    if (formation_calc_cce && functional_cce.empty()) {
      formation_calc_cce = false;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " formation_calc_cce=" << formation_calc_cce << endl;
    }
    // CO20200624 STOP - CCE

    // LDAU ?
    if (FORMATION_CALC) { // CO20200624 - we can correct PBE now with CCE
      //[CO20200624 - we can correct PBE now with CCE]if(formation_calc_U) FORMATION_CALC=false; // no formation for LDAU
      if (formation_calc_U && !(formation_calc_cce && functional_cce == "PBE+U:ICSD")) {
        FORMATION_CALC = false;
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " FORMATION_CALC=" << FORMATION_CALC << endl;
    }

    // PP AVAILABLE ?
    // line XXX1
    if (FORMATION_CALC) {
      for (uint i = 0; i < (uint) data.nspecies; i++) {
        if (!xPOTCAR_EnthalpyReference_AUID(data.vspecies_pp_AUID.at(i), data.METAGGA)) {
          FORMATION_CALC = false; // WORKS WITH SCANN TOO
          if (AFLOWLIB_VERBOSE) {
            cout << MESSAGE << " REFERENCE NOT AVAILABLE species_pp=" << xstr.species_pp_version.at(i) << "  species_pp_AUID=" << data.vspecies_pp_AUID.at(i)
                 << (data.METAGGA.empty() ? string("") : string("  METAGGA=[" + data.METAGGA + "]")) << "   Href=nothing" << endl;
          }
        }
      }
    }

    // OPERATE FORMATION
    // line XXX2
    uint i = 0;
    if (FORMATION_CALC) { // no LDAU yet
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [FCALC=1]" << endl;
      }
      vector<double> venthalpy_atom_ref;
      double enthalpy_atom_ref = data.enthalpy_atom; // if there is 1 then there is only one
      string aus_gs_structure;
      double aus_gs_atom = 0.0;
      double aus_volume_atom = 0.0;
      double aus_spin_atom = 0.0;
      for (i = 0; i < (uint) data.nspecies; i++) {
        //      string pseudopotential,string type,vector<double> LDAU
        enthalpy_atom_ref = data.enthalpy_atom; // if there is 1 then there is only one
        aus_gs_structure = "";
        // NEW STYLE
        xPOTCAR_EnthalpyReference_AUID(data.vspecies_pp_AUID.at(i), data.METAGGA, aus_gs_structure, aus_gs_atom, aus_volume_atom, aus_spin_atom);
        enthalpy_atom_ref = aus_gs_atom;
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " REFERENCE species=" << xstr.species.at(i) << " species_pp_AUID=" << data.vspecies_pp_AUID.at(i) << (data.METAGGA.empty() ? string("") : string("  METAGGA=[" + data.METAGGA + "]"))
               << "   Href=" << enthalpy_atom_ref << endl;
        }
        venthalpy_atom_ref.push_back(enthalpy_atom_ref);
      }

      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [FCALC=2]" << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [FCALC=2] data.nspecies=" << data.nspecies << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [FCALC=2] venthalpy_atom_ref.size()=" << venthalpy_atom_ref.size() << endl;
      }

      // calculation of REF
      data.enthalpy_formation_cell = data.enthalpy_cell;
      data.enthalpy_formation_atom = data.enthalpy_atom;

      double data_natoms = 0.0; // needs to be double for pocc
      for (i = 0; i < xstr.comp_each_type.size(); i++) {
        data_natoms += xstr.comp_each_type[i];
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " data_natoms=" << data_natoms << endl;
      }

      for (i = 0; i < (uint) data.nspecies; i++) {
        data.enthalpy_formation_cell = data.enthalpy_formation_cell - (double(venthalpy_atom_ref.at(i) * xstr.comp_each_type.at(i)));
      }
      data.enthalpy_formation_atom = data.enthalpy_formation_cell / data_natoms;

      // CO20200624 START - adding cce variants
      if (formation_calc_cce && !functional_cce.empty()) {
        vector<double> enthalpy_formation_cell_corrections_cce;
        try {
          enthalpy_formation_cell_corrections_cce = cce::calculate_corrections(xstr, functional_cce);
        } catch (aurostd::xerror& excpt) {
          pflow::logger(excpt.whereFileName(), excpt.whereFunction(), excpt.buildMessageString(), cout, _LOGGER_ERROR_);
          formation_calc_cce = false;
        }
        if (formation_calc_cce && enthalpy_formation_cell_corrections_cce.size() == 2) { // the first is at 300K, the second at 0K
          data.enthalpy_formation_cce_300K_cell = data.enthalpy_formation_cce_0K_cell = data.enthalpy_formation_cell;
          data.enthalpy_formation_cce_300K_cell -= enthalpy_formation_cell_corrections_cce[0];
          data.enthalpy_formation_cce_0K_cell -= enthalpy_formation_cell_corrections_cce[1];
          data.enthalpy_formation_cce_300K_atom = data.enthalpy_formation_cce_300K_cell / data_natoms;
          data.enthalpy_formation_cce_0K_atom = data.enthalpy_formation_cce_0K_cell / data_natoms;
        }
        if (formation_calc_cce == false && formation_calc_U == true) { // CO20210828 - if CCE fails for +U calc, clear out data.enthalpy_formation_cell and _atom
          if (AFLOWLIB_VERBOSE) {
            cout << MESSAGE << " CCE failed for +U calculation, clearing out data.enthalpy_formation_cell and _atom" << endl;
          }
          data.enthalpy_formation_cell = AUROSTD_NAN;
          data.enthalpy_formation_atom = AUROSTD_NAN;
          FORMATION_CALC = false;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " data.enthalpy_formation_cell=" << data.enthalpy_formation_cell << endl;
      cerr << __AFLOW_FUNC__ << " data.enthalpy_formation_atom=" << data.enthalpy_formation_atom << endl;
      cerr << __AFLOW_FUNC__ << " data.enthalpy_formation_cce_300K_cell=" << data.enthalpy_formation_cce_300K_cell << endl;
      cerr << __AFLOW_FUNC__ << " data.enthalpy_formation_cce_0K_cell=" << data.enthalpy_formation_cce_0K_cell << endl;
      cerr << __AFLOW_FUNC__ << " data.enthalpy_formation_cce_300K_atom=" << data.enthalpy_formation_cce_300K_atom << endl;
      cerr << __AFLOW_FUNC__ << " data.enthalpy_formation_cce_0K_atom=" << data.enthalpy_formation_cce_0K_atom << endl;
    }
    // CO20200624 STOP - adding cce variants

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [FCALC=3]" << endl;
    }

    if (FORMATION_CALC) {
      data.entropic_temperature = 0;
      if (data.vstoichiometry.size() > 1) {
        for (i = 0; i < (size_t) data.vstoichiometry.size(); i++) {
          if (data.vstoichiometry[i] > _EPSILON_COMPOSITION_ && data.vstoichiometry[i] < 1 - _EPSILON_COMPOSITION_) {
            data.entropic_temperature += data.vstoichiometry[i] * logl(data.vstoichiometry[i]);
          }
        }
        data.entropic_temperature = data.enthalpy_formation_atom / (data.entropic_temperature * KBOLTZEV);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [FCALC=4]" << endl;
    }
    return FORMATION_CALC;
  }
} // namespace aflowlib

namespace aflowlib {
  /// @brief LIB2RAW subroutine for thermodynamics results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @param LOCAL whether to run locally
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_Thermodynamics(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE, bool LOCAL) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream messagestream; // dummy stringstream, can output to cout, but not used right now

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [-1]" << endl;
    }
    // ZIP-AGNOSTIC
    if (XHOST.vext.size() != XHOST.vzip.size()) {
      messagestream << "XHOST.vext.size()!=XHOST.vzip.size()";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _INDEX_MISMATCH_);
    }
    // CO+DX START 20170713 - adding symmetry output to RAW
    _aflags aflags;
    aflags.Directory = directory_RAW;
    ofstream FileMESSAGE; // dummy ofstream, not really used
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [-2]" << endl;
    }
    //[CO20200829 - inside data]string system_name=KBIN::ExtractSystemName(directory_LIB);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [-3]" << endl;
    }
    // CO+DX STOP 20170713 - adding symmetry output to RAW
    deque<string> deq_species;
    aurostd::string2tokens(data.species, deq_species, ","); // DX20190620
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "directory_LIB=\"" << directory_LIB << "\"" << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "directory_RAW=\"" << directory_RAW << "\"" << endl;
    }
    //[CO20200624 - MOVING UP]if(directory_LIB.at(directory_LIB.size()-1)=='/')  directory_LIB=directory_LIB.substr(0,directory_LIB.size()-1);
    //[CO20200624 - MOVING UP]if(directory_RAW.at(directory_RAW.size()-1)=='/')  directory_RAW=directory_RAW.substr(0,directory_RAW.size()-1);
    if (LDEBUG) {
      cerr << XPID << "directory_LIB=\"" << directory_LIB << "\"" << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "directory_RAW=\"" << directory_RAW << "\"" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - species = " << data.vspecies.size() << endl;
    }
    const uint relax_max = RELAX_MAX;
    stringstream command;
    command.clear();
    command.str(std::string());
    // directories must exist already
    bool flag_EDATA_ORIG_ = true;
    bool flag_EDATA_RELAX_ = true;
    bool flag_EDATA_BANDS_ = false;
    bool flag_DATA_ORIG_ = false;
    bool flag_DATA_RELAX_ = false;
    bool flag_DATA_BANDS_ = false;
    const bool flag_TIMING = true;
    const bool flag_ENERGY1 = false;
    const bool flag_MAGNETIC = true;
    bool flag_SG1 = true;
    bool flag_SG2 = true;
    bool flag_ICSD = false;
    bool flag_LIB0 = false;
    const bool flag_LIB1 = false;
    bool flag_LIB2 = false;
    bool flag_ERROR = false;

    struct file {
      string lib;
      string raw;
    };
    std::unordered_map<string, file> filecar_map; // car-->(lib,raw)
    string tmp_key;
    string FileName_OUTCAR_relax;

    if (LOCAL) {
      flag_EDATA_BANDS_ = true;
      flag_DATA_ORIG_ = true;
      flag_DATA_RELAX_ = true;
      flag_DATA_BANDS_ = true;
    } else {
      if (aurostd::substring2bool(directory_LIB, "ICSD")) {
        flag_ICSD = true;
      }
      if (aurostd::substring2bool(directory_LIB, "LIB0")) {
        flag_EDATA_ORIG_ = false;
        flag_EDATA_RELAX_ = false;
        flag_SG1 = false;
        flag_SG2 = false;
      }
      if (aurostd::substring2bool(directory_LIB, "LIB2") || aurostd::substring2bool(directory_LIB, "LIBRARYX")) {
        flag_LIB2 = true;
      }
    }

    flag_LIB0 = aurostd::substring2bool(directory_LIB, "LIB0");

    // check for flag_EDATA_BANDS_
    // check for flag_DATA_BANDS_
    tmp_key = "POSCAR";
    if ((flag_EDATA_ORIG_ && flag_EDATA_RELAX_) || (flag_DATA_ORIG_ && flag_DATA_RELAX_)) {
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        filecar_map[tmp_key].lib = string(directory_LIB).append("/").append(tmp_key).append(".bands").append(XHOST.vext[iext]);
        filecar_map[tmp_key].raw = string(directory_RAW).append("/").append(tmp_key).append(".bands").append(XHOST.vext[iext]);
        if (aurostd::FileExist(filecar_map[tmp_key].lib)) {
          flag_EDATA_BANDS_ = flag_EDATA_ORIG_ && flag_EDATA_RELAX_;
          flag_DATA_BANDS_ = flag_DATA_ORIG_ && flag_DATA_RELAX_;
          aurostd::CopyFile(filecar_map[tmp_key].lib, filecar_map[tmp_key].raw);
          aurostd::DecompressFile(filecar_map[tmp_key].raw);
          vfiles.emplace_back(tmp_key + ".bands");
        }
      }
    }
    if (flag_ICSD) {
      ;
    } // dummy load
    if (flag_MAGNETIC) {
      ;
    } // dummy load
    if (flag_LIB0) {
      ;
    } // dummy load
    if (flag_LIB1) {
      ;
    } // dummy load

    xstructure str_orig;
    xstructure str_relax1;
    xstructure str_relax;

    vector<string> tokens;

    double data1_energy_cell = 0.0;
    double data1_energy_atom = 0.0;
    string data_sg1_pre;
    string data_sg1_mid;
    string data_sg1_post;
    string data_sg2_pre;
    string data_sg2_mid;
    string data_sg2_post;
    xvector<double> data_abcabc;

    xOUTCAR xOUT;
    // xDOSCAR xDOS;
    // xEIGENVAL xEIG;
    xKPOINTS kpoints;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [2]" << endl;
    }
    // copy _AFLOWIN_ LOCK
    // _AFLOWIN_
    data.vloop.emplace_back("thermodynamics");
    // star

    if (!aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "LOCK file does not exist", _FILE_NOT_FOUND_);
    }

    // FILES
    aflowlib::LIB2RAW_FileNeeded(directory_LIB, _AFLOWIN_, directory_RAW, _AFLOWIN_, vfiles, MESSAGE); // _AFLOWIN_
    aflowlib::LIB2RAW_FileNeeded(directory_LIB, _AFLOWLOCK_, directory_RAW, _AFLOWLOCK_, vfiles, MESSAGE); // LOCK

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [3]" << endl;
    }
    if (true || flag_DATA_ORIG_ || flag_EDATA_ORIG_ || flag_SG1 || flag_SG2) { // POSCAR.orig.EXT
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " [3.1]" << endl;
      }
      // if(flag_ORIG==false) { //[CO20200106 - close bracket for indenting]}
      bool found = false;
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " BUILDING POSCAR.orig from POSCAR.orig.EXT" << endl;
        }
        if (!found && aurostd::FileExist(directory_LIB + "/POSCAR.orig" + XHOST.vext[iext])) {
          found = true;
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " BUILDING POSCAR.orig from POSCAR.orig" << XHOST.vext[iext] << endl;
          }
          aflowlib::LIB2RAW_FileNeeded(directory_LIB, "POSCAR.orig", directory_RAW, "POSCAR.orig", vfiles, MESSAGE); // POSCAR.orig
        }
      }
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " BUILDING POSCAR.orig from POSCAR.relax1.EXT" << endl;
        }
        if (!found && aurostd::FileExist(directory_LIB + "/POSCAR.relax1" + XHOST.vext[iext])) {
          found = true;
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " BUILDING POSCAR.orig from POSCAR.relax1" << XHOST.vext[iext] << endl;
          }
          aflowlib::LIB2RAW_FileNeeded(directory_LIB, "POSCAR.relax1", directory_RAW, "POSCAR.orig", vfiles, MESSAGE); // POSCAR.orig
        }
      }
      if (!found) {
        found = true;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " BUILDING POSCAR.orig from " << _AFLOWIN_ << "" << endl;
        }
        const string out = aurostd::execute2string(string("cat ") + "\"" + directory_RAW + "/" + _AFLOWIN_ + "\"" + R"( | aflow --justbetween="[VASP_POSCAR_MODE_EXPLICIT]START","[VASP_POSCAR_MODE_EXPLICIT]STOP")");
        aurostd::string2file(out, directory_RAW + "/POSCAR.orig", aurostd::compression_type::None, "WRITE");
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [4]" << endl;
    }
    if (!flag_LIB0) { // no LIB0
      if (true || flag_DATA_RELAX_ || flag_EDATA_RELAX_ || flag_SG1 || flag_SG2) { // CONTCAR.relax.EXT
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          for (uint i = 1; i <= relax_max; i++) {
            const vector<string> vcarkeys{"CONTCAR", "OUTCAR", "KPOINTS", "INCAR"};
            for (const auto& key : vcarkeys) {
              filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".relax").append(aurostd::utype2string(i)).append(XHOST.vext[iext]));
              filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax").append(XHOST.vext[iext]));
            }
            tmp_key = "CONTCAR";
            if (aurostd::FileExist(filecar_map[tmp_key].lib)) {
              aurostd::CopyFile(filecar_map[tmp_key].lib, filecar_map[tmp_key].raw);
              aurostd::DecompressFile(filecar_map[tmp_key].raw);
              vfiles.emplace_back(tmp_key + ".relax");
            }
            tmp_key = "OUTCAR";
            if (aurostd::FileExist(filecar_map[tmp_key].lib)) {
              aurostd::CopyFile(filecar_map[tmp_key].lib, filecar_map[tmp_key].raw);
              aurostd::DecompressFile(filecar_map[tmp_key].raw);
              vfiles.emplace_back(tmp_key + ".relax");
              FileName_OUTCAR_relax = filecar_map[tmp_key].lib;
            }
          }
        }
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          if (!aurostd::FileExist(directory_RAW + "/CONTCAR.relax") && aurostd::FileExist(directory_LIB + "/CONTCAR.static" + XHOST.vext[iext]) && !aurostd::FileExist(directory_RAW + "/OUTCAR.relax") &&
              aurostd::FileExist(directory_LIB + "/OUTCAR.static" + XHOST.vext[iext]) && !aurostd::FileExist(directory_RAW + "/KPOINTS.relax") && aurostd::FileExist(directory_LIB + "/KPOINTS.static" + XHOST.vext[iext])) {
            const vector<string> vcarkeys{"CONTCAR", "OUTCAR", "KPOINTS"};
            for (const auto& key : vcarkeys) {
              filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".static").append(XHOST.vext[iext]));
              filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax").append(XHOST.vext[iext]));
              cout << MESSAGE << " WARNING - PATCHING " << key << ".relax with " << key << ".static " << filecar_map[key].lib << endl;
              if (aurostd::FileExist(filecar_map[key].lib)) {
                aurostd::CopyFile(filecar_map[key].lib, filecar_map[key].raw);
                aurostd::DecompressFile(filecar_map[key].raw);
                vfiles.emplace_back(key + ".relax");
                if (key == "OUTCAR") {
                  FileName_OUTCAR_relax = filecar_map[key].lib;
                }
              }
            }
          }
        }
        vector<string> vcarkeys{"CONTCAR", "OUTCAR", "KPOINTS", "INCAR"};
        for (size_t i = 0; i < vcarkeys.size(); i++) {
          const string& key = vcarkeys[i];
          if (aurostd::FileExist(filecar_map[key].raw)) {
            messagestream << MESSAGE << " [" << i + 1 << "] - file not prepared " << filecar_map[key].lib;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _FILE_ERROR_);
          }
        }
      }
    } // no LIB0

    if (flag_LIB0) {
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        if (aurostd::FileExist(directory_LIB + "/CONTCAR.static" + XHOST.vext[iext]) && aurostd::FileExist(directory_LIB + "/OUTCAR.static" + XHOST.vext[iext])) {
          const vector<string> vcarkeys{"CONTCAR", "OUTCAR", "KPOINTS"};
          for (const auto& key : vcarkeys) {
            filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".static").append(XHOST.vext[iext]));
            filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax").append(XHOST.vext[iext]));
            cout << MESSAGE << " WARNING - PATCHING " << key << ".relax with " << key << ".static " << filecar_map[key].lib << endl;
            if (aurostd::FileExist(filecar_map[key].lib)) {
              aurostd::CopyFile(filecar_map[key].lib, filecar_map[key].raw);
              aurostd::DecompressFile(filecar_map[key].raw);
              vfiles.emplace_back(key + ".relax");
              if (key == "OUTCAR") {
                FileName_OUTCAR_relax = filecar_map[key].lib;
              }
            }
          }
        }
      }
      vector<string> vcarkeys{"CONTCAR", "OUTCAR", "KPOINTS", "INCAR"};
      for (size_t i = 0; i < vcarkeys.size(); i++) {
        const string& key = vcarkeys[i];
        if (aurostd::FileExist(filecar_map[key].raw)) {
          messagestream << MESSAGE << " [" << i + 1 << "] - file not prepared " << filecar_map[key].lib;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _FILE_ERROR_);
        }
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [5]" << endl;
    }
    if (flag_SG1 || flag_SG2) { // CONTCAR.relax1
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        const vector<string> vcarkeys{"CONTCAR", "OUTCAR"};
        for (const auto& key : vcarkeys) {
          filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".relax1").append(XHOST.vext[iext]));
          filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax1").append(XHOST.vext[iext]));
          if (aurostd::FileExist(filecar_map[key].lib)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB, key + ".relax1", directory_RAW, key + ".relax1", vfiles, MESSAGE);
          }
        }
        if (!aurostd::FileExist(directory_RAW + "/CONTCAR.relax1") && aurostd::FileExist(directory_LIB + "/CONTCAR.static" + XHOST.vext[iext]) && !aurostd::FileExist(directory_RAW + "/OUTCAR.relax1") &&
            aurostd::FileExist(directory_LIB + "/OUTCAR.static" + XHOST.vext[iext])) {
          // loop cars again because the && in above if-check, combining into single loop would be equivalent to using || for the two cars
          for (const auto& key : vcarkeys) {
            filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".static").append(XHOST.vext[iext]));
            filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax1").append(XHOST.vext[iext]));
            cout << MESSAGE << " WARNING - PATCHING " << key << ".relax1 with " << key << ".static " << filecar_map[key].lib << endl;
            if (aurostd::FileExist(filecar_map[key].lib)) {
              aurostd::CopyFile(filecar_map[key].lib, filecar_map[key].raw);
              aurostd::DecompressFile(filecar_map[key].raw);
              vfiles.emplace_back(key + ".relax1");
            }
          }
        }
        if (aurostd::FileExist(filecar_map["CONTCAR"].raw)) {
          messagestream << MESSAGE << " [4] - file not prepared " << filecar_map["CONTCAR"].lib;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _FILE_ERROR_);
        }
        if (aurostd::FileExist(filecar_map["OUTCAR"].raw)) {
          messagestream << MESSAGE << " [5] - file not prepared " << filecar_map["OUTCAR"].lib;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _FILE_ERROR_);
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [6]" << endl;
    }
    if (flag_DATA_RELAX_ || flag_EDATA_RELAX_ || true) { // OSZICAR.relax.EXT
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) {
        // SKIP uncompressed
        // trying to get an OUTCAR
        if (aurostd::FileExist(directory_LIB + "/OUTCAR.relax1" + XHOST.vext[iext])) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB, "OUTCAR.relax1", directory_RAW, "OUTCAR.relax1", vfiles, MESSAGE); // _AFLOWIN_
        }
        if (aurostd::FileExist(directory_LIB + "/CONTCAR.relax1" + XHOST.vext[iext])) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB, "CONTCAR.relax1", directory_RAW, "CONTCAR.relax1", vfiles, MESSAGE); // _AFLOWIN_
        }
        const vector<string> vcarkeys{"OUTCAR", "CONTCAR", "KPOINTS", "INCAR"};
        for (const auto& key : vcarkeys) {
          filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".relax2").append(XHOST.vext[iext]));
          filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax").append(XHOST.vext[iext]));
        }
        if (aurostd::FileExist(filecar_map["OUTCAR"].lib) && aurostd::FileExist(filecar_map["CONTCAR"].lib) && aurostd::FileExist(filecar_map["KPOINTS"].lib) && aurostd::FileExist(filecar_map["INCAR"].lib)) {
          const vector<string> vcarkeys2{"OUTCAR", "CONTCAR", "KPOINTS"};
          for (const auto& key2 : vcarkeys2) {
            aurostd::CopyFile(filecar_map[key2].lib, filecar_map[key2].raw);
            aurostd::DecompressFile(filecar_map[key2].raw);
            vfiles.emplace_back(key2 + ".relax");
            if (key2 == "OUTCAR") {
              FileName_OUTCAR_relax = filecar_map[key2].lib;
            }
          }
          // nocopy INCAR
        } else {
          for (const auto& key : vcarkeys) {
            filecar_map[key].lib = aurostd::CleanFileName(string(directory_LIB).append("/").append(key).append(".relax1").append(XHOST.vext[iext]));
            filecar_map[key].raw = aurostd::CleanFileName(string(directory_RAW).append("/").append(key).append(".relax").append(XHOST.vext[iext]));
          }
          if (aurostd::FileExist(filecar_map["OUTCAR"].lib) && aurostd::FileExist(filecar_map["CONTCAR"].lib) && aurostd::FileExist(filecar_map["KPOINTS"].lib) && aurostd::FileExist(filecar_map["INCAR"].lib)) {
            const vector<string> vcarkeys2{"OUTCAR", "CONTCAR", "KPOINTS"};
            for (const auto& key2 : vcarkeys2) {
              aurostd::CopyFile(filecar_map[key2].lib, filecar_map[key2].raw);
              aurostd::DecompressFile(filecar_map[key2].raw);
              vfiles.emplace_back(key2 + ".relax");
              if (key2 == "OUTCAR") {
                FileName_OUTCAR_relax = filecar_map[key2].lib;
              }
            }
            // nocopy INCAR
          }
        }
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [6.9]" << endl;
    }
    // get code
    data.code = "nan";
    vector<string> vlines;
    aurostd::file2vectorstring(directory_RAW + "/OUTCAR.relax", vlines);
    for (auto iter = vlines.rbegin(); iter != vlines.rend(); ++iter) {
      if (aurostd::substring2bool(*iter, "vasp")) {
        aurostd::string2tokens(*iter, tokens, " ");
        break;
      }
    }
    if (tokens.size() > 1) {
      data.code = tokens[0];
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " CODE = " << data.code << endl;
    }

    // create structures
    // CO20171021 - and write xstr_json
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [7]" << endl;
    }

    if (str_orig.num_each_type.empty() && aurostd::FileExist(directory_RAW + "/POSCAR.orig")) {
      const xstructure _str_orig(directory_RAW + "/POSCAR.orig", IOVASP_AUTO);
      str_orig = _str_orig;
      str_orig.SetSpecies(deq_species); // DX20190620 - add species to xstructure
      str_orig.ReScale(1.0);
    } // CO20171025

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [7.1]" << endl;
    }

    if (str_relax.num_each_type.empty() && aurostd::FileExist(directory_RAW + "/CONTCAR.relax")) {
      const xstructure _str_relax(directory_RAW + "/CONTCAR.relax", IOVASP_AUTO);
      str_relax = _str_relax;
      str_relax.SetSpecies(deq_species); // DX20190620 - add species to xstructure
      str_relax.ReScale(1.0);
      xstructure2json(str_relax).saveFile(directory_RAW + "/" + data.system_name + "_structure_relax.json");
    } // CO20171025

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [7.2]" << endl;
    }

    if (str_relax.num_each_type.empty() && aurostd::FileExist(directory_RAW + "/CONTCAR.static")) {
      const xstructure _str_relax(directory_RAW + "/CONTCAR.static", IOVASP_AUTO);
      str_relax = _str_relax;
      str_relax.SetSpecies(deq_species); // DX20190620 - add species to xstructure
      str_relax.ReScale(1.0);
      xstructure2json(str_relax).saveFile(directory_RAW + "/" + data.system_name + "_structure_relax.json");
    } // CO20171025

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [7.3]" << endl;
    }

    if (aurostd::FileExist(directory_RAW + "/CONTCAR.relax1")) {
      const xstructure _str_relax1(directory_RAW + "/CONTCAR.relax1", IOVASP_AUTO);
      str_relax1 = _str_relax1;
      str_relax1.SetSpecies(deq_species); // DX20190620 - add species to xstructure
      str_relax1.ReScale(1.0);
      xstructure2json(str_relax1).saveFile(directory_RAW + "/" + data.system_name + "_structure_relax1.json");
    } // CO20171025
    // do the extractions

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [8]" << endl;
    }

    // LOAD STRUCTURES
    data.nspecies = str_relax.num_each_type.size();
    data.natoms = str_relax.atoms.size();
    data.volume_cell = str_relax.GetVolume();
    data.volume_atom = str_relax.GetVolume() / (double) data.natoms;
    data_abcabc = Getabc_angles(str_relax.lattice, DEGREES);
    data.vgeometry = aurostd::xvector2vector(data_abcabc);

    // DX20190124 - add original crystal info - START
    //  LOAD ORIGINAL STRUCTURE
    data.natoms_orig = str_orig.atoms.size();
    data.volume_cell_orig = str_orig.GetVolume();
    data.volume_atom_orig = str_orig.GetVolume() / (double) data.natoms_orig;
    data_abcabc = Getabc_angles(str_orig.lattice, DEGREES);
    data.vgeometry_orig = aurostd::xvector2vector(data_abcabc);

    // DX20190124 - add original crystal info - END

    // CO, get fpos now, cpos comes from outcar later (either way works)
    vector<string> fpos_strings;
    vector<string> fpos_strings_combined;
    data.vpositions_fractional.clear();
    for (size_t i = 0; i < str_relax.atoms.size(); i++) {
      data.vpositions_fractional.emplace_back(str_relax.atoms[i].fpos);
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " POSITIONS_FRACTIONAL = " << data.vpositions_fractional[i][1] << "," << data.vpositions_fractional[i][2] << "," << data.vpositions_fractional[i][3] << "," << endl;
      }
      // prepare for string variant
      for (uint j = 1; j < (uint) str_relax.atoms[i].fpos.rows + 1; j++) {
        fpos_strings.emplace_back(aurostd::utype2string(str_relax.atoms[i].fpos[j], 8));
      }
      fpos_strings_combined.emplace_back(aurostd::joinWDelimiter(fpos_strings, ","));
      fpos_strings.clear();
    }
    data.positions_fractional = aurostd::joinWDelimiter(fpos_strings_combined, ";");

    data.geometry = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vgeometry, _AFLOWLIB_DATA_GEOMETRY_PREC_), ",");
    data.geometry_orig = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vgeometry_orig, _AFLOWLIB_DATA_GEOMETRY_PREC_), ",");
    // DX20190124 - add original crystal info - END

    data.vstoichiometry = aurostd::deque2vector(str_relax.stoich_each_type);
    // CO20200624 START - mimic stoich from PrintData1(): aflow_pflow_print.cpp, this is really obsolete
    stringstream stoich_ss;
    stoich_ss.precision(4);
    for (size_t it = 0; it < str_relax.stoich_each_type.size(); it++) {
      stoich_ss << setw(8) << str_relax.stoich_each_type[it] << " ";
    }
    data.stoich = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(stoich_ss.str());
    // CO20200624 END - mimic stoich from PrintData1(): aflow_pflow_print.cpp, this is really obsolete
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " NSPECIES = " << data.nspecies << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " NATOMS = " << data.natoms << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " VOLUME (A^3) = " << data.volume_cell << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " VOLUME_ATOM (A^3) = " << data.volume_atom << "   " << directory_LIB << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " GEOMETRY (A,A,A,deg,deg,deg) = " << data.geometry << endl;
    }

    // LOAD ENERGY DATA1
    if (flag_ENERGY1) {
      xOUT.GetPropertiesFile(directory_RAW + "/OUTCAR.relax1", data.natoms, true);
      //   ExtractDataOSZICAR(directory_RAW+"/OSZICAR.relax1",data.natoms,data1_dE,data1_dEN,data1_mag,data1_mag_atom);
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENERGY1 total E0 (eV) = " << (data1_energy_cell = xOUT.energy_cell) << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENERGY1 per atom E0/N (eV) = " << (data1_energy_atom = xOUT.energy_atom) << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENTHALPY1 total E0 (eV) = " << xOUT.enthalpy_cell << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENTHALPY1 per atom E0/N (eV) = " << xOUT.enthalpy_atom << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " SPIN1 mag (\\mu) = " << xOUT.mag_cell << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " SPIN1 per atom mag/N (\\mu) = " << xOUT.mag_atom << endl;
      }
    }
    // LOAD ENERGY DATA
    xOUT.GetPropertiesFile(directory_RAW + "/OUTCAR.relax", data.natoms, true);
    kpoints.GetPropertiesFile(directory_RAW + "/KPOINTS.relax", true);

    data.pressure = xOUT.pressure;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PRESSURE (kB) = " << data.pressure << endl;
    }
    data.vstress_tensor.clear();
    data.vstress_tensor.emplace_back(xOUT.stress(1, 1));
    data.vstress_tensor.emplace_back(xOUT.stress(1, 2));
    data.vstress_tensor.emplace_back(xOUT.stress(1, 3));
    data.vstress_tensor.emplace_back(xOUT.stress(2, 1));
    data.vstress_tensor.emplace_back(xOUT.stress(2, 2));
    data.vstress_tensor.emplace_back(xOUT.stress(2, 3));
    data.vstress_tensor.emplace_back(xOUT.stress(3, 1));
    data.vstress_tensor.emplace_back(xOUT.stress(3, 2));
    data.vstress_tensor.emplace_back(xOUT.stress(3, 3));
    data.stress_tensor = aurostd::vecDouble2String(data.vstress_tensor, 7);
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " STRESS_TENSOR (kB) = " << data.stress_tensor << endl;
    }
    data.pressure_residual = xOUT.pressure_residual;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PRESSURE_RESIDUAL (kB) = " << data.pressure_residual << endl;
    }
    data.Pulay_stress = xOUT.Pulay_stress;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PULAY_STRESS (kB) = " << data.Pulay_stress << endl;
    }
    data.energy_cell = xOUT.energy_cell;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ENERGY total E0 (eV) = " << data.energy_cell << endl;
    }
    data.energy_atom = xOUT.energy_atom;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ENERGY per atom E0/N (eV) = " << data.energy_atom << endl;
    }
    data.enthalpy_cell = xOUT.enthalpy_cell;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ENTHALPY total E0 (eV) = " << data.enthalpy_cell << endl;
    }
    data.enthalpy_atom = xOUT.enthalpy_atom;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ENTHALPY per atom E0/N (eV) = " << data.enthalpy_atom << "   " << directory_LIB << endl;
    }
    data.eentropy_cell = xOUT.eentropy_cell;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " E-ENTROPY total E0 (eV) = " << data.eentropy_cell << endl;
    }
    data.eentropy_atom = xOUT.eentropy_atom;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " E-ENTROPY per atom E0/N (eV) = " << data.eentropy_atom << endl;
    }
    data.PV_cell = xOUT.PV_cell;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PV total E0 (eV) = " << data.PV_cell << endl;
    }
    data.PV_atom = xOUT.PV_atom;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PV per atom E0/N (eV) = " << data.PV_atom << endl;
    }
    data.spin_cell = xOUT.mag_cell;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " SPIN mag (\\mu) = " << data.spin_cell << endl;
    }
    data.spin_atom = xOUT.mag_atom;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " SPIN per atom mag/N (\\mu) = " << data.spin_atom << "   " << directory_LIB << endl;
    }
    // CO20180130 START
    // moving FROM magnetic loop so we keep spin/cell, spin/atom, and spinD all together
    data.spinD = aurostd::vecDouble2String(xOUT.vmag, 5);
    data.vspinD.clear();
    if (!xOUT.vmag.empty()) {
      data.spinD = aurostd::vecDouble2String(xOUT.vmag, 5);
      data.vspinD = xOUT.vmag;
    } else {
      const size_t s = static_cast<size_t>(xOUT.natoms);
      data.spinD = aurostd::vecDouble2String(vector<double>(s, 0.0));
      data.vspinD.assign(s, 0.0);
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " SPIND (\\mu) = " << data.spinD << "   " << directory_LIB << endl;
    }
    // CO20180130 STOP
    data.energy_cutoff = xOUT.ENCUT;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ENERGY_CUTOFF (eV) = " << data.energy_cutoff << endl;
    }
    data.delta_electronic_energy_convergence = xOUT.total_energy_change;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " DELTA_ELECTRONIC_ENERGY_CONVERGENCE (eV) = " << data.delta_electronic_energy_convergence << endl; // CORMAC
    }
    data.delta_electronic_energy_threshold = xOUT.EDIFF;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " DELTA_ELECTRONIC_ENERGY_THRESHOLD (eV) = " << data.delta_electronic_energy_threshold << endl; // CORMAC
    }
    data.nkpoints = kpoints.nkpoints;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " NKPOINTS (from KPOINTS) = " << data.nkpoints << endl;
    }
    data.kpoints_nnn_relax = kpoints.nnn_kpoints;
    data.kpoints = aurostd::joinWDelimiter(kpoints.nnn_kpoints, ",");
    data.nkpoints_irreducible = xOUT.nkpoints_irreducible;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " NKPOINTS_IRREDUCIBLE (from OUTCAR) = " << data.nkpoints_irreducible << endl;
    }
    data.kppra = data.natoms * kpoints.nkpoints;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " NKPPRA = " << data.kppra << endl;
    }
    data.vforces = xOUT.vforces;
    data.vpositions_cartesian = xOUT.vpositions_cartesian;
    int precfp = _DOUBLE_PRECISION_; // DX20190320 - changed from uint to int, otherwise breaks
    vector<string> tmp;
    auto xvec_joiner = [&](const string& assignment_message) {
      // lamba parameterized by the verbose message it sends
      return [&](const xvector<double>& pos) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " " << assignment_message << " = " << aurostd::xvecDouble2String(pos) << endl;
        }
        return aurostd::xvecDouble2String(pos, precfp);
      };
    };
    std::transform(data.vpositions_cartesian.begin(), data.vpositions_cartesian.end(), std::back_inserter(tmp), xvec_joiner("POSITIONS_CARTESIAN"));
    data.positions_cartesian = aurostd::joinWDelimiter(tmp, ';');
    tmp.clear();
    std::transform(data.forces.begin(), data.forces.end(), std::back_inserter(tmp), xvec_joiner("FORCES(eV/Angst)"));
    data.forces = aurostd::joinWDelimiter(tmp, ';');
    tmp.clear();

    // LOAD SPECIES
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " " << data.nspecies << " " << str_relax.species.size() << " " << str_relax.species_pp.size() << " " << str_relax.species_pp_type.size() << " "
           << str_relax.species_pp_version.size() << " " << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_volume.size() << " " << str_relax.species_mass.size() << endl;
    }

    str_relax.species.clear();
    str_relax.species_pp.clear();
    str_relax.species_pp_type.clear();
    str_relax.species_pp_version.clear();
    str_relax.species_pp_ZVAL.clear();
    str_relax.species_pp_vLDAU.clear();
    str_relax.species_volume.clear();
    str_relax.species_mass.clear();

    // try OUTCARs
    data.vspecies.clear();
    data.vspecies_pp.clear();
    data.vspecies_pp_version.clear();
    data.vspecies_pp_ZVAL.clear();
    data.vspecies_pp_AUID.clear();
    data.METAGGA = "";

    vector<string> stries{"relax1", "relax2", "relax3", "static", "bands", "dielectric"};
    tmp_key = "OUTCAR";
    for (size_t i = 0; i < stries.size() && str_relax.species.empty(); i++) { // ONLY OUTCAR as POTCARS are obsolete
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        filecar_map[tmp_key].lib = string(directory_LIB).append("/").append(tmp_key).append(".").append(stries[i]).append(XHOST.vext[iext]);
        if (aurostd::FileExist(filecar_map[tmp_key].lib)) {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " fileE_LIB=" << filecar_map[tmp_key].lib << endl;
          }
          xOUT.GetPropertiesFile(filecar_map[tmp_key].lib);
          // DEBUG	for(size_t i=0;i<xOUT.species.size();i++) cerr << XPID << "xOUT.species.at(i)=" << xOUT.species.at(i) << endl;
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " xOUT.species.size()=" << xOUT.species.size() << endl;
          }
          str_relax.species.assign(xOUT.species.begin(), xOUT.species.end());
          str_relax.species_pp.assign(xOUT.species_pp.begin(), xOUT.species_pp.end());
          str_relax.species_pp_type.assign(xOUT.species_pp_type.begin(), xOUT.species_pp_type.end());
          str_relax.species_pp_version.assign(xOUT.species_pp_version.begin(), xOUT.species_pp_version.end());
          str_relax.species_pp_ZVAL.assign(xOUT.vZVAL.begin(), xOUT.vZVAL.end());
          data.vspecies.insert(data.vspecies.end(), xOUT.species.begin(), xOUT.species.end());
          data.vspecies_pp.insert(data.vspecies_pp.end(), xOUT.species_pp.begin(), xOUT.species_pp.end());
          data.vspecies_pp_version.insert(data.vspecies_pp_version.end(), xOUT.species_pp_version.begin(), xOUT.species_pp_version.end());
          data.vspecies_pp_ZVAL.insert(data.vspecies_pp_ZVAL.end(), xOUT.vZVAL.begin(), xOUT.vZVAL.end());
          data.vspecies_pp_AUID.assign(xOUT.species_pp_AUID.begin(), xOUT.species_pp_AUID.end());
          data.dft_type = xOUT.pp_type;
          // CO20210213 - check types are all the same, if not issue warning/error (mixing is not advisable)
          for (size_t j = 0; j < xOUT.species_pp_type.size(); j++) {
            if (xOUT.species_pp_type[j] != xOUT.pp_type) {
              pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Mismatch in species_pp_types (" + xOUT.species_pp_type[j] + " vs. " + xOUT.pp_type + ")", _LOGGER_WARNING_);
            }
          }
          data.vdft_type = {xOUT.pp_type}; // CO, this is technically a vector (RESTAPI paper)
          str_relax.species_pp_vLDAU = xOUT.species_pp_vLDAU;
          data.ldau_TLUJ = xOUT.string_LDAU;
          //[CO+ME20210713 - keep legacy behavior, only print when non-zero]if(data.ldau_TLUJ.empty()){data.ldau_TLUJ=aurostd::utype2string(0);} //CO20210713 - no +U
          data.METAGGA = xOUT.METAGGA;

          // ME20190124 BEGIN - Store LDAU information individually
          //  Note that the vector here has the species in the columns, not the
          //  rows because this is closer to the format in the out and json files.
          if (!xOUT.species_pp_vLDAU.empty()) {
            data.vLDAU.resize(xOUT.species_pp_vLDAU[0].size());
          } // CO20200731
          else { // CO20210713 - set ldau_type=0
            data.vLDAU.resize(4);
            for (size_t j = 0; j < xOUT.species.size(); j++) {
              data.vLDAU[0].push_back(0);
            }
          }
          for (size_t k = 0; k < xOUT.species_pp_vLDAU.size(); k++) {
            for (size_t j = 0; j < xOUT.species_pp_vLDAU[k].size(); j++) { // CO20200731 - this WILL break if xOUT.species_pp_vLDAU[i].size()==0 //4
              data.vLDAU[j].emplace_back(xOUT.species_pp_vLDAU[k][j]);
            }
          }
          // ME20190124 END
          if (AFLOWLIB_VERBOSE && !data.ldau_TLUJ.empty()) {
            cout << MESSAGE << " LDAU_string=" << data.ldau_TLUJ << endl;
          }
        }
      }
    }

    if (str_relax.species.empty()) {
      messagestream << "OUTCAR/POTCAR not FOUND in " << directory_LIB;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _FILE_NOT_FOUND_);
    }

    for (size_t j = 0; j < str_relax.species.size(); j++) {
      str_relax.species_volume.emplace_back(GetAtomVolume(str_relax.species[j]));
      str_relax.species_mass.emplace_back(GetAtomMass(str_relax.species[j]));
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " " << data.nspecies << " " << str_relax.species.size() << " " << str_relax.species_pp.size() << " " << str_relax.species_pp_type.size() << " " << str_relax.species_pp_version.size()
           << " " << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_pp_vLDAU.size() << " " << str_relax.species_volume.size() << " " << str_relax.species_mass.size() << endl;
    }
    vector<std::pair<size_t, string>> size_checks{
        {           str_relax.species.size(),            "str_relax.species.size()"},
        {        str_relax.species_pp.size(),         "str_relax.species_pp.size()"},
        {   str_relax.species_pp_type.size(),    "str_relax.species_pp_type.size()"},
        {str_relax.species_pp_version.size(), "str_relax.species_pp_version.size()"},
        {   str_relax.species_pp_ZVAL.size(),    "str_relax.species_pp_ZVAL.size()"},
        {       data.vspecies_pp_AUID.size(),        "data.vspecies_pp_AUID.size()"},
        {    str_relax.species_volume.size(),     "str_relax.species_volume.size()"},
        {      str_relax.species_mass.size(),       "str_relax.species_mass.size()"}
    };
    for (size_t i = 0; i < size_checks.size(); i++) {
      if (data.nspecies != size_checks[i].first) {
        messagestream << MESSAGE << " [" << i + 1 << "] - data.nspecies[" << data.nspecies << "]!=" << size_checks[i].second << "]";
        if (i == 0) {
          messagestream << endl << str_relax;
        }
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _INDEX_MISMATCH_);
      }
    }

    data.compound.clear();
    data.density = 0.0;
    data.valence_cell_iupac = 0.0;
    data.valence_cell_std = 0.0;

    data.composition = aurostd::joinWDelimiter(str_relax.num_each_type, ',');
    data.species = aurostd::joinWDelimiter(str_relax.species, ',');
    data.species_pp = aurostd::joinWDelimiter(str_relax.species_pp, ',');
    data.species_pp_version = aurostd::joinWDelimiter(str_relax.species_pp_version, ',');
    data.species_pp_AUID = aurostd::joinWDelimiter(data.vspecies_pp_AUID, ',');
    data.species_pp_ZVAL = aurostd::vecDouble2String(str_relax.species_pp_ZVAL);
    data.stoichiometry = aurostd::vecDouble2String(data.vstoichiometry, _AFLOWLIB_STOICH_PRECISION_);
    data.vcomposition = vector<double>(str_relax.num_each_type.begin(), str_relax.num_each_type.end());

    for (uint i = 0; i < data.nspecies; i++) {
      data.compound += str_relax.species[i];
      data.compound += aurostd::utype2string(str_relax.num_each_type[i]);
      data.density += (double) str_relax.num_each_type[i] * GetAtomMass(str_relax.species.at(i));
      data.valence_cell_iupac += str_relax.num_each_type[i] * GetAtomValenceIupac(str_relax.species[i]);
      data.valence_cell_std += str_relax.num_each_type[i] * GetAtomValenceStd(str_relax.species[i]);
    }

    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " VALENCE_IUPAC = " << data.valence_cell_iupac << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " VALENCE_STD = " << data.valence_cell_std << endl;
    }

    // density
    data.density /= data.volume_cell;
    data.density *= 1000.0; // grams instead of kilos
    data.density *= 1e8 * 1e8 * 1e8; // cm^3 instead of A^3
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " DENSITY (grams/cm^3) = " << data.density << endl;
    }

    // DX20190124 - original density info - START
    data.density_orig = 0.0;
    for (uint i = 0; i < (uint) data.nspecies; i++) {
      data.density_orig += (double) str_orig.num_each_type.at(i) * GetAtomMass(str_orig.species.at(i));
    }
    // density
    data.density_orig /= data.volume_cell_orig;
    data.density_orig *= 1000.0; // grams instead of kilos
    data.density_orig *= 1e8 * 1e8 * 1e8; // cm^3 instead of A^3
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " DENSITY_ORIG (grams/cm^3) = " << data.density_orig << endl;
    }
    // DX20190124 - original density info - START

    // scintillation_attenuation_length
    data.scintillation_attenuation_length = 0.0;
    data.scintillation_attenuation_length = GetCompoundAttenuationLength(str_relax.species, str_relax.num_each_type, data.density);
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " SCINTILLATION_ATTENUATION_LENGTH (cm) = " << data.scintillation_attenuation_length << endl;
    }

    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PSEUDOPOTENTIAL dft_type=" << data.dft_type << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_version = " << data.species_pp_version << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_ZVAL = " << data.species_pp_ZVAL << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_AUID = " << data.species_pp_AUID << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " PSEUDOPOTENTIAL METAGGA = [" << data.METAGGA << "]" << endl;
    }

    const bool FORMATION_CALC = LIB2RAW_Calculate_FormationEnthalpy(data, str_relax, MESSAGE);

    if (FORMATION_CALC == true) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENTHALPY FORMATION total E0 (eV) = " << data.enthalpy_formation_cell << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENTHALPY FORMATION per atom E0/N (eV) = " << data.enthalpy_formation_atom << "   " << directory_LIB << endl;
      }
      // CO20200624 START - CCE
      if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_300K_cell != AUROSTD_NAN) {
        cout << MESSAGE << " ENTHALPY FORMATION CCE total E(300K) (eV) = " << data.enthalpy_formation_cce_300K_cell << endl;
      }
      if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_300K_atom != AUROSTD_NAN) {
        cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(300K)/N (eV) = " << data.enthalpy_formation_cce_300K_atom << "   " << directory_LIB << endl;
      }
      if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_0K_cell != AUROSTD_NAN) {
        cout << MESSAGE << " ENTHALPY FORMATION CCE total E(0K) (eV) = " << data.enthalpy_formation_cce_0K_cell << endl;
      }
      if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_0K_atom != AUROSTD_NAN) {
        cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(0K)/N (eV) = " << data.enthalpy_formation_cce_0K_atom << "   " << directory_LIB << endl;
      }
      // CO20200624 END - CCE
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENTROPIC_TEMPERATURE (eV) = " << data.entropic_temperature * KBOLTZEV << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " ENTROPIC_TEMPERATURE (K) = " << data.entropic_temperature << "   " << directory_LIB << endl;
      }
    }
    // DONE WITH THERMO

    // do the TIMING
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [9]" << endl;
    }
    data.calculation_cores = 1;
    data.calculation_time = 0.0;
    data.calculation_memory = 0.0;
    if (flag_TIMING) { // OUTCAR.relax.EXT
      for (uint i = 1; i <= relax_max; i++) {
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          filecar_map["OUTCAR"].lib = directory_LIB + "/OUTCAR.relax" + aurostd::utype2string(i) + XHOST.vext[iext];
          if (aurostd::FileExist(filecar_map["OUTCAR"].lib)) {
            xOUTCAR outcar_tmp;
            outcar_tmp.GetPropertiesFile(filecar_map["OUTCAR"].lib); // OK
            data.calculation_cores = aurostd::max(data.calculation_cores, (int) outcar_tmp.calculation_cores);
            data.calculation_time += outcar_tmp.calculation_time; // will multiply after
            data.calculation_memory = aurostd::max(data.calculation_memory, outcar_tmp.calculation_memory);
          }
        }
      }
      const vector<string> post_relax = {"static", "bands", "dielectric"};
      for (const string& run_type : post_relax) {
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
          filecar_map["OUTCAR"].lib = directory_LIB + "/OUTCAR." + run_type + XHOST.vext[iext];
          if (aurostd::FileExist(filecar_map["OUTCAR"].lib)) {
            xOUTCAR outcar_tmp;
            outcar_tmp.GetPropertiesFile(filecar_map["OUTCAR"].lib); // OK
            data.calculation_cores = aurostd::max(data.calculation_cores, (int) outcar_tmp.calculation_cores);
            data.calculation_time += outcar_tmp.calculation_time; // will multiply after
            data.calculation_memory = aurostd::max(data.calculation_memory, outcar_tmp.calculation_memory);
          }
        }
      }

      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " CALCULATION time (sec) = " << data.calculation_time << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " CALCULATION mem (MB) = " << data.calculation_memory << endl;
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " CALCULATION cores = " << data.calculation_cores << endl;
      }
    }
    // do the SG
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [10]" << endl;
    }
    if (flag_SG1 || flag_SG2) { // POSCAR.orig CONTCAR.relax1 CONTCAR.relax
#ifdef USE_PLATON_SG
      string space_group_calculator_sg1_cmd = " | aflow --platonSG=3.0,0.5,0.5,0.5";
      string space_group_calculator_sg2_cmd = " | aflow --platonSG=1.5,0.25,0.25,0.25";

      string ilattice_cmd = " | aflow --ilattice 2.0";
      if (flag_SG1) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg1_cmd << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg1_cmd << endl;
        ssfile << directory_RAW << endl;
        // I run outside so a segfault+core does not block the aflow
        // INSIDE no more inside to avoid BOMBING
        // str_orig.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_pre=str_orig.spacegroup;     //  aflow --platonSG 3.0 0.5 0.5 0.5
        // str_relax1.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_mid=str_relax1.spacegroup; //  aflow --platonSG 3.0 0.5 0.5 0.5
        // str_relax.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_post=str_relax.spacegroup;  //  aflow --platonSG 3.0 0.5 0.5 0.5
        // OUTSIDE
        // sg1_pre
        data_sg1_pre = aurostd::execute2string("cat \"" + directory_RAW + "/POSCAR.orig" + "\"" + space_group_calculator_sg1_cmd);
        aurostd::string2tokens(data_sg1_pre, tokens, "#");
        if (tokens.size() != 2) {
          data_sg1_pre = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg1_pre = NOSG;
          } else {
            aurostd::StringSubst(data_sg1_pre, "\n", "");
          }
        }
        if (data_sg1_pre == NOSG) {
          if (AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_pre" << endl;
          data_sg1_pre = aurostd::execute2string("cat \"" + directory_RAW + "/POSCAR.orig" + "\"" + ilattice_cmd + space_group_calculator_sg1_cmd);
        }
        aurostd::string2tokens(data_sg1_pre, tokens, "#");
        if (tokens.size() != 2) {
          data_sg1_pre = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg1_pre = NOSG;
          } else {
            aurostd::StringSubst(data_sg1_pre, "\n", "");
          }
        }
        if (aurostd::substring2bool(data_sg1_pre, "SymbolnotKnown")) data_sg1_pre = NOSG; // give up
        // sg1_mid
        data_sg1_mid = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax1" + "\"" + space_group_calculator_sg1_cmd);
        aurostd::string2tokens(data_sg1_mid, tokens, "#");
        if (tokens.size() != 2) {
          data_sg1_mid = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg1_mid = NOSG;
          } else {
            aurostd::StringSubst(data_sg1_mid, "\n", "");
          }
        }
        if (data_sg1_mid == NOSG) {
          if (AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_mid" << endl;
          data_sg1_mid = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax1" + "\"" + ilattice_cmd + space_group_calculator_sg1_cmd);
        }
        aurostd::string2tokens(data_sg1_mid, tokens, "#");
        if (tokens.size() != 2) {
          data_sg1_mid = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg1_mid = NOSG;
          } else {
            aurostd::StringSubst(data_sg1_mid, "\n", "");
          }
        }
        if (aurostd::substring2bool(data_sg1_mid, "SymbolnotKnown")) data_sg1_mid = NOSG; // give up
        // sg1_mid
        data_sg1_post = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax" + "\"" + space_group_calculator_sg1_cmd);
        aurostd::string2tokens(data_sg1_post, tokens, "#");
        if (tokens.size() != 2) {
          data_sg1_post = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg1_post = NOSG;
          } else {
            aurostd::StringSubst(data_sg1_post, "\n", "");
          }
        }
        if (data_sg1_post == NOSG) {
          if (AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_post" << endl;
          data_sg1_post = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax" + "\"" + ilattice_cmd + space_group_calculator_sg1_cmd);
        }
        aurostd::string2tokens(data_sg1_post, tokens, "#");
        if (tokens.size() != 2) {
          data_sg1_post = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg1_post = NOSG;
          } else {
            aurostd::StringSubst(data_sg1_post, "\n", "");
          }
        }
        if (aurostd::substring2bool(data_sg1_post, "SymbolnotKnown")) data_sg1_post = NOSG; // give up
        // DONE
        ssfile << "PRE  " << data_sg1_pre << endl;
        ssfile << "MID  " << data_sg1_mid << endl;
        ssfile << "POST " << data_sg1_post << endl;
        //      aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP1_OUT);
        data.sg = data_sg1_pre + "," + data_sg1_mid + "," + data_sg1_post;
        data.vsg.clear();
        data.vsg.emplace_back(data_sg1_pre);
        data.vsg.emplace_back(data_sg1_mid);
        data.vsg.emplace_back(data_sg1_post); // CO20171202
        if (AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP1 = " << data.sg << endl;
      }
      if (flag_SG2) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg2_cmd << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg2_cmd << endl;
        ssfile << directory_RAW << endl;
        // I run outside so a segfault+core does not block the aflow
        // INSIDE
        // str_orig.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_pre=str_orig.spacegroup;     //  aflow --platonSG 1.5 0.25 0.25 0.25
        // str_relax1.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_mid=str_relax1.spacegroup; //  aflow --platonSG 1.5 0.25 0.25 0.25
        // str_relax.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_post=str_relax.spacegroup;  //  aflow --platonSG 1.5 0.25 0.25 0.25
        // OUTSIDE
        // sg2_pre
        data_sg2_pre = aurostd::execute2string("cat \"" + directory_RAW + "/POSCAR.orig" + "\"" + space_group_calculator_sg2_cmd);
        aurostd::string2tokens(data_sg2_pre, tokens, "#");
        if (tokens.size() != 2) {
          data_sg2_pre = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg2_pre = NOSG;
          } else {
            aurostd::StringSubst(data_sg2_pre, "\n", "");
          }
        }
        if (data_sg2_pre == NOSG) {
          if (AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_pre" << endl;
          data_sg2_pre = aurostd::execute2string("cat \"" + directory_RAW + "/POSCAR.orig" + "\"" + ilattice_cmd + space_group_calculator_sg2_cmd);
        }
        aurostd::string2tokens(data_sg2_pre, tokens, "#");
        if (tokens.size() != 2) {
          data_sg2_pre = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg2_pre = NOSG;
          } else {
            aurostd::StringSubst(data_sg2_pre, "\n", "");
          }
        }
        if (aurostd::substring2bool(data_sg2_pre, "SymbolnotKnown")) data_sg2_pre = NOSG; // give up
        // sg2_mid
        data_sg2_mid = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax1" + "\"" + space_group_calculator_sg2_cmd);
        aurostd::string2tokens(data_sg2_mid, tokens, "#");
        if (tokens.size() != 2) {
          data_sg2_mid = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg2_mid = NOSG;
          } else {
            aurostd::StringSubst(data_sg2_mid, "\n", "");
          }
        }
        if (data_sg2_mid == NOSG) {
          if (AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_mid" << endl;
          data_sg2_mid = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax1" + "\"" + ilattice_cmd + space_group_calculator_sg2_cmd);
        }
        aurostd::string2tokens(data_sg2_mid, tokens, "#");
        if (tokens.size() != 2) {
          data_sg2_mid = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg2_mid = NOSG;
          } else {
            aurostd::StringSubst(data_sg2_mid, "\n", "");
          }
        }
        if (aurostd::substring2bool(data_sg2_mid, "SymbolnotKnown")) data_sg2_mid = NOSG; // give up
        // sg2_mid
        data_sg2_post = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax" + "\"" + space_group_calculator_sg2_cmd);
        aurostd::string2tokens(data_sg2_post, tokens, "#");
        if (tokens.size() != 2) {
          data_sg2_post = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg2_post = NOSG;
          } else {
            aurostd::StringSubst(data_sg2_post, "\n", "");
          }
        }
        if (data_sg2_post == NOSG) {
          if (AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_post" << endl;
          data_sg2_post = aurostd::execute2string("cat \"" + directory_RAW + "/CONTCAR.relax" + "\"" + ilattice_cmd + space_group_calculator_sg2_cmd);
        }
        aurostd::string2tokens(data_sg2_post, tokens, "#");
        if (tokens.size() != 2) {
          data_sg2_post = NOSG;
        } else {
          if (aurostd::string2utype<uint>(tokens[1]) == 0) {
            data_sg2_post = NOSG;
          } else {
            aurostd::StringSubst(data_sg2_post, "\n", "");
          }
        }
        if (aurostd::substring2bool(data_sg2_post, "SymbolnotKnown")) data_sg2_post = NOSG; // give up
        // DONE
        ssfile << "PRE  " << data_sg2_pre << endl;
        ssfile << "MID  " << data_sg2_mid << endl;
        ssfile << "POST " << data_sg2_post << endl;
        // aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP2_OUT);
        data.sg2 = data_sg2_pre + "," + data_sg2_mid + "," + data_sg2_post;
        data.vsg2.clear();
        data.vsg2.emplace_back(data_sg2_pre);
        data.vsg2.emplace_back(data_sg2_mid);
        data.vsg2.emplace_back(data_sg2_post); // CO20171202
        if (AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP2 = " << data.sg2 << endl;
      }
#endif
      // DX20190319 - moved endif to end of external functions

      // DX20190319 - now call aflowSG internally - START
#ifdef USE_AFLOW_SG
      // DX+CO START

      const string space_group_calculator_sg1_function = "aflowSG(loose)";
      const string space_group_calculator_sg2_function = "aflowSG(tight)";
      bool no_scan = false;

      // pre
      const xstructure str_orig(directory_RAW + "/POSCAR.orig", IOAFLOW_AUTO);
      const double pre_default_tolerance = SYM::defaultTolerance(str_orig);
      // mid
      const xstructure str_mid(directory_RAW + "/CONTCAR.relax1", IOAFLOW_AUTO);
      const double mid_default_tolerance = SYM::defaultTolerance(str_mid);
      // post
      const xstructure str_post(directory_RAW + "/CONTCAR.relax", IOAFLOW_AUTO);
      const double post_default_tolerance = SYM::defaultTolerance(str_post);

      if (flag_SG1) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg1_function << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg1_function << endl;
        ssfile << directory_RAW << endl;

        // sg1_pre
        xstructure str_sg1_pre = str_orig;
        double pre_loose_tolerance = pre_default_tolerance * 10.0;
        const uint sgroup_sg1_pre = str_sg1_pre.SpaceGroup_ITC(pre_loose_tolerance, no_scan);
        if (sgroup_sg1_pre < 1 || sgroup_sg1_pre > 230) {
          data_sg1_pre = NOSG;
        } else {
          data_sg1_pre = GetSpaceGroupName(sgroup_sg1_pre, str_sg1_pre.directory) + " #" + aurostd::utype2string(sgroup_sg1_pre);
        }

        // sg1_mid
        xstructure str_sg1_mid = str_mid;
        double mid_loose_tolerance = mid_default_tolerance * 10.0;
        const uint sgroup_sg1_mid = str_sg1_mid.SpaceGroup_ITC(mid_loose_tolerance, no_scan);
        if (sgroup_sg1_mid < 1 || sgroup_sg1_mid > 230) {
          data_sg1_mid = NOSG;
        } else {
          data_sg1_mid = GetSpaceGroupName(sgroup_sg1_mid, str_sg1_mid.directory) + " #" + aurostd::utype2string(sgroup_sg1_mid);
        }

        // sg1_post
        xstructure str_sg1_post = str_post;
        double post_loose_tolerance = post_default_tolerance * 10.0;
        const uint sgroup_sg1_post = str_sg1_post.SpaceGroup_ITC(post_loose_tolerance, no_scan);
        if (sgroup_sg1_post < 1 || sgroup_sg1_post > 230) {
          data_sg1_post = NOSG;
        } else {
          data_sg1_post = GetSpaceGroupName(sgroup_sg1_post, str_sg1_post.directory) + " #" + aurostd::utype2string(sgroup_sg1_post);
        }

        // DONE
        ssfile << "PRE  " << data_sg1_pre << endl;
        ssfile << "MID  " << data_sg1_mid << endl;
        ssfile << "POST " << data_sg1_post << endl;
        data.sg = data_sg1_pre + "," + data_sg1_mid + "," + data_sg1_post;
        data.vsg = {data_sg1_pre, data_sg1_mid, data_sg1_post}; // CO20171202
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " SPACEGROUP1 = " << data.sg << endl;
        }
      }
      if (flag_SG2) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg2_function << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg2_function << endl;
        ssfile << directory_RAW << endl;

        // sg2_pre
        xstructure str_sg2_pre = str_orig;
        double pre_tight_tolerance = pre_default_tolerance;
        const uint sgroup_sg2_pre = str_sg2_pre.SpaceGroup_ITC(pre_tight_tolerance, no_scan);
        if (sgroup_sg2_pre < 1 || sgroup_sg2_pre > 230) {
          data_sg2_pre = NOSG;
        } else {
          data_sg2_pre = GetSpaceGroupName(sgroup_sg2_pre, str_sg2_pre.directory) + " #" + aurostd::utype2string(sgroup_sg2_pre);
        }

        // sg2_mid
        xstructure str_sg2_mid = str_mid;
        double mid_tight_tolerance = mid_default_tolerance;
        const uint sgroup_sg2_mid = str_sg2_mid.SpaceGroup_ITC(mid_tight_tolerance, no_scan);
        if (sgroup_sg2_mid < 1 || sgroup_sg2_mid > 230) {
          data_sg2_mid = NOSG;
        } else {
          data_sg2_mid = GetSpaceGroupName(sgroup_sg2_mid, str_sg2_mid.directory) + " #" + aurostd::utype2string(sgroup_sg2_mid);
        }

        // sg2_post
        xstructure str_sg2_post = str_post;
        double post_tight_tolerance = post_default_tolerance;
        const uint sgroup_sg2_post = str_sg2_post.SpaceGroup_ITC(post_tight_tolerance, no_scan);
        if (sgroup_sg2_post < 1 || sgroup_sg2_post > 230) {
          data_sg2_post = NOSG;
        } else {
          data_sg2_post = GetSpaceGroupName(sgroup_sg2_post, str_sg2_post.directory) + " #" + aurostd::utype2string(sgroup_sg2_post);
        }

        // DONE
        ssfile << "PRE  " << data_sg2_pre << endl;
        ssfile << "MID  " << data_sg2_mid << endl;
        ssfile << "POST " << data_sg2_post << endl;
        data.sg2 = data_sg2_pre + "," + data_sg2_mid + "," + data_sg2_post;
        data.vsg2 = {data_sg2_pre, data_sg2_mid, data_sg2_post}; // CO20171202
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " SPACEGROUP2 = " << data.sg2 << endl;
        }
      }

#endif
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [13]" << endl;
    }
    // check for ERRORS
    if (flag_ENERGY1) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " [flag_ENERGY1]" << endl;
      }
      if (aurostd::abs(data.energy_atom - data1_energy_atom) > ENERGY_ATOM_ERROR_meV / 1000.0) {
        if (data.energy_atom < 0.0 && data1_energy_atom > 0.0) {
          flag_ERROR = true;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [14]" << endl;
    }
    // done now WRITE aflowlib.out
    if (flag_ERROR) {
      data.error_status = "ERROR";
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [15]" << endl;
    }

    // PERFORM EDATA DATA AND CIF STEPS -------------------------------------------------------------------------
    xstructure str;
    xstructure str_sp;
    xstructure str_sc;
    vector<xstructure> vcif;

    const vector<std::tuple<char, string, string>> groups_outs_json{
        {           _PGROUP_,            DEFAULT_AFLOW_PGROUP_OUT,            DEFAULT_AFLOW_PGROUP_JSON},
        {          _PGROUPK_,           DEFAULT_AFLOW_PGROUPK_OUT,           DEFAULT_AFLOW_PGROUPK_JSON},
        {           _FGROUP_,            DEFAULT_AFLOW_FGROUP_OUT,            DEFAULT_AFLOW_FGROUP_JSON},
        {      _PGROUP_XTAL_,       DEFAULT_AFLOW_PGROUP_XTAL_OUT,       DEFAULT_AFLOW_PGROUP_XTAL_JSON},
        {     _PGROUPK_XTAL_,      DEFAULT_AFLOW_PGROUPK_XTAL_OUT,      DEFAULT_AFLOW_PGROUPK_XTAL_JSON},
        {_PGROUPK_PATTERSON_, DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT, DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON},
        {           _IATOMS_,            DEFAULT_AFLOW_IATOMS_OUT,            DEFAULT_AFLOW_IATOMS_JSON},
        {           _AGROUP_,            DEFAULT_AFLOW_AGROUP_OUT,            DEFAULT_AFLOW_AGROUP_JSON}
    };

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [16]" << endl;
    }
    // PERFORM EDATA STEP
    if (flag_EDATA_ORIG_ || flag_EDATA_RELAX_) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " EDATA start: " << directory_RAW << endl;
      }
    }
    if (flag_EDATA_ORIG_) { // ORIG
      if (!aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_EDATA_ORIG_OUT) && aurostd::FileExist(directory_RAW + "/POSCAR.orig")) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " EDATA doing orig (POSCAR.orig) text format: " << directory_RAW << endl;
        }
        str = xstructure(directory_RAW + "/POSCAR.orig", IOAFLOW_AUTO);
        str_sp.clear();
        str_sc.clear(); // DX20191220 - uppercase to lowercase clear
        // DX START
        // DX END
        stringstream sss;
        sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym = str; // CO20171027 //DX20180226 - set equal str
        aurostd::xoption vpflow_edata_orig; // DX20180823 - added xoption
        sss << pflow::PrintData(str, str_sym, str_sp, str_sc, vpflow_edata_orig, "EDATA", txt_ft, false); // 1=EDATA //CO20171025 //CO20171027 //DX20210301 - void to string output, "txt" to txt_ft
        aurostd::stringstream2file(sss, directory_RAW + "/" + DEFAULT_FILE_EDATA_ORIG_OUT);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " EDATA doing orig (POSCAR.orig) json format: " << directory_RAW << endl;
        }
        stringstream jjj; // CO20171025
        jjj << pflow::PrintData(str_sym, str_sym, str_sp, str_sc, vpflow_edata_orig, "EDATA", json_ft, true); // 1=EDATA, already_calculated!  //CO20171025 //CO20171027  //DX20210301 - void to string output, "json" to json_ft
        aurostd::stringstream2file(jjj, directory_RAW + "/" + DEFAULT_FILE_EDATA_ORIG_JSON); // CO20171025
        vcif = {str, str_sp, str_sc};
        // CO+DX START 20170713 - adding symmetry output to RAW
        string new_sym_file;
        for (const auto& [group, out, json] : groups_outs_json) {
          const string out_path = string(directory_RAW).append("/").append(out);
          const string json_path = string(directory_RAW).append("/").append(json);
          // txt variants
          KBIN_SymmetryWrite(FileMESSAGE, str_sp, aflags, group, false, messagestream, "txt");
          if (aurostd::FileExist(out_path)) {
            AddFileNameBeforeExtension(out_path, "orig", new_sym_file);
            aurostd::file2file(out_path, new_sym_file);
          }
          // json variants
          KBIN_SymmetryWrite(FileMESSAGE, str_sp, aflags, group, false, messagestream, "json");
          if (aurostd::FileExist(json_path)) {
            AddFileNameBeforeExtension(json_path, "orig", new_sym_file);
            aurostd::file2file(json_path, new_sym_file);
          }
        }
        // cout << messagestream;
        // CO+DX STOP 20170713 - adding symmetry output to RAW
        //  now extract info
        // DX20180823 - extract info from xoption - START
        if (vpflow_edata_orig.flag("EDATA::CALCULATED")) {
          if (data.Bravais_lattice_orig.empty()) { // Bravais_Lattice_orig
            data.Bravais_lattice_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_TYPE");
          }
          if (data.lattice_variation_orig.empty()) { // Bravais_Lattice_Variation_orig
            data.lattice_variation_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE");
          }
          if (data.lattice_system_orig.empty()) { // Lattice_System_orig
            data.lattice_system_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_SYSTEM");
          }
          if (data.Pearson_symbol_orig.empty()) { // Pearson_orig
            data.Pearson_symbol_orig = vpflow_edata_orig.getattachedscheme("EDATA::PEARSON_SYMBOL");
          }
          // DX20190124 - extract additional info from xoption - START
          if (data.crystal_system_orig.empty()) { // crystal_system_orig
            data.crystal_system_orig = vpflow_edata_orig.getattachedscheme("EDATA::CRYSTAL_SYSTEM");
          }
          if (data.crystal_family_orig.empty()) { // crystal_family_orig
            data.crystal_family_orig = vpflow_edata_orig.getattachedscheme("EDATA::CRYSTAL_FAMILY");
          }
          if (data.point_group_Hermann_Mauguin_orig.empty()) { // point_group_Hermann_Mauguin_orig
            data.point_group_Hermann_Mauguin_orig = vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_HERMANN_MAUGUIN");
          }
          if (data.crystal_class_orig.empty()) { // crystal_class_orig
            data.crystal_class_orig = vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_CRYSTAL_CLASS");
          }
          if (data.point_group_Schoenflies_orig.empty()) { // point_group_Schoenflies_orig
            data.point_group_Schoenflies_orig = vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_SCHOENFLIES");
          }
          if (data.point_group_orbifold_orig.empty()) { // point_group_orbifold_orig
            data.point_group_orbifold_orig = vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_ORBIFOLD");
          }
          if (data.point_group_type_orig.empty()) { // point_group_type_orig
            data.point_group_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_TYPE");
          }
          if (data.point_group_order_orig == AUROSTD_NAN) { // point_group_order_orig
            data.point_group_order_orig = vpflow_edata_orig.getattachedutype<uint>("EDATA::POINT_GROUP_ORDER");
          }
          if (data.point_group_structure_orig.empty()) { // point_group_structure_orig
            data.point_group_structure_orig = vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_STRUCTURE");
          }
          if (data.Bravais_lattice_lattice_type_orig.empty()) { // Bravais_lattice_lattice_type_orig
            data.Bravais_lattice_lattice_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE");
          }
          if (data.Bravais_lattice_lattice_variation_type_orig.empty()) { // Bravais_lattice_lattice_variation_type_orig
            data.Bravais_lattice_lattice_variation_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE");
          }
          if (data.Bravais_lattice_lattice_system_orig.empty()) { // Bravais_lattice_lattice_system_orig
            data.Bravais_lattice_lattice_system_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM");
          }
          if (data.Bravais_superlattice_lattice_type_orig.empty()) { // Bravais_superlattice_lattice_type_orig
            data.Bravais_superlattice_lattice_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE");
          }
          if (data.Bravais_superlattice_lattice_variation_type_orig.empty()) { // Bravais_superlattice_lattice_variation_type_orig
            data.Bravais_superlattice_lattice_variation_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE");
          }
          if (data.Bravais_superlattice_lattice_system_orig.empty()) { // Bravais_superlattice_lattice_system_orig
            data.Bravais_superlattice_lattice_system_orig = vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM");
          }
          if (data.Pearson_symbol_superlattice_orig.empty()) { // Pearson_symbol_superlattice_orig
            data.Pearson_symbol_superlattice_orig = vpflow_edata_orig.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE");
          }
          if (data.reciprocal_geometry_orig.empty()) { // reciprocal_geometry_orig
            data.reciprocal_geometry_orig = vpflow_edata_orig.getattachedscheme("EDATA::RECIPROCAL_LATTICE_PARAMETERS");
            vector<string> ktokens;
            aurostd::string2tokens(data.reciprocal_geometry_orig, ktokens, ",");
            for (size_t t = 0; t < ktokens.size(); t++) {
              data.vreciprocal_geometry_orig.emplace_back(aurostd::string2utype<double>(ktokens[t]));
            }
          }
          if (data.reciprocal_volume_cell_orig == AUROSTD_NAN) { // reciprocal_volume_cell_orig
            data.reciprocal_volume_cell_orig = vpflow_edata_orig.getattachedutype<double>("EDATA::RECIPROCAL_SPACE_VOLUME");
          }
          if (data.reciprocal_lattice_type_orig.empty()) { // reciprocal_lattice_type_orig
            data.reciprocal_lattice_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE");
          }
          if (data.reciprocal_lattice_variation_type_orig.empty()) { // reciprocal_lattice_variation_type_orig
            data.reciprocal_lattice_variation_type_orig = vpflow_edata_orig.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE");
          }
          // DX20190131 - use self-consistent space group orig - START
          if (data.spacegroup_orig == AUROSTD_NAN) { // CO20201111
            data.spacegroup_orig = vpflow_edata_orig.getattachedutype<uint>("SGDATA::SPACE_GROUP_NUMBER"); // CO20201111
            if (AFLOWLIB_VERBOSE) {
              cout << MESSAGE << " SPACEGROUP_ORIG = " << data.spacegroup_orig << endl;
            }
          }
          // DX20190131 - use self-consistent space group orig - END
          if (data.Wyckoff_letters_orig.empty()) { // Wyckoff_letters_orig
            data.Wyckoff_letters_orig = vpflow_edata_orig.getattachedscheme("SGDATA::WYCKOFF_LETTERS");
          }
          if (data.Wyckoff_multiplicities_orig.empty()) { // Wyckoff_multiplicities_orig
            data.Wyckoff_multiplicities_orig = vpflow_edata_orig.getattachedscheme("SGDATA::WYCKOFF_MULTIPLICITIES");
          }
          if (data.Wyckoff_site_symmetries_orig.empty()) { // Wyckoff_site_symmetries_orig
            data.Wyckoff_site_symmetries_orig = vpflow_edata_orig.getattachedscheme("SGDATA::WYCKOFF_SITE_SYMMETRIES");
          }
        }
        // DX20180823 - extract info from xoption - END
        else { // DX20180823 - if no xoption, read from file (safety)
          vector<string> vline_edata;
          aurostd::string2vectorstring(sss.str(), vline_edata);
          for (size_t iline = 0; iline < vline_edata.size(); iline++) {
            if (data.Bravais_lattice_orig.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.Bravais_lattice_orig = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.lattice_variation_orig.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Lattice Variation")) { // Bravais_Lattice_Variation
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.lattice_variation_orig = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.lattice_system_orig.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Lattice System")) { // Lattice_System
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.lattice_system_orig = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.Pearson_symbol_orig.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata[iline], "Superlattice")) { // Pearson
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.Pearson_symbol_orig = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
          }
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_orig << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_orig << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_orig << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_orig << endl;
        }
        // DX20190208 - add AFLOW label/parameters/parameter values - START
        const double anrl_symmetry_tolerance = str_sym.sym_eps;
        xstructure str_anrl = str;
        const uint setting = SG_SETTING_ANRL;
        anrl::structure2anrl(str_anrl, anrl_symmetry_tolerance, setting);
        if (data.aflow_prototype_label_orig.empty()) { // aflow label
          data.aflow_prototype_label_orig = str_anrl.prototype;
        }
        if (data.aflow_prototype_params_list_orig.empty()) { // aflow parameter list
          data.aflow_prototype_params_list_orig = aurostd::joinWDelimiter(str_anrl.prototype_parameter_list, ",");
        }
        if (data.aflow_prototype_params_values_orig.empty()) { // anrl parameter values
          data.aflow_prototype_params_values_orig = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(str_anrl.prototype_parameter_values, _AFLOWLIB_DATA_DOUBLE_PREC_), ",");
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] AFLOW Label = " << data.aflow_prototype_label_orig << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] AFLOW parameter list = " << data.aflow_prototype_params_list_orig << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.ORIG] AFLOW parameter values = " << data.aflow_prototype_params_values_orig << endl;
        }
        // DX20190208 - add AFLOW label/parameters/parameter values - END
      }
    }

    // ONLY TO GET REFERENCE ENERGIES AND INTERCEPT PROTOTYPES
    if (data.nspecies == 1) {
      const string s = str_relax.species.at(0);
      string TXT1;
      string TXT2;
      TXT1 = "XXX if(AUID==\"" + data.vspecies_pp_AUID.at(0);
      if (!aurostd::substring2bool(str_relax.species_pp_version.at(0), "SCAN")) {
        TXT1 += "\" && nKIN) {found=true;groundstate_structure=\"";
      } else {
        TXT1 += "\" && SCAN) {found=true;groundstate_structure=\"";
      }
      TXT2 = "\";groundstate_energy=" + aurostd::utype2string<double>(data.enthalpy_atom, 7) + ";volume_atom=" + aurostd::utype2string<double>(data.volume_atom, 7);
      if (aurostd::abs(data.spin_atom) < 0.1) {
        TXT2 += ";spin_atom=0.0;} // ";
      } else {
        TXT2 += ";spin_atom=" + aurostd::utype2string<double>(data.spin_atom, 7) + ";} // ";
      }
      TXT2 += str_relax.species_pp_version.at(0);
      if (!data.aflow_prototype_label_orig.empty()) {
        cout << data.aflow_prototype_label_orig << endl;
      }

      // A1
      if (data.aflow_prototype_label_orig == "A_cF4_225_a" && (s == "Ac" || s == "Ag" || s == "Al" || s == "Au" || s == "Ca" || s == "Ce" || s == "Cu" || s == "Ir" || s == "La" || s == "Ni" || s == "Pb" ||
                                                               s == "Pd" || s == "Pt" || s == "Rh" || s == "Sr" || s == "Yb" || s == "Ar" || s == "Ne" || s == "Xe" || s == "Kr")) { // A1
        cout << TXT1 << "A1" << TXT2 << endl;
      }
      // A2
      if (data.aflow_prototype_label_orig == "A_cI2_229_a" &&
          (s == "Ba" || s == "Cr" || s == "Fe" || s == "K" || s == "Li" || s == "Mo" || s == "Na" || s == "Nb" || s == "Ta" || s == "V" || s == "W" || s == "Cs" || s == "Eu")) { // A2
        cout << TXT1 << "A2" << TXT2 << endl;
      }
      // A3
      if (data.aflow_prototype_label_orig == "A_hP2_194_c" && (s == "Be" || s == "Cd" || s == "Co" || s == "Dy" || s == "Hf" || s == "Hg" || s == "Ho" || s == "Mg" || s == "Os" || s == "Re" || s == "Ru" ||
                                                               s == "Sc" || s == "Tc" || s == "Ti" || s == "Tl" || s == "Y" || s == "Zn" || s == "Zr" || s == "He")) { // A3
        cout << TXT1 << "A3" << TXT2 << endl;
      }
      // A4
      if (data.aflow_prototype_label_orig == "A_cF8_227_a" && (s == "Ge" || s == "Si")) { // A4
        cout << TXT1 << "A4" << TXT2 << endl;
      }
      // A5
      if (data.aflow_prototype_label_orig == "A_tI4_141_a" && (s == "Sn")) { // A5
        cout << TXT1 << "A5" << TXT2 << endl;
      }
      // A6
      if (data.aflow_prototype_label_orig == "A_tI2_139_a" && (s == "In")) { // A6
        cout << TXT1 << "A6" << TXT2 << endl;
      }
      // A7
      if (data.aflow_prototype_label_orig == "A_hR2_166_c" && (s == "As" || s == "Bi" || s == "Sb" || s == "P")) { // A7 P??
        cout << TXT1 << "A7" << TXT2 << endl;
      }
      // A8
      if (data.aflow_prototype_label_orig == "A_hP3_152_a" && (s == "Se" || s == "Te")) { // A8
        cout << TXT1 << "A8" << TXT2 << endl;
      }
      // A9
      if (data.aflow_prototype_label_orig == "A_hP4_194_bc" && (s == "C")) { // A9
        cout << TXT1 << "A9" << TXT2 << endl;
      }
      // A10
      if (data.aflow_prototype_label_orig == "A_hR1_166_a" && (s == "Hg")) { // A10
        cout << TXT1 << "A10" << TXT2 << endl;
      }
      // A11
      if (data.aflow_prototype_label_orig == "A_oC8_64_f" && (s == "Ga" || s == "Br")) { // A11
        cout << TXT1 << "A11" << TXT2 << endl;
      }
      // A12
      if (data.aflow_prototype_label_orig == "A_cI58_217_ac2g" && (s == "Mn")) { // A12
        cout << TXT1 << "A12" << TXT2 << endl;
      }
      // diatom (A_tP2_123_g) just RICO choice
      if (data.aflow_prototype_label_orig == "A_tP2_123_g" && (s == "O" || s == "N" || s == "F" || s == "H" || s == "Cl")) { // diatom //  cat /tmp/xscrubber_ppAUID.LIB1 | grep diatom | grep '/O'
        cout << TXT1 << "diatom" << TXT2 << endl;
      }
      // A14
      if (data.aflow_prototype_label_orig == "A_oC8_64_f" && (s == "I")) { // A14
        cout << TXT1 << "A14" << TXT2 << endl;
      }
      // A16
      if (data.aflow_prototype_label_orig == "A_oF128_70_4h" && (s == "S")) { // A16
        cout << TXT1 << "A16" << TXT2 << endl;
      }
      if (data.aflow_prototype_label_orig == "A_hR12_166_2h" && (s == "B")) { // ICSD_56992
        cout << TXT1 << "ICSD_56992" << TXT2 << endl;
      }
      if (data.aflow_prototype_label_orig == "A_hR3_166_ac" && (s == "Sm")) { // C19
        cout << TXT1 << "C19" << TXT2 << endl;
      }
      if (data.aflow_prototype_label_orig.empty() && (s == "Xe")) { // ISOLATED
        cout << TXT1 << "isolated" << TXT2 << endl;
      }
    }

    if (flag_EDATA_RELAX_) { // RELAX
      if (!aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_EDATA_RELAX_OUT) && aurostd::FileExist(directory_RAW + "/CONTCAR.relax")) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " EDATA doing relax (CONTCAR.relax) text format: " << directory_RAW << endl;
        }
        str = xstructure(directory_RAW + "/CONTCAR.relax", IOAFLOW_AUTO);
        str_sp.clear();
        str_sc.clear(); // DX20191220 - uppercase to lowercase clear
        stringstream sss;
        sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym = str; // CO20171027 //DX20180226 - set equal str
        aurostd::xoption vpflow_edata_relax; // DX20180823 - added xoption
        sss << pflow::PrintData(str, str_sym, str_sp, str_sc, vpflow_edata_relax, "EDATA", txt_ft, false); // EDATA //CO20171025 //CO20171027 //DX20210301 - void to string output, "txt" to txt_ft
        aurostd::stringstream2file(sss, directory_RAW + "/" + DEFAULT_FILE_EDATA_RELAX_OUT);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " EDATA doing relax (CONTCAR.relax) json format: " << directory_RAW << endl;
        }
        stringstream jjj; // CO20171025
        jjj << pflow::PrintData(str_sym, str_sym, str_sp, str_sc, vpflow_edata_relax, "EDATA", json_ft, true); // EDATA already_calculated! //CO20171025 //CO20171027 //DX20210301 - void to string output, "json" to json_ft
        aurostd::stringstream2file(jjj, directory_RAW + "/" + DEFAULT_FILE_EDATA_RELAX_JSON); // CO20171025
        vcif = {str, str_sp, str_sc};
        // CO+DX START 20170713 - adding symmetry output to RAW
        string new_sym_file;
        for (const auto& [group, out, json] : groups_outs_json) {
          const string out_path = string(directory_RAW).append("/").append(out);
          const string json_path = string(directory_RAW).append("/").append(json);
          // txt variants
          KBIN_SymmetryWrite(FileMESSAGE, str_sp, aflags, group, false, messagestream, "txt");
          if (aurostd::FileExist(out_path)) {
            AddFileNameBeforeExtension(out_path, "relax", new_sym_file);
            aurostd::file2file(out_path, new_sym_file);
          }
          // json variants
          KBIN_SymmetryWrite(FileMESSAGE, str_sp, aflags, group, false, messagestream, "json");
          if (aurostd::FileExist(json_path)) {
            AddFileNameBeforeExtension(json_path, "relax", new_sym_file);
            aurostd::file2file(json_path, new_sym_file);
          }
        }
        // cout << messagestream;
        // CO+DX STOP 20170713 - adding symmetry output to RAW
        //  now extract info
        // DX20180823 - extract info from xoption - START
        if (vpflow_edata_relax.flag("EDATA::CALCULATED")) {
          if (data.Bravais_lattice_relax.empty()) { // Bravais_Lattice
            data.Bravais_lattice_relax = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_TYPE");
          }
          if (data.lattice_variation_relax.empty()) { // Bravais_Lattice_Variation
            data.lattice_variation_relax = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE");
          }
          if (data.lattice_system_relax.empty()) { // Lattice_System
            data.lattice_system_relax = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_SYSTEM");
          }
          if (data.Pearson_symbol_relax.empty()) { // Pearson
            data.Pearson_symbol_relax = vpflow_edata_relax.getattachedscheme("EDATA::PEARSON_SYMBOL");
          }
          if (data.crystal_system.empty()) { // crystal_system
            data.crystal_system = vpflow_edata_relax.getattachedscheme("EDATA::CRYSTAL_SYSTEM");
          }
          if (data.crystal_family.empty()) { // crystal_family
            data.crystal_family = vpflow_edata_relax.getattachedscheme("EDATA::CRYSTAL_FAMILY");
          }
          // DX20190115 - typo in Hermann-Mauguin and crystal class entries - START
          if (data.point_group_Hermann_Mauguin.empty()) { // point_group_Hermann_Mauguin
            data.point_group_Hermann_Mauguin = vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_HERMANN_MAUGUIN");
          }
          if (data.crystal_class.empty()) { // crystal_class
            data.crystal_class = vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_CRYSTAL_CLASS");
          }
          // DX20190115 - typo in Hermann-Mauguin and crystal class entries - END
          if (data.point_group_Schoenflies.empty()) { // point_group_Schoenflies
            data.point_group_Schoenflies = vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_SCHOENFLIES");
          }
          if (data.point_group_orbifold.empty()) { // point_group_orbifold
            data.point_group_orbifold = vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_ORBIFOLD");
          }
          if (data.point_group_type.empty()) { // point_group_type
            data.point_group_type = vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_TYPE");
          }
          if (data.point_group_order == AUROSTD_NAN) { // point_group_order
            data.point_group_order = vpflow_edata_relax.getattachedutype<uint>("EDATA::POINT_GROUP_ORDER");
          }
          if (data.point_group_structure.empty()) { // point_group_structure
            data.point_group_structure = vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_STRUCTURE");
          }
          if (data.Bravais_lattice_lattice_type.empty()) { // Bravais_lattice_lattice_type
            data.Bravais_lattice_lattice_type = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE");
          }
          // DX20190115 - typo, missing "variation" - START
          if (data.Bravais_lattice_lattice_variation_type.empty()) { // Bravais_lattice_lattice_variation_type
            data.Bravais_lattice_lattice_variation_type = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE");
          }
          // DX20190115 - typo, missing "variation" - END
          if (data.Bravais_lattice_lattice_system.empty()) { // Bravais_lattice_lattice_system
            data.Bravais_lattice_lattice_system = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM");
          }
          if (data.Bravais_superlattice_lattice_type.empty()) { // Bravais_superlattice_lattice_type
            data.Bravais_superlattice_lattice_type = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE");
          }
          if (data.Bravais_superlattice_lattice_variation_type.empty()) { // Bravais_superlattice_lattice_variation_type
            data.Bravais_superlattice_lattice_variation_type = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE");
          }
          if (data.Bravais_superlattice_lattice_system.empty()) { // Bravais_superlattice_lattice_system
            data.Bravais_superlattice_lattice_system = vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM");
          }
          if (data.Pearson_symbol_superlattice.empty()) { // Pearson_symbol_superlattice
            data.Pearson_symbol_superlattice = vpflow_edata_relax.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE");
          }
          if (data.reciprocal_geometry_relax.empty()) { // reciprocal_geometry_relax //CO20220719 _relax
            data.reciprocal_geometry_relax = vpflow_edata_relax.getattachedscheme("EDATA::RECIPROCAL_LATTICE_PARAMETERS"); // CO20220719 _relax
            vector<string> ktokens;
            aurostd::string2tokens(data.reciprocal_geometry_relax, ktokens, ","); // CO20220719 _relax
            for (size_t t = 0; t < ktokens.size(); t++) {
              data.vreciprocal_geometry_relax.emplace_back(aurostd::string2utype<double>(ktokens[t]));
            } // CO20220719 _relax
          }
          if (data.reciprocal_volume_cell == AUROSTD_NAN) { // reciprocal_volume_cell
            data.reciprocal_volume_cell = vpflow_edata_relax.getattachedutype<double>("EDATA::RECIPROCAL_SPACE_VOLUME");
          }
          if (data.reciprocal_lattice_type.empty()) { // reciprocal_lattice_type
            data.reciprocal_lattice_type = vpflow_edata_relax.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE");
          }
          if (data.reciprocal_lattice_variation_type.empty()) { // reciprocal_lattice_variation_type
            data.reciprocal_lattice_variation_type = vpflow_edata_relax.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE");
          }
          // DX20190131 - use self-consistent space group relax - START
          if (data.spacegroup_relax == AUROSTD_NAN) { // CO20201111
            data.spacegroup_relax = vpflow_edata_relax.getattachedutype<uint>("SGDATA::SPACE_GROUP_NUMBER"); // CO20201111
            if (AFLOWLIB_VERBOSE) {
              cout << MESSAGE << " SPACEGROUP_RELAX = " << data.spacegroup_relax << endl;
            }
          }
          // DX20190131 - use self-consistent space group orig - END
          if (data.Wyckoff_letters.empty()) { // Wyckoff_letters
            data.Wyckoff_letters = vpflow_edata_relax.getattachedscheme("SGDATA::WYCKOFF_LETTERS");
          }
          if (data.Wyckoff_multiplicities.empty()) { // Wyckoff_multiplicities
            data.Wyckoff_multiplicities = vpflow_edata_relax.getattachedscheme("SGDATA::WYCKOFF_MULTIPLICITIES");
          }
          if (data.Wyckoff_site_symmetries.empty()) { // Wyckoff_site_symmetries
            data.Wyckoff_site_symmetries = vpflow_edata_relax.getattachedscheme("SGDATA::WYCKOFF_SITE_SYMMETRIES");
          }
        }
        // DX20180823 - extract info from xoption - END
        else { // DX20180823 - if no xoption, read from file (safety)
          vector<string> vline_edata;
          aurostd::string2vectorstring(sss.str(), vline_edata);
          for (size_t iline = 0; iline < vline_edata.size(); iline++) {
            if (data.Bravais_lattice_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.Bravais_lattice_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.lattice_variation_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Lattice Variation")) { // Bravais_Lattice_Variation
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.lattice_variation_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.lattice_system_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Lattice System")) { // Lattice_System
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.lattice_system_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.Pearson_symbol_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata[iline], "Superlattice")) { // Pearson
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.Pearson_symbol_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
          }
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_relax << endl;
        }
        // DX20190208 - add AFLOW label/parameters/parameter values - START
        const double anrl_symmetry_tolerance = str_sym.sym_eps;
        xstructure str_anrl = str;
        const uint setting = SG_SETTING_ANRL;
        anrl::structure2anrl(str_anrl, anrl_symmetry_tolerance, setting);
        if (data.aflow_prototype_label_relax.empty()) { // aflow label
          data.aflow_prototype_label_relax = str_anrl.prototype;
        }
        if (data.aflow_prototype_params_list_relax.empty()) { // aflow parameter list
          data.aflow_prototype_params_list_relax = aurostd::joinWDelimiter(str_anrl.prototype_parameter_list, ",");
        }
        if (data.aflow_prototype_params_values_relax.empty()) { // aflow parameter values
          data.aflow_prototype_params_values_relax = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(str_anrl.prototype_parameter_values, _AFLOWLIB_DATA_DOUBLE_PREC_), ",");
        }
        // DX20190208 - add AFLOW label/parameters/parameter values - END
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] AFLOW Label = " << data.aflow_prototype_label_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] AFLOW parameter list = " << data.aflow_prototype_params_list_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.RELAX] AFLOW parameter values = " << data.aflow_prototype_params_values_relax << endl;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [17]" << endl;
    }
    if (flag_EDATA_BANDS_) { // BANDS
      if (!aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_EDATA_BANDS_OUT) && aurostd::FileExist(directory_RAW + "/POSCAR.bands")) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " EDATA doing bands (POSCAR.bands) text format: " << directory_RAW << endl;
        }
        str = xstructure(directory_RAW + "/POSCAR.bands", IOAFLOW_AUTO);
        str_sp.clear();
        str_sc.clear(); // DX20191220 - uppercase to lowercase clear
        stringstream sss;
        sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym = str; // CO20171027 //DX20180226 - set equal str
        aurostd::xoption vpflow_edata_bands; // DX20180823 - added xoption
        sss << pflow::PrintData(str, str_sym, str_sp, str_sc, vpflow_edata_bands, "EDATA", txt_ft, false); // EDATA //CO20171025 //CO20171027 //DX20210301 - void to string output, "txt" to txt_ft
        aurostd::stringstream2file(sss, directory_RAW + "/" + DEFAULT_FILE_EDATA_BANDS_OUT);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " EDATA doing bands (POSCAR.bands) json format: " << directory_RAW << endl;
        }
        stringstream jjj; // CO20171025
        jjj << pflow::PrintData(str_sym, str_sym, str_sp, str_sc, vpflow_edata_bands, "EDATA", json_ft, true); // EDATA already_calculated! //CO20171025 //CO20171027 //DX20210301 - void to string output, "json" to json_ft
        aurostd::stringstream2file(jjj, directory_RAW + "/" + DEFAULT_FILE_EDATA_BANDS_JSON); // CO20171025
        vcif = {str, str_sp, str_sc};
        // CO+DX START 20170713 - adding symmetry output to RAW
        string new_sym_file;
        for (const auto& [group, out, json] : groups_outs_json) {
          const string out_path = string(directory_RAW).append("/").append(out);
          const string json_path = string(directory_RAW).append("/").append(json);
          // txt variants
          KBIN_SymmetryWrite(FileMESSAGE, str_sp, aflags, group, false, messagestream, "txt");
          if (aurostd::FileExist(out_path)) {
            AddFileNameBeforeExtension(out_path, "bands", new_sym_file);
            aurostd::file2file(out_path, new_sym_file);
          }
          // json variants
          KBIN_SymmetryWrite(FileMESSAGE, str_sp, aflags, group, false, messagestream, "json");
          if (aurostd::FileExist(json_path)) {
            AddFileNameBeforeExtension(json_path, "bands", new_sym_file);
            aurostd::file2file(json_path, new_sym_file);
          }
        }
        // cout << messagestream;
        // CO+DX STOP 20170713 - adding symmetry output to RAW
        //  now extract info
        // DX20180823 - extract info from xoption - START
        if (vpflow_edata_bands.flag("EDATA::CALCULATED")) {
          if (data.Bravais_lattice_relax.empty()) { // Bravais_Lattice
            data.Bravais_lattice_relax = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_TYPE");
          }
          if (data.lattice_variation_relax.empty()) { // Bravais_Lattice_Variation
            data.lattice_variation_relax = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE");
          }
          if (data.lattice_system_relax.empty()) { // Lattice_System
            data.lattice_system_relax = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_SYSTEM");
          }
          if (data.Pearson_symbol_relax.empty()) { // Pearson
            data.Pearson_symbol_relax = vpflow_edata_bands.getattachedscheme("EDATA::PEARSON_SYMBOL");
          }
          if (data.crystal_system.empty()) { // crystal_system
            data.crystal_system = vpflow_edata_bands.getattachedscheme("EDATA::CRYSTAL_SYSTEM");
          }
          if (data.crystal_family.empty()) { // crystal_family
            data.crystal_family = vpflow_edata_bands.getattachedscheme("EDATA::CRYSTAL_FAMILY");
          }
          // DX20190115 - typo in Hermann-Mauguin and crystal class entries - START
          if (data.point_group_Hermann_Mauguin.empty()) { // point_group_Hermann_Mauguin
            data.point_group_Hermann_Mauguin = vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_HERMANN_MAUGUIN");
          }
          if (data.crystal_class.empty()) { // crystal_class
            data.crystal_class = vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_CRYSTAL_CLASS");
          }
          // DX20190115 - typo in Hermann-Mauguin and crystal class entries - END
          if (data.point_group_Schoenflies.empty()) { // point_group_Schoenflies
            data.point_group_Schoenflies = vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_SCHOENFLIES");
          }
          if (data.point_group_orbifold.empty()) { // point_group_orbifold
            data.point_group_orbifold = vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_ORBIFOLD");
          }
          if (data.point_group_type.empty()) { // point_group_type
            data.point_group_type = vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_TYPE");
          }
          if (data.point_group_order == AUROSTD_NAN) { // point_group_order
            data.point_group_order = vpflow_edata_bands.getattachedutype<uint>("EDATA::POINT_GROUP_ORDER");
          }
          if (data.point_group_structure.empty()) { // point_group_structure
            data.point_group_structure = vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_STRUCTURE");
          }
          if (data.Bravais_lattice_lattice_type.empty()) { // Bravais_lattice_lattice_type
            data.Bravais_lattice_lattice_type = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE");
          }
          // DX20190115 - typo, missing "variation" - START
          if (data.Bravais_lattice_lattice_variation_type.empty()) { // Bravais_lattice_lattice_variation_type
            data.Bravais_lattice_lattice_variation_type = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE");
          }
          // DX20190115 - typo, missing "variation" - END
          if (data.Bravais_lattice_lattice_system.empty()) { // Bravais_lattice_lattice_system
            data.Bravais_lattice_lattice_system = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM");
          }
          if (data.Bravais_superlattice_lattice_type.empty()) { // Bravais_superlattice_lattice_type
            data.Bravais_superlattice_lattice_type = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE");
          }
          if (data.Bravais_superlattice_lattice_variation_type.empty()) { // Bravais_superlattice_lattice_variation_type
            data.Bravais_superlattice_lattice_variation_type = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE");
          }
          if (data.Bravais_superlattice_lattice_system.empty()) { // Bravais_superlattice_lattice_system
            data.Bravais_superlattice_lattice_system = vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM");
          }
          if (data.Pearson_symbol_superlattice.empty()) { // Pearson_symbol_superlattice
            data.Pearson_symbol_superlattice = vpflow_edata_bands.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE");
          }
          if (data.reciprocal_geometry_relax.empty()) { // reciprocal_geometry_relax //CO20220719 _relax
            data.reciprocal_geometry_relax = vpflow_edata_bands.getattachedscheme("EDATA::RECIPROCAL_LATTICE_PARAMETERS"); // CO20220719 _relax
            vector<string> ktokens;
            aurostd::string2tokens(data.reciprocal_geometry_relax, ktokens, ","); // CO20220719 _relax
            for (size_t t = 0; t < ktokens.size(); t++) {
              data.vreciprocal_geometry_relax.emplace_back(aurostd::string2utype<double>(ktokens[t]));
            } // CO20220719 _relax
          }
          if (data.reciprocal_volume_cell == AUROSTD_NAN) { // reciprocal_volume_cell
            data.reciprocal_volume_cell = vpflow_edata_bands.getattachedutype<double>("EDATA::RECIPROCAL_SPACE_VOLUME");
          }
          if (data.reciprocal_lattice_type.empty()) { // reciprocal_lattice_type
            data.reciprocal_lattice_type = vpflow_edata_bands.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE");
          }
          if (data.reciprocal_lattice_variation_type.empty()) { // reciprocal_lattice_variation_type
            data.reciprocal_lattice_variation_type = vpflow_edata_bands.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE");
          }
          // DX20190131 - use self-consistent space group relax - START
          if (data.spacegroup_relax == AUROSTD_NAN) { // CO20201111
            data.spacegroup_relax = vpflow_edata_bands.getattachedutype<uint>("SGDATA::SPACE_GROUP_NUMBER"); // CO20201111
            if (AFLOWLIB_VERBOSE) {
              cout << MESSAGE << " SPACEGROUP_RELAX = " << data.spacegroup_relax << endl;
            }
          }
          // DX20190131 - use self-consistent space group orig - END
          if (data.Wyckoff_letters.empty()) { // Wyckoff_letters
            data.Wyckoff_letters = vpflow_edata_bands.getattachedscheme("SGDATA::WYCKOFF_LETTERS");
          }
          if (data.Wyckoff_multiplicities.empty()) { // Wyckoff_multiplicities
            data.Wyckoff_multiplicities = vpflow_edata_bands.getattachedscheme("SGDATA::WYCKOFF_MULTIPLICITIES");
          }
          if (data.Wyckoff_site_symmetries.empty()) { // Wyckoff_site_symmetries
            data.Wyckoff_site_symmetries = vpflow_edata_bands.getattachedscheme("SGDATA::WYCKOFF_SITE_SYMMETRIES");
          }
        }
        // DX20180823 - extract info from xoption - END
        else { // DX20180823 - if no xoption, read from file (safety)
          vector<string> vline_edata;
          aurostd::string2vectorstring(sss.str(), vline_edata);
          for (size_t iline = 0; iline < vline_edata.size(); iline++) {
            if (data.Bravais_lattice_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.Bravais_lattice_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.lattice_variation_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Lattice Variation")) { // Bravais_Lattice_Variation
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.lattice_variation_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.lattice_system_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Lattice System")) { // Lattice_System
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.lattice_system_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
            if (data.Pearson_symbol_relax.empty() && aurostd::substring2bool(vline_edata[iline], "Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata[iline], "Superlattice")) { // Pearson
              aurostd::string2tokens(vline_edata[iline], tokens, "=");
              data.Pearson_symbol_relax = aurostd::RemoveWhiteSpaces(tokens.at(tokens.size() - 1));
            }
          }
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_relax << endl;
        }
        // DX20190208 - add AFLOW label/parameters/parameter values - START
        const double anrl_symmetry_tolerance = str_sym.sym_eps;
        xstructure str_anrl = str;
        const uint setting = SG_SETTING_ANRL;
        anrl::structure2anrl(str_anrl, anrl_symmetry_tolerance, setting);
        if (data.aflow_prototype_label_relax.empty()) { // anrl label
          data.aflow_prototype_label_relax = str_anrl.prototype;
        }
        if (data.aflow_prototype_params_list_relax.empty()) { // anrl parameter list
          data.aflow_prototype_params_list_relax = aurostd::joinWDelimiter(str_anrl.prototype_parameter_list, ",");
        }
        if (data.aflow_prototype_params_values_relax.empty()) { // anrl parameter values
          data.aflow_prototype_params_values_relax = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(str_anrl.prototype_parameter_values, _AFLOWLIB_DATA_DOUBLE_PREC_), ",");
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] AFLOW Label = " << data.aflow_prototype_label_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] AFLOW parameter list = " << data.aflow_prototype_params_list_relax << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " [EDATA.BANDS] AFLOW parameter values = " << data.aflow_prototype_params_values_relax << endl;
        }
        // DX20190208 - add AFLOW label/parameters/parameter values - END
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [18]" << endl;
    }
    // PERFORM DATA STEP
    if (flag_DATA_ORIG_ || flag_DATA_RELAX_) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " DATA start: " << directory_RAW << endl;
      }
    }
    if (flag_DATA_ORIG_) { // ORIG
      if (!aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_DATA_ORIG_OUT) && aurostd::FileExist(directory_RAW + "/POSCAR.orig")) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " DATA doing orig (POSCAR.orig) text format: " << directory_RAW << endl;
        }
        str = xstructure(directory_RAW + "/POSCAR.orig", IOAFLOW_AUTO);
        str_sp.clear();
        str_sc.clear(); // DX20191220 - uppercase to lowercase clear
        stringstream sss;
        sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym = str; // CO20171027 //DX20180226 - set equal str
        sss << pflow::PrintData(str, str_sym, str_sp, str_sc, "DATA", txt_ft, false); // DATA //CO20171025 //CO20171027 //DX20210301 - void to string output, "txt" to txt_ft
        aurostd::stringstream2file(sss, directory_RAW + "/" + DEFAULT_FILE_DATA_ORIG_OUT);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " DATA doing orig (POSCAR.orig) json format: " << directory_RAW << endl;
        }
        stringstream jjj; // CO20171025
        jjj << pflow::PrintData(str_sym, str_sym, str_sp, str_sc, "DATA", json_ft, true); // DATA already_calculated! //CO20171025 //CO20171027 //DX20210301 - void to string output, "json" to json_ft
        aurostd::stringstream2file(jjj, directory_RAW + "/" + DEFAULT_FILE_DATA_ORIG_JSON); // CO20171025
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [19]" << endl;
    }
    if (flag_DATA_RELAX_) { // RELAX
      if (!aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_DATA_RELAX_OUT) && aurostd::FileExist(directory_RAW + "/CONTCAR.relax")) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " DATA doing relax (CONTCAR.relax) text format: " << directory_RAW << endl;
        }
        str = xstructure(directory_RAW + "/CONTCAR.relax", IOAFLOW_AUTO);
        str_sp.clear();
        str_sc.clear(); // DX20191220 - uppercase to lowercase clear
        stringstream sss;
        sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym = str; // CO20171027 //DX20180226 - set equal str
        sss << pflow::PrintData(str, str_sym, str_sp, str_sc, "DATA", txt_ft, false); // DATA //CO20171025 //CO20171027 //DX20210301 - void to string output, "txt" to txt_ft
        aurostd::stringstream2file(sss, directory_RAW + "/" + DEFAULT_FILE_DATA_RELAX_OUT);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " DATA doing relax (CONTCAR.relax) json format: " << directory_RAW << endl;
        }
        stringstream jjj; // CO20171025
        jjj << pflow::PrintData(str_sym, str_sym, str_sp, str_sc, "DATA", json_ft, true); // DATA already_calculated! //CO20171025 //CO20171027 //DX20210301 - void to string output, "json" to json_ft
        aurostd::stringstream2file(jjj, directory_RAW + "/" + DEFAULT_FILE_DATA_RELAX_JSON); // CO20171025
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [20]" << endl;
    }
    if (flag_DATA_BANDS_) { // BANDS
      if (!aurostd::FileExist(directory_RAW + "/" + DEFAULT_FILE_DATA_BANDS_OUT) && aurostd::FileExist(directory_RAW + "/POSCAR.bands")) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " DATA doing relax (POSCAR.bands) text format: " << directory_RAW << endl;
        }
        str = xstructure(directory_RAW + "/POSCAR.bands", IOAFLOW_AUTO);
        str_sp.clear();
        str_sc.clear(); // DX20191220 - uppercase to lowercase clear
        stringstream sss;
        sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym = str; // CO20171027 //DX20180226 - set equal str
        sss << pflow::PrintData(str, str_sym, str_sp, str_sc, "DATA", txt_ft, false); // DATA //CO20171025 //CO20171027 //DX20210301 - void to string output, "txt" to txt_ft
        aurostd::stringstream2file(sss, directory_RAW + "/" + DEFAULT_FILE_DATA_BANDS_OUT);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " DATA doing relax (POSCAR.bands) json format: " << directory_RAW << endl;
        }
        stringstream jjj; // CO20171025
        jjj << pflow::PrintData(str_sym, str_sym, str_sp, str_sc, "DATA", json_ft, true); // DATA already_calculated! //CO20171025 //CO20171027 //DX20210301 - void to string output, "json" to json_ft
        aurostd::stringstream2file(jjj, directory_RAW + "/" + DEFAULT_FILE_DATA_BANDS_JSON); // CO20171025
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [21]" << endl;
    }
    // CREATING CIFS
    if (vcif.empty() && aurostd::FileExist(directory_RAW + "/POSCAR.bands")) {
      vcif.emplace_back(directory_RAW + "/POSCAR.bands", IOVASP_AUTO);
    }
    const xvector<double> nvec(3);
    nvec(1) = 1;
    nvec(2) = 1;
    nvec(3) = 1;
    const double angle = 45;

    for (size_t j = 0; j < vcif.size(); j++) {
      vcif[j] = GetLTFVCell(nvec, angle, vcif[j]);
      if (j == 0) {
        cout << MESSAGE << " CIF creation: data.spacegroup_relax=" << data.spacegroup_relax << endl;
      }
      if (j == 0) {
        cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal" << endl;
      }
      if (j == 1) {
        cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim" << endl;
      }
      if (j == 2) {
        cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv" << endl;
      }
      stringstream oss;
      //    for(size_t i=0;i<data.vspecies.size();i++) cerr << XPID << data.vspecies.at(i) << endl;
      vcif[j].species.assign(data.vspecies.begin(), data.vspecies.end());
      vcif[j].species_pp.assign(data.vspecies.begin(), data.vspecies.end());
      vcif[j].species_pp_type.assign(data.vspecies.size(), "");
      vcif[j].species_pp_version.assign(data.vspecies.size(), "");
      vcif[j].species_pp_ZVAL.assign(data.vspecies.size(), 0.0);
      vcif[j].species_pp_vLDAU.assign(data.vspecies.size(), deque<double>());
      for (size_t i = 0; i < vcif[j].atoms.size(); i++) {
        vcif[j].atoms[i].name = vcif[j].species.at(vcif[j].atoms[i].type);
      }
      for (size_t i = 0; i < vcif[j].atoms.size(); i++) {
        vcif[j].atoms[i].cleanname = vcif[j].species.at(vcif[j].atoms[i].type);
      }
      pflow::PrintCIF(oss, vcif[j], 1); // aurostd::string2utype<int>(data.spacegroup_relax));
      if (j == 0) {
        aurostd::stringstream2file(oss, directory_RAW + "/" + KBIN::ExtractSystemName(directory_LIB) + ".cif");
      }
      if (j == 1) {
        aurostd::stringstream2file(oss, directory_RAW + "/" + KBIN::ExtractSystemName(directory_LIB) + "_sprim.cif");
      }
      if (j == 2) {
        aurostd::stringstream2file(oss, directory_RAW + "/" + KBIN::ExtractSystemName(directory_LIB) + "_sconv.cif");
      }
      vcif[j].AddCorners();
      oss.clear();
      oss.str("");
      pflow::PrintCIF(oss, vcif[j], 1); // aurostd::string2utype<int>(data.spacegroup_relax));
      if (j == 0) {
        cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal_corner" << endl;
      }
      if (j == 1) {
        cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim_corner" << endl;
      }
      if (j == 2) {
        cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv_corner" << endl;
      }
      if (j == 0) {
        aurostd::stringstream2file(oss, directory_RAW + "/" + KBIN::ExtractSystemName(directory_LIB) + "_corner.cif");
      }
      if (j == 1) {
        aurostd::stringstream2file(oss, directory_RAW + "/" + KBIN::ExtractSystemName(directory_LIB) + "_sprim_corner.cif");
      }
      if (j == 2) {
        aurostd::stringstream2file(oss, directory_RAW + "/" + KBIN::ExtractSystemName(directory_LIB) + "_sconv_corner.cif");
      }
    }

    // CREATING GEOMETRY
    if (aurostd::FileExist(directory_RAW + "/CONTCAR.relax")) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " LOAD SPECIES=" << data.species << endl;
      }
      xstructure xstr(directory_RAW + "/CONTCAR.relax", IOAFLOW_AUTO);
      xstr.SetSpecies(deq_species);
      const vector<std::pair<string, int>> file_io{
          {  "vasp",   IOVASP_AUTO},
          {    "qe",     IOQE_GEOM},
          {"abinit", IOABINIT_GEOM},
          {  "aims",   IOAIMS_GEOM}
      };
      for (const auto& [file_type, iomode] : file_io) {
        stringstream temp_output;
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " WRITE relax positions for " << file_type << endl;
        }
        if (file_type != "vasp") {
          xstr.ReScale(1.0);
          xstr.neg_scale = false;
          xstr.coord_flag = _COORDS_FRACTIONAL_;
          if (xstr.title.empty()) {
            xstr.buildGenericTitle();
          }
        }
        xstr.iomode = iomode;
        temp_output << xstr;
        aurostd::stringstream2file(temp_output, directory_RAW + "/CONTCAR.relax." + file_type);
      }
    }

    if (flag_LIB2 == true) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " PATCHING for LIB2: " << directory_RAW << endl;
      }
      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        filecar_map["aflow"].lib = directory_LIB + "/" + DEFAULT_AFLOW_QMVASP_OUT + XHOST.vext[iext];
        filecar_map["aflow"].raw = directory_RAW + "/" + DEFAULT_AFLOW_QMVASP_OUT + XHOST.vext[iext];
        if (aurostd::FileExist(filecar_map["aflow"].lib)) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB, DEFAULT_AFLOW_QMVASP_OUT + XHOST.vext[iext], directory_RAW, DEFAULT_AFLOW_QMVASP_OUT + XHOST.vext[iext], vfiles, MESSAGE);
        }
      }

      for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
        filecar_map["POSCAR"].lib = directory_LIB + "/POSCAR.relax2" + XHOST.vext[iext];
        filecar_map["POSCAR"].raw = directory_RAW + "/POSCAR.relax2" + XHOST.vext[iext];
        if (aurostd::FileExist(filecar_map["POSCAR"].lib)) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB, "POSCAR.relax2" + XHOST.vext[iext], directory_RAW, "POSCAR.relax2" + XHOST.vext[iext], vfiles, MESSAGE);
        }
      }
      aurostd::CopyFile(directory_RAW + "/CONTCAR.relax", directory_RAW + "/CONTCAR.relax2");
      vfiles.emplace_back("CONTCAR.relax2");
    }

    // LINKING
    if (!FileName_OUTCAR_relax.empty()) {
      const string& FROM = FileName_OUTCAR_relax;
      string TO = FileName_OUTCAR_relax;
      aurostd::StringSubstInPlace(TO, "LIB/", "RAW/");
      for (const auto& entry : {"relax1", "relax2", "relax3", "static", "bands"}) {
        aurostd::StringSubstInPlace(TO, entry, "relax");
      }
      cout << MESSAGE << " linking " << FROM << "->" << TO << endl;
      aurostd::LinkFile(FROM, TO);
    }

    // DONE
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Magnetic
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for magnetic results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_Magnetic(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    // CO20180130 - note that here we extract STATIC/BANDS properties only
    // some spin properties (spin/cell, spin/atom, spinD) can be extracted from relax2
    // we do this in the thermo loop and DO NOT attempt to redo here
    // spinF is nonzero for magnetic metals, so we need to determine Egap here as well (need BANDS)
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("magnetic");
    xOUTCAR outcar_bands;
    xOUTCAR outcar_static;
    xDOSCAR doscar;

    const double MAG_EPS = 1e-6;

    // try to grab STATIC properties first, Efermi and spin/cell, spin/atom, spinD
    // NB, spin properties from relax2 are okay (lower kppra), we grab it in thermo loop, so we don't attempt to here
    double EFERMI = AUROSTD_NAN;
    outcar_static.clear();
    if (aurostd::FileExist(directory_LIB + "/" + "OUTCAR.static") || aurostd::CompressFileExist(directory_LIB + "/" + "OUTCAR.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "OUTCAR.static", directory_RAW, "OUTCAR.static", vfiles, MESSAGE); // OUTCAR.static
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "OUTCAR.static") << endl;
      }
      if (outcar_static.GetPropertiesFile(directory_RAW + "/" + "OUTCAR.static")) {
        EFERMI = outcar_static.Efermi;
        data.spin_cell = outcar_static.mag_cell;
        data.spin_atom = outcar_static.mag_atom;
        data.spinD = "";
        data.vspinD.clear();
        if (!outcar_static.vmag.empty()) {
          for (size_t i = 0; i < outcar_static.vmag.size(); i++) {
            data.spinD += aurostd::utype2string<double>(outcar_static.vmag[i], 5) + (i < outcar_static.vmag.size() - 1 ? "," : "");
            data.vspinD.emplace_back(outcar_static.vmag[i]);
          }
        } else {
          for (uint i = 0; i < outcar_static.natoms; i++) { // use outcar_static.natoms as there can be a primitivization between relax and static
            data.spinD += aurostd::utype2string<double>(0) + (i < outcar_static.natoms - 1 ? "," : "");
            data.vspinD.push_back(0.0);
          }
        }
      } else {
        cout << MESSAGE << " ERROR OUTCAR.static properties cannot be extracted" << endl;
      } // CO20200404 - removed xOUTCAR.ERROR
    } else {
      cout << MESSAGE << " MISSING OUTCAR.static" << endl;
    }
    if (EFERMI == AUROSTD_NAN) {
      cout << MESSAGE << " unable to load OUTCAR.static, using Efermi from bands" << endl;
    }

    // ideally we grab both bands and static
    // we want static because we trust this Efermi (self-consistent)
    // however, we fall back on Efermi from bands

    // BANDS CALCULATION
    data.spinF = 0.0; // DEFAULT
    outcar_bands.clear();
    doscar.clear();
    if (aurostd::FileExist(directory_LIB + "/" + "OUTCAR.bands") || aurostd::CompressFileExist(directory_LIB + "/" + "OUTCAR.bands")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "OUTCAR.bands", directory_RAW, "OUTCAR.bands", vfiles, MESSAGE); // OUTCAR.bands
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "OUTCAR.bands") << endl;
      }
      if (outcar_bands.GetPropertiesFile(directory_RAW + "/" + "OUTCAR.bands")) {
        if (outcar_bands.GetBandGap(EFERMI)) {
          data.Egap = outcar_bands.Egap_net;
          data.Egap_fit = outcar_bands.Egap_fit_net;
          data.Egap_type = outcar_bands.Egap_type_net;
          if (aurostd::substring2bool(data.Egap_type, "metal")) {
            data.Egap = 0.0; // half-metal
          }
          if (aurostd::substring2bool(data.Egap_type, "metal")) {
            data.Egap_fit = 0.0; // half-metal
          }

          // SPIN POLARIZATION AT FERMI LEVEL
          if (data.Egap < MAG_EPS && aurostd::abs(data.spin_cell) > MAG_EPS) { // must be metal and magnetic
            if (doscar.content.empty() && (aurostd::FileExist(directory_LIB + "/" + "DOSCAR.static") || aurostd::CompressFileExist(directory_LIB + "/" + "DOSCAR.static"))) {
              aflowlib::LIB2RAW_FileNeeded(directory_LIB, "DOSCAR.static", directory_RAW, "DOSCAR.static", vfiles, MESSAGE); // DOSCAR.static
              if (AFLOWLIB_VERBOSE) {
                cout << MESSAGE << " loading " << string(directory_RAW + "/" + "DOSCAR.static") << endl;
              }
              doscar.GetPropertiesFile(directory_RAW + "/" + "DOSCAR.static");
            }
            if (doscar.content.empty() && (aurostd::FileExist(directory_LIB + "/" + "DOSCAR.relax2") || aurostd::CompressFileExist(directory_LIB + "/" + "DOSCAR.relax2"))) {
              aflowlib::LIB2RAW_FileNeeded(directory_LIB, "DOSCAR.relax2", directory_RAW, "DOSCAR.relax", vfiles, MESSAGE); // DOSCAR.relax2
              if (AFLOWLIB_VERBOSE) {
                cout << MESSAGE << " loading " << string(directory_RAW + "/" + "DOSCAR.relax") << endl;
              }
              doscar.GetPropertiesFile(directory_RAW + "/" + "DOSCAR.relax");
            }
            if (!doscar.content.empty()) {
              data.spinF = doscar.spinF;
            } else {
              cout << MESSAGE << " MISSING DOSCAR.static and DOSCAR.relax[2]" << endl;
            }
          }
        } else {
          cout << MESSAGE << " ERROR OUTCAR.bands BandGap() cannot be extracted" << endl;
        } // CO20181129  //CO20200404 - removed xOUTCAR.ERROR
      } else {
        cout << MESSAGE << " ERROR OUTCAR.bands properties cannot be extracted" << endl;
      } // CO20181129 //CO20200404 - removed xOUTCAR.ERROR
    } else {
      cout << MESSAGE << " MISSING OUTCAR.bands" << endl;
    } // CO20181129

    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " spin/cell = " << data.spin_cell << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " spin/atom = " << data.spin_atom << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " spinD     = " << data.spinD << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " spinF     = " << ((data.spinF != AUROSTD_NAN) ? aurostd::utype2string(data.spinF, 5) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " Egap (eV)     = " << ((data.Egap != AUROSTD_NAN) ? aurostd::utype2string(data.Egap, 5) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " Egap_fit (eV) = " << ((data.Egap_fit != AUROSTD_NAN) ? aurostd::utype2string(data.Egap_fit, 5) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " Egap_type     = " << ((!data.Egap_type.empty()) ? data.Egap_type : "unavailable") << endl;
    }

    // done
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bader  //CO WORK HERE
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for bader charge results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_Bader(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("bader");

    if (aurostd::FileExist(directory_LIB + "/" + "CHGCAR.static") || aurostd::CompressFileExist(directory_LIB + "/" + "CHGCAR.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "CHGCAR.static", directory_RAW, "CHGCAR.static", vfiles, MESSAGE); // CHGCAR.static
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "CHGCAR.static") << endl;
      }
    } else {
      cout << MESSAGE << " MISSING CHGCAR.static" << endl;
      return;
    }
    if (aurostd::FileExist(directory_LIB + "/" + "AECCAR0.static") || aurostd::CompressFileExist(directory_LIB + "/" + "AECCAR0.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "AECCAR0.static", directory_RAW, "AECCAR0.static", vfiles, MESSAGE); // AECCAR0.static
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "AECCAR0.static") << endl;
      }
    } else {
      cout << MESSAGE << " MISSING AECCAR0.static" << endl;
      return;
    }
    if (aurostd::FileExist(directory_LIB + "/" + "AECCAR2.static") || aurostd::CompressFileExist(directory_LIB + "/" + "AECCAR2.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "AECCAR2.static", directory_RAW, "AECCAR2.static", vfiles, MESSAGE); // AECCAR2.static
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "AECCAR2.static") << endl;
      }
    } else {
      cout << MESSAGE << " MISSING AECCAR2.static" << endl;
      return;
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [2]" << endl;
    }

    vector<string> vspecies;
    aurostd::string2tokens(data.species, vspecies, ",");
    for (size_t i = 0; i < vspecies.size(); i++) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " " << vspecies[i] << endl;
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [3]" << endl;
    }

    vector<double> vspecies_pp_ZVAL = data.vspecies_pp_ZVAL;
    for (size_t i = 0; i < vspecies_pp_ZVAL.size(); i++) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " " << vspecies_pp_ZVAL[i] << endl;
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [4]" << endl;
    }

    deque<int> num_each_type;
    xstructure xstr_relaxed = KBIN::GetMostRelaxedStructure(directory_LIB);
    for (size_t i = 0; i < xstr_relaxed.num_each_type.size(); i++) {
      num_each_type.emplace_back(xstr_relaxed.num_each_type[i]);
    } // vector vs. deque
    if (LDEBUG) {
      for (size_t i = 0; i < num_each_type.size(); i++) {
        cerr << XPID << __AFLOW_FUNC__ << " num_each_type[" << i << "]=" << num_each_type[i] << endl;
      }
    }

    vector<double> vZVAL;
    for (size_t i = 0; i < vspecies_pp_ZVAL.size(); i++) {
      for (uint j = 0; j < (uint) num_each_type.at(i); j++) {
        vZVAL.emplace_back(vspecies_pp_ZVAL[i]);
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [5]" << endl;
    }

    // flags
    aurostd::xoption bader_flags;
    // DX+CO START
    bader_flags.flag("BADER::AFLOWLIB_LIBRARY", true); // net charge file
    bader_flags.flag("BADER::JVXL_ALL_SPECIES", true); // for flag automation, no need to specify cutoffs,downsample_ratios here
    bader_flags.flag("BADER::JVXL_CYCLIC", true); // CYCLIC MODE
    // bader_flags.flag("BADER::JVXL_CYCLIC",false);   //SETS MODE
    // DX+CO END

    const vector<double> cutoffs = {0.10, 0.20, 0.25, 0.30, 0.333333333333333, 0.40, 0.50, 0.60, 0.666666666666666, 0.75, 0.80, 0.90};
    vector<int> downsample_ratios;
    downsample_ratios.push_back(2); // for all

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [6]" << endl;
    }

    string bader_options;
    ostringstream oss;
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " beginning BADER calculation--please be patient" << endl;
    }
    bader_functions::Flags2BaderCommands(bader_flags, bader_options, oss);
    // DX+CO START
    bader_functions::BaderCalc(bader_flags, bader_options, KBIN::ExtractSystemName(directory_LIB), vspecies, num_each_type, vZVAL, cutoffs, downsample_ratios, directory_RAW,
                               oss); // CO+ME20200601 - no data.prototype, use ExtractSystemName() //data.prototype
    // DX+CO END
    if (AFLOWLIB_VERBOSE) {
      cout << oss.str();
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [7]" << endl;
    }

    vector<string> vline;
    vector<string> tokens;
    stringstream abader_ss;
    string abader_file = data.system_name + "_abader.out"; // CO20200731 - not prototype, title
    abader_ss.clear();
    if (aurostd::CompressFileExist(abader_file, abader_file)) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(abader_file) << endl;
      }
      abader_ss.str(aurostd::substring2string(aurostd::compressfile2string(abader_file), "[BADER_RESULTS]START", "[BADER_RESULTS]STOP", 1));
      aurostd::stream2vectorstring(abader_ss, vline);
      for (size_t i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "bader_net_charges") {
            data.bader_net_charges = tokens[1];
          }
          if (tokens[0] == "bader_atomic_volumes") {
            data.bader_atomic_volumes = tokens[1];
          }
        }
      }
      aurostd::string2tokens<double>(data.bader_net_charges, data.vbader_net_charges, ","); // be careful, sets precision to 20
      aurostd::string2tokens<double>(data.bader_atomic_volumes, data.vbader_atomic_volumes, ","); // be careful, sets precision to 20
      if (data.vbader_net_charges.empty()) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " data.vbader_net_charges was not extracted correctly" << endl;
        }
        return;
      }
      if (data.vbader_net_charges.size() != data.vbader_atomic_volumes.size()) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " length of data.vbader_net_charges does not match length of data.vbader_atomic_volumes" << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " there was a problem with the bader data extraction" << endl;
        }
        return;
      }
      stringstream num_prec;
      if (AFLOWLIB_VERBOSE) {
        for (size_t i = 0; i < data.vbader_net_charges.size(); i++) {
          num_prec << std::fixed << setprecision(4) << data.vbader_net_charges[i];
          cout << MESSAGE << " bader net charge for atom " << i + 1 << " = " << num_prec.str() << endl;
          num_prec.str("");
        }
        for (size_t i = 0; i < data.vbader_atomic_volumes.size(); i++) {
          num_prec << std::fixed << setprecision(4) << data.vbader_atomic_volumes[i];
          cout << MESSAGE << " bader atomic volume for atom " << i + 1 << " = " << num_prec.str() << endl;
          num_prec.str("");
        }
      }
    } else {
      return;
    }
    // done
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_AGL  // CORMAC
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for AGL results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_AGL(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("agl");
    vector<string> vline;
    vector<string> tokens;
    stringstream aflow_agl_out;
    string aflow_agl_raw_file = directory_RAW + "/" + "aflow.agl.out";

    // AFLOW AGL
    aflow_agl_out.clear();
    if (aurostd::CompressFileExist(directory_LIB + "/" + "aflow.agl.out")) {
      const vector<string> agl_file_list{"aflow.agl.out",
                                         "AGL.out",
                                         "AGL_energies_temperature.out",
                                         "AGL_thermal_properties_temperature.out",
                                         "AGL_edos_gap_pressure.out",
                                         "AGL_edos_gap_pressure.json",
                                         "AGL_energy.json",
                                         "AGL_energy_structures.json",
                                         "AGL_energy_volume.out",
                                         "AGL_gibbs_energy_pT.out",
                                         "AGL_Hugoniot.out"};
      for (const auto& file : agl_file_list) {
        aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "aflow.agl.out") << endl;
      }
      if (!(aurostd::CompressFileExist(aflow_agl_raw_file, aflow_agl_raw_file))) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File " + aflow_agl_raw_file + " not found", _FILE_NOT_FOUND_);
      }
      aflow_agl_out.str(aurostd::substring2string(aurostd::compressfile2string(aflow_agl_raw_file), "[AGL_RESULTS]START", "[AGL_RESULTS]STOP", 1));
      aurostd::stream2vectorstring(aflow_agl_out, vline);
      for (size_t i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "agl_thermal_conductivity_300K") {
            data.agl_thermal_conductivity_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_debye") {
            data.agl_debye = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_acoustic_debye") {
            data.agl_acoustic_debye = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_gruneisen") {
            data.agl_gruneisen = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_heat_capacity_Cv_300K") {
            data.agl_heat_capacity_Cv_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_heat_capacity_Cp_300K") {
            data.agl_heat_capacity_Cp_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_thermal_expansion_300K") {
            data.agl_thermal_expansion_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_bulk_modulus_static_300K") {
            data.agl_bulk_modulus_static_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_bulk_modulus_isothermal_300K") {
            data.agl_bulk_modulus_isothermal_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "agl_poisson_ratio_source") {
            data.agl_poisson_ratio_source = tokens[1]; // CT20181212
          }
          if (tokens[0] == "agl_vibrational_free_energy_300K_cell") {
            data.agl_vibrational_free_energy_300K_cell = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "agl_vibrational_free_energy_300K_atom") {
            data.agl_vibrational_free_energy_300K_atom = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "agl_vibrational_entropy_300K_cell") {
            data.agl_vibrational_entropy_300K_cell = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "agl_vibrational_entropy_300K_atom") {
            data.agl_vibrational_entropy_300K_atom = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
        }
      }
    } else {
      return;
    }

    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_thermal_conductivity_300K (W/m*K) = " << ((data.agl_thermal_conductivity_300K != AUROSTD_NAN) ? aurostd::utype2string(data.agl_thermal_conductivity_300K, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_debye (K) = " << ((data.agl_debye != AUROSTD_NAN) ? aurostd::utype2string(data.agl_debye, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_acoustic_debye (K) = " << ((data.agl_acoustic_debye != AUROSTD_NAN) ? aurostd::utype2string(data.agl_acoustic_debye, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_gruneisen = " << ((data.agl_gruneisen != AUROSTD_NAN) ? aurostd::utype2string(data.agl_gruneisen, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_heat_capacity_Cv_300K (kB/cell) = " << ((data.agl_heat_capacity_Cv_300K != AUROSTD_NAN) ? aurostd::utype2string(data.agl_heat_capacity_Cv_300K, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_heat_capacity_Cp_300K (kB/cell) = " << ((data.agl_heat_capacity_Cp_300K != AUROSTD_NAN) ? aurostd::utype2string(data.agl_heat_capacity_Cp_300K, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_thermal_expansion_300K (1/K) = " << ((data.agl_thermal_expansion_300K != AUROSTD_NAN) ? aurostd::utype2string(data.agl_thermal_expansion_300K, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_bulk_modulus_static_300K (GPa) = " << ((data.agl_bulk_modulus_static_300K != AUROSTD_NAN) ? aurostd::utype2string(data.agl_bulk_modulus_static_300K, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_bulk_modulus_isothermal_300K (GPa) = " << ((data.agl_bulk_modulus_isothermal_300K != AUROSTD_NAN) ? aurostd::utype2string(data.agl_bulk_modulus_isothermal_300K, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_poisson_ratio_source = " << ((!data.agl_poisson_ratio_source.empty()) ? (data.agl_poisson_ratio_source) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_vibrational_free_energy_300K_cell (meV/cell) = "
           << ((data.agl_vibrational_free_energy_300K_cell != AUROSTD_NAN) ? aurostd::utype2string(data.agl_vibrational_free_energy_300K_cell, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_vibrational_free_energy_300K_atom (meV/atom) = "
           << ((data.agl_vibrational_free_energy_300K_atom != AUROSTD_NAN) ? aurostd::utype2string(data.agl_vibrational_free_energy_300K_atom, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_vibrational_entropy_300K_cell (meV/cell*K) = " << ((data.agl_vibrational_entropy_300K_cell != AUROSTD_NAN) ? aurostd::utype2string(data.agl_vibrational_entropy_300K_cell, 10) : "unavailable")
           << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " agl_vibrational_entropy_300K_atom (meV/atom*K) = " << ((data.agl_vibrational_entropy_300K_atom != AUROSTD_NAN) ? aurostd::utype2string(data.agl_vibrational_entropy_300K_atom, 10) : "unavailable")
           << endl; // CT20181212
    }
    // done
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_AEL  // CORMAC
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for AEL results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change siganture\, optimize\, cleanup}
  void LIB2RAW_Loop_AEL(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("ael");

    vector<string> vline;
    vector<string> tokens;
    stringstream aflow_ael_out;

    // AFLOW AEL
    if (aurostd::FileExist(directory_LIB + "/" + "aflow.ael.out") || aurostd::CompressFileExist(directory_LIB + "/" + "aflow.ael.out")) {
      const vector<string> ael_file_list{
          "aflow.ael.out", "AEL_Elastic_constants.out", "AEL_Compliance_tensor.out", "AEL_elastic_tensor.json", "AEL_energy_structures.json",
      };
      for (const auto& file : ael_file_list) {
        aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "aflow.ael.out") << endl;
      }
      // ME20191105
      const string ael_out_str = aurostd::compressfile2string(directory_RAW + "/aflow.ael.out");
      aflow_ael_out.str(aurostd::substring2string(ael_out_str, "[AEL_RESULTS]START", "[AEL_RESULTS]STOP", 1));
      aurostd::stream2vectorstring(aflow_ael_out, vline);
      for (size_t i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "ael_poisson_ratio") {
            data.ael_poisson_ratio = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_bulk_modulus_voigt" || tokens[0] == "ael_bulk_modulus_voight") {
            data.ael_bulk_modulus_voigt = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_bulk_modulus_reuss") {
            data.ael_bulk_modulus_reuss = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_shear_modulus_voigt" || tokens[0] == "ael_shear_modulus_voigth") {
            data.ael_shear_modulus_voigt = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_shear_modulus_reuss") {
            data.ael_shear_modulus_reuss = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_bulk_modulus_vrh") {
            data.ael_bulk_modulus_vrh = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_shear_modulus_vrh") {
            data.ael_shear_modulus_vrh = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "ael_elastic_anisotropy") {
            data.ael_elastic_anisotropy = aurostd::string2utype<double>(tokens[1]); // CO20181129
          }
          if (tokens[0] == "ael_youngs_modulus_vrh") {
            data.ael_youngs_modulus_vrh = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_speed_sound_transverse") {
            data.ael_speed_sound_transverse = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_speed_sound_longitudinal") {
            data.ael_speed_sound_longitudinal = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_speed_sound_average") {
            data.ael_speed_sound_average = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_pughs_modulus_ratio") {
            data.ael_pughs_modulus_ratio = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_debye_temperature") {
            data.ael_debye_temperature = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_applied_pressure") {
            data.ael_applied_pressure = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "ael_average_external_pressure") {
            data.ael_average_external_pressure = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
        }
      }

      // ME20191105 BEGIN
      xmatrix<double> tensor(6, 6);
      vector<double> row;
      aflow_ael_out.str(aurostd::substring2string(ael_out_str, "[AEL_STIFFNESS_TENSOR]START", "[AEL_STIFFNESS_TENSOR]STOP", 1));
      aurostd::stream2vectorstring(aflow_ael_out, vline);
      if (!vline.empty()) {
        vline.pop_back(); // Remove empty line at the end
        try {
          if (vline.size() != 6) {
            stringstream message;
            message << "Could not read stiffness tensor: wrong number of lines"
                    << " (found " << vline.size() << ", need 6).";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          for (int i = 0; i < 6; i++) {
            aurostd::string2tokens(vline[i], row);
            if (row.size() != 6) {
              stringstream message;
              message << "Could not read stiffness tensor."
                      << " Wrong number of columns in line " << (i + 1) << " (found " << row.size() << ", need 6).";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            for (int j = 0; j < 6; j++) {
              tensor[i + 1][j + 1] = row[j];
            }
          }
          data.ael_stiffness_tensor = tensor;
        } catch (aurostd::xerror& e) {
          std::cout << MESSAGE << " ERROR - " << e.buildMessageString() << std::endl;
        }
      } else {
        std::cout << MESSAGE << " WARNING - No stiffness tensor found in aflow.ael.out." << std::endl;
      }

      tensor.clear();
      aflow_ael_out.str(aurostd::substring2string(ael_out_str, "[AEL_COMPLIANCE_TENSOR]START", "[AEL_COMPLIANCE_TENSOR]STOP", 1));
      aurostd::stream2vectorstring(aflow_ael_out, vline);
      if (!vline.empty()) {
        vline.pop_back(); // Remove empty line at the end
        try {
          if (vline.size() != 6) {
            stringstream message;
            message << "Could not read compliance tensor: wrong number of lines"
                    << " (found " << vline.size() << ", need 6).";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
          }
          for (int i = 0; i < 6; i++) {
            aurostd::string2tokens(vline[i], row);
            if (row.size() != 6) {
              stringstream message;
              message << "Could not read compliance tensor:"
                      << " wrong number of columns in line " << (i + 1) << " (found " << row.size() << ", need 6).";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
            }
            for (int j = 0; j < 6; j++) {
              tensor[i + 1][j + 1] = row[j];
            }
          }
          data.ael_compliance_tensor = tensor;
        } catch (aurostd::xerror& e) {
          std::cout << MESSAGE << " ERROR - " << e.buildMessageString() << std::endl;
        }
      } else {
        std::cout << MESSAGE << " WARNING - No compliance tensor found in aflow.ael.out." << std::endl;
      }
      // ME20191105 END
    } else {
      return;
    }

    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_poisson_ratio = " << ((data.ael_poisson_ratio != AUROSTD_NAN) ? aurostd::utype2string(data.ael_poisson_ratio, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_bulk_modulus_voigt (GPa) = " << ((data.ael_bulk_modulus_voigt != AUROSTD_NAN) ? aurostd::utype2string(data.ael_bulk_modulus_voigt, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_bulk_modulus_reuss (GPa) = " << ((data.ael_bulk_modulus_reuss != AUROSTD_NAN) ? aurostd::utype2string(data.ael_bulk_modulus_reuss, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_shear_modulus_voigt (GPa) = " << ((data.ael_shear_modulus_voigt != AUROSTD_NAN) ? aurostd::utype2string(data.ael_shear_modulus_voigt, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_shear_modulus_reuss (GPa) = " << ((data.ael_shear_modulus_reuss != AUROSTD_NAN) ? aurostd::utype2string(data.ael_shear_modulus_reuss, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_bulk_modulus_vrh (GPa) = " << ((data.ael_bulk_modulus_vrh != AUROSTD_NAN) ? aurostd::utype2string(data.ael_bulk_modulus_vrh, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_shear_modulus_vrh (GPa) = " << ((data.ael_shear_modulus_vrh != AUROSTD_NAN) ? aurostd::utype2string(data.ael_shear_modulus_vrh, 10) : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_elastic_anisotropy = " << ((data.ael_elastic_anisotropy != AUROSTD_NAN) ? aurostd::utype2string(data.ael_elastic_anisotropy, 10) : "unavailable") << endl; // CO20181129
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_youngs_modulus_vrh (GPa) = " << ((data.ael_youngs_modulus_vrh != AUROSTD_NAN) ? aurostd::utype2string(data.ael_youngs_modulus_vrh, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_speed_sound_transverse (m/s) = " << ((data.ael_speed_sound_transverse != AUROSTD_NAN) ? aurostd::utype2string(data.ael_speed_sound_transverse, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_speed_sound_longitudinal (m/s) = " << ((data.ael_speed_sound_longitudinal != AUROSTD_NAN) ? aurostd::utype2string(data.ael_speed_sound_longitudinal, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_speed_sound_average (m/s) = " << ((data.ael_speed_sound_average != AUROSTD_NAN) ? aurostd::utype2string(data.ael_speed_sound_average, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_pughs_modulus_ratio = " << ((data.ael_pughs_modulus_ratio != AUROSTD_NAN) ? aurostd::utype2string(data.ael_pughs_modulus_ratio, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_debye_temperature (K) = " << ((data.ael_debye_temperature != AUROSTD_NAN) ? aurostd::utype2string(data.ael_debye_temperature, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_applied_pressure (GPa) = " << ((data.ael_applied_pressure != AUROSTD_NAN) ? aurostd::utype2string(data.ael_applied_pressure, 10) : "unavailable") << endl; // CT20181212
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " ael_average_external_pressure (GPa) = " << ((data.ael_average_external_pressure != AUROSTD_NAN) ? aurostd::utype2string(data.ael_average_external_pressure, 10) : "unavailable") << endl; // CT20181212
    }
    // done
    if (AFLOWLIB_VERBOSE) { // ME20191105
      std::cout << MESSAGE << " ael_stiffness_tensor = ";
      if ((data.ael_stiffness_tensor.rows != 6) || (data.ael_stiffness_tensor.cols != 6)) {
        std::cout << "unavailable" << std::endl;
      } else {
        std::cout << std::endl << data.ael_stiffness_tensor << std::endl;
      }

      std::cout << MESSAGE << " ael_compliance_tensor = ";
      if ((data.ael_compliance_tensor.rows != 6) || (data.ael_compliance_tensor.cols != 6)) {
        std::cout << "unavailable" << std::endl;
      } else {
        std::cout << std::endl << data.ael_compliance_tensor << std::endl;
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// BEGIN ME20210901
// ***************************************************************************
// aflowlib::LIB2RAW_Loop_APL
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for APL results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_APL(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " [1]" << std::endl;
    }
    if (AFLOWLIB_VERBOSE) {
      std::cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << std::endl;
    }

    data.vloop.emplace_back("apl");

    vector<string> vlines;
    vector<string> tokens;
    string file;
    string lines;

    // Always need PHPOSCAR
    file = DEFAULT_APL_PHPOSCAR_FILE;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);

    // aflow.apl.out
    file = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_OUT_FILE;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    file = directory_RAW + "/" + DEFAULT_APL_PHPOSCAR_FILE;
    if (!(aurostd::CompressFileExist(file, file))) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File " + file + " not found", _FILE_NOT_FOUND_);
    }
    lines = aurostd::substring2string(aurostd::compressfile2string(file), "[APL_THERMO_RESULTS]START", "[APL_THERMO_RESULTS]STOP", 1);
    aurostd::string2vectorstring(lines, vlines);
    for (size_t i = 0; i < vlines.size(); i++) {
      aurostd::StringSubstInPlace(vlines[i], "=", " ");
      aurostd::string2tokens(vlines[i], tokens, " ");
      if (tokens.size() >= 2) {
        if (tokens[0] == "energy_free_vibrational_cell_apl_300K") {
          data.energy_free_vibrational_cell_apl_300K = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "energy_free_vibrational_atom_apl_300K") {
          data.energy_free_vibrational_atom_apl_300K = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "entropy_vibrational_cell_apl_300K") {
          data.entropy_vibrational_cell_apl_300K = 1000.0 * KBOLTZEV * aurostd::string2utype<double>(tokens[1]); // Convert to meV/K to be consistent with AGL
        } else if (tokens[0] == "entropy_vibrational_atom_apl_300K") {
          data.entropy_vibrational_atom_apl_300K = 1000.0 * KBOLTZEV * aurostd::string2utype<double>(tokens[1]); // Convert to meV/K to be consistent with AGL
        } else if (tokens[0] == "energy_internal_vibrational_cell_apl_300K") {
          data.energy_internal_vibrational_cell_apl_300K = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "energy_internal_vibrational_atom_apl_300K") {
          data.energy_internal_vibrational_atom_apl_300K = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "energy_zero_point_cell_apl") {
          data.energy_zero_point_cell_apl = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "energy_zero_point_atom_apl") {
          data.energy_zero_point_atom_apl = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "heat_capacity_Cv_cell_apl_300K") {
          data.heat_capacity_Cv_cell_apl_300K = aurostd::string2utype<double>(tokens[1]);
        } else if (tokens[0] == "heat_capacity_Cv_atom_apl_300K") {
          data.heat_capacity_Cv_atom_apl_300K = aurostd::string2utype<double>(tokens[1]);
        }
      }
    }

    if (AFLOWLIB_VERBOSE) {
      std::cout << MESSAGE << " energy_free_vibrational_cell_apl_300K = "
                << ((data.energy_free_vibrational_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_free_vibrational_cell_apl_300K) : "unavailable") << std::endl;
      std::cout << MESSAGE << " energy_free_vibrational_atom_apl_300K = "
                << ((data.energy_free_vibrational_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_free_vibrational_atom_apl_300K) : "unavailable") << std::endl;
      std::cout << MESSAGE << " entropy_vibrational_cell_apl_300K = " << ((data.entropy_vibrational_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.entropy_vibrational_cell_apl_300K) : "unavailable")
                << std::endl;
      std::cout << MESSAGE << " entropy_vibrational_atom_apl_300K = " << ((data.entropy_vibrational_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.entropy_vibrational_atom_apl_300K) : "unavailable")
                << std::endl;
      std::cout << MESSAGE << " energy_internal_vibrational_cell_apl_300K = "
                << ((data.energy_internal_vibrational_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_internal_vibrational_cell_apl_300K) : "unavailable") << std::endl;
      std::cout << MESSAGE << " energy_internal_vibrational_atom_apl_300K = "
                << ((data.energy_internal_vibrational_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_internal_vibrational_atom_apl_300K) : "unavailable") << std::endl;
      std::cout << MESSAGE << " energy_zero_point_cell_apl = " << ((data.energy_zero_point_cell_apl != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_zero_point_cell_apl) : "unavailable") << std::endl;
      std::cout << MESSAGE << " energy_zero_point_atom_apl = " << ((data.energy_zero_point_atom_apl != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_zero_point_atom_apl) : "unavailable") << std::endl;
      std::cout << MESSAGE << " heat_capacity_Cv_cell_apl_300K = " << ((data.heat_capacity_Cv_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.heat_capacity_Cv_cell_apl_300K) : "unavailable")
                << std::endl;
      std::cout << MESSAGE << " heat_capacity_Cv_atom_apl_300K = " << ((data.heat_capacity_Cv_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.heat_capacity_Cv_atom_apl_300K) : "unavailable")
                << std::endl;
    }

    // Thermo json file
    file = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_JSON;
    if (aurostd::CompressFileExist(directory_LIB + "/" + file)) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    }

    // Mean square displacement file
    file = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_MSQRDISP_FILE;
    if (aurostd::CompressFileExist(directory_LIB + "/" + file)) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    }

    // Group velocities file
    file = DEFAULT_APL_FILE_PREFIX + DEFAULT_AAPL_GVEL_FILE;
    if (aurostd::CompressFileExist(directory_LIB + "/" + file)) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
    }

    // Plot dispersions and/or DOS
    const bool plot_disp = (aurostd::CompressFileExist(directory_LIB + "/" + DEFAULT_APL_PHEIGENVAL_FILE) && aurostd::CompressFileExist(directory_LIB + "/" + DEFAULT_APL_PHKPOINTS_FILE));
    const bool plot_dos = aurostd::CompressFileExist(directory_LIB + "/" + DEFAULT_APL_PHDOSCAR_FILE);
    if (plot_disp || plot_dos) {
      if (AFLOWLIB_VERBOSE) {
        std::cout << MESSAGE << " Plotting phonon dispersions and/or DOS." << std::endl;
      }
      if (plot_dos) {
        aflowlib::LIB2RAW_FileNeeded(directory_LIB, DEFAULT_APL_PHDOSCAR_FILE, directory_RAW, DEFAULT_APL_PHDOSCAR_FILE, vfiles, MESSAGE);
      }
      if (plot_disp) {
        aflowlib::LIB2RAW_FileNeeded(directory_LIB, DEFAULT_APL_PHEIGENVAL_FILE, directory_RAW, DEFAULT_APL_PHEIGENVAL_FILE, vfiles, MESSAGE);
        aflowlib::LIB2RAW_FileNeeded(directory_LIB, DEFAULT_APL_PHKPOINTS_FILE, directory_RAW, DEFAULT_APL_PHKPOINTS_FILE, vfiles, MESSAGE);
      }
      aurostd::xoption opts;
      aurostd::xoption plotoptions;
      string plottype;
      if (plot_disp && plot_dos) {
        plottype = "PLOT_PHDISPDOS";
      } else if (plot_disp) {
        plottype = "PLOT_PHDISP";
      } else {
        plottype = "PLOT_PHDOS";
      }
      opts.push_attached(plottype, directory_RAW);
      opts.push_attached("PLOTTER::PRINT", "png");
      if (plot_dos) {
        opts.push_attached("PLOTTER::PROJECTION", "atoms");
      }
      plotoptions = plotter::getPlotOptionsPhonons(opts, plottype);
      if (plot_disp && plot_dos) {
        plotter::PLOT_PHDISPDOS(plotoptions);
        // Convert to json for web
        const xKPOINTS xkpts(directory_LIB + "/" + DEFAULT_APL_PHKPOINTS_FILE);
        const xEIGENVAL xeigen(directory_LIB + "/" + DEFAULT_APL_PHEIGENVAL_FILE);
        const xDOSCAR xdos(directory_LIB + "/" + DEFAULT_APL_PHDOSCAR_FILE);
        xoption jsonoptions;
        jsonoptions.push_attached("DIRECTORY", directory_RAW);
        jsonoptions.flag("NOSHIFT", true);
        ofstream dummy;
        const string json = plotter::bandsDOS2JSON(xdos, xeigen, xkpts, jsonoptions, dummy).toString();
        const string filename = directory_RAW + "/" + data.system_name + "_phdosdata.json";
        aurostd::string2file(json, filename);
      } else if (plot_disp) {
        plotter::PLOT_PHDISP(plotoptions);
      } else {
        plotter::PLOT_DOS(plotoptions);
      }
    }
  }
} // namespace aflowlib
// END ME20210901
// BEGIN AS20200831
// ***************************************************************************
// aflowlib::LIB2RAW_Loop_QHA  // SMOLYANYUK
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for QHA results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_QHA(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("qha");
    vector<string> vline;
    vector<string> tokens;
    stringstream aflow_qha_out;
    string qha_raw_file = directory_RAW + "/" + DEFAULT_QHA_FILE_PREFIX + "out";

    if (aurostd::CompressFileExist(directory_LIB + "/" + DEFAULT_QHA_FILE_PREFIX + "out")) { //"aflow.qha.out"
      const vector<string> qha_file_list{"aflow_qha.in",
                                         DEFAULT_QHA_FILE_PREFIX + "out",
                                         DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_THERMO_FILE,
                                         DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_THERMO_FILE,
                                         DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_FVT_FILE,
                                         DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_PDIS_FILE + ".T300K.out",
                                         DEFAULT_APL_PHPOSCAR_FILE,
                                         DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_KPOINTS_FILE};
      for (const auto& file : qha_file_list) {
        aflowlib::LIB2RAW_FileNeeded(directory_LIB, file, directory_RAW, file, vfiles, MESSAGE);
      }

      // read QHA data from the aflow.qha.out file
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(qha_raw_file) << endl;
      }
      if (!(aurostd::CompressFileExist(qha_raw_file, qha_raw_file))) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File " + qha_raw_file + " not found", _FILE_NOT_FOUND_);
      }

      aflow_qha_out.str(aurostd::substring2string(aurostd::compressfile2string(qha_raw_file), "[QHA_RESULTS]START", "[QHA_RESULTS]STOP", 1));
      aurostd::stream2vectorstring(aflow_qha_out, vline);
      for (size_t i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "gruneisen_qha") {
            data.gruneisen_qha = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "gruneisen_qha_300K") {
            data.gruneisen_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "thermal_expansion_qha_300K") {
            data.thermal_expansion_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "modulus_bulk_qha_300K") {
            data.modulus_bulk_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "modulus_bulk_derivative_pressure_qha_300K") {
            data.modulus_bulk_derivative_pressure_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "heat_capacity_Cv_atom_qha_300K") {
            data.heat_capacity_Cv_atom_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "heat_capacity_Cv_cell_qha_300K") {
            data.heat_capacity_Cv_cell_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "heat_capacity_Cp_atom_qha_300K") {
            data.heat_capacity_Cp_atom_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "heat_capacity_Cp_cell_qha_300K") {
            data.heat_capacity_Cp_cell_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "volume_atom_qha_300K") {
            data.volume_atom_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "energy_free_atom_qha_300K") {
            data.energy_free_atom_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "energy_free_cell_qha_300K") {
            data.energy_free_cell_qha_300K = aurostd::string2utype<double>(tokens[1]);
          }
        }
      }

      // plot thermodynamic
      if (aurostd::CompressFileExist(directory_LIB + "/" + DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_THERMO_FILE)) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " plotting QHA thermodynamic data " << endl;
        }
        aurostd::xoption opt;
        opt.flag("PLOT_THERMO_QHA", true);
        opt.addattachedscheme("PLOT_THERMO_QHA", directory_RAW, true);
        opt.push_attached("PLOTTER::PRINT", "png");
        aurostd::xoption plotopts = plotter::getPlotOptionsQHAthermo(opt, "PLOT_THERMO_QHA"); // CO+AS20210713
        plotter::PLOT_THERMO_QHA(plotopts);
      }

      // convert T-dependent phonon dispersions to JSON file
      if (aurostd::CompressFileExist(directory_RAW + "/" + DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_PDIS_FILE + ".T300K.out") && aurostd::CompressFileExist(directory_RAW + "/" + DEFAULT_APL_PHPOSCAR_FILE) &&
          aurostd::CompressFileExist(directory_RAW + "/" + DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_KPOINTS_FILE)) {
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " converting T-dependent phonon dispersions to JSON format " << endl;
        }
        stringstream json;
        const xstructure xstr(directory_RAW + "/" + DEFAULT_APL_PHPOSCAR_FILE);
        const xKPOINTS xkpts(directory_RAW + "/" + DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_KPOINTS_FILE);
        const xEIGENVAL xeig(directory_RAW + "/" + DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_PDIS_FILE + ".T300K.out");
        xoption xopt;
        xopt.push_attached("EFERMI", "0.0");
        xopt.push_attached("OUTPUT_FORMAT", "JSON");
        xopt.push_attached("DIRECTORY", directory_RAW);
        xopt.flag("NOSHIFT", true);
        plotter::generateBandPlot(json, xeig, xkpts, xstr, xopt);
        aurostd::stringstream2file(json, directory_RAW + "/" + DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_PDIS_FILE + ".T300K.json");
      }
    } else {
      return;
    }

    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " gruneisen_qha = " << ((data.gruneisen_qha != AUROSTD_NAN) ? aurostd::utype2string(data.gruneisen_qha, 10) : "unavailable") << endl;
      cout << MESSAGE << " gruneisen_qha_300K = " << ((data.gruneisen_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.gruneisen_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " thermal_expansion_qha_300K (10^-5/K) = " << ((data.thermal_expansion_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.thermal_expansion_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " modulus_bulk_qha_300K (GPa) = " << ((data.modulus_bulk_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.modulus_bulk_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " modulus_bulk_derivative_pressure_qha_300K = "
           << ((data.modulus_bulk_derivative_pressure_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.modulus_bulk_derivative_pressure_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " heat_capacity_Cv_atom_qha_300K = " << ((data.heat_capacity_Cv_atom_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.heat_capacity_Cv_atom_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " heat_capacity_Cv_cell_qha_300K = " << ((data.heat_capacity_Cv_cell_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.heat_capacity_Cv_cell_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " heat_capacity_Cp_atom_qha_300K = " << ((data.heat_capacity_Cp_atom_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.heat_capacity_Cp_atom_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " heat_capacity_Cp_cell_qha_300K = " << ((data.heat_capacity_Cp_cell_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.heat_capacity_Cp_cell_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " volume_atom_qha_300K = " << ((data.volume_atom_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.volume_atom_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " energy_free_atom_qha_300K = " << ((data.energy_free_atom_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.energy_free_atom_qha_300K, 10) : "unavailable") << endl;
      cout << MESSAGE << " energy_free_cell_qha_300K = " << ((data.energy_free_cell_qha_300K != AUROSTD_NAN) ? aurostd::utype2string(data.energy_free_cell_qha_300K, 10) : "unavailable") << endl;
    }
    // done
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib
// END AS20200831

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_POCC - CO20200624
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for POCC results
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_POCC(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    stringstream message;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("pocc");

    // CMo_pvNb_svTa_pvV_svW_pv:PAW_PBE.AB_cF8_225_a_b.AB:POCC_P0-1xA_P1-0.2xB-0.2xC-0.2xD-0.2xE-0.2xF
    data.pocc_parameters = data.system_name;
    string::size_type loc;
    loc = data.pocc_parameters.find(TAG_TITLE_POCC);
    data.pocc_parameters = data.pocc_parameters.substr(loc + TAG_TITLE_POCC.size(), string::npos);
    if (data.pocc_parameters.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No pocc_parameters found", _INPUT_ILLEGAL_);
    }
    data.pocc_parameters = pocc::addDefaultPOCCTOL2string(data.pocc_parameters);
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " pocc_parameters=" << data.pocc_parameters << endl;
    }

    vector<string> vline;
    vector<string> tokens;
    stringstream aflow_pocc_out;
    stringstream aflow_pocc_agl_out;
    uint i = 0;
    uint j = 0;
    vector<double> v_dgs;
    vector<double> v_energies;
    vector<string> vdirfiles;
    string pocc_out_lib_file = directory_LIB + "/" + POCC_FILE_PREFIX + POCC_OUT_FILE;
    string pocc_agl_lib_file = directory_LIB + "/" + "aflow.pocc_agl.out";
    string pocc_agl_raw_file = directory_RAW + "/" + "aflow.pocc_agl.out";
    string pocc_agl_raw_file_content;

    if (aurostd::CompressFileExist(pocc_out_lib_file, pocc_out_lib_file)) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "aflow.in", directory_RAW, "aflow.in", vfiles, MESSAGE); // aflow.in, needed for ExtractSystemName() in plotter
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, POCC_FILE_PREFIX + POCC_OUT_FILE, directory_RAW, POCC_FILE_PREFIX + POCC_OUT_FILE, vfiles, MESSAGE); // aflow.pocc.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, POCC_FILE_PREFIX + POCC_UNIQUE_SUPERCELLS_FILE, directory_RAW, POCC_FILE_PREFIX + POCC_UNIQUE_SUPERCELLS_FILE, vfiles, MESSAGE); // aflow.pocc.structures_unique.out
      // get DOSCAR.pocc + png's
      aurostd::DirectoryLS(directory_LIB, vdirfiles);
      for (i = 0; i < vdirfiles.size(); i++) {
        if (vdirfiles[i].find(POCC_DOSCAR_FILE) != string::npos) {
          cout << directory_LIB + "/" + vdirfiles[i] << endl;
          aurostd::LinkFile(directory_LIB + "/" + vdirfiles[i], directory_RAW); // link the file, no need to de-compress
        }
        const vector<string> doss{"_dos_orbitals_", "_dos_species_", "_dos_atoms_", "_phdos_"};
        for (const auto& dos : doss) {
          if ((vdirfiles[i].find(dos) != string::npos) && (vdirfiles[i].find(".png") != string::npos)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB, vdirfiles[i], directory_RAW, vdirfiles[i], vfiles, MESSAGE);
          }
        }
      }
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_LIB + "/" + POCC_FILE_PREFIX + POCC_OUT_FILE) << endl;
      }

      aflow_pocc_out.str(aurostd::substring2string(aurostd::compressfile2string(pocc_out_lib_file), "[AFLOW_POCC]START_TEMPERATURE=ALL", "[AFLOW_POCC]STOP_TEMPERATURE=ALL", 1));
      aurostd::stream2vectorstring(aflow_pocc_out, vline);
      for (i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "enthalpy_mix_atom") {
            data.energy_atom = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "entropy_forming_ability") {
            data.entropy_forming_ability = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0].find("degeneracy_supercell_") != string::npos) {
            v_dgs.emplace_back(aurostd::string2utype<double>(tokens[1]));
          }
          if (tokens[0].find("enthalpy_atom_supercell_") != string::npos) {
            v_energies.emplace_back(aurostd::string2utype<double>(tokens[1]));
          }
        }
      }
      if (AFLOWLIB_VERBOSE && data.energy_atom != AUROSTD_NAN) {
        cout << MESSAGE << " energy_atom=" << data.energy_atom << endl;
      }
      if (AFLOWLIB_VERBOSE && data.entropy_forming_ability != AUROSTD_NAN) {
        cout << MESSAGE << " entropy_forming_ability=" << data.entropy_forming_ability << endl;
      }
    }

    if (v_dgs.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No v_dgs found", _FILE_CORRUPT_);
    }
    if (v_energies.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No v_energies found", _FILE_CORRUPT_);
    }
    if (v_dgs.size() != v_energies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_dgs.size()!=v_energies.size()", _FILE_CORRUPT_);
    }

    if (aurostd::CompressFileExist(directory_LIB + "/" + POCC_FILE_PREFIX + POCC_ALL_HNF_MATRICES_FILE)) { // old pocc doesn't print this
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, POCC_FILE_PREFIX + POCC_ALL_HNF_MATRICES_FILE, directory_RAW, POCC_FILE_PREFIX + POCC_ALL_HNF_MATRICES_FILE, vfiles, MESSAGE); // aflow.pocc.hnf_matrices.out
    }
    if (aurostd::CompressFileExist(directory_LIB + "/" + POCC_FILE_PREFIX + POCC_ALL_SITE_CONFIGURATIONS_FILE)) { // old pocc doesn't print this
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, POCC_FILE_PREFIX + POCC_ALL_SITE_CONFIGURATIONS_FILE, directory_RAW, POCC_FILE_PREFIX + POCC_ALL_SITE_CONFIGURATIONS_FILE, vfiles, MESSAGE); // aflow.pocc.site_configurations.out
    }

    if (aurostd::CompressFileExist(pocc_agl_lib_file, pocc_agl_lib_file)) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, "aflow.pocc_agl.out", directory_RAW, "aflow.pocc_agl.out", vfiles, MESSAGE); // aflow.pocc_agl.out
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " loading " << string(directory_RAW + "/" + "aflow.pocc_agl.out") << endl;
      }
      // 0K START
      if (!(aurostd::CompressFileExist(pocc_agl_raw_file, pocc_agl_raw_file))) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File " + pocc_agl_raw_file + " not found", _FILE_NOT_FOUND_);
      }
      aurostd::compressfile2string(pocc_agl_raw_file, pocc_agl_raw_file_content);
      aflow_pocc_agl_out.str(aurostd::substring2string(pocc_agl_raw_file_content, "[POCC_AGL_RESULTS]START_TEMPERATURE=0000_K", "[POCC_AGL_RESULTS]STOP_TEMPERATURE=0000_K", 1));
      aurostd::stream2vectorstring(aflow_pocc_agl_out, vline);
      for (i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "pocc_agl_debye") {
            data.agl_debye = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_acoustic_debye") {
            data.agl_acoustic_debye = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_gruneisen") {
            data.agl_gruneisen = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_poisson_ratio_source") {
            data.agl_poisson_ratio_source = tokens[1]; // CT20181212
          }
          // what to do about:
          // pocc_agl_gibbs_energy_atom_ave
          // pocc_agl_vibrational_energy_atom_ave
        }
      }
      if (AFLOWLIB_VERBOSE && data.agl_debye != AUROSTD_NAN) {
        cout << MESSAGE << " agl_debye=" << data.agl_debye << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_acoustic_debye != AUROSTD_NAN) {
        cout << MESSAGE << " agl_acoustic_debye=" << data.agl_acoustic_debye << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_gruneisen != AUROSTD_NAN) {
        cout << MESSAGE << " agl_gruneisen=" << data.agl_gruneisen << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.agl_poisson_ratio_source.empty()) {
        cout << MESSAGE << " agl_poisson_ratio_source=" << data.agl_poisson_ratio_source << endl;
      }
      // 0K STOP
      // 300K START
      aflow_pocc_agl_out.str(aurostd::substring2string(pocc_agl_raw_file_content, "[POCC_AGL_RESULTS]START_TEMPERATURE=0300_K", "[POCC_AGL_RESULTS]STOP_TEMPERATURE=0300_K", 1));
      aurostd::stream2vectorstring(aflow_pocc_agl_out, vline);
      for (i = 0; i < vline.size(); i++) {
        aurostd::StringSubstInPlace(vline[i], "=", " ");
        aurostd::string2tokens(vline[i], tokens, " ");
        if (tokens.size() >= 2) {
          if (tokens[0] == "pocc_agl_heat_capacity_Cv_300K") {
            data.agl_heat_capacity_Cv_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_heat_capacity_Cp_300K") {
            data.agl_heat_capacity_Cp_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_thermal_expansion_300K") {
            data.agl_thermal_expansion_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_bulk_modulus_static_300K") {
            data.agl_bulk_modulus_static_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_bulk_modulus_isothermal_300K") {
            data.agl_bulk_modulus_isothermal_300K = aurostd::string2utype<double>(tokens[1]);
          }
          if (tokens[0] == "pocc_agl_vibrational_free_energy_300K_cell") {
            data.agl_vibrational_free_energy_300K_cell = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "pocc_agl_vibrational_free_energy_300K_atom") {
            data.agl_vibrational_free_energy_300K_atom = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "pocc_agl_vibrational_entropy_300K_cell") {
            data.agl_vibrational_entropy_300K_cell = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "pocc_agl_vibrational_entropy_300K_atom") {
            data.agl_vibrational_entropy_300K_atom = aurostd::string2utype<double>(tokens[1]); // CT20181212
          }
          if (tokens[0] == "pocc_agl_thermal_conductivity_300K") {
            data.agl_thermal_conductivity_300K = aurostd::string2utype<double>(tokens[1]);
          }
          // what to do about:
          // pocc_agl_gibbs_energy_atom_ave
          // pocc_agl_vibrational_energy_atom_ave
        }
      }
      if (AFLOWLIB_VERBOSE && data.agl_heat_capacity_Cv_300K != AUROSTD_NAN) {
        cout << MESSAGE << " agl_heat_capacity_Cv_300K=" << data.agl_heat_capacity_Cv_300K << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_heat_capacity_Cp_300K != AUROSTD_NAN) {
        cout << MESSAGE << " agl_heat_capacity_Cp_300K=" << data.agl_heat_capacity_Cp_300K << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_thermal_expansion_300K != AUROSTD_NAN) {
        cout << MESSAGE << " agl_thermal_expansion_300K=" << data.agl_thermal_expansion_300K << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_bulk_modulus_static_300K != AUROSTD_NAN) {
        cout << MESSAGE << " agl_bulk_modulus_static_300K=" << data.agl_bulk_modulus_static_300K << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_bulk_modulus_isothermal_300K != AUROSTD_NAN) {
        cout << MESSAGE << " agl_bulk_modulus_isothermal_300K=" << data.agl_bulk_modulus_isothermal_300K << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_vibrational_free_energy_300K_cell != AUROSTD_NAN) {
        cout << MESSAGE << " agl_vibrational_free_energy_300K_cell=" << data.agl_vibrational_free_energy_300K_cell << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_vibrational_free_energy_300K_atom != AUROSTD_NAN) {
        cout << MESSAGE << " agl_vibrational_free_energy_300K_atom=" << data.agl_vibrational_free_energy_300K_atom << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_vibrational_entropy_300K_cell != AUROSTD_NAN) {
        cout << MESSAGE << " agl_vibrational_entropy_300K_cell=" << data.agl_vibrational_entropy_300K_cell << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_vibrational_entropy_300K_atom != AUROSTD_NAN) {
        cout << MESSAGE << " agl_vibrational_entropy_300K_atom=" << data.agl_vibrational_entropy_300K_atom << endl;
      }
      if (AFLOWLIB_VERBOSE && data.agl_thermal_conductivity_300K != AUROSTD_NAN) {
        cout << MESSAGE << " agl_thermal_conductivity_300K=" << data.agl_thermal_conductivity_300K << endl;
      }
      // 300K STOP
    }

    // ME20210927 - APL
    const string aplout = POCC_FILE_PREFIX + POCC_APL_OUT_FILE;
    string apl_out_raw_file = directory_RAW + "/" + aplout;
    if (aurostd::CompressFileExist(directory_LIB + "/" + aplout)) {
      // Link PHPOSCAR for plotting
      for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
        if (aurostd::FileExist(directory_LIB + "/" + DEFAULT_APL_PHPOSCAR_FILE + XHOST.vext[iext])) {
          aurostd::LinkFile(directory_LIB + "/" + DEFAULT_APL_PHPOSCAR_FILE + XHOST.vext[iext], directory_RAW);
          break;
        }
      }
      aflowlib::LIB2RAW_FileNeeded(directory_LIB, aplout, directory_RAW, aplout, vfiles, MESSAGE);
      if (AFLOWLIB_VERBOSE) {
        std::cout << MESSAGE << " loading " << directory_RAW << "/" << aplout << std::endl;
      }
      if (!(aurostd::CompressFileExist(apl_out_raw_file, apl_out_raw_file))) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "File " + apl_out_raw_file + " not found", _FILE_NOT_FOUND_);
      }
      const string lines = aurostd::substring2string(aurostd::compressfile2string(apl_out_raw_file), "[POCC_APL_RESULTS]START_TEMPERATURE=0300_K", "[POCC_APL_RESULTS]STOP_TEMPERATURE=0300_K", 1);
      if (!lines.empty()) {
        const string lines_300K = aurostd::substring2string(lines, "[APL_THERMO_RESULTS]START", "[APL_THERMO_RESULTS]STOP", 1);
        vector<string> vlines;
        aurostd::string2vectorstring(lines_300K, vlines);
        for (size_t k = 0; k < vlines.size(); k++) {
          aurostd::StringSubstInPlace(vlines[k], "=", " ");
          aurostd::string2tokens(vlines[k], tokens, " ");
          if (tokens.size() >= 2) {
            if (tokens[0] == "energy_free_vibrational_cell_apl_300K") {
              data.energy_free_vibrational_cell_apl_300K = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "energy_free_vibrational_atom_apl_300K") {
              data.energy_free_vibrational_atom_apl_300K = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "entropy_vibrational_cell_apl_300K") {
              data.entropy_vibrational_cell_apl_300K = 1000.0 * KBOLTZEV * aurostd::string2utype<double>(tokens[1]); // Convert to meV/K to be consistent with AGL
            } else if (tokens[0] == "entropy_vibrational_atom_apl_300K") {
              data.entropy_vibrational_atom_apl_300K = 1000.0 * KBOLTZEV * aurostd::string2utype<double>(tokens[1]); // Convert to meV/K to be consistent with AGL
            } else if (tokens[0] == "energy_internal_vibrational_cell_apl_300K") {
              data.energy_internal_vibrational_cell_apl_300K = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "energy_internal_vibrational_atom_apl_300K") {
              data.energy_internal_vibrational_atom_apl_300K = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "energy_zero_point_cell_apl") {
              data.energy_zero_point_cell_apl = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "energy_zero_point_atom_apl") {
              data.energy_zero_point_atom_apl = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "heat_capacity_Cv_cell_apl_300K") {
              data.heat_capacity_Cv_cell_apl_300K = aurostd::string2utype<double>(tokens[1]);
            } else if (tokens[0] == "heat_capacity_Cv_atom_apl_300K") {
              data.heat_capacity_Cv_atom_apl_300K = aurostd::string2utype<double>(tokens[1]);
            }
          }
        }
      }
      if (AFLOWLIB_VERBOSE) {
        std::cout << MESSAGE << " energy_free_vibrational_cell_apl_300K = "
                  << ((data.energy_free_vibrational_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_free_vibrational_cell_apl_300K) : "unavailable") << std::endl;
        std::cout << MESSAGE << " energy_free_vibrational_atom_apl_300K = "
                  << ((data.energy_free_vibrational_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_free_vibrational_atom_apl_300K) : "unavailable") << std::endl;
        std::cout << MESSAGE << " entropy_vibrational_cell_apl_300K = " << ((data.entropy_vibrational_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.entropy_vibrational_cell_apl_300K) : "unavailable")
                  << std::endl;
        std::cout << MESSAGE << " entropy_vibrational_atom_apl_300K = " << ((data.entropy_vibrational_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.entropy_vibrational_atom_apl_300K) : "unavailable")
                  << std::endl;
        std::cout << MESSAGE << " energy_internal_vibrational_cell_apl_300K = "
                  << ((data.energy_internal_vibrational_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_internal_vibrational_cell_apl_300K) : "unavailable") << std::endl;
        std::cout << MESSAGE << " energy_internal_vibrational_atom_apl_300K = "
                  << ((data.energy_internal_vibrational_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_internal_vibrational_atom_apl_300K) : "unavailable") << std::endl;
        std::cout << MESSAGE << " energy_zero_point_cell_apl = " << ((data.energy_zero_point_cell_apl != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_zero_point_cell_apl) : "unavailable") << std::endl;
        std::cout << MESSAGE << " energy_zero_point_atom_apl = " << ((data.energy_zero_point_atom_apl != AUROSTD_NAN) ? aurostd::utype2string<double>(data.energy_zero_point_atom_apl) : "unavailable") << std::endl;
        std::cout << MESSAGE << " heat_capacity_Cv_cell_apl_300K = " << ((data.heat_capacity_Cv_cell_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.heat_capacity_Cv_cell_apl_300K) : "unavailable")
                  << std::endl;
        std::cout << MESSAGE << " heat_capacity_Cv_atom_apl_300K = " << ((data.heat_capacity_Cv_atom_apl_300K != AUROSTD_NAN) ? aurostd::utype2string<double>(data.heat_capacity_Cv_atom_apl_300K) : "unavailable")
                  << std::endl;
      }
    }

    string AflowIn_file;
    string AflowIn;
    KBIN::getAflowInFromDirectory(directory_LIB, AflowIn_file, AflowIn);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " loaded aflow.in from dir=" << directory_LIB << endl;
    }
    xstructure xstr_pocc = pocc::extractPARTCAR(AflowIn);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " loaded PARTCAR" << endl;
      cerr << xstr_pocc << endl;
    }
    //
    data.volume_cell_orig = xstr_pocc.GetVolume();
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " volume_cell_orig=" << data.volume_cell_orig << endl;
    }
    //
    xvector<double> data_abcabc;
    data_abcabc = Getabc_angles(xstr_pocc.scale * xstr_pocc.lattice, DEGREES);
    data.vgeometry_orig = aurostd::xvector2vector(data_abcabc);
    data.geometry_orig = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vgeometry_orig, _AFLOWLIB_DATA_GEOMETRY_PREC_), ",");
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " geometry_orig=" << data.geometry_orig << endl;
    }
    //
    data.compound = getGenericTitleXStructure(xstr_pocc, false);
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " compound=" << data.compound << endl;
    }
    //
    data.vstoichiometry = aurostd::deque2vector(xstr_pocc.stoich_each_type);
    data.stoichiometry = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vstoichiometry, _AFLOWLIB_STOICH_PRECISION_), ",");
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " stoichiometry=" << data.stoichiometry << endl;
    }
    // CO20200624 START - mimic stoich from PrintData1(): aflow_pflow_print.cpp, this is really obsolete
    stringstream stoich_ss;
    stoich_ss.precision(4);
    for (i = 0; i < xstr_pocc.stoich_each_type.size(); i++) {
      stoich_ss << setw(8) << xstr_pocc.stoich_each_type[i] << " ";
    }
    data.stoich = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(stoich_ss.str());
    // CO20200624 END - mimic stoich from PrintData1(): aflow_pflow_print.cpp, this is really obsolete

    // load properties for ARUN.POCC_0
    // only properties that MUST be true for all systems
    // this does NOT include cell properties
    // remember: pocc right now may have cell sizes the same, but this may not always be true
    _aflags aflags;
    aflags.Directory = directory_LIB;
    pocc::POccCalculator pcalc(aflags);
    pcalc.loadDataIntoCalculator();
    pcalc.setTemperatureStringParameters(); // needed for DOSCAR plots
    if (pcalc.m_ARUN_directories.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No ARUN.POCC_* runs found", _FILE_CORRUPT_);
    }
    if (pcalc.m_ARUN_directories.size() != v_dgs.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "pcalc.m_ARUN_directories.size()!=v_dgs.size()", _FILE_CORRUPT_);
    }
    if (pcalc.m_ARUN_directories.size() != v_energies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "pcalc.m_ARUN_directories.size()!=v_energies.size()", _FILE_CORRUPT_);
    }
    const xvector<double> xv_dgs = aurostd::vector2xvector(v_dgs); // for aurostd::meanWeighted()

    // do DOSCAR plots
    vdirfiles.clear();
    aurostd::DirectoryLS(directory_RAW, vdirfiles);
    std::sort(vdirfiles.begin(), vdirfiles.end()); // get in order
    for (i = 0; i < vdirfiles.size(); i++) {
      if ((vdirfiles[i].find(POCC_DOSCAR_FILE) != string::npos) || (vdirfiles[i].find(POCC_PHDOSCAR_FILE) != string::npos)) { // ME20210927 - added PHDOSCAR
        // need to grab POSCAR from ARUN.POCC_01
        // inside plotter we change '/RAW/' to '/LIB/', everything in RAW must be self-contained
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " plotting " << vdirfiles[i] << endl;
        }
        pcalc.plotAvgDOSCAR(directory_RAW + "/" + vdirfiles[i], directory_RAW);
      }
    }

    // xOUTCAR
    string filename;
    if (!aurostd::CompressFileExist(directory_LIB + "/" + pcalc.m_ARUN_directories[0] + "/OUTCAR.relax2", filename)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUTCAR cannot be extracted", _FILE_CORRUPT_);
    }
    xOUTCAR xOUT;
    xOUT.GetPropertiesFile(filename);
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xOUT.species.size()=" << xOUT.species.size() << endl;
    }
    if (xOUT.species.size() != data.vspecies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUT.species.size()!=data.vspecies.size()", _FILE_CORRUPT_);
    }
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xOUT.species_pp.size()=" << xOUT.species_pp.size() << endl;
    }
    if (xOUT.species_pp.size() != data.vspecies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUT.species_pp.size()!=data.vspecies.size()", _FILE_CORRUPT_);
    }
    if (xOUT.species_pp.size() != xstr_pocc.species.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUT.species_pp.size()!=xstr_pocc.species.size()", _FILE_CORRUPT_);
    }
    //
    data.vspecies_pp = xOUT.species_pp; // for aflowlib_libraries.cpp
    xstr_pocc.species_pp.assign(xOUT.species_pp.begin(), xOUT.species_pp.end()); // for LIB2RAW_Calculate_FormationEnthalpy
    data.species_pp = aurostd::joinWDelimiter(data.vspecies_pp, ",");
    if (AFLOWLIB_VERBOSE && !data.species_pp.empty()) {
      cout << MESSAGE << " species_pp=" << data.species_pp << endl;
    }
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xOUT.species_pp_version.size()=" << xOUT.species_pp_version.size() << endl;
    }
    if (xOUT.species_pp_version.size() != data.vspecies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUT.species_pp_version.size()!=data.vspecies.size()", _FILE_CORRUPT_);
    }
    data.vspecies_pp_version = xOUT.species_pp_version; // for aflowlib_libraries.cpp
    xstr_pocc.species_pp_version.assign(xOUT.species_pp_version.begin(), xOUT.species_pp_version.end()); // for LIB2RAW_Calculate_FormationEnthalpy
    data.species_pp_version = aurostd::joinWDelimiter(data.vspecies_pp_version, ",");
    if (AFLOWLIB_VERBOSE && !data.species_pp_version.empty()) {
      cout << MESSAGE << " species_pp_version=" << data.species_pp_version << endl;
    }
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xOUT.vZVAL.size()=" << xOUT.vZVAL.size() << endl;
    }
    if (xOUT.vZVAL.size() != data.vspecies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUT.vZVAL.size()!=data.vspecies.size()", _FILE_CORRUPT_);
    }
    data.vspecies_pp_ZVAL = xOUT.vZVAL; // for aflowlib_libraries.cpp
    xstr_pocc.species_pp_ZVAL.assign(xOUT.vZVAL.begin(), xOUT.vZVAL.end()); // for LIB2RAW_Calculate_FormationEnthalpy
    data.species_pp_ZVAL = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vspecies_pp_ZVAL), ",");
    if (AFLOWLIB_VERBOSE && !data.species_pp_ZVAL.empty()) {
      cout << MESSAGE << " species_pp_ZVAL=" << data.species_pp_ZVAL << endl;
    }
    //
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xOUT.species_pp_AUID.size()=" << xOUT.species_pp_AUID.size() << endl;
    }
    if (xOUT.species_pp_AUID.size() != data.vspecies.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xOUT.species_pp_AUID.size()!=data.vspecies.size()", _FILE_CORRUPT_);
    }
    data.vspecies_pp_AUID = xOUT.species_pp_AUID; // for aflowlib_libraries.cpp
    data.species_pp_AUID = aurostd::joinWDelimiter(data.vspecies_pp_AUID, ",");
    if (AFLOWLIB_VERBOSE && !data.species_pp_AUID.empty()) {
      cout << MESSAGE << " species_pp_AUID=" << data.species_pp_AUID << endl;
    }
    //
    data.vdft_type = {xOUT.pp_type}; // CO, this is technically a vector (RESTAPI paper)
    data.dft_type = aurostd::joinWDelimiter(data.vdft_type, ",");
    if (AFLOWLIB_VERBOSE && !data.dft_type.empty()) {
      cout << MESSAGE << " dft_type=" << data.dft_type << endl;
    }
    data.ldau_TLUJ = xOUT.string_LDAU;
    //[CO+ME20210713 - keep legacy behavior, only print when non-zero]if(data.ldau_TLUJ.empty()){data.ldau_TLUJ=aurostd::utype2string(0);} //CO20210713 - no +U
    if (AFLOWLIB_VERBOSE && !data.ldau_TLUJ.empty()) {
      cout << MESSAGE << " ldau_TLUJ=" << data.ldau_TLUJ << endl;
    }
    // ME20190124 BEGIN - Store LDAU information individually
    //  Note that the vector here has the species in the columns, not the
    //  rows because this is closer to the format in the out and json files.
    xstr_pocc.species_pp_vLDAU.clear();
    if (!xOUT.species_pp_vLDAU.empty()) {
      data.vLDAU.resize(xOUT.species_pp_vLDAU[0].size());
    } // CO20200731
    else { // CO20210713 - set ldau_type=0
      data.vLDAU.resize(4);
      for (size_t k = 0; k < xOUT.species.size(); k++) {
        data.vLDAU[0].push_back(0);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xOUT.species_pp_vLDAU.size()=" << xOUT.species_pp_vLDAU.size() << endl;
      for (i = 0; i < xOUT.species_pp_vLDAU.size(); i++) {
        cerr << __AFLOW_FUNC__ << " xOUT.species_pp_vLDAU[i=" << i << "].size()=" << xOUT.species_pp_vLDAU[i].size() << endl;
        for (j = 0; j < xOUT.species_pp_vLDAU[i].size(); j++) {
          cerr << __AFLOW_FUNC__ << " xOUT.species_pp_vLDAU[i=" << i << "][j=" << j << "]=" << xOUT.species_pp_vLDAU[i][j] << endl;
        }
      }
    }
    for (i = 0; i < xOUT.species_pp_vLDAU.size(); i++) {
      xstr_pocc.species_pp_vLDAU.emplace_back(xOUT.species_pp_vLDAU[i]); // keep the same structure // for LIB2RAW_Calculate_FormationEnthalpy
      for (j = 0; j < xOUT.species_pp_vLDAU[i].size(); j++) {
        data.vLDAU[j].emplace_back(xOUT.species_pp_vLDAU[i][j]);
      }
    }
    // ME20190124 END
    data.METAGGA = xOUT.METAGGA;
    if (AFLOWLIB_VERBOSE && !data.METAGGA.empty()) {
      cout << MESSAGE << " METAGGA=" << data.METAGGA << endl;
    }
    data.energy_cutoff = xOUT.ENCUT;
    if (AFLOWLIB_VERBOSE && data.energy_cutoff != AUROSTD_NAN) {
      cout << MESSAGE << " energy_cutoff=" << data.energy_cutoff << endl;
    }
    if (std::abs(xOUT.PV_atom) < PRESSURE_ZERO_ENTHALPY_ENERGY) { // only set if it's zero because that means ALL ARUNs are zero (as input), otherwise it's ARUN-specific
      data.pressure = xOUT.pressure;
      if (AFLOWLIB_VERBOSE && data.pressure != AUROSTD_NAN) {
        cout << MESSAGE << " pressure=" << data.pressure << endl;
      }
      data.PV_cell = xOUT.PV_cell;
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " PV total E0 (eV) = " << data.PV_cell << endl;
      }
      data.PV_atom = xOUT.PV_atom;
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " PV per atom E0/N (eV) = " << data.PV_atom << endl;
      }
      if (data.energy_atom != AUROSTD_NAN) {
        data.enthalpy_atom = data.energy_atom;
      }
      if (AFLOWLIB_VERBOSE && data.enthalpy_atom != AUROSTD_NAN) {
        cout << MESSAGE << " enthalpy_atom=" << data.enthalpy_atom << endl;
      }
    } // do else later, see aflow_ovasp
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " data.nspecies=" << data.nspecies;
      cerr << " xstr_pocc.species.size()=" << xstr_pocc.species.size();
      cerr << " xstr_pocc.species_pp.size()=" << xstr_pocc.species_pp.size();
      cerr << " xstr_pocc.species_pp_type.size()=" << xstr_pocc.species_pp_type.size();
      cerr << " xstr_pocc.species_pp_version.size()=" << xstr_pocc.species_pp_version.size();
      cerr << " xstr_pocc.species_pp_ZVAL.size()=" << xstr_pocc.species_pp_ZVAL.size();
      cerr << " xstr_pocc.species_pp_vLDAU.size()=" << xstr_pocc.species_pp_vLDAU.size();
      cerr << " xstr_pocc.species_volume.size()=" << xstr_pocc.species_volume.size();
      cerr << " xstr_pocc.species_mass.size()=" << xstr_pocc.species_mass.size();
      cerr << endl;
    }
    if (data.nspecies != xstr_pocc.species.size()) {
      message << MESSAGE << " [1] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species.size()[" << xstr_pocc.species.size() << "]" << endl << xstr_pocc;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != xstr_pocc.species_pp.size()) {
      message << MESSAGE << " [2] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species_pp.size()[" << xstr_pocc.species_pp.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != xstr_pocc.species_pp_type.size()) {
      message << MESSAGE << " [3] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species_pp_type.size()[" << xstr_pocc.species_pp_type.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != xstr_pocc.species_pp_version.size()) {
      message << MESSAGE << " [4] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species_pp_version.size()[" << xstr_pocc.species_pp_version.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != xstr_pocc.species_pp_ZVAL.size()) {
      message << MESSAGE << " [5] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species_pp_ZVAL.size()[" << xstr_pocc.species_pp_ZVAL.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != data.vspecies_pp_AUID.size()) {
      message << MESSAGE << " [5] - data.nspecies[" << data.nspecies << "]!=data.vspecies_pp_AUID.size()[" << data.vspecies_pp_AUID.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != xstr_pocc.species_volume.size()) {
      message << MESSAGE << " [6] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species_volume.size()[" << xstr_pocc.species_volume.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }
    if (data.nspecies != xstr_pocc.species_mass.size()) {
      message << MESSAGE << " [7] - data.nspecies[" << data.nspecies << "]!=xstr_pocc.species_mass.size()[" << xstr_pocc.species_mass.size() << "]";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
    }

    // parent structure
    if (pcalc.xstr_sym.atoms.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "pcalc.xstr_sym was not found", _RUNTIME_ERROR_);
    }
    const xstructure& xstr_pocc_parent = pcalc.xstr_sym;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xstr_pocc_parent=" << endl;
      cerr << xstr_pocc_parent << endl;
    }
    xstructure xstr_sp;
    xstructure xstr_sc;
    stringstream sss;
    aurostd::xoption vpflow_edata; // DX20180823 - added xoption
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " starting EDATA analysis" << endl;
    }
    sss << pflow::PrintData(xstr_pocc_parent, xstr_sp, xstr_sc, vpflow_edata, "EDATA"); // 1=EDATA //CO20171025 //CO20171027 //DX20210301 void to string output
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " EDATA analysis finished" << endl;
    }
    xmatrix<double> klattice;
    // do NOT write out file, not all the properties are applicable
    if (vpflow_edata.flag("EDATA::CALCULATED")) {
      data.Bravais_lattice_lattice_type_orig = vpflow_edata.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE"); // Bravais_lattice_lattice_type_orig
      data.Bravais_lattice_lattice_variation_type_orig = vpflow_edata.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE"); // Bravais_lattice_lattice_variation_type_orig
      data.Bravais_lattice_lattice_system_orig = vpflow_edata.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM"); // Bravais_lattice_lattice_system_orig
      data.Bravais_superlattice_lattice_type_orig = vpflow_edata.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE"); // Bravais_superlattice_lattice_type_orig
      data.Bravais_superlattice_lattice_variation_type_orig = vpflow_edata.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE"); // Bravais_superlattice_lattice_variation_type_orig
      data.Bravais_superlattice_lattice_system_orig = vpflow_edata.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM"); // Bravais_superlattice_lattice_system_orig
      data.Pearson_symbol_superlattice_orig = vpflow_edata.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE"); // Pearson_symbol_superlattice_orig
      // reciprocal_geometry_orig
      klattice = ReciprocalLattice(xstr_pocc.lattice, xstr_pocc.scale);
      data_abcabc = Getabc_angles(klattice, DEGREES);
      data.vreciprocal_geometry_orig = aurostd::xvector2vector(data_abcabc);
      data.reciprocal_geometry_orig = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vreciprocal_geometry_orig, _AFLOWLIB_DATA_GEOMETRY_PREC_), ",");
      data.reciprocal_volume_cell_orig = det(klattice);
      //
      data.reciprocal_lattice_type_orig = vpflow_edata.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE"); // reciprocal_lattice_type_orig
      data.reciprocal_lattice_variation_type_orig = vpflow_edata.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE"); // reciprocal_lattice_variation_type_orig
    }
    // print
    if (AFLOWLIB_VERBOSE && !data.Bravais_lattice_lattice_type_orig.empty()) {
      cout << MESSAGE << " Bravais_lattice_lattice_type_orig=" << data.Bravais_lattice_lattice_type_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.Bravais_lattice_lattice_variation_type_orig.empty()) {
      cout << MESSAGE << " Bravais_lattice_lattice_variation_type_orig=" << data.Bravais_lattice_lattice_variation_type_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.Bravais_lattice_lattice_system_orig.empty()) {
      cout << MESSAGE << " Bravais_lattice_lattice_system_orig=" << data.Bravais_lattice_lattice_system_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.Bravais_superlattice_lattice_type_orig.empty()) {
      cout << MESSAGE << " Bravais_superlattice_lattice_type_orig=" << data.Bravais_superlattice_lattice_type_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.Bravais_superlattice_lattice_variation_type_orig.empty()) {
      cout << MESSAGE << " Bravais_superlattice_lattice_variation_type_orig=" << data.Bravais_superlattice_lattice_variation_type_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.Bravais_superlattice_lattice_system_orig.empty()) {
      cout << MESSAGE << " Bravais_superlattice_lattice_system_orig=" << data.Bravais_superlattice_lattice_system_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.Pearson_symbol_superlattice_orig.empty()) {
      cout << MESSAGE << " Pearson_symbol_superlattice_orig=" << data.Pearson_symbol_superlattice_orig << endl;
    }
    //
    if (AFLOWLIB_VERBOSE && !data.reciprocal_geometry_orig.empty()) {
      cout << MESSAGE << " reciprocal_geometry_orig=" << data.reciprocal_geometry_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && data.reciprocal_volume_cell_orig != AUROSTD_NAN) {
      cout << MESSAGE << " reciprocal_volume_cell_orig=" << data.reciprocal_volume_cell_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.reciprocal_lattice_type_orig.empty()) {
      cout << MESSAGE << " reciprocal_lattice_type_orig=" << data.reciprocal_lattice_type_orig << endl;
    }
    if (AFLOWLIB_VERBOSE && !data.reciprocal_lattice_variation_type_orig.empty()) {
      cout << MESSAGE << " reciprocal_lattice_variation_type_orig=" << data.reciprocal_lattice_variation_type_orig << endl;
    }

    if (false) {
      // there is no relaxed pocc material: so copy over the original properties
      // this may change in the future
      // set
      data.Bravais_lattice_lattice_type = data.Bravais_lattice_lattice_type_orig;
      data.Bravais_lattice_lattice_variation_type = data.Bravais_lattice_lattice_variation_type_orig;
      data.Bravais_lattice_lattice_system = data.Bravais_lattice_lattice_system_orig;
      data.Bravais_superlattice_lattice_type = data.Bravais_superlattice_lattice_type_orig;
      data.Bravais_superlattice_lattice_variation_type = data.Bravais_superlattice_lattice_variation_type_orig;
      data.Bravais_superlattice_lattice_system = data.Bravais_superlattice_lattice_system_orig;
      data.Pearson_symbol_superlattice = data.Pearson_symbol_superlattice_orig;
      // reciprocal_geometry_relax //CO20220719 _relax
      data.vreciprocal_geometry_relax = data.vreciprocal_geometry_orig; // CO20220719 _relax
      data.reciprocal_geometry_relax = data.reciprocal_geometry_orig; // CO20220719 _relax
      data.reciprocal_volume_cell = data.reciprocal_volume_cell_orig;
      //
      data.reciprocal_lattice_type = data.reciprocal_lattice_type_orig;
      data.reciprocal_lattice_variation_type = data.reciprocal_lattice_variation_type_orig;
      // print
      if (AFLOWLIB_VERBOSE && !data.Bravais_lattice_lattice_type.empty()) {
        cout << MESSAGE << " Bravais_lattice_lattice_type=" << data.Bravais_lattice_lattice_type << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.Bravais_lattice_lattice_variation_type.empty()) {
        cout << MESSAGE << " Bravais_lattice_lattice_variation_type=" << data.Bravais_lattice_lattice_variation_type << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.Bravais_lattice_lattice_system.empty()) {
        cout << MESSAGE << " Bravais_lattice_lattice_system=" << data.Bravais_lattice_lattice_system << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.Bravais_superlattice_lattice_type.empty()) {
        cout << MESSAGE << " Bravais_superlattice_lattice_type=" << data.Bravais_superlattice_lattice_type << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.Bravais_superlattice_lattice_variation_type.empty()) {
        cout << MESSAGE << " Bravais_superlattice_lattice_variation_type=" << data.Bravais_superlattice_lattice_variation_type << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.Bravais_superlattice_lattice_system.empty()) {
        cout << MESSAGE << " Bravais_superlattice_lattice_system=" << data.Bravais_superlattice_lattice_system << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.Pearson_symbol_superlattice.empty()) {
        cout << MESSAGE << " Pearson_symbol_superlattice=" << data.Pearson_symbol_superlattice << endl;
      }
      //
      if (AFLOWLIB_VERBOSE && !data.reciprocal_geometry_orig.empty()) {
        cout << MESSAGE << " reciprocal_geometry_orig=" << data.reciprocal_geometry_orig << endl;
      }
      if (AFLOWLIB_VERBOSE && data.reciprocal_volume_cell_orig != AUROSTD_NAN) {
        cout << MESSAGE << " reciprocal_volume_cell_orig=" << data.reciprocal_volume_cell_orig << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.reciprocal_lattice_type.empty()) {
        cout << MESSAGE << " reciprocal_lattice_type=" << data.reciprocal_lattice_type << endl;
      }
      if (AFLOWLIB_VERBOSE && !data.reciprocal_lattice_variation_type.empty()) {
        cout << MESSAGE << " reciprocal_lattice_variation_type=" << data.reciprocal_lattice_variation_type << endl;
      }
    }

    // anrl properties
    xstructure xstr_anrl = xstr_pocc_parent; // make a copy so functions don't pollute structure
    anrl::structure2anrl(xstr_anrl, xstr_anrl.sym_eps, SG_SETTING_ANRL);
    data.aflow_prototype_label_orig = xstr_anrl.prototype;
    if (AFLOWLIB_VERBOSE && !data.aflow_prototype_label_orig.empty()) {
      cout << MESSAGE << " aflow_prototype_label_orig=" << data.aflow_prototype_label_orig << endl;
    }
    data.aflow_prototype_params_list_orig = aurostd::joinWDelimiter(xstr_anrl.prototype_parameter_list, ",");
    if (AFLOWLIB_VERBOSE && !data.aflow_prototype_params_list_orig.empty()) {
      cout << MESSAGE << " aflow_prototype_params_list_orig=" << data.aflow_prototype_params_list_orig << endl;
    }
    // build data.aflow_prototype_params_values_orig from scratch
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " prototype_parameter_values(PARENT)=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(xstr_anrl.prototype_parameter_values, _AFLOWLIB_DATA_DOUBLE_PREC_), ",") << endl;
    }
    vector<double> prototype_parameter_values = xstr_anrl.prototype_parameter_values;
    if (prototype_parameter_values.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "AFLOW parameter values builder failed: prototype_parameter_values.empty()", _RUNTIME_ERROR_);
    }
    prototype_parameter_values[0] *= (xstr_pocc.GetVolume() / xstr_pocc_parent.GetVolume()); // scale for the volume of the pocc structure, remember: anrl is based on std_conv, not the primitive
    if (prototype_parameter_values.size() != xstr_anrl.prototype_parameter_values.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "AFLOW parameter values builder failed: prototype_parameter_values.size()!=xstr_anrl.prototype_parameter_values.size()", _RUNTIME_ERROR_);
    }
    data.aflow_prototype_params_values_orig = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(prototype_parameter_values, _AFLOWLIB_DATA_DOUBLE_PREC_), ",");
    if (AFLOWLIB_VERBOSE && !data.aflow_prototype_params_values_orig.empty()) {
      cout << MESSAGE << " aflow_prototype_params_values_orig=" << data.aflow_prototype_params_values_orig << endl;
    }

    // get other _atom/_cell properties
    double data_natoms = 0.0; // needs to be double for pocc
    for (i = 0; i < xstr_pocc.comp_each_type.size(); i++) {
      data_natoms += xstr_pocc.comp_each_type[i];
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " natoms=" << data_natoms << endl;
    }
    //
    if (data.volume_cell_orig != AUROSTD_NAN) {
      data.volume_atom_orig = data.volume_cell_orig / data_natoms;
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " volume_atom_orig=" << data.volume_atom_orig << endl;
      }
    }
    //
    if (data.energy_atom != AUROSTD_NAN) {
      data.energy_cell = data.energy_atom * data_natoms;
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " energy_cell=" << data.energy_cell << endl;
      }
    }
    //
    if (data.enthalpy_atom != AUROSTD_NAN) {
      data.enthalpy_cell = data.enthalpy_atom * data_natoms;
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " enthalpy_cell=" << data.enthalpy_cell << endl;
      }
    }
    //
    if (data.volume_cell_orig != AUROSTD_NAN) {
      data.density_orig = 0.0;
      for (i = 0; i < xstr_pocc.comp_each_type.size(); i++) {
        data.density_orig += xstr_pocc.comp_each_type[i] * GetAtomMass(xstr_pocc.species[i]);
      }
      data.density_orig /= data.volume_cell_orig;
      data.density_orig *= 1000.0; // grams instead of kilos
      data.density_orig *= 1e8 * 1e8 * 1e8; // cm^3 instead of A^3
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " DENSITY_ORIG (grams/cm^3) = " << data.density_orig << endl;
      }
    }

    // enthalpy_formation_atom
    if (data.enthalpy_atom != AUROSTD_NAN && data.enthalpy_cell != AUROSTD_NAN) {
      // the enthalpy_mix could be corrected directly, but cce requires knowledge of each structure
      // so we grab ALL the enthalpies and structures for the enthalpy_formation (+cce variants) calculation
      // then average at the end
      xstructure xstr_derivative;
      stringstream xstr_ss;
      bool FORMATION_CALC = false;
      aflowlib::_aflowlib_entry data_derivative;
      vector<double> v_Hfs;
      vector<double> v_Hfs_cce_300K;
      vector<double> v_Hfs_cce_0K;
      vector<double> v_Ts; // all _atom
      for (i = 0; i < pcalc.m_ARUN_directories.size(); i++) {
        xstr_derivative.clear();
        // if you want to grab from the derivative structure, do as below
        // following algorithm as SC does above (see fileX_LIB)
        // these approaches should be consolidated into a single algorithm so the MOST relaxed structure is always extracted
        // see, e.g., KBIN::GetMostRelaxedStructure() in aflow_ivasp
        // we avoid for now: do not fix what is not broken (until necessary)
        if (xstr_derivative.atoms.empty() && aurostd::CompressFileExist(directory_LIB + "/" + pcalc.m_ARUN_directories[i] + "/CONTCAR.static", filename)) {
          aurostd::compressfile2stringstream(filename, xstr_ss);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " found CONTCAR.static:" << endl;
            cerr << xstr_ss.str() << endl;
          }
          xstr_ss >> xstr_derivative;
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " loaded CONTCAR.static" << endl;
          }
        }
        if (xstr_derivative.atoms.empty()) {
          for (j = RELAX_MAX; j > 0 && xstr_derivative.atoms.empty(); j--) { // no relax0
            if (xstr_derivative.atoms.empty() && aurostd::CompressFileExist(directory_LIB + "/" + pcalc.m_ARUN_directories[i] + "/CONTCAR.relax" + aurostd::utype2string(j), filename)) {
              aurostd::compressfile2stringstream(filename, xstr_ss);
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " found CONTCAR.relax" + aurostd::utype2string(j) + ":" << endl;
                cerr << xstr_ss.str() << endl;
              }
              xstr_ss >> xstr_derivative;
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " loaded CONTCAR.relax" + aurostd::utype2string(j) << endl;
              }
            }
          }
        }
        if (xstr_derivative.atoms.empty()) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xstr_derivative cannot be extracted", _FILE_CORRUPT_);
        }
        // load in xstr_pocc (from xOUT)
        xstr_derivative.species = xstr_pocc.species; // for LIB2RAW_Calculate_FormationEnthalpy
        xstr_derivative.species_pp = xstr_pocc.species_pp; // for LIB2RAW_Calculate_FormationEnthalpy
        xstr_derivative.species_pp_version = xstr_pocc.species_pp_version; // for LIB2RAW_Calculate_FormationEnthalpy
        xstr_derivative.species_pp_ZVAL = xstr_pocc.species_pp_ZVAL; // for LIB2RAW_Calculate_FormationEnthalpy
        xstr_derivative.species_pp_vLDAU = xstr_pocc.species_pp_vLDAU; // keep the same structure // for LIB2RAW_Calculate_FormationEnthalpy
        pflow::fixEmptyAtomNames(xstr_derivative, true); // force fix
        // copy data over
        data_derivative = data;
        data_derivative.enthalpy_atom = v_energies[i];
        data_derivative.enthalpy_cell = v_energies[i] * xstr_derivative.atoms.size();
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTHALPY total E0 (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_cell << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTHALPY per atom E0/N (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_atom << "   " << directory_LIB << endl;
        }
        //
        FORMATION_CALC = LIB2RAW_Calculate_FormationEnthalpy(data_derivative, xstr_derivative, MESSAGE);
        if (FORMATION_CALC == false) {
          break;
        }
        // print out data for derivative structure START
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTHALPY FORMATION total E0 (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_formation_cell << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTHALPY FORMATION per atom E0/N (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_formation_atom << "   " << directory_LIB << endl;
        }
        // CO20200624 START - CCE
        if (AFLOWLIB_VERBOSE && data_derivative.enthalpy_formation_cce_300K_cell != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE total E(300K) (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_formation_cce_300K_cell << endl;
        }
        if (AFLOWLIB_VERBOSE && data_derivative.enthalpy_formation_cce_300K_atom != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(300K)/N (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_formation_cce_300K_atom << "   " << directory_LIB << endl;
        }
        if (AFLOWLIB_VERBOSE && data_derivative.enthalpy_formation_cce_0K_cell != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE total E(0K) (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_formation_cce_0K_cell << endl;
        }
        if (AFLOWLIB_VERBOSE && data_derivative.enthalpy_formation_cce_0K_atom != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(0K)/N (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.enthalpy_formation_cce_0K_atom << "   " << directory_LIB << endl;
        }
        // CO20200624 END - CCE
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTROPIC_TEMPERATURE (eV) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.entropic_temperature * KBOLTZEV << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTROPIC_TEMPERATURE (K) [" << pcalc.m_ARUN_directories[i] << "] = " << data_derivative.entropic_temperature << "   " << directory_LIB << endl;
        }
        // print out data for derivative structure END
        v_Hfs.push_back(data_derivative.enthalpy_formation_atom);
        if (data_derivative.enthalpy_formation_cce_300K_atom != AUROSTD_NAN) {
          v_Hfs_cce_300K.push_back(data_derivative.enthalpy_formation_cce_300K_atom);
        }
        if (data_derivative.enthalpy_formation_cce_0K_atom != AUROSTD_NAN) {
          v_Hfs_cce_0K.push_back(data_derivative.enthalpy_formation_cce_0K_atom);
        }
        v_Ts.push_back(data_derivative.entropic_temperature);
      }
      if (FORMATION_CALC == true) { // they must ALL be true
        // we gather v_Hfs as vector and convert to xvector for the function
        // keep as is: we don't want to create an xvector of defined size
        // the size of the vector tells us whether the values were extracted correctly
        if (v_dgs.size() != v_Hfs.size()) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_dgs.size()!=v_Hfs.size()", _FILE_CORRUPT_);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " xv_dgs=" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(xv_dgs, 5), ",") << endl;
          cerr << __AFLOW_FUNC__ << " v_Hfs=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(v_Hfs, 5), ",") << endl;
        }
        data.enthalpy_formation_atom = aurostd::meanWeighted(aurostd::vector2xvector(v_Hfs), xv_dgs);
        data.enthalpy_formation_cell = data.enthalpy_formation_atom * data_natoms;
        if (v_dgs.size() == v_Hfs_cce_300K.size()) { // they all worked
          data.enthalpy_formation_cce_300K_atom = aurostd::meanWeighted(aurostd::vector2xvector(v_Hfs_cce_300K), xv_dgs);
          data.enthalpy_formation_cce_300K_cell = data.enthalpy_formation_cce_300K_atom * data_natoms;
        }
        if (v_dgs.size() == v_Hfs_cce_0K.size()) { // they all worked
          data.enthalpy_formation_cce_0K_atom = aurostd::meanWeighted(aurostd::vector2xvector(v_Hfs_cce_0K), xv_dgs);
          data.enthalpy_formation_cce_0K_cell = data.enthalpy_formation_cce_0K_atom * data_natoms;
        }
        if (v_dgs.size() != v_Ts.size()) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "v_dgs.size()!=v_Ts.size()", _FILE_CORRUPT_);
        }
        data.entropic_temperature = aurostd::meanWeighted(aurostd::vector2xvector(v_Ts), xv_dgs);
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTHALPY FORMATION total E0 (eV) = " << data.enthalpy_formation_cell << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTHALPY FORMATION per atom E0/N (eV) = " << data.enthalpy_formation_atom << "   " << directory_LIB << endl;
        }
        // CO20200624 START - CCE
        if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_300K_cell != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE total E(300K) (eV) = " << data.enthalpy_formation_cce_300K_cell << endl;
        }
        if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_300K_atom != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(300K)/N (eV) = " << data.enthalpy_formation_cce_300K_atom << "   " << directory_LIB << endl;
        }
        if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_0K_cell != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE total E(0K) (eV) = " << data.enthalpy_formation_cce_0K_cell << endl;
        }
        if (AFLOWLIB_VERBOSE && data.enthalpy_formation_cce_0K_atom != AUROSTD_NAN) {
          cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(0K)/N (eV) = " << data.enthalpy_formation_cce_0K_atom << "   " << directory_LIB << endl;
        }
        // CO20200624 END - CCE
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTROPIC_TEMPERATURE (eV) = " << data.entropic_temperature * KBOLTZEV << endl;
        }
        if (AFLOWLIB_VERBOSE) {
          cout << MESSAGE << " ENTROPIC_TEMPERATURE (K) = " << data.entropic_temperature << "   " << directory_LIB << endl;
        }
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_LOCK
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for LOCK data
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_LOCK(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    // Stefano Curtarolo 2009-2010-2011-2012
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    data.vloop.emplace_back("lock");

    vector<string> vlock;
    vector<string> vtokens;
    string tmp;
    string::size_type loc;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB, _AFLOWLOCK_, directory_RAW, _AFLOWLOCK_, vfiles, MESSAGE); // OUTCAR.static
    aurostd::file2vectorstring(directory_RAW + "/" + _AFLOWLOCK_, vlock);
    _XHOST aus_XHOST;
    // ---------------------------------------------------------------
    data.vaflowlib_date.clear(); // clear here,
    for (size_t iline = 0; iline < vlock.size(); iline++) { // CO20200624 - adding lock date to vaflowlib_date  //grab both FIRST and LAST dates - data.vaflowlib_date.empty()
      if (aurostd::substring2bool(vlock[iline], "date=") && aurostd::substring2bool(vlock[iline], "[") && aurostd::substring2bool(vlock[iline], "]")) {
        loc = vlock[iline].find("date=");
        tmp = vlock[iline].substr(loc, string::npos);
        loc = tmp.find(']');
        tmp = tmp.substr(0, loc);
        tmp = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tmp);
        aurostd::StringSubstInPlace(tmp, "[", "");
        aurostd::StringSubstInPlace(tmp, "]", ""); // just in case
        aurostd::StringSubstInPlace(tmp, "date=", ""); // just in case
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " FOUND LOCK date = " << tmp << endl;
        }
        tmp = aflow_convert_time_ctime2aurostd(tmp); // CO20210626 - utc_offset already included inside
        if (!tmp.empty()) {
          if (data.vaflowlib_date.empty()) {
            data.vaflowlib_date.push_back(tmp);
          } // get first date
          else { // get last date
            if (data.vaflowlib_date.size() > 1) {
              data.vaflowlib_date.pop_back();
            }
            data.vaflowlib_date.push_back(tmp);
          }
        }
      }
    }
    if (data.vaflowlib_date.empty()) { // CO20210213
      // try looking for old style, should be the last item in lines containing -
      // look that very last item in line is between 2000 and current year
      double dtmp = 0.0;
      for (size_t iline = 0; iline < vlock.size(); iline++) { // CO20200624 - adding lock date to vaflowlib_date  //grab both FIRST and LAST dates - data.vaflowlib_date.empty()
        aurostd::string2tokens(vlock[iline], vtokens, " ");
        if (!vtokens.empty() && aurostd::isfloat(vtokens.back())) {
          dtmp = aurostd::string2utype<double>(vtokens.back());
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " dtmp=" << dtmp << endl;
          }
          if (aurostd::isinteger(dtmp) && (dtmp > 2000 && dtmp < (double) aurostd::get_year())) {
            aurostd::string2tokens(vlock[iline], vtokens, "-");
            if (!vtokens.empty()) {
              // START taken from above
              tmp = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vtokens.back());
              aurostd::StringSubstInPlace(tmp, "[", "");
              aurostd::StringSubstInPlace(tmp, "]", ""); // just in case
              aurostd::StringSubstInPlace(tmp, "date=", ""); // just in case
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " FOUND LOCK date = " << tmp << endl;
              }
              tmp = aflow_convert_time_ctime2aurostd(tmp); // CO20210626 - utc_offset already included inside
              if (!tmp.empty()) {
                if (data.vaflowlib_date.empty()) {
                  data.vaflowlib_date.push_back(tmp);
                } // get first date
                else { // get last date
                  if (data.vaflowlib_date.size() > 1) {
                    data.vaflowlib_date.pop_back();
                  }
                  data.vaflowlib_date.push_back(tmp);
                }
              }
              // STOP taken from above
            }
          }
        }
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " lock_date (START) = " << ((!data.vaflowlib_date.empty()) ? data.vaflowlib_date[0] : "unavailable") << endl;
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " lock_date (END) = " << ((data.vaflowlib_date.size() > 1) ? data.vaflowlib_date[1] : "unavailable") << endl;
    }
    // ---------------------------------------------------------------
    for (size_t iline = 0; iline < vlock.size() && data.aflow_version.empty(); iline++) {
      if (aurostd::substring2bool(vlock[iline], "NFS") && aurostd::substring2bool(vlock[iline], "(") && aurostd::substring2bool(vlock[iline], ")")) {
        aurostd::string2tokens(vlock[iline], vtokens);
        data.aflow_version = aurostd::RemoveWhiteSpaces(vtokens.at(vtokens.size() - 1));
        aurostd::StringSubstInPlace(data.aflow_version, "(", "");
        aurostd::StringSubstInPlace(data.aflow_version, ")", "");
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " aflow_version = " << ((!data.aflow_version.empty()) ? data.aflow_version : "unavailable") << endl;
    }
    // XHOST.CPU_Model ---------------------------------------------------------------
    for (size_t iline = 0; iline < vlock.size() && aus_XHOST.CPU_Model.empty(); iline++) {
      if (aurostd::substring2bool(vlock[iline], "XHOST.CPU_Model") && aurostd::substring2bool(vlock[iline], ":")) {
        aurostd::string2tokens(vlock[iline], vtokens, ":");
        aus_XHOST.CPU_Model = aurostd::RemoveWhiteSpaces(vtokens.at(vtokens.size() - 1));
        data.node_CPU_Model = aus_XHOST.CPU_Model;
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " XHOST.CPU_Model = " << ((!aus_XHOST.CPU_Model.empty()) ? aus_XHOST.CPU_Model : "unavailable") << endl;
    }
    // XHOST.CPU_Cores ---------------------------------------------------------------
    for (size_t iline = 0; iline < vlock.size() && aus_XHOST.CPU_Cores == 0; iline++) {
      if (aurostd::substring2bool(vlock[iline], "XHOST.CPU_Cores") && aurostd::substring2bool(vlock[iline], ":")) {
        aurostd::string2tokens(vlock[iline], vtokens, ":");
        aus_XHOST.CPU_Cores = aurostd::string2utype<int>(aurostd::RemoveWhiteSpaces(vtokens.at(vtokens.size() - 1)));
        aus_XHOST.CPU_Cores = aurostd::string2utype<int>(aurostd::RemoveWhiteSpaces(vtokens.at(vtokens.size() - 1)));
        data.node_CPU_Cores = aus_XHOST.CPU_Cores;
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " XHOST.CPU_Cores = " << ((aus_XHOST.CPU_Cores) ? aurostd::utype2string<int>(aus_XHOST.CPU_Cores) : "unavailable") << endl;
    }
    // XHOST.CPU_MHz ---------------------------------------------------------------
    for (size_t iline = 0; iline < vlock.size() && aus_XHOST.CPU_MHz.empty(); iline++) {
      if (aurostd::substring2bool(vlock[iline], "XHOST.CPU_MHz") && aurostd::substring2bool(vlock[iline], ":")) {
        aurostd::string2tokens(vlock[iline], vtokens, ":");
        aus_XHOST.CPU_MHz = aurostd::RemoveWhiteSpaces(vtokens.at(vtokens.size() - 1));
        data.node_CPU_MHz = ceil(aurostd::string2utype<double>(aus_XHOST.CPU_MHz));
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " XHOST.CPU_MHz = " << ((!aus_XHOST.CPU_MHz.empty()) ? aus_XHOST.CPU_MHz : "unavailable") << endl;
    }
    // XHOST.RAM_GB ---------------------------------------------------------------
    for (size_t iline = 0; iline < vlock.size() && aus_XHOST.RAM_GB < 0.001; iline++) {
      if (aurostd::substring2bool(vlock[iline], "XHOST.RAM_GB") && aurostd::substring2bool(vlock[iline], ":")) {
        aurostd::string2tokens(vlock[iline], vtokens, ":");
        aus_XHOST.RAM_GB = ceil(aurostd::string2utype<double>(aurostd::RemoveWhiteSpaces(vtokens.at(vtokens.size() - 1))));
        data.node_RAM_GB = aus_XHOST.RAM_GB;
      }
    }
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " XHOST.RAM_GB = " << ((aus_XHOST.RAM_GB) ? aurostd::utype2string<double>(aus_XHOST.RAM_GB) : "unavailable") << endl;
    }
    // REMOVING  ---------------------------------------------------------------
    const vector<string> vremove = {"aflow.in~",     "WAVECAR",        "REPORT",         "POTCAR.relax1",  "POTCAR.relax2",  "POTCAR.relax3",  "POTCAR.static",
                                    "POTCAR.bands",  "AECCAR0",        "AECCAR1",        "AECCAR2",        "AECCAR1.static", "AECCAR0.bands",  "AECCAR1.bands",
                                    "AECCAR2.bands", "AECCAR0.relax1", "AECCAR1.relax1", "AECCAR2.relax1", "AECCAR0.relax2", "AECCAR1.relax2", "AECCAR2.relax2"};
    for (const auto& dir_entry : std::filesystem::recursive_directory_iterator(directory_LIB)) {
      // iterate over directory, then try to match files
      for (const auto& remove : vremove) {
        for (const auto& ext : XHOST.vext) {
          if (dir_entry.path().filename().string() == remove + ext && dir_entry.exists()) {
            std::filesystem::remove(dir_entry.path());
            if (AFLOWLIB_VERBOSE) {
              cout << MESSAGE << " removing file: " << remove + ext << endl;
            }
            goto dir_entry_iterate1; // double break
          }
        } // ext
      } // remove
    dir_entry_iterate1:;
    }

    vector<string> vremoveERROR = {"aflow.error.rotmat",      "aflow.error.nbands",       "aflow.error.symprec", "aflow.error.read_kpoints_rd_sym",
                                   "aflow.error.ibzkpt",      "aflow.error.gamma_shift",  "aflow.error.mpich11", "aflow.error.mpich139",
                                   "aflow.error.nkxyz_ikptd", "aflow.error.eddrmm",       "aflow.error.lreal",   "aflow.error.exccor",
                                   "aflow.error.brmix",       "aflow.error.dav",          "aflow.error.edddav",  "aflow.error.efield_pead",
                                   "aflow.error.zpotrf",      "aflow.error.zpotrf_potim", "aflow.error.natoms",  "aflow.error.psmaxn",
                                   "aflow.error.npar",        "aflow.error.nparc",        "aflow.error.nparn",   "aflow.error.npar_remove",
                                   "aflow.error.csloshing",   "aflow.error.dentet"};
    for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
      for (size_t iremove = 0; iremove < vremoveERROR.size(); iremove++) {
        if (aurostd::FileExist(directory_LIB + "/" + vremoveERROR[iremove] + XHOST.vext[iext])) {
          aurostd::RemoveFile(directory_LIB + "/" + vremoveERROR[iremove] + XHOST.vext[iext]);
        }
      } // iremove
    } // iext

    // done   ---------------------------------------------------------------
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_PATCH
// ***************************************************************************
namespace aflowlib {
  /// @brief LIB2RAW subroutine for patching
  /// @param directory_LIB directory to use for LIB
  /// @param directory_RAW directory to use for RAW
  /// @param vfiles list to keep track of added files
  /// @param data data structure for the library entry
  /// @param MESSAGE message to use for logging errors
  /// @authors
  /// @mod{ST,20241022,created doxy\, change signature\, optimize\, cleanup}
  void LIB2RAW_Loop_PATCH(const string& directory_LIB, const string& directory_RAW, vector<string>& vfiles, aflowlib::_aflowlib_entry& data, const string& MESSAGE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " [1]" << endl;
    }
    // Stefano Curtarolo 2009-2010-2011-2012
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_LIB << endl;
    }
    if (AFLOWLIB_VERBOSE && false) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << directory_RAW << endl; // to avoid unused
    }
    if (AFLOWLIB_VERBOSE && false) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - begin " << vfiles.size() << endl; // to avoid unused
    }

    // --------------------------------------------------------------------
    // CHECK existence of DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT
    string file;
    for (size_t iext = 1; iext < XHOST.vext.size() && file.empty(); iext++) { // SKIM uncompressed
      if (aurostd::FileExist(directory_LIB + "/" + DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT + XHOST.vext[iext])) {
        file = directory_LIB + "/" + DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT + XHOST.vext[iext];
      }
    }

    if (!file.empty()) {
      // if(AFLOWLIB_VERBOSE)
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " FOUND " << file << endl;
      }
    } else {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " WRITING " << directory_LIB << "/" << DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << endl;
      }
      KBIN::VASP_Write_ppAUID_FILE(directory_LIB, data.vspecies_pp_AUID, data.vspecies);
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " COMPRESSING " << directory_LIB << "/" << DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << endl;
      }
      aurostd::CompressFile(directory_LIB + "/" + DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT);
    }
    // --------------------------------------------------------------------
    // CHECK existence of VASP_POTCAR_AUID in aflow.in
    vector<string> vVASP_POTCAR_AUID;
    aurostd::file2vectorstring(directory_LIB + "/" + _AFLOWIN_, vVASP_POTCAR_AUID);
    auto predicate = [](const string& line) { return !aurostd::substring2bool(line, "VASP_POTCAR_AUID"); };
    vVASP_POTCAR_AUID.erase(std::remove_if(vVASP_POTCAR_AUID.begin(), vVASP_POTCAR_AUID.end(), predicate), vVASP_POTCAR_AUID.end());
    // todo: need to make a function for this grep style filtering
    if (LDEBUG) {
      cout << __AFLOW_FUNC__ << " vVASP_POTCAR_AUID.size()=" << vVASP_POTCAR_AUID.size() << endl;
    }

    if (!vVASP_POTCAR_AUID.empty()) {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " FOUND [VASP_POTCAR_AUID] in " << _AFLOWIN_ << ": " << vVASP_POTCAR_AUID[0] << endl;
      }
    } else {
      if (AFLOWLIB_VERBOSE) {
        cout << MESSAGE << " FIXING " << directory_LIB << "/" << _AFLOWIN_ << endl;
      }
      KBIN::VASP_Write_ppAUID_AFLOWIN(directory_LIB, data.vspecies_pp_AUID, data.vspecies);
    }
    // done   ---------------------------------------------------------------
    if (AFLOWLIB_VERBOSE) {
      cout << MESSAGE << " " << __AFLOW_FUNC__ << " - end " << directory_LIB << endl;
    }
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::XPLUG
// ***************************************************************************
namespace aflowlib {
  /// @brief zip directories of completed AFLOW runs
  /// @param xplug_opts command line options
  /// @authors
  /// @mod{CO,20200501,created function}
  /// @mod{SD,20240416,rewritten using filesystem and libarchive}
  /// @mod{HE,20240907,changing to relative path zipping}
  void XPLUG(const aurostd::xoption& xplug_opts) {
    const string directory = aurostd::CleanFileName(xplug_opts.getattachedscheme("XPLUG::DIRECTORY"));
    if (directory.empty()) {
      const string message = "XPLUG::DIRECTORY cannot be empty";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
    }
    if (!aurostd::FileExist(directory)) {
      const string message = "XPLUG::DIRECTORY does not exist";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    const bool clean = xplug_opts.flag("XPLUG::CLEAN");
    double last_mod_days = 3.0;
    if (!xplug_opts.getattachedscheme("XPLUG::LAST_MOD_DAYS").empty()) {
      last_mod_days = aurostd::string2utype<double>(xplug_opts.getattachedscheme("XPLUG::LAST_MOD_DAYS"));
    }
    int zip_size_gb = 9;
    if (!xplug_opts.getattachedscheme("XPLUG::ZIP_SIZE_GB").empty()) {
      zip_size_gb = aurostd::string2utype<int>(xplug_opts.getattachedscheme("XPLUG::ZIP_SIZE_GB"));
    }
    std::string prefix;
    if (!xplug_opts.getattachedscheme("XPLUG::PREFIX").empty()) {
      prefix = xplug_opts.getattachedscheme("XPLUG::PREFIX");
    }
    cerr << "prefix: " << prefix << endl;
    std::filesystem::path relative_path;
    if (xplug_opts.getattachedscheme("XPLUG::RELATIVE").empty()) {
      relative_path = fs::current_path();
    } else {
      const std::string relative_path_raw = xplug_opts.getattachedscheme("XPLUG::RELATIVE");
      relative_path = (aurostd::CleanFileName(relative_path_raw));
    }
    const std::filesystem::path root_path(directory);
    vector<std::filesystem::path> vdir_all;
    vector<std::filesystem::path> vdir_zip;
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::recursive_directory_iterator(root_path)) {
      if (dir_entry.is_directory()) {
        vdir_all.emplace_back(dir_entry.path());
      }
    }
    if (vdir_all.empty()) {
      vdir_all = {root_path};
    }
    XPLUG_Check(vdir_all, vdir_zip, clean, last_mod_days);
    XPLUG_Zip(vdir_zip, zip_size_gb, relative_path, prefix);
  }

  /// @brief finds valid directories for zipping
  /// @param vdir_all directories to check
  /// @param vdir_zip directories to zip
  /// @param clean invalid directories
  /// @param last_mod_days time, in days, after which to clean incomplete runs
  /// @authors
  /// @mod{SD,20240416,created function}
  void XPLUG_Check(const vector<std::filesystem::path>& vdir_all, vector<std::filesystem::path>& vdir_zip, const bool clean, double last_mod_days) {
    const std::set<string> required_singular = {_AFLOWLOCK_,          "aflow.agroup.out.xz", "aflow.fgroup.out.xz", "aflow.iatoms.out.xz", "aflow.in", "aflow.pgroup.out.xz", "aflow.pseudopotential_auid.out.xz",
                                                "aflow.qmvasp.out.xz"};
    const vector<std::regex> required_plural = {std::regex("CHGCAR.*.xz"),   std::regex("CHG.*.xz"),         std::regex("CONTCAR.*.xz"), std::regex("DOSCAR.*.xz"), std::regex("EIGENVAL.*.xz"),
                                                std::regex("IBZKPT.*.xz"),   std::regex("INCAR.*.xz"),       std::regex("KPOINTS.*.xz"), std::regex("PCDAT.*.xz"),  std::regex("POSCAR.*.xz"),
                                                std::regex("vasp.out.*.xz"), std::regex("vasprun.xml.*.xz"), std::regex("XDATCAR.*.xz")};
    const std::regex regex_oszicar("OSZICAR.*.xz");
    const std::regex regex_outcar("OUTCAR.*.xz");

    vector<std::filesystem::path> vfile;
    vector<string> voutcar;
    vector<string> voszicar;
    std::set<string> sfile;
    bool valid = true;
    bool found = false;
    for (const std::filesystem::path& dir : vdir_all) {
      valid = true;
      vfile.clear();
      voutcar.clear();
      voszicar.clear();
      sfile.clear();

      // Collect the files as a set (for singluar) and vector (for plural)
      for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(dir)) {
        if (dir_entry.is_regular_file()) {
          sfile.insert(dir_entry.path().filename().string());
          vfile.emplace_back(dir_entry.path());
        }
      }

      // Check if LOCK file exists, otherwise skip directory
      if (!sfile.count(_AFLOWLOCK_)) {
        continue;
      }

      // Check if all required singular files exist, otherwise directory is invalid
      for (const string& required_file : required_singular) {
        if (!sfile.count(required_file)) {
          valid = false;
          break;
        }
      }

      // If directory is invalid, we can clean it
      if (!valid) {
        if (!clean) {
          continue;
        }
        const std::set<string>::iterator it_vasp = sfile.find("vasp.out");
        if (it_vasp != sfile.end()) {
          const double last_mod_vasp = (long double) aurostd::SecondsSinceFileModified(vfile[std::distance(sfile.begin(), it_vasp)]) / 3600;
          if (last_mod_vasp > last_mod_days) {
            KBIN::Clean(dir.string());
          }
        }
        continue;
      }

      // Check if all required plural files exist, otherwise directory is invalid
      for (const std::regex& required_file : required_plural) {
        found = false;
        for (const std::filesystem::path& file : vfile) {
          if (std::regex_match(file.filename().string(), required_file)) {
            found = true;
            break;
          }
        }
        if (!found) {
          valid = false;
          break;
        }
      }

      // If directory is invalid, we can clean it
      if (!valid) {
        if (!clean) {
          continue;
        }
        KBIN::Clean(dir.string());
        continue;
      }

      // Check that the run has been converged, otherwise directory is invalid
      for (const std::filesystem::path& file : vfile) {
        if (std::regex_match(file.filename().string(), regex_oszicar)) {
          voszicar.emplace_back(file.string());
        } else if (std::regex_match(file.filename().string(), regex_outcar)) {
          voutcar.emplace_back(file.string());
        }
      }
      if (!voszicar.empty() && voszicar.size() == voutcar.size()) {
        for (size_t icar = 0; icar < voszicar.size() && valid; icar++) {
          if (KBIN::VASP_OSZICARUnconverged(voszicar[icar], voutcar[icar]) && KBIN::VASP_getNSTEPS(voszicar[icar]) < AFLOWRC_MAX_VASP_NELM) {
            valid = false;
          }
        }
      } else {
        valid = false;
      }

      // If directory is invalid, we can clean it
      if (!valid) {
        if (!clean) {
          continue;
        }
        KBIN::Clean(dir.string());
        continue;
      }

      // Directory is valid
      vdir_zip.push_back(dir);
    }
  }

  /// @brief zip valid directories
  /// @param vdir_zip directories to zip
  /// @param zip_size_gb maximum size, in GB, of each zip file
  /// @param relative_path store paths in zip relative to this one
  /// @param prefix to add to the archive names
  /// @authors
  /// @mod{SD,20240416,created function}
  /// @mod{HE,20240906,changing to relative path zipping}
  void XPLUG_Zip(const vector<std::filesystem::path>& vdir_zip, const int zip_size_gb, const std::filesystem::path& relative_path, const std::string& prefix) {
    if (vdir_zip.empty()) {
      return;
    }
    string archive_name_base = "update_" + prefix;
    if (!prefix.empty()) {
      archive_name_base += "_";
    }
    archive_name_base = archive_name_base + aurostd::get_datetime_formatted("", true, "-", "") + "_" + XHOST.hostname.substr(0, XHOST.hostname.find('.'));
    string archive_name;
    const unsigned long long int zip_size = 1E9 * zip_size_gb;
    unsigned long long int file_size = 0;
    const vector<std::filesystem::path> vdir;
    vector<std::filesystem::path> vzip_tmp;
    const vector<long long int> vsize;
    vector<vector<std::filesystem::path>> vzip;
    vector<string> vfile;
    vector<fs::path> vname;

    // Lambda expression to ignore certain files
    const std::function<bool(string)> ignore_file = [](const string& name) {
      if (name == "core" || std::regex_match(name, std::regex("POTCAR.*")) || std::regex_match(name, std::regex("AECCAR.*"))) {
        return true;
      }
      return false;
    };

    // Collect directories based on file size
    for (const std::filesystem::path& dir : vdir_zip) {
      for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(dir)) {
        if (dir_entry.is_regular_file() && !ignore_file(dir_entry.path().filename().string())) {
          file_size += dir_entry.file_size();
        }
      }
      if (file_size < zip_size) {
        vzip_tmp.push_back(dir);
      } else {
        vzip.push_back(vzip_tmp);
        vzip_tmp.clear();
        file_size = 0;
      }
    }
    vzip.push_back(vzip_tmp);

    // Compress files
    for (size_t iz = 0; iz < vzip.size(); iz++) {
      archive_name = archive_name_base + "_" + aurostd::utype2string(iz + 1) + "_of_" + aurostd::utype2string(vzip.size());
      vfile.clear();
      for (const std::filesystem::path& dir : vzip[iz]) {
        for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(dir)) {
          if (dir_entry.is_regular_file() && !ignore_file(dir_entry.path().filename().string())) {
            vfile.emplace_back(dir_entry.path().string());
          }
        }
      }
      vname.push_back(aurostd::CompressFiles(vfile, relative_path, archive_name, aurostd::compression_type::ZIP, true));
    }

    // Append hash to archive names
    for (const fs::path& original_path : vname) {
      // HE20240908 SHA1 is a fast hash on modern hardware; use sha1sum to check file integrity
      const std::string hash_type = "SHA1";
      const string hash = aurostd::file2hash(original_path, hash_type);
      const std::string new_name = original_path.stem().string() + "_" + hash_type + "_" + hash + original_path.extension().string();
      fs::rename(original_path, original_path.parent_path() / new_name);
    }

    // Delete the zipped directories
    for (const std::filesystem::path& dir : vdir_zip) {
      std::filesystem::remove_all(dir);
    }
  }
} // namespace aflowlib

// ***************************************************************************
// void aflowlib::VaspFileExist(const string& str_dir, const string& FILE)
// ***************************************************************************
namespace aflowlib {
  /// @brief check if vasp files exist
  /// @param str_dir directory to check
  /// @param FILE base file to check
  /// @return true any available vasp files exist
  /// @authors
  /// @mod{ST,20241022,created doxy}
  bool VaspFileExist(const string& str_dir, const string& FILE) {
    bool RUN_FLAG = false;
    if (!RUN_FLAG && aurostd::FileExist(str_dir + "/" + FILE)) {
      RUN_FLAG = true;
    }
    if (!RUN_FLAG && aurostd::CompressFileExist(str_dir + "/" + FILE)) {
      RUN_FLAG = true;
    }
    if (!RUN_FLAG && aurostd::FileExist(str_dir + "/" + FILE + ".static")) {
      RUN_FLAG = true;
    }
    if (!RUN_FLAG && aurostd::CompressFileExist(str_dir + "/" + FILE + ".static")) {
      RUN_FLAG = true;
    }
    if (!RUN_FLAG && aurostd::FileExist(str_dir + "/" + FILE + ".bands")) {
      RUN_FLAG = true;
    }
    if (!RUN_FLAG && aurostd::CompressFileExist(str_dir + "/" + FILE + ".bands")) {
      RUN_FLAG = true;
    }
    if (!RUN_FLAG) {
      cerr << FILE + " or " + FILE + ".static/bands or " + FILE + ".static/bands.EXT not found in the directory!" << endl;
    }
    return RUN_FLAG;
  }
} // namespace aflowlib

// ***************************************************************************
// void aflowlib::vaspfile2stringstream
// ***************************************************************************
namespace aflowlib {
  string vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& sss) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    // if XXXX.bands exists, you may also use this function by setting FILE=vaspfile.bands
    sss.clear();
    sss.str(std::string(""));
    bool gfound = false;
    if (!gfound && (FILE == "EIGENVAL" || FILE == "KPOINTS")) {
      gfound = true;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FILE=" << FILE << endl;
      }
      deque<string> vtype = {"", ".bz2", ".gz", ".xz", ".bands", ".bands.bz2", ".bands.gz", ".bands.xz"}; // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for (size_t i = 0; i < vtype.size(); i++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " TESTING FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
        }
        if (aurostd::FileExist(str_dir + "/" + FILE + vtype[i])) {
          aurostd::compressfile2stringstream(str_dir + "/" + FILE + vtype[i], sss);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " FOUND FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
          }
          return aurostd::CleanFileName(str_dir + "/" + FILE + vtype[i]);
        }
      }
      // ME20190627 BEGIN
      const string message = FILE + " or " + FILE + ".bands or " + FILE + ".bands.EXT not found in the directory, aborting!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      // ME20190627 END
    }
    if (!gfound && (FILE == "DOSCAR")) {
      gfound = true;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FILE=" << FILE << endl;
      }
      deque<string> vtype = {"", ".bz2", ".gz", ".xz", ".static", ".static.bz2", ".static.gz", ".static.xz"}; // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for (size_t i = 0; i < vtype.size(); i++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " TESTING FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
        }
        if (aurostd::FileExist(str_dir + "/" + FILE + vtype[i])) {
          aurostd::compressfile2stringstream(str_dir + "/" + FILE + vtype[i], sss);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " FOUND FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
          }
          return aurostd::CleanFileName(str_dir + "/" + FILE + vtype[i]);
        }
      }
      // ME20190627 BEGIN
      const string message = FILE + " or " + FILE + ".static or " + FILE + ".static.EXT not found in the directory, aborting!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      // ME20190627 END
    }
    if (!gfound && (FILE == "POSCAR" || FILE == "CONTCAR")) { // RF20200409
      gfound = true;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FILE=" << FILE << endl;
      }
      deque<string> vtype = {"",           ".bz2",    ".gz",         ".xz",        ".bands",     ".bands.bz2", ".bands.gz",   ".bands.xz",  ".static",   ".static.bz2", ".static.gz",
                             ".static.xz", ".relax2", ".relax2.bz2", ".relax2.gz", ".relax2.xz", ".relax1",    ".relax1.bz2", ".relax1.gz", ".relax1.xz"}; // have to add emptyness to vtype at the beginning //RF20200409
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for (size_t i = 0; i < vtype.size(); i++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " TESTING FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
        }
        if (aurostd::FileExist(str_dir + "/" + FILE + vtype[i])) {
          aurostd::compressfile2stringstream(str_dir + "/" + FILE + vtype[i], sss);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " FOUND FILE=" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
          }
          return aurostd::CleanFileName(str_dir + "/" + FILE + vtype[i]);
        }
      }
      // ME20190627 BEGIN
      const string message = FILE + " or " + FILE + ".bands/static/relax or " + FILE + ".bands./static/relax.EXT not found in the directory, aborting!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      // ME20190627 END
    }
    if (!gfound) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FILE=" << FILE << endl;
      }
      deque<string> vtype = {"",        ".bz2",        ".gz",        ".xz",        ".static", ".static.bz2", ".static.gz", ".static.xz",
                             ".relax1", ".relax1.bz2", ".relax1.gz", ".relax1.xz", ".bands",  ".bands.bz2",  ".bands.gz",  ".bands.xz"}; // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for (size_t i = 0; i < vtype.size(); i++) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " TESTING FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
        }
        if (aurostd::FileExist(str_dir + "/" + FILE + vtype[i])) {
          aurostd::compressfile2stringstream(str_dir + "/" + FILE + vtype[i], sss);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " FOUND FILE=[" << str_dir + "/" + FILE + vtype[i] << "]" << endl;
          }
          return aurostd::CleanFileName(str_dir + "/" + FILE + vtype[i]);
        }
      }
      // ME20190627 BEGIN
      const string message = FILE + " or " + FILE + ".static/relax1/bands or " + FILE + ".static/relax1/bands.EXT not found in the directory, aborting!";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      // ME20190627 END
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return string("");
  }
} // namespace aflowlib

namespace aflowlib {
  string vaspfile2stringstream(const string& str_dir, const string& FILE) {
    stringstream sss;
    return vaspfile2stringstream(str_dir, FILE, sss);
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::LIB2LIB
// ***************************************************************************
// Added by CT20181212
namespace aflowlib {
  /// @brief LIB2RAW subroutine for archiving runs
  /// @param options options in the form of dir or all[,dir]
  /// @param flag_FORCE whether to force if target directory already exists
  /// @param LOCAL perform database operations locally only, useful for testing
  /// @return true on succesful completion, otherwise false
  /// @authors
  /// @mod{ST,20241022,created doxy\, optimize\, cleanup}
  bool LIB2LIB(const string& options, bool flag_FORCE, bool LOCAL) {
    const bool LDEBUG = (false || XHOST.DEBUG); // CO20200624
    stringstream message;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " options=" << options << endl;
    }
    vector<string> tokens;

    if (!aurostd::substring2bool(options, "POCC")) {
      aurostd::string2tokens(options, tokens, ","); // QUICK FIX
    } else {
      tokens.clear();
      tokens.push_back(options);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " tokens.size()=" << tokens.size() << endl;
    }
    if (tokens.empty() || tokens.size() > 2) {
      init::ErrorOption(options, __AFLOW_FUNC__, aurostd::liststring2string("aflow --lib2lib=directory", "aflow --lib2lib=all[,dir]"));
    }

    if (!tokens.empty()) {
      if (tokens[0] == "all") {
        XHOST.sensors_allowed = false;
        XHOST.vflag_pflow.flag("MULTI=SH", false);
        return LIB2RAW_ALL(options, flag_FORCE);
        XHOST.sensors_allowed = true;
      }
    }

    string directory_LIB;
    string PROJECT_LIBRARY;
    if (LOCAL) {
      flag_FORCE = true;
      string directory = options;
      CleanDirectoryLIB(directory);
      directory_LIB = directory;
      PROJECT_LIBRARY = directory_LIB;
    } else { // normal run
      CheckMaterialServer(__AFLOW_FUNC__); // must be in AFLOW_MATERIALS_SERVER
      string directory = options;
      CleanDirectoryLIB(directory);

      PROJECT_LIBRARY = aflowlib::LIB2RAW_CheckProjectFromDirectory(directory);
      cout << __AFLOW_FUNC__ << " AFLOW V" << string(AFLOW_VERSION) << endl;
      cout << __AFLOW_FUNC__ << " directory=" << directory << endl;
      cout << __AFLOW_FUNC__ << " hostname=" << XHOST.hostname << endl;
      directory_LIB = directory;
      // wait for temperature (to avoid overheat)
      // WAIT FOR TEMPERATURE
      vector<string> vmessage;
      vmessage.push_back(XHOST.sPID);
      vmessage.emplace_back("       dir=" + directory_LIB);
      init::WaitTEMP(AFLOW_CORE_TEMPERATURE_LIB2RAW, cout, true, vmessage); // fixed at define at the beginning
    }

    bool run_directory = false;
    bool agl_aflowin_found = false;
    bool ael_aflowin_found = false;
    bool apl_aflowin_found = false; // ME20210927
    bool qha_aflowin_found = false; // AS20200901
    string AflowInName = _AFLOWIN_;
    string FileLockName = _AFLOWLOCK_;
    string stmp;
    vector<string> vAflowInName;
    vector<string> vFileLockName;

    if (pocc::structuresGenerated(directory_LIB)) { // CO20200624
      run_directory = true;
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " FOUND POCC" << endl;
      }
      if (!aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_ + ".pocc.preprocessing")) {
        if (aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_ + ".old")) {
          aurostd::file2file(directory_LIB + "/" + _AFLOWLOCK_ + ".old", directory_LIB + "/" + _AFLOWLOCK_ + ".pocc.preprocessing");
        } else if (aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_ + ".OLD")) {
          aurostd::file2file(directory_LIB + "/" + _AFLOWLOCK_ + ".OLD", directory_LIB + "/" + _AFLOWLOCK_ + ".pocc.preprocessing");
        } else if (aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_)) {
          aurostd::file2file(directory_LIB + "/" + _AFLOWLOCK_, directory_LIB + "/" + _AFLOWLOCK_ + ".pocc.preprocessing");
        }
      }
      AflowInName = _AFLOWIN_;
      if (aurostd::CompressFileExist(directory_LIB + "/" + AflowInName, stmp) && aurostd::IsCompressed(stmp)) {
        aurostd::DecompressFile(stmp);
      } // CO20210204 - fix aflow.in.xz
      vAflowInName.push_back(AflowInName); // AS20200915
      FileLockName = _AFLOWLOCK_;
      if (aurostd::CompressFileExist(directory_LIB + "/" + FileLockName, stmp) && aurostd::IsCompressed(stmp)) {
        aurostd::DecompressFile(stmp);
      } // CO20210204 - fix LOCK.xz
      vFileLockName.push_back(FileLockName); // AS20200915
      // ME20211008 - Check for AFLOW modules
      // APL
      apl_aflowin_found = apl::APL_Get_AflowInName(AflowInName, directory_LIB);
      if (apl_aflowin_found) { // CO20211125
        aurostd::RemoveFile(directory_LIB + "/" + POCC_FILE_PREFIX + POCC_APL_OUT_FILE);
        FileLockName = "LOCK.apl";
        if (aurostd::CompressFileExist(directory_LIB + "/" + FileLockName, stmp) && aurostd::IsCompressed(stmp)) {
          aurostd::DecompressFile(stmp);
        }
        vAflowInName.push_back(AflowInName);
        vFileLockName.push_back(FileLockName);
      }
    } else {
      AGL_functions::AGL_Get_AflowInName(AflowInName, directory_LIB, agl_aflowin_found); // CT20200713 Call function to find correct aflow.in file name  //CO20210204 - fix aflow.in.xz inside
      if (agl_aflowin_found) {
        run_directory = true;
        const vector<string> agl_and_ael_file_list{"aflow.agl.out",
                                                   "AGL.out",
                                                   "AGL_edos_gap_pressure.out",
                                                   "AGL_edos_gap_pressure.json",
                                                   "AGL_energies_temperature.out",
                                                   "AGL_energy.json",
                                                   "AGL_energy_structures.json",
                                                   "AGL_energy_volume.out",
                                                   "AGL_gibbs_energy_pT.out",
                                                   "AGL_Hugoniot.out",
                                                   "AGL_thermal_properties_temperature.out",
                                                   "AGL_THERMO",

                                                   "aflow.ael.out",
                                                   "AEL_Compliance_tensor.out",
                                                   "AEL_Elastic_constants.out",
                                                   "AEL_elastic_tensor.json",
                                                   "AEL_energy_structures.json"};
        for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
          for (const auto& file : agl_and_ael_file_list) {
            const string path = string(directory_LIB).append("/").append(file).append(XHOST.vext[iext]);
            aurostd::RemoveFile(path);
          }
        }
        // CORMAC FIX BELOW I NEED TO COMMENT TO RUN THE _AFLOWIN_AGL_DEFAULT_ containing only statics.
        FileLockName = "agl.LOCK"; // CO20210204 - reset from before
        if (aurostd::CompressFileExist(directory_LIB + "/" + FileLockName, stmp) && aurostd::IsCompressed(stmp)) {
          aurostd::DecompressFile(stmp);
        } // CO20210204 - fix LOCK.xz

        // AS20200904
        // save for later since we will need to loop among all possible submodules, i.e.
        // AGL, QHA,...
        vAflowInName.push_back(AflowInName); // AS20200904
        vFileLockName.push_back(FileLockName); // AS20200904
      } else {
        // Check for AEL input file
        AEL_functions::AEL_Get_AflowInName(AflowInName, directory_LIB, ael_aflowin_found); // CT20200715 Call function to find correct aflow.in file name  //CO20210204 - fix aflow.in.xz inside
        if (ael_aflowin_found) {
          run_directory = true;
          const vector<string> ael_file_list{"aflow.ael.out", "AEL_Compliance_tensor.out", "AEL_Elastic_constants.out", "AEL_elastic_tensor.json", "AEL_energy_structures.json"};
          for (size_t iext = 0; iext < XHOST.vext.size(); iext++) {
            for (const auto& file : ael_file_list) {
              const string path = string(directory_LIB).append("/").append(file).append(XHOST.vext[iext]);
              aurostd::RemoveFile(path);
            }
          }
          FileLockName = "ael.LOCK"; // CO20210204 - reset from before
          if (aurostd::CompressFileExist(directory_LIB + "/" + FileLockName, stmp) && aurostd::IsCompressed(stmp)) {
            aurostd::DecompressFile(stmp);
          } // CO20210204 - fix LOCK.xz

          // AS20200904
          // save for later since we will need to loop among all possible submodules,
          // i.e. AGL, QHA,...
          vAflowInName.push_back(AflowInName); // AS20200904
          vFileLockName.push_back(FileLockName); // AS20200904
        }
      }

      // ME20210901 BEGIN
      // Check for APL aflow.in
      apl_aflowin_found = apl::APL_Get_AflowInName(AflowInName, directory_LIB);
      if (apl_aflowin_found) {
        // Abort if there is an APL aflow.in file, but no PHPOSCAR (i.e., directory
        // was never run). Otherwise, lib2raw will run VASP.
        if (!aurostd::CompressFileExist(directory_LIB + "/" + DEFAULT_APL_PHPOSCAR_FILE)) {
          const string message = "PHPOSCAR not found. Cannot run directory " + directory_LIB + " for post-processing.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
        }
        run_directory = true;
        const vector<string> apl_file_list{DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_OUT_FILE,
                                           DEFAULT_APL_PHDOSCAR_FILE,
                                           DEFAULT_APL_PHKPOINTS_FILE,
                                           DEFAULT_APL_PHEIGENVAL_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_MSQRDISP_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_AAPL_GVEL_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HARMIFC_FILE,
                                           DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_POLAR_FILE};
        for (const auto& file : apl_file_list) {
          const string path = string(directory_LIB).append("/").append(file);
          aurostd::RemoveFile(path);
        }

        FileLockName = "LOCK.apl";
        if (aurostd::CompressFileExist(directory_LIB + "/" + FileLockName, stmp) && aurostd::IsCompressed(stmp)) {
          aurostd::DecompressFile(stmp);
        }
        // save for later since we will need to loop among all possible submodules
        vAflowInName.push_back(AflowInName);
        vFileLockName.push_back(FileLockName);
      }
      // ME20210901 END

      // AS20200902 BEGIN
      // Check for QHA input file
      qha_aflowin_found = apl::QHA_Get_AflowInName(AflowInName, directory_LIB);
      if (qha_aflowin_found) {
        run_directory = true;

        // clean QHA output files
        const vector<string> qha_file_regexs{DEFAULT_QHA_FILE_PREFIX + ".*", DEFAULT_QHA3P_FILE_PREFIX + ".*", DEFAULT_QHANP_FILE_PREFIX + ".*", DEFAULT_SCQHA_FILE_PREFIX + ".*", ".*qha.*png.*"};
        for (const auto& re : qha_file_regexs) {
          aurostd::RemoveFile(directory_LIB, std::regex(re));
        }

        FileLockName = "LOCK.qha"; // CO20210204 - reset from before
        if (aurostd::CompressFileExist(directory_LIB + "/" + FileLockName, stmp) && aurostd::IsCompressed(stmp)) {
          aurostd::DecompressFile(stmp);
        } // CO20210204 - fix LOCK.xz

        // AS20200904
        // save for later since we will need to loop among all possible submodules,
        // i.e. AGL, QHA,...
        vAflowInName.push_back(AflowInName); // AS20200904
        vFileLockName.push_back(FileLockName); // AS20200904
      }
      // AS20200902 END
    }
    // CT20200624 Calls functions to run AEL and AGL in postprocessing mode instead of executing an AFLOW run
    // CT20200624 This should help prevent VASP from running when performing postprocessing, since we go direct to AEL/AGL routines

    // KBIN::RUN_Directory() should be used instead of individual post-processing functions
    // otherwise we will have to program every single instance here
    // this also ensures we compress everything at the end
    if (run_directory) {
      _aflags aflags;
      aflags.Directory = directory_LIB;
      aflags.AFLOW_FORCE_RUN = true; // CO20201111 - force module run
      string aid_file; // CO20201111
      if (aurostd::FileExist(aflags.Directory + "/ALREADY_IN_DATABASE", aid_file) || aurostd::CompressFileExist(aflags.Directory + "/ALREADY_IN_DATABASE", aid_file)) { // CO20201111 - fix some broken in database
        aurostd::RemoveFile(aid_file);
      }

      // save originals
      const string _AFLOWIN_orig = _AFLOWIN_;
      const string _AFLOWLOCK_orig = _AFLOWLOCK_;

      // a sanity check
      if (vAflowInName.size() != vFileLockName.size()) {
        message << "vAflowInName.size()!=vFileLockName.size()";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }

      // AS20200904 loop over all detected earlier submodules
      for (size_t i = 0; i < vAflowInName.size(); i++) {
        // set env for RUN_Directory()
        // AS20200904 BEGIN
        //_AFLOWIN_=AflowInName;
        //_AFLOWLOCK_=FileLockName;

        _AFLOWIN_ = vAflowInName[i];
        _AFLOWLOCK_ = vFileLockName[i];
        cout << __AFLOW_FUNC__ << " Running KBIN::RUN_Directory() with aflow.in=";
        cout << vAflowInName[i] << " and LOCK=" << vFileLockName[i] << endl;
        // AS20200904 END

        // CO20200829 - because of LOCK and agl.LOCK in the same directory, sometimes we see LOCK.xz, we need to decompress
        // otherwise aflow can't run in the directory
        if (XHOST.vext.size() != XHOST.vzip.size()) { // CO20200829 - check for LOCK.xz and decompress first
          message << "XHOST.vext.size()!=XHOST.vzip.size()";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
        }
        for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // CO20200829 - check for LOCK.xz and decompress first // SKIP uncompressed
          if (aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_ + XHOST.vext[iext])) {
            aurostd::DecompressFile(directory_LIB + "/" + _AFLOWLOCK_ + XHOST.vext[iext]);
          }
        }
        if (aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_)) {
          aurostd::file2file(directory_LIB + "/" + _AFLOWLOCK_, directory_LIB + "/" + _AFLOWLOCK_ + ".run"); // keep original LOCK
        }
        KBIN::RUN_Directory(aflags);
      }

      // return to original
      _AFLOWIN_ = _AFLOWIN_orig;
      _AFLOWLOCK_ = _AFLOWLOCK_orig;

      // CO20210304 - patch, there are some problematic AEL/AGL runs having LOCK.run and no LOCK
      // these runs should have moved agl.LOCK, not LOCK
      // having no LOCK breaks all of the LOCK-reading functionality of lib2raw
      if (!aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_) && aurostd::FileExist(directory_LIB + "/" + _AFLOWLOCK_ + ".run")) {
        aurostd::file2file(directory_LIB + "/" + _AFLOWLOCK_ + ".run", directory_LIB + "/" + _AFLOWLOCK_); // LOCK.run->LOCK
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    return true;
  }
} // namespace aflowlib

//////////////////////////////////////////////////////////////////////////////////////
// CO20201111
#define _DEBUG_STOICH_FEATURES_ false
namespace aflowlib {
  void insertStoichStats(const vector<string> vstats, const xvector<double>& nspecies_xv, const xvector<double>& stoich_xv, vector<double>& vfeatures) {
    const bool LDEBUG = (false || _DEBUG_STOICH_FEATURES_ || XHOST.DEBUG);

    uint k = 0;
    uint l = 0;
    int index = 0;
    double d_tmp;
    vector<uint> vi_tmp;

    // check for NNN or AUROSTD_NAN
    bool has_NaN = false;
    for (index = nspecies_xv.lrows; index <= nspecies_xv.urows && !has_NaN; index++) {
      if (aurostd::isNaN(nspecies_xv[index])) {
        has_NaN = true;
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " nspecies_xv=" << nspecies_xv << endl;
    }
    for (k = 0; k < vstats.size(); k++) {
      if (vstats[k] == "min") {
        if (has_NaN) {
          vfeatures.push_back(NNN);
          continue;
        }
        vfeatures.emplace_back(aurostd::min(nspecies_xv));
      } else if (vstats[k] == "max") {
        if (has_NaN) {
          vfeatures.push_back(NNN);
          continue;
        }
        vfeatures.emplace_back(aurostd::max(nspecies_xv));
      } else if (vstats[k] == "range") {
        if (has_NaN) {
          vfeatures.push_back(NNN);
          continue;
        }
        vfeatures.emplace_back(aurostd::max(nspecies_xv) - aurostd::min(nspecies_xv));
      } else if (vstats[k] == "mean") {
        if (has_NaN) {
          vfeatures.push_back(NNN);
          continue;
        }
        vfeatures.emplace_back(aurostd::scalar_product(stoich_xv, nspecies_xv));
      } else if (vstats[k] == "dev") {
        if (has_NaN) {
          vfeatures.push_back(NNN);
          continue;
        }
        d_tmp = aurostd::scalar_product(stoich_xv, nspecies_xv); // mean
        vfeatures.emplace_back(aurostd::scalar_product(stoich_xv, aurostd::abs(nspecies_xv - d_tmp)));
      } else if (vstats[k] == "mode") { // property of most promiment species
        if (has_NaN) {
          vfeatures.push_back(NNN);
          continue;
        }
        d_tmp = aurostd::max(stoich_xv); // stoich_max
        vi_tmp.clear();
        for (index = stoich_xv.lrows; index <= stoich_xv.urows; index++) {
          if (aurostd::isequal(stoich_xv[index], d_tmp)) {
            vi_tmp.push_back(index);
          }
        }
        if (vi_tmp.size() == 1) {
          vfeatures.emplace_back(nspecies_xv[vi_tmp[0]]);
        } // easy case
        else {
          // take average
          d_tmp = 0;
          for (l = 0; l < vi_tmp.size(); l++) {
            d_tmp += nspecies_xv[vi_tmp[l]];
          }
          d_tmp /= (double) vi_tmp.size();
          vfeatures.push_back(d_tmp);
        }
      } else {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown statistic type: " + vstats[k], _RUNTIME_ERROR_);
      }
    }
  }

  void _aflowlib_entry::getStoichFeatures(vector<string>& vheaders, const string& e_props) {
    vector<double> vfeatures; // dummy
    return getStoichFeatures(vheaders, vfeatures, true, e_props);
  }
  void _aflowlib_entry::getStoichFeatures(vector<string>& vheaders, vector<double>& vfeatures, bool vheaders_only, const string& e_props) {
    // follows supplementary of 10.1038/npjcompumats.2016.28
    const bool LDEBUG = (false || _DEBUG_STOICH_FEATURES_ || XHOST.DEBUG);
    vheaders.clear();
    vfeatures.clear();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // headers
    stringstream tmp_ss;
    xelement::xelement xel;
    vector<xelement::xelement> vxel;
    uint i = 0;
    uint j = 0;
    uint k = 0;

    // L^p norms
    vector<uint> vp = {0, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    for (i = 0; i < vp.size(); i++) {
      aurostd::StringstreamClean(tmp_ss);
      tmp_ss << "stoich_norm_p_" << std::setfill('0') << std::setw(2) << vp[i];
      vheaders.emplace_back(tmp_ss.str());
    }

    // element-property-based
    // get which properties to average
    xel.populate(1); // dummy to get properties
    vector<string> vproperties_full;
    vector<string> vproperties;
    aurostd::string2tokens(e_props, vproperties_full, ",");
    vector<string> vstats = {"min", "max", "range", "mean", "dev", "mode"};
    // load up vheaders
    for (i = 0; i < vproperties_full.size(); i++) {
      if (vproperties_full[i] == "oxidation_states") {
        continue;
      } // skip this
      if (xel.getType(vproperties_full[i]) == "number" || xel.getType(vproperties_full[i]) == "numbers") {
        vproperties.emplace_back(vproperties_full[i]);

        if (xel.getType(vproperties.back()) == "number") {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " " << vproperties.back() << " is a number" << endl;
          }
          for (j = 0; j < vstats.size(); j++) {
            vheaders.emplace_back(vproperties.back() + "_stoich_" + vstats[j]);
          }
        } else if (xel.getType(vproperties.back()) == "numbers") {
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << " " << vproperties.back() << " are numbers" << endl;
          }
          if (vproperties.back() == "lattice_constants") {
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_a_stoich_" + vstats[j]);
            }
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_b_stoich_" + vstats[j]);
            }
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_c_stoich_" + vstats[j]);
            }
          } else if (vproperties.back() == "lattice_angles") {
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_alpha_stoich_" + vstats[j]);
            }
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_beta_stoich_" + vstats[j]);
            }
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_gamma_stoich_" + vstats[j]);
            }
          } else if (vproperties.back() == "oxidation_states_preferred") {
            for (j = 0; j < vstats.size(); j++) {
              vheaders.emplace_back(vproperties.back() + "_stoich_" + vstats[j]); // only use 0th oxidation_state_preferred
            }
          } else if (vproperties.back() == "energies_ionization") {
            for (k = 0; k < _ENERGIES_IONIZATION_MAX_AFLOWMACHL_; k++) {
              for (j = 0; j < vstats.size(); j++) {
                vheaders.emplace_back(vproperties.back() + "_" + aurostd::utype2string(k + 1) + "_stoich_" + vstats[j]);
              }
            }
          } else {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown numbers type: " + vproperties.back(), _RUNTIME_ERROR_);
          }
        }
      }
    }

    // valence (un)occupation
    vector<string> vorbitals = {"s", "p", "d", "f"};
    for (i = 0; i < vorbitals.size(); i++) {
      vheaders.emplace_back("valence_fraction_occupied_" + vorbitals[i]);
      vheaders.emplace_back("valence_fraction_unoccupied_" + vorbitals[i]);
    }

    // ionic character
    vheaders.emplace_back("formability_ionic");
    vector<string> vEN;
    for (i = 0; i < vproperties.size(); i++) {
      if (vproperties[i].find("electronegativity") != string::npos && xel.getUnits(vproperties[i]).empty()) { // must have NO units (goes in exp)
        vEN.emplace_back(vproperties[i]);
      }
    }
    for (i = 0; i < vEN.size(); i++) {
      vheaders.emplace_back("character_ionic_" + vEN[i] + "_max");
      vheaders.emplace_back("character_ionic_" + vEN[i] + "_mean");
    }

    if (LDEBUG) {
      for (i = 0; i < vheaders.size(); i++) {
        cerr << __AFLOW_FUNC__ << " vheaders[i=" << i << "]=\"" << vheaders[i] << "\"" << endl;
      }
    }

    if (vheaders_only) {
      return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // features

    const uint nspecies = vspecies.size();
    if (nspecies == 0) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "nspecies==0", _RUNTIME_ERROR_);
    }
    const xvector<double> nspecies_xv(nspecies);

    // L^p norms
    const xvector<double> stoich_xv(nspecies);
    for (j = 0; j < nspecies; j++) {
      stoich_xv[stoich_xv.lrows + j] = vcomposition[j] / natoms;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " stoich_xv=" << stoich_xv << endl;
    }
    for (i = 0; i < vp.size(); i++) {
      for (j = 0; j < nspecies; j++) {
        nspecies_xv[nspecies_xv.lrows + j] = std::pow(stoich_xv[stoich_xv.lrows + j], (double) vp[i]);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " nspecies_xv[\"stoich_norm_p_" + aurostd::utype2string(vp[i]) + "\"]=" << nspecies_xv << endl;
      }
      vfeatures.emplace_back(std::pow(sum(nspecies_xv), (vp[i] == 0 ? 1.0 : 1.0 / vp[i])));
    }

    // element-property-based
    // load up vxel
    int index = 0;
    int index_min = 0;
    int index_max = 0;
    for (j = 0; j < nspecies; j++) {
      vxel.emplace_back(vspecies[j]);
    }
    for (i = 0; i < vproperties.size(); i++) {
      if (xel.getType(vproperties[i]) == "number") {
        for (j = 0; j < nspecies; j++) {
          nspecies_xv[nspecies_xv.lrows + j] = vxel[j].getPropertyDouble(vproperties[i]);
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " nspecies_xv[\"" + vproperties[i] + "\"]=" << nspecies_xv << endl;
        }
        insertStoichStats(vstats, nspecies_xv, stoich_xv, vfeatures);
      }
      if (xel.getType(vproperties[i]) == "numbers") {
        if (vproperties[i] == "lattice_constants" || vproperties[i] == "lattice_angles") {
          index_min = 1;
          index_max = 3;
          for (index = index_min; index <= index_max; index++) {
            for (j = 0; j < nspecies; j++) {
              const xvector<double>& xvec = vxel[j].getPropertyXVectorDouble(vproperties[i]);
              nspecies_xv[nspecies_xv.lrows + j] = xvec[index];
            }
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " nspecies_xv[\"" + vproperties[i] + "_index_" + aurostd::utype2string(index) + "\"]=" << nspecies_xv << endl;
            }
            insertStoichStats(vstats, nspecies_xv, stoich_xv, vfeatures);
          }
        } else if (vproperties[i] == "oxidation_states_preferred" || vproperties[i] == "energies_ionization") {
          if (vproperties[i] == "oxidation_states_preferred") {
            index_min = 0;
            index_max = 0;
          } else if (vproperties[i] == "energies_ionization") {
            index_min = 0;
            index_max = _ENERGIES_IONIZATION_MAX_AFLOWMACHL_ - 1;
          } else {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown numbers property (vector): " + vproperties[i], _RUNTIME_ERROR_);
          }
          for (index = index_min; index <= index_max; index++) {
            for (j = 0; j < nspecies; j++) {
              const vector<double>& vec = vxel[j].getPropertyVectorDouble(vproperties[i]);
              nspecies_xv[nspecies_xv.lrows + j] = (index < (int) vec.size() ? vec[index] : NNN);
            }
            if (LDEBUG) {
              cerr << __AFLOW_FUNC__ << " nspecies_xv[\"" + vproperties[i] + "_index=" + aurostd::utype2string(index) + "\"]=" << nspecies_xv << endl;
            }
            insertStoichStats(vstats, nspecies_xv, stoich_xv, vfeatures);
          }
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown numbers property: " + vproperties[i], _RUNTIME_ERROR_);
        }
      }
    }

    // valence (un)occupation
    for (j = 0; j < nspecies; j++) {
      nspecies_xv[nspecies_xv.lrows + j] = vxel[j].getPropertyDouble("valence_std");
    }
    const double denom = aurostd::scalar_product(stoich_xv, nspecies_xv); // same for all quantities
    vector<double> vval_total = {2, 6, 10, 14};
    for (i = 0; i < vorbitals.size(); i++) {
      for (j = 0; j < nspecies; j++) {
        nspecies_xv[nspecies_xv.lrows + j] = vxel[j].getPropertyDouble("valence_" + vorbitals[i]);
      } // populate with orbital occupation
      vfeatures.emplace_back(aurostd::scalar_product(stoich_xv, nspecies_xv) / denom); // occupied
      for (j = 0; j < nspecies; j++) {
        nspecies_xv[nspecies_xv.lrows + j] = vval_total[i] - nspecies_xv[nspecies_xv.lrows + j];
      } // has occupied inside already
      vfeatures.emplace_back(aurostd::scalar_product(stoich_xv, nspecies_xv) / denom); // unoccupied
    }

    // ionic character
    // ionic formability
    bool formability_ionic = false;
    bool has_NaN = false;
    vector<int> nspecies_v;
    aurostd::xcombos xc;
    for (j = 0; j < nspecies && !has_NaN; j++) {
      const vector<double>& oxidation_states = vxel[j].getPropertyVectorDouble("oxidation_states");
      if (oxidation_states.empty()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No oxidation states found for element: " + vxel[j].symbol, _RUNTIME_ERROR_);
      } // should have NNN
      if (oxidation_states.size() == 1 && aurostd::isNaN(oxidation_states[0])) {
        has_NaN = true;
      }
      nspecies_v.push_back((int) oxidation_states.size());
    }
    if (has_NaN) {
      cerr << __AFLOW_FUNC__ << " has NaN" << endl;
    }
    if (!has_NaN) {
      const xvector<double> natoms_xv(natoms);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " oxidation_states_count=" << aurostd::joinWDelimiter(nspecies_v, ",") << endl;
      }
      xc.reset(nspecies_v, 'E');
      while (xc.increment() && !formability_ionic) {
        const vector<int>& indices = xc.getCombo();
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " indices=" << aurostd::joinWDelimiter(indices, ",") << endl;
        }
        i = 0;
        for (j = 0; j < nspecies && !has_NaN; j++) {
          const vector<double>& oxidation_states = vxel[j].getPropertyVectorDouble("oxidation_states");
          for (k = 0; k < vcomposition[j]; k++) {
            if (aurostd::isNaN(oxidation_states[indices[j]])) {
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Found NaN among populated oxidation_states", _RUNTIME_ERROR_);
            }
            natoms_xv[natoms_xv.lrows + (i++)] = oxidation_states[indices[j]];
          }
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " natoms_xv[\"oxidation_states\"]=" << natoms_xv << endl;
        }
        if (aurostd::isequal(aurostd::sum(natoms_xv), 0.0)) {
          formability_ionic = true;
        }
      }
    }
    vfeatures.push_back((!has_NaN && formability_ionic ? 1 : 0));
    // character_ionic_max and _mean
    const xvector<double> pairs_xv(aurostd::nCk((int) nspecies, 2)); // electronegativities
    const xvector<double> pairs2_xv(aurostd::nCk((int) nspecies, 2)); // stoich
    for (i = 0; i < vEN.size(); i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " EN=" << vEN[i] << endl;
      }
      xc.reset(nspecies, 2);
      k = 0;
      while (xc.increment()) {
        const vector<int>& indices = xc.getIndices();
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " indices=" << aurostd::joinWDelimiter(indices, ",") << endl;
        }
        pairs_xv[pairs_xv.lrows + k] = 1.0 - std::exp(-0.25 * (std::pow(vxel[indices[0]].getPropertyDouble(vEN[i]) - vxel[indices[1]].getPropertyDouble(vEN[i]), 2.0)));
        pairs2_xv[pairs2_xv.lrows + k] = stoich_xv[stoich_xv.lrows + indices[0]] * stoich_xv[stoich_xv.lrows + indices[1]];
        k++;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " pairs_xv[\"" + vEN[i] + "\"]=" << pairs_xv << endl;
      }
      vfeatures.emplace_back(aurostd::max(pairs_xv));
      vfeatures.emplace_back(aurostd::scalar_product(pairs_xv, pairs2_xv));
    }

    if (vheaders.size() != vfeatures.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "vheaders.size()!=vfeatures.size()", _RUNTIME_ERROR_);
    }

    if (LDEBUG) {
      for (i = 0; i < vheaders.size(); i++) {
        cerr << __AFLOW_FUNC__ << " " << vheaders[i] << "=" << vfeatures[i] << endl;
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }
} // namespace aflowlib

#endif //  _AFLOWLIB_LIBRARIES_CPP_
// ***************************************************************************
