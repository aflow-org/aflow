// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#include "aflow_xhost.h"

#include <cstddef>
#include <string>
#include <vector>

#include <pthread.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"

using std::string;
using std::vector;

// ***************************************************************************
// ***************************************************************************
// _XHOST
// look into aflow.h for the definitions
// constructors

_XHOST::_XHOST() {  // constructor PUBLIC
  PID = 0;
  TID = 0;  // CO20200502 - threadID
  ostrPID.clear();
  ostrPID.str(string(""));
  ostrTID.clear();
  ostrTID.str(string(""));  // CO20200502 - threadID
  sPID = "";
  sTID = "";
  showPID = false;
  showTID = false;
  QUIET = false;
  QUIET_GLOBAL = false; // CO20220630
  QUIET_CERR = false; // extra quiet SC20210617
  QUIET_COUT = false; // extra quiet SC20210617
  LOGGER_WHITELIST.clear();  // HE+ME20220305
  LOGGER_BLACKLIST.clear();  // HE+ME20220305
  TEST = false;
  DEBUG = false;
  MPI = false;
  GENERATE_AFLOWIN_ONLY = false;  // CT20180719
  POSTPROCESS = false;  // CO20200624
  ARUN_POSTPROCESS = false;  // CT20181212
  AVOID_RUNNING_VASP = false;  // CO20200624
  PSEUDOPOTENTIAL_GENERATOR = false; // SC20200327
  hostname = "";
  machine_type = "";
  tmpfs = "";
  user = "";
  group = "";
  home = "";
  shell = "";
  progname = "aflow";
  Find_Parameters = "";
  sensors_allowed = true;
  argv.clear();
  AFLOW_MATERIALS_SERVER = "";
  AFLOW_WEB_SERVER = "";
  RAM = 0.0;
  RAM_MB = 0.0;
  RAM_GB = 0.0;
  CPU_Cores = 0;
  CPU_active = 1;  // ME20220130 - count main thread
  CPU_Model = "";
  CPU_MHz = "";
  vTemperatureCore.clear();
  Time_starting = 0.0;
  Time_now = 0.0;
  Date = 0;
  Day = "";
  Month = "";
  Year = "";
  Copyright_Years = "";
  // PTHREADS_FLAG=false;
  // PTHREADS_MAX=0;
  // PTHREADS_RUNNING=0;
  // thread.clear();
  // iret.clear();
  // thread_busy.clear();
  vcmd.clear();
  maxmem = 100.00;
  AFLOW_RUNDIRflag = false;
  AFLOW_MULTIflag = false;
  AFLOW_RUNXflag = false;
  AFLOW_RUNXnumber = 0;
  is_PBS = false;
  PBS_NUM_PPN = 0;
  PBS_NNODES = 0;
  is_SLURM = false;
  SLURM_CPUS_ON_NODE = 0;
  SLURM_NNODES = 0;
  SLURM_NTASKS = 0;
  is_MACHINE_FULTON_MARYLOU = false;
  vGlobal_uint.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vGlobal_uint.push_back(0);
  }
  vGlobal_string.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vGlobal_string.emplace_back("");
  }
  vvGlobal_string.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vvGlobal_string.emplace_back();
  }
  vvLIBS.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vvLIBS.emplace_back();
  }
  vflag_aflow.clear();
  vflag_pflow.clear();
  vflag_outreach.clear();
  vflag_control.clear();
  vschema.clear();
  vschema_internal.clear(); // ME20220208
  XHOST_LIBRARY_LIB0 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB1 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB2 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB3 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB4 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB5 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB6 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB7 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB8 = LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB9 = LIBRARY_NOTHING;
  XHOST_LIBRARY_ICSD = LIBRARY_NOTHING;
  XHOST_vLibrary_ICSD.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    XHOST_vLibrary_ICSD.emplace_back("");  // needs some initialization
  }
  // extensions
  vcat.clear();
  vcat.emplace_back("cat");
  vcat.emplace_back("bzcat");
  vcat.emplace_back("xzcat");
  vcat.emplace_back("gzcat");
  vext.clear();
  vext.emplace_back("");
  vext.emplace_back(".bz2");
  vext.emplace_back(".xz");
  vext.emplace_back(".gz");
  vzip.clear();
  vzip.emplace_back("");
  vzip.emplace_back("bzip2");
  vzip.emplace_back("xz");
  vzip.emplace_back("gzip");
  // AFLOWRC
  aflowrc_filename = "";    // AFLOWRC
  aflowrc_content = "";    // AFLOWRC
  vaflowrc.clear();   // AFLOWRC
  adefault.clear();    // AFLOWRC
  // AFLOWSYM
  SKEW_TEST = false; // DX20171019
  SKEW_TOL = AUROSTD_NAN; // DX20171019
  // xstructure
  READ_SPIN_FROM_ATOMLABEL = false; // SD20220316 - spin and pp label can conflict with one another
  // WEB
  //[CO20200404 - overload with --www]WEB_MODE=false; //CO20190402
};

_XHOST::~_XHOST() { // destructor PUBLIC
  free();
}

void _XHOST::copy(const _XHOST& b) { // copy PRIVATE
  PID = b.PID;
  TID = b.TID;  // CO20200502 - threadID
  ostrPID.clear();
  ostrPID.str(string(""));
  ostrPID << b.ostrPID.str();
  ostrTID.clear();
  ostrTID.str(string(""));
  ostrTID << b.ostrTID.str(); // CO20200502 - threadID
  sPID = b.sPID;
  sTID = b.sTID;
  showPID = b.showPID;
  showTID = b.showTID;
  QUIET = b.QUIET;
  QUIET_GLOBAL = b.QUIET_GLOBAL;  // CO20220630
  QUIET_CERR = b.QUIET_CERR; // extra quiet SC20210617
  QUIET_COUT = b.QUIET_COUT; // extra quiet SC20210617
  LOGGER_WHITELIST = b.LOGGER_WHITELIST;  // HE+ME20220305
  LOGGER_BLACKLIST = b.LOGGER_BLACKLIST;  // HE+ME20220305
  TEST = b.TEST;
  DEBUG = b.DEBUG;
  MPI = b.MPI;
  GENERATE_AFLOWIN_ONLY = b.GENERATE_AFLOWIN_ONLY;  // CT20180719
  POSTPROCESS = b.POSTPROCESS;  // CO20200624
  ARUN_POSTPROCESS = b.ARUN_POSTPROCESS;  // CT20181212
  AVOID_RUNNING_VASP = b.AVOID_RUNNING_VASP;  // CO20200624
  PSEUDOPOTENTIAL_GENERATOR = b.PSEUDOPOTENTIAL_GENERATOR; // SC20200327
  hostname = b.hostname;
  machine_type = b.machine_type;
  tmpfs = b.tmpfs;
  user = b.user;
  group = b.group;
  home = b.home;
  progname = b.progname;
  Find_Parameters = b.Find_Parameters;
  sensors_allowed = b.sensors_allowed;
  argv.clear();
  for (size_t i = 0; i < b.argv.size(); i++) {
    argv.push_back(b.argv[i]);
  }
  AFLOW_MATERIALS_SERVER = b.AFLOW_MATERIALS_SERVER;
  AFLOW_WEB_SERVER = b.AFLOW_WEB_SERVER;
  RAM = b.RAM;
  RAM_MB = b.RAM_MB;
  RAM_GB = b.RAM_GB;
  CPU_Cores = b.CPU_Cores;
  CPU_active = b.CPU_active;  // ME20220130
  CPU_Model = b.CPU_Model;
  CPU_MHz = b.CPU_MHz;
  vTemperatureCore.clear();
  for (size_t i = 0; i < b.vTemperatureCore.size(); i++) {
    vTemperatureCore.push_back(b.vTemperatureCore[i]);
  }
  Time_starting = b.Time_starting;
  Time_now = b.Time_now;
  Date = b.Date;
  Day = b.Day;
  Month = b.Month;
  Year = b.Year;
  Copyright_Years = b.Copyright_Years;
  // PTHREADS_FLAG=b.PTHREADS_FLAG;
  // PTHREADS_MAX=b.PTHREADS_MAX;
  // PTHREADS_RUNNING=b.PTHREADS_RUNNING;
  // thread.clear();for(size_t i=0;i<b.thread.size();i++) thread.push_back(b.thread.at(i));
  // iret.clear();for(size_t i=0;i<b.iret.size();i++) iret.push_back(b.iret.at(i));
  // thread_busy.clear();for(size_t i=0;i<b.thread.size();i++) thread.push_back(b.thread.at(i));
  vcmd.clear();
  for (size_t i = 0; i < b.vcmd.size(); i++) {
    vcmd.push_back(b.vcmd[i]);
  }
  maxmem = b.maxmem;
  AFLOW_RUNDIRflag = b.AFLOW_RUNDIRflag;
  AFLOW_MULTIflag = b.AFLOW_MULTIflag;
  AFLOW_RUNXflag = b.AFLOW_RUNXflag;
  AFLOW_RUNXnumber = b.AFLOW_RUNXnumber;
  is_PBS = b.is_PBS;
  PBS_NUM_PPN = b.PBS_NUM_PPN;
  PBS_NNODES = b.PBS_NNODES;
  is_SLURM = b.is_SLURM;
  SLURM_CPUS_ON_NODE = b.SLURM_CPUS_ON_NODE;
  SLURM_NNODES = b.SLURM_NNODES;
  SLURM_NTASKS = b.SLURM_NTASKS;
  is_MACHINE_FULTON_MARYLOU = b.is_MACHINE_FULTON_MARYLOU;
  vGlobal_uint.clear();
  for (size_t i = 0; i < b.vGlobal_uint.size(); i++) {
    vGlobal_uint.push_back(b.vGlobal_uint[i]);
  }
  vGlobal_string.clear();
  for (size_t i = 0; i < b.vGlobal_string.size(); i++) {
    vGlobal_string.push_back(b.vGlobal_string[i]);
  }
  vvGlobal_string.clear();
  for (size_t i = 0; i < b.vvGlobal_string.size(); i++) {
    vvGlobal_string.push_back(b.vvGlobal_string[i]);
  }
  vvLIBS.clear();
  for (size_t i = 0; i < b.vvLIBS.size(); i++) {
    vvLIBS.push_back(b.vvLIBS[i]);
  }
  vflag_aflow = b.vflag_aflow;
  vflag_pflow = b.vflag_pflow;
  vflag_outreach = b.vflag_outreach;
  vflag_control = b.vflag_control;
  vschema = b.vschema;
  vschema_internal = b.vschema_internal;  // ME20220208
  // extensions
  vcat.clear();
  for (size_t i = 0; i < b.vcat.size(); i++) {
    vcat.push_back(b.vcat[i]);
  }
  vext.clear();
  for (size_t i = 0; i < b.vext.size(); i++) {
    vext.push_back(b.vext[i]);
  }
  vzip.clear();
  for (size_t i = 0; i < b.vzip.size(); i++) {
    vzip.push_back(b.vzip[i]);
  }
  // AFLOWRC
  aflowrc_filename = b.aflowrc_filename;    // AFLOWRC
  aflowrc_content = b.aflowrc_content;    // AFLOWRC
  vaflowrc.clear();
  for (size_t i = 0; i < b.vaflowrc.size(); i++) {
    vaflowrc.push_back(b.vaflowrc[i]);   // AFLOWRC
  }
  adefault.clear();
  adefault = b.adefault; // AFLOWRC
  // AFLOWSYM
  SKEW_TEST = b.SKEW_TEST; // DX20171019
  SKEW_TOL = b.SKEW_TOL; // DX20171019
  // xstructure
  READ_SPIN_FROM_ATOMLABEL = b.READ_SPIN_FROM_ATOMLABEL; // SD20220316
  // WEB
  //[CO20200404 - overload with --www]WEB_MODE=b.WEB_MODE;  //CO20190402
}

const _XHOST& _XHOST::operator=(const _XHOST& b) {  // operator= PUBLIC
  if (this != &b) {
    free();
    copy(b);
  }
  return *this;
}

//_XHOST::_XHOST(const _XHOST& b) { // copy PUBLIC
////  free();*this=b;
// copy(b);
// }

void _XHOST::free() { // free PRIVATE
  ostrPID.clear();
  ostrPID.str(string(""));
  ostrTID.clear();
  ostrTID.str(string(""));  // CO20200502 - threadID
  vTemperatureCore.clear();
  // thread.clear();
  // iret.clear();
  // thread_busy.clear();
  vcmd.clear();
  vGlobal_uint.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vGlobal_uint.push_back(0);
  }
  vGlobal_string.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vGlobal_string.emplace_back("");
  }
  vvGlobal_string.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vvGlobal_string.emplace_back();
  }
  vvLIBS.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    vvLIBS.emplace_back();
  }
  XHOST_vLibrary_ICSD.clear();
  for (uint i = 0; i < XHOST_vGlobal_MAX; i++) {
    XHOST_vLibrary_ICSD.emplace_back("");  // needs some initialization
  }
  vflag_aflow.clear();
  vflag_pflow.clear();
  vflag_outreach.clear();
  vflag_control.clear();
  vschema.clear();
  vschema_internal.clear();  // ME20220208
  // extensions
  vcat.clear();
  vcat.emplace_back("cat");
  vcat.emplace_back("bzcat");
  vcat.emplace_back("xzcat");
  vcat.emplace_back("gzcat");
  vext.clear();
  vext.emplace_back("");
  vext.emplace_back(".bz2");
  vext.emplace_back(".xz");
  vext.emplace_back(".gz");
  vzip.clear();
  vzip.emplace_back("");
  vzip.emplace_back("bzip2");
  vzip.emplace_back("xz");
  vzip.emplace_back("gzip");
  // AFLOWRC
  aflowrc_filename.clear();    // AFLOWRC
  aflowrc_content.clear();    // AFLOWRC
  vaflowrc.clear();   // AFLOWRC
  adefault.clear();  // AFLOWRC
}

void _XHOST::clear() {  // clear PRIVATE
  _XHOST const XHOST_temp;
  copy(XHOST_temp);
}

// other public functions of XHOST

pthread_mutex_t mutex_XAFLOW_XHOST = PTHREAD_MUTEX_INITIALIZER;

std::string _XHOST::command(const string& command) {
  string _command = command;
#ifdef _MACOSX_
  if (command == "beep") return string("echo -ne '\007'");
#endif
  // first check for EXACT match in vcmds  //CO20180705
  if (command == "aflow_data" && aurostd::FileExist("./aflow_data")) {
    return "./aflow_data";
  } // CO20180705 - hack for developers, prefer ./aflow_data over aflow_data in PATH, as it is probably newer (development)
  for (size_t i = 0; i < vcmd.size(); i++) {
    if (vcmd[i] == command) {
      return vcmd[i];
    }// found before.. only == otherwise cat gets confused with bzcat
  }
  // next check if we can find the command in PATH or somewher else common //CO20180705
  if (aurostd::IsCommandAvailableModify(_command)) {
    pthread_mutex_lock(&mutex_XAFLOW_XHOST);
    // pthread_mutex_unlock(&mutex_XAFLOW_XHOST);
    vcmd.push_back(_command);
    pthread_mutex_unlock(&mutex_XAFLOW_XHOST);
    return _command;
  } // found and added
  // CO20180705 START
  // requested command is "aflow_data", and we have ./aflow_data in our vcmd's
  // we need to strip vcmd's to basename and check
  vector<string> path_parts;
  string basename;
  for (size_t i = 0; i < vcmd.size(); i++) {
    aurostd::string2tokens(vcmd[i], path_parts, "/", true);  // consecutive
    if (!path_parts.empty()) {
      basename = path_parts[path_parts.size() - 1];
      if (basename == command) {
        return vcmd[i];
      }
    }
  }
  // CO20180705 STOP
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "command=" + command + " not found", _INPUT_MISSING_);
  return string();
}

bool _XHOST::is_command(const string& command) {
  try {
    const string path = _XHOST::command(command);
    return !path.empty();
  }    // CO20190629 - using xerror now
  catch (aurostd::xerror& excpt) {
    return false;
  }  // CO20190629 - using xerror now
  //[CO20190629 - using xerror now]return !_XHOST::command(command).empty(); //CO20180705
}
