// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AFLOW_CLASSES_CPP
#define _AFLOW_CLASSES_CPP
#include "aflow.h"

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _XHOST
// look into aflow.h for the definitions
// constructors

_XHOST::_XHOST() {  // constructor PUBLIC
  PID=0;
  TID=0;  //CO20200502 - threadID
  ostrPID.clear();ostrPID.str(string(""));
  ostrTID.clear();ostrTID.str(string(""));  //CO20200502 - threadID
  sPID="";
  sTID="";
  showPID=FALSE;
  showTID=FALSE;
  QUIET=FALSE;
  QUIET_GLOBAL=FALSE; //CO20220630
  QUIET_CERR=FALSE; // extra quiet SC20210617
  QUIET_COUT=FALSE; // extra quiet SC20210617
  LOGGER_WHITELIST.clear();  //HE+ME20220305
  LOGGER_BLACKLIST.clear();  //HE+ME20220305
  TEST=FALSE;
  DEBUG=FALSE;
  MPI=FALSE;
  GENERATE_AFLOWIN_ONLY=FALSE;  //CT20180719
  POSTPROCESS=FALSE;  //CO20200624
  ARUN_POSTPROCESS=FALSE;  //CT20181212
  AVOID_RUNNING_VASP=FALSE;  //CO20200624
  PSEUDOPOTENTIAL_GENERATOR=FALSE; //SC20200327
  hostname="";
  machine_type="";
  tmpfs="";
  user="";
  group="";
  home="";
  shell="";
  progname="aflow";
  Find_Parameters="";
  sensors_allowed=TRUE;
  argv.clear();
  AFLOW_MATERIALS_SERVER="";
  AFLOW_WEB_SERVER="";
  RAM=0.0;
  RAM_MB=0.0;
  RAM_GB=0.0;
  CPU_Cores=0;
  CPU_active=1;  //ME20220130 - count main thread
  CPU_Model="";
  CPU_MHz="";
  vTemperatureCore.clear();
  Time_starting=0.0;
  Time_now=0.0;
  Date=0;
  Day="";
  Month="";
  Year="";
  Copyright_Years="";
  // PTHREADS_FLAG=FALSE;
  // PTHREADS_MAX=0;
  // PTHREADS_RUNNING=0;
  // thread.clear();
  // iret.clear();
  // thread_busy.clear();
  vcmd.clear();
  maxmem=100.00;
  AFLOW_RUNDIRflag=FALSE;
  AFLOW_MULTIflag=FALSE;
  AFLOW_RUNXflag=FALSE;
  AFLOW_RUNXnumber=0;
  is_PBS=FALSE;
  PBS_NUM_PPN=0;
  PBS_NNODES=0;
  is_SLURM=FALSE;
  SLURM_CPUS_ON_NODE=0;
  SLURM_NNODES=0;
  SLURM_NTASKS=0;
  is_MACHINE_FULTON_MARYLOU=FALSE;
  APENNSY_USE_SERVER=FALSE;
  APENNSY_USE_LIBRARY=FALSE;
  APENNSY_SERVER_AFLOWLIB_ORG=FALSE;
  vGlobal_uint.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vGlobal_uint.push_back(0); 
  vGlobal_string.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vGlobal_string.push_back(""); 
  vvGlobal_string.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vvGlobal_string.push_back(vector<string>()); 
  vvLIBS.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vvLIBS.push_back(vector<string>()); 
  // [OBSOLETE] vLibrary_ICSD.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vLibrary_ICSD.push_back(""); 
  // [OBSOLETE] vLibrary_ICSD_ALL.clear(); // for(uint i=0;i<XHOST_vGlobal_MAX;i++) vLibrary_ICSD_ALL.push_back(""); 
  // [OBSOLETE] Library_ICSD_ALL="";
  vflag_aflow.clear();
  vflag_pflow.clear();
  vflag_apennsy.clear();
  vflag_outreach.clear();
  vflag_control.clear();
  vschema.clear();
  vschema_internal.clear(); //ME20220208
  XHOST_LIBRARY_LIB0=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB1=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB2=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB3=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB4=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB5=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB6=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB7=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB8=LIBRARY_NOTHING;
  XHOST_LIBRARY_LIB9=LIBRARY_NOTHING;
  XHOST_LIBRARY_ICSD=LIBRARY_NOTHING;
  XHOST_vLibrary_ICSD.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) XHOST_vLibrary_ICSD.push_back("");  // needs some initialization
  // extensions
  vcat.clear();vcat.push_back("cat");vcat.push_back("bzcat"); vcat.push_back("xzcat");vcat.push_back("gzcat");
  vext.clear();vext.push_back("");   vext.push_back(".bz2");  vext.push_back(".xz");  vext.push_back(".gz");
  vzip.clear();vzip.push_back("");   vzip.push_back("bzip2"); vzip.push_back("xz");   vzip.push_back("gzip");
  // AFLOWRC
  aflowrc_filename="";    // AFLOWRC
  aflowrc_content="";    // AFLOWRC
  vaflowrc.clear();   // AFLOWRC
  adefault.clear();    // AFLOWRC
  // AFLOWSYM
  SKEW_TEST=FALSE; //DX20171019
  SKEW_TOL=AUROSTD_NAN; //DX20171019
  // xstructure
  READ_SPIN_FROM_ATOMLABEL=FALSE; //SD20220316 - spin and pp label can conflict with one another
  // WEB
  //[CO20200404 - overload with --www]WEB_MODE=FALSE; //CO20190402
};

_XHOST::~_XHOST() { // destructor PUBLIC
  free();
}

void _XHOST::copy(const _XHOST& b) { // copy PRIVATE
  PID=b.PID;
  TID=b.TID;  //CO20200502 - threadID
  ostrPID.clear();ostrPID.str(string(""));ostrPID << b.ostrPID.str();
  ostrTID.clear();ostrTID.str(string(""));ostrTID << b.ostrTID.str(); //CO20200502 - threadID
  sPID=b.sPID;
  sTID=b.sTID;
  showPID=b.showPID;
  showTID=b.showTID;
  QUIET=b.QUIET;
  QUIET_GLOBAL=b.QUIET_GLOBAL;  //CO20220630
  QUIET_CERR=b.QUIET_CERR; // extra quiet SC20210617
  QUIET_COUT=b.QUIET_COUT; // extra quiet SC20210617
  LOGGER_WHITELIST = b.LOGGER_WHITELIST;  //HE+ME20220305
  LOGGER_BLACKLIST = b.LOGGER_BLACKLIST;  //HE+ME20220305
  TEST=b.TEST;
  DEBUG=b.DEBUG;
  MPI=b.MPI;
  GENERATE_AFLOWIN_ONLY=b.GENERATE_AFLOWIN_ONLY;  //CT20180719
  POSTPROCESS=b.POSTPROCESS;  //CO20200624
  ARUN_POSTPROCESS=b.ARUN_POSTPROCESS;  //CT20181212
  AVOID_RUNNING_VASP=b.AVOID_RUNNING_VASP;  //CO20200624
  PSEUDOPOTENTIAL_GENERATOR=b.PSEUDOPOTENTIAL_GENERATOR; //SC20200327
  hostname=b.hostname;
  machine_type=b.machine_type;
  tmpfs=b.tmpfs;
  user=b.user;
  group=b.group;
  home=b.home;
  progname=b.progname;
  Find_Parameters=b.Find_Parameters;
  sensors_allowed=b.sensors_allowed;
  argv.clear();for(uint i=0;i<b.argv.size();i++) argv.push_back(b.argv.at(i));
  AFLOW_MATERIALS_SERVER=b.AFLOW_MATERIALS_SERVER;
  AFLOW_WEB_SERVER=b.AFLOW_WEB_SERVER;
  RAM=b.RAM;
  RAM_MB=b.RAM_MB;
  RAM_GB=b.RAM_GB;
  CPU_Cores=b.CPU_Cores;
  CPU_active=b.CPU_active;  //ME20220130
  CPU_Model=b.CPU_Model;
  CPU_MHz=b.CPU_MHz;
  vTemperatureCore.clear();for(uint i=0;i<b.vTemperatureCore.size();i++) vTemperatureCore.push_back(b.vTemperatureCore.at(i));
  Time_starting=b.Time_starting;
  Time_now=b.Time_now;
  Date=b.Date;
  Day=b.Day;
  Month=b.Month;
  Year=b.Year;
  Copyright_Years=b.Copyright_Years;
  // PTHREADS_FLAG=b.PTHREADS_FLAG;
  // PTHREADS_MAX=b.PTHREADS_MAX;
  // PTHREADS_RUNNING=b.PTHREADS_RUNNING;
  // thread.clear();for(uint i=0;i<b.thread.size();i++) thread.push_back(b.thread.at(i));
  // iret.clear();for(uint i=0;i<b.iret.size();i++) iret.push_back(b.iret.at(i));
  // thread_busy.clear();for(uint i=0;i<b.thread.size();i++) thread.push_back(b.thread.at(i));
  vcmd.clear();for(uint i=0;i<b.vcmd.size();i++) vcmd.push_back(b.vcmd.at(i));
  maxmem=b.maxmem;
  AFLOW_RUNDIRflag=b.AFLOW_RUNDIRflag;
  AFLOW_MULTIflag=b.AFLOW_MULTIflag;
  AFLOW_RUNXflag=b.AFLOW_RUNXflag;
  AFLOW_RUNXnumber=b.AFLOW_RUNXnumber;
  is_PBS=b.is_PBS;
  PBS_NUM_PPN=b.PBS_NUM_PPN;
  PBS_NNODES=b.PBS_NNODES;
  is_SLURM=b.is_SLURM;
  SLURM_CPUS_ON_NODE=b.SLURM_CPUS_ON_NODE;
  SLURM_NNODES=b.SLURM_NNODES;
  SLURM_NTASKS=b.SLURM_NTASKS;
  is_MACHINE_FULTON_MARYLOU=b.is_MACHINE_FULTON_MARYLOU;
  APENNSY_USE_SERVER=b.APENNSY_USE_SERVER;
  APENNSY_USE_LIBRARY=b.APENNSY_USE_LIBRARY;
  APENNSY_SERVER_AFLOWLIB_ORG=b.APENNSY_SERVER_AFLOWLIB_ORG;
  vGlobal_uint.clear();for(uint i=0;i<b.vGlobal_uint.size();i++) vGlobal_uint.push_back(b.vGlobal_uint.at(i));
  vGlobal_string.clear();for(uint i=0;i<b.vGlobal_string.size();i++) vGlobal_string.push_back(b.vGlobal_string.at(i));
  vvGlobal_string.clear();for(uint i=0;i<b.vvGlobal_string.size();i++) vvGlobal_string.push_back(b.vvGlobal_string.at(i));
  vvLIBS.clear();for(uint i=0;i<b.vvLIBS.size();i++) vvLIBS.push_back(b.vvLIBS.at(i));
  // [OBSOLETE] vLibrary_ICSD.clear();for(uint i=0;i<b.vLibrary_ICSD.size();i++) vLibrary_ICSD.push_back(b.vLibrary_ICSD.at(i));
  // [OBSOLETE] vLibrary_ICSD_ALL.clear();for(uint i=0;i<b.vLibrary_ICSD_ALL.size();i++) vLibrary_ICSD_ALL.push_back(b.vLibrary_ICSD_ALL.at(i));
  // [OBSOLETE] Library_ICSD_ALL=b.Library_ICSD_ALL;
  vflag_aflow=b.vflag_aflow;
  vflag_pflow=b.vflag_pflow;
  vflag_apennsy=b.vflag_apennsy;
  vflag_outreach=b.vflag_outreach;
  vflag_control=b.vflag_control;
  vschema=b.vschema;
  vschema_internal = b.vschema_internal;  //ME20220208
  // extensions
  vcat.clear();for(uint i=0;i<b.vcat.size();i++) vcat.push_back(b.vcat.at(i));
  vext.clear();for(uint i=0;i<b.vext.size();i++) vext.push_back(b.vext.at(i));
  vzip.clear();for(uint i=0;i<b.vzip.size();i++) vzip.push_back(b.vzip.at(i));
  // AFLOWRC
  aflowrc_filename=b.aflowrc_filename;    // AFLOWRC
  aflowrc_content=b.aflowrc_content;    // AFLOWRC
  vaflowrc.clear();for(uint i=0;i<b.vaflowrc.size();i++) vaflowrc.push_back(b.vaflowrc.at(i));   // AFLOWRC
  adefault.clear(); adefault=b.adefault; // AFLOWRC
  // AFLOWSYM
  SKEW_TEST=b.SKEW_TEST; //DX20171019
  SKEW_TOL=b.SKEW_TOL; //DX20171019
  // xstructure
  READ_SPIN_FROM_ATOMLABEL=b.READ_SPIN_FROM_ATOMLABEL; //SD20220316
  // WEB
  //[CO20200404 - overload with --www]WEB_MODE=b.WEB_MODE;  //CO20190402
}

const _XHOST& _XHOST::operator=(const _XHOST& b) {  // operator= PUBLIC
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

//_XHOST::_XHOST(const _XHOST& b) { // copy PUBLIC
////  free();*this=b;
//copy(b);
//}

void _XHOST::free() { // free PRIVATE
  ostrPID.clear();ostrPID.str(string(""));
  ostrTID.clear();ostrTID.str(string(""));  //CO20200502 - threadID
  vTemperatureCore.clear();
  // thread.clear();
  // iret.clear();
  // thread_busy.clear();
  vcmd.clear();
  vGlobal_uint.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vGlobal_uint.push_back(0); 
  vGlobal_string.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vGlobal_string.push_back(""); 
  vvGlobal_string.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vvGlobal_string.push_back(vector<string>()); 
  vvLIBS.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vvLIBS.push_back(vector<string>()); 
  // [OBSOLETE] vLibrary_ICSD.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) vLibrary_ICSD.push_back(""); 
  // [OBSOLETE] vLibrary_ICSD_ALL.clear();
  XHOST_vLibrary_ICSD.clear();for(uint i=0;i<XHOST_vGlobal_MAX;i++) XHOST_vLibrary_ICSD.push_back("");  // needs some initialization
  vflag_aflow.clear();
  vflag_pflow.clear();
  vflag_apennsy.clear();
  vflag_outreach.clear();
  vflag_control.clear();
  vschema.clear();
  vschema_internal.clear();  //ME20220208
  // extensions
  vcat.clear();vcat.push_back("cat");vcat.push_back("bzcat"); vcat.push_back("xzcat");vcat.push_back("gzcat");
  vext.clear();vext.push_back("");   vext.push_back(".bz2");  vext.push_back(".xz");  vext.push_back(".gz");
  vzip.clear();vzip.push_back("");   vzip.push_back("bzip2"); vzip.push_back("xz");   vzip.push_back("gzip");
  // AFLOWRC
  aflowrc_filename.clear();    // AFLOWRC
  aflowrc_content.clear();    // AFLOWRC
  vaflowrc.clear();   // AFLOWRC
  adefault.clear();  // AFLOWRC
}

void _XHOST::clear() {  // clear PRIVATE
  _XHOST XHOST_temp;
  copy(XHOST_temp);
}

// other public functions of XHOST

pthread_mutex_t mutex_XAFLOW_XHOST=PTHREAD_MUTEX_INITIALIZER;

std::string _XHOST::command(const string& command) {
  string _command=command;
#ifdef _MACOSX_
  if(command=="beep") return string("echo -ne '\007'");
#endif
  //first check for EXACT match in vcmds  //CO20180705
  if(command=="aflow_data"&&aurostd::FileExist("./aflow_data")) {return "./aflow_data";} //CO20180705 - hack for developers, prefer ./aflow_data over aflow_data in PATH, as it is probably newer (development)
  for(uint i=0;i<vcmd.size();i++) {
    if(vcmd.at(i)==command) {
      return vcmd.at(i);
    }// found before.. only == otherwise cat gets confused with bzcat
  }
  //next check if we can find the command in PATH or somewher else common //CO20180705
  if(aurostd::IsCommandAvailableModify(_command)) {
    pthread_mutex_lock(&mutex_XAFLOW_XHOST);
    // pthread_mutex_unlock(&mutex_XAFLOW_XHOST);
    vcmd.push_back(_command);
    pthread_mutex_unlock(&mutex_XAFLOW_XHOST);
    return _command; } // found and added
  //CO20180705 START
  //requested command is "aflow_data", and we have ./aflow_data in our vcmd's
  //we need to strip vcmd's to basename and check
  vector<string> path_parts;
  string basename;
  for(uint i=0;i<vcmd.size();i++) {
    aurostd::string2tokens(vcmd[i],path_parts,"/",true);  //consecutive
    if(path_parts.size()) {
      basename=path_parts[path_parts.size()-1];
      if(basename==command) {
        return vcmd[i];
      }
    }
  }
  //CO20180705 STOP
  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"command="+command+" not found",_INPUT_MISSING_);
  return string();
}

bool _XHOST::is_command(const string& command) {
  try{string path=_XHOST::command(command);return !path.empty();}    //CO20190629 - using xerror now
  catch(aurostd::xerror& excpt) {return false;}  //CO20190629 - using xerror now
  //[CO20190629 - using xerror now]return !_XHOST::command(command).empty(); //CO20180705
  //[CO20180705 OBSOLETE]string _command=command;
  //[CO20180705 OBSOLETE]#ifdef _MACOSX_
  //[CO20180705 OBSOLETE]if(command=="beep") return TRUE;
  //[CO20180705 OBSOLETE]#endif
  //[CO20180705 OBSOLETE]for(uint i=0;i<vcmd.size();i++) if(aurostd::substring2bool(vcmd.at(i),_command)) return TRUE; // found before..
  //[CO20180705 OBSOLETE]if(aurostd::IsCommandAvailableModify(_command)) return TRUE;
  //[CO20180705 OBSOLETE]return FALSE; // not found
}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// _AFLAGS
// look into aflow.h for the definitions

// constructors
_aflags::_aflags() {
  AFLOW_PTHREADS_NUMBER                 = 0;
  LocalDirectory                        = ""; // should do PWD
  Directory                             = "";
  QUIET                                 = FALSE;
  AFLOW_FORCE_RUN                       = FALSE;
  AFLOW_PERFORM_DIRECTORY               = FALSE;
  AFLOW_PERFORM_FILE                    = FALSE;
  AFLOW_PERFORM_ORDER_SORT              = FALSE;
  AFLOW_PERFORM_ORDER_REVERSE           = FALSE;
  AFLOW_PERFORM_ORDER_RANDOM            = FALSE;
  AFLOW_MODE_GENERATE                   = FALSE;
  AFLOW_MODE_QSUB_MODE1                 = FALSE;
  AFLOW_MODE_QSUB_MODE2                 = FALSE;
  AFLOW_MODE_QSUB_MODE3                 = FALSE;
  KBIN_RUN_AFLOWIN                      = FALSE;
  KBIN_GEN_GENERAL                      = FALSE;  //CO20180405 - general --generate flag, we need to read aflow.in to determine whether vasp/sym/aims/etc.
  KBIN_GEN_VASP_FROM_AFLOWIN            = FALSE;
  KBIN_GEN_AIMS_FROM_AFLOWIN            = FALSE;
  KBIN_GEN_AFLOWIN_FROM_VASP            = FALSE;
  //DX START
  KBIN_GEN_SYMMETRY_OF_AFLOWIN          = FALSE;
  //DX END
  KBIN_DELETE_AFLOWIN                   = FALSE;
  AFLOW_FORCE_MPI                       = FALSE;
  AFLOW_FORCE_SERIAL                    = FALSE;
  AFLOW_GLOBAL_NCPUS                    = -1;
  // Perform TASKS
  AFLOW_PERFORM_CLEAN                   = FALSE;
  // host related things
  AFLOW_MACHINE_GLOBAL.clear();
  AFLOW_MACHINE_LOCAL.clear();
  // APENNSY Flags
  // [OBSOLETE] APENNSY_VERBOSE_flag                  = FALSE;   // APENNSY flags
  // [OBSOLETE] APENNSY_LIST_flag                     = FALSE;   // APENNSY flags
  // [OBSOLETE] APENNSY_SERVER_AFLOWLIB_ORG_flag      = FALSE;   // APENNSY flags
  APENNSY_LATTICE_flag                  = "ALL";   // APENNSY flags
  APENNSY_GNUPLOT_FONT_str              = DEFAULT_GNUPLOT_EPS_FONT;  // APENNSY content
  APENNSY_GNUPLOT_FONT_BOLD_str         = DEFAULT_GNUPLOT_EPS_FONT_BOLD; // APENNSY content
  APENNSY_GNUPLOT_FONT_ITALICS_str      = DEFAULT_GNUPLOT_EPS_FONT_ITALICS; // APENNSY content
  APENNSY_NEGLECT_STRUCTURES_vstrs.clear();        // APENNSY content
  vflag.clear();
}

// destructor
_aflags::~_aflags() {
  free();
}

void _aflags::free() {
}

void _aflags::copy(const _aflags& b) {
  AFLOW_PTHREADS_NUMBER                 = b.AFLOW_PTHREADS_NUMBER;
  LocalDirectory                        = b.LocalDirectory;
  Directory                             = b.Directory;
  QUIET                                 = b.QUIET;
  AFLOW_FORCE_RUN                       = b.AFLOW_FORCE_RUN;
  AFLOW_PERFORM_DIRECTORY               = b.AFLOW_PERFORM_DIRECTORY;
  AFLOW_PERFORM_FILE                    = b.AFLOW_PERFORM_FILE;
  AFLOW_PERFORM_ORDER_SORT              = b.AFLOW_PERFORM_ORDER_SORT;
  AFLOW_PERFORM_ORDER_REVERSE           = b.AFLOW_PERFORM_ORDER_REVERSE;
  AFLOW_PERFORM_ORDER_RANDOM            = b.AFLOW_PERFORM_ORDER_RANDOM;
  AFLOW_MODE_GENERATE                   = b.AFLOW_MODE_GENERATE;
  AFLOW_MODE_QSUB_MODE1                 = b.AFLOW_MODE_QSUB_MODE1;
  AFLOW_MODE_QSUB_MODE2                 = b.AFLOW_MODE_QSUB_MODE2;
  AFLOW_MODE_QSUB_MODE3                 = b.AFLOW_MODE_QSUB_MODE3;
  KBIN_RUN_AFLOWIN                      = b.KBIN_RUN_AFLOWIN;
  KBIN_GEN_VASP_FROM_AFLOWIN            = b.KBIN_GEN_VASP_FROM_AFLOWIN;
  KBIN_GEN_AIMS_FROM_AFLOWIN            = b.KBIN_GEN_AIMS_FROM_AFLOWIN;
  KBIN_GEN_AFLOWIN_FROM_VASP            = b.KBIN_GEN_AFLOWIN_FROM_VASP;
  //DX START
  KBIN_GEN_SYMMETRY_OF_AFLOWIN          = b.KBIN_GEN_SYMMETRY_OF_AFLOWIN;
  //DX END
  KBIN_DELETE_AFLOWIN                   = b.KBIN_DELETE_AFLOWIN;
  AFLOW_FORCE_MPI                       = b.AFLOW_FORCE_MPI;
  AFLOW_FORCE_SERIAL                    = b.AFLOW_FORCE_SERIAL;
  AFLOW_GLOBAL_NCPUS                    = b.AFLOW_GLOBAL_NCPUS;
  // Perform TASKS
  AFLOW_PERFORM_CLEAN                   = b.AFLOW_PERFORM_CLEAN;
  // host related things
  AFLOW_MACHINE_GLOBAL                  = b.AFLOW_MACHINE_GLOBAL;
  AFLOW_MACHINE_LOCAL                   = b.AFLOW_MACHINE_LOCAL;
  // APENNSY Flags
  // [OBSOLETE] APENNSY_VERBOSE_flag                  = b.APENNSY_VERBOSE_flag;                          // APENNSY flags
  // [OBSOLETE] APENNSY_LIST_flag                     = b.APENNSY_LIST_flag;                            // APENNSY flags
  // [OBSOLETE] APENNSY_SERVER_AFLOWLIB_ORG_flag      = b.APENNSY_SERVER_AFLOWLIB_ORG_flag;                     // APENNSY flags
  APENNSY_LATTICE_flag                  = b.APENNSY_LATTICE_flag;                          // APENNSY flags
  APENNSY_GNUPLOT_FONT_str              = b.APENNSY_GNUPLOT_FONT_str;                      // APENNSY flags
  APENNSY_GNUPLOT_FONT_BOLD_str         = b.APENNSY_GNUPLOT_FONT_BOLD_str;                 // APENNSY flags
  APENNSY_GNUPLOT_FONT_ITALICS_str      = b.APENNSY_GNUPLOT_FONT_ITALICS_str;              // APENNSY flags
  APENNSY_NEGLECT_STRUCTURES_vstrs.clear();                                                // APENNSY content
  for(uint i=0;i<b.APENNSY_NEGLECT_STRUCTURES_vstrs.size();i++)                            // APENNSY content
    APENNSY_NEGLECT_STRUCTURES_vstrs.push_back(b.APENNSY_NEGLECT_STRUCTURES_vstrs.at(i));  // APENNSY content
  vflag.clear();vflag=b.vflag;
}

// copy
_aflags::_aflags(const _aflags& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _aflags& _aflags::operator=(const _aflags& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _aflags::clear() {
  _aflags aflags_temp;
  copy(aflags_temp);
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// _KFLAGS
// look into aflow.h for the definitions

// constructors
_kflags::_kflags() {
  AFLOW_MODE_ALIEN                                 = FALSE;
  AFLOW_MODE_MATLAB                                = FALSE;
  AFLOW_MATLAB_MODE_EXPLICIT                       = FALSE;
  AFLOW_MATLAB_MODE_EXPLICIT_START_STOP            = FALSE;
  AFLOW_MATLAB_MODE_IMPLICIT                       = FALSE;
  AFLOW_MATLAB_MODE_EXTERNAL                       = FALSE;
  AFLOW_MATLAB_FILE                                = FALSE;
  AFLOW_MATLAB_FILE_FILE                           = FALSE;
  AFLOW_MATLAB_FILE_COMMAND                        = FALSE;
  AFLOW_MODE_VASP                                  = FALSE;
  AFLOW_MODE_AIMS                                  = FALSE;
  AFLOW_MODE_PRESCRIPT_EXPLICIT                    = FALSE;
  AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP         = FALSE;
  AFLOW_MODE_PRESCRIPT.str(std::string());
  AFLOW_MODE_POSTSCRIPT_EXPLICIT                   = FALSE;
  AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP        = FALSE;
  AFLOW_MODE_POSTSCRIPT.str(std::string());
  AFLOW_MODE_EMAIL                                 = FALSE;
  KBIN_BIN                                         = "";
  KBIN_SERIAL_BIN                                  = ""; //ME20190107
  KZIP_BIN                                         = "";
  KZIP_COMPRESS                                    = FALSE;
  KBIN_MPI                                         = FALSE;
  KBIN_MPI_NCPUS                                   = 0;
  KBIN_MPI_NCPUS_STRING                            = ""; //ME20181216
  KBIN_MPI_NCPUS_ORIG                              = 0; //CO20210804
  KBIN_MPI_START                                   = "";
  KBIN_MPI_STOP                                    = "";
  KBIN_MPI_COMMAND                                 = "";
  KBIN_MPI_AUTOTUNE                                = FALSE;
  KBIN_MPI_BIN                                     = "";
  KBIN_MPI_OPTIONS                                 = "";
  KBIN_QSUB                                        = FALSE;
  KBIN_QSUB_MODE1                                  = FALSE;
  KBIN_QSUB_MODE2                                  = FALSE;
  KBIN_QSUB_MODE3                                  = FALSE;
  KBIN_QSUB_COMMAND                                = "";
  KBIN_QSUB_PARAMS                                 = "";
  KBIN_QSUB_MODE_EXPLICIT                          = FALSE;
  KBIN_QSUB_MODE_EXPLICIT_START_STOP               = FALSE;
  KBIN_QSUB_MODE_IMPLICIT                          = FALSE;
  KBIN_QSUB_FILE                                   = FALSE;
  KBIN_SYMMETRY_CALCULATION                        = FALSE;
  //DX START
  KBIN_SYMMETRY_NO_SCAN                            = FALSE;
  KBIN_SYMMETRY_EPS                                = AUROSTD_NAN;
  KBIN_SYMMETRY_CALCULATE_PGROUP                   = TRUE; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUPK                  = TRUE; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_FGROUP                   = TRUE; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL              = TRUE; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL             = TRUE; //DX20200423 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON        = TRUE; //DX20200423 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_IATOMS                   = TRUE; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_AGROUP                   = TRUE; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_SGROUP                   = TRUE; //DX20170814 - Specify what to calculate/verify
  //DX END
  KBIN_SYMMETRY_PGROUP_WRITE                       = FALSE;
  KBIN_SYMMETRY_PGROUP_XTAL_WRITE                  = FALSE;
  KBIN_SYMMETRY_PGROUPK_WRITE                      = FALSE;
  KBIN_SYMMETRY_PGROUPK_XTAL_WRITE                 = FALSE; //DX20200423
  KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE            = FALSE; //DX20200423
  KBIN_SYMMETRY_FGROUP_WRITE                       = FALSE;
  KBIN_SYMMETRY_SGROUP_WRITE                       = FALSE;
  KBIN_SYMMETRY_AGROUP_WRITE                       = FALSE;
  KBIN_SYMMETRY_IATOMS_WRITE                       = FALSE;
  KBIN_SYMMETRY_SGROUP_RADIUS                      = 0.0;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_CALCULATION                       = FALSE;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_WRITE                             = FALSE;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_RADIUS                            = 0.0;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_DRADIUS                           = 0.0;
  KBIN_POCC                                        = FALSE;
  KBIN_POCC_CALCULATION                            = FALSE;
  KBIN_POCC_TEMPERATURE_STRING                     = "";  //CO20191110
  KBIN_POCC_ARUNS2SKIP_STRING                      = "";  //CO20200624
  KBIN_POCC_EXCLUDE_UNSTABLE                       = FALSE;  //ME20210927
  KBIN_FROZSL                                      = FALSE;
  KBIN_FROZSL_DOWNLOAD                             = FALSE;
  KBIN_FROZSL_FILE                                 = FALSE;
  KBIN_FROZSL_FILE_NAME                            = DEFAULT_AFLOW_FROZSL_INPUT_OUT;
  // [OBSOLETE]  KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT              = FALSE;
  // [OBSOLETE]  KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT_START_STOP   = FALSE;
  // [OBSOLETE]  KBIN_FROZSL_PRESCRIPT_STRING                     = "";
  // [OBSOLETE]  KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT             = FALSE;
  // [OBSOLETE]  KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT_START_STOP  = FALSE;
  // [OBSOLETE]  KBIN_FROZSL_POSTSCRIPT_STRING                    = "";
  KBIN_FROZSL_STRUCTURE_MODE_FILE                  = FALSE;
  KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP   = FALSE;
  KBIN_FROZSL_STRUCTURE_STRING                     = "";
  KBIN_FROZSL_DIELECTRIC_MODE_FILE                 = FALSE;
  KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP  = FALSE;
  KBIN_FROZSL_DIELECTRIC_ZEFF                      = FALSE;
  KBIN_FROZSL_DIELECTRIC_STRING                    = "";
  KBIN_PHONONS_CALCULATION_APL                     = FALSE;
  KBIN_PHONONS_CALCULATION_QHA                     = FALSE; //CO20170601
  KBIN_PHONONS_CALCULATION_AAPL                    = FALSE; //CO20170601
  KBIN_PHONONS_CALCULATION_AGL                     = FALSE;
  KBIN_PHONONS_CALCULATION_AEL                     = FALSE;
  KBIN_PHONONS_CALCULATION_FROZSL                  = FALSE;
  KBIN_PHONONS_CALCULATION_FROZSL_output           = "";
  KBIN_PHONONS_CALCULATION_FROZSL_poscars          = "";
}

// destructor
_kflags::~_kflags() {
  free();
}

void _kflags::free() {
}

void _kflags::copy(const _kflags& b) {
  AFLOW_MODE_ALIEN                                 = b.AFLOW_MODE_ALIEN;
  AFLOW_MODE_MATLAB                                = b.AFLOW_MODE_MATLAB;
  AFLOW_MATLAB_MODE_EXPLICIT                       = b.AFLOW_MATLAB_MODE_EXPLICIT;
  AFLOW_MATLAB_MODE_EXPLICIT_START_STOP            = b.AFLOW_MATLAB_MODE_EXPLICIT_START_STOP;
  AFLOW_MATLAB_MODE_IMPLICIT                       = b.AFLOW_MATLAB_MODE_IMPLICIT;
  AFLOW_MATLAB_MODE_EXTERNAL                       = b.AFLOW_MATLAB_MODE_EXTERNAL;
  AFLOW_MATLAB_FILE                                = b.AFLOW_MATLAB_FILE;
  AFLOW_MATLAB_FILE_FILE                           = b.AFLOW_MATLAB_FILE_FILE;
  AFLOW_MATLAB_FILE_COMMAND                        = b.AFLOW_MATLAB_FILE_COMMAND;
  AFLOW_MODE_VASP                                  = b.AFLOW_MODE_VASP;
  AFLOW_MODE_AIMS                                  = b.AFLOW_MODE_AIMS;
  AFLOW_MODE_PRESCRIPT_EXPLICIT                    = b.AFLOW_MODE_PRESCRIPT_EXPLICIT;
  AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP         = b.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP;
  AFLOW_MODE_PRESCRIPT.str(std::string()); AFLOW_MODE_PRESCRIPT << b.AFLOW_MODE_PRESCRIPT.str();
  AFLOW_MODE_POSTSCRIPT_EXPLICIT                   = b.AFLOW_MODE_POSTSCRIPT_EXPLICIT;
  AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP        = b.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP;
  AFLOW_MODE_POSTSCRIPT.str(std::string()); AFLOW_MODE_POSTSCRIPT << b.AFLOW_MODE_POSTSCRIPT.str();
  AFLOW_MODE_EMAIL                                 = b.AFLOW_MODE_EMAIL;
  KBIN_BIN                                         = b.KBIN_BIN;
  KBIN_SERIAL_BIN                                  = b.KBIN_SERIAL_BIN;  //ME20190107
  KZIP_BIN                                         = b.KZIP_BIN;
  KZIP_COMPRESS                                    = b.KZIP_COMPRESS;
  KBIN_MPI                                         = b.KBIN_MPI;
  KBIN_MPI_NCPUS                                   = b.KBIN_MPI_NCPUS;
  KBIN_MPI_NCPUS_STRING                            = b.KBIN_MPI_NCPUS_STRING;	//ME20181216
  KBIN_MPI_NCPUS_ORIG                              = b.KBIN_MPI_NCPUS_ORIG; //CO20210804
  KBIN_MPI_START                                   = b.KBIN_MPI_START;
  KBIN_MPI_STOP                                    = b.KBIN_MPI_STOP;
  KBIN_MPI_COMMAND                                 = b.KBIN_MPI_COMMAND;
  KBIN_MPI_AUTOTUNE                                = b.KBIN_MPI_AUTOTUNE;
  KBIN_MPI_BIN                                     = b.KBIN_MPI_BIN;
  KBIN_MPI_OPTIONS                                 = b.KBIN_MPI_OPTIONS;
  KBIN_QSUB                                        = b.KBIN_QSUB;
  KBIN_QSUB_MODE1                                  = b.KBIN_QSUB_MODE1;
  KBIN_QSUB_MODE2                                  = b.KBIN_QSUB_MODE2;
  KBIN_QSUB_MODE3                                  = b.KBIN_QSUB_MODE3;
  KBIN_QSUB_COMMAND                                = b.KBIN_QSUB_COMMAND;
  KBIN_QSUB_PARAMS                                 = b.KBIN_QSUB_PARAMS;
  KBIN_QSUB_MODE_EXPLICIT                          = b.KBIN_QSUB_MODE_EXPLICIT;
  KBIN_QSUB_MODE_EXPLICIT_START_STOP               = b.KBIN_QSUB_MODE_EXPLICIT_START_STOP;
  KBIN_QSUB_MODE_IMPLICIT                          = b.KBIN_QSUB_MODE_IMPLICIT;
  KBIN_QSUB_FILE                                   = b.KBIN_QSUB_FILE;
  KBIN_SYMMETRY_CALCULATION                        = b.KBIN_SYMMETRY_CALCULATION;
  //DX START
  KBIN_SYMMETRY_NO_SCAN                            = b.KBIN_SYMMETRY_NO_SCAN;
  KBIN_SYMMETRY_EPS                                = b.KBIN_SYMMETRY_EPS;
  KBIN_SYMMETRY_CALCULATE_PGROUP                   = b.KBIN_SYMMETRY_CALCULATE_PGROUP; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUPK                  = b.KBIN_SYMMETRY_CALCULATE_PGROUPK; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_FGROUP                   = b.KBIN_SYMMETRY_CALCULATE_FGROUP; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL              = b.KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL             = b.KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL; //DX20200423 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON        = b.KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON; //DX20200423 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_IATOMS                   = b.KBIN_SYMMETRY_CALCULATE_IATOMS; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_AGROUP                   = b.KBIN_SYMMETRY_CALCULATE_AGROUP; //DX20170814 - Specify what to calculate/verify
  KBIN_SYMMETRY_CALCULATE_SGROUP                   = b.KBIN_SYMMETRY_CALCULATE_SGROUP; //DX20170814 - Specify what to calculate/verify
  //DX END
  KBIN_SYMMETRY_PGROUP_WRITE                       = b.KBIN_SYMMETRY_PGROUP_WRITE;
  KBIN_SYMMETRY_PGROUP_XTAL_WRITE                  = b.KBIN_SYMMETRY_PGROUP_XTAL_WRITE;
  KBIN_SYMMETRY_PGROUPK_WRITE                      = b.KBIN_SYMMETRY_PGROUPK_WRITE;
  KBIN_SYMMETRY_PGROUPK_XTAL_WRITE                 = b.KBIN_SYMMETRY_PGROUPK_XTAL_WRITE; //DX20200423
  KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE            = b.KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE; //DX20200423
  KBIN_SYMMETRY_FGROUP_WRITE                       = b.KBIN_SYMMETRY_FGROUP_WRITE;
  KBIN_SYMMETRY_SGROUP_WRITE                       = b.KBIN_SYMMETRY_SGROUP_WRITE;
  KBIN_SYMMETRY_AGROUP_WRITE                       = b.KBIN_SYMMETRY_AGROUP_WRITE;
  KBIN_SYMMETRY_IATOMS_WRITE                       = b.KBIN_SYMMETRY_IATOMS_WRITE;
  KBIN_SYMMETRY_SGROUP_RADIUS                      = b.KBIN_SYMMETRY_SGROUP_RADIUS;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_CALCULATION                       = b.KBIN_NEIGHBORS_CALCULATION;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_WRITE                             = b.KBIN_NEIGHBORS_WRITE;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_RADIUS                            = b.KBIN_NEIGHBORS_RADIUS;
  //DX20210122 [OBSOLETE] KBIN_NEIGHBORS_DRADIUS                           = b.KBIN_NEIGHBORS_DRADIUS;
  KBIN_POCC                                        = b.KBIN_POCC;
  KBIN_POCC_CALCULATION                            = b.KBIN_POCC_CALCULATION; //CO20191110
  KBIN_POCC_TEMPERATURE_STRING                     = b.KBIN_POCC_TEMPERATURE_STRING; //CO20191110
  KBIN_POCC_ARUNS2SKIP_STRING                      = b.KBIN_POCC_ARUNS2SKIP_STRING; //CO20200627
  KBIN_POCC_EXCLUDE_UNSTABLE                       = b.KBIN_POCC_EXCLUDE_UNSTABLE;  //ME20210927
  KBIN_FROZSL                                      = b.KBIN_FROZSL;
  KBIN_FROZSL_DOWNLOAD                             = b.KBIN_FROZSL_DOWNLOAD;
  KBIN_FROZSL_FILE                                 = b.KBIN_FROZSL_FILE;
  KBIN_FROZSL_FILE_NAME                            = b.KBIN_FROZSL_FILE_NAME;
  // [OBSOLETE]  KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT              = b.KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT;
  // [OBSOLETE]  KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT_START_STOP   = b.KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT_START_STOP;
  // [OBSOLETE]  KBIN_FROZSL_PRESCRIPT_STRING                     = b.KBIN_FROZSL_PRESCRIPT_STRING;
  // [OBSOLETE]  KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT             = b.KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT;
  // [OBSOLETE]  KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT_START_STOP  = b.KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT_START_STOP;
  // [OBSOLETE]  KBIN_FROZSL_POSTSCRIPT_STRING                    = b.KBIN_FROZSL_POSTSCRIPT_STRING;
  KBIN_FROZSL_STRUCTURE_MODE_FILE                  = b.KBIN_FROZSL_STRUCTURE_MODE_FILE;
  KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP   = b.KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP;
  KBIN_FROZSL_STRUCTURE_STRING                     = b.KBIN_FROZSL_STRUCTURE_STRING;
  KBIN_FROZSL_DIELECTRIC_MODE_FILE                 = b.KBIN_FROZSL_DIELECTRIC_MODE_FILE;
  KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP  = b.KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP;
  KBIN_FROZSL_DIELECTRIC_STRING                    = b.KBIN_FROZSL_DIELECTRIC_STRING;
  KBIN_PHONONS_CALCULATION_APL                     = b.KBIN_PHONONS_CALCULATION_APL;
  KBIN_PHONONS_CALCULATION_QHA                     = b.KBIN_PHONONS_CALCULATION_QHA;  //CO20170601
  KBIN_PHONONS_CALCULATION_AAPL                    = b.KBIN_PHONONS_CALCULATION_AAPL; //CO20170601
  KBIN_PHONONS_CALCULATION_AGL                     = b.KBIN_PHONONS_CALCULATION_AGL;
  KBIN_PHONONS_CALCULATION_AEL                     = b.KBIN_PHONONS_CALCULATION_AEL;
  KBIN_PHONONS_CALCULATION_FROZSL                  = b.KBIN_PHONONS_CALCULATION_FROZSL;
  KBIN_PHONONS_CALCULATION_FROZSL_output           = b.KBIN_PHONONS_CALCULATION_FROZSL_output;
  KBIN_PHONONS_CALCULATION_FROZSL_poscars          = b.KBIN_PHONONS_CALCULATION_FROZSL_poscars;
  KBIN_MODULE_OPTIONS                              = b.KBIN_MODULE_OPTIONS;  //ME20181027
}

// copy
_kflags::_kflags(const _kflags& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _kflags& _kflags::operator=(const _kflags& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _kflags::clear() {
  _kflags kflags_temp;
  copy(kflags_temp);
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// _VFLAGS
// look into aflow.h for the definitions

// constructors
_vflags::_vflags() {free();}

// destructor
_vflags::~_vflags() {free();}

void _vflags::free() {
  // SYSTEM
  AFLOW_SYSTEM.clear(); //ME20181121

  KBIN_VASP_RUN_NRELAX                                           = 0;
  //  RUN
  KBIN_VASP_RUN.clear();                                         // RUN
  KBIN_VASP_RUN.push("");                                   // RUN
  //  REPEAT
  KBIN_VASP_REPEAT.clear();                                      // REPEAT
  KBIN_VASP_REPEAT.push("");                                // REPEAT

  KBIN_VASP_FORCE_OPTION_NOTUNE.clear();                         // NOTUNE
  KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO.clear();                    // SYSTEM_AUTO

  //  RELAX_MODE
  KBIN_VASP_FORCE_OPTION_RELAX_MODE.clear();                     // RELAX_MODE forces/energy     
  KBIN_VASP_FORCE_OPTION_RELAX_MODE.push(DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME); // RELAX_MODE forces/energy    

  //  RELAX_TYPE
  KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();                     // RELAX_TYPE geometry   
  KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("");               // RELAX_TYPE geometry    

  KBIN_VASP_FORCE_OPTION_PREC.clear();                           // PREC
  KBIN_VASP_FORCE_OPTION_PREC.push(DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME); // PREC
  KBIN_VASP_FORCE_OPTION_ALGO.clear();                           // ALGO
  KBIN_VASP_FORCE_OPTION_ALGO.push(DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME); // ALGO
  KBIN_VASP_FORCE_OPTION_METAGGA.clear();                        // METAGGA
  KBIN_VASP_FORCE_OPTION_METAGGA.push(DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME); // METAGGA
  KBIN_VASP_FORCE_OPTION_IVDW.clear();                           // IVDW
  KBIN_VASP_FORCE_OPTION_IVDW.push(DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME); // IVDW
  KBIN_VASP_FORCE_OPTION_ABMIX.clear();                          // ABMIX
  KBIN_VASP_FORCE_OPTION_ABMIX.push(DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME); // ABMIX
  KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.clear();          // AUTO_PSEUDOPOTENTIALS
  KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.push(DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE); // AUTO_PSEUDOPOTENTIALS

  // ENMAX_MULTIPLY
  KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.clear();
  //  KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.push("1.4");

  // NBANDS
  KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry                     = FALSE;
  // NBANDS
  KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.clear();
  //  KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.push("0");
  // POTIM
  KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.clear();
  KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_PREC_POTIM));
  // PSTRESS
  KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.clear();
  // KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.push("0.0");
  // EDIFFG
  KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.clear();
  KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG));
  // NELM  //CO20200624
  KBIN_VASP_FORCE_OPTION_NELM_EQUAL.clear(); //CO20200624
  KBIN_VASP_FORCE_OPTION_NELM_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_NELM));  //default //CO20200624
  // NELM_STATIC  //CO20200624
  KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.clear(); //CO20200624
  KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_NELM_STATIC));  //default //CO20200624
  // ISMEAR  //CO20181128
  KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.clear(); //CO20181128
  KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME));  //default //CO20181128
  // SIGMA  //CO20181128
  KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.clear(); //CO20181128
  KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_SIGMA));  //default //CO20181128
  // ISMEAR_STATIC  //CO20181128
  KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.clear(); //CO20181128
  KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME));  //default //CO20181128
  // SIGMA_STATIC  //CO20181128
  KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.clear(); //CO20181128
  KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC));  //default //CO20181128
  // ISMEAR_BANDS  //CO20181128
  KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.clear(); //CO20181128
  KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME));  //default //CO20181128
  // SIGMA_BANDS  //CO20181128
  KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.clear(); //CO20181128
  KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.push(aurostd::utype2string(DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS));  //default //CO20181128
  // RWIGS
  KBIN_VASP_FORCE_OPTION_RWIGS_STATIC                            = FALSE;

  KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.clear();                     // SKIP_NOMIX

  //  KBIN_VASP_FORCE_OPTION_TRISTATE.clear();                       // 

  KBIN_VASP_FORCE_OPTION_SPIN.clear();                           // SPIN
  KBIN_VASP_FORCE_OPTION_SPIN.option                             = DEFAULT_VASP_FORCE_OPTION_SPIN;  // SPIN 
  KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1                     = DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;  // SPIN_REMOVE_RELAX_1
  KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2                     = DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2;  // SPIN_REMOVE_RELAX_2
  KBIN_VASP_FORCE_OPTION_BADER.clear();                          // BADER
  KBIN_VASP_FORCE_OPTION_BADER.option                            = DEFAULT_VASP_FORCE_OPTION_BADER;  // BADER
  KBIN_VASP_FORCE_OPTION_ELF.clear();                            // ELF
  KBIN_VASP_FORCE_OPTION_ELF.option                              = DEFAULT_VASP_FORCE_OPTION_ELF;  // ELF
  KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.clear();                    // AUTO_MAGMOM
  KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option                      = DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM; // AUTO_MAGMOM
  KBIN_VASP_FORCE_OPTION_SYM.clear();                            // SYM
  KBIN_VASP_FORCE_OPTION_SYM.option                              = DEFAULT_VASP_FORCE_OPTION_SYM; // SYM
  KBIN_VASP_FORCE_OPTION_WAVECAR.clear();                        // WAVECAR
  KBIN_VASP_FORCE_OPTION_WAVECAR.option                          = DEFAULT_VASP_FORCE_OPTION_WAVECAR; // WAVECAR
  KBIN_VASP_FORCE_OPTION_CHGCAR.clear();                         // CHGCAR
  KBIN_VASP_FORCE_OPTION_CHGCAR.option                           = DEFAULT_VASP_FORCE_OPTION_CHGCAR; // CHGCAR
  KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.clear();                    //ME20191028
  KBIN_VASP_FORCE_OPTION_LSCOUPLING.clear();                     // LSCOUPLING
  KBIN_VASP_FORCE_OPTION_LSCOUPLING.option                       = DEFAULT_VASP_FORCE_OPTION_LSCOUPLING; // LSCOUPLING

  KBIN_VASP_FORCE_OPTION_LDAU0.clear();                          // LDAU0
  KBIN_VASP_FORCE_OPTION_LDAU1.clear();                          // LDAU1
  KBIN_VASP_FORCE_OPTION_LDAU2.clear();                          // LDAU2
  KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.clear();                 // LDAU_ADIABATIC
  KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int=0;           // LDAU_ADIABATIC
  KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.clear();                    // LDAU_CUTOFF
  KBIN_VASP_LDAU_SPECIES                                         = "";
  KBIN_VASP_LDAU_PARAMETERS                                      = "";
  KBIN_VASP_LDAU_AFLOW_AUTO_flag                                 = TRUE;
  // FORCE_OPTIONS
  KBIN_VASP_FORCE_OPTION_TYPE.clear();                           // TYPE
  KBIN_VASP_FORCE_OPTION_TYPE.push(DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME); // TYPE

  KBIN_VASP_FORCE_OPTION_NSW_EQUAL                               = FALSE;
  KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE                         = 0;

  // AFIX
  KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.clear();                    // AFIX
  // CONVERT
  KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.clear(); 
  //  KBIN_VASP_FORCE_OPTION_VOLUME
  KBIN_VASP_FORCE_OPTION_VOLUME.clear();                         // RUN
  KBIN_VASP_FORCE_OPTION_VOLUME.push("");                   // RUN

  KBIN_VASP_FORCE_OPTION_KPOINTS.clear();                        // KPOINTS

  // KBIN_VASP_INCAR_MODE
  KBIN_VASP_INCAR_MODE.clear();                                  // all false
  // RELAX
  // KBIN_VASP_KPOINTS_MODE
  KBIN_VASP_KPOINTS_MODE.clear();                                // all false

  KBIN_VASP_KPOINTS_KMODE.clear();                               // KMODE
  KBIN_VASP_KPOINTS_KMODE.push("0");                        // KMODE
  KBIN_VASP_KPOINTS_KPPRA.clear();                               // KPPRA
  KBIN_VASP_KPOINTS_KPPRA.push("100");                      // KPPRA
  KBIN_VASP_KPOINTS_KSCHEME.clear();                             // KSCHEME
  KBIN_VASP_KPOINTS_KSCHEME.push(DEFAULT_KSCHEME);//"Monkhorst-Pack"; // KSCHEME
  KBIN_VASP_KPOINTS_KSHIFT.clear();                              // KSCHEME
  KBIN_VASP_KPOINTS_KSHIFT.push("0 0 0");                   // DEFAULT FIX
  // STATIC
  KBIN_VASP_KPOINTS_STATIC_KMODE.clear();                        // KMODE
  KBIN_VASP_KPOINTS_STATIC_KMODE.push("0");                 // KMODE
  KBIN_VASP_KPOINTS_STATIC_KPPRA.clear();                        // KPPRA
  KBIN_VASP_KPOINTS_STATIC_KPPRA.push("1");                 // KPPRA
  KBIN_VASP_KPOINTS_STATIC_KSCHEME.clear();                      // KSCHEME
  KBIN_VASP_KPOINTS_STATIC_KSCHEME.push(DEFAULT_STATIC_KSCHEME);//"Monkhorst-Pack"); // KSCHEME
  KBIN_VASP_KPOINTS_STATIC_KSHIFT.clear();                       // KSCHEME
  KBIN_VASP_KPOINTS_STATIC_KSHIFT.push("0 0 0"); // DEFAULT FIX
  // PHONONS
  KBIN_VASP_KPOINTS_PHONONS_KPPRA.clear();                       // KPPRA
  KBIN_VASP_KPOINTS_PHONONS_KPPRA.push("1");                // KPPRA
  KBIN_VASP_KPOINTS_PHONONS_KSCHEME.clear();                     // KSCHEME
  KBIN_VASP_KPOINTS_PHONONS_KSCHEME.push(DEFAULT_PHONONS_KSCHEME);  //"Gamma");           // KSCHEME

  KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.clear();         // PARITY
  KBIN_VASP_KPOINTS_PHONONS_GRID.clear();                        //ME20200427
  KBIN_VASP_KPOINTS_PHONONS_SHIFT.clear();                       //ME20200427
  // BANDS
  // KBIN_VASP_KPOINTS_BANDS_LATTICE
  KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();                       // LATTICE
  KBIN_VASP_KPOINTS_BANDS_LATTICE.push(DEFAULT_BANDS_LATTICE);             // some default
  // [OBSOLETE] KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG                           = FALSE;
  // [OBSOLETE] KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE                          = "";

  //[CO20210805 - OBSOLETE]KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG                      = FALSE;
  KBIN_VASP_KPOINTS_BANDS_GRID.clear();                         //CO20210805
  KBIN_VASP_KPOINTS_BANDS_GRID.push(aurostd::utype2string(DEFAULT_BANDS_GRID));        //CO20210805
  //[CO20210805 - OBSOLETE]KBIN_VASP_KPOINTS_BANDS_GRID_FLAG                              = FALSE;
  //[CO20210805 - OBSOLETE]KBIN_VASP_KPOINTS_BANDS_GRID_VALUE                             = 0;
  KBIN_VASP_WRITE_KPOINTS                                        = FALSE;
  // KBIN_VASP_POSCAR_MODE
  KBIN_VASP_POSCAR_MODE.clear();                                 // all false
  KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.clear();
  KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.clear();
  // KBIN_VASP_INCAR_FILE
  KBIN_VASP_INCAR_VERBOSE                                        = TRUE; // VERBOSITY IS GOOD !!!
  KBIN_VASP_INCAR_FILE.clear();                                  // all false
  KBIN_VASP_INCAR_EXPLICIT.clear(); //ME20181226
  KBIN_VASP_INCAR_EXPLICIT_START_STOP.str(""); //ME20181226 //CO20190401 - clear ss with .str("")
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_KEYWORD                     = FALSE;
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_SYSTEM_AUTO                 = FALSE;
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_FILE                        = FALSE;
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_COMMAND                     = FALSE;

  // KBIN_VASP_KPOINTS_FILE
  KBIN_VASP_KPOINTS_FILE.clear();                                // all false
  KBIN_VASP_KPOINTS_EXPLICIT.clear(); //ME20181226
  KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(""); //ME20181226 //CO20190401 - clear ss with .str("")
  // [OBSOLETE] KBIN_VASP_KPOINTS_FILE_KEYWORD                   = FALSE;
  // [OBSOLETE] KBIN_VASP_KPOINTS_FILE_FILE                      = FALSE;
  // [OBSOLETE] KBIN_VASP_KPOINTS_FILE_COMMAND                   = FALSE;

  // KBIN_VASP_POSCAR_FILE
  KBIN_VASP_POSCAR_FILE.clear();                                 // all false
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_KEYWORD                    = FALSE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_PROTOTYPE                  = FALSE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_FILE                       = FALSE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_COMMAND                    = FALSE;

  //  KBIN_VASP_POSCAR_FILE_VOLUME
  KBIN_VASP_POSCAR_FILE_VOLUME.clear();                          // RUN
  KBIN_VASP_POSCAR_FILE_VOLUME.push("");                    // RUN

  // KBIN_VASP_POTCAR_MODE  
  KBIN_VASP_POTCAR_MODE.clear();                                 // all false
  // KBIN_VASP_POTCAR_FILE
  KBIN_VASP_POTCAR_FILE.clear();                                 // all false
  KBIN_VASP_POTCAR_EXPLICIT.clear(); //ME20181226
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_KEYWORD                    = FALSE;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_SYSTEM_AUTO                = FALSE;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_PREFIX                     = FALSE;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_SUFFIX                     = FALSE;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_FILE                       = FALSE;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_COMMAND                    = FALSE;
}

void _vflags::copy(const _vflags& b) {
  AFLOW_SYSTEM                                                   = b.AFLOW_SYSTEM; //ME20181226
  KBIN_VASP_RUN_NRELAX                                           = b.KBIN_VASP_RUN_NRELAX;
  // RUN
  KBIN_VASP_RUN                                                  = b.KBIN_VASP_RUN;  
  // REPEAT
  KBIN_VASP_REPEAT                                               = b.KBIN_VASP_REPEAT;

  KBIN_VASP_FORCE_OPTION_NOTUNE                                  = b.KBIN_VASP_FORCE_OPTION_NOTUNE;  // NOTUNE
  KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO                             = b.KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO;  // SYSTEM_AUTO

  // RELAX_MODE
  KBIN_VASP_FORCE_OPTION_RELAX_MODE                              = b.KBIN_VASP_FORCE_OPTION_RELAX_MODE; // RELAX_MODE

  // RELAX_TYPE
  KBIN_VASP_FORCE_OPTION_RELAX_TYPE                              = b.KBIN_VASP_FORCE_OPTION_RELAX_TYPE; // RELAX_MODE

  KBIN_VASP_FORCE_OPTION_PREC                                    = b.KBIN_VASP_FORCE_OPTION_PREC; // PREC
  KBIN_VASP_FORCE_OPTION_ALGO                                    = b.KBIN_VASP_FORCE_OPTION_ALGO; // ALGO
  KBIN_VASP_FORCE_OPTION_METAGGA                                 = b.KBIN_VASP_FORCE_OPTION_METAGGA; // METAGGA
  KBIN_VASP_FORCE_OPTION_IVDW                                    = b.KBIN_VASP_FORCE_OPTION_IVDW; // IVDW
  KBIN_VASP_FORCE_OPTION_ABMIX                                   = b.KBIN_VASP_FORCE_OPTION_ABMIX; // ABMIX
  KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS                   = b.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS; // AUTO_PSEUDOPOTENTIALS

  KBIN_VASP_FORCE_OPTION_SKIP_NOMIX                              = b.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX; // SKIP_NOMIX
  // ENMAX_MULTIPLY
  KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL                    = b.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL;
  // NBANDS
  KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry                     = b.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry;
  KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL                            = b.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL;
  // POTIM
  KBIN_VASP_FORCE_OPTION_POTIM_EQUAL                             = b.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL;
  // PSTRESS
  KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL                           = b.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL;
  // EDIFFG
  KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL                            = b.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL;
  // NELM
  KBIN_VASP_FORCE_OPTION_NELM_EQUAL                              = b.KBIN_VASP_FORCE_OPTION_NELM_EQUAL; //CO20200624
  // NELM_STATIC
  KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL                       = b.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL; //CO20200624
  // ISMEAR
  KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL                            = b.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL; //CO20181128
  // SIGMA
  KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL                             = b.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL; //CO20181128
  // ISMEAR_STATIC
  KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL                     = b.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL; //CO20181128
  // SIGMA_STATIC
  KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL                      = b.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL; //CO20181128
  // ISMEAR_BANDS
  KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL                      = b.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL; //CO20181128
  // SIGMA_BANDS
  KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL                       = b.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL; //CO20181128
  // RWIGS
  KBIN_VASP_FORCE_OPTION_RWIGS_STATIC                            = b.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC;

  // KBIN_VASP_FORCE_OPTION_TRISTATE                                = b.KBIN_VASP_FORCE_OPTION_TRISTATE; // TRISTATE

  KBIN_VASP_FORCE_OPTION_SPIN                                    = b.KBIN_VASP_FORCE_OPTION_SPIN; // SPIN
  KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1                     = b.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;
  KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2                     = b.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2; 
  KBIN_VASP_FORCE_OPTION_BADER                                   = b.KBIN_VASP_FORCE_OPTION_BADER; // BADER
  KBIN_VASP_FORCE_OPTION_ELF                                     = b.KBIN_VASP_FORCE_OPTION_ELF; // ELF
  KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM                             = b.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM; // AUTO_MAGMOM
  KBIN_VASP_FORCE_OPTION_SYM                                     = b.KBIN_VASP_FORCE_OPTION_SYM; // SYM
  KBIN_VASP_FORCE_OPTION_WAVECAR                                 = b.KBIN_VASP_FORCE_OPTION_WAVECAR; // WAVECAR
  KBIN_VASP_FORCE_OPTION_CHGCAR                                  = b.KBIN_VASP_FORCE_OPTION_CHGCAR; // CHGCAR
  KBIN_VASP_FORCE_OPTION_CHGCAR_FILE                             = b.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE;  //ME20191028
  KBIN_VASP_FORCE_OPTION_LSCOUPLING                              = b.KBIN_VASP_FORCE_OPTION_LSCOUPLING; // LSCOUPLING
  KBIN_VASP_FORCE_OPTION_LDAU0                                   = b.KBIN_VASP_FORCE_OPTION_LDAU0;  // LDAU0
  KBIN_VASP_FORCE_OPTION_LDAU1                                   = b.KBIN_VASP_FORCE_OPTION_LDAU1;  // LDAU1
  KBIN_VASP_FORCE_OPTION_LDAU2                                   = b.KBIN_VASP_FORCE_OPTION_LDAU2;  // LDAU2
  KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC                          = b.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC; // LDAU_ADIABATIC
  KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF                             = b.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF; // LDAU_CUTOFF

  KBIN_VASP_LDAU_SPECIES                                         = b.KBIN_VASP_LDAU_SPECIES;
  KBIN_VASP_LDAU_PARAMETERS                                      = b.KBIN_VASP_LDAU_PARAMETERS;
  KBIN_VASP_LDAU_AFLOW_AUTO_flag                                 = b.KBIN_VASP_LDAU_AFLOW_AUTO_flag;
  // FORCE_OPTION
  KBIN_VASP_FORCE_OPTION_TYPE                                    = b.KBIN_VASP_FORCE_OPTION_TYPE; // TYPE

  KBIN_VASP_FORCE_OPTION_NSW_EQUAL                               = b.KBIN_VASP_FORCE_OPTION_NSW_EQUAL;
  KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE                         = b.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE;

  KBIN_VASP_FORCE_OPTION_KPOINTS                                 = b.KBIN_VASP_FORCE_OPTION_KPOINTS;  // KPOINTS
  KBIN_VASP_FORCE_OPTION_IGNORE_AFIX                             = b.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX; // IGNORE_AFIX
  KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL                       = b.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL; // CONVERT_UNIT_CELL

  //  KBIN_VASP_FORCE_OPTION_VOLUME
  KBIN_VASP_FORCE_OPTION_VOLUME                                  = b.KBIN_VASP_FORCE_OPTION_VOLUME;
  // KBIN_VASP_INCAR_MODE
  KBIN_VASP_INCAR_MODE                                           = b.KBIN_VASP_INCAR_MODE;
  // KBIN_VASP_KPOINTS_MODE
  KBIN_VASP_KPOINTS_MODE                                         = b.KBIN_VASP_KPOINTS_MODE;

  // [OBSOLETE] KBIN_VASP_INCAR_MODE_EXPLICIT                    = b.KBIN_VASP_INCAR_MODE_EXPLICIT;
  // [OBSOLETE] KBIN_VASP_INCAR_MODE_EXPLICIT_START_STOP         = b.KBIN_VASP_INCAR_MODE_EXPLICIT_START_STOP;
  // [OBSOLETE] KBIN_VASP_INCAR_MODE_IMPLICIT                    = b.KBIN_VASP_INCAR_MODE_IMPLICIT;
  // [OBSOLETE] KBIN_VASP_INCAR_MODE_EXTERNAL                    = b.KBIN_VASP_INCAR_MODE_EXTERNAL;
  // RELAX
  // [OBSOLETE] KBIN_VASP_KPOINTS_MODE_EXPLICIT                  = b.KBIN_VASP_KPOINTS_MODE_EXPLICIT;
  // [OBSOLETE] KBIN_VASP_KPOINTS_MODE_EXPLICIT_START_STOP       = b.KBIN_VASP_KPOINTS_MODE_EXPLICIT_START_STOP;
  // [OBSOLETE] KBIN_VASP_KPOINTS_MODE_IMPLICIT                  = b.KBIN_VASP_KPOINTS_MODE_IMPLICIT;
  // [OBSOLETE] KBIN_VASP_KPOINTS_MODE_EXTERNAL                  = b.KBIN_VASP_KPOINTS_MODE_EXTERNAL;
  // RELAX
  KBIN_VASP_KPOINTS_KMODE                                        = b.KBIN_VASP_KPOINTS_KMODE;
  KBIN_VASP_KPOINTS_KPPRA                                        = b.KBIN_VASP_KPOINTS_KPPRA;
  KBIN_VASP_KPOINTS_KSCHEME                                      = b.KBIN_VASP_KPOINTS_KSCHEME;
  KBIN_VASP_KPOINTS_KSHIFT                                       = b.KBIN_VASP_KPOINTS_KSHIFT;
  // STATIC
  KBIN_VASP_KPOINTS_STATIC_KMODE                                 = b.KBIN_VASP_KPOINTS_STATIC_KMODE;
  KBIN_VASP_KPOINTS_STATIC_KPPRA                                 = b.KBIN_VASP_KPOINTS_STATIC_KPPRA;
  KBIN_VASP_KPOINTS_STATIC_KSCHEME                               = b.KBIN_VASP_KPOINTS_STATIC_KSCHEME;
  KBIN_VASP_KPOINTS_STATIC_KSHIFT                                = b.KBIN_VASP_KPOINTS_STATIC_KSHIFT;
  // PHONONS
  KBIN_VASP_KPOINTS_PHONONS_KPPRA                                = b.KBIN_VASP_KPOINTS_PHONONS_KPPRA;
  KBIN_VASP_KPOINTS_PHONONS_KSCHEME                              = b.KBIN_VASP_KPOINTS_PHONONS_KSCHEME;
  KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY                  = b.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY;
  KBIN_VASP_KPOINTS_PHONONS_GRID                                 = b.KBIN_VASP_KPOINTS_PHONONS_GRID;  //ME20200427
  KBIN_VASP_KPOINTS_PHONONS_SHIFT                                = b.KBIN_VASP_KPOINTS_PHONONS_SHIFT; //ME20200427
  // BANDS
  KBIN_VASP_KPOINTS_BANDS_LATTICE                                = b.KBIN_VASP_KPOINTS_BANDS_LATTICE;
  // [OBSOLETE] KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG                           = b.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG;
  // [OBSOLETE] KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE                          = b.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE;
  //[CO20210805 - OBSOLETE]KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG                      = b.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG;
  KBIN_VASP_KPOINTS_BANDS_GRID                                   = b.KBIN_VASP_KPOINTS_BANDS_GRID;  //CO20210805
  //[CO20210805 - OBSOLETE]KBIN_VASP_KPOINTS_BANDS_GRID_FLAG                              = b.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG;
  //[CO20210805 - OBSOLETE]KBIN_VASP_KPOINTS_BANDS_GRID_VALUE                             = b.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE;
  KBIN_VASP_WRITE_KPOINTS                                        = b.KBIN_VASP_WRITE_KPOINTS;
  // KBIN_VASP_POSCAR_MODE
  KBIN_VASP_POSCAR_MODE                                          = b.KBIN_VASP_POSCAR_MODE;
  KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.clear();
  for(uint i=0;i<b.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++)
    KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.push_back(b.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i));
  KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.clear();
  for(uint i=0;i<b.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size();i++)
    KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.push_back(b.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(i));
  // KBIN_VASP_POTCAR_MODE  
  KBIN_VASP_POTCAR_MODE                                          = b.KBIN_VASP_POTCAR_MODE;
  // KBIN_VASP_INCAR_FILE
  KBIN_VASP_INCAR_FILE                                           = b.KBIN_VASP_INCAR_FILE;
  KBIN_VASP_INCAR_EXPLICIT.clear(); KBIN_VASP_INCAR_EXPLICIT << b.KBIN_VASP_INCAR_EXPLICIT.str(); //ME20181226
  KBIN_VASP_INCAR_EXPLICIT_START_STOP.str(""); KBIN_VASP_INCAR_EXPLICIT_START_STOP << b.KBIN_VASP_INCAR_EXPLICIT_START_STOP.str(); //ME20181226 //CO20190401 - clear ss with .str("")
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_KEYWORD                     = b.KBIN_VASP_INCAR_FILE_KEYWORD;
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_SYSTEM_AUTO                 = b.KBIN_VASP_INCAR_FILE_SYSTEM_AUTO;
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_FILE                        = b.KBIN_VASP_INCAR_FILE_FILE;
  // [OBSOLETE] KBIN_VASP_INCAR_FILE_COMMAND                     = b.KBIN_VASP_INCAR_FILE_COMMAND;
  KBIN_VASP_INCAR_VERBOSE                                        = b.KBIN_VASP_INCAR_VERBOSE;

  // KBIN_VASP_KPOINTS_FILE
  KBIN_VASP_KPOINTS_FILE                                         = b.KBIN_VASP_KPOINTS_FILE;
  KBIN_VASP_KPOINTS_EXPLICIT.clear(); KBIN_VASP_KPOINTS_EXPLICIT << b.KBIN_VASP_KPOINTS_EXPLICIT.str(); //ME20181226
  KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(""); KBIN_VASP_KPOINTS_EXPLICIT_START_STOP << b.KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(); //ME20181226 //CO20190401 - clear ss with .str("")
  // [OBSOLETE] KBIN_VASP_KPOINTS_FILE_KEYWORD                   = b.KBIN_VASP_KPOINTS_FILE_KEYWORD;
  // [OBSOLETE] KBIN_VASP_KPOINTS_FILE_FILE                      = b.KBIN_VASP_KPOINTS_FILE_FILE;
  // [OBSOLETE] KBIN_VASP_KPOINTS_FILE_COMMAND                   = b.KBIN_VASP_KPOINTS_FILE_COMMAND;

  // KBIN_VASP_POSCAR_FILE
  KBIN_VASP_POSCAR_FILE                                          = b.KBIN_VASP_POSCAR_FILE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_KEYWORD                    = b.KBIN_VASP_POSCAR_FILE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_PROTOTYPE                  = b.KBIN_VASP_POSCAR_FILE_PROTOTYPE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_FILE                       = b.KBIN_VASP_POSCAR_FILE_FILE;
  // [OBSOLETE] KBIN_VASP_POSCAR_FILE_COMMAND                    = b.KBIN_VASP_POSCAR_FILE_COMMAND;

  // KBIN_VASP_POSCAR_FILE_VOLUME
  KBIN_VASP_POSCAR_FILE_VOLUME                                   = b.KBIN_VASP_POSCAR_FILE_VOLUME;

  // KBIN_VASP_POTCAR_FILE
  KBIN_VASP_POTCAR_FILE                                          = b.KBIN_VASP_POTCAR_FILE;
  KBIN_VASP_POTCAR_EXPLICIT.clear(); KBIN_VASP_POTCAR_EXPLICIT << b.KBIN_VASP_POTCAR_EXPLICIT.str(); //ME20181226
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_KEYWORD                    = b.KBIN_VASP_POTCAR_FILE_KEYWORD;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_SYSTEM_AUTO                = b.KBIN_VASP_POTCAR_FILE_SYSTEM_AUTO;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_PREFIX                     = b.KBIN_VASP_POTCAR_FILE_PREFIX;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_SUFFIX                     = b.KBIN_VASP_POTCAR_FILE_SUFFIX;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_FILE                       = b.KBIN_VASP_POTCAR_FILE_FILE;
  // [OBSOLETE] KBIN_VASP_POTCAR_FILE_COMMAND                    = b.KBIN_VASP_POTCAR_FILE_COMMAND;
}

// copy
_vflags::_vflags(const _vflags& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _vflags& _vflags::operator=(const _vflags& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _vflags::clear() {free();}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// XVASP
// look into aflow.h for the definitions

// constructors
_xvasp::_xvasp() {
  str                         = xstructure();
  Directory                   = "";
  AnalyzeLabel                = "";
  xqsub                       = _xqsub();
  aopts.clear();
  // VASP INPUT
  aopts.flag("AFLOWIN_FLAG::VASP",TRUE);     // standard VASP
  // [OBSOLETE] AFLOWIN_FLAG_VASP            = TRUE; // standard VASP
  aopts.flag("AFLOWIN_FLAG::QE",FALSE);      // standard non QE
  // [OBSOLETE] AFLOWIN_FLAG_QE              = FALSE; // standard non QE
  aopts.flag("AFLOWIN_FLAG::ABINIT",FALSE);  // standard non ABINIT
  // [OBSOLETE] AFLOWIN_FLAG_ABINIT          = FALSE; // standard non ABINIT
  aopts.flag("AFLOWIN_FLAG::AIMS",FALSE);  // standard non AIMS
  // [OBSOLETE] AFLOWIN_FLAG_AIMS          = FALSE; // standard non AIMS
  POSCAR.str(std::string());
  POSCAR_orig.str(std::string());
  aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
  // [OBSOLETE] POSCAR_generated            = FALSE;
  aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);
  // [OBSOLETE] POSCAR_changed              = FALSE;
  POSCAR_index                = 0;
  INCAR.str(std::string());
  INCAR_orig.str(std::string());
  aopts.flag("FLAG::XVASP_INCAR_generated",FALSE);
  // [OBSOLETE] INCAR_generated             = FALSE;
  aopts.flag("FLAG::XVASP_INCAR_changed",FALSE);
  // [OBSOLETE] INCAR_changed               = FALSE;
  KPOINTS.str(std::string());
  KPOINTS_orig.str(std::string());
  aopts.flag("FLAG::XVASP_KPOINTS_generated",FALSE);
  // [OBSOLETE] KPOINTS_generated           = FALSE;
  aopts.flag("FLAG::XVASP_KPOINTS_changed",FALSE);
  // [OBSOLETE] KPOINTS_changed             = FALSE;
  POTCAR.str(std::string());
  POTCAR_orig.str(std::string());
  aopts.flag("FLAG::XVASP_POTCAR_generated",FALSE);
  // [OBSOLETE] POTCAR_generated            = FALSE;
  aopts.flag("FLAG::XVASP_POTCAR_changed",FALSE);
  // [OBSOLETE] POTCAR_changed              = FALSE;
  POTCAR_TYPE                 = "";
  POTCAR_TYPE_DATE_PRINT_flag = FALSE;
  POTCAR_TYPE_PRINT_flag      = FALSE;  //CO20181226
  POTCAR_ENMAX                = 0.0;
  POTCAR_ENMIN                = 0.0;
  POTCAR_PAW                  = FALSE;
  POTCAR_POTENTIALS.str(std::string());
  POTCAR_AUID.clear();        // md5sums
  NCPUS                       = 1;
  NRELAX                      = 0;
  NRELAXING                   = 0;
  // VASP OUTPUT
  OUTCAR.str(std::string());
  CONTCAR.str(std::string());
  OSZICAR.str(std::string());
  // QE INPUT
  QE_GEOM.str(std::string());
  QE_GEOM_orig.str(std::string());
  QE_GEOM_generated            = FALSE;
  QE_GEOM_changed              = FALSE;
  QE_GEOM_index                = 0;
  // ABINIT INPUT
  // AIMS INPUT
  // AVASP
  //aopts.clear();                                  // default  //CO20180420 - clears AFLOWIN_FLAG::VASP default (TRUE)
  AVASP_aflowin_only_if_missing=FALSE;  // DEFAULT VALUES
  AVASP_arun = false;         //ME20181019
  AVASP_arun_mode = "";       //ME20181019
  AVASP_arun_runname = "";    //ME20181019
  aplopts.clear();            //ME20181025
  aaplopts.clear();           //ME20181025
  qhaopts.clear();            //AS20200302
  AVASP_dirbase="";
  AVASP_libbase="";
  AVASP_label="";
  AVASP_parameters="";
  AVASP_pocc_parameters=""; //CO20181226
  AVASP_pocc_tol=""; //CO20181226
  AVASP_volume_in=0.0;
  AVASP_alpha_fix=FALSE;                          // DEFAULT VALUES
  AVASP_prototype_mode=LIBRARY_MODE_HTQC;         // DEFAULT VALUES
  AVASP_prototype_from_library_=TRUE;             // DEFAULT VALUES
  AVASP_directory_from_library_=TRUE;             // DEFAULT VALUES
  AVASP_potential=DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE;                   // DEFAULT VALUES
  // [OBSOLETE] AVASP_flag_AFLOW_WRITE.clear();          // EMPTY
  aopts.flag("FLAG::AVASP_AUTO_PSEUDOPOTENTIALS",TRUE);  // DEFAULT VALUES
  aopts.flag("AFLOWIN_FLAG::ENMAX_MULTIPLY",FALSE);      // DEFAULT VALUES
  aopts.flag("AFLOWIN_FLAG::NBANDS",TRUE);               // DEFAULT VALUES
  aopts.flag("AFLOWIN_FLAG::POTIM",FALSE);               // DEFAULT VALUES - nothing specified
  aopts.flag("AFLOWIN_FLAG::PSTRESS",FALSE);             // DEFAULT VALUES - nothing specified
  aopts.flag("AFLOWIN_FLAG::EDIFFG",FALSE);              // DEFAULT VALUES - nothing specified
  aopts.flag("AFLOWIN_FLAG::METAGGA",FALSE);              // DEFAULT VALUES - nothing specified
  aopts.flag("AFLOWIN_FLAG::IVDW",FALSE);              // DEFAULT VALUES - nothing specified
  aopts.flag("FLAG::AVASP_SKIP_NOMIX",FALSE);            // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_WAVECAR",DEFAULT_VASP_FORCE_OPTION_WAVECAR);                           // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_CHGCAR",DEFAULT_VASP_FORCE_OPTION_CHGCAR);                             // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_SPIN",DEFAULT_VASP_FORCE_OPTION_SPIN);                                 // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_1",DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1);   // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_SPIN_REMOVE_RELAX_2",DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2);   // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_BADER",DEFAULT_VASP_FORCE_OPTION_BADER);                               // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_ELF",DEFAULT_VASP_FORCE_OPTION_ELF);                                   // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_LSCOUPLING",DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);                     // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_AUTO_MAGMOM",DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);                   // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_SYM",DEFAULT_VASP_FORCE_OPTION_SYM);                   // DEFAULT VALUES //ME20181226
  aopts.flag("FLAG::AVASP_RELAX_FORCES",FALSE);          // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_KPPRA",FALSE);                 // DEFAULT VALUES
  AVASP_value_KPPRA=DEFAULT_KPPRA;                       // DEFAULT VALUES
  AVASP_KSCHEME=DEFAULT_KSCHEME;                         // DEFAULT VALUES //ME20181226
  AVASP_value_KPPRA_STATIC=DEFAULT_KPPRA_STATIC;         // DEFAULT VALUES //ME20181226
  AVASP_STATIC_KSCHEME=DEFAULT_STATIC_KSCHEME;           // DEFAULT VALUES //ME20181226
  AVASP_KPOINTS=""; //ME20181226
  AVASP_value_NSW=0;                                     // DEFAULT VALUES
  aopts.flag("FLAG::PRECISION_SET",FALSE);               // DEFAULT VALUES no PRECISION specified
  AVASP_flag_PRECISION_scheme="N";                       // DEFAULT VALUES Normal PRECISION
  aopts.flag("FLAG::PRECISION_PRESERVED",FALSE);         // DEFAULT VALUES no preservation
  aopts.flag("FLAG::ALGO_SET",FALSE);                    // DEFAULT VALUES no ALGO specified
  AVASP_flag_ALGO_scheme="N";                            // DEFAULT VALUES Normal ALGO
  aopts.flag("FLAG::ALGO_PRESERVED",FALSE);              // DEFAULT VALUES no preservation
  aopts.flag("FLAG::ABMIX_SET",FALSE);                   // DEFAULT VALUES no ABMIX specified
  AVASP_flag_ABMIX_scheme=DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME; // DEFAULT VALUES Normal MIX
  AVASP_flag_TYPE.clear();                               // DEFAULT VALUES
  AVASP_flag_TYPE.xscheme=                                "DEFAULT";
  aopts.flag("FLAG::AVASP_FORCE_LDAU",FALSE);            // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_FORCE_NOLDAU",FALSE);          // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_LDAU1",FALSE);                 // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_LDAU2",FALSE);                 // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_LDAU_ADIABATIC",FALSE);        // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_LDAU_CUTOFF",FALSE);           // DEFAULT VALUES
  AVASP_LDAU_PARAMETERS_STRING="";                       // DEFAULT VALUES
  AVASP_LDAU_PARAMETERS_UJSUM=0.0;                       // DEFAULT VALUES
  aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_STANDARD_PRIMITIVE",TRUE);
  aopts.flag("AVASP_flag_CONVERT_UNIT_CELL_MINKOWSKI",TRUE);
  aopts.flag("FLAG::VOLUME_PRESERVED",FALSE);            // DEFAULT VALUES
  aopts.flag("FLAG::EXTRA_INCAR",FALSE);                 // DEFAULT VALUES
  AVASP_EXTRA_INCAR.str("");                             // DEFAULT VALUES  //ME20190614  clear ss with .str("")
  AVASP_INCAR_KEYWORD.clear();                           // DEFAULT VALUES  //ME20180601
  AVASP_INCAR_EXPLICIT_START_STOP.str("");               // DEFAULT VALUES  //ME20180601  //CO20190401 - clear ss with .str("")
  AVASP_KPOINTS_KEYWORD.clear();                         // DEFAULT VALUES  //ME20180601
  AVASP_KPOINTS_EXPLICIT_START_STOP.str("");             // DEFAULT VALUES  //ME20180601  //CO20190401 - clear ss with .str("")
  AVASP_POTCAR_KEYWORD.clear();                          // DEFAULT VALUES  //CO20170601
  AVASP_flag_MPI=FALSE;                                  // DEFAULT VALUES
  AVASP_flag_RUN_RELAX=TRUE;                             // DEFAULT VALUES
  AVASP_flag_RUN_RELAX_STATIC=FALSE;                     // DEFAULT VALUES
  AVASP_flag_RUN_RELAX_STATIC_BANDS=FALSE;               // DEFAULT VALUES
  AVASP_flag_RUN_STATIC_BANDS=FALSE;                     // DEFAULT VALUES
  AVASP_path_BANDS="";                                   // DEFAULT VALUES
  AVASP_value_BANDS_GRID=0;                              // DEFAULT VALUES
  AVASP_flag_RUN_STATIC=FALSE;                           // DEFAULT VALUES
  AVASP_flag_GENERATE=FALSE;                             // DEFAULT VALUES
  aopts.flag("FLAG::POSCAR_PRESERVED",FALSE);            // DEFAULT VALUES
  aopts.flag("FLAG::CHGCAR_PRESERVED",FALSE);            // DEFAULT VALUES
  aopts.flag("FLAG::KPOINTS_PRESERVED",FALSE);           // DEFAULT VALUES
  aopts.flag("FLAG::WAVECAR_PRESERVED",FALSE);           // DEFAULT VALUES
  aopts.flag("FLAG::WAVEDER_PRESERVED",FALSE);           // DEFAULT VALUES
  aopts.flag("FLAG::AVASP_SYMMETRY=OFF",FALSE);          // DEFAULT VALUES   //CO20170601
  //DX20210122 [OBSOLETE] aopts.flag("FLAG::AVASP_NEIGHBORS=OFF",FALSE);         // DEFAULT VALUES  //CO20170601
  aopts.flag("FLAG::AVASP_APL=OFF",FALSE);               // DEFAULT VALUES  //CO20170601
  aopts.flag("FLAG::AVASP_QHA=OFF",FALSE);               // DEFAULT VALUES  //CO20170601
  aopts.flag("FLAG::AVASP_AAPL=OFF",FALSE);              // DEFAULT VALUES  //CO20170601
  aopts.flag("FLAG::AVASP_AGL=OFF",FALSE);               // DEFAULT VALUES  //CO20170601
  aopts.flag("FLAG::AVASP_AEL=OFF",FALSE);               // DEFAULT VALUES  //CO20170601
}

// destructor
_xvasp::~_xvasp() {
  free();
}

void _xvasp::free() {
}

void _xvasp::clear() {
  _xvasp tmp;
  copy(tmp);
}

void _xvasp::copy(const _xvasp& b) {
  str                                      = b.str;
  Directory                                = b.Directory;
  AnalyzeLabel                             = b.AnalyzeLabel;
  xqsub                                    = b.xqsub;
  aopts                                    = b.aopts;
  // VASP INPUT
  // [OBSOLETE] AFLOWIN_FLAG_VASP                        = b.AFLOWIN_FLAG_VASP;
  POSCAR.str(std::string()); POSCAR << b.POSCAR.str();
  POSCAR_orig.str(std::string()); POSCAR_orig << b.POSCAR_orig.str();
  // [OBSOLETE] POSCAR_generated                         = b.POSCAR_generated;
  // [OBSOLETE] POSCAR_changed                           = b.POSCAR_changed;
  POSCAR_index                             = b.POSCAR_index;
  INCAR.str(std::string()); INCAR << b.INCAR.str();
  INCAR_orig.str(std::string()); INCAR_orig << b.INCAR_orig.str();
  // [OBSOLETE] INCAR_generated                          = b.INCAR_generated;
  // [OBSOLETE] INCAR_changed                            = b.INCAR_changed;
  KPOINTS.str(std::string()); KPOINTS << b.KPOINTS.str();
  KPOINTS_orig.str(std::string()); KPOINTS_orig << b.KPOINTS_orig.str();
  // [OBSOLETE] KPOINTS_generated                        = b.KPOINTS_generated;
  // [OBSOLETE] KPOINTS_changed                          = b.KPOINTS_changed;
  POTCAR.str(std::string()); POTCAR << b.POTCAR.str();
  POTCAR_orig.str(std::string()); POTCAR_orig << b.POTCAR_orig.str();
  // [OBSOLETE] POTCAR_generated                         = b.POTCAR_generated;
  // [OBSOLETE] POTCAR_changed                           = b.POTCAR_changed;
  POTCAR_TYPE                              = b.POTCAR_TYPE;
  POTCAR_TYPE_DATE_PRINT_flag              = b.POTCAR_TYPE_DATE_PRINT_flag;
  POTCAR_TYPE_PRINT_flag                   = b.POTCAR_TYPE_PRINT_flag; //CO20181226
  POTCAR_ENMAX                             = b.POTCAR_ENMAX;
  POTCAR_ENMIN                             = b.POTCAR_ENMIN;
  POTCAR_PAW                               = b.POTCAR_PAW;
  POTCAR_POTENTIALS.str(std::string()); POTCAR_POTENTIALS << b.POTCAR_POTENTIALS.str();
  POTCAR_AUID.clear();for(uint i=0;i<b.POTCAR_AUID.size();i++) POTCAR_AUID.push_back(b.POTCAR_AUID.at(i));
  NCPUS                                    = b.NCPUS;
  NRELAX                                   = b.NRELAX;
  NRELAXING                                = b.NRELAXING;
  // VASP OUTPUT
  OUTCAR.str(std::string()); OUTCAR << b.OUTCAR.str();
  CONTCAR.str(std::string()); CONTCAR << b.CONTCAR.str();
  OSZICAR.str(std::string()); OSZICAR << b.OSZICAR.str();
  // QE INPUT
  // [OBSOLETE] AFLOWIN_FLAG_QE                           = b.AFLOWIN_FLAG_QE;
  QE_GEOM.str(std::string()); QE_GEOM << b.QE_GEOM.str();
  QE_GEOM_orig.str(std::string()); QE_GEOM_orig << b.QE_GEOM_orig.str();
  QE_GEOM_generated                        = b.QE_GEOM_generated;
  QE_GEOM_changed                          = b.QE_GEOM_changed;
  QE_GEOM_index                            = b.QE_GEOM_index;
  // ABINIT INPUT
  // [OBSOLETE] AFLOWIN_FLAG_ABINIT                    = b.AFLOWIN_FLAG_ABINIT;
  // AIMS INPUT
  // [OBSOLETE] AFLOWIN_FLAG_AIMS                      = b.AFLOWIN_FLAG_AIMS;
  // AVASP
  AVASP_aflowin_only_if_missing            = b.AVASP_aflowin_only_if_missing;
  AVASP_arun                               = b.AVASP_arun;            //ME20181019
  AVASP_arun_mode                          = b.AVASP_arun_mode;       //ME20181019
  AVASP_arun_runname                       = b.AVASP_arun_runname;    //ME20181019
  aplopts                                  = b.aplopts;               //ME20181025
  aaplopts                                 = b.aaplopts;              //ME20181025
  qhaopts                                  = b.qhaopts;               //AS20200302
  AVASP_dirbase                            = b.AVASP_dirbase;
  AVASP_libbase                            = b.AVASP_libbase;
  AVASP_label                              = b.AVASP_label;
  AVASP_parameters                         = b.AVASP_parameters;
  AVASP_pocc_parameters                    = b.AVASP_pocc_parameters; //CO20181226
  AVASP_pocc_tol                           = b.AVASP_pocc_tol; //CO20181226
  AVASP_volume_in                          = b.AVASP_volume_in;
  AVASP_potential                          = b.AVASP_potential;
  AVASP_alpha_fix                          = b.AVASP_alpha_fix;
  AVASP_prototype_mode                     = b.AVASP_prototype_mode;
  AVASP_prototype_from_library_            = b.AVASP_prototype_from_library_;
  AVASP_directory_from_library_            = b.AVASP_directory_from_library_;
  AVASP_value_KPPRA                        = b.AVASP_value_KPPRA;
  AVASP_KSCHEME                            = b.AVASP_KSCHEME;
  AVASP_value_KPPRA_STATIC                 = b.AVASP_value_KPPRA_STATIC;
  AVASP_STATIC_KSCHEME                     = b.AVASP_STATIC_KSCHEME;
  AVASP_KPOINTS                            = b.AVASP_KPOINTS;  //ME20180601
  AVASP_value_NSW                          = b.AVASP_value_NSW;
  // [OBSOLETE] AVASP_flag_PRECISION_flag                = b.AVASP_flag_PRECISION_flag;
  AVASP_flag_PRECISION_scheme              = b.AVASP_flag_PRECISION_scheme;
  // [OBSOLETE] AVASP_flag_PRECISION_preserved           = b.AVASP_flag_PRECISION_preserved;
  // [OBSOLETE] AVASP_flag_ALGO_flag                     = b.AVASP_flag_ALGO_flag;
  AVASP_flag_ALGO_scheme                   = b.AVASP_flag_ALGO_scheme;
  // [OBSOLETE] AVASP_flag_ALGO_preserved                = b.AVASP_flag_ALGO_preserved;
  // [OBSOLETE] AVASP_flag_ABMIX_flag                    = b.AVASP_flag_ABMIX_flag;
  AVASP_flag_ABMIX_scheme                  = b.AVASP_flag_ABMIX_scheme;
  AVASP_flag_TYPE                          = b.AVASP_flag_TYPE;
  // [OBSOLETE] AVASP_flag_forceLDAU                     = b.AVASP_flag_forceLDAU;
  // [OBSOLETE] AVASP_flag_forceNOLDAU                   = b.AVASP_flag_forceNOLDAU;
  // [OBSOLETE] AVASP_flag_LDAU1                         = b.AVASP_flag_LDAU1;
  // [OBSOLETE] AVASP_flag_LDAU2                         = b.AVASP_flag_LDAU2;
  // [OBSOLETE] AVASP_flag_LDAU_ADIABATIC                = b.AVASP_flag_LDAU_ADIABATIC;
  // [OBSOLETE] AVASP_flag_LDAU_CUTOFF                   = b.AVASP_flag_LDAU_CUTOFF;
  AVASP_LDAU_PARAMETERS_STRING             = b.AVASP_LDAU_PARAMETERS_STRING;
  AVASP_LDAU_PARAMETERS_UJSUM              = b.AVASP_LDAU_PARAMETERS_UJSUM;
  // [OBSOLETE] AVASP_flag_CONVERT_UNIT_CELL             = b.AVASP_flag_CONVERT_UNIT_CELL;
  // [OBSOLETE] AVASP_flag_PRESERVE_VOLUME               = b.AVASP_flag_PRESERVE_VOLUME;
  // [OBSOLETE] AVASP_flag_EXTRA_INCAR                   = b.AVASP_flag_EXTRA_INCAR;
  AVASP_EXTRA_INCAR.clear(); AVASP_EXTRA_INCAR << b.AVASP_EXTRA_INCAR.str();
  AVASP_INCAR_KEYWORD                      = b.AVASP_INCAR_KEYWORD;  //ME20180601
  AVASP_INCAR_EXPLICIT_START_STOP.str(""); AVASP_INCAR_EXPLICIT_START_STOP << b.AVASP_INCAR_EXPLICIT_START_STOP.str(); //ME20180601 //CO20190401 - clear ss with .str("")
  AVASP_KPOINTS_KEYWORD                    = b.AVASP_KPOINTS_KEYWORD; //ME20180601
  AVASP_KPOINTS_EXPLICIT_START_STOP.str(""); AVASP_KPOINTS_EXPLICIT_START_STOP << b.AVASP_KPOINTS_EXPLICIT_START_STOP.str(); //ME20180601 //CO20190401 - clear ss with .str("")
  AVASP_POTCAR_KEYWORD                    = b.AVASP_POTCAR_KEYWORD; //CO20170601
  // [OBSOLETE] AVASP_flag_LSDAU                         = b.AVASP_flag_LSDAU;
  // [OBSOLETE] AVASP_Get_LDAU2_ParametersU_MODE         = b.AVASP_Get_LDAU2_ParametersU_MODE;
  AVASP_flag_MPI                           = b.AVASP_flag_MPI;
  AVASP_flag_RUN_RELAX                     = b.AVASP_flag_RUN_RELAX;
  AVASP_flag_RUN_RELAX_STATIC              = b.AVASP_flag_RUN_RELAX_STATIC;
  AVASP_flag_RUN_RELAX_STATIC_BANDS        = b.AVASP_flag_RUN_RELAX_STATIC_BANDS;
  AVASP_flag_RUN_STATIC_BANDS              = b.AVASP_flag_RUN_STATIC_BANDS;
  AVASP_path_BANDS                         = b.AVASP_path_BANDS;
  AVASP_value_BANDS_GRID                   = b.AVASP_value_BANDS_GRID;
  AVASP_flag_RUN_STATIC                    = b.AVASP_flag_RUN_STATIC;
  AVASP_flag_GENERATE                      = b.AVASP_flag_GENERATE;
}

// copy
_xvasp::_xvasp(const _xvasp& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _xvasp& _xvasp::operator=(const _xvasp& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

// public functions

double _GetZVAL(const stringstream& sss,vector<double>& vZVAL) {return GetZVAL(sss,vZVAL);} // wrapper
double _GetZVAL(const _xvasp& xvasp,vector<double>& vZVAL) {return GetZVAL(xvasp,vZVAL);} // wrapper
double _GetZVAL(const string& directory,vector<double>& vZVAL) {return GetZVAL(directory,vZVAL);} // wrapper
double _GetCellAtomZVAL(const string& directory,vector<double>& vZVAL,vector<double>& sZVAL,string mode) {return GetCellAtomZVAL(directory,vZVAL,sZVAL,mode);} // wrapper
double _GetCellAtomZVAL(const stringstream& zvalCAR,vector<double>& vZVAL,const stringstream& xstructureCAR,vector<double>& sZVAL,string mode) {return GetCellAtomZVAL(zvalCAR,vZVAL,xstructureCAR,sZVAL,mode);} // wrapper

double _xvasp::GetZVAL(void) {
  vector<double> vZVAL;
  // return GetZVAL(this->POTCAR,vZVAL);
  return _GetZVAL(POTCAR,vZVAL);
}

double _xvasp::GetCellAtomZVAL(string mode) {
  vector<double> vZVAL,sZVAL;
  return _GetCellAtomZVAL(POTCAR,vZVAL,POSCAR,sZVAL,mode);
}

double _GetPOMASS(const stringstream& sss,vector<double>& vPOMASS) {return GetPOMASS(sss,vPOMASS);} // wrapper
double _GetPOMASS(const _xvasp& xvasp,vector<double>& vPOMASS) {return GetPOMASS(xvasp,vPOMASS);} // wrapper
double _GetPOMASS(const string& directory,vector<double>& vPOMASS) {return GetPOMASS(directory,vPOMASS);} // wrapper
double _GetCellAtomPOMASS(const string& directory,vector<double>& vPOMASS,vector<double>& sPOMASS,string mode) {return GetCellAtomPOMASS(directory,vPOMASS,sPOMASS,mode);} // wrapper
double _GetCellAtomPOMASS(const stringstream& pomassCAR,vector<double>& vPOMASS,const stringstream& xstructureCAR,vector<double>& sPOMASS,string mode) {return GetCellAtomPOMASS(pomassCAR,vPOMASS,xstructureCAR,sPOMASS,mode);} // wrapper

double _xvasp::GetPOMASS(void) {
  vector<double> vPOMASS;
  return _GetPOMASS(POTCAR,vPOMASS);
}

double _xvasp::GetCellAtomPOMASS(string mode) {
  vector<double> vPOMASS,sPOMASS;
  return _GetCellAtomPOMASS(POTCAR,vPOMASS,POSCAR,sPOMASS,mode);
}


// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// _AIMS stuff
// look into aflow.h for the definitions
_aimsflags::_aimsflags() {
  KBIN_AIMS_FORCE_OPTION_NOTUNE.clear();
  KBIN_AIMS_FORCE_OPTION_NOTUNE.push("");
  KBIN_AIMS_RUN.clear();
  KBIN_AIMS_RUN.push("");
  KBIN_AIMS_GEOM_MODE.clear();
  KBIN_AIMS_GEOM_MODE.push("");
  KBIN_AIMS_GEOM_FILE.clear();
  KBIN_AIMS_GEOM_FILE.push("");
  KBIN_AIMS_GEOM_FILE_VOLUME.clear();
  KBIN_AIMS_GEOM_FILE_VOLUME.push("");
  KBIN_AIMS_FORCE_OPTION_VOLUME.clear();
  KBIN_AIMS_FORCE_OPTION_VOLUME.push("");
  KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.clear();
  KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL.push("");
  KBIN_AIMS_CONTROL_MODE.clear();
  KBIN_AIMS_CONTROL_MODE.push("");
  KBIN_AIMS_CONTROL_FILE.clear();
  KBIN_AIMS_CONTROL_FILE.push("");
  KBIN_AIMS_CONTROL_VERBOSE=false;

  free();
}

_aimsflags::~_aimsflags() {
  free();
}

void _aimsflags::free() {
  KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.clear();
  for(uint i=0;i<KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.size();i++) {KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE[i].clear();} //DX20191220 - uppercase to lowercase clear
  KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.clear();
}

void _aimsflags::copy(const _aimsflags& b) {
  KBIN_AIMS_FORCE_OPTION_NOTUNE=b.KBIN_AIMS_FORCE_OPTION_NOTUNE;
  KBIN_AIMS_RUN=b.KBIN_AIMS_RUN;
  KBIN_AIMS_GEOM_MODE=b.KBIN_AIMS_GEOM_MODE;
  KBIN_AIMS_GEOM_FILE=b.KBIN_AIMS_GEOM_FILE;
  KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.clear(); for(uint i=0;i<b.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.size();i++) {KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING.push_back(b.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING[i]);}
  KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.clear(); for(uint i=0;i<b.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.size();i++) {KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE.push_back(b.KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE[i]);}
  KBIN_AIMS_GEOM_FILE_VOLUME=b.KBIN_AIMS_GEOM_FILE_VOLUME;
  KBIN_AIMS_FORCE_OPTION_VOLUME=b.KBIN_AIMS_FORCE_OPTION_VOLUME;
  KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL=b.KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL;
  KBIN_AIMS_CONTROL_MODE=b.KBIN_AIMS_CONTROL_MODE;
  KBIN_AIMS_CONTROL_FILE=b.KBIN_AIMS_CONTROL_FILE;
  KBIN_AIMS_CONTROL_VERBOSE=b.KBIN_AIMS_CONTROL_VERBOSE;
}

_aimsflags::_aimsflags(const _aimsflags& b) {
  copy(b);
}

const _aimsflags& _aimsflags::operator=(const _aimsflags& b) {
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _aimsflags::clear() {
  _aimsflags aimsflags_temp;
  copy(aimsflags_temp);
}

///////

_xaims::_xaims() {
  GEOM_index=0;
  Directory="";
  aopts.clear();
  aopts.push("");
  NCPUS=0;
  CONTROL.str("");
  CONTROL_orig.str("");
  CONTROL_generated=false;
  CONTROL_changed=false;
  CONTROL_FILE_NAME="control.in";
  GEOM.str("");
  GEOM_orig.str("");
  GEOM_generated=false;
  GEOM_changed=false;
  GEOM_FILE_NAME="geometry.in";
  OUTPUT_FILE_NAME="aims.out";
  free();
}

_xaims::~_xaims() {
  free();
}

void _xaims::free() {
  str.clear(); //DX20191220 - uppercase to lowercase clear
  xqsub.clear();
}

void _xaims::copy(const _xaims& b) {
  GEOM_index=b.GEOM_index;
  str=b.str;
  Directory=b.Directory;
  xqsub=b.xqsub;
  aopts=b.aopts;
  NCPUS=b.NCPUS;
  CONTROL.str(std::string()); CONTROL << b.CONTROL.str();
  CONTROL_orig.str(std::string()); CONTROL_orig << b.CONTROL_orig.str();
  CONTROL_generated=b.CONTROL_generated;
  CONTROL_changed=b.CONTROL_changed;
  CONTROL_FILE_NAME=b.CONTROL_FILE_NAME;
  GEOM.str(std::string()); GEOM << b.GEOM.str();
  GEOM_orig.str(std::string()); GEOM_orig << b.GEOM_orig.str();
  GEOM_generated=b.GEOM_generated;
  GEOM_changed=b.GEOM_changed;
  GEOM_FILE_NAME=b.GEOM_FILE_NAME;
  OUTPUT_FILE_NAME=b.OUTPUT_FILE_NAME;
}

_xaims::_xaims(const _xaims& b) {
  copy(b);
}

const _xaims& _xaims::operator=(const _xaims& b) {
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _xaims::clear() {
  _xaims aimsflags_temp;
  copy(aimsflags_temp);
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// _AFLIENFLAGS
// look into aflow.h for the definitions

// constructors
_alienflags::_alienflags() {
  KBIN_ALIEN_COMMAND_BINARY_FLAG            = FALSE;
  KBIN_ALIEN_COMMAND_BINARY_VALUE           = DEFAULT_KBIN_ALIEN_BIN;
  KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG = FALSE;
  KBIN_ALIEN_FORCE_OPTION_NOTUNE            = FALSE;
  KBIN_ALIEN_FORCE_OPTION_SOMETHING         = FALSE;
  KBIN_ALIEN_INPUT_MODE_EXPLICIT            = FALSE;
  KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP = FALSE;
  KBIN_ALIEN_INPUT_MODE_IMPLICIT            = FALSE;
  KBIN_ALIEN_INPUT_MODE_EXTERNAL            = FALSE;
  KBIN_ALIEN_INPUT_FILE                     = FALSE;
  KBIN_ALIEN_INPUT_FILE_FILE_FLAG           = FALSE;
  KBIN_ALIEN_INPUT_FILE_FILE_VALUE          = "";
  KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG        = FALSE;
  KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE       = "";
  KBIN_ALIEN_INPUT_MODE_INPUT_FLAG          = FALSE;
  KBIN_ALIEN_INPUT_MODE_INPUT_VALUE         = ALIEN_INPUT_FILE_NAME_DEFAULT;
  KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG        = FALSE;
  KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE       = ALIEN_OUTPUT_FILE_NAME_DEFAULT;
}


// destructor
_alienflags::~_alienflags() {
  free();
}

void _alienflags::free() {
}

void _alienflags::copy(const _alienflags& b) {
  KBIN_ALIEN_COMMAND_BINARY_FLAG            = b.KBIN_ALIEN_COMMAND_BINARY_FLAG;
  KBIN_ALIEN_COMMAND_BINARY_VALUE           = b.KBIN_ALIEN_COMMAND_BINARY_VALUE;
  KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG = b.KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG;
  KBIN_ALIEN_FORCE_OPTION_NOTUNE            = b.KBIN_ALIEN_FORCE_OPTION_NOTUNE;
  KBIN_ALIEN_FORCE_OPTION_SOMETHING         = b.KBIN_ALIEN_FORCE_OPTION_SOMETHING;
  KBIN_ALIEN_INPUT_MODE_EXPLICIT            = b.KBIN_ALIEN_INPUT_MODE_EXPLICIT;
  KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP = b.KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP;
  KBIN_ALIEN_INPUT_MODE_IMPLICIT            = b.KBIN_ALIEN_INPUT_MODE_IMPLICIT;
  KBIN_ALIEN_INPUT_MODE_EXTERNAL            = b.KBIN_ALIEN_INPUT_MODE_EXTERNAL;
  KBIN_ALIEN_INPUT_FILE                     = b.KBIN_ALIEN_INPUT_FILE;
  KBIN_ALIEN_INPUT_FILE_FILE_FLAG           = b.KBIN_ALIEN_INPUT_FILE_FILE_FLAG;
  KBIN_ALIEN_INPUT_FILE_FILE_VALUE          = b.KBIN_ALIEN_INPUT_FILE_FILE_VALUE;
  KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG        = b.KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG;
  KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE       = b.KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE;
  KBIN_ALIEN_INPUT_MODE_INPUT_FLAG          = b.KBIN_ALIEN_INPUT_MODE_INPUT_FLAG;
  KBIN_ALIEN_INPUT_MODE_INPUT_VALUE         = b.KBIN_ALIEN_INPUT_MODE_INPUT_VALUE;
  KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG        = b.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG;
  KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE       = b.KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE;
}

// copy
_alienflags::_alienflags(const _alienflags& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _alienflags& _alienflags::operator=(const _alienflags& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _alienflags::clear() {
  _alienflags alienflags_temp;
  copy(alienflags_temp);
}


// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// XALIEN
// look into aflow.h for the definitions

// constructors
_xalien::_xalien() {
  Directory               = "";
  xqsub                   = _xqsub();
  INPUT.str(std::string());
  INPUT_orig.str(std::string());
  INPUT_generated         = FALSE;
  INPUT_changed           = FALSE;
  INPUT_FILE_NAME         = ALIEN_INPUT_FILE_NAME_DEFAULT;
  OUTPUT_FILE_NAME        = ALIEN_OUTPUT_FILE_NAME_DEFAULT;
  NCPUS                   = 1;
  NRELAX                  = 0;
}

// destructor
_xalien::~_xalien() {
  free();
}

void _xalien::free() {
}

void _xalien::clear() {
  _xalien xalien_temp;
  copy(xalien_temp);
}

void _xalien::copy(const _xalien& b) {
  Directory                    = b.Directory;
  xqsub                        = b.xqsub;
  INPUT.str(std::string());               INPUT << b.INPUT.str();
  INPUT_orig.str(std::string());          INPUT_orig << b.INPUT_orig.str();
  INPUT_generated              = b.INPUT_generated;
  INPUT_changed                = b.INPUT_changed;
  INPUT_FILE_NAME              = b.INPUT_FILE_NAME;
  OUTPUT_FILE_NAME             = b.OUTPUT_FILE_NAME;
  NCPUS                        = b.NCPUS;
  NRELAX                       = b.NRELAX;
}

// copy
_xalien::_xalien(const _xalien& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _xalien& _xalien::operator=(const _xalien& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}


// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// generic classes (for vasp vs. aims vs. alien, etc.)
// look into aflow.h for the definitions
_xflags::_xflags() {
  AFLOW_MODE_VASP=false;
  AFLOW_MODE_AIMS=false;
  AFLOW_MODE_ALIEN=false;
  free();
}
_xflags::_xflags(_vflags& vflags) {setVFlags(vflags);}
_xflags::_xflags(_aimsflags& aimsflags) {setAIMSFlags(aimsflags);}
_xflags::_xflags(_alienflags& alienflags) {setALIENFlags(alienflags);}

_xflags::~_xflags() {
  free();
}

void _xflags::free() {
  vflags.clear();
  aimsflags.clear();
  alienflags.clear();
}

void _xflags::copy(const _xflags& b) {
  AFLOW_MODE_VASP=b.AFLOW_MODE_VASP;
  vflags=b.vflags;
  AFLOW_MODE_AIMS=b.AFLOW_MODE_AIMS;
  aimsflags=b.aimsflags;
  AFLOW_MODE_ALIEN=b.AFLOW_MODE_ALIEN;
  alienflags=b.alienflags;
}

_xflags::_xflags(const _xflags& b) {
  copy(b);
}

const _xflags& _xflags::operator=(const _xflags& b) {
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _xflags::clear() {
  _xflags xflags_temp;
  copy(xflags_temp);
}

void _xflags::setVFlags(_vflags& in_vflags) {
  vflags=in_vflags;
  AFLOW_MODE_VASP=true;
  AFLOW_MODE_AIMS=false;
  AFLOW_MODE_ALIEN=false;
}

void _xflags::setAIMSFlags(_aimsflags& in_aimsflags) {
  aimsflags=in_aimsflags;
  AFLOW_MODE_VASP=false;
  AFLOW_MODE_AIMS=true;
  AFLOW_MODE_ALIEN=false;
}

void _xflags::setALIENFlags(_alienflags& in_alienflags) {
  alienflags=in_alienflags;
  AFLOW_MODE_VASP=false;
  AFLOW_MODE_AIMS=false;
  AFLOW_MODE_ALIEN=true;
}

/////

_xinput::_xinput() {
  AFLOW_MODE_VASP=false;
  AFLOW_MODE_AIMS=false;
  AFLOW_MODE_ALIEN=false;
  free();
}

_xinput::_xinput(_xvasp& xvasp) {setXVASP(xvasp);}
_xinput::_xinput(_xaims& xaims)  {setXAIMS(xaims);}
_xinput::_xinput(_xalien& xalien) {setXALIEN(xalien);}

_xinput::~_xinput() {
  free();
}

void _xinput::free() {
  xvasp.clear();
  xaims.clear();
  xalien.clear();
}

void _xinput::copy(const _xinput& b) {
  AFLOW_MODE_VASP=b.AFLOW_MODE_VASP;
  xvasp=b.xvasp;
  AFLOW_MODE_AIMS=b.AFLOW_MODE_AIMS;
  xaims=b.xaims;
  AFLOW_MODE_ALIEN=b.AFLOW_MODE_ALIEN;
  xalien=b.xalien;
}

_xinput::_xinput(const _xinput& b) {
  copy(b);
}

const _xinput& _xinput::operator=(const _xinput& b) {
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _xinput::clear() {
  _xinput xinput_temp;
  copy(xinput_temp);
}

void _xinput::setXVASP(_xvasp& in_xvasp) {
  xvasp=in_xvasp;
  AFLOW_MODE_VASP=true;
  AFLOW_MODE_AIMS=false;
  AFLOW_MODE_ALIEN=false;
}

void _xinput::setXAIMS(_xaims& in_xaims) {
  xaims=in_xaims;
  AFLOW_MODE_VASP=false;
  AFLOW_MODE_AIMS=true;
  AFLOW_MODE_ALIEN=false;
}

void _xinput::setXALIEN(_xalien& in_xalien) {
  xalien=in_xalien;
  AFLOW_MODE_VASP=false;
  AFLOW_MODE_AIMS=false;
  AFLOW_MODE_ALIEN=true;
}

xstructure& _xinput::getXStr() {
  if(!(AFLOW_MODE_VASP || AFLOW_MODE_AIMS)) {
    throw aurostd::xerror(__AFLOW_FILE__,XPID+"_xinput::getXStr():","No structure available.",_INPUT_MISSING_);
  }
  if(AFLOW_MODE_VASP) {return xvasp.str;}
  if(AFLOW_MODE_AIMS) {return xaims.str;}
  return xvasp.str; //otherwise, return garbage
}

string& _xinput::getDirectory() {
  if(!(AFLOW_MODE_VASP || AFLOW_MODE_AIMS || AFLOW_MODE_ALIEN)) {
    throw aurostd::xerror(__AFLOW_FILE__,XPID+"_xinput::getDirectory():","No directory available.",_INPUT_MISSING_);
  }
  if(AFLOW_MODE_VASP) {return xvasp.Directory;}
  if(AFLOW_MODE_AIMS) {return xaims.Directory;}
  if(AFLOW_MODE_ALIEN) {return xalien.Directory;}
  return xvasp.Directory; //otherwise, return garbage
}

void _xinput::setXStr(const xstructure& str,bool set_all) {
  if(!(AFLOW_MODE_VASP || AFLOW_MODE_AIMS)) {
    throw aurostd::xerror(__AFLOW_FILE__,XPID+"_xinput::setXStr():","No structure available.",_INPUT_MISSING_);
  }
  if(AFLOW_MODE_VASP || set_all) {xvasp.str=str;}
  if(AFLOW_MODE_AIMS || set_all) {xaims.str=str;}
}

void _xinput::setDirectory(string Directory,bool set_all) {
  if(!(AFLOW_MODE_VASP || AFLOW_MODE_AIMS || AFLOW_MODE_ALIEN)) {
    throw aurostd::xerror(__AFLOW_FILE__,XPID+"_xinput::setDirectory():","No directory available.",_INPUT_MISSING_);
  }
  if(AFLOW_MODE_VASP  || set_all) {xvasp.Directory=  Directory;}
  if(AFLOW_MODE_AIMS  || set_all) {xaims.Directory=  Directory;}
  if(AFLOW_MODE_ALIEN || set_all) {xalien.Directory= Directory;}
}

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// XQSUB
// look into aflow.h for the definitions

// constructors
_xqsub::_xqsub() {
  QSUB.str(std::string());
  QSUB_orig.str(std::string());
  QSUB_generated        = FALSE;
  QSUB_changed          = FALSE;
}

// destructor
_xqsub::~_xqsub() {
  free();
}

void _xqsub::free() {
}

void _xqsub::clear() {
  _xqsub xqsub_temp;
  copy(xqsub_temp);
}

void _xqsub::copy(const _xqsub& b) {
  QSUB.str(std::string());              QSUB << b.QSUB.str();
  QSUB_orig.str(std::string());         QSUB_orig << b.QSUB_orig.str();
  QSUB_generated             = b.QSUB_generated;
  QSUB_changed               = b.QSUB_changed;
}

// copy
_xqsub::_xqsub(const _xqsub& b) {
  //  free();
  // *this=b;
  copy(b);
}

// copy operator b=a
const _xqsub& _xqsub::operator=(const _xqsub& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}


// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************

// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
// XSTREAM - Corey Oses - 20180420
// look into aflow.h for the definitions

// constructors
xStream::xStream(ostream& oss) : p_FileMESSAGE(NULL),f_new_ofstream(false) {initialize(oss);} //{free();}
xStream::xStream(ofstream& ofs,ostream& oss) : p_FileMESSAGE(NULL),f_new_ofstream(false) {initialize(ofs,oss);} //{free();}
xStream::~xStream() {free();} //CO20190318
void xStream::free() {
  if(f_new_ofstream) {delete p_FileMESSAGE;}  //first delete, then set to null
  p_FileMESSAGE=NULL;
  f_new_ofstream=false;
  p_oss=NULL;
}
void xStream::copy(const xStream& b) {
  //must handle this very carefully
  //first, since we are interested in copying the stream from b, we must delete the new pointer of self (if it's new)
  //if b is new, don't copy b's pointer, as it doesn't belong to self (and should not be deleted by self)
  //simply create a new one, and declare it as such
  //p_FileMESSAGE=b.p_FileMESSAGE;
  free();
  if(b.f_new_ofstream) {p_FileMESSAGE=(new ofstream());}
  else {setOFStream(*b.p_FileMESSAGE);}  //p_FileMESSAGE=b.p_FileMESSAGE;
  f_new_ofstream=b.f_new_ofstream;  //very important! seg faults otherwise
  setOSS(*b.p_oss); //p_oss=b.p_oss;
}
ostream* xStream::getOSS() const {return p_oss;} //CO20191110
ofstream* xStream::getOFStream() const {return p_FileMESSAGE;} //CO20191110
void xStream::setOFStream(ofstream& FileMESSAGE) {p_FileMESSAGE=&FileMESSAGE;}
void xStream::setOSS(ostream& oss) {p_oss=&oss;}
//ME20200427 START - Initializer functions
void xStream::initialize(ostream& oss) {
  free();
  p_FileMESSAGE=new ofstream();f_new_ofstream=true;
  initialize(*p_FileMESSAGE,oss);
  f_new_ofstream=true;  // override
}
void xStream::initialize(ofstream& ofs,ostream& oss) {
  setOFStream(ofs);f_new_ofstream=false;
  setOSS(oss);
}
//ME20200427 STOP - Initializer functions

#endif  // _AFLOW_CLASSES_CPP


// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
