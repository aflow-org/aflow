// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_INIT_CPP_
#define _AFLOW_INIT_CPP_
#include "aflow.h"

// ***************************************************************************

string _AFLOWIN_; 
string _AFLOWLOCK_; 
const string _LOCK_LINK_SUFFIX_=".init"; //SD20220224
string _STOPFLOW_;  //CO20210315

// THREADS
namespace AFLOW_PTHREADS {
  extern bool FLAG;        // run pthread YES/NO
  extern int MAX;         // how many MAX threads I can use  default or --np
  extern int RUNNING;      // how many threads are actually running
  extern pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  extern int viret[MAX_ALLOCATABLE_PTHREADS];          // the thread runnings
  extern bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
}

#define CPU_File     string("/proc/cpuinfo")

_XHOST XHOST;  // GLOBAL

// vector<string> vVASP_POTCAR_DIRECTORIES;
// vector<string> vAFLOW_LIBRARY_DIRECTORIES;
// vector<string> vAFLOW_PROJECTS_DIRECTORIES;

// uint LIBRARY_LIB2=LIBRARY_NOTHING;
// uint LIBRARY_ICSD=LIBRARY_NOTHING;
// uint LIBRARY_LIB3=LIBRARY_NOTHING;
// uint LIBRARY_LIB4=LIBRARY_NOTHING;
// uint LIBRARY_LIB5=LIBRARY_NOTHING;
// uint LIBRARY_LIB6=LIBRARY_NOTHING;
// uint LIBRARY_LIB7=LIBRARY_NOTHING;
// uint LIBRARY_LIB8=LIBRARY_NOTHING;
// uint LIBRARY_LIB9=LIBRARY_NOTHING;
// uint LIBRARY_LIB1=LIBRARY_NOTHING;
// uint LIBRARY_LIB0=LIBRARY_NOTHING;
// uint LIBRARY_AUID=LIBRARY_NOTHING;

// ***************************************************************************
// init::InitMachine
// ***************************************************************************
namespace init {
  int GetCPUCores() { //CO20180124
    int ncpus=sysconf(_SC_NPROCESSORS_ONLN);
    if(ncpus==96) ncpus=48; // fix the hyperthreading lie
    if(ncpus<1) ncpus=1;
    return ncpus;
  }
  //ME20200724 - returns an int now to remove exits. -1 means that aflow can
  //continue, otherwise terminate with exit code.
  int InitMachine(bool INIT_VERBOSE,vector<string>& argv,vector<string>& cmds,std::ostream& oss) {
    // DECLARATIONS
    bool LDEBUG=(FALSE || XHOST.DEBUG),found=FALSE;
    string message = "";
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [BEGIN]" << endl;

    //CO20220906 WARNING - do not run ANYTHING before the aflowrc read/write/load below
    //XHOST is largely unset before this
    //aflowrc loads TMPDIR needed for many AUROSTD functions

    // AFLOWRC LOAD DEFAULTS FROM AFLOWRC.
    // XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    // XHOST.vflag_control.flag("AFLOWRC::OVERWRITE",aurostd::args2flag(XHOST.argv,cmds,"--aflowrc=overwrite|--aflowrc_overwrite"));
    // getenv returns a char pointer, that can be null if the environment variable is not set
    // assigning a nullptr to a string is undefined behavior
    // always use aurostd::getenv2string that contains the necessary check HE20220701 //SD20220701
    XHOST.home = aurostd::getenv2string("HOME"); //AS SOON AS POSSIBLE
    XHOST.user = aurostd::getenv2string("USER");  //AS SOON AS POSSIBLE
    //CO20220906 START - chicken and egg problem with /tmp file, need a temporary solution for XHOST.tmpfs as we need to write out tmp files immediately
    XHOST.argv.clear();for(uint i=0;i<argv.size();i++){XHOST.argv.push_back(argv.at(i));} //AS SOON AS POSSIBLE
    vector<string> tokens;
    string tmpfs_str=AFLOWRC_DEFAULT_TMPFS_DIRECTORIES; //CO20220906 - always defined, does not require reading aflow.rc
    const string tmpfs_str_input=aurostd::args2attachedstring(XHOST.argv,"--use_tmpfs=","");
    if(!tmpfs_str_input.empty()){tmpfs_str=tmpfs_str_input;}  //override with input from command line, always takes precedence
    tmpfs_str+=",./";    //CO20220906 - hack, add current directory at the end, should work unless there's a real problem
    aurostd::string2tokens(tmpfs_str,tokens,",");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " tokens[tmpfs_str]=\"" << aurostd::joinWDelimiter(tokens,",") << "\"" << endl;}
    for(uint i=0;i<tokens.size() && XHOST.tmpfs.empty();i++){
      if(aurostd::DirectoryWritable(tokens[i])){XHOST.tmpfs=tokens[i];}
    }
    if(XHOST.tmpfs.empty()){
      message="tmp directories are not writable";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }
    XHOST.tmpfs=aurostd::CleanFileName(XHOST.tmpfs+"/");
    //
    if(!aflowrc::is_available(oss,INIT_VERBOSE || XHOST.DEBUG)) aflowrc::write_default(oss,INIT_VERBOSE || XHOST.DEBUG);
    aflowrc::read(oss,INIT_VERBOSE || XHOST.DEBUG);
    XHOST.vflag_control.flag("AFLOWRC::READ",aurostd::args2flag(XHOST.argv,cmds,"--aflowrc=read|--aflowrc_read"));
    if(XHOST.vflag_control.flag("AFLOWRC::READ")) {aflowrc::print_aflowrc(oss,TRUE);return false;}
    //
    // SD20220223 - read AFLOWRC to determine the temporary directory, which is needed for execute2string
    // SD20220223 - check to make sure the temporary directory is writable
    tmpfs_str=DEFAULT_TMPFS_DIRECTORIES;  //CO20220906 - redefine with user aflow.rc parameters
    if(!tmpfs_str_input.empty()){tmpfs_str=tmpfs_str_input;}  //override with input from command line (already snatched before), always takes precedence
    aurostd::string2tokens(tmpfs_str,tokens,",");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " tokens[tmpfs_str]=\"" << aurostd::joinWDelimiter(tokens,",") << "\"" << endl;}
    for(uint i=0;i<tokens.size() && XHOST.tmpfs.empty();i++){
      if(aurostd::DirectoryWritable(tokens[i])){XHOST.tmpfs=tokens[i];}
    }
    if(XHOST.tmpfs.empty()) {
      message="tmp directories are not writable";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }
    XHOST.tmpfs=aurostd::CleanFileName(XHOST.tmpfs+"/");
    //[SD20220223 - OBSOLETE]XHOST.tmpfs=aurostd::args2attachedstring(argv,"--use_tmpfs","/tmp");
    tokens.clear();
    // SD20220223 - try alternatives
    if(XHOST.home.empty()){XHOST.home=aurostd::execute2string("cd && pwd");}
    if(XHOST.home.empty()){XHOST.home="~";}
    if(XHOST.user.empty()){XHOST.user=aurostd::execute2string("whoami");}
    if(XHOST.user.empty()){XHOST.user="UNKNOWN_USER";}

    int depth_short=20,depth_long=45;
    string position;
    if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
    if(INIT_VERBOSE) oss << "* AFLOW V=" << string(AFLOW_VERSION) << " - machine information " << endl;
    if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;

    // IMMEDIATELY
    XHOST.QUIET=XHOST.QUIET_GLOBAL=aurostd::args2flag(argv,cmds,"--quiet|-q");  //CO20220630 - global allows headless server, do not modify inside
    XHOST.QUIET_CERR=aurostd::args2flag(argv,cmds,"--quiet=cerr"); // extra quiet SC20210617
    XHOST.QUIET_COUT=aurostd::args2flag(argv,cmds,"--quiet=cout"); // extra quiet SC20210617
    XHOST.DEBUG=aurostd::args2flag(argv,cmds,"--debug");
    XHOST.TEST=aurostd::args2flag(argv,cmds,"--test|-test");
    XHOST.SKEW_TEST=aurostd::args2flag(argv,cmds,"--skew_test"); //DX20171025
    XHOST.READ_SPIN_FROM_ATOMLABEL=aurostd::args2flag(argv,cmds,"--read_spin_from_atomlabel"); //SD20220316
    //[CO20200404 - overload with --www]XHOST.WEB_MODE=aurostd::args2flag(argv,cmds,"--web_mode"); //CO20190402
    XHOST.MPI=aurostd::args2flag(argv,"--MPI|--mpi");

    XHOST.GENERATE_AFLOWIN_ONLY=aurostd::args2flag(argv,cmds,"--generate_aflowin_only");  //CT20180719
    XHOST.POSTPROCESS=aurostd::args2attachedflag(argv,cmds,"--lib2raw=|--lib2lib=");  //CO20200624
    XHOST.ARUN_POSTPROCESS=aurostd::args2flag(argv,cmds,"--postprocess");  //CT20181212

    // IMMEDIATELY GET PIDS
    // get PID
    XHOST.PID=getpid();    // PID number
    aurostd::StringstreamClean(XHOST.ostrPID);
    XHOST.ostrPID<<XHOST.PID;  // PID as stringstream
    XHOST.sPID="";
    XHOST.showPID=aurostd::args2flag(argv,cmds,"--showPID");
    if(XHOST.showPID) XHOST.sPID="[PID="+aurostd::PaddedPRE(XHOST.ostrPID.str(),7)+"] "; // PID as a comment
    aurostd::xerror_PID=XHOST.sPID;

    // get PGID
    XHOST.PGID=getpgrp();    // PGID number
    aurostd::StringstreamClean(XHOST.ostrPGID);
    XHOST.ostrPGID<<XHOST.PGID;  // PGID as stringstream
    XHOST.sPGID="";
    XHOST.showPGID=aurostd::args2flag(argv,cmds,"--showPGID");
    if(XHOST.showPGID) XHOST.sPGID="[PGID="+aurostd::PaddedPRE(XHOST.ostrPGID.str(),7)+"] "; // PGID as a comment

    // get TID
    XHOST.TID=getpid();    // TID number
    aurostd::StringstreamClean(XHOST.ostrTID);
    XHOST.ostrTID<<XHOST.TID;  // TID as stringstream
    XHOST.sTID="";
    XHOST.showTID=aurostd::args2flag(argv,cmds,"--showTID");
    if(XHOST.showTID) XHOST.sTID="[TID="+aurostd::PaddedPRE(XHOST.ostrTID.str(),7)+"] "; // TID as a comment

    //    if(XHOST.showPID) XHOST.sPID="[TID="+aurostd::PaddedPRE(XHOST.ostrTID.str(),7)+"] "; // TID as a comment  //CO20200502 - threadID
    //   XHOST.sPID="[PID="+aurostd::PaddedPRE(XHOST.ostrPID.str(),7)+"] "; // for the time being (LIB4)

    // DO THREADS IMMEDIATELY
    // std::vector<pthread_t> _thread(MAX_ALLOCATABLE_PTHREADS);AFLOW_PTHREADS::vpthread=_thread;
    // std::vector<int>    _iret(MAX_ALLOCATABLE_PTHREADS);AFLOW_PTHREADS::viret=_iret;
    // std::vector<bool>   _thread_busy(MAX_ALLOCATABLE_PTHREADS);AFLOW_PTHREADS::vpthread_busy=_thread_busy;
    // AFLOW_PTHREADS::vpthread.clear();AFLOW_PTHREADS::vpthread=*(new vector<pthread_t>(MAX_ALLOCATABLE_PTHREADS));      // they are static now
    // AFLOW_PTHREADS::viret.clear();AFLOW_PTHREADS::viret=*(new vector<int>(MAX_ALLOCATABLE_PTHREADS));                  // they are static now
    // AFLOW_PTHREADS::vpthread_busy.clear();AFLOW_PTHREADS::vpthread_busy=*(new vector<bool>(MAX_ALLOCATABLE_PTHREADS)); // they are static now
    // CPUS
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [1]" << endl;
    XHOST.CPU_Cores=GetCPUCores();//sysconf(_SC_NPROCESSORS_ONLN);
    //if(XHOST.CPU_Cores==96) XHOST.CPU_Cores=48; // fix the hyperthreading lie
    //if(XHOST.CPU_Cores<1) XHOST.CPU_Cores=1;
    XHOST.CPU_Model="nan";
    XHOST.CPU_MHz="nan";
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [2]" << endl;
    // #ifndef _MACOSX_
    if(aurostd::FileExist(CPU_File)) {  // LINUX SYSTEMS
      // XHOST.CPU_Model
      aurostd::string2vectorstring(aurostd::execute2string("cat "+CPU_File+" | grep \"model name\""),tokens);
      if(tokens.size()>0) {
        aurostd::StringSubst(tokens.at(0),": ",":");aurostd::string2tokens(string(tokens.at(0)),tokens,":");
        if(tokens.size()>1) {
          aurostd::StringSubst(tokens.at(1)," ","_");aurostd::StringSubst(tokens.at(1),"__","_");
          aurostd::StringSubst(tokens.at(1),"__","_");aurostd::StringSubst(tokens.at(1),"__","_");
          XHOST.CPU_Model=tokens.at(1);
        }
      }
      // XHOST.CPU_MHz
      aurostd::string2vectorstring(aurostd::execute2string("cat "+CPU_File+" | grep \"cpu MHz\""),tokens);
      if(tokens.size()>0) {
        aurostd::StringSubst(tokens.at(0),": ",":");aurostd::string2tokens(string(tokens.at(0)),tokens,":");
        if(tokens.size()>1) {
          aurostd::StringSubst(tokens.at(1)," ","_");XHOST.CPU_MHz=aurostd::utype2string(ceil(aurostd::string2utype<double>(tokens.at(1))));
        }
      }
    }
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [3]" << endl;
    // #endif
    // if(INIT_VERBOSE) {
    // oss << "--- MACHINE ------------ " << endl;
    // oss << aurostd::PaddedPOST("NCPUS = ",depth_short) << XHOST.CPU_Cores << endl;
    // oss << aurostd::PaddedPOST("random_seed = ",depth_short) << random_seed << endl;
    //}
    // MEMORY
    XHOST.RAM=init::GetRAM();
    XHOST.RAM_MB=XHOST.RAM/1024/1024;
    XHOST.RAM_GB=ceil(XHOST.RAM/1024/1024/1024);
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [4]" << endl;
    // RANDOM NUMBER GENERATOR
    long int random_seed=aurostd::_random_initialize();
    // oss << floor(10000*aurostd::ran0()) << " " << time << endl;
    // TMP FILE system
    //    XHOST.tmpfs="/tmp/"; if(aurostd::FileExist("/run/shm/")) if(aurostd::DirectoryWritable("/run/shm/")) XHOST.tmpfs="/run/shm/";
    if(INIT_VERBOSE) oss << "XHOST.tmpfs=" << XHOST.tmpfs << endl;
    // AFLOW_TIME
    XHOST.Day=aurostd::utype2string(aurostd::get_day());
    XHOST.Month=aurostd::utype2string(aurostd::get_month());
    XHOST.Year=aurostd::utype2string(aurostd::get_year());
    XHOST.Copyright_Years="2003-"+aurostd::utype2string(aurostd::get_year());
    XHOST.Time_starting=0.0;
    XHOST.Time_now=0.0;
    XHOST.Time_starting=aurostd::get_seconds();
    XHOST.Time_now=XHOST.Time_starting;
    // AFLOW_DATE
    XHOST.Date=aurostd::get_date();
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [5]" << endl;
    // some verbose
    if(INIT_VERBOSE) {
      oss << aurostd::PaddedPOST("XHOST.ostrPID = ",depth_short) << XHOST.ostrPID.str() << endl;
      oss << aurostd::PaddedPOST("XHOST.ostrPGID = ",depth_short) << XHOST.ostrPGID.str() << endl;
      oss << aurostd::PaddedPOST("XHOST.ostrTID = ",depth_short) << XHOST.ostrTID.str() << endl;  //CO20200502 - threadID
      oss << aurostd::PaddedPOST("DEFAULT_KZIP_BIN = ",depth_short) << DEFAULT_KZIP_BIN << endl;
      oss << aurostd::PaddedPOST("DEFAULT_KZIP_EXT = ",depth_short) << DEFAULT_KZIP_EXT << endl;
      oss << "--- MACHINE ------------ " << endl;
      oss << aurostd::PaddedPOST("CPU_Model = ",depth_short) << XHOST.CPU_Model << endl;
      oss << aurostd::PaddedPOST("CPU_MHz = ",depth_short) << XHOST.CPU_MHz << endl;
      oss << aurostd::PaddedPOST("CPU_Cores = ",depth_short) << XHOST.CPU_Cores << endl;
      oss << aurostd::PaddedPOST("RAM = ",depth_short) << (long int) XHOST.RAM << endl;
      oss << aurostd::PaddedPOST("RAM_MB = ",depth_short) << XHOST.RAM_MB << " (" << XHOST.RAM_GB << "GB)" << endl;
      oss << aurostd::PaddedPOST("random_seed = ",depth_short) << random_seed << endl;
    }

    // NAME and OS
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [6]" << endl;
    static struct utsname os;
    if((uname(&os)) < 0) { return FALSE;} // NULL
    XHOST.hostname=os.nodename;//aurostd::execute2string(string("hostname"));
    if(XHOST.hostname=="nietzsche") XHOST.hostname="nietzsche.mems.duke.edu";
    if(XHOST.hostname=="materials") XHOST.hostname="materials.duke.edu";
    if(XHOST.hostname=="aflowlib") XHOST.hostname="aflowlib.duke.edu";
    if(XHOST.hostname=="quser") XHOST.hostname="quser.materials.duke.edu";  //CO20200526
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("hostname = ",depth_short) << XHOST.hostname << endl;
    if(AFLOW_BlackList(XHOST.hostname)) {
      message = "HOSTNAME BLACKLISTED = " + XHOST.hostname;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
    }

    // MACHINE TYPE
    XHOST.machine_type="linux";
    if(aurostd::substring2bool(os.sysname,"Darwin")) XHOST.machine_type="macosx";
    if(aurostd::substring2bool(os.sysname,"Linux")) XHOST.machine_type="linux";
    if(aurostd::substring2bool(os.sysname,"OSF1")) XHOST.machine_type="alpha";
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("machinetype = ",depth_short) << XHOST.machine_type << endl;
    // SERVER STUFF
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [7]" << endl;
    XHOST.AFLOW_MATERIALS_SERVER=AFLOW_MATERIALS_SERVER_DEFAULT;XHOST.AFLOW_WEB_SERVER=AFLOW_WEB_SERVER_DEFAULT; // DEFAULT
    if(aurostd::substring2bool(XHOST.hostname,"nietzsche")) {XHOST.AFLOW_MATERIALS_SERVER=AFLOW_MATERIALS_SERVER_DEFAULT;XHOST.AFLOW_WEB_SERVER=AFLOW_WEB_SERVER_DEFAULT;}
    if(aurostd::substring2bool(XHOST.hostname,"aflowlib")) {XHOST.AFLOW_MATERIALS_SERVER="aflowlib.duke.edu";XHOST.AFLOW_WEB_SERVER="aflowlib.duke.edu";}
    if(INIT_VERBOSE) {
      oss << "--- SERVER ------------------ " << endl;
      oss << aurostd::PaddedPOST("XHOST.AFLOW_MATERIALS_SERVER = ",depth_short) << XHOST.AFLOW_MATERIALS_SERVER << endl;
      oss << aurostd::PaddedPOST("XHOST.AFLOW_WEB_SERVER = ",depth_short) << XHOST.AFLOW_WEB_SERVER << endl;
    }
    // fix FIND --noleaf or not noleaf
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [8]" << endl;
    XHOST.Find_Parameters=DEFAULT_AFLOW_FIND_PARAMETERS_NOLEAF;
    if(XHOST.machine_type=="macosx") XHOST.Find_Parameters=DEFAULT_AFLOW_FIND_PARAMETERS_NORMAL;
    if(INIT_VERBOSE) {
      oss << "--- OS ------------------ " << endl;
      oss << aurostd::PaddedPOST("os.sysname = ",depth_short) << os.sysname << endl;
      oss << aurostd::PaddedPOST("os.release = ",depth_short) << os.release << endl;
      oss << aurostd::PaddedPOST("os.version = ",depth_short) << os.version << endl;
      oss << aurostd::PaddedPOST("os.machine = ",depth_short) << os.machine << endl;
      oss << aurostd::PaddedPOST("os.nodename = ",depth_short) << os.nodename << endl;
    }
    // USER //  uid_t uid;uid=geteuid();
    // OLD
    // struct passwd *pw; // pw=getpwuid(uid);
    // XHOST.user=string(pw->pw_name);
    // NEW
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [9]" << endl;
    // XHOST.user="boot";
    aurostd::StringSubst(XHOST.user,"\n","");
    // USER DONE
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("username = ",depth_short) << XHOST.user << endl;
    // GROUP
    aurostd::string2tokens(aurostd::execute2string("groups 2> /dev/null"),tokens);
    XHOST.group="none";
    if(tokens.size()>0) XHOST.group=tokens.at(0);
    aurostd::StringSubst(XHOST.group,"\n","");
    // OLD
    // struct group *grp; // grp=getgrgid(pw->pw_gid);
    // XHOST.group=string(grp->gr_name);
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("groupname = ",depth_short) << XHOST.group << endl;
    // HOME DONE
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("home = ",depth_short) << XHOST.home << endl;
    // SHELL
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [10]" << endl;
    XHOST.shell=aurostd::execute2string(string("echo $SHELL"));
    aurostd::string2tokens(XHOST.shell,tokens,"/");
    if(tokens.size()>0) XHOST.shell=tokens.at(tokens.size()-1);
    aurostd::StringSubst(XHOST.shell,"\n","");
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("shell = ",depth_short) << XHOST.shell << endl;
    // PROGNAME
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [11]" << endl;
    XHOST.progname="aflow";
    if(aurostd::substring2bool(XHOST.argv.at(0),"aconvasp") || aurostd::substring2bool(XHOST.argv.at(0),"convasp")) XHOST.progname="aconvasp";
    if(aurostd::substring2bool(XHOST.argv.at(0),"aflow") || aurostd::substring2bool(XHOST.argv.at(0),"aflowd")) XHOST.progname="aflow";
    if(aurostd::substring2bool(XHOST.argv.at(0),"aflow1") || aurostd::substring2bool(XHOST.argv.at(0),"aflowd1")) XHOST.progname="aflow1";
    if(aurostd::substring2bool(XHOST.argv.at(0),"aflow2") || aurostd::substring2bool(XHOST.argv.at(0),"aflowd2")) XHOST.progname="aflow2";
    if(aurostd::substring2bool(XHOST.argv.at(0),"apennsy") || aurostd::substring2bool(XHOST.argv.at(0),"apennsy")) XHOST.progname="apennsy";
    if(INIT_VERBOSE) oss << "--- PROGNAME ------------ " << endl;
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("progname = ",depth_short) << XHOST.progname << endl;

    // IP
    // ifconfig | grep inet | grep -v 127 | grep -v inet6
    // GET PROGRAMS
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [12]" << endl;
    if(INIT_VERBOSE) oss << "--- COMPILER ------------ " << endl;
    if(INIT_VERBOSE) oss << aurostd::PaddedPOST("G++ version = ",depth_short) << GCC_VERSION << endl;
    if(INIT_VERBOSE) oss << "--- BINARIES @ " << XHOST.hostname << " --- " << endl;
    // GET AFLOW_DATA
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [13]" << endl;
    found=FALSE;
    //CO20180706 - above we override Progname ./aflow with aflow, so it will never look in ./aflow
    //I am assuming this is for a good reason, but it screws things up here
    //so we need to rely on XHOST.argv.at(0) as above
    if(!found&&(XHOST.argv[0]=="aflow"||XHOST.argv[0]=="aflowd"||XHOST.argv[0]=="aconvasp"||XHOST.argv[0]=="apennsy")) {
      string aflow_data_command="aflow_data";
      //CO20180706 - note, IsCommandAvailableModify() simply overwrites aflow_data_command with true path from `which'
      if(aurostd::IsCommandAvailableModify(aflow_data_command)) {XHOST.vcmd.push_back(aflow_data_command);found=TRUE;}  //CO20180703
    }
    if(!found&&(XHOST.argv[0]=="./aflow"||XHOST.argv[0]=="./aflowd"||XHOST.argv[0]=="./aconvasp"||XHOST.argv[0]=="./apennsy")) {
      if(aurostd::FileExist("./aflow_data")) {XHOST.vcmd.push_back("./aflow_data");found=TRUE;}
    }
    if(!found&&(XHOST.argv[0]=="/usr/local/bin/aflow"||XHOST.argv[0]=="/usr/local/bin/aflowd"||XHOST.argv[0]=="/usr/local/bin/aconvasp"||XHOST.argv[0]=="/usr/local/bin/apennsy")) {
      if(aurostd::FileExist("/usr/local/bin/aflow_data")) {XHOST.vcmd.push_back("/usr/local/bin/aflow_data");found=TRUE;}
    }
    //[OBSOLETE CO20180706]if(!found&&(XHOST.progname=="aflow"||XHOST.progname=="aflowd"||XHOST.progname=="aconvasp"||XHOST.progname=="apennsy")) {XHOST.vcmd.push_back("aflow_data");found=TRUE;}
    //[OBSOLETE CO20180706]if(!found&&(XHOST.progname=="./aflow"||XHOST.progname=="./aflowd"||XHOST.progname=="./aconvasp"||XHOST.progname=="./apennsy")) {XHOST.vcmd.push_back("./aflow_data");found=TRUE;}
    //[OBSOLETE CO20180706]if(!found&&(XHOST.progname=="/usr/local/bin/aflow"||XHOST.progname=="/usr/local/bin/aflowd"||XHOST.progname=="/usr/local/bin/aconvasp"||XHOST.progname=="/usr/local/bin/apennsy")) {XHOST.vcmd.push_back("/usr/local/bin/aflow_data");found=TRUE;}
    if(!found) {
      aurostd::string2tokens(XHOST.argv.at(0),tokens,"/");  //CO20180703 - note, string2tokens without consecutive keeps beginning / with first entry
      string aflow_data="";
      for(uint i=0;i<tokens.size()-1;i++) {aflow_data+=tokens.at(i)+"/";} aflow_data+=string("aflow_data");
      if(aurostd::FileExist(aflow_data)) {
        XHOST.vcmd.push_back(aflow_data);
        for(uint i=0;i<XHOST.vcmd.size();i++)
          if(aurostd::substring2bool(XHOST.vcmd.at(i),"aflow_data")) 
            XHOST.vcmd.at(i)=aflow_data;
      }
    }
    // search for updates and proxies
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [14]" << endl;

    if(XHOST.hostname=="qrats.materials.duke.edu" && aurostd::FileExist("/usr/local/maui/bin/showq")){XHOST.vcmd.push_back("/usr/local/maui/bin/showq");}  //CO20200526

    // if(XHOST.is_command("wget")) {aurostd::execute("wget -q http://materials.duke.edu/aflow_update/"+XHOST.user+"/"+XHOST.hostname);};

    // SOME LOADING UP
    if(INIT_VERBOSE) {
      if(XHOST.is_command("aflow_data")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"aflow_data\")=TRUE",depth_long) << "[" << XHOST.command("aflow_data") << "]" << endl;} else {oss << "XHOST.is_command(\"aflow_data\")=FALSE" << endl;}
      if(XHOST.is_command("beep")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"beep\")=TRUE",depth_long) << "[" << XHOST.command("beep") << "]" << endl;} else {oss << "XHOST.is_command(\"beep\")=FALSE" << endl;}
      if(XHOST.is_command("bkill")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"bkill\")=TRUE",depth_long) << "[" << XHOST.command("bkill") << "]" << endl;} else {oss << "XHOST.is_command(\"bkill\")=FALSE" << endl;}
      if(XHOST.is_command("bsub")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"bsub\")=TRUE",depth_long) << "[" << XHOST.command("bsub") << "]" << endl;} else {oss << "XHOST.is_command(\"bsub\")=FALSE" << endl;}
      if(XHOST.is_command("bzip2")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"bzip2\")=TRUE",depth_long) << "[" << "bzip2" << "]" << endl;} else {oss << "XHOST.is_command(\"bzip2\")=FALSE" << endl;}
      if(XHOST.is_command("compress")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"compress\")=TRUE",depth_long) << "[" << XHOST.command("compress") << "]" << endl;} else {oss << "XHOST.is_command(\"compress\")=FALSE" << endl;}
      if(XHOST.is_command("convert")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"convert\")=TRUE",depth_long) << "[" << XHOST.command("convert") << "]" << endl;} else {oss << "XHOST.is_command(\"convert\")=FALSE" << endl;}
      if(XHOST.is_command("curl")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"curl\")=TRUE",depth_long) << "[" << XHOST.command("curl") << "]" << endl;} else {oss << "XHOST.is_command(\"curl\")=FALSE" << endl;}
      if(XHOST.is_command("dvipdf")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"dvipdf\")=TRUE",depth_long) << "[" << XHOST.command("dvipdf") << "]" << endl;} else {oss << "XHOST.is_command(\"dvipdf\")=FALSE" << endl;}
      if(XHOST.is_command("dvips")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"dvips\")=TRUE",depth_long) << "[" << XHOST.command("dvips") << "]" << endl;} else {oss << "XHOST.is_command(\"dvips\")=FALSE" << endl;}
      if(XHOST.is_command("findsym")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"findsym\")=TRUE",depth_long) << "[" << XHOST.command("findsym") << "]" << endl;} else {oss << "XHOST.is_command(\"findsym\")=FALSE" << endl;}
      if(XHOST.is_command("frozsl")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"frozsl\")=TRUE",depth_long) << "[" << XHOST.command("frozsl") << "]" << endl;} else {oss << "XHOST.is_command(\"frozsl\")=FALSE" << endl;}
      if(XHOST.is_command("frozsl_init")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"frozsl_init\")=TRUE",depth_long) << "[" << XHOST.command("frozsl_init") << "]" << endl;} else {oss << "XHOST.is_command(\"frozsl_init\")=FALSE" << endl;}
      if(XHOST.is_command("gnuplot")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"gnuplot\")=TRUE",depth_long) << "[" << XHOST.command("gnuplot") << "]" << endl;} else {oss << "XHOST.is_command(\"gnuplot\")=FALSE" << endl;}
      if(XHOST.is_command("gzip")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"gzip\")=TRUE",depth_long) << "[" << "gzip" << "]" << endl;} else {oss << "XHOST.is_command(\"gzip\")=FALSE" << endl;}
      if(XHOST.is_command("halt")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"halt\")=TRUE",depth_long) << "[" << XHOST.command("halt") << "]" << endl;} else {oss << "XHOST.is_command(\"halt\")=FALSE" << endl;}
      if(XHOST.is_command("kill")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"kill\")=TRUE",depth_long) << "[" << XHOST.command("kill") << "]" << endl;} else {oss << "XHOST.is_command(\"kill\")=FALSE" << endl;}
      if(XHOST.is_command("ifconfig")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"ifconfig\")=TRUE",depth_long) << "[" << XHOST.command("ifconfig") << "]" << endl;} else {oss << "XHOST.is_command(\"ifconfig\")=FALSE" << endl;}
      if(XHOST.is_command("java")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"java\")=TRUE",depth_long) << "[" << XHOST.command("java") << "]" << endl;} else {oss << "XHOST.is_command(\"java\")=FALSE" << endl;}
      if(XHOST.is_command("jmol")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"jmol\")=TRUE",depth_long) << "[" << XHOST.command("jmol") << "]" << endl;} else {oss << "XHOST.is_command(\"jmol\")=FALSE" << endl;}
      if(XHOST.is_command("latex")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"latex\")=TRUE",depth_long) << "[" << XHOST.command("latex") << "]" << endl;} else {oss << "XHOST.is_command(\"latex\")=FALSE" << endl;}
      if(XHOST.is_command("matlab")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"matlab\")=TRUE",depth_long) << "[" << XHOST.command("matlab") << "]" << endl;} else {oss << "XHOST.is_command(\"matlab\")=FALSE" << endl;}
      if(XHOST.is_command("mpivasp46s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp46s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp46s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp46s\")=FALSE" << endl;}
      if(XHOST.is_command("mpivasp52s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp52s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp52s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp52s\")=FALSE" << endl;}
      if(XHOST.is_command("mpivasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp54s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp54s\")=FALSE" << endl;}
      if(XHOST.is_command("mpivasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"mpivasp54s\")=TRUE",depth_long) << "[" << XHOST.command("mpivasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"mpivasp54s\")=FALSE" << endl;}
      if(XHOST.is_command("pbsnodes")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"pbsnodes\")=TRUE",depth_long) << "[" << XHOST.command("pbsnodes") << "]" << endl;} else {oss << "XHOST.is_command(\"pbsnodes\")=FALSE" << endl;} //CO20200526
      if(XHOST.is_command("pdflatex")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"pdflatex\")=TRUE",depth_long) << "[" << XHOST.command("pdflatex") << "]" << endl;} else {oss << "XHOST.is_command(\"pdflatex\")=FALSE" << endl;}
      if(XHOST.is_command("platon")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"platon\")=TRUE",depth_long) << "[" << XHOST.command("platon") << "]" << endl;} else {oss << "XHOST.is_command(\"platon\")=FALSE" << endl;}
      if(XHOST.is_command("ps2pdf")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"ps2pdf\")=TRUE",depth_long) << "[" << XHOST.command("ps2pdf") << "]" << endl;} else {oss << "XHOST.is_command(\"ps2pdf\")=FALSE" << endl;}
      if(XHOST.is_command("pwd")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"pwd\")=TRUE",depth_long) << "[" << XHOST.command("pwd") << "]" << endl;} else {oss << "XHOST.is_command(\"pwd\")=FALSE" << endl;}
      if(XHOST.is_command("qconvex")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qconvex\")=TRUE",depth_long) << "[" << XHOST.command("qconvex") << "]" << endl;} else {oss << "XHOST.is_command(\"qconvex\")=FALSE" << endl;}
      if(XHOST.is_command("qdel")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qdel\")=TRUE",depth_long) << "[" << XHOST.command("qdel") << "]" << endl;} else {oss << "XHOST.is_command(\"qdel\")=FALSE" << endl;}
      if(XHOST.is_command("qhull")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qhull\")=TRUE",depth_long) << "[" << XHOST.command("qhull") << "]" << endl;} else {oss << "XHOST.is_command(\"qhull\")=FALSE" << endl;}
      if(XHOST.is_command("qstat")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qstat\")=TRUE",depth_long) << "[" << XHOST.command("qstat") << "]" << endl;} else {oss << "XHOST.is_command(\"qstat\")=FALSE" << endl;} //CO20200526
      if(XHOST.is_command("qsub")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"qsub\")=TRUE",depth_long) << "[" << XHOST.command("qsub") << "]" << endl;} else {oss << "XHOST.is_command(\"qsub\")=FALSE" << endl;}
      if(XHOST.is_command("rasmol")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"rasmol\")=TRUE",depth_long) << "[" << XHOST.command("rasmol") << "]" << endl;} else {oss << "XHOST.is_command(\"rasmol\")=FALSE" << endl;}
      if(XHOST.is_command("rsync")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"rsync\")=TRUE",depth_long) << "[" << XHOST.command("rsync") << "]" << endl;} else {oss << "XHOST.is_command(\"rsync\")=FALSE" << endl;}
      if(XHOST.is_command("sbatch")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"sbatch\")=TRUE",depth_long) << "[" << XHOST.command("sbatch") << "]" << endl;} else {oss << "XHOST.is_command(\"sbatch\")=FALSE" << endl;}
      if(XHOST.is_command("showq")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"showq\")=TRUE",depth_long) << "[" << XHOST.command("showq") << "]" << endl;} else {oss << "XHOST.is_command(\"showq\")=FALSE" << endl;} //CO20200526
      if(XHOST.is_command("sinfo")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"sinfo\")=TRUE",depth_long) << "[" << XHOST.command("sinfo") << "]" << endl;} else {oss << "XHOST.is_command(\"sinfo\")=FALSE" << endl;} //CO20200526
      if(XHOST.is_command("squeue")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"squeue\")=TRUE",depth_long) << "[" << XHOST.command("squeue") << "]" << endl;} else {oss << "XHOST.is_command(\"squeue\")=FALSE" << endl;} //CO20200526
      if(XHOST.is_command("scancel")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"scancel\")=TRUE",depth_long) << "[" << XHOST.command("scancel") << "]" << endl;} else {oss << "XHOST.is_command(\"scancel\")=FALSE" << endl;}
      if(XHOST.is_command("sensors")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"sensors\")=TRUE",depth_long) << "[" << XHOST.command("sensors") << "]" << endl;} else {oss << "XHOST.is_command(\"sensors\")=FALSE" << endl;}
      if(XHOST.is_command("unzip")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"unzip\")=TRUE",depth_long) << "[" << XHOST.command("unzip") << "]" << endl;} else {oss << "XHOST.is_command(\"unzip\")=FALSE" << endl;}
      if(XHOST.is_command("vasp46s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp46s\")=TRUE",depth_long) << "[" << XHOST.command("vasp46s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp46s\")=FALSE" << endl;}
      if(XHOST.is_command("vasp52s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp52s\")=TRUE",depth_long) << "[" << XHOST.command("vasp52s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp52s\")=FALSE" << endl;}
      if(XHOST.is_command("vasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp54s\")=TRUE",depth_long) << "[" << XHOST.command("vasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp54s\")=FALSE" << endl;}
      if(XHOST.is_command("vasp54s")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp54s\")=TRUE",depth_long) << "[" << XHOST.command("vasp54s") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp54s\")=FALSE" << endl;}
      if(XHOST.is_command("vasp_std")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"vasp_std\")=TRUE",depth_long) << "[" << XHOST.command("vasp_std") << "]" << endl;} else {oss << "XHOST.is_command(\"vasp_std\")=FALSE" << endl;} //CO20200624
      if(XHOST.is_command("wget")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"wget\")=TRUE",depth_long) << "[" << XHOST.command("wget") << "]" << endl;} else {oss << "XHOST.is_command(\"wget\")=FALSE" << endl;}
      if(XHOST.is_command("zip")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"zip\")=TRUE",depth_long) << "[" << XHOST.command("zip") << "]" << endl;} else {oss << "XHOST.is_command(\"zip\")=FALSE" << endl;}
      if(XHOST.is_command("xz")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"xz\")=TRUE",depth_long) << "[" << "xz" << "]" << endl;} else {oss << "XHOST.is_command(\"xz\")=FALSE" << endl;}
      if(XHOST.is_command("xxx")) {oss << aurostd::PaddedPOST("XHOST.is_command(\"xxx\")=TRUE",depth_long) << "[" << XHOST.command("xxx") << "]" << endl;} else {oss << "XHOST.is_command(\"xxx\")=FALSE" << endl;}
    }
    // #include "aflow_init_aus.cpp"
    // TEMPERATURE
    XHOST.sensors_allowed=TRUE;
    if(0)  if(XHOST.is_command("sensors")) {
      init::GetTEMP();
      if(INIT_VERBOSE) oss << "XHOST.vTemperatureCore.size()=" << XHOST.vTemperatureCore.size() << endl;
      if(INIT_VERBOSE && XHOST.vTemperatureCore.size()) {
        oss << "--- TEMPERATURES --------------- " << endl;
        for(uint i=0;i<XHOST.vTemperatureCore.size();i++) {oss << (i==0?"TEMP(C)=[":"") << XHOST.vTemperatureCore.at(i) << (i<XHOST.vTemperatureCore.size()-1?",":"]\n");}
        //;if(i<XHOST.vTemperatureCore.size()-1) oss << ","; else oss << "]" << endl; //CO20200106 - patching for auto-indenting
      }
    }
    // QUEUE STUFF
    XHOST.is_PBS=FALSE;
    XHOST.PBS_NUM_PPN=aurostd::getenv2uint("PBS_NUM_PPN");
    XHOST.PBS_NNODES=aurostd::getenv2uint("PBS_NNODES");
    if(XHOST.PBS_NUM_PPN!=0 || XHOST.PBS_NNODES!=0) XHOST.is_PBS=TRUE;
    XHOST.is_SLURM=FALSE;
    XHOST.SLURM_CPUS_ON_NODE=aurostd::getenv2int("SLURM_CPUS_ON_NODE");
    XHOST.SLURM_NNODES=aurostd::getenv2int("SLURM_NNODES");          
    XHOST.SLURM_NTASKS=aurostd::getenv2int("SLURM_NTASKS");          
    if(XHOST.SLURM_CPUS_ON_NODE!=0 || XHOST.SLURM_NNODES!=0 || XHOST.SLURM_NTASKS!=0) XHOST.is_SLURM=TRUE;
    if(INIT_VERBOSE) {
      oss << "--- QUEUES --------------- " << endl;
      oss << "is_PBS=" << XHOST.is_PBS << endl;
      oss << "PBS_NUM_PPN=" << XHOST.PBS_NUM_PPN << endl;
      oss << "PBS_NNODES=" << XHOST.PBS_NNODES << endl;
      oss << "is_SLURM=" << XHOST.is_SLURM << endl;
      oss << "SLURM_CPUS_ON_NODE=" << XHOST.SLURM_CPUS_ON_NODE << endl;
      oss << "SLURM_NNODES=" << XHOST.SLURM_NNODES << endl;
      oss << "SLURM_NTASKS=" << XHOST.SLURM_NTASKS << endl;
    }
    // maxmem
    XHOST.maxmem=aurostd::args2attachedutype<double>(XHOST.argv,"--mem=|--maxmem=",101.0);

    // MPIs
    if(INIT_VERBOSE) {
      oss << "--- DATES --------------- " << endl;
      oss << "XHOST.Date = " << XHOST.Date << endl;
      oss << "XHOST.maxmem = " << XHOST.maxmem << endl;
      oss << "--- MPI COMMANDS --- " << endl;
      oss << "MPI_OPTIONS_DUKE_BETA_MPICH=" << MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_BETA_MPICH=" << MPI_COMMAND_DUKE_BETA_MPICH << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_BETA_MPICH=" << MPI_BINARY_DIR_DUKE_BETA_MPICH << "\"" << endl;
      oss << "MPI_OPTIONS_DUKE_BETA_OPENMPI=" << MPI_OPTIONS_DUKE_BETA_OPENMPI << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_BETA_OPENMPI=" << MPI_COMMAND_DUKE_BETA_OPENMPI << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_BETA_OPENMPI=" << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << "\"" << endl;
      oss << "MPI_OPTIONS_DUKE_QRATS_MPICH=" << MPI_OPTIONS_DUKE_QRATS_MPICH << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_QRATS_MPICH=" << MPI_COMMAND_DUKE_QRATS_MPICH << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_QRATS_MPICH=" << MPI_BINARY_DIR_DUKE_QRATS_MPICH << "\"" << endl;
      //CO START
      oss << "MPI_OPTIONS_DUKE_QFLOW_OPENMPI=" << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_QFLOW_OPENMPI=" << MPI_COMMAND_DUKE_QFLOW_OPENMPI << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI=" << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << "\"" << endl;
      //CO20201220 X START
      oss << "MPI_OPTIONS_DUKE_X_X=" << MPI_OPTIONS_DUKE_X_X << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_X_X=" << MPI_COMMAND_DUKE_X_X << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_X_X=" << MPI_BINARY_DIR_DUKE_X_X << "\"" << endl;
      oss << "MPI_OPTIONS_DUKE_X_CRAY=" << MPI_OPTIONS_DUKE_X_CRAY << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_X_CRAY=" << MPI_COMMAND_DUKE_X_CRAY << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_X_CRAY=" << MPI_BINARY_DIR_DUKE_X_CRAY << "\"" << endl;
      oss << "MPI_OPTIONS_DUKE_X_OLDCRAY=" << MPI_OPTIONS_DUKE_X_OLDCRAY << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_X_OLDCRAY=" << MPI_COMMAND_DUKE_X_OLDCRAY << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_X_OLDCRAY=" << MPI_BINARY_DIR_DUKE_X_OLDCRAY << "\"" << endl;
      oss << "MPI_OPTIONS_DUKE_X_SMB=" << MPI_OPTIONS_DUKE_X_SMB << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_X_SMB=" << MPI_COMMAND_DUKE_X_SMB << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_X_SMB=" << MPI_BINARY_DIR_DUKE_X_SMB << "\"" << endl;
      //CO20201220 X STOP
      //CO20220818 JHU_ROCKFISH START
      oss << "MPI_OPTIONS_JHU_ROCKFISH=" << MPI_OPTIONS_JHU_ROCKFISH << "\"" << endl;
      oss << "MPI_COMMAND_JHU_ROCKFISH=" << MPI_COMMAND_JHU_ROCKFISH << "\"" << endl;
      oss << "MPI_BINARY_DIR_JHU_ROCKFISH=" << MPI_BINARY_DIR_JHU_ROCKFISH << "\"" << endl;
      //CO20220818 JHU_ROCKFISH STOP
      oss << "MPI_OPTIONS_MPCDF_EOS=" << MPI_OPTIONS_MPCDF_EOS << "\"" << endl;
      oss << "MPI_COMMAND_MPCDF_EOS=" << MPI_COMMAND_MPCDF_EOS << "\"" << endl;
      oss << "MPI_NCPUS_MPCDF_EOS=" << MPI_NCPUS_MPCDF_EOS << "\"" << endl;
      oss << "MPI_HYPERTHREADING_MPCDF_EOS=" << MPI_HYPERTHREADING_MPCDF_EOS << "\"" << endl;
      oss << "MPI_BINARY_DIR_MPCDF_EOS=" << MPI_BINARY_DIR_MPCDF_EOS << "\"" << endl;
      oss << "MPI_OPTIONS_MPCDF_DRACO=" << MPI_OPTIONS_MPCDF_DRACO << "\"" << endl;
      oss << "MPI_COMMAND_MPCDF_DRACO=" << MPI_COMMAND_MPCDF_DRACO << "\"" << endl;
      oss << "MPI_NCPUS_MPCDF_DRACO=" << MPI_NCPUS_MPCDF_DRACO << "\"" << endl;
      oss << "MPI_HYPERTHREADING_MPCDF_DRACO=" << MPI_HYPERTHREADING_MPCDF_DRACO << "\"" << endl;
      oss << "MPI_BINARY_DIR_MPCDF_DRACO=" << MPI_BINARY_DIR_MPCDF_DRACO << "\"" << endl;
      oss << "MPI_OPTIONS_MPCDF_COBRA=" << MPI_OPTIONS_MPCDF_COBRA << "\"" << endl;
      oss << "MPI_COMMAND_MPCDF_COBRA=" << MPI_COMMAND_MPCDF_COBRA << "\"" << endl;
      oss << "MPI_NCPUS_MPCDF_COBRA=" << MPI_NCPUS_MPCDF_COBRA << "\"" << endl;
      oss << "MPI_HYPERTHREADING_MPCDF_COBRA=" << MPI_HYPERTHREADING_MPCDF_COBRA << "\"" << endl;
      oss << "MPI_BINARY_DIR_MPCDF_COBRA=" << MPI_BINARY_DIR_MPCDF_COBRA << "\"" << endl;
      oss << "MPI_OPTIONS_MPCDF_HYDRA=" << MPI_OPTIONS_MPCDF_HYDRA << "\"" << endl;
      oss << "MPI_COMMAND_MPCDF_HYDRA=" << MPI_COMMAND_MPCDF_HYDRA << "\"" << endl;
      oss << "MPI_NCPUS_MPCDF_HYDRA=" << MPI_NCPUS_MPCDF_HYDRA << "\"" << endl;
      oss << "MPI_HYPERTHREADING_MPCDF_HYDRA=" << MPI_HYPERTHREADING_MPCDF_HYDRA << "\"" << endl;
      oss << "MPI_BINARY_DIR_MPCDF_HYDRA=" << MPI_BINARY_DIR_MPCDF_HYDRA << "\"" << endl;
      //CO END
      //DX20190509 - MACHINE001 - START
      oss << "MPI_OPTIONS_MACHINE001=" << MPI_OPTIONS_MACHINE001 << "\"" << endl;
      oss << "MPI_COMMAND_MACHINE001=" << MPI_COMMAND_MACHINE001 << "\"" << endl;
      oss << "MPI_BINARY_DIR_MACHINE001=" << MPI_BINARY_DIR_MACHINE001 << "\"" << endl;
      //DX20190509 - MACHINE001 - END
      //DX20190509 - MACHINE002 - START
      oss << "MPI_OPTIONS_MACHINE002=" << MPI_OPTIONS_MACHINE002 << "\"" << endl;
      oss << "MPI_COMMAND_MACHINE002=" << MPI_COMMAND_MACHINE002 << "\"" << endl;
      oss << "MPI_BINARY_DIR_MACHINE002=" << MPI_BINARY_DIR_MACHINE002 << "\"" << endl;
      //DX20190509 - MACHINE002 - END
      //DX20201005 - MACHINE003 - START
      oss << "MPI_OPTIONS_MACHINE003=" << MPI_OPTIONS_MACHINE003 << "\"" << endl;
      oss << "MPI_COMMAND_MACHINE003=" << MPI_COMMAND_MACHINE003 << "\"" << endl;
      oss << "MPI_BINARY_DIR_MACHINE003=" << MPI_BINARY_DIR_MACHINE003 << "\"" << endl;
      //DX20201005 - MACHINE003 - END
      //DX20190107 - CMU EULER - START
      oss << "MPI_OPTIONS_CMU_EULER=" << MPI_OPTIONS_CMU_EULER << "\"" << endl;
      oss << "MPI_COMMAND_CMU_EULER=" << MPI_COMMAND_CMU_EULER << "\"" << endl;
      oss << "MPI_BINARY_DIR_CMU_EULER=" << MPI_BINARY_DIR_CMU_EULER << "\"" << endl;
      //DX20190107 - CMU EULER - END
      oss << "MPI_OPTIONS_DUKE_MATERIALS=" << MPI_OPTIONS_DUKE_MATERIALS << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_MATERIALS=" << MPI_COMMAND_DUKE_MATERIALS << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_MATERIALS=" << MPI_BINARY_DIR_DUKE_MATERIALS << "\"" << endl;
      oss << "MPI_OPTIONS_DUKE_AFLOWLIB=" << MPI_OPTIONS_DUKE_AFLOWLIB << "\"" << endl;
      oss << "MPI_COMMAND_DUKE_AFLOWLIB=" << MPI_COMMAND_DUKE_AFLOWLIB << "\"" << endl;
      oss << "MPI_BINARY_DIR_DUKE_AFLOWLIB=" << MPI_BINARY_DIR_DUKE_AFLOWLIB << "\"" << endl;
      oss << "MPI_OPTIONS_FULTON_MARYLOU=" << MPI_OPTIONS_FULTON_MARYLOU << "\"" << endl;
      oss << "MPI_COMMAND_FULTON_MARYLOU=" << MPI_COMMAND_FULTON_MARYLOU << "\"" << endl;
      oss << "MPI_BINARY_DIR_FULTON_MARYLOU=" << MPI_BINARY_DIR_FULTON_MARYLOU << "\"" << endl;
      oss << "MPI_OPTIONS_MACHINE1=" << MPI_OPTIONS_MACHINE1 << "\"" << endl;
      oss << "MPI_COMMAND_MACHINE1=" << MPI_COMMAND_MACHINE1 << "\"" << endl;
      oss << "MPI_BINARY_DIR_MACHINE1=" << MPI_BINARY_DIR_MACHINE1 << "\"" << endl;
      oss << "MPI_OPTIONS_MACHINE2=" << MPI_OPTIONS_MACHINE2 << "\"" << endl;
      oss << "MPI_COMMAND_MACHINE2=" << MPI_COMMAND_MACHINE2 << "\"" << endl;
      oss << "MPI_BINARY_DIR_MACHINE2=" << MPI_BINARY_DIR_MACHINE2 << "\"" << endl;
      //   oss << endl;
    }
    // DO LISRS
    vector<string> vstrs;
    // DO VARIABLES
    aurostd::string2tokens(DEFAULT_VASP_POTCAR_DIRECTORIES,vVASP_POTCAR_DIRECTORIES,",");// vVASP_POTCAR_DIRECTORIES;
    // LIBRARIES
    aurostd::string2tokens(DEFAULT_AFLOW_LIBRARY_DIRECTORIES,vAFLOW_LIBRARY_DIRECTORIES,",");// vAFLOW_LIBRARY_DIRECTORIES;
    // PROJECTS
    vAFLOW_PROJECTS_DIRECTORIES.clear();
    // XHOST_LIBRARY_LIB2=LIBRARY_NOTHING;
    aurostd::string2tokens(DEFAULT_AFLOW_PROJECTS_DIRECTORIES,vstrs,",");// vAFLOW_PROJECTS_DIRECTORIES;
    // DEBUG  cerr << "vAFLOW_PROJECTS_DIRECTORIES.size()=" << vAFLOW_PROJECTS_DIRECTORIES.size() << endl; 
    for(uint i=0;i<vstrs.size();i++) { //  cerr << vstrs.at(i) << endl;
      if(aurostd::FileExist(vstrs.at(i))) {
        if(aurostd::substring2bool(vstrs.at(i),"AUID")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_AUID=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
        if(aurostd::FileExist(vstrs.at(i)+"/LIB")) {
          if(aurostd::substring2bool(vstrs.at(i),"ICSD")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_ICSD=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB0")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB0=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB1")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB1=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB2")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB2=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB3")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB3=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB4")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB4=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB5")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB5=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB6")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB6=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB7")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB7=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB8")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB8=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
          if(aurostd::substring2bool(vstrs.at(i),"LIB9")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB9=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}
        }else{  //CO20200624 - patch for LIB7 which has no LIB
          if(aurostd::substring2bool(vstrs.at(i),"LIB7")) {vAFLOW_PROJECTS_DIRECTORIES.push_back(vstrs.at(i)); XHOST_LIBRARY_LIB7=vAFLOW_PROJECTS_DIRECTORIES.size()-1;}  //CO20200624
        }
      }
    }

    // get position JSONL
    if(aurostd::EFileExist(DEFAULT_AFLOW_DB_DATA_PATH+"/aflow:00.jsonl") && aurostd::EFileExist(DEFAULT_AFLOW_DB_DATA_PATH+"/aflow:ff.jsonl")) {XHOST_LIBRARY_JSONL=DEFAULT_AFLOW_DB_DATA_PATH;} // check 00 and ff for being sure
    // DEBUG cerr << "vAFLOW_PROJECTS_DIRECTORIES.size()=" << vAFLOW_PROJECTS_DIRECTORIES.size() << endl;
    // DEBUG   for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++) cerr << "vAFLOW_PROJECTS_DIRECTORIES.at(i)=" << vAFLOW_PROJECTS_DIRECTORIES.at(i) << endl;

    if(INIT_VERBOSE) {
      oss << "--- PROJECTS @ " << XHOST.hostname << " --- " << endl;
      oss << "DEFAULT_AFLOW_DB_DATA_PATH=" << DEFAULT_AFLOW_DB_DATA_PATH << endl;
      oss << "XHOST_LIBRARY_JSONL=" << XHOST_LIBRARY_JSONL << endl;
      if(XHOST_LIBRARY_AUID!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_AUID=" << XHOST_LIBRARY_AUID << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_AUID << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_AUID) << endl;
      if(XHOST_LIBRARY_ICSD!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_ICSD=" << XHOST_LIBRARY_ICSD << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_ICSD << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD) << endl;
      if(XHOST_LIBRARY_LIB0!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB0=" << XHOST_LIBRARY_LIB0 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB0 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB0) << endl;
      if(XHOST_LIBRARY_LIB1!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB1=" << XHOST_LIBRARY_LIB1 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB1 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1) << endl;
      if(XHOST_LIBRARY_LIB2!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB2=" << XHOST_LIBRARY_LIB2 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB2 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2) << endl;
      if(XHOST_LIBRARY_LIB3!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB3=" << XHOST_LIBRARY_LIB3 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB3 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3) << endl;
      if(XHOST_LIBRARY_LIB4!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB4=" << XHOST_LIBRARY_LIB4 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB4 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4) << endl;
      if(XHOST_LIBRARY_LIB5!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB5=" << XHOST_LIBRARY_LIB5 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB5 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5) << endl;
      if(XHOST_LIBRARY_LIB6!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB6=" << XHOST_LIBRARY_LIB6 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB6 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6) << endl;
      if(XHOST_LIBRARY_LIB7!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB7=" << XHOST_LIBRARY_LIB7 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB7 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7) << endl;
      if(XHOST_LIBRARY_LIB8!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB8=" << XHOST_LIBRARY_LIB8 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB8 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8) << endl;
      if(XHOST_LIBRARY_LIB9!=LIBRARY_NOTHING) oss << "XHOST_LIBRARY_LIB9=" << XHOST_LIBRARY_LIB9 << " : vAFLOW_PROJECTS_DIRECTORIES.at(" << XHOST_LIBRARY_LIB9 << ")=" << vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9) << endl;

      init::InitLoadString("vLIBS");
      if(XHOST_LIBRARY_ICSD!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_ICSD_LIB.size()=" << XHOST_Library_CALCULATED_ICSD_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB0!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB0_LIB.size()=" << XHOST_Library_CALCULATED_LIB0_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB1!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB1_LIB.size()=" << XHOST_Library_CALCULATED_LIB1_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB2!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB2_LIB.size()=" << XHOST_Library_CALCULATED_LIB2_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB3!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB3_LIB.size()=" << XHOST_Library_CALCULATED_LIB3_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB4!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB4_LIB.size()=" << XHOST_Library_CALCULATED_LIB4_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB5!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB5_LIB.size()=" << XHOST_Library_CALCULATED_LIB5_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB6!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB6_LIB.size()=" << XHOST_Library_CALCULATED_LIB6_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB7!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB7_LIB.size()=" << XHOST_Library_CALCULATED_LIB7_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB8!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB8_LIB.size()=" << XHOST_Library_CALCULATED_LIB8_LIB.size() << endl; }
      if(XHOST_LIBRARY_LIB9!=LIBRARY_NOTHING) { oss << "Library_CALCULATED_LIB9_LIB.size()=" << XHOST_Library_CALCULATED_LIB9_LIB.size() << endl; }
    }

    // OLD aurostd::string2tokens(string(AFLOW_PROJECTS_DIRECTORIES),vAFLOW_PROJECTS_DIRECTORIES,",");// vAFLOW_PROJECTS_DIRECTORIES;

    // check for MACHINES MARYLOU
    XHOST.is_MACHINE_FULTON_MARYLOU=FALSE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(XHOST.user=="fslcollab8" || XHOST.user=="glh43" || XHOST.user=="legoses") XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(aurostd::substring2bool(XHOST.group,"fslcollab")) XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(aurostd::substring2bool(XHOST.group,"fslg")) XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    if(!XHOST.is_MACHINE_FULTON_MARYLOU) if(aurostd::substring2bool(XHOST.group,"glh43")) XHOST.is_MACHINE_FULTON_MARYLOU=TRUE;
    // some other technique to get MARYLOU

    // check for APENNSY_USE_SERVER/AFLOWLIB
    if(XHOST.hostname=="nietzsche.mems.duke.edu" || XHOST.hostname=="materials.duke.edu" || XHOST.hostname=="aflowlib.duke.edu") {
      XHOST.APENNSY_USE_SERVER=TRUE;XHOST.APENNSY_USE_LIBRARY=FALSE;XHOST.APENNSY_SERVER_AFLOWLIB_ORG=FALSE;
    } else {
      XHOST.APENNSY_USE_SERVER=FALSE;XHOST.APENNSY_USE_LIBRARY=FALSE;XHOST.APENNSY_SERVER_AFLOWLIB_ORG=TRUE;
    }
    if(INIT_VERBOSE) {
      oss << "--- APENNSY_USE_*** @ " << XHOST.hostname << " --- " << endl;
      oss << "XHOST.APENNSY_USE_SERVER=" << XHOST.APENNSY_USE_SERVER << endl;
      oss << "XHOST.APENNSY_USE_LIBRARY=" << XHOST.APENNSY_USE_LIBRARY << endl;
      oss << "XHOST.APENNSY_SERVER_AFLOWLIB_ORG=" << XHOST.APENNSY_SERVER_AFLOWLIB_ORG << endl;
    }    
    // DO aflow.in and LOCK
    if(INIT_VERBOSE) oss << "--- LOADING @ _AFLOWIN_ and _AFLOWLOCK_ --- " << endl;
    _AFLOWIN_=aurostd::args2attachedstring(XHOST.argv,"--use_aflow.in=",_AFLOWIN_DEFAULT_);
    _AFLOWLOCK_=aurostd::args2attachedstring(XHOST.argv,"--use_LOCK=",_AFLOWLOCK_DEFAULT_);
    _STOPFLOW_=aurostd::args2attachedstring(XHOST.argv,"--use_stop_file=","STOPFLOW");  //CO20210315
    if(INIT_VERBOSE) oss << "_AFLOWIN_=" << _AFLOWIN_ << endl;
    if(INIT_VERBOSE) oss << "_AFLOWLOCK_=" << _AFLOWLOCK_ << endl;
    if(INIT_VERBOSE) oss << "_STOPFLOW_=" << _STOPFLOW_ << endl;  //CO20210315

    // LOAD control stuff
    if(INIT_VERBOSE) oss << "--- LOADING @ control --- " << endl;
    XHOST.vflag_control.flag("MACHINE",aurostd::args2flag(argv,cmds,"--machine|--machine="));
    XHOST.vflag_control.flag("MULTI=SH",aurostd::args2flag(argv,cmds,"--multi=sh"));
    XHOST.vflag_control.flag("MULTI=BZIP2",aurostd::args2flag(argv,cmds,"--multi=bzip2"));
    XHOST.vflag_control.flag("MULTI=BUNZIP2",aurostd::args2flag(argv,cmds,"--multi=bunzip2"));
    XHOST.vflag_control.flag("MULTI=GZIP",aurostd::args2flag(argv,cmds,"--multi=gzip"));
    XHOST.vflag_control.flag("MULTI=GUNZIP",aurostd::args2flag(argv,cmds,"--multi=gunzip"));
    XHOST.vflag_control.flag("MULTI=XZIP",aurostd::args2flag(argv,cmds,"--multi=xz|--multi=xzip"));
    XHOST.vflag_control.flag("MULTI=XUNZIP",aurostd::args2flag(argv,cmds,"--multi=xunzip"));
    XHOST.vflag_control.flag("MULTI=BZ2XZ",aurostd::args2flag(argv,cmds,"--multi=bz2xz"));
    XHOST.vflag_control.flag("MULTI=GZ2XZ",aurostd::args2flag(argv,cmds,"--multi=gz2xz"));
    XHOST.vflag_control.flag("MULTI=ZIP",aurostd::args2flag(argv,cmds,"--multi=zip"));
    XHOST.vflag_control.flag("MONITOR",aurostd::args2flag(argv,cmds,"--monitor"));
    XHOST.vflag_control.flag("MONITOR_VASP",aurostd::args2flag(argv,cmds,"--monitor_vasp"));
    XHOST.vflag_control.flag("AFLOW_STARTUP_SCRIPT",aurostd::args2attachedflag(argv,cmds,"--aflow_startup_script=")); //SD20220502
    XHOST.vflag_control.push_attached("AFLOW_STARTUP_SCRIPT",aurostd::args2attachedstring(argv,"--aflow_startup_script=","")); //SD20220502
    //[SD20220402 - OBSOLETE]XHOST.vflag_control.flag("KILL_VASP_ALL",aurostd::args2flag(argv,cmds,"--kill_vasp_all|--killvaspall"));  //CO20210315 - issue non-specific killall vasp command
    XHOST.vflag_control.flag("KILL_VASP_OOM",aurostd::args2flag(argv,cmds,"--kill_vasp_oom|--killvaspoom"));  //CO20210315 - kill vasp if approaching OOM
    XHOST.vflag_control.flag("GETTEMP",aurostd::args2flag(argv,cmds,"--getTEMP|--getTEMPS|--getTEMPs|--gettemp|--gettemps"));
    XHOST.vflag_control.flag("SWITCH_AFLOW",
        aurostd::args2flag(argv,cmds,"--run|--clean|--xclean|--multi|--generate") ||
        aurostd::args2attachedflag(argv,cmds,"--run=") ||
        aurostd::args2flag(argv,cmds,"--generate_vasp_from_aflowin|--generate_aflowin_from_vasp"));
    //DX START
    XHOST.vflag_control.flag("AFLOWIN_SYM",aurostd::args2flag(argv,cmds,"--generate_symmetry|--generate_sym")); //DX
    //DX END
    // [OBSOLETE]    XHOST.vflag_control.flag("SWITCH_AFLOWLIB",aurostd::args2flag(argv,cmds,"--aflowlib") || aurostd::args2attachedflag(argv,cmds,"--aflowlib="));
    XHOST.vflag_control.flag("SWITCH_APENNSY1",aurostd::args2flag(argv,cmds,"--apennsy|--lib2|--lib2u|--lib2pgm|--LIB2|--LIB2U|--LIB2PGM|--libraryX|--libraryU|--libraryPGM|--alloy"));
    XHOST.vflag_control.flag("SWITCH_APENNSY2",aurostd::args2flag(argv,cmds,"--apool|--apool_private|--apool_test|--library2|-simpls|--simpls|--VASPIN|--energy|--psenergy"));

    // DIRECTORY NEEDS A SPECIAL TREATMENT
    found=FALSE;

    string dir_default="./",dir=dir_default;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--DIRECTORY|--directory|--D|--d"))) {dir=aurostd::args2string(argv,"--DIRECTORY|--directory|--D|--d",dir_default);if(INIT_VERBOSE) cerr << "--DIRECTORY " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"./"))) {dir="./";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"."))) {dir=".";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"../"))) {dir="../";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"~/"))) {dir="~/";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2flag(argv,"~"))) {dir="~";if(INIT_VERBOSE) cerr << "taking: " << dir << " " << endl;}
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--DIRECTORY=|--D=|--d="))) {dir=aurostd::args2attachedstring(argv,"--DIRECTORY=|--D=|--d=",dir_default);if(INIT_VERBOSE) cerr << "--DIRECTORY=" << dir << " " << endl;}
    XHOST.vflag_control.flag("DIRECTORY",found);  // if found
    if(XHOST.vflag_control.flag("DIRECTORY")) XHOST.vflag_control.push_attached("DIRECTORY",dir);
    if(XHOST.vflag_control.flag("DIRECTORY")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"DIR\")=[" << XHOST.vflag_control.getattachedscheme("DIRECTORY") << "]" << endl; 
    //CO20190402 START - cleaning directory, giving us something to print with logger
    string directory_clean="";
    if(XHOST.vflag_control.flag("DIRECTORY")) {directory_clean=XHOST.vflag_control.getattachedscheme("DIRECTORY");}
    aurostd::StringSubst(directory_clean,"//","/");  //clean
    directory_clean=aurostd::RemoveWhiteSpaces(directory_clean);
    if(directory_clean.empty() || directory_clean=="./" || directory_clean==".") {directory_clean=aurostd::getPWD()+"/";}  //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
    if(!directory_clean.empty()) {XHOST.vflag_control.flag("DIRECTORY_CLEAN",TRUE);XHOST.vflag_control.push_attached("DIRECTORY_CLEAN",directory_clean);}
    if(LDEBUG) {cerr << directory_clean << endl;}
    //CO20190402 STOP - cleaning directory, giving us something to print with logger

    // LOGICS to intercept --file= => XHOST.vflag_control.flag("FILE") XHOST.vflag_control.getattachedscheme("FILE")
    // FILE NEEDS A SPECIAL TREATMENT
    string file_default="xxxx",file=file_default;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--FILE|--file|--F|--f"))) {file=aurostd::args2string(argv,"--FILE|--file|--F|--f",file_default);if(INIT_VERBOSE) cerr << "--FILE " << file << " " << endl;}
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--FILE=|--file=|--F=|--f="))) {file=aurostd::args2attachedstring(argv,"--FILE=|--file=|--F=|--f=",file_default);if(INIT_VERBOSE) cerr << "--FILE=" << file << " " << endl;}
    XHOST.vflag_control.flag("FILE",found);  // if found
    if(XHOST.vflag_control.flag("FILE")) XHOST.vflag_control.push_attached("FILE",file);
    if(XHOST.vflag_control.flag("FILE")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"FILE\")=[" << XHOST.vflag_control.getattachedscheme("FILE") << "]" << endl;


    // VFILES NEEDS A SPECIAL TREATMENT
    string dirs_default="xxxx",dirs=dirs_default;
    vector<string> vdir;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--DIRECTORY|--directory|--D|--d"))) {dirs="";vdir=aurostd::args2vectorstring(argv,"--DIRECTORY|--directory|--D|--d",dirs_default);if(INIT_VERBOSE) cerr << "--DIRECTORY " << vdir.size() << " " << endl;}
    if(found) for(uint i=0;i<vdir.size();i++) dirs+=vdir.at(i)+(i<vdir.size()-1?",":""); // glue
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--DIRECTORY=|--D=|--d="))) {dirs=aurostd::args2attachedstring(argv,"--DIRECTORY=|--D=|--d=",dirs_default);if(INIT_VERBOSE) cerr << "--DIRECTORY=" << dirs << " " << endl;}
    XHOST.vflag_control.flag("VDIR",found);  // if found
    if(XHOST.vflag_control.flag("VDIR")) XHOST.vflag_control.push_attached("VDIR",dirs); 
    if(XHOST.vflag_control.flag("VDIR")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"VDIR\")=[" << XHOST.vflag_control.getattachedscheme("VDIR") << "]" << endl; 
    // VFILES NEEDS A SPECIAL TREATMENT
    string files_default="xxxx",files=files_default;
    vector<string> vfile;
    found=FALSE;
    if(!found&&(found=aurostd::args2flag(argv,"--FILE|--file|--F|--f"))) {files="";vfile=aurostd::args2vectorstring(argv,"--FILE|--file|--F|--f",files_default);if(INIT_VERBOSE) cerr << "--FILE " << vfile.size() << " " << endl;}
    if(found) for(uint i=0;i<vfile.size();i++) files+=vfile.at(i)+(i<vfile.size()-1?",":""); // glue
    if(!found&&(found=aurostd::args2attachedflag(argv,cmds,"--FILE=|--file=|--F=|--f="))) {files=aurostd::args2attachedstring(argv,"--FILE=|--file=|--F=|--f=",files_default);if(INIT_VERBOSE) cerr << "--FILE=" << files << " " << endl;}
    XHOST.vflag_control.flag("VFILES",found);  // if found
    if(XHOST.vflag_control.flag("VFILES")) XHOST.vflag_control.push_attached("VFILES",files); 
    if(XHOST.vflag_control.flag("VFILES")) if(INIT_VERBOSE) cerr << "XHOST.vflag_control.flag(\"VFILES\")=[" << XHOST.vflag_control.getattachedscheme("VFILES") << "]" << endl; 

    XHOST.vflag_control.flag("AFLOW_HELP",aurostd::args2flag(argv,cmds,"-h|--help"));
    XHOST.vflag_control.flag("AFLOW_EXCEPTIONS", aurostd::args2flag(argv, cmds, "-e|--errors|--exceptions"));  //ME20180531
    XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3",aurostd::args2flag(argv,cmds,"-l|--license"));
    XHOST.vflag_control.flag("README_AFLOW",aurostd::args2flag(argv,cmds,"--readme=run|--readme=aflow|--readme_aflow"));
    XHOST.vflag_control.flag("README_AFLOW_PFLOW",aurostd::args2flag(argv,cmds,"--readme=pflow|--readme=processor|--readme=aconvasp|--readme_aconvasp"));
    XHOST.vflag_control.flag("README_FROZSL",aurostd::args2flag(argv,cmds,"--readme=frozsl|--readme_frozsl"));
    XHOST.vflag_control.flag("README_APL",aurostd::args2flag(argv,cmds,"--readme=apl|--readme_apl"));
    XHOST.vflag_control.flag("README_QHA",aurostd::args2flag(argv,cmds,"--readme=qha|--readme_qha|--readme=qha3p|--readme_qha3p|--readme=scqha|--readme_scqha"));
    XHOST.vflag_control.flag("README_AAPL",aurostd::args2flag(argv,cmds,"--readme=aapl|--readme_aapl"));
    XHOST.vflag_control.flag("README_AGL",aurostd::args2flag(argv,cmds,"--readme=agl|--readme_agl"));
    XHOST.vflag_control.flag("README_AEL",aurostd::args2flag(argv,cmds,"--readme=ael|--readme_ael"));
    XHOST.vflag_control.flag("README_ANRL",aurostd::args2flag(argv,cmds,"--readme=prototypes|--readme_prototypes|--readme=anrl|--readme_anrl")); //DX20210422 - added prototypes variants
    XHOST.vflag_control.flag("README_COMPARE",aurostd::args2flag(argv,cmds,"--readme=xtalfinder|--readme_xtalfinder|--readme=compare|--readme_compare")); //DX20210422 - added xtalfinder variants
    XHOST.vflag_control.flag("README_GFA",aurostd::args2flag(argv,cmds,"--readme=gfa|--readme_gfa")); //CO20190401
    XHOST.vflag_control.flag("README_SYMMETRY",aurostd::args2flag(argv,cmds,"--readme=symmetry|--readme_symmetry"));
    XHOST.vflag_control.flag("README_CCE",aurostd::args2flag(argv,cmds,"--readme=cce|--readme_cce")); //CO20190620
    XHOST.vflag_control.flag("README_CHULL",aurostd::args2flag(argv,cmds,"--readme=chull|--readme_chull")); //CO20190620
    XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION",aurostd::args2flag(argv,cmds,"--readme=partial_occupation|--readme=pocc|--readme_pocc"));
    XHOST.vflag_control.flag("README_APENNSY",aurostd::args2flag(argv,cmds,"--readme=apennsy|--readme_apennsy"));
    XHOST.vflag_control.flag("README_SCRIPTING",aurostd::args2flag(argv,cmds,"--readme=scripting|--readme_scripting|--readme=script|--readme_script"));
    XHOST.vflag_control.flag("README_EXCEPTIONS", aurostd::args2flag(argv, cmds, "--readme=errors|--readme_errors|--readme=exceptions|--readme_exceptions"));  //ME20180531
    XHOST.vflag_control.flag("README_HTRESOURCES",aurostd::args2flag(argv,cmds,"--readme=htresources|--readme_htresources|--readme=resources|--readme_resources"));
    XHOST.vflag_control.flag("README_XAFLOW",aurostd::args2flag(argv,cmds,"--readme=xaflow|--readme_xaflow"));
    XHOST.vflag_control.flag("README_AFLOWRC",aurostd::args2flag(argv,cmds,"--readme=aflowrc|--readme_aflowrc"));
    if(!(XHOST.vflag_control.flag("AFLOW_HELP") || 
          XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3") ||
          XHOST.vflag_control.flag("README_AFLOW") ||
          XHOST.vflag_control.flag("README_AFLOW_PFLOW") ||
          XHOST.vflag_control.flag("README_FROZSL") ||
          XHOST.vflag_control.flag("README_APL") ||
          XHOST.vflag_control.flag("README_QHA") ||
          XHOST.vflag_control.flag("README_AAPL") ||
          XHOST.vflag_control.flag("README_AGL") ||
          XHOST.vflag_control.flag("README_AEL") ||
          XHOST.vflag_control.flag("README_ANRL") ||
          XHOST.vflag_control.flag("README_COMPARE") ||
          XHOST.vflag_control.flag("README_GFA") || //CO20190401
          XHOST.vflag_control.flag("README_SYMMETRY") ||
          XHOST.vflag_control.flag("README_CCE") || //CO20190620
          XHOST.vflag_control.flag("README_CHULL") || //CO20190620
          XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION") ||
          XHOST.vflag_control.flag("README_APENNSY") ||
          XHOST.vflag_control.flag("README_SCRIPTING") ||
          XHOST.vflag_control.flag("README_EXCEPTIONS") ||  //ME20180531
          XHOST.vflag_control.flag("README_HTRESOURCES") ||
          XHOST.vflag_control.flag("README_XAFLOW") ||
          XHOST.vflag_control.flag("README_AFLOWRC") ||
          FALSE)){  //CO20180306 - need to catch --readme such that it doesn't interfere with --readme=
            XHOST.vflag_control.flag("AFLOW_HELP",aurostd::args2flag(argv,cmds,"--readme"));  //CO20180306
          }
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AFLOW_HELP\")=" << XHOST.vflag_control.flag("AFLOW_HELP") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOW_LICENSE_GPL3\")=" << XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOW\")=" << XHOST.vflag_control.flag("README_AFLOW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOW_PFLOW\")=" << XHOST.vflag_control.flag("README_AFLOW_PFLOW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_FROZSL\")=" << XHOST.vflag_control.flag("README_FROZSL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_APL\")=" << XHOST.vflag_control.flag("README_APL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AGL\")=" << XHOST.vflag_control.flag("README_AGL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AEL\")=" << XHOST.vflag_control.flag("README_AEL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_ANRL\")=" << XHOST.vflag_control.flag("README_ANRL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_COMPARE\")=" << XHOST.vflag_control.flag("README_COMPARE") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_GFA\")=" << XHOST.vflag_control.flag("README_GFA") << endl;  //CO20190401
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_SYMMETRY\")=" << XHOST.vflag_control.flag("README_SYMMETRY") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_CCE\")=" << XHOST.vflag_control.flag("README_CCE") << endl;  //CO20190620
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_CHULL\")=" << XHOST.vflag_control.flag("README_CHULL") << endl;  //CO20190620
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_PARTIAL_OCCUPATION\")=" << XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_APENNSY\")=" << XHOST.vflag_control.flag("README_APENNSY") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_SCRIPTING\")=" << XHOST.vflag_control.flag("README_SCRIPTING") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_HTRESOURCES\")=" << XHOST.vflag_control.flag("README_HTRESOURCES") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_XAFLOW\")=" << XHOST.vflag_control.flag("README_XAFLOW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"README_AFLOWRC\")=" << XHOST.vflag_control.flag("README_AFLOWRC") << endl;

    // arguments
    string keep=aurostd::args2attachedstring(XHOST.argv,"--keep=","");
    XHOST.vflag_control.flag("KEEP::TEX",aurostd::substring2bool(keep,"tex") || aurostd::substring2bool(keep,"TEX"));
    XHOST.vflag_control.flag("KEEP::DVI",aurostd::substring2bool(keep,"dvi") || aurostd::substring2bool(keep,"DVI"));
    XHOST.vflag_control.flag("KEEP::TOC",aurostd::substring2bool(keep,"toc") || aurostd::substring2bool(keep,"TOC"));
    XHOST.vflag_control.flag("KEEP::EPS",aurostd::substring2bool(keep,"eps") || aurostd::substring2bool(keep,"EPS"));
    XHOST.vflag_control.flag("KEEP::PDF",aurostd::substring2bool(keep,"pdf") || aurostd::substring2bool(keep,"PDF"));
    XHOST.vflag_control.flag("KEEP::JPG",aurostd::substring2bool(keep,"jpg") || aurostd::substring2bool(keep,"JPG"));
    XHOST.vflag_control.flag("KEEP::PNG",aurostd::substring2bool(keep,"png") || aurostd::substring2bool(keep,"PNG"));
    XHOST.vflag_control.flag("KEEP::GIF",aurostd::substring2bool(keep,"gif") || aurostd::substring2bool(keep,"GIF"));
    XHOST.vflag_control.flag("KEEP::GPL",aurostd::substring2bool(keep,"gpl") || aurostd::substring2bool(keep,"GPL") || aurostd::substring2bool(keep,"gnuplot"));
    XHOST.vflag_control.flag("KEEP::MAT",aurostd::substring2bool(keep,"mat") || aurostd::substring2bool(keep,"matlab"));
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::TEX\")=" << XHOST.vflag_control.flag("KEEP::TEX") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::DVI\")=" << XHOST.vflag_control.flag("KEEP::DVI") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::TOC\")=" << XHOST.vflag_control.flag("KEEP::TOC") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::EPS\")=" << XHOST.vflag_control.flag("KEEP::EPS") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::PDF\")=" << XHOST.vflag_control.flag("KEEP::PDF") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::GPL\")=" << XHOST.vflag_control.flag("KEEP::GPL") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::MAT\")=" << XHOST.vflag_control.flag("KEEP::MAT") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::JPG\")=" << XHOST.vflag_control.flag("KEEP::JPG") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::PNG\")=" << XHOST.vflag_control.flag("KEEP::PNG") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"KEEP::GIF\")=" << XHOST.vflag_control.flag("KEEP::GIF") << endl;

    XHOST.vflag_control.flag("PRINT_MODE::HTML",aurostd::args2flag(XHOST.argv,cmds,"--print=html|--print_html")); 
    XHOST.vflag_control.flag("PRINT_MODE::TXT",aurostd::args2flag(XHOST.argv,cmds,"--print=txt|--print_txt")); 
    XHOST.vflag_control.flag("PRINT_MODE::JSON",aurostd::args2flag(XHOST.argv,cmds,"--print=json|--print_json")); //DX20170907 - Add json
    XHOST.vflag_control.flag("PRINT_MODE::PYTHON",aurostd::args2flag(XHOST.argv,cmds,"--print=python|--print_python")); //DX20201228 - add Python
    XHOST.vflag_control.flag("PRINT_MODE::LATEX",aurostd::args2flag(XHOST.argv,cmds,"--print=latex|--print_latex"));
    XHOST.vflag_control.flag("PRINT_MODE::YEAR",aurostd::args2flag(XHOST.argv,cmds,"--print=year|--print_year"));
    XHOST.vflag_control.flag("PRINT_MODE::DOI",aurostd::args2flag(XHOST.argv,cmds,"--print=doi|--print_doi"));
    XHOST.vflag_control.flag("PRINT_MODE::BIBTEX",aurostd::args2flag(XHOST.argv,cmds,"--print=bibtex|--print_bibtex"));
    XHOST.vflag_control.flag("PRINT_MODE::EXTRA",aurostd::args2flag(argv,cmds,"--print=extra|--print=vextra_html|--print_vextra_html"));
    XHOST.vflag_control.flag("PRINT_MODE::NUMBER",aurostd::args2flag(XHOST.argv,cmds,"--print=number|--print_number"));
    XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS",aurostd::args2flag(XHOST.argv,cmds,"--print=hyperlinks|--print_hyperlinks|--print=hyperlink|--print_hyperlink"));
    XHOST.vflag_control.flag("PRINT_MODE::NOTE",aurostd::args2flag(XHOST.argv,cmds,"--print=note|--print_note|--print=notes|--print_notes"));
    XHOST.vflag_control.flag("PRINT_MODE::NEW",aurostd::args2flag(XHOST.argv,cmds,"--print=new|--print_new"));
    XHOST.vflag_control.flag("PRINT_MODE::DATE",aurostd::args2flag(XHOST.argv,cmds,"--print=date|--print_date"));
    XHOST.vflag_control.flag("APENNSY::LATEX_SNAPSHOT",aurostd::args2flag(argv,cmds,"--snapshot"));
    XHOST.vflag_control.flag("APENNSY::LATEX_OUTPUT",!aurostd::args2flag(argv,cmds,"--NOLATEX|--nolatex"));
    XHOST.vflag_control.flag("APENNSY::LATEX_CITE",aurostd::args2flag(argv,cmds,"--cite"));

    XHOST.vflag_control.flag("OSS::COUT",aurostd::args2flag(XHOST.argv,cmds,"--oss=cout|--oss_cout|--COUT|--cout"));
    XHOST.vflag_control.flag("OSS::CERR",aurostd::args2flag(XHOST.argv,cmds,"--oss=cerr|--oss_cerr|--CERR|--cerr"));
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::TXT\")=" << XHOST.vflag_control.flag("PRINT_MODE::TXT") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::JSON\")=" << XHOST.vflag_control.flag("PRINT_MODE::JSON") << endl; //DX20170907 - Add json
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::PYTHON\")=" << XHOST.vflag_control.flag("PRINT_MODE::PYTHON") << endl; //DX20170907 - add Python
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::LATEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::LATEX") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::YEAR\")=" << XHOST.vflag_control.flag("PRINT_MODE::YEAR") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DOI\")=" << XHOST.vflag_control.flag("PRINT_MODE::DOI") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::BIBTEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EXTRA\")=" << XHOST.vflag_control.flag("PRINT_MODE::EXTRA") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NUMBER\")=" << XHOST.vflag_control.flag("PRINT_MODE::NUMBER") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NOTE\")=" << XHOST.vflag_control.flag("PRINT_MODE::NOTE") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NEW\")=" << XHOST.vflag_control.flag("PRINT_MODE::NEW") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DATE\")=" << XHOST.vflag_control.flag("PRINT_MODE::DATE") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"OSS::COUT\")=" << XHOST.vflag_control.flag("OSS::COUT") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"OSS::CERR\")=" << XHOST.vflag_control.flag("OSS::CERR") << endl;
    //  INIT_VERBOSE=FALSE;

    XHOST.vflag_control.flag("BEEP",aurostd::args2flag(XHOST.argv,cmds,"--beep")); 
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"BEEP\")=" << XHOST.vflag_control.flag("BEEP") << endl;
    XHOST.vflag_control.flag("REMOVE",aurostd::args2flag(XHOST.argv,cmds,"--remove")); 
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"REMOVE\")=" << XHOST.vflag_control.flag("REMOVE") << endl;
    XHOST.vflag_control.flag("ZIP",aurostd::args2flag(XHOST.argv,cmds,"--zip")); 
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"ZIP\")=" << XHOST.vflag_control.flag("ZIP") << endl;

    //  XHOST.vflag_control.flag("PRINT_MODE::EPS",aurostd::args2flag(XHOST.argv,cmds,"--print=eps|--print=eps"));
    XHOST.vflag_control.flag("PRINT_MODE::EPS",TRUE); // default
    XHOST.vflag_control.flag("PRINT_MODE::PDF",aurostd::args2flag(XHOST.argv,cmds,"--print=pdf|--print_pdf")); if(XHOST.vflag_control.flag("PRINT_MODE::PDF")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::GIF",aurostd::args2flag(XHOST.argv,cmds,"--print=gif|--print_gif")); if(XHOST.vflag_control.flag("PRINT_MODE::GIF")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::JPG",aurostd::args2flag(XHOST.argv,cmds,"--print=jpg|--print_jpg")); if(XHOST.vflag_control.flag("PRINT_MODE::JPG")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    XHOST.vflag_control.flag("PRINT_MODE::PNG",aurostd::args2flag(XHOST.argv,cmds,"--print=png|--print_png")); if(XHOST.vflag_control.flag("PRINT_MODE::PNG")) XHOST.vflag_control.flag("PRINT_MODE::EPS",FALSE);
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EPS\")=" << XHOST.vflag_control.flag("PRINT_MODE::EPS") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::PDF\")=" << XHOST.vflag_control.flag("PRINT_MODE::PDF") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::GIF\")=" << XHOST.vflag_control.flag("PRINT_MODE::GIF") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::JPG\")=" << XHOST.vflag_control.flag("PRINT_MODE::JPG") << endl;
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::PNG\")=" << XHOST.vflag_control.flag("PRINT_MODE::PNG") << endl;

    // Allow the preselection of the EntryLoader source (will automatically fall back to the next best source if given source is not accessible)
    XHOST.vflag_control.flag("ENTRY_LOADER::SOURCE",aurostd::args2attachedflag(argv,"--entry_loader_source=")); //HE20220826
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"ENTRY_LOADER::SOURCE\")=" << XHOST.vflag_control.flag("ENTRY_LOADER::SOURCE") << endl;  //HE20220826
    if(XHOST.vflag_control.flag("ENTRY_LOADER::SOURCE")) XHOST.vflag_control.push_attached("ENTRY_LOADER::SOURCE",aurostd::args2attachedstring(argv,"--entry_loader_source=","")); //HE20220826
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"ENTRY_LOADER::SOURCE\")=" << XHOST.vflag_control.getattachedscheme("ENTRY_LOADER::SOURCE") << endl;  //HE20220826


    //[CO20191110]run pocc post-processing for particular temperatures from command line
    XHOST.vflag_control.flag("CALCULATION_TEMPERATURE",aurostd::args2attachedflag(argv,"--temperature=|--temp="));  //CO20191110
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"CALCULATION_TEMPERATURE\")=" << XHOST.vflag_control.flag("CALCULATION_TEMPERATURE") << endl;  //CO20191110
    if(XHOST.vflag_control.flag("CALCULATION_TEMPERATURE")) XHOST.vflag_control.push_attached("CALCULATION_TEMPERATURE",aurostd::args2attachedstring(argv,"--temperature=|--temp=","300")); //CO20191110
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"CALCULATION_TEMPERATURE\")=" << XHOST.vflag_control.getattachedscheme("CALCULATION_TEMPERATURE") << endl;  //CO20191110
    //[CO20200624]run pocc post-processing to skip bad aruns
    XHOST.vflag_control.flag("ARUNS2SKIP",aurostd::args2attachedflag(argv,"--aruns2skip=|--arun2skip="));  //CO20200624
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"ARUNS2SKIP\")=" << XHOST.vflag_control.flag("ARUNS2SKIP") << endl;  //CO20200624
    if(XHOST.vflag_control.flag("ARUNS2SKIP")) XHOST.vflag_control.push_attached("ARUNS2SKIP",aurostd::args2attachedstring(argv,"--aruns2skip=|--arun2skip=","")); //CO20200624
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"ARUNS2SKIP\")=" << XHOST.vflag_control.getattachedscheme("ARUNS2SKIP") << endl;  //CO20200624
    XHOST.vflag_control.flag("NEGLECT_CCE",aurostd::args2attachedflag(argv,"--neglect_cce|--no_cce|--neglectcce|--nocce")); //CO20210115
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"NEGLECT_CCE\")=" << XHOST.vflag_control.flag("NEGLECT_CCE") << endl;  //CO20210115
    XHOST.vflag_control.flag("FORCE_POCC",aurostd::args2attachedflag(argv,"--force_pocc")); //CO20210115
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"FORCE_POCC\")=" << XHOST.vflag_control.flag("FORCE_POCC") << endl;  //CO20210115

    // [CT20200320] run full AEL post-processing for POCC
    XHOST.vflag_control.flag("AEL_RUN_POSTPROCESSING",aurostd::args2flag(XHOST.argv,cmds,"--ael_run_postprocessing"));  //CT20200320
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AEL_RUN_POSTPROCESSING\")=" << XHOST.vflag_control.flag("AEL_RUN_POSTPROCESSING") << endl;  //CT20200320
    // [CT20200320] write full results for POCC+AEL
    XHOST.vflag_control.flag("AEL_WRITE_FULL_RESULTS",aurostd::args2flag(XHOST.argv,cmds,"--ael_write_full_results"));  //CT20200320
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AEL_WRITE_FULL_RESULTS\")=" << XHOST.vflag_control.flag("AEL_FULL_RESULTS") << endl;  //CT20200320
    // [CT20200323] run full AGL post-processing for POCC
    XHOST.vflag_control.flag("AGL_RUN_POSTPROCESSING",aurostd::args2flag(XHOST.argv,cmds,"--agl_run_postprocessing"));  //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AGL_RUN_POSTPROCESSING\")=" << XHOST.vflag_control.flag("AGL_RUN_POSTPROCESSING") << endl;  //CT20200323
    // [CT20200323] write full results for POCC+AGL
    XHOST.vflag_control.flag("AGL_WRITE_FULL_RESULTS",aurostd::args2flag(XHOST.argv,cmds,"--agl_write_full_results"));  //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AGL_WRITE_FULL_RESULTS\")=" << XHOST.vflag_control.flag("AGL_FULL_RESULTS") << endl;  //CT20200323
    //[CT20200323] set number of AGL temperatures from command line
    XHOST.vflag_control.flag("AGL_NTEMPERATURE",aurostd::args2attachedflag(argv,"--agl_ntemperature="));  //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AGL_NTEMPERATURE\")=" << XHOST.vflag_control.flag("AGL_NTEMPERATURE") << endl;  //CT20200323
    if(XHOST.vflag_control.flag("AGL_NTEMPERATURE")) XHOST.vflag_control.push_attached("AGL_NTEMPERATURE",aurostd::args2attachedstring(argv,"--agl_ntemperature=","201")); //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"AGL_NTEMPERATURE\")=" << XHOST.vflag_control.getattachedscheme("AGL_NTEMPERATURE") << endl;  //CT20200323
    //[CT20200323] set AGL temperature step size from command line
    XHOST.vflag_control.flag("AGL_STEMPERATURE",aurostd::args2attachedflag(argv,"--agl_stemperature="));  //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AGL_STEMPERATURE\")=" << XHOST.vflag_control.flag("AGL_STEMPERATURE") << endl;  //CT20200323
    if(XHOST.vflag_control.flag("AGL_STEMPERATURE")) XHOST.vflag_control.push_attached("AGL_STEMPERATURE",aurostd::args2attachedstring(argv,"--agl_stemperature=","10.0")); //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"AGL_NTEMPERATURE\")=" << XHOST.vflag_control.getattachedscheme("AGL_NTEMPERATURE") << endl;  //CT20200323
    //[CT20200323] set number of AGL pressures from command line
    XHOST.vflag_control.flag("AGL_NPRESSURE",aurostd::args2attachedflag(argv,"--agl_npressure="));  //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AGL_NPRESSURE\")=" << XHOST.vflag_control.flag("AGL_NPRESSURE") << endl;  //CT20200323
    if(XHOST.vflag_control.flag("AGL_NPRESSURE")) XHOST.vflag_control.push_attached("AGL_NPRESSURE",aurostd::args2attachedstring(argv,"--agl_npressure=","101")); //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"AGL_NPRESSURE\")=" << XHOST.vflag_control.getattachedscheme("AGL_NPRESSURE") << endl;  //CT20200323
    //[CT20200323] set AGL pressure step size from command line
    XHOST.vflag_control.flag("AGL_SPRESSURE",aurostd::args2attachedflag(argv,"--agl_spressure="));  //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.flag(\"AGL_SPRESSURE\")=" << XHOST.vflag_control.flag("AGL_SPRESSURE") << endl;  //CT20200323
    if(XHOST.vflag_control.flag("AGL_SPRESSURE")) XHOST.vflag_control.push_attached("AGL_SPRESSURE",aurostd::args2attachedstring(argv,"--agl_spressure=","1.0")); //CT20200323
    if(INIT_VERBOSE) oss << "XHOST.vflag_control.getattachedscheme(\"AGL_SPRESSURE\")=" << XHOST.vflag_control.getattachedscheme("AGL_SPRESSURE") << endl;  //CT20200323

    XHOST.AVOID_RUNNING_VASP=aurostd::args2attachedflag(argv,cmds,"--avoid_running_vasp|--no_vasp|--novasp");  //CO20200624 - VERY important, prevents VASP from running
    if( XHOST.GENERATE_AFLOWIN_ONLY ||
        XHOST.POSTPROCESS || 
        XHOST.ARUN_POSTPROCESS ||
        XHOST.vflag_control.flag("AEL_RUN_POSTPROCESSING") || //CT20200722
        XHOST.vflag_control.flag("AGL_RUN_POSTPROCESSING") || //CT20200722
        FALSE) XHOST.AVOID_RUNNING_VASP=TRUE;  //CO20200624
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " XHOST.AVOID_RUNNING_VASP=" << XHOST.AVOID_RUNNING_VASP << endl;}

    XHOST.vflag_control.flag("XPLUG_DO_CLEAN",aurostd::args2flag(XHOST.argv,cmds,"--doclean"));
    XHOST.vflag_control.flag("XPLUG_DO_ADD",aurostd::args2flag(argv,"--add"));

    XHOST.vflag_control.flag("XPLUG_PREFIX",aurostd::args2attachedflag(argv,"--prefix="));
    if(XHOST.vflag_control.flag("XPLUG_PREFIX")) XHOST.vflag_control.push_attached("XPLUG_PREFIX",aurostd::args2attachedstring(argv,"--prefix=",""));
    XHOST.vflag_control.flag("XPLUG_NUM_ZIP",aurostd::args2attachedflag(argv,"--nzip="));
    if(XHOST.vflag_control.flag("XPLUG_NUM_ZIP")) XHOST.vflag_control.push_attached("XPLUG_NUM_ZIP",aurostd::args2attachedstring(argv,"--nzip=",""));
    if(!XHOST.vflag_control.flag("XPLUG_NUM_ZIP")) XHOST.vflag_control.push_attached("XPLUG_NUM_ZIP",aurostd::utype2string(1));
    XHOST.vflag_control.flag("XPLUG_NUM_SIZE",aurostd::args2attachedflag(argv,"--nsize="));
    if(XHOST.vflag_control.flag("XPLUG_NUM_SIZE")) XHOST.vflag_control.push_attached("XPLUG_NUM_SIZE",aurostd::args2attachedstring(argv,"--nsize=",""));
    if(!XHOST.vflag_control.flag("XPLUG_NUM_SIZE")) XHOST.vflag_control.push_attached("XPLUG_NUM_SIZE",aurostd::utype2string(128));
    XHOST.vflag_control.flag("XPLUG_NUM_THREADS",aurostd::args2attachedflag(argv,"--np="));
    XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX",aurostd::args2attachedflag(argv,"--npmax")); //CO20180124
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")||XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX")){
      XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS",aurostd::args2attachedstring(argv,"--np=","1")); //SC20200319
      if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX")||aurostd::toupper(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"))=="MAX"){ //so both --npmax and --np=max work
        XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS","MAX"); //ME20181113
      }
    }
    //  else {XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS",aurostd::utype2string(XHOST.CPU_Cores/2));  OBSOLETE ME20181113  //[CO20200106 - close bracket for indenting]}
    
    // USEFUL shortcuts //SC20200319
    if(!aurostd::args2attachedflag(argv,"--np=")) {
      deque<string> vshort; 
      for(uint ishort=0;ishort<=128;ishort++)
        vshort.push_back(aurostd::utype2string(ishort));
      for(uint ishort=0;ishort<vshort.size();ishort++) {
        if(aurostd::args2flag(argv,cmds,"--multi="+vshort.at(ishort))) {  //SC20200319
          XHOST.vflag_control.flag("MULTI=SH",TRUE);
          XHOST.vflag_control.flag("XPLUG_NUM_THREADS",TRUE);
          XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS",vshort.at(ishort));
          //	  if(INIT_VERBOSE)
          cerr << XPID << "init::InitMachine: FOUND MULTI=SH with np=" << XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS") << endl;
        }
      }
    }
    // USEFUL shortcuts //SC20200323 
    if(!XHOST.vflag_control.flag("FILE")) { // not specified file
      deque<string> vshort; aurostd::string2tokens("0,1,2,3,4,5,6,7,8,9",vshort,",");
      for(uint ishort=0;ishort<vshort.size();ishort++) {
        if(aurostd::args2flag(argv,cmds,"--multi=x.lib"+vshort.at(ishort)) || aurostd::args2flag(argv,cmds,"--multi=./x.lib"+vshort.at(ishort))) {  //SC20200323
          XHOST.vflag_control.flag("MULTI=SH",TRUE);
          XHOST.vflag_control.flag("XPLUG_NUM_THREADS",TRUE);
          XHOST.vflag_control.push_attached("XPLUG_NUM_THREADS","16");
          XHOST.vflag_control.flag("FILE",TRUE);  // if found
          XHOST.vflag_control.push_attached("FILE","x.lib"+vshort.at(ishort));
          cerr << XPID << "init::InitMachine: FOUND FILE=" << XHOST.vflag_control.getattachedscheme("FILE") << endl;
          cerr << XPID << "init::InitMachine: FOUND MULTI=SH with np=" << XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS") << endl;
        }
      }
    }

    //ME20181103 - set MPI when number of threads is larger than 1
    if (XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX") || 
        (aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS")) > 1)) {
      XHOST.MPI = true;
    }

    // AFLOWLIB_SERVER
    XHOST.vflag_control.flag("AFLOWLIB_SERVER",aurostd::args2attachedflag(argv,"--server="));
    if(XHOST.vflag_control.flag("AFLOWLIB_SERVER")) XHOST.vflag_control.push_attached("AFLOWLIB_SERVER",aurostd::args2attachedstring(argv,"--server=",""));
    if(XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="default" || XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="aflowlib") {
      XHOST.vflag_control.pop_attached("AFLOWLIB_SERVER");
      XHOST.vflag_control.push_attached("AFLOWLIB_SERVER","aflowlib.duke.edu");
    }
    if(XHOST.vflag_control.flag("AFLOWLIB_SERVER") &&
        !(XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="aflowlib.duke.edu" || XHOST.vflag_control.getattachedscheme("AFLOWLIB_SERVER")=="materials.duke.edu")) {
      message = XPID + "\"--server=\" can only be \"aflowlib.duke.edu\" or \"materials.duke.edu\"";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
    // LOAD options
    if(INIT_VERBOSE) oss << "--- LOADING @ options --- " << endl;
    //if(INIT_VERBOSE) XHOST.DEBUG=TRUE;
    if(INIT_VERBOSE) oss << "--- LOADING @ aflow options --- " << endl;
    XHOST.vflag_aflow.flag("LOOP",aurostd::args2flag(XHOST.argv,cmds,"-loop|--loop"));
    XHOST.vflag_aflow.flag("CLEAN",aurostd::args2flag(XHOST.argv,cmds,"-c|-clean|--CLEAN|--clean"));
    XHOST.vflag_aflow.args2addattachedscheme(XHOST.argv,cmds,"XCLEAN","--xclean=","");

    XHOST.AFLOW_RUNDIRflag=aurostd::args2flag(XHOST.argv,"--run|-run");
    XHOST.AFLOW_MULTIflag=aurostd::args2flag(XHOST.argv,"--run=multi|-run=multi|--multi|-multi");	
    XHOST.AFLOW_RUNXflag=!XHOST.AFLOW_MULTIflag && (aurostd::args2attachedflag(XHOST.argv,"--run=") || aurostd::args2attachedflag(XHOST.argv,"-run="));

    XHOST.AFLOW_RUNXnumber=0;
    XHOST.vflag_pflow.clear(); 
    XHOST.vflag_apennsy.clear(); 
    XHOST.vflag_outreach.clear(); 

    // ME20210206 - Load --web_mode before parsing arguments
    // LOADING ANRL WEB
    XHOST.vflag_control.flag("WWW",aurostd::args2flag(argv,cmds,"--www|--web|--web_mode|--php|--html|-www|-web|-web_mode|-php|-html"));  //CO20200404
    if(XHOST.user=="www-data"){XHOST.vflag_control.flag("WWW",true);} //CO20201215

    if(INIT_VERBOSE) oss << "--- LOADING @ aconvasp options --- " << endl;
    PflowARGs(XHOST.argv,cmds,XHOST.vflag_pflow);
    if(INIT_VERBOSE) oss << "--- LOADING @ apennsy options --- " << endl;
    ApennsyARGs(XHOST.argv,cmds,XHOST.vflag_apennsy); 
    // LOADING OUTREACH
    if(INIT_VERBOSE) oss << "--- LOADING @ outreach options --- " << endl;   
    XHOST.vflag_control.flag("CV::PUBS",aurostd::args2flag(argv,cmds,"--cv=pubs"));
    XHOST.vflag_control.flag("CV::ITALKS",aurostd::args2flag(argv,cmds,"--cv=italks|--cv=talks|--presentations"));
    XHOST.vflag_control.flag("CV::ACADEMIC",aurostd::args2flag(argv,cmds,"--cv=academic"));
    XHOST.vflag_control.flag("CV::RESEARCH",aurostd::args2flag(argv,cmds,"--cv=research"));
    XHOST.vflag_control.flag("CV::EDUCATION",aurostd::args2flag(argv,cmds,"--cv=education"));
    XHOST.vflag_control.flag("CV::TEACHING",aurostd::args2flag(argv,cmds,"--cv=teaching"));
    XHOST.vflag_control.flag("CV::ADVISING",aurostd::args2flag(argv,cmds,"--cv=advising"));
    XHOST.vflag_control.flag("CV::AWARDS",aurostd::args2flag(argv,cmds,"--cv=awards"));
    XHOST.vflag_control.flag("CV::PRESS",aurostd::args2flag(argv,cmds,"--cv=press"));
    XHOST.vflag_control.flag("CV::PATENTS",aurostd::args2flag(argv,cmds,"--cv=patents"));
    XHOST.vflag_control.flag("CV::SERVICE_OUTSIDE",aurostd::args2flag(argv,cmds,"--cv=service_outside|--cv=serviceoutside"));
    XHOST.vflag_control.flag("CV::SERVICE_INSIDE",aurostd::args2flag(argv,cmds,"--cv=service_inside|--cv=serviceinside"));
    XHOST.vflag_control.flag("PHP::CENTER_MISSION",aurostd::args2flag(argv,cmds,"--php_center_mission|--center_mission|--center-mission"));
    XHOST.vflag_control.flag("CV::AUTHOR",aurostd::args2attachedflag(argv,"--author="));
    if(XHOST.vflag_control.flag("CV::AUTHOR")) XHOST.vflag_control.push_attached("CV::AUTHOR",aurostd::args2attachedstring(argv,"--author=",""));
    XHOST.vflag_control.flag("PHP::PUBS_ALLOY",aurostd::args2attachedflag(argv,"--php_pubs_alloy="));
    if(XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) XHOST.vflag_control.push_attached("PHP::PUBS_ALLOY",aurostd::args2attachedstring(argv,"--php_pubs_alloy=",""));
    XHOST.vflag_control.flag("PHP::PUBS_KEYWORD",aurostd::args2attachedflag(argv,"--php_pubs_keyword="));
    if(XHOST.vflag_control.flag("PHP::PUBS_KEYWORD")) XHOST.vflag_control.push_attached("PHP::PUBS_KEYWORD",aurostd::args2attachedstring(argv,"--php_pubs_keyword=",""));
    XHOST.vflag_control.flag("GRANTS",aurostd::args2flag(argv,cmds,"--grant|--grants"));
    if(XHOST.vflag_control.flag("GRANTS")) XHOST.vflag_control.push_attached("GRANTS",aurostd::args2attachedstring(argv,"--grant=|--grants=",""));
    if(LDEBUG) cout << "OUTREACH OPTIONS: xscheme=" << XHOST.vflag_control.xscheme << endl;
    //    if(LDEBUG) cout << "OUTREACH OPTIONS: vxscheme.size()=" << XHOST.vflag_control.vxscheme.size() << endl;  OBSOLETE ME20181102
    if(LDEBUG) cout << "OUTREACH OPTIONS: vxsghost.size()=" << XHOST.vflag_control.vxsghost.size() << endl;
    if(LDEBUG) cout << "OUTREACH OPTIONS: argv.size()=" << argv.size() << endl;

    //FANCY_PRINT
    XHOST.vflag_control.flag("NO_FANCY_PRINT",aurostd::args2flag(argv,cmds,"--no_fancy_print|--nofancyprint"));  //CO20200404

    // DEFAULT options
    if(INIT_VERBOSE) oss << "--- DEFAULTSs --- " << endl;
    if(INIT_VERBOSE) aflowrc::print_aflowrc(oss,INIT_VERBOSE || XHOST.DEBUG);
    //if(INIT_VERBOSE) XHOST.DEBUG=TRUE;

    // FINISHED
    if(INIT_VERBOSE) oss << endl;
    if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
    if(INIT_VERBOSE) oss << "* AFLOW V=" << string(AFLOW_VERSION) << " - machine information " << endl;
    if(INIT_VERBOSE) oss << "*********************************************************************************" << endl;
    if(INIT_VERBOSE) return 0;
    // CHECK CRC
    // aurostd::crc64_main();

    // NOW LOAD schema
    if (init::InitSchema(INIT_VERBOSE) == 0) return 0;
    init::InitSchemaInternal(INIT_VERBOSE);

    // DONE
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitMachine: [END]" << endl;

    return -1;
  }
} // namespace init

// ***************************************************************************
// long init::GetRAM(void)
// ***************************************************************************
namespace init {
#ifndef _MACOSX_
#include <sys/sysinfo.h>
  long GetRAM(void) {
    long pages=sysconf(_SC_PHYS_PAGES);
    long page_size=sysconf(_SC_PAGE_SIZE);
    return pages*page_size;
  }
  long _GetRAM(void) {
    struct sysinfo s;
    if(sysinfo(&s)!=0) {
      string message = "sysinfo error";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return s.totalram;
  }
#endif
#ifdef _MACOSX_
#include <sys/sysctl.h>
  long GetRAM(void) {
    int mib[2]={CTL_HW,HW_MEMSIZE};
    u_int namelen=sizeof(mib)/sizeof(mib[0]);
    uint64_t size;
    size_t len=sizeof(size);
    if(sysctl(mib,namelen,&size,&len,NULL,0)<0) {
      string message = "sysctl returned an error";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return (long) size;
  }
#endif
} // namespace init

// ***************************************************************************
// init::InitLoadString
// ***************************************************************************
namespace init {
  string InitLoadString(string str2load,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || LVERBOSE);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " str2load=" << str2load << endl; 

    if((str2load=="vLIBS" || str2load=="XHOST_vLIBS") && XHOST_vLIBS.size()==3) return ""; // intercept before it reloads it again

    if(!XHOST.is_command("aflow_data")) {
      string message = "aflow_data is not in the path.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } 
    if(LDEBUG) cerr << "00000  MESSAGE AFLOW INIT Loading data = [" << str2load << "]";
    if(LDEBUG) cerr.flush();

    string out;
    string aflow_data_path=aurostd::args2attachedstring(XHOST.argv,"--aflow_data_path=",(string) "");
    if(aflow_data_path=="") {
      if(LDEBUG) cerr << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST.hostname=" << XHOST.hostname << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST.user=" << XHOST.user << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;

      if(XHOST.hostname=="nietzsche.mems.duke.edu"&&XHOST.user=="auro"&&aurostd::FileExist(XHOST.home+"/work/AFLOW3/aflow_data")) {  //CO, special SC
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " FOUND " << XHOST.home << "/work/AFLOW3/aflow_data" << endl;
        //	if(LDEBUG) cerr << __AFLOW_FUNC__ << " out.length()=" << out.length() << endl;
        out=aurostd::execute2string(XHOST.home+"/work/AFLOW3/aflow_data"+" "+str2load);
        if(LDEBUG) cerr << __AFLOW_FUNC__ << " out.length()=" << out.length() << endl;
        // 	vector<string> vout;aurostd::string2vectorstring(out,vout);
        // 	if(LDEBUG) cerr << __AFLOW_FUNC__ << " vout.size()=" << vout.size() << endl;
        // 	for(uint i=0;i<vout.size();i+=3) {
        // 	  if(aurostd::substring2bool(vout.at(i),"LIB6") || aurostd::substring2bool(vout.at(i+1),"LIB6") || aurostd::substring2bool(vout.at(i+2),"LIB6")) 
        // 	    cout << vout.at(i) << "    " << vout.at(i+1) << "    " << vout.at(i+2) << endl;
        // 	}
        // 	//       if(LDEBUG) cerr << __AFLOW_FUNC__ << " out=" << out << endl; 
        // if(LDEBUG) cerr << __AFLOW_FUNC__ << " str2load=" << str2load << endl;
      } else {
        //	cerr <<  __AFLOW_FUNC__ << " [2] " << endl;
        if(LDEBUG) {cerr << __AFLOW_FUNC__ << " issuing command: " << XHOST.command("aflow_data") << " " << str2load << endl;}
        out=aurostd::execute2string(XHOST.command("aflow_data")+" "+str2load);
      }
    } else { // cerr << string(aflow_data_path+"/"+XHOST.command("aflow_data")) << endl;
      //      cerr <<  __AFLOW_FUNC__ << " [3] " << endl;
      out=aurostd::execute2string(aflow_data_path+"/"+XHOST.command("aflow_data")+" "+str2load);
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " out.length()=" << out.length() << endl;
    if(LDEBUG) cerr.flush();
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_vLIBS.size()=" << XHOST_vLIBS.size() << " [1]" << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " str2load=" << str2load << " [1]" << endl;

    if((str2load=="vLIBS" || str2load=="XHOST_vLIBS")) { // && XHOST_vLIBS.size()!=3)
      if(XHOST_vLIBS.size()) for(uint i=0;i<XHOST_vLIBS.size();i++) XHOST_vLIBS.at(i).clear();
      XHOST_vLIBS.clear();
      XHOST_vLIBS.push_back(vector<string>()); // AURL
      XHOST_vLIBS.push_back(vector<string>()); // AUID
      XHOST_vLIBS.push_back(vector<string>()); // LOOP
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_vLIBS.size()=" << XHOST_vLIBS.size() << " [2]" << endl;
      vector<string> vout;
      string aurl,auid,loop;
      bool found=FALSE;
      aurostd::string2vectorstring(out,vout);
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " vout.size()=" << vout.size() << endl;
      // make some checks if it  is divisible by three

      for(uint i=0;i<vout.size();) {
        aurl=vout.at(i++);XHOST_vLIBS.at(0).push_back(aurl); // AURL
        auid=vout.at(i++);XHOST_vLIBS.at(1).push_back(auid); // AUID
        loop=vout.at(i++);XHOST_vLIBS.at(2).push_back(loop); // LOOP
        aurostd::StringSubst(aurl,"aflowlib.duke.edu:","");
        aurostd::StringSubst(aurl,"materials.duke.edu:","");
        found=FALSE;
        // do LIB3 first to accelerate
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB3")) { // XHOST_Library_CALCULATED_LIB3
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB3_RAW/","");
          XHOST_Library_CALCULATED_LIB3_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB3_RAW.push_back(aurl);
        }
        // do LIB4 second to accelerate
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB4")) { // XHOST_Library_CALCULATED_LIB4
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB4_RAW/","");
          XHOST_Library_CALCULATED_LIB4_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB4_RAW.push_back(aurl);
        }
        // do LIB2 third to accelerate
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB2")) { // XHOST_Library_CALCULATED_LIB2
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB2_RAW/","");
          XHOST_Library_CALCULATED_LIB2_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB2_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/ICSD")) { // XHOST_Library_CALCULATED_ICSD
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/ICSD_RAW/","");
          aurostd::StringSubst(aurl,"AFLOWDATA/ICSD_WEB/","");
          XHOST_Library_CALCULATED_ICSD_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_ICSD_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB0")) { // XHOST_Library_CALCULATED_LIB0
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB0_RAW/","");
          XHOST_Library_CALCULATED_LIB0_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB0_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB1")) { // XHOST_Library_CALCULATED_LIB1
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB1_RAW/","");
          XHOST_Library_CALCULATED_LIB1_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB1_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB5")) { // XHOST_Library_CALCULATED_LIB5
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB5_RAW/","");
          XHOST_Library_CALCULATED_LIB5_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB5_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB6")) { // XHOST_Library_CALCULATED_LIB6
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB6_RAW/","");
          XHOST_Library_CALCULATED_LIB6_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB6_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB7")) { // XHOST_Library_CALCULATED_LIB7
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB7_RAW/","");
          XHOST_Library_CALCULATED_LIB7_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB7_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB8")) { // XHOST_Library_CALCULATED_LIB8
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB8_RAW/","");
          XHOST_Library_CALCULATED_LIB8_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB8_RAW.push_back(aurl);
        }
        if(!found) if(aurostd::substring2bool(aurl,"AFLOWDATA/LIB9")) { // XHOST_Library_CALCULATED_LIB9
          found=TRUE;
          aurostd::StringSubst(aurl,"AFLOWDATA/LIB9_RAW/","");
          XHOST_Library_CALCULATED_LIB9_LIB.push_back(aurl);
          XHOST_Library_CALCULATED_LIB9_RAW.push_back(aurl);
        }
      }
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_vLIBS.at(0).size()=" << XHOST_vLIBS.at(0).size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_vLIBS.at(1).size()=" << XHOST_vLIBS.at(1).size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_vLIBS.at(2).size()=" << XHOST_vLIBS.at(2).size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_ICSD_LIB.size()=" << XHOST_Library_CALCULATED_ICSD_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB0_LIB.size()=" << XHOST_Library_CALCULATED_LIB0_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB1_LIB.size()=" << XHOST_Library_CALCULATED_LIB1_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB2_LIB.size()=" << XHOST_Library_CALCULATED_LIB2_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB3_LIB.size()=" << XHOST_Library_CALCULATED_LIB3_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB4_LIB.size()=" << XHOST_Library_CALCULATED_LIB4_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB5_LIB.size()=" << XHOST_Library_CALCULATED_LIB5_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB6_LIB.size()=" << XHOST_Library_CALCULATED_LIB6_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB7_LIB.size()=" << XHOST_Library_CALCULATED_LIB7_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB8_LIB.size()=" << XHOST_Library_CALCULATED_LIB8_LIB.size() << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " XHOST_Library_CALCULATED_LIB9_LIB.size()=" << XHOST_Library_CALCULATED_LIB9_LIB.size() << endl;
    }
    return out; 
  }
} // namespace init

// [OBSOLETE] bool vAURL_cutout(vector<string>& vvAURL,string cutout,bool LVERBOSE) {
// [OBSOLETE]   bool LDEBUG=(FALSE || XHOST.DEBUG || LVERBOSE);
// [OBSOLETE]   stringstream sss;
// [OBSOLETE]   if(LDEBUG) cerr << "vAUID_cutout: XHOST_vAUID.size()=" << XHOST_vAUID.size() << endl;
// [OBSOLETE]   if(LDEBUG) cerr << "vAURL_cutout: XHOST_vAURL.size()=" << XHOST_vAURL.size() << endl;
// [OBSOLETE]   if(LDEBUG) cerr << "vLOOP_cutout: XHOST_vLOOP.size()=" << XHOST_vLOOP.size() << endl;
// [OBSOLETE]   if(LDEBUG) cerr << "vLIBS_cutout: XHOST_vLIBS.size()=" << XHOST_vLIBS.size() << endl;
// [OBSOLETE]   vector<string> tokens;
// [OBSOLETE]   aurostd::string2tokens(cutout,tokens,"|");
// [OBSOLETE]   string aus;
// [OBSOLETE]   if(!XHOST_vAURL.size()) init::InitLoadString("vLIBS",LVERBOSE);
// [OBSOLETE]   for(uint i=0;i<tokens.size();i++) {
// [OBSOLETE]     for(uint j=0;j<XHOST_vAURL.size();j++) {
// [OBSOLETE]       aus=XHOST_vAURL.at(j);
// [OBSOLETE]       if(aurostd::substring2bool(aus,tokens.at(i))) {
// [OBSOLETE]         aurostd::StringSubst(aus,tokens.at(i),"");
// [OBSOLETE]         aurostd::StringSubst(aus,"aflowlib.duke.edu:","");
// [OBSOLETE]         aurostd::StringSubst(aus,"materials.duke.edu:","");
// [OBSOLETE]         vvAURL.push_back(aus);
// [OBSOLETE]       }
// [OBSOLETE]     }
// [OBSOLETE]   }
// [OBSOLETE]   return true;
// [OBSOLETE] }
// [OBSOLETE] 
// [OBSOLETE] string vAURL_cutout(string cutout,bool LVERBOSE) {
// [OBSOLETE]   bool LDEBUG=(FALSE || XHOST.DEBUG || LVERBOSE);
// [OBSOLETE]   stringstream sss;
// [OBSOLETE]   if(LDEBUG) cerr << "vAUID_cutout: XHOST_vAUID.size()=" << XHOST_vAUID.size() << endl;
// [OBSOLETE]   if(LDEBUG) cerr << "vAURL_cutout: XHOST_vAURL.size()=" << XHOST_vAURL.size() << endl;
// [OBSOLETE]   if(LDEBUG) cerr << "vLOOP_cutout: XHOST_vLOOP.size()=" << XHOST_vLOOP.size() << endl;
// [OBSOLETE]   if(LDEBUG) cerr << "vLIBS_cutout: XHOST_vLIBS.size()=" << XHOST_vLIBS.size() << endl;
// [OBSOLETE]   vector<string> tokens;
// [OBSOLETE]   aurostd::string2tokens(cutout,tokens,"|");
// [OBSOLETE]   string aus;
// [OBSOLETE]   if(!XHOST_vAURL.size()) init::InitLoadString("vLIBS",LVERBOSE);
// [OBSOLETE]   for(uint i=0;i<tokens.size();i++) {
// [OBSOLETE]     for(uint j=0;j<XHOST_vAURL.size();j++) {
// [OBSOLETE]       aus=XHOST_vAURL.at(j);
// [OBSOLETE]       if(aurostd::substring2bool(aus,tokens.at(i))) {
// [OBSOLETE]         aurostd::StringSubst(aus,tokens.at(i),"");
// [OBSOLETE]         aurostd::StringSubst(aus,"aflowlib.duke.edu:","");
// [OBSOLETE]         aurostd::StringSubst(aus,"materials.duke.edu:","");
// [OBSOLETE] 	sss << aus << endl;
// [OBSOLETE]       }
// [OBSOLETE]     }
// [OBSOLETE]   }
// [OBSOLETE]   return sss.str();
// [OBSOLETE] }

// ***************************************************************************
// init::InitGlobalObject
// ***************************************************************************
namespace init {
  string InitGlobalObject(string str,string grep,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // || LVERBOSE;
    string out="";
    string FileLibrary="";

    if(str=="Library_HTQC") { 
      if(XHOST_Library_HTQC.empty()) {
        return XHOST_Library_HTQC=init::InitLoadString(str,LVERBOSE);
      } else { 
        return XHOST_Library_HTQC;
      }
    } // FIX
    // FILES CALCULATED
    if(str=="Library_CALCULATED_ICSD_LIB" || str=="Library_CALCULATED_ICSD_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB0_LIB" || str=="Library_CALCULATED_LIB0_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB1_LIB" || str=="Library_CALCULATED_LIB1_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB2_LIB" || str=="Library_CALCULATED_LIB2_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB3_LIB" || str=="Library_CALCULATED_LIB3_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB4_LIB" || str=="Library_CALCULATED_LIB4_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB5_LIB" || str=="Library_CALCULATED_LIB5_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB6_LIB" || str=="Library_CALCULATED_LIB6_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB7_LIB" || str=="Library_CALCULATED_LIB7_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB8_LIB" || str=="Library_CALCULATED_LIB8_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    if(str=="Library_CALCULATED_LIB9_LIB" || str=="Library_CALCULATED_LIB9_RAW") { init::InitLoadString("vLIBS",LVERBOSE); }
    // AUID AURL LOOP LIBS
    if(str=="vLIBS" || str=="XHOST_vLIBS") { init::InitLoadString("vLIBS",LVERBOSE);} // just make them all
    // AFLOWLIB THINGS
    // LOAD
    // if(str=="aflowlib_lib0") { if(XHOST_aflowlib_lib0.empty()) { return XHOST_aflowlib_lib0=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib0;}} // 
    // if(str=="aflowlib_lib1") { if(XHOST_aflowlib_lib1.empty()) { return XHOST_aflowlib_lib1=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib1;}} //
    // if(str=="aflowlib_lib2") { if(XHOST_aflowlib_lib2.empty()) { return XHOST_aflowlib_lib2=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib2;}} //
    // if(str=="aflowlib_lib3") { if(XHOST_aflowlib_lib3.empty()) { return XHOST_aflowlib_lib3=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib3;}} // 
    // if(str=="aflowlib_lib4") { if(XHOST_aflowlib_lib4.empty()) { return XHOST_aflowlib_lib4=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib4;}} // 
    // if(str=="aflowlib_lib5") { if(XHOST_aflowlib_lib5.empty()) { return XHOST_aflowlib_lib5=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib5;}} // 
    // if(str=="aflowlib_lib6") { if(XHOST_aflowlib_lib6.empty()) { return XHOST_aflowlib_lib6=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib6;}} //  
    // if(str=="aflowlib_lib7") { if(XHOST_aflowlib_lib7.empty()) { return XHOST_aflowlib_lib7=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib7;}} //  
    // if(str=="aflowlib_lib8") { if(XHOST_aflowlib_lib8.empty()) { return XHOST_aflowlib_lib8=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib8;}} //  
    // if(str=="aflowlib_lib9") { if(XHOST_aflowlib_lib9.empty()) { return XHOST_aflowlib_lib9=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_lib9;}} //  
    // if(str=="aflowlib_icsd") { if(XHOST_aflowlib_icsd.empty()) { return XHOST_aflowlib_icsd=init::InitLoadString(str,LVERBOSE);} else { return XHOST_aflowlib_icsd;}} // 
    // LOAD
    // STOKES THINGS
    if(str=="FINDSYM_data_space_txt") { if(XHOST_FINDSYM_data_space_txt.empty()) { return XHOST_FINDSYM_data_space_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FINDSYM_data_space_txt;}} // LOADED TXTS
    if(str=="FINDSYM_data_wyckoff_txt") { if(XHOST_FINDSYM_data_wyckoff_txt.empty()) { return XHOST_FINDSYM_data_wyckoff_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FINDSYM_data_wyckoff_txt;}} // LOADED TXTS
    if(str=="FROZSL_data_space_txt") { if(XHOST_FROZSL_data_space_txt.empty()) { return XHOST_FROZSL_data_space_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_space_txt;}} // LOADED TXTS
    if(str=="FROZSL_data_wyckoff_txt") { if(XHOST_FROZSL_data_wyckoff_txt.empty()) { return XHOST_FROZSL_data_wyckoff_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_wyckoff_txt;}} // LOADED TXTS
    if(str=="FROZSL_data_images_txt") { if(XHOST_FROZSL_data_images_txt.empty()) { return XHOST_FROZSL_data_images_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_images_txt;}} // LOADED TXTS
    if(str=="FROZSL_data_irreps_txt") { if(XHOST_FROZSL_data_irreps_txt.empty()) { return XHOST_FROZSL_data_irreps_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_irreps_txt;}} // LOADED TXTS
    if(str=="FROZSL_data_isotropy_txt") { if(XHOST_FROZSL_data_isotropy_txt.empty()) { return XHOST_FROZSL_data_isotropy_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_isotropy_txt;}} // LOADED TXTS
    if(str=="FROZSL_data_little_txt") { if(XHOST_FROZSL_data_little_txt.empty()) { return XHOST_FROZSL_data_little_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_data_little_txt;}} // LOADED TXTS
    if(str=="FROZSL_symmetry2_dat") { if(XHOST_FROZSL_symmetry2_dat.empty()) { return XHOST_FROZSL_symmetry2_dat=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_symmetry2_dat;}} // LOADED TXTS
    if(str=="FROZSL_const_dat") { if(XHOST_FROZSL_const_dat.empty()) { return XHOST_FROZSL_const_dat=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_const_dat;}} // LOADED TXTS
    if(str=="FROZSL_phvaspsetup_AFLOW") { if(XHOST_FROZSL_phvaspsetup_AFLOW.empty()) { return XHOST_FROZSL_phvaspsetup_AFLOW=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_phvaspsetup_AFLOW;}} // LOADED TXTS
    if(str=="FROZSL_phvaspsetup_POSCAR") { if(XHOST_FROZSL_phvaspsetup_POSCAR.empty()) { return XHOST_FROZSL_phvaspsetup_POSCAR=init::InitLoadString(str,LVERBOSE);} else { return XHOST_FROZSL_phvaspsetup_POSCAR;}} // LOADED TXTS
    // README THINGS
    if(str=="README_AFLOW_LICENSE_GPL3_TXT") { if(XHOST_README_AFLOW_LICENSE_GPL3_TXT.empty()) { return XHOST_README_AFLOW_LICENSE_GPL3_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_LICENSE_GPL3_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_TXT") { if(XHOST_README_AFLOW_TXT.empty()) { return XHOST_README_AFLOW_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_PFLOW_TXT") { if(XHOST_README_AFLOW_PFLOW_TXT.empty()) { return XHOST_README_AFLOW_PFLOW_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_PFLOW_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_APENNSY_TXT") { if(XHOST_README_AFLOW_APENNSY_TXT.empty()) { return XHOST_README_AFLOW_APENNSY_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_APENNSY_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_SCRIPTING_TXT") { if(XHOST_README_AFLOW_SCRIPTING_TXT.empty()) { return XHOST_README_AFLOW_SCRIPTING_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_SCRIPTING_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_FROZSL_TXT") { if(XHOST_README_AFLOW_FROZSL_TXT.empty()) { return XHOST_README_AFLOW_FROZSL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_FROZSL_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_POCC_TXT") { if(XHOST_README_AFLOW_POCC_TXT.empty()) { return XHOST_README_AFLOW_POCC_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_POCC_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_APL_TXT") { if(XHOST_README_AFLOW_APL_TXT.empty()) { return XHOST_README_AFLOW_APL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_APL_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_AGL_TXT") { if(XHOST_README_AFLOW_AGL_TXT.empty()) { return XHOST_README_AFLOW_AGL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_AGL_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_AEL_TXT") { if(XHOST_README_AFLOW_AEL_TXT.empty()) { return XHOST_README_AFLOW_AEL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_AEL_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_ANRL_TXT") { if(XHOST_README_AFLOW_ANRL_TXT.empty()) { return XHOST_README_AFLOW_ANRL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_ANRL_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_COMPARE_TXT") { if(XHOST_README_AFLOW_COMPARE_TXT.empty()) { return XHOST_README_AFLOW_COMPARE_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_COMPARE_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_GFA_TXT") { if(XHOST_README_AFLOW_GFA_TXT.empty()) { return XHOST_README_AFLOW_GFA_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_GFA_TXT;}} // LOADED TXTS  //CO20190401
    if(str=="README_AFLOW_SYM_TXT") { if(XHOST_README_AFLOW_SYM_TXT.empty()) { return XHOST_README_AFLOW_SYM_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_SYM_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_CCE_TXT") { if(XHOST_README_AFLOW_CCE_TXT.empty()) { return XHOST_README_AFLOW_CCE_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_CCE_TXT;}} // LOADED TXTS  //CO20190620
    if(str=="README_AFLOW_CHULL_TXT") { if(XHOST_README_AFLOW_CHULL_TXT.empty()) { return XHOST_README_AFLOW_CHULL_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_CHULL_TXT;}} // LOADED TXTS  //CO20190620
    if(str=="README_AFLOW_EXCEPTIONS_TXT") {if(XHOST_README_AFLOW_EXCEPTIONS_TXT.empty()){ return XHOST_README_AFLOW_EXCEPTIONS_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_EXCEPTIONS_TXT;}}  //ME20180531
    if(str=="README_PROTO_TXT") { if(XHOST_README_PROTO_TXT.empty()) { return XHOST_README_PROTO_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_PROTO_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_XAFLOW_TXT") { if(XHOST_README_AFLOW_XAFLOW_TXT.empty()) { return XHOST_README_AFLOW_XAFLOW_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_XAFLOW_TXT;}} // LOADED TXTS
    if(str=="README_AFLOW_AFLOWRC_TXT") { if(XHOST_README_AFLOW_AFLOWRC_TXT.empty()) { return XHOST_README_AFLOW_AFLOWRC_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_README_AFLOW_AFLOWRC_TXT;}} // LOADED TXTS
    // SCINTILLATION THINGS
    if(str=="ElectronStoppingPower_txt") { if(XHOST_ElectronStoppingPower_txt.empty()) { return XHOST_ElectronStoppingPower_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_ElectronStoppingPower_txt;}} // LOADED TXTS
    if(str=="PhotonCrossSection_txt") { if(XHOST_PhotonCrossSection_txt.empty()) { return XHOST_PhotonCrossSection_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_PhotonCrossSection_txt;}} // LOADED TXTS
    if(str=="PhotonStoppingPower_txt") { if(XHOST_PhotonStoppingPower_txt.empty()) { return XHOST_PhotonStoppingPower_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_PhotonStoppingPower_txt;}} // LOADED TXTS
    if(str=="ICSD_List_txt") { if(XHOST_ICSD_List_txt.empty()) { return XHOST_ICSD_List_txt=init::InitLoadString(str,LVERBOSE);} else { return XHOST_ICSD_List_txt;}} // LOADED TXTS
    if(str=="AFLOW_PSEUDOPOTENTIALS") { if(XHOST_AFLOW_PSEUDOPOTENTIALS.empty()) { return XHOST_AFLOW_PSEUDOPOTENTIALS=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_PSEUDOPOTENTIALS;}} // LOADED TXTS
    if(str=="AFLOW_PSEUDOPOTENTIALS_TXT") { if(XHOST_AFLOW_PSEUDOPOTENTIALS_TXT.empty()) { return XHOST_AFLOW_PSEUDOPOTENTIALS_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_PSEUDOPOTENTIALS_TXT;}} // LOADED TXTS
    if(str=="AFLOW_PSEUDOPOTENTIALS_LIST_TXT") { if(XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT.empty()) { return XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT=init::InitLoadString(str,LVERBOSE);} else { return XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT;}} // LOADED TXTS
    if(str=="f144468a7ccc2d3a72ba44000715efdb") {
      if(XHOST_f144468a7ccc2d3a72ba44000715efdb.empty()) {
        //	cerr << "init::InitGlobalObject" << " [1a]" << endl;
        XHOST_f144468a7ccc2d3a72ba44000715efdb=init::InitLoadString(str,LVERBOSE);
        //	cerr << "init::InitGlobalObject" << " [1b]" << endl;
        return XHOST_f144468a7ccc2d3a72ba44000715efdb;
      } else {
        //	cerr << "init::InitGlobalObject" << " [2]" << endl;
        return XHOST_f144468a7ccc2d3a72ba44000715efdb;
      }
    } // LOADED TXTS
    // [OBSOLETE] if(str=="d0f1b0e47f178ae627a388d3bf65d2d2") { if(XHOST_d0f1b0e47f178ae627a388d3bf65d2d2.empty()) { return XHOST_d0f1b0e47f178ae627a388d3bf65d2d2=init::InitLoadString(str,LVERBOSE);} else { return XHOST_d0f1b0e47f178ae627a388d3bf65d2d2;}} // LOADED TXTS
    // [OBSOLETE] if(str=="decf00ca3ad2fe494eea8e543e929068") { if(XHOST_decf00ca3ad2fe494eea8e543e929068.empty()) { return XHOST_decf00ca3ad2fe494eea8e543e929068=init::InitLoadString(str,LVERBOSE);} else { return XHOST_decf00ca3ad2fe494eea8e543e929068;}} // LOADED TXTS

    // SEARCH IN AFLOW_DATA AND IF NOT FROM LIBRARIES
    // pure and auro are inside aflow_data
    if(str=="Library_ICSD" || str=="aflowlib_lib0" || str=="aflowlib_lib1" || str=="aflowlib_lib2" || str=="aflowlib_lib3" || str=="aflowlib_lib4" || str=="aflowlib_lib5" || str=="aflowlib_lib6" || str=="aflowlib_icsd") {
      //  cerr << "(*vLibrary).length()=" << (*vLibrary).length() << endl;
      // if(LDEBUG)

      string *vLibrary = NULL;
      if(str=="Library_ICSD")  vLibrary=&XHOST_Library_ICSD_ALL;
      if(str=="aflowlib_icsd") vLibrary=&XHOST_aflowlib_icsd;
      if(str=="aflowlib_lib0") vLibrary=&XHOST_aflowlib_lib1;
      if(str=="aflowlib_lib1") vLibrary=&XHOST_aflowlib_lib1;
      if(str=="aflowlib_lib2") vLibrary=&XHOST_aflowlib_lib2;
      if(str=="aflowlib_lib3") vLibrary=&XHOST_aflowlib_lib3;
      if(str=="aflowlib_lib4") vLibrary=&XHOST_aflowlib_lib4;
      if(str=="aflowlib_lib5") vLibrary=&XHOST_aflowlib_lib5;
      if(str=="aflowlib_lib6") vLibrary=&XHOST_aflowlib_lib6;
      if(str=="aflowlib_lib7") vLibrary=&XHOST_aflowlib_lib7;
      if(str=="aflowlib_lib8") vLibrary=&XHOST_aflowlib_lib8;
      if(str=="aflowlib_lib9") vLibrary=&XHOST_aflowlib_lib9;
      // check SHORTCUTS
      cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitGlobalObject: Loading library \"" << str << "\" from \"ram\"" << endl;
      if(!(*vLibrary).empty()) {return (*vLibrary);}
      // THEN try aflow_data
      cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitGlobalObject: Loading library \"" << str << "\" from \"aflow_data\"" << endl;
      (*vLibrary)=init::InitLoadString(str,LVERBOSE);
      if(!(*vLibrary).empty()) {return (*vLibrary);}
      // THEN try aflow libraries
      cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitGlobalObject: Loading library \"" << str << "\" from \"aflow libs\"" << endl;

      // check if available
      if((*vLibrary).empty()) {   // find and LOAD
        string str2search=str;
        aurostd::StringSubst(str2search,"Library_ICSD","aflow_library_icsd");
        (*vLibrary)="";
        for(uint j=0;j<vAFLOW_LIBRARY_DIRECTORIES.size() && (*vLibrary).empty();j++) {   // cycle through possible directories
          FileLibrary=aurostd::CleanFileName(vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/"+str2search+".dat");
          if(aurostd::FileExist(FileLibrary) && aurostd::FileNotEmpty(FileLibrary)) {
            if(LDEBUG || LVERBOSE) cerr << "00000  AFLOW LIBRARY  (" << j << ")  found=" <<  FileLibrary << endl;
            if(LDEBUG || LVERBOSE) cerr << "loading... ";
            if(LDEBUG || LVERBOSE) cerr.flush();
            if(grep=="") {
              aurostd::file2string(FileLibrary,(*vLibrary));
            } else {
              (*vLibrary)=aurostd::execute2string("cat "+FileLibrary+" | grep -E '"+grep+"'");
            }
            if(LDEBUG || LVERBOSE) cerr << "length=" << (*vLibrary).size();// << " " << endl;
            if(LDEBUG || LVERBOSE) cerr.flush();
          }
        } // cycle through possible directories
        if((*vLibrary).empty()) {
          cerr << "WARNING - init::InitGlobalObject: " << str << " not found! " << endl;
          return "";
        }
        out=(*vLibrary);
      }
    } 

    return out;
  }
} // namespace init

// ***************************************************************************
// init::InitLibraryObject
// ***************************************************************************
namespace init {
  string InitLibraryObject(string str,bool LVERBOSE) {
    bool LDEBUG=FALSE;
    // Search LIBRARY
    string str2search=str;
    aurostd::StringSubst(str2search,"Library_ICSD","aflow_library_icsd.dat");
    if(str=="Library_ICSD") { 
      if(XHOST_Library_ICSD_ALL.empty()) { 
        return XHOST_Library_ICSD_ALL=init::InitLoadString(str,LVERBOSE);
      } else { 
        return XHOST_Library_ICSD_ALL;
      }
    } // FIX

    string FileLibrary,out="";
    for(uint j=0;j<vAFLOW_LIBRARY_DIRECTORIES.size() && out.empty();j++) {   // cycle through possible directories
      FileLibrary=vAFLOW_LIBRARY_DIRECTORIES.at(j)+"/"+str2search;
      if(LDEBUG || LVERBOSE) cerr << "DDDDD  InitLibraryObject: (" << j << ")  " <<  FileLibrary << endl;
      if(aurostd::FileExist(FileLibrary) && aurostd::FileNotEmpty(FileLibrary))
        out=aurostd::file2string(FileLibrary,out);
    } // cycle through possible directories
    if(!out.empty()) {
      if(LDEBUG || LVERBOSE) cerr << "00000  MESSAGE InitLibraryObject: AFLOW LIBRARY  Found library file = [" << FileLibrary << "]" << endl;
    } else {
      cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") initLibraryObject: AFLOW_LIBRARY not found! " << endl;
    }

    return out;
  }
} // namespace init



// ***************************************************************************
// init::AFLOW_Projects_Directories
// ***************************************************************************
namespace init {
  string AFLOW_Projects_Directories(string lib) {
    bool LDEBUG=FALSE;
    if(LDEBUG) {;} //CO20190906 - keep LDEBUG busy
    string out="";
    //ME20200707 - The LIBRARY_NOTHING check is important or this function
    //breaks when the LIB directory does not exist
    if((XHOST_LIBRARY_AUID != LIBRARY_NOTHING) && (aurostd::toupper(lib)=="AUID")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_AUID);
    if((XHOST_LIBRARY_ICSD != LIBRARY_NOTHING) && (aurostd::toupper(lib)=="ICSD")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD);
    if((XHOST_LIBRARY_LIB0 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB0") ||  lib=="0")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB0);
    if((XHOST_LIBRARY_LIB1 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB1") ||  lib=="1")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1);
    if((XHOST_LIBRARY_LIB2 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB2") ||  lib=="2")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2);
    if((XHOST_LIBRARY_LIB3 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB3") ||  lib=="3")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3);
    if((XHOST_LIBRARY_LIB4 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB4") ||  lib=="4")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4);
    if((XHOST_LIBRARY_LIB5 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB5") ||  lib=="5")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5);
    if((XHOST_LIBRARY_LIB6 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB6") ||  lib=="6")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6);
    if((XHOST_LIBRARY_LIB7 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB7") ||  lib=="7")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7);
    if((XHOST_LIBRARY_LIB8 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB8") ||  lib=="8")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8);
    if((XHOST_LIBRARY_LIB9 != LIBRARY_NOTHING) && ((aurostd::toupper(lib)=="LIB9") ||  lib=="9")) out=vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9);
    return out;

    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_AUID)" "init::AFLOW_Projects_Directories(\"AUID\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_ICSD)" "init::AFLOW_Projects_Directories(\"ICSD\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB0)" "init::AFLOW_Projects_Directories(\"LIB0\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB1)" "init::AFLOW_Projects_Directories(\"LIB1\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB2)" "init::AFLOW_Projects_Directories(\"LIB2\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB3)" "init::AFLOW_Projects_Directories(\"LIB3\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB4)" "init::AFLOW_Projects_Directories(\"LIB4\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB5)" "init::AFLOW_Projects_Directories(\"LIB5\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB6)" "init::AFLOW_Projects_Directories(\"LIB6\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB7)" "init::AFLOW_Projects_Directories(\"LIB7\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB8)" "init::AFLOW_Projects_Directories(\"LIB8\")" *cpp
    //subst "vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_LIB9)" "init::AFLOW_Projects_Directories(\"LIB9\")" *cpp

  }
} // namespace init

// ***************************************************************************
// uint init::GetTEMP(void) // need sensors package
// ***************************************************************************
namespace init {
  pthread_mutex_t mutex_INIT_GetTEMP=PTHREAD_MUTEX_INITIALIZER;
  uint GetTEMP(void) {
    pthread_mutex_lock(&mutex_INIT_GetTEMP);
    // pthread_mutex_unlock(&mutex_INIT_GetTEMP);

    // if(aurostd::execute2string("ps aux | grep sensors | grep -v sensorsd | grep -v grep")!="") LOCAL_is_sensor=FALSE; // must postpone
    XHOST.vTemperatureCore.clear();
    if(XHOST.sensors_allowed) {
      bool LOCAL_is_sensor=XHOST.is_command("sensors");
      // test sensors
      if(LOCAL_is_sensor)
        if(aurostd::execute2utype<int>("bash -c \"sensors 2>&1 2> /dev/null\" | grep -c temp")==0)
          LOCAL_is_sensor=FALSE;
      // now check
      if(LOCAL_is_sensor) { // need sensors package
        vector<string> vline_temp1,vline_temp2,tokens,tokens2; string sj=" ";
        aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("sensors")+" | grep temp1 | grep -v low | grep high | head -8 "),vline_temp1);
        aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("sensors")+" | grep temp2 | grep -v low | grep high | head -8 "),vline_temp2);
        for(uint i=0;i<vline_temp1.size();i++) {
          aurostd::StringSubst(vline_temp1.at(i),"temp1"," ");
          aurostd::StringSubst(vline_temp1.at(i),":"," ");aurostd::StringSubst(vline_temp1.at(i),"="," ");
          aurostd::StringSubst(vline_temp1.at(i),"("," ");aurostd::StringSubst(vline_temp1.at(i),")"," ");
          //aurostd::StringSubst(vline_temp1.at(i),""," ");aurostd::StringSubst(vline_temp1.at(i),"C"," "); //CO, PREVIOUSLY DEGREE SYMBOL, removed to get rid of warnings on mac
          aurostd::StringSubst(vline_temp1.at(i),"\u00B0"," ");aurostd::StringSubst(vline_temp1.at(i),"C"," "); //\u00B0:  http://www.fileformat.info/info/unicode/char/b0/index.htm
          aurostd::StringSubst(vline_temp1.at(i),"+"," ");aurostd::StringSubst(vline_temp1.at(i),"-"," ");
          for(unsigned char j=1;j<255;j++)
            if((j>=1 && j<=45 && j!=32) || (j==47) || (j>=58 && j<255)) {sj[0]=j;aurostd::StringSubst(vline_temp1.at(i),sj," ");}
          //   cerr << vline_temp1.at(i) << endl;
          aurostd::string2tokens(vline_temp1.at(i),tokens," ");
          // cerr << vline_temp1.at(i) << endl;
          for(uint j=0;j<tokens.size();j++) {
            if(aurostd::substring2bool(tokens.at(j),".")) {
              // aurostd::string2tokens(tokens.at(j),tokens2,".");
              tokens2=tokens;
              if(tokens2.size()>0) {XHOST.vTemperatureCore.push_back(aurostd::string2utype<double>(tokens2.at(0)));}
              break;
            }
          }
          // check for combinations
        }
        //  while(vline_temp2.size()>0) {vline_temp1.pop_back();vline_temp2.pop_back();} // remove the temp2&temp1 stuff
        if(vline_temp1.size()==9) vline_temp1.pop_back();
      } else {
        XHOST.vTemperatureCore.clear();
      }
    }
    pthread_mutex_unlock(&mutex_INIT_GetTEMP);
    return XHOST.vTemperatureCore.size();
  }
} // namespace init

// ***************************************************************************
// uint init::GetTEMP(void) // need sensors package
// ***************************************************************************
namespace init {
  double WaitTEMP(double TRESHOLD,ostream& oss,bool LVERBOSE,vector<string> vmessage) {
    bool _tmp_=XHOST.sensors_allowed;
    XHOST.sensors_allowed=TRUE;
    string message_PRE="",message_POST="";

    if(vmessage.size()>0) message_PRE=vmessage.at(0);
    if(vmessage.size()>1) message_POST=vmessage.at(1);

    stringstream sss;
    sss.setf(std::ios::fixed,std::ios::floatfield);
    sss.precision(1);

    init::GetTEMP();
    while (aurostd::max(XHOST.vTemperatureCore)>TRESHOLD) {
      int sleep=20+aurostd::abs(100.0*(TRESHOLD-aurostd::max(XHOST.vTemperatureCore))*aurostd::ran0());
      sss.clear();sss.str("");sss << aurostd::max(XHOST.vTemperatureCore) << " >  " << TRESHOLD;
      if(LVERBOSE) { oss << message_PRE << "init::WaitTEMP: max(TEMP) " << sss.str() << "   ... waiting, sleeping secs = " << sleep << message_POST << endl;oss.flush();}
      aurostd::execute("sleep "+aurostd::utype2string<int>(sleep));
      init::GetTEMP();
    }
    sss.clear();sss.str("");sss << aurostd::max(XHOST.vTemperatureCore) << " <= " << TRESHOLD;
    if(LVERBOSE) { oss << message_PRE << "init::WaitTEMP: max(TEMP) " << sss.str() << message_POST << endl;oss.flush();}
    XHOST.sensors_allowed=_tmp_;
    return aurostd::max(XHOST.vTemperatureCore);
  }
}


// ***************************************************************************
// AFLOW_getTEMP
// ***************************************************************************
uint AFLOW_getTEMP(const vector<string>& argv) {
  bool LDEBUG=(TRUE || XHOST.DEBUG);
  bool RUNBAR=aurostd::args2flag(argv,"--runbar|--RUNBAR|--runBAR|--bar|--BAR");
  bool RUNSTAT=aurostd::args2flag(argv,"--runstat|--RUNSTAT|--runSTAT|--stat|--STAT");
  string WRITE=aurostd::args2attachedstring(argv,"--write=","");
  double maxmem=aurostd::args2attachedutype<double>(argv,"--mem=|--maxmem=",XHOST.maxmem);
  double refresh=aurostd::args2attachedutype<double>(argv,"--refresh=",AFLOW_CORE_TEMPERATURE_REFRESH);
  double warning_beep=aurostd::args2attachedutype<double>(argv,"--warning_beep=",AFLOW_CORE_TEMPERATURE_BEEP);
  double warning_halt=aurostd::args2attachedutype<double>(argv,"--warning_halt=",AFLOW_CORE_TEMPERATURE_HALT);

  if(LDEBUG) cerr << "AFLOW_getTEMP: RUNBAR=" << RUNBAR << endl;
  if(LDEBUG) cerr << "AFLOW_getTEMP: RUNSTAT=" << RUNSTAT << endl;
  if(LDEBUG) cerr << "AFLOW_getTEMP: write=" << WRITE << endl;
  if(LDEBUG) cerr << "AFLOW_getTEMP: maxmem=" << maxmem << endl;
  if(LDEBUG) cerr << "AFLOW_getTEMP: refresh=" << refresh << "   -   AFLOW_CORE_TEMPERATURE_REFRESH=" << AFLOW_CORE_TEMPERATURE_REFRESH << endl;
  if(LDEBUG) cerr << "AFLOW_getTEMP: warning_beep=" << warning_beep << "   -   AFLOW_CORE_TEMPERATURE_BEEP=" << AFLOW_CORE_TEMPERATURE_BEEP << endl;
  if(LDEBUG) cerr << "AFLOW_getTEMP: warning_halt=" << warning_halt << "   -   AFLOW_CORE_TEMPERATURE_HALT=" << AFLOW_CORE_TEMPERATURE_HALT << endl;
  if(WRITE!="") aurostd::RemoveFile(WRITE);

  while(init::GetTEMP()) {
    stringstream oss;
    double Tmax=aurostd::max(XHOST.vTemperatureCore);
    double Tmin=aurostd::min(XHOST.vTemperatureCore);
    double Tzero=30.0;
    oss << "00000  MESSAGE " << aurostd::get_time() << " ";
    if(RUNSTAT || (!RUNSTAT && !RUNBAR)) {
      string soss="- [temp(C)=";
      for(uint i=0;i<XHOST.vTemperatureCore.size();i++) {soss+=aurostd::utype2string(XHOST.vTemperatureCore.at(i),3)+(i<XHOST.vTemperatureCore.size()-1?",":"]");}
      oss << aurostd::PaddedPOST(soss,XHOST.vTemperatureCore.size()*5+11," ");
      soss=" - ["+aurostd::utype2string(Tmin,3)+","+aurostd::utype2string(Tmax,3)+"]";
      oss << aurostd::PaddedPOST(soss,14," ") << "  beep=" << warning_beep << "  halt=" << warning_halt << " ";
    }
    if(RUNBAR) {
      for(double i=0;i<2*(Tmax-30.0);i++) oss << "*";
      oss << " " << Tmax;
    }
    if(maxmem<100.0) oss << " - [mem=" << aurostd::utype2string<double>(AFLOW_checkMEMORY("vasp",maxmem),4) << " (" << maxmem << ")]";
    if(Tmax>=warning_beep) oss << "  (beep) MAX>" << warning_beep;// << endl;
    if(Tmax>=warning_halt) oss << "  (halt) SHUTDOWN>" << warning_halt;// << endl;

    if(WRITE!="") {
      stringstream aus;vector<string> vlines;
      //      if(aurostd::FileExist(WRITE))
      {aurostd::file2vectorstring(WRITE,vlines);}
      aus<<"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">"<<endl;
      aus<<"<html> <head>"<<endl;
      aus<<"<META HTTP-EQUIV=\"expires\" CONTENT=\"0\"> <META NAME=\"robots\" CONTENT=\"none\"> <META NAME=\"robots\" CONTENT=\"noindex,nofollow\"> <META NAME=\"robots\" CONTENT=\"noarchive\">"<<endl;
      aus<<"<?php $page=$_SERVER['PHP_SELF']; $sec=\"" << int(refresh) << "\"; header(\"Refresh: $sec; url=$page\"); ?>"<<endl;
      aus << "</head> <!?php print strftime('%c'); ?> <pre>"<<endl;
      aus <<  oss.str() << endl;
      for(uint i=0;i<vlines.size();i++)
        if(i>4 && i<vlines.size()-1) aus << vlines.at(i) << endl;
      aus << "</pre> </body> </html>"<<endl;
      aurostd::stringstream2file(aus,WRITE);
    }

    if(Tmax>=warning_beep) { aurostd::execute(XHOST.command("beep")+" -l 100 -f "+aurostd::utype2string<double>(50*(Tmax-Tzero)));}
    if(Tmax>=warning_halt) {
      aurostd::execute(XHOST.command("beep")+" -f 1000");aurostd::execute(XHOST.command("beep")+" -f 1500");
      aurostd::execute(XHOST.command("halt"));aurostd::execute(XHOST.command("beep")+" -f 2000");
      aurostd::execute(XHOST.command("beep")+" -f 2500");}

    //   if(maxmem>0.0 && maxmem<100.0) oss << aurostd::utype2string<double>(AFLOW_checkMEMORY("vasp",maxmem),4);
    //   if(maxmem>0.0 && maxmem<100.0) AFLOW_checkMEMORY("vasp",maxmem);
    //    if(maxmem>0.0 && maxmem<100.0) AFLOW_checkMEMORY("aflow",maxmem);
    // if(maxmem>0.0 && maxmem<100.0) AFLOW_checkMEMORY("clamd",maxmem);

    cout << oss.str();// cerr << oss.str();
    oss.str(std::string());
    cout << endl;cout.flush();//cerr << endl;cerr.flush();

    if(!RUNSTAT && !RUNBAR) break;
    sleep(refresh*(1.0-(Tmax-Tzero)/40.0));
  }
  return 1;
}

// ***************************************************************************
// AFLOW_monitor
// ***************************************************************************
uint AFLOW_monitor(const vector<string>& argv) {
  cout << "MMMMM  Aflow: starting AFLOW_monitor" << endl;
  cerr << "MMMMM  Aflow: starting AFLOW_monitor" << endl;

  // double maxmem=aurostd::args2attachedutype<double>(argv,"--mem=","--maxmem=",double (95.0/XHOST.CPU_Cores));
  double maxmem=aurostd::args2attachedutype<double>(argv,"--mem=|--maxmem=",double (3.0));
  maxmem=round(10*maxmem)/10; // to get easy numbers.

  vector<string> argv_local;
  argv_local.push_back(argv.at(0));
  argv_local.push_back("--runstat");
  argv_local.push_back("--refresh=180");
  argv_local.push_back("--mem="+aurostd::utype2string<double>(maxmem));
  argv_local.push_back("--warning_beep=57");

  //  argv_local.push_back();

  return AFLOW_getTEMP(argv_local);
}

// ***************************************************************************
// CheckAFLOWLIBMaterialServer
// ***************************************************************************
bool CheckMaterialServer(void) {return CheckMaterialServer("");}
bool CheckMaterialServer(const string& message) { //CO20200624
  if(XHOST.hostname==XHOST.AFLOW_MATERIALS_SERVER) return TRUE;
  if(XHOST.hostname==XHOST.AFLOW_WEB_SERVER) return TRUE;
  if(XHOST.hostname=="habana") return TRUE;
  if(XHOST.hostname=="aflowlib") return TRUE;
  stringstream messagestream;
  messagestream << "Your machine is \"" << XHOST.hostname << "\". ";
  if(message.length()>0) messagestream << "Command \"" << message << "\" can run only on \"" << XHOST.AFLOW_MATERIALS_SERVER << "\" or \"" << XHOST.AFLOW_WEB_SERVER << "\"." << endl;
  else messagestream << "The procedure can run only on \"" << XHOST.AFLOW_MATERIALS_SERVER << "\" or \"" << XHOST.AFLOW_WEB_SERVER << "\".";
  throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, messagestream, _RUNTIME_ERROR_);
  return FALSE;
}

// ***************************************************************************
// aflow_get_time_string
// ***************************************************************************
string aflow_get_time_string(void) {
  //OUTPUT: http://www.cplusplus.com/reference/ctime/ctime/
  //Www Mmm dd hh:mm:ss yyyy
  //Where Www is the weekday, Mmm the month (in letters), dd the day of the month, hh:mm:ss the time, and yyyy the year.
  //The string is followed by a new-line character ('\n') and terminated with a null-character.
  //[CO20200624 - OBSOLETE]#ifdef ALPHA
  //[CO20200624 - OBSOLETE]  ostringstream aus;
  //[CO20200624 - OBSOLETE]  string OUT;
  //[CO20200624 - OBSOLETE]  aus<<"date | sed \"s/ /_/g\" > "+XHOST.tmpfs+"/date."<< XHOST.ostrPID.str() << "." << XHOST.ostrTID.str() << " " <<endl;  //CO20200502 - threadID
  //[CO20200624 - OBSOLETE]  system(aus.str().c_str());
  //[CO20200624 - OBSOLETE]  ifstream FileAUS;
  //[CO20200624 - OBSOLETE]  string FileNameAUS=XHOST.tmpfs+"/date."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str();  //CO20200502 - threadID
  //[CO20200624 - OBSOLETE]  FileAUS.open(FileNameAUS.c_str(),std::ios::in);
  //[CO20200624 - OBSOLETE]  FileAUS >> OUT;
  //[CO20200624 - OBSOLETE]  FileAUS.clear();FileAUS.close();
  //[CO20200624 - OBSOLETE]  // return (char*) OUT.c_str();
  //[CO20200624 - OBSOLETE]  return string("NotAvailable \n");
  //[CO20200624 - OBSOLETE]#else
  long ltime=time(NULL);

  string date=string(ctime(&ltime));
  if(date.length()>0)
    if(date.at(date.length()-1)=='\n')
      date.erase(date.length()-1);
  return date;
  //[CO20200624 - OBSOLETE]#endif
}
string aflow_convert_time_ctime2aurostd(const string& time_LOCK){ //CO20200624
  //refer to aflow_get_time_string()
  //convert Www Mmm dd hh:mm:ss yyyy style to aurostd::get_datetime() one
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  vector<string> tokens;
  aurostd::string2tokens(time_LOCK,tokens);
  if(tokens.size()!=5){return "";}

  //https://en.cppreference.com/w/c/chrono/strftime
  //'Www Mmm dd hh:mm:ss yyyy' === '%a %b %d %H:%M:%S %Y'
  //https://stackoverflow.com/questions/19524720/using-strptime-converting-string-to-time-but-getting-garbage
  tm tstruct;
  if(!strptime(time_LOCK.c_str(),"%a %b %d %H:%M:%S %Y",&tstruct)){return "";}
  tstruct.tm_isdst=-1;  //let computer figure it out
  std::mktime(&tstruct);  //get is_dst

  if(LDEBUG){
    char buffer[30];
    strftime(buffer,30,"%F %T %Z",&tstruct);
    cerr << __AFLOW_FUNC__ << " tstruct=" << buffer << endl;
  }

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}

  bool include_utc_offset=true;
  return aurostd::get_datetime(tstruct,include_utc_offset);
}

// ***************************************************************************
// aflow_get_time_string_short
// ***************************************************************************
string aflow_get_time_string_short(void) {
  string date;
  vector<string> tokens;
  aurostd::string2tokens(aflow_get_time_string(),tokens);
  date="na";
  if(tokens.size()>4) {
    if(tokens.at(2).size()>1)
      date=tokens.at(4).substr(2,2)+tokens.at(1)+tokens.at(2);
    else
      date=tokens.at(4).substr(2,2)+tokens.at(1)+"0"+tokens.at(2);
    aurostd::StringSubst(date,"Jan","01");aurostd::StringSubst(date,"Feb","02");aurostd::StringSubst(date,"Mar","03");
    aurostd::StringSubst(date,"Apr","04");aurostd::StringSubst(date,"May","05");aurostd::StringSubst(date,"Jun","06");
    aurostd::StringSubst(date,"Jul","07");aurostd::StringSubst(date,"Aug","08");aurostd::StringSubst(date,"Sep","09");
    aurostd::StringSubst(date,"Oct","10");aurostd::StringSubst(date,"Nov","11");aurostd::StringSubst(date,"Dec","12");
    date=date.substr(0,6);
  }
  if(date.length()>0)
    if(date.at(date.length()-1)=='\n')
      date.erase(date.length()-1);
  return date;
}

// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] // strPID
// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] string strPID(void) {
// [OBSOLETE]   int PID=getpid();
// [OBSOLETE]   ostringstream oss;
// [OBSOLETE]   oss << PID;
// [OBSOLETE]   return (string) oss.str();
// [OBSOLETE] }

// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] // strTID
// [OBSOLETE] // ***************************************************************************
// [OBSOLETE] string strTID(void) { //CO20200502 - threadID
// [OBSOLETE]   int TID=aurostd::getTID();
// [OBSOLETE]   ostringstream oss;
// [OBSOLETE]   oss << TID;
// [OBSOLETE]   return (string) oss.str();
// [OBSOLETE] }

// ***************************************************************************
// Messages
// ***************************************************************************
double AFLOW_checkMEMORY(const string& progname,double memory) {
  vector<string> vps,tokens;string command;
  double maxmem=0.0;
  if(progname.empty()) aurostd::string2vectorstring(aurostd::execute2string("ps aux | grep -v \" 0.0  0.0 \" | grep "+XHOST.user),vps);
  else aurostd::string2vectorstring(aurostd::execute2string("ps aux | grep \""+progname+"\" | grep -v \" 0.0  0.0 \" | grep "+XHOST.user),vps);
  for(uint i=0;i<vps.size();i++) {
    aurostd::string2tokens(vps.at(i),tokens);
    if(tokens.size()>4) {
      if(aurostd::string2utype<double>(tokens.at(3))>maxmem) maxmem=aurostd::string2utype<double>(tokens.at(3));
      if(memory>0.0 && memory<100.0) {
        if(aurostd::string2utype<double>(tokens.at(3))>memory) {
          command=string(XHOST.command("kill")+" -9 "+tokens.at(1));
          aurostd::execute(command);
          //	  cerr << endl << "AFLOW_checkMEMORY: killing(" << memory << ") = " << vps.at(i) << endl;
          cout << endl << "AFLOW_checkMEMORY [date=" << aflow_get_time_string() << "]: kill(" << tokens.at(3) << ">" << aurostd::utype2string<double>(memory,4) << ") = [" << vps.at(i) << "]" << endl;
        }
      }
    }
  }
  return maxmem;
}

// ***************************************************************************
// GetVASPBinaryFromLOCK
// ***************************************************************************
bool GetVASPBinaryFromLOCK(const string& directory,string& vasp_bin){  //CO20210315
  int ncpus=0;
  return GetVASPBinaryFromLOCK(directory,vasp_bin,ncpus);
}
bool GetVASPBinaryFromLOCK(const string& directory,string& vasp_bin,int& ncpus){  //CO20210315
  bool LDEBUG=(FALSE || XHOST.DEBUG);

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  //reset
  vasp_bin="";
  ncpus=0;

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking for " << directory+"/"+_AFLOWLOCK_ << endl;}
  if(!aurostd::FileExist(directory+"/"+_AFLOWLOCK_)){return false;}
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " FOUND " << directory+"/"+_AFLOWLOCK_ << endl;}

  uint i=0,j=0,vlines_size=0,vtokens_size=0;
  vector<string> vlines,vtokens;
  vlines_size=aurostd::file2vectorstring(directory+"/"+_AFLOWLOCK_,vlines);
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << directory+"/"+_AFLOWLOCK_ << " contains " << vlines_size << " lines" << endl;}

  for(i=vlines_size-1;i<vlines_size;i--){ //go backwards
    if(vlines[i].find(VASP_KEYWORD_EXECUTION)==string::npos){continue;} //look for 'Executing:' line
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found line containing execution command: \"" << vlines[i] << "\"" << endl;}
    vtokens_size=aurostd::string2tokens(vlines[i],vtokens," ");
    for(j=2;j<vtokens_size;j++){
      if(vtokens[j]==DEFAULT_VASP_OUT||vtokens[j]==DEFAULT_VASP_OUT+";"){ //looking for sequence VASP_BIN >> VASP_OUT //if j goes below 0, then it goes to max uint //CO20210315 - adding semicolon as NO HOST mode used to print it
        vasp_bin=vtokens[j-2];
        ncpus=1;  //set default
        if(j>2 && aurostd::isfloat(vtokens[j-3])){ncpus=aurostd::string2utype<int>(vtokens[j-3]);} //grab if available, it will be just before the bin
        return true;
      }
    }
  }
  return false;
}

// ***************************************************************************
// processFlagsFromLOCK
// ***************************************************************************
void processFlagsFromLOCK(_xvasp& xvasp,_vflags& vflags,aurostd::xoption& xfixed){  //CO20210315
  bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || XHOST.DEBUG);

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  if(!aurostd::FileExist(xvasp.Directory+"/"+_AFLOWLOCK_)){return;}

  uint i=0;

  //clear out past options
  //[needs to be targeted, keep FLAG::AFIX_DRYRUN]xvasp.aopts.clear();
  vector<string> flags;
  flags=xvasp.aopts.vxscheme;
  for(i=0;i<flags.size();i++){ //capture all _PRESERVED flags (KPOINTS, POSCAR, etc.)
    const string& flag=flags[i];
    if(flag.find("FLAG::")!=string::npos && flag.find("_PRESERVED")!=string::npos){xvasp.aopts.flag(flag,false);}
  }
  vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.clear();
  vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved=false;
  xfixed.clear();

  uint vlines_size=0;
  vector<string> vlines,vtokens;
  vlines_size=aurostd::file2vectorstring(xvasp.Directory+"/"+_AFLOWLOCK_,vlines);

  string str_xvasp_aopts_start="xvasp.aopts.flag(\"";
  string str_vflags_ignore_afix_start="vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"";
  string str_xfixed_start="xfixed.flag(\"";
  string str_flag_end="\")=1";
  string scheme="";

  string::size_type loc_start=0,loc_end=0;

  uint nexecuting=0;
  for(i=vlines_size-1;i<vlines_size;i--){ //go backwards
    const string& line=vlines[i];
    if(line.find(VASP_KEYWORD_EXECUTION)!=string::npos){nexecuting++;} //look for 'Executing:' line
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " vlines[i=" << i << "]=" << line << endl;
      cerr << __AFLOW_FUNC__ << " nexecuting=" << nexecuting << endl;
    }
    if(nexecuting>1){break;}

    loc_start=line.find(str_xvasp_aopts_start);
    loc_end=line.find(str_flag_end);
    if(loc_start!=string::npos && loc_end!=string::npos){
      scheme=line.substr(loc_start+str_xvasp_aopts_start.length(),loc_end-(loc_start+str_xvasp_aopts_start.length()));
      xvasp.aopts.flag(scheme,true);
    }

    loc_start=line.find(str_vflags_ignore_afix_start);
    loc_end=line.find(str_flag_end);
    if(loc_start!=string::npos && loc_end!=string::npos){
      scheme=line.substr(loc_start+str_vflags_ignore_afix_start.length(),loc_end-(loc_start+str_vflags_ignore_afix_start.length()));
      vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(scheme,true);
    }

    if(line.find("vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved=1")!=string::npos){vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved=true;}

    loc_start=line.find(str_xfixed_start);
    loc_end=line.find(str_flag_end);
    if(loc_start!=string::npos && loc_end!=string::npos){
      scheme=line.substr(loc_start+str_xfixed_start.length(),loc_end-(loc_start+str_xfixed_start.length()));
      xfixed.flag(scheme,true);
    }

  }

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " printing xvasp.aopts" << endl;
    for(uint i=0;i<xvasp.aopts.vxscheme.size();i++){cerr << __AFLOW_FUNC__ << " xvasp.aopts.flag(\"" << xvasp.aopts.vxscheme[i] << "\")=" << xvasp.aopts.flag(xvasp.aopts.vxscheme[i]) << endl;}
    cerr << __AFLOW_FUNC__ << " printing xfixed" << endl;
    for(uint i=0;i<xfixed.vxscheme.size();i++){cerr << __AFLOW_FUNC__ << " xfixed.flag(\"" << xfixed.vxscheme[i] << "\")=" << xfixed.flag(xfixed.vxscheme[i]) << endl;}
  }
}

// ***************************************************************************
// AFLOW_VASP_instance_running
// ***************************************************************************
bool AFLOW_VASP_instance_running(){ //CO20210315
  //[SD20220402 - OBSOLETE]//this needs to become more complicated as we add options other than --kill_vasp_all
  //[SD20220402 - OBSOLETE]if(XHOST.vflag_control.flag("MONITOR_VASP") && XHOST.vflag_control.flag("KILL_VASP_ALL")==false){
  //[SD20220402 - OBSOLETE]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"a targeted kill command for VASP not available yet: try --kill_vasp_all",_INPUT_ILLEGAL_);
  //[SD20220402 - OBSOLETE]}
  return (aurostd::ProcessPIDs("aflow").size()>1);  //check that the instance of aflow running vasp is running
}
bool AFLOW_VASP_instance_running(const string& vasp_pgid){
  //SD20220329 - only check instances of aflow with the same PGID as the monitor
  string output_syscall="";
  return (aurostd::ProcessPIDs("aflow",vasp_pgid,output_syscall).size()>1);  //check that the instance of aflow running vasp is running
}


// ***************************************************************************
// AFLOW_MONITOR_instance_running
// ***************************************************************************
bool AFLOW_MONITOR_instance_running(const _aflags& aflags){ //CO20210315
  //[SD20220402 - OBSOLETE]//this needs to become more complicated as we add options other than --kill_vasp_all
  //[SD20220402 - OBSOLETE]if(XHOST.vflag_control.flag("MONITOR_VASP") && XHOST.vflag_control.flag("KILL_VASP_ALL")==false){ //CO20210315 - this check doesn't really apply, we wouldn't call this function if we were running with --monitor_vasp, but it's a good reminder that this code needs to become smarter in the future with flags other than --kill_vasp_all
  //[SD20220402 - OBSOLETE]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"a targeted kill command for VASP not available yet: try --kill_vasp_all",_INPUT_ILLEGAL_);
  //[SD20220402 - OBSOLETE]}
  return aurostd::FileExist(aflags.Directory+"/"+_AFLOWLOCK_+"."+FILE_VASP_MONITOR);
}

// ***************************************************************************
// VASP_instance_running
// ***************************************************************************
bool VASP_instance_running(const string& vasp_bin){ //CO20210315
  //[SD20220402 - OBSOLETE]//this needs to become more complicated as we add options other than --kill_vasp_all
  //[SD20220402 - OBSOLETE]if(XHOST.vflag_control.flag("MONITOR_VASP") && XHOST.vflag_control.flag("KILL_VASP_ALL")==false){
  //[SD20220402 - OBSOLETE]  throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"a targeted kill command for VASP not available yet: try --kill_vasp_all",_INPUT_ILLEGAL_);
  //[SD20220402 - OBSOLETE]}
  return aurostd::ProcessRunning(vasp_bin);
}
bool VASP_instance_running(const string& vasp_bin,const string& vasp_pgid){
  //SD20220329 - only check instances of vasp with the same PGID as the monitor
  return aurostd::ProcessRunning(vasp_bin,vasp_pgid);
}

// ***************************************************************************
// AFLOW_monitor_VASP
// ***************************************************************************
#define NCOUNTS_WAIT_MONITOR 10 //wait no more than 10*sleep_seconds (should be 10 minutes)
void AFLOW_monitor_VASP(){  //CO20210601
  //AFLOW_monitor_VASP() with no input arguments will look for FILE/DIRECTORY input from XHOST (--FILE or --D)
  //then it will pass the path to AFLOW_monitor_VASP(const string& directory)
  //this function will keep reading the --FILE input for new directories
  bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || XHOST.DEBUG);

  if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

  string directory="";

  if(XHOST.vflag_control.flag("FILE")){
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found FILE input: " << XHOST.vflag_control.getattachedscheme("FILE") << endl;}
    uint nfailures=0,nsuccesses=0;  //nfailures waits for file to be written, nsuccesses makes sure we didn't run out of directories
    string file=XHOST.vflag_control.getattachedscheme("FILE");
    string file_dir=aurostd::dirname(file);
    vector<string> vlines,vtokens;
    uint i=0,j=0,vlines_size=0,vtokens_size=0;
    bool found=false;
    string search_str="[dir=";
    string directory_old="";
    while(nfailures<NCOUNTS_WAIT_MONITOR && nsuccesses<NCOUNTS_WAIT_MONITOR){  //10 minutes, keeps us from wasting walltime
      if(!aurostd::FileExist(file)){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " FILE not yet found, waiting...: FILE=" << file << endl;}
        nfailures++;  //increment, only wait 10 minutes for a file to be written
        aurostd::Sleep(SECONDS_SLEEP_VASP_MONITOR);
        continue;
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " FILE found: FILE=" << file << endl;}
      nfailures=0;  //reset
      vlines_size=aurostd::file2vectorstring(file,vlines);
      found=false;
      for(i=vlines_size-1;i<vlines_size&&!found;i--){ //go backwards
        //looking for "... - [dir=.] - [user=aflow] ..."
        if(vlines[i].find(search_str)==string::npos){continue;}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found line containing \"" << search_str << "\": " << vlines[i] << endl;}
        vtokens_size=aurostd::string2tokens(vlines[i],vtokens," ");
        for(j=0;j<vtokens_size&&!found;j++){
          if(vtokens[j].find(search_str)==string::npos){continue;}
          directory=vtokens[j];
          aurostd::StringSubst(directory,search_str,"");aurostd::StringSubst(directory,"]","");
          found=true;
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " dir=" << directory << endl;}
          if(directory==directory_old){ //increment, only wait 10 minutes for a new directory to be found
            if(LDEBUG){cerr << __AFLOW_FUNC__ << " already ran this directory, waiting...: dir=" << directory << endl;}
            nsuccesses++;
            continue;
          }
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " found directory to run: dir=" << directory << endl;}
          directory_old=directory;
          nsuccesses=0;
          AFLOW_monitor_VASP(directory);
        }
      }
      //if --FILE=LOCK, this will be useful
      if(aurostd::EFileExist(file_dir+"/"+DEFAULT_AFLOW_END_OUT)){break;}
      if(aurostd::EFileExist(file_dir+"/"+_STOPFLOW_)){aurostd::RemoveFile(file_dir+"/"+_STOPFLOW_);break;}

      aurostd::Sleep(SECONDS_SLEEP_VASP_MONITOR);

      //if --FILE=LOCK, this will be useful
      if(aurostd::EFileExist(file_dir+"/"+DEFAULT_AFLOW_END_OUT)){break;}
      if(aurostd::EFileExist(file_dir+"/"+_STOPFLOW_)){aurostd::RemoveFile(file_dir+"/"+_STOPFLOW_);break;}
    }
  }
  else if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " found DIRECTORY input: " << XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN") << endl;}
    directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");
    return AFLOW_monitor_VASP(directory);
  }
  else{throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no directory input provided",_INPUT_ILLEGAL_);}
  if(LDEBUG){cerr << __AFLOW_FUNC__ << " END" << endl;}
}

void AFLOW_monitor_VASP(const string& directory){ //CO20210601
  //AFLOW_monitor_VASP() with no input arguments will look for FILE/DIRECTORY input from XHOST (--FILE or --D)
  //then it will pass the path to AFLOW_monitor_VASP(const string& directory)
  //this function will finish once the calculation in the particular directory finishes
  stringstream message;

  //aflags
  _aflags aflags;
  aflags.Directory=directory;
  if(directory.empty()){return;} //throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"directory input is empty",_INPUT_ILLEGAL_);
  if(!aurostd::FileExist(aflags.Directory+"/"+_AFLOWIN_)){return;}

  //output objects
  ofstream FileMESSAGE,FileMESSAGE_devnull;
  string FileNameLOCK=aflags.Directory+"/"+_AFLOWLOCK_+"."+FILE_VASP_MONITOR;
  FileMESSAGE.open(FileNameLOCK.c_str(),std::ios::out);
  bool oss_silent=true;
  ostream& oss=cout;if(oss_silent){oss.setstate(std::ios_base::badbit);}  //like NULL - cannot print to cout with two instances of aflow running

  message << aflow::Banner("BANNER_NORMAL");pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);

  //aflow.in
  string AflowIn_file="",AflowIn="";
  try{KBIN::getAflowInFromAFlags(aflags,AflowIn_file,AflowIn,FileMESSAGE,oss);}
  catch(aurostd::xerror& err){
    pflow::logger(err.whereFileName(), err.whereFunction(), err.what(), aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
    return;
  }

  //other flags
  _kflags kflags=KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,oss);
  _vflags vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags,oss);

  _xvasp xvasp;
  xvasp.Directory=aflags.Directory; //arun_directory;
  message << "START        - " << xvasp.Directory;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_); //include directory so noticeable space remains
  AVASP_populateXVASP(aflags,kflags,vflags,xvasp);
  xvasp.aopts.flag("FLAG::AFIX_DRYRUN",true); //VERY IMPORTANT

  //a note about xmonitor: it only controls the amount of output to the LOCK file, it does NOT affect performance

  int nloop=1;
  uint i=0;
  uint nexecuting=0,nexecuting_old=0;
  uint n_not_running=0,n_running=0;
  uint sleep_seconds=SECONDS_SLEEP_VASP_MONITOR;
  uint sleep_seconds_afterkill=sleep_seconds;
  aurostd::xoption xmessage,xwarning,xmonitor,xfixed;
  bool VERBOSE=(FALSE || VERBOSE_MONITOR_VASP);
  bool vasp_running=false;
  uint vlines_lock_size=0;
  vector<string> vlines_lock;

  //there are issues getting the correct vasp binary since this is an entirely different aflow instance
  //we might run the other aflow instance with --mpi or --machine flags that affect which vasp binary we use
  //therefore, the most robust way to define the binary is to search the LOCK file
  //[CO20210315 - OBSOLETE]string& vasp_bin=kflags.KBIN_MPI_BIN;
  //[CO20210315 - OBSOLETE]if(!(kflags.KBIN_MPI==true||XHOST.MPI==true)){vasp_bin=kflags.KBIN_BIN;}

  //CO20210315 - when we generalize this code to run for targeted instances of vasp, we also need to target the right instances of aflow

  string vasp_bin="",vasp_pgid=XHOST.ostrPGID.str();
  int ncpus=0;
  nloop=0;
  string memory_string="";
  double usage_percentage_ram=0.0,usage_percentage_swap=0.0;

  if(vasp_pgid.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no vasp PGID found",_RUNTIME_ERROR_);}
  if(VERBOSE){message << "vasp_pgid=\"" << vasp_pgid << "\" (from XHOST)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}
  while((AFLOW_VASP_instance_running(vasp_pgid) || (nloop++)<NCOUNTS_WAIT_MONITOR) && vasp_bin.empty()){  //wait no more than 10 minutes for vasp bin to start up
    GetVASPBinaryFromLOCK(xvasp.Directory,vasp_bin,ncpus);
    vasp_bin=aurostd::basename(vasp_bin); //remove directory stuff
    if(vasp_bin.empty()){
      if(VERBOSE){message << "sleeping for " << sleep_seconds_afterkill << " seconds, waiting for VASP binary to start running and " << _AFLOWLOCK_ << " to be written";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}
      aurostd::Sleep(sleep_seconds_afterkill); //sleep at least a minute to let aflow start up
      continue;
    }
  }
  if(vasp_bin.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"no vasp binary found in "+xvasp.Directory+"/"+_AFLOWLOCK_,_RUNTIME_ERROR_);}
  message << "vasp_binary=\"" << vasp_bin << "\" (from " << _AFLOWLOCK_ << ")";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  //SD20220602 - check that an OUTCAR is produced on the initial launch of vasp
  if(XHOST.vflag_control.flag("AFLOW_STARTUP_SCRIPT") && !aurostd::FileExist(xvasp.Directory+"/OUTCAR")){
    string startup_script_name=XHOST.vflag_control.getattachedscheme("AFLOW_STARTUP_SCRIPT");
    aurostd::ProcessKill(startup_script_name,vasp_pgid,true,5); //send SIGTRAP rather than SIGKILL so we can trap it
    while (aurostd::ProcessRunning(startup_script_name,vasp_pgid) && (n_running++)<NCOUNTS_WAIT_MONITOR){ //try killing the script for no more than 10 minutes
      message << "initial launch of vasp did not produce an OUTCAR file; killing the whole AFLOW startup script \"" << startup_script_name << "\"";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      aurostd::ProcessKill(startup_script_name,vasp_pgid); //send SIGKILL
      aurostd::Sleep(sleep_seconds_afterkill); //sleep to wait for the script to get killed
    }
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Initial launch of vasp did not produce an OUTCAR file, tried killing the whole AFLOW startup script",_RUNTIME_ERROR_);
  }

  nloop=0;
  n_not_running=0;
  while(AFLOW_VASP_instance_running(vasp_pgid) || n_not_running<NCOUNTS_WAIT_MONITOR){  //wait no more than 10 minutes for vasp bin to start up (again)
    if(!aurostd::FileExist(xvasp.Directory+"/"+_AFLOWLOCK_)){break;} //we needed it above to get the vasp_bin
    if(aurostd::EFileExist(xvasp.Directory+"/"+DEFAULT_AFLOW_END_OUT)){break;}  //check before continue below
    if(aurostd::EFileExist(xvasp.Directory+"/"+_STOPFLOW_)){aurostd::RemoveFile(xvasp.Directory+"/"+_STOPFLOW_);break;} //check before continue below

    vasp_running=VASP_instance_running(vasp_bin,vasp_pgid);
    if(VERBOSE){
      message << "nloop=" << (nloop++);pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      message << "program \"" << vasp_bin << "\" is " << (vasp_running?"":"NOT ") << "running";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    }
    if(vasp_running==false){
      n_not_running++;
      if(VERBOSE){message << "sleeping for " << sleep_seconds_afterkill << " seconds, waiting for \"" << vasp_bin << "\" to start running";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}
      aurostd::Sleep(sleep_seconds_afterkill); //sleep at least a minute to let aflow sleep since OUTCAR is incomplete
      continue;
    }
    n_not_running=0;  //reset

    //determine whether we need to clear xmonitor with new vasp instance (relax1->relax2)
    nexecuting=0;
    vlines_lock_size=aurostd::file2vectorstring(xvasp.Directory+"/"+_AFLOWLOCK_,vlines_lock);  //we already checked above that it exists
    for(i=0;i<vlines_lock_size;i++){
      if(vlines_lock[i].find(VASP_KEYWORD_EXECUTION)==string::npos){continue;}
      nexecuting++;
    }

    //SD20220602 - check that an OUTCAR is produced on the new vasp instance
    if(XHOST.vflag_control.flag("AFLOW_STARTUP_SCRIPT") && !aurostd::FileExist(xvasp.Directory+"/OUTCAR")){
      string startup_script_name=XHOST.vflag_control.getattachedscheme("AFLOW_STARTUP_SCRIPT");
      aurostd::ProcessKill(startup_script_name,vasp_pgid,true,5); //send SIGTRAP rather than SIGKILL so we can trap it
      while (aurostd::ProcessRunning(startup_script_name,vasp_pgid) && (n_running++)<NCOUNTS_WAIT_MONITOR){ //try killing the script for no more than 10 minutes
        message << "new vasp instance did not produce an OUTCAR file, killing the whole AFLOW startup script \"" << startup_script_name << "\"";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
        aurostd::ProcessKill(startup_script_name,vasp_pgid); //send SIGKILL
        aurostd::Sleep(sleep_seconds_afterkill); //sleep to wait for the script to get killed
      }
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"New vasp instance did not produce an OUTCAR file, tried killing the whole AFLOW startup script",_RUNTIME_ERROR_);
    }

    if(nexecuting!=nexecuting_old){
      message << "found new instance of \"" << vasp_bin << "\"";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      xmonitor.clear(); //new run, clear IGNORE_WARNINGS //we've transitioned from relax1->relax2, etc.
      nexecuting_old=nexecuting;
    }

    //check vasp output file here
    KBIN::VASP_ProcessWarnings(xvasp,aflags,kflags,xmessage,xwarning,xmonitor,FileMESSAGE);

    //check memory again, it's possible it floated above the threshold only for a second
    usage_percentage_ram=0.0;usage_percentage_swap=0.0;
    if(0){  //do not turn off MEMORY because it fails GetMemory(), it's possible MEMORY was triggered for other reasons (e.g., FROZEN_CALC)
      if(xwarning.flag("MEMORY")){
        bool ignore_memory=false;
        if(aurostd::GetMemoryUsagePercentage(usage_percentage_ram,usage_percentage_swap)==false){
          message << "ignoring xwarning.flag(\"MEMORY\"), could not retrieve memory status";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
          ignore_memory=true;
        }else{
          if(usage_percentage_ram<MEMORY_MAX_USAGE_RAM && usage_percentage_swap<MEMORY_MAX_USAGE_SWAP){
            message << "ignoring xwarning.flag(\"MEMORY\"), memory usage dropped below threshold ("+aurostd::utype2string(usage_percentage_ram,2,FIXED_STREAM)+"% ram usage,"+aurostd::utype2string(usage_percentage_swap,2,FIXED_STREAM)+"% swap usage)";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
            ignore_memory=true;
          }
        }
        if(ignore_memory){xwarning.flag("MEMORY",false);}
      }
    }

    bool kill_vasp=false;
    if(xwarning.flag()){  //if any flag is on
      kill_vasp=true;
      //the --monitor_vasp instance will not have the right ncpus set, so grab it from the LOCK
      GetVASPBinaryFromLOCK(xvasp.Directory,vasp_bin,ncpus);  //grab ncpus and set to kflags.KBIN_MPI_NCPUS
      vasp_bin=aurostd::basename(vasp_bin); //remove directory stuff
      if(ncpus>0){kflags.KBIN_MPI_NCPUS=ncpus;}
      //HERE: plug in exceptions from xfixed, etc. to turn OFF kill_vasp
      //read LOCK to see what has been issued already
      processFlagsFromLOCK(xvasp,vflags,xfixed);
      xfixed.flag("ALL",false); //otherwise new attempts cannot be tried
      if(VERBOSE){message << "ls (pre)" << endl << aurostd::execute2string("ls -l") << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}
      if(!KBIN::VASP_FixErrors(xvasp,xmessage,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE)){
        if(xwarning.flag("CALC_FROZEN")==false && xwarning.flag("OUTPUT_LARGE")==false){  //kill vasp if the calc is frozen (no solution)
          kill_vasp=false;
          message << "kill_vasp=" << kill_vasp << ": all fixes exhausted" << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
        }
      }
      if(VERBOSE){message << "ls (post)" << endl << aurostd::execute2string("ls -l") << endl;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}
    }

    //it's possible KBIN::VASP_ProcessWarnings() takes a long time (big vasp.out)
    //double check that there isn't a new instance of vasp running
    nexecuting=0;
    vlines_lock_size=aurostd::file2vectorstring(xvasp.Directory+"/"+_AFLOWLOCK_,vlines_lock);  //we already checked above that it exists
    for(i=0;i<vlines_lock_size;i++){
      if(vlines_lock[i].find(VASP_KEYWORD_EXECUTION)==string::npos){continue;}
      nexecuting++;
    }
    if(nexecuting!=nexecuting_old){
      message << "found new instance of \"" << vasp_bin << "\"";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      kill_vasp=false;
      nexecuting_old=nexecuting;
    }

    if(kill_vasp){
      vasp_running=VASP_instance_running(vasp_bin,vasp_pgid);
      if(vasp_running){
        //special case for MEMORY, the error will be triggered in the --monitor_vasp instance, and not in the --run one
        //so write out "AFLOW ERROR: AFLOW_MEMORY" so it gets caught in the --run instance
        if(xwarning.flag("MEMORY")){
          memory_string=" "+string(AFLOW_MEMORY_TAG); //pre-pending space to match formatting
          if(aurostd::GetMemoryUsagePercentage(usage_percentage_ram,usage_percentage_swap)){
            memory_string+=" ("+aurostd::utype2string(usage_percentage_ram,2,FIXED_STREAM)+"% ram usage,"+aurostd::utype2string(usage_percentage_swap,2,FIXED_STREAM)+"% swap usage)";  //CO20210315 - use fixed stream here, since we'll have two sig figs before decimal, and two after (99.99%)
          }
          memory_string+="\n";
          aurostd::string2file(memory_string,xvasp.Directory+"/"+DEFAULT_VASP_OUT,"APPEND");
        }
        //write BEFORE issuing the kill, the other instance of aflow will start to act as soon as the process is dead
        message << "issuing kill command for: \"" << vasp_bin << "\"";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
        if(0){  //super debug
          string output_syscall="";
          vector<string> vpids=aurostd::ProcessPIDs(vasp_bin,vasp_pgid,output_syscall);
          message << "output_syscall=";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
          message << output_syscall;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);
          message << "PIDs2kill="+aurostd::joinWDelimiter(vpids,",");pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
        }
        aurostd::ProcessKill(vasp_bin,vasp_pgid); //send SIGKILL
      }else{
        message << "\"" << vasp_bin << "\" has died before the kill command could be issued";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      }
      xmonitor.clear(); //new run, clear IGNORE_WARNINGS
      message << "sleeping for " << sleep_seconds_afterkill << " seconds after kill command";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
      aurostd::Sleep(sleep_seconds_afterkill); //there are two aflow instances. the aflow-monitor should sleep at least a minute, as (in the worst case) aflow-vasp needs 15-30 seconds to sleep while it waits for an incomplete OUTCAR to finish writing. you don't want aflow-vasp to pick up AFTER aflow-monitor.
    }

    if(aurostd::EFileExist(xvasp.Directory+"/"+DEFAULT_AFLOW_END_OUT)){break;}
    if(aurostd::EFileExist(xvasp.Directory+"/"+_STOPFLOW_)){aurostd::RemoveFile(xvasp.Directory+"/"+_STOPFLOW_);break;}

    if(VERBOSE){message << "sleeping for " << sleep_seconds << " seconds";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);}
    aurostd::Sleep(sleep_seconds);

    if(aurostd::EFileExist(xvasp.Directory+"/"+DEFAULT_AFLOW_END_OUT)){break;}
    if(aurostd::EFileExist(xvasp.Directory+"/"+_STOPFLOW_)){aurostd::RemoveFile(xvasp.Directory+"/"+_STOPFLOW_);break;}
  }
  message << "END        - " << xvasp.Directory;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_); //include directory so noticeable space remains
  FileMESSAGE.flush();FileMESSAGE.clear();FileMESSAGE.close();
}

// ***************************************************************************
// Messages
// ***************************************************************************
pthread_mutex_t mutex_INIT_Message=PTHREAD_MUTEX_INITIALIZER;
string Message(const string& filename,const string& list2print){
  pthread_mutex_lock(&mutex_INIT_Message);
  // pthread_mutex_unlock(&mutex_INIT_Message);
  stringstream strout;
  string LIST2PRINT=aurostd::toupper(list2print); //CO+DX20200825
  if(aurostd::substring2bool(LIST2PRINT,"USER")) strout << " - [user=" << XHOST.user << "]";
  if(aurostd::substring2bool(LIST2PRINT,"GROUP")) strout << " - [group=" << XHOST.group << "]";
  if(aurostd::substring2bool(LIST2PRINT,"HOST") || aurostd::substring2bool(LIST2PRINT,"HOSTNAME") || aurostd::substring2bool(LIST2PRINT,"MACHINE")) strout << " - [host=" << XHOST.hostname << "]";
  if(aurostd::substring2bool(LIST2PRINT,"TEMPERATURE")) if(init::GetTEMP()) for(uint i=0;i<XHOST.vTemperatureCore.size();i++) {strout << (i==0?" - [temp(C)=":"") << XHOST.vTemperatureCore.at(i) << (i<XHOST.vTemperatureCore.size()-1?",":"]");}
  if(aurostd::substring2bool(LIST2PRINT,"PID")) strout << " - [PID=" << XHOST.PID << "]";  //CO20200502
  if(aurostd::substring2bool(LIST2PRINT,"PGID")) strout << " - [PGID=" << XHOST.PGID << "]";  //SD20220330
  if(aurostd::substring2bool(LIST2PRINT,"TID")) strout << " - [TID=" << XHOST.TID << "]";  //CO20200502
  if(LIST2PRINT.empty() || aurostd::substring2bool(LIST2PRINT,"TIME") || aurostd::substring2bool(LIST2PRINT,"DATE")) strout << " - [date=" << aflow_get_time_string() << "]";   //CO20200624
  if(aurostd::substring2bool(LIST2PRINT,"MEMORY") && (XHOST.maxmem>0.0 && XHOST.maxmem<100)) strout << " - [mem=" << aurostd::utype2string<double>(AFLOW_checkMEMORY("vasp",XHOST.maxmem),4) << " (" << XHOST.maxmem << ")]"; //CO20170628 - slow otherwise!!!
  if(XHOST.vTemperatureCore.size()>0) if(max(XHOST.vTemperatureCore)>AFLOW_CORE_TEMPERATURE_BEEP) strout << " - [ERROR_TEMPERATURE=" << max(XHOST.vTemperatureCore) << ">" << AFLOW_CORE_TEMPERATURE_BEEP << "@ host=" << XHOST.hostname<< "]";
  // strout << endl;
  pthread_mutex_unlock(&mutex_INIT_Message);
  // do some killing
  //if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0) AFLOW_checkMEMORY("vasp",XHOST.maxmem);  //CO20170628 - this is already run above, very slow
  // if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0) AFLOW_checkMEMORY("aflow",XHOST.maxmem);
  // if(XHOST.maxmem>0.0 && XHOST.maxmem<100.0) AFLOW_checkMEMORY("clamd",XHOST.maxmem);
  if(!filename.empty()) strout << " - ["  <<  filename << "]";  //CO20200713
  return strout.str();
}
string Message(const string& filename,const _aflags& aflags,const string& list2print){ //CO20200713
  stringstream strout;  //CO20200713
  if(!list2print.empty()) strout << " - [dir=" << aflags.Directory << "]";  //CO20200713
  if(AFLOW_PTHREADS::FLAG) strout << " - [thread=" << aurostd::utype2string(aflags.AFLOW_PTHREADS_NUMBER) << "/" << aurostd::utype2string(AFLOW_PTHREADS::MAX_PTHREADS) << "]"; //CO20200713
  strout << Message(filename,list2print); //CO20200713
  return strout.str();  //CO20200713
}

// ***************************************************************************
// AFLOW_BlackList
// ***************************************************************************
bool AFLOW_BlackList(const string& h) { //CO20200713
  // cerr << h << endl;
  // if(h=="nietzsche" || h=="nietzsche.mems.duke.edu" || h=="material.duke.edu") return TRUE;
  if(h=="blacklisted_hostname") return TRUE;
  //  if(h=="m6-11-6") return TRUE;
  //  if(h=="m6-2-2") return TRUE;
  return FALSE;
}

// ***************************************************************************
// init::ErrorOptions
// ***************************************************************************
namespace init {
  void MessageOption(const string& options, const string& routine,vector<string> vusage) {  //CO20200624 //DX20200724 - bool to void
    ostream& oss=cerr;
    string usage="      "+routine+" Usage: ";  //CO20200624
    for(uint i=0;i<vusage.size();i++) {
      if(aurostd::substring2bool(vusage.at(i),"options:")) usage="              ";
      if(vusage.at(i)!="") oss << usage << vusage.at(i) << endl;
    }
    oss << "       options=[" << options << "]" << endl;
    //DX20200724 [OBSOLETE] return true;
  }
  void MessageOption(const string& options, const string& routine,string usage) { //CO20200624 //DX20200724 - bool to void
    vector<string> vusage;
    aurostd::string2vectorstring(usage,vusage);
    MessageOption(options,routine,vusage); //DX20200724 - removed return
  }
  void ErrorOption(const string& options, const string& routine,vector<string> vusage) { //DX20200724 - bool to void
    stringstream message;

    vector<string> tokens_options;
    aurostd::string2tokens(options,tokens_options,",");

    message << "Routine " << routine << ": Wrong number/type of input parameters (" << tokens_options.size();pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_ERROR_);
    //[CO20200624 - OBSOLETE]oss << "ERROR: " << routine << ":" << endl;
    //[CO20200624 - OBSOLETE]oss << "       Wrong number/type of input parameters! (" << tokens_options.size() << ")" << endl;

    //DX20200724 [OBSOLETE] return MessageOption(options,routine,vusage);
    MessageOption(options,routine,vusage); //DX20200724
    throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Incorrect usage.",_RUNTIME_ERROR_); //DX20200724
  }
  void ErrorOption(const string& options, const string& routine,string usage) { //DX20200725 - bool to void
    vector<string> vusage;
    aurostd::string2vectorstring(usage,vusage);
    ErrorOption(options,routine,vusage); //DX20200724 - removed return
  }
}

// ***************************************************************************
// init::SchemaFixName
// ***************************************************************************
namespace init {
  void SchemaFixName(string s1, string s2, string s3) {
    string line1="// schema is CAPITAL, content is not necessarily";
    string line2="XHOST.vschema.push_attached(\"SCHEMA::NAME:"+aurostd::toupper(s1)+"\",\""+s1+"\");";
    string line3="XHOST.vschema.push_attached(\"SCHEMA::UNIT:"+aurostd::toupper(s1)+"\",\""+s2+"\");";
    string line4="XHOST.vschema.push_attached(\"SCHEMA::TYPE:"+aurostd::toupper(s1)+"\",\""+s3+"\");";
    string line5="nschema++;";
    cout << line1 << endl;
    cout << line2 << endl;
    cout << line3 << endl;
    cout << line4 << endl;
    cout << line5 << endl;
    cout << endl;
    //  cout << aurostd::PaddedPOST(line2,105) << " // schema is CAPITAL" << endl;
    //  cout << aurostd::PaddedPOST(line3,105) << " // schema is CAPITAL" << endl;
    //  cout << aurostd::PaddedPOST(line4,105) << " // schema is CAPITAL" << endl;
    //  cout << "nschema++" << endl << endl;
  }
}

// ***************************************************************************
// init::InitSchema
// ***************************************************************************
namespace init {
  uint InitSchema(bool INIT_VERBOSE) {
    // DECLARATIONS
    bool LDEBUG=(FALSE || XHOST.DEBUG || INIT_VERBOSE);
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitSchema: [BEGIN]" << endl;

    uint nschema=0;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_APPLIED_PRESSURE","ael_applied_pressure");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_APPLIED_PRESSURE","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_APPLIED_PRESSURE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_AVERAGE_EXTERNAL_PRESSURE","ael_average_external_pressure");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_AVERAGE_EXTERNAL_PRESSURE","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_AVERAGE_EXTERNAL_PRESSURE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_BULK_MODULUS_REUSS","ael_bulk_modulus_reuss");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_BULK_MODULUS_REUSS","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_BULK_MODULUS_REUSS","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_BULK_MODULUS_VOIGT","ael_bulk_modulus_voigt");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_BULK_MODULUS_VOIGT","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_BULK_MODULUS_VOIGT","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_BULK_MODULUS_VRH","ael_bulk_modulus_vrh");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_BULK_MODULUS_VRH","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_BULK_MODULUS_VRH","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_COMPLIANCE_TENSOR","ael_compliance_tensor");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_COMPLIANCE_TENSOR","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_COMPLIANCE_TENSOR","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_DEBYE_TEMPERATURE","ael_debye_temperature");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_DEBYE_TEMPERATURE","K");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_DEBYE_TEMPERATURE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_ELASTIC_ANISOTROPY","ael_elastic_anisotropy");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_ELASTIC_ANISOTROPY","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_ELASTIC_ANISOTROPY","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_POISSON_RATIO","ael_poisson_ratio");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_POISSON_RATIO","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_POISSON_RATIO","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_PUGHS_MODULUS_RATIO","ael_pughs_modulus_ratio");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_PUGHS_MODULUS_RATIO","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_PUGHS_MODULUS_RATIO","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_SHEAR_MODULUS_REUSS","ael_shear_modulus_reuss");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_SHEAR_MODULUS_REUSS","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_SHEAR_MODULUS_REUSS","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_SHEAR_MODULUS_VOIGT","ael_shear_modulus_voigt");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_SHEAR_MODULUS_VOIGT","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_SHEAR_MODULUS_VOIGT","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_SHEAR_MODULUS_VRH","ael_shear_modulus_vrh");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_SHEAR_MODULUS_VRH","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_SHEAR_MODULUS_VRH","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_SPEED_SOUND_AVERAGE","ael_speed_sound_average");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_SPEED_SOUND_AVERAGE","m/s");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_SPEED_SOUND_AVERAGE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_SPEED_SOUND_LONGITUDINAL","ael_speed_sound_longitudinal");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_SPEED_SOUND_LONGITUDINAL","m/s");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_SPEED_SOUND_LONGITUDINAL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_SPEED_SOUND_TRANSVERSE","ael_speed_sound_transverse");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_SPEED_SOUND_TRANSVERSE","m/s");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_SPEED_SOUND_TRANSVERSE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_STIFFNESS_TENSOR","ael_stiffness_tensor");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_STIFFNESS_TENSOR","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_STIFFNESS_TENSOR","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AEL_YOUNGS_MODULUS_VRH","ael_youngs_modulus_vrh");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AEL_YOUNGS_MODULUS_VRH","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AEL_YOUNGS_MODULUS_VRH","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_VERSION","aflow_version");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_VERSION","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_VERSION","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOWLIB_DATE","aflowlib_date");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOWLIB_DATE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOWLIB_DATE","strings");  //CO+ME20200624
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOWLIB_VERSION","aflowlib_version");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOWLIB_VERSION","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOWLIB_VERSION","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_ACOUSTIC_DEBYE","agl_acoustic_debye");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_ACOUSTIC_DEBYE","K");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_ACOUSTIC_DEBYE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_BULK_MODULUS_ISOTHERMAL_300K","agl_bulk_modulus_isothermal_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_BULK_MODULUS_ISOTHERMAL_300K","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_BULK_MODULUS_ISOTHERMAL_300K","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_BULK_MODULUS_STATIC_300K","agl_bulk_modulus_static_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_BULK_MODULUS_STATIC_300K","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_BULK_MODULUS_STATIC_300K","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_DEBYE","agl_debye");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_DEBYE","K");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_DEBYE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_GRUNEISEN","agl_gruneisen");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_GRUNEISEN","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_GRUNEISEN","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_HEAT_CAPACITY_CP_300K","agl_heat_capacity_Cp_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_HEAT_CAPACITY_CP_300K","kB/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_HEAT_CAPACITY_CP_300K","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_HEAT_CAPACITY_CV_300K","agl_heat_capacity_Cv_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_HEAT_CAPACITY_CV_300K","kB/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_HEAT_CAPACITY_CV_300K","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_POISSON_RATIO_SOURCE","agl_poisson_ratio_source");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_POISSON_RATIO_SOURCE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_POISSON_RATIO_SOURCE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_THERMAL_CONDUCTIVITY_300K","agl_thermal_conductivity_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_THERMAL_CONDUCTIVITY_300K","W/(m K)");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_THERMAL_CONDUCTIVITY_300K","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_THERMAL_EXPANSION_300K","agl_thermal_expansion_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_THERMAL_EXPANSION_300K","K^-1");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_THERMAL_EXPANSION_300K","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_VIBRATIONAL_ENTROPY_300K_ATOM","agl_vibrational_entropy_300K_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_VIBRATIONAL_ENTROPY_300K_ATOM","meV/(K atom)");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_VIBRATIONAL_ENTROPY_300K_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_VIBRATIONAL_ENTROPY_300K_CELL","agl_vibrational_entropy_300K_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_VIBRATIONAL_ENTROPY_300K_CELL","meV/(K cell)");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_VIBRATIONAL_ENTROPY_300K_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_VIBRATIONAL_FREE_ENERGY_300K_ATOM","agl_vibrational_free_energy_300K_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_VIBRATIONAL_FREE_ENERGY_300K_ATOM","meV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_VIBRATIONAL_FREE_ENERGY_300K_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AGL_VIBRATIONAL_FREE_ENERGY_300K_CELL","agl_vibrational_free_energy_300K_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AGL_VIBRATIONAL_FREE_ENERGY_300K_CELL","meV/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AGL_VIBRATIONAL_FREE_ENERGY_300K_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily //DX20200902
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_PROTOTYPE_LABEL_ORIG","aflow_prototype_label_orig"); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_PROTOTYPE_LABEL_ORIG",""); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_PROTOTYPE_LABEL_ORIG","string"); //DX20201001 - updated keyword
    nschema++;

    // schema is CAPITAL, content is not necessarily //DX20200902
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_PROTOTYPE_PARAMS_LIST_ORIG","aflow_prototype_params_list_orig"); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_PROTOTYPE_PARAMS_LIST_ORIG",""); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_PROTOTYPE_PARAMS_LIST_ORIG","strings"); //DX20201001 - updated keyword
    nschema++;

    // schema is CAPITAL, content is not necessarily //DX20200902
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_PROTOTYPE_PARAMS_VALUES_ORIG","aflow_prototype_params_values_orig"); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_PROTOTYPE_PARAMS_VALUES_ORIG",""); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_PROTOTYPE_PARAMS_VALUES_ORIG","numbers"); //DX20201001 - updated keyword
    nschema++;

    // schema is CAPITAL, content is not necessarily //DX20200902
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_PROTOTYPE_LABEL_RELAX","aflow_prototype_label_relax"); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_PROTOTYPE_LABEL_RELAX",""); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_PROTOTYPE_LABEL_RELAX","string"); //DX20201001 - updated keyword
    nschema++;

    // schema is CAPITAL, content is not necessarily //DX20200902
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_PROTOTYPE_PARAMS_LIST_RELAX","aflow_prototype_params_list_relax"); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_PROTOTYPE_PARAMS_LIST_RELAX",""); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_PROTOTYPE_PARAMS_LIST_RELAX","strings"); //DX20201001 - updated keyword
    nschema++;

    // schema is CAPITAL, content is not necessarily //DX20200902
    XHOST.vschema.push_attached("SCHEMA::NAME:AFLOW_PROTOTYPE_PARAMS_VALUES_RELAX","aflow_prototype_params_values_relax"); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::UNIT:AFLOW_PROTOTYPE_PARAMS_VALUES_RELAX",""); //DX20201001 - updated keyword
    XHOST.vschema.push_attached("SCHEMA::TYPE:AFLOW_PROTOTYPE_PARAMS_VALUES_RELAX","numbers"); //DX20201001 - updated keyword
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AUID","auid");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AUID","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AUID","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:AURL","aurl");
    XHOST.vschema.push_attached("SCHEMA::UNIT:AURL","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:AURL","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BADER_ATOMIC_VOLUMES","bader_atomic_volumes");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BADER_ATOMIC_VOLUMES","A^3");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BADER_ATOMIC_VOLUMES","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BADER_NET_CHARGES","bader_net_charges");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BADER_NET_CHARGES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BADER_NET_CHARGES","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_LATTICE_SYSTEM","Bravais_lattice_lattice_system");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_LATTICE_SYSTEM","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_LATTICE_SYSTEM","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_LATTICE_SYSTEM_ORIG","Bravais_lattice_lattice_system_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_LATTICE_SYSTEM_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_LATTICE_SYSTEM_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_LATTICE_TYPE","Bravais_lattice_lattice_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_LATTICE_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_LATTICE_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_LATTICE_TYPE_ORIG","Bravais_lattice_lattice_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_LATTICE_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_LATTICE_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE","Bravais_lattice_lattice_variation_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE_ORIG","Bravais_lattice_lattice_variation_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_ORIG","Bravais_lattice_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_LATTICE_RELAX","Bravais_lattice_relax");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_LATTICE_RELAX","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_LATTICE_RELAX","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_SUPERLATTICE_LATTICE_SYSTEM","Bravais_superlattice_lattice_system");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_SUPERLATTICE_LATTICE_SYSTEM","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_SUPERLATTICE_LATTICE_SYSTEM","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_SUPERLATTICE_LATTICE_SYSTEM_ORIG","Bravais_superlattice_lattice_system_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_SUPERLATTICE_LATTICE_SYSTEM_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_SUPERLATTICE_LATTICE_SYSTEM_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_SUPERLATTICE_LATTICE_TYPE","Bravais_superlattice_lattice_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_SUPERLATTICE_LATTICE_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_SUPERLATTICE_LATTICE_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_SUPERLATTICE_LATTICE_TYPE_ORIG","Bravais_superlattice_lattice_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_SUPERLATTICE_LATTICE_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_SUPERLATTICE_LATTICE_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_SUPERLATTICE_LATTICE_VARIATION_TYPE","Bravais_superlattice_lattice_variation_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_SUPERLATTICE_LATTICE_VARIATION_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_SUPERLATTICE_LATTICE_VARIATION_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:BRAVAIS_SUPERLATTICE_LATTICE_VARIATION_TYPE_ORIG","Bravais_superlattice_lattice_variation_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:BRAVAIS_SUPERLATTICE_LATTICE_VARIATION_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:BRAVAIS_SUPERLATTICE_LATTICE_VARIATION_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CALCULATION_CORES","calculation_cores");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CALCULATION_CORES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CALCULATION_CORES","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CALCULATION_MEMORY","calculation_memory");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CALCULATION_MEMORY","MB");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CALCULATION_MEMORY","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CALCULATION_TIME","calculation_time");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CALCULATION_TIME","seconds");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CALCULATION_TIME","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CATALOG","catalog");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CATALOG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CATALOG","string");
    nschema++;

    //CO20200829 START
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:COEFFICIENT_ENTROPY_STABILIZATION","coefficient_entropy_stabilization");
    XHOST.vschema.push_attached("SCHEMA::UNIT:COEFFICIENT_ENTROPY_STABILIZATION","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:COEFFICIENT_ENTROPY_STABILIZATION","number");
    nschema++;
    //CO20200829 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CODE","code");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CODE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CODE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:COMPOSITION","composition");
    XHOST.vschema.push_attached("SCHEMA::UNIT:COMPOSITION","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:COMPOSITION","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:COMPOUND","compound");
    XHOST.vschema.push_attached("SCHEMA::UNIT:COMPOUND","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:COMPOUND","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CRYSTAL_CLASS","crystal_class");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CRYSTAL_CLASS","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CRYSTAL_CLASS","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CRYSTAL_CLASS_ORIG","crystal_class_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CRYSTAL_CLASS_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CRYSTAL_CLASS_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CRYSTAL_FAMILY","crystal_family");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CRYSTAL_FAMILY","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CRYSTAL_FAMILY","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CRYSTAL_FAMILY_ORIG","crystal_family_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CRYSTAL_FAMILY_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CRYSTAL_FAMILY_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CRYSTAL_SYSTEM","crystal_system");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CRYSTAL_SYSTEM","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CRYSTAL_SYSTEM","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:CRYSTAL_SYSTEM_ORIG","crystal_system_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:CRYSTAL_SYSTEM_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:CRYSTAL_SYSTEM_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DATA_API","data_api");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DATA_API","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DATA_API","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DATA_SOURCE","data_source");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DATA_SOURCE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DATA_SOURCE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DELTA_ELECTRONIC_ENERGY_CONVERGENCE","delta_electronic_energy_convergence");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DELTA_ELECTRONIC_ENERGY_CONVERGENCE","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DELTA_ELECTRONIC_ENERGY_CONVERGENCE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DELTA_ELECTRONIC_ENERGY_THRESHOLD","delta_electronic_energy_threshold");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DELTA_ELECTRONIC_ENERGY_THRESHOLD","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DELTA_ELECTRONIC_ENERGY_THRESHOLD","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DENSITY","density");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DENSITY","g/cm^3");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DENSITY","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DFT_TYPE","dft_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DFT_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DFT_TYPE","string");
    nschema++;

    //CO20200829 START
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:DISTANCE_HULL","distance_hull");
    XHOST.vschema.push_attached("SCHEMA::UNIT:DISTANCE_HULL","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:DISTANCE_HULL","number");
    nschema++;
    //CO20200829 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:EENTROPY_ATOM","eentropy_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:EENTROPY_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:EENTROPY_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:EENTROPY_CELL","eentropy_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:EENTROPY_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:EENTROPY_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:EGAP","Egap");
    XHOST.vschema.push_attached("SCHEMA::UNIT:EGAP","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:EGAP","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:EGAP_FIT","Egap_fit");
    XHOST.vschema.push_attached("SCHEMA::UNIT:EGAP_FIT","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:EGAP_FIT","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:EGAP_TYPE","Egap_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:EGAP_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:EGAP_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_ATOM","energy_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_CELL","energy_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_CUTOFF","energy_cutoff");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_CUTOFF","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_CUTOFF","number");
    nschema++;

    //AS20201008 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_FREE_ATOM_QHA_300K","energy_free_atom_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_FREE_ATOM_QHA_300K","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_FREE_ATOM_QHA_300K","number");
    //AS20201008 END
    //
    //AS20201207 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_FREE_CELL_QHA_300K","energy_free_cell_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_FREE_CELL_QHA_300K","eV/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_FREE_CELL_QHA_300K","number");
    //AS20201207 END

    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_FREE_VIBRATIONAL_ATOM_APL_300K","energy_free_vibrational_atom_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_FREE_VIBRATIONAL_ATOM_APL_300K","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_FREE_VIBRATIONAL_ATOM_APL_300K","number");
    //ME20210927 END
    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_FREE_VIBRATIONAL_CELL_APL_300K","energy_free_vibrational_cell_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_FREE_VIBRATIONAL_CELL_APL_300K","eV/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_FREE_VIBRATIONAL_CELL_APL_300K","number");
    //ME20210927 END

    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_INTERNAL_VIBRATIONAL_ATOM_APL_300K","energy_internal_vibrational_atom_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_INTERNAL_VIBRATIONAL_ATOM_APL_300K","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_INTERNAL_VIBRATIONAL_ATOM_APL_300K","number");
    //ME20210927 END
    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_INTERNAL_VIBRATIONAL_CELL_APL_300K","energy_internal_vibrational_cell_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_INTERNAL_VIBRATIONAL_CELL_APL_300K","eV/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_INTERNAL_VIBRATIONAL_CELL_APL_300K","number");
    //ME20210927 END

    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_ZERO_POINT_ATOM_APL","energy_zero_point_atom_apl");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_ZERO_POINT_ATOM_APL","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_ZERO_POINT_ATOM_APL","number");
    //ME20210927 END
    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENERGY_ZERO_POINT_CELL_APL","energy_zero_point_cell_apl");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENERGY_ZERO_POINT_CELL_APL","eV/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENERGY_ZERO_POINT_CELL_APL","number");
    //ME20210927 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_ATOM","enthalpy_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_CELL","enthalpy_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_FORMATION_ATOM","enthalpy_formation_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_FORMATION_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_FORMATION_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_FORMATION_CELL","enthalpy_formation_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_FORMATION_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_FORMATION_CELL","number");
    nschema++;

    //CO20200829 START
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_FORMATION_CCE_300K_CELL","enthalpy_formation_cce_300K_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_FORMATION_CCE_300K_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_FORMATION_CCE_300K_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_FORMATION_CCE_300K_ATOM","enthalpy_formation_cce_300K_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_FORMATION_CCE_300K_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_FORMATION_CCE_300K_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_FORMATION_CCE_0K_CELL","enthalpy_formation_cce_0K_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_FORMATION_CCE_0K_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_FORMATION_CCE_0K_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTHALPY_FORMATION_CCE_0K_ATOM","enthalpy_formation_cce_0K_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTHALPY_FORMATION_CCE_0K_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTHALPY_FORMATION_CCE_0K_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTROPY_FORMING_ABILITY","entropy_forming_ability");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTROPY_FORMING_ABILITY","(eV/atom)^{-1}");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTROPY_FORMING_ABILITY","number");
    nschema++;
    //CO20200829 END

    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTROPY_VIBRATIONAL_ATOM_APL_300K","entropy_vibrational_atom_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTROPY_VIBRATIONAL_ATOM_APL_300K","meV/(K atom)");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTROPY_VIBRATIONAL_ATOM_APL_300K","number");
    //ME20210927 END
    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTROPY_VIBRATIONAL_CELL_APL_300K","entropy_vibrational_cell_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTROPY_VIBRATIONAL_CELL_APL_300K","meV/(K cell)");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTROPY_VIBRATIONAL_CELL_APL_300K","number");
    //ME20210927 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:ENTROPIC_TEMPERATURE","entropic_temperature");
    XHOST.vschema.push_attached("SCHEMA::UNIT:ENTROPIC_TEMPERATURE","K");
    XHOST.vschema.push_attached("SCHEMA::TYPE:ENTROPIC_TEMPERATURE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:FILES","files");
    XHOST.vschema.push_attached("SCHEMA::UNIT:FILES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:FILES","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:FORCES","forces");
    XHOST.vschema.push_attached("SCHEMA::UNIT:FORCES","eV/A");
    XHOST.vschema.push_attached("SCHEMA::TYPE:FORCES","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:GEOMETRY","geometry");
    XHOST.vschema.push_attached("SCHEMA::UNIT:GEOMETRY","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:GEOMETRY","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:GEOMETRY_ORIG","geometry_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:GEOMETRY_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:GEOMETRY_ORIG","numbers");
    nschema++;

    //CO20200829 START
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:GROUND_STATE","ground_state");
    XHOST.vschema.push_attached("SCHEMA::UNIT:GROUND_STATE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:GROUND_STATE","bool");
    nschema++;
    //CO20200829 END

    //AS20200915
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:GRUNEISEN_QHA","gruneisen_qha");
    XHOST.vschema.push_attached("SCHEMA::UNIT:GRUNEISEN_QHA","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:GRUNEISEN_QHA","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:GRUNEISEN_QHA_300K","gruneisen_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:GRUNEISEN_QHA_300K","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:GRUNEISEN_QHA_300K","number");
    nschema++;
    //AS20200915 END

    //AS20201008 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:HEAT_CAPACITY_CP_ATOM_QHA_300K","heat_capacity_Cv_atom_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:HEAT_CAPACITY_CP_ATOM_QHA_300K","kB/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:HEAT_CAPACITY_CP_ATOM_QHA_300K","number");

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:HEAT_CAPACITY_CV_ATOM_QHA_300K","heat_capacity_Cp_atom_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:HEAT_CAPACITY_CV_ATOM_QHA_300K","kB/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:HEAT_CAPACITY_CV_ATOM_QHA_300K","number");
    //AS20201008 END

    //AS20201207 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:HEAT_CAPACITY_CP_CELL_QHA_300K","heat_capacity_Cv_cell_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:HEAT_CAPACITY_CP_CELL_QHA_300K","kB/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:HEAT_CAPACITY_CP_CELL_QHA_300K","number");

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:HEAT_CAPACITY_CV_CELL_QHA_300K","heat_capacity_Cp_cell_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:HEAT_CAPACITY_CV_CELL_QHA_300K","kB/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:HEAT_CAPACITY_CV_CELL_QHA_300K","number");
    //AS20201207 END

    //ME20210927 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:HEAT_CAPACITY_CV_ATOM_APL_300K","heat_capacity_Cv_atom_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:HEAT_CAPACITY_CV_ATOM_APL_300K","kB/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:HEAT_CAPACITY_CV_ATOM_APL_300K","number");
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:HEAT_CAPACITY_CV_CELL_APL_300K","heat_capacity_Cv_cell_apl_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:HEAT_CAPACITY_CV_CELL_APL_300K","kB/cell");
    XHOST.vschema.push_attached("SCHEMA::TYPE:HEAT_CAPACITY_CV_CELL_APL_300K","number");
    //ME20210927 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:KPOINTS","kpoints");
    XHOST.vschema.push_attached("SCHEMA::UNIT:KPOINTS","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:KPOINTS","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:KPOINTS_BANDS_NKPTS","kpoints_bands_nkpts");
    XHOST.vschema.push_attached("SCHEMA::UNIT:KPOINTS_BANDS_NKPTS","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:KPOINTS_BANDS_NKPTS","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:KPOINTS_BANDS_PATH","kpoints_bands_path");
    XHOST.vschema.push_attached("SCHEMA::UNIT:KPOINTS_BANDS_PATH","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:KPOINTS_BANDS_PATH","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:KPOINTS_RELAX","kpoints_relax");
    XHOST.vschema.push_attached("SCHEMA::UNIT:KPOINTS_RELAX","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:KPOINTS_RELAX","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:KPOINTS_STATIC","kpoints_static");
    XHOST.vschema.push_attached("SCHEMA::UNIT:KPOINTS_STATIC","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:KPOINTS_STATIC","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LATTICE_SYSTEM_ORIG","lattice_system_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LATTICE_SYSTEM_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LATTICE_SYSTEM_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LATTICE_SYSTEM_RELAX","lattice_system_relax");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LATTICE_SYSTEM_RELAX","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LATTICE_SYSTEM_RELAX","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LATTICE_VARIATION_ORIG","lattice_variation_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LATTICE_VARIATION_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LATTICE_VARIATION_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LATTICE_VARIATION_RELAX","lattice_variation_relax");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LATTICE_VARIATION_RELAX","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LATTICE_VARIATION_RELAX","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LDAU_J","ldau_j");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LDAU_J","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LDAU_J","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LDAU_L","ldau_l");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LDAU_L","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LDAU_L","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LDAU_TLUJ","ldau_TLUJ");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LDAU_TLUJ","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LDAU_TLUJ","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LDAU_TYPE","ldau_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LDAU_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LDAU_TYPE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LDAU_U","ldau_u");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LDAU_U","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LDAU_U","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:LOOP","loop");
    XHOST.vschema.push_attached("SCHEMA::UNIT:LOOP","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:LOOP","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:METAGGA","metagga");
    XHOST.vschema.push_attached("SCHEMA::UNIT:METAGGA","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:METAGGA","string");
    nschema++;

    //AS20200915 BEGIN
    XHOST.vschema.push_attached("SCHEMA::NAME:MODULUS_BULK_QHA_300K","modulus_bulk_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:MODULUS_BULK_QHA_300K","GPa");
    XHOST.vschema.push_attached("SCHEMA::TYPE:MODULUS_BULK_QHA_300K","number");
    nschema++;
    //AS20200915 END

    //AS20201008 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:MODULUS_BULK_DERIVATIVE_PRESSURE_QHA_300K","modulus_bulk_derivative_pressure_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:MODULUS_BULK_DERIVATIVE_PRESSURE_QHA_300K","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:MODULUS_BULK_DERIVATIVE_PRESSURE_QHA_300K","number");
    nschema++;
    //AS20201008 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NATOMS","natoms");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NATOMS","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NATOMS","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NBONDXX","nbondxx");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NBONDXX","A");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NBONDXX","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NODE_CPU_CORES","node_CPU_Cores");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NODE_CPU_CORES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NODE_CPU_CORES","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NODE_CPU_MHZ","node_CPU_MHz");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NODE_CPU_MHZ","MHz");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NODE_CPU_MHZ","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NODE_CPU_MODEL","node_CPU_Model");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NODE_CPU_MODEL","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NODE_CPU_MODEL","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NODE_RAM_GB","node_RAM_GB");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NODE_RAM_GB","GB");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NODE_RAM_GB","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:NSPECIES","nspecies");
    XHOST.vschema.push_attached("SCHEMA::UNIT:NSPECIES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:NSPECIES","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PEARSON_SYMBOL_ORIG","Pearson_symbol_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PEARSON_SYMBOL_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PEARSON_SYMBOL_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PEARSON_SYMBOL_RELAX","Pearson_symbol_relax");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PEARSON_SYMBOL_RELAX","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PEARSON_SYMBOL_RELAX","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PEARSON_SYMBOL_SUPERLATTICE","Pearson_symbol_superlattice");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PEARSON_SYMBOL_SUPERLATTICE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PEARSON_SYMBOL_SUPERLATTICE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PEARSON_SYMBOL_SUPERLATTICE_ORIG","Pearson_symbol_superlattice_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PEARSON_SYMBOL_SUPERLATTICE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PEARSON_SYMBOL_SUPERLATTICE_ORIG","string");
    nschema++;

    //CO20200829 START
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POCC_PARAMETERS","pocc_parameters");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POCC_PARAMETERS","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POCC_PARAMETERS","string");
    nschema++;
    //CO20200829 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_HERMANN_MAUGUIN","point_group_Hermann_Mauguin");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_HERMANN_MAUGUIN","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_HERMANN_MAUGUIN","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_HERMANN_MAUGUIN_ORIG","point_group_Hermann_Mauguin_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_HERMANN_MAUGUIN_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_HERMANN_MAUGUIN_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_ORBIFOLD","point_group_orbifold");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_ORBIFOLD","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_ORBIFOLD","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_ORBIFOLD_ORIG","point_group_orbifold_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_ORBIFOLD_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_ORBIFOLD_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_ORDER","point_group_order");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_ORDER","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_ORDER","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_ORDER_ORIG","point_group_order_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_ORDER_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_ORDER_ORIG","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_SCHOENFLIES","point_group_Schoenflies");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_SCHOENFLIES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_SCHOENFLIES","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_SCHOENFLIES_ORIG","point_group_Schoenflies_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_SCHOENFLIES_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_SCHOENFLIES_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_STRUCTURE","point_group_structure");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_STRUCTURE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_STRUCTURE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_STRUCTURE_ORIG","point_group_structure_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_STRUCTURE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_STRUCTURE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_TYPE","point_group_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POINT_GROUP_TYPE_ORIG","point_group_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POINT_GROUP_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POINT_GROUP_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POSITIONS_CARTESIAN","positions_cartesian");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POSITIONS_CARTESIAN","A");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POSITIONS_CARTESIAN","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:POSITIONS_FRACTIONAL","positions_fractional");
    XHOST.vschema.push_attached("SCHEMA::UNIT:POSITIONS_FRACTIONAL","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:POSITIONS_FRACTIONAL","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PRESSURE","pressure");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PRESSURE","kbar");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PRESSURE","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PRESSURE_FINAL","pressure_final");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PRESSURE_FINAL","kbar");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PRESSURE_FINAL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PRESSURE_RESIDUAL","pressure_residual");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PRESSURE_RESIDUAL","kbar");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PRESSURE_RESIDUAL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PULAY_STRESS","Pulay_stress");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PULAY_STRESS","kbar");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PULAY_STRESS","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PV_ATOM","PV_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PV_ATOM","eV/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PV_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PV_CELL","PV_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PV_CELL","eV");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PV_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:PROTOTYPE","prototype");
    XHOST.vschema.push_attached("SCHEMA::UNIT:PROTOTYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:PROTOTYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_GEOMETRY_RELAX","reciprocal_geometry_relax");  //CO20220719 _relax
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_GEOMETRY_RELAX",""); //CO20220719 _relax
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_GEOMETRY_RELAX","numbers");  //CO20220719 _relax
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_GEOMETRY_ORIG","reciprocal_geometry_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_GEOMETRY_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_GEOMETRY_ORIG","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_LATTICE_TYPE","reciprocal_lattice_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_LATTICE_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_LATTICE_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_LATTICE_TYPE_ORIG","reciprocal_lattice_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_LATTICE_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_LATTICE_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_LATTICE_VARIATION_TYPE","reciprocal_lattice_variation_type");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_LATTICE_VARIATION_TYPE","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_LATTICE_VARIATION_TYPE","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_LATTICE_VARIATION_TYPE_ORIG","reciprocal_lattice_variation_type_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_LATTICE_VARIATION_TYPE_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_LATTICE_VARIATION_TYPE_ORIG","string");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_VOLUME_CELL","reciprocal_volume_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_VOLUME_CELL","A^-3");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_VOLUME_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:RECIPROCAL_VOLUME_CELL_ORIG","reciprocal_volume_cell_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:RECIPROCAL_VOLUME_CELL_ORIG","A^-3/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:RECIPROCAL_VOLUME_CELL_ORIG","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SCINTILLATION_ATTENUATION_LENGTH","scintillation_attenuation_length");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SCINTILLATION_ATTENUATION_LENGTH","cm");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SCINTILLATION_ATTENUATION_LENGTH","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SG","sg");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SG","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SG2","sg2");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SG2","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SG2","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPACEGROUP_ORIG","spacegroup_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPACEGROUP_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPACEGROUP_ORIG","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPACEGROUP_RELAX","spacegroup_relax");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPACEGROUP_RELAX","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPACEGROUP_RELAX","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPECIES","species");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPECIES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPECIES","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPECIES_PP","species_pp");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPECIES_PP","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPECIES_PP","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPECIES_PP_AUID","species_pp_AUID");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPECIES_PP_AUID","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPECIES_PP_AUID","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPECIES_PP_ZVAL","species_pp_ZVAL");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPECIES_PP_ZVAL","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPECIES_PP_ZVAL","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPECIES_PP_VERSION","species_pp_version");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPECIES_PP_VERSION","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPECIES_PP_VERSION","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPIND","spinD");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPIND","uB");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPIND","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPINF","spinF");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPINF","uB");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPINF","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPIN_ATOM","spin_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPIN_ATOM","uB/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPIN_ATOM","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:SPIN_CELL","spin_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:SPIN_CELL","uB");
    XHOST.vschema.push_attached("SCHEMA::TYPE:SPIN_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:STOICHIOMETRY","stoichiometry");
    XHOST.vschema.push_attached("SCHEMA::UNIT:STOICHIOMETRY","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:STOICHIOMETRY","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:STRESS_TENSOR","stress_tensor");
    XHOST.vschema.push_attached("SCHEMA::UNIT:STRESS_TENSOR","kbar");
    XHOST.vschema.push_attached("SCHEMA::TYPE:STRESS_TENSOR","numbers");
    nschema++;

    //AS20200915 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:THERMAL_EXPANSION_QHA_300K","thermal_expansion_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:THERMAL_EXPANSION_QHA_300K","K^-1");
    XHOST.vschema.push_attached("SCHEMA::TYPE:THERMAL_EXPANSION_QHA_300K","number");
    nschema++;
    //AS20200915 END

    // schema is CAPITAL, content is not necessarily
    // OBSOLETE ME20200829
    //XHOST.vschema.push_attached("SCHEMA::NAME:TITLE","title");
    //XHOST.vschema.push_attached("SCHEMA::UNIT:TITLE","");
    //XHOST.vschema.push_attached("SCHEMA::TYPE:TITLE","string");
    //nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:VALENCE_CELL_IUPAC","valence_cell_iupac");
    XHOST.vschema.push_attached("SCHEMA::UNIT:VALENCE_CELL_IUPAC","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:VALENCE_CELL_IUPAC","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:VALENCE_CELL_STD","valence_cell_std");
    XHOST.vschema.push_attached("SCHEMA::UNIT:VALENCE_CELL_STD","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:VALENCE_CELL_STD","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:VOLUME_ATOM","volume_atom");
    XHOST.vschema.push_attached("SCHEMA::UNIT:VOLUME_ATOM","A^3/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:VOLUME_ATOM","number");
    nschema++;

    //AS20201008 BEGIN
    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:VOLUME_ATOM_QHA_300K","volume_atom_qha_300K");
    XHOST.vschema.push_attached("SCHEMA::UNIT:VOLUME_ATOM_QHA_300K","A^3/atom");
    XHOST.vschema.push_attached("SCHEMA::TYPE:VOLUME_ATOM_QHA_300K","number");
    //AS20201008 END

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:VOLUME_CELL","volume_cell");
    XHOST.vschema.push_attached("SCHEMA::UNIT:VOLUME_CELL","A^3");
    XHOST.vschema.push_attached("SCHEMA::TYPE:VOLUME_CELL","number");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:WYCKOFF_LETTERS","Wyckoff_letters");
    XHOST.vschema.push_attached("SCHEMA::UNIT:WYCKOFF_LETTERS","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:WYCKOFF_LETTERS","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:WYCKOFF_LETTERS_ORIG","Wyckoff_letters_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:WYCKOFF_LETTERS_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:WYCKOFF_LETTERS_ORIG","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:WYCKOFF_MULTIPLICITIES","Wyckoff_multiplicities");
    XHOST.vschema.push_attached("SCHEMA::UNIT:WYCKOFF_MULTIPLICITIES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:WYCKOFF_MULTIPLICITIES","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:WYCKOFF_MULTIPLICITIES_ORIG","Wyckoff_multiplicities_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:WYCKOFF_MULTIPLICITIES_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:WYCKOFF_MULTIPLICITIES_ORIG","numbers");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:WYCKOFF_SITE_SYMMETRIES","Wyckoff_site_symmetries");
    XHOST.vschema.push_attached("SCHEMA::UNIT:WYCKOFF_SITE_SYMMETRIES","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:WYCKOFF_SITE_SYMMETRIES","strings");
    nschema++;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema.push_attached("SCHEMA::NAME:WYCKOFF_SITE_SYMMETRIES_ORIG","Wyckoff_site_symmetries_orig");
    XHOST.vschema.push_attached("SCHEMA::UNIT:WYCKOFF_SITE_SYMMETRIES_ORIG","");
    XHOST.vschema.push_attached("SCHEMA::TYPE:WYCKOFF_SITE_SYMMETRIES_ORIG","strings");
    nschema++;

    // read them as
    LDEBUG=0; //DX+CO put LDEBUG=1 so you seen how they answer
    if(LDEBUG) cerr << "XHOST.vschema.isdefined(\"SCHEMA::NAME:AEL_BULK_MODULUS_REUSS\")=" << XHOST.vschema.isscheme("SCHEMA::NAME:AEL_BULK_MODULUS_REUSS") << endl; // set or not set
    if(LDEBUG) cerr << "XHOST.vschema.getattachedscheme(\"SCHEMA::NAME:AEL_BULK_MODULUS_REUSS\")=" << XHOST.vschema.getattachedscheme("SCHEMA::NAME:AEL_BULK_MODULUS_REUSS") << endl;  // content
    if(LDEBUG) cerr << "XHOST.vschema.isdefined(\"SCHEMA::UNIT:AEL_BULK_MODULUS_REUSS\")=" << XHOST.vschema.isscheme("SCHEMA::UNIT:AEL_BULK_MODULUS_REUSS") << endl; // set or not set
    if(LDEBUG) cerr << "XHOST.vschema.getattachedscheme(\"SCHEMA::UNIT:AEL_BULK_MODULUS_REUSS\")=" << XHOST.vschema.getattachedscheme("SCHEMA::UNIT:AEL_BULK_MODULUS_REUSS") << endl;  // content
    if(LDEBUG) cerr << "XHOST.vschema.isdefined(\"SCHEMA::TYPE:AEL_BULK_MODULUS_REUSS\")=" << XHOST.vschema.isscheme("SCHEMA::TYPE:AEL_BULK_MODULUS_REUSS") << endl; // set or not set
    if(LDEBUG) cerr << "XHOST.vschema.getattachedscheme(\"SCHEMA::TYPE:AEL_BULK_MODULUS_REUSS\")=" << XHOST.vschema.getattachedscheme("SCHEMA::TYPE:AEL_BULK_MODULUS_REUSS") << endl;  // content
    if(LDEBUG) cerr << "XHOST.vschema.isdefined(\"SCHEMA::TYPE:AEL_BULK_MODULUS_OSES\")=" << XHOST.vschema.isscheme("SCHEMA::TYPE:AEL_BULK_MODULUS_OSES") << endl; // set or not set
    if(LDEBUG) cerr << "XHOST.vschema.getattachedscheme(\"SCHEMA::TYPE:AEL_BULK_MODULUS_OSES\")=" << XHOST.vschema.getattachedscheme("SCHEMA::TYPE:AEL_BULK_MODULUS_OSES") << endl;  // content

    if(LDEBUG) cerr << "nschema=" << nschema << endl;
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitSchema: [END]" << endl;
    if(LDEBUG) return 0;

    return nschema;

    // since DX gave me stuff in small letters, I used this to fix the text.. you can use it too if you have a bunch of them

    //SchemaFixName("ael_bulk_modulus_reuss","GPa","number");
    //SchemaFixName("ael_bulk_modulus_voigt","GPa","number");
    //SchemaFixName("ael_bulk_modulus_vrh","GPa","number");
    //SchemaFixName("ael_compliance_tensor","","numbers");
    //SchemaFixName("ael_elastic_anisotropy","","number");
    //SchemaFixName("ael_poisson_ratio","","number");
    //SchemaFixName("ael_shear_modulus_reuss","GPa","number");
    //SchemaFixName("ael_shear_modulus_voigt","GPa","number");
    //SchemaFixName("ael_shear_modulus_vrh","GPa","number");
    //SchemaFixName("ael_stiffness_tensor","","numbers");
    //SchemaFixName("aflow_version","","string");
    //SchemaFixName("aflowlib_date","","string");
    //SchemaFixName("aflowlib_version","","string");
    //SchemaFixName("agl_acoustic_debye","K","number");
    //SchemaFixName("agl_bulk_modulus_isothermal_300K","GPa","number");
    //SchemaFixName("agl_bulk_modulus_static_300K","GPa","number");
    //SchemaFixName("agl_debye","K","number");
    //SchemaFixName("agl_gruneisen","","number");
    //SchemaFixName("agl_heat_capacity_Cp_300K","kB/cell","number");
    //SchemaFixName("agl_heat_capacity_Cv_300K","kB/cell","number");
    //SchemaFixName("agl_thermal_conductivity_300K","W/(m K)","number");
    //SchemaFixName("agl_thermal_expansion_300K","K^-1","number");
    //SchemaFixName("auid","","string");
    //SchemaFixName("aurl","","string");
    //SchemaFixName("bader_atomic_volumes","A^3","numbers");
    //SchemaFixName("bader_net_charges","","numbers");
    //SchemaFixName("Bravais_lattice_lattice_system","","string");
    //SchemaFixName("Bravais_lattice_lattice_system_orig","","string");
    //SchemaFixName("Bravais_lattice_lattice_type","","string");
    //SchemaFixName("Bravais_lattice_lattice_type_orig","","string");
    //SchemaFixName("Bravais_lattice_lattice_variation_type","","string");
    //SchemaFixName("Bravais_lattice_lattice_variation_type_orig","","string");
    //SchemaFixName("Bravais_lattice_orig","","string");
    //SchemaFixName("Bravais_lattice_relax","","string");
    //SchemaFixName("Bravais_superlattice_lattice_system","","string");
    //SchemaFixName("Bravais_superlattice_lattice_system_orig","","string");
    //SchemaFixName("Bravais_superlattice_lattice_type","","string");
    //SchemaFixName("Bravais_superlattice_lattice_type_orig","","string");
    //SchemaFixName("Bravais_superlattice_lattice_variation_type","","string");
    //SchemaFixName("Bravais_superlattice_lattice_variation_type_orig","","string");
    //SchemaFixName("calculation_cores","","number");
    //SchemaFixName("calculation_memory","MB","number");
    //SchemaFixName("calculation_time","seconds","number");
    //SchemaFixName("catalog","","string");
    //SchemaFixName("code","","string");
    //SchemaFixName("composition","","numbers");
    //SchemaFixName("compound","","string");
    //SchemaFixName("crystal_class","","string");
    //SchemaFixName("crystal_class_orig","","string");
    //SchemaFixName("crystal_family","","string");
    //SchemaFixName("crystal_family_orig","","string");
    //SchemaFixName("crystal_system","","string");
    //SchemaFixName("crystal_system_orig","","string");
    //SchemaFixName("data_api","","string");
    //SchemaFixName("data_source","","string");
    //SchemaFixName("delta_electronic_energy_convergence","eV","number");
    //SchemaFixName("delta_electronic_energy_threshold","eV","number");
    //SchemaFixName("density","g/cm^3","number");
    //SchemaFixName("dft_type","","string");
    //SchemaFixName("eentropy_atom","eV/atom","number");
    //SchemaFixName("eentropy_cell","eV","number");
    //SchemaFixName("Egap","eV","number");
    //SchemaFixName("Egap_fit","eV","number");
    //SchemaFixName("Egap_type","","string");
    //SchemaFixName("energy_atom","eV/atom","number");
    //SchemaFixName("energy_cell","eV","number");
    //SchemaFixName("energy_cutoff","eV","number");
    //SchemaFixName("enthalpy_atom","eV/atom","number");
    //SchemaFixName("enthalpy_cell","eV","number");
    //SchemaFixName("enthalpy_formation_atom","eV/atom","number");
    //SchemaFixName("enthalpy_formation_cell","eV","number");
    //SchemaFixName("entropic_temperature","K","number");
    //SchemaFixName("files","","strings");
    //SchemaFixName("forces","eV/A","numbers");
    //SchemaFixName("geometry","","numbers");
    //SchemaFixName("geometry_orig","","numbers");
    //SchemaFixName("kpoints","","strings");
    //SchemaFixName("kpoints_bands_nkpts","","number");
    //SchemaFixName("kpoints_bands_path","","strings");
    //SchemaFixName("kpoints_relax","","numbers");
    //SchemaFixName("kpoints_static","","numbers");
    //SchemaFixName("lattice_system_orig","","string");
    //SchemaFixName("lattice_system_relax","","string");
    //SchemaFixName("lattice_variation_orig","","string");
    //SchemaFixName("lattice_variation_relax","","string");
    //SchemaFixName("ldau_j","eV","numbers");
    //SchemaFixName("ldau_l","","numbers");
    //SchemaFixName("ldau_TLUJ","","numbers");
    //SchemaFixName("ldau_type","","number");
    //SchemaFixName("ldau_u","eV","numbers");
    //SchemaFixName("loop","","strings");
    //SchemaFixName("natoms","","number");
    //SchemaFixName("nbondxx","A","numbers");
    //SchemaFixName("node_CPU_Cores","","number");
    //SchemaFixName("node_CPU_MHz","MHz","number");
    //SchemaFixName("node_CPU_Model","","string");
    //SchemaFixName("node_RAM_GB","GB","number");
    //SchemaFixName("nspecies","","number");
    //SchemaFixName("Pearson_symbol_orig","","string");
    //SchemaFixName("Pearson_symbol_relax","","string");
    //SchemaFixName("Pearson_symbol_superlattice","","string");
    //SchemaFixName("Pearson_symbol_superlattice_orig","","string");
    //SchemaFixName("point_group_Hermann_Mauguin","","string");
    //SchemaFixName("point_group_Hermann_Mauguin_orig","","string");
    //SchemaFixName("point_group_orbifold","","string");
    //SchemaFixName("point_group_orbifold_orig","","string");
    //SchemaFixName("point_group_order","","number");
    //SchemaFixName("point_group_order_orig","","number");
    //SchemaFixName("point_group_Schoenflies","","string");
    //SchemaFixName("point_group_Schoenflies_orig","","string");
    //SchemaFixName("point_group_structure","","string");
    //SchemaFixName("point_group_structure_orig","","string");
    //SchemaFixName("point_group_type","","string");
    //SchemaFixName("point_group_type_orig","","string");
    //SchemaFixName("positions_cartesian","A","numbers");
    //SchemaFixName("positions_fractional","","numbers");
    //SchemaFixName("pressure","kbar" ,"number");
    //SchemaFixName("pressure_final","kbar" ,"number");
    //SchemaFixName("pressure_residual","kbar" ,"number");
    //SchemaFixName("Pulay_stress","kbar","number");
    //SchemaFixName("PV_atom","eV/atom","number");
    //SchemaFixName("PV_cell","eV","number");
    //SchemaFixName("prototype","","string");
    //SchemaFixName("reciprocal_geometry","","numbers");
    //SchemaFixName("reciprocal_geometry_orig","","numbers");
    //SchemaFixName("reciprocal_lattice_type","","string");
    //SchemaFixName("reciprocal_lattice_type_orig","","string");
    //SchemaFixName("reciprocal_lattice_variation_type","","string");
    //SchemaFixName("reciprocal_lattice_variation_type_orig","","string");
    //SchemaFixName("reciprocal_volume_cell","A^-3","number");
    //SchemaFixName("reciprocal_volume_cell_orig","A^-3/atom","number");
    //SchemaFixName("scintillation_attenuation_length","cm","number");
    //SchemaFixName("sg","","strings");
    //SchemaFixName("sg2","","strings");
    //SchemaFixName("spacegroup_orig","","number");
    //SchemaFixName("spacegroup_relax","","number");
    //SchemaFixName("species","","strings");
    //SchemaFixName("species_pp","","strings");
    //SchemaFixName("species_pp_ZVAL","","numbers");
    //SchemaFixName("species_pp_version","","strings");
    //SchemaFixName("spinD","uB","numbers");
    //SchemaFixName("spinF","uB","number");
    //SchemaFixName("spin_atom","uB/atom","number");
    //SchemaFixName("spin_cell","uB","number");
    //SchemaFixName("stoichiometry","","numbers");
    //SchemaFixName("stress_tensor","kbar","numbers");
    //SchemaFixName("title","","string");
    //SchemaFixName("valence_cell_iupac","","number");
    //SchemaFixName("valence_cell_std","","number");
    //SchemaFixName("volume_atom","A^3/atom","number");
    //SchemaFixName("volume_cell","A^3","number");
    //SchemaFixName("Wyckoff_letters","","strings");
    //SchemaFixName("Wyckoff_letters_orig","","strings");
    //SchemaFixName("Wyckoff_multiplicities","","numbers");
    //SchemaFixName("Wyckoff_multiplicities_orig","","numbers");
    //SchemaFixName("Wyckoff_site_symmetries","","strings");
    //SchemaFixName("Wyckoff_site_symmetries_orig","","strings");
  }
}

namespace init {

  //ME20220208 - Initialize internal schema, which contain keywords that
  //are inside the SQLite database, but are not served to the public.
  uint InitSchemaInternal(bool INIT_VERBOSE) {
    // DECLARATIONS
    bool LDEBUG=(FALSE || XHOST.DEBUG || INIT_VERBOSE);
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitSchemaInternal: [BEGIN]" << endl;

    uint nschema = 0;

    // schema is CAPITAL, content is not necessarily
    XHOST.vschema_internal.push_attached("SCHEMA::NAME:ALLOY", "alloy");
    XHOST.vschema_internal.push_attached("SCHEMA::UNIT:ALLOY", "");
    XHOST.vschema_internal.push_attached("SCHEMA::TYPE:ALLOY", "string");
    nschema++;

    if(LDEBUG) cerr << "nschema=" << nschema << endl;
    if(LDEBUG) cerr << "AFLOW V(" << string(AFLOW_VERSION) << ") init::InitSchemaInternal: [END]" << endl;

    return nschema;
  }

}

namespace init {

  //getSchemaKeys/////////////////////////////////////////////////////////////
  // Returns the keys from the AFLOW schema.
  // Adapted from AflowDB
  vector<string> getSchemaKeys(const aurostd::xoption& vschema) {
    vector<string> keys;
    string key = "";
    for (uint i = 0, n = vschema.vxsghost.size(); i < n; i += 2) {
      if(vschema.vxsghost[i].find("SCHEMA::NAME:") != string::npos) {  //CO20200520
        key = aurostd::RemoveSubString(vschema.vxsghost[i], "SCHEMA::NAME:");
        // schema keys are upper case
        keys.push_back(aurostd::toupper(key));
      }
    }
    return keys;
  }

  //getSchemaNames/////////////////////////////////////////////////////////////
  // CO20200520
  // Returns the names of the AFLOW schema keys
  vector<string> getSchemaNames(const aurostd::xoption& vschema) {
    vector<string> keys;
    for (uint i = 0, n = vschema.vxsghost.size(); i < n; i += 2) {
      if(vschema.vxsghost[i].find("SCHEMA::NAME:") != string::npos) {
        keys.push_back(aurostd::RemoveSubString(vschema.vxsghost[i + 1], "SCHEMA::NAME:"));
      }
    }
    return keys;
  }

  //getSchemaTypes/////////////////////////////////////////////////////////////
  // Gets the data types of the schema keys.
  // Adapted from AflowDB
  vector<string> getSchemaTypes(const aurostd::xoption& vschema) {
    return getSchemaTypes(vschema, getSchemaKeys(vschema));
  }

  vector<string> getSchemaTypes(const aurostd::xoption& vschema, const vector<string>& keys) {
    uint nkeys = keys.size();
    vector<string> types(nkeys);
    string type = "";
    for (uint k = 0; k < nkeys; k++) {
      types.push_back(vschema.getattachedscheme("SCHEMA::TYPE:" + aurostd::toupper(keys[k])));
    }
    return types;
  }
}


// **************************************************************************

#endif

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
