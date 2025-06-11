// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef AUROSTD_XHOST_H
#define AUROSTD_XHOST_H

#include <sstream>
#include <string>
#include <sys/types.h>
#include <vector>

#include "AUROSTD/aurostd_xoption.h"

#define XPID XHOST.sPID
#define XPGID XHOST.sPGID
#define XTID XHOST.sTID

#define XHOST_vGlobal_MAX                              256
#define XHOST_aflowlib_icsd                            XHOST.vGlobal_string.at(23)
#define XHOST_aflowlib_lib0                            XHOST.vGlobal_string.at(24)
#define XHOST_aflowlib_lib1                            XHOST.vGlobal_string.at(25)
#define XHOST_aflowlib_lib2                            XHOST.vGlobal_string.at(26)
#define XHOST_aflowlib_lib3                            XHOST.vGlobal_string.at(27)
#define XHOST_aflowlib_lib4                            XHOST.vGlobal_string.at(28)
#define XHOST_aflowlib_lib5                            XHOST.vGlobal_string.at(29)
#define XHOST_aflowlib_lib6                            XHOST.vGlobal_string.at(30)
#define XHOST_aflowlib_lib7                            XHOST.vGlobal_string.at(31)
#define XHOST_aflowlib_lib8                            XHOST.vGlobal_string.at(32)
#define XHOST_aflowlib_lib9                            XHOST.vGlobal_string.at(33)
//#define XHOST_AUID                                     XHOST.vGlobal_string.at(34)
//#define XHOST_AURL                                     XHOST.vGlobal_string.at(35)
//#define XHOST_LOOP                                     XHOST.vGlobal_string.at(36)
#define XHOST_Library_ICSD_ALL                         XHOST.vGlobal_string.at(38)
//  string Library_ICSD_ALL;          // the complete library

#define XHOST_vLIBS XHOST.vvLIBS
#define XHOST_vAURL XHOST.vvLIBS.at(0)
#define XHOST_vAUID XHOST.vvLIBS.at(1)
#define XHOST_vLOOP XHOST.vvLIBS.at(2)
#define XHOST_LIBRARY_JSONL                            XHOST.vGlobal_string.at(3)

#define vVASP_POTCAR_DIRECTORIES                       XHOST.vvGlobal_string.at(4)
#define vAFLOW_LIBRARY_DIRECTORIES                     XHOST.vvGlobal_string.at(5)
#define vAFLOW_PROJECTS_DIRECTORIES                    XHOST.vvGlobal_string.at(6)
#define XHOST_vLibrary_ICSD                            XHOST.vvGlobal_string.at(7)
#define XHOST_vLibrary_ICSD_ALL                        XHOST.vvGlobal_string.at(8)
#define XHOST_Library_CALCULATED_ICSD_LIB              XHOST.vvGlobal_string.at(9)
#define XHOST_Library_CALCULATED_ICSD_RAW              XHOST.vvGlobal_string.at(10)
#define XHOST_Library_CALCULATED_LIB0_LIB              XHOST.vvGlobal_string.at(11)
#define XHOST_Library_CALCULATED_LIB0_RAW              XHOST.vvGlobal_string.at(12)
#define XHOST_Library_CALCULATED_LIB1_LIB              XHOST.vvGlobal_string.at(13)
#define XHOST_Library_CALCULATED_LIB1_RAW              XHOST.vvGlobal_string.at(14)
#define XHOST_Library_CALCULATED_LIB2_LIB              XHOST.vvGlobal_string.at(15)
#define XHOST_Library_CALCULATED_LIB2_RAW              XHOST.vvGlobal_string.at(16)
#define XHOST_Library_CALCULATED_LIB3_LIB              XHOST.vvGlobal_string.at(17)
#define XHOST_Library_CALCULATED_LIB3_RAW              XHOST.vvGlobal_string.at(18)
#define XHOST_Library_CALCULATED_LIB4_LIB              XHOST.vvGlobal_string.at(19)
#define XHOST_Library_CALCULATED_LIB4_RAW              XHOST.vvGlobal_string.at(20)
#define XHOST_Library_CALCULATED_LIB5_LIB              XHOST.vvGlobal_string.at(21)
#define XHOST_Library_CALCULATED_LIB5_RAW              XHOST.vvGlobal_string.at(22)
#define XHOST_Library_CALCULATED_LIB6_LIB              XHOST.vvGlobal_string.at(23)
#define XHOST_Library_CALCULATED_LIB6_RAW              XHOST.vvGlobal_string.at(24)
#define XHOST_Library_CALCULATED_LIB7_LIB              XHOST.vvGlobal_string.at(25)
#define XHOST_Library_CALCULATED_LIB7_RAW              XHOST.vvGlobal_string.at(26)
#define XHOST_Library_CALCULATED_LIB8_LIB              XHOST.vvGlobal_string.at(27)
#define XHOST_Library_CALCULATED_LIB8_RAW              XHOST.vvGlobal_string.at(28)
#define XHOST_Library_CALCULATED_LIB9_LIB              XHOST.vvGlobal_string.at(29)
#define XHOST_Library_CALCULATED_LIB9_RAW              XHOST.vvGlobal_string.at(30)

//  vector<string> vLibrary_ICSD;     // ordered by #species
//  vector<string> vLibrary_ICSD_ALL; // line by line

// all the at(N) need to be sequetial !!!
#define XHOST_ElectronStoppingPower_txt                XHOST.vGlobal_string.at(82)
#define XHOST_PhotonCrossSection_txt                   XHOST.vGlobal_string.at(83)
#define XHOST_PhotonStoppingPower_txt                  XHOST.vGlobal_string.at(84)
#define XHOST_ICSD_List_txt                            XHOST.vGlobal_string.at(85)
#define XHOST_AFLOW_PSEUDOPOTENTIALS                   XHOST.vGlobal_string.at(86)
#define XHOST_AFLOW_PSEUDOPOTENTIALS_TXT               XHOST.vGlobal_string.at(87)
#define XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT          XHOST.vGlobal_string.at(88)
#define XHOST_f144468a7ccc2d3a72ba44000715efdb         XHOST.vGlobal_string.at(90)

// LOADENTRIES DEFAULTS
#define _AFLOW_LIB_MAX_ 10                             //LIB11 does not exist yet, modify accordingly

#define XHOST_LIBRARY_LIB0                             XHOST.vGlobal_uint.at(0)
#define XHOST_LIBRARY_LIB1                             XHOST.vGlobal_uint.at(1)
#define XHOST_LIBRARY_LIB2                             XHOST.vGlobal_uint.at(2)
#define XHOST_LIBRARY_LIB3                             XHOST.vGlobal_uint.at(3)
#define XHOST_LIBRARY_LIB4                             XHOST.vGlobal_uint.at(4)
#define XHOST_LIBRARY_LIB5                             XHOST.vGlobal_uint.at(5)
#define XHOST_LIBRARY_LIB6                             XHOST.vGlobal_uint.at(6)
#define XHOST_LIBRARY_LIB7                             XHOST.vGlobal_uint.at(7)
#define XHOST_LIBRARY_LIB8                             XHOST.vGlobal_uint.at(8)
#define XHOST_LIBRARY_LIB9                             XHOST.vGlobal_uint.at(9)
#define XHOST_LIBRARY_ICSD                             XHOST.vGlobal_uint.at(10)
#define XHOST_LIBRARY_AUID                             XHOST.vGlobal_uint.at(11)

// --------------------------------------------------------------------------
// this is a container of general global choices
class _XHOST {
  public:
    // constructor destructor                         // constructor/destructor
    _XHOST();                                         // default, just allocate
    ~_XHOST();                                        // kill everything
    // _XHOST(const _XHOST& b);                          // constructor copy
    const _XHOST& operator=(const _XHOST &b);         // copy
    // BOOT
    int PGID,PID,TID;                // aflow_init.cpp  PID/TID number  //CO20200508 //SD20220329 PGID number
    std::ostringstream ostrPGID,ostrPID,ostrTID; // aflow_init.cpp  PID/TID in ostringstream... //CO20200508
    std::string sPGID,sPID,sTID;           // aflow_init.cpp  [PID=12345678]  [TID=12345678]
    bool showPGID,showPID,showTID;       // aflow_init.cpp  check if --showPID
    // machinery
    bool QUIET;           //CO20220630 - can be overridden by LOGGER_WHITELIST/_BLACKLIST, modifiable within functions
    bool QUIET_GLOBAL;    //CO20220630 - exclusively for --quiet (headless server), do not set/unset inside code (GLOBAL SILENCE)
    bool QUIET_CERR;      //CO20220630 - silences cerr exclusively  // extra quiet SC20210617
    bool QUIET_COUT;      //CO20220630 - silences cout exclusively  // extra quiet SC20210617
    bool TEST,DEBUG,MPI;
    std::vector<std::string> LOGGER_WHITELIST;  //HE+ME20220305 - for logging
    std::vector<std::string> LOGGER_BLACKLIST;  //HE+ME20220305 - for logging
    bool GENERATE_AFLOWIN_ONLY; //CT20180719
    bool POSTPROCESS; //CO20200624 - generic postprocessing, including --lib2raw and --lib2lib
    bool ARUN_POSTPROCESS; //CT20181212 - this is for the --postprocess flag needed for AEL/AGL, can be extended to other modules too
    bool AVOID_RUNNING_VASP; //CO20200624
    bool PSEUDOPOTENTIAL_GENERATOR; //SC20200327
    // HARDWARE/SOFTWARE
    std::string hostname,machine_type,tmpfs,user,group,home,shell,progname;
    std::string Find_Parameters;
    bool sensors_allowed;
    // ARGUMENTS
    std::vector<std::string> argv;          // argv of line command
    // SERVERS
    std::string AFLOW_MATERIALS_SERVER,AFLOW_WEB_SERVER;
    long double RAM,RAM_MB,RAM_GB;
    int CPU_Cores;
    int CPU_active;  //ME20220130
    std::string CPU_Model;
    std::string CPU_MHz;
    std::vector<double> vTemperatureCore;
    long double Time_starting,Time_now;
    long int Date;
    std::string Day,Month,Year;
    std::string Copyright_Years; // =string("2003-YEAR_FROM_DATE");
    // MULTHREADS
    // bool PTHREADS_FLAG;        // run pthread YES/NO
    // int  PTHREADS_MAX;         // how many MAX threads I can use  default or --np
    // int PTHREADS_RUNNING;      // how many threads are actually running
    // vector<pthread_t> thread;  // the actual thread
    // vector<int> iret;          // the thread runnings
    // vector<bool> thread_busy;  // is the thread busy
    // COMMANDS
    std::vector<std::string> vcmd;
    // RAM CHECK
    double maxmem;
    // FUNCTIONS
    std::string command(const std::string& command);
    bool is_command(const std::string& command);
    // AFLOW STUFF
    // vflag_aflow.flag("LOOP");
    // vflag_aflow.flag("CLEAN");
    // vflag_aflow.isflag*"XCLEAN");
    bool AFLOW_RUNDIRflag;
    bool AFLOW_MULTIflag;
    bool AFLOW_RUNXflag;
    uint AFLOW_RUNXnumber;
    // QUEQUE STUFF
    bool is_PBS;int PBS_NUM_PPN,PBS_NNODES;
    bool is_SLURM;int SLURM_CPUS_ON_NODE,SLURM_NNODES,SLURM_NTASKS;
    bool is_MACHINE_FULTON_MARYLOU; // some flags
    // Library_CALCULATED*
    std::vector<uint>   vGlobal_uint;      // parameters uint
    std::vector<std::string> vGlobal_string;    // parameters as strings
    std::vector<std::vector<std::string> > vvGlobal_string; // parameters as vector strings
    std::vector<std::vector<std::string> > vvLIBS; // parameters as vector strings
    // vector<string> vLibrary_ICSD;     // ordered by #species (needs to be allocated)
    // vector<string> vLibrary_ICSD_ALL; // line by line
    // string Library_ICSD_ALL;          // the complete library
    // vector<string> vVASP_POTCAR_DIRECTORIES;
    // vector<string> vAFLOW_LIBRARY_DIRECTORIES;
    // vector<string> vAFLOW_PROJECTS_DIRECTORIES;
    // AFLOW flags/options
    aurostd::xoption vflag_aflow;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_pflow;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_outreach;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_control;  // argv/argc options following the xoption structure
    aurostd::xoption vschema;        // keywords, names, units etc etc
    aurostd::xoption vschema_internal;  //ME20220208
    // USUAL COMMANDS
    std::vector<std::string> vcat; //     cat, bzcat, xzcat, gzcat
    std::vector<std::string> vext; //      "",  .bz2,   .xz,   .gz
    std::vector<std::string> vzip; //      "", bzip2,    xz,  gzip
    // AFLOWRC
    std::string aflowrc_filename;
    std::string aflowrc_content;
    std::vector<std::string> vaflowrc;
    aurostd::xoption adefault;            // default  xoption
    // AFLOWSYM
    bool SKEW_TEST; //DX20171019
    double SKEW_TOL; //DX20171019
    // xstructure
    bool READ_SPIN_FROM_ATOMLABEL; //SD20220316
    // WEB MODE
    //[CO20200404 - overload with --www]bool WEB_MODE;  //CO20190401
  private:                                                //
    void free();                                           // free space
    void copy(const _XHOST& b);                            //
    void clear();                                          // free space
};

// max is 128
extern _XHOST XHOST; // this will be global

#endif //AUROSTD_XHOST_H
