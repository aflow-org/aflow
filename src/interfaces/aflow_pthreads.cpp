// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo
// contains routines to run PTHREADS

#ifndef _AFLOW_PTHREADS_CPP
#define _AFLOW_PTHREADS_CPP

#include "interfaces/aflow_pthreads.h"

#include <cstddef>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <pthread.h>
#include <unistd.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_argv.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"

#include "aflow.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "aflow_xthread.h"
#include "flow/aflow_xclasses.h"

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

// #define  AFLOW_PTHREADS_MULTISH_PREEMPTIVE_
#ifndef AFLOW_MULTITHREADS_ENABLE
#define AFLOW_PTHREADS_MULTISH_TIMESHARING_
#endif
// #define  AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_
// #define  AFLOW_PTHREADS::MULTISH_TIMESHARING_CONCURRENT_

#define _KBIN_SEEK_THREAD_SLEEP_ 33
#define _KBIN_START_THREAD_SLEEP_ 2
#define _KBIN_FLUSH_THREAD_SLEEP_ 1

// #define _PTHREAD_FLUSH_TIME_ 1

namespace AFLOW_PTHREADS {
  bool FLAG;                                     // run pthread YES/NO
  int MAX_PTHREADS;                              // how many MAX threads I can use  default or --np
  int RUNNING;                                   // how many threads are actually running
  pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  int viret[MAX_ALLOCATABLE_PTHREADS];           // the thread runnings
  bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
  bool MULTISH_TIMESHARING_SEQUENTIAL_ = false;
  bool MULTISH_TIMESHARING_CONCURRENT_ = true;
} // namespace AFLOW_PTHREADS

#define CPU_File string("/proc/cpuinfo")
#define CPU_String string("cpu MHz")

// ***************************************************************************
// AFLOW_PTHREADS::GetTotalCPUs
// ***************************************************************************
// This function returns the max number of CPUS by interrogating as
// cat /proc/cpuinfo | grep -c "cpu MHz"
// if the file is not found, it returns CPUs=1
namespace AFLOW_PTHREADS {
  int GetTotalCPUs() {
    int CPU_Cores = sysconf(_SC_NPROCESSORS_ONLN);
    if (CPU_Cores < 1) {
      CPU_Cores = 1;
    }
    return CPU_Cores;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Check_Threads
// **************************************************************************
// This function checks the input argv and set up the proper multithread
// parameters (SC Dec07)
namespace AFLOW_PTHREADS {
  bool Check_Threads(vector<string> argv, const bool& VERBOSE) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (VERBOSE) {
      ;
    } // dummy load
    AFLOW_PTHREADS::FLAG = false;
    AFLOW_PTHREADS::MAX_PTHREADS = 1;
    const bool fnp = XHOST.vflag_control.flag("NUM_THREADS"); // aurostd::args2attachedflag(argv,"--np=");

    const bool fnpmax = aurostd::args2flag(argv, "--npmax");
    const bool multi_sh = aurostd::args2flag(argv, "--multi=sh|--multi=sh");
    if (!fnp && !fnpmax) {
      AFLOW_PTHREADS::MAX_PTHREADS = 1;
    }
    if (fnp && !fnpmax) {
      AFLOW_PTHREADS::MAX_PTHREADS = aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("NUM_THREADS"));
    } // aurostd::args2attachedutype<int>(argv,"--np=",0);; //SC20200319
    if (!fnp && fnpmax) {
      AFLOW_PTHREADS::MAX_PTHREADS = AFLOW_PTHREADS::GetTotalCPUs();
    }
    if (fnp && fnpmax) {
      AFLOW_PTHREADS::MAX_PTHREADS = aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("NUM_THREADS"));
    } // aurostd::args2attachedutype<int>(argv,"--np=",0);; //SC20200319
    if (multi_sh && !fnp && !fnpmax) {
      AFLOW_PTHREADS::MAX_PTHREADS = AFLOW_PTHREADS::GetTotalCPUs();
    }

    if (AFLOW_PTHREADS::MAX_PTHREADS > 1) {
      AFLOW_PTHREADS::FLAG = true;
      if (AFLOW_PTHREADS::MAX_PTHREADS > MAX_ALLOCATABLE_PTHREADS) {
        AFLOW_PTHREADS::MAX_PTHREADS = MAX_ALLOCATABLE_PTHREADS;
      }
      if (LDEBUG) {
        cerr << "AAAAA  AFLOW THREADED VERSION  threads=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
      }
    }
    if (AFLOW_PTHREADS::MAX_PTHREADS <= 1) {
      AFLOW_PTHREADS::FLAG = false;
      AFLOW_PTHREADS::MAX_PTHREADS = 1;
      if (LDEBUG) {
        cerr << "AAAAA  AFLOW SERIAL VERSION threads=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
      }
    }
    //  AFLOW_PTHREADS::FLAG=true;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " fnp=" << fnp << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " fnpmax=" << fnpmax << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " multi_sh=" << multi_sh << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << endl;
    }
    //  if(LDEBUG) throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_); // for debug
    return AFLOW_PTHREADS::FLAG;
  }
} // namespace AFLOW_PTHREADS

namespace AFLOW_PTHREADS {
  bool Check_Threads_WrapperNP(vector<string> argv, uint size_to_check, const bool& VERBOSE) {
    ostringstream aus;
    AFLOW_PTHREADS::FLAG = true;
    AFLOW_PTHREADS::Check_Threads(argv, VERBOSE);
    if (AFLOW_PTHREADS::MAX_PTHREADS > (int) size_to_check) {
      if (VERBOSE) {
        aus << "AFLOW multi_XXXX (" << string(AFLOW_VERSION) << "): WARNING (threads>commands) => threads=commands" << endl;
      }
      AFLOW_PTHREADS::MAX_PTHREADS = (int) size_to_check;
    }
    if (AFLOW_PTHREADS::MAX_PTHREADS >= 2) {
      if (VERBOSE) {
        aus << "AFLOW multi_XXXX (" << string(AFLOW_VERSION) << "): threads=" << AFLOW_PTHREADS::MAX_PTHREADS << "  -  commands=" << size_to_check << endl;
      }
    }
    if (AFLOW_PTHREADS::MAX_PTHREADS <= 1) {
      if (VERBOSE) {
        aus << "AFLOW multi_XXXX (" << string(AFLOW_VERSION) << "): serial  -  commands=" << size_to_check << endl;
      }
    }
    //  if(VERBOSE)
    aurostd::PrintMessageStream(aus, XHOST.QUIET);

    return AFLOW_PTHREADS::FLAG;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Clean_Threads
// **************************************************************************
// This function clears the busy-ness of all the pthreads flags
namespace AFLOW_PTHREADS {
  void Clean_Threads() {
    for (uint ithread = 0; ithread < MAX_ALLOCATABLE_PTHREADS; ithread++) {         // clean threads
      AFLOW_PTHREADS::vpthread_busy[ithread] = false;          // clean threads
    }
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::No_Threads
// **************************************************************************
// This function removes threads
namespace AFLOW_PTHREADS {
  void No_Threads() {
    AFLOW_PTHREADS::FLAG = false;
    AFLOW_PTHREADS::MAX_PTHREADS = 1;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Available_Free_Threads
// **************************************************************************
// This function return true and a free PTHREAD index if a free thread is
// available
namespace AFLOW_PTHREADS {
  bool Available_Free_Threads(int& fthread) {
    bool free_thread = false;
    fthread = -1;
    for (int ithread = 0; ithread < AFLOW_PTHREADS::MAX_PTHREADS; ithread++) {
      if (AFLOW_PTHREADS::vpthread_busy[ithread] == false) {
        free_thread = true;
        fthread = ithread;
      }
    }
    return free_thread;
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// AFLOW_PTHREADS::Wait_Available_Free_Threads
// **************************************************************************
// This function return true and a free PTHREAD index if a free thread is
// available
namespace AFLOW_PTHREADS {
  bool Wait_Available_Free_Threads(int& fthread, const double& pthread_wait, const bool& VERBOSE) {
    ostringstream aus;
    bool free_thread = false;
    while (!free_thread) {            // waiting for free thread to start
      free_thread = AFLOW_PTHREADS::Available_Free_Threads(fthread);
      if (!free_thread) {              // do something else !
        if (VERBOSE) {
          aus << "MMMMM  Aflow: MULTI-THREADED: Waiting for free threads: " << pthread_wait << " seconds " << endl;
        }
        if (VERBOSE) {
          aurostd::PrintMessageStream(aus, false);
        }
        aurostd::Sleep((int) pthread_wait);
      }
      if (free_thread) {
        if (VERBOSE) {
          aus << "MMMMM  Aflow: MULTI-THREADED: Found free thread  fthread=" << fthread << " - " << endl;
        }
        if (VERBOSE) {
          aurostd::PrintMessageStream(aus, false);
        }
      }
    }
    return free_thread;
  }
} // namespace AFLOW_PTHREADS

namespace AFLOW_PTHREADS {
  bool Wait_Available_Free_Threads(int& fthread, const bool& VERBOSE) {
    return AFLOW_PTHREADS::Wait_Available_Free_Threads(fthread, _KBIN_SEEK_THREAD_SLEEP_, VERBOSE);
  }
} // namespace AFLOW_PTHREADS

// **************************************************************************
// KBIN::RUN_Directory_PTHREADS
// **************************************************************************
// Interfaces for KBIN_RUN_Directory

namespace KBIN {
  typedef struct {
    _aflags* paflags;     // FOR KBIN (ALL)
    int itbusy;      // FOR KBIN (ALL)
    bool VERBOSE;     // FOR KBIN (ALL)
    string command;     // FOR MULTISH PREEMPTIVE
    int ITHREAD;     // FOR MULTISH_TIMESHARING
    int THREADS_MAX; // FOR MULTISH_TIMESHARING
    deque<string>* dcmds; // FOR MULTISH_TIMESHARING
  } _threaded_params;
} // namespace KBIN

KBIN::_threaded_params params[MAX_ALLOCATABLE_PTHREADS];
_aflags taflags[MAX_ALLOCATABLE_PTHREADS];
pthread_mutex_t mutex_PTHREAD = PTHREAD_MUTEX_INITIALIZER;

namespace KBIN {
  void RUN_Directory_PTHREADS(_aflags& aflags) {
    stringstream message;

    const int ithread = aflags.AFLOW_PTHREADS_NUMBER;
    if (ithread < 0) {
      message << "ithread<0  ithread=" << ithread;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    if (ithread >= MAX_ALLOCATABLE_PTHREADS) {
      message << "ithread>=MAX_ALLOCATABLE_PTHREADS  ithread=" << ithread;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    if (ithread >= AFLOW_PTHREADS::MAX_PTHREADS) {
      message << "ithread>=AFLOW_PTHREADS::MAX_PTHREADS  ithread=" << ithread;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    taflags[ithread] = aflags;
    params[ithread].paflags = &taflags[ithread];
    params[ithread].command = "";// command;
    params[ithread].itbusy = ithread;
    AFLOW_PTHREADS::viret[ithread] = pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]), nullptr, KBIN::_threaded_interface_RUN_Directory, (void*) &params[ithread]);
    aurostd::Sleep(_KBIN_START_THREAD_SLEEP_);
  }
} // namespace KBIN

namespace KBIN {
  void* _threaded_interface_RUN_Directory(void* ptr) {
    KBIN::_threaded_params* pparams;
    pparams = (KBIN::_threaded_params*) ptr;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = true;
    AFLOW_PTHREADS::RUNNING++;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP_);
    aurostd::execute(XHOST.command("aflow") + " --run=1 --DIRECTORY=" + (*pparams->paflags).Directory);  // run it OUTSIDE
    //  KBIN_RUN_Directory((*pparams->paflags)); // RUN IT INSIDE
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = false;
    AFLOW_PTHREADS::RUNNING--;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP_);
    return nullptr;
  }
} // namespace KBIN

// **************************************************************************
// FUNCTION AFLOW_PTHREADS_MULTISH_PREEMPTIVE_
// **************************************************************************
// Interfaces for AFLOW_PTHREADS_MULTISH_PREEMPTIVE_

#ifdef AFLOW_PTHREADS_MULTISH_PREEMPTIVE_
#warning "aflow_pthreads.cpp with AFLOW_PTHREADS_MULTISH_PREEMPTIVE_"

namespace KBIN {
  void* _threaded_interface_MULTIRUN_sh(void* ptr) {
    KBIN::_threaded_params* pparams;
    pparams = (KBIN::_threaded_params*) ptr;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = true;
    AFLOW_PTHREADS::RUNNING++;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP);
    cout << pparams->itbusy << "  - " << pparams->command << endl;
    aurostd::execute(pparams->command);
    // KBIN_RUN_Directory((*pparams->paflags));
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = false;
    AFLOW_PTHREADS::RUNNING--;
    aurostd::Sleep(_KBIN_FLUSH_THREAD_SLEEP);
    return NULL;
  }
} // namespace KBIN

namespace KBIN {
  void MULTIRUN_PTHREADS(_aflags& aflags, string command) {
    stringstream message;

    int ithread = aflags.AFLOW_PTHREADS_NUMBER;
    if (ithread < 0) {
      message << "ithread<0  ithread=" << ithread;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    if (ithread >= MAX_ALLOCATABLE_PTHREADS) {
      message << "ithread>=MAX_ALLOCATABLE_PTHREADS  ithread=" << ithread;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    if (ithread >= AFLOW_PTHREADS::MAX_PTHREADS) {
      message << "ithread>=AFLOW_PTHREADS::MAX_PTHREADS  ithread=" << ithread;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_BOUNDS_);
    }
    taflags[ithread] = aflags;
    params[ithread].paflags = &taflags[ithread];
    params[ithread].command = command;
    params[ithread].itbusy = ithread;
    AFLOW_PTHREADS::viret[ithread] = pthread_create(&(AFLOW_PTHREADS::vpthread[ithread]), NULL, KBIN::_threaded_interface_MULTIRUN_sh, (void*) &params[ithread]);
    aurostd::Sleep(_KBIN_START_THREAD_SLEEP_);
  }
} // namespace KBIN

namespace AFLOW_PTHREADS {
  bool MULTI_sh(const vector<string>& argv) {
    stringstream message;
    ostringstream aus;
    _aflags aflags;
    string file_name = XHOST.vflag_control.getattachedscheme("FILE");
    bool free_thread;
    int ithread = 0;
    bool VERBOSE = false;

    AFLOW_PTHREADS::FLAG = true;
    AFLOW_PTHREADS::MAX_PTHREADS = PTHREAD_DEFAULT; // safety...

    if (!aurostd::FileExist(file_name)) {
      message << "FILE_NOT_FOUND = " << file_name;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    if (aurostd::FileEmpty(file_name)) {
      message << "FILE_EMPTY = " << file_name;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    aus << "MMMMM Loading File = " << file_name << endl;
    aurostd::PrintMessageStream(aus, XHOST.QUIET);
    vector<string> vcmds;
    vcmds.clear();
    aurostd::file2vectorstring(file_name, vcmds);
    aus << "MMMMM Loaded Lines = " << vcmds.size() << endl;
    aurostd::PrintMessageStream(aus, XHOST.QUIET);

    AFLOW_PTHREADS::Clean_Threads();                                  // clean threads

    for (size_t i = 0; i < vcmds.size(); i++) {
      // loop this
      free_thread = AFLOW_PTHREADS::Wait_Available_Free_Threads(ithread, 0, VERBOSE);        // WAIT A WHILE !!
      if (free_thread) {
        //  aus << "MMMMM Aflow: Found subdirectory to run " << aflags.Directory<< " - " << XHOST.hostname << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
        aus << "MMMMM Aflow: MULTI-THREADED: Starting  pthread_free=" << ithread << "  pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " vline=" << i << " - " << XHOST.hostname << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        aflags.AFLOW_PTHREADS_NUMBER = ithread;
        AFLOW_PTHREADS::vpthread_busy[ithread] = true;
        KBIN::MULTIRUN_PTHREADS(aflags, vcmds[i]);
        AFLOW_PTHREADS::vpthread_busy[ithread] = false;
      }
      //    cout << free_thread << " " << ithread << endl;
    }

    aus << "MMMMM  Aflow: MULTI-THREADED: FLUSHING PTHREADS - " << XHOST.hostname << endl;
    aurostd::PrintMessageStream(aus, XHOST.QUIET);
    for (ithread = 0; ithread < AFLOW_PTHREADS::MAX_PTHREADS; ithread++)
      if (AFLOW_PTHREADS::vpthread_busy[ithread] == true) {
        aus << "MMMMM  Aflow: MULTI-THREADED: Flushing   pthread=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << XHOST.hostname << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        pthread_join(thread[ithread], NULL);
      }
    return true;
  }
} // namespace AFLOW_PTHREADS

#endif //  AFLOW_PTHREADS_MULTISH_PREEMPTIVE_

// **************************************************************************
// FUNCTION AFLOW_PTHREADS_MULTISH_TIMESHARING_
// **************************************************************************
// Interfaces for AFLOW_PTHREADS_MULTISH_TIMESHARING_

#ifdef AFLOW_PTHREADS_MULTISH_TIMESHARING_
#warning "aflow_pthreads.cpp with AFLOW_PTHREADS::MULTISH_TIMESHARING_"

namespace AFLOW_PTHREADS {
  void* _threaded_COMMANDS(void* ptr) {
    stringstream message;
    if (AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_) {
      // cerr << XPID << "AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_ AFLOW_PTHREADS::_threaded_COMMANDS" << endl;
      KBIN::_threaded_params* pparams = (KBIN::_threaded_params*) ptr;
      string command;
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = true;
      AFLOW_PTHREADS::RUNNING++;
      if ((pparams->VERBOSE)) {
        pthread_mutex_lock(&mutex_PTHREAD);
        cout << __AFLOW_FUNC__ << " " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;
        pthread_mutex_unlock(&mutex_PTHREAD);
      }
      for (size_t ithread = (pparams->ITHREAD); ithread < (*pparams->dcmds).size(); ithread += (pparams->THREADS_MAX)) {
        if ((pparams->VERBOSE)) {
          pthread_mutex_lock(&mutex_PTHREAD);
          cout << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << " " << ithread << " " << (*pparams->dcmds).at(ithread) << endl;
          pthread_mutex_unlock(&mutex_PTHREAD);
        }
        command = (*pparams->dcmds).at(ithread);
        aurostd::execute(command);
        pthread_mutex_lock(&mutex_PTHREAD);
        cout << command << endl;
        cout.flush();
        pthread_mutex_unlock(&mutex_PTHREAD);
        // if((pparams->VERBOSE)) {cout << char('A'+(pparams->ITHREAD));cout.flush();}
      }
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = false;
      AFLOW_PTHREADS::RUNNING--;
      // aurostd::Sleep(_KBIN_START_THREAD_SLEEP_);
      return NULL;
    }
    if (AFLOW_PTHREADS::MULTISH_TIMESHARING_CONCURRENT_) { // it bombs and I do not know why...
      bool FRONT = true;
      KBIN::_threaded_params* pparams = (KBIN::_threaded_params*) ptr;
      string command;
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = true;
      AFLOW_PTHREADS::RUNNING++;
      while ((*pparams->dcmds).size() > 0) {
        pthread_mutex_lock(&mutex_PTHREAD);
        if (FRONT) {
          command = (*pparams->dcmds).at(0);
          (*pparams->dcmds).pop_front();
        }  // from the front
        if (!FRONT) {
          command = (*pparams->dcmds).at((*pparams->dcmds).size() - 1);
          (*pparams->dcmds).pop_back();
        }  // from the back
        pthread_mutex_unlock(&mutex_PTHREAD);
        if ((pparams->VERBOSE)) {
          pthread_mutex_lock(&mutex_PTHREAD);
          cout << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << ": " << command << endl;
          pthread_mutex_unlock(&mutex_PTHREAD);
        }
        //    aurostd::Sleep((int) _KBIN_START_THREAD_SLEEP_);
        //    cerr << "[1]" << endl;
        //   cout << command << endl;cout.flush();
        aurostd::execute_thread_safe(command); // HE20240219 changed to execute_thread_safe
        // if((pparams->VERBOSE)) {cout << char('A'+(pparams->ITHREAD));cout.flush();}
      }
      AFLOW_PTHREADS::vpthread_busy[pparams->itbusy] = false;
      AFLOW_PTHREADS::RUNNING--;
      // cerr << "[2]" << endl;
      return NULL;
    }
    if (!AFLOW_PTHREADS::MULTISH_TIMESHARING_SEQUENTIAL_ && !AFLOW_PTHREADS::MULTISH_TIMESHARING_CONCURRENT_) {
      message << "You must specify MULTISH_TIMESHARING_SEQUENTIAL_ of MULTISH_TIMESHARING_CONCURRENT_";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_MISSING_);
    }
    return NULL;
  }
} // namespace AFLOW_PTHREADS

namespace AFLOW_PTHREADS {
  bool MULTI_sh(const vector<string>& argv) {
    stringstream message;
    ostringstream aus;
    _aflags aflags;
    string file_name = XHOST.vflag_control.getattachedscheme("FILE");
    if (file_name.empty() || file_name == "--f") file_name = argv.at(argv.size() - 1);
    //  bool free_thread;int ithread=0;
    bool VERBOSE = false;

    if (!aurostd::FileExist(file_name)) {
      message << "FILE_NOT_FOUND = " << file_name;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    if (aurostd::FileEmpty(file_name)) {
      message << "FILE_EMPTY = " << file_name;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    if (VERBOSE) {
      aus << "MMMMM  Loading File = " << file_name << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
    }
    vector<string> vcmds;
    vcmds.clear();
    aurostd::file2vectorstring(file_name, vcmds);

    AFLOW_PTHREADS::Check_Threads_WrapperNP(argv, vcmds.size(), VERBOSE);   // check treads from NP
    aurostd::multithread_execute(vcmds, AFLOW_PTHREADS::MAX_PTHREADS, VERBOSE);
    cout << endl;
    return true;
  }
} // namespace AFLOW_PTHREADS

namespace aurostd {
  bool multithread_execute(const deque<string>& dcmds_in, int _NUM_THREADS, bool VERBOSE) {
    deque<string> dcmds = dcmds_in;
    bool LDEBUG = true;
    int NUM_THREADS = _NUM_THREADS;                                          // SAFETY
    if ((int) dcmds.size() <= NUM_THREADS) NUM_THREADS = (uint) dcmds.size();   // SAFETY

    if (NUM_THREADS <= 1) {                                                   // run singular
      for (size_t i = 0; i < dcmds.size(); i++)                                     // run singular
        aurostd::execute(dcmds[i]);                                     // run singular
    }
    if (NUM_THREADS >= 2) {                                                   // multithread
      AFLOW_PTHREADS::FLAG = true;
      AFLOW_PTHREADS::MAX_PTHREADS = NUM_THREADS;  // prepare
      if (AFLOW_PTHREADS::MAX_PTHREADS > MAX_ALLOCATABLE_PTHREADS) AFLOW_PTHREADS::MAX_PTHREADS = MAX_ALLOCATABLE_PTHREADS; // check max

      if (LDEBUG) cerr << __AFLOW_FUNC__ << " _NUM_THREADS=" << _NUM_THREADS << endl;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " NUM_THREADS=" << NUM_THREADS << endl;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " MAX_ALLOCATABLE_PTHREADS=" << MAX_ALLOCATABLE_PTHREADS << endl;
      if (LDEBUG) cerr << __AFLOW_FUNC__ << " AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;

      _aflags aflags;
      AFLOW_PTHREADS::Clean_Threads();                                     // clean threads
      KBIN::_threaded_params params[MAX_ALLOCATABLE_PTHREADS];             // prepare
      for (int ithread = 0; ithread < AFLOW_PTHREADS::MAX_PTHREADS; ithread++) {  // prepare loop
        params[ithread].paflags = &aflags;                                   // prepare params
        params[ithread].ITHREAD = ithread;                                   // prepare params
        params[ithread].THREADS_MAX = AFLOW_PTHREADS::MAX_PTHREADS;                    // prepare params
        //    cerr << AFLOW_PTHREADS::MAX_PTHREADS << endl;
        params[ithread].dcmds = &dcmds;                                      // prepare params
        params[ithread].itbusy = ithread;                                    // prepare params
        params[ithread].VERBOSE = VERBOSE;                                   // prepare params
      }
      for (int ithread = 0; ithread < AFLOW_PTHREADS::MAX_PTHREADS; ithread++)                          // run threads
        AFLOW_PTHREADS::viret[ithread] = pthread_create(&AFLOW_PTHREADS::vpthread[ithread], NULL, AFLOW_PTHREADS::_threaded_COMMANDS, (void*) &params[ithread]); // run threads
      for (int ithread = 0; ithread < AFLOW_PTHREADS::MAX_PTHREADS; ithread++)                          // flush
        pthread_join(AFLOW_PTHREADS::vpthread[ithread], NULL);                                    // flush
    }
    return true;
  }
} // namespace aurostd

#endif //  AFLOW_PTHREADS_MULTISH_TIMESHARING_

// ***************************************************************************
// MultiThread Execute vectors/deque of Strings
// ***************************************************************************
// adding something to aurostd
namespace aurostd {
  bool multithread_execute(const deque<string>& cmds, int _NUM_THREADS, bool VERBOSE) {
    const bool LDEBUG = false;
    if (VERBOSE) {}
    if (LDEBUG) {
      std::cerr << "Commands to run:\n" << aurostd::joinWDelimiter(cmds, "\n") << std::endl;
    }

    std::function<void(const deque<string>::const_iterator&)> fn = [&](const deque<string>::const_iterator& it) { aurostd::execute_thread_safe(*it); }; // HE20240220 changed to execute_thread_safe
    xthread::xThread xt(_NUM_THREADS);
    xt.run(cmds, fn);
    return true;
  }
  bool multithread_execute(const vector<string>& vcmds, int _NUM_THREADS, bool VERBOSE) {
    int NUM_THREADS = _NUM_THREADS;                                          // SAFETY
    if ((int) vcmds.size() <= NUM_THREADS) {
      NUM_THREADS = (uint) vcmds.size();   // SAFETY
    }

    deque<string> dcmds;
    dcmds.clear();
    for (size_t i = 0; i < vcmds.size(); i++) {
      dcmds.push_back(vcmds[i]);         // copy
    }
    return aurostd::multithread_execute(dcmds, NUM_THREADS, VERBOSE);
  }
  bool multithread_execute(const deque<string>& vcommand, int NUM_THREADS) {
    return multithread_execute(vcommand, NUM_THREADS, false);
  }
  bool multithread_execute(const vector<string>& vcommand, int NUM_THREADS) {
    return multithread_execute(vcommand, NUM_THREADS, false);
  }
  bool multithread_execute(const deque<string>& vcommand) {
    AFLOW_PTHREADS::MAX_PTHREADS = AFLOW_PTHREADS::GetTotalCPUs();
    return multithread_execute(vcommand, AFLOW_PTHREADS::MAX_PTHREADS, false);
  }
  bool multithread_execute(const vector<string>& vcommand) {
    AFLOW_PTHREADS::MAX_PTHREADS = AFLOW_PTHREADS::GetTotalCPUs();
    return multithread_execute(vcommand, AFLOW_PTHREADS::MAX_PTHREADS, false);
  }
} // namespace aurostd

#ifdef AFLOW_MULTITHREADS_ENABLE

namespace AFLOW_PTHREADS {
  bool MULTI_sh(const vector<string>& argv) {
    stringstream message;
    ostringstream aus;
    const _aflags aflags;
    string file_name = XHOST.vflag_control.getattachedscheme("FILE");
    if (file_name.empty() || file_name == "--f") {
      file_name = argv.at(argv.size() - 1);
    }
    const bool VERBOSE = false;

    if (!aurostd::FileExist(file_name)) {
      message << "FILE_NOT_FOUND = " << file_name;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    if (aurostd::FileEmpty(file_name)) {
      message << "FILE_EMPTY = " << file_name;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    aus << "MMMMM Loading File = " << file_name << endl;
    aurostd::PrintMessageStream(aus, XHOST.QUIET);
    vector<string> vcmds;
    vcmds.clear();
    aurostd::file2vectorstring(file_name, vcmds);
    aus << "MMMMM Loaded Lines = " << vcmds.size() << endl;
    aurostd::PrintMessageStream(aus, XHOST.QUIET);
    return aurostd::multithread_execute(vcmds, KBIN::get_NCPUS(), VERBOSE);
  }
} // namespace AFLOW_PTHREADS

#endif

// **************************************************************************
// NOT MULTITHREAD BUT GOOD ENOUGH....

// ***************************************************************************
// sflow::KILL
// ***************************************************************************
namespace sflow {
  void KILL(string options) {
    vector<int> jobs;
    aurostd::StringCommasColumsVectorInt(options, jobs);
    vector<string> vcmds;

    stringstream aus;
    for (size_t i = 0; i < jobs.size(); i++) {
      if (XHOST.is_command("kill")) {
        aus.clear();
        aus.str(std::string());
        aus << XHOST.command("kill") << " -9 " << jobs[i] << endl;  // command = kill
        vcmds.push_back(aus.str());
      }
    }
    for (size_t i = 0; i < vcmds.size(); i++) {
      cout << "EXECUTING: " << vcmds[i];// << endl;
      aurostd::execute(vcmds[i]);
    }
  }
} // namespace sflow

// ***************************************************************************
// sflow::JUST
// ***************************************************************************
namespace sflow {
  void JUST(string options, istream& input, string mode) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    string strout;
    ostringstream osstemp;
    osstemp << input.rdbuf();
    strout = osstemp.str();

    if (mode == "JUSTAFTER" || mode == "AFTER") {
      if (tokens.size() != 1) {
        init::ErrorOption(options, "sflow::JUST", "aflow --justafter=string < something");
      }
      vector<string> vline;
      aurostd::string2vectorstring(strout, vline);
      bool found = false;
      for (size_t i = 0; i < vline.size(); i++) {
        if (found) {
          cout << vline[i] << endl;
        }
        if (!found) {
          found = aurostd::substring2bool(vline[i], tokens.at(0));
        }
      }
    }
    if (mode == "JUSTBEFORE" || mode == "BEFORE") {
      if (tokens.size() != 1) {
        init::ErrorOption(options, "sflow::JUST", "aflow --justbefore=string < something");
      }
      if (strout.find(tokens.at(0)) == string::npos) {
        cout << strout;
      }
      strout = strout.substr(0, strout.find(tokens.at(0)));
      strout = strout.substr(0, strout.find_last_of("\n") + 1);
      cout << strout;
    }

    if (mode == "JUSTBETWEEN" || mode == "BETWEEN") {
      if (tokens.size() > 2) {
        init::ErrorOption(options, "sflow::JUST", "aflow --justbetween=string_start[,string_stop] < something");
      }

      string strfind_from;
      string strfind_to;
      string avoid1_strfind_from;
      string avoid1_strfind_to;
      if (tokens.size() == 1) {
        strfind_from = "START." + tokens.at(0);
        strfind_to = "STOP." + tokens.at(0);
      }
      if (tokens.size() == 2) {
        strfind_from = tokens.at(0);
        strfind_to = tokens.at(1);
      }
      avoid1_strfind_from = "#" + strfind_from;
      avoid1_strfind_to = "#" + strfind_to;

      vector<string> vstranalyze;
      uint istart = 0;
      uint istop = 0;
      aurostd::string2vectorstring(osstemp.str(), vstranalyze);
      for (size_t i = 0; i < vstranalyze.size(); i++) {
        if (aurostd::substring2bool(vstranalyze[i], strfind_from)) {
          istart = i;
        }
        if (aurostd::substring2bool(vstranalyze[i], strfind_to)) {
          istop = i;
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " istart=" << istart << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " istop=" << istop << endl;
      }

      for (size_t i = 0; i < vstranalyze.size(); i++) {
        if (i > istart && i < istop) {
          cout << vstranalyze[i] << endl;
        }
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " END" << endl;
    }
    // throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Throw for debugging purposes.",_GENERIC_ERROR_);
  }
} // namespace sflow

// ***************************************************************************
// sflow::QDEL qdel bkill scancel
// ***************************************************************************
namespace sflow {
  void QDEL(string options, string cmd) {
    vector<int> jobs;
    aurostd::StringCommasColumsVectorInt(options, jobs);
    vector<string> vcmds;
    stringstream aus;
    for (size_t i = 0; i < jobs.size(); i++) {
      aus.clear();
      aus.str(std::string());
      aus << cmd << " " << jobs[i] << endl;  // cmd = qdel OR bkill
      vcmds.push_back(aus.str());
    }
    for (size_t i = 0; i < vcmds.size(); i++) {
      cout << "EXECUTING: " << vcmds[i];// << endl;
      aurostd::execute(vcmds[i]);
    }
  }
} // namespace sflow

namespace sflow {
  void QDEL(string options) {
    vector<int> jobs;
    aurostd::StringCommasColumsVectorInt(options, jobs);
    vector<string> vcmds;
    // if(aurostd::args2flag(argv,"--scancel")) {XHOST.is_command("qdel")=false;XHOST.is_command("bkill")=false;} // force
    // if(aurostd::args2flag(argv,"--bkill")) {XHOST.is_command("scancel")=false;XHOST.is_command("bkill")=false;} // force
    // if(XHOST.is_command("qdel")) {XHOST.is_command("scancel")=false;XHOST.is_command("bkill")=false;} // some priority

    stringstream aus;
    for (size_t i = 0; i < jobs.size(); i++) {
      if (XHOST.is_command("scancel")) {
        aus.clear();
        aus.str(std::string());
        aus << XHOST.command("scancel") << " " << jobs[i] << endl;  // cmd = qdel OR bkill
        vcmds.push_back(aus.str());
      } else {
        if (XHOST.is_command("qdel")) {
          aus.clear();
          aus.str(std::string());
          aus << XHOST.command("qdel") << " " << jobs[i] << endl;  // cmd = qdel OR bkill
          vcmds.push_back(aus.str());
        } else {
          if (XHOST.is_command("bkill")) {
            aus.clear();
            aus.str(std::string());
            aus << XHOST.command("bkill") << " " << jobs[i] << endl;  // cmd = qdel OR bkill
            vcmds.push_back(aus.str());
          }
        }
      }
    }
    for (size_t i = 0; i < vcmds.size(); i++) {
      cout << "EXECUTING: " << vcmds[i];// << endl;
      aurostd::execute(vcmds[i]);
    }
  }
} // namespace sflow

// ***************************************************************************
// sflow::QSUB qsub bsub sbatch
// ***************************************************************************
namespace sflow {
  void QSUB(string options, string cmd) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << XPID << "sflow::QSUB: BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    if (tokens.size() != 2) {
      init::ErrorOption(options, "sflow::QSUB", "aflow --qsub=N,file");
    }

    vector<string> vcmds;
    stringstream aus;
    for (uint i = 0; i < (uint) aurostd::string2utype<int>(tokens.at(0)); i++) {
      aus.clear();
      aus.str(std::string());
      aus << cmd << " " << tokens.at(1) << endl;  // cmd = sbatch OR bsub <
      vcmds.push_back(aus.str());
    }
    for (size_t i = 0; i < vcmds.size(); i++) {
      cout << "EXECUTING: " << vcmds[i];// << endl;
      aurostd::execute(vcmds[i]);
    }
    if (LDEBUG) {
      cerr << XPID << "sflow::QSUB: END" << endl;
    }
  }
} // namespace sflow

namespace sflow {
  void QSUB(string options) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << XPID << "sflow::QSUB: BEGIN" << endl;
    }
    vector<string> tokens;
    aurostd::string2tokens(options, tokens, ",");

    if (tokens.size() != 2) {
      init::ErrorOption(options, "sflow::QSUB", "aflow --qsub=N,file");
    }

    vector<string> vcmds;
    stringstream aus;
    for (uint i = 0; i < (uint) aurostd::string2utype<int>(tokens.at(0)); i++) {
      if (XHOST.is_command("sbatch")) {
        aus.clear();
        aus.str(std::string());
        if (XHOST.is_MACHINE_FULTON_MARYLOU) {
          aus << XHOST.command("sbatch") << " -C beta " << tokens.at(1) << endl;  // cmd = sbatch
        } else {
          aus << XHOST.command("sbatch") << "  " << tokens.at(1) << endl;  // cmd = sbatch
        }
        vcmds.push_back(aus.str());
      } else {
        if (XHOST.is_command("bsub")) {
          aus.clear();
          aus.str(std::string());
          aus << XHOST.command("bsub") << " <" << " " << tokens.at(1) << endl;  // cmd = bsub <
          vcmds.push_back(aus.str());
        } else {
          if (XHOST.is_command("qsub")) {
            aus.clear();
            aus.str(std::string());
            aus << XHOST.command("qsub") << " " << tokens.at(1) << endl;  // cmd = qsub
            vcmds.push_back(aus.str());
          }
        }
      }
    }
    for (size_t i = 0; i < vcmds.size(); i++) {
      cout << "EXECUTING: " << vcmds[i];// << endl;
      aurostd::execute(vcmds[i]);
    }
    if (LDEBUG) {
      cerr << XPID << "sflow::QSUB: END" << endl;
    }
  }
} // namespace sflow

// **************************************************************************

#endif  // _PTHREADS_IMPLEMENTATIONS_

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
