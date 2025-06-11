
#ifndef AFLOW_PTHREADS_H
#define AFLOW_PTHREADS_H

#include <deque>
#include <istream>
#include <string>
#include <vector>

#include "aflow_defs.h"

namespace AFLOW_PTHREADS {
  extern bool FLAG;         // run pthread YES/NO
  extern int MAX_PTHREADS;  // how many MAX threads I can use  default or --np
  extern int RUNNING;       // how many threads are actually running
  extern pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  extern int viret[MAX_ALLOCATABLE_PTHREADS];          // the thread runnings
  extern bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
} // namespace AFLOW_PTHREADS

namespace AFLOW_PTHREADS {
  int GetTotalCPUs();
  bool Check_Threads(std::vector<std::string> argv, const bool& VERBOSE);
  void Clean_Threads();
  void No_Threads();
  bool Available_Free_Threads(int& fthread);
  bool Wait_Available_Free_Threads(int& fthread, const double& pthread_wait, const bool& VERBOSE);
  bool Wait_Available_Free_Threads(int& fthread, const bool& VERBOSE);

  bool MULTI_sh(const std::vector<std::string>& argv);
} // namespace AFLOW_PTHREADS

namespace aurostd { // Multithreaded add on to aurostd
  bool multithread_execute(const std::deque<std::string>& vcommand, int NUM_THREADS, bool VERBOSE);
  bool multithread_execute(const std::deque<std::string>& vcommand, int NUM_THREADS);
  bool multithread_execute(const std::deque<std::string>& vcommand);
  bool multithread_execute(const std::vector<std::string>& vcommand, int NUM_THREADS, bool VERBOSE);
  bool multithread_execute(const std::vector<std::string>& vcommand, int NUM_THREADS);
  bool multithread_execute(const std::vector<std::string>& vcommand);
} // namespace aurostd

namespace sflow {
  void KILL(std::string options);
  void JUST(std::string options, std::istream& input, std::string mode);
  void QSUB(std::string options);
  void QSUB(std::string options, std::string cmd);
  void QDEL(std::string options);
  void QDEL(std::string options, std::string cmd);
} // namespace sflow

#endif // AFLOW_PTHREADS_H
