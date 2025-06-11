
#ifndef AFLOW_INIT_H
#define AFLOW_INIT_H

#include <iostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"

#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "flow/aflow_xclasses.h"

// --------------------------------------------------------------------------
// aflow_init.cpp
namespace init {
  int GetCPUCores();
  int InitMachine(bool INIT_VERBOSE, std::vector<std::string>& argv, std::vector<std::string>& cmds, std::ostream& outf);  // ME20200724 - changed to int
  std::string AFLOW_Projects_Directories(std::string string2load);
  long GetRAM();
  uint GetTEMP();
  double WaitTEMP(double TRESHOLD = AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT, std::ostream& oss = std::cout, bool LVERBOSE = false, std::vector<std::string> vmessage = std::vector<std::string>(0));
  uint InitSchema(bool INIT_VERBOSE);
  uint InitSchemaInternal(bool INIT_VERBOSE);  // ME20220208
  std::vector<std::string> getSchemaKeys(const aurostd::xoption& vschema);  // ME20220223
  std::vector<std::string> getSchemaNames(const aurostd::xoption& vschema);  // CO20200520
  std::vector<std::string> getSchemaTypes(const aurostd::xoption& vschema);  // ME20220223
  std::vector<std::string> getSchemaTypes(const aurostd::xoption& vschema, const std::vector<std::string>& keys);  // ME20220223
} // namespace init

uint AFLOW_getTEMP(const std::vector<std::string>& argv);
uint AFLOW_monitor(const std::vector<std::string>& argv);
double AFLOW_checkMEMORY(const std::string& progname = "", double = 102.0);
bool CheckMaterialServer(const std::string& message);  // CO20200624
bool CheckMaterialServer();

class _xvasp; // forward declaration

bool GetVASPBinaryFromLOCK(const std::string& directory, std::string& vasp_bin);  // CO20210315
bool GetVASPBinaryFromLOCK(const std::string& directory, std::string& vasp_bin, int& ncpus);  // CO20210315
void processFlagsFromLOCK(_xvasp& xvasp, _vflags& vflags, aurostd::xoption& xfixed);  // CO20210315
bool AFLOW_VASP_instance_running(); // CO20210315
bool AFLOW_VASP_instance_running(const std::string& pgid); // SD20220330
bool AFLOW_MONITOR_instance_running(const _aflags& aflags); // CO20210315
bool VASP_instance_running(const std::string& vasp_bin); // CO20210315
bool VASP_instance_running(const std::string& vasp_bin, const std::string& pgid); // SD20220330
void AFLOW_monitor_VASP();  // CO20210315
void AFLOW_monitor_VASP(const std::string& directory);  // CO20210315

std::string Message(const std::string& filename, const std::string& list2print = _AFLOW_MESSAGE_DEFAULTS_);  // CO20200713
std::string Message(const std::string& filename, const _aflags& aflags, const std::string& list2print = _AFLOW_MESSAGE_DEFAULTS_);  // CO20200713
bool AFLOW_BlackList(const std::string& h);  // CO20200713
namespace init {
  void MessageOption(const std::string& options, const std::string& routine, std::vector<std::string> vusage);  // CO20200624 - should go to cerr for web //DX20200724 - bool to void
  void MessageOption(const std::string& options, const std::string& routine, std::string vusage);  // CO20200624 - should go to cerr for web //DX20200724 - bool to void
  void ErrorOption(const std::string& options, const std::string& routine, std::vector<std::string> vusage);  // CO20200624 - should go to cerr for web //DX20200724 - bool to void
  void ErrorOption(const std::string& options, const std::string& routine, std::string vusage);  // CO20200624 - should go to cerr for web //DX20200724 - bool to void
} // namespace init

#endif // AFLOW_INIT_H
