#ifndef AFLOW_SQLITE_H
#define AFLOW_SQLITE_H

#include <string>
#include <vector>

#include "extern/SQLITE/sqlite3.h"

namespace sql {
  void SQLexecuteCommand(sqlite3*, const std::string&);
  std::string SQLexecuteCommandSCALAR(sqlite3*, const std::string&);
  std::vector<std::string> SQLexecuteCommandVECTOR(sqlite3*, const std::string&);
  std::vector<std::vector<std::string>> SQLexecuteCommand2DVECTOR(sqlite3*, const std::string&);
  int SQLcallback(void*, int, char**, char**);
  int SQLcallbackSCALAR(void*, int, char**, char**);
  int SQLcallbackVECTOR(void*, int, char**, char**);
  int SQLcallback2DVECTOR(void*, int, char**, char**);
} // namespace sql

#endif
