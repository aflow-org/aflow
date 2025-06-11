//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
//****************************************************************************

// These functions provide a direct interface to SQLite. They execute commands
// and handle callbacks.

#include "interfaces/aflow_sqlite.h"

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"

#include "aflow_xhost.h"
#include "extern/SQLITE/sqlite3.h"

#define _AFLOW_SQL_DEBUG_ false
#define _SQL_COMMAND_DEBUG_ false  // debug SQL commands that are sent - verbose output
#define _SQL_CALLBACK_DEBUG_ false  // debug SQL callbacks - extremely verbose output

using std::string;
using std::vector;

namespace sql {

  // Execute command functions -----------------------------------------------
  // These functions provide a framework to execute SQLite commands (including
  // exception handling). The return types should all be void or string-typed
  // to keep the number of functions to a minimum. Conversion to other data
  // types should be handled outside these functions.

  void SQLexecuteCommand(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) {
      std::cerr << XPID << "sql::SQLexecuteCommand(): command = " << command << std::endl;
    }
    char* sqlErrMsg = nullptr;
    const int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback, nullptr, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      const std::string function = "sql::SQLexecuteCommand():";
      const std::string message = string(sqlErrMsg) + " in command " + command;
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    }
  }

  string SQLexecuteCommandSCALAR(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) {
      std::cerr << XPID << "sql::SQLexecuteCommandSCALAR(): command = " << command << std::endl;
    }
    char* sqlErrMsg = nullptr;
    string returnstring;
    const int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackSCALAR, &returnstring, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      const string function = XPID + "sql::SQLexecuteCommandSCALAR():";
      string message = string(sqlErrMsg) + " in command " + command;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    } else {
      return returnstring;
    }
  }

  vector<string> SQLexecuteCommandVECTOR(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) {
      std::cerr << XPID << "sql::SQLexecuteCommandVECTOR(): command = " << command << std::endl;
    }
    char* sqlErrMsg = nullptr;
    vector<string> returnvector;
    const int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackVECTOR, &returnvector, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      const string function = XPID + "sql::SQLexecuteCommandVECTOR():";
      string message = string(sqlErrMsg) + " in command " + command;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    } else {
      return returnvector;
    }
  }

  vector<vector<string>> SQLexecuteCommand2DVECTOR(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) {
      std::cerr << XPID << "sql::SQLexecuteCommand2DVECTOR(): command = " << command << std::endl;
    }
    char* sqlErrMsg = nullptr;
    vector<vector<string>> returnvector;
    const int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback2DVECTOR, &returnvector, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      const string function = XPID + "sql::SQLexecuteCommand2DVECTOR():";
      string message = string(sqlErrMsg) + " in command " + command;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    } else {
      return returnvector;
    }
  }

  // Callback functions ------------------------------------------------------
  // The following functions are the callback functions passed into
  // sqlite_exec. Each executeCommand function should have its own callback
  // function.

  int SQLcallback(void* data, int argc, char** argv, char** col) {
    (void) data;  // To suppress compiler warnings
    if (_SQL_CALLBACK_DEBUG_) {
      for (int i = 0; i < argc; i++) {
        std::cerr << XPID << "sql::SQLcallback()[" << i << "]: " << col[i] << " = " << (argv[i] ? argv[i] : "NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    return 0;
  }

  int SQLcallbackSCALAR(void* data, int argc, char** argv, char** col) {
    if (_SQL_CALLBACK_DEBUG_) {
      for (int i = 0; i < argc; i++) {
        std::cerr << XPID << "sql::SQLcallbackSCALAR()[" << i << "]: " << col[i] << " = " << (argv[i] ? argv[i] : "NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    string* val = static_cast<string*>(data);
    if (argc == 1) {
      if (argv[0] == nullptr) {
        *val = "";
      } else {
        *val = std::string(argv[0]);
      }
      return 0;
    } else {
      return 1;
    }
  }

  int SQLcallbackVECTOR(void* data, int argc, char** argv, char** col) {
    if (_SQL_CALLBACK_DEBUG_) {
      for (int i = 0; i < argc; i++) {
        std::cerr << XPID << "sql::SQLcallbackVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i] ? argv[i] : "NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    vector<string>* vec = static_cast<vector<string>*>(data);
    for (int i = 0; i < argc; i++) {
      if (argv[i] != nullptr) {
        vec->emplace_back(argv[i]);
      } else {
        vec->emplace_back("");
      }
    }
    return 0;
  }

  int SQLcallback2DVECTOR(void* data, int argc, char** argv, char** col) {
    if (_SQL_CALLBACK_DEBUG_) {
      for (int i = 0; i < argc; i++) {
        std::cerr << XPID << "sql::SQLcallback2DVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i] ? argv[i] : "NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    vector<vector<string>>* vec2d = static_cast<vector<vector<string>>*>(data);
    vector<string> vec(argc, "");
    for (int i = 0; i < argc; i++) {
      if (argv[i] != nullptr) {
        vec[i] = string(argv[i]);
      }
    }
    vec2d->push_back(vec);
    return 0;
  }

}  // namespace sql

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
//****************************************************************************
