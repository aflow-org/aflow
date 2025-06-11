//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
//****************************************************************************

// These functions provide a direct interface to SQLite. They execute commands
// and handle callbacks.

#include "aflow_sqlite.h"

#define _AFLOW_SQL_DEBUG_    false
#define _SQL_COMMAND_DEBUG_  false  // debug SQL commands that are sent - verbose output
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
    if (_SQL_COMMAND_DEBUG_) std::cerr << XPID << "sql::SQLexecuteCommand(): command = " << command << std::endl;
    char* sqlErrMsg = 0;
    int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback, 0, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      std::string function = "sql::SQLexecuteCommand():";
      std::string message = string(sqlErrMsg) + " in command " + command;
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    }
  }

  string SQLexecuteCommandSCALAR(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) std::cerr << XPID << "sql::SQLexecuteCommandSCALAR(): command = " << command << std::endl;
    char* sqlErrMsg = 0;
    string returnstring = "";
    int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackSCALAR, &returnstring, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      string function = XPID + "sql::SQLexecuteCommandSCALAR():";
      string message = string(sqlErrMsg) + " in command " + command;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    } else {
      return returnstring;
    }
  }

  vector<string> SQLexecuteCommandVECTOR(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) std::cerr << XPID << "sql::SQLexecuteCommandVECTOR(): command = " << command << std::endl;
    char *sqlErrMsg = 0;
    vector<string> returnvector;
    int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackVECTOR, &returnvector, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      string function = XPID + "sql::SQLexecuteCommandVECTOR():";
      string message = string(sqlErrMsg) + " in command " + command;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, function, message, _RUNTIME_SQL_);
    } else {
      return returnvector;
    }
  }

  vector<vector<string> > SQLexecuteCommand2DVECTOR(sqlite3* cursor, const string& command) {
    if (_SQL_COMMAND_DEBUG_) std::cerr << XPID << "sql::SQLexecuteCommand2DVECTOR(): command = " << command << std::endl;
    char *sqlErrMsg = 0;
    vector<vector<string> > returnvector;
    int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback2DVECTOR, &returnvector, &sqlErrMsg);
    if (sql_code != SQLITE_OK) {
      string function = XPID + "sql::SQLexecuteCommand2DVECTOR():";
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
        std::cerr << XPID << "sql::SQLcallback()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    return 0;
  }

  int SQLcallbackSCALAR(void* data, int argc, char** argv, char** col) {
    if (_SQL_CALLBACK_DEBUG_) {
      for (int i = 0; i < argc; i++) {
        std::cerr << XPID << "sql::SQLcallbackSCALAR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    string* val = static_cast<string*>(data);
    if (argc == 1) {
      if (argv[0] == NULL) {
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
        std::cerr << XPID << "sql::SQLcallbackVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    vector<string>* vec = static_cast<vector<string>*>(data);
    for (int i = 0; i < argc; i++) {
      if (argv[i] != NULL) vec->push_back(string(argv[i]));
      else vec->push_back("");
    }
    return 0;
  }

  int SQLcallback2DVECTOR(void* data, int argc, char** argv, char** col) {
    if (_SQL_CALLBACK_DEBUG_) {
      for (int i = 0; i < argc; i++) {
        std::cerr << XPID << "sql::SQLcallback2DVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
      }
      std::cerr << std::endl;
    }
    vector<vector<string> >* vec2d = static_cast<vector<vector<string> >*>(data);
    vector<string> vec(argc, "");
    for (int i = 0; i < argc; i++) {
      if (argv[i] != NULL) vec[i] = string(argv[i]);
    }
    vec2d->push_back(vec);
    return 0;
  }

}  // namespace sql

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
//****************************************************************************
