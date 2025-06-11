#include "../AUROSTD/aurostd.h"
#include "sqlite3.h"

namespace sql {
  void SQLexecuteCommand(sqlite3*, const string&);
  string SQLexecuteCommandSCALAR(sqlite3*, const string&);
  vector<string> SQLexecuteCommandVECTOR(sqlite3*, const string&);
  vector<vector<string> > SQLexecuteCommand2DVECTOR(sqlite3*, const string&);
  int SQLcallback(void*, int, char**, char**);
  int SQLcallbackSCALAR(void*, int, char**, char**);
  int SQLcallbackVECTOR(void*, int, char**, char**);
  int SQLcallback2DVECTOR(void*, int, char**, char**);
}
