//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
//****************************************************************************

// Class for the AFLOW database, based on SQLite. It provides the framework
// for rebuilding and interacting with the database file using data in JSON
// format. This file is divided into the following sections:
//
// * Constructor
// * Temp file handling
// * Database builder
// * Database analyzer
// * User-level functions to interact with the database
//
// ---------------------------------- USAGE ----------------------------------
//
// To update the database, i.e. to rebuild the database only when new entries
// are available, use: aflow --update_database
//
// To force a database rebuild, use: aflow --rebuild_database
//
// Database statistics will be collected in each case. To only perform the
// analysis, use: aflow --analyze_database
//
// ----------------------------- REBUILD PROCESS -----------------------------
//
// The database rebuild process has the following steps:
//   * Check if the database needs to be rebuild. This can be triggered either
//     by the user, by new entries in the schema, or by new/updated dat files.
//     The dat files are collections of aflowlib.json files for a set of
//     AUIDs (i.e. aflow:00.jsonl, aflow:01.jsonl, etc.).
//   * Create a temporary database file and populate with data from these JSONs.
//     Rebuilding from scratch instead of incrementally adding into the existing
//     database protects the database from corruption and injection attacks.
//     Using a temporary database file prevents that the active database file
//     get corrupted due to an error in the build process.
//   * Compare temporary and current database file and copy over if necessary.
//     The database file only gets copied after successful completion of the
//     rebuild process. A rebuild is considered successful when the number of
//     entries and properties of the new database file is greater than or equal
//     to that of the old file.
//   * Analyze the database and output the database statistics.
//
// ------------------------------ RETURN CODES -------------------------------
//
// The database rebuild uses return codes to convey to outside processes why
// it terminated. The codes are listed below.

// Return code definitions

// Success
#define _AFLOW_DB_SUCCESS_            0  // Database successfully updated

// Codes that do not trigger a database analysis
#define _AFLOW_DB_NOT_UPDATED_       100  // Nothing to update or patch
#define _AFLOW_DB_NOT_COPIED_        101  // Temp database file was not copied
#define _AFLOW_DB_NOT_PATCHED_       102  // No patches applied even though files were given

// Codes that trigger a database analysis
#define _AFLOW_DB_PATCH_SKIPPED_     200  // Returned when rebuild was successful but no patch was applied
#define _AFLOW_DB_PATCH_INCOMPLETE_  201  // At least one patch was incomplete

#include "aflowlib.h"

#define _AFLOW_DB_DEBUG_ false

using std::string;
using std::vector;
using namespace std::placeholders;

static const uint _DEFAULT_SET_LIMIT_ = 16; //DX20200317 - int -> uint
static const int _N_AUID_TABLES_ = 256;

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace aflowlib {

  // Open the database for read access
  AflowDB::AflowDB(const string& db_file, ostream& oss) : xStream(oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " reading database" << std::endl;
    free();
    aurostd::xoption dummy;
    initialize(db_file, "", "", SQLITE_OPEN_READONLY, dummy, dummy);
  }

  AflowDB::AflowDB(const string& db_file, const aurostd::xoption& vschema_in, ostream& oss) : xStream(oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " reading database" << std::endl;
    free();
    aurostd::xoption dummy;
    initialize(db_file, "", "", SQLITE_OPEN_READONLY, vschema_in, dummy);
  }

  AflowDB::AflowDB(const string& db_file, const aurostd::xoption& vschema_in, const aurostd::xoption& vschema_internal_in, ostream& oss) : xStream(oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " reading database" << std::endl;
    free();
    initialize(db_file, "", "", SQLITE_OPEN_READONLY, vschema_in, vschema_internal_in);
  }

  // Open the database for write access
  AflowDB::AflowDB(const string& db_file, const string& dt_path, const string& lck_file,
      const aurostd::xoption& vschema_in, ostream& oss) : xStream (oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    free();
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Database file: " << database_file << std::endl;
      std::cerr << __AFLOW_FUNC__ << " Data path: " << data_path << std::endl;
      std::cerr << __AFLOW_FUNC__ << " Lock file: " << lock_file << std::endl;
    }
    aurostd::xoption dummy;
    initialize(db_file, dt_path, lck_file, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, vschema_in, dummy);
  }

  AflowDB::AflowDB(const string& db_file, const string& dt_path, const string& lck_file,
      const aurostd::xoption& vschema_in, const aurostd::xoption& vschema_internal_in, ostream& oss) : xStream (oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    free();
    if (LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " Database file: " << database_file << std::endl;
      std::cerr << __AFLOW_FUNC__ << " Data path: " << data_path << std::endl;
      std::cerr << __AFLOW_FUNC__ << " Lock file: " << lock_file << std::endl;
    }
    initialize(db_file, dt_path, lck_file, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, vschema_in, vschema_internal_in);
  }

  void AflowDB::initialize(const string& db_file, const string& dt_path, const string& lck_file,
    int open_flags, const aurostd::xoption& schema_in, const aurostd::xoption& vschema_internal_in) {
    data_path = dt_path;
    database_file = db_file;
    lock_file = lck_file;
    vschema = schema_in;
    vschema_internal = vschema_internal_in;
    open(open_flags);
  }

  // Copy constructors
  AflowDB::AflowDB(const AflowDB& that) : xStream(*that.getOFStream(), *that.getOSS()) {
    copy(that);
  }

  AflowDB& AflowDB::operator=(const AflowDB& that) {
    copy(that);
    return *this;
  }

  void AflowDB::copy(const AflowDB& that) {
    if (this == &that) return;
    xStream::copy(that);
    data_path = that.data_path;
    database_file = that.database_file;
    lock_file = that.lock_file;
    vschema_internal = that.vschema_internal;
    // Two databases should never operate on a temporary
    // file or have write access at the same time.
    is_tmp = false;
    open(SQLITE_OPEN_READONLY);
  }

  // Destructor
  AflowDB::~AflowDB() {
    close();
    free();
  }

  void AflowDB::free() {
    data_path = "";
    database_file = "";
    lock_file = "";
    vschema.clear();
    vschema_internal.clear();
    is_tmp = false;
  }

  void AflowDB::clear() {
    close();
    free();
  }

  // Opens the database file and creates the main cursor
  void AflowDB::open(int open_flags) {
    string message = "Opening " + database_file + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    int sql_code = sqlite3_open_v2(database_file.c_str(), &db, open_flags, NULL); //DX20200319 - nullptr -> NULL
    if (sql_code != SQLITE_OK) {
      message = "Could not open database file " + database_file;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  // Closes the database file and removes the main cursor
  void AflowDB::close() {
    if (isTMP()) closeTmpFile(false, true, true);
    string message = "Closing " + database_file + ".";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    int sql_code = sqlite3_close(db);
    if (sql_code != SQLITE_OK) {
      message = "Could not close database file " + database_file;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

}  // namespace aflowlib

/********************************* TMP FILE **********************************/

namespace aflowlib {

  //openTmpFile///////////////////////////////////////////////////////////////
  // Opens a temporary database file 
  void AflowDB::openTmpFile(int open_flags, bool copy_original) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    stringstream message;

    string tmp_file = database_file + ".tmp";
    message << "Opening " << tmp_file;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);

    // Create a symbolic link to lock the tmp file creation process. Since creating
    // symbolic links is an atomic process, it can be used to prevent race conditions.
    string lock_link = lock_file + ".lnk";
    if (symlink(lock_file.c_str(), lock_link.c_str())) {
      message << "Could not create symbolic link to lock file " + lock_file + ": ";
      if (errno == EEXIST) {
        message << "Another process has created it already.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_BUSY_);
      } else {
        message << "Process exited with errno " << errno << ".";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }

    // If there is a temporary database file, a rebuild process is either in
    // progress or failed.
    if (aurostd::FileExist(tmp_file)) {
      long int tm_tmp = aurostd::GetTimestampModified(tmp_file);
      time_t t = std::time(NULL); //DX20200319 - nullptr -> NULL
      long int tm_curr = (long int) t;
      int pid = -1;
      bool del_tmp = false;
      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " Found temporary database file with time stamp " << tm_tmp
          << " (current time: " << tm_curr << ")." << std::endl;
      }
      // Check if the rebuild process that is saved in the lock file
      // is still active.
      if (aurostd::FileExist(lock_file) && !aurostd::FileEmpty(lock_file)) {
        pid = aurostd::string2utype<int>(aurostd::file2string(lock_file));
        if (kill(pid, 0)) pid = -1;
      }
      if (tm_curr - tm_tmp < DEFAULT_AFLOW_DB_STALE_THRESHOLD) {
        if (pid > -1) del_tmp = false;
        else del_tmp = true;
      } else {
        del_tmp = true;
      }

      if (del_tmp) {
        if (LDEBUG) {
          message << "Temporary database file already exists, but process has become stale."
            << " Removing temporary files and killing outstanding process.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
        }
        if (pid > -1) {
          int killreturn = kill(pid, 9);
          if (killreturn) {
            aurostd::RemoveFile(lock_link);
            message << "Could not kill active process " << pid << "(errno =  " << errno << ").";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }
        }
        aurostd::RemoveFile(tmp_file);
        aurostd::RemoveFile(tmp_file + "-journal");
      } else {
        aurostd::RemoveFile(lock_link);
        message << "Could not create temporary database file. File already exists and rebuild process is active.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_BUSY_);
      }
    }

    aurostd::string2file(aurostd::utype2string<int>(getpid()), lock_file);
    int sql_code = sqlite3_close(db);
    if (sql_code != SQLITE_OK) {
      aurostd::RemoveFile(lock_link);
      message << "Could not close main database file " + database_file << " (SQL code " << sql_code << ").";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    if (copy_original) {
      bool copied = aurostd::CopyFile(database_file, tmp_file);
      // Make sure the file can be modified
      if (copied) {
        aurostd::ChmodFile("664", tmp_file);
      } else {
        aurostd::RemoveFile(lock_link);
        message << "Unable to copy original database file.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
    }

    sql_code = sqlite3_open_v2(tmp_file.c_str(), &db, open_flags, NULL); //DX20200319 - nullptr -> NULL
    if (sql_code != SQLITE_OK) {
      message << "Could not open tmp database file " + tmp_file << " (SQL code " << sql_code << ").";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    is_tmp = true;
    aurostd::RemoveFile(lock_link);
  }

  //closeTmpFile//////////////////////////////////////////////////////////////
  // Closes a temporary database file and overwrites the original if the
  // build was successful. Unless --rebuild_database was selected by the user,
  // this function tests whether the new database has less entries than the
  // old one. To avoid data loss, the original file should not be overwritten.
  // nocopy: just close the tmp file without copying
  bool AflowDB::closeTmpFile(bool force_copy, bool keep, bool nocopy) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    string tmp_file = database_file + ".tmp";
    stringstream message;
    message << "Closing " << tmp_file;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);

    // To throw
    bool found_error = false;
    stringstream message_error;

    // If the tmp file is empty, the build failed, so all other tests will fail as well.
    if (aurostd::FileEmpty(tmp_file)) {
      message_error << "Temporary database file empty. Database will not be overwritten.";
      nocopy = true;
      found_error = true;
    }
    if (!found_error && !nocopy) {
      if ((int) getTables().size() != _N_AUID_TABLES_) {
        message_error << "Temporary database file has the wrong number of tables."
          << " Database file will not be copied.";
        nocopy = true;
        found_error = true;
      }
    }

    // If the user does not force a rebuild, test that the new database does not
    // result in less data than the old one. To build a database with less
    // properties, aflow needs to be run with the --rebuild_database option.
    uint nentries_tmp = 0, ncols_tmp = 0, nentries_old = 0, ncols_old = 0;
    if (!found_error && !force_copy && !nocopy) {
      vector<string> props, tables;
      tables = getTables();
      // Get number of properties
      props = getColumnNames(tables[0]);  // All tables have the same columns, so any table is good
      ncols_tmp = props.size();

      // Get number of entries
      props = getDatabasePropertyMultiTables("COUNT", tables, "*");
      for (int i = 0; i < _N_AUID_TABLES_; i++) nentries_tmp += aurostd::string2utype<uint>(props[i]);
    }

    int sql_code = sqlite3_close(db);
    if (sql_code != SQLITE_OK) {
      message_error << "Could not close tmp database file " + tmp_file
        << " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message_error, _FILE_ERROR_);
    }
    is_tmp = false;

    bool copied = false;
    if (found_error || nocopy) {
      copied = false;
    } else if (force_copy) {
      message << "Force copy selected. Database will be overwritten.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
      copied = aurostd::CopyFile(tmp_file, database_file);
    } else {
      if (!aurostd::FileEmpty(database_file)) {
        message << "Old database file found. Determining number of entries and properties to prevent data loss.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);

        sqlite3* main_db;
        sql_code = sqlite3_open_v2(database_file.c_str(), &main_db, SQLITE_OPEN_READONLY, NULL);
        if (sql_code != SQLITE_OK) {
          message_error << "Could not close main database file " + database_file
            << " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message_error, _FILE_ERROR_);
        }
        vector<string> props, tables;
        tables = getTables(main_db);

        // Get number of properties
        props = getColumnNames(main_db, tables[0]);
        ncols_old = props.size();

        // Get number of entries
        props = getDatabasePropertyMultiTables(main_db, "COUNT", tables, "*");
        for (int i = 0; i < _N_AUID_TABLES_; i++) nentries_old += aurostd::string2utype<uint>(props[i]);

        sql_code = sqlite3_close(main_db);
        if (sql_code != SQLITE_OK) {
          message_error << "Could not close main database file " + database_file
            << " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message_error, _FILE_ERROR_);
        }
      } else {
        if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " No old datbase found." << std::endl;
        ncols_old = 0;
        nentries_old = 0;
      }

      if (LDEBUG) {
        std::cerr << __AFLOW_FUNC__ << " Number of entries: " << nentries_tmp << " (current database: " << nentries_old << ")." << std::endl;
        std::cerr << __AFLOW_FUNC__ << " Number of properties: " << ncols_tmp << " (current database: " << ncols_old << ")." << std::endl;
      }

      if (nentries_tmp < nentries_old) {
        message_error << "The rebuild process resulted in less entries"
          << " than in the current database. To prevent accidental"
          << " data loss, the temporary database will not be copied."
          << " Rerun as aflow --rebuild_database if the database should"
          << " be updated regardless. The temporary database file will be"
          << " kept to allow debugging.";
        keep = true;
        copied = false;
        found_error = true;
      } else if (ncols_tmp < ncols_old) {
        message_error << "The rebuild process resulted in less properties"
          << " than in the current database. To prevent accidental"
          << " data loss, the temporary database will not be copied."
          << " Rerun as aflow --rebuild_database if the database should"
          << " be updated regardless. The temporary database file will be"
          << " kept to allow debugging.";
        keep = true;
        copied = false;
        found_error = true;
      } else if (!found_error) {
        // No problems found, so replace old database file with new file
        copied = aurostd::CopyFile(tmp_file, database_file);
      }
      message << "Database file " << (copied?"successfully":"not") << " copied.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    }
    if (!found_error && !keep) {
      aurostd::RemoveFile(tmp_file);
      message << "Temporary database file removed.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    } else if (LDEBUG) {
      message << "Temporary database file " << tmp_file << " will not be deleted.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    }

    if (found_error) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    open();

    aurostd::string2file("", lock_file);
    if (copied) aurostd::ChmodFile("444", database_file);
    return copied;
  }

  //isTMP/////////////////////////////////////////////////////////////////////
  // Returns whether or not the cursor points to a temporary file or not
  bool AflowDB::isTMP() {
    return is_tmp;
  }

}  // namespace aflowlib

/***************************** REBUILD DATABASE *****************************/

namespace aflowlib {

  //rebuildDatabase///////////////////////////////////////////////////////////
  // This function first checks if a rebuild is necessary and then initiates
  // the rebuilding functions. Returns true if the database has been
  // successfully rebuilt.
  int AflowDB::rebuildDatabase(bool force_rebuild) {
    vector<string> patch_placeholder;
    return rebuildDatabase(patch_placeholder, force_rebuild);
  }

  int AflowDB::rebuildDatabase(const string& patch_files, bool force_rebuild) {
    vector<string> patch_files_input;
    if (!patch_files.empty()) aurostd::string2tokens(patch_files, patch_files_input, ",");
    return rebuildDatabase(patch_files_input, force_rebuild);
  }

  int AflowDB::rebuildDatabase(const vector<string>& patch_files_input, bool force_rebuild) {
    stringstream message;
    bool rebuild_db = false;

    // Always rebuild when the user wants the rebuild.
    rebuild_db = force_rebuild;
    if (rebuild_db) {
      message << "Rebuilding database (user rebuild).";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    }

    // Always rebuild if the database file doesn't exist. Note: When the
    // database is opened and the database file doesn't exist, SQLite creates
    // an empty file, so check if the file is empty and not if it exists.
    if (!rebuild_db) {
      rebuild_db = aurostd::FileEmpty(database_file); // file_empty needed for LDEBUG
      if (rebuild_db) {
        message << "Rebuilding database (file not found or empty).";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
      }
    }

    // Rebuild when the schema has been updated
    if (!rebuild_db) {
      vector<string> columns;
      string table = getTables()[0];  // All tables have the same columns, so any table is good
      columns = getColumnNames(table);

      // Properties
      vector<string> keys_schema = getAllSchemaKeys();
      uint nkeys = keys_schema.size();
      if (nkeys > columns.size()) {
        rebuild_db = true;
      } else if (nkeys < columns.size()) {
        message << "The number of properties (" << nkeys << ") is smaller than in the"
          << " current database (" << columns.size() << "). To prevent accidental"
          << " data loss, AFLOW will abort the update process. Rerun as"
          << " aflow --rebuild_database if the database should be updated regardless.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }

      // Data types
      if (!rebuild_db) {
        vector<string> types_db = getColumnTypes(table);
        vector<string> types_schema = getDataTypes(keys_schema, false);
        int index = -1;
        string key = "";
        uint k = 0;
        for (k = 0; k < nkeys; k++) {
          key = vschema.getattachedscheme("SCHEMA::NAME:" + keys_schema[k]);
          types_schema[k] = aurostd::RemoveSubString(types_schema[k], " COLLATE NOCASE");  // TEXT may contain directive to be case insensitive
          if (key.empty()) key = vschema_internal.getattachedscheme("SCHEMA::NAME:" + keys_schema[k]);
          if (!aurostd::WithinList(columns, key, index)) break;
          if (types_db[index] != types_schema[k]) break;
        }
        rebuild_db = (k != nkeys);
      }
      if (rebuild_db) {
        message << " Rebuilding database (schema updated).";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
      }
    }

    // Check if any relevant files are newer than the database.
    if (!rebuild_db) {
      vector<string> json_files(_N_AUID_TABLES_);
      for (int i = 0; i < _N_AUID_TABLES_; i++) {
        stringstream t;
        t << std::setfill('0') << std::setw(2) << std::hex << i;
        json_files[i] = aurostd::CleanFileName(data_path + "/aflow:" + t.str() + ".jsonl");
        if (!aurostd::EFileExist(json_files[i]) && !aurostd::FileExist(json_files[i])) {
          message << data_path << " is not a valid data path. Missing file for aflow:" << t.str() << ".";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
        }
      }

      long int tm_db = aurostd::GetTimestampModified(database_file);
      int i = 0;
      for (i = 0; i < _N_AUID_TABLES_; i++) {
        if (aurostd::GetTimestampModified(json_files[i]) > tm_db) break;
      }
      rebuild_db = (i != _N_AUID_TABLES_);

      if (rebuild_db) message << "Rebuilding database (updated files found).";
      else message << "No updated files found. Database will not be rebuilt.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    }

    // Check if the database may need to be patched
    uint npatch_input = patch_files_input.size();

    if (rebuild_db) {
      openTmpFile();
      rebuildDB();
      bool close_tmp = closeTmpFile(force_rebuild);
      // If the tmp database was not copied, something went wrong - return 1
      if (!close_tmp) return _AFLOW_DB_NOT_COPIED_;

      // No patches needed, so return
      if (npatch_input == 0) return (close_tmp?_AFLOW_DB_SUCCESS_:_AFLOW_DB_NOT_COPIED_);
    }

    if (npatch_input > 0) {
      // Do not check timestamps after a rebuild because patches must be applied
      int patch_code = patchDatabase(patch_files_input, !rebuild_db);
      if (rebuild_db) {
        if (patch_code == _AFLOW_DB_NOT_PATCHED_) return _AFLOW_DB_PATCH_SKIPPED_; // Database was rebuilt but not patched
        else return patch_code;
      } else if (patch_code == _AFLOW_DB_NOT_PATCHED_) {
        return _AFLOW_DB_NOT_UPDATED_;  // Database not rebuilt and not patched
      } else {
        return patch_code;
      }
    }

    return _AFLOW_DB_NOT_UPDATED_;
  }

  //patchDatabase/////////////////////////////////////////////////////////////
  // Patches the database using jsonl files. This can be used for emergency
  // patches before a lib2raw run or to add data that lib2raw does not cover.
  int AflowDB::patchDatabase(const string& patch_files, bool check_timestamps) {
    vector<string> patch_files_input;
    if (!patch_files.empty()) aurostd::string2tokens(patch_files, patch_files_input, ",");
    return patchDatabase(patch_files_input, check_timestamps);
  }

  int AflowDB::patchDatabase(const vector<string>& patch_files_input, bool check_timestamps) {
    stringstream message;
    int patch_code = 0;

    vector<string> patch_files;
    uint npatch_input = patch_files_input.size();
    if (npatch_input == 0) return 0;  // No files, so nothing to patch

    long int tm_db = 0;
    if (check_timestamps) tm_db = aurostd::GetTimestampModified(database_file);

    // Check files that need to be patched
    for (uint i = 0; i < npatch_input; i++) {
      string filename = patch_files_input[i];
      if (!aurostd::EFileExist(filename)) {
        // Check if the file is in the data path
        if (aurostd::EFileExist(data_path + "/" + filename)) {
          filename = aurostd::CleanFileName(data_path + "/" + filename);
        } else {
          message << "File " << filename << " not found. Skipping.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);

          // If a file is missing, then the patch is incomplete
          patch_code = _AFLOW_DB_PATCH_INCOMPLETE_;
          continue;
        }
      }

      if (!check_timestamps || (aurostd::GetTimestampModified(filename) > tm_db)) {
        patch_files.push_back(filename);
        message << "Adding file " << filename << " to patch list.";
      } else {
        message << "Skipping file " << filename << ". File is older than the database.";
      }
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
    }

    uint npatch = patch_files.size();
    if (npatch == 0) {  // Nothing to patch
      if (check_timestamps) return _AFLOW_DB_NOT_PATCHED_;  // No new files
      else return _AFLOW_DB_PATCH_INCOMPLETE_;  // All (new) files missing
    }

    // Start patching
    uint patches_applied = 0;
    uint entries_patched = 0;
    for (uint i = 0; i < npatch; i++) {
      message << "Patching database using " << patch_files[i] << ".";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
      openTmpFile(SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, true);  // Copy original

      // Do not constantly synchronize the database with file on disk (increases
      // performance significantly).
      sql::SQLexecuteCommand(db, "PRAGMA synchronous = OFF");

      vector<string> data;
      aurostd::efile2vectorstring(patch_files[i], data);
      entries_patched = applyPatchFromJsonl(data);
      if (entries_patched > 0) {
        if (closeTmpFile(false)) {  // Do not force copy or a bad patch may break the database
          message << "Patched " << entries_patched << "/" << data.size() << " entries.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
          // Patch incomplete, so adjust code (unless an error was found before)
          if ((patch_code < 200) && (entries_patched < data.size())) patch_code = _AFLOW_DB_PATCH_INCOMPLETE_;
          patches_applied++;
        } else {
          message << "Could not close temporary database file.";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        }
      } else {
        closeTmpFile(false, false, true);  // Do not copy tmp file
        message << "No patches applied.";
        // Patch incomplete, so adjust code (unless an error was found before)
        if ((patch_code < 200) && (entries_patched < data.size())) patch_code = _AFLOW_DB_PATCH_INCOMPLETE_;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
    }

    if (patches_applied == 0) patch_code = _AFLOW_DB_NOT_PATCHED_;

    return patch_code;
  }

  // Rebuild -----------------------------------------------------------------

  //rebuildDB/////////////////////////////////////////////////////////////////
  // Rebuilds the database from scratch.
  void AflowDB::rebuildDB() {
    stringstream message;

    message << "Starting rebuild.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);

    // Do not constantly synchronize the database with file on disk (increases
    // performance significantly).
    sql::SQLexecuteCommand(db, "PRAGMA synchronous = OFF");

    // Get columns and types from schema
    vector<string> keys = getAllSchemaKeys();
    uint nkeys = keys.size();
    vector<string> columns(nkeys);
    for (uint k = 0; k < nkeys; k++) {
      columns[k] = vschema.getattachedscheme("SCHEMA::NAME:" + keys[k]);
      if (columns[k].empty()) columns[k] = vschema_internal.getattachedscheme("SCHEMA::NAME:" + keys[k]);
    }
    vector<string> types = getDataTypes(columns, true);

    // Rebuild
#ifdef AFLOW_MULTITHREADS_ENABLE
    int ncpus = init::GetCPUCores();
    int max_cpus = 16;
    if (ncpus < 1) ncpus = 1;
    if (ncpus > max_cpus) ncpus = max_cpus;
    xthread::xThread xt(ncpus);
    std::function<void(int, const vector<string>&, const vector<string>&)> fn = std::bind(&AflowDB::buildTable, this, _1, _2, _3);
    xt.run(_N_AUID_TABLES_, fn, columns, types);
#else
    for (uint i = 0; i < _N_AUID_TABLES_; i++) buildTable(i, columns, types);
#endif
    message << __AFLOW_FUNC__ << " Finished rebuild.";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);
  }

  //buildTables///////////////////////////////////////////////////////////////
  // Reads a .jsonl files and processes the JSONs for the database writer.
  void AflowDB::buildTable(int i, const vector<string>& columns, const vector<string>& types) {
    stringstream t;
    t << std::setfill('0') << std::setw(2) << std::hex << i;
    string table = "auid_" + t.str();
    string auid_prefix = "\"aflow:" + t.str();

    string jsonfile = aurostd::CleanFileName(data_path + "/aflow:" + t.str() + ".jsonl");
    vector<string> data;
    aurostd::efile2vectorstring(jsonfile, data);
    vector<vector<string> > values;
    uint ndata = data.size();
    string aurl = "", auid = "";
    for (uint d = 0; d < ndata; d++) {
      // Filter auids that are in the wrong file
      auid = aurostd::extractJsonValueAflow(data[d], "auid");
      if (auid.rfind(auid_prefix, 0) != 0) continue;
      // Filter non-POCC ARUNs
      aurl = aurostd::extractJsonValueAflow(data[d], "aurl");
      if ((aurl.find("ARUN") != string::npos)
          && (aurl.find("ARUN.POCC") == string::npos)) {
        continue;
      }
      values.push_back(getDataValues(data[d], columns, types));
    }
    populateTable(table, columns, types, values);
  }

  //populateTable/////////////////////////////////////////////////////////////
  // Populates the database tables and creates indexes. This function uses a
  // mutex to make the writing thread-safe. While SQLITE does not allow
  // concurrent writing anyway, multiple threads may open a transaction, which
  // causes the build to fail (hence the mutex). The mutex is also the reason
  // why this function should never do any processing.
  void AflowDB::populateTable(const string& table, const vector<string>& columns, const vector<string>& types, const vector<vector<string> >& values) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::unique_lock<std::mutex> lk(write_mutex);
#endif
    createTable(table, columns, types);
    int chunk_size = 1000, count = 0;
    transaction(true);
    for (uint v = 0, nvals = values.size(); v < nvals; v++) {
      insertValues(table, columns, values[v]);
      if (++count % chunk_size == 0 || v == nvals - 1) {
        transaction(false);
        transaction(true);
        count = 0;
      }
    }

    // Create indexes on important database properties
    vector<string> index_cols;
    string indexes = "alloy,auid,aurl,catalog,compound,nspecies,Pearson_symbol_relax,spacegroup_relax";
    aurostd::string2tokens(indexes, index_cols, ",");
    string index, index_expression;
    for (uint i = 0; i < index_cols.size(); i++) {
      index = "index_" + table + "_" + index_cols[i];
      createIndex(index, table, index_cols[i]);
    }
    // Create special indexes on the species strings to accelerate database queries
    for (int e = 1; e < NUM_ELEMENTS; e++) {
      index = "index_" + table + "_" + vatom_symbol[e];
      index_expression = "INSTR(species, '\"" + vatom_symbol[e] + "\"')";
      createIndex(index, table, index_expression);
    }
    transaction(false);
    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " Finished building table " << table << std::endl;
  }

  // Patch -------------------------------------------------------------------

  //applyPatchFromJsonl///////////////////////////////////////////////////////
  // Applies database patches from data in jsonl format. Returns the number
  // of entries that were patched.
  uint AflowDB::applyPatchFromJsonl(const vector<string>& data) {
    stringstream message;

    uint nlines = data.size();
    if (nlines == 0) return 0;  // Nothing to do

    uint chunk_size = 1000;

    vector<string> schema_keys = getAllSchemaKeys();

    // Patch
    vector<string> keys, vals, cols, types;
    string auid = "";
    uint npatches = 0, nkeys = 0;
    transaction(true);
    for (uint l = 0; l < nlines; l++) {
      cols.clear();
      keys = aurostd::extractJsonKeysAflow(data[l]);
      nkeys = keys.size();

      // Check
      for (uint k = 0; k < nkeys; k++) {
        if (keys[k] == "auid") {
          auid = aurostd::extractJsonValueAflow(data[l], keys[k]);
          auid = auid.substr(1, auid.size() - 2);  // Remove quotes
        } else if (aurostd::WithinList(schema_keys, keys[k])) {
          cols.push_back(keys[k]);
        }
      }

      if (auid.empty()) {
        message << "Skipping line " << l << ". No AUID found in JSON.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }

      if (!auidInDatabase(auid)) {
        message << "Skipping line " << l << ". AUID " << auid << " not in database.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }

      if (cols.size() == 0) {
        message << "Skipping line " << l << ". No keys match properties in schema.";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }

      types = getDataTypes(cols, false);
      vals = getDataValues(data[l], cols, types);

      // Update
      try {
        updateEntry(auid, cols, vals);
        npatches++;
      } catch (aurostd::xerror& e) {
        message << "Failed to patch using line " << l << " (SQL error): " << e.buildMessageString() << ".";
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        continue;
      }

      if ((npatches + 1) % chunk_size == 0) {
        transaction(false);
        transaction(true);
      }
    }

    transaction(false);

    return npatches;
  }

  //auidInDatabase////////////////////////////////////////////////////////////
  // Test if an AUID is in the database by retrieving a value of the row (here
  // AUID since it is indexed).
  bool AflowDB::auidInDatabase(const string& auid) {
    string table = "auid_" + auid.substr(6, 2);
    string where = "auid='\"" + auid + "\"'";
    string result = getValue(table, "auid", where);
    return (!result.empty());
  }

  //updateEntry///////////////////////////////////////////////////////////////
  // Updates columns of an entry with a specific AUID.
  void AflowDB::updateEntry(const string& auid, const vector<string>& cols, const vector<string>& vals) {
    string table = "auid_" + auid.substr(6, 2);
    string where = "AUID='\"" + auid + "\"'";
    updateRow(table, cols, vals, where);
  }

  //getAllSchemaKeys//////////////////////////////////////////////////////////
  // Returns the keys from the schema and the internal schema
  vector<string> AflowDB::getAllSchemaKeys() {
    vector<string> keys = init::getSchemaKeys(vschema);
    vector<string> keys_internal = init::getSchemaKeys(vschema_internal);
    for (const string& key : keys_internal) keys.push_back(key);
    return keys;
  }

  // Data --------------------------------------------------------------------


  //getDataNames//////////////////////////////////////////////////////////////
  // Gets the data names of the schema keys and converts them into SQLite
  // types. Note that SQLite does not recognize arrays or Booleans, so they
  // will be stored as text.
  vector<string> getDataNames() { //CO20200520
    vector<string> keys;
    for (uint i = 0, n = XHOST.vschema.vxsghost.size(); i < n; i += 2) {
      if(XHOST.vschema.vxsghost[i].find("SCHEMA::NAME:")!=string::npos) {
        const string& key=XHOST.vschema.vxsghost[i+1];
        keys.push_back(key);
      }
    }
    return keys;
  }

  //getDataTypes//////////////////////////////////////////////////////////////
  // Gets the data types of the schema keys and converts them into SQLite
  // types. Note that SQLite does not recognize arrays or Booleans, so they
  // will be stored as text.
  vector<string> AflowDB::getDataTypes(const vector<string>& keys, bool unique) {
    uint nkeys = keys.size();
    vector<string> types(nkeys);
    string type = "";
    for (uint k = 0; k < nkeys; k++) {
      // AUID has to be unique
      type = vschema.getattachedscheme("SCHEMA::TYPE:" + aurostd::toupper(keys[k]));
      if (type.empty()) type = vschema_internal.getattachedscheme("SCHEMA::TYPE:" + aurostd::toupper(keys[k]));
      if (unique && (keys[k] == "AUID")) {
        types[k] = "TEXT UNIQUE NOT NULL COLLATE NOCASE";  // Make string search not case sensitive
      } else if (type == "number") {
        types[k] = "REAL";
      } else {
        types[k] = "TEXT COLLATE NOCASE";  // Make string search not case sensitive
      }
    }
    return types;
  }

  //getDataValues/////////////////////////////////////////////////////////////
  // Retrieves the values of each property from the aflowlib.json file.
  vector<string> AflowDB::getDataValues(const string& entry, const vector<string>& cols, const vector<string>& types) {
    string value = "", id = "";
    uint ncols = cols.size();
    vector<string> values(ncols, "NULL");
    string species = "";  // For special key handling
    for (uint c = 0; c < ncols; c++) {
      value = aurostd::extractJsonValueAflow(entry, cols[c]);

      // Store for later
      if (cols[c] == "species") species = value;
      else if (cols[c] == "aurl") {
        aurostd::StringSubst(value, "_RAW/", "_WEB/");
        aurostd::StringSubst(value, "_LIB/", "_WEB/");
      }

      // Check for synonyms for changed parameter names
      if (value.empty()) {
        if (cols[c] == "ldau_type") {
          value = "0"; // if no LDAU type in json, assume it's 0 (no LDAU)
        } else if (cols[c] == "aflow_prototype_label_relax") {
          value = aurostd::extractJsonValueAflow(entry, "anrl_label_relax"); //replace anrl_label_relax with aflow_prototype_label_relax - DO NOT TOUCH
        } else if (cols[c] == "aflow_prototype_label_orig") {
          value = aurostd::extractJsonValueAflow(entry, "anrl_label_orig"); //replace anrl_label_orig with aflow_prototype_label_orig - DO NOT TOUCH
        } else if (cols[c] == "aflow_prototype_params_list_relax") {
          value = aurostd::extractJsonValueAflow(entry, "anrl_parameter_list_relax"); //replace anrl_parameter_list_relax with aflow_prototype_params_list_relax - DO NOT TOUCH
          if (value.empty()) value = aurostd::extractJsonValueAflow(entry, "aflow_prototype_parameter_list_relax"); //replace anrl_prototype_parameter_list_relax with aflow_prototype_params_list_relax - DO NOT TOUCH
        } else if (cols[c] == "aflow_prototype_params_list_orig") {
          value = aurostd::extractJsonValueAflow(entry, "anrl_parameter_list_orig"); //replace anrl_parameter_list_orig with aflow_prototype_params_list_orig - DO NOT TOUCH
          if (value.empty()) value = aurostd::extractJsonValueAflow(entry, "aflow_prototype_parameter_list_orig"); //replace anrl_prototype_parameter_list_orig with aflow_prototype_params_list_orig - DO NOT TOUCH
        } else if (cols[c] == "aflow_prototype_params_values_relax") {
          value = aurostd::extractJsonValueAflow(entry, "anrl_parameter_values_relax"); //replace anrl_parameter_values_relax with aflow_prototype_params_values_relax - DO NOT TOUCH
          if (value.empty()) value = aurostd::extractJsonValueAflow(entry, "aflow_prototype_parameter_values_relax"); //replace aflow_prototype_parameter_values_relax with aflow_prototype_params_values_relax - DO NOT TOUCH
        } else if (cols[c] == "aflow_prototype_params_values_orig") {
          value = aurostd::extractJsonValueAflow(entry, "anrl_parameter_values_orig"); //replace anrl_parameter_values_orig with aflow_prototype_params_values_orig - DO NOT TOUCH
          if (value.empty()) value = aurostd::extractJsonValueAflow(entry, "aflow_prototype_parameter_values_orig"); //replace aflow_prototype_parameter_values_orig with aflow_prototype_params_values_orig - DO NOT TOUCH
        }
      }

      // If not found in the json, Check if column is part of the extra schema
      // Only check if not found in case the json is part of a patch file
      if (value.empty()) {
        if ((cols[c] == "alloy") && (!species.empty())) {
          for (uint i = 0; i < species.size(); i++) {
            if (isalpha(species[i])) value += species[i];
          }
        }
      }
      if (!value.empty() && (aurostd::toupper(value) != "NULL")) {
        if (types[c] != "REAL") values[c] = "'" + value + "'";
        else values[c] = value;
      }
    }
    return values;
  }

}  // namespace aflowlib

/**************************** DATABASE ANALYSIS *****************************/

namespace aflowlib {

  //analyzeDatabase///////////////////////////////////////////////////////////
  // Provides analytics for the database in JSON format.
  void AflowDB::analyzeDatabase(const string& outfile) {

    if (aurostd::FileEmpty(database_file)) {
      string message = "Cannot analyze database. File empty.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }

    // Get properties and tables for which statistics need to be collected
    // Since all tables have the same columns, only one table needs to be searched
    vector<string> tables = getTables();

    vector<string> catalogs = getSetMultiTables(tables, "catalog", true);
    uint ncatalogs = catalogs.size();

    vector<string> loop_entries, loops;
    loop_entries = getSetMultiTables(tables, "loop", true);
    loops = getUniqueFromJsonArrays(loop_entries);

    vector<DBStats> db_stats(ncatalogs + 1);
    DBStats& total_stats = db_stats[0];  // Declare to make code more legible
    total_stats = initDBStats("\"total\"", loops);

    for (uint c = 0; c < ncatalogs; c++) {
      DBStats& catalog_stats = db_stats[c+1];  // Declare to make code more legible
      catalog_stats = getCatalogStats(catalogs[c], tables, loops);

      // Add to totals
      total_stats.nentries += catalog_stats.nentries;
      // Since we cannot determine which systems are in the ICSD but not in
      // other catalogs, skip it for the system count. The error should be small.
      if (catalog_stats.catalog != "ICSD") total_stats.nsystems += catalog_stats.nsystems;

      // Columns
      for (uint i = 0; i < total_stats.columns.size(); i++) {
        total_stats.count[i][0] += catalog_stats.count[i][0];
        total_stats.count[i][1] += catalog_stats.count[i][1];
        if ((catalog_stats.types[i] != "bool") && (catalog_stats.count[i][0] + catalog_stats.count[i][1] > 0)) {
          if (total_stats.types[i] == "number") {
            if (total_stats.max[i].empty()) total_stats.max[i] = catalog_stats.max[i];
            else if (aurostd::string2utype<double>(catalog_stats.max[i]) > aurostd::string2utype<double>(total_stats.max[i])) total_stats.max[i] = catalog_stats.max[i];

            if (total_stats.min[i].empty()) total_stats.min[i] = catalog_stats.min[i];
            else if (aurostd::string2utype<double>(total_stats.min[i]) > aurostd::string2utype<double>(catalog_stats.min[i])) total_stats.min[i] = catalog_stats.min[i];
          }

          // Just append to the sets for now and sort later
          uint nset = total_stats.set[i].size();
          for (uint s = 0; s < catalog_stats.set[i].size(); s++) {
            if ((nset <= _DEFAULT_SET_LIMIT_) && !aurostd::WithinList(total_stats.set[i], catalog_stats.set[i][s])) {
              total_stats.set[i].push_back(catalog_stats.set[i][s]);
              nset++;
            }
          }
        }
      }

      // Loops
      for (std::map<string, uint>::iterator it = total_stats.loop_counts.begin(); it != total_stats.loop_counts.end(); ++it) {
        const string& key = (*it).first;
        (*it).second += catalog_stats.loop_counts[key];
      }

      // Species
      for (uint s = 0; s < catalog_stats.species.size(); s++) {
        if (!aurostd::WithinList(total_stats.species, catalog_stats.species[s])) {
          total_stats.species.push_back(catalog_stats.species[s]);
        }
      }
    }

    // Get the sets for the totals
    for (uint i = 0; i < total_stats.columns.size(); i++) {
      uint nset = total_stats.set[i].size();
      if (nset <= _DEFAULT_SET_LIMIT_) {
        if (total_stats.types[i] == "bool") {
          total_stats.set[i].clear();
          total_stats.set[i].push_back("true");
          total_stats.set[i].push_back("false");
        } else if (total_stats.types[i] == "number") {
          vector<double> set_dbl(nset);
          for (uint s = 0; s < nset; s++) set_dbl[s] = aurostd::string2utype<double>(total_stats.set[i][s]);
          aurostd::sort(set_dbl, total_stats.set[i]);
        } else {
          std::sort(total_stats.set[i].begin(), total_stats.set[i].end());
        }
      }
    }
    std::sort(total_stats.species.begin(), total_stats.species.end());

    // Write output
    ncatalogs = ncatalogs + 1; // Include totals

    std::stringstream json;
    json << "{\"Aflow_DBs\":{";

    for (uint c = 0; c < ncatalogs; c++) {
      json << db_stats[c].catalog << ":" << stats2json(db_stats[c]) << ((c < ncatalogs - 1)?",":"");
    }
    json << "}}" << std::endl;
    aurostd::stringstream2file(json, outfile);
  }

  DBStats AflowDB::initDBStats(const string& catalog, const vector<string>& loops) {
    vector<string> keys = getAllSchemaKeys();
    vector<string> excluded_properties, cols;
    string exclude = "alloy";
    aurostd::string2tokens(exclude, excluded_properties, ",");
    string key = "";
    for (uint i = 0, nkeys = keys.size(); i < nkeys; i++) {
      key = vschema.getattachedscheme("SCHEMA::NAME:" + keys[i]);
      if (key.empty()) key = vschema_internal.getattachedscheme("SCHEMA::NAME:" + keys[i]);
      if (!key.empty() && !aurostd::WithinList(excluded_properties, key)) cols.push_back(key);
    }
    uint ncols = cols.size();

    DBStats stats;
    stats.catalog = catalog;
    stats.columns = cols;
    stats.count.assign(ncols, vector<int>(2, 0));
    stats.max.assign(ncols, "");
    stats.min.assign(ncols, "");
    stats.nentries = 0;
    stats.nsystems = 0;
    stats.set.resize(ncols);
    for (uint i = 0; i < loops.size(); i++) {
      stats.loop_counts[loops[i]] = 0;
    }

    // Get types for post-processing
    vector<string> types(ncols);
    for (uint c = 0; c < ncols; c++) {
      types[c] = vschema.getattachedscheme("SCHEMA::TYPE:" + aurostd::toupper(cols[c]));
    }
    stats.types = types;

    return stats;
  }

  //getCatalogStats///////////////////////////////////////////////////////////
  // Gets the statistics for all properties in the catalog.
  DBStats AflowDB::getCatalogStats(const string& catalog, const vector<string>& tables, const vector<string>& loops) {
    stringstream message;

    DBStats stats = initDBStats(catalog, loops);
    vector<DBStats> colstats(tables.size(), stats);

    string where = "catalog='" + catalog + "'";
    vector<string> entries = getDatabasePropertyMultiTables("COUNT", tables, "*", where);
    uint ntables = tables.size();
    for (uint t = 0; t < ntables; t++) stats.nentries += aurostd::string2utype<int>(entries[t]);
    message << "Starting analysis for catalog " << catalog << " (" << stats.nentries << " entries).";
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, *p_FileMESSAGE, *p_oss);

    if (stats.nentries > 0) {
#ifdef AFLOW_MULTITHREADS_ENABLE
      int ncpus = init::GetCPUCores();
      if (ncpus < 1) ncpus = 1;
      // The maximum number of CPUs are empirically found values that appear
      // to result in the shortest run times. Further testing may be necessary.
      int max_cpus = 32;
      if (stats.nentries < 10000) max_cpus = 4;
      else if (stats.nentries < 50000)  max_cpus = 8;
      else if (stats.nentries < 100000) max_cpus = 16;
      if (ncpus > max_cpus) ncpus = max_cpus;
      xthread::xThread xt(ncpus);
      std::function<void(uint, int, const vector<string>&, vector<DBStats>&)> fn = std::bind(&AflowDB::getColStats, this, _1, _2, _3, _4);
      xt.runPredistributed(ntables, fn, tables, colstats);
#else
      getColStats(0, ntables, tables, colstats);
#endif

      // Properties: count, max, min, set
      vector<string> set;
      string max = "", min = "";
      uint nset = 0, n = 0;
      uint ncols = stats.columns.size();
      const vector<string>& types = stats.types;
      for (uint c = 0; c < ncols; c++) {
        for (uint t = 0; t < ntables; t++) {
          stats.count[c][0] += colstats[t].count[c][0];
          stats.count[c][1] += colstats[t].count[c][1];
        }
        if (stats.count[c][0] + stats.count[c][1] > 0) {
          set.clear();
          max = ""; min = "";
          nset = 0; n = 0;
          if (types[c] != "bool") {  // No max, min, or set for bool
            for (uint t = 0; t < ntables; t++) {
              const DBStats& cstats = colstats[t];
              if (cstats.count[c][0] > 0) {
                if (types[c] == "number") {
                  if (max.empty()) max = cstats.max[c];
                  else if (aurostd::string2utype<double>(cstats.max[c]) > aurostd::string2utype<double>(max)) max = cstats.max[c];

                  if (min.empty()) min = cstats.min[c];
                  else if (aurostd::string2utype<double>(cstats.min[c]) < aurostd::string2utype<double>(min)) min = cstats.min[c];
                }
                if (nset <= _DEFAULT_SET_LIMIT_) {
                  n = cstats.set[c].size();
                  if (n > _DEFAULT_SET_LIMIT_) {
                    set = cstats.set[c];
                    nset = n;
                  } else {
                    for (const string& s : cstats.set[c]) {
                      if (!aurostd::WithinList(set, s)) {
                        set.push_back(s);
                        nset++;
                      }
                      if (nset > _DEFAULT_SET_LIMIT_) break;
                    }
                  }
                }
              }
            }
          }
          stats.max[c] = max;
          stats.min[c] = min;
          if (types[c] == "bool") {
            set.clear();
            set.push_back("true");
            set.push_back("false");
          } else if (nset <= _DEFAULT_SET_LIMIT_) {
            if (types[c] == "number") {
              vector<double> set_dbl(nset);
              for (uint i = 0; i < nset; i++) set_dbl[i] = aurostd::string2utype<double>(set[i]);
              aurostd::sort(set_dbl, set);
            } else {
              std::sort(set.begin(), set.end());
            }
          }
          stats.set[c] = set;
        }
      }

      // Species, systems
      vector<string> species = getSetMultiTables(tables, "species", true, where);
      stats.nsystems = species.size();
      stats.species = getUniqueFromJsonArrays(species);

      // Loop counts
      for (std::map<string, uint>::iterator it = stats.loop_counts.begin(); it != stats.loop_counts.end(); ++it) {
         const string& key = (*it).first;
         for (uint t = 0; t < ntables; t++) (*it).second += colstats[t].loop_counts[key];
      }
    }

    return stats;
  }

  //getColStats///////////////////////////////////////////////////////////////
  // Retrieves the statistics for each database property and the loops.
  void AflowDB::getColStats(int startIndex, int endIndex, const vector<string>& tables, vector<DBStats>& colstats) {
    sqlite3* cursor;
    string message = "";
    int sql_code = sqlite3_open_v2(database_file.c_str(), &cursor, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, NULL); //DX20200319 - nullptr -> NULL
    if (sql_code != SQLITE_OK) {
      message = "Could not open cursor on database file " + database_file + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    string where = "";
    for (int i = startIndex; i < endIndex; i++) {
      // Writing directly to the colstats vector may not be thread-safe,
      // so create copy first and merge later.
      DBStats cstats = colstats[i];
      const string& catalog = cstats.catalog;
      const vector<string>& cols = cstats.columns;
      const vector<string>& types = cstats.types;
      uint ncols = cols.size();
      for (uint c = 0; c < ncols; c++) {
        if (types[c] == "bool") {
          where = "catalog='" + catalog + "' AND " + cols[c] + "=";
          cstats.count[c][0] = aurostd::string2utype<int>(getDatabaseProperty(cursor, "COUNT", tables[i], cols[c], where + "'true'"));
          cstats.count[c][1] = aurostd::string2utype<int>(getDatabaseProperty(cursor, "COUNT", tables[i], cols[c], where + "'false'"));
        } else {
          where = "catalog='" + catalog + "' AND " + cols[c] + " NOT NULL";
          cstats.count[c][0] = aurostd::string2utype<int>(getDatabaseProperty(cursor, "COUNT", tables[i], cols[c], where));
        }
        // No need to determine max, min, or set for bool
        if ((types[c] != "bool") && (cstats.count[c][0] > 0)) {
          // Max and  min only make sense for numbers
          if (types[c] == "number") {
            cstats.max[c] = getDatabaseProperty(cursor, "MAX", tables[i], cols[c], where);
            cstats.min[c] = getDatabaseProperty(cursor, "MIN", tables[i], cols[c], where);
          }
          cstats.set[c] = getSet(cursor, tables[i], cols[c], true, where, _DEFAULT_SET_LIMIT_ + 1);
        }
      }
      for (std::map<string, uint>::iterator it = cstats.loop_counts.begin(); it != cstats.loop_counts.end(); ++it) {
        const string& key = (*it).first;
        where = "catalog='" + catalog + "' AND loop LIKE '%\"" + key + "\"%'";
        (*it).second = aurostd::string2utype<int>(getDatabaseProperty("COUNT", tables[i], "loop", where));
      }
#ifdef AFLOW_MULTITHREADS_ENABLE
      std::lock_guard<std::mutex> lk(write_mutex);
#endif
      colstats[i] = std::move(cstats);
    }

    sql_code = sqlite3_close(cursor);
    if (sql_code != SQLITE_OK) {
      message = "Could not close cursor on database file " + database_file + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  //getUniqueFromJsonArrays///////////////////////////////////////////////////
  // Determines the unique array elements in a set of 1D-array strings.
  vector<string> AflowDB::getUniqueFromJsonArrays(const vector<string>& arrays) {
    vector<string> unique, tokens;
    string arr = "";
    int nunique = 0, u = 0;
    for (uint a = 0; a < arrays.size(); a++) {
      vector<string> tokens;
      arr = arrays[a];
      arr=aurostd::RemoveSubString(arr, "[");
      arr=aurostd::RemoveSubString(arr, "]");
      arr=aurostd::RemoveSubString(arr, "\"");
      aurostd::string2tokens(arr, tokens, ", ");
      for (uint t = 0; t < tokens.size(); t++) {
        if (nunique == 0) {
          unique.push_back(tokens[t]);
          nunique = 1;
        } else {
          for (u = 0; u < nunique; u++) {
            if (unique[u] == tokens[t]) break;
          }
          if (u == nunique) {
            unique.push_back(tokens[t]);
            nunique++;
          }
        }
      }
    }
    return unique;
  }

  //writeStatsToJson//////////////////////////////////////////////////////////
  // Writes the database statistics into a JSON-formatted string(stream).
  string AflowDB::stats2json(const DBStats& db_stats) {
    stringstream json;
    json << "{";
    json << "\"count\":" << db_stats.nentries << ",";
    json << "\"systems\":" <<  db_stats.nsystems << ",";
    for (std::map<string, uint>::const_iterator it = db_stats.loop_counts.begin(); it != db_stats.loop_counts.end(); ++it) {
      json << "\"" << (*it).first << "\":" << (*it).second << ",";
    }
    json << "\"columns\":{";
    uint ncols = db_stats.columns.size();
    for (uint c = 0; c < ncols; c++) {
      json << "\"" << db_stats.columns[c] << "\":{";
      json << "\"count\":";
      if (db_stats.types[c] == "bool") {
        json << "{\"true\":" << db_stats.count[c][0] << ",\"false\":" << db_stats.count[c][1] << "},";
      } else {
        json << aurostd::utype2string<int>(db_stats.count[c][0]) << ",";
      }
      json << "\"min\":" << (db_stats.min[c].empty()?"null":db_stats.min[c])<< ",";
      json << "\"max\":" << (db_stats.max[c].empty()?"null":db_stats.max[c]) << ",";

      // Write set
      json << "\"set\":";
      if (db_stats.set[c].size() > _DEFAULT_SET_LIMIT_) json << "null";
      else json << "[" << aurostd::joinWDelimiter(db_stats.set[c], ",") << "]";

      json << "}" << ((c < ncols - 1)?",":"");
    }
    json << "},";

    json << "\"species\":";
    if (db_stats.species.size() > 0) {
      json << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(db_stats.species, "\"", "\""), ",") << "]";
    } else {
      json << "null";
    }

    json << "}";
    return json.str();
  }

}  // namespace aflowlib

/******************************* GET ENTRIES ********************************/

namespace aflowlib {

  //getEntry//////////////////////////////////////////////////////////////////
  // Returns all data for a specific AUID in the specified format
  string AflowDB::getEntry(const string& auid, filetype ft) {
    string entry = "";
    if (auidInDatabase(auid)) {
      string table = "auid_" + auid.substr(6,2);
      string where = "auid='\"" + auid + "\"'";
      vector<string> values = getRows(table, where)[0];
      if (values.size() == 0) return "";  // Nothing found
      vector<string> keys = getColumnNames(table);
      uint nkeys = keys.size();
      vector<string> keyval;
      switch (ft) {
        case json_ft:
          for (uint i = 0; i < nkeys; i++) {
            if (!values[i].empty()) {
              keyval.push_back("\"" + keys[i] + "\":" + values[i]);
            }
          }
          entry = "{" + aurostd::joinWDelimiter(keyval, ",") + "}";
          break;
        case aflow_ft:
        default:
          for (uint i = 0; i < nkeys; i++) {
            if (!values[i].empty()) {
              aurostd::StringSubst(values[i], "],[", ";");
              aurostd::StringSubst(values[i], "[", "");
              aurostd::StringSubst(values[i], "]", "");
              aurostd::StringSubst(values[i], "\"", "");
              keyval.push_back(keys[i] + "=" + values[i]);
            }
          }
          entry = aurostd::joinWDelimiter(keyval, "|");
      }
    }
    return entry;
  }


  //getEntryAentry////////////////////////////////////////////////////////////
  // Returns the data of a specific auid as an _aflowlib_entry object
  _aflowlib_entry AflowDB::getEntryAentry(const string& auid) {
    _aflowlib_entry aentry;
    string aflowout = getEntry(auid, aflow_ft);
    if (!aflowout.empty()) aentry.Load(aflowout, *p_oss);
    return aentry;
  }

  //getEntrySet///////////////////////////////////////////////////////////////
  // Returns all data of a set of auid based on a search condition
  vector<string > AflowDB::getEntrySet(const string& where, filetype ft) {
    vector<vector<string> > data = getRowsMultiTables(where);
    uint ndata = data.size();
    vector<string> entries(ndata);
    if (ndata > 0) {
      vector<string> keys = getColumnNames("auid_00");
      uint nkeys = keys.size();
      vector<string> keyval;
      for (uint i = 0; i < ndata; i++) {
        keyval.clear();
        switch (ft) {
          case json_ft:
            for (uint j = 0; j < nkeys; j++) {
              if (!data[i][j].empty()) {
                keyval.push_back("\"" + keys[j] + "\":" + data[i][j]);
              }
            }
            entries[i] = "{" + aurostd::joinWDelimiter(keyval, ",") + "}";
            break;
          case aflow_ft:
          default:
            for (uint j = 0; j < nkeys; j++) {
              if (!data[i][j].empty()) {
                aurostd::StringSubst(data[i][j], "],[", ";");
                aurostd::StringSubst(data[i][j], "[", "");
                aurostd::StringSubst(data[i][j], "]", "");
                aurostd::StringSubst(data[i][j], "\"", "");
                keyval.push_back(keys[j] + "=" + data[i][j]);
              }
            }
            entries[i] = aurostd::joinWDelimiter(keyval, "|");
        }
      }
    }
    return entries;
  }

  vector<_aflowlib_entry> AflowDB::getEntrySetAentry(const string& where) {
    vector<string> aflowouts = getEntrySet(where, aflow_ft);
    uint nentries = aflowouts.size();
    vector<_aflowlib_entry> aentries(nentries);
    for (uint i = 0; i < nentries; i++) aentries[i].Load(aflowouts[i], *p_oss);
    return aentries;
  }

}

/***************************** SQLite FUNCTIONS *****************************/

// These functions are higher level SQLite functions that call functions of
// the SQLite interface. They are essentially syntactic sugar to make the code
// easier to read and to facilitate the implementation of the AflowDB class
// into other parts of AFLOW.
//
// Getter functions are overloaded to take a cursor other than the default
// cursor. Since SQLite does not allow concurrent writing, the writer
// functions should not be overloaded.

namespace aflowlib {

  // INDEX -------------------------------------------------------------------

  //createIndex///////////////////////////////////////////////////////////////
  // Creates an index.
  void AflowDB::createIndex(const string& index, const string& table, const string& column) {
    string command = "CREATE INDEX " + index + " ON " + table + "(" + column +  ")";
    sql::SQLexecuteCommand(db, command);
  }

  //dropIndex/////////////////////////////////////////////////////////////////
  // Removes an index.
  void AflowDB::dropIndex(const string& index) {
    string command = "DROP INDEX " + index;
    sql::SQLexecuteCommand(db, command);
  }

  // TRANSACTION -------------------------------------------------------------

  //transaction///////////////////////////////////////////////////////////////
  // Begings (begin == true) or ends a database transaction.
  void AflowDB::transaction(bool begin) {
    string command = string(begin?"BEGIN":"END") + " TRANSACTION";
    sql::SQLexecuteCommand(db, command);
  }

  // TABLE -------------------------------------------------------------------

  //dropTable/////////////////////////////////////////////////////////////////
  // Deletes a table from the database.
  void AflowDB::dropTable(const string& table) {
    string command = "DROP TABLE IF EXISTS " + table;
    sql::SQLexecuteCommand(db, command);
  }

  //getTables/////////////////////////////////////////////////////////////////
  // Retrieves a set of tables. If where is empty, all tables in the database
  // will be returned.
  vector<string> AflowDB::getTables(const string& where) {
    return getTables(db, where);
  }

  vector<string> AflowDB::getTables(sqlite3* cursor, const string& where) {
    string command = "SELECT name FROM sqlite_master WHERE type='table'";
    if (!where.empty()) command += " AND (" + where + ")";
    return sql::SQLexecuteCommandVECTOR(cursor, command);
  }

  //createTable///////////////////////////////////////////////////////////////
  // Creates a table where all columns have the same type.
  void AflowDB::createTable(const string& table, const vector<string>& cols, const string& type) {
    vector<string> types(cols.size(), type);
    createTable(table, cols, types);
  }

  // Creates a table where each column is assigned its own type
  void AflowDB::createTable(const string& table, const vector<string>& cols, const vector<string>& types) {
    uint ncols = cols.size();
    if (ncols != types.size()) {
      string message = "Could not create table. ";
      message += "Number of columns and number of types do not match.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      string command = "CREATE TABLE IF NOT EXISTS " + table + " (";
      for (uint c = 0; c < ncols; c++) {
        command += cols[c] + " " + types[c];
        if (c < ncols - 1) command += ", ";
      }
      command += ")";
      sql::SQLexecuteCommand(db, command);
    }
  }

  // INSERT ------------------------------------------------------------------

  //insertValues//////////////////////////////////////////////////////////////
  // Inserts a set of values into a table.
  void AflowDB::insertValues(const string& table, const vector<string>& vals) {
    vector<string> cols;
    insertValues(table, cols, vals);
  }

  // Inserts a set of values into a table (if cols is empty) or into specific
  // columns of the table.
  void AflowDB::insertValues(const string& table, const vector<string>& cols,
      const vector<string>& vals) {
    uint ncols = cols.size();
    uint nvals = vals.size();
    if ((ncols > 0) && (ncols != nvals)) {
      string message = "Could not insert values. ";
      message += "Number of columns and number of values do not match.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      string command = "INSERT INTO " + table;
      if (ncols > 0) command += " (" + aurostd::joinWDelimiter(cols, ", ") + ")";
      command += " VALUES(" + aurostd::joinWDelimiter(vals, ", ") + ")";
      sql::SQLexecuteCommand(db, command);
    }
  }

  // UPDATE ------------------------------------------------------------------
  // Performs an UPDATE command.
  // WARNING: if WHERE is empty, all rows will updated! Never default to an
  // empty string!
  void AflowDB::updateRow(const string& table, const vector<string>& cols,
      const vector<string>& vals, const string& where) {
    uint ncols = cols.size();
    uint nvals = vals.size();
    if (ncols == 0) {
      string message = "Could not update row. ";
      message += "No columns selected.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else if (ncols != nvals) {
      string message = "Could not update row. ";
      message += "Number of columns and number of values do not match.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    } else {
      string command = "UPDATE " + table + " SET ";
      for (uint i = 0; i < ncols; i++) {
        command += cols[i] + "=" + vals[i];
        if (i < ncols - 1) command += ",";
      }
      if (!where.empty()) command += " WHERE " + where;
      sql::SQLexecuteCommand(db, command);
    }
  }


  // GET ---------------------------------------------------------------------

  //getColumnNames////////////////////////////////////////////////////////////
  // Returns the names of all columns in a specific table.
  vector<string> AflowDB::getColumnNames(const string& table) {
    return getColumnNames(db, table);
  }

  vector<string> AflowDB::getColumnNames(sqlite3* cursor, const string& table) {
    string command = "PRAGMA table_info(" + table + ")";
    vector<vector<string> > pragma_results = sql::SQLexecuteCommand2DVECTOR(cursor, command);
    uint nresults = pragma_results.size();
    vector<string> columns(nresults);
    for (uint r = 0; r < nresults; r++) {
      columns[r] = pragma_results[r][1];
    }
    return columns;
  }

  //getColumnTypes////////////////////////////////////////////////////////////
  // Returns the data types of all columns in a specific table.
  vector<string> AflowDB::getColumnTypes(const string& table) {
    return getColumnTypes(db, table);
  }

  vector<string> AflowDB::getColumnTypes(sqlite3* cursor, const string& table) {
    string command = "PRAGMA table_info(" + table + ")";
    vector<vector<string> > pragma_results = sql::SQLexecuteCommand2DVECTOR(cursor, command);
    uint nresults = pragma_results.size();
    vector<string> columns(nresults);
    for (uint r = 0; r < nresults; r++) {
      columns[r] = pragma_results[r][2];
    }
    return columns;
  }

  //getRows///////////////////////////////////////////////////////////////////
  // Returns the rows for a specific where condition. If the where condition is
  // ambiguous, the first row that fits the criteria will be returned.
  vector<vector<string> > AflowDB::getRows(const string& table, const string& where) {
    return getRows(db, table, where);
  }

  vector<vector<string> > AflowDB::getRows(sqlite3* cursor, const string& table, const string& where) {
    string command = prepareSELECT(table, "", "*", where, 0, "");
    return sql::SQLexecuteCommand2DVECTOR(cursor, command);
  }

  vector<vector<string> > AflowDB::getRowsMultiTables(const string& where) {
    vector<string> tables = getTables();
    return getRowsMultiTables(db, tables, where);
  }

  vector<vector<string> > AflowDB::getRowsMultiTables(sqlite3* cursor, const string& where) {
    vector<string> tables = getTables();
    return getRowsMultiTables(cursor, tables, where);
  }

  vector<vector<string> > AflowDB::getRowsMultiTables(const vector<string>& tables, const string& where) {
    return getRowsMultiTables(db, tables, where);
  }

  vector<vector<string> > AflowDB::getRowsMultiTables(sqlite3* cursor, const vector<string>& tables, const string& where) {
    uint ntables = tables.size();
    vector<string> commands(ntables);
    for (uint t = 0; t < ntables; t++) commands[t] = prepareSELECT(tables[t], "", "*", where);
    return sql::SQLexecuteCommand2DVECTOR(cursor, aurostd::joinWDelimiter(commands, " UNION ALL "));
  }

  //getValue//////////////////////////////////////////////////////////////////
  // Gets a value from a specific column. The row must be specified in the
  // where condition, or else it just takes the first value.
  string AflowDB::getValue(const string& table, const string& col, const string& where) {
    return getValue(db, table, col, where);
  }

  string AflowDB::getValue(sqlite3* cursor, const string& table, const string& col, const string& where) {
    string command = prepareSELECT(table, "", col, where, 0, "");
    return sql::SQLexecuteCommandSCALAR(cursor, command);
  }

  // Collects the values from a single column from all tables. //HE20220405
  std::vector<std::string> AflowDB::getValuesMultiTable(const std::string & col, const std::string & where){
    vector<string> tables = getTables();
    uint ntables = tables.size();
    vector<string> commands(ntables);
    for (uint t = 0; t < ntables; t++) commands[t] = prepareSELECT(tables[t], "", col, where);
    return sql::SQLexecuteCommandVECTOR(db, aurostd::joinWDelimiter(commands, " UNION ALL "));
  }


  //getDatabaseProperty///////////////////////////////////////////////////////
  // Gets a database property for a specific column.
  string AflowDB::getDatabaseProperty(const string& property, const string& table,
      const string& col, const string& where) {
    return getDatabaseProperty(db, property, table, col, where);
  }
  string AflowDB::getDatabaseProperty(sqlite3* cursor, const string& property, const string& table,
      const string& col, const string& where) {
    string command = prepareSELECT(table, property, col, where, 0, "");
    return sql::SQLexecuteCommandSCALAR(cursor, command);
  }

  //getDatabasePropertyMultiTables////////////////////////////////////////////
  // Gets a database property for a specific column across multiple tables.
  vector<string> AflowDB::getDatabasePropertyMultiTables(const string& property, const vector<string>& tables,
      const string& col, const string& where) {
    return getDatabasePropertyMultiTables(db, property, tables, col, where);
  }

  vector<string> AflowDB::getDatabasePropertyMultiTables(sqlite3* cursor, const string& property, const vector<string>& tables,
      const string& col, const string& where) {
    uint ntables = tables.size();
    vector<string> commands(ntables);
    for (uint t = 0; t < ntables; t++) commands[t] = prepareSELECT(tables[t], property, col, where);
    return sql::SQLexecuteCommandVECTOR(cursor, aurostd::joinWDelimiter(commands, " UNION ALL "));
  }

  //getSet////////////////////////////////////////////////////////////////////
  // Retrieves a (distinct) set from a single column.
  vector<string> AflowDB::getSet(const string& table, const string& col, bool distinct,
      const string& where, int limit, const string& order_by) {
    return getSet(db, table, col, distinct, where, limit, order_by);
  }

  vector<string> AflowDB::getSet(sqlite3* cursor, const string& table, const string& col, bool distinct,
      const string& where, int limit, const string& order_by) {
    string property = string((distinct?"DISTINCT":""));
    string command = prepareSELECT(table, property, col, where, limit, order_by);
    return sql::SQLexecuteCommandVECTOR(cursor, command);
  }

  //getSetMulitTables/////////////////////////////////////////////////////////
  // Retrieves a (distinct) set from a single column across multiple tables.
  // The result is sorted already, so there is not need for order_by.
  vector<string> AflowDB::getSetMultiTables(const vector<string>& tables, const string& col,
      bool distinct, const string& where, int limit) {
    return getSetMultiTables(db, tables, col, distinct, where, limit);
  }
  vector<string> AflowDB::getSetMultiTables(sqlite3* cursor, const vector<string>& tables, const string& col,
      bool distinct, const string& where, int limit) {
    string property = string((distinct?"DISTINCT":""));
    uint ntables = tables.size();
    vector<string> commands(ntables);
    for (uint t = 0; t < ntables; t++) commands[t] = prepareSELECT(tables[t], property, col, where, 0);
    string union_string = " UNION ";
    if (!distinct) union_string += "ALL ";
    string command = aurostd::joinWDelimiter(commands, union_string);
    if (limit > 0) command += " LIMIT " + aurostd::utype2string<int>(limit);
    return sql::SQLexecuteCommandVECTOR(cursor, command);
  }

  // SELECT ------------------------------------------------------------------

  //prepateSELECT/////////////////////////////////////////////////////////////
  // Lower level function to prepare a SELECT statement for all GET functions.
  string AflowDB::prepareSELECT(const string& table, const string& property, const string& cols,
      const string& where, int limit, const string& order_by) {
    stringstream command;
    command << "SELECT ";
    if (!property.empty()) command << property << ((property == "DISTINCT")?" ":"(");
    command << cols;
    if (!property.empty()) command << ((property == "DISTINCT")?"":")");
    command << " FROM " << table;
    if (!where.empty()) command << " WHERE (" << where << ")";
    if (!order_by.empty()) command << " ORDER BY " << order_by;
    if (limit > 0) command << " LIMIT " << limit;
    return command.str();
  }

  // Wrapper function for a vector representation of the columns.
  string AflowDB::prepareSELECT(const string& table, const string& property, const vector<string>& cols,
      const string& where, int limit, const string& order_by) {
    return prepareSELECT(table, property, aurostd::joinWDelimiter(cols, ", "), where, limit, order_by);
  }

}  // namespace aflowlib

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
//****************************************************************************
