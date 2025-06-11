
#ifndef AFLOWLIB_DATABASE_H
#define AFLOWLIB_DATABASE_H

#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"

#include "aflowlib/aflowlib_web_interface.h"
#include "extern/SQLITE/sqlite3.h"
#include "flow/aflow_support_types.h"
#include "flow/aflow_xclasses.h"

namespace aflowlib {
  struct DBStats {
    std::vector<std::string> columns;
    std::vector<std::vector<int>> count; // 2D to accommodate bool
    std::map<std::string, uint> loop_counts;
    std::vector<std::string> max;
    std::vector<std::string> min;
    int nentries;
    int nsystems;
    std::vector<std::vector<std::string>> set;
    std::vector<std::string> species;
    std::vector<std::string> types;
    std::string catalog;
  };

  class AflowDB : public xStream {
  public:
    AflowDB(const std::string&, std::ostream& oss = std::cout);
    AflowDB(const std::string&, const aurostd::xoption&, std::ostream& oss = std::cout);
    AflowDB(const std::string&, const aurostd::xoption&, const aurostd::xoption&, std::ostream& oss = std::cout);
    AflowDB(const std::string&, const std::string&, const std::string&, const aurostd::xoption&, std::ostream& oss = std::cout);
    AflowDB(const std::string&, const std::string&, const std::string&, const aurostd::xoption&, const aurostd::xoption&, std::ostream& oss = std::cout);
    AflowDB(const AflowDB&);
    AflowDB& operator=(const AflowDB&);
    ~AflowDB();
    void clear();

    [[nodiscard]] bool isTMP() const;

    int rebuildDatabase(bool force_rebuild = false);
    int rebuildDatabase(const std::string&, bool force_rebuild = false);
    int rebuildDatabase(const std::vector<std::string>&, bool force_rebuild = false);
    int patchDatabase(const std::string&, bool check_timestamps = false);
    int patchDatabase(const std::vector<std::string>&, bool check_timestamps = false);
    void analyzeDatabase(const std::string&);

    std::string getEntry(const std::string&, filetype);
    _aflowlib_entry getEntryAentry(const std::string&);
    std::vector<std::string> getEntrySet(const std::string&, filetype);
    std::vector<_aflowlib_entry> getEntrySetAentry(const std::string&);

    std::vector<std::vector<std::string>> getEntrySetData(const std::string&);

    std::vector<std::string> getTables(const std::string& where = "");
    std::vector<std::string> getTables(sqlite3*, const std::string& where = "");

    std::vector<std::string> getColumnNames(const std::string&);
    std::vector<std::string> getColumnNames(sqlite3*, const std::string&);
    std::vector<std::string> getColumnTypes(const std::string&);
    std::vector<std::string> getColumnTypes(sqlite3*, const std::string&);

    std::vector<std::vector<std::string>> getRows(const std::string&, const std::string& where = "");
    std::vector<std::vector<std::string>> getRows(sqlite3*, const std::string&, const std::string& where = "");
    std::vector<std::vector<std::string>> getRowsMultiTables(const std::string& where = "");
    std::vector<std::vector<std::string>> getRowsMultiTables(sqlite3*, const std::string& where = "");
    std::vector<std::vector<std::string>> getRowsMultiTables(const std::vector<std::string>&, const std::string& where = "");
    std::vector<std::vector<std::string>> getRowsMultiTables(sqlite3*, const std::vector<std::string>&, const std::string& where = "");
    std::string getValue(const std::string&, const std::string&, const std::string& where = "");
    std::string getValue(sqlite3*, const std::string&, const std::string&, const std::string& where = "");
    std::vector<std::string> getValuesMultiTable(const std::string& col, const std::string& where); // HE20220405
    std::string getDatabaseProperty(const std::string&, const std::string&, const std::string&, const std::string& where = "");
    std::string getDatabaseProperty(sqlite3*, const std::string&, const std::string&, const std::string&, const std::string& where = "");
    std::vector<std::string> getDatabasePropertyMultiTables(const std::string&, const std::vector<std::string>&, const std::string&, const std::string& where = "");
    std::vector<std::string> getDatabasePropertyMultiTables(sqlite3*, const std::string&, const std::vector<std::string>&, const std::string&, const std::string& where = "");
    std::vector<std::string> getSet(const std::string&, const std::string&, bool distinct = false, const std::string& where = "", int limit = 0, const std::string& order_by = "");
    std::vector<std::string> getSet(sqlite3*, const std::string&, const std::string&, bool distinct = false, const std::string& where = "", int limit = 0, const std::string& order_by = "");
    std::vector<std::string> getSetMultiTables(const std::vector<std::string>&, const std::string&, bool distinct = false, const std::string& where = "", int limit = 0);
    std::vector<std::string> getSetMultiTables(sqlite3*, const std::vector<std::string>&, const std::string&, bool distinct = false, const std::string& where = "", int limit = 0);

    void transaction(bool);

  private:
    void free();
    void open(int = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX);
    void close();
    void copy(const AflowDB&);
    void initialize(const std::string& db_file, const std::string& dt_path, const std::string& lck_file, int open_flags, const aurostd::xoption& schema_in, const aurostd::xoption& schema_internal_in);

    sqlite3* db;
    bool is_tmp;
    aurostd::xoption vschema_internal;
    aurostd::xoption vschema;
    std::string data_path;
    std::string database_file;
    std::string lock_file;
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex write_mutex;
#endif

    void openTmpFile(int open_flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, bool copy_original = false);
    bool closeTmpFile(bool force_copy = false, bool keep = false, bool nocopy = false);

    void rebuildDB();
    void buildTable(int, const std::vector<std::string>&, const std::vector<std::string>&);
    void populateTable(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<std::vector<std::string>>&);

    uint applyPatchFromJsonl(const std::vector<std::string>&);
    bool auidInDatabase(const std::string&);
    void updateEntry(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&);

    std::vector<std::string> getAllSchemaKeys();
    std::vector<std::string> getDataTypes(const std::vector<std::string>&, bool);
    std::vector<std::string> getDataValues(const aurostd::JSON::object&, const std::vector<std::string>&, const std::vector<std::string>&);

    DBStats initDBStats(const std::string&, const std::vector<std::string>&);
    DBStats getCatalogStats(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&);
    void getColStats(int, int, const std::vector<std::string>&, std::vector<DBStats>&);
    std::vector<std::string> getUniqueFromJsonArrays(const std::vector<std::string>&);
    std::string stats2json(const DBStats&);

    void createIndex(const std::string&, const std::string&, const std::string&);
    void dropIndex(const std::string&);
    void dropTable(const std::string&);
    void createTable(const std::string&, const std::vector<std::string>&, const std::string&);
    void createTable(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&);
    void insertValues(const std::string&, const std::vector<std::string>&);
    void insertValues(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&);
    void updateRow(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&, const std::string&);
    std::string prepareSELECT(const std::string&, const std::string&, const std::string&, const std::string& where = "", int limit = 0, const std::string& order_by = "");
    std::string prepareSELECT(const std::string&, const std::string&, const std::vector<std::string>&, const std::string& where = "", int limit = 0, const std::string& order_by = "");
  };

  std::vector<std::string> getSchemaKeys();

} // namespace aflowlib

#endif // AFLOWLIB_DATABASE_H
