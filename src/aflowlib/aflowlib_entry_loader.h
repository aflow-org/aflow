
#ifndef AFLOWLIB_ENTRY_LOADER_H
#define AFLOWLIB_ENTRY_LOADER_H

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"

#include "aflowlib/aflowlib_database.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_xclasses.h"

namespace aflowlib {
  class AflowDB;

  // HE20220420 START
  class EntryLoader : public xStream {
  public:
      // NECESSARY PUBLIC CLASS METHODS - START
      // constructors - START
    EntryLoader(std::ostream& oss = std::cout);
    EntryLoader(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    EntryLoader(const EntryLoader& b);
      // constructors - STOP
    ~EntryLoader();
    const EntryLoader& operator=(const EntryLoader& other);
    void clear();
      // NECESSARY PUBLIC CLASS METHODS - STOP

      /// @brief available sources to load lib entries
    enum class Source {
      SQLITE,         ///< internal AFLUX database
      AFLUX,          ///< public AFLUX API
      FILESYSTEM,     ///< search by public alloy DB, loaded from internal filesystem
      FILESYSTEM_RAW, ///< internal filesystem
      RESTAPI,        ///< search by public alloy DB, loaded from public REST API
      RESTAPI_RAW,    ///< public REST API
      NONE,           ///< no source set
      FAILED          ///< setting a source failed
    };

      // Settings
    bool m_out_silent; ///< silents all output
    bool m_out_debug;  ///< verbose output (overwrites #m_out_silent)

    bool m_xstructure_relaxed;  ///< add the relaxed structure into the lib entry
    bool m_xstructure_original; ///< add the original structure into the lib entry

      /// sources for the relaxed structure (in order of priority)
    std::vector<std::string> m_xstructure_final_file_name;

    std::string m_sqlite_file;          ///< location of the internal AFLUX SQLITE DB file
    std::string m_sqlite_alloy_file;    ///< location of the public alloy SQLITE DB file
    std::string m_sqlite_collection;    ///< collection to use for queries

    std::string m_aflux_server;         ///< server that provides the AFLUX API
    std::string m_aflux_path;           ///< base path including AFLUX API version
    std::string m_aflux_collection;     ///< collection to use for queries
      /// default AFLUX API directives to set format and disable paging
    std::map<std::string, std::string> m_aflux_directives;

    std::string m_restapi_server;        ///< server that provides the AFLOW RESTAPI
    std::string m_restapi_path;          ///< AFLOW RESTAPI base path
    std::string m_restapi_directives;    ///< directive to receive an entry (typically "?format=text" or "aflowlib.out")
    std::string m_restapi_listing_dirs;  ///< directive to list sub entries
    std::string m_restapi_listing_files; ///< directive to list files available for an entry
    std::string m_restapi_collection;    ///< collection to use for queries

    std::string m_filesystem_outfile;    ///< file name of the AFLOW lib entry
    std::string m_filesystem_path;       ///< base path to the internal filesystem
    std::string m_filesystem_collection; ///< collection to use for queries

      // Attributes
    std::unique_ptr<aflowlib::AflowDB> m_sqlite_db_ptr;       ///< pointer to an instance of the internal AFLUX SQLITE DB
    std::unique_ptr<aflowlib::AflowDB> m_sqlite_alloy_db_ptr; ///< pointer to an instance of the public alloy SQLITE DB
      /// @brief set that contains all loaded AUIDs
      /// @note `std::set` is stored sorted and is therefore faster at finding entries compared to `std::vector`
    std::vector<std::string> m_auid_list;

      // Data views
      /// @brief flat collection of loaded lib entries
    std::shared_ptr<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>> m_entries_flat;
      /// @brief tiered collection of loaded lib entries (number of elements -> alloy -> entries)
      /// @note creating a copy of this smart pointer is the most efficient way to use
      ///       the loaded entries after the EntryLoader class goes out of scope
    std::shared_ptr<std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>> m_entries_layered_map;

      // Generic entry loaders
    void loadAUID(std::string AUID);
    void loadAUID(const std::vector<std::string>& AUID);

    void loadAURL(std::string AURL);
    void loadAURL(const std::vector<std::string>& AURL);

    void loadAlloy(const std::string& alloy, bool recursive = true);
    void loadAlloy(const std::vector<std::string>& alloy, bool recursive = true);

      // Direct entry loader
    void loadSqliteWhere(const std::string& where);
    void loadAFLUXQuery(const std::string& query);
    void loadAFLUXMatchbook(const std::map<std::string, std::string>& matchbook);
    void loadRestAPIQueries(const std::vector<std::string>& queries, bool full_url = false);
    void loadFiles(const std::vector<std::string>& files);
    void loadFiles(const std::string& file_path);
    void loadFolders(const std::vector<std::string>& folders);
    void loadFolders(const std::string& folder_path);
    void loadText(const std::vector<std::string>& raw_data_lines, const std::string& source = "");
    void loadVector(const std::vector<std::string>& keys, const std::vector<std::vector<std::string>>& content, const std::string& source = "");

      // Source setter and getter
    bool setSource(EntryLoader::Source new_source);
    EntryLoader::Source getSource() const;
    std::string getSourceString() const;

      // Getter for raw data
    std::vector<std::string> getRawSqliteWhere(const std::string& where) const;
    std::vector<std::string> getRawAFLUXMatchbook(const std::map<std::string, std::string>& matchbook);
    std::vector<std::string> getRawAFLUXQuery(const std::string& query);
    std::string getRawRestAPIQuery(const std::string& query, bool full_url = false);

      // xstructure loaders
    void addXstructure(aflowlib::_aflowlib_entry& entry, bool orig = false);
    bool loadXstructureFile(const aflowlib::_aflowlib_entry& entry, xstructure& new_structure, std::vector<std::string> possible_files = {});
    bool loadXstructureAflowIn(const aflowlib::_aflowlib_entry& entry, xstructure& new_structure, const int index = -1);

      // getter for views on the loaded aflowlib entries (using zero copy shared_ptr)
    void getEntriesViewFlat(std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>& result) const;
    void getEntriesViewTwoLayer(std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>& result) const;
    void getEntriesViewThreeLayer(std::vector<std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>& result) const;
    void getEntriesViewMap(std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>& result) const;

      // getter that copies the loaded aflowlib entries into new structure
    void getEntriesFlat(std::vector<aflowlib::_aflowlib_entry>& result) const;
    void getEntriesTwoLayer(std::vector<std::vector<aflowlib::_aflowlib_entry>>& result) const;
    void getEntriesThreeLayer(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& result) const;

    static bool cleanAUID(std::string& AUID);
    static bool cleanAURL(std::string& AURL);

  private:
      // NECESSARY private CLASS METHODS - START
    void init();
    void free();
    void copy(const EntryLoader& b);

      // NECESSARY END CLASS METHODS - END

      /// @brief Determines which source is used to load data
      /// @note should just be set by setSource()
    Source m_current_source;
    bool m_filesystem_available; ///< helper to avoid repeatedly testing if the filesystem is available
    std::stringstream m_logger_message;  ///< reusable stringstream for logging
    bool m_out_super_silent;     ///< mute all output - used for private checks

    void selectSource();
    bool getCommandlineSource();
    void listRestAPI(std::string url, std::vector<std::string>& result, bool directories = true);
    void getAlloyAUIDList(const std::vector<std::string>& alloy_list, std::vector<std::string>& auid_list);
    void loadAlloySearchFSR(const std::vector<std::string>& alloy_list, uint lib_max); // FILESYSTEM_RAW
    void loadAlloySearchRR(const std::vector<std::string>& alloy_list, uint lib_max); // RESTAPI_RAW
    std::string buildAFLUXQuery(const std::map<std::string, std::string>& matchbook) const;
    std::string extractAlloy(std::string name, char lib_type) const;

      // Logging helper
    void outInfo(const std::string& function_name);
    void outDebug(const std::string& function_name);
    void outError(const std::string& function_name, int line_number);
    void outHardError(const std::string& function_name, int line_number, int error_type);
  };
  // HE20220420 END
} // namespace aflowlib

#endif // AFLOWLIB_ENTRY_LOADER_H
