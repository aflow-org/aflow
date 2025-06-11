// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// This EntryLoader class is the combination/evolution of different prior solutions to integrate with the AFLOW database.
// 2016 - Corey Oses (CHULL)
// 2018 - Marco Esters (database design/integration)
// 2018 - Frisco Rose (API requests)
// 2020 - David Hicks (structural comparisons)
// 2022 - Hagen Eckert (unified EntryLoader for all available database sources)
// hagen.eckert@duke.edu


#ifndef _AFLOWLIB_ENTRY_LOADER_CPP_
#define _AFLOWLIB_ENTRY_LOADER_CPP_

#include "aflow.h"
#include "aflowlib.h"

#define _DEBUG_ENTRY_LOADER_ false


namespace aflowlib {

  /// @class EntryLoader
  /// @brief unified interface to gather AFLOW lib entries from different sources
  ///
  /// @authors
  /// @mod{CO,2016,CHULL}
  /// @mod{ME,2018,database design/integration}
  /// @mod{FR,2018,API requests}
  /// @mod{DX,2020,structural comparisons}
  /// @mod{HE,20220819,created EntryLoader}
  ///
  /// Basic usage
  /// @code
  /// std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> results;
  /// { // ensures that all resources are released as soon as possible
  ///   aflowlib::EntryLoader el;
  ///   el.m_out_debug = true; // change the behavior of the EntryLoader class
  ///   el.loadAlloy("NiCa"); // the load functions can be called multiple times
  ///   el.loadAUID("aflow:7dd846bc04c764e8"); // no duplicates will be stored
  ///   el.getEntriesViewFlat(results);
  /// }
  /// for (std::shared_ptr<aflowlib::_aflowlib_entry>> entry : results) std::cout << entry->auid << std::endl;
  /// @endcode

  // class constructor
  EntryLoader::EntryLoader(std::ostream &oss) : xStream(oss) { init(); }

  EntryLoader::EntryLoader(std::ofstream &FileMESSAGE, std::ostream &oss) : xStream(FileMESSAGE, oss) { init(); }

  /// @brief create a copy at initialization
  /// @note DB connections are reused and not recreated
  /// @note a new view on the data is created, but the same lib entries are reused
  EntryLoader::EntryLoader(const EntryLoader &b) : xStream(*b.getOFStream(), *b.getOSS()) { copy(b); } // copy PUBLIC

  // class de-constructor
  EntryLoader::~EntryLoader() { xStream::free(); }

  /// @brief create a copy through assigning
  /// @note DB connections are reused and not recreated
  /// @note a new view on the data is created, but the same lib entries are reused
  const EntryLoader &EntryLoader::operator=(const EntryLoader &other) {
    if (this != &other) { copy(other); }
    return *this;
  }

  /// @brief initialize the class (private)
  /// create shared pointer used to store the data views and read default values from flags
  void EntryLoader::init() {
    m_out_silent = true;
    m_out_debug = (XHOST.DEBUG || _DEBUG_ENTRY_LOADER_);
    m_xstructure_relaxed = false;
    m_xstructure_original = false;
    m_xstructure_final_file_name = {"CONTCAR.relax2", "CONTCAR.relax", "POSCAR.static", "POSCAR.bands", "CONTCAR.static", "CONTCAR.bands"};

    m_sqlite_file = DEFAULT_AFLOW_DB_FILE;
    m_sqlite_alloy_file = DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE;
    m_sqlite_collection = "RAW";

    m_aflux_server = DEFAULT_ENTRY_LOADER_AFLUX_SERVER;
    m_aflux_path = DEFAULT_ENTRY_LOADER_AFLUX_PATH;
    m_aflux_collection = "RAW";
    m_aflux_directives = {{"format", "aflow"}, {"paging", "0"}};

    m_restapi_server = DEFAULT_ENTRY_LOADER_RESTAPI_SERVER;
    m_restapi_path = DEFAULT_ENTRY_LOADER_RESTAPI_PATH;
    m_restapi_directives = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;
    m_restapi_listing_dirs = "?aflowlib_entries";
    m_restapi_listing_files = "?files";
    m_restapi_collection = "WEB";
    m_filesystem_outfile = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;
    m_filesystem_path = DEFAULT_ENTRY_LOADER_FS_PATH;
    m_filesystem_collection = "RAW";

    m_entries_flat = std::make_shared<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>();
    m_entries_layered_map = std::make_shared<std::map<short,
                                             std::map<std::string,
                                             std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>>();

    m_current_source = Source::NONE;
    m_filesystem_available = false;
    m_out_super_silent = false;
  }

  /// @brief clean a EntryLoader object
  void EntryLoader::clear() { free(); }

  /// @brief reset a EntryLoader object
  void EntryLoader::free() { *this = {}; } // calling the default constructor to ensure a clean instance of EntryLoader

  /// @brief create a copy (privat)
  void EntryLoader::copy(const EntryLoader &b) {  //copy PRIVATE
    m_out_debug = b.m_out_debug;
    m_out_silent = b.m_out_silent;
    m_xstructure_relaxed = b.m_xstructure_relaxed;
    m_xstructure_original = b.m_xstructure_original;
    m_xstructure_final_file_name = b.m_xstructure_final_file_name;
    m_sqlite_file = b.m_sqlite_file;
    m_sqlite_alloy_file = b.m_sqlite_alloy_file;
    m_sqlite_collection = b.m_sqlite_collection;
    m_aflux_server = b.m_aflux_server;
    m_aflux_path = b.m_aflux_path;
    m_aflux_collection = b.m_aflux_collection;
    m_aflux_directives = b.m_aflux_directives;
    m_restapi_server = b.m_restapi_server;
    m_restapi_path = b.m_restapi_path;
    m_restapi_directives = b.m_restapi_directives;
    m_restapi_listing_dirs = b.m_restapi_listing_dirs;
    m_restapi_listing_files = b.m_restapi_listing_files;
    m_restapi_collection = b.m_restapi_collection;
    m_filesystem_outfile = b.m_filesystem_outfile;
    m_filesystem_path = b.m_filesystem_path;
    m_filesystem_collection = b.m_filesystem_collection;
    m_current_source = b.m_current_source;
    m_filesystem_available = b.m_filesystem_available;
    m_auid_list = b.m_auid_list;
    aurostd::StringstreamClean(m_logger_message);

    m_sqlite_db_ptr = nullptr;
    m_sqlite_alloy_db_ptr = nullptr;

    m_entries_flat = std::make_shared < std::vector < std::shared_ptr < aflowlib::_aflowlib_entry>>>(*b.m_entries_flat);
    m_entries_layered_map = std::make_shared < std::map < short,
        std::map < std::string,
        std::vector < std::shared_ptr < aflowlib::_aflowlib_entry >>>>>(*b.m_entries_layered_map);
  }

  /// @brief load a single entry from a AUID
  /// @param AUID AFLOW unique ID
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note following AUIDs would be equivalent:
  /// @note `aflow:d912e209c81aeb94`
  /// @note `auid:d912e209c81aeb94`
  /// @note `d912e209c81aeb94`
  void EntryLoader::loadAUID(std::string AUID) {
    selectSource();

    m_logger_message << "Try loading: " << AUID;
    outInfo(__AFLOW_FUNC__);
    if (!cleanAUID(AUID)) {
      m_logger_message << "AUID cleaning failed";
      outError(__AFLOW_FUNC__, __LINE__);
      return;
    }

    switch (m_current_source) {

      case Source::SQLITE: {
        std::string where = "auid='\"" + AUID + "\"'";
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        std::map <std::string, std::string> matchbook{{"*",    ""},
                                                      {"auid", "'" + AUID + "'"}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        loadAUID(std::vector<std::string>({AUID}));
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        loadAUID(std::vector<std::string>({AUID}));
        break;
      }

      default:
        break;
    }
  }

  /// @brief load multiple entries from a vector of AUIDs
  /// @param AUID list of AURLs
  /// @authors
  /// @mod{HE,20220216,created}
  /// @mod{HE,20221117,ensure that aflux URL is not longer than 8000 characters}
  void EntryLoader::loadAUID(const std::vector<std::string> &AUID) {
    selectSource();
    std::vector <std::string> clean_AUID;

    m_logger_message << "Try loading " << AUID.size() <<" AUIDs";
    outInfo(__AFLOW_FUNC__);

    for (std::vector<std::string>::const_iterator AUID_single = AUID.begin(); AUID_single != AUID.end(); AUID_single++) {
      std::string temp = *AUID_single;
      if (cleanAUID(temp)) clean_AUID.push_back(temp);
    }

    if (clean_AUID.size()!=AUID.size()) {
      m_logger_message << "cleaning of " << AUID.size() - clean_AUID.size()  << " AUIDs failed";
      outError(__AFLOW_FUNC__, __LINE__);
    }

    switch (m_current_source) {

      case Source::AFLUX: {
        std::string AUID_combined = "";
        bool first = true;
        for (size_t i=0; i<clean_AUID.size(); i++){
          if ((AUID_combined+ clean_AUID[i]).size()  < 5000) {
            if (first) {
              AUID_combined=clean_AUID[i];
              first=false;
            }
            else AUID_combined+="':'" + clean_AUID[i];
          }
          else {
            AUID_combined = "'" + AUID_combined + "'";
            loadAFLUXMatchbook({{"*",    ""},
                                {"auid", AUID_combined}});
            AUID_combined = clean_AUID[i];
          }
        }

        AUID_combined = "'" + AUID_combined + "'";
        loadAFLUXMatchbook({{"*",    ""},
                            {"auid", AUID_combined}});
        break;
      }

      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(clean_AUID, "\"','\"");
        where = "auid IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        std::vector <std::string> queries;
        for (std::vector<std::string>::const_iterator AUID_single = clean_AUID.begin(); AUID_single != clean_AUID.end(); AUID_single++) {
          std::string rest_query = "AUID/aflow:";
          for (uint part = 6; part <= 20; part+= 2) rest_query += AUID_single->substr(part, 2) + "/";
          rest_query += "WEB/";
          queries.push_back(rest_query);
        }
        loadRestAPIQueries(queries);
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        std::vector <std::string> files;
        for (std::vector<std::string>::const_iterator AUID_single = clean_AUID.begin(); AUID_single != clean_AUID.end(); AUID_single++) {
          std::string file_path = m_filesystem_path + "AUID/aflow:";
          for (uint part = 6; part <= 20; part+= 2) file_path += AUID_single->substr(part, 2) + "/";
          file_path += m_filesystem_collection + "/" + m_filesystem_outfile;
          files.push_back(file_path);
        }
        loadFiles(files);
        break;
      }

      default:
        break;
    }
  }

  /// @brief load a single entry from an AURL
  /// @param AURL AFLOW unique URL
  /// @authors
  /// @mod{HE,20220217,created}
  /// @note following AURLs would be equivalent:
  /// @note `aflowlib.duke.edu:AFLOWDATA/LIB2_RAW/Ca_svCu_pv/138`
  /// @note `AFLOWDATA/LIB2_RAW/Ca_svCu_pv/138`
  /// @note `LIB2_RAW/Ca_svCu_pv/138`
  /// @note `WEB`, `RAW`, and `LIB` will be replaced by the selected collection
  ///       (#m_sqlite_collection, #m_aflux_collection, #m_restapi_collection, #m_filesystem_collection)
  void EntryLoader::loadAURL(std::string AURL) {
    selectSource();

    if (!cleanAURL(AURL)) return;

    m_logger_message << "Try loading " << AURL;
    outInfo(__AFLOW_FUNC__);
    if (!cleanAURL(AURL)) {
      m_logger_message << "AURL cleaning failed";
      outError(__AFLOW_FUNC__, __LINE__);
      return;
    }

    switch (m_current_source) {

      case Source::SQLITE: {
        std::string where = "aurl='\"" + AURL + "\"'";
        where = std::regex_replace(where, aurostd::regex_aurl2file, "$1_" + m_sqlite_collection + "/");
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        AURL = std::regex_replace(AURL, aurostd::regex_aurl2file, "$1_" + m_aflux_collection + "/");
        std::map <std::string, std::string> matchbook{{"*",    ""},
                                                      {"aurl", "'" + AURL + "'"}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        loadAURL(std::vector<std::string>({AURL}));
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        loadAURL(std::vector<std::string>({AURL}));
        break;
      }

      default:
        break;
    }
  }

  /// @brief load multiple entries from a vector of AURLs
  /// @param AURL list of AURLs
  /// @authors
  /// @mod{HE,20220322,created}
  /// @mod{HE,20221117,ensure that aflux URL is not longer than 8000 characters}
  void EntryLoader::loadAURL(const std::vector <std::string> &AURL) {
    selectSource();

    m_logger_message << "Try loading " << AURL.size() <<" AURLs";
    outInfo(__AFLOW_FUNC__);

    std::vector <std::string> clean_AURL;
    for (std::vector<std::string>::const_iterator AURL_single = AURL.begin(); AURL_single != AURL.end(); AURL_single++) {
      std::string temp = *AURL_single;
      if (cleanAURL(temp)) clean_AURL.push_back(temp);
    }

    if (clean_AURL.size()!=AURL.size()) {
      m_logger_message << "cleaning of " << AURL.size() - clean_AURL.size()  << " AURLs failed";
      outError(__AFLOW_FUNC__, __LINE__);
    }

    switch (m_current_source) {

      case Source::AFLUX: {
        std::string AURL_combined = "";
        bool first = true;
        for (size_t i=0; i<clean_AURL.size(); i++){
          if ((AURL_combined+ clean_AURL[i]).size()  < 5000) {
            if (first) {
              AURL_combined=clean_AURL[i];
              first=false;
            }
            else AURL_combined+="':'" + clean_AURL[i];
          }
          else {
            AURL_combined = "'" + AURL_combined + "'";
            AURL_combined = std::regex_replace(AURL_combined, aurostd::regex_aurl2file, "$1_" + m_aflux_collection + "/");
            loadAFLUXMatchbook({{"*",    ""},
                                {"aurl", AURL_combined}});
            AURL_combined = clean_AURL[i];
          }
        }

        AURL_combined = "'" + AURL_combined + "'";
        AURL_combined = std::regex_replace(AURL_combined, aurostd::regex_aurl2file, "$1_" + m_aflux_collection + "/");
        loadAFLUXMatchbook({{"*",    ""},
                            {"aurl", AURL_combined}});
        break;
      }

      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(clean_AURL, "\"','\"");
        where = std::regex_replace(where, aurostd::regex_aurl2file, "$1_" + m_sqlite_collection + "/");
        where = "aurl IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        std::vector <std::string> queries;
        for (std::vector<std::string>::const_iterator AURL_single = clean_AURL.begin(); AURL_single != clean_AURL.end(); AURL_single++) {
          std::string rest_query = std::regex_replace(AURL_single->substr(28), aurostd::regex_aurl2file, "$1_" + m_sqlite_collection + "/") + "/";
          queries.push_back(rest_query);
        }
        loadRestAPIQueries(queries);
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        std::vector <std::string> files;
        for (std::vector<std::string>::const_iterator AURL_single = clean_AURL.begin(); AURL_single != clean_AURL.end(); AURL_single++) {
          std::string file_path = m_filesystem_path + AURL_single->substr(28) + "/" + m_filesystem_outfile;
          file_path = std::regex_replace(file_path, aurostd::regex_aurl2file, "$1/" + m_filesystem_collection + "/");
          files.push_back(file_path);
        }
        loadFiles(files);
        break;
      }

      default:
        break;
    }
  }

  /// @brief load all entries that are available for a alloy set
  /// @param alloy alloy descriptor
  /// @param recursive include sub alloys
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note an alloy can be given as comma separated list `"Cu, Au, Fe"`
  ///       or as combined string `"CuAuFe"` (elements don't need to be sorted)
  /// @note `recursive=true` (default) will generate sub-alloys so `CuAuFe`
  ///       will load `AuCuFe`, `AuCu`, `AuFe`, `CuFe`, `Au`, `Cu`, and `Fe`
  void EntryLoader::loadAlloy(const std::string &alloy, bool recursive) {
    std::vector <std::string> alloy_elements;
    if (alloy.find(", ") != std::string::npos) aurostd::string2tokens(alloy, alloy_elements, ", ");
    else if (alloy.find(',') != std::string::npos) aurostd::string2tokens(alloy, alloy_elements, ",");
    else alloy_elements = aurostd::getElements(alloy);
    loadAlloy(alloy_elements, recursive);
  }

  /// @brief load all entries that are available for a set alloy
  /// @param alloy split alloy elements
  /// @param recursive include sub alloys
  /// @authors
  /// @mod{HE,20220216,created}
  void EntryLoader::loadAlloy(const std::vector <std::string> &alloy, bool recursive) {
    selectSource();
    std::vector <std::string> alloy_clean;

    // check that all parts of alloy are actually elements
    for (std::vector<std::string>::const_iterator element = alloy.begin(); element != alloy.end(); element++) {
      if (xelement::xelement::isElement(*element) == 0) {
        m_logger_message << *element << " is not an element";
        outError(__AFLOW_FUNC__, __LINE__);
      } else alloy_clean.emplace_back(*element);
    }

    // build list of all needed sub-alloys
    std::vector<std::vector<std::string>> alloy_sub_list;
    if (recursive) {
      aurostd::xcombos xc;
      for (size_t i = 1; i <= alloy_clean.size(); i++) {
        xc.reset(alloy_clean.size(), i);
        std::vector<int> choice;
        while (xc.increment()) {
          std::vector <std::string> alloy_tmp;
          choice = xc.getCombo();
          for (size_t idx = 0; idx < choice.size(); ++idx) {
            if (choice[idx]) alloy_tmp.push_back(alloy_clean[idx]);
          }
          alloy_sub_list.push_back(alloy_tmp);
        }
      }
    } else {
      alloy_sub_list.push_back(alloy_clean);
    }

    // build sorted alloy strings
    std::vector <std::string> final_alloy_list;
    for (std::vector<std::vector<std::string>>::iterator a = alloy_sub_list.begin(); a != alloy_sub_list.end(); a++) {
      std::string alloy_tmp = "";
      std::sort(a->begin(), a->end());
      for (std::vector<std::string>::iterator e = a->begin(); e != a->end(); e++) alloy_tmp += *e;
      final_alloy_list.push_back(alloy_tmp);
    }
    m_logger_message << "Try loading " << final_alloy_list.size() << " systems: ";
    for (std::vector<std::string>::const_iterator i_alloy = final_alloy_list.begin(); i_alloy != final_alloy_list.end(); i_alloy++) {
      m_logger_message << "\"" << *i_alloy << "\" ";
    }
    outInfo(__AFLOW_FUNC__);

    // load the alloys
    switch (m_current_source) {
      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(final_alloy_list, "','");
        where = "alloy IN ('" + where + "')";
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        std::string alloy_match = aurostd::joinWDelimiter(final_alloy_list, ":");
        std::map <std::string, std::string> matchbook{{"*",     ""},
                                                      {"alloy", alloy_match}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI: {
        std::vector <std::string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }

      case Source::RESTAPI_RAW: {
        loadAlloySearchRR(final_alloy_list, alloy_clean.size());
        break;
      }

      case Source::FILESYSTEM: {
        std::vector <std::string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }

      case Source::FILESYSTEM_RAW: {
        loadAlloySearchFSR(final_alloy_list, alloy_clean.size());
        break;
      }

      default:
        break;
    }

  }

  /// @brief load entries based on a custom query against the internal AFLUX SQLITE DB
  /// @param where part of the SQLITE query after a `WHERE` keyword
  /// @authors
  /// @mod{HE,20220422,created}
  /// @note use setSource() to change the current source to Source::SQLITE to ensure that the database is connected
  void EntryLoader::loadSqliteWhere(const std::string &where) {
//    std::vector <std::string> raw_lines = getRawSqliteWhere(where);
    vector<vector<string>> content = m_sqlite_db_ptr->getRowsMultiTables(where);
    vector<string> keys = m_sqlite_db_ptr->getColumnNames("auid_00");

    m_logger_message << content.size() << " entries found in DB for " << where;
    outDebug(__AFLOW_FUNC__);

    size_t start_size = m_entries_flat->size();
    loadVector(keys, content, where);

    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__AFLOW_FUNC__);

  }

  /// @brief load entries from a costum AFLUX query
  /// @param query AFLUX query
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note #m_aflux_server and #m_aflux_path will be added
  void EntryLoader::loadAFLUXQuery(const std::string &query) {
    size_t start_size = m_entries_flat->size();
    loadText(getRawAFLUXQuery(query), query);
    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__AFLOW_FUNC__);
  }

  /// @brief load entries from an AFLUX matchbook
  /// @param matchbook map with keyword and modifier
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note a pair creates `<keyword>(<modifier>)`
  /// @note if modifier is empty just the keyword is used
  void EntryLoader::loadAFLUXMatchbook(const std::map <std::string, std::string> &matchbook) {
    loadAFLUXQuery(buildAFLUXQuery(matchbook));
  }

  /// @brief load entries from a list of REST API queries
  /// @authors
  /// @mod{HE,20220216,created}
  /// @param queries vector of queries
  /// @param full_url if `true` don't add #m_restapi_server, #m_restapi_path, and #m_restapi_directives
  void EntryLoader::loadRestAPIQueries(const std::vector <std::string> &queries, bool full_url) {
    size_t start_size = m_entries_flat->size();
    size_t done_downloads = 0;
    for (std::vector<std::string>::const_iterator query = queries.begin(); query != queries.end(); query++) {
      loadText({getRawRestAPIQuery(*query, full_url)}, *query);
      done_downloads++;
      if (done_downloads%100 == 0){
        m_logger_message << "Loaded " << done_downloads << " of " << queries.size();
        outInfo(__AFLOW_FUNC__);
      }
    }
    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__AFLOW_FUNC__);
  }

  /// @brief load entries from a list of file paths
  /// @param files list of file paths
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note doesn't add #m_filesystem_path or #m_filesystem_collection
  void EntryLoader::loadFiles(const std::vector <std::string> &files) {
    size_t start_size = m_entries_flat->size();
    Source last_source = m_current_source;
    m_current_source = Source::FILESYSTEM_RAW;
    for (std::vector<std::string>::const_iterator file_path = files.begin(); file_path != files.end(); file_path++) {
       std::string file_content;
       if (aurostd::file2string(*file_path, file_content) > 0) {
         loadText({file_content}, *file_path);
       }
    }
    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__AFLOW_FUNC__);
    m_current_source = last_source;
  }

  /// @brief load a entry from a file path
  /// @param file_path file paths
  /// @authors
  /// @mod{HE,20220913,created}
  /// @note doesn't add #m_filesystem_path or #m_filesystem_collection
  void EntryLoader::loadFiles(const std::string &file_path) {
    loadFiles((vector<std::string>){file_path});
  }

  /// @brief load all entries in the given folders and their sub-folders
  /// @param folders list of folder paths
  /// @authors
  /// @mod{HE,20220913,created}
  /// @note EntryLoader::m_filesystem_outfile (default aflowlib.out) is used to identify entries
  void EntryLoader::loadFolders(const std::vector<std::string> &folders) {
    std::vector <std::string> found_entries;

    char **paths = new char*[folders.size()+1](); // initialize the array with nullptr, +1 ensures it is null terminated
    for (size_t list_idx = 0; list_idx < folders.size(); list_idx++) {
      paths[list_idx] = (char *) folders[list_idx].c_str();
    }

    size_t scanned = 0;
    size_t found = 0;
    struct stat file_stats{};

    // initialize a file tree scan
    // https://man7.org/linux/man-pages/man3/fts.3.html
    // FTS_PHYSICAL - don't follow symlinks
    // FTS_NOCHDIR - don't change the workdir of the program
    // FTS_XDEV - don't descend into folders that are on a different device
    FTS *tree = fts_open(paths, FTS_PHYSICAL | FTS_NOCHDIR | FTS_XDEV, nullptr);
    FTSENT *node;
    if (tree == nullptr) {
      m_logger_message << "Failed to initialize the file tree used to search for alloys!";
      outHardError(__AFLOW_FUNC__, __LINE__,_FILE_ERROR_);
      return;
    }

    // Iterate over all entries found in the folders listed in paths
    while ((node = fts_read(tree))) {
      scanned += 1;
      if (m_out_debug && ((scanned-1) % 1000 == 0)) {
        m_logger_message << (scanned-1) << " objects scanned; " << found << " entries found; next scan: " << node->fts_path;
        outDebug(__AFLOW_FUNC__);
      }
      if ((node->fts_info & FTS_D)) {
        std::string base_path = node->fts_path;
        std::string full_path = base_path + "/" + m_filesystem_outfile;
        if (stat(full_path.c_str(), &file_stats) == 0) {
          found += 1;
          found_entries.emplace_back(full_path);
        }
      }
    }
    fts_close(tree);
    delete[] paths;

    m_logger_message << "Finishing search in the filesystem after scanning " << scanned << " objects";
    outInfo(__AFLOW_FUNC__);
    Source last_source = m_current_source;
    m_current_source = Source::FILESYSTEM_RAW;
    loadFiles(found_entries);
    m_current_source = last_source;
  }

  /// @brief load all entries recursive in the given folder and its sub-folders
  /// @param folder_path folder path
  /// @authors
  /// @mod{HE,20220913,created}
  /// @note EntryLoader::m_filesystem_outfile (default aflowlib.out) is used to identify entries
  void EntryLoader::loadFolders(const std::string &folder_path) {
    loadFolders((vector<std::string>){folder_path});
  }

  /// @brief load entries from a list of strings
  /// @param raw_data_lines list of data strings
  /// @param source string describing the source, like a SQL query or a path
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note each entry in raw_data_lines should correspond to one AFLOW lib entry
  void EntryLoader::loadText(const std::vector <std::string> &raw_data_lines, const std::string &source) {
    for (std::vector<std::string>::const_iterator line = raw_data_lines.begin(); line != raw_data_lines.end(); line++) {
      std::shared_ptr <aflowlib::_aflowlib_entry> entry = std::make_shared<aflowlib::_aflowlib_entry>();
      entry->Load(*line, *p_oss);
      entry->el_source_type = getSourceString();
      entry->el_source = source;
      if (!entry->auid.empty() && (std::find(m_auid_list.begin(),m_auid_list.end(), entry->auid) == m_auid_list.end())) {
        m_entries_flat->push_back(entry);
        (*m_entries_layered_map)[entry->nspecies][entry->species].push_back(entry);
        m_auid_list.emplace_back(entry->auid);
        if (m_xstructure_original) addXstructure(*entry, true);
        if (m_xstructure_relaxed) addXstructure(*entry);
      }
    }
  }

  /// @brief load entries from vectors
  /// @param keys list of keys
  /// @param content list of keys
  /// @param source string describing the source, like a SQL query or a path
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note each entry in content should correspond to one AFLOW lib entry
  void EntryLoader::loadVector(const std::vector<std::string> &keys, const std::vector<std::vector<std::string>> & content, const std::string & source) {
    std::vector<uint64_t> hash_list;
    for (std::vector<std::string>::const_iterator key = keys.begin(); key != keys.end(); key++) hash_list.emplace_back(aurostd::crc64(*key));
    for (std::vector<std::vector<std::string>>::const_iterator row=content.begin(); row!=content.end(); row++) {
      std::shared_ptr <aflowlib::_aflowlib_entry> entry = std::make_shared<aflowlib::_aflowlib_entry>();
      entry->Load(hash_list, *row);
      entry->el_source_type = getSourceString();
      entry->el_source = source;
      if (!entry->auid.empty() && (std::find(m_auid_list.begin(),m_auid_list.end(), entry->auid) == m_auid_list.end())) {
        m_entries_flat->push_back(entry);
        (*m_entries_layered_map)[entry->nspecies][entry->species].push_back(entry);
        m_auid_list.emplace_back(entry->auid);
        if (m_xstructure_original) addXstructure(*entry, true);
        if (m_xstructure_relaxed) addXstructure(*entry);
      }
    }
  }

  /// @brief returns the active data source
  EntryLoader::Source EntryLoader::getSource() const{
    return m_current_source;
  }

  /// @brief returns the active data source as string
  std::string EntryLoader::getSourceString() const {
    switch(m_current_source) {
      case Source::SQLITE: return "SQLITE";
      case Source::AFLUX: return "AFLUX";
      case Source::FILESYSTEM: return "FILESYSTEM";
      case Source::FILESYSTEM_RAW: return "FILESYSTEM_RAW";
      case Source::RESTAPI: return "RESTAPI";
      case Source::RESTAPI_RAW: return "RESTAPI_RAW";
      case Source::NONE: return "NONE";
      case Source::FAILED: return "FAILED";
      default: return "NONE";
    }
  }

  /// @brief change the currently used source and prepares them
  /// @param new_source source to change to
  /// @param silent silent all logging when called from selectSource() (default false)
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note if the public alloy SQLITE DB is not found Source::FILESYSTEM and Source::RESTAPI
  ///       are mapped to Source::FILESYSTEM_RAW and Source::RESTAPI_RAW respectively
  bool EntryLoader::setSource(EntryLoader::Source new_source) {
    if (new_source == m_current_source) return true;
    m_current_source = Source::FAILED;
    m_filesystem_available = aurostd::IsDirectory(m_filesystem_path + "AUID/");

    switch (new_source) {
      case Source::SQLITE: {
        if (aurostd::FileExist(m_sqlite_file)) {
          if (m_sqlite_db_ptr == nullptr) {
            std::unique_ptr<aflowlib::AflowDB> tmpPtr(new aflowlib::AflowDB(m_sqlite_file));
            m_sqlite_db_ptr = std::move(tmpPtr);
          }
          m_current_source = Source::SQLITE;
          return true;
        } else {
          m_logger_message << "Internal AFLUX SQLITE DB not found at " << m_sqlite_file;
          outError(__AFLOW_FUNC__, __LINE__);
        }
        break;
      }

      case Source::AFLUX: {
        if ("AFLUXtest" == aurostd::httpGet(m_aflux_server+"/test/?echo=AFLUXtest")) {
          m_current_source = Source::AFLUX;
          return true;
        } else {
          m_logger_message << "AFLUX API could not be reached at " << m_aflux_server;
          outError(__AFLOW_FUNC__, __LINE__);
        }
        break;
      }

      case Source::FILESYSTEM: {
        if (m_filesystem_available) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) { //
            if (m_sqlite_alloy_db_ptr == nullptr) {
              std::unique_ptr<aflowlib::AflowDB> tmpPtr(new aflowlib::AflowDB(m_sqlite_alloy_file));
              m_sqlite_alloy_db_ptr = std::move(tmpPtr);
            }
            m_current_source = Source::FILESYSTEM;
            return true;
          } else {
            m_logger_message << "Could not find public alloy SQLITE DB at " << m_sqlite_alloy_file << "; switching to Source::FILESYSTEM_RAW";
            outInfo(__AFLOW_FUNC__);
            return setSource(Source::FILESYSTEM_RAW);
          }
        } else {
          m_logger_message << "Could not find AFLOW entries in the filesystem at " << m_filesystem_path;
          outError(__AFLOW_FUNC__, __LINE__);
        }
        break;
      }

      case Source::FILESYSTEM_RAW: {
        if (m_filesystem_available) { //
          if (aurostd::FileExist(m_sqlite_alloy_file)) { // use the alloy DB to optimize the RAW results
            if (m_sqlite_alloy_db_ptr == nullptr) {
              std::unique_ptr<aflowlib::AflowDB> tmpPtr(new aflowlib::AflowDB(m_sqlite_alloy_file));
              m_sqlite_alloy_db_ptr = std::move(tmpPtr);
            }
          }
          m_current_source = Source::FILESYSTEM_RAW;
          return true;
        } else {
          m_logger_message << "Could not find AFLOW entries in the filesystem at " << m_filesystem_path;
          outError(__AFLOW_FUNC__, __LINE__);
        }
        break;
      }

      case Source::RESTAPI: {
        if (200 == aurostd::httpGetStatus(m_restapi_server+m_restapi_path+"ICSD_WEB/")) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) {
            if (m_sqlite_alloy_db_ptr == nullptr) {
              std::unique_ptr<aflowlib::AflowDB> tmpPtr(new aflowlib::AflowDB(m_sqlite_alloy_file));
              m_sqlite_alloy_db_ptr = std::move(tmpPtr);
            }
            m_current_source = Source::RESTAPI;
            return true;
          } else {
            m_logger_message << "Could not find public alloy SQLITE DB at " << m_sqlite_alloy_file << "; switching to Source::RESTAPI_RAW";
            outInfo(__AFLOW_FUNC__);
            return setSource(Source::RESTAPI_RAW);
          }
        } else {
          m_logger_message << "AFLOW REST API could not be reached at " << m_aflux_server+m_restapi_path;
          outError(__AFLOW_FUNC__, __LINE__);
        }
        break;
      }

      case Source::RESTAPI_RAW: {
        if (200 == aurostd::httpGetStatus(m_restapi_server+m_restapi_path+"ICSD_WEB/")) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) {
            if (m_sqlite_alloy_db_ptr == nullptr) {
              std::unique_ptr<aflowlib::AflowDB> tmpPtr(new aflowlib::AflowDB(m_sqlite_alloy_file));
              m_sqlite_alloy_db_ptr = std::move(tmpPtr);
            }
          }
          m_current_source = Source::RESTAPI_RAW;
          return true;
        }
        else {
          m_logger_message << "AFLOW REST API could not be reached at " << m_aflux_server+m_restapi_path;
          outError(__AFLOW_FUNC__, __LINE__);
        }
        break;
      }

      case Source::NONE: {
        m_current_source = Source::NONE;
        return true;
      }

      case Source::FAILED: {
        m_current_source = Source::FAILED;
        return true;
      }

    }
    return false;
  }

  /// @brief query the internal AFLUX SQLITE DB
  /// @param where part of the SQLITE query after a `WHERE` keyword
  /// @authors
  /// @mod{HE,20220216,created}
  /// @return list of result query strings
  std::vector <std::string> EntryLoader::getRawSqliteWhere(const std::string &where) const{
    if (m_sqlite_db_ptr != nullptr) return m_sqlite_db_ptr->getEntrySet(where, aflow_ft);
    else return {};
  }

  /// @brief query AFLUX using a matchbook
  /// @param matchbook map with keyword and modifier
  /// @authors
  /// @mod{HE,20220216,created}
  /// @return list of query result strings
  std::vector <std::string> EntryLoader::getRawAFLUXMatchbook(const std::map <std::string, std::string> &matchbook){
    return getRawAFLUXQuery(buildAFLUXQuery(matchbook));
  }

  /// @brief query AFLUX using a raw query
  /// @param matchbook map with keyword and modifier
  /// @authors
  /// @mod{HE,20220216,created}
  /// @return list of query result strings
  /// @note #m_aflux_server and #m_aflux_path will be added
  std::vector <std::string> EntryLoader::getRawAFLUXQuery(const std::string &query){
    std::string output = "";
    std::vector <std::string> raw_lines;
    if (200 == aurostd::httpGetStatus(m_aflux_server, m_aflux_path, query, output)){
      aurostd::string2vectorstring(output, raw_lines);
    } else {
      m_logger_message << "Failed to get AFLUX query ";
      m_logger_message << "(" << m_aflux_server << " | " << m_aflux_path << " | " << query << ")";
      outError(__AFLOW_FUNC__, __LINE__);
    }
    return raw_lines;
  }

  /// @brief query the REST API
  /// @param query REST API query
  /// @param full_url if `true` don't add #m_restapi_server, #m_restapi_path, and #m_restapi_directives
  /// @authors
  /// @mod{HE,20220216,created}
  /// @return query result string
  std::string EntryLoader::getRawRestAPIQuery(const std::string &query, bool full_url) {
    std::string output = "";
    std::string url;
    if (full_url) url = query;
    else url = m_restapi_server + m_restapi_path + query + m_restapi_directives;
    if (200 == aurostd::httpGetStatus(url, output)){
      return output;
    } else {
      m_logger_message << "Failed to get REST API query " << "(" << url << ")";
      outError(__AFLOW_FUNC__, __LINE__);
      return "";
    }
  }

  /// @brief add an xstructure to an AFLOW lib entry
  /// @param entry AFLOW lib entry
  /// @param orig if `true` load the original structure otherwise the final structure (default)
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note use this function only if you want the xstructure for a subset of entries;
  ///       if you want xstructures for all entries set #m_xstructure_relaxed or #m_xstructure_original to `true`
  /// @note adds the structure to entry.vstr
  void EntryLoader::addXstructure(aflowlib::_aflowlib_entry &entry, bool orig) {
    if (orig) {
      xstructure new_structure;
      if (loadXstructureAflowIn(entry, new_structure, 1)) { // load from aflow.in
        entry.vstr.push_back(new_structure);
      } else {
        m_logger_message << "Failed to add original structure to " << entry.auid << " (" << entry.aurl << ")";
        outError(__AFLOW_FUNC__, __LINE__);
      }
    } else {
      xstructure new_structure;
      if (pflow::loadXstructureLibEntry(entry, new_structure)) { // load directly from the entry (AFLUX, SQLITE)
        entry.vstr.push_back(new_structure);
      } else if (loadXstructureFile(entry, new_structure)){
        entry.vstr.push_back(new_structure);
      } else if (loadXstructureAflowIn(entry, new_structure, -1)) { // load from aflow.in
        entry.vstr.push_back(new_structure);
      } else {
        m_logger_message << "Failed to add relaxed structure to " << entry.auid << " (" << entry.aurl << ")";
        outError(__AFLOW_FUNC__, __LINE__);
      }
    }
  }

  /// @brief load a structure from a structure file
  /// @param entry AFLOW lib entry
  /// @param new_structure save-to structure
  /// @param possible_files load the structure from the first exiting sources given (if list is empty use m_xstructure_final_file_name)
  /// @return xstructure
  /// @note does not add the structure to entry.vstr
  bool EntryLoader::loadXstructureFile(const aflowlib::_aflowlib_entry &entry, xstructure &new_structure, std::vector <std::string> possible_files) {
    std::string base_url = m_restapi_server + m_restapi_path + entry.aurl.substr(28) + "/";
    std::string base_folder;
    bool direct_fs = entry.el_source_type.substr(0, 10) == "FILESYSTEM"; // catch FILESYSTEM and FILESYSTEM_RAW
    if (direct_fs){
      base_folder = entry.el_source.substr(0, entry.el_source.size()-m_filesystem_outfile.size());
    } else {
      base_folder = m_filesystem_path + entry.aurl.substr(28) + "/";
      base_folder = std::regex_replace(base_folder, aurostd::regex_aurl2file, "$1/" + m_filesystem_collection + "/");
    }
    std::string poscar;
    if (entry.catalog =="LIB0" && !m_filesystem_available && entry.aurl.substr(entry.aurl.size() - 2)=="/0"){
      return false; // no entries in RESTAPI for WEB0 /0
    }

    // if no file names given use the class defaults
    if (possible_files.empty()) possible_files = m_xstructure_final_file_name;

    std::vector <std::string> available_files;
    if (m_filesystem_available || direct_fs) {
      aurostd::DirectoryLS(base_folder, available_files);
    } else {
      listRestAPI(base_url, available_files, false);
    }

    std::string selected_file;
    for (std::vector<std::string>::const_iterator file_name = possible_files.begin(); file_name != possible_files.end(); file_name++) {
      if (aurostd::EWithinList(available_files, *file_name, selected_file)) {
        if (m_filesystem_available || direct_fs) aurostd::efile2string(base_folder + selected_file, poscar);
        else poscar = getRawRestAPIQuery(base_url + selected_file, true);
        if (!poscar.empty()) { // load from aflow.in
          new_structure = xstructure((std::stringstream) poscar, IOVASP_AUTO);
          return true;
        }
      }
    }
    return false;
  }

  /// @brief load the first structure in an `aflow.in` file
  /// @param entry AFLOW lib entry
  /// @param new_structure save-to structure
  /// @authors
  /// @mod{HE,20220216,created}
  /// @return xstructure
  /// @note this is always the original structure
  /// @note does not add the structure to entry.vstr
  bool EntryLoader::loadXstructureAflowIn(const aflowlib::_aflowlib_entry &entry, xstructure &new_structure, const int index) {
    std::string base_folder;
    bool direct_fs = entry.el_source_type.substr(0, 10) == "FILESYSTEM"; // catch FILESYSTEM and FILESYSTEM_RAW
    if (direct_fs){
      base_folder = entry.el_source.substr(0, entry.el_source.size()-m_filesystem_outfile.size());
    } else {
      base_folder = m_filesystem_path + entry.aurl.substr(28) + "/";
      base_folder = std::regex_replace(base_folder, aurostd::regex_aurl2file, "$1/" + m_filesystem_collection + "/");
    }
    std::string base_url = m_restapi_server + m_restapi_path + entry.aurl.substr(28) + "/";
    std::string aflowin_content;

    if (m_filesystem_available || direct_fs) {
      std::stringstream buffer;
      std::ifstream open_file(base_folder + "aflow.in");
      buffer << open_file.rdbuf();
      aflowin_content = buffer.str();
    } else { // load form REST API
      m_out_super_silent = true;
      aflowin_content = getRawRestAPIQuery(base_url + "aflow.in", true);
      m_out_super_silent = false;
    }
    std::stringstream poscar;
    bool poscar_loaded =  KBIN::ExtractPOSCARStringStreamFromAFLOWIN(aflowin_content, poscar, index);
    if (poscar_loaded) { // load from aflow.in
      new_structure = xstructure((std::stringstream) poscar.str(), IOVASP_AUTO);
      return true;
    }
    return false;
  }

  /// @brief save the shared pointers to the AFLOW lib entries into an external vector
  /// @param result save-to vector
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note creating a copy the smart pointer #m_entries_flat is a bit more efficient way to use
  ///       the loaded entries after the EntryLoader class goes out of scope
  void EntryLoader::getEntriesViewFlat(std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> &result) const {
    result = *m_entries_flat;
  }

  /// @brief save the shared pointers to the AFLOW lib entries into a two layer vector
  /// @param result save-to vector
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note the entries are grouped by number of elements (smallest to largest)
  void EntryLoader::getEntriesViewTwoLayer(vector<vector<std::shared_ptr < aflowlib::_aflowlib_entry>> > &result) const{
  for (std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>::iterator layer1 = m_entries_layered_map->begin(); layer1 != m_entries_layered_map->end(); layer1++) {
  std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> collected_entries;
    for (std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>::iterator layer2 = layer1->second.begin(); layer2 != layer1->second.end(); layer2++) {
      for (std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>::iterator entry = layer2->second.begin(); entry != layer2->second.end(); entry++) collected_entries.push_back(*entry);
      }
      result.push_back(collected_entries);
    }
  }

  /// @brief save the shared pointers to the AFLOW lib entries into a three layer vector
  /// @param result save-to vector
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note the entries are grouped first by number of elements (smallest to largest) and then by alloy
  void EntryLoader::getEntriesViewThreeLayer(std::vector<std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> &result) const{
    for (std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>::iterator layer1 = m_entries_layered_map->begin(); layer1 != m_entries_layered_map->end(); layer1++) {
      std::vector < std::vector < std::shared_ptr < aflowlib::_aflowlib_entry>>> collected_entries_l1;
      for (std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>::iterator layer2 = layer1->second.begin(); layer2 != layer1->second.end(); layer2++) {
        std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> collected_entries_l2;
        for (std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>::iterator entry = layer2->second.begin(); entry != layer2->second.end(); entry++) collected_entries_l2.push_back(*entry);
        collected_entries_l1.push_back(collected_entries_l2);
      }
      result.push_back(collected_entries_l1);
    }
  }

  /// @brief save the shared pointers to the AFLOW lib entries into a three layer map
  /// @param result save-to map
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note the entries are grouped first by number of elements and then by alloy
  /// @note creating a copy of the smart pointer #m_entries_layered_map is a bit more efficient way to use
  ///       the loaded entries after the EntryLoader class goes out of scope
  void EntryLoader::getEntriesViewMap(std::map < short, std::map < std::string,
                                      std::vector < std::shared_ptr < aflowlib::_aflowlib_entry >> >> &result) const {
    result = *m_entries_layered_map;
  }

  /// @brief copy the AFLOW lib entries into a vector
  /// @param result save-to vector
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries are copied into a continuous chunk of memory
  void EntryLoader::getEntriesFlat(std::vector <aflowlib::_aflowlib_entry> &result) const {
    for (std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>::iterator entry = m_entries_flat->begin(); entry != m_entries_flat->end(); entry++){
      result.emplace_back(*(*entry));
    }
  }

  /// @brief copy the AFLOW lib entries into a two layer vector
  /// @param result save-to vector
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries are copied into a continuous chunk of memory
  /// @note the entries are grouped by number of elements (smallest to largest)
  void EntryLoader::getEntriesTwoLayer(std::vector <std::vector<aflowlib::_aflowlib_entry>> &result) const {
    for (std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>::iterator layer1 = m_entries_layered_map->begin(); layer1 != m_entries_layered_map->end(); layer1++) {
      std::vector <aflowlib::_aflowlib_entry> collected_entries;
      for (std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>::iterator layer2 = layer1->second.begin(); layer2 != layer1->second.end(); layer2++) {
        for (std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>::iterator entry = layer2->second.begin(); entry != layer2->second.end(); entry++) collected_entries.push_back(*(*entry));
      }
      result.push_back(collected_entries);
    }
  }


  /// @brief copy the AFLOW lib entries into a three layer vector
  /// @param result save-to vector
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the underlying entries are copied into a continuous chunk of memory
  /// @note the entries are grouped first by number of elements (smallest to largest) and then by alloy
  void EntryLoader::getEntriesThreeLayer(std::vector<std::vector<vector<aflowlib::_aflowlib_entry>>> & result) const {
    for (std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>::iterator layer1 = m_entries_layered_map->begin(); layer1 != m_entries_layered_map->end(); layer1++) {
      std::vector <std::vector<aflowlib::_aflowlib_entry>> collected_entries_l1;
      for (std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>::iterator layer2 = layer1->second.begin(); layer2 != layer1->second.end(); layer2++) {
        std::vector <aflowlib::_aflowlib_entry> collected_entries_l2;
        for (std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>::iterator entry = layer2->second.begin(); entry != layer2->second.end(); entry++) collected_entries_l2.push_back(*(*entry));
        collected_entries_l1.push_back(collected_entries_l2);
      }
      result.push_back(collected_entries_l1);
    }
  }

  // private functions
  /// @brief find the best available source
  /// @authors
  /// @mod{HE,20220216,created}
  void EntryLoader::selectSource() {
    if (m_current_source == Source::NONE){
      EntryLoader::getCommandlineSource();
    }

    if (m_current_source == Source::NONE || m_current_source == Source::FAILED) {
      m_out_super_silent = true;
      if (EntryLoader::setSource(Source::SQLITE)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected Source::SQLITE";
        outInfo(__AFLOW_FUNC__);
        return;
      }
      if (EntryLoader::setSource(Source::AFLUX)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected Source::AFLUX";
        outInfo(__AFLOW_FUNC__);
        return;
      }
      if (EntryLoader::setSource(Source::FILESYSTEM)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected";
        if (m_current_source==Source::FILESYSTEM) m_logger_message << " Source::FILESYSTEM";
        else if (m_current_source==Source::FILESYSTEM_RAW) m_logger_message << " Source::FILESYSTEM_RAW";
        outInfo(__AFLOW_FUNC__);
        return;
      }
      if (EntryLoader::setSource(Source::RESTAPI)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected";
        if (m_current_source==Source::RESTAPI) m_logger_message << " Source::RESTAPI";
        else if (m_current_source==Source::RESTAPI_RAW) m_logger_message << " Source::RESTAPI_RAW";
        outInfo(__AFLOW_FUNC__);
        return;
      }
      EntryLoader::setSource(Source::FAILED);
      m_out_super_silent = false;
      m_logger_message << "Failed to find a working source!";
      outHardError(__AFLOW_FUNC__, __LINE__, _GENERIC_ERROR_);
    }
  }

  bool EntryLoader::getCommandlineSource(){
    string raw_source = XHOST.vflag_control.getattachedscheme("ENTRY_LOADER::SOURCE");
    raw_source = aurostd::tolower(raw_source);
    if (raw_source == "aflux") return EntryLoader::setSource(Source::AFLUX);
    else if (raw_source == "sqlite") return EntryLoader::setSource(Source::SQLITE);
    else if (raw_source == "file" || raw_source == "filesystem") return EntryLoader::setSource(Source::FILESYSTEM);
    else if (raw_source == "rest" || raw_source == "restapi") return EntryLoader::setSource(Source::RESTAPI);
    else return false;
  }

  /// @brief list directories or files for the REST API
  /// @authors
  /// @mod{HE,20220216,created}
  /// @param url entry url
  /// @param result save-to vector
  /// @param directories list directories if `true` (default)
  void EntryLoader::listRestAPI(std::string url, std::vector<std::string> &result, bool directories) {
    result.clear();
    if (directories) url += m_restapi_listing_dirs;
    else url += m_restapi_listing_files;
    std::string output;
    if (200 == aurostd::httpGetStatus(url, output)) {
      if (output[output.size()-1]=='\n') output.erase(output.size()-1);
      aurostd::string2tokens(output, result, ",");
    } else {
      m_logger_message << "Could not list content for " << url;
      outInfo(__AFLOW_FUNC__);
    };
  }

  /// @brief generate a list of AUID for a given list of alloys by querying the public alloy SQLITE DB
  /// @param alloy_list list of cleaned and sorted alloys
  /// @param auid_list resulting AUIDs
  /// @authors
  /// @mod{HE,20220216,created}
  void EntryLoader::getAlloyAUIDList(const std::vector<std::string> & alloy_list, std::vector<std::string> &auid_list) {
    if (m_sqlite_alloy_db_ptr == nullptr) {
      m_logger_message << "Public alloy SQLITE DB is not ready";
      outHardError(__AFLOW_FUNC__, __LINE__, _RUNTIME_ERROR_);
      return;
    }
    auid_list.clear();
    std::string where = aurostd::joinWDelimiter(alloy_list, "','");
    where = "alloy IN ('" + where + "')";
    auid_list = m_sqlite_alloy_db_ptr->getValuesMultiTable("auid", where);
  }

  /// @brief load AFLOW lib entries for a given alloy directly from the filesystem (source FILESYSTEM_RAW)
  /// @param alloy_list list of cleaned and sorted alloys
  /// @param lib_max number of elements in the largest alloy
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note if the public alloy SQLITE DB is available the results are expanded with missed entries
  void EntryLoader::loadAlloySearchFSR(const std::vector <std::string> &alloy_list, uint lib_max) {
    std::vector <std::string> found_entries;
    std::vector <std::string> search_path_list;
    search_path_list.emplace_back(m_filesystem_path + "ICSD/" + m_filesystem_collection + "/");
    if (alloy_list.size()>1) {
      for (uint lib_idx = 1; lib_idx <= lib_max; lib_idx++) {
        search_path_list.emplace_back(
            m_filesystem_path + "LIB" + std::to_string(lib_idx) + "/" + m_filesystem_collection + "/");
      }
    } else if (alloy_list.size()==1){
      search_path_list.emplace_back(
          m_filesystem_path + "LIB" + std::to_string(lib_max) + "/" + m_filesystem_collection + "/");
    } else {
      return;
    }

    char **paths = new char*[search_path_list.size()+1](); // initialize the array with nullptr, +1 ensures it is null terminated
    for (size_t list_idx = 0; list_idx < search_path_list.size(); list_idx++) {
      paths[list_idx] = (char *) search_path_list[list_idx].c_str();
    }

    size_t scanned = 0;
    size_t found = 0;
    size_t check_idx = m_filesystem_path.size();
    struct stat file_stats{};

    // initialize a file tree scan
    // https://man7.org/linux/man-pages/man3/fts.3.html
    // FTS_PHYSICAL - don't follow symlinks
    // FTS_NOCHDIR - don't change the workdir of the program
    // FTS_XDEV - don't descend into folders that are on a different device
    FTS *tree = fts_open(paths, FTS_PHYSICAL | FTS_NOCHDIR | FTS_XDEV, nullptr);
    FTSENT *node;
    if (tree == nullptr) {
      m_logger_message << "Failed to initialize the file tree used to search for alloys!";
      outHardError(__AFLOW_FUNC__, __LINE__,_FILE_ERROR_);
      return;
    }

    // Iterate over all entries found in the folders listed in paths
    while ((node = fts_read(tree))) {
      scanned += 1;
      if ((scanned-1) % 10000 == 0) {
        m_logger_message << (scanned-1) << " objects scanned; " << found << " entries found; next scan: " << node->fts_path;
        outDebug(__AFLOW_FUNC__);
      }
      if ((node->fts_info & FTS_D)) {
        if (node->fts_path[check_idx] == 'I') { // ICSD
          if (node->fts_level == 2) {
            if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(node->fts_name, 'I')) == alloy_list.end()) {
              fts_set(tree, node, FTS_SKIP);
              continue;
            } else {
              std::string base_path = node->fts_path;
              std::string full_path = base_path + "/" + m_filesystem_outfile;
              if (stat(full_path.c_str(), &file_stats) == 0) {
                found += 1;
                found_entries.emplace_back(full_path);
                if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
              }
            }
          }
        } else if (node->fts_path[check_idx] == 'L') { // LIBX
          if (node->fts_level == 1) {
            if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(node->fts_name, 'L')) == alloy_list.end()) {
              fts_set(tree, node, FTS_SKIP);
              continue;
            } else {
              std::string base_path = node->fts_path;
              std::string full_path = base_path + "/" + m_filesystem_outfile;
              if (stat(full_path.c_str(), &file_stats) == 0) {
                found += 1;
                found_entries.emplace_back(full_path);
                if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
              }
            }
          } else if (node->fts_level > 1) {
            std::string base_path = node->fts_path;
            std::string full_path = base_path + "/" + m_filesystem_outfile;
            if (stat(full_path.c_str(), &file_stats) == 0) {
              found += 1;
              found_entries.emplace_back(full_path);
              if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
            }
          }
        }
      }
    }
    fts_close(tree);
    delete[] paths;


    m_logger_message << "Finishing search in the filesystem after scanning " << scanned << " objects";
    outInfo(__AFLOW_FUNC__);
    loadFiles(found_entries);
    if (m_sqlite_alloy_db_ptr != nullptr) {
      std::vector <std::string> known_AUID_list;
      getAlloyAUIDList(alloy_list, known_AUID_list);
      std::vector <std::string> missing_AUID;
      for (std::vector<std::string>::const_iterator AUID = known_AUID_list.begin(); AUID != known_AUID_list.end(); AUID++) {
        if (std::find(m_auid_list.begin(), m_auid_list.end(), *AUID) == m_auid_list.end()) missing_AUID.emplace_back(*AUID);
      }
      m_logger_message << "Found " << missing_AUID.size() << " entries in the public alloy SQLITE DB that where not found in the filesystem search.";
      outDebug(__AFLOW_FUNC__);
      loadAUID(missing_AUID);
    }
  }

  /// @brief load AFLOW lib entries for a given alloy set directly from the AFLOW REST API (source RESTAPI_RAW)
  /// @param alloy_list list of cleaned and sorted alloys
  /// @param lib_max number of elements in the largest alloy
  /// @authors
  /// @mod{HE,20220216,created}
  void EntryLoader::loadAlloySearchRR(const std::vector <std::string> & alloy_list, uint lib_max) {
    std::vector<std::string> rest_api_queries;
    std::vector<std::string> icsd_search_path_list;
    std::vector<std::string> libx_search_path_list;
    std::vector<std::string> libx_crawl_list;
    for (std::vector<std::string>::const_iterator sym = BRAVAIS_LATTICES.begin(); sym != BRAVAIS_LATTICES.end(); sym++) {
      icsd_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "ICSD" + "_" + m_restapi_collection + "/" + *sym + "/");
    }
    if (alloy_list.size()>1) {
      for (uint lib_idx = 1; lib_idx <= lib_max; lib_idx++) {
        libx_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "LIB" + std::to_string(lib_idx) + "_" + m_restapi_collection + "/" );
      }
    } else if (alloy_list.size()==1){
      libx_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "LIB" + std::to_string(lib_max) + "_" + m_restapi_collection + "/" );
    } else {
      return;
    }
    std::vector<std::string> listing;
    for (std::vector<std::string>::const_iterator url = icsd_search_path_list.begin(); url != icsd_search_path_list.end(); url++) {
      listRestAPI(*url, listing);
      for (std::vector<std::string>::const_iterator name = listing.begin(); name != listing.end(); name++) {
        if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(*name, 'I')) != alloy_list.end()){
          rest_api_queries.emplace_back(*url + *name + "/" + m_restapi_directives);
        }
      }
    }

    m_logger_message << "Found " << rest_api_queries.size() << " ICSD entries in REST API search";
    outDebug(__AFLOW_FUNC__);
    for (std::vector<std::string>::const_iterator url = libx_search_path_list.begin(); url != libx_search_path_list.end(); url++) {
      listRestAPI(*url, listing);
      for (std::vector<std::string>::const_iterator name = listing.begin(); name != listing.end(); name++) {
        if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(*name, 'L')) != alloy_list.end()){
          libx_crawl_list.emplace_back(*url + *name + "/" );
        }
      }
    }

    m_logger_message << "Found " << libx_crawl_list.size() << " LIBX base folders in the REST API search";
    outDebug(__AFLOW_FUNC__);

    uint scanned=0;
    uint found=0;
    while(!libx_crawl_list.empty()){
      scanned += 1;
      std::string url = libx_crawl_list.back();
      libx_crawl_list.pop_back();
      if ((scanned-1)%100 == 0) {
        m_logger_message << (scanned-1) << " folders scanned; " << found << " entries found; next scan: " << url;
        outDebug(__AFLOW_FUNC__);
      }
      listRestAPI(url, listing);
      if (listing.empty()) {
        rest_api_queries.emplace_back(url + m_restapi_directives);
        found += 1;
      } else {
        for (std::vector<std::string>::const_iterator name = listing.begin(); name != listing.end(); name++) {
          libx_crawl_list.emplace_back(url + *name + "/");
        }
      }
    }


    loadRestAPIQueries(rest_api_queries, true);
    if (m_sqlite_alloy_db_ptr != nullptr) {
      std::vector <std::string> known_AUID_list;
      getAlloyAUIDList(alloy_list, known_AUID_list);
      std::vector <std::string> missing_AUID;
      for (std::vector<std::string>::const_iterator AUID = known_AUID_list.begin(); AUID != known_AUID_list.end(); AUID++) {
        if (std::find(m_auid_list.begin(), m_auid_list.end(), *AUID)== m_auid_list.end()) missing_AUID.emplace_back(*AUID);
      }
      m_logger_message << "Found " << missing_AUID.size() << " entries in the public alloy SQLITE DB that where not found in the REST API search.";
      outDebug(__AFLOW_FUNC__);
      loadAUID(missing_AUID);
    }

  }

  /// @brief map different AUID styles to a default
  /// @param AUID AUID to clean
  /// @return `false` if AUID is not valid
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the AUID is cleaned in place
  bool EntryLoader::cleanAUID(std::string & AUID) {
    if (AUID.find('\"')!=std::string::npos) aurostd::StringSubst(AUID, "\"", "");
    if (AUID.substr(0, 5) == "auid:") AUID="aflow:" + AUID.substr(5);
    else if (AUID.substr(0, 6) != "aflow:") AUID="aflow:" + AUID;
    if (aurostd::_ishex(AUID.substr(6)) && AUID.size()==22) return true;
    else return false;
  }

  /// @brief map different AURL styles to a default
  /// @param AUID AURL to clean
  /// @return `false` if AURL is not valid
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note AURL is not validated
  /// @note the AURL is cleaned in place
  bool EntryLoader::cleanAURL(std::string & AURL) {
    if (AURL.substr(0, 5) == "aurl:")  AURL = AURL.substr(5);
    if (AURL.substr(0, 28)== "aflowlib.duke.edu:AFLOWDATA/") return true;
    if (AURL.substr(0, 9) == "AFLOWDATA") {
      AURL = "aflowlib.duke.edu:" + AURL;
      return true;
    }
    if ((AURL.substr(0,3) == "LIB" || AURL.substr(0,4) == "ICSD")) {
      AURL = "aflowlib.duke.edu:AFLOWDATA/" + AURL;
      return true;
    }
    return false;
  }

  /// @brief create a AFLUX query string from a matchbook
  /// @param matchbook map with keyword and modifier
  /// @returns AFLUX query
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note does not include #m_aflux_server or #m_aflux_path
  /// @note the string is percent encoded
  std::string EntryLoader::buildAFLUXQuery(const std::map<std::string, std::string> & matchbook) const {
    std::string query = "";
    std::vector<std::map<std::string, std::string>> joined_directives = {matchbook, m_aflux_directives};
    for (std::vector<std::map<std::string, std::string>>::const_iterator part = joined_directives.begin(); part != joined_directives.end(); part++) {
      for (std::map<std::string, std::string>::const_iterator it = part->begin(); it != part->end(); it++) {
        query += "," + it->first;
        if (!it->second.empty()) query += "(" + it->second + ")";
      }
    }
    query.erase(0,1);
    return "?" + aurostd::httpPercentEncodingFull(query);
  }

  /// @brief extract alloy string from AFLOW run name
  /// @param name AFLOW run name (from path structure)
  /// @param lib_type indicate the LIB type | `I` for `ICSD` or `L` for 'LIBx'
  /// @return alloy string
  /// @authors
  /// @mod{HE,20220216,created}
  /// @note the alloy sting is sorted
  std::string EntryLoader::extractAlloy(std::string name, char lib_type) const {
    if (lib_type == 'I') { //ICSD
      name.erase(std::min(name.find("_ICSD_"), name.size()));
    } else if (lib_type == 'L') { //LIBX
      name.erase(std::min(name.find(':'), name.size()));
      name = std::regex_replace(name, aurostd::regex_ppclean,"");
    } else {
      return "";
    }
    std::vector <std::string> element_match(std::sregex_token_iterator(name.begin(), name.end(), aurostd::regex_elements),
                                            std::sregex_token_iterator());
    std::sort(element_match.begin(), element_match.end());
    return aurostd::joinWDelimiter(element_match, "");
  }

  void EntryLoader::outInfo(const std::string & function_name) {
    if ((!m_out_silent || m_out_debug) && !m_out_super_silent){
      pflow::logger(__FILE__, function_name, m_logger_message, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
    }
    aurostd::StringstreamClean(m_logger_message);
  }

  void EntryLoader::outDebug(const std::string & function_name) {
    if (m_out_debug && !m_out_super_silent){
      pflow::logger(__FILE__, function_name, m_logger_message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }
    aurostd::StringstreamClean(m_logger_message);
  }

  void EntryLoader::outError(const std::string &function_name, int line_number) {
    if (!m_out_super_silent) {
      pflow::logger((std::string) __FILE__ + ":" + std::to_string(line_number),
                    function_name, m_logger_message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    aurostd::StringstreamClean(m_logger_message);
  }

  void EntryLoader::outHardError(const std::string &function_name, int line_number, int error_type) {
    throw aurostd::xerror((std::string) __FILE__ + ":" + std::to_string(line_number), function_name, m_logger_message, error_type);
    aurostd::StringstreamClean(m_logger_message);
  }


}

#endif  // _AFLOWLIB_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
