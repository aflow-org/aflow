// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by hagen.eckert@duke.edu

#ifndef _AUROSTD_DATA_CPP_
#define _AUROSTD_DATA_CPP_

#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <map>
#include <string>
#include <utility>

#include <archive.h>
#include <archive_entry.h>

#include "aurostd.h"
#include "aurostd_incbin_defs.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"
#include "aurostd_xparser_json.h"

#ifndef CMAKE_SOURCE_DIR
#define CMAKE_SOURCE_DIR ""
#endif

namespace aurostd::EmbData {
// Include files either as text (INCTXT) or binary (INCBIN) into the aflow binary
// see AUROSTD/aurostd_incbin_defs.h
  extern "C" {
    // strings used for unittests
  INCBIN(FINDSYM, CMAKE_SOURCE_DIR "DATA/FINDSYM.tar.xz");
  INCBIN(FROZSL, CMAKE_SOURCE_DIR "DATA/FROZSL.tar.xz");
  INCBIN(PHVASP, CMAKE_SOURCE_DIR "DATA/PHVASP.tar.xz");
  INCBIN(ReadMe, CMAKE_SOURCE_DIR "DATA/README.tar.xz");
  INCBIN(Test, CMAKE_SOURCE_DIR "DATA/TEST.tar.xz");
  INCBIN(Images, CMAKE_SOURCE_DIR "DATA/IMAGES.tar.xz");
  INCBIN(Scripts, CMAKE_SOURCE_DIR "DATA/SCRIPTS.tar.xz");
  INCBIN(Prototypes, CMAKE_SOURCE_DIR "DATA/PROTO.tar.xz");
  }

  static const std::map<std::string, std::pair<char *, unsigned int>> aflow_data_collections = {
      { "README",         {(char *) &afdataReadMeData[0], afdataReadMeSize}},
      {"FINDSYM",       {(char *) &afdataFINDSYMData[0], afdataFINDSYMSize}},
      { "FROZSL",         {(char *) &afdataFROZSLData[0], afdataFROZSLSize}},
      { "PHVASP",         {(char *) &afdataPHVASPData[0], afdataPHVASPSize}},
      {   "TEST",             {(char *) &afdataTestData[0], afdataTestSize}},
      { "IMAGES",         {(char *) &afdataImagesData[0], afdataImagesSize}},
      {"SCRIPTS",       {(char *) &afdataScriptsData[0], afdataScriptsSize}},
      {  "PROTO", {(char *) &afdataPrototypesData[0], afdataPrototypesSize}},
  };

  std::string read_data(struct archive *ar) {
    std::string content;
    const size_t buff_size = 16384;
    char buff[buff_size + 1];
    size_t size;
    for (;;) {
      size = archive_read_data(ar, &buff, buff_size);
      if (size == 0) {
        break;
      }
      buff[size] = '\0';
      content += buff;
    }
    return content;
  }

  void write_collection_folder(const std::pair<char *, unsigned int> raw_data, const std::string &path, const std::string &folder_out_raw) {
    if (path.back() != '/') {
      const std::string message = "The requested path (" + path + ") is not a folder - use write_collection_file";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    struct archive *a = archive_read_new();
    struct archive_entry *entry;
    static char buff[16384];
    bool found = false;
    size_t len;
    int r;
    std::filesystem::path file_out;
    std::filesystem::path raw_out;
    const std::filesystem::path folder_out = std::filesystem::path(aurostd::CleanFileName(folder_out_raw));
    archive_read_support_format_tar(a);
    archive_read_support_filter_xz(a);

    if (archive_read_open_memory(a, raw_data.first, raw_data.second)) {
      const std::string message = "Embedded archive could not be open while trying to read \"" + path + "\". " + archive_error_string(a);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }

    for (;;) {
      r = archive_read_next_header(a, &entry);
      if (r == ARCHIVE_EOF) {
        break;
      }
      if (r != ARCHIVE_OK) {
        const std::string message = "Failed to open entries in archive while trying to read \"" + path + "\". " + archive_error_string(a);
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }

      // check if the beginning of the entries filename equals the path
      if (strncmp(archive_entry_pathname(entry), path.c_str(), path.size()) == 0) {
        // build the new file path relative to folder_out
        raw_out = std::filesystem::path(archive_entry_pathname(entry));
        file_out = folder_out / std::filesystem::proximate(raw_out, path);
        // ensure that the folder exist
        std::filesystem::create_directories(file_out.parent_path());
        // write the data
        std::fstream fs(file_out, std::ios::out);
        len = archive_read_data(a, buff, sizeof(buff));
        while (len > 0) { // read until archive is empty
          fs.write(buff, len);
          len = archive_read_data(a, buff, sizeof(buff));
        }
        fs.close();
        found = true;
      } else {
        archive_read_data_skip(a);
      }
    }
    if (not found) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Failed to find any entry at \"" + path + "\" in archive.", _FILE_NOT_FOUND_);
    }
    archive_read_free(a);
  }

  void write_collection_file(const std::pair<char *, unsigned int> raw_data, const std::string &filename, const std::string &outfile_raw) {
    string outfile(CleanFileName(outfile_raw));
    struct archive *a = archive_read_new();
    struct archive_entry *entry;
    static char buff[16384];
    bool found = false;
    size_t len;

    const std::filesystem::path path = std::filesystem::path(outfile);
    archive_read_support_format_tar(a);
    archive_read_support_filter_xz(a);

    int r;
    if (archive_read_open_memory(a, raw_data.first, raw_data.second)) {
      const std::string message = "Embedded archive could not be open while trying to read \"" + filename + "\". " + archive_error_string(a);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    for (;;) {
      r = archive_read_next_header(a, &entry);
      if (r == ARCHIVE_EOF) {
        break;
      }
      if (r != ARCHIVE_OK) {
        const std::string message = "Failed to open entries in archive while trying to read \"" + filename + "\". " + archive_error_string(a);
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }

      if (archive_entry_pathname(entry) == filename) {
        // open the new file and write it in chunks
        if (outfile.empty()) {
          outfile = path.parent_path().string() + "/" + path.stem().string();
        }
        std::fstream fs(outfile, std::ios::out);
        len = archive_read_data(a, buff, sizeof(buff));
        while (len > 0) { // read until archive is empty
          fs.write(buff, len);
          len = archive_read_data(a, buff, sizeof(buff));
        }
        fs.close();
        found = true;
        break;
      } else {
        archive_read_data_skip(a);
      }
    }
    if (not found) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Failed to find \"" + filename + "\" in archive.", _FILE_NOT_FOUND_);
    }
    archive_read_free(a);
  }

  std::string get_collection_text(const std::pair<char *, unsigned int> raw_data, const std::string &filename) {
    struct archive *a = archive_read_new();
    struct archive_entry *entry;
    bool found = false;
    archive_read_support_format_tar(a);
    archive_read_support_filter_xz(a);
    std::string content;
    int r;

    if (archive_read_open_memory(a, raw_data.first, raw_data.second)) {
      const std::string message = "Embedded archive could not be open while trying to read \"" + filename + "\". " + archive_error_string(a);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
    for (;;) {
      r = archive_read_next_header(a, &entry);
      if (r == ARCHIVE_EOF) {
        break;
      }
      if (r != ARCHIVE_OK) {
        const std::string message = "Failed to open entries in archive while trying to read \"" + filename + "\". " + archive_error_string(a);
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
      }
      if (archive_entry_pathname(entry) == filename) {
        content = read_data(a);
        found = true;
        break;
      } else {
        archive_read_data_skip(a);
      }
    }
    if (not found) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Failed to find \"" + filename + "\" in archive.", _FILE_NOT_FOUND_);
    }
    archive_read_free(a);
    return content;
  }

  /// @brief return data for unittests
  /// @return content of `DATA/TEST/unit_test.json`
  aurostd::JSON::object get_unit_test() {
    return aurostd::JSON::loadString(get_content("unit_test.json", "TEST"));
  }

  /// @brief return data for unittests
  /// @return content of DATA/TEST/<file>`
  string get_test_file(string file) {
    return get_content(file, "TEST");
  }

  /// @brief get the content of an embedded file
  /// @param filename name of the file
  /// @param collection name of the collection
  /// @return content of the embedded file
  std::string get_content(const std::string &filename, const std::string &collection) {
    if (aflow_data_collections.count(collection)) {
      return get_collection_text(aflow_data_collections.at(collection), filename);
    } else {
      const string message = "Collection \"" + collection + " \" was not embedded";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
  }

  /// @brief write a embedded file into the filesystem
  /// @param filename name of the file (as in the DATA folder)
  /// @param target_path target location
  /// @authors
  /// @mod{HE,20230413,create function}
  void save_to_file(const std::string &filename, const std::string &collection, const std::string &target_path) {
    if (aflow_data_collections.count(collection)) {
      write_collection_file(aflow_data_collections.at(collection), filename, target_path);
    } else {
      const string message = "Collection \"" + collection + " \" was not embedded";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
  }

  /// @brief write a embedded file into the filesystem
  /// @param path path a folder in the archive (must end in '/')
  /// @param target_path target location
  /// @authors
  /// @mod{HE,20240403,create function}
  void save_to_folder(const std::string &path, const std::string &collection, const std::string &target_path) {
    if (aflow_data_collections.count(collection)) {
      write_collection_folder(aflow_data_collections.at(collection), path, target_path);
    } else {
      const string message = "Collection \"" + collection + " \" was not embedded";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
  }

} // end namespace aurostd::EmbData

#endif  //_AUROSTD_DATA_CPP_
