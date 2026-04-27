// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// This file hold function dealing with files
// hagen.eckert@duke.edu

#include "aurostd_xfile.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <archive.h>
#include <archive_entry.h>

#include "aurostd.h"
#include "aurostd_defs.h"
#include "aurostd_time.h"
#include "aurostd_xerror.h"
#include "aurostd_xrandom.h"

#include "aflow_xhost.h" // todo required for XPID use and XHOST.DEBUG use

using std::cerr;
using std::deque;
using std::endl;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

namespace aurostd {
  /// @brief create a directory
  /// @param directory_raw to be created
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool DirectoryMake(const string& directory_raw) {
    const string directory = CleanFileName(directory_raw);
    if (!std::filesystem::is_directory(directory)) {
      return std::filesystem::create_directories(directory);
    }
    return true;
  }

  /// @brief create a directory on a remote computer
  /// @note need a working public key setup
  /// @authors
  /// @mod{SC,,created}
  [[deprecated("Should be rewritten without using system calls")]]
  bool SSH_DirectoryMake(const string& user, const string& machine, const string& directory_raw) {  // "" compliant SC20190401
    const bool LDEBUG = (false || XHOST.DEBUG);
    const string Directory(aurostd::CleanFileName(directory_raw));
    std::deque<string> tokens;
    string dir_tmp;
    string command;
    aurostd::string2tokens(Directory, tokens, "/");
    if (Directory.at(0) == '/') {
      tokens[0] = "/" + tokens[0]; // fix the root thing
    }
    for (const auto& token : tokens) {
      dir_tmp += token + "/";
      command = "ssh " + user + "@" + machine + " mkdir \"" + dir_tmp + "\"";
      if (LDEBUG) {
        cerr << "SSH_DirectoryMake creating directory=" << command << endl;
      }
      aurostd::execute(command);
    }
    return true;
  }

  /// @brief returns all the subdirectories of a directory, including itself
  /// @param directory_raw to inspect
  /// @param vsubd list of directories created by this function
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool SubDirectoryLS(const string& directory_raw, vector<string>& vsubd) {
    const string directory = CleanFileName(directory_raw);
    if (directory.empty()) {
      return false;
    }
    vsubd.clear();
    vsubd.push_back(directory);
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(directory)) {
      if (dir_entry.is_directory()) {
        vsubd.push_back(dir_entry.path().filename().string());
      }
    }
    return true;
  }

  /// @brief returns the content of a directory
  /// @param directory_raw to inspect
  /// @param vfiles list of content
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool DirectoryLS(const string& directory_raw, vector<string>& vfiles) {
    const string directory = CleanFileName(directory_raw);
    if (directory.empty()) {
      return false;
    }
    vfiles.clear();
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(directory)) {
      vfiles.push_back(dir_entry.path().filename().string());
    }
    return true;
  }

  bool DirectoryLS(const string& directory_raw, deque<string>& vfiles) {
    vector<string> _vfiles;
    if (DirectoryLS(directory_raw, _vfiles)) {
      vfiles = aurostd::vector2deque(_vfiles);
      return true;
    }
    return false;
  }

  /// @brief returns the current directory
  /// @authors
  /// @mod{CO,20191112,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  string getPWD() {
    return std::filesystem::absolute(std::filesystem::current_path()).string();
  }

  /// @brief returns the dirname of file
  /// @param file of interest
  /// @authors
  /// @mod{CO,20210315,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  string dirname(const string& file) {
    return std::filesystem::path(aurostd::CleanFileName(file)).parent_path();
  }

  /// @brief returns the basename of file
  /// @param file of interest
  /// @authors
  /// @mod{CO,20210315,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  string basename(const string& file) {
    return std::filesystem::path(aurostd::CleanFileName(file)).stem();
  }

  /// @brief returns the file extension
  /// @authors
  /// @mod{SD,20240326,created function}
  /// @note legacy function to work with strings rather than filesystem objects directly
  string GetFileExtension(const string& FileName) {
    return std::filesystem::path(FileName).extension().string();
  }

  /// @brief returns a relative path object based on the first occurrence of a given folder name
  /// @param starting_path path where the relative search starts
  /// @param folder_name folder name to search for
  /// @note starting_path `/A/C/B/C/D/E.txt` + folder_name `C` results ib `D/E.txt`
  /// @return relative path, or empty fs::path if nothing is found
  /// @authors
  /// @mod{HE,20250717,created}
  fs::path GetRelativePath(const fs::path& starting_path, const string& folder_name) {
    fs::path current_path = fs::absolute(starting_path);
    while (current_path.has_parent_path()) {
      current_path = current_path.parent_path();
      if (current_path.filename().string() == folder_name) {
        return fs::relative(starting_path, current_path);
      }
      if (current_path == current_path.parent_path()) { // found root directory, stop search
        break;
      }
    }
    return {};
  }

  /// @brief cleans file names from obvious things
  /// @note replaces ? and * - which could be less obvious
  /// @authors
  /// @mod{SC,,created}
  /// @mod{ME,20200922,remove all //}
  /// @mod{ME,20211001,clean low ASCII characters}
  /// @mod{SD,20241126,replace ~ with HOME}
  string CleanFileName(const string& fileIN) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (fileIN.empty()) {
      return fileIN;
    }
    // ME20211001
    // Remove any control characters (below ASCII 32) while copying. This is useful
    // when fileIN is read from a file, which can have all sorts of junk and causes
    // FileExist to break. We cannot use RemoveControlCodeCharactersFromString()
    // because it keeps tabs and linebreaks and cannot use CleanStringASCII because
    // it keeps control characters.
    string fileOUT;
    for (const char c : fileIN) {
      if (c > 31) {
        fileOUT += c;
      }
    }
    if (fileOUT.find('~') != string::npos) {
      aurostd::StringSubstInPlace(fileOUT, "~", getenv("HOME"));
    }
    // ME20200922 - Cleaning // must be in a while loop, or it won't clean e.g. ///
    while (fileOUT.find("//") != string::npos) {
      aurostd::StringSubstInPlace(fileOUT, "//", "/");
    }
    aurostd::StringSubstInPlace(fileOUT, "/./", "/");

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " " << fileOUT << endl;
    }
    return fileOUT;
  }
  /// @brief overload of CleanFileName to support fs::path
  /// @note replaces ? and * - which could be less obvious
  /// @authors
  /// @mod{HE,20241217,created}
  string CleanFileName(const fs::path& fileIN) {
    return CleanFileName(fileIN.string());
  }

  /// @brief check the existence of a LOCK file
  bool DirectoryLocked(const string& directory, const string& LOCK) {
    if (FileExist(directory + "/" + LOCK)) {
      return true;
    }
    if (FileExist(directory + "/" + LOCK + _LOCK_LINK_SUFFIX_)) {
      return true;
    }
    if (FileExist(directory + "/" + LOCK + ".xz")) {
      return true;
    }
    if (FileExist(directory + "/" + LOCK + ".gz")) {
      return true;
    }
    if (FileExist(directory + "/" + LOCK + ".bz2")) {
      return true;
    }
    return false;
  }

  bool DirectoryLocked(const string& directory) {
    return DirectoryLocked(directory, "LOCK");
  }

  /// @brief check the existence of a SKIP file
  bool DirectorySkipped(const string& directory) {
    if (FileExist(directory + "/SKIP")) {
      return true;
    }
    if (FileExist(directory + "/SKIP.xz")) {
      return true;
    }
    if (FileExist(directory + "/SKIP.gz")) {
      return true;
    }
    if (FileExist(directory + "/SKIP.bz2")) {
      return true;
    }
    return false;
  }

  /// @brief create a temporary file to check if a directory is writable
  bool DirectoryWritable(const string& directory_raw) {  // "" compliant SC20190401
    const string Directory(aurostd::CleanFileName(directory_raw));
    const string filename = aurostd::TmpFileCreate("DirectoryWritable", Directory, true); // SD20220223 - uses TmpFileCreate
    const bool writable = aurostd::string2file("DirectoryWritable", filename);
    if (!writable || !FileExist(filename)) {
      return false;
    }
    aurostd::RemoveFile(filename);
    return true;
  }

  /// @brief check if path is a directory
  /// @param path to check
  /// @authors
  /// @mod{CO,20200531,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool IsDirectory(const string& path) {
    return std::filesystem::is_directory(aurostd::CleanFileName(path));
  }

  /// @brief check if a file is compressed
  /// @param FileNameIn to be checked
  /// @authors
  /// @mod{CO,20180220,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  //***************************************************************************//
  // aurostd::IsCompressed
  //***************************************************************************//
  bool IsCompressed(const string& FileNameIn) {
    const std::filesystem::path path(FileNameIn);
    archive* a;
    archive_entry* entry;
    bool check = aurostd::FileNotEmpty(FileNameIn);
    int r;

    // initialized the new archive
    a = archive_read_new();
    if (path.extension() == ".zip") {
      archive_read_support_format_zip(a);
    } else {
      archive_read_support_format_raw(a);
      archive_read_support_format_empty(a);
      archive_read_support_filter_all(a);
    }

    // open the archive
    r = archive_read_open_filename(a, FileNameIn.c_str(), 16384);
    if (r != ARCHIVE_OK) {
      check = false;
    }

    if (check) {
      // check if the content is valid (only quick check - file could be corrupt later)
      archive_read_next_header(a, &entry);
      if ((static_cast<std::string>(archive_filter_name(a, r)) == "none") && (static_cast<std::string>(archive_format_name(a)) == "raw")) {
        check = false;
      }
    }
    archive_read_free(a);
    return check;
  }

  /// @brief check if a file is compressed
  /// @param FileNameIn to be checked
  /// @param FileNameOut name of the decompressed file
  /// @authors
  /// @mod{CO,20180220,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  bool IsCompressed(const string& FileNameIn, string& FileNameOut) {
    const std::filesystem::path path(FileNameIn);
    FileNameOut = path.parent_path().string() + "/" + path.stem().string();
    return IsCompressed(FileNameIn);
  }

  /// @brief check if file exists
  /// @param FileName to check
  /// @authors
  /// @mod{SC,20080101,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool FileExist(const string& FileName) {
    try {
      return std::filesystem::exists(aurostd::CleanFileName(FileName));
    } catch (std::filesystem::filesystem_error& e) {
      return false;
    }
  }

  /// @brief check if file exists
  /// @param FileName to check
  /// @param FileNameOut overwrite with FileName if it exists
  /// @authors
  /// @mod{SC,20080101,created function}
  bool FileExist(const string& FileName, string& FileNameOut) {
    if (FileExist(FileName)) {
      FileNameOut = FileName;
      return true;
    }
    return false;
  }

  /// @brief check if a file or its compressed variant exits
  bool CompressFileExist(const string& FileName) {
    string FileNameOut;
    return CompressFileExist(FileName, FileNameOut);
  }

  /// @brief check if file or its compressed variant exits
  /// @param FileNameRaw to check
  /// @param FileNameOut overwrite with (possibly compressed) FileName if it exists
  /// @authors
  /// @mod{SC,20080101,created function}
  /// @mod{CO,20191110,add CleanFileName}
  bool CompressFileExist(const string& FileNameRaw, string& FileNameOut) {
    const string FileName = aurostd::CleanFileName(FileNameRaw);
    if (FileExist(FileName)) {
      FileNameOut = FileName;
      return true;
    }
    if (FileExist(FileName + ".xz")) {
      FileNameOut = FileName + ".xz";
      return true;
    }
    if (FileExist(FileName + ".bz2")) {
      FileNameOut = FileName + ".bz2";
      return true;
    }
    if (FileExist(FileName + ".gz")) {
      FileNameOut = FileName + ".gz";
      return true;
    }
    if (FileExist(FileName + ".zip")) {
      FileNameOut = FileName + ".zip";
      return true;
    }
    return false;
  }

  /// @brief Determines whether a compressed file is within the list
  /// @param list of files
  /// @param input filename
  /// @return boolean if the input was found in the list or not
  /// @authors
  /// @mod{CO,20200223,created function}
  bool CompressFileWithinList(const vector<string>& list, const string& input) {
    string output;
    return CompressFileWithinList(list, input, output);
  }

  /// @brief Determines whether a compressed file is within the list
  /// @param list of files
  /// @param input filename
  /// @param output filename
  /// @return boolean if the input was found in the list or not
  /// @authors
  /// @mod{CO,20200223,created function}
  /// @mod{HE,20241217,change to loop for all supported compression types}
  bool CompressFileWithinList(const vector<string>& list, const string& input, string& output) {
    output = "";
    for (const string& entry : list) {
      for (const string& suffix : compression_suffix) {
        if (entry == input + suffix) {
          output = input + suffix;
          return true;
        }
      }
    }
    return false;
  }

  /// @brief check if file is empty
  /// @param FileNameRaw to check
  /// @authors
  /// @mod{SC,20080101,created function}
  /// @mod{SD,20240312,rewritten using traits_type::eof}
  bool FileEmpty(const string& FileNameRaw) {
    const string FileName(aurostd::CleanFileName(FileNameRaw));
    if (!FileExist(FileName)) {
      return true;  // does not exist hence empty
    }
    std::ifstream FileStream(FileName);
    return (FileStream.peek() == std::ifstream::traits_type::eof());
  }

  /// @brief check if file is not empty
  bool FileNotEmpty(const string& FileName) {
    return !FileEmpty(FileName);
  }

  /// @brief gets modification time
  /// @param FileNameRaw to check
  /// @return modified timestamp (in seconds since epoch)
  /// @authors
  /// @mod{ME,20180712,created}
  /// @mod{HE,20260114,change to std::filesystem calls}
  long int GetTimestampModified(const string& FileNameRaw) {
    const fs::path p = fs::path(aurostd::CleanFileName(FileNameRaw));
    const std::filesystem::file_time_type ftime = std::filesystem::last_write_time(p);
    return std::chrono::duration_cast<std::chrono::seconds>(ftime.time_since_epoch()).count();
  }

  /// @brief gets duration since last file modification
  /// @param FileNameRaw to check
  /// @return duration in seconds
  /// @authors
  /// @mod{CO,20210315,created}
  /// @mod{ME,20210315,created}
  /// @mod{HE,20260114,change to std::filesystem calls}
  long int SecondsSinceFileModified(const string& FileNameRaw) {
    const fs::path p = fs::path(aurostd::CleanFileName(FileNameRaw));
    const auto file_ts = std::filesystem::last_write_time(p);
    const auto now_ts = std::filesystem::file_time_type::clock::now();
    const long int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now_ts - file_ts).count();
    if (elapsed_seconds > 0) {
      return elapsed_seconds;
    }
    return 0;
  }

  /// @brief gets duration since the folder or its content was last modified
  /// @param FileNameRaw to check
  /// @return duration in seconds
  /// @authors
  /// @mod{HE,20260114,created}
  long int SecondsSinceDirectoryContentModified(const string& FileNameRaw) {
    const fs::path p = fs::path(aurostd::CleanFileName(FileNameRaw));
    if (fs::is_directory(p)) {
      // use the last modification of the folder itself as a start point
      std::filesystem::file_time_type last_ts = std::filesystem::last_write_time(p);
      for (const auto& folder_entry : fs::directory_iterator(p)) {
        const std::filesystem::file_time_type ftime = std::filesystem::last_write_time(folder_entry);
        if (ftime > last_ts) {
          last_ts = ftime;
        }
      }
      const auto now_ts = std::filesystem::file_time_type::clock::now();
      const long int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now_ts - last_ts).count();
      if (elapsed_seconds > 0) {
        return elapsed_seconds;
      }
      return 0;
    }
    return SecondsSinceFileModified(FileNameRaw);
  }

  /// @brief returns in bytes the size of a file
  /// @param FileName to check
  /// @authors
  /// @mod{SC,20080101,created function}
  /// @mod{ME,20191001,changed to unsigned long long int to accommodate large files}
  /// @mod{CO,20210315,made faster}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  unsigned long long int FileSize(const string& FileName) {
    return std::filesystem::file_size(aurostd::CleanFileName(FileName));
  }

  /// @brief checks if a filestream is empty
  /// @todo should be replaces by aurostd::FileEmpty() in code and then removed
  /// @authors
  /// @mod{DM,,created}
  void InFileExistCheck(const string& routine, const string& FileNameRaw, const std::ifstream& file_to_check) {
    const string FileName(aurostd::CleanFileName(FileNameRaw));
    if (!file_to_check) {
      const string message = "In routine " + routine + ". Cannot open file " + FileName + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
  }

  /// @brief create a string pointing to the location of a unique temp file or directory
  /// @param identifier_raw overwrite the identifier (default: tmp)
  /// @param tmpdir_raw overwrite the tmpdir (default: `XHOST.tmpfs`)
  /// @param hidden hide the file
  /// @param directory create it for a directory (no suffix)
  /// @note the file or folder will not be created
  /// @return path
  string TmpStrCreate(const string& identifier_raw, const string& tmpdir_raw, const bool hidden, const bool directory) {
    string identifier = identifier_raw;
    if (identifier.empty()) {
      identifier = "tmp";
    } // CO20210624
    string tmpdir = tmpdir_raw;
    if (tmpdir.empty()) {
      tmpdir = XHOST.tmpfs;
    } // CO20210315
    string str = tmpdir
               + "/"
               + (hidden ? "." : "")
               + "_aflow_"
               + identifier
               + "."
               + XHOST.user
               + ".pid"
               + XHOST.ostrPID.str()
               + ".tid"
               + XHOST.ostrTID.str()
               + ".a"
               + AFLOW_VERSION
               + ".rnd"
               + aurostd::utype2string(uint(std::floor((double) 100000 * aurostd::ran0())))
               + ".u"
               + aurostd::utype2string(uint((double) aurostd::get_useconds()))
               + (directory ? "_" : ".")
               + "tmp"; // CO20200502 - threadID
    str = aurostd::CleanFileName(str);
    return str;
  }

  /// @brief create a string pointing to the location of a unique temp file
  /// @see
  /// @xlink{aurostd::TmpStrCreate}
  string TmpFileCreate(const string& identifier, const string& tmpdir, const bool hidden) { // CO20210315
    return TmpStrCreate(identifier, tmpdir, hidden, false); // directory=false -> creating a file
  }

  /// @brief create a string pointing to the location of a unique temp directory
  /// @note creates the folder
  /// @see
  /// @xlink{aurostd::TmpStrCreate}
  string TmpDirectoryCreate(const string& identifier, const string& tmpdir, const bool hidden) { // CO20210315
    string dir = TmpStrCreate(identifier, tmpdir, hidden, true); // directory=true -> creating a directory
    DirectoryMake(dir);
    return dir;
  }

  /// @brief copy a file
  /// @param file_from source
  /// @param file_to target
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @mod{HE,20240815,catch identical files}
  /// @mod{ST,20251208,catch errors ENODATA EIO ENOSYS EINVAL and try safer function}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool CopyFile(const string& file_from, const string& file_to) {
    const std::filesystem::path from_path = CleanFileName(file_from);
    const std::filesystem::path to_path = CleanFileName(file_to);

    if (!std::filesystem::is_regular_file(from_path)) {
      return false;
    }
    if (std::filesystem::is_regular_file(to_path)) { // equivalence check can only be performed if both paths point to a file
      if (std::filesystem::equivalent(from_path, to_path)) {
        return false;
      }
    }
    // now try to do the copy
    try {
      return std::filesystem::copy_file(from_path, to_path, std::filesystem::copy_options::overwrite_existing);
    } catch (const std::filesystem::filesystem_error& e) {
      // uh oh, check the error codes and try to use a different pathway if the error might be sendfile related
      switch (e.code().value()) {
        case ENODATA: // occurs from bug on some Lustre based file systems even where there is data
        case EIO: // in case the fs driver/implementation does not support the system call but reports so incorrectly, or there is an issue on a networked filesystem
        case ENOSYS: // suggested by sendfile(2) man pages, this would mean the system call does not exist
        case EINVAL: // suggested by sendfile(2) man pages, this would mean the kernel does not want to use the system call on our files
          return CopyFile_safer(from_path, to_path);
        default: return false;
      }
    }
  }

  /// @brief alternative copy file function that tries to use a user-space pathway
  /// This function uses iostream system to try to copy via a user-space read/write
  /// rather than the sendfile kernel-space pathway that std::filesystem::copy_file will use.
  /// While sendfile is theoretically faster, it might choke if there are cache or network
  /// issues with the supporting filesystem. Read and write are more stable system calls.
  /// While we don't use read/write syscalls here explicitly, the streams should pass the
  /// data through user-space and compile down to read/write syscalls.
  /// @param file_from source
  /// @param file_to target
  /// @authors
  /// @mod{ST,20251208,created function}
  /// @note Do not use this, use CopyFile instead
  bool CopyFile_safer(const std::string& file_from, const std::string& file_to) {
    std::ifstream in_file{file_from, std::ios::in | std::ios::binary};
    if (!in_file.is_open()) {
      return false;
    }
    std::ofstream out_file{file_to, std::ios::out | std::ios::binary};
    if (!out_file.is_open()) {
      in_file.close();
      return false;
    }

    out_file << in_file.rdbuf();
    if (!out_file.good()) {
      in_file.close();
      out_file.close();
      return false;
    }

    try {
      std::filesystem::permissions(file_to, std::filesystem::status(file_from).permissions());
    } catch (...) {
      // do nothing
    }

    in_file.close();
    out_file.close();
    return true;
  }

  /// @brief copy a file
  /// @param directory_from source
  /// @param directory_to target
  /// @authors
  /// @mod{HE,20260129,created function}
  /// @note legacy function to work with strings rather than filesystem objects directly
  /// @note enables workarounds for problematic filesystems
  bool CopyDirectory(const std::string& directory_from, const std::string& directory_to) {
    const std::filesystem::path from_path = CleanFileName(directory_from);
    const std::filesystem::path to_path = CleanFileName(directory_to);

    if (!std::filesystem::is_directory(from_path)) {
      return false;
    }

    if (std::filesystem::is_regular_file(to_path)) {
      return false;
    }

    if (std::filesystem::is_directory(to_path)) {
      if (std::filesystem::equivalent(from_path, to_path)) {
        return false;
      }
    }
    try {
      fs::create_directories(to_path.parent_path());
      fs::copy(from_path, to_path, fs::copy_options::recursive);
      return true;
    } catch (const std::filesystem::filesystem_error& e) {
      // uh oh, check the error codes and try to use a different pathway if the error might be sendfile related
      switch (e.code().value()) {
        case ENODATA: // occurs from bug on some Lustre based file systems even where there is data
        case EIO: // in case the fs driver/implementation does not support the system call but reports so incorrectly, or there is an issue on a networked filesystem
        case ENOSYS: // suggested by sendfile(2) man pages, this would mean the system call does not exist
        case EINVAL: // suggested by sendfile(2) man pages, this would mean the kernel does not want to use the system call on our files
          return CopyDirectory_safer(from_path, to_path);
        default: return false;
      }
    }
  }

  /// @brief alternative copy directory function that uses CopyFile_safer
  /// This function similar to CopyFile_safer is a fallback option for problematic filesystems.
  /// It will copy only regular files and directories recursively.
  /// @param directory_from source
  /// @param directory_to target
  /// @authors
  /// @mod{HE,20260129,created function}
  /// @note Do not use this, use CopyDirectory instead
  bool CopyDirectory_safer(const std::string& directory_from, const std::string& directory_to) {
    const std::filesystem::path from_path = CleanFileName(directory_from);
    const std::filesystem::path to_path = CleanFileName(directory_to);
    fs::create_directories(to_path.parent_path()); // ensure parent exist
    fs::create_directory(to_path, from_path); // make new folder with the same attributes as old one
    for (const std::filesystem::directory_entry& x : std::filesystem::directory_iterator(from_path)) {
      if (x.is_regular_file()) {
        if (!CopyFile(x.path(), directory_to / x.path().filename())) {
          return false;
        }
      } else if (x.is_directory()) {
        if (!CopyDirectory_safer(x.path(), directory_to / x.path().filename())) {
          return false;
        }
      }
    }
    return true;
  }

  /// @brief link a file
  /// @param file_from source
  /// @param file_to target
  /// @param soft create a soft-link instead of a hard-link
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool LinkFile(const string& file_from, const string& file_to, const bool soft) {
    const string from_clean = CleanFileName(file_from);
    const string to_clean = CleanFileName(file_to);
    const std::filesystem::path from_path = from_clean;
    const std::filesystem::path to_path = to_clean;
    if (from_path.empty() || to_path.empty() || !std::filesystem::exists(from_path) || std::filesystem::exists(to_path) || !std::filesystem::exists(to_path.parent_path())) {
      return false;
    }
    if (soft) {
      std::filesystem::create_symlink(from_clean, to_clean);
      return std::filesystem::is_symlink(to_clean);
    }
    std::filesystem::create_hard_link(from_clean, to_clean);
    return std::filesystem::is_regular_file(to_clean);
  }

  /// @brief move file to directory
  /// @param file to move
  /// @param destination name of directory
  /// @authors
  /// @mod{CO,20171025,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool file2directory(const string& file, const string& destination) {
    const std::filesystem::path old_path = CleanFileName(file);
    const std::filesystem::path dir_path = CleanFileName(destination);
    if (!std::filesystem::is_regular_file(old_path)) {
      return false;
    }
    if (!std::filesystem::is_directory(dir_path)) {
      return false;
    }
    try {
      std::filesystem::rename(old_path, dir_path / old_path.filename());
    } catch (std::filesystem::filesystem_error& e) {
      std::filesystem::copy(old_path, dir_path / old_path.filename(), std::filesystem::copy_options::overwrite_existing);
      std::filesystem::remove(old_path);
    }
    return true;
  }

  /// @brief move a list of files to directory
  /// @see
  /// @xlink{aurostd::file2directory}
  bool file2directory(const vector<string>& files, const string& destination) {
    for (const auto& file : files) {
      if (!file2directory(file, destination)) {
        return false;
      }
    }

    return true;
  }

  // ***************************************************************************
  // Function SplitFileName
  // ***************************************************************************
  /// @brief WARNING: this function is deprecated; users should use the std::filesystem::path functions instead
  /// @param filePath
  /// @return std::array<std::string,3> {parent, name, suffix}
  std::array<std::string, 3> splitFilePath(const string& filePath) {
    size_t type_sep = filePath.rfind(".");
    size_t folder_sep = filePath.rfind("/");
    std::string suffix;
    std::string name;
    std::string parent;
    if (folder_sep != std::string::npos) {
      parent = filePath.substr(0, folder_sep + 1);
      if (type_sep != std::string::npos) {
        if (type_sep > folder_sep) {
          suffix = filePath.substr(type_sep);
          name = filePath.substr(folder_sep + 1, type_sep - folder_sep - 1);
        } else {
          name = filePath.substr(folder_sep + 1);
        }
      } else {
        name = filePath.substr(folder_sep + 1);
      }
    } else {
      if (type_sep != std::string::npos) {
        suffix = filePath.substr(type_sep);
        name = filePath.substr(0, type_sep);
      } else {
        name = filePath;
      }
    }
    return {parent, name, suffix};
  }

  /// @brief move file to another file
  /// @param file to move
  /// @param destination name of file
  /// @authors
  /// @mod{CO,20171025,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool file2file(const string& file, const string& destination) {
    const std::filesystem::path old_path = CleanFileName(file);
    const std::filesystem::path new_path = CleanFileName(destination);
    if (!std::filesystem::is_regular_file(old_path)) {
      return false;
    }
    try {
      std::filesystem::rename(old_path, new_path);
    } catch (std::filesystem::filesystem_error& e) {
      std::filesystem::copy(old_path, new_path, std::filesystem::copy_options::overwrite_existing);
      std::filesystem::remove(old_path);
    }
    return true;
  }

  /// @brief remove a file
  /// @param file to be removed
  /// @authors
  /// @mod{SC,20080101,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool RemoveFile(const string& file) {
    const std::filesystem::path path(aurostd::CleanFileName(file));
    if (std::filesystem::is_regular_file(path)) {
      return std::filesystem::remove(path);
    }
    return false;
  }

  /// @brief remove files that match a regex pattern
  /// @param directory to search
  /// @param regex_pattern to match
  /// @authors
  /// @mod{SD,20240312,created function}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool RemoveFile(const string& directory, const std::regex& regex_pattern) {
    const std::filesystem::path path(aurostd::CleanFileName(directory));
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(std::filesystem::path(aurostd::CleanFileName(directory)))) {
      if (dir_entry.is_regular_file() && std::regex_match(dir_entry.path().filename().string(), regex_pattern)) {
        std::filesystem::remove(dir_entry.path());
      }
    }
    return true;
  }

  /// @brief removes all regular files inside a directory, not recursive
  /// @param directory the directory in which files are removed
  /// @authors
  /// @mod{ST,20241016,created}
  bool RemoveFiles(const string& directory) {
    for (const auto& dir_entry : std::filesystem::directory_iterator(aurostd::CleanFileName(directory))) {
      if (dir_entry.is_regular_file()) {
        std::filesystem::remove(dir_entry.path());
      }
    }
    return true;
  }

  /// @brief remove a list of files to directory
  /// @see
  /// @xlink{aurostd::RemoveFile}
  bool RemoveFile(const vector<string>& files) {
    for (const auto& file : files) {
      if (!RemoveFile(file)) {
        return false;
      }
    }
    return true;
  }

  /// @brief remove a file
  /// @param directory to be removed
  /// @authors
  /// @mod{CO,20171025,created function}
  /// @mod{SD,20240312,rewritten using filesystem}
  /// @note legacy function to work with strings rather than filesystem objects directly
  bool RemoveDirectory(const string& directory) {
    const std::filesystem::path path(aurostd::CleanFileName(directory));
    if (std::filesystem::is_directory(path)) {
      std::filesystem::remove_all(path);
    }
    return false;
  }

  compression_type GetCompressionType(const fs::path& path) {
    const string extension = aurostd::GetFileExtension(path);
    const auto it = std::find(compression_suffix.begin(), compression_suffix.end(), extension);
    if (it == compression_suffix.end()) {
      return compression_type::None;
    }
    const compression_type ct = static_cast<compression_type>(it - compression_suffix.begin());
    return ct;
  }

  /// @brief read the content of a file into a string
  /// @authors
  /// @mod{SC,,created}
  /// @mod{HE,20220221,avoids reading the file char by char}
  /// @mod{HE,20241216,add transparent compression support}
  size_t file2string(const string& FileNameIN, string& StringIN) {
    const string FileNameINClean = aurostd::CleanFileName(FileNameIN);
    if (GetCompressionType(FileNameINClean)) {
      return compressfile2string(FileNameINClean, StringIN);
    }
    const std::ifstream open_file(FileNameINClean);
    std::stringstream buffer;
    buffer << open_file.rdbuf();
    StringIN = buffer.str();
    return StringIN.length();
  }

  /// @brief return content file content as string
  string file2string(const string& FileNameIN) {
    string StringIN;
    file2string(FileNameIN, StringIN);
    return StringIN;
  }

  /// @brief compressed file to string
  /// @param FileNameIN to be read
  /// @param content of the file
  /// @authors
  /// @mod{CO,20210624,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  size_t compressfile2string(const string& FileNameIN, string& content) {
    const string file(aurostd::CleanFileName(FileNameIN));
    std::stringstream ss;
    compressfile2stringstream(file, ss);
    content = ss.str();
    return content.length();
  }

  /// @brief return content of a compressed file as string
  string compressfile2string(const string& FileNameIN) {
    string StringIN;
    compressfile2string(FileNameIN, StringIN);
    return StringIN;
  }

  /// @brief write the content of a file into a vector of strings
  size_t file2vectorstring(const string& FileNameIN, vector<string>& vline, bool consecutive, bool trim_edges) {
    return aurostd::string2vectorstring(file2string(aurostd::CleanFileName(FileNameIN)), vline, consecutive, trim_edges);
  }

  /// @brief write the content of a file into a deque of strings
  size_t file2vectorstring(const string& FileNameIN, deque<string>& vline) {
    return aurostd::string2dequestring(file2string(aurostd::CleanFileName(FileNameIN)), vline);
  }

  /// @brief write the content of a compressed file into a vector of strings
  /// @todo check usage - file2vectorstring should cover this use case now
  size_t compressfile2vectorstring(const string& FileNameIN, vector<string>& vline, bool consecutive, bool trim_edges) {
    return aurostd::string2vectorstring(compressfile2string(aurostd::CleanFileName(FileNameIN)), vline, consecutive, trim_edges);
  }

  /// @brief write the content of a compressed file into a deque of strings
  /// @todo check usage - file2vectorstring should cover this use case now
  size_t compressfile2vectorstring(const string& FileNameIN, deque<string>& vline) {
    return aurostd::string2dequestring(compressfile2string(aurostd::CleanFileName(FileNameIN)), vline);
  }

  /// @brief write the content of a file into a stringstream
  /// @mod{HE,20241216,add transparent compression support}
  void file2stringstream(const string& FileNameIN, stringstream& StringstreamIN) {
    const string FileNameINClean = aurostd::CleanFileName(FileNameIN);
    if (GetCompressionType(FileNameINClean)) {
      compressfile2stringstream(FileNameINClean, StringstreamIN);
    } else {
      const std::ifstream open_file(FileNameINClean);
      StringstreamIN << open_file.rdbuf();
    }
  }

  /// @brief write the content of a compressed file into a stringstream
  /// @param file_raw to be read
  /// @param content of the file
  /// @authors
  /// @mod{CO,20210624,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  /// @mod{HE,20250728,add auto suffix detection}
  void compressfile2stringstream(const string& file_raw, stringstream& content) {
    string file(aurostd::CleanFileName(file_raw));
    // find the actual file
    for (const std::string& suffix : compression_suffix) {
      if (const string file_compressed = file + suffix; std::filesystem::exists(file_compressed)) {
        file = file_compressed;
        break;
      }
    }

    static constexpr size_t buff_len = 16384;
    thread_local char buff[buff_len + 1];
    archive* a = archive_read_new();
    archive_entry* entry;

    // initialized the archive to be read

    if (aurostd::GetFileExtension(file) == ".zip") {
      archive_read_support_format_zip(a);
    } else {
      archive_read_support_format_raw(a);
      archive_read_support_format_empty(a);
      archive_read_support_filter_all(a);
    }

    // open the archive
    if (const int r = archive_read_open_filename(a, file.c_str(), buff_len); r != ARCHIVE_OK) {
      const std::string message = "Problem loading output: ";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message + archive_error_string(a), _FILE_NOT_FOUND_);
    }
    archive_read_next_header(a, &entry);

    // prepare the stringstream and write the data
    aurostd::StringstreamClean(content);
    size_t len = archive_read_data(a, buff, buff_len);
    while (len > 0) { // read until archive is empty
      buff[len] = '\0'; // ensure proper termination of the c style string
      content << buff;
      len = archive_read_data(a, buff, buff_len);
    }

    // free objects created in this function
    archive_read_free(a);
  }

  /// @brief write string into a file
  /// @param StringOUTPUT content to write
  /// @param FileNameRawOUTPUT file to write to
  /// @param ct compression type (default None)
  /// @param mode how to write the data (APPEND, PRE, WRITE) - default WRITE
  /// @return success
  /// @note this functions ensures that the full path is available
  /// @authors
  /// @mod{SC,20190401,created}
  /// @mod{HE,20241216,add compression}
  bool string2file(const string& StringOUTPUT, const std::string& FileNameRawOUTPUT, const compression_type ct, const string& mode) {
    const fs::path FileNameOUTPUT = CleanFileName(FileNameRawOUTPUT);
    try {
      fs::create_directories(fs::absolute(FileNameOUTPUT).parent_path());
    } catch (fs::filesystem_error& e) {
      return false;
    }
    bool writable = true; // CO20190808 - captures whether we can open/write file

    if (mode == "POST" || mode == "APPEND" || mode == "PRE") {
      string FileINPUT;
      if (fs::exists(FileNameOUTPUT)) {
        aurostd::file2string(FileNameOUTPUT, FileINPUT);
      }
      if (ct) {
        if (mode == "POST" || mode == "APPEND") {
          return string2compressfile(StringOUTPUT + FileINPUT, FileNameOUTPUT, ct);
        }
        return string2compressfile(FileINPUT + StringOUTPUT, FileNameOUTPUT, ct);
      }
      std::ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(), std::ios::out);
      writable = FileOUTPUT.is_open(); // CO20190808 - captures whether we can open/write file
      if (mode == "POST" || mode == "APPEND") {
        FileOUTPUT << FileINPUT;
      }
      FileOUTPUT << StringOUTPUT;
      if (mode == "PRE") {
        FileOUTPUT << FileINPUT;
      }
      FileOUTPUT.flush();
      FileOUTPUT.clear();
      FileOUTPUT.close();
      return writable; // return false if something got messed up //CO20190808 - captures whether we can open/write file
    }

    if (mode == "WRITE" || mode.empty()) {
      if (ct) {
        return string2compressfile(StringOUTPUT, FileNameOUTPUT, ct);
      }
      std::ofstream FileOUTPUT;
      FileOUTPUT.open(FileNameOUTPUT.c_str(), std::ios::out);
      writable = FileOUTPUT.is_open(); // CO20190808 - captures whether we can open/write file
      FileOUTPUT << StringOUTPUT;
      FileOUTPUT.flush();
      FileOUTPUT.clear();
      FileOUTPUT.close();
      return writable; // return false if something got messed up //CO20190808 - captures whether we can open/write file
    }
    return false;
  }

  /// @brief string to compressed file
  /// @param content of the string
  /// @param file_raw file to save data to
  /// @param ct compression type
  /// @note this function creates the parent folder if necessary
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  bool string2compressfile(const string& content, const string& file_raw, const compression_type ct) {
    fs::path file(aurostd::CleanFileName(file_raw));
    try {
      fs::create_directories(fs::absolute(file).parent_path());
    } catch (fs::filesystem_error& e) {
      return false;
    }

    archive* a;
    archive_entry* entry;
    int r;

    // initialized the new archive
    a = archive_write_new();
    entry = archive_entry_new();

    // select compression and archive format
    aurostd::CompressInit(a, ct, true);

    // open the new archive
    if (const std::string ext = file.extension(); ext != compression_suffix[ct]) {
      file = file.replace_extension(ext + compression_suffix[ct]);
    }

    r = archive_write_open_filename(a, file.c_str());
    if (r != ARCHIVE_OK) {
      const std::string message = "Problem preparing output: ";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message + archive_error_string(a), _FILE_NOT_FOUND_);
    }

    // configure new entry
    archive_entry_set_pathname(entry, file.c_str());
    archive_entry_set_size(entry, content.size());
    archive_entry_set_filetype(entry, AE_IFREG);

    // save file information to archive
    archive_write_header(a, entry);

    // write string to archive
    archive_write_data(a, content.c_str(), content.size());

    // free objects created in this function
    archive_entry_free(entry);
    archive_write_close(a);
    archive_write_free(a);
    return true;
  }

  /// @brief write stringstream into a file
  bool stringstream2file(const stringstream& StringstreamOUTPUT, const string& FileNameOUTPUT, const compression_type ct, const string& mode) {
    return string2file(StringstreamOUTPUT.str(), FileNameOUTPUT, ct, mode);
  }

  /// @brief stringstream to compressed file
  /// @param content of the string
  /// @param file name
  /// @param ct compression type
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  bool stringstream2compressfile(const stringstream& content, const string& file, const compression_type ct) {
    return string2compressfile(content.str(), file, ct);
  }

  /// @brief write a list of strings to file
  /// @tparam vtype container type (vector, deque)
  /// @param vline container of lists
  /// @param FileNameOUT name of file to write
  /// @note this function creates the parent folder if necessary
  /// @return success
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{HE,20240328,changed to template approach}
  template <template <class...> class vtype> bool vectorstring2file(const vtype<string>& vline, const string& FileNameOUT, const compression_type ct) {
    const fs::path file = aurostd::CleanFileName(FileNameOUT);

    try {
      fs::create_directories(fs::absolute(file).parent_path());
    } catch (fs::filesystem_error& e) {
      return false;
    }

    if (ct) {
      std::stringstream full_string;
      for (const auto& iline : vline) {
        full_string << iline << endl;
      }
      return stringstream2compressfile(full_string, file, ct);
    }
    std::ofstream FileOUT;
    FileOUT.open(file.c_str(), std::ios::out);
    const bool writable = FileOUT.is_open(); // CO20190808 - captures whether we can open/write file
    for (const auto& iline : vline) {
      FileOUT << iline << endl;
    }
    FileOUT.flush();
    FileOUT.clear();
    FileOUT.close();
    return writable;
  }
  template bool vectorstring2file(const vector<string>& vline, const string& FileNameOUT, compression_type ct);
  template bool vectorstring2file(const deque<string>& vline, const string& FileNameOUT, compression_type ct);

  ///@brief check compression calls
  void CompressionErrorHandling(archive* a, const int r, const std::string& message) {
    if (r != ARCHIVE_OK) {
      if (r == ARCHIVE_WARN) {
        const string full_message = "LibArchive Warning during " + message + ": ";
        if (XHOST.DEBUG) {
          cerr << full_message << archive_error_string(a) << " (not a fatal error)" << endl;
        }
      } else {
        const string full_message = "LibArchive Error during " + message + ": ";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, full_message + archive_error_string(a), _RUNTIME_INIT_);
      }
    }
  }

  /// @brief initialize compression
  /// @param a archive
  /// @param ct compression type (default zt)
  /// @param single if true only compress a single file directly without any folder structure
  /// @authors
  /// @mod{SD,20240331,created function}
  /// @mod{HE,20241217,change to use compression_type with}
  void CompressInit(archive* a, const compression_type ct, const bool single) {
    int r = ARCHIVE_OK;
    const std::string error_message = "initialization";
    switch (ct) {
      case BZ2:
        if (single) {
          r = archive_write_set_format_raw(a);
        } else {
          r = archive_write_set_format_pax_restricted(a);
        }
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_add_filter_bzip2(a);
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_set_options(a, "compression-level=9");
        CompressionErrorHandling(a, r, error_message);
        break;
      case GZ:
        if (single) {
          r = archive_write_set_format_raw(a);
        } else {
          r = archive_write_set_format_pax_restricted(a);
        }
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_add_filter_gzip(a);
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_set_options(a, "compression-level=9");
        CompressionErrorHandling(a, r, error_message);
        break;
      case XZ:
        if (single) {
          r = archive_write_set_format_raw(a);
        } else {
          r = archive_write_set_format_pax_restricted(a);
        }
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_add_filter_xz(a);
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_set_options(a, "compression-level=9,threads=1");
        CompressionErrorHandling(a, r, error_message);
        break;
      case ZIP:
        r = archive_write_set_format_zip(a);
        CompressionErrorHandling(a, r, error_message);
        if (single) {
          r = archive_write_set_options(a, "compression=deflate,zip64");
        } else {
          r = archive_write_set_options(a, "compression=store,zip64");
        }
        CompressionErrorHandling(a, r, error_message);
        break;
      case ZSTD:
        if (single) {
          r = archive_write_set_format_raw(a);
        } else {
          r = archive_write_set_format_pax_restricted(a);
        }
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_add_filter_zstd(a);
        CompressionErrorHandling(a, r, error_message);
        r = archive_write_set_options(a, "compression-level=20");
        CompressionErrorHandling(a, r, error_message);
        break;
      default: throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown archive format", _INPUT_UNKNOWN_);
    }
  }

  /// @brief decompress file
  /// @param file_raw to be decompressed
  /// @param outfile_raw location to save the content of the compressed file
  /// @param keep the compressed file
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240325,rewritten using libarchive}
  void DecompressFile(const string& file_raw, const string& outfile_raw, const bool keep) {
    const string file(aurostd::CleanFileName(file_raw));
    string outfile(aurostd::CleanFileName(outfile_raw));
    const std::filesystem::path path(file);
    thread_local char buff[16384];
    archive* a;
    archive_entry* entry;
    size_t len;
    int r = ARCHIVE_OK;
    const std::string error_message = "decompression of a single file";

    // initialized the new archive and the file entries
    a = archive_read_new();
    if (path.extension() == ".zip") {
      r = archive_read_support_format_zip(a);
      CompressionErrorHandling(a, r, error_message);
    } else {
      r = archive_read_support_format_raw(a);
      CompressionErrorHandling(a, r, error_message);
      r = archive_read_support_format_empty(a);
      CompressionErrorHandling(a, r, error_message);
      r = archive_read_support_filter_all(a);
      CompressionErrorHandling(a, r, error_message);
    }

    // open the archive
    r = archive_read_open_filename(a, file.c_str(), sizeof(buff));
    CompressionErrorHandling(a, r, error_message);
    if (outfile.empty()) {
      outfile = path.parent_path().string() + "/" + path.stem().string();
    }
    r = archive_read_next_header(a, &entry);
    if (r == ARCHIVE_EOF) { // create empty file if archive has no data
      std::fstream fs(outfile, std::ios::out);
      fs.close();
    } else {
      CompressionErrorHandling(a, r, error_message);
      // open the new file and write it in chunks
      std::fstream fs(outfile, std::ios::out);
      len = archive_read_data(a, buff, sizeof(buff));
      while (len > 0) { // read until archive is empty
        fs.write(buff, len);
        len = archive_read_data(a, buff, sizeof(buff));
      }
      fs.close();
    }

    // free objects created in this function
    if (!keep) {
      aurostd::RemoveFile(file);
    }
    r = archive_read_free(a);
    CompressionErrorHandling(a, r, error_message);
  }

  void DecompressFile(const string& file_raw, const bool keep) {
    DecompressFile(file_raw, "", keep);
  }

  /// @brief decompress files
  /// @param archive_file to be decompressed
  /// @param directory to decompress files
  /// @param keep the compressed file
  /// @authors
  /// @mod{SD,20240403,created function}
  void DecompressFiles(const string& archive_file, const string& directory, const bool keep) {
    const string archive_file_clean(aurostd::CleanFileName(archive_file));
    string directory_clean(aurostd::CleanFileName(directory));
    std::filesystem::path path(archive_file_clean);
    thread_local char buff[16384];
    archive* a;
    archive_entry* entry;
    size_t len;
    int r = ARCHIVE_OK;
    const std::string error_message = "decompression of multiple files";

    // initialized the new archive and the file entries
    a = archive_read_new();
    if (path.extension() == ".zip") {
      r = archive_read_support_format_zip(a);
      CompressionErrorHandling(a, r, error_message);
    } else {
      r = archive_read_support_format_tar(a);
      CompressionErrorHandling(a, r, error_message);
      r = archive_read_support_filter_all(a);
      CompressionErrorHandling(a, r, error_message);
    }

    // open the archive
    r = archive_read_open_filename(a, archive_file_clean.c_str(), sizeof(buff));
    CompressionErrorHandling(a, r, error_message);

    // open the new files and write it in chunks
    if (directory_clean.empty()) {
      directory_clean = aurostd::getPWD();
    }
    while (archive_read_next_header(a, &entry) == ARCHIVE_OK) {
      if ((static_cast<std::string>(archive_filter_name(a, r)) == "none") && (static_cast<std::string>(archive_format_name(a)) == "raw")) {
        continue;
      }
      path = std::filesystem::path(directory_clean + "/" + archive_entry_pathname(entry));
      if (!std::filesystem::is_directory(path.parent_path())) {
        std::filesystem::create_directories(path.parent_path());
      }
      std::fstream fs(path.string(), std::ios::out);
      len = archive_read_data(a, buff, sizeof(buff));
      while (len > 0) { // read until archive is empty
        fs.write(buff, len);
        len = archive_read_data(a, buff, sizeof(buff));
      }
      fs.close();
    }

    // free objects created in this function
    if (!keep) {
      aurostd::RemoveFile(archive_file_clean);
    }
    r = archive_read_free(a);
    CompressionErrorHandling(a, r, error_message);
  }

  void DecompressFiles(const string& archive, const bool keep) {
    DecompressFiles(archive, "", keep);
  }

  /// @brief compress file
  /// @param file_raw to be compressed
  /// @param ct compression type
  /// @param keep the file
  /// @authors
  /// @mod{SC,20190401,created function}
  /// @mod{SD,20240325,rewritten using libarchive}
  void CompressFile(const string& file_raw, const compression_type ct, const bool keep) {
    const fs::path file(aurostd::CleanFileName(file_raw));

    thread_local char buff[16384];
    archive* a;
    archive* disk;
    archive_entry* entry;
    int r = ARCHIVE_OK;
    const std::string error_message = "compression of a single file";

    // initialized the new archive, the disk reader and the file entries
    a = archive_write_new();
    entry = archive_entry_new();
    disk = archive_read_disk_new();
    r = archive_read_disk_set_standard_lookup(disk);
    CompressionErrorHandling(a, r, error_message);

    // select compression and archive format
    aurostd::CompressInit(a, ct, true);

    // open the disk reader
    r = archive_read_disk_open(disk, file.c_str());
    CompressionErrorHandling(a, r, error_message);

    // open the new archive
    fs::path outfile = file;
    outfile.replace_extension(static_cast<string>(file.extension()) + compression_suffix[ct]);

    r = archive_write_open_filename(a, outfile.c_str());
    CompressionErrorHandling(a, r, error_message);

    // load information of file from disk
    r = archive_read_next_header2(disk, entry);
    CompressionErrorHandling(a, r, error_message);
    r = archive_read_disk_descend(disk);
    CompressionErrorHandling(a, r, error_message);

    // save file information to archive
    r = archive_write_header(a, entry);
    CompressionErrorHandling(a, r, error_message);

    // open the input file and write it in chunks to new archive
    std::ifstream is(archive_entry_sourcepath(entry), std::ifstream::binary);
    while (is) {
      is.read(buff, sizeof(buff));
      archive_write_data(a, buff, is.gcount());
    }
    is.close();

    // free objects created in this function
    if (!keep) {
      aurostd::RemoveFile(file);
    }
    archive_entry_free(entry);
    r = archive_read_free(disk);
    CompressionErrorHandling(a, r, error_message);
    r = archive_write_free(a);
    CompressionErrorHandling(a, r, error_message);
  }

  /// @brief compress files
  /// @param files_raw files to be compressed
  /// @param relative_to_path base path to construct filestructure in the archive
  /// @param archive_name of the archive
  /// @param ct compression type
  /// @param keep the files
  /// @authors
  /// @mod{SD,20240331,created function}
  /// @mod{HE,20240908,adding relative paths}
  fs::path CompressFiles(const vector<string>& files_raw, const fs::path& relative_to_path, const string& archive_name, const compression_type ct, const bool keep) {
    thread_local char buff[16384];
    archive* a;
    archive* disk;
    archive_entry* entry;
    int r = ARCHIVE_OK;
    const std::string error_message = "compression of multiple files";
    vector<string> files;
    files.reserve(files_raw.size());
    for (const auto& f : files_raw) {
      files.push_back(aurostd::CleanFileName(f));
    }

    // initialized the new archive and the disk reader
    a = archive_write_new();

    // select compression and archive format
    aurostd::CompressInit(a, ct, false);

    // open the new archive
    string tar_suffix = ".tar";
    if (ct == compression_type::ZIP) {
      tar_suffix = "";
    }
    fs::path outfile = archive_name + tar_suffix + compression_suffix[ct];
    try {
      fs::create_directories(fs::absolute(outfile).parent_path());
    } catch (fs::filesystem_error& e) {
      const string error = e.what();
      const string message = "Parent folder could not be created: " + static_cast<string>(outfile.parent_path()) + "\n" + error;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    r = archive_write_open_filename(a, outfile.c_str());
    CompressionErrorHandling(a, r, error_message);

    for (const auto& file : files) {
      // initialized the file entries
      entry = archive_entry_new();

      // initialized the disk reader
      disk = archive_read_disk_new();
      r = archive_read_disk_set_standard_lookup(disk);
      CompressionErrorHandling(a, r, error_message);

      // open the disk reader
      r = archive_read_disk_open(disk, file.c_str());
      CompressionErrorHandling(a, r, error_message);

      // load information of file from disk
      r = archive_read_next_header2(disk, entry);
      CompressionErrorHandling(a, r, error_message);
      r = archive_read_disk_descend(disk);
      CompressionErrorHandling(a, r, error_message);

      // HE20240907 change path information to be relative
      const std::string original_path_raw = archive_entry_pathname(entry);
      const fs::path original_path(original_path_raw);
      // using proximate to avoid empty paths in case `relative_to_path` and `original_path` are not compatible
      const fs::path new_path = fs::proximate(original_path, relative_to_path);
      const std::string new_path_raw = new_path.string();
      archive_entry_set_pathname(entry, new_path_raw.c_str());

      // save file information to archive
      r = archive_write_header(a, entry);
      CompressionErrorHandling(a, r, error_message);

      // open the input file and write it in chunks to new archive
      std::ifstream is(archive_entry_sourcepath(entry), std::ifstream::binary);
      while (is) {
        is.read(buff, sizeof(buff));
        archive_write_data(a, buff, is.gcount());
      }
      is.close();

      // free read objects
      if (!keep) {
        aurostd::RemoveFile(file);
      }
      archive_entry_free(entry);
      r = archive_read_free(disk);
      CompressionErrorHandling(a, r, error_message);
    }

    // free write objects
    r = archive_write_close(a);
    CompressionErrorHandling(a, r, error_message);
    r = archive_write_free(a);
    CompressionErrorHandling(a, r, error_message);

    return outfile;
  }

  /// @brief compress a file based on the compression of another filed
  /// @param CompressedFileName compressed filed
  /// @param FileNameOUT file to compress
  /// @authors
  /// @mod{CO,20190401,created function}
  /// @mod{SD,20240326,updated using filesystem}
  bool MatchCompressed(const string& CompressedFileName, const string& FileNameOUT) {
    if (!aurostd::IsCompressed(CompressedFileName) && !aurostd::IsCompressed(FileNameOUT)) {
      return true;
    }
    if (aurostd::IsCompressed(FileNameOUT)) {
      if (aurostd::GetFileExtension(CompressedFileName) == aurostd::GetFileExtension(FileNameOUT)) {
        return true;
      }
      return false;
    }
    const compression_type ct = aurostd::GetCompressionType(CompressedFileName);
    if (ct == compression_type::None) {
      return false;
    }
    aurostd::CompressFile(FileNameOUT, ct);
    return true;
  }

  /// @brief decompress a file to a temp file
  /// @param FileNameIN to be decompressed
  /// @param FileNameOUT name of temp file
  /// @authors
  /// @mod{CO,20210623,created function}
  /// @mod{SD,20240326,rewritten using libarchive}
  bool compressfile2tempfile(const string& FileNameIN, string& FileNameOUT) {
    bool tempfile_created = false;
    return compressfile2tempfile(FileNameIN, FileNameOUT, tempfile_created);
  }

  /// @brief decompress a file to a specific temp file
  bool compressfile2tempfile(const string& FileNameRawIN, string& FileNameOUT, bool& tempfile_created) {
    string FileNameIN;
    tempfile_created = false;
    if (!(aurostd::FileExist(FileNameRawIN, FileNameIN) || aurostd::CompressFileExist(FileNameRawIN, FileNameIN))) {
      return false;
    }
    if (!aurostd::IsCompressed(FileNameIN)) {
      FileNameOUT = FileNameIN;
      return true;
    }
    tempfile_created = true;
    FileNameOUT = aurostd::TmpStrCreate();
    aurostd::DecompressFile(FileNameIN, FileNameOUT, true);
    return true;
  }

  /// @brief convert a compressed file to another compression
  /// @param file to convert
  /// @param ct compression type
  /// @param keep the original compressed file
  /// @authors
  /// @mod{SD,20240329,created function}
  void compressfile2compressfile(const string& file, const compression_type ct, const bool keep) {
    const std::filesystem::path path(file);
    thread_local char buff[16384];
    archive* a;
    archive* a_new;
    archive_entry* entry;
    size_t len;
    int r = ARCHIVE_OK;
    const std::string error_message = "re-compressing of a single file";

    // initialized the new archive, the disk and the file entries
    a = archive_read_new();
    if (path.extension() == ".zip") {
      r = archive_read_support_format_zip(a);
      CompressionErrorHandling(a, r, error_message);
    } else {
      r = archive_read_support_format_raw(a);
      CompressionErrorHandling(a, r, error_message);
      r = archive_read_support_format_empty(a);
      CompressionErrorHandling(a, r, error_message);
      r = archive_read_support_filter_all(a);
      CompressionErrorHandling(a, r, error_message);
    }
    a_new = archive_write_new();

    // select compression and archive format
    aurostd::CompressInit(a_new, ct, true);

    // open the archive
    r = archive_read_open_filename(a, file.c_str(), sizeof(buff));
    CompressionErrorHandling(a, r, error_message);
    r = archive_read_next_header(a, &entry);
    CompressionErrorHandling(a, r, error_message);

    // open the new archive
    const string outfile = path.parent_path().string() + "/" + path.stem().string() + compression_suffix[ct];
    r = archive_write_open_filename(a_new, outfile.c_str());
    CompressionErrorHandling(a, r, error_message);

    // save file information to archive
    r = archive_write_header(a_new, entry);
    CompressionErrorHandling(a, r, error_message);

    // open the input archive and write it in chunks to new archive
    len = archive_read_data(a, buff, sizeof(buff));
    while (len > 0) { // read until archive is empty
      archive_write_data(a_new, buff, len);
      len = archive_read_data(a, buff, sizeof(buff));
    }

    // free objects created in this function
    if (!keep) {
      aurostd::RemoveFile(file);
    }
    r = archive_write_free(a_new);
    CompressionErrorHandling(a, r, error_message);
    r = archive_write_free(a);
    CompressionErrorHandling(a, r, error_message);
  }

  /// @brief convert compressed files in a directory to another compression
  /// @param directory to search for compressed files
  /// @param ct compression type
  /// @param keep the original compressed files
  /// @authors
  /// @mod{SD,20240329,created function}
  void compressfiles2compressfiles(const string& directory, const compression_type ct, const bool keep) {
    const string directory_clean = aurostd::CleanFileName(directory);
    if (directory_clean.empty()) {
      return;
    }
    const string& ext_new = compression_suffix[ct];
    string ext = ext_new;
    for (const std::filesystem::directory_entry& dir_entry : std::filesystem::directory_iterator(directory_clean)) {
      ext = dir_entry.path().extension().string();
      if ((ext == ".bz2" || ext == ".gz" || ext == ".xz" || ext == ".zip" || ext == ".zstd") && (ext != ext_new)) {
        aurostd::compressfile2compressfile(dir_entry.path().string(), ct, keep);
      }
    }
  }

} // namespace aurostd
