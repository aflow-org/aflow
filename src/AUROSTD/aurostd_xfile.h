// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XFILE_H_
#define _AUROSTD_XFILE_H_

#include <cstddef>
#include <deque>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <archive.h>

namespace aurostd {
  namespace fs = std::filesystem;

  // compresion helper //HE20241216
  enum compression_type { None, XZ, ZIP, BZ2, GZ, ZSTD };
  const std::vector<std::string> compression_suffix = {"", ".xz", ".zip", ".bz2", ".gz", ".zst"};
  compression_type GetCompressionType(const fs::path& path);

  // make directories
  bool DirectoryMake(const std::string& directory_raw);
  bool SSH_DirectoryMake(const std::string& user, std::string machine, std::string Directory);

  // list directory content
  bool SubDirectoryLS(const std::string& directory_raw, std::vector<std::string>& vsubd);  // CO20200731
  bool DirectoryLS(const std::string& directory_raw, std::vector<std::string>& vfiles);
  bool DirectoryLS(const std::string& directory_raw, std::deque<std::string>& vfiles);

  // filename helper
  std::string dirname(const std::string& file);  // CO20210315
  std::string basename(const std::string& file); // CO20210315
  std::string GetFileExtension(const std::string& FileName);
  std::string CleanFileName(const std::string& fileIN);

  // aflow directory status
  bool DirectoryLocked(const std::string& directory, const std::string& LOCK);
  bool DirectoryLocked(const std::string& directory);
  bool DirectorySkipped(const std::string& directory);

  // general file/directory status
  bool DirectoryWritable(const std::string& directory_raw);
  bool IsDirectory(const std::string& path);
  bool IsCompressed(const std::string& FileNameIn, std::string& FileNameOut);
  bool IsCompressed(const std::string& FileNameIn);
  bool FileExist(const std::string& FileName);
  bool FileExist(const std::string& FileName, std::string& FileNameOut);
  bool CompressFileExist(const std::string& FileName);
  bool CompressFileExist(const std::string& FileNameRaw, std::string& FileNameOut);
  bool CompressFileWithinList(const std::vector<std::string>& list, const std::string& input);
  bool CompressFileWithinList(const std::vector<std::string>& list, const std::string& input, std::string& output);
  bool FileEmpty(const std::string& FileNameRaw);
  bool FileNotEmpty(const std::string& FileName);
  long int GetTimestampModified(const std::string&);  // ME20180712
  long int SecondsSinceFileModified(const std::string&);  // CO20210315
  unsigned long long int FileSize(const std::string& FileName);  // ME20191001
  void InFileExistCheck(const std::string& routine, const std::string& FileNameRaw, const std::ifstream& file_to_check);

  // temporary files
  std::string TmpStrCreate(const std::string& identifier_raw = "", const std::string& tmpdir_raw = "", bool hidden = false, bool directory = false);  // CO20210624
  std::string TmpFileCreate(const std::string& identifier = "", const std::string& tmpdir = "", bool hidden = false);  // CO20210315 - empty tmpdir means use XHOST.tmpfs
  std::string TmpDirectoryCreate(const std::string& identifier = "", const std::string& tmpdir = "", bool hidden = false);  // CO20210315 - empty tmpdir means use XHOST.tmpfs

  // copy link move
  bool CopyFile(const std::string& file_from, const std::string& file_to);
  bool LinkFile(const std::string& file_from, const std::string& file_to, bool soft = true);
  bool file2directory(const std::string& file, const std::string& destination);
  bool file2directory(const std::vector<std::string>& files, const std::string& destination);
  bool file2file(const std::string& file, const std::string& destination); // CO20171025

  // delete
  bool RemoveFile(const std::string& file);
  bool RemoveFile(const std::string& directory, const std::regex& regex_pattern);
  bool RemoveFiles(const std::string& directory);
  bool RemoveFile(const std::vector<std::string>& files);
  bool RemoveDirectory(const std::string& directory);

  // read file to std::string
  std::string file2string(const std::string& FileNameIN);
  size_t file2string(const std::string& FileNameIN, std::string& StringIN);
  std::string compressfile2string(const std::string& FileNameIN);
  size_t compressfile2string(const std::string& FileNameIN, std::string& content);

  // read file to vector of std::string
  size_t file2vectorstring(const std::string& FileNameIN, std::vector<std::string>& vline, bool consecutive = false, bool trim_edges = false);
  size_t file2vectorstring(const std::string& FileNameIN, std::deque<std::string>& vline);
  size_t compressfile2vectorstring(const std::string& FileNameIN, std::vector<std::string>& vline, bool consecutive = false, bool trim_edges = false);
  size_t compressfile2vectorstring(const std::string& FileNameIN, std::deque<std::string>& vline);

  // read file into a stringstream
  void file2stringstream(const std::string& FileNameIN, std::stringstream& StringstreamIN);
  void compressfile2stringstream(const std::string& file_raw, std::stringstream& content);

  // write files
  bool string2file(const std::string& StringOUTPUT, const std::string& FileNameRawOUTPUT, compression_type ct = None, const std::string& mode = "");
  bool string2compressfile(const std::string& content, const std::string& file_raw, compression_type ct = XZ);
  bool stringstream2file(const std::stringstream& StringstreamOUTPUT, const std::string& FileNameOUTPUT, compression_type ct = None, const std::string& mode = "");
  bool stringstream2compressfile(const std::stringstream& content, const std::string& file, compression_type ct = XZ);
  template <template <class...> class vtype> bool vectorstring2file(const vtype<std::string>& vline, const std::string& FileNameOUT, compression_type ct = None);

  // compression helpers

  void CompressInit(archive* a, compression_type ct, bool single);
  void DecompressFile(const std::string& file_raw, const std::string& outfile_raw = "", bool keep = false);
  void DecompressFile(const std::string& file_raw, bool keep);
  void DecompressFiles(const std::string& archive_file, const std::string& directory = "", bool keep = false);
  void DecompressFiles(const std::string& archive_file, bool keep);
  void CompressFile(const std::string& file_raw, compression_type ct = XZ, bool keep = false); // with default
  fs::path CompressFiles(const std::vector<std::string>& files, const fs::path& relative_to_path, const std::string& archive_name, compression_type ct = XZ, bool keep = false);
  bool MatchCompressed(const std::string& CompressedFileName, const std::string& FileNameOUT);
  bool compressfile2tempfile(const std::string& FileNameIN, std::string& FileNameOUT);
  bool compressfile2tempfile(const std::string& FileNameRawIN, std::string& FileNameOUT, bool& tempfile_created);
  void compressfile2compressfile(const std::string& file, compression_type ct, bool keep = false);
  void compressfiles2compressfiles(const std::string& directory, compression_type ct, bool keep = false);

  namespace JSON { // forward declaration //HE20241216
    struct object;
  }
  // Helper functions to retrieve embedded files //HE20230413
  namespace EmbData {
    JSON::object get_unit_test();
    std::string get_test_file(std::string file);
    std::string get_content(const std::string& filename, const std::string& collection);
    void save_to_file(const std::string& filename, const std::string& collection, const std::string& target_path);
    void save_to_folder(const std::string& path, const std::string& collection, const std::string& target_path);
  } // namespace EmbData
} // namespace aurostd
#endif  //_AUROSTD_XFILE_H_
