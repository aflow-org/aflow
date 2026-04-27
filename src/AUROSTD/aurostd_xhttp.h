// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// Written by Frisco Rose in 2018
// Complete rewrite by Hagen Eckert in 2022
// hagen.eckert@duke.edu

#ifndef AFLOW_SRC_AUROSTD_XHTTP_H
#define AFLOW_SRC_AUROSTD_XHTTP_H

#include <cstddef>
#include <deque>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace aurostd {
  void httpGet(const std::string& url, std::string& output, long& response_code, std::map<std::string, std::string>& header);

  std::string httpGet(const std::string& url);
  std::string httpGet(const std::string& url, long& response_code);
  std::string httpGet(const std::string& url, long& response_code, std::map<std::string, std::string>& header);

  long httpGetStatus(const std::string& url);
  long httpGetStatus(const std::string& url, std::string& output);
  long httpGetStatus(const std::string& url, std::string& output, std::map<std::string, std::string>& header);

  long httpGetStatus(const std::string& host, const std::string& path, const std::string& query, std::string& output);
  long httpGetStatus(const std::string& host, const std::string& path, const std::string& query, std::string& output, std::map<std::string, std::string>& header);

  void httpGetFile(const std::string& url, const std::string& filename, long& response_code, std::map<std::string, std::string>& header);

  long httpGetFileStatus(const std::string& url, const std::string& filename);
  long httpGetFileStatus(const std::string& url, const std::string& filename, std::map<std::string, std::string>& header);

  std::string httpGetCompressedFileContent(const std::string& url);
  std::string httpGetCompressedFileContent(const std::string& url, long& response_code);
  std::string httpGetCompressedFileContent(const std::string& url, long& response_code, std::map<std::string, std::string>& header);

  template <typename utype> size_t httpGetTokens(const std::string& url, std::vector<utype>& tokens, const std::string& delimiters = " ");

  std::string httpPercentEncodingFull(std::string work_str);
} // namespace aurostd

#endif // AFLOW_SRC_AUROSTD_XHTTP_H
