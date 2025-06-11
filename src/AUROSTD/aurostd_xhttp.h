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

namespace aurostd{
  struct xURL {
    std::string scheme;
    std::string user;
    std::string host;
    unsigned int port;
    std::string path;
    std::string query;
  };

  int httpGetStatus(const std::string &url);
  int httpGetStatus(const std::string &url, std::string &output);
  int httpGetStatus(const std::string &url, std::string &output, std::map<std::string, std::string> &header);

  int httpGetStatus(const std::string &host, const std::string &path, const std::string &query, std::string &output);
  int httpGetStatus(const std::string &host, const std::string &path, const std::string &query, std::string &output, std::map<std::string, std::string> &header);

  std::string httpGet(const std::string &url);
  std::string httpGet(const std::string &url, int &status_code);
  std::string httpGet(const std::string &url, int &status_code, std::map<std::string, std::string> &header);

  std::string httpGet(const std::string &host, const std::string &path, const std::string &query);
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query, int &status_code);
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query, int &status_code, std::map<std::string, std::string> &header);

  std::string httpPercentEncodingFull(std::string work_str);

  xURL httpParseURL(const std::string &url, const bool strict = false);
}

#endif //AFLOW_SRC_AUROSTD_XHTTP_H
