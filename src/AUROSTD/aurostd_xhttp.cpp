// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// First written by Frisco Rose in 2018
// Complete rewrite by Hagen Eckert in 2022
// Changed to use CURL by Hagen Eckert in 2025
// hagen.eckert@duke.edu

#include "aurostd_xhttp.h"

#include "config.h"

#include <array>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <curl/curl.h>
#include <curl/easy.h>
#include <curl/header.h>

#include "aurostd.h"
#include "aurostd_automatic_template.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"

#include "aflow_xhost.h" // todo required for XPID use and XHOST.DEBUG use

using std::cerr;
using std::endl;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

#define _DEBUG_XHTTP_ false

namespace {
  // anonymous namespace for functions that are only used locally, providing basic implementation

  /// @brief callback function to save data collected by curl into a std::string
  /// @param ptr delivered data
  /// @param size always 1
  /// @param nmemb size of the data
  /// @param data pointer to string that is used to save the information
  /// @return size written
  /// @xlink{"CURL_WRITEFUNCTION definition",https://curl.se/libcurl/c/CURLOPT_WRITEFUNCTION.html}
  /// @authors
  /// @mod{HE,20250602,created}
  size_t writeData2string(void* ptr, size_t size, size_t nmemb, void* data) {
    static_cast<std::string*>(data)->append(static_cast<char*>(ptr), size * nmemb);
    return size * nmemb;
  }

  /// @brief callback function to save data collected by curl into a file
  /// @param ptr delivered data
  /// @param size always 1
  /// @param nmemb size of the data
  /// @param stream pointer a file descriptor FILE*
  /// @return size written
  /// @xlink{"CURL_WRITEFUNCTION definition",https://curl.se/libcurl/c/CURLOPT_WRITEFUNCTION.html}
  /// @authors
  /// @mod{HE,20250602,created}
  size_t writeData2File(void* ptr, size_t size, size_t nmemb, void* stream) {
    const size_t written = fwrite(ptr, size, nmemb, static_cast<FILE*>(stream));
    if (written != nmemb * size) {
      std::stringstream message;
      message << "error " << errno << ": " << strerror(errno);
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return written;
  }

  /// @brief function to check for curl errors and throw if needed
  /// @param status curl status
  /// @param errbuf curl error buffer
  /// @xlink{"list of CURL errors",https://curl.se/libcurl/c/libcurl-errors.html}
  /// @authors
  /// @mod{HE,20250602,created}
  void checkCURLStatus(const CURLcode& status, const std::array<char, CURL_ERROR_SIZE>& errbuf) {
    if (status == CURLE_OK) {
      return;
    }
    const string message(errbuf.data());
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }

  /// @brief wrapper around CURL lib calls to get the content of a webpage
  /// @param url web address
  /// @param filename save to this location if not empty
  /// @param output save to this string if filename in empty
  /// @param response_code response code of the server
  /// @param header map of returned headers
  /// @xlink{"CURL lib",https://curl.se/libcurl/c/}
  /// @note this function is not designed to be used directly
  /// @authors
  /// @mod{HE,20250602,created}
  void httpGetBase(const string& url, const std::string& filename, std::string& output, long& response_code, std::map<std::string, std::string>& header) {
    const bool LDEBUG = (false || XHOST.DEBUG);

    // ensure the initial response code is set
    response_code = -1;
    // define a user agent string
    const std::string agent = "AFLOW/" + string(AFLOW_VERSION);
    // create the basic curl handle - cleaned up with curl_easy_cleanup at the end of the function
    CURL* handle = curl_easy_init();
    // define two header structs to loop through all headers
    struct curl_header* raw_header = nullptr;
    struct curl_header* prev_header = nullptr;
    // reserve space to capture human friendly error information and make sure it is empty
    std::array<char, CURL_ERROR_SIZE> errbuf;
    errbuf[0] = 0;
    // initialize the curl status to be OK before the first use
    CURLcode status = CURLE_OK;

    // set up the CURL library
    // set error buffer
    curl_easy_setopt(handle, CURLOPT_ERRORBUFFER, errbuf.data());
    // set url (the URL is not parsed until curl_easy_perform)
    status = curl_easy_setopt(handle, CURLOPT_URL, url.c_str());
    checkCURLStatus(status, errbuf);
    // allways follow redirections (30X)
    status = curl_easy_setopt(handle, CURLOPT_FOLLOWLOCATION, 1L); // starting with LIBCURL 8.13.0 1L could be replaced with CURLFOLLOW_ALL
    checkCURLStatus(status, errbuf);
    // set the user agent
    status = curl_easy_setopt(handle, CURLOPT_USERAGENT, agent.c_str());
    checkCURLStatus(status, errbuf);
    if constexpr (_DEBUG_XHTTP_) {
      // be very verbose during debugging
      status = curl_easy_setopt(handle, CURLOPT_VERBOSE, 1L);
      checkCURLStatus(status, errbuf);
    }
    if (LDEBUG) {
      // show progressbar during debugging
      status = curl_easy_setopt(handle, CURLOPT_NOPROGRESS, 0L);
      checkCURLStatus(status, errbuf);
    } else {
      // disable progressbar
      status = curl_easy_setopt(handle, CURLOPT_NOPROGRESS, 1L);
      checkCURLStatus(status, errbuf);
    }
    if (filename.empty()) {
      // use the write to string callback function
      status = curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, writeData2string);
      checkCURLStatus(status, errbuf);
      // set the output location to the output string
      status = curl_easy_setopt(handle, CURLOPT_WRITEDATA, &output);
      checkCURLStatus(status, errbuf);
      // perform the actual interaction with the webpage
      status = curl_easy_perform(handle);
      checkCURLStatus(status, errbuf);
    } else {
      // use the write to file callback function
      status = curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, writeData2File);
      checkCURLStatus(status, errbuf);
      // clean the filename and open it
      const std::string filename_clean = aurostd::CleanFileName(filename);
      if (FILE* output_file = fopen(filename_clean.c_str(), "wb")) {
        // set the file descriptor as the write target if it was created successful
        status = curl_easy_setopt(handle, CURLOPT_WRITEDATA, output_file);
        checkCURLStatus(status, errbuf);
        // perform the actual interaction with the webpage
        status = curl_easy_perform(handle);
        checkCURLStatus(status, errbuf);
        // close the file after writing
        if (const int r = fclose(output_file); r != 0) {
          const string message = "Unable to close file " + filename_clean;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
      } else {
        // throw an error if the file could not be opened
        const string message = "Unable to open file " + filename_clean;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
    // retrieve the response code
    status = curl_easy_getinfo(handle, CURLINFO_RESPONSE_CODE, &response_code);
    checkCURLStatus(status, errbuf);
    // collect all relevant headers
    // CURLH_HEADER The header arrived as a header from the server.
    // CURLH_1XX The header arrived in an HTTP 1xx response. A 1xx response is an "intermediate" response that might happen before the "real" response.
    // CURLH_TRAILER The header arrived as a trailer. A header that arrives after the body.
    constexpr unsigned int header_origin = CURLH_HEADER | CURLH_1XX | CURLH_TRAILER;
    while ((raw_header = curl_easy_nextheader(handle, header_origin, -1, prev_header)) != nullptr) {
      header.insert({raw_header->name, raw_header->value});
      prev_header = raw_header;
    }
    // clean up curl lib
    curl_easy_cleanup(handle);
  }
} // namespace

namespace aurostd {
  /// @brief get a web resource as string
  /// @param url web address
  /// @param output save to this string
  /// @param response_code response code of the server
  /// @param header map of returned headers
  void httpGet(const string& url, std::string& output, long& response_code, std::map<std::string, std::string>& header) {
    httpGetBase(url, "", output, response_code, header);
  }

  /// @brief get a web resource as file
  /// @param url web address
  /// @param filename save to this location
  /// @param response_code response code of the server
  /// @param header map of returned headers
  void httpGetFile(const string& url, const std::string& filename, long& response_code, std::map<std::string, std::string>& header) {
    std::string discarded_output;
    httpGetBase(url, filename, discarded_output, response_code, header);
  }

  /// @brief Retrieve data from an url
  /// @param url content url
  /// @return HTTP response code (-1 on failure)
  long httpGetStatus(const std::string& url) {
    std::string output;
    long response_code = -1;
    std::map<std::string, std::string> header;
    httpGet(url, output, response_code, header);
    return response_code;
  }

  /// @brief Retrieve data from an url
  /// @param url content url
  /// @param output message body
  /// @return HTTP response code (-1 on failure)
  long httpGetStatus(const std::string& url, std::string& output) {
    long response_code = -1;
    std::map<std::string, std::string> header;
    httpGet(url, output, response_code, header);
    return response_code;
  }

  /// @brief Retrieve data from an url
  /// @param url content url
  /// @param output message body
  /// @param header response header
  /// @return HTTP response code (-1 on failure)
  long httpGetStatus(const std::string& url, std::string& output, std::map<std::string, std::string>& header) {
    long response_code = -1;
    httpGet(url, output, response_code, header);
    return response_code;
  }

  /// @brief Retrieve data from url components
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @param output message body
  /// @return HTTP response code (-1 on failure)
  /// @note using http:// as this function was designed for AFLOW API calls, so encryption is not helpful
  long httpGetStatus(const std::string& host, const std::string& path, const std::string& query, std::string& output) {
    const std::string url = "http://" + host + "/" + path + query;
    long response_code = -1;
    std::map<std::string, std::string> header;
    httpGet(url, output, response_code, header);
    return response_code;
  }

  /// @brief Retrieve data from url components
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @param output message body
  /// @param header response header
  /// @return HTTP response code (-1 on failure)
  long httpGetStatus(const std::string& host, const std::string& path, const std::string& query, std::string& output, std::map<std::string, std::string>& header) {
    const std::string url = "http://" + host + "/" + path + query;
    long response_code = -1;
    httpGet(url, output, response_code, header);
    return response_code;
  }

  /// @brief Retrieve data from an url
  /// @param url content url
  /// @return message body
  std::string httpGet(const std::string& url) {
    std::string output;
    long response_code = -1;
    std::map<std::string, std::string> header;

    httpGet(url, output, response_code, header);
    return output;
  }

  /// @brief Retrieve data from an url
  /// @param url content url
  /// @param response_code HTTP response code (-1 on failure)
  /// @return message body
  std::string httpGet(const std::string& url, long& response_code) {
    std::string output;
    std::map<std::string, std::string> header;

    response_code = -1;
    httpGet(url, output, response_code, header);
    return output;
  }

  /// @brief Retrieve data from an url
  /// @param url content url
  /// @param response_code HTTP response code (-1 on failure)
  /// @param header response header
  /// @return message body
  std::string httpGet(const std::string& url, long& response_code, std::map<std::string, std::string>& header) {
    std::string output;
    response_code = -1;
    httpGet(url, output, response_code, header);
    return output;
  }

  /// @brief Download data as file
  /// @param url content url
  /// @param filename save to this location
  /// @return HTTP response code (-1 on failure)
  long httpGetFileStatus(const std::string& url, const std::string& filename) {
    long response_code = -1;
    std::map<std::string, std::string> header;
    httpGetFile(url, filename, response_code, header);
    return response_code;
  }

  /// @brief download data as file
  /// @param url content url
  /// @param filename save to this location
  /// @param header response header
  /// @return HTTP response code (-1 on failure)
  long httpGetFileStatus(const std::string& url, const std::string& filename, std::map<std::string, std::string>& header) {
    long response_code = -1;
    httpGetFile(url, filename, response_code, header);
    return response_code;
  }

  /// @brief get the data split into tokens
  /// @param url content url
  /// @param tokens data
  /// @param delimiters separators to use for splitting
  /// @return number of tokens
  template <typename utype> size_t httpGetTokens(const string& url, vector<utype>& tokens, const string& delimiters) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Loading url=" << url << endl;
    }
    tokens.clear();
    string content;
    const long status = httpGetStatus(url, content);
    if (status != 200 || content.empty()) {
      return 0;
    }

    vector<string> stokens;
    aurostd::string2tokens(content, stokens, delimiters);
    for (const string& stoken : stokens) {
      if (!stoken.empty()) {
        if constexpr (std::is_convertible_v<utype, std::string_view>) {
          tokens.push_back(stoken);
        } else {
          tokens.push_back(aurostd::string2utype<utype>(stoken));
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "Loaded " << tokens.size() << " tokens from " << url << endl;
    }
    return tokens.size();
  }

#define AST_TEMPLATE(utype) template size_t httpGetTokens(const string&, std::vector<utype>&, const string&);
  AST_GEN_1(AST_UTYPE_NUM)
  AST_GEN_1(AST_UTYPE_STRING)
#undef AST_TEMPLATE

  /// @brief get the content of a file from the web, decompress locally if needed
  /// @param url content url
  /// @param response_code HTTP response code (-1 on failure)
  /// @param header response header
  /// @return message body
  /// @note could be optimized by avoiding temporary files, but this function is used rarely
  std::string httpGetCompressedFileContent(const string& url, long& response_code, std::map<std::string, std::string>& header) {
    std::string content;
    const string ext = GetFileExtension(url);
    if (ext.empty()) {
      return httpGet(url, response_code, header);
    }
    const string temp_file = aurostd::TmpFileCreate("eurl2string") + ext;
    httpGetFile(url, temp_file, response_code, header);
    compressfile2string(temp_file, content);
    return content;
  }

  /// @brief get the content of a file from the web, decompress locally if needed
  /// @param url content url
  /// @param response_code HTTP response code (-1 on failure)
  /// @return message body
  /// @note could be optimized by avoiding temporary files, but this function is used rarely
  std::string httpGetCompressedFileContent(const string& url, long& response_code) {
    std::map<std::string, std::string> header;
    return httpGetCompressedFileContent(url, response_code, header);
  }

  /// @brief get the content of a file from the web, decompress locally if needed
  /// @param url content url
  /// @return message body
  /// @note could be optimized by avoiding temporary files, but this function is used rarely
  std::string httpGetCompressedFileContent(const string& url) {
    std::map<std::string, std::string> header;
    long response_code = -1;
    return httpGetCompressedFileContent(url, response_code, header);
  }

  /// @brief Fully percent encode a string
  /// @param work_str sting to escape
  /// @return escaped string
  /// @note https://www.rfc-editor.org/rfc/rfc3986#section-2.1
  /// @note just leave unreserved characters (ALPHA / DIGIT / "-" / "." / "_" / "~")
  string httpPercentEncodingFull(string work_str) {
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XHTTP_);

    const char* allowed =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789"
        "-_.~";

    size_t pos = 0;
    int to_replace = 0;
    std::stringstream output;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Escaping '" << work_str << "'" << endl;
    }

    while (!work_str.empty()) {
      pos = std::strspn(work_str.c_str(), allowed);
      to_replace = work_str[pos];
      if (to_replace < 0) {
        to_replace += 256;
      }
      output << work_str.substr(0, pos) << "%" << std::uppercase << std::hex << std::setfill('0') << std::setw(2) << to_replace;
      if (LDEBUG) {
        cerr << " Match '" << work_str[pos] << "' (%" << std::uppercase << std::hex << std::setfill('0') << to_replace << std::dec << ")" << endl;
      }
      work_str.erase(0, pos + 1);
    }
    return output.str();
  }
} // namespace aurostd
