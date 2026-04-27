// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#include "aurostd_hash.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <sstream>
#include <string>
#include <vector>

#include <openssl/evp.h>
#include <openssl/types.h>

#include "aurostd.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"

// ***************************************************************************
// CRC64
// ***************************************************************************

/* Redis uses the CRC64 variant with "Jones" coefficients and init value of 0.
 *
 * Specification of this CRC64 variant follows:
 * Name: crc-64-jones
 * Width: 64 bites
 * Poly: 0xad93d23594c935a9
 * Reflected In: True
 * Xor_In: 0xffffffffffffffff
 * Reflected_Out: True
 * Xor_Out: 0x0
 * Check("123456789"): 0xe9c6d914c4b8d9ca
 *
 * Copyright (c) 2012, Salvatore Sanfilippo <antirez at gmail dot com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * * Neither the name of Redis nor the names of its contributors may be used
 * to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE. */
// #include <stdio.h>
// #include <stdint.h>

// typedef unsigned long long int uint64_t;
#ifndef UINT64_C
#define UINT64_C(c) c##ULL
#endif

using std::ifstream;
using std::string;
using std::stringstream;
using std::vector;

namespace aurostd {

  uint64_t crc64(uint64_t crc, const unsigned char* s, uint64_t l) {
    for (uint64_t j = 0; j < l; j++) {
      const uint8_t byte = s[j];
      crc = crc64_tab[(uint8_t) crc ^ byte] ^ (crc >> 8);
    }
    return crc;
  }

  uint64_t crc64(uint64_t crc, const string& s) {
    return crc64(crc, (unsigned char*) s.c_str(), (uint64_t) s.length());
  }

  uint64_t crc64(const string& s) { // HE20220404 runtime partner function to aurostd::ctcrc64 starting at crc 0
    return crc64(0, (unsigned char*) s.c_str(), (uint64_t) s.length());
  }

  string crc2string(uint64_t crc) {
    stringstream oss;
    oss << std::hex << (unsigned long long) crc;
    const string sss = oss.str();
    return aurostd::PaddedPRE(sss, 16, "0");
  }

  /// @brief return the CRC64 checksum of a file
  /// @authors
  /// @mod{SC,20200326,created function}
  string file2auid(const string& file) {
    const vector<string> vtokens;
    if (aurostd::FileExist(file)) {
      uint64_t crc = 0;
      crc = aurostd::crc64(crc, aurostd::compressfile2string(file)); // DON'T TOUCH THIS
      return aurostd::crc2string(crc);
    }
    return "";
  }

  /// @brief Generates the checksum of a file
  /// @note Taken from old APL/apl_hroutines
  /// @todo reimplement using small chunks instead of reading the whole file at once
  /// authors
  /// @mod{ME,20190219,moved to aurostd}
  unsigned int getFileCheckSum(const string& filename, const string& algo) {
    ifstream infile(filename.c_str(), std::ios::in | std::ios::binary);
    if (!infile.is_open()) {
      const string message = "Cannot open file " + filename + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }

    // Get file length
    infile.seekg(0, std::ios::end);
    unsigned long length = infile.tellg();
    infile.seekg(0, std::ios::beg);

    // Setup read buffer (for whole file)
    if (length % 2 != 0) {
      length++;
    }
    char* buffer = new char[length];
    buffer[length - 1] = 0x00;

    // Read it in!
    infile.read(buffer, length);
    infile.close();

    // Get checksum
    unsigned int checksum;
    if (algo == "Fletcher32") {
      checksum = getFletcher32((unsigned short*) buffer, length >> 1);
    } else {
      checksum = 0;
    }
    delete[] buffer;

    // Return value
    return checksum;
  }

  /// @brief Generates the 32 bit checksum of a string based on Fletcher's algorithm.
  /// @note Taken from old APL/apl_hroutines
  /// @xlink{Wikipedia, http://en.wikipedia.org/wiki/Fletcher%27s_checksum}
  /// authors
  /// @mod{ME,20190219,moved to aurostd}
  unsigned int getFletcher32(unsigned short* data, size_t len) {
    unsigned int sum1 = 0xffff;
    unsigned int sum2 = 0xffff;

    while (len) {
      unsigned tlen = len > 360 ? 360 : len;
      len -= tlen;
      do {
        sum1 += *data++;
        sum2 += sum1;
      } while (--tlen);
      sum1 = (sum1 & 0xffff) + (sum1 >> 16);
      sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }

    // Second reduction step to reduce sums to 16 bits
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);

    return sum2 << 16 | sum1;
  }

  /// @brief convert a binary char array into a hex-encoded version
  /// @param value to convert
  /// @param data_length length to read
  /// @return hex-encoded string
  /// @warning as `value` is only a pointer the code can not know the length of it at this point; extra care needs to be applyed to the passed value to avoid an out of bound read
  std::string binary2hex(const unsigned char* value, const uint data_length) {
    // convert into hex representation
    char hex[data_length * 2 + 1];
    for (int i = 0; i < data_length; ++i) {
      snprintf(hex + i * 2, 3, "%02x", value[i]);
    }
    hex[data_length * 2] = '\0';
    return hex;
  }

  /// @brief create a hex-encoded hash (message digest) of a string
  /// @param input string to be hashed
  /// @param hash_type choosen hash varient (e.g. `MD5`, `SHA1`, `SHA256`)
  /// @return hash as hex-encoded string
  /// @see file2hash
  std::string string2hash(const std::string& input, const std::string& hash_type) {
    // variables to store the temp hash content
    unsigned int md_len;
    unsigned char md_value[EVP_MAX_MD_SIZE];

    // initialize hash context mamager
    EVP_MD_CTX* hash_context = EVP_MD_CTX_new();
    ;

    // try to load get hash function
    const EVP_MD* md = EVP_get_digestbyname(hash_type.c_str());
    if (md == nullptr) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, hash_type + " is not a valid message digest", _VALUE_ILLEGAL_);
    }

    // create a context manager and load hash variant
    if (!EVP_DigestInit_ex2(hash_context, md, nullptr)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Can't create hashing context", _RUNTIME_ERROR_);
    }

    // update the hash based on the string
    EVP_DigestUpdate(hash_context, input.c_str(), strlen(input.c_str()));

    // save result and clean up
    EVP_DigestFinal_ex(hash_context, md_value, &md_len);
    EVP_MD_CTX_free(hash_context);

    // convert into hex representation
    return binary2hex(md_value, md_len);
  }

  /// @brief create a hex-encoded hash (message digest) of a file
  /// @param file_path file to be hashed
  /// @param hash_type choosen hash varient (e.g. `MD5`, `SHA1`, `SHA256`)
  /// @return hash as hex-encoded string
  /// @note to list the local available digests use `openssl list -digest-algorithms`
  /// @note by using openSSL, hardware accaleration is utalized making common cryptographic hashes faster than simple ones
  /// @note to find the fastest hash variant use `openssl speed md5 sha1 sha256` or `openssl speed` for oll
  /// @note on Apple M1 MAX: MD5 0.51 Gb/s, SHA1 1.83 Gb/s, SHA256 1.86 Gb/s
  /// @mod{HE,20240907,function created}
  std::string file2hash(const fs::path& file_path, const std::string& hash_type) {
    // variables to store the temp hash content
    unsigned int md_len;
    unsigned char md_value[EVP_MAX_MD_SIZE];

    // initialize hash context mamager
    EVP_MD_CTX* hash_context = EVP_MD_CTX_new();
    ;

    // try to load get hash function
    const EVP_MD* md = EVP_get_digestbyname(hash_type.c_str());
    if (md == nullptr) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, hash_type + " is not a valid message digest", _VALUE_ILLEGAL_);
    }

    // open the file for hashing
    std::ifstream file(file_path, std::ifstream::binary);
    if (!file) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Can't open file for hashing: " + file_path.string(), _FILE_NOT_FOUND_);
    }

    // create a context manager and load hash variant
    if (!EVP_DigestInit_ex2(hash_context, md, nullptr)) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Can't create hashing context", _RUNTIME_ERROR_);
    }

    // loop over the file in small chuncks and update the hash
    while (file.good()) {
      char buf[1024 * 16];
      file.read(buf, sizeof(buf));
      EVP_DigestUpdate(hash_context, buf, file.gcount());
    }

    // save result and clean up
    EVP_DigestFinal_ex(hash_context, md_value, &md_len);
    EVP_MD_CTX_free(hash_context);

    // convert into hex representation
    return binary2hex(md_value, md_len);
  }

  /// @brief convert a crc64 hash (uint64_t) to a compact human readable string
  /// @param crc hash result (or other uint64_t)
  /// @param width number of characters
  /// @note if the given crc is too small the result is padded with the first entry in aurostd::human_alphanum_choices
  /// @note if width is zero the max length string is returned without padding
  /// @authors
  /// @mod{HE,20230221,created}
  string crc2human(const uint64_t crc, const uint width) {
    string result;
    uint64_t work = crc;
    while (work > 0) {
      result += human_alphanum_choices[work % human_alphanum_choices.size()];
      work = work / human_alphanum_choices.size();
      if (width > 0 && result.size() >= width) {
        break;
      }
    }
    if (result.size() < width) {
      result.insert(result.size(), width - result.size(), human_alphanum_choices[0]);
    }
    return result;
  }

  /// @brief uses crc64 to create a compact human readable hash string
  /// @param input string to be hashed
  /// @param width number of characters
  /// @note if width is zero the max length string is returned without padding
  string crc2human(const string& input, const uint width) {
    return crc2human(crc64(input), width);
  }

} // namespace aurostd

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
