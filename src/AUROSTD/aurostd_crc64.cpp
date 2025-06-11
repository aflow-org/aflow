// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AUROSTD_CRC64_CPP_
#define _AUROSTD_CRC64_CPP_
#include "aurostd.h"

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
//#include <stdio.h>
//#include <stdint.h>

//typedef unsigned long long int uint64_t;
#ifndef UINT64_C
#define UINT64_C(c) c ## ULL
#endif
//
//#include "aurostd_crc64.h"

namespace aurostd {

  uint64_t crc64(uint64_t crc, const unsigned char *s, uint64_t l) {
    for (uint64_t j = 0; j < l; j++) {
      uint8_t byte = s[j];
      crc = crc64_tab[(uint8_t)crc ^ byte] ^ (crc >> 8);
    }
    return crc;
  }

  uint64_t crc64(uint64_t crc, const string s) {
    return crc64(crc,(unsigned char*) s.c_str(),(uint64_t) s.length());
  }

  uint64_t crc64(const string s) { //HE20220404 runtime partner function to aurostd::ctcrc64 starting at crc 0
    return crc64(0,(unsigned char*) s.c_str(),(uint64_t) s.length());
  }
  
  string crc2string(uint64_t crc) {
    stringstream oss; oss << std::hex << (unsigned long long) crc; 
    string sss=oss.str();
    return aurostd::PaddedPRE(sss,16,"0");
  }

  /// @brief convert a crc64 hash (uint64_t) to a compact human readable string
  /// @param crc hash result (or other uint64_t)
  /// @param width number of characters
  /// @note if the given crc is too small the result is padded with the first entry in aurostd::human_alphanum_choices
  /// @note if width is zero the max length string is returned without padding
  /// @authors
  /// @mod{HE,20230221,created}
  string crc2human(const uint64_t crc, const uint width){
    string result = "";
    uint64_t work = crc;
    while (work > 0) {
      result += human_alphanum_choices[work % human_alphanum_choices.size()];
      work = work / human_alphanum_choices.size();
      if (width>0 && result.size()>=width) break;
    }
    if (result.size() < width) result.insert(result.size(), width - result.size(), human_alphanum_choices[0]);
    return result;
  }

  /// @brief uses crc64 to create a compact human readable hash string
  /// @param crc hash result (or other uint64_t)
  /// @param width number of characters
  /// @note if width is zero the max length string is returned without padding
  string crc2human(const string & input, const uint width) {
    return crc2human(crc64(input), width);
  }

  /* Test main */
  int crc64_main(void) {
    //  printf("e9c6d914c4b8d9ca == %016llx\n",(unsigned long long) aurostd::crc64(0,(unsigned char*)"123456789",9));
    cerr << "TYPE1: e9c6d914c4b8d9ca == " << std::hex << (unsigned long long) aurostd::crc64(0,(unsigned char*)"123456789",9) << endl;
    cerr << "TYPE2: e9c6d914c4b8d9ca == " << std::hex << (unsigned long long) aurostd::crc64(0,(string)"123456789") << endl;
    cerr << "TYPE3: e9c6d914c4b8d9ca == " << aurostd::crc2string(crc64(0,(string)"123456789")) << endl;
    return 0;
  }
  
}

#endif  // _AURO_IMPLEMENTATIONS_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

