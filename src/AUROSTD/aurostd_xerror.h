//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
#ifndef _AUROSTD_XERROR_H_
#define _AUROSTD_XERROR_H_

namespace aurostd {
  // BEGIN CONSTANTS
  // Definitions for the named error code constant (see README, section 2.1.).

#define _GENERIC_ERROR_              1
#define _ILLEGAL_CODE_               2
#define _INPUT_ERROR_               10
#define _INPUT_UNKNOWN_             11
#define _INPUT_MISSING_             12
#define _INPUT_AMBIGUOUS_           13
#define _INPUT_ILLEGAL_             14
#define _INPUT_NUMBER_              15
#define _FILE_ERROR_                20
#define _FILE_NOT_FOUND_            21
#define _FILE_WRONG_FORMAT_         22
#define _FILE_CORRUPT_              23
#define _VALUE_ERROR_               30
#define _VALUE_ILLEGAL_             31
#define _VALUE_RANGE_               32
#define _INDEX_ERROR_               40
#define _INDEX_ILLEGAL_             41
#define _INDEX_BOUNDS_              42
#define _INDEX_MISMATCH_            43
#define _RUNTIME_ERROR_             50
#define _RUNTIME_INIT_              51
#define _RUNTIME_SQL_               52
#define _RUNTIME_BUSY_              53
#define _RUNTIME_EXTERNAL_MISS_     54 //CO20200531
#define _RUNTIME_EXTERNAL_FAIL_     55 //CO20200531
#define _RUNTIME_HTTP_              56 //ME20220426
#define _ALLOC_ERROR_               60
#define _ALLOC_ALLOCATE_            61
#define _ALLOC_INSUFFICIENT_        62

  extern string xerror_PID; //SC20200508
  // END CONSTANTS

  class xerror {
    public:
      xerror(const std::string&, const std::string&, const std::string&, int = 1);
      xerror(const std::string&, const std::string&, const std::stringstream&, int = 1);
      xerror(const std::string&, const std::stringstream&, const std::stringstream&, int = 1);
      xerror(const std::string&, const std::stringstream&, const std::string&, int = 1);
      ~xerror() throw() {};
      std::string what();
      int whatCode();
      std::string whereFunction();  //CO20191201
      std::string whereFileName();  //CO20191201
      std::string buildMessageString();
    private:
      int error_code;
      int error_number;
      int error_type;
      std::string file_name;
      std::string function_name;
      std::string error_message;

      void buildException(const std::string&, const std::string&, const std::string&, const int&);
      bool codeValid();
      std::string error_string();
  };
} // namespace aurostd

#endif
//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
