//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
// Class to handle AFLOW exceptions.
// 
// Since AFLOW should return 0 when it terminates normally, 0 cannot be used
// as an error code. See README_AFLOW_EXCEPTIONS.TXT for details and examples.
// If exceptions are added, please update the README as well!

#ifndef _AUROSTD_XERROR_CPP_
#define _AUROSTD_XERROR_CPP_

#ifndef _AUROSTD_XERROR_H_
#include "aurostd_xerror.h"
#endif

namespace aurostd {
  // List of error types and errors ////////////////////////////////////////////
#define _AFLOW_NUM_ERR_TYPES_ 7  // Number of error types

  string xerror_PID;  //SC20200508

  static const std::string error_types[_AFLOW_NUM_ERR_TYPES_] = 
  {"", "Input Error", "File Error", "Value Error",
    "Index Error", "Runtime Error", "Allocation Error"};

  //Errors associated with those types
  static const std::string errors[_AFLOW_NUM_ERR_TYPES_][9] =
  {{"Generic error", "Illegal error code", "", "", "", "", "", "", ""},
    {"unknown flag", "missing flag", "input ambiguous", "illegal parameter", "number of parameters", "", "", "", ""},
    {"file not found", "wrong format", "file corrupt", "", "", "", "", "", ""},
    {"illegal value", "out of range", "", "", "", "", "", "", ""},
    {"illegal value", "out of bounds", "mismatch", "", "", "", "", "", ""},
    {"not initialized", "SQL error", "busy", "external command not found", "external command failed", "HTTP error", "", "", ""}, //CO20200531
    {"could not allocate", "insufficient memory", "", "", "", "", "", "", ""}};

  //Constructors////////////////////////////////////////////////////////////////
  xerror::xerror(const std::string& filename, const std::string& function, const std::string& msg, int code) {
    buildException(filename, function, msg, code);
  }

  xerror::xerror(const std::string& filename, const std::string& function, const std::stringstream& msg, int code) {
    buildException(filename, function, msg.str(), code);
  }

  xerror::xerror(const std::string& filename, const std::stringstream& function, const std::string& msg, int code) {
    buildException(filename, function.str(), msg, code);
  }

  xerror::xerror(const std::string& filename, const std::stringstream& function, const std::stringstream& msg, int code) {
    buildException(filename, function.str(), msg.str(), code);
  }

  //buildExeption///////////////////////////////////////////////////////////////
  // builds the xerror object.
  void xerror::buildException(const std::string& filename, const std::string& fname, const std::string& msg, const int& code) {
    file_name = filename;
    function_name = fname;
    error_message = msg;
    error_code = code;
    error_type = error_code/10;
    error_number = error_code%10;
    if (!codeValid()) {
      error_code = 2;
      error_type = 0;
      error_number = 2;
    }
  }

  //codeValid///////////////////////////////////////////////////////////////////
  // Checks if a valid error code was provided
  bool xerror::codeValid() {
    if (error_code < 0) return false;
    if (error_type > _AFLOW_NUM_ERR_TYPES_) return false;
    if (error_type == 0 && error_number == 0) return false;  // "Error code" 0 is reserved
    if (error_number > 0) {  // an error value of 0 is the generic error, so it's always valid
      if (errors[error_type][error_number -1] == "") return false;
    }
    return true;
  }

  /******************* Functions to build the message string *******************/
  //buildMessageString/////////////////////////////////////////////////////////
  std::string xerror::buildMessageString() {
    std::stringstream msgstr;
    msgstr << xerror_PID; //SC20200508
    msgstr << "ERROR " << error_code << " in ";
    msgstr << whereFunction() << ": ";  // function name
    msgstr << error_string() << endl; //<< " - ";  // error type  //CO20181226
    if (error_code == 2) {
      //   msgstr << xerror_PID;  //SC20200508
      msgstr << "There was an error, but the supplied error code is invalid. Please contact the developers. ";
      msgstr << "Supplied error message: ";
    }
    msgstr << what();  // detailed error message
    //[CO20211125 - OBSOLETE]if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){  //CO20200624
    //[CO20211125 - OBSOLETE]  msgstr << " [dir=" << XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN") << "]";
    //[CO20211125 - OBSOLETE]}
    return msgstr.str();
  }

  std::string xerror::whereFunction() { //CO20191201
    return function_name;
  }

  std::string xerror::whereFileName() { //CO20191201
    return file_name;
  }

  std::string xerror::error_string() {
    std::stringstream errorstr;
    if (error_type > 0) {
      errorstr << error_types[error_type];
      errorstr << " (";
      if (error_number == 0) {  // reserve 0 for errors that cannot be specified otherwise
        errorstr << "generic";
      } else {
        errorstr << errors[error_type][error_number - 1];
      }
      errorstr << ")";
    } else {
      errorstr << errors[0][error_number - 1];
    }
    return errorstr.str();
  }

  std::string xerror::what() {
    return error_message;
  }

  int xerror::whatCode() {
    return error_code;
  }
} // namespace aurostd
#endif

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2018-2021               *
// *                                                                         *
//****************************************************************************
