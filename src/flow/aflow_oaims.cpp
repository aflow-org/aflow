// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - Corey Oses

#include <cstddef>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_xhost.h"

using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

using aurostd::xvector;

//---------------------------------------------------------------------------------
// class xAIMSOUT
//---------------------------------------------------------------------------------
bool xAIMSOUT::GetProperties(const string& stringIN, bool QUIET) {
  stringstream sss;
  sss.str(stringIN);
  if (filename.empty()) {
    filename = "string";
  }
  return xAIMSOUT::GetProperties(sss, QUIET);
}

bool xAIMSOUT::GetPropertiesFile(const string& fileIN, bool QUIET) {
  stringstream sss;
  aurostd::compressfile2stringstream(fileIN, sss);
  if (filename.empty()) {
    filename = fileIN;
  }
  return xAIMSOUT::GetProperties(sss, QUIET);
}

bool xAIMSOUT::GetPropertiesFile(const string& fileIN, uint natoms_check, bool QUIET) {
  const bool flag = GetPropertiesFile(fileIN, QUIET);
  if (aurostd::abs(natoms_check - natoms) > 0.1) {
    stringstream message;
    message << "natoms_check(" << natoms_check << ")!= (int) natoms(" << natoms << ")";
    throw aurostd::xerror(__AFLOW_FILE__, "xAIMSOUT::GetPropertiesFile():", message, _FILE_CORRUPT_); // CO20200624
  }
  return flag;
}

bool xAIMSOUT::GetPropertiesUrlFile(const string& url, const string& file, bool VERBOSE) {
  const string tmpfile = aurostd::TmpFileCreate("xAIMSOUT_GetProperties"); // CO20200502 - threadID
  aurostd::httpGetFileStatus(url + "/" + file, tmpfile);
  const bool out = GetPropertiesFile(tmpfile);
  aurostd::RemoveFile(tmpfile);
  return out;
}

bool xAIMSOUT::GetProperties(const stringstream& stringstreamIN, bool QUIET) {
  const bool LVERBOSE = (false || XHOST.DEBUG || !QUIET);
  bool ERROR_flag = false;
  clear();
  stringstream sss;
  sss.str(stringstreamIN.str());
  content = stringstreamIN.str();
  vcontent.clear();
  const vector<string> vline;
  vector<string> tokens;
  aurostd::string2vectorstring(content, vcontent);
  string line;
  if (filename.empty()) {
    filename = "stringstream";
  }
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: ---------------------------------" << endl;
  }
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: BEGIN" << endl;
  }
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: filename=[" << filename << "]" << endl;
  }
  if (LVERBOSE) {
    cerr.precision(12);
  }
  // ----------------------------------------------------------------------
  // grab natoms first
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: ---------------------------------" << endl;
  }
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: LOAD NATOMS" << endl;
  }
  line = "";
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (aurostd::substring2bool(vcontent.at(iline), "Number of atoms")) { // VASP
      aurostd::string2tokens(vcontent.at(iline), tokens, ":");
      if (tokens.size() == 2) {
        natoms = aurostd::string2utype<double>(tokens[1]);
      }
    }
  }
  if (natoms == 0) {
    if (LVERBOSE) {
      cerr << "WARNING - xAIMSOUT::GetProperties:" << " no natoms tag found";
    }
    ERROR_flag = true;
  }
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: natoms=" << natoms << endl;
  }
  // ----------------------------------------------------------------------
  // now grab vforces
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: ---------------------------------" << endl;
  }
  if (LVERBOSE) {
    cerr << "xAIMSOUT::GetProperties: LOAD VFORCES" << endl;
  }
  int iforce;
  for (int iline = (int) vcontent.size() - 1; iline >= 0; iline--) {  // NEW - FROM THE BACK
    if (vforces.empty() && aurostd::substring2bool(vcontent.at(iline), "Total atomic forces")) { // VASP
      for (size_t i = 1; i < natoms + 1 && iline + i < vcontent.size() - 1; i++) {
        line = aurostd::RemoveCharacter(vcontent[iline + i], '|');  // remove pesky | in the beginning
        aurostd::string2tokens(line, tokens, " ");
        if (tokens.size() != 4) {
          if (LVERBOSE) {
            cerr << "WARNING - xAIMSOUT::GetProperties: ";
            cerr << "ill-written forces line, tokens.size()=" << tokens.size() << ";";
            cerr << "should be 4" << endl;
          }
          ERROR_flag = true;
        } else {
          iforce = aurostd::string2utype<int>(tokens[0]);
          if ((iforce != (int) i)) {
            if (LVERBOSE) {
              cerr << "WARNING - xAIMSOUT::GetProperties: ";
              cerr << "missing force " << i << "; ";
              cerr << "found force " << iforce << " instead" << endl;
            }
            ERROR_flag = true;
          } else {
            xvector<double> force(3);
            force[1] = aurostd::string2utype<double>(tokens[1]);
            force[2] = aurostd::string2utype<double>(tokens[2]);
            force[3] = aurostd::string2utype<double>(tokens[3]);
            vforces.push_back(force);
            if (LVERBOSE) {
              cerr << "xAIMSOUT::GetProperties: found new force[" << vforces.size() << "]=" << force << endl;
            }
          }
        }
      }
    }
  }
  if (vforces.empty() || vforces.size() != natoms) {
    if (LVERBOSE) {
      cerr << "WARNING - xAIMSOUT::GetProperties: not all forces found, vforces.size()=" << vforces.size() << endl;
    }
  }
  // DONE NOW RETURN
  if (ERROR_flag && LVERBOSE) {
    cerr << "WARNING - xAIMSOUT::GetProperties: ERROR_flag set in xAIMSOUT" << endl;
  }
  if (ERROR_flag) {
    return false;
  }
  return true;
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
