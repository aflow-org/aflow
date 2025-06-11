// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AEL_GET_STRESS_CPP
#define _AFLOW_AEL_GET_STRESS_CPP

#include <cstddef>
#include <deque>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "flow/aflow_avasp.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_kvasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/AEL/aflow_ael_elasticity.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

// ###############################################################################
//                  AFLOW Automatic Elasticity Library (AEL) (2014-2021)
// ###############################################################################
//
// Uses strain-stress calculations to obtain elastic constants of materials
// Based on original Python program written by M. de Jong et al.
// See Scientific Data 2, 150009 (2015) for details of original program
// See Phys. Rev. Materials 1, 015401 (2017) for details of this implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AEL
//

// *****************************************************************************************************************
// The following functions are for setting up AEL inputs for postprocessing runs called from other parts of AFLOW
// *****************************************************************************************************************
namespace AEL_functions {
  uint AEL_xvasp_flags_populate(_xvasp& xvasp, string& AflowIn, string& AflowInName, string& FileLockName, const string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) {
    ifstream FileAFLOWIN;
    string FileNameAFLOWIN;
    ostringstream aus;
    const vector<string> vAflowInCheck;
    bool ael_aflowin_found = false;
    bool Krun = true;
    const bool load_POSCAR_from_xvasp = false;
    aurostd::xoption USER_AEL_POISSON_RATIO;
    USER_AEL_POISSON_RATIO.option = false;
    // Set aflags
    aflags.Directory = directory_LIB;
    aurostd::StringstreamClean(aus);
    aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
    aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
    aus << _AELSTR_MESSAGE_ << "AflowInName = " << AflowInName << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    if (aflags.Directory[0] != '/' && aflags.Directory[0] != '.' && aflags.Directory[0] != ' ') {
      aflags.Directory = "./" + aflags.Directory;
    }
    aflags.KBIN_RUN_AFLOWIN = true;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN = false;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP = false;
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN = false;
    aflags.KBIN_DELETE_AFLOWIN = false;
    // CO20200502 START - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future
    // CT20200720 Moved finding aflow.in file to separate function

    // CT20200721 Calls function AEL_Get_AflowInName to find correct aflow.in filename
    AEL_functions::AEL_Get_AflowInName(AflowInName, directory_LIB, ael_aflowin_found);
    FileNameAFLOWIN = directory_LIB + "/" + AflowInName;

    // CO20200502 STOP - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future

    if (ael_aflowin_found) {
      // Set FileMESSAGE name
      // CT20200624 Moved down so LOCK file is only renamed if this is actually an AEL main directory and not an ARUN.AEL directory
      if (!FileLockName.empty()) {
        if (aurostd::FileExist(directory_LIB + "/" + FileLockName)) {
          aurostd::file2file(aurostd::CleanFileName(directory_LIB + "/" + FileLockName), aurostd::CleanFileName(directory_LIB + "/" + FileLockName + ".run"));
        }
        const string FileNameMessage = directory_LIB + "/" + FileLockName;
        FileMESSAGE.open(FileNameMessage.c_str(), std::ios::app);
      } else {
        if (aurostd::FileExist(directory_LIB + "/ael.LOCK")) {
          aurostd::file2file(aurostd::CleanFileName(directory_LIB + "/ael.LOCK"), aurostd::CleanFileName(directory_LIB + "/ael.LOCK.run"));
          const string FileNameMessage = directory_LIB + "/ael.LOCK";
          FileMESSAGE.open(FileNameMessage.c_str(), std::ios::app);
        } else if (aurostd::FileExist(directory_LIB + "/agl.LOCK")) {
          aurostd::file2file(aurostd::CleanFileName(directory_LIB + "/agl.LOCK"), aurostd::CleanFileName(directory_LIB + "/ael.LOCK.run"));
          const string FileNameMessage = directory_LIB + "/agl.LOCK";
          FileMESSAGE.open(FileNameMessage.c_str(), std::ios::app);
        } else {
          const string FileNameMessage = directory_LIB + "/ael.LOCK";
          FileMESSAGE.open(FileNameMessage.c_str(), std::ios::app);
        }
      }
      // Search for AEL aflow.in filename
      FileAFLOWIN.open(FileNameAFLOWIN.c_str(), std::ios::in);
      FileAFLOWIN.clear();
      FileAFLOWIN.seekg(0);
      AflowIn = "";
      char c;
      // READ _AFLOWIN_ and put into AflowInCheck
      while (FileAFLOWIN.get(c)) {
        AflowIn += c;
      }
      FileAFLOWIN.clear();
      FileAFLOWIN.seekg(0);
      AflowIn = aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
      vector<string> vAflowIn;
      aurostd::string2vectorstring(AflowIn, vAflowIn);
      // Set kflags
      kflags.KBIN_MPI = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI]");
      kflags.AFLOW_MODE_VASP = aurostd::substring2bool(AflowIn, "[AFLOW_MODE=VASP]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_VASP]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE]VASP"); // check VASP
      if (kflags.AFLOW_MODE_VASP && !aflags.KBIN_GEN_VASP_FROM_AFLOWIN) {
        aflags.KBIN_GEN_VASP_FROM_AFLOWIN = true;
      } // do vasp last, default
      kflags.KBIN_SYMMETRY_CALCULATION = aurostd::substring2bool(AflowIn, "[AFLOW_SYMMETRY]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_SYMMETRY]CALC", true);
      kflags.KBIN_SYMMETRY_NO_SCAN = aurostd::substring2bool(AflowIn, "[AFLOW_SYMMETRY]NO_SCAN", true);
      if (aurostd::substring2bool(AflowIn, "[AFLOW_SYMMETRY]SYM_EPS=", true)) {
        kflags.KBIN_SYMMETRY_EPS = aurostd::substring2utype<double>(AflowIn, "[AFLOW_SYMMETRY]SYM_EPS=", true);
      }
      // parameters for zip/compression
      kflags.KZIP_COMPRESS = true;
      aurostd::StringstreamClean(aus);
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP=none]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP=NONE]") || !aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP")) {
        kflags.KZIP_COMPRESS = false;
        for (int i = 0; i < 1; i++) {
          aus << "WWWWW  Warning no compression of output files... " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintWarningStream(FileMESSAGE, aus, XHOST.QUIET);
        }
      } else {
        if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP")) { // "[AFLOW_MODE_ZIP=" not found
          kflags.KZIP_BIN = DEFAULT_KZIP_BIN; // take default
          aus << "00000  MESSAGE Taking DEFAULT KZIP_BIN=\"" << kflags.KZIP_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP]")) { // "[AFLOW_MODE_ZIP]" not found
          kflags.KZIP_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_ZIP]");
          aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP=")) { // "[AFLOW_MODE_ZIP=" found
          kflags.KZIP_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_ZIP="), ']');
          aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
      }
      // parameters for AAPL
      aurostd::xoption KBIN_PHONONS_CALCULATION_AAPL;
      KBIN_PHONONS_CALCULATION_AAPL.option = false;
      KBIN_PHONONS_CALCULATION_AAPL.options2entry(AflowIn, string("[AFLOW_AAPL]KAPPA=|[AFLOW_PHONONS]KAPPA="), KBIN_PHONONS_CALCULATION_AAPL.option, KBIN_PHONONS_CALCULATION_AAPL.xscheme);
      KBIN_PHONONS_CALCULATION_AAPL.option |= aurostd::substring2bool(AflowIn, "[AFLOW_AAPL]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_AAPL]CALC", true); // legacy
      kflags.KBIN_PHONONS_CALCULATION_AAPL = KBIN_PHONONS_CALCULATION_AAPL.option;
      // Parameters for QHA-APL
      if (!kflags.KBIN_PHONONS_CALCULATION_AAPL) {
        kflags.KBIN_PHONONS_CALCULATION_QHA = aurostd::substring2bool(AflowIn, "[AFLOW_QHA]CALC", true) || aurostd::substring2bool(AflowIn, "VASP_QHA]CALC", true);
      }
      // parameters for APL
      if (!(kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA)) { // mutually exclusive
        kflags.KBIN_PHONONS_CALCULATION_APL = aurostd::substring2bool(AflowIn, "[AFLOW_APL]CALC", true) || aurostd::substring2bool(AflowIn, "[AFLOW_PHONONS]CALC", true) ||
                                              aurostd::substring2bool(AflowIn, "[VASP_PHONONS]CALC", true);
      }
      // parameters for AGL (Debye Model)
      for (size_t i = 0; i < vAflowIn.size() && !kflags.KBIN_PHONONS_CALCULATION_AGL; i++) {
        if ((aurostd::substring2bool(vAflowIn[i], "[AFLOW_AGL]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_AGL]CALC", true)) &&
            !(aurostd::substring2bool(vAflowIn[i], "[AFLOW_AGL]CALC_", true) || aurostd::substring2bool(vAflowIn[i], "[VASP_AGL]CALC_", true) || aurostd::substring2bool(vAflowIn[i], "[AFLOW_AGL]CALCS", true) ||
              aurostd::substring2bool(vAflowIn[i], "[VASP_AGL]CALCS", true) || false)) {
          kflags.KBIN_PHONONS_CALCULATION_AGL = true;
        }
      }
      // parameters for AEL (Elastic constants)
      for (size_t i = 0; i < vAflowIn.size() && !kflags.KBIN_PHONONS_CALCULATION_AEL; i++) {
        if ((aurostd::substring2bool(vAflowIn[i], "[AFLOW_AEL]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_AEL]CALC", true)) &&
            !(aurostd::substring2bool(vAflowIn[i], "[AFLOW_AEL]CALC_", true) || aurostd::substring2bool(vAflowIn[i], "[VASP_AEL]CALC_", true) || aurostd::substring2bool(vAflowIn[i], "[AFLOW_AEL]CALCS", true) ||
              aurostd::substring2bool(vAflowIn[i], "[VASP_AEL]CALCS", true) || false)) {
          kflags.KBIN_PHONONS_CALCULATION_AEL = true;
        }
      }
      // parameters for POCC CALCULATIONS
      kflags.KBIN_POCC = false;
      kflags.KBIN_POCC_CALCULATION = aurostd::substring2bool(AflowIn, "[AFLOW_POCC]CALC", true) &&
                                     (aurostd::substring2bool(AflowIn, "[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE", true) && aurostd::substring2bool(AflowIn, "[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE", true)); // CO20180419
      if (kflags.KBIN_POCC_CALCULATION) {
        kflags.KBIN_POCC = true;
      }
      // parameters for FROZSL
      kflags.KBIN_FROZSL = false;
      kflags.KBIN_PHONONS_CALCULATION_FROZSL = aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]CALC", true);
      kflags.KBIN_FROZSL_DOWNLOAD = (aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]DOWN", true) || aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]DOWNLOAD", true));
      kflags.KBIN_FROZSL_FILE = aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]FILE", true);
      if (kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_FROZSL_DOWNLOAD || kflags.KBIN_FROZSL_FILE) {
        kflags.KBIN_FROZSL = true;
      }
      // Set KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_VASP_BIN;
      KBIN::MPI_Extract(AflowIn, FileMESSAGE, aflags, kflags);
      // Set vflags from AflowIN
      vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, FileMESSAGE, aflags, kflags);

      // Set-up xvasp
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      xvasp.clear();
      const uint ixvasp = 0;
      xvasp.POSCAR_index = ixvasp;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      KBIN::readModulesFromAflowIn(AflowIn, kflags, xvasp);
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      xvasp.Directory = aflags.Directory;
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      if (Krun) {
        Krun = (Krun && KBIN::VASP_Produce_INPUT(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags, load_POSCAR_from_xvasp));
      }
      if (Krun) {
        Krun = (Krun && KBIN::VASP_Modify_INPUT(xvasp, FileMESSAGE, aflags, kflags, vflags));
      }
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ << "xvasp.Directory = " << xvasp.Directory << endl;
      aus << _AELSTR_MESSAGE_ << "aflags.Directory = " << aflags.Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      // Fix blank species
      if (!xvasp.str.species.empty()) {
        if (xvasp.str.species[0].empty()) {
          pflow::fixEmptyAtomNames(xvasp.str);
        }
      }
      if (Krun) {
        return 0;
      } else {
        return 1;
      }
    } else {
      cerr << _AELSTR_MESSAGE_ << "AEL input file not found!" << endl; // CT20200624 If no AEL file present, then this is not an AEL main directory: write to cerr and return 2 to inform calling function
      cerr << _AELSTR_MESSAGE_ << "Not an AEL main directory" << endl;
      return 2;
    }
  }
} // namespace AEL_functions

// *****************************************************************************************************************
// Finds aflow.in file for AEL calculation (if it exists) //CT20200715
// *****************************************************************************************************************
namespace AEL_functions {
  uint AEL_Get_AflowInName(string& AflowInName, const string& directory_LIB, bool& ael_aflowin_found) {
    ifstream FileAFLOWINcheck;
    string FileNameAFLOWIN;
    string FileNameAFLOWINcheck;
    string AflowInCheck;
    string aflowinname;
    string stmp;
    uint aelerror = 0;
    uint filelength = 0;
    vector<string> vAflowInCheck;
    ael_aflowin_found = false;
    aurostd::xoption USER_AEL_POISSON_RATIO;
    USER_AEL_POISSON_RATIO.option = false;
    vector<string> vaflowins;
    if (!AflowInName.empty()) {
      vaflowins.push_back(AflowInName);
    } // Check if AflowInName exists
    if (!_AFLOWIN_.empty()) {
      vaflowins.push_back(_AFLOWIN_);
    } // Otherwise, check if _AFLOWIN_ file is AEL input file
    vaflowins.push_back(_AFLOWIN_AEL_DEFAULT_); // Otherwise, check for other commonly used names for AEL aflow.in file
    vaflowins.push_back(_AFLOWIN_AGL_DEFAULT_);
    for (size_t iaf = 0; iaf < vaflowins.size() && !ael_aflowin_found; iaf++) {
      aflowinname = vaflowins[iaf];
      if (aurostd::CompressFileExist(directory_LIB + "/" + aflowinname, stmp) && aurostd::IsCompressed(stmp)) {
        aurostd::DecompressFile(stmp);
      } // CO20210204 - fix aflow.in.xz
      if ((!ael_aflowin_found) && (aurostd::FileExist(directory_LIB + "/" + aflowinname))) {
        FileNameAFLOWINcheck = directory_LIB + "/" + aflowinname;
        AflowInCheck = "";
        // READ aflowinname and put into AflowInCheck
        filelength = aurostd::file2string(FileNameAFLOWINcheck, AflowInCheck);
        if (filelength > 0) {
          aelerror = 0;
        } else {
          aelerror = 1;
        }
        AflowInCheck = aurostd::RemoveComments(AflowInCheck); // NOW Clean AFLOWIN
        vAflowInCheck.clear();
        aurostd::string2vectorstring(AflowInCheck, vAflowInCheck);
        // Check if aflowinname contains command to run AEL
        for (size_t i = 0; i < vAflowInCheck.size() && !ael_aflowin_found; i++) {
          if ((aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AEL]CALC", true) || aurostd::substring2bool(AflowInCheck, "[VASP_AEL]CALC", true)) &&
              !(aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AEL]CALC_", true) || aurostd::substring2bool(vAflowInCheck[i], "[VASP_AEL]CALC_", true) ||
                aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AEL]CALCS", true) || aurostd::substring2bool(vAflowInCheck[i], "[VASP_AEL]CALCS", true) || false)) {
            FileNameAFLOWIN = FileNameAFLOWINcheck;
            ael_aflowin_found = true;
            AflowInName = aflowinname;
          } else if ((aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AGL]CALC", true) || aurostd::substring2bool(AflowInCheck, "[VASP_AGL]CALC", true)) &&
                     !(aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AGL]CALC_", true) || aurostd::substring2bool(vAflowInCheck[i], "[VASP_AGL]CALC_", true) ||
                       aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AGL]CALCS", true) || aurostd::substring2bool(vAflowInCheck[i], "[VASP_AGL]CALCS", true) || false)) {
            if (aurostd::substring2bool(AflowInCheck, "[AFLOW_AGL]AEL_POISSON_RATIO=", true)) {
              USER_AEL_POISSON_RATIO.options2entry(AflowInCheck, "[AFLOW_AGL]AEL_POISSON_RATIO=", USER_AEL_POISSON_RATIO.option);
            } else if (aurostd::substring2bool(AflowInCheck, "[AFLOW_AGL]AELPOISSONRATIO=", true)) {
              USER_AEL_POISSON_RATIO.options2entry(AflowInCheck, "[AFLOW_AGL]AELPOISSONRATIO=", USER_AEL_POISSON_RATIO.option);
            }
            if (USER_AEL_POISSON_RATIO.option) {
              FileNameAFLOWIN = FileNameAFLOWINcheck;
              ael_aflowin_found = true;
              AflowInName = aflowinname;
            }
          }
        }
        FileAFLOWINcheck.close();
      }
    }
    return aelerror;
  }
} // namespace AEL_functions

// *******************************************************************************
// The following functions are for generating _AFLOWIN_ files
// *******************************************************************************

// ***************************************************************************
// AEL_functions::aelvaspflags
// ***************************************************************************
namespace AEL_functions {
  //
  // Function to assign values for VASP input flags from aflow.in file to vaspRun _xvasp class
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint aelvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, string& dirrunname, _AEL_data& AEL_data, ofstream& FileMESSAGE) {
    ostringstream aus;
    vector<string> vfile;
    string vfilename;
    bool vfileexist = false;
    if (AEL_data.relax_static || AEL_data.static_only) {
      aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"), vfile, ",");
      for (size_t ij = 0; ij < vfile.size(); ij++) {
        if (aurostd::FileExist(dirrunname + "/" + vfile[ij])) {
          vfilename = vfile[ij];
          vfileexist = true;
        }
      }
    } else {
      aurostd::string2tokens(string("OUTCAR.relax2.bz2,OUTCAR.relax2.gz,OUTCAR.relax2.xz,OUTCAR.relax2"), vfile, ",");
      for (size_t ij = 0; ij < vfile.size(); ij++) {
        if (aurostd::FileExist(dirrunname + "/" + vfile[ij])) {
          vfilename = vfile[ij];
          vfileexist = true;
        }
      }
    }
    // SOME WARNINGS: check existence of LOCK and OUTCAR.relax2 files
    if (!(aurostd::FileExist(dirrunname + "/" + _AFLOWLOCK_) ||
          ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/ael.LOCK") || aurostd::FileExist(dirrunname + "/LOCK")))) &&
        (vfileexist)) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "found " << vfilename << " but no " << _AFLOWLOCK_ << " in " << dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 1;
    }

    if ((aurostd::FileExist(dirrunname + "/" + _AFLOWLOCK_) ||
         ((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/ael.LOCK") || aurostd::FileExist(dirrunname + "/LOCK")))) &&
        !(vfileexist)) {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_WARNING_ + "found " << _AFLOWLOCK_ << " but no OUTCAR in " << dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 1;
    }

    // Switch off autotune
    kbinFlags.KBIN_MPI_AUTOTUNE = true;

    std::map<std::string, bool*> RUN_map = {
        {               "GENERATE",                    &vaspRun.AVASP_flag_GENERATE},
        {                  "RELAX",                   &vaspRun.AVASP_flag_RUN_RELAX},
        {           "RELAX_STATIC",            &vaspRun.AVASP_flag_RUN_RELAX_STATIC},
        {     "RELAX_STATIC_BANDS",      &vaspRun.AVASP_flag_RUN_RELAX_STATIC_BANDS},
        {"RELAX_STATIC_DIELECTRIC", &vaspRun.AVASP_flag_RUN_RELAX_STATIC_DIELECTRIC},
        {                 "STATIC",                  &vaspRun.AVASP_flag_RUN_STATIC},
        {           "STATIC_BANDS",            &vaspRun.AVASP_flag_RUN_STATIC_BANDS},
        {      "STATIC_DIELECTRIC",       &vaspRun.AVASP_flag_RUN_STATIC_DIELECTRIC}
    };
    for (auto& [name, flag] : RUN_map) {
      *flag = false;
    }

    if (AEL_data.relax_static) {
      vaspRun.AVASP_flag_RUN_RELAX_STATIC = true;
      *RUN_map["RELAX_STATIC"] = true;
      vaspFlags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry = true;
      vaspFlags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.content_string = "IONS";
      vaspRun.aopts.flag("FLAG::VOLUME_PRESERVED", true);
    } else {
      vaspRun.AVASP_flag_RUN_STATIC = true;
      *RUN_map["STATIC"] = true;
    }
    for (const auto& [name, flag] : RUN_map) {
      vaspFlags.KBIN_VASP_RUN.flag(name, *flag);
    }

    // Set unit cell conversion to "PRESERVE"
    setPreserveUnitCell(vaspRun);
    vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry = true;
    vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.pop("STANDARD_PRIMITIVE");
    vaspFlags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.push("PRESERVE");

    if (AEL_data.precaccalgonorm) {
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.isentry = true;
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.content_string = "ACCURATE";
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.isentry = true;
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.content_string = "NORMAL";
    }

    // Switch off VASP symmetry - this can help when applied strains break the symmetry of the primitive cell
    if (!AEL_data.vasp_symmetry) {
      vaspFlags.KBIN_VASP_FORCE_OPTION_SYM.option = false;
    }

    // Change format of POSCAR
    if ((!kbinFlags.KBIN_MPI && (kbinFlags.KBIN_BIN.find("46") != string::npos)) || (kbinFlags.KBIN_MPI && (kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
      vaspRun.str.is_vasp5_poscar_format = false;
    }

    return 0;
  }
} // namespace AEL_functions

// ************************************************************************************************
// This set of functions extract stress tensor data from VASP runs
// ************************************************************************************************

// ***************************************************************************
// AEL_functions::extractstress
// ***************************************************************************
namespace AEL_functions {
  //
  // extractstress: Extract stress tensors from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractstress(vector<_xvasp>& vaspRuns, _AEL_data& AEL_data, vector<string>& dirrunname, ofstream& FileMESSAGE) {
    const bool LVERBOSE = (false || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile;
    vector<string> dfile;
    xOUTCAR outcar;
    xVASPRUNXML vasprunxml;
    string vfilename;
    string dfilename;
    string ffilename;
    bool skipdir = false;
    aurostd::xmatrix<double> stress_tensor(3, 3);
    double pressure_val = 0.0;
    double energy_cell_val = 0.0;
    AEL_data.normal_stress.clear();
    AEL_data.shear_stress.clear();
    AEL_data.normal_stress.resize(AEL_data.normal_strain.size());
    AEL_data.shear_stress.resize(AEL_data.shear_strain.size());
    AEL_data.normal_deformations_complete.resize(AEL_data.normal_strain.size());
    AEL_data.shear_deformations_complete.resize(AEL_data.shear_strain.size());
    AEL_data.energycalculated.clear();
    AEL_data.pressurecalculated.clear();
    AEL_data.stresscalculated.clear();
    AEL_data.structurecalculated.clear();
    if (AEL_data.relax_static || AEL_data.static_only) {
      if (AEL_data.vasprunxmlstress) {
        aurostd::string2tokens(string("vasprun.xml.static.bz2,vasprun.xml.static.gz,vasprun.xml.static.xz,vasprun.xml.static"), vfile, ",");
      } else {
        aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"), vfile, ",");
      }
    } else if (AEL_data.relax_only) {
      if (AEL_data.vasprunxmlstress) {
        aurostd::string2tokens(string("vasprun.xml.relax2.bz2,vasprun.xml.relax2.gz,vasprun.xml.relax2.xz,vasprun.xml.relax2"), vfile, ",");
      } else {
        aurostd::string2tokens(string("OUTCAR.relax2.bz2,OUTCAR.relax2.gz,OUTCAR.relax2.xz,OUTCAR.relax2"), vfile, ",");
      }
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_ERROR_ + "No run type selected" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 1;
    }

    // Loop over normal strains and read stress tensors
    uint idVaspRun = 0;
    for (size_t i = 1; i <= AEL_data.normal_strain.size(); i++) {
      for (size_t j = 0; j < AEL_data.normal_deformations.size(); j++) {
        skipdir = false;
        if (idVaspRun > vaspRuns.size()) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          return 1;
        }
        const double strainfactor = 1.0 + AEL_data.normal_deformations[j];
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns[idVaspRun].Directory << endl;
        aus << _AELSTR_MESSAGE_ + "Normal stress: i = " << i << ", j = " << j << ", strain factor = " << strainfactor << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        aurostd::string2tokens(vaspRuns[idVaspRun].Directory, dfile, "/");
        dfilename = dfile.at(dfile.size() - 1);
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        // Check if structure is on list of failed runs to be skipped
        // If so, then skip reading and continue to next structure
        for (size_t ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
          ffilename = AEL_data.failed_arun_list[ij];
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
            aus << _AELSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
          if (aurostd::substring2bool(dfilename, ffilename, true)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            skipdir = true;
          }
        }
        if (skipdir) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          idVaspRun++;
          skipdir = false;
          continue;
        }
        // if(skipdir) continue;
        // continue;

        // If tarred and compressed directory exists...
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.bz2")) {
          aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.bz2");
        } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.bz2")) {
          aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.bz2");
        } // Extract all...
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.gz")) {
          aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.gz");
        } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.gz")) {
          aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.gz");
        } // Extract all...
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.xz")) {
          aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.xz");
        } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.xz")) {
          aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.xz");
        } // Extract all...

        // If the LOCK file is missing, then it is probably a corrupted run
        // Do not accept it and wait for the new run
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "LOCK file path = " << dirrunname[idVaspRun] + "/" + _AFLOWLOCK_ << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        if (!aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + _AFLOWLOCK_) && !aurostd::FileExist(dirrunname[idVaspRun] + "/" + _AFLOWLOCK_) &&
            !(((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && ((aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/agl.LOCK")) || (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/LOCK")) ||
                                                                    (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/ael.LOCK")) || (aurostd::FileExist(dirrunname[idVaspRun] + "/agl.LOCK")) ||
                                                                    (aurostd::FileExist(dirrunname[idVaspRun] + "/LOCK")) || (aurostd::FileExist(dirrunname[idVaspRun] + "/ael.LOCK")))))) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (AEL_data.autoskipfailedaruns) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            idVaspRun++;
            continue;
          } else {
            throw AELStageBreak();
          }
        }
        if (AEL_data.vasprunxmlstress) {
          for (size_t ij = 0; ij < vfile.size() && (vasprunxml.content.empty()); ij++) {
            if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns[idVaspRun].Directory + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[ij]);
              vfilename = vfile[ij];
            } else if (aurostd::FileExist(dirrunname[idVaspRun] + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname[idVaspRun] + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(dirrunname[idVaspRun] + "/" + vfile[ij]);
              vfilename = vfile[ij];
            }
          }
          if (vasprunxml.content.empty()) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -vasprunxml.stress;
          vasprunxml.clear();
        } else {
          for (size_t ij = 0; ij < vfile.size() && (outcar.content.empty()); ij++) {
            if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns[idVaspRun].Directory + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              outcar.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[ij]);
              vfilename = vfile[ij];
            } else if (aurostd::FileExist(dirrunname[idVaspRun] + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname[idVaspRun] + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              outcar.GetPropertiesFile(dirrunname[idVaspRun] + "/" + vfile[ij]);
              vfilename = vfile[ij];
            }
          }
          if (outcar.content.empty()) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -outcar.stress;
          energy_cell_val = outcar.energy_cell;
          pressure_val = outcar.pressure_residual;
          outcar.clear();
        }
        AEL_data.normal_stress[i - 1].push_back(stress_tensor);
        AEL_data.normal_deformations_complete[i - 1].push_back(AEL_data.normal_deformations[j]);
        AEL_data.energycalculated.push_back(energy_cell_val);
        AEL_data.pressurecalculated.push_back(pressure_val);
        AEL_data.stresscalculated.push_back(stress_tensor);
        AEL_data.structurecalculated.push_back(idVaspRun);
        idVaspRun++;
        // Print out stress tensor
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Normal stress tensor = " << AEL_data.normal_stress[i - 1].at(AEL_data.normal_stress[i - 1].size() - 1) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
    }

    // Loop over shear strains and read stress tensors
    for (size_t i = 1; i <= AEL_data.shear_strain.size(); i++) {
      for (size_t j = 0; j < AEL_data.shear_deformations.size(); j++) {
        skipdir = false;
        if (idVaspRun > vaspRuns.size()) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          return 1;
        }
        const double strainfactor = 1.0 + AEL_data.shear_deformations[j];
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns[idVaspRun].Directory << endl;
        aus << _AELSTR_MESSAGE_ + "Shear stress: i = " << i << ", j = " << j << ", strain factor = " << strainfactor << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        aurostd::string2tokens(vaspRuns[idVaspRun].Directory, dfile, "/");
        dfilename = dfile.at(dfile.size() - 1);
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        // Check if structure is on list of failed runs to be skipped
        // If so, then skip reading and continue to next structure
        for (size_t ij = 0; ij < AEL_data.failed_arun_list.size(); ij++) {
          ffilename = AEL_data.failed_arun_list[ij];
          if (LVERBOSE) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
            aus << _AELSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
          if (aurostd::substring2bool(dfilename, ffilename, true)) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            skipdir = true;
          }
        }
        if (skipdir) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          skipdir = false;
          idVaspRun++;
          continue;
        }

        // If tarred and compressed directory exists...
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.bz2")) {
          aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.bz2");
        } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.bz2")) {
          aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.bz2");
        } // Extract all...
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.gz")) {
          aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.gz");
        } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.gz")) {
          aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.gz");
        } // Extract all...
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.xz")) {
          aurostd::execute(string("tar -xf") + vaspRuns[idVaspRun].Directory + ".tar.xz");
        } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.xz")) {
          aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.xz");
        } // Extract all...

        // If the LOCK file is missing, then it is probably a corrupted run
        // Do not accept it and wait for the new run
        if (!aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + _AFLOWLOCK_) && !aurostd::FileExist(dirrunname[idVaspRun] + "/" + _AFLOWLOCK_) &&
            !(((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && ((aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/agl.LOCK")) || (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/LOCK")) ||
                                                                    (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/ael.LOCK")) || (aurostd::FileExist(dirrunname[idVaspRun] + "/agl.LOCK")) ||
                                                                    (aurostd::FileExist(dirrunname[idVaspRun] + "/LOCK")) || (aurostd::FileExist(dirrunname[idVaspRun] + "/ael.LOCK")))))) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (AEL_data.autoskipfailedaruns) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            idVaspRun++;
            continue;
          } else {
            throw AELStageBreak();
          }
        }
        if (AEL_data.vasprunxmlstress) {
          for (size_t ij = 0; ij < vfile.size() && (vasprunxml.content.empty()); ij++) {
            if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns[idVaspRun].Directory + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[ij]);
              vfilename = vfile[ij];
            } else if (aurostd::FileExist(dirrunname[idVaspRun] + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname[idVaspRun] + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              vasprunxml.GetPropertiesFile(dirrunname[idVaspRun] + "/" + vfile[ij]);
              vfilename = vfile[ij];
            }
          }
          if (vasprunxml.content.empty()) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -vasprunxml.stress;
          vasprunxml.clear();
        } else {
          for (size_t ij = 0; ij < vfile.size() && (outcar.content.empty()); ij++) {
            if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns[idVaspRun].Directory + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              outcar.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[ij]);
              vfilename = vfile[ij];
            } else if (aurostd::FileExist(dirrunname[idVaspRun] + "/" + vfile[ij])) {
              if (LVERBOSE) {
                aurostd::StringstreamClean(aus);
                aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
                aus << _AELSTR_MESSAGE_ + "file = " << dirrunname[idVaspRun] + "/" + vfile[ij] << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              outcar.GetPropertiesFile(dirrunname[idVaspRun] + "/" + vfile[ij]);
              vfilename = vfile[ij];
            }
          }
          if (outcar.content.empty()) {
            aurostd::StringstreamClean(aus);
            aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (AEL_data.autoskipfailedaruns) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              idVaspRun++;
              continue;
            } else {
              throw AELStageBreak();
            }
          }
          stress_tensor = -outcar.stress;
          energy_cell_val = outcar.energy_cell;
          pressure_val = outcar.pressure_residual;
          outcar.clear();
        }
        AEL_data.shear_stress[i - 1].push_back(stress_tensor);
        AEL_data.shear_deformations_complete[i - 1].push_back(AEL_data.shear_deformations[j]);
        AEL_data.energycalculated.push_back(energy_cell_val);
        AEL_data.pressurecalculated.push_back(pressure_val);
        AEL_data.stresscalculated.push_back(stress_tensor);
        AEL_data.structurecalculated.push_back(idVaspRun);
        // Print out stress tensor
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Shear stress tensor = " << AEL_data.shear_stress[i - 1].at(AEL_data.shear_stress[i - 1].size() - 1) << endl;
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", energy = " << outcar.energy_cell << endl;
        aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", pressure = " << outcar.pressure_residual << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        idVaspRun++;
      }
    }
    if (AEL_data.calcstrainorigin) {
      if (idVaspRun > vaspRuns.size()) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "idVaspRun = " << idVaspRun << " > vaspRuns.size() = " << vaspRuns.size() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        return 1;
      }
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Directory = " << vaspRuns[idVaspRun].Directory << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

      // If tarred and compressed directory exists...
      if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.bz2");
      } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.bz2");
      } // Extract all...
      if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.gz")) {
        aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.gz");
      } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.bz2");
      } // Extract all...
      if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.xz")) {
        aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.xz");
      } else if (aurostd::FileExist(dirrunname[idVaspRun] + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + dirrunname[idVaspRun] + ".tar.bz2");
      } // Extract all...

      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if (!aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + _AFLOWLOCK_) && !aurostd::FileExist(dirrunname[idVaspRun] + "/" + _AFLOWLOCK_) &&
          !(((XHOST.ARUN_POSTPROCESS || AEL_data.postprocess) && ((aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/agl.LOCK")) || (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/LOCK")) ||
                                                                  (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/ael.LOCK")) || (aurostd::FileExist(dirrunname[idVaspRun] + "/agl.LOCK")) ||
                                                                  (aurostd::FileExist(dirrunname[idVaspRun] + "/LOCK")) || (aurostd::FileExist(dirrunname[idVaspRun] + "/ael.LOCK")))))) {
        aurostd::StringstreamClean(aus);
        aus << _AELSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        throw AELStageBreak();
      }
      if (AEL_data.vasprunxmlstress) {
        for (size_t ij = 0; ij < vfile.size() && (vasprunxml.content.empty()); ij++) {
          if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[ij])) {
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns[idVaspRun].Directory + "/" + vfile[ij] << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            vasprunxml.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[ij]);
            vfilename = vfile[ij];
          } else if (aurostd::FileExist(dirrunname[idVaspRun] + "/" + vfile[ij])) {
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << dirrunname[idVaspRun] + "/" + vfile[ij] << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            vasprunxml.GetPropertiesFile(dirrunname[idVaspRun] + "/" + vfile[ij]);
            vfilename = vfile[ij];
          }
        }
        if (vasprunxml.content.empty()) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -vasprunxml.stress;
        vasprunxml.clear();
      } else {
        for (size_t ij = 0; ij < vfile.size() && (outcar.content.empty()); ij++) {
          if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[ij])) {
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << vaspRuns[idVaspRun].Directory + "/" + vfile[ij] << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            outcar.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[ij]);
            vfilename = vfile[ij];
          } else if (aurostd::FileExist(dirrunname[idVaspRun] + "/" + vfile[ij])) {
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << dirrunname[idVaspRun] + "/" + vfile[ij] << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            outcar.GetPropertiesFile(dirrunname[idVaspRun] + "/" + vfile[ij]);
            vfilename = vfile[ij];
          }
        }
        if (outcar.content.empty()) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -outcar.stress;
        energy_cell_val = outcar.energy_cell;
        pressure_val = outcar.pressure_residual;
        outcar.clear();
      }
      AEL_data.origin_stress.push_back(stress_tensor);
      AEL_data.energycalculated.push_back(energy_cell_val);
      AEL_data.pressurecalculated.push_back(pressure_val);
      AEL_data.stresscalculated.push_back(stress_tensor);
      AEL_data.structurecalculated.push_back(idVaspRun);
      // Print out stress tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Origin stress tensor = " << AEL_data.origin_stress[0] << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      idVaspRun++;
    } else if (AEL_data.fitrelaxedstruct) {
      const string relaxedstructdirname = AEL_data.dirpathname + "/";
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "Initial relaxed structure directory = " << relaxedstructdirname << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

      // If tarred and compressed directory exists...
      if (aurostd::FileExist(relaxedstructdirname + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + relaxedstructdirname + ".tar.bz2");
      } // Extract all...
      if (aurostd::FileExist(relaxedstructdirname + ".tar.gz")) {
        aurostd::execute(string("tar -xf ") + relaxedstructdirname + ".tar.gz");
      } // Extract all...
      if (aurostd::FileExist(relaxedstructdirname + ".tar.xz")) {
        aurostd::execute(string("tar -xf ") + relaxedstructdirname + ".tar.xz");
      } // Extract all...

      if (AEL_data.vasprunxmlstress) {
        for (size_t ij = 0; ij < vfile.size() && (vasprunxml.content.empty()); ij++) {
          if (aurostd::FileExist(relaxedstructdirname + "/" + vfile[ij])) {
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << relaxedstructdirname + "/" + vfile[ij] << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            vasprunxml.GetPropertiesFile(relaxedstructdirname + "/" + vfile[ij]);
            vfilename = vfile[ij];
          }
        }
        if (vasprunxml.content.empty()) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << relaxedstructdirname << " directory is missing." << endl;
          aus << _AELSTR_ERROR_ + "The flag [AFLOW_AEL]FITRELAXEDSTRUCT=ON is set in the input file." << endl;
          aus << _AELSTR_ERROR_ + "This requires the results of a relaxed calculation to be present in the directory " << relaxedstructdirname << endl;
          aus << _AELSTR_ERROR_ + "Either set the flag to OFF or run AEL in the appropriate directory." << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -vasprunxml.stress;
        vasprunxml.clear();
      } else {
        for (size_t ij = 0; ij < vfile.size() && (outcar.content.empty()); ij++) {
          if (aurostd::FileExist(relaxedstructdirname + "/" + vfile[ij])) {
            if (LVERBOSE) {
              aurostd::StringstreamClean(aus);
              aus << _AELSTR_MESSAGE_ + "vfile = " << vfile[ij] << endl;
              aus << _AELSTR_MESSAGE_ + "file = " << relaxedstructdirname + "/" + vfile[ij] << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            outcar.GetPropertiesFile(relaxedstructdirname + "/" + vfile[ij]);
            vfilename = vfile[ij];
          }
        }
        if (outcar.content.empty()) {
          aurostd::StringstreamClean(aus);
          aus << _AELSTR_WARNING_ + "The " << vfilename << " file in " << relaxedstructdirname << " directory is missing." << endl;
          aus << _AELSTR_ERROR_ + "The flag [AFLOW_AEL]FITRELAXEDSTRUCT=ON is set in the input file." << endl;
          aus << _AELSTR_ERROR_ + "This requires the results of a relaxed calculation to be present in the directory " << relaxedstructdirname << endl;
          aus << _AELSTR_ERROR_ + "Either set the flag to OFF or run AEL in the appropriate directory." << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          throw AELStageBreak();
        }
        stress_tensor = -outcar.stress;
        energy_cell_val = outcar.energy_cell;
        pressure_val = outcar.pressure_residual;
        outcar.clear();
      }
      AEL_data.origin_stress.push_back(stress_tensor);
      AEL_data.energycalculated.push_back(energy_cell_val);
      AEL_data.pressurecalculated.push_back(pressure_val);
      AEL_data.stresscalculated.push_back(stress_tensor);
      AEL_data.structurecalculated.push_back(idVaspRun);
      // Print out stress tensor
      aurostd::StringstreamClean(aus);
      aus << _AELSTR_MESSAGE_ + "System number = " << idVaspRun << ", Origin stress tensor = " << AEL_data.origin_stress[0] << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      idVaspRun++;
    }
    return 0;
  }
} // namespace AEL_functions

// **************************************************************************
//  End of AFLOW AEL set-up and extract stress-strain data
// **************************************************************************

#endif // _AFLOW_AEL_GET_STRESS_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
