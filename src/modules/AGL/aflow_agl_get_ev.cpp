// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu
#ifndef _AFLOW_AGL_GET_EV_CPP
#define _AFLOW_AGL_GET_EV_CPP

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
#include "flow/aflow_ivasp.h"
#include "flow/aflow_kvasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/AGL/aflow_agl_debye.h"

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
//                  AFLOW Automatic GIBBS Library (AGL) (2013-2021)
// ###############################################################################
//
// Uses quasi-harmonic Debye model to obtain thermodynamic properties of materials
// Based on original Fortran program written by M. A. Blanco et al.
// See Computer Physics Communications 158, 57-72 (2004) and Journal of Molecular Structure (Theochem) 368, 245-255 (1996) for details of original GIBBS program
// See C. Toher et al., Phys. Rev. B 90, 174107 (2014), Phys. Rev. 1, 015401 (2017) and references therein for description of this AGL implementation
// Please cite these works in addition to the general AFLOW papers if you use results generated using AGL
//

// *****************************************************************************************************************
// The following functions are for setting up AGL inputs for postprocessing runs called from other parts of AFLOW
// *****************************************************************************************************************
namespace AGL_functions {
  uint AGL_xvasp_flags_populate(_xvasp& xvasp, string& AflowIn, string& AflowInName, string& FileLockName, const string& directory_LIB, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) {
    ifstream FileAFLOWIN;
    string FileNameAFLOWIN;
    ostringstream aus;
    const vector<string> vAflowInCheck;
    bool agl_aflowin_found = false;
    bool Krun = true;
    const bool load_POSCAR_from_xvasp = false;
    // Set aflags
    aflags.Directory = directory_LIB;
    if (aflags.Directory.at(0) != '/' && aflags.Directory.at(0) != '.' && aflags.Directory.at(0) != ' ') {
      aflags.Directory = "./" + aflags.Directory;
    }
    aflags.KBIN_RUN_AFLOWIN = true;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN = false;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP = false;
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN = false;
    aflags.KBIN_DELETE_AFLOWIN = false;

    // CT20200721 Calls function AGL_Get_AflowInName to find correct aflow.in filename
    AGL_functions::AGL_Get_AflowInName(AflowInName, directory_LIB, agl_aflowin_found);
    FileNameAFLOWIN = directory_LIB + "/" + AflowInName;

    // CO20200502 STOP - CT, I am consolidating the following code with an outer loop, it should make it easier to patch in the future
    if (agl_aflowin_found) {
      // Set FileMESSAGE name
      // CT20200624 Moved down so LOCK file is only renamed if this is actually an AGL main directory and not an ARUN.AGL directory
      if (!FileLockName.empty()) {
        if (aurostd::FileExist(directory_LIB + "/" + FileLockName)) {
          aurostd::file2file(aurostd::CleanFileName(directory_LIB + "/" + FileLockName), aurostd::CleanFileName(directory_LIB + "/" + FileLockName + ".run"));
        }
        const string FileNameMessage = directory_LIB + "/" + FileLockName;
        FileMESSAGE.open(FileNameMessage.c_str(), std::ios::app);
      } else {
        if (aurostd::FileExist(directory_LIB + "/agl.LOCK")) {
          aurostd::file2file(aurostd::CleanFileName(directory_LIB + "/agl.LOCK"), aurostd::CleanFileName(directory_LIB + "/agl.LOCK.run"));
        }
        const string FileNameMessage = directory_LIB + "/agl.LOCK";
        FileMESSAGE.open(FileNameMessage.c_str(), std::ios::app);
      }
      // Search for AGL aflow.in filename
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "AFLOW Input file name = " << FileNameAFLOWIN << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
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
      xvasp.clear();
      const uint ixvasp = 0;
      xvasp.POSCAR_index = ixvasp;
      KBIN::readModulesFromAflowIn(AflowIn, kflags, xvasp);
      xvasp.Directory = aflags.Directory;
      if (Krun) {
        Krun = (Krun && KBIN::VASP_Produce_INPUT(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags, load_POSCAR_from_xvasp));
      }
      if (Krun) {
        Krun = (Krun && KBIN::VASP_Modify_INPUT(xvasp, FileMESSAGE, aflags, kflags, vflags));
      }
      // Fix blank species
      if (!xvasp.str.species.empty()) {
        if (xvasp.str.species.at(0).empty()) {
          pflow::fixEmptyAtomNames(xvasp.str);
        }
      }
      if (Krun) {
        return 0;
      } else {
        return 1;
      }
    } else {
      cerr << _AGLSTR_MESSAGE_ << "AGL input file not found!" << endl; // CT20200624 If no AGL file present, then this is not an AGL main directory: write to cerr and return 2 to inform calling function
      cerr << _AGLSTR_MESSAGE_ << "Not an AGL main directory" << endl;
      return 2;
    }
  }
} // namespace AGL_functions

// *****************************************************************************************************************
// Finds aflow.in file for AGL calculation (if it exists) //CT20200713
// *****************************************************************************************************************
namespace AGL_functions {
  uint AGL_Get_AflowInName(string& AflowInName, const string& directory_LIB, bool& agl_aflowin_found) {
    ifstream FileAFLOWINcheck;
    string FileNameAFLOWIN;
    string FileNameAFLOWINcheck;
    string AflowInCheck;
    string aflowinname;
    string stmp;
    uint aglerror = 0;
    uint filelength = 0;
    vector<string> vAflowInCheck;
    agl_aflowin_found = false;
    vector<string> vaflowins;
    if (!AflowInName.empty()) {
      vaflowins.push_back(AflowInName);
    } // Check if AflowInName exists
    if (!_AFLOWIN_.empty()) {
      vaflowins.push_back(_AFLOWIN_);
    } // Otherwise, check if _AFLOWIN_ file is AGL input file
    vaflowins.push_back(_AFLOWIN_AGL_DEFAULT_); // Otherwise, check for other commonly used names for AGL aflow.in file
    for (size_t iaf = 0; iaf < vaflowins.size() && !agl_aflowin_found; iaf++) {
      aflowinname = vaflowins.at(iaf);
      if (aurostd::CompressFileExist(directory_LIB + "/" + aflowinname, stmp) && aurostd::IsCompressed(stmp)) {
        aurostd::DecompressFile(stmp);
      } // CO20210204 - fix aflow.in.xz
      if ((!agl_aflowin_found) && (aurostd::FileExist(directory_LIB + "/" + aflowinname))) {
        FileNameAFLOWINcheck = directory_LIB + "/" + aflowinname;
        AflowInCheck = "";
        // READ aflowinname and put into AflowInCheck
        filelength = aurostd::file2string(FileNameAFLOWINcheck, AflowInCheck);
        if (filelength > 0) {
          aglerror = 0;
        } else {
          aglerror = 1;
        }
        AflowInCheck = aurostd::RemoveComments(AflowInCheck); // NOW Clean AFLOWIN
        vAflowInCheck.clear();
        aurostd::string2vectorstring(AflowInCheck, vAflowInCheck);
        // Check if aflowinname contains command to run AGL
        for (size_t i = 0; i < vAflowInCheck.size() && !agl_aflowin_found; i++) {
          if ((aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AGL]CALC", true) || aurostd::substring2bool(AflowInCheck, "[VASP_AGL]CALC", true)) &&
              !(aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AGL]CALC_", true) || aurostd::substring2bool(vAflowInCheck[i], "[VASP_AGL]CALC_", true) ||
                aurostd::substring2bool(vAflowInCheck[i], "[AFLOW_AGL]CALCS", true) || aurostd::substring2bool(vAflowInCheck[i], "[VASP_AGL]CALCS", true) || false)) {
            FileNameAFLOWIN = FileNameAFLOWINcheck;
            agl_aflowin_found = true;
            AflowInName = aflowinname;
          }
        }
        FileAFLOWINcheck.close();
      }
    }
    return aglerror;
  }
} // namespace AGL_functions

// *******************************************************************************
// The following functions are for generating _AFLOWIN_ files
// *******************************************************************************

// ***************************************************************************
// AGL_functions::aglvaspflags
// ***************************************************************************
namespace AGL_functions {
  //
  // Function to assign values for VASP input flags from _AFLOWIN_ file to vaspRun _xvasp class
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint aglvaspflags(_xvasp& vaspRun, _vflags& vaspFlags, _kflags& kbinFlags, string& dirrunname, _AGL_data& AGL_data, ofstream& FileMESSAGE) {
    const bool LVERBOSE = (false || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile;
    vector<string> dfile;
    string vfilename;
    string dfilename;
    string ffilename;
    bool vfileexist = false;
    bool skipdir = false;
    aurostd::string2tokens(dirrunname, dfile, "/");
    dfilename = dfile.at(dfile.size() - 1);
    if (AGL_data.relax_static || AGL_data.static_only) {
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
    // SOME WARNINGS: check existence of LOCK and OUTCAR.static files
    if (!(aurostd::FileExist(dirrunname + "/" + _AFLOWLOCK_) || ((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/LOC"
                                                                                                                                                                                                      "K")))) &&
        (vfileexist)) {
      aurostd::StringstreamClean(aus);
      // [OBOLSETE] aus << _AGLSTR_WARNING_ + "found OUTCAR.static but no LOCK in " <<  vaspRun.Directory << endl;
      aus << _AGLSTR_WARNING_ + "found " << vfilename << " but no " << _AFLOWLOCK_ << " in " << dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (size_t ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list[ij];
        if (LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (aurostd::substring2bool(dfilename, ffilename, true)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          skipdir = true;
        }
      }
      if (skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Directory: " << dfilename << " will be skipped." << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        skipdir = false;
      } else {
        return 1;
      }
    }

    if ((aurostd::FileExist(dirrunname + "/" + _AFLOWLOCK_) || ((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) && (aurostd::FileExist(dirrunname + "/agl.LOCK") || aurostd::FileExist(dirrunname + "/LOC"
                                                                                                                                                                                                     "K")))) &&
        !(vfileexist)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "found " << _AFLOWLOCK_ << " but no OUTCAR in " << dirrunname << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (size_t ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list[ij];
        if (LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (aurostd::substring2bool(dfilename, ffilename, true)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          skipdir = true;
        }
      }
      if (skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Directory: " << dfilename << " will be skipped." << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        skipdir = false;
      } else {
        return 1;
      }
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

    if (AGL_data.relax_static) {
      vaspRun.AVASP_flag_RUN_RELAX_STATIC = true;
      *RUN_map["RELAX_STATIC"] = true;
      vaspRun.aopts.flag("FLAG::VOLUME_PRESERVED", true);
    } else {
      vaspRun.AVASP_flag_RUN_STATIC = true;
      *RUN_map["STATIC"] = true;
      // Increase DOSCAR density for STATIC-only runs
      vaspRun.aopts.flag("FLAG::EXTRA_INCAR", true);
      vaspRun.AVASP_EXTRA_INCAR << "# Added by [AFLOW_AGL] begin" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "EMIN= -30.0    # For finer DOS grid" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "EMAX=  45.0    # For finer DOS grid" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "NEDOS= 5000    # For finer DOS grid" << std::endl;
      vaspRun.AVASP_EXTRA_INCAR << "# Added by [AFLOW_AGL] end" << std::endl;
    }
    for (const auto& [name, flag] : RUN_map) {
      vaspFlags.KBIN_VASP_RUN.flag(name, *flag);
    }

    if (AGL_data.precaccalgonorm) {
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.isentry = true;
      vaspFlags.KBIN_VASP_FORCE_OPTION_PREC.content_string = "ACCURATE";
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.clear();
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.isentry = true;
      vaspFlags.KBIN_VASP_FORCE_OPTION_ALGO.content_string = "NORMAL";
    }

    // Change format of POSCAR
    if ((!kbinFlags.KBIN_MPI && (kbinFlags.KBIN_BIN.find("46") != string::npos)) || (kbinFlags.KBIN_MPI && (kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
      vaspRun.str.is_vasp5_poscar_format = false;
    }
    return 0;
  }
} // namespace AGL_functions

// ************************************************************************************************
// This set of functions extract, sort and check (E, V) data from VASP runs
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::extractenerg
// ***************************************************************************
namespace AGL_functions {
  //
  // extractenerg: Extract final energies from the completed VASP calculations
  // Adapted from section of AFLOW APL function DirectMethodPC::runVASPCalculations()
  //
  uint extractenerg(vector<_xvasp>& vaspRuns, _AGL_data& AGL_data, vector<string>& dirrunname, ofstream& FileMESSAGE) {
    const bool LVERBOSE = (false || XHOST.DEBUG);
    ostringstream aus;
    vector<string> vfile;
    vector<string> dfile;
    aurostd::string2tokens(string("OUTCAR.static.bz2,OUTCAR.static.gz,OUTCAR.static.xz,OUTCAR.static"), vfile, ",");
    xOUTCAR outcar;
    string dfilename;
    string ffilename;
    bool skipdir = false;
    aurostd::xmatrix<double> stress_tensor(3, 3);
    AGL_data.volumeinput.clear();
    AGL_data.energyinput.clear();
    AGL_data.pressurecalculated.clear();
    AGL_data.stresscalculated.clear();
    for (size_t idVaspRun = 0; idVaspRun < vaspRuns.size(); idVaspRun++) {
      skipdir = false;
      aurostd::StringstreamClean(aus);
      // Print out total energy
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      aurostd::string2tokens(vaspRuns[idVaspRun].Directory, dfile, "/");
      dfilename = dfile.at(dfile.size() - 1);
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "Directory name = " << dfilename << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      // Check if structure is on list of failed runs to be skipped
      // If so, then skip reading and continue to next structure
      for (size_t ij = 0; ij < AGL_data.failed_arun_list.size(); ij++) {
        ffilename = AGL_data.failed_arun_list[ij];
        if (LVERBOSE) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "dfilename = " << dfilename << endl;
          aus << _AGLSTR_MESSAGE_ + "ffilename = " << ffilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (aurostd::substring2bool(dfilename, ffilename, true)) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Found directory in to-skip list: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          skipdir = true;
        }
      }
      if (skipdir) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        skipdir = false;
        continue;
      }

      // If tarred and compressed directory exists...
      if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.bz2");
      } else if (aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.bz2")) {
        aurostd::execute(string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.bz2");
      } // Extract all...
      if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.gz")) {
        aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.gz");
      } else if (aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.gz")) {
        aurostd::execute(string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.gz");
      } // Extract all...
      if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + ".tar.xz")) {
        aurostd::execute(string("tar -xf ") + vaspRuns[idVaspRun].Directory + ".tar.xz");
      } else if (aurostd::FileExist(dirrunname.at(idVaspRun) + ".tar.xz")) {
        aurostd::execute(string("tar -xf ") + dirrunname.at(idVaspRun) + ".tar.xz");
      } // Extract all...

      // If the LOCK file is missing, then it is probably a corrupted run
      // Do not accept it and wait for the new run
      if (!aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + _AFLOWLOCK_) && !aurostd::FileExist(dirrunname.at(idVaspRun) + "/" + _AFLOWLOCK_) &&
          !((XHOST.ARUN_POSTPROCESS || AGL_data.postprocess) &&
            ((aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/agl.LOCK")) || (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/LOCK")) || (aurostd::FileExist(dirrunname.at(idVaspRun) + "/agl.LOCK")) ||
             (aurostd::FileExist(dirrunname.at(idVaspRun) + "/LOCK"))))) { // CT20200625 Modify to check for other LOCK file names if postprocessing run
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "The " << _AFLOWLOCK_ << " file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        if (AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }

      // for(size_t i=0;i<vfile.size()&&(outcar.outcar=="");i++)
      for (size_t i = 0; i < vfile.size() && (outcar.content.empty()); i++) {
        if (aurostd::FileExist(vaspRuns[idVaspRun].Directory + "/" + vfile[i])) {
          outcar.GetPropertiesFile(vaspRuns[idVaspRun].Directory + "/" + vfile[i]);
        } else if (aurostd::FileExist(dirrunname.at(idVaspRun) + "/" + vfile[i])) {
          outcar.GetPropertiesFile(dirrunname.at(idVaspRun) + "/" + vfile[i]);
        }
      }
      // if(outcar.outcar=="")
      if (outcar.content.empty()) {
        aurostd::StringstreamClean(aus);
        aus << _AGLSTR_WARNING_ + "The OUTCAR.static file in " << vaspRuns[idVaspRun].Directory << " directory is missing." << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        if (AGL_data.autoskipfailedaruns) {
          aurostd::StringstreamClean(aus);
          aus << _AGLSTR_MESSAGE_ + "Skipping directory: " << dfilename << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          continue;
        } else {
          throw AGLStageBreak();
        }
      }

      AGL_data.energyinput.push_back(outcar.energy_cell);
      AGL_data.volumeinput.push_back(vaspRuns[idVaspRun].str.Volume());
      AGL_data.pressurecalculated.push_back(outcar.pressure_residual);
      stress_tensor = -outcar.stress;
      AGL_data.stresscalculated.push_back(stress_tensor);
      AGL_data.structurecalculated.push_back(idVaspRun);

      aurostd::StringstreamClean(aus);
      // Print out total energy, volume and calculated residual pressure
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Energy (eV) = " << AGL_data.energyinput.at(AGL_data.energyinput.size() - 1) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Volume (Ang^3) = " << AGL_data.volumeinput.at(AGL_data.volumeinput.size() - 1) << endl;
      aus << _AGLSTR_MESSAGE_ << "System number = " << idVaspRun << ", Pressure (kB) = " << AGL_data.pressurecalculated.at(AGL_data.pressurecalculated.size() - 1) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      outcar.clear();
    }
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::checkmin
// ***************************************************************************
namespace AGL_functions {
  //
  // checkmin: checks position of lowest energy in (E, V)
  // If lowest energy corresponds to smallest volume, sets cmerr = 1
  // If lowest energy corresponds to largest volume, sets cmerr = 2
  // Otherwise, sets cmerr = 0
  //
  uint checkmin(_AGL_data& AGL_data, int& cmerr, ofstream& FileMESSAGE) {
    ostringstream aus;
    vector<double> energy(AGL_data.energyinput.size());
    vector<double> volume(AGL_data.volumeinput.size());
    cmerr = 0;
    uint aglerror = 0;

    for (size_t i = 0; i < AGL_data.energyinput.size(); i++) {
      energy.at(i) = AGL_data.energyinput[i];
      volume.at(i) = AGL_data.volumeinput.at(i);
    }
    // Sorts (E, V) data into order of increasing volume
    // AGL_functions::qcksort (vol, idx, 0, is-1, FileMESSAGE);
    aglerror = AGL_functions::qcksortev(volume, energy, FileMESSAGE);
    if (aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V)" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return aglerror;
    }
    // Find global minimum of energy data points
    double etref = energy.at(0);
    uint itmin = 0;
    for (size_t i = 0; i < AGL_data.energyinput.size(); i++) {
      if (energy.at(i) < etref) {
        etref = energy.at(i);
        itmin = i;
      }
    }
    // Check that the minimum energy does not correspond to the largest or smallest volume
    // If it does, then this suggests that a larger or smaller volume may have a lower energy
    if (itmin == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Minimum energy is for smallest volume" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      cmerr = 1;
      return aglerror;
    } else if (itmin == (AGL_data.energyinput.size() - 1)) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Minimum energy is for largest volume" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      cmerr = 2;
      return aglerror;
    } else {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "Energy minimum is contained in (E, V) data" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      cmerr = 0;
      return aglerror;
    }
    energy.clear();
    volume.clear();
    return aglerror;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::checkconcav
// ***************************************************************************
namespace AGL_functions {
  //
  // checkconcav: checks concavity of (E, V) data
  // If data around global minimum of supplied (E, V) data is not concave, sets ccerr = 1
  // Otherwise, sets ccerr = 0
  //
  uint checkconcav(_AGL_data& AGL_data, int& ccerr, ofstream& FileMESSAGE) {
    ostringstream aus;
    uint j = 1;
    vector<double> energy(AGL_data.energyinput.size());
    vector<double> volume(AGL_data.volumeinput.size());
    ccerr = 0;
    uint aglerror = 0;
    for (size_t i = 0; i < AGL_data.energyinput.size(); i++) {
      energy.at(i) = AGL_data.energyinput[i];
      volume.at(i) = AGL_data.volumeinput.at(i);
    }
    // Sort (E, V) data in order of increasing volume
    aglerror = AGL_functions::qcksortev(volume, energy, FileMESSAGE);
    if (aglerror != 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_ERROR_ + "Failed to sort E(V)" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return aglerror;
    }
    // Find global minimum of energy data points
    double etref = energy.at(0);
    uint itmin = 0;
    for (size_t i = 0; i < AGL_data.energyinput.size(); i++) {
      if (energy.at(i) < etref) {
        etref = energy.at(i);
        itmin = i;
      }
    }
    // Finds first acceptable point (first point where E-V curve is concave)
    while (((energy.at(j) - energy.at(j - 1)) / (volume.at(j) - volume.at(j - 1)) >= (energy.at(j) - energy.at(j + 1)) / (volume.at(j) - volume.at(j + 1))) && j < AGL_data.energyinput.size() - 2) {
      j = j + 1;
    }
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "Concavity check: j = " << j << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    // If only the last three points are concave, then the (E, V) data cannot be used for GIBBS
    // A warning is given and the signal is given to rerun the VASP calculations with more k-points
    if (j == AGL_data.energyinput.size() - 1) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "All points show convex patterns" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    // If the first accepted point is already past the global minimum, the (E, V) could cause problems for GIBBS
    // Gives a warning and sends the signal to rerun the VASP calculations with more k-points
    if (j >= itmin) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "First concave point is already passed the global minimum" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    j = j + 1;
    uint jtmax = 2;

    // Point j marks the last accepted point, i the new trial point
    for (size_t i = j + 1; i <= AGL_data.energyinput.size() - 1; i++) {
      if ((energy.at(j) - energy.at(j - 1)) / (volume.at(j) - volume.at(j - 1)) < (energy.at(j) - energy.at(i)) / (volume.at(j) - volume.at(i))) {
        j = j + 1;
        jtmax = i;
      }
    }

    // If the global minimum lies outside of the range of the accepted points, then there will be problems using the (E, V) data for GIBBS
    // Gives a warning and then gives the signal to rerun the VASP calculations with more k-points
    // Problems with noise in the (E, V) data are often caused by insufficient K-points or basis set
    if (jtmax <= itmin) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_WARNING_ + "Global minimum of (E, V) data lies outside of initial range of points accepted for concavity" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      ccerr = 1;
      return aglerror;
    }
    // If the data is concave and includes the global minimum, then data should be good for GIBBS
    // Gives success message and sends signal to continue to GIBBS method
    aurostd::StringstreamClean(aus);
    aus << _AGLSTR_MESSAGE_ << "(E, V) data is concave" << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    ccerr = 0;
    return aglerror;
  }
} // namespace AGL_functions

// ************************************************************************************************
//  This set of functions sorts the (E, V) data, fits it by a polynomial, and finds the minimum
// ************************************************************************************************

// ***************************************************************************
// AGL_functions::qcksortev
// ***************************************************************************
namespace AGL_functions {
  //
  // Quick-sort algorithm: sorts the elements of array in ascending order
  // Calls aurostd::quicksort to actually sort data
  // Returns vol and energ sorted in order of increasing vol
  //
  uint qcksortev(vector<double>& vol, vector<double>& energ, ofstream& FileMESSAGE) {
    int icheck = 0;
    ostringstream aus;
    // Check if data is already in correct order
    for (size_t i = 0; i < (vol.size() - 1); i++) {
      if (vol.at(i + 1) < vol.at(i)) {
        icheck = 1;
      }
    }
    // If data is already in correct order, exits function without sorting
    if (icheck == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ + "qcksort: Data is already arranged in increasing order" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 0;
    }
    // If data is not in the correct order, calls quicksort2 from aurostd library to sort it
    // Sorts both vol and energ so that they are in the order of increasing vol
    aurostd::sort(vol, energ);
    return 0;
  }
} // namespace AGL_functions

// ***************************************************************************
// AGL_functions::qcksortevt
// ***************************************************************************
namespace AGL_functions {
  //
  // Quick-sort algorithm: sorts the elements of array in ascending order
  // Calls aurostd::quicksort to actually sort data
  // Returns vol, energ and tdebye sorted in order of increasing vol
  //
  uint qcksortevt(vector<double>& vol, vector<double>& energ, vector<double>& tdebye, ofstream& FileMESSAGE) {
    int icheck = 0;
    ostringstream aus;
    // Check if data is already in correct order
    for (size_t i = 0; i < (vol.size() - 1); i++) {
      if (vol.at(i + 1) < vol.at(i)) {
        icheck = 1;
      }
    }
    // If data is already in correct order, exits function without sorting
    if (icheck == 0) {
      aurostd::StringstreamClean(aus);
      aus << _AGLSTR_MESSAGE_ << "qcksort: Data is already arranged in increasing order" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      return 0;
    }
    // If data is not in the correct order, calls quicksort2 from aurostd library to sort it
    // Sorts both vol and energ so that they are in the order of increasing vol
    aurostd::sort(vol, energ, tdebye);
    return 0;
  }
} // namespace AGL_functions

// **************************************************************************
//  End of AFLOW AGL set-up, extract, sort and check (E, V) data
// **************************************************************************

#endif // _AFLOW_AGL_GET_EV_CPP
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                Aflow CORMAC TOHER - Duke University 2013-2021           *
// *                                                                         *
// ***************************************************************************
