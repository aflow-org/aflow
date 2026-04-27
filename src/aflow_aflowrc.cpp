
#include "aflow_aflowrc.h"

#include "config.h"

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xfile.h"

#include "aflow_init.h" // NOLINT // ST: Needed for stuff behind aflow_aflowrc.h macros. We can't put aflow_init.h in it due to cyclic include problems.
#include "aflow_xhost.h"
#include "flow/aflow_pflow.h"

using std::endl;
using std::ifstream;
using std::iostream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

// ***************************************************************************
// aflowrc::load_default
// ***************************************************************************
namespace aflowrc {
  bool load_default(string schema, string schema_default) {
    bool found = false;
    string aus;
    string string_to_add = schema_default;
    vector<string> tokens;
    for (size_t i = 0; i < XHOST.vaflowrc.size() && !found; i++) {
      aurostd::string2tokens(XHOST.vaflowrc.at(i), tokens, "=");
      if (!tokens.empty() && !found) {
        if (aurostd::RemoveWhiteSpaces(tokens.at(0)) == schema && !found) {
          // CO20181226 - it is possible to have '=' inside value: MPI_START_DEFAULT="export OMP_NUM_THREADS=1"
          // treat tokens[0] as special
          tokens.erase(tokens.begin()); // remove key
          found = true;
          aus = aurostd::RemoveWhiteSpacesFromTheBack(aurostd::joinWDelimiter(tokens, "=")); // CO20180705 - if there are spaces between ="" and //, then the value is set to the spaces (not an empty string!) //CO20181226 - join again by '='
          aurostd::string2tokens(aus, tokens, "\""); //	    if(tokens.size()>0) cerr << tokens.at(0) << endl;
          if (!tokens.empty()) {
            string_to_add = tokens.at(0);
          }
        }
      }
    }
    // fix ~/ with XHOST.user
    if (aurostd::substring2bool(string_to_add, "~/")) {
      aurostd::StringSubst(string_to_add, "~/", XHOST.home + "/");
    }
    XHOST.adefault.push_attached(schema, string_to_add); // add what is present or the default if not present
    return found;
  }
  template <class utype> bool load_default(string schema, utype schema_default) {
    const bool found = XHOST.adefault.args2addattachedscheme(XHOST.vaflowrc, schema, string(schema + "="), ""); // add what is present
    if (!found) {
      XHOST.adefault.push_attached(schema, aurostd::utype2string<utype>(schema_default)); // add default if not present
    }
    return found;
  }
} // namespace aflowrc

// ***************************************************************************
// aflowrc::is_available
// ***************************************************************************
namespace aflowrc {
  bool is_available(std::ostream& oss, bool AFLOWRC_VERBOSE) {
    const bool LDEBUG = (false || XHOST.DEBUG || AFLOWRC_VERBOSE);
    bool aflowrc_local = false;
    bool aflowrc_global = false;
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;
    }

    // If a user installs AFLOW using snap https://snapcraft.io/aflow, it cannot access hidden files in the home folder.
    // Therefore, we need to verify if we are running in a snap container and modify the location of the configuration file accordingly.
    // The configuration file should be located under `~/snap/aflow/common/aflow.rc`.
    // Note that it is not a hidden file because it is a specific folder dedicated to the aflow snap.
    if (const char* snap_p = std::getenv("SNAP_USER_COMMON")) {
      XHOST.aflowrc_filename = std::string(snap_p) + "/aflow.rc";
      return aurostd::FileExist(XHOST.aflowrc_filename);
    }

    // TESTING LOCAL OR USER BASED
    if (XHOST.aflowrc_filename.empty()) {
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_LOCAL;
    }
    aflowrc_local = aurostd::FileExist(AFLOWRC_FILENAME_LOCAL);
    aflowrc_global = aurostd::FileExist(AFLOWRC_FILENAME_GLOBAL);

    // LOCAL=true && GLOBAL=true => take LOCAL
    if (aflowrc_local && aflowrc_global) {
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " LOCAL=true && GLOBAL=true => LOCAL " << endl;
      }
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_LOCAL;
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " END" << endl;
      }
      return true;
    }
    // LOCAL=true && GLOBAL=false => take LOCAL
    if (aflowrc_local && !aflowrc_global) {
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " LOCAL=true && GLOBAL=false => LOCAL " << endl;
      }
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_LOCAL;
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " END" << endl;
      }
      return true;
    }
    // LOCAL=false && GLOBAL=true => take GLOBAL
    if (!aflowrc_local && aflowrc_global) {
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " LOCAL=false && GLOBAL=true => GLOBAL " << endl;
      }
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_GLOBAL;
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " END" << endl;
      }
      return true;
    }
    // LOCAL=false && GLOBAL=false => take NOTHING AND REWRITE
    if (!aflowrc_local && !aflowrc_global) {
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " LOCAL=false && GLOBAL=false => NOTHING " << endl;
      }
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_LOCAL; // because it is going to write it
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      }
      if (LDEBUG) {
        oss << __AFLOW_FUNC__ << " END" << endl;
      }
      return false;
    }

    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " END" << endl;
    }
    return false;
  }
} // namespace aflowrc

// ***************************************************************************
// aflowrc::read
// ***************************************************************************
namespace aflowrc {
  bool read(std::ostream& oss, bool AFLOWRC_VERBOSE) {
    const bool LDEBUG = (false || XHOST.DEBUG || AFLOWRC_VERBOSE);
    stringstream message; // CO20200404
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;
    }
    if (XHOST.aflowrc_filename.empty()) {
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_LOCAL;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
    }

    if (!aflowrc::is_available(oss, AFLOWRC_VERBOSE)) {
      if (!XHOST.vflag_control.flag("WWW")) { // CO20200404 - new web flag
        if (!(aurostd::substring2bool(XHOST.aflowrc_filename, "/mnt/MAIN") || aurostd::substring2bool(XHOST.aflowrc_filename, "/mnt/uMAIN"))) { // CO20200404 - patching for new disk
          message << XHOST.aflowrc_filename << " not found, loading DEFAULT values";
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_MESSAGE_); // CO20200404 - LEAVE std::cerr here, FR needs this for web
        }
      }
    }

    aurostd::file2string(XHOST.aflowrc_filename, XHOST.aflowrc_content);
    // oss << "BEGIN" << endl << XHOST.aflowrc_content << "END" << endl;
    // XHOST.aflowrc_content=aurostd::RemoveComments(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    //   XHOST.aflowrc_content=aurostd::RemoveWhiteSpaces(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    XHOST.aflowrc_content = aurostd::RemoveComments(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    // oss << "BEGIN" << endl << XHOST.aflowrc_content << "END" << endl;
    aurostd::string2vectorstring(XHOST.aflowrc_content, XHOST.vaflowrc); // vectorize

    // DEFAULT DEFINITIONS
    aflowrc::load_default("DEFAULT_KZIP_BIN", AFLOWRC_DEFAULT_KZIP_BIN);
    aflowrc::load_default("DEFAULT_KZIP_EXT", AFLOWRC_DEFAULT_KZIP_EXT);
    aflowrc::load_default("DEFAULT_TMPFS_DIRECTORIES", AFLOWRC_DEFAULT_TMPFS_DIRECTORIES);

    // HE20220218 START
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE", AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_AFLUX_SERVER", AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_AFLUX_PATH", AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_RESTAPI_SERVER", AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_RESTAPI_PATH", AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_FS_PATH", AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH);
    // HE20220218 STOP

    // ME20191001 START
    //  AFLOW database files
    aflowrc::load_default("DEFAULT_AFLOW_DB_FILE", AFLOWRC_DEFAULT_AFLOW_DB_FILE);
    aflowrc::load_default("DEFAULT_AFLOW_DB_STATS_FILE", AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE);
    aflowrc::load_default("DEFAULT_AFLOW_DB_DATA_PATH", AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH);
    aflowrc::load_default("DEFAULT_AFLOW_DB_LOCK_FILE", AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE);
    aflowrc::load_default("DEFAULT_AFLOW_DB_STALE_THRESHOLD", AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD);
    // ME20191001 END
    //  FILENAMES FOR AFLOW.ORG ANALYSIS
    aflowrc::load_default("DEFAULT_FILE_AFLOWLIB_ENTRY_OUT", AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
    aflowrc::load_default("DEFAULT_FILE_AFLOWLIB_ENTRY_JSON", AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_JSON);
    aflowrc::load_default("DEFAULT_FILE_EDATA_ORIG_OUT", AFLOWRC_DEFAULT_FILE_EDATA_ORIG_OUT);
    aflowrc::load_default("DEFAULT_FILE_EDATA_RELAX_OUT", AFLOWRC_DEFAULT_FILE_EDATA_RELAX_OUT);
    aflowrc::load_default("DEFAULT_FILE_EDATA_BANDS_OUT", AFLOWRC_DEFAULT_FILE_EDATA_BANDS_OUT);
    aflowrc::load_default("DEFAULT_FILE_DATA_ORIG_OUT", AFLOWRC_DEFAULT_FILE_DATA_ORIG_OUT);
    aflowrc::load_default("DEFAULT_FILE_DATA_RELAX_OUT", AFLOWRC_DEFAULT_FILE_DATA_RELAX_OUT);
    aflowrc::load_default("DEFAULT_FILE_DATA_BANDS_OUT", AFLOWRC_DEFAULT_FILE_DATA_BANDS_OUT);
    aflowrc::load_default("DEFAULT_FILE_EDATA_ORIG_JSON", AFLOWRC_DEFAULT_FILE_EDATA_ORIG_JSON);
    aflowrc::load_default("DEFAULT_FILE_EDATA_RELAX_JSON", AFLOWRC_DEFAULT_FILE_EDATA_RELAX_JSON);
    aflowrc::load_default("DEFAULT_FILE_EDATA_BANDS_JSON", AFLOWRC_DEFAULT_FILE_EDATA_BANDS_JSON);
    aflowrc::load_default("DEFAULT_FILE_DATA_ORIG_JSON", AFLOWRC_DEFAULT_FILE_DATA_ORIG_JSON);
    aflowrc::load_default("DEFAULT_FILE_DATA_RELAX_JSON", AFLOWRC_DEFAULT_FILE_DATA_RELAX_JSON);
    aflowrc::load_default("DEFAULT_FILE_DATA_BANDS_JSON", AFLOWRC_DEFAULT_FILE_DATA_BANDS_JSON);
    aflowrc::load_default("DEFAULT_FILE_TIME_OUT", AFLOWRC_DEFAULT_FILE_TIME_OUT);
    aflowrc::load_default("DEFAULT_FILE_SPACEGROUP1_OUT", AFLOWRC_DEFAULT_FILE_SPACEGROUP1_OUT);
    aflowrc::load_default("DEFAULT_FILE_SPACEGROUP2_OUT", AFLOWRC_DEFAULT_FILE_SPACEGROUP2_OUT);
    aflowrc::load_default("DEFAULT_FILE_VOLDISTPARAMS_OUT", AFLOWRC_DEFAULT_FILE_VOLDISTPARAMS_OUT);
    aflowrc::load_default("DEFAULT_FILE_VOLDISTEVOLUTION_OUT", AFLOWRC_DEFAULT_FILE_VOLDISTEVOLUTION_OUT);

    // FILENAMES FOR AFLOW OPERATION
    aflowrc::load_default("DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT", AFLOWRC_DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PRESCRIPT_OUT", AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PRESCRIPT_COMMAND", AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_COMMAND);
    aflowrc::load_default("DEFAULT_AFLOW_POSTSCRIPT_OUT", AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_POSTSCRIPT_COMMAND", AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_COMMAND);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_OUT", AFLOWRC_DEFAULT_AFLOW_PGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_JSON", AFLOWRC_DEFAULT_AFLOW_PGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_XTAL_OUT", AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_XTAL_JSON", AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT", AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT); // DX20200129
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON", AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON); // DX20200129
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_OUT", AFLOWRC_DEFAULT_AFLOW_PGROUPK_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_JSON", AFLOWRC_DEFAULT_AFLOW_PGROUPK_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_XTAL_OUT", AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_XTAL_JSON", AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_FGROUP_OUT", AFLOWRC_DEFAULT_AFLOW_FGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FGROUP_JSON", AFLOWRC_DEFAULT_AFLOW_FGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_SGROUP_OUT", AFLOWRC_DEFAULT_AFLOW_SGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_SGROUP_JSON", AFLOWRC_DEFAULT_AFLOW_SGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_AGROUP_OUT", AFLOWRC_DEFAULT_AFLOW_AGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_AGROUP_JSON", AFLOWRC_DEFAULT_AFLOW_AGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_IATOMS_OUT", AFLOWRC_DEFAULT_AFLOW_IATOMS_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_IATOMS_JSON", AFLOWRC_DEFAULT_AFLOW_IATOMS_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_ICAGES_OUT", AFLOWRC_DEFAULT_AFLOW_ICAGES_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_SURFACE_OUT", AFLOWRC_DEFAULT_AFLOW_SURFACE_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_QMVASP_OUT", AFLOWRC_DEFAULT_AFLOW_QMVASP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_ERVASP_OUT", AFLOWRC_DEFAULT_AFLOW_ERVASP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_IMMISCIBILITY_OUT", AFLOWRC_DEFAULT_AFLOW_IMMISCIBILITY_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_MEMORY_OUT", AFLOWRC_DEFAULT_AFLOW_MEMORY_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_INPUT_OUT", AFLOWRC_DEFAULT_AFLOW_FROZSL_INPUT_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_POSCAR_OUT", AFLOWRC_DEFAULT_AFLOW_FROZSL_POSCAR_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_MODES_OUT", AFLOWRC_DEFAULT_AFLOW_FROZSL_MODES_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_EIGEN_OUT", AFLOWRC_DEFAULT_AFLOW_FROZSL_EIGEN_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_END_OUT", AFLOWRC_DEFAULT_AFLOW_END_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_DIELECTRIC_FILE", AFLOWRC_DEFAULT_AFLOW_DIELECTRIC_FILE);

    // DEFAULT GENERIC MPI
    aflowrc::load_default("MPI_START_DEFAULT", AFLOWRC_MPI_START_DEFAULT);
    aflowrc::load_default("MPI_STOP_DEFAULT", AFLOWRC_MPI_STOP_DEFAULT);
    aflowrc::load_default("MPI_COMMAND_DEFAULT", AFLOWRC_MPI_COMMAND_DEFAULT);
    aflowrc::load_default("MPI_NCPUS_DEFAULT", AFLOWRC_MPI_NCPUS_DEFAULT);
    aflowrc::load_default("MPI_NCPUS_MAX", AFLOWRC_MPI_NCPUS_MAX);

    // BINARY VASP
    aflowrc::load_default("DEFAULT_VASP_GAMMA_BIN", AFLOWRC_DEFAULT_VASP_GAMMA_BIN);
    aflowrc::load_default("DEFAULT_VASP_GAMMA_MPI_BIN", AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN);
    aflowrc::load_default("DEFAULT_VASP_BIN", AFLOWRC_DEFAULT_VASP_BIN);
    aflowrc::load_default("DEFAULT_VASP_MPI_BIN", AFLOWRC_DEFAULT_VASP_MPI_BIN);
    aflowrc::load_default("DEFAULT_VASP5_BIN", AFLOWRC_DEFAULT_VASP5_BIN);
    aflowrc::load_default("DEFAULT_VASP5_MPI_BIN", AFLOWRC_DEFAULT_VASP5_MPI_BIN);
    // BINARY AIMS
    aflowrc::load_default("DEFAULT_AIMS_BIN", AFLOWRC_DEFAULT_AIMS_BIN);

    // POTCARS
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIRECTORIES", AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE", AFLOWRC_DEFAULT_VASP_POTCAR_DATE);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_SUFFIX", AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE_POT_LDA", AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE_POT_GGA", AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POT_LDA", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_LDA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POT_GGA", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_GGA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POT_PBE", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_PBE);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN", AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN);

    // DEFAULT KPOINTS/DOS
    aflowrc::load_default("DEFAULT_BANDS_GRID", AFLOWRC_DEFAULT_BANDS_GRID);
    aflowrc::load_default("DEFAULT_BANDS_LATTICE", AFLOWRC_DEFAULT_BANDS_LATTICE);
    aflowrc::load_default("DEFAULT_KSCHEME", AFLOWRC_DEFAULT_KSCHEME);
    aflowrc::load_default("DEFAULT_KPPRA", AFLOWRC_DEFAULT_KPPRA);
    aflowrc::load_default("DEFAULT_STATIC_KSCHEME", AFLOWRC_DEFAULT_STATIC_KSCHEME);
    aflowrc::load_default("DEFAULT_KPPRA_STATIC", AFLOWRC_DEFAULT_KPPRA_STATIC);
    aflowrc::load_default("DEFAULT_KPPRA_ICSD", AFLOWRC_DEFAULT_KPPRA_ICSD);
    aflowrc::load_default("DEFAULT_UNARY_BANDS_GRID", AFLOWRC_DEFAULT_UNARY_BANDS_GRID);
    aflowrc::load_default("DEFAULT_UNARY_KPPRA", AFLOWRC_DEFAULT_UNARY_KPPRA);
    aflowrc::load_default("DEFAULT_UNARY_KPPRA_STATIC", AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC);
    aflowrc::load_default("DEFAULT_UNARY_KPPRA_DIELECTRIC", AFLOWRC_DEFAULT_UNARY_KPPRA_DIELECTRIC);
    aflowrc::load_default("DEFAULT_PHONONS_KSCHEME", AFLOWRC_DEFAULT_PHONONS_KSCHEME);
    aflowrc::load_default("DEFAULT_PHONONS_KPPRA", AFLOWRC_DEFAULT_PHONONS_KPPRA);
    aflowrc::load_default("DEFAULT_DIELECTRIC_KSCHEME", AFLOWRC_DEFAULT_DIELECTRIC_KSCHEME);
    aflowrc::load_default("DEFAULT_KPPRA_DIELECTRIC", AFLOWRC_DEFAULT_KPPRA_DIELECTRIC);
    aflowrc::load_default("DEFAULT_DOS_EMIN", AFLOWRC_DEFAULT_DOS_EMIN);
    aflowrc::load_default("DEFAULT_DOS_EMAX", AFLOWRC_DEFAULT_DOS_EMAX);
    aflowrc::load_default("DEFAULT_DOS_SCALE", AFLOWRC_DEFAULT_DOS_SCALE);

    // PRECISION
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_LOW", AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_MEDIUM", AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_NORMAL", AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_HIGH", AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_ACCURATE", AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE);
    aflowrc::load_default("DEFAULT_VASP_ENMAX_MINIMUM", AFLOWRC_DEFAULT_VASP_ENMAX_MINIMUM);
    aflowrc::load_default("DEFAULT_VASP_SPIN_REMOVE_CUTOFF", AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF);
    aflowrc::load_default("DEFAULT_VASP_PREC_POTIM", AFLOWRC_DEFAULT_VASP_PREC_POTIM);
    aflowrc::load_default("DEFAULT_VASP_PREC_EDIFFG", AFLOWRC_DEFAULT_VASP_PREC_EDIFFG);

    // OPTIONS
    aflowrc::load_default("DEFAULT_VASP_OUT", AFLOWRC_DEFAULT_VASP_OUT);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_INCAR", AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_POSCAR", AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_POTCAR", AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_KPOINTS", AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS);
    aflowrc::load_default("DEFAULT_AIMS_EXTERNAL_CONTROL", AFLOWRC_DEFAULT_AIMS_EXTERNAL_CONTROL);
    aflowrc::load_default("DEFAULT_AIMS_EXTERNAL_GEOM", AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM);
    aflowrc::load_default("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE", AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SIGMA", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_NELM", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM); // CO20200624
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_NELM_STATIC", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM_STATIC); // CO20200624
    aflowrc::load_default("MAX_VASP_NELM", AFLOWRC_MAX_VASP_NELM); // CO20200624
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SYM", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_BADER", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_BADER_STATIC", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER_STATIC);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ELF", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_WAVECAR", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_CHGCAR", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_LSCOUPLING", AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);

    // AFLOW_LIBRARY AFLOW_PROJECT
    aflowrc::load_default("DEFAULT_AFLOW_LIBRARY_DIRECTORIES", AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES);
    aflowrc::load_default("DEFAULT_AFLOW_PROJECTS_DIRECTORIES", AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES);
    aflowrc::load_default("DEFAULT_AFLOWDATA_WEB_DIRECTORY", AFLOWRC_DEFAULT_AFLOWDATA_WEB_DIRECTORY); // CO+ME20200731

    // DEFAULT PLATON/FINDSYM
    aflowrc::load_default("DEFAULT_PLATON_P_EQUAL", AFLOWRC_DEFAULT_PLATON_P_EQUAL);
    aflowrc::load_default("DEFAULT_PLATON_P_EXACT", AFLOWRC_DEFAULT_PLATON_P_EXACT);
    aflowrc::load_default("DEFAULT_PLATON_P_ANG", AFLOWRC_DEFAULT_PLATON_P_ANG);
    aflowrc::load_default("DEFAULT_PLATON_P_D1", AFLOWRC_DEFAULT_PLATON_P_D1);
    aflowrc::load_default("DEFAULT_PLATON_P_D2", AFLOWRC_DEFAULT_PLATON_P_D2);
    aflowrc::load_default("DEFAULT_PLATON_P_D3", AFLOWRC_DEFAULT_PLATON_P_D3);
    aflowrc::load_default("DEFAULT_FINDSYM_TOL", AFLOWRC_DEFAULT_FINDSYM_TOL);

    // DEFAULT GNUPLOT
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT", AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT);
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_BOLD", AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD);
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_ITALICS", AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS);
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS", AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS);
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT", AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT);
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_BOLD", AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD);
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_ITALICS", AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS);
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS", AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS);
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT", AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT);
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT_BOLD", AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD);
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT_ITALICS", AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_ITALICS);
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS", AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS);

    // DEFAULT NHULL
    aflowrc::load_default("DEFAULT_NHULL_ALLOWED_DFT_TYPES", AFLOWRC_DEFAULT_NHULL_ALLOWED_DFT_TYPES);
    aflowrc::load_default("DEFAULT_NHULL_ALLOW_ALL_FORMATION_ENERGIES", AFLOWRC_DEFAULT_NHULL_ALLOW_ALL_FORMATION_ENERGIES);
    aflowrc::load_default("DEFAULT_NHULL_COUNT_THRESHOLD_BINARIES", AFLOWRC_DEFAULT_NHULL_COUNT_THRESHOLD_BINARIES);
    aflowrc::load_default("DEFAULT_NHULL_PERFORM_OUTLIER_ANALYSIS", AFLOWRC_DEFAULT_NHULL_PERFORM_OUTLIER_ANALYSIS);
    aflowrc::load_default("DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES", AFLOWRC_DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES);
    aflowrc::load_default("DEFAULT_NHULL_OUTLIER_MULTIPLIER", AFLOWRC_DEFAULT_NHULL_OUTLIER_MULTIPLIER);
    aflowrc::load_default("DEFAULT_NHULL_IGNORE_KNOWN_ILL_CONVERGED", AFLOWRC_DEFAULT_NHULL_IGNORE_KNOWN_ILL_CONVERGED);

    // DEFAULT GFA  //CO20190628
    aflowrc::load_default("DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF", AFLOWRC_DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF); // CO20190628

    // DEFAULT ARUN
    aflowrc::load_default("ARUN_DIRECTORY_PREFIX", AFLOWRC_ARUN_DIRECTORY_PREFIX);

    // DEFAULT POCC
    aflowrc::load_default("DEFAULT_POCC_STRUCTURE_GENERATION_ALGO", AFLOWRC_DEFAULT_POCC_STRUCTURE_GENERATION_ALGO);
    aflowrc::load_default("DEFAULT_POCC_TEMPERATURE_STRING", AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING);
    aflowrc::load_default("DEFAULT_POCC_EXCLUDE_UNSTABLE", AFLOWRC_DEFAULT_POCC_EXCLUDE_UNSTABLE); // ME20210927
    aflowrc::load_default("DEFAULT_POCC_SITE_TOL", AFLOWRC_DEFAULT_POCC_SITE_TOL);
    aflowrc::load_default("DEFAULT_POCC_STOICH_TOL", AFLOWRC_DEFAULT_POCC_STOICH_TOL);
    aflowrc::load_default("DEFAULT_UFF_BONDING_DISTANCE", AFLOWRC_DEFAULT_UFF_BONDING_DISTANCE);
    aflowrc::load_default("DEFAULT_UFF_ENERGY_TOLERANCE", AFLOWRC_DEFAULT_UFF_ENERGY_TOLERANCE);
    aflowrc::load_default("DEFAULT_UFF_CLUSTER_RADIUS", AFLOWRC_DEFAULT_UFF_CLUSTER_RADIUS);
    aflowrc::load_default("DEFAULT_POCC_RDF_RMAX", AFLOWRC_DEFAULT_POCC_RDF_RMAX);
    aflowrc::load_default("DEFAULT_POCC_RDF_NBINS", AFLOWRC_DEFAULT_POCC_RDF_NBINS);
    aflowrc::load_default("DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON", AFLOWRC_DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON);
    aflowrc::load_default("DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS", AFLOWRC_DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS);
    aflowrc::load_default("POCC_FILE_PREFIX", AFLOWRC_POCC_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_POCC_JSON", AFLOWRC_DEFAULT_POCC_JSON);
    aflowrc::load_default("POCC_APL_OUT_FILE", AFLOWRC_POCC_APL_OUT_FILE); // ME20210927
    aflowrc::load_default("POCC_ALL_SUPERCELLS_FILE", AFLOWRC_POCC_ALL_SUPERCELLS_FILE);
    aflowrc::load_default("POCC_UNIQUE_SUPERCELLS_FILE", AFLOWRC_POCC_UNIQUE_SUPERCELLS_FILE);
    aflowrc::load_default("POCC_ALL_HNF_MATRICES_FILE", AFLOWRC_POCC_ALL_HNF_MATRICES_FILE);
    aflowrc::load_default("POCC_ALL_SITE_CONFIGURATIONS_FILE", AFLOWRC_POCC_ALL_SITE_CONFIGURATIONS_FILE);
    aflowrc::load_default("POCC_DOSCAR_FILE", AFLOWRC_POCC_DOSCAR_FILE);
    aflowrc::load_default("POCC_PHDOSCAR_FILE", AFLOWRC_POCC_PHDOSCAR_FILE); // ME20210927
    aflowrc::load_default("POCC_ANIONS_LIST", AFLOWRC_POCC_ANIONS_LIST);

    // DEFAULT APL
    //// DEFAULT APL SUPERCELL
    aflowrc::load_default("DEFAULT_APL_PREC", AFLOWRC_DEFAULT_APL_PREC);
    aflowrc::load_default("DEFAULT_APL_ENGINE", AFLOWRC_DEFAULT_APL_ENGINE);
    aflowrc::load_default("DEFAULT_APL_HIBERNATE", AFLOWRC_DEFAULT_APL_HIBERNATE);
    aflowrc::load_default("DEFAULT_APL_MINSHELL", AFLOWRC_DEFAULT_APL_MINSHELL);
    aflowrc::load_default("DEFAULT_APL_MINATOMS", AFLOWRC_DEFAULT_APL_MINATOMS);
    aflowrc::load_default("DEFAULT_APL_POLAR", AFLOWRC_DEFAULT_APL_POLAR);
    aflowrc::load_default("DEFAULT_APL_DMAG", AFLOWRC_DEFAULT_APL_DMAG);
    aflowrc::load_default("DEFAULT_APL_DXYZONLY", AFLOWRC_DEFAULT_APL_DXYZONLY);
    aflowrc::load_default("DEFAULT_APL_DSYMMETRIZE", AFLOWRC_DEFAULT_APL_DSYMMETRIZE); // CO20181226
    aflowrc::load_default("DEFAULT_APL_DINEQUIV_ONLY", AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY); // CO20181226
    aflowrc::load_default("DEFAULT_APL_DPM", AFLOWRC_DEFAULT_APL_DPM);
    aflowrc::load_default("DEFAULT_APL_RELAX", AFLOWRC_DEFAULT_APL_RELAX);
    aflowrc::load_default("DEFAULT_APL_RELAX_COMMENSURATE", AFLOWRC_DEFAULT_APL_RELAX_COMMENSURATE); // ME20200427
    aflowrc::load_default("DEFAULT_APL_ZEROSTATE", AFLOWRC_DEFAULT_APL_ZEROSTATE);
    aflowrc::load_default("DEFAULT_APL_ZEROSTATE_CHGCAR", AFLOWRC_DEFAULT_APL_ZEROSTATE_CHGCAR); // ME20191029
    aflowrc::load_default("DEFAULT_APL_USE_LEPSILON", AFLOWRC_DEFAULT_APL_USE_LEPSILON);

    //// DEFAULT APL PHONON PROPERTIES
    aflowrc::load_default("DEFAULT_APL_FREQFORMAT", AFLOWRC_DEFAULT_APL_FREQFORMAT);
    aflowrc::load_default("DEFAULT_APL_DC", AFLOWRC_DEFAULT_APL_DC);
    aflowrc::load_default("DEFAULT_APL_DCPATH", AFLOWRC_DEFAULT_APL_DCPATH);
    aflowrc::load_default("DEFAULT_APL_DCPOINTS", AFLOWRC_DEFAULT_APL_DCPOINTS);
    aflowrc::load_default("DEFAULT_APL_DOS", AFLOWRC_DEFAULT_APL_DOS);
    aflowrc::load_default("DEFAULT_APL_DOSMETHOD", AFLOWRC_DEFAULT_APL_DOSMETHOD);
    aflowrc::load_default("DEFAULT_APL_DOSMESH", AFLOWRC_DEFAULT_APL_DOSMESH);
    aflowrc::load_default("DEFAULT_APL_DOSPOINTS", AFLOWRC_DEFAULT_APL_DOSPOINTS);
    aflowrc::load_default("DEFAULT_APL_DOSSMEAR", AFLOWRC_DEFAULT_APL_DOSSMEAR);
    aflowrc::load_default("DEFAULT_APL_DOS_PROJECT", AFLOWRC_DEFAULT_APL_DOS_PROJECT); // ME20200213
    aflowrc::load_default("DEFAULT_APL_TP", AFLOWRC_DEFAULT_APL_TP);
    aflowrc::load_default("DEFAULT_APL_DISPLACEMENTS", AFLOWRC_DEFAULT_APL_DISPLACEMENTS); // ME20200421
    aflowrc::load_default("DEFAULT_APL_TPT", AFLOWRC_DEFAULT_APL_TPT);
    aflowrc::load_default("DEFAULT_APL_GVEL", AFLOWRC_DEFAULT_APL_GVEL); // ME20200517

    //// DEFAULT APL FILES
    aflowrc::load_default("DEFAULT_APL_FILE_PREFIX", AFLOWRC_DEFAULT_APL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_APL_OUT_FILE", AFLOWRC_DEFAULT_APL_OUT_FILE); // ME20210927
    aflowrc::load_default("DEFAULT_APL_PDIS_FILE", AFLOWRC_DEFAULT_APL_PDIS_FILE);
    aflowrc::load_default("DEFAULT_APL_PDOS_FILE", AFLOWRC_DEFAULT_APL_PDOS_FILE);
    aflowrc::load_default("DEFAULT_APL_THERMO_FILE", AFLOWRC_DEFAULT_APL_THERMO_FILE);
    aflowrc::load_default("DEFAULT_APL_THERMO_JSON", AFLOWRC_DEFAULT_APL_THERMO_JSON); // ME20211019
    aflowrc::load_default("DEFAULT_APL_DYNMAT_FILE", AFLOWRC_DEFAULT_APL_DYNMAT_FILE);
    aflowrc::load_default("DEFAULT_APL_HARMIFC_FILE", AFLOWRC_DEFAULT_APL_HARMIFC_FILE);
    aflowrc::load_default("DEFAULT_APL_POLAR_FILE", AFLOWRC_DEFAULT_APL_POLAR_FILE); // ME20200415
    aflowrc::load_default("DEFAULT_APL_HSKPTS_FILE", AFLOWRC_DEFAULT_APL_HSKPTS_FILE);
    aflowrc::load_default("DEFAULT_APL_MSQRDISP_FILE", AFLOWRC_DEFAULT_APL_MSQRDISP_FILE); // ME20200329
    aflowrc::load_default("DEFAULT_APL_GVEL_FILE", AFLOWRC_DEFAULT_APL_GVEL_FILE); // ME20200517
    // ME20190614 BEGIN
    aflowrc::load_default("DEFAULT_APL_PHDOSCAR_FILE", AFLOWRC_DEFAULT_APL_PHDOSCAR_FILE);
    aflowrc::load_default("DEFAULT_APL_PHPOSCAR_FILE", AFLOWRC_DEFAULT_APL_PHPOSCAR_FILE);
    aflowrc::load_default("DEFAULT_APL_PHKPOINTS_FILE", AFLOWRC_DEFAULT_APL_PHKPOINTS_FILE);
    aflowrc::load_default("DEFAULT_APL_PHEIGENVAL_FILE", AFLOWRC_DEFAULT_APL_PHEIGENVAL_FILE);
    // ME20190614 END
    aflowrc::load_default("DEFAULT_APL_STATE_FILE", AFLOWRC_DEFAULT_APL_STATE_FILE); // ME20200224
    // ME20200329 BEGIN
    aflowrc::load_default("DEFAULT_APL_ADISP_SCENE_FORMAT", AFLOWRC_DEFAULT_APL_ADISP_SCENE_FORMAT);
    aflowrc::load_default("DEFAULT_APL_ADISP_AMPLITUDE", AFLOWRC_DEFAULT_APL_ADISP_AMPLITUDE);
    aflowrc::load_default("DEFAULT_APL_ADISP_NSTEPS", AFLOWRC_DEFAULT_APL_ADISP_NSTEPS);
    aflowrc::load_default("DEFAULT_APL_ADISP_NPERIODS", AFLOWRC_DEFAULT_APL_ADISP_NPERIODS);
    // ME20200329 END

    // DEFAULT QHA
    //// DEFAULT QHA VALUES
    aflowrc::load_default("DEFAULT_QHA_MODE", AFLOWRC_DEFAULT_QHA_MODE);
    aflowrc::load_default("DEFAULT_QHA_EOS", AFLOWRC_DEFAULT_QHA_EOS);
    aflowrc::load_default("DEFAULT_QHA_EOS_DISTORTION_RANGE", AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE);
    aflowrc::load_default("DEFAULT_QHA_EOS_MODEL", AFLOWRC_DEFAULT_QHA_EOS_MODEL); // AS20200818
    aflowrc::load_default("DEFAULT_QHA_GP_DISTORTION", AFLOWRC_DEFAULT_QHA_GP_DISTORTION);
    aflowrc::load_default("DEFAULT_QHA_TAYLOR_EXPANSION_ORDER", AFLOWRC_DEFAULT_QHA_TAYLOR_EXPANSION_ORDER); // AS20200602
    aflowrc::load_default("DEFAULT_QHA_INCLUDE_ELEC_CONTRIB", AFLOWRC_DEFAULT_QHA_INCLUDE_ELEC_CONTRIB);
    aflowrc::load_default("DEFAULT_QHA_SOMMERFELD_EXPANSION", AFLOWRC_DEFAULT_QHA_SOMMERFELD_EXPANSION); // AS20200528
    aflowrc::load_default("DEFAULT_QHA_PDIS_T", AFLOWRC_DEFAULT_QHA_PDIS_T);
    // AS20200508 BEGIN
    aflowrc::load_default("DEFAULT_QHA_GP_FINITE_DIFF", AFLOWRC_DEFAULT_QHA_GP_FINITE_DIFF);
    aflowrc::load_default("DEFAULT_QHA_IGNORE_IMAGINARY", AFLOWRC_DEFAULT_QHA_IGNORE_IMAGINARY);
    aflowrc::load_default("DEFAULT_QHA_RELAX_IONS_CELL", AFLOWRC_DEFAULT_QHA_RELAX_IONS_CELL); // AS20201123
    //// DEFAULT QHA FILES
    aflowrc::load_default("DEFAULT_QHA_FILE_PREFIX", AFLOWRC_DEFAULT_QHA_FILE_PREFIX);
    // AS20200709 BEGIN
    aflowrc::load_default("DEFAULT_QHA3P_FILE_PREFIX", AFLOWRC_DEFAULT_QHA3P_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_QHANP_FILE_PREFIX", AFLOWRC_DEFAULT_QHANP_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_SCQHA_FILE_PREFIX", AFLOWRC_DEFAULT_SCQHA_FILE_PREFIX);
    // AS20200709 END
    aflowrc::load_default("DEFAULT_QHA_GP_PATH_FILE", AFLOWRC_DEFAULT_QHA_GP_PATH_FILE);
    aflowrc::load_default("DEFAULT_QHA_GP_MESH_FILE", AFLOWRC_DEFAULT_QHA_GP_MESH_FILE);
    aflowrc::load_default("DEFAULT_QHA_GP_AVG_FILE", AFLOWRC_DEFAULT_QHA_GP_AVG_FILE);
    aflowrc::load_default("DEFAULT_QHA_THERMO_FILE", AFLOWRC_DEFAULT_QHA_THERMO_FILE);
    aflowrc::load_default("DEFAULT_QHA_FREQS_FILE", AFLOWRC_DEFAULT_QHA_FREQS_FILE);
    aflowrc::load_default("DEFAULT_QHA_FVT_FILE", AFLOWRC_DEFAULT_QHA_FVT_FILE);
    // AS20200508 END
    aflowrc::load_default("DEFAULT_QHA_COEFF_FILE", AFLOWRC_DEFAULT_QHA_COEFF_FILE); // AS20210517
    aflowrc::load_default("DEFAULT_QHA_IMAG_FILE", AFLOWRC_DEFAULT_QHA_IMAG_FILE); // AS20210517
    aflowrc::load_default("DEFAULT_QHA_PDIS_FILE", AFLOWRC_DEFAULT_QHA_PDIS_FILE); // AS20201022
    aflowrc::load_default("DEFAULT_QHA_PDOS_FILE", AFLOWRC_DEFAULT_QHA_PDOS_FILE); // AS20201201
    aflowrc::load_default("DEFAULT_QHA_KPOINTS_FILE", AFLOWRC_DEFAULT_QHA_KPOINTS_FILE); // AS20201112
    // AS20210914 BEGIN
    aflowrc::load_default("DEFAULT_POCC_QHA_THERMO_FILE", AFLOWRC_DEFAULT_POCC_QHA_THERMO_FILE);
    aflowrc::load_default("DEFAULT_POCC_QHA_AVGTHERMO_FILE", AFLOWRC_DEFAULT_POCC_QHA_AVGTHERMO_FILE);
    // AS20210914 END

    // DEFAULT AAPL
    //// DEFAULT AAPL VALUES
    aflowrc::load_default("DEFAULT_AAPL_BTE", AFLOWRC_DEFAULT_AAPL_BTE);
    //[ME20181226]aflowrc::load_default("DEFAULT_AAPL_BZMETHOD",AFLOWRC_DEFAULT_AAPL_BZMETHOD);
    aflowrc::load_default("DEFAULT_AAPL_FOURTH_ORDER", AFLOWRC_DEFAULT_AAPL_FOURTH_ORDER);
    aflowrc::load_default("DEFAULT_AAPL_CUT_RAD", AFLOWRC_DEFAULT_AAPL_CUT_RAD);
    aflowrc::load_default("DEFAULT_AAPL_CUT_SHELL", AFLOWRC_DEFAULT_AAPL_CUT_SHELL);
    aflowrc::load_default("DEFAULT_AAPL_THERMALGRID", AFLOWRC_DEFAULT_AAPL_THERMALGRID);
    aflowrc::load_default("DEFAULT_AAPL_TCT", AFLOWRC_DEFAULT_AAPL_TCT);
    aflowrc::load_default("DEFAULT_AAPL_SUMRULE", AFLOWRC_DEFAULT_AAPL_SUMRULE);
    aflowrc::load_default("DEFAULT_AAPL_SUMRULE_MAX_ITER", AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER);
    aflowrc::load_default("DEFAULT_AAPL_MIXING_COEFFICIENT", AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT);
    aflowrc::load_default("DEFAULT_AAPL_ISOTOPE", AFLOWRC_DEFAULT_AAPL_ISOTOPE);
    aflowrc::load_default("DEFAULT_AAPL_BOUNDARY", AFLOWRC_DEFAULT_AAPL_BOUNDARY);
    aflowrc::load_default("DEFAULT_AAPL_CUMULATIVEK", AFLOWRC_DEFAULT_AAPL_CUMULATIVEK);
    aflowrc::load_default("DEFAULT_AAPL_NANO_SIZE", AFLOWRC_DEFAULT_AAPL_NANO_SIZE);
    //// DEFAULT AAPL FILES
    aflowrc::load_default("DEFAULT_AAPL_FILE_PREFIX", AFLOWRC_DEFAULT_AAPL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_AAPL_IRRQPTS_FILE", AFLOWRC_DEFAULT_AAPL_IRRQPTS_FILE);
    aflowrc::load_default("DEFAULT_AAPL_GVEL_FILE", AFLOWRC_DEFAULT_AAPL_GVEL_FILE);
    aflowrc::load_default("DEFAULT_AAPL_PS_FILE", AFLOWRC_DEFAULT_AAPL_PS_FILE); // ME20191104
    aflowrc::load_default("DEFAULT_AAPL_GRUENEISEN_FILE", AFLOWRC_DEFAULT_AAPL_GRUENEISEN_FILE); // ME20191104
    aflowrc::load_default("DEFAULT_AAPL_RATES_FILE", AFLOWRC_DEFAULT_AAPL_RATES_FILE);
    aflowrc::load_default("DEFAULT_AAPL_RATES_3RD_FILE", AFLOWRC_DEFAULT_AAPL_RATES_3RD_FILE);
    aflowrc::load_default("DEFAULT_AAPL_RATES_4TH_FILE", AFLOWRC_DEFAULT_AAPL_RATES_4TH_FILE);
    aflowrc::load_default("DEFAULT_AAPL_ISOTOPE_FILE", AFLOWRC_DEFAULT_AAPL_ISOTOPE_FILE);
    aflowrc::load_default("DEFAULT_AAPL_BOUNDARY_FILE", AFLOWRC_DEFAULT_AAPL_BOUNDARY_FILE);
    aflowrc::load_default("DEFAULT_AAPL_TCOND_FILE", AFLOWRC_DEFAULT_AAPL_TCOND_FILE);

    // DEFAULT AEL
    //// DEFAULT AEL STRAIN CALCS
    aflowrc::load_default("DEFAULT_AEL_STRAIN_SYMMETRY", AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY);
    aflowrc::load_default("DEFAULT_AEL_NNORMAL_STRAINS", AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS);
    aflowrc::load_default("DEFAULT_AEL_NSHEAR_STRAINS", AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS);
    aflowrc::load_default("DEFAULT_AEL_NORMAL_STRAIN_STEP", AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP);
    aflowrc::load_default("DEFAULT_AEL_SHEAR_STRAIN_STEP", AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP);
    aflowrc::load_default("DEFAULT_AEL_ORIGIN_STRAIN_CALC", AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC);
    aflowrc::load_default("DEFAULT_AEL_ORIGIN_STRAIN_FIT", AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT);
    aflowrc::load_default("DEFAULT_AEL_RELAXED_STRUCT_FIT", AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT);
    aflowrc::load_default("DEFAULT_AEL_NEG_STRAINS", AFLOWRC_DEFAULT_AEL_NEG_STRAINS);
    aflowrc::load_default("DEFAULT_AEL_NIND_STRAIN_DIRS", AFLOWRC_DEFAULT_AEL_VASPSYM);
    aflowrc::load_default("DEFAULT_AEL_VASPSYM", AFLOWRC_DEFAULT_AEL_VASPSYM);
    aflowrc::load_default("DEFAULT_AEL_PRECACC_ALGONORM", AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM);
    aflowrc::load_default("DEFAULT_AEL_VASPRUNXML_STRESS", AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS);
    aflowrc::load_default("DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS", AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS);
    aflowrc::load_default("DEFAULT_AEL_SKIP_ARUNS_MAX", AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX);

    //// DEFAULT AEL CHECKS AND PROCESSING
    aflowrc::load_default("DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY", AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY);
    aflowrc::load_default("DEFAULT_AEL_SYMMETRIZE", AFLOWRC_DEFAULT_AEL_SYMMETRIZE);

    //// DEFAULT AEL OUTPUT FILES
    aflowrc::load_default("DEFAULT_AEL_FILE_PREFIX", AFLOWRC_DEFAULT_AEL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_AEL_WRITE_FULL_RESULTS", AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS);
    aflowrc::load_default("DEFAULT_AEL_DIRNAME_ARUN", AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN);

    // DEFAULT AGL
    //// DEFAULT AGL STRAIN CALCS
    aflowrc::load_default("DEFAULT_AGL_AEL_POISSON_RATIO", AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO);
    aflowrc::load_default("DEFAULT_AGL_NSTRUCTURES", AFLOWRC_DEFAULT_AGL_NSTRUCTURES);
    aflowrc::load_default("DEFAULT_AGL_STRAIN_STEP", AFLOWRC_DEFAULT_AGL_STRAIN_STEP);
    aflowrc::load_default("DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS", AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS);
    aflowrc::load_default("DEFAULT_AGL_SKIP_ARUNS_MAX", AFLOWRC_DEFAULT_AGL_SKIP_ARUNS_MAX);

    //// DEFAULT AGL CHECKS AND PROCESSING
    aflowrc::load_default("DEFAULT_AGL_NTEMPERATURE", AFLOWRC_DEFAULT_AGL_NTEMPERATURE);
    aflowrc::load_default("DEFAULT_AGL_STEMPERATURE", AFLOWRC_DEFAULT_AGL_STEMPERATURE);
    aflowrc::load_default("DEFAULT_AGL_NPRESSURE", AFLOWRC_DEFAULT_AGL_NPRESSURE);
    aflowrc::load_default("DEFAULT_AGL_SPRESSURE", AFLOWRC_DEFAULT_AGL_SPRESSURE);
    aflowrc::load_default("DEFAULT_AGL_POISSON_RATIO", AFLOWRC_DEFAULT_AGL_POISSON_RATIO);
    aflowrc::load_default("DEFAULT_AGL_IEOS", AFLOWRC_DEFAULT_AGL_IEOS);
    aflowrc::load_default("DEFAULT_AGL_IDEBYE", AFLOWRC_DEFAULT_AGL_IDEBYE);
    aflowrc::load_default("DEFAULT_AGL_FIT_TYPE", AFLOWRC_DEFAULT_AGL_FIT_TYPE);
    aflowrc::load_default("DEFAULT_AGL_CHECK_EV_CONCAVITY", AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY);
    aflowrc::load_default("DEFAULT_AGL_CHECK_EV_MIN", AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN);
    aflowrc::load_default("DEFAULT_AGL_HUGONIOT_CALC", AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC);
    aflowrc::load_default("DEFAULT_AGL_HUGONIOT_EXTRAPOLATE", AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE);
    aflowrc::load_default("DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE", AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE);

    //// DEFAULT AGL OUTPUT FILES
    aflowrc::load_default("DEFAULT_AGL_FILE_PREFIX", AFLOWRC_DEFAULT_AGL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_AGL_WRITE_FULL_RESULTS", AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS);
    aflowrc::load_default("DEFAULT_AGL_DIRNAME_ARUN", AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN);
    aflowrc::load_default("DEFAULT_AGL_WRITE_GIBBS_INPUT", AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT);
    aflowrc::load_default("DEFAULT_AGL_PLOT_RESULTS", AFLOWRC_DEFAULT_AGL_PLOT_RESULTS);

    // DEFAULT QCA
    aflowrc::load_default("DEFAULT_QCA_MIN_SLEEP_SECONDS", AFLOWRC_DEFAULT_QCA_MIN_SLEEP_SECONDS);
    aflowrc::load_default("DEFAULT_QCA_MAX_NUM_ATOMS", AFLOWRC_DEFAULT_QCA_MAX_NUM_ATOMS);
    aflowrc::load_default("DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS", AFLOWRC_DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS);
    aflowrc::load_default("DEFAULT_QCA_CV_CUTOFF", AFLOWRC_DEFAULT_QCA_CV_CUTOFF);
    aflowrc::load_default("DEFAULT_QCA_CONC_NPTS", AFLOWRC_DEFAULT_QCA_CONC_NPTS);
    aflowrc::load_default("DEFAULT_QCA_TEMP_NPTS", AFLOWRC_DEFAULT_QCA_TEMP_NPTS);
    aflowrc::load_default("DEFAULT_QCA_TEMP_MIN", AFLOWRC_DEFAULT_QCA_TEMP_MIN);
    aflowrc::load_default("DEFAULT_QCA_TEMP_MAX", AFLOWRC_DEFAULT_QCA_TEMP_MAX);
    aflowrc::load_default("DEFAULT_QCA_TEMP_MIN_LIMIT", AFLOWRC_DEFAULT_QCA_TEMP_MIN_LIMIT);
    aflowrc::load_default("DEFAULT_QCA_PRINT", AFLOWRC_DEFAULT_QCA_PRINT);

    // RF20200413 START
    //  DEFAULT CCE
    aflowrc::load_default("DEFAULT_CCE_OX_METHOD", AFLOWRC_DEFAULT_CCE_OX_METHOD);
    aflowrc::load_default("DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION", AFLOWRC_DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION);
    aflowrc::load_default("DEFAULT_CCE_OX_TOL", AFLOWRC_DEFAULT_CCE_OX_TOL);
    aflowrc::load_default("DEFAULT_CCE_PEROX_CUTOFF", AFLOWRC_DEFAULT_CCE_PEROX_CUTOFF);
    aflowrc::load_default("DEFAULT_CCE_SUPEROX_CUTOFF", AFLOWRC_DEFAULT_CCE_SUPEROX_CUTOFF);
    aflowrc::load_default("DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF", AFLOWRC_DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF);
    aflowrc::load_default("DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF", AFLOWRC_DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF);
    // RF20200413 END

    // DEFAULT XTALFINDER
    aflowrc::load_default("DEFAULT_XTALFINDER_MISFIT_MATCH", AFLOWRC_DEFAULT_XTALFINDER_MISFIT_MATCH); // DX20201118
    aflowrc::load_default("DEFAULT_XTALFINDER_MISFIT_FAMILY", AFLOWRC_DEFAULT_XTALFINDER_MISFIT_FAMILY); // DX20201118
    aflowrc::load_default("DEFAULT_XTALFINDER_SUPERCELL_METHOD", AFLOWRC_DEFAULT_XTALFINDER_SUPERCELL_METHOD); // DX20201223
    aflowrc::load_default("DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING", AFLOWRC_DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING); // DX20200709
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_MATERIAL", AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL); // DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_STRUCTURE", AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE); // DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_DUPLICATE", AFLOWRC_DEFAULT_XTALFINDER_FILE_DUPLICATE); // DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE", AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE); // DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE", AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE); // DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE", AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE); // DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE", AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE); // DX20201228

    // DX20200720 - START
    //  DEFAULT ANRL
    aflowrc::load_default("DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL", AFLOWRC_DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL);
    // DX20200720 - END

    // DEFAULT CORE
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_BEEP", AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP);
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_HALT", AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT);
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_REFRESH", AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH);

    // VASP MACHINE SETTINGS
    aflowrc::load_default("SECONDS_SLEEP_VASP_COMPLETION", AFLOWRC_SECONDS_SLEEP_VASP_COMPLETION); // CO20201111
    aflowrc::load_default("SECONDS_SLEEP_VASP_MONITOR", AFLOWRC_SECONDS_SLEEP_VASP_MONITOR); // CO20201111
    aflowrc::load_default("SECONDS_STALE_OUTCAR", AFLOWRC_SECONDS_STALE_OUTCAR); // CO20201111
    aflowrc::load_default("BYTES_MAX_VASP_OUT", AFLOWRC_BYTES_MAX_VASP_OUT); // CO20201111
    aflowrc::load_default("MEMORY_MAX_USAGE_RAM", AFLOWRC_MEMORY_MAX_USAGE_RAM); // CO20201111
    aflowrc::load_default("MEMORY_MAX_USAGE_SWAP", AFLOWRC_MEMORY_MAX_USAGE_SWAP); // CO20201111
    aflowrc::load_default("FILE_VASP_MONITOR", AFLOWRC_FILE_VASP_MONITOR); // CO20201111
    aflowrc::load_default("INTEL_COMPILER_PATHS", AFLOWRC_INTEL_COMPILER_PATHS); // CO20201111

    // DEFAULT MACHINE DEPENDENT MPI
    aflowrc::load_default("MPI_OPTIONS_DUKE_BETA_MPICH", AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH);
    aflowrc::load_default("MPI_COMMAND_DUKE_BETA_MPICH", AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_BETA_MPICH", AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH);

    aflowrc::load_default("MPI_OPTIONS_DUKE_BETA_OPENMPI", AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI);
    aflowrc::load_default("MPI_COMMAND_DUKE_BETA_OPENMPI", AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_BETA_OPENMPI", AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI);

    aflowrc::load_default("MPI_OPTIONS_DUKE_MATERIALS", AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS);
    aflowrc::load_default("MPI_COMMAND_DUKE_MATERIALS", AFLOWRC_MPI_COMMAND_DUKE_MATERIALS);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_MATERIALS", AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS);

    aflowrc::load_default("MPI_OPTIONS_DUKE_AFLOWLIB", AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB);
    aflowrc::load_default("MPI_COMMAND_DUKE_AFLOWLIB", AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_AFLOWLIB", AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB);

    aflowrc::load_default("MPI_OPTIONS_DUKE_HABANA", AFLOWRC_MPI_OPTIONS_DUKE_HABANA);
    aflowrc::load_default("MPI_COMMAND_DUKE_HABANA", AFLOWRC_MPI_COMMAND_DUKE_HABANA);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_HABANA", AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA);

    aflowrc::load_default("MPI_OPTIONS_DUKE_QRATS_MPICH", AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH);
    aflowrc::load_default("MPI_COMMAND_DUKE_QRATS_MPICH", AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_QRATS_MPICH", AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH);

    aflowrc::load_default("MPI_OPTIONS_DUKE_QFLOW_OPENMPI", AFLOWRC_MPI_OPTIONS_DUKE_QFLOW_OPENMPI);
    aflowrc::load_default("MPI_COMMAND_DUKE_QFLOW_OPENMPI", AFLOWRC_MPI_COMMAND_DUKE_QFLOW_OPENMPI);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI", AFLOWRC_MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI);

    // CO20201220 X START
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_X", AFLOWRC_MPI_OPTIONS_DUKE_X_X);
    aflowrc::load_default("MPI_COMMAND_DUKE_X_X", AFLOWRC_MPI_COMMAND_DUKE_X_X);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_X", AFLOWRC_MPI_BINARY_DIR_DUKE_X_X);
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_CRAY", AFLOWRC_MPI_OPTIONS_DUKE_X_CRAY);
    aflowrc::load_default("MPI_COMMAND_DUKE_X_CRAY", AFLOWRC_MPI_COMMAND_DUKE_X_CRAY);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_CRAY", AFLOWRC_MPI_BINARY_DIR_DUKE_X_CRAY);
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_OLDCRAY", AFLOWRC_MPI_OPTIONS_DUKE_X_OLDCRAY);
    aflowrc::load_default("MPI_COMMAND_DUKE_X_OLDCRAY", AFLOWRC_MPI_COMMAND_DUKE_X_OLDCRAY);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_OLDCRAY", AFLOWRC_MPI_BINARY_DIR_DUKE_X_OLDCRAY);
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_SMB", AFLOWRC_MPI_OPTIONS_DUKE_X_SMB);
    aflowrc::load_default("MPI_COMMAND_DUKE_X_SMB", AFLOWRC_MPI_COMMAND_DUKE_X_SMB);
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_SMB", AFLOWRC_MPI_BINARY_DIR_DUKE_X_SMB);
    // CO20201220 X STOP

    // CO20220818 JHU_ROCKFISH START
    aflowrc::load_default("MPI_OPTIONS_JHU_ROCKFISH", AFLOWRC_MPI_OPTIONS_JHU_ROCKFISH);
    aflowrc::load_default("MPI_COMMAND_JHU_ROCKFISH", AFLOWRC_MPI_COMMAND_JHU_ROCKFISH);
    aflowrc::load_default("MPI_BINARY_DIR_JHU_ROCKFISH", AFLOWRC_MPI_BINARY_DIR_JHU_ROCKFISH);
    // CO20220818 JHU_ROCKFISH STOP

    // DX20190509 - MACHINE001 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE001", AFLOWRC_MPI_OPTIONS_MACHINE001);
    aflowrc::load_default("MPI_COMMAND_MACHINE001", AFLOWRC_MPI_COMMAND_MACHINE001);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE001", AFLOWRC_MPI_BINARY_DIR_MACHINE001);
    // DX20190509 - MACHINE001 - END

    // DX20190509 - MACHINE002 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE002", AFLOWRC_MPI_OPTIONS_MACHINE002);
    aflowrc::load_default("MPI_COMMAND_MACHINE002", AFLOWRC_MPI_COMMAND_MACHINE002);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE002", AFLOWRC_MPI_BINARY_DIR_MACHINE002);
    // DX20190509 - MACHINE002 - END

    // DX20201005 - MACHINE003 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE003", AFLOWRC_MPI_OPTIONS_MACHINE003);
    aflowrc::load_default("MPI_COMMAND_MACHINE003", AFLOWRC_MPI_COMMAND_MACHINE003);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE003", AFLOWRC_MPI_BINARY_DIR_MACHINE003);
    // DX20201005 - MACHINE003 - END

    // DX20211011 - MACHINE004 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE004", AFLOWRC_MPI_OPTIONS_MACHINE004);
    aflowrc::load_default("MPI_COMMAND_MACHINE004", AFLOWRC_MPI_COMMAND_MACHINE004);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE004", AFLOWRC_MPI_BINARY_DIR_MACHINE004);
    // DX20211011 - MACHINE004 - END

    // DX20190107 - CMU EULER - START
    aflowrc::load_default("MPI_OPTIONS_CMU_EULER", AFLOWRC_MPI_OPTIONS_CMU_EULER);
    aflowrc::load_default("MPI_COMMAND_CMU_EULER", AFLOWRC_MPI_COMMAND_CMU_EULER);
    aflowrc::load_default("MPI_BINARY_DIR_CMU_EULER", AFLOWRC_MPI_BINARY_DIR_CMU_EULER);
    // DX20190107 - CMU EULER - END

    aflowrc::load_default("MPI_OPTIONS_MPCDF_EOS", AFLOWRC_MPI_OPTIONS_MPCDF_EOS);
    aflowrc::load_default("MPI_COMMAND_MPCDF_EOS", AFLOWRC_MPI_COMMAND_MPCDF_EOS);
    aflowrc::load_default("MPI_NCPUS_MPCDF_EOS", AFLOWRC_MPI_NCPUS_MPCDF_EOS);
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_EOS", AFLOWRC_MPI_HYPERTHREADING_MPCDF_EOS);
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_EOS", AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS);

    aflowrc::load_default("MPI_OPTIONS_MPCDF_DRACO", AFLOWRC_MPI_OPTIONS_MPCDF_DRACO);
    aflowrc::load_default("MPI_COMMAND_MPCDF_DRACO", AFLOWRC_MPI_COMMAND_MPCDF_DRACO);
    aflowrc::load_default("MPI_NCPUS_MPCDF_DRACO", AFLOWRC_MPI_NCPUS_MPCDF_DRACO);
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_DRACO", AFLOWRC_MPI_HYPERTHREADING_MPCDF_DRACO);
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_DRACO", AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO);

    aflowrc::load_default("MPI_OPTIONS_MPCDF_COBRA", AFLOWRC_MPI_OPTIONS_MPCDF_COBRA);
    aflowrc::load_default("MPI_COMMAND_MPCDF_COBRA", AFLOWRC_MPI_COMMAND_MPCDF_COBRA);
    aflowrc::load_default("MPI_NCPUS_MPCDF_COBRA", AFLOWRC_MPI_NCPUS_MPCDF_COBRA);
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_COBRA", AFLOWRC_MPI_HYPERTHREADING_MPCDF_COBRA);
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_COBRA", AFLOWRC_MPI_BINARY_DIR_MPCDF_COBRA);

    aflowrc::load_default("MPI_OPTIONS_MPCDF_HYDRA", AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA);
    aflowrc::load_default("MPI_COMMAND_MPCDF_HYDRA", AFLOWRC_MPI_COMMAND_MPCDF_HYDRA);
    aflowrc::load_default("MPI_NCPUS_MPCDF_HYDRA", AFLOWRC_MPI_NCPUS_MPCDF_HYDRA);
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_HYDRA", AFLOWRC_MPI_HYPERTHREADING_MPCDF_HYDRA);
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_HYDRA", AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA);

    aflowrc::load_default("MPI_OPTIONS_FULTON_MARYLOU", AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU);
    aflowrc::load_default("MPI_COMMAND_FULTON_MARYLOU", AFLOWRC_MPI_COMMAND_FULTON_MARYLOU);
    aflowrc::load_default("MPI_BINARY_DIR_FULTON_MARYLOU", AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU);

    aflowrc::load_default("MPI_OPTIONS_MACHINE1", AFLOWRC_MPI_OPTIONS_MACHINE1);
    aflowrc::load_default("MPI_COMMAND_MACHINE1", AFLOWRC_MPI_COMMAND_MACHINE1);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE1", AFLOWRC_MPI_BINARY_DIR_MACHINE1);

    aflowrc::load_default("MPI_OPTIONS_MACHINE2", AFLOWRC_MPI_OPTIONS_MACHINE2);
    aflowrc::load_default("MPI_COMMAND_MACHINE2", AFLOWRC_MPI_COMMAND_MACHINE2);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE2", AFLOWRC_MPI_BINARY_DIR_MACHINE2);

    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " END" << endl;
    }

    return true;
  }
} // namespace aflowrc

// ***************************************************************************
// aflowrc::write_default
// ***************************************************************************
namespace aflowrc {
  bool write_default(std::ostream& oss, bool AFLOWRC_VERBOSE) {
    const bool LDEBUG = (false || XHOST.DEBUG || AFLOWRC_VERBOSE);
    stringstream message;
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " BEGIN" << endl;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;
    }
    if (XHOST.aflowrc_filename.empty()) {
      XHOST.aflowrc_filename = AFLOWRC_FILENAME_LOCAL;
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
    }

    stringstream aflowrc("");
    aflowrc << "// ****************************************************************************************************" << endl;
    aflowrc << "// *                                                                                                  *" << endl;
    aflowrc << "// *                          aflow - Automatic-FLOW for materials discovery                          *" << endl;
    aflowrc << "// *                aflow.org consortium - High-Throughput ab-initio Computing Project                *" << endl;
    aflowrc << "// *                                                                                                  *" << endl;
    aflowrc << "// ****************************************************************************************************" << endl;
    aflowrc << "// DEFAULT .aflow.rc generated by AFLOW V" << string(AFLOW_VERSION) << endl;
    aflowrc << "// comments with // ignored... " << endl;
    aflowrc << "// strings are with=\"...\" " << endl;

    aflowrc << " " << endl;
    aflowrc << "AFLOWRC=\"" << AFLOWRC_AFLOWRC << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULT DEFINITIONS" << endl;
    aflowrc << "DEFAULT_KZIP_BIN=\"" << AFLOWRC_DEFAULT_KZIP_BIN << "\"" << endl;
    aflowrc << "DEFAULT_KZIP_EXT=\"" << AFLOWRC_DEFAULT_KZIP_EXT << "\"" << endl;
    aflowrc << "DEFAULT_TMPFS_DIRECTORIES=\"" << AFLOWRC_DEFAULT_TMPFS_DIRECTORIES << "\"" << endl;

    aflowrc << " " << endl;

    // HE20220218 START
    aflowrc << "// DEFAULTS ENTRY LOADER" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_AFLUX_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_AFLUX_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_RESTAPI_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_RESTAPI_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_FS_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH << "\"" << endl;
    aflowrc << " " << endl;
    // HE20220218 STOP

    // ME20191001 START
    aflowrc << "// DEFAULT AFLOW DATABASE" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_STATS_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_DATA_PATH=\"" << AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_LOCK_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_STALE_THRESHOLD=" << AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD << endl;
    aflowrc << " " << endl;
    // ME20191001 STOP
    aflowrc << "// FILENAMES FOR AFLOW.ORG ANALYSIS" << endl;
    aflowrc << "DEFAULT_FILE_AFLOWLIB_ENTRY_OUT=\"" << AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_AFLOWLIB_ENTRY_JSON=\"" << AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_EDATA_ORIG_OUT=\"" << AFLOWRC_DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_EDATA_RELAX_OUT=\"" << AFLOWRC_DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_EDATA_BANDS_OUT=\"" << AFLOWRC_DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_DATA_ORIG_OUT=\"" << AFLOWRC_DEFAULT_FILE_DATA_ORIG_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_DATA_RELAX_OUT=\"" << AFLOWRC_DEFAULT_FILE_DATA_RELAX_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_DATA_BANDS_OUT=\"" << AFLOWRC_DEFAULT_FILE_DATA_BANDS_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_EDATA_ORIG_JSON=\"" << AFLOWRC_DEFAULT_FILE_EDATA_ORIG_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_EDATA_RELAX_JSON=\"" << AFLOWRC_DEFAULT_FILE_EDATA_RELAX_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_EDATA_BANDS_JSON=\"" << AFLOWRC_DEFAULT_FILE_EDATA_BANDS_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_DATA_ORIG_JSON=\"" << AFLOWRC_DEFAULT_FILE_DATA_ORIG_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_DATA_RELAX_JSON=\"" << AFLOWRC_DEFAULT_FILE_DATA_RELAX_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_DATA_BANDS_JSON=\"" << AFLOWRC_DEFAULT_FILE_DATA_BANDS_JSON << "\"" << endl;
    aflowrc << "DEFAULT_FILE_TIME_OUT=\"" << AFLOWRC_DEFAULT_FILE_TIME_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_SPACEGROUP1_OUT=\"" << AFLOWRC_DEFAULT_FILE_SPACEGROUP1_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_SPACEGROUP2_OUT=\"" << AFLOWRC_DEFAULT_FILE_SPACEGROUP2_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_VOLDISTPARAMS_OUT=\"" << AFLOWRC_DEFAULT_FILE_VOLDISTPARAMS_OUT << "\"" << endl;
    aflowrc << "DEFAULT_FILE_VOLDISTEVOLUTION_OUT=\"" << AFLOWRC_DEFAULT_FILE_VOLDISTEVOLUTION_OUT << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// FILENAMES FOR AFLOW OPERATION" << endl;
    aflowrc << "DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PRESCRIPT_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PRESCRIPT_COMMAND=\"" << AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_COMMAND << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_POSTSCRIPT_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_POSTSCRIPT_COMMAND=\"" << AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_COMMAND << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUP_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUP_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUP_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUP_XTAL_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUP_XTAL_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT << "\"" << endl; // DX20200129
    aflowrc << "DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON << "\"" << endl; // DX20200129
    aflowrc << "DEFAULT_AFLOW_PGROUPK_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUPK_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUPK_XTAL_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PGROUPK_XTAL_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_FGROUP_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_FGROUP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_FGROUP_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_FGROUP_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_SGROUP_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_SGROUP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_SGROUP_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_SGROUP_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_AGROUP_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_AGROUP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_AGROUP_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_AGROUP_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_IATOMS_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_IATOMS_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_IATOMS_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_IATOMS_JSON << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_ICAGES_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_ICAGES_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_SURFACE_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_SURFACE_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_QMVASP_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_QMVASP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_ERVASP_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_ERVASP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_IMMISCIBILITY_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_IMMISCIBILITY_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_MEMORY_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_MEMORY_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_FROZSL_INPUT_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_FROZSL_INPUT_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_FROZSL_POSCAR_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_FROZSL_POSCAR_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_FROZSL_MODES_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_FROZSL_MODES_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_FROZSL_EIGEN_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_FROZSL_EIGEN_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_END_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_END_OUT << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DIELECTRIC_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DIELECTRIC_FILE << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULT GENERIC MPI " << endl;
    aflowrc << "MPI_START_DEFAULT=\"" << AFLOWRC_MPI_START_DEFAULT << "\"" << endl;
    aflowrc << "MPI_STOP_DEFAULT=\"" << AFLOWRC_MPI_STOP_DEFAULT << "\"" << endl;
    aflowrc << "MPI_COMMAND_DEFAULT=\"" << AFLOWRC_MPI_COMMAND_DEFAULT << "\"" << endl;
    aflowrc << "MPI_NCPUS_DEFAULT=" << AFLOWRC_MPI_NCPUS_DEFAULT << endl;
    aflowrc << "MPI_NCPUS_MAX=" << AFLOWRC_MPI_NCPUS_MAX << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS BINARY" << endl;
    aflowrc << "DEFAULT_VASP_GAMMA_BIN=\"" << AFLOWRC_DEFAULT_VASP_GAMMA_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_GAMMA_MPI_BIN=\"" << AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_BIN=\"" << AFLOWRC_DEFAULT_VASP_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_MPI_BIN=\"" << AFLOWRC_DEFAULT_VASP_MPI_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP5_BIN=\"" << AFLOWRC_DEFAULT_VASP5_BIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP5_MPI_BIN=\"" << AFLOWRC_DEFAULT_VASP5_MPI_BIN << "\"" << endl;
    aflowrc << "DEFAULT_AIMS_BIN=\"" << AFLOWRC_DEFAULT_AIMS_BIN << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS POTCARS" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIRECTORIES=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DATE=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DATE << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_SUFFIX=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DATE_POT_LDA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DATE_POT_GGA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA << "\"" << endl;

    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POT_LDA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_LDA << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POT_GGA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_GGA << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POT_PBE=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_PBE << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN << "\"" << endl;
    aflowrc << "DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN=\"" << AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS KPOINTS/DOS" << endl;
    aflowrc << "DEFAULT_BANDS_GRID=" << AFLOWRC_DEFAULT_BANDS_GRID << endl;
    aflowrc << "DEFAULT_BANDS_LATTICE=\"" << AFLOWRC_DEFAULT_BANDS_LATTICE << "\"" << endl;
    aflowrc << "DEFAULT_KSCHEME=\"" << AFLOWRC_DEFAULT_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_KPPRA=" << AFLOWRC_DEFAULT_KPPRA << endl;
    aflowrc << "DEFAULT_STATIC_KSCHEME=\"" << AFLOWRC_DEFAULT_STATIC_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_KPPRA_STATIC=" << AFLOWRC_DEFAULT_KPPRA_STATIC << endl;
    aflowrc << "DEFAULT_KPPRA_ICSD=" << AFLOWRC_DEFAULT_KPPRA_ICSD << endl;
    aflowrc << "DEFAULT_UNARY_BANDS_GRID=" << AFLOWRC_DEFAULT_UNARY_BANDS_GRID << endl;
    aflowrc << "DEFAULT_UNARY_KPPRA=" << AFLOWRC_DEFAULT_UNARY_KPPRA << endl;
    aflowrc << "DEFAULT_UNARY_KPPRA_STATIC=" << AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC << endl;
    aflowrc << "DEFAULT_UNARY_KPPRA_DIELECTRIC=" << AFLOWRC_DEFAULT_UNARY_KPPRA_DIELECTRIC << endl;
    aflowrc << "DEFAULT_PHONONS_KSCHEME=\"" << AFLOWRC_DEFAULT_PHONONS_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_PHONONS_KPPRA=" << AFLOWRC_DEFAULT_PHONONS_KPPRA << endl;
    aflowrc << "DEFAULT_DIELECTRIC_KSCHEME=\"" << AFLOWRC_DEFAULT_DIELECTRIC_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_KPPRA_DIELECTRIC=" << AFLOWRC_DEFAULT_KPPRA_DIELECTRIC << endl;
    aflowrc << "DEFAULT_DOS_EMIN=" << AFLOWRC_DEFAULT_DOS_EMIN << endl;
    aflowrc << "DEFAULT_DOS_EMAX=" << AFLOWRC_DEFAULT_DOS_EMAX << endl;
    aflowrc << "DEFAULT_DOS_SCALE=" << AFLOWRC_DEFAULT_DOS_SCALE << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS PRECISION" << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_LOW=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_MEDIUM=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_NORMAL=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_HIGH=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH << endl;
    aflowrc << "DEFAULT_VASP_PREC_ENMAX_ACCURATE=" << AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE << endl;
    aflowrc << "DEFAULT_VASP_ENMAX_MINIMUM=" << AFLOWRC_DEFAULT_VASP_ENMAX_MINIMUM << endl;
    aflowrc << "DEFAULT_VASP_SPIN_REMOVE_CUTOFF=" << AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF << endl;
    aflowrc << "DEFAULT_VASP_PREC_POTIM=" << AFLOWRC_DEFAULT_VASP_PREC_POTIM << endl;
    aflowrc << "DEFAULT_VASP_PREC_EDIFFG=" << AFLOWRC_DEFAULT_VASP_PREC_EDIFFG << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS OPTIONS " << endl;
    aflowrc << "DEFAULT_VASP_OUT=\"" << AFLOWRC_DEFAULT_VASP_OUT << "\"" << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_INCAR=\"" << AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR << "\"" << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_POSCAR=\"" << AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR << "\"" << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_POTCAR=\"" << AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR << "\"" << endl;
    aflowrc << "DEFAULT_VASP_EXTERNAL_KPOINTS=\"" << AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS << "\"" << endl;
    aflowrc << "DEFAULT_AIMS_EXTERNAL_CONTROL=\"" << AFLOWRC_DEFAULT_AIMS_EXTERNAL_CONTROL << "\"" << endl;
    aflowrc << "DEFAULT_AIMS_EXTERNAL_GEOM=\"" << AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM << "\"" << endl;
    aflowrc << "DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE=\"" << AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE << "\"" << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SIGMA=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_NELM=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM << endl; // CO20200624
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_NELM_STATIC=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM_STATIC << endl; // CO20200624
    aflowrc << "MAX_VASP_NELM=" << AFLOWRC_MAX_VASP_NELM << endl; // CO20200624
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SYM=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SPIN=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_BADER=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_BADER_STATIC=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER_STATIC << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_ELF=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_WAVECAR=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_CHGCAR=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR << endl;
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_LSCOUPLING=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING << endl;

    aflowrc << " " << endl;
    aflowrc << "// AFLOW_LIBRARY AFLOW_PROJECT" << endl;
    aflowrc << "DEFAULT_AFLOW_LIBRARY_DIRECTORIES=\"" << AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_PROJECTS_DIRECTORIES=\"" << AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES << "\"" << endl;
    aflowrc << "DEFAULT_AFLOWDATA_WEB_DIRECTORY=\"" << AFLOWRC_DEFAULT_AFLOWDATA_WEB_DIRECTORY << "\"" << endl; // CO+ME20200731

    aflowrc << " " << endl;
    aflowrc << "// DEFAULT PLATON/FINDSYM" << endl;
    aflowrc << "DEFAULT_PLATON_P_EQUAL=" << AFLOWRC_DEFAULT_PLATON_P_EQUAL << endl;
    aflowrc << "DEFAULT_PLATON_P_EXACT=" << AFLOWRC_DEFAULT_PLATON_P_EXACT << endl;
    aflowrc << "DEFAULT_PLATON_P_ANG=" << AFLOWRC_DEFAULT_PLATON_P_ANG << endl;
    aflowrc << "DEFAULT_PLATON_P_D1=" << AFLOWRC_DEFAULT_PLATON_P_D1 << endl;
    aflowrc << "DEFAULT_PLATON_P_D2=" << AFLOWRC_DEFAULT_PLATON_P_D2 << endl;
    aflowrc << "DEFAULT_PLATON_P_D3=" << AFLOWRC_DEFAULT_PLATON_P_D3 << endl;
    aflowrc << "DEFAULT_FINDSYM_TOL=" << AFLOWRC_DEFAULT_FINDSYM_TOL << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS GNUPLOT" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT_BOLD=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT_BOLD=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_GREEK_FONT=\"" << AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_GREEK_FONT_BOLD=\"" << AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_GREEK_FONT_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_ITALICS << "\"" << endl;
    aflowrc << "DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS=\"" << AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS NHULL" << endl;
    aflowrc << "DEFAULT_NHULL_ALLOWED_DFT_TYPES=\"" << AFLOWRC_DEFAULT_NHULL_ALLOWED_DFT_TYPES << "\"  // comma-separated list of dft_types to include (string match)" << endl;
    aflowrc << "DEFAULT_NHULL_ALLOW_ALL_FORMATION_ENERGIES=" << AFLOWRC_DEFAULT_NHULL_ALLOW_ALL_FORMATION_ENERGIES << " // 0 - false, 1 - true" << endl;
    aflowrc << "DEFAULT_NHULL_COUNT_THRESHOLD_BINARIES=" << AFLOWRC_DEFAULT_NHULL_COUNT_THRESHOLD_BINARIES << " // INT" << endl;
    aflowrc << "DEFAULT_NHULL_PERFORM_OUTLIER_ANALYSIS=" << AFLOWRC_DEFAULT_NHULL_PERFORM_OUTLIER_ANALYSIS << " // 0 - false, 1 - true" << endl;
    aflowrc << "DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES=" << AFLOWRC_DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES << " // INT" << endl;
    aflowrc << "DEFAULT_NHULL_OUTLIER_MULTIPLIER=" << AFLOWRC_DEFAULT_NHULL_OUTLIER_MULTIPLIER << " // DOUBLE" << endl;
    aflowrc << "DEFAULT_NHULL_IGNORE_KNOWN_ILL_CONVERGED=" << AFLOWRC_DEFAULT_NHULL_IGNORE_KNOWN_ILL_CONVERGED << " // 0 - false (NOT recommended), 1 - true" << endl;

    aflowrc << " " << endl; // CO20190628
    aflowrc << "// DEFAULTS GFA" << endl; // CO20190628
    aflowrc << "DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF=" << AFLOWRC_DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF << " // DOUBLE in eV" << endl; // CO20190628

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS ARUN" << endl;
    aflowrc << "ARUN_DIRECTORY_PREFIX=\"" << AFLOWRC_ARUN_DIRECTORY_PREFIX << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS POCC" << endl;
    aflowrc << "DEFAULT_POCC_STRUCTURE_GENERATION_ALGO=\"" << AFLOWRC_DEFAULT_POCC_STRUCTURE_GENERATION_ALGO << "\"" << " // UFF" << endl;
    aflowrc << "DEFAULT_POCC_TEMPERATURE_STRING=\"" << AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING << "\"" << endl;
    aflowrc << "DEFAULT_POCC_EXCLUDE_UNSTABLE=" << AFLOWRC_DEFAULT_POCC_EXCLUDE_UNSTABLE << endl; // ME20210927
    aflowrc << "DEFAULT_POCC_SITE_TOL=" << AFLOWRC_DEFAULT_POCC_SITE_TOL << endl;
    aflowrc << "DEFAULT_POCC_STOICH_TOL=" << AFLOWRC_DEFAULT_POCC_STOICH_TOL << endl;
    aflowrc << "DEFAULT_UFF_BONDING_DISTANCE=" << AFLOWRC_DEFAULT_UFF_BONDING_DISTANCE << endl;
    aflowrc << "DEFAULT_UFF_ENERGY_TOLERANCE=" << AFLOWRC_DEFAULT_UFF_ENERGY_TOLERANCE << endl;
    aflowrc << "DEFAULT_UFF_CLUSTER_RADIUS=" << AFLOWRC_DEFAULT_UFF_CLUSTER_RADIUS << endl;
    aflowrc << "DEFAULT_POCC_RDF_RMAX=" << AFLOWRC_DEFAULT_POCC_RDF_RMAX << endl;
    aflowrc << "DEFAULT_POCC_RDF_NBINS=" << AFLOWRC_DEFAULT_POCC_RDF_NBINS << endl;
    aflowrc << "DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON=" << AFLOWRC_DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON << endl;
    aflowrc << "DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS=" << AFLOWRC_DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS << endl;
    aflowrc << "POCC_FILE_PREFIX=\"" << AFLOWRC_POCC_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_POCC_JSON=\"" << AFLOWRC_DEFAULT_POCC_JSON << "\"" << endl;
    aflowrc << "POCC_APL_OUT_FILE=\"" << AFLOWRC_POCC_APL_OUT_FILE << "\"" << endl; // ME20210927
    aflowrc << "POCC_ALL_SUPERCELLS_FILE=\"" << AFLOWRC_POCC_ALL_SUPERCELLS_FILE << "\"" << endl;
    aflowrc << "POCC_UNIQUE_SUPERCELLS_FILE=\"" << AFLOWRC_POCC_UNIQUE_SUPERCELLS_FILE << "\"" << endl;
    aflowrc << "POCC_ALL_HNF_MATRICES_FILE=\"" << AFLOWRC_POCC_ALL_HNF_MATRICES_FILE << "\"" << endl;
    aflowrc << "POCC_ALL_SITE_CONFIGURATIONS_FILE=\"" << AFLOWRC_POCC_ALL_SITE_CONFIGURATIONS_FILE << "\"" << endl;
    aflowrc << "POCC_DOSCAR_FILE=\"" << AFLOWRC_POCC_DOSCAR_FILE << "\"" << endl;
    aflowrc << "POCC_PHDOSCAR_FILE=\"" << AFLOWRC_POCC_PHDOSCAR_FILE << "\"" << endl; // ME20210927
    aflowrc << "POCC_ANIONS_LIST=\"" << AFLOWRC_POCC_ANIONS_LIST << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS APL" << endl;
    aflowrc << "DEFAULT_APL_PREC=\"" << AFLOWRC_DEFAULT_APL_PREC << "\"" << endl;
    aflowrc << "DEFAULT_APL_ENGINE=\"" << AFLOWRC_DEFAULT_APL_ENGINE << "\"" << endl;
    aflowrc << "DEFAULT_APL_HIBERNATE=" << AFLOWRC_DEFAULT_APL_HIBERNATE << endl;
    aflowrc << "DEFAULT_APL_MINSHELL=" << AFLOWRC_DEFAULT_APL_MINSHELL << endl;
    aflowrc << "DEFAULT_APL_MINATOMS=" << AFLOWRC_DEFAULT_APL_MINATOMS << endl;
    aflowrc << "DEFAULT_APL_POLAR=" << AFLOWRC_DEFAULT_APL_POLAR << endl;
    aflowrc << "DEFAULT_APL_DMAG=" << AFLOWRC_DEFAULT_APL_DMAG << endl;
    aflowrc << "DEFAULT_APL_DXYZONLY=" << AFLOWRC_DEFAULT_APL_DXYZONLY << endl;
    aflowrc << "DEFAULT_APL_DSYMMETRIZE=" << AFLOWRC_DEFAULT_APL_DSYMMETRIZE << endl;
    aflowrc << "DEFAULT_APL_DINEQUIV_ONLY=" << AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY << endl;
    aflowrc << "DEFAULT_APL_DPM=\"" << AFLOWRC_DEFAULT_APL_DPM << "\"" << endl;
    aflowrc << "DEFAULT_APL_RELAX=" << AFLOWRC_DEFAULT_APL_RELAX << endl;
    aflowrc << "DEFAULT_APL_RELAX_COMMENSURATE=" << AFLOWRC_DEFAULT_APL_RELAX_COMMENSURATE << endl; // ME20200427
    aflowrc << "DEFAULT_APL_ZEROSTATE=" << AFLOWRC_DEFAULT_APL_ZEROSTATE << endl;
    aflowrc << "DEFAULT_APL_ZEROSTATE_CHGCAR=" << AFLOWRC_DEFAULT_APL_ZEROSTATE_CHGCAR << endl; // ME20191029
    aflowrc << "DEFAULT_APL_USE_LEPSILON=" << AFLOWRC_DEFAULT_APL_USE_LEPSILON << endl;
    aflowrc << "DEFAULT_APL_FREQFORMAT=\"" << AFLOWRC_DEFAULT_APL_FREQFORMAT << "\"" << endl;
    aflowrc << "DEFAULT_APL_DC=" << AFLOWRC_DEFAULT_APL_DC << endl;
    aflowrc << "DEFAULT_APL_DCPATH=\"" << AFLOWRC_DEFAULT_APL_DCPATH << "\"" << endl;
    aflowrc << "DEFAULT_APL_DCPOINTS=" << AFLOWRC_DEFAULT_APL_DCPOINTS << endl; // CO20181226
    aflowrc << "DEFAULT_APL_DOS=" << AFLOWRC_DEFAULT_APL_DOS << endl;
    aflowrc << "DEFAULT_APL_DOSMETHOD=\"" << AFLOWRC_DEFAULT_APL_DOSMETHOD << "\"" << endl;
    aflowrc << "DEFAULT_APL_DOSMESH=\"" << AFLOWRC_DEFAULT_APL_DOSMESH << "\"" << endl;
    aflowrc << "DEFAULT_APL_DOSPOINTS=" << AFLOWRC_DEFAULT_APL_DOSPOINTS << endl;
    aflowrc << "DEFAULT_APL_DOSSMEAR=" << AFLOWRC_DEFAULT_APL_DOSSMEAR << endl;
    aflowrc << "DEFAULT_APL_DOS_PROJECT=" << AFLOWRC_DEFAULT_APL_DOS_PROJECT << endl; // ME20200421
    aflowrc << "DEFAULT_APL_TP=" << AFLOWRC_DEFAULT_APL_TP << endl;
    aflowrc << "DEFAULT_APL_DISPLACEMENTS=" << AFLOWRC_DEFAULT_APL_DISPLACEMENTS << endl; // ME20200421
    aflowrc << "DEFAULT_APL_TPT=\"" << AFLOWRC_DEFAULT_APL_TPT << "\"" << endl;
    aflowrc << "DEFAULT_APL_GVEL=" << AFLOWRC_DEFAULT_APL_GVEL << endl; // ME20200517
    aflowrc << "DEFAULT_APL_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_APL_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_APL_OUT_FILE=\"" << AFLOWRC_DEFAULT_APL_OUT_FILE << "\"" << endl; // ME20210927
    aflowrc << "DEFAULT_APL_PDIS_FILE=\"" << AFLOWRC_DEFAULT_APL_PDIS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PDOS_FILE=\"" << AFLOWRC_DEFAULT_APL_PDOS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_THERMO_FILE=\"" << AFLOWRC_DEFAULT_APL_THERMO_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_THERMO_JSON=\"" << AFLOWRC_DEFAULT_APL_THERMO_JSON << "\"" << endl; // ME20211019
    aflowrc << "DEFAULT_APL_DYNMAT_FILE=\"" << AFLOWRC_DEFAULT_APL_DYNMAT_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_HARMIFC_FILE=\"" << AFLOWRC_DEFAULT_APL_HARMIFC_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_POLAR_FILE=\"" << AFLOWRC_DEFAULT_APL_POLAR_FILE << "\"" << endl; // ME20200415
    aflowrc << "DEFAULT_APL_HSKPTS_FILE=\"" << AFLOWRC_DEFAULT_APL_HSKPTS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_MSQRDISP_FILE=\"" << AFLOWRC_DEFAULT_APL_MSQRDISP_FILE << "\"" << endl; // ME20200329
    aflowrc << "DEFAULT_APL_GVEL_FILE=\"" << AFLOWRC_DEFAULT_APL_GVEL_FILE << "\"" << endl; // ME20200517
    // ME20190614 BEGIN
    aflowrc << "DEFAULT_APL_PHDOSCAR_FILE=\"" << AFLOWRC_DEFAULT_APL_PHDOSCAR_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PHPOSCAR_FILE=\"" << AFLOWRC_DEFAULT_APL_PHPOSCAR_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PHKPOINTS_FILE=\"" << AFLOWRC_DEFAULT_APL_PHKPOINTS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PHEIGENVAL_FILE=\"" << AFLOWRC_DEFAULT_APL_PHEIGENVAL_FILE << "\"" << endl;
    // ME20190614 END
    aflowrc << "DEFAULT_APL_STATE_FILE=\"" << AFLOWRC_DEFAULT_APL_STATE_FILE << "\"" << endl; // ME20200224
    // ME20200329 BEGIN
    aflowrc << "DEFAULT_APL_ADISP_SCENE_FORMAT=\"" << AFLOWRC_DEFAULT_APL_ADISP_SCENE_FORMAT << "\"" << endl;
    aflowrc << "DEFAULT_APL_ADISP_AMPLITUDE=" << AFLOWRC_DEFAULT_APL_ADISP_AMPLITUDE << endl;
    aflowrc << "DEFAULT_APL_ADISP_NSTEPS=" << AFLOWRC_DEFAULT_APL_ADISP_NSTEPS << endl;
    aflowrc << "DEFAULT_APL_ADISP_NPERIODS=" << AFLOWRC_DEFAULT_APL_ADISP_NPERIODS << endl;
    // ME20200329 END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS QHA" << endl;
    aflowrc << "DEFAULT_QHA_MODE=\"" << AFLOWRC_DEFAULT_QHA_MODE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_EOS=" << AFLOWRC_DEFAULT_QHA_EOS << endl;
    aflowrc << "DEFAULT_QHA_EOS_DISTORTION_RANGE=\"" << AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_EOS_MODEL=\"" << AFLOWRC_DEFAULT_QHA_EOS_MODEL << "\"" << endl; // AS20200818
    aflowrc << "DEFAULT_QHA_GP_DISTORTION=" << AFLOWRC_DEFAULT_QHA_GP_DISTORTION << endl;
    aflowrc << "DEFAULT_QHA_TAYLOR_EXPANSION_ORDER=" << AFLOWRC_DEFAULT_QHA_TAYLOR_EXPANSION_ORDER << endl; // AS20200602
    aflowrc << "DEFAULT_QHA_INCLUDE_ELEC_CONTRIB=" << AFLOWRC_DEFAULT_QHA_INCLUDE_ELEC_CONTRIB << endl;
    aflowrc << "DEFAULT_QHA_SOMMERFELD_EXPANSION=" << AFLOWRC_DEFAULT_QHA_SOMMERFELD_EXPANSION << endl; // AS20200528
    aflowrc << "DEFAULT_QHA_PDIS_T=\"" << AFLOWRC_DEFAULT_QHA_PDIS_T << "\"" << endl;
    // AS20200508 BEGIN
    aflowrc << "DEFAULT_QHA_GP_FINITE_DIFF=" << AFLOWRC_DEFAULT_QHA_GP_FINITE_DIFF << endl;
    aflowrc << "DEFAULT_QHA_IGNORE_IMAGINARY=" << AFLOWRC_DEFAULT_QHA_IGNORE_IMAGINARY << endl;
    aflowrc << "DEFAULT_QHA_RELAX_IONS_CELL=" << AFLOWRC_DEFAULT_QHA_RELAX_IONS_CELL << endl; // AS20201123
    aflowrc << "DEFAULT_QHA_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_QHA_FILE_PREFIX << "\"" << endl;
    // AS20200709 BEGIN
    aflowrc << "DEFAULT_QHA3P_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_QHA3P_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_QHANP_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_QHANP_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_SCQHA_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_SCQHA_FILE_PREFIX << "\"" << endl;
    // AS20200709 END
    aflowrc << "DEFAULT_QHA_GP_PATH_FILE=\"" << AFLOWRC_DEFAULT_QHA_GP_PATH_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_GP_MESH_FILE=\"" << AFLOWRC_DEFAULT_QHA_GP_MESH_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_GP_AVG_FILE=\"" << AFLOWRC_DEFAULT_QHA_GP_AVG_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_THERMO_FILE=\"" << AFLOWRC_DEFAULT_QHA_THERMO_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_FREQS_FILE=\"" << AFLOWRC_DEFAULT_QHA_FREQS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_FVT_FILE=\"" << AFLOWRC_DEFAULT_QHA_FVT_FILE << "\"" << endl;
    // AS20200508 END
    aflowrc << "DEFAULT_QHA_COEFF_FILE=\"" << AFLOWRC_DEFAULT_QHA_COEFF_FILE << "\"" << endl; // AS20210517
    aflowrc << "DEFAULT_QHA_IMAG_FILE=\"" << AFLOWRC_DEFAULT_QHA_IMAG_FILE << "\"" << endl; // AS20210517
    aflowrc << "DEFAULT_QHA_PDIS_FILE=\"" << AFLOWRC_DEFAULT_QHA_PDIS_FILE << "\"" << endl; // AS20201022
    aflowrc << "DEFAULT_QHA_PDOS_FILE=\"" << AFLOWRC_DEFAULT_QHA_PDOS_FILE << "\"" << endl; // AS20201201
    aflowrc << "DEFAULT_QHA_KPOINTS_FILE=\"" << AFLOWRC_DEFAULT_QHA_KPOINTS_FILE << "\"" << endl; // AS20201112
    // AS20210914 BEGIN
    aflowrc << "DEFAULT_POCC_QHA_THERMO_FILE=\"" << AFLOWRC_DEFAULT_POCC_QHA_THERMO_FILE << "\"" << endl;
    aflowrc << "DEFAULT_POCC_QHA_AVGTHERMO_FILE=\"" << AFLOWRC_DEFAULT_POCC_QHA_AVGTHERMO_FILE << "\"" << endl;
    // AS20210914 END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS AAPL" << endl;
    aflowrc << "DEFAULT_AAPL_BTE=\"" << AFLOWRC_DEFAULT_AAPL_BTE << "\"" << endl;
    //[ME20181226]aflowrc << "DEFAULT_AAPL_BZMETHOD=\"" << AFLOWRC_DEFAULT_AAPL_BZMETHOD << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_FOURTH_ORDER=" << AFLOWRC_DEFAULT_AAPL_FOURTH_ORDER << endl;
    aflowrc << "DEFAULT_AAPL_CUT_RAD=\"" << AFLOWRC_DEFAULT_AAPL_CUT_RAD << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_CUT_SHELL=\"" << AFLOWRC_DEFAULT_AAPL_CUT_SHELL << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_THERMALGRID=\"" << AFLOWRC_DEFAULT_AAPL_THERMALGRID << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_TCT=\"" << AFLOWRC_DEFAULT_AAPL_TCT << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_SUMRULE=" << AFLOWRC_DEFAULT_AAPL_SUMRULE << endl;
    aflowrc << "DEFAULT_AAPL_SUMRULE_MAX_ITER=" << AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER << endl;
    aflowrc << "DEFAULT_AAPL_MIXING_COEFFICIENT=" << AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT << endl;
    aflowrc << "DEFAULT_AAPL_ISOTOPE=" << AFLOWRC_DEFAULT_AAPL_ISOTOPE << endl;
    aflowrc << "DEFAULT_AAPL_BOUNDARY=" << AFLOWRC_DEFAULT_AAPL_BOUNDARY << endl;
    aflowrc << "DEFAULT_AAPL_CUMULATIVE=" << AFLOWRC_DEFAULT_AAPL_CUMULATIVEK << endl;
    aflowrc << "DEFAULT_AAPL_NANO_SIZE=" << AFLOWRC_DEFAULT_AAPL_NANO_SIZE << endl;
    aflowrc << "DEFAULT_AAPL_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_AAPL_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_IRRQPTS_FILE=\"" << AFLOWRC_DEFAULT_AAPL_IRRQPTS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_GVEL_FILE=\"" << AFLOWRC_DEFAULT_AAPL_GVEL_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_PS_FILE=\"" << AFLOWRC_DEFAULT_AAPL_PS_FILE << "\"" << endl; // ME20191104
    aflowrc << "DEFAULT_AAPL_GRUENEISEN_FILE=\"" << AFLOWRC_DEFAULT_AAPL_GRUENEISEN_FILE << "\"" << endl; // ME20191104
    aflowrc << "DEFAULT_AAPL_RATES_FILE=\"" << AFLOWRC_DEFAULT_AAPL_RATES_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_RATES_3RD_FILE=\"" << AFLOWRC_DEFAULT_AAPL_RATES_3RD_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_RATES_4TH_FILE=\"" << AFLOWRC_DEFAULT_AAPL_RATES_4TH_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_ISOTOPE_FILE=\"" << AFLOWRC_DEFAULT_AAPL_ISOTOPE_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_BOUNDARY_FILE=\"" << AFLOWRC_DEFAULT_AAPL_BOUNDARY_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_TCOND_FILE=\"" << AFLOWRC_DEFAULT_AAPL_TCOND_FILE << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS AEL" << endl;
    aflowrc << "DEFAULT_AEL_STRAIN_SYMMETRY=" << AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY << endl;
    aflowrc << "DEFAULT_AEL_NNORMAL_STRAINS=" << AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS << endl;
    aflowrc << "DEFAULT_AEL_NSHEAR_STRAINS=" << AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS << endl;
    aflowrc << "DEFAULT_AEL_NORMAL_STRAIN_STEP=" << AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP << endl;
    aflowrc << "DEFAULT_AEL_SHEAR_STRAIN_STEP=" << AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP << endl;
    aflowrc << "DEFAULT_AEL_ORIGIN_STRAIN_CALC=" << AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC << endl;
    aflowrc << "DEFAULT_AEL_ORIGIN_STRAIN_FIT=" << AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT << endl;
    aflowrc << "DEFAULT_AEL_RELAXED_STRUCT_FIT=" << AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT << endl;
    aflowrc << "DEFAULT_AEL_NEG_STRAINS=" << AFLOWRC_DEFAULT_AEL_NEG_STRAINS << endl;
    aflowrc << "DEFAULT_AEL_NIND_STRAIN_DIRS=" << AFLOWRC_DEFAULT_AEL_NIND_STRAIN_DIRS << endl;
    aflowrc << "DEFAULT_AEL_VASPSYM=" << AFLOWRC_DEFAULT_AEL_VASPSYM << endl;
    aflowrc << "DEFAULT_AEL_PRECACC_ALGONORM=" << AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM << endl;
    aflowrc << "DEFAULT_AEL_VASPRUNXML_STRESS=" << AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS << endl;
    aflowrc << "DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS=" << AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS << endl;
    aflowrc << "DEFAULT_AEL_SKIP_ARUNS_MAX=" << AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX << endl;
    aflowrc << "DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY=" << AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY << endl;
    aflowrc << "DEFAULT_AEL_SYMMETRIZE=" << AFLOWRC_DEFAULT_AEL_SYMMETRIZE << endl;
    aflowrc << "DEFAULT_AEL_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_AEL_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_AEL_WRITE_FULL_RESULTS=" << AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS << endl;
    aflowrc << "DEFAULT_AEL_DIRNAME_ARUN=" << AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS AGL" << endl;
    aflowrc << "DEFAULT_AGL_AEL_POISSON_RATIO=" << AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO << endl;
    aflowrc << "DEFAULT_AGL_NSTRUCTURES=" << AFLOWRC_DEFAULT_AGL_NSTRUCTURES << endl;
    aflowrc << "DEFAULT_AGL_STRAIN_STEP=" << AFLOWRC_DEFAULT_AGL_STRAIN_STEP << endl;
    aflowrc << "DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS=" << AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS << endl;
    aflowrc << "DEFAULT_AGL_SKIP_ARUNS_MAX=" << AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX << endl;
    aflowrc << "DEFAULT_AGL_NTEMPERATURE=" << AFLOWRC_DEFAULT_AGL_NTEMPERATURE << endl;
    aflowrc << "DEFAULT_AGL_STEMPERATURE=" << AFLOWRC_DEFAULT_AGL_STEMPERATURE << endl;
    aflowrc << "DEFAULT_AGL_NPRESSURE=" << AFLOWRC_DEFAULT_AGL_NPRESSURE << endl;
    aflowrc << "DEFAULT_AGL_SPRESSURE=" << AFLOWRC_DEFAULT_AGL_SPRESSURE << endl;
    aflowrc << "DEFAULT_AGL_POISSON_RATIO=" << AFLOWRC_DEFAULT_AGL_POISSON_RATIO << endl;
    aflowrc << "DEFAULT_AGL_IEOS=" << AFLOWRC_DEFAULT_AGL_IEOS << endl;
    aflowrc << "DEFAULT_AGL_IDEBYE=" << AFLOWRC_DEFAULT_AGL_IDEBYE << endl;
    aflowrc << "DEFAULT_AGL_FIT_TYPE=" << AFLOWRC_DEFAULT_AGL_FIT_TYPE << endl;
    aflowrc << "DEFAULT_AGL_CHECK_EV_CONCAVITY=" << AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY << endl;
    aflowrc << "DEFAULT_AGL_CHECK_EV_MIN=" << AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN << endl;
    aflowrc << "DEFAULT_AGL_HUGONIOT_CALC=" << AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC << endl;
    aflowrc << "DEFAULT_AGL_HUGONIOT_EXTRAPOLATE=" << AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE << endl;
    aflowrc << "DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE=" << AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE << endl;
    aflowrc << "DEFAULT_AGL_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_AGL_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_AGL_WRITE_FULL_RESULTS=" << AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS << endl;
    aflowrc << "DEFAULT_AGL_DIRNAME_ARUN=" << AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN << endl;
    aflowrc << "DEFAULT_AGL_WRITE_GIBBS_INPUT=" << AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT << endl;
    aflowrc << "DEFAULT_AGL_PLOT_RESULTS=" << AFLOWRC_DEFAULT_AGL_PLOT_RESULTS << endl;

    // SD20220323 - QCA START
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS QCA " << endl;
    aflowrc << "DEFAULT_QCA_MIN_SLEEP_SECONDS" << AFLOWRC_DEFAULT_QCA_MIN_SLEEP_SECONDS << endl;
    aflowrc << "DEFAULT_QCA_MAX_NUM_ATOMS" << AFLOWRC_DEFAULT_QCA_MAX_NUM_ATOMS << endl;
    aflowrc << "DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS" << AFLOWRC_DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS << endl;
    aflowrc << "DEFAULT_QCA_CV_CUTOFF" << AFLOWRC_DEFAULT_QCA_CV_CUTOFF << endl;
    aflowrc << "DEFAULT_QCA_CONC_NPTS" << AFLOWRC_DEFAULT_QCA_CONC_NPTS << endl;
    aflowrc << "DEFAULT_QCA_TEMP_NPTS" << AFLOWRC_DEFAULT_QCA_TEMP_NPTS << endl;
    aflowrc << "DEFAULT_QCA_TEMP_MIN" << AFLOWRC_DEFAULT_QCA_TEMP_MIN << endl;
    aflowrc << "DEFAULT_QCA_TEMP_MAX" << AFLOWRC_DEFAULT_QCA_TEMP_MAX << endl;
    aflowrc << "DEFAULT_QCA_TEMP_MIN_LIMIT" << AFLOWRC_DEFAULT_QCA_TEMP_MIN_LIMIT << endl;
    aflowrc << "DEFAULT_QCA_PRINT" << AFLOWRC_DEFAULT_QCA_PRINT << endl;
    // SD20220323 - QCA END

    // RF20200413 START
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS CCE" << endl;
    aflowrc << "DEFAULT_CCE_OX_METHOD=" << AFLOWRC_DEFAULT_CCE_OX_METHOD << "" << "  // 1 - ELECTRONEGATIVITY_ALLEN, 2 - BADER" << endl;
    aflowrc << "DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION=" << AFLOWRC_DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION << "" << "  // 0.4 Ang tolerance between shortest and longest bonds when testing for multi-anion compound" << endl;
    aflowrc << "DEFAULT_CCE_OX_TOL=" << AFLOWRC_DEFAULT_CCE_OX_TOL << "" << "  // sum of oxidation states might not be exactly zero due to numerics" << endl;
    aflowrc << "DEFAULT_CCE_PEROX_CUTOFF=" << AFLOWRC_DEFAULT_CCE_PEROX_CUTOFF << "" << "  // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang" << endl;
    aflowrc << "DEFAULT_CCE_SUPEROX_CUTOFF=" << AFLOWRC_DEFAULT_CCE_SUPEROX_CUTOFF << "" << "  // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang" << endl;
    aflowrc << "DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF=" << AFLOWRC_DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF << "" << "  // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    aflowrc << "DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF=" << AFLOWRC_DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF << "" << "  // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    // RF20200413 END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS XTALFINDER" << endl;
    aflowrc << "DEFAULT_XTALFINDER_MISFIT_MATCH=" << AFLOWRC_DEFAULT_XTALFINDER_MISFIT_MATCH << " // values below this threshold: similar structures have similar properties" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_MISFIT_FAMILY=" << AFLOWRC_DEFAULT_XTALFINDER_MISFIT_FAMILY << " // values above this threshold: matched structures do not have similar properties" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_SUPERCELL_METHOD=" << AFLOWRC_DEFAULT_XTALFINDER_SUPERCELL_METHOD << " // // supercell method for comparing (robust, but slow, superceded by transformation method)" << endl; // DX20201223
    aflowrc << "DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING=" << AFLOWRC_DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING << " // factor that divides minimum interatomic distance" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_MATERIAL=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL << " // results file prefix" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_STRUCTURE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE << " // results file prefix" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_DUPLICATE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_DUPLICATE << " // results file prefix" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE << " // results file prefix" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE << " // results file prefix" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE << " // results file prefix" << endl; // DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE << " // results file prefix" << endl; // DX20201118

    // DX20200720 - START
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS ANRL" << endl;
    aflowrc << "DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL=" << AFLOWRC_DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL << " // tolerance for equivalent Wyckoff coordinates" << endl;
    // DX20200720 - END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS CORE" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_BEEP=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP << " // Celsius" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_HALT=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT << " // Celsius" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_REFRESH=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH << " // seconds" << endl;

    aflowrc << "SECONDS_SLEEP_VASP_COMPLETION=" << AFLOWRC_SECONDS_SLEEP_VASP_COMPLETION << " // seconds" << endl; // CO20201111
    aflowrc << "SECONDS_SLEEP_VASP_MONITOR=" << AFLOWRC_SECONDS_SLEEP_VASP_MONITOR << " // seconds" << endl; // CO20201111
    aflowrc << "SECONDS_STALE_OUTCAR=" << AFLOWRC_SECONDS_STALE_OUTCAR << " // seconds" << endl; // CO20201111
    aflowrc << "BYTES_MAX_VASP_OUT=" << AFLOWRC_BYTES_MAX_VASP_OUT << " // bytes" << endl; // CO20201111
    aflowrc << "MEMORY_MAX_USAGE_RAM=" << AFLOWRC_MEMORY_MAX_USAGE_RAM << " // bytes" << endl; // CO20201111
    aflowrc << "MEMORY_MAX_USAGE_SWAP=" << AFLOWRC_MEMORY_MAX_USAGE_SWAP << " // bytes" << endl; // CO20201111
    aflowrc << "FILE_VASP_MONITOR=" << AFLOWRC_FILE_VASP_MONITOR << " // monitor file postfix" << endl; // CO20201111
    aflowrc << "INTEL_COMPILER_PATHS=" << AFLOWRC_INTEL_COMPILER_PATHS << " // comma-separated paths to search (for sourcing)" << endl; // CO20201111

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS MACHINE DEPENDENT MPI" << endl;
    aflowrc << "MPI_OPTIONS_DUKE_BETA_MPICH=\"" << AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << "  // DUKE_BETA_MPICH" << endl;
    aflowrc << "MPI_COMMAND_DUKE_BETA_MPICH=\"" << AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH << "\"" << "  // DUKE_BETA_MPICH" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_BETA_MPICH=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH << "\"" << "  // DUKE_BETA_MPICH" << endl;

    aflowrc << "MPI_OPTIONS_DUKE_BETA_OPENMPI=\"" << AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI << "\"" << "  // DUKE_BETA_OPENMPI" << endl;
    aflowrc << "MPI_COMMAND_DUKE_BETA_OPENMPI=\"" << AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI << "\"" << "  // DUKE_BETA_OPENMPI" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_BETA_OPENMPI=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI << "\"" << "  // DUKE_BETA_OPENMPI" << endl;

    aflowrc << "MPI_OPTIONS_DUKE_MATERIALS=\"" << AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS << "\"" << "  // DUKE_MATERIALS" << endl;
    aflowrc << "MPI_COMMAND_DUKE_MATERIALS=\"" << AFLOWRC_MPI_COMMAND_DUKE_MATERIALS << "\"" << "  // DUKE_MATERIALS" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_MATERIALS=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS << "\"" << "  // DUKE_MATERIALS" << endl;

    aflowrc << "MPI_OPTIONS_DUKE_AFLOWLIB=\"" << AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB << "\"" << "  // DUKE_AFLOWLIB" << endl;
    aflowrc << "MPI_COMMAND_DUKE_AFLOWLIB=\"" << AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB << "\"" << "  // DUKE_AFLOWLIB" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_AFLOWLIB=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB << "\"" << "  // DUKE_AFLOWLIB" << endl;

    aflowrc << "MPI_OPTIONS_DUKE_HABANA=\"" << AFLOWRC_MPI_OPTIONS_DUKE_HABANA << "\"" << "  // DUKE_HABANA" << endl;
    aflowrc << "MPI_COMMAND_DUKE_HABANA=\"" << AFLOWRC_MPI_COMMAND_DUKE_HABANA << "\"" << "  // DUKE_HABANA" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_HABANA=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA << "\"" << "  // DUKE_HABANA" << endl;

    aflowrc << "MPI_OPTIONS_DUKE_QRATS_MPICH=\"" << AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH << "\"" << "  // DUKE_QRATS_MPICH" << endl;
    aflowrc << "MPI_COMMAND_DUKE_QRATS_MPICH=\"" << AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH << "\"" << "  // DUKE_QRATS_MPICH" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_QRATS_MPICH=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH << "\"" << "  // DUKE_QRATS_MPICH" << endl;

    aflowrc << "MPI_OPTIONS_DUKE_QFLOW_OPENMPI=\"" << AFLOWRC_MPI_OPTIONS_DUKE_QFLOW_OPENMPI << "\"" << "  // DUKE_QFLOW_OPENMPI" << endl;
    aflowrc << "MPI_COMMAND_DUKE_QFLOW_OPENMPI=\"" << AFLOWRC_MPI_COMMAND_DUKE_QFLOW_OPENMPI << "\"" << "  // DUKE_QFLOW_OPENMPI" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << "\"" << "  // DUKE_QFLOW_OPENMPI" << endl;

    // CO20201220 X START
    aflowrc << "MPI_OPTIONS_DUKE_X_X=\"" << AFLOWRC_MPI_OPTIONS_DUKE_X_X << "\"" << "  // DUKE_X_X" << endl;
    aflowrc << "MPI_COMMAND_DUKE_X_X=\"" << AFLOWRC_MPI_COMMAND_DUKE_X_X << "\"" << "  // DUKE_X_X" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_X_X=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_X_X << "\"" << "  // DUKE_X_X" << endl;
    aflowrc << "MPI_OPTIONS_DUKE_X_CRAY=\"" << AFLOWRC_MPI_OPTIONS_DUKE_X_CRAY << "\"" << "  // DUKE_X_CRAY" << endl;
    aflowrc << "MPI_COMMAND_DUKE_X_CRAY=\"" << AFLOWRC_MPI_COMMAND_DUKE_X_CRAY << "\"" << "  // DUKE_X_CRAY" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_X_CRAY=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_X_CRAY << "\"" << "  // DUKE_X_CRAY" << endl;
    aflowrc << "MPI_OPTIONS_DUKE_X_OLDCRAY=\"" << AFLOWRC_MPI_OPTIONS_DUKE_X_OLDCRAY << "\"" << "  // DUKE_X_OLDCRAY" << endl;
    aflowrc << "MPI_COMMAND_DUKE_X_OLDCRAY=\"" << AFLOWRC_MPI_COMMAND_DUKE_X_OLDCRAY << "\"" << "  // DUKE_X_OLDCRAY" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_X_OLDCRAY=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_X_OLDCRAY << "\"" << "  // DUKE_X_OLDCRAY" << endl;
    aflowrc << "MPI_OPTIONS_DUKE_X_SMB=\"" << AFLOWRC_MPI_OPTIONS_DUKE_X_SMB << "\"" << "  // DUKE_X_SMB" << endl;
    aflowrc << "MPI_COMMAND_DUKE_X_SMB=\"" << AFLOWRC_MPI_COMMAND_DUKE_X_SMB << "\"" << "  // DUKE_X_SMB" << endl;
    aflowrc << "MPI_BINARY_DIR_DUKE_X_SMB=\"" << AFLOWRC_MPI_BINARY_DIR_DUKE_X_SMB << "\"" << "  // DUKE_X_SMB" << endl;
    // CO20201220 X STOP

    // CO20220818 JHU_ROCKFISH START
    aflowrc << "MPI_OPTIONS_JHU_ROCKFISH=\"" << AFLOWRC_MPI_OPTIONS_JHU_ROCKFISH << "\"" << "  // JHU_ROCKFISH" << endl;
    aflowrc << "MPI_COMMAND_JHU_ROCKFISH=\"" << AFLOWRC_MPI_COMMAND_JHU_ROCKFISH << "\"" << "  // JHU_ROCKFISH" << endl;
    aflowrc << "MPI_BINARY_DIR_JHU_ROCKFISH=\"" << AFLOWRC_MPI_BINARY_DIR_JHU_ROCKFISH << "\"" << "  // JHU_ROCKFISH" << endl;
    // CO20220818 JHU_ROCKFISH STOP

    // DX20190509 - MACHINE001 - START
    aflowrc << "MPI_OPTIONS_MACHINE001=\"" << AFLOWRC_MPI_OPTIONS_MACHINE001 << "\"" << "// MACHINE001" << endl;
    aflowrc << "MPI_COMMAND_MACHINE001=\"" << AFLOWRC_MPI_COMMAND_MACHINE001 << "\"" << "// MACHINE001" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE001=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE001 << "\"" << "// MACHINE001" << endl;
    // DX20190509 - MACHINE001 - END

    // DX20190509 - MACHINE002 - START
    aflowrc << "MPI_OPTIONS_MACHINE002=\"" << AFLOWRC_MPI_OPTIONS_MACHINE002 << "\"" << "// MACHINE002" << endl;
    aflowrc << "MPI_COMMAND_MACHINE002=\"" << AFLOWRC_MPI_COMMAND_MACHINE002 << "\"" << "// MACHINE002" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE002=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE002 << "\"" << "// MACHINE002" << endl;
    // DX20190509 - MACHINE002 - START

    // DX20201005 - MACHINE003 - START
    aflowrc << "MPI_OPTIONS_MACHINE003=\"" << AFLOWRC_MPI_OPTIONS_MACHINE003 << "\"" << "// MACHINE003" << endl;
    aflowrc << "MPI_COMMAND_MACHINE003=\"" << AFLOWRC_MPI_COMMAND_MACHINE003 << "\"" << "// MACHINE003" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE003=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE003 << "\"" << "// MACHINE003" << endl;
    // DX20201005 - MACHINE003 - START

    // DX20211011 - MACHINE004 - START
    aflowrc << "MPI_OPTIONS_MACHINE004=\"" << AFLOWRC_MPI_OPTIONS_MACHINE004 << "\"" << "// MACHINE004" << endl;
    aflowrc << "MPI_COMMAND_MACHINE004=\"" << AFLOWRC_MPI_COMMAND_MACHINE004 << "\"" << "// MACHINE004" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE004=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE004 << "\"" << "// MACHINE004" << endl;
    // DX20211011 - MACHINE004 - START

    // DX20190107 - CMU EULER - START
    aflowrc << "MPI_OPTIONS_CMU_EULER=\"" << AFLOWRC_MPI_OPTIONS_CMU_EULER << "\"" << "// CMU_EULER" << endl;
    aflowrc << "MPI_COMMAND_CMU_EULER=\"" << AFLOWRC_MPI_COMMAND_CMU_EULER << "\"" << "// CMU_EULER" << endl;
    aflowrc << "MPI_BINARY_DIR_CMU_EULER=\"" << AFLOWRC_MPI_BINARY_DIR_CMU_EULER << "\"" << "// CMU_EULER" << endl;
    // DX20190107 - CMU EULER - END

    aflowrc << "MPI_OPTIONS_MPCDF_EOS=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_EOS=\"" << AFLOWRC_MPI_COMMAND_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_EOS=\"" << AFLOWRC_MPI_NCPUS_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_EOS=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI        // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_EOS=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;

    aflowrc << "MPI_OPTIONS_MPCDF_DRACO=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_DRACO=\"" << AFLOWRC_MPI_COMMAND_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_DRACO=\"" << AFLOWRC_MPI_NCPUS_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_DRACO=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI        // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_DRACO=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;

    aflowrc << "MPI_OPTIONS_MPCDF_COBRA=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_COBRA=\"" << AFLOWRC_MPI_COMMAND_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_COBRA=\"" << AFLOWRC_MPI_NCPUS_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_COBRA=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI        // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_COBRA=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;

    aflowrc << "MPI_OPTIONS_MPCDF_HYDRA=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_HYDRA=\"" << AFLOWRC_MPI_COMMAND_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_HYDRA=\"" << AFLOWRC_MPI_NCPUS_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_HYDRA=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI        // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_HYDRA=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;

    aflowrc << "MPI_OPTIONS_FULTON_MARYLOU=\"" << AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU << "\"" << "  // FULTON_MARYLOU" << endl;
    aflowrc << "MPI_COMMAND_FULTON_MARYLOU=\"" << AFLOWRC_MPI_COMMAND_FULTON_MARYLOU << "\"" << "  // FULTON_MARYLOU" << endl;
    aflowrc << "MPI_BINARY_DIR_FULTON_MARYLOU=\"" << AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU << "\"" << "  // FULTON_MARYLOU" << endl;

    aflowrc << "MPI_OPTIONS_MACHINE1=\"" << AFLOWRC_MPI_OPTIONS_MACHINE1 << "\"" << "  // MACHINE1" << endl;
    aflowrc << "MPI_COMMAND_MACHINE1=\"" << AFLOWRC_MPI_COMMAND_MACHINE1 << "\"" << "  // MACHINE1" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE1=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE1 << "\"" << "  // MACHINE1" << endl;

    aflowrc << "MPI_OPTIONS_MACHINE2=\"" << AFLOWRC_MPI_OPTIONS_MACHINE2 << "\"" << "  // MACHINE2" << endl;
    aflowrc << "MPI_COMMAND_MACHINE2=\"" << AFLOWRC_MPI_COMMAND_MACHINE2 << "\"" << "  // MACHINE2" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE2=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE2 << "\"" << "  // MACHINE2" << endl;

    aflowrc << " " << endl;
    aflowrc << "// ****************************************************************************************************" << endl;

    //   XHOST.DEBUG=true;
    //[CO20190808 - issue this ONLY if it was written, should fix www-data]cerr << "WARNING: aflowrc::write_default: WRITING default " << XHOST.aflowrc_filename << endl;
    if (aurostd::stringstream2file(aflowrc, XHOST.aflowrc_filename) && aurostd::FileExist(XHOST.aflowrc_filename)) {
      if (!XHOST.vflag_control.flag("WWW")) { // CO20200404 - new web flag
        message << "WRITING default " << XHOST.aflowrc_filename;
        pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, std::cerr, _LOGGER_MESSAGE_); // CO20200404 - LEAVE std::cerr here, FR needs this for web
      }
    }
    if (LDEBUG) {
      oss << __AFLOW_FUNC__ << " END" << endl;
    }
    return true;
  }
} // namespace aflowrc

// ***************************************************************************
// aflowrc::print_aflowrc
// ***************************************************************************
namespace aflowrc {
  bool print_aflowrc(std::ostream& oss, bool AFLOWRC_VERBOSE) {
    const bool LDEBUG = (false || XHOST.DEBUG || AFLOWRC_VERBOSE);
    if (LDEBUG) {
      oss << "aflowrc::print_aflowrc: BEGIN" << endl;
    }
    if (LDEBUG) {
      oss << "aflowrc::print_aflowrc: XHOST.home=" << XHOST.home << endl;
    }
    if (LDEBUG) {
      oss << "aflowrc::print_aflowrc: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
    }

    // HE20220218 START
    if (LDEBUG) {
      oss << "// DEFAULTS ENTRY LOADER" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_ENTRY_LOADER_AFLUX_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_ENTRY_LOADER_AFLUX_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_ENTRY_LOADER_RESTAPI_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_ENTRY_LOADER_RESTAPI_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_ENTRY_LOADER_FS_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH << "\"" << endl;
    }
    // HE20220218 STOP

    // ME20191001 START
    if (LDEBUG) {
      oss << "// DEFAULT AFLOW DATABASE" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_AFLOW_DB_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_AFLOW_DB_STATS_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_AFLOW_DB_DATA_PATH=\"" << AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_AFLOW_DB_LOCK_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "DEFAULT_AFLOW_DB_STALE_THRESHOLD=" << AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD << endl;
    }
    // ME20191001 STOP
    if (LDEBUG) {
      oss << "// DEFAULT DEFINITIONS" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_KZIP_BIN")=")" << DEFAULT_KZIP_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_KZIP_EXT")=")" << DEFAULT_KZIP_EXT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_TMPFS_DIRECTORIES")=")" << DEFAULT_TMPFS_DIRECTORIES << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// FILENAMES FOR AFLOW.ORG ANALYSIS" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_AFLOWLIB_ENTRY_OUT")=")" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_AFLOWLIB_ENTRY_JSON")=")" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_ORIG_OUT")=")" << DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_RELAX_OUT")=")" << DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_BANDS_OUT")=")" << DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_ORIG_OUT")=")" << DEFAULT_FILE_DATA_ORIG_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_RELAX_OUT")=")" << DEFAULT_FILE_DATA_RELAX_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_BANDS_OUT")=")" << DEFAULT_FILE_DATA_BANDS_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_ORIG_JSON")=")" << DEFAULT_FILE_EDATA_ORIG_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_RELAX_JSON")=")" << DEFAULT_FILE_EDATA_RELAX_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_BANDS_JSON")=")" << DEFAULT_FILE_EDATA_BANDS_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_ORIG_JSON")=")" << DEFAULT_FILE_DATA_ORIG_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_RELAX_JSON")=")" << DEFAULT_FILE_DATA_RELAX_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_BANDS_JSON")=")" << DEFAULT_FILE_DATA_BANDS_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_TIME_OUT")=")" << DEFAULT_FILE_TIME_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_SPACEGROUP1_OUT")=")" << DEFAULT_FILE_SPACEGROUP1_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_SPACEGROUP2_OUT")=")" << DEFAULT_FILE_SPACEGROUP2_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_VOLDISTPARAMS_OUT")=")" << DEFAULT_FILE_VOLDISTPARAMS_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_FILE_VOLDISTEVOLUTION_OUT")=")" << DEFAULT_FILE_VOLDISTEVOLUTION_OUT << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// FILENAMES FOR AFLOW OPERATION" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT")=")" << DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PRESCRIPT_OUT")=")" << DEFAULT_AFLOW_PRESCRIPT_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PRESCRIPT_COMMAND")=")" << DEFAULT_AFLOW_PRESCRIPT_COMMAND << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_POSTSCRIPT_OUT")=")" << DEFAULT_AFLOW_POSTSCRIPT_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_POSTSCRIPT_COMMAND")=")" << DEFAULT_AFLOW_POSTSCRIPT_COMMAND << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_OUT")=")" << DEFAULT_AFLOW_PGROUP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_JSON")=")" << DEFAULT_AFLOW_PGROUP_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_XTAL_OUT")=")" << DEFAULT_AFLOW_PGROUP_XTAL_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_XTAL_JSON")=")" << DEFAULT_AFLOW_PGROUP_XTAL_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT")=")" << DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT << "\"" << endl; // DX20200129
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON")=")" << DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON << "\"" << endl; // DX20200129
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_OUT")=")" << DEFAULT_AFLOW_PGROUPK_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_JSON")=")" << DEFAULT_AFLOW_PGROUPK_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_XTAL_OUT")=")" << DEFAULT_AFLOW_PGROUPK_XTAL_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_XTAL_JSON")=")" << DEFAULT_AFLOW_PGROUPK_XTAL_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FGROUP_OUT")=")" << DEFAULT_AFLOW_FGROUP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FGROUP_JSON")=")" << DEFAULT_AFLOW_FGROUP_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SGROUP_OUT")=")" << DEFAULT_AFLOW_SGROUP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SGROUP_JSON")=")" << DEFAULT_AFLOW_SGROUP_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_AGROUP_OUT")=")" << DEFAULT_AFLOW_AGROUP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_AGROUP_JSON")=")" << DEFAULT_AFLOW_AGROUP_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IATOMS_OUT")=")" << DEFAULT_AFLOW_IATOMS_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IATOMS_JSON")=")" << DEFAULT_AFLOW_IATOMS_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_ICAGES_OUT")=")" << DEFAULT_AFLOW_ICAGES_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SURFACE_OUT")=")" << DEFAULT_AFLOW_SURFACE_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_QMVASP_OUT")=")" << DEFAULT_AFLOW_QMVASP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_ERVASP_OUT")=")" << DEFAULT_AFLOW_ERVASP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IMMISCIBILITY_OUT")=")" << DEFAULT_AFLOW_IMMISCIBILITY_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_MEMORY_OUT")=")" << DEFAULT_AFLOW_MEMORY_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_INPUT_OUT")=")" << DEFAULT_AFLOW_FROZSL_INPUT_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_POSCAR_OUT")=")" << DEFAULT_AFLOW_FROZSL_POSCAR_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_MODES_OUT")=")" << DEFAULT_AFLOW_FROZSL_MODES_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_EIGEN_OUT")=")" << DEFAULT_AFLOW_FROZSL_EIGEN_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_END_OUT")=")" << DEFAULT_AFLOW_END_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DIELECTRIC_FILE")=")" << DEFAULT_AFLOW_DIELECTRIC_FILE << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULT GENERIC MPI" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_START_DEFAULT")=")" << MPI_START_DEFAULT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_STOP_DEFAULT")=")" << MPI_STOP_DEFAULT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DEFAULT")=")" << MPI_COMMAND_DEFAULT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_DEFAULT\")=" << MPI_NCPUS_DEFAULT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MAX\")=" << MPI_NCPUS_MAX << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS BINARY" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_BIN")=")" << DEFAULT_VASP_GAMMA_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_MPI_BIN")=")" << DEFAULT_VASP_GAMMA_MPI_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_BIN")=")" << DEFAULT_VASP_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_MPI_BIN")=")" << DEFAULT_VASP_MPI_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP5_BIN")=")" << DEFAULT_VASP5_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP5_MPI_BIN")=")" << DEFAULT_VASP5_MPI_BIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AIMS_BIN")=")" << DEFAULT_AIMS_BIN << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS POTCARS" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIRECTORIES")=")" << DEFAULT_VASP_POTCAR_DIRECTORIES << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE")=")" << DEFAULT_VASP_POTCAR_DATE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_SUFFIX")=")" << DEFAULT_VASP_POTCAR_SUFFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_LDA")=")" << DEFAULT_VASP_POTCAR_DATE_POT_LDA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_GGA")=")" << DEFAULT_VASP_POTCAR_DATE_POT_GGA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_LDA")=")" << DEFAULT_VASP_POTCAR_DIR_POT_LDA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_GGA")=")" << DEFAULT_VASP_POTCAR_DIR_POT_GGA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_PBE")=")" << DEFAULT_VASP_POTCAR_DIR_POT_PBE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA")=")" << DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA")=")" << DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE")=")" << DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN")=")" << DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("VASP_PSEUDOPOTENTIAL_DIRECTORY_POTPAW_PBE_KIN_")=")" << DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULT KPOINTS/DOS" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_BANDS_GRID\")=" << DEFAULT_BANDS_GRID << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_BANDS_LATTICE")=")" << DEFAULT_BANDS_LATTICE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_KSCHEME")=")" << DEFAULT_KSCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA\")=" << DEFAULT_KPPRA << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_STATIC_KSCHEME")=")" << DEFAULT_STATIC_KSCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_STATIC\")=" << DEFAULT_KPPRA_STATIC << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_ICSD\")=" << DEFAULT_KPPRA_ICSD << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_BANDS_GRID\")=" << DEFAULT_UNARY_BANDS_GRID << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA\")=" << DEFAULT_UNARY_KPPRA << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA_STATIC\")=" << DEFAULT_UNARY_KPPRA_STATIC << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA_DIELECTRIC\")=" << DEFAULT_UNARY_KPPRA_DIELECTRIC << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_PHONONS_KSCHEME")=")" << DEFAULT_PHONONS_KSCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_PHONONS_KPPRA\")=" << DEFAULT_PHONONS_KPPRA << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_DIELECTRIC_KSCHEME")=")" << DEFAULT_DIELECTRIC_KSCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_DIELECTRIC\")=" << DEFAULT_KPPRA_DIELECTRIC << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_EMIN\")=" << DEFAULT_DOS_EMIN << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_EMAX\")=" << DEFAULT_DOS_EMAX << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_SCALE\")=" << DEFAULT_DOS_SCALE << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULT PRECISION" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_LOW\")=" << DEFAULT_VASP_PREC_ENMAX_LOW << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_MEDIUM\")=" << DEFAULT_VASP_PREC_ENMAX_MEDIUM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_NORMAL\")=" << DEFAULT_VASP_PREC_ENMAX_NORMAL << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_HIGH\")=" << DEFAULT_VASP_PREC_ENMAX_HIGH << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_ACCURATE\")=" << DEFAULT_VASP_PREC_ENMAX_ACCURATE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_ENMAX_MINIMUM\")=" << DEFAULT_VASP_ENMAX_MINIMUM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_SPIN_REMOVE_CUTOFF\")=" << DEFAULT_VASP_SPIN_REMOVE_CUTOFF << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_POTIM\")=" << DEFAULT_VASP_PREC_POTIM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_EDIFFG\")=" << DEFAULT_VASP_PREC_EDIFFG << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS OPTIONS " << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_OUT")=")" << DEFAULT_VASP_OUT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_INCAR")=")" << DEFAULT_VASP_EXTERNAL_INCAR << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POSCAR")=")" << DEFAULT_VASP_EXTERNAL_POSCAR << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POTCAR")=")" << DEFAULT_VASP_EXTERNAL_POTCAR << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_KPOINT")=")" << DEFAULT_VASP_EXTERNAL_KPOINTS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_CONTROL")=")" << DEFAULT_AIMS_EXTERNAL_CONTROL << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_GEOM")=")" << DEFAULT_AIMS_EXTERNAL_GEOM << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE")=")" << DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT\")=" << DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_FORCE_OPTION_SIGMA")=")" << DEFAULT_VASP_FORCE_OPTION_SIGMA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC")=")" << DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS")=")" << DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_NELM\")=" << DEFAULT_VASP_FORCE_OPTION_NELM << endl; // CO20200624
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_NELM_STATIC\")=" << DEFAULT_VASP_FORCE_OPTION_NELM_STATIC << endl; // CO20200624
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"MAX_VASP_NELM\")=" << MAX_VASP_NELM << endl; // CO20200624
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME")=")" << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SYM\")=" << DEFAULT_VASP_FORCE_OPTION_SYM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_BADER\")=" << DEFAULT_VASP_FORCE_OPTION_BADER << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_BADER_STATIC\")=" << DEFAULT_VASP_FORCE_OPTION_BADER_STATIC << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_ELF\")=" << DEFAULT_VASP_FORCE_OPTION_ELF << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM\")=" << DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_WAVECAR\")=" << DEFAULT_VASP_FORCE_OPTION_WAVECAR << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_CHGCAR\")=" << DEFAULT_VASP_FORCE_OPTION_CHGCAR << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_LSCOUPLING\")=" << DEFAULT_VASP_FORCE_OPTION_LSCOUPLING << endl;
    }

    if (LDEBUG) {
      oss << "// AFLOW_LIBRARY AFLOW_PROJECT" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_LIBRARY_DIRECTORIES")=")" << DEFAULT_AFLOW_LIBRARY_DIRECTORIES << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PROJECTS_DIRECTORIES")=")" << DEFAULT_AFLOW_PROJECTS_DIRECTORIES << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AFLOWDATA_WEB_DIRECTORY")=")" << DEFAULT_AFLOWDATA_WEB_DIRECTORY << "\"" << endl; // CO+ME20200731
    }

    if (LDEBUG) {
      oss << "// DEFAULT PLATON/FINDSYM" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_PLATON_P_EQUAL\")=" << DEFAULT_PLATON_P_EQUAL << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_PLATON_P_EXACT\")=" << DEFAULT_PLATON_P_EXACT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_ANG\")=" << DEFAULT_PLATON_P_ANG << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D1\")=" << DEFAULT_PLATON_P_D1 << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D2\")=" << DEFAULT_PLATON_P_D2 << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D3\")=" << DEFAULT_PLATON_P_D3 << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_FINDSYM_TOL\")=" << DEFAULT_FINDSYM_TOL << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULT GNUPLOT" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT")=")" << DEFAULT_GNUPLOT_EPS_FONT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD")=")" << DEFAULT_GNUPLOT_EPS_FONT_BOLD << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_ITALICS")=")" << DEFAULT_GNUPLOT_EPS_FONT_ITALICS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS")=")" << DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT")=")" << DEFAULT_GNUPLOT_PNG_FONT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD")=")" << DEFAULT_GNUPLOT_PNG_FONT_BOLD << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_ITALICS")=")" << DEFAULT_GNUPLOT_PNG_FONT_ITALICS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS")=")" << DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT")=")" << DEFAULT_GNUPLOT_GREEK_FONT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_BOLD")=")" << DEFAULT_GNUPLOT_GREEK_FONT_BOLD << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_ITALICS")=")" << DEFAULT_GNUPLOT_GREEK_FONT_ITALICS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS")=")" << DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULT NHULL" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_NHULL_ALLOWED_DFT_TYPES\")=" << DEFAULT_NHULL_ALLOWED_DFT_TYPES << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_NHULL_ALLOW_ALL_FORMATION_ENERGIES\")=" << DEFAULT_NHULL_ALLOW_ALL_FORMATION_ENERGIES << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_NHULL_COUNT_THRESHOLD_BINARIES\")=" << DEFAULT_NHULL_COUNT_THRESHOLD_BINARIES << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_NHULL_PERFORM_OUTLIER_ANALYSIS\")=" << DEFAULT_NHULL_PERFORM_OUTLIER_ANALYSIS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES\")=" << DEFAULT_NHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_NHULL_OUTLIER_MULTIPLIER\")=" << DEFAULT_NHULL_OUTLIER_MULTIPLIER << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_NHULL_IGNORE_KNOWN_ILL_CONVERGED\")=" << DEFAULT_NHULL_IGNORE_KNOWN_ILL_CONVERGED << endl;
    }
    if (LDEBUG) {
      oss << "// DEFAULTS GFA" << endl; // CO20190628
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF\")=" << DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF << endl; // CO20190628
    }

    if (LDEBUG) {
      oss << "// DEFAULTS ARUN" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("ARUN_DIRECTORY_PREFIX")=")" << ARUN_DIRECTORY_PREFIX << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS POCC" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_POCC_STRUCTURE_GENERATION_ALGO")=")" << DEFAULT_POCC_STRUCTURE_GENERATION_ALGO << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_POCC_TEMPERATURE_STRING")=")" << DEFAULT_POCC_TEMPERATURE_STRING << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_EXCLUDE_UNSTABLE\")=" << DEFAULT_POCC_EXCLUDE_UNSTABLE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_SITE_TOL\")=" << DEFAULT_POCC_SITE_TOL << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_STOICH_TOL\")=" << DEFAULT_POCC_STOICH_TOL << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_UFF_BONDING_DISTANCE\")=" << DEFAULT_UFF_BONDING_DISTANCE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_UFF_ENERGY_TOLERANCE\")=" << DEFAULT_UFF_ENERGY_TOLERANCE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_UFF_CLUSTER_RADIUS\")=" << DEFAULT_UFF_CLUSTER_RADIUS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_RDF_RMAX\")=" << DEFAULT_POCC_RDF_RMAX << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_RDF_NBINS\")=" << DEFAULT_POCC_RDF_NBINS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON\")=" << DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS\")=" << DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_FILE_PREFIX")=")" << POCC_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_POCC_JSON")=")" << DEFAULT_POCC_JSON << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_APL_OUT_FILE")=")" << POCC_APL_OUT_FILE << "\"" << endl; // ME20210927
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_ALL_SUPERCELLS_FILE")=")" << POCC_ALL_SUPERCELLS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_UNIQUE_SUPERCELLS_FILE")=")" << POCC_UNIQUE_SUPERCELLS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_ALL_HNF_MATRICES_FILE")=")" << POCC_ALL_HNF_MATRICES_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_ALL_SITE_CONFIGURATIONS_FILE")=")" << POCC_ALL_SITE_CONFIGURATIONS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_DOSCAR_FILE")=")" << POCC_DOSCAR_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_PHDOSCAR_FILE")=")" << POCC_PHDOSCAR_FILE << "\"" << endl; // ME20210927
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("POCC_ANIONS_LIST")=")" << POCC_ANIONS_LIST << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS APL" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PREC")=")" << DEFAULT_APL_PREC << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_ENGINE")=")" << DEFAULT_APL_ENGINE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_HIBERNATE\")=" << DEFAULT_APL_HIBERNATE << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_MINSHELL\")=" << DEFAULT_APL_MINSHELL << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_MINATOMS\")=" << DEFAULT_APL_MINATOMS << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_POLAR\")=" << DEFAULT_APL_POLAR << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DMAG\")=" << DEFAULT_APL_DMAG << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DXYZONLY\")=" << DEFAULT_APL_DXYZONLY << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DSYMMETRIZE\")=" << DEFAULT_APL_DSYMMETRIZE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DINEQUIV_ONLY\")=" << DEFAULT_APL_DINEQUIV_ONLY << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_DPM")=")" << DEFAULT_APL_DPM << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_RELAX\")=" << DEFAULT_APL_RELAX << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_RELAX_COMMENSURATE\")=" << DEFAULT_APL_RELAX_COMMENSURATE << endl; // ME20200427
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ZEROSTATE\")=" << DEFAULT_APL_ZEROSTATE << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ZEROSTATE_CHGCAR\")=" << DEFAULT_APL_ZEROSTATE_CHGCAR << endl; // ME20191029
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_USE_LEPSILON\")=" << DEFAULT_APL_USE_LEPSILON << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_FREQFORMAT")=")" << DEFAULT_APL_FREQFORMAT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DC\")=" << DEFAULT_APL_DC << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_DCPATH")=")" << DEFAULT_APL_DCPATH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DCPOINTS\")=" << DEFAULT_APL_DCPOINTS << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOS\")=" << DEFAULT_APL_DOS << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_DOSMETHOD")=")" << DEFAULT_APL_DOSMETHOD << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_DOSMESH")=")" << DEFAULT_APL_DOSMESH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOSPOINTS\")=" << DEFAULT_APL_DOSPOINTS << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOSSMEAR\")=" << DEFAULT_APL_DOSSMEAR << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOS_PROJECT\")=" << DEFAULT_APL_DOS_PROJECT << endl; // ME20200213
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_TP\")=" << DEFAULT_APL_TP << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DISPLACEMENTS\")=" << DEFAULT_APL_DISPLACEMENTS << endl; // ME20200421
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_TPT")=")" << DEFAULT_APL_TPT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_GVEL\")=" << DEFAULT_APL_GVEL << endl; // ME20200517
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_FILE_PREFIX")=")" << DEFAULT_APL_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_OUT_FILE")=")" << DEFAULT_APL_OUT_FILE << "\"" << endl; // ME20210927
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PDIS_FILE")=")" << DEFAULT_APL_PDIS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PDOS_FILE")=")" << DEFAULT_APL_PDOS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_THERMO_FILE")=")" << DEFAULT_APL_THERMO_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_THERMO_JSON")=")" << DEFAULT_APL_THERMO_JSON << "\"" << endl; // ME20211019
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_DYNMAT_FILE")=")" << DEFAULT_APL_DYNMAT_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_HARMIFC_FILE")=")" << DEFAULT_APL_HARMIFC_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_POLAR_FILE")=")" << DEFAULT_APL_POLAR_FILE << "\"" << endl; // ME20200415
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_HSKPTS_FILE")=")" << DEFAULT_APL_HSKPTS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_MSQRDISP_FILE")=")" << DEFAULT_APL_MSQRDISP_FILE << "\"" << endl; // ME20200329
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_GVEL_FILE")=")" << DEFAULT_APL_GVEL_FILE << "\"" << endl; // ME20200517
    }
    // ME20190614 BEGIN
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PHDOSCAR_FILE")=")" << DEFAULT_APL_PHDOSCAR_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PHPOSCAR_FILE")=")" << DEFAULT_APL_PHPOSCAR_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PHKPOINTS_FILE")=")" << DEFAULT_APL_PHKPOINTS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_PHEIGENVAL_FILE")=")" << DEFAULT_APL_PHEIGENVAL_FILE << "\"" << endl;
    }
    // ME20190614 END
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_STATE_FILE")=")" << DEFAULT_APL_STATE_FILE << "\"" << endl; // ME20200224
    }
    // ME20200329 BEGIN
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_APL_ADISP_SCENE_FORMAT")=")" << DEFAULT_APL_ADISP_SCENE_FORMAT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_AMPLITUDE\")=" << DEFAULT_APL_ADISP_AMPLITUDE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_NSTEPS\")=" << DEFAULT_APL_ADISP_NSTEPS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_NPERIODS\")=" << DEFAULT_APL_ADISP_NPERIODS << endl;
    }
    // ME20200329 END

    if (LDEBUG) {
      oss << "// DEFAULTS QHA" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_MODE")=")" << DEFAULT_QHA_MODE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_EOS\")=" << DEFAULT_QHA_EOS << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_EOS_DISTORTION_RANGE")=")" << DEFAULT_QHA_EOS_DISTORTION_RANGE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_EOS_MODEL")=")" << DEFAULT_QHA_EOS_MODEL << "\"" << endl; // AS20200818
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_DISTORTION\")=" << DEFAULT_QHA_GP_DISTORTION << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_TAYLOR_EXPANSION_ORDER\")=" << DEFAULT_QHA_TAYLOR_EXPANSION_ORDER << endl; // AS20200602
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_INCLUDE_ELEC_CONTRIB\")=" << DEFAULT_QHA_INCLUDE_ELEC_CONTRIB << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_SOMMERFELD_EXPANSION\")=" << DEFAULT_QHA_SOMMERFELD_EXPANSION << endl; // AS20200528
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDIS_T")=")" << DEFAULT_QHA_PDIS_T << "\"" << endl;
    }
    // AS20200508 BEGIN
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_FINITE_DIFF\")=" << DEFAULT_QHA_GP_FINITE_DIFF << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_IGNORE_IMAGINARY\")=" << DEFAULT_QHA_IGNORE_IMAGINARY << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_RELAX_IONS_CELL\")=" << DEFAULT_QHA_IGNORE_IMAGINARY << endl; // AS20201123
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_FILE_PREFIX")=")" << DEFAULT_QHA_FILE_PREFIX << "\"" << endl;
    }
    // AS20200709 BEGIN
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA3P_FILE_PREFIX")=")" << DEFAULT_QHA3P_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHANP_FILE_PREFIX")=")" << DEFAULT_QHANP_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_SCQHA_FILE_PREFIX")=")" << DEFAULT_SCQHA_FILE_PREFIX << "\"" << endl;
    }
    // AS20200709 END
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_PATH_FILE")=")" << DEFAULT_QHA_GP_PATH_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_MESH_FILE")=")" << DEFAULT_QHA_GP_MESH_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_AVG_FILE")=")" << DEFAULT_QHA_GP_AVG_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_THERMO_FILE")=")" << DEFAULT_QHA_THERMO_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_FREQS_FILE")=")" << DEFAULT_QHA_FREQS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_FVT_FILE")=")" << DEFAULT_QHA_FVT_FILE << "\"" << endl;
    }
    // AS20200508 END
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_COEFF_FILE")=")" << DEFAULT_QHA_COEFF_FILE << "\"" << endl; // AS20210517
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_IMAG_FILE")=")" << DEFAULT_QHA_IMAG_FILE << "\"" << endl; // AS20210517
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDIS_FILE")=")" << DEFAULT_QHA_PDIS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDOS_FILE")=")" << DEFAULT_QHA_PDOS_FILE << "\"" << endl; // AS20201201
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_QHA_KPOINTS_FILE")=")" << DEFAULT_QHA_KPOINTS_FILE << "\"" << endl; // AS20201112
    }
    // AS20210914 BEGIN
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_POCC_QHA_THERMO_FILE")=")" << DEFAULT_POCC_QHA_THERMO_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_POCC_QHA_AVGTHERMO_FILE")=")" << DEFAULT_POCC_QHA_AVGTHERMO_FILE << "\"" << endl;
    }
    // AS20210914 END
    if (LDEBUG) {
      oss << "// DEFAULTS AAPL" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_BTE")=")" << DEFAULT_AAPL_BTE << "\"" << endl;
    }
    //[ME20181226]if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_BZMETHOD\")=\"" << DEFAULT_AAPL_BZMETHOD << "\"" << endl;
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_FOURTH_ORDER\")=" << DEFAULT_AAPL_FOURTH_ORDER << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_CUT_RAD")=")" << DEFAULT_AAPL_CUT_RAD << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_CUT_SHELL")=")" << DEFAULT_AAPL_CUT_SHELL << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_THERMALGRID")=")" << DEFAULT_AAPL_THERMALGRID << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_TCT")=")" << DEFAULT_AAPL_TCT << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_SUMRULE\")=" << DEFAULT_AAPL_SUMRULE << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_SUMRULE_MAX_ITER\")=" << DEFAULT_AAPL_SUMRULE_MAX_ITER << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_MIXING_COEFFICIENT\")=" << DEFAULT_AAPL_MIXING_COEFFICIENT << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_ISOTOPE\")=" << DEFAULT_AAPL_ISOTOPE << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_BOUNDARY\")=" << DEFAULT_AAPL_BOUNDARY << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_CUMULATIVEK\")=" << DEFAULT_AAPL_CUMULATIVEK << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_NANO_SIZE\")=" << DEFAULT_AAPL_NANO_SIZE << endl; // ME20190112
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_FILE_PREFIX")=")" << DEFAULT_AAPL_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_IRRQPTS_FILE")=")" << DEFAULT_AAPL_IRRQPTS_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_GVEL_FILE")=")" << DEFAULT_AAPL_GVEL_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_PS_FILE")=")" << DEFAULT_AAPL_PS_FILE << "\"" << endl; // ME20191104
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_GRUENEISEN_FILE")=")" << DEFAULT_AAPL_GRUENEISEN_FILE << "\"" << endl; // ME20191104
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_FILE")=")" << DEFAULT_AAPL_RATES_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_3RD_FILE")=")" << DEFAULT_AAPL_RATES_3RD_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_4TH_FILE")=")" << DEFAULT_AAPL_RATES_4TH_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_ISOTOPE_FILE")=")" << DEFAULT_AAPL_ISOTOPE_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_BOUNDARY_FILE")=")" << DEFAULT_AAPL_BOUNDARY_FILE << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AAPL_TCOND_FILE")=")" << DEFAULT_AAPL_TCOND_FILE << "\"" << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS AEL" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_STRAIN_SYMMETRY\")=" << AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NNORMAL_STRAINS\")=" << AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NSHEAR_STRAINS\")=" << AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NORMAL_STRAIN_STEP\")=" << AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_SHEAR_STRAIN_STEP\")=" << AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_ORIGIN_STRAIN_CALC\")=" << AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_ORIGIN_STRAIN_FIT\")=" << AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_RELAXED_STRUCT_FIT\")=" << AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NEG_STRAINS\")=" << AFLOWRC_DEFAULT_AEL_NEG_STRAINS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NIND_STRAIN_DIRS\")=" << AFLOWRC_DEFAULT_AEL_NIND_STRAIN_DIRS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_VASPSYM\")=" << AFLOWRC_DEFAULT_AEL_VASPSYM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_PRECACC_ALGONORM\")=" << AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_VASPRUNXML_STRESS\")=" << AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS\")=" << AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_SKIP_ARUNS_MAX\")=" << AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY\")=" << AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_SYMMETRIZE\")=" << AFLOWRC_DEFAULT_AEL_SYMMETRIZE << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AEL_FILE_PREFIX")=")" << AFLOWRC_DEFAULT_AEL_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_WRITE_FULL_RESULTS\")=" << AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_DIRNAME_ARUN\")=" << AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN << endl;
    }

    if (LDEBUG) {
      oss << "// DEFAULTS AGL" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_AEL_POISSON_RATIO\")=" << AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_NSTRUCTURES\")=" << AFLOWRC_DEFAULT_AGL_NSTRUCTURES << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_STRAIN_STEP\")=" << AFLOWRC_DEFAULT_AGL_STRAIN_STEP << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS\")=" << AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_SKIP_ARUNS_MAX\")=" << AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_NTEMPERATURE\")=" << AFLOWRC_DEFAULT_AGL_NTEMPERATURE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_STEMPERATURE\")=" << AFLOWRC_DEFAULT_AGL_STEMPERATURE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_NPRESSURE\")=" << AFLOWRC_DEFAULT_AGL_NPRESSURE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_SPRESSURE\")=" << AFLOWRC_DEFAULT_AGL_SPRESSURE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_POISSON_RATIO\")=" << AFLOWRC_DEFAULT_AGL_POISSON_RATIO << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_IEOS\")=" << AFLOWRC_DEFAULT_AGL_IEOS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_IDEBYE\")=" << AFLOWRC_DEFAULT_AGL_IDEBYE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_FIT_TYPE\")=" << AFLOWRC_DEFAULT_AGL_FIT_TYPE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_CHECK_EV_CONCAVITY\")=" << AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_CHECK_EV_MIN\")=" << AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_HUGONIOT_CALC\")=" << AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_HUGONIOT_EXTRAPOLATE\")=" << AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE\")=" << AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_AGL_FILE_PREFIX")=")" << AFLOWRC_DEFAULT_AGL_FILE_PREFIX << "\"" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_WRITE_FULL_RESULTS\")=" << AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_DIRNAME_ARUN\")=" << AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_WRITE_GIBBS_INPUT\")=" << AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_PLOT_RESULTS\")=" << AFLOWRC_DEFAULT_AGL_PLOT_RESULTS << endl;
    }

    // SD20220323 - QCA START
    if (LDEBUG) {
      oss << "// DEFAULTS QCA " << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_MIN_SLEEP_SECONDS\")=" << AFLOWRC_DEFAULT_QCA_MIN_SLEEP_SECONDS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_MAX_NUM_ATOMS\")=" << AFLOWRC_DEFAULT_QCA_MAX_NUM_ATOMS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS\")=" << AFLOWRC_DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_CV_CUTOFF\")=" << AFLOWRC_DEFAULT_QCA_CV_CUTOFF << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_CONC_NPTS\")=" << AFLOWRC_DEFAULT_QCA_CONC_NPTS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_NPTS\")=" << AFLOWRC_DEFAULT_QCA_TEMP_NPTS << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_MIN\")=" << AFLOWRC_DEFAULT_QCA_TEMP_MIN << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_MAX\")=" << AFLOWRC_DEFAULT_QCA_TEMP_MAX << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_MIN_LIMIT\")=" << AFLOWRC_DEFAULT_QCA_TEMP_MIN_LIMIT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_PRINT\")=" << AFLOWRC_DEFAULT_QCA_PRINT << endl;
    }
    // SD20220323 - QCA END

    // RF20200413 START
    if (LDEBUG) {
      oss << "// DEFAULTS CCE" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_OX_METHOD")=")" << DEFAULT_CCE_OX_METHOD << "\"" << "               // 1 - ELECTRONEGATIVITY_ALLEN, 2 - BADER" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION")=")" << DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION << "\""
          << "               // 0.4 Ang tolerance between shortest and longest bonds when testing for multi-anion compound" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_OX_TOL")=")" << DEFAULT_CCE_OX_TOL << "\"" << "               // sum of oxidation states might not be exactly zero due to numerics" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_PEROX_CUTOFF")=")" << DEFAULT_CCE_PEROX_CUTOFF << "\""
          << "               // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_SUPEROX_CUTOFF")=")" << DEFAULT_CCE_SUPEROX_CUTOFF << "\""
          << "               // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF")=")" << DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF << "\"" << "               // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF")=")" << DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF << "\"" << "               // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    }
    // RF20200413 END

    // DX20200708 - START
    if (LDEBUG) {
      oss << "// DEFAULTS XTALFINDER" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_MISFIT_MATCH\")=" << DEFAULT_XTALFINDER_MISFIT_MATCH << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_MISFIT_FAMILY\")=" << DEFAULT_XTALFINDER_MISFIT_FAMILY << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_SUPERCELL_METHOD\")=" << DEFAULT_XTALFINDER_SUPERCELL_METHOD << endl; // DX20201223
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING\")=" << DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_MATERIAL\")=" << DEFAULT_XTALFINDER_FILE_MATERIAL << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_STRUCTURE\"=" << DEFAULT_XTALFINDER_FILE_STRUCTURE << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_DUPLICATE\")=" << DEFAULT_XTALFINDER_FILE_DUPLICATE << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE\")=" << DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE\")=" << DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE\")=" << DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE << endl; // DX20201118
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE\")=" << DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE << endl; // DX20201118
    }
    // DX20200708 - END

    // DX20200720 - START
    if (LDEBUG) {
      oss << "// DEFAULTS ANRL" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL\")=" << DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL << endl;
    }
    // DX20200708 - END

    if (LDEBUG) {
      oss << "// DEFAULT CORE" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_BEEP\")=" << AFLOW_CORE_TEMPERATURE_BEEP << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_HALT\")=" << AFLOW_CORE_TEMPERATURE_HALT << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_REFRESH\")=" << AFLOW_CORE_TEMPERATURE_REFRESH << endl;
    }

    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"SECONDS_SLEEP_VASP_COMPLETION\")=" << SECONDS_SLEEP_VASP_COMPLETION << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"SECONDS_SLEEP_VASP_MONITOR\")=" << SECONDS_SLEEP_VASP_MONITOR << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"SECONDS_STALE_OUTCAR\")=" << SECONDS_STALE_OUTCAR << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<unsigned long long int>(\"BYTES_MAX_VASP_OUT\")=" << BYTES_MAX_VASP_OUT << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"MEMORY_MAX_USAGE_RAM\")=" << MEMORY_MAX_USAGE_RAM << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedutype<double>(\"MEMORY_MAX_USAGE_SWAP\")=" << MEMORY_MAX_USAGE_SWAP << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"FILE_VASP_MONITOR\")=" << FILE_VASP_MONITOR << endl; // CO20201111
    }
    if (LDEBUG) {
      oss << "XHOST.adefault.getattachedscheme(\"INTEL_COMPILER_PATHS\")=" << INTEL_COMPILER_PATHS << endl; // CO20201111
    }

    if (LDEBUG) {
      oss << "// DEFAULT MACHINE DEPENDENT MPI" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_MPICH")=")" << MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_MPICH")=")" << MPI_COMMAND_DUKE_BETA_MPICH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_MPICH")=")" << MPI_BINARY_DIR_DUKE_BETA_MPICH << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_OPENMPI")=")" << MPI_OPTIONS_DUKE_BETA_OPENMPI << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_OPENMPI")=")" << MPI_COMMAND_DUKE_BETA_OPENMPI << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_OPENMPI")=")" << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_MATERIALS")=")" << MPI_OPTIONS_DUKE_MATERIALS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_MATERIALS")=")" << MPI_COMMAND_DUKE_MATERIALS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_MATERIALS")=")" << MPI_BINARY_DIR_DUKE_MATERIALS << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_AFLOWLIB")=")" << MPI_OPTIONS_DUKE_AFLOWLIB << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_AFLOWLIB")=")" << MPI_COMMAND_DUKE_AFLOWLIB << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_AFLOWLIB")=")" << MPI_BINARY_DIR_DUKE_AFLOWLIB << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_HABANA")=")" << MPI_OPTIONS_DUKE_HABANA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_HABANA")=")" << MPI_COMMAND_DUKE_HABANA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_HABANA")=")" << MPI_BINARY_DIR_DUKE_HABANA << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QRATS_MPICH")=")" << MPI_OPTIONS_DUKE_QRATS_MPICH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QRATS_MPICH")=")" << MPI_COMMAND_DUKE_QRATS_MPICH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QRATS_MPICH")=")" << MPI_BINARY_DIR_DUKE_QRATS_MPICH << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QFLOW_OPENMPI")=")" << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QFLOW_OPENMPI")=")" << MPI_COMMAND_DUKE_QFLOW_OPENMPI << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI")=")" << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << "\"" << endl;
    }

    // CO20201220 X START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_X")=")" << MPI_OPTIONS_DUKE_X_X << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_X")=")" << MPI_COMMAND_DUKE_X_X << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_X")=")" << MPI_BINARY_DIR_DUKE_X_X << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_CRAY")=")" << MPI_OPTIONS_DUKE_X_CRAY << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_CRAY")=")" << MPI_COMMAND_DUKE_X_CRAY << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_CRAY")=")" << MPI_BINARY_DIR_DUKE_X_CRAY << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_OLDCRAY")=")" << MPI_OPTIONS_DUKE_X_OLDCRAY << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_OLDCRAY")=")" << MPI_COMMAND_DUKE_X_OLDCRAY << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_OLDCRAY")=")" << MPI_BINARY_DIR_DUKE_X_OLDCRAY << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_SMB")=")" << MPI_OPTIONS_DUKE_X_SMB << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_SMB")=")" << MPI_COMMAND_DUKE_X_SMB << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_SMB")=")" << MPI_BINARY_DIR_DUKE_X_SMB << "\"" << endl;
    }
    // CO20201220 X STOP

    // CO20220818 JHU_ROCKFISH START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_JHU_ROCKFISH")=")" << MPI_OPTIONS_JHU_ROCKFISH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_JHU_ROCKFISH")=")" << MPI_COMMAND_JHU_ROCKFISH << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_JHU_ROCKFISH")=")" << MPI_BINARY_DIR_JHU_ROCKFISH << "\"" << endl;
    }
    // CO20220818 JHU_ROCKFISH STOP

    // DX20190509 - MACHINE001 - START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE001")=")" << MPI_OPTIONS_MACHINE001 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE001")=")" << MPI_COMMAND_MACHINE001 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE001")=")" << MPI_BINARY_DIR_MACHINE001 << "\"" << endl;
    }
    // DX20190509 - MACHINE001 - END

    // DX20190509 - MACHINE002 - START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE002")=")" << MPI_OPTIONS_MACHINE002 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE002")=")" << MPI_COMMAND_MACHINE002 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE002")=")" << MPI_BINARY_DIR_MACHINE002 << "\"" << endl;
    }
    // DX20190509 - MACHINE002 - END

    // DX20201005 - MACHINE003 - START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE003")=")" << MPI_OPTIONS_MACHINE003 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE003")=")" << MPI_COMMAND_MACHINE003 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE003")=")" << MPI_BINARY_DIR_MACHINE003 << "\"" << endl;
    }
    // DX20201005 - MACHINE003 - END

    // DX20211011 - MACHINE004 - START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE004")=")" << MPI_OPTIONS_MACHINE004 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE004")=")" << MPI_COMMAND_MACHINE004 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE004")=")" << MPI_BINARY_DIR_MACHINE004 << "\"" << endl;
    }
    // DX2021101 - MACHINE004 - END

    // DX20190107 - CMU EULER - START
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_CMU_EULER")=")" << MPI_OPTIONS_CMU_EULER << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_CMU_EULER")=")" << MPI_COMMAND_CMU_EULER << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_CMU_EULER")=")" << MPI_BINARY_DIR_CMU_EULER << "\"" << endl;
    }
    // DX20190107 - CMU EULER - END

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_EOS")=")" << MPI_OPTIONS_MPCDF_EOS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_EOS")=")" << MPI_COMMAND_MPCDF_EOS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_EOS")=")" << MPI_NCPUS_MPCDF_EOS << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype("MPI_HYPERTHREADING_MPCDF_EOS")=")" << MPI_HYPERTHREADING_MPCDF_EOS << "\"" << "            // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_EOS")=")" << MPI_BINARY_DIR_MPCDF_EOS << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_DRACO")=")" << MPI_OPTIONS_MPCDF_DRACO << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_DRACO")=")" << MPI_COMMAND_MPCDF_DRACO << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_DRACO")=")" << MPI_NCPUS_MPCDF_DRACO << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype("MPI_HYPERTHREADING_MPCDF_DRACO")=")" << MPI_HYPERTHREADING_MPCDF_DRACO << "\"" << "            // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_DRACO")=")" << MPI_BINARY_DIR_MPCDF_DRACO << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_COBRA")=")" << MPI_OPTIONS_MPCDF_COBRA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_COBRA")=")" << MPI_COMMAND_MPCDF_COBRA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_COBRA")=")" << MPI_NCPUS_MPCDF_COBRA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype("MPI_HYPERTHREADING_MPCDF_COBRA")=")" << MPI_HYPERTHREADING_MPCDF_COBRA << "\"" << "            // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_COBRA")=")" << MPI_BINARY_DIR_MPCDF_COBRA << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_HYDRA")=")" << MPI_OPTIONS_MPCDF_HYDRA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_HYDRA")=")" << MPI_COMMAND_MPCDF_HYDRA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_HYDRA")=")" << MPI_NCPUS_MPCDF_HYDRA << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedutype("MPI_HYPERTHREADING_MPCDF_HYDRA")=")" << MPI_HYPERTHREADING_MPCDF_HYDRA << "\"" << "            // false/OFF, IGNORE/NEGLECT, true/ON " << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_HYDRA")=")" << MPI_BINARY_DIR_MPCDF_HYDRA << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_FULTON_MARYLOU")=")" << MPI_OPTIONS_FULTON_MARYLOU << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_FULTON_MARYLOU")=")" << MPI_COMMAND_FULTON_MARYLOU << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_FULTON_MARYLOU")=")" << MPI_BINARY_DIR_FULTON_MARYLOU << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE1")=")" << MPI_OPTIONS_MACHINE1 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE1")=")" << MPI_COMMAND_MACHINE1 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE1")=")" << MPI_BINARY_DIR_MACHINE1 << "\"" << endl;
    }

    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE2")=")" << MPI_OPTIONS_MACHINE2 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE2")=")" << MPI_COMMAND_MACHINE2 << "\"" << endl;
    }
    if (LDEBUG) {
      oss << R"(XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE2")=")" << MPI_BINARY_DIR_MACHINE2 << "\"" << endl;
    }

    //   if(LDEBUG) oss << "XHOST.adefault.content=" << XHOST.adefault.content_string << endl;

    if (LDEBUG) {
      oss << "aflowrc::print_aflowrc: END" << endl;
    }

    oss.flush();
    return false;
  }
} // namespace aflowrc

// **************************************************************************
// **************************************************************************
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
