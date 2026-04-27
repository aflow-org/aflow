// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to run VASP in KBIN
// Stefano Curtarolo - 2007-2012 Duke
//   GENERATE, STATIC, KPOINTS, RELAX, RELAX_STATIC, RELAX_STATIC_BANDS, RELAX_STATIC_DIELECTRIC, STATIC_BANDS, STATIC_DIELECTRIC

#include "flow/aflow_kvasp.h"

#include "config.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pthread.h>
#include <unistd.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_argv.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_time.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_init.h"
#include "aflow_xhost.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/SYM/aflow_symmetry.h"
#include "structure/aflow_lattice.h"

using aurostd::RemoveWhiteSpaces;
using std::cerr;
using std::cout;
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

using aurostd::xmatrix;
using aurostd::xoption;
using aurostd::xvector;

#define _incarpad_ 26
#define _KVASP_VASP_SLEEP_ 2
#define _KVASP_WAIT_SLEEP_ 10
#define KBIN_WRONG_ENTRY_STRING string("WRONG_ENTRY")
#define KBIN_WRONG_ENTRY_NUMBER -123

#define _DEBUG_KVASP_ false  // CO20190116

#define _VASP_CONTCAR_SAVE_ true

#define _STROPT_ string("[VASP_FORCE_OPTION]")
#define LDAU_ADIABATIC_RELAX_DEFAULT 6

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

namespace KBIN {
  _kflags VASP_Get_Kflags_from_AflowIN(const string& AflowIn, _aflags& aflags, ostream& oss) { // CO20200624
    ofstream FileMESSAGE("/dev/null");
    return KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn, FileMESSAGE, aflags, oss);
  }
} // namespace KBIN

namespace KBIN {
  _kflags VASP_Get_Kflags_from_AflowIN(const string& _AflowIn, ofstream& FileMESSAGE, _aflags& aflags, ostream& oss) {  // CO20200624
    const bool LDEBUG = (false || _DEBUG_KVASP_ || XHOST.DEBUG);
    _kflags kflags;
    const string AflowIn = aurostd::RemoveComments(_AflowIn); // for safety //CO20180502
    vector<string> vAflowIn;
    aurostd::string2vectorstring(AflowIn, vAflowIn);
    const string BflowIn = AflowIn;
    ostringstream aus;

    if (LDEBUG) {
      cerr << "DEBUG: " << __AFLOW_FUNC__ << " (START)" << endl;
    }

    if (aflags.Directory.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "aflags.Directory not set", _INPUT_MISSING_);
    }  // CO20200624 - prevent segfault
    if (aflags.Directory.at(0) != '/' && aflags.Directory.at(0) != '.' && aflags.Directory.at(0) != ' ') {
      aflags.Directory = "./" + aflags.Directory;
    }

    // ***************************************************************************
    // FIND MPI
    kflags.KBIN_MPI = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MPI]");  // search for MPI string
    // ***************************************************************************
    // FIND HOST
    // duke_beta
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]BETA") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_BETA")) { // check DUKE_BETA
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // duke_beta_openmpi
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]BETA_OPENMPI")
        || // check DUKE_BETA_OPENMPI
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_BETA_OPENMPI")) { // check DUKE_BETA_OPENMPI
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // duke_qrats
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]QRATS") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_QRATS")) { // check DUKE_QRATS
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // duke_qflow
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]QFLOW")
        || // backwards compatible //CO20180409
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_QFLOW")
        || // backwards compatible //CO20180409
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]QUSER")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_QUSER")) { // check DUKE_QFLOW
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // CO20201220 X START
    //  duke_x
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_X")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]X_X")
        || // backwards compatible //CO20180409
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_X_X")) { // check DUKE_X_X //CO20180409
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_X")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_CRAY") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]X_CRAY") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_X_CRAY")) { // check DUKE_X_CRAY //SD20221006
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_CRAY")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_OLDCRAY") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]X_OLDCRAY") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_X_OLDCRAY")) { // check DUKE_X_OLDCRAY //SD20221006
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_OLDCRAY")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_SMB") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]X_SMB") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_X_SMB")) { // check DUKE_X_SMB //SD20221006
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_SMB")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // CO20201220 X STOP
    // CO20220818 JHU_ROCKFISH START
    //  jhu_rockfish
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::JHU_ROCKFISH")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]ROCKFISH")
        || // backwards compatible //CO20220830
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]JHU_ROCKFISH")) { // check JHU_ROCKFISH //CO20220830
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::JHU_ROCKFISH")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // CO20220818 JHU_ROCKFISH STOP
    //  mpcdf_eos
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]EOS") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MPCDF_EOS")) { // check MPCDF_EOS
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // mpcdf_draco
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DRACO") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MPCDF_DRACO")) { // check MPCDF_DRACO
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // mpcdf_cobra
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_COBRA") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]COBRA") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MPCDF_COBRA")) { // check MPCDF_COBRA
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // mpcdf_hydra
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]HYDRA") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MPCDF_HYDRA")) { // check MPCDF_HYDRA
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20190509 - MACHINE001 - START
    //  machine001
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE001") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE001")) { // check MACHINE001
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20190509 - MACHINE001 - END
    // DX20190509 - MACHINE002 - START
    //  machine002
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE002") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE002")) { // check MACHINE002
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20190509 - MACHINE002 - END
    // DX20201005 - MACHINE003 - START
    //  machine003
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE003") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE003") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE003")) { // check MACHINE003
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE003")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20201005 - MACHINE003 - END
    // DX20211011 - MACHINE004 - START
    //  machine004
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE004") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE004") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE004")) { // check MACHINE004
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE004")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20211011 - MACHINE004 - END
    //  duke_materials
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_MATERIALS")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MATERIALS")
        || // check DUKE_MATERIALS
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_MATERIALS")) { // check DUKE_MATERIALS
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // duke_aflowlib
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_AFLOWLIB")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]AFLOWLIB")
        || // check DUKE_AFLOWLIB
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_AFLOWLIB")) { // check DUKE_AFLOWLIB
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // duke_habana
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_HABANA")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]HABANA")
        || // check DUKE_HABANA
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]DUKE_HABANA")) { // check DUKE_HABANA
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // fulton_marylou
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::FULTON_MARYLOU")
        || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MARYLOU")
        || // check FULTON_MARYLOU
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]FULTON_MARYLOU")) { // check FULTON_MARYLOU
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // OL
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::OHAD")
        || // CO20181113
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE2")
        || // check MACHINE2
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE2")) { // check MACHINE2
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) { // CO20181113
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // host1
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::HOST1")
        || // CO20181113
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE1")
        || // check MACHINE1
        aurostd::substring2bool(AflowIn, "[AFLOW_HOST]MACHINE1")) { // check MACHINE1
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) { // CO20181113
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20190107 - CMU EULER - START
    //  cmu_euler
    if (aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::CMU_EULER") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]CMU_EULER") || aurostd::substring2bool(AflowIn, "[AFLOW_HOST]CMU_EULER")) { // check CMU_EULER
      aflags.AFLOW_MACHINE_LOCAL = aflags.AFLOW_MACHINE_GLOBAL;
    }
    if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.KBIN_MPI = true; // overrides the MPI for machines
    }
    // DX20190107 - CMU EULER - END

    // ***************************************************************************
    // OTHER CHECKS FOR MPI
    // machines are done withing the VASP/ALIEN stuff, if necessary
    if (aflags.AFLOW_FORCE_MPI) {
      kflags.KBIN_MPI = true; // forcing
    }
    if (aflags.AFLOW_FORCE_SERIAL) {
      kflags.KBIN_MPI = false; // forcing
    }

    kflags.KBIN_QSUB = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_QSUB]") && !aurostd::substring2bool(AflowIn, "[AFLOW_MODE_QSUB]MODE"); // search for QSUB string
    kflags.KBIN_QSUB_MODE1 = aflags.AFLOW_MODE_QSUB_MODE1 || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_QSUB]MODE1"); // search for QSUB string mode1
    kflags.KBIN_QSUB_MODE2 = aflags.AFLOW_MODE_QSUB_MODE2 || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_QSUB]MODE2"); // search for QSUB string mode2
    kflags.KBIN_QSUB_MODE3 = aflags.AFLOW_MODE_QSUB_MODE3 || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_QSUB]MODE3"); // search for QSUB string mode3
    kflags.AFLOW_MODE_ALIEN = // check ALIEN
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE=ALIEN]")
        || // check ALIEN
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ALIEN]")
        || // check ALIEN
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE]ALIEN"); // check ALIEN
    kflags.AFLOW_MODE_MATLAB = // check MATLAB
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE=MATLAB]")
        || // check MATLAB
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE_MATLAB]")
        || // check MATLAB
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE]MATLAB"); // check MATLAB
    kflags.AFLOW_MODE_VASP = // check VASP
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE=VASP]")
        || // check VASP
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE_VASP]")
        || // check VASP
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE]VASP"); // check VASP
    kflags.AFLOW_MODE_AIMS = // check AIMS
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE=AIMS]")
        || // check AIMS
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE_AIMS]")
        || // check AIMS
        aurostd::substring2bool(AflowIn, "[AFLOW_MODE]AIMS"); // check AIMS
    // CO20180406 - fix generate flags
    if (aflags.KBIN_GEN_GENERAL) {
      if (kflags.AFLOW_MODE_AIMS && !aflags.KBIN_GEN_AIMS_FROM_AFLOWIN) {
        aflags.KBIN_GEN_AIMS_FROM_AFLOWIN = true;
      } // very safe
      if (kflags.AFLOW_MODE_VASP && !aflags.KBIN_GEN_VASP_FROM_AFLOWIN) {
        aflags.KBIN_GEN_VASP_FROM_AFLOWIN = true;
      } // do vasp last, default
    }
    kflags.KBIN_SYMMETRY_CALCULATION = aurostd::substring2bool(AflowIn, "[AFLOW_SYMMETRY]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_SYMMETRY]CALC", true);
    // DX START
    kflags.KBIN_SYMMETRY_NO_SCAN = aurostd::substring2bool(AflowIn, "[AFLOW_SYMMETRY]NO_SCAN", true);
    if (aurostd::substring2bool(AflowIn, "[AFLOW_SYMMETRY]SYM_EPS=", true)) {
      kflags.KBIN_SYMMETRY_EPS = aurostd::substring2utype<double>(AflowIn, "[AFLOW_SYMMETRY]SYM_EPS=", true);
    }
    // DX END
    //  ---------------------------------------------------------
    //  parameters for AAPL - CO20170601
    //  to make backwards compatible, we need to not only look for substring, but need to see if "KAPPA=y"
    //  start with AAPL first, then QHA, then APL, they are mutually exclusive
    aurostd::xoption KBIN_PHONONS_CALCULATION_AAPL;
    KBIN_PHONONS_CALCULATION_AAPL.option = false;
    KBIN_PHONONS_CALCULATION_AAPL.options2entry(AflowIn, string("[AFLOW_AAPL]KAPPA=|[AFLOW_PHONONS]KAPPA="), KBIN_PHONONS_CALCULATION_AAPL.option, KBIN_PHONONS_CALCULATION_AAPL.xscheme); // CO20170601
    KBIN_PHONONS_CALCULATION_AAPL.option |= aurostd::substring2bool(AflowIn, "[AFLOW_AAPL]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_AAPL]CALC", true); // legacy
    kflags.KBIN_PHONONS_CALCULATION_AAPL = KBIN_PHONONS_CALCULATION_AAPL.option;
    // ---------------------------------------------------------
    // parameters for QHA - CO20170601
    // to make backwards compatible, we need to not only look for substring, but need to see if "[AFLOW_QHA]CALC"
    if (!kflags.KBIN_PHONONS_CALCULATION_AAPL) { // mutually exclusive
      kflags.KBIN_PHONONS_CALCULATION_QHA = aurostd::substring2bool(AflowIn, "[AFLOW_QHA]CALC", true) || aurostd::substring2bool(AflowIn, "VASP_QHA]CALC", true);
      /////////////////////////////
      // aurostd::xoption KBIN_PHONONS_CALCULATION_QHA; //PN20180705
      // KBIN_PHONONS_CALCULATION_QHA.option=false; //PN20180705
      // KBIN_PHONONS_CALCULATION_QHA.options2entry(AflowIn, string("[AFLOW_QHA]GRUNEISEN=|[AFLOW_PHONONS]GRUNEISEN="), KBIN_PHONONS_CALCULATION_QHA.option, KBIN_PHONONS_CALCULATION_QHA.xscheme); //CO20170601 //PN20180705
      // KBIN_PHONONS_CALCULATION_QHA.option |= aurostd::substring2bool(AflowIn,"[AFLOW_QHA]CALC",true) || aurostd::substring2bool(AflowIn,"[VASP_QHA]CALC",true); //legacy //PN20180705
      // kflags.KBIN_PHONONS_CALCULATION_QHA  = KBIN_PHONONS_CALCULATION_QHA.option; //PN20180705
    }
    // ---------------------------------------------------------
    // parameters for APL
    // if(LDEBUG) cout << XPID << "KBIN::RUN_Directory: kflags.KBIN_PHONONS_CALCULATION_APL=" << kflags.KBIN_PHONONS_CALCULATION_APL << endl;
    if (!(kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA)) { // mutually exclusive
      kflags.KBIN_PHONONS_CALCULATION_APL = aurostd::substring2bool(AflowIn, "[AFLOW_APL]CALC", true)
                                         || aurostd::substring2bool(AflowIn, "[AFLOW_PHONONS]CALC", true)
                                         || aurostd::substring2bool(AflowIn, "[VASP_PHONONS]CALC", true);
    }
    // if(LDEBUG) cout << XPID << "KBIN::RUN_Directory: kflags.KBIN_PHONONS_CALCULATION_APL=" << kflags.KBIN_PHONONS_CALCULATION_APL << endl;
    // ---------------------------------------------------------
    // parameters for AGL (Debye Model)
    // Cormac created CALCSTRAINORIGIN, so we need to check [AFLOW_AEL]CALC vs. [AFLOW_AEL]CALCSTRAINORIGIN
    // kflags.KBIN_PHONONS_CALCULATION_AGL  = aurostd::substring2bool(AflowIn,"[AFLOW_AGL]CALC",true) || aurostd::substring2bool(AflowIn,"[VASP_AGL]CALC",true) || aurostd::substring2bool(AflowIn,"[AFLOW_GIBBS]CALC",true) || aurostd::substring2bool(AflowIn,"[VASP_GIBBS]CALC",true);
    for (size_t i = 0; i < vAflowIn.size() && !kflags.KBIN_PHONONS_CALCULATION_AGL; i++) {
      if ((aurostd::substring2bool(vAflowIn[i], "[AFLOW_AGL]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_AGL]CALC", true))
          && !(aurostd::substring2bool(vAflowIn[i], "[AFLOW_AGL]CALC_", true)
               || aurostd::substring2bool(vAflowIn[i], "[VASP_AGL]CALC_", true)
               || aurostd::substring2bool(vAflowIn[i], "[AFLOW_AGL]CALCS", true)
               || aurostd::substring2bool(vAflowIn[i], "[VASP_AGL]CALCS", true)
               || false)) {
        kflags.KBIN_PHONONS_CALCULATION_AGL = true;
      }
    }
    // ---------------------------------------------------------
    // parameters for AEL (Elastic constants)
    // Cormac created CALCSTRAINORIGIN, so we need to check [AFLOW_AEL]CALC vs. [AFLOW_AEL]CALCSTRAINORIGIN
    // kflags.KBIN_PHONONS_CALCULATION_AEL  = aurostd::substring2bool(AflowIn,"[AFLOW_AEL]CALC",true) || aurostd::substring2bool(AflowIn,"[VASP_AEL]CALC",true);
    for (size_t i = 0; i < vAflowIn.size() && !kflags.KBIN_PHONONS_CALCULATION_AEL; i++) {
      if ((aurostd::substring2bool(vAflowIn[i], "[AFLOW_AEL]CALC", true) || aurostd::substring2bool(AflowIn, "[VASP_AEL]CALC", true))
          && !(aurostd::substring2bool(vAflowIn[i], "[AFLOW_AEL]CALC_", true)
               || aurostd::substring2bool(vAflowIn[i], "[VASP_AEL]CALC_", true)
               || aurostd::substring2bool(vAflowIn[i], "[AFLOW_AEL]CALCS", true)
               || aurostd::substring2bool(vAflowIn[i], "[VASP_AEL]CALCS", true)
               || false)) {
        kflags.KBIN_PHONONS_CALCULATION_AEL = true;
      }
    }
    // ---------------------------------------------------------
    // Warn user if both APL/AAPL and AEL/AGL flags are set, since they are mutually exclusive
    if ((kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA) && (kflags.KBIN_PHONONS_CALCULATION_AGL || kflags.KBIN_PHONONS_CALCULATION_AEL)) {
      aus << "WWWWW  WARNING: APL/AAPL/QHA and AEL/AGL flags both set" << endl;
      aus << "WWWWW  WARNING: These runs are mutually exclusive" << endl;
      aus << "WWWWW  WARNING: APL/AAPL/QHA runs will take priority" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    } // CT20200520 - added warning
    // ---------------------------------------------------------
    // parameters for NEIGHBORS
    // ---------------------------------------------------------
    // parameters for POCC CALCULATIONS, KESONG YANG
    kflags.KBIN_POCC = false;
    kflags.KBIN_POCC_CALCULATION = aurostd::substring2bool(AflowIn, "[AFLOW_POCC]CALC", true)
                                && (aurostd::substring2bool(AflowIn, "[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE", true) && aurostd::substring2bool(AflowIn, "[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE", true)); // CO20180419
    if (kflags.KBIN_POCC_CALCULATION) {
      aus << "00000  MESSAGE POCC_CALCULATION " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    }
    if (kflags.KBIN_POCC_CALCULATION) {
      kflags.KBIN_POCC = true;
    } // CO20180419
    if (kflags.KBIN_POCC) { // CO20191110
      kflags.KBIN_POCC_TEMPERATURE_STRING = aurostd::substring2string(AflowIn, "[AFLOW_POCC]TEMPERATURE="); // CO20191110
      if (kflags.KBIN_POCC_TEMPERATURE_STRING.empty()) { // CO20191110
        kflags.KBIN_POCC_TEMPERATURE_STRING = DEFAULT_POCC_TEMPERATURE_STRING; // CO20191110
      }
      aus << "00000  MESSAGE POCC_TEMPERATURE_STRING=" << kflags.KBIN_POCC_TEMPERATURE_STRING << Message(__AFLOW_FILE__, aflags) << endl; // CO20191110
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss); // CO20191110
      kflags.KBIN_POCC_ARUNS2SKIP_STRING = aurostd::substring2string(AflowIn, "[AFLOW_POCC]ARUNS2SKIP="); // CO20200624
      if (!kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) { // CO20200624
        aus << "00000  MESSAGE POCC_ARUNS2SKIP_STRING=" << kflags.KBIN_POCC_ARUNS2SKIP_STRING << Message(__AFLOW_FILE__, aflags) << endl; // CO20200624
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss); // CO20200624
      }
      // ME20210927 - EXCLUDE_UNSTABLE
      string exclude = aurostd::toupper(aurostd::substring2string(AflowIn, "[AFLOW_POCC]EXCLUDE_UNSTABLE="));
      if (exclude.empty()) {
        kflags.KBIN_POCC_EXCLUDE_UNSTABLE = DEFAULT_POCC_EXCLUDE_UNSTABLE;
      } else {
        if (exclude[0] == 'T') {
          kflags.KBIN_POCC_EXCLUDE_UNSTABLE = true;
        } else if (exclude[0] == 'F') {
          kflags.KBIN_POCC_EXCLUDE_UNSTABLE = false;
        } else {
          kflags.KBIN_POCC_EXCLUDE_UNSTABLE = DEFAULT_POCC_EXCLUDE_UNSTABLE;
        }
        aus << "00000  MESSAGE POCC_EXCLUDE_UNSTABLE=" << (kflags.KBIN_POCC_EXCLUDE_UNSTABLE ? "true" : "false") << Message(__AFLOW_FILE__, aflags) << endl;
      }
    }
    // ---------------------------------------------------------
    // parameters for FROZSL
    kflags.KBIN_FROZSL = false;
    kflags.KBIN_PHONONS_CALCULATION_FROZSL = aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]CALC", true);
    if (kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
      aus << "00000  MESSAGE FROZSL_CALCULATION " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    }
    kflags.KBIN_FROZSL_DOWNLOAD = (aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]DOWN", true) || aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]DOWNLOAD", true));
    if (kflags.KBIN_FROZSL_DOWNLOAD) {
      aus << "00000  MESSAGE FROZSL_DOWNLOAD " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    }
    if (kflags.KBIN_FROZSL_FILE) {
      aus << "00000  MESSAGE FROZSL_FILE " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    }
    kflags.KBIN_FROZSL_FILE = aurostd::substring2bool(AflowIn, "[AFLOW_FROZSL]FILE", true); // load of file somewhere else
    if (kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_FROZSL_DOWNLOAD || kflags.KBIN_FROZSL_FILE) {
      kflags.KBIN_FROZSL = true;
    }
    // ---------------------------------------------------------
    // the rest of symmetry stuff is sought inside ivasp or
    if (kflags.AFLOW_MODE_ALIEN) {
      kflags.AFLOW_MODE_MATLAB = false; // fix PRIORITY
      kflags.AFLOW_MODE_VASP = false; // fix PRIORITY
      kflags.KBIN_MPI = false; // fix PRIORITY
    }
    if (kflags.AFLOW_MODE_MATLAB) {
      kflags.AFLOW_MODE_VASP = false; // fix PRIORITY
      kflags.KBIN_MPI = false;
    }
    if (LDEBUG) {
      cout << "DEBUG kflags.AFLOW_MODE_ALIEN=" << kflags.AFLOW_MODE_ALIEN << endl;
    }
    if (LDEBUG) {
      cout << "DEBUG kflags.AFLOW_MODE_MATLAB=" << kflags.AFLOW_MODE_MATLAB << endl;
    }
    if (LDEBUG) {
      cout << "DEBUG kflags.AFLOW_MODE_VASP=" << kflags.AFLOW_MODE_VASP << endl;
    }
    // ***************************************************************************
    // ZIP COMPRESS
    // ***************************************************************************
    kflags.KZIP_COMPRESS = true;
    aurostd::StringstreamClean(aus);
    if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP=none]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP=NONE]") || !aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP")) {
      kflags.KZIP_COMPRESS = false;
      for (int i = 0; i < 1; i++) {
        aus << "WWWWW  WARNING no compression of output files..." << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintWarningStream(FileMESSAGE, aus, XHOST.QUIET);
      }
    } else {
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP")) { // "[AFLOW_MODE_ZIP=" not found
        kflags.KZIP_BIN = DEFAULT_KZIP_BIN; // take default
        aus << "00000  MESSAGE Taking DEFAULT KZIP_BIN=\"" << kflags.KZIP_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP]")) { // "[AFLOW_MODE_ZIP]" not found
        kflags.KZIP_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_ZIP]");
        aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_ZIP=")) { // "[AFLOW_MODE_ZIP=" found
        kflags.KZIP_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_ZIP="), ']');
        aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
    }
    // ************************************************************************************************************************************
    // Get the KZIP_BIN name - moved inside EACH mode
    // ************************************************************************************************************************************
    // LOAD PREFIX POSTFIX
    KBIN::StartStopCheck(AflowIn, "[AFLOW_MODE_PRESCRIPT]", kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT, kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP);
    KBIN::StartStopCheck(AflowIn, "[AFLOW_MODE_POSTSCRIPT]", kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT, kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP);
    if (kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT) { // [AFLOW_MODE_PRESCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_PRESCRIPT_COMMAND << " file from " << _AFLOWIN_ << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.AFLOW_MODE_PRESCRIPT.str(aurostd::substring2string(AflowIn, "[AFLOW_MODE_PRESCRIPT]", 0));
    }
    if (kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) { // [AFLOW_MODE_PRESCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_PRESCRIPT_COMMAND << " file from START/STOP " << _AFLOWIN_ << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.AFLOW_MODE_PRESCRIPT.str(aurostd::substring2string(AflowIn, "[AFLOW_MODE_PRESCRIPT]START", "[AFLOW_MODE_PRESCRIPT]STOP", 1));
    }
    if (kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT) { // [AFLOW_MODE_POSTSCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_POSTSCRIPT_COMMAND << " file from " << _AFLOWIN_ << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.AFLOW_MODE_POSTSCRIPT.str(aurostd::substring2string(AflowIn, "[AFLOW_MODE_POSTSCRIPT]", 0));
    }
    if (kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) { // [AFLOW_MODE_POSTSCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_POSTSCRIPT_COMMAND << " file from START/STOP " << _AFLOWIN_ << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      kflags.AFLOW_MODE_POSTSCRIPT.str(aurostd::substring2string(AflowIn, "[AFLOW_MODE_POSTSCRIPT]START", "[AFLOW_MODE_POSTSCRIPT]STOP", 1));
    }
    // ************************************************************************************************************************************
    // ALIEN MODE
    if (kflags.AFLOW_MODE_ALIEN) {
      aus << XPID << "00000  MESSAGE [AFLOW_MODE=ALIEN] found in " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      // ***************************************************************************
      // Getting KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_KBIN_ALIEN_BIN; // take default
      aurostd::StringstreamClean(aus);
      if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
        kflags.KBIN_BIN = DEFAULT_KBIN_ALIEN_BIN; // take default
        aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
        kflags.KBIN_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY]");
        aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
        kflags.KBIN_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY="), ']');
        aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
      // ME20190107 - Grab the serial binary to propagate into child aflow.in files
      kflags.KBIN_SERIAL_BIN = kflags.KBIN_BIN;
      // ***************************************************************************
      // ALIEN MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_EMAIL]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE]EMAIL");
      // ***************************************************************************
    }
    // ************************************************************************************************************************************
    // MATLAB MODE
    if (kflags.AFLOW_MODE_MATLAB) {
      aus << XPID << "00000  MESSAGE [AFLOW_MODE=MATLAB] found in " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      // ***************************************************************************
      // MATLAB MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_EMAIL]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE]EMAIL");
      // ***************************************************************************
    }
    // ************************************************************************************************************************************
    // AIMS MODE
    if (kflags.AFLOW_MODE_AIMS) {
      aus << XPID << "00000  MESSAGE [AFLOW_MODE=AIMS] found in " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      aurostd::StringstreamClean(aus);
      // ***************************************************************************
      // Getting KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_AIMS_BIN; // take default  dont touch MPI as it has already been dealt by  KBIN::MPI_Extract
      if (kflags.KBIN_MPI == false) { // only if no MPI is specified
        if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
          kflags.KBIN_BIN = DEFAULT_AIMS_BIN; // take default
          aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
          kflags.KBIN_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY]");
          aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
          kflags.KBIN_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY="), ']');
          aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
        }
        // ME20190107 - Grab the serial binary to propagate into child aflow.in files
        kflags.KBIN_SERIAL_BIN = kflags.KBIN_BIN;
      } else {
        kflags.KBIN_BIN = kflags.KBIN_MPI_BIN;
        aus << "00000  MESSAGE Taking KBIN_BIN=KBIN_MPI_BIN=\"" << kflags.KBIN_MPI_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }

      // ***************************************************************************
      // AIMS MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_EMAIL]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE]EMAIL");
      // ***************************************************************************
    }
    // ************************************************************************************************************************************
    // MPI SWTICHES
    if (kflags.KBIN_MPI) {
      KBIN::MPI_Extract(AflowIn, FileMESSAGE, aflags, kflags);
    }
    // ************************************************************************************************************************************
    // ************************************************************************************************************************************
    // VASP MODE
    if (kflags.AFLOW_MODE_VASP) {
      aus << XPID << "00000  MESSAGE [AFLOW_MODE=VASP] found in " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      // ***************************************************************************
      // Getting KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_VASP_BIN; // take default  dont touch MPI as it has already been dealt by  KBIN::MPI_Extract
      aurostd::StringstreamClean(aus);
      if (kflags.KBIN_MPI == false) { // only if no MPI is specified
        if (!aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
          kflags.KBIN_BIN = DEFAULT_VASP_BIN; // take default
          aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
          kflags.KBIN_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY]");
          aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
        }
        if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
          kflags.KBIN_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY="), ']');
          aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
        }
        // ME20190107 - Grab the serial binary to propagate into child aflow.in files
        kflags.KBIN_SERIAL_BIN = kflags.KBIN_BIN;
      } else {
        kflags.KBIN_BIN = kflags.KBIN_MPI_BIN;
        aus << "00000  MESSAGE Taking KBIN_BIN=KBIN_MPI_BIN=\"" << kflags.KBIN_MPI_BIN << "\" " << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
      }
      // ***************************************************************************
      // VASP MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL = aurostd::substring2bool(AflowIn, "[AFLOW_MODE_EMAIL]") || aurostd::substring2bool(AflowIn, "[AFLOW_MODE]EMAIL");
    }
    // ***************************************************************************
    // ************************************************************************************************************************************
    // MATLAB MODE
    if (kflags.KBIN_PHONONS_CALCULATION_FROZSL && !kflags.AFLOW_MODE_VASP) {
      aus << XPID << "00000  MESSAGE [AFLOW_FROZSL]CALC found in " << _AFLOWIN_ << " " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    }
    // ************************************************************************************************************************************
    // NO MODE MODE
    if (!kflags.AFLOW_MODE_VASP && !kflags.AFLOW_MODE_AIMS && !kflags.AFLOW_MODE_MATLAB && !kflags.AFLOW_MODE_ALIEN && !kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
      aus << "EEEEE  [AFLOW_MODE=????] invalid found in     " << Message(__AFLOW_FILE__, aflags) << endl;
      aus << "EEEEE  [AFLOW_MODE=ALIEN]        is supported " << Message(__AFLOW_FILE__, aflags) << endl;
      aus << "EEEEE  [AFLOW_MODE=MATLAB]       is supported " << Message(__AFLOW_FILE__, aflags) << endl;
      aus << "EEEEE  [AFLOW_MODE=VASP]         is supported " << Message(__AFLOW_FILE__, aflags) << endl;
      aus << "EEEEE  [AFLOW_FROZSL]CALC        is supported " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    // ***************************************************************************
    // FINALIZE LOCK
    aus << "XXXXX  KBIN DIRECTORY END (aflow" << string(AFLOW_VERSION) << ")  " << Message(__AFLOW_FILE__, aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    // ***************************************************************************
    // PREPARE MESSAGE FOR LOG TO BE INTERCEPTED IN COMPRESSION
    aus << "XXXXX  KBIN_DIRECTORY_END " << aflags.Directory << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET, oss);
    // ***************************************************************************

    if (LDEBUG) {
      cerr << "DEBUG: " << __AFLOW_FUNC__ << " (END)" << endl;
    }

    return kflags;
  }
} // namespace KBIN

namespace KBIN {
  _vflags VASP_Get_Vflags_from_AflowIN(const string& AflowIn, _aflags& aflags, _kflags& kflags, ostream& oss) {
    ofstream FileMESSAGE("/dev/null");
    return KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, FileMESSAGE, aflags, kflags, oss);
  }
} // namespace KBIN

namespace KBIN {
  /// @brief Search the aflow.in for flags related to VASP runs
  /// @param _AflowIn aflow.in to search
  /// @param FileMESSAGE file stream for logging
  /// @param aflags general aflow flags
  /// @param kflags general flags for kbinary
  /// @param oss out stream for logging
  /// @return VASP flags
  _vflags VASP_Get_Vflags_from_AflowIN(const string& _AflowIn, ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, ostream& oss) {
    const bool LDEBUG = (false || _DEBUG_KVASP_ || XHOST.DEBUG);
    string message;
    _vflags vflags;
    const string AflowIn = aurostd::RemoveComments(_AflowIn); // for safety //CO20180502
    vector<string> vAflowIn;
    aurostd::string2vectorstring(AflowIn, vAflowIn);
    string BflowIn = AflowIn;

    if (LDEBUG) {
      cerr << "DEBUG: " << __AFLOW_FUNC__ << " (START)" << endl;
    }
    // HOW TO RUN
    vflags.KBIN_VASP_RUN_NRELAX = 0;
    vflags.KBIN_VASP_RUN.clear();

    std::vector<std::pair<std::string, bool>> RUN_ordered_vec = {
        {               "GENERATE", false},
        {     "RELAX_STATIC_BANDS", false},
        {"RELAX_STATIC_DIELECTRIC", false},
        {           "RELAX_STATIC", false},
        {           "STATIC_BANDS", false},
        {      "STATIC_DIELECTRIC", false},
        {                 "STATIC", false},
        {                "KPOINTS", false},
        {                  "RELAX", false}
    };
    std::vector<std::pair<std::string, bool>> REPEAT_ordered_vec = {
        {           "REPEAT_STATIC", false},
        {     "REPEAT_STATIC_BANDS", false},
        {"REPEAT_STATIC_DIELECTRIC", false},
        {            "REPEAT_BANDS", false},
        {       "REPEAT_DIELECTRIC", false},
        {           "REPEAT_DELSOL", false}
    };

    std::smatch match;
    for (auto& pair : RUN_ordered_vec) {
      const std::regex re(R"((^|\n)\s*?\[VASP_RUN(_|\]))" + pair.first + "(=|[^_]|\\s)(\\d+|)");
      if (std::regex_search(AflowIn, match, re)) {
        pair.second = true;
        if (pair.first.find("RELAX") != std::string::npos) {
          vflags.KBIN_VASP_RUN_NRELAX = aurostd::string2utype<int>(match[4].str());
        }
      }
    }
    if (aflags.KBIN_GEN_VASP_FROM_AFLOWIN) {
      vflags.KBIN_VASP_RUN.push("GENERATE");
    }

    if (!vflags.KBIN_VASP_RUN.xscheme.empty()) {
      vflags.KBIN_VASP_RUN.isentry = true;
    }

    vflags.KBIN_VASP_REPEAT.clear();
    for (auto& pair : REPEAT_ordered_vec) {
      const std::regex re(R"((^|\n)\s*?\[VASP_RUN(_|\]))" + pair.first);
      std::smatch match;
      if (std::regex_search(AflowIn, match, re) || aurostd::FileExist(aflags.Directory + "/" + pair.first)) {
        pair.second = true;
      }
    }

    // priorities about RUN
    for (const auto& pair1 : RUN_ordered_vec) {
      if (pair1.second) {
        vflags.KBIN_VASP_RUN.push(pair1.first);
        vflags.KBIN_VASP_RUN.flag(pair1.first, true);
        for (const auto& pair2 : RUN_ordered_vec) {
          if (pair1.first != pair2.first) {
            vflags.KBIN_VASP_RUN.flag(pair2.first, false);
          }
        }
        break;
      }
    }

    // priorities about REPEAT
    for (const auto& pair1 : REPEAT_ordered_vec) {
      if (pair1.second) {
        vflags.KBIN_VASP_REPEAT.push(pair1.first);
        vflags.KBIN_VASP_REPEAT.flag(pair1.first, true);
        for (const auto& pair2 : REPEAT_ordered_vec) {
          if (pair1.first != pair2.first) {
            vflags.KBIN_VASP_REPEAT.flag(pair2.first, false);
          }
        }
        for (const auto& pair3 : RUN_ordered_vec) {
          vflags.KBIN_VASP_RUN.flag(pair3.first, false);
        }
        break;
      }
    }

    if (kflags.KBIN_PHONONS_CALCULATION_APL
        || kflags.KBIN_PHONONS_CALCULATION_QHA
        || kflags.KBIN_PHONONS_CALCULATION_AAPL
        || kflags.KBIN_PHONONS_CALCULATION_FROZSL
        || kflags.KBIN_PHONONS_CALCULATION_AGL
        || kflags.KBIN_PHONONS_CALCULATION_AEL) {
      for (const auto& pair : RUN_ordered_vec) {
        vflags.KBIN_VASP_RUN.flag(pair.first, false);
      }
      for (const auto& pair : REPEAT_ordered_vec) {
        vflags.KBIN_VASP_REPEAT.flag(pair.first, false);
      }
      kflags.KBIN_SYMMETRY_CALCULATION = false;
    }

    // RELAX_MODE AND PRIORITIES  // ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY (default ENERGY) "
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.options2entry(AflowIn, _STROPT_ + "RELAX_MODE=", false, DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme = KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY", "ENERGY");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCES", "FORCES");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCE", "FORCES");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY_FORCES", "ENERGY_FORCES");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY_FORCE", "ENERGY_FORCES");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCES_ENERGY", "FORCES_ENERGY");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCE_ENERGY", "FORCES_ENERGY");
    if (vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.isentry && vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // FORCE OPTIONS
    vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.options2entry(AflowIn, _STROPT_ + "NOTUNE");

    // FORCE OPTIONS SYSTEM_AUTO
    vflags.KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO.options2entry(AflowIn, _STROPT_ + "SYSTEM_AUTO");
    vflags.AFLOW_SYSTEM.options2entry(AflowIn, "[AFLOW]SYSTEM=", false, ""); // ME20181121

    // FORCE OPTIONS STATIC RELAX_ALL RELAX_IONS RELAX CELL_VOLUME
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "STATIC", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("STATIC");
    }
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_ALL", true) || aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("ALL");
    }
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_IONS", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS");
    }
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_CELL_SHAPE", true) || aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_SHAPE", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_SHAPE");
    }
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_CELL_VOLUME", true) || aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_VOLUME", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_VOLUME");
    }
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_IONS_CELL_VOLUME", true) || aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_IONS_VOLUME", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_VOLUME");
    }
    // AS20201123 BEGIN
    if (aurostd::substring2bool(AflowIn, _STROPT_ + "RELAX_IONS_CELL_SHAPE", true)) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_SHAPE");
    }
    // AS20201123 END
    if (!vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.xscheme.empty()) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry = true;
    }

    if (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS", false);
    }
    if (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL")
        && (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC")
            || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")
            || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")
            || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")
            || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL", false);
    }
    if (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC")) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("STATIC");
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry = true;
    }

    // PRECISION AND PRIORITIES // (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED
    vflags.KBIN_VASP_FORCE_OPTION_PREC.options2entry(AflowIn, _STROPT_ + "PREC=", false, DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme = KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('L', "LOW");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('M', "MEDIUM");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('N', "NORMAL");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('H', "HIGH");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('A', "ACCURATE");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('P', "PHONONS"); // JJPR Modification
    if (vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry && vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_PREC.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_PREC.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ALGO AND PRIORITIES // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.options2entry(AflowIn, _STROPT_ + "ALGO=", false, DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme = KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('N', "NORMAL");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('V', "VERYFAST");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('F', "FAST");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('A', "ALL");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('D', "DAMPED");
    if (vflags.KBIN_VASP_FORCE_OPTION_ALGO.isentry && vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_ALGO.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_ALGO.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved = vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved || aurostd::substring2bool(AflowIn, _STROPT_ + "ALGO_PRESERVED", true); // FIX ALGO_PRESERVED

    // ABMIX AND PRIORITIES // empty | [AUTO | US | PAW | #AMIX,#BMIX[,#AMIX_MAG,#BMIX_MAG]]
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.options2entry(AflowIn, _STROPT_ + "ABMIX=", false, DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme = KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('A', "AUTO");
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('U', "US");
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('L', "US"); // LDA
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('G', "US"); // GGA
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('P', "PAW"); // something with PAW..
    if (vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry && vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_ABMIX.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_ABMIX.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // METAGGA AND PRIORITIES // TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " METAGGA" << endl;
    }
    vflags.KBIN_VASP_FORCE_OPTION_METAGGA.options2entry(AflowIn, _STROPT_ + "METAGGA=", false, DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME);
    if (vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry && vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_METAGGA.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_METAGGA.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " METAGGA vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " METAGGA vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << endl;
    }

    // IVDW AND PRIORITIES // [number_for_VASP_see_manual_for_IVDW | 0]
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " IVDW" << endl;
    }
    vflags.KBIN_VASP_FORCE_OPTION_IVDW.options2entry(AflowIn, _STROPT_ + "IVDW=", false, DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME);
    if (vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry && vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_IVDW.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_IVDW.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " IVDW vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " IVDW vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << endl;
    }

    // NEGLECT_NOMIX
    vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.options2entry(AflowIn, string(_STROPT_ + "NEGLECT_IMMISCIBLE" + "|" + _STROPT_ + "NEGLECT_NOMIX" + "|" + _STROPT_ + "SKIP_NOMIX"));

    // AUTO_PSEUDOPOTENTIALS and AUTO_PSEUDOPOTENTIALS_TYPE
    vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.options2entry(AflowIn, _STROPT_ + "AUTO_PSEUDOPOTENTIALS=", false, DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE);

    // POTIM
    vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.options2entry(AflowIn, _STROPT_ + "POTIM=", true, vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_VASP_PREC_POTIM"

    // PSTRESS
    vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.options2entry(AflowIn, _STROPT_ + "PSTRESS=", true, vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.0"

    // EDIFFG
    vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.options2entry(AflowIn, _STROPT_ + "EDIFFG=", true, vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_VASP_PREC_EDIFFG"

    // NELM //CO20200624
    vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.options2entry(AflowIn, _STROPT_ + "NELM=", false, vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "60" - default //CO20200624

    // NELM_STATIC //CO20200624
    vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.options2entry(AflowIn, _STROPT_ + "NELM_STATIC=", false,
                                                                  vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "120" - default  //CO20200624

    // ISMEAR
    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.options2entry(AflowIn, _STROPT_ + "ISMEAR=", false, vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1" - default //CO20181128

    // SIGMA
    vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.options2entry(AflowIn, _STROPT_ + "SIGMA=", false, vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.1" - default //CO20181128

    // ISMEAR_STATIC  //CO20210315
    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.options2entry(AflowIn, _STROPT_ + "ISMEAR_STATIC=", false,
                                                                    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1" - default  //CO20181128

    // SIGMA_STATIC //CO20210315
    vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.options2entry(AflowIn, _STROPT_ + "SIGMA_STATIC=", false,
                                                                   vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.1" - default  //CO20181128

    // ISMEAR_BANDS  //CO20210315
    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.options2entry(AflowIn, _STROPT_ + "ISMEAR_BANDS=" + "|" + _STROPT_ + "ISMEAR_STATIC_BANDS=", false,
                                                                   vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1" - default  //CO20181128  //CO20210624 - backwards compatibility with ISMEAR_STATIC_BANDS

    // SIGMA_BANDS //CO20210315
    vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.options2entry(AflowIn, _STROPT_ + "SIGMA_BANDS=" + "|" + _STROPT_ + "SIGMA_STATIC_BANDS=", false,
                                                                  vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.1" - default  //CO20181128  //CO20210624 - backwards compatibility with SIGMA_STATIC_BANDS

    // NBANDS and/or NBANDS=
    vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry = aurostd::substring2bool(AflowIn, _STROPT_ + "NBANDS", true);
    vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.options2entry(AflowIn, _STROPT_ + "NBANDS=", true, vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"
    if (vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry) {
      vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry = false;
    }

    vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.options2entry(AflowIn, _STROPT_ + "ENMAX_MULTIPLY=", true, vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.0"

    // RWIGS_STATIC
    vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC = aurostd::substring2bool(AflowIn, _STROPT_ + "RWIGS_STATIC", true);

    // SPIN AND PRIORITIES // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1; // DEFAULT
    vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2; // DEFAULT
    vflags.KBIN_VASP_FORCE_OPTION_SPIN.options2entry(AflowIn, _STROPT_ + "SPIN=", DEFAULT_VASP_FORCE_OPTION_SPIN);
    if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) { // ME+RF20200225; fixes bug that SPIN was switched OFF in static calc. when SPIN=ON in aflow.in and REMOVE_RELAX was set by default
      if (!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {
        if (aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string, "REMOVE_RELAX_1") || aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string, "REMOVE_RELAX_2")) {
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "SPIN is OFF. REMOVE_RELAX_1/2 will be switched off.", aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
          vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = false;
          vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = false;
        }
      } else {
        vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = false;
        vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = false;
        if (aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string, "REMOVE_RELAX_1")) {
          vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = true;
        }
        if (aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string, "REMOVE_RELAX_2")) {
          vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = true;
        }
      }
    }
    if (!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {
      vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = false; // nothing to remove
    }
    if (!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {
      vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = false; // nothing to remove
    }

    // BADER AND PRIORITIES // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_BADER.options2entry(AflowIn, _STROPT_ + "BADER=", DEFAULT_VASP_FORCE_OPTION_BADER);
    if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {
      vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry = true; // DEFAULT
    }
    if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {
      vflags.KBIN_VASP_FORCE_OPTION_BADER.option = true; // DEFAULT
    }

    // ELF AND PRIORITIES // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_ELF.options2entry(AflowIn, _STROPT_ + "ELF=", DEFAULT_VASP_FORCE_OPTION_ELF);

    // AUTO_MAGMOM AND PRIORITIES  // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.options2entry(AflowIn, _STROPT_ + "AUTO_MAGMOM=", DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);
    // LSCOUPLING AND PRIORITIES  // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.options2entry(AflowIn, _STROPT_ + "LSCOUPLING=", DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);
    if (vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option) {
      if (!aurostd::substring2bool(kflags.KBIN_BIN, "LS") && !aurostd::substring2bool(kflags.KBIN_BIN, "ls")) {
        kflags.KBIN_BIN += "LS";
      }
      if (!aurostd::substring2bool(kflags.KBIN_MPI_BIN, "LS") && !aurostd::substring2bool(kflags.KBIN_MPI_BIN, "ls")) {
        kflags.KBIN_MPI_BIN += "LS";
      }
      kflags.KBIN_BIN = aurostd::RemoveCharacter(kflags.KBIN_BIN, ' '); // if there is junk
      kflags.KBIN_MPI_BIN = aurostd::RemoveCharacter(kflags.KBIN_MPI_BIN, ' '); // if there is junk
    }
    // SYM AND PRIORITIES
    vflags.KBIN_VASP_FORCE_OPTION_SYM.options2entry(AflowIn, _STROPT_ + "SYM=", DEFAULT_VASP_FORCE_OPTION_SYM);
    // WAVECAR AND PRIORITIES
    vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.options2entry(AflowIn, _STROPT_ + "WAVECAR=", DEFAULT_VASP_FORCE_OPTION_WAVECAR);
    // CHGCAR AND PRIORITIES
    vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.options2entry(AflowIn, _STROPT_ + "CHGCAR=", DEFAULT_VASP_FORCE_OPTION_CHGCAR);
    // ME20191028 - specify CHGCAR file to use
    vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.options2entry(AflowIn, _STROPT_ + "CHGCAR_FILE=", 0, "");

    // LDAU2 AND PRIORITIES
    vflags.KBIN_VASP_LDAU_SPECIES = "";
    vflags.KBIN_VASP_LDAU_PARAMETERS = "";
    vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag = true;

    BflowIn = AflowIn;
    aurostd::StringSubstInPlace(BflowIn, "LDAU1=", "LDAU=");
    aurostd::StringSubstInPlace(BflowIn, "LDAU2=", "LDAU=");
    vflags.KBIN_VASP_FORCE_OPTION_LDAU0.options2entry(BflowIn, string(_STROPT_ + "LDAU=OFF" + "|" + _STROPT_ + "LDAU=0" + "|" + _STROPT_ + "LDAU=N" + "|" + _STROPT_ + "LDAU=false"));
    vflags.KBIN_VASP_FORCE_OPTION_LDAU1
        .options2entry(AflowIn, string(_STROPT_ + "LDAU1=ON" + "|" + _STROPT_ + "LDAU1=1" + "|" + "LDAU1=Y" + "|" + _STROPT_ + "LDAU1=true" + "|" + _STROPT_ + "LDAU1=ADIABATIC" + "|" + _STROPT_ + "LDAU1=CUTOFF"));
    vflags.KBIN_VASP_FORCE_OPTION_LDAU2
        .options2entry(AflowIn, string(_STROPT_ + "LDAU2=ON" + "|" + _STROPT_ + "LDAU2=1" + "|" + "LDAU2=Y" + "|" + _STROPT_ + "LDAU2=true" + "|" + _STROPT_ + "LDAU2=ADIABATIC" + "|" + _STROPT_ + "LDAU2=CUTOFF"));
    if (vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
      vflags.KBIN_VASP_FORCE_OPTION_LDAU0.isentry = false;
    }
    if (vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
      if (aurostd::substring2bool(AflowIn, _STROPT_ + "LDAU_SPECIES=", true)) {
        vflags.KBIN_VASP_LDAU_SPECIES = aurostd::substring2string(AflowIn, _STROPT_ + "LDAU_SPECIES=", 1, false);
      }
      if (aurostd::substring2bool(AflowIn, _STROPT_ + "LDAU1_SPECIES=", true)) {
        vflags.KBIN_VASP_LDAU_SPECIES = aurostd::substring2string(AflowIn, _STROPT_ + "LDAU1_SPECIES=", 1, false);
      }
      if (aurostd::substring2bool(AflowIn, _STROPT_ + "LDAU2_SPECIES=", true)) {
        vflags.KBIN_VASP_LDAU_SPECIES = aurostd::substring2string(AflowIn, _STROPT_ + "LDAU2_SPECIES=", 1, false);
      }
      if (aurostd::substring2bool(AflowIn, _STROPT_ + "LDAU_PARAMETERS=", true)) {
        vflags.KBIN_VASP_LDAU_PARAMETERS = RemoveWhiteSpaces(aurostd::substring2string(AflowIn, _STROPT_ + "LDAU_PARAMETERS=", 1, false));
      }
      if (!vflags.KBIN_VASP_LDAU_SPECIES.empty()) {
        vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag = true;
      }
      if (!vflags.KBIN_VASP_LDAU_PARAMETERS.empty()) {
        vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag = false;
      }
    }
    // ADIABATIC
    vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.options2entry(AflowIn, string(_STROPT_ + "LDAU1=ADIABATIC" + "|" + _STROPT_ + "LDAU2=ADIABATIC" + "|" + _STROPT_ + "LDAU=ADIABATIC"));
    if (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry) {
      if (vflags.KBIN_VASP_RUN_NRELAX < LDAU_ADIABATIC_RELAX_DEFAULT) {
        vflags.KBIN_VASP_RUN_NRELAX = LDAU_ADIABATIC_RELAX_DEFAULT;
      }
      vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int = vflags.KBIN_VASP_RUN_NRELAX;
    }
    // CUTOFF
    vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.options2entry(AflowIn, string(_STROPT_ + "LDAU1=CUTOFF" + "|" + _STROPT_ + "LDAU2=CUTOFF" + "|" + _STROPT_ + "LDAU=CUTOFF"));
    if (vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry) {
      vflags.KBIN_VASP_RUN_NRELAX++;
    }
    // KPOINTS
    BflowIn = AflowIn;
    aurostd::StringSubstInPlace(BflowIn, "=", "_");
    aurostd::StringSubstInPlace(BflowIn, "KPOINTS_", "KPOINTS="); // bypass for getting all "_"
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.options2entry(BflowIn, string(_STROPT_ + "KPOINTS="), aurostd_xoptionMULTI, ""); // stack them all

    // TYPE AND PRIORITIES // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.options2entry(AflowIn, _STROPT_ + "TYPE=", false, DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme = KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('M', "METAL");
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('I', "INSULATOR");
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('S', "SEMICONDUCTOR");
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('D', "DEFAULT");
    if (vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry && vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme == KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_TYPE.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_TYPE.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // PARAMETERS FOR INCAR
    // NSW=
    vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL = aurostd::substring2bool(AflowIn, _STROPT_ + "NSW=", true);
    if (vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL) {
      vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE = aurostd::substring2utype<int>(AflowIn, _STROPT_ + "NSW=");
    } else {
      vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE = 0;
    }

    // IGNORE_AFIX stuff
    BflowIn = AflowIn;
    aurostd::StringSubstInPlace(BflowIn, "=", "_");
    aurostd::StringSubstInPlace(BflowIn, "IGNORE_AFIX_", "IGNORE_AFIX="); // bypass for getting all "_"
    vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.options2entry(BflowIn, string(_STROPT_ + "IGNORE_AFIX="), aurostd_xoptionMULTI, ""); // stack them all

    // INPUT FILES

    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_FILE]")) {
      vflags.KBIN_VASP_INCAR_FILE.push("KEYWORD");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_FILE]SYSTEM_AUTO", true)) {
      vflags.KBIN_VASP_INCAR_FILE.push("SYSTEM_AUTO");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_FILE]FILE=", true)) { // ME20181113
      const string file = aurostd::substring2string(AflowIn, "[VASP_INCAR_FILE]FILE=", 1, true);
      vflags.KBIN_VASP_INCAR_FILE.push_attached("FILE", file);
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_FILE]COMMAND=", true)) { // ME20181113
      const string command = aurostd::substring2string(AflowIn, "[VASP_INCAR_FILE]COMMAND=", 1, true);
      vflags.KBIN_VASP_INCAR_FILE.push_attached("COMMAND", command);
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_MODE_EXPLICIT]")) {
      vflags.KBIN_VASP_INCAR_MODE.push("EXPLICIT");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn, "[VASP_INCAR_MODE_EXPLICIT]STOP")) {
      vflags.KBIN_VASP_INCAR_MODE.push("EXPLICIT_START_STOP");
    }
    if (vflags.KBIN_VASP_INCAR_FILE.flag("KEYWORD")) { // ME20181113
      vflags.KBIN_VASP_INCAR_EXPLICIT.str(aurostd::substring2string(AflowIn, "[VASP_INCAR_FILE]", 0));
    }
    if (vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT_START_STOP")) { // ME20181113
      vflags.KBIN_VASP_INCAR_EXPLICIT_START_STOP.str(aurostd::substring2string(AflowIn, "[VASP_INCAR_MODE_EXPLICIT]START", "[VASP_INCAR_MODE_EXPLICIT]STOP", 1));
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_MODE_IMPLICIT]")) {
      vflags.KBIN_VASP_INCAR_MODE.push("IMPLICIT");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_INCAR_MODE_EXTERNAL]")) {
      vflags.KBIN_VASP_INCAR_MODE.push("EXTERNAL");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_FILE]")) {
      vflags.KBIN_VASP_KPOINTS_FILE.push("KEYWORD");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_MODE_EXPLICIT]")) {
      vflags.KBIN_VASP_KPOINTS_MODE.push("EXPLICIT");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_MODE_EXPLICIT]STOP")) {
      vflags.KBIN_VASP_KPOINTS_MODE.push("EXPLICIT_START_STOP");
    }
    if (vflags.KBIN_VASP_KPOINTS_FILE.flag("KEYWORD")) { // ME20181113
      vflags.KBIN_VASP_KPOINTS_EXPLICIT.str(aurostd::substring2string(AflowIn, "[VASP_KPOINTS_FILE]", 0));
    }
    if (vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")) { // ME20181113
      vflags.KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(aurostd::substring2string(AflowIn, "[VASP_KPOINTS_MODE_EXPLICIT]START", "[VASP_KPOINTS_MODE_EXPLICIT]STOP", 1));
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_MODE_IMPLICIT]")) {
      vflags.KBIN_VASP_KPOINTS_MODE.push("IMPLICIT");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_MODE_EXTERNAL]")) {
      vflags.KBIN_VASP_KPOINTS_MODE.push("EXTERNAL");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_FILE]FILE=", true)) { // ME20181113
      const string file = aurostd::substring2string(AflowIn, "[VASP_KPOINTS_FILE]FILE=", 1, true);
      vflags.KBIN_VASP_KPOINTS_FILE.push_attached("FILE", file);
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_KPOINTS_FILE]COMMAND=", true)) { // ME20181113
      const string command = aurostd::substring2string(AflowIn, "[VASP_KPOINTS_FILE]COMMAND=", 1, true);
      vflags.KBIN_VASP_KPOINTS_FILE.push_attached("COMMAND", command);
    }

    // KPOINTS FOR RELAX
    vflags.KBIN_VASP_KPOINTS_KMODE.options2entry(AflowIn, "[VASP_KPOINTS_FILE]KMODE=", false, vflags.KBIN_VASP_KPOINTS_KMODE.xscheme);
    vflags.KBIN_VASP_KPOINTS_KPPRA.options2entry(AflowIn, "[VASP_KPOINTS_FILE]KPPRA=", false, vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme);
    if (vflags.KBIN_VASP_KPOINTS_KPPRA.isentry == false) {
      vflags.KBIN_VASP_KPOINTS_KPPRA.clear();
      vflags.KBIN_VASP_KPOINTS_KPPRA.isentry = true;
      vflags.KBIN_VASP_KPOINTS_KPPRA.push("100");
    }
    vflags.KBIN_VASP_KPOINTS_KSCHEME.options2entry(AflowIn, "[VASP_KPOINTS_FILE]KSCHEME=", false, vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme);
    if (vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry == false) {
      vflags.KBIN_VASP_KPOINTS_KSCHEME.clear();
      vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry = true;
      vflags.KBIN_VASP_KPOINTS_KSCHEME.push(DEFAULT_KSCHEME);
    }
    vflags.KBIN_VASP_KPOINTS_KSHIFT.options2entry(AflowIn, "[VASP_KPOINTS_FILE]KSHIFT=", false, vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme);

    // KPOINTS FOR STATIC
    vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.options2entry(AflowIn, "[VASP_KPOINTS_FILE]STATIC_KMODE=", false, vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme);
    vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.options2entry(AflowIn, "[VASP_KPOINTS_FILE]STATIC_KPPRA=", false, vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme);
    vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.options2entry(AflowIn, "[VASP_KPOINTS_FILE]STATIC_KSCHEME=", false, vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme);
    vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.options2entry(AflowIn, "[VASP_KPOINTS_FILE]STATIC_KSHIFT=", false, vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme);

    // KPOINTS FOR DIELECTRIC
    vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KMODE.options2entry(AflowIn, "[VASP_KPOINTS_FILE]DIELECTRIC_KMODE=", false, vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KMODE.xscheme);
    vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KPPRA.options2entry(AflowIn, "[VASP_KPOINTS_FILE]DIELECTRIC_KPPRA=", false, vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KPPRA.xscheme);
    vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KSCHEME.options2entry(AflowIn, "[VASP_KPOINTS_FILE]DIELECTRIC_KSCHEME=", false, vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KSCHEME.xscheme);
    vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KSHIFT.options2entry(AflowIn, "[VASP_KPOINTS_FILE]DIELECTRIC_KSHIFT=", false, vflags.KBIN_VASP_KPOINTS_DIELECTRIC_KSHIFT.xscheme);

    vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.options2entry(AflowIn, "[VASP_KPOINTS_FILE]BANDS_LATTICE=", false, vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.xscheme); // scheme already loaded in aflow_xclasses.cpp is AUTO
    if (!vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry) {
      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();
      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(DEFAULT_BANDS_LATTICE);
    }

    if ((vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) && !vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry) {
      cerr << R"(WARNING: if you use vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") or vflags.KBIN_VASP_RUN.flag("STATIC_BANDS"), you must specify KBIN_VASP_KPOINTS_BANDS_LATTICE)" << endl;
      cerr << "         Taking defauls vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
    }

    vflags.KBIN_VASP_KPOINTS_BANDS_GRID.options2entry(AflowIn, "[VASP_KPOINTS_FILE]BANDS_GRID=", false, vflags.KBIN_VASP_KPOINTS_BANDS_GRID.xscheme); // scheme already loaded in aflow_xclasses.cpp is 20
    if (!vflags.KBIN_VASP_KPOINTS_BANDS_GRID.isentry) { // CO20210805
      vflags.KBIN_VASP_KPOINTS_BANDS_GRID.clear();
      vflags.KBIN_VASP_KPOINTS_BANDS_GRID.push(aurostd::utype2string(DEFAULT_BANDS_GRID));
    }
    if ((vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("REPEAT_BANDS"))
        && !vflags.KBIN_VASP_KPOINTS_BANDS_GRID.isentry) {
      message = "";
      message += "if you run RELAX_STATIC_BANDS, STATIC_BANDS, REPEAT_STATIC_BANDS, or REPEAT_BANDS, you must specify KBIN_VASP_KPOINTS_BANDS_GRID";
      message += "taking default KBIN_VASP_KPOINTS_BANDS_GRID=" + vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_string;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, _LOGGER_WARNING_);
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]")) {
      vflags.KBIN_VASP_POSCAR_FILE.push("KEYWORD");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_MODE_EXPLICIT]")) {
      vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT");
    }

    if (aurostd::substring2bool(AflowIn, _VASP_POSCAR_MODE_EXPLICIT_START_) && aurostd::substring2bool(AflowIn, _VASP_POSCAR_MODE_EXPLICIT_STOP_)) {
      vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP");
    }

    if ((aurostd::substring2bool(AflowIn, _VASP_POSCAR_MODE_EXPLICIT_START_P_) && aurostd::substring2bool(AflowIn, _VASP_POSCAR_MODE_EXPLICIT_STOP_P_))) {
      vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP_POINT"); // CO20200624
    }

    if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") && !kflags.KBIN_FROZSL) { // NO FROZSL
      if (LDEBUG) {
        cerr << "DEBUG: vflags.KBIN_VASP_POSCAR_MODE.flag(\"EXPLICIT_START_STOP_POINT\")=" << vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") << endl;
      }
      if (LDEBUG) {
        cerr << "DEBUG: kflags.KBIN_PHONONS_CALCULATION_FROZSL=" << kflags.KBIN_PHONONS_CALCULATION_FROZSL << endl;
      }
      if (LDEBUG) {
        cerr << "DEBUG: kflags.KBIN_FROZSL_DOWNLOAD=" << kflags.KBIN_FROZSL_DOWNLOAD << endl;
      }
      if (LDEBUG) {
        cerr << "DEBUG: kflags.KBIN_FROZSL_FILE=" << kflags.KBIN_FROZSL_FILE << endl;
      }
      stringstream input_file;
      input_file.clear();
      // loading
      if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        input_file.str(AflowIn);
      }
      if (kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
        FROZSL::Setup_frozsl_init_input(AflowIn, FileMESSAGE, input_file, aflags, kflags);
        FROZSL::Extract_INPUT(AflowIn, FileMESSAGE, input_file, aflags, kflags);
      }
      if (kflags.KBIN_FROZSL_DOWNLOAD) {
        FROZSL::Setup_frozsl_init_input(AflowIn, FileMESSAGE, input_file, aflags, kflags);
      }
      if (kflags.KBIN_FROZSL_FILE) {
        FROZSL::File_INPUT(AflowIn, FileMESSAGE, input_file, aflags, kflags);
      }

      vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP_POINT");
      // done loading now load structures up
      vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP", false); // some default
      aurostd::substring2strings(input_file.str(), vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING, _VASP_POSCAR_MODE_EXPLICIT_START_P_); // CO20200624
      // some verbose
      for (size_t i = 0; i < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size(); i++) {
        if (LDEBUG) {
          cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i] << endl;
        }
      }
      // load up the structures
      for (size_t i = 0; i < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size(); i++) {
        const string START = _VASP_POSCAR_MODE_EXPLICIT_START_P_ + vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i]; // CO20200624
        const string STOP = _VASP_POSCAR_MODE_EXPLICIT_STOP_P_ + vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i]; // CO20200624
        stringstream POSCAR;
        POSCAR.str(aurostd::substring2string(input_file.str(), START, STOP, -1));
        if (!POSCAR.str().empty()) {
          vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.emplace_back(POSCAR, IOVASP_AUTO);
        }
      }
      if (LDEBUG) {
        cerr << "DEBUG " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << endl;
      }
      if (LDEBUG) {
        cerr << "DEBUG " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size() << endl;
      }
      if (vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() != vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size()) {
        message = "IN " + _AFLOWIN_ + " in Directory=" + aflags.Directory + '\n';
        message += "vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()=" + aurostd::utype2string<uint>(vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) + '\n';
        message += "vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size()=" + aurostd::utype2string<uint>(vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size());
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      for (size_t i = 0; i < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size(); i++) {
        if (LDEBUG) {
          cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(i) << endl;
        }
      }
    } else {
      vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.clear();
      vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.clear();
    }
    // the rest for POSCAR
    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_MODE_IMPLICIT]")) {
      vflags.KBIN_VASP_POSCAR_MODE.push("IMPLICIT");
    }

    // vflags.KBIN_VASP_POSCAR_FILE_SYSTEM_AUTO                =   aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]SYSTEM_AUTO",true);
    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]PROTOTYPE=", true)) {
      vflags.KBIN_VASP_POSCAR_FILE.push("PROTOTYPE");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_MODE_EXTERNAL]")) {
      vflags.KBIN_VASP_POSCAR_MODE.push("EXTERNAL");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]FILE=", true)) {
      vflags.KBIN_VASP_POSCAR_FILE.push("FILE");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]COMMAND=", true)) {
      vflags.KBIN_VASP_POSCAR_FILE.push("COMMAND");
    }

    // VOLUMES
    vflags.KBIN_VASP_POSCAR_FILE_VOLUME.clear();
    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]VOLUME=", true)) {
      vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("EQUAL_EQUAL");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]VOLUME+=", true)) {
      vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("PLUS_EQUAL");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POSCAR_FILE]VOLUME*=", true)) {
      vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("MULTIPLY_EQUAL");
    }
    if (!vflags.KBIN_VASP_POSCAR_FILE_VOLUME.xscheme.empty()) {
      vflags.KBIN_VASP_POSCAR_FILE_VOLUME.isentry = true;
    }

    // CONVERT_UNIT_CELL stuff
    BflowIn = AflowIn;
    aurostd::StringSubstInPlace(BflowIn, "=", "_");
    aurostd::StringSubstInPlace(BflowIn, "CONVERT_UNIT_CELL_", "CONVERT_UNIT_CELL="); // bypass for getting all "_"
    vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.options2entry(BflowIn, string(_STROPT_ + "CONVERT_UNIT_CELL="), aurostd_xoptionMULTI, ""); // stack them all
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " BEFORE vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme << endl; // ME20181113
    }
    // // PRIORITIES
    if (vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") || vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI", false);
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL", false);
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT", false);
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ", false);
    } // some PRIORITIES
    if (vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL", false);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " AFTER vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme << endl; // ME20181113
    }

    // DEBUG
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"MINKOWSKI\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;
    }

    // VOLUMES
    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.clear();

    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn, "EQUAL_EQUAL", _STROPT_ + "VOLUME=", "");
    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn, "PLUS_EQUAL", _STROPT_ + "VOLUME+=", "");
    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn, "MULTIPLY_EQUAL", _STROPT_ + "VOLUME*=", "");
    if (LDEBUG) {
      cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("EQUAL_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"PLUS_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"PLUS_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("PLUS_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    }
    if (LDEBUG) {
      cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("MULTIPLY_EQUAL") << endl;
    }

    if (!vflags.KBIN_VASP_FORCE_OPTION_VOLUME.xscheme.empty()) {
      vflags.KBIN_VASP_FORCE_OPTION_VOLUME.isentry = true;
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]")) {
      vflags.KBIN_VASP_POTCAR_FILE.push("KEYWORD");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]SYSTEM_AUTO", true)) {
      vflags.KBIN_VASP_POTCAR_FILE.push("SYSTEM_AUTO");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]PREFIX=", true)) {
      vflags.KBIN_VASP_POTCAR_FILE.push("PREFIX");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]SUFFIX=", true)) {
      vflags.KBIN_VASP_POTCAR_FILE.push("SUFFIX");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]FILE=", true)) {
      vflags.KBIN_VASP_POTCAR_FILE.push("FILE");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]COMMAND=", true)) {
      vflags.KBIN_VASP_POTCAR_FILE.push("COMMAND");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_MODE_EXPLICIT]")) {
      vflags.KBIN_VASP_POTCAR_MODE.push("EXPLICIT");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_MODE_IMPLICIT]")) {
      vflags.KBIN_VASP_POTCAR_MODE.push("IMPLICIT");
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_MODE_EXTERNAL]")) {
      vflags.KBIN_VASP_POTCAR_MODE.push("EXTERNAL");
    }

    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]PREFIX=", true)) { // CO20181113
      const string file = aurostd::substring2string(AflowIn, "[VASP_POTCAR_FILE]PREFIX=", 1, true);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("PREFIX", file);
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]SUFFIX=", true)) { // CO20181113
      const string file = aurostd::substring2string(AflowIn, "[VASP_POTCAR_FILE]SUFFIX=", 1, true);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("SUFFIX", file);
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]FILE=", true)) { // CO20181113
      const string file = aurostd::substring2string(AflowIn, "[VASP_POTCAR_FILE]FILE=", 1, true);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("FILE", file);
    }
    if (aurostd::substring2bool(AflowIn, "[VASP_POTCAR_FILE]COMMAND=", true)) { // CO20181113
      const string command = aurostd::substring2string(AflowIn, "[VASP_POTCAR_FILE]COMMAND=", 1, true);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("COMMAND", command);
    }

    if (vflags.KBIN_VASP_POTCAR_FILE.flag("KEYWORD")) { // CO20181113
      vflags.KBIN_VASP_POTCAR_EXPLICIT.str(aurostd::substring2string(AflowIn, "[VASP_POTCAR_FILE]", 0));
    }

    // APL ENTRIES
    if (LDEBUG) {
      cerr << "DEBUG: " << __AFLOW_FUNC__ << " (APL)" << endl;
    }

    // CO20170601 START
    // CO make backwards and forwards compatible with all possible workflows
    vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.options2entry(AflowIn, "[AFLOW_APL]KPPRA=|[AFLOW_QHA]KPPRA=|[AFLOW_AAPL]KPPRA=|[AFLOW_PHONONS]KPPRA=", false,
                                                         vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
    vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.options2entry(AflowIn, "[AFLOW_APL]KSCHEME=|[AFLOW_QHA]KSCHEME=|[AFLOW_AAPL]KSCHEME=|[AFLOW_PHONONS]KSCHEME=", false,
                                                           vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_SCHEME"

    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.options2entry(AflowIn, "[AFLOW_APL]KPOINTS=|[AFLOW_QHA]KPOINTS=|[AFLOW_AAPL]KPOINTS=", false, vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.xscheme);
    // ME202020427 - APL k-point handling needs to be moved to modules eventually
    //  This is a non-standard feature and should not be defaulted
    vflags.KBIN_VASP_KPOINTS_PHONONS_GRID.options2entry(AflowIn, "[AFLOW_APL]KPOINTS_GRID=|[AFLOW_QHA]KPOINTS_GRID=|[AFLOW_AAPL]KPOINTS_GRID=", false, "");
    vflags.KBIN_VASP_KPOINTS_PHONONS_SHIFT.options2entry(AflowIn, "[AFLOW_APL]KPOINTS_SHIFT=|[AFLOW_QHA]KPOINTS_SHIFT=|[AFLOW_AAPL]KPOINTS_SHIFT=", false, "");
    // CO20170601 END

    // FROZSL ENTRIES
    if (LDEBUG) {
      cerr << "DEBUG: " << __AFLOW_FUNC__ << " (FROZSL)" << endl;
    }

    if (LDEBUG) {
      cerr << "DEBUG: " << __AFLOW_FUNC__ << " (STOP)" << endl;
    }

    return vflags;
  }
} // namespace KBIN

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
namespace KBIN {
  bool VASP_ExtractNGF(string OUTCAR, int& NGXF, int& NGYF, int& NGZF);
} // namespace KBIN

namespace KBIN {
  /// @brief processes and runs the aflow routine on the files in the active directory
  bool VASP_Directory(ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    const bool LDEBUG = (false || _DEBUG_KVASP_ || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << "DEBUG: KBIN::VASP_Directory (BEGIN)" << endl;
    }
    ostringstream aus;
    const string::iterator pos;

    ifstream FileAFLOWIN;
    string FileNameAFLOWIN;
    string AflowIn;
    FileNameAFLOWIN = aflags.Directory + "/" + _AFLOWIN_;
    FileAFLOWIN.open(FileNameAFLOWIN.c_str(), std::ios::in);
    FileAFLOWIN.clear();
    FileAFLOWIN.seekg(0);
    AflowIn = "";
    char c;
    while (FileAFLOWIN.get(c)) {
      AflowIn += c; // READ _AFLOWIN_ and put into AflowIn
    }
    FileAFLOWIN.clear();
    FileAFLOWIN.seekg(0);
    AflowIn = aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    if (!FileAFLOWIN) { // ******* _AFLOWIN_ does not exist
      aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(aus, XHOST.QUIET);
      return false;
    }
    aflags.QUIET = false;
    _vflags vflags;
    vflags = KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn, FileMESSAGE, aflags, kflags);

    // *********************************************************************************************************************
    // OPERATIONS related to PARTICULAR MACHINES ***************************************************************************

    if (LDEBUG) {
      cerr << "[DEBUG] aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL.getattachedscheme("NAME") << endl;
    }
    if (LDEBUG) {
      cerr << "[DEBUG] aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << endl; // HE20220309 use machine name
    }

    // ***************************************************************************
    if (LDEBUG) {
      cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") << endl; // CO20210315
    }
    if (LDEBUG) {
      cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") << endl;
    }
    if (LDEBUG) {
      cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_DIELECTRIC\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_DIELECTRIC") << endl;
    }
    if (LDEBUG) {
      cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS") << endl;
    }
    if (LDEBUG) {
      cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_DIELECTRIC\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_DIELECTRIC") << endl;
    }
    if (LDEBUG) {
      cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_DELSOL\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL") << endl;
    }

    // ***************************************************************************
    // Get the KBIN_BIN name
    aurostd::StringstreamClean(aus);
    aus << "00000  MESSAGE KBIN::VASP_Directory Running KBIN_BIN=\"" << kflags.KBIN_BIN << "\"" << Message(__AFLOW_FILE__, aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    // ***************************************************************************
    // Some verbose
    if (kflags.KBIN_POCC) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running POCC_CALCULATION" << Message(__AFLOW_FILE__, aflags) << endl;
    } // CO20180419 //POCC is special needs to run first because there is NO poscar defined yet
    else if (kflags.KBIN_PHONONS_CALCULATION_APL) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_APL" << Message(__AFLOW_FILE__, aflags) << endl;
    } else if (kflags.KBIN_PHONONS_CALCULATION_QHA) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_QHA" << Message(__AFLOW_FILE__, aflags) << endl; // CO20170601
    } else if (kflags.KBIN_PHONONS_CALCULATION_AAPL) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AAPL" << Message(__AFLOW_FILE__, aflags) << endl; // CO20170601
    } else if (kflags.KBIN_PHONONS_CALCULATION_AGL) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AGL (Debye Model)" << Message(__AFLOW_FILE__, aflags) << endl;
    } else if (kflags.KBIN_PHONONS_CALCULATION_AEL) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AEL (Elastic constants)" << Message(__AFLOW_FILE__, aflags) << endl;
    } else if (kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
      aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_FROZSL" << Message(__AFLOW_FILE__, aflags) << endl;
    } else {
      if (vflags.KBIN_VASP_RUN.flag("GENERATE")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_GENERATE" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("STATIC")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("KPOINTS")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_KPOINTS" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("RELAX")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC_BANDS" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("STATIC_DIELECTRIC")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC_DIELECTRIC" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC_BANDS" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_DIELECTRIC")) {
        aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC_DIELECTRIC" << Message(__AFLOW_FILE__, aflags) << endl;
      }
    }
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    // ***************************************************************************
    size_t ntasks = 0;
    ntasks = 1; // default
    if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
      ntasks = vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() + 1; // CO20200624 - include head directory as well (at the end)
      aus << "00000  MESSAGE Loaded ntasks = " << ntasks << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      size_t i = 0;
      for (i = 0; i < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size(); i++) {
        aus << "00000  MESSAGE task " << i + 1 << "/" << ntasks << " in subdirectory " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING[i] << endl; // CO20200624 - +1
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      aus << "00000  MESSAGE task " << ((i++) + 1) << "/" << ntasks << " in main directory " << aflags.Directory << endl; // CO20200624 - include head directory as well (at the end)
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    // ***************************************************************************
    // start the loop !
    _aflags aflags_backup;
    aflags_backup = aflags;
    _kflags kflags_backup;
    kflags_backup = kflags;

    for (size_t ixvasp = 0; ixvasp < ntasks; ixvasp++) { // LOOP ixvasp
      // declarations
      _xvasp xvasp;
      xvasp.clear();
      xvasp.POSCAR_index = ixvasp;
      aflags = aflags_backup;
      kflags = kflags_backup; // load it up
      readModulesFromAflowIn(AflowIn, kflags, xvasp); // ME20181027
      // some verbose
      if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") && ixvasp < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) { // CO20200624 - include head directory as well (at the end)
        aus << "00000  MESSAGE START loop " << xvasp.POSCAR_index + 1 << "/" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << Message(__AFLOW_FILE__, aflags) << endl; // CO20200624 - +1
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      if (LDEBUG) {
        cerr << XPID << "KBIN::VASP_Directory: [1]" << xvasp.str << endl;
      }
      // ------------------------------------------
      // now start for each xvasp
      xvasp.Directory = aflags.Directory;
      if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") && ixvasp < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) { // CO20200624 - include head directory as well (at the end)
        xvasp.Directory = aflags.Directory + "/" + KBIN_SUBDIRECTORIES + vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(xvasp.POSCAR_index);
        aus << "00000  MESSAGE Taking loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") && ixvasp < vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) { // CO20200624 - include head directory as well (at the end)
        if (aurostd::FileExist(xvasp.Directory)) {
          aus << "00000  MESSAGE Skipping loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          return false; // avoid rerunning
        } else {
          // before making it, check it again... NFS problem... check LOCK again
          if (aurostd::FileExist(xvasp.Directory + "/" + _AFLOWLOCK_)
              || aurostd::CompressFileExist(xvasp.Directory + "/" + _AFLOWLOCK_)
              || aurostd::FileExist(xvasp.Directory + "/LLOCK")
              || aurostd::CompressFileExist(xvasp.Directory + "/LLOCK")) {
            return false; // to fight against NFS cache
          }
          aurostd::DirectoryMake(xvasp.Directory);
          aus << "00000  MESSAGE Creating loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          aurostd::string2file("NNNNN  KBIN LLOCK ASAP for NFS concurrent jobs (aflow" + std::string(AFLOW_VERSION) + ")", xvasp.Directory + "/LLOCK");
        }
      }

      aflags.Directory = xvasp.Directory; // so we are set ! since there are plenty of routines with aflags.Directory inside
      aus << "00000  MESSAGE Performing loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      // ------------------------------------------
      // do the flags
      if (LDEBUG) {
        cerr << XPID << "KBIN::VASP_Directory: [2]" << xvasp.str << endl;
      }
      vflags.KBIN_VASP_INCAR_VERBOSE = true; // ALWAYS

      if (vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option) { // CO20211217 - prevents VASP_Backup() from deleting and recycles for next run (VASP_RecycleExtraFile())
        xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED", vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option); // passing flag to xvasp
        aus << "00000  MESSAGE Saving WAVECAR files" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }

      // CO, come back here at some point
      // think about taking pocc structure and reducing first if possible (and requested)
      // need to first convert to non-pocc structure (already programmed in aflow_pocc.cpp)
      // plug this in as xvasp.str (also create xvasp.str_pocc then)
      // reduce as requested and re-pocc the structure
      // think about if we need a separate flag for reducing pocc vs. reducing derivative structures
      //  produce BEFORE NOMIX
      if (!(kflags.KBIN_POCC || kflags.KBIN_PHONONS_CALCULATION_FROZSL)) { // CO20180419, do NOT produce POSCAR for POCC
        if (!KBIN::VASP_Produce_INPUT(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags)) {
          return false;
        }
        if (!KBIN::VASP_Modify_INPUT(xvasp, FileMESSAGE, aflags, kflags, vflags)) {
          return false;
        }
        if (kflags.KBIN_QSUB && !KBIN::QSUB_Extract(xvasp.xqsub, AflowIn, FileAFLOWIN, FileMESSAGE, aflags, kflags)) {
          return false;
        }
        if (kflags.KBIN_QSUB_MODE1 && !KBIN::QSUB_Extract_Mode1(xvasp.xqsub, FileMESSAGE, aflags, kflags)) {
          return false;
        }
        if (kflags.KBIN_QSUB_MODE2 && !KBIN::QSUB_Extract_Mode2(xvasp.xqsub, FileMESSAGE, aflags, kflags)) {
          return false;
        }
        if (kflags.KBIN_QSUB_MODE3 && !KBIN::QSUB_Extract_Mode3(xvasp.xqsub, FileMESSAGE, aflags, kflags)) {
          return false;
        }
      }
      if (!XHOST.GENERATE_AFLOWIN_ONLY && vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.isentry) { // ME20210709 - For now, skip check if generate_aflowin_only. In the future, the elements (not the pseudopotentials) need to be grabbed from the aflow.in
        const string potentials = xvasp.POTCAR_POTENTIALS.str();
        if (!aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/1/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/2/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/3/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/58/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/59/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/60/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/115/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/116/")
            && !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory + "/"), "/117/")) {
          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SKIP_NOMIX (NEGLECT_NOMIX, NEGLECT_IMMISCIBLE)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (MiscibilityCheck(potentials) == MISCIBILITY_SYSTEM_NOMIX) {
            aus << "00000  MESSAGE Skipping system: " << aurostd::VASP_PseudoPotential_CleanName(potentials) << " is known to be immiscible (aflow_nomix.cpp)" << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            stringstream command("");
            command << "cat " << xvasp.Directory << "/" << _AFLOWLOCK_ << " > " << xvasp.Directory << "/" << DEFAULT_AFLOW_IMMISCIBILITY_OUT << endl;
            aurostd::execute(command);
            return false;
          }
        }
        if (MiscibilityCheck(potentials) == MISCIBILITY_SYSTEM_MISCIBLE) {
          aus << "00000  MESSAGE Running system: " << aurostd::VASP_PseudoPotential_CleanName(potentials) << " is known to be miscible (aflow_nomix.cpp)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (MiscibilityCheck(potentials) == MISCIBILITY_SYSTEM_UNKNOWN) {
          aus << "00000  MESSAGE Running system: " << aurostd::VASP_PseudoPotential_CleanName(potentials) << " is unknown (aflow_nomix.cpp)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
      }

      // ***************************************************************************
      // READY TO RUN
      if (LDEBUG) {
        cerr << XPID << "KBIN::VASP_Directory: [3]" << endl;
      }
      if (LDEBUG) {
        cerr << xvasp.str << endl;
      }
      xvasp.NRELAX = 0;
      ostringstream aus;
      const bool PAWGGA2 = false;
      // ***************************************************************************
      // directory check
      ifstream DirectoryStream;
      DirectoryStream.open(xvasp.Directory.c_str(), std::ios::in);
      if (!DirectoryStream) {
        aus << "XXXXX  MAKING DIRECTORY = " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(aus, XHOST.QUIET);
        aurostd::DirectoryMake(xvasp.Directory);
      }
      // ***************************************************************************
      // DO THE SYMMETRY NEIGHBORS CALCULATION
      // DX

      if (!(kflags.KBIN_POCC || kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_FROZSL)
          || aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN) { // CO, do internally
        // DX
        if (aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN) {
          return KBIN_StepSymmetryPerform(xvasp.str, AflowIn, FileMESSAGE, aflags, kflags, true, cout);
        }
        if (!KBIN_StepSymmetryPerform(xvasp.str, AflowIn, FileMESSAGE, aflags, kflags, true, cout)) {
          return false;
        } // DO THE SYMMETRY CALCULATION
        // DX
      }
      // VASP VASP WRITE
      // ***************************************************************************
      // VASP INPUT FILES ARE DONE, NOW WE CAN USE OR MODYFYING THEM
      if (vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NOTUNE, no tuning xCARs - ";
        aus << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - ";
        aus << XHOST.hostname << " - " << aflow_get_time_string() << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      // ***************************************************************************
      // VASP HOW TO RUN ??
      // GENERATE ONLY -------------------------------------------------------------
      if (vflags.KBIN_VASP_RUN.flag("GENERATE")) {
        KBIN::VASP_Write_INPUT(xvasp, vflags); // VASP VASP WRITE
        aus << "00000  MESSAGE VASP generation files ONLY" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = 0;
        return false;
      }
      // CO20180420 - added if/else-if for workflows that need to PRECEDE relax/static/etc.
      //  RUN SOMETHING
      if (kflags.KBIN_POCC) { // RUN POCC ------------------------  //CO20180419 //POCC is special, run as priority
        aus << "00000  MESSAGE PERFORMING POCC_CALCULATION" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
      } else if (kflags.KBIN_PHONONS_CALCULATION_APL) { // RUN PHONONS APL ------------------------
        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_APL" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
        // CO20170601 START
      } else if (kflags.KBIN_PHONONS_CALCULATION_QHA) { // RUN PHONONS QHA ------------------------
        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_QHA" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
      } else if (kflags.KBIN_PHONONS_CALCULATION_AAPL) { // RUN PHONONS AAPL ------------------------
        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AAPL" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
        // CO20170601 END
      } else if (kflags.KBIN_PHONONS_CALCULATION_AGL) { // RUN PHONONS AGL ------------------------
        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AGL (Debye Model)" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
      } else if (kflags.KBIN_PHONONS_CALCULATION_AEL) { // RUN PHONONS AEL ------------------------
        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AEL (Elastic constants)" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
      } else if (kflags.KBIN_PHONONS_CALCULATION_FROZSL) { // RUN PHONONS FROZSL ------------------------
        aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_FROZSL" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xvasp.NRELAX = -3;
      } else {
        if (vflags.KBIN_VASP_RUN.flag("STATIC")) { // RUN STATIC ------------------------
          aus << "00000  MESSAGE Performing Static RUN" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          xvasp.NRELAX = -1;
        }
        if (vflags.KBIN_VASP_RUN.flag("KPOINTS")) { // RUN KPOINTS ------------------------
          aus << "00000  MESSAGE Running KPOINTS swap" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          xvasp.NRELAX = -2;
        }
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) { // RUN RELAX_STATIC_BANDS ------------------------
          xvasp.NRELAX = vflags.KBIN_VASP_RUN_NRELAX;
          if (xvasp.NRELAX <= 0) {
            aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintErrorStream(FileMESSAGE, aus, XHOST.QUIET);
            xvasp.NRELAX = 0;
            return false;
          }
          aus << "00000  MESSAGE RELAX_STATIC_BANDS Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_DIELECTRIC")) { // RUN RELAX_STATIC_DIELECTRIC ------------------------
          xvasp.NRELAX = vflags.KBIN_VASP_RUN_NRELAX;
          if (xvasp.NRELAX <= 0) {
            aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintErrorStream(FileMESSAGE, aus, XHOST.QUIET);
            xvasp.NRELAX = 0;
            return false;
          }
          aus << "00000  MESSAGE RELAX_STATIC_DIELECTRIC Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) { // RUN RELAX_STATIC ------------------------
          xvasp.NRELAX = vflags.KBIN_VASP_RUN_NRELAX;
          if (xvasp.NRELAX <= 0) {
            aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintErrorStream(FileMESSAGE, aus, XHOST.QUIET);
            xvasp.NRELAX = 0;
            return false;
          }
          aus << "00000  MESSAGE RELAX_STATIC Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) { // RUN STATIC_BANDS ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE STATIC_BANDS Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_RUN.flag("STATIC_DIELECTRIC")) { // RUN STATIC_DIELECTRIC ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE STATIC_DIELECTRIC Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_RUN.flag("RELAX")) { // RUN RELAX ------------------------
          xvasp.NRELAX = vflags.KBIN_VASP_RUN_NRELAX;
          if (xvasp.NRELAX <= 0) {
            aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintErrorStream(FileMESSAGE, aus, XHOST.QUIET);
            xvasp.NRELAX = 0;
            return false;
          }
          aus << "00000  MESSAGE RELAX Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) { // RUN REPEAT_STATIC ------------------------ //CO20210315
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE REPEAT_STATIC Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) { // RUN REPEAT_STATIC_BANDS ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE REPEAT_STATIC_BANDS Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_DIELECTRIC")) { // RUN REPEAT_STATIC_DIELECTRIC ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE REPEAT_STATIC_DIELECTRIC Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) { // RUN REPEAT_BANDS ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE REPEAT_BANDS Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_DIELECTRIC")) { // RUN REPEAT_DIELECTRIC ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE REPEAT_DIELECTRIC Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) { // RUN REPEAT_DELSOL ------------------------
          xvasp.NRELAX = -1;
          aus << "00000  MESSAGE REPEAT_DELSOL Running " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
      }

      // ***************************************************************************
      // READY TO RUN
      // ***************************************************************************
      // START
      if (LDEBUG) {
        cerr << XPID << "KBIN::VASP_Directory: [4]" << xvasp.str << endl;
      }
      // ***************************************************************************
      // FIX BLANK SPECIES
      if (!xvasp.str.species.empty()) {
        // CO, fixing for RHT routines
        if (xvasp.str.species.at(0).empty()) {
          pflow::fixEmptyAtomNames(xvasp.str);
        }
      }
      // ***************************************************************************
      // PRESCRIPT
      if (kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) {
        KBIN::RUN_DirectoryScript(aflags, DEFAULT_AFLOW_PRESCRIPT_COMMAND, DEFAULT_AFLOW_PRESCRIPT_OUT);
      }
      // ***************************************************************************
      // CO20180419 - POCC always comes first (NO POSCAR), need to convert PARTCAR -> POSCARs
      // other workflows follow, all of these precede relaxation/static/etc.
      if (kflags.KBIN_POCC) {
        KBIN::VASP_RunPOCC(xvasp, AflowIn, aflags, kflags, vflags, FileMESSAGE);
      } // CO20180419
      else if (kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL) {
        // ME20200107 - Wrap in a try statement so that faulty APL runs don't kill other post-processing
        try {
          KBIN::VASP_RunPhonons_APL(xvasp, AflowIn, aflags, kflags, vflags, FileMESSAGE); // PHONONIC PHONONIC PHONONIC //CO20170601
        } catch (aurostd::xerror e) {
          pflow::logger(e.whereFileName(), e.whereFunction(), e.buildMessageString(), aflags.Directory, FileMESSAGE, std::cout, _LOGGER_ERROR_);
        }
      } else if (kflags.KBIN_PHONONS_CALCULATION_AGL) {
        KBIN::VASP_RunPhonons_AGL(xvasp, AflowIn, aflags, kflags, vflags, FileMESSAGE);
      } else if (kflags.KBIN_PHONONS_CALCULATION_AEL) {
        KBIN::VASP_RunPhonons_AEL(xvasp, AflowIn, aflags, kflags, vflags, FileMESSAGE);
      } else if (kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
        KBIN::VASP_RunPhonons_FROZSL(xvasp, AflowIn, aflags, kflags, vflags, FileMESSAGE);
      } else {
        if (LDEBUG) {
          cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species.size()=" << xvasp.str.species.size() << endl;
        }
        if (LDEBUG) {
          for (size_t i = 0; i < xvasp.str.species.size(); i++) {
            cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species[i]=[" << xvasp.str.species[i] << "]" << endl;
          }
        }
        if (LDEBUG) {
          cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species_pp.size()=" << xvasp.str.species_pp.size() << endl;
        }
        if (LDEBUG) {
          for (size_t i = 0; i < xvasp.str.species_pp.size(); i++) {
            cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species_pp[i]=[" << xvasp.str.species_pp[i] << "]" << endl;
          }
        }
        if (LDEBUG) {
          cerr << XPID << "KBIN::VASP_Directory: [6]" << xvasp.str << endl;
        }
        // --------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------
        // RELAX RELAX RELAX
        if (vflags.KBIN_VASP_RUN.flag("RELAX")) { // xvasp.RELAX>0
          KBIN::VASP_Write_INPUT(xvasp, vflags); // VASP VASP WRITE
          if (PAWGGA2) { // WAS A BUG IN PAW MAYBE IT IS FIXED
            // STEP 1
            aus << "11111  RELAXATION - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax2paw_gga", true, FileMESSAGE)) {
              KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [PAWGGA2 REL]");
              return false;
            }
            aus << "22222  END        - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          } else {
            if (xvasp.NRELAX == 0) {
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "STATIC RUN FIX INCAR: should not be here", _RUNTIME_ERROR_);
            } else { // DYNAMIC RUN
              for (xvasp.NRELAXING = 1; xvasp.NRELAXING <= xvasp.NRELAX; xvasp.NRELAXING++) {
                aus
                    << 11111 * xvasp.NRELAXING
                    << "  RELAXATION - "
                    << xvasp.Directory
                    << " - K=["
                    << xvasp.str.kpoints_k1
                    << " "
                    << xvasp.str.kpoints_k2
                    << " "
                    << xvasp.str.kpoints_k3
                    << "]"
                    << " - "
                    << kflags.KBIN_BIN
                    << Message(__AFLOW_FILE__, aflags)
                    << endl;
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
                if (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int > 0) {
                  KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp, xvasp.NRELAXING); // ADIABATIC
                }
                if (xvasp.NRELAXING < xvasp.NRELAX) {
                  if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string(xvasp.NRELAXING), "relax" + aurostd::utype2string(xvasp.NRELAXING), true, FileMESSAGE)) {
                    KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [RELAXATION<]");
                    return false;
                  }
                  KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp, aflags, vflags, xvasp.NRELAXING, true, FileMESSAGE); // check if it is the case of turning off spin //CO20210315 - always write_incar here
                  KBIN::XVASP_KPOINTS_IBZKPT_UPDATE(xvasp, aflags, vflags, xvasp.NRELAXING, true, FileMESSAGE); // check if it is the case of updating IBZKPT  //CO20210315 - always write_incar here
                  // ME20190301 BEGIN
                  //  CHGCAR/WAVECAR needs to be recycled if CHGCAR/WAVECAR=ON or VASP
                  //  won't be able to read the files. Bug found by Rico Friedrich
                  if (vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option) {
                    // ME20191031 - Only recycle when ICHARG was found
                    const string extra_incar = xvasp.AVASP_EXTRA_INCAR.str();
                    const size_t nlines = aurostd::GetNLinesString(extra_incar);
                    size_t l;
                    string line;
                    for (l = 1; l <= nlines; l++) {
                      line = aurostd::RemoveWhiteSpaces(aurostd::GetLineString(extra_incar, l));
                      if (aurostd::substring2bool(line, "ICHARG=1", true)) {
                        break;
                      }
                    }
                    if (l <= nlines) {
                      KBIN::VASP_RecycleExtraFile(xvasp, "CHGCAR", "relax" + aurostd::utype2string<int>(xvasp.NRELAXING));
                    }
                  }
                  // CO+MM20211217
                  // In general it's not a good idea to use the WAVECAR from relax1 to start relax2
                  // a big change in the unit cell size/shape is going drastically change the plane waves that should be used
                  // we can add another option to do this in the future
                  //[CO+MM20211217]if(vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option) KBIN::VASP_RecycleExtraFile(xvasp, "WAVECAR", "relax"+aurostd::utype2string<int>(xvasp.NRELAXING));
                  // ME20190301 END
                }
                if (xvasp.NRELAXING == xvasp.NRELAX) {
                  if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string(xvasp.NRELAXING), true, FileMESSAGE)) {
                    KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [RELAXATION=]");
                    return false;
                  }
                  KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp, aflags, vflags, xvasp.NRELAXING, false, FileMESSAGE); // ME20190610 - or else SPIN_REMOVE_RELAX_2 won't work  //CO20210315 - never write_incar here (last step for sure)
                }
                KBIN::XVASP_INCAR_ADJUST_ICHARG(xvasp, vflags, aflags, xvasp.NRELAXING, (xvasp.NRELAXING < xvasp.NRELAX), FileMESSAGE); // ME20191028 //CO20210315 - only write_incar if it's not the last relaxation (last step)
              }
              xvasp.NRELAXING = xvasp.NRELAX;
              xvasp.NRELAXING++;
              aus
                  << 11111 * xvasp.NRELAXING
                  << "  END        - "
                  << xvasp.Directory
                  << " - K=["
                  << xvasp.str.kpoints_k1
                  << " "
                  << xvasp.str.kpoints_k2
                  << " "
                  << xvasp.str.kpoints_k3
                  << "]"
                  << " - "
                  << kflags.KBIN_BIN
                  << Message(__AFLOW_FILE__, aflags)
                  << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
          }
        }
        // --------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------
        // --------------------------------------------------------------------------------------------------------------------
        // STATIC STATIC STATIC
        if (vflags.KBIN_VASP_RUN.flag("STATIC")) { // xvasp.RELAX=-1
          xvasp.aopts.flag("FLAG::POSCAR_PRESERVED", true); // in case of errors it is not lost but recycled
          KBIN::VASP_Write_INPUT(xvasp, vflags); // VASP VASP WRITE
          aus << 11111 << "  STATIC - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [STATIC]");
            return false;
          }
          if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error (OUTCAR_INCOMPLETE) [STATIC]");
            return false;
          }
          KBIN::VASP_Backup(xvasp, true, string("static"));
        }

        vector<double> xvasp_spin_evolution;
        xmatrix<double> rlattice(xvasp.str.lattice);

        string run_type_str;
        const std::vector<std::string> run_type_vec = {"RELAX_STATIC_BANDS",  "RELAX_STATIC_DIELECTRIC",  "RELAX_STATIC", "STATIC_BANDS",      "STATIC_DIELECTRIC", "STATIC", "KPOINTS", "RELAX", "REPEAT_STATIC",
                                                       "REPEAT_STATIC_BANDS", "REPEAT_STATIC_DIELECTRIC", "REPEAT_BANDS", "REPEAT_DIELECTRIC", "REPEAT_DELSOL"};
        for (auto& str : run_type_vec) {
          if (vflags.KBIN_VASP_RUN.flag(str) || vflags.KBIN_VASP_REPEAT.flag(str)) {
            run_type_str = str;
          }
        }
        aus << "00000  MESSAGE MODE= (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

        // DO THE RELAX PART (IF ANY)
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_DIELECTRIC") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {
          KBIN::VASP_Write_INPUT(xvasp, vflags); // VASP VASP WRITE
          for (xvasp.NRELAXING = 1; xvasp.NRELAXING <= xvasp.NRELAX; xvasp.NRELAXING++) {
            aus
                << 11111 * xvasp.NRELAXING
                << "  RELAXATION ("
                << run_type_str
                << ") - "
                << xvasp.Directory
                << " - K=["
                << xvasp.str.kpoints_k1
                << " "
                << xvasp.str.kpoints_k2
                << " "
                << xvasp.str.kpoints_k3
                << "]"
                << " - "
                << kflags.KBIN_BIN
                << Message(__AFLOW_FILE__, aflags)
                << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int > 0) {
              KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp, xvasp.NRELAXING); // ADIABATIC
            }
            if (xvasp.NRELAXING < xvasp.NRELAX) {
              if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string(xvasp.NRELAXING), "relax" + aurostd::utype2string(xvasp.NRELAXING), true, FileMESSAGE)) {
                KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " RELAXATION<]");
                return false;
              }
              // ME20190301 BEGIN
              //  CHGCAR/WAVECAR needs to be recycled if CHGCAR/WAVECAR=ON or VASP
              //  won't be able to read the files. Bug found by Rico Friedrich
              if (vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option) {
                // ME20191031 - Only recycle when ICHARG was found
                const string extra_incar = xvasp.AVASP_EXTRA_INCAR.str();
                const size_t nlines = aurostd::GetNLinesString(extra_incar);
                size_t l;
                string line;
                for (l = 1; l <= nlines; l++) {
                  line = aurostd::RemoveWhiteSpaces(aurostd::GetLineString(extra_incar, l));
                  if (aurostd::substring2bool(line, "ICHARG=1", true)) {
                    break;
                  }
                }
                if (l <= nlines) {
                  KBIN::VASP_RecycleExtraFile(xvasp, "CHGCAR", "relax" + aurostd::utype2string<int>(xvasp.NRELAXING));
                }
              }
              // CO+MM20211217
              // In general it's not a good idea to use the WAVECAR from relax1 to start relax2
              // a big change in the unit cell size/shape is going drastically change the plane waves that should be used
              // we can add another option to do this in the future
              //[CO+MM20211217]if(vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option) KBIN::VASP_RecycleExtraFile(xvasp, "WAVECAR", "relax"+aurostd::utype2string<int>(xvasp.NRELAXING));
              // ME20190301 END
            }
            if (xvasp.NRELAXING == xvasp.NRELAX && !KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string(xvasp.NRELAXING), true, FileMESSAGE)) {
              KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " RELAXATION=]");
              return false;
            }
            KBIN::XVASP_INCAR_ADJUST_ICHARG(xvasp, vflags, aflags, xvasp.NRELAXING, true, FileMESSAGE); // ME20191028 //CO20210315 - always write_incar, there's a STATIC that follows (at least)
            xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
            aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size() - 1) << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            //[CO20210315 - made robust enough to work in all cases]if(xvasp.NRELAXING<xvasp.NRELAX)
            KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp, aflags, vflags, xvasp.NRELAXING, true, FileMESSAGE); // check if it is the case of turning off spin  //CO20210315 - always write_incar, there's a STATIC that follows (at least)
          }
          if (xvasp.NRELAX > 0) {
            KBIN::VASP_Recycle(xvasp, "relax" + aurostd::utype2string(xvasp.NRELAX)); // bring back the stuff
          }
          //[CO20210315 - not sure why only if NRELAX==2, we could have NRELAX==1 and this still might apply]if(xvasp.NRELAX==2)
          KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp, aflags, vflags, xvasp.NRELAX, true,
                                              FileMESSAGE); // check if it is the case of turning off spin  //CO20210315 - always write_incar, there's a STATIC that follows (at least) //CO20210315 - no longer necessary per above (we check after every relaxation, but it doesn't hurt
          xvasp.NRELAXING = xvasp.NRELAX;
        }

        if (vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_DIELECTRIC")) {
          KBIN::VASP_Write_INPUT(xvasp, vflags); // VASP VASP WRITE
          aus << "00000  NO RELAXATION IN (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          xvasp.NRELAX = 0;
          xvasp.NRELAXING = xvasp.NRELAX;
        }

        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_DIELECTRIC")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_DIELECTRIC")) {
          // LOAD FORMER LOCK
          if (aurostd::FileExist(xvasp.Directory + "/" + run_type_str)) {
            stringstream lock_recycled;
            aurostd::file2stringstream(xvasp.Directory + "/" + run_type_str, lock_recycled);
            aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
            aus << "XXXXX FORMER LOCK BEGIN, recycled (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
            aus << "XXXXX FORMER LOCK END, recycled (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
            aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
          // clean up old files
          aurostd::RemoveFile(xvasp.Directory, std::regex(DEFAULT_AFLOW_END_OUT + ".*"));
          if (run_type_str.find("STATIC") != std::string::npos) {
            aurostd::RemoveFile(xvasp.Directory, std::regex(".*static.*"));
          }
          if (run_type_str.find("BANDS") != std::string::npos) {
            aurostd::RemoveFile(xvasp.Directory, std::regex(".*bands.*"));
          }
          if (run_type_str.find("DIELECTRIC") != std::string::npos) {
            aurostd::RemoveFile(xvasp.Directory, std::regex(".*dielectric.*"));
          }
          // unzip directory
          KBIN::DecompressDirectory(xvasp.Directory);
          // recycle files
          if (run_type_str.find("STATIC") != std::string::npos) {
            if (aurostd::FileExist(xvasp.Directory + string("/POSCAR.relax2"))) {
              KBIN::VASP_Recycle(xvasp, "relax2");
            } else {
              if (aurostd::FileExist(xvasp.Directory + string("/POSCAR.relax1"))) {
                KBIN::VASP_Recycle(xvasp, "relax1");
              } else {
                aus << run_type_str << ": RELAX2 or RELAX1 must be present.";
                throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
              }
            }
          } else {
            if (aurostd::FileExist(xvasp.Directory + string("/POSCAR.static"))) {
              KBIN::VASP_Recycle(xvasp, "static");
            } else {
              aus << run_type_str << ": STATIC must be present.";
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
            }
          }
        }

        // STATIC PART ----------------------------------------------------------------------------
        // NOW DO THE STATIC PATCHING POSCAR
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")
            || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_DIELECTRIC")
            || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")
            || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")
            || vflags.KBIN_VASP_RUN.flag("STATIC_DIELECTRIC")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_DIELECTRIC")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")
            || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
          aus << "00000  MESSAGE Patching POSCAR " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          // LOAD THE RELAXED STRUCTURE WHICH WILL BE USED FOR THE DIRECTIONS
          stringstream straus;
          aurostd::file2stringstream(xvasp.Directory + "/POSCAR", straus);
          xvasp.str = xstructure(straus, IOVASP_AUTO);
          xvasp.str.FixLattices();
          rlattice = xvasp.str.lattice; // in rlattice I`ve always the final structure
          const bool STATIC_DEBUG = false; // true;
          if (run_type_str != "STATIC") {
            if (STATIC_DEBUG) {
              aus << "STATIC_DEBUG: " << endl;
              aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry << endl;
              aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
              aus << "STATIC_DEBUG: " << xvasp.str << endl;
              aus << "STATIC_DEBUG: " << rlattice << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            string stringBZ;
            bool foundBZ = false;
            if (vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string == "AUTO") { // CO20210805
              // AS20200724 BEGIN
              //  when KBIN_VASP_KPOINTS_BANDS_LATTICE=AUTO we need to retrieve
              //  the lattice type. For example, if CONVERT_UNIT_CELL=PRES the
              //  call to KPOINTS_Directions() will lead to an error since
              //  vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string is empty
              xvasp.str.GetLatticeType();
              vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();
              vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(xvasp.str.bravais_lattice_variation_type);
              // AS20200724 END
            }

            // always recalculate standardization
            if (vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE")
                || vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT")
                || vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")
                || vflags.KBIN_VASP_KPOINTS_MODE.flag("EXTERNAL")) {
              // nothing
              aus << "00000  MESSAGE PRESERVING ORIGINAL STRUCTURE" << Message(__AFLOW_FILE__, aflags) << endl;
              aus << "00000  MESSAGE ORIGINAL: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
              aus << "00000  MESSAGE ORIGINAL: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            } else {
              aus << "00000  MESSAGE WARNING RECALCULATING STANDARD STRUCTURE" << Message(__AFLOW_FILE__, aflags) << endl;
              aus << "00000  MESSAGE BEFORE: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
              aus << "00000  MESSAGE BEFORE: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // reshuffle the structure
              if (vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
                xvasp.str.Standard_Conventional_UnitCellForm();
              } else {
                xvasp.str.Standard_Primitive_UnitCellForm();
              }
              xvasp.POSCAR.str(std::string());
              xvasp.POSCAR.clear();
              xvasp.POSCAR << xvasp.str;
              xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated", true);
              xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed", true);
              aurostd::stringstream2file(xvasp.POSCAR, string(xvasp.Directory + "/POSCAR"));
              xvasp.str.FixLattices();
              rlattice = xvasp.str.lattice; // in rlattice I`ve always the final structure
              // CO20210805 - both Standard_Conventional_UnitCellForm() and Standard_Primitive_UnitCellForm()
              // return back updated bravais_lattice_variation_type, no need to recalculate with GetLatticeType()
              vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();
              vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(xvasp.str.bravais_lattice_variation_type); // WSETYAWAN mod
              aus << "00000  MESSAGE AFTER: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
              aus << "00000  MESSAGE AFTER: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            stringBZ = LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string, rlattice, vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int, xvasp.str.iomode, foundBZ); // rlattice = updated structure
            if (STATIC_DEBUG) {
              aus << "STATIC_DEBUG: " << endl;
              aus << "STATIC_DEBUG: " << stringBZ << endl;
              aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;
              aus << "STATIC_DEBUG: " << xvasp.str << endl;
              aus << "STATIC_DEBUG: " << rlattice << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            if (foundBZ == false) {
              aus << "Unrecoverable error, lattice not found:" << std::endl;
              aus << xvasp.str;
              throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
            }
          }
          // done with the fixing
          xvasp.str.FixLattices();
          rlattice = xvasp.str.lattice; // in rlattice I`ve always the final structure
          xvasp.aopts.flag("FLAG::POSCAR_PRESERVED", true); // CO20210315 - correct placement // in case of errors it is not lost but recycled
          // NOW DO THE STATIC PATCHING KPOINTS
          aus << "00000  MESSAGE Patching KPOINTS " << Message(__AFLOW_FILE__, aflags) << endl;
          vflags.KBIN_VASP_RUN.push("STATIC"); // use STATIC KPOINTS
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          KBIN::VASP_Produce_KPOINTS(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
          KBIN::VASP_Modify_KPOINTS(xvasp, FileMESSAGE, aflags, vflags);
          aurostd::stringstream2file(xvasp.KPOINTS, string(xvasp.Directory + "/KPOINTS"));
          // NOW DO THE STATIC PATCHING INCAR
          aus << "00000  MESSAGE [" << run_type_str << "] Patching INCAR (static_patching)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          KBIN::VASP_Reread_INCAR(xvasp, FileMESSAGE, aflags); // REREAD IT
          KBIN::XVASP_INCAR_Relax_Static_ON(xvasp, vflags); // FIX
          // do the RWIGS ON
          if (vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) {
            KBIN::XVASP_INCAR_RWIGS_Static(xvasp, vflags, FileMESSAGE, true);
          }
          // done write INCAR
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory + "/INCAR"));
          // NOW DO THE STATIC RUN
          if (vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_DIELECTRIC") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
            xvasp.NRELAXING = xvasp.NRELAX; // 0;
          }
          xvasp.NRELAXING++;
          aus
              << aurostd::PaddedPRE(aurostd::utype2string(11111 * xvasp.NRELAXING), 5, "0")
              << "  STATIC ("
              << run_type_str
              << ") - "
              << xvasp.Directory
              << " - K=["
              << xvasp.str.kpoints_k1
              << " "
              << xvasp.str.kpoints_k2
              << " "
              << xvasp.str.kpoints_k3
              << "]"
              << " - "
              << kflags.KBIN_BIN
              << Message(__AFLOW_FILE__, aflags)
              << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + "]");
            return false;
          }
          if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error (OUTCAR_INCOMPLETE) [" + run_type_str + "]");
            return false;
          }
          KBIN::VASP_Backup(xvasp, true, string("static"));
          xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
          aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size() - 1) << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
        // BANDS PART ----------------------------------------------------------------------------
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
          // NOW DO THE BANDS PATCHING KPOINTS (if necessary...)
          KBIN::VASP_Recycle(xvasp, "static"); // bring back the stuff
          KBIN::VASP_RecycleExtraFile(xvasp, "CHGCAR", "static"); // bring back the stuff
          xvasp.aopts.flag("FLAG::POSCAR_PRESERVED", true); // CO20210315 - correct placement // in case of errors it is not lost but recycled
          xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED", true); // in case of errors it is not lost but recycled
          aus << "00000  MESSAGE Patching KPOINTS with BANDS LATTICE = \"" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << "\"" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          // poscar was already conventionalized in the static part
          xvasp.KPOINTS.clear();
          xvasp.KPOINTS.str(std::string());
          string stringBZ;
          bool foundBZ = false;
          stringBZ = LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string, rlattice, vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int, xvasp.str.iomode, foundBZ); // rlattice = updated structure
          if (foundBZ == false) { // CO20210805
            aus << "Unrecoverable error, lattice not found:" << std::endl;
            aus << xvasp.str;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
          }
          // removed stuff BELOW
          xvasp.KPOINTS << stringBZ;
          aurostd::stringstream2file(xvasp.KPOINTS, string(xvasp.Directory + "/KPOINTS"));
          xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED", true); // don't touch kpoints if there are flaws
          // NOW DO THE BANDS PATCHING INCAR
          aus << "00000  MESSAGE [" << run_type_str << "] Patching INCAR (bands_patching)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          KBIN::VASP_Reread_INCAR(xvasp, FileMESSAGE, aflags); // REREAD IT
          KBIN::XVASP_INCAR_Bands_ON(xvasp, vflags); // FIX
          // do the RWIGS OFF
          if (vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) {
            KBIN::XVASP_INCAR_RWIGS_Static(xvasp, vflags, FileMESSAGE, false);
          }
          // done write INCAR
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory + "/INCAR"));
          // NOW DO THE BANDS RUN
          xvasp.NRELAXING++;
          aus
              << 11111 * xvasp.NRELAXING
              << "  BANDS ("
              << run_type_str
              << ") - "
              << xvasp.Directory
              << " - K=["
              << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string
              << ","
              << vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int
              << "]"
              << Message(__AFLOW_FILE__, aflags)
              << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " BANDS]");
            return false;
          }
          if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error (OUTCAR_INCOMPLETE) [" + run_type_str + " BANDS]");
            return false;
          }
          KBIN::VASP_Backup(xvasp, true, string("bands"));
          xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED", false); // bands are done... I can refix the KPOINTS
          xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
          aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size() - 1) << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }

        // DIELECTRIC PART ----------------------------------------------------------------------------
        if (vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_DIELECTRIC") || vflags.KBIN_VASP_RUN.flag("STATIC_DIELECTRIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_DIELECTRIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_DIELECTRIC")) {
          // Keep static files and old CHGCAR
          KBIN::VASP_Recycle(xvasp, "static");
          KBIN::VASP_RecycleExtraFile(xvasp, "CHGCAR", "static");
          KBIN::VASP_RecycleExtraFile(xvasp, "POSCAR", "static");
          xvasp.aopts.flag("FLAG::POSCAR_PRESERVED", true);
          xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED", true);
          // NOW DO THE DIELECTRIC PATCHING KPOINTS
          aus << "00000  MESSAGE Patching KPOINTS " << Message(__AFLOW_FILE__, aflags) << endl;
          vflags.KBIN_VASP_RUN.push("DIELECTRIC"); // use DIELECTRIC KPOINTS
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          KBIN::VASP_Produce_KPOINTS(xvasp, AflowIn, FileMESSAGE, aflags, kflags, vflags);
          KBIN::VASP_Modify_KPOINTS(xvasp, FileMESSAGE, aflags, vflags);
          aurostd::stringstream2file(xvasp.KPOINTS, string(xvasp.Directory + "/KPOINTS"));
          // NOW DO THE DIELECTRIC PATCHING INCAR
          aus << "00000  MESSAGE [" << run_type_str << "] Patching INCAR (dielectric_patching)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          KBIN::VASP_Reread_INCAR(xvasp, FileMESSAGE, aflags); // REREAD IT
          KBIN::XVASP_INCAR_Dielectric_ON(xvasp, vflags, aflags); // FIX
          if (vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) {
            KBIN::XVASP_INCAR_RWIGS_Static(xvasp, vflags, FileMESSAGE, false);
          }
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory + "/INCAR"));
          // NOW DO THE DIELECTRIC RUN
          xvasp.NRELAXING++;
          aus
              << 11111 * xvasp.NRELAXING
              << "  DIELECTRIC ("
              << run_type_str
              << ") - "
              << xvasp.Directory
              << " - K=["
              << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string
              << ","
              << vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int
              << "]"
              << Message(__AFLOW_FILE__, aflags)
              << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " DIELECTRIC]");
            return false;
          }
          if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error (OUTCAR_INCOMPLETE) [" + run_type_str + " DIELECTRIC]");
            return false;
          }
          KBIN::VASP_Backup(xvasp, true, string("dielectric"));
          xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
          aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size() - 1) << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }

        // FINISHED
        xvasp.NRELAXING++;
        aus << 11111 * xvasp.NRELAXING << "  END (" << run_type_str << ")        - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        vflags.KBIN_VASP_RUN.flag("STATIC", false); // put back the options

        // CLEAN-UP
        if (!vflags.KBIN_VASP_RUN.flag("RELAX")) {
          aurostd::RemoveFile(xvasp.Directory, std::regex("CHG.relax.*"));
          aurostd::RemoveFile(xvasp.Directory, std::regex("CHGCAR.relax.*"));
        }
        aurostd::RemoveFile(xvasp.Directory, std::regex("CHG.bands.*"));
        aurostd::RemoveFile(xvasp.Directory, std::regex("CHGCAR.bands.*"));
        aurostd::RemoveFile(xvasp.Directory, std::regex("CHG.dielectric.*"));
        aurostd::RemoveFile(xvasp.Directory, std::regex("CHGCAR.dielectric.*"));

        // --------------------------------------------------------------------------------------------------------------------
        // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL
        // PRL 105, 196403 (2010)
        if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) {
          bool delsol_d = false;
          float NELECT = 0;
          float Nr = 0;
          const xmatrix<double> rlattice(3, 3);
          string stmp;
          const string fnamedelsol = xvasp.Directory + string("/delsol.tmp");
          stringstream command;
          stringstream strdelsol;
          const ifstream fdelsol;
          if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) {
            run_type_str = "REPEAT_DELSOL";
          }
          aus << "00000  MESSAGE MODE= (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          // LOAD FORMER LOCK
          if (aurostd::FileExist(xvasp.Directory + string("/REPEAT_DELSOL"))) {
            stringstream lock_recycled;
            aurostd::file2stringstream(xvasp.Directory + "/REPEAT_DELSOL", lock_recycled);
            aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
            aus << "XXXXX FORMER LOCK BEGIN, recycled (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
            aus << lock_recycled.str();
            aus << "XXXXX FORMER LOCK END, recycled (" << run_type_str << ") - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
            aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
          // UNZIP EVERYTHING
          KBIN::DecompressDirectory(xvasp.Directory);

          // copy INCAR, POSCAR, KPOINTS, POTCAR from *.static
          KBIN::VASP_Recycle(xvasp, "static"); // bring back the stuff
          // Scanning whether it is sp or spd from the POTCAR
          command.clear();
          command.str(std::string());
          command << "cd " << xvasp.Directory << endl;
          command << R"(grep VRHFIN POTCAR.static | sed 's/:/\n/g' | grep -v VRHFIN > delsol.tmp)" << endl;
          aurostd::execute(command);
          strdelsol.clear();
          strdelsol.str(std::string());
          aurostd::file2stringstream(xvasp.Directory + "/delsol.tmp", strdelsol);
          aurostd::RemoveFile(xvasp.Directory + "/delsol.tmp");
          delsol_d = false;
          if ((aurostd::substring2bool(strdelsol.str(), "d"))) {
            delsol_d = true;
          }
          // Scanning NELECT from OUTCAR.static
          command.clear();
          command.str(std::string());
          command << "cd " << xvasp.Directory << endl;
          command << R"(grep NELECT OUTCAR.static | sed 's/=/\n/g' | grep -v NELECT > delsol.tmp)" << endl;
          aurostd::execute(command);
          strdelsol.clear();
          strdelsol.str(std::string());
          aurostd::file2stringstream(xvasp.Directory + "/delsol.tmp", strdelsol);
          aurostd::RemoveFile(xvasp.Directory + "/delsol.tmp");
          strdelsol >> NELECT;
          // if(NELECT<1.0) ;//need to add error handling here
          Nr = NELECT / 68.0;
          if (delsol_d) {
            Nr = NELECT / 72.0;
          }
          aus << "DELSOL: NELECT N0=" << NELECT << endl << "DELSOL: NELECT n=" << Nr << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

          // NOW MODIFY THE INCAR
          aus << "00000  MESSAGE [" << run_type_str << "] modifying INCAR (delsol_patching)" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          KBIN::VASP_Reread_INCAR(xvasp, FileMESSAGE, aflags); // REREAD IT
          stmp = "NELECT=" + aurostd::utype2string(NELECT + Nr);
          xvasp.INCAR << aurostd::PaddedPOST(stmp, _incarpad_) << " # NELECT = N0 + n for DELSOL plus" << endl;
          // do the RWIGS OFF
          if (vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) {
            KBIN::XVASP_INCAR_RWIGS_Static(xvasp, vflags, FileMESSAGE, false);
          }
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory + "/INCAR"));

          KBIN::VASP_RecycleExtraFile(xvasp, "POSCAR", "static"); // bring back the stuff
          KBIN::VASP_RecycleExtraFile(xvasp, "KPOINTS", "static"); // bring back the stuff

          // NOW RUN DELSOL plus
          uint vrelax = 7;
          aus << 11111 * vrelax << "  DELSOL plus (" << run_type_str << ") - " << xvasp.Directory << " - " << stmp << " " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          vrelax++;
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + "]");
            return false;
          }
          if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error (OUTCAR_INCOMPLETE) [" + run_type_str + "]");
            return false;
          }
          KBIN::VASP_Backup(xvasp, true, string("dsolp"));

          // NOW DO DELSOL minus
          KBIN::VASP_Recycle(xvasp, "static");
          KBIN::VASP_Reread_INCAR(xvasp, FileMESSAGE, aflags);
          stmp = "NELECT=" + aurostd::utype2string(NELECT - Nr);
          xvasp.INCAR << aurostd::PaddedPOST(stmp, _incarpad_) << " # NELECT = N0 - n for DELSOL minus" << endl;
          // do the RWIGS OFF
          if (vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) {
            KBIN::XVASP_INCAR_RWIGS_Static(xvasp, vflags, FileMESSAGE, false);
          }
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory + "/INCAR"));

          KBIN::VASP_RecycleExtraFile(xvasp, "POSCAR", "static"); // bring back the stuff
          KBIN::VASP_RecycleExtraFile(xvasp, "KPOINTS", "static"); // bring back the stuff
          // NOW RUN DELSOL minus
          aus << 11111 * vrelax << "  DELSOL minus (" << run_type_str << ") - " << xvasp.Directory << " - " << stmp << " " << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          vrelax++;
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " minus]");
            return false;
          }
          if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error (OUTCAR_INCOMPLETE) [" + run_type_str + " minus]");
            return false;
          }
          KBIN::VASP_Backup(xvasp, false, string("dsolm"));
          // FINISHED
          aus << 11111 * vrelax << "  END (" << run_type_str << ")        - " << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          vflags.KBIN_VASP_RUN.flag("STATIC", false); // put back the options
          // CLEAN-UP BY WSETYAWAN
          if (vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
            aurostd::RemoveFile(xvasp.Directory, std::regex("CHG.dsol.*"));
            aurostd::RemoveFile(xvasp.Directory, std::regex("CHGCAR.dsol.*"));
          }
        }
        // --------------------------------------------------------------------------------------------------------------------
        // KPOINTS KPOINTS KPOINTS
        if (vflags.KBIN_VASP_RUN.flag("KPOINTS")) { // NON THREADED
          KBIN::VASP_Write_INPUT(xvasp, vflags); // VASP VASP WRITE
          xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed", true); // BACKUP KPOINTS
          aurostd::stringstream2file(xvasp.KPOINTS_orig, string(xvasp.Directory + "/KPOINTS.orig")); // BACKUP KPOINTS
          const int kbak_k1 = xvasp.str.kpoints_k1;
          const int kbak_k2 = xvasp.str.kpoints_k2;
          const int kbak_k3 = xvasp.str.kpoints_k3;
          int kk1;
          int kk2;
          int kk3;
          string relax;
          string relaxfile;
          // 1st step: 1/2
          kk1 = (kbak_k1 + 1) / 2;
          kk2 = (kbak_k2 + 1) / 2;
          kk3 = (kbak_k3 + 1) / 2;
          relax = "11111a ";
          relaxfile = "relax1";
          aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          xvasp.str.kpoints_k1 = aurostd::min(kk1, kbak_k1);
          xvasp.str.kpoints_k2 = aurostd::min(kk2, kbak_k2);
          xvasp.str.kpoints_k3 = aurostd::min(kk3, kbak_k3);
          xvasp.str.kpoints_kmax = aurostd::max(xvasp.str.kpoints_k1, xvasp.str.kpoints_k2, xvasp.str.kpoints_k3);
          KBIN::XVASP_KPOINTS_numbers2string(xvasp);
          aurostd::stringstream2file(xvasp.KPOINTS, string(xvasp.Directory + "/KPOINTS")); // BACKUP KPOINTS
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, relaxfile, relaxfile, false, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " 1]");
            return false;
          }
          // 2nd step: 1/2
          kk1 = 3 * (kbak_k1 + 1) / 4;
          kk2 = 3 * (kbak_k2 + 1) / 4;
          kk3 = 3 * (kbak_k3 + 1) / 4;
          relax = "11111b ";
          relaxfile = "relax1";
          aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          xvasp.str.kpoints_k1 = aurostd::min(kk1, kbak_k1);
          xvasp.str.kpoints_k2 = aurostd::min(kk2, kbak_k2);
          xvasp.str.kpoints_k3 = aurostd::min(kk3, kbak_k3);
          xvasp.str.kpoints_kmax = aurostd::max(xvasp.str.kpoints_k1, xvasp.str.kpoints_k2, xvasp.str.kpoints_k3);
          KBIN::XVASP_KPOINTS_numbers2string(xvasp);
          aurostd::stringstream2file(xvasp.KPOINTS, string(xvasp.Directory + "/KPOINTS")); // BACKUP KPOINTS
          if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, relaxfile, relaxfile, false, FileMESSAGE)) {
            KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " 2]");
            return false;
          }
          // 3rd step: 1/2
          kk1 = kbak_k1;
          kk2 = kbak_k2;
          kk3 = kbak_k3;
          relax = "11111c ";
          relaxfile = "relax1";
          aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          xvasp.str.kpoints_k1 = aurostd::min(kk1, kbak_k1);
          xvasp.str.kpoints_k2 = aurostd::min(kk2, kbak_k2);
          xvasp.str.kpoints_k3 = aurostd::min(kk3, kbak_k3);
          xvasp.str.kpoints_kmax = aurostd::max(xvasp.str.kpoints_k1, xvasp.str.kpoints_k2, xvasp.str.kpoints_k3);
          KBIN::XVASP_KPOINTS_numbers2string(xvasp);
          aurostd::stringstream2file(xvasp.KPOINTS, string(xvasp.Directory + "/KPOINTS")); // BACKUP KPOINTS
          // ----------------------------------------------------------
          // with PAWGGA2
          if (PAWGGA2) {
            // STEP 1
            aus << "11111  RELAXATION - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax2paw_gga", false, FileMESSAGE)) {
              KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " PAWGGA2]");
              return false;
            }
            aus << "22222  END        - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
          // ----------------------------------------------------------
          // norma, without PAWGGA2
          if (!PAWGGA2) {
            int vrelax = 1;
            for (int i = 1; i <= xvasp.NRELAX; i++) {
              aus << 11111 * vrelax << "   RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              if (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int > 0) {
                KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp, i); // ADIABATIC
              }
              if (i < xvasp.NRELAX && !KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string(vrelax), "relax" + aurostd::utype2string(vrelax), true, FileMESSAGE)) {
                KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " 4]");
                return false;
              }
              if (i == xvasp.NRELAX && !KBIN::VASP_Run(xvasp, aflags, kflags, vflags, "relax" + aurostd::utype2string(vrelax), true, FileMESSAGE)) {
                KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  runtime error [" + run_type_str + " 5]");
                return false;
              }
              vrelax++;
            }
            aus << 11111 * vrelax << "   END        - " << xvasp.Directory << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
        } // KPOINTS KPOINTS KPOINTS
        // ***************************************************************************
        // POSTSCRIPT
        if (!vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
          if (kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) {
            KBIN::RUN_DirectoryScript(aflags, DEFAULT_AFLOW_POSTSCRIPT_COMMAND, DEFAULT_AFLOW_POSTSCRIPT_OUT);
          }
        }
        // ***************************************************************************
      }
      // **********
      // some verbose
      if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus
            << "00000  MESSAGE END loop "
            << xvasp.POSCAR_index + 1
            << "/"
            << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() // CO20200624 - +1
            << Message(__AFLOW_FILE__, aflags)
            << endl;
        aus << "00000  MESSAGE END loop in directory =" << xvasp.Directory << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        // compress the subdirectories
        if (kflags.KZIP_COMPRESS && !KBIN::CompressDirectory(aflags)) {
          return false;
        }
      }
      aflags = aflags_backup;
      kflags = kflags_backup; // RESTORE
    } // LOOP ixvasp
    // ***************************************************************************
    aflags = aflags_backup;
    kflags = kflags_backup; // RESTORE
    // POSTSCRIPT
    if (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
      if (kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) {
        KBIN::RUN_DirectoryScript(aflags, DEFAULT_AFLOW_POSTSCRIPT_COMMAND, DEFAULT_AFLOW_POSTSCRIPT_OUT);
      }
    }
    // ***************************************************************************
    FileAFLOWIN.clear();
    FileAFLOWIN.close();
    return true;
  }
} // namespace KBIN

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

int CheckStringInFile(string FileIn, string str, int PID, int TID) { // CO20200502 - threadID
  int out;
  ostringstream aus_exec;
  aus_exec << "cat " << FileIn << " | grep -c \"" << str << "\" > aflow.out." << PID << "." << TID << endl; // CO20200502 - threadID
  aurostd::execute(aus_exec);
  ifstream FileCHECK;
  stringstream FineNameAflowTmpPIDTID; // CO20200502 - threadID
  FineNameAflowTmpPIDTID << "aflow.out." << PID << "." << TID; // CO20200502 - threadID
  FileCHECK.open(FineNameAflowTmpPIDTID.str().c_str(), std::ios::in); // CO20200502 - threadID
  FileCHECK >> out;
  FileCHECK.clear();
  FileCHECK.close();
  aurostd::execute(aus_exec);
  aurostd::RemoveFile("aflow.out." + aurostd::utype2string(PID) + "." + aurostd::utype2string(TID));
  return out;
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// FLAG Class for KBIN_VASP_Run

namespace KBIN {
  bool ReachedAccuracy2bool(const string& scheme, const aurostd::xoption& xRequiresAccuracy, const aurostd::xoption& xmessage, bool vasp_still_running) { // CO20210315
    const bool LDEBUG = (false || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " xRequiresAccuracy.flag(\"" << scheme << "\")=" << xRequiresAccuracy.flag(scheme) << endl;
      cerr << __AFLOW_FUNC__ << " vasp_still_running=" << vasp_still_running << endl;
      cerr << __AFLOW_FUNC__ << " xmessage.flag(\"REACHED_ACCURACY\")=" << xmessage.flag("REACHED_ACCURACY") << endl;
    }

    if (xRequiresAccuracy.flag(scheme) == false) {
      return true;
    } // this scheme does not require xmessage.flag("REACHED_ACCURACY"), return true for '&&' logic

    // CO20210601 - I am decoupling vasp_still_running and xmessage.flag("REACHED_ACCURACY") FOR RELAXATIONS ONLY
    // on qrats, I noticed that a calculation can reach accuracy, but not finish (incomplete OUTCAR), thus it is frozen
    // a good solution would be to take the CONTCAR as the new input and restart (RELAXATIONS ONLY)
    // to diagnose and treat this problem correctly, we need to consider vasp_still_running and xmessage.flag("REACHED_ACCURACY") independently (for this case only)
    //[CO20210602 - must use, otherwise NUM_PROB appears early]if(vasp_still_running==false){;}  //keep busy
    return (vasp_still_running == false && xmessage.flag("REACHED_ACCURACY") == false); // vasp_still_running==false && - xmessage.flag("REACHED_ACCURACY") should already understand whether to consider vasp_still_running
  }
  void VASP_InitWarnings(std::unordered_map<string, bool>& vasp_warnings) {
    vasp_warnings = {
        {                                             "reached required accuracy", false},
        {                    "VERY BAD NEWS! internal error in subroutine SGRCON", false},
        {                     "Found some non-integer element in rotation matrix", false},
        {                    "VERY BAD NEWS! internal error in subroutine IBZKPT", false},
        {"Reciprocal lattice and k-lattice belong to different class of lattices", false},
        {                                                             "NKX>IKPTD", false},
        {                                                             "NKY>IKPTD", false},
        {                                                             "NKZ>IKPTD", false},
        {                    "VERY BAD NEWS! internal error in subroutine INVGRP", false},
        {           "inverse of rotation matrix was not found (increase SYMPREC)", false},
        {                               "WARNING in EDDRMM: call to ZHEGV failed", false},
        {                                                              "num prob", false},
        {                                          "ZBRENT: can't locate minimum", false},
        {                                     "ZBRENT: fatal error in bracketing", false},
        {                                            "we suggest increasing NELM", false},
        {                                          "BRMIX: very serious problems", false},
        {                     "WARNING: Sub-Space-Matrix is not hermitian in DAV", false},
        {                      "WARNING: DENTET: can't reach specified precision", false},
        {                                    "Error EDDDAV: Call to ZHEGV failed", false},
        {                                              "EFIELD_PEAD is too large", false},
        {                                 "EFIELD_PEAD are too large for comfort", false},
        {                      "ERROR FEXCF: supplied exchange-correlation table", false},
        {                      "ERROR FEXCP: supplied Exchange-correletion table", false},
        {                                              "shift your grid to Gamma", false},
        {                             "LRF_COMMUTATOR internal error: the vector", false},
        {                                                        AFLOW_MEMORY_TAG, false},
        {                  "BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES", false},
        {                                                         "EXIT CODE: 11", false},
        {                                                  "KILLED BY SIGNAL: 11", false},
        {                                                        "EXIT CODE: 139", false},
        {                                                 "KILLED BY SIGNAL: 139", false},
        {                                                        "EXIT CODE: 174", false},
        {                                                 "KILLED BY SIGNAL: 174", false},
        {                              "distance between some ions is very small", false},
        {                                                                "NBANDS", false},
        {               "number of bands is not sufficient to hold all electrons", false},
        {                    "Number of bands NBANDS too small to hold electrons", false},
        {                        "Your highest band is occupied at some k-points", false},
        {             "number of bands has been changed from the values supplied", false},
        {                                                        "now  NBANDS  =", false},
        {                                               "please rerun with NPAR=", false},
        {                                                              "NPAR = 4", false},
        {                                                  "NPAR=number of cores", false},
        {                                                  "NPAR=number of nodes", false},
        {                        "Please remove the tag NPAR from the INCAR file", false},
        {                     "WARNING: PSMAXN for non-local potential too small", false},
        {                                              "REAL_OPT: internal ERROR", false},
        {                                           "ERROR in RE_READ_KPOINTS_RD", false},
        {                                                   "switch off symmetry", false},
        {                                              "REAL_OPT: internal ERROR", false},
        {                                       "REAL_OPTLAY: internal error (1)", false},
        {                                         "LAPACK: Routine ZPOTRF failed", false},
    };
  }
  void VASP_ProcessWarnings(_xvasp& xvasp, _aflags& aflags, _kflags& kflags, size_t& offset, std::unordered_map<string, bool>& vasp_warnings, aurostd::xoption& xmessage, aurostd::xoption& xwarning, ofstream& FileMESSAGE) { // CO20210315
    aurostd::xoption xmonitor;
    return VASP_ProcessWarnings(xvasp, aflags, kflags, offset, vasp_warnings, xmessage, xwarning, xmonitor, FileMESSAGE);
  }
  void VASP_ProcessWarnings(_xvasp& xvasp, _aflags& aflags, _kflags& kflags, size_t& offset, std::unordered_map<string, bool>& vasp_warnings, aurostd::xoption& xmessage, aurostd::xoption& xwarning, aurostd::xoption& xmonitor, ofstream& FileMESSAGE) { // CO20210315
    const bool LDEBUG = (false || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    stringstream aus;
    const bool VERBOSE = (false || XHOST.vflag_control.flag("MONITOR_VASP") == false || LDEBUG);

    if (!aurostd::FileExist(xvasp.Directory + "/" + DEFAULT_VASP_OUT)) {
      return;
    }
    if (!aurostd::FileExist(xvasp.Directory + "/INCAR")) {
      return;
    }
    if (!aurostd::FileExist(aflags.Directory + "/" + _AFLOWLOCK_)) {
      return;
    } // we needed it above to get the vasp_bin
    const bool vasp_monitor_running = AFLOW_MONITOR_instance_running(aflags);

    const long int tmod_outcar = aurostd::SecondsSinceFileModified(xvasp.Directory + "/" + "OUTCAR"); // better to look at OUTCAR than vasp.out, when vasp is killed you get errors in vasp.out, resetting the time
    const unsigned long long int fsize_vaspout = aurostd::FileSize(xvasp.Directory + "/" + DEFAULT_VASP_OUT);
    if (VERBOSE) {
      aus << "00000  MESSAGE time since " << "OUTCAR" << " last modified: " << tmod_outcar << " seconds (max=" << SECONDS_STALE_OUTCAR << " seconds)" << Message(__AFLOW_FILE__, aflags) << endl;
      if (LDEBUG) {
        cerr << aus.str();
      }
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      aus << "00000  MESSAGE size of " << DEFAULT_VASP_OUT << ": " << fsize_vaspout << " bytes (max=" << BYTES_MAX_VASP_OUT << " bytes)" << Message(__AFLOW_FILE__, aflags) << endl;
      if (LDEBUG) {
        cerr << aus.str();
      }
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    // do memory check
    double usage_percentage_ram = 0.0;
    double usage_percentage_swap = 0.0;
    bool approaching_oom = false;
    if (aurostd::GetMemoryUsagePercentage(usage_percentage_ram, usage_percentage_swap)) {
      approaching_oom = (usage_percentage_ram >= MEMORY_MAX_USAGE_RAM && usage_percentage_swap >= MEMORY_MAX_USAGE_SWAP);
    }
    if (approaching_oom) { // might be a quick memory spike, try again
      if (VERBOSE) {
        aus << "00000  MESSAGE ram memory used: " << aurostd::utype2string(usage_percentage_ram, 2, FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_RAM << "%)" << Message(__AFLOW_FILE__, aflags) << endl;
        aus << "00000  MESSAGE swap memory used: " << aurostd::utype2string(usage_percentage_swap, 2, FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_SWAP << "%)" << Message(__AFLOW_FILE__, aflags) << endl;
        aus << "00000  MESSAGE reading memory again after " << SECONDS_SLEEP_VASP_MONITOR << " second sleep" << endl;
        if (LDEBUG) {
          cerr << aus.str();
        }
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      aurostd::Sleep(SECONDS_SLEEP_VASP_MONITOR);
      if (aurostd::GetMemoryUsagePercentage(usage_percentage_ram, usage_percentage_swap)) {
        approaching_oom = (usage_percentage_ram >= MEMORY_MAX_USAGE_RAM && usage_percentage_swap >= MEMORY_MAX_USAGE_SWAP);
      }
    }
    if (VERBOSE || approaching_oom) {
      aus << "00000  MESSAGE ram memory used: " << aurostd::utype2string(usage_percentage_ram, 2, FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_RAM << "%)" << Message(__AFLOW_FILE__, aflags) << endl;
      aus << "00000  MESSAGE swap memory used: " << aurostd::utype2string(usage_percentage_swap, 2, FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_SWAP << "%)" << Message(__AFLOW_FILE__, aflags) << endl;
      if (LDEBUG) {
        cerr << aus.str();
      }
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    // get INCAR
    VASP_Reread_INCAR(xvasp); // preload incar

    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " [1]" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    // get whether relaxing or not
    bool relaxing = false;
    if (aurostd::kvpair2bool(xvasp.INCAR, "IBRION", "=")) {
      if (aurostd::kvpair2bool(xvasp.INCAR, "NSW", "=") && aurostd::kvpair2utype<int>(xvasp.INCAR, "NSW", "=") > 0) {
        relaxing = true;
      }
    }

    // need to get ISYM and ISPIN (ISPIND too)
    // get ISYM
    int isym_current = 2; // CO20200624 - VASP default for non-USPP runs, add check later for USPP //https://www.vasp.at/wiki/index.php/ISYM
    if (aurostd::substring2bool(xvasp.INCAR, "ISYM=")) {
      isym_current = aurostd::substring2utype<int>(xvasp.INCAR, "ISYM=");
    }
    // get ISPIN
    int ispin_current = 1; // CO20200624 - VASP default  //https://www.vasp.at/wiki/index.php/ISPIN
    if (aurostd::substring2bool(xvasp.INCAR, "ISPIN=")) {
      ispin_current = aurostd::substring2utype<int>(xvasp.INCAR, "ISPIN=");
    }

    // there are issues getting the correct vasp binary since this is an entirely different aflow instance
    // we might run the other aflow instance with --mpi or --machine flags that affect which vasp binary we use
    // therefore, the most robust way to define the binary is to search the LOCK file
    string vasp_bin;
    string vasp_pgid;
    GetVASPBinaryFromLOCK(xvasp.Directory, vasp_bin);
    vasp_bin = aurostd::basename(vasp_bin); // remove directory stuff
    vasp_pgid = aurostd::utype2string(getpgrp()); // SD20220406 - need PGID for VASP_instance_running
    if (vasp_bin.empty()) { // rely on defaults here in case we're not running --monitor_vasp
      vasp_bin = kflags.KBIN_MPI_BIN;
      if (!(kflags.KBIN_MPI == true || XHOST.MPI == true)) {
        vasp_bin = kflags.KBIN_BIN;
      }
    }

    // CO20210315 - determine if vasp is still running
    const bool vasp_still_running = VASP_instance_running(vasp_bin, vasp_pgid); // SD20220406 - we are not using the --kill_vasp_all flag anymore
    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " vasp_still_running=" << vasp_still_running << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    uint i = 0;

    // determine which schemes require xmessage.flag("REACHED_ACCURACY")
    // CO+AS202010315 - considering DENTET for this list, patches seem to be working. See: LIB3/CTa_pvTi_sv/ABC2_tP8_123_h_h_abg-001.ABC
    aurostd::xoption xRequiresAccuracy;
    const vector<string> vtokens{"DAV", "EDDRMM", "IBZKPT", "NUM_PROB", "ZBRENT"}; // DENTET
    for (i = 0; i < vtokens.size(); i++) {
      xRequiresAccuracy.flag(vtokens[i], true);
    }

    xwarning.clear(); // CO20210315 - very important to clear!

    const bool vasp_oszicar_unconverged = KBIN::VASP_OSZICARUnconverged(xvasp.Directory + "/OSZICAR", xvasp.Directory + "/OUTCAR");

    // WARNINGS START
    aurostd::substrings_map_present_file(xvasp.Directory + "/" + DEFAULT_VASP_OUT, vasp_warnings, offset);
    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " checking warnings" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    xwarning.flag("OUTCAR_INCOMPLETE", vasp_still_running == false && !KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, false)); // CO20201111
    // CO20210601 - I am decoupling vasp_still_running and xmessage.flag("REACHED_ACCURACY") FOR RELAXATIONS ONLY
    // on qrats, I noticed that a calculation can reach accuracy, but not finish (incomplete OUTCAR), thus it is frozen
    // a good solution would be to take the CONTCAR as the new input and restart (RELAXATIONS ONLY)
    // to diagnose and treat this problem correctly, we need to consider vasp_still_running and xmessage.flag("REACHED_ACCURACY") independently (for this case only)
    const bool reached_accuracy_relaxing = (relaxing == true && vasp_warnings["reached required accuracy"]); //[CO20210601 - only check for static, not relax (might be frozen, and we can fix by appending CONTCAR)]vasp_still_running==false &&
    const bool reached_accuracy_static = (relaxing == false && vasp_still_running == false && vasp_oszicar_unconverged == false && xwarning.flag("OUTCAR_INCOMPLETE") == false); // CO20210315 - "reached accuracy" for static/bands calculation is a converged electronic scf, need to also check OUTCAR_INCOMPLETE, as a converged OSZICAR might actually be a run that ended because of an error
    xmessage.flag("REACHED_ACCURACY", (reached_accuracy_relaxing || reached_accuracy_static));

    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " relaxing=" << relaxing << endl;
      aus << __AFLOW_FUNC__ << " vasp_oszicar_unconverged=" << vasp_oszicar_unconverged << endl;
      aus << __AFLOW_FUNC__ << " xwarning.flag(\"OUTCAR_INCOMPLETE\")=" << xwarning.flag("OUTCAR_INCOMPLETE") << endl;
      aus << __AFLOW_FUNC__ << " xmessage.flag(\"REACHED_ACCURACY\")=" << xmessage.flag("REACHED_ACCURACY") << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    string scheme;
    bool found_warning = false;
    // VASP's internal symmetry routines START
    // CO20200624 - these are all related to VASP's internal symmetry routines
    // they would all benefit from similar fixes (except NKXYZ_IKPTD which requires KPOINTS to be reduced)
    // SGRCON+NIRMAT: https://dannyvanpoucke.be/vasp-errors-en/
    // IBZKPT+KKSYM: http://www.error.wiki/VERY_BAD_NEWS!_internal_error_in_subroutine_IBZKPT
    // IBZKPT+KKSYM: https://www.vasp.at/forum/viewtopic.php?f=3&t=7811
    // INVGRP+SYMPREC: https://www.vasp.at/forum/viewtopic.php?t=486
    // INVGRP+SYMPREC: https://www.vasp.at/forum/viewtopic.php?f=3&t=9435&p=9473
    //
    scheme = "SGRCON"; // usually goes with NIRMAT
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["VERY BAD NEWS! internal error in subroutine SGRCON"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NIRMAT"; // usually goes with SGRCON
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["Found some non-integer element in rotation matrix"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "IBZKPT"; // usually goes with KKSYM or NKXYZ_IKPTD
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["VERY BAD NEWS! internal error in subroutine IBZKPT"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "KKSYM"; // usually goes with IBZKPT
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["Reciprocal lattice and k-lattice belong to different class of lattices"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NKXYZ_IKPTD"; // usually goes with IBZKPT
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    bool found_nkxyz_ikptd = false; // break up for readability
    found_nkxyz_ikptd = (found_nkxyz_ikptd || vasp_warnings["NKX>IKPTD"]);
    found_nkxyz_ikptd = (found_nkxyz_ikptd || vasp_warnings["NKY>IKPTD"]);
    found_nkxyz_ikptd = (found_nkxyz_ikptd || vasp_warnings["NKZ>IKPTD"]);
    found_warning = (found_warning && found_nkxyz_ikptd);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "INVGRP"; // usually goes with SYMPREC
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["VERY BAD NEWS! internal error in subroutine INVGRP"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "SYMPREC"; // usually goes with INVGR
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["inverse of rotation matrix was not found (increase SYMPREC)"]);
    xwarning.flag(scheme, found_warning);
    // VASP's internal symmetry routines END

    // VASP issues with RMM-DIIS START
    scheme = "EDDRMM"; // CO20210315 - look here https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["WARNING in EDDRMM: call to ZHEGV failed"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NUM_PROB"; // CO20210315 - look here https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS  //CO20210315 - num prob can be a big problem for the DOS: https://www.vasp.at/forum/viewtopic.php?f=3&t=18028
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["num prob"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "ZBRENT"; // CO20210315 - look here https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings["ZBRENT: can't locate minimum"] || vasp_warnings["ZBRENT: fatal error in bracketing"]));
    xwarning.flag(scheme, found_warning);
    // VASP issues with RMM-DIIS END

    // CSLOSHING and NELM START
    // CSLOSHING and NELM warnings are similar, CSLOSHING will apply a fix before VASP finishes running, while NELM only cares about the LAST iteration
    // do not check turn on CSLOSHING with NELM, caused collisions with EDDDAV. it bounces back and forth between ALGO=VERFAST and ALGO=NORMAL
    // however, if CSLOSHING, turn on NELM, as the NELM patches will work for CSLOSHING too
    // the check for xwarning.flag("OUTCAR_INCOMPLETE")==false ensures we don't flag a run that was killed by --monitor_vasp, not xmessage.flag("REACHED_ACCURACY") (must be unconverged)
    scheme = "NELM"; // CO20210315
    found_warning = (vasp_still_running == false && xwarning.flag("OUTCAR_INCOMPLETE") == false && (vasp_oszicar_unconverged || vasp_warnings["we suggest increasing NELM"])); // check from OSZICAR
    xwarning.flag(scheme, found_warning);
    //
    scheme = "CSLOSHING";
    found_warning = (KBIN::VASP_OSZICARUnconverging(xvasp.Directory)); // check from OSZICAR //xwarning.flag("NELM")
    xwarning.flag(scheme, found_warning);
    if (xwarning.flag("CSLOSHING")) {
      xwarning.flag("NELM", true);
    }
    // CSLOSHING and NELM END

    // ALL OTHERS (in alphabetic order) START
    scheme = "BRMIX";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["BRMIX: very serious problems"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "CALC_FROZEN"; // CO20210315
    found_warning = (tmod_outcar >= SECONDS_STALE_OUTCAR);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "DAV";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["WARNING: Sub-Space-Matrix is not hermitian in DAV"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "DENTET";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["WARNING: DENTET: can't reach specified precision"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "EDDDAV";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["Error EDDDAV: Call to ZHEGV failed"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "EFIELD_PEAD";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings["EFIELD_PEAD is too large"] || vasp_warnings["EFIELD_PEAD are too large for comfort"])); // 20190704 - new VASP
    xwarning.flag(scheme, found_warning);
    //
    scheme = "EXCCOR";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings["ERROR FEXCF: supplied exchange-correlation table"] || vasp_warnings["ERROR FEXCP: supplied Exchange-correletion table"])); // CO20210315 - some formatting changes between versions (also some bad spelling)
    xwarning.flag(scheme, found_warning);
    //
    scheme = "GAMMA_SHIFT";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["shift your grid to Gamma"]); // CO20190704 - captures both old/new versions of VASP
    xwarning.flag(scheme, found_warning);
    //
    scheme = "LRF_COMMUTATOR"; // GET ALL TIMES
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["LRF_COMMUTATOR internal error: the vector"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "MEMORY";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings[AFLOW_MEMORY_TAG] || (XHOST.vflag_control.flag("KILL_VASP_OOM") && approaching_oom)));
    xwarning.flag(scheme, found_warning);
    //
    // on qrats, we see this
    //===================================================================================
    //=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
    //=   PID 27264 RUNNING AT prod-0004
    //=   EXIT CODE: 9
    //=   CLEANING UP REMAINING PROCESSES
    //=   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
    //===================================================================================
    // YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Killed (signal 9)
    // This typically refers to a problem with your application.
    // Please see the FAQ page for debugging suggestions
    //
    // on quser, we see nothing...
    //
    // on x, we see this
    //===================================================================================
    //=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
    //=   RANK 58 PID 63228 RUNNING AT node6
    //=   KILLED BY SIGNAL: 9 (Killed)
    //===================================================================================
    //
    const bool found_bad_termination = vasp_warnings["BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES"];
    bool found_exit_code = false; // specific exit code
    //
    scheme = "MPICH0"; // 0 just means that it is generic, maybe fix name later
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && found_bad_termination);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "MPICH11";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_exit_code = false;
    found_exit_code = (found_exit_code || vasp_warnings["EXIT CODE: 11"]);
    found_exit_code = (found_exit_code || vasp_warnings["KILLED BY SIGNAL: 11"]);
    found_warning = (found_warning && (found_bad_termination && found_exit_code));
    xwarning.flag(scheme, found_warning);
    //
    scheme = "MPICH139";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_exit_code = false;
    found_exit_code = (found_exit_code || vasp_warnings["EXIT CODE: 139"]);
    found_exit_code = (found_exit_code || vasp_warnings["KILLED BY SIGNAL: 139"]);
    found_warning = (found_warning && (found_bad_termination && found_exit_code));
    xwarning.flag(scheme, found_warning);
    //
    scheme = "MPICH174"; // CO20210315
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_exit_code = false;
    found_exit_code = (found_exit_code || vasp_warnings["EXIT CODE: 174"]);
    found_exit_code = (found_exit_code || vasp_warnings["KILLED BY SIGNAL: 174"]);
    found_warning = (found_warning && (found_bad_termination && found_exit_code));
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NATOMS"; // look for problem for distance
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["distance between some ions is very small"]);
    xwarning.flag(scheme, found_warning);
    //
    // CO+ME20210315 - simplified
    // ME20190620 - Avoid changing NBANDS in the aflow.in file just because VASP throws the warning that NBANDS is changed because of NPAR.
    // However, if you have that warning AND the error that the number of bands is not sufficient, aflow needs to act.
    scheme = "NBANDS";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning
                     && (vasp_warnings["NBANDS"]
                         && (vasp_warnings["number of bands is not sufficient to hold all electrons"]
                             || vasp_warnings["Number of bands NBANDS too small to hold electrons"]
                             || vasp_warnings["Your highest band is occupied at some k-points"]))); // The NBANDS warning due to NPAR is not an error we want to fix, so set to false if found
    const bool vasp_corrected = (vasp_warnings["number of bands has been changed from the values supplied"] || vasp_warnings["now  NBANDS  ="]); // Need explicit check or else the NPAR warning prevents this NBANDS error from being corrected
    found_warning = (found_warning && vasp_corrected == false);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NPAR";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["please rerun with NPAR="]); // not only npar==1
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NPARC";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings["NPAR = 4"] && vasp_warnings["NPAR=number of cores"])); // fix with NPAR=cores in MPI
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NPARN";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings["NPAR = 4"] && vasp_warnings["NPAR=number of nodes"])); // fix with NPAR=nodes in MPI
    xwarning.flag(scheme, found_warning);
    //
    scheme = "NPAR_REMOVE";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["Please remove the tag NPAR from the INCAR file"]); // and restart the
    xwarning.flag(scheme, found_warning);
    //
    scheme = "PSMAXN";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    // found_warning=(found_warning && vasp_warnings["WARNING: PSMAXN for non-local potential too small"]);
    found_warning = (found_warning && vasp_warnings["REAL_OPT: internal ERROR"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "READ_KPOINTS_RD_SYM";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && (vasp_warnings["ERROR in RE_READ_KPOINTS_RD"] && vasp_warnings["switch off symmetry"]));
    xwarning.flag(scheme, found_warning);
    //
    scheme = "REAL_OPT";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["REAL_OPT: internal ERROR"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "REAL_OPTLAY_1";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["REAL_OPTLAY: internal error (1)"]);
    xwarning.flag(scheme, found_warning);
    //
    scheme = "ZPOTRF";
    found_warning = ReachedAccuracy2bool(scheme, xRequiresAccuracy, xmessage, vasp_still_running);
    found_warning = (found_warning && vasp_warnings["LAPACK: Routine ZPOTRF failed"]);
    xwarning.flag(scheme, found_warning);
    // ALL OTHERS (in alphabetic order) END

    // CO20210315 - this bit must be done before wdebug
    // CO20210315 - only ignore KKSYM if OUTCAR is not registered as incomplete
    if ((xwarning.flag("KKSYM")) && (xwarning.flag("OUTCAR_INCOMPLETE")) && ((ispin_current == 2 && isym_current == -1) || (ispin_current == 1 && isym_current == 0))) { // CO20200624 - needs to change if we do magnetic systems //ispind_current==2 &&
      xmonitor.flag("IGNORING_WARNINGS:KKSYM", false);
    }

    const bool wdebug = (false && LDEBUG);
    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " printing warnings" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    if (!xmonitor.flag("IGNORING_MESSAGES:REACHED_ACCURACY") && (wdebug || xmessage.flag("REACHED_ACCURACY"))) {
      aus << "MMMMM  MESSAGE xmessage.flag(\"REACHED_ACCURACY\")=" << xmessage.flag("REACHED_ACCURACY") << endl;
    }
    if (wdebug) {
      aus << "MMMMM  MESSAGE VASP_release=" << xmessage.getattachedscheme("SVERSION") << endl;
    }
    if (wdebug) {
      aus << "MMMMM  MESSAGE VASP_version=" << xmessage.getattachedscheme("DVERSION") << endl;
    }
    if (wdebug) {
      aus << "MMMMM  MESSAGE AFLOW_version=" << AFLOW_VERSION << endl;
    }
    if (wdebug || xwarning.flag("OUTCAR_INCOMPLETE")) {
      aus << "WWWWW  WARNING xwarning.flag(\"OUTCAR_INCOMPLETE\")=" << xwarning.flag("OUTCAR_INCOMPLETE") << endl; // CO20201111
    }
    if (wdebug || xwarning.flag("BRMIX")) {
      aus << "WWWWW  WARNING xwarning.flag(\"BRMIX\")=" << xwarning.flag("BRMIX") << endl;
    }
    if (wdebug || xwarning.flag("CSLOSHING")) {
      aus << "WWWWW  WARNING xwarning.flag(\"CSLOSHING\")=" << xwarning.flag("CSLOSHING") << endl;
    }
    if (wdebug || xwarning.flag("CALC_FROZEN")) {
      aus << "WWWWW  WARNING xwarning.flag(\"CALC_FROZEN\")=" << xwarning.flag("CALC_FROZEN") << endl; // CO20210135
    }
    if (wdebug || xwarning.flag("DAV")) {
      aus << "WWWWW  WARNING xwarning.flag(\"DAV\")=" << xwarning.flag("DAV") << endl;
    }
    if (wdebug || xwarning.flag("DENTET")) {
      aus << "WWWWW  WARNING xwarning.flag(\"DENTET\")=" << xwarning.flag("DENTET") << endl;
    }
    if (wdebug || xwarning.flag("EDDDAV")) {
      aus << "WWWWW  WARNING xwarning.flag(\"EDDDAV\")=" << xwarning.flag("EDDDAV") << endl;
    }
    if (wdebug || xwarning.flag("EDDRMM")) {
      aus << "WWWWW  WARNING xwarning.flag(\"EDDRMM\")=" << xwarning.flag("EDDRMM") << endl;
    }
    if (wdebug || xwarning.flag("EFIELD_PEAD")) {
      aus << "WWWWW  WARNING xwarning.flag(\"EFIELD_PEAD\")=" << xwarning.flag("EFIELD_PEAD") << endl;
    }
    if (wdebug || xwarning.flag("EXCCOR")) {
      aus << "WWWWW  WARNING xwarning.flag(\"EXCCOR\")=" << xwarning.flag("EXCCOR") << endl;
    }
    if (wdebug || xwarning.flag("GAMMA_SHIFT")) {
      aus << "WWWWW  WARNING xwarning.flag(\"GAMMA_SHIFT\")=" << xwarning.flag("GAMMA_SHIFT") << endl;
    }
    if (!xmonitor.flag("IGNORING_WARNINGS:KKSYM") && !xmonitor.flag("IGNORING_WARNINGS:IBZKPT") && (wdebug || xwarning.flag("IBZKPT"))) {
      aus << "WWWWW  WARNING xwarning.flag(\"IBZKPT\")=" << xwarning.flag("IBZKPT") << endl;
    }
    if (wdebug || xwarning.flag("INVGRP")) {
      aus << "WWWWW  WARNING xwarning.flag(\"INVGRP\")=" << xwarning.flag("INVGRP") << endl;
    }
    if (!xmonitor.flag("IGNORING_WARNINGS:KKSYM") && (wdebug || xwarning.flag("KKSYM"))) {
      aus << "WWWWW  WARNING xwarning.flag(\"KKSYM\")=" << xwarning.flag("KKSYM") << endl;
    }
    if (wdebug || xwarning.flag("LRF_COMMUTATOR")) {
      aus << "WWWWW  WARNING xwarning.flag(\"LRF_COMMUTATOR\")=" << xwarning.flag("LRF_COMMUTATOR") << endl;
    }
    if (wdebug || xwarning.flag("MEMORY")) {
      aus << "WWWWW  WARNING xwarning.flag(\"MEMORY\")=" << xwarning.flag("MEMORY") << endl;
    }
    if (wdebug || xwarning.flag("MPICH0")) {
      aus << "WWWWW  WARNING xwarning.flag(\"MPICH0\")=" << xwarning.flag("MPICH0") << endl;
    }
    if (wdebug || xwarning.flag("MPICH11")) {
      aus << "WWWWW  WARNING xwarning.flag(\"MPICH11\")=" << xwarning.flag("MPICH11") << endl;
    }
    if (wdebug || xwarning.flag("MPICH139")) {
      aus << "WWWWW  WARNING xwarning.flag(\"MPICH139\")=" << xwarning.flag("MPICH139") << endl;
    }
    if (wdebug || xwarning.flag("MPICH174")) {
      aus << "WWWWW  WARNING xwarning.flag(\"MPICH174\")=" << xwarning.flag("MPICH174") << endl;
    }
    if (wdebug || xwarning.flag("NATOMS")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NATOMS\")=" << xwarning.flag("NATOMS") << endl;
    }
    if (wdebug || xwarning.flag("NBANDS")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NBANDS\")=" << xwarning.flag("NBANDS") << endl;
    }
    if (wdebug || xwarning.flag("NELM")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NELM\")=" << xwarning.flag("NELM") << endl;
    }
    if (wdebug || xwarning.flag("NIRMAT")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NIRMAT\")=" << xwarning.flag("NIRMAT") << endl;
    }
    if (wdebug || xwarning.flag("NKXYZ_IKPTD")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NKXYZ_IKPTD\")=" << xwarning.flag("NKXYZ_IKPTD") << endl;
    }
    if (wdebug || xwarning.flag("NPAR")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NPAR\")=" << xwarning.flag("NPAR") << endl;
    }
    if (!xmonitor.flag("IGNORING_WARNINGS:NPARC") && (wdebug || xwarning.flag("NPARC"))) {
      aus << "WWWWW  WARNING xwarning.flag(\"NPARC\")=" << xwarning.flag("NPARC") << endl;
    }
    if (!xmonitor.flag("IGNORING_WARNINGS:NPARN") && (wdebug || xwarning.flag("NPARN"))) {
      aus << "WWWWW  WARNING xwarning.flag(\"NPARN\")=" << xwarning.flag("NPARN") << endl;
    }
    if (wdebug || xwarning.flag("NPAR_REMOVE")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NPAR_REMOVE\")=" << xwarning.flag("NPAR_REMOVE") << endl;
    }
    if (wdebug || xwarning.flag("NUM_PROB")) {
      aus << "WWWWW  WARNING xwarning.flag(\"NUM_PROB\")=" << xwarning.flag("NUM_PROB") << endl; // CO20210315
    }
    if (wdebug || xwarning.flag("PSMAXN")) {
      aus << "WWWWW  WARNING xwarning.flag(\"PSMAXN\")=" << xwarning.flag("PSMAXN") << endl;
    }
    if (wdebug || xwarning.flag("READ_KPOINTS_RD_SYM")) {
      aus << "WWWWW  WARNING xwarning.flag(\"READ_KPOINTS_RD_SYM\")=" << xwarning.flag("READ_KPOINTS_RD_SYM") << endl;
    }
    if (wdebug || xwarning.flag("REAL_OPT")) {
      aus << "WWWWW  WARNING xwarning.flag(\"REAL_OPT\")=" << xwarning.flag("REAL_OPT") << endl;
    }
    if (wdebug || xwarning.flag("REAL_OPTLAY_1")) {
      aus << "WWWWW  WARNING xwarning.flag(\"REAL_OPTLAY_1\")=" << xwarning.flag("REAL_OPTLAY_1") << endl;
    }
    if (wdebug || xwarning.flag("SGRCON")) {
      aus << "WWWWW  WARNING xwarning.flag(\"SGRCON\")=" << xwarning.flag("SGRCON") << endl;
    }
    if (wdebug || xwarning.flag("SYMPREC")) {
      aus << "WWWWW  WARNING xwarning.flag(\"SYMPREC\")=" << xwarning.flag("SYMPREC") << endl;
    }
    if (wdebug || xwarning.flag("ZBRENT")) {
      aus << "WWWWW  WARNING xwarning.flag(\"ZBRENT\")=" << xwarning.flag("ZBRENT") << endl; // CO20210315
    }
    if (wdebug || xwarning.flag("ZPOTRF")) {
      aus << "WWWWW  WARNING xwarning.flag(\"ZPOTRF\")=" << xwarning.flag("ZPOTRF") << endl;
    }
    if (LDEBUG) {
      cerr << aus.str();
    }
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

    // CO20210601 - only print reached accuracy message once, if calc is frozen this can be printed out many times
    if (xmessage.flag("REACHED_ACCURACY")) {
      xmonitor.flag("IGNORING_MESSAGES:REACHED_ACCURACY", true);
    }

    //[CO20210315 - doesn't work]//CO20210315 - this appears to fix ICSD/LIB/CUB/Ag1Sm1_ICSD_604546
    //[CO20210315 - doesn't work]if(xwarning.flag("EDDRMM") && xwarning.flag("ZPOTRF")){ // fix EDDRMM first
    //[CO20210315 - doesn't work]  aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"ZPOTRF\"): prioritizing xwarning.flag(\"EDDRMM\")" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    //[CO20210315 - doesn't work]  xwarning.flag("ZPOTRF",false);
    //[CO20210315 - doesn't work]  //we don't need an xmonitor here, this is only for prioritizing errors
    //[CO20210315 - doesn't work]}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // must do before RMM_DIIS and ROTMAT
    if (true || vasp_monitor_running) { // CO20210315 - might consider running this always, come back and test
      // vasp_monitor will kill vasp prematurely, triggering warnings that require xmessage.flag("REACHED_ACCURACY") (false positives)
      // prioritize the other warnings
      // might be a false positive without vasp_monitor_running, vasp can also trigger its own premature exiting
      // in the future, if these must be turned off, be careful of switching between DAV (EDDDAV) and RMM-DIIS problems
      // one requires ALGO=NORMAL, the other =VERYFAST and they are mutually exclusive
      // therefore, create an xwarnings_fixed which stores which warnings have been fixed previously
      // if warnings like RMM-DIIS have been fixed before, are not problems now, and we encounter CSLOSHING, we should NOT try =NORMAL
      uint n_require_accuracy = 0;
      for (i = 0; i < xwarning.vxscheme.size(); i++) {
        if (xRequiresAccuracy.flag(xwarning.vxscheme[i])) {
          n_require_accuracy++;
        }
      }
      uint n_derivative = 0;
      // these warnings are derivative: e.g., an incomplete outcar could be the result of many errors, including those requiring xmessage.flag("REACHED_ACCURACY")
      const vector<string> vwarnings_derivative{"OUTCAR_INCOMPLETE", "CALC_FROZEN"};
      for (i = 0; i < vwarnings_derivative.size(); i++) {
        if (xwarning.flag(vwarnings_derivative[i])) {
          n_derivative++;
        }
      }
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " xwarning.vxscheme=" << aurostd::joinWDelimiter(xwarning.vxscheme, ",") << endl;
        aus << __AFLOW_FUNC__ << " xwarning.vxscheme.size()=" << xwarning.vxscheme.size() << endl;
        aus << __AFLOW_FUNC__ << " n_require_accuracy=" << n_require_accuracy << endl;
        aus << __AFLOW_FUNC__ << " n_derivative=" << n_derivative << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      if (xwarning.vxscheme.size() > (n_require_accuracy + n_derivative)) { // this means we have some real errors inside
        vector<string> xwarning_vxscheme = xwarning.vxscheme; // make a copy since we're deleting entries of the vector
        for (i = xwarning_vxscheme.size() - 1; i < xwarning_vxscheme.size(); i--) { // go backwards since we're removing entries
          if (xRequiresAccuracy.flag(xwarning_vxscheme[i])) {
            aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"" + xwarning_vxscheme[i] + R"("): prioritizing other warnings first (requires xmessage.flag("REACHED_ACCURACY")))" << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET); //; possible false positive
            xwarning.flag(xwarning_vxscheme[i], false);
            // we don't need an xmonitor here, this is only for prioritizing errors
          }
        }
      }
    }
    // NUM_PROB and ZBRENT are soft errors, only fix if the OUTCAR is incomplete (and not if accuracy not reached)
    // the ZBRENT patch (changing IBRION) may be beneficial
    // conjugate gradient might fail too close to the minimum, so this patch might be good for relax2
    // however, it has shown to steer other calculations in bad directions, leading to other warnings that cannot be fixed
    // better not to over-correct
    // better to check for convergence of the calculation later with --xplug
    if (xwarning.flag("OUTCAR_INCOMPLETE") == false) {
      if (xwarning.flag("NUM_PROB")) {
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NUM_PROB\"): OUTCAR is complete (soft warning)" << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xwarning.flag("NUM_PROB", false);
      }
      if (xwarning.flag("ZBRENT")) {
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"ZBRENT\"): OUTCAR is complete (soft warning)" << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xwarning.flag("ZBRENT", false);
      }
    }

    // CO20210315 - only ignore KKSYM if OUTCAR is not registered as incomplete
    if ((xwarning.flag("KKSYM")) && ((ispin_current == 2 && isym_current == -1) || (ispin_current == 1 && isym_current == 0))) { // CO20200624 - needs to change if we do magnetic systems //ispind_current==2 &&
      if (!xmonitor.flag("IGNORING_WARNINGS:KKSYM")) {
        aus << "MMMMM  MESSAGE ignoring KKSYM warning: ISYM=" << isym_current << " ISPIN=" << ispin_current << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xmonitor.flag("IGNORING_WARNINGS:KKSYM", true); // so we don't clog output files
      }
      xwarning.flag("KKSYM", false);
      if (xwarning.flag("IBZKPT") && xwarning.flag("OUTCAR_INCOMPLETE") == false) { // goes with KKSYM
        aus << "MMMMM  MESSAGE ignoring IBZKPT warning associated with KKSYM: ISYM=" << isym_current << " ISPIN=" << ispin_current << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xwarning.flag("IBZKPT", false);
      }
    }

    // do just before RMM_DIIS and ROTMAT
    // we generally treat CALC_FROZEN as a out-of-memory error UNLESS the calculation has reached its accuracy (finished) and the OUTCAR is incomplete (odd occurrence, noticed on qrats so far)
    // only convert CALC_FROZEN to MEMORY if it's the only warning (OUTCAR_INCOMPLETE is derivative)
    if (xwarning.flag("CALC_FROZEN") && xmessage.flag("REACHED_ACCURACY") == false) {
      if ((xwarning.vxscheme.size() == 1) || (xwarning.flag("OUTCAR_INCOMPLETE") && xwarning.vxscheme.size() == 2)) {
        xwarning.flag("CALC_FROZEN", false);
        xwarning.flag("MEMORY", true);
        aus << R"(MMMMM  MESSAGE treating xwarning.flag("CALC_FROZEN") as xwarning.flag("MEMORY"))" << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET); //; possible false positive
      }
    }

    // the memory might ramp up too quickly, in which the OOM killer will kill vasp
    // if we are lucky, we'll get the "BAD TERMINATION..." warning in the vasp.out (works on qrats)
    if (xwarning.flag("MPICH0") && xmessage.flag("REACHED_ACCURACY") == false) {
      if ((xwarning.vxscheme.size() == 1) || (xwarning.flag("OUTCAR_INCOMPLETE") && xwarning.vxscheme.size() == 2)) {
        xwarning.flag("MPICH0", false);
        xwarning.flag("MEMORY", true);
        aus << R"(MMMMM  MESSAGE treating xwarning.flag("MPICH0") as xwarning.flag("MEMORY"))" << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET); //; possible false positive
      }
    }

    // do last
    bool rmm_diis_warning = false;
    rmm_diis_warning = (rmm_diis_warning || xwarning.flag("EDDRMM"));
    rmm_diis_warning = (rmm_diis_warning || xwarning.flag("NUM_PROB"));
    rmm_diis_warning = (rmm_diis_warning || xwarning.flag("ZBRENT"));
    // CO20210315 - can probably add others to this list as well
    xwarning.flag("RMM_DIIS", rmm_diis_warning);

    bool rotmat_warning = false;
    rotmat_warning = (rotmat_warning || xwarning.flag("SGRCON"));
    rotmat_warning = (rotmat_warning || xwarning.flag("NIRMAT"));
    rotmat_warning = (rotmat_warning || xwarning.flag("IBZKPT"));
    rotmat_warning = (rotmat_warning || xwarning.flag("KKSYM"));
    rotmat_warning = (rotmat_warning || xwarning.flag("INVGRP"));
    rotmat_warning = (rotmat_warning || xwarning.flag("SYMPREC"));
    // CO20210315 - can probably add others to this list as well
    xwarning.flag("ROTMAT", rotmat_warning);
    ///////////////////////////////////////////////////////////////////////////////////////////////

    if (wdebug || xwarning.flag("RMM_DIIS")) {
      aus << "WWWWW  WARNING xwarning.flag(\"RMM_DIIS\")=" << xwarning.flag("RMM_DIIS") << endl; // CO20210315
    }
    if (wdebug || xwarning.flag("ROTMAT")) {
      aus << "WWWWW  WARNING xwarning.flag(\"ROTMAT\")=" << xwarning.flag("ROTMAT") << endl;
    }
    if (LDEBUG) {
      cerr << aus.str();
    }
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

    // WARNINGS STOP

    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " [3]" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    // decide which warnings to ignore

    if (xwarning.flag("MPICH11") && xwarning.flag("NBANDS")) { // fix MPICH11 first
      aus << R"(MMMMM  MESSAGE ignoring xwarning.flag("NBANDS"): prioritizing xwarning.flag("MPICH11"))" << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      xwarning.flag("NBANDS", false);
      // we don't need an xmonitor here, this is only for prioritizing errors
    }
    if (xwarning.flag("NPARC")
        && (aurostd::kvpair2utype<int>(xvasp.INCAR, "NPAR", "=") == 2
            || // dont bother for small NPAR
            aurostd::kvpair2string(xvasp.INCAR, "LRPA", "=") == ".true."
            || // CO20210315 - would be better to check if .true. or .false.
            aurostd::kvpair2string(xvasp.INCAR, "LEPSILON", "=") == ".true."
            || // CO20210315 - would be better to check if .true. or .false.
            aurostd::kvpair2string(xvasp.INCAR, "LOPTICS", "=") == ".true.")) { // dont touch NPARC if LRPA or LEPSILON or LOPTICS necessary
      if (!xmonitor.flag("IGNORING_WARNINGS:NPARC")) {
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NPARC\"): found either LRPA, LEPSILON, or LOPTICS" << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xmonitor.flag("IGNORING_WARNINGS:NPARC", true); // so we don't clog output files
      }
      xwarning.flag("NPARC", false);
    }
    if (xwarning.flag("NPARN")
        && (aurostd::kvpair2string(xvasp.INCAR, "LRPA", "=") == ".true."
            || // CO20210315 - would be better to check if .true. or .false.
            aurostd::kvpair2string(xvasp.INCAR, "LEPSILON", "=") == ".true."
            || // CO20210315 - would be better to check if .true. or .false.
            aurostd::kvpair2string(xvasp.INCAR, "LOPTICS", "=") == ".true.")) { // dont touch NPARN if LRPA or LEPSILON or LOPTICS necessary
      if (!xmonitor.flag("IGNORING_WARNINGS:NPARN")) {
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NPARN\"): found either LRPA, LEPSILON, or LOPTICS" << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        xmonitor.flag("IGNORING_WARNINGS:NPARN", true); // so we don't clog output files
      }
      xwarning.flag("NPARN", false);
    }
    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " [4]" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Error2Fix(const string& error, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    int submode = 0; // default
    const bool try_last_ditch_efforts = true; // default
    return VASP_Error2Fix(error, error, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, const string& mode, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    int submode = 0; // default
    const bool try_last_ditch_efforts = true; // default
    return VASP_Error2Fix(error, mode, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, int& submode, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    const bool try_last_ditch_efforts = true; // default
    return VASP_Error2Fix(error, error, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, const string& mode, int& submode, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    const bool try_last_ditch_efforts = true; // default
    return VASP_Error2Fix(error, mode, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, bool try_last_ditch_efforts, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    int submode = 0; // default
    return VASP_Error2Fix(error, error, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, const string& mode, bool try_last_ditch_efforts, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    int submode = 0; // default
    return VASP_Error2Fix(error, mode, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, int& submode, bool try_last_ditch_efforts, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    return VASP_Error2Fix(error, error, submode, try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error, const string& mode, int& submode, bool try_last_ditch_efforts, _xvasp& xvasp, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // CO20210315
    const bool LDEBUG = (false || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    stringstream aus;

    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " [CHECK " << error << " PROBLEMS]" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    bool apply_fix = xwarning.flag(error);
    const bool VERBOSE = (false || XHOST.vflag_control.flag("MONITOR_VASP") == false || LDEBUG);
    if (apply_fix && xfixed.flag("ALL")) {
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " xfixed.flag(\"ALL\")==true: skipping " << error << " fix" << Message(__AFLOW_FILE__, aflags) << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      apply_fix = false;
    }
    if (apply_fix && vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:ALL")) {
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:ALL\")==true: skipping " << error << " fix" << Message(__AFLOW_FILE__, aflags) << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      apply_fix = false;
    }
    if (apply_fix && vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:" + error)) {
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:" << error << "\")==true: skipping " << error << " fix" << Message(__AFLOW_FILE__, aflags) << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      apply_fix = false;
    }
    if (apply_fix && VERBOSE) {
      aus << "MMMMM  MESSAGE attempting to fix ERROR=\"" << error << "\" (mode=\"" << mode << "\")" << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET); // CO20210315 - do not reference submode after KBIN::XVASP_Afix(), (submode>=0?" (SUBMODE="+aurostd::utype2string(submode)+")":"")
    }
    // do not reference submode below KBIN::XVASP_Afix(), as it will have been incremented (perhaps by 2)
    // if KBIN::XVASP_Afix() fails, submode will be restored to original, so it is okay to reference for the LDEBUG statement
    if (apply_fix && !KBIN::XVASP_Afix(mode, submode, try_last_ditch_efforts, xfixed, xvasp, kflags, vflags, aflags, FileMESSAGE)) { // CO20210315
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " KBIN::XVASP_Afix(mode=\"" << mode << "\"" << (submode >= 0 ? ",submode=" + aurostd::utype2string(submode) : "") << ") failed" << Message(__AFLOW_FILE__, aflags) << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      } // if KBIN::XVASP_Afix() fails, submode will be restored to original, so it is okay to reference for the LDEBUG statement
      apply_fix = false;
    }
    if (apply_fix && VERBOSE) {
      aus << "MMMMM  MESSAGE fix applied for ERROR=\"" << error << "\" (mode=\"" << mode << "\")" << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET); // CO20210315 - do not reference submode after KBIN::XVASP_Afix(), (submode>=0?" (SUBMODE="+aurostd::utype2string(submode)+")":"")
    }
    if (apply_fix) {
      if (xvasp.aopts.flag("FLAG::AFIX_DRYRUN") == false) {
        KBIN::VASP_Error(xvasp, "WWWWW  ERROR " + __AFLOW_FUNC__ + " \"" + error + "\" problems" + Message(__AFLOW_FILE__, aflags));
        //[CO20210315 - old style]xfixed.flag(error,true);
      }
      xfixed.flag("ALL", true);
    }
    return apply_fix;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_FixErrors(_xvasp& xvasp, aurostd::xoption& xmessage, aurostd::xoption& xwarning, aurostd::xoption& xfixed, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) {
    // a note about the fixes below
    // they generally compound, which I believe is the right approach
    // however, there might be some options which conflict
    // add KBIN::XVASP_INCAR_REMOVE_ENTRY() as necessary
    // check also submode fixes

    bool fixed_applied = false;
    bool try_last_ditch_efforts = true;
    uint i = 0;

    for (i = 0; i < 2 && fixed_applied == false; i++) { // for loop goes twice, once with try_last_ditch_efforts==false, then again with ==true
      try_last_ditch_efforts = (i == 1);

      // check NBANDS/LRF_COMMUTATOR problems immediately
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NBANDS", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      //[CO20210315 - fix previously removed]KBIN::VASP_Error2Fix("LRF_COMMUTATOR",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);

      // fix symmetry issues next
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("GAMMA_SHIFT", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      // ********* APPLY PREFERRED SYMMETRY FIXES ******************  //all of these must come before ROTMAT
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NKXYZ_IKPTD", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE)); // must come before IBZKPT
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("KKSYM", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE)); // must come before IBZKPT
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("IBZKPT", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("SYMPREC", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE)); // must come before INVGRP
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("INVGRP", "SYMPREC", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("SGRCON", "SYMPREC", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      // ********* APPLY GENERIC SYMMETRY FIXES ******************
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("ROTMAT", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));

      // fix MPI/NPAR problems next
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("MPICH11", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("MPICH139", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied
                       || KBIN::VASP_Error2Fix("MPICH174", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags,
                                               FileMESSAGE)); // CO20210315 - testing, exit code 174 looks like an error on the node, basically try rerunning with more memory
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NPAR", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NPARC", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NPARN", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NPAR_REMOVE", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));

      // all other fixes, no priority here (alphabetic order)
      //  ********* APPLY OTHER FIXES ******************
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("BRMIX", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("CSLOSHING", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE)); // must come before NELM
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("DAV", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("DENTET", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("EDDDAV", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("EDDRMM", "RMM_DIIS", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("EFIELD_PEAD", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("EXCCOR", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("MEMORY", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NATOMS", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("PSMAXN", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("REAL_OPTLAY_1", "LREAL", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("REAL_OPT", "LREAL", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("ZPOTRF", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));

      // patch only if above warnings are not patched first
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("NELM", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));

      // apply these last if no other fixes worked
      //  ********* APPLY PREFERRED RMM-DIIS FIXES ******************  //all of these must come before RMM-DIIS
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("ZBRENT", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
      // ********* APPLY GENERIC RMM-DIIS FIXES ******************
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("RMM_DIIS", try_last_ditch_efforts, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));

      //[CO20210315 - do not apply patches for frozen calc]//CO20210315 - do last, fixes assume out-of-memory error
      //[CO20210315 - do not apply patches for frozen calc]fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","MEMORY",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
    }

    // sometimes VASP will die after it prints the REACHED_ACCURACY state, but before the OUTCAR is finished (not sure why)
    // this is rare...
    // might be a threading/NFS issue
    // this means any errors inside that require REACHED_ACCURACY will not be triggered
    // in this case, try restarting the calculation from CONTCAR
    // lowering NCPUS has been shown to work, indicating that this is indeed a threading/mpi issue
    // try from the most relaxed CONTCAR to save time
    // this is NOT a magic bullet, it looks like the threading solution works for some structures and not others
    // I am leaving "THREADS" vs. going to "MEMORY" solutions which will change NBANDS, KPOINTS, etc.
    // better to run on another machine/different binary
    if (fixed_applied == false && xwarning.flag("CALC_FROZEN") && xmessage.flag("REACHED_ACCURACY") && xwarning.flag("OUTCAR_INCOMPLETE")) {
      //[CO20210315 - not shown to work]fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","RESTART_CALC",false,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      //[CO20210621 - not shown to work (alone)]fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","RECYCLE_CONTCAR",false,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied = (fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN", "THREADS", false, xvasp, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE));
    }

    // print out all xfixed BEFORE adding "ALL"
    std::sort(xfixed.vxscheme.begin(), xfixed.vxscheme.end()); // sort for printing
    // print ALL first
    stringstream aus;
    if (xfixed.flag("ALL")) {
      aus << "MMMMM  MESSAGE xfixed.flag(\"ALL\")=" << xfixed.flag("ALL") << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    if (xfixed.flag("ALL") || XHOST.vflag_control.flag("MONITOR_VASP") == false) { // very important that we print all xfixed even if not "ALL" for LOCK file, --monitor_vasp reads the LOCK looking here. otherwise, only print if "ALL"
      for (size_t i = 0; i < xfixed.vxscheme.size(); i++) {
        if (xfixed.vxscheme[i] == "ALL") {
          continue;
        } // already done above
        aus << "MMMMM  MESSAGE xfixed.flag(\"" + xfixed.vxscheme[i] + "\")=" << xfixed.flag(xfixed.vxscheme[i]) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
    }

    return fixed_applied;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Run(_xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMESSAGE) { // AFLOW_FUNCTION_IMPLEMENTATION
    const bool LDEBUG = (false || _DEBUG_KVASP_ || XHOST.DEBUG);
    const string function = "KBIN::VASP_Run";
    ostringstream aus_exec;
    ostringstream aus;

    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " BEGIN" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    if (XHOST.AVOID_RUNNING_VASP) { // CO20200624
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "VASP should NOT be running", _INPUT_ILLEGAL_); // better to throw to avoid VASP_Backup(), etc.
      // return false;
    }

    xoption xwarning;
    xoption xfixed;
    xoption xmessage;
    bool vasp_start = true;
    aurostd::StringstreamClean(aus_exec);
    aurostd::StringstreamClean(aus);
    int nrun = 0;
    const int maxrun = 100; // CO20210315 - increase from 15 to 100 //NBANDS can take a lot of iterations to reach goal

    // get CPUS from PBS/SLURM
    // string ausenv;
    aus << "DDDDD  PBS_NUM_PPN=" << XHOST.PBS_NUM_PPN << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  PBS_NNODES=" << XHOST.PBS_NNODES << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  SLURM_CPUS_ON_NODE=" << XHOST.SLURM_CPUS_ON_NODE << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  SLURM_NNODES=" << XHOST.SLURM_NNODES << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  SLURM_NTASKS=" << XHOST.SLURM_NTASKS << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  kflags.KBIN_MPI_NCPUS=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  XHOST.CPU_Cores=" << XHOST.CPU_Cores << Message(__AFLOW_FILE__, aflags) << endl;
    aus << "DDDDD  aflags.AFLOW_GLOBAL_NCPUS=" << aflags.AFLOW_GLOBAL_NCPUS << Message(__AFLOW_FILE__, aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);

    // for reducint CPUs on the fly
    if (aflags.AFLOW_GLOBAL_NCPUS < 0) {
      kflags.KBIN_MPI_NCPUS = -aflags.AFLOW_GLOBAL_NCPUS; // this to force things on reducing CPUS
    }

    aus << "DDDDD  kflags.KBIN_MPI_NCPUS=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__, aflags) << endl;

    // for for LS coupling
    if (vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option) {
      if (!aurostd::substring2bool(kflags.KBIN_BIN, VASPLS_BIN_POSTFIX_DEFAULT)) {
        kflags.KBIN_BIN = kflags.KBIN_BIN + VASPLS_BIN_POSTFIX_DEFAULT; // standard LS
      }
      if (!aurostd::substring2bool(kflags.KBIN_MPI_BIN, VASPLS_BIN_POSTFIX_DEFAULT)) {
        kflags.KBIN_MPI_BIN = kflags.KBIN_MPI_BIN + VASPLS_BIN_POSTFIX_DEFAULT; // standard LS
      }
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      aus << "00000  MESSAGE SPIN-ORBIT TYPE CALCULATIONS , adding " << VASPLS_BIN_POSTFIX_DEFAULT << " to BIN" << Message(__AFLOW_FILE__, aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }
    if (kflags.KBIN_MPI) {
      kflags.KBIN_BIN = kflags.KBIN_MPI_BIN; // forcing, no matter what
    }

    uint xvasp_aopts_vxscheme_size = 0;
    uint vflags_KBIN_VASP_FORCE_OPTION_IGNORE_AFIX_vxscheme_size = 0;

    while (vasp_start) {
      // ********* RUN VASP
      { // ERRORS
        bool error = false;
        if (aurostd::FileEmpty(xvasp.Directory + "/INCAR")) {
          KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR " + function + ": Empty INCAR" + Message(__AFLOW_FILE__, aflags));
          error = true;
          return false;
        }
        if (aurostd::FileEmpty(xvasp.Directory + "/POSCAR")) {
          KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR " + function + ": Empty POSCAR" + Message(__AFLOW_FILE__, aflags));
          error = true;
          return false;
        }
        if (aurostd::FileEmpty(xvasp.Directory + "/KPOINTS")) {
          KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR " + function + ": Empty KPOINTS" + Message(__AFLOW_FILE__, aflags));
          error = true;
          return false;
        }
        if (aurostd::FileEmpty(xvasp.Directory + "/POTCAR")
            && (KBIN::VASP_Produce_POTCAR(xvasp, aurostd::file2string(xvasp.Directory + "/" + _AFLOWIN_), FileMESSAGE, aflags, kflags, vflags) ? !aurostd::stringstream2file(xvasp.POTCAR, xvasp.Directory + "/POTCAR") : true)) {
          KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR " + function + ": Empty POTCAR" + Message(__AFLOW_FILE__, aflags));
          error = true;
          return false;
        }
        if (error) {
          return false;
        }

        if (LDEBUG) {
          aus << __AFLOW_FUNC__ << " [1]" << Message(__AFLOW_FILE__, aflags) << endl;
          cerr << aus.str();
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }

        // FIX INCAR if alternating
        if (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")) {
          if (aurostd::_isodd(xvasp.NRELAXING)) {
            aus << "00000  MESSAGE Alternating: RELAX_CELL_VOLUME" << Message(__AFLOW_FILE__, aflags) << endl;
          }
          if (aurostd::_iseven(xvasp.NRELAXING)) {
            aus << "00000  MESSAGE Alternating: RELAX_IONS" << Message(__AFLOW_FILE__, aflags) << endl;
          }
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC", false);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL", false);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS", false);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE", false);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME", false);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_VOLUME");
          KBIN::XVASP_INCAR_Relax_ON(xvasp, vflags, xvasp.NRELAXING);
        }

        // CO20210315
        // print out these schemes so they can be picked up by the vasp monitor
        xvasp_aopts_vxscheme_size = xvasp.aopts.vxscheme.size();
        for (uint i = 0; i < xvasp_aopts_vxscheme_size; i++) {
          const string& flag = xvasp.aopts.vxscheme[i];
          if (flag.find("FLAG::") != string::npos && flag.find("_PRESERVED") != string::npos) {
            if (xvasp.aopts.flag(flag)) {
              aus << "MMMMM  MESSAGE xvasp.aopts.flag(\"" << flag << "\")=" << xvasp.aopts.flag(flag) << endl;
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
          }
        }
        vflags_KBIN_VASP_FORCE_OPTION_IGNORE_AFIX_vxscheme_size = vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme.size();
        for (uint i = 0; i < vflags_KBIN_VASP_FORCE_OPTION_IGNORE_AFIX_vxscheme_size; i++) {
          const string& flag = vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme[i];
          if (vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(flag)) {
            aus << "MMMMM  MESSAGE vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"" << flag << "\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(flag) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
        }
        if (vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved) {
          aus << "MMMMM  MESSAGE vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved=" << vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }

        // RUN VASP NON QUEUE ------------------------------------------------------------------------
        if (kflags.KBIN_QSUB == false) {
          nrun++;
          aus_exec << "cd " << xvasp.Directory << endl;
          aurostd::RemoveFile(xvasp.Directory + "/" + DEFAULT_VASP_OUT);
          if (kflags.KBIN_MPI == false) {
            aus << "00000  MESSAGE SERIAL job - [" << xvasp.str.atoms.size() << "atoms]" << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            aus_exec << kflags.KBIN_BIN << " > " << DEFAULT_VASP_OUT << endl;
            aus << "00000  MESSAGE" << VASP_KEYWORD_EXECUTION << kflags.KBIN_BIN << " > " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory") << endl; // CO20170628 - SLOW WITH MEMORY
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            aurostd::execute(aus_exec);
            aurostd::Sleep(_KVASP_VASP_SLEEP_);
          } else {
            aus << "00000  MESSAGE MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__, aflags) << endl;
            if (!kflags.KBIN_MPI_OPTIONS.empty()) {
              aus << "00000  MESSAGE MPI OPTIONS=[" << kflags.KBIN_MPI_OPTIONS << "]" << Message(__AFLOW_FILE__, aflags) << endl;
            }
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            if (LDEBUG) {
              aus << __AFLOW_FUNC__ << " aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl;
              aus << __AFLOW_FUNC__ << " aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
              cerr << aus.str();
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            }
            // NO HOST ------------------------------------------------------------------------
            if (!aflags.AFLOW_MACHINE_LOCAL.flag()) {
              aus << "00000  MESSAGE" << VASP_KEYWORD_EXECUTION;
              if (!kflags.KBIN_MPI_OPTIONS.empty()) {
                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                //[CO20210315 - do not print, will confuse vasp monitor]aus << kflags.KBIN_MPI_OPTIONS << "; ";
              }
              if (!kflags.KBIN_MPI_START.empty()) {
                aus_exec << kflags.KBIN_MPI_START << " > " << DEFAULT_VASP_OUT << endl;
                //[CO20210315 - do not print, will confuse vasp monitor]aus << kflags.KBIN_MPI_START << " > " << DEFAULT_VASP_OUT << "; ";
              }
              aus_exec << kflags.KBIN_MPI_COMMAND << " " << kflags.KBIN_MPI_NCPUS << " " << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aus << kflags.KBIN_MPI_COMMAND << " " << kflags.KBIN_MPI_NCPUS << " " << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT; //[CO20210315 - do not print, will confuse vasp monitor]<< "; ";
              if (!kflags.KBIN_MPI_STOP.empty()) {
                aus_exec << kflags.KBIN_MPI_STOP << " >> " << DEFAULT_VASP_OUT << endl;
                //[CO20210315 - do not print, will confuse vasp monitor]aus << kflags.KBIN_MPI_STOP << " >> " << DEFAULT_VASP_OUT;
              }
              aus << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory") << endl; // CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_BETA_MPICH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_BETA_MPICH
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_BETA_MPICH
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_BETA_MPICH << endl;
              aus_exec << MPI_COMMAND_DUKE_BETA_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_MPICH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_BETA_OPENMPI ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
              if (!aurostd::substring2bool(kflags.KBIN_MPI_BIN, "_openmpi")) {
                kflags.KBIN_MPI_BIN = kflags.KBIN_MPI_BIN + "_openmpi"; // fix the OPENMPI
              }
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_BETA_OPENMPI
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_BETA_OPENMPI
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_BETA_OPENMPI << endl;
              aus_exec << MPI_COMMAND_DUKE_BETA_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_QRATS_MPICH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_QRATS_MPICH
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_QRATS_MPICH
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_QRATS_MPICH << endl;
              aus_exec << MPI_COMMAND_DUKE_QRATS_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QRATS_MPICH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_QFLOW_OPENMPI ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_QFLOW_OPENMPI
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << endl;
              aus_exec << MPI_COMMAND_DUKE_QFLOW_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // CO20201220 X START
            //  HOST DUKE_X_X ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_X")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_X_X
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_X_X
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_X << endl;
              aus_exec << MPI_COMMAND_DUKE_X_X << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_X << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_CRAY")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_X_CRAY
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_X_CRAY
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_CRAY << endl;
              aus_exec << MPI_COMMAND_DUKE_X_CRAY << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_CRAY << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //        aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_OLDCRAY")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_X_OLDCRAY
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_X_OLDCRAY
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_OLDCRAY << endl;
              aus_exec << MPI_COMMAND_DUKE_X_OLDCRAY << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_OLDCRAY << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //        aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_SMB")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_X_SMB
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_X_SMB
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_SMB << endl;
              aus_exec << MPI_COMMAND_DUKE_X_SMB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_SMB << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //        aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // CO20201220 X STOP
            // CO20220818 JHU_ROCKFISH START
            //  HOST JHU_ROCKFISH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::JHU_ROCKFISH")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_JHU_ROCKFISH
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_JHU_ROCKFISH
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_JHU_ROCKFISH << endl;
              aus_exec << MPI_COMMAND_JHU_ROCKFISH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_JHU_ROCKFISH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // CO20220818 JHU_ROCKFISH STOP
            //  HOST MPCDF_EOS_MPI ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {
              // verbosization
              int local_NCPUS = kflags.KBIN_MPI_NCPUS;
              if (MPI_NCPUS_MPCDF_EOS > 0) {
                local_NCPUS = MPI_NCPUS_MPCDF_EOS;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_EOS" << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << local_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MPCDF_EOS
                  << " "
                  << local_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MPCDF_EOS
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_EOS << endl;
              aus_exec << MPI_COMMAND_MPCDF_EOS << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_EOS << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST MPCDF_DRACO_MPI ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {
              // verbosization
              int local_NCPUS = kflags.KBIN_MPI_NCPUS;
              if (MPI_NCPUS_MPCDF_DRACO > 0) {
                local_NCPUS = MPI_NCPUS_MPCDF_DRACO;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_DRACO" << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << local_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MPCDF_DRACO
                  << " "
                  << local_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MPCDF_DRACO
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_DRACO << endl;
              aus_exec << MPI_COMMAND_MPCDF_DRACO << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_DRACO << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST MPCDF_COBRA_MPI ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {
              // verbosization
              int local_NCPUS = kflags.KBIN_MPI_NCPUS;
              if (MPI_NCPUS_MPCDF_COBRA > 0) {
                local_NCPUS = MPI_NCPUS_MPCDF_COBRA;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_COBRA" << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << local_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MPCDF_COBRA
                  << " "
                  << local_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MPCDF_COBRA
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_COBRA << endl;
              aus_exec << MPI_COMMAND_MPCDF_COBRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_COBRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST MPCDF_HYDRA_MPI ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {
              // verbosization
              int local_NCPUS = kflags.KBIN_MPI_NCPUS;
              if (MPI_NCPUS_MPCDF_HYDRA > 0) {
                local_NCPUS = MPI_NCPUS_MPCDF_HYDRA;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_HYDRA" << Message(__AFLOW_FILE__, aflags) << endl; // HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              }
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << local_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              //	      aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MPCDF_HYDRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MPCDF_HYDRA
                  << " "
                  << MPI_BINARY_DIR_MPCDF_HYDRA
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name   // poe not MPI run
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_HYDRA << endl;
              //	      aus_exec << MPI_COMMAND_MPCDF_HYDRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aus_exec << MPI_COMMAND_MPCDF_HYDRA << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl; // poe not MPI run
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // DX20190509 - MACHINE001 - START
            //  HOST MACHINE001_MPICH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MACHINE001
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MACHINE001
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE001 << endl;
              aus_exec << MPI_COMMAND_MACHINE001 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE001 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // DX20190509 - MACHINE001 - END
            // DX20190509 - MACHINE002 - START
            //  HOST MACHINE002_MPICH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MACHINE002
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MACHINE002
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE002 << endl;
              aus_exec << MPI_COMMAND_MACHINE002 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE002 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // DX20190509 - MACHINE002 - END
            // DX20201005 - MACHINE003 - START
            //  HOST MACHINE003_MPICH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE003")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MACHINE003
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MACHINE003
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE003 << endl;
              aus_exec << MPI_COMMAND_MACHINE003 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE003 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // DX20201005 - MACHINE003 - END
            // DX20211011 - MACHINE004 - START
            //  HOST MACHINE004_MPICH ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE004")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs  "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MACHINE004
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_MACHINE004
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE004 << endl;
              aus_exec << MPI_COMMAND_MACHINE004 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE004 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // DX20211011 - MACHINE004 - END
            //  HOST DUKE_MATERIALS ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_MATERIALS
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_MATERIALS
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_MATERIALS << endl;
              aus_exec << MPI_COMMAND_DUKE_MATERIALS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_MATERIALS << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_AFLOWLIB ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_AFLOWLIB
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_AFLOWLIB
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_AFLOWLIB << endl;
              aus_exec << MPI_COMMAND_DUKE_AFLOWLIB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_AFLOWLIB << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_HABANA ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_DUKE_HABANA
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_DUKE_HABANA
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_HABANA << endl;
              aus_exec << MPI_COMMAND_DUKE_HABANA << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_HABANA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST FULTON_MARYLOU ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              //	      aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_FULTON_MARYLOU << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_FULTON_MARYLOU
                  << " "
                  << MPI_BINARY_DIR_FULTON_MARYLOU
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_FULTON_MARYLOU << endl;
              aus_exec << MPI_COMMAND_FULTON_MARYLOU << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              // aus_exec << MPI_COMMAND_FULTON_MARYLOU << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;  // with --np
              aurostd::execute(aus_exec);
            }
            // HOST CMU_EULER ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_CMU_EULER
                  << " "
                  << kflags.KBIN_MPI_NCPUS
                  << " "
                  << MPI_BINARY_DIR_CMU_EULER
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_CMU_EULER << endl;
              aus_exec << MPI_COMMAND_CMU_EULER << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_CMU_EULER << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST OL ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MACHINE2
                  << " "
                  << MPI_BINARY_DIR_MACHINE2
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE2 << endl; // CO20181226
              aus_exec << MPI_COMMAND_MACHINE2 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE2 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl; // CO20181226 - adding kflags.KBIN_MPI_NCPUS
              aurostd::execute(aus_exec);
            }
            // HOST HOST1 ------------------------------------------------------------------------
            if (aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) {
              // verbosization
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << "  MPI PARALLEL job - ["
                  << xvasp.str.atoms.size()
                  << "atoms] - "
                  << " MPI="
                  << kflags.KBIN_MPI_NCPUS
                  << "CPUs "
                  << Message(__AFLOW_FILE__, aflags)
                  << endl; // HE20220309 use machine name
              aus
                  << "00000  MESSAGE HOST="
                  << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME")
                  << " "
                  << VASP_KEYWORD_EXECUTION
                  << MPI_COMMAND_MACHINE1
                  << " "
                  << MPI_BINARY_DIR_MACHINE1
                  << kflags.KBIN_MPI_BIN
                  << " >> "
                  << DEFAULT_VASP_OUT
                  << Message(__AFLOW_FILE__, aflags, string(_AFLOW_MESSAGE_DEFAULTS_) + ",memory")
                  << endl; // HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE1 << endl; // CO20181226
              aus_exec << MPI_COMMAND_MACHINE1 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE1 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl; // CO20181226 - adding kflags.KBIN_MPI_NCPUS
              aurostd::execute(aus_exec);
            }
            // DONE ------------------------------------------------------------------------
          }
          aurostd::Sleep(_KVASP_VASP_SLEEP_);
          vasp_start = false;
        }
        // RUN VASP QUEUED ------------------------------------------------------------------------
        if (kflags.KBIN_QSUB) {
          nrun++;
          aus_exec << "cd " << xvasp.Directory << endl;
          if (kflags.KBIN_MPI == false) {
            aus << "00000  MESSAGE QUEUED SERIAL job - [" << xvasp.str.atoms.size() << "atoms]" << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            aurostd::RemoveFile(string(xvasp.Directory + "/aflow.qsub.done"));
            aurostd::stringstream2file(xvasp.xqsub.QSUB, string(xvasp.Directory + "/aflow.qsub.run"));
            aurostd::Chmod(0755, string(xvasp.Directory + "/aflow.qsub.run"));
            aus_exec << kflags.KBIN_QSUB_COMMAND << " " << kflags.KBIN_QSUB_PARAMS << " " << "./aflow.qsub.run &" << endl;
            aurostd::execute(aus_exec);
            KBIN::QSUB_WaitFinished(aflags, FileMESSAGE, false);
            aurostd::RemoveFile(string(xvasp.Directory + "/aflow.qsub.done"));
          } else {
            aus << "00000  MESSAGE QUEUED MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
            aurostd::RemoveFile(string(xvasp.Directory + "/aflow.qsub.done"));
            aurostd::stringstream2file(xvasp.xqsub.QSUB, string(xvasp.Directory + "/aflow.qsub.run"));
            aurostd::Chmod(0755, string(xvasp.Directory + "/aflow.qsub.run"));
            aus_exec << kflags.KBIN_QSUB_COMMAND << " " << kflags.KBIN_QSUB_PARAMS << " " << "./aflow.qsub.run &" << endl;
            aurostd::execute(aus_exec);
            KBIN::QSUB_WaitFinished(aflags, FileMESSAGE, false);
            aurostd::RemoveFile(string(xvasp.Directory + "/aflow.qsub.done"));
          }
          aurostd::Sleep(_KVASP_VASP_SLEEP_);
          vasp_start = false;
        }
      }
      KBIN::WaitFinished(xvasp, aflags, FileMESSAGE, 2, false); // CO20201111 - try twice and NO verbose, we verbose in bigger VASP_Run() loop
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " [2]" << Message(__AFLOW_FILE__, aflags) << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }

      if (aurostd::FileEmpty(xvasp.Directory + "/" + DEFAULT_VASP_OUT)) {
        KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR " + function + ": Empty " + DEFAULT_VASP_OUT + Message(__AFLOW_FILE__, aflags));
        return false;
      }
      if (aurostd::FileEmpty(xvasp.Directory + "/OUTCAR")) {
        KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR " + function + ": Empty OUTCAR" + Message(__AFLOW_FILE__, aflags));
        return false;
      }
      // DONT CHECK CONTCAR it can be empty
      // DONT CHECK OSZICAR it can be empty

      // update kpoints table

      // ***************** CHECK FOR ERRORS *********
      if (LDEBUG) {
        aus << __AFLOW_FUNC__ << " [3a]  nrun=" << nrun << Message(__AFLOW_FILE__, aflags) << endl;
        aus << __AFLOW_FUNC__ << " [3b]  maxrun=" << maxrun << Message(__AFLOW_FILE__, aflags) << endl;
        cerr << aus.str();
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }

      xmessage.clear(); // CO20210315

      // check VASP version
      const string SVERSION = KBIN::OUTCAR2VASPVersion(xvasp.Directory + "/OUTCAR");
      const double DVERSION = aurostd::VersionString2Double(SVERSION);
      xmessage.push_attached("SVERSION", SVERSION); // CO20210315 - put to xmessage
      xmessage.push_attached("DVERSION", aurostd::utype2string(DVERSION)); // CO20210315 - put to xmessage

      // get algo_current - START
      KBIN::VASP_Reread_INCAR(xvasp);
      string algo_current = "NORMAL"; // vasp default: https://www.vasp.at/wiki/index.php/ALGO
      if (aurostd::substring2bool(xvasp.INCAR, "ALGO=", true)) {
        algo_current = aurostd::toupper(aurostd::kvpair2string(xvasp.INCAR, "ALGO", "="));
      } // CO20210315 - remove whitespaces
      else if (aurostd::substring2bool(xvasp.INCAR, "IALGO=", true)) { // aflow also prints IALGO sometimes, need to check, remove whitespaces
        const string algo_current_tmp = aurostd::toupper(KBIN::INCAR_IALGO2ALGO(aurostd::kvpair2utype<int>(xvasp.INCAR, "IALGO", "=")));
        if (!algo_current_tmp.empty()) {
          algo_current = algo_current_tmp;
        }
      }
      if (!algo_current.empty()) { // so we don't retry later as a fix
        xfixed.flag("ALGO=" + algo_current, true);
        aus << "MMMMM  MESSAGE adding current \"ALGO=" << algo_current << "\" to xfixed" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      // get algo_current - END

      if (nrun < maxrun) {
        if (LDEBUG) {
          aus << __AFLOW_FUNC__ << " [4]" << Message(__AFLOW_FILE__, aflags) << endl;
          cerr << aus.str();
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }

        size_t offset = 0;
        std::unordered_map<string, bool> vasp_warnings;
        VASP_InitWarnings(vasp_warnings);
        VASP_ProcessWarnings(xvasp, aflags, kflags, offset, vasp_warnings, xmessage, xwarning, FileMESSAGE);

        xfixed.flag("ALL", false);
        vasp_start = false;

        if (LDEBUG) {
          aus << __AFLOW_FUNC__ << " [5]" << Message(__AFLOW_FILE__, aflags) << endl;
          cerr << aus.str();
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }

        KBIN::VASP_FixErrors(xvasp, xmessage, xwarning, xfixed, aflags, kflags, vflags, FileMESSAGE);

        // come back - should we check if any xwarning() flag is still on?
        // CO20210315 - NO, aflow patches if possible, otherwise let it run. we'll catch BIG issues with --xplug

        // ********* VASP TO BE RESTARTED *********
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " [DONE WITH CHECKS]" << Message(__AFLOW_FILE__, aflags) << endl;
        }
        if (xfixed.flag("ALL")) {
          vasp_start = true;
        }
        if (vasp_start) {
          if (LDEBUG) {
            aus << __AFLOW_FUNC__ << " [VASP TO BE RESTARTED]" << Message(__AFLOW_FILE__, aflags) << endl;
            cerr << aus.str();
            aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
          }
          aus << "00000  RESTART VASP" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        }
      }
    }

    if (LDEBUG) {
      aus << "MMMMM  MESSAGE tested all the errors" << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    if (!aurostd::FileExist(xvasp.Directory + "/" + DEFAULT_VASP_OUT)) {
      KBIN::VASP_Error(xvasp, "EEEEE  file does not exist=" + DEFAULT_VASP_OUT);
      return false;
    }
    if (aurostd::FileEmpty(xvasp.Directory + "/" + DEFAULT_VASP_OUT)) {
      KBIN::VASP_Error(xvasp, "EEEEE  file empty=" + DEFAULT_VASP_OUT);
      return false;
    }
    if (!aurostd::FileExist(xvasp.Directory + "/OUTCAR")) {
      KBIN::VASP_Error(xvasp, "EEEEE  file does not exist=OUTCAR");
      return false;
    }
    if (aurostd::FileEmpty(xvasp.Directory + "/OUTCAR")) {
      KBIN::VASP_Error(xvasp, "EEEEE  file empty=OUTCAR");
      return false;
    }
    if (!aurostd::FileExist(xvasp.Directory + "/CONTCAR")) {
      KBIN::VASP_Error(xvasp, "EEEEE  file does not exist=CONTCAR");
      return false;
    }
    if (aurostd::FileEmpty(xvasp.Directory + "/CONTCAR")) {
      KBIN::VASP_Error(xvasp, "EEEEE  file is empty=CONTCAR");
      return false;
    }

    if (LDEBUG) {
      aus << __AFLOW_FUNC__ << " END" << Message(__AFLOW_FILE__, aflags) << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    return true;
    // ********* FINISH
    //  return 1;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Run(_xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, string relax, bool qmwrite, ofstream& FileMESSAGE) { // AFLOW_FUNCTION_IMPLEMENTATION
    if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
      KBIN::VASP_Error(xvasp, "EEEEE  ERROR KBIN::VASP_Run" + Message(__AFLOW_FILE__, aflags));
      return false;
    }
    if (_VASP_CONTCAR_SAVE_) {
      KBIN::VASP_CONTCAR_Save(xvasp, string(relax));
    }
    if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
      KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR KBIN::VASP_RunFinished: OUTCAR is incomplete" + Message(__AFLOW_FILE__, aflags)); // CO20201111  //AFTER CONTCAR_SAVE_
      return false;
    }
    KBIN::VASP_Backup(xvasp, qmwrite, relax);
    return true;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Run(_xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, string relaxA, string relaxB, bool qmwrite, ofstream& FileMESSAGE) { // AFLOW_FUNCTION_IMPLEMENTATION
    if (relaxA != relaxB) {
      const string message = "relaxA (" + relaxA + ") != relaxB (" + relaxB + ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if (!KBIN::VASP_Run(xvasp, aflags, kflags, vflags, FileMESSAGE)) {
      KBIN::VASP_Error(xvasp, "EEEEE  ERROR KBIN::VASP_Run" + Message(__AFLOW_FILE__, aflags));
      return false;
    }
    if (_VASP_CONTCAR_SAVE_) {
      KBIN::VASP_CONTCAR_Save(xvasp, string(relaxA));
    }
    if (!KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, true)) {
      KBIN::VASP_Error(xvasp, FileMESSAGE, "EEEEE  ERROR KBIN::VASP_RunFinished: OUTCAR is incomplete" + Message(__AFLOW_FILE__, aflags)); // CO20201111 //AFTER CONTCAR_SAVE_
      return false;
    }
    KBIN::VASP_Backup(xvasp, qmwrite, relaxA);
    KBIN::VASP_Recycle(xvasp, relaxB);
    return true;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_RunFinished(_xvasp& xvasp, _aflags& aflags, ofstream& FileMESSAGE, bool verbose) {
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    // if(verbose) aus << "00000  MESSAGE RUN CHECK FINISHED :" << Message(__AFLOW_FILE__,aflags) << endl;
    // if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OUTCAR DOES NOT EXIST
    if (!aurostd::FileExist(xvasp.Directory + "/OUTCAR")) {
      if (verbose) {
        aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR does not exist) :" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (verbose) {
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      return false;
    }
    // OUTCAR EXISTS BUT EMPTY
    if (aurostd::FileEmpty(xvasp.Directory + "/OUTCAR")) {
      if (verbose) {
        aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR is empty) :" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (verbose) {
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      return false;
    }
    // OUTCAR EXISTS
    const bool RemoveWS = true;
    const bool case_insensitive = true;
    const bool expect_near_end = true;
    if (aurostd::substring_present_file(xvasp.Directory + "/OUTCAR", "Total CPU time used (sec)")) { // CO20210601
      if (verbose) {
        aus << "00000  MESSAGE RUN FINISHED (OUTCAR is complete) :" << Message(__AFLOW_FILE__, aflags) << endl;
      }
      if (verbose) {
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
      }
      return true;
    }
    if (verbose) {
      aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR is incomplete)" << Message(__AFLOW_FILE__, aflags) << endl;
    }
    if (verbose) {
      aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
    }

    return false;
  }
} // namespace KBIN

namespace KBIN {
  void WaitFinished(_xvasp& xvasp, _aflags& aflags, ofstream& FileMESSAGE, uint max_count, bool verbose) {
    uint i = 0;
    uint sleep_seconds = SECONDS_SLEEP_VASP_COMPLETION;
    if ((max_count * sleep_seconds) > SECONDS_SLEEP_VASP_MONITOR) { // safety for --monitor_vasp
      sleep_seconds = (uint) (aurostd::max(3.0, ((double) SECONDS_SLEEP_VASP_MONITOR / (double) max_count) - 5.0)); // the max ensures we don't go below 0 (if SECONDS_SLEEP_VASP_MONITOR is too low)
    }
    while ((i++) < max_count && !KBIN::VASP_RunFinished(xvasp, aflags, FileMESSAGE, verbose)) {
      aurostd::Sleep(sleep_seconds); // CO20201111
    }
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Error(const _xvasp& xvasp, const string& message1, const string& message2, const string& message3) { // CO20210315 - cleaned up
    ofstream FileXVASP;
    return VASP_Error(xvasp, FileXVASP, message1, message2, message3);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Error(const _xvasp& xvasp, ofstream& FileMESSAGE, const string& message1, const string& message2, const string& message3) { // CO20210315 - cleaned up
    const string FileNameXVASP = xvasp.Directory + "/" + DEFAULT_AFLOW_ERVASP_OUT;
    stringstream message;
    message << message1 << message2 << message3 << endl;
    aurostd::stringstream2file(message, FileNameXVASP, aurostd::compression_type::None, "APPEND");
    aurostd::PrintMessageStream(FileMESSAGE, message, XHOST.QUIET);
  }
} // namespace KBIN

namespace KBIN {
  string VASP_Analyze(_xvasp& xvasp, bool qmwrite) { // AFLOW_FUNCTION_IMPLEMENTATION
    // CHECK ERRORS
    bool error = false;
    if (aurostd::FileEmpty(xvasp.Directory + "/CONTCAR")) {
      KBIN::VASP_Error(xvasp, "EEEEE  ERROR " + __AFLOW_FUNC__ + ": Empty CONTCAR");
      error = true;
    }
    if (aurostd::FileEmpty(xvasp.Directory + "/OUTCAR") || !aurostd::substring_present_file(xvasp.Directory + "/OUTCAR", "Total CPU time used (sec)")) {
      KBIN::VASP_Error(xvasp, "EEEEE  ERROR " + __AFLOW_FUNC__ + ": Empty OUTCAR");
      error = true;
    } // SD20221019 - OUTCAR does not exist or is incomplete
    if (aurostd::FileEmpty(xvasp.Directory + "/INCAR")) {
      KBIN::VASP_Error(xvasp, "EEEEE  ERROR " + __AFLOW_FUNC__ + ": Empty INCAR");
      error = true;
    }
    if (aurostd::FileEmpty(xvasp.Directory + "/vasprun.xml")) {
      KBIN::VASP_Error(xvasp, "EEEEE  ERROR " + __AFLOW_FUNC__ + ": Empty vasprun.xml");
      error = true;
    }
    if (error) {
      return "";
    }

    xvasp.str.qm_clear();
    xvasp.str.qm_load(xvasp.Directory);
    // LOAD OUTPUTS
    stringstream strstream;
    strstream.clear();
    strstream.str(std::string());
    strstream.setf(std::ios::fixed, std::ios::floatfield);
    // OUTCAR OPERATIONS ---------------------------------------------------------------
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream << "[KBIN_ANALYZE]START_" << xvasp.AnalyzeLabel << endl;
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream << "# POSITION                                       TOTAL-FORCE (eV/Angst)               " << endl;
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream.precision(_DOUBLE_WRITE_PRECISION_); // CO20200731 - 12
    for (size_t i = 0; i < xvasp.str.atoms.size(); i++) { // clear        (from the previous step)
      for (uint j = 1; j <= 3; j++) {
        if (std::abs(xvasp.str.qm_positions.at(i)[j]) < 10.0) {
          strstream << " ";
        }
        if (xvasp.str.qm_positions.at(i)[j] >= 0.0) {
          strstream << " ";
        }
        strstream << "   " << xvasp.str.qm_positions.at(i)[j] << " ";
      }
      for (uint j = 1; j <= 3; j++) {
        if (std::abs(xvasp.str.qm_forces.at(i)[j]) < 10.0) {
          strstream << " ";
        }
        if (xvasp.str.qm_forces.at(i)[j] >= 0.0) {
          strstream << " ";
        }
        strstream << "   " << xvasp.str.qm_forces.at(i)[j] << " ";
      }
      strstream << endl;
    }
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    // OUTCAR OPERATIONS ---------------------------------------------------------------
    strstream.setf(std::ios::scientific, std::ios::floatfield);
    strstream.setf(std::ios::left, std::ios::adjustfield);
    strstream << "E_cell=" << xvasp.str.qm_E_cell << "  (eV/cell)" << endl;
    strstream << "E_atom=" << xvasp.str.qm_E_atom << "  (eV/at)" << endl;
    strstream << "H_cell=" << xvasp.str.qm_H_cell << "  (eV/cell)" << endl;
    strstream << "H_atom=" << xvasp.str.qm_H_atom << "  (eV/at)" << endl;
    strstream << "PV_cell=" << xvasp.str.qm_PV_cell << "  (eV/cell)" << endl;
    strstream << "PV_atom=" << xvasp.str.qm_PV_atom << "  (eV/at)" << endl;
    strstream << "mag_cell=" << xvasp.str.qm_mag_cell << "  (mu/cell)" << endl;
    strstream << "mag_atom=" << xvasp.str.qm_mag_atom << "  (mu/at)" << endl;
    xstructure qm_str(xvasp.str); // suck it in !
    qm_str.qm_recycle();
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream << "[VASP_POSCAR_MODE_EXPLICIT]START " << endl;
    strstream << qm_str;
    strstream << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << endl;

    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream << "[KBIN_ANALYZE]STOP_" << xvasp.AnalyzeLabel << endl;
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;

    if (qmwrite) {
      const string FileNameXVASP = xvasp.Directory + "/" + DEFAULT_AFLOW_QMVASP_OUT;
      stringstream FileXVASPout;
      // ME20200304 - do not overwrite prior runs
      string FileNameXVASPfull;
      if (aurostd::CompressFileExist(FileNameXVASP, FileNameXVASPfull) && aurostd::IsCompressed(FileNameXVASPfull)) {
        aurostd::DecompressFile(FileNameXVASPfull);
      }
      if (aurostd::FileExist(FileNameXVASP)) { // RECYCLE PREVIOUS STUFF
        stringstream FileXVASPin;
        aurostd::file2stringstream(FileNameXVASP, FileXVASPin);
        FileXVASPout << FileXVASPin.str();
      }
      FileXVASPout << strstream.str();
      aurostd::stringstream2file(FileXVASPout, xvasp.Directory + "/" + DEFAULT_AFLOW_QMVASP_OUT);
    }
    xvasp.str.qm_calculated = true;
    return strstream.str();
  }
} // namespace KBIN

namespace KBIN {
  void GenerateAflowinFromVASPDirectory(_aflags& aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ifstream FileSUBDIR;
    string FileNameSUBDIR;
    FileNameSUBDIR = aflags.Directory;
    FileSUBDIR.open(FileNameSUBDIR.c_str(), std::ios::in);
    FileSUBDIR.clear();
    FileSUBDIR.close();
    ostringstream aus;

    if (aflags.Directory.at(0) != '/' && aflags.Directory.at(0) != '.' && aflags.Directory.at(0) != ' ') {
      aflags.Directory = "./" + aflags.Directory;
    }

    if (!FileSUBDIR) { // ******* Directory is non existent
      aus << "Directory not found";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _FILE_NOT_FOUND_);
    } else { // ******* Directory EXISTS
      // Check LOCK again
      // ifstream FileLOCK0;string FileNameLOCK0=aflags.Directory+"/"+_AFLOWLOCK_;    FileLOCK0.open(FileNameLOCK0.c_str(),std::ios::in);FileLOCK0.close();
      // ifstream FileLOCK1;string FileNameLOCK1=aflags.Directory+"/"+_AFLOWLOCK_+".gz"; FileLOCK1.open(FileNameLOCK1.c_str(),std::ios::in);FileLOCK1.close();
      // ifstream FileLOCK2;string FileNameLOCK2=aflags.Directory+"/"+_AFLOWLOCK_+".bz2";FileLOCK2.open(FileNameLOCK2.c_str(),std::ios::in);FileLOCK2.close();
      // ifstream FileSKIP0;string FileNameSKIP0=aflags.Directory+"/SKIP";    FileSKIP0.open(FileNameSKIP0.c_str(),std::ios::in);FileSKIP0.close();
      // ifstream FileSKIP1;string FileNameSKIP1=aflags.Directory+"/SKIP.gz"; FileSKIP1.open(FileNameSKIP1.c_str(),std::ios::in);FileSKIP1.close();
      // ifstream FileSKIP2;string FileNameSKIP2=aflags.Directory+"/SKIP.bz2";FileSKIP2.open(FileNameSKIP2.c_str(),std::ios::in);FileSKIP2.close();
      // // CHECK FOR LOCK
      // if(FileLOCK0 || FileLOCK1 || FileLOCK2) {                                                                 // ******* Directory is locked
      // 	// LOCK exist, then RUN already RUN
      // 	aus << "EEEEE  LOCKED" << Message(__AFLOW_FILE__,aflags) << endl;
      // 	aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // }
      // if(FileSKIP0 || FileSKIP1 || FileSKIP2) {                                                                 // ******* Directory is skipped
      // 	// SKIP exist, then RUN already RUN
      // 	aus << "EEEEE  SKIPPED" << Message(__AFLOW_FILE__,aflags) << endl;
      // 	aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // }
      if (aurostd::FileExist(aflags.Directory + "/" + _AFLOWLOCK_) || aurostd::CompressFileExist(aflags.Directory + "/" + _AFLOWLOCK_)) { // ******* Directory is locked
        // LOCK exist, then RUN already RUN
        aus << "Directory LOCKED";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
      }
      if (aurostd::FileExist(aflags.Directory + "/SKIP") || aurostd::CompressFileExist(aflags.Directory + "/SKIP")) { // ******* Directory is skipped
        // SKIP exist, then RUN already RUN
        aus << "Directory SKIPPED";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
      }

      // ******* Directory is un locked/skipped
      /// ******************************************************************
      // RESET LOCK
      ofstream FileLOCK;
      const string FileNameLOCK = aflags.Directory + "/" + _AFLOWLOCK_;
      FileLOCK.open(FileNameLOCK.c_str(), std::ios::out);
      /// ******************************************************************
      // CHECK FOR INCAR KPOINTS POSCAR POTCAR
      ifstream FileINCAR;
      const string FileNameINCAR = aflags.Directory + "/INCAR";
      FileINCAR.open(FileNameINCAR.c_str(), std::ios::in);
      if (!FileINCAR) { // ******* INCAR does not exist
        aus << "EEEEE  INCAR ABSENT  =" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK, aus, XHOST.QUIET);
      }
      ifstream FileKPOINTS;
      const string FileNameKPOINTS = aflags.Directory + "/KPOINTS";
      FileKPOINTS.open(FileNameKPOINTS.c_str(), std::ios::in);
      if (!FileKPOINTS) { // ******* KPOINTS does not exist
        aus << "EEEEE  KPOINTS ABSENT  =" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK, aus, XHOST.QUIET);
      }
      ifstream FilePOSCAR;
      const string FileNamePOSCAR = aflags.Directory + "/POSCAR";
      FilePOSCAR.open(FileNamePOSCAR.c_str(), std::ios::in);
      if (!FilePOSCAR) { // ******* POSCAR does not exist
        aus << "EEEEE  POSCAR ABSENT  =" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK, aus, XHOST.QUIET);
      }
      ifstream FilePOTCAR;
      const string FileNamePOTCAR = aflags.Directory + "/POTCAR";
      FilePOTCAR.open(FileNamePOTCAR.c_str(), std::ios::in);
      if (!FilePOTCAR) { // ******* POTCAR does not exist
        aus << "EEEEE  POTCAR ABSENT  =" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK, aus, XHOST.QUIET);
      }
      // ----------------------------------------------------------------------------------------------------
      if (FileINCAR && FileKPOINTS && FilePOSCAR && FilePOTCAR) {
        // VASP INCAR KPOINTS POSCAR POTCAR ARE PRESENT
        /// ******************************************************************
        // WRITE LOCK
        aus << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow" << Message(__AFLOW_FILE__, aflags) << endl;
        aus << "MMMMM  (C) " << XHOST.Copyright_Years << ", Stefano Curtarolo - Duke University " << Message(__AFLOW_FILE__, aflags) << endl;
        aus << "MMMMM  High-Throughput ab-initio Computing" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileLOCK, aus, XHOST.QUIET);
        aus << "00000  MESSAGE GENERATING " << _AFLOWIN_ << " from VASP-xCARs files" << Message(__AFLOW_FILE__, aflags) << endl;
        aurostd::PrintMessageStream(FileLOCK, aus, XHOST.QUIET);
        /// ******************************************************************
        // RESET AFLOWIN
        ofstream FileAFLOWIN;
        const string FileNameAFLOWIN = aflags.Directory + "/" + _AFLOWIN_;
        FileAFLOWIN.open(FileNameAFLOWIN.c_str(), std::ios::out);
        /// ******************************************************************
        // WRITE AFLOWIN
        // WRITE TITLE
        string str1;
        string str2;
        getline(FileINCAR, str1);
        FileINCAR.seekg(0);
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[AFLOW] Automatically generated from XCARS by aflow/aflowd " << string(AFLOW_VERSION) << endl;
        FileAFLOWIN << "[AFLOW] Automatic-Flow" << Message(__AFLOW_FILE__, aflags) << endl;
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        // WRITE HEADER
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[AFLOW] " << str1 << endl;
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[AFLOW] input file for aflow " << endl;
        FileAFLOWIN << "[AFLOW] comments with label " << endl;
        FileAFLOWIN << "[AFLOW] separating with __ the options makes them ignored " << endl;
        FileAFLOWIN << "[AFLOW_MODE=VASP] " << endl;
        FileAFLOWIN << "[VASP] *************************************************** " << endl;
        for (int i = 0; i < (int) XHOST.argv.size() - 1; i++) {
          str1 = XHOST.argv[i];
          str2 = XHOST.argv.at(i + 1);
          if (str1 == "--set" && str2.at(0) == '[') {
            aus << "00000  MESSAGE Adding " << str2 << " to " << _AFLOWIN_ << Message(__AFLOW_FILE__, aflags) << endl;
            aurostd::PrintMessageStream(FileLOCK, aus, XHOST.QUIET);
            FileAFLOWIN << str2 << endl;
          }
        }
        // WRITE INCAR
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
        while (getline(FileINCAR, str1)) {
          FileAFLOWIN << "[VASP_INCAR_FILE]" << str1 << endl;
        }
        // WRITE KPOINTS
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
        while (getline(FileKPOINTS, str1)) {
          FileAFLOWIN << "[VASP_KPOINTS_FILE]" << str1 << endl;
        }
        // WRITE POSCAR
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
        while (getline(FilePOSCAR, str1)) {
          FileAFLOWIN << "[VASP_POSCAR_FILE]" << str1 << endl;
        }
        // WRITE POTCAR
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_POTCAR_MODE_EXPLICIT]" << endl;
        while (getline(FilePOTCAR, str1)) {
          FileAFLOWIN << str1 << endl;
        }
        // close everything.
        FileINCAR.clear();
        FileINCAR.close();
        FileKPOINTS.clear();
        FileKPOINTS.close();
        FilePOSCAR.clear();
        FilePOSCAR.close();
        FilePOTCAR.clear();
        FilePOTCAR.close();
        FileAFLOWIN.flush();
        FileAFLOWIN.clear();
        FileAFLOWIN.close();
        /// ******************************************************************
        // everything is done. check if we nned to delete VASP FILES
        if (aurostd::args2flag(XHOST.argv, "--delete_xcars")) {
          aus << "00000  MESSAGE Removing vasp files in" << Message(__AFLOW_FILE__, aflags) << endl;
          aurostd::PrintMessageStream(FileLOCK, aus, XHOST.QUIET);
          aus << "cd " << aflags.Directory << endl;
          aus << "rm -f `ls | grep -v " << _AFLOWIN_ << " | grep -v LOCK ` " << endl;
          aurostd::execute(aus);
        }
      }
      FileLOCK.flush();
      FileLOCK.clear();
      FileLOCK.close();
    }
  }
} // namespace KBIN

namespace KBIN {
  /// @brief Runs a backup of VASP files for an iteration of aflow, e.g. OUTCAR -> OUTCAR.relax1
  /// @param xvasp the xvasp settings object
  /// @param qmwrite
  /// @param ext the extension to save to e.g. 'relax1'
  /// @authors
  /// @mod{ST,20250618,cleaned and added LOCPOT}
  void VASP_Backup(_xvasp& xvasp, bool qmwrite, const string& ext) {
    xvasp.AnalyzeLabel = ext;
    KBIN::VASP_Analyze(xvasp, qmwrite);
    const ostringstream aus;
    const std::filesystem::path dir = xvasp.Directory;

    for (const string& cext : aurostd::compression_suffix) {
      aurostd::RemoveFile(dir / ("core" + cext));
    }

    vector<string> files_to_remove{
        "aflow.qsub.run",
        "aflow.qsub.out",
    };
    vector<string> files_to_save{
        "AECCAR0", "AECCAR1", "AECCAR2", "CHG",   "CHGCAR", "CONTCAR", "DYNMAT", "DOSCAR",  "ELFCAR",      "EIGENVAL",  "IBZKPT",     "LOCPOT",      "INCAR",
        "KPOINTS", "OSZICAR", "OUTCAR",  "PCDAT", "POSCAR", "POTCAR",  "PROCAR", "XDATCAR", "vasprun.xml", "vaspin.h5", "vaspout.h5", "vaspwave.h5", DEFAULT_VASP_OUT,
    };

    if (xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED")) {
      files_to_save.emplace_back("WAVECAR");
    } else {
      files_to_remove.emplace_back("WAVECAR");
    }
    if (xvasp.aopts.flag("FLAG::WAVEDER_PRESERVED")) {
      files_to_save.emplace_back("WAVEDER");
    } else {
      files_to_remove.emplace_back("WAVEDER");
    }

    for (const string& file : files_to_remove) {
      aurostd::RemoveFile(dir / file);
    }
    for (const string& file : files_to_save) {
      aurostd::file2file(dir / file, dir / (file + "." + ext));
    }
  }
} // namespace KBIN

namespace KBIN {
  void VASP_CONTCAR_Save(const _xvasp& xvasp, const string& ext) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    return VASP_CONTCAR_Save(xvasp.Directory, ext); // CO20210716
  }
  void VASP_CONTCAR_Save(const string& directory, const string& ext) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    const string function = "KBIN::VASP_CONTCAR_Save";
    const string operation = function + "(" + ext + ")";
    const string aflowin = directory + "/" + _AFLOWIN_;
    const string contcar = directory + string("/CONTCAR");
    if (aurostd::FileExist(aflowin) && aurostd::FileExist(contcar) && !aurostd::FileEmpty(contcar)) {
      const xstructure xstr(contcar, IOAFLOW_AUTO);
      ostringstream aus;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      aus << "[AFLOW] SELF-MODIFICATION" << endl;
      aus << "[AFLOW] Recycling CONTCAR of " << ext << endl;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      aus << _VASP_POSCAR_MODE_EXPLICIT_START_ << endl;
      aus << xstr;
      aus << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << endl;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      KBIN::AFLOWIN_ADD(aflowin, aus, "");
      KBIN::AFLOWIN_REMOVE(aflowin, "[VASP_FORCE_OPTION]VOLUME", operation); // CO20210315
    }
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Recycle(const _xvasp& xvasp, const string& ext) { // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210315
    aurostd::CopyFile(xvasp.Directory + "/CONTCAR." + ext, xvasp.Directory + "/POSCAR");
    aurostd::CopyFile(xvasp.Directory + "/INCAR." + ext, xvasp.Directory + "/INCAR");
    aurostd::CopyFile(xvasp.Directory + "/KPOINTS." + ext, xvasp.Directory + "/KPOINTS");
    aurostd::CopyFile(xvasp.Directory + "/POTCAR." + ext, xvasp.Directory + "/POTCAR");
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Recycle(const _xvasp& xvasp, int relax_number) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      aurostd::execute(XHOST.vzip[iext] + " -dqf " + aurostd::CleanFileName(xvasp.Directory + "/*" + XHOST.vext[iext]));
    }
    aurostd::CopyFile(xvasp.Directory + "/CONTCAR.relax" + aurostd::utype2string<int>(relax_number), xvasp.Directory + "/POSCAR");
    aurostd::CopyFile(xvasp.Directory + "/INCAR.relax" + aurostd::utype2string<int>(relax_number), xvasp.Directory + "/INCAR");
    aurostd::CopyFile(xvasp.Directory + "/KPOINTS.relax" + aurostd::utype2string<int>(relax_number), xvasp.Directory + "/KPOINTS");
    aurostd::CopyFile(xvasp.Directory + "/POTCAR.relax" + aurostd::utype2string<int>(relax_number), xvasp.Directory + "/POTCAR");
  }
} // namespace KBIN

namespace KBIN {
  void VASP_RecycleExtraFile(const _xvasp& xvasp, const string& xfile, const string& relax) { // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    aurostd::CopyFile(xvasp.Directory + "/" + xfile + "." + relax, xvasp.Directory + "/" + xfile);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_RecycleExtraFile(const _xvasp& xvasp, const string& xfile, int relax_number) { // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210315
    for (size_t iext = 1; iext < XHOST.vext.size(); iext++) { // SKIP uncompressed
      aurostd::execute(XHOST.vzip[iext] + " -dqf " + aurostd::CleanFileName(xvasp.Directory + "/" + xfile + XHOST.vext[iext]));
    }
    aurostd::CopyFile(xvasp.Directory + "/" + xfile + ".relax" + aurostd::utype2string<int>(relax_number), xvasp.Directory + "/" + xfile);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_BackupOriginal(_xvasp xvasp) { // AFLOW_FUNCTION_IMPLEMENTATION
    aurostd::CopyFile(xvasp.Directory + "/KPOINTS", xvasp.Directory + "/KPOINTS.orig");
    aurostd::CopyFile(xvasp.Directory + "/INCAR", xvasp.Directory + "/INCAR.orig");
    aurostd::CopyFile(xvasp.Directory + "/POSCAR", xvasp.Directory + "/POSCAR.orig");
  }
} // namespace KBIN

namespace KBIN {
  // SD20221026 - eliminated subshell
  // CO20210315 - reconsider rewriting these functions to eliminate subshell
  // NELM comes from xOUTCAR
  // NSTEPS comes from xOSZICAR (to be created)
  uint VASP_getNELM(const string& outcar) { // CO20200624
    const bool LDEBUG = (false || _DEBUG_KVASP_ || XHOST.DEBUG);
    ifstream FileOUTCAR;
    FileOUTCAR.open(outcar.c_str(), std::ios::in);
    const string tmp = aurostd::kvpair2string(FileOUTCAR, "NELM", "=");
    FileOUTCAR.close();
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " " << outcar << " NELM grep response=\"" << tmp << "\"" << endl;
    }
    int NELM = 60; // VASP default
    if (!tmp.empty() && aurostd::isfloat(tmp)) {
      NELM = aurostd::string2utype<int>(tmp);
    }
    return NELM;
  }
  uint VASP_getNSTEPS(const string& oszicar) { // CO20200624
    const bool LDEBUG = (false || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    ifstream FileOSZICAR;
    FileOSZICAR.open(oszicar.c_str(), std::ios::in);
    string tmp = aurostd::kvpair2string(FileOSZICAR, "DAV", ":", -1);
    FileOSZICAR.close();
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " " << oszicar << " NSTEPS grep response=\"" << tmp << "\"" << endl;
    }
    int NSTEPS = 0; // VASP default
    if (!tmp.empty()) {
      if (aurostd::isfloat(tmp)) {
        NSTEPS = aurostd::string2utype<int>(tmp);
      } else if (tmp.find("*") != string::npos) {
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " found number bigger than 999" << endl;
        }
        vector<string> vlines;
        vector<string> tokens;
        aurostd::compressfile2vectorstring(oszicar, vlines);
        bool found_F_line = false;
        uint i = 0;
        uint j = 0;
        uint nsteps = 0;
        for (i = vlines.size() - 1; i < vlines.size() && NSTEPS == 0; i--) { // go backwards
          if (found_F_line) {
            aurostd::string2tokens(vlines[i], tokens, " ");
            if (tokens.size() < 2) {
              continue;
            }
            if (tokens[0].find(":") == string::npos) {
              continue;
            } // safety check, should be "DAV:" or "RMM:"
            tmp = tokens[1];
            if (aurostd::isfloat(tmp)) {
              nsteps = aurostd::string2utype<uint>(tmp);
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " last countable nstep: " << nsteps << endl;
                cerr << __AFLOW_FUNC__ << " steps after last counterable NSTEP: " << j << endl;
              }
              NSTEPS = nsteps + j;
              if (LDEBUG) {
                cerr << __AFLOW_FUNC__ << " NSTEPS=" << NSTEPS << endl;
              }
            } else {
              j++;
            }
          }
          if (vlines[i].find("F=") != string::npos) {
            found_F_line = true;
            continue;
          }
        }
      }
    }
    return NSTEPS;
  }
  bool VASP_OSZICARUnconverging(const string& dir, uint cutoff) { // CO20210601
    // this function will read the whole OSZICAR looking for electronic sloshing issues
    // if there are $cutoff ionic steps that did not converge electronically, return true
    const bool LDEBUG = (false || _DEBUG_KVASP_ || XHOST.DEBUG);
    vector<string> vlines;
    vector<string> vrelax;
    vector<string> tokens;
    aurostd::file2vectorstring(dir + "/OSZICAR", vlines);
    uint i = 0;
    for (i = 1; i < vlines.size(); i++) { // start at 1 so i-1 isn't a problem
      if (vlines[i].find("F=") != string::npos) {
        vrelax.push_back(vlines[i - 1]);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " nionic=" << vrelax.size() << endl;
    }
    if (vrelax.size() < cutoff) {
      return false; // no problem
    }
    // otherwise check for issues.
    const uint NELM = KBIN::VASP_getNELM(dir + "/OUTCAR");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " NELM=" << NELM << endl;
    }
    uint nsteps = 0;
    uint nissues = 0;
    for (i = 0; i < vrelax.size() && nissues < cutoff; i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vrelax[i=" << i << "]=\"" << vrelax[i] << "\"" << endl;
      }
      aurostd::string2tokens(vrelax[i], tokens, " ");
      if (tokens.size() > 1) {
        nsteps = aurostd::string2utype<uint>(tokens[1]);
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " nsteps=" << nsteps << endl;
        }
        if (nsteps != 0 && nsteps >= NELM) {
          nissues++;
        }
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " nissues=" << nissues << endl;
    }
    if (nissues == cutoff) {
      return true;
    }
    return false;
  }
  bool VASP_OSZICARUnconverged(const string& oszicar, const string& outcar) { // CO20210601
    // this function only looks at the last electronic SC step (different than VASP_OSZICARUnconverging, good for STATIC calcs)
    // if it is unconverged, return true
    const bool LDEBUG = (false || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    const uint NELM = KBIN::VASP_getNELM(outcar);
    const uint NSTEPS = KBIN::VASP_getNSTEPS(oszicar);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " NELM=" << NELM << endl;
      cerr << __AFLOW_FUNC__ << " NSTEPS=" << NSTEPS << endl;
    }
    if (NELM != 0 && NSTEPS != 0 && NSTEPS >= NELM) {
      return true;
    }
    return false;
  }
} // namespace KBIN

// ***************************************************************************
// functions written by CAMILO CALDERON
// 2013: camilo.calderon@duke.edu

// todo:
// Finish the DYNADIEL tag
// OUTCAR file & type as a separate subroutine
// Add more options to the statdiel tag (various dielectric tensor types)

// ***************************************************************************
namespace KBIN {
  void GetStatDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi) { // loop GetStatDiel
    string message;
    string outcarpath;
    const string outcarpath_tmp = aurostd::TmpFileCreate("OUTCARc1.tmp");
    vector<string> outcarlines;
    vector<string> endline;
    vector<string> startline;
    vector<string> vasptoken;
    xmatrix<double> statdiel(3, 3);
    const xmatrix<double> eigenvec(3, 3);
    const double eps = 1.0E-5; // need to define this more rigorously
    const string work_dir = aurostd::getPWD(); // CO20191112

    if (!aurostd::FileExist(outcar)) {
      message = "check filename || file missing";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    } else {
      outcarpath = "/" + outcar;
      outcarpath = work_dir + outcarpath;
      vector<string> outcardata;
      aurostd::string2tokens(outcarpath, outcardata, ".");
      if (outcardata[outcardata.size() - 1] == "bz2" || outcardata[outcardata.size() - 1] == "xz" || outcardata[outcardata.size() - 1] == "gz") {
        aurostd::compressfile2vectorstring(outcarpath, outcarlines);
      } else {
        aurostd::file2vectorstring(outcarpath, outcarlines);
      }
    }
    // check the loaded OUTCAR
    aurostd::string2tokens(outcarlines.at(0), startline, " ");
    aurostd::string2tokens(startline.at(0), vasptoken, ".");
    aurostd::string2tokens(outcarlines.at(outcarlines.size() - 1), endline, " ");
    if (vasptoken.at(0) != "vasp" || endline.at(0) != "Voluntary") { // first and last line check
      message = "OUTCAR file is probably corrupt";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    uint sec_count = 0;
    for (size_t ii = outcarlines.size() - 12; ii < outcarlines.size(); ii++) { // presence timing information check
      vector<string> timetoken;
      aurostd::string2tokens(outcarlines.at(ii), timetoken, " ");
      if (!timetoken.empty()) {
        for (size_t jj = 0; jj < timetoken.size(); jj++) {
          if (timetoken[jj] == "(sec):") {
            sec_count += 1;
          }
        }
      }
    }
    if (sec_count != 4) { // first and last line check
      message = "OUTCAR file is probably corrupt";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    // OUTCAR is now in memory, now parse the info
    vector<string> words_line;
    vector<string> vec1;
    vector<string> vec2;
    vector<string> vec3;
    bool check_digit = false;
    uint refline = 0;
    for (size_t ii = outcarlines.size() - 1; ii > 1; ii--) { // line contents
      for (size_t jj = 0; jj < words_line.size(); jj++) {
        const string search_term = "MACROSCOPIC";
        const string test_word = words_line[jj];
        if (test_word == search_term) { // start of dielectric tensor
          refline = ii + 2;
          check_digit = true;
        }
      }
      if (check_digit) { // put the tensor info into the string vectors
        aurostd::string2tokens(outcarlines.at(refline + 0), vec1, " ");
        aurostd::string2tokens(outcarlines.at(refline + 1), vec2, " ");
        aurostd::string2tokens(outcarlines.at(refline + 2), vec3, " ");
        for (uint jj = 1; jj <= 3; jj++) { // string to double, 3x3 matrix, be careful with array bounds
          statdiel(1, jj) = atof(vec1.at(jj - 1).c_str());
          statdiel(2, jj) = atof(vec2.at(jj - 1).c_str());
          statdiel(3, jj) = atof(vec3.at(jj - 1).c_str());
        }
        break;
      }
    }
    if (!check_digit) {
      message = outcar + " lacks MACROSCOPIC statement";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    } // DONE PARSING //
    bool matcheck = false;
    for (uint ii = 1; ii <= 3; ii++) { // clean up spuriously small values: e.g. "-0.000001"
      if (std::abs(statdiel(1, ii)) < eps) {
        statdiel(1, ii) = 0.00;
      }
      if (std::abs(statdiel(2, ii)) < eps) {
        statdiel(2, ii) = 0.00;
      }
      if (std::abs(statdiel(3, ii)) < eps) {
        statdiel(3, ii) = 0.00;
      }
    }
    for (uint ii = 1; ii <= 3; ii++) { // check if it is asymmetric & if large off-diags exist
      for (uint jj = 1; jj <= 3; jj++) {
        const double testdiff = statdiel[ii][jj] - statdiel[jj][ii];
        if (testdiff >= eps) { // eps is a bit arbitrary right now ..
          // serious issues with VASP calculation here:
          message = "asymmetric dielectric tensor";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        } else { // only if small
          statdiel(ii, jj) = statdiel(jj, ii);
        }
        if (ii != jj) {
          if (std::abs(statdiel(ii, jj)) > 0 || std::abs(statdiel(jj, ii)) > 0) {
            matcheck = true;
            break;
          }
        }
      }
      if (matcheck) {
        break;
      }
    }
    matcheck = true;
    if (matcheck) { // diagonalize the 3x3 matrix
      aurostd::eigen(statdiel, eigr, eigi);
    }
  } // loop GetStatDiel
} // namespace KBIN

// ***************************************************************************

namespace KBIN {
  /// @brief gets a string containing the VASP version from binary
  ///
  /// @param binfile path to the binary file
  ///
  /// @return string containing the VASP version
  ///
  /// @authors
  /// @mod{SD,20250218,updated for any VASP version using `vasp --version` and regex}
  /// @mod{SD,20220923,updated for vasp6.3}
  /// @mod{SD,20220401,created function based on getVASPVersionString}
  string BIN2VASPVersion(const string& binfile) {
    const std::regex regex_vasp("vasp\\.((\\.*?[0-9]{1,4}){3})");
    std::smatch match;
    const string version_info = aurostd::execute2string(binfile + " --version");
    if (!std::regex_search(version_info, match, regex_vasp)) {
      return "";
    }
    return match[1];
  }

  /// @brief gets a string containing the VASP version from OUTCAR
  ///
  /// @param outcar path to the OUTCAR
  ///
  /// @return string containing the VASP version
  ///
  /// @authors
  /// @mod{SD,20250218,updated using getline and regex}
  /// @mod{CO,20210315,created function}
  string OUTCAR2VASPVersion(const string& outcar) {
    const std::regex regex_vasp("vasp\\.((\\.*?[0-9]{1,4}){3})");
    std::smatch match;
    std::ifstream file(outcar);
    string line;
    if (file.is_open()) {
      if (!std::getline(file, line)) {
        return "";
      }
      file.close();
    }
    if (!std::regex_search(line, match, regex_vasp)) {
      return "";
    }
    return match[1];
  }
} // namespace KBIN

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
