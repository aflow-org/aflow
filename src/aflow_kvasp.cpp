// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to run VASP in KBIN
// Stefano Curtarolo - 2007-2012 Duke
//   GENERATE, STATIC, KPOINTS, RELAX, RELAX_STATIC, RELAX_STATIC_BANDS, STATIC_BANDS, DIELECTRIC_STATIC, DIELECTRIC_DYNAMIC, DSCF

#ifndef _AFLOW_KVASP_CPP
#define _AFLOW_KVASP_CPP
#define _incarpad_ 26
#include "aflow.h"

#define _KVASP_VASP_SLEEP_   2
#define _KVASP_WAIT_SLEEP_   10
//[CO20201111 created rc parameter SECONDS_SLEEP_VASP_COMPLETION]#define _KVASP_CHECK_SLEEP_  30 //CO20201111  //60   //10
#define KBIN_WRONG_ENTRY_STRING string("WRONG_ENTRY")
#define KBIN_WRONG_ENTRY_NUMBER -123

#define _DEBUG_KVASP_ false  //CO20190116

using aurostd::RemoveWhiteSpaces;
using aurostd::RemoveWhiteSpacesFromTheBack;
using aurostd::FileExist;

pthread_mutex_t mutex_KVASP=PTHREAD_MUTEX_INITIALIZER;

#define _VASP_CONTCAR_SAVE_  TRUE

#define _STROPT_ string("[VASP_FORCE_OPTION]")
#define LDAU_ADIABATIC_RELAX_DEFAULT 6
#define DUKE_MATERIALS_VASP5_CORES_DIELECTRIC 16
#define AFLOWLIB_VASP5_CORES_DIELECTRIC 16
#define DUKE_BETA_VASP5_CORES_DIELECTRIC 16


// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

namespace KBIN {
  _kflags VASP_Get_Kflags_from_AflowIN(const string &AflowIn,_aflags &aflags,ostream& oss) { //CO20200624
    ofstream FileMESSAGE("/dev/null");
    return KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,oss);
  }
} // namespace KBIN


namespace KBIN {
  _kflags VASP_Get_Kflags_from_AflowIN(const string &_AflowIn,ofstream &FileMESSAGE,_aflags &aflags,ostream& oss) {  //CO20200624
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    _kflags kflags;
    string AflowIn=aurostd::RemoveComments(_AflowIn); // for safety //CO20180502
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
    string BflowIn=AflowIn;
    ostringstream aus;

    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (START)" << endl;

    if(aflags.Directory.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"aflags.Directory not set",_INPUT_MISSING_);}  //CO20200624 - prevent segfault
    if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;

    // ***************************************************************************
    // FIND MPI	
    kflags.KBIN_MPI= aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI]");  // search for MPI string
    // ***************************************************************************
    // FIND HOST
    // duke_beta	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]BETA") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_BETA"))   // check DUKE_BETA
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // duke_beta_openmpi	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]BETA_OPENMPI") ||    // check DUKE_BETA_OPENMPI
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_BETA_OPENMPI"))   // check DUKE_BETA_OPENMPI
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // duke_qrats	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]QRATS") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_QRATS"))   // check DUKE_QRATS
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // duke_qflow
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QFLOW_OPENMPI") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]QFLOW") ||  //backwards compatible //CO20180409
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_QFLOW") || //backwards compatible //CO20180409
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]QUSER") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_QUSER"))   // check DUKE_QFLOW
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //CO20201220 X START
    // duke_x
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_X") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]X_X") ||  //backwards compatible //CO20180409
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_X_X"))  //check DUKE_X_X //CO20180409
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_X")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_CRAY") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]X_CRAY") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_X_CRAY"))  //check DUKE_X_CRAY //SD20221006
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_CRAY")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_OLDCRAY") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]X_OLDCRAY") ||  
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_X_OLDCRAY"))  //check DUKE_X_OLDCRAY //SD20221006
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_OLDCRAY")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_X_SMB") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]X_SMB") ||  
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_X_SMB"))  //check DUKE_X_SMB //SD20221006
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_SMB")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //CO20201220 X STOP
    //CO20220818 JHU_ROCKFISH START
    // jhu_rockfish
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::JHU_ROCKFISH") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]ROCKFISH") ||  //backwards compatible //CO20220830
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]JHU_ROCKFISH"))  //check JHU_ROCKFISH //CO20220830
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::JHU_ROCKFISH")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //CO20220818 JHU_ROCKFISH STOP
    // mpcdf_eos	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]EOS") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_EOS"))   // check MPCDF_EOS
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // mpcdf_draco	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DRACO") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_DRACO"))   // check MPCDF_DRACO
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // mpcdf_cobra	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_COBRA") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]COBRA") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_COBRA"))   // check MPCDF_COBRA
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // mpcdf_hydra	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]HYDRA") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MPCDF_HYDRA"))   // check MPCDF_HYDRA
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20190509 - MACHINE001 - START
    // machine001	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE001") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE001"))   // check MACHINE001
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20190509 - MACHINE001 - END
    //DX20190509 - MACHINE002 - START
    // machine002
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE002") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE002"))   // check MACHINE002
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20190509 - MACHINE002 - END
    //DX20201005 - MACHINE003 - START
    // machine003
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE003") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE003") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE003"))   // check MACHINE003
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE003")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20201005 - MACHINE003 - END
    //DX20211011 - MACHINE004 - START
    // machine004
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE004") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE004") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE004"))   // check MACHINE004
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE004")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20211011 - MACHINE004 - END
    // duke_materials	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_MATERIALS") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MATERIALS") ||    // check DUKE_MATERIALS
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_MATERIALS"))   // check DUKE_MATERIALS
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // duke_aflowlib	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_AFLOWLIB") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]AFLOWLIB") ||    // check DUKE_AFLOWLIB
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_AFLOWLIB"))   // check DUKE_AFLOWLIB
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // duke_habana	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_HABANA") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]HABANA") ||    // check DUKE_HABANA
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]DUKE_HABANA"))   // check DUKE_HABANA
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // fulton_marylou	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::FULTON_MARYLOU") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MARYLOU") ||    // check FULTON_MARYLOU
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]FULTON_MARYLOU"))   // check FULTON_MARYLOU
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //OL	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::OHAD") || //CO20181113
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE2") ||  // check MACHINE2
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE2"))   // check MACHINE2
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) { //CO20181113
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    // host1	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::HOST1") || //CO20181113
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE1") ||  // check MACHINE1
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]MACHINE1"))   // check MACHINE1
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) { //CO20181113
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20190107 - CMU EULER - START
    // cmu_euler	
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::CMU_EULER") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]CMU_EULER") || 
        aurostd::substring2bool(AflowIn,"[AFLOW_HOST]CMU_EULER"))   // check CMU_EULER
      aflags.AFLOW_MACHINE_LOCAL=aflags.AFLOW_MACHINE_GLOBAL;
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) {
      aus << "00000  MESSAGE Taking HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
    }
    //DX20190107 - CMU EULER - END

    // ***************************************************************************
    // OTHER CHECKS FOR MPI
    // machines are done withing the VASP/ALIEN stuff, if necessary
    if(aflags.AFLOW_FORCE_MPI) kflags.KBIN_MPI=TRUE;      // forcing
    if(aflags.AFLOW_FORCE_SERIAL) kflags.KBIN_MPI=FALSE;  // forcing

    kflags.KBIN_QSUB= aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]") && !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE");  // search for QSUB string
    kflags.KBIN_QSUB_MODE1=aflags.AFLOW_MODE_QSUB_MODE1 || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE1"); // search for QSUB string mode1
    kflags.KBIN_QSUB_MODE2=aflags.AFLOW_MODE_QSUB_MODE2 || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE2"); // search for QSUB string mode2
    kflags.KBIN_QSUB_MODE3=aflags.AFLOW_MODE_QSUB_MODE3 || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_QSUB]MODE3"); // search for QSUB string mode3
    kflags.AFLOW_MODE_ALIEN=                                               // check ALIEN
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE=ALIEN]") ||             // check ALIEN
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ALIEN]") ||             // check ALIEN
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE]ALIEN");                // check ALIEN
    kflags.AFLOW_MODE_MATLAB=                                              // check MATLAB
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE=MATLAB]") ||            // check MATLAB
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MATLAB]") ||            // check MATLAB
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE]MATLAB");               // check MATLAB
    kflags.AFLOW_MODE_VASP=                                                // check VASP
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE=VASP]") ||              // check VASP
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE_VASP]") ||              // check VASP
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE]VASP");                 // check VASP
    kflags.AFLOW_MODE_AIMS=                                                // check AIMS
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE=AIMS]") ||              // check AIMS
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE_AIMS]") ||              // check AIMS
      aurostd::substring2bool(AflowIn,"[AFLOW_MODE]AIMS");                 // check AIMS
    //CO20180406 - fix generate flags
    if(aflags.KBIN_GEN_GENERAL){
      if(kflags.AFLOW_MODE_AIMS && !aflags.KBIN_GEN_AIMS_FROM_AFLOWIN){aflags.KBIN_GEN_AIMS_FROM_AFLOWIN=true;} //very safe
      if(kflags.AFLOW_MODE_VASP && !aflags.KBIN_GEN_VASP_FROM_AFLOWIN){aflags.KBIN_GEN_VASP_FROM_AFLOWIN=true;} //do vasp last, default
    }
    kflags.KBIN_SYMMETRY_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_SYMMETRY]CALC",TRUE);
    //DX START
    kflags.KBIN_SYMMETRY_NO_SCAN  = aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]NO_SCAN",TRUE);
    //cerr << kflags.KBIN_SYMMETRY_EPS << endl;
    if(aurostd::substring2bool(AflowIn,"[AFLOW_SYMMETRY]SYM_EPS=",TRUE)){
      kflags.KBIN_SYMMETRY_EPS      = aurostd::substring2utype<double>(AflowIn,"[AFLOW_SYMMETRY]SYM_EPS=",TRUE);
    }
    //DX END
    // ---------------------------------------------------------
    // parameters for AAPL - CO20170601
    // to make backwards compatible, we need to not only look for substring, but need to see if "KAPPA=y"
    // start with AAPL first, then QHA, then APL, they are mutually exclusive
    aurostd::xoption KBIN_PHONONS_CALCULATION_AAPL;
    KBIN_PHONONS_CALCULATION_AAPL.option=false;
    KBIN_PHONONS_CALCULATION_AAPL.options2entry(AflowIn, string("[AFLOW_AAPL]KAPPA=|[AFLOW_PHONONS]KAPPA="), KBIN_PHONONS_CALCULATION_AAPL.option, KBIN_PHONONS_CALCULATION_AAPL.xscheme); //CO20170601
    KBIN_PHONONS_CALCULATION_AAPL.option |= aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AAPL]CALC",TRUE);  //legacy
    kflags.KBIN_PHONONS_CALCULATION_AAPL  = KBIN_PHONONS_CALCULATION_AAPL.option;
    // ---------------------------------------------------------
    // parameters for QHA - CO20170601
    // to make backwards compatible, we need to not only look for substring, but need to see if "[AFLOW_QHA]CALC"
    if(!kflags.KBIN_PHONONS_CALCULATION_AAPL){  //mutually exclusive
      kflags.KBIN_PHONONS_CALCULATION_QHA  = aurostd::substring2bool(AflowIn,"[AFLOW_QHA]CALC",TRUE) || aurostd::substring2bool(AflowIn,"VASP_QHA]CALC",TRUE);
      /////////////////////////////
      //aurostd::xoption KBIN_PHONONS_CALCULATION_QHA; //PN20180705
      //KBIN_PHONONS_CALCULATION_QHA.option=false; //PN20180705
      //KBIN_PHONONS_CALCULATION_QHA.options2entry(AflowIn, string("[AFLOW_QHA]GRUNEISEN=|[AFLOW_PHONONS]GRUNEISEN="), KBIN_PHONONS_CALCULATION_QHA.option, KBIN_PHONONS_CALCULATION_QHA.xscheme); //CO20170601 //PN20180705
      //KBIN_PHONONS_CALCULATION_QHA.option |= aurostd::substring2bool(AflowIn,"[AFLOW_QHA]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_QHA]CALC",TRUE); //legacy //PN20180705
      //kflags.KBIN_PHONONS_CALCULATION_QHA  = KBIN_PHONONS_CALCULATION_QHA.option; //PN20180705
    }
    // ---------------------------------------------------------
    // parameters for APL
    // if(LDEBUG) cout << XPID << "KBIN::RUN_Directory: kflags.KBIN_PHONONS_CALCULATION_APL=" << kflags.KBIN_PHONONS_CALCULATION_APL << endl;
    if(!(kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA)){ //mutually exclusive
      kflags.KBIN_PHONONS_CALCULATION_APL  = aurostd::substring2bool(AflowIn,"[AFLOW_APL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_PHONONS]CALC",TRUE);
    }
    // if(LDEBUG) cout << XPID << "KBIN::RUN_Directory: kflags.KBIN_PHONONS_CALCULATION_APL=" << kflags.KBIN_PHONONS_CALCULATION_APL << endl;
    // ---------------------------------------------------------
    // parameters for AGL (Debye Model)
    //Cormac created CALCSTRAINORIGIN, so we need to check [AFLOW_AEL]CALC vs. [AFLOW_AEL]CALCSTRAINORIGIN
    //kflags.KBIN_PHONONS_CALCULATION_AGL  = aurostd::substring2bool(AflowIn,"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AGL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_GIBBS]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_GIBBS]CALC",TRUE);
    for(uint i=0;i<vAflowIn.size()&&!kflags.KBIN_PHONONS_CALCULATION_AGL;i++){
      if((aurostd::substring2bool(vAflowIn[i],"[AFLOW_AGL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AGL]CALC",TRUE)) 
          && !(
            aurostd::substring2bool(vAflowIn[i],"[AFLOW_AGL]CALC_",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AGL]CALC_",TRUE) ||
            aurostd::substring2bool(vAflowIn[i],"[AFLOW_AGL]CALCS",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AGL]CALCS",TRUE) ||
            FALSE)){
        kflags.KBIN_PHONONS_CALCULATION_AGL=true;
      }
    }
    // ---------------------------------------------------------
    // parameters for AEL (Elastic constants)
    //Cormac created CALCSTRAINORIGIN, so we need to check [AFLOW_AEL]CALC vs. [AFLOW_AEL]CALCSTRAINORIGIN
    //kflags.KBIN_PHONONS_CALCULATION_AEL  = aurostd::substring2bool(AflowIn,"[AFLOW_AEL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AEL]CALC",TRUE);
    for(uint i=0;i<vAflowIn.size()&&!kflags.KBIN_PHONONS_CALCULATION_AEL;i++){
      if((aurostd::substring2bool(vAflowIn[i],"[AFLOW_AEL]CALC",TRUE) || aurostd::substring2bool(AflowIn,"[VASP_AEL]CALC",TRUE)) 
          && !(
            aurostd::substring2bool(vAflowIn[i],"[AFLOW_AEL]CALC_",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AEL]CALC_",TRUE) ||
            aurostd::substring2bool(vAflowIn[i],"[AFLOW_AEL]CALCS",TRUE) || aurostd::substring2bool(vAflowIn[i],"[VASP_AEL]CALCS",TRUE) ||
            FALSE)){
        kflags.KBIN_PHONONS_CALCULATION_AEL=true;
      }
    }
    // ---------------------------------------------------------
    // Warn user if both APL/AAPL and AEL/AGL flags are set, since they are mutually exclusive
    if ((kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_QHA) && (kflags.KBIN_PHONONS_CALCULATION_AGL || kflags.KBIN_PHONONS_CALCULATION_AEL)) {
      aus << "WWWWW  WARNING: APL/AAPL/QHA and AEL/AGL flags both set" << endl;
      aus << "WWWWW  WARNING: These runs are mutually exclusive" << endl;
      aus << "WWWWW  WARNING: APL/AAPL/QHA runs will take priority" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    } //CT20200520 - added warning
    // ---------------------------------------------------------
    // parameters for NEIGHBORS
    //DX20210122 [OBSOLETE] kflags.KBIN_NEIGHBORS_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_NEIGHBOURS]CALC",TRUE) || aurostd::substring2bool(AflowIn, "[AFLOW_NEIGHBORS]CALC");  //ME20190107 - Added American spelling
    // ---------------------------------------------------------
    // parameters for POCC CALCULATIONS, KESONG YANG
    kflags.KBIN_POCC=FALSE;
    kflags.KBIN_POCC_CALCULATION  = aurostd::substring2bool(AflowIn,"[AFLOW_POCC]CALC",TRUE) && (aurostd::substring2bool(AflowIn,"[POCC_MODE_EXPLICIT]START.POCC_STRUCTURE",TRUE) && aurostd::substring2bool(AflowIn,"[POCC_MODE_EXPLICIT]STOP.POCC_STRUCTURE",TRUE)); //CO20180419
    if(kflags.KBIN_POCC_CALCULATION) {
      aus << "00000  MESSAGE POCC_CALCULATION "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    }
    if(kflags.KBIN_POCC_CALCULATION) {kflags.KBIN_POCC=TRUE;} //CO20180419
    if(kflags.KBIN_POCC){ //CO20191110
      kflags.KBIN_POCC_TEMPERATURE_STRING=aurostd::substring2string(AflowIn,"[AFLOW_POCC]TEMPERATURE=");  //CO20191110
      if(kflags.KBIN_POCC_TEMPERATURE_STRING.empty()){  //CO20191110
        kflags.KBIN_POCC_TEMPERATURE_STRING=DEFAULT_POCC_TEMPERATURE_STRING;  //CO20191110
      }
      aus << "00000  MESSAGE POCC_TEMPERATURE_STRING=" << kflags.KBIN_POCC_TEMPERATURE_STRING << Message(__AFLOW_FILE__,aflags) << endl;  //CO20191110
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);  //CO20191110
      kflags.KBIN_POCC_ARUNS2SKIP_STRING=aurostd::substring2string(AflowIn,"[AFLOW_POCC]ARUNS2SKIP=");  //CO20200624
      if(!kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()){  //CO20200624
        aus << "00000  MESSAGE POCC_ARUNS2SKIP_STRING=" << kflags.KBIN_POCC_ARUNS2SKIP_STRING << Message(__AFLOW_FILE__,aflags) << endl;  //CO20200624
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);  //CO20200624
      }
      //ME20210927 - EXCLUDE_UNSTABLE
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
        aus << "00000  MESSAGE POCC_EXCLUDE_UNSTABLE=" << (kflags.KBIN_POCC_EXCLUDE_UNSTABLE?"TRUE":"FALSE") << Message(__AFLOW_FILE__,aflags) << endl;
      }
    }
    // ---------------------------------------------------------
    // parameters for FROZSL
    kflags.KBIN_FROZSL=FALSE;
    kflags.KBIN_PHONONS_CALCULATION_FROZSL  = aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]CALC",TRUE);
    if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
      aus << "00000  MESSAGE FROZSL_CALCULATION "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    }
    kflags.KBIN_FROZSL_DOWNLOAD     = (aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]DOWN",TRUE) ||
        aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]DOWNLOAD",TRUE));
    if(kflags.KBIN_FROZSL_DOWNLOAD) {
      aus << "00000  MESSAGE FROZSL_DOWNLOAD "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    }
    if(kflags.KBIN_FROZSL_FILE) {
      aus << "00000  MESSAGE FROZSL_FILE "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    }
    kflags.KBIN_FROZSL_FILE  = aurostd::substring2bool(AflowIn,"[AFLOW_FROZSL]FILE",TRUE); // load of file somewhere else
    if(kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_FROZSL_DOWNLOAD|| kflags.KBIN_FROZSL_FILE) kflags.KBIN_FROZSL=TRUE;
    // ---------------------------------------------------------
    // the rest of symmetry stuff is sought inside ivasp or
    if(kflags.AFLOW_MODE_ALIEN) {
      kflags.AFLOW_MODE_MATLAB=FALSE;                  // fix PRIORITY
      kflags.AFLOW_MODE_VASP=FALSE;                    // fix PRIORITY
      kflags.KBIN_MPI=FALSE;                           // fix PRIORITY
    }
    if(kflags.AFLOW_MODE_MATLAB) {
      kflags.AFLOW_MODE_VASP=FALSE;                    // fix PRIORITY
      kflags.KBIN_MPI=FALSE;
    }
    if(LDEBUG) cout << "DEBUG kflags.AFLOW_MODE_ALIEN=" << kflags.AFLOW_MODE_ALIEN << endl;
    if(LDEBUG) cout << "DEBUG kflags.AFLOW_MODE_MATLAB=" << kflags.AFLOW_MODE_MATLAB << endl;
    if(LDEBUG) cout << "DEBUG kflags.AFLOW_MODE_VASP=" << kflags.AFLOW_MODE_VASP << endl;
    // ***************************************************************************
    // ZIP COMPRESS
    // ***************************************************************************
    kflags.KZIP_COMPRESS=TRUE;
    aurostd::StringstreamClean(aus);
    if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=none]") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=NONE]") ||
        !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP")) {
      kflags.KZIP_COMPRESS=FALSE;
      for(int i=0;i<1;i++) {
        aus << "WWWWW  WARNING no compression of output files..." << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintWarningStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    } else {
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP")) { // "[AFLOW_MODE_ZIP=" not found
        kflags.KZIP_BIN=DEFAULT_KZIP_BIN;  // take default
        aus << "00000  MESSAGE Taking DEFAULT KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP]")) { // "[AFLOW_MODE_ZIP]" not found
        kflags.KZIP_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_ZIP]");
        aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_ZIP=")) { // "[AFLOW_MODE_ZIP=" found
        kflags.KZIP_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_ZIP="),']');
        aus << "00000  MESSAGE Taking KZIP_BIN=\"" << kflags.KZIP_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
    }
    // ************************************************************************************************************************************
    // Get the KZIP_BIN name - moved inside EACH mode
    // ************************************************************************************************************************************
    // LOAD PREFIX POSTFIX
    KBIN::StartStopCheck(AflowIn,"[AFLOW_MODE_PRESCRIPT]",kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT,kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP);
    KBIN::StartStopCheck(AflowIn,"[AFLOW_MODE_POSTSCRIPT]",kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT,kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP);
    if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT) {  // [AFLOW_MODE_PRESCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_PRESCRIPT_COMMAND << " file from " << _AFLOWIN_ << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.AFLOW_MODE_PRESCRIPT.str(aurostd::substring2string(AflowIn,"[AFLOW_MODE_PRESCRIPT]",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,kflags.AFLOW_MODE_PRESCRIPT,"[AFLOW_MODE_PRESCRIPT]"); //CO20200624 - FileAFLOWIN->AflowIn
    }
    if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_PRESCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_PRESCRIPT_COMMAND << " file from START/STOP " << _AFLOWIN_ << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.AFLOW_MODE_PRESCRIPT.str(aurostd::substring2string(AflowIn,"[AFLOW_MODE_PRESCRIPT]START","[AFLOW_MODE_PRESCRIPT]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,kflags.AFLOW_MODE_PRESCRIPT,"[AFLOW_MODE_PRESCRIPT]START","[AFLOW_MODE_PRESCRIPT]STOP"); //CO20200624 - FileAFLOWIN->AflowIn
    }
    if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT) {  // [AFLOW_MODE_POSTSCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_POSTSCRIPT_COMMAND << " file from " << _AFLOWIN_ << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.AFLOW_MODE_POSTSCRIPT.str(aurostd::substring2string(AflowIn,"[AFLOW_MODE_POSTSCRIPT]",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,kflags.AFLOW_MODE_POSTSCRIPT,"[AFLOW_MODE_POSTSCRIPT]"); //CO20200624 - FileAFLOWIN->AflowIn
    }
    if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_POSTSCRIPT] construction
      aus << "00000  MESSAGE Generating " << DEFAULT_AFLOW_POSTSCRIPT_COMMAND << " file from START/STOP " << _AFLOWIN_ << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      kflags.AFLOW_MODE_POSTSCRIPT.str(aurostd::substring2string(AflowIn,"[AFLOW_MODE_POSTSCRIPT]START","[AFLOW_MODE_POSTSCRIPT]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn,kflags.AFLOW_MODE_POSTSCRIPT,"[AFLOW_MODE_POSTSCRIPT]START","[AFLOW_MODE_POSTSCRIPT]STOP");  //CO20200624 - FileAFLOWIN->AflowIn
    }
    // ************************************************************************************************************************************
    // ALIEN MODE
    if(kflags.AFLOW_MODE_ALIEN) {
      aus      << XPID << "00000  MESSAGE [AFLOW_MODE=ALIEN] found in " << _AFLOWIN_ << " "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      // ***************************************************************************
      // Getting KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_KBIN_ALIEN_BIN;  // take default
      aurostd::StringstreamClean(aus);
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
        kflags.KBIN_BIN=DEFAULT_KBIN_ALIEN_BIN;  // take default
        aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
        kflags.KBIN_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY]");
        aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
      if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
        kflags.KBIN_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
        aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
      //ME20190107 - Grab the serial binary to propagate into child aflow.in files
      kflags.KBIN_SERIAL_BIN = kflags.KBIN_BIN;
      // ***************************************************************************
      // ALIEN MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL            =
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL") ;
      // ***************************************************************************
    }
    // ************************************************************************************************************************************
    // MATLAB MODE
    if(kflags.AFLOW_MODE_MATLAB) {
      aus      << XPID << "00000  MESSAGE [AFLOW_MODE=MATLAB] found in " << _AFLOWIN_ << " "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      // ***************************************************************************
      // MATLAB MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL            =
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL") ;
      // ***************************************************************************
    }
    // ************************************************************************************************************************************
    // AIMS MODE
    if(kflags.AFLOW_MODE_AIMS) {
      aus      << XPID << "00000  MESSAGE [AFLOW_MODE=AIMS] found in " << _AFLOWIN_ << " "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      aurostd::StringstreamClean(aus);
      if(1){  //no support yet
        // ***************************************************************************
        // Getting KBIN_BIN
        kflags.KBIN_BIN = DEFAULT_AIMS_BIN;  // take default  dont touch MPI as it has already been dealt by  KBIN::MPI_Extract
        if(kflags.KBIN_MPI==FALSE) { // only if no MPI is specified
          if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
            kflags.KBIN_BIN=DEFAULT_AIMS_BIN;  // take default
            aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
          }
          if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
            kflags.KBIN_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY]");
            aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
          }
          if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
            kflags.KBIN_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
            aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
          }
          //ME20190107 - Grab the serial binary to propagate into child aflow.in files
          kflags.KBIN_SERIAL_BIN = kflags.KBIN_BIN;
        } else {
          kflags.KBIN_BIN=kflags.KBIN_MPI_BIN;
          aus << "00000  MESSAGE Taking KBIN_BIN=KBIN_MPI_BIN=\"" << kflags.KBIN_MPI_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
        }
      }
      // ***************************************************************************
      // AIMS MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL            =
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL");
      // ***************************************************************************
    }
    // ************************************************************************************************************************************
    // MPI SWTICHES
    if(kflags.KBIN_MPI) KBIN::MPI_Extract(AflowIn,FileMESSAGE,aflags,kflags);
    // ************************************************************************************************************************************
    // ************************************************************************************************************************************
    // VASP MODE
    if(kflags.AFLOW_MODE_VASP) {
      aus      << XPID << "00000  MESSAGE [AFLOW_MODE=VASP] found in " << _AFLOWIN_ << " "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      // ***************************************************************************
      // Getting KBIN_BIN
      kflags.KBIN_BIN = DEFAULT_VASP_BIN;  // take default  dont touch MPI as it has already been dealt by  KBIN::MPI_Extract
      aurostd::StringstreamClean(aus);
      // old Get BIN
      // 	  if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" not found
      // 	    aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\Message(__AFLOW_FILE__,aflags) << endl;
      // 	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      // 	    //   cerr << "take KBIN=" << kflags.KBIN_BIN << endl;
      // 	  } else {
      // 	    kflags.KBIN_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
      // 	    aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\Message(__AFLOW_FILE__,aflags) << endl;
      // 	    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      // 	  }
      if(kflags.KBIN_MPI==FALSE) { // only if no MPI is specified
        if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY")) { // "[AFLOW_MODE_BINARY=" not found
          kflags.KBIN_BIN=DEFAULT_VASP_BIN;  // take default
          aus << "00000  MESSAGE Taking DEFAULT KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
        }
        if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY]")) { // "[AFLOW_MODE_BINARY]" not found
          kflags.KBIN_BIN=aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY]");
          aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
        }
        if(aurostd::substring2bool(AflowIn,"[AFLOW_MODE_BINARY=")) { // "[AFLOW_MODE_BINARY=" found
          kflags.KBIN_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_BINARY="),']');
          aus << "00000  MESSAGE Taking KBIN_BIN=\"" << kflags.KBIN_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
        }
        //ME20190107 - Grab the serial binary to propagate into child aflow.in files
        kflags.KBIN_SERIAL_BIN = kflags.KBIN_BIN;
      } else {
        kflags.KBIN_BIN=kflags.KBIN_MPI_BIN;
        aus << "00000  MESSAGE Taking KBIN_BIN=KBIN_MPI_BIN=\"" << kflags.KBIN_MPI_BIN << "\" "  << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
      }
      // ***************************************************************************
      // VASP MODE  // must contain EMAIL perform
      kflags.AFLOW_MODE_EMAIL            =
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE_EMAIL]") ||
        aurostd::substring2bool(AflowIn,"[AFLOW_MODE]EMAIL");
    }
    // ***************************************************************************
    // ************************************************************************************************************************************
    // MATLAB MODE
    if(kflags.KBIN_PHONONS_CALCULATION_FROZSL && !kflags.AFLOW_MODE_VASP) {
      aus      << XPID << "00000  MESSAGE [AFLOW_FROZSL]CALC found in " << _AFLOWIN_ << " "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    }
    // ************************************************************************************************************************************
    // NO MODE MODE
    if(!kflags.AFLOW_MODE_VASP && !kflags.AFLOW_MODE_AIMS && !kflags.AFLOW_MODE_MATLAB && !kflags.AFLOW_MODE_ALIEN && !kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
      aus << "EEEEE  [AFLOW_MODE=????] invalid found in     "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [AFLOW_MODE=ALIEN]        is supported "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [AFLOW_MODE=MATLAB]       is supported "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [AFLOW_MODE=VASP]         is supported "  << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "EEEEE  [AFLOW_FROZSL]CALC        is supported "  << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // ***************************************************************************
    // FINALIZE LOCK
    aus << "XXXXX  KBIN DIRECTORY END (aflow" << string(AFLOW_VERSION) << ")  "  << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    // ***************************************************************************
    // PREPARE MESSAGE FOR LOG TO BE INTERCEPTED IN COMPRESSION
    aus << "XXXXX  KBIN_DIRECTORY_END " << aflags.Directory << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET,oss);
    // ***************************************************************************

    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (END)" << endl;

    return kflags;
  }
}

namespace KBIN {
  _vflags VASP_Get_Vflags_from_AflowIN(const string &AflowIn,_aflags &aflags,_kflags &kflags,ostream& oss) {
    ofstream FileMESSAGE("/dev/null");
    return KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags,oss);
  }
} // namespace KBIN


namespace KBIN {
  _vflags VASP_Get_Vflags_from_AflowIN(const string &_AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,ostream& oss) {
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    string message = "";
    _vflags vflags;
    string AflowIn=aurostd::RemoveComments(_AflowIn); // for safety //CO20180502
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
    string BflowIn=AflowIn;

    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (START)" << endl;
    // HOW TO RUN
    vflags.KBIN_VASP_RUN_NRELAX=0;
    // [OBSOLETE]  vflags.KBIN_VASP_RUN_GENERATE           =(aurostd::substring2bool(AflowIn,"[VASP_RUN_GENERATE]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]GENERATE")) || aflags.KBIN_GEN_VASP_FROM_AFLOWIN;
    // [OBSOLETE] vflags.KBIN_VASP_RUN_STATIC              =(aurostd::substring2bool(AflowIn,"[VASP_RUN_STATIC]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]STATIC"));
    // [OBSOLETE] vflags.KBIN_VASP_RUN_KPOINTS             =(aurostd::substring2bool(AflowIn,"[VASP_RUN_KPOINTS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]KPOINTS"));

    vflags.KBIN_VASP_RUN.clear();
    if((aurostd::substring2bool(AflowIn,"[VASP_RUN_GENERATE]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]GENERATE")) || aflags.KBIN_GEN_VASP_FROM_AFLOWIN) 
      vflags.KBIN_VASP_RUN.push("GENERATE");
    if((aurostd::substring2bool(AflowIn,"[VASP_RUN_STATIC]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]STATIC")))
      vflags.KBIN_VASP_RUN.push("STATIC");
    if((aurostd::substring2bool(AflowIn,"[VASP_RUN_KPOINTS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]KPOINTS"))) 
      vflags.KBIN_VASP_RUN.push("KPOINTS");

    for(uint i=0;i<vAflowIn.size();i++) {
      if(aurostd::substring2bool(vAflowIn.at(i),"VASP_RUN")) {
        string vasp_run_string=vAflowIn.at(i);
        if(vasp_run_string.find("#")!=string::npos) vasp_run_string=vasp_run_string.substr(0,vasp_run_string.find("#"));
        if(vasp_run_string.find("//")!=string::npos) vasp_run_string=vasp_run_string.substr(0,vasp_run_string.find("//"));
        if(vasp_run_string.find("!")!=string::npos) vasp_run_string=vasp_run_string.substr(0,vasp_run_string.find("!"));

        //      cout << vasp_run_string << endl;
        vector<string> aflowin_tokens;
        aurostd::string2tokens(vasp_run_string,aflowin_tokens,",");
        if(aflowin_tokens.size()>0) {
          // [OBSOLETE] vflags.KBIN_VASP_RUN_RELAX               =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX="));
          if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX="))) vflags.KBIN_VASP_RUN.push("RELAX");
          if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN_RELAX=");
          if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN]RELAX=");
          // [OBSOLETE] vflags.KBIN_VASP_RUN_RELAX_STATIC        =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC="));
          if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC="))) vflags.KBIN_VASP_RUN.push("RELAX_STATIC");
          if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC=");
          if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC=");
          // [OBSOLETE] vflags.KBIN_VASP_RUN_RELAX_STATIC_BANDS  =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS="));
          if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS="))) vflags.KBIN_VASP_RUN.push("RELAX_STATIC_BANDS");
          if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN_RELAX_STATIC_BANDS=");
          if(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS=")) vflags.KBIN_VASP_RUN_NRELAX=aurostd::substring2utype<int>(aflowin_tokens.at(0),"[VASP_RUN]RELAX_STATIC_BANDS=");
          // [OBSOLETE] vflags.KBIN_VASP_RUN_STATIC_BANDS        =(aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_STATIC_BANDS]") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]STATIC_BANDS"));
          if((aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN_STATIC_BANDS]") || aurostd::substring2bool(aflowin_tokens.at(0),"[VASP_RUN]STATIC_BANDS"))) vflags.KBIN_VASP_RUN.push("STATIC_BANDS");
        }

        for(uint j=1;j<aflowin_tokens.size();j++) {
          // [OBSOLETE] vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC   =(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || aurostd::substring2bool(aflowin_tokens.at(j),"DS") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_STATIC"));
          if((vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC") || aurostd::substring2bool(aflowin_tokens.at(j),"DS") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_STATIC"))) vflags.KBIN_VASP_RUN.push("DIELECTRIC_STATIC");
          // [OBSOLETE] vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC  =(vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC || aurostd::substring2bool(aflowin_tokens.at(j),"DD") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_DYNAMIC"));
          if((vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC") || aurostd::substring2bool(aflowin_tokens.at(j),"DD") || aurostd::substring2bool(aflowin_tokens.at(j),"DIELECTRIC_DYNAMIC"))) vflags.KBIN_VASP_RUN.push("DIELECTRIC_DYNAMIC");
          // [OBSOLETE] vflags.KBIN_VASP_RUN_DSCF                =(vflags.KBIN_VASP_RUN_DSCF || aurostd::substring2bool(aflowin_tokens.at(j),"DSCF"));
          if((vflags.KBIN_VASP_RUN.flag("DSCF") || aurostd::substring2bool(aflowin_tokens.at(j),"DSCF"))) vflags.KBIN_VASP_RUN.push("DSCF");
        }
        if(vflags.KBIN_VASP_RUN.flag("DSCF")) vflags.KBIN_VASP_RUN.push("DIELECTRIC_DYNAMIC");
        if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC")) vflags.KBIN_VASP_RUN.push("DIELECTRIC_STATIC");
      }
    }
    if(vflags.KBIN_VASP_RUN.xscheme!="") vflags.KBIN_VASP_RUN.isentry=TRUE;

    // if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) cout << "vflags.KBIN_VASP_RUN.flag(\"DIELECTRIC_STATIC\")" << endl;
    // if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC")) cout << "vflags.KBIN_VASP_RUN.flag(\"DIELECTRIC_DYNAMIC\")" << endl;
    // if(vflags.KBIN_VASP_RUN.flag("DSCF")) cout << "vflags.KBIN_VASP_RUN.flag(\"DSCF\")" << endl;

    // [OBSOLETE] vflags.KBIN_VASP_REPEAT_STATIC        = aurostd::FileExist(aflags.Directory+string("/REPEAT_STATIC")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_STATIC]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_STATIC]");
    // [OBSOLETE] vflags.KBIN_VASP_REPEAT_STATIC_BANDS = aurostd::FileExist(aflags.Directory+string("/REPEAT_STATIC_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_STATIC_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_STATIC_BANDS]");
    // [OBSOLETE] vflags.KBIN_VASP_REPEAT_BANDS        = aurostd::FileExist(aflags.Directory+string("/REPEAT_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_BANDS]");
    // [OBSOLETE] vflags.KBIN_VASP_REPEAT_DELSOL       = aurostd::FileExist(aflags.Directory+string("/REPEAT_DELSOL")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_DELSOL]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_DELSOL]");

    vflags.KBIN_VASP_REPEAT.clear();
    if(aurostd::FileExist(aflags.Directory+string("/REPEAT_STATIC_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_STATIC_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_STATIC_BANDS")) vflags.KBIN_VASP_REPEAT.push("REPEAT_STATIC_BANDS"); //CO20210315 - fixing typo  //CO20210315 - needs to go BEFORE REPEAT_STATIC
    else if(aurostd::FileExist(aflags.Directory+string("/REPEAT_STATIC")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_STATIC]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_STATIC")) vflags.KBIN_VASP_REPEAT.push("REPEAT_STATIC"); //CO20210315 - needs to go AFTER REPEAT_STATIC_BANDS
    else if(aurostd::FileExist(aflags.Directory+string("/REPEAT_BANDS")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_BANDS]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_BANDS")) vflags.KBIN_VASP_REPEAT.push("REPEAT_BANDS"); //CO20210315 - fixing typo
    else if(aurostd::FileExist(aflags.Directory+string("/REPEAT_DELSOL")) || aurostd::substring2bool(AflowIn,"[VASP_RUN_REPEAT_DELSOL]") || aurostd::substring2bool(AflowIn,"[VASP_RUN]REPEAT_DELSOL")) vflags.KBIN_VASP_REPEAT.push("REPEAT_DELSOL"); //CO20210315 - fixing typo

    if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) cout << "vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC\")" << endl;
    if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) cout << "vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")" << endl;

    // priorities about RUN
    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {  // RELAX_STATIC_BANDS
      //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC_BANDS\")==TRUE" << endl;
      vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
      vflags.KBIN_VASP_RUN.push("RELAX_STATIC_BANDS");
    } else {
      if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {  // RELAX_STATIC
        //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC\")==TRUE" << endl;
        vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
        vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
        vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
        vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
        vflags.KBIN_VASP_RUN.push("RELAX_STATIC");
        vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
      } else {
        if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) {  // STATIC_BANDS
          //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\")==TRUE" << endl;
          vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
          vflags.KBIN_VASP_RUN.push("STATIC_BANDS");
          vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
          vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
          vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
          vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
        } else {                                  
          if(vflags.KBIN_VASP_RUN.flag("STATIC")) {  // STATIC
            //	  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"STATIC\")==TRUE" << endl;
            vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
            vflags.KBIN_VASP_RUN.push("STATIC");
            vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
            vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
          } else {                            
            if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) {  // KPOINTS
              //  cerr << "[DEBUG] vflags.KBIN_VASP_RUN.flag(\"KPOINTS\")==TRUE" << endl;
              vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
              vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
              vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
              vflags.KBIN_VASP_RUN.push("KPOINTS");
              vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
              vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
            } else {
              //   cerr << "[DEBUG] DEFAULT" << endl;
              vflags.KBIN_VASP_RUN.push("RELAX");
              vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
              vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
              vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
              vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
              vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
            }
          }
        }
      }
    }
    // priorities about REPEAT
    if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) { //CO20210315
      vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL",FALSE);  //CO20210315
    }
    else if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
      vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE);
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL",FALSE);  //CO20210315
    }
    else if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
      vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL",FALSE);  //CO20210315
    }
    else if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) {
      vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE);
    }

    if(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL || kflags.KBIN_PHONONS_CALCULATION_FROZSL || kflags.KBIN_PHONONS_CALCULATION_AGL || kflags.KBIN_PHONONS_CALCULATION_AEL) {  //CO20170601
      vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_RUN.flag("KPOINTS",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC",FALSE);
      vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC",FALSE);  //CO20210315
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS",FALSE);
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS",FALSE); 
      vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL",FALSE);
      kflags.KBIN_SYMMETRY_CALCULATION=FALSE;
    }

    // RELAX_MODE AND PRIORITIES  // ENERGY | FORCES | ENERGY_FORCES | FORCES_ENERGY (default ENERGY) "
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.options2entry(AflowIn,_STROPT_+"RELAX_MODE=",FALSE,DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY","ENERGY");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCES","FORCES");vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCE","FORCES");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY_FORCES","ENERGY_FORCES");vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("ENERGY_FORCE","ENERGY_FORCES");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCES_ENERGY","FORCES_ENERGY");vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.scheme2scheme("FORCE_ENERGY","FORCES_ENERGY");
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.isentry && vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // FORCE OPTIONS
    vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.options2entry(AflowIn,_STROPT_+"NOTUNE");

    // FORCE OPTIONS SYSTEM_AUTO
    vflags.KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO.options2entry(AflowIn,_STROPT_+"SYSTEM_AUTO");
    vflags.AFLOW_SYSTEM.options2entry(AflowIn,"[AFLOW]SYSTEM=", false, "");  //ME20181121

    // FORCE OPTIONS STATIC RELAX_ALL RELAX_IONS RELAX CELL_VOLUME 
    // cerr << "vflags.KBIN_VASP_FORCE_OPTION_STATIC= " << vflags.KBIN_VASP_FORCE_OPTION_STATIC << endl;
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
    if(aurostd::substring2bool(AflowIn,_STROPT_+"STATIC",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("STATIC");
    if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_ALL",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX",TRUE))  vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("ALL");
    if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS");
    if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_SHAPE",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_SHAPE",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_SHAPE");
    if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_VOLUME",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_VOLUME");
    if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_VOLUME",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_VOLUME");
    //AS20201123 BEGIN
    if(aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_CELL_SHAPE",TRUE)){
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_SHAPE");
    }
    //AS20201123 END
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.xscheme!="") vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry=TRUE;

    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_STATIC=aurostd::substring2bool(AflowIn,_STROPT_+"STATIC",TRUE);
    // [OBSOLETE] cerr << aflow_aconvasp_main.cpp "vflags.KBIN_VASP_FORCE_OPTION_STATIC= " << vflags.KBIN_VASP_FORCE_OPTION_STATIC << endl;
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_RELAX_ALL = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_ALL",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS",TRUE);
    // [OBSOLETE] vvflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_SHAPE = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_SHAPE",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_SHAPE",TRUE);
    // [OBSOLETE] vvflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_VOLUME = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_VOLUME",TRUE);
    // [OBSOLETE] vvflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS_CELL_VOLUME = aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_CELL_VOLUME",TRUE) || aurostd::substring2bool(AflowIn,_STROPT_+"RELAX_IONS_VOLUME",TRUE);
    // [OBSOLETE] if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS_CELL_VOLUME) vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS=FALSE;
    // [OBSOLETE] if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_ALL && (vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS || vflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_SHAPE ||
    // [OBSOLETE]  vflags.KBIN_VASP_FORCE_OPTION_RELAX_CELL_VOLUME || vflags.KBIN_VASP_FORCE_OPTION_RELAX_IONS_CELL_VOLUME)) vflags.KBIN_VASP_FORCE_OPTION_RELAX_ALL=FALSE;

    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS",FALSE);
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL") && (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC") || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS") || 
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE") || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME") ||
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")))    
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL",FALSE);
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC")) {
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.clear();
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("STATIC");
      vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry=TRUE;}

    // PRECISION AND PRIORITIES // (LOW | MEDIUM | NORMAL | HIGH | ACCURATE), PRESERVED
    vflags.KBIN_VASP_FORCE_OPTION_PREC.options2entry(AflowIn,_STROPT_+"PREC=",FALSE,DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('L',"LOW");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('M',"MEDIUM");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('N',"NORMAL");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('H',"HIGH");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('A',"ACCURATE");
    vflags.KBIN_VASP_FORCE_OPTION_PREC.scheme2scheme('P',"PHONONS"); //JJPR Modification
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry && vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_PREC.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_PREC.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ALGO AND PRIORITIES // (NORMAL | VERYFAST | FAST | ALL | DAMPED), PRESERVED
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.options2entry(AflowIn,_STROPT_+"ALGO=",FALSE,DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('N',"NORMAL");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('V',"VERYFAST");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('F',"FAST");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('A',"ALL");
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme('D',"DAMPED");
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.isentry && vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_ALGO.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_ALGO.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved= vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved || aurostd::substring2bool(AflowIn,_STROPT_+"ALGO_PRESERVED",TRUE); // FIX ALGO_PRESERVED

    // ABMIX AND PRIORITIES // empty | [AUTO | US | PAW | #AMIX,#BMIX[,#AMIX_MAG,#BMIX_MAG]]
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.options2entry(AflowIn,_STROPT_+"ABMIX=",FALSE,DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme=KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('A',"AUTO");
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('U',"US");
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('L',"US"); // LDA
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('G',"US"); // GGA
    vflags.KBIN_VASP_FORCE_OPTION_ABMIX.scheme2scheme('P',"PAW");  // something with PAW..
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry && vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message =  "vflags.KBIN_VASP_FORCE_OPTION_ABMIX.content_string="  + vflags.KBIN_VASP_FORCE_OPTION_ABMIX.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // METAGGA AND PRIORITIES // TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " METAGGA" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_METAGGA.options2entry(AflowIn,_STROPT_+"METAGGA=",FALSE,DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME);
    // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=KBIN_WRONG_ENTRY_STRING;
    // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('T',"TPSS");
    // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('R',"RTPSS");
    // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('S',"SCAN");
    // vflags.KBIN_VASP_FORCE_OPTION_METAGGA.scheme2scheme('N',"NONE");
    if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry && vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_METAGGA.content_string="  +  vflags.KBIN_VASP_FORCE_OPTION_METAGGA.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    } 
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " METAGGA vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " METAGGA vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << endl;

    // IVDW AND PRIORITIES // [number_for_VASP_see_manual_for_IVDW | 0] 
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " IVDW" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_IVDW.options2entry(AflowIn,_STROPT_+"IVDW=",FALSE,DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME);
    if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry && vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_IVDW.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_IVDW.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    } 
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " IVDW vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " IVDW vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << endl;

    // NEGLECT_NOMIX
    vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.options2entry(AflowIn,string(_STROPT_+"NEGLECT_IMMISCIBLE"+"|"+_STROPT_+"NEGLECT_NOMIX"+"|"+_STROPT_+"SKIP_NOMIX"));

    // AUTO_PSEUDOPOTENTIALS and AUTO_PSEUDOPOTENTIALS_TYPE
    vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.options2entry(AflowIn,_STROPT_+"AUTO_PSEUDOPOTENTIALS=",FALSE,DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE);

    // POTIM
    // cerr << "POTIM" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.options2entry(AflowIn,_STROPT_+"POTIM=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_VASP_PREC_POTIM"

    // PSTRESS
    // cerr << "PSTRESS" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.options2entry(AflowIn,_STROPT_+"PSTRESS=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.0"  

    // EDIFFG
    // cerr << "EDIFFG" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.options2entry(AflowIn,_STROPT_+"EDIFFG=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_VASP_PREC_EDIFFG"  

    // NELM //CO20200624
    // cerr << "NELM" << endl;  //CO20200624
    vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.options2entry(AflowIn,_STROPT_+"NELM=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_NELM_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "60" - default  //CO20200624

    // NELM_STATIC //CO20200624
    // cerr << "NELM_STATIC" << endl; //CO20200624
    vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.options2entry(AflowIn,_STROPT_+"NELM_STATIC=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "120" - default  //CO20200624

    // ISMEAR
    // cerr << "ISMEAR" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.options2entry(AflowIn,_STROPT_+"ISMEAR=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1" - default  //CO20181128

    // SIGMA
    // cerr << "SIGMA" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.options2entry(AflowIn,_STROPT_+"SIGMA=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.1" - default  //CO20181128

    // ISMEAR_STATIC  //CO20210315
    // cerr << "ISMEAR_STATIC" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.options2entry(AflowIn,_STROPT_+"ISMEAR_STATIC=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1" - default  //CO20181128

    // SIGMA_STATIC //CO20210315
    // cerr << "SIGMA_STATIC" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.options2entry(AflowIn,_STROPT_+"SIGMA_STATIC=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.1" - default  //CO20181128

    // ISMEAR_BANDS  //CO20210315
    // cerr << "ISMEAR_BANDS" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.options2entry(AflowIn,_STROPT_+"ISMEAR_BANDS="+"|"+_STROPT_+"ISMEAR_STATIC_BANDS=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1" - default  //CO20181128  //CO20210624 - backwards compatibility with ISMEAR_STATIC_BANDS

    // SIGMA_BANDS //CO20210315
    // cerr << "SIGMA_BANDS" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.options2entry(AflowIn,_STROPT_+"SIGMA_BANDS="+"|"+_STROPT_+"SIGMA_STATIC_BANDS=",FALSE,vflags.KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.1" - default  //CO20181128  //CO20210624 - backwards compatibility with SIGMA_STATIC_BANDS

    // NBANDS and/or NBANDS=
    //  cerr << "NBANDS_AUTO" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry = aurostd::substring2bool(AflowIn,_STROPT_+"NBANDS",TRUE);
    // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry << endl;
    // cerr << "NBANDS_EQUAL" << endl;
    // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme << endl;
    vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.options2entry(AflowIn,_STROPT_+"NBANDS=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"  
    // [OBSOLETE]  vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL_isentry  = aurostd::substring2bool(AflowIn,_STROPT_+"NBANDS=",TRUE);
    if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry) vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry=FALSE;
    // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry << endl;
    // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry << endl;
    // cerr << "vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.xscheme << endl;

    // cerr << "ENMAX_MULTIPLY" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.options2entry(AflowIn,_STROPT_+"ENMAX_MULTIPLY=",TRUE,vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0.0"  

    // RWIGS_STATIC
    // cerr << "RWIGS" << endl;
    vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC   =
      aurostd::substring2bool(AflowIn,_STROPT_+"RWIGS_STATIC",TRUE);

    // SPIN AND PRIORITIES // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;   // DEFAULT
    vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2;   // DEFAULT
    vflags.KBIN_VASP_FORCE_OPTION_SPIN.options2entry(AflowIn,_STROPT_+"SPIN=",DEFAULT_VASP_FORCE_OPTION_SPIN);
    if (vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) { //ME+RF20200225; fixes bug that SPIN was switched OFF in static calc. when SPIN=ON in aflow.in and REMOVE_RELAX was set by default
      if(!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option){
        if(aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string,"REMOVE_RELAX_1") || aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string,"REMOVE_RELAX_2")){
          pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "SPIN is OFF. REMOVE_RELAX_1/2 will be switched off.", aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
          vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=FALSE;
          vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=FALSE;
        }
      } else {
        vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=FALSE;
        vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=FALSE;
        if(aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string,"REMOVE_RELAX_1")){vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=TRUE;}
        if(aurostd::substring2bool(vflags.KBIN_VASP_FORCE_OPTION_SPIN.content_string,"REMOVE_RELAX_2")){vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=TRUE;}
      }
    }
    if(!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1=FALSE; // nothing to remove
    if(!vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2=FALSE; // nothing to remove

    // BADER AND PRIORITIES // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_BADER.options2entry(AflowIn,_STROPT_+"BADER=",DEFAULT_VASP_FORCE_OPTION_BADER);
    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry=TRUE; // DEFAULT 
    if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_RUN.flag("STATIC")) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_DIELECTRIC_STATIC) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_DIELECTRIC_DYNAMIC) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_DSCF) vflags.KBIN_VASP_FORCE_OPTION_BADER.option=TRUE; // DEFAULT 

    // ELF AND PRIORITIES // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_ELF.options2entry(AflowIn,_STROPT_+"ELF=",DEFAULT_VASP_FORCE_OPTION_ELF);
    //  if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry=TRUE; // DEFAULT 
    // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_FORCE_OPTION_ELF.option=TRUE; // DEFAULT 

    // AUTO_MAGMOM AND PRIORITIES  // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.options2entry(AflowIn,_STROPT_+"AUTO_MAGMOM=",DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);
    // LSCOUPLING AND PRIORITIES  // ON | OFF
    vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.options2entry(AflowIn,_STROPT_+"LSCOUPLING=",DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);
    if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option) {
      if(!aurostd::substring2bool(kflags.KBIN_BIN,"LS") && !aurostd::substring2bool(kflags.KBIN_BIN,"ls")) kflags.KBIN_BIN+="LS";
      if(!aurostd::substring2bool(kflags.KBIN_MPI_BIN,"LS") && !aurostd::substring2bool(kflags.KBIN_MPI_BIN,"ls")) kflags.KBIN_MPI_BIN+="LS";
      kflags.KBIN_BIN=aurostd::RemoveCharacter(kflags.KBIN_BIN,' ');            // if there is junk
      kflags.KBIN_MPI_BIN=aurostd::RemoveCharacter(kflags.KBIN_MPI_BIN,' ');    // if there is junk
    }
    // SYM AND PRIORITIES
    vflags.KBIN_VASP_FORCE_OPTION_SYM.options2entry(AflowIn,_STROPT_+"SYM=",DEFAULT_VASP_FORCE_OPTION_SYM);
    // WAVECAR AND PRIORITIES
    vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.options2entry(AflowIn,_STROPT_+"WAVECAR=",DEFAULT_VASP_FORCE_OPTION_WAVECAR);
    // CHGCAR AND PRIORITIES
    vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.options2entry(AflowIn,_STROPT_+"CHGCAR=",DEFAULT_VASP_FORCE_OPTION_CHGCAR);
    //ME20191028 - specify CHGCAR file to use
    vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.options2entry(AflowIn,_STROPT_+"CHGCAR_FILE=",0,"");

    // LDAU2 AND PRIORITIES
    vflags.KBIN_VASP_LDAU_SPECIES="";
    vflags.KBIN_VASP_LDAU_PARAMETERS="";
    vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag=TRUE;

    BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"LDAU1=","LDAU=");aurostd::StringSubst(BflowIn,"LDAU2=","LDAU=");
    vflags.KBIN_VASP_FORCE_OPTION_LDAU0.options2entry(BflowIn,string(_STROPT_+"LDAU=OFF"+"|"+_STROPT_+"LDAU=0"+"|"+_STROPT_+"LDAU=N"+"|"+_STROPT_+"LDAU=FALSE"));
    vflags.KBIN_VASP_FORCE_OPTION_LDAU1.options2entry(AflowIn,string(_STROPT_+"LDAU1=ON"+"|"+_STROPT_+"LDAU1=1"+"|"+"LDAU1=Y"+"|"+_STROPT_+"LDAU1=TRUE"+"|"+_STROPT_+"LDAU1=ADIABATIC"+"|"+_STROPT_+"LDAU1=CUTOFF"));
    vflags.KBIN_VASP_FORCE_OPTION_LDAU2.options2entry(AflowIn,string(_STROPT_+"LDAU2=ON"+"|"+_STROPT_+"LDAU2=1"+"|"+"LDAU2=Y"+"|"+_STROPT_+"LDAU2=TRUE"+"|"+_STROPT_+"LDAU2=ADIABATIC"+"|"+_STROPT_+"LDAU2=CUTOFF"));
    if(vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry)  vflags.KBIN_VASP_FORCE_OPTION_LDAU0.isentry=FALSE;
    if(vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {
      if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU_SPECIES=",TRUE))
        vflags.KBIN_VASP_LDAU_SPECIES=aurostd::substring2string(AflowIn,_STROPT_+"LDAU_SPECIES=",1,FALSE);
      if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU1_SPECIES=",TRUE))
        vflags.KBIN_VASP_LDAU_SPECIES=aurostd::substring2string(AflowIn,_STROPT_+"LDAU1_SPECIES=",1,FALSE);
      if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU2_SPECIES=",TRUE))
        vflags.KBIN_VASP_LDAU_SPECIES=aurostd::substring2string(AflowIn,_STROPT_+"LDAU2_SPECIES=",1,FALSE);
      if(aurostd::substring2bool(AflowIn,_STROPT_+"LDAU_PARAMETERS=",TRUE)) 
        vflags.KBIN_VASP_LDAU_PARAMETERS=RemoveWhiteSpaces(aurostd::substring2string(AflowIn,_STROPT_+"LDAU_PARAMETERS=",1,FALSE));
      if(vflags.KBIN_VASP_LDAU_SPECIES!="") 
        vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag=TRUE;
      if(vflags.KBIN_VASP_LDAU_PARAMETERS!="") 
        vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag=FALSE;
    }
    // ADIABATIC
    vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.options2entry(AflowIn,string(_STROPT_+"LDAU1=ADIABATIC"+"|"+_STROPT_+"LDAU2=ADIABATIC"+"|"+_STROPT_+"LDAU=ADIABATIC"));
    if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry) {
      if(vflags.KBIN_VASP_RUN_NRELAX<LDAU_ADIABATIC_RELAX_DEFAULT)
        vflags.KBIN_VASP_RUN_NRELAX=LDAU_ADIABATIC_RELAX_DEFAULT;
      vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int=vflags.KBIN_VASP_RUN_NRELAX;
    }
    // CUTOFF
    vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.options2entry(AflowIn,string(_STROPT_+"LDAU1=CUTOFF"+"|"+_STROPT_+"LDAU2=CUTOFF"+"|"+_STROPT_+"LDAU=CUTOFF"));
    if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry) {
      vflags.KBIN_VASP_RUN_NRELAX++;
    }
    // KPOINTS
    BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"KPOINTS_","KPOINTS="); // bypass for getting all "_"
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.options2entry(BflowIn,string(_STROPT_+"KPOINTS="),aurostd_xoptionMULTI,""); // stack them all
    if(0) {
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.content_string=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.content_string << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KEEPK\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KEEPK") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"EVEN\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("EVEN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"ODD\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("ODD") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSHIFT_GAMMA_EVEN\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_EVEN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSHIFT_GAMMA_ODD\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_ODD") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"GAMMA\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("GAMMA") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSCHEME_MONKHORST_PACK\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_MONKHORST_PACK") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSCHEME_GAMMA\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_GAMMA") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"KSCHEME_AUTO\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_AUTO") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag(\"IBZKPT\")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT") << endl;
    }

    // TYPE AND PRIORITIES // METAL | INSULATOR | SEMICONDUCTOR | DEFAULT
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.options2entry(AflowIn,_STROPT_+"TYPE=",FALSE,DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME);
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme=KBIN_WRONG_ENTRY_STRING;
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('M',"METAL");
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('I',"INSULATOR");
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('S',"SEMICONDUCTOR");
    vflags.KBIN_VASP_FORCE_OPTION_TYPE.scheme2scheme('D',"DEFAULT");
    if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry && vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme==KBIN_WRONG_ENTRY_STRING) {
      message = "vflags.KBIN_VASP_FORCE_OPTION_TYPE.content_string=" + vflags.KBIN_VASP_FORCE_OPTION_TYPE.content_string;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // PARAMETERS FOR INCAR
    // NSW=
    vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL =
      aurostd::substring2bool(AflowIn,_STROPT_+"NSW=",TRUE);
    if(vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL)
      vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE=aurostd::substring2utype<int>(AflowIn,_STROPT_+"NSW=");
    else
      vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE=0;

    // IGNORE_AFIX stuff
    BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"IGNORE_AFIX_","IGNORE_AFIX="); // bypass for getting all "_"
    vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.options2entry(BflowIn,string(_STROPT_+"IGNORE_AFIX="),aurostd_xoptionMULTI,""); // stack them all
    if(0) {
      //ERRORS
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:ALL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:ALL") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:BRMIX\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:BRMIX") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:CSLOSHING\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:CSLOSHING") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:CALC_FROZEN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:CALC_FROZEN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:DAV\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:DAV") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:DENTET\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:DENTET") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:EDDDAV\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:EDDDAV") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:EDDRMM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:EDDRMM") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:EFIELD_PEAD\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:EFIELD_PEAD") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:EXCCOR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:EXCCOR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:GAMMA_SHIFT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:GAMMA_SHIFT") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:IBZKPT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:IBZKPT") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:INVGRP\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:INVGRP") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:KKSYM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:KKSYM") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:LRF_COMMUTATOR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:LRF_COMMUTATOR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:MEMORY\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:MEMORY") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:MPICH0\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:MPICH0") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:MPICH11\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:MPICH11") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:MPICH139\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:MPICH139") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:MPICH174\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:MPICH174") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NATOMS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NATOMS") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NBANDS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NBANDS") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NELM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NELM") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NIRMAT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NIRMAT") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NKXYZ_IKPTD\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NKXYZ_IKPTD") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NPAR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NPAR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NPARC\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NPARC") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NPARN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NPARN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NPAR_REMOVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NPAR_REMOVE") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:NUM_PROB\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:NUM_PROB") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:OUTPUT_LARGE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:OUTPUT_LARGE") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:PSMAXN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:PSMAXN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:READ_KPOINTS_RD_SYM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:READ_KPOINTS_RD_SYM") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:REAL_OPT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:REAL_OPT") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:REAL_OPTLAY_1\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:REAL_OPTLAY_1") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:RMM_DIIS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:RMM_DIIS") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:ROTMAT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:ROTMAT") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:SGRCON\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:SGRCON") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:SYMPREC\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:SYMPREC") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:ZBRENT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:ZBRENT") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:ZPOTRF\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:ZPOTRF") << endl;
      //FIXES
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ALL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ALL") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ALGO\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ALGO") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ALGO=FAST\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ALGO=FAST") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ALGO=NORMAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ALGO=NORMAL") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ALGO=VERYFAST\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ALGO=VERYFAST") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:AMIN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:AMIN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:AMIN=0.01\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:AMIN=0.01") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:AMIX\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:AMIX") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:AMIX=0.1\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:AMIX=0.1") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:BMIX\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:BMIX") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:BMIX=0.01\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:BMIX=0.01") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:BMIX=3\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:BMIX=3") << endl;  //CO20210315 - what to do about 3 vs. 3.0?
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:EDIFF\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:EDIFF") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:EFIELD_PEAD\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:EFIELD_PEAD") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ENMAX\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ENMAX") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ISMEAR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ISMEAR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ISMEAR=2\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ISMEAR=2") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ISMEAR=-1\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ISMEAR=-1") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS=GAMMA\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS=GAMMA") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS=GAMMA_EVEN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS=GAMMA_EVEN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS=GAMMA_ODD\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS=GAMMA_ODD") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS=KMAX\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS=KMAX") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS--\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS--") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:KPOINTS-=2\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:KPOINTS-=2") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:LREAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:LREAL") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NBANDS\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NBANDS") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NBANDS++\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NBANDS++") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NBANDS--\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NBANDS--") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NPAR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NPAR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NPAR=1\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NPAR=1") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NPAR=4\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NPAR=4") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:NPAR_REMOVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:NPAR_REMOVE") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR_SCALE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR_SCALE") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR_SCALE*=1.2\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR_SCALE*=1.2") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR_VOLUME\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR_VOLUME") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR_VOLUME*=1.05\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR_VOLUME*=1.05") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR_VOLUME*=2\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR_VOLUME*=2") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POSCAR=STANDARD_CONVENTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POSCAR=STANDARD_CONVENTIONAL") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:POTIM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:POTIM") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:RECYCLE_CONTCAR\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:RECYCLE_CONTCAR") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:RELAX_MODE\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:RELAX_MODE") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:RELAX_MODE=FORCES\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:RELAX_MODE=FORCES") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:SKIP_RUN\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:SKIP_RUN") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:SYM\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:SYM") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:SYM=OFF\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:SYM=OFF") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:SYMPREC\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:SYMPREC") << endl;
      cerr << "vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"FIX:ULIMIT\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("FIX:ULIMIT") << endl;
    }

    // INPUT FILES
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_KEYWORD                =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]");      
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_SYSTEM_AUTO            =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]SYSTEM_AUTO",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_FILE                   =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]FILE=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_FILE_COMMAND                =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]COMMAND=",TRUE);

    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]")) vflags.KBIN_VASP_INCAR_FILE.push("KEYWORD");
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]SYSTEM_AUTO",TRUE)) vflags.KBIN_VASP_INCAR_FILE.push("SYSTEM_AUTO");
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]FILE=",TRUE)) { //ME20181113
      string file = aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]FILE=", 1, TRUE);
      vflags.KBIN_VASP_INCAR_FILE.push_attached("FILE", file);
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_FILE]COMMAND=",TRUE)) { //ME20181113
      string command = aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]COMMAND=", 1, TRUE);
      vflags.KBIN_VASP_INCAR_FILE.push_attached("COMMAND", command);
    }

    // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_EXPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]");
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_EXPLICIT_START_STOP  =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]STOP");
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_IMPLICIT]");
    // [OBSOLETE] vflags.KBIN_VASP_INCAR_MODE_EXTERNAL             =    aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXTERNAL]");
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]")) vflags.KBIN_VASP_INCAR_MODE.push("EXPLICIT");
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]STOP")) vflags.KBIN_VASP_INCAR_MODE.push("EXPLICIT_START_STOP");
    if (vflags.KBIN_VASP_INCAR_FILE.flag("KEYWORD")) { //ME20181113
      vflags.KBIN_VASP_INCAR_EXPLICIT.str(aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn, vflags.KBIN_VASP_INCAR_EXPLICIT, "[VASP_INCAR_FILE]");
    }
    if (vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT_START_STOP")) { //ME20181113
      vflags.KBIN_VASP_INCAR_EXPLICIT_START_STOP.str(aurostd::substring2string(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START", "[VASP_INCAR_MODE_EXPLICIT]STOP",0));
      //[TB20221206 - FIXED]vflags.KBIN_VASP_INCAR_EXPLICIT.str(aurostd::substring2string(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START", "[VASP_INCAR_MODE_EXPLICIT]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn, vflags.KBIN_VASP_INCAR_EXPLICIT_START_STOP, "[VASP_INCAR_MODE_EXPLICIT]START", "[VASP_INCAR_MODE_EXPLICIT]STOP");
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_IMPLICIT]")) vflags.KBIN_VASP_INCAR_MODE.push("IMPLICIT");
    if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXTERNAL]")) vflags.KBIN_VASP_INCAR_MODE.push("EXTERNAL");

    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_FILE_KEYWORD                      =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]");      
    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]")) vflags.KBIN_VASP_KPOINTS_FILE.push("KEYWORD");
    //  vflags.KBIN_VASP_KPOINTS_FILE_SYSTEM_AUTO          =  aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]SYSTEM_AUTO",TRUE);

    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_EXPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]");
    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_EXPLICIT_START_STOP  =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]STOP");
    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_IMPLICIT]");
    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_MODE_EXTERNAL             =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXTERNAL]");
    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]")) vflags.KBIN_VASP_KPOINTS_MODE.push("EXPLICIT");
    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]STOP")) vflags.KBIN_VASP_KPOINTS_MODE.push("EXPLICIT_START_STOP");
    if (vflags.KBIN_VASP_KPOINTS_FILE.flag("KEYWORD")) { //ME20181113
      vflags.KBIN_VASP_KPOINTS_EXPLICIT.str(aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn, vflags.KBIN_VASP_KPOINTS_EXPLICIT, "[VASP_KPOINTS_FILE]");
    }
    if (vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")) { //ME20181113
      vflags.KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(aurostd::substring2string(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]START", "[VASP_KPOINTS_MODE_EXPLICIT]STOP",0));
      //[TB20221206 - FIXED]vflags.KBIN_VASP_KPOINTS_EXPLICIT.str(aurostd::substring2string(AflowIn,"[VASP_KPOINTS_MODE_EXPLICIT]START", "[VASP_KPOINTS_MODE_EXPLICIT]STOP",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn, vflags.KBIN_VASP_KPOINTS_EXPLICIT_START_STOP, "[VASP_KPOINTS_MODE_EXPLICIT]START", "[VASP_KPOINTS_MODE_EXPLICIT]STOP");
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_IMPLICIT]")) vflags.KBIN_VASP_KPOINTS_MODE.push("IMPLICIT");
    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_MODE_EXTERNAL]")) vflags.KBIN_VASP_KPOINTS_MODE.push("EXTERNAL");
    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_FILE_FILE                 =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]FILE=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_FILE_COMMAND              =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=",TRUE);

    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]FILE=",TRUE)) { //ME20181113
      string file = aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]FILE=", 1, TRUE);
      vflags.KBIN_VASP_KPOINTS_FILE.push_attached("FILE", file);
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=",TRUE)) { //ME20181113
      string command = aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=", 1, TRUE);
      vflags.KBIN_VASP_KPOINTS_FILE.push_attached("COMMAND", command);
    }

    // KPOINTS FOR RELAX
    vflags.KBIN_VASP_KPOINTS_KMODE.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KMODE=",FALSE,vflags.KBIN_VASP_KPOINTS_KMODE.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"  
    // cerr << "vflags.KBIN_VASP_KPOINTS_KMODE.isentry=" << vflags.KBIN_VASP_KPOINTS_KMODE.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KMODE.xscheme=" << vflags.KBIN_VASP_KPOINTS_KMODE.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_KPPRA.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KPPRA=",FALSE,vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
    if(vflags.KBIN_VASP_KPOINTS_KPPRA.isentry==FALSE) {vflags.KBIN_VASP_KPOINTS_KPPRA.clear();vflags.KBIN_VASP_KPOINTS_KPPRA.isentry=TRUE;vflags.KBIN_VASP_KPOINTS_KPPRA.push("100");}
    // cerr << "vflags.KBIN_VASP_KPOINTS_KPPRA.isentry=" << vflags.KBIN_VASP_KPOINTS_KPPRA.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme=" << vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_KSCHEME.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KSCHEME=",FALSE,vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
    if(vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry==FALSE) {vflags.KBIN_VASP_KPOINTS_KSCHEME.clear();vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry=TRUE;vflags.KBIN_VASP_KPOINTS_KSCHEME.push(DEFAULT_KSCHEME);}
    // cerr << "vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry=" << vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme=" << vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_KSHIFT.options2entry(AflowIn,"[VASP_KPOINTS_FILE]KSHIFT=",FALSE,vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0 0 0"
    // cerr << "vflags.KBIN_VASP_KPOINTS_KSHIFT.isentry=" << vflags.KBIN_VASP_KPOINTS_KSHIFT.isentry << endl << "vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme=" << vflags.KBIN_VASP_KPOINTS_KSHIFT.xscheme << endl;

    // KPOINTS FOR STATIC
    vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KMODE=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0"
    // cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KPPRA=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
    // cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KSCHEME=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "Monkhorst-Pack"
    // cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.options2entry(AflowIn,"[VASP_KPOINTS_FILE]STATIC_KSHIFT=",FALSE,vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme); // scheme already loaded in aflow_xclasses.cpp is "0 0 0"
    //cerr << "vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry << endl << "vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.xscheme << endl;

    // [OBSOLETE] vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG        =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=",TRUE);
    // [OBSOLETE] if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG) vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE=aurostd::RemoveWhiteSpaces(aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=",1,TRUE));
    // [OBSOLETE]  else {vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE=DEFAULT_BANDS_LATTICE;vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG=TRUE;} // DEFAULT FIX

    vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.options2entry(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=",FALSE,vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.xscheme); // scheme already loaded in aflow_xclasses.cpp is AUTO
    if(!vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry){
      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();
      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(DEFAULT_BANDS_LATTICE);
    }

    if((vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) && !vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry) {
      cerr << "WARNING: if you use vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC_BANDS\") or vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\"), you must specify KBIN_VASP_KPOINTS_BANDS_LATTICE" << endl;
      cerr << "         Taking defauls vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
    }
    //[CO20210805 - OBSOLETE]vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG    =    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=AUTO",TRUE) ||    aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_LATTICE=auto",TRUE);
    // [OBSOLETE]  if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG) vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG=FALSE;
    //[CO20210805 - OBSOLETE]if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=="AUTO") vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG=TRUE;
    //[CO20210805 - OBSOLETE]if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG) vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();

    vflags.KBIN_VASP_KPOINTS_BANDS_GRID.options2entry(AflowIn,"[VASP_KPOINTS_FILE]BANDS_GRID=",FALSE,vflags.KBIN_VASP_KPOINTS_BANDS_GRID.xscheme); // scheme already loaded in aflow_xclasses.cpp is 20
    //[CO20210805 - OBSOLETE]vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG    =     aurostd::substring2bool(AflowIn,"[VASP_KPOINTS_FILE]BANDS_GRID=",TRUE);
    //[CO20210805 - OBSOLETE]if(vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG)
    //[CO20210805 - OBSOLETE]  vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=aurostd::substring2utype<int>(AflowIn,"[VASP_KPOINTS_FILE]BANDS_GRID=",TRUE);
    //[CO20210805 - OBSOLETE]else vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=16;vflags.KBIN_VASP_KPOINTS_BANDS_GRID_FLAG=TRUE;
    if(!vflags.KBIN_VASP_KPOINTS_BANDS_GRID.isentry){  //CO20210805
      vflags.KBIN_VASP_KPOINTS_BANDS_GRID.clear();
      vflags.KBIN_VASP_KPOINTS_BANDS_GRID.push(aurostd::utype2string(DEFAULT_BANDS_GRID));
    }
    //    cerr << "WARNING: if you use VASP_RUN_RELAX_STATIC_BANDS or vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\"), you must specify KBIN_VASP_KPOINTS_BANDS_GRID_FLAG" << endl;
    //   cerr << "         Taking defauls vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE=" << vflags.KBIN_VASP_KPOINTS_BANDS_GRID_VALUE << endl;
    if((vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || 
          vflags.KBIN_VASP_RUN.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("REPEAT_BANDS")) && 
        !vflags.KBIN_VASP_KPOINTS_BANDS_GRID.isentry) {
      message="";
      message+="if you run RELAX_STATIC_BANDS, STATIC_BANDS, REPEAT_STATIC_BANDS, or REPEAT_BANDS, you must specify KBIN_VASP_KPOINTS_BANDS_GRID";
      message+="taking default KBIN_VASP_KPOINTS_BANDS_GRID="+vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_string;
      pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,_LOGGER_WARNING_);
    }

    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_KEYWORD                 =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]"); 
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]"))  vflags.KBIN_VASP_POSCAR_FILE.push("KEYWORD");

    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT                  =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]");
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]")) vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT");

    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_START_STOP       =    aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_START_) &&  aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_STOP_);
    if(aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_START_) &&  aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_STOP_)) vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP");

    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_START_STOP_POINT =    (aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_START_P_) && aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_STOP_P_)); //CO20200624
    if((aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_START_P_) && aurostd::substring2bool(AflowIn,_VASP_POSCAR_MODE_EXPLICIT_STOP_P_))) vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP_POINT");  //CO20200624

    if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") && !kflags.KBIN_FROZSL) {  // NO FROZSL
      if(LDEBUG) cerr << "DEBUG: vflags.KBIN_VASP_POSCAR_MODE.flag(\"EXPLICIT_START_STOP_POINT\")=" << vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT") << endl;
      if(LDEBUG) cerr << "DEBUG: kflags.KBIN_PHONONS_CALCULATION_FROZSL=" << kflags.KBIN_PHONONS_CALCULATION_FROZSL << endl;
      if(LDEBUG) cerr << "DEBUG: kflags.KBIN_FROZSL_DOWNLOAD=" << kflags.KBIN_FROZSL_DOWNLOAD << endl;
      if(LDEBUG) cerr << "DEBUG: kflags.KBIN_FROZSL_FILE=" << kflags.KBIN_FROZSL_FILE << endl;
      stringstream input_file;
      input_file.clear();
      // loading
      if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) input_file.str(AflowIn);
      if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {
        FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags);
        FROZSL::Extract_INPUT(AflowIn,FileMESSAGE,input_file,aflags,kflags);
      }
      if(kflags.KBIN_FROZSL_DOWNLOAD)    FROZSL::Setup_frozsl_init_input(AflowIn,FileMESSAGE,input_file,aflags,kflags);
      if(kflags.KBIN_FROZSL_FILE)        FROZSL::File_INPUT(AflowIn,FileMESSAGE,input_file,aflags,kflags);

      vflags.KBIN_VASP_POSCAR_MODE.push("EXPLICIT_START_STOP_POINT");
      // done loading now load structures up
      vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP",FALSE); // some default
      aurostd::substring2strings(input_file.str(),vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING,_VASP_POSCAR_MODE_EXPLICIT_START_P_);  //CO20200624
      // some verbose
      for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++)
        if(LDEBUG) cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i) << endl;
      // load up the structures
      for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
        string START=_VASP_POSCAR_MODE_EXPLICIT_START_P_+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i); //CO20200624
        string STOP=_VASP_POSCAR_MODE_EXPLICIT_STOP_P_+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i); //CO20200624
        stringstream POSCAR;
        POSCAR.str(aurostd::substring2string(input_file.str(),START,STOP,-1));
        //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(input_file.str(),START) && aurostd::substring2bool(input_file.str(),STOP))
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(input_file.str(),POSCAR,START,STOP);
        if(!POSCAR.str().empty()) vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.push_back(xstructure(POSCAR,IOVASP_AUTO));
      }
      if(LDEBUG) cerr << "DEBUG " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << endl;
      if(LDEBUG) cerr << "DEBUG " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size() << endl;
      if(vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() != vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size()) {
        message = "IN " + _AFLOWIN_ + " in Directory=" + aflags.Directory + '\n';
        message += "vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()=" + aurostd::utype2string<uint>(vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) + '\n';
        message += "vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size()=" + aurostd::utype2string<uint>(vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.size());
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INDEX_MISMATCH_);
      }
      for(uint i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++)
        if(LDEBUG) cerr << "DEBUG= " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(i) << endl;
    } else {
      vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.clear();
      vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.clear();
    }
    // the rest for POSCAR
    // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_IMPLICIT]");
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_IMPLICIT]")) vflags.KBIN_VASP_POSCAR_MODE.push("IMPLICIT");

    // vflags.KBIN_VASP_POSCAR_FILE_SYSTEM_AUTO                =   aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]SYSTEM_AUTO",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_PROTOTYPE                     =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]PROTOTYPE=",TRUE);
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]PROTOTYPE=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE.push("PROTOTYPE");

    // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_MODE_EXTERNAL                      =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXTERNAL]");
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXTERNAL]")) vflags.KBIN_VASP_POSCAR_MODE.push("EXTERNAL");

    // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_FILE_FILE                          =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]FILE=",TRUE);
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]FILE=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE.push("FILE");
    // [OBSOLETE]  vflags.KBIN_VASP_POSCAR_FILE_COMMAND                       =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]COMMAND=",TRUE);
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]COMMAND=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE.push("COMMAND");

    // VOLUMES
    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_EQUAL_EQUAL            =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_PLUS_EQUAL             =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME+=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_MINUS_EQUAL            =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME-=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_MULTIPLY_EQUAL         =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME*=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POSCAR_FILE_VOLUME_DIVIDE_EQUAL           =    aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME/=",TRUE);
    vflags.KBIN_VASP_POSCAR_FILE_VOLUME.clear();
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("EQUAL_EQUAL");
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME+=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("PLUS_EQUAL");
    // [OBSOLETE] if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME-=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("MINUS_EQUAL");
    if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME*=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("MULTIPLY_EQUAL");
    // [OBSOLETE] if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_FILE]VOLUME/=",TRUE)) vflags.KBIN_VASP_POSCAR_FILE_VOLUME.push("DIVIDE_EQUAL");
    if(vflags.KBIN_VASP_POSCAR_FILE_VOLUME.xscheme!="") vflags.KBIN_VASP_POSCAR_FILE_VOLUME.isentry=TRUE;

    // CONVERT_UNIT_CELL stuff
    BflowIn=AflowIn;aurostd::StringSubst(BflowIn,"=","_");aurostd::StringSubst(BflowIn,"CONVERT_UNIT_CELL_","CONVERT_UNIT_CELL="); // bypass for getting all "_"
    vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.options2entry(BflowIn,string(_STROPT_+"CONVERT_UNIT_CELL="),aurostd_xoptionMULTI,""); // stack them all
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEFORE vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme << endl; //ME20181113
    // // PRIORITIES
    if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") || vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI",FALSE);
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL",FALSE);
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT",FALSE);
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ",FALSE);
    } // some PRIORITIES
    if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {
      vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL",FALSE);
    }
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " AFTER vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme << endl; //ME20181113

    // DEBUG
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string  << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"MINKOWSKI\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" <<  vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;

    // VOLUMES
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_EQUAL_EQUAL      =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_PLUS_EQUAL       =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME+=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_MINUS_EQUAL      =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME-=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_MULTIPLY_EQUAL   =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME*=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME_DIVIDE_EQUAL     =    aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME/=",TRUE);
    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.clear();

    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("EQUAL_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME+=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("PLUS_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME-=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("MINUS_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME*=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("MULTIPLY_EQUAL");
    // [OBSOLETE]  if(aurostd::substring2bool(AflowIn,_STROPT_+"VOLUME/=",TRUE)) vflags.KBIN_VASP_FORCE_OPTION_VOLUME.push("DIVIDE_EQUAL");

    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"EQUAL_EQUAL",_STROPT_+"VOLUME=","");
    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"PLUS_EQUAL",_STROPT_+"VOLUME+=","");
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"MINUS_EQUAL",_STROPT_+"VOLUME-=","");
    vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"MULTIPLY_EQUAL",_STROPT_+"VOLUME*=","");
    // [OBSOLETE] vflags.KBIN_VASP_FORCE_OPTION_VOLUME.args2addattachedscheme(vAflowIn,"DIVIDE_EQUAL",_STROPT_+"VOLUME/=","");
    if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("EQUAL_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"PLUS_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"PLUS_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("PLUS_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(LDEBUG) cerr << "CORMAC STUFF  " << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedscheme("MULTIPLY_EQUAL") << endl;

    if(vflags.KBIN_VASP_FORCE_OPTION_VOLUME.xscheme!="") vflags.KBIN_VASP_FORCE_OPTION_VOLUME.isentry=TRUE;

    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_KEYWORD              =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]");      
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_SYSTEM_AUTO          =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SYSTEM_AUTO",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_PREFIX               =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_SUFFIX               =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_FILE                 =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]FILE=",TRUE);
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_FILE_COMMAND              =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",TRUE);

    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]")) vflags.KBIN_VASP_POTCAR_FILE.push("KEYWORD");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SYSTEM_AUTO",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("SYSTEM_AUTO");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("PREFIX");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("SUFFIX");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]FILE=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("FILE");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",TRUE)) vflags.KBIN_VASP_POTCAR_FILE.push("COMMAND");

    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_MODE_EXPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXPLICIT]");
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_MODE_IMPLICIT             =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_IMPLICIT]");
    // [OBSOLETE] vflags.KBIN_VASP_POTCAR_MODE_EXTERNAL             =    aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXTERNAL]");

    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXPLICIT]")) vflags.KBIN_VASP_POTCAR_MODE.push("EXPLICIT");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_IMPLICIT]")) vflags.KBIN_VASP_POTCAR_MODE.push("IMPLICIT");
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_MODE_EXTERNAL]")) vflags.KBIN_VASP_POTCAR_MODE.push("EXTERNAL");

    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",TRUE)) { //CO20181113
      string file = aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]PREFIX=", 1, TRUE);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("PREFIX", file);
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",TRUE)) { //CO20181113
      string file = aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=", 1, TRUE);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("SUFFIX", file);
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]FILE=",TRUE)) { //CO20181113
      string file = aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]FILE=", 1, TRUE);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("FILE", file);
    }
    if(aurostd::substring2bool(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",TRUE)) { //CO20181113
      string command = aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]COMMAND=", 1, TRUE);
      vflags.KBIN_VASP_POTCAR_FILE.push_attached("COMMAND", command);
    }

    if (vflags.KBIN_VASP_POTCAR_FILE.flag("KEYWORD")) { //CO20181113
      vflags.KBIN_VASP_POTCAR_EXPLICIT.str(aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]",0));
      //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(AflowIn, vflags.KBIN_VASP_POTCAR_EXPLICIT, "[VASP_POTCAR_FILE]");
    }


    // APL ENTRIES
    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (APL)" << endl;

    //CO20170601 START
    //CO make backwards and forwards compatible with all possible workflows
    vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.options2entry(AflowIn,"[AFLOW_APL]KPPRA=|[AFLOW_QHA]KPPRA=|[AFLOW_AAPL]KPPRA=|[AFLOW_PHONONS]KPPRA=",
        FALSE,vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme); // scheme already loaded in aflow_xclasses.cpp is "1"
    // cerr << "vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry << endl << "vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.xscheme << endl;
    vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.options2entry(AflowIn,"[AFLOW_APL]KSCHEME=|[AFLOW_QHA]KSCHEME=|[AFLOW_AAPL]KSCHEME=|[AFLOW_PHONONS]KSCHEME=",
        FALSE,vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme); // scheme already loaded in aflow_xclasses.cpp is "DEFAULT_SCHEME"
    // cerr << "vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry << endl << "vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme=" << vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.xscheme << endl;


    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.options2entry(AflowIn,"[AFLOW_APL]KPOINTS=|[AFLOW_QHA]KPOINTS=|[AFLOW_AAPL]KPOINTS=",
        FALSE,vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.xscheme);
    //ME202020427 - APL k-point handling needs to be moved to modules eventually
    // This is a non-standard feature and should not be defaulted
    vflags.KBIN_VASP_KPOINTS_PHONONS_GRID.options2entry(AflowIn,"[AFLOW_APL]KPOINTS_GRID=|[AFLOW_QHA]KPOINTS_GRID=|[AFLOW_AAPL]KPOINTS_GRID=", FALSE,"");
    vflags.KBIN_VASP_KPOINTS_PHONONS_SHIFT.options2entry(AflowIn,"[AFLOW_APL]KPOINTS_SHIFT=|[AFLOW_QHA]KPOINTS_SHIFT=|[AFLOW_AAPL]KPOINTS_SHIFT=", FALSE,"");
    //[ME20181216]    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.clear();
    //[ME20181216]    if(aurostd::substring2bool(AflowIn,"[AFLOW_APL]KPOINTS=EVEN",TRUE) || 
    //[ME20181216]        aurostd::substring2bool(AflowIn,"[AFLOW_QHA]KPOINTS=EVEN",TRUE) ||
    //[ME20181216]        aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]KPOINTS=EVEN",TRUE) ||
    //[ME20181216]        aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_EVEN",TRUE)) 
    //[ME20181216]      {vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.push("EVEN");}
    //[ME20181216]    if(aurostd::substring2bool(AflowIn,"[AFLOW_APL]KPOINTS=ODD",TRUE) || 
    //[ME20181216]        aurostd::substring2bool(AflowIn,"[AFLOW_QHA]KPOINTS=ODD",TRUE) ||
    //[ME20181216]        aurostd::substring2bool(AflowIn,"[AFLOW_AAPL]KPOINTS=ODD",TRUE) ||
    //[ME20181216]        aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_ODD",TRUE)) 
    //[ME20181216]      {vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.push("ODD");}
    //[ME20181216]    //  vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY_EVEN=aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS=EVEN",TRUE)||aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_EVEN",TRUE);
    //[ME20181216]    // vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY_ODD=aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS=ODD",TRUE)||aurostd::substring2bool(AflowIn,"[AFLOW_PHONONS]KPOINTS_ODD",TRUE);
    //CO20170601 END

    // FROZSL ENTRIES
    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (FROZSL)" << endl;

    if(LDEBUG) cerr << "DEBUG: " << __AFLOW_FUNC__ << " (STOP)" << endl;

    return vflags;
  }
} // namespace KBIN

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
namespace KBIN {
  bool VASP_ExtractNGF(string OUTCAR,int &NGXF,int &NGYF,int &NGZF);
} // namespace KBIN

namespace KBIN {
  bool VASP_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    if(LDEBUG) cerr << "DEBUG: KBIN::VASP_Directory (BEGIN)" << endl;
    //  bool KBIN_MPI_LOCAL;KBIN_MPI_LOCAL=MPI;
    // bool KBIN_VASP_WRITE_KPOINTS;
    // string::size_type sub_size1,sub_size2;
    string subS,subS1,subS2;
    ostringstream aus;
    string::iterator pos;
    bool Krun=TRUE;

    ifstream FileAFLOWIN;
    string FileNameAFLOWIN;
    string AflowIn;
    FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
    FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
    FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
    AflowIn="";char c; while (FileAFLOWIN.get(c)) AflowIn+=c;  // READ _AFLOWIN_ and put into AflowIn
    FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
    AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
      aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      return FALSE;
    }
    aflags.QUIET=FALSE;
    _vflags vflags;
    vflags=KBIN::VASP_Get_Vflags_from_AflowIN(AflowIn,FileMESSAGE,aflags,kflags);

    // *********************************************************************************************************************
    // OPERATIONS related to PARTICULAR MACHINES ***************************************************************************

    if(LDEBUG) cerr << "[DEBUG] aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL.getattachedscheme("NAME") << endl;
    if(LDEBUG) cerr << "[DEBUG] aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << endl; //HE20220309 use machine name

    // ***************************************************************************
    if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") << endl; //CO20210315
    if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") << endl;
    if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS") << endl;
    if(LDEBUG) cerr << "[DEBUG] vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_DELSOL\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL") << endl;

    // ***************************************************************************
    // Get the KBIN_BIN name
    aurostd::StringstreamClean(aus);
    aus << "00000  MESSAGE KBIN::VASP_Directory Running KBIN_BIN=\"" << kflags.KBIN_BIN << "\"" << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // ***************************************************************************
    // Some verbose
    if(kflags.KBIN_POCC) {aus << "00000  MESSAGE KBIN::VASP_Directory Running POCC_CALCULATION" << Message(__AFLOW_FILE__,aflags) << endl;} //CO20180419 //POCC is special needs to run first because there is NO poscar defined yet
    else if(kflags.KBIN_PHONONS_CALCULATION_APL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_APL" << Message(__AFLOW_FILE__,aflags) << endl;
    else if(kflags.KBIN_PHONONS_CALCULATION_QHA) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_QHA" << Message(__AFLOW_FILE__,aflags) << endl;   //CO20170601
    else if(kflags.KBIN_PHONONS_CALCULATION_AAPL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AAPL" << Message(__AFLOW_FILE__,aflags) << endl; //CO20170601
    else if(kflags.KBIN_PHONONS_CALCULATION_AGL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AGL (Debye Model)" << Message(__AFLOW_FILE__,aflags) << endl;
    else if(kflags.KBIN_PHONONS_CALCULATION_AEL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_AEL (Elastic constants)" << Message(__AFLOW_FILE__,aflags) << endl;
    else if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) aus << "00000  MESSAGE KBIN::VASP_Directory Running PHONONS_CALCULATION_FROZSL" << Message(__AFLOW_FILE__,aflags) << endl;
    else {
      if(vflags.KBIN_VASP_RUN.flag("GENERATE")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_GENERATE" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("STATIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_KPOINTS" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("RELAX")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_STATIC_BANDS" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_RELAX_STATIC_BANDS" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_DIELECTRIC_STATIC" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_DIELECTRIC_DYNAMIC" << Message(__AFLOW_FILE__,aflags) << endl;
      if(vflags.KBIN_VASP_RUN.flag("DSCF")) aus << "00000  MESSAGE KBIN::VASP_Directory Running RUN_DSCF" << Message(__AFLOW_FILE__,aflags) << endl;
    }
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // ***************************************************************************
    uint ntasks=0;
    ntasks=1; // default
    if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
      ntasks=vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()+1;  //CO20200624 - include head directory as well (at the end)
      aus << "00000  MESSAGE Loaded ntasks = " << ntasks << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      uint i=0;
      for(i=0;i<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size();i++) {
        aus << "00000  MESSAGE task " << i+1 << "/" << ntasks << " in subdirectory " << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(i) << endl; //CO20200624 - +1
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      aus << "00000  MESSAGE task " << ((i++)+1) << "/" << ntasks << " in main directory " << aflags.Directory << endl;  //CO20200624 - include head directory as well (at the end)
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // ***************************************************************************
    // start the loop !
    _aflags aflags_backup;aflags_backup=aflags;
    _kflags kflags_backup;kflags_backup=kflags;
    //  _vflags vflags_backup;vflags_backup=vflags;


    for(uint ixvasp=0;ixvasp<ntasks;ixvasp++) {  // LOOP ixvasp
      // declarations
      _xvasp xvasp;xvasp.clear();
      xvasp.POSCAR_index=ixvasp;
      aflags=aflags_backup;kflags=kflags_backup; // load it up
      readModulesFromAflowIn(AflowIn, kflags, xvasp);  //ME20181027
      // some verbose
      if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")&&ixvasp<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) { //CO20200624 - include head directory as well (at the end)
        aus << "00000  MESSAGE START loop " << xvasp.POSCAR_index+1 << "/" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size() << Message(__AFLOW_FILE__,aflags) << endl;  //CO20200624 - +1
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [1]" << xvasp.str << endl; 
      // ------------------------------------------
      // now start for each xvasp
      Krun=TRUE;  // guess everything is intelligent !
      xvasp.Directory=aflags.Directory;
      if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")&&ixvasp<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) {  //CO20200624 - include head directory as well (at the end)
        xvasp.Directory=aflags.Directory+"/"+KBIN_SUBDIRECTORIES+vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.at(xvasp.POSCAR_index);
        aus << "00000  MESSAGE Taking loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // check for directory KY CHECK THIS (if Krun=FALSE, everything stops).
      if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")&&ixvasp<vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()) { //CO20200624 - include head directory as well (at the end)
        if(aurostd::FileExist(xvasp.Directory)) {
          Krun=FALSE; // avoid rerunning
          aus << "00000  MESSAGE Skipping loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          // before making it, check it again... NFS problem... check LOCK again
          if(Krun && aurostd::FileExist(xvasp.Directory+"/"+_AFLOWLOCK_)) Krun=FALSE;    // to fight against NFS cache
          if(Krun && aurostd::EFileExist(xvasp.Directory+"/"+_AFLOWLOCK_)) Krun=FALSE;     // to fight against NFS cache
          if(Krun && aurostd::FileExist(xvasp.Directory+"/LLOCK")) Krun=FALSE;     // to fight against NFS cache
          if(Krun && aurostd::EFileExist(xvasp.Directory+"/LLOCK")) Krun=FALSE;     // to fight against NFS cache
          if(Krun) {
            aurostd::DirectoryMake(xvasp.Directory);
            aus << "00000  MESSAGE Creating loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            aurostd::execute("echo \"NNNNN  KBIN LLOCK ASAP for NFS concurrent jobs (aflow"+string(AFLOW_VERSION)+")\" >> "+xvasp.Directory+"/LLOCK");
          }
        }
      }


      if(Krun) {
        aflags.Directory=xvasp.Directory; // so we are set ! since there are plenty of routines with aflags.Directory inside
        aus << "00000  MESSAGE Performing loop directory = " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // ------------------------------------------
      // do the flags
      if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [2]" << xvasp.str << endl;
      vflags.KBIN_VASP_INCAR_VERBOSE=TRUE; // ALWAYS

      if(0){  //CO20210315 - better to keep verbosity ON
        if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
        if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
        if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) vflags.KBIN_VASP_INCAR_VERBOSE=FALSE; // TURN OFF VERBOSITY
      }

      if(Krun && vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option) {  //CO20211217 - prevents VASP_Backup() from deleting and recycles for next run (VASP_RecycleExtraFile())
        xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED",vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option); //passing flag to xvasp
        aus << "00000  MESSAGE Saving WAVECAR files" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }

      //CO, come back here at some point
      //think about taking pocc structure and reducing first if possible (and requested)
      //need to first convert to non-pocc structure (already programmed in aflow_pocc.cpp)
      //plug this in as xvasp.str (also create xvasp.str_pocc then)
      //reduce as requested and re-pocc the structure
      //think about if we need a separate flag for reducing pocc vs. reducing derivative structures
      // produce BEFORE NOMIX
      if(!(
            kflags.KBIN_POCC ||
            //kflags.KBIN_PHONONS_CALCULATION_APL ||    //CO20180515 - KEEP APL/QHA/AAPL to grab structure
            //kflags.KBIN_PHONONS_CALCULATION_QHA ||    //CO20180515 - KEEP APL/QHA/AAPL to grab structure
            //kflags.KBIN_PHONONS_CALCULATION_AAPL ||   //CO20180515 - KEEP APL/QHA/AAPL to grab structure
            kflags.KBIN_PHONONS_CALCULATION_FROZSL ||   //CO20180503 - KEEP AEL/AGL stuff running per normal
            //kflags.KBIN_PHONONS_CALCULATION_AGL ||    //CO20180503 - KEEP AEL/AGL stuff running per normal
            //kflags.KBIN_PHONONS_CALCULATION_AEL ||    //CO20180503 - KEEP AEL/AGL stuff running per normal
            FALSE)  //identity
        ) { //CO20180419, do NOT produce POSCAR for POCC
        if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
        if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
        if(Krun && kflags.KBIN_QSUB) Krun=(Krun && KBIN::QSUB_Extract(xvasp.xqsub,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags));
        if(Krun && kflags.KBIN_QSUB_MODE1) Krun=(Krun && KBIN::QSUB_Extract_Mode1(xvasp.xqsub,FileMESSAGE,aflags,kflags));
        if(Krun && kflags.KBIN_QSUB_MODE2) Krun=(Krun && KBIN::QSUB_Extract_Mode2(xvasp.xqsub,FileMESSAGE,aflags,kflags));
        if(Krun && kflags.KBIN_QSUB_MODE3) Krun=(Krun && KBIN::QSUB_Extract_Mode3(xvasp.xqsub,FileMESSAGE,aflags,kflags));
      }
      //ME20210709 [OBSOLETE] if(Krun && vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.isentry) //DX+ME20210709 - Throws an error with --genertate_aflowin_only because it does not populate the POTCARs, so the miscibility check fails.
      if(Krun && !XHOST.GENERATE_AFLOWIN_ONLY && vflags.KBIN_VASP_FORCE_OPTION_SKIP_NOMIX.isentry) { //ME20210709 - For now, skip check if generate_aflowin_only. In the future, the elements (not the pseudopotentials) need to be grabbed from the aflow.in
        string potentials=xvasp.POTCAR_POTENTIALS.str();
        if(!aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/1/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/2/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/3/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/58/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/59/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/60/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/115/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/116/") &&
            !aurostd::substring2bool(aurostd::CleanFileName(xvasp.Directory+"/"),"/117/")
          ) {
          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SKIP_NOMIX (NEGLECT_NOMIX, NEGLECT_IMMISCIBLE)" << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          //	cerr << "potentials=" << potentials << endl;
          if(MiscibilityCheck(potentials)==MISCIBILITY_SYSTEM_NOMIX) {
            aus << "00000  MESSAGE Skipping system: " << KBIN::VASP_PseudoPotential_CleanName(potentials) << " is known to be immiscible (aflow_nomix.cpp)" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            stringstream command("");
            command << "cat " << xvasp.Directory << "/" << _AFLOWLOCK_ << " > " << xvasp.Directory << "/" << DEFAULT_AFLOW_IMMISCIBILITY_OUT << endl;
            aurostd::execute(command);
            Krun=FALSE;
          }
        }
        if(MiscibilityCheck(potentials)==MISCIBILITY_SYSTEM_MISCIBLE) {
          aus << "00000  MESSAGE Running system: " << KBIN::VASP_PseudoPotential_CleanName(potentials) << " is known to be miscible (aflow_nomix.cpp)" << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=TRUE;
        }
        if(MiscibilityCheck(potentials)==MISCIBILITY_SYSTEM_UNKNOWN) {
          aus << "00000  MESSAGE Running system: " << KBIN::VASP_PseudoPotential_CleanName(potentials) << " is unknown (aflow_nomix.cpp)" << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=TRUE;
        }
      } 

      // produce AFTER NOMIX
      // if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags));
      // if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
      // if(Krun && kflags.KBIN_QSUB) Krun=(Krun && KBIN::QSUB_Extract(xvasp.xqsub,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags));
      // if(Krun && kflags.KBIN_QSUB_MODE1) Krun=(Krun && KBIN::QSUB_Extract_Mode1(xvasp.xqsub,FileMESSAGE,aflags,kflags));
      // if(Krun && kflags.KBIN_QSUB_MODE2) Krun=(Krun && KBIN::QSUB_Extract_Mode2(xvasp.xqsub,FileMESSAGE,aflags,kflags));
      // if(Krun && kflags.KBIN_QSUB_MODE3) Krun=(Krun && KBIN::QSUB_Extract_Mode3(xvasp.xqsub,FileMESSAGE,aflags,kflags));


      // ***************************************************************************
      // READY TO RUN
      if(Krun) {
        if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [3]" << endl;
        if(LDEBUG) cerr << xvasp.str << endl;
        xvasp.NRELAX=0;
        bool Krun=true;
        ostringstream aus;
        bool PAWGGA2=FALSE;
        // ***************************************************************************
        // directory check
        ifstream DirectoryStream;
        DirectoryStream.open(xvasp.Directory.c_str(),std::ios::in);
        if(!DirectoryStream) {
          //   aus << "EEEEE  DIRECTORY_NOT_FOUND = " << Message(__AFLOW_FILE__,aflags) << endl;
          aus << "XXXXX  MAKING DIRECTORY = " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
          string str="mkdir "+xvasp.Directory;
          system(str.c_str());
        }
        // some check
        if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
          //    aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = " << Message(__AFLOW_FILE__,aflags) << endl;
          //    aurostd::PrintMessageStream(aus,XHOST.QUIET);
          //    return FALSE;
        }
        // ***************************************************************************
        // DO THE SYMMETRY NEIGHBORS CALCULATION
        //if(!kflags.KBIN_PHONONS_CALCULATION_FROZSL) { //[CO20200106 - close bracket for indenting]}
        //DX

        if(!(
              kflags.KBIN_POCC ||
              kflags.KBIN_PHONONS_CALCULATION_APL ||
              kflags.KBIN_PHONONS_CALCULATION_QHA || 
              kflags.KBIN_PHONONS_CALCULATION_AAPL || 
              kflags.KBIN_PHONONS_CALCULATION_FROZSL ||     //CO20180503 - KEEP AEL/AGL stuff running per normal
              //kflags.KBIN_PHONONS_CALCULATION_AGL ||      //CO20180503 - KEEP AEL/AGL stuff running per normal
              //kflags.KBIN_PHONONS_CALCULATION_AEL ||      //CO20180503 - KEEP AEL/AGL stuff running per normal
              FALSE) || //CO20170601 //identity
            aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN
          ) {  //CO, do internally
          //DX
          if(Krun) Krun=KBIN_StepSymmetryPerform(xvasp.str,AflowIn,FileMESSAGE,aflags,kflags,TRUE,cout); // DO THE SYMMETRY CALCULATION
          //DX20210122 [OBSOLETE - function doesn't calculate anything, removed] if(Krun) Krun=StepNeighborsPerform(xvasp.str,AflowIn,FileMESSAGE,aflags,kflags); // DO THE NEIGHBORS CALCULATION
          //DX
          //cerr << "KBIN GEN SYMMETRY OF AFLOWIN: " << aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN << endl;
          if(aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN){
            return Krun;
          }
          //DX
        }
        // VASP VASP WRITE
        //   if(Krun) Krun=(Krun && KBIN::VASP_Write_INPUT(xvasp,vflags));
        // ***************************************************************************
        // VASP INPUT FILES ARE DONE, NOW WE CAN USE OR MODYFYING THEM
        if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry) {
          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NOTUNE, no tuning xCARs - ";
          aus << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - ";
          aus << XHOST.hostname << " - " << aflow_get_time_string() << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }	  
        // ***************************************************************************
        // VASP HOW TO RUN ??
        // GENERATE ONLY -------------------------------------------------------------
        if(vflags.KBIN_VASP_RUN.flag("GENERATE")) {
          KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
          aus << "00000  MESSAGE VASP generation files ONLY" << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          xvasp.NRELAX=0;
        } else {
          //CO20180420 - added if/else-if for workflows that need to PRECEDE relax/static/etc.
          // RUN SOMETHING
          if(kflags.KBIN_POCC) {  // RUN POCC ------------------------  //CO20180419 //POCC is special, run as priority
            aus << "00000  MESSAGE PERFORMING POCC_CALCULATION" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
          }else if(kflags.KBIN_PHONONS_CALCULATION_APL) {  // RUN PHONONS APL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_APL" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
            //CO20170601 START
          }else if(kflags.KBIN_PHONONS_CALCULATION_QHA) {  // RUN PHONONS QHA ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_QHA" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
          }else if(kflags.KBIN_PHONONS_CALCULATION_AAPL) {  // RUN PHONONS AAPL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AAPL" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
            //CO20170601 END
          }else if(kflags.KBIN_PHONONS_CALCULATION_AGL) {  // RUN PHONONS AGL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AGL (Debye Model)" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
          }else if(kflags.KBIN_PHONONS_CALCULATION_AEL) {  // RUN PHONONS AEL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_AEL (Elastic constants)" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
          }else if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {  // RUN PHONONS FROZSL ------------------------
            aus << "00000  MESSAGE PERFORMING PHONONS_CALCULATION_FROZSL" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.NRELAX=-3;
          }else {
            if(vflags.KBIN_VASP_RUN.flag("STATIC")) {  // RUN STATIC ------------------------
              aus << "00000  MESSAGE Performing Static RUN" << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              //	if(vflags.KBIN_VASP_KPOINTS_KMODE_isentry==TRUE || vflags.KBIN_VASP_KPOINTS_KSCHEME_isentry==TRUE || vflags.KBIN_VASP_KPOINTS_KPPRA_isentry==TRUE || vflags.KBIN_VASP_KPOINTS_KSHIFT_isentry) {
              //	  aus << "00000  MESSAGE Patching KPOINT for the Static RUN" << Message(__AFLOW_FILE__,aflags) << endl;
              //	  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              //	}
              xvasp.NRELAX=-1;
            }
            if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) {  // RUN KPOINTS ------------------------
              aus << "00000  MESSAGE Running KPOINTS swap" << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              xvasp.NRELAX=-2;
            }
            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) {  // RUN RELAX_STATIC_BANDS ------------------------
              xvasp.NRELAX=vflags.KBIN_VASP_RUN_NRELAX; //  aurostd::substring2utype<int>(AflowIn,"[VASP_RUN_RELAX_STATIC_BANDS=");
              if(xvasp.NRELAX<0)  {
                aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
                Krun=FALSE;	   //	  FileINPUT.clear();FileINPUT.close();FileMESSAGE.clear();FileMESSAGE.close();
                xvasp.NRELAX=0;
              }
              //	if(xvasp.NRELAX>1 && xvasp.NRELAX!=2)
              {
                aus << "00000  MESSAGE RELAX_STATIC_BANDS Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
              }
            }
            if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {  // RUN RELAX_STATIC ------------------------
              xvasp.NRELAX=vflags.KBIN_VASP_RUN_NRELAX; //  aurostd::substring2utype<int>(AflowIn,"[VASP_RUN_RELAX_STATIC=");
              if(xvasp.NRELAX<0)  {
                aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
                Krun=FALSE;	   //	  FileINPUT.clear();FileINPUT.close();FileMESSAGE.clear();FileMESSAGE.close();
                xvasp.NRELAX=0;
              }
              //	if(xvasp.NRELAX>1 && xvasp.NRELAX!=2)
              {
                aus << "00000  MESSAGE RELAX_STATIC Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
              }
            }
            if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) { // RUN STATIC_BANDS ------------------------
              xvasp.NRELAX=-1;	
              aus << "00000  MESSAGE STATIC_BANDS Running " << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
              //	  xvasp.NRELAX=0;//nrelax;	
            }
            if(vflags.KBIN_VASP_RUN.flag("RELAX")) { // RUN RELAX ------------------------
              if(!(aurostd::substring2bool(AflowIn,"[VASP_RUN_RELAX=") || aurostd::substring2bool(AflowIn,"[VASP_RUN]RELAX="))) {
                xvasp.NRELAX=2;
                aus << "00000  MESSAGE Running DEFAULT [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              } else  { 	  
                xvasp.NRELAX=vflags.KBIN_VASP_RUN_NRELAX; //  aurostd::substring2utype<int>(AflowIn,"[VASP_RUN_RELAX=");
              }
              if(xvasp.NRELAX==0 || xvasp.NRELAX<0)  {
                aus << "EEEEE  No relaxation to run or nrelax<0 [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);    
                Krun=FALSE;	   //	  FileINPUT.clear();FileINPUT.close();FileMESSAGE.clear();FileMESSAGE.close();
                xvasp.NRELAX=0;
              }
              aus << "00000  MESSAGE RELAX Running [nrelax=" << xvasp.NRELAX << "] " << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
            }
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) { // RUN REPEAT_STATIC ------------------------ //CO20210315
              xvasp.NRELAX=-1;	
              aus << "00000  MESSAGE REPEAT_STATIC Running " << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
            }
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) { // RUN REPEAT_STATIC_BANDS ------------------------
              xvasp.NRELAX=-1;	
              aus << "00000  MESSAGE REPEAT_STATIC_BANDS Running " << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
            }
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) { // RUN REPEAT_BANDS ------------------------
              xvasp.NRELAX=-1;	
              aus << "00000  MESSAGE REPEAT_BANDS Running " << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
            }
            if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) { // RUN REPEAT_DELSOL ------------------------
              xvasp.NRELAX=-1;	
              aus << "00000  MESSAGE REPEAT_DELSOL Running " << Message(__AFLOW_FILE__,aflags) << endl;
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
            }
          }

          // ***************************************************************************
          // READY TO RUN
          if(Krun) {   // survived all troubles
            // ***************************************************************************
            // START
            if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [4]" << xvasp.str << endl;
            // ***************************************************************************
            // FIX BLANC SPECIES
            if(xvasp.str.species.size()>0) {
              //[OBSOLETE] CO, fixing for RHT routines, FIXED INSIDE RHT
              //      if(xvasp.str.species.at(0)=="A") {
              //	      for(uint itype=0;itype<xvasp.str.species.size();itype++) {
              //    xvasp.str.species.at(itype)="";
              //  }
              //}
              //CO, fixing for RHT routines
              if(xvasp.str.species.at(0)=="") {
                pflow::fixEmptyAtomNames(xvasp.str);  //CO moved to pflow
                //  for(uint itype=0;itype<xvasp.str.species.size();itype++) {
                //    if(xvasp.str.species.size()==xvasp.str.species_pp.size()) {
                //      if((xvasp.str.species.at(itype)=="") && xvasp.str.species_pp.at(itype)!="") 
                //        xvasp.str.species.at(itype)=KBIN::VASP_PseudoPotential_CleanName(xvasp.str.species_pp.at(itype));
                //    }
                //  }  // cormac I`ll write a short pflow for this stuff
                //  int iatom=0;
                //  for(uint itype=0;itype<xvasp.str.num_each_type.size();itype++) {
                //    string species=string(xvasp.str.species.at(itype));
                //    xvasp.str.species.at(itype)=species;
                //    for(int j=0;j<xvasp.str.num_each_type.at(itype);j++) {
                //      xvasp.str.atoms.at(iatom).name=species;    // CONVASP_MODE
                //      xvasp.str.atoms.at(iatom).CleanName();
                //      xvasp.str.atoms.at(iatom).CleanSpin();
                //      xvasp.str.atoms.at(iatom).name_is_given=TRUE;
                //      iatom++;
                //    }
                //  }
              }
            }
            // ***************************************************************************
            // PRESCRIPT
            if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
              KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_PRESCRIPT_COMMAND,DEFAULT_AFLOW_PRESCRIPT_OUT);
            // ***************************************************************************
            //CO20180419 - POCC always comes first (NO POSCAR), need to convert PARTCAR -> POSCARs
            //other workflows follow, all of these precede relaxation/static/etc.
            if(kflags.KBIN_POCC){
              //[CO20200624 - OBSOLETE]if (kflags.KBIN_PHONONS_CALCULATION_AEL) xvasp.aopts.push_attached("AFLOWIN_FLAG::MODULE", "AEL");  //CT20200319
              //[CO20200624 - OBSOLETE]if (kflags.KBIN_PHONONS_CALCULATION_AGL) xvasp.aopts.push_attached("AFLOWIN_FLAG::MODULE", "AGL");  //CT20200319
              KBIN::VASP_RunPOCC(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);
            } //CO20180419
            else if(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL) {
              //ME20200107 - Wrap in a try statement so that faulty APL runs don't kill other post-processing
              try {
                KBIN::VASP_RunPhonons_APL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE); // PHONONIC PHONONIC PHONONIC //CO20170601
              } catch (aurostd::xerror e) {
                pflow::logger(e.whereFileName(), e.whereFunction(), e.buildMessageString(), aflags.Directory, FileMESSAGE, std::cout, _LOGGER_ERROR_);
              }
            }
            else if(kflags.KBIN_PHONONS_CALCULATION_AGL==TRUE) {KBIN::VASP_RunPhonons_AGL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);}
            else if(kflags.KBIN_PHONONS_CALCULATION_AEL==TRUE) {KBIN::VASP_RunPhonons_AEL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);}
            else if(kflags.KBIN_PHONONS_CALCULATION_FROZSL) {KBIN::VASP_RunPhonons_FROZSL(xvasp,AflowIn,aflags,kflags,vflags,FileMESSAGE);}
            else {
              if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species.size()=" << xvasp.str.species.size() << endl;
              if(LDEBUG) for(uint i=0;i<xvasp.str.species.size();i++) cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species.at(i)=[" << xvasp.str.species.at(i) << "]" << endl;
              if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species_pp.size()=" << xvasp.str.species_pp.size() << endl;
              if(LDEBUG) for(uint i=0;i<xvasp.str.species_pp.size();i++) cerr << XPID << "KBIN::VASP_Directory: [5] xvasp.str.species_pp.at(i)=[" << xvasp.str.species_pp.at(i) << "]" << endl;
              //	    KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
              //	    cerr << xvasp.POTCAR.str() << endl;
              if(LDEBUG) cerr << XPID << "KBIN::VASP_Directory: [6]" << xvasp.str << endl;
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // STATIC STATIC STATIC
              if(vflags.KBIN_VASP_RUN.flag("STATIC")) {    // xvasp.RELAX=-1
                xvasp.aopts.flag("FLAG::POSCAR_PRESERVED",TRUE); // in case of errors it is not lost but recycled
                KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                aus << 11111 << "  STATIC - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [STATIC]");return Krun;}
                //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("static"));
                Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [STATIC]");return Krun;} //CO20201111 //AFTER CONTCAR_SAVE_
                bool qmwrite=TRUE;
                KBIN::VASP_Backup(xvasp,qmwrite,string("static"));
              }		
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // RELAX RELAX RELAX
              if(vflags.KBIN_VASP_RUN.flag("RELAX")) {    // xvasp.RELAX>0
                KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                if(PAWGGA2) {  // WAS A BUG IN PAW MAYBE IT IS FIXED
                  // STEP 1
                  aus << "11111  RELAXATION - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax2paw_gga",TRUE,FileMESSAGE);
                  if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [PAWGGA2 REL]");return Krun;}
                  //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                  //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [PAWGGA2 REL]");return Krun;} //CO20201111
                  aus << "22222  END        - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                } else {
                  if(xvasp.NRELAX==0) {
                    cerr << "STATIC RUN FIX INCAR: should not be here" << endl;
                    return FALSE;
                  } else { // DYNAMIC RUN
                    for(xvasp.NRELAXING=1;xvasp.NRELAXING<=xvasp.NRELAX;xvasp.NRELAXING++) {
                      aus << 11111*xvasp.NRELAXING << "  RELAXATION - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                      if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int>0) KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp,xvasp.NRELAXING); // ADIABATIC
                      if(xvasp.NRELAXING<xvasp.NRELAX)  {
                        Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                        if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAXATION<]");return Krun;}
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RELAXATION<]");return Krun;} //CO20201111
                        KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAXING,true,FileMESSAGE);         // check if it is the case of turning off spin //CO20210315 - always write_incar here
                        KBIN::XVASP_KPOINTS_IBZKPT_UPDATE(xvasp,aflags,vflags,xvasp.NRELAXING,true,FileMESSAGE);           // check if it is the case of updating IBZKPT  //CO20210315 - always write_incar here
                        //ME20190301 BEGIN
                        // CHGCAR/WAVECAR needs to be recycled if CHGCAR/WAVECAR=ON or VASP
                        // won't be able to read the files. Bug found by Rico Friedrich
                        if(vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option) {
                          //ME20191031 - Only recycle when ICHARG was found
                          string extra_incar = xvasp.AVASP_EXTRA_INCAR.str();
                          int nlines = aurostd::GetNLinesString(extra_incar);
                          int l;
                          string line;
                          for (l = 1; l <= nlines; l++) {
                            line = aurostd::RemoveWhiteSpaces(aurostd::GetLineString(extra_incar, l));
                            if (aurostd::substring2bool(line, "ICHARG=1", true)) break;
                          }
                          if (l <= nlines) KBIN::VASP_RecycleExtraFile(xvasp, "CHGCAR", "relax"+aurostd::utype2string<int>(xvasp.NRELAXING));
                        }
                        //CO+MM20211217
                        //In general it's not a good idea to use the WAVECAR from relax1 to start relax2
                        //a big change in the unit cell size/shape is going drastically change the plane waves that should be used
                        //we can add another option to do this in the future
                        //[CO+MM20211217]if(vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option) KBIN::VASP_RecycleExtraFile(xvasp, "WAVECAR", "relax"+aurostd::utype2string<int>(xvasp.NRELAXING));
                        //ME20190301 END
                      }
                      if(xvasp.NRELAXING==xvasp.NRELAX) {
                        Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                        if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAXATION=]");return Krun;}
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RELAXATION=]");return Krun;} //CO20201111
                        KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAXING,false,FileMESSAGE);  //ME20190610 - or else SPIN_REMOVE_RELAX_2 won't work  //CO20210315 - never write_incar here (last step for sure)
                      }
                      KBIN::XVASP_INCAR_ADJUST_ICHARG(xvasp, vflags, aflags, xvasp.NRELAXING, (xvasp.NRELAXING<xvasp.NRELAX), FileMESSAGE);  //ME20191028 //CO20210315 - only write_incar if it's not the last relaxation (last step)
                    }
                    xvasp.NRELAXING=xvasp.NRELAX;
                    xvasp.NRELAXING++;
                    aus << 11111*xvasp.NRELAXING << "  END        - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 

                    // if((vflags.KBIN_VASP_FORCE_OPTION_LDAU0.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) && vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry) {
                    //   aus << 11111*xvasp.NRELAXING << "  EXTRA vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF" << endl;
                    //   aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); 
                    // }
                  }
                }
              }
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // RELAX_STATIC_BANDS RELAX_STATIC_BANDS RELAX_STATIC_BANDS REPEAT_STATIC REPEAT_STATIC_BANDS REPEAT_BANDS
              // STATIC_BANDS STATIC_BANDS STATIC_BANDS 	
              // RELAX_STATIC RELAX_STATIC RELAX_STATIC
              // REPEAT_STATIC REPEAT_STATIC REPEAT_STATIC
              // REPEAT_STATIC_BANDS REPEAT_STATIC_BANDS
              // REPEAT_BANDS REPEAT_BANDS REPEAT_BANDS
              if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {    // xvasp.RELAX>0 //CO20210315
                vector<double> xvasp_spin_evolution;
                xmatrix<double> rlattice(xvasp.str.lattice);

                string STRING_TO_SHOW="";
                if(vflags.KBIN_VASP_RUN.flag("STATIC")) STRING_TO_SHOW="STATIC";
                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) STRING_TO_SHOW="RELAX_STATIC";
                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) STRING_TO_SHOW="RELAX_STATIC_BANDS";
                if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) STRING_TO_SHOW="STATIC_BANDS";
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) STRING_TO_SHOW="REPEAT_STATIC"; //CO20210315
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) STRING_TO_SHOW="REPEAT_STATIC_BANDS";
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) STRING_TO_SHOW="REPEAT_BANDS";
                aus << "00000  MESSAGE MODE= (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    

                //[CO20210315 - wrong placement, change this with relaxations]xvasp.aopts.flag("FLAG::POSCAR_PRESERVED",TRUE); // in case of errors it is not lost but recycled
                //[CO20210315 - wrong placement, change this with relaxations]xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED",TRUE); // in case of errors it is not lost but recycled

                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
                  // DO THE RELAX PART (IF ANY)
                  KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                  if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC")) {
                    for(xvasp.NRELAXING=1;xvasp.NRELAXING<=xvasp.NRELAX;xvasp.NRELAXING++) {
                      aus << 11111*xvasp.NRELAXING << "  RELAXATION (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                      if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int>0) KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp,xvasp.NRELAXING); // ADIABATIC
                      if(xvasp.NRELAXING<xvasp.NRELAX)  {
                        Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                        if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS RELAXATION<]");return Krun;}
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RELAX_STATIC_BANDS RELAXATION<]");return Krun;} //CO20201111
                        //ME20190301 BEGIN
                        // CHGCAR/WAVECAR needs to be recycled if CHGCAR/WAVECAR=ON or VASP
                        // won't be able to read the files. Bug found by Rico Friedrich
                        if(vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option) {
                          //ME20191031 - Only recycle when ICHARG was found
                          string extra_incar = xvasp.AVASP_EXTRA_INCAR.str();
                          int nlines = aurostd::GetNLinesString(extra_incar);
                          int l;
                          string line;
                          for (l = 1; l <= nlines; l++) {
                            line = aurostd::RemoveWhiteSpaces(aurostd::GetLineString(extra_incar, l));
                            if (aurostd::substring2bool(line, "ICHARG=1", true)) break;
                          }
                          if (l <= nlines) KBIN::VASP_RecycleExtraFile(xvasp, "CHGCAR", "relax"+aurostd::utype2string<int>(xvasp.NRELAXING));
                        }
                        //CO+MM20211217
                        //In general it's not a good idea to use the WAVECAR from relax1 to start relax2
                        //a big change in the unit cell size/shape is going drastically change the plane waves that should be used
                        //we can add another option to do this in the future
                        //[CO+MM20211217]if(vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option) KBIN::VASP_RecycleExtraFile(xvasp, "WAVECAR", "relax"+aurostd::utype2string<int>(xvasp.NRELAXING));
                        //ME20190301 END
                      }
                      if(xvasp.NRELAXING==xvasp.NRELAX) {
                        Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(xvasp.NRELAXING),TRUE,FileMESSAGE);
                        if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS RELAXATION=]");return Krun;}
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                        //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RELAX_STATIC_BANDS RELAXATION=]");return Krun;} //CO20201111
                      }
                      KBIN::XVASP_INCAR_ADJUST_ICHARG(xvasp, vflags, aflags, xvasp.NRELAXING, true, FileMESSAGE);  //ME20191028 //CO20210315 - always write_incar, there's a STATIC that follows (at least)
                      xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
                      aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size()-1) << Message(__AFLOW_FILE__,aflags) << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                      //[CO20210315 - made robust enough to work in all cases]if(xvasp.NRELAXING<xvasp.NRELAX) 
                      KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAXING,true,FileMESSAGE); 	// check if it is the case of turning off spin  //CO20210315 - always write_incar, there's a STATIC that follows (at least)
                    }
                    if(xvasp.NRELAX>0) KBIN::VASP_Recycle(xvasp,"relax"+aurostd::utype2string(xvasp.NRELAX));  // bring back the stuff
                    //[CO20210315 - not sure why only if NRELAX==2, we could have NRELAX==1 and this still might apply]if(xvasp.NRELAX==2) 
                    KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp,aflags,vflags,xvasp.NRELAX,true,FileMESSAGE); 	// check if it is the case of turning off spin  //CO20210315 - always write_incar, there's a STATIC that follows (at least) //CO20210315 - no longer necessary per above (we check after every relaxation, but it doesn't hurt
                  }
                  if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) { //CO20210315 - looks like a typo to me //RELAX_STATIC
                    aus << "00000  NO RELAXATION IN (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                    xvasp.NRELAX=0; //CO20210315 - mimicking STATIC below
                  }
                  if(vflags.KBIN_VASP_RUN.flag("STATIC")) {
                    aus << "00000  NO RELAXATION IN (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    xvasp.NRELAX=0;
                  }
                  xvasp.NRELAXING=xvasp.NRELAX;
                } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_RUN.flag("STATIC")
                // REPEAT_STATIC ----------------------------------------------------------------------------
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")) {  //CO20210315
                  // LOAD FORMER LOCK
                  if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_STATIC"))) {
                    stringstream lock_recycled;
                    aurostd::file2stringstream(xvasp.Directory+"/REPEAT_STATIC",lock_recycled);
                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                    aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    //	aus << lock_recycled.str();
                    aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  }
                  // clean up worthless stuff
                  aurostd::execute("cd "+xvasp.Directory+" && rm -f *static* "+DEFAULT_AFLOW_END_OUT+"*");  //CO20210315 - delete aflow.end.out for --monitor_vasp
                  // UNZIP EVERYTHING
                  for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed 
                    // aurostd::execute("cd "+xvasp.Directory+" && "+"bzip2 -dfq *bz2 "); // ORIGINAL
                    aus << "00000  MESSAGE attempting UNZIP=" << XHOST.vzip[iext] << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                    aurostd::execute(XHOST.vzip[iext]+" -dfq `find \""+aurostd::CleanFileName(xvasp.Directory)+"\" -name \"*"+XHOST.vext[iext]+"\"`");
                  }		
                  if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.relax2"))) {
                    KBIN::VASP_Recycle(xvasp,"relax2");
                  } else {
                    if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.relax1"))) {
                      KBIN::VASP_Recycle(xvasp,"relax1");
                    } else {
                      aus << "REPEAT_STATIC: RELAX2 or RELAX1 must be present.";
                      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
                    }
                  }
                  //[CO20210315 - do before unzip]// clean up worthless stuff
                  //[CO20210315 - do before unzip]aurostd::execute("cd "+xvasp.Directory+" && rm -f *static*");
                } // vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC")
                // REPEAT_STATIC_BANDS PART ----------------------------------------------------------------------------
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {  //CO20210315
                  // LOAD FORMER LOCK
                  if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_STATIC_BANDS"))) {
                    stringstream lock_recycled;
                    aurostd::file2stringstream(xvasp.Directory+"/REPEAT_STATIC_BANDS",lock_recycled);
                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                    aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    //	aus << lock_recycled.str();
                    aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  }
                  // clean up worthless stuff
                  aurostd::execute("cd "+xvasp.Directory+" && rm -f *static* *bands* "+DEFAULT_AFLOW_END_OUT+"*");  //CO20210315 - delete aflow.end.out for --monitor_vasp
                  // UNZIP EVERYTHING
                  for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed 
                    // aurostd::execute("cd "+xvasp.Directory+" && "+"bzip2 -dfq *bz2 "); // ORIGINAL
                    aus << "00000  MESSAGE attempting UNZIP=" << XHOST.vzip[iext] << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                    aurostd::execute(XHOST.vzip[iext]+" -dfq `find \""+aurostd::CleanFileName(xvasp.Directory)+"\" -name \"*"+XHOST.vext[iext]+"\"`");
                  }		
                  if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.relax2"))) {
                    KBIN::VASP_Recycle(xvasp,"relax2");
                  } else {
                    if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.relax1"))) {
                      KBIN::VASP_Recycle(xvasp,"relax1");
                    } else {
                      aus << "REPEAT_STATIC_BANDS: RELAX2 or RELAX1 must be present.";
                      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
                    }
                  }
                  //[CO20210315 - do before unzip]// clean up worthless stuff
                  //[CO20210315 - do before unzip]aurostd::execute("cd "+xvasp.Directory+" && rm -f *static* *bands*");
                } // vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")
                // REPEAT_BANDS PART ----------------------------------------------------------------------------
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
                  // LOAD FORMER LOCK
                  if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_BANDS"))) {
                    stringstream lock_recycled;
                    aurostd::file2stringstream(xvasp.Directory+"/REPEAT_BANDS",lock_recycled);
                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                    aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    //	aus << lock_recycled.str();
                    aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  }
                  // clean up worthless stuff
                  aurostd::execute("cd "+xvasp.Directory+" && rm -f *bands* "+DEFAULT_AFLOW_END_OUT+"*"); //CO20210315 - delete aflow.end.out for --monitor_vasp
                  // UNZIP EVERYTHING
                  for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed 
                    // aurostd::execute("cd "+xvasp.Directory+" && "+"bzip2 -dfq *bz2 "); // ORIGINAL
                    aus << "00000  MESSAGE attempting UNZIP=" << XHOST.vzip[iext] << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                    aurostd::execute(XHOST.vzip[iext]+" -dfq `find \""+aurostd::CleanFileName(xvasp.Directory)+"\" -name \"*"+XHOST.vext[iext]+"\"`");
                  }		

                  if(aurostd::FileExist(xvasp.Directory+string("/POSCAR.static"))) {
                    KBIN::VASP_Recycle(xvasp,"static");
                  } else {
                    aus << "REPEAT_BANDS: STATIC must be present.";
                    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
                  }
                  //[CO20210315 - do before unzip]// clean up worthless stuff
                  //[CO20210315 - do before unzip]aurostd::execute("cd "+xvasp.Directory+" && rm -f *bands* ");
                } // vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
                // STATIC PART ----------------------------------------------------------------------------
                // STATIC PART ----------------------------------------------------------------------------
                // STATIC PART ----------------------------------------------------------------------------
                // STATIC PART ----------------------------------------------------------------------------
                // NOW DO THE STATIC PATCHING POSCAR
                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
                  aus << "00000  MESSAGE Patching POSCAR " << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  vflags.KBIN_VASP_RUN.push("STATIC");        // force to suck them up from STATIC_KPPRA....
                  // LOAD THE RELAXED STRUCTURE WHICH WILL BE USED FOR THE DIRECTIONS
                  stringstream straus;
                  aurostd::file2stringstream(xvasp.Directory+"/POSCAR",straus);
                  xvasp.str=xstructure(straus,IOVASP_AUTO);
                  xvasp.str.FixLattices();
                  rlattice=xvasp.str.lattice; // in rlattice I`ve always the final structure
                  bool STATIC_DEBUG=FALSE;//TRUE;
                  // RECREATE CONVENTIONAL OR PRIMITIVE
                  if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
                    // shall we restandardize ?
                  }
                  //WSETYAWAN_LOOK
                  //    STATIC_DEBUG=TRUE;
                  if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_RUN.flag(\"RELAX_STATIC_BANDS\")=" << vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_RUN.flag(\"STATIC_BANDS\")=" << vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_STATIC_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_REPEAT.flag(\"REPEAT_BANDS\")=" << vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS") << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.isentry << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;}
                    //[CO20210805 - OBSOLETE]if(STATIC_DEBUG) {aus << "STATIC_DEBUG: vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG=" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG << endl;}
                    if(STATIC_DEBUG) {aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    string stringBZ="";
                    bool foundBZ=FALSE;
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << xvasp.str << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << rlattice << endl;}
                    //[CO20210805 - OBSOLETE]if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG==FALSE)
                    if(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string=="AUTO") { //CO20210805
                      foundBZ=FALSE;
                      //AS20200724 BEGIN
                      // when KBIN_VASP_KPOINTS_BANDS_LATTICE=AUTO we need to retrieve
                      // the lattice type. For example, if CONVERT_UNIT_CELL=PRES the
                      // call to KPOINTS_Directions() will lead to an error since 
                      // vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string is empty
                      xvasp.str.GetLatticeType();
                      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();
                      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(xvasp.str.bravais_lattice_variation_type);
                      //AS20200724 END
                    }
                    if(0){  //CO20210805 - no need to run this code twice, run once after you standardize the structure
                      stringBZ=LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string,rlattice,vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int,xvasp.str.iomode,foundBZ); // rlattice = updated structure
                      if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                      if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << stringBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                      if(STATIC_DEBUG) {aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                      if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    }

                    // always recalculate standardization
                    if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE")==TRUE) {
                      // nothing
                      aus << "00000  MESSAGE PRESERVING ORIGINAL STRUCTURE" << Message(__AFLOW_FILE__,aflags) << endl;
                      aus << "00000  MESSAGE ORIGINAL: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                      aus << "00000  MESSAGE ORIGINAL: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    }else{
                      aus << "00000  MESSAGE WARNING RECALCULATING STANDARD STRUCTURE" << Message(__AFLOW_FILE__,aflags) << endl;
                      aus << "00000  MESSAGE BEFORE: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                      aus << "00000  MESSAGE BEFORE: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                      // reshuffle the structure
                      if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {xvasp.str.Standard_Conventional_UnitCellForm();}
                      else {xvasp.str.Standard_Primitive_UnitCellForm();}
                      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
                      xvasp.POSCAR << xvasp.str;
                      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
                      aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR"));
                      xvasp.str.FixLattices();
                      rlattice=xvasp.str.lattice; // in rlattice I`ve always the final structure
                      //CO20210805 - both Standard_Conventional_UnitCellForm() and Standard_Primitive_UnitCellForm() 
                      //return back updated bravais_lattice_variation_type, no need to recalculate with GetLatticeType()
                      vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.clear();vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.push(xvasp.str.bravais_lattice_variation_type); //WSETYAWAN mod
                      aus << "00000  MESSAGE AFTER: a,b,c,alpha,beta,gamma " << xvasp.str.a << "," << xvasp.str.b << "," << xvasp.str.c << "," << xvasp.str.alpha << "," << xvasp.str.beta << "," << xvasp.str.gamma << endl;
                      aus << "00000  MESSAGE AFTER: lattice: " << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    }
                    stringBZ=LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string,rlattice,vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int,xvasp.str.iomode,foundBZ); // rlattice = updated structure
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << stringBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: foundBZ=" << foundBZ << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << xvasp.str << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << rlattice << endl;}
                    if(STATIC_DEBUG) {aus << "STATIC_DEBUG: " << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
                    if(foundBZ==FALSE) {
                      aus << "Unrecoverable error, lattice not found:" << std::endl;
                      aus << xvasp.str;
                      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
                    }
                  }
                  // done with the fixing
                  xvasp.str.FixLattices();
                  rlattice=xvasp.str.lattice; // in rlattice I`ve always the final structure
                  xvasp.aopts.flag("FLAG::POSCAR_PRESERVED",TRUE); //CO20210315 - correct placement // in case of errors it is not lost but recycled
                  // NOW DO THE STATIC PATCHING KPOINTS
                  aus << "00000  MESSAGE Patching KPOINTS " << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  //
                  // [OBSOLETE]	      KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
                  KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags);
                  KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags);
                  aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
                  // NOW DO THE STATIC PATCHING INCAR
                  aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] Patching INCAR (static_patching)" << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                  // KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags); // BETTER than produce, SHOULD reread it
                  KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                  // KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags);  // MODIFY ACCORDINGLY
                  KBIN::XVASP_INCAR_Relax_Static_ON(xvasp,vflags);     // FIX
                  // do the RWIGS ON
                  if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC) KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,ON);
                  // done write INCAR
                  aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                  // NOW DO THE STATIC RUN
                  if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS")) xvasp.NRELAXING=xvasp.NRELAX; //0;
                  if(vflags.KBIN_VASP_RUN.flag("STATIC")) xvasp.NRELAXING=xvasp.NRELAX; // 0;
                  xvasp.NRELAXING++;
                  aus << aurostd::PaddedPRE(aurostd::utype2string(11111*xvasp.NRELAXING),5,"0") << "  STATIC (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory 
                    << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                  if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS STATIC]");return Krun;}
                  //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("static"));
                  Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                  if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RELAX_STATIC_BANDS STATIC]");return Krun;} //CO20201111  //AFTER CONTCAR_SAVE_
                  bool qmwrite=TRUE;
                  KBIN::VASP_Backup(xvasp,qmwrite,string("static"));
                  xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
                  aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size()-1) << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC")
                // BANDS PART ----------------------------------------------------------------------------
                // BANDS PART ----------------------------------------------------------------------------
                // BANDS PART ----------------------------------------------------------------------------
                // BANDS PART ----------------------------------------------------------------------------
                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")) {
                  // NOW DO THE BANDS PATCHING KPOINTS (if necessary...)
                  KBIN::VASP_Recycle(xvasp,"static");  // bring back the stuff
                  KBIN::VASP_RecycleExtraFile(xvasp,"CHGCAR","static");  // bring back the stuff
                  xvasp.aopts.flag("FLAG::POSCAR_PRESERVED",TRUE); //CO20210315 - correct placement // in case of errors it is not lost but recycled
                  xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED",TRUE); // in case of errors it is not lost but recycled
                  aus << "00000  MESSAGE Patching KPOINTS with BANDS LATTICE = \"" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << "\"" << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  // KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags);
                  // KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags);
                  // poscar was already conventionalized in the static part
                  xvasp.KPOINTS.clear();xvasp.KPOINTS.str(std::string());
                  //	      xvasp.KPOINTS <<
                  string stringBZ="";
                  bool foundBZ=false;
                  stringBZ=LATTICE::KPOINTS_Directions(vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string,rlattice,vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int,xvasp.str.iomode,foundBZ); // rlattice = updated structure
                  if(foundBZ==FALSE) {  //CO20210805
                    aus << "Unrecoverable error, lattice not found:" << std::endl;
                    aus << xvasp.str;
                    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
                  }
                  // removed stuff BELOW
                  xvasp.KPOINTS << stringBZ;
                  aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
                  xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED",TRUE); // don't touch kpoints if there are flaws
                  // NOW DO THE BANDS PATCHING INCAR
                  aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] Patching INCAR (bands_patching)" << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                  //  KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileAFLOWIN,FileMESSAGE,aflags,kflags,vflags); // BETTER than produce, SHOULD reread it
                  KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                  // KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags); // MODIFY ACCORDINGLY
                  KBIN::XVASP_INCAR_Bands_ON(xvasp,vflags);     // FIX
                  // do the RWIGS OFF
                  if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC)
                    KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,OFF);
                  // done write INCAR
                  aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                  if(0)  {
                    stringstream command;
                    command << "cd \"" <<  xvasp.Directory << "\"" << endl;
                    command << "cat INCAR | grep -v NGXF | grep -v NGYF | grep -v NGZF > INCAR.new" << endl;
                    command << "cat OUTCAR.static | grep NGXF | grep dimension | sed \"s/NG/\nNG/g\" | grep -v dimension | sed \"s/ //g\" >> INCAR.new" << endl;
                    command << "mv INCAR.new INCAR " << endl;
                    aurostd::execute(command);
                  }
                  // NOW DO THE BANDS RUN
                  xvasp.NRELAXING++;
                  aus << 11111*xvasp.NRELAXING << "  BANDS (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - K=[" << vflags.KBIN_VASP_KPOINTS_BANDS_LATTICE.content_string << "," << vflags.KBIN_VASP_KPOINTS_BANDS_GRID.content_int << "]" << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                  //[CO20210315 - OBSOLETE]vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("ROTMAT");	// dont mess up KPOINTS in bands
                  //[CO20210315 - OBSOLETE]vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("IBZKPT");    // dont mess up KPOINTS in bands
                  //[CO20210315 - OBSOLETE]vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("EDDRMM");	// dont mess up KPOINTS in bands
                  //[CO20210315 - OBSOLETE as above: xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED",TRUE)]vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.push("FIX:KPOINTS");	//CO20210315 - this should also be obsolete, Afix knows not to fix a non-automesh, also covered by xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED",TRUE) above // dont mess up KPOINTS in bands
                  Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                  if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RELAX_STATIC_BANDS BANDS]");return Krun;}
                  //  if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("bands"));
                  Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                  if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RELAX_STATIC_BANDS BANDS]");return Krun;} //CO20201111 //AFTER CONTCAR_SAVE_
                  bool qmwrite=FALSE;
                  KBIN::VASP_Backup(xvasp,qmwrite,string("bands"));
                  xvasp_spin_evolution.push_back(xvasp.str.qm_mag_atom); // keep track of spins
                  aus << "00000  MESSAGE RESULT SPIN=" << xvasp_spin_evolution.at(xvasp_spin_evolution.size()-1) << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
                // DIELECTRIC PART - DSCF ----------------------------------------------------------------
                // DIELECTRIC PART -----------------------------------------------------------------------
                // DIELECTRIC PART -----------------------------------------------------------------------
                // DIELECTRIC PART -----------------------------------------------------------------------
                xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED",FALSE); // bands are done... I can refix the KPOINTS
                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC")) {
                  // have static
                  if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) {  // check for DIELECTRIC STATIC
                    // check VASP version
                    string sversion=KBIN::OUTCAR2VASPVersionNumber(xvasp.Directory+"/OUTCAR.static"); //CO20210315
                    double dversion=KBIN::VASPVersionString2Double(sversion); //CO20210315
                    //[CO20210315 - not needed]vector<string> tokens; aurostd::string2tokensAdd(sversion,tokens,".");
                    //[CO20210315 - not needed]if(tokens.size()>0) dversion+=aurostd::string2utype<double>(tokens.at(0));
                    //[CO20210315 - not needed]if(tokens.size()>1) dversion+=aurostd::string2utype<double>(tokens.at(1))/10.0;
                    aus << "00000  MESSAGE Found VASP version=" << sversion << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    if(dversion<5.2) { // cant do it
                      aus << "EEEEE  ERROR: Dielectric calculations need VASP >=5.2" << Message(__AFLOW_FILE__,aflags) << endl;
                      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                      Krun=FALSE;return Krun;
                    }
                    // PROCEED
                    xvasp.NRELAXING++;
                    aus << 11111*xvasp.NRELAXING << "  RUN_DIELECTRIC_STATIC (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS") && kflags.KBIN_MPI_NCPUS==24) {
                      uint ncpus_before=kflags.KBIN_MPI_NCPUS;
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) kflags.KBIN_MPI_NCPUS=DUKE_MATERIALS_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      aflags.AFLOW_GLOBAL_NCPUS=-kflags.KBIN_MPI_NCPUS;
                      aus << "00000  MESSAGE Running RUN_DIELECTRIC_STATIC fixing mpivasp5 with " << ncpus_before << "-AMD cores to " << kflags.KBIN_MPI_NCPUS << "-AMD cores" << Message(__AFLOW_FILE__,aflags) << endl;
                    }	
                    if((aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH") || aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI") || aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB"))
                        && (kflags.KBIN_MPI_NCPUS==64 || kflags.KBIN_MPI_NCPUS==48 || kflags.KBIN_MPI_NCPUS==32)) {
                      uint ncpus_before=kflags.KBIN_MPI_NCPUS;
                      kflags.KBIN_MPI_NCPUS=DUKE_BETA_VASP5_CORES_DIELECTRIC; // something
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) kflags.KBIN_MPI_NCPUS=DUKE_BETA_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) kflags.KBIN_MPI_NCPUS=DUKE_BETA_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_X")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //CO20201220 X_X
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_CRAY")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //SD20221006 X_X
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_OLDCRAY")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //SD20221006 X_X
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_SMB")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //SD20221006 X_X
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::JHU_ROCKFISH")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //CO20220818 ROCKFISH
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //DX20190509 - MACHINE001
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //DX20190509 - MACHINE002
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE003")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //DX20201005 - MACHINE003
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE004")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //DX20211011 - MACHINE004
                      if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) kflags.KBIN_MPI_NCPUS=AFLOWLIB_VASP5_CORES_DIELECTRIC;  // bug in mpivasp5 //DX20190107 - CMU EULER
                      aflags.AFLOW_GLOBAL_NCPUS=-kflags.KBIN_MPI_NCPUS;
                      aus << "00000  MESSAGE Running RUN_DIELECTRIC_STATIC fixing mpivasp5 with " << ncpus_before << "-AMD cores to " << kflags.KBIN_MPI_NCPUS << "-AMD cores" << Message(__AFLOW_FILE__,aflags) << endl;
                    }	
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    // a. Reuse the INCAR.static  NELM = 0, and remove the NBANDS parameter.
                    // a. Set the k-point DKGRID < 0.10 (I've been using "aflow -k" for this).
                    // b. Retain the following static run entries and their values: ALGO, LREAL, NSIM, ISYM, IBRION, NSW, NELM, NELMIN, ENMAX, ISPIN, ISMEAR, SIGMA, and everything LDA+U related.
                    // c. Set NBANDS to a value that is around 10x the VASP default obtained in STEP 00.
                    // d. Eliminate PSTRESS, EMIN, EMAX, LORBIT, ISIF, NEDOS.
                    aus << "00000  MESSAGE Running RUN_DIELECTRIC_STATIC recycling static" << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    KBIN::VASP_Recycle(xvasp,"static");  // bring back the stuff
                    KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                    KBIN::VASP_Reread_KPOINTS(xvasp,FileMESSAGE,aflags); // REREAD IT
                    KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET(xvasp,kflags,vflags,"STATIC");     // FIX
                    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
                    aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
                    aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RUN_DIELECTRIC_STATIC]");return Krun;}
                    Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RUN_DIELECTRIC_STATIC]");return Krun;} //CO20201111
                    xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED",TRUE); // WAVECAR.dielectric_static
                    bool qmwrite=TRUE;
                    KBIN::VASP_Backup(xvasp,qmwrite,string("dielectric_static"));
                  }
                  if(vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC") && vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) {  // check for DIELECTRIC DYNAMIC
                    xvasp.NRELAXING++;
                    aus << 11111*xvasp.NRELAXING << "  RUN_DIELECTRIC_DYNAMIC (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    // a. Reuse the STEP 01 WAVECAR + ALGO=EXACT  NELM=1 LOPTICS=.TRUE. CSHIFT=0.15 OMEGAMAX=25 NEDOS=12500  Remove LEPSILON and LRPA
                    aus << "00000  MESSAGE Running RUN_DIELECTRIC_DYNAMIC recycling dielectric_static" << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                    KBIN::VASP_Recycle(xvasp,"dielectric_static"); // bring back the stuff
                    KBIN::VASP_RecycleExtraFile(xvasp,"WAVECAR","dielectric_static");  // bring back the stuff
                    KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                    KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET(xvasp,kflags,vflags,"DYNAMIC");   // FIX
                    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
                    aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
                    Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [RUN_DIELECTRIC_DYNAMIC]");return Krun;}
                    Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                    if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [RUN_DIELECTRIC_DYNAMIC]");return Krun;} //CO20201111
                    aurostd::execute("rm -f "+xvasp.Directory+"/WAVECAR.dielectric_static");
                    xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED",FALSE); // all gone
                    bool qmwrite=TRUE;
                    KBIN::VASP_Backup(xvasp,qmwrite,string("dielectric_dynamic"));
                  }
                  if(vflags.KBIN_VASP_RUN.flag("DSCF") && vflags.KBIN_VASP_RUN.flag("DIELECTRIC_DYNAMIC") && vflags.KBIN_VASP_RUN.flag("DIELECTRIC_STATIC")) {  // check for DIELECTRIC DYNAMIC
                    string message = "Dielectric calculations not implemented.";
                    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
                  }
                }
                // FINISHED
                xvasp.NRELAXING++;
                aus << 11111*xvasp.NRELAXING << "  END (" << STRING_TO_SHOW << ")        - " <<  xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);  // put back the options

                // CLEAN-UP BY WSETYAWAN
                ostringstream xaus;
                xaus << "cd " << xvasp.Directory << endl;
                // if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) xaus << "rm -f CHG.relax* CHGCAR.relax* POTCAR.* CHGCAR.bands CHG.bands" << endl;
                // if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS"))  xaus << "rm -f POTCAR.* CHGCAR.bands CHG.bands" << endl;
                if(vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS")) xaus << "rm -f CHG.relax* CHGCAR.relax* CHGCAR.bands CHG.bands" << endl;
                if(vflags.KBIN_VASP_RUN.flag("STATIC_BANDS"))  xaus << "rm -f CHGCAR.bands CHG.bands" << endl;
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS"))  xaus << "rm -f CHGCAR.bands CHG.bands" << endl;
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS")) {;}
                aurostd::execute(xaus);
                // done ....
                //WSETYAWAN, you might ask something more here
              } // vflags.KBIN_VASP_RUN.flag("RELAX_STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("STATIC_BANDS") || vflags.KBIN_VASP_RUN.flag("RELAX_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_STATIC_BANDS") || vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS")
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // --------------------------------------------------------------------------------------------------------------------
              // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL 	
              // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL 	
              // REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL REPEAT_DELSOL
              // PRL 105, 196403 (2010)
              if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) {
                bool delsol_d=FALSE;
                float NELECT=0,Nr=0;
                xmatrix<double> rlattice(3,3);
                string STRING_TO_SHOW="",stmp="";
                string fnamedelsol=xvasp.Directory+string("/delsol.tmp");
                stringstream command,strdelsol;
                ifstream fdelsol;
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")) STRING_TO_SHOW="REPEAT_DELSOL";
                aus << "00000  MESSAGE MODE= (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                // LOAD FORMER LOCK
                if(aurostd::FileExist(xvasp.Directory+string("/REPEAT_DELSOL"))) {
                  stringstream lock_recycled;
                  aurostd::file2stringstream(xvasp.Directory+"/REPEAT_DELSOL",lock_recycled);
                  aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                  aus << "XXXXX FORMER LOCK BEGIN, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                  aus << lock_recycled.str();
                  aus << "XXXXX FORMER LOCK END, recycled (" << STRING_TO_SHOW << ") - " << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                  aus << "XXXXX ---------------------------------------------------------------------------------------------- " << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                }
                // UNZIP EVERYTHING	      
                for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed 
                  // command.clear();command.str(std::string());  // ORIGINAL
                  // command << "cd " <<  xvasp.Directory << endl; command << "bzip2" << " -dfq *bz2 " << endl;  // ORIGINAL
                  // aurostd::execute(command);  // ORIGINAL
                  aus << "00000  MESSAGE attempting UNZIP=" << XHOST.vzip[iext] << endl; aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);		  
                  aurostd::execute(XHOST.vzip[iext]+" -dfq `find \""+aurostd::CleanFileName(xvasp.Directory)+"\" -name \"*"+XHOST.vext[iext]+"\"`");
                }		

                // copy INCAR, POSCAR, KPOINTS, POTCAR from *.static
                KBIN::VASP_Recycle(xvasp,"static");  // bring back the stuff
                //Scanning whether it is sp or spd from the POTCAR
                command.clear();command.str(std::string());
                command << "cd " <<  xvasp.Directory << endl;
                command << "grep VRHFIN POTCAR.static | sed \'s/:/\\n/g\' | grep -v VRHFIN > delsol.tmp" << endl;
                aurostd::execute(command);	
                strdelsol.clear();strdelsol.str(std::string());
                aurostd::file2stringstream(xvasp.Directory+"/delsol.tmp",strdelsol);
                command.clear();command.str(std::string());
                command << "rm -f " << xvasp.Directory << "/delsol.tmp" << endl;
                aurostd::execute(command);
                delsol_d=FALSE;
                if((aurostd::substring2bool(strdelsol.str(),"d"))) delsol_d=TRUE;
                //Scanning NELECT from OUTCAR.static
                command.clear();command.str(std::string());
                command << "cd " << xvasp.Directory << endl;
                command << "grep NELECT OUTCAR.static | sed \'s/=/\\n/g\' | grep -v NELECT > delsol.tmp" << endl;
                aurostd::execute(command);
                strdelsol.clear();strdelsol.str(std::string());
                aurostd::file2stringstream(xvasp.Directory+"/delsol.tmp",strdelsol);
                command.clear();command.str(std::string());
                command << "rm -f " << xvasp.Directory << "/delsol.tmp" << endl;
                aurostd::execute(command);
                strdelsol>>NELECT;
                //if(NELECT<1.0) ;//need to add error handling here
                Nr=NELECT/68.0;
                if(delsol_d) Nr=NELECT/72.0;
                aus << "DELSOL: NELECT N0=" << NELECT << endl << "DELSOL: NELECT n=" << Nr << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    

                // NOW MODIFY THE INCAR
                aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] modifying INCAR (delsol_patching)" << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                stmp="NELECT="+aurostd::utype2string(NELECT+Nr);
                xvasp.INCAR << aurostd::PaddedPOST(stmp,_incarpad_) << " # NELECT = N0 + n for DELSOL plus" << endl;
                // do the RWIGS OFF
                if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC)
                  KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,OFF);
                aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));

                ////Reread POSCAR
                //aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] rereading POSCAR (delsol)" << Message(__AFLOW_FILE__,aflags) << endl;
                //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                //KBIN::VASP_Reread_POSCAR(xvasp,FileMESSAGE,aflags); // REREAD IT
                ////Reread KPOINTS
                //aus << "00000  MESSAGE [" << STRING_TO_SHOW << "] rereading KPOINTS (delsol)" << Message(__AFLOW_FILE__,aflags) << endl;
                //aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                //KBIN::VASP_Reread_KPOINTS(xvasp,FileMESSAGE,aflags); // REREAD IT

                KBIN::VASP_RecycleExtraFile(xvasp,"POSCAR","static");  // bring back the stuff
                KBIN::VASP_RecycleExtraFile(xvasp,"KPOINTS","static");  // bring back the stuff

                // NOW RUN DELSOL plus
                uint vrelax=7;
                aus << 11111*vrelax << "  DELSOL plus (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << stmp << " "<< Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                vrelax++;
                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [DELSOL]");return Krun;}
                //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("dsolp"));
                Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [DELSOL]");return Krun;} //CO20201111 //AFTER CONTCAR_SAVE_
                bool qmwrite=FALSE;
                KBIN::VASP_Backup(xvasp,qmwrite,string("dsolp"));

                //NOW DO DELSOL minus
                KBIN::VASP_Recycle(xvasp,"static");
                KBIN::VASP_Reread_INCAR(xvasp,FileMESSAGE,aflags);
                stmp="NELECT="+aurostd::utype2string(NELECT-Nr);
                xvasp.INCAR << aurostd::PaddedPOST(stmp,_incarpad_) << " # NELECT = N0 - n for DELSOL minus" << endl;
                // do the RWIGS OFF
                if(vflags.KBIN_VASP_FORCE_OPTION_RWIGS_STATIC)
                  KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,FileMESSAGE,OFF);
                aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));

                KBIN::VASP_RecycleExtraFile(xvasp,"POSCAR","static");  // bring back the stuff
                KBIN::VASP_RecycleExtraFile(xvasp,"KPOINTS","static");  // bring back the stuff
                // NOW RUN DELSOL minus
                aus << 11111*vrelax << "  DELSOL minus (" << STRING_TO_SHOW << ") - " <<  xvasp.Directory << " - " << stmp << " "<< Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                vrelax++;
                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE);
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [DELSOL minus]");return Krun;}
                //	    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string("dsolm"));
                Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [DELSOL minus]");return Krun;} //CO20201111 //AFTER CONTCAR_SAVE_
                qmwrite=FALSE;
                KBIN::VASP_Backup(xvasp,qmwrite,string("dsolm"));		
                // FINISHED
                aus << 11111*vrelax << "  END (" << STRING_TO_SHOW << ")        - " <<  xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
                vflags.KBIN_VASP_RUN.flag("STATIC",FALSE);  // put back the options
                // CLEAN-UP BY WSETYAWAN
                ostringstream xaus;
                xaus << "cd " << xvasp.Directory << endl;
                if(vflags.KBIN_VASP_REPEAT.flag("REPEAT_BANDS"))  xaus << "rm -f CHGCAR.dsol* CHG.dsol*" << endl;
                aurostd::execute(xaus);
                // done ....
                //WSETYAWAN, you might ask something more here
              } // (vflags.KBIN_VASP_REPEAT.flag("REPEAT_DELSOL")
              //************* END OF REPEAT_DELSOL ************************
              //************* END OF REPEAT_DELSOL ************************
              //************* END OF REPEAT_DELSOL ************************
              //************* END OF REPEAT_DELSOL ************************
              //************* END OF REPEAT_DELSOL ************************
              //************* END OF REPEAT_DELSOL ************************
              // --------------------------------------------------------------------------------------------------------------------
              // KPOINTS KPOINTS KPOINTS
              if(vflags.KBIN_VASP_RUN.flag("KPOINTS")) {            // NON THREADED
                KBIN::VASP_Write_INPUT(xvasp,vflags); // VASP VASP WRITE
                xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);                                                                // BACKUP KPOINTS
                aurostd::stringstream2file(xvasp.KPOINTS_orig,string(xvasp.Directory+"/KPOINTS.orig"));   // BACKUP KPOINTS
                int kbak_k1=xvasp.str.kpoints_k1;
                int kbak_k2=xvasp.str.kpoints_k2;
                int kbak_k3=xvasp.str.kpoints_k3;
                //	    int kbak_kmax=xvasp.str.kpoints_kmax;   kbak_kmax=max(kbak_k1,kbak_k2,kbak_k3);
                int kk1,kk2,kk3;
                string relax,relaxfile;
                // 1st step: 1/2
                kk1=(kbak_k1+1)/2;kk2=(kbak_k2+1)/2;kk3=(kbak_k3+1)/2;
                relax="11111a ";relaxfile="relax1";
                aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                xvasp.str.kpoints_k1=min(kk1,kbak_k1);xvasp.str.kpoints_k2=min(kk2,kbak_k2);xvasp.str.kpoints_k3=min(kk3,kbak_k3);
                xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
                KBIN::XVASP_KPOINTS_numbers2string(xvasp);
                aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));   // BACKUP KPOINTS
                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,relaxfile,relaxfile,FALSE,FileMESSAGE);
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS 1]");return Krun;}
                //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [KPOINTS 1]");return Krun;} //CO20201111
                // 2nd step: 1/2
                kk1=3*(kbak_k1+1)/4;kk2=3*(kbak_k2+1)/4;kk3=3*(kbak_k3+1)/4;
                relax="11111b ";relaxfile="relax1";
                aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                xvasp.str.kpoints_k1=min(kk1,kbak_k1);xvasp.str.kpoints_k2=min(kk2,kbak_k2);xvasp.str.kpoints_k3=min(kk3,kbak_k3);
                xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
                KBIN::XVASP_KPOINTS_numbers2string(xvasp);
                aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));   // BACKUP KPOINTS
                Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,relaxfile,relaxfile,FALSE,FileMESSAGE);
                if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS 2]");return Krun;}
                //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [KPOINTS 2]");return Krun;} //CO20201111
                // 3rd step: 1/2
                kk1=kbak_k1;kk2=kbak_k2;kk3=kbak_k3;
                relax="11111c ";relaxfile="relax1";
                aus << relax << "RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3 << "]" << Message(__AFLOW_FILE__,aflags) << endl;
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                xvasp.str.kpoints_k1=min(kk1,kbak_k1);xvasp.str.kpoints_k2=min(kk2,kbak_k2);xvasp.str.kpoints_k3=min(kk3,kbak_k3);
                xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
                KBIN::XVASP_KPOINTS_numbers2string(xvasp);
                aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));   // BACKUP KPOINTS
                // ----------------------------------------------------------
                // with PAWGGA2
                if(PAWGGA2) {
                  // STEP 1
                  aus << "11111  RELAXATION - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                  Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax2paw_gga",FALSE,FileMESSAGE);
                  if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS PAWGGA2]");return Krun;}
                  //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                  //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [KPOINTS PAWGGA2]");return Krun;} //CO20201111
                  aus << "22222  END        - " <<  xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                }
                // ----------------------------------------------------------
                // norma, without PAWGGA2
                if(!PAWGGA2) {
                  int vrelax=1;
                  for(int i=1;i<=xvasp.NRELAX;i++) {
                    aus << 11111*vrelax << "   RELAXATION - " << xvasp.Directory << " - K=[" << kk1 << " " << kk2 << " " << kk3
                      << "]" << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                    if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int>0) KBIN::XVASP_INCAR_LDAU_ADIABATIC(xvasp,i); // ADIABATIC
                    if(i<xvasp.NRELAX)  {
                      Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(vrelax),"relax"+aurostd::utype2string(vrelax),TRUE,FileMESSAGE);
                      if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINT 4]");return Krun;}
                      //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                      //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [KPOINT 4]");return Krun;} //CO20201111
                    }
                    if(i==xvasp.NRELAX) {
                      Krun=KBIN::VASP_Run(xvasp,aflags,kflags,vflags,"relax"+aurostd::utype2string(vrelax),TRUE,FileMESSAGE);
                      if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error [KPOINTS 5]");return Krun;}
                      //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]Krun=KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true); //CO20201111
                      //[CO20210104 - OUTCAR has already been moved to OUTCAR.RELAX, check is inside VASP_RUN()]if(!Krun) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  runtime error (OUTCAR_INCOMPLETE) [KPOINTS 5]");return Krun;} //CO20201111
                    }
                    vrelax++;
                  }
                  aus << 11111*vrelax << "   END        - " << xvasp.Directory << " - " << kflags.KBIN_BIN << Message(__AFLOW_FILE__,aflags) << endl;
                  aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
                }
              } // KPOINTS KPOINTS KPOINTS
              // ***************************************************************************
              // POSTSCRIPT
              if(!vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT"))
                if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
                  KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
              // ***************************************************************************
            }
          }
        }
      }
      // **********
      // some verbose
      if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus << "00000  MESSAGE END loop " << xvasp.POSCAR_index+1 << "/" << vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING.size()  //CO20200624 - +1
          << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE END loop in directory =" << xvasp.Directory << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // compress the subdirectories
        if(Krun && kflags.KZIP_COMPRESS) Krun=(Krun && KBIN::CompressDirectory(aflags,kflags));
      }
      aflags=aflags_backup;kflags=kflags_backup; // RESTORE
    } // LOOP ixvasp
    // ***************************************************************************
    aflags=aflags_backup;kflags=kflags_backup; // RESTORE
    // POSTSCRIPT
    if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT"))
      if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
        KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
    // ***************************************************************************
    FileAFLOWIN.clear();FileAFLOWIN.close();
    return Krun;
  }
} // namespace

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

int CheckStringInFile(string FileIn,string str,int PID,int TID) { //CO20200502 - threadID
  int out;
  ostringstream aus_exec;
  aus_exec << "cat " << FileIn << " | grep -c \"" << str << "\" > aflow.out." << PID << "." << TID << endl; //CO20200502 - threadID
  aurostd::execute(aus_exec);
  ifstream FileCHECK;
  stringstream FineNameAflowTmpPIDTID;  //CO20200502 - threadID
  FineNameAflowTmpPIDTID << "aflow.out." << PID << "." << TID;  //CO20200502 - threadID
  FileCHECK.open(FineNameAflowTmpPIDTID.str().c_str(),std::ios::in);  //CO20200502 - threadID
  FileCHECK >> out;
  FileCHECK.clear();FileCHECK.close();
  aus_exec << "rm -f aflow.out." << PID << "." << TID << endl;  //CO20200502 - threadID
  aurostd::execute(aus_exec);
  return out;
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// FLAG Class for KBIN_VASP_Run

namespace KBIN {
  bool ReachedAccuracy2bool(const string& scheme,const aurostd::xoption& xRequiresAccuracy,const aurostd::xoption& xmessage,bool vasp_still_running){ //CO20210315
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " xRequiresAccuracy.flag(\"" << scheme << "\")=" << xRequiresAccuracy.flag(scheme) << endl;
      cerr << __AFLOW_FUNC__ << " vasp_still_running=" << vasp_still_running << endl;
      cerr << __AFLOW_FUNC__ << " xmessage.flag(\"REACHED_ACCURACY\")=" << xmessage.flag("REACHED_ACCURACY") << endl;
    }

    if(xRequiresAccuracy.flag(scheme)==false){return true;}  //this scheme does not require xmessage.flag("REACHED_ACCURACY"), return true for '&&' logic

    //CO20210601 - I am decoupling vasp_still_running and xmessage.flag("REACHED_ACCURACY") FOR RELAXATIONS ONLY
    //on qrats, I noticed that a calculation can reach accuracy, but not finish (incomplete OUTCAR), thus it is frozen
    //a good solution would be to take the CONTCAR as the new input and restart (RELAXATIONS ONLY)
    //to diagnose and treat this problem correctly, we need to consider vasp_still_running and xmessage.flag("REACHED_ACCURACY") independently (for this case only)
    //[CO20210602 - must use, otherwise NUM_PROB appears early]if(vasp_still_running==false){;}  //keep busy
    return (vasp_still_running==false && xmessage.flag("REACHED_ACCURACY")==false);  //vasp_still_running==false && - xmessage.flag("REACHED_ACCURACY") should already understand whether to consider vasp_still_running
  }
  void VASP_ProcessWarnings(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,aurostd::xoption& xmessage,aurostd::xoption& xwarning,ofstream &FileMESSAGE) { //CO20210315
    aurostd::xoption xmonitor;
    return VASP_ProcessWarnings(xvasp,aflags,kflags,xmessage,xwarning,xmonitor,FileMESSAGE);
  }
  void VASP_ProcessWarnings(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,aurostd::xoption& xmessage,aurostd::xoption& xwarning,aurostd::xoption& xmonitor,ofstream &FileMESSAGE) { //CO20210315
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    stringstream aus;
    bool VERBOSE=(FALSE || XHOST.vflag_control.flag("MONITOR_VASP")==false || LDEBUG);

    if(!aurostd::FileExist(xvasp.Directory+"/"+DEFAULT_VASP_OUT)){return;}
    if(!aurostd::FileExist(xvasp.Directory+"/INCAR")){return;}
    if(!aurostd::FileExist(aflags.Directory+"/"+_AFLOWLOCK_)){return;} //we needed it above to get the vasp_bin
    bool vasp_monitor_running=AFLOW_MONITOR_instance_running(aflags);

    long int tmod_outcar=aurostd::SecondsSinceFileModified(xvasp.Directory+"/"+"OUTCAR"); //better to look at OUTCAR than vasp.out, when vasp is killed you get errors in vasp.out, resetting the time
    unsigned long long int fsize_vaspout=aurostd::FileSize(xvasp.Directory+"/"+DEFAULT_VASP_OUT);
    if(VERBOSE){
      aus << "00000  MESSAGE time since " << "OUTCAR" << " last modified: " << tmod_outcar << " seconds (max=" << SECONDS_STALE_OUTCAR << " seconds)" << Message(__AFLOW_FILE__,aflags) << endl;
      if(LDEBUG){cerr << aus.str();}
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE size of " << DEFAULT_VASP_OUT << ": " << fsize_vaspout << " bytes (max=" << BYTES_MAX_VASP_OUT << " bytes)" << Message(__AFLOW_FILE__,aflags) << endl;
      if(LDEBUG){cerr << aus.str();}
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    //CO20210315 - reading the full vasp.out does not work
    //vasp can spit out so many warnings that the file can be >50GB, killing aflow's memory
    //the subprocess is really the fastest way to grep (regex support)
    //even better than reading the file line by line and processing string (toupper, remove white spaces, etc.)
    //CO20210315 - opening up subshells for grep (substring_present_file_FAST) is very expensive, especially many times
    //better to read file in once
    //[CO20210315 - does not work for big files]string content_vasp_out=aurostd::file2string(xvasp.Directory+"/"+DEFAULT_VASP_OUT);  //no "comments" to remove in this output file
    //[CO20210315 - does not work for big files]content_vasp_out=aurostd::RemoveWhiteSpaces(content_vasp_out);  //remove whitespaces
    //[CO20210315 - does not work for big files]content_vasp_out=aurostd::toupper(content_vasp_out);  //put toupper to eliminate case-sensitivity 

    //do memory check
    double usage_percentage_ram=0.0,usage_percentage_swap=0.0;
    bool approaching_oom=false;
    if(aurostd::GetMemoryUsagePercentage(usage_percentage_ram,usage_percentage_swap)){approaching_oom=(usage_percentage_ram>=MEMORY_MAX_USAGE_RAM && usage_percentage_swap>=MEMORY_MAX_USAGE_SWAP);}
    if(approaching_oom){  //might be a quick memory spike, try again
      if(VERBOSE){
        aus << "00000  MESSAGE ram memory used: " << aurostd::utype2string(usage_percentage_ram,2,FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_RAM << "%)" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE swap memory used: " << aurostd::utype2string(usage_percentage_swap,2,FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_SWAP << "%)" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "00000  MESSAGE reading memory again after " << SECONDS_SLEEP_VASP_MONITOR << " second sleep" << endl;
        if(LDEBUG){cerr << aus.str();}
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      aurostd::Sleep(SECONDS_SLEEP_VASP_MONITOR);
      if(aurostd::GetMemoryUsagePercentage(usage_percentage_ram,usage_percentage_swap)){approaching_oom=(usage_percentage_ram>=MEMORY_MAX_USAGE_RAM && usage_percentage_swap>=MEMORY_MAX_USAGE_SWAP);}
    }
    if(VERBOSE||approaching_oom){
      aus << "00000  MESSAGE ram memory used: " << aurostd::utype2string(usage_percentage_ram,2,FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_RAM << "%)" << Message(__AFLOW_FILE__,aflags) << endl;
      aus << "00000  MESSAGE swap memory used: " << aurostd::utype2string(usage_percentage_swap,2,FIXED_STREAM) << "% (max=" << MEMORY_MAX_USAGE_SWAP << "%)" << Message(__AFLOW_FILE__,aflags) << endl;
      if(LDEBUG){cerr << aus.str();}
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    //get INCAR
    VASP_Reread_INCAR(xvasp);  //preload incar

    if(LDEBUG){aus << __AFLOW_FUNC__ << " [1]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    //get whether relaxing or not
    bool relaxing=false;
    if(aurostd::kvpair2bool(xvasp.INCAR,"IBRION","=")){
      if(aurostd::kvpair2bool(xvasp.INCAR,"NSW","=") && aurostd::kvpair2utype<int>(xvasp.INCAR,"NSW","=")>0){relaxing=true;}
    }

    //need to get ISYM and ISPIN (ISPIND too)
    //get ISYM
    int isym_current=2; //CO20200624 - VASP default for non-USPP runs, add check later for USPP //https://www.vasp.at/wiki/index.php/ISYM
    if(aurostd::substring2bool(xvasp.INCAR,"ISYM=")){isym_current=aurostd::substring2utype<int>(xvasp.INCAR,"ISYM=");}
    //[CO20210315 - ISPIND not necessary]//get ISPIND
    //[CO20210315 - ISPIND not necessary]//seems ISPIND is not read from the INCAR: https://www.vasp.at/forum/viewtopic.php?f=3&t=3037
    //[CO20210315 - ISPIND not necessary]int ispind_current=1; //CO20200624 - VASP default //https://cms.mpi.univie.ac.at/vasp/guide/node87.html
    //[CO20210315 - ISPIND not necessary]if(aurostd::substring2bool(xvasp.INCAR,"ISPIND=")){ispind_current=aurostd::substring2utype<int>(xvasp.INCAR,"ISPIND=");}
    //get ISPIN
    int ispin_current=1; //CO20200624 - VASP default  //https://www.vasp.at/wiki/index.php/ISPIN
    if(aurostd::substring2bool(xvasp.INCAR,"ISPIN=")){ispin_current=aurostd::substring2utype<int>(xvasp.INCAR,"ISPIN=");}
    //[CO20210315 - ISPIND not necessary]if(ispin_current==2 && ispind_current==1){ispind_current=2;}  //in case ISPIND is not written but spin is on

    //there are issues getting the correct vasp binary since this is an entirely different aflow instance
    //we might run the other aflow instance with --mpi or --machine flags that affect which vasp binary we use
    //therefore, the most robust way to define the binary is to search the LOCK file
    //[CO20210315 - OBSOLETE]string& vasp_bin=kflags.KBIN_MPI_BIN;
    //[CO20210315 - OBSOLETE]if(!(kflags.KBIN_MPI==true||XHOST.MPI==true)){vasp_bin=kflags.KBIN_BIN;}
    string vasp_bin="",vasp_pgid="";
    GetVASPBinaryFromLOCK(xvasp.Directory,vasp_bin);
    vasp_bin=aurostd::basename(vasp_bin); //remove directory stuff
    vasp_pgid=aurostd::utype2string(getpgrp()); //SD20220406 - need PGID for VASP_instance_running
    if(vasp_bin.empty()){ //rely on defaults here in case we're not running --monitor_vasp
      vasp_bin=kflags.KBIN_MPI_BIN;
      if(!(kflags.KBIN_MPI==true||XHOST.MPI==true)){vasp_bin=kflags.KBIN_BIN;}
    }

    //CO20210315 - determine if vasp is still running
    //[SD20220406 - OBSOLETE]bool vasp_still_running=false; //CO20210315 - look at wording, very careful, this bool implies vasp WAS running, and now is not
    //[SD20220406 - OBSOLETE]if(XHOST.vflag_control.flag("KILL_VASP_ALL")){
    //[SD20220406 - OBSOLETE]  //KILL_VASP_ALL allows us to check for ANY instance of vasp running (assumes exclusivity in node environment)
    //[SD20220406 - OBSOLETE]  //we will develop more precise methods for tracking parent/child processes of aflow to target specific vasp instances in the future
    //[SD20220406 - OBSOLETE]  vasp_still_running=VASP_instance_running(vasp_bin);
    //[SD20220406 - OBSOLETE]}
    bool vasp_still_running=VASP_instance_running(vasp_bin,vasp_pgid); //SD20220406 - we are not using the --kill_vasp_all flag anymore
    if(LDEBUG){aus << __AFLOW_FUNC__ << " vasp_still_running=" << vasp_still_running << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    //vasp can spit out many warnings quickly, going from ~10MB to 1GB within 2 minutes
    //we will pass BYTES_MAX_VASP_OUT into substring_present_file_FAST as a stop condition
    //then OUTPUT_LARGE will be triggered at the bottom
    unsigned long long int grep_stop_condition=AUROSTD_MAX_ULLINT;
    if(vasp_monitor_running && vasp_still_running){grep_stop_condition=BYTES_MAX_VASP_OUT;} //no stop condition when vasp is not running, we must process the errors

    //CO20210315 - might consider "renice 20 -p VASP_PIDs" to give greps below priority

    bool renice=false;
    if(vasp_still_running && aurostd::ReniceAvailable() && fsize_vaspout>=10000000){renice=aurostd::ProcessRenice(vasp_bin,15);} //renice vasp if vasp.out>=10Mb so grep can work faster //CO20221028 - renice up/down is not always possible, now we check first

    uint i=0;
    vector<string> vtokens;

    //determine which schemes require xmessage.flag("REACHED_ACCURACY")
    //CO+AS202010315 - considering DENTET for this list, patches seem to be working. See: LIB3/CTa_pvTi_sv/ABC2_tP8_123_h_h_abg-001.ABC
    aurostd::xoption xRequiresAccuracy;
    aurostd::string2tokens("DAV,EDDRMM,IBZKPT,NUM_PROB,ZBRENT",vtokens,","); //DENTET,
    for(i=0;i<vtokens.size();i++){xRequiresAccuracy.flag(vtokens[i],true);}

    xwarning.clear(); //CO20210315 - very important to clear!

    bool vasp_oszicar_unconverged=KBIN::VASP_OSZICARUnconverged(xvasp.Directory+"/OSZICAR",xvasp.Directory+"/OUTCAR");

    bool RemoveWS=true; //remove whitespaces
    bool case_insensitive=true;
    bool expect_near_end=true;  //cat vs. tac

    //WARNINGS START
    if(LDEBUG){aus << __AFLOW_FUNC__ << " checking warnings" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    xwarning.flag("OUTCAR_INCOMPLETE",vasp_still_running==false && !KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,false));  //CO20201111
    //CO20210601 - I am decoupling vasp_still_running and xmessage.flag("REACHED_ACCURACY") FOR RELAXATIONS ONLY
    //on qrats, I noticed that a calculation can reach accuracy, but not finish (incomplete OUTCAR), thus it is frozen
    //a good solution would be to take the CONTCAR as the new input and restart (RELAXATIONS ONLY)
    //to diagnose and treat this problem correctly, we need to consider vasp_still_running and xmessage.flag("REACHED_ACCURACY") independently (for this case only)
    bool reached_accuracy_relaxing=(relaxing==true  && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"reached required accuracy",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );  //[CO20210601 - only check for static, not relax (might be frozen, and we can fix by appending CONTCAR)]vasp_still_running==false &&
    bool reached_accuracy_static  =(relaxing==false && vasp_still_running==false && vasp_oszicar_unconverged==false && xwarning.flag("OUTCAR_INCOMPLETE")==false); //CO20210315 - "reached accuracy" for static/bands calculation is a converged electronic scf, need to also check OUTCAR_INCOMPLETE, as a converged OSZICAR might actually be a run that ended because of an error
    xmessage.flag("REACHED_ACCURACY",(reached_accuracy_relaxing || reached_accuracy_static));

    if(LDEBUG){
      aus << __AFLOW_FUNC__ << " relaxing=" << relaxing << endl;
      aus << __AFLOW_FUNC__ << " vasp_oszicar_unconverged=" << vasp_oszicar_unconverged << endl;
      aus << __AFLOW_FUNC__ << " xwarning.flag(\"OUTCAR_INCOMPLETE\")=" << xwarning.flag("OUTCAR_INCOMPLETE") << endl;
      aus << __AFLOW_FUNC__ << " xmessage.flag(\"REACHED_ACCURACY\")=" << xmessage.flag("REACHED_ACCURACY") << endl;
      cerr << aus.str();
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    string scheme="";
    bool found_warning=false;
    //VASP's internal symmetry routines START
    //CO20200624 - these are all related to VASP's internal symmetry routines
    //they would all benefit from similar fixes (except NKXYZ_IKPTD which requires KPOINTS to be reduced)
    //SGRCON+NIRMAT: https://dannyvanpoucke.be/vasp-errors-en/
    //IBZKPT+KKSYM: http://www.error.wiki/VERY_BAD_NEWS!_internal_error_in_subroutine_IBZKPT
    //IBZKPT+KKSYM: https://www.vasp.at/forum/viewtopic.php?f=3&t=7811
    //INVGRP+SYMPREC: https://www.vasp.at/forum/viewtopic.php?t=486
    //INVGRP+SYMPREC: https://www.vasp.at/forum/viewtopic.php?f=3&t=9435&p=9473
    //
    scheme="SGRCON";  //usually goes with NIRMAT
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"VERY BAD NEWS! internal error in subroutine SGRCON",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="NIRMAT"; //usually goes with SGRCON
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"Found some non-integer element in rotation matrix",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="IBZKPT";  //usually goes with KKSYM or NKXYZ_IKPTD
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"VERY BAD NEWS! internal error in subroutine IBZKPT",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="KKSYM"; //usually goes with IBZKPT
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"Reciprocal lattice and k-lattice belong to different class of lattices",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="NKXYZ_IKPTD"; //usually goes with IBZKPT
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    bool found_nkxyz_ikptd=false; //break up for readability
    found_nkxyz_ikptd=(found_nkxyz_ikptd || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NKX>IKPTD",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_nkxyz_ikptd=(found_nkxyz_ikptd || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NKY>IKPTD",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_nkxyz_ikptd=(found_nkxyz_ikptd || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NKZ>IKPTD",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_warning=(found_warning && found_nkxyz_ikptd);
    xwarning.flag(scheme,found_warning);
    //
    scheme="INVGRP";  //usually goes with SYMPREC
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"VERY BAD NEWS! internal error in subroutine INVGRP",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="SYMPREC"; //usually goes with INVGR
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"inverse of rotation matrix was not found (increase SYMPREC)",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //VASP's internal symmetry routines END

    //VASP issues with RMM-DIIS START
    scheme="EDDRMM"; //CO20210315 - look here https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"WARNING in EDDRMM: call to ZHEGV failed",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="NUM_PROB";  //CO20210315 - look here https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS  //CO20210315 - num prob can be a big problem for the DOS: https://www.vasp.at/forum/viewtopic.php?f=3&t=18028
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"num prob",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="ZBRENT";  //CO20210315 - look here https://www.vasp.at/wiki/index.php/IALGO#RMM-DIIS
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"ZBRENT: can't locate minimum",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //VASP issues with RMM-DIIS END

    //CSLOSHING and NELM START
    //CSLOSHING and NELM warnings are similar, CSLOSHING will apply a fix before VASP finishes running, while NELM only cares about the LAST iteration
    //do not check turn on CSLOSHING with NELM, caused collisions with EDDDAV. it bounces back and forth between ALGO=VERFAST and ALGO=NORMAL
    //however, if CSLOSHING, turn on NELM, as the NELM patches will work for CSLOSHING too
    //[CO20210315 - OBSOLETE]//if there is a patch to be applied for the error (CSLOSHING does), then check both when vasp running and when it's not running
    //[CO20210315 - OBSOLETE]//check NELM first, and set CSLOSHING on if NELM, that way CSLOSHING patches are applied first (work for both errors)
    //the check for xwarning.flag("OUTCAR_INCOMPLETE")==false ensures we don't flag a run that was killed by --monitor_vasp, not xmessage.flag("REACHED_ACCURACY") (must be unconverged)
    scheme="NELM";  //CO20210315
    found_warning=(vasp_still_running==false && xwarning.flag("OUTCAR_INCOMPLETE")==false && vasp_oszicar_unconverged);  // check from OSZICAR
    xwarning.flag(scheme,found_warning);
    //
    scheme="CSLOSHING";
    found_warning=(KBIN::VASP_OSZICARUnconverging(xvasp.Directory)); // check from OSZICAR //xwarning.flag("NELM")
    xwarning.flag(scheme,found_warning);
    if(1){if(xwarning.flag("CSLOSHING")){xwarning.flag("NELM",true);}}
    //CSLOSHING and NELM END

    //ALL OTHERS (in alphabetic order) START
    scheme="BRMIX";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"BRMIX: very serious problems",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="CALC_FROZEN"; //CO20210315
    found_warning=(tmod_outcar>=SECONDS_STALE_OUTCAR);
    xwarning.flag(scheme,found_warning);
    //
    scheme="DAV";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"WARNING: Sub-Space-Matrix is not hermitian in DAV",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="DENTET";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"WARNING: DENTET: can't reach specified precision",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="EDDDAV";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"Error EDDDAV: Call to ZHEGV failed",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="EFIELD_PEAD";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"EFIELD_PEAD is too large",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"EFIELD_PEAD are too large for comfort",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) )); //20190704 - new VASP
    xwarning.flag(scheme,found_warning);
    //
    scheme="EXCCOR";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"ERROR FEXCF: supplied exchange-correlation table",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"ERROR FEXCP: supplied Exchange-correletion table",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) ));  //CO20210315 - some formatting changes between versions (also some bad spelling)
    xwarning.flag(scheme,found_warning);
    //
    scheme="GAMMA_SHIFT";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"shift your grid to Gamma",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));  //CO20190704 - captures both old/new versions of VASP
    xwarning.flag(scheme,found_warning);
    //
    scheme="LRF_COMMUTATOR";   // GET ALL TIMES
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"LRF_COMMUTATOR internal error: the vector",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="MEMORY";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,AFLOW_MEMORY_TAG,RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) || (XHOST.vflag_control.flag("KILL_VASP_OOM") && approaching_oom)));
    xwarning.flag(scheme,found_warning);
    //
    //on qrats, we see this
    //===================================================================================
    //=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
    //=   PID 27264 RUNNING AT prod-0004
    //=   EXIT CODE: 9
    //=   CLEANING UP REMAINING PROCESSES
    //=   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
    //===================================================================================
    //YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Killed (signal 9)
    //This typically refers to a problem with your application.
    //Please see the FAQ page for debugging suggestions
    //
    //on quser, we see nothing...
    //
    //on x, we see this
    //===================================================================================
    //=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
    //=   RANK 58 PID 63228 RUNNING AT node6
    //=   KILLED BY SIGNAL: 9 (Killed)
    //===================================================================================
    //
    bool found_bad_termination=aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition);
    bool found_exit_code=false; //specific exit code
    //
    scheme="MPICH0";  //0 just means that it is generic, maybe fix name later
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && found_bad_termination);
    xwarning.flag(scheme,found_warning);
    //
    scheme="MPICH11";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_exit_code=false;
    found_exit_code=(found_exit_code || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"EXIT CODE: 11",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_exit_code=(found_exit_code || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"KILLED BY SIGNAL: 11",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_warning=(found_warning && (found_bad_termination && found_exit_code));
    xwarning.flag(scheme,found_warning);
    //
    scheme="MPICH139";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_exit_code=false;
    found_exit_code=(found_exit_code || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"EXIT CODE: 139",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_exit_code=(found_exit_code || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"KILLED BY SIGNAL: 139",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_warning=(found_warning && (found_bad_termination && found_exit_code));
    xwarning.flag(scheme,found_warning);
    //
    scheme="MPICH174";  //CO20210315
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_exit_code=false;
    found_exit_code=(found_exit_code || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"EXIT CODE: 174",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_exit_code=(found_exit_code || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"KILLED BY SIGNAL: 174",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) );
    found_warning=(found_warning && (found_bad_termination && found_exit_code));
    xwarning.flag(scheme,found_warning);
    //
    scheme="NATOMS";  //look for problem for distance
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"distance between some ions is very small",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    //CO+ME20210315 - simplified
    //ME20190620 - Avoid changing NBANDS in the aflow.in file just because VASP throws the warning that NBANDS is changed because of NPAR. 
    //However, if you have that warning AND the error that the number of bands is not sufficient, aflow needs to act.
    scheme="NBANDS";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,scheme,RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"number of bands is not sufficient to hold all electrons",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"Number of bands NBANDS too small to hold electrons",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"Your highest band is occupied at some k-points",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition)) ));  // The NBANDS warning due to NPAR is not an error we want to fix, so set to false if found
    bool vasp_corrected=(aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"number of bands has been changed from the values supplied",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) || aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"now  NBANDS  =",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));  // Need explicit check or else the NPAR warning prevents this NBANDS error from being corrected
    found_warning=(found_warning && vasp_corrected==false);
    xwarning.flag(scheme,found_warning);
    //
    scheme="NPAR";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"please rerun with NPAR=",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition)); //not only npar==1
    xwarning.flag(scheme,found_warning);
    //
    scheme="NPARC";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NPAR = 4",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NPAR=number of cores",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition)));  // fix with NPAR=cores in MPI
    xwarning.flag(scheme,found_warning);
    //
    scheme="NPARN";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NPAR = 4",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"NPAR=number of nodes",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition)));  // fix with NPAR=nodes in MPI
    xwarning.flag(scheme,found_warning);
    //
    scheme="NPAR_REMOVE";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"Please remove the tag NPAR from the INCAR file",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));  //and restart the
    xwarning.flag(scheme,found_warning);
    //
    scheme="PSMAXN";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    //found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"WARNING: PSMAXN for non-local potential too small",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"REAL_OPT: internal ERROR",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="READ_KPOINTS_RD_SYM";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && (aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"ERROR in RE_READ_KPOINTS_RD",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"switch off symmetry",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition) ));
    xwarning.flag(scheme,found_warning);
    //
    scheme="REAL_OPT";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"REAL_OPT: internal ERROR",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="REAL_OPTLAY_1";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"REAL_OPTLAY: internal error (1)",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //
    scheme="ZPOTRF";
    found_warning=ReachedAccuracy2bool(scheme,xRequiresAccuracy,xmessage,vasp_still_running);
    found_warning=(found_warning && aurostd::substring_present_file_FAST(xvasp.Directory+"/"+DEFAULT_VASP_OUT,"LAPACK: Routine ZPOTRF failed",RemoveWS,case_insensitive,expect_near_end,grep_stop_condition));
    xwarning.flag(scheme,found_warning);
    //ALL OTHERS (in alphabetic order) END

    //check at the end (only if other warnings are found and can be corrected)
    //this is actually a critical check, what appeared to be vasp stuck was actually aflow stuck reading (cat/sed/grep) a VERY long vasp.out full of warnings
    //a single check for errors could take hours
    //files this large will put a lot of stress on cat/grep for warnings
    //vasp.out should NEVER get this large, means its filled with warnings/errors
    //do not look for other warnings, some might depend on xmessage.flag(\"REACHED_ACCURACY\")
    //make sure BYTES_MAX_VASP_OUT is as LARGE as possible, we want vasp to exit gracefully whenever possible. real errors might appear only after long list of warnings.
    //should only be an error when vasp is running, there is no treatment for this warning except to kill vasp
    scheme="OUTPUT_LARGE";  //CO20210315
    found_warning=(vasp_still_running==true && fsize_vaspout>=BYTES_MAX_VASP_OUT);  //100MB, make aflowrc parameter, 1GB is too large for NFS mounted nodes (NFS has to push the entire file over cable), when set to 100MB it still took 3 hours //CO20210315 - these stats were pre-renice, the processing time is much better now
    xwarning.flag(scheme,found_warning);

    if(LDEBUG){aus << __AFLOW_FUNC__ << " [2]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    if(renice){aurostd::ProcessRenice(vasp_bin,0);} //renice vasp back to normal after grep

    //CO20210315 - this bit must be done before wdebug
    //CO20210315 - only ignore KKSYM if OUTCAR is not registered as incomplete
    if( ( xwarning.flag("KKSYM") ) && ( xwarning.flag("OUTCAR_INCOMPLETE") ) &&
        ((ispin_current==2 && isym_current==-1) || (ispin_current==1 && isym_current==0))){  //CO20200624 - needs to change if we do magnetic systems //ispind_current==2 &&
      xmonitor.flag("IGNORING_WARNINGS:KKSYM",FALSE);
    }

    bool wdebug=(FALSE && LDEBUG);
    if(1) {
      if(LDEBUG){aus << __AFLOW_FUNC__ << " printing warnings" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      if(!xmonitor.flag("IGNORING_MESSAGES:REACHED_ACCURACY") && (wdebug || xmessage.flag("REACHED_ACCURACY"))) aus << "MMMMM  MESSAGE xmessage.flag(\"REACHED_ACCURACY\")=" << xmessage.flag("REACHED_ACCURACY") << endl;
      if(wdebug) aus << "MMMMM  MESSAGE VASP_release=" << xmessage.getattachedscheme("SVERSION") << endl;
      if(wdebug) aus << "MMMMM  MESSAGE VASP_version=" << xmessage.getattachedscheme("DVERSION") << endl;
      if(wdebug) aus << "MMMMM  MESSAGE AFLOW_version=" << AFLOW_VERSION << endl;
      if(wdebug || xwarning.flag("OUTCAR_INCOMPLETE")) aus << "WWWWW  WARNING xwarning.flag(\"OUTCAR_INCOMPLETE\")=" << xwarning.flag("OUTCAR_INCOMPLETE") << endl; //CO20201111
      //
      if(wdebug || xwarning.flag("BRMIX")) aus << "WWWWW  WARNING xwarning.flag(\"BRMIX\")=" << xwarning.flag("BRMIX") << endl;
      if(wdebug || xwarning.flag("CSLOSHING")) aus << "WWWWW  WARNING xwarning.flag(\"CSLOSHING\")=" << xwarning.flag("CSLOSHING") << endl;
      if(wdebug || xwarning.flag("CALC_FROZEN")) aus << "WWWWW  WARNING xwarning.flag(\"CALC_FROZEN\")=" << xwarning.flag("CALC_FROZEN") << endl; //CO20210135
      if(wdebug || xwarning.flag("DAV")) aus << "WWWWW  WARNING xwarning.flag(\"DAV\")=" << xwarning.flag("DAV") << endl;
      if(wdebug || xwarning.flag("DENTET")) aus << "WWWWW  WARNING xwarning.flag(\"DENTET\")=" << xwarning.flag("DENTET") << endl;
      if(wdebug || xwarning.flag("EDDDAV")) aus << "WWWWW  WARNING xwarning.flag(\"EDDDAV\")=" << xwarning.flag("EDDDAV") << endl;
      if(wdebug || xwarning.flag("EDDRMM")) aus << "WWWWW  WARNING xwarning.flag(\"EDDRMM\")=" << xwarning.flag("EDDRMM") << endl;
      if(wdebug || xwarning.flag("EFIELD_PEAD")) aus << "WWWWW  WARNING xwarning.flag(\"EFIELD_PEAD\")=" << xwarning.flag("EFIELD_PEAD") << endl;
      if(wdebug || xwarning.flag("EXCCOR")) aus << "WWWWW  WARNING xwarning.flag(\"EXCCOR\")=" << xwarning.flag("EXCCOR") << endl;
      if(wdebug || xwarning.flag("GAMMA_SHIFT")) aus << "WWWWW  WARNING xwarning.flag(\"GAMMA_SHIFT\")=" << xwarning.flag("GAMMA_SHIFT") << endl;
      if(!xmonitor.flag("IGNORING_WARNINGS:KKSYM") && !xmonitor.flag("IGNORING_WARNINGS:IBZKPT") && (wdebug || xwarning.flag("IBZKPT"))) aus << "WWWWW  WARNING xwarning.flag(\"IBZKPT\")=" << xwarning.flag("IBZKPT") << endl;
      if(wdebug || xwarning.flag("INVGRP")) aus << "WWWWW  WARNING xwarning.flag(\"INVGRP\")=" << xwarning.flag("INVGRP") << endl;
      if(!xmonitor.flag("IGNORING_WARNINGS:KKSYM") && (wdebug || xwarning.flag("KKSYM"))) aus << "WWWWW  WARNING xwarning.flag(\"KKSYM\")=" << xwarning.flag("KKSYM") << endl;
      if(wdebug || xwarning.flag("LRF_COMMUTATOR")) aus << "WWWWW  WARNING xwarning.flag(\"LRF_COMMUTATOR\")=" << xwarning.flag("LRF_COMMUTATOR") << endl;
      if(wdebug || xwarning.flag("MEMORY")) aus << "WWWWW  WARNING xwarning.flag(\"MEMORY\")=" << xwarning.flag("MEMORY") << endl;
      if(wdebug || xwarning.flag("MPICH0")) aus << "WWWWW  WARNING xwarning.flag(\"MPICH0\")=" << xwarning.flag("MPICH0") << endl;
      if(wdebug || xwarning.flag("MPICH11")) aus << "WWWWW  WARNING xwarning.flag(\"MPICH11\")=" << xwarning.flag("MPICH11") << endl;
      if(wdebug || xwarning.flag("MPICH139")) aus << "WWWWW  WARNING xwarning.flag(\"MPICH139\")=" << xwarning.flag("MPICH139") << endl;
      if(wdebug || xwarning.flag("MPICH174")) aus << "WWWWW  WARNING xwarning.flag(\"MPICH174\")=" << xwarning.flag("MPICH174") << endl;
      if(wdebug || xwarning.flag("NATOMS")) aus << "WWWWW  WARNING xwarning.flag(\"NATOMS\")=" << xwarning.flag("NATOMS") << endl;
      if(wdebug || xwarning.flag("NBANDS")) aus << "WWWWW  WARNING xwarning.flag(\"NBANDS\")=" << xwarning.flag("NBANDS") << endl;
      if(wdebug || xwarning.flag("NELM")) aus << "WWWWW  WARNING xwarning.flag(\"NELM\")=" << xwarning.flag("NELM") << endl;
      if(wdebug || xwarning.flag("NIRMAT")) aus << "WWWWW  WARNING xwarning.flag(\"NIRMAT\")=" << xwarning.flag("NIRMAT") << endl;
      if(wdebug || xwarning.flag("NKXYZ_IKPTD")) aus << "WWWWW  WARNING xwarning.flag(\"NKXYZ_IKPTD\")=" << xwarning.flag("NKXYZ_IKPTD") << endl;
      if(wdebug || xwarning.flag("NPAR")) aus << "WWWWW  WARNING xwarning.flag(\"NPAR\")=" << xwarning.flag("NPAR") << endl;
      if(!xmonitor.flag("IGNORING_WARNINGS:NPARC") && (wdebug || xwarning.flag("NPARC"))) aus << "WWWWW  WARNING xwarning.flag(\"NPARC\")=" << xwarning.flag("NPARC") << endl;
      if(!xmonitor.flag("IGNORING_WARNINGS:NPARN") && (wdebug || xwarning.flag("NPARN"))) aus << "WWWWW  WARNING xwarning.flag(\"NPARN\")=" << xwarning.flag("NPARN") << endl;
      if(wdebug || xwarning.flag("NPAR_REMOVE")) aus << "WWWWW  WARNING xwarning.flag(\"NPAR_REMOVE\")=" << xwarning.flag("NPAR_REMOVE") << endl;
      if(wdebug || xwarning.flag("NUM_PROB")) aus << "WWWWW  WARNING xwarning.flag(\"NUM_PROB\")=" << xwarning.flag("NUM_PROB") << endl;  //CO20210315
      if(wdebug || xwarning.flag("OUTPUT_LARGE")) aus << "WWWWW  WARNING xwarning.flag(\"OUTPUT_LARGE\")=" << xwarning.flag("OUTPUT_LARGE") << endl;
      if(wdebug || xwarning.flag("PSMAXN")) aus << "WWWWW  WARNING xwarning.flag(\"PSMAXN\")=" << xwarning.flag("PSMAXN") << endl;
      if(wdebug || xwarning.flag("READ_KPOINTS_RD_SYM")) aus << "WWWWW  WARNING xwarning.flag(\"READ_KPOINTS_RD_SYM\")=" << xwarning.flag("READ_KPOINTS_RD_SYM") << endl;
      if(wdebug || xwarning.flag("REAL_OPT")) aus << "WWWWW  WARNING xwarning.flag(\"REAL_OPT\")=" << xwarning.flag("REAL_OPT") << endl;
      if(wdebug || xwarning.flag("REAL_OPTLAY_1")) aus << "WWWWW  WARNING xwarning.flag(\"REAL_OPTLAY_1\")=" << xwarning.flag("REAL_OPTLAY_1") << endl;
      if(wdebug || xwarning.flag("SGRCON")) aus << "WWWWW  WARNING xwarning.flag(\"SGRCON\")=" << xwarning.flag("SGRCON") << endl;
      if(wdebug || xwarning.flag("SYMPREC")) aus << "WWWWW  WARNING xwarning.flag(\"SYMPREC\")=" << xwarning.flag("SYMPREC") << endl;
      if(wdebug || xwarning.flag("ZBRENT")) aus << "WWWWW  WARNING xwarning.flag(\"ZBRENT\")=" << xwarning.flag("ZBRENT") << endl;  //CO20210315
      if(wdebug || xwarning.flag("ZPOTRF")) aus << "WWWWW  WARNING xwarning.flag(\"ZPOTRF\")=" << xwarning.flag("ZPOTRF") << endl;
      if(LDEBUG){cerr << aus.str();}
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    //CO20210601 - only print reached accuracy message once, if calc is frozen this can be printed out many times
    if(xmessage.flag("REACHED_ACCURACY")){xmonitor.flag("IGNORING_MESSAGES:REACHED_ACCURACY",TRUE);} 

    //[CO20210315 - doesn't work]//CO20210315 - this appears to fix ICSD/LIB/CUB/Ag1Sm1_ICSD_604546
    //[CO20210315 - doesn't work]if(xwarning.flag("EDDRMM") && xwarning.flag("ZPOTRF")){ // fix EDDRMM first
    //[CO20210315 - doesn't work]  aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"ZPOTRF\"): prioritizing xwarning.flag(\"EDDRMM\")" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    //[CO20210315 - doesn't work]  xwarning.flag("ZPOTRF",FALSE);
    //[CO20210315 - doesn't work]  //we don't need an xmonitor here, this is only for prioritizing errors
    //[CO20210315 - doesn't work]}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //must do before RMM_DIIS and ROTMAT
    if(1||vasp_monitor_running){ //CO20210315 - might consider running this always, come back and test
      //vasp_monitor will kill vasp prematurely, triggering warnings that require xmessage.flag("REACHED_ACCURACY") (false positives)
      //prioritize the other warnings
      //might be a false positive without vasp_monitor_running, vasp can also trigger its own premature exiting
      //in the future, if these must be turned off, be careful of switching between DAV (EDDDAV) and RMM-DIIS problems
      //one requires ALGO=NORMAL, the other =VERYFAST and they are mutually exclusive
      //therefore, create an xwarnings_fixed which stores which warnings have been fixed previously
      //if warnings like RMM-DIIS have been fixed before, are not problems now, and we encounter CSLOSHING, we should NOT try =NORMAL
      uint n_require_accuracy=0;
      for(i=0;i<xwarning.vxscheme.size();i++){
        if(xRequiresAccuracy.flag(xwarning.vxscheme[i])){n_require_accuracy++;}
      }
      uint n_derivative=0;
      vector<string> vwarnings_derivative; //these warnings are derivative: e.g., an incomplete outcar could be the result of many errors, including those requiring xmessage.flag("REACHED_ACCURACY")
      aurostd::string2tokens("OUTCAR_INCOMPLETE,CALC_FROZEN,OUTPUT_LARGE",vwarnings_derivative,",");
      for(i=0;i<vwarnings_derivative.size();i++){
        if(xwarning.flag(vwarnings_derivative[i])){n_derivative++;}
      }
      if(LDEBUG){
        aus << __AFLOW_FUNC__ << " xwarning.vxscheme=" << aurostd::joinWDelimiter(xwarning.vxscheme,",") << endl;
        aus << __AFLOW_FUNC__ << " xwarning.vxscheme.size()=" << xwarning.vxscheme.size() << endl;
        aus << __AFLOW_FUNC__ << " n_require_accuracy=" << n_require_accuracy << endl;
        aus << __AFLOW_FUNC__ << " n_derivative=" << n_derivative << endl;
        cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(xwarning.vxscheme.size()>(n_require_accuracy+n_derivative)){ //this means we have some real errors inside
        vector<string> xwarning_vxscheme=xwarning.vxscheme; //make a copy since we're deleting entries of the vector
        for(i=xwarning_vxscheme.size()-1;i<xwarning_vxscheme.size();i--){  //go backwards since we're removing entries
          if(xRequiresAccuracy.flag(xwarning_vxscheme[i])){
            aus << "MMMMM  MESSAGE ignoring xwarning.flag(\""+xwarning_vxscheme[i]+"\"): prioritizing other warnings first (requires xmessage.flag(\"REACHED_ACCURACY\"))" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //; possible false positive
            xwarning.flag(xwarning_vxscheme[i],FALSE);
            //we don't need an xmonitor here, this is only for prioritizing errors
          }
        }
      }
    }
    //NUM_PROB and ZBRENT are soft errors, only fix if the OUTCAR is incomplete (and not if accuracy not reached)
    //the ZBRENT patch (changing IBRION) may be beneficial
    //conjugate gradient might fail too close to the minimum, so this patch might be good for relax2
    //however, it has shown to steer other calculations in bad directions, leading to other warnings that cannot be fixed
    //better not to over-correct
    //better to check for convergence of the calculation later with --xplug
    if(xwarning.flag("OUTCAR_INCOMPLETE")==false){
      if(xwarning.flag("NUM_PROB")){
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NUM_PROB\"): OUTCAR is complete (soft warning)" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xwarning.flag("NUM_PROB",FALSE);
      }
      if(xwarning.flag("ZBRENT")){
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"ZBRENT\"): OUTCAR is complete (soft warning)" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xwarning.flag("ZBRENT",FALSE);
      }
    }
    //[CO20210315 - OBSOLETE]//ZBRENT is a soft error, fix anything else first
    //[CO20210315 - OBSOLETE]if(xwarning.flag("ZBRENT") && 
    //[CO20210315 - OBSOLETE]    ((xwarning.flag("OUTCAR_INCOMPLETE")==false && xwarning.vxscheme.size()>1) || (xwarning.flag("OUTCAR_INCOMPLETE") && xwarning.vxscheme.size()>2))){
    //[CO20210315 - OBSOLETE]  aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"ZBRENT\"): prioritizing other warnings first" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    //[CO20210315 - OBSOLETE]  xwarning.flag("ZBRENT",FALSE);
    //[CO20210315 - OBSOLETE]  //we don't need an xmonitor here, this is only for prioritizing errors
    //[CO20210315 - OBSOLETE]}

    //CO20210315 - only ignore KKSYM if OUTCAR is not registered as incomplete
    if( ( xwarning.flag("KKSYM") ) && 
        ((ispin_current==2 && isym_current==-1) || (ispin_current==1 && isym_current==0))){  //CO20200624 - needs to change if we do magnetic systems //ispind_current==2 &&
      if(!xmonitor.flag("IGNORING_WARNINGS:KKSYM")){
        aus << "MMMMM  MESSAGE ignoring KKSYM warning: ISYM=" << isym_current << " ISPIN=" << ispin_current << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xmonitor.flag("IGNORING_WARNINGS:KKSYM",TRUE);  //so we don't clog output files
      }
      xwarning.flag("KKSYM",FALSE);
      if(xwarning.flag("IBZKPT") && xwarning.flag("OUTCAR_INCOMPLETE")==false){  //goes with KKSYM
        aus << "MMMMM  MESSAGE ignoring IBZKPT warning associated with KKSYM: ISYM=" << isym_current << " ISPIN=" << ispin_current << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xwarning.flag("IBZKPT",FALSE);
      }
    }

    //do just before RMM_DIIS and ROTMAT
    //we generally treat CALC_FROZEN as a out-of-memory error UNLESS the calculation has reached its accuracy (finished) and the OUTCAR is incomplete (odd occurrence, noticed on qrats so far)
    //only convert CALC_FROZEN to MEMORY if it's the only warning (OUTCAR_INCOMPLETE is derivative)
    if(xwarning.flag("CALC_FROZEN") && xmessage.flag("REACHED_ACCURACY")==false){
      if((xwarning.vxscheme.size()==1)||(xwarning.flag("OUTCAR_INCOMPLETE") && xwarning.vxscheme.size()==2)){
        xwarning.flag("CALC_FROZEN",false);xwarning.flag("MEMORY",true);
        aus << "MMMMM  MESSAGE treating xwarning.flag(\"CALC_FROZEN\") as xwarning.flag(\"MEMORY\")" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //; possible false positive
      }
    }

    //the memory might ramp up too quickly, in which the OOM killer will kill vasp
    //if we are lucky, we'll get the "BAD TERMINATION..." warning in the vasp.out (works on qrats)
    if(xwarning.flag("MPICH0") && xmessage.flag("REACHED_ACCURACY")==false){
      if((xwarning.vxscheme.size()==1)||(xwarning.flag("OUTCAR_INCOMPLETE") && xwarning.vxscheme.size()==2)){
        xwarning.flag("MPICH0",false);xwarning.flag("MEMORY",true);
        aus << "MMMMM  MESSAGE treating xwarning.flag(\"MPICH0\") as xwarning.flag(\"MEMORY\")" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //; possible false positive
      }
    }

    //do last
    bool rmm_diis_warning=false;
    rmm_diis_warning=(rmm_diis_warning || xwarning.flag("EDDRMM"));
    rmm_diis_warning=(rmm_diis_warning || xwarning.flag("NUM_PROB"));
    rmm_diis_warning=(rmm_diis_warning || xwarning.flag("ZBRENT"));
    //CO20210315 - can probably add others to this list as well
    xwarning.flag("RMM_DIIS",rmm_diis_warning);

    bool rotmat_warning=false;
    rotmat_warning=(rotmat_warning || xwarning.flag("SGRCON"));
    rotmat_warning=(rotmat_warning || xwarning.flag("NIRMAT"));
    rotmat_warning=(rotmat_warning || xwarning.flag("IBZKPT"));
    rotmat_warning=(rotmat_warning || xwarning.flag("KKSYM"));
    rotmat_warning=(rotmat_warning || xwarning.flag("INVGRP"));
    rotmat_warning=(rotmat_warning || xwarning.flag("SYMPREC"));
    //CO20210315 - can probably add others to this list as well
    xwarning.flag("ROTMAT",rotmat_warning);
    ///////////////////////////////////////////////////////////////////////////////////////////////

    if(1){
      if(wdebug || xwarning.flag("RMM_DIIS")) aus << "WWWWW  WARNING xwarning.flag(\"RMM_DIIS\")=" << xwarning.flag("RMM_DIIS") << endl;  //CO20210315
      if(wdebug || xwarning.flag("ROTMAT")) aus << "WWWWW  WARNING xwarning.flag(\"ROTMAT\")=" << xwarning.flag("ROTMAT") << endl;
      if(LDEBUG){cerr << aus.str();}
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    //[CO20210315 - CSLOSHING already picks this up]xOUTCAR xout(xvasp.Directory+"/OUTCAR",true);  //quiet, there might be issues with halfway-written OUTCARs
    //[CO20210315 - CSLOSHING already picks this up]int NBANDS=xout.NBANDS;
    //[CO20210315 - CSLOSHING already picks this up]int NELM=xout.NELM;
    //[CO20210315 - CSLOSHING already picks this up]int NSTEPS=KBIN::VASP_getNSTEPS(xvasp.Directory+"/OSZICAR");
    //[CO20210315 - CSLOSHING already picks this up]aus << "MMMMM  MESSAGE NBANDS=" << NBANDS << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);aus.str(std::string());aus.clear();
    //[CO20210315 - CSLOSHING already picks this up]aus << "MMMMM  MESSAGE NELM=" << NELM << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);aus.str(std::string());aus.clear();
    //[CO20210315 - CSLOSHING already picks this up]aus << "MMMMM  MESSAGE NSTEPS=" << NSTEPS << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);aus.str(std::string());aus.clear();
    //[CO20210315 - APL can go way above]if(xwarning.flag("NBANDS") && NBANDS>1000) xwarning.flag("NBANDS",FALSE); // for safety
    //[CO20210315 - these flags are NOT found INCAR, look instead in vflags.KBIN_VASP_RUN.flag(), someone needs to recheck these operations]if(xwarning.flag("NBANDS") && aurostd::substring2bool(xvasp.INCAR,"DIELECTRIC_STATIC") && NBANDS>1000) xwarning.flag("NBANDS",FALSE); // for safety
    //[CO20210315 - these flags are NOT found INCAR, look instead in vflags.KBIN_VASP_RUN.flag(), someone needs to recheck these operations]if(xwarning.flag("NBANDS") && aurostd::substring2bool(xvasp.INCAR,"DIELECTRIC_DYNAMIC") && NBANDS>1000) xwarning.flag("NBANDS",FALSE); // for safety
    //[CO20210315 - these flags are NOT found INCAR, look instead in vflags.KBIN_VASP_RUN.flag(), someone needs to recheck these operations]if(xwarning.flag("NBANDS") && aurostd::substring2bool(xvasp.INCAR,"DSCF") && NBANDS>1000) xwarning.flag("NBANDS",FALSE); // for safety
    //[CO20210315 - CSLOSHING already picks this up]if(NELM!=0 && NSTEPS!=0 && NSTEPS>=NELM) { xwarning.flag("NELM",TRUE); } else { xwarning.flag("NELM",FALSE); }

    //WARNINGS STOP

    if(LDEBUG){aus << __AFLOW_FUNC__ << " [3]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    //decide which warnings to ignore

    if(xwarning.flag("MPICH11") && xwarning.flag("NBANDS")){ // fix MPICH11 first
      aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NBANDS\"): prioritizing xwarning.flag(\"MPICH11\")" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xwarning.flag("NBANDS",FALSE);
      //we don't need an xmonitor here, this is only for prioritizing errors
    }
    if(xwarning.flag("NPARC") && (aurostd::kvpair2utype<int>(xvasp.INCAR,"NPAR","=")==2 || // dont bother for small NPAR
          aurostd::kvpair2string(xvasp.INCAR,"LRPA","=")==".TRUE." ||  //CO20210315 - would be better to check if .TRUE. or .FALSE.
          aurostd::kvpair2string(xvasp.INCAR,"LEPSILON","=")==".TRUE." ||  //CO20210315 - would be better to check if .TRUE. or .FALSE.
          aurostd::kvpair2string(xvasp.INCAR,"LOPTICS","=")==".TRUE.")){  // dont touch NPARC if LRPA or LEPSILON or LOPTICS necessary
      if(!xmonitor.flag("IGNORING_WARNINGS:NPARC")){
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NPARC\"): found either LRPA, LEPSILON, or LOPTICS" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xmonitor.flag("IGNORING_WARNINGS:NPARC",TRUE);  //so we don't clog output files
      }
      xwarning.flag("NPARC",FALSE);
    }
    if(xwarning.flag("NPARN") && (aurostd::kvpair2string(xvasp.INCAR,"LRPA","=")==".TRUE." ||  //CO20210315 - would be better to check if .TRUE. or .FALSE.
          aurostd::kvpair2string(xvasp.INCAR,"LEPSILON","=")==".TRUE." ||  //CO20210315 - would be better to check if .TRUE. or .FALSE.
          aurostd::kvpair2string(xvasp.INCAR,"LOPTICS","=")==".TRUE.")){  // dont touch NPARN if LRPA or LEPSILON or LOPTICS necessary
      if(!xmonitor.flag("IGNORING_WARNINGS:NPARN")){
        aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"NPARN\"): found either LRPA, LEPSILON, or LOPTICS" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xmonitor.flag("IGNORING_WARNINGS:NPARN",TRUE);  //so we don't clog output files
      }
      xwarning.flag("NPARN",FALSE);
    }
    //[CO20210315 - flag would not be turned on above]if(xmessage.flag("REACHED_ACCURACY") && xwarning.flag("IBZKPT")){  //CO20210315 - I guess it's not an issue then?
    //[CO20210315 - flag would not be turned on above]  if(!xmonitor.flag("IGNORING_WARNINGS:IBZKPT")){
    //[CO20210315 - flag would not be turned on above]    aus << "MMMMM  MESSAGE ignoring xwarning.flag(\"IBZKPT\"): VASP calculation achieved required accuracy anyway" << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    //[CO20210315 - flag would not be turned on above]    xmonitor.flag("IGNORING_WARNINGS:IBZKPT",TRUE); //so we don't clog output files
    //[CO20210315 - flag would not be turned on above]  }
    //[CO20210315 - flag would not be turned on above]  xwarning.flag("IBZKPT",FALSE);
    //[CO20210315 - flag would not be turned on above]}

    //[CO20210315 - OBSOLETE, we prioritize the order below]if(xwarning.flag("NKXYZ_IKPTD")){xwarning.flag("IBZKPT",FALSE);} // priority
    //[try NIRMAT first]if(xwarning.flag("NIRMAT") && xwarning.flag("SGRCON")) xwarning.flag("SGRCON",FALSE);
    //[no must fix the LATTICE]if(xwarning.flag("EDDRMM")) xwarning.flag("ZPOTRF",FALSE);

    if(LDEBUG){aus << __AFLOW_FUNC__ << " [4]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Error2Fix(const string& error,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) { //CO20210315
    int submode=0; //default
    bool try_last_ditch_efforts=true; //default
    return VASP_Error2Fix(error,error,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,const string& mode,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) {  //CO20210315
    int submode=0; //default
    bool try_last_ditch_efforts=true; //default
    return VASP_Error2Fix(error,mode,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,int& submode,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) { //CO20210315
    bool try_last_ditch_efforts=true; //default
    return VASP_Error2Fix(error,error,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,const string& mode,int& submode,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) {  //CO20210315
    bool try_last_ditch_efforts=true; //default
    return VASP_Error2Fix(error,mode,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,bool try_last_ditch_efforts,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) { //CO20210315
    int submode=0; //default
    return VASP_Error2Fix(error,error,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,const string& mode,bool try_last_ditch_efforts,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) {  //CO20210315
    int submode=0; //default
    return VASP_Error2Fix(error,mode,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,int& submode,bool try_last_ditch_efforts,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) { //CO20210315
    return VASP_Error2Fix(error,error,submode,try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);
  }
  bool VASP_Error2Fix(const string& error,const string& mode,int& submode,bool try_last_ditch_efforts,_xvasp &xvasp,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) {  //CO20210315
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    stringstream aus;

    if(LDEBUG){aus << __AFLOW_FUNC__ << " [CHECK " << error << " PROBLEMS]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    bool apply_fix=xwarning.flag(error);
    bool VERBOSE=(FALSE || XHOST.vflag_control.flag("MONITOR_VASP")==false || LDEBUG);
    if(apply_fix && xfixed.flag("ALL")){
      if(LDEBUG){aus << __AFLOW_FUNC__ << " xfixed.flag(\"ALL\")==TRUE: skipping " << error << " fix" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      apply_fix=false;
    }
    if(apply_fix && vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:ALL")){
      if(LDEBUG){aus << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:ALL\")==TRUE: skipping " << error << " fix" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      apply_fix=false;
    }
    if(apply_fix && vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag("ERROR:"+error)){
      if(LDEBUG){aus << __AFLOW_FUNC__ << " vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"ERROR:" << error << "\")==TRUE: skipping " << error << " fix" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
      apply_fix=false;
    }
    if(apply_fix && VERBOSE){
      aus << "MMMMM  MESSAGE attempting to fix ERROR=\"" << error << "\" (mode=\"" << mode << "\")" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  //CO20210315 - do not reference submode after KBIN::XVASP_Afix(), (submode>=0?" (SUBMODE="+aurostd::utype2string(submode)+")":"")
    }
    //do not reference submode below KBIN::XVASP_Afix(), as it will have been incremented (perhaps by 2)
    //if KBIN::XVASP_Afix() fails, submode will be restored to original, so it is okay to reference for the LDEBUG statement
    if(apply_fix && !KBIN::XVASP_Afix(mode,submode,try_last_ditch_efforts,xfixed,xvasp,kflags,vflags,aflags,FileMESSAGE)){   //CO20210315
      if(LDEBUG){aus << __AFLOW_FUNC__ << " KBIN::XVASP_Afix(mode=\"" << mode << "\"" << (submode>=0?",submode="+aurostd::utype2string(submode):"") << ") failed" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}  //if KBIN::XVASP_Afix() fails, submode will be restored to original, so it is okay to reference for the LDEBUG statement
      apply_fix=false;
    }
    if(apply_fix && VERBOSE){
      aus << "MMMMM  MESSAGE fix applied for ERROR=\"" << error << "\" (mode=\"" << mode << "\")" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  //CO20210315 - do not reference submode after KBIN::XVASP_Afix(), (submode>=0?" (SUBMODE="+aurostd::utype2string(submode)+")":"")
    }
    if(apply_fix){
      if(xvasp.aopts.flag("FLAG::AFIX_DRYRUN")==false){
        KBIN::VASP_Error(xvasp,"WWWWW  ERROR "+__AFLOW_FUNC__+" \""+error+"\" problems"+Message(__AFLOW_FILE__,aflags));
        //[CO20210315 - old style]xfixed.flag(error,TRUE);
      }
      xfixed.flag("ALL",TRUE);
    }
    return apply_fix;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_FixErrors(_xvasp &xvasp,aurostd::xoption& xmessage,aurostd::xoption& xwarning,aurostd::xoption& xfixed,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE){
    //a note about the fixes below
    //they generally compound, which I believe is the right approach
    //however, there might be some options which conflict
    //add KBIN::XVASP_INCAR_REMOVE_ENTRY() as necessary
    //check also submode fixes

    bool fixed_applied=false;
    bool try_last_ditch_efforts=true;
    uint i=0;

    for(i=0;i<2&&fixed_applied==false;i++){ //for loop goes twice, once with try_last_ditch_efforts==false, then again with ==true
      try_last_ditch_efforts=(i==1);

      //check NBANDS/LRF_COMMUTATOR problems immediately
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NBANDS",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      //[CO20210315 - fix previously removed]KBIN::VASP_Error2Fix("LRF_COMMUTATOR",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);

      //fix symmetry issues next
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("GAMMA_SHIFT",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      // ********* APPLY PREFERRED SYMMETRY FIXES ******************  //all of these must come before ROTMAT
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NKXYZ_IKPTD",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));  //must come before IBZKPT
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("KKSYM",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));  //must come before IBZKPT
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("IBZKPT",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("SYMPREC",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));  //must come before INVGRP
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("INVGRP","SYMPREC",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("SGRCON","SYMPREC",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      // ********* APPLY GENERIC SYMMETRY FIXES ******************
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("ROTMAT",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));

      //fix MPI/NPAR problems next
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("MPICH11",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("MPICH139",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("MPICH174",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE)); //CO20210315 - testing, exit code 174 looks like an error on the node, basically try rerunning with more memory
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NPAR",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NPARC",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NPARN",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NPAR_REMOVE",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));

      //all other fixes, no priority here (alphabetic order)
      // ********* APPLY OTHER FIXES ******************
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("BRMIX",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CSLOSHING",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE)); //must come before NELM
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("DAV",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("DENTET",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("EDDDAV",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("EDDRMM","RMM_DIIS",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("EFIELD_PEAD",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("EXCCOR",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("MEMORY",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NATOMS",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("PSMAXN",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("REAL_OPTLAY_1","LREAL",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("REAL_OPT","LREAL",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("ZPOTRF",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));

      //patch only if above warnings are not patched first
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("NELM",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));

      //apply these last if no other fixes worked
      // ********* APPLY PREFERRED RMM-DIIS FIXES ******************  //all of these must come before RMM-DIIS
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("ZBRENT",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      // ********* APPLY GENERIC RMM-DIIS FIXES ******************
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("RMM_DIIS",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));

      //[CO20210315 - do not apply patches for frozen calc]//CO20210315 - do last, fixes assume out-of-memory error
      //[CO20210315 - do not apply patches for frozen calc]fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","MEMORY",try_last_ditch_efforts,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
    }

    //sometimes VASP will die after it prints the REACHED_ACCURACY state, but before the OUTCAR is finished (not sure why)
    //this is rare...
    //might be a threading/NFS issue
    //this means any errors inside that require REACHED_ACCURACY will not be triggered
    //in this case, try restarting the calculation from CONTCAR
    //lowering NCPUS has been shown to work, indicating that this is indeed a threading/mpi issue
    //try from the most relaxed CONTCAR to save time
    //this is NOT a magic bullet, it looks like the threading solution works for some structures and not others
    //I am leaving "THREADS" vs. going to "MEMORY" solutions which will change NBANDS, KPOINTS, etc.
    //better to run on another machine/different binary
    if(fixed_applied==false && xwarning.flag("CALC_FROZEN") && xmessage.flag("REACHED_ACCURACY") && xwarning.flag("OUTCAR_INCOMPLETE")){
      //[CO20210315 - not shown to work]fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","RESTART_CALC",false,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      //[CO20210621 - not shown to work (alone)]fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","RECYCLE_CONTCAR",false,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
      fixed_applied=(fixed_applied || KBIN::VASP_Error2Fix("CALC_FROZEN","THREADS",false,xvasp,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE));
    }

    //print out all xfixed BEFORE adding "ALL"
    std::sort(xfixed.vxscheme.begin(),xfixed.vxscheme.end()); //sort for printing
    //print ALL first
    stringstream aus;
    if(xfixed.flag("ALL")){aus << "MMMMM  MESSAGE xfixed.flag(\"ALL\")=" << xfixed.flag("ALL") << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
    if(xfixed.flag("ALL") || XHOST.vflag_control.flag("MONITOR_VASP")==false){  //very important that we print all xfixed even if not "ALL" for LOCK file, --monitor_vasp reads the LOCK looking here. otherwise, only print if "ALL"
      for(uint i=0;i<xfixed.vxscheme.size();i++){
        if(xfixed.vxscheme[i]=="ALL"){continue;}  //already done above
        aus << "MMMMM  MESSAGE xfixed.flag(\""+xfixed.vxscheme[i]+"\")=" << xfixed.flag(xfixed.vxscheme[i]) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    return fixed_applied;
  }
}

namespace KBIN {
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    string function="KBIN::VASP_Run";
    ostringstream aus_exec,aus;

    if(LDEBUG){aus << __AFLOW_FUNC__ << " BEGIN" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    if(XHOST.AVOID_RUNNING_VASP){  //CO20200624
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"VASP should NOT be running",_INPUT_ILLEGAL_);  //better to throw to avoid VASP_Backup(), etc.
      //return false;
    }

    xoption xwarning,xfixed,xmessage;
    bool vasp_start=TRUE;
    aurostd::StringstreamClean(aus_exec);
    aurostd::StringstreamClean(aus);
    int nrun=0,maxrun=100; //CO20210315 - increase from 15 to 100 //NBANDS can take a lot of iterations to reach goal

    // get CPUS from PBS/SLURM
    // string ausenv;
    aus << "DDDDD  PBS_NUM_PPN=" << XHOST.PBS_NUM_PPN << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  PBS_NNODES=" << XHOST.PBS_NNODES << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  SLURM_CPUS_ON_NODE=" << XHOST.SLURM_CPUS_ON_NODE << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  SLURM_NNODES=" << XHOST.SLURM_NNODES << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  SLURM_NTASKS=" << XHOST.SLURM_NTASKS << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  kflags.KBIN_MPI_NCPUS=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  XHOST.CPU_Cores=" << XHOST.CPU_Cores << Message(__AFLOW_FILE__,aflags) << endl;
    aus << "DDDDD  aflags.AFLOW_GLOBAL_NCPUS=" << aflags.AFLOW_GLOBAL_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) { 	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
    //   kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;
    // }
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    // //CO20201220 X START
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    // //CO20201220 X STOP
    // //CO20220818 JHU_ROCKFISH START
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::JHU_ROCKFISH")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    // //CO20220818 JHU_ROCKFISH STOP
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;
    // }
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;
    // }
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;
    // }
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {	//CO
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_CPUS_ON_NODE;
    // }
    // //DX20190509 - MACHINE001 - START
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) {
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    // //DX20190509 - MACHINE001 - END
    // //DX20190509 - MACHINE002 - START
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) {
    //   if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    //DX20190509 - MACHINE002 - END
    //DX20190107 - CMU EULER - START
    // if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) {
    //  if(kflags.KBIN_MPI_NCPUS==0) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;
    // }
    //DX20180502 - CMU EULER - END

    // for reducint CPUs on the fly
    if(aflags.AFLOW_GLOBAL_NCPUS<0) kflags.KBIN_MPI_NCPUS=-aflags.AFLOW_GLOBAL_NCPUS; // this to force things on reducing CPUS

    aus << "DDDDD  kflags.KBIN_MPI_NCPUS=" << kflags.KBIN_MPI_NCPUS << Message(__AFLOW_FILE__,aflags) << endl;

    // for for LS coupling
    if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option) {
      if(!aurostd::substring2bool(kflags.KBIN_BIN,VASPLS_BIN_POSTFIX_DEFAULT)) kflags.KBIN_BIN=kflags.KBIN_BIN+VASPLS_BIN_POSTFIX_DEFAULT; // standard LS
      if(!aurostd::substring2bool(kflags.KBIN_MPI_BIN,VASPLS_BIN_POSTFIX_DEFAULT)) kflags.KBIN_MPI_BIN=kflags.KBIN_MPI_BIN+VASPLS_BIN_POSTFIX_DEFAULT; // standard LS
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
      aus << "00000  MESSAGE SPIN-ORBIT TYPE CALCULATIONS , adding " << VASPLS_BIN_POSTFIX_DEFAULT << " to BIN" << Message(__AFLOW_FILE__,aflags) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(kflags.KBIN_MPI) kflags.KBIN_BIN=kflags.KBIN_MPI_BIN; // forcing, no matter what

    uint xvasp_aopts_vxscheme_size=0;
    uint vflags_KBIN_VASP_FORCE_OPTION_IGNORE_AFIX_vxscheme_size=0;

    while(vasp_start) {
      // ********* RUN VASP                
      { // ERRORS
        bool error=FALSE;
        if(aurostd::FileEmpty(xvasp.Directory+"/INCAR"))   {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR "+function+": Empty INCAR"+Message(__AFLOW_FILE__,aflags));error=TRUE;return FALSE;}
        if(aurostd::FileEmpty(xvasp.Directory+"/POSCAR"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR "+function+": Empty POSCAR"+Message(__AFLOW_FILE__,aflags));error=TRUE;return FALSE;}
        if(aurostd::FileEmpty(xvasp.Directory+"/KPOINTS")) {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR "+function+": Empty KPOINTS"+Message(__AFLOW_FILE__,aflags));error=TRUE;return FALSE;}
        if(aurostd::FileEmpty(xvasp.Directory+"/POTCAR"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR "+function+": Empty POTCAR"+Message(__AFLOW_FILE__,aflags));error=TRUE;return FALSE;}
        if(error) return FALSE;

        if(LDEBUG){aus << __AFLOW_FUNC__ << " [1]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

        // FIX INCAR if alternating
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")) {
          if(aurostd::_isodd(xvasp.NRELAXING))  aus << "00000  MESSAGE Alternating: RELAX_CELL_VOLUME" << Message(__AFLOW_FILE__,aflags) << endl;
          if(aurostd::_iseven(xvasp.NRELAXING)) aus << "00000  MESSAGE Alternating: RELAX_IONS" << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC",FALSE);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL",FALSE);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS",FALSE);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE",FALSE);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME",FALSE);
          vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS_CELL_VOLUME");
          KBIN::XVASP_INCAR_Relax_ON(xvasp,vflags,xvasp.NRELAXING);
        }

        //CO20210315
        //print out these schemes so they can be picked up by the vasp monitor
        xvasp_aopts_vxscheme_size=xvasp.aopts.vxscheme.size();
        for(uint i=0;i<xvasp_aopts_vxscheme_size;i++){
          const string& flag=xvasp.aopts.vxscheme[i];
          if(flag.find("FLAG::")!=string::npos && flag.find("_PRESERVED")!=string::npos){
            if(xvasp.aopts.flag(flag)){aus << "MMMMM  MESSAGE xvasp.aopts.flag(\"" << flag << "\")=" << xvasp.aopts.flag(flag) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
          }
        }
        vflags_KBIN_VASP_FORCE_OPTION_IGNORE_AFIX_vxscheme_size=vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme.size();
        for(uint i=0;i<vflags_KBIN_VASP_FORCE_OPTION_IGNORE_AFIX_vxscheme_size;i++){
          const string& flag=vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme[i];
          if(vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(flag)){aus << "MMMMM  MESSAGE vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(\"" << flag << "\")=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.flag(flag) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
        }
        if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved){aus << "MMMMM  MESSAGE vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved=" << vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

        // RUN VASP NON QUEUE ------------------------------------------------------------------------
        if(kflags.KBIN_QSUB==FALSE) {
          nrun++;
          aus_exec << "cd " << xvasp.Directory << endl;
          aus_exec << "rm -f " << DEFAULT_VASP_OUT << endl;
          if(kflags.KBIN_MPI==FALSE) {
            aus << "00000  MESSAGE SERIAL job - [" << xvasp.str.atoms.size() << "atoms]" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            aus_exec << kflags.KBIN_BIN << " > " << DEFAULT_VASP_OUT << endl;
            aus << "00000  MESSAGE" << VASP_KEYWORD_EXECUTION << kflags.KBIN_BIN << " > " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl;  //CO20170628 - SLOW WITH MEMORY
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            aurostd::execute(aus_exec);
            aurostd::Sleep(_KVASP_VASP_SLEEP_);
          } else {
            aus << "00000  MESSAGE MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl;
            if(kflags.KBIN_MPI_OPTIONS!="") aus << "00000  MESSAGE MPI OPTIONS=[" << kflags.KBIN_MPI_OPTIONS << "]" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
            if(LDEBUG){
              aus << __AFLOW_FUNC__ << " aflags.AFLOW_MACHINE_GLOBAL=" << aflags.AFLOW_MACHINE_GLOBAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl;
              aus << __AFLOW_FUNC__ << " aflags.AFLOW_MACHINE_LOCAL=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            }
            // NO HOST ------------------------------------------------------------------------
            if(!aflags.AFLOW_MACHINE_LOCAL.flag()) {
              aus << "00000  MESSAGE" << VASP_KEYWORD_EXECUTION;
              if(!kflags.KBIN_MPI_OPTIONS.empty()){
                aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
                //[CO20210315 - do not print, will confuse vasp monitor]aus << kflags.KBIN_MPI_OPTIONS << "; ";
              }
              if(!kflags.KBIN_MPI_START.empty()){
                aus_exec << kflags.KBIN_MPI_START << " > " << DEFAULT_VASP_OUT << endl;
                //[CO20210315 - do not print, will confuse vasp monitor]aus << kflags.KBIN_MPI_START << " > " << DEFAULT_VASP_OUT << "; ";
              }
              aus_exec << kflags.KBIN_MPI_COMMAND << " " << kflags.KBIN_MPI_NCPUS << " " << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aus << kflags.KBIN_MPI_COMMAND << " " << kflags.KBIN_MPI_NCPUS << " " << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT; //[CO20210315 - do not print, will confuse vasp monitor]<< "; ";
              if(!kflags.KBIN_MPI_STOP.empty()){
                aus_exec << kflags.KBIN_MPI_STOP << " >> " << DEFAULT_VASP_OUT << endl;
                //[CO20210315 - do not print, will confuse vasp monitor]aus << kflags.KBIN_MPI_STOP << " >> " << DEFAULT_VASP_OUT;
              }
              aus << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_BETA_MPICH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_BETA_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_MPICH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_BETA_MPICH << endl;
              aus_exec << MPI_COMMAND_DUKE_BETA_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_MPICH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_BETA_OPENMPI ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) {
              if(!aurostd::substring2bool(kflags.KBIN_MPI_BIN,"_openmpi")) kflags.KBIN_MPI_BIN=kflags.KBIN_MPI_BIN+"_openmpi"; // fix the OPENMPI
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_BETA_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_BETA_OPENMPI << endl;
              aus_exec << MPI_COMMAND_DUKE_BETA_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_QRATS_MPICH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_QRATS_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QRATS_MPICH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_QRATS_MPICH << endl;
              aus_exec << MPI_COMMAND_DUKE_QRATS_MPICH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QRATS_MPICH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_QFLOW_OPENMPI ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_QFLOW_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << endl;
              aus_exec << MPI_COMMAND_DUKE_QFLOW_OPENMPI << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            //CO20201220 X START
            // HOST DUKE_X_X ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_X")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_X_X << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_X << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_X << endl;
              aus_exec << MPI_COMMAND_DUKE_X_X << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_X << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_CRAY")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_X_CRAY << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_CRAY << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_CRAY << endl;
              aus_exec << MPI_COMMAND_DUKE_X_CRAY << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_CRAY << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //        aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_OLDCRAY")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_X_OLDCRAY << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_OLDCRAY << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_OLDCRAY << endl;
              aus_exec << MPI_COMMAND_DUKE_X_OLDCRAY << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_OLDCRAY << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //        aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_X_SMB")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_X_SMB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_SMB << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_X_SMB << endl;
              aus_exec << MPI_COMMAND_DUKE_X_SMB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_X_SMB << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //        aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            //CO20201220 X STOP
            //CO20220818 JHU_ROCKFISH START
            // HOST JHU_ROCKFISH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::JHU_ROCKFISH")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_JHU_ROCKFISH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_JHU_ROCKFISH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_JHU_ROCKFISH << endl;
              aus_exec << MPI_COMMAND_JHU_ROCKFISH << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_JHU_ROCKFISH << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            //CO20220818 JHU_ROCKFISH STOP
            // HOST MPCDF_EOS_MPI ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) {
              // verbosization
              int local_NCPUS=kflags.KBIN_MPI_NCPUS;
              if(MPI_NCPUS_MPCDF_EOS>0) {
                local_NCPUS=MPI_NCPUS_MPCDF_EOS;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_EOS" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              // [OBSOLETE] the HT should come out from the ntasks and not intercepted anymore here
              // [OBSOLETE]	      if(MPI_HYPERTHREADING_MPCDF_EOS=="FALSE" || MPI_HYPERTHREADING_MPCDF_EOS=="OFF") {
              // [OBSOLETE]  local_NCPUS=local_NCPUS/2;
              // [OBSOLETE]	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = OFF" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE]	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_EOS=="TRUE" || MPI_HYPERTHREADING_MPCDF_EOS=="ON") {
              // [OBSOLETE]	local_NCPUS=local_NCPUS*2;
              // [OBSOLETE]	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = ON" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE]	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MPCDF_EOS << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_EOS << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_EOS << endl;
              aus_exec << MPI_COMMAND_MPCDF_EOS << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_EOS << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST MPCDF_DRACO_MPI ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) {
              // verbosization
              int local_NCPUS=kflags.KBIN_MPI_NCPUS;
              if(MPI_NCPUS_MPCDF_DRACO>0) {
                local_NCPUS=MPI_NCPUS_MPCDF_DRACO;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_DRACO" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              // [OBSOLETE] the HT should come out from the ntasks and not intercepted anymore here
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_DRACO=="FALSE" || MPI_HYPERTHREADING_MPCDF_DRACO=="OFF") {
              // [OBSOLETE] 	local_NCPUS=local_NCPUS/2;
              // [OBSOLETE] 	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = OFF" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE] 	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_DRACO=="TRUE" || MPI_HYPERTHREADING_MPCDF_DRACO=="ON") {
              // [OBSOLETE] 	local_NCPUS=local_NCPUS*2;
              // [OBSOLETE] 	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = ON" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE] 	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MPCDF_DRACO << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_DRACO << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_DRACO << endl;
              aus_exec << MPI_COMMAND_MPCDF_DRACO << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_DRACO << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST MPCDF_COBRA_MPI ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) {
              // verbosization 
              int local_NCPUS=kflags.KBIN_MPI_NCPUS;
              if(MPI_NCPUS_MPCDF_COBRA>0) {
                local_NCPUS=MPI_NCPUS_MPCDF_COBRA;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_COBRA" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              // [OBSOLETE] the HT should come out from the ntasks and not intercepted anymore here
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_COBRA=="FALSE" || MPI_HYPERTHREADING_MPCDF_COBRA=="OFF") {
              // [OBSOLETE] 	local_NCPUS=local_NCPUS/2;
              // [OBSOLETE] 	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = OFF" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE] 	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_COBRA=="TRUE" || MPI_HYPERTHREADING_MPCDF_COBRA=="ON") {
              // [OBSOLETE] 	local_NCPUS=local_NCPUS*2;
              // [OBSOLETE] 	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = ON" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE] 	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MPCDF_COBRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_COBRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_COBRA << endl;
              aus_exec << MPI_COMMAND_MPCDF_COBRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_COBRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            // HOST MPCDF_HYDRA_MPI ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) {
              // verbosization 
              int local_NCPUS=kflags.KBIN_MPI_NCPUS;
              if(MPI_NCPUS_MPCDF_HYDRA>0) {
                local_NCPUS=MPI_NCPUS_MPCDF_HYDRA;
                aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: kflags.KBIN_MPI_NCPUS=MPI_NCPUS_MPCDF_HYDRA" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
                aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              }
              // [OBSOLETE] the HT should come out from the ntasks and not intercepted anymore here
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_HYDRA=="FALSE" || MPI_HYPERTHREADING_MPCDF_HYDRA=="OFF") {
              // [OBSOLETE] 	local_NCPUS=local_NCPUS/2;
              // [OBSOLETE] 	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = OFF" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE] 	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              // [OBSOLETE] if(MPI_HYPERTHREADING_MPCDF_HYDRA=="TRUE" || MPI_HYPERTHREADING_MPCDF_HYDRA=="ON") {
              // [OBSOLETE] 	local_NCPUS=local_NCPUS*2;
              // [OBSOLETE] 	aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  Forcing: HYPERTHREADING = ON" << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              // [OBSOLETE] 	aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // [OBSOLETE] }
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << local_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              //	      aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MPCDF_HYDRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MPCDF_HYDRA << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name   // poe not MPI run
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MPCDF_HYDRA << endl;
              //	      aus_exec << MPI_COMMAND_MPCDF_HYDRA << " " << local_NCPUS << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aus_exec << MPI_COMMAND_MPCDF_HYDRA << " " << MPI_BINARY_DIR_MPCDF_HYDRA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;  // poe not MPI run
              //	    aurostd::PrintMessageStream(FileMESSAGE,aus_exec,XHOST.QUIET);
              aurostd::execute(aus_exec);
            }
            //DX20190509 - MACHINE001 - START
            // HOST MACHINE001_MPICH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MACHINE001 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE001 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE001 << endl;
              aus_exec << MPI_COMMAND_MACHINE001 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE001 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            //DX20190509 - MACHINE001 - END
            //DX20190509 - MACHINE002 - START
            // HOST MACHINE002_MPICH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MACHINE002 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE002 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE002 << endl;
              aus_exec << MPI_COMMAND_MACHINE002 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE002 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            //DX20190509 - MACHINE002 - END
            //DX20201005 - MACHINE003 - START
            // HOST MACHINE003_MPICH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE003")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MACHINE003 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE003 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE003 << endl;
              aus_exec << MPI_COMMAND_MACHINE003 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE003 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            //DX20201005 - MACHINE003 - END
            //DX20211011 - MACHINE004 - START
            // HOST MACHINE004_MPICH ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE004")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs  " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MACHINE004 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE004 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE004 << endl;
              aus_exec << MPI_COMMAND_MACHINE004 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE004 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            //DX20211011 - MACHINE004 - END
            // HOST DUKE_MATERIALS ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_MATERIALS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_MATERIALS << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_MATERIALS << endl;
              aus_exec << MPI_COMMAND_DUKE_MATERIALS << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_MATERIALS << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_AFLOWLIB ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_AFLOWLIB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_AFLOWLIB << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_AFLOWLIB << endl;
              aus_exec << MPI_COMMAND_DUKE_AFLOWLIB << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_AFLOWLIB << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST DUKE_HABANA ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_DUKE_HABANA << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_HABANA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_DUKE_HABANA << endl;
              aus_exec << MPI_COMMAND_DUKE_HABANA << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_DUKE_HABANA << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST FULTON_MARYLOU ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              //	      aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_FULTON_MARYLOU << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_FULTON_MARYLOU << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_FULTON_MARYLOU << endl;
              aus_exec << MPI_COMMAND_FULTON_MARYLOU << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              // aus_exec << MPI_COMMAND_FULTON_MARYLOU << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_FULTON_MARYLOU << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;  // with --np
              aurostd::execute(aus_exec);
            }
            // HOST CMU_EULER ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_CMU_EULER << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_CMU_EULER << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_CMU_EULER << endl;
              aus_exec << MPI_COMMAND_CMU_EULER << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_CMU_EULER << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;
              aurostd::execute(aus_exec);
            }
            // HOST OL ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MACHINE2 << " " << MPI_BINARY_DIR_MACHINE2 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE2 << endl;	//CO20181226
              aus_exec << MPI_COMMAND_MACHINE2 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE2 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;	//CO20181226 - adding kflags.KBIN_MPI_NCPUS
              aurostd::execute(aus_exec);
            }
            // HOST HOST1 ------------------------------------------------------------------------
            if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) {
              // verbosization
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << "  MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl; //HE20220309 use machine name
              aus << "00000  MESSAGE HOST=" << aflags.AFLOW_MACHINE_LOCAL.getattachedscheme("NAME") << " " << VASP_KEYWORD_EXECUTION << MPI_COMMAND_MACHINE1 << " " << MPI_BINARY_DIR_MACHINE1 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << Message(__AFLOW_FILE__,aflags,string(_AFLOW_MESSAGE_DEFAULTS_)+",memory") << endl; //HE20220309 use machine name  //CO20170628 - SLOW WITH MEMORY
              aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
              // run
              aus_exec << kflags.KBIN_MPI_OPTIONS << endl;
              aus_exec << MPI_OPTIONS_MACHINE1 << endl;	//CO20181226
              aus_exec << MPI_COMMAND_MACHINE1 << " " << kflags.KBIN_MPI_NCPUS << " " << MPI_BINARY_DIR_MACHINE1 << kflags.KBIN_MPI_BIN << " >> " << DEFAULT_VASP_OUT << endl;	//CO20181226 - adding kflags.KBIN_MPI_NCPUS
              aurostd::execute(aus_exec);
            }
            // DONE ------------------------------------------------------------------------
          }
          aurostd::Sleep(_KVASP_VASP_SLEEP_);
          vasp_start=FALSE;
        }
        // RUN VASP QUEUED ------------------------------------------------------------------------
        if(kflags.KBIN_QSUB) {
          nrun++;
          aus_exec << "cd " << xvasp.Directory << endl;
          if(kflags.KBIN_MPI==FALSE) {
            aus << "00000  MESSAGE QUEUED SERIAL job - [" << xvasp.str.atoms.size() << "atoms]" << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
            aurostd::stringstream2file(xvasp.xqsub.QSUB,string(xvasp.Directory+"/aflow.qsub.run"));
            aurostd::ChmodFile("755",string(xvasp.Directory+"/aflow.qsub.run"));
            aus_exec << kflags.KBIN_QSUB_COMMAND << " " << kflags.KBIN_QSUB_PARAMS << " " << "./aflow.qsub.run &" << endl;
            aurostd::execute(aus_exec);
            KBIN::QSUB_WaitFinished(aflags,FileMESSAGE,FALSE);
            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
          } else {
            aus << "00000  MESSAGE QUEUED MPI PARALLEL job - [" << xvasp.str.atoms.size() << "atoms] - " << " MPI=" << kflags.KBIN_MPI_NCPUS << "CPUs " << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
            aurostd::stringstream2file(xvasp.xqsub.QSUB,string(xvasp.Directory+"/aflow.qsub.run"));
            aurostd::ChmodFile("755",string(xvasp.Directory+"/aflow.qsub.run"));
            aus_exec << kflags.KBIN_QSUB_COMMAND << " " << kflags.KBIN_QSUB_PARAMS << " " << "./aflow.qsub.run &" << endl;
            aurostd::execute(aus_exec);
            KBIN::QSUB_WaitFinished(aflags,FileMESSAGE,FALSE);
            aurostd::RemoveFile(string(xvasp.Directory+"/aflow.qsub.done"));
          }	
          aurostd::Sleep(_KVASP_VASP_SLEEP_);
          vasp_start=FALSE;
        }
      }
      KBIN::WaitFinished(xvasp,aflags,FileMESSAGE,2,false);  //CO20201111 - try twice and NO verbose, we verbose in bigger VASP_Run() loop
      if(LDEBUG){aus << __AFLOW_FUNC__ << " [2]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

      if(aurostd::FileEmpty(xvasp.Directory+"/"+DEFAULT_VASP_OUT))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR "+function+": Empty "+DEFAULT_VASP_OUT+Message(__AFLOW_FILE__,aflags));return FALSE;}
      if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR"))  {KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR "+function+": Empty OUTCAR"+Message(__AFLOW_FILE__,aflags));return FALSE;}
      // DONT CHECK CONTCAR it can be empty
      // DONT CHECK OSZICAR it can be empty

      // update kpoints table

      // ***************** CHECK FOR ERRORS *********
      if(LDEBUG){
        aus << __AFLOW_FUNC__ << " [3a]  nrun=" << nrun << Message(__AFLOW_FILE__,aflags) << endl;
        aus << __AFLOW_FUNC__ << " [3b]  maxrun=" << maxrun << Message(__AFLOW_FILE__,aflags) << endl;
        cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }

      xmessage.clear(); //CO20210315

      // check VASP version
      string SVERSION=KBIN::OUTCAR2VASPVersionNumber(xvasp.Directory+"/OUTCAR"); //CO20210315
      double DVERSION=KBIN::VASPVersionString2Double(SVERSION); //CO20210315
      xmessage.push_attached("SVERSION",SVERSION);  //CO20210315 - put to xmessage
      xmessage.push_attached("DVERSION",aurostd::utype2string(DVERSION)); //CO20210315 - put to xmessage

      //get algo_current - START
      KBIN::VASP_Reread_INCAR(xvasp);
      string algo_current="NORMAL"; //vasp default: https://www.vasp.at/wiki/index.php/ALGO
      if(aurostd::substring2bool(xvasp.INCAR,"ALGO=",true)){algo_current=aurostd::toupper(aurostd::kvpair2string(xvasp.INCAR,"ALGO","="));}  //CO20210315 - remove whitespaces
      else if(aurostd::substring2bool(xvasp.INCAR,"IALGO=",true)){ //aflow also prints IALGO sometimes, need to check, remove whitespaces
        string algo_current_tmp=aurostd::toupper(KBIN::INCAR_IALGO2ALGO(aurostd::kvpair2utype<int>(xvasp.INCAR,"IALGO","=")));
        if(!algo_current_tmp.empty()){algo_current=algo_current_tmp;}
      }
      if(!algo_current.empty()){  //so we don't retry later as a fix
        xfixed.flag("ALGO="+algo_current,true);
        aus << "MMMMM  MESSAGE adding current \"ALGO=" << algo_current << "\" to xfixed" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      //get algo_current - END

      if(nrun<maxrun) {

        if(LDEBUG){aus << __AFLOW_FUNC__ << " [4]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

        VASP_ProcessWarnings(xvasp,aflags,kflags,xmessage,xwarning,FileMESSAGE);

        xfixed.flag("ALL",FALSE);
        vasp_start=FALSE;

        if(LDEBUG){aus << __AFLOW_FUNC__ << " [5]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

        KBIN::VASP_FixErrors(xvasp,xmessage,xwarning,xfixed,aflags,kflags,vflags,FileMESSAGE);

        //come back - should we check if any xwarning() flag is still on?
        //CO20210315 - NO, aflow patches if possible, otherwise let it run. we'll catch BIG issues with --xplug

        // ********* VASP TO BE RESTARTED *********
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " [DONE WITH CHECKS]" << Message(__AFLOW_FILE__,aflags) << endl;}
        if(xfixed.flag("ALL")) vasp_start=TRUE;
        if(vasp_start) {
          if(LDEBUG){aus << __AFLOW_FUNC__ << " [VASP TO BE RESTARTED]" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}
          aus << "00000  RESTART VASP" << Message(__AFLOW_FILE__,aflags) << endl;aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }
    }

    if(LDEBUG){aus << "MMMMM  MESSAGE tested all the errors" << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    bool Krun=TRUE;
    if(!aurostd::FileExist(xvasp.Directory+"/"+DEFAULT_VASP_OUT)) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file does not exist="+DEFAULT_VASP_OUT);}
    if(aurostd::FileEmpty(xvasp.Directory+"/"+DEFAULT_VASP_OUT)) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file empty="+DEFAULT_VASP_OUT);}
    if(!aurostd::FileExist(xvasp.Directory+"/OUTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file does not exist=OUTCAR");}
    if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file empty=OUTCAR");}
    if(!aurostd::FileExist(xvasp.Directory+"/CONTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file does not exist=CONTCAR");}
    if(aurostd::FileEmpty(xvasp.Directory+"/CONTCAR")) {Krun=FALSE;KBIN::VASP_Error(xvasp,"EEEEE  file is empty=CONTCAR");}

    if(LDEBUG) {aus << __AFLOW_FUNC__ << " Krun=" << Krun << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    if(LDEBUG) {aus << __AFLOW_FUNC__ << " END" << Message(__AFLOW_FILE__,aflags) << endl;cerr << aus.str();aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);}

    return Krun;
    // ********* FINISH
    //  return 1;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relax,bool qmwrite,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool Krun=TRUE;
    if(KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE)){Krun=(Krun&&true);}
    else{
      KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Run"+Message(__AFLOW_FILE__,aflags));
      Krun=false;
    }
    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string(relax));
    if(KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true)){Krun=(Krun&&true);} //CO20201111
    else{
      KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_RunFinished: OUTCAR is incomplete"+Message(__AFLOW_FILE__,aflags)); //CO20201111  //AFTER CONTCAR_SAVE_
      Krun=false;
    }
    KBIN::VASP_Backup(xvasp,qmwrite,relax);
    return Krun;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,string relaxB,bool qmwrite,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool Krun=TRUE;
    if(relaxA!=relaxB) {
      string message = "relaxA (" + relaxA + ") != relaxB (" + relaxB + ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    if(KBIN::VASP_Run(xvasp,aflags,kflags,vflags,FileMESSAGE)){Krun=(Krun&&true);}
    else{
      KBIN::VASP_Error(xvasp,"EEEEE  ERROR KBIN::VASP_Run"+Message(__AFLOW_FILE__,aflags));
      Krun=false;
    }
    if(_VASP_CONTCAR_SAVE_) KBIN::VASP_CONTCAR_Save(xvasp,string(relaxA));
    if(KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,true)){Krun=(Krun&&true);} //CO20201111
    else{
      KBIN::VASP_Error(xvasp,FileMESSAGE,"EEEEE  ERROR KBIN::VASP_RunFinished: OUTCAR is incomplete"+Message(__AFLOW_FILE__,aflags)); //CO20201111 //AFTER CONTCAR_SAVE_
      Krun=false;
    }
    KBIN::VASP_Backup(xvasp,qmwrite,relaxA);
    KBIN::VASP_Recycle(xvasp,relaxB);
    return Krun;
  }
} // namespace KBIN

namespace KBIN {
  bool VASP_RunFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool verbose) {
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    // if(verbose) aus << "00000  MESSAGE RUN CHECK FINISHED :" << Message(__AFLOW_FILE__,aflags) << endl;
    // if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // OUTCAR DOES NOT EXIST
    if(!aurostd::FileExist(xvasp.Directory+"/OUTCAR")) {
      if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR does not exist) :" << Message(__AFLOW_FILE__,aflags) << endl;
      if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return FALSE;
    }
    // OUTCAR EXISTS BUT EMPTY
    if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR")) {
      if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR is empty) :" << Message(__AFLOW_FILE__,aflags) << endl;
      if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return FALSE;
    }
    // OUTCAR EXISTS
    bool RemoveWS=true,case_insensitive=true,expect_near_end=true;
    if(aurostd::substring_present_file_FAST(xvasp.Directory+"/OUTCAR","Total CPU time used (sec)",RemoveWS,case_insensitive,expect_near_end)){  //CO20210601
      if(verbose) aus << "00000  MESSAGE RUN FINISHED (OUTCAR is complete) :" << Message(__AFLOW_FILE__,aflags) << endl;
      if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      return TRUE;
    }
    if(0){  //faster approach above
      ostringstream aus_exec;
      aurostd::StringstreamClean(aus_exec);
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat OUTCAR | grep CPU > aflow.check_outcar.tmp " << endl;
      aurostd::execute(aus_exec);
      aurostd::StringstreamClean(aus_exec);
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "rm -f aflow.check_outcar.tmp " << endl;
      //CO20201111 - not super efficient with cat/grep into file, but it's safe for the spaces, may change later
      if(aurostd::substring_present_file(xvasp.Directory+"/aflow.check_outcar.tmp",aurostd::RemoveWhiteSpaces("Total CPU time used (sec)"),TRUE)) {
        if(verbose) aus << "00000  MESSAGE RUN FINISHED (OUTCAR is complete) :" << Message(__AFLOW_FILE__,aflags) << endl;
        if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aurostd::execute(aus_exec);
        return TRUE;
      }
      aurostd::execute(aus_exec);
    }
    if(verbose) aus << "00000  MESSAGE RUN NOT FINISHED (OUTCAR is incomplete)" << Message(__AFLOW_FILE__,aflags) << endl;
    if(verbose) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    return FALSE;
  }
} // namespace KBIN

namespace KBIN {
  void WaitFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,uint max_count,bool verbose) {
    uint i=0;
    uint sleep_seconds=SECONDS_SLEEP_VASP_COMPLETION;
    if((max_count*sleep_seconds)>SECONDS_SLEEP_VASP_MONITOR){ //safety for --monitor_vasp
      sleep_seconds=(uint)(max(3.0,((double)SECONDS_SLEEP_VASP_MONITOR/(double)max_count)-5.0)); //the max ensures we don't go below 0 (if SECONDS_SLEEP_VASP_MONITOR is too low)
    }
    while((i++)<max_count && !KBIN::VASP_RunFinished(xvasp,aflags,FileMESSAGE,verbose)) {
      aurostd::Sleep(sleep_seconds); //CO20201111
    }
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Error(const _xvasp& xvasp,const string& message1,const string& message2,const string& message3) { //CO20210315 - cleaned up
    ofstream FileXVASP;
    return VASP_Error(xvasp,FileXVASP,message1,message2,message3);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Error(const _xvasp& xvasp,ofstream &FileMESSAGE,const string& message1,const string& message2,const string& message3) { //CO20210315 - cleaned up
    string FileNameXVASP=xvasp.Directory+"/"+DEFAULT_AFLOW_ERVASP_OUT;
    stringstream message;message << message1 << message2 << message3 << endl;
    aurostd::stringstream2file(message,FileNameXVASP,"APPEND");
    aurostd::PrintMessageStream(FileMESSAGE,message,XHOST.QUIET);
  }
} // namespace KBIN

namespace KBIN {
  string VASP_Analyze(_xvasp &xvasp,bool qmwrite) {       // AFLOW_FUNCTION_IMPLEMENTATION
    string function="KBIN::VASP_Analyze";
    // CHECK ERRORS
    bool error=FALSE;
    // if(aurostd::FileEmpty(xvasp.Directory+"/EIGENVAL")) {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty EIGENVAL");error=TRUE;}
    // if(aurostd::FileEmpty(xvasp.Directory+"/CHG"))      {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty CHG");error=TRUE;}
    // if(aurostd::FileEmpty(xvasp.Directory+"/CHGCAR"))   {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty CHGCAR");error=TRUE;}
    // if(aurostd::FileEmpty(xvasp.Directory+"/DOSCAR"))   {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty DOSCAR");error=TRUE;}
    if(aurostd::FileEmpty(xvasp.Directory+"/CONTCAR"))  {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty CONTCAR");error=TRUE;}
    if(aurostd::FileEmpty(xvasp.Directory+"/OUTCAR") || !aurostd::substring_present_file_FAST(xvasp.Directory+"/OUTCAR","Total CPU time used (sec)",true,true,true))   {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty OUTCAR");error=TRUE;} // SD20221019 - OUTCAR does not exist or is incomplete
    if(aurostd::FileEmpty(xvasp.Directory+"/INCAR"))    {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty INCAR");error=TRUE;}
    if(aurostd::FileEmpty(xvasp.Directory+"/vasprun.xml"))    {KBIN::VASP_Error(xvasp,"EEEEE  ERROR "+function+": Empty vasprun.xml");error=TRUE;}
    if(error) return "";

    // cerr << "# KBIN::VASP_Analyze BEGIN" << endl;
    xvasp.str.qm_clear();
    xvasp.str.qm_load(xvasp.Directory);
    // LOAD OUTPUTS
    stringstream strstream;
    strstream.clear();
    strstream.str(std::string());
    strstream.setf(std::ios::fixed,std::ios::floatfield);
    // OUTCAR OPERATIONS ---------------------------------------------------------------
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream << "[KBIN_ANALYZE]START_" << xvasp.AnalyzeLabel << endl;
    if(xvasp.AnalyzeLabel!="dielectric_static" && xvasp.AnalyzeLabel!="dielectric_dynamic") {
      strstream << "[AFLOW] **************************************************************************************************************************" << endl;
      strstream << "# POSITION                                       TOTAL-FORCE (eV/Angst)               " << endl;
      strstream << "[AFLOW] **************************************************************************************************************************" << endl;
      strstream.precision(_DOUBLE_WRITE_PRECISION_);  //CO20200731 - 12
      for(uint i=0;i<xvasp.str.atoms.size();i++) {   // clear        (from the previous step)
        for(uint j=1;j<=3;j++) {if(abs(xvasp.str.qm_positions.at(i)[j])<10.0) strstream << " ";if(xvasp.str.qm_positions.at(i)[j]>=0.0) strstream << " "; strstream << "   " << xvasp.str.qm_positions.at(i)[j] << " ";}
        for(uint j=1;j<=3;j++) {if(abs(xvasp.str.qm_forces.at(i)[j])<10.0) strstream << " ";if(xvasp.str.qm_forces.at(i)[j]>=0.0) strstream << " "; strstream << "   " << xvasp.str.qm_forces.at(i)[j] << " ";}
        strstream << endl;
      }
      strstream << "[AFLOW] **************************************************************************************************************************" << endl;
      // OUTCAR OPERATIONS ---------------------------------------------------------------
      strstream.setf(std::ios::scientific,std::ios::floatfield);
      strstream.setf(std::ios::left,std::ios::adjustfield);
      strstream << "E_cell=" << xvasp.str.qm_E_cell << "  (eV/cell)" << endl; 
      strstream << "E_atom=" << xvasp.str.qm_E_atom << "  (eV/at)" << endl;
      strstream << "H_cell=" << xvasp.str.qm_H_cell << "  (eV/cell)" << endl; 
      strstream << "H_atom=" << xvasp.str.qm_H_atom << "  (eV/at)" << endl;
      strstream << "PV_cell=" << xvasp.str.qm_PV_cell << "  (eV/cell)" << endl; 
      strstream << "PV_atom=" << xvasp.str.qm_PV_atom << "  (eV/at)" << endl;
      strstream << "mag_cell=" << xvasp.str.qm_mag_cell << "  (mu/cell)" << endl;
      strstream << "mag_atom="<< xvasp.str.qm_mag_atom << "  (mu/at)" << endl;
      xstructure qm_str(xvasp.str);    // suck it in !
      // qm_str=xvasp.str;
      qm_str.qm_recycle();
      strstream << "[AFLOW] **************************************************************************************************************************" << endl;
      strstream << "[VASP_POSCAR_MODE_EXPLICIT]START " << endl;
      strstream << qm_str;
      strstream << "[VASP_POSCAR_MODE_EXPLICIT]STOP " << endl;
    }

    //  if(aurostd::substring2bool(aurostd::execute2string("grep LEPSILON "+xvasp.Directory+"/OUTCAR"),"LEPSILON=T",TRUE))
    if(xvasp.AnalyzeLabel=="dielectric_static") {
      strstream << "[AFLOW] **************************************************************************************************************************" << endl;
      strstream << "[KBIN_ANALYZE]START_DIELECTRIC_STATIC" << endl;
      vector<string> vlines,tokens;
      aurostd::string2vectorstring(aurostd::execute2string("grep -A 4 \"MACROSCOPIC STATIC DIELECTRIC TENSOR\" "+xvasp.Directory+"/OUTCAR  | tail -n 3"),vlines);
      xmatrix<double> epsilon(3,3);
      if(vlines.size()==3) {
        for(uint i=1;i<=3;i++) {
          aurostd::string2tokens(vlines.at(i-1),tokens," ");
          if(tokens.size()==3) {
            for(uint j=1;j<=3;j++)
              epsilon(i,j)=aurostd::string2utype<double>(tokens.at(j-1));
          }
        }
      }
      strstream << " epsilon = " << endl;
      strstream << "  " << epsilon(1,1) << "  " << epsilon(1,2) << "  " << epsilon(1,3) << "  " << endl;
      strstream << "  " << epsilon(2,1) << "  " << epsilon(2,2) << "  " << epsilon(2,3) << "  " << endl;
      strstream << "  " << epsilon(3,1) << "  " << epsilon(3,2) << "  " << epsilon(3,3) << "  " << endl;
      strstream << "[KBIN_ANALYZE]STOP_DIELECTRIC_STATIC" << endl;
    }
    //  if(aurostd::substring2bool(aurostd::execute2string("grep LOPTICS "+xvasp.Directory+"/OUTCAR"),"LOPTICS=T",TRUE)) {  //[CO20200106 - close bracket for indenting]}
    if(xvasp.AnalyzeLabel=="dielectric_dynamic") {
      strstream << "[AFLOW] **************************************************************************************************************************" << endl;
      strstream << " DIELECTRIC DYNAMIC " << endl;  
      vector<string> vlines,tokens;string line;
      vector<xvector<double> > vepsilonIMAG,vepsilonREAL;
      aurostd::file2vectorstring(xvasp.Directory+"/OUTCAR",vlines);
      uint IMAG_start=0,IMAG_end=0,REAL_start=0,REAL_end=0;
      for(uint i=0;i<vlines.size();i++) {
        if(aurostd::substring2bool(vlines.at(i),"IMAGINARY DIELECTRIC FUNCTION")) IMAG_start=i;
        if(aurostd::substring2bool(vlines.at(i),"REAL DIELECTRIC FUNCTION")) {REAL_start=i;IMAG_end=i-2;REAL_end=(IMAG_end-IMAG_start)+REAL_start;}
      }
      // WRITING DIELECTRIC_IMAGINARY
      strstream << "[KBIN_ANALYZE]START_DIELECTRIC_IMAGINARY" << endl;
      for(uint i=IMAG_start;i<=IMAG_end;i++) {
        line=vlines.at(i);
        aurostd::StringSubst(line,"-0.000000","0");aurostd::StringSubst(line,"0.000000","0");aurostd::StringSubst(line,"\t"," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");
        strstream << line << endl;
        if(i>=IMAG_start+3) {
          xvector<double> epsilonIMAG(7);
          aurostd::string2tokens(vlines.at(i),tokens," ");
          if(tokens.size()==7) for(uint j=0;j<tokens.size();j++) epsilonIMAG(j+1)=aurostd::string2utype<double>(tokens.at(j));
          vepsilonIMAG.push_back(epsilonIMAG);
        }
      }
      strstream << "[KBIN_ANALYZE]STOPT_DIELECTRIC_IMAGINARY" << endl;
      // WRITING DIELECTRIC_REAL
      strstream << "[KBIN_ANALYZE]START_DIELECTRIC_REAL" << endl;
      for(uint i=REAL_start;i<=REAL_end;i++) {
        line=vlines.at(i);
        aurostd::StringSubst(line,"-0.000000","0");aurostd::StringSubst(line,"0.000000","0");aurostd::StringSubst(line,"\t"," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");aurostd::StringSubst(line,"  "," ");
        strstream << line << endl;
        if(i>=REAL_start+3) {
          xvector<double> epsilonREAL(7);
          aurostd::string2tokens(vlines.at(i),tokens," ");
          if(tokens.size()==7) for(uint j=0;j<tokens.size();j++) epsilonREAL(j+1)=aurostd::string2utype<double>(tokens.at(j));
          vepsilonREAL.push_back(epsilonREAL);
        }
      }
      strstream << "[KBIN_ANALYZE]STOPT_DIELECTRIC_REAL" << endl;
    }  
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;
    strstream << "[KBIN_ANALYZE]STOP_" << xvasp.AnalyzeLabel << endl;
    strstream << "[AFLOW] **************************************************************************************************************************" << endl;

    if(qmwrite) {
      string FileNameXVASP=xvasp.Directory+"/"+DEFAULT_AFLOW_QMVASP_OUT;
      stringstream FileXVASPout;
      //ME20200304 - do not overwrite prior runs
      string FileNameXVASPfull = "";
      if (aurostd::EFileExist(FileNameXVASP, FileNameXVASPfull)) aurostd::UncompressFile(FileNameXVASPfull);
      if(aurostd::FileExist(FileNameXVASP)) { //RECYCLE PREVIOUS STUFF
        stringstream FileXVASPin;
        aurostd::file2stringstream(FileNameXVASP, FileXVASPin);
        FileXVASPout << FileXVASPin.str();
      }
      FileXVASPout << strstream.str();
      aurostd::stringstream2file(FileXVASPout,xvasp.Directory+"/"+DEFAULT_AFLOW_QMVASP_OUT);
    }
    xvasp.str.qm_calculated=TRUE;
    //  cerr << "# KBIN::VASP_Analyze END" << endl;
    return strstream.str();
  }
}  // namespace KBIN

namespace KBIN {
  void GenerateAflowinFromVASPDirectory(_aflags &aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ifstream FileSUBDIR;string FileNameSUBDIR;
    FileNameSUBDIR=aflags.Directory;
    FileSUBDIR.open(FileNameSUBDIR.c_str(),std::ios::in);
    FileSUBDIR.clear();FileSUBDIR.close();
    ostringstream aus;

    if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;

    if(!FileSUBDIR) {                                                                                           // ******* Directory is non existent
      aus << "Directory not found";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _FILE_NOT_FOUND_);
    } else {                                                                                                    // ******* Directory EXISTS
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
      if(aurostd::FileExist(aflags.Directory+"/"+_AFLOWLOCK_) || aurostd::EFileExist(aflags.Directory+"/"+_AFLOWLOCK_))	{ // ******* Directory is locked
        // LOCK exist, then RUN already RUN
        aus << "Directory LOCKED";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
      }
      if(aurostd::FileExist(aflags.Directory+"/SKIP") || aurostd::EFileExist(aflags.Directory+"/SKIP")) {	// ******* Directory is skipped
        // SKIP exist, then RUN already RUN
        aus << "Directory SKIPPED";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, aus.str(), _RUNTIME_ERROR_);
      }

      // ******* Directory is un locked/skipped
      /// ******************************************************************
      // RESET LOCK
      ofstream FileLOCK;
      string FileNameLOCK=aflags.Directory+"/"+_AFLOWLOCK_;
      FileLOCK.open(FileNameLOCK.c_str(),std::ios::out);
      /// ******************************************************************
      // CHECK FOR INCAR KPOINTS POSCAR POTCAR
      ifstream FileINCAR;string FileNameINCAR=aflags.Directory+"/INCAR";FileINCAR.open(FileNameINCAR.c_str(),std::ios::in);
      if(!FileINCAR)  {                                                                                        // ******* INCAR does not exist
        aus << "EEEEE  INCAR ABSENT  =" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
      }
      ifstream FileKPOINTS;string FileNameKPOINTS=aflags.Directory+"/KPOINTS";FileKPOINTS.open(FileNameKPOINTS.c_str(),std::ios::in);
      if(!FileKPOINTS)  {                                                                                        // ******* KPOINTS does not exist
        aus << "EEEEE  KPOINTS ABSENT  =" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
      }
      ifstream FilePOSCAR;string FileNamePOSCAR=aflags.Directory+"/POSCAR";FilePOSCAR.open(FileNamePOSCAR.c_str(),std::ios::in);
      if(!FilePOSCAR)  {                                                                                        // ******* POSCAR does not exist
        aus << "EEEEE  POSCAR ABSENT  =" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
      }
      ifstream FilePOTCAR;string FileNamePOTCAR=aflags.Directory+"/POTCAR";FilePOTCAR.open(FileNamePOTCAR.c_str(),std::ios::in);
      if(!FilePOTCAR)  {                                                                                        // ******* POTCAR does not exist
        aus << "EEEEE  POTCAR ABSENT  =" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintErrorStream(FileLOCK,aus,XHOST.QUIET);
      }
      // ----------------------------------------------------------------------------------------------------
      if(FileINCAR && FileKPOINTS && FilePOSCAR && FilePOTCAR) {
        // VASP INCAR KPOINTS POSCAR POTCAR ARE PRESENT
        /// ******************************************************************
        // WRITE LOCK
        aus << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow" << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "MMMMM  (C) " << XHOST.Copyright_Years << ", Stefano Curtarolo - Duke University " << Message(__AFLOW_FILE__,aflags) << endl;
        aus << "MMMMM  High-Throughput ab-initio Computing" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
        aus << "00000  MESSAGE GENERATING " << _AFLOWIN_ << " from VASP-xCARs files" << Message(__AFLOW_FILE__,aflags) << endl;
        aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
        /// ******************************************************************
        // RESET AFLOWIN
        ofstream FileAFLOWIN;
        string FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
        FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::out);
        /// ******************************************************************
        // WRITE AFLOWIN
        // WRITE TITLE
        string str1,str2;
        getline(FileINCAR,str1);
        FileINCAR.seekg(0);
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[AFLOW] Automatically generated from XCARS by aflow/aflowd " << string(AFLOW_VERSION) << endl;
        FileAFLOWIN << "[AFLOW] Automatic-Flow" << Message(__AFLOW_FILE__,aflags) << endl;
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
        for(int i=0;i<(int) XHOST.argv.size()-1;i++) {
          str1=XHOST.argv.at(i);
          str2=XHOST.argv.at(i+1);
          if(str1=="--set" && str2.at(0)=='[') {
            aus << "00000  MESSAGE Adding " << str2 << " to " << _AFLOWIN_ << Message(__AFLOW_FILE__,aflags) << endl;
            aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
            FileAFLOWIN << str2 << endl;
          }
        }
        // WRITE INCAR
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
        while (getline(FileINCAR,str1)) FileAFLOWIN << "[VASP_INCAR_FILE]" << str1 << endl;
        // WRITE KPOINTS
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
        while (getline(FileKPOINTS,str1)) FileAFLOWIN << "[VASP_KPOINTS_FILE]" << str1 << endl;
        // WRITE POSCAR
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
        while (getline(FilePOSCAR,str1)) FileAFLOWIN << "[VASP_POSCAR_FILE]" << str1 << endl;
        // WRITE POTCAR
        FileAFLOWIN << AFLOWIN_SEPARATION_LINE << endl;
        FileAFLOWIN << "[VASP_POTCAR_MODE_EXPLICIT]" << endl;
        while (getline(FilePOTCAR,str1)) FileAFLOWIN << str1 << endl;
        // close everything.
        FileINCAR.clear();FileINCAR.close();
        FileKPOINTS.clear();FileKPOINTS.close();
        FilePOSCAR.clear();FilePOSCAR.close();
        FilePOTCAR.clear();FilePOTCAR.close();
        FileAFLOWIN.flush();FileAFLOWIN.clear();FileAFLOWIN.close();
        /// ******************************************************************
        // everything is done. check if we nned to delete VASP FILES
        if(aurostd::args2flag(XHOST.argv,"--delete_xcars")) {
          aus << "00000  MESSAGE Removing vasp files in" << Message(__AFLOW_FILE__,aflags) << endl;
          aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
          aus << "cd " << aflags.Directory << endl;
          aus << "rm -f `ls | grep -v " << _AFLOWIN_ << " | grep -v LOCK ` " << endl;
          aurostd::execute(aus);
        }
      }
      FileLOCK.flush();FileLOCK.clear();FileLOCK.close();
    }
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Backup(_xvasp& xvasp,bool qmwrite,const string& ext) {        // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210315
    xvasp.AnalyzeLabel=ext;
    KBIN::VASP_Analyze(xvasp,qmwrite);
    ostringstream aus;

    for(uint iext=0;iext<XHOST.vext.size();iext++) { 
      if(aurostd::FileExist(xvasp.Directory+"/core"+XHOST.vext[iext]))
        aurostd::execute("rm -f "+xvasp.Directory+"/core"+XHOST.vext[iext]);
    }
    if(!xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVECAR")) aurostd::RemoveFile(xvasp.Directory+"/WAVECAR");
    if(!xvasp.aopts.flag("FLAG::WAVEDER_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVEDER")) aurostd::RemoveFile(xvasp.Directory+"/WAVEDER");
    if(aurostd::FileExist(xvasp.Directory+"/aflow.qsub.run")) aurostd::RemoveFile(xvasp.Directory+"/aflow.qsub.run");
    if(aurostd::FileExist(xvasp.Directory+"/aflow.qsub.out")) aurostd::RemoveFile(xvasp.Directory+"/aflow.qsub.out");
    if(aurostd::FileExist(xvasp.Directory+"/AECCAR0")) aurostd::file2file(xvasp.Directory+"/AECCAR0",xvasp.Directory+"/AECCAR0."+ext);  // BADER
    if(aurostd::FileExist(xvasp.Directory+"/AECCAR1")) aurostd::file2file(xvasp.Directory+"/AECCAR1",xvasp.Directory+"/AECCAR1."+ext);  // BADER
    if(aurostd::FileExist(xvasp.Directory+"/AECCAR2")) aurostd::file2file(xvasp.Directory+"/AECCAR2",xvasp.Directory+"/AECCAR2."+ext);  // BADER
    if(aurostd::FileExist(xvasp.Directory+"/CHG")) aurostd::file2file(xvasp.Directory+"/CHG",xvasp.Directory+"/CHG."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/CHGCAR")) aurostd::file2file(xvasp.Directory+"/CHGCAR",xvasp.Directory+"/CHGCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/CONTCAR")) aurostd::file2file(xvasp.Directory+"/CONTCAR",xvasp.Directory+"/CONTCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/DYNMAT")) aurostd::file2file(xvasp.Directory+"/DYNMAT",xvasp.Directory+"/DYNMAT."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/DOSCAR")) aurostd::file2file(xvasp.Directory+"/DOSCAR",xvasp.Directory+"/DOSCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/ELFCAR")) aurostd::file2file(xvasp.Directory+"/ELFCAR",xvasp.Directory+"/ELFCAR."+ext);  // ELF
    if(aurostd::FileExist(xvasp.Directory+"/EIGENVAL")) aurostd::file2file(xvasp.Directory+"/EIGENVAL",xvasp.Directory+"/EIGENVAL."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/IBZKPT")) aurostd::file2file(xvasp.Directory+"/IBZKPT",xvasp.Directory+"/IBZKPT."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/INCAR")) aurostd::file2file(xvasp.Directory+"/INCAR",xvasp.Directory+"/INCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/KPOINTS")) aurostd::file2file(xvasp.Directory+"/KPOINTS",xvasp.Directory+"/KPOINTS."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/OSZICAR")) aurostd::file2file(xvasp.Directory+"/OSZICAR",xvasp.Directory+"/OSZICAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/OUTCAR")) aurostd::file2file(xvasp.Directory+"/OUTCAR",xvasp.Directory+"/OUTCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/PCDAT")) aurostd::file2file(xvasp.Directory+"/PCDAT",xvasp.Directory+"/PCDAT."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/POSCAR")) aurostd::file2file(xvasp.Directory+"/POSCAR",xvasp.Directory+"/POSCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/POTCAR")) aurostd::file2file(xvasp.Directory+"/POTCAR",xvasp.Directory+"/POTCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/PROCAR")) aurostd::file2file(xvasp.Directory+"/PROCAR",xvasp.Directory+"/PROCAR."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/XDATCAR")) aurostd::file2file(xvasp.Directory+"/XDATCAR",xvasp.Directory+"/XDATCAR."+ext);
    if(xvasp.aopts.flag("FLAG::WAVECAR_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVECAR")) aurostd::file2file(xvasp.Directory+"/WAVECAR",xvasp.Directory+"/WAVECAR."+ext);
    if(xvasp.aopts.flag("FLAG::WAVEDER_PRESERVED") && aurostd::FileExist(xvasp.Directory+"/WAVEDER")) aurostd::file2file(xvasp.Directory+"/WAVEDER",xvasp.Directory+"/WAVEDER."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/"+DEFAULT_VASP_OUT)) aurostd::file2file(xvasp.Directory+"/"+DEFAULT_VASP_OUT,xvasp.Directory+"/"+DEFAULT_VASP_OUT+"."+ext);
    if(aurostd::FileExist(xvasp.Directory+"/vasprun.xml")) aurostd::file2file(xvasp.Directory+"/vasprun.xml",xvasp.Directory+"/vasprun.xml."+ext);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_CONTCAR_Save(const _xvasp& xvasp,const string& ext) {        // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    return VASP_CONTCAR_Save(xvasp.Directory,ext);  //CO20210716
  }
  void VASP_CONTCAR_Save(const string& directory,const string& ext) {        // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    string function="KBIN::VASP_CONTCAR_Save";
    string operation=function+"("+ext+")";
    string aflowin=directory+"/"+_AFLOWIN_;
    string contcar=directory+string("/CONTCAR");
    if(aurostd::FileExist(aflowin) && aurostd::FileExist(contcar) && !aurostd::FileEmpty(contcar)) {
      xstructure xstr(contcar,IOAFLOW_AUTO);
      ostringstream aus;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      aus << "[AFLOW] SELF-MODIFICATION" << endl;
      aus << "[AFLOW] Recycling CONTCAR of " << ext << endl;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      aus << _VASP_POSCAR_MODE_EXPLICIT_START_ << endl;
      aus << xstr;
      aus << _VASP_POSCAR_MODE_EXPLICIT_STOP_ << endl;
      aus << AFLOWIN_SEPARATION_LINE << endl;
      KBIN::AFLOWIN_ADD(aflowin,aus,"");
      KBIN::AFLOWIN_REMOVE(aflowin,"[VASP_FORCE_OPTION]VOLUME",operation); //CO20210315
    }
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Recycle(const _xvasp& xvasp,const string& ext) {        // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210315
    aurostd::CopyFile(xvasp.Directory+"/CONTCAR."+ext,xvasp.Directory+"/POSCAR");
    aurostd::CopyFile(xvasp.Directory+"/INCAR."+ext,xvasp.Directory+"/INCAR");
    aurostd::CopyFile(xvasp.Directory+"/KPOINTS."+ext,xvasp.Directory+"/KPOINTS");
    aurostd::CopyFile(xvasp.Directory+"/POTCAR."+ext,xvasp.Directory+"/POTCAR");
  }
} // namespace KBIN

namespace KBIN {
  void VASP_Recycle(const _xvasp& xvasp,int relax_number) {        // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
      aurostd::execute(XHOST.vzip[iext]+" -dqf "+aurostd::CleanFileName(xvasp.Directory+"/*"+XHOST.vext[iext]));
    }
    aurostd::CopyFile(xvasp.Directory+"/CONTCAR.relax"+aurostd::utype2string<int>(relax_number),xvasp.Directory+"/POSCAR");
    aurostd::CopyFile(xvasp.Directory+"/INCAR.relax"+aurostd::utype2string<int>(relax_number),xvasp.Directory+"/INCAR");
    aurostd::CopyFile(xvasp.Directory+"/KPOINTS.relax"+aurostd::utype2string<int>(relax_number),xvasp.Directory+"/KPOINTS");
    aurostd::CopyFile(xvasp.Directory+"/POTCAR.relax"+aurostd::utype2string<int>(relax_number),xvasp.Directory+"/POTCAR");
  }
} // namespace KBIN

namespace KBIN {
  void VASP_RecycleExtraFile(const _xvasp& xvasp,const string& xfile,const string& relax) {        // AFLOW_FUNCTION_IMPLEMENTATION //CO20210315
    aurostd::CopyFile(xvasp.Directory+"/"+xfile+"."+relax,xvasp.Directory+"/"+xfile);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_RecycleExtraFile(const _xvasp& xvasp,const string& xfile,int relax_number) {        // AFLOW_FUNCTION_IMPLEMENTATION  //CO20210315
    for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed 
      aurostd::execute(XHOST.vzip[iext]+" -dqf "+aurostd::CleanFileName(xvasp.Directory+"/"+xfile+XHOST.vext[iext]));
    }
    aurostd::CopyFile(xvasp.Directory+"/"+xfile+".relax"+aurostd::utype2string<int>(relax_number),xvasp.Directory+"/"+xfile);
  }
} // namespace KBIN

namespace KBIN {
  void VASP_BackupOriginal(_xvasp xvasp) {        // AFLOW_FUNCTION_IMPLEMENTATION
    aurostd::CopyFile(xvasp.Directory+"/KPOINTS",xvasp.Directory+"/KPOINTS.orig");
    aurostd::CopyFile(xvasp.Directory+"/INCAR",xvasp.Directory+"/INCAR.orig");
    aurostd::CopyFile(xvasp.Directory+"/POSCAR",xvasp.Directory+"/POSCAR.orig");
  }
} // namespace KBIN

namespace KBIN {
  //SD20221026 - eliminated subshell
  //CO20210315 - reconsider rewriting these functions to eliminate subshell
  //NELM comes from xOUTCAR
  //NSTEPS comes from xOSZICAR (to be created)
  uint VASP_getNELM(const string& outcar){ //CO20200624
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    ifstream FileOUTCAR;
    FileOUTCAR.open(outcar.c_str(),std::ios::in);
    string tmp=aurostd::kvpair2string(FileOUTCAR,"NELM","=");
    FileOUTCAR.close();
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << outcar << " NELM grep response=\"" << tmp << "\"" << endl;}
    int NELM=60;  //VASP default
    if(!tmp.empty() && aurostd::isfloat(tmp)){NELM=aurostd::string2utype<int>(tmp);}
    return NELM;
  }
  uint VASP_getNSTEPS(const string& oszicar){  //CO20200624
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    ifstream FileOSZICAR;
    FileOSZICAR.open(oszicar.c_str(),std::ios::in);
    string tmp=aurostd::kvpair2string(FileOSZICAR,"DAV",":",-1);
    FileOSZICAR.close();
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " " << oszicar << " NSTEPS grep response=\"" << tmp << "\"" << endl;}
    int NSTEPS=0;  //VASP default
    if(!tmp.empty()){
      if(aurostd::isfloat(tmp)){NSTEPS=aurostd::string2utype<int>(tmp);}
      else if(tmp.find("*")!=string::npos){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " found number bigger than 999" << endl;}
        vector<string> vlines,tokens;
        aurostd::efile2vectorstring(oszicar,vlines);
        bool found_F_line=false;
        uint i=0,j=0,nsteps=0;
        for(i=vlines.size()-1;i<vlines.size()&&NSTEPS==0;i--){ //go backwards
          if(found_F_line){
            aurostd::string2tokens(vlines[i],tokens," ");
            if(tokens.size()<2){continue;}
            if(tokens[0].find(":")==string::npos){continue;}  //safety check, should be "DAV:" or "RMM:"
            tmp=tokens[1];
            if(aurostd::isfloat(tmp)){
              nsteps=aurostd::string2utype<uint>(tmp);
              if(LDEBUG){
                cerr << __AFLOW_FUNC__ << " last countable nstep: " << nsteps << endl;
                cerr << __AFLOW_FUNC__ << " steps after last counterable NSTEP: " << j << endl;
              }
              NSTEPS=nsteps+j;
              if(LDEBUG){cerr << __AFLOW_FUNC__ << " NSTEPS=" << NSTEPS << endl;}
            }
            else{j++;}
          }
          if(vlines[i].find("F=")!=string::npos){found_F_line=true;continue;}
        }
      }
    }
    return NSTEPS;
  }
  bool VASP_OSZICARUnconverging(const string& dir,uint cutoff) {  //CO20210601
    //this function will read the whole OSZICAR looking for electronic sloshing issues
    //if there are $cutoff ionic steps that did not converge electronically, return true
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    vector<string> vlines,vrelax,tokens;
    aurostd::file2vectorstring(dir+"/OSZICAR",vlines);
    uint i=0;
    for(i=1;i<vlines.size();i++){ //start at 1 so i-1 isn't a problem
      if(vlines[i].find("F=")!=string::npos){vrelax.push_back(vlines[i-1]);}
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " nionic=" << vrelax.size() << endl;}
    if(vrelax.size()<cutoff) return FALSE; // no problem
    // otherwise check for issues.
    uint NELM=KBIN::VASP_getNELM(dir+"/OUTCAR");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " NELM=" << NELM << endl;}
    uint nsteps=0,nissues=0;
    for(i=0;i<vrelax.size()&&nissues<cutoff;i++) {
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vrelax[i=" << i << "]=\"" << vrelax[i] << "\"" << endl;}
      aurostd::string2tokens(vrelax[i],tokens," ");
      if(tokens.size()>1){
        nsteps=aurostd::string2utype<uint>(tokens[1]);
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " nsteps=" << nsteps << endl;}
        if(nsteps!=0 && nsteps>=NELM){nissues++;}
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " nissues=" << nissues << endl;}
    if(nissues==cutoff) return true;
    return false;
  }
  bool VASP_OSZICARUnconverged(const string& oszicar,const string& outcar) {  //CO20210601
    //this function only looks at the last electronic SC step (different than VASP_OSZICARUnconverging, good for STATIC calcs)
    //if it is unconverged, return true
    bool LDEBUG=(FALSE || VERBOSE_MONITOR_VASP || _DEBUG_KVASP_ || XHOST.DEBUG);
    uint NELM=KBIN::VASP_getNELM(outcar);
    uint NSTEPS=KBIN::VASP_getNSTEPS(oszicar);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " NELM=" << NELM << endl;
      cerr << __AFLOW_FUNC__ << " NSTEPS=" << NSTEPS << endl;
    }
    if(NELM!=0 && NSTEPS!=0 && NSTEPS>=NELM){return true;}
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
    //[CO20191112 - OBSOLETE]int PATH_LENGTH_MAX = 1024 ;
    //[CO20191112 - OBSOLETE]char work_dir[PATH_LENGTH_MAX] ;
    string message = "";
    string outcarfile, outcarpath ;
    string outcarpath_tmp = aurostd::TmpFileCreate("OUTCARc1.tmp") ;
    vector<string> outcarlines, endline, startline, vasptoken ;
    xmatrix<double> statdiel(3,3), eigenvec(3,3) ;
    double eps = 1.0E-5 ; // need to define this more rigorously
    //[CO20191112 - OBSOLETE]getcwd(work_dir, PATH_LENGTH_MAX) ;
    string work_dir=aurostd::getPWD();  //CO20191112

    if(!aurostd::FileExist(outcar)) {
      message = "check filename || file missing";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    } else {
      outcarpath = "/" + outcar ;
      outcarpath = work_dir + outcarpath ;
      vector<string> outcardata ;
      aurostd::string2tokens(outcarpath, outcardata, ".") ;
      if(outcardata.at(outcardata.size()-1) == "bz2") { // compressed option
        aurostd::execute("bzcat " + outcarpath + " > " + outcarpath_tmp) ;
        aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
      } else if(outcardata.at(outcardata.size()-1) == "xz") { // compressed option
        aurostd::execute("xzcat " + outcarpath + " > " + outcarpath_tmp) ;
        aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
      } else if(outcardata.at(outcardata.size()-1) == "gz") { // compressed option
        aurostd::execute("gzcat " + outcarpath + " > " + outcarpath_tmp) ;
        aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
      } else { // plain text option
        aurostd::execute("cat " + outcarpath + " > " + outcarpath_tmp) ;
        aurostd::file2vectorstring(outcarpath_tmp, outcarlines) ;
      }
    }
    // check the loaded OUTCAR
    aurostd::string2tokens(outcarlines.at(0),startline," ");
    aurostd::string2tokens(startline.at(0),vasptoken,".");
    aurostd::string2tokens(outcarlines.at(outcarlines.size()-1),endline," ");
    if(vasptoken.at(0) != "vasp" || endline.at(0) != "Voluntary") { // first and last line check
      message =  "OUTCAR file is probably corrupt";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    uint sec_count = 0 ;
    for (uint ii=outcarlines.size()-12 ; ii<outcarlines.size() ; ii++) { // presence timing information check
      vector<string> timetoken ;
      aurostd::string2tokens(outcarlines.at(ii),timetoken," ") ;
      if(timetoken.size() > 0) {
        for (uint jj=0 ; jj<timetoken.size() ; jj++)
        { if(timetoken.at(jj) == "(sec):") sec_count+=1 ; }
      }
    }
    if(sec_count != 4) { // first and last line check
      message =  "OUTCAR file is probably corrupt";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    // OUTCAR is now in memory, now parse the info
    vector<string> words_line ;
    vector<string> vec1, vec2, vec3 ;
    bool  check_digit = false ;
    uint  refline = 0;
    for (uint ii=outcarlines.size()-1 ; ii > 1 ; ii--) { // line contents
      for (uint jj=0 ; jj < words_line.size() ; jj++) {
        string search_term = "MACROSCOPIC" ;
        string test_word = words_line.at(jj) ;
        if(test_word == search_term) { // start of dielectric tensor
          refline = ii + 2 ;
          check_digit = true ;
        }
      }
      if(check_digit) { // put the tensor info into the string vectors
        aurostd::string2tokens(outcarlines.at(refline+0),vec1," ") ;
        aurostd::string2tokens(outcarlines.at(refline+1),vec2," ") ;
        aurostd::string2tokens(outcarlines.at(refline+2),vec3," ") ;
        for (uint jj=1 ; jj <= 3 ; jj++) { // string to double, 3x3 matrix, be careful with array bounds
          statdiel(1,jj) = atof(vec1.at(jj-1).c_str()) ;
          statdiel(2,jj) = atof(vec2.at(jj-1).c_str()) ;
          statdiel(3,jj) = atof(vec3.at(jj-1).c_str()) ;
        }
        break ;
      }
    }
    if(!check_digit) {
      message = outcar + " lacks MACROSCOPIC statement";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    } // DONE PARSING //
    bool matcheck = false ;
    for (uint ii = 1 ; ii <= 3 ; ii++) { // clean up spuriously small values: e.g. "-0.000001"
      if( abs(statdiel(1,ii)) < eps ) statdiel(1,ii) = 0.00 ;
      if( abs(statdiel(2,ii)) < eps ) statdiel(2,ii) = 0.00 ;
      if( abs(statdiel(3,ii)) < eps ) statdiel(3,ii) = 0.00 ;
    }
    for (uint ii = 1 ; ii <= 3 ; ii++) { // check if it is asymmetric & if large off-diags exist
      for (uint jj = 1 ; jj <= 3 ; jj++) {
        double testdiff = statdiel[ii][jj] - statdiel[jj][ii] ;
        if(testdiff >= eps) { // eps is a bit arbitrary right now ..
          // serious issues with VASP calculation here: 
          message = "asymmetric dielectric tensor";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        } else { // only if small
          statdiel(ii,jj) = statdiel(jj,ii) ;
        }
        if(ii != jj) {
          if(abs(statdiel(ii,jj)) > 0 || abs(statdiel(jj,ii)) > 0) {
            matcheck = true ;
            break ;
          }
        }
      }
      if(matcheck) break ;
    }
    matcheck = true ;
    if(matcheck)
    { // diagonalize the 3x3 matrix
      aurostd::eigen(statdiel,eigr,eigi) ;
    }
  } // loop GetStatDiel
} // namespace KBIN

// ***************************************************************************

namespace KBIN {
  /// @brief gets a string containing the VASP version from binary
  ///
  /// @param binfile absolute path to the binary file
  ///
  /// @return string containing the VASP version
  ///
  /// @note To get the actual clean VASP version, one must still pass this function
  /// to VASPVersionString2Number or VASPVersionString2Double
  ///
  /// @authors
  /// @mod{SD,20220923,updated for vasp6.3}
  /// @mod{SD,20220401,created function based on getVASPVersionString}
  string BIN2VASPVersion(const string& binfile){
    ifstream infile(binfile.c_str(), std::ios::in | std::ios::binary);
    if (!infile.is_open()) {return "";}
    int bufferSize = 1024, i;
    char buffer[bufferSize];
    string vaspVersion = "", buffer_str = "";
    bool found_vasp = false;
    while (vaspVersion.empty() && !infile.eof()) {
      if (!infile.read(buffer, bufferSize)) {bufferSize = infile.gcount();}
      //SD20220401 - need to search for multiple keywords to find the correct line with the VASP version; this could change in the future with new VASP releases
      for (i = 0; !found_vasp && i < bufferSize - 5; i++) { // search for "vasp."  or "vasp(" in the buffer
        if ((buffer[i] == 'v') &&
            (buffer[i + 1] == 'a') &&
            (buffer[i + 2] == 's') &&
            (buffer[i + 3] == 'p') &&
            ((buffer[i + 4] == '.') || (buffer[i + 4] == '('))) {
          found_vasp = true;
          break;
        }
      }
      if (found_vasp) {
        buffer_str += std::string(reinterpret_cast<char*>(buffer), bufferSize); // avoid null-terminator
        if (i != 0) { // find("vasp.")
          buffer_str = buffer_str.substr(i); // get the buffer string starting from "vasp."
          i = 0;
        }
        buffer_str.replace(4, 1, "."); // replaces buffer string to always be "vasp."
        if (buffer_str.find("complex") != string::npos) { // second keyword
          buffer_str = buffer_str.substr(0, buffer_str.find("complex")); // get the buffer string up to "complex"
          if (buffer_str.find("\n") != string::npos || buffer_str.find("\\n") != string::npos) {
            buffer_str = "";
            found_vasp = false; // if there are newlines between the keywords then it is the wrong line
          }
          else { // no newlines, so this is the correct line
            vaspVersion = buffer_str.substr(0, buffer_str.find(" "));
            break;
          }
        }
        else if (buffer_str.find("scpc") != string::npos) { // second keyword
          buffer_str = buffer_str.substr(0, buffer_str.find("scpc")); // get the buffer string up to "scpc"
          if (buffer_str.find("\n") != string::npos || buffer_str.find("\\n") != string::npos) {
            buffer_str = "";
            found_vasp = false; // if there are newlines between the keywords then it is the wrong line
          }
          else { // no newlines, so this is the correct line
            vaspVersion = buffer_str.substr(0, buffer_str.find(" "));
            break;
          }
        }
        else if (buffer_str.find("\n") != string::npos || buffer_str.find("\\n") != string::npos) {
          buffer_str = "";
          found_vasp = false; // if there are newlines in the buffer then it is the wrong line
        }
      }
      infile.seekg(-10, std::ios::cur); // shift cursor to avoid the case where "vasp." is on the boundary of two buffers
    }
    return vaspVersion;
  }

  string BIN2VASPVersionNumber(const string& binfile){  //SD20220331
    return VASPVersionString2Number(BIN2VASPVersion(binfile));
  }

  double BIN2VASPVersionDouble(const string& binfile){  //SD20220331
    return VASPVersionString2Double(BIN2VASPVersion(binfile));
  }

  string OUTCAR2VASPVersion(const string& outcar){  //CO20210315
    //outcar -> vasp.4.6.35
    //outcar -> vasp.5.4.4.18Apr17-6-g9f103f2a35
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " outcar=" << outcar << endl;}
    if(!aurostd::FileExist(outcar)){return "";}
    vector<string> vlines;
    uint vlines_size=aurostd::file2vectorstring(outcar,vlines);
    vector<string> tokens;
    uint tokens_size=0;
    for(uint iline=0;iline<vlines_size;iline++){
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " vlines[iline]=\"" << vlines[iline] << "\"" << endl;}
      if(vlines[iline].find("vasp.")!=string::npos){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " FOUND 'vasp.' line" << endl;}
        tokens_size=aurostd::string2tokens(vlines[iline],tokens," ");
        for(uint i=0;i<tokens_size;i++){
          if(tokens[i].find("vasp.")!=string::npos){
            return tokens[i];
          }
        }
      }
    }
    return "";
  }

  string OUTCAR2VASPVersionNumber(const string& outcar){  //CO20210315
    //outcar -> 4.6.35
    //outcar -> 5.4.4
    return VASPVersionString2Number(OUTCAR2VASPVersion(outcar));
  }

  double OUTCAR2VASPVersionDouble(const string& outcar){  //CO20210315
    //outcar -> 4.635
    //outcar -> 5.44
    return VASPVersionString2Double(OUTCAR2VASPVersion(outcar));
  }

  string VASPVersionString2Number(const string& vasp_version){  //CO20210315
    //vasp.4.6.35 -> 4.6.35
    //vasp.5.4.4.18Apr17-6-g9f103f2a35 -> 5.4.4
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    string version_str=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vasp_version);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " version_str=\"" << version_str << "\"" << endl;}
    if(version_str.empty()){return "";}
    aurostd::StringSubst(version_str,"vasp.",""); //remove 'vasp.'
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " version_str=\"" << version_str << "\"" << endl;}
    //isfloat() does not work here: "35 3Apr08" is considered float: 35
    vector<string> vtokens;
    aurostd::string2tokens(version_str,vtokens,".");  //split by '.', check if pieces are all digits
    string version_str_num="";
    uint i=0,j=0;
    bool all_digits=true;
    for(i=0;i<vtokens.size();i++){
      all_digits=true;
      for(j=0;j<vtokens[i].size()&&all_digits;j++){
        if(!isdigit(vtokens[i][j])){all_digits=false;}
      }
      if(all_digits){
        if(version_str_num.empty()){version_str_num+=vtokens[i];}
        else{version_str_num+="."+vtokens[i];}
      }else{break;} //stop at 18Apr17...
    }
    //[CO20210315 - not good, will catch 18 (date) in vasp.5.4.4.18Apr17-6-g9f103f2a35]for(uint i=0;i<version_str.size();i++){
    //[CO20210315 - not good, will catch 18 (date) in vasp.5.4.4.18Apr17-6-g9f103f2a35]  if(isdigit(version_str[i]) || version_str[i]=='.'){
    //[CO20210315 - not good, will catch 18 (date) in vasp.5.4.4.18Apr17-6-g9f103f2a35]    version_str_num+=version_str[i];
    //[CO20210315 - not good, will catch 18 (date) in vasp.5.4.4.18Apr17-6-g9f103f2a35]  }else{break;}
    //[CO20210315 - not good, will catch 18 (date) in vasp.5.4.4.18Apr17-6-g9f103f2a35]}
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " version_str_num=\"" << version_str_num << "\"" << endl;}
    if(version_str_num.empty()){return "";}  //repetita iuvant
    return version_str_num;
  }

  double VASPVersionString2Double(const string& vasp_version){  //CO20210315
    //SD20220331 - Changed how the double is returned, so we can do version comparison by comparing doubles,
    //for example now: 4.1.311 > 4.1.4 since 4.001311 > 4.001004
    //vasp.4.6.35 -> 4.006035
    //vasp.5.4.4.18Apr17-6-g9f103f2a35 -> 5.004004
    //differs from 2Number in that it returns a double
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    string version_str=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vasp_version);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " version_str=\"" << version_str << "\"" << endl;}
    return aurostd::VersionString2Double(KBIN::VASPVersionString2Number(version_str));
  }

  /// @brief gets a string containing the VASP version
  ///
  /// @param binfile binary filename
  ///
  /// @return string containing the VASP version
  ///
  /// @note Versions 6.3+ of VASP can only be compiled in parallel. Attemping to run VASP in serial
  /// will cause VASP, as well as AFLOW, to hang (not crash), wasting precious CPU hours.
  /// Therefore, either one would need to always call VASP using MPI commands or read the version
  /// from the binary, which would need to be updated periodically. We have chosen to go with the
  /// latter, as the worst case is an unknown VASP version.
  ///
  /// @authors
  /// @mod{SD,20220923,updated for vasp6.3}
  /// @mod{CO,20210315,updated for vasp5}
  /// @mod{ME,20200114,return empty string instead of throwing xerror when the binary
  /// is not found or not a valid VASP binary}
  /// @mod{ME,20190219,created function}

  //ME20190219 - getVASPVersionString
  // Retrives the VASP version of a binary file.
  // Taken from old APL/apl_hroutines
  //ME20200114 - Return empty string instead of throwing xerror when the binary
  // is not found or not a valid VASP binary. Throwing errors would kill aflow
  // when aflow.in files are moved between machines and the VASP binary files
  // have different names. This is not desirable when VASP does not need to be
  // run (e.g. for post-processing).
  string getVASPVersion(const string& binfile) {
    // /home/bin/vasp_std -> vasp.4.6.35
    // /home/bin/vasp_std -> vasp.5.4.4.18Apr17-6-g9f103f2a35
    bool LDEBUG=(FALSE || _DEBUG_KVASP_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " binfile=" << binfile << endl;}
    if (!XHOST.is_command(binfile)) return "";
    // Get the full path to the binary
    string fullPathBinaryName = XHOST.command(binfile);
    if (fullPathBinaryName.empty()) return "";
    string vaspVersion = KBIN::BIN2VASPVersion(fullPathBinaryName);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " vaspVersion from BIN=" << vaspVersion << endl;}
    return vaspVersion;
  }

  string getVASPVersionNumber(const string& binfile) {  //CO20200610
    // /home/bin/vasp_std -> 4.6.35
    // /home/bin/vasp_std -> 5.4.4
    return VASPVersionString2Number(getVASPVersion(binfile));
  }

  double getVASPVersionDouble(const string& binfile) {  //CO20200610
    // /home/bin/vasp_std -> 4.635
    // /home/bin/vasp_std -> 5.44
    return VASPVersionString2Double(getVASPVersion(binfile));
  }
}  // namespace KBIN

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
