// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by SC 2017-2018

#include "aflow.h"

// ***************************************************************************
#ifdef _AFLOW_AFLOWRC_H_

#define AFLOWRC_FILENAME_LOCAL   XHOST.home+"/.aflow.rc"
#define AFLOWRC_FILENAME_GLOBAL  "/etc/aflow.conf"

// DEFAULT DEFINITIONS
#define AFLOWRC_AFLOWRC string(AFLOW_VERSION)

// DEFAULT DEFINITIONS
#define AFLOWRC_DEFAULT_KZIP_BIN                        string("xz")
#define         DEFAULT_KZIP_BIN                        XHOST.adefault.getattachedscheme("DEFAULT_KZIP_BIN")
#define AFLOWRC_DEFAULT_KZIP_EXT                        string(".xz") 
#define         DEFAULT_KZIP_EXT                        XHOST.adefault.getattachedscheme("DEFAULT_KZIP_EXT")
#define AFLOWRC_DEFAULT_TMPFS_DIRECTORIES               string("/tmp/,/run/shm/,/dev/shm/")
#define         DEFAULT_TMPFS_DIRECTORIES               XHOST.adefault.getattachedscheme("DEFAULT_TMPFS_DIRECTORIES")


//HE20220218 START
// DEFAULTS ENTRY LOADER
#define AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE      string("~/.aflow/aflowlib_alloy.db")
#define         DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE      XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER       string("aflowlib.duke.edu")
#define         DEFAULT_ENTRY_LOADER_AFLUX_SERVER       XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_AFLUX_SERVER")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH         string("/API/aflux/")
#define         DEFAULT_ENTRY_LOADER_AFLUX_PATH         XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_AFLUX_PATH")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER     string("aflowlib.duke.edu")
#define         DEFAULT_ENTRY_LOADER_RESTAPI_SERVER     XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_RESTAPI_SERVER")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH       string("/AFLOWDATA/")
#define         DEFAULT_ENTRY_LOADER_RESTAPI_PATH       XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_RESTAPI_PATH")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH            string("/common/")
#define         DEFAULT_ENTRY_LOADER_FS_PATH            XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_FS_PATH")
//HE20220218 STOP

//ME20191001 START
// DEFAULTS AFLOW DATABASE
#define AFLOWRC_DEFAULT_AFLOW_DB_FILE                   string("/var/cache/aflow_data/AFLOWDB/aflowlib.db")
#define         DEFAULT_AFLOW_DB_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_FILE")
#define AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE             string("/var/cache/aflow_data/AFLOWDB/aflowlib.json")
#define         DEFAULT_AFLOW_DB_STATS_FILE             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_STATS_FILE")
#define AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH              string("/common/AFLOW/LIBS/")
#define         DEFAULT_AFLOW_DB_DATA_PATH              XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_DATA_PATH")
#define AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE              string("/var/cache/aflow_data/AFLOWDB/ADB_Idle.lock")
#define         DEFAULT_AFLOW_DB_LOCK_FILE              XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_LOCK_FILE")
#define AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD        3*3600
#define         DEFAULT_AFLOW_DB_STALE_THRESHOLD        XHOST.adefault.getattachedutype<long int>("DEFAULT_AFLOW_DB_STALE_THRESHOLD")
//ME20191001 STOP

// FILENAMES FOR AFLOW.ORG ANALYSIS
#define AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_OUT         string("aflowlib.out")
#define         DEFAULT_FILE_AFLOWLIB_ENTRY_OUT         XHOST.adefault.getattachedscheme("DEFAULT_FILE_AFLOWLIB_ENTRY_OUT")
#define AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_JSON        string("aflowlib.json")
#define         DEFAULT_FILE_AFLOWLIB_ENTRY_JSON        XHOST.adefault.getattachedscheme("DEFAULT_FILE_AFLOWLIB_ENTRY_JSON")
#define AFLOWRC_DEFAULT_FILE_EDATA_ORIG_OUT             string("edata.orig.out")
#define         DEFAULT_FILE_EDATA_ORIG_OUT             XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_ORIG_OUT")
#define AFLOWRC_DEFAULT_FILE_EDATA_RELAX_OUT            string("edata.relax.out")
#define         DEFAULT_FILE_EDATA_RELAX_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_RELAX_OUT")
#define AFLOWRC_DEFAULT_FILE_EDATA_BANDS_OUT            string("edata.bands.out")
#define         DEFAULT_FILE_EDATA_BANDS_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_BANDS_OUT")
#define AFLOWRC_DEFAULT_FILE_DATA_ORIG_OUT              string("data.orig.out")
#define         DEFAULT_FILE_DATA_ORIG_OUT              XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_ORIG_OUT")
#define AFLOWRC_DEFAULT_FILE_DATA_RELAX_OUT             string("data.relax.out")
#define         DEFAULT_FILE_DATA_RELAX_OUT             XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_RELAX_OUT")
#define AFLOWRC_DEFAULT_FILE_DATA_BANDS_OUT             string("data.bands.out")
#define         DEFAULT_FILE_DATA_BANDS_OUT             XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_BANDS_OUT")
#define AFLOWRC_DEFAULT_FILE_EDATA_ORIG_JSON            string("edata.orig.json")
#define         DEFAULT_FILE_EDATA_ORIG_JSON            XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_ORIG_JSON")
#define AFLOWRC_DEFAULT_FILE_EDATA_RELAX_JSON           string("edata.relax.json")
#define         DEFAULT_FILE_EDATA_RELAX_JSON           XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_RELAX_JSON")
#define AFLOWRC_DEFAULT_FILE_EDATA_BANDS_JSON           string("edata.bands.json")
#define         DEFAULT_FILE_EDATA_BANDS_JSON           XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_BANDS_JSON")
#define AFLOWRC_DEFAULT_FILE_DATA_ORIG_JSON             string("data.orig.json")
#define         DEFAULT_FILE_DATA_ORIG_JSON             XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_ORIG_JSON")
#define AFLOWRC_DEFAULT_FILE_DATA_RELAX_JSON            string("data.relax.json")
#define         DEFAULT_FILE_DATA_RELAX_JSON            XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_RELAX_JSON")
#define AFLOWRC_DEFAULT_FILE_DATA_BANDS_JSON            string("data.bands.json")
#define         DEFAULT_FILE_DATA_BANDS_JSON            XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_BANDS_JSON")
#define AFLOWRC_DEFAULT_FILE_TIME_OUT                   string("time")
#define         DEFAULT_FILE_TIME_OUT                   XHOST.adefault.getattachedscheme("DEFAULT_FILE_TIME_OUT")
#define AFLOWRC_DEFAULT_FILE_SPACEGROUP1_OUT            string("SpaceGroup")
#define         DEFAULT_FILE_SPACEGROUP1_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_SPACEGROUP1_OUT")
#define AFLOWRC_DEFAULT_FILE_SPACEGROUP2_OUT            string("SpaceGroup2")
#define         DEFAULT_FILE_SPACEGROUP2_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_SPACEGROUP2_OUT")
#define AFLOWRC_DEFAULT_FILE_VOLDISTPARAMS_OUT          string("VOLDISTParams")
#define         DEFAULT_FILE_VOLDISTPARAMS_OUT          XHOST.adefault.getattachedscheme("DEFAULT_FILE_VOLDISTPARAMS_OUT")
#define AFLOWRC_DEFAULT_FILE_VOLDISTEVOLUTION_OUT       string("VOLDISTEvolution")
#define         DEFAULT_FILE_VOLDISTEVOLUTION_OUT       XHOST.adefault.getattachedscheme("DEFAULT_FILE_VOLDISTEVOLUTION_OUT")

// FILENAMES FOR AFLOW OPERATION
#define AFLOWRC_DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT  string("aflow.pseudopotential_auid.out") 
#define         DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT  XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_OUT             string("aflow.prescript.out") 
#define         DEFAULT_AFLOW_PRESCRIPT_OUT             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PRESCRIPT_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_COMMAND         string("aflow.prescript.command") 
#define         DEFAULT_AFLOW_PRESCRIPT_COMMAND         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PRESCRIPT_COMMAND")
#define AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_OUT            string("aflow.postscript.out") 
#define         DEFAULT_AFLOW_POSTSCRIPT_OUT            XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_POSTSCRIPT_OUT")
#define AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_COMMAND        string("aflow.postscript.command") 
#define         DEFAULT_AFLOW_POSTSCRIPT_COMMAND        XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_POSTSCRIPT_COMMAND")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_OUT                string("aflow.pgroup.out")
#define         DEFAULT_AFLOW_PGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_JSON               string("aflow.pgroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_PGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_OUT           string("aflow.pgroup_xtal.out")
#define         DEFAULT_AFLOW_PGROUP_XTAL_OUT           XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_XTAL_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_JSON          string("aflow.pgroup_xtal.json") //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_PGROUP_XTAL_JSON          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_XTAL_JSON")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT     string("aflow.pgroupk_Patterson.out") //DX20200129
#define         DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT     XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT") //DX20200129
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON    string("aflow.pgroupk_Patterson.json") //DX20200129
#define         DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON    XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON") //DX20200129
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_OUT               string("aflow.pgroupk.out")
#define         DEFAULT_AFLOW_PGROUPK_OUT               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_JSON              string("aflow.pgroupk.json")     //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_PGROUPK_JSON              XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_JSON")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_OUT          string("aflow.pgroupk_xtal.out") //DX20171205 - Added pgroupk_xtal
#define         DEFAULT_AFLOW_PGROUPK_XTAL_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_XTAL_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_JSON         string("aflow.pgroupk_xtal.json")//DX20170802 - Add JSON //DX20171205 - Added pgroupk_xtal
#define         DEFAULT_AFLOW_PGROUPK_XTAL_JSON         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_XTAL_JSON")
#define AFLOWRC_DEFAULT_AFLOW_FGROUP_OUT                string("aflow.fgroup.out")
#define         DEFAULT_AFLOW_FGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FGROUP_JSON               string("aflow.fgroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_FGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_SGROUP_OUT                string("aflow.sgroup.out")
#define         DEFAULT_AFLOW_SGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_SGROUP_JSON               string("aflow.sgroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_SGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_AGROUP_OUT                string("aflow.agroup.out")
#define         DEFAULT_AFLOW_AGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_AGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_AGROUP_JSON               string("aflow.agroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_AGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_AGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_IATOMS_OUT                string("aflow.iatoms.out")
#define         DEFAULT_AFLOW_IATOMS_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IATOMS_OUT")
#define AFLOWRC_DEFAULT_AFLOW_IATOMS_JSON               string("aflow.iatoms.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_IATOMS_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IATOMS_JSON")
#define AFLOWRC_DEFAULT_AFLOW_ICAGES_OUT                string("aflow.icages.out")
#define         DEFAULT_AFLOW_ICAGES_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_ICAGES_OUT")
#define AFLOWRC_DEFAULT_AFLOW_SURFACE_OUT               string("aflow.surface.out")
#define         DEFAULT_AFLOW_SURFACE_OUT               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SURFACE_OUT")
#define AFLOWRC_DEFAULT_AFLOW_QMVASP_OUT                string("aflow.qmvasp.out")
#define         DEFAULT_AFLOW_QMVASP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_QMVASP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_ERVASP_OUT                string("aflow.error.out")
#define         DEFAULT_AFLOW_ERVASP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_ERVASP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_IMMISCIBILITY_OUT         string("aflow.immiscibility.out")
#define         DEFAULT_AFLOW_IMMISCIBILITY_OUT         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IMMISCIBILITY_OUT")
#define AFLOWRC_DEFAULT_AFLOW_MEMORY_OUT                string("aflow.memory.out")
#define         DEFAULT_AFLOW_MEMORY_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_MEMORY_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_INPUT_OUT          string("aflow.frozsl_input.out")
#define         DEFAULT_AFLOW_FROZSL_INPUT_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_INPUT_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_POSCAR_OUT         string("aflow.frozsl_poscar.out")
#define         DEFAULT_AFLOW_FROZSL_POSCAR_OUT         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_POSCAR_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_MODES_OUT          string("aflow.frozsl_energies.out")
#define         DEFAULT_AFLOW_FROZSL_MODES_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_MODES_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_EIGEN_OUT          string("aflow.frozsl_eigen.out")
#define         DEFAULT_AFLOW_FROZSL_EIGEN_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_EIGEN_OUT")
#define AFLOWRC_DEFAULT_AFLOW_END_OUT                   string("aflow.end.out")
#define         DEFAULT_AFLOW_END_OUT                   XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_END_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PLASMONICS_FILE           string("aflow.plasmonics_eps")
#define         DEFAULT_AFLOW_PLASMONICS_FILE           XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PLASMONICS_FILE")

// GENERIC MPI   // DONE
#define AFLOWRC_MPI_START_DEFAULT                       string("")
#define         MPI_START_DEFAULT                       XHOST.adefault.getattachedscheme("MPI_START_DEFAULT")
#define AFLOWRC_MPI_STOP_DEFAULT                        string("")
#define         MPI_STOP_DEFAULT                        XHOST.adefault.getattachedscheme("MPI_STOP_DEFAULT")
#define AFLOWRC_MPI_COMMAND_DEFAULT                     string("mpirun -np")
#define         MPI_COMMAND_DEFAULT                     XHOST.adefault.getattachedscheme("MPI_COMMAND_DEFAULT")
#define AFLOWRC_MPI_NCPUS_DEFAULT                       4
#define         MPI_NCPUS_DEFAULT                       XHOST.adefault.getattachedutype<int>("MPI_NCPUS_DEFAULT")
#define AFLOWRC_MPI_NCPUS_MAX                           init::GetCPUCores() //16  //CO20180124
#define         MPI_NCPUS_MAX                           XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MAX")

// BINARY    // DONE
#define AFLOWRC_DEFAULT_VASP_GAMMA_BIN                  string("vasp46s_gamma")
#define         DEFAULT_VASP_GAMMA_BIN                  XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_BIN")
#define AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN              string("mpivasp46s_gamma")
#define         DEFAULT_VASP_GAMMA_MPI_BIN              XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_MPI_BIN")
#define AFLOWRC_DEFAULT_VASP_BIN                        string("vasp46s")
#define         DEFAULT_VASP_BIN                        XHOST.adefault.getattachedscheme("DEFAULT_VASP_BIN")
#define AFLOWRC_DEFAULT_VASP_MPI_BIN                    string("mpivasp46s")
#define         DEFAULT_VASP_MPI_BIN                    XHOST.adefault.getattachedscheme("DEFAULT_VASP_MPI_BIN")
#define AFLOWRC_DEFAULT_VASP5_BIN                       string("vasp54s")
#define         DEFAULT_VASP5_BIN                       XHOST.adefault.getattachedscheme("DEFAULT_VASP5_BIN")
#define AFLOWRC_DEFAULT_VASP5_MPI_BIN                   string("mpivasp54s")
#define         DEFAULT_VASP5_MPI_BIN                   XHOST.adefault.getattachedscheme("DEFAULT_VASP5_MPI_BIN")
//aims
#define AFLOWRC_DEFAULT_AIMS_BIN                        string("aims")
#define         DEFAULT_AIMS_BIN                        XHOST.adefault.getattachedscheme("DEFAULT_AIMS_BIN")

// POTCARS // DONE
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES               string("/common/VASP,/common/AFLOW/VASP,/home/aflow/common/AFLOW/VASP,/fslhome/fslcollab8/group/VASP,/fslhome/glh43/src/,/share/home/00470/tg457283/common/AFLOW/VASP/,/share/home/00457/tg457357/common/AFLOW/VASP/,/home/mehl/bin/AFLOW/VASP/,~/common/VASP/,~/common/AFLOW/VASP/,/home/aflow/common/VASP/,/nics/a/proj/aflow/common/AFLOW/VASP/,/home/users/aflow/common/VASP,/share/apps/AFLOW3/VASP,/share/apps/vasp/PP,/projects/kyang-group/common/VASP,/home/Tools/src/vasp/,/somewhere/")  // first is default, tokenized with "," //DX20190107 - added CMU path
#define         DEFAULT_VASP_POTCAR_DIRECTORIES               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIRECTORIES")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE                      string("current")
#define         DEFAULT_VASP_POTCAR_DATE                      XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX                    string("/POTCAR")
#define         DEFAULT_VASP_POTCAR_SUFFIX                    XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_SUFFIX")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA              string("01Apr2000")   // when no date is given for pot_LDA
#define         DEFAULT_VASP_POTCAR_DATE_POT_LDA              XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA              string("01Apr2000")   // when no date is given for pot_GGA
#define         DEFAULT_VASP_POTCAR_DATE_POT_GGA              XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_GGA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_LDA               string("pot_LDA")
#define         DEFAULT_VASP_POTCAR_DIR_POT_LDA               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_GGA               string("pot_GGA")
#define         DEFAULT_VASP_POTCAR_DIR_POT_GGA               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_GGA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_PBE               string("pot_PBE")
#define         DEFAULT_VASP_POTCAR_DIR_POT_PBE               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_PBE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA            string("potpaw_LDA")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA            string("potpaw_GGA")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE            string("potpaw_PBE")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN        string("potpaw_LDA.54")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN        XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN        string("potpaw_PBE.54")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN        XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN")

// KPOINTS/DOS // DONE
#define AFLOWRC_DEFAULT_BANDS_GRID                            20
#define         DEFAULT_BANDS_GRID                            XHOST.adefault.getattachedutype<int>("DEFAULT_BANDS_GRID") 
#define AFLOWRC_DEFAULT_BANDS_LATTICE                         string("AUTO")
#define         DEFAULT_BANDS_LATTICE                         XHOST.adefault.getattachedscheme("DEFAULT_BANDS_LATTICE")
#define AFLOWRC_DEFAULT_KSCHEME                               string("AUTO")  //WSETYAWAN DEFAULT  //CO20200329  // WAS M but ivasp fixes
#define         DEFAULT_KSCHEME                               XHOST.adefault.getattachedscheme("DEFAULT_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA                                 6000
#define         DEFAULT_KPPRA                                 XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA")
#define AFLOWRC_DEFAULT_KPPRA_STATIC                          10000
#define         DEFAULT_KPPRA_STATIC                          XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_STATIC")
#define AFLOWRC_DEFAULT_STATIC_KSCHEME                        string("AUTO")  //WSETYAWAN DEFAULT  //CO20200329  // WAS M but ivasp fixes
#define         DEFAULT_STATIC_KSCHEME                        XHOST.adefault.getattachedscheme("DEFAULT_STATIC_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA_ICSD                            8000
#define         DEFAULT_KPPRA_ICSD                            XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_ICSD")
#define AFLOWRC_DEFAULT_UNARY_BANDS_GRID                      128
#define         DEFAULT_UNARY_BANDS_GRID                      XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_BANDS_GRID")
#define AFLOWRC_DEFAULT_UNARY_KPPRA                           8000 // 32768 // 27000
#define         DEFAULT_UNARY_KPPRA                           XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA")
#define AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC                    8000 // 32768 // 27000
#define         DEFAULT_UNARY_KPPRA_STATIC                    XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA_STATIC")
#define AFLOWRC_DEFAULT_PHONONS_KSCHEME                       string("G")
#define         DEFAULT_PHONONS_KSCHEME                       XHOST.adefault.getattachedscheme("DEFAULT_PHONONS_KSCHEME")
#define AFLOWRC_DEFAULT_PHONONS_KPPRA                         3000  //CO20181226  //ME20190205 - 8000 uses too much memory e.g. for NaF - 2000 appears sufficient; ME20200108 - 3000 minimum for metals
#define         DEFAULT_PHONONS_KPPRA                         XHOST.adefault.getattachedutype<int>("DEFAULT_PHONONS_KPPRA") //CO20181226 //ME20190112
#define AFLOWRC_DEFAULT_DOS_EMIN                              -10.0
#define         DEFAULT_DOS_EMIN                              XHOST.adefault.getattachedutype<double>("DEFAULT_DOS_EMIN")
#define AFLOWRC_DEFAULT_DOS_EMAX                              10.0
#define         DEFAULT_DOS_EMAX                              XHOST.adefault.getattachedutype<double>("DEFAULT_DOS_EMAX")
#define AFLOWRC_DEFAULT_DOS_SCALE                             1.2
#define         DEFAULT_DOS_SCALE                             XHOST.adefault.getattachedutype<double>("DEFAULT_DOS_SCALE")

// PRECISION // DONE
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW                   1.0
#define         DEFAULT_VASP_PREC_ENMAX_LOW                   XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_LOW")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM                1.3
#define         DEFAULT_VASP_PREC_ENMAX_MEDIUM                XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_MEDIUM")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL                1.3
#define         DEFAULT_VASP_PREC_ENMAX_NORMAL                XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_NORMAL")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH                  1.4
#define         DEFAULT_VASP_PREC_ENMAX_HIGH                  XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_HIGH")
#define AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE              1.4
#define         DEFAULT_VASP_PREC_ENMAX_ACCURATE              XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_ENMAX_ACCURATE")
#define AFLOWRC_DEFAULT_VASP_ENMAX_MINIMUM                    0.25
#define         DEFAULT_VASP_ENMAX_MINIMUM                    XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_ENMAX_MINIMUM")
#define AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF               0.05
#define         DEFAULT_VASP_SPIN_REMOVE_CUTOFF               XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_SPIN_REMOVE_CUTOFF")
#define AFLOWRC_DEFAULT_VASP_PREC_POTIM                       0.5
#define         DEFAULT_VASP_PREC_POTIM                       XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_POTIM")
#define AFLOWRC_DEFAULT_VASP_PREC_EDIFFG                      -1E-3
#define         DEFAULT_VASP_PREC_EDIFFG                      XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_PREC_EDIFFG")

// OPTIONS // DONE
#define AFLOWRC_DEFAULT_VASP_OUT                                      string("vasp.out")
#define         DEFAULT_VASP_OUT                                      XHOST.adefault.getattachedscheme("DEFAULT_VASP_OUT")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR                           string("./INCAR")
#define         DEFAULT_VASP_EXTERNAL_INCAR                           XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_INCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR                          string("./POSCAR")
#define         DEFAULT_VASP_EXTERNAL_POSCAR                          XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POSCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR                          string("./POTCAR")
#define         DEFAULT_VASP_EXTERNAL_POTCAR                          XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POTCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS                         string("./KPOINTS")
#define         DEFAULT_VASP_EXTERNAL_KPOINTS                         XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_KPOINTS")
#define AFLOWRC_DEFAULT_AIMS_EXTERNAL_CONTROL                         string("./control.in")
#define         DEFAULT_AIMS_EXTERNAL_CONTROL                         XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_CONTROL")
#define AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM                            string("./geometry.in")
#define         DEFAULT_AIMS_EXTERNAL_GEOM                            XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_GEOM")
#define AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE                     string("potpaw_PBE")
#define         DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE                     XHOST.adefault.getattachedscheme("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME           string("ENERGY")
#define         DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME           XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT                 2 //CO20181226
#define         DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT                 XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT") //CO20181226
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME                 string("ACCURATE")
#define         DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME                 string("NORMAL")
#define         DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME              string("NONE")
#define         DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME              XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME                 string("0")  
#define         DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME                 string("DEFAULT")
#define         DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME               1 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME               XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME") //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME        -5 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME        XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME") //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME         0 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME         XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME") //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA                       0.1 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_SIGMA                       XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_FORCE_OPTION_SIGMA")  //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC                0.05 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC                XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC")  //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS                 0.05 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS                 XHOST.adefault.getattachedutype<double>("DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS")  //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM                        60  //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_NELM                        XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_NELM") //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM_STATIC                 120 //CO20200624
#define         DEFAULT_VASP_FORCE_OPTION_NELM_STATIC                 XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_NELM_STATIC")  //CO20200624
#define AFLOWRC_MAX_VASP_NELM                                         300  //CO20200624
#define         MAX_VASP_NELM                                         XHOST.adefault.getattachedutype<int>("MAX_VASP_NELM") //CO20200624
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME                string("DEFAULT")
#define         DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME                XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM                         TRUE
#define         DEFAULT_VASP_FORCE_OPTION_SYM                         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SYM") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN                        TRUE
#define         DEFAULT_VASP_FORCE_OPTION_SPIN                        XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1         FALSE
#define         DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2         TRUE  //ME20190308 - remove spin after two relaxations if zero
#define         DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER                       FALSE
#define         DEFAULT_VASP_FORCE_OPTION_BADER                       XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_BADER") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER_STATIC                TRUE
#define         DEFAULT_VASP_FORCE_OPTION_BADER_STATIC                XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_BADER_STATIC") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF                         FALSE
#define         DEFAULT_VASP_FORCE_OPTION_ELF                         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_ELF") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM                 FALSE   // TRUE
#define         DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM                 XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR                     FALSE
#define         DEFAULT_VASP_FORCE_OPTION_WAVECAR                     XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_WAVECAR") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR                      TRUE
#define         DEFAULT_VASP_FORCE_OPTION_CHGCAR                      XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_CHGCAR") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING                  FALSE
#define         DEFAULT_VASP_FORCE_OPTION_LSCOUPLING                  XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_LSCOUPLING") 

// AFLOW_LIBRARY AFLOW_PROJECT // DONE
#define AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES             string("/common/AFLOW/LIBS/,/home/aflow/common/AFLOW/LIBS/,/fslhome/glh43/src/,/usr/local/bin/,/fslhome/fslcollab8/group/bin/,/home/auro/work/AFLOW3/,~/common/AFLOW/LIBS/,./,/nics/a/proj/aflow/common/AFLOW/LIBS/,/home/users/aflow/common/AFLOW/LIBS,/home/junkai/PROTO_DATABASE/,/projects/kyang-group/common/LIBS,/somewhere/")  // first is default, tokenized with ","
#define         DEFAULT_AFLOW_LIBRARY_DIRECTORIES             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_LIBRARY_DIRECTORIES")
#define AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES            string("/common/AUID,/common/ICSD,/common/LIB0,/common/LIB1,/common/LIB2,/common/LIB3,/common/LIB4,/common/LIB5,/common/LIB6,/common/LIB7,/common/LIB8,/common/LIB9")  // first is default, tokenized with ","
#define         DEFAULT_AFLOW_PROJECTS_DIRECTORIES            XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PROJECTS_DIRECTORIES")
#define AFLOWRC_DEFAULT_AFLOWDATA_WEB_DIRECTORY               string("/www/AFLOWDATA")  //CO+ME20200731
#define         DEFAULT_AFLOWDATA_WEB_DIRECTORY               XHOST.adefault.getattachedscheme("DEFAULT_AFLOWDATA_WEB_DIRECTORY") //CO+ME20200731

// PLATON/FINDSYM // DONE
#define AFLOWRC_DEFAULT_PLATON_P_EQUAL                        FALSE
#define         DEFAULT_PLATON_P_EQUAL                        XHOST.adefault.getattachedutype<bool>("DEFAULT_PLATON_P_EQUAL") 
#define AFLOWRC_DEFAULT_PLATON_P_EXACT                        FALSE
#define         DEFAULT_PLATON_P_EXACT                        XHOST.adefault.getattachedutype<bool>("DEFAULT_PLATON_P_EXACT") 
#define AFLOWRC_DEFAULT_PLATON_P_ANG                          1.0
#define         DEFAULT_PLATON_P_ANG                          XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_ANG") 
#define AFLOWRC_DEFAULT_PLATON_P_D1                           0.25
#define         DEFAULT_PLATON_P_D1                           XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_D1") 
#define AFLOWRC_DEFAULT_PLATON_P_D2                           0.25
#define         DEFAULT_PLATON_P_D2                           XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_D2") 
#define AFLOWRC_DEFAULT_PLATON_P_D3                           0.25
#define         DEFAULT_PLATON_P_D3                           XHOST.adefault.getattachedutype<double>("DEFAULT_PLATON_P_D3") 
#define AFLOWRC_DEFAULT_FINDSYM_TOL                           1.0e-3
#define         DEFAULT_FINDSYM_TOL                           XHOST.adefault.getattachedutype<double>("DEFAULT_FINDSYM_TOL") 

// GNUPLOT // DONE
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT                      string("Helvetica")              // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT                      XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD                 string("Helvetica-Bold")         // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_BOLD                 XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS              string("Helvetica-Oblique")      // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_ITALICS              XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS         string("Helvetica-BoldOblique")  // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS         XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT                      string("Arial")                  // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT                      XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD                 string("Arial_Bold")             // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_BOLD                 XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS              string("Arial_Italic")           // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_ITALICS              XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS         string("Arial_BoldItalic")       // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS         XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT                    string("Symbol") 
#define         DEFAULT_GNUPLOT_GREEK_FONT                    XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD               string("Symbol-Bold")
#define         DEFAULT_GNUPLOT_GREEK_FONT_BOLD               XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_ITALICS            string("Symbol-Oblique")
#define         DEFAULT_GNUPLOT_GREEK_FONT_ITALICS            XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS       string("Symbol-BoldOblique")
#define         DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS       XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS")

// DEFAULT CHULL
#define AFLOWRC_DEFAULT_CHULL_ALLOWED_DFT_TYPES                           string("PAW_PBE,PAW_PBE_KIN")
#define         DEFAULT_CHULL_ALLOWED_DFT_TYPES                           XHOST.adefault.getattachedscheme("DEFAULT_CHULL_ALLOWED_DFT_TYPES")
#define AFLOWRC_DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES                FALSE
#define         DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES                XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES")
#define AFLOWRC_DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES                    200
#define         DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES                    XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES")
#define AFLOWRC_DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS                    TRUE
#define         DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS                    XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS")
#define AFLOWRC_DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES   75
#define         DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES   XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES")
#define AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER                          3.25
#define         DEFAULT_CHULL_OUTLIER_MULTIPLIER                          XHOST.adefault.getattachedutype<double>("DEFAULT_CHULL_OUTLIER_MULTIPLIER")
#define AFLOWRC_DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED                  TRUE
#define         DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED                  XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED")
#define AFLOWRC_DEFAULT_CHULL_LATEX_BANNER                                2 //CO20180827
#define         DEFAULT_CHULL_LATEX_BANNER                                XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_BANNER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN                      FALSE
#define         DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN")
#define AFLOWRC_DEFAULT_CHULL_LATEX_STOICH_HEADER                         FALSE
#define         DEFAULT_CHULL_LATEX_STOICH_HEADER                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_STOICH_HEADER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNARIES                          FALSE
#define         DEFAULT_CHULL_LATEX_PLOT_UNARIES                          XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_PLOT_UNARIES")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL                         -1
#define         DEFAULT_CHULL_LATEX_PLOT_OFF_HULL                         XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_PLOT_OFF_HULL")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE                         FALSE
#define         DEFAULT_CHULL_LATEX_PLOT_UNSTABLE                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_PLOT_UNSTABLE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME                         string("")
#define         DEFAULT_CHULL_LATEX_FILTER_SCHEME                         XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_FILTER_SCHEME")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE                          0.0
#define         DEFAULT_CHULL_LATEX_FILTER_VALUE                          XHOST.adefault.getattachedutype<double>("DEFAULT_CHULL_LATEX_FILTER_VALUE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR                             TRUE
#define         DEFAULT_CHULL_LATEX_COLOR_BAR                             XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COLOR_BAR")
#define AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP                              TRUE
#define         DEFAULT_CHULL_LATEX_HEAT_MAP                              XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_HEAT_MAP")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT                        TRUE
#define         DEFAULT_CHULL_LATEX_COLOR_GRADIENT                        XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COLOR_GRADIENT")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP                             string("")
#define         DEFAULT_CHULL_LATEX_COLOR_MAP                             XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_COLOR_MAP")
#define AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR                   string("")
#define         DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR                   XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR")
#define AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS                          FALSE
#define         DEFAULT_CHULL_LATEX_REVERSE_AXIS                          XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_REVERSE_AXIS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW                FALSE
#define         DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW                XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LINKS                                 1
#define         DEFAULT_CHULL_LATEX_LINKS                                 XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_LINKS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME                            string("")
#define         DEFAULT_CHULL_LATEX_LABEL_NAME                            XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_LABEL_NAME")
#define AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS                           FALSE
#define         DEFAULT_CHULL_LATEX_META_LABELS                           XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_META_LABELS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL                       FALSE
#define         DEFAULT_CHULL_LATEX_LABELS_OFF_HULL                       XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_LABELS_OFF_HULL")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION              -1
#define         DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION              XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION")
#define AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT                        FALSE
#define         DEFAULT_CHULL_LATEX_HELVETICA_FONT                        XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_HELVETICA_FONT")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE                             string("")
#define         DEFAULT_CHULL_LATEX_FONT_SIZE                             XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_FONT_SIZE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS                         TRUE
#define         DEFAULT_CHULL_LATEX_ROTATE_LABELS                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_ROTATE_LABELS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS                           -1
#define         DEFAULT_CHULL_LATEX_BOLD_LABELS                           XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_BOLD_LABELS")
#define AFLOWRC_DEFAULT_CHULL_PNG_RESOLUTION                              300
#define         DEFAULT_CHULL_PNG_RESOLUTION                              XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_PNG_RESOLUTION")

// DEFAULT GFA
#define AFLOWRC_DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF                     0.05  //CO20190628
#define         DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF                     XHOST.adefault.getattachedutype<double>("DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF")  //CO20190628

// DEFAULT ARUN
#define AFLOWRC_ARUN_DIRECTORY_PREFIX                         string("ARUN.")
#define         ARUN_DIRECTORY_PREFIX                         XHOST.adefault.getattachedscheme("ARUN_DIRECTORY_PREFIX")

//DEFAULT POCC //CO20181226
#define AFLOWRC_DEFAULT_POCC_STRUCTURE_GENERATION_ALGO            string("UFF")
#define         DEFAULT_POCC_STRUCTURE_GENERATION_ALGO            XHOST.adefault.getattachedscheme("DEFAULT_POCC_STRUCTURE_GENERATION_ALGO")
#define AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING                   string("0:2400:300")
#define         DEFAULT_POCC_TEMPERATURE_STRING                   XHOST.adefault.getattachedscheme("DEFAULT_POCC_TEMPERATURE_STRING")
#define AFLOWRC_DEFAULT_POCC_EXCLUDE_UNSTABLE                     true  //ME20210927
#define         DEFAULT_POCC_EXCLUDE_UNSTABLE                     XHOST.adefault.getattachedutype<bool>("DEFAULT_POCC_EXCLUDE_UNSTABLE")  //ME20210927
#define AFLOWRC_DEFAULT_POCC_SITE_TOL                             0.001
#define         DEFAULT_POCC_SITE_TOL                             XHOST.adefault.getattachedutype<double>("DEFAULT_POCC_SITE_TOL")
#define AFLOWRC_DEFAULT_POCC_STOICH_TOL                           0.001
#define         DEFAULT_POCC_STOICH_TOL                           XHOST.adefault.getattachedutype<double>("DEFAULT_POCC_STOICH_TOL")
#define AFLOWRC_DEFAULT_UFF_BONDING_DISTANCE                      0.5
#define         DEFAULT_UFF_BONDING_DISTANCE                      XHOST.adefault.getattachedutype<double>("DEFAULT_UFF_BONDING_DISTANCE")
#define AFLOWRC_DEFAULT_UFF_ENERGY_TOLERANCE                      1e-6
#define         DEFAULT_UFF_ENERGY_TOLERANCE                      XHOST.adefault.getattachedutype<double>("DEFAULT_UFF_ENERGY_TOLERANCE")
#define AFLOWRC_DEFAULT_UFF_CLUSTER_RADIUS                        10
#define         DEFAULT_UFF_CLUSTER_RADIUS                        XHOST.adefault.getattachedutype<double>("DEFAULT_UFF_CLUSTER_RADIUS")
#define AFLOWRC_DEFAULT_POCC_RDF_RMAX                             50
#define         DEFAULT_POCC_RDF_RMAX                             XHOST.adefault.getattachedutype<double>("DEFAULT_POCC_RDF_RMAX")
#define AFLOWRC_DEFAULT_POCC_RDF_NBINS                            50
#define         DEFAULT_POCC_RDF_NBINS                            XHOST.adefault.getattachedutype<double>("DEFAULT_POCC_RDF_NBINS")
#define AFLOWRC_DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON  FALSE
#define         DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON  XHOST.adefault.getattachedutype<bool>("DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON")
#define AFLOWRC_DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS             TRUE
#define         DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS             XHOST.adefault.getattachedutype<bool>("DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS")
#define AFLOWRC_POCC_FILE_PREFIX                                  string("aflow.pocc.")
#define         POCC_FILE_PREFIX                                  XHOST.adefault.getattachedscheme("POCC_FILE_PREFIX")
#define AFLOWRC_POCC_OUT_FILE                                     string("out")
#define         POCC_OUT_FILE                                     XHOST.adefault.getattachedscheme("POCC_OUT_FILE")
#define AFLOWRC_POCC_APL_OUT_FILE                                 string("apl.out")  //ME20210927
#define         POCC_APL_OUT_FILE                                 XHOST.adefault.getattachedscheme("POCC_APL_OUT_FILE")  //ME20210927
#define AFLOWRC_POCC_ALL_SUPERCELLS_FILE                          string("structures_all.out")
#define         POCC_ALL_SUPERCELLS_FILE                          XHOST.adefault.getattachedscheme("POCC_ALL_SUPERCELLS_FILE")
#define AFLOWRC_POCC_UNIQUE_SUPERCELLS_FILE                       string("structures_unique.out")
#define         POCC_UNIQUE_SUPERCELLS_FILE                       XHOST.adefault.getattachedscheme("POCC_UNIQUE_SUPERCELLS_FILE")
#define AFLOWRC_POCC_ALL_HNF_MATRICES_FILE                        string("hnf_matrices.out")
#define         POCC_ALL_HNF_MATRICES_FILE                        XHOST.adefault.getattachedscheme("POCC_ALL_HNF_MATRICES_FILE")
#define AFLOWRC_POCC_ALL_SITE_CONFIGURATIONS_FILE                 string("site_configurations.out")
#define         POCC_ALL_SITE_CONFIGURATIONS_FILE                 XHOST.adefault.getattachedscheme("POCC_ALL_SITE_CONFIGURATIONS_FILE")
#define AFLOWRC_POCC_DOSCAR_FILE                                  string("DOSCAR.pocc")
#define         POCC_DOSCAR_FILE                                  XHOST.adefault.getattachedscheme("POCC_DOSCAR_FILE")
#define AFLOWRC_POCC_PHDOSCAR_FILE                                string("PHDOSCAR.pocc")  //ME20210927
#define         POCC_PHDOSCAR_FILE                                XHOST.adefault.getattachedscheme("POCC_PHDOSCAR_FILE")  //ME20210927
#define AFLOWRC_POCC_ANIONS_LIST                                  string("B,C,N,O")
#define         POCC_ANIONS_LIST                                  XHOST.adefault.getattachedscheme("POCC_ANIONS_LIST")

// DEFAULT APL
//// DEFAULT APL SUPERCELL
#define AFLOWRC_DEFAULT_APL_PREC                              string("PHONONS")
#define         DEFAULT_APL_PREC                              XHOST.adefault.getattachedscheme("DEFAULT_APL_PREC")
#define AFLOWRC_DEFAULT_APL_ENGINE                            string("DM")
#define         DEFAULT_APL_ENGINE                            XHOST.adefault.getattachedscheme("DEFAULT_APL_ENGINE")
#define AFLOWRC_DEFAULT_APL_HIBERNATE                         TRUE
#define         DEFAULT_APL_HIBERNATE                         XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_HIBERNATE")
#define AFLOWRC_DEFAULT_APL_MINSHELL                          6
#define         DEFAULT_APL_MINSHELL                          XHOST.adefault.getattachedutype<int>("DEFAULT_APL_MINSHELL")
#define AFLOWRC_DEFAULT_APL_MINATOMS                          175  //ME20190301
#define         DEFAULT_APL_MINATOMS                          XHOST.adefault.getattachedutype<int>("DEFAULT_APL_MINATOMS")
#define AFLOWRC_DEFAULT_APL_POLAR                             TRUE  //CO20181226
#define         DEFAULT_APL_POLAR                             XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_POLAR")
#define AFLOWRC_DEFAULT_APL_DMAG                              0.015
#define         DEFAULT_APL_DMAG                              XHOST.adefault.getattachedutype<double>("DEFAULT_APL_DMAG")
#define AFLOWRC_DEFAULT_APL_DXYZONLY                          FALSE
#define         DEFAULT_APL_DXYZONLY                          XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DXYZONLY")
#define AFLOWRC_DEFAULT_APL_DSYMMETRIZE                       TRUE
#define         DEFAULT_APL_DSYMMETRIZE                       XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DSYMMETRIZE")
#define AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY                     TRUE
#define         DEFAULT_APL_DINEQUIV_ONLY                     XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DINEQUIV_ONLY")
#define AFLOWRC_DEFAULT_APL_DPM                               string("ON")
#define         DEFAULT_APL_DPM                               XHOST.adefault.getattachedscheme("DEFAULT_APL_DPM")
#define AFLOWRC_DEFAULT_APL_RELAX                             TRUE
#define         DEFAULT_APL_RELAX                             XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_RELAX")
#define AFLOWRC_DEFAULT_APL_RELAX_COMMENSURATE                TRUE  //ME20200427
#define         DEFAULT_APL_RELAX_COMMENSURATE                XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_RELAX_COMMENSURATE")  //ME20200427
#define AFLOWRC_DEFAULT_APL_ZEROSTATE                         FALSE  //CO2018121  //ME20220415 - ZEROSTATE=ON and DPM=ON is unnecessary
#define         DEFAULT_APL_ZEROSTATE                         XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_ZEROSTATE")
#define AFLOWRC_DEFAULT_APL_ZEROSTATE_CHGCAR                  FALSE  //ME20191029
#define         DEFAULT_APL_ZEROSTATE_CHGCAR                  XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_ZEROSTATE_CHGCAR")  //ME20191029
#define AFLOWRC_DEFAULT_APL_USE_LEPSILON                      TRUE
#define         DEFAULT_APL_USE_LEPSILON                      XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_USE_LEPSILON")

//// DEFAULT APL PHONON PROPERTIES
#define AFLOWRC_DEFAULT_APL_FREQFORMAT                        string("THZ|ALLOW_NEGATIVE")  //CO20181226 - no spaces!
#define         DEFAULT_APL_FREQFORMAT                        XHOST.adefault.getattachedscheme("DEFAULT_APL_FREQFORMAT")
#define AFLOWRC_DEFAULT_APL_DC                                TRUE
#define         DEFAULT_APL_DC                                XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DC")
#define AFLOWRC_DEFAULT_APL_DCPATH                            string("LATTICE")
#define         DEFAULT_APL_DCPATH                            XHOST.adefault.getattachedscheme("DEFAULT_APL_DCPATH")
#define AFLOWRC_DEFAULT_APL_DCPOINTS                          100
#define         DEFAULT_APL_DCPOINTS                          XHOST.adefault.getattachedutype<int>("DEFAULT_APL_DCPOINTS")
#define AFLOWRC_DEFAULT_APL_DOS                               TRUE
#define         DEFAULT_APL_DOS                               XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DOS")
#define AFLOWRC_DEFAULT_APL_DOSMETHOD                         string("LT")
#define         DEFAULT_APL_DOSMETHOD                         XHOST.adefault.getattachedscheme("DEFAULT_APL_DOSMETHOD")
#define AFLOWRC_DEFAULT_APL_DOSMESH                           string("21x21x21")
#define         DEFAULT_APL_DOSMESH                           XHOST.adefault.getattachedscheme("DEFAULT_APL_DOSMESH")
#define AFLOWRC_DEFAULT_APL_DOSPOINTS                         2000
#define         DEFAULT_APL_DOSPOINTS                         XHOST.adefault.getattachedutype<int>("DEFAULT_APL_DOSPOINTS")
#define AFLOWRC_DEFAULT_APL_DOSSMEAR                          0.0
#define         DEFAULT_APL_DOSSMEAR                          XHOST.adefault.getattachedutype<double>("DEFAULT_APL_DOSSMEAR")
#define AFLOWRC_DEFAULT_APL_DOS_PROJECT                       FALSE  //ME20200213
#define         DEFAULT_APL_DOS_PROJECT                       XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DOS_PROJECT")  //ME20200213
#define AFLOWRC_DEFAULT_APL_TP                                TRUE
#define         DEFAULT_APL_TP                                XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_TP")
#define AFLOWRC_DEFAULT_APL_DISPLACEMENTS                     TRUE  //ME20200421
#define         DEFAULT_APL_DISPLACEMENTS                     XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DISPLACEMENTS")  //ME20200421
#define AFLOWRC_DEFAULT_APL_TPT                               string("0:2000:10")
#define         DEFAULT_APL_TPT                               XHOST.adefault.getattachedscheme("DEFAULT_APL_TPT")
#define AFLOWRC_DEFAULT_APL_GVEL                              TRUE  // ME20200517
#define         DEFAULT_APL_GVEL                              XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DISPLACEMENTS")  // ME20200517

//// DEFAULT APL FILES
#define AFLOWRC_DEFAULT_APL_FILE_PREFIX                       string("aflow.apl.")
#define         DEFAULT_APL_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_APL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_APL_OUT_FILE                          string("out")  // ME20210927
#define         DEFAULT_APL_OUT_FILE                          XHOST.adefault.getattachedscheme("DEFAULT_APL_OUT_FILE")  // ME20210927
#define AFLOWRC_DEFAULT_APL_PDIS_FILE                         string("phonon_dispersion.out")
#define         DEFAULT_APL_PDIS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_APL_PDIS_FILE")
#define AFLOWRC_DEFAULT_APL_PDOS_FILE                         string("phonon_dos.out")
#define         DEFAULT_APL_PDOS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_APL_PDOS_FILE")
#define AFLOWRC_DEFAULT_APL_THERMO_FILE                       string("thermodynamic_properties.out")
#define         DEFAULT_APL_THERMO_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_APL_THERMO_FILE")
#define AFLOWRC_DEFAULT_APL_THERMO_JSON                       string("thermodynamic_properties.json")  //ME20211019
#define         DEFAULT_APL_THERMO_JSON                       XHOST.adefault.getattachedscheme("DEFAULT_APL_THERMO_JSON")  //ME20211019
#define AFLOWRC_DEFAULT_APL_DYNMAT_FILE                       string("DYNMAT.out")
#define         DEFAULT_APL_DYNMAT_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_APL_DYNMAT_FILE")
#define AFLOWRC_DEFAULT_APL_HARMIFC_FILE                      string("harmonicIFCs.xml")
#define         DEFAULT_APL_HARMIFC_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_APL_HARMIFC_FILE")
//ME20200415
#define AFLOWRC_DEFAULT_APL_POLAR_FILE                        string("polar.xml")
#define         DEFAULT_APL_POLAR_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_APL_POLAR_FILE")
#define AFLOWRC_DEFAULT_APL_HSKPTS_FILE                       string("hskpoints.out")
#define         DEFAULT_APL_HSKPTS_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_APL_HSKPTS_FILE")
#define AFLOWRC_DEFAULT_APL_MSQRDISP_FILE                     string("displacements.out")  //ME20200329
#define         DEFAULT_APL_MSQRDISP_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_APL_MSQRDISP_FILE")  // ME20200329
#define AFLOWRC_DEFAULT_APL_GVEL_FILE                         string("group_velocities.out")  //ME202000517
#define         DEFAULT_APL_GVEL_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_APL_GVEL_FILE")  //ME20200517
//ME20190614 BEGIN
#define AFLOWRC_DEFAULT_APL_PHDOSCAR_FILE                     string("PHDOSCAR")
#define         DEFAULT_APL_PHDOSCAR_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_APL_PHDOSCAR_FILE")
#define AFLOWRC_DEFAULT_APL_PHPOSCAR_FILE                     string("PHPOSCAR")
#define         DEFAULT_APL_PHPOSCAR_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_APL_PHPOSCAR_FILE")
#define AFLOWRC_DEFAULT_APL_PHKPOINTS_FILE                    string("PHKPOINTS")
#define         DEFAULT_APL_PHKPOINTS_FILE                    XHOST.adefault.getattachedscheme("DEFAULT_APL_PHKPOINTS_FILE")
#define AFLOWRC_DEFAULT_APL_PHEIGENVAL_FILE                   string("PHEIGENVAL")
#define         DEFAULT_APL_PHEIGENVAL_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_APL_PHEIGENVAL_FILE")
//ME20190614 END
#define AFLOWRC_DEFAULT_APL_STATE_FILE                        string("fccalc_state.out")
#define         DEFAULT_APL_STATE_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_APL_STATE_FILE")  //ME20200224

//ME20200329 BEGIN
#define AFLOWRC_DEFAULT_APL_ADISP_SCENE_FORMAT                string("XCRYSDEN")
#define         DEFAULT_APL_ADISP_SCENE_FORMAT                XHOST.adefault.getattachedscheme("DEFAULT_APL_ADISP_SCENE_FORMAT")
#define AFLOWRC_DEFAULT_APL_ADISP_AMPLITUDE                   0.25
#define         DEFAULT_APL_ADISP_AMPLITUDE                   XHOST.adefault.getattachedutype<double>("DEFAULT_APL_ADISP_AMPLITUDE")
#define AFLOWRC_DEFAULT_APL_ADISP_NSTEPS                      20
#define         DEFAULT_APL_ADISP_NSTEPS                      XHOST.adefault.getattachedutype<int>("DEFAULT_APL_ADISP_NSTEPS")
#define AFLOWRC_DEFAULT_APL_ADISP_NPERIODS                    1
#define         DEFAULT_APL_ADISP_NPERIODS                    XHOST.adefault.getattachedutype<int>("DEFAULT_APL_ADISP_NPERIODS")
//ME20190614 END

// DEFAULT QHA
//// DEFAULT QHA VALUES
#define AFLOWRC_DEFAULT_QHA_MODE                              string("QHA")
#define         DEFAULT_QHA_MODE                              XHOST.adefault.getattachedscheme("DEFAULT_QHA_MODE")
#define AFLOWRC_DEFAULT_QHA_EOS                               TRUE
#define         DEFAULT_QHA_EOS                               XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_EOS")
#define AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE              string("-12:16:3")
#define         DEFAULT_QHA_EOS_DISTORTION_RANGE              XHOST.adefault.getattachedscheme("DEFAULT_QHA_EOS_DISTORTION_RANGE")
//AS20200818 BEGIN
#define AFLOWRC_DEFAULT_QHA_EOS_MODEL                         string("SJ")
#define         DEFAULT_QHA_EOS_MODEL                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_EOS_MODEL")
//AS20200818 END
#define AFLOWRC_DEFAULT_QHA_GP_DISTORTION                     1.0
#define         DEFAULT_QHA_GP_DISTORTION                     XHOST.adefault.getattachedutype<double>("DEFAULT_QHA_GP_DISTORTION")
//AS20200602 BEGIN
#define AFLOWRC_DEFAULT_QHA_TAYLOR_EXPANSION_ORDER            2
#define         DEFAULT_QHA_TAYLOR_EXPANSION_ORDER            XHOST.adefault.getattachedutype<double>("DEFAULT_QHA_TAYLOR_EXPANSION_ORDER")
//AS20200602 END
#define AFLOWRC_DEFAULT_QHA_INCLUDE_ELEC_CONTRIB              FALSE
#define         DEFAULT_QHA_INCLUDE_ELEC_CONTRIB              XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_INCLUDE_ELEC_CONTRIB")
//AS20200528 BEGIN
#define AFLOWRC_DEFAULT_QHA_SOMMERFELD_EXPANSION              FALSE
#define         DEFAULT_QHA_SOMMERFELD_EXPANSION              XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_SOMMERFELD_EXPANSION")
//AS20200528 END
#define AFLOWRC_DEFAULT_QHA_PDIS_T                            string("300")
#define         DEFAULT_QHA_PDIS_T                            XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDIS_T")
//AS20200508 BEGIN
#define AFLOWRC_DEFAULT_QHA_GP_FINITE_DIFF                    FALSE
#define         DEFAULT_QHA_GP_FINITE_DIFF                    XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_GP_FINITE_DIFF")
#define AFLOWRC_DEFAULT_QHA_IGNORE_IMAGINARY                  FALSE
#define         DEFAULT_QHA_IGNORE_IMAGINARY                  XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_IGNORE_IMAGINARY")
//AS20201123 BEGIN
#define AFLOWRC_DEFAULT_QHA_RELAX_IONS_CELL                   FALSE
#define         DEFAULT_QHA_RELAX_IONS_CELL                   XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_RELAX_IONS_CELL")
//AS20201123 END

//// DEFAULT QHA FILES
#define AFLOWRC_DEFAULT_QHA_FILE_PREFIX                       string("aflow.qha.")
#define         DEFAULT_QHA_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_QHA_FILE_PREFIX")
//AS20200709 BEGIN
#define AFLOWRC_DEFAULT_QHA3P_FILE_PREFIX                     string("aflow.qha3p.")
#define         DEFAULT_QHA3P_FILE_PREFIX                     XHOST.adefault.getattachedscheme("DEFAULT_QHA3P_FILE_PREFIX")
#define AFLOWRC_DEFAULT_QHANP_FILE_PREFIX                     string("aflow.qhanp.")
#define         DEFAULT_QHANP_FILE_PREFIX                     XHOST.adefault.getattachedscheme("DEFAULT_QHANP_FILE_PREFIX")
#define AFLOWRC_DEFAULT_SCQHA_FILE_PREFIX                     string("aflow.scqha.")
#define         DEFAULT_SCQHA_FILE_PREFIX                     XHOST.adefault.getattachedscheme("DEFAULT_SCQHA_FILE_PREFIX")
//AS20200709 END
#define AFLOWRC_DEFAULT_QHA_GP_PATH_FILE                      string("gp.disp.out")
#define         DEFAULT_QHA_GP_PATH_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_PATH_FILE")
#define AFLOWRC_DEFAULT_QHA_GP_MESH_FILE                      string("gp.mesh.out")
#define         DEFAULT_QHA_GP_MESH_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_MESH_FILE")
#define AFLOWRC_DEFAULT_QHA_GP_AVG_FILE                       string("gp.avg.out")
#define         DEFAULT_QHA_GP_AVG_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_AVG_FILE")
#define AFLOWRC_DEFAULT_QHA_THERMO_FILE                       string("thermo.out")
#define         DEFAULT_QHA_THERMO_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_QHA_THERMO_FILE")
#define AFLOWRC_DEFAULT_QHA_FREQS_FILE                        string("frequencies.out")
#define         DEFAULT_QHA_FREQS_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_QHA_FREQS_FILE")
#define AFLOWRC_DEFAULT_QHA_FVT_FILE                          string("FVT.out")
#define         DEFAULT_QHA_FVT_FILE                          XHOST.adefault.getattachedscheme("DEFAULT_QHA_FVT_FILE")
//AS20200508 END
//AS20210517 BEGIN
#define AFLOWRC_DEFAULT_QHA_COEFF_FILE                        string("coeff.out")
#define         DEFAULT_QHA_COEFF_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_QHA_COEFF_FILE")
#define AFLOWRC_DEFAULT_QHA_IMAG_FILE                         string("imag.out")
#define         DEFAULT_QHA_IMAG_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_IMAG_FILE")
//AS20210517 END
//AS20201022 BEGIN
#define AFLOWRC_DEFAULT_QHA_PDIS_FILE                         string("dispersion_phonon")
#define         DEFAULT_QHA_PDIS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDIS_FILE")
//AS20201022 END
//AS20201201 BEGIN
#define AFLOWRC_DEFAULT_QHA_PDOS_FILE                         string("dos_phonon")
#define         DEFAULT_QHA_PDOS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDOS_FILE")
//AS20201201 END
//AS20201112 BEGIN
#define AFLOWRC_DEFAULT_QHA_KPOINTS_FILE                      string("kpoints.out")
#define         DEFAULT_QHA_KPOINTS_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_QHA_KPOINTS_FILE")
//AS20201112 
//AS20210914 BEGIN
#define AFLOWRC_DEFAULT_POCC_QHA_THERMO_FILE                  string("qha.thermo.out")
#define         DEFAULT_POCC_QHA_THERMO_FILE                  XHOST.adefault.getattachedscheme("DEFAULT_POCC_QHA_THERMO_FILE")
#define AFLOWRC_DEFAULT_POCC_QHA_AVGTHERMO_FILE               string("qha.avgthermo.out")
#define         DEFAULT_POCC_QHA_AVGTHERMO_FILE               XHOST.adefault.getattachedscheme("DEFAULT_POCC_QHA_AVGTHERMO_FILE")
//AS20210914 END

// DEFAULT AAPL
//// DEFAULT AAPL VALUES
#define AFLOWRC_DEFAULT_AAPL_BTE                              string("FULL")
#define         DEFAULT_AAPL_BTE                              XHOST.adefault.getattachedscheme("DEFAULT_AAPL_BTE")
#define AFLOWRC_DEFAULT_AAPL_FOURTH_ORDER                     FALSE
#define         DEFAULT_AAPL_FOURTH_ORDER                     XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_FOURTH_ORDER")
#define AFLOWRC_DEFAULT_AAPL_CUT_RAD                          string("0.0") //ME20190308 - use CUT_SHELL by default //ME20191029
#define         DEFAULT_AAPL_CUT_RAD                          XHOST.adefault.getattachedscheme("DEFAULT_AAPL_CUT_RAD")
#define AFLOWRC_DEFAULT_AAPL_CUT_SHELL                        string("6")  //ME20190301  //ME20190408  //ME20191029
#define         DEFAULT_AAPL_CUT_SHELL                        XHOST.adefault.getattachedscheme("DEFAULT_AAPL_CUT_SHELL")
#define AFLOWRC_DEFAULT_AAPL_THERMALGRID                      string("21x21x21")  //ME20200110 - 21x21x21 more than enough for tetrahedron method; odd is preferred (Gamma-centered by construction)
#define         DEFAULT_AAPL_THERMALGRID                      XHOST.adefault.getattachedscheme("DEFAULT_AAPL_THERMALGRID")
#define AFLOWRC_DEFAULT_AAPL_TCT                              string("50:550:50")
#define         DEFAULT_AAPL_TCT                              XHOST.adefault.getattachedscheme("DEFAULT_AAPL_TCT")
#define AFLOWRC_DEFAULT_AAPL_SUMRULE                          1e-7
#define         DEFAULT_AAPL_SUMRULE                          XHOST.adefault.getattachedutype<double>("DEFAULT_AAPL_SUMRULE")
#define AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER                 2000
#define         DEFAULT_AAPL_SUMRULE_MAX_ITER                 XHOST.adefault.getattachedutype<int>("DEFAULT_AAPL_SUMRULE_MAX_ITER")
#define AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT               0.0
#define         DEFAULT_AAPL_MIXING_COEFFICIENT               XHOST.adefault.getattachedutype<double>("DEFAULT_AAPL_MIXING_COEFFICIENT")
#define AFLOWRC_DEFAULT_AAPL_ISOTOPE                          TRUE
#define         DEFAULT_AAPL_ISOTOPE                          XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_ISOTOPE")
#define AFLOWRC_DEFAULT_AAPL_BOUNDARY                         FALSE
#define         DEFAULT_AAPL_BOUNDARY                         XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_BOUNDARY")
#define AFLOWRC_DEFAULT_AAPL_CUMULATIVEK                      FALSE
#define         DEFAULT_AAPL_CUMULATIVEK                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_CUMULATIVEK")
#define AFLOWRC_DEFAULT_AAPL_NANO_SIZE                        100.0
#define         DEFAULT_AAPL_NANO_SIZE                        XHOST.adefault.getattachedutype<double>("DEFAULT_AAPL_NANO_SIZE")

//// DEFAULT AAPL FILES
#define AFLOWRC_DEFAULT_AAPL_FILE_PREFIX                      string("aflow.aapl.")
#define         DEFAULT_AAPL_FILE_PREFIX                      XHOST.adefault.getattachedscheme("DEFAULT_AAPL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_AAPL_IRRQPTS_FILE                     string("irred_qpoints.out")
#define         DEFAULT_AAPL_IRRQPTS_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_AAPL_IRRQPTS_FILE")
#define AFLOWRC_DEFAULT_AAPL_GVEL_FILE                        string("group_velocities.out")
#define         DEFAULT_AAPL_GVEL_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_AAPL_GVEL_FILE")
#define AFLOWRC_DEFAULT_AAPL_PS_FILE                          string("phase_space.out")  //ME20191104
#define         DEFAULT_AAPL_PS_FILE                          XHOST.adefault.getattachedscheme("DEFAULT_AAPL_PS_FILE")  //ME20191104
#define AFLOWRC_DEFAULT_AAPL_GRUENEISEN_FILE                  string("grueneisen.out")  //ME20191104
#define         DEFAULT_AAPL_GRUENEISEN_FILE                  XHOST.adefault.getattachedscheme("DEFAULT_AAPL_GRUENEISEN_FILE")  //ME20191104
#define AFLOWRC_DEFAULT_AAPL_RATES_FILE                       string("scattering_rates_total.out")
#define         DEFAULT_AAPL_RATES_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_FILE")
#define AFLOWRC_DEFAULT_AAPL_RATES_3RD_FILE                   string("scattering_rates_anharmonic_3rd.out")
#define         DEFAULT_AAPL_RATES_3RD_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_3RD_FILE")
#define AFLOWRC_DEFAULT_AAPL_RATES_4TH_FILE                   DEFAULT_AAPL_FILE_PREFIX+string("scattering_rates_anharmonic_4th.out")
#define         DEFAULT_AAPL_RATES_4TH_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_4TH_FILE")
#define AFLOWRC_DEFAULT_AAPL_ISOTOPE_FILE                     string("scattering_rates_isotope.out")
#define         DEFAULT_AAPL_ISOTOPE_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_AAPL_ISOTOPE_FILE")
#define AFLOWRC_DEFAULT_AAPL_BOUNDARY_FILE                    string("scattering_rates_boundary.out")
#define         DEFAULT_AAPL_BOUNDARY_FILE                    XHOST.adefault.getattachedscheme("DEFAULT_AAPL_BOUNDARY_FILE")
#define AFLOWRC_DEFAULT_AAPL_TCOND_FILE                       string("thermal_conductivity.out")
#define         DEFAULT_AAPL_TCOND_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_AAPL_TCOND_FILE")


// DEFAULT AEL
//// DEFAULT AEL STRAIN CALCS
#define AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY                   TRUE
#define         DEFAULT_AEL_STRAIN_SYMMETRY                   XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_STRAIN_SYMMETRY")
#define AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS                   4
#define         DEFAULT_AEL_NNORMAL_STRAINS                   XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_NNORMAL_STRAINS")
#define AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS                    4
#define         DEFAULT_AEL_NSHEAR_STRAINS                    XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_NSHEAR_STRAINS")
#define AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP                0.005
#define         DEFAULT_AEL_NORMAL_STRAIN_STEP                XHOST.adefault.getattachedutype<double>("DEFAULT_AEL_NORMAL_STRAIN_STEP")
#define AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP                 0.005
#define         DEFAULT_AEL_SHEAR_STRAIN_STEP                 XHOST.adefault.getattachedutype<double>("DEFAULT_AEL_SHEAR_STRAIN_STEP")
#define AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC                FALSE
#define         DEFAULT_AEL_ORIGIN_STRAIN_CALC                XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_ORIGIN_STRAIN_CALC")
#define AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT                 FALSE
#define         DEFAULT_AEL_ORIGIN_STRAIN_FIT                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_ORIGIN_STRAIN_FIT")
#define AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT                FALSE
#define         DEFAULT_AEL_RELAXED_STRUCT_FIT                XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_RELAXED_STRUCT_FIT")
#define AFLOWRC_DEFAULT_AEL_NEG_STRAINS                       TRUE
#define         DEFAULT_AEL_NEG_STRAINS                       XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_NEG_STRAINS")
#define AFLOWRC_DEFAULT_AEL_NIND_STRAIN_DIRS                  3
#define         DEFAULT_AEL_NIND_STRAIN_DIRS                  XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_NIND_STRAIN_DIRS")
#define AFLOWRC_DEFAULT_AEL_VASPSYM                           FALSE
#define         DEFAULT_AEL_VASPSYM                           XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_VASPSYM")
#define AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM                  FALSE
#define         DEFAULT_AEL_PRECACC_ALGONORM                  XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_PRECACC_ALGONORM")
#define AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS                 FALSE
#define         DEFAULT_AEL_VASPRUNXML_STRESS                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_VASPRUNXML_STRESS")
#define AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS             FALSE
#define         DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS             XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS")
#define AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX                    1
#define         DEFAULT_AEL_SKIP_ARUNS_MAX                    XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_SKIP_ARUNS_MAX")

//// DEFAULT AEL CHECKS AND PROCESSING
#define AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY            TRUE
#define         DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY            XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY")
#define AFLOWRC_DEFAULT_AEL_SYMMETRIZE                        FALSE
#define         DEFAULT_AEL_SYMMETRIZE                        XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_SYMMETRIZE")

//// DEFAULT AEL OUTPUT FILES
#define AFLOWRC_DEFAULT_AEL_FILE_PREFIX                       string("aflow.ael.")
#define         DEFAULT_AEL_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_AEL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS                FALSE
#define         DEFAULT_AEL_WRITE_FULL_RESULTS                XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_WRITE_FULL_RESULTS")
#define AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN                      TRUE
#define         DEFAULT_AEL_DIRNAME_ARUN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_DIRNAME_ARUN")

// DEFAULT AGL
//// DEFAULT AGL STRAIN CALCS
#define AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO                 TRUE
#define         DEFAULT_AGL_AEL_POISSON_RATIO                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_AEL_POISSON_RATIO")
#define AFLOWRC_DEFAULT_AGL_NSTRUCTURES                       28
#define         DEFAULT_AGL_NSTRUCTURES                       XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_NSTRUCTURES")
#define AFLOWRC_DEFAULT_AGL_STRAIN_STEP                       0.01
#define         DEFAULT_AGL_STRAIN_STEP                       XHOST.adefault.getattachedutype<double>("DEFAULT_AGL_STRAIN_STEP")
#define AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS             FALSE
#define         DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS             XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS")
#define AFLOWRC_DEFAULT_AGL_SKIP_ARUNS_MAX                    7
#define         DEFAULT_AGL_SKIP_ARUNS_MAX                    XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_SKIP_ARUNS_MAX")

//// DEFAULT AGL CHECKS AND PROCESSING
#define AFLOWRC_DEFAULT_AGL_NTEMPERATURE                      201
#define         DEFAULT_AGL_NTEMPERATURE                      XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_NTEMPERATURE")
#define AFLOWRC_DEFAULT_AGL_STEMPERATURE                      10.0
#define         DEFAULT_AGL_STEMPERATURE                      XHOST.adefault.getattachedutype<double>("DEFAULT_AGL_STEMPERATURE")
#define AFLOWRC_DEFAULT_AGL_NPRESSURE                         101
#define         DEFAULT_AGL_NPRESSURE                         XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_NPRESSURE")
#define AFLOWRC_DEFAULT_AGL_SPRESSURE                         1.0
#define         DEFAULT_AGL_SPRESSURE                         XHOST.adefault.getattachedutype<double>("DEFAULT_AGL_SPRESSURE")
#define AFLOWRC_DEFAULT_AGL_POISSON_RATIO                     0.25
#define         DEFAULT_AGL_POISSON_RATIO                     XHOST.adefault.getattachedutype<double>("DEFAULT_AGL_POISSON_RATIO")
#define AFLOWRC_DEFAULT_AGL_IEOS                              0
#define         DEFAULT_AGL_IEOS                              XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_IEOS")
#define AFLOWRC_DEFAULT_AGL_IDEBYE                            0
#define         DEFAULT_AGL_IDEBYE                            XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_IDEBYE")
#define AFLOWRC_DEFAULT_AGL_FIT_TYPE                          0
#define         DEFAULT_AGL_FIT_TYPE                          XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_FIT_TYPE")
#define AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY                FALSE
#define         DEFAULT_AGL_CHECK_EV_CONCAVITY                XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_CHECK_EV_CONCAVITY")
#define AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN                      FALSE
#define         DEFAULT_AGL_CHECK_EV_MIN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_CHECK_EV_MIN")
#define AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC                     TRUE
#define         DEFAULT_AGL_HUGONIOT_CALC                     XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_HUGONIOT_CALC")
#define AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE              FALSE
#define         DEFAULT_AGL_HUGONIOT_EXTRAPOLATE              XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_HUGONIOT_EXTRAPOLATE")
#define AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE      FALSE
#define         DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE")

//// DEFAULT AGL OUTPUT FILES
#define AFLOWRC_DEFAULT_AGL_FILE_PREFIX                       string("aflow.agl.")
#define         DEFAULT_AGL_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_AGL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS                FALSE
#define         DEFAULT_AGL_WRITE_FULL_RESULTS                XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_WRITE_FULL_RESULTS")
#define AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN                      TRUE
#define         DEFAULT_AGL_DIRNAME_ARUN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_DIRNAME_ARUN")
#define AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT                 FALSE
#define         DEFAULT_AGL_WRITE_GIBBS_INPUT                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_WRITE_GIBBS_INPUT")
#define AFLOWRC_DEFAULT_AGL_PLOT_RESULTS                      FALSE
#define         DEFAULT_AGL_PLOT_RESULTS                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_PLOT_RESULTS")

//// DEFAULT QCA
#define AFLOWRC_DEFAULT_QCA_MIN_SLEEP_SECONDS                 60 // seconds
#define         DEFAULT_QCA_MIN_SLEEP_SECONDS                 XHOST.adefault.getattachedutype<int>("DEFAULT_QCA_MIN_SLEEP_SECONDS")
#define AFLOWRC_DEFAULT_QCA_MAX_NUM_ATOMS                     8
#define         DEFAULT_QCA_MAX_NUM_ATOMS                     XHOST.adefault.getattachedutype<int>("DEFAULT_QCA_MAX_NUM_ATOMS")
#define AFLOWRC_DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS               4
#define         DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS               XHOST.adefault.getattachedutype<int>("DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS")
#define AFLOWRC_DEFAULT_QCA_CV_CUTOFF                         0.05
#define         DEFAULT_QCA_CV_CUTOFF                         XHOST.adefault.getattachedutype<double>("DEFAULT_QCA_CV_CUTOFF")
#define AFLOWRC_DEFAULT_QCA_CONC_NPTS                         20
#define         DEFAULT_QCA_CONC_NPTS                         XHOST.adefault.getattachedutype<double>("DEFAULT_QCA_CONC_NPTS")
#define AFLOWRC_DEFAULT_QCA_TEMP_NPTS                         150
#define         DEFAULT_QCA_TEMP_NPTS                         XHOST.adefault.getattachedutype<double>("DEFAULT_QCA_TEMP_NPTS")
#define AFLOWRC_DEFAULT_QCA_TEMP_MIN                          300 // K
#define         DEFAULT_QCA_TEMP_MIN                          XHOST.adefault.getattachedutype<double>("DEFAULT_QCA_TEMP_MIN")
#define AFLOWRC_DEFAULT_QCA_TEMP_MAX                          5000 // K
#define         DEFAULT_QCA_TEMP_MAX                          XHOST.adefault.getattachedutype<double>("DEFAULT_QCA_TEMP_MAX")
#define AFLOWRC_DEFAULT_QCA_TEMP_MIN_LIMIT                    10000 // K
#define         DEFAULT_QCA_TEMP_MIN_LIMIT                    XHOST.adefault.getattachedutype<double>("DEFAULT_QCA_TEMP_MIN_LIMIT")
#define AFLOWRC_DEFAULT_QCA_PRINT                             string("txt")
#define         DEFAULT_QCA_PRINT                             XHOST.adefault.getattachedscheme("DEFAULT_QCA_PRINT")

//RF20200413 START
// DEFAULT CCE
#define AFLOWRC_DEFAULT_CCE_OX_METHOD                         1
#define         DEFAULT_CCE_OX_METHOD                         XHOST.adefault.getattachedutype<int>("DEFAULT_CCE_OX_METHOD")
#define AFLOWRC_DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION           0.4
#define         DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION           XHOST.adefault.getattachedutype<double>("DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION")
#define AFLOWRC_DEFAULT_CCE_OX_TOL                            0.001
#define         DEFAULT_CCE_OX_TOL                            XHOST.adefault.getattachedutype<double>("DEFAULT_CCE_OX_TOL")
#define AFLOWRC_DEFAULT_CCE_PEROX_CUTOFF                      1.6
#define         DEFAULT_CCE_PEROX_CUTOFF                      XHOST.adefault.getattachedutype<double>("DEFAULT_CCE_PEROX_CUTOFF")
#define AFLOWRC_DEFAULT_CCE_SUPEROX_CUTOFF                    1.4
#define         DEFAULT_CCE_SUPEROX_CUTOFF                    XHOST.adefault.getattachedutype<double>("DEFAULT_CCE_SUPEROX_CUTOFF")
#define AFLOWRC_DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF          1.3
#define         DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF          XHOST.adefault.getattachedutype<double>("DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF")
#define AFLOWRC_DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF          1.2
#define         DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF          XHOST.adefault.getattachedutype<double>("DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF")
//RF20200413 END

// DEFAULT XTALFINDER
#define AFLOWRC_DEFAULT_XTALFINDER_MISFIT_MATCH               0.1 // values below this threshold: similar structures have similar properties // DX20201118
#define         DEFAULT_XTALFINDER_MISFIT_MATCH               XHOST.adefault.getattachedutype<double>("DEFAULT_XTALFINDER_MISFIT_MATCH") //DX20201118
#define AFLOWRC_DEFAULT_XTALFINDER_MISFIT_FAMILY              0.2 // values above this threshold: matched structures do not have similar properties //DX20201118
#define         DEFAULT_XTALFINDER_MISFIT_FAMILY              XHOST.adefault.getattachedutype<double>("DEFAULT_XTALFINDER_MISFIT_FAMILY") //DX20201118
#define AFLOWRC_DEFAULT_XTALFINDER_SUPERCELL_METHOD           FALSE // supercell method for comparing (robust, but slow, superceded by transformation method)
#define         DEFAULT_XTALFINDER_SUPERCELL_METHOD           XHOST.adefault.getattachedutype<bool>("DEFAULT_XTALFINDER_SUPERCELL_METHOD") //DX20201223
//DX20200709 - START
#define AFLOWRC_DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING    4.0 // factor that divides minimum interatomic distance
#define         DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING    XHOST.adefault.getattachedutype<double>("DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING")
//DX20200709 - END
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL                    string("material_comparison_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_MATERIAL                    XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_MATERIAL") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE                   string("structure_comparison_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_STRUCTURE                   XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_STRUCTURE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_DUPLICATE                   string("duplicate_compounds_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_DUPLICATE                   XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_DUPLICATE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE   string("material_comparison_compare2database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE   XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE  string("structure_comparison_compare2database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE  XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE           string("material_comparison_database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE           XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE          string("structure_comparison_database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE          XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE") //DX20201228

//DX20200720 - START
// DEFAULT ANRL
#define AFLOWRC_DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL           1e-6 // tolerance for equivalent Wyckoff coordinates
#define         DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL           XHOST.adefault.getattachedutype<double>("DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL")
//DX20200720 - END

// CORES // DONE
#define AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP                   56.0    // Celsius
#define         AFLOW_CORE_TEMPERATURE_BEEP                   XHOST.adefault.getattachedutype<double>("AFLOW_CORE_TEMPERATURE_BEEP") 
#define AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT                   65.0    // Celsius, you need to run aflow as root to shutdown
#define         AFLOW_CORE_TEMPERATURE_HALT                   XHOST.adefault.getattachedutype<double>("AFLOW_CORE_TEMPERATURE_HALT") 
#define AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH                5.0    // seconds
#define         AFLOW_CORE_TEMPERATURE_REFRESH                XHOST.adefault.getattachedutype<double>("AFLOW_CORE_TEMPERATURE_REFRESH") 

// VASP MACHINE SETTINGS
#define AFLOWRC_SECONDS_SLEEP_VASP_COMPLETION                 15    // seconds
#define         SECONDS_SLEEP_VASP_COMPLETION                 XHOST.adefault.getattachedutype<double>("SECONDS_SLEEP_VASP_COMPLETION") 
#define AFLOWRC_SECONDS_SLEEP_VASP_MONITOR                    60    // seconds
#define         SECONDS_SLEEP_VASP_MONITOR                    XHOST.adefault.getattachedutype<double>("SECONDS_SLEEP_VASP_MONITOR") 
#define AFLOWRC_SECONDS_STALE_OUTCAR                          21600    // seconds
#define         SECONDS_STALE_OUTCAR                          XHOST.adefault.getattachedutype<double>("SECONDS_STALE_OUTCAR") 
#define AFLOWRC_BYTES_MAX_VASP_OUT                            20000000000    // bytes
#define         BYTES_MAX_VASP_OUT                            XHOST.adefault.getattachedutype<unsigned long long int>("BYTES_MAX_VASP_OUT") 
#define AFLOWRC_MEMORY_MAX_USAGE_RAM                          98    // percent
#define         MEMORY_MAX_USAGE_RAM                          XHOST.adefault.getattachedutype<double>("MEMORY_MAX_USAGE_RAM") 
#define AFLOWRC_MEMORY_MAX_USAGE_SWAP                         45    // percent  //shouldn't go above 50, sometimes it ramps up quickly, so set to 45 to be safe
#define         MEMORY_MAX_USAGE_SWAP                         XHOST.adefault.getattachedutype<double>("MEMORY_MAX_USAGE_SWAP") 
#define AFLOWRC_FILE_VASP_MONITOR                             string("monitor_vasp")
#define         FILE_VASP_MONITOR                             XHOST.adefault.getattachedscheme("FILE_VASP_MONITOR")
#define AFLOWRC_INTEL_COMPILER_PATHS                          string("/opt/intel/bin/compilervars.sh,/opt/intel/bin/compilervars.csh,/app/intel/parallel_studio_xe_2020_update1/bin/compilervars.sh")
#define         INTEL_COMPILER_PATHS                          XHOST.adefault.getattachedscheme("INTEL_COMPILER_PATHS")

// MACHINE DEPENDENT MPI
#define AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH                   string("ulimit -s unlimited ") // DUKE_BETA_MPICH
#define         MPI_OPTIONS_DUKE_BETA_MPICH                   XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_MPICH")
#define AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH                   string("/usr/bin/mpiexec -np") // DUKE_BETA_MPICH
#define         MPI_COMMAND_DUKE_BETA_MPICH                   XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_MPICH")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH                string("/usr/local/bin/") // DUKE_BETA_MPICH
#define         MPI_BINARY_DIR_DUKE_BETA_MPICH                XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_MPICH")

#define AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI                 string("ulimit -s unlimited ") // DUKE_BETA_OPENMPI
#define         MPI_OPTIONS_DUKE_BETA_OPENMPI                 XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_OPENMPI")
#define AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI                 string("/usr/bin/mpirun.openmpi -np") // DUKE_BETA_OPENMPI
#define         MPI_COMMAND_DUKE_BETA_OPENMPI                 XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_OPENMPI")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI              string("/usr/local/bin/") // DUKE_BETA_OPENMPI
#define         MPI_BINARY_DIR_DUKE_BETA_OPENMPI              XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_OPENMPI")

#define AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS                    string("ulimit -s unlimited ") // DUKE_MATERIALS
#define         MPI_OPTIONS_DUKE_MATERIALS                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_MATERIALS")
#define AFLOWRC_MPI_COMMAND_DUKE_MATERIALS                    string("/usr/bin/mpiexec -np") // DUKE_MATERIALS
#define         MPI_COMMAND_DUKE_MATERIALS                    XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_MATERIALS")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS                 string("/usr/local/bin/")  // DUKE_MATERIALS
#define         MPI_BINARY_DIR_DUKE_MATERIALS                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_MATERIALS")

#define AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB                     string("ulimit -s unlimited ") // DUKE_AFLOWLIB
#define         MPI_OPTIONS_DUKE_AFLOWLIB                     XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_AFLOWLIB")
#define AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB                     string("/usr/bin/mpiexec -np") // DUKE_AFLOWLIB
#define         MPI_COMMAND_DUKE_AFLOWLIB                     XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_AFLOWLIB")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB                  string("/usr/local/bin/") // DUKE_AFLOWLIB
#define         MPI_BINARY_DIR_DUKE_AFLOWLIB                  XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_AFLOWLIB")

#define AFLOWRC_MPI_OPTIONS_DUKE_HABANA                       string("ulimit -s unlimited ") // DUKE_HABANA
#define         MPI_OPTIONS_DUKE_HABANA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_HABANA")
#define AFLOWRC_MPI_COMMAND_DUKE_HABANA                       string("/usr/bin/mpiexec -np") // DUKE_HABANA
#define         MPI_COMMAND_DUKE_HABANA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_HABANA")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA                    string("/usr/local/bin/") // DUKE_HABANA
#define         MPI_BINARY_DIR_DUKE_HABANA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_HABANA")

#define AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH                  string("ulimit -s unlimited ") // DUKE_QRATS_MPICH
#define         MPI_OPTIONS_DUKE_QRATS_MPICH                  XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QRATS_MPICH")
#define AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH                  string("/MAIN/bin/MPICH/bin/mpirun -np") // DUKE_QRATS_MPICH
#define         MPI_COMMAND_DUKE_QRATS_MPICH                  XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QRATS_MPICH")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH               string("/usr/local/bin/") // DUKE_QRATS_MPICH
#define         MPI_BINARY_DIR_DUKE_QRATS_MPICH               XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QRATS_MPICH")

#define AFLOWRC_MPI_OPTIONS_DUKE_QFLOW_OPENMPI                string("ulimit -s unlimited ") // DUKE_QFLOW_MPICH
#define         MPI_OPTIONS_DUKE_QFLOW_OPENMPI                XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QFLOW_OPENMPI")
#define AFLOWRC_MPI_COMMAND_DUKE_QFLOW_OPENMPI                string("/home/bin/local/bin/mpirun -n") // DUKE_QFLOW_MPICH
#define         MPI_COMMAND_DUKE_QFLOW_OPENMPI                XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QFLOW_OPENMPI")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI             string("/home/bin/") // DUKE_QFLOW_MPICH
#define         MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI             XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI")

//CO20201220 X START
#define AFLOWRC_MPI_OPTIONS_DUKE_X_X                          string("ulimit -s unlimited ") // DUKE_X_X_MPICH
#define         MPI_OPTIONS_DUKE_X_X                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_X")
#define AFLOWRC_MPI_COMMAND_DUKE_X_X                          string("srun --mpi=pmix --cpus-per-task") // DUKE_X_X_MPICH
#define         MPI_COMMAND_DUKE_X_X                          XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_X")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_X                       string("~/bin/x/") // DUKE_X_X_MPICH
#define         MPI_BINARY_DIR_DUKE_X_X                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_X")
#define AFLOWRC_MPI_OPTIONS_DUKE_X_CRAY                       string("ulimit -s unlimited ") // DUKE_X_CRAY_MPICH
#define         MPI_OPTIONS_DUKE_X_CRAY                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_CRAY")
#define AFLOWRC_MPI_COMMAND_DUKE_X_CRAY                       string("srun --mpi=pmix --cpus-per-task") // DUKE_X_CRAY_MPICH
#define         MPI_COMMAND_DUKE_X_CRAY                       XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_CRAY")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_CRAY                    string("~/bin/cray/") // DUKE_X_CRAY_MPICH
#define         MPI_BINARY_DIR_DUKE_X_CRAY                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_CRAY")
#define AFLOWRC_MPI_OPTIONS_DUKE_X_OLDCRAY                    string("ulimit -s unlimited ") // DUKE_X_OLDCRAY_MPICH
#define         MPI_OPTIONS_DUKE_X_OLDCRAY                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_OLDCRAY")
#define AFLOWRC_MPI_COMMAND_DUKE_X_OLDCRAY                    string("srun --mpi=pmix --cpus-per-task") // DUKE_X_OLDCRAY_MPICH
#define         MPI_COMMAND_DUKE_X_OLDCRAY                    XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_OLDCRAY")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_OLDCRAY                 string("~/bin/oldcray/") // DUKE_X_OLDCRAY_MPICH
#define         MPI_BINARY_DIR_DUKE_X_OLDCRAY                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_OLDCRAY")
#define AFLOWRC_MPI_OPTIONS_DUKE_X_SMB                        string("ulimit -s unlimited ") // DUKE_X_SMB_MPICH
#define         MPI_OPTIONS_DUKE_X_SMB                        XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_SMB")
#define AFLOWRC_MPI_COMMAND_DUKE_X_SMB                        string("srun --mpi=pmix --cpus-per-task") // DUKE_X_SMB_MPICH
#define         MPI_COMMAND_DUKE_X_SMB                        XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_SMB")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_SMB                     string("~/bin/smb/") // DUKE_X_SMB_MPICH
#define         MPI_BINARY_DIR_DUKE_X_SMB                     XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_SMB")
//CO20201220 X STOP

//CO20220818 JHU_ROCKFISH START
#define AFLOWRC_MPI_OPTIONS_JHU_ROCKFISH                      string("ulimit -s unlimited ") // JHU_ROCKFISH_MPICH
#define         MPI_OPTIONS_JHU_ROCKFISH                      XHOST.adefault.getattachedscheme("MPI_OPTIONS_JHU_ROCKFISH")
#define AFLOWRC_MPI_COMMAND_JHU_ROCKFISH                      string("mpirun -n") // JHU_ROCKFISH_MPICH
#define         MPI_COMMAND_JHU_ROCKFISH                      XHOST.adefault.getattachedscheme("MPI_COMMAND_JHU_ROCKFISH")
#define AFLOWRC_MPI_BINARY_DIR_JHU_ROCKFISH                   string("~/bin/") // JHU_ROCKFISH_MPICH
#define         MPI_BINARY_DIR_JHU_ROCKFISH                   XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_JHU_ROCKFISH")
//CO20220818 JHU_ROCKFISH STOP

//DX20190509 - MACHINE001 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE001                        string("") // MACHINE001
#define         MPI_OPTIONS_MACHINE001                        XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE001")
#define AFLOWRC_MPI_COMMAND_MACHINE001                        string("aprun -n") // MACHINE001
#define         MPI_COMMAND_MACHINE001                        XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE001")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE001                     string("~/bin/") // MACHINE001
#define         MPI_BINARY_DIR_MACHINE001                     XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE001")
//DX20190509 - MACHINE001 - END

//DX20190509 - MACHINE002 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE002                       string("") // MACHINE002
#define         MPI_OPTIONS_MACHINE002                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE002")
#define AFLOWRC_MPI_COMMAND_MACHINE002                       string("mpirun -np") // MACHINE002
#define         MPI_COMMAND_MACHINE002                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE002")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE002                    string("~/bin/") // MACHINE002
#define         MPI_BINARY_DIR_MACHINE002                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE002")
//DX20190509 - MACHINE002 - END

//DX20201005 - MACHINE003 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE003                       string("") // MACHINE003
#define         MPI_OPTIONS_MACHINE003                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE003")
#define AFLOWRC_MPI_COMMAND_MACHINE003                       string("mpirun -np") // MACHINE003
#define         MPI_COMMAND_MACHINE003                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE003")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE003                    string("~/bin/") // MACHINE003
#define         MPI_BINARY_DIR_MACHINE003                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE003")
//DX20201005 - MACHINE003 - END

//DX20211011 - MACHINE004 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE004                       string("") // MACHINE004
#define         MPI_OPTIONS_MACHINE004                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE004")
#define AFLOWRC_MPI_COMMAND_MACHINE004                       string("mpirun -np") // MACHINE004
#define         MPI_COMMAND_MACHINE004                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE004")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE004                    string("~/bin/") // MACHINE004
#define         MPI_BINARY_DIR_MACHINE004                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE004")
//DX20211011 - MACHINE004 - END

#define AFLOWRC_MPI_OPTIONS_MPCDF_EOS                         string("ulimit -s unlimited ") // MPCDF_EOS_MPICH
#define         MPI_OPTIONS_MPCDF_EOS                         XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_EOS")
#define AFLOWRC_MPI_COMMAND_MPCDF_EOS                         string("/usr/bin/srun -n") // MPCDF_EOS_MPICH
#define         MPI_COMMAND_MPCDF_EOS                         XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_EOS")
#define AFLOWRC_MPI_NCPUS_MPCDF_EOS                           32 // 32 // MPCDF_EOS_MPICH
#define         MPI_NCPUS_MPCDF_EOS                           XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_EOS")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_EOS                  string("NEGLECT")  // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON
#define         MPI_HYPERTHREADING_MPCDF_EOS                  XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_EOS") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS                      string("~/bin/") // MPCDF_EOS_MPICH
#define         MPI_BINARY_DIR_MPCDF_EOS                      XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_EOS")

#define AFLOWRC_MPI_OPTIONS_MPCDF_DRACO                       string("ulimit -s unlimited ") // MPCDF_DRACO_MPICH  // FIX_DRACO
#define         MPI_OPTIONS_MPCDF_DRACO                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_DRACO")  // FIX_DRACO
#define AFLOWRC_MPI_COMMAND_MPCDF_DRACO                       string("/usr/bin/srun -n") // MPCDF_DRACO_MPICH // FIX_DRACO
#define         MPI_COMMAND_MPCDF_DRACO                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_DRACO") // FIX_DRACO
#define AFLOWRC_MPI_NCPUS_MPCDF_DRACO                         0 // 32 // MPCDF_DRACO_MPICH
#define         MPI_NCPUS_MPCDF_DRACO                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_DRACO")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_DRACO                string("OFF")  // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON
#define         MPI_HYPERTHREADING_MPCDF_DRACO                XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_DRACO") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO                    string("~/bin/") // MPCDF_DRACO_MPICH // FIX_DRACO
#define         MPI_BINARY_DIR_MPCDF_DRACO                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_DRACO") // FIX_DRACO

#define AFLOWRC_MPI_OPTIONS_MPCDF_COBRA                       string("ulimit -s unlimited ") // MPCDF_COBRA_MPICH  // FIX_COBRA
#define         MPI_OPTIONS_MPCDF_COBRA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_COBRA")  // FIX_COBRA
#define AFLOWRC_MPI_COMMAND_MPCDF_COBRA                       string("/usr/bin/srun -n") // MPCDF_COBRA_MPICH // FIX_COBRA
#define         MPI_COMMAND_MPCDF_COBRA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_COBRA") // FIX_COBRA
#define AFLOWRC_MPI_NCPUS_MPCDF_COBRA                         0 // 40 // MPCDF_COBRA_MPICH
#define         MPI_NCPUS_MPCDF_COBRA                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_COBRA")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_COBRA                string("OFF")  // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON
#define         MPI_HYPERTHREADING_MPCDF_COBRA                XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_COBRA") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_COBRA                    string("~/bin/") // MPCDF_COBRA_MPICH // FIX_COBRA
#define         MPI_BINARY_DIR_MPCDF_COBRA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_COBRA") // FIX_COBRA

#define AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA                       string("ulimit -s unlimited ") // MPCDF_HYDRA_MPICH
#define         MPI_OPTIONS_MPCDF_HYDRA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_HYDRA")
#define AFLOWRC_MPI_COMMAND_MPCDF_HYDRA                       string("poe ") // MPCDF_HYDRA_MPICH
#define         MPI_COMMAND_MPCDF_HYDRA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_HYDRA")
#define AFLOWRC_MPI_NCPUS_MPCDF_HYDRA                         0 // 24 // MPCDF_HYDRA_MPICH
#define         MPI_NCPUS_MPCDF_HYDRA                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_HYDRA")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_HYDRA                string("OFF")  // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON
#define         MPI_HYPERTHREADING_MPCDF_HYDRA                XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_HYDRA") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA                    string("~/bin/") // MPCDF_HYDRA_MPICH
#define         MPI_BINARY_DIR_MPCDF_HYDRA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_HYDRA")

//define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    string("module purge") // FULTON_MARYLOU
#define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    string("export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE") // FULTON_MARYLOU
//#define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    string("export OMP_NUM_THREADS=$SLURM_NTASKS") // FULTON_MARYLOU
#define         MPI_OPTIONS_FULTON_MARYLOU                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_FULTON_MARYLOU")
#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    string("srun ") // FULTON_MARYLOU  
//#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    string("mpiexec") // FULTON_MARYLOU  
//#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    string("mpiexec -np") // FULTON_MARYLOU WITH NP
#define         MPI_COMMAND_FULTON_MARYLOU                    XHOST.adefault.getattachedscheme("MPI_COMMAND_FULTON_MARYLOU")
#define AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU                 string("/fslgroup/fslg_datamining/bin/") // FULTON_MARYLOU
#define         MPI_BINARY_DIR_FULTON_MARYLOU                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_FULTON_MARYLOU")

//DX - CMU EULER - START
#define AFLOWRC_MPI_OPTIONS_CMU_EULER                         string("") // CMU EULER
#define         MPI_OPTIONS_CMU_EULER                         XHOST.adefault.getattachedscheme("MPI_OPTIONS_CMU_EULER")
#define AFLOWRC_MPI_COMMAND_CMU_EULER                         string("mpirun -np") // CMU_EULER
#define         MPI_COMMAND_CMU_EULER                         XHOST.adefault.getattachedscheme("MPI_COMMAND_CMU_EULER")
#define AFLOWRC_MPI_BINARY_DIR_CMU_EULER                      string("/home/Tools/bin/") // CMU_EULER
#define         MPI_BINARY_DIR_CMU_EULER                      XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_CMU_EULER")
//DX - CMU EULER - END

#define AFLOWRC_MPI_OPTIONS_MACHINE1                          string("") // future expansions
#define         MPI_OPTIONS_MACHINE1                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE1")
#define AFLOWRC_MPI_COMMAND_MACHINE1                          string("...something ...")  // future expansions
#define         MPI_COMMAND_MACHINE1                          XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE1")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE1                       string("/somewhere/")  // future expansions
#define         MPI_BINARY_DIR_MACHINE1                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE1")

#define AFLOWRC_MPI_OPTIONS_MACHINE2                          string("") // future expansions
#define         MPI_OPTIONS_MACHINE2                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE2")
#define AFLOWRC_MPI_COMMAND_MACHINE2                          string("stub not used")  // future expansions
#define         MPI_COMMAND_MACHINE2                          XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE2")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE2                       string("/home/aflow/bin/")  // future expansions
#define         MPI_BINARY_DIR_MACHINE2                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE2")

#endif // _AFLOW_AFLOWRC_H_

// POCC STUFF
// defaults go in 4 positions; Here with #define AFLOWRC_DEFAULT, in read(), in write_default(), and in print_aflowrc()...
// I coded strings (without spaces), <int>.. you can do <doubles> just like the <int>
// for strings with spaces I need to fix the code. Dont add them now. Then you need to go around the whole code and fix the use of DEFAULTS, possibly also in the READMEs if they are specified.
// STRING     string blablabla=DEFAULT_STRING =>  string blablabla=XHOST.adefault.getattachedscheme("DEFAULT_STRING")
// INT        int blablabla=DEFAULT_INT       =>  int blablabla=XHOST.adefault.getattachedscheme<int>("DEFAULT_INT")
// DOUBLE     double blablabla=DEFAULT_DOUBLE =>  double blablabla=XHOST.adefault.getattachedscheme<double>("DEFAULT_DOUBLE")
// ./aflow --machine to check them out


#ifndef _AFLOW_AFLOWRC_CPP_
#define _AFLOW_AFLOWRC_CPP_

// ***************************************************************************
// aflowrc::load_default
// ***************************************************************************
namespace aflowrc {
  bool load_default(string schema,string schema_default) {
    bool found=FALSE;
    string aus,string_to_add=schema_default;
    vector<string> tokens;
    for(uint i=0;i<XHOST.vaflowrc.size()&&!found;i++) {
      aurostd::string2tokens(XHOST.vaflowrc.at(i),tokens,"=");
      if(tokens.size()>0&&!found) {
        if(aurostd::RemoveWhiteSpaces(tokens.at(0))==schema&&!found) {
          //CO20181226 - it is possible to have '=' inside value: MPI_START_DEFAULT="export OMP_NUM_THREADS=1"
          //treat tokens[0] as special
          tokens.erase(tokens.begin()); //remove key
          found=TRUE; aus=aurostd::RemoveWhiteSpacesFromTheBack(aurostd::joinWDelimiter(tokens,"="));  //CO20180705 - if there are spaces between ="" and //, then the value is set to the spaces (not an empty string!)  //CO20181226 - join again by '='
          aurostd::string2tokens(aus,tokens,"\""); //	    if(tokens.size()>0) cerr << tokens.at(0) << endl;
          if(tokens.size()>0) string_to_add=tokens.at(0);
        }
      }
    }
    // fix ~/ with XHOST.user
    if(aurostd::substring2bool(string_to_add,"~/")) aurostd::StringSubst(string_to_add,"~/",XHOST.home+"/");
    XHOST.adefault.push_attached(schema,string_to_add); // add what is present or the default if not present
    return found;
  }
  template<class utype> 
    bool load_default(string schema,utype schema_default) {
      bool found=XHOST.adefault.args2addattachedscheme(XHOST.vaflowrc,schema,string(schema+"="),""); // add what is present
      if(!found) XHOST.adefault.push_attached(schema,aurostd::utype2string<utype>(schema_default));  // add default if not present
      return found;
    }
}

// ***************************************************************************
// aflowrc::is_available
// ***************************************************************************
namespace aflowrc {
  bool is_available(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    bool aflowrc_local=FALSE;
    bool aflowrc_global=FALSE;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " BEGIN" << endl;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;
    // TESTING LOCAL OR USER BASED
    if(XHOST.aflowrc_filename.empty()) XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    aflowrc_local=aurostd::FileExist(AFLOWRC_FILENAME_LOCAL);
    aflowrc_global=aurostd::FileExist(AFLOWRC_FILENAME_GLOBAL);

    // LOCAL=TRUE && GLOBAL=TRUE => take LOCAL
    if(aflowrc_local && aflowrc_global) {
      if(LDEBUG) oss << __AFLOW_FUNC__ << " LOCAL=TRUE && GLOBAL=TRUE => LOCAL " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
      if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;
      return TRUE;
    }
    // LOCAL=TRUE && GLOBAL=FALSE => take LOCAL
    if(aflowrc_local && !aflowrc_global) {
      if(LDEBUG) oss << __AFLOW_FUNC__ << " LOCAL=TRUE && GLOBAL=FALSE => LOCAL " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL; 
      if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;
      return TRUE;
    }
    // LOCAL=FALSE && GLOBAL=TRUE => take GLOBAL
    if(!aflowrc_local && aflowrc_global) {
      if(LDEBUG) oss << __AFLOW_FUNC__ << " LOCAL=FALSE && GLOBAL=TRUE => GLOBAL " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_GLOBAL;
      if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;
      return TRUE;
    }
    // LOCAL=FALSE && GLOBAL=FALSE => take NOTHING AND REWRITE
    if(!aflowrc_local && !aflowrc_global) {
      if(LDEBUG) oss << __AFLOW_FUNC__ << " LOCAL=FALSE && GLOBAL=FALSE => NOTHING " << endl;
      XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL; // because it is going to write it
      if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;
      if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;
      return FALSE;
    }

    if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;
    return FALSE;
  }
} // namespace aflowrc


// ***************************************************************************
// aflowrc::read
// ***************************************************************************
namespace aflowrc {
  bool read(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    stringstream message; //CO20200404
    if(LDEBUG) oss << __AFLOW_FUNC__ << " BEGIN" << endl;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;
    if(XHOST.aflowrc_filename.empty()) XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;

    if(!aflowrc::is_available(oss,AFLOWRC_VERBOSE)){
      if(!XHOST.vflag_control.flag("WWW")){ //CO20200404 - new web flag
        if(!(aurostd::substring2bool(XHOST.aflowrc_filename,"/mnt/MAIN") || aurostd::substring2bool(XHOST.aflowrc_filename,"/mnt/uMAIN"))){ //CO20200404 - patching for new disk
          //[CO20200404 - OBSOLETE]cout << "WARNING: aflowrc::read: " << XHOST.aflowrc_filename << " not found, loading DEFAULT values" << endl;
          message << XHOST.aflowrc_filename << " not found, loading DEFAULT values";pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,std::cerr,_LOGGER_MESSAGE_);  //CO20200404 - LEAVE std::cerr here, FR needs this for web
        }
      }
    }

    aurostd::file2string(XHOST.aflowrc_filename,XHOST.aflowrc_content);
    // oss << "BEGIN" << endl << XHOST.aflowrc_content << "END" << endl;
    // XHOST.aflowrc_content=aurostd::RemoveComments(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    //   XHOST.aflowrc_content=aurostd::RemoveWhiteSpaces(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    XHOST.aflowrc_content=aurostd::RemoveComments(XHOST.aflowrc_content); // NOW Clean XHOST.aflowrc_content
    // oss << "BEGIN" << endl << XHOST.aflowrc_content << "END" << endl;
    aurostd::string2vectorstring(XHOST.aflowrc_content,XHOST.vaflowrc); // vectorize

    // DEFAULT DEFINITIONS
    aflowrc::load_default("DEFAULT_KZIP_BIN",AFLOWRC_DEFAULT_KZIP_BIN);
    aflowrc::load_default("DEFAULT_KZIP_EXT",AFLOWRC_DEFAULT_KZIP_EXT);
    aflowrc::load_default("DEFAULT_TMPFS_DIRECTORIES",AFLOWRC_DEFAULT_TMPFS_DIRECTORIES);

    //HE20220218 START
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE", AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_AFLUX_SERVER", AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_AFLUX_PATH", AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_RESTAPI_SERVER", AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_RESTAPI_PATH", AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH);
    aflowrc::load_default("DEFAULT_ENTRY_LOADER_FS_PATH", AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH);
    //HE20220218 STOP

    //ME20191001 START
    // AFLOW database files
    aflowrc::load_default("DEFAULT_AFLOW_DB_FILE", AFLOWRC_DEFAULT_AFLOW_DB_FILE);
    aflowrc::load_default("DEFAULT_AFLOW_DB_STATS_FILE", AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE);
    aflowrc::load_default("DEFAULT_AFLOW_DB_DATA_PATH", AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH);
    aflowrc::load_default("DEFAULT_AFLOW_DB_LOCK_FILE", AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE);
    aflowrc::load_default("DEFAULT_AFLOW_DB_STALE_THRESHOLD", AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD);
    //ME20191001 END
    // FILENAMES FOR AFLOW.ORG ANALYSIS
    aflowrc::load_default("DEFAULT_FILE_AFLOWLIB_ENTRY_OUT",AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
    aflowrc::load_default("DEFAULT_FILE_AFLOWLIB_ENTRY_JSON",AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_JSON);
    aflowrc::load_default("DEFAULT_FILE_EDATA_ORIG_OUT",AFLOWRC_DEFAULT_FILE_EDATA_ORIG_OUT);
    aflowrc::load_default("DEFAULT_FILE_EDATA_RELAX_OUT",AFLOWRC_DEFAULT_FILE_EDATA_RELAX_OUT);
    aflowrc::load_default("DEFAULT_FILE_EDATA_BANDS_OUT",AFLOWRC_DEFAULT_FILE_EDATA_BANDS_OUT);
    aflowrc::load_default("DEFAULT_FILE_DATA_ORIG_OUT",AFLOWRC_DEFAULT_FILE_DATA_ORIG_OUT);
    aflowrc::load_default("DEFAULT_FILE_DATA_RELAX_OUT",AFLOWRC_DEFAULT_FILE_DATA_RELAX_OUT);
    aflowrc::load_default("DEFAULT_FILE_DATA_BANDS_OUT",AFLOWRC_DEFAULT_FILE_DATA_BANDS_OUT);
    aflowrc::load_default("DEFAULT_FILE_EDATA_ORIG_JSON",AFLOWRC_DEFAULT_FILE_EDATA_ORIG_JSON);
    aflowrc::load_default("DEFAULT_FILE_EDATA_RELAX_JSON",AFLOWRC_DEFAULT_FILE_EDATA_RELAX_JSON);
    aflowrc::load_default("DEFAULT_FILE_EDATA_BANDS_JSON",AFLOWRC_DEFAULT_FILE_EDATA_BANDS_JSON);
    aflowrc::load_default("DEFAULT_FILE_DATA_ORIG_JSON",AFLOWRC_DEFAULT_FILE_DATA_ORIG_JSON);
    aflowrc::load_default("DEFAULT_FILE_DATA_RELAX_JSON",AFLOWRC_DEFAULT_FILE_DATA_RELAX_JSON);
    aflowrc::load_default("DEFAULT_FILE_DATA_BANDS_JSON",AFLOWRC_DEFAULT_FILE_DATA_BANDS_JSON);
    aflowrc::load_default("DEFAULT_FILE_TIME_OUT",AFLOWRC_DEFAULT_FILE_TIME_OUT);
    aflowrc::load_default("DEFAULT_FILE_SPACEGROUP1_OUT",AFLOWRC_DEFAULT_FILE_SPACEGROUP1_OUT);
    aflowrc::load_default("DEFAULT_FILE_SPACEGROUP2_OUT",AFLOWRC_DEFAULT_FILE_SPACEGROUP2_OUT);
    aflowrc::load_default("DEFAULT_FILE_VOLDISTPARAMS_OUT",AFLOWRC_DEFAULT_FILE_VOLDISTPARAMS_OUT);
    aflowrc::load_default("DEFAULT_FILE_VOLDISTEVOLUTION_OUT",AFLOWRC_DEFAULT_FILE_VOLDISTEVOLUTION_OUT);

    // FILENAMES FOR AFLOW OPERATION
    aflowrc::load_default("DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT",AFLOWRC_DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PRESCRIPT_OUT",AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PRESCRIPT_COMMAND",AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_COMMAND);
    aflowrc::load_default("DEFAULT_AFLOW_POSTSCRIPT_OUT",AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_POSTSCRIPT_COMMAND",AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_COMMAND);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_OUT",AFLOWRC_DEFAULT_AFLOW_PGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_JSON",AFLOWRC_DEFAULT_AFLOW_PGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_XTAL_OUT",AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUP_XTAL_JSON",AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT",AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT); //DX20200129
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON",AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON); //DX20200129
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_OUT",AFLOWRC_DEFAULT_AFLOW_PGROUPK_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_JSON",AFLOWRC_DEFAULT_AFLOW_PGROUPK_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_XTAL_OUT",AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PGROUPK_XTAL_JSON",AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_JSON);  
    aflowrc::load_default("DEFAULT_AFLOW_FGROUP_OUT",AFLOWRC_DEFAULT_AFLOW_FGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FGROUP_JSON",AFLOWRC_DEFAULT_AFLOW_FGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_SGROUP_OUT",AFLOWRC_DEFAULT_AFLOW_SGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_SGROUP_JSON",AFLOWRC_DEFAULT_AFLOW_SGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_AGROUP_OUT",AFLOWRC_DEFAULT_AFLOW_AGROUP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_AGROUP_JSON",AFLOWRC_DEFAULT_AFLOW_AGROUP_JSON);
    aflowrc::load_default("DEFAULT_AFLOW_IATOMS_OUT",AFLOWRC_DEFAULT_AFLOW_IATOMS_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_IATOMS_JSON",AFLOWRC_DEFAULT_AFLOW_IATOMS_JSON);  
    aflowrc::load_default("DEFAULT_AFLOW_ICAGES_OUT",AFLOWRC_DEFAULT_AFLOW_ICAGES_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_SURFACE_OUT",AFLOWRC_DEFAULT_AFLOW_SURFACE_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_QMVASP_OUT",AFLOWRC_DEFAULT_AFLOW_QMVASP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_ERVASP_OUT",AFLOWRC_DEFAULT_AFLOW_ERVASP_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_IMMISCIBILITY_OUT",AFLOWRC_DEFAULT_AFLOW_IMMISCIBILITY_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_MEMORY_OUT",AFLOWRC_DEFAULT_AFLOW_MEMORY_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_INPUT_OUT",AFLOWRC_DEFAULT_AFLOW_FROZSL_INPUT_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_POSCAR_OUT",AFLOWRC_DEFAULT_AFLOW_FROZSL_POSCAR_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_MODES_OUT",AFLOWRC_DEFAULT_AFLOW_FROZSL_MODES_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_FROZSL_EIGEN_OUT",AFLOWRC_DEFAULT_AFLOW_FROZSL_EIGEN_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_END_OUT",AFLOWRC_DEFAULT_AFLOW_END_OUT);
    aflowrc::load_default("DEFAULT_AFLOW_PLASMONICS_FILE",AFLOWRC_DEFAULT_AFLOW_PLASMONICS_FILE);

    // DEFAULT GENERIC MPI
    aflowrc::load_default("MPI_START_DEFAULT",AFLOWRC_MPI_START_DEFAULT); 
    aflowrc::load_default("MPI_STOP_DEFAULT",AFLOWRC_MPI_STOP_DEFAULT); 
    aflowrc::load_default("MPI_COMMAND_DEFAULT",AFLOWRC_MPI_COMMAND_DEFAULT); 
    aflowrc::load_default("MPI_NCPUS_DEFAULT",AFLOWRC_MPI_NCPUS_DEFAULT); 
    aflowrc::load_default("MPI_NCPUS_MAX",AFLOWRC_MPI_NCPUS_MAX); 

    // BINARY VASP
    aflowrc::load_default("DEFAULT_VASP_GAMMA_BIN",AFLOWRC_DEFAULT_VASP_GAMMA_BIN); 
    aflowrc::load_default("DEFAULT_VASP_GAMMA_MPI_BIN",AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN); 
    aflowrc::load_default("DEFAULT_VASP_BIN",AFLOWRC_DEFAULT_VASP_BIN); 
    aflowrc::load_default("DEFAULT_VASP_MPI_BIN",AFLOWRC_DEFAULT_VASP_MPI_BIN); 
    aflowrc::load_default("DEFAULT_VASP5_BIN",AFLOWRC_DEFAULT_VASP5_BIN); 
    aflowrc::load_default("DEFAULT_VASP5_MPI_BIN",AFLOWRC_DEFAULT_VASP5_MPI_BIN); 
    // BINARY AIMS
    aflowrc::load_default("DEFAULT_AIMS_BIN",AFLOWRC_DEFAULT_AIMS_BIN); 

    // POTCARS
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIRECTORIES",AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE",AFLOWRC_DEFAULT_VASP_POTCAR_DATE); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_SUFFIX",AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE_POT_LDA",AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DATE_POT_GGA",AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA); 
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POT_LDA",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_LDA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POT_GGA",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_GGA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POT_PBE",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_PBE);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN);
    aflowrc::load_default("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN",AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN);

    // DEFAULT KPOINTS/DOS
    aflowrc::load_default("DEFAULT_BANDS_GRID",AFLOWRC_DEFAULT_BANDS_GRID); 
    aflowrc::load_default("DEFAULT_BANDS_LATTICE",AFLOWRC_DEFAULT_BANDS_LATTICE); 
    aflowrc::load_default("DEFAULT_KSCHEME",AFLOWRC_DEFAULT_KSCHEME); 
    aflowrc::load_default("DEFAULT_KPPRA",AFLOWRC_DEFAULT_KPPRA); 
    aflowrc::load_default("DEFAULT_STATIC_KSCHEME",AFLOWRC_DEFAULT_STATIC_KSCHEME); 
    aflowrc::load_default("DEFAULT_KPPRA_STATIC",AFLOWRC_DEFAULT_KPPRA_STATIC); 
    aflowrc::load_default("DEFAULT_KPPRA_ICSD",AFLOWRC_DEFAULT_KPPRA_ICSD); 
    aflowrc::load_default("DEFAULT_UNARY_BANDS_GRID",AFLOWRC_DEFAULT_UNARY_BANDS_GRID); 
    aflowrc::load_default("DEFAULT_UNARY_KPPRA",AFLOWRC_DEFAULT_UNARY_KPPRA); 
    aflowrc::load_default("DEFAULT_UNARY_KPPRA_STATIC",AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC); 
    aflowrc::load_default("DEFAULT_PHONONS_KSCHEME",AFLOWRC_DEFAULT_PHONONS_KSCHEME); 
    aflowrc::load_default("DEFAULT_PHONONS_KPPRA",AFLOWRC_DEFAULT_PHONONS_KPPRA); 
    aflowrc::load_default("DEFAULT_DOS_EMIN",AFLOWRC_DEFAULT_DOS_EMIN); 
    aflowrc::load_default("DEFAULT_DOS_EMAX",AFLOWRC_DEFAULT_DOS_EMAX); 
    aflowrc::load_default("DEFAULT_DOS_SCALE",AFLOWRC_DEFAULT_DOS_SCALE); 

    // PRECISION
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_LOW",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_LOW);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_MEDIUM",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_MEDIUM);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_NORMAL",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_NORMAL);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_HIGH",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_HIGH);
    aflowrc::load_default("DEFAULT_VASP_PREC_ENMAX_ACCURATE",AFLOWRC_DEFAULT_VASP_PREC_ENMAX_ACCURATE);
    aflowrc::load_default("DEFAULT_VASP_ENMAX_MINIMUM",AFLOWRC_DEFAULT_VASP_ENMAX_MINIMUM);
    aflowrc::load_default("DEFAULT_VASP_SPIN_REMOVE_CUTOFF",AFLOWRC_DEFAULT_VASP_SPIN_REMOVE_CUTOFF);
    aflowrc::load_default("DEFAULT_VASP_PREC_POTIM",AFLOWRC_DEFAULT_VASP_PREC_POTIM);
    aflowrc::load_default("DEFAULT_VASP_PREC_EDIFFG",AFLOWRC_DEFAULT_VASP_PREC_EDIFFG);

    // OPTIONS
    aflowrc::load_default("DEFAULT_VASP_OUT",AFLOWRC_DEFAULT_VASP_OUT);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_INCAR",AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_POSCAR",AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_POTCAR",AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR);
    aflowrc::load_default("DEFAULT_VASP_EXTERNAL_KPOINTS",AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS);
    aflowrc::load_default("DEFAULT_AIMS_EXTERNAL_CONTROL",AFLOWRC_DEFAULT_AIMS_EXTERNAL_CONTROL);
    aflowrc::load_default("DEFAULT_AIMS_EXTERNAL_GEOM",AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM);
    aflowrc::load_default("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE",AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SIGMA",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_NELM",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM); //CO20200624
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_NELM_STATIC",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM_STATIC); //CO20200624
    aflowrc::load_default("MAX_VASP_NELM",AFLOWRC_MAX_VASP_NELM); //CO20200624
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SYM",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_BADER",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_BADER_STATIC",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER_STATIC);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_ELF",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_WAVECAR",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_CHGCAR",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR);
    aflowrc::load_default("DEFAULT_VASP_FORCE_OPTION_LSCOUPLING",AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING);

    // AFLOW_LIBRARY AFLOW_PROJECT
    aflowrc::load_default("DEFAULT_AFLOW_LIBRARY_DIRECTORIES",AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES);
    aflowrc::load_default("DEFAULT_AFLOW_PROJECTS_DIRECTORIES",AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES);
    aflowrc::load_default("DEFAULT_AFLOWDATA_WEB_DIRECTORY",AFLOWRC_DEFAULT_AFLOWDATA_WEB_DIRECTORY); //CO+ME20200731

    // DEFAULT PLATON/FINDSYM
    aflowrc::load_default("DEFAULT_PLATON_P_EQUAL",AFLOWRC_DEFAULT_PLATON_P_EQUAL);
    aflowrc::load_default("DEFAULT_PLATON_P_EXACT",AFLOWRC_DEFAULT_PLATON_P_EXACT);
    aflowrc::load_default("DEFAULT_PLATON_P_ANG",AFLOWRC_DEFAULT_PLATON_P_ANG);
    aflowrc::load_default("DEFAULT_PLATON_P_D1",AFLOWRC_DEFAULT_PLATON_P_D1);
    aflowrc::load_default("DEFAULT_PLATON_P_D2",AFLOWRC_DEFAULT_PLATON_P_D2);
    aflowrc::load_default("DEFAULT_PLATON_P_D3",AFLOWRC_DEFAULT_PLATON_P_D3);
    aflowrc::load_default("DEFAULT_FINDSYM_TOL",AFLOWRC_DEFAULT_FINDSYM_TOL);

    // DEFAULT GNUPLOT
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT);
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_BOLD",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD); 
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_BOLD",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT",AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT); 
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT_BOLD",AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD); 
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_ITALICS); 
    aflowrc::load_default("DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS",AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS); 

    // DEFAULT CHULL
    aflowrc::load_default("DEFAULT_CHULL_ALLOWED_DFT_TYPES",AFLOWRC_DEFAULT_CHULL_ALLOWED_DFT_TYPES); 
    aflowrc::load_default("DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES",AFLOWRC_DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES); 
    aflowrc::load_default("DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES",AFLOWRC_DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES); 
    aflowrc::load_default("DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS",AFLOWRC_DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS); 
    aflowrc::load_default("DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES",AFLOWRC_DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES); 
    aflowrc::load_default("DEFAULT_CHULL_OUTLIER_MULTIPLIER",AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER); 
    aflowrc::load_default("DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED",AFLOWRC_DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_BANNER",AFLOWRC_DEFAULT_CHULL_LATEX_BANNER); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN",AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_STOICH_HEADER",AFLOWRC_DEFAULT_CHULL_LATEX_STOICH_HEADER); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_UNARIES",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNARIES); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_OFF_HULL",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_UNSTABLE",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FILTER_SCHEME",AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FILTER_VALUE",AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COLOR_BAR",AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_HEAT_MAP",AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COLOR_GRADIENT",AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_COLOR_MAP",AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR",AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_REVERSE_AXIS",AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW",AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_LINKS",AFLOWRC_DEFAULT_CHULL_LATEX_LINKS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_LABEL_NAME",AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_META_LABELS",AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_LABELS_OFF_HULL",AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION",AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_HELVETICA_FONT",AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_FONT_SIZE",AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_ROTATE_LABELS",AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS); 
    aflowrc::load_default("DEFAULT_CHULL_LATEX_BOLD_LABELS",AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS); 
    aflowrc::load_default("DEFAULT_CHULL_PNG_RESOLUTION",AFLOWRC_DEFAULT_CHULL_PNG_RESOLUTION); 

    // DEFAULT GFA  //CO20190628
    aflowrc::load_default("DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF",AFLOWRC_DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF); //CO20190628

    // DEFAULT ARUN
    aflowrc::load_default("ARUN_DIRECTORY_PREFIX",AFLOWRC_ARUN_DIRECTORY_PREFIX);

    // DEFAULT POCC
    aflowrc::load_default("DEFAULT_POCC_STRUCTURE_GENERATION_ALGO",AFLOWRC_DEFAULT_POCC_STRUCTURE_GENERATION_ALGO);
    aflowrc::load_default("DEFAULT_POCC_TEMPERATURE_STRING",AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING);
    aflowrc::load_default("DEFAULT_POCC_EXCLUDE_UNSTABLE",AFLOWRC_DEFAULT_POCC_EXCLUDE_UNSTABLE);  //ME20210927
    aflowrc::load_default("DEFAULT_POCC_SITE_TOL",AFLOWRC_DEFAULT_POCC_SITE_TOL);
    aflowrc::load_default("DEFAULT_POCC_STOICH_TOL",AFLOWRC_DEFAULT_POCC_STOICH_TOL);
    aflowrc::load_default("DEFAULT_UFF_BONDING_DISTANCE",AFLOWRC_DEFAULT_UFF_BONDING_DISTANCE);
    aflowrc::load_default("DEFAULT_UFF_ENERGY_TOLERANCE",AFLOWRC_DEFAULT_UFF_ENERGY_TOLERANCE);
    aflowrc::load_default("DEFAULT_UFF_CLUSTER_RADIUS",AFLOWRC_DEFAULT_UFF_CLUSTER_RADIUS);
    aflowrc::load_default("DEFAULT_POCC_RDF_RMAX",AFLOWRC_DEFAULT_POCC_RDF_RMAX);
    aflowrc::load_default("DEFAULT_POCC_RDF_NBINS",AFLOWRC_DEFAULT_POCC_RDF_NBINS);
    aflowrc::load_default("DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON",AFLOWRC_DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON);
    aflowrc::load_default("DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS",AFLOWRC_DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS);
    aflowrc::load_default("POCC_FILE_PREFIX",AFLOWRC_POCC_FILE_PREFIX);
    aflowrc::load_default("POCC_OUT_FILE",AFLOWRC_POCC_OUT_FILE);
    aflowrc::load_default("POCC_APL_OUT_FILE",AFLOWRC_POCC_APL_OUT_FILE);  //ME20210927
    aflowrc::load_default("POCC_ALL_SUPERCELLS_FILE",AFLOWRC_POCC_ALL_SUPERCELLS_FILE);
    aflowrc::load_default("POCC_UNIQUE_SUPERCELLS_FILE",AFLOWRC_POCC_UNIQUE_SUPERCELLS_FILE);
    aflowrc::load_default("POCC_ALL_HNF_MATRICES_FILE",AFLOWRC_POCC_ALL_HNF_MATRICES_FILE);
    aflowrc::load_default("POCC_ALL_SITE_CONFIGURATIONS_FILE",AFLOWRC_POCC_ALL_SITE_CONFIGURATIONS_FILE);
    aflowrc::load_default("POCC_DOSCAR_FILE",AFLOWRC_POCC_DOSCAR_FILE);
    aflowrc::load_default("POCC_PHDOSCAR_FILE",AFLOWRC_POCC_PHDOSCAR_FILE);  //ME20210927
    aflowrc::load_default("POCC_ANIONS_LIST",AFLOWRC_POCC_ANIONS_LIST);

    // DEFAULT APL
    //// DEFAULT APL SUPERCELL
    aflowrc::load_default("DEFAULT_APL_PREC",AFLOWRC_DEFAULT_APL_PREC);
    aflowrc::load_default("DEFAULT_APL_ENGINE",AFLOWRC_DEFAULT_APL_ENGINE);
    aflowrc::load_default("DEFAULT_APL_HIBERNATE",AFLOWRC_DEFAULT_APL_HIBERNATE);
    aflowrc::load_default("DEFAULT_APL_MINSHELL",AFLOWRC_DEFAULT_APL_MINSHELL);
    aflowrc::load_default("DEFAULT_APL_MINATOMS",AFLOWRC_DEFAULT_APL_MINATOMS);
    aflowrc::load_default("DEFAULT_APL_POLAR",AFLOWRC_DEFAULT_APL_POLAR);
    aflowrc::load_default("DEFAULT_APL_DMAG",AFLOWRC_DEFAULT_APL_DMAG);
    aflowrc::load_default("DEFAULT_APL_DXYZONLY",AFLOWRC_DEFAULT_APL_DXYZONLY);
    aflowrc::load_default("DEFAULT_APL_DSYMMETRIZE",AFLOWRC_DEFAULT_APL_DSYMMETRIZE); //CO20181226
    aflowrc::load_default("DEFAULT_APL_DINEQUIV_ONLY",AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY); //CO20181226
    aflowrc::load_default("DEFAULT_APL_DPM",AFLOWRC_DEFAULT_APL_DPM);
    aflowrc::load_default("DEFAULT_APL_RELAX",AFLOWRC_DEFAULT_APL_RELAX);
    aflowrc::load_default("DEFAULT_APL_RELAX_COMMENSURATE",AFLOWRC_DEFAULT_APL_RELAX_COMMENSURATE);  //ME20200427
    aflowrc::load_default("DEFAULT_APL_ZEROSTATE",AFLOWRC_DEFAULT_APL_ZEROSTATE);
    aflowrc::load_default("DEFAULT_APL_ZEROSTATE_CHGCAR",AFLOWRC_DEFAULT_APL_ZEROSTATE_CHGCAR);  //ME20191029
    aflowrc::load_default("DEFAULT_APL_USE_LEPSILON",AFLOWRC_DEFAULT_APL_USE_LEPSILON);

    //// DEFAULT APL PHONON PROPERTIES
    aflowrc::load_default("DEFAULT_APL_FREQFORMAT",AFLOWRC_DEFAULT_APL_FREQFORMAT);
    aflowrc::load_default("DEFAULT_APL_DC",AFLOWRC_DEFAULT_APL_DC);
    aflowrc::load_default("DEFAULT_APL_DCPATH",AFLOWRC_DEFAULT_APL_DCPATH);
    aflowrc::load_default("DEFAULT_APL_DCPOINTS",AFLOWRC_DEFAULT_APL_DCPOINTS);
    aflowrc::load_default("DEFAULT_APL_DOS",AFLOWRC_DEFAULT_APL_DOS);
    aflowrc::load_default("DEFAULT_APL_DOSMETHOD",AFLOWRC_DEFAULT_APL_DOSMETHOD);
    aflowrc::load_default("DEFAULT_APL_DOSMESH",AFLOWRC_DEFAULT_APL_DOSMESH);
    aflowrc::load_default("DEFAULT_APL_DOSPOINTS",AFLOWRC_DEFAULT_APL_DOSPOINTS);
    aflowrc::load_default("DEFAULT_APL_DOSSMEAR",AFLOWRC_DEFAULT_APL_DOSSMEAR);
    aflowrc::load_default("DEFAULT_APL_DOS_PROJECT",AFLOWRC_DEFAULT_APL_DOS_PROJECT);  //ME20200213
    aflowrc::load_default("DEFAULT_APL_TP",AFLOWRC_DEFAULT_APL_TP);
    aflowrc::load_default("DEFAULT_APL_DISPLACEMENTS",AFLOWRC_DEFAULT_APL_DISPLACEMENTS);  //ME20200421
    aflowrc::load_default("DEFAULT_APL_TPT",AFLOWRC_DEFAULT_APL_TPT);
    aflowrc::load_default("DEFAULT_APL_GVEL",AFLOWRC_DEFAULT_APL_GVEL);  //ME20200517

    //// DEFAULT APL FILES
    aflowrc::load_default("DEFAULT_APL_FILE_PREFIX",AFLOWRC_DEFAULT_APL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_APL_OUT_FILE",AFLOWRC_DEFAULT_APL_OUT_FILE);  //ME20210927
    aflowrc::load_default("DEFAULT_APL_PDIS_FILE",AFLOWRC_DEFAULT_APL_PDIS_FILE);
    aflowrc::load_default("DEFAULT_APL_PDOS_FILE",AFLOWRC_DEFAULT_APL_PDOS_FILE);
    aflowrc::load_default("DEFAULT_APL_THERMO_FILE",AFLOWRC_DEFAULT_APL_THERMO_FILE);
    aflowrc::load_default("DEFAULT_APL_THERMO_JSON",AFLOWRC_DEFAULT_APL_THERMO_JSON);  //ME20211019
    aflowrc::load_default("DEFAULT_APL_DYNMAT_FILE",AFLOWRC_DEFAULT_APL_DYNMAT_FILE);
    aflowrc::load_default("DEFAULT_APL_HARMIFC_FILE",AFLOWRC_DEFAULT_APL_HARMIFC_FILE);
    aflowrc::load_default("DEFAULT_APL_POLAR_FILE",AFLOWRC_DEFAULT_APL_POLAR_FILE);  //ME20200415
    aflowrc::load_default("DEFAULT_APL_HSKPTS_FILE",AFLOWRC_DEFAULT_APL_HSKPTS_FILE);
    aflowrc::load_default("DEFAULT_APL_MSQRDISP_FILE", AFLOWRC_DEFAULT_APL_MSQRDISP_FILE);  //ME20200329
    aflowrc::load_default("DEFAULT_APL_GVEL_FILE", AFLOWRC_DEFAULT_APL_GVEL_FILE);  // ME20200517
    //ME20190614 BEGIN
    aflowrc::load_default("DEFAULT_APL_PHDOSCAR_FILE",AFLOWRC_DEFAULT_APL_PHDOSCAR_FILE);
    aflowrc::load_default("DEFAULT_APL_PHPOSCAR_FILE",AFLOWRC_DEFAULT_APL_PHPOSCAR_FILE);
    aflowrc::load_default("DEFAULT_APL_PHKPOINTS_FILE",AFLOWRC_DEFAULT_APL_PHKPOINTS_FILE);
    aflowrc::load_default("DEFAULT_APL_PHEIGENVAL_FILE",AFLOWRC_DEFAULT_APL_PHEIGENVAL_FILE);
    //ME20190614 END
    aflowrc::load_default("DEFAULT_APL_STATE_FILE",AFLOWRC_DEFAULT_APL_STATE_FILE);  //ME20200224
    //ME20200329 BEGIN
    aflowrc::load_default("DEFAULT_APL_ADISP_SCENE_FORMAT",AFLOWRC_DEFAULT_APL_ADISP_SCENE_FORMAT);
    aflowrc::load_default("DEFAULT_APL_ADISP_AMPLITUDE",AFLOWRC_DEFAULT_APL_ADISP_AMPLITUDE);
    aflowrc::load_default("DEFAULT_APL_ADISP_NSTEPS",AFLOWRC_DEFAULT_APL_ADISP_NSTEPS);
    aflowrc::load_default("DEFAULT_APL_ADISP_NPERIODS",AFLOWRC_DEFAULT_APL_ADISP_NPERIODS);
    //ME20200329 END

    // DEFAULT QHA
    //// DEFAULT QHA VALUES
    aflowrc::load_default("DEFAULT_QHA_MODE", AFLOWRC_DEFAULT_QHA_MODE);
    aflowrc::load_default("DEFAULT_QHA_EOS", AFLOWRC_DEFAULT_QHA_EOS);
    aflowrc::load_default("DEFAULT_QHA_EOS_DISTORTION_RANGE", AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE);
    aflowrc::load_default("DEFAULT_QHA_EOS_MODEL", AFLOWRC_DEFAULT_QHA_EOS_MODEL);//AS20200818
    aflowrc::load_default("DEFAULT_QHA_GP_DISTORTION", AFLOWRC_DEFAULT_QHA_GP_DISTORTION);
    aflowrc::load_default("DEFAULT_QHA_TAYLOR_EXPANSION_ORDER", AFLOWRC_DEFAULT_QHA_TAYLOR_EXPANSION_ORDER);//AS20200602
    aflowrc::load_default("DEFAULT_QHA_INCLUDE_ELEC_CONTRIB", AFLOWRC_DEFAULT_QHA_INCLUDE_ELEC_CONTRIB);
    aflowrc::load_default("DEFAULT_QHA_SOMMERFELD_EXPANSION", AFLOWRC_DEFAULT_QHA_SOMMERFELD_EXPANSION);//AS20200528
    aflowrc::load_default("DEFAULT_QHA_PDIS_T", AFLOWRC_DEFAULT_QHA_PDIS_T);
    //AS20200508 BEGIN
    aflowrc::load_default("DEFAULT_QHA_GP_FINITE_DIFF", AFLOWRC_DEFAULT_QHA_GP_FINITE_DIFF);
    aflowrc::load_default("DEFAULT_QHA_IGNORE_IMAGINARY", AFLOWRC_DEFAULT_QHA_IGNORE_IMAGINARY);
    aflowrc::load_default("DEFAULT_QHA_RELAX_IONS_CELL", AFLOWRC_DEFAULT_QHA_RELAX_IONS_CELL);//AS20201123
    //// DEFAULT QHA FILES
    aflowrc::load_default("DEFAULT_QHA_FILE_PREFIX", AFLOWRC_DEFAULT_QHA_FILE_PREFIX);
    //AS20200709 BEGIN
    aflowrc::load_default("DEFAULT_QHA3P_FILE_PREFIX", AFLOWRC_DEFAULT_QHA3P_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_QHANP_FILE_PREFIX", AFLOWRC_DEFAULT_QHANP_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_SCQHA_FILE_PREFIX", AFLOWRC_DEFAULT_SCQHA_FILE_PREFIX);
    //AS20200709 END
    aflowrc::load_default("DEFAULT_QHA_GP_PATH_FILE", AFLOWRC_DEFAULT_QHA_GP_PATH_FILE);
    aflowrc::load_default("DEFAULT_QHA_GP_MESH_FILE", AFLOWRC_DEFAULT_QHA_GP_MESH_FILE);
    aflowrc::load_default("DEFAULT_QHA_GP_AVG_FILE", AFLOWRC_DEFAULT_QHA_GP_AVG_FILE);
    aflowrc::load_default("DEFAULT_QHA_THERMO_FILE", AFLOWRC_DEFAULT_QHA_THERMO_FILE);
    aflowrc::load_default("DEFAULT_QHA_FREQS_FILE", AFLOWRC_DEFAULT_QHA_FREQS_FILE);
    aflowrc::load_default("DEFAULT_QHA_FVT_FILE", AFLOWRC_DEFAULT_QHA_FVT_FILE);
    //AS20200508 END
    aflowrc::load_default("DEFAULT_QHA_COEFF_FILE", AFLOWRC_DEFAULT_QHA_COEFF_FILE);//AS20210517
    aflowrc::load_default("DEFAULT_QHA_IMAG_FILE", AFLOWRC_DEFAULT_QHA_IMAG_FILE);//AS20210517
    aflowrc::load_default("DEFAULT_QHA_PDIS_FILE", AFLOWRC_DEFAULT_QHA_PDIS_FILE);//AS20201022
    aflowrc::load_default("DEFAULT_QHA_PDOS_FILE", AFLOWRC_DEFAULT_QHA_PDOS_FILE);//AS20201201
    aflowrc::load_default("DEFAULT_QHA_KPOINTS_FILE", AFLOWRC_DEFAULT_QHA_KPOINTS_FILE);//AS20201112
    //AS20210914 BEGIN
    aflowrc::load_default("DEFAULT_POCC_QHA_THERMO_FILE", AFLOWRC_DEFAULT_POCC_QHA_THERMO_FILE);
    aflowrc::load_default("DEFAULT_POCC_QHA_AVGTHERMO_FILE", AFLOWRC_DEFAULT_POCC_QHA_AVGTHERMO_FILE);
    //AS20210914 END

    // DEFAULT AAPL
    //// DEFAULT AAPL VALUES
    aflowrc::load_default("DEFAULT_AAPL_BTE",AFLOWRC_DEFAULT_AAPL_BTE);
    //[ME20181226]aflowrc::load_default("DEFAULT_AAPL_BZMETHOD",AFLOWRC_DEFAULT_AAPL_BZMETHOD);
    aflowrc::load_default("DEFAULT_AAPL_FOURTH_ORDER",AFLOWRC_DEFAULT_AAPL_FOURTH_ORDER);
    aflowrc::load_default("DEFAULT_AAPL_CUT_RAD",AFLOWRC_DEFAULT_AAPL_CUT_RAD);
    aflowrc::load_default("DEFAULT_AAPL_CUT_SHELL",AFLOWRC_DEFAULT_AAPL_CUT_SHELL);
    aflowrc::load_default("DEFAULT_AAPL_THERMALGRID",AFLOWRC_DEFAULT_AAPL_THERMALGRID);
    aflowrc::load_default("DEFAULT_AAPL_TCT",AFLOWRC_DEFAULT_AAPL_TCT);
    aflowrc::load_default("DEFAULT_AAPL_SUMRULE",AFLOWRC_DEFAULT_AAPL_SUMRULE);
    aflowrc::load_default("DEFAULT_AAPL_SUMRULE_MAX_ITER",AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER);
    aflowrc::load_default("DEFAULT_AAPL_MIXING_COEFFICIENT",AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT);
    aflowrc::load_default("DEFAULT_AAPL_ISOTOPE",AFLOWRC_DEFAULT_AAPL_ISOTOPE);
    aflowrc::load_default("DEFAULT_AAPL_BOUNDARY",AFLOWRC_DEFAULT_AAPL_BOUNDARY);
    aflowrc::load_default("DEFAULT_AAPL_CUMULATIVEK",AFLOWRC_DEFAULT_AAPL_CUMULATIVEK);
    aflowrc::load_default("DEFAULT_AAPL_NANO_SIZE",AFLOWRC_DEFAULT_AAPL_NANO_SIZE);
    //// DEFAULT AAPL FILES
    aflowrc::load_default("DEFAULT_AAPL_FILE_PREFIX",AFLOWRC_DEFAULT_AAPL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_AAPL_IRRQPTS_FILE",AFLOWRC_DEFAULT_AAPL_IRRQPTS_FILE);
    aflowrc::load_default("DEFAULT_AAPL_GVEL_FILE",AFLOWRC_DEFAULT_AAPL_GVEL_FILE);
    aflowrc::load_default("DEFAULT_AAPL_PS_FILE",AFLOWRC_DEFAULT_AAPL_PS_FILE);  //ME20191104
    aflowrc::load_default("DEFAULT_AAPL_GRUENEISEN_FILE",AFLOWRC_DEFAULT_AAPL_GRUENEISEN_FILE);  //ME20191104
    aflowrc::load_default("DEFAULT_AAPL_RATES_FILE",AFLOWRC_DEFAULT_AAPL_RATES_FILE);
    aflowrc::load_default("DEFAULT_AAPL_RATES_3RD_FILE",AFLOWRC_DEFAULT_AAPL_RATES_3RD_FILE);
    aflowrc::load_default("DEFAULT_AAPL_RATES_4TH_FILE",AFLOWRC_DEFAULT_AAPL_RATES_4TH_FILE);
    aflowrc::load_default("DEFAULT_AAPL_ISOTOPE_FILE",AFLOWRC_DEFAULT_AAPL_ISOTOPE_FILE);
    aflowrc::load_default("DEFAULT_AAPL_BOUNDARY_FILE",AFLOWRC_DEFAULT_AAPL_BOUNDARY_FILE);
    aflowrc::load_default("DEFAULT_AAPL_TCOND_FILE",AFLOWRC_DEFAULT_AAPL_TCOND_FILE);
    //aflowrc::load_default("DEFAULT_AAPL_TCOND_PLOT_FILE",AFLOWRC_DEFAULT_AAPL_TCOND_PLOT_FILE);  OBSOLETE ME20191104

    // DEFAULT AEL
    //// DEFAULT AEL STRAIN CALCS
    aflowrc::load_default("DEFAULT_AEL_STRAIN_SYMMETRY",AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY);
    aflowrc::load_default("DEFAULT_AEL_NNORMAL_STRAINS",AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS);
    aflowrc::load_default("DEFAULT_AEL_NSHEAR_STRAINS",AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS);
    aflowrc::load_default("DEFAULT_AEL_NORMAL_STRAIN_STEP",AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP);
    aflowrc::load_default("DEFAULT_AEL_SHEAR_STRAIN_STEP",AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP);
    aflowrc::load_default("DEFAULT_AEL_ORIGIN_STRAIN_CALC",AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC);
    aflowrc::load_default("DEFAULT_AEL_ORIGIN_STRAIN_FIT",AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT);
    aflowrc::load_default("DEFAULT_AEL_RELAXED_STRUCT_FIT",AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT);
    aflowrc::load_default("DEFAULT_AEL_NEG_STRAINS",AFLOWRC_DEFAULT_AEL_NEG_STRAINS); 
    aflowrc::load_default("DEFAULT_AEL_NIND_STRAIN_DIRS",AFLOWRC_DEFAULT_AEL_VASPSYM); 
    aflowrc::load_default("DEFAULT_AEL_VASPSYM",AFLOWRC_DEFAULT_AEL_VASPSYM);
    aflowrc::load_default("DEFAULT_AEL_PRECACC_ALGONORM",AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM);
    aflowrc::load_default("DEFAULT_AEL_VASPRUNXML_STRESS",AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS);
    aflowrc::load_default("DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS",AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS);
    aflowrc::load_default("DEFAULT_AEL_SKIP_ARUNS_MAX",AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX);

    //// DEFAULT AEL CHECKS AND PROCESSING
    aflowrc::load_default("DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY",AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY);
    aflowrc::load_default("DEFAULT_AEL_SYMMETRIZE",AFLOWRC_DEFAULT_AEL_SYMMETRIZE);

    //// DEFAULT AEL OUTPUT FILES
    aflowrc::load_default("DEFAULT_AEL_FILE_PREFIX",AFLOWRC_DEFAULT_AEL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_AEL_WRITE_FULL_RESULTS",AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS);
    aflowrc::load_default("DEFAULT_AEL_DIRNAME_ARUN",AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN);

    // DEFAULT AGL
    //// DEFAULT AGL STRAIN CALCS
    aflowrc::load_default("DEFAULT_AGL_AEL_POISSON_RATIO",AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO);
    aflowrc::load_default("DEFAULT_AGL_NSTRUCTURES",AFLOWRC_DEFAULT_AGL_NSTRUCTURES);
    aflowrc::load_default("DEFAULT_AGL_STRAIN_STEP",AFLOWRC_DEFAULT_AGL_STRAIN_STEP);
    aflowrc::load_default("DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS",AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS);
    aflowrc::load_default("DEFAULT_AGL_SKIP_ARUNS_MAX",AFLOWRC_DEFAULT_AGL_SKIP_ARUNS_MAX);

    //// DEFAULT AGL CHECKS AND PROCESSING
    aflowrc::load_default("DEFAULT_AGL_NTEMPERATURE",AFLOWRC_DEFAULT_AGL_NTEMPERATURE);
    aflowrc::load_default("DEFAULT_AGL_STEMPERATURE",AFLOWRC_DEFAULT_AGL_STEMPERATURE);
    aflowrc::load_default("DEFAULT_AGL_NPRESSURE",AFLOWRC_DEFAULT_AGL_NPRESSURE);
    aflowrc::load_default("DEFAULT_AGL_SPRESSURE",AFLOWRC_DEFAULT_AGL_SPRESSURE);    
    aflowrc::load_default("DEFAULT_AGL_POISSON_RATIO",AFLOWRC_DEFAULT_AGL_POISSON_RATIO);
    aflowrc::load_default("DEFAULT_AGL_IEOS",AFLOWRC_DEFAULT_AGL_IEOS);
    aflowrc::load_default("DEFAULT_AGL_IDEBYE",AFLOWRC_DEFAULT_AGL_IDEBYE);
    aflowrc::load_default("DEFAULT_AGL_FIT_TYPE",AFLOWRC_DEFAULT_AGL_FIT_TYPE);
    aflowrc::load_default("DEFAULT_AGL_CHECK_EV_CONCAVITY",AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY);
    aflowrc::load_default("DEFAULT_AGL_CHECK_EV_MIN",AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN);    
    aflowrc::load_default("DEFAULT_AGL_HUGONIOT_CALC",AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC);
    aflowrc::load_default("DEFAULT_AGL_HUGONIOT_EXTRAPOLATE",AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE);
    aflowrc::load_default("DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE",AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE);

    //// DEFAULT AGL OUTPUT FILES
    aflowrc::load_default("DEFAULT_AGL_FILE_PREFIX",AFLOWRC_DEFAULT_AGL_FILE_PREFIX);
    aflowrc::load_default("DEFAULT_AGL_WRITE_FULL_RESULTS",AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS);
    aflowrc::load_default("DEFAULT_AGL_DIRNAME_ARUN",AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN);
    aflowrc::load_default("DEFAULT_AGL_WRITE_GIBBS_INPUT",AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT);
    aflowrc::load_default("DEFAULT_AGL_PLOT_RESULTS",AFLOWRC_DEFAULT_AGL_PLOT_RESULTS);    

    // DEFAULT QCA
    aflowrc::load_default("DEFAULT_QCA_MIN_SLEEP_SECONDS",AFLOWRC_DEFAULT_QCA_MIN_SLEEP_SECONDS);
    aflowrc::load_default("DEFAULT_QCA_MAX_NUM_ATOMS",AFLOWRC_DEFAULT_QCA_MAX_NUM_ATOMS);
    aflowrc::load_default("DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS",AFLOWRC_DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS);
    aflowrc::load_default("DEFAULT_QCA_CV_CUTOFF",AFLOWRC_DEFAULT_QCA_CV_CUTOFF);
    aflowrc::load_default("DEFAULT_QCA_CONC_NPTS",AFLOWRC_DEFAULT_QCA_CONC_NPTS);
    aflowrc::load_default("DEFAULT_QCA_TEMP_NPTS",AFLOWRC_DEFAULT_QCA_TEMP_NPTS);
    aflowrc::load_default("DEFAULT_QCA_TEMP_MIN",AFLOWRC_DEFAULT_QCA_TEMP_MIN);
    aflowrc::load_default("DEFAULT_QCA_TEMP_MAX",AFLOWRC_DEFAULT_QCA_TEMP_MAX);
    aflowrc::load_default("DEFAULT_QCA_TEMP_MIN_LIMIT",AFLOWRC_DEFAULT_QCA_TEMP_MIN_LIMIT);
    aflowrc::load_default("DEFAULT_QCA_PRINT",AFLOWRC_DEFAULT_QCA_PRINT);

    //RF20200413 START
    // DEFAULT CCE
    aflowrc::load_default("DEFAULT_CCE_OX_METHOD",AFLOWRC_DEFAULT_CCE_OX_METHOD);
    aflowrc::load_default("DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION",AFLOWRC_DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION);
    aflowrc::load_default("DEFAULT_CCE_OX_TOL",AFLOWRC_DEFAULT_CCE_OX_TOL);
    aflowrc::load_default("DEFAULT_CCE_PEROX_CUTOFF",AFLOWRC_DEFAULT_CCE_PEROX_CUTOFF);
    aflowrc::load_default("DEFAULT_CCE_SUPEROX_CUTOFF",AFLOWRC_DEFAULT_CCE_SUPEROX_CUTOFF);
    aflowrc::load_default("DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF",AFLOWRC_DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF);
    aflowrc::load_default("DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF",AFLOWRC_DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF);
    //RF20200413 END

    // DEFAULT XTALFINDER
    aflowrc::load_default("DEFAULT_XTALFINDER_MISFIT_MATCH",AFLOWRC_DEFAULT_XTALFINDER_MISFIT_MATCH); //DX20201118
    aflowrc::load_default("DEFAULT_XTALFINDER_MISFIT_FAMILY",AFLOWRC_DEFAULT_XTALFINDER_MISFIT_FAMILY); //DX20201118
    aflowrc::load_default("DEFAULT_XTALFINDER_SUPERCELL_METHOD",AFLOWRC_DEFAULT_XTALFINDER_SUPERCELL_METHOD); //DX20201223
    aflowrc::load_default("DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING",AFLOWRC_DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING); //DX20200709
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_MATERIAL",AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL); //DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_STRUCTURE",AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE); //DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_DUPLICATE",AFLOWRC_DEFAULT_XTALFINDER_FILE_DUPLICATE); //DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE",AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE); //DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE",AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE); //DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE",AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE); //DX20201228
    aflowrc::load_default("DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE",AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE); //DX20201228

    //DX20200720 - START
    // DEFAULT ANRL
    aflowrc::load_default("DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL",AFLOWRC_DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL);
    //DX20200720 - END

    // DEFAULT CORE
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_BEEP",AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP);
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_HALT",AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT);
    aflowrc::load_default("AFLOW_CORE_TEMPERATURE_REFRESH",AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH);

    // VASP MACHINE SETTINGS
    aflowrc::load_default("SECONDS_SLEEP_VASP_COMPLETION",AFLOWRC_SECONDS_SLEEP_VASP_COMPLETION); //CO20201111
    aflowrc::load_default("SECONDS_SLEEP_VASP_MONITOR",AFLOWRC_SECONDS_SLEEP_VASP_MONITOR); //CO20201111
    aflowrc::load_default("SECONDS_STALE_OUTCAR",AFLOWRC_SECONDS_STALE_OUTCAR); //CO20201111
    aflowrc::load_default("BYTES_MAX_VASP_OUT",AFLOWRC_BYTES_MAX_VASP_OUT); //CO20201111
    aflowrc::load_default("MEMORY_MAX_USAGE_RAM",AFLOWRC_MEMORY_MAX_USAGE_RAM); //CO20201111
    aflowrc::load_default("MEMORY_MAX_USAGE_SWAP",AFLOWRC_MEMORY_MAX_USAGE_SWAP); //CO20201111
    aflowrc::load_default("FILE_VASP_MONITOR",AFLOWRC_FILE_VASP_MONITOR); //CO20201111
    aflowrc::load_default("INTEL_COMPILER_PATHS",AFLOWRC_INTEL_COMPILER_PATHS); //CO20201111

    // DEFAULT MACHINE DEPENDENT MPI
    aflowrc::load_default("MPI_OPTIONS_DUKE_BETA_MPICH",AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH); 
    aflowrc::load_default("MPI_COMMAND_DUKE_BETA_MPICH",AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_BETA_MPICH",AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_BETA_OPENMPI",AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI); 
    aflowrc::load_default("MPI_COMMAND_DUKE_BETA_OPENMPI",AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_BETA_OPENMPI",AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_MATERIALS",AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS); 
    aflowrc::load_default("MPI_COMMAND_DUKE_MATERIALS",AFLOWRC_MPI_COMMAND_DUKE_MATERIALS); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_MATERIALS",AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_AFLOWLIB",AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB); 
    aflowrc::load_default("MPI_COMMAND_DUKE_AFLOWLIB",AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_AFLOWLIB",AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_HABANA",AFLOWRC_MPI_OPTIONS_DUKE_HABANA); 
    aflowrc::load_default("MPI_COMMAND_DUKE_HABANA",AFLOWRC_MPI_COMMAND_DUKE_HABANA); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_HABANA",AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_QRATS_MPICH",AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH); 
    aflowrc::load_default("MPI_COMMAND_DUKE_QRATS_MPICH",AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_QRATS_MPICH",AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH); 

    aflowrc::load_default("MPI_OPTIONS_DUKE_QFLOW_OPENMPI",AFLOWRC_MPI_OPTIONS_DUKE_QFLOW_OPENMPI); 
    aflowrc::load_default("MPI_COMMAND_DUKE_QFLOW_OPENMPI",AFLOWRC_MPI_COMMAND_DUKE_QFLOW_OPENMPI); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI",AFLOWRC_MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI); 

    //CO20201220 X START
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_X",AFLOWRC_MPI_OPTIONS_DUKE_X_X); 
    aflowrc::load_default("MPI_COMMAND_DUKE_X_X",AFLOWRC_MPI_COMMAND_DUKE_X_X); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_X",AFLOWRC_MPI_BINARY_DIR_DUKE_X_X); 
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_CRAY",AFLOWRC_MPI_OPTIONS_DUKE_X_CRAY); 
    aflowrc::load_default("MPI_COMMAND_DUKE_X_CRAY",AFLOWRC_MPI_COMMAND_DUKE_X_CRAY); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_CRAY",AFLOWRC_MPI_BINARY_DIR_DUKE_X_CRAY); 
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_OLDCRAY",AFLOWRC_MPI_OPTIONS_DUKE_X_OLDCRAY); 
    aflowrc::load_default("MPI_COMMAND_DUKE_X_OLDCRAY",AFLOWRC_MPI_COMMAND_DUKE_X_OLDCRAY); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_OLDCRAY",AFLOWRC_MPI_BINARY_DIR_DUKE_X_OLDCRAY); 
    aflowrc::load_default("MPI_OPTIONS_DUKE_X_SMB",AFLOWRC_MPI_OPTIONS_DUKE_X_SMB); 
    aflowrc::load_default("MPI_COMMAND_DUKE_X_SMB",AFLOWRC_MPI_COMMAND_DUKE_X_SMB); 
    aflowrc::load_default("MPI_BINARY_DIR_DUKE_X_SMB",AFLOWRC_MPI_BINARY_DIR_DUKE_X_SMB); 
    //CO20201220 X STOP
    
    //CO20220818 JHU_ROCKFISH START
    aflowrc::load_default("MPI_OPTIONS_JHU_ROCKFISH",AFLOWRC_MPI_OPTIONS_JHU_ROCKFISH); 
    aflowrc::load_default("MPI_COMMAND_JHU_ROCKFISH",AFLOWRC_MPI_COMMAND_JHU_ROCKFISH); 
    aflowrc::load_default("MPI_BINARY_DIR_JHU_ROCKFISH",AFLOWRC_MPI_BINARY_DIR_JHU_ROCKFISH); 
    //CO20220818 JHU_ROCKFISH STOP

    //DX20190509 - MACHINE001 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE001",AFLOWRC_MPI_OPTIONS_MACHINE001); 
    aflowrc::load_default("MPI_COMMAND_MACHINE001",AFLOWRC_MPI_COMMAND_MACHINE001); 
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE001",AFLOWRC_MPI_BINARY_DIR_MACHINE001); 
    //DX20190509 - MACHINE001 - END

    //DX20190509 - MACHINE002 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE002",AFLOWRC_MPI_OPTIONS_MACHINE002);
    aflowrc::load_default("MPI_COMMAND_MACHINE002",AFLOWRC_MPI_COMMAND_MACHINE002);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE002",AFLOWRC_MPI_BINARY_DIR_MACHINE002);
    //DX20190509 - MACHINE002 - END

    //DX20201005 - MACHINE003 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE003",AFLOWRC_MPI_OPTIONS_MACHINE003);
    aflowrc::load_default("MPI_COMMAND_MACHINE003",AFLOWRC_MPI_COMMAND_MACHINE003);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE003",AFLOWRC_MPI_BINARY_DIR_MACHINE003);
    //DX20201005 - MACHINE003 - END

    //DX20211011 - MACHINE004 - START
    aflowrc::load_default("MPI_OPTIONS_MACHINE004",AFLOWRC_MPI_OPTIONS_MACHINE004);
    aflowrc::load_default("MPI_COMMAND_MACHINE004",AFLOWRC_MPI_COMMAND_MACHINE004);
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE004",AFLOWRC_MPI_BINARY_DIR_MACHINE004);
    //DX20211011 - MACHINE004 - END

    //DX20190107 - CMU EULER - START
    aflowrc::load_default("MPI_OPTIONS_CMU_EULER",AFLOWRC_MPI_OPTIONS_CMU_EULER); 
    aflowrc::load_default("MPI_COMMAND_CMU_EULER",AFLOWRC_MPI_COMMAND_CMU_EULER); 
    aflowrc::load_default("MPI_BINARY_DIR_CMU_EULER",AFLOWRC_MPI_BINARY_DIR_CMU_EULER); 
    //DX20190107 - CMU EULER - END

    aflowrc::load_default("MPI_OPTIONS_MPCDF_EOS",AFLOWRC_MPI_OPTIONS_MPCDF_EOS); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_EOS",AFLOWRC_MPI_COMMAND_MPCDF_EOS); 
    aflowrc::load_default("MPI_NCPUS_MPCDF_EOS",AFLOWRC_MPI_NCPUS_MPCDF_EOS); 
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_EOS",AFLOWRC_MPI_HYPERTHREADING_MPCDF_EOS); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_EOS",AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS); 

    aflowrc::load_default("MPI_OPTIONS_MPCDF_DRACO",AFLOWRC_MPI_OPTIONS_MPCDF_DRACO); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_DRACO",AFLOWRC_MPI_COMMAND_MPCDF_DRACO); 
    aflowrc::load_default("MPI_NCPUS_MPCDF_DRACO",AFLOWRC_MPI_NCPUS_MPCDF_DRACO); 
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_DRACO",AFLOWRC_MPI_HYPERTHREADING_MPCDF_DRACO); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_DRACO",AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO); 

    aflowrc::load_default("MPI_OPTIONS_MPCDF_COBRA",AFLOWRC_MPI_OPTIONS_MPCDF_COBRA); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_COBRA",AFLOWRC_MPI_COMMAND_MPCDF_COBRA); 
    aflowrc::load_default("MPI_NCPUS_MPCDF_COBRA",AFLOWRC_MPI_NCPUS_MPCDF_COBRA); 
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_COBRA",AFLOWRC_MPI_HYPERTHREADING_MPCDF_COBRA); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_COBRA",AFLOWRC_MPI_BINARY_DIR_MPCDF_COBRA); 

    aflowrc::load_default("MPI_OPTIONS_MPCDF_HYDRA",AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA); 
    aflowrc::load_default("MPI_COMMAND_MPCDF_HYDRA",AFLOWRC_MPI_COMMAND_MPCDF_HYDRA); 
    aflowrc::load_default("MPI_NCPUS_MPCDF_HYDRA",AFLOWRC_MPI_NCPUS_MPCDF_HYDRA); 
    aflowrc::load_default("MPI_HYPERTHREADING_MPCDF_HYDRA",AFLOWRC_MPI_HYPERTHREADING_MPCDF_HYDRA); 
    aflowrc::load_default("MPI_BINARY_DIR_MPCDF_HYDRA",AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA); 

    aflowrc::load_default("MPI_OPTIONS_FULTON_MARYLOU",AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU); 
    aflowrc::load_default("MPI_COMMAND_FULTON_MARYLOU",AFLOWRC_MPI_COMMAND_FULTON_MARYLOU); 
    aflowrc::load_default("MPI_BINARY_DIR_FULTON_MARYLOU",AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU); 

    aflowrc::load_default("MPI_OPTIONS_MACHINE1",AFLOWRC_MPI_OPTIONS_MACHINE1); 
    aflowrc::load_default("MPI_COMMAND_MACHINE1",AFLOWRC_MPI_COMMAND_MACHINE1); 
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE1",AFLOWRC_MPI_BINARY_DIR_MACHINE1); 

    aflowrc::load_default("MPI_OPTIONS_MACHINE2",AFLOWRC_MPI_OPTIONS_MACHINE2); 
    aflowrc::load_default("MPI_COMMAND_MACHINE2",AFLOWRC_MPI_COMMAND_MACHINE2); 
    aflowrc::load_default("MPI_BINARY_DIR_MACHINE2",AFLOWRC_MPI_BINARY_DIR_MACHINE2); 

    if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;

    return TRUE;
  }
} // namespace aflowrc


// ***************************************************************************
// aflowrc::write_default
// ***************************************************************************
namespace aflowrc {
  bool write_default(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    stringstream message;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " BEGIN" << endl;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.home=" << XHOST.home << endl;
    if(XHOST.aflowrc_filename.empty()) XHOST.aflowrc_filename=AFLOWRC_FILENAME_LOCAL;
    if(LDEBUG) oss << __AFLOW_FUNC__ << " XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;

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

    //HE20220218 START
    aflowrc << "// DEFAULTS ENTRY LOADER" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_AFLUX_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_AFLUX_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_RESTAPI_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_RESTAPI_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH << "\"" << endl;
    aflowrc << "DEFAULT_ENTRY_LOADER_FS_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH << "\"" << endl;
    aflowrc << " " << endl;
    //HE20220218 STOP

    //ME20191001 START
    aflowrc << "// DEFAULT AFLOW DATABASE" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_STATS_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_DATA_PATH=\"" << AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_LOCK_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AFLOW_DB_STALE_THRESHOLD=" << AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD << endl;
    aflowrc << " " << endl;
    //ME20191001 STOP
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
    aflowrc << "DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT << "\"" << endl; //DX20200129
    aflowrc << "DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON=\"" << AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON << "\"" << endl; //DX20200129
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
    aflowrc << "DEFAULT_AFLOW_PLASMONICS_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_PLASMONICS_FILE << "\"" << endl;

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
    aflowrc << "DEFAULT_UNARY_KPPRA=" << AFLOWRC_DEFAULT_UNARY_KPPRA << " // 32768 // 27000" << endl;
    aflowrc << "DEFAULT_UNARY_KPPRA_STATIC=" << AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC << "// 32768 // 27000" << endl;
    aflowrc << "DEFAULT_PHONONS_KSCHEME=\"" << AFLOWRC_DEFAULT_PHONONS_KSCHEME << "\"" << endl;
    aflowrc << "DEFAULT_PHONONS_KPPRA=" << AFLOWRC_DEFAULT_PHONONS_KPPRA << endl;
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
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_NELM=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM << endl; //CO20200624
    aflowrc << "DEFAULT_VASP_FORCE_OPTION_NELM_STATIC=" << AFLOWRC_DEFAULT_VASP_FORCE_OPTION_NELM_STATIC << endl; //CO20200624
    aflowrc << "MAX_VASP_NELM=" << AFLOWRC_MAX_VASP_NELM << endl; //CO20200624
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
    aflowrc << "DEFAULT_AFLOWDATA_WEB_DIRECTORY=\"" << AFLOWRC_DEFAULT_AFLOWDATA_WEB_DIRECTORY << "\"" << endl; //CO+ME20200731

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
    aflowrc << "// DEFAULTS CHULL" << endl;
    aflowrc << "DEFAULT_CHULL_ALLOWED_DFT_TYPES=\"" << AFLOWRC_DEFAULT_CHULL_ALLOWED_DFT_TYPES << "\"  // comma-separated list of dft_types to include (string match)" << endl;
    aflowrc << "DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES=" << AFLOWRC_DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES=" << AFLOWRC_DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES << " // INT" << endl;
    aflowrc << "DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS=" << AFLOWRC_DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES=" << AFLOWRC_DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES << " // INT" << endl;
    aflowrc << "DEFAULT_CHULL_OUTLIER_MULTIPLIER=" << AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER << " // DOUBLE" << endl;
    aflowrc << "DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED=" << AFLOWRC_DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED << " // 0 - FALSE (NOT recommended), 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_BANNER=" << AFLOWRC_DEFAULT_CHULL_LATEX_BANNER << " // 0 - no banner, 1 - full banner, 2 - small banner" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN=" << AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_STOICH_HEADER=" << AFLOWRC_DEFAULT_CHULL_LATEX_STOICH_HEADER << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_UNARIES=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNARIES << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_OFF_HULL=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL << " // -1 - default based on dimension/filter_by settings, 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_UNSTABLE=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FILTER_SCHEME=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME << "\"" << "  // Z-axis (also Energy-axis) or Distance; only reads first letter (case-insensitive)" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FILTER_VALUE=" << AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE << " // DOUBLE (filter scheme must be set)" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COLOR_BAR=" << AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_HEAT_MAP=" << AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COLOR_GRADIENT=" << AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_COLOR_MAP=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP << "\"" << "  // default - rgb(0pt)=(0,0,1); rgb(63pt)=(1,0.644,0) (latex pgfplots color maps)" << endl; //: http://www.phy.ntnu.edu.tw/demolab/doc/texlive-pictures-doc/latex/pgfplots/pgfplots.pdf page 58)" << endl;  //the url confuses RemoveComment(), perhaps we recursive remove comments? (might be dangerous)
    aflowrc << "DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR << "\"" << "  // white, black, red, green, blue, cyan, magenta, or yellow" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_REVERSE_AXIS=" << AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW=" << AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_LINKS=" << AFLOWRC_DEFAULT_CHULL_LATEX_LINKS << " // 0 - no links whatsoever, 1 - internal and external links, 2 - external links only, 3 - internal links only" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_LABEL_NAME=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME << "\"" << "  // Compound, Prototype, Both (Compound and Prototype), ICSD, or None; only reads first letter (case-insensitive)" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_META_LABELS=" << AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_LABELS_OFF_HULL=" << AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION=" << AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION << " // -1 - default based on dimension/label settings, 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_HELVETICA_FONT=" << AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_FONT_SIZE=\"" << AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE << "\"" << "// tiny, scriptsize, footnotesize, small, normalsize, large (default), Large (default large), LARGE, huge (default large Helvetica), or Huge; fontsize{5}{7}\\\\selectfont also works" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_ROTATE_LABELS=" << AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS << " // 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_LATEX_BOLD_LABELS=" << AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS << " // -1 - default: no bold unless the compound is a ternary, 0 - FALSE, 1 - TRUE" << endl;
    aflowrc << "DEFAULT_CHULL_PNG_RESOLUTION=" << AFLOWRC_DEFAULT_CHULL_PNG_RESOLUTION << " // INT" << endl;

    aflowrc << " " << endl; //CO20190628
    aflowrc << "// DEFAULTS GFA" << endl; //CO20190628
    aflowrc << "DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF=" << AFLOWRC_DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF << " // DOUBLE in eV" << endl; //CO20190628

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS ARUN" << endl;
    aflowrc << "ARUN_DIRECTORY_PREFIX=\"" << AFLOWRC_ARUN_DIRECTORY_PREFIX << "\"" << endl;

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS POCC" << endl;
    aflowrc << "DEFAULT_POCC_STRUCTURE_GENERATION_ALGO=\"" << AFLOWRC_DEFAULT_POCC_STRUCTURE_GENERATION_ALGO << "\"" << " // UFF" << endl;
    aflowrc << "DEFAULT_POCC_TEMPERATURE_STRING=\"" << AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING << "\"" << endl;
    aflowrc << "DEFAULT_POCC_EXCLUDE_UNSTABLE=" << AFLOWRC_DEFAULT_POCC_EXCLUDE_UNSTABLE << endl;  //ME20210927
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
    aflowrc << "POCC_OUT_FILE=\"" << AFLOWRC_POCC_OUT_FILE << "\"" << endl;
    aflowrc << "POCC_APL_OUT_FILE=\"" << AFLOWRC_POCC_APL_OUT_FILE << "\"" << endl;  //ME20210927
    aflowrc << "POCC_ALL_SUPERCELLS_FILE=\"" << AFLOWRC_POCC_ALL_SUPERCELLS_FILE << "\"" << endl;
    aflowrc << "POCC_UNIQUE_SUPERCELLS_FILE=\"" << AFLOWRC_POCC_UNIQUE_SUPERCELLS_FILE << "\"" << endl;
    aflowrc << "POCC_ALL_HNF_MATRICES_FILE=\"" << AFLOWRC_POCC_ALL_HNF_MATRICES_FILE << "\"" << endl;
    aflowrc << "POCC_ALL_SITE_CONFIGURATIONS_FILE=\"" << AFLOWRC_POCC_ALL_SITE_CONFIGURATIONS_FILE << "\"" << endl;
    aflowrc << "POCC_DOSCAR_FILE=\"" << AFLOWRC_POCC_DOSCAR_FILE << "\"" << endl;
    aflowrc << "POCC_PHDOSCAR_FILE=\"" << AFLOWRC_POCC_PHDOSCAR_FILE << "\"" << endl;  //ME20210927
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
    aflowrc << "DEFAULT_APL_RELAX_COMMENSURATE=" << AFLOWRC_DEFAULT_APL_RELAX_COMMENSURATE << endl;  //ME20200427
    aflowrc << "DEFAULT_APL_ZEROSTATE=" << AFLOWRC_DEFAULT_APL_ZEROSTATE << endl;
    aflowrc << "DEFAULT_APL_ZEROSTATE_CHGCAR=" << AFLOWRC_DEFAULT_APL_ZEROSTATE_CHGCAR << endl;  //ME20191029
    aflowrc << "DEFAULT_APL_USE_LEPSILON=" << AFLOWRC_DEFAULT_APL_USE_LEPSILON << endl;
    aflowrc << "DEFAULT_APL_FREQFORMAT=\"" << AFLOWRC_DEFAULT_APL_FREQFORMAT << "\"" << endl;
    aflowrc << "DEFAULT_APL_DC=" << AFLOWRC_DEFAULT_APL_DC << endl;
    aflowrc << "DEFAULT_APL_DCPATH=\"" << AFLOWRC_DEFAULT_APL_DCPATH << "\"" << endl;
    aflowrc << "DEFAULT_APL_DCPOINTS=" << AFLOWRC_DEFAULT_APL_DCPOINTS << endl; //CO20181226
    aflowrc << "DEFAULT_APL_DOS=" << AFLOWRC_DEFAULT_APL_DOS << endl;
    aflowrc << "DEFAULT_APL_DOSMETHOD=\"" << AFLOWRC_DEFAULT_APL_DOSMETHOD << "\"" << endl;
    aflowrc << "DEFAULT_APL_DOSMESH=\"" << AFLOWRC_DEFAULT_APL_DOSMESH << "\"" << endl;
    aflowrc << "DEFAULT_APL_DOSPOINTS=" << AFLOWRC_DEFAULT_APL_DOSPOINTS << endl;
    aflowrc << "DEFAULT_APL_DOSSMEAR=" << AFLOWRC_DEFAULT_APL_DOSSMEAR << endl;
    aflowrc << "DEFAULT_APL_DOS_PROJECT=" << AFLOWRC_DEFAULT_APL_DOS_PROJECT << endl; //ME20200421
    aflowrc << "DEFAULT_APL_TP=" << AFLOWRC_DEFAULT_APL_TP << endl;
    aflowrc << "DEFAULT_APL_DISPLACEMENTS=" << AFLOWRC_DEFAULT_APL_DISPLACEMENTS << endl;  //ME20200421
    aflowrc << "DEFAULT_APL_TPT=\"" << AFLOWRC_DEFAULT_APL_TPT << "\"" << endl;
    aflowrc << "DEFAULT_APL_GVEL=" << AFLOWRC_DEFAULT_APL_GVEL << endl;  //ME20200517
    aflowrc << "DEFAULT_APL_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_APL_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_APL_OUT_FILE=\"" << AFLOWRC_DEFAULT_APL_OUT_FILE << "\"" << endl;  //ME20210927
    aflowrc << "DEFAULT_APL_PDIS_FILE=\"" << AFLOWRC_DEFAULT_APL_PDIS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PDOS_FILE=\"" << AFLOWRC_DEFAULT_APL_PDOS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_THERMO_FILE=\"" << AFLOWRC_DEFAULT_APL_THERMO_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_THERMO_JSON=\"" << AFLOWRC_DEFAULT_APL_THERMO_JSON << "\"" << endl;  //ME20211019
    aflowrc << "DEFAULT_APL_DYNMAT_FILE=\"" << AFLOWRC_DEFAULT_APL_DYNMAT_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_HARMIFC_FILE=\"" << AFLOWRC_DEFAULT_APL_HARMIFC_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_POLAR_FILE=\"" << AFLOWRC_DEFAULT_APL_POLAR_FILE << "\"" << endl;  //ME20200415
    aflowrc << "DEFAULT_APL_HSKPTS_FILE=\"" << AFLOWRC_DEFAULT_APL_HSKPTS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_MSQRDISP_FILE=\"" << AFLOWRC_DEFAULT_APL_MSQRDISP_FILE << "\"" << endl;  //ME20200329
    aflowrc << "DEFAULT_APL_GVEL_FILE=\"" << AFLOWRC_DEFAULT_APL_GVEL_FILE << "\"" << endl;  // ME20200517
    //ME20190614 BEGIN
    aflowrc << "DEFAULT_APL_PHDOSCAR_FILE=\"" << AFLOWRC_DEFAULT_APL_PHDOSCAR_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PHPOSCAR_FILE=\"" << AFLOWRC_DEFAULT_APL_PHPOSCAR_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PHKPOINTS_FILE=\"" << AFLOWRC_DEFAULT_APL_PHKPOINTS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_APL_PHEIGENVAL_FILE=\"" << AFLOWRC_DEFAULT_APL_PHEIGENVAL_FILE << "\"" << endl;
    //ME20190614 END
    aflowrc << "DEFAULT_APL_STATE_FILE=\"" << AFLOWRC_DEFAULT_APL_STATE_FILE << "\"" << endl;  //ME20200224
    //ME20200329 BEGIN
    aflowrc << "DEFAULT_APL_ADISP_SCENE_FORMAT=\"" << AFLOWRC_DEFAULT_APL_ADISP_SCENE_FORMAT << "\"" << endl;
    aflowrc << "DEFAULT_APL_ADISP_AMPLITUDE=" << AFLOWRC_DEFAULT_APL_ADISP_AMPLITUDE << endl;
    aflowrc << "DEFAULT_APL_ADISP_NSTEPS=" << AFLOWRC_DEFAULT_APL_ADISP_NSTEPS << endl;
    aflowrc << "DEFAULT_APL_ADISP_NPERIODS=" << AFLOWRC_DEFAULT_APL_ADISP_NPERIODS << endl;
    //ME20200329 END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS QHA" << endl;
    aflowrc << "DEFAULT_QHA_MODE=\"" << AFLOWRC_DEFAULT_QHA_MODE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_EOS=" << AFLOWRC_DEFAULT_QHA_EOS  << endl;
    aflowrc << "DEFAULT_QHA_EOS_DISTORTION_RANGE=\"" << AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_EOS_MODEL=\"" << AFLOWRC_DEFAULT_QHA_EOS_MODEL << "\"" << endl;//AS20200818
    aflowrc << "DEFAULT_QHA_GP_DISTORTION=" << AFLOWRC_DEFAULT_QHA_GP_DISTORTION  << endl;
    aflowrc << "DEFAULT_QHA_TAYLOR_EXPANSION_ORDER=" << AFLOWRC_DEFAULT_QHA_TAYLOR_EXPANSION_ORDER  << endl;//AS20200602
    aflowrc << "DEFAULT_QHA_INCLUDE_ELEC_CONTRIB=" << AFLOWRC_DEFAULT_QHA_INCLUDE_ELEC_CONTRIB  << endl;
    aflowrc << "DEFAULT_QHA_SOMMERFELD_EXPANSION=" << AFLOWRC_DEFAULT_QHA_SOMMERFELD_EXPANSION  << endl;//AS20200528
    aflowrc << "DEFAULT_QHA_PDIS_T=\"" << AFLOWRC_DEFAULT_QHA_PDIS_T << "\"" << endl;
    //AS20200508 BEGIN
    aflowrc << "DEFAULT_QHA_GP_FINITE_DIFF=" << AFLOWRC_DEFAULT_QHA_GP_FINITE_DIFF  << endl;
    aflowrc << "DEFAULT_QHA_IGNORE_IMAGINARY=" << AFLOWRC_DEFAULT_QHA_IGNORE_IMAGINARY  << endl;
    aflowrc << "DEFAULT_QHA_RELAX_IONS_CELL=" << AFLOWRC_DEFAULT_QHA_RELAX_IONS_CELL  << endl;//AS20201123
    aflowrc << "DEFAULT_QHA_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_QHA_FILE_PREFIX << "\"" << endl;
    //AS20200709 BEGIN
    aflowrc << "DEFAULT_QHA3P_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_QHA3P_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_QHANP_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_QHANP_FILE_PREFIX << "\"" << endl;
    aflowrc << "DEFAULT_SCQHA_FILE_PREFIX=\"" << AFLOWRC_DEFAULT_SCQHA_FILE_PREFIX << "\"" << endl;
    //AS20200709 END
    aflowrc << "DEFAULT_QHA_GP_PATH_FILE=\"" << AFLOWRC_DEFAULT_QHA_GP_PATH_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_GP_MESH_FILE=\"" << AFLOWRC_DEFAULT_QHA_GP_MESH_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_GP_AVG_FILE=\"" << AFLOWRC_DEFAULT_QHA_GP_AVG_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_THERMO_FILE=\"" << AFLOWRC_DEFAULT_QHA_THERMO_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_FREQS_FILE=\"" << AFLOWRC_DEFAULT_QHA_FREQS_FILE << "\"" << endl;
    aflowrc << "DEFAULT_QHA_FVT_FILE=\"" << AFLOWRC_DEFAULT_QHA_FVT_FILE << "\"" << endl;
    //AS20200508 END
    aflowrc << "DEFAULT_QHA_COEFF_FILE=\"" << AFLOWRC_DEFAULT_QHA_COEFF_FILE << "\"" << endl;//AS20210517
    aflowrc << "DEFAULT_QHA_IMAG_FILE=\"" << AFLOWRC_DEFAULT_QHA_IMAG_FILE << "\"" << endl;//AS20210517
    aflowrc << "DEFAULT_QHA_PDIS_FILE=\"" << AFLOWRC_DEFAULT_QHA_PDIS_FILE << "\"" << endl;//AS20201022
    aflowrc << "DEFAULT_QHA_PDOS_FILE=\"" << AFLOWRC_DEFAULT_QHA_PDOS_FILE << "\"" << endl;//AS20201201
    aflowrc << "DEFAULT_QHA_KPOINTS_FILE=\"" << AFLOWRC_DEFAULT_QHA_KPOINTS_FILE << "\"" << endl;//AS20201112
    //AS20210914 BEGIN
    aflowrc << "DEFAULT_POCC_QHA_THERMO_FILE=\"" << AFLOWRC_DEFAULT_POCC_QHA_THERMO_FILE << "\"" << endl;
    aflowrc << "DEFAULT_POCC_QHA_AVGTHERMO_FILE=\"" << AFLOWRC_DEFAULT_POCC_QHA_AVGTHERMO_FILE << "\"" << endl;
    //AS20210914 END

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
    aflowrc << "DEFAULT_AAPL_PS_FILE=\"" << AFLOWRC_DEFAULT_AAPL_PS_FILE << "\"" << endl;  //ME20191104
    aflowrc << "DEFAULT_AAPL_GRUENEISEN_FILE=\"" << AFLOWRC_DEFAULT_AAPL_GRUENEISEN_FILE << "\"" << endl;  //ME20191104
    aflowrc << "DEFAULT_AAPL_RATES_FILE=\"" << AFLOWRC_DEFAULT_AAPL_RATES_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_RATES_3RD_FILE=\"" << AFLOWRC_DEFAULT_AAPL_RATES_3RD_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_RATES_4TH_FILE=\"" << AFLOWRC_DEFAULT_AAPL_RATES_4TH_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_ISOTOPE_FILE=\"" << AFLOWRC_DEFAULT_AAPL_ISOTOPE_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_BOUNDARY_FILE=\"" << AFLOWRC_DEFAULT_AAPL_BOUNDARY_FILE << "\"" << endl;
    aflowrc << "DEFAULT_AAPL_TCOND_FILE=\"" << AFLOWRC_DEFAULT_AAPL_TCOND_FILE << "\"" << endl;
    //aflowrc << "DEFAULT_AAPL_TCOND_PLOT_FILE=\"" << AFLOWRC_DEFAULT_AAPL_TCOND_PLOT_FILE << "\"" << endl;  OBSOLETE ME20191104

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

    //SD20220323 - QCA START
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
    //SD20220323 - QCA END

    //RF20200413 START
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS CCE" << endl;
    aflowrc << "DEFAULT_CCE_OX_METHOD=" << AFLOWRC_DEFAULT_CCE_OX_METHOD << "" << "  // 1 - ELECTRONEGATIVITY_ALLEN, 2 - BADER" << endl;
    aflowrc << "DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION=" << AFLOWRC_DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION << "" << "  // 0.4 Ang tolerance between shortest and longest bonds when testing for multi-anion compound" << endl;
    aflowrc << "DEFAULT_CCE_OX_TOL=" << AFLOWRC_DEFAULT_CCE_OX_TOL << "" << "  // sum of oxidation states might not be exactly zero due to numerics" << endl;
    aflowrc << "DEFAULT_CCE_PEROX_CUTOFF=" << AFLOWRC_DEFAULT_CCE_PEROX_CUTOFF << "" << "  // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang" << endl;
    aflowrc << "DEFAULT_CCE_SUPEROX_CUTOFF=" << AFLOWRC_DEFAULT_CCE_SUPEROX_CUTOFF << "" << "  // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang" << endl;
    aflowrc << "DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF=" << AFLOWRC_DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF << "" << "  // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    aflowrc << "DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF=" << AFLOWRC_DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF << "" << "  // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    //RF20200413 END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS XTALFINDER" << endl;
    aflowrc << "DEFAULT_XTALFINDER_MISFIT_MATCH=" << AFLOWRC_DEFAULT_XTALFINDER_MISFIT_MATCH << " // values below this threshold: similar structures have similar properties" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_MISFIT_FAMILY=" << AFLOWRC_DEFAULT_XTALFINDER_MISFIT_FAMILY << " // values above this threshold: matched structures do not have similar properties" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_SUPERCELL_METHOD=" << AFLOWRC_DEFAULT_XTALFINDER_SUPERCELL_METHOD << " // // supercell method for comparing (robust, but slow, superceded by transformation method)" << endl; //DX20201223
    aflowrc << "DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING=" << AFLOWRC_DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING << " // factor that divides minimum interatomic distance" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_MATERIAL=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL << " // results file prefix" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_STRUCTURE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE << " // results file prefix" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_DUPLICATE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_DUPLICATE << " // results file prefix" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE << " // results file prefix" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE << " // results file prefix" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE << " // results file prefix" << endl; //DX20201118
    aflowrc << "DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE=" << AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE << " // results file prefix" << endl; //DX20201118

    //DX20200720 - START
    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS ANRL" << endl;
    aflowrc << "DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL=" << AFLOWRC_DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL << " // tolerance for equivalent Wyckoff coordinates" << endl;
    //DX20200720 - END

    aflowrc << " " << endl;
    aflowrc << "// DEFAULTS CORE" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_BEEP=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_BEEP << " // Celsius" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_HALT=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT << " // Celsius" << endl;
    aflowrc << "AFLOW_CORE_TEMPERATURE_REFRESH=" << AFLOWRC_AFLOW_CORE_TEMPERATURE_REFRESH << " // seconds"   << endl;

    aflowrc << "SECONDS_SLEEP_VASP_COMPLETION=" << AFLOWRC_SECONDS_SLEEP_VASP_COMPLETION << " // seconds"   << endl;  //CO20201111
    aflowrc << "SECONDS_SLEEP_VASP_MONITOR=" << AFLOWRC_SECONDS_SLEEP_VASP_MONITOR << " // seconds"   << endl;  //CO20201111
    aflowrc << "SECONDS_STALE_OUTCAR=" << AFLOWRC_SECONDS_STALE_OUTCAR << " // seconds"   << endl;  //CO20201111
    aflowrc << "BYTES_MAX_VASP_OUT=" << AFLOWRC_BYTES_MAX_VASP_OUT << " // bytes"   << endl;  //CO20201111
    aflowrc << "MEMORY_MAX_USAGE_RAM=" << AFLOWRC_MEMORY_MAX_USAGE_RAM << " // bytes"   << endl;  //CO20201111
    aflowrc << "MEMORY_MAX_USAGE_SWAP=" << AFLOWRC_MEMORY_MAX_USAGE_SWAP << " // bytes"   << endl;  //CO20201111
    aflowrc << "FILE_VASP_MONITOR=" << AFLOWRC_FILE_VASP_MONITOR << " // monitor file postfix"   << endl;  //CO20201111
    aflowrc << "INTEL_COMPILER_PATHS=" << AFLOWRC_INTEL_COMPILER_PATHS << " // comma-separated paths to search (for sourcing)"   << endl;  //CO20201111

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

    //CO20201220 X START
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
    //CO20201220 X STOP
    
    //CO20220818 JHU_ROCKFISH START
    aflowrc << "MPI_OPTIONS_JHU_ROCKFISH=\"" << AFLOWRC_MPI_OPTIONS_JHU_ROCKFISH << "\"" << "  // JHU_ROCKFISH" << endl;
    aflowrc << "MPI_COMMAND_JHU_ROCKFISH=\"" << AFLOWRC_MPI_COMMAND_JHU_ROCKFISH << "\"" << "  // JHU_ROCKFISH" << endl;
    aflowrc << "MPI_BINARY_DIR_JHU_ROCKFISH=\"" << AFLOWRC_MPI_BINARY_DIR_JHU_ROCKFISH << "\"" << "  // JHU_ROCKFISH" << endl; 
    //CO20220818 JHU_ROCKFISH STOP

    //DX20190509 - MACHINE001 - START
    aflowrc << "MPI_OPTIONS_MACHINE001=\"" << AFLOWRC_MPI_OPTIONS_MACHINE001 << "\"" << "// MACHINE001" << endl;
    aflowrc << "MPI_COMMAND_MACHINE001=\"" << AFLOWRC_MPI_COMMAND_MACHINE001 << "\"" << "// MACHINE001" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE001=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE001 << "\"" << "// MACHINE001" << endl; 
    //DX20190509 - MACHINE001 - END

    //DX20190509 - MACHINE002 - START
    aflowrc << "MPI_OPTIONS_MACHINE002=\"" << AFLOWRC_MPI_OPTIONS_MACHINE002 << "\"" << "// MACHINE002" << endl;
    aflowrc << "MPI_COMMAND_MACHINE002=\"" << AFLOWRC_MPI_COMMAND_MACHINE002 << "\"" << "// MACHINE002" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE002=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE002 << "\"" << "// MACHINE002" << endl;
    //DX20190509 - MACHINE002 - START

    //DX20201005 - MACHINE003 - START
    aflowrc << "MPI_OPTIONS_MACHINE003=\"" << AFLOWRC_MPI_OPTIONS_MACHINE003 << "\"" << "// MACHINE003" << endl;
    aflowrc << "MPI_COMMAND_MACHINE003=\"" << AFLOWRC_MPI_COMMAND_MACHINE003 << "\"" << "// MACHINE003" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE003=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE003 << "\"" << "// MACHINE003" << endl;
    //DX20201005 - MACHINE003 - START

    //DX20211011 - MACHINE004 - START
    aflowrc << "MPI_OPTIONS_MACHINE004=\"" << AFLOWRC_MPI_OPTIONS_MACHINE004 << "\"" << "// MACHINE004" << endl;
    aflowrc << "MPI_COMMAND_MACHINE004=\"" << AFLOWRC_MPI_COMMAND_MACHINE004 << "\"" << "// MACHINE004" << endl;
    aflowrc << "MPI_BINARY_DIR_MACHINE004=\"" << AFLOWRC_MPI_BINARY_DIR_MACHINE004 << "\"" << "// MACHINE004" << endl;
    //DX20211011 - MACHINE004 - START

    //DX20190107 - CMU EULER - START
    aflowrc << "MPI_OPTIONS_CMU_EULER=\"" << AFLOWRC_MPI_OPTIONS_CMU_EULER << "\"" << "// CMU_EULER" << endl;
    aflowrc << "MPI_COMMAND_CMU_EULER=\"" << AFLOWRC_MPI_COMMAND_CMU_EULER << "\"" << "// CMU_EULER" << endl;
    aflowrc << "MPI_BINARY_DIR_CMU_EULER=\"" << AFLOWRC_MPI_BINARY_DIR_CMU_EULER << "\"" << "// CMU_EULER" << endl; 
    //DX20190107 - CMU EULER - END

    aflowrc << "MPI_OPTIONS_MPCDF_EOS=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_EOS=\"" << AFLOWRC_MPI_COMMAND_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_EOS=\"" << AFLOWRC_MPI_NCPUS_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_EOS=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI        // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_EOS=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS << "\"" << "  // MPCDF_EOS_MPI" << endl; 

    aflowrc << "MPI_OPTIONS_MPCDF_DRACO=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_DRACO=\"" << AFLOWRC_MPI_COMMAND_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_DRACO=\"" << AFLOWRC_MPI_NCPUS_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_DRACO=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI        // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_DRACO=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO << "\"" << "  // MPCDF_DRACO_MPI" << endl; 

    aflowrc << "MPI_OPTIONS_MPCDF_COBRA=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_COBRA=\"" << AFLOWRC_MPI_COMMAND_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_COBRA=\"" << AFLOWRC_MPI_NCPUS_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_COBRA=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI        // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    aflowrc << "MPI_BINARY_DIR_MPCDF_COBRA=\"" << AFLOWRC_MPI_BINARY_DIR_MPCDF_COBRA << "\"" << "  // MPCDF_COBRA_MPI" << endl; 

    aflowrc << "MPI_OPTIONS_MPCDF_HYDRA=\"" << AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;
    aflowrc << "MPI_COMMAND_MPCDF_HYDRA=\"" << AFLOWRC_MPI_COMMAND_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;
    aflowrc << "MPI_NCPUS_MPCDF_HYDRA=\"" << AFLOWRC_MPI_NCPUS_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI" << endl;
    aflowrc << "MPI_HYPERTHREADING_MPCDF_HYDRA=\"" << AFLOWRC_MPI_HYPERTHREADING_MPCDF_HYDRA << "\"" << "  // MPCDF_HYDRA_MPI        // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
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

    //   XHOST.DEBUG=TRUE;
    //[CO20190808 - issue this ONLY if it was written, should fix www-data]cerr << "WARNING: aflowrc::write_default: WRITING default " << XHOST.aflowrc_filename << endl;
    if(aurostd::stringstream2file(aflowrc,XHOST.aflowrc_filename) && aurostd::FileExist(XHOST.aflowrc_filename)){
      if(!XHOST.vflag_control.flag("WWW")){ //CO20200404 - new web flag
        //[CO20200404 - OBSOLETE]cerr << "WARNING: aflowrc::write_default: WRITING default " << XHOST.aflowrc_filename << endl;  //CO20190808 - issue this ONLY if it was written, should fix www-data
        message << "WRITING default " << XHOST.aflowrc_filename;pflow::logger(__AFLOW_FILE__,__AFLOW_FUNC__,message,std::cerr,_LOGGER_MESSAGE_);  //CO20200404 - LEAVE std::cerr here, FR needs this for web
      }
    }
    if(LDEBUG) oss << __AFLOW_FUNC__ << " END" << endl;
    return TRUE;
  }
} // namespace aflowrc

// ***************************************************************************
// aflowrc::print_aflowrc
// ***************************************************************************
namespace aflowrc {
  bool print_aflowrc(std::ostream& oss,bool AFLOWRC_VERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG || AFLOWRC_VERBOSE);   
    if(LDEBUG) oss << "aflowrc::print_aflowrc: BEGIN" << endl;
    if(LDEBUG) oss << "aflowrc::print_aflowrc: XHOST.home=" << XHOST.home << endl;
    if(LDEBUG) oss << "aflowrc::print_aflowrc: XHOST.aflowrc_filename=" << XHOST.aflowrc_filename << endl;

    //HE20220218 START
    if(LDEBUG) oss << "// DEFAULTS ENTRY LOADER" << endl;
    if(LDEBUG) oss << "DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_ENTRY_LOADER_AFLUX_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_ENTRY_LOADER_AFLUX_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_ENTRY_LOADER_RESTAPI_SERVER=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_ENTRY_LOADER_RESTAPI_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_ENTRY_LOADER_FS_PATH=\"" << AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH << "\"" << endl;
    //HE20220218 STOP

    //ME20191001 START
    if(LDEBUG) oss << "// DEFAULT AFLOW DATABASE" << endl;
    if(LDEBUG) oss << "DEFAULT_AFLOW_DB_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_FILE << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_AFLOW_DB_STATS_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_AFLOW_DB_DATA_PATH=\"" << AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_AFLOW_DB_LOCK_FILE=\"" << AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE << "\"" << endl;
    if(LDEBUG) oss << "DEFAULT_AFLOW_DB_STALE_THRESHOLD=" << AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD << endl;
    //ME20191001 STOP
    if(LDEBUG) oss << "// DEFAULT DEFINITIONS" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_KZIP_BIN\")=\"" << DEFAULT_KZIP_BIN << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_KZIP_EXT\")=\"" << DEFAULT_KZIP_EXT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_TMPFS_DIRECTORIES\")=\"" << DEFAULT_TMPFS_DIRECTORIES << "\"" << endl;

    if(LDEBUG) oss << "// FILENAMES FOR AFLOW.ORG ANALYSIS" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_AFLOWLIB_ENTRY_OUT\")=\"" << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_AFLOWLIB_ENTRY_JSON\")=\"" << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_EDATA_ORIG_OUT\")=\"" << DEFAULT_FILE_EDATA_ORIG_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_EDATA_RELAX_OUT\")=\"" << DEFAULT_FILE_EDATA_RELAX_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_EDATA_BANDS_OUT\")=\"" << DEFAULT_FILE_EDATA_BANDS_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_DATA_ORIG_OUT\")=\"" << DEFAULT_FILE_DATA_ORIG_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_DATA_RELAX_OUT\")=\"" << DEFAULT_FILE_DATA_RELAX_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_DATA_BANDS_OUT\")=\"" << DEFAULT_FILE_DATA_BANDS_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_EDATA_ORIG_JSON\")=\"" << DEFAULT_FILE_EDATA_ORIG_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_EDATA_RELAX_JSON\")=\"" << DEFAULT_FILE_EDATA_RELAX_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_EDATA_BANDS_JSON\")=\"" << DEFAULT_FILE_EDATA_BANDS_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_DATA_ORIG_JSON\")=\"" << DEFAULT_FILE_DATA_ORIG_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_DATA_RELAX_JSON\")=\"" << DEFAULT_FILE_DATA_RELAX_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_DATA_BANDS_JSON\")=\"" << DEFAULT_FILE_DATA_BANDS_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_TIME_OUT\")=\"" << DEFAULT_FILE_TIME_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_SPACEGROUP1_OUT\")=\"" << DEFAULT_FILE_SPACEGROUP1_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_SPACEGROUP2_OUT\")=\"" << DEFAULT_FILE_SPACEGROUP2_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_VOLDISTPARAMS_OUT\")=\"" << DEFAULT_FILE_VOLDISTPARAMS_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_FILE_VOLDISTEVOLUTION_OUT\")=\"" << DEFAULT_FILE_VOLDISTEVOLUTION_OUT << "\"" << endl;

    if(LDEBUG) oss << "// FILENAMES FOR AFLOW OPERATION" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT\")=\"" << DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PRESCRIPT_OUT\")=\"" << DEFAULT_AFLOW_PRESCRIPT_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PRESCRIPT_COMMAND\")=\"" << DEFAULT_AFLOW_PRESCRIPT_COMMAND << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_POSTSCRIPT_OUT\")=\"" << DEFAULT_AFLOW_POSTSCRIPT_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_POSTSCRIPT_COMMAND\")=\"" << DEFAULT_AFLOW_POSTSCRIPT_COMMAND << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUP_OUT\")=\"" << DEFAULT_AFLOW_PGROUP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUP_JSON\")=\"" << DEFAULT_AFLOW_PGROUP_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUP_XTAL_OUT\")=\"" << DEFAULT_AFLOW_PGROUP_XTAL_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUP_XTAL_JSON\")=\"" << DEFAULT_AFLOW_PGROUP_XTAL_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT\")=\"" << DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT << "\"" << endl; //DX20200129
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON\")=\"" << DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON << "\"" << endl; //DX20200129
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUPK_OUT\")=\"" << DEFAULT_AFLOW_PGROUPK_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUPK_JSON\")=\"" << DEFAULT_AFLOW_PGROUPK_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUPK_XTAL_OUT\")=\"" << DEFAULT_AFLOW_PGROUPK_XTAL_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PGROUPK_XTAL_JSON\")=\"" << DEFAULT_AFLOW_PGROUPK_XTAL_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_FGROUP_OUT\")=\"" << DEFAULT_AFLOW_FGROUP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_FGROUP_JSON\")=\"" << DEFAULT_AFLOW_FGROUP_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_SGROUP_OUT\")=\"" << DEFAULT_AFLOW_SGROUP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_SGROUP_JSON\")=\"" << DEFAULT_AFLOW_SGROUP_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_AGROUP_OUT\")=\"" << DEFAULT_AFLOW_AGROUP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_AGROUP_JSON\")=\"" << DEFAULT_AFLOW_AGROUP_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_IATOMS_OUT\")=\"" << DEFAULT_AFLOW_IATOMS_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_IATOMS_JSON\")=\"" << DEFAULT_AFLOW_IATOMS_JSON << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_ICAGES_OUT\")=\"" << DEFAULT_AFLOW_ICAGES_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_SURFACE_OUT\")=\"" << DEFAULT_AFLOW_SURFACE_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_QMVASP_OUT\")=\"" << DEFAULT_AFLOW_QMVASP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_ERVASP_OUT\")=\"" << DEFAULT_AFLOW_ERVASP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_IMMISCIBILITY_OUT\")=\"" << DEFAULT_AFLOW_IMMISCIBILITY_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_MEMORY_OUT\")=\"" << DEFAULT_AFLOW_MEMORY_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_FROZSL_INPUT_OUT\")=\"" << DEFAULT_AFLOW_FROZSL_INPUT_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_FROZSL_POSCAR_OUT\")=\"" << DEFAULT_AFLOW_FROZSL_POSCAR_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_FROZSL_MODES_OUT\")=\"" << DEFAULT_AFLOW_FROZSL_MODES_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_FROZSL_EIGEN_OUT\")=\"" << DEFAULT_AFLOW_FROZSL_EIGEN_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_END_OUT\")=\"" << DEFAULT_AFLOW_END_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PLASMONICS_FILE\")=\"" << DEFAULT_AFLOW_PLASMONICS_FILE << "\"" << endl;

    if(LDEBUG) oss << "// DEFAULT GENERIC MPI" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_START_DEFAULT\")=\"" << MPI_START_DEFAULT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_STOP_DEFAULT\")=\"" << MPI_STOP_DEFAULT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DEFAULT\")=\"" << MPI_COMMAND_DEFAULT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_DEFAULT\")=" << MPI_NCPUS_DEFAULT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MAX\")=" << MPI_NCPUS_MAX << endl;

    if(LDEBUG) oss << "// DEFAULTS BINARY" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_GAMMA_BIN\")=\"" << DEFAULT_VASP_GAMMA_BIN << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_GAMMA_MPI_BIN\")=\"" << DEFAULT_VASP_GAMMA_MPI_BIN << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_BIN\")=\"" << DEFAULT_VASP_BIN << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_MPI_BIN\")=\"" << DEFAULT_VASP_MPI_BIN << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP5_BIN\")=\"" << DEFAULT_VASP5_BIN << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP5_MPI_BIN\")=\"" << DEFAULT_VASP5_MPI_BIN << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AIMS_BIN\")=\"" << DEFAULT_AIMS_BIN << "\""   << endl;

    if(LDEBUG) oss << "// DEFAULTS POTCARS" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIRECTORIES\")=\"" << DEFAULT_VASP_POTCAR_DIRECTORIES << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DATE\")=\"" << DEFAULT_VASP_POTCAR_DATE << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_SUFFIX\")=\"" << DEFAULT_VASP_POTCAR_SUFFIX << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DATE_POT_LDA\")=\"" << DEFAULT_VASP_POTCAR_DATE_POT_LDA << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DATE_POT_GGA\")=\"" << DEFAULT_VASP_POTCAR_DATE_POT_GGA << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POT_LDA\")=\"" << DEFAULT_VASP_POTCAR_DIR_POT_LDA << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POT_GGA\")=\"" << DEFAULT_VASP_POTCAR_DIR_POT_GGA << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POT_PBE\")=\"" << DEFAULT_VASP_POTCAR_DIR_POT_PBE << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA\")=\"" << DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA\")=\"" << DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE\")=\"" << DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN\")=\"" << DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"VASP_PSEUDOPOTENTIAL_DIRECTORY_POTPAW_PBE_KIN_\")=\"" << DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN << "\"" << endl;

    if(LDEBUG) oss << "// DEFAULT KPOINTS/DOS" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_BANDS_GRID\")=" << DEFAULT_BANDS_GRID << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_BANDS_LATTICE\")=\"" << DEFAULT_BANDS_LATTICE << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_KSCHEME\")=\"" << DEFAULT_KSCHEME << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA\")=" << DEFAULT_KPPRA << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_STATIC_KSCHEME\")=\"" << DEFAULT_STATIC_KSCHEME << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_STATIC\")=" << DEFAULT_KPPRA_STATIC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_KPPRA_ICSD\")=" << DEFAULT_KPPRA_ICSD << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_BANDS_GRID\")=" << DEFAULT_UNARY_BANDS_GRID << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA\")=" << DEFAULT_UNARY_KPPRA << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_UNARY_KPPRA_STATIC\")=" << DEFAULT_UNARY_KPPRA_STATIC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_PHONONS_KSCHEME\")=\"" << DEFAULT_PHONONS_KSCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_PHONONS_KPPRA\")=" << DEFAULT_PHONONS_KPPRA << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_EMIN\")=" << DEFAULT_DOS_EMIN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_EMAX\")=" << DEFAULT_DOS_EMAX << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_DOS_SCALE\")=" << DEFAULT_DOS_SCALE << endl;

    if(LDEBUG) oss << "// DEFAULT PRECISION" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_LOW\")=" << DEFAULT_VASP_PREC_ENMAX_LOW << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_MEDIUM\")=" << DEFAULT_VASP_PREC_ENMAX_MEDIUM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_NORMAL\")=" << DEFAULT_VASP_PREC_ENMAX_NORMAL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_HIGH\")=" << DEFAULT_VASP_PREC_ENMAX_HIGH << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_ENMAX_ACCURATE\")=" << DEFAULT_VASP_PREC_ENMAX_ACCURATE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_ENMAX_MINIMUM\")=" << DEFAULT_VASP_ENMAX_MINIMUM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_SPIN_REMOVE_CUTOFF\")=" << DEFAULT_VASP_SPIN_REMOVE_CUTOFF << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_POTIM\")=" << DEFAULT_VASP_PREC_POTIM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_PREC_EDIFFG\")=" << DEFAULT_VASP_PREC_EDIFFG << endl;

    if(LDEBUG) oss << "// DEFAULTS OPTIONS " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_OUT\")=\"" << DEFAULT_VASP_OUT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_INCAR\")=\"" << DEFAULT_VASP_EXTERNAL_INCAR << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_POSCAR\")=\"" << DEFAULT_VASP_EXTERNAL_POSCAR << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_POTCAR\")=\"" << DEFAULT_VASP_EXTERNAL_POTCAR << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_EXTERNAL_KPOINT\")=\"" << DEFAULT_VASP_EXTERNAL_KPOINTS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AIMS_EXTERNAL_CONTROL\")=\"" << DEFAULT_AIMS_EXTERNAL_CONTROL << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AIMS_EXTERNAL_GEOM\")=\"" << DEFAULT_AIMS_EXTERNAL_GEOM << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE\")=\"" << DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT\")=" << DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_FORCE_OPTION_SIGMA\")=\"" << DEFAULT_VASP_FORCE_OPTION_SIGMA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC\")=\"" << DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS\")=\"" << DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_NELM\")=" << DEFAULT_VASP_FORCE_OPTION_NELM << endl; //CO20200624
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_VASP_FORCE_OPTION_NELM_STATIC\")=" << DEFAULT_VASP_FORCE_OPTION_NELM_STATIC << endl; //CO20200624
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MAX_VASP_NELM\")=" << MAX_VASP_NELM << endl; //CO20200624
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME\")=\"" << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SYM\")=" << DEFAULT_VASP_FORCE_OPTION_SYM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2\")=" << DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_BADER\")=" << DEFAULT_VASP_FORCE_OPTION_BADER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_BADER_STATIC\")=" << DEFAULT_VASP_FORCE_OPTION_BADER_STATIC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_ELF\")=" << DEFAULT_VASP_FORCE_OPTION_ELF << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM\")=" << DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_WAVECAR\")=" << DEFAULT_VASP_FORCE_OPTION_WAVECAR << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_CHGCAR\")=" << DEFAULT_VASP_FORCE_OPTION_CHGCAR << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_VASP_FORCE_OPTION_LSCOUPLING\")=" << DEFAULT_VASP_FORCE_OPTION_LSCOUPLING << endl;

    if(LDEBUG) oss << "// AFLOW_LIBRARY AFLOW_PROJECT" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_LIBRARY_DIRECTORIES\")=\"" << DEFAULT_AFLOW_LIBRARY_DIRECTORIES << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOW_PROJECTS_DIRECTORIES\")=\"" << DEFAULT_AFLOW_PROJECTS_DIRECTORIES << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AFLOWDATA_WEB_DIRECTORY\")=\"" << DEFAULT_AFLOWDATA_WEB_DIRECTORY << "\"" << endl;  //CO+ME20200731

    if(LDEBUG) oss << "// DEFAULT PLATON/FINDSYM" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_PLATON_P_EQUAL\")=" << DEFAULT_PLATON_P_EQUAL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_PLATON_P_EXACT\")=" << DEFAULT_PLATON_P_EXACT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_ANG\")=" << DEFAULT_PLATON_P_ANG << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D1\")=" << DEFAULT_PLATON_P_D1 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D2\")=" << DEFAULT_PLATON_P_D2 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_PLATON_P_D3\")=" << DEFAULT_PLATON_P_D3 << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_FINDSYM_TOL\")=" << DEFAULT_FINDSYM_TOL << endl;

    if(LDEBUG) oss << "// DEFAULT GNUPLOT" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT\")=\"" << DEFAULT_GNUPLOT_EPS_FONT << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT_BOLD\")=\"" << DEFAULT_GNUPLOT_EPS_FONT_BOLD << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT_ITALICS\")=\"" << DEFAULT_GNUPLOT_EPS_FONT_ITALICS << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS\")=\"" << DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT\")=\"" << DEFAULT_GNUPLOT_PNG_FONT << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_BOLD\")=\"" << DEFAULT_GNUPLOT_PNG_FONT_BOLD << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_ITALICS\")=\"" << DEFAULT_GNUPLOT_PNG_FONT_ITALICS << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS\")=\"" << DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_GREEK_FONT\")=\"" << DEFAULT_GNUPLOT_GREEK_FONT << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_GREEK_FONT_BOLD\")=\"" << DEFAULT_GNUPLOT_GREEK_FONT_BOLD << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_ITALICS\")=\"" << DEFAULT_GNUPLOT_GREEK_FONT_ITALICS << "\""   << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS\")=\"" << DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS << "\""   << endl;

    if(LDEBUG) oss << "// DEFAULT CHULL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CHULL_ALLOWED_DFT_TYPES\")=" << DEFAULT_CHULL_ALLOWED_DFT_TYPES << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES\")=" << DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES\")=" << DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS\")=" << DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES\")=" << DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_OUTLIER_MULTIPLIER\")=" << DEFAULT_CHULL_OUTLIER_MULTIPLIER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED\")=" << DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_LATEX_BANNER\")=" << DEFAULT_CHULL_LATEX_BANNER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN\")=" << DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_STOICH_HEADER\")=" << DEFAULT_CHULL_LATEX_STOICH_HEADER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_PLOT_UNARIES\")=" << DEFAULT_CHULL_LATEX_PLOT_UNARIES << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_LATEX_PLOT_OFF_HULL\")=" << DEFAULT_CHULL_LATEX_PLOT_OFF_HULL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_PLOT_UNSTABLE\")=" << DEFAULT_CHULL_LATEX_PLOT_UNSTABLE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CHULL_LATEX_FILTER_SCHEME\")=\"" << DEFAULT_CHULL_LATEX_FILTER_SCHEME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"DEFAULT_CHULL_LATEX_FILTER_VALUE\")=" << DEFAULT_CHULL_LATEX_FILTER_VALUE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_COLOR_BAR\")=" << DEFAULT_CHULL_LATEX_COLOR_BAR << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_HEAT_MAP\")=" << DEFAULT_CHULL_LATEX_HEAT_MAP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_COLOR_GRADIENT\")=" << DEFAULT_CHULL_LATEX_COLOR_GRADIENT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CHULL_LATEX_COLOR_MAP\")=\"" << DEFAULT_CHULL_LATEX_COLOR_MAP << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR\")=\"" << DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_REVERSE_AXIS\")=" << DEFAULT_CHULL_LATEX_REVERSE_AXIS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW\")=" << DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_LATEX_LINKS\")=" << DEFAULT_CHULL_LATEX_LINKS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CHULL_LATEX_LABEL_NAME\")=\"" << DEFAULT_CHULL_LATEX_LABEL_NAME << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_META_LABELS\")=" << DEFAULT_CHULL_LATEX_META_LABELS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_LABELS_OFF_HULL\")=" << DEFAULT_CHULL_LATEX_LABELS_OFF_HULL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION\")=" << DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_HELVETICA_FONT\")=" << DEFAULT_CHULL_LATEX_HELVETICA_FONT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CHULL_LATEX_FONT_SIZE\")=\"" << DEFAULT_CHULL_LATEX_FONT_SIZE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<bool>(\"DEFAULT_CHULL_LATEX_ROTATE_LABELS\")=" << DEFAULT_CHULL_LATEX_ROTATE_LABELS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_LATEX_BOLD_LABELS\")=" << DEFAULT_CHULL_LATEX_BOLD_LABELS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"DEFAULT_CHULL_PNG_RESOLUTION\")=" << DEFAULT_CHULL_PNG_RESOLUTION << endl;

    if(LDEBUG) oss << "// DEFAULTS GFA" << endl;  //CO20190628
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF\")=" << DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF << endl;  //CO20190628

    if(LDEBUG) oss << "// DEFAULTS ARUN" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"ARUN_DIRECTORY_PREFIX\")=\"" << ARUN_DIRECTORY_PREFIX << "\"" << endl;

    if(LDEBUG) oss << "// DEFAULTS POCC" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_STRUCTURE_GENERATION_ALGO\")=\"" << DEFAULT_POCC_STRUCTURE_GENERATION_ALGO << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_TEMPERATURE_STRING\")=\"" << DEFAULT_POCC_TEMPERATURE_STRING << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_EXCLUDE_UNSTABLE\")=" << DEFAULT_POCC_EXCLUDE_UNSTABLE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_SITE_TOL\")=" << DEFAULT_POCC_SITE_TOL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_STOICH_TOL\")=" << DEFAULT_POCC_STOICH_TOL << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_UFF_BONDING_DISTANCE\")=" << DEFAULT_UFF_BONDING_DISTANCE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_UFF_ENERGY_TOLERANCE\")=" << DEFAULT_UFF_ENERGY_TOLERANCE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_UFF_CLUSTER_RADIUS\")=" << DEFAULT_UFF_CLUSTER_RADIUS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_RDF_RMAX\")=" << DEFAULT_POCC_RDF_RMAX << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_RDF_NBINS\")=" << DEFAULT_POCC_RDF_NBINS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON\")=" << DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS\")=" << DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_FILE_PREFIX\")=\"" << POCC_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_OUT_FILE\")=\"" << POCC_OUT_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_APL_OUT_FILE\")=\"" << POCC_APL_OUT_FILE << "\"" << endl;  //ME20210927
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_ALL_SUPERCELLS_FILE\")=\"" << POCC_ALL_SUPERCELLS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_UNIQUE_SUPERCELLS_FILE\")=\"" << POCC_UNIQUE_SUPERCELLS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_ALL_HNF_MATRICES_FILE\")=\"" << POCC_ALL_HNF_MATRICES_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_ALL_SITE_CONFIGURATIONS_FILE\")=\"" << POCC_ALL_SITE_CONFIGURATIONS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_DOSCAR_FILE\")=\"" << POCC_DOSCAR_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_PHDOSCAR_FILE\")=\"" << POCC_PHDOSCAR_FILE << "\"" << endl;  // ME20210927
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"POCC_ANIONS_LIST\")=\"" << POCC_ANIONS_LIST << "\"" << endl;

    if(LDEBUG) oss << "// DEFAULTS APL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PREC\")=\"" << DEFAULT_APL_PREC << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ENGINE\")=\"" << DEFAULT_APL_ENGINE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_HIBERNATE\")=" << DEFAULT_APL_HIBERNATE << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_MINSHELL\")=" << DEFAULT_APL_MINSHELL << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_MINATOMS\")=" << DEFAULT_APL_MINATOMS << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_POLAR\")=" << DEFAULT_APL_POLAR << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DMAG\")=" << DEFAULT_APL_DMAG << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DXYZONLY\")=" << DEFAULT_APL_DXYZONLY << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DSYMMETRIZE\")=" << DEFAULT_APL_DSYMMETRIZE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DINEQUIV_ONLY\")=" << DEFAULT_APL_DINEQUIV_ONLY << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DPM\")=\"" << DEFAULT_APL_DPM << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_RELAX\")=" << DEFAULT_APL_RELAX << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_RELAX_COMMENSURATE\")=" << DEFAULT_APL_RELAX_COMMENSURATE << endl;  //ME20200427
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ZEROSTATE\")=" << DEFAULT_APL_ZEROSTATE << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ZEROSTATE_CHGCAR\")=" << DEFAULT_APL_ZEROSTATE_CHGCAR << endl;  //ME20191029
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_USE_LEPSILON\")=" << DEFAULT_APL_USE_LEPSILON << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_FREQFORMAT\")=\"" << DEFAULT_APL_FREQFORMAT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DC\")=" << DEFAULT_APL_DC << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DCPATH\")=\"" << DEFAULT_APL_DCPATH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DCPOINTS\")=" << DEFAULT_APL_DCPOINTS << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOS\")=" << DEFAULT_APL_DOS << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOSMETHOD\")=\"" << DEFAULT_APL_DOSMETHOD << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOSMESH\")=\"" << DEFAULT_APL_DOSMESH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOSPOINTS\")=" << DEFAULT_APL_DOSPOINTS << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOSSMEAR\")=" << DEFAULT_APL_DOSSMEAR << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DOS_PROJECT\")=" << DEFAULT_APL_DOS_PROJECT << endl;  //ME20200213
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_TP\")=" << DEFAULT_APL_TP << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DISPLACEMENTS\")=" << DEFAULT_APL_DISPLACEMENTS << endl;  //ME20200421
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_TPT\")=\"" << DEFAULT_APL_TPT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_GVEL\")=" << DEFAULT_APL_GVEL << endl;  //ME20200517
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_FILE_PREFIX\")=\"" << DEFAULT_APL_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_OUT_FILE\")=\"" << DEFAULT_APL_OUT_FILE << "\"" << endl;  // ME20210927
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PDIS_FILE\")=\"" << DEFAULT_APL_PDIS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PDOS_FILE\")=\"" << DEFAULT_APL_PDOS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_THERMO_FILE\")=\"" << DEFAULT_APL_THERMO_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_THERMO_JSON\")=\"" << DEFAULT_APL_THERMO_JSON << "\"" << endl;  //ME20211019
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_DYNMAT_FILE\")=\"" << DEFAULT_APL_DYNMAT_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_HARMIFC_FILE\")=\"" << DEFAULT_APL_HARMIFC_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_POLAR_FILE\")=\"" << DEFAULT_APL_POLAR_FILE << "\"" << endl;  //ME20200415
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_HSKPTS_FILE\")=\"" << DEFAULT_APL_HSKPTS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_MSQRDISP_FILE\")=\"" << DEFAULT_APL_MSQRDISP_FILE << "\"" << endl;  //ME20200329
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_GVEL_FILE\")=\"" << DEFAULT_APL_GVEL_FILE << "\"" << endl;  // ME20200517
    //ME20190614 BEGIN
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PHDOSCAR_FILE\")=\"" << DEFAULT_APL_PHDOSCAR_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PHPOSCAR_FILE\")=\"" << DEFAULT_APL_PHPOSCAR_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PHKPOINTS_FILE\")=\"" << DEFAULT_APL_PHKPOINTS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_PHEIGENVAL_FILE\")=\"" << DEFAULT_APL_PHEIGENVAL_FILE << "\"" << endl;
    //ME20190614 END
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_STATE_FILE\")=\"" << DEFAULT_APL_STATE_FILE << "\"" << endl;  //ME20200224
    //ME20200329 BEGIN
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_SCENE_FORMAT\")=\"" << DEFAULT_APL_ADISP_SCENE_FORMAT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_AMPLITUDE\")=" << DEFAULT_APL_ADISP_AMPLITUDE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_NSTEPS\")=" << DEFAULT_APL_ADISP_NSTEPS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_APL_ADISP_NPERIODS\")=" << DEFAULT_APL_ADISP_NPERIODS << endl;
    //ME20200329 END

    if(LDEBUG) oss << "// DEFAULTS QHA" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_MODE\")=\"" << DEFAULT_QHA_MODE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_EOS\")=" << DEFAULT_QHA_EOS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_EOS_DISTORTION_RANGE\")=\"" << DEFAULT_QHA_EOS_DISTORTION_RANGE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_EOS_MODEL\")=\"" << DEFAULT_QHA_EOS_MODEL << "\"" << endl;//AS20200818
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_DISTORTION\")=" << DEFAULT_QHA_GP_DISTORTION << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_TAYLOR_EXPANSION_ORDER\")=" << DEFAULT_QHA_TAYLOR_EXPANSION_ORDER << endl;//AS20200602
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_INCLUDE_ELEC_CONTRIB\")=" << DEFAULT_QHA_INCLUDE_ELEC_CONTRIB << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_SOMMERFELD_EXPANSION\")=" << DEFAULT_QHA_SOMMERFELD_EXPANSION << endl;//AS20200528
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_PDIS_T\")=\"" << DEFAULT_QHA_PDIS_T << "\"" << endl;
    //AS20200508 BEGIN
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_FINITE_DIFF\")=" << DEFAULT_QHA_GP_FINITE_DIFF << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_IGNORE_IMAGINARY\")=" << DEFAULT_QHA_IGNORE_IMAGINARY << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_RELAX_IONS_CELL\")=" << DEFAULT_QHA_IGNORE_IMAGINARY << endl;//AS20201123
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_FILE_PREFIX\")=\"" << DEFAULT_QHA_FILE_PREFIX << "\"" << endl;
    //AS20200709 BEGIN
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA3P_FILE_PREFIX\")=\"" << DEFAULT_QHA3P_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHANP_FILE_PREFIX\")=\"" << DEFAULT_QHANP_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_SCQHA_FILE_PREFIX\")=\"" << DEFAULT_SCQHA_FILE_PREFIX << "\"" << endl;
    //AS20200709 END
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_PATH_FILE\")=\"" << DEFAULT_QHA_GP_PATH_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_MESH_FILE\")=\"" << DEFAULT_QHA_GP_MESH_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_GP_AVG_FILE\")=\"" << DEFAULT_QHA_GP_AVG_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_THERMO_FILE\")=\"" << DEFAULT_QHA_THERMO_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_FREQS_FILE\")=\"" << DEFAULT_QHA_FREQS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_FVT_FILE\")=\"" << DEFAULT_QHA_FVT_FILE << "\"" << endl;
    //AS20200508 END
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_COEFF_FILE\")=\"" << DEFAULT_QHA_COEFF_FILE << "\"" << endl;//AS20210517
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_IMAG_FILE\")=\"" << DEFAULT_QHA_IMAG_FILE << "\"" << endl;//AS20210517
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_PDIS_FILE\")=\"" << DEFAULT_QHA_PDIS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_PDOS_FILE\")=\"" << DEFAULT_QHA_PDOS_FILE << "\"" << endl;//AS20201201
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QHA_KPOINTS_FILE\")=\"" << DEFAULT_QHA_KPOINTS_FILE << "\"" << endl;//AS20201112
    //AS20210914 BEGIN
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_QHA_THERMO_FILE\")=\"" << DEFAULT_POCC_QHA_THERMO_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_POCC_QHA_AVGTHERMO_FILE\")=\"" << DEFAULT_POCC_QHA_AVGTHERMO_FILE << "\"" << endl;
    //AS20210914 END
    if(LDEBUG) oss << "// DEFAULTS AAPL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_BTE\")=\"" << DEFAULT_AAPL_BTE << "\"" << endl;
    //[ME20181226]if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_BZMETHOD\")=\"" << DEFAULT_AAPL_BZMETHOD << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_FOURTH_ORDER\")=" << DEFAULT_AAPL_FOURTH_ORDER << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_CUT_RAD\")=\"" << DEFAULT_AAPL_CUT_RAD << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_CUT_SHELL\")=\"" << DEFAULT_AAPL_CUT_SHELL << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_THERMALGRID\")=\"" << DEFAULT_AAPL_THERMALGRID << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_TCT\")=\"" << DEFAULT_AAPL_TCT << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_SUMRULE\")=" << DEFAULT_AAPL_SUMRULE << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_SUMRULE_MAX_ITER\")=" << DEFAULT_AAPL_SUMRULE_MAX_ITER << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_MIXING_COEFFICIENT\")=" << DEFAULT_AAPL_MIXING_COEFFICIENT << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_ISOTOPE\")=" << DEFAULT_AAPL_ISOTOPE << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_BOUNDARY\")=" << DEFAULT_AAPL_BOUNDARY << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_CUMULATIVEK\")=" << DEFAULT_AAPL_CUMULATIVEK << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_NANO_SIZE\")=" << DEFAULT_AAPL_NANO_SIZE << endl;  //ME20190112
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_FILE_PREFIX\")=\"" << DEFAULT_AAPL_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_IRRQPTS_FILE\")=\"" << DEFAULT_AAPL_IRRQPTS_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_GVEL_FILE\")=\"" << DEFAULT_AAPL_GVEL_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_PS_FILE\")=\"" << DEFAULT_AAPL_PS_FILE << "\"" << endl;  //ME20191104
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_GRUENEISEN_FILE\")=\"" << DEFAULT_AAPL_GRUENEISEN_FILE << "\"" << endl;  //ME20191104
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_RATES_FILE\")=\"" << DEFAULT_AAPL_RATES_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_RATES_3RD_FILE\")=\"" << DEFAULT_AAPL_RATES_3RD_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_RATES_4TH_FILE\")=\"" << DEFAULT_AAPL_RATES_4TH_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_ISOTOPE_FILE\")=\"" << DEFAULT_AAPL_ISOTOPE_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_BOUNDARY_FILE\")=\"" << DEFAULT_AAPL_BOUNDARY_FILE << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_TCOND_FILE\")=\"" << DEFAULT_AAPL_TCOND_FILE << "\"" << endl;
    //if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AAPL_TCOND_PLOT_FILE\")=\"" << DEFAULT_AAPL_TCOND_PLOT_FILE << "\"" << endl;  OBSOLETE ME20190411

    if(LDEBUG) oss << "// DEFAULTS AEL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_STRAIN_SYMMETRY\")=" << AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NNORMAL_STRAINS\")=" << AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NSHEAR_STRAINS\")=" << AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NORMAL_STRAIN_STEP\")=" << AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_SHEAR_STRAIN_STEP\")=" << AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_ORIGIN_STRAIN_CALC\")=" << AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_ORIGIN_STRAIN_FIT\")=" << AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_RELAXED_STRUCT_FIT\")=" << AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NEG_STRAINS\")=" << AFLOWRC_DEFAULT_AEL_NEG_STRAINS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_NIND_STRAIN_DIRS\")=" << AFLOWRC_DEFAULT_AEL_NIND_STRAIN_DIRS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_VASPSYM\")=" << AFLOWRC_DEFAULT_AEL_VASPSYM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_PRECACC_ALGONORM\")=" << AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_VASPRUNXML_STRESS\")=" << AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS\")=" << AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_SKIP_ARUNS_MAX\")=" << AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY\")=" << AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_SYMMETRIZE\")=" << AFLOWRC_DEFAULT_AEL_SYMMETRIZE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_FILE_PREFIX\")=\"" << AFLOWRC_DEFAULT_AEL_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_WRITE_FULL_RESULTS\")=" << AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AEL_DIRNAME_ARUN\")=" << AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN << endl;

    if(LDEBUG) oss << "// DEFAULTS AGL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_AEL_POISSON_RATIO\")=" << AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_NSTRUCTURES\")=" << AFLOWRC_DEFAULT_AGL_NSTRUCTURES << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_STRAIN_STEP\")=" << AFLOWRC_DEFAULT_AGL_STRAIN_STEP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS\")=" << AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_SKIP_ARUNS_MAX\")=" << AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX << endl;    
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_NTEMPERATURE\")=" << AFLOWRC_DEFAULT_AGL_NTEMPERATURE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_STEMPERATURE\")=" << AFLOWRC_DEFAULT_AGL_STEMPERATURE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_NPRESSURE\")=" << AFLOWRC_DEFAULT_AGL_NPRESSURE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_SPRESSURE\")=" << AFLOWRC_DEFAULT_AGL_SPRESSURE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_POISSON_RATIO\")=" << AFLOWRC_DEFAULT_AGL_POISSON_RATIO << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_IEOS\")=" << AFLOWRC_DEFAULT_AGL_IEOS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_IDEBYE\")=" << AFLOWRC_DEFAULT_AGL_IDEBYE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_FIT_TYPE\")=" << AFLOWRC_DEFAULT_AGL_FIT_TYPE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_CHECK_EV_CONCAVITY\")=" << AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_CHECK_EV_MIN\")=" << AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_HUGONIOT_CALC\")=" << AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_HUGONIOT_EXTRAPOLATE\")=" << AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE\")=" << AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE << endl;    
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_FILE_PREFIX\")=\"" << AFLOWRC_DEFAULT_AGL_FILE_PREFIX << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_WRITE_FULL_RESULTS\")=" << AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_DIRNAME_ARUN\")=" << AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_WRITE_GIBBS_INPUT\")=" << AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_AGL_PLOT_RESULTS\")=" << AFLOWRC_DEFAULT_AGL_PLOT_RESULTS << endl;    

    //SD20220323 - QCA START
    if(LDEBUG) oss << "// DEFAULTS QCA " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_MIN_SLEEP_SECONDS\")=" << AFLOWRC_DEFAULT_QCA_MIN_SLEEP_SECONDS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_MAX_NUM_ATOMS\")=" << AFLOWRC_DEFAULT_QCA_MAX_NUM_ATOMS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS\")=" << AFLOWRC_DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_CV_CUTOFF\")=" << AFLOWRC_DEFAULT_QCA_CV_CUTOFF << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_CONC_NPTS\")=" << AFLOWRC_DEFAULT_QCA_CONC_NPTS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_NPTS\")=" << AFLOWRC_DEFAULT_QCA_TEMP_NPTS << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_MIN\")=" << AFLOWRC_DEFAULT_QCA_TEMP_MIN << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_MAX\")=" << AFLOWRC_DEFAULT_QCA_TEMP_MAX << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_TEMP_MIN_LIMIT\")=" << AFLOWRC_DEFAULT_QCA_TEMP_MIN_LIMIT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_QCA_PRINT\")=" << AFLOWRC_DEFAULT_QCA_PRINT << endl;
    //SD20220323 - QCA END

    //RF20200413 START
    if(LDEBUG) oss << "// DEFAULTS CCE" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_OX_METHOD\")=\"" << DEFAULT_CCE_OX_METHOD << "\"" << "               // 1 - ELECTRONEGATIVITY_ALLEN, 2 - BADER" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION\")=\"" << DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION << "\"" << "               // 0.4 Ang tolerance between shortest and longest bonds when testing for multi-anion compound" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_OX_TOL\")=\"" << DEFAULT_CCE_OX_TOL << "\"" << "               // sum of oxidation states might not be exactly zero due to numerics" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_PEROX_CUTOFF\")=\"" << DEFAULT_CCE_PEROX_CUTOFF << "\"" << "               // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_SUPEROX_CUTOFF\")=\"" << DEFAULT_CCE_SUPEROX_CUTOFF << "\"" << "               // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF\")=\"" << DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF << "\"" << "               // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF\")=\"" << DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF << "\"" << "               // O-O bonds in the O2 molecule is about 1.21 Ang." << endl;
    //RF20200413 END

    //DX20200708 - START
    if(LDEBUG) oss << "// DEFAULTS XTALFINDER" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_MISFIT_MATCH\")=" << DEFAULT_XTALFINDER_MISFIT_MATCH << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_MISFIT_FAMILY\")=" << DEFAULT_XTALFINDER_MISFIT_FAMILY << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_SUPERCELL_METHOD\")=" << DEFAULT_XTALFINDER_SUPERCELL_METHOD << endl; //DX20201223
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING\")=" << DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_MATERIAL\")=" << DEFAULT_XTALFINDER_FILE_MATERIAL << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_STRUCTURE\"=" << DEFAULT_XTALFINDER_FILE_STRUCTURE << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_DUPLICATE\")=" << DEFAULT_XTALFINDER_FILE_DUPLICATE << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE\")=" << DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE\")=" << DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE\")=" << DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE << endl; //DX20201118
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE\")=" << DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE << endl; //DX20201118
    //DX20200708 - END

    //DX20200720 - START
    if(LDEBUG) oss << "// DEFAULTS ANRL" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL\")=" << DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL << endl;
    //DX20200708 - END

    if(LDEBUG) oss << "// DEFAULT CORE" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_BEEP\")=" << AFLOW_CORE_TEMPERATURE_BEEP << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_HALT\")=" << AFLOW_CORE_TEMPERATURE_HALT << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"AFLOW_CORE_TEMPERATURE_REFRESH\")=" << AFLOW_CORE_TEMPERATURE_REFRESH << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"SECONDS_SLEEP_VASP_COMPLETION\")=" << SECONDS_SLEEP_VASP_COMPLETION << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"SECONDS_SLEEP_VASP_MONITOR\")=" << SECONDS_SLEEP_VASP_MONITOR << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"SECONDS_STALE_OUTCAR\")=" << SECONDS_STALE_OUTCAR << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<unsigned long long int>(\"BYTES_MAX_VASP_OUT\")=" << BYTES_MAX_VASP_OUT << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"MEMORY_MAX_USAGE_RAM\")=" << MEMORY_MAX_USAGE_RAM << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<double>(\"MEMORY_MAX_USAGE_SWAP\")=" << MEMORY_MAX_USAGE_SWAP << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"FILE_VASP_MONITOR\")=" << FILE_VASP_MONITOR << endl; //CO20201111
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"INTEL_COMPILER_PATHS\")=" << INTEL_COMPILER_PATHS << endl; //CO20201111

    if(LDEBUG) oss << "// DEFAULT MACHINE DEPENDENT MPI" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_BETA_MPICH\")=\"" << MPI_OPTIONS_DUKE_BETA_MPICH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_BETA_MPICH\")=\"" << MPI_COMMAND_DUKE_BETA_MPICH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_BETA_MPICH\")=\"" << MPI_BINARY_DIR_DUKE_BETA_MPICH << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_BETA_OPENMPI\")=\"" << MPI_OPTIONS_DUKE_BETA_OPENMPI << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_BETA_OPENMPI\")=\"" << MPI_COMMAND_DUKE_BETA_OPENMPI << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_BETA_OPENMPI\")=\"" << MPI_BINARY_DIR_DUKE_BETA_OPENMPI << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_MATERIALS\")=\"" << MPI_OPTIONS_DUKE_MATERIALS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_MATERIALS\")=\"" << MPI_COMMAND_DUKE_MATERIALS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_MATERIALS\")=\"" << MPI_BINARY_DIR_DUKE_MATERIALS << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_AFLOWLIB\")=\"" << MPI_OPTIONS_DUKE_AFLOWLIB << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_AFLOWLIB\")=\"" << MPI_COMMAND_DUKE_AFLOWLIB << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_AFLOWLIB\")=\"" << MPI_BINARY_DIR_DUKE_AFLOWLIB << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_HABANA\")=\"" << MPI_OPTIONS_DUKE_HABANA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_HABANA\")=\"" << MPI_COMMAND_DUKE_HABANA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_HABANA\")=\"" << MPI_BINARY_DIR_DUKE_HABANA << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_QRATS_MPICH\")=\"" << MPI_OPTIONS_DUKE_QRATS_MPICH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_QRATS_MPICH\")=\"" << MPI_COMMAND_DUKE_QRATS_MPICH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_QRATS_MPICH\")=\"" << MPI_BINARY_DIR_DUKE_QRATS_MPICH << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_QFLOW_OPENMPI\")=\"" << MPI_OPTIONS_DUKE_QFLOW_OPENMPI << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_QFLOW_OPENMPI\")=\"" << MPI_COMMAND_DUKE_QFLOW_OPENMPI << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI\")=\"" << MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI << "\"" << endl;

    //CO20201220 X START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_X_X\")=\"" << MPI_OPTIONS_DUKE_X_X << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_X_X\")=\"" << MPI_COMMAND_DUKE_X_X << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_X_X\")=\"" << MPI_BINARY_DIR_DUKE_X_X << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_X_CRAY\")=\"" << MPI_OPTIONS_DUKE_X_CRAY << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_X_CRAY\")=\"" << MPI_COMMAND_DUKE_X_CRAY << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_X_CRAY\")=\"" << MPI_BINARY_DIR_DUKE_X_CRAY << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_X_OLDCRAY\")=\"" << MPI_OPTIONS_DUKE_X_OLDCRAY << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_X_OLDCRAY\")=\"" << MPI_COMMAND_DUKE_X_OLDCRAY << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_X_OLDCRAY\")=\"" << MPI_BINARY_DIR_DUKE_X_OLDCRAY << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_DUKE_X_SMB\")=\"" << MPI_OPTIONS_DUKE_X_SMB << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_DUKE_X_SMB\")=\"" << MPI_COMMAND_DUKE_X_SMB << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_DUKE_X_SMB\")=\"" << MPI_BINARY_DIR_DUKE_X_SMB << "\"" << endl;
    //CO20201220 X STOP
    
    //CO20220818 JHU_ROCKFISH START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_JHU_ROCKFISH\")=\"" << MPI_OPTIONS_JHU_ROCKFISH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_JHU_ROCKFISH\")=\"" << MPI_COMMAND_JHU_ROCKFISH << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_JHU_ROCKFISH\")=\"" << MPI_BINARY_DIR_JHU_ROCKFISH << "\"" << endl;
    //CO20220818 JHU_ROCKFISH STOP

    //DX20190509 - MACHINE001 - START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE001\")=\"" << MPI_OPTIONS_MACHINE001 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE001\")=\"" << MPI_COMMAND_MACHINE001 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE001\")=\"" << MPI_BINARY_DIR_MACHINE001 << "\"" << endl;
    //DX20190509 - MACHINE001 - END

    //DX20190509 - MACHINE002 - START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE002\")=\"" << MPI_OPTIONS_MACHINE002 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE002\")=\"" << MPI_COMMAND_MACHINE002 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE002\")=\"" << MPI_BINARY_DIR_MACHINE002 << "\"" << endl;
    //DX20190509 - MACHINE002 - END

    //DX20201005 - MACHINE003 - START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE003\")=\"" << MPI_OPTIONS_MACHINE003 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE003\")=\"" << MPI_COMMAND_MACHINE003 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE003\")=\"" << MPI_BINARY_DIR_MACHINE003 << "\"" << endl;
    //DX20201005 - MACHINE003 - END

    //DX20211011 - MACHINE004 - START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE004\")=\"" << MPI_OPTIONS_MACHINE004 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE004\")=\"" << MPI_COMMAND_MACHINE004 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE004\")=\"" << MPI_BINARY_DIR_MACHINE004 << "\"" << endl;
    //DX2021101 - MACHINE004 - END

    //DX20190107 - CMU EULER - START
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_CMU_EULER\")=\"" << MPI_OPTIONS_CMU_EULER << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_CMU_EULER\")=\"" << MPI_COMMAND_CMU_EULER << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_CMU_EULER\")=\"" << MPI_BINARY_DIR_CMU_EULER << "\"" << endl;
    //DX20190107 - CMU EULER - END

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_EOS\")=\"" << MPI_OPTIONS_MPCDF_EOS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_EOS\")=\"" << MPI_COMMAND_MPCDF_EOS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MPCDF_EOS\")=\"" << MPI_NCPUS_MPCDF_EOS << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype(\"MPI_HYPERTHREADING_MPCDF_EOS\")=\"" << MPI_HYPERTHREADING_MPCDF_EOS << "\"" << "            // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_EOS\")=\"" << MPI_BINARY_DIR_MPCDF_EOS << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_DRACO\")=\"" << MPI_OPTIONS_MPCDF_DRACO << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_DRACO\")=\"" << MPI_COMMAND_MPCDF_DRACO << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MPCDF_DRACO\")=\"" << MPI_NCPUS_MPCDF_DRACO << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype(\"MPI_HYPERTHREADING_MPCDF_DRACO\")=\"" << MPI_HYPERTHREADING_MPCDF_DRACO << "\"" << "            // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_DRACO\")=\"" << MPI_BINARY_DIR_MPCDF_DRACO << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_COBRA\")=\"" << MPI_OPTIONS_MPCDF_COBRA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_COBRA\")=\"" << MPI_COMMAND_MPCDF_COBRA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MPCDF_COBRA\")=\"" << MPI_NCPUS_MPCDF_COBRA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype(\"MPI_HYPERTHREADING_MPCDF_COBRA\")=\"" << MPI_HYPERTHREADING_MPCDF_COBRA << "\"" << "            // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_COBRA\")=\"" << MPI_BINARY_DIR_MPCDF_COBRA << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MPCDF_HYDRA\")=\"" << MPI_OPTIONS_MPCDF_HYDRA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MPCDF_HYDRA\")=\"" << MPI_COMMAND_MPCDF_HYDRA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype<int>(\"MPI_NCPUS_MPCDF_HYDRA\")=\"" << MPI_NCPUS_MPCDF_HYDRA << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedutype(\"MPI_HYPERTHREADING_MPCDF_HYDRA\")=\"" << MPI_HYPERTHREADING_MPCDF_HYDRA << "\"" << "            // FALSE/OFF, IGNORE/NEGLECT, TRUE/ON " << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MPCDF_HYDRA\")=\"" << MPI_BINARY_DIR_MPCDF_HYDRA << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_FULTON_MARYLOU\")=\"" << MPI_OPTIONS_FULTON_MARYLOU << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_FULTON_MARYLOU\")=\"" << MPI_COMMAND_FULTON_MARYLOU << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_FULTON_MARYLOU\")=\"" << MPI_BINARY_DIR_FULTON_MARYLOU << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE1\")=\"" << MPI_OPTIONS_MACHINE1 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE1\")=\"" << MPI_COMMAND_MACHINE1 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE1\")=\"" << MPI_BINARY_DIR_MACHINE1 << "\"" << endl;

    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_OPTIONS_MACHINE2\")=\"" << MPI_OPTIONS_MACHINE2 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_COMMAND_MACHINE2\")=\"" << MPI_COMMAND_MACHINE2 << "\"" << endl;
    if(LDEBUG) oss << "XHOST.adefault.getattachedscheme(\"MPI_BINARY_DIR_MACHINE2\")=\"" << MPI_BINARY_DIR_MACHINE2 << "\"" << endl;

    //   if(LDEBUG) oss << "XHOST.adefault.content=" << XHOST.adefault.content_string << endl;

    if(LDEBUG) oss << "aflowrc::print_aflowrc: END" << endl;

    oss.flush();
    return FALSE;
  }
} // namespace aflowrc

// **************************************************************************
// **************************************************************************

#endif // _AFLOW_AFLOWRC_CPP_

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2023              *
// *                                                                        *
// **************************************************************************
