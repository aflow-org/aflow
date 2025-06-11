// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by SC 2017-2018

// ***************************************************************************
#ifndef _AFLOW_AFLOWRC_H_
#define _AFLOW_AFLOWRC_H_

#include <string>
#include <iosfwd>

#include "aflow_xhost.h"
// #include "aflow_init.h"

#define AFLOWRC_FILENAME_LOCAL   XHOST.home+"/.aflow.rc"
#define AFLOWRC_FILENAME_GLOBAL  "/etc/aflow.conf"

// DEFAULT DEFINITIONS
#define AFLOWRC_AFLOWRC std::string(AFLOW_VERSION)

// DEFAULT DEFINITIONS
#define AFLOWRC_DEFAULT_KZIP_BIN                        std::string("xz")
#define         DEFAULT_KZIP_BIN                        XHOST.adefault.getattachedscheme("DEFAULT_KZIP_BIN")
#define AFLOWRC_DEFAULT_KZIP_EXT                        std::string(".xz")
#define         DEFAULT_KZIP_EXT                        XHOST.adefault.getattachedscheme("DEFAULT_KZIP_EXT")
#define AFLOWRC_DEFAULT_TMPFS_DIRECTORIES               std::string("/tmp/,/run/shm/,/dev/shm/")
#define         DEFAULT_TMPFS_DIRECTORIES               XHOST.adefault.getattachedscheme("DEFAULT_TMPFS_DIRECTORIES")


//HE20220218 START
// DEFAULTS ENTRY LOADER
#define AFLOWRC_DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE      std::string("~/.aflow/aflowlib_alloy.db")
#define         DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE      XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_SERVER       std::string("aflowlib.duke.edu")
#define         DEFAULT_ENTRY_LOADER_AFLUX_SERVER       XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_AFLUX_SERVER")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_AFLUX_PATH         std::string("/API/aflux/")
#define         DEFAULT_ENTRY_LOADER_AFLUX_PATH         XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_AFLUX_PATH")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_SERVER     std::string("aflowlib.duke.edu")
#define         DEFAULT_ENTRY_LOADER_RESTAPI_SERVER     XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_RESTAPI_SERVER")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_RESTAPI_PATH       std::string("/AFLOWDATA/")
#define         DEFAULT_ENTRY_LOADER_RESTAPI_PATH       XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_RESTAPI_PATH")
#define AFLOWRC_DEFAULT_ENTRY_LOADER_FS_PATH            std::string("/common/")
#define         DEFAULT_ENTRY_LOADER_FS_PATH            XHOST.adefault.getattachedscheme("DEFAULT_ENTRY_LOADER_FS_PATH")
//HE20220218 STOP

//ME20191001 START
// DEFAULTS AFLOW DATABASE
#define AFLOWRC_DEFAULT_AFLOW_DB_FILE                   std::string("/var/cache/aflow_data/AFLOWDB/aflowlib.db")
#define         DEFAULT_AFLOW_DB_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_FILE")
#define AFLOWRC_DEFAULT_AFLOW_DB_STATS_FILE             std::string("/var/cache/aflow_data/AFLOWDB/aflowlib.json")
#define         DEFAULT_AFLOW_DB_STATS_FILE             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_STATS_FILE")
#define AFLOWRC_DEFAULT_AFLOW_DB_DATA_PATH              std::string("/common/AFLOW/LIBS/")
#define         DEFAULT_AFLOW_DB_DATA_PATH              XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_DATA_PATH")
#define AFLOWRC_DEFAULT_AFLOW_DB_LOCK_FILE              std::string("/var/cache/aflow_data/AFLOWDB/ADB_Idle.lock")
#define         DEFAULT_AFLOW_DB_LOCK_FILE              XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_DB_LOCK_FILE")
#define AFLOWRC_DEFAULT_AFLOW_DB_STALE_THRESHOLD        3*3600
#define         DEFAULT_AFLOW_DB_STALE_THRESHOLD        XHOST.adefault.getattachedutype<long int>("DEFAULT_AFLOW_DB_STALE_THRESHOLD")
//ME20191001 STOP

// FILENAMES FOR AFLOW.ORG ANALYSIS
#define AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_OUT         std::string("aflowlib.out")
#define         DEFAULT_FILE_AFLOWLIB_ENTRY_OUT         XHOST.adefault.getattachedscheme("DEFAULT_FILE_AFLOWLIB_ENTRY_OUT")
#define AFLOWRC_DEFAULT_FILE_AFLOWLIB_ENTRY_JSON        std::string("aflowlib.json")
#define         DEFAULT_FILE_AFLOWLIB_ENTRY_JSON        XHOST.adefault.getattachedscheme("DEFAULT_FILE_AFLOWLIB_ENTRY_JSON")
#define AFLOWRC_DEFAULT_FILE_EDATA_ORIG_OUT             std::string("edata.orig.out")
#define         DEFAULT_FILE_EDATA_ORIG_OUT             XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_ORIG_OUT")
#define AFLOWRC_DEFAULT_FILE_EDATA_RELAX_OUT            std::string("edata.relax.out")
#define         DEFAULT_FILE_EDATA_RELAX_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_RELAX_OUT")
#define AFLOWRC_DEFAULT_FILE_EDATA_BANDS_OUT            std::string("edata.bands.out")
#define         DEFAULT_FILE_EDATA_BANDS_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_BANDS_OUT")
#define AFLOWRC_DEFAULT_FILE_DATA_ORIG_OUT              std::string("data.orig.out")
#define         DEFAULT_FILE_DATA_ORIG_OUT              XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_ORIG_OUT")
#define AFLOWRC_DEFAULT_FILE_DATA_RELAX_OUT             std::string("data.relax.out")
#define         DEFAULT_FILE_DATA_RELAX_OUT             XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_RELAX_OUT")
#define AFLOWRC_DEFAULT_FILE_DATA_BANDS_OUT             std::string("data.bands.out")
#define         DEFAULT_FILE_DATA_BANDS_OUT             XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_BANDS_OUT")
#define AFLOWRC_DEFAULT_FILE_EDATA_ORIG_JSON            std::string("edata.orig.json")
#define         DEFAULT_FILE_EDATA_ORIG_JSON            XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_ORIG_JSON")
#define AFLOWRC_DEFAULT_FILE_EDATA_RELAX_JSON           std::string("edata.relax.json")
#define         DEFAULT_FILE_EDATA_RELAX_JSON           XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_RELAX_JSON")
#define AFLOWRC_DEFAULT_FILE_EDATA_BANDS_JSON           std::string("edata.bands.json")
#define         DEFAULT_FILE_EDATA_BANDS_JSON           XHOST.adefault.getattachedscheme("DEFAULT_FILE_EDATA_BANDS_JSON")
#define AFLOWRC_DEFAULT_FILE_DATA_ORIG_JSON             std::string("data.orig.json")
#define         DEFAULT_FILE_DATA_ORIG_JSON             XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_ORIG_JSON")
#define AFLOWRC_DEFAULT_FILE_DATA_RELAX_JSON            std::string("data.relax.json")
#define         DEFAULT_FILE_DATA_RELAX_JSON            XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_RELAX_JSON")
#define AFLOWRC_DEFAULT_FILE_DATA_BANDS_JSON            std::string("data.bands.json")
#define         DEFAULT_FILE_DATA_BANDS_JSON            XHOST.adefault.getattachedscheme("DEFAULT_FILE_DATA_BANDS_JSON")
#define AFLOWRC_DEFAULT_FILE_TIME_OUT                   std::string("time")
#define         DEFAULT_FILE_TIME_OUT                   XHOST.adefault.getattachedscheme("DEFAULT_FILE_TIME_OUT")
#define AFLOWRC_DEFAULT_FILE_SPACEGROUP1_OUT            std::string("SpaceGroup")
#define         DEFAULT_FILE_SPACEGROUP1_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_SPACEGROUP1_OUT")
#define AFLOWRC_DEFAULT_FILE_SPACEGROUP2_OUT            std::string("SpaceGroup2")
#define         DEFAULT_FILE_SPACEGROUP2_OUT            XHOST.adefault.getattachedscheme("DEFAULT_FILE_SPACEGROUP2_OUT")
#define AFLOWRC_DEFAULT_FILE_VOLDISTPARAMS_OUT          std::string("VOLDISTParams")
#define         DEFAULT_FILE_VOLDISTPARAMS_OUT          XHOST.adefault.getattachedscheme("DEFAULT_FILE_VOLDISTPARAMS_OUT")
#define AFLOWRC_DEFAULT_FILE_VOLDISTEVOLUTION_OUT       std::string("VOLDISTEvolution")
#define         DEFAULT_FILE_VOLDISTEVOLUTION_OUT       XHOST.adefault.getattachedscheme("DEFAULT_FILE_VOLDISTEVOLUTION_OUT")

// FILENAMES FOR AFLOW OPERATION
#define AFLOWRC_DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT  std::string("aflow.pseudopotential_auid.out")
#define         DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT  XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_OUT             std::string("aflow.prescript.out")
#define         DEFAULT_AFLOW_PRESCRIPT_OUT             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PRESCRIPT_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PRESCRIPT_COMMAND         std::string("aflow.prescript.command")
#define         DEFAULT_AFLOW_PRESCRIPT_COMMAND         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PRESCRIPT_COMMAND")
#define AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_OUT            std::string("aflow.postscript.out")
#define         DEFAULT_AFLOW_POSTSCRIPT_OUT            XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_POSTSCRIPT_OUT")
#define AFLOWRC_DEFAULT_AFLOW_POSTSCRIPT_COMMAND        std::string("aflow.postscript.command")
#define         DEFAULT_AFLOW_POSTSCRIPT_COMMAND        XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_POSTSCRIPT_COMMAND")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_OUT                std::string("aflow.pgroup.out")
#define         DEFAULT_AFLOW_PGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_JSON               std::string("aflow.pgroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_PGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_OUT           std::string("aflow.pgroup_xtal.out")
#define         DEFAULT_AFLOW_PGROUP_XTAL_OUT           XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_XTAL_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUP_XTAL_JSON          std::string("aflow.pgroup_xtal.json") //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_PGROUP_XTAL_JSON          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUP_XTAL_JSON")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT     std::string("aflow.pgroupk_Patterson.out") //DX20200129
#define         DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT     XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT") //DX20200129
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON    std::string("aflow.pgroupk_Patterson.json") //DX20200129
#define         DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON    XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON") //DX20200129
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_OUT               std::string("aflow.pgroupk.out")
#define         DEFAULT_AFLOW_PGROUPK_OUT               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_JSON              std::string("aflow.pgroupk.json")     //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_PGROUPK_JSON              XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_JSON")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_OUT          std::string("aflow.pgroupk_xtal.out") //DX20171205 - Added pgroupk_xtal
#define         DEFAULT_AFLOW_PGROUPK_XTAL_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_XTAL_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PGROUPK_XTAL_JSON         std::string("aflow.pgroupk_xtal.json")//DX20170802 - Add JSON //DX20171205 - Added pgroupk_xtal
#define         DEFAULT_AFLOW_PGROUPK_XTAL_JSON         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PGROUPK_XTAL_JSON")
#define AFLOWRC_DEFAULT_AFLOW_FGROUP_OUT                std::string("aflow.fgroup.out")
#define         DEFAULT_AFLOW_FGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FGROUP_JSON               std::string("aflow.fgroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_FGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_SGROUP_OUT                std::string("aflow.sgroup.out")
#define         DEFAULT_AFLOW_SGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_SGROUP_JSON               std::string("aflow.sgroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_SGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_AGROUP_OUT                std::string("aflow.agroup.out")
#define         DEFAULT_AFLOW_AGROUP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_AGROUP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_AGROUP_JSON               std::string("aflow.agroup.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_AGROUP_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_AGROUP_JSON")
#define AFLOWRC_DEFAULT_AFLOW_IATOMS_OUT                std::string("aflow.iatoms.out")
#define         DEFAULT_AFLOW_IATOMS_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IATOMS_OUT")
#define AFLOWRC_DEFAULT_AFLOW_IATOMS_JSON               std::string("aflow.iatoms.json")      //DX20170802 - Add JSON
#define         DEFAULT_AFLOW_IATOMS_JSON               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IATOMS_JSON")
#define AFLOWRC_DEFAULT_AFLOW_ICAGES_OUT                std::string("aflow.icages.out")
#define         DEFAULT_AFLOW_ICAGES_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_ICAGES_OUT")
#define AFLOWRC_DEFAULT_AFLOW_SURFACE_OUT               std::string("aflow.surface.out")
#define         DEFAULT_AFLOW_SURFACE_OUT               XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_SURFACE_OUT")
#define AFLOWRC_DEFAULT_AFLOW_QMVASP_OUT                std::string("aflow.qmvasp.out")
#define         DEFAULT_AFLOW_QMVASP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_QMVASP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_ERVASP_OUT                std::string("aflow.error.out")
#define         DEFAULT_AFLOW_ERVASP_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_ERVASP_OUT")
#define AFLOWRC_DEFAULT_AFLOW_IMMISCIBILITY_OUT         std::string("aflow.immiscibility.out")
#define         DEFAULT_AFLOW_IMMISCIBILITY_OUT         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_IMMISCIBILITY_OUT")
#define AFLOWRC_DEFAULT_AFLOW_MEMORY_OUT                std::string("aflow.memory.out")
#define         DEFAULT_AFLOW_MEMORY_OUT                XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_MEMORY_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_INPUT_OUT          std::string("aflow.frozsl_input.out")
#define         DEFAULT_AFLOW_FROZSL_INPUT_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_INPUT_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_POSCAR_OUT         std::string("aflow.frozsl_poscar.out")
#define         DEFAULT_AFLOW_FROZSL_POSCAR_OUT         XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_POSCAR_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_MODES_OUT          std::string("aflow.frozsl_energies.out")
#define         DEFAULT_AFLOW_FROZSL_MODES_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_MODES_OUT")
#define AFLOWRC_DEFAULT_AFLOW_FROZSL_EIGEN_OUT          std::string("aflow.frozsl_eigen.out")
#define         DEFAULT_AFLOW_FROZSL_EIGEN_OUT          XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_FROZSL_EIGEN_OUT")
#define AFLOWRC_DEFAULT_AFLOW_END_OUT                   std::string("aflow.end.out")
#define         DEFAULT_AFLOW_END_OUT                   XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_END_OUT")
#define AFLOWRC_DEFAULT_AFLOW_PLASMONICS_FILE           std::string("aflow.plasmonics_eps")
#define         DEFAULT_AFLOW_PLASMONICS_FILE           XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PLASMONICS_FILE")

// GENERIC MPI   // DONE
#define AFLOWRC_MPI_START_DEFAULT                       std::string("")
#define         MPI_START_DEFAULT                       XHOST.adefault.getattachedscheme("MPI_START_DEFAULT")
#define AFLOWRC_MPI_STOP_DEFAULT                        std::string("")
#define         MPI_STOP_DEFAULT                        XHOST.adefault.getattachedscheme("MPI_STOP_DEFAULT")
#define AFLOWRC_MPI_COMMAND_DEFAULT                     std::string("mpirun -np")
#define         MPI_COMMAND_DEFAULT                     XHOST.adefault.getattachedscheme("MPI_COMMAND_DEFAULT")
#define AFLOWRC_MPI_NCPUS_DEFAULT                       4
#define         MPI_NCPUS_DEFAULT                       XHOST.adefault.getattachedutype<int>("MPI_NCPUS_DEFAULT")
#define AFLOWRC_MPI_NCPUS_MAX                           init::GetCPUCores() //16  //CO20180124
#define         MPI_NCPUS_MAX                           XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MAX")

// BINARY    // DONE
#define AFLOWRC_DEFAULT_VASP_GAMMA_BIN                  std::string("vasp_gam")
#define         DEFAULT_VASP_GAMMA_BIN                  XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_BIN")
#define AFLOWRC_DEFAULT_VASP_GAMMA_MPI_BIN              std::string("vasp_gam")
#define         DEFAULT_VASP_GAMMA_MPI_BIN              XHOST.adefault.getattachedscheme("DEFAULT_VASP_GAMMA_MPI_BIN")
#define AFLOWRC_DEFAULT_VASP_BIN                        std::string("vasp_std")
#define         DEFAULT_VASP_BIN                        XHOST.adefault.getattachedscheme("DEFAULT_VASP_BIN")
#define AFLOWRC_DEFAULT_VASP_MPI_BIN                    std::string("vasp_std")
#define         DEFAULT_VASP_MPI_BIN                    XHOST.adefault.getattachedscheme("DEFAULT_VASP_MPI_BIN")
#define AFLOWRC_DEFAULT_VASP5_BIN                       std::string("vasp_std")
#define         DEFAULT_VASP5_BIN                       XHOST.adefault.getattachedscheme("DEFAULT_VASP5_BIN")
#define AFLOWRC_DEFAULT_VASP5_MPI_BIN                   std::string("vasp_std")
#define         DEFAULT_VASP5_MPI_BIN                   XHOST.adefault.getattachedscheme("DEFAULT_VASP5_MPI_BIN")
//aims
#define AFLOWRC_DEFAULT_AIMS_BIN                        std::string("aims")
#define         DEFAULT_AIMS_BIN                        XHOST.adefault.getattachedscheme("DEFAULT_AIMS_BIN")

// POTCARS // DONE
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIRECTORIES               std::string("/common/VASP,/common/AFLOW/VASP,/home/aflow/common/AFLOW/VASP,/fslhome/fslcollab8/group/VASP,/fslhome/glh43/src/,/share/home/00470/tg457283/common/AFLOW/VASP/,/share/home/00457/tg457357/common/AFLOW/VASP/,/home/mehl/bin/AFLOW/VASP/,~/common/VASP/,~/common/AFLOW/VASP/,/home/aflow/common/VASP/,/nics/a/proj/aflow/common/AFLOW/VASP/,/home/users/aflow/common/VASP,/share/apps/AFLOW3/VASP,/share/apps/vasp/PP,/projects/kyang-group/common/VASP,/home/Tools/src/vasp/,/somewhere/")  // first is default, tokenized with "," //DX20190107 - added CMU path
#define         DEFAULT_VASP_POTCAR_DIRECTORIES               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIRECTORIES")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE                      std::string("current")
#define         DEFAULT_VASP_POTCAR_DATE                      XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_SUFFIX                    std::string("/POTCAR")
#define         DEFAULT_VASP_POTCAR_SUFFIX                    XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_SUFFIX")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_LDA              std::string("01Apr2000")   // when no date is given for pot_LDA
#define         DEFAULT_VASP_POTCAR_DATE_POT_LDA              XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DATE_POT_GGA              std::string("01Apr2000")   // when no date is given for pot_GGA
#define         DEFAULT_VASP_POTCAR_DATE_POT_GGA              XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DATE_POT_GGA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_LDA               std::string("pot_LDA")
#define         DEFAULT_VASP_POTCAR_DIR_POT_LDA               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_GGA               std::string("pot_GGA")
#define         DEFAULT_VASP_POTCAR_DIR_POT_GGA               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_GGA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POT_PBE               std::string("pot_PBE")
#define         DEFAULT_VASP_POTCAR_DIR_POT_PBE               XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POT_PBE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA            std::string("potpaw_LDA")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA            std::string("potpaw_GGA")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE            std::string("potpaw_PBE")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE            XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN        std::string("potpaw_LDA.54")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN        XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN")
#define AFLOWRC_DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN        std::string("potpaw_PBE.54")
#define         DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN        XHOST.adefault.getattachedscheme("DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN")

// KPOINTS/DOS // DONE
#define AFLOWRC_DEFAULT_BANDS_GRID                            20
#define         DEFAULT_BANDS_GRID                            XHOST.adefault.getattachedutype<int>("DEFAULT_BANDS_GRID") 
#define AFLOWRC_DEFAULT_BANDS_LATTICE                         std::string("AUTO")
#define         DEFAULT_BANDS_LATTICE                         XHOST.adefault.getattachedscheme("DEFAULT_BANDS_LATTICE")
#define AFLOWRC_DEFAULT_KSCHEME                               std::string("AUTO")
#define         DEFAULT_KSCHEME                               XHOST.adefault.getattachedscheme("DEFAULT_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA                                 6000
#define         DEFAULT_KPPRA                                 XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA")
#define AFLOWRC_DEFAULT_KPPRA_STATIC                          10000
#define         DEFAULT_KPPRA_STATIC                          XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_STATIC")
#define AFLOWRC_DEFAULT_STATIC_KSCHEME                        std::string("AUTO")
#define         DEFAULT_STATIC_KSCHEME                        XHOST.adefault.getattachedscheme("DEFAULT_STATIC_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA_DIELECTRIC                      12500 
#define         DEFAULT_KPPRA_DIELECTRIC                      XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_DIELECTRIC")
#define AFLOWRC_DEFAULT_DIELECTRIC_KSCHEME                    std::string("AUTO")
#define         DEFAULT_DIELECTRIC_KSCHEME                    XHOST.adefault.getattachedscheme("DEFAULT_DIELECTRIC_KSCHEME")
#define AFLOWRC_DEFAULT_KPPRA_ICSD                            8000
#define         DEFAULT_KPPRA_ICSD                            XHOST.adefault.getattachedutype<int>("DEFAULT_KPPRA_ICSD")
#define AFLOWRC_DEFAULT_UNARY_BANDS_GRID                      128
#define         DEFAULT_UNARY_BANDS_GRID                      XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_BANDS_GRID")
#define AFLOWRC_DEFAULT_UNARY_KPPRA                           8000
#define         DEFAULT_UNARY_KPPRA                           XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA")
#define AFLOWRC_DEFAULT_UNARY_KPPRA_STATIC                    8000
#define         DEFAULT_UNARY_KPPRA_STATIC                    XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA_STATIC")
#define AFLOWRC_DEFAULT_UNARY_KPPRA_DIELECTRIC                8000
#define         DEFAULT_UNARY_KPPRA_DIELECTRIC                XHOST.adefault.getattachedutype<int>("DEFAULT_UNARY_KPPRA_DIELECTRIC")
#define AFLOWRC_DEFAULT_PHONONS_KSCHEME                       std::string("G")
#define         DEFAULT_PHONONS_KSCHEME                       XHOST.adefault.getattachedscheme("DEFAULT_PHONONS_KSCHEME")
#define AFLOWRC_DEFAULT_PHONONS_KPPRA                         3000  //CO20181226  //ME20190205 - 8000 uses too much memory e.g. for NaF - 2000 appears sufficient; ME20200108 - 3000 minimum for metals
#define         DEFAULT_PHONONS_KPPRA                         XHOST.adefault.getattachedutype<int>("DEFAULT_PHONONS_KPPRA")
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
#define AFLOWRC_DEFAULT_VASP_OUT                                      std::string("vasp.out")
#define         DEFAULT_VASP_OUT                                      XHOST.adefault.getattachedscheme("DEFAULT_VASP_OUT")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_INCAR                           std::string("./INCAR")
#define         DEFAULT_VASP_EXTERNAL_INCAR                           XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_INCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_POSCAR                          std::string("./POSCAR")
#define         DEFAULT_VASP_EXTERNAL_POSCAR                          XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POSCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_POTCAR                          std::string("./POTCAR")
#define         DEFAULT_VASP_EXTERNAL_POTCAR                          XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_POTCAR")
#define AFLOWRC_DEFAULT_VASP_EXTERNAL_KPOINTS                         std::string("./KPOINTS")
#define         DEFAULT_VASP_EXTERNAL_KPOINTS                         XHOST.adefault.getattachedscheme("DEFAULT_VASP_EXTERNAL_KPOINTS")
#define AFLOWRC_DEFAULT_AIMS_EXTERNAL_CONTROL                         std::string("./control.in")
#define         DEFAULT_AIMS_EXTERNAL_CONTROL                         XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_CONTROL")
#define AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM                            std::string("./geometry.in")
#define         DEFAULT_AIMS_EXTERNAL_GEOM                            XHOST.adefault.getattachedscheme("DEFAULT_AIMS_EXTERNAL_GEOM")
#define AFLOWRC_DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE                     std::string("potpaw_PBE")
#define         DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE                     XHOST.adefault.getattachedscheme("DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME           std::string("ENERGY")
#define         DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME           XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT                 2 //CO20181226
#define         DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT                 XHOST.adefault.getattachedutype<int>("DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT") //CO20181226
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME                 std::string("ACCURATE")
#define         DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME                 std::string("NORMAL")
#define         DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME              std::string("NONE")
#define         DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME              XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME                 std::string("0")
#define         DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME                 XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME                 std::string("DEFAULT")
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
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME                std::string("DEFAULT")
#define         DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME                XHOST.adefault.getattachedscheme("DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME")
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SYM                         true
#define         DEFAULT_VASP_FORCE_OPTION_SYM                         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SYM") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN                        true
#define         DEFAULT_VASP_FORCE_OPTION_SPIN                        XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1         false
#define         DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2         true  //ME20190308 - remove spin after two relaxations if zero
#define         DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER                       false
#define         DEFAULT_VASP_FORCE_OPTION_BADER                       XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_BADER") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_BADER_STATIC                true
#define         DEFAULT_VASP_FORCE_OPTION_BADER_STATIC                XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_BADER_STATIC") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_ELF                         false
#define         DEFAULT_VASP_FORCE_OPTION_ELF                         XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_ELF") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM                 false   // true
#define         DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM                 XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_WAVECAR                     false
#define         DEFAULT_VASP_FORCE_OPTION_WAVECAR                     XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_WAVECAR") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_CHGCAR                      true
#define         DEFAULT_VASP_FORCE_OPTION_CHGCAR                      XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_CHGCAR") 
#define AFLOWRC_DEFAULT_VASP_FORCE_OPTION_LSCOUPLING                  false
#define         DEFAULT_VASP_FORCE_OPTION_LSCOUPLING                  XHOST.adefault.getattachedutype<bool>("DEFAULT_VASP_FORCE_OPTION_LSCOUPLING") 

// AFLOW_LIBRARY AFLOW_PROJECT // DONE
#define AFLOWRC_DEFAULT_AFLOW_LIBRARY_DIRECTORIES             std::string("/common/AFLOW/LIBS/,/home/aflow/common/AFLOW/LIBS/,/fslhome/glh43/src/,/usr/local/bin/,/fslhome/fslcollab8/group/bin/,/home/auro/work/AFLOW3/,~/common/AFLOW/LIBS/,./,/nics/a/proj/aflow/common/AFLOW/LIBS/,/home/users/aflow/common/AFLOW/LIBS,/home/junkai/PROTO_DATABASE/,/projects/kyang-group/common/LIBS,/somewhere/")  // first is default, tokenized with ","
#define         DEFAULT_AFLOW_LIBRARY_DIRECTORIES             XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_LIBRARY_DIRECTORIES")
#define AFLOWRC_DEFAULT_AFLOW_PROJECTS_DIRECTORIES            std::string("/common/AUID,/common/ICSD,/common/LIB0,/common/LIB1,/common/LIB2,/common/LIB3,/common/LIB4,/common/LIB5,/common/LIB6,/common/LIB7,/common/LIB8,/common/LIB9")  // first is default, tokenized with ","
#define         DEFAULT_AFLOW_PROJECTS_DIRECTORIES            XHOST.adefault.getattachedscheme("DEFAULT_AFLOW_PROJECTS_DIRECTORIES")
#define AFLOWRC_DEFAULT_AFLOWDATA_WEB_DIRECTORY               std::string("/www/AFLOWDATA")  //CO+ME20200731
#define         DEFAULT_AFLOWDATA_WEB_DIRECTORY               XHOST.adefault.getattachedscheme("DEFAULT_AFLOWDATA_WEB_DIRECTORY") //CO+ME20200731

// PLATON/FINDSYM // DONE
#define AFLOWRC_DEFAULT_PLATON_P_EQUAL                        false
#define         DEFAULT_PLATON_P_EQUAL                        XHOST.adefault.getattachedutype<bool>("DEFAULT_PLATON_P_EQUAL") 
#define AFLOWRC_DEFAULT_PLATON_P_EXACT                        false
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
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT                      std::string("Helvetica")              // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT                      XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD                 std::string("Helvetica-Bold")         // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_BOLD                 XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_ITALICS              std::string("Helvetica-Oblique")      // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_ITALICS              XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS         std::string("Helvetica-BoldOblique")  // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS         XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT                      std::string("Arial")                  // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT                      XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD                 std::string("Arial_Bold")             // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_BOLD                 XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_ITALICS              std::string("Arial_Italic")           // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_ITALICS              XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS         std::string("Arial_BoldItalic")       // we do not have time to waste with ttf
#define         DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS         XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT                    std::string("Symbol")
#define         DEFAULT_GNUPLOT_GREEK_FONT                    XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD               std::string("Symbol-Bold")
#define         DEFAULT_GNUPLOT_GREEK_FONT_BOLD               XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_BOLD")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_ITALICS            std::string("Symbol-Oblique")
#define         DEFAULT_GNUPLOT_GREEK_FONT_ITALICS            XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_ITALICS")
#define AFLOWRC_DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS       std::string("Symbol-BoldOblique")
#define         DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS       XHOST.adefault.getattachedscheme("DEFAULT_GNUPLOT_GREEK_FONT_BOLD_ITALICS")

// DEFAULT CHULL
#define AFLOWRC_DEFAULT_CHULL_ALLOWED_DFT_TYPES                           std::string("PAW_PBE,PAW_PBE_KIN")
#define         DEFAULT_CHULL_ALLOWED_DFT_TYPES                           XHOST.adefault.getattachedscheme("DEFAULT_CHULL_ALLOWED_DFT_TYPES")
#define AFLOWRC_DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES                false
#define         DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES                XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES")
#define AFLOWRC_DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES                    200
#define         DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES                    XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES")
#define AFLOWRC_DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS                    true
#define         DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS                    XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS")
#define AFLOWRC_DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES   75
#define         DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES   XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES")
#define AFLOWRC_DEFAULT_CHULL_OUTLIER_MULTIPLIER                          3.25
#define         DEFAULT_CHULL_OUTLIER_MULTIPLIER                          XHOST.adefault.getattachedutype<double>("DEFAULT_CHULL_OUTLIER_MULTIPLIER")
#define AFLOWRC_DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED                  true
#define         DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED                  XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED")
#define AFLOWRC_DEFAULT_CHULL_LATEX_BANNER                                2 //CO20180827
#define         DEFAULT_CHULL_LATEX_BANNER                                XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_BANNER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN                      false
#define         DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN")
#define AFLOWRC_DEFAULT_CHULL_LATEX_STOICH_HEADER                         false
#define         DEFAULT_CHULL_LATEX_STOICH_HEADER                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_STOICH_HEADER")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNARIES                          false
#define         DEFAULT_CHULL_LATEX_PLOT_UNARIES                          XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_PLOT_UNARIES")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_OFF_HULL                         -1
#define         DEFAULT_CHULL_LATEX_PLOT_OFF_HULL                         XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_PLOT_OFF_HULL")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_UNSTABLE                         false
#define         DEFAULT_CHULL_LATEX_PLOT_UNSTABLE                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_PLOT_UNSTABLE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_SCHEME                         std::string("")
#define         DEFAULT_CHULL_LATEX_FILTER_SCHEME                         XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_FILTER_SCHEME")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FILTER_VALUE                          0.0
#define         DEFAULT_CHULL_LATEX_FILTER_VALUE                          XHOST.adefault.getattachedutype<double>("DEFAULT_CHULL_LATEX_FILTER_VALUE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_BAR                             true
#define         DEFAULT_CHULL_LATEX_COLOR_BAR                             XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COLOR_BAR")
#define AFLOWRC_DEFAULT_CHULL_LATEX_HEAT_MAP                              true
#define         DEFAULT_CHULL_LATEX_HEAT_MAP                              XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_HEAT_MAP")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_GRADIENT                        true
#define         DEFAULT_CHULL_LATEX_COLOR_GRADIENT                        XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_COLOR_GRADIENT")
#define AFLOWRC_DEFAULT_CHULL_LATEX_COLOR_MAP                             std::string("")
#define         DEFAULT_CHULL_LATEX_COLOR_MAP                             XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_COLOR_MAP")
#define AFLOWRC_DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR                   std::string("")
#define         DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR                   XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR")
#define AFLOWRC_DEFAULT_CHULL_LATEX_REVERSE_AXIS                          false
#define         DEFAULT_CHULL_LATEX_REVERSE_AXIS                          XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_REVERSE_AXIS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW                false
#define         DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW                XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LINKS                                 1
#define         DEFAULT_CHULL_LATEX_LINKS                                 XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_LINKS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LABEL_NAME                            std::string("")
#define         DEFAULT_CHULL_LATEX_LABEL_NAME                            XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_LABEL_NAME")
#define AFLOWRC_DEFAULT_CHULL_LATEX_META_LABELS                           false
#define         DEFAULT_CHULL_LATEX_META_LABELS                           XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_META_LABELS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_LABELS_OFF_HULL                       false
#define         DEFAULT_CHULL_LATEX_LABELS_OFF_HULL                       XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_LABELS_OFF_HULL")
#define AFLOWRC_DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION              -1
#define         DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION              XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_PLOT_REDUCED_COMPOSITION")
#define AFLOWRC_DEFAULT_CHULL_LATEX_HELVETICA_FONT                        false
#define         DEFAULT_CHULL_LATEX_HELVETICA_FONT                        XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_HELVETICA_FONT")
#define AFLOWRC_DEFAULT_CHULL_LATEX_FONT_SIZE                             std::string("")
#define         DEFAULT_CHULL_LATEX_FONT_SIZE                             XHOST.adefault.getattachedscheme("DEFAULT_CHULL_LATEX_FONT_SIZE")
#define AFLOWRC_DEFAULT_CHULL_LATEX_ROTATE_LABELS                         true
#define         DEFAULT_CHULL_LATEX_ROTATE_LABELS                         XHOST.adefault.getattachedutype<bool>("DEFAULT_CHULL_LATEX_ROTATE_LABELS")
#define AFLOWRC_DEFAULT_CHULL_LATEX_BOLD_LABELS                           -1
#define         DEFAULT_CHULL_LATEX_BOLD_LABELS                           XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_LATEX_BOLD_LABELS")
#define AFLOWRC_DEFAULT_CHULL_PNG_RESOLUTION                              300
#define         DEFAULT_CHULL_PNG_RESOLUTION                              XHOST.adefault.getattachedutype<int>("DEFAULT_CHULL_PNG_RESOLUTION")

// DEFAULT GFA
#define AFLOWRC_DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF                     0.05  //CO20190628
#define         DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF                     XHOST.adefault.getattachedutype<double>("DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF")  //CO20190628

// DEFAULT ARUN
#define AFLOWRC_ARUN_DIRECTORY_PREFIX                         std::string("ARUN.")
#define         ARUN_DIRECTORY_PREFIX                         XHOST.adefault.getattachedscheme("ARUN_DIRECTORY_PREFIX")

//DEFAULT POCC //CO20181226
#define AFLOWRC_DEFAULT_POCC_STRUCTURE_GENERATION_ALGO            std::string("UFF")
#define         DEFAULT_POCC_STRUCTURE_GENERATION_ALGO            XHOST.adefault.getattachedscheme("DEFAULT_POCC_STRUCTURE_GENERATION_ALGO")
#define AFLOWRC_DEFAULT_POCC_TEMPERATURE_STRING                   std::string("0:2400:300")
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
#define AFLOWRC_DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON  false
#define         DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON  XHOST.adefault.getattachedutype<bool>("DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON")
#define AFLOWRC_DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS             false
#define         DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS             XHOST.adefault.getattachedutype<bool>("DEFAULT_POCC_WRITE_OUT_ALL_SUPERCELLS")
#define AFLOWRC_POCC_FILE_PREFIX                                  std::string("aflow.pocc.")
#define         POCC_FILE_PREFIX                                  XHOST.adefault.getattachedscheme("POCC_FILE_PREFIX")
#define AFLOWRC_POCC_OUT_FILE                                     std::string("out")
#define         POCC_OUT_FILE                                     XHOST.adefault.getattachedscheme("POCC_OUT_FILE")
#define AFLOWRC_POCC_APL_OUT_FILE                                 std::string("apl.out")  //ME20210927
#define         POCC_APL_OUT_FILE                                 XHOST.adefault.getattachedscheme("POCC_APL_OUT_FILE")  //ME20210927
#define AFLOWRC_POCC_ALL_SUPERCELLS_FILE                          std::string("structures_all.out")
#define         POCC_ALL_SUPERCELLS_FILE                          XHOST.adefault.getattachedscheme("POCC_ALL_SUPERCELLS_FILE")
#define AFLOWRC_POCC_UNIQUE_SUPERCELLS_FILE                       std::string("structures_unique.out")
#define         POCC_UNIQUE_SUPERCELLS_FILE                       XHOST.adefault.getattachedscheme("POCC_UNIQUE_SUPERCELLS_FILE")
#define AFLOWRC_POCC_ALL_HNF_MATRICES_FILE                        std::string("hnf_matrices.out")
#define         POCC_ALL_HNF_MATRICES_FILE                        XHOST.adefault.getattachedscheme("POCC_ALL_HNF_MATRICES_FILE")
#define AFLOWRC_POCC_ALL_SITE_CONFIGURATIONS_FILE                 std::string("site_configurations.out")
#define         POCC_ALL_SITE_CONFIGURATIONS_FILE                 XHOST.adefault.getattachedscheme("POCC_ALL_SITE_CONFIGURATIONS_FILE")
#define AFLOWRC_POCC_DOSCAR_FILE                                  std::string("DOSCAR.pocc")
#define         POCC_DOSCAR_FILE                                  XHOST.adefault.getattachedscheme("POCC_DOSCAR_FILE")
#define AFLOWRC_POCC_PHDOSCAR_FILE                                std::string("PHDOSCAR.pocc")  //ME20210927
#define         POCC_PHDOSCAR_FILE                                XHOST.adefault.getattachedscheme("POCC_PHDOSCAR_FILE")  //ME20210927
#define AFLOWRC_POCC_ANIONS_LIST                                  std::string("B,C,N,O")
#define         POCC_ANIONS_LIST                                  XHOST.adefault.getattachedscheme("POCC_ANIONS_LIST")

// DEFAULT APL
//// DEFAULT APL SUPERCELL
#define AFLOWRC_DEFAULT_APL_PREC                              std::string("PHONONS")
#define         DEFAULT_APL_PREC                              XHOST.adefault.getattachedscheme("DEFAULT_APL_PREC")
#define AFLOWRC_DEFAULT_APL_ENGINE                            std::string("DM")
#define         DEFAULT_APL_ENGINE                            XHOST.adefault.getattachedscheme("DEFAULT_APL_ENGINE")
#define AFLOWRC_DEFAULT_APL_HIBERNATE                         true
#define         DEFAULT_APL_HIBERNATE                         XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_HIBERNATE")
#define AFLOWRC_DEFAULT_APL_MINSHELL                          6
#define         DEFAULT_APL_MINSHELL                          XHOST.adefault.getattachedutype<int>("DEFAULT_APL_MINSHELL")
#define AFLOWRC_DEFAULT_APL_MINATOMS                          175  //ME20190301
#define         DEFAULT_APL_MINATOMS                          XHOST.adefault.getattachedutype<int>("DEFAULT_APL_MINATOMS")
#define AFLOWRC_DEFAULT_APL_POLAR                             true  //CO20181226
#define         DEFAULT_APL_POLAR                             XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_POLAR")
#define AFLOWRC_DEFAULT_APL_DMAG                              0.015
#define         DEFAULT_APL_DMAG                              XHOST.adefault.getattachedutype<double>("DEFAULT_APL_DMAG")
#define AFLOWRC_DEFAULT_APL_DXYZONLY                          false
#define         DEFAULT_APL_DXYZONLY                          XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DXYZONLY")
#define AFLOWRC_DEFAULT_APL_DSYMMETRIZE                       true
#define         DEFAULT_APL_DSYMMETRIZE                       XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DSYMMETRIZE")
#define AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY                     true
#define         DEFAULT_APL_DINEQUIV_ONLY                     XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DINEQUIV_ONLY")
#define AFLOWRC_DEFAULT_APL_DPM                               std::string("ON")
#define         DEFAULT_APL_DPM                               XHOST.adefault.getattachedscheme("DEFAULT_APL_DPM")
#define AFLOWRC_DEFAULT_APL_RELAX                             true
#define         DEFAULT_APL_RELAX                             XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_RELAX")
#define AFLOWRC_DEFAULT_APL_RELAX_COMMENSURATE                true  //ME20200427
#define         DEFAULT_APL_RELAX_COMMENSURATE                XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_RELAX_COMMENSURATE")  //ME20200427
#define AFLOWRC_DEFAULT_APL_ZEROSTATE                         false  //CO2018121  //ME20220415 - ZEROSTATE=ON and DPM=ON is unnecessary
#define         DEFAULT_APL_ZEROSTATE                         XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_ZEROSTATE")
#define AFLOWRC_DEFAULT_APL_ZEROSTATE_CHGCAR                  false  //ME20191029
#define         DEFAULT_APL_ZEROSTATE_CHGCAR                  XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_ZEROSTATE_CHGCAR")  //ME20191029
#define AFLOWRC_DEFAULT_APL_USE_LEPSILON                      true
#define         DEFAULT_APL_USE_LEPSILON                      XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_USE_LEPSILON")

//// DEFAULT APL PHONON PROPERTIES
#define AFLOWRC_DEFAULT_APL_FREQFORMAT                        std::string("THZ|ALLOW_NEGATIVE")  //CO20181226 - no spaces!
#define         DEFAULT_APL_FREQFORMAT                        XHOST.adefault.getattachedscheme("DEFAULT_APL_FREQFORMAT")
#define AFLOWRC_DEFAULT_APL_DC                                true
#define         DEFAULT_APL_DC                                XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DC")
#define AFLOWRC_DEFAULT_APL_DCPATH                            std::string("LATTICE")
#define         DEFAULT_APL_DCPATH                            XHOST.adefault.getattachedscheme("DEFAULT_APL_DCPATH")
#define AFLOWRC_DEFAULT_APL_DCPOINTS                          100
#define         DEFAULT_APL_DCPOINTS                          XHOST.adefault.getattachedutype<int>("DEFAULT_APL_DCPOINTS")
#define AFLOWRC_DEFAULT_APL_DOS                               true
#define         DEFAULT_APL_DOS                               XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DOS")
#define AFLOWRC_DEFAULT_APL_DOSMETHOD                         std::string("LT")
#define         DEFAULT_APL_DOSMETHOD                         XHOST.adefault.getattachedscheme("DEFAULT_APL_DOSMETHOD")
#define AFLOWRC_DEFAULT_APL_DOSMESH                           std::string("21x21x21")
#define         DEFAULT_APL_DOSMESH                           XHOST.adefault.getattachedscheme("DEFAULT_APL_DOSMESH")
#define AFLOWRC_DEFAULT_APL_DOSPOINTS                         2000
#define         DEFAULT_APL_DOSPOINTS                         XHOST.adefault.getattachedutype<int>("DEFAULT_APL_DOSPOINTS")
#define AFLOWRC_DEFAULT_APL_DOSSMEAR                          0.0
#define         DEFAULT_APL_DOSSMEAR                          XHOST.adefault.getattachedutype<double>("DEFAULT_APL_DOSSMEAR")
#define AFLOWRC_DEFAULT_APL_DOS_PROJECT                       false  //ME20200213
#define         DEFAULT_APL_DOS_PROJECT                       XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DOS_PROJECT")  //ME20200213
#define AFLOWRC_DEFAULT_APL_TP                                true
#define         DEFAULT_APL_TP                                XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_TP")
#define AFLOWRC_DEFAULT_APL_DISPLACEMENTS                     true  //ME20200421
#define         DEFAULT_APL_DISPLACEMENTS                     XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DISPLACEMENTS")  //ME20200421
#define AFLOWRC_DEFAULT_APL_TPT                               std::string("0:2000:10")
#define         DEFAULT_APL_TPT                               XHOST.adefault.getattachedscheme("DEFAULT_APL_TPT")
#define AFLOWRC_DEFAULT_APL_GVEL                              true  // ME20200517
#define         DEFAULT_APL_GVEL                              XHOST.adefault.getattachedutype<bool>("DEFAULT_APL_DISPLACEMENTS")  // ME20200517

//// DEFAULT APL FILES
#define AFLOWRC_DEFAULT_APL_FILE_PREFIX                       std::string("aflow.apl.")
#define         DEFAULT_APL_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_APL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_APL_OUT_FILE                          std::string("out")  // ME20210927
#define         DEFAULT_APL_OUT_FILE                          XHOST.adefault.getattachedscheme("DEFAULT_APL_OUT_FILE")  // ME20210927
#define AFLOWRC_DEFAULT_APL_PDIS_FILE                         std::string("phonon_dispersion.out")
#define         DEFAULT_APL_PDIS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_APL_PDIS_FILE")
#define AFLOWRC_DEFAULT_APL_PDOS_FILE                         std::string("phonon_dos.out")
#define         DEFAULT_APL_PDOS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_APL_PDOS_FILE")
#define AFLOWRC_DEFAULT_APL_THERMO_FILE                       std::string("thermodynamic_properties.out")
#define         DEFAULT_APL_THERMO_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_APL_THERMO_FILE")
#define AFLOWRC_DEFAULT_APL_THERMO_JSON                       std::string("thermodynamic_properties.json")  //ME20211019
#define         DEFAULT_APL_THERMO_JSON                       XHOST.adefault.getattachedscheme("DEFAULT_APL_THERMO_JSON")  //ME20211019
#define AFLOWRC_DEFAULT_APL_DYNMAT_FILE                       std::string("DYNMAT.out")
#define         DEFAULT_APL_DYNMAT_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_APL_DYNMAT_FILE")
#define AFLOWRC_DEFAULT_APL_HARMIFC_FILE                      std::string("harmonicIFCs.xml")
#define         DEFAULT_APL_HARMIFC_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_APL_HARMIFC_FILE")
//ME20200415
#define AFLOWRC_DEFAULT_APL_POLAR_FILE                        std::string("polar.xml")
#define         DEFAULT_APL_POLAR_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_APL_POLAR_FILE")
#define AFLOWRC_DEFAULT_APL_HSKPTS_FILE                       std::string("hskpoints.out")
#define         DEFAULT_APL_HSKPTS_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_APL_HSKPTS_FILE")
#define AFLOWRC_DEFAULT_APL_MSQRDISP_FILE                     std::string("displacements.out")  //ME20200329
#define         DEFAULT_APL_MSQRDISP_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_APL_MSQRDISP_FILE")  // ME20200329
#define AFLOWRC_DEFAULT_APL_GVEL_FILE                         std::string("group_velocities.out")  //ME202000517
#define         DEFAULT_APL_GVEL_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_APL_GVEL_FILE")  //ME20200517
//ME20190614 BEGIN
#define AFLOWRC_DEFAULT_APL_PHDOSCAR_FILE                     std::string("PHDOSCAR")
#define         DEFAULT_APL_PHDOSCAR_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_APL_PHDOSCAR_FILE")
#define AFLOWRC_DEFAULT_APL_PHPOSCAR_FILE                     std::string("PHPOSCAR")
#define         DEFAULT_APL_PHPOSCAR_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_APL_PHPOSCAR_FILE")
#define AFLOWRC_DEFAULT_APL_PHKPOINTS_FILE                    std::string("PHKPOINTS")
#define         DEFAULT_APL_PHKPOINTS_FILE                    XHOST.adefault.getattachedscheme("DEFAULT_APL_PHKPOINTS_FILE")
#define AFLOWRC_DEFAULT_APL_PHEIGENVAL_FILE                   std::string("PHEIGENVAL")
#define         DEFAULT_APL_PHEIGENVAL_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_APL_PHEIGENVAL_FILE")
//ME20190614 END
#define AFLOWRC_DEFAULT_APL_STATE_FILE                        std::string("fccalc_state.out")
#define         DEFAULT_APL_STATE_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_APL_STATE_FILE")  //ME20200224

//ME20200329 BEGIN
#define AFLOWRC_DEFAULT_APL_ADISP_SCENE_FORMAT                std::string("XCRYSDEN")
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
#define AFLOWRC_DEFAULT_QHA_MODE                              std::string("QHA")
#define         DEFAULT_QHA_MODE                              XHOST.adefault.getattachedscheme("DEFAULT_QHA_MODE")
#define AFLOWRC_DEFAULT_QHA_EOS                               true
#define         DEFAULT_QHA_EOS                               XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_EOS")
#define AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE              std::string("-12:16:3")
#define         DEFAULT_QHA_EOS_DISTORTION_RANGE              XHOST.adefault.getattachedscheme("DEFAULT_QHA_EOS_DISTORTION_RANGE")
//AS20200818 BEGIN
#define AFLOWRC_DEFAULT_QHA_EOS_MODEL                         std::string("SJ")
#define         DEFAULT_QHA_EOS_MODEL                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_EOS_MODEL")
//AS20200818 END
#define AFLOWRC_DEFAULT_QHA_GP_DISTORTION                     1.0
#define         DEFAULT_QHA_GP_DISTORTION                     XHOST.adefault.getattachedutype<double>("DEFAULT_QHA_GP_DISTORTION")
//AS20200602 BEGIN
#define AFLOWRC_DEFAULT_QHA_TAYLOR_EXPANSION_ORDER            2
#define         DEFAULT_QHA_TAYLOR_EXPANSION_ORDER            XHOST.adefault.getattachedutype<double>("DEFAULT_QHA_TAYLOR_EXPANSION_ORDER")
//AS20200602 END
#define AFLOWRC_DEFAULT_QHA_INCLUDE_ELEC_CONTRIB              false
#define         DEFAULT_QHA_INCLUDE_ELEC_CONTRIB              XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_INCLUDE_ELEC_CONTRIB")
//AS20200528 BEGIN
#define AFLOWRC_DEFAULT_QHA_SOMMERFELD_EXPANSION              false
#define         DEFAULT_QHA_SOMMERFELD_EXPANSION              XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_SOMMERFELD_EXPANSION")
//AS20200528 END
#define AFLOWRC_DEFAULT_QHA_PDIS_T                            std::string("300")
#define         DEFAULT_QHA_PDIS_T                            XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDIS_T")
//AS20200508 BEGIN
#define AFLOWRC_DEFAULT_QHA_GP_FINITE_DIFF                    false
#define         DEFAULT_QHA_GP_FINITE_DIFF                    XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_GP_FINITE_DIFF")
#define AFLOWRC_DEFAULT_QHA_IGNORE_IMAGINARY                  false
#define         DEFAULT_QHA_IGNORE_IMAGINARY                  XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_IGNORE_IMAGINARY")
//AS20201123 BEGIN
#define AFLOWRC_DEFAULT_QHA_RELAX_IONS_CELL                   false
#define         DEFAULT_QHA_RELAX_IONS_CELL                   XHOST.adefault.getattachedutype<bool>("DEFAULT_QHA_RELAX_IONS_CELL")
//AS20201123 END

//// DEFAULT QHA FILES
#define AFLOWRC_DEFAULT_QHA_FILE_PREFIX                       std::string("aflow.qha.")
#define         DEFAULT_QHA_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_QHA_FILE_PREFIX")
//AS20200709 BEGIN
#define AFLOWRC_DEFAULT_QHA3P_FILE_PREFIX                     std::string("aflow.qha3p.")
#define         DEFAULT_QHA3P_FILE_PREFIX                     XHOST.adefault.getattachedscheme("DEFAULT_QHA3P_FILE_PREFIX")
#define AFLOWRC_DEFAULT_QHANP_FILE_PREFIX                     std::string("aflow.qhanp.")
#define         DEFAULT_QHANP_FILE_PREFIX                     XHOST.adefault.getattachedscheme("DEFAULT_QHANP_FILE_PREFIX")
#define AFLOWRC_DEFAULT_SCQHA_FILE_PREFIX                     std::string("aflow.scqha.")
#define         DEFAULT_SCQHA_FILE_PREFIX                     XHOST.adefault.getattachedscheme("DEFAULT_SCQHA_FILE_PREFIX")
//AS20200709 END
#define AFLOWRC_DEFAULT_QHA_GP_PATH_FILE                      std::string("gp.disp.out")
#define         DEFAULT_QHA_GP_PATH_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_PATH_FILE")
#define AFLOWRC_DEFAULT_QHA_GP_MESH_FILE                      std::string("gp.mesh.out")
#define         DEFAULT_QHA_GP_MESH_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_MESH_FILE")
#define AFLOWRC_DEFAULT_QHA_GP_AVG_FILE                       std::string("gp.avg.out")
#define         DEFAULT_QHA_GP_AVG_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_QHA_GP_AVG_FILE")
#define AFLOWRC_DEFAULT_QHA_THERMO_FILE                       std::string("thermo.out")
#define         DEFAULT_QHA_THERMO_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_QHA_THERMO_FILE")
#define AFLOWRC_DEFAULT_QHA_FREQS_FILE                        std::string("frequencies.out")
#define         DEFAULT_QHA_FREQS_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_QHA_FREQS_FILE")
#define AFLOWRC_DEFAULT_QHA_FVT_FILE                          std::string("FVT.out")
#define         DEFAULT_QHA_FVT_FILE                          XHOST.adefault.getattachedscheme("DEFAULT_QHA_FVT_FILE")
//AS20200508 END
//AS20210517 BEGIN
#define AFLOWRC_DEFAULT_QHA_COEFF_FILE                        std::string("coeff.out")
#define         DEFAULT_QHA_COEFF_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_QHA_COEFF_FILE")
#define AFLOWRC_DEFAULT_QHA_IMAG_FILE                         std::string("imag.out")
#define         DEFAULT_QHA_IMAG_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_IMAG_FILE")
//AS20210517 END
//AS20201022 BEGIN
#define AFLOWRC_DEFAULT_QHA_PDIS_FILE                         std::string("dispersion_phonon")
#define         DEFAULT_QHA_PDIS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDIS_FILE")
//AS20201022 END
//AS20201201 BEGIN
#define AFLOWRC_DEFAULT_QHA_PDOS_FILE                         std::string("dos_phonon")
#define         DEFAULT_QHA_PDOS_FILE                         XHOST.adefault.getattachedscheme("DEFAULT_QHA_PDOS_FILE")
//AS20201201 END
//AS20201112 BEGIN
#define AFLOWRC_DEFAULT_QHA_KPOINTS_FILE                      std::string("kpoints.out")
#define         DEFAULT_QHA_KPOINTS_FILE                      XHOST.adefault.getattachedscheme("DEFAULT_QHA_KPOINTS_FILE")
//AS20201112 
//AS20210914 BEGIN
#define AFLOWRC_DEFAULT_POCC_QHA_THERMO_FILE                  std::string("qha.thermo.out")
#define         DEFAULT_POCC_QHA_THERMO_FILE                  XHOST.adefault.getattachedscheme("DEFAULT_POCC_QHA_THERMO_FILE")
#define AFLOWRC_DEFAULT_POCC_QHA_AVGTHERMO_FILE               std::string("qha.avgthermo.out")
#define         DEFAULT_POCC_QHA_AVGTHERMO_FILE               XHOST.adefault.getattachedscheme("DEFAULT_POCC_QHA_AVGTHERMO_FILE")
//AS20210914 END

// DEFAULT AAPL
//// DEFAULT AAPL VALUES
#define AFLOWRC_DEFAULT_AAPL_BTE                              std::string("FULL")
#define         DEFAULT_AAPL_BTE                              XHOST.adefault.getattachedscheme("DEFAULT_AAPL_BTE")
#define AFLOWRC_DEFAULT_AAPL_FOURTH_ORDER                     false
#define         DEFAULT_AAPL_FOURTH_ORDER                     XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_FOURTH_ORDER")
#define AFLOWRC_DEFAULT_AAPL_CUT_RAD                          std::string("0.0") //ME20190308 - use CUT_SHELL by default //ME20191029
#define         DEFAULT_AAPL_CUT_RAD                          XHOST.adefault.getattachedscheme("DEFAULT_AAPL_CUT_RAD")
#define AFLOWRC_DEFAULT_AAPL_CUT_SHELL                        std::string("6")  //ME20190301  //ME20190408  //ME20191029
#define         DEFAULT_AAPL_CUT_SHELL                        XHOST.adefault.getattachedscheme("DEFAULT_AAPL_CUT_SHELL")
#define AFLOWRC_DEFAULT_AAPL_THERMALGRID                      std::string("21x21x21")  //ME20200110 - 21x21x21 more than enough for tetrahedron method; odd is preferred (Gamma-centered by construction)
#define         DEFAULT_AAPL_THERMALGRID                      XHOST.adefault.getattachedscheme("DEFAULT_AAPL_THERMALGRID")
#define AFLOWRC_DEFAULT_AAPL_TCT                              std::string("50:550:50")
#define         DEFAULT_AAPL_TCT                              XHOST.adefault.getattachedscheme("DEFAULT_AAPL_TCT")
#define AFLOWRC_DEFAULT_AAPL_SUMRULE                          1e-7
#define         DEFAULT_AAPL_SUMRULE                          XHOST.adefault.getattachedutype<double>("DEFAULT_AAPL_SUMRULE")
#define AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER                 2000
#define         DEFAULT_AAPL_SUMRULE_MAX_ITER                 XHOST.adefault.getattachedutype<int>("DEFAULT_AAPL_SUMRULE_MAX_ITER")
#define AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT               0.0
#define         DEFAULT_AAPL_MIXING_COEFFICIENT               XHOST.adefault.getattachedutype<double>("DEFAULT_AAPL_MIXING_COEFFICIENT")
#define AFLOWRC_DEFAULT_AAPL_ISOTOPE                          true
#define         DEFAULT_AAPL_ISOTOPE                          XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_ISOTOPE")
#define AFLOWRC_DEFAULT_AAPL_BOUNDARY                         false
#define         DEFAULT_AAPL_BOUNDARY                         XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_BOUNDARY")
#define AFLOWRC_DEFAULT_AAPL_CUMULATIVEK                      false
#define         DEFAULT_AAPL_CUMULATIVEK                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AAPL_CUMULATIVEK")
#define AFLOWRC_DEFAULT_AAPL_NANO_SIZE                        100.0
#define         DEFAULT_AAPL_NANO_SIZE                        XHOST.adefault.getattachedutype<double>("DEFAULT_AAPL_NANO_SIZE")

//// DEFAULT AAPL FILES
#define AFLOWRC_DEFAULT_AAPL_FILE_PREFIX                      std::string("aflow.aapl.")
#define         DEFAULT_AAPL_FILE_PREFIX                      XHOST.adefault.getattachedscheme("DEFAULT_AAPL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_AAPL_IRRQPTS_FILE                     std::string("irred_qpoints.out")
#define         DEFAULT_AAPL_IRRQPTS_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_AAPL_IRRQPTS_FILE")
#define AFLOWRC_DEFAULT_AAPL_GVEL_FILE                        std::string("group_velocities.out")
#define         DEFAULT_AAPL_GVEL_FILE                        XHOST.adefault.getattachedscheme("DEFAULT_AAPL_GVEL_FILE")
#define AFLOWRC_DEFAULT_AAPL_PS_FILE                          std::string("phase_space.out")  //ME20191104
#define         DEFAULT_AAPL_PS_FILE                          XHOST.adefault.getattachedscheme("DEFAULT_AAPL_PS_FILE")  //ME20191104
#define AFLOWRC_DEFAULT_AAPL_GRUENEISEN_FILE                  std::string("grueneisen.out")  //ME20191104
#define         DEFAULT_AAPL_GRUENEISEN_FILE                  XHOST.adefault.getattachedscheme("DEFAULT_AAPL_GRUENEISEN_FILE")  //ME20191104
#define AFLOWRC_DEFAULT_AAPL_RATES_FILE                       std::string("scattering_rates_total.out")
#define         DEFAULT_AAPL_RATES_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_FILE")
#define AFLOWRC_DEFAULT_AAPL_RATES_3RD_FILE                   std::string("scattering_rates_anharmonic_3rd.out")
#define         DEFAULT_AAPL_RATES_3RD_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_3RD_FILE")
#define AFLOWRC_DEFAULT_AAPL_RATES_4TH_FILE                   DEFAULT_AAPL_FILE_PREFIX+std::string("scattering_rates_anharmonic_4th.out")
#define         DEFAULT_AAPL_RATES_4TH_FILE                   XHOST.adefault.getattachedscheme("DEFAULT_AAPL_RATES_4TH_FILE")
#define AFLOWRC_DEFAULT_AAPL_ISOTOPE_FILE                     std::string("scattering_rates_isotope.out")
#define         DEFAULT_AAPL_ISOTOPE_FILE                     XHOST.adefault.getattachedscheme("DEFAULT_AAPL_ISOTOPE_FILE")
#define AFLOWRC_DEFAULT_AAPL_BOUNDARY_FILE                    std::string("scattering_rates_boundary.out")
#define         DEFAULT_AAPL_BOUNDARY_FILE                    XHOST.adefault.getattachedscheme("DEFAULT_AAPL_BOUNDARY_FILE")
#define AFLOWRC_DEFAULT_AAPL_TCOND_FILE                       std::string("thermal_conductivity.out")
#define         DEFAULT_AAPL_TCOND_FILE                       XHOST.adefault.getattachedscheme("DEFAULT_AAPL_TCOND_FILE")


// DEFAULT AEL
//// DEFAULT AEL STRAIN CALCS
#define AFLOWRC_DEFAULT_AEL_STRAIN_SYMMETRY                   true
#define         DEFAULT_AEL_STRAIN_SYMMETRY                   XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_STRAIN_SYMMETRY")
#define AFLOWRC_DEFAULT_AEL_NNORMAL_STRAINS                   4
#define         DEFAULT_AEL_NNORMAL_STRAINS                   XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_NNORMAL_STRAINS")
#define AFLOWRC_DEFAULT_AEL_NSHEAR_STRAINS                    4
#define         DEFAULT_AEL_NSHEAR_STRAINS                    XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_NSHEAR_STRAINS")
#define AFLOWRC_DEFAULT_AEL_NORMAL_STRAIN_STEP                0.005
#define         DEFAULT_AEL_NORMAL_STRAIN_STEP                XHOST.adefault.getattachedutype<double>("DEFAULT_AEL_NORMAL_STRAIN_STEP")
#define AFLOWRC_DEFAULT_AEL_SHEAR_STRAIN_STEP                 0.005
#define         DEFAULT_AEL_SHEAR_STRAIN_STEP                 XHOST.adefault.getattachedutype<double>("DEFAULT_AEL_SHEAR_STRAIN_STEP")
#define AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_CALC                false
#define         DEFAULT_AEL_ORIGIN_STRAIN_CALC                XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_ORIGIN_STRAIN_CALC")
#define AFLOWRC_DEFAULT_AEL_ORIGIN_STRAIN_FIT                 false
#define         DEFAULT_AEL_ORIGIN_STRAIN_FIT                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_ORIGIN_STRAIN_FIT")
#define AFLOWRC_DEFAULT_AEL_RELAXED_STRUCT_FIT                false
#define         DEFAULT_AEL_RELAXED_STRUCT_FIT                XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_RELAXED_STRUCT_FIT")
#define AFLOWRC_DEFAULT_AEL_NEG_STRAINS                       true
#define         DEFAULT_AEL_NEG_STRAINS                       XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_NEG_STRAINS")
#define AFLOWRC_DEFAULT_AEL_NIND_STRAIN_DIRS                  3
#define         DEFAULT_AEL_NIND_STRAIN_DIRS                  XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_NIND_STRAIN_DIRS")
#define AFLOWRC_DEFAULT_AEL_VASPSYM                           false
#define         DEFAULT_AEL_VASPSYM                           XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_VASPSYM")
#define AFLOWRC_DEFAULT_AEL_PRECACC_ALGONORM                  false
#define         DEFAULT_AEL_PRECACC_ALGONORM                  XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_PRECACC_ALGONORM")
#define AFLOWRC_DEFAULT_AEL_VASPRUNXML_STRESS                 false
#define         DEFAULT_AEL_VASPRUNXML_STRESS                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_VASPRUNXML_STRESS")
#define AFLOWRC_DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS             false
#define         DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS             XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_AUTOSKIP_FAILED_ARUNS")
#define AFLOWRC_DEFAULT_AEL_SKIP_ARUNS_MAX                    1
#define         DEFAULT_AEL_SKIP_ARUNS_MAX                    XHOST.adefault.getattachedutype<int>("DEFAULT_AEL_SKIP_ARUNS_MAX")

//// DEFAULT AEL CHECKS AND PROCESSING
#define AFLOWRC_DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY            true
#define         DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY            XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_CHECK_ELASTIC_SYMMETRY")
#define AFLOWRC_DEFAULT_AEL_SYMMETRIZE                        false
#define         DEFAULT_AEL_SYMMETRIZE                        XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_SYMMETRIZE")

//// DEFAULT AEL OUTPUT FILES
#define AFLOWRC_DEFAULT_AEL_FILE_PREFIX                       std::string("aflow.ael.")
#define         DEFAULT_AEL_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_AEL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_AEL_WRITE_FULL_RESULTS                false
#define         DEFAULT_AEL_WRITE_FULL_RESULTS                XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_WRITE_FULL_RESULTS")
#define AFLOWRC_DEFAULT_AEL_DIRNAME_ARUN                      true
#define         DEFAULT_AEL_DIRNAME_ARUN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AEL_DIRNAME_ARUN")

// DEFAULT AGL
//// DEFAULT AGL STRAIN CALCS
#define AFLOWRC_DEFAULT_AGL_AEL_POISSON_RATIO                 true
#define         DEFAULT_AGL_AEL_POISSON_RATIO                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_AEL_POISSON_RATIO")
#define AFLOWRC_DEFAULT_AGL_NSTRUCTURES                       28
#define         DEFAULT_AGL_NSTRUCTURES                       XHOST.adefault.getattachedutype<int>("DEFAULT_AGL_NSTRUCTURES")
#define AFLOWRC_DEFAULT_AGL_STRAIN_STEP                       0.01
#define         DEFAULT_AGL_STRAIN_STEP                       XHOST.adefault.getattachedutype<double>("DEFAULT_AGL_STRAIN_STEP")
#define AFLOWRC_DEFAULT_AGL_AUTOSKIP_FAILED_ARUNS             false
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
#define AFLOWRC_DEFAULT_AGL_CHECK_EV_CONCAVITY                false
#define         DEFAULT_AGL_CHECK_EV_CONCAVITY                XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_CHECK_EV_CONCAVITY")
#define AFLOWRC_DEFAULT_AGL_CHECK_EV_MIN                      false
#define         DEFAULT_AGL_CHECK_EV_MIN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_CHECK_EV_MIN")
#define AFLOWRC_DEFAULT_AGL_HUGONIOT_CALC                     true
#define         DEFAULT_AGL_HUGONIOT_CALC                     XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_HUGONIOT_CALC")
#define AFLOWRC_DEFAULT_AGL_HUGONIOT_EXTRAPOLATE              false
#define         DEFAULT_AGL_HUGONIOT_EXTRAPOLATE              XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_HUGONIOT_EXTRAPOLATE")
#define AFLOWRC_DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE      false
#define         DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_RUN_ALL_PRESSURE_TEMPERATURE")

//// DEFAULT AGL OUTPUT FILES
#define AFLOWRC_DEFAULT_AGL_FILE_PREFIX                       std::string("aflow.agl.")
#define         DEFAULT_AGL_FILE_PREFIX                       XHOST.adefault.getattachedscheme("DEFAULT_AGL_FILE_PREFIX")
#define AFLOWRC_DEFAULT_AGL_WRITE_FULL_RESULTS                false
#define         DEFAULT_AGL_WRITE_FULL_RESULTS                XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_WRITE_FULL_RESULTS")
#define AFLOWRC_DEFAULT_AGL_DIRNAME_ARUN                      true
#define         DEFAULT_AGL_DIRNAME_ARUN                      XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_DIRNAME_ARUN")
#define AFLOWRC_DEFAULT_AGL_WRITE_GIBBS_INPUT                 false
#define         DEFAULT_AGL_WRITE_GIBBS_INPUT                 XHOST.adefault.getattachedutype<bool>("DEFAULT_AGL_WRITE_GIBBS_INPUT")
#define AFLOWRC_DEFAULT_AGL_PLOT_RESULTS                      false
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
#define AFLOWRC_DEFAULT_QCA_PRINT                             std::string("txt")
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
#define AFLOWRC_DEFAULT_XTALFINDER_SUPERCELL_METHOD           false // supercell method for comparing (robust, but slow, superceded by transformation method)
#define         DEFAULT_XTALFINDER_SUPERCELL_METHOD           XHOST.adefault.getattachedutype<bool>("DEFAULT_XTALFINDER_SUPERCELL_METHOD") //DX20201223
//DX20200709 - START
#define AFLOWRC_DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING    4.0 // factor that divides minimum interatomic distance
#define         DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING    XHOST.adefault.getattachedutype<double>("DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING")
//DX20200709 - END
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL                    std::string("material_comparison_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_MATERIAL                    XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_MATERIAL") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE                   std::string("structure_comparison_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_STRUCTURE                   XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_STRUCTURE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_DUPLICATE                   std::string("duplicate_compounds_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_DUPLICATE                   XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_DUPLICATE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE   std::string("material_comparison_compare2database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE   XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE  std::string("structure_comparison_compare2database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE  XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE           std::string("material_comparison_database_output") // results file prefix
#define         DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE           XHOST.adefault.getattachedscheme("DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE") //DX20201228
#define AFLOWRC_DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE          std::string("structure_comparison_database_output") // results file prefix
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
#define AFLOWRC_FILE_VASP_MONITOR                             std::string("monitor_vasp")
#define         FILE_VASP_MONITOR                             XHOST.adefault.getattachedscheme("FILE_VASP_MONITOR")
#define AFLOWRC_INTEL_COMPILER_PATHS                          std::string("/opt/intel/bin/compilervars.sh,/opt/intel/bin/compilervars.csh,/app/intel/parallel_studio_xe_2020_update1/bin/compilervars.sh")
#define         INTEL_COMPILER_PATHS                          XHOST.adefault.getattachedscheme("INTEL_COMPILER_PATHS")

// MACHINE DEPENDENT MPI
#define AFLOWRC_MPI_OPTIONS_DUKE_BETA_MPICH                   std::string("ulimit -s unlimited ") // DUKE_BETA_MPICH
#define         MPI_OPTIONS_DUKE_BETA_MPICH                   XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_MPICH")
#define AFLOWRC_MPI_COMMAND_DUKE_BETA_MPICH                   std::string("/usr/bin/mpiexec -np") // DUKE_BETA_MPICH
#define         MPI_COMMAND_DUKE_BETA_MPICH                   XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_MPICH")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_MPICH                std::string("/usr/local/bin/") // DUKE_BETA_MPICH
#define         MPI_BINARY_DIR_DUKE_BETA_MPICH                XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_MPICH")

#define AFLOWRC_MPI_OPTIONS_DUKE_BETA_OPENMPI                 std::string("ulimit -s unlimited ") // DUKE_BETA_OPENMPI
#define         MPI_OPTIONS_DUKE_BETA_OPENMPI                 XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_BETA_OPENMPI")
#define AFLOWRC_MPI_COMMAND_DUKE_BETA_OPENMPI                 std::string("/usr/bin/mpirun.openmpi -np") // DUKE_BETA_OPENMPI
#define         MPI_COMMAND_DUKE_BETA_OPENMPI                 XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_BETA_OPENMPI")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_BETA_OPENMPI              std::string("/usr/local/bin/") // DUKE_BETA_OPENMPI
#define         MPI_BINARY_DIR_DUKE_BETA_OPENMPI              XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_BETA_OPENMPI")

#define AFLOWRC_MPI_OPTIONS_DUKE_MATERIALS                    std::string("ulimit -s unlimited ") // DUKE_MATERIALS
#define         MPI_OPTIONS_DUKE_MATERIALS                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_MATERIALS")
#define AFLOWRC_MPI_COMMAND_DUKE_MATERIALS                    std::string("/usr/bin/mpiexec -np") // DUKE_MATERIALS
#define         MPI_COMMAND_DUKE_MATERIALS                    XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_MATERIALS")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_MATERIALS                 std::string("/usr/local/bin/")  // DUKE_MATERIALS
#define         MPI_BINARY_DIR_DUKE_MATERIALS                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_MATERIALS")

#define AFLOWRC_MPI_OPTIONS_DUKE_AFLOWLIB                     std::string("ulimit -s unlimited ") // DUKE_AFLOWLIB
#define         MPI_OPTIONS_DUKE_AFLOWLIB                     XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_AFLOWLIB")
#define AFLOWRC_MPI_COMMAND_DUKE_AFLOWLIB                     std::string("/usr/bin/mpiexec -np") // DUKE_AFLOWLIB
#define         MPI_COMMAND_DUKE_AFLOWLIB                     XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_AFLOWLIB")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_AFLOWLIB                  std::string("/usr/local/bin/") // DUKE_AFLOWLIB
#define         MPI_BINARY_DIR_DUKE_AFLOWLIB                  XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_AFLOWLIB")

#define AFLOWRC_MPI_OPTIONS_DUKE_HABANA                       std::string("ulimit -s unlimited ") // DUKE_HABANA
#define         MPI_OPTIONS_DUKE_HABANA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_HABANA")
#define AFLOWRC_MPI_COMMAND_DUKE_HABANA                       std::string("/usr/bin/mpiexec -np") // DUKE_HABANA
#define         MPI_COMMAND_DUKE_HABANA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_HABANA")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_HABANA                    std::string("/usr/local/bin/") // DUKE_HABANA
#define         MPI_BINARY_DIR_DUKE_HABANA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_HABANA")

#define AFLOWRC_MPI_OPTIONS_DUKE_QRATS_MPICH                  std::string("ulimit -s unlimited ") // DUKE_QRATS_MPICH
#define         MPI_OPTIONS_DUKE_QRATS_MPICH                  XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QRATS_MPICH")
#define AFLOWRC_MPI_COMMAND_DUKE_QRATS_MPICH                  std::string("/MAIN/bin/MPICH/bin/mpirun -np") // DUKE_QRATS_MPICH
#define         MPI_COMMAND_DUKE_QRATS_MPICH                  XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QRATS_MPICH")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_QRATS_MPICH               std::string("/usr/local/bin/") // DUKE_QRATS_MPICH
#define         MPI_BINARY_DIR_DUKE_QRATS_MPICH               XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QRATS_MPICH")

#define AFLOWRC_MPI_OPTIONS_DUKE_QFLOW_OPENMPI                std::string("ulimit -s unlimited ") // DUKE_QFLOW_MPICH
#define         MPI_OPTIONS_DUKE_QFLOW_OPENMPI                XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_QFLOW_OPENMPI")
#define AFLOWRC_MPI_COMMAND_DUKE_QFLOW_OPENMPI                std::string("/home/bin/local/bin/mpirun -n") // DUKE_QFLOW_MPICH
#define         MPI_COMMAND_DUKE_QFLOW_OPENMPI                XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_QFLOW_OPENMPI")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI             std::string("/home/bin/") // DUKE_QFLOW_MPICH
#define         MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI             XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI")

//CO20201220 X START
#define AFLOWRC_MPI_OPTIONS_DUKE_X_X                          std::string("ulimit -s unlimited ") // DUKE_X_X_MPICH
#define         MPI_OPTIONS_DUKE_X_X                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_X")
#define AFLOWRC_MPI_COMMAND_DUKE_X_X                          std::string("srun --mpi=pmix --cpus-per-task") // DUKE_X_X_MPICH
#define         MPI_COMMAND_DUKE_X_X                          XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_X")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_X                       std::string("~/bin/x/") // DUKE_X_X_MPICH
#define         MPI_BINARY_DIR_DUKE_X_X                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_X")
#define AFLOWRC_MPI_OPTIONS_DUKE_X_CRAY                       std::string("ulimit -s unlimited ") // DUKE_X_CRAY_MPICH
#define         MPI_OPTIONS_DUKE_X_CRAY                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_CRAY")
#define AFLOWRC_MPI_COMMAND_DUKE_X_CRAY                       std::string("srun --mpi=pmix --cpus-per-task") // DUKE_X_CRAY_MPICH
#define         MPI_COMMAND_DUKE_X_CRAY                       XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_CRAY")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_CRAY                    std::string("~/bin/cray/") // DUKE_X_CRAY_MPICH
#define         MPI_BINARY_DIR_DUKE_X_CRAY                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_CRAY")
#define AFLOWRC_MPI_OPTIONS_DUKE_X_OLDCRAY                    std::string("ulimit -s unlimited ") // DUKE_X_OLDCRAY_MPICH
#define         MPI_OPTIONS_DUKE_X_OLDCRAY                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_OLDCRAY")
#define AFLOWRC_MPI_COMMAND_DUKE_X_OLDCRAY                    std::string("srun --mpi=pmix --cpus-per-task") // DUKE_X_OLDCRAY_MPICH
#define         MPI_COMMAND_DUKE_X_OLDCRAY                    XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_OLDCRAY")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_OLDCRAY                 std::string("~/bin/oldcray/") // DUKE_X_OLDCRAY_MPICH
#define         MPI_BINARY_DIR_DUKE_X_OLDCRAY                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_OLDCRAY")
#define AFLOWRC_MPI_OPTIONS_DUKE_X_SMB                        std::string("ulimit -s unlimited ") // DUKE_X_SMB_MPICH
#define         MPI_OPTIONS_DUKE_X_SMB                        XHOST.adefault.getattachedscheme("MPI_OPTIONS_DUKE_X_SMB")
#define AFLOWRC_MPI_COMMAND_DUKE_X_SMB                        std::string("srun --mpi=pmix --cpus-per-task") // DUKE_X_SMB_MPICH
#define         MPI_COMMAND_DUKE_X_SMB                        XHOST.adefault.getattachedscheme("MPI_COMMAND_DUKE_X_SMB")
#define AFLOWRC_MPI_BINARY_DIR_DUKE_X_SMB                     std::string("~/bin/smb/") // DUKE_X_SMB_MPICH
#define         MPI_BINARY_DIR_DUKE_X_SMB                     XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_DUKE_X_SMB")
//CO20201220 X STOP

//CO20220818 JHU_ROCKFISH START
#define AFLOWRC_MPI_OPTIONS_JHU_ROCKFISH                      std::string("ulimit -s unlimited ") // JHU_ROCKFISH_MPICH
#define         MPI_OPTIONS_JHU_ROCKFISH                      XHOST.adefault.getattachedscheme("MPI_OPTIONS_JHU_ROCKFISH")
#define AFLOWRC_MPI_COMMAND_JHU_ROCKFISH                      std::string("mpirun -n") // JHU_ROCKFISH_MPICH
#define         MPI_COMMAND_JHU_ROCKFISH                      XHOST.adefault.getattachedscheme("MPI_COMMAND_JHU_ROCKFISH")
#define AFLOWRC_MPI_BINARY_DIR_JHU_ROCKFISH                   std::string("~/bin/") // JHU_ROCKFISH_MPICH
#define         MPI_BINARY_DIR_JHU_ROCKFISH                   XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_JHU_ROCKFISH")
//CO20220818 JHU_ROCKFISH STOP

//DX20190509 - MACHINE001 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE001                        std::string("") // MACHINE001
#define         MPI_OPTIONS_MACHINE001                        XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE001")
#define AFLOWRC_MPI_COMMAND_MACHINE001                        std::string("aprun -n") // MACHINE001
#define         MPI_COMMAND_MACHINE001                        XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE001")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE001                     std::string("~/bin/") // MACHINE001
#define         MPI_BINARY_DIR_MACHINE001                     XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE001")
//DX20190509 - MACHINE001 - END

//DX20190509 - MACHINE002 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE002                       std::string("") // MACHINE002
#define         MPI_OPTIONS_MACHINE002                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE002")
#define AFLOWRC_MPI_COMMAND_MACHINE002                       std::string("mpirun -np") // MACHINE002
#define         MPI_COMMAND_MACHINE002                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE002")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE002                    std::string("~/bin/") // MACHINE002
#define         MPI_BINARY_DIR_MACHINE002                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE002")
//DX20190509 - MACHINE002 - END

//DX20201005 - MACHINE003 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE003                       std::string("") // MACHINE003
#define         MPI_OPTIONS_MACHINE003                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE003")
#define AFLOWRC_MPI_COMMAND_MACHINE003                       std::string("mpiexec -n") // MACHINE003
#define         MPI_COMMAND_MACHINE003                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE003")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE003                    std::string("~/bin/") // MACHINE003
#define         MPI_BINARY_DIR_MACHINE003                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE003")
//DX20201005 - MACHINE003 - END

//DX20211011 - MACHINE004 - START
#define AFLOWRC_MPI_OPTIONS_MACHINE004                       std::string("") // MACHINE004
#define         MPI_OPTIONS_MACHINE004                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE004")
#define AFLOWRC_MPI_COMMAND_MACHINE004                       std::string("mpirun --mca pml ob1 --mca btl_ofi_mode 2 -np") // MACHINE004
#define         MPI_COMMAND_MACHINE004                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE004")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE004                    std::string("~/bin/") // MACHINE004
#define         MPI_BINARY_DIR_MACHINE004                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE004")
//DX20211011 - MACHINE004 - END

#define AFLOWRC_MPI_OPTIONS_MPCDF_EOS                         std::string("ulimit -s unlimited ") // MPCDF_EOS_MPICH
#define         MPI_OPTIONS_MPCDF_EOS                         XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_EOS")
#define AFLOWRC_MPI_COMMAND_MPCDF_EOS                         std::string("/usr/bin/srun -n") // MPCDF_EOS_MPICH
#define         MPI_COMMAND_MPCDF_EOS                         XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_EOS")
#define AFLOWRC_MPI_NCPUS_MPCDF_EOS                           32 // 32 // MPCDF_EOS_MPICH
#define         MPI_NCPUS_MPCDF_EOS                           XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_EOS")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_EOS                  std::string("NEGLECT")  // false/OFF, IGNORE/NEGLECT, true/ON
#define         MPI_HYPERTHREADING_MPCDF_EOS                  XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_EOS") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_EOS                      std::string("~/bin/") // MPCDF_EOS_MPICH
#define         MPI_BINARY_DIR_MPCDF_EOS                      XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_EOS")

#define AFLOWRC_MPI_OPTIONS_MPCDF_DRACO                       std::string("ulimit -s unlimited ") // MPCDF_DRACO_MPICH  // FIX_DRACO
#define         MPI_OPTIONS_MPCDF_DRACO                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_DRACO")  // FIX_DRACO
#define AFLOWRC_MPI_COMMAND_MPCDF_DRACO                       std::string("/usr/bin/srun -n") // MPCDF_DRACO_MPICH // FIX_DRACO
#define         MPI_COMMAND_MPCDF_DRACO                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_DRACO") // FIX_DRACO
#define AFLOWRC_MPI_NCPUS_MPCDF_DRACO                         0 // 32 // MPCDF_DRACO_MPICH
#define         MPI_NCPUS_MPCDF_DRACO                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_DRACO")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_DRACO                std::string("OFF")  // false/OFF, IGNORE/NEGLECT, true/ON
#define         MPI_HYPERTHREADING_MPCDF_DRACO                XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_DRACO") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_DRACO                    std::string("~/bin/") // MPCDF_DRACO_MPICH // FIX_DRACO
#define         MPI_BINARY_DIR_MPCDF_DRACO                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_DRACO") // FIX_DRACO

#define AFLOWRC_MPI_OPTIONS_MPCDF_COBRA                       std::string("ulimit -s unlimited ") // MPCDF_COBRA_MPICH  // FIX_COBRA
#define         MPI_OPTIONS_MPCDF_COBRA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_COBRA")  // FIX_COBRA
#define AFLOWRC_MPI_COMMAND_MPCDF_COBRA                       std::string("/usr/bin/srun -n") // MPCDF_COBRA_MPICH // FIX_COBRA
#define         MPI_COMMAND_MPCDF_COBRA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_COBRA") // FIX_COBRA
#define AFLOWRC_MPI_NCPUS_MPCDF_COBRA                         0 // 40 // MPCDF_COBRA_MPICH
#define         MPI_NCPUS_MPCDF_COBRA                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_COBRA")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_COBRA                std::string("OFF")  // false/OFF, IGNORE/NEGLECT, true/ON
#define         MPI_HYPERTHREADING_MPCDF_COBRA                XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_COBRA") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_COBRA                    std::string("~/bin/") // MPCDF_COBRA_MPICH // FIX_COBRA
#define         MPI_BINARY_DIR_MPCDF_COBRA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_COBRA") // FIX_COBRA

#define AFLOWRC_MPI_OPTIONS_MPCDF_HYDRA                       std::string("ulimit -s unlimited ") // MPCDF_HYDRA_MPICH
#define         MPI_OPTIONS_MPCDF_HYDRA                       XHOST.adefault.getattachedscheme("MPI_OPTIONS_MPCDF_HYDRA")
#define AFLOWRC_MPI_COMMAND_MPCDF_HYDRA                       std::string("poe ") // MPCDF_HYDRA_MPICH
#define         MPI_COMMAND_MPCDF_HYDRA                       XHOST.adefault.getattachedscheme("MPI_COMMAND_MPCDF_HYDRA")
#define AFLOWRC_MPI_NCPUS_MPCDF_HYDRA                         0 // 24 // MPCDF_HYDRA_MPICH
#define         MPI_NCPUS_MPCDF_HYDRA                         XHOST.adefault.getattachedutype<int>("MPI_NCPUS_MPCDF_HYDRA")
#define AFLOWRC_MPI_HYPERTHREADING_MPCDF_HYDRA                std::string("OFF")  // false/OFF, IGNORE/NEGLECT, true/ON
#define         MPI_HYPERTHREADING_MPCDF_HYDRA                XHOST.adefault.getattachedscheme("MPI_HYPERTHREADING_MPCDF_HYDRA") 
#define AFLOWRC_MPI_BINARY_DIR_MPCDF_HYDRA                    std::string("~/bin/") // MPCDF_HYDRA_MPICH
#define         MPI_BINARY_DIR_MPCDF_HYDRA                    XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MPCDF_HYDRA")

//define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    std::string("module purge") // FULTON_MARYLOU
#define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    std::string("export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE") // FULTON_MARYLOU
//#define AFLOWRC_MPI_OPTIONS_FULTON_MARYLOU                    std::string("export OMP_NUM_THREADS=$SLURM_NTASKS") // FULTON_MARYLOU
#define         MPI_OPTIONS_FULTON_MARYLOU                    XHOST.adefault.getattachedscheme("MPI_OPTIONS_FULTON_MARYLOU")
#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    std::string("srun ") // FULTON_MARYLOU
//#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    std::string("mpiexec") // FULTON_MARYLOU
//#define AFLOWRC_MPI_COMMAND_FULTON_MARYLOU                    std::string("mpiexec -np") // FULTON_MARYLOU WITH NP
#define         MPI_COMMAND_FULTON_MARYLOU                    XHOST.adefault.getattachedscheme("MPI_COMMAND_FULTON_MARYLOU")
#define AFLOWRC_MPI_BINARY_DIR_FULTON_MARYLOU                 std::string("/fslgroup/fslg_datamining/bin/") // FULTON_MARYLOU
#define         MPI_BINARY_DIR_FULTON_MARYLOU                 XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_FULTON_MARYLOU")

//DX - CMU EULER - START
#define AFLOWRC_MPI_OPTIONS_CMU_EULER                         std::string("") // CMU EULER
#define         MPI_OPTIONS_CMU_EULER                         XHOST.adefault.getattachedscheme("MPI_OPTIONS_CMU_EULER")
#define AFLOWRC_MPI_COMMAND_CMU_EULER                         std::string("mpirun -np") // CMU_EULER
#define         MPI_COMMAND_CMU_EULER                         XHOST.adefault.getattachedscheme("MPI_COMMAND_CMU_EULER")
#define AFLOWRC_MPI_BINARY_DIR_CMU_EULER                      std::string("/home/Tools/bin/") // CMU_EULER
#define         MPI_BINARY_DIR_CMU_EULER                      XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_CMU_EULER")
//DX - CMU EULER - END

#define AFLOWRC_MPI_OPTIONS_MACHINE1                          std::string("") // future expansions
#define         MPI_OPTIONS_MACHINE1                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE1")
#define AFLOWRC_MPI_COMMAND_MACHINE1                          std::string("...something ...")  // future expansions
#define         MPI_COMMAND_MACHINE1                          XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE1")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE1                       std::string("/somewhere/")  // future expansions
#define         MPI_BINARY_DIR_MACHINE1                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE1")

#define AFLOWRC_MPI_OPTIONS_MACHINE2                          std::string("") // future expansions
#define         MPI_OPTIONS_MACHINE2                          XHOST.adefault.getattachedscheme("MPI_OPTIONS_MACHINE2")
#define AFLOWRC_MPI_COMMAND_MACHINE2                          std::string("stub not used")  // future expansions
#define         MPI_COMMAND_MACHINE2                          XHOST.adefault.getattachedscheme("MPI_COMMAND_MACHINE2")
#define AFLOWRC_MPI_BINARY_DIR_MACHINE2                       std::string("/home/aflow/bin/")  // future expansions
#define         MPI_BINARY_DIR_MACHINE2                       XHOST.adefault.getattachedscheme("MPI_BINARY_DIR_MACHINE2")

// --------------------------------------------------------------------------
// aflow_aflowrc.cpp
namespace aflowrc {
  bool is_available(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool read(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool write_default(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool print_aflowrc(std::ostream& oss,bool AFLOWRC_VERBOSE);
} // namespace aflowrc


#endif // _AFLOW_AFLOWRC_H_

// POCC STUFF
// defaults go in 4 positions; Here with #define AFLOWRC_DEFAULT, in read(), in write_default(), and in print_aflowrc()...
// I coded strings (without spaces), <int>.. you can do <doubles> just like the <int>
// for strings with spaces I need to fix the code. Dont add them now. Then you need to go around the whole code and fix the use of DEFAULTS, possibly also in the READMEs if they are specified.
// STRING     string blablabla=DEFAULT_STRING =>  string blablabla=XHOST.adefault.getattachedscheme("DEFAULT_STRING")
// INT        int blablabla=DEFAULT_INT       =>  int blablabla=XHOST.adefault.getattachedscheme<int>("DEFAULT_INT")
// DOUBLE     double blablabla=DEFAULT_DOUBLE =>  double blablabla=XHOST.adefault.getattachedscheme<double>("DEFAULT_DOUBLE")
// ./aflow --machine to check them out
