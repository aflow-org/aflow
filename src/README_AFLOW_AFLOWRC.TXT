****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
LATEST VERSION OF THE FILE:                     materials.duke.edu/AFLOW/README_AFLOW_AFLOWRC.TXT
****************************************************************************************************

AFLOWRC README:

README written by: Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
GENERAL OVERVIEW

The .aflow.rc file is a general configuration file for AFLOW. 
Included in the file are machine-specific paths (VASP executables, pseudopotential files, 
database locations) and general run defaults.
Many of these settings can be overridden by tags in the aflow.in, and even by command line
flags when running AFLOW.
The general priority is as follows: command-line flags, aflow.in tags, and .aflow.rc settings.
The .aflow.rc file is located in your $HOME directory (~/.aflow.rc).
In the event that the file does not exist, AFLOW will generate a new .aflow.rc.
It is a good practice to keep the .aflow.rc file up to date with every new version.
Consider moving ~/.aflow.rc to ~/.aflow.rc_OLD and melding your system modifications.

****************************************************************************************************
LIST OF AFLOWRC SETTINGS:

AFLOWRC
  description: Indicates the version of AFLOW used to generate the .aflow.rc.
  type: string
  default: "X.X.XXX"

####################################################################################################
DEFAULT DEFINITIONS

DEFAULT_KZIP_BIN
  description: Compression format executable for files after the run has completed.
  related aflow.in tag: [AFLOW_MODE_ZIP=xz]
  type: string
  default: "xz"

DEFAULT_KZIP_EXT
  description: Compression format extension for files after the run has completed.
  type: string
  default: ".xz"

DEFAULT_TMPFS_DIRECTORIES
  description: Path(s) to search for the temporary directory 
  type: string
  default: "/tmp/,/run/shm/,/dev/shm/"

####################################################################################################
DEFAULTS FOR THE ENTRY LOADER

DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE
  description: location of the public alloy SQLITE DB file
  type: string
  default: "~/.aflow/aflowlib_alloy.db"

DEFAULT_ENTRY_LOADER_AFLUX_SERVER
  description: server that provides the AFLUX API
  type: string
  default: "aflowlib.duke.edu"

DEFAULT_ENTRY_LOADER_AFLUX_PATH
  description: base path including AFLUX API version
  type: string
  default: "/API/aflux/v1.0/"

DEFAULT_ENTRY_LOADER_RESTAPI_SERVER
  description: server that provides the AFLOW RESTAPI
  type: string
  default: "aflowlib.duke.edu"

DEFAULT_ENTRY_LOADER_RESTAPI_PATH
  description: AFLOW RESTAPI base path
  type: string
  default: "/AFLOWDATA/"

DEFAULT_ENTRY_LOADER_FS_PATH
  description: base path to the internal filesystem
  type: string
  default: "/common/"

####################################################################################################
DEFAULTS FOR THE AFLOW DATABASE

DEFAULT_AFLOW_DB_FILE
  description: The location of the AFLOW database file.
  type: string
  default: "/mnt/ramd/aflowlib.db"

DEFAULT_AFLOW_DB_STATS_FILE
  description: The location of the AFLOW database statistics file.
  type: string
  default: "/mnt/ramd/aflowlib.json"

DEFAULT_AFLOW_DB_DATA_PATH
  description: The location of the AFLOW database JSON files.
  type: string
  default: "/common/AFLOW/LIBS/"

DEFAULT_AFLOW_DB_LOCK_FILE
  description: The location of the AFLOW database lock file.
  type: string
  default: "/mnt/ramd/ADB_Idle.lock"

DEFAULT_AFLOW_DB_STALE_THRESHOLD
  description: The time in seconds at which the AFLOW database rebuild process is considered stale.
  type: double
  default: 10800 (3 hours)

####################################################################################################
FILENAMES FOR AFLOW.ORG ANALYSIS

DEFAULT_FILE_AFLOWLIB_ENTRY_OUT
  description: File name for the entry summary properties in plain text (--lib2raw).
  type: string
  default: "aflowlib.out"

DEFAULT_FILE_AFLOWLIB_ENTRY_JSON
  description: File name for the entry summary properties in JSON format (--lib2raw).  
  type: string
  default: "aflowlib.json"

DEFAULT_FILE_EDATA_ORIG_OUT
  description: File name for the extended crystallographic data of the original structure in plain text (--edata).
  type: string
  default: "edata.orig.out"

DEFAULT_FILE_EDATA_RELAX_OUT
  description: File name for the extended crystallographic data of the relaxed structure in plain text (--edata).
  type: string
  default: "edata.relax.out"

DEFAULT_FILE_EDATA_BANDS_OUT
  description: File name for the extended crystallographic data of the primitivized relaxed structure in plain text (--edata).
  type: string
  default: "edata.bands.out"

DEFAULT_FILE_DATA_ORIG_OUT
  description: File name for the crystallographic data of the original structure in plain text (--data).
  type: string
  default: "data.orig.out"

DEFAULT_FILE_DATA_RELAX_OUT
  description: File name for the crystallographic data of the relaxed structure in plain text (--data).
  type: string
  default: "data.relax.out"

DEFAULT_FILE_DATA_BANDS_OUT
  description: File name for the crystallographic data of the primitivized relaxed structure in plain text (--data).
  type: string
  default: "data.bands.out"

DEFAULT_FILE_EDATA_ORIG_JSON
  description: File name for the extended crystallographic data of the original structure in JSON format (--edata --print=json).
  type: string
  default: "edata.orig.json"

DEFAULT_FILE_EDATA_RELAX_JSON
  description: File name for the extended crystallographic data of the relaxed structure in JSON format (--edata --print=json).
  type: string
  default: "edata.relax.json"

DEFAULT_FILE_EDATA_BANDS_JSON
  description: File name for the extended crystallographic data of the primitivized relaxed structure in JSON format (--edata --print=json).
  type: string
  default: "edata.bands.json"

DEFAULT_FILE_DATA_ORIG_JSON
  description: File name for the crystallographic data of the original structure in plain text (--data --print=json).
  type: string
  default: "data.orig.json"

DEFAULT_FILE_DATA_RELAX_JSON
  description: File name for the crystallographic data of the relaxed structure in plain text (--data --print=json).
  type: string
  default: "data.relax.json"

DEFAULT_FILE_DATA_BANDS_JSON
  description: File name for the crystallographic data of the primitivized relaxed structure in plain text (--data --print=json).
  type: string
  default: "data.bands.json"

DEFAULT_FILE_TIME_OUT
  description: File name for the entry calculation time (--lib2raw).
  type: string
  default: "time"

DEFAULT_FILE_SPACEGROUP1_OUT
  description: File name for the entry space group pre/mid/post relaxation (--lib2raw --aflowSG=loose).
  type: string
  default: "SpaceGroup"

DEFAULT_FILE_SPACEGROUP2_OUT
  description: File name for the entry space group pre/mid/post relaxation (--lib2raw --aflowSG=tight).
  type: string
  default: "SpaceGroup2"

DEFAULT_FILE_VOLDISTPARAMS_OUT
  description: File name for the entry volume/nearest-neighbor-distances post relaxation (--lib2raw).
  type: string
  default: "VOLDISTParams"

DEFAULT_FILE_VOLDISTEVOLUTION_OUT
  description: File name for the entry volume/nearest-neighbor-distances pre/mid/post relaxation (--lib2raw).
  type: string
  default: "VOLDISTEvolution"

####################################################################################################
FILENAMES FOR AFLOW OPERATION

DEFAULT_AFLOW_PRESCRIPT_OUT
  description: File name for the output of the script to be run before AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.prescript.out"

DEFAULT_AFLOW_PRESCRIPT_COMMAND
  description: Script to be run before AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.prescript.command"

DEFAULT_AFLOW_POSTSCRIPT_OUT
  description: File name for the output of the script to be run after AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.postscript.out"

DEFAULT_AFLOW_POSTSCRIPT_COMMAND
  description: Script to be run after AFLOW_MODE_BINARY simulations.
  type: string
  default: "aflow.postscript.command"

DEFAULT_AFLOW_PGROUP_OUT
  description: File name for the point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroup.out"

DEFAULT_AFLOW_PGROUP_JSON
  description: File name for the point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroup.json"

DEFAULT_AFLOW_PGROUP_XTAL_OUT
  description: File name for the crystallographic point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroup_xtal.out"

DEFAULT_AFLOW_PGROUP_XTAL_JSON
  description: File name for the crystallographic point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroup_xtal.json"

DEFAULT_AFLOW_PGROUPK_OUT
  description: File name for the reciprocal point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroupk.out"

DEFAULT_AFLOW_PGROUPK_JSON
  description: File name for the reciprocal point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroupk.json"

DEFAULT_AFLOW_PGROUPK_XTAL_OUT
  description: File name for the dual of the crystallographic point group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.pgroupk_xtal.out"

DEFAULT_AFLOW_PGROUPK_XTAL_JSON
  description: File name for the dual of the crystallographic point group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.pgroupk_xtal.json"

DEFAULT_AFLOW_FGROUP_OUT
  description: File name for the factor group (coset representative) operations in plain text (--aflow-sym).
  type: string
  default: "aflow.fgroup.out"

DEFAULT_AFLOW_FGROUP_JSON
  description: File name for the factor group (coset representative) operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.fgroup.json"

DEFAULT_AFLOW_SGROUP_OUT
  description: File name for the full space group operations in plain text (--aflow-sym).
  type: string
  default: "aflow.sgroup.out"

DEFAULT_AFLOW_SGROUP_JSON
  description: File name for the full space group operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.sgroup.json"

DEFAULT_AFLOW_AGROUP_OUT
  description: File name for the site symmetry operations in plain text (--aflow-sym).
  type: string
  default: "aflow.agroup.out"

DEFAULT_AFLOW_AGROUP_JSON
  description: File name for the site symmetry operations in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.agroup.json"

DEFAULT_AFLOW_IATOMS_OUT
  description: File name for the inequivalent atoms calculation in plain text (--aflow-sym).
  type: string
  default: "aflow.iatoms.out"

DEFAULT_AFLOW_IATOMS_JSON
  description: File name for the inequivalent atoms calculation in JSON format (--aflow-sym --print=json).
  type: string
  default: "aflow.iatoms.json"

DEFAULT_AFLOW_ICAGES_OUT
  description: File name for the interstitial defects (cages) calculation in plain text (--cages).
  type: string
  default: "aflow.icages.out"

DEFAULT_AFLOW_SURFACE_OUT
  description: File name for the surface (hkl) calculation in plain text (--hkl_search).
  type: string
  default: "aflow.surface.out"

DEFAULT_AFLOW_QMVASP_OUT
  description: File name for the VASP calculation summary in plain text (--qmvasp).
  type: string
  default: "aflow.qmvasp.out"

DEFAULT_AFLOW_ERVASP_OUT
  description: File name for the VASP calculation errors in plain text.
  type: string
  default: "aflow.error.out"

DEFAULT_AFLOW_IMMISCIBILITY_OUT
  description: File name for the miscibility determination in plain text (--miscibility).
  type: string
  default: "aflow.immiscibility.out"

DEFAULT_AFLOW_MEMORY_OUT
  description: File name for the memory check in plain text.
  type: string
  default: "aflow.memory.out"

DEFAULT_AFLOW_FROZSL_INPUT_OUT
  description: File name for the input settings of the FROZSL calculation.
  type: string
  default: "aflow.frozsl_input.out"

DEFAULT_AFLOW_FROZSL_POSCAR_OUT
  description: File name for the structure used in the FROZSL calculation.
  type: string
  default: "aflow.frozsl_poscar.out"

DEFAULT_AFLOW_FROZSL_MODES_OUT
  description: File name for the phonon eigen-energies from the FROZSL calculation.
  type: string
  default: "aflow.frozsl_energies.out"

DEFAULT_AFLOW_FROZSL_EIGEN_OUT
  description: File name for the phonon eigen-modes from the FROZSL calculation.
  type: string
  default: "aflow.frozsl_eigen.out"

DEFAULT_AFLOW_END_OUT
  description: File name for the output file indicating AFLOW is finished.
  type: string
  default: "aflow.end.out"

DEFAULT_AFLOW_PLASMONICS_OUT
  description: File name for the plasmonics data output file.
  type: string
  default: "aflow.plasmonics_eps_0.01.out"

####################################################################################################
DEFAULT GENERIC MPI

MPI_START_DEFAULT
  description: Commands to run before running the MPI program.
  type: string
  default: ""

MPI_STOP_DEFAULT
  description: Commands to run after running the MPI program.
  type: string
  default: ""

MPI_COMMAND_DEFAULT
  description: MPI program command.
  type: string
  default: "mpirun -np"

MPI_NCPUS_DEFAULT
  description: Number of cores to run with the MPI program.
  type: integer
  default: 4

MPI_NCPUS_MAX
  description: Maximum number of cores that AFLOW can consume.
  type: integer
  default: 64

####################################################################################################
DEFAULTS BINARY

DEFAULT_VASP_GAMMA_BIN
  description: Name for serial Gamma-Point-only VASP executable.
  type: string
  default: "vasp46s_gamma"

DEFAULT_VASP_GAMMA_MPI_BIN
  description: Name for MPI Gamma-Point-only VASP executable.
  type: string
  default: "mpivasp46s_gamma"

DEFAULT_VASP_BIN
  description: Name for serial VASP executable.
  type: string
  default: "vasp46s"

DEFAULT_VASP_MPI_BIN
  description: Name for MPI VASP executable.
  type: string
  default: "mpivasp46s"

DEFAULT_VASP5_BIN
  description: Name for serial VASP 5 executable.
  type: string
  default: "vasp54s"

DEFAULT_VASP5_MPI_BIN
  description: Name for MPI VASP 5 executable.
  type: string
  default: "mpivasp54s"

DEFAULT_AIMS_BIN
  description: Name for FHI-AIMS executable.
  type: string
  default: "aims"

####################################################################################################
DEFAULTS POTCARS

DEFAULT_VASP_POTCAR_DIRECTORIES
  description: Path(s) to search for the VASP pseudopotential files, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "/common/VASP,/common/AFLOW/VASP,/home/aflow/common/AFLOW/VASP,/fslhome/fslcollab8/group/VASP,/fslhome/glh43/src/,/share/home/00470/tg457283/common/AFLOW/VASP/,/share/home/00457/tg457357/common/AFLOW/VASP/,/home/mehl/bin/AFLOW/VASP/,~/common/VASP/,~/common/AFLOW/VASP/,/nics/a/proj/aflow/common/AFLOW/VASP/,/home/users/aflow/common/VASP,/share/apps/AFLOW3/VASP,/projects/kyang-group/common/VASP,/somewhere/"

DEFAULT_VASP_POTCAR_DATE
  description: Date component of VASP pseudopotential files path, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "current"

DEFAULT_VASP_POTCAR_SUFFIX
  description: Last component of VASP pseudopotential files path, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "/POTCAR"

DEFAULT_VASP_POTCAR_DATE_POT_LDA
  description: Date component of VASP LDA pseudopotential files to be added to the --aflow_proto directory, e.g., V_pv:LDA:01Apr2000.
  type: string
  default: "01Apr2000"

DEFAULT_VASP_POTCAR_DATE_POT_GGA
  description: Date component of VASP GGA pseudopotential files to be added to the --aflow_proto directory, e.g., V_pv:LDA:01Apr2000.
  type: string
  default: "01Apr2000"

DEFAULT_VASP_POTCAR_DIR_POT_LDA
  description: Path component of VASP LDA pseudopotential files, e.g., /common/VASP/pot_LDA/current/Mn_pv/POTCAR.
  type: string
  default: "pot_LDA"

DEFAULT_VASP_POTCAR_DIR_POT_GGA
  description: Path component of VASP GGA pseudopotential files, e.g., /common/VASP/pot_GGA/current/Mn/POTCAR.
  type: string
  default: "pot_GGA"

DEFAULT_VASP_POTCAR_DIR_POT_PBE
  description: Path component of VASP PBE pseudopotential files, e.g., /common/VASP/pot_PBE/current/Mn/POTCAR.
  type: string
  default: "pot_PBE"

DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA
  description: Path component of VASP PAW LDA pseudopotential files, e.g., /common/VASP/potpaw_LDA/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_LDA"

DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA
  description: Path component of VASP PAW GGA pseudopotential files, e.g., /common/VASP/potpaw_GGA/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_GGA"

DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE
  description: Path component of VASP PAW PBE pseudopotential files, e.g., /common/VASP/potpaw_PBE/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_PBE"

DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN
  description: Path component of VASP PAW LDA pseudopotential files containing information on the kinetic energy density of the core-electrons, e.g., /common/VASP/potpaw_PBE.54/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_LDA.54"

DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN
  description: Path component of VASP PAW PBE pseudopotential files containing information on the kinetic energy density of the core-electrons, e.g., /common/VASP/potpaw_PBE.54/current/Mn_pv/POTCAR.
  type: string
  default: "potpaw_PBE.54"

####################################################################################################
DEFAULTS KPOINTS/DOS

DEFAULT_BANDS_GRID
  description: Number of points in which the energy dispersion is resolved between each leg of the high-symmetry path.
  type: integer
  default: 20

DEFAULT_BANDS_LATTICE
  description: Crystal lattice of the structure specifying the high-symmetry path.
  type: string
  choices: "AUTO", "BCC", "BCT1", "BCT2", "CUB", "FCC", "HEX", "MCL", "MCLC1", "MCLC2", "MCLC3", "MCLC4", "MCLC5", "ORC", "ORCC", "ORCF1", "ORCF2", "ORCF3", "ORCI", "RHL1", "RHL2", "TET", "TRI1a", "TRI1b", "TRI2a", "TRI2b"
  default: "AUTO"
  doi: 10.1016/j.commatsci.2010.05.010

DEFAULT_KSCHEME
  description: K-points mesh scheme.
  type: string
  choices: "M", "G"
  default: "M"

DEFAULT_KPPRA
  description: Number of k-points per reciprocal atom.
  type: integer
  default: 6000

DEFAULT_STATIC_KSCHEME
  description: K-points mesh scheme for the static (energy/DOS) run.
  type: string
  choices: "M", "G"
  default: "M"

DEFAULT_KPPRA_STATIC
  description: Number of k-points per reciprocal atom for the static (energy/DOS) run.
  type: integer
  default: 10000

DEFAULT_KPPRA_ICSD
  description: Number of k-points per reciprocal atom for structures from the ICSD catalog.
  type: integer
  default: 8000

DEFAULT_UNARY_BANDS_GRID
  description: Number of points in which the energy dispersion is resolved between each leg of the high-symmetry path for unary structures.
  type: integer
  default: 128

DEFAULT_UNARY_KPPRA
  description: Number of k-points per reciprocal atom for unary structures.
  type: integer
  default: 8000

DEFAULT_UNARY_KPPRA_STATIC
  description: Number of k-points per reciprocal atom for unary structures and the static (energy/DOS) run..
  type: integer
  default: 8000

DEFAULT_PHONONS_KSCHEME
  description: K-points mesh scheme for phonon calculations.
  type: string
  choices: "M", "G"
  default: "G"

DEFAULT_PHONONS_KPPRA
  description: Number of k-points per reciprocal atom for phonon calculations.
  type: integer
  default: 3000

DEFAULT_DOS_EMIN
  description: Minimum of the energy range over which the DOS is plotted.
  type: integer
  default: -10

DEFAULT_DOS_EMAX
  description: Maximum of the energy range over which the DOS is plotted.
  type: integer
  default: 10

DEFAULT_DOS_SCALE
  description: Interval for the energy axis over which the DOS is plotted.
  type: double
  default: 1.2

####################################################################################################
DEFAULTS PRECISION

DEFAULT_VASP_PREC_ENMAX_LOW
  description: Low-level (default) multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1

DEFAULT_VASP_PREC_ENMAX_MEDIUM
  description: Medium-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.3

DEFAULT_VASP_PREC_ENMAX_NORMAL
  description: Normal-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.3

DEFAULT_VASP_PREC_ENMAX_HIGH
  description: High-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.4

DEFAULT_VASP_PREC_ENMAX_ACCURATE
  description: Accurate-level multiplier for the max energy cutoff of the pseudopotentials.
  type: double
  default: 1.4

DEFAULT_VASP_ENMAX_MINIMUM
  description: Minimum multiplier for the max energy cutoff of the pseudopotentials (if PSMAXN error is encountered).
  type: double
  default: 0.25

DEFAULT_VASP_SPIN_REMOVE_CUTOFF
  description: Magnetic moment threshold below which the spin-polarization setting will be turned off.
  type: double
  default: 0.05

DEFAULT_VASP_PREC_POTIM
  description: Sets starting value for POTIM in VASP (scaling constant for the forces). AFLOW may adjust as errors occur.
  type: double
  default: 0.5

DEFAULT_VASP_PREC_EDIFFG
  description: Sets EDIFFG for VASP. Note that a positive value defines the break condition for the ionic relaxation loop, while a negative value defines the break condition based on the forces.
  type: double
  default: -0.001

####################################################################################################
DEFAULTS OPTIONS

DEFAULT_VASP_OUT
  description: Filename for vasp output.
  type: string
  default: vasp.out

DEFAULT_VASP_EXTERNAL_INCAR
  description: File to load in as default VASP INCAR.
  type: string
  default: ./INCAR

DEFAULT_VASP_EXTERNAL_POSCAR
  description: File to load in as default VASP POSCAR.
  type: string
  default: ./POSCAR

DEFAULT_VASP_EXTERNAL_POTCAR
  description: File to load in as default VASP POTCAR.
  type: string
  default: ./POTCAR

DEFAULT_VASP_EXTERNAL_KPOINTS
  description: File to load in as default VASP KPOINTS.
  type: string
  default: ./KPOINTS

DEFAULT_AIMS_EXTERNAL_CONTROL
  description: File to load in as default AIMS control.in.
  type: string
  default: ./control.in

DEFAULT_AIMS_EXTERNAL_GEOM
  description: File to load in as default AIMS geom.in.
  type: string
  default: ./geometry.in

DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE
  description: Default pseudopotential type.
  type: string
  default: potpaw_PBE

DEFAULT_VASP_FORCE_OPTION_RELAX_MODE_SCHEME
  description: Default relax mode setting.
  type: string
  choices: ENERGY, FORCES, ENERGY_FORCES, FORCES_ENERGY
  default: ENERGY

DEFAULT_VASP_FORCE_OPTION_RELAX_COUNT
  description: Default number of relaxations.
  type: integer
  default: 2

DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME
  description: Default VASP PREC setting.
  type: string
  choices: LOW, MEDIUM, NORMAL, HIGH, ACCURATE
  default: ACCURATE

DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME
  description: Default VASP ALGO setting.
  type: string
  choices: NORMAL, VERYFAST, FAST, ALL, DAMPED
  default: NORMAL

DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME
  description: Default VASP METAGGA setting.
  type: string
  choices: NONE, TPSS, RTPSS, M06L, MS0, MS1 MS2, MBJ, SCAN
  default: NONE

DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME
  description: Default VASP IVDW setting.
  type: bool
  choices: 0, 1, 10, 11, 12, 2, 20, 21, 202, 4
  default: 0

DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME
  description: Default VASP integration scheme.
  type: string
  choices: DEFAULT, METAL, INSULATOR, SEMICONDUCTOR
  default: DEFAULT

DEFAULT_VASP_FORCE_OPTION_ISMEAR_SCHEME
  description: Default VASP ISMEAR scheme.
  type: int
  default: 1

DEFAULT_VASP_FORCE_OPTION_ISMEAR_STATIC_SCHEME
  description: Default VASP ISMEAR scheme for STATIC runs.
  type: int
  default: -5

DEFAULT_VASP_FORCE_OPTION_ISMEAR_BANDS_SCHEME
  description: Default VASP ISMEAR scheme for BANDS runs.
  type: int
  default: 0

DEFAULT_VASP_FORCE_OPTION_SIGMA
  description: Default VASP SIGMA value.
  type: double
  default: 0.1

DEFAULT_VASP_FORCE_OPTION_SIGMA_STATIC
  description: Default VASP SIGMA value for STATIC runs.
  type: double
  default: 0.05

DEFAULT_VASP_FORCE_OPTION_SIGMA_BANDS
  description: Default VASP SIGMA value for BANDS runs.
  type: double
  default: 0.05

DEFAULT_VASP_NELM
  description: Sets the maximum number of electronic self-consistency steps (NELM) for VASP runs.
  type: int
  default: 60

DEFAULT_VASP_NELM_STATIC
  description: Sets the maximum number of electronic self-consistency steps (NELM) for VASP STATIC runs.
  type: int
  default: 120

MAX_VASP_NELM
  description: Sets the maximum number of electronic self-consistency steps (NELM) for VASP should DEFAULT_VASP_FORCE_OPTION_NELM or DEFAULT_VASP_FORCE_OPTION_NELM_STATIC not be enough for convergence.
  type: int
  default: 300

DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME
  description: Default VASP Kerker mixing scheme settings. Defines AMIX/BMIX/AMIX_MAG/BMIX_MAG.
  type: string
  choices: DEFAULT (unspecified), AUTO, US, PAW
  default: DEFAULT

DEFAULT_VASP_FORCE_OPTION_SYM
  description: Default VASP ISYM option.
  type: integer
  choices: 0, 1, 2
  default: 1

DEFAULT_VASP_FORCE_OPTION_SPIN
  description: Adapts INCAR for VASP spin-polarized calculations.
  type: bool
  default: 1

DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1
  description: Turns off spin-polarized calculation after first relaxation if the magnetic moment is less than VASP_SPIN_REMOVE_CUTOFF.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2
  description: Turns off spin-polarized calculation after second relaxation if the magnetic moment is less than VASP_SPIN_REMOVE_CUTOFF.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_BADER
  description: Adapts INCAR for VASP to retrieve core charge (AECCARx) for Bader partial charge calculation. Turning this setting on will write these files out for any calculation type (RELAX, STATIC, STATIC_BANDS, etc.).
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_BADER_STATIC
  description: Adapts INCAR for VASP to retrieve core charge (AECCARx) for Bader partial charge calculation. Turning this setting on will write these files out if the calculation type includes STATIC (RELAX_STATIC, RELAX_STATIC_BANDS, STATIC_BANDS, etc.).
  type: bool
  default: 1

DEFAULT_VASP_FORCE_OPTION_ELF
  description: Adapts INCAR for VASP to calculate the electron localization function.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM
  description: Adapts INCAR for VASP to initialize magnetic moment for all species. It defaults to a ferromagnetic solution.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_WAVECAR
  description: Adapts INCAR for VASP to save the WAVECAR.
  type: bool
  default: 0

DEFAULT_VASP_FORCE_OPTION_CHGCAR
  description: Adapts INCAR for VASP to save the CHGCAR.
  type: bool
  default: 1

DEFAULT_VASP_FORCE_OPTION_LSCOUPLING
  description: Adapts INCAR for VASP spin-orbit coupling calculations.
  type: bool
  default: 0

####################################################################################################
AFLOW_LIBRARY AFLOW_PROJECT

DEFAULT_AFLOW_LIBRARY_DIRECTORIES
  description: Path(s) to search for the library of materials calculations.
  type: string
  default: "/common/AFLOW/LIBS/,/home/aflow/common/AFLOW/LIBS/,/fslhome/glh43/src/,/usr/local/bin/,/fslhome/fslcollab8/group/bin/,/home/auro/work/AFLOW3/,~/common/AFLOW/LIBS/,./,/nics/a/proj/aflow/common/AFLOW/LIBS/,/home/users/aflow/common/AFLOW/LIBS,/home/junkai/PROTO_DATABASE/,/projects/kyang-group/common/LIBS,/somewhere/"

DEFAULT_AFLOW_PROJECTS_DIRECTORIES
  description: Path(s) of the various library projects.
  type: string
  default: "/common/AUID,/common/ICSD,/common/LIB0,/common/LIB1,/common/LIB2,/common/LIB3,/common/LIB4,/common/LIB5,/common/LIB6,/common/LIB7,/common/LIB8,/common/LIB9"

####################################################################################################
DEFAULT PLATON/FINDSYM

DEFAULT_PLATON_P_EQUAL
  description: Sets PLATON EQUAL tag (incorporate atom types).
  type: bool
  default: 0

DEFAULT_PLATON_P_EXACT
  description: Sets PLATON EXACT tag (match all atoms).
  type: bool
  default: 0

DEFAULT_PLATON_P_ANG
  description: Sets PLATON ang (angle criterion) for metrical symmetry of lattice.
  type: bool
  default: 1

DEFAULT_PLATON_P_D1
  description: Sets PLATON d1 (distance criterion) for non-inversion symmetry.
  type: double
  default: 0.25

DEFAULT_PLATON_P_D2
  description: Sets PLATON d2 (distance criterion) for inversion symmetry.
  type: double
  default: 0.25

DEFAULT_PLATON_P_D3
  description: Sets PLATON d3 (distance criterion) for translation symmetry.
  type: double
  default: 0.25

DEFAULT_FINDSYM_TOL
  description: Sets FINDSYM tolerance.
  type: double
  default: 0.001

####################################################################################################
DEFAULTS GNUPLOT

DEFAULT_GNUPLOT_EPS_FONT
  description: Sets gnuplot eps terminal font.
  type: string
  default: "Helvetica"

DEFAULT_GNUPLOT_EPS_FONT_BOLD
  description: Sets gnuplot eps terminal bold font.
  type: string
  default: "Helvetica-Bold"

DEFAULT_GNUPLOT_EPS_FONT_ITALICS
  description: Sets gnuplot eps terminal italicized font.
  type: string
  default: "Helvetica-Oblique"

DEFAULT_GNUPLOT_EPS_FONT_BOLD_ITALICS
  description: Sets gnuplot eps terminal bold italicized font.
  type: string
  default: "Helvetica-BoldOblique"

DEFAULT_GNUPLOT_PNG_FONT
  description: Sets gnuplot png terminal font.
  type: string
  default: "Arial"

DEFAULT_GNUPLOT_PNG_FONT_BOLD
  description: Sets gnuplot png terminal bold font.
  type: string
  default: "Arial_Bold"

DEFAULT_GNUPLOT_PNG_FONT_ITALICS
  description: Sets gnuplot png terminal italicized font.
  type: string
  default: "Arial_Italic"

DEFAULT_GNUPLOT_PNG_FONT_BOLD_ITALICS
  description: Sets gnuplot png terminal bold italicized font.
  type: string
  default: "Arial_BoldItalic"

####################################################################################################
DEFAULTS CHULL

DEFAULT_CHULL_ALLOWED_DFT_TYPES
  description: Defines the allowed entries based on the DFT calculation type.
  type: string
  choices: "US", "GGA", "PAW_LDA", "PAW_GGA", "PAW_PBE", "PAW_PBE_KIN", "GW", "HSE06"
  default: "PAW_PBE,PAW_PBE_KIN"

DEFAULT_CHULL_ALLOW_ALL_FORMATION_ENERGIES
  description: Allows all entries independent of DFT calculation type.
  type: bool
  default: 0

DEFAULT_CHULL_COUNT_THRESHOLD_BINARIES
  description: Defines the minimum number of entries for a reliable binary hull.
  type: integer
  default: 200

DEFAULT_CHULL_PERFORM_OUTLIER_ANALYSIS
  description: Enables determination of outliers.
  type: bool
  default: 1

DEFAULT_CHULL_OUTLIER_ANALYSIS_COUNT_THRESHOLD_BINARIES
  description: Defines the minimum number of entries for a reliable outlier analysis. Only phases stable with respect to their end-members are considered for the outlier analysis (below the zero H_f tie-line).
  type: integer
  default: 75

DEFAULT_CHULL_OUTLIER_MULTIPLIER
  description: Defines the bounds beyond the interquartile range for which points are considered outliers.
  type: double
  default: 3.25

DEFAULT_CHULL_IGNORE_KNOWN_ILL_CONVERGED
  description: AFLOW maintains a list of (older) prototypes known to have converged poorly. These entries are likely outliers. If this flag is on (1), then the entries are removed automatically. Turning this flag off (0) is not recommended.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_BANNER
  description: Defines banner type.
  type: integer
  choices: 0 (no banner), 1 (normal banner), 2 (small banner)
  default: 1

DEFAULT_CHULL_LATEX_COMPOUNDS_COLUMN
  description: Adds an additional column to the PDF table output with the entry compound label.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_STOICH_HEADER
  description: Toggles between a reduced composition header (default) or a stoichiometry header for the entries table.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_PLOT_UNARIES
  description: Incorporates the end-members in the convex hull illustration.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_PLOT_OFF_HULL
  description: Incorporates unstable phases in the convex hull illustration, but excludes phases unstable with respect to their end-members (above the zero H_f tie-line).
  type: integer
  choices: -1 (default: true for 2-dimensional systems, false for 3-dimensional systems), 0 (false), 1 (true).
  default: -1

DEFAULT_CHULL_LATEX_PLOT_UNSTABLE
  description: Incorporates all unstable phases in the convex hull illustration.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_FILTER_SCHEME
  description: Defines exclusion scheme for the convex hull illustration. In contrast to --neglect, this scheme is limited only to the illustration, points are still included in the analysis/report. The following strings are accepted: Z-axis (also Energy-axis) or Distance. Z-axis refers to a scheme that eliminates structures from the illustration based on their formation enthalpy. On the other hand, Distance refers to a scheme that eliminates structures from the illustration based on their distances from the hull. The criteria (value) for elimination is defined by DEFAULT_CHULL_LATEX_FILTER_VALUE.
  type: string
  choices: "z-axis", "distance"
  default: ""

DEFAULT_CHULL_LATEX_FILTER_VALUE
  description: Defines the value beyond which points are excluded per the scheme defined with DEFAULT_CHULL_LATEX_FILTER_SCHEME.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_COLOR_BAR
  description: Defines whether to show the color bar graphic. Colors can still be incorporated without the color bar graphic.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_HEAT_MAP
  description: Defines whether to color the facets with heat maps illustrating their depth (3-dimensional illustration only).
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_COLOR_GRADIENT
  description: Defines whether to incorporate a color scheme at all in the illustration. Turning this flag off will also turn off DEFAULT_CHULL_LATEX_COLOR_BAR and DEFAULT_CHULL_LATEX_HEAT_MAP.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_COLOR_MAP
  description: Defines the color map, options are presented in the PGFPLOTS Package Manual. Default is "rgb(0pt)=(0.035,0.270,0.809); rgb(63pt)=(1,0.644,0)".
  type: string
  default: ""

DEFAULT_CHULL_LATEX_TERNARY_LABEL_COLOR
  description: Defines the color used for ternary compounds. Default is white.
  type: string
  default: ""

DEFAULT_CHULL_LATEX_REVERSE_AXIS
  description: Flips the axes of binary convex hulls.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_FACET_LINE_DROP_SHADOW
  description: Adds a drop shadow underneath facet lines.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_LINKS
  description: Defines the links scheme. True/false, i.e., 0/1, will toggle all links on/off. 2 enables external hyperlinks only (no links to other sections of the PDF). 3 enables internal links only (no links to external pages).
  type: bool
  choices: 0 (false), 1 (true), 2 (external-only), or 3 (internal-only)
  default: 1

DEFAULT_CHULL_LATEX_LABEL_NAME
  description: Defines the labeling scheme for phases shown on the convex hull. By default, the compound labels are shown, while the prototype label can also be specified. "icsd" shows the ICSD entry number designation (lowest for multiple equivalent ground-state structures reflecting icsd_canonical_auid) if appropriate. Also acceptable: both and none.
  type: string
  choices: "compound", "prototype", "both", "icsd", "none"
  default: ""

DEFAULT_CHULL_LATEX_META_LABELS
  description: Enables verbose labels, including compound, prototype, H_f, T_S, and Delta H_f. Warning, significant overlap of labels should be expected.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_LABELS_OFF_HULL
  description: Enables labels for unstable points.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_PLOT_REDUCE_COMPOSITION
  description: Dictates whether compound labels on the plot are reduced. Only three values are accepted: -1 (default: true unless a binary is requested and an explicit "compounds" option is set for DEFAULT_CHULL_LATEX_LABEL_NAME), 0 (false), and 1 (true).
  type: integer
  choices: -1 (default), 0 (false), 1 (true)
  default: -1

DEFAULT_CHULL_LATEX_HELVETICA_FONT
  description: Toggles the font scheme from Computer Modern (default) to Helvetica.
  type: bool
  default: 0

DEFAULT_CHULL_LATEX_FONT_SIZE
  description: Defines the font size of the labels on the convex hull illustration. Warning, other settings may override this default.
  type: string
  choices: "tiny", "scriptsize", "footnotesize", "small", "normalsize", "large" (default), "Large", "LARGE", "huge", "Huge".
  default: ""

DEFAULT_CHULL_LATEX_ROTATE_LABELS
  description: Toggles whether labels are appropriately rotated.
  type: bool
  default: 1

DEFAULT_CHULL_LATEX_BOLD_LABELS
  description: Toggles whether labels are bolded. Only three values are accepted: -1 (default: false unless phase is a ternary), 0 (false), 1 (true).
  type: integer
  default: -1
  choices -1 (default), 0 (false), or 1 (true)

DEFAULT_CHULL_PNG_RESOLUTION
  description: Default resolution for conversion pdf to png (convert via ImageMagick).
  type: integer
  default: 300

####################################################################################################
DEFAULTS CCE

DEFAULT_CCE_OX_METHOD
  description: Default method to determine oxidation numbers.
  type: int
  choices: 1 (ELECTRONEGATIVITY_ALLEN), 2 (BADER)
  default: 1

DEFAULT_CCE_NN_DIST_TOL_MULTI_ANION
  description: Tolerance between shortest and longest bonds when testing for multi-anion compound in Ang.
  type: double
  default: 0.4

DEFAULT_CCE_OX_TOL
  description: Tolerance for sum of oxidation numbers to be considered as zero.
  type: double
  default: 0.001

DEFAULT_CCE_SELF_DIST_TOL
  description: Distance tolerance in Ang for neighbor screening to savely exclude the central atom itself.
  type: double
  default: 0.001

DEFAULT_CCE_PEROX_CUTOFF
  description: Upper bound for O-O bond length in peroxides in Ang.
  type: double
  default: 1.6

DEFAULT_CCE_SUPEROX_CUTOFF
  description: Upper bound for O-O bond length in superoxides in Ang.
  type: double
  default: 1.4

DEFAULT_CCE_O2_MOLECULE_UPPER_CUTOFF
  description: Upper bound for O-O bond length for molecular O2 in Ang.
  type: double
  default: 1.3

DEFAULT_CCE_O2_MOLECULE_LOWER_CUTOFF
  description: Lower bound for O-O bond length for molecular O2 in Ang.
  type: double
  default: 1.2

####################################################################################################
DEFAULTS GFA

DEFAULT_GFA_FORMATION_ENTHALPY_CUTOFF
  description: Default cutoff formation enthalpy (in eV) for systems included in the GFA analysis.
  type: double
  default: 0.05

####################################################################################################
ARUN_DIRECTORY_PREFIX
  description: Default prefix for ARUN directories.
  type: string
  default: "ARUN."

####################################################################################################
DEFAULTS POCC

DEFAULT_POCC_TEMPERATURE_STRING
  description: Temperature range over which POCC will perform post-processing. Scan 0 to 1200 (K) in 300 (K) increments. 0 is converted to _ZERO_TOL_ since T=0 exactly is undefined for the Boltzmann distribution.
  type: string
  default: "0:1200:300"

DEFAULT_POCC_EXCLUDE_UNSTABLE
  description: Decided whether to exclude unstabale structures from the ensemble average
  type: bool
  default: 1

DEFAULT_POCC_SITE_TOL
  description: Default site tolerance for HNF calculation.
  type: double
  default: 0.001

DEFAULT_POCC_STOICH_TOL
  description: Default stoichiometry tolerance for HNF calculation.
  type: double
  default: 0.001

DEFAULT_UFF_BONDING_DISTANCE
  description: Default distance cutoff for bonded atoms.
  type: double
  default: 0.5 Angstroms

DEFAULT_UFF_ENERGY_TOLERANCE
  description: Default energy tolerance for equivalent structures.
  type: double
  default: 1e-6 eV

DEFAULT_UFF_CLUSTER RADIUS
  description: Default cluster radius for UFF energy determination.
  type: double
  default: 10 Angstroms

DEFAULT_POCC_PERFORM_ROBUST_STRUCTURE_COMPARISON
  description: Default option for performing robust structure comparison (AFLOW-XtalMatch) IN ADDITION TO UFF energy analysis.
  type: bool
  default: 0

DEFAULT_POCC_WRITE_OUT_ALL_DERIVATIVE_STRUCTURES
  description: Default option for writing out all derivative structures to file (named POCC_ALL_DERIVATIVE_STRUCTURES_FILE).
  type: bool
  default: 1

POCC_FILE_PREFIX
  description: Prefix for files written by POCC module.
  type: string
  default: "aflow.pocc."

POCC_ALL_DERIVATIVE_STRUCTURES_FILE
  description: Name for all derivative structures file.
  type: string
  default: "all_structures.out"

POCC_UNIQUE_DERIVATIVE_STRUCTURES_FILE
  description: Name for unique derivative structures file.
  type: string
  default: "unique_structures.out"

POCC_ALL_HNF_MATRICES_FILE
  description: Name for HNF matrices file.
  type: string
  default: "hnf_matrices.out"

POCC_ALL_SITE_CONFIGURATIONS_FILE
  description: Name for site configurations file.
  type: string
  default: "site_configurations.out"

POCC_OUT_FILE
  description: Name for POCC post-processing properties file.
  type: string
  default: "out"

POCC_APL_OUT_FILE
  description: Name for the POCC-APL post-processing properties file.
  type: string
  default: "apl.out"

POCC_DOSCAR_FILE
  description: Name for ensemble average DOSCAR.
  type: string
  default: "DOSCAR.pocc"

POCC_PHDOSCAR_FILE
  description: Name for ensemble average PHDOSCAR.
  type: string
  default: "PHDOSCAR.pocc"

POCC_ANIONS_LIST
  description: List of elements that could be considered anions.
  type: string
  default: "B,C,N,O"

####################################################################################################
DEFAULTS QCA

DEFAULT_QCA_MIN_SLEEP_SECONDS
  description: Minimum time to sleep between waiting for ATAT to finish.
  type: integer
  default: 60

DEFAULT_QCA_MAX_NUM_ATOMS
  description: Maximum number of atoms in the unit cell for derivative structures created by the ATAT structure generator.
  type: integer
  default: 8

DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS
  description: Maximum number of atoms in the unit cell for derivative structures found in the AFLOW calculations.
  type: integer
  default: 4

DEFAULT_QCA_CV_CUTOFF
  description: Cross-validation score cut-off (in eV) for the cluster expansion.
  type: double
  default: 0.05

DEFAULT_QCA_CONC_NPTS
  description: Number of points to sub-divide the concentration profile.
  type: integer
  default: 20

DEFAULT_QCA_TEMP_NPTS
  description: Number of points to sub-divide the temperature profile.
  type: integer
  default: 150

DEFAULT_QCA_TEMP_MIN
  description: Minimum temperature (in K) of the temperature profile.
  type: double
  default: 300

DEFAULT_QCA_TEMP_MAX
  description: Maximum temperature (in K) of the temperature profile.
  type: double
  default: 5000

DEFAULT_QCA_TEMP_MIN_LIMIT
  description: Temperature limit (in K) at which to end the calculation if a solution cannot be found.
  type: double
  default: 10000

DEFAULT_QCA_FORMAT_DATA
  description: Format of the output data.
  type: string
  default: "txt"

DEFAULT_QCA_FORMAT_PLOT
  description: Format of the plotted image.
  type: string
  default: "png"

####################################################################################################
DEFAULTS APL

DEFAULT_APL_PREC
  description: PREC value for phonon calculations, including the relaxation calculation.
  type: string
  default: PHONONS

DEFAULT_APL_ENGINE
  description: Sets whether to use the direct (DM) or linear response (LR) method for phonon calculations.
  type: string
  choices: "DM", "LR"
  default: "DM"

DEFAULT_APL_HIBERNATE
  description: Use the hibernation functionality.
  type: bool
  default: 1

DEFAULT_APL_MINATOMS
  description: Minimum number of atoms to use for the supercell generation.
  type: int
  default: 100

DEFAULT_APL_MINSHELL
  description: Minimum number of coordination shells to use for the supercell generation.
  type: int
  default: 6

DEFAULT_APL_POLAR
  description: Use polar corrections for phonon calculations.
  type: bool
  default: 1

DEFAULT_APL_DMAG
  description: The distortion magnitude for the direct method. Also used for AAPL.
  type: double
  default: 0.015

DEFAULT_APL_DXYZONLY
  description: Use distortions along the crystallographic axes only.
  type: bool
  default: 0

DEFAULT_APL_DSYMMETRIZE
  description: Use site symmetry to reduce the number of distortions (VASP calculations) needed.
  type: bool
  default: 1

DEFAULT_APL_DPM
  description: Defines whether the negative distortion directions should be taken into account for phonon calculations.
  type: string
  choices: "ON", "OFF", "AUTO"
  default: "AUTO"

DEFAULT_APL_RELAX
  description: Do high-precision structural relaxation before phonon calculations.
  type: bool
  default: 1

DEFAULT_APL_ZEROSTATE
  description: Calculate the forces of the undistorted supercell.
  type: bool
  default: 1

DEFAULT_APL_ZEROSTATE_CHGCAR
  description: Use the CHGCAR of the ZEROSTATE calculation to be read by other APL and AAPL calculations.
  type: bool
  default: 0

DEFAULT_APL_USE_LEPSILON
  description: Use LEPSILON for the calculation of the Born effective charge tensor. Otherwise, use LCALCEPS.
  type: bool
  default: 1

DEFAULT_APL_FREQFORMAT
  description: Format for the phonon frequencies.
  type: string
  default: "THZ | ALLOW_NEGATIVE"

DEFAULT_APL_DC
  description: Calculate phonon dispersion curves.
  type: bool
  default: 1

DEFAULT_APL_DCPATH
  description: The q-point path of the phonon dispersion curve.
  type: string
  default: "LATTICE"

DEFAULT_APL_DCPOINTS
  description: Number of points between two high-symmetry points in the phonon dispersion curve.
  type: int
  default: 100

DEFAULT_APL_DOS
  description: Calculate the phonon density of states.
  type: bool
  default: 1

DEFAULT_APL_DOSMETHOD
  description: Determines whether to use the linear tetrahedron (LT) or root sampling (RS) method for the phonon DOS.
  type: string
  choices: "LT", "RS"
  default:

DEFAULT_APL_DOSPOINTS
  description: Number of points for the phonon DOS.
  type: int
  default: 100

DEFAULT_APL_DOSSMEAR
  description: Smearing value for the root sampling method for the phonon DOS.
  type: double
  default: 0.0

DEFAULT_APL_TP
  description: Calculate thermodynamic properties.
  type: bool
  default: 1

DEFAULT_APL_TPT
  description: Temperature range (Tstart, Tend) and step size (Tstep) for the thermodynamic properties calculation.
  type: string
  default: "0:2000:10"

DEFAULT_APL_FILE_PREFIX
  description: APL file prefix. All APL output files start with this prefix.
  type: string
  default: "aflow.apl."

DEFAULT_APL_FILE_PREFIX
  description: The APL summary file.
  type: string
  default: ".out"

DEFAULT_APL_PDIS_FILE
  description: Output file for the phonon dispersion curve.
  type: string
  default: "phonon_dispersion.out"

DEFAULT_APL_PDOS_FILE
  description: Output file for the phonon DOS.
  type: string
  default: "phonon_dos.out"

DEFAULT_APL_THERMO_FILE
  description: Output file for the thermodynamic properties.
  type: string
  default: "thermodynamic_properties.out"

DEFAULT_APL_THERMO_JSON
  description: Output file for the thermodynamic properties in json format.
  type: string
  default: "thermodynamic_properties.json"

DEFAULT_APL_HARMIFC_FILE
  description: Output file for the harmonic force constants. Also functions as the APL hibernation file.
  type: string
  default: "harmonicIFCs.xml"

DEFAULT_APL_POLAR_FILE
  description: Output file for the Born effective charge tensors and the dielectric tensor.
  type: string
  default: "polar.xml"

DEFAULT_APL_HSKPTS_FILE
  description: Output file for the high-symmetry k-points in APL.
  type: string
  default: "hskpoints.out"

DEFAULT_APL_MSQRDISP_FILE
  description: Output file for the mean square displacements.
  type: string
  default: "displacements.out"

DEFAULT_APL_GVEL_FILE
  description: Output file for the group velocities.
  type: string
  default: "group_velocities.out"

DEFAULT_APL_PHDOSCAR_FILE
  description: Output file for the phonon DOS in VASP DOSCAR format.
  type: string
  default: "PHDOSCAR"

DEFAULT_APL_PHPOSCAR_FILE
  description: Output file for the unit cell of the phonon calcultions in VASP5 POSCAR format.
  type: string
  default: "PHPOSCAR"

DEFAULT_APL_PHKPOINTS_FILE
  description: Output file for the phonon dispersion path in VASP KPOINTS format.
  type: string
  default: "PHKPOINTS"

DEFAULT_APL_PHEIGENVAL_FILE
  description: Output file for the phonon dispersion in VASP EIGENVAL format.
  type: string
  default: "PHEIGENVAL"

DEFAULT_APL_STATE_FILE
  description: Output file to save the state of the force constant calculators.
  type: string
  default: "fccalc_state.out"

DEFAULT_APL_ADISP_SCENE_FORMAT
  description: Default format to visualize atomic displacements for phonon modes.
  type: string
  default: "XCRYSDEN"

DEFAULT_APL_ADISP_AMPLITUDE
  description: Default amplitude for atomic displacements.
  type: double
  default: 0.25

DEFAULT_APL_ADISP_NSTEPS
  description: The number of steps per period in an XCRYSDEN animated structure file.
  type: integer
  default: 20

DEFAULT_APL_ADISP_NPERIODS
  description: The number of periods in an XCRYSDEN animated structure file.
  type: integer
  default: 1

####################################################################################################
DEFAULTS QHA

DEFAULT_QHA_MODE
  description: Defines QHA mode to run (it is possible to run several modes simultaneously).
  type: string
  choices: "QHA", "QHA3P", "SCHQA" or their comma-separated combination
  default: "QHA"

DEFAULT_QHA_EOS
  description: Runs equation of state (EOS) calculation (includes several APL runs) to obtain thermodynamic properties of material.
  type: bool
  default: 1

DEFAULT_QHA_EOS_DISTORTION_RANGE
  description: Sets the range of volumes and the volume increment for the set of EOS calculations.
  type: string
  default: "-12:16:3"

DEFAULT_QHA_EOS_MODEL
  description: Defines what equation of state model(s) to use for the E(V) relation fitting.
  type: string
  choices: "SJ", "BM3", "BM4", "M" or their comma-separated combination
  default: "SJ,BM4"

DEFAULT_QHA_GP_DISTORTION
  description: Sets the value of volume change (dV) for the calculation of the Grueneisen parameter using the finite difference method.
  type: double
  default: 1.0

DEFAULT_QHA_TAYLOR_EXPANSION_ORDER
  description: Sets the order of Taylor expansion that is used by QHANP method.
  type: int
  default: 2

DEFAULT_QHA_INCLUDE_ELEC_CONTRIB
  description: Defines if electronic contribution to the free energy is included.
  type: bool
  default: 0

DEFAULT_QHA_SOMMERFELD_EXPANSION
  description: Defines if electronic contribution to the free energy is calculated using the Sommerfeld expansion.
  type: bool
  default: 0

DEFAULT_QHA_PDIS_T
  description: Defines the set of temperatures to calculate the temperature-dependent phonon band structure with the SCQHA method.
  type: string
  default: "50,100,600"

DEFAULT_QHA_GP_FINITE_DIFF
  description: Defines if Grueneisen parameter is calculated using finite difference method.
  type: bool
  default: 0

DEFAULT_QHA_IGNORE_IMAGINARY
  description: Defines if imaginary frequencies in phonon calculations are ignored.
  type: bool
  default: 0

DEFAULT_QHA_RELAX_IONS_CELL
  description: Defines if ions and cell shape will be relaxed in QHA-related subcalculations.
  type: bool
  default: 0

DEFAULT_QHA_FILE_PREFIX
  description: QHA file prefix. All QHA output files start with this prefix.
  type: string
  default: "aflow.qha."

DEFAULT_QHA3P_FILE_PREFIX
  description: QHA3P file prefix. All QHA3P output files start with this prefix.
  type: string
  default: "aflow.qha3p."

DEFAULT_QHANP_FILE_PREFIX
  description: QHANP file prefix. All QHANP output files start with this prefix.
  type: string
  default: "aflow.qhanp."

DEFAULT_SCQHA_FILE_PREFIX
  description: SCQHA file prefix. All SCQHA output files start with this prefix.
  type: string
  default: "aflow.scqha."

DEFAULT_QHA_GP_PATH_FILE
  description: Output file for mode-resolved Grueneisen parameter dispersion curve.
  type: string
  default: "gp.disp.out"

DEFAULT_QHA_GP_MESH_FILE
  description: Output file for mode-resolved Grueneisen parameter calculated on q-points mesh.
  type: string
  default: "gp.mesh.out"

DEFAULT_QHA_GP_AVG_FILE
  description: Output file for temperature-dependent average Grueneisen parameter calculated using finite differences.
  type: string
  default: "gp.avg.out"

DEFAULT_QHA_THERMO_FILE
  description: Output file for thermal properties (equilibrium volume, free energy, bulk modulus, isobaric and isochoric specific heat, average Grueneisen parameter).
  type: string
  default: "gp.thermo.out"

DEFAULT_QHA_FREQS_FILE
  description: Output file for volume-dependent phonon frequencies.
  type: string
  default: "frequencies.out"

DEFAULT_QHA_FVT_FILE
  description: Output file for temperature- and volume-dependent free energy.
  type: string
  default: "FVT.out"

DEFAULT_QHA_PDIS_FILE
  description: Output file for temperature-dependent phonon dispersion relations.
  type: string
  default: "dispersion_phonon"

DEFAULT_QHA_PDOS_FILE
  description: Output file for temperature-dependent phonon DOS.
  type: string
  default: "dos_phonon"

DEFAULT_QHA_KPOINTS_FILE
  description: Output file for path in the BZ used to calculated temperature-dependent phonon dispersion relations.
  type: string
  default: "kpoints.out"

####################################################################################################
DEFAULTS AAPL

DEFAULT_AAPL_BTE
  description: Describes how to solve the Boltzmann Transport Equation.
  type: string
  choices: "FULL", "RTA"
  default: "FULL"

DEFAULT_AAPL_CUT_RAD
  description: Cut-off radius for anharmonic force constants.
  type: string
  format: "rad3"
  default: "0.0"

DEFAULT_AAPL_CUT_SHELL
  description: Cut-off shell for anharmonic force constants.
  type: string
  format: "shell3"
  default: "6"

DEFAULT_AAPL_THERMALGRID
  description: q-point grid for thermal conductivity calculations.
  type: string
  format: "N1xN2xN3"
  default: "21x21x21"

DEFAULT_AAPL_TCT
  description: Temperature range (Tstart, Tend) and step size (Tstep) for the thermal conductivity calculations.
  type: string
  format: "Tstart:Tend:Tstep"
  default: "50:550:50"

DEFAULT_AAPL_SUMRULE
  description: Convergence criterion for the sumrule in the determination of anharmonic force constants.
  type: double
  default: 1e-7

DEFAULT_AAPL_MIXING_COEFFICIENT
  description: The mixing coefficient in the iteration procedure for anharmonic force constants.
  type: double
  default: 0.0

DEFAULT_AAPL_ISOTOPE
  description: Includes isotope scattering processes for thermal conductivity calculations.
  type: bool
  default: 1

DEFAULT_AAPL_BOUNDARY
  description: Includes boundary corrections for thermal conductivity calculations.
  type: bool
  default: 0

DEFAULT_AAPL_CUMULATIVEK
  description: Uses the free mean path for boundary corrections.
  type: bool
  default: 0

DEFAULT_AAPL_NANO_SIZE
  description: Grain boundary/mean free path size for boundary corrections to the thermal conductivity.
  type: double
  default: 100.0

DEFAULT_AAPL_FILE_PREFIX
  description: AAPL file prefix. All AAPL output files start with this prefix.
  type: string
  default: "aflow.aapl."

DEFAULT_AAPL_IRRQPTS_FILE
  description: Output file for the irreducible q-points used for thermal conductivity calculations.
  type: string
  default: "irred_qpoints.out"

DEFAULT_AAPL_GVEL_FILE
  description: Output file for the group velocities of the phonon modes and frequencies used for thermal conductivity calculations.
  type: string
  default: "group_velocities.out"

DEFAULT_AAPL_PS_FILE
  description: Output file for the phase space of the phonon modes obtained from thermal conductivity calculations.
  type: string
  default: "phase_space.out"

DEFAULT_AAPL_GRUENEISEN_FILE
  description: Output file for the average and mode Grueneisen parameters obtained from thermal conductivity calculations.
  type: string
  default: "grueneisen.out"

DEFAULT_AAPL_RATES_FILE
  description: Output file for the total phonon scattering rates.
  type: string
  default: "scattering_rates_total.out"

DEFAULT_AAPL_RATES_3RD_FILE
  description: Output file for the scattering rates from three-phonon scattering processes.
  type: string
  default: "scattering_rates_anharmonic_3rd.out"

DEFAULT_AAPL_RATES_4TH_FILE
  description: Output file for the scattering rates from four-phonon scattering processes.
  type: string
  default: "scattering_rates_anharmonic_4th.out"

DEFAULT_AAPL_ISOTOPE_FILE
  description: Output file for the scattering rates from isotope scattering processes.
  type: string
  default: "scattering_rates_isotope.out"

DEFAULT_AAPL_BOUNDARY_FILE
  description: Output file for the scattering rates from grain boundary scattering processes.
  type: string
  default: "scattering_rates_boundary.out"

DEFAULT_AAPL_TCOND_FILE
  description: Output file for the thermal conductivity values.
  type: string
  default: "thermal_conductivity.out"

DEFAULT_AAPL_TCOND_PLOT_FILE
  description: Gnuplot file for the thermal conductivity plot.
  type: string
  default: "thermal_conductivity.plt"

####################################################################################################
DEFAULTS XTALFINDER

DEFAULT_XTALFINDER_MISFIT_MATCH
  description: Similarity threshold for similar structures.
  type: double
  default: 0.1

DEFAULT_XTALFINDER_MISFIT_FAMILY
  description: Similarity threshold for structures in the same family.
  type: double
  default: 0.2

DEFAULT_XTALFINDER_MISFIT_SUPERCELL_METHOD
  description: Supercell method for converting crystal structures into a new representation (well tested, but slow). The alternative (default) is the transformation method (fast).
  type: bool
  default: 0

DEFAULT_XTALFINDER_SAFE_ATOM_MATCH_SCALING
  description: Scaling factor that divides the minimum interatomic distance, establishing when mapped atoms are sufficiently close and other mappings do not need to be explored.
  type: double
  default: 4.0

DEFAULT_XTALFINDER_FILE_MATERIAL
  description: Output file prefix for comparing materials (map alike atomic species).
  type: string
  default: "material_comparison_output"

DEFAULT_XTALFINDER_FILE_STRUCTURE
  description: Output file prefix for comparing structures (ignore atomic species).
  type: string
  default: "structure_comparison_output"

DEFAULT_XTALFINDER_FILE_DUPLICATE
  description: Output file prefix for finding duplicate compounds (map alike atomic species and do not scale volumes).
  type: string
  default: "duplicate_comparison_output"

DEFAULT_XTALFINDER_FILE_MATERIAL_COMPARE2DATABASE
  description: Output file prefix for comparing an input material to the AFLOW.org repository (map alike atomic species).
  type: string
  default: "material_comparison_compare2database_output"

DEFAULT_XTALFINDER_FILE_STRUCTURE_COMPARE2DATABASE
  description: Output file prefix for comparing an input structure to the AFLOW.org repository (ignore atomic species).
  type: string
  default: "structure_comparison_compare2database_output"

DEFAULT_XTALFINDER_FILE_MATERIAL_DATABASE
  description: Output file prefix for comparing materials in the AFLOW.org repository (map alike atomic species).
  type: string
  default: "material_comparison_database_output"

DEFAULT_XTALFINDER_FILE_STRUCTURE_DATABASE
  description: Output file prefix for comparing structures in the AFLOW.org repository (ignore atomic species).
  type: string
  default: "structure_comparison_database_output"

####################################################################################################
DEFAULTS ANRL

DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL
  description: Fractional tolerance value for identifying equivalent Wyckoff coordinates.
  type: double
  default: 1e-6

####################################################################################################
DEFAULTS CORE

AFLOW_CORE_TEMPERATURE_BEEP
  description: Threshold computer temperature (Celsius) above which the machine beeps.
  type: integer
  default: 56

AFLOW_CORE_TEMPERATURE_HALT
  description: Threshold computer temperature (Celsius) above which AFLOW stops running.
  type: integer
  default: 65

AFLOW_CORE_TEMPERATURE_REFRESH
  description: Refresh rate (seconds) for computer temperature readings.
  type: integer
  default: 5

####################################################################################################
DEFAULTS VASP MACHINE

SECONDS_SLEEP_VASP_COMPLETION
  description: Time (seconds) to wait for VASP to finish writing output files. Useful on clusters with NFS delays. Twice this number should be less than SECONDS_SLEEP_VASP_MONITOR to run with --monitor_vasp.
  type: integer
  default: 15

SECONDS_SLEEP_VASP_MONITOR
  description: Time (seconds) to wait between monitoring loops for --monitor_vasp.
  type: integer
  default: 60

SECONDS_STALE_VASP_OUT
  description: Time (seconds) at which the DEFAULT_VASP_OUT file is considered stale (since last modified).
  type: integer
  default: 21600 (6 hours)

BYTES_MAX_VASP_OUT
  description: Maximum size of the DEFAULT_VASP_OUT file. The larger the size, the longer the delays in applying patches. Keep as large as possible, it is better to let vasp finish gracefully.
  type: integer
  default: 20000000000 (20GB)

MEMORY_MAX_USAGE_RAM
  description: The ram memory limit (in percent) before a MEMORY error is issued by --monitor_vasp.
  type: integer
  default: 98

MEMORY_MAX_USAGE_SWAP
  description: The swap memory limit (in percent) before a MEMORY error is issued by --monitor_vasp.
  type: integer
  default: 45

FILE_VASP_MONITOR
  description: File name for the aflow vasp monitor (--monitor_vasp). In the default case, the name would be LOCK.monitor_vasp.
  type: string
  default: "monitor_vasp"

INTEL_COMPILER_PATHS
  description: Comma-separated paths to search for sourcing the Intel compiler environmental variables.
  type: string
  default: "/opt/intel/bin/compilervars.sh,/opt/intel/bin/compilervars.csh,/app/intel/parallel_studio_xe_2020_update1/bin/compilervars.sh"

####################################################################################################
DEFAULTS MACHINE DEPENDENT MPI

MPI_OPTIONS_DUKE_BETA_MPICH
  description: MPI options for --machine=beta.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_BETA_MPICH
  description: MPI command for --machine=beta.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_BETA_MPICH
  description: MPI executable directory for --machine=beta.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_BETA_OPENMPI
  description: MPI options for --machine=beta_openmpi.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_BETA_OPENMPI
  description: MPI command for --machine=beta_openmpi.
  type: string
  default: "/usr/bin/mpirun.openmpi -np"

MPI_BINARY_DIR_DUKE_BETA_OPENMPI
  description: MPI executable directory for --machine=beta_openmpi.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_MATERIALS
  description: MPI options for --machine=materials.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_MATERIALS
  description: MPI command for --machine=materials.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_MATERIALS
  description: MPI executable directory for --machine=materials.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_AFLOWLIB
  description: MPI options for --machine=aflowlib.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_AFLOWLIB
  description: MPI command for --machine=aflowlib.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_AFLOWLIB
  description: MPI executable directory for --machine=aflowlib.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_HABANA
  description: MPI options for --machine=habana.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_HABANA
  description: MPI command for --machine=habana.
  type: string
  default: "/usr/bin/mpiexec -np"

MPI_BINARY_DIR_DUKE_HABANA
  description: MPI executable directory for --machine=habana.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_QRATS_MPICH
  description: MPI options for --machine=qrats.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_QRATS_MPICH
  description: MPI command for --machine=qrats.
  type: string
  default: "/MAIN/bin/MPICH/bin/mpirun -np"

MPI_BINARY_DIR_DUKE_QRATS_MPICH
  description: MPI executable directory for --machine=qrats.
  type: string
  default: "/usr/local/bin/"

MPI_OPTIONS_DUKE_QFLOW_OPENMPI
  description: MPI options for --machine=qflow.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_QFLOW_OPENMPI
  description: MPI command for --machine=qflow.
  type: string
  default: "/usr/bin/mpirun -n"

MPI_BINARY_DIR_DUKE_QFLOW_OPENMPI
  description: MPI executable directory for --machine=qflow.
  type: string
  default: "/home/bin/"

MPI_OPTIONS_DUKE_X
  description: MPI options for --machine=x.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_DUKE_X
  description: MPI command for --machine=x.
  type: string
  default: "mpirun -n"

MPI_BINARY_DIR_DUKE_X
  description: MPI executable directory for --machine=x.
  type: string
  default: "/home/bin/"

MPI_OPTIONS_MACHINE001
  description: MPI options for --machine=machine001.
  type: string
  default: ""

MPI_COMMAND_MACHINE001
  description: MPI command for --machine=machine001.
  type: string
  default: "aprun -n"

MPI_BINARY_DIR_MACHINE001
  description: MPI executable directory for --machine=machine001.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MACHINE002
  description: MPI options for --machine=machine002.
  type: string
  default: ""

MPI_COMMAND_MACHINE002
  description: MPI command for --machine=machine002.
  type: string
  default: "/p/app/intel/parallel_studio_xe_2017_update4/impi/2017.3.196/intel64/bin/mpirun -np"

MPI_BINARY_DIR_MACHINE002
  description: MPI executable directory for --machine=machine002.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_EOS
  description: MPI options for --machine=eos.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_EOS
  description: MPI command for --machine=eos.
  type: string
  default: "/usr/bin/srun -n"

MPI_NCPUS_MPCDF_EOS
  description: MPI core limit for --machine=eos.
  type: string
  default: "32"

MPI_HYPERTHREADING_MPCDF_EOS
  description: Adjusts NCPUS for hyperthreading with --machine=eos. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "NEGLECT"

MPI_BINARY_DIR_MPCDF_EOS
  description: MPI executable directory for --machine=eos.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_DRACO
  description: MPI options for --machine=draco.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_DRACO
  description: MPI command for --machine=draco.
  type: string
  default: "/usr/bin/srun -n"

MPI_NCPUS_MPCDF_DRACO
  description: MPI core limit for --machine=draco.
  type: string
  default: "0"

MPI_HYPERTHREADING_MPCDF_DRACO
  description: Adjusts NCPUS for hyperthreading with --machine=draco. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "OFF"

MPI_BINARY_DIR_MPCDF_DRACO
  description: MPI executable directory for --machine=draco.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_COBRA
  description: MPI options for --machine=cobra.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_COBRA
  description: MPI command for --machine=cobra.
  type: string
  default: "/usr/bin/srun -n"

MPI_NCPUS_MPCDF_COBRA
  description: MPI core limit for --machine=cobra.
  type: string
  default: "0"

MPI_HYPERTHREADING_MPCDF_COBRA
  description: Adjusts NCPUS for hyperthreading with --machine=cobra. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "OFF"

MPI_BINARY_DIR_MPCDF_COBRA
  description: MPI executable directory for --machine=cobra.
  type: string
  default: "~/bin/"

MPI_OPTIONS_MPCDF_HYDRA
  description: MPI options for --machine=hydra.
  type: string
  default: "ulimit -s unlimited "

MPI_COMMAND_MPCDF_HYDRA
  description: MPI command for --machine=hydra.
  type: string
  default: "poe "

MPI_NCPUS_MPCDF_HYDRA
  description: MPI core limit for --machine=hydra.
  type: string
  default: "0"

MPI_HYPERTHREADING_MPCDF_HYDRA
  description: Adjusts NCPUS for hyperthreading with --machine=hydra. "NEGLECT" does nothing. "OFF" divides NCPUS in half. "ON" doubles NCPUS.
  type: string
  choices: "NEGLECT", "OFF", "ON"
  default: "OFF"

MPI_BINARY_DIR_MPCDF_HYDRA
  description: MPI executable directory for --machine=hydra.
  type: string
  default: "~/bin/"

MPI_OPTIONS_FULTON_MARYLOU
  description: MPI options for --machine=marylou.
  type: string
  default: "export OMP_NUM_THREADS

MPI_COMMAND_FULTON_MARYLOU
  description: MPI command for --machine=marylou.
  type: string
  default: "mpiexec"

MPI_BINARY_DIR_FULTON_MARYLOU
  description: MPI executable directory for --machine=marylou.
  type: string
  default: "/fslgroup/fslg_datamining/bin/"

MPI_OPTIONS_MACHINE1
  description: MPI options for --machine=host1.
  type: string
  default: ""

MPI_COMMAND_MACHINE1
  description: MPI command for --machine=host1.
  type: string
  default: "...something ..."

MPI_BINARY_DIR_MACHINE1
  description: MPI executable directory for --machine=host1.
  type: string
  default: "/somewhere/"

MPI_OPTIONS_MACHINE2
  description: MPI options for --machine=ohad.
  type: string
  default: ""

MPI_COMMAND_MACHINE2
  description: MPI command for --machine=ohad.
  type: string
  default: "stub not used"

MPI_BINARY_DIR_MACHINE2
  description: MPI executable directory for --machine=ohad.
  type: string
  default: "/home/aflow/bin/"


****************************************************************************************************

For additional information contact: Stefano Curtarolo (stefano@duke.edu) or 
                                    Corey Oses (corey.oses@duke.edu)

****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
