#ifndef AFLOW_XCLASSES_H
#define AFLOW_XCLASSES_H

#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xserialization.h"

#include "structure/aflow_xstructure.h"

class xstructure;

// ME20181026 - Container for APL options
struct _moduleOptions {
  // APL
  std::vector<aurostd::xoption> aplflags;

  // AAPL
  std::vector<aurostd::xoption> aaplflags;

  // QHA
  std::vector<aurostd::xoption> qhaflags;

  // AEL
  std::vector<aurostd::xoption> aelflags;

  // AGL
  std::vector<aurostd::xoption> aglflags;
};

// CO20180420
// for stream management with objects
class xStream {
public:
  // NECESSARY PUBLIC CLASS METHODS - START
  // constructors - START
  xStream(std::ostream& oss = std::cout);
  xStream(std::ofstream& ofs, std::ostream& oss = std::cout);
  // constructors - STOP
  ~xStream();
  // NECESSARY PUBLIC CLASS METHODS - END

  // initializers
  void initialize(std::ostream& oss = std::cout); // ME20200427
  void initialize(std::ofstream& ofs, std::ostream& oss = std::cout); // ME20200427

  // getters
  [[nodiscard]] std::ostream* getOSS() const; // CO20191110
  [[nodiscard]] std::ofstream* getOFStream() const; // CO20191110
protected:
  // NECESSARY private CLASS METHODS - START
  void free();
  void copy(const xStream& b);
  // NECESSARY END CLASS METHODS - END
  // logger variables
  std::ostream* p_oss;
  std::ofstream* p_FileMESSAGE;
  bool f_new_ofstream; // for deletion later
  // general setters
  void setOFStream(std::ofstream& FileMESSAGE);
  void setOSS(std::ostream& oss);
};

// --------------------------------------------------------------------------
// general flags to run aflow
class _aflags : public JsonSerializable<_aflags> {
public:
  // trivial constructurs/destuctors/operators
  _aflags(); // default, just allocate
  ~_aflags(); // kill everything
  _aflags(const _aflags& b); // constructor copy
  const _aflags& operator=(const _aflags& b); // copy
  void clear(); // clear
  // CONTENT
  bool QUIET;
  int AFLOW_PTHREADS_NUMBER; // cant be GLOBAL as this is a local run stuff
  // particular
  std::string LocalDirectory; // where is aflow now
  std::string Directory; // where aflow must run
  bool AFLOW_FORCE_RUN; // Force run also in database
  bool AFLOW_PERFORM_DIRECTORY; // Directory is specified (sometimes it is useful).
  bool AFLOW_PERFORM_FILE; // File is specified (sometimes it is useful).
  bool AFLOW_PERFORM_ORDER_SORT; // Sorts the _AFLOWIN_ in the list
  bool AFLOW_PERFORM_ORDER_REVERSE; // Reverse the _AFLOWIN_ in the list
  bool AFLOW_PERFORM_ORDER_RANDOM; // Randomize the _AFLOWIN_ in the list
  bool AFLOW_MODE_GENERATE; // TODO OVVERRIDE all _AFLOWIN_
  bool AFLOW_MODE_QSUB_MODE1; // TODO OVVERRIDE all _AFLOWIN_
  bool AFLOW_MODE_QSUB_MODE2; // TODO OVVERRIDE all _AFLOWIN_
  bool AFLOW_MODE_QSUB_MODE3; // TODO OVVERRIDE all _AFLOWIN_
  // general flags to operate in the directory
  bool KBIN_RUN_AFLOWIN;
  bool KBIN_GEN_GENERAL; // CO20180409
  bool KBIN_GEN_VASP_FROM_AFLOWIN;
  bool KBIN_GEN_AIMS_FROM_AFLOWIN;
  bool KBIN_GEN_AFLOWIN_FROM_VASP;
  // DX START
  bool KBIN_GEN_SYMMETRY_OF_AFLOWIN;
  // DX END
  bool KBIN_DELETE_AFLOWIN;
  bool AFLOW_FORCE_MPI; // not yet implemented
  bool AFLOW_FORCE_SERIAL; // not yet implemented
  int AFLOW_GLOBAL_NCPUS; // Forced CPUS
  // Perform TASKS
  bool AFLOW_PERFORM_CLEAN; // to clean a directory
  // host related things
  aurostd::xoption AFLOW_MACHINE_GLOBAL;
  aurostd::xoption AFLOW_MACHINE_LOCAL; // flag for duke_beta_mpich
  aurostd::xoption vflag;

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  _aflags deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "_aflags"; }

private: //
  void free(); // free space
  void copy(const _aflags& b); //
};

// SERIALIZATION MEMBERS, ignoring xoption members
#define JSON_aflags_MEMBERS                                                                                                                                                                                 \
  QUIET, AFLOW_PTHREADS_NUMBER, LocalDirectory, Directory, AFLOW_FORCE_RUN, AFLOW_PERFORM_DIRECTORY, AFLOW_PERFORM_FILE, AFLOW_PERFORM_ORDER_SORT, AFLOW_PERFORM_ORDER_REVERSE, AFLOW_PERFORM_ORDER_RANDOM, \
      AFLOW_MODE_GENERATE, AFLOW_MODE_QSUB_MODE1, AFLOW_MODE_QSUB_MODE2, AFLOW_MODE_QSUB_MODE3, KBIN_RUN_AFLOWIN, KBIN_GEN_GENERAL, KBIN_GEN_VASP_FROM_AFLOWIN, KBIN_GEN_AIMS_FROM_AFLOWIN,                 \
      KBIN_GEN_AFLOWIN_FROM_VASP, KBIN_GEN_SYMMETRY_OF_AFLOWIN, KBIN_DELETE_AFLOWIN, AFLOW_FORCE_MPI, AFLOW_FORCE_SERIAL, AFLOW_GLOBAL_NCPUS, AFLOW_PERFORM_CLEAN

// --------------------------------------------------------------------------
// general flags for kbinary (all)
class _kflags : public JsonSerializable<_kflags> {
public:
  // trivial constructurs/destuctors/operators
  _kflags(); // default, just allocate
  ~_kflags(); // kill everything
  _kflags(const _kflags& b); // constructor copy
  const _kflags& operator=(const _kflags& b); // copy
  void clear(); // clear
  // CONTENT
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in KBIN MODE
  //
  bool AFLOW_MODE_ALIEN;
  //
  bool AFLOW_MODE_MATLAB;
  bool AFLOW_MATLAB_MODE_EXPLICIT;
  bool AFLOW_MATLAB_MODE_EXPLICIT_START_STOP;
  bool AFLOW_MATLAB_MODE_IMPLICIT;
  bool AFLOW_MATLAB_MODE_EXTERNAL;
  bool AFLOW_MATLAB_FILE;
  bool AFLOW_MATLAB_FILE_FILE;
  bool AFLOW_MATLAB_FILE_COMMAND;
  //
  bool AFLOW_MODE_VASP;
  bool AFLOW_MODE_AIMS;
  //
  bool AFLOW_MODE_PRESCRIPT_EXPLICIT;
  bool AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP;
  std::stringstream AFLOW_MODE_PRESCRIPT;
  bool AFLOW_MODE_POSTSCRIPT_EXPLICIT;
  bool AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP;
  std::stringstream AFLOW_MODE_POSTSCRIPT;
  //
  bool AFLOW_MODE_EMAIL;
  // normal binary
  std::string KBIN_BIN;
  std::string KBIN_SERIAL_BIN; // ME20190107
  std::string KZIP_BIN;
  bool KZIP_COMPRESS;
  // MPI binaries and flags
  bool KBIN_MPI;
  int KBIN_MPI_NCPUS;
  std::string KBIN_MPI_NCPUS_STRING; // ME20181216
  int KBIN_MPI_NCPUS_ORIG; // CO20210804 - repurposing
  std::string KBIN_MPI_START;
  std::string KBIN_MPI_STOP;
  std::string KBIN_MPI_COMMAND;
  bool KBIN_MPI_AUTOTUNE;
  std::string KBIN_MPI_BIN;
  std::string KBIN_MPI_OPTIONS;
  // QSUB
  bool KBIN_QSUB;
  bool KBIN_QSUB_MODE1;
  bool KBIN_QSUB_MODE2;
  bool KBIN_QSUB_MODE3;
  std::string KBIN_QSUB_COMMAND;
  std::string KBIN_QSUB_PARAMS;
  bool KBIN_QSUB_MODE_EXPLICIT;
  bool KBIN_QSUB_MODE_EXPLICIT_START_STOP;
  bool KBIN_QSUB_MODE_IMPLICIT;
  bool KBIN_QSUB_FILE;
  // symmetry operation lists
  bool KBIN_SYMMETRY_CALCULATION;
  // DX START
  bool KBIN_SYMMETRY_NO_SCAN;
  double KBIN_SYMMETRY_EPS;
  bool KBIN_SYMMETRY_CALCULATE_PGROUP; // DX20170814 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_PGROUPK; // DX20170814 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_FGROUP; // DX20170814 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL; // DX20170814 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL; // DX20171205 - Specify what to calculate/verify; Added pgroupk_xtal
  bool KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON; // DX20200129 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_IATOMS; // DX20170814 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_AGROUP; // DX20170814 - Specify what to calculate/verify
  bool KBIN_SYMMETRY_CALCULATE_SGROUP; // DX20170814 - Specify what to calculate/verify
  // DX END
  bool KBIN_SYMMETRY_PGROUP_WRITE; // taken true by default
  bool KBIN_SYMMETRY_PGROUPK_WRITE; // taken true by default
  bool KBIN_SYMMETRY_PGROUP_XTAL_WRITE; // taken true by default
  bool KBIN_SYMMETRY_PGROUPK_XTAL_WRITE; // DX20171205 - Added pgroupk_xtal
  bool KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE; // DX20200129 - taken true by default
  bool KBIN_SYMMETRY_FGROUP_WRITE; // taken true by default
  bool KBIN_SYMMETRY_SGROUP_WRITE;
  bool KBIN_SYMMETRY_AGROUP_WRITE; // taken true by default
  bool KBIN_SYMMETRY_IATOMS_WRITE; // taken true by default
  double KBIN_SYMMETRY_SGROUP_RADIUS;
  // pocc operation lists
  bool KBIN_POCC;
  bool KBIN_POCC_CALCULATION;
  std::string KBIN_POCC_TEMPERATURE_STRING; // CO20191114
  std::string KBIN_POCC_ARUNS2SKIP_STRING; // CO20200627
  bool KBIN_POCC_EXCLUDE_UNSTABLE; // ME20210927
  // frozsl operation lists
  bool KBIN_FROZSL;
  bool KBIN_FROZSL_DOWNLOAD;
  bool KBIN_FROZSL_FILE;
  std::string KBIN_FROZSL_FILE_NAME;
  bool KBIN_FROZSL_STRUCTURE_MODE_FILE;
  bool KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP;
  std::string KBIN_FROZSL_STRUCTURE_STRING;
  bool KBIN_FROZSL_DIELECTRIC_MODE_FILE;
  bool KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP;
  bool KBIN_FROZSL_DIELECTRIC_ZEFF;
  std::string KBIN_FROZSL_DIELECTRIC_STRING;
  // phonons operation lists
  bool KBIN_PHONONS_CALCULATION_APL;
  bool KBIN_PHONONS_CALCULATION_QHA; // CO20170601
  bool KBIN_PHONONS_CALCULATION_AAPL; // CO20170601
  bool KBIN_PHONONS_CALCULATION_AGL;
  bool KBIN_PHONONS_CALCULATION_AEL;
  bool KBIN_PHONONS_CALCULATION_FROZSL;
  std::string KBIN_PHONONS_CALCULATION_FROZSL_output;
  std::string KBIN_PHONONS_CALCULATION_FROZSL_poscars;
  _moduleOptions KBIN_MODULE_OPTIONS; // ME20181027

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  _kflags deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "_kflags"; }

private: //
  void free(); // free space
  void copy(const _kflags& b); //

    // SERIALIZATION MEMBERS Avoiding following: AFLOW_MODE_PRESCRIPT, AFLOW_MODE_POSTSCRIPT, KBIN_MODULE_OPTIONS
#define JSON_kflags_MEMBERS                                                                                                                                                                                      \
  AFLOW_MODE_ALIEN, AFLOW_MODE_MATLAB, AFLOW_MATLAB_MODE_EXPLICIT, AFLOW_MATLAB_MODE_EXPLICIT_START_STOP, AFLOW_MATLAB_MODE_IMPLICIT, AFLOW_MATLAB_MODE_EXTERNAL, AFLOW_MATLAB_FILE, AFLOW_MATLAB_FILE_FILE,     \
      AFLOW_MATLAB_FILE_COMMAND, AFLOW_MODE_VASP, AFLOW_MODE_AIMS, AFLOW_MODE_PRESCRIPT_EXPLICIT, AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP, AFLOW_MODE_POSTSCRIPT_EXPLICIT,                                      \
      AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP, AFLOW_MODE_EMAIL, KBIN_BIN, KBIN_SERIAL_BIN, KZIP_BIN, KZIP_COMPRESS, KBIN_MPI, KBIN_MPI_NCPUS, KBIN_MPI_NCPUS_STRING, KBIN_MPI_NCPUS_ORIG, KBIN_MPI_START,     \
      KBIN_MPI_STOP, KBIN_MPI_COMMAND, KBIN_MPI_AUTOTUNE, KBIN_MPI_BIN, KBIN_MPI_OPTIONS, KBIN_QSUB, KBIN_QSUB_MODE1, KBIN_QSUB_MODE2, KBIN_QSUB_MODE3, KBIN_QSUB_COMMAND, KBIN_QSUB_PARAMS,                     \
      KBIN_QSUB_MODE_EXPLICIT, KBIN_QSUB_MODE_EXPLICIT_START_STOP, KBIN_QSUB_MODE_IMPLICIT, KBIN_QSUB_FILE, KBIN_SYMMETRY_CALCULATION, KBIN_SYMMETRY_NO_SCAN, KBIN_SYMMETRY_EPS, KBIN_SYMMETRY_CALCULATE_PGROUP, \
      KBIN_SYMMETRY_CALCULATE_PGROUPK, KBIN_SYMMETRY_CALCULATE_FGROUP, KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL, KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL, KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON,                     \
      KBIN_SYMMETRY_CALCULATE_IATOMS, KBIN_SYMMETRY_CALCULATE_AGROUP, KBIN_SYMMETRY_CALCULATE_SGROUP, KBIN_SYMMETRY_PGROUP_WRITE, KBIN_SYMMETRY_PGROUPK_WRITE, KBIN_SYMMETRY_PGROUP_XTAL_WRITE,                  \
      KBIN_SYMMETRY_PGROUPK_XTAL_WRITE, KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE, KBIN_SYMMETRY_FGROUP_WRITE, KBIN_SYMMETRY_SGROUP_WRITE, KBIN_SYMMETRY_AGROUP_WRITE, KBIN_SYMMETRY_IATOMS_WRITE,                   \
      KBIN_SYMMETRY_SGROUP_RADIUS, KBIN_POCC, KBIN_POCC_CALCULATION, KBIN_POCC_TEMPERATURE_STRING, KBIN_POCC_ARUNS2SKIP_STRING, KBIN_POCC_EXCLUDE_UNSTABLE, KBIN_FROZSL, KBIN_FROZSL_DOWNLOAD, KBIN_FROZSL_FILE, \
      KBIN_FROZSL_FILE_NAME, KBIN_FROZSL_STRUCTURE_MODE_FILE, KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP, KBIN_FROZSL_STRUCTURE_STRING, KBIN_FROZSL_DIELECTRIC_MODE_FILE,                                    \
      KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP, KBIN_FROZSL_DIELECTRIC_ZEFF, KBIN_FROZSL_DIELECTRIC_STRING, KBIN_PHONONS_CALCULATION_APL, KBIN_PHONONS_CALCULATION_QHA, KBIN_PHONONS_CALCULATION_AAPL,    \
      KBIN_PHONONS_CALCULATION_AGL, KBIN_PHONONS_CALCULATION_AEL, KBIN_PHONONS_CALCULATION_FROZSL, KBIN_PHONONS_CALCULATION_FROZSL_output, KBIN_PHONONS_CALCULATION_FROZSL_poscars
};

// --------------------------------------------------------------------------
// general flags for vasp mode
class xstructure; // prototype of structure, just to compile
class _vflags : public JsonSerializable<_vflags> {
public:
  // trivial constructurs/destuctors/operators
  _vflags(); // default, just allocate
  ~_vflags(); // kill everything
  _vflags(const _vflags& b); // constructor copy
  const _vflags& operator=(const _vflags& b); // copy
  void clear(); // clear
  // CONTENT
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in VASP MODE
  aurostd::xoption AFLOW_SYSTEM; // ME20181121
  int KBIN_VASP_RUN_NRELAX;
  aurostd::xoption KBIN_VASP_RUN;
  // GENERATE, STATIC, KPOINTS, RELAX, RELAX_STATIC, RELAX_STATIC_BANDS, RELAX_STATIC_DIELECTRIC, STATIC_BANDS, STATIC_DIELECTRIC
  aurostd::xoption KBIN_VASP_REPEAT; // REPEAT_STATIC REPEAT_STATIC_BANDS REPEAT_BANDS REPEAT_DIELECTRIC REPEAT_DELSOL
  aurostd::xoption KBIN_VASP_FORCE_OPTION_NOTUNE; // NOTUNE
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO; // SYSTEM_AUTO
  aurostd::xoption KBIN_VASP_FORCE_OPTION_RELAX_MODE; // RELAX_MODE  forces/energy
  aurostd::xoption KBIN_VASP_FORCE_OPTION_RELAX_TYPE;
  // RELAX_TYPE  STATIC, ALL, IONS, CELL_SHAPE, CELL_VOLUME, IONS_CELL_VOLUME
  aurostd::xoption KBIN_VASP_FORCE_OPTION_PREC; // PREC
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ALGO; // ALGO
  aurostd::xoption KBIN_VASP_FORCE_OPTION_METAGGA; // METAGGA
  aurostd::xoption KBIN_VASP_FORCE_OPTION_IVDW; // IVDW
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ABMIX; // ABMIX
  aurostd::xoption KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS; // AUTO_PSEUDOPOTENTIALS
  // ENMAX_MULTIPLY
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL; // isentry and content_int
  // NBANDS
  aurostd::xoption KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL; // isentry and content_int
  bool KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry;
  // POTIM
  aurostd::xoption KBIN_VASP_FORCE_OPTION_POTIM_EQUAL; // isentry and content_double
  // PSTRESS
  aurostd::xoption KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL; // isentry and content_double
  // EDIFFG
  aurostd::xoption KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL; // isentry and content_double
  // NELM
  aurostd::xoption KBIN_VASP_FORCE_OPTION_NELM_EQUAL; // isentry and content_double //CO20200624
  // NELM_STATIC
  aurostd::xoption KBIN_VASP_FORCE_OPTION_NELM_STATIC_EQUAL; // isentry and content_double //CO20200624
  // ISMEAR
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL; // isentry and content_double //CO20181129
  // SIGMA
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL; // isentry and content_double //CO20181129
  // ISMEAR_STATIC
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ISMEAR_STATIC_EQUAL; // isentry and content_double //CO20181129
  // SIGMA_STATIC
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SIGMA_STATIC_EQUAL; // isentry and content_double //CO20181129
  // ISMEAR_BANDS
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ISMEAR_BANDS_EQUAL; // isentry and content_double //CO20181129
  // SIGMA_BANDS
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SIGMA_BANDS_EQUAL; // isentry and content_double //CO20181129
  // RWIGS
  bool KBIN_VASP_FORCE_OPTION_RWIGS_STATIC;
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SKIP_NOMIX; // SKIP_NOMIX
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SPIN; // SPIN
  bool KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;
  bool KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2;
  // xoption KBIN_VASP_FORCE_OPTION_TRISTATE;      //  SYM
  aurostd::xoption KBIN_VASP_FORCE_OPTION_BADER; // BADER=ON | OFF | NONE
  aurostd::xoption KBIN_VASP_FORCE_OPTION_ELF; // ELF=ON | OFF | NONE
  aurostd::xoption KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM; // AUTO_MAGMOM
  aurostd::xoption KBIN_VASP_FORCE_OPTION_SYM; // SYM
  aurostd::xoption KBIN_VASP_FORCE_OPTION_WAVECAR; // WAVECAR
  aurostd::xoption KBIN_VASP_FORCE_OPTION_CHGCAR; // CHGCAR
  aurostd::xoption KBIN_VASP_FORCE_OPTION_CHGCAR_FILE; // ME20191028
  aurostd::xoption KBIN_VASP_FORCE_OPTION_LSCOUPLING; // LSCOUPLING
  aurostd::xoption KBIN_VASP_FORCE_OPTION_LDAU0; // LDAU0
  aurostd::xoption KBIN_VASP_FORCE_OPTION_LDAU1; // LDAU1
  aurostd::xoption KBIN_VASP_FORCE_OPTION_LDAU2; // LDAU2
  aurostd::xoption KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC; // LDAU_ADIABATIC
  aurostd::xoption KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF; // LDAU_CUTOFF
  std::string KBIN_VASP_LDAU_SPECIES;
  std::string KBIN_VASP_LDAU_PARAMETERS;
  bool KBIN_VASP_LDAU_AFLOW_AUTO_flag;
  // FORCE_OPTION
  aurostd::xoption KBIN_VASP_FORCE_OPTION_TYPE; // TYPE
  bool KBIN_VASP_FORCE_OPTION_NSW_EQUAL;
  int KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE;

  aurostd::xoption KBIN_VASP_FORCE_OPTION_IGNORE_AFIX; // AFIX
  // xoption kopts;
  aurostd::xoption KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL; // CONVERT_UNIT_CELL
  aurostd::xoption KBIN_VASP_FORCE_OPTION_VOLUME; // EQUAL_EQUAL, MULTIPLY_EQUAL,PLUS_EQUAL
  aurostd::xoption KBIN_VASP_FORCE_OPTION_KPOINTS; // KPOINTS
  aurostd::xoption KBIN_VASP_INCAR_MODE; // EXPLICIT, EXPLICIT_START_STOP, IMPLICIT, EXTERNAL;
  // RELAX
  aurostd::xoption KBIN_VASP_KPOINTS_MODE; // EXPLICIT, EXPLICIT_START_STOP, IMPLICIT, EXTERNAL;
  aurostd::xoption KBIN_VASP_KPOINTS_KMODE; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_KPPRA; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_KSCHEME; // isentry and content_string
  aurostd::xoption KBIN_VASP_KPOINTS_KSHIFT; // isentry and content_string
  // STATIC
  aurostd::xoption KBIN_VASP_KPOINTS_STATIC_KMODE; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_STATIC_KPPRA; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_STATIC_KSCHEME; // isentry and content_string
  aurostd::xoption KBIN_VASP_KPOINTS_STATIC_KSHIFT; // isentry and content_string
  // PHONONS
  aurostd::xoption KBIN_VASP_KPOINTS_PHONONS_KPPRA; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_PHONONS_KSCHEME; // isentry and content_string
  aurostd::xoption KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY; // EVEN ODD
  aurostd::xoption KBIN_VASP_KPOINTS_PHONONS_GRID; // ME20200427
  aurostd::xoption KBIN_VASP_KPOINTS_PHONONS_SHIFT; // ME20200427
  // BANDS
  aurostd::xoption KBIN_VASP_KPOINTS_BANDS_LATTICE;
  aurostd::xoption KBIN_VASP_KPOINTS_BANDS_GRID; // CO20210805
  // DIELECTRIC
  aurostd::xoption KBIN_VASP_KPOINTS_DIELECTRIC_KMODE; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_DIELECTRIC_KPPRA; // isentry and content_int
  aurostd::xoption KBIN_VASP_KPOINTS_DIELECTRIC_KSCHEME; // isentry and content_string
  aurostd::xoption KBIN_VASP_KPOINTS_DIELECTRIC_KSHIFT; // isentry and content_string
  bool KBIN_VASP_WRITE_KPOINTS;

  aurostd::xoption KBIN_VASP_POSCAR_MODE;
  // EXPLICIT, EXPLICIT_START_STOP, EXPLICIT_START_STOP_POINT, IMPLICIT, EXTERNAL;
  std::vector<std::string> KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING;
  std::vector<xstructure> KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE;
  bool KBIN_VASP_INCAR_VERBOSE; // VERBOSITY
  aurostd::xoption KBIN_VASP_INCAR_FILE; // KEYWORD, SYSTEM_AUTO, FILE, COMMAND
  std::stringstream KBIN_VASP_INCAR_EXPLICIT; // ME20181127
  std::stringstream KBIN_VASP_INCAR_EXPLICIT_START_STOP; // ME20181127
  aurostd::xoption KBIN_VASP_KPOINTS_FILE; // KEYWORD, FILE, COMMAND
  std::stringstream KBIN_VASP_KPOINTS_EXPLICIT; // ME20181127
  std::stringstream KBIN_VASP_KPOINTS_EXPLICIT_START_STOP; // ME20181127
  aurostd::xoption KBIN_VASP_POSCAR_FILE; // KEYWORD, PROTOTYPE, FILE, COMMAND
  aurostd::xoption KBIN_VASP_POSCAR_FILE_VOLUME; // EQUAL_EQUAL, MULTIPLY_EQUAL PLUS_EQUAL
  aurostd::xoption KBIN_VASP_POTCAR_MODE; // EXPLICIT, IMPLICIT, EXTERNAL;
  aurostd::xoption KBIN_VASP_POTCAR_FILE; // KEYWORD, SYSTEM_AUTO, PREFIX, SUFFIX, FILE, COMMAND, WRITE
  std::stringstream KBIN_VASP_POTCAR_EXPLICIT; // CO20181226

  // JSON load/dump
protected:
  [[nodiscard]] aurostd::JSON::object serialize() const override;
  _vflags deserialize(const aurostd::JSON::object& jo) override;
  [[nodiscard]] std::string getJsonID() const override { return "_vflags"; }

private: //
  void free(); // free space
  void copy(const _vflags& b); //
};

// SERIALIZATION MEMBERS ignoring all xoption, stringstream members
#define JSON_vflags_MEMBERS                                                                                                                                                                      \
  KBIN_VASP_RUN_NRELAX, KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry, KBIN_VASP_FORCE_OPTION_RWIGS_STATIC, KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1, KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2, \
      KBIN_VASP_LDAU_SPECIES, KBIN_VASP_LDAU_PARAMETERS, KBIN_VASP_LDAU_AFLOW_AUTO_flag, KBIN_VASP_FORCE_OPTION_NSW_EQUAL, KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE, KBIN_VASP_WRITE_KPOINTS,      \
      KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING, KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE, KBIN_VASP_INCAR_VERBOSE

// --------------------------------------------------------------------------
// general flags for aims mode
class _aimsflags {
public:
  _aimsflags();
  ~_aimsflags();
  _aimsflags(const _aimsflags& b);
  const _aimsflags& operator=(const _aimsflags& b);
  void clear();
  // CONTENT
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in AIMS MODE
  aurostd::xoption KBIN_AIMS_FORCE_OPTION_NOTUNE;
  aurostd::xoption KBIN_AIMS_RUN;
  aurostd::xoption KBIN_AIMS_GEOM_MODE;
  aurostd::xoption KBIN_AIMS_GEOM_FILE;
  std::vector<std::string> KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING;
  std::vector<xstructure> KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE;
  aurostd::xoption KBIN_AIMS_GEOM_FILE_VOLUME;
  aurostd::xoption KBIN_AIMS_FORCE_OPTION_VOLUME;
  aurostd::xoption KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL;
  aurostd::xoption KBIN_AIMS_CONTROL_MODE;
  aurostd::xoption KBIN_AIMS_CONTROL_FILE;
  bool KBIN_AIMS_CONTROL_VERBOSE; // VERBOSITY
private:
  void free(); // free space
  void copy(const _aimsflags& b); //
};

// --------------------------------------------------------------------------
// general flags for alien mode
class _alienflags {
public:
  // trivial constructurs/destuctors/operators
  _alienflags(); // default, just allocate
  ~_alienflags(); // kill everything
  _alienflags(const _alienflags& b); // constructor copy
  const _alienflags& operator=(const _alienflags& b); // copy
  void clear(); // clear
  // CONTENT
  bool KBIN_ALIEN_COMMAND_BINARY_FLAG;
  std::string KBIN_ALIEN_COMMAND_BINARY_VALUE;
  bool KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG;
  // in this struct we put all the flags which are used on LOCAL DIRECTORIES in ALIEN MODE
  bool KBIN_ALIEN_FORCE_OPTION_NOTUNE;
  bool KBIN_ALIEN_FORCE_OPTION_SOMETHING; // SOMETHING

  bool KBIN_ALIEN_INPUT_MODE_EXPLICIT;
  bool KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP;
  bool KBIN_ALIEN_INPUT_MODE_IMPLICIT;
  bool KBIN_ALIEN_INPUT_MODE_EXTERNAL;
  bool KBIN_ALIEN_INPUT_FILE;
  bool KBIN_ALIEN_INPUT_FILE_FILE_FLAG;
  std::string KBIN_ALIEN_INPUT_FILE_FILE_VALUE;
  bool KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG;
  std::string KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE;
  bool KBIN_ALIEN_INPUT_MODE_INPUT_FLAG;
  std::string KBIN_ALIEN_INPUT_MODE_INPUT_VALUE;
  bool KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG;
  std::string KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE;

private: //
  void free(); // free space
  void copy(const _alienflags& b); //
};

// --------------------------------------------------------------------------
// general container for any set of flags
class _xflags {
public:
  _xflags();
  _xflags(_vflags& vflags);
  _xflags(_aimsflags& aimsflags);
  _xflags(_alienflags& alienflags);
  ~_xflags();
  _xflags(const _xflags& b);
  const _xflags& operator=(const _xflags& b);
  void clear();
  bool AFLOW_MODE_VASP;
  _vflags vflags;
  bool AFLOW_MODE_AIMS;
  _aimsflags aimsflags;
  bool AFLOW_MODE_ALIEN;
  _alienflags alienflags;
  void setVFlags(_vflags& vflags);
  void setAIMSFlags(_aimsflags& aimsflags);
  // add qe and others here eventually
  void setALIENFlags(_alienflags& alienflags);

private:
  void free(); // free space
  void copy(const _xflags& b); //
};

// for queue
class _xqsub {
public:
  // trivial constructurs/destuctors/operators
  _xqsub(); // default, just allocate
  ~_xqsub(); // kill everything
  _xqsub(const _xqsub& b); // constructor copy
  const _xqsub& operator=(const _xqsub& b); // copy
  void clear(); // clear
  std::stringstream QSUB; //
  std::stringstream QSUB_orig; //
  bool QSUB_generated; //
  bool QSUB_changed; //
private: //
  void free(); // free space
  void copy(const _xqsub& b); //
};

// for a vasp run
class _xvasp {
public:
  // trivial constructurs/destuctors/operators
  _xvasp(); // default, just allocate
  ~_xvasp(); // kill everything
  _xvasp(const _xvasp& b); // constructor copy
  const _xvasp& operator=(const _xvasp& b); // copy
  void clear(); // clear
  // aflow_xatom.cpp contains the code
  // CONTENT
  xstructure str;
  std::string Directory;
  std::string AnalyzeLabel;
  _xqsub xqsub;
  aurostd::xoption aopts;
  // --------------------------------
  // VASP INPUT CONTENT
  std::stringstream POSCAR;
  std::stringstream POSCAR_orig;
  uint POSCAR_index;
  std::stringstream INCAR;
  std::stringstream INCAR_orig;
  std::stringstream KPOINTS;
  std::stringstream KPOINTS_orig;
  std::stringstream POTCAR;
  std::stringstream POTCAR_orig;
  std::string POTCAR_TYPE;
  bool POTCAR_TYPE_DATE_PRINT_flag; // to add Zr:POT_PAW:01Apr2001 to directory...
  bool POTCAR_TYPE_PRINT_flag; // to add Zr:POT_PAW to directory... (no date) //CO20181226
  double POTCAR_ENMAX;
  double POTCAR_ENMIN;
  bool POTCAR_PAW;
  std::stringstream POTCAR_POTENTIALS;
  std::deque<std::string> POTCAR_AUID; // md5sum one line by one line of the POTCAR
  // --------------------------------
  // VASP INPUT CONTENT
  std::stringstream OUTCAR; // OUTPUT
  std::stringstream CONTCAR; // OUTPUT
  std::stringstream OSZICAR; // OUTPUT
  // --------------------------------
  // QE INPUT CONTENT
  std::stringstream QE_GEOM; // FUTURE
  std::stringstream QE_GEOM_orig; // FUTURE
  bool QE_GEOM_generated; // FUTURE
  bool QE_GEOM_changed; // FUTURE
  uint QE_GEOM_index; // FUTURE
  // --------------------------------
  // ABINIT INPUT CONTENT
  // --------------------------------
  // AIMS INPUT CONTENT
  // --------------------------------
  int NCPUS;
  int NRELAX; // job to do         // -1 (static) 0(error) 1,2,3,4.... (relaxes)  -2(run kpoints)
  int NRELAXING; // job doing (to monitor odd/even)
  // --------------------------------
  // content for AVASP generation
  bool AVASP_aflowin_only_if_missing; //
  bool AVASP_arun; // ME20181019
  std::string AVASP_arun_mode; // ME20181019
  std::string AVASP_arun_runname; // ME20181019
  aurostd::xoption aplopts; // ME20181025
  aurostd::xoption aaplopts; // ME20181025
  aurostd::xoption qhaopts; // AS20200302
  std::string AVASP_dirbase;
  std::string AVASP_libbase;
  std::string AVASP_label;
  std::string AVASP_parameters;
  std::string AVASP_pocc_parameters; // CO20181226
  std::string AVASP_pocc_tol; // CO20181226
  double AVASP_volume_in;
  std::string AVASP_potential;
  bool AVASP_alpha_fix;
  uint AVASP_prototype_mode;
  bool AVASP_prototype_from_library_;
  bool AVASP_directory_from_library_;
  int AVASP_value_NSW;
  int AVASP_value_KPPRA;
  std::string AVASP_KSCHEME;
  int AVASP_value_KPPRA_STATIC;
  int AVASP_value_KPPRA_DIELECTRIC;
  std::string AVASP_STATIC_KSCHEME;
  std::string AVASP_DIELECTRIC_KSCHEME;
  std::string AVASP_KPOINTS; // ME20181226
  std::string AVASP_flag_PRECISION_scheme;
  std::string AVASP_flag_ALGO_scheme;
  std::string AVASP_flag_ABMIX_scheme;
  aurostd::xoption AVASP_flag_TYPE; // TYPE
  std::string AVASP_LDAU_PARAMETERS_STRING;
  double AVASP_LDAU_PARAMETERS_UJSUM;
  std::stringstream AVASP_EXTRA_INCAR;
  std::vector<std::string> AVASP_INCAR_KEYWORD; // ME20181226
  std::stringstream AVASP_INCAR_EXPLICIT_START_STOP; // ME20181226
  std::vector<std::string> AVASP_KPOINTS_KEYWORD; // ME20181226
  std::stringstream AVASP_KPOINTS_EXPLICIT_START_STOP; // ME20181226
  std::vector<std::string> AVASP_POTCAR_KEYWORD; // ME20181226
  bool AVASP_flag_MPI;
  bool AVASP_flag_RUN_RELAX;
  bool AVASP_flag_RUN_RELAX_STATIC;
  bool AVASP_flag_RUN_RELAX_STATIC_BANDS;
  bool AVASP_flag_RUN_RELAX_STATIC_DIELECTRIC;
  bool AVASP_flag_RUN_STATIC_BANDS;
  bool AVASP_flag_RUN_STATIC_DIELECTRIC;
  std::string AVASP_path_BANDS;
  uint AVASP_value_BANDS_GRID;
  bool AVASP_flag_RUN_STATIC;
  bool AVASP_flag_GENERATE;
  // -------------------------------- FUNCTIONS
  double GetZVAL();
  double GetCellAtomZVAL(std::string mode); // CELL ATOM
  double GetPOMASS();
  double GetCellAtomPOMASS(std::string mode); // CELL ATOM
private: //
  void free(); // free space
  void copy(const _xvasp& b); //
};

// for an aims run
class _xaims {
public:
  _xaims();
  ~_xaims();
  _xaims(const _xaims& b);
  const _xaims& operator=(const _xaims& b);
  void clear();
  uint GEOM_index;
  xstructure str;
  std::string Directory;
  _xqsub xqsub;
  aurostd::xoption aopts;
  int NCPUS;
  // --------------------------------
  // AIMS INPUT CONTENT
  std::stringstream CONTROL;
  std::stringstream CONTROL_orig;
  bool CONTROL_generated;
  bool CONTROL_changed;
  std::string CONTROL_FILE_NAME;
  std::stringstream GEOM;
  std::stringstream GEOM_orig;
  bool GEOM_generated;
  bool GEOM_changed;
  std::string GEOM_FILE_NAME;
  std::string OUTPUT_FILE_NAME;

private: //
  void free(); // free space
  void copy(const _xaims& b); //
};

// for a alien run
class _xalien {
public:
  // trivial constructurs/destuctors/operators
  _xalien(); // default, just allocate
  ~_xalien(); // kill everything
  _xalien(const _xalien& b); // constructor copy
  const _xalien& operator=(const _xalien& b); // copy
  void clear(); // clear
  // aflow_xatom.cpp contains the code
  // CONTENT
  std::string Directory;
  _xqsub xqsub;
  std::stringstream INPUT;
  std::stringstream INPUT_orig;
  bool INPUT_generated;
  bool INPUT_changed;
  std::string INPUT_FILE_NAME;
  std::string OUTPUT_FILE_NAME;
  // ----------------
  int NCPUS;
  int NRELAX; // -1 (static) 0(error) 1,2,3,4.... (relaxes)  -2(run kpoints)
private: //
  void free(); // free space
  void copy(const _xalien& b); //
};

// for a generic run
class _xinput {
public:
  _xinput();
  _xinput(_xvasp& xvasp);
  _xinput(_xaims& xaims);
  _xinput(_xalien& xalien);
  ~_xinput();
  _xinput(const _xinput& b);
  const _xinput& operator=(const _xinput& b);
  void clear();
  bool AFLOW_MODE_VASP;
  _xvasp xvasp;
  bool AFLOW_MODE_AIMS;
  _xaims xaims;
  bool AFLOW_MODE_ALIEN;
  _xalien xalien;
  void setXVASP(_xvasp& xvasp);
  void setXAIMS(_xaims& xaims);
  void setXALIEN(_xalien& xalien);
  xstructure& getXStr();
  std::string& getDirectory();
  void setXStr(const xstructure& str, bool set_all = false);
  void setDirectory(const std::string Directory, bool set_all = false);

private:
  void free();
  void copy(const _xinput& b);
};

#endif // AFLOW_XCLASSES_H
