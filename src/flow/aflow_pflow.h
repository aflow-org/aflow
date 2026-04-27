// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_PFLOW_H_
#define _AFLOW_PFLOW_H_

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow_defs.h"
#include "aflow_support_types.h"
#include "aflowlib/aflowlib_web_interface.h"
#include "flow/aflow_xclasses.h"
#include "structure/aflow_xstructure.h"

// pflow prototypes
// aflow_pflow.cpp
#define kwarning "[WARNING] aflow: "
#define kintro " "

// LOGGER MODES
static const char _LOGGER_ERROR_ = 'E';
static const char _LOGGER_WARNING_ = 'W';
static const char _LOGGER_COMPLETE_ = 'C';
static const char _LOGGER_OPTION_ = 'O';
static const char _LOGGER_RAW_ = 'R';
static const char _LOGGER_MESSAGE_ = 'M';
static const char _LOGGER_NOTICE_ = 'N';

#define XRAY_THETA_TOL 1e-5;                                  // CO20190409

//[MOVED to aflow.h]// LOADENTRIES DEFAULTS
//[MOVED to aflow.h]static const uint _AFLOW_LIB_MAX_ = 10;  //LIB11 does not exist yet, modify accordingly

// aflow_pflow_main.cpp
namespace pflow {
  int main(std::vector<std::string>& argv, std::vector<std::string>& cmds);
}

namespace pflow {
  bool CheckCommands(std::vector<std::string>, const std::vector<std::string>& cmds);
  void CheckIntegritiy();
  std::string Intro_pflow(std::string x);
} // namespace pflow
std::string AFLOW_PrototypesIcsdReadme();
namespace pflow {
  xstructure XXX(std::istream& input, const int& iatom);
  bool XXX(std::vector<std::string>, std::istream& input);
} // namespace pflow
namespace pflow {
  xstructure ABCCAR(std::istream& input);
  void ACE(std::istream& input);
  bool AddSpinToXstructure(xstructure& a, std::vector<double>& vmag); // DX20170927 - Magnetic symmetry
  bool AddSpinToXstructure(xstructure& a, std::vector<aurostd::xvector<double>>& vmag_noncoll); // DX20171205 - Magnetic symmetry (non-collinear)
  void AGROUP(_aflags& aflags, std::istream& input);
  void AGROUP2(std::istream& input);
  void AGROUP2m(std::istream& input);
  xstructure ALPHABETIC(std::istream& input);
  std::string ALPHACompound(const std::string& options);
  std::string ALPHASpecies(const std::string& options);
  std::string AFLOWIN(std::istream& input);
  void ANGLES(const std::string& options, std::istream& input);
  std::string ATOMSMAX(const std::string& options, std::istream& input);
  void BANDS(const std::string& options, std::istream& input);
  void BANDGAP(aurostd::xoption& vpflow, std::ostream& oss = std::cout); // CAMILO  //CO20171006
  void BANDGAP_DOS(aurostd::xoption& vpflow, std::ostream& oss = std::cout); // CAMILO  //CO20171006  //CO20191110
  void BANDSTRUCTURE(_aflags& aflags);
  std::string BZDirectionsLATTICE(const std::string& options);
  std::string BZDirectionsSTRUCTURE(std::istream& input, aurostd::xoption& vpflow); // DX20181102 - add options
  void CAGES(_aflags& aflags, const std::string& options, std::istream& input);
  // DX+CO START
  bool PerformFullSymmetry(xstructure& a);
  bool PerformFullSymmetry(xstructure& a, std::ofstream& FileMESSAGE, const std::string& directory, _kflags& kflags, bool osswrite, std::ostream& oss, std::string format = "txt");  // ME20200224
  bool PerformFullSymmetry(xstructure& a, std::ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, bool osswrite, std::ostream& oss, std::string format = "txt");
  bool PerformFullSymmetry(xstructure& a,
                           double& tolerance,
                           bool no_scan,
                           bool force_perform,
                           std::ofstream& FileMESSAGE,
                           _aflags& aflags,
                           _kflags& kflags,
                           bool osswrite,
                           std::ostream& oss,
                           std::string format = "txt");
  void ProcessAndAddSpinToXstructure(xstructure& a, const std::string& magmom_info); // DX20190801
  void defaultKFlags4SymWrite(_kflags& kflags, bool write = true);
  void defaultKFlags4SymCalc(_kflags& kflags, bool calc = true);
  bool CalculateFullSymmetry(std::istream& input, aurostd::xoption& vpflow, std::ostream& oss = std::cout);
  bool CalculateFullSymmetry(_aflags& aflags, _kflags& kflags, xstructure& a, aurostd::xoption& vpflow, bool osswrite, std::ostream& oss = std::cout);
  bool fixEmptyAtomNames(xstructure& xstr, bool force_fix = false);  // force_fix=true if you want to override what is already in species
  // DX+CO END
  xstructure CART(std::istream& input);
  xstructure CORNERS(std::istream& input);
  void ChangeSuffix(const std::string& options);
  std::string CHGDIFF(aurostd::xoption vpflow);
  bool CHGDIFF(const std::string& chgcar1_file, const std::string& chgcar2_file, const std::string& output_File, std::ostream& oss = std::cout);
  // DX+CO START
  void CHGINT(std::vector<std::string>);
  // DX+CO END
  std::string CHGSUM(aurostd::xoption vpflow);
  bool CHGSUM(const std::string& chgcar_in1, const std::string& chgcar_in2, std::ostream& oss = std::cout);
  bool CHGSUM(std::string& species_header, const std::string& chgcar_in1, const std::string& chgcar_in2, const std::string& output_file, std::ostream& oss = std::cout);
  bool CHGSUM(const std::string& chgcar_in1, const std::string& chgcar_in2, const std::string& output_file, std::ostream& oss = std::cout);
  bool CHGSUM(const std::vector<std::string>& chgcar_files, std::ostream& oss = std::cout);
  bool CHGSUM(const std::vector<std::string>& chgcar_files, const std::string& output, std::ostream& oss = std::cout);
  bool CHGSUM(std::string& species_header, const std::vector<std::string>& chgcar_files, std::ostream& oss = std::cout);
  bool CHGSUM(std::string& species_header, const std::vector<std::string>& chgcar_files, const std::string& output_file, std::ostream& oss = std::cout);
  void CIF(std::istream& input, aurostd::xoption& vpflow);
  void getCIFOptions(xstructure& a, aurostd::xoption& vpflow); // DX20250319
  void CLAT(const std::string& options);
  void CLEAN(std::vector<std::string>);
  void CLEANALL(std::istream& input);
  void CMPSTR(std::vector<std::string>);
  void COMPARE(const std::string& options);
  bool DATA(std::istream& input, const aurostd::xoption& vpflow, const std::string& smode = "DATA", std::ostream& oss = std::cout); // DX20170901 - SGDATA + JSON //DX20210302 - added const to vpflow
  void DEBYE(const std::string& options);
  void DISP(const std::string& options, std::istream& input);
  void DIST(const std::string& options, std::istream& input);
  void DYNADIEL(std::vector<std::string>& argv); // CAMILO
  void EDOS(std::vector<std::string>);
  void EFFMASS(std::vector<std::string>& argv, std::ostream& oss = std::cout); // CAMILO
  void EIGCURV(const std::string& options, std::ostream& oss = std::cout); // CAMILO
  std::string EQUIVALENT(_aflags& aflags, std::istream& input, aurostd::xoption& vpflow);
  void EWALD(const std::string& options, std::istream& input);
  std::string EXTRACT_xcar(_aflags& aflags, std::vector<std::string>, std::string, std::string);
  std::string EXTRACT_Symmetry(_aflags& aflags, std::vector<std::string>);
  void FGROUP(_aflags& aflags, std::istream& input);
  bool FIXBANDS(_aflags& aflags, std::string opts);
  void FINDSYM(aurostd::xoption& vpflow, uint mode, std::istream& input);
  xstructure FRAC(std::istream& input);
  std::string FROZSL_VASPSETUP(std::vector<std::string> argv, int mode);
  std::string FROZSL_ANALYZE(std::istream& input);
  std::string FROZSL_INPUT();
  std::string FROZSL_OUTPUT();
  std::string GEOMETRY(std::istream& input); // CO20191110
  bool GetCollinearMagneticInfo(uint num_atoms, const std::string& magmom_info, std::vector<double>& vmag); // DX20170927 - Magnetic symmetry //DX20191107 - int to uint
  bool GetNonCollinearMagneticInfo(uint num_atoms, const std::string& magmom_info, std::vector<aurostd::xvector<double>>& vmag_noncoll); // DX20171205 - Magnetic symmetry non-collinear //DX20191107 - int to uint
  void GLASS_FORMING_ABILITY(aurostd::xoption& vpflow); // DF20190329
  void ATOMIC_ENVIRONMENT(const aurostd::xoption& vpflow); // HE20210331
  void GULP(std::istream& input);
  void HKL(const std::string& options, _aflags& aflags, std::istream& input);
  void HKLSearch(const std::string& options, _aflags& aflags, std::istream& input, const std::string& smode);
  bool setPOCCTOL(xstructure& xstr, const std::string& pocc_tol_string); // CO20181226
  bool POCC_COMMAND_LINE(aurostd::xoption& vpflow, std::istream& input, std::ostream& oss = std::cout); // CO20181226
  void ICSD_2WYCK(std::istream& input, bool SOF);
  void ICSD(std::vector<std::string> argv, std::istream& input);
  xstructure IDENTICAL(std::istream& input);
  xstructure INCELL(std::istream& input);
  xstructure INCOMPACT(std::istream& input);
  void INTPOL(const std::string& options);
  xstructure INWS(std::istream& input);
  void JMOL(const std::string& options, std::istream& input);
  void KBAND(std::vector<std::string>);
  xstructure INFLATE_LATTICE(const std::string& options, std::istream& input);
  xstructure INFLATE_VOLUME(const std::string& options, std::istream& input);
  void KPATH(std::istream& input, bool WWW);
  void KPATH(std::istream& input, double grid, bool WWW);
  xstructure KPOINTS(const std::string& options, std::istream& input, std::ostream& oss = std::cout);
  xstructure KPOINTS_DELTA(aurostd::xoption& vpflow, std::istream& input, std::ostream& oss = std::cout);
  void JOINSTRLIST(std::vector<std::string>);
  void MAKESTRLIST(std::vector<std::string>);
  xstructure LATTICEREDUCTION(std::istream& input);
  std::string LATTICE_TYPE(std::istream& input, aurostd::xoption& vpflow); // DX20200820 - added vpflow
  std::string LATTICE_LATTICE_TYPE(std::istream& input, aurostd::xoption& vpflow); // DX20200820
  std::string listPrototypeLabels(aurostd::xoption& vpflow); // DX20181004
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // START - all relevent functions for loading entries here
  // Added by Corey Oses - May 2017
  // load entries is heavily overloaded, mostly to accommodate entries separated as
  // std::vector<std::vector<std::vector<> > > entries (unaries vs. binaries, then species-specific, good for convex hull),
  // std::vector<std::vector<> > entries (unaries vs. binaries), OR
  // std::vector<> entries (all together)
  ////////////////////////////////////////////////////////////////////////////////
  std::string arity_string(uint arity, bool capital = false, bool plural = false); // CO20180329
  // loadEntries
  bool loadEntries(std::vector<std::string>& velements, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow,
                   std::vector<std::string>& velements,
                   std::string server,
                   std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& entries,
                   std::ofstream& FileMESSAGE,
                   std::ostream& oss = std::cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadEntries
  bool loadEntries(std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadEntries
  bool loadEntries(std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadEntries(std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadEntries(aurostd::xoption& vpflow, std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  ////////////////////////////////////////////////////////////////////////////////
  bool loadFromCommon(aurostd::xoption& vpflow);
  ////////////////////////////////////////////////////////////////////////////////
  // load and merging LIBX
  bool loadAndMergeLIBX(std::string combination, std::string LIB, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(std::string _combination, std::string LIB, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(std::vector<std::string>& combination, std::string LIB, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(std::vector<std::string>& combination, std::string LIB, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, std::string combination, std::string LIB, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow,
                        std::string _combination,
                        std::string LIB,
                        std::string server,
                        std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries,
                        std::ofstream& FileMESSAGE,
                        std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, std::vector<std::string>& combination, std::string LIB, std::string server, std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries, std::ostream& oss = std::cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow,
                        std::vector<std::string>& combination,
                        std::string LIB,
                        std::string server,
                        std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>>& naries,
                        std::ofstream& FileMESSAGE,
                        std::ostream& oss = std::cout);
  ////////////////////////////////////////////////////////////////////////////////
  uint SubLayersRestAPILS(const std::string& url, std::vector<std::string>& vsuburl); // CO20200731
  ////////////////////////////////////////////////////////////////////////////////
  // loadLIBX std::string elements
  bool loadLIBX(std::string LIB, std::string elements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::string elements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::string elements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::string elements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  // loadLIBX std::vector elements
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<aflowlib::_aflowlib_entry>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadLIBX std::string elements, organized by -naries
  bool loadLIBX(std::string LIB, std::string elements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::string elements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::string elements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::string elements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::string elements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  // loadLIBX_nested std::vector elements
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ostream& oss = std::cout);
  bool loadLIBX(aurostd::xoption& vpflow, std::string LIB, std::vector<std::string>& velements, std::string server, std::vector<std::vector<aflowlib::_aflowlib_entry>>& entries, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  ////////////////////////////////////////////////////////////////////////////////
  // get elemental combinations (recursively)
  std::vector<std::vector<std::string>> elementalCombinations(const std::vector<std::string>& velements, uint nary);
  ////////////////////////////////////////////////////////////////////////////////
  // easy way to think about it:  do compounds belong to the hull?
  bool compoundsBelong(const std::vector<std::string>& velements,
                       const std::string& input,
                       std::ostream& oss = std::cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const std::vector<std::string>& velements,
                       const std::string& input,
                       std::vector<std::string>& input_velements,
                       std::vector<double>& input_vcomposition,
                       std::ostream& oss = std::cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const std::vector<std::string>& velements,
                       const std::string& input,
                       std::ofstream& FileMESSAGE,
                       std::ostream& oss = std::cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const std::vector<std::string>& velements,
                       const std::string& input,
                       std::vector<std::string>& input_velements,
                       std::vector<double>& input_vcomposition,
                       std::ofstream& FileMESSAGE,
                       std::ostream& oss = std::cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const std::vector<std::string>& velements, const std::vector<std::string>& elements, std::ostream& oss = std::cout, bool sort_elements = false);
  bool compoundsBelong(const std::vector<std::string>& velements, const std::vector<std::string>& elements, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool sort_elements = false);
  ////////////////////////////////////////////////////////////////////////////////
  // loads xstructures
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, std::ostream& oss = std::cout, bool relaxed_only = true, std::string path = "", bool is_url_path = false);
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool relaxed_only = true, std::string path = "", bool is_url_path = false);
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, std::vector<std::string>& structure_files, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, bool relaxed_only = true, std::string path = "", bool is_url_path = false); // DX20200224
  bool loadXstructureLibEntry(const aflowlib::_aflowlib_entry& entry, xstructure& new_structure); // HE20220606
  ////////////////////////////////////////////////////////////////////////////////
  // sets default flags
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, std::ostream& oss = std::cout, std::string input = "A", bool entry_output = true, bool silent = false);
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout, std::string input = "A", bool entry_output = true, bool silent = false);
  ////////////////////////////////////////////////////////////////////////////////
  bool prototypeMatch(std::string proto_database, std::string proto_search); // smarter than == for prototype matches, deals with 549 vs. 549.bis
  ////////////////////////////////////////////////////////////////////////////////
  // Added by Corey Oses - May 2017
  // END - all relevent functions for loading entries here
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // START - added by Corey Oses - May 2017
  // effectively logs EVERYTHING, deals with std::cout and logger
  void updateProgressBar(unsigned long long int current, unsigned long long int end, std::ostream& oss = std::cout);
  void logger(const std::string& filename, const std::string& function_name, std::stringstream& message, const char& type, std::ostream& oss = std::cout, bool silent = false, const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              std::stringstream& message,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              std::stringstream& message,
              std::ofstream& FileMESSAGE,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              std::stringstream& message,
              const std::string& directory,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              std::stringstream& message,
              const std::string& directory,
              std::ofstream& FileMESSAGE,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              std::stringstream& message,
              const _aflags& aflags,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              std::stringstream& message,
              const _aflags& aflags,
              std::ofstream& FileMESSAGE,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              const char& type = _LOGGER_MESSAGE_,
              std::ostream& oss = std::cout,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              std::ofstream& FileMESSAGE,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              const std::string& directory,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              const _aflags& aflags,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              const std::string& directory,
              std::ofstream& FileMESSAGE,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // main std::function
  void logger(const std::string& filename,
              const std::string& function_name,
              const std::string& _message,
              const _aflags& aflags,
              std::ofstream& FileMESSAGE,
              std::ostream& oss = std::cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const std::string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // main function
  // END - added by Corey Oses - May 2017
  xstructure LTCELL(const std::string& options, std::istream& input);
  void MagneticParameters(std::string _directory, std::ostream& oss = std::cout);
  xstructure MINKOWSKIBASISREDUCTION(std::istream& input);
  std::string MISCIBILITY(std::vector<std::string> argv);
  void MOM(std::istream& input);
  void MSI(std::istream& input);
  uint NATOMS(std::istream& input);
  //[CO20190520 - not in pflow]xmatrix<double> GetDistMatrix(const xstructure& a); //CO20171025
  std::string NBONDXX(std::istream& input, bool aflowlib_legacy_format = false); // CO20171025
  std::string NBONDXX(const xstructure& a, bool aflowlib_legacy_format = false); // CO20171025
  xstructure NAMES(std::vector<std::string>, std::istream& input);
  xstructure NANOPARTICLE(std::istream& input, const aurostd::xvector<double>& iparams);
  xstructure NIGGLI(std::istream& input);
  void NDATA(std::istream& input);
  double NNDIST(std::istream& input);
  xstructure NOORDERPARAMETER(std::istream& input);
  xstructure NOSD(std::istream& input);
  xstructure NUMNAMES(std::vector<std::string>, std::istream& input);
  uint NSPECIES(std::istream& input);
  bool OPARAMETER(std::vector<std::string>, std::istream& input);
  bool SHIRLEY(std::vector<std::string>, std::istream& input);
  bool SHIRLEY2(std::vector<std::string>, std::istream& input);
  std::string PEARSON_SYMBOL(std::istream& input);
  std::string PEARSON_SYMBOL(std::istream& input, const aurostd::xoption& vpflow); // DX20210611 - added vpflow
  bool POCCUPATION(std::vector<std::string>, std::istream& input);
  void PDB(std::istream& input);
  void PDOS(std::vector<std::string>);
  void PHONONS(_aflags& aflags, std::istream& input, const double& radius);
  void PGROUP(_aflags& aflags, std::istream& input);
  void PGROUPXTAL(_aflags& aflags, std::istream& input);
  void PGROUPK(_aflags& aflags, std::istream& input);
  void PLANEDENS(std::vector<std::string>);
  std::string PLATON(const std::string& options, std::istream& input);
  std::string SG(aurostd::xoption& vpflow, std::istream& input, std::string mode, std::string print);
  void STATDIEL(std::vector<std::string>& argv); // CAMILO
  bool SYMMETRY_GROUPS(_aflags& aflags, std::istream& input, aurostd::xoption& vpflow, std::ostream& oss = std::cout); // DX20170818 - Add no_scan option to all symmetry Xgroups
  void POCC(std::vector<std::string>);
  std::string POSCAR2AFLOWIN(std::istream& input, const std::string& = ""); // Modified ME20181113
  void POSCAR2WYCKOFF(std::istream& input);
  std::vector<std::string> GENERATE_CERAMICS(const std::vector<std::string>& nonmetals, const std::vector<std::string>& metals, uint metal_arity); // CO20200731
  std::vector<std::string> GENERATE_CERAMICS(const aurostd::xoption& vpflow); // CO20200731
  std::string GENERATE_CERAMICS_PRINT(const aurostd::xoption& vpflow); // CO20200731
  bool PSEUDOPOTENTIALS_CHECK(const aurostd::xoption& vpflow, const std::string& file, std::ostream& oss = std::cout); // SC20200330
  void PYTHON_MODULES(const std::string& modules, std::ostream& oss = std::cout); // ME20211103
  void PYTHON_MODULES(const std::string& modules, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // ME20211103
  void PYTHON_MODULES(const std::vector<std::string>& vmodules, std::ostream& oss = std::cout); // ME20211103
  void PYTHON_MODULES(const std::vector<std::string>& vmodules, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // ME20211103
  bool QMVASP(aurostd::xoption& vpflow); // std::vector<std::string> argv); //CO20180703
  xstructure POSCAR(std::istream& input);
  void RENDER(std::istream& input, const std::filesystem::path& output);
  aurostd::xmatrix<double> QE_ibrav2lattice(const int& ibrav, const aurostd::xvector<double>& parameters, const bool& isabc); // DX20180123 - added more robust QE reader
} // namespace pflow
bool RequestedAlphabeticLabeling(std::string& label);
bool AlphabetizePrototypeLabelSpecies(std::deque<std::string>& species, std::deque<std::string>& species_pp, std::deque<double>& volumes, std::deque<double>& masses, std::string& label);
bool AlphabetizePrototypeLabelSpecies(std::deque<std::string>& species, std::deque<std::string>& species_pp, std::string& label);
bool AlphabetizePrototypeLabelSpecies(std::deque<std::string>& species, std::deque<double>& volumes, std::string& label);
bool AlphabetizePrototypeLabelSpecies(std::deque<std::string>& species, std::string& label);
std::string AlphabetizePrototypeLabelSpeciesArgv(std::vector<std::string>& argv);
namespace pflow {
  bool PROTO_PARSE_INPUT(const std::vector<std::string>& params, std::vector<std::vector<std::string>>& vstr, std::vector<std::vector<double>>& vnum, bool ignore_label = false, bool reverse = false); // CO20181226
  bool PROTO_TEST_INPUT(const std::vector<std::vector<std::string>>& vvstr, const std::vector<std::vector<double>>& vvnum, uint& nspeciesHTQC, bool patch_nspecies = false); // CO20181226
  bool sortPOCCSites(const std::string& p1, const std::string& p2); // CO20181226
  bool sortPOCCOccs(const std::string& occ1, const std::string& occ2); // CO20181226
  bool FIX_PRECISION_POCC(const std::string& occ, std::string& new_occ); // CO20181226
  void FIX_POCC_PARAMS(const xstructure& xstr, std::string& pocc_params); // CO20181226
  bool checkAnionSublattice(const xstructure& xstr); // CO20210201
  bool convertXStr2POCC(xstructure& xstr, const std::string& pocc_params, const std::vector<std::string>& vspecies, const std::vector<double>& vvolumes); // CO20181226
  bool POccInputs2Xstr(const std::string& pocc_input, aurostd::xoption& pocc_settings, xstructure& xstr, std::ostream& oss); // CO20211130
  bool POccInputs2Xstr(const std::string& pocc_input, aurostd::xoption& pocc_settings, xstructure& xstr, std::ofstream& FileMESSAGE, std::ostream& oss); // CO20211130
  xstructure PROTO_LIBRARIES(aurostd::xoption vpflow);
  bool PROTO_AFLOW(aurostd::xoption vpflow, bool flag_REVERSE); // too many options
  bool PROTO_ICSD_AFLOWIN(std::vector<std::string>& argv);
  xstructure PRIM(std::istream& input, uint mode);
  void RASMOL(const std::string& options, std::istream& input);
  void RBANAL(std::vector<std::string>);
  void RBDIST(std::vector<std::string>);
  xstructure RMATOM(std::istream& input, const int& iatom);
  xstructure RMCOPIES(std::istream& input);
  void RAYTRACE(std::vector<std::string>);
  xstructure SCALE(const std::string& options, std::istream& input);
  void RDF(const std::string& options, std::istream& input, bool raw_counts = false); // CO20220627
  void RDFCMP(const std::string& options);
  void RSM(std::vector<std::string>, std::istream& input);
  xstructure SD(std::vector<std::string>, std::istream& input);
  xstructure SETCM(std::istream& input, const aurostd::xvector<double>& cm);
  xstructure SETORIGIN(std::istream& input, const aurostd::xvector<double>& origin);
  xstructure SETORIGIN(std::istream& input, const int& natom);
  void SEWALD(std::vector<std::string>, std::istream& input);
  void SG(std::istream& input);
  void SGROUP(_aflags& aflags, std::istream& input, double radius);
  void SHELL(const std::string& options, std::istream& input);
  std::string SPECIES(std::istream& input);
  xstructure SHIFT(const std::string& options, std::istream& input);
  void SPLINE(std::vector<std::string>);
  void SUMPDOS(std::vector<std::string>);
  xstructure SUPERCELL(const std::string& options, std::istream& input);
  void SUPERCELLSTRLIST(const std::string& options);
  xstructure xstrSWAP(std::vector<std::string>, std::istream& input);
  xstructure VOLUME(const std::string& options, std::istream& input);
  std::string WyckoffPositions(aurostd::xoption& vpflow, std::istream& input); // DX20180807 - added wyccar to pflow //DX20210525 - changed name and generalized std::function
  xstructure WYCKOFF(std::vector<std::string>, std::istream& input);
  void XRAY(const std::string& options, std::istream& input);
  void XRAY_PEAKS(const aurostd::xoption& vpflow, std::istream& input); // CO20190409
  void READ_XRAY_DATA(const std::string& filename, std::vector<double>& v_twotheta, std::vector<double>& intensity); // CO20190620
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, std::istream& input); // CO20190409
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, const xstructure& str); // CO20190409
  void PRINT_XRAY_DATA_PLOT(std::istream& input, double lambda = XRAY_RADIATION_COPPER_Kalpha, const std::string& directory = ""); // CO20190409
  void PRINT_XRAY_DATA_PLOT(const xstructure& str, double lambda = XRAY_RADIATION_COPPER_Kalpha, const std::string& directory = ""); // CO20190409
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, const std::string& directory = ""); // CO20190620
  void PRINT_XRAY_DATA_PLOT(const std::string& filename, const std::string& directory = ""); // CO20190620
  void PRINT_XRAY_DATA_PLOT(const std::vector<double>& v_twotheta, const std::vector<double>& v_intensity, const std::string& directory = ""); // CO20190620
  void PLOT_XRAY(const aurostd::xoption& vpflow, std::istream& input); // CO20190409
  void PLOT_XRAY(const aurostd::xoption& vpflow, const xstructure& str); // CO20190409
  void PLOT_XRAY(std::istream& input, double lambda = XRAY_RADIATION_COPPER_Kalpha, const std::string& directory = "", bool keep_gp = false, bool force_generic_title = false); // CO20190409
  void PLOT_XRAY(const xstructure& str, double lambda = XRAY_RADIATION_COPPER_Kalpha, const std::string& directory = "", bool keep_gp = false, bool force_generic_title = false); // CO20190409
  void PLOT_XRAY(const aurostd::xoption& vpflow, const std::string& title = "", const std::string& directory = "", bool keep_gp = false); // CO20190620
  void PLOT_XRAY(const std::string& filename, const std::string& title = "", const std::string& directory = "", bool keep_gp = false); // CO20190620
  void PLOT_XRAY(const std::vector<double>& v_twotheta, const std::vector<double>& v_intensity, const std::string& title = "", const std::string& directory = "", bool keep_gp = false); // CO20190620
  void XYZ(const std::string& options, std::istream& input);
  void XYZINSPHERE(std::istream& input, double radius);
  void XYZWS(std::istream& input);
  void XelementPrint(const std::string& options, std::ostream& oss = std::cout);
  void ZVAL(const std::string& options);
} // namespace pflow

// aflow_pflow_print.cpp
namespace pflow {
  void PrintACE(const xstructure&, std::ostream& oss = std::cout);
  void PrintAngles(xstructure str, const double& cutoff, std::ostream& oss = std::cout);
  class projdata;
  void PrintBands(const pflow::projdata& pd);
  bool PrintCHGCAR(const xstructure& str,
                   const std::stringstream& chgcar_header,
                   const std::vector<int>& ngrid,
                   const std::vector<int>& format_dim,
                   const std::vector<double>& chg_tot,
                   const std::vector<double>& chg_diff,
                   const std::string& output_name,
                   std::ostream& oss = std::cout);
  void PrintChgInt(std::vector<aurostd::matrix<double>>& rad_chg_int, aurostd::matrix<double>& vor_chg_int, std::ostream& oss = std::cout); // CO20200404 pflow::matrix()->aurostd::matrix()
  void PrintCIF(std::ostream& oss, const xstructure&, int = 1, int = 1); // DX20180806 - added setting default
  void PrintClat(const aurostd::xvector<double>& data, std::ostream& oss = std::cout);
  void PrintCmpStr(const xstructure& str1, const xstructure& str2, const double& rcut, std::ostream& oss = std::cout);
  std::string PrintData(const xstructure& xstr, const std::string& smode = "DATA", filetype ftype = txt_ft, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1); // DX20210301
  std::string PrintData(const xstructure& xstr,
                        aurostd::xoption& vpflow,
                        const std::string& smode = "DATA",
                        filetype ftype = txt_ft,
                        bool already_calculated = false,
                        double sym_eps = AUROSTD_MAX_DOUBLE,
                        bool no_scan = false,
                        int setting = 1); // DX20210301
  std::string PrintData(const xstructure& xstr,
                        xstructure& str_sp,
                        xstructure& str_sc,
                        aurostd::xoption& vpflow,
                        const std::string& smode = "DATA",
                        filetype ftype = txt_ft,
                        bool already_calculated = false,
                        double sym_eps = AUROSTD_MAX_DOUBLE,
                        bool no_scan = false,
                        int setting = 1); // DX20210301
  std::string PrintData(const xstructure& xstr,
                        xstructure& str_sym,
                        xstructure& str_sp,
                        xstructure& str_sc,
                        const std::string& smode = "DATA",
                        filetype ftype = txt_ft,
                        bool already_calculated = false,
                        double sym_eps = AUROSTD_MAX_DOUBLE,
                        bool no_scan = false,
                        int setting = 1); // DX20210301
  std::string PrintData(const xstructure& xstr,
                        xstructure& str_sym,
                        xstructure& str_sp,
                        xstructure& str_sc,
                        aurostd::xoption& vpflow,
                        const std::string& smode = "DATA",
                        filetype ftype = txt_ft,
                        bool already_calculated = false,
                        double sym_eps = AUROSTD_MAX_DOUBLE,
                        bool no_scan = false,
                        int setting = 1); // DX20210301
  std::string PrintRealLatticeData(const xstructure& xstr, const std::string& smode = "DATA", filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  std::string PrintRealLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, const std::string& smode = "DATA", filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  std::string PrintRealLatticeData(const xstructure& xstr,
                                   aurostd::xoption& vpflow,
                                   aurostd::JSON::object& json,
                                   const std::string& smode = "DATA",
                                   filetype ftype = txt_ft,
                                   bool standalone = true,
                                   bool already_calculated = false,
                                   double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  std::string PrintLatticeLatticeData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  std::string PrintLatticeLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  std::string PrintLatticeLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  std::string PrintCrystalPointGroupData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  std::string PrintCrystalPointGroupData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  std::string PrintCrystalPointGroupData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  std::string PrintReciprocalLatticeData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  std::string PrintReciprocalLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  std::string PrintReciprocalLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  std::string PrintSuperlatticeData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  std::string PrintSuperlatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  std::string PrintSuperlatticeData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209

  void PrintDisplacements(xstructure str, const double cutoff, std::ostream& oss = std::cout);
  void PrintDistances(xstructure str, const double cutoff, std::ostream& oss = std::cout);
  void PrintEwald(const xstructure& in_str, double& epoint, double& ereal, double& erecip, double& eewald, double& eta, const double& SUMTOL, std::ostream& oss = std::cout);
  void PrintGulp(const xstructure&, std::ostream& oss = std::cout);
  std::string PrintSGData(xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1, bool supress_Wyckoff = false); // DX20210211
  std::string PrintSGData(xstructure& xstr,
                          aurostd::xoption& vpflow,
                          filetype ftype = txt_ft,
                          bool standalone = true,
                          bool already_calculated = false,
                          double sym_eps = AUROSTD_MAX_DOUBLE,
                          bool no_scan = false,
                          int setting = 1,
                          bool suppress_Wyckoff = false); // DX20210211
  std::string PrintSGData(xstructure& xstr,
                          aurostd::xoption& vpflow,
                          aurostd::JSON::object& json,
                          filetype ftype = txt_ft,
                          bool standalone = true,
                          bool already_calculated = false,
                          double sym_eps = AUROSTD_MAX_DOUBLE,
                          bool no_scan = false,
                          int setting = 1,
                          bool suppress_Wyckoff = false); // DX20210211

  std::string PrintWyckoffData(xstructure& xstr, filetype ftype = txt_ft, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1); // DX20210610
} // namespace pflow
void PrintKmesh(const aurostd::xmatrix<double>& kmesh, std::ostream& oss = std::cout); // HERE
void PrintImages(xstructure strA, xstructure strB, const int& ni, const std::string& path_flag);
void PrintMSI(const xstructure&, std::ostream& oss = std::cout);
void PrintNdata(const xstructure&, std::ostream& oss = std::cout);
// void PrintNeatProj(projdata& pd);
void PrintPDB(const xstructure&, std::ostream& oss = std::cout);
void platon2print(xstructure, bool P_EQUAL, bool P_EXACT, double P_ang, double P_d1, double P_d2, double P_d3, std::ostream& sout);
std::string RDF2string(const xstructure& xstr, const double rmax, const int nbins, const aurostd::xmatrix<double>& rdf_all); // CO20220627
void PrintRDF(const xstructure& xstr, const double rmax, const int nbins, const aurostd::xmatrix<double>& rdf_all, std::ostream& oss = std::cout); // CO20220627
void PrintRDFCmp(const xstructure& str_A,
                 const xstructure& str_B,
                 const double& rmax,
                 const int nbins,
                 const double& smooth_width,
                 const int nsh,
                 const aurostd::matrix<double>& rdfsh_all_A,
                 const aurostd::matrix<double>& rdfsh_all_B,
                 const std::vector<int>& best_match,
                 const aurostd::matrix<double>& rms_mat,
                 std::ostream& oss = std::cout); // CO20200404 pflow::matrix()->aurostd::matrix()
void PrintRSM(const xstructure&, std::ostream& oss = std::cout);
void PrintShell(const xstructure& str, const int& ns, const double& rmin, const double& rmax, const std::string& sname, const int lin_dens, std::ostream& oss = std::cout);
double CorrectionFactor(const double& th);
void PrintXray(const xstructure& str, double l, std::ostream& oss = std::cout); // CO20190520
void PrintXYZ(const xstructure& a, const aurostd::xvector<int>& n, std::ostream& oss = std::cout);
void PrintXYZws(const xstructure& a, std::ostream& oss = std::cout);
void PrintXYZInSphere(const xstructure& a, const double& radius, std::ostream& oss = std::cout);

// aflow_pflow_funcs.cpp
double DebyeWallerFactor(double theta, double temp, double debye_temp, double mass, double lambda = XRAY_RADIATION_COPPER_Kalpha);
std::string getGenericTitleXStructure(const xstructure& xstr, bool latex = false); // CO20190520
aurostd::xvector<double> balanceChemicalEquation(const std::vector<aurostd::xvector<double>>& _lhs, const std::vector<aurostd::xvector<double>>& _rhs, bool normalize, double tol); // CO20180817
aurostd::xvector<double> balanceChemicalEquation(const aurostd::xmatrix<double>& _composition_matrix, bool normalize, double tol);
void ParseChemFormula(std::string& ChemFormula, std::vector<std::string>& ChemName, std::vector<float>& ChemConc);
void ParseChemFormulaIndividual(uint nchar, std::string& ChemFormula, std::string& AtomSymbol, float& AtomConc);

namespace pflow { // CO20190601
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout); // CO20190321

  bool findClosedPackingPlane(std::istream& input); // CO20191110
  bool findClosedPackingPlane(const xstructure& xstr); // CO20191110
} // namespace pflow

namespace pflow {
  void GetXray2ThetaIntensity(const xstructure& str, std::vector<double>& v_twotheta, std::vector<double>& v_intensity, double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190520
  std::vector<uint> GetXrayPeaks(const xstructure& str, std::vector<double>& v_twotheta, std::vector<double>& v_intensity, std::vector<double>& v_intensity_smooth, double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190520 //CO20190620 - v_peaks_amplitude not needed
  std::vector<uint> GetXrayPeaks(const std::vector<double>& v_twotheta, const std::vector<double>& v_intensity, std::vector<double>& v_intensity_smooth); // CO20190520  //CO20190620 - v_peaks_amplitude not needed
  void GetXray(const xstructure& str, std::vector<double>& dist, std::vector<double>& sf, std::vector<double>& scatt_fact, std::vector<double>& mass, std::vector<double>& twoB_vec, double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190520
  void GetXrayData(const xstructure& str,
                   std::vector<double>& dist,
                   std::vector<double>& sf,
                   std::vector<double>& scatt_fact,
                   std::vector<double>& mass,
                   std::vector<double>& twoB_vec,
                   std::vector<std::vector<double>>& ids,
                   aurostd::matrix<double>& data,
                   double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190409  //CO20190620 - intmax can be grabbed later  //CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<std::vector<int>> getvvitypesRDF(const xstructure& xstr); // CO20220627
  void GetRDF(const xstructure& xstr, aurostd::xmatrix<double>& rdf_all, const double rmax = 5.0, const int nbins = 25, bool raw_counts = false, const double sigma = 0.0, const int window_gaussian = 0); // CO20220627 - rewritten
  void GetRDFShells(const xstructure& str, const double& rmax, const int& nbins, const int& smooth_width, const aurostd::matrix<double>& rdf, aurostd::matrix<double>& rdfsh, aurostd::matrix<double>& rdfsh_loc); // CO20200404 pflow::matrix()->aurostd::matrix()
  double RdfSh_RMS(const int iaA, const int iaB, const int nsh_max, const int nt, const aurostd::matrix<double>& rdfsh_all_A, const aurostd::matrix<double>& rdfsh_all_B); // CO20200404 pflow::matrix()->aurostd::matrix()
  void CmpRDFShells(const xstructure& str_A,
                    const xstructure& str_B,
                    const aurostd::matrix<double>& rdfsh_all_A,
                    const aurostd::matrix<double>& rdfsh_all_B,
                    const int nsh,
                    std::vector<int>& best_match,
                    aurostd::matrix<double>& rms_mat); // CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> GetSmoothRDF(const aurostd::matrix<double>& rdf, const double& sigma); // CO20200404 pflow::matrix()->aurostd::matrix()
  void CmpStrDist(xstructure str1, xstructure str2, const double& cutoff, aurostd::matrix<double>& dist1, aurostd::matrix<double>& dist2, aurostd::matrix<double>& dist_diff, aurostd::matrix<double>& dist_diff_n); // CO20200404 pflow::matrix()->aurostd::matrix()
} // namespace pflow

// aflow_pflow.cpp
int SignNoZero(const double& x);
int Nint(const double& x);
int Sign(const double& x);

// ---------------------------------------------------------------------------
// PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDOSDATA PDO
namespace pflow {
  class pdosdata {
  public:
    // Constructors
    pdosdata(); // default
    void PrintParams(std::ostream& oss, const std::vector<std::string>& Ltotnames);
    void PrintPDOS(std::ostream& oss, const int& sp);
    // variables
    std::string PDOSinfile;
    aurostd::matrix<int> pdos_at; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<int> pdos_k; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<int> pdos_b; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<int> pdos_lm; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> pdos; // CO20200404 pflow::matrix()->aurostd::matrix()
    double emin;
    double emax;
    double efermi;
    double smooth_sigma;
    int spin;
    int nlm;
    int natoms;
    int print_params;
    int nbins;
  };
} // namespace pflow

// ---------------------------------------------------------------------------
// RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY TRACING RAY
namespace pflow {
  class rtparams {
  public:
    void free();
    void copy(const rtparams& b);
    // Constructors
    rtparams(); // default
    rtparams(const rtparams& b); // default
    // Operators
    const rtparams& operator=(const rtparams& b);
    // Accessors
    void Print() const;
    // variables
    double resx;
    double resy;
    std::vector<double> viewdir;
    int viewdir_s;
    std::vector<double> updir;
    int updir_s;
    double zoom;
    double aspectratio;
    double antialiasing;
    double raydepth;
    std::vector<double> center;
    std::vector<double> center_guide;
    int center_s;
    std::vector<double> background;
    aurostd::matrix<double> lightcenter; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<double> lightrad;
    aurostd::matrix<double> lightcolor; // CO20200404 pflow::matrix()->aurostd::matrix()
    // Sphere texture variables (ambient, diffuse, specular, opacity)
    aurostd::matrix<double> sphtex_tex; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<double> sphtex_tex_def;
    aurostd::matrix<double> sphtex_color; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<double> sphtex_color_def;
    aurostd::matrix<double> sphtex_phong; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<double> sphtex_phong_def;
    std::vector<std::string> sphtex_names;
    std::vector<double> sph_rad;
    // Plane variables
    int plane;
    int plane_s;
    std::vector<double> plane_center;
    std::vector<double> plane_normal;
    std::vector<double> plane_color;
    std::vector<double> planetex_tex;
    std::vector<double> plane_center_def;
    std::vector<double> plane_normal_def;
    std::vector<double> plane_color_def;
    std::vector<double> planetex_tex_def;
    int plane_center_s;
    int plane_normal_s;
    int plane_color_s;
    int planetex_tex_s;
    double sph_rad_def;
    std::string shading;
    std::string outfile;
    aurostd::matrix<double> sc; // CO20200404 pflow::matrix()->aurostd::matrix()
    int sc_s;
    int calc_type;
    std::vector<std::string> input_files;
    int first_set;
    std::string insert_file;
    std::vector<double> rotation;
    std::vector<double> struct_origin;
    int struct_origin_s;
  };

  void SetRTParams(xstructure& str, pflow::rtparams& rtp);
  std::vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(std::ifstream& outcar_inf, std::ifstream& xdatcar_inf);
  void GetDatFromOutcar(std::vector<aurostd::matrix<double>>& lat_vec, std::deque<int>& num_each_type, std::ifstream& outcar_inf); // CO20200404 pflow::matrix()->aurostd::matrix()
  void GetDatFromXdatcar(std::vector<aurostd::matrix<double>>& fpos_vec, std::ifstream& xdatcar_inf); // CO20200404 pflow::matrix()->aurostd::matrix()
  std::vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(std::ifstream& outcar_inf, std::ifstream& xdatcar_inf);
  void PrintStrVec(const std::vector<xstructure>& str_vec, std::ostream& outf);
  void ReadInRTParams(std::ifstream& rtinfile, pflow::rtparams& rtp);
  void ReadInStrVec(std::vector<xstructure>& str_vec, std::ifstream& strlist_inf);
  void JoinStrVec(std::vector<xstructure> str_vec_1, std::vector<xstructure> str_vec_2, std::vector<xstructure>& str_vec_tot);
  void SetStrFromRTParams(xstructure& str, pflow::rtparams& rtp);
  void SuperCellStrVec(std::vector<xstructure>& str_vec, const aurostd::matrix<double>& sc); // CO20200404 pflow::matrix()->aurostd::matrix()
  void UpDateRTParams(pflow::rtparams& rtp, const int& istr, int nstr);
  void SetRTParams(xstructure& str, pflow::rtparams& rtp);
  void GetRTDatFile(xstructure str, const pflow::rtparams& rtp, std::ostringstream& rtdat_file);
  std::string PrintRTDatFile(std::ostringstream& rtdat_file, const pflow::rtparams& rt_params);
  std::string CreateRTtgaFile(const std::string& datfile, const pflow::rtparams& rt_params);
  std::string CreateRTjpgFile(const std::string& tgafile, const pflow::rtparams& rt_params);
  void GetRTencFile(const pflow::rtparams& rtp, const int nim, std::ostringstream& os);
  std::string PrintRTencFile(const pflow::rtparams& rt_params, std::ostringstream& rtenc_file);
  std::string CreateRTmpgFile(const pflow::rtparams& rt_params, const std::string& encfile);
  void RayTraceManager(std::vector<std::string>);

} // namespace pflow

// ---------------------------------------------------------------------------
// PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PRO

namespace pflow {
  class projdata {
  public:
    // Constructors
    projdata(); // default
    void Print(std::ostream& outf);
    // variables
    int nl_max; // 4 for s,p,d,f orbitals
    int nlm_max; // 16 for s,p,d,f orbitals
    int nlmtot_max; // 20 for s,p,d,f orbitals + p,d,f,all totals
    int nl; // 3 for spd
    int nlm; // 9 for spd
    int nlmtot; // 9+Psum+Dsum+Allsum=12 (for spd)
    int nkpts;
    int nbands;
    int nions;
    int ntypes;
    std::vector<int> num_each_type;
    aurostd::matrix<double> wfermi_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> wfermi_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<double> wkpt;
    std::vector<aurostd::matrix<aurostd::matrix<std::complex<double>>>> pdat_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<aurostd::matrix<std::complex<double>>>> pdat_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lmtot_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lmtot_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_kpt_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_kpt_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_kpt_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_kpt_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_bnd_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_bnd_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_bnd_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<aurostd::matrix<double>> occ_vs_ion_bnd_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> ener_k_b_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> ener_k_b_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    int sp;
    int rspin;
    aurostd::matrix<double> kpts; // CO20200404 pflow::matrix()->aurostd::matrix()
    std::vector<std::string> LMnames;
    std::vector<std::string> Lnames;
    std::vector<std::string> LLMnames;
    std::string PROOUTinfile;
    aurostd::matrix<double> lat; // CO20200404 pflow::matrix()->aurostd::matrix()
  };
} // namespace pflow

// ---------------------------------------------------------------------------
// PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJFUNCS PROJF

// in aflow_pflow_funcs.cpp
namespace pflow {
  std::complex<double> ProcessProjection(const std::complex<double>& proj);
  void ReadInProj(projdata& pd);
  void CalcNeatProj(projdata& pd, int only_occ);
  void ReadInPDOSData(const projdata& prd, pdosdata& pdd);
  void CalcPDOS(const projdata& prd, pdosdata& pdd);
  void SmoothPDOS(const projdata& prd, pdosdata& pdd);
  void AtomCntError(const std::string& tok, const int tokcnt, const int atom_cnt);
} // namespace pflow

// in aflow_pflow_print.cpp
void PrintNeatProj(pflow::projdata& pd, std::ostream& outf);

// ---------------------------------------------------------------------------
// SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUNCS SUMPDOSFUN

// in aflow_pflow_funcs.cpp
namespace pflow {
  void ReadSumDOSParams(std::ifstream& infile, pflow::pdosdata& pdd);
  void ReadInPDOSData(aurostd::matrix<aurostd::matrix<double>>& allpdos, pflow::pdosdata& pdd, std::ifstream& infile); // CO20200404 pflow::matrix()->aurostd::matrix()
  void SumPDOS(const aurostd::matrix<aurostd::matrix<double>>& allpdos, pflow::pdosdata& pdd); // CO20200404 pflow::matrix()->aurostd::matrix()
} // namespace pflow

// in pflow_print
void PrintSumPDOS(pflow::pdosdata& pdd, std::ostream& outf);

// ---------------------------------------------------------------------------
// RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBFUNCS RBF

// in aflow_pflow_funcs.cpp
namespace pflow {
  double TotalAtomDist(xstructure str, xstructure str00, const std::string& path_flag);
  std::vector<std::string> GetRBDir(const int& nim);
  std::vector<double> GetRBEner(const int& nim);
  std::vector<xstructure> GetRBStruct(const int& nim);
  std::vector<double> GetRBDistCum(const std::vector<xstructure>& str_vec, const std::string& path_flag);
  std::vector<double> GetRBDistFromStrI(const std::vector<xstructure>& str_vec, const xstructure& strI, const std::string& path_flag);
  void RBPoscarDisp(const xstructure& str1in, const xstructure& str2in, xstructure& diffstr, double& totdist, aurostd::matrix<double>& cm, const std::string& path_flag); // CO20200404 pflow::matrix()->aurostd::matrix()
} // namespace pflow

// in aflow_pflow_print.cpp
void PrintRBAnal(const int& nim, const std::string& path_flag, std::ostream& outf);
void PrintRBPoscarDisp(const xstructure& diffstr, double& totdist, aurostd::matrix<double>& cm, const std::string& path_flag, std::ostream& outf); // CO20200404 pflow::matrix()->aurostd::matrix()

// ---------------------------------------------------------------------------
// CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUNCS CHARGE FUN

namespace pflow {
  class pd_params {
  public:
    std::string type;
    double scale;
    aurostd::matrix<double> pts; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> dpts; // CO20200404 pflow::matrix()->aurostd::matrix()
    int Nx, Ny;
    std::string orig_loc;
    std::string ortho;
    void Print(std::ostream& outf) const;
  };

  bool ReadCHGCAR(xstructure& str,
                  std::stringstream& chgcar_header,
                  std::vector<int>& ngrid,
                  std::vector<int>& format_dim,
                  std::vector<double>& chg_tot,
                  std::vector<double>& chg_diff,
                  std::stringstream& chgcar_ss,
                  std::ostream& oss = std::cout);
  bool ReadChg(xstructure& str, std::vector<int>& ngrid, std::vector<double>& chg_tot, std::vector<double>& chg_diff, std::istream& chgfile);
  void GetChgInt(std::vector<aurostd::matrix<double>>& rad_chg_int,
                 aurostd::matrix<double>& vor_chg_int, // CO20200404 pflow::matrix()->aurostd::matrix()
                 xstructure& str,
                 std::vector<int>& ngrid,
                 std::vector<double>& chg_tot,
                 std::vector<double>& chg_diff);
  void ReadPlaneDensParams(const xstructure& str, pd_params& pdp, std::istream& infile);
  void GetPlaneDens(const pd_params& pdp, std::vector<double>& dens2d_tot, std::vector<double>& dens2d_diff, const xstructure& str, const std::vector<int>& ngrid, const std::vector<double>& chg_tot, const std::vector<double>& chg_diff);
  void PrintPlaneDens(const pd_params& pdp, const std::vector<double>& dens2d_tot, const std::vector<double>& dens2d_diff, const xstructure& str);
} // namespace pflow

// ---------------------------------------------------------------------------
// EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWALD FUNCS EWA
namespace pflow {
  void Ewald(const xstructure& in_str, double& epoint, double& ereal, double& erecip, double& eewald, double& eta, const double& SUMTOL);
  double GetEta(const int& natoms, const double& vol);
  double GetPointEner(const double& rteta, const std::vector<double>& atchg, const double& vol);
  double GetRecipEner(const double& eta, const std::vector<double>& atchg, const double& vol, const aurostd::matrix<double>& rlat, const aurostd::matrix<double>& fpos, const double& SUMTOL); // CO20200404 pflow::matrix()->aurostd::matrix()
  double GetRealEner(const double& eta, const std::vector<double>& atchg, const double& vol, const aurostd::matrix<double>& lat, const aurostd::matrix<double>& fpos, const double& SUMTOL); // CO20200404 pflow::matrix()->aurostd::matrix()
  double GetScreenedESEner();
  double ScreenedESEner(const xstructure& in_str, const double& Ks, const double& SUMTOL);
} // namespace pflow

// Output help information of an option
void helpIndividualOption(std::vector<std::string>& argv);

// ---------------------------------------------------------------------------
// FORMER WAHYU.H

void AConvaspBandgap(std::vector<std::string>& bandsdir);
void AConvaspBandgaps(std::istream& bandsdir, std::ostream& oss = std::cout);
void AConvaspBandgaps(std::istream& bandsdir, std::ostringstream& oss);
void AConvaspBandgapFromDOS(std::istream& doscar);
void AConvaspBandgapListFromDOS(std::istream& doscar);
namespace pflow {
  void ICSD(std::vector<std::string> argv, std::istream& input);
  void ICSD_2POSCAR(std::istream& input);
  void ICSD_2PROTO(std::istream& input);
  void ICSD_2WYCK(std::istream& input, bool SOF);
  void ICSD_ListMetals();
} // namespace pflow
float GetBandGap_WAHYU(std::stringstream& straus, float Efermi, char& gaptype);
float GetBandgapFromDOS(std::ifstream& doscar);
float GetBandgapFromDOS(std::istream& doscar);
std::vector<std::string> GetMetalsList(bool v);
bool IsMetal(const std::string s);
void ParseChemicalFormula(std::string Formula, std::vector<std::string>& Zname, std::vector<float>& Zconc);
std::string RemoveCharFromString(const std::string s, char c);
std::string RemoveStringFromTo(const std::string s, char cstart, char cstop);
int StringCrop(std::string s, std::vector<std::string>& vstring);
std::string StringCropAt(const std::string s, int icrop);
std::vector<float> SortFloat(std::vector<float> v, int mode);
aurostd::xvector<double> cross(const aurostd::xvector<double> a, const aurostd::xvector<double> b);

// ---------------------------------------------------------------------------
// FORMER RICHARD.H
namespace pflow {
  double GetAtomicPlaneDist(const std::string& options, std::istream& input);
  double frac2dbl(std::string str);
  bool havechar(std::string str_in, char c);
  int whereischar(std::string str, char c);
  void cleanupstring(std::string& str);
} // namespace pflow

// from KY's old files
namespace pflow {
  void BZMAX(std::istream& input);
}

// ME20190628
namespace pflow {
  // Precision for pretty printing
  const int COEF_PRECISION = 4;

  std::string prettyPrintCompound(const std::string&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20190629
  std::string prettyPrintCompound(const std::vector<std::string>&, const std::vector<uint>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //DX20200727 - uint std::variant
  std::string prettyPrintCompound(const std::vector<std::string>&, const std::vector<double>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20190629
  std::string prettyPrintCompound(const std::vector<std::string>&, const aurostd::xvector<uint>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20200727 - uint std::variant
  std::string prettyPrintCompound(const std::vector<std::string>&, const aurostd::xvector<double>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20190629

} // namespace pflow

//[CO20200526 - EASY TEMPLATE CLASS]namespace pflow {
//[CO20200526 - EASY TEMPLATE CLASS]  class AQueue : public xStream {
//[CO20200526 - EASY TEMPLATE CLASS]    public:
//[CO20200526 - EASY TEMPLATE CLASS]      //NECESSARY PUBLIC CLASS METHODS - START
//[CO20200526 - EASY TEMPLATE CLASS]      //constructors - START
//[CO20200526 - EASY TEMPLATE CLASS]      AQueue(ostream& oss=cout);
//[CO20200526 - EASY TEMPLATE CLASS]      AQueue(ofstream& FileMESSAGE,ostream& oss=cout);
//[CO20200526 - EASY TEMPLATE CLASS]      AQueue(const AQueue& b);
//[CO20200526 - EASY TEMPLATE CLASS]      //constructors - STOP
//[CO20200526 - EASY TEMPLATE CLASS]      ~AQueue();
//[CO20200526 - EASY TEMPLATE CLASS]      const AQueue& operator=(const AQueue& other);
//[CO20200526 - EASY TEMPLATE CLASS]      void clear();
//[CO20200526 - EASY TEMPLATE CLASS]      //NECESSARY PUBLIC CLASS METHODS - STOP
//[CO20200526 - EASY TEMPLATE CLASS]
//[CO20200526 - EASY TEMPLATE CLASS]      //general attributes
//[CO20200526 - EASY TEMPLATE CLASS]      bool m_initialized;
//[CO20200526 - EASY TEMPLATE CLASS]
//[CO20200526 - EASY TEMPLATE CLASS]      //initialization methods
//[CO20200526 - EASY TEMPLATE CLASS]      bool initialize(ostream& oss);
//[CO20200526 - EASY TEMPLATE CLASS]      bool initialize(ofstream& FilMESSAGE,ostream& oss);
//[CO20200526 - EASY TEMPLATE CLASS]    private:
//[CO20200526 - EASY TEMPLATE CLASS]      //NECESSARY private CLASS METHODS - START
//[CO20200526 - EASY TEMPLATE CLASS]      void free();
//[CO20200526 - EASY TEMPLATE CLASS]      void copy(const AQueue& b);
//[CO20200526 - EASY TEMPLATE CLASS]      //NECESSARY END CLASS METHODS - END
//[CO20200526 - EASY TEMPLATE CLASS]  };
//[CO20200526 - EASY TEMPLATE CLASS]}

enum job_status { // CO20200526
  JOB_RUNNING,
  JOB_QUEUED,
  JOB_HELD,
  JOB_DONE
};

enum node_status { // CO20200526
  NODE_FREE,
  NODE_OCCUPIED,
  NODE_FULL,
  NODE_DOWN,
  NODE_OFFLINE,
  NODE_OPERATIONAL, // NOT ASSIGNED - this is an aggregate of free+occupied+full
  NODE_NONOPERATIONAL, // NOT ASSIGNED - this is an aggregate of down+offline
};

enum cpus_status { // CO20200526
  CPUS_FREE,
  CPUS_OCCUPIED,
  CPUS_TOTAL,
};

enum queue_system { // CO20200526
  QUEUE_SLURM,
  QUEUE_TORQUE
};

namespace pflow {
  // AJob stays a struct until we need more than just free
  struct AJob { // CO20200526
    uint m_index; // reflection to m_jobs
    uint m_id;
    std::string m_user;
    job_status m_status;
    uint m_ncpus; // this is a "total" ncpus for the job (NOT an index)
    std::vector<uint> m_vinodes;
    std::vector<uint> m_vncpus; // this is ncpus split across nodes (NOT an index)
    std::vector<uint> m_vipartitions;
    void free();
  };
  // ANode stays a struct until we need more than just free
  struct ANode { // CO20200526
    uint m_index; // reflection to m_nodes
    std::string m_name;
    node_status m_status;
    uint m_ncpus;
    uint m_ncpus_occupied; // if we need to collect job information later, then this should become a getter based on job count
    std::string m_properties; // needed to match with queues
    std::vector<uint> m_vijobs;
    std::vector<uint> m_vipartitions;
    void free();
    [[nodiscard]] bool isStatus(const node_status& status) const;
  };
  // APartition stays a struct until we need more than just free
  struct APartition { // CO20200526
    uint m_index; // reflection to m_partitions
    std::string m_name;
    std::string m_properties_node; // needed to match with queues //also seems to be available ONLY to root user, so we hack for QRATS //http://docs.adaptivecomputing.com/torque/4-2-8/Content/topics/4-serverPolicies/mappingQueueToRes.htm
    std::vector<uint> m_inodes;
    std::vector<uint> m_vijobs;
    void free();
  };
} // namespace pflow

// CO20200526 - queueing class
namespace pflow {
  uint getTORQUEIDFromString(const std::string& torqueid_str);
  class AQueue : public xStream {
  public:
    AQueue(std::ostream& oss = std::cout) : xStream(oss) {}
    AQueue(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss) {}
    AQueue(const aurostd::xoption& vpflow, std::ostream& oss = std::cout) : xStream(oss), m_flags(vpflow) {}
    AQueue(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss), m_flags(vpflow) {}

    void clear() { *this = AQueue(); }

    // general attributes
    bool m_initialized = false;
    aurostd::xoption m_flags;
    queue_system m_qsys = QUEUE_SLURM;
    std::vector<APartition> m_partitions;
    std::vector<ANode> m_nodes;
    std::vector<AJob> m_jobs;

    // setters
    void setFlags(const aurostd::xoption& vpflow);

    // getters
    [[nodiscard]] uint getNNodes() const;
    [[nodiscard]] uint getNCPUS() const;
    [[nodiscard]] uint getNNodes(const APartition& partition) const;
    [[nodiscard]] uint getNCPUS(const APartition& partition) const;
    [[nodiscard]] uint getNNodes(const APartition& partition, const node_status& status) const;
    [[nodiscard]] uint getNCPUS(const APartition& partition, const node_status& status_node, const cpus_status& status_cpus = CPUS_TOTAL) const;
    [[nodiscard]] uint getNCPUS(const std::string& user, const std::string& partition, const job_status& status) const;
    [[nodiscard]] uint getNCPUS(const std::string& user, const APartition& partition, const job_status& status) const;
    [[nodiscard]] double getPercentage(const std::string& user, const std::string& partition, const job_status& status) const;
    [[nodiscard]] double getPercentage(const std::string& user, const APartition& partition, const job_status& status) const;
    [[nodiscard]] uint nodeName2Index(const std::string& name) const;
    [[nodiscard]] uint partitionName2Index(const std::string& name) const;

    // methods
    void getQueue(); // wrapper around processQueue() with try's for failed external calls
  private:
    void freeQueue();
    void processQueue(); // main processer for external queue commands (pbsnodes, qstat, squeue, etc.)

    void readNodesPartitionsSLURM();
    void readJobsSLURM();
    void readPartitionsTORQUE();
    void readNodesJobsTORQUE();
    void readJobsTORQUE();

    bool addJob(const AJob& _job);
    bool addPartition(const APartition& _partition);
    bool addNode(const ANode& _node);
    void nodePartitionMapping(ANode& node);
    void jobMapping(AJob& job);
  };
} // namespace pflow

// CO20200526 - queueing class
namespace pflow {
  std::string getQueueStatus(const aurostd::xoption& vpflow);
}

namespace pflow {
  std::vector<std::string> getFakeElements(uint nspecies); // DX20200728
  bool hasRealElements(const xstructure& xstr); // DX20210113
  double getSymmetryTolerance(const xstructure& xstr, const std::string& tolerance_string);
  std::vector<double> getSymmetryToleranceSpectrum(const std::string& tolerance_range_string);
  uint getSpaceGroupSetting(const std::string& setting_string, uint mode_default = 0); // DX20210420 - mode_default=0: unspecified, AFLOW will determine
} // namespace pflow

namespace pflow {
  void outputAtomicEnvironment(const std::string& auid, uint aeMode = 1, double radius = 4.0);
}

#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
