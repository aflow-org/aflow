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
  using std::cout;
  using std::deque;
  using std::ofstream;
  using std::ostream;
  using std::string;
  using std::stringstream;
  using std::vector;
  xstructure ABCCAR(std::istream& input);
  void ACE(std::istream& input);
  bool AddSpinToXstructure(xstructure& a, vector<double>& vmag); // DX20170927 - Magnetic symmetry
  bool AddSpinToXstructure(xstructure& a, vector<aurostd::xvector<double>>& vmag_noncoll); // DX20171205 - Magnetic symmetry (non-collinear)
  void AGROUP(_aflags& aflags, std::istream& input);
  void AGROUP2(std::istream& input);
  void AGROUP2m(std::istream& input);
  xstructure ALPHABETIC(std::istream& input);
  string ALPHACompound(const string& options);
  string ALPHASpecies(const string& options);
  string AFLOWIN(std::istream& input);
  void ANGLES(const string& options, std::istream& input);
  string ATOMSMAX(const string& options, std::istream& input);
  void BANDS(const string& options, std::istream& input);
  void BANDGAP(aurostd::xoption& vpflow, ostream& oss = cout); // CAMILO  //CO20171006
  void BANDGAP_DOS(aurostd::xoption& vpflow, ostream& oss = cout); // CAMILO  //CO20171006  //CO20191110
  void BANDSTRUCTURE(_aflags& aflags);
  string BZDirectionsLATTICE(const string& options);
  string BZDirectionsSTRUCTURE(std::istream& input, aurostd::xoption& vpflow); // DX20181102 - add options
  void CAGES(_aflags& aflags, const string& options, std::istream& input);
  // DX+CO START
  bool PerformFullSymmetry(xstructure& a);
  bool PerformFullSymmetry(xstructure& a, ofstream& FileMESSAGE, const string& directory, _kflags& kflags, bool osswrite, ostream& oss, string format = "txt");  // ME20200224
  bool PerformFullSymmetry(xstructure& a, ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, bool osswrite, ostream& oss, string format = "txt");
  bool PerformFullSymmetry(xstructure& a, double& tolerance, bool no_scan, bool force_perform, ofstream& FileMESSAGE, _aflags& aflags, _kflags& kflags, bool osswrite, ostream& oss, string format = "txt");
  void ProcessAndAddSpinToXstructure(xstructure& a, const string& magmom_info); // DX20190801
  void defaultKFlags4SymWrite(_kflags& kflags, bool write = true);
  void defaultKFlags4SymCalc(_kflags& kflags, bool calc = true);
  bool CalculateFullSymmetry(std::istream& input, aurostd::xoption& vpflow, ostream& oss = cout);
  bool CalculateFullSymmetry(_aflags& aflags, _kflags& kflags, xstructure& a, aurostd::xoption& vpflow, bool osswrite, ostream& oss = cout);
  bool fixEmptyAtomNames(xstructure& xstr, bool force_fix = false);  // force_fix=true if you want to override what is already in species
  // DX+CO END
  xstructure CART(std::istream& input);
  xstructure CORNERS(std::istream& input);
  void ChangeSuffix(const string& options);
  string CHGDIFF(aurostd::xoption vpflow);
  bool CHGDIFF(const string& chgcar1_file, const string& chgcar2_file, const string& output_File, ostream& oss = cout);
  // DX+CO START
  void CHGINT(vector<string>);
  // DX+CO END
  string CHGSUM(aurostd::xoption vpflow);
  bool CHGSUM(const string& chgcar_in1, const string& chgcar_in2, ostream& oss = cout);
  bool CHGSUM(string& species_header, const string& chgcar_in1, const string& chgcar_in2, const string& output_file, ostream& oss = cout);
  bool CHGSUM(const string& chgcar_in1, const string& chgcar_in2, const string& output_file, ostream& oss = cout);
  bool CHGSUM(const vector<string>& chgcar_files, ostream& oss = cout);
  bool CHGSUM(const vector<string>& chgcar_files, const string& output, ostream& oss = cout);
  bool CHGSUM(string& species_header, const vector<string>& chgcar_files, ostream& oss = cout);
  bool CHGSUM(string& species_header, const vector<string>& chgcar_files, const string& output_file, ostream& oss = cout);
  void CIF(std::istream& input, aurostd::xoption& vpflow);
  void getCIFOptions(xstructure& a, aurostd::xoption& vpflow); // DX20250319
  void CLAT(const string& options);
  void CLEAN(vector<string>);
  void CLEANALL(std::istream& input);
  void CMPSTR(vector<string>);
  void COMPARE(const string& options);
  bool DATA(std::istream& input, const aurostd::xoption& vpflow, const string& smode = "DATA", ostream& oss = cout); // DX20170901 - SGDATA + JSON //DX20210302 - added const to vpflow
  void DEBYE(const string& options);
  void DISP(const string& options, std::istream& input);
  void DIST(const string& options, std::istream& input);
  void DYNADIEL(vector<string>& argv); // CAMILO
  void EDOS(vector<string>);
  void EFFMASS(vector<string>& argv, ostream& oss = cout); // CAMILO
  void EIGCURV(const string& options, ostream& oss = cout); // CAMILO
  string EQUIVALENT(_aflags& aflags, std::istream& input, aurostd::xoption& vpflow);
  void EWALD(const string& options, std::istream& input);
  string EXTRACT_xcar(_aflags& aflags, vector<string>, string, string);
  string EXTRACT_Symmetry(_aflags& aflags, vector<string>);
  void FGROUP(_aflags& aflags, std::istream& input);
  bool FIXBANDS(_aflags& aflags, string opts);
  void FINDSYM(aurostd::xoption& vpflow, uint mode, std::istream& input);
  xstructure FRAC(std::istream& input);
  string FROZSL_VASPSETUP(vector<string> argv, int mode);
  string FROZSL_ANALYZE(std::istream& input);
  string FROZSL_INPUT();
  string FROZSL_OUTPUT();
  string GEOMETRY(std::istream& input); // CO20191110
  bool GetCollinearMagneticInfo(uint num_atoms, const string& magmom_info, vector<double>& vmag); // DX20170927 - Magnetic symmetry //DX20191107 - int to uint
  bool GetNonCollinearMagneticInfo(uint num_atoms, const string& magmom_info, vector<aurostd::xvector<double>>& vmag_noncoll); // DX20171205 - Magnetic symmetry non-collinear //DX20191107 - int to uint
  void GLASS_FORMING_ABILITY(aurostd::xoption& vpflow); // DF20190329
  void ATOMIC_ENVIRONMENT(const aurostd::xoption& vpflow); // HE20210331
  void GULP(std::istream& input);
  void HKL(const string& options, _aflags& aflags, std::istream& input);
  void HKLSearch(const string& options, _aflags& aflags, std::istream& input, const string& smode);
  bool setPOCCTOL(xstructure& xstr, const string& pocc_tol_string); // CO20181226
  bool POCC_COMMAND_LINE(aurostd::xoption& vpflow, std::istream& input, ostream& oss = cout); // CO20181226
  void ICSD_2WYCK(std::istream& input, bool SOF);
  void ICSD(vector<string> argv, std::istream& input);
  xstructure IDENTICAL(std::istream& input);
  xstructure INCELL(std::istream& input);
  xstructure INCOMPACT(std::istream& input);
  void INTPOL(const string& options);
  xstructure INWS(std::istream& input);
  void JMOL(const string& options, std::istream& input);
  void KBAND(vector<string>);
  xstructure INFLATE_LATTICE(const string& options, std::istream& input);
  xstructure INFLATE_VOLUME(const string& options, std::istream& input);
  void KPATH(std::istream& input, bool WWW);
  void KPATH(std::istream& input, double grid, bool WWW);
  xstructure KPOINTS(const string& options, std::istream& input, ostream& oss = cout);
  xstructure KPOINTS_DELTA(aurostd::xoption& vpflow, std::istream& input, ostream& oss = cout);
  void JOINSTRLIST(vector<string>);
  void MAKESTRLIST(vector<string>);
  xstructure LATTICEREDUCTION(std::istream& input);
  string LATTICE_TYPE(std::istream& input, aurostd::xoption& vpflow); // DX20200820 - added vpflow
  string LATTICE_LATTICE_TYPE(std::istream& input, aurostd::xoption& vpflow); // DX20200820
  string listPrototypeLabels(aurostd::xoption& vpflow); // DX20181004
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // START - all relevent functions for loading entries here
  // Added by Corey Oses - May 2017
  // load entries is heavily overloaded, mostly to accommodate entries separated as
  // vector<vector<vector<> > > entries (unaries vs. binaries, then species-specific, good for convex hull),
  // vector<vector<> > entries (unaries vs. binaries), OR
  // vector<> entries (all together)
  ////////////////////////////////////////////////////////////////////////////////
  string arity_string(uint arity, bool capital = false, bool plural = false); // CO20180329
  // loadEntries
  bool loadEntries(vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadEntries
  bool loadEntries(vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadEntries
  bool loadEntries(vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadEntries(vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadEntries(aurostd::xoption& vpflow, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  ////////////////////////////////////////////////////////////////////////////////
  bool loadFromCommon(aurostd::xoption& vpflow);
  ////////////////////////////////////////////////////////////////////////////////
  // load and merging LIBX
  bool loadAndMergeLIBX(string combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss = cout);
  bool loadAndMergeLIBX(string _combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss = cout);
  bool loadAndMergeLIBX(vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, string combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss = cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, string _combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ostream& oss = cout);
  bool loadAndMergeLIBX(aurostd::xoption& vpflow, vector<string>& combination, string LIB, string server, vector<vector<vector<aflowlib::_aflowlib_entry>>>& naries, ofstream& FileMESSAGE, ostream& oss = cout);
  ////////////////////////////////////////////////////////////////////////////////
  uint SubLayersRestAPILS(const string& url, vector<string>& vsuburl); // CO20200731
  ////////////////////////////////////////////////////////////////////////////////
  // loadLIBX string elements
  bool loadLIBX(string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  // loadLIBX vector elements
  bool loadLIBX(string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<aflowlib::_aflowlib_entry>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  ////////////////////////////////////////////////////////////////////////////////
  // loadLIBX string elements, organized by -naries
  bool loadLIBX(string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, string elements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  // loadLIBX_nested vector elements
  bool loadLIBX(string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ostream& oss = cout);
  bool loadLIBX(aurostd::xoption& vpflow, string LIB, vector<string>& velements, string server, vector<vector<aflowlib::_aflowlib_entry>>& entries, ofstream& FileMESSAGE, ostream& oss = cout);
  ////////////////////////////////////////////////////////////////////////////////
  // get elemental combinations (recursively)
  vector<vector<string>> elementalCombinations(const vector<string>& velements, uint nary);
  ////////////////////////////////////////////////////////////////////////////////
  // easy way to think about it:  do compounds belong to the hull?
  bool compoundsBelong(const vector<string>& velements, const string& input, ostream& oss = cout, bool clean = true, bool sort_elements = false, elements_string_type e_str_type = composition_string, bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const vector<string>& velements,
                       const string& input,
                       vector<string>& input_velements,
                       vector<double>& input_vcomposition,
                       ostream& oss = cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const vector<string>& velements,
                       const string& input,
                       ofstream& FileMESSAGE,
                       ostream& oss = cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const vector<string>& velements,
                       const string& input,
                       vector<string>& input_velements,
                       vector<double>& input_vcomposition,
                       ofstream& FileMESSAGE,
                       ostream& oss = cout,
                       bool clean = true,
                       bool sort_elements = false,
                       elements_string_type e_str_type = composition_string,
                       bool shortcut_pp_string_AFLOW_database = false);
  bool compoundsBelong(const vector<string>& velements, const vector<string>& elements, ostream& oss = cout, bool sort_elements = false);
  bool compoundsBelong(const vector<string>& velements, const vector<string>& elements, ofstream& FileMESSAGE, ostream& oss = cout, bool sort_elements = false);
  ////////////////////////////////////////////////////////////////////////////////
  // loads xstructures
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, ostream& oss = cout, bool relaxed_only = true, string path = "", bool is_url_path = false);
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, ofstream& FileMESSAGE, ostream& oss = cout, bool relaxed_only = true, string path = "", bool is_url_path = false);
  bool loadXstructures(aflowlib::_aflowlib_entry& entry, vector<string>& structure_files, ofstream& FileMESSAGE, ostream& oss = cout, bool relaxed_only = true, string path = "", bool is_url_path = false); // DX20200224
  bool loadXstructureLibEntry(const aflowlib::_aflowlib_entry& entry, xstructure& new_structure); // HE20220606
  ////////////////////////////////////////////////////////////////////////////////
  // sets default flags
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ostream& oss = cout, string input = "A", bool entry_output = true, bool silent = false);
  void defaultLoadEntriesFlags(aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& oss = cout, string input = "A", bool entry_output = true, bool silent = false);
  ////////////////////////////////////////////////////////////////////////////////
  bool prototypeMatch(string proto_database, string proto_search); // smarter than == for prototype matches, deals with 549 vs. 549.bis
  ////////////////////////////////////////////////////////////////////////////////
  // Added by Corey Oses - May 2017
  // END - all relevent functions for loading entries here
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // START - added by Corey Oses - May 2017
  // effectively logs EVERYTHING, deals with cout and logger
  void updateProgressBar(unsigned long long int current, unsigned long long int end, ostream& oss = std::cout);
  void logger(const string& filename, const string& function_name, stringstream& message, const char& type, ostream& oss = cout, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, stringstream& message, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, stringstream& message, ofstream& FileMESSAGE, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, stringstream& message, const string& directory, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename,
              const string& function_name,
              stringstream& message,
              const string& directory,
              ofstream& FileMESSAGE,
              ostream& oss = cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, stringstream& message, const _aflags& aflags, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename,
              const string& function_name,
              stringstream& message,
              const _aflags& aflags,
              ofstream& FileMESSAGE,
              ostream& oss = cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, const string& _message, const char& type = _LOGGER_MESSAGE_, ostream& oss = cout, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, const string& _message, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, const string& _message, ofstream& FileMESSAGE, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, const string& _message, const string& directory, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename, const string& function_name, const string& _message, const _aflags& aflags, ostream& oss = cout, const char& type = _LOGGER_MESSAGE_, bool silent = false, const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // overload
  void logger(const string& filename,
              const string& function_name,
              const string& _message,
              const string& directory,
              ofstream& FileMESSAGE,
              ostream& oss = cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // main function
  void logger(const string& filename,
              const string& function_name,
              const string& _message,
              const _aflags& aflags,
              ofstream& FileMESSAGE,
              ostream& oss = cout,
              const char& type = _LOGGER_MESSAGE_,
              bool silent = false,
              const string& message_metadata = _AFLOW_MESSAGE_DEFAULTS_); // main function
  // END - added by Corey Oses - May 2017
  xstructure LTCELL(const string& options, std::istream& input);
  void MagneticParameters(string _directory, ostream& oss = cout);
  xstructure MINKOWSKIBASISREDUCTION(std::istream& input);
  string MISCIBILITY(vector<string> argv);
  void MOM(std::istream& input);
  void MSI(std::istream& input);
  uint NATOMS(std::istream& input);
  //[CO20190520 - not in pflow]xmatrix<double> GetDistMatrix(const xstructure& a); //CO20171025
  string NBONDXX(std::istream& input, bool aflowlib_legacy_format = false); // CO20171025
  string NBONDXX(const xstructure& a, bool aflowlib_legacy_format = false); // CO20171025
  xstructure NAMES(vector<string>, std::istream& input);
  xstructure NANOPARTICLE(std::istream& input, const aurostd::xvector<double>& iparams);
  xstructure NIGGLI(std::istream& input);
  void NDATA(std::istream& input);
  double NNDIST(std::istream& input);
  xstructure NOORDERPARAMETER(std::istream& input);
  xstructure NOSD(std::istream& input);
  xstructure NUMNAMES(vector<string>, std::istream& input);
  uint NSPECIES(std::istream& input);
  bool OPARAMETER(vector<string>, std::istream& input);
  bool SHIRLEY(vector<string>, std::istream& input);
  bool SHIRLEY2(vector<string>, std::istream& input);
  string PEARSON_SYMBOL(std::istream& input);
  string PEARSON_SYMBOL(std::istream& input, const aurostd::xoption& vpflow); // DX20210611 - added vpflow
  bool POCCUPATION(vector<string>, std::istream& input);
  void PDB(std::istream& input);
  void PDOS(vector<string>);
  void PHONONS(_aflags& aflags, std::istream& input, const double& radius);
  void PGROUP(_aflags& aflags, std::istream& input);
  void PGROUPXTAL(_aflags& aflags, std::istream& input);
  void PGROUPK(_aflags& aflags, std::istream& input);
  void PLANEDENS(vector<string>);
  string PLATON(const string& options, std::istream& input);
  string SG(aurostd::xoption& vpflow, std::istream& input, string mode, string print);
  void STATDIEL(vector<string>& argv); // CAMILO
  bool SYMMETRY_GROUPS(_aflags& aflags, std::istream& input, aurostd::xoption& vpflow, ostream& oss = cout); // DX20170818 - Add no_scan option to all symmetry Xgroups
  void POCC(vector<string>);
  string POSCAR2AFLOWIN(std::istream& input, const string& = ""); // Modified ME20181113
  void POSCAR2WYCKOFF(std::istream& input);
  vector<string> GENERATE_CERAMICS(const vector<string>& nonmetals, const vector<string>& metals, uint metal_arity); // CO20200731
  vector<string> GENERATE_CERAMICS(const aurostd::xoption& vpflow); // CO20200731
  string GENERATE_CERAMICS_PRINT(const aurostd::xoption& vpflow); // CO20200731
  bool PSEUDOPOTENTIALS_CHECK(const aurostd::xoption& vpflow, const string& file, ostream& oss = cout); // SC20200330
  void PYTHON_MODULES(const string& modules, ostream& oss = std::cout); // ME20211103
  void PYTHON_MODULES(const string& modules, ofstream& FileMESSAGE, ostream& oss = std::cout); // ME20211103
  void PYTHON_MODULES(const vector<string>& vmodules, ostream& oss = std::cout); // ME20211103
  void PYTHON_MODULES(const vector<string>& vmodules, ofstream& FileMESSAGE, ostream& oss = std::cout); // ME20211103
  bool QMVASP(aurostd::xoption& vpflow); // vector<string> argv); //CO20180703
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
  using std::cout;
  using std::deque;
  using std::ofstream;
  using std::string;
  using std::stringstream;
  using std::vector;
  bool PROTO_PARSE_INPUT(const vector<string>& params, vector<vector<string>>& vstr, vector<vector<double>>& vnum, bool ignore_label = false, bool reverse = false); // CO20181226
  bool PROTO_TEST_INPUT(const vector<vector<string>>& vvstr, const vector<vector<double>>& vvnum, uint& nspeciesHTQC, bool patch_nspecies = false); // CO20181226
  bool sortPOCCSites(const string& p1, const string& p2); // CO20181226
  bool sortPOCCOccs(const string& occ1, const string& occ2); // CO20181226
  bool FIX_PRECISION_POCC(const string& occ, string& new_occ); // CO20181226
  void FIX_POCC_PARAMS(const xstructure& xstr, string& pocc_params); // CO20181226
  bool checkAnionSublattice(const xstructure& xstr); // CO20210201
  bool convertXStr2POCC(xstructure& xstr, const string& pocc_params, const vector<string>& vspecies, const vector<double>& vvolumes); // CO20181226
  bool POccInputs2Xstr(const string& pocc_input, aurostd::xoption& pocc_settings, xstructure& xstr, ostream& oss); // CO20211130
  bool POccInputs2Xstr(const string& pocc_input, aurostd::xoption& pocc_settings, xstructure& xstr, ofstream& FileMESSAGE, ostream& oss); // CO20211130
  xstructure PROTO_LIBRARIES(aurostd::xoption vpflow);
  bool PROTO_AFLOW(aurostd::xoption vpflow, bool flag_REVERSE); // too many options
  bool PROTO_ICSD_AFLOWIN(vector<string>& argv);
  xstructure PRIM(std::istream& input, uint mode);
  void RASMOL(const string& options, std::istream& input);
  void RBANAL(vector<string>);
  void RBDIST(vector<string>);
  xstructure RMATOM(std::istream& input, const int& iatom);
  xstructure RMCOPIES(std::istream& input);
  void RAYTRACE(vector<string>);
  xstructure SCALE(const string& options, std::istream& input);
  void RDF(const string& options, std::istream& input, bool raw_counts = false); // CO20220627
  void RDFCMP(const string& options);
  void RSM(vector<string>, std::istream& input);
  xstructure SD(vector<string>, std::istream& input);
  xstructure SETCM(std::istream& input, const aurostd::xvector<double>& cm);
  xstructure SETORIGIN(std::istream& input, const aurostd::xvector<double>& origin);
  xstructure SETORIGIN(std::istream& input, const int& natom);
  void SEWALD(vector<string>, std::istream& input);
  void SG(std::istream& input);
  void SGROUP(_aflags& aflags, std::istream& input, double radius);
  void SHELL(const string& options, std::istream& input);
  string SPECIES(std::istream& input);
  xstructure SHIFT(const string& options, std::istream& input);
  void SPLINE(vector<string>);
  void SUMPDOS(vector<string>);
  xstructure SUPERCELL(const string& options, std::istream& input);
  void SUPERCELLSTRLIST(const string& options);
  xstructure xstrSWAP(vector<string>, std::istream& input);
  xstructure VOLUME(const string& options, std::istream& input);
  string WyckoffPositions(aurostd::xoption& vpflow, std::istream& input); // DX20180807 - added wyccar to pflow //DX20210525 - changed name and generalized function
  xstructure WYCKOFF(vector<string>, std::istream& input);
  void XRAY(const string& options, std::istream& input);
  void XRAY_PEAKS(const aurostd::xoption& vpflow, std::istream& input); // CO20190409
  void READ_XRAY_DATA(const string& filename, vector<double>& v_twotheta, vector<double>& intensity); // CO20190620
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, std::istream& input); // CO20190409
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, const xstructure& str); // CO20190409
  void PRINT_XRAY_DATA_PLOT(std::istream& input, double lambda = XRAY_RADIATION_COPPER_Kalpha, const string& directory = ""); // CO20190409
  void PRINT_XRAY_DATA_PLOT(const xstructure& str, double lambda = XRAY_RADIATION_COPPER_Kalpha, const string& directory = ""); // CO20190409
  void PRINT_XRAY_DATA_PLOT(const aurostd::xoption& vpflow, const string& directory = ""); // CO20190620
  void PRINT_XRAY_DATA_PLOT(const string& filename, const string& directory = ""); // CO20190620
  void PRINT_XRAY_DATA_PLOT(const vector<double>& v_twotheta, const vector<double>& v_intensity, const string& directory = ""); // CO20190620
  void PLOT_XRAY(const aurostd::xoption& vpflow, std::istream& input); // CO20190409
  void PLOT_XRAY(const aurostd::xoption& vpflow, const xstructure& str); // CO20190409
  void PLOT_XRAY(std::istream& input, double lambda = XRAY_RADIATION_COPPER_Kalpha, const string& directory = "", bool keep_gp = false, bool force_generic_title = false); // CO20190409
  void PLOT_XRAY(const xstructure& str, double lambda = XRAY_RADIATION_COPPER_Kalpha, const string& directory = "", bool keep_gp = false, bool force_generic_title = false); // CO20190409
  void PLOT_XRAY(const aurostd::xoption& vpflow, const string& title = "", const string& directory = "", bool keep_gp = false); // CO20190620
  void PLOT_XRAY(const string& filename, const string& title = "", const string& directory = "", bool keep_gp = false); // CO20190620
  void PLOT_XRAY(const vector<double>& v_twotheta, const vector<double>& v_intensity, const string& title = "", const string& directory = "", bool keep_gp = false); // CO20190620
  void XYZ(const string& options, std::istream& input);
  void XYZINSPHERE(std::istream& input, double radius);
  void XYZWS(std::istream& input);
  void XelementPrint(const string& options, ostream& oss = cout);
  void ZVAL(const string& options);
} // namespace pflow

// aflow_pflow_print.cpp
namespace pflow {
  using std::cout;
  using std::deque;
  using std::ofstream;
  using std::string;
  using std::stringstream;
  using std::vector;
  void PrintACE(const xstructure&, ostream& oss = cout);
  void PrintAngles(xstructure str, const double& cutoff, ostream& oss = cout);
  class projdata;
  void PrintBands(const pflow::projdata& pd);
  bool PrintCHGCAR(const xstructure& str, const stringstream& chgcar_header, const vector<int>& ngrid, const vector<int>& format_dim, const vector<double>& chg_tot, const vector<double>& chg_diff, const string& output_name, ostream& oss = cout);
  void PrintChgInt(vector<aurostd::matrix<double>>& rad_chg_int, aurostd::matrix<double>& vor_chg_int, ostream& oss = cout); // CO20200404 pflow::matrix()->aurostd::matrix()
  void PrintCIF(ostream& oss, const xstructure&, int = 1, int = 1); // DX20180806 - added setting default
  void PrintClat(const aurostd::xvector<double>& data, ostream& oss = cout);
  void PrintCmpStr(const xstructure& str1, const xstructure& str2, const double& rcut, ostream& oss = cout);
  string PrintData(const xstructure& xstr, const string& smode = "DATA", filetype ftype = txt_ft, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1); // DX20210301
  string PrintData(const xstructure& xstr, aurostd::xoption& vpflow, const string& smode = "DATA", filetype ftype = txt_ft, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1); // DX20210301
  string PrintData(const xstructure& xstr,
                   xstructure& str_sp,
                   xstructure& str_sc,
                   aurostd::xoption& vpflow,
                   const string& smode = "DATA",
                   filetype ftype = txt_ft,
                   bool already_calculated = false,
                   double sym_eps = AUROSTD_MAX_DOUBLE,
                   bool no_scan = false,
                   int setting = 1); // DX20210301
  string PrintData(const xstructure& xstr,
                   xstructure& str_sym,
                   xstructure& str_sp,
                   xstructure& str_sc,
                   const string& smode = "DATA",
                   filetype ftype = txt_ft,
                   bool already_calculated = false,
                   double sym_eps = AUROSTD_MAX_DOUBLE,
                   bool no_scan = false,
                   int setting = 1); // DX20210301
  string PrintData(const xstructure& xstr,
                   xstructure& str_sym,
                   xstructure& str_sp,
                   xstructure& str_sc,
                   aurostd::xoption& vpflow,
                   const string& smode = "DATA",
                   filetype ftype = txt_ft,
                   bool already_calculated = false,
                   double sym_eps = AUROSTD_MAX_DOUBLE,
                   bool no_scan = false,
                   int setting = 1); // DX20210301
  string PrintRealLatticeData(const xstructure& xstr, const string& smode = "DATA", filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  string PrintRealLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, const string& smode = "DATA", filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  string PrintRealLatticeData(
      const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, const string& smode = "DATA", filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  string PrintLatticeLatticeData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  string PrintLatticeLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  string PrintLatticeLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  string PrintCrystalPointGroupData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  string PrintCrystalPointGroupData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210211
  string PrintCrystalPointGroupData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  string PrintReciprocalLatticeData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  string PrintReciprocalLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  string PrintReciprocalLatticeData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // HE20240221

  string PrintSuperlatticeData(const xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  string PrintSuperlatticeData(const xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209
  string PrintSuperlatticeData(const xstructure& xstr, aurostd::xoption& vpflow, aurostd::JSON::object& json, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE); // DX20210209

  void PrintDisplacements(xstructure str, const double cutoff, ostream& oss = cout);
  void PrintDistances(xstructure str, const double cutoff, ostream& oss = cout);
  void PrintEwald(const xstructure& in_str, double& epoint, double& ereal, double& erecip, double& eewald, double& eta, const double& SUMTOL, ostream& oss = cout);
  void PrintGulp(const xstructure&, ostream& oss = cout);
  string PrintSGData(xstructure& xstr, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1, bool supress_Wyckoff = false); // DX20210211
  string PrintSGData(xstructure& xstr, aurostd::xoption& vpflow, filetype ftype = txt_ft, bool standalone = true, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1, bool suppress_Wyckoff = false); // DX20210211
  string PrintSGData(xstructure& xstr,
                     aurostd::xoption& vpflow,
                     aurostd::JSON::object& json,
                     filetype ftype = txt_ft,
                     bool standalone = true,
                     bool already_calculated = false,
                     double sym_eps = AUROSTD_MAX_DOUBLE,
                     bool no_scan = false,
                     int setting = 1,
                     bool suppress_Wyckoff = false); // DX20210211

  string PrintWyckoffData(xstructure& xstr, filetype ftype = txt_ft, bool already_calculated = false, double sym_eps = AUROSTD_MAX_DOUBLE, bool no_scan = false, int setting = 1); // DX20210610
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
  using std::cout;
  using std::deque;
  using std::ofstream;
  using std::string;
  using std::stringstream;
  using std::vector;
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void GeneralizedStackingFaultEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, std::istream& input, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321
  void CleavageEnergyCalculation(const aurostd::xoption& vpflow, const xstructure& xstr_in, const _aflags& aflags, const _kflags& kflags, const _vflags& vflags, ofstream& FileMESSAGE, ostream& oss = cout); // CO20190321

  bool findClosedPackingPlane(std::istream& input); // CO20191110
  bool findClosedPackingPlane(const xstructure& xstr); // CO20191110
} // namespace pflow

namespace pflow {
  using std::cout;
  using std::deque;
  using std::ofstream;
  using std::string;
  using std::stringstream;
  using std::vector;
  void GetXray2ThetaIntensity(const xstructure& str, vector<double>& v_twotheta, vector<double>& v_intensity, double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190520
  vector<uint> GetXrayPeaks(const xstructure& str, vector<double>& v_twotheta, vector<double>& v_intensity, vector<double>& v_intensity_smooth, double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190520 //CO20190620 - v_peaks_amplitude not needed
  vector<uint> GetXrayPeaks(const vector<double>& v_twotheta, const vector<double>& v_intensity, vector<double>& v_intensity_smooth); // CO20190520  //CO20190620 - v_peaks_amplitude not needed
  void GetXray(const xstructure& str, vector<double>& dist, vector<double>& sf, vector<double>& scatt_fact, vector<double>& mass, vector<double>& twoB_vec, double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190520
  void GetXrayData(const xstructure& str,
                   vector<double>& dist,
                   vector<double>& sf,
                   vector<double>& scatt_fact,
                   vector<double>& mass,
                   vector<double>& twoB_vec,
                   vector<vector<double>>& ids,
                   aurostd::matrix<double>& data,
                   double lambda = XRAY_RADIATION_COPPER_Kalpha); // CO20190409  //CO20190620 - intmax can be grabbed later  //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<vector<int>> getvvitypesRDF(const xstructure& xstr); // CO20220627
  void GetRDF(const xstructure& xstr, aurostd::xmatrix<double>& rdf_all, const double rmax = 5.0, const int nbins = 25, bool raw_counts = false, const double sigma = 0.0, const int window_gaussian = 0); // CO20220627 - rewritten
  void GetRDFShells(const xstructure& str, const double& rmax, const int& nbins, const int& smooth_width, const aurostd::matrix<double>& rdf, aurostd::matrix<double>& rdfsh, aurostd::matrix<double>& rdfsh_loc); // CO20200404 pflow::matrix()->aurostd::matrix()
  double RdfSh_RMS(const int iaA, const int iaB, const int nsh_max, const int nt, const aurostd::matrix<double>& rdfsh_all_A, const aurostd::matrix<double>& rdfsh_all_B); // CO20200404 pflow::matrix()->aurostd::matrix()
  void CmpRDFShells(const xstructure& str_A,
                    const xstructure& str_B,
                    const aurostd::matrix<double>& rdfsh_all_A,
                    const aurostd::matrix<double>& rdfsh_all_B,
                    const int nsh,
                    vector<int>& best_match,
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
  using std::cout;
  using std::ofstream;
  using std::ostream;
  using std::string;
  using std::vector;
  class pdosdata {
  public:
    // Constructors
    pdosdata(); // default
    void PrintParams(ostream& oss, const vector<string>& Ltotnames);
    void PrintPDOS(ostream& oss, const int& sp);
    // variables
    string PDOSinfile;
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
  using std::cout;
  using std::ifstream;
  using std::ofstream;
  using std::ostream;
  using std::string;
  using std::vector;
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
    vector<double> viewdir;
    int viewdir_s;
    vector<double> updir;
    int updir_s;
    double zoom;
    double aspectratio;
    double antialiasing;
    double raydepth;
    vector<double> center;
    vector<double> center_guide;
    int center_s;
    vector<double> background;
    aurostd::matrix<double> lightcenter; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> lightrad;
    aurostd::matrix<double> lightcolor; // CO20200404 pflow::matrix()->aurostd::matrix()
    // Sphere texture variables (ambient, diffuse, specular, opacity)
    aurostd::matrix<double> sphtex_tex; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> sphtex_tex_def;
    aurostd::matrix<double> sphtex_color; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> sphtex_color_def;
    aurostd::matrix<double> sphtex_phong; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> sphtex_phong_def;
    vector<string> sphtex_names;
    vector<double> sph_rad;
    // Plane variables
    int plane;
    int plane_s;
    vector<double> plane_center;
    vector<double> plane_normal;
    vector<double> plane_color;
    vector<double> planetex_tex;
    vector<double> plane_center_def;
    vector<double> plane_normal_def;
    vector<double> plane_color_def;
    vector<double> planetex_tex_def;
    int plane_center_s;
    int plane_normal_s;
    int plane_color_s;
    int planetex_tex_s;
    double sph_rad_def;
    string shading;
    string outfile;
    aurostd::matrix<double> sc; // CO20200404 pflow::matrix()->aurostd::matrix()
    int sc_s;
    int calc_type;
    vector<string> input_files;
    int first_set;
    string insert_file;
    vector<double> rotation;
    vector<double> struct_origin;
    int struct_origin_s;
  };

  void SetRTParams(xstructure& str, pflow::rtparams& rtp);
  vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(ifstream& outcar_inf, ifstream& xdatcar_inf);
  void GetDatFromOutcar(vector<aurostd::matrix<double>>& lat_vec, deque<int>& num_each_type, ifstream& outcar_inf); // CO20200404 pflow::matrix()->aurostd::matrix()
  void GetDatFromXdatcar(vector<aurostd::matrix<double>>& fpos_vec, ifstream& xdatcar_inf); // CO20200404 pflow::matrix()->aurostd::matrix()
  vector<xstructure> GetStrVecFromOUTCAR_XDATCAR(ifstream& outcar_inf, ifstream& xdatcar_inf);
  void PrintStrVec(const vector<xstructure>& str_vec, ostream& outf);
  void ReadInRTParams(ifstream& rtinfile, pflow::rtparams& rtp);
  void ReadInStrVec(vector<xstructure>& str_vec, ifstream& strlist_inf);
  void JoinStrVec(vector<xstructure> str_vec_1, vector<xstructure> str_vec_2, vector<xstructure>& str_vec_tot);
  void SetStrFromRTParams(xstructure& str, pflow::rtparams& rtp);
  void SuperCellStrVec(vector<xstructure>& str_vec, const aurostd::matrix<double>& sc); // CO20200404 pflow::matrix()->aurostd::matrix()
  void UpDateRTParams(pflow::rtparams& rtp, const int& istr, int nstr);
  void SetRTParams(xstructure& str, pflow::rtparams& rtp);
  void GetRTDatFile(xstructure str, const pflow::rtparams& rtp, std::ostringstream& rtdat_file);
  string PrintRTDatFile(std::ostringstream& rtdat_file, const pflow::rtparams& rt_params);
  string CreateRTtgaFile(const string& datfile, const pflow::rtparams& rt_params);
  string CreateRTjpgFile(const string& tgafile, const pflow::rtparams& rt_params);
  void GetRTencFile(const pflow::rtparams& rtp, const int nim, std::ostringstream& os);
  string PrintRTencFile(const pflow::rtparams& rt_params, std::ostringstream& rtenc_file);
  string CreateRTmpgFile(const pflow::rtparams& rt_params, const string& encfile);
  void RayTraceManager(vector<string>);

} // namespace pflow

// ---------------------------------------------------------------------------
// PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PROJDATA PRO

namespace pflow {
  using std::cout;
  using std::ofstream;
  using std::ostream;
  using std::string;
  using std::vector;
  class projdata {
  public:
    // Constructors
    projdata(); // default
    void Print(ostream& outf);
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
    vector<int> num_each_type;
    aurostd::matrix<double> wfermi_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> wfermi_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<double> wkpt;
    vector<aurostd::matrix<aurostd::matrix<std::complex<double>>>> pdat_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<aurostd::matrix<std::complex<double>>>> pdat_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lmtot_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<aurostd::matrix<double>> occ_vs_ion_kpt_bnd_lmtot_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_kpt_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_kpt_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_kpt_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_kpt_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_bnd_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_bnd_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_bnd_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<aurostd::matrix<double>> occ_vs_ion_bnd_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_lm_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_lm_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_l_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> occ_vs_ion_l_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> ener_k_b_u; // CO20200404 pflow::matrix()->aurostd::matrix()
    aurostd::matrix<double> ener_k_b_d; // CO20200404 pflow::matrix()->aurostd::matrix()
    int sp;
    int rspin;
    aurostd::matrix<double> kpts; // CO20200404 pflow::matrix()->aurostd::matrix()
    vector<string> LMnames;
    vector<string> Lnames;
    vector<string> LLMnames;
    string PROOUTinfile;
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

  bool ReadCHGCAR(xstructure& str, std::stringstream& chgcar_header, std::vector<int>& ngrid, std::vector<int>& format_dim, std::vector<double>& chg_tot, std::vector<double>& chg_diff, std::stringstream& chgcar_ss, std::ostream& oss = std::cout);
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

// ME20190628 - prettyPrintCompound from CHULL
namespace pflow {
  // Precision for pretty printing
  const int COEF_PRECISION = 4;

  std::string prettyPrintCompound(const std::string&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20190629
  std::string prettyPrintCompound(const std::vector<std::string>&, const std::vector<uint>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //DX20200727 - uint variant
  std::string prettyPrintCompound(const std::vector<std::string>&, const std::vector<double>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20190629
  std::string prettyPrintCompound(const std::vector<std::string>&, const aurostd::xvector<uint>&, vector_reduction_type vred = gcd_vrt, bool = true, filetype ftype = latex_ft); // char=_latex_  //CO20200727 - uint variant
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
    // NECESSARY PUBLIC CLASS METHODS - START
    // constructors - START
    AQueue(std::ostream& oss = std::cout);
    AQueue(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    AQueue(const aurostd::xoption& vpflow, std::ostream& oss = std::cout);
    AQueue(const aurostd::xoption& vpflow, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    AQueue(const AQueue& b);
    // constructors - STOP
    ~AQueue();
    const AQueue& operator=(const AQueue& other);
    void clear();
    // NECESSARY PUBLIC CLASS METHODS - STOP

    // general attributes
    bool m_initialized;
    aurostd::xoption m_flags;
    queue_system m_qsys;
    std::vector<APartition> m_partitions;
    std::vector<ANode> m_nodes;
    std::vector<AJob> m_jobs;

    // initialization methods
    bool initialize(std::ostream& oss);
    bool initialize(std::ofstream& FilMESSAGE, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, std::ostream& oss);
    bool initialize(const aurostd::xoption& vpflow, std::ofstream& FilMESSAGE, std::ostream& oss);
    bool initialize();
    bool initialize(const aurostd::xoption& vpflow);

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
    // NECESSARY private CLASS METHODS - START
    void free();
    void copy(const AQueue& b);
    // NECESSARY END CLASS METHODS - END

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
