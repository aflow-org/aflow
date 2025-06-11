// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_ANRL_H_
#define _AFLOW_ANRL_H_

#include <deque>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"

#include "aflow.h"
#include "aflow_defs.h"
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

using aurostd::xmatrix;
using std::deque;
using std::string;
using std::vector;

// printing modes
#define _PROTO_GENERATOR_GEOMETRY_FILE_ 0
#define _PROTO_GENERATOR_EQUATIONS_ONLY_ 1
#define _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_ 2

// ---------------------------------------------------------------------------
// below are functions limited to the aflow_anrl.cpp file, if you want to
// include functions into other parts of aflow, put the functions in aflow.h
// (since the SYMBOLICC++ is all header files, i.e., inline definitions, we
// cannot import the aflow_anrl.h file into other files)
namespace anrl {
  using std::cout;
  using std::deque;
  using std::ofstream;
  using std::ostream;
  using std::string;
  using std::vector;
  // ---------------------------------------------------------------------------

  /// @brief provides the prototype data and avoid multiple parsings of the JSON source
  /// @authors
  /// @mod{HE,20250515,created}
  struct ProtoData {
    aurostd::JSON::object lookup;
    aurostd::JSON::object content;
    static ProtoData& get() {
      static ProtoData instance;
      return instance;
    }

  private:
    ProtoData() {
      lookup = aurostd::JSON::loadString(aurostd::EmbData::get_content("lookup.json", "PROTO"));
      content = aurostd::JSON::loadString(aurostd::EmbData::get_content("data.json", "PROTO"));
    }
  };

  std::string getPrototypeUID(const std::string& search_string, bool allow_legacy = true);

  bool vproto2tokens(string proto, string& label, uint& nspecies, uint& natoms, uint& spacegroup, uint& nunderscores, uint& nparameters, string& Pearson_symbol, string& params, string& Strukturbericht, string& prototype, string& dialect);
  // ---------------------------------------------------------------------------
  // functions to determine atomic positions from Wyckoff and parameters
  vector<uint> extractStoichiometry(string& anrl_label);
  // ---------------------------------------------------------------------------
  // checking functions
  bool PrototypeANRL_Consistency(
      uint vparameters_size, uint proto_nparameters, string proto_prototype, string proto_label, string proto_Strukturbericht, string proto_Pearson_symbol, uint proto_spacegroup, string proto_params, uint print_mode); // DX20180710 - added print_mode //DX20200207 - oss no longer needed
  // ---------------------------------------------------------------------------
  // helper functions to determine label and internal degrees of freedom
  string groupedWyckoffPosition2ANRLString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize);
  vector<string> getANRLLatticeParameterString(char& lattice_type);
  vector<double> getANRLLatticeParameterValuesFromWyccar(const vector<string>& wyccar_ITC, char lattice_type, char lattice_centering, uint setting); // DX20191031
  vector<double> getANRLLatticeParameterValuesFromABCAngles(const xstructure& xstr, char lattice_type, char lattice_centering, uint setting); // DX20191031
  vector<double> getANRLLatticeParameterValues(const vector<double>& all_lattice_parameters, char lattice_type, char lattice_centering, uint setting); // DX20191031
  uint getANRLSettingChoice(int spacegroup); // DX20191031 - removed reference
  // ---------------------------------------------------------------------------
  // map structure to label and internal degrees of freedom
  string structure2anrl(std::istream& input, aurostd::xoption& vpflow); // xoption
  string structure2anrl(xstructure& xstr, bool recalculate_symmetry = true); // use default options //DX20191031 - added recalculate_symmetry
  string structure2anrl(xstructure& xstr, double tolerance); // specify symmetry tolerance //CO20190520 - removed pointers for bools and doubles, added const where possible
  string structure2anrl(xstructure& xstr, uint setting); // specify setting
  string structure2anrl(xstructure& xstr,
                        double tolerance,
                        uint setting,
                        bool recalculate_symmetry = true,
                        bool print_element_names = false,
                        bool print_atomic_numbers = false); // main function //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190829 - added recalculate_symmetry //DX20191031 - removed reference //DX20210622 - added printing options
  // ---------------------------------------------------------------------------
  // generic prototype generator (main function)
  xstructure PrototypeANRL_Generator(string& label, string& parameters, deque<string>& vatomX, deque<double>& vvolumeX, ostream& logstream = cout, bool silence_logger = true); // DX20200528 - command line = no logger
  xstructure PrototypeANRL_Generator(string& label, string& parameters, deque<string>& vatomX, deque<double>& vvolumeX, ofstream& FileMESSAGE, ostream& logstream = cout, bool silence_logger = false); // DX20200528 - internal = logger
  xstructure PrototypeANRL_Generator_UID(const string& uid, string& parameters, deque<string>& vatomX, deque<double>& vvolumeX);
  // ---------------------------------------------------------------------------
  // [OLD] hard-coded generator (requires ANRL/ subdirectory)
  xstructure PrototypeANRL(ostream& oss, string label, string parameters, deque<string>& vatomX, deque<double>& vvolumeX, double volume_in, int mode, bool flip_option);
  // ---------------------------------------------------------------------------
  // get lattice functions (primitive/conventional) - ITC/ANRL standard
  xmatrix<double> getLattice(const string& lattice_and_centering, const char& space_group_letter, const vector<double>& lattice_parameter_values, uint mode = 0);
  xmatrix<double> getTriclinicLattice(const vector<double>& lattice_parameter_values, uint mode = 0);
  xmatrix<double> getMonoclinicLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode = 0);
  xmatrix<double> getOrthorhombicLattice(const string& lattice_and_centering, const char& space_group_letter, const vector<double>& lattice_parameter_values, uint mode = 0);
  xmatrix<double> getTetragonalLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode = 0);
  xmatrix<double> getHexagonalLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode = 0);
  xmatrix<double> getCubicLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode = 0);
  xstructure rhl2hex(const xstructure& str, double& a, double& c);

  // ---------------------------------------------------------------------------
  // functions to determine atomic positions from Wyckoff and parameters
  deque<_atom> getAtomsFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_positions, const xmatrix<double>& lattice_conventional);
  vector<string> determineWyckoffVariables(vector<wyckoffsite_ITC>& Wyckoff_positions);
  void applyWyckoffValues(const vector<double>& Wyckoff_parameter_values, vector<wyckoffsite_ITC>& Wyckoff_positions); // perhaps add this as a method to wyckoffsite_ITC?
  bool containsDuplicateWyckoffCoordinate(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered = false);
  vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(const vector<string>& Wyckoff_tokens, const vector<string>& species, uint space_group_number, int setting = SG_SETTING_ANRL);

  // ---------------------------------------------------------------------------
  // checking functions
  double specialCaseSymmetryTolerances(const string& label_input);
  bool isSpecialCaseEquivalentPrototypes(const vector<string>& labels_matched);
  bool structureAndLabelConsistent(const xstructure& _xstr,
                                   const string& label_input,
                                   string& label_and_params_calculated,
                                   double tolerance_sym_input = AUROSTD_MAX_DOUBLE); // DX20201105 - added tolerance input
} // namespace anrl

// Symbolic functions are defined in aflow_symbolic.cpp

#endif
