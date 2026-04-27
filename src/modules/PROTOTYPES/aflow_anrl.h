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

  bool vproto2tokens(std::string proto,
                     std::string& label,
                     uint& nspecies,
                     uint& natoms,
                     uint& spacegroup,
                     uint& nunderscores,
                     uint& nparameters,
                     std::string& Pearson_symbol,
                     std::string& params,
                     std::string& Strukturbericht,
                     std::string& prototype,
                     std::string& dialect);
  // ---------------------------------------------------------------------------
  // functions to determine atomic positions from Wyckoff and parameters
  std::vector<uint> extractStoichiometry(std::string& anrl_label);
  // ---------------------------------------------------------------------------
  // checking functions
  bool PrototypeANRL_Consistency(uint vparameters_size,
                                 uint proto_nparameters,
                                 std::string proto_prototype,
                                 std::string proto_label,
                                 std::string proto_Strukturbericht,
                                 std::string proto_Pearson_symbol,
                                 uint proto_spacegroup,
                                 std::string proto_params,
                                 uint print_mode); // DX20180710 - added print_mode //DX20200207 - oss no longer needed
  // ---------------------------------------------------------------------------
  // helper functions to determine label and internal degrees of freedom
  std::string groupedWyckoffPosition2ANRLString(const std::vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize);
  std::vector<std::string> getANRLLatticeParameterString(char& lattice_type);
  std::vector<double> getANRLLatticeParameterValuesFromWyccar(const std::vector<std::string>& wyccar_ITC, char lattice_type, char lattice_centering, uint setting); // DX20191031
  std::vector<double> getANRLLatticeParameterValuesFromABCAngles(const xstructure& xstr, char lattice_type, char lattice_centering, uint setting); // DX20191031
  std::vector<double> getANRLLatticeParameterValues(const std::vector<double>& all_lattice_parameters, char lattice_type, char lattice_centering, uint setting); // DX20191031
  uint getANRLSettingChoice(int spacegroup); // DX20191031 - removed reference
  // ---------------------------------------------------------------------------
  // std::map structure to label and internal degrees of freedom
  std::string structure2anrl(std::istream& input, aurostd::xoption& vpflow); // xoption
  std::string structure2anrl(xstructure& xstr, bool recalculate_symmetry = true); // use default options //DX20191031 - added recalculate_symmetry
  std::string structure2anrl(xstructure& xstr, double tolerance); // specify symmetry tolerance //CO20190520 - removed pointers for bools and doubles, added const where possible
  std::string structure2anrl(xstructure& xstr, uint setting); // specify setting
  std::string structure2anrl(xstructure& xstr,
                             double tolerance,
                             uint setting,
                             bool recalculate_symmetry = true,
                             bool print_element_names = false,
                             bool print_atomic_numbers = false); // main std::function //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190829 - added recalculate_symmetry //DX20191031 - removed reference //DX20210622 - added printing options
  // ---------------------------------------------------------------------------
  // generic prototype generator (main std::function)
  xstructure PrototypeANRL_Generator(std::string& label, std::string& parameters, std::deque<std::string>& vatomX, std::deque<double>& vvolumeX, std::ostream& logstream = std::cout, bool silence_logger = true); // DX20200528 - command line = no logger
  xstructure PrototypeANRL_Generator(std::string& label, std::string& parameters, std::deque<std::string>& vatomX, std::deque<double>& vvolumeX, std::ofstream& FileMESSAGE, std::ostream& logstream = std::cout, bool silence_logger = false); // DX20200528 - internal = logger
  xstructure PrototypeANRL_Generator_UID(const std::string& uid, std::string& parameters, const std::string& permutation, std::deque<std::string>& vatomX, std::deque<double>& vvolumeX);
  // ---------------------------------------------------------------------------
  // [OLD] hard-coded generator (requires ANRL/ subdirectory)
  xstructure PrototypeANRL(std::ostream& oss, std::string label, std::string parameters, std::deque<std::string>& vatomX, std::deque<double>& vvolumeX, double volume_in, int mode, bool flip_option);
  // ---------------------------------------------------------------------------
  // get lattice functions (primitive/conventional) - ITC/ANRL standard
  aurostd::xmatrix<double> getLattice(const std::string& lattice_and_centering, const char& space_group_letter, const std::vector<double>& lattice_parameter_values, uint mode = 0);
  aurostd::xmatrix<double> getTriclinicLattice(const std::vector<double>& lattice_parameter_values, uint mode = 0);
  aurostd::xmatrix<double> getMonoclinicLattice(const std::string& lattice_and_centering, const std::vector<double>& lattice_parameter_values, uint mode = 0);
  aurostd::xmatrix<double> getOrthorhombicLattice(const std::string& lattice_and_centering, const char& space_group_letter, const std::vector<double>& lattice_parameter_values, uint mode = 0);
  aurostd::xmatrix<double> getTetragonalLattice(const std::string& lattice_and_centering, const std::vector<double>& lattice_parameter_values, uint mode = 0);
  aurostd::xmatrix<double> getHexagonalLattice(const std::string& lattice_and_centering, const std::vector<double>& lattice_parameter_values, uint mode = 0);
  aurostd::xmatrix<double> getCubicLattice(const std::string& lattice_and_centering, const std::vector<double>& lattice_parameter_values, uint mode = 0);
  xstructure rhl2hex(const xstructure& str, double& a, double& c);

  // ---------------------------------------------------------------------------
  // functions to determine atomic positions from Wyckoff and parameters
  std::deque<_atom> getAtomsFromWyckoff(const std::vector<wyckoffsite_ITC>& Wyckoff_positions, const aurostd::xmatrix<double>& lattice_conventional);
  std::vector<std::string> determineWyckoffVariables(std::vector<wyckoffsite_ITC>& Wyckoff_positions);
  void applyWyckoffValues(const std::vector<double>& Wyckoff_parameter_values, std::vector<wyckoffsite_ITC>& Wyckoff_positions); // perhaps add this as a method to wyckoffsite_ITC?
  bool containsDuplicateWyckoffCoordinate(const std::vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered = false);
  std::vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(const std::vector<std::string>& Wyckoff_tokens, const std::vector<std::string>& species, uint space_group_number, int setting = SG_SETTING_ANRL);

  // ---------------------------------------------------------------------------
  // checking functions
  double specialCaseSymmetryTolerances(const std::string& label_input);
  bool isSpecialCaseEquivalentPrototypes(const std::vector<std::string>& labels_matched);
  bool structureAndLabelConsistent(const xstructure& _xstr,
                                   const std::string& label_input,
                                   std::string& label_and_params_calculated,
                                   double tolerance_sym_input = AUROSTD_MAX_DOUBLE); // DX20201105 - added tolerance input
} // namespace anrl

// Symbolic functions are defined in aflow_symbolic.cpp

#endif
