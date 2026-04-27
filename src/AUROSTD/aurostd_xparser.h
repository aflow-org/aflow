// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XPARSER_H_
#define _AUROSTD_XPARSER_H_

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "aurostd.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xvector.h"

// compound specification is how a compound is specified
// composition (Mn2Pt3) is ORTHOGONAL to pseudopotential string (Mn_pvPt)
// for instance, H1.25 can be a pseudopotential and NOT a composition
enum elements_string_type {
  composition_string,
  pp_string,
};

// CO20190712 - see VASP_PseudoPotential_CleanName_InPlace() in aflow_ivasp.cpp
const std::string CAPITAL_LETTERS_PP_LIST =
    "_GW2"    // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
    ",_GW"    // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/As_GW
    ",_ZORA"  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pt_ZORA
    ",_LDApU" // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Zn_sv_LDApU
    ",_AE"    // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
    ",_NC2"   // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_NC2
    ",_200eV"
    "";

namespace aurostd {
  void VASP_PseudoPotential_CleanName_InPlace(std::string& species, bool capital_letters_only = false, bool remove_floats = true); // CO20190712  //CO20210623 - added remove_floats
  std::string VASP_PseudoPotential_CleanName(const std::string& speciesIN);
  ////////////////////////////////////////////////////////////////////////////////
  void elementsFromCompositionString(const std::string& input);  // CO20190712
  template <class utype> void elementsFromCompositionString(const std::string& input, std::vector<std::string>& velements, std::vector<utype>& vcomposition); // CO20190712
  void elementsFromPPString(const std::string& input, std::vector<std::string>& velements, bool keep_pp = false); // CO20190712
  ////////////////////////////////////////////////////////////////////////////////
  // returns UNSORTED vector<string> from string
  std::vector<std::string> getElements(const std::string& input); // CO20190712
  std::vector<std::string> getElements(const std::string& input, elements_string_type e_str_type, bool clean = true,
    bool sort_elements = false, bool keep_pp = false, std::ostream& oss = std::cout);
  template <class utype>
  std::vector<std::string> getElements(const std::string& input, std::vector<utype>& vcomposition, bool clean = true,
    bool sort_elements = false, bool keep_pp = false, std::ostream& oss = std::cout);
  std::vector<std::string> getElements(const std::string& input, elements_string_type e_str_type, std::ofstream& FileMESSAGE,
    bool clean = true, bool sort_elements = false, bool keep_pp = false, std::ostream& oss = std::cout);
  template <class utype>
  std::vector<std::string> getElements(const std::string& input, std::vector<utype>& vcomposition, elements_string_type e_str_type,
    bool clean = true, bool sort_elements = false, bool keep_pp = false, std::ostream& oss = std::cout);
  template <class utype>
  std::vector<std::string> getElements(
      const std::string& input, std::vector<utype>& vcomposition, elements_string_type e_str_type, std::ofstream& FileMESSAGE,
      bool clean = true, bool sort_elements = false, bool keep_pp = false, std::ostream& oss = std::cout);
} // namespace aurostd

namespace aurostd {
  std::vector<std::string> extractJsonKeysAflow(const std::string& json);
  std::string extractJsonValueAflow(const std::string& json, std::string key);
  std::vector<std::string> extractJsonVectorAflow(const std::string& json, std::string key); // SD20220504
  std::vector<std::vector<std::string>> extractJsonMatrixAflow(const std::string& json, std::string key); // SD20220504
} // namespace aurostd

#endif // _AUROSTD_XPARSER_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
