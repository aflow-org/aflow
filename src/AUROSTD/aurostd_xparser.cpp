// ***************************************************************************
// *                                                                         *
// *                  Aflow  - Duke University 2003-2024                     *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_XPARSER_CPP_
#define _AUROSTD_XPARSER_CPP_

#include "aurostd_xparser.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "aurostd.h"
#include "aurostd_automatic_template.h"
#include "aurostd_defs.h"
#include "aurostd_time.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xscalar.h"
#include "aurostd_xvector.h"

#include "aflow_aflowrc.h" // todo required for some defs
#include "aflow_xhost.h" // todo required for XHOST.DEBUG use
#include "flow/aflow_pflow.h" // todo required for some logging

using std::cerr;
using std::endl;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::string;
using std::stringstream;
using std::vector;

namespace aurostd {

  void VASP_PseudoPotential_CleanName_InPlace(string& species, bool capital_letters_only, bool remove_floats) { // CO20190712  //CO20210623 - added remove_floats
    // WARNING: to anyone adding to this list, BE CAREFUL to avoid adding entries that contain capital letters
    // they must be added to CAPITAL_LETTERS_PP_LIST in aurostd.h
    // these pp suffixes cause problems when parsing compounds (capital letters)
    //  todo perhaps this and the elements function should be moved to aflow

    vector<string> vCAPITAL_LETTERS_PP;
    aurostd::string2tokens(CAPITAL_LETTERS_PP_LIST, vCAPITAL_LETTERS_PP, ",");
    for (size_t i = 0; i < vCAPITAL_LETTERS_PP.size(); i++) {// capital letter ones to watch out for when parsing compounds
      aurostd::RemoveSubStringInPlace(species, vCAPITAL_LETTERS_PP[i]);
    }

    if (capital_letters_only == false) {
      // from AFLOW.org database //CO20210315 - must come before .5 (removed below)
      aurostd::RemoveSubStringInPlace(species, "pot_LDA/");
      aurostd::RemoveSubStringInPlace(species, "pot_GGA/");
      aurostd::RemoveSubStringInPlace(species, "pot_PBE/");
      aurostd::RemoveSubStringInPlace(species, "potpaw_LDA/");
      aurostd::RemoveSubStringInPlace(species, "potpaw_GGA/");
      aurostd::RemoveSubStringInPlace(species, "potpaw_PBE/");
      aurostd::RemoveSubStringInPlace(species, "potpaw_LDA.54/");
      aurostd::RemoveSubStringInPlace(species, "potpaw_PBE.54/");

      // general database  //CO20210315 - must come before .5 (removed below)
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POT_LDA + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POT_GGA + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POT_PBE + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN + "/");
      aurostd::RemoveSubStringInPlace(species, DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN + "/");

      aurostd::RemoveSubStringInPlace(species, "_old");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Si_h_old
      aurostd::RemoveSubStringInPlace(species, ".old");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Mg_pv.old
      aurostd::RemoveSubStringInPlace(species, "_vnew");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pd_vnew
      aurostd::RemoveSubStringInPlace(species, "_new2");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Ti_sv_new2
      aurostd::RemoveSubStringInPlace(species, "_new");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Au_new

      aurostd::RemoveSubStringInPlace(species, "_pvf");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Cu_pvf
      aurostd::RemoveSubStringInPlace(species, "_rel");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pb_d_rel
      aurostd::RemoveSubStringInPlace(species, "_ref");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Ge_d_GW_ref
      aurostd::RemoveSubStringInPlace(species, "_local");  // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/C_local
      aurostd::RemoveSubStringInPlace(species, "_nopc");  // CO20190712 - potpaw_LDA/potpaw_PBE.20100505/Si_nopc
      aurostd::RemoveSubStringInPlace(species, ".nrel");  // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Ga_pv_GW.nrel
      aurostd::RemoveSubStringInPlace(species, "_nr");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/C_h_nr
      aurostd::RemoveSubStringInPlace(species, "_nc");  // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/H_nc_GW
      aurostd::RemoveSubStringInPlace(species, "_n");  // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_GW_n
      aurostd::RemoveSubStringInPlace(species, "_parsv");  // CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Mg_pv_parsv_GW
      aurostd::RemoveSubStringInPlace(species, "_sv2");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Li_sv2
      aurostd::RemoveSubStringInPlace(species, "_sv");
      aurostd::RemoveSubStringInPlace(species, "_vs"); // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/N_vs
      aurostd::RemoveSubStringInPlace(species, "_pv");
      aurostd::RemoveSubStringInPlace(species, "_dr");  // CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.20100505/Pb_dr
      aurostd::RemoveSubStringInPlace(species, "_d3");  // CO20190712 - BEFORE _d //potpaw_PBE/potpaw_PBE.20100506/Ge_d3
      aurostd::RemoveSubStringInPlace(species, "_d2");  // CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.05May2010/As_d2_GW
      aurostd::RemoveSubStringInPlace(species, "_d");
      aurostd::RemoveSubStringInPlace(species, "_soft");  // CO20190712 - BEFORE _s
      aurostd::RemoveSubStringInPlace(species, "_s");
      aurostd::RemoveSubStringInPlace(species, "_h");
      aurostd::RemoveSubStringInPlace(species, "_f");  // CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Cu_f
      aurostd::RemoveSubStringInPlace(species, "_af"); // CO20191110 - SHACHAR aflow pp

      aurostd::RemoveSubStringInPlace(species, "_1");
      aurostd::RemoveSubStringInPlace(species, "_2");
      aurostd::RemoveSubStringInPlace(species, "_3");

      // CO20210623 - selectively remove floats, this might interfere with extracting composition
      if (remove_floats) {
        aurostd::RemoveSubStringInPlace(species, "1.75"); // CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.75
        aurostd::RemoveSubStringInPlace(species, "1.66"); // CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.66
        aurostd::RemoveSubStringInPlace(species, "1.33"); // CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H1.33
        aurostd::RemoveSubStringInPlace(species, "1.25"); // CO20190712 - before all other decimal numbers
        aurostd::RemoveSubStringInPlace(species, "1.5"); // CO20190712 - potpaw_PBE/potpaw_PBE.06May2010/H1.5
        aurostd::RemoveSubStringInPlace(species, ".75");  // CO20190712 - before 0.5
        aurostd::RemoveSubStringInPlace(species, ".25");  // CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.25
        aurostd::RemoveSubStringInPlace(species, ".66"); // CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.66
        aurostd::RemoveSubStringInPlace(species, ".33"); // CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.33
        aurostd::RemoveSubStringInPlace(species, ".42"); // CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.42
        aurostd::RemoveSubStringInPlace(species, ".58"); // CO20190712 - before 0.5 //potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.58
        aurostd::RemoveSubStringInPlace(species, ".5");
      }

      aurostd::RemoveSubStringInPlace(species, "+1");
      aurostd::RemoveSubStringInPlace(species, "+3");
      aurostd::RemoveSubStringInPlace(species, "+5");
      aurostd::RemoveSubStringInPlace(species, "+7");
      aurostd::RemoveSubStringInPlace(species, "-1");
      aurostd::RemoveSubStringInPlace(species, "-3");
      aurostd::RemoveSubStringInPlace(species, "-5");
      aurostd::RemoveSubStringInPlace(species, "-7");

      aurostd::RemoveSubStringInPlace(species, "__"); // CO20190712 - BEFORE _ - potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW__
      aurostd::RemoveSubStringInPlace(species, "_");  // CO20190712  //potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW_
    }
  }

  // use only as a supplement for getElements(), do NOT use outside
  // this assumes a very simple Mn2Pd5, non-stoich is ok (e.g., Mn2.5Pd5)
  // no pseudo potential specification
  // no junk at the end (_ICSD_, :LDAU2, :PAW_PBE, .OLD, etc.), pre-process before
  // this is FASTER than getElements(), but not as robust for general input (specialized)
  void elementsFromCompositionString(const string& input, vector<string>& velements) {
    vector<double> vcomposition;
    return elementsFromCompositionString(input, velements, vcomposition);
  }  // CO20190712
  template <class utype> void elementsFromCompositionString(const string& input, vector<string>& velements, vector<utype>& vcomposition) { // CO20190712
    const bool LDEBUG = (false || XHOST.DEBUG);
    velements.clear();
    vcomposition.clear();  // ME20190628

    //////////////////////////////////////////////////////////////////////////////
    // START Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " original input=" << input << endl;
    }

    // CO20180409 - running through input twice, no need, simply check at the end
    // uint numberOfElements = 0;
    // for (size_t i = 0; i < input.size(); i++) {
    //   if(isupper(input[i])) {
    //     numberOfElements++;
    //   }
    // }
    // if(numberOfElements == 0) {
    //   pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    //   return velements;
    // }

    //////////////////////////////////////////////////////////////////////////////
    // END Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Parsing input
    //////////////////////////////////////////////////////////////////////////////

    // CO20180316 - fixed this function to be simpler, too complicated before
    string auxstr;
    for (size_t i = 0; i < input.size(); i++) {
      if (isupper(input[i])) {
        auxstr.clear();
        auxstr += input[i++];
        while (((i < input.size()) && (input[i] >= 'a' && input[i] <= 'z'))) {
          auxstr += input[i++];
        }
        i--;  // very important since we increase at the top of the loop (equivalent to i+j-1)
        // while (((i < input.size()) && isalpha(input[i]) && !isupper(input[i]))){auxstr+=input[i++];}
        // while(!(clean && !isalpha(input[i]))) //(input[i]=='_' || input[i]==':' || input[i]=='.' || isdigit(input[i]))))
        // isalpha() saves us from all issues with VASP_PseudoPotential_CleanName() except, e.g., potpaw_PBE/Na, we took care of that above
        // if(clean)
        //{ //CO20200106 - patching for auto-indenting
        //   auxstr = KBIN::VASP_PseudoPotential_CleanName(auxstr);  //fix vasp pp
        //   //CO20180409 - again, no need to run through essentially a third time, we already cut at these characters
        //   //look for bad characters and cut the string
        //   //for(size_t j=1;j<auxstr.size();j++){
        //   //  if(auxstr[j]=='_' || auxstr[j]==':' || isdigit(auxstr[j])){auxstr=auxstr.substr(0,j);break;}  //fix aflow stuff like ':'
        //   //}
        // }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " element found: " << auxstr << endl;
        }
        velements.push_back(auxstr);
        // ME20190628 - get composition, too
      } else if ((input[i] >= '0' && input[i] <= '9') || (input[i] == '.')) {  // CO20190712 - just in case we have H.25 (not good form but try to catch anyway, never produced by aflow automatically)
        auxstr.clear();
        auxstr += input[i++];
        while ((i < input.size()) && ((input[i] >= '0' && input[i] <= '9') || (input[i] == '.'))) {
          auxstr += input[i++];
        }
        i--;
        if (LDEBUG) {
          std::cerr << __AFLOW_FUNC__ << " found element count: " << auxstr << " of element " << (velements.size() - 1) << ".";
          if (vcomposition.size() != velements.size()) {
            std::cerr << " Will add ones to elements " << vcomposition.size() << " to " << (velements.size() - 2) << ".";
          }
          std::cerr << std::endl;
        }
        // Add implicit ones
        for (size_t i = vcomposition.size(); i < velements.size() - 1; i++) {
          vcomposition.push_back((utype) 1.0);
        }
        vcomposition.push_back(aurostd::string2utype<utype>(auxstr));
      }
    }
    // Add implicit ones
    for (size_t i = vcomposition.size(); i < velements.size(); i++) {
      vcomposition.push_back((utype) 1.0);
    }
  }

  // use only as a supplement for getElements(), do NOT use outside
  // this assumes Mn_pvPt
  // no composition information
  // no junk at the end (_ICSD_, :LDAU2, :PAW_PBE, .OLD, etc.), pre-process before
  void elementsFromPPString(const string& input, vector<string>& velements, bool keep_pp) { // CO20190712
    const bool LDEBUG = (false || XHOST.DEBUG);
    velements = getElements(input);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements, "\""), ",") << endl;
    }
    if (keep_pp == false) {
      return;
    }

    // copy info into vspecies and clear velements
    vector<string> vspecies;
    for (size_t i = 0; i < velements.size(); i++) {
      vspecies.push_back(velements[i]);
    }
    velements.clear();

    // parse string around these elements
    string::size_type loc1 = 0;
    string::size_type loc2 = string::npos;
    vector<string> vCAPITAL_LETTERS_PP;
    aurostd::string2tokens(CAPITAL_LETTERS_PP_LIST, vCAPITAL_LETTERS_PP, ",");
    bool found_CAPITAL_LETTERS_PP = false;
    bool found_CAPITAL_LETTERS = false;
    for (size_t i = 0; i < vspecies.size(); i++) {
      if ((i + 1) >= vspecies.size()) {
        loc2 = string::npos;
      } else {
        loc2 = input.find(vspecies[i + 1], loc1);
      }
      while (loc2 != string::npos) {
        found_CAPITAL_LETTERS_PP = false;
        for (size_t j = 0; j < vCAPITAL_LETTERS_PP.size() && found_CAPITAL_LETTERS_PP == false; j++) {
          if ((loc2 - (vCAPITAL_LETTERS_PP[j].size() - 1)) < input.size()) {
            continue;
          }
          found_CAPITAL_LETTERS = true;
          for (size_t k = 0; k < vCAPITAL_LETTERS_PP[j].size() && found_CAPITAL_LETTERS == true; k++) {
            if (input[loc2 - k] != vCAPITAL_LETTERS_PP[j][vCAPITAL_LETTERS_PP[j].size() - k - 1]) {
              found_CAPITAL_LETTERS = false;
            }
          }
          if (found_CAPITAL_LETTERS) {
            found_CAPITAL_LETTERS_PP = true;
          }
        }
        if (found_CAPITAL_LETTERS_PP == false) {
          break;
        }
        loc2 = input.find(vspecies[i], loc2 + 1);
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " loc1=" << loc1 << ", loc2=" << loc2 << endl;
      }
      if (loc2 == string::npos) { // CO20210315
        velements.push_back(input.substr(loc1));
        break;
      } else {
        velements.push_back(input.substr(loc1, loc2 - loc1));  // loc2-loc1 because it is the distance
        loc1 = loc2;
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements, "\""), ",") << endl;
    }
  }

  // ***************************************************************************
  // aurostd::getElements(string input,ostream&
  // oss,ofstream& FileMESSAGE)
  // ***************************************************************************
  // returns UNSORTED vector<string> from string
  vector<string> getElements(const string& input) { // CO20190712 //borrowed from XATOM_SplitAlloySpecies() //slow since we create many strings, but definitely works
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " original input=\"" << input << "\"" << endl;
    }
    string alloy = input;
    //[CO20190712 - no need for multiple passes anymore]for(uint i=1;i<=2;i++){alloy=KBIN::VASP_PseudoPotential_CleanName(alloy);} //be certain you clean everything, especially _GW (worst offender)
    aurostd::VASP_PseudoPotential_CleanName_InPlace(alloy); // be certain you clean everything, especially _GW (worst offender)
    aurostd::RemoveNumbersInPlace(alloy);              // remove composition
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " cleaned input=\"" << alloy << "\"" << endl;
    }
    vector<string> vspecies;
    for (uint i = 0; i < alloy.length(); i++) {
      if (alloy[i] >= 'A' && alloy[i] <= 'Z') {
        vspecies.emplace_back("");
      }
      vspecies.back() += alloy[i];
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vspecies pre ASCII clean=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies, "\""), ",") << endl;
    }
    for (size_t i = 0; i < vspecies.size(); i++) {
      aurostd::CleanStringASCII_InPlace(vspecies[i]);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " vspecies post ASCII clean=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies, "\""), ",") << endl;
    }
    return vspecies;
  }
  vector<string> getElements(const string& input, elements_string_type e_str_type, bool clean, bool sort_elements, bool keep_pp, ostream& oss) {  // overload
    ofstream FileMESSAGE;
    return getElements(input, e_str_type, FileMESSAGE, clean, sort_elements, keep_pp, oss);
  }
  // ME20190628 - added variant that also determines the composition
  template <class utype> vector<string> getElements(const string& input, vector<utype>& vcomposition, bool clean, bool sort_elements, bool keep_pp, ostream& oss) {
    ofstream FileMESSAGE;
    return getElements(input, vcomposition, composition_string, FileMESSAGE, clean, sort_elements, keep_pp, oss);  // this gets composition_string by default, pp_string has no composition
  }
  // cannot deduce utype from this construction
  vector<string> getElements(const string& input, elements_string_type e_str_type, ofstream& FileMESSAGE, bool clean, bool sort_elements, bool keep_pp, ostream& oss) {  // overload
    vector<double> vcomposition;
    return getElements(input, vcomposition, e_str_type, FileMESSAGE, clean, sort_elements, keep_pp, oss);
  }
  template <class utype> vector<string> getElements(const string& input, vector<utype>& vcomposition, elements_string_type e_str_type, bool clean, bool sort_elements, bool keep_pp, ostream& oss) { // overload
    ofstream FileMESSAGE;
    return getElements(input, vcomposition, e_str_type, FileMESSAGE, clean, sort_elements, keep_pp, oss);  // this gets composition_string by default, pp_string has no composition
  }
  template <class utype>
  vector<string> getElements(const string& _input, vector<utype>& vcomposition, elements_string_type e_str_type, ofstream& FileMESSAGE, bool clean, bool sort_elements, bool keep_pp, ostream& oss) { // main function
    const bool LDEBUG = (false || XHOST.DEBUG);
    vector<string> velements;
    vcomposition.clear(); // ME20190628

    //////////////////////////////////////////////////////////////////////////////
    // START Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " original input=" << _input << endl;
    }

    if (_input.empty()) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Empty input", FileMESSAGE, oss, _LOGGER_ERROR_);
      return velements;
    }

    string input = _input;

    if (clean && (e_str_type == composition_string || (e_str_type == pp_string && keep_pp == false))) {
      // in case we run into potpaw_PBE/Na, but only works for single elements, must be before check for isupper(input[0])
      const bool capital_letters_only = false; // default
      const bool remove_floats = false; // CO20210623 - explicitly KEEP floats for composition
      aurostd::VASP_PseudoPotential_CleanName_InPlace(input, capital_letters_only, remove_floats);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " checking input [1] =" << input << endl;
    }

    if (!isupper(input[0])) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Elements must be properly capitalized (input=" + input + ")", FileMESSAGE, oss, _LOGGER_ERROR_);
      return velements;
    }

    // we have a LIB1 problem... grab first everything before :
    // this is safe, as aflow generally introduces : in prototype, e.g., :LDAU2
    // this is safe anyway because elements would be BEFORE :
    if (clean) {
      // FAST
      string::size_type loc;
      //:
      loc = input.find(':');
      input = input.substr(0, loc);
      //_ICSD_
      loc = input.find("_ICSD_");
      input = input.substr(0, loc);
      // SLOW
      // vector<string> tokens;
      // aurostd::string2tokens(input,tokens,":");
      // input=tokens[0];
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " checking input [2] =" << input << endl;
    }

    // CO20180409 - running through input twice, no need, simply check at the end
    // uint numberOfElements = 0;
    // for (size_t i = 0; i < input.size(); i++) {
    //   if(isupper(input[i])) {
    //     numberOfElements++;
    //   }
    // }
    // if(numberOfElements == 0) {
    //   pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    //   return velements;
    // }

    //////////////////////////////////////////////////////////////////////////////
    // END Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Parsing input
    //////////////////////////////////////////////////////////////////////////////

    if (e_str_type == composition_string) {
      elementsFromCompositionString(input, velements, vcomposition);
    } else if (e_str_type == pp_string) {
      elementsFromPPString(input, velements, keep_pp);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Unknown compound designation", _INPUT_ILLEGAL_);
    }

    if (clean) {
      for (size_t i = 0; i < velements.size(); i++) {
        aurostd::CleanStringASCII_InPlace(velements[i]);
      } // CO20190712 - extra cleaning from XATOM_SplitAlloySpecies
    }

    // Add implicit ones
    for (size_t i = vcomposition.size(); i < velements.size(); i++) {
      vcomposition.push_back((utype) 1.0);
    }

    //////////////////////////////////////////////////////////////////////////////
    // END Parsing input
    //////////////////////////////////////////////////////////////////////////////

    if (velements.empty()) {
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, "No elements found", FileMESSAGE, oss, _LOGGER_ERROR_);
    }

    if (sort_elements && velements.size() > 1) {
      // this is MORE efficient that std::swap which has a copy constructor inside
      // http://www.cplusplus.com/reference/algorithm/swap/
      string etmp;
      utype ctmp = (utype) 0.0;
      for (size_t i = 0; i < velements.size() - 1; i++) {
        for (size_t j = i + 1; j < velements.size(); j++) {
          if (velements[i] > velements[j]) {
            etmp = velements[j]; // fix old j
            velements[j] = velements[i]; // swap
            velements[i] = etmp; // set i to old j
            ctmp = vcomposition[j]; // fix old j
            vcomposition[j] = vcomposition[i]; // swap
            vcomposition[i] = ctmp; // set i to old j
          }
        }
      }
    }

    return velements;
  }
#define AST_TEMPLATE(utype) template vector<string> getElements(const string&, vector<utype>&, bool, bool, bool, ostream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template vector<string> getElements(const string&, vector<utype>&, elements_string_type, bool, bool, bool, ostream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE
#define AST_TEMPLATE(utype) template vector<string> getElements(const string&, vector<utype>&, elements_string_type, ofstream&, bool, bool, bool, ostream&);
  AST_GEN_1(AST_UTYPE_NUM)
#undef AST_TEMPLATE

  // extractJsonKeysAflow//////////////////////////////////////////////////////
  // This function extracts keys from an aflowlib.json file. It is much
  // faster than using SQLite's JSON extension, but was designed to only
  // work for the aflowlib.json. It cannot handle nested JSONs!
  vector<string> extractJsonKeysAflow(const string& json) {
    vector<string> keys;
    string substring;
    string::size_type pos = 0;
    string::size_type lastPos = 0;
    string::size_type dpos = 0;
    string::size_type quote1 = 0;
    string::size_type quote2 = 0;
    string::size_type colon = 0;

    // Find the first comma - this is either the end of the key-value pair
    // or part of an array. Either way, the key is inside.
    pos = json.find(",");
    lastPos = 1; // First character is a curly brace, so skip
    dpos = pos - lastPos;
    while ((pos != string::npos) || (lastPos != string::npos)) {
      // A comma could be separating a key-value pair an array
      // or numbers or strings
      substring = json.substr(lastPos, dpos);

      // Find the colon - if there is no colon, it cannot be a key-value pair
      colon = substring.find(":");
      if (colon != string::npos) {
        // A key is enclosed in quotes, so there must be at least two of them
        quote1 = substring.find("\"");
        if (quote1 != string::npos) {
          quote2 = substring.find("\"", quote1 + 1);
          // Most non-keys are filtered out by now. There could still be array
          // elements left. In that case, however, the colon is between the quotes,
          // so make sure that the first two quotes appear before the colon and
          // take everything in-between as the key. This breaks if quotes, colons,
          // and commas are inside a string in the right sequence, but should not
          // be the case in AFLOW's JSON files.
          if ((quote2 != string::npos) && (quote1 < colon) && (quote2 < colon)) {
            substring = substring.substr(quote1 + 1, quote2 - quote1 - 1);
            if (!substring.empty()) {
              keys.push_back(substring);
            }
          }
        }
      }
      // Move on to the next comma
      lastPos = json.find_first_not_of(",", pos);
      pos = json.find(",", lastPos);
      dpos = pos - lastPos;
    }
    return keys;
  }

  // extractJsonValueAflow/////////////////////////////////////////////////////
  //  This function extracts values from an aflowlib.json file. It is much
  //  faster than using SQLite's JSON extension, but has was designed to only
  //  work for the aflowlib.json. It cannot handle nested JSONs!
  string extractJsonValueAflow(const string& json, string key) {
    string value;
    key = "\"" + key + "\":";
    string::size_type start = 0;
    string::size_type end = 0;
    start = json.find(key);
    if (start != string::npos) {
      start += key.length();
      end = json.find("\":", start);
      if (end != string::npos) {
        // In case there is any white space between key and value
        value = aurostd::RemoveWhiteSpacesFromTheFront(json.substr(start, end - start));
        // If we have a nested object, "value" should only be '{' + white space by now.
        if (value[0] == '{') {
          const string message = "JSON parser cannot read nested objects.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
        end = value.find_last_of(",");
        value = value.substr(0, end);
      } else {
        end = json.find("}", start);
        // In case there is any white space between key and value
        value = aurostd::RemoveWhiteSpacesFromTheFront(json.substr(start, end - start));
        // If we have a nested object, it should start with '{'
        if (value[0] == '{') {
          const string message = "JSON parser cannot read nested objects.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
        }
      }
    }
    return value;
  }

  // extractJsonVectorAflow/////////////////////////////////////////////////////
  vector<string> extractJsonVectorAflow(const string& json, string key) { // SD20220504
    const string value = extractJsonValueAflow(json, key);
    const string::size_type start = value.find("[");
    const string::size_type stop = value.rfind("]");
    vector<string> vec;
    aurostd::string2tokens(value.substr(start + 1, stop - 1), vec, ",");
    return vec;
  }

  // extractJsonMatrixAflow/////////////////////////////////////////////////////
  vector<vector<string>> extractJsonMatrixAflow(const string& json, string key) { // SD20220504
    const string value = extractJsonValueAflow(json, key);
    const string::size_type start = value.find("[[");
    const string::size_type stop = value.rfind("]]");
    vector<vector<string>> mat;
    vector<string> vec;
    vector<string> tokens;
    aurostd::string2tokensByDelimiter(value.substr(start + 2, stop - 2), tokens, "],[");
    for (size_t i = 0; i < tokens.size(); i++) {
      aurostd::string2tokens(tokens[i], vec, ",");
      mat.push_back(vec);
    }
    return mat;
  }

  /// @class x3DWriter
  /// @brief tools to create 3D scence and export them in different
  ///
  /// @authors
  /// @mod{HE,20220922,created x3Dwriter}
  /// @mod{HE,20240530,implemented in AFLOW4}
  ///
  /// Basic usage
  /// @code
  /// // load a structure to render
  /// const std::string structure_file = "example.vasp";
  /// xstructure work_structure(structure_file, IOAFLOW_AUTO);
  /// work_structure.ReScale(1.0);
  /// // prepare x3DWriter
  /// aurostd::x3DWriter w;
  /// w.scene_center = (work_structure.lattice.getcol(1) + work_structure.lattice.getcol(2) + work_structure.lattice.getcol(3)) / 3.0;
  /// w.addLatticeBox(work_structure.lattice, 0.1);
  /// for (auto &a: work_structure.atoms)
  ///   w.addSphere({a.cpos[1], a.cpos[2], a.cpos[3]}, 0.3, "bm_blue");
  /// aurostd::string2file(w.toHTML(), "/Users/nathan/Projects/AFLOW4/test_struc/AB_cP2_221_b_a.html");
  /// @endcode

  /// @brief reconstruct / clear
  void x3DWriter::clear() {
    *this = {};
  } // calling the constructor

  /// @brief create a copy of x3DWriter
  void x3DWriter::copy(const x3DWriter& x3w) {
    objects = x3w.objects;
    max_distance = x3w.max_distance;
    tachyon_zoom = x3w.tachyon_zoom;
    tachyon_camera_position = x3w.tachyon_camera_position;
    tachyon_camera_theta = x3w.tachyon_camera_theta;
    tachyon_camera_phi = x3w.tachyon_camera_phi;
    ani_type = x3w.ani_type;
    meta = x3w.meta;
    scene_center = x3w.scene_center;
    join_threshold = x3w.join_threshold;
  }

  /// @brief class constractor
  x3DWriter::x3DWriter() {
    meta.emplace_back("reference", "https://aflow.org");
    meta.emplace_back("generator", "AFLOW " + static_cast<std::string>(AFLOW_VERSION));
    meta.emplace_back("created", get_datetime_formatted("-", true, " ", ":"));

    // set up basic materials
    addMaterial("bm_black", {0.0, 0.0, 0.0});
    addMaterial("bm_white", {1.0, 1.0, 1.0});
    addMaterial("bm_red", {0.58, 0.07, 0.0});
    addMaterial("bm_blue", {0.01, 0.10, 0.58});
    addMaterial("bm_green", {0.0, 0.56, 0.0});
    addMaterial("bm_grey", {0.3, 0.3, 0.3});
    // set up grey glass
    Material matGreyGlass;
    matGreyGlass.name = "bm_grey_glass";
    matGreyGlass.color = {0.3, 0.3, 0.3};
    matGreyGlass.opacity = 0.3;
    matGreyGlass.specular = 0.0;
    matGreyGlass.ambient = 0.1;
    matGreyGlass.diffuse = 0.7;
    addMaterial(matGreyGlass);

    scene_center = {0.0, 0.0, 0.0};
  }

  /// @brief class copy constructor
  x3DWriter::x3DWriter(const x3DWriter& x3w) {
    copy(x3w);
  }

  /// @brief default class de-constructor
  x3DWriter::~x3DWriter() = default;

  /// @brief assignment operator
  x3DWriter& x3DWriter::operator=(const x3DWriter& x3w) {
    if (this == &x3w) {
      return *this;
    }
    copy(x3w);
    return *this;
  }

  /// @brief prepare a scene to show a lattice by setting the center and a view collection
  void x3DWriter::prepareSceneLattice(const xmatrix<double>& lattice) {
    const xvector<double> endpoint = lattice(1) + lattice(2) + lattice(3);
    scene_center = endpoint / 2.0;
    constexpr std::array<std::pair<int, int>, 3> view_set({
        {{1, 2}, {2, 3}, {3, 1}}
    });
    tachyon_lattice_views.clear();
    for (const auto& [main, secondary] : view_set) {
      const xvector<double> plane_normal = aurostd::vector_product(lattice(main), lattice(secondary));
      const xvector<double> up_dir = std::cos(pi * 0.5) * lattice(main) + std::sin(pi * 0.5) * plane_normal;
      tachyon_lattice_views.emplace_back(lattice(main) * 2, aurostd::normalizeSumToOne(up_dir));
    }
  }

  /// @brief add a Sphere to the scene
  void x3DWriter::addSphere(const xvector<double>& center, double radius, const std::string& material) {
    const std::shared_ptr<x3DWriter::Sphere> newSphere = std::make_shared<x3DWriter::Sphere>();
    newSphere->center = center - scene_center;
    newSphere->radius = radius;
    newSphere->material = material;
    max_distance = max(max_distance, aurostd::modulus(center));
    x3DWriter::storage_object so;
    so.obj = newSphere; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::SPHERE;
    objects.emplace_back(so);
  }

  /// @brief add a set of Sphere with a list of center and constant radius and materials
  void x3DWriter::addSpheres(const vector<xvector<double>>& center, const double radius, const std::string& material) {
    for (const xvector<double>& c : center) {
      addSphere(c, radius, material);
    }
  }

  /// @brief add a box forming a unit cell to the scene
  /// @param lattice 3x3 matrix like xstructure.lattice
  /// @param radius radius of the cylinders forming the box
  /// @param material box material
  void x3DWriter::addLatticeBox(const xmatrix<double>& lattice, double radius, const std::string& material, bool axis) {
    vector<xvector<double>> points;
    points.push_back({0, 0, 0});
    for (int i = lattice.lrows; i <= lattice.urows; i++) {
      points.push_back(lattice(i));
    }
    points.push_back(points[1] + points[2]);
    points.push_back(points[1] + points[3]);
    points.push_back(points[2] + points[3]);
    points.push_back(points[1] + points[2] + points[3]);
    addSpheres(points, radius, material);
    vector<std::pair<uint, uint>> connections = {
        {1, 4},
        {1, 5},
        {2, 4},
        {2, 6},
        {3, 5},
        {3, 6},
        {7, 4},
        {7, 5},
        {7, 6}
    };
    ;
    if (!axis) {
      connections.emplace_back(std::pair<uint, uint>({0, 1}));
      connections.emplace_back(std::pair<uint, uint>({0, 2}));
      connections.emplace_back(std::pair<uint, uint>({0, 3}));
    }

    for (auto [base, apex] : connections) {
      addOpenCylinder(points[base], points[apex], radius, material);
    }
    if (axis) {
      addOpenCylinder(points[0], points[1], radius, "bm_red");
      addOpenCylinder(points[0], points[2], radius, "bm_green");
      addOpenCylinder(points[0], points[3], radius, "bm_blue");
    }
  }

  /// @brief add a OpenCylinder to the scene
  /// @param base start point
  /// @param apex end point
  /// @param radius cylinder radius
  /// @param material cylinder material
  void x3DWriter::addOpenCylinder(const xvector<double>& base, const xvector<double>& apex, double radius, const std::string& material) {
    const std::shared_ptr<x3DWriter::OpenCylinder> newCylinder = std::make_shared<x3DWriter::OpenCylinder>();
    const xvector<double> shift_base = base - scene_center;
    const xvector<double> shift_apex = apex - scene_center;
    const xvector<double> axis = shift_base - shift_apex;
    newCylinder->base = shift_base;
    newCylinder->apex = shift_apex;
    newCylinder->center = shift_apex + (axis / 2);
    newCylinder->radius = radius;
    newCylinder->material = material;
    newCylinder->height = aurostd::modulus(axis);
    newCylinder->theta = pi / 2 + std::acos(axis(1) / newCylinder->height);
    newCylinder->phi = std::atan2(axis(3), axis(2));
    max_distance = max(max_distance, aurostd::modulus(shift_base));
    max_distance = max(max_distance, aurostd::modulus(shift_apex));
    x3DWriter::storage_object so;
    so.obj = newCylinder; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::OPEN_CYLINDER;
    objects.emplace_back(so);
  }

  /// @brief add ConvexFacets to the Scene
  /// @param vertexes facet corners
  /// @param facets vertex index that form the facets
  /// @param material facet materials
  /// @param shift
  void x3DWriter::addConvexFacets(const vector<xvector<double>>& vertexes, const vector<vector<uint>>& facets, const std::string& material, const xvector<double>& shift) {
    const std::shared_ptr<x3DWriter::ConvexFacets> newFacet = std::make_shared<x3DWriter::ConvexFacets>();
    newFacet->material = material;

    vector<xvector<double>> shifted_vertexes;
    shifted_vertexes.reserve(vertexes.size());
    for (const xvector<double>& vertex : vertexes) {
      shifted_vertexes.push_back(vertex + shift - scene_center);
    }
    newFacet->vertexes = shifted_vertexes;

    newFacet->facets = facets;
    for (const auto& vertex : vertexes) {
      max_distance = max(max_distance, aurostd::modulus(vertex));
    }
    x3DWriter::storage_object so;
    so.obj = newFacet; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::FACET;
    objects.emplace_back(so);
  }

  /// @brief add ConvexFacets with a grey glass material to the Scene
  /// @param vertexes facet corners
  /// @param facets vertex index that form the facets
  /// @param shift
  void x3DWriter::addConvexFacets(const vector<xvector<double>>& vertexes, const vector<vector<uint>>& facets, const xvector<double>& shift) {
    addConvexFacets(vertexes, facets, "bm_grey_glass", shift);
  }

  /// @brief add a new basic Material
  /// @param name name to reference the Material
  /// @param color as R,G,B xvector
  void x3DWriter::addMaterial(const std::string& name, const xvector<double>& color) {
    Material newMat;
    newMat.name = name;
    newMat.color = color;
    addMaterial(newMat);
  }

  /// @brief add a new Material
  /// @param newMaterial
  void x3DWriter::addMaterial(const Material& newMaterial) {
    const std::shared_ptr<x3DWriter::Material> newMat = std::make_shared<x3DWriter::Material>(newMaterial);
    x3DWriter::storage_object so;
    so.obj = newMat; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::MATERIAL;
    objects.emplace_back(so);
  }

  /// @brief create a set of materials using the turbo colormap
  /// @param count number of materials to create
  ///
  /// the materials will be named "auto_cs_{idx}"
  vector<std::string> x3DWriter::addColorSpreadMaterial(const uint count) {
    vector<std::string> created_materials;
    const uint color_step = 256 / count;
    for (uint color_index = 0; color_index < count; color_index++) {
      const xvector<double> color;
      color(1) = turbo_cmap[color_index * color_step][0];
      color(2) = turbo_cmap[color_index * color_step][1];
      color(3) = turbo_cmap[color_index * color_step][2];
      created_materials.emplace_back("auto_cs_" + std::to_string(color_index));
      addMaterial(created_materials[color_index], color);
    }
    return created_materials;
  }
  template <typename Container> std::map<std::string, std::string> x3DWriter::addNamedColorSpreadMaterial(const Container& names) {
    std::vector<string> color_spread = addColorSpreadMaterial(names.size());
    std::map<std::string, std::string> mapped_colors;
    std::transform(names.begin(), names.end(), color_spread.begin(), std::inserter(mapped_colors, mapped_colors.end()), std::make_pair<const string&, const string&>);
    return mapped_colors;
  }
  template std::map<std::string, std::string> x3DWriter::addNamedColorSpreadMaterial(const std::vector<std::string>&);
  template std::map<std::string, std::string> x3DWriter::addNamedColorSpreadMaterial(const std::deque<std::string>&);

  /// @brief combine overlapping vertexes and facets
  /// @param vertexes originial vertexes positions
  /// @param facets list of vertex indexes forming the facets
  /// @param new_vertexes updated vertexes positions
  /// @param new_facets updated list of vertex indexes forming the facets
  void x3DWriter::joinFacets(vector<xvector<double>>& vertexes, vector<vector<uint>>& facets, const vector<xvector<double>>& new_vertexes, const vector<vector<uint>>& new_facets) const {
    if (vertexes.empty() || facets.empty()) {
      vertexes = new_vertexes;
      facets = new_facets;
      return;
    }
    bool is_set = false;

    vector<uint> vertex_idx_map;
    const size_t original_vertex_count = vertexes.size();
    vector<uint> prepared_facet;
    for (const auto& new_vertexe : new_vertexes) {
      is_set = false;
      for (uint ov_id = 0; ov_id < original_vertex_count; ov_id++) {
        if (distance(vertexes[ov_id], new_vertexe) <= join_threshold) {
          vertex_idx_map.emplace_back(ov_id);
          is_set = true;
          break;
        }
      }
      if (is_set) {
        continue;
      }
      vertex_idx_map.emplace_back(vertexes.size());
      vertexes.emplace_back(new_vertexe);
    }

    for (const auto& new_facet : new_facets) {
      prepared_facet.clear();
      for (const uint v_id : new_facet) {
        prepared_facet.emplace_back(vertex_idx_map[v_id]);
      }
      facets.emplace_back(prepared_facet);
    }

    std::set<uint> facets_to_remove;
    const size_t facet_count = facets.size();
    bool overlap = true;
    // remove double facets to avoid z-fighting
    for (size_t row = 0; row < facet_count; row++) {
      for (size_t col = row + 1; col < facet_count; col++) {
        if (facets_to_remove.find(row) != facets_to_remove.end() || facets_to_remove.find(col) != facets_to_remove.end()) {
          continue;
        }
        overlap = true;
        if (facets[row].size() < facets[col].size()) {
          for (size_t vert_idx = 0; vert_idx < facets[row].size(); vert_idx++) {
            overlap = std::find(facets[col].begin(), facets[col].end(), facets[row][vert_idx]) != facets[col].end();
            if (!overlap) {
              break;
            }
          }
          if (overlap) {
            facets_to_remove.insert(row);
            break;
          }
        } else {
          for (size_t vert_idx = 0; vert_idx < facets[col].size(); vert_idx++) {
            overlap = std::find(facets[row].begin(), facets[row].end(), facets[col][vert_idx]) != facets[row].end();
            if (!overlap) {
              break;
            }
          }
          if (overlap) {
            facets_to_remove.insert(col);
          }
        }
      }
    }
    vector<vector<uint>> cleaned_facets;
    for (uint of_idx = 0; of_idx < facets.size(); of_idx++) {
      if (facets_to_remove.find(of_idx) == facets_to_remove.end()) {
        cleaned_facets.emplace_back(facets[of_idx]);
      }
    }
    facets = cleaned_facets;
  }

  /// @brief prepare an animation
  /// @param duration animation time in second
  /// @param out_folder folder to store the generated files
  /// @param fps frames per second
  /// @param lr move left to right
  ///
  /// In the out_folder a render.sh script that will create a file with the stem _Animation;
  /// for rendering [ffmpeg](https://ffmpeg.org/) and [tachyon](http://jedi.ks.uiuc.edu/~johns/tachyon/) is needed
  void x3DWriter::animate(float duration, const std::filesystem::path& out_folder, uint fps, bool lr) {
    std::filesystem::path file_path;
    std::error_code err;
    if (!std::filesystem::create_directories(out_folder, err)) {
      if (err.value() != 0) // exist is okay
      {
        const string message = "Can not create result folder " + out_folder.string() + ". Error: " + err.message();
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
    }
    const uint frame_count = duration * fps;
    std::stringstream command_content;
    const double angle_step = 360.0 / frame_count;
    for (uint frame_index = 0; frame_index < frame_count; frame_index++) {
      std::stringstream filename;
      filename << "frame_" << std::setw(6) << std::setfill('0') << frame_index;
      file_path = out_folder / filename.str();
      if (lr) {
        tachyon_camera_phi = frame_index * angle_step;
      } else {
        tachyon_camera_theta = frame_index * angle_step;
      }
      aurostd::string2file(toTachyon(), file_path.string() + ".dat");
      command_content << "tachyon " << file_path.filename() << ".dat -o " << file_path.filename() << ".ppm -format PPM -auto_skylight 0.7 -fullshade -numthreads 8" << endl;
    }

    switch (ani_type) {
      case (animation_format::MP4): {
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -vcodec libx264 -crf 20 -pix_fmt yuv420p  -an _Animation.mp4" << endl;
        break;
      }
      case (animation_format::GIF): {
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -filter_complex [0:v] fps=" << fps << ", split [a][b];[a] palettegen [p];[b][p] paletteuse _Animation.gif" << endl;
        break;
      }
      case (animation_format::WEBM): {
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -c:v libvpx-vp9 -b:v 0 -crf 20 -pix_fmt yuv420p -pass 1 -an -f null /dev/null" << endl;
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -c:v libvpx-vp9 -b:v 0 -crf 20 -pix_fmt yuv420p -pass 2 -an _Animation.webm" << endl;
        break;
      }
    }
    aurostd::string2file(command_content.str(), (out_folder / "render.sh"));
  }

  /// convert a material to x3d
  /// @return x3d version of a material (xml)
  std::string x3DWriter::x3d_material(const std::shared_ptr<x3DWriter::Material>& material) {
    std::stringstream mat;
    mat << "<Material";
    mat << " diffuseColor='" << material->color << "'";
    mat << " ambientIntensity='" << material->ambient << "'";
    mat << " shininess='" << material->specular << "'";
    mat << " transparency='" << 1.0 - material->opacity << "'/>";
    return mat.str();
  }

  /// @brief calculate the camera position (used fot tachyon)
  void x3DWriter::tachyon_calculate_camera() {
    if (tachyon_camera_orthographic) {
      tachyon_zoom = 1.0 / (max_distance * 1.05);
    } else {
      tachyon_zoom = 2.0;
    }

    const double camera_distance = std::tan((180.0 - (tachyon_camera_angle / tachyon_zoom)) / (180 * 2.0) * pi) * max_distance;
    tachyon_camera_position(1) = camera_distance * std::sin(tachyon_camera_phi / 180.0 * pi) * std::cos(tachyon_camera_theta / 180.0 * pi);
    tachyon_camera_position(2) = camera_distance * std::cos(tachyon_camera_phi / 180.0 * pi) * std::sin(tachyon_camera_theta / 180.0 * pi);
    tachyon_camera_position(3) = camera_distance * std::cos(tachyon_camera_phi / 180.0 * pi);
  }

  /// @brief save scene for the tachyon renderer
  /// use [tachyon](http://jedi.ks.uiuc.edu/~johns/tachyon/) to render the generated file
  std::string x3DWriter::toTachyon() {
    tachyon_calculate_camera();
    stringstream content;
    for (const auto& [name, meta_content] : meta) {
      content << "# " << name << ": " << meta_content << endl;
    }
    content << endl;
    content << "BEGIN_SCENE" << endl;
    content << "  RESOLUTION 720 720" << endl << endl;
    content << "CAMERA" << endl;
    if (tachyon_camera_orthographic) {
      content << "  PROJECTION ORTHOGRAPHIC" << endl;
    }
    content << "  ZOOM " << tachyon_zoom << "  ASPECTRATIO -1.0 ANTIALIASING 2 RAYDEPTH 12" << endl;
    if (tachyon_lattice_views_idx == -1) {
      content << "  CENTER " << tachyon_camera_position << endl;
      content << "  VIEWDIR " << -tachyon_camera_position / max_distance << endl;
      content << "  UPDIR 0 1 0" << endl;
    } else {
      const auto& [center, up] = tachyon_lattice_views[tachyon_lattice_views_idx];
      content << "CENTER " << center << endl;
      content << "VIEWDIR " << -aurostd::normalizeSumToOne(center) << endl;
      content << "UPDIR " << up << endl << endl;
    }
    content << "END_CAMERA" << endl << endl;
    content << "BACKGROUND 1.0 1.0 1.0" << endl;

    for (const auto& [type, obj] : objects) {
      if (type == x3DWriter::object_types::MATERIAL) {
        const std::shared_ptr<x3DWriter::Material> material = std::static_pointer_cast<x3DWriter::Material>(obj);
        content << "TEXDEF " << material->name << endl;
        content << "  AMBIENT " << material->ambient;
        content << "  DIFFUSE " << material->diffuse;
        content << "  SPECULAR " << material->specular;
        content << "  OPACITY " << material->opacity << endl;
        content << "  COLOR " << material->color;
        content << "  TEXFUNC 0" << endl;
      }
    }

    for (const auto& [type, obj] : objects) {
      switch (type) {
        case object_types::MATERIAL: {
          // written to the start of the file
          break;
        }
        case object_types::SPHERE: {
          const std::shared_ptr<x3DWriter::Sphere> sphere = std::static_pointer_cast<x3DWriter::Sphere>(obj);
          content << "SPHERE" << endl;
          content << "  CENTER " << sphere->center << endl;
          content << "  RAD " << sphere->radius << endl;
          content << "  " << sphere->material << endl << endl;
          break;
        }
        case object_types::FACET: {
          const std::shared_ptr<x3DWriter::ConvexFacets> facet_container = std::static_pointer_cast<x3DWriter::ConvexFacets>(obj);
          for (auto facet : facet_container->facets) {
            for (size_t tri_index = 1; tri_index < facet.size() - 1; tri_index++) {
              content << "TRI" << endl;
              content << "  V0 " << facet_container->vertexes[facet[0]] << endl;
              content << "  V1 " << facet_container->vertexes[facet[tri_index]] << endl;
              content << "  V2 " << facet_container->vertexes[facet[tri_index + 1]] << endl;
              content << "  " << facet_container->material << endl << endl;
            }
          }
          break;
        }
        case object_types::OPEN_CYLINDER: {
          const std::shared_ptr<x3DWriter::OpenCylinder> cylinder = std::static_pointer_cast<x3DWriter::OpenCylinder>(obj);
          content << "FCylinder" << endl;
          content << "  BASE " << cylinder->base << endl;
          content << "  APEX " << cylinder->apex << endl;
          content << "  RAD " << cylinder->radius << endl;
          content << "  " << cylinder->material << endl << endl;
          break;
        }
      }
    }

    content << "END_SCENE" << endl;

    return content.str();
  }

  /// @brief save scene as x3d (xml)
  /// @param include_xml create as stand alone xml
  /// @param replace_material switch if Materials are exported
  std::string x3DWriter::toX3D(const bool include_xml, const bool replace_material) {
    stringstream x3d_content;
    if (include_xml) {
      x3d_content << R"(<X3D profile="Immersive" version="3.0">)" << endl << endl;
      x3d_content << "<head>" << endl;
      for (const auto& [key, value] : meta) {
        x3d_content << "  <meta name=\"" << key << "\" content=\"" << value << "\"/>" << endl;
      }
      x3d_content << "</head>" << endl << endl;
    }

    x3d_content << "<Scene>" << endl;
    std::map<std::string, std::string> material_lookup;

    for (const auto& [type, obj] : objects) {
      if (type == x3DWriter::object_types::MATERIAL) {
        const std::shared_ptr<x3DWriter::Material> material = std::static_pointer_cast<x3DWriter::Material>(obj);
        material_lookup[material->name] = x3d_material(material);
      }
    }
    if (!replace_material) {
      for (const auto& [name, content] : material_lookup) {
        x3d_content << "<Appearance DEF='" << name << "'>";
        x3d_content << content;
        x3d_content << "</Appearance>" << endl;
      }
    }

    for (const auto& [type, obj] : objects) {
      switch (type) {
        case object_types::MATERIAL: {
          // written to the start of the file
          break;
        }
        case object_types::SPHERE: {
          const std::shared_ptr<x3DWriter::Sphere> sphere = std::static_pointer_cast<x3DWriter::Sphere>(obj);
          x3d_content << "<Transform translation='" << sphere->center << "'><Shape>" << endl;
          if (replace_material) {
            x3d_content << "  <Appearance>" << material_lookup[sphere->material] << "</Appearance>" << endl;
          } else {
            x3d_content << "  <Appearance USE='" << sphere->material << "'/>" << endl;
          }
          x3d_content << "  <Sphere radius='" << sphere->radius << "'/>" << endl;
          x3d_content << "</Shape></Transform>" << endl << endl;
          break;
        }
        case object_types::FACET: {
          const std::shared_ptr<x3DWriter::ConvexFacets> facet_container = std::static_pointer_cast<x3DWriter::ConvexFacets>(obj);
          x3d_content << "<Shape>" << endl;
          if (replace_material) {
            x3d_content << "  <Appearance>" << material_lookup[facet_container->material] << "</Appearance>" << endl;
          } else {
            x3d_content << "  <Appearance USE='" << facet_container->material << "'/>" << endl;
          }
          x3d_content << R"(  <IndexedFaceSet solid="false" colorPerVertex="false" normalPerVertex="false" coordIndex=)" << endl;
          for (const auto& facet : facet_container->facets) {
            x3d_content << "    ";
            for (const auto point_index : facet) {
              x3d_content << point_index << " ";
            }
            x3d_content << "-1" << endl; // end a facet
          }
          x3d_content << "  \">" << endl;
          x3d_content << "  <Coordinate point=\"" << endl;
          for (const auto& coordinate : facet_container->vertexes) {
            x3d_content << "    " << coordinate << endl;
          }
          x3d_content << "  \"/>" << endl;
          x3d_content << "  </IndexedFaceSet>" << endl;
          x3d_content << "</Shape>" << endl << endl;
          break;
        }
        case object_types::OPEN_CYLINDER: {
          const std::shared_ptr<x3DWriter::OpenCylinder> cylinder = std::static_pointer_cast<x3DWriter::OpenCylinder>(obj);
          x3d_content << "<Transform translation='" << cylinder->center << "'>" << endl;
          x3d_content << "<Transform rotation='1 0 0 " << cylinder->phi << "'><Transform rotation='0 0 1 " << cylinder->theta << "'><Shape>" << endl;
          if (replace_material) {
            x3d_content << "  <Appearance>" << material_lookup[cylinder->material] << "</Appearance>" << endl;
          } else {
            x3d_content << "  <Appearance USE='" << cylinder->material << "'/>" << endl;
          }
          x3d_content << "  <Cylinder radius='" << cylinder->radius << "' height='" << cylinder->height << "' top='FALSE' bottom='FALSE'/>" << endl;
          x3d_content << "</Shape></Transform></Transform></Transform>" << endl << endl;
          break;
        }
      }
    }

    x3d_content << "</Scene>" << endl;
    if (include_xml) {
      x3d_content << "</X3D>" << endl;
    }

    return x3d_content.str();
  }

  /// @brief save scene as x3d and embed it in a html file
  std::string x3DWriter::toHTML() {
    stringstream content;

    content << "<!DOCTYPE html>\n<html lang=\"en\">" << endl << endl;
    content << "<head>" << endl;
    content << "<meta charset=\"UTF-8\">" << endl;
    content << R"(<meta name="viewport" content="width=device-width, initial-scale=1.0">)" << endl;
    for (const auto& [name, meta_content] : meta) {
      content << "  <meta name=\"" << name << "\" content=\"" << meta_content << "\"/>" << endl;
    }
    content << "<title>AFLOW X3D</title>" << endl;
    content << "<script src='https://www.x3dom.org/download/x3dom.js'> </script>" << endl;
    // CDN makes it directly usable
    content << "<link rel='stylesheet' type='text/css' href='https://www.x3dom.org/download/x3dom.css'/>" << endl;
    content << "</head>" << endl << endl;
    content << "<body>" << endl;
    content << "<h1>AFLOW 3D viewer</h1>" << endl << endl;

    content << "<x3d width='1200px' height='1200px'>" << endl;
    content << toX3D(false, true) << endl << endl;
    content << "</x3d>" << endl;
    content << "</body>\n</html>" << endl;

    return content.str();
  }

} // namespace aurostd

#endif // _AUROSTD_XPARSER_CPP_
