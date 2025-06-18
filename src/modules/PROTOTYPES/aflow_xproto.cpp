// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2008-2015
// this is one of the most complicate parts of aflow
// fix Ag1Te3_ICSD_37186

/// home/auro/work/AFLOW3/aflow --potential=potpaw_GGA --ldau2 --potential_complete --aflow_proto=ICSD_52482.A,ICSD_52483.A,ICSD_76031.A,ICSD_652630.A,ICSD_652632.A,ICSD_652633.A,ICSD_652635.A,ICSD_652637.A,ICSD_652640.A,ICSD_108746.A,ICSD_246657.A:Sm_2,Sm_3

// /home/auro/work/AFLOW3/aflow --potential=potpaw_GGA  --potential_complete --aflow_proto=A3,ICSD_94429.A,ICSD_108026.A,ICSD_165132.A,ICSD_240995.A,ICSD_22300,A,ICSD_240995.A,ICSD_14288.A,ICSD_43431.A:B,B_s,B_h

#ifndef _AFLOW_XPROTO_CPP
#define _AFLOW_XPROTO_CPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <iterator>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xrandom.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_xhost.h"
#include "modules/COMPARE/aflow_compare_structure.h" //DX20181009
#include "modules/PROTOTYPES/aflow_anrl.h"
#include "modules/SYM/aflow_symmetry_spacegroup.h" //DX20181009
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_xatom.h"

#define _EPS_ 0.02

using std::deque;
using std::ostream;
using std::vector;

string *LOAD_Library_ICSD(string file);

namespace aurostd {
  template <class utype> void swap(utype &a, utype &b) {
    utype temp = a;
    a = b;
    b = temp;
  }
} // namespace aurostd

// ***************************************************************************
// Function extra operator << for vector
template <class utype> // operator <<  vector<>
std::ostream &operator<<(std::ostream &buf, const std::vector<utype> &x) {
  for (size_t i = 0; i < x.size(); i++) {
    buf << x[i] << " ";
  }
  return buf;
}
// ***************************************************************************
// Function extra operator << for deque
template <class utype> // operator <<  deque<>
std::ostream &operator<<(std::ostream &buf, const std::deque<utype> &x) {
  for (size_t i = 0; i < x.size(); i++) {
    buf << x[i] << " ";
  }
  return buf;
}

// ***************************************************************************
// ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL
// ***************************************************************************

// ***************************************************************************
namespace aflowlib {
  string PrototypeCleanLatticeString(const string &latticeIN) {
    string lattice = latticeIN;
    lattice = aurostd::RemoveSubStringFirst(lattice, "fcc");
    lattice = aurostd::RemoveSubStringFirst(lattice, "FCC");
    lattice = aurostd::RemoveSubStringFirst(lattice, "bcc");
    lattice = aurostd::RemoveSubStringFirst(lattice, "BCC");
    lattice = aurostd::RemoveSubStringFirst(lattice, "hcp");
    lattice = aurostd::RemoveSubStringFirst(lattice, "HCP");
    lattice = aurostd::RemoveSubStringFirst(lattice, "sc");
    lattice = aurostd::RemoveSubStringFirst(lattice, "SC");
    lattice = aurostd::RemoveSubStringFirst(lattice, "f");
    lattice = aurostd::RemoveSubStringFirst(lattice, "F");
    lattice = aurostd::RemoveSubStringFirst(lattice, "b");
    lattice = aurostd::RemoveSubStringFirst(lattice, "B");
    lattice = aurostd::RemoveSubStringFirst(lattice, "h");
    lattice = aurostd::RemoveSubStringFirst(lattice, "H");
    lattice = aurostd::RemoveSubStringFirst(lattice, "c");
    lattice = aurostd::RemoveSubStringFirst(lattice, "C");
    lattice = aurostd::RemoveSubStringFirst(lattice, "p");
    lattice = aurostd::RemoveSubStringFirst(lattice, "P");
    lattice = aurostd::RemoveSubStringFirst(lattice, "s");
    lattice = aurostd::RemoveSubStringFirst(lattice, "S");
    lattice = aurostd::RemoveSubStringFirst(lattice, " ");
    return lattice;
  }
} // namespace aflowlib

// ***************************************************************************
// ***************************************************************************
namespace aflowlib {
  std::set<string> GetPrototypesByStoichiometry(const vector<uint> &stoichiometry, const string &library) {
    // return uid list
    std::vector<uint> stoichiometry_reduced;
    aurostd::reduceByGCD(stoichiometry, stoichiometry_reduced);
    std::sort(stoichiometry_reduced.begin(), stoichiometry_reduced.end(), std::greater<>());
    const std::string stoichiometry_key = aurostd::joinWDelimiter(stoichiometry_reduced, ":");
    const anrl::ProtoData pd = anrl::ProtoData::get();
    if (pd.lookup["stoichiometry"].find(stoichiometry_key) == pd.lookup["stoichiometry"].end()) {
      return {};
    }
    return pd.lookup["stoichiometry"][stoichiometry_key];
  }
} // namespace aflowlib

// ***************************************************************************
// ***************************************************************************
namespace aflowlib {
  std::set<string> GetPrototypesBySymmetry(const vector<uint> &stoichiometry, const uint space_group_number, const vector<GroupedWyckoffPosition> &grouped_Wyckoff_positions, const uint setting, const string &library) {
    const anrl::ProtoData pd = anrl::ProtoData::get();
    const std::set<std::string> vuid_stoichiometry = GetPrototypesByStoichiometry(stoichiometry, library);
    if (pd.lookup["space_group_number"].find(std::to_string(space_group_number)) == pd.lookup["space_group_number"].end()) {
       return {};
    }
    const std::set<std::string> vuid_space_group = pd.lookup["space_group_number"][std::to_string(space_group_number)];
    std::vector<std::string> vuid_symmetry_temp;
    std::set_intersection(vuid_stoichiometry.begin(), vuid_stoichiometry.end(), vuid_space_group.begin(), vuid_space_group.end(), std::back_inserter(vuid_symmetry_temp));
    std::set<std::string> vuid_symmetry(vuid_symmetry_temp.begin(), vuid_symmetry_temp.end());
    if (grouped_Wyckoff_positions.empty()) {
      vector<GroupedWyckoffPosition> prototype_grouped_Wyckoffs;
      for (const auto &uid : vuid_symmetry) {
        const vector<vector<string>> prototype_grouped_Wyckoff_letters = compare::convertANRLWyckoffString2GroupedPositions(static_cast<string>(pd.content[uid]["label"]));
        compare::groupWyckoffPositionsFromGroupedString(static_cast<uint>(pd.content[uid]["space_group_number"]), setting, prototype_grouped_Wyckoff_letters, prototype_grouped_Wyckoffs);
        if (!compare::matchableWyckoffPositions(grouped_Wyckoff_positions, prototype_grouped_Wyckoffs, false)) {
          // remove uid if not Wyckoff match
          vuid_symmetry.erase(uid);
        }
      }
    }

    return vuid_symmetry;
  }
} // namespace aflowlib

// ***************************************************************************
// ***************************************************************************
namespace aflowlib {
  /// @brief return the number of species for a prototype
  /// @param _label user provided label to check
  /// @note searches first in the database and then assumes standard label format
  /// @return number of species
  /// @authors
  /// @mod{CO,20181226,created}
  /// @mod{HE,20250519,updated to new prototype DB}
  uint PrototypeLibrariesSpeciesNumber(const string &_label) { // CO20181226
    if (_label.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "label empty [1]", _VALUE_ILLEGAL_);
    }
    // in avasp, input here can be comma-separated series of labels
    // assume user intelligent and that all labels provided are same nspecies count
    // take first one
    vector<string> labeltokens;
    aurostd::string2tokens(_label, labeltokens, ",");
    if (labeltokens.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "label empty [2]", _VALUE_ILLEGAL_);
    }
    const string label = labeltokens[0];

    // check the database
    const anrl::ProtoData pd = anrl::ProtoData::get();
    const std::string uid = anrl::getPrototypeUID(label);
    if (!uid.empty()) {
      return static_cast<uint>(pd.content[uid]["number_of_species"]);
    }

    // check ANRL
    vector<string> tokens;
    aurostd::string2tokens(label, tokens, "_");
    if (tokens.size() >= 4) {
      return static_cast<uint>(tokens.size()) - 3;
    }

    // otherwise binary
    return (uint) 2;
  }
} // namespace aflowlib

// ***************************************************************************
// ***************************************************************************
namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss, const string &label, const string &parameters, int mode) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,int mode)" << endl;
    }
    deque<string> atomX;
    deque<double> volumeX;
    const uint nspeciesHTQC = aflowlib::PrototypeLibrariesSpeciesNumber(label);
    for (uint i = 0; i < nspeciesHTQC; i++) {
      string string_tmp = "A";
      string_tmp[0] += i;
      atomX.push_back(string_tmp);
      volumeX.push_back(1.0);
    }
    return aflowlib::PrototypeLibraries(oss, label, parameters, atomX, volumeX, -2.0, mode, false);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss, const string &label, const string &parameters, deque<string> &atomX, int mode) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX)" << endl;
    }
    deque<double> volumeX;
    for (size_t i = 0; i < atomX.size(); i++) {
      if (LDEBUG) { // CO20181106
        cerr << __AFLOW_FUNC__ << " atomX.at(" << i << ")=" << atomX[i] << endl;
        cerr << __AFLOW_FUNC__ << " GetAtomVolume()=" << GetAtomVolume(atomX[i]) << endl;
      }
      volumeX.push_back(GetAtomVolume(atomX[i]));
      //[CO20181106]volumeX.push_back(GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomX[i])));
    }
    return aflowlib::PrototypeLibraries(oss, label, parameters, atomX, volumeX, -3.0, mode, false);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss, const string &label, const string &parameters, deque<string> &atomX, deque<double> &volumeX, int mode) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,deque<double> &volumeX)" << endl;
    }
    return aflowlib::PrototypeLibraries(oss, label, parameters, atomX, volumeX, -4.0, mode, false);
  }
} // namespace aflowlib

namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss, const string &label, const string &parameters, deque<string> &atomX, deque<double> &volumeX, double volume_in, int mode) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &atomX,deque<double> &volumeX,double volume_in)" << endl;
    }
    return aflowlib::PrototypeLibraries(oss, label, parameters, atomX, volumeX, volume_in, mode, false);
  }
} // namespace aflowlib

// double NearestNeighbor(const xstructure &str_in);

// ***************************************************************************
// the mother of all the prototypes
namespace aflowlib {
  xstructure PrototypeLibraries(ostream &oss, string label, string parameters, deque<string> &vatomX, deque<double> &vvolumeX, double volume_in, int mode, bool flip_option) { // COMPLETE ONE
    // XHOST.DEBUG=true;
    const bool LDEBUG = (false || XHOST.DEBUG);

    const std::string uid = anrl::getPrototypeUID(label);
    const anrl::ProtoData pd = anrl::ProtoData::get();

    stringstream message;
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option)" << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: label=" << label << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: parameters=" << parameters << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: vatomX.size()=" << vatomX.size() << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: vatomX =";
      for (size_t i = 0; i < vatomX.size(); i++) {
        cerr << " " << vatomX[i];
      }
      cerr << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: vvolumeX.size()=" << vvolumeX.size() << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: vvolumeX =";
      for (size_t i = 0; i < vvolumeX.size(); i++) {
        cerr << " " << vvolumeX[i];
      }
      cerr << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: volume_in=" << volume_in << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: mode=" << mode << endl;
    }
    if (LDEBUG) {
      cerr << XPID << "aflowlib::PrototypeLibraries: flip_option=" << flip_option << endl;
    }
    if (LDEBUG) {
      if (uid.empty()) {
        cerr << XPID << "aflowlib::PrototypeLibraries: can't find uid" << endl;
      } else {
        cerr << XPID << "aflowlib::PrototypeLibraries: found uid=" << uid << endl;
      }
    }
    // check for ANRL
    if (!uid.empty()) {
      return anrl::PrototypeANRL_Generator_UID(uid, parameters, vatomX, vvolumeX);
    }

    vector<string> vlabel_ANRL;
    if (aurostd::string2tokens(label, vlabel_ANRL, "_") >= 4) {
      if (LDEBUG) {
        cerr << XPID << "aflowlib::PrototypeLibraries: ANRL=" << 1 << endl;
      }
      // DX20190708 START
      //  add default permutation to label if not included in input //DX20190708
      vector<string> perm_tokens;
      if (aurostd::string2tokens(label, perm_tokens, ".") == 1) {
        const string default_permutation = aurostd::RemoveNumbers(vlabel_ANRL[0]);
        label += "." + default_permutation;
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " added default permutation designation to ANRL label; label=" << label << endl;
        }
      }
      return anrl::PrototypeANRL_Generator(label, parameters, vatomX, vvolumeX); // DX20200423
    }

    message << "Could not find or create a prototype with the label " << label << endl;
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
  }
} // namespace aflowlib

// DX20210114 [MOVED NearestNeighbor() to aflow_xatom.cpp]

// ***************************************************************************
// aflowlib::PrototypesHelp
// ***************************************************************************
namespace aflowlib {
  string PrototypesHelp() {
    const anrl::ProtoData pd = anrl::ProtoData::get();
    stringstream result;
    // intro(strstream);
    result << endl;
    //  strstream << aflow::Banner("BANNER_TINY") << endl;
    result << aflow::Banner("BANNER_BIG") << endl;
    vector<string> number_of_species = pd.lookup["number_of_species"].keys();
    std::sort(number_of_species.begin(), number_of_species.end());
    for (const std::string &number : number_of_species) {
      result << "PROTOTYPES " << number << "-COMPONENTS (UID | LINK | LABEL)" << endl;
      for (const std::string &uid : static_cast<std::vector<std::string>>(pd.lookup["number_of_species"][number])) {
        if (pd.content[uid]["legacy"]) {
          result << "  " << uid << " | - | " << static_cast<string>(pd.content[uid]["label"]) << endl;
        } else {
          result << "  " << uid << " | https://aflow.org/p/" << uid << " | " << static_cast<string>(pd.content[uid]["label"]) << endl;
        }
      }
    }
    return result.str();
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::PrototypesIcsdHelp
// ***************************************************************************
namespace aflowlib {
  /// @brief write all available prototypes that correspond to an ICSD number with a certain number of elements
  /// @param options number of elements as comma seperated list, or empty for all available
  /// @return
  string PrototypesIcsdHelp(const string &options) {
    stringstream result;
    vector<string> number_of_species;
    const anrl::ProtoData pd = anrl::ProtoData::get();

    if (options.empty()) {
      number_of_species = pd.lookup["number_of_species"].keys();
    } else {
      vector<string> tokens;
      aurostd::string2tokens(options, tokens, ",");
      for (const auto &n : tokens) {
        if (pd.lookup["number_of_species"].find(n) != pd.lookup["number_of_species"].end()) {
          number_of_species.push_back(n);
        } else {
          cerr << "Ignoring " << n << " - not found in the database." << endl;
        }
      }
    }
    std::sort(number_of_species.begin(), number_of_species.end());

    for (const std::string &number : number_of_species) {
      result << number << "-COMPONENTS (ICSD | UID | LINK | LABEL)" << endl;
      for (const std::string &uid : static_cast<std::vector<std::string>>(pd.lookup["number_of_species"][number])) {
        if (pd.content[uid]["icsd"]) {
          result << "  ICSD_" << static_cast<string>(pd.content[uid]["icsd"]) << " | " << uid << " | https://aflow.org/p/" << uid << " | " << static_cast<string>(pd.content[uid]["label"]) << endl;
        }
      }
    }

    return result.str();
  }
} // namespace aflowlib

// ***************************************************************************
// aflowlib::CALCULATED
// ***************************************************************************
namespace aflowlib {
  /// @brief gather the current statistic on available AFLUX entries
  /// @return statistic summary
  /// @authors
  /// @mod{HE,20240326,changed to use AFLUX}
  string CALCULATED() {
    const std::string statistic_raw = aurostd::httpGet("http://aflow.org/API/aflowlib_stats/");
    const aurostd::JSON::object aflux_statistic = aurostd::JSON::loadString(statistic_raw);
    stringstream out;
    out << "Calculations available on AFLUX " << TODAY << endl;
    out << "Library_ICSD_CALCULATED = " << (long) aflux_statistic["Aflow_DBs"]["ICSD"]["count"] << endl;
    out << "Library_LIB0_CALCULATED = " << (long) aflux_statistic["Aflow_DBs"]["LIB0"]["count"] << endl;
    out << "Library_LIB1_CALCULATED = " << (long) aflux_statistic["Aflow_DBs"]["LIB1"]["count"] << endl;
    out << "Library_LIB2_CALCULATED = " << (long) aflux_statistic["Aflow_DBs"]["LIB2"]["count"] << endl;
    out << "Library_LIB3_CALCULATED = " << (long) aflux_statistic["Aflow_DBs"]["LIB3"]["count"] << endl;
    out << "Library_LIB4_CALCULATED = " << (long) aflux_statistic["Aflow_DBs"]["LIB4"]["count"] << endl;
    out << "Total = " << (long) aflux_statistic["Aflow_DBs"]["total"]["count"] << endl;
    return out.str();
  }
} // namespace aflowlib

// ***************************************************************************
// Aflowlib_CALCULATED_ICSD_RANDOM
// ***************************************************************************
namespace aflowlib {
  /// @brief return the AURL of a random entry from the ICSD database
  /// @return AURL
  /// @authors
  /// @mod{HE,20240326,changed to use AFLUX}
  string CALCULATED_ICSD_RANDOM() {
    const aurostd::JSON::object aflux_statistic = aurostd::JSON::loadString(aurostd::httpGet("http://aflow.org/API/aflowlib_stats/"));

    const uint rnd_choice = (uint) floor((double) aflux_statistic["Aflow_DBs"]["ICSD"]["count"] * aurostd::ran0());
    const aurostd::JSON::object selected_system = aurostd::JSON::loadString(aurostd::httpGet("http://aflow.org/API/aflux/?$catalog(ICSD),$paging(" + std::to_string(rnd_choice) + ",1)"));

    return (string) selected_system[(size_t) 0]["aurl"] + "\n";
  }
} // namespace aflowlib

// ***************************************************************************

// ***************************************************************************

#endif // _AFLOW_XPROTO_CPP
// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
