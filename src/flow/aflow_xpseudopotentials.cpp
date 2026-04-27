// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2021

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_automatic_template.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xscalar.h"

#include "aflow.h"
#include "aflow_xhost.h"

using std::cerr;
using std::cout;
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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XPSEUDOPOTENTIAL_PROTOTYPES_
/*
   ./aflow --scrub=POTCAR --FILE /common/VASP/potpaw_PBE/current/Mo_pv/POTCAR
   ./aflow --scrub=OUTCAR --FILE /common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz
   ./aflow --pseudopotentials_check=/tmp/POTCAR1
   ./aflow --pseudopotentials_check=/tmp/POTCAR2
   ./aflow --pseudopotentials_check=/tmp/OUTCAR.relax2
   ./aflow --pseudopotentials_check=/common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz

   ./aflow --use_aflow.in=aflow.in --beep --force --showPID --lib2raw="/common/LIB3/LIB/TeW_pvY_sv/TFCC001.ABC"

#!/bin/sh
#echo "$1"
#echo "$2"
STR1=`cat "/common/LIB1/LIB/$1/A1/aflow.in" | grep AUID | head -1 | sed "s/\[VASP_POTCAR_AUID\]/if(AUID==\""/g`
STR2="\") {"$2"}   //   "$1
echo $STR1$STR2

*/

#define PSEUDOPOTENTIAL_GENERATOR_pad 70

namespace /* anonymous */ {

  /// @brief Returns the pspot data from the packaged in-memory json file as a vector of xPOTCAR
  /// @authors
  /// @mod{ST,20260203,created}
  std::vector<xPOTCAR> load_pspot_data() {
    const aurostd::JSON::object potcar_list = aurostd::JSON::loadString(aurostd::EmbData::get_content("aflow_xpseudopotentials_data.json", "PSPOTS"));
    std::vector<xPOTCAR> vxpseudopotential;
    vxpseudopotential.reserve(potcar_list.size());

#define JSON_SETTER_HERE(var) xpotcar.var = static_cast<decltype(xpotcar.var)>(obj[#var]);
#define MEMBERS_HERE \
  filename, AUID, vTITEL, species, species_Z, species_pp, species_pp_type, species_pp_version, species_pp_AUID, species_pp_groundstate_energy, species_pp_groundstate_structure, vEATOM, vRMAX, vLEXCH
    const auto& objs = aurostd::JSON::List(potcar_list);
    for (const auto& obj : objs) {
      xPOTCAR xpotcar;
      AST_JSON_ACCESSOR(JSON_SETTER_HERE, MEMBERS_HERE)
      vxpseudopotential.push_back(xpotcar);
    }
#undef JSON_SETTER_HERE
#undef MEMBERS_HERE

    return vxpseudopotential;
  }

  /// @brief Returns the pspot reference data from the packaged in-memory json file as a json object
  /// @authors
  /// @mod{ST,20260203,created}
  aurostd::JSON::object load_pspot_refs() {
    return aurostd::JSON::loadString(aurostd::EmbData::get_content("pseudopotential_references.json", "PSPOTS"));
  }

  /// @brief Returns the pspot data from the packaged in-memory json file as a vector of xPOTCAR using static cache
  /// @authors
  /// @mod{ST,20260203,created}
  std::vector<xPOTCAR> get_pspot_data() {
    static const std::vector<xPOTCAR> data = load_pspot_data();
    return data;
  }

  /// @brief Returns the pspot reference data from the packaged in-memory json file as a json object using static cache
  /// @authors
  /// @mod{ST,20260203,created}
  aurostd::JSON::object get_pspot_refs() {
    static const aurostd::JSON::object obj = load_pspot_refs();
    return obj;
  }

} // namespace

/// @brief Loads the pspot data from static cache. Builds from json on first call.
/// @returns pspot data
/// @authors
/// @mod{ST,20260203,created}
std::vector<xPOTCAR> get_pseudopotential_data() {
  return get_pspot_data();
}

/// @brief Loads the pspot reference data from static cache. Builds from json on first call.
/// @returns pspot reference data
/// @authors
/// @mod{ST,20260203,created}
aurostd::JSON::object get_pseudopotential_enthalpy_references() {
  return get_pspot_refs();
}

bool xPOTCAR_FixBoot(xPOTCAR& xPOT) {
  bool fix = false;
  if (xPOT.vENMAX.empty()) {
    fix = true;
    xPOT.vENMAX.push_back(NNN);
  }                                    // if not identified for BOOT STRAP
  if (xPOT.vENMIN.empty()) {
    fix = true;
    xPOT.vENMIN.push_back(NNN);
  }                                    // if not identified for BOOT STRAP
  if (xPOT.vPOMASS.empty()) {
    fix = true;
    xPOT.vPOMASS.push_back(NNN);
  }                                  // if not identified for BOOT STRAP
  if (xPOT.vZVAL.empty()) {
    fix = true;
    xPOT.vZVAL.push_back(NNN);
  }                                      // if not identified for BOOT STRAP
  if (xPOT.vEATOM.empty()) {
    fix = true;
    xPOT.vEATOM.push_back(NNN);
  }                                    // if not identified for BOOT STRAP
  if (xPOT.vRCORE.empty()) {
    fix = true;
    xPOT.vRCORE.push_back(NNN);
  }                                    // if not identified for BOOT STRAP
  if (xPOT.vRWIGS.empty()) {
    fix = true;
    xPOT.vRWIGS.push_back(NNN);
  }                                    // if not identified for BOOT STRAP
  if (xPOT.vEAUG.empty()) {
    fix = true;
    xPOT.vEAUG.push_back(NNN);
  }                                      // if not identified for BOOT STRAP
  if (xPOT.vRAUG.empty()) {
    fix = true;
    xPOT.vRAUG.push_back(NNN);
  }                                      // if not identified for BOOT STRAP
  if (xPOT.vRMAX.empty()) {
    fix = true;
    xPOT.vRMAX.push_back(NNN);
  }                                      // if not identified for BOOT STRAP
  if (xPOT.vTITEL.empty()) {
    fix = true;
    xPOT.vTITEL.emplace_back("N/A");
  }                                  // if not identified for BOOT STRAP
  if (xPOT.vLEXCH.empty()) {
    fix = true;
    xPOT.vLEXCH.emplace_back("N/A");
  }                                  // if not identified for BOOT STRAP
  if (xPOT.species.empty()) {
    fix = true;
    xPOT.species.emplace_back("N/A");
  }                                // if not identified for BOOT STRAP
  if (xPOT.species_Z.empty()) {
    fix = true;
    xPOT.species_Z.push_back(0);
  }                                // if not identified for BOOT STRAP
  if (xPOT.species_pp.empty()) {
    fix = true;
    xPOT.species_pp.emplace_back("N/A");
  }                          // if not identified for BOOT STRAP
  if (xPOT.species_pp_type.empty()) {
    fix = true;
    xPOT.species_pp_type.emplace_back("N/A");
  }                // if not identified for BOOT STRAP
  if (xPOT.species_pp_version.empty()) {
    fix = true;
    xPOT.species_pp_version.emplace_back("N/A");
  }          // if not identified for BOOT STRAP
  if (xPOT.species_pp_AUID.empty()) {
    fix = true;
    xPOT.species_pp_AUID.emplace_back("N/A");
  }                // if not identified for BOOT STRAP
  if (xPOT.species_pp_groundstate_energy.empty()) {
    fix = true;
    xPOT.species_pp_groundstate_energy.push_back(NNN);
  }          // if not identified for BOOT STRAP
  if (xPOT.species_pp_groundstate_structure.empty()) {
    fix = true;
    xPOT.species_pp_groundstate_structure.emplace_back("N/A");
  }  // if not identified for BOOT STRAP
  return fix;
}

/// @brief Locates the first potcar data corresponding to given values
/// @param species_pp_AUID[out] AUID of matching species
/// @param species_pp_AUID_collisions[out] AUIDs of colliding matching species - not used - no collisions
/// @param TITEL[in] corresponding potcar value
/// @param LEXCH[in] corresponding potcar value
/// @param EATOM[in] corresponding potcar value
/// @param RMAX[in] corresponding potcar value
/// @param LVERBOSE[in] true to log verbosely
/// @return Matching xPOTCAR
xPOTCAR xPOTCAR_Finder(vector<string>& species_pp_AUID, vector<string>& species_pp_AUID_collisions, const string& TITEL, const string& LEXCH, const double& EATOM, const double& RMAX, bool LVERBOSE) {
  xPOTCAR xPOT;
  bool found = false;
  const auto vxpseudopotential = get_pseudopotential_data();
  for (size_t ipp = 0; ipp < vxpseudopotential.size(); ipp++) {
    const bool test = (TITEL == vxpseudopotential[ipp].vTITEL.at(0))
                   && (LEXCH == vxpseudopotential[ipp].vLEXCH.at(0))
                   && (aurostd::abs(EATOM - vxpseudopotential[ipp].vEATOM.at(0)) < 0.00001)
                   && (aurostd::abs(RMAX - vxpseudopotential[ipp].vRMAX.at(0)) < 0.0001);
    if (test) {
      found = true;
      if (LVERBOSE) {
        cerr << XPID << "xPOTCAR::xPOTCAR_Finder: FOUND: POTCAR=" << vxpseudopotential[ipp].filename << endl;
      }
      species_pp_AUID.push_back(vxpseudopotential[ipp].AUID);
      xPOT = vxpseudopotential[ipp];
      break;
    }
  }
  if (!found) {
    if (!vxpseudopotential.empty()) {
      cerr << XPID << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: TITEL=" << TITEL << " LEXCH=" << LEXCH << " EATOM=" << EATOM << " RMAX =" << RMAX << endl;
    }
    species_pp_AUID.emplace_back("N/A");
  }

  xPOTCAR_FixBoot(xPOT);
  return xPOT;
}

/// @brief Locates the first potcar data with the given AUID
/// @param AUID AUID to find
/// @param LVERBOSE true to log verbosely
/// @return Matching xPOTCAR
xPOTCAR xPOTCAR_Finder(const string& AUID, bool LVERBOSE) {
  xPOTCAR xPOT;
  bool found = false;
  const auto vxpseudopotential = get_pseudopotential_data();
  for (size_t ipp = 0; ipp < vxpseudopotential.size(); ipp++) {
    if (AUID == vxpseudopotential[ipp].AUID) {
      found = true;
      if (LVERBOSE) {
        cerr << XPID << "xPOTCAR::xPOTCAR_Finder: FOUND: POTCAR=" << vxpseudopotential[ipp].filename << endl;
      }
      xPOT = vxpseudopotential[ipp];
      break;
    }
  }
  if (!found) {
    if (!vxpseudopotential.empty()) {
      cerr << XPID << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: AUID=" << AUID << endl;
    }
  }

  xPOTCAR_FixBoot(xPOT);
  return xPOT;
}

bool xPOTCAR_PURE_Printer(xPOTCAR& xPOT, ostream& oss, bool LVERBOSE) {
  if (XHOST.PSEUDOPOTENTIAL_GENERATOR && xPOT.species.size() == 1) {  // SC20200326
    string comment;

    double groundstate_energy = NNN;
    string groundstate_structure = "N/A_" + xPOT.AUID;

    double volume_atom;
    double spin_atom;

    const bool found = xPOTCAR_EnthalpyReference_AUID(xPOT.AUID, "", groundstate_structure, groundstate_energy, volume_atom, spin_atom);

    if (found) {
      cerr << XPID << "xPOTCAR_PURE_Printer:     FOUND AUID=" << xPOT.AUID << endl;
    }
    if (!found) {
      cerr << XPID << "xPOTCAR_PURE_Printer: NOT_FOUND AUID=" << xPOT.AUID << endl;
    }

    comment = xPOT.species_pp_version.at(0);
    xPOT.vTITEL.at(0) = aurostd::RemoveWhiteSpaces(xPOT.vTITEL.at(0));
    xPOT.vLEXCH.at(0) = aurostd::RemoveWhiteSpaces(xPOT.vLEXCH.at(0));
    if (LVERBOSE) {
      cerr << XPID << "xPOTCAR_PURE_Printer: vTITEL.at(0)=" << xPOT.vTITEL.at(0) << endl;
    }   // SC20200326
    if (LVERBOSE) {
      cerr << XPID << "xPOTCAR_PURE_Printer: species.at(0)=" << xPOT.species.at(0) << endl;
    }   // SC20200326
    if (LVERBOSE) {
      cerr << XPID << "xPOTCAR_PURE_Printer: species_Z.at(0)=" << xPOT.species_Z.at(0) << endl;
    }    // SC20200326
    oss << "  " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl;    // SC20200326
    oss << "  " << aurostd::PaddedPOST("{", PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("xPOTCAR x;", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.filename=\"" + xPOT.filename + "\";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.AUID=\"" + xPOT.AUID + "\";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vTITEL.push_back(\"" + xPOT.vTITEL.at(0) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species.push_back(\"" + xPOT.species.at(0) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_Z.push_back(" + aurostd::utype2string<int>(xPOT.species_Z.at(0)) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp.push_back(\"" + xPOT.species_pp.at(0) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_type.push_back(\"" + xPOT.species_pp_type.at(0) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_version.push_back(\"" + xPOT.species_pp_version.at(0) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_AUID.push_back(\"" + xPOT.AUID + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    if (!found) {
      oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back(" + aurostd::utype2string<double>(groundstate_energy, 10) + ");//" + xPOT.AUID, PSEUDOPOTENTIAL_GENERATOR_pad) << "    // "
          << comment << endl;    // SC20200326
    }
    if (found) {
      oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back(" + aurostd::utype2string<double>(groundstate_energy, 10) + ");//", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    }
    oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_structure.push_back(\"" + groundstate_structure + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vEATOM.push_back(" + aurostd::utype2string<double>(xPOT.vEATOM.at(0), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vRMAX.push_back(" + aurostd::utype2string<double>(xPOT.vRMAX.at(0), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vLEXCH.push_back(\"" + xPOT.vLEXCH.at(0) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vxpseudopotential.push_back(x);", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "  " << aurostd::PaddedPOST("}", PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl; // SC20200326
    oss << "  // " << xPOT.filename << endl; // SC20200326
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << endl; // SC20200326
    return true;
  }
  return false;
}

ostream& operator<<(ostream& oss, const xPOTCAR& xPOT) {
  oss.setf(std::ios::fixed, std::ios::floatfield);
  oss.precision(10);
  for (size_t i = 0; i < xPOT.species.size(); i++) {
    const string comment = xPOT.species_pp_version.at(i);
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << "  // [AFLOW]START=" << comment << " " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("filename=\"" + xPOT.filename + "\";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("AUID=\"" + xPOT.AUID + "\";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vTITEL.push_back(\"" + xPOT.vTITEL.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("pp_type=\"" + xPOT.pp_type + "\";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species.push_back(\"" + xPOT.species[i] + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_Z.push_back(" + aurostd::utype2string<int>(xPOT.species_Z.at(i)) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp.push_back(\"" + xPOT.species_pp.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_type.push_back(\"" + xPOT.species_pp_type.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_version.push_back(\"" + xPOT.species_pp_version.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_AUID.push_back(\"" + xPOT.species_pp_AUID.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_energy.push_back(" + aurostd::utype2string<double>(xPOT.species_pp_groundstate_energy.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad)
        << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_structure.push_back(\"" + xPOT.species_pp_groundstate_structure.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_PAW=" + aurostd::bool2string(xPOT.POTCAR_PAW) + ";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_TYPE=\"" + xPOT.POTCAR_TYPE + "\";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_KINETIC=" + aurostd::bool2string(xPOT.POTCAR_KINETIC) + ";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_GW=" + aurostd::bool2string(xPOT.POTCAR_GW) + ";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_AE=" + aurostd::bool2string(xPOT.POTCAR_AE) + ";", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMAX.push_back(" + aurostd::utype2string<double>(xPOT.vENMAX.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMIN.push_back(" + aurostd::utype2string<double>(xPOT.vENMIN.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vPOMASS.push_back(" + aurostd::utype2string<double>(xPOT.vPOMASS.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vZVAL.push_back(" + aurostd::utype2string<double>(xPOT.vZVAL.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vEATOM.push_back(" + aurostd::utype2string<double>(xPOT.vEATOM.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRCORE.push_back(" + aurostd::utype2string<double>(xPOT.vRCORE.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRWIGS.push_back(" + aurostd::utype2string<double>(xPOT.vRWIGS.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vEAUG.push_back(" + aurostd::utype2string<double>(xPOT.vEAUG.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRAUG.push_back(" + aurostd::utype2string<double>(xPOT.vRAUG.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRMAX.push_back(" + aurostd::utype2string<double>(xPOT.vRMAX.at(i), 10) + ");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "    " << aurostd::PaddedPOST("vLEXCH.push_back(\"" + xPOT.vLEXCH.at(i) + "\");", PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl; // SC20200326
    oss << "  // " << xPOT.filename << endl; // SC20200326
    oss << "  // [AFLOW]STOP=" << comment << " " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
  }
  return oss;
}

bool xPOTCAR_EnthalpyReference_AUID(string AUID, string METAGGA) {
  string groundstate_structure;
  double groundstate_energy = 0.0;
  double volume_atom = 0.0;
  double spin_atom = 0.0;
  return xPOTCAR_EnthalpyReference_AUID(AUID, METAGGA, groundstate_structure, groundstate_energy, volume_atom, spin_atom);
}

bool xPOTCAR_EnthalpyReference_AUID(string AUID, string METAGGA, string& groundstate_structure, double& groundstate_energy, double& volume_atom, double& spin_atom) {
  const bool LDEBUG = (false || XHOST.DEBUG);
  const bool VERBOSE = false; // true;
  if (LDEBUG) {
    cerr << XPID << "xPOTCAR_EnthalpyReference_AUID: [BEGIN]" << endl;
  }
  bool found = false;
  if (LDEBUG) {
    cout << "ERROR (xPOTCAR_EnthalpyReference_AUID): AUID=[" << AUID << "]" << endl; // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  }
  if (LDEBUG) {
    cout << "ERROR (xPOTCAR_EnthalpyReference_AUID): METAGGA=[" << METAGGA << "]" << endl; // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  }

  bool nKIN = false;
  bool SCAN = false;
  if (METAGGA.empty() || METAGGA == "none" || METAGGA == "NONE") {
    nKIN = true;
    SCAN = false;
  }
  if (METAGGA == "SCAN" || METAGGA == "scan") {
    nKIN = false;
    SCAN = true;
  }

  if (!XHOST.PSEUDOPOTENTIAL_GENERATOR) {
    if (VERBOSE) {
      cout << "xPOTCAR_EnthalpyReference_AUID: AUID=[" << AUID << "]  METAGGA=[" << METAGGA << "]" << endl;
    }
  }

  const aurostd::JSON::object pseudopotential_refs = get_pseudopotential_enthalpy_references();

  const auto& objs = aurostd::JSON::List(pseudopotential_refs);
  for (const auto& obj : objs) {
    if (static_cast<string>(obj["AUID"]) == AUID && static_cast<string>(obj["METAGGA"]) == METAGGA) {
      found = true;
      groundstate_structure = static_cast<string>(obj["groundstate_structure"]);
      groundstate_energy = static_cast<double>(obj["groundstate_energy"]);
      volume_atom = static_cast<double>(obj["volume_atom"]);
      spin_atom = static_cast<double>(obj["spin_atom"]);
    }
  }

  // ./xgo Sm_3:PAW_GGA:11May2000 "found=true;groundstate_structure=\"ICSD_246657\";groundstate_energy=-4.621400;volume_atom=33.447633;spin_atom=0.0;"// FIX
  // ./xgo Sm_3:PAW_GGA:11May2000 && 0 "found=true;groundstate_structure=\"ICSD_652637\";groundstate_energy=-4.64136;volume_atom=33.5075;spin_atom=0.0;"// IT HAS LDAU
  // ./xgo Sm_3:PAW_PBE:07Sep2000 && 0 "found=true;groundstate_structure=\"A1\";groundstate_energy=-4.7062;volume_atom=33.8339;spin_atom=0.0;"// IT HAS LDAU
  // ./xgo Ce "found=true;groundstate_structure=\"A1\";groundstate_energy=-5.92998;volume_atom=26.0579;spin_atom=0.0;"
  // ./xgo Ce "found=true;groundstate_structure=\"ICSD_2284-mS4\";groundstate_energy=-5.93013;volume_atom=26.0697;spin_atom=0.0;"
  // ./xgo Cl_h:PAW_PBE:08Apr2002 "found=true;groundstate_structure=\"A11\";groundstate_energy=-1.8156;volume_atom=37.3299;spin_atom=0.0;"WAITING

  if (!found) {
    volume_atom = 999999, spin_atom = 999999;
  } // some defaults
  //  if(!found) cerr <<"ERROR (xPOTCAR_EnthalpyReference_AUID): NOT FOUND: AUID=" << AUID << endl;// throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  if (LDEBUG && !found) {
    cout << "ERROR (xPOTCAR_EnthalpyReference_AUID): NOT FOUND: AUID=" << AUID << endl; // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  }
  if (LDEBUG && found) {
    cout << "ERROR (xPOTCAR_EnthalpyReference_AUID): FOUND: AUID=" << AUID << endl; // throw aurostd::xerror(__AFLOW_FILE__,XPID+"xPOTCAR_EnthalpyReference_AUID():","Throw for debugging purposes.",_GENERIC_ERROR_);
  }
  if (LDEBUG) {
    cerr << XPID << "xPOTCAR_EnthalpyReference_AUID: [END]" << endl;
  }
  return found;
}

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
