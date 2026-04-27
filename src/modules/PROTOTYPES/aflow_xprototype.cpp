// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2019
#include <cstddef>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"

#include "aflow.h"
#include "aflow_defs.h"

using std::endl;
using std::ifstream;
using std::iostream;
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
// _XPROTOTYPE_

namespace xprototype {
  // constructors
  xprototype::xprototype() {
    free();
  }
  // destructor
  xprototype::~xprototype() {
    free();
  }
  // free

  void xprototype::free() {
    verbose = false;
    // add all others
    catalog = "";
    volume = NNN;
    label = "";
    parameter_list.clear();
    parameter_values.clear();
    parameter_set_id = "";
    weblink = "";
    stoichiometry.clear();
    Pearson_symbol = "";
    space_group_number = NNN;
    space_group_symbol_H_M = "";
    space_group_symbol_Hall = "";
    space_group_symbol_Schoenflies = "";
    for (size_t i = 0; i < Wyckoff_letters.size(); i++) {
      Wyckoff_letters[i].clear();
    }
    Wyckoff_letters.clear();
    for (size_t i = 0; i < Wyckoff_site_symmetries.size(); i++) {
      Wyckoff_site_symmetries[i].clear();
    }
    Wyckoff_site_symmetries.clear();
    for (size_t i = 0; i < Wyckoff_multiplicities.size(); i++) {
      Wyckoff_multiplicities[i].clear();
    }
    Wyckoff_multiplicities.clear();
    prototype_material = "";
    common_name = "";
    mineral_name = "";
    phase = "";
    strukturbericht = "";
    similar_materials.clear();
    comments.clear();
    title = "";
  }

  const xprototype& xprototype::operator=(const xprototype& b) {      // operator=
    if (this != &b) {
      free();
      copy(b);
    }
    return *this;
  }

  void xprototype::copy(const xprototype& b) {  // copy PRIVATE
    // will populate
    verbose = b.verbose;
    catalog = b.catalog;
    volume = b.volume;
    label = b.label;
    parameter_list.clear();
    for (size_t i = 0; i < b.parameter_list.size(); i++) {
      parameter_list.push_back(b.parameter_list[i]);
    }
    parameter_values.clear();
    for (size_t i = 0; i < b.parameter_values.size(); i++) {
      parameter_values.push_back(b.parameter_values[i]);
    }
    parameter_set_id = b.parameter_set_id;
    weblink = b.weblink;
    stoichiometry.clear();
    for (size_t i = 0; i < b.stoichiometry.size(); i++) {
      stoichiometry.push_back(b.stoichiometry[i]);
    }
    Pearson_symbol = b.Pearson_symbol;
    space_group_number = b.space_group_number;
    space_group_symbol_H_M = b.space_group_symbol_H_M;
    space_group_symbol_Hall = b.space_group_symbol_Hall;
    space_group_symbol_Schoenflies = b.space_group_symbol_Schoenflies;
    for (size_t i = 0; i < Wyckoff_letters.size(); i++) {
      Wyckoff_letters[i].clear();
    }
    Wyckoff_letters.clear();
    for (size_t i = 0; i < b.Wyckoff_letters.size(); i++) {
      Wyckoff_letters.push_back(b.Wyckoff_letters[i]);
    }
    for (size_t i = 0; i < Wyckoff_site_symmetries.size(); i++) {
      Wyckoff_site_symmetries[i].clear();
    }
    Wyckoff_site_symmetries.clear();
    for (size_t i = 0; i < b.Wyckoff_site_symmetries.size(); i++) {
      Wyckoff_site_symmetries.push_back(b.Wyckoff_site_symmetries[i]);
    }
    for (size_t i = 0; i < Wyckoff_multiplicities.size(); i++) {
      Wyckoff_multiplicities[i].clear();
    }
    Wyckoff_multiplicities.clear();
    for (size_t i = 0; i < b.Wyckoff_multiplicities.size(); i++) {
      Wyckoff_multiplicities.push_back(b.Wyckoff_multiplicities[i]);
    }
    prototype_material = b.prototype_material;
    common_name = b.common_name;
    mineral_name = b.mineral_name;
    phase = b.phase;
    strukturbericht = b.strukturbericht;
    similar_materials.clear();
    for (size_t i = 0; i < b.similar_materials.size(); i++) {
      similar_materials.push_back(b.similar_materials[i]);
    }
    comments.clear();
    for (size_t i = 0; i < b.comments.size(); i++) {
      comments.push_back(b.comments[i]);
    }
    title = b.title;
  }

  void xprototype::clear() {
    const xprototype a;
    (*this) = a;
  }

  ostream& operator<<(ostream& oss, const xprototype& prototype) {
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss.precision(10);
    oss << "verbose=" << prototype.verbose << endl;
    oss << "catalog=" << prototype.catalog << endl;
    oss << "volume=" << prototype.volume << endl;
    oss << "label=" << prototype.label << endl;
    oss << "parameter_list=";
    for (size_t i = 0; i < prototype.parameter_list.size(); i++) {
      oss << prototype.parameter_list[i] << " ";
    }
    oss << endl;
    oss << "parameter_values=";
    for (size_t i = 0; i < prototype.parameter_values.size(); i++) {
      oss << prototype.parameter_values[i] << " ";
    }
    oss << endl;
    oss << "parameter_set_id=" << prototype.parameter_set_id << endl;
    oss << "weblink=" << prototype.weblink << endl;
    oss << "stoichiometry=";
    for (size_t i = 0; i < prototype.stoichiometry.size(); i++) {
      oss << prototype.stoichiometry[i] << " ";
    }
    oss << endl;
    oss << "Pearson_symbol=" << prototype.Pearson_symbol << endl;
    oss << "space_group_number=" << prototype.space_group_number << endl;
    oss << "space_group_symbol_H_M=" << prototype.space_group_symbol_H_M << endl;
    oss << "space_group_symbol_Hall=" << prototype.space_group_symbol_Hall << endl;
    oss << "space_group_symbol_Schoenflies=" << prototype.space_group_symbol_Schoenflies << endl;
    oss << "Wyckoff_letters=";
    for (size_t i = 0; i < prototype.Wyckoff_letters.size(); i++) {
      for (size_t j = 0; i < prototype.Wyckoff_letters.at(i).size(); j++) {
        oss << prototype.Wyckoff_letters.at(i).at(j) << " ";
      }
      oss << endl;
    }
    oss << "Wyckoff_site_symmetries=";
    for (size_t i = 0; i < prototype.Wyckoff_site_symmetries.size(); i++) {
      for (size_t j = 0; i < prototype.Wyckoff_site_symmetries.at(i).size(); j++) {
        oss << prototype.Wyckoff_site_symmetries.at(i).at(j) << " ";
      }
      oss << endl;
    }
    oss << "Wyckoff_multiplicities=";
    for (size_t i = 0; i < prototype.Wyckoff_multiplicities.size(); i++) {
      for (size_t j = 0; i < prototype.Wyckoff_multiplicities.at(i).size(); j++) {
        oss << prototype.Wyckoff_multiplicities.at(i).at(j) << " ";
      }
      oss << endl;
    }
    oss << "prototype_material=" << prototype.prototype_material << endl;
    oss << "common_name=" << prototype.common_name << endl;
    oss << "mineral_name=" << prototype.mineral_name << endl;
    oss << "phase=" << prototype.phase << endl;
    oss << "strukturbericht=" << prototype.strukturbericht << endl;
    oss << "similar_materials=";
    for (size_t i = 0; i < prototype.similar_materials.size(); i++) {
      oss << prototype.similar_materials[i] << " ";
    }
    oss << endl;
    oss << "comments=";
    for (size_t i = 0; i < prototype.comments.size(); i++) {
      oss << prototype.comments[i] << " ";
    }
    oss << endl;
    oss << "title=" << prototype.title << endl;
    return oss;
  }

  // ********************************************************************************************************************************************************

  void xprototype::populate(const string& prototype) {
    free();
    // DEFAULT
    verbose = false;

    // OFFSET
    if (prototype.empty()) {
      const xprototype a;
      (*this) = a;
      return;
    }
    // ********************************************************************************************************************************************************
    // A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b A2BC4D_tI16_121_d_a_i_b
    else if (prototype == "A2BC4D_tI16_121_d_a_i_b" || aurostd::toupper(prototype) == aurostd::toupper("stannite")) { // A2BC4D_tI16_121_d_a_i_b
      // stannite example
      catalog = "anrl";
      volume = 1;
      label = "A2BC4D_tI16_121_d_a_i_b";
      parameter_list.emplace_back("a");
      parameter_list.emplace_back("c/a");
      parameter_list.emplace_back("x4");
      parameter_list.emplace_back("z4");
      parameter_values.push_back(5.46);
      parameter_values.push_back(1.96428571429);
      parameter_values.push_back(0.245);
      parameter_values.push_back(0.132);
      parameter_set_id = "001";
      weblink = _AFLOW_PROTOTYPE_ENCYCLOPEDIA_ + "A2BC4D_tI16_121_d_a_i_b.html";
      stoichiometry.push_back(2);
      stoichiometry.push_back(1);
      stoichiometry.push_back(4);
      stoichiometry.push_back(1);
      // symmetry
      Pearson_symbol = "tI16";
      space_group_number = 121;
      space_group_symbol_H_M = "I-42m"; // or lookup table
      space_group_symbol_Hall = "I -4 2"; // or lookup table
      space_group_symbol_Schoenflies = "D_{2d}^{11}"; // or lookup table
      vector<string> Wyckoff_letters_tmp;
      // letters
      Wyckoff_letters_tmp.emplace_back("d");
      Wyckoff_letters.push_back(Wyckoff_letters_tmp);
      Wyckoff_letters_tmp.clear();
      Wyckoff_letters_tmp.emplace_back("a");
      Wyckoff_letters.push_back(Wyckoff_letters_tmp);
      Wyckoff_letters_tmp.clear();
      Wyckoff_letters_tmp.emplace_back("i");
      Wyckoff_letters.push_back(Wyckoff_letters_tmp);
      Wyckoff_letters_tmp.clear();
      Wyckoff_letters_tmp.emplace_back("b");
      Wyckoff_letters.push_back(Wyckoff_letters_tmp);
      Wyckoff_letters_tmp.clear();
      // site symmetries
      vector<string> Wyckoff_site_symmetries_tmp;
      Wyckoff_site_symmetries_tmp.emplace_back("-4..");
      Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);
      Wyckoff_site_symmetries_tmp.clear();
      Wyckoff_site_symmetries_tmp.emplace_back("-42m");
      Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);
      Wyckoff_site_symmetries_tmp.clear();
      Wyckoff_site_symmetries_tmp.emplace_back("..m");
      Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);
      Wyckoff_site_symmetries_tmp.clear();
      Wyckoff_site_symmetries_tmp.emplace_back("-42m");
      Wyckoff_site_symmetries.push_back(Wyckoff_site_symmetries_tmp);
      Wyckoff_site_symmetries_tmp.clear();
      vector<uint> Wyckoff_multiplicities_tmp;
      Wyckoff_multiplicities_tmp.push_back(4);
      Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);
      Wyckoff_multiplicities_tmp.clear();
      Wyckoff_multiplicities_tmp.push_back(2);
      Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);
      Wyckoff_multiplicities_tmp.clear();
      Wyckoff_multiplicities_tmp.push_back(8);
      Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);
      Wyckoff_multiplicities_tmp.clear();
      Wyckoff_multiplicities_tmp.push_back(2);
      Wyckoff_multiplicities.push_back(Wyckoff_multiplicities_tmp);
      Wyckoff_multiplicities_tmp.clear();
      // designations
      prototype_material = "Cu2FeS4Sn";
      common_name = "";
      mineral_name = "stannite";
      phase = "";
      strukturbericht = "H2_{6}";
      similar_materials.emplace_back("Cu2CdSe4Sn");
      similar_materials.emplace_back("CoCu2S4Sn");
      similar_materials.emplace_back("Cu2GeHgS4");
      similar_materials.emplace_back("Cu2HgS4Sn");
      similar_materials.emplace_back("Ag2FeS4Sn");
      comments.emplace_back(
          "If $c=2a$, $x=1/4$, and $z=3/8$, the atoms are on the sites of the diamond ($A4$) structure. If, in addition, the Cu, Fe, and Sn atoms are replaced by a single atom type, the crystal reduces to the "
          "zincblende ($B3$) structure.");
      title = "Stannite (Cu2FeS4Sn, H2_{6}) Structure";
      return;
    }
    // ********************************************************************************************************************************************************

    throw aurostd::xerror(__AFLOW_FILE__, "xprototype::xprototype():", "Prototype does not exist: " + prototype, _FILE_NOT_FOUND_);
  }
} // namespace xprototype

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2024              *
// *                                                                        *
// **************************************************************************
