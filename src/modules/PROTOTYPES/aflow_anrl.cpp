// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2016
// Written by David Hicks (DX) - 2016-2021 (generic prototype generator)

#ifndef _AFLOW_ANRL_CPP
#define _AFLOW_ANRL_CPP

#include "aflow_anrl.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <deque>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser.h"
#include "AUROSTD/aurostd_xscalar.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "aflow_aflowrc.h"
#include "aflow_defs.h"
#include "aflow_xhost.h"
#include "extern/SYMBOLICCPLUSPLUS/symbolic.h"
#include "flow/aflow_ivasp.h"
#include "flow/aflow_pflow.h"
#include "flow/aflow_support_types.h"
#include "interfaces/aflow_symbolic.h"
#include "modules/COMPARE/aflow_compare_structure.h"
#include "modules/SYM/aflow_symmetry.h"
#include "modules/SYM/aflow_symmetry_spacegroup.h"
#include "modules/SYM/aflow_wyckoff.h"
#include "structure/aflow_lattice.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

#define _DEBUG_ANRL_ false // DX20200625

// ***************************************************************************
// anrl::PrototypeANRL_Consistency()
// ***************************************************************************
namespace anrl {

  bool PrototypeANRL_Consistency(uint vparameters_size, uint nparameters, string prototype, string label, string Strukturbericht, string Pearson_symbol, uint spacegroup, string params, uint print_mode) { // DX20180710 - added print_mode info //DX20200207 - oss no longer needed

    if (vparameters_size != nparameters && print_mode != 1) { // DX20180710 - if equations only (print_mode==1), we do not need the parameters
      stringstream message;
      message << "anrl::PrototypeANRL" << endl;
      message << " Prototype                   : " << prototype << endl;
      message << " AFLOW prototype label       : " << label << endl;
      message << " Strukturbericht Designation : " << Strukturbericht << endl;
      message << " Pearson Symbol              : " << Pearson_symbol << endl;
      message << " Space group number          : " << GetSpaceGroupName(spacegroup) << endl;
      message << " Space group symbol          : " << spacegroup << endl;
      message << " AFLOW prototype command     : aflow --proto=" << label << endl;
      message << "                                     --params=" << params << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_); // DX20200207 - return false -> throw (throw here instead of in another function)
    }
    return true;
  }
} // namespace anrl

// ***************************************************************************
// anrl::vproto2tokens()
// ***************************************************************************
namespace anrl {
  bool vproto2tokens(string proto_line, string& label, uint& nspecies, uint& natoms, uint& spacegroup, uint& nunderscores, uint& nparameters, string& Pearson_symbol, string& params, string& Strukturbericht, string& prototype, string& dialect) {
    vector<string> tokens;
    uint j = 0;

    aurostd::string2tokens(proto_line, tokens, ";");
    label = tokens.at(j++);
    nspecies = aurostd::string2utype<uint>(tokens.at(j++));
    natoms = aurostd::string2utype<uint>(tokens.at(j++));
    spacegroup = aurostd::string2utype<uint>(tokens.at(j++));
    nunderscores = aurostd::string2utype<uint>(tokens.at(j++));
    nparameters = aurostd::string2utype<uint>(tokens.at(j++));
    Pearson_symbol = tokens.at(j++);
    params = tokens.at(j++);
    Strukturbericht = tokens.at(j++);
    prototype = tokens.at(j++);
    dialect = tokens.at(j++);

    return true;
  }
} // namespace anrl

// DX20190208 - for making ANRL label -START
//  ***************************************************************************
//  anrl::extractStoichiometry()
//  ***************************************************************************
namespace anrl {
  vector<uint> extractStoichiometry(string& anrl_label) {
    vector<uint> stoichiometry;

    vector<string> tokens;
    aurostd::string2tokens(anrl_label, tokens, "_");
    string stoichiometry_string = tokens[0];

    bool is_previous_alpha = false;
    for (size_t i = 0; i < stoichiometry_string.size(); i++) {
      if (is_previous_alpha && isalpha(stoichiometry_string[i])) {
        stoichiometry.push_back(1);
      } else if (isdigit(stoichiometry_string[i])) {
        stringstream tmp;
        tmp << stoichiometry_string[i];
        stoichiometry.push_back(aurostd::string2utype<uint>(tmp.str()));
      }
      is_previous_alpha = isalpha(stoichiometry_string[i]);
    }
    if (is_previous_alpha) {
      stoichiometry.push_back(1);
    }
    return stoichiometry;
  }
} // namespace anrl
// DX20190208 - for making ANRL label - END

// ***************************************************************************
// anrl::rhl2hex()
// ***************************************************************************
namespace anrl {
  xstructure rhl2hex(const xstructure& str, double& a, double& c) {
    // RHL to HEX transformation

    xstructure hex_str; // make new xstructure object
    hex_str = str;

    hex_str.atoms.clear();

    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    const xmatrix<double> rhl_lattice;
    const xmatrix<double> hex_lattice;
    const xmatrix<double> rtransf;
    const xmatrix<double> htransf;

    // HEX lattice
    a1 = (1.0 / 2.0) * a * xn - (sqrt(3.0) / 2.0) * a * yn;
    a2 = (1.0 / 2.0) * a * xn + (sqrt(3.0) / 2.0) * a * yn;
    a3 = c * zn;
    hex_str.lattice(1, 1) = a1(1);
    hex_str.lattice(1, 2) = a1(2);
    hex_str.lattice(1, 3) = a1(3);
    hex_str.lattice(2, 1) = a2(1);
    hex_str.lattice(2, 2) = a2(2);
    hex_str.lattice(2, 3) = a2(3);
    hex_str.lattice(3, 1) = a3(1);
    hex_str.lattice(3, 2) = a3(2);
    hex_str.lattice(3, 3) = a3(3);

    hex_str.FixLattices(); // Reciprocal/f2c/c2f

    // RHL Transformation matrix
    rtransf(1, 1) = (1.0 / 2.0);
    rtransf(1, 2) = -(1.0 / (2.0 * sqrt(3.0)));
    rtransf(1, 3) = (1.0 / 3.0);
    rtransf(2, 1) = 0.0;
    rtransf(2, 2) = (1.0 / sqrt(3.0));
    rtransf(2, 3) = (1.0 / 3.0);
    rtransf(3, 1) = -(1.0 / 2.0);
    rtransf(3, 2) = -(1.0 / (2.0 * sqrt(3.0)));
    rtransf(3, 3) = (1.0 / 3.0);

    // HEX Transformtion matrix
    htransf(1, 1) = (1.0 / 2.0);
    htransf(1, 2) = -(sqrt(3.0) / 2.0);
    htransf(1, 3) = 0.0;
    htransf(2, 1) = (1.0 / 2.0);
    htransf(2, 2) = (sqrt(3.0) / 2.0);
    htransf(2, 3) = 0.0;
    htransf(3, 1) = 0.0;
    htransf(3, 2) = 0.0;
    htransf(3, 3) = 1.0;

    const xvector<double> c1(3);
    c1(1) = (2.0 / 3.0);
    c1(2) = (1.0 / 3.0);
    c1(3) = (1.0 / 3.0); // centering translation
    const xvector<double> c2(3);
    c2(1) = (1.0 / 3.0);
    c2(2) = (2.0 / 3.0);
    c2(3) = (2.0 / 3.0); // centering translation

    _atom atom_tmp; // DX20200907
    for (size_t a = 0; a < str.atoms.size(); a++) {
      atom_tmp.clear(); // DX20200907
      atom_tmp.name = str.atoms[a].name;
      atom_tmp.type = str.atoms[a].type;
      atom_tmp.basis = str.atoms[a].basis;
      xvector<double> center_pos;
      center_pos = trasp(inverse(htransf)) * (trasp(rtransf) * str.atoms[a].fpos); // Method for transforming RHL to HEX
      atom_tmp.fpos = center_pos;
      hex_str.comp_each_type.at(atom_tmp.type) += 1.0;
      hex_str.atoms.push_back(atom_tmp);
      // add centering c1
      atom_tmp.fpos = center_pos + c1;
      hex_str.comp_each_type.at(atom_tmp.type) += 1.0;
      hex_str.atoms.push_back(atom_tmp);
      // add centering c2
      atom_tmp.fpos = center_pos + c2;
      hex_str.comp_each_type.at(atom_tmp.type) += 1.0;
      hex_str.atoms.push_back(atom_tmp);
    }
    return hex_str;
  }
} // namespace anrl

// DX20190208 - for making ANRL label -START
//  ***************************************************************************
//  anrl::groupedWyckoffPosition2ANRLString()
//  ***************************************************************************
namespace anrl {
  string groupedWyckoffPosition2ANRLString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize) {
    vector<string> all_Wyckoff_sets;
    for (size_t i = 0; i < grouped_positions.size(); i++) {
      vector<string> Wyckoff_letters = grouped_positions[i].letters;
      if (alphabetize) {
        std::sort(Wyckoff_letters.begin(), Wyckoff_letters.end());
      }
      vector<string> Wyckoff_set;
      vector<uint> Wyckoff_set_count;
      for (size_t j = 0; j < Wyckoff_letters.size(); j++) {
        bool letter_stored = false;
        for (size_t k = 0; k < Wyckoff_set.size(); k++) {
          if (Wyckoff_letters[j] == Wyckoff_set[k]) {
            Wyckoff_set_count[k] = Wyckoff_set_count[k] + 1;
            letter_stored = true;
            break;
          }
        }
        if (!letter_stored) {
          Wyckoff_set.push_back(Wyckoff_letters[j]);
          Wyckoff_set_count.push_back(1);
        }
      }

      string tmp;
      for (size_t j = 0; j < Wyckoff_set.size(); j++) {
        if (Wyckoff_set_count[j] == 1) {
          tmp += Wyckoff_set[j];
        } else {
          tmp += aurostd::utype2string<uint>(Wyckoff_set_count[j]) + Wyckoff_set[j];
        }
      }
      all_Wyckoff_sets.push_back(tmp);
    }
    return aurostd::joinWDelimiter(all_Wyckoff_sets, "_");
  }
} // namespace anrl

// ***************************************************************************
// anrl::getANRLLatticeParameterString()
// ***************************************************************************
namespace anrl {
  vector<string> getANRLLatticeParameterString(char& lattice_type) {
    // get lattice parameter degrees of freedom based on lattice type

    vector<string> lattice_parameter_list;
    // triclinic
    if (lattice_type == 'a') {
      lattice_parameter_list.emplace_back("a");
      lattice_parameter_list.emplace_back("b/a");
      lattice_parameter_list.emplace_back("c/a");
      lattice_parameter_list.emplace_back("alpha");
      lattice_parameter_list.emplace_back("beta");
      lattice_parameter_list.emplace_back("gamma");
    }
    // monoclinic
    else if (lattice_type == 'm') {
      lattice_parameter_list.emplace_back("a");
      lattice_parameter_list.emplace_back("b/a");
      lattice_parameter_list.emplace_back("c/a");
      lattice_parameter_list.emplace_back("beta");
    }
    // orthorhombic
    else if (lattice_type == 'o') {
      lattice_parameter_list.emplace_back("a");
      lattice_parameter_list.emplace_back("b/a");
      lattice_parameter_list.emplace_back("c/a");
    }
    // tetragonal/trigonal/hexagonal
    else if (lattice_type == 't' || lattice_type == 'h') {
      lattice_parameter_list.emplace_back("a");
      lattice_parameter_list.emplace_back("c/a");
    }
    // cubic
    else if (lattice_type == 'c') {
      lattice_parameter_list.emplace_back("a");
    }

    return lattice_parameter_list;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getANRLLatticeParameterValuesFromWyccar()
// ***************************************************************************
namespace anrl {
  vector<double> getANRLLatticeParameterValuesFromWyccar(const vector<string>& wyccar_ITC, char lattice_type, char lattice_centering, uint setting) {
    // get lattice parameter values from the WYCCAR based on the degrees of freedom given the lattice type

    // extract lattice parameters from wyccar
    const vector<double> all_lattice_parameters = SYM::ExtractLatticeParametersFromWyccar(wyccar_ITC);

    return getANRLLatticeParameterValues(all_lattice_parameters, lattice_type, lattice_centering, setting);
  }

  vector<double> getANRLLatticeParameterValuesFromABCAngles(const xstructure& xstr, char lattice_type, char lattice_centering, uint setting) {
    // get lattice parameter values from the XSTRUCTURE

    // extract lattice parameters from xstructure
    vector<double> all_lattice_parameters;
    all_lattice_parameters.push_back(xstr.a);
    all_lattice_parameters.push_back(xstr.b);
    all_lattice_parameters.push_back(xstr.c);
    all_lattice_parameters.push_back(xstr.alpha);
    all_lattice_parameters.push_back(xstr.beta);
    all_lattice_parameters.push_back(xstr.gamma);

    // ensure all lattice parameters have been set (i.e., not zero or negative)
    for (size_t i = 0; i < all_lattice_parameters.size(); i++) {
      if (all_lattice_parameters[i] <= _ZERO_TOL_) {
        stringstream message;
        message << "The " << i << "th lattice parameter is ill-defined. Please check input.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
    }

    return getANRLLatticeParameterValues(all_lattice_parameters, lattice_type, lattice_centering, setting);
  }

  vector<double> getANRLLatticeParameterValues(const vector<double>& all_lattice_parameters, char lattice_type, char lattice_centering, uint setting) {
    // store only relevant parameter values
    vector<double> lattice_parameter_values;
    // triclinic
    if (lattice_type == 'a') {
      lattice_parameter_values.push_back(all_lattice_parameters[0]); // a
      lattice_parameter_values.push_back(all_lattice_parameters[1] / all_lattice_parameters[0]); // b/a
      lattice_parameter_values.push_back(all_lattice_parameters[2] / all_lattice_parameters[0]); // c/a
      lattice_parameter_values.push_back(all_lattice_parameters[3]); // alpha
      lattice_parameter_values.push_back(all_lattice_parameters[4]); // beta
      lattice_parameter_values.push_back(all_lattice_parameters[5]); // gamma
    }
    // monoclinic
    else if (lattice_type == 'm') {
      lattice_parameter_values.push_back(all_lattice_parameters[0]); // a
      lattice_parameter_values.push_back(all_lattice_parameters[1] / all_lattice_parameters[0]); // b/a
      lattice_parameter_values.push_back(all_lattice_parameters[2] / all_lattice_parameters[0]); // c/a
      lattice_parameter_values.push_back(all_lattice_parameters[4]); // beta
    }
    // orthorhombic
    else if (lattice_type == 'o') {
      lattice_parameter_values.push_back(all_lattice_parameters[0]); // a
      lattice_parameter_values.push_back(all_lattice_parameters[1] / all_lattice_parameters[0]); // b/a
      lattice_parameter_values.push_back(all_lattice_parameters[2] / all_lattice_parameters[0]); // c/a
    }
    // tetragonal/trigonal/hexagonal
    else if (lattice_type == 't' || lattice_type == 'h') {
      // if rhl setting, need to get correct a and c from a' and alpha'
      if (lattice_type == 'h' && lattice_centering == 'R' && (setting == SG_SETTING_1 || setting == SG_SETTING_ANRL)) {
        const double a_prime = all_lattice_parameters[0]; // a'
        const double alpha_prime = all_lattice_parameters[3] * deg2rad; // alpha'
        // see ITC (5th edition) pg. 16 for conversion
        const double c = a_prime * aurostd::sqrt(3.0) * aurostd::sqrt(1.0 + 2.0 * cos(alpha_prime));
        const double a = a_prime * aurostd::sqrt(2.0) * aurostd::sqrt(1.0 - cos(alpha_prime));
        lattice_parameter_values.push_back(a); // a
        lattice_parameter_values.push_back(c / a); // c/a
      } else {
        lattice_parameter_values.push_back(all_lattice_parameters[0]); // a
        lattice_parameter_values.push_back(all_lattice_parameters[2] / all_lattice_parameters[0]); // c/a
      }
    }
    // cubic
    else if (lattice_type == 'c') {
      lattice_parameter_values.push_back(all_lattice_parameters[0]); // a
    }

    return lattice_parameter_values;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getANRLSettingChoice()
// ***************************************************************************
namespace anrl {
  uint getANRLSettingChoice(int spacegroup) { // DX20191031 - remove reference

    // ANRL setting choice
    // rhl: rhombohedral setting: setting=1
    // monoclinic: unique axis-b: setting=1
    // centrosymmetric: origin on inversion site: setting=2

    uint anrl_setting = 1;

    // check for centrosymmetric cases
    if (spacegroup == 48 || spacegroup == 50 || spacegroup == 59 || spacegroup == 68 || spacegroup == 70 || spacegroup == 85 || spacegroup == 86 || spacegroup == 88 || spacegroup == 125 || spacegroup == 126 ||
        spacegroup == 129 || spacegroup == 130 || spacegroup == 133 || spacegroup == 134 || spacegroup == 137 || spacegroup == 138 || spacegroup == 141 || spacegroup == 142 || spacegroup == 201 ||
        spacegroup == 203 || spacegroup == 222 || spacegroup == 224 || spacegroup == 227 || spacegroup == 228) {
      anrl_setting = 2;
    }
    return anrl_setting;
  }
} // namespace anrl

// ***************************************************************************
// anrl::structure2anrl() [FROM COMMAND-LINE]
// ***************************************************************************
namespace anrl {
  string structure2anrl(istream& input, aurostd::xoption& vpflow) {
    // determine anrl label, parameters, and parameter values of the input structure

    bool recalculate_symmetry = true; // DX20191030
    uint setting = SG_SETTING_ANRL; // anrl setting choice is default

    const string usage = "aflow --prototype < POSCAR"; // DX20200721

    // load structure
    xstructure xstr(input, IOAFLOW_AUTO);

    // ---------------------------------------------------------------------------
    // ensure structure is alphabetic, otherwise the prototype convention breaks //DX20210706
    xstr.SpeciesPutAlphabetic();
    std::stable_sort(xstr.atoms.begin(), xstr.atoms.end(), sortAtomsNames);
    xstr.MakeBasis();
    xstr.MakeTypes();

    // DX20191217 START
    //  ---------------------------------------------------------------------------
    //  print format
    string format = "text";
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      format = "text";
    }
    if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      format = "json";
    }
    // DX20191217 END

    // DX20191030 - add force Wyckoff choice option - START
    //  check if forcing certain Wyckoff convention
    //  Wyckoff positions must be provided (either in CIF, Wyccar, or in Wyckoff object)
    if (vpflow.flag("STRUCTURE2ANRL::FORCE_WYCKOFF")) {
      // check if Wyckoff information is available
      if (xstr.wyckoff_sites_ITC.empty()) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot use --force_Wyckoff option, Wyckoff positions must be given.", _INPUT_ILLEGAL_);
      }
      if (vpflow.flag("STRUCTURE2ANRL::SETTING")) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot use options --setting and --force_Wyckoff together.", _INPUT_AMBIGUOUS_);
      }
      setting = xstr.spacegroupnumberoption;
      recalculate_symmetry = false;
    }
    // DX20191030 - add force Wyckoff choice option - END

    // symmetry tolerance
    const double tolerance = pflow::getSymmetryTolerance(xstr, vpflow.getattachedscheme("STRUCTURE2ANRL::TOLERANCE")); // DX20200820 - consolidated setting tolerance into a function

    // space group setting
    setting = pflow::getSpaceGroupSetting(vpflow.getattachedscheme("STRUCTURE2ANRL::SETTING"), SG_SETTING_ANRL); // DX20210421 - consolidated space group setting into function

    // print element names //DX20210622
    bool print_element_names = false;
    if (vpflow.flag("STRUCTURE2ANRL::PRINT_ELEMENT_NAMES")) {
      print_element_names = true;
    }

    // print atomic number //DX20210622
    bool print_atomic_numbers = false;
    if (vpflow.flag("STRUCTURE2ANRL::PRINT_ATOMIC_NUMBERS")) {
      print_atomic_numbers = true;
    }

    return structure2anrl(xstr, tolerance, setting, recalculate_symmetry, print_element_names, print_atomic_numbers);
  }
} // namespace anrl

// ***************************************************************************
// anrl::structure2anrl() OVERLOADS
// ***************************************************************************
namespace anrl {
  string structure2anrl(xstructure& xstr, bool recalculate_symmetry) { // DX20190829 - added recalculate_symmetry
    // determine anrl label, parameters, and parameter values of the input structure
    const double default_tolerance = SYM::defaultTolerance(xstr);
    const uint setting = SG_SETTING_ANRL; // anrl setting choice is default
    return structure2anrl(xstr, default_tolerance, setting, recalculate_symmetry); // DX20190829 - added recalculate_symmetry
  }
} // namespace anrl

// ***************************************************************************
namespace anrl {
  string structure2anrl(xstructure& xstr, double tolerance) { // CO20190520 - removed pointers for bools and doubles, added const where possible
    // determine anrl label, parameters, and parameter values of the input structure
    const uint setting = SG_SETTING_ANRL; // anrl setting choice is default
    return structure2anrl(xstr, tolerance, setting);
  }
} // namespace anrl

// ***************************************************************************
namespace anrl {
  string structure2anrl(xstructure& xstr, uint setting) { // DX20191031 - removed reference
    // determine anrl label, parameters, and parameter values of the input structure
    const double default_tolerance = SYM::defaultTolerance(xstr);
    return structure2anrl(xstr, default_tolerance, setting);
  }
} // namespace anrl

// ***************************************************************************
// anrl::structure2anrl() [MAIN]
// ***************************************************************************
namespace anrl {
  string structure2anrl(xstructure& xstr, double tolerance, uint input_setting, bool recalculate_symmetry, bool print_element_names, bool print_atomic_numbers) { // CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190829 - added recalculate_symmetry //DX20191031 - removed reference
    // determine anrl label, parameters, and parameter values of the input structure
    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);

    ostringstream oss;
    const ofstream FileMESSAGE;

    uint setting = input_setting; // DX20191230 - if symmetry already calculated, we want to store true setting

    // Calculate symmetry
    uint space_group_number = 0;
    if ((xstr.space_group_ITC == 0 && xstr.spacegroupnumber == 0) || recalculate_symmetry) { // DX20190829 - added if-statement; don't recalculate, it is faster
      space_group_number = xstr.SpaceGroup_ITC(tolerance, input_setting);
    } else if (xstr.space_group_ITC >= 1 && xstr.space_group_ITC <= 230) {
      space_group_number = xstr.space_group_ITC;
      setting = xstr.setting_ITC; // DX20191230
    } else if (xstr.spacegroupnumber >= 1 && xstr.spacegroupnumber <= 230) {
      space_group_number = xstr.spacegroupnumber;
      setting = xstr.spacegroupnumberoption; // DX20191230
    }

    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions);

    // ===== Determine ANRL label ===== //
    // stoichiometry
    xstructure tmp_xstr = xstr; // make new one so we do not override atom names
    ReScale(tmp_xstr, 1.0); // DX20191031
    const bool remove_ones = true;
    tmp_xstr.DecorateWithFakeElements(); // DX20200728 - fakeAtomNames() -> DecorateWithFakeElements()
    const vector<uint> reduced_composition = tmp_xstr.GetReducedComposition(false); // numerical sort is false
    const vector<string> fake_elements = tmp_xstr.GetElements();
    const string reduced_stoichiometry = pflow::prettyPrintCompound(fake_elements, reduced_composition, no_vrt, remove_ones, txt_ft); // remove ones is true

    // Pearson symbol (quick)
    uint conventional_cell_atoms_count = 0;
    for (size_t i = 0; i < xstr.wyckoff_sites_ITC.size(); i++) {
      conventional_cell_atoms_count += xstr.wyckoff_sites_ITC[i].multiplicity;
    }
    string lattice_type_and_centering = LATTICE::SpaceGroup2LatticeTypeAndCentering(space_group_number); // DX20191031
    char lattice_type = lattice_type_and_centering[0];
    const char lattice_centering = lattice_type_and_centering[1];

    // rhl fixes
    if (lattice_centering == 'R' && setting == SG_SETTING_2) {
      conventional_cell_atoms_count /= 3;
    } // for rhl, number in Pearson symbol is given wrt to rhl cell

    stringstream tmp;
    tmp << lattice_type << lattice_centering << conventional_cell_atoms_count;
    const string Pearson_symbol = tmp.str();

    // space group number
    const string space_group_number_str = aurostd::utype2string<uint>(space_group_number);

    // Wyckoff positions
    const string Wyckoff_string = anrl::groupedWyckoffPosition2ANRLString(grouped_Wyckoff_positions, true); // alphabetize=true

    // combine stoichiometry + Pearson + space group number + Wyckoff positions to make label
    const string aflow_label = reduced_stoichiometry + "_" + Pearson_symbol + "_" + space_group_number_str + "_" + Wyckoff_string;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << ":: AFLOW ANRL label = " << aflow_label << endl;
    }

    // ===== Determine parameters ===== //
    vector<string> parameter_list;
    vector<string> lattice_parameter_list;
    vector<string> Wyckoff_parameter_list;
    vector<double> parameter_values;
    vector<double> lattice_parameter_values;
    vector<double> Wyckoff_parameter_values;

    // lattice parameters
    lattice_parameter_list = anrl::getANRLLatticeParameterString(lattice_type);
    if (!xstr.wyccar_ITC.empty()) {
      lattice_parameter_values = anrl::getANRLLatticeParameterValuesFromWyccar(xstr.wyccar_ITC, lattice_type, lattice_centering, setting);
    } else {
      lattice_parameter_values = anrl::getANRLLatticeParameterValuesFromABCAngles(xstr, lattice_type, lattice_centering, setting);
    }

    // Wyckoff parameters
    vector<wyckoffsite_ITC> ordered_Wyckoff_sites_ITC = xstr.wyckoff_sites_ITC;
    // reorder Wyckoff positions alphabetically by Wyckoff letter, then by species
    std::sort(ordered_Wyckoff_sites_ITC.begin(), ordered_Wyckoff_sites_ITC.end());

    if (LDEBUG) {
      for (size_t i = 0; i < ordered_Wyckoff_sites_ITC.size(); i++) {
        cerr << __AFLOW_FUNC__ << ":Ordered Wyckoff site: " << ordered_Wyckoff_sites_ITC[i] << endl;
      }
    }

    // determine degrees of freedom in Wyckoff positions
    for (size_t i = 0; i < ordered_Wyckoff_sites_ITC.size(); i++) {
      bool contains_x = false;
      bool contains_y = false;
      bool contains_z = false;
      if (!ordered_Wyckoff_sites_ITC[i].equations.empty()) {
        for (size_t j = 0; j < ordered_Wyckoff_sites_ITC[i].equations[0].size(); j++) { // DX20190311 - used ordered Wyckoff variable
          if (ordered_Wyckoff_sites_ITC[i].equations[0][j].find("x") != std::string::npos) {
            contains_x = true;
          } // DX20190311 - used ordered Wyckoff variable
          if (ordered_Wyckoff_sites_ITC[i].equations[0][j].find("y") != std::string::npos) {
            contains_y = true;
          } // DX20190311 - used ordered Wyckoff variable
          if (ordered_Wyckoff_sites_ITC[i].equations[0][j].find("z") != std::string::npos) {
            contains_z = true;
          } // DX20190311 - used ordered Wyckoff variable
        }
        // store
        string variable_designation;
        if (contains_x) {
          variable_designation = "x" + aurostd::utype2string<uint>(i + 1);
          Wyckoff_parameter_list.push_back(variable_designation);
          Wyckoff_parameter_values.push_back(ordered_Wyckoff_sites_ITC[i].coord(1));
        }
        if (contains_y) {
          variable_designation = "y" + aurostd::utype2string<uint>(i + 1);
          Wyckoff_parameter_list.push_back(variable_designation);
          Wyckoff_parameter_values.push_back(ordered_Wyckoff_sites_ITC[i].coord(2));
        }
        if (contains_z) {
          variable_designation = "z" + aurostd::utype2string<uint>(i + 1);
          Wyckoff_parameter_list.push_back(variable_designation);
          Wyckoff_parameter_values.push_back(ordered_Wyckoff_sites_ITC[i].coord(3));
        }
      } else {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "The equations for site " + aurostd::utype2string(i) + "are not provided. Check symmetry", _INPUT_MISSING_); // CO20200624
      }
    }

    // combine parameter vectors
    parameter_list = lattice_parameter_list;
    parameter_list.insert(parameter_list.end(), Wyckoff_parameter_list.begin(), Wyckoff_parameter_list.end());

    parameter_values = lattice_parameter_values;
    parameter_values.insert(parameter_values.end(), Wyckoff_parameter_values.begin(), Wyckoff_parameter_values.end());

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << "::ANRL parameters:" << endl;
      for (size_t i = 0; i < parameter_list.size(); i++) {
        cerr << parameter_list[i] << "=" << parameter_values[i] << endl;
      }
    }

    // store label/params/params values/etc. in xstructure
    xstr.prototype = aflow_label;
    xstr.prototype_parameter_list = parameter_list;
    xstr.prototype_parameter_values = parameter_values;
    xstr.num_parameters = parameter_list.size();
    xstr.num_lattice_parameters = lattice_parameter_list.size();

    // print label/params/params values
    string format = "text";
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      format = "text";
    } else if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      format = "json";
    }

    // store atomic number //DX20210622
    vector<uint> atomic_numbers;
    if (print_atomic_numbers && pflow::hasRealElements(xstr)) {
      for (size_t i = 0; i < xstr.species.size(); i++) {
        atomic_numbers.push_back(xelement::symbol2Z(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])));
      }
    }

    if (format == "json") {
      const string eendl;
      const bool roff = true; // round off
      stringstream sscontent_json;
      vector<string> vcontent_json;

      sscontent_json << R"("aflow_prototype_label":")" << aflow_label << "\"" << eendl;
      vcontent_json.push_back(sscontent_json.str());
      sscontent_json.str("");
      sscontent_json << "\"aflow_prototype_params_list\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(parameter_list, "\""), ",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str());
      sscontent_json.str("");
      sscontent_json << "\"aflow_prototype_params_values\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(parameter_values, 8, roff), ",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str());
      sscontent_json.str("");
      if (print_element_names) { // DX20210622
        sscontent_json << "\"element_names\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(xstr.species, "\""), ",") << "]" << eendl;
        vcontent_json.push_back(sscontent_json.str());
        sscontent_json.str("");
      }
      if (print_atomic_numbers) { // DX20210622
        sscontent_json << "\"atomic_numbers\":[" << aurostd::joinWDelimiter(atomic_numbers, ",") << "]" << eendl;
        vcontent_json.push_back(sscontent_json.str());
        sscontent_json.str("");
      }

      oss << "{" << aurostd::joinWDelimiter(vcontent_json, ",") << "}";
    } else {
      oss << "AFLOW label    : " << aflow_label << endl;
      oss << "params         : " << aurostd::joinWDelimiter(parameter_list, ",") << endl;
      oss << "params values  : " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(parameter_values, 6), ",") << endl;
      if (print_element_names) {
        oss << "element names  : " << aurostd::joinWDelimiter(xstr.species, ",") << endl;
      } // DX20210622
      if (print_atomic_numbers) {
        oss << "atomic numbers : " << aurostd::joinWDelimiter(atomic_numbers, ",") << endl;
      } // DX20210622
    }

    return oss.str();
  }
} // namespace anrl
// DX20190208 - for making ANRL label - END

// ***************************************************************************
// anrl::getLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getLattice(const string& lattice_and_centering, const char& space_group_letter, const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the lattice in the ANRL convention and populates with the
    // relevant lattice parameters.
    // space_group_letter : needed to differentiate between the A, and C
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " Lattice mode=" << mode << endl;
    }

    // ---------------------------------------------------------------------------
    // triclinic crystal system
    if (lattice_and_centering == "aP") {
      return getTriclinicLattice(lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    // ---------------------------------------------------------------------------
    // monoclinic crystal system
    else if (lattice_and_centering == "mP" || lattice_and_centering == "mC") {
      return getMonoclinicLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    // ---------------------------------------------------------------------------
    // orthorhombic crystal system
    else if (lattice_and_centering == "oP" || lattice_and_centering == "oC" || lattice_and_centering == "oI" || lattice_and_centering == "oF") {
      return getOrthorhombicLattice(lattice_and_centering, space_group_letter, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    // ---------------------------------------------------------------------------
    // tetragonal crystal system
    else if (lattice_and_centering == "tP" || lattice_and_centering == "tI") {
      return getTetragonalLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    // ---------------------------------------------------------------------------
    // hexagonal crystal system (includes hex + rhl)
    else if (lattice_and_centering == "hP" || lattice_and_centering == "hR") {
      return getHexagonalLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    // ---------------------------------------------------------------------------
    // cubic crystal system
    else if (lattice_and_centering == "cP" || lattice_and_centering == "cF" || lattice_and_centering == "cI") {
      return getCubicLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    else {
      message << "Lattice type and centering are not possible: " << lattice_and_centering;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ILLEGAL_);
    }
  }
} // namespace anrl

// ***************************************************************************
// anrl::getTriclinicLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getTriclinicLattice(const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the triclinic lattice in the ANRL convention and populates
    // with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice
    // NOTE : the lattice_and_centering input is not required; triclinic
    //        systems only have one centering option (P)

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs
    if (lattice_parameter_values.size() != 6) {
      message << "There needs to be 6 lattice parameters to build the triclinic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }

    // ---------------------------------------------------------------------------
    // main variables
    const xmatrix<double> lattice;
    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    uint i = 0;
    const double a = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a=" << a << endl;
    }
    const double bovera = lattice_parameter_values[i++];
    const double b = bovera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " b=" << b << " (b/a=" << bovera << ")" << endl;
    }
    const double covera = lattice_parameter_values[i++];
    const double c = covera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " c=" << c << " (c/a=" << covera << ")" << endl;
    }
    const double alpha = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " alpha=" << alpha << endl;
    }
    const double beta = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " beta=" << beta << endl;
    }
    const double gamma = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " gamma=" << gamma << endl;
    }

    // ---------------------------------------------------------------------------
    // triclinic - tri (aP)
    if (mode == 0 || mode == 1) { // primitive and conventional cells are the same

      const double cx = c * cos(deg2rad * beta);
      const double cy = c * (cos(deg2rad * alpha) - cos(deg2rad * beta) * cos(deg2rad * gamma)) / sin(deg2rad * gamma);
      const double cz = sqrt(pow(c, 2.0) - pow(cx, 2.0) - pow(cy, 2.0));

      a1 = a * xn;
      a2 = b * cos(deg2rad * gamma) * xn + b * sin(deg2rad * gamma) * yn;
      a3 = cx * xn + cy * yn + cz * zn;

      // ---------------------------------------------------------------------------
      // build lattice
      lattice(1, 1) = a1(1);
      lattice(1, 2) = a1(2);
      lattice(1, 3) = a1(3);
      lattice(2, 1) = a2(1);
      lattice(2, 2) = a2(2);
      lattice(2, 3) = a2(3);
      lattice(3, 1) = a3(1);
      lattice(3, 2) = a3(2);
      lattice(3, 3) = a3(3);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice = " << lattice << endl;
    }

    return lattice;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getMonoclinicLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getMonoclinicLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the monoclinic lattice in the ANRL convention and populates
    // with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs
    if (lattice_parameter_values.size() != 4) {
      message << "There needs to be 4 lattice parameters to build the monoclinic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }

    // ---------------------------------------------------------------------------
    // main variables
    const xmatrix<double> lattice;
    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    uint i = 0;
    const double a = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a=" << a << endl;
    }
    const double bovera = lattice_parameter_values[i++];
    const double b = bovera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " b=" << b << " (b/a=" << bovera << ")" << endl;
    }
    const double covera = lattice_parameter_values[i++];
    const double c = covera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " c=" << c << " (c/a=" << covera << ")" << endl;
    }
    const double beta = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " beta=" << beta << endl;
    }

    // ---------------------------------------------------------------------------
    // primitive lattices
    if (mode == 0) {
      // ---------------------------------------------------------------------------
      // simple monoclinic - mcl (mP)
      if (lattice_and_centering == "mP") {
        a1 = a * xn;
        a2 = b * yn;
        a3 = c * cos(deg2rad * beta) * xn + c * sin(deg2rad * beta) * zn;
      }

      // ---------------------------------------------------------------------------
      // base-centered monoclinic - mclc (mC)
      if (lattice_and_centering == "mC") {
        a1 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * b * yn;
        a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
        a3 = c * cos(deg2rad * beta) * xn + c * sin(deg2rad * beta) * zn;
      }
    }

    // ---------------------------------------------------------------------------
    // conventional lattice
    else if (mode == 1) {
      a1 = a * xn;
      a2 = b * yn;
      a3 = c * cos(deg2rad * beta) * xn + c * sin(deg2rad * beta) * zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice
    lattice(1, 1) = a1(1);
    lattice(1, 2) = a1(2);
    lattice(1, 3) = a1(3);
    lattice(2, 1) = a2(1);
    lattice(2, 2) = a2(2);
    lattice(2, 3) = a2(3);
    lattice(3, 1) = a3(1);
    lattice(3, 2) = a3(2);
    lattice(3, 3) = a3(3);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice = " << lattice << endl;
    }

    return lattice;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getOrthorhombicLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getOrthorhombicLattice(const string& lattice_and_centering, const char& space_group_letter, const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the orthorhombic lattice in the ANRL convention and populates
    // with the relevant lattice parameters.
    // space_group_letter : needed to differentiate between the A, and C
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs
    if (lattice_parameter_values.size() != 3) {
      message << "There needs to be 3 lattice parameters to build the orthorhombic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }

    // ---------------------------------------------------------------------------
    // main variables
    const xmatrix<double> lattice;
    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    uint i = 0;
    const double a = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a=" << a << endl;
    }
    const double bovera = lattice_parameter_values[i++];
    const double b = bovera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " b=" << b << " (b/a=" << bovera << ")" << endl;
    }
    const double covera = lattice_parameter_values[i++];
    const double c = covera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " c=" << c << " (c/a=" << covera << ")" << endl;
    }

    // ---------------------------------------------------------------------------
    // primitive lattices
    if (mode == 0) {
      // ---------------------------------------------------------------------------
      // simple orthorhombic - orc (oP)
      if (lattice_and_centering == "oP") {
        a1 = a * xn;
        a2 = b * yn;
        a3 = c * zn;
      }

      // ---------------------------------------------------------------------------
      // base-centered orthorhombic - orcc (oC)
      if (lattice_and_centering == "oC") {
        // A-centered
        if (space_group_letter == 'A') {
          a1 = a * xn;
          a2 = (1.0 / 2.0) * b * yn - (1.0 / 2.0) * c * zn;
          a3 = (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
        }
        // C-centered
        if (space_group_letter == 'C') {
          a1 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * b * yn;
          a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
          a3 = c * zn;
        }
      }

      // ---------------------------------------------------------------------------
      // body-centered orthorhombic - orci (oI)
      if (lattice_and_centering == "oI") {
        a1 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
        a2 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
        a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn - (1.0 / 2.0) * c * zn;
      }

      // ---------------------------------------------------------------------------
      // face-centered orthorhombic - orcf (oF)
      if (lattice_and_centering == "oF") {
        a1 = (1.0 / 2.0) * b * yn + (1.0 / 2.0) * c * zn;
        a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * c * zn;
        a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * b * yn;
      }
    }

    // ---------------------------------------------------------------------------
    // conventional lattice
    else if (mode == 1) {
      a1 = a * xn;
      a2 = b * yn;
      a3 = c * zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice
    lattice(1, 1) = a1(1);
    lattice(1, 2) = a1(2);
    lattice(1, 3) = a1(3);
    lattice(2, 1) = a2(1);
    lattice(2, 2) = a2(2);
    lattice(2, 3) = a2(3);
    lattice(3, 1) = a3(1);
    lattice(3, 2) = a3(2);
    lattice(3, 3) = a3(3);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice = " << lattice << endl;
    }

    return lattice;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getTetragonaLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getTetragonalLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the tetragonal lattice in the ANRL convention and populates
    // with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs
    if (lattice_parameter_values.size() != 2) {
      message << "There needs to be 2 lattice parameters to build the tetragonal lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }

    // ---------------------------------------------------------------------------
    // main variables
    const xmatrix<double> lattice;
    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    uint i = 0;
    const double a = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a=" << a << endl;
    }
    const double covera = lattice_parameter_values[i++];
    const double c = covera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " c=" << c << " (c/a=" << covera << ")" << endl;
    }

    // ---------------------------------------------------------------------------
    // primitive lattices
    if (mode == 0) {
      // ---------------------------------------------------------------------------
      // simple tetragonal - tet (tP)
      if (lattice_and_centering == "tP") {
        a1 = a * xn;
        a2 = a * yn;
        a3 = c * zn;
      }

      // ---------------------------------------------------------------------------
      // body-centered tegtragonal - bct (tI)
      if (lattice_and_centering == "tI") {
        a1 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn + (1.0 / 2.0) * c * zn;
        a2 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * a * yn + (1.0 / 2.0) * c * zn;
        a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn - (1.0 / 2.0) * c * zn;
      }
    }
    // ---------------------------------------------------------------------------
    // conventional lattice
    else if (mode == 1) {
      a1 = a * xn;
      a2 = a * yn;
      a3 = c * zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice
    lattice(1, 1) = a1(1);
    lattice(1, 2) = a1(2);
    lattice(1, 3) = a1(3);
    lattice(2, 1) = a2(1);
    lattice(2, 2) = a2(2);
    lattice(2, 3) = a2(3);
    lattice(3, 1) = a3(1);
    lattice(3, 2) = a3(2);
    lattice(3, 3) = a3(3);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice = " << lattice << endl;
    }

    return lattice;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getHexagonalLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getHexagonalLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the hexagonal/rhombohedral lattice in the ANRL convention and
    // populates with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check number of inputs
    if (lattice_parameter_values.size() != 2) {
      message << "There needs to be 2 lattice parameters to build the hexagonal/rhombohedral lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }

    // ---------------------------------------------------------------------------
    // main variables
    const xmatrix<double> lattice;
    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    uint i = 0;
    const double a = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a=" << a << endl;
    }
    const double covera = lattice_parameter_values[i++];
    const double c = covera * a;
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " c=" << c << " (c/a=" << covera << ")" << endl;
    }

    // ---------------------------------------------------------------------------
    // primitive lattices
    if (mode == 0) {
      // ---------------------------------------------------------------------------
      // hexagonal - hex (hP)
      if (lattice_and_centering == "hP") {
        a1 = (1.0 / 2.0) * a * xn - (sqrt(3.0) / 2.0) * a * yn;
        a2 = (1.0 / 2.0) * a * xn + (sqrt(3.0) / 2.0) * a * yn;
        a3 = c * zn;
      }

      // ---------------------------------------------------------------------------
      // rhombohedral - rhl (hR)
      if (lattice_and_centering == "hR") {
        a1 = (1.0 / 2.0) * a * xn - (1.0 / (2.0 * sqrt(3.0))) * a * yn + (1.0 / 3.0) * c * zn;
        a2 = (1.0 / sqrt(3.0)) * a * yn + (1.0 / 3.0) * c * zn;
        a3 = -(1.0 / 2.0) * a * xn - (1.0 / (2.0 * sqrt(3.0))) * a * yn + (1.0 / 3.0) * c * zn;
      }
    }
    // ---------------------------------------------------------------------------
    // conventional lattice
    else if (mode == 1) {
      a1 = (1.0 / 2.0) * a * xn - (sqrt(3.0) / 2.0) * a * yn;
      a2 = (1.0 / 2.0) * a * xn + (sqrt(3.0) / 2.0) * a * yn;
      a3 = c * zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice
    lattice(1, 1) = a1(1);
    lattice(1, 2) = a1(2);
    lattice(1, 3) = a1(3);
    lattice(2, 1) = a2(1);
    lattice(2, 2) = a2(2);
    lattice(2, 3) = a2(3);
    lattice(3, 1) = a3(1);
    lattice(3, 2) = a3(2);
    lattice(3, 3) = a3(3);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice = " << lattice << endl;
    }

    return lattice;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getCubicLattice()
// ***************************************************************************
namespace anrl {
  xmatrix<double> getCubicLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode) {
    // Returns the cubic lattice in the ANRL convention and
    // populates with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // check number of inputs
    if (lattice_parameter_values.size() != 1) {
      message << "There needs to be 1 lattice parameters to build the cubic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
    }

    // ---------------------------------------------------------------------------
    // main variables
    const xmatrix<double> lattice;
    const xvector<double> xn(3);
    xn(1) = 1.0;
    xn(2) = 0.0;
    xn(3) = 0.0;
    const xvector<double> yn(3);
    yn(1) = 0.0;
    yn(2) = 1.0;
    yn(3) = 0.0;
    const xvector<double> zn(3);
    zn(1) = 0.0;
    zn(2) = 0.0;
    zn(3) = 1.0;
    xvector<double> a1(3);
    xvector<double> a2(3);
    xvector<double> a3(3);

    uint i = 0;
    const double a = lattice_parameter_values[i++];
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " a=" << a << endl;
    }

    // ---------------------------------------------------------------------------
    // primitive lattices
    if (mode == 0) {
      // ---------------------------------------------------------------------------
      // simple cubic - cub (cP)
      if (lattice_and_centering == "cP") {
        a1 = a * xn;
        a2 = a * yn;
        a3 = a * zn;
      }

      // ---------------------------------------------------------------------------
      // face-centered cubic - fcc (cF)
      if (lattice_and_centering == "cF") {
        a1 = (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
        a2 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * zn;
        a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn;
      }

      // ---------------------------------------------------------------------------
      // body-centered cubic - bcc (cI)
      if (lattice_and_centering == "cI") {
        a1 = -(1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
        a2 = (1.0 / 2.0) * a * xn - (1.0 / 2.0) * a * yn + (1.0 / 2.0) * a * zn;
        a3 = (1.0 / 2.0) * a * xn + (1.0 / 2.0) * a * yn - (1.0 / 2.0) * a * zn;
      }
    }

    // ---------------------------------------------------------------------------
    // conventional lattice
    else if (mode == 1) {
      a1 = a * xn;
      a2 = a * yn;
      a3 = a * zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice
    lattice(1, 1) = a1(1);
    lattice(1, 2) = a1(2);
    lattice(1, 3) = a1(3);
    lattice(2, 1) = a2(1);
    lattice(2, 2) = a2(2);
    lattice(2, 3) = a2(3);
    lattice(3, 1) = a3(1);
    lattice(3, 2) = a3(2);
    lattice(3, 3) = a3(3);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " lattice = " << lattice << endl;
    }

    return lattice;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getAtomsFromWyckoff()
// ***************************************************************************
namespace anrl {
  deque<_atom> getAtomsFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_positions, const xmatrix<double>& lattice_conventional) {
    // create atoms (deque<_atom>) from the Wyckoff positions by
    // plugging in the parameter values into the Wyckoff equations

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // variables
    deque<_atom> atoms_conventional_cell;
    _atom atom_tmp;

    // ---------------------------------------------------------------------------
    // create f2c once (efficiency)
    const xmatrix<double> f2c = trasp(lattice_conventional);

    // ---------------------------------------------------------------------------
    // create an _atom for each Wyckoff position by plugging in the relevant parameters
    for (size_t i = 0; i < Wyckoff_positions.size(); i++) {
      // get x, y, and z coordinates from the Wyckoff object
      // added format=FIXED_STREAM since SYM::simplify has trouble with scientific notation
      const string x_value_string = aurostd::utype2string<double>(Wyckoff_positions[i].coord(1), AUROSTD_DEFAULT_PRECISION, FIXED_STREAM); // DX20201028 - added precision and format
      const string y_value_string = aurostd::utype2string<double>(Wyckoff_positions[i].coord(2), AUROSTD_DEFAULT_PRECISION, FIXED_STREAM); // DX20201028 - added precision and format
      const string z_value_string = aurostd::utype2string<double>(Wyckoff_positions[i].coord(3), AUROSTD_DEFAULT_PRECISION, FIXED_STREAM); // DX20201028 - added precision and format
      for (size_t j = 0; j < Wyckoff_positions[i].equations.size(); j++) {
        vector<string> coordinate_vstring = Wyckoff_positions[i].equations[j];
        const xvector<double> coordinate(3);
        for (size_t k = 0; k < coordinate_vstring.size(); k++) {
          // substitute variable with value
          aurostd::StringSubstInPlace(coordinate_vstring[k], "x", x_value_string);
          aurostd::StringSubstInPlace(coordinate_vstring[k], "y", y_value_string);
          aurostd::StringSubstInPlace(coordinate_vstring[k], "z", z_value_string);
          // simplify string
          const vector<SYM::sdouble> component_tmp = SYM::simplify(coordinate_vstring[k]);
          const string component_string = SYM::formatWyckoffPosition(component_tmp, false);
          // ensure the string is numeric
          if (!aurostd::isfloat(component_string)) {
            message << "There are non-numeric characters in the string after variable-value substitution: component " << k << "=" << component_string << endl;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _VALUE_ERROR_);
          }
          coordinate(k + 1) = aurostd::string2utype<double>(component_string);
        }
        // store atom info
        atom_tmp.fpos = coordinate; // coordinate from ITC is fpos
        atom_tmp.cpos = f2c * coordinate;
        atom_tmp.type = Wyckoff_positions[i].index;
        atom_tmp.name = Wyckoff_positions[i].type;
        atoms_conventional_cell.push_back(atom_tmp);
      }
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " atoms in the conventional cell:" << endl;
      for (size_t i = 0; i < atoms_conventional_cell.size(); i++) {
        cerr << "atoms: " << atoms_conventional_cell[i] << " " << atoms_conventional_cell[i].name << endl;
      }
    }
    return atoms_conventional_cell;
  }
} // namespace anrl

// ***************************************************************************
// anrl::determineWyckoffVariables()
// ***************************************************************************
namespace anrl {
  vector<string> determineWyckoffVariables(vector<wyckoffsite_ITC>& Wyckoff_positions) {
    // Determines the variable coordinates in the Wyckoff positions and
    // returns a vector of Wyckoff variables (x, y, or z) that need to be
    // specified. The Wyckoff positions should be ordered by Wyckoff letter
    // (alphabetic) in accordance with the ANRL convention.
    // The format is : x1, y1, z1, x2, y2, z2, x3, ...
    // Note: Wyckoff_positions are updated, i.e., the corresponding parameter
    // index is assigned so we can substitute values later

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    // ---------------------------------------------------------------------------
    // variables
    vector<string> Wyckoff_parameter_list;

    for (size_t i = 0; i < Wyckoff_positions.size(); i++) {
      bool contains_x = false;
      bool contains_y = false;
      bool contains_z = false;
      if (!Wyckoff_positions[i].equations.empty()) {
        // ---------------------------------------------------------------------------
        // look at the representative Wyckoff position (index 0) and see if it
        // x, y, and/or z
        for (size_t j = 0; j < Wyckoff_positions[i].equations[0].size(); j++) {
          if (Wyckoff_positions[i].equations[0][j].find("x") != std::string::npos) {
            contains_x = true;
          }
          if (Wyckoff_positions[i].equations[0][j].find("y") != std::string::npos) {
            contains_y = true;
          }
          if (Wyckoff_positions[i].equations[0][j].find("z") != std::string::npos) {
            contains_z = true;
          }
        }

        // ---------------------------------------------------------------------------
        // check for x-coordinate
        if (contains_x) {
          const string variable = "x";
          const string variable_name = variable + aurostd::utype2string<uint>(i + 1);
          Wyckoff_parameter_list.push_back(variable_name);
          Wyckoff_positions[i].parameter_index = i + 1;
        }
        // ---------------------------------------------------------------------------
        // check for y-coordinate
        if (contains_y) {
          const string variable = "y";
          const string variable_name = variable + aurostd::utype2string<uint>(i + 1);
          Wyckoff_parameter_list.push_back(variable_name);
          Wyckoff_positions[i].parameter_index = i + 1;
        }
        // ---------------------------------------------------------------------------
        // check for z-coordinate
        if (contains_z) {
          const string variable = "z";
          const string variable_name = variable + aurostd::utype2string<uint>(i + 1);
          Wyckoff_parameter_list.push_back(variable_name);
          Wyckoff_positions[i].parameter_index = i + 1;
        }
      } else {
        message << "The equations for site " << i << "are not provided. Check symmetry.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " parameters=" << aurostd::joinWDelimiter(Wyckoff_parameter_list, ",") << endl;
    }

    return Wyckoff_parameter_list;
  }
} // namespace anrl

// ***************************************************************************
// anrl::applyWyckoffValues()
// ***************************************************************************
namespace anrl {
  void applyWyckoffValues(const vector<double>& Wyckoff_parameter_values, vector<wyckoffsite_ITC>& Wyckoff_positions) {
    // Applies the Wyckoff parameter values to the Wyckoff object
    // i.e., the coord attribute, indicating the degree of freedom (x, y, or z)
    // The Wyckoff positions should be ordered by Wyckoff letter
    // (alphabetic) in accordance with the ANRL convention.
    // The format is : x1, y1, z1, x2, y2, z2, x3, ...

    stringstream message;

    // ---------------------------------------------------------------------------
    // variables
    const vector<string> Wyckoff_parameter_list;
    uint w = 0; // Wyckoff counter

    for (size_t i = 0; i < Wyckoff_positions.size(); i++) {
      bool contains_x = false;
      bool contains_y = false;
      bool contains_z = false;
      if (!Wyckoff_positions[i].equations.empty()) {
        // ---------------------------------------------------------------------------
        // look at the representative Wyckoff position (index 0) and see if it
        // x, y, and/or z
        for (size_t j = 0; j < Wyckoff_positions[i].equations[0].size(); j++) {
          if (Wyckoff_positions[i].equations[0][j].find("x") != std::string::npos) {
            contains_x = true;
          }
          if (Wyckoff_positions[i].equations[0][j].find("y") != std::string::npos) {
            contains_y = true;
          }
          if (Wyckoff_positions[i].equations[0][j].find("z") != std::string::npos) {
            contains_z = true;
          }
        }

        // ---------------------------------------------------------------------------
        // check for x-coordinate
        if (contains_x) {
          // store parameter value
          if (w < Wyckoff_parameter_values.size()) {
            Wyckoff_positions[i].coord(1) = Wyckoff_parameter_values[w++];
          } else {
            message << "There are too few input parameters; could not populate the x-coordinate for Wyckoff position " << i;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
          }
        }
        // ---------------------------------------------------------------------------
        // check for y-coordinate
        if (contains_y) {
          // store parameter value
          if (w < Wyckoff_parameter_values.size()) {
            Wyckoff_positions[i].coord(2) = Wyckoff_parameter_values[w++];
          } else {
            message << "There are too few input parameters; could not populate the y-coordinate for Wyckoff position " << i;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
          }
        }
        // ---------------------------------------------------------------------------
        // check for z-coordinate
        if (contains_z) {
          // store parameter value
          if (w < Wyckoff_parameter_values.size()) {
            Wyckoff_positions[i].coord(3) = Wyckoff_parameter_values[w++];
          } else {
            message << "There are too few input parameters; could not populate the z-coordinate for Wyckoff position " << i;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
          }
        }
      } else {
        message << "The equations for site " << i << "are not provided. Check symmetry.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
  }
} // namespace anrl

// ***************************************************************************
// anrl::containsDuplicateWyckoffCoordinate()
// ***************************************************************************
namespace anrl {
  bool containsDuplicateWyckoffCoordinate(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered) {
    // Checks if Wyckoff positions occur in the structure multiple times
    // Two cases:
    //   1) multiple instances of a fixed Wyckoff position (i.e., no variables)
    //   2) multiple instances of a variable Wyckoff position with the same
    //      parameters

    stringstream message;

    // ---------------------------------------------------------------------------
    // tolerance indicating if Wyckoff coordinates are the same (heuristic)
    const double _WYCKOFF_FRACTIONAL_TOL_ = DEFAULT_ANRL_WYCKOFF_FRACTIONAL_TOL; // default = 1e-6

    vector<wyckoffsite_ITC> ordered_Wyckoff_sites = wyckoff_sites_ITC;
    if (!already_ordered) {
      std::sort(ordered_Wyckoff_sites.begin(), ordered_Wyckoff_sites.end(), sortWyckoffByLetter);
    }

    for (size_t i = 0; i < ordered_Wyckoff_sites.size(); i++) {
      for (size_t j = i + 1; j < ordered_Wyckoff_sites.size(); j++) {
        if (ordered_Wyckoff_sites[i].letter == ordered_Wyckoff_sites[j].letter) {
          // ---------------------------------------------------------------------------
          // case 1: no variables in representative Wyckoff positions (first equation)
          // means we should only have one instance of this Wyckoff position
          bool contains_variable = false;
          for (size_t k = 0; k < ordered_Wyckoff_sites[i].equations[0].size(); k++) {
            if (ordered_Wyckoff_sites[i].equations[0][k].find("x") != std::string::npos) {
              contains_variable = true;
              break;
            }
            if (ordered_Wyckoff_sites[i].equations[0][k].find("y") != std::string::npos) {
              contains_variable = true;
              break;
            }
            if (ordered_Wyckoff_sites[i].equations[0][k].find("z") != std::string::npos) {
              contains_variable = true;
              break;
            }
          }
          if (!contains_variable) {
            message << "Contains multiple static (i.e., no variable) Wyckoff positions: " << ordered_Wyckoff_sites[i].letter;
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
            return true;
          }
          // ---------------------------------------------------------------------------
          // case 2: contains variables, but values are the same; this is a quick check
          // (i.e., use fractional tol), more rigorous check done later
          else if (aurostd::isequal(ordered_Wyckoff_sites[i].coord(1), ordered_Wyckoff_sites[j].coord(1), _WYCKOFF_FRACTIONAL_TOL_) &&
                   aurostd::isequal(ordered_Wyckoff_sites[i].coord(2), ordered_Wyckoff_sites[j].coord(2), _WYCKOFF_FRACTIONAL_TOL_) &&
                   aurostd::isequal(ordered_Wyckoff_sites[i].coord(3), ordered_Wyckoff_sites[j].coord(3), _WYCKOFF_FRACTIONAL_TOL_)) {
            message << "Contains duplicate Wyckoff letters with the same degrees of freedom: " << aurostd::joinWDelimiter(xvecDouble2vecString(ordered_Wyckoff_sites[i].coord), ",");
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
            return true;
          }
        }
        // ---------------------------------------------------------------------------
        // since ordered by Wyckoff letter, we can skip the rest (break) and start
        // where we left off in the second loop (i=j)
        else {
          i = j;
          break;
        }
      }
    }
    return false;
  }
} // namespace anrl

// ***************************************************************************
// anrl::getWyckoffSitesFromANRL()
// ***************************************************************************
namespace anrl {
  vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(const vector<string>& Wyckoff_tokens, const vector<string>& species, uint space_group_number, int setting) {
    // Get Wyckoff positions/sites from the ANRL Wyckoff designation
    // i.e., given Wyckoff letter and number of times they are used
    // (e.g., 4a, 2b, c) get the Wyckoff letter, multiplicity, site symmetry,
    // and equations

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);

    vector<wyckoffsite_ITC> wyckoff_sites_ITC;
    for (size_t i = 0; i < Wyckoff_tokens.size(); i++) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " Wyckoff designation=" << Wyckoff_tokens[i] << endl;
      }

      uint Wyckoff_multiplication_factor = 1;
      stringstream ss_Wyckoff_letter;
      stringstream ss_factor;

      for (size_t j = 0; j < Wyckoff_tokens[i].size(); j++) {
        // ---------------------------------------------------------------------------
        // extract the prefactor (if it exists; prefactor=1 is not usually given)
        if (isdigit(Wyckoff_tokens[i][j])) {
          ss_factor << Wyckoff_tokens[i][j];
          continue;
        }
        // ---------------------------------------------------------------------------
        // extract the Wyckoff letter
        else {
          ss_Wyckoff_letter << Wyckoff_tokens[i][j];
          if (!ss_factor.str().empty()) {
            Wyckoff_multiplication_factor = aurostd::string2utype<uint>(ss_factor.str());
          }
        }
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " " << Wyckoff_multiplication_factor << " Wyckoff position(s) with letter" << ss_Wyckoff_letter.str() << endl;
        }

        // ---------------------------------------------------------------------------
        // populates Wyckoff with corresponding letter, multiplicity, equations, etc.
        wyckoffsite_ITC Wyckoff_tmp;
        Wyckoff_tmp.getWyckoffFromLetter(space_group_number, ss_Wyckoff_letter.str(), setting);
        Wyckoff_tmp.index = i;
        Wyckoff_tmp.type = species[i];
        if (LDEBUG) {
          cerr << __AFLOW_FUNC__ << " extracted Wyckoff position:" << Wyckoff_tmp << endl;
        }

        // ---------------------------------------------------------------------------
        // store the Wyckoff position as indicated by the prefactor,
        // e.g., 4b-> four Wyckoff positions with letter b
        for (uint m = 0; m < Wyckoff_multiplication_factor; m++) {
          wyckoff_sites_ITC.push_back(Wyckoff_tmp);
        }
        // reset
        Wyckoff_multiplication_factor = 1;
        ss_Wyckoff_letter.str("");
        ss_factor.str("");
      }
    }

    return wyckoff_sites_ITC;
  }
} // namespace anrl

// ***************************************************************************
// anrl::specialCaseSymmetryTolerances
// ***************************************************************************
namespace anrl {
  double specialCaseSymmetryTolerances(const string& label_input) {
    // symmetry tolerances for specific prototypes
    // some parameter values can be "close" to higher symmetry points,
    // causing structures to fall into higher symmetries when analyzed with
    // certain tolerances; occurs for certain AFLOW Prototype Encyclopedia
    // structures

    // ---------------------------------------------------------------------------
    // A2B7C2_oF88_22_k_bdefghij_k-001 (Predicted Phase IV Cd2Re2O7)
    // see comments in http://aflow.org/prototype-encyclopedia/A2B7C2_oF88_22_k_bdefghij_k.html
    if (label_input == "A2B7C2_oF88_22_k_bdefghij_k-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // AB_oP8_33_a_ai-001 (Modderite)
    // see comments in http://aflow.org/prototype-encyclopedia/AB_oP8_33_a_a.html
    if (label_input == "AB_oP8_33_a_a-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // A2B_oC12_38_de_ab-001 (Au2V)
    // see comments in http://aflow.org/prototype-encyclopedia/A2B_oC12_38_de_ab.html
    else if (label_input == "A2B_oC12_38_de_ab-001") {
      return 0.0001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // AB4_oC20_41_a_2b-001 (PtSn4, Struk: D1_{c})
    // see comments in http://aflow.org/prototype-encyclopedia/AB4_oC20_41_a_2b.html
    else if (label_input == "AB4_oC20_41_a_2b-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // AB2_oC24_41_2a_2b-001 (PdSn2, Struk: C_{e})
    else if (label_input == "AB2_oC24_41_2a_2b-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // A5B2_oP14_49_dehq_ab-001 (beta-Ta2O5)
    else if (label_input == "A5B2_oP14_49_dehq_ab-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // A10B3C4_oP68_55_2e2fgh2i_adef_2e2f-001 (Orthorhombic Sr4Ru3O10, part 3)
    // see comments in http://aflow.org/prototype-encyclopedia/A10B3C4_oP68_55_2e2fgh2i_adef_2e2f.html
    else if (label_input == "A10B3C4_oP68_55_2e2fgh2i_adef_2e2f-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // AB_oC8_67_a_g-001 (alpha-FeSe)
    else if (label_input == "AB_oC8_67_a_g-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // AB_oC8_67_a_g-002 (alpha-PbO)
    else if (label_input == "AB_oC8_67_a_g-002") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // AB_tP8_111_n_n-001 (VN, low-temperature)
    // see comments in http://aflow.org/prototype-encyclopedia/AB_tP8_111_n_n.html
    else if (label_input == "AB_tP8_111_n_n-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // A4B6C_hP11_143_bd_2d_a-001 (ScRh6P4)
    // see comments in http://aflow.org/prototype-encyclopedia/A4B6C_hP11_143_bd_2d_a.html
    else if (label_input == "A4B6C_hP11_143_bd_2d_a-001") {
      return 0.001; // symmetry tolerance
    }
    // ---------------------------------------------------------------------------
    // A12BC4_cP34_195_2j_ab_2e-001 (PrRu4P12)
    // see comments in http://aflow.org/prototype-encyclopedia/A12BC4_cP34_195_2j_ab_2e.html
    else if (label_input == "A12BC4_cP34_195_2j_ab_2e-001") {
      return 0.001; // symmetry tolerance
    }
    return AUROSTD_MAX_DOUBLE;
  }
} // namespace anrl

// ***************************************************************************
// anrl::isSpecialCaseEquivalentPrototypes() //DX20210421
// ***************************************************************************
namespace anrl {
  bool isSpecialCaseEquivalentPrototypes(const vector<string>& labels_matched) {
    // Check if prototypes are expected to be duplicates of one another.
    // In "make check", we ensure newly added prototypes do not match with
    // existing ones.
    // However, some prototypes can match for the following reasons:
    //  1) prototypes are enantiomorphs (i.e., they are duplicates by
    //     construction, to help represent all 230 space groups)
    //  2) for historical reasons; due to structure refinement, unique
    //     Strukturbericht labeling, or other significance described in
    //     literature (these are usually explained in the comments of the
    //     prototype encyclopedia)
    //  3) improvements to XtalFinder reveal structures match (in general, the
    //     misfits will be just below the default threshold of 0.1; if they
    //     are not, then we have problems)
    // New prototype-matches should be investigated with XtalFinder and
    // only reported here if we wish to keep the equivalent prototypes.

    const uint nlabels_matched = labels_matched.size();

    // ---------------------------------------------------------------------------
    // list of 2 labels matching
    if (nlabels_matched == 2) {
      // ---------------------------------------------------------------------------
      // A2B_mP12_14_2e_e-001 (ZrO2, Baddeleyite, Struk: C43) == A2B_mP12_14_2e_e-009 (ZrO2, ICSD #659226)
      // misfit=0.0998
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      if (aurostd::WithinList(labels_matched, "A2B_mP12_14_2e_e-001") && aurostd::WithinList(labels_matched, "A2B_mP12_14_2e_e-009")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB3C4_oP16_31_a_ab_2ab-001 (AsCu3S4, Enargite, Struk:H2_5) == A3B4C_oP16_31_ab_2ab_a-001 (Li3O4V1, ICSD #19002)
      // misfit=0.0884
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "AB3C4_oP16_31_a_ab_2ab-001") && aurostd::WithinList(labels_matched, "A3B4C_oP16_31_ab_2ab_a-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB_oP8_62_c_c-002 (MnP, Struk:B31) == AB_oP8_62_c_c-005 (FeAs, Westerveldite, Struk:B14)
      // misfit=0.0442
      // REASON FOR DUPLICATE: historical; different Strukturbericht designations
      // see comments in http://aflow.org/prototype-encyclopedia/AB_oP8_62_c_c.FeAs.html
      else if (aurostd::WithinList(labels_matched, "AB_oP8_62_c_c-002") && aurostd::WithinList(labels_matched, "AB_oP8_62_c_c-005")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A_tP30_136_bf2ij-001 (beta-U, Struk:A_d) == sigma_tP30_136_bf2ij-001 (sigma-CrFe, Struk:D8_b)
      // misfit=0
      // REASON FOR DUPLICATE: historical; different Strukturbericht designations
      // see comments in http://aflow.org/prototype-encyclopedia/sigma_tP30_136_bf2ij.html
      else if (aurostd::WithinList(labels_matched, "A_tP30_136_bf2ij-001") && aurostd::WithinList(labels_matched, "sigma_tP30_136_bf2ij-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B_hP24_151_3c_2a-001 (CrCl3, Struk:D0_4) == A3B_hP24_153_3c_2b-001 (CrCl3, enantiomorph)
      // misfit=0
      // REASON FOR DUPLICATE: use enantiomorph to represent SG #153
      // see comments in http://aflow.org/prototype-encyclopedia/A3B_hP24_153_3c_2b.html
      else if (aurostd::WithinList(labels_matched, "A3B_hP24_151_3c_2a-001") && aurostd::WithinList(labels_matched, "A3B_hP24_153_3c_2b-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // ABC3_hR10_161_a_a_b-002 (Na1Nb1O3, ICSD #9645) == ABC3_hR10_161_a_a_b-004 (Ga1La1O3, ICSD #51036)
      // misfit=0.0963
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "ABC3_hR10_161_a_a_b-002") && aurostd::WithinList(labels_matched, "ABC3_hR10_161_a_a_b-004")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A2B3_hR10_167_c_e-001 (Al2O3, Corundum, Struk:D5_1) == A2B3_hR10_167_c_e-002 (O3V2 binary oxide)
      // misfit=0.0889
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A2B3_hR10_167_c_e-001") && aurostd::WithinList(labels_matched, "A2B3_hR10_167_c_e-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A2B_hP9_180_j_c-001 (beta-Quartz, Struk:C8) == A2B_hP9_181_j_c-001 (beta-SiO2, enantiomorph)
      // misfit=0
      // REASON FOR DUPLICATE: use enantiomorph to represent SG #181
      // see comments in http://aflow.org/prototype-encyclopedia/A2B_hP9_181_j_c.html
      else if (aurostd::WithinList(labels_matched, "A2B_hP9_180_j_c-001") && aurostd::WithinList(labels_matched, "A2B_hP9_181_j_c-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // ABC_tP24_91_d_d_d-001 (ThBC) == ABC_tP24_95_d_d_d-001 (ThBC, enantiomorph)
      // misfit=0
      // REASON FOR DUPLICATE: use enantiomorph to represent SG #95
      // see comments in http://aflow.org/prototype-encyclopedia/ABC_tP24_95_d_d_d.html
      else if (aurostd::WithinList(labels_matched, "ABC_tP24_91_d_d_d-001") && aurostd::WithinList(labels_matched, "ABC_tP24_95_d_d_d-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A2B3_hP30_169_2a_3a-001 (alpha-Al2S3) == A2B3_hP30_170_2a_3a-001 (alpha-Al2S3, enantiomorph)
      // misfit=0
      // REASON FOR DUPLICATE: use enantiomorph to represent SG #170
      // see comments in http://aflow.org/prototype-encyclopedia/A2B3_hP30_170_2a_3a.html
      else if (aurostd::WithinList(labels_matched, "A2B3_hP30_169_2a_3a-001") && aurostd::WithinList(labels_matched, "A2B3_hP30_170_2a_3a-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB3_hP24_178_b_ac-001 (AuF3) == AB3_hP24_179_b_ac-001 (AuF3)
      // misfit=0
      // REASON FOR DUPLICATE: use enantiomorph to represent SG #179
      // see comments in http://aflow.org/prototype-encyclopedia/AB3_hP24_179_b_ac.html
      else if (aurostd::WithinList(labels_matched, "AB3_hP24_178_b_ac-001") && aurostd::WithinList(labels_matched, "AB3_hP24_179_b_ac-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B_hP24_185_ab2c_c-001 (Cu3P) == AB3_hP24_185_c_ab2c-001 (Na3As, Struk:D0_18)
      // misfit=0.0753
      // REASON FOR DUPLICATE: historical, discrepancies in literature
      // see comments in http://aflow.org/prototype-encyclopedia/A3B_hP24_185_ab2c_c.html
      else if (aurostd::WithinList(labels_matched, "A3B_hP24_185_ab2c_c-001") && aurostd::WithinList(labels_matched, "AB3_hP24_185_c_ab2c-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB4C_mP12_13_f_2g_e-001 (MgO4W, ICSD #67903) == AB4C_mP12_13_f_2g_e-004 (CuO4W, ICSD #182751)
      // misfit=0.0822
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "AB4C_mP12_13_f_2g_e-001") && aurostd::WithinList(labels_matched, "AB4C_mP12_13_f_2g_e-004")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B_mC8_12_di_a-001 (N3Na1, ICSD #29370) == A3B_mC8_12_di_a-002 (N3Na1, ICSD #29376)
      // misfit=0.098383
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A3B_mC8_12_di_a-001") && aurostd::WithinList(labels_matched, "A3B_mC8_12_di_a-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B2_mC30_12_a4i_3i-001 (Ca3N2, ICSD #162794) == A3B2_mC30_12_a4i_3i-002 (Ca3N2, ICSD #169726)
      // misfit=0.0989
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A3B2_mC30_12_a4i_3i-001") && aurostd::WithinList(labels_matched, "A3B2_mC30_12_a4i_3i-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B2C2_mC14_12_ai_i_i-003 (C3Ho2Mo2, ICSD #88511) == A3B2C2_mC14_12_ai_i_i-004 (C3Ce2Mo2, ICSD #417827)
      // misfit=0.0921
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A3B2C2_mC14_12_ai_i_i-003") && aurostd::WithinList(labels_matched, "A3B2C2_mC14_12_ai_i_i-004")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A4B3C4_oI22_71_n_af_eh-001 (Ag4Dy3Sn4, ICSD #156968) == A3B4C4_oI22_71_af_eh_n-001 (La3Pd4Zn4, ICSD #182774)
      // misfit=0.0879
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A4B3C4_oI22_71_n_af_eh-001") && aurostd::WithinList(labels_matched, "A3B4C4_oI22_71_af_eh_n-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB5C2_tI32_140_a_cl_h-001 (Bi1Er5Pt2, ICSD #107217) == A2BC5_tI32_140_h_a_cl-001 (Au2Bi1Tb5, ICSD #156956)
      // misfit=0.0914
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "AB5C2_tI32_140_a_cl_h-001") && aurostd::WithinList(labels_matched, "A2BC5_tI32_140_h_a_cl-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A2BC_hR12_166_h_bc_ac-001 (Al2Cu1Yb1, ICSD #604213) == AB2C_hR12_166_bc_h_ac-001 (Ag1Al2Pr1, ICSD #604688)
      // misfit=0.0995
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A2BC_hR12_166_h_bc_ac-001") && aurostd::WithinList(labels_matched, "AB2C_hR12_166_bc_h_ac-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B2C_hP6_191_g_c_a-001 (Ag3Al2La1, ICSD #57329) == A2BC3_hP6_191_c_a_g-002 (Al2Ce1Pt3, ICSD #658142)
      // misfit=0.0989
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A3B2C_hP6_191_g_c_a-001") && aurostd::WithinList(labels_matched, "A2BC3_hP6_191_c_a_g-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A4B2C5_mC22_12_2i_i_aj-001 (B4La2Ni5, ICSD #63501) == A4B2C5_mC22_12_2i_i_aj-002 (B4La2Ni5, ICSD #170618)
      // misfit=0.0677
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A4B2C5_mC22_12_2i_i_aj-001") && aurostd::WithinList(labels_matched, "A4B2C5_mC22_12_2i_i_aj-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // ABC2_mC8_12_c_a_i-001 (metal-oxide; Na1Ni1O2, ICSD #26609) == ABC2_mC8_12_a_c_i-002 (Mn1Na1O2, ICSD #16270)
      // misfit=0.0629
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "ABC2_mC8_12_c_a_i-001") && aurostd::WithinList(labels_matched, "ABC2_mC8_12_a_c_i-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // ABC2_mC8_12_c_a_i-001 (metal-oxide; Na1Ni1O2, ICSD #26609) == AB2C_mC8_12_a_i_c-001 (Na1O2V1, ICSD #420138)
      // misfit=0.0741
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "ABC2_mC8_12_c_a_i-001") && aurostd::WithinList(labels_matched, "AB2C_mC8_12_a_i_c-001")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB2C_mC16_15_e_f_e-001 (Bi1O2Rb1, ICSD #407208) == ABC2_mC16_15_e_e_f-003 (Bi1Cs1O2, ICSD #406564)
      // misfit=0.0801
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "AB2C_mC16_15_e_f_e-001") && aurostd::WithinList(labels_matched, "ABC2_mC16_15_e_e_f-003")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // A3B9C4_hR16_146_3a_3b_4a-001 (Ba3O9Yb4, ICSD #33239) == A3B9C4_hR16_146_3a_3b_4a-002 (Ba3Ho4O9, ICSD #33807)
      // misfit=0.0999
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "A3B9C4_hR16_146_3a_3b_4a-001") && aurostd::WithinList(labels_matched, "A3B9C4_hR16_146_3a_3b_4a-002")) {
        return true;
      }
      // ---------------------------------------------------------------------------
      // AB4C_oC24_63_a_fg_c-001 (MgSO4, anrl part 1) == ABC4_oC24_63_c_a_fg-002 (Cr1Hg1O4, ICSD #416147)
      // misfit=0.0924
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      else if (aurostd::WithinList(labels_matched, "AB4C_oC24_63_a_fg_c-001") && aurostd::WithinList(labels_matched, "ABC4_oC24_63_c_a_fg-002")) {
        return true;
      }
    }
    // ---------------------------------------------------------------------------
    // list of 3 labels matching
    else if (nlabels_matched == 3) {
      // ---------------------------------------------------------------------------
      // ABC2_mC8_12_c_a_i-001 (metal-oxide; Na1Ni1O2, ICSD #26609) == AB2C_mC8_12_a_i_c-001 (Na1O2V1, ICSD #420138) == ABC2_mC8_12_a_c_i-002 (Mn1Na1O2, ICSD #16270)
      // misfits=0.0746 and 0.0639
      // REASON FOR DUPLICATE: Improvement to XtalFinder code, reduces misfit
      if (aurostd::WithinList(labels_matched, "ABC2_mC8_12_c_a_i-001") && aurostd::WithinList(labels_matched, "AB2C_mC8_12_a_i_c-001") && aurostd::WithinList(labels_matched, "ABC2_mC8_12_a_c_i-002")) {
        return true;
      }
    }

    return false; // not a special case
  }
} // namespace anrl

// ***************************************************************************
// anrl::structureAndLabelConsistent()
// ***************************************************************************
namespace anrl {
  bool structureAndLabelConsistent(const xstructure& _xstr, const string& label_input, string& label_and_params_calculated,
                                   double tolerance_sym_input) { // DX20201105

    // Checks if the created structure is consistent with the label;
    // it is possible that the provided parameters elevate the structure
    // to a higher symmetry

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);

    xstructure xstr = _xstr; // copy

    // ---------------------------------------------------------------------------
    // set symmetry tolerance
    double tolerance_sym = tolerance_sym_input;
    if (tolerance_sym == AUROSTD_MAX_DOUBLE) {
      tolerance_sym = SYM::defaultTolerance(xstr);
    }

    // ---------------------------------------------------------------------------
    // determine label from structure (reverse process)
    label_and_params_calculated = structure2anrl(xstr, tolerance_sym); // DX20201105 - pass in symmetry tolerance

    // cannot do a strict string comparison of labels, symmetry analysis may
    // change origin (i.e., Wyckoff letters); need to check if labels are
    // isopointal (check SG and Wyckoff multiplicities and site symmetries)

    vector<string> label_fields;
    aurostd::string2tokens(label_input, label_fields, "_");

    // ---------------------------------------------------------------------------
    // check space groups
    const uint space_group_in_label = aurostd::string2utype<uint>(label_fields[2]);
    if (!compare::matchableSpaceGroups(xstr.space_group_ITC, space_group_in_label)) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " the calculated and label-designated space groups are incommensurate: "
             << "calculated=" << xstr.space_group_ITC << " vs "
             << "label= " << space_group_in_label << endl;
      }
      return false;
    }

    // ---------------------------------------------------------------------------
    // check Wyckoff positions

    // get Wyckoff information from label and format to compare
    const vector<vector<string>> Wyckoff_fields = compare::convertANRLWyckoffString2GroupedPositions(label_input);
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions_label;
    compare::groupWyckoffPositionsFromGroupedString(space_group_in_label, xstr.setting_ITC, Wyckoff_fields, grouped_Wyckoff_positions_label);

    // get Wyckoff information from xstructure and format to compare
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions_structure;
    compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions_structure);
    const string Wyckoff_string_structure = anrl::groupedWyckoffPosition2ANRLString(grouped_Wyckoff_positions_structure, true);

    // print grouped Wyckoff sequences
    if (LDEBUG) {
      print(grouped_Wyckoff_positions_label);
      cerr << "-------------------------" << endl;
      print(grouped_Wyckoff_positions_structure);
    }

    if (!compare::matchableWyckoffPositions(grouped_Wyckoff_positions_label, grouped_Wyckoff_positions_structure,
                                            false)) { // same_species=false - the structure MAY be decorated, but the label is NOT
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " the calculated and label-designated Wyckoff positions are incommensurate: "
             << "calculated=" << Wyckoff_string_structure << " vs "
             << "label= " << label_input << endl;
      }
      return false;
    }

    // ---------------------------------------------------------------------------
    // all tests passed; the structure and label are commensurate
    return true;
  }
} // namespace anrl

namespace anrl {
  std::string getPrototypeUID(const std::string& search_string, bool allow_legacy) {
    const bool LDEBUG = (false || XHOST.DEBUG);
    std::string uid;
    const ProtoData pd = ProtoData::get();
    std::string clean_search_string;
    // ensure extra information is cut away
    clean_search_string = search_string.substr(0, search_string.find(':'));

    // Prototype ID
    if (pd.content.find(clean_search_string) != pd.content.end()) {
      uid = clean_search_string;
    }
    // One to One aliases
    if (uid.empty()) {
      for (const std::string& id_type : {"label", "alias", "icsd", "ccdc"}) {
        if (pd.lookup[id_type].find(clean_search_string) != pd.lookup[id_type].end()) {
          uid = static_cast<string>(pd.lookup[id_type][clean_search_string]);
          break;
        }
      }
    }
    // Multi value aliases
    if (uid.empty()) {
      for (const std::string& id_type : {"prototype"}) {
        if (pd.lookup[id_type].find(clean_search_string) != pd.lookup[id_type].end()) {
          if (pd.lookup[id_type][clean_search_string].size() == 1) {
            uid = static_cast<string>(pd.lookup[id_type][clean_search_string][static_cast<size_t>(0)]);
          } else {
            cerr << "Found multiple possible protoypes:" << endl;
            for (const std::string& possible_uid : static_cast<vector<std::string>>(pd.lookup[id_type][clean_search_string])) {
              if (pd.content[possible_uid]["legacy"]) {
                cerr << "  " << possible_uid << "(legacy prototype " << static_cast<string>(pd.content[possible_uid]["label"]) << ")" << endl;
              } else {
                cerr << "  " << possible_uid << " | https://aflow.org/p/" << possible_uid << endl;
              }
            }
            stringstream message;
            message << "Could not find a unique prototype for the search string: '" << clean_search_string << "'";
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, clean_search_string, _INPUT_ILLEGAL_);
          }
          break;
        }
      }
    }
    if (LDEBUG) {
      if (uid.empty()) {
        cerr << "Search string (" << search_string << ") not found in the database" << endl;
      } else {
        cerr << "Found " << uid << " in the database";
        if (!pd.content[uid]["legacy"]) {
          cerr << " - https://aflow.org/p/" << uid;
        }
        cerr << endl;
      }
    }
    if (!uid.empty() && !allow_legacy && pd.content[uid]["legacy"]) {
      uid.clear();
      if (LDEBUG) {
        cerr << "Legacy prototypes not enabled, returning no result";
      }
    }
    return uid;
  }
} // namespace anrl

// ***************************************************************************
// anrl::PrototypeANRL_Generator()
// ***************************************************************************
// Returns a ANRL prototype structure based on the label and internal
// degrees of freedom.
// The function is generic and will build ANY prototype as long as:
//   1) the label and parameters are valid (the function has many checks) AND
//   2) the structure is a crystal (i.e., built from Wyckoff positions)
// A symbolic representations of the crystal can be returned in terms of:
// lattice variables: a, b, c, alpha, beta, gamma AND
// Wyckoff variables: x, y, and z
namespace anrl {

  xstructure PrototypeANRL_Generator_UID(const string& uid, string& parameters, deque<string>& vatomX, deque<double>& vvolumeX) {
    const anrl::ProtoData pd = anrl::ProtoData::get();
    std::string label = static_cast<string>(pd.content[uid]["label"]);
    if (parameters.empty()) {
      parameters = static_cast<string>(pd.content[uid]["parameter_values"]);
    }
    ofstream FileMESSAGE;
    return PrototypeANRL_Generator(label, parameters, vatomX, vvolumeX, FileMESSAGE, cout, true);
  }

  xstructure PrototypeANRL_Generator(string& label, string& parameters, deque<string>& vatomX, deque<double>& vvolumeX, ostream& logstream, bool silence_logger) {
    // command line version (no need for FileMESSAGE or logger)

    ofstream FileMESSAGE;

    const xstructure prototype = PrototypeANRL_Generator(label, parameters, vatomX, vvolumeX, FileMESSAGE, logstream, silence_logger);

    return prototype;
  }
} // namespace anrl

namespace anrl {
  xstructure PrototypeANRL_Generator(string& label, string& parameters, deque<string>& vatomX, deque<double>& vvolumeX, ofstream& FileMESSAGE, ostream& logstream, bool silence_logger) {
    // main version

    const bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_ANRL_);
    stringstream message;

    xstructure str;
    const ProtoData pd = ProtoData::get();
    aurostd::JSON::object loaded_prototype;

    // ---------------------------------------------------------------------------
    // determine print mode
    uint print_mode = _PROTO_GENERATOR_GEOMETRY_FILE_; // no equations
    if (XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
      print_mode = _PROTO_GENERATOR_EQUATIONS_ONLY_; // equations only
      str.symbolic_math_representation_only = true; // DX20180618 print symbolic math representation only
      message << "Printing the symbolic equations only";
    } else if (XHOST.vflag_pflow.flag("PROTO::ADD_EQUATIONS")) {
      print_mode = _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_; // equations + parameters
      str.constrained_symmetry_calculation = true; // DX20180618 appends information to geometry file for calculation
      message << "Printing the geometry file and the symbolic equations";
    } else {
      message << "Printing geometry file only";
    }
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_, silence_logger);

    // ---------------------------------------------------------------------------
    // declare variables

    vector<string> tokens;
    string label_anrl;
    string number_id = "001"; // for predefined anrls //DX20191207
    string label_permutations;
    deque<uint> vpermutation;

    // ---------------------------------------------------------------------------
    // handle corner cases //DX20200929
    if (label.find("sigma_tP30_136_bf2ij") != std::string::npos) {
      aurostd::StringSubstInPlace(label, "sigma_tP30_136_bf2ij", "A_tP30_136_bf2ij"); // label
      aurostd::StringSubstInPlace(label, ".sigma", ".A"); // permutation
    }

    // ---------------------------------------------------------------------------
    // search for label_permutations
    aurostd::string2tokens(label, tokens, ".");
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << ": tokens.size()=" << tokens.size() << endl;
    }

    if (tokens.empty()) {
      label_anrl = label;
    }
    if (tokens.size() == 1) {
      label_anrl = tokens.at(0);
    }
    if (tokens.size() == 2) {
      label_anrl = tokens.at(0);
      label_permutations = tokens.at(1);
    }
    for (size_t i = 0; i < label_permutations.size(); i++) {
      vpermutation.push_back(aurostd::mod(label_permutations.at(i) - 65, 32));
    }

    // ---------------------------------------------------------------------------
    // check if preset suffix is included with label, e.g., A_hR2_166_c-001
    if (label_anrl.find('-') != std::string::npos) {
      tokens.clear();
      aurostd::string2tokens(label_anrl, tokens, "-");
      if (tokens.size() == 2) {
        label_anrl = tokens[0];
        number_id = tokens[1];
      }
    }
    message << "The input label is " << label_anrl;
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_, silence_logger);

    if (!number_id.empty()) {
      message << "The preset parameters " << number_id << " will be extracted (if they exist)";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_, silence_logger);
    }

    const std::string search_label = label_anrl + "-" + number_id;
    for (const std::string& id_type : {"label", "alias", "icsd", "prototype"}) {
      auto result = pd.lookup[id_type].find(search_label);
      if (result != pd.lookup[id_type].end()) {
        loaded_prototype = pd.content[static_cast<std::string>(result->second)];
        break;
      }
    }

    // ---------------------------------------------------------------------------
    // not found, new label
    if (!static_cast<bool>(loaded_prototype) and !XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
      message << "This label (" << label_anrl << ") does not currently exist in the AFLOW library. Consider adding it to the library.";
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_, silence_logger);
    }

    // -------------------------------------------------------------------------
    // check if using original anrl lattice parameter value when using the
    // preset parameter functionality
    bool keep_anrl_lattice_parameter = false;
    bool scale_volume_by_species = false; // DX20201104 - default should be false
    if (parameters == "use_anrl_lattice_param") {
      keep_anrl_lattice_parameter = true;
      scale_volume_by_species = false;
      parameters = ""; // clear the hack
    }

    // -------------------------------------------------------------------------
    // if no parameters given
    if (parameters.empty()) {
      // -------------------------------------------------------------------------
      // extract values from library
      if (static_cast<bool>(loaded_prototype)) {
        std::vector<double> parameter_values = loaded_prototype["parameter_values"];
        if (!keep_anrl_lattice_parameter) {
          parameter_values[0] = -1;
        }
        parameters = aurostd::joinWDelimiter(parameter_values, ",");
      }
      message << "Extracted the following parameters (internal degrees of freedom) from the AFLOW parameters= " << parameters;
      pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_, silence_logger);
    }

    // ---------------------------------------------------------------------------
    // split label into fields
    aurostd::string2tokens(label_anrl, tokens, "_");

    message << "The AFLOW label has been partitioned into " << tokens.size() << " fields : " << aurostd::joinWDelimiter(tokens, " ");
    pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_, silence_logger);

    if (tokens.size() < 4) {
      message << "Number of fields in label is too small, should be 4 or more: label" << label_anrl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // check stoichometry and get species
    const string compound_string = tokens[0];
    vector<uint> stoichiometry;
    vector<string> species = aurostd::getElements(compound_string, stoichiometry); // DX20200724
    vector<uint> reduced_stoich;
    aurostd::reduceByGCD(stoichiometry, reduced_stoich);
    if (!compare::sameStoichiometry(stoichiometry, reduced_stoich)) {
      message << "The input stoichiometry (first field in label=" << compound_string << ") is not reduced, it should be: " << aurostd::joinWDelimiter(reduced_stoich, ":");
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    for (size_t i = 0; i < species.size(); i++) { // number of species
      str.num_each_type.push_back(0);
      str.comp_each_type.push_back(0.0);
      str.species.emplace_back("");
      str.species_pp.emplace_back("");
      str.species_pp_type.emplace_back("");
      str.species_pp_version.emplace_back("");
      str.species_pp_ZVAL.push_back(0.0);
      str.species_pp_vLDAU.emplace_back();
      str.species_volume.push_back(0.0);
      str.species_mass.push_back(0.0);
    }

    // ---------------------------------------------------------------------------
    // get Pearson symbol
    string Pearson_symbol = tokens[1];
    char lattice_type = Pearson_symbol[0];
    char lattice_centering = Pearson_symbol[1];
    const uint number_of_atoms_conventional = aurostd::string2utype<uint>(Pearson_symbol.substr(2, Pearson_symbol.size()));

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " # atoms in conventional cell (from Pearson): " << number_of_atoms_conventional << endl;
    }

    // ---------------------------------------------------------------------------
    // get space group number
    const uint space_group_number = aurostd::string2utype<uint>(tokens[2]);

    if (space_group_number < 1 || space_group_number > 230) {
      message << "The space group number is invalid; it must be between 1-230: spacegroup=" << space_group_number;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // check if Pearson symbol and space group match
    const string lattice_and_centering_from_Pearson = Pearson_symbol.substr(0, 2);
    const string lattice_and_centering_from_sg = SYM::spacegroup2latticeAndCentering(space_group_number);
    if (lattice_and_centering_from_Pearson != lattice_and_centering_from_sg) {
      message << "Pearson symbol and space group number are incommensurate; the lattice centerings do not match:";
      message << "Pearson=" << Pearson_symbol << " (centering=" << lattice_and_centering_from_Pearson << ") vs ";
      message << "SG=" << space_group_number << "(centering=" << lattice_and_centering_from_sg << ")";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // get space group information
    string space_group_symbol = GetSpaceGroupName(space_group_number);
    const char space_group_letter = space_group_symbol[0];

    // ---------------------------------------------------------------------------
    // get Wyckoff positions (the remaining fields after the first three)
    vector<string> Wyckoff_tokens;
    Wyckoff_tokens.insert(Wyckoff_tokens.end(), tokens.begin() + 3, tokens.end()); // get Wyckoff tokens

    // ---------------------------------------------------------------------------
    // check if number of Wyckoff positions match the number of species
    if (Wyckoff_tokens.size() != species.size()) {
      message << "The number of species does not match the number of Wyckoff species: # species=" << species.size() << " vs ";
      message << "# Wyckoff species=" << Wyckoff_tokens.size() << " (input label=" << label << ")" << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " the Wyckoff sequences associated with each species are:" << endl;
      for (size_t i = 0; i < species.size(); i++) {
        cerr << species[i] << ": " << Wyckoff_tokens[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // initialize ITC space group/Wyckoff position object
    const uint setting = SG_SETTING_ANRL;

    vector<wyckoffsite_ITC> wyckoff_sites_ITC = anrl::getWyckoffSitesFromANRL(Wyckoff_tokens, species, space_group_number, setting);
    const vector<uint> number_of_each_type = SYM::numberEachTypeFromWyckoff(wyckoff_sites_ITC);

    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " the Wyckoff positions are:" << endl;
      print(wyckoff_sites_ITC);
    }

    // ---------------------------------------------------------------------------
    // check Wyckoff positions and stoichiometry
    // vector<uint> Wyckoff_reduced_stoich = compare::gcdStoich(number_of_each_type);
    vector<uint> Wyckoff_reduced_stoich;
    aurostd::reduceByGCD(number_of_each_type, Wyckoff_reduced_stoich);
    if (!compare::sameStoichiometry(stoichiometry, Wyckoff_reduced_stoich)) {
      message << "The input composition and Wyckoff positions yield different stoichiometries: composition=" << aurostd::joinWDelimiter(stoichiometry, ":");
      message << ", Wyckoff=" << aurostd::joinWDelimiter(Wyckoff_reduced_stoich, ":") << endl;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // check Wyckoff multiplicity and number of atoms in the conventional cell
    // (should make this into a function)
    uint Wyckoff_multiplicity_sum = 0;
    for (size_t i = 0; i < wyckoff_sites_ITC.size(); i++) {
      Wyckoff_multiplicity_sum += wyckoff_sites_ITC[i].multiplicity;
    }
    if (Wyckoff_multiplicity_sum != number_of_atoms_conventional) {
      message << "The sum of the Wyckoff multiplicity does not add up to the number of atoms in the conventional cell (from Pearson symbol); bad prototype label: ";
      message << " Wyckoff multiplicity sum = " << Wyckoff_multiplicity_sum;
      message << ", Pearson symbol = " << Pearson_symbol;
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // determine parameters
    vector<string> parameter_list;
    vector<string> lattice_parameter_list;
    vector<string> Wyckoff_parameter_list;
    vector<double> lattice_parameter_values;
    vector<double> Wyckoff_parameter_values;

    // -------------------------------------------------------------------------
    // lattice parameters
    lattice_parameter_list = getANRLLatticeParameterString(lattice_type);

    // ---------------------------------------------------------------------------
    // reorder Wyckoff positions alphabetically by Wyckoff letter, then by species
    vector<wyckoffsite_ITC> ordered_Wyckoff_sites_ITC = wyckoff_sites_ITC;
    std::sort(ordered_Wyckoff_sites_ITC.begin(), ordered_Wyckoff_sites_ITC.end(), sortWyckoffByLetter);

    if (LDEBUG) {
      for (size_t i = 0; i < ordered_Wyckoff_sites_ITC.size(); i++) {
        cerr << __AFLOW_FUNC__ << ":Ordered Wyckoff site: " << ordered_Wyckoff_sites_ITC[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // determine degrees of freedom in Wyckoff positions
    Wyckoff_parameter_list = determineWyckoffVariables(ordered_Wyckoff_sites_ITC);

    // ---------------------------------------------------------------------------
    // combine parameter vectors
    parameter_list = lattice_parameter_list;
    parameter_list.insert(parameter_list.end(), Wyckoff_parameter_list.begin(), Wyckoff_parameter_list.end());
    if (XHOST.vflag_pflow.flag("PROTO::SYMBOLS_ONLY")) {
      cout << aurostd::joinWDelimiter(parameter_list, ",") << endl;
      return str;
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " parameters=" << aurostd::joinWDelimiter(parameter_list, ",") << endl;
    }

    // -------------------------------------------------------------------------
    // if no parameters are provided and more than one parameter is needed,
    // we throw an error; this is a new prototype
    if (!XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
      if (parameters.empty() && parameter_list.size() != 1) {
        message << "No parameters provided. Since this is a new prototype label with more than one degree of freedom,";
        message << "you must add parameter values with --params=..." << endl;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      // if only one parameter is needed, we can generate the structure,
      // i.e., there are no degrees of freedom (other than the lattice parameter)
      // for this label, and it does not require an enumeration suffix
      else if (parameters.empty() && parameter_list.size() == 1) {
        parameters = "1.0";
      }

      // ---------------------------------------------------------------------------
      // partition in parameter values
      vector<string> vparameters_temp;
      aurostd::string2tokens(parameters, vparameters_temp, ",");
      if (aurostd::string2utype<double>(vparameters_temp[0]) < _ZERO_TOL_) {
        // DX20201104 - was missing
        scale_volume_by_species = true;
        vparameters_temp[0] = "1.0"; // fix
        parameters = aurostd::joinWDelimiter(vparameters_temp, ",");
      }
      vector<double> vparameters = aurostd::vectorstring2vectorutype<double>(vparameters_temp);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " parameter_values=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vparameters, AUROSTD_DEFAULT_PRECISION, FIXED_STREAM), ",") << endl;
      }

      // ---------------------------------------------------------------------------
      // check for automatic volume scaling (i.e., first parameter is negative)
      if (vparameters[0] <= 0.0) {
        // CO20181226 forget signbit, also include 0
        vparameters[0] = 1.0; // fix
        vparameters_temp[0] = "1.0"; // fix
        parameters = aurostd::joinWDelimiter(vparameters_temp, ",");
      }

      if (vparameters.size() != parameter_list.size()) {
        message << "The number of input parameters does not match the number required by the lattice and Wyckoff positions: ";
        message << "input parameters=" << parameters << " vs ";
        message << "parameter_list=" << aurostd::joinWDelimiter(parameter_list, ",") << endl;
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_NUMBER_);
      }

      // ---------------------------------------------------------------------------
      // populate degree of freedom values
      lattice_parameter_values.insert(lattice_parameter_values.end(), vparameters.begin(), vparameters.begin() + lattice_parameter_list.size());
      Wyckoff_parameter_values.insert(Wyckoff_parameter_values.end(), vparameters.begin() + lattice_parameter_list.size(), vparameters.end());

      anrl::applyWyckoffValues(Wyckoff_parameter_values, ordered_Wyckoff_sites_ITC);

      // ---------------------------------------------------------------------------
      // check to ensure no duplicate
      if (containsDuplicateWyckoffCoordinate(ordered_Wyckoff_sites_ITC, true)) {
        message << "Contains duplicate Wyckoff letters with the same degrees of freedom.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }

      // ---------------------------------------------------------------------------
      // generate lattice (in ANRL convention) based on symmetry and lattice
      // parameter values

      // primitive (0)
      const xmatrix<double> lattice_primitive = getLattice(lattice_and_centering_from_Pearson, space_group_letter, lattice_parameter_values, 0);
      // conventional (1)
      const xmatrix<double> lattice_conventional = getLattice(lattice_and_centering_from_Pearson, space_group_letter, lattice_parameter_values, 1);

      // ---------------------------------------------------------------------------
      // generate atoms based Wyckoff equations and Wyckoff parameter values
      const deque<_atom> atoms_conventional_cell = getAtomsFromWyckoff(ordered_Wyckoff_sites_ITC, lattice_conventional);

      // ---------------------------------------------------------------------------
      // get interatomic distance to find a good "fold-in" tolerance //DX20201021
      xstructure str_conv;
      str_conv.lattice = lattice_conventional;
      str_conv.atoms = atoms_conventional_cell;
      str_conv.sym_eps = SYM::defaultTolerance(str_conv);

      deque<_atom> atoms_primitive_cell;
      // special case: if using the rhombohedral setting, then the Wyckoff positions
      // are already given wrt to the primitive cell; no need to perform conversion
      if (setting == SG_SETTING_ANRL && lattice_centering == 'R') {
        atoms_primitive_cell = atoms_conventional_cell;
        // need to sort atoms alphabetically
        std::stable_sort(atoms_primitive_cell.begin(), atoms_primitive_cell.end(), sortAtomsNames);
      }
      // generic case: convert conventional to primitive
      else {
        atoms_primitive_cell = foldAtomsInCell(atoms_conventional_cell, lattice_conventional, lattice_primitive, false,
                                               str_conv.sym_eps, // DX20201019 - use sym_eps instead of 1e-6
                                               false);
      }

      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " atoms in the primitive cell (" << atoms_primitive_cell.size() << "):" << endl;
        for (size_t i = 0; i < atoms_primitive_cell.size(); i++) {
          cerr << atoms_primitive_cell[i] << " " << atoms_primitive_cell[i].name << endl;
        }
      }

      // ---------------------------------------------------------------------------
      // check ratio between conventional and primitive atoms
      if (atoms_conventional_cell.size() % atoms_primitive_cell.size() != 0) {
        message << "The ratio of atoms between the conventional cell and primitive cell is not an integer; check the tolerance: #conventional=" << atoms_conventional_cell.size()
                << " vs #primitive=" << atoms_primitive_cell.size();
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      const uint ratio_calculated = atoms_conventional_cell.size() / atoms_primitive_cell.size();
      const uint ratio_conventional2primitive = LATTICE::Conventional2PrimitiveRatio(lattice_centering);
      if (!aurostd::isequal(ratio_calculated, ratio_conventional2primitive)) {
        if (!(ratio_calculated == 1 && lattice_centering == 'R' && setting == SG_SETTING_ANRL)) {
          message << "The calculated ratio and the expected conventional2primtive ratio do not match: calculated=" << ratio_calculated << " vs expected=" << ratio_conventional2primitive;
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
      }

      // ---------------------------------------------------------------------------
      // create xstructure
      str.iomode = IOVASP_AUTO;
      str.title = label + " params=" + parameters + " SG=" + aurostd::utype2string(space_group_number) + DOI_ANRL;
      // CO20190520
      str.scale = 1.0;
      str.lattice = lattice_primitive;
      str.atoms = atoms_primitive_cell;
      str.dist_nn_min = SYM::minimumDistance(str);
      // DX20210114 - calculate dist_nn_min before default tolerance, defaultTolerance will use dist_nn_min
      str.sym_eps = SYM::defaultTolerance(str);
      // need sym_eps for AddAtom later (otherwise it breaks for systems like A12B6C_cF608_210_4h_2h_e)
      str.sym_eps_calculated = true; // DX20200929

      // ---------------------------------------------------------------------------
      // add ANRL info to xstructure
      str.num_parameters = parameter_list.size();
      str.num_lattice_parameters = lattice_parameter_list.size();
      str.prototype_parameter_list = parameter_list;
      str.prototype_parameter_values = {};
      str.setting_ITC = setting;

      // ---------------------------------------------------------------------------
      // convert RHL to HEX setting ([--hex] option)
      if (XHOST.vflag_pflow.flag("PROTO::HEX") && lattice_centering == 'R') {
        vector<double> vparameters;
        aurostd::string2tokens(parameters, vparameters, ",");
        uint i = 0;
        double a = vparameters.at(i++);
        const double covera = vparameters.at(i++);
        double c = covera * a;
        str = rhl2hex(str, a, c);
      }

      for (size_t iat = 0; iat < str.atoms.size(); iat++) {
        str.atoms[iat].name_is_given = true;
        //[CO20200130 // reference position for convasp
        str.atoms[iat].basis = iat; // position in the basis
        if (print_mode != _PROTO_GENERATOR_EQUATIONS_ONLY_) {
          // equations only //DX20180618
          str.atoms[iat].cpos = F2C(str.lattice, str.atoms[iat].fpos);
        }
        str.num_each_type.at(str.atoms[iat].type)++;
        //     str.comp_each_type.at(str.atoms[iat].type)+=1.0; inside code
        str.species.at(str.atoms[iat].type) = str.atoms[iat].name;
      }
    }

    // ---------------------------------------------------------------------------
    // symbolic representation of prototypes
    if (print_mode == _PROTO_GENERATOR_EQUATIONS_ONLY_ || print_mode == _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_) {
#if USE_SYMBOLIC_SOURCE // DX20200831 - defined in aflow.h
      // ---------------------------------------------------------------------------
      // get symbolic lattice
      const symbolic::Symbolic lattice_symbolic = SymbolicANRLPrimitiveLattices(lattice_and_centering_from_Pearson, space_group_letter);

      // ---------------------------------------------------------------------------
      // order alphabetically by species //DX20210217
      vector<wyckoffsite_ITC> Wyckoff_sites_alphabetic = ordered_Wyckoff_sites_ITC;
      std::sort(Wyckoff_sites_alphabetic.begin(), Wyckoff_sites_alphabetic.end(), sortWyckoffByType);

      // ---------------------------------------------------------------------------
      // convert Wyckoff site into symbolic notation
      vector<SymbolicWyckoffSite> Wyckoff_sites_symbolic;
      for (size_t i = 0; i < ordered_Wyckoff_sites_ITC.size(); i++) {
        Wyckoff_sites_symbolic.push_back(initializeSymbolicWyckoffSite(Wyckoff_sites_alphabetic[i]));
      }

      // ---------------------------------------------------------------------------
      // transform to ANRL primitive cell
      for (size_t i = 0; i < Wyckoff_sites_symbolic.size(); i++) {
        Wyckoff_sites_symbolic[i].equations = convertEquations2FractionalEquations(lattice_and_centering_from_Pearson, lattice_symbolic, Wyckoff_sites_symbolic[i].equations);
      }

      // ---------------------------------------------------------------------------
      // convert generic variable to the parameter designation, e.g., x -> x2
      substituteVariableWithParameterDesignation(Wyckoff_sites_symbolic);
      vector<symbolic::Symbolic> symbolic_equations;
      for (size_t i = 0; i < Wyckoff_sites_symbolic.size(); i++) {
        symbolic_equations.insert(symbolic_equations.end(), Wyckoff_sites_symbolic[i].equations.begin(), Wyckoff_sites_symbolic[i].equations.end());
      }

      // This shouldn't affect anything except web page generation.
      // It may add some extraneous output in some cases.
      // This outputs the basis equations and wyckoff positions for use in
      // generating the web pages: Originally AZ2023, Added to AFLOW4 by MJM20250521
      str.symbolic_math_lattice = symbolic::matrix2VectorVectorString(lattice_symbolic);
      if (XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
        char parmeq[32];
        char parteq[16];
        cout << "BEGIN EQUATIONS WITH WYCKOFF" << endl;
        for (uint i = 0; i < Wyckoff_sites_symbolic.size(); i++) {
          // We don't know how many atom types there will be, so format big: MJM 20230131
          snprintf(parmeq, 32, "%4.4i", Wyckoff_sites_symbolic[i].parameter_index);
          for (uint j = 0; j < Wyckoff_sites_symbolic[i].equations.size(); j++) {
            snprintf(parteq, 16, "%2.2i", j);
            // <em>Try</em> to put the Wyckoff symbol at the front as part of my index
            // <em>Maybe</em> this will help MJM 230201
            cout << Wyckoff_sites_symbolic[i].letter << "." << parmeq << "." << parteq << " " << Wyckoff_sites_symbolic[i].equations[j].row(0) << " " << Wyckoff_sites_symbolic[i].equations[j].row(1) << " "
                 << Wyckoff_sites_symbolic[i].equations[j].row(2) << " " << Wyckoff_sites_symbolic[i].type << " " << Wyckoff_sites_symbolic[i].multiplicity << " " << Wyckoff_sites_symbolic[i].letter << endl;
          }
        }
        cout << "END EQUATIONS WITH WYCKOFF" << endl;

        cout << "BEGIN SYMBOLIC LATTICE" << endl;
        for (const vector<std::string>& row : str.symbolic_math_lattice) {
          cout << aurostd::joinWDelimiter(row, "  ") << endl;
        }
        cout << "END SYMBOLIC LATTICE" << endl;
        return str;
      }
      // ---------------------------------------------------------------------------
      // convert to vector<string> and add to _atom
      if (!XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
        addSymbolicEquation2Atoms(symbolic_equations, str.atoms);
      }
#else
      // ---------------------------------------------------------------------------
      // if the SymbolicC++ code is not compiled
      message << "The SymbolicC++ source code has not been compiled, symbolic equations cannot be printed";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
#endif
    }

    // ---------------------------------------------------------------------------
    // DONE
    if (print_mode != _PROTO_GENERATOR_EQUATIONS_ONLY_) { // equations only //DX20180618
      xvector<double> data(6);
      data = Getabc_angles(str.lattice, DEGREES);
      str.a = data[1];
      str.b = data[2];
      str.c = data[3];
      str.alpha = data[4];
      str.beta = data[5];
      str.gamma = data[6];
      clear(str.origin);
    }
    //  if(vpflow.flag("STDPRIMCELL")) {cout << "EUREKA"<< endl;} //cout << GetStandardPrimitive(xstructure(cin,IOAFLOW_AUTO));

    // ---------------------------------------------------------------------------
    // NOW PLAY WITH PERMUTATIONS and ATOMX
    if (!vpermutation.empty() || !vatomX.empty()) {
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " PERMUTATIONS" << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vpermutation.size()=" << vpermutation.size() << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vpermutation =";
        for (size_t i = 0; i < vpermutation.size(); i++) {
          cerr << " " << vpermutation.at(i);
        }
        cerr << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " ATOMX" << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vatomX.size()=" << vatomX.size() << endl;
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " vatomX =";
        for (size_t i = 0; i < vatomX.size(); i++) {
          cerr << " " << vatomX.at(i);
        }
        cerr << endl;
      }
      if (print_mode != _PROTO_GENERATOR_EQUATIONS_ONLY_) { // equations only //DX20180618
        std::deque<_atom> atoms;
        atoms = str.atoms;
        // STRIP ALL ATOMS
        while (!str.atoms.empty()) {
          str.RemoveAtom(0);
        }
        // ADD MODIFIED ATOMS
        for (size_t i = 0; i < atoms.size(); i++) {
          const uint type = atoms[i].type;
          if (!vpermutation.empty()) {
            atoms[i].type = vpermutation.at(type);
          } // PERMUTATIONS
          if (!vpermutation.empty() || !vatomX.empty()) {
            atoms[i].name = vatomX.at(atoms[i].type);
          } // PERMUTATIONS AND ATOMX
          //	atoms[i].name=aurostd::mod(label_permutations.at(type)-65,32)+65;
          str.AddAtom(atoms[i], false); // CO20230319 - add by type
          // DX20181205 - Volume scaling by atomic species - START
          // ---------------------------------------------------------------------------
          // if a=1.0 for prototype (i.e., no scaling factor), use atomic species to get volume
          if (scale_volume_by_species == true) {
            double volume = 0.0;
            for (size_t i = 0; i < str.num_each_type.size(); i++) {
              for (uint j = 0; j < (uint) str.num_each_type[i]; j++) {
                volume += vvolumeX[i];
                if (LDEBUG) {
                  cerr << __AFLOW_FUNC__ << " volume=" << volume << "  (" << vvolumeX[i] << ")" << endl;
                }
              }
            }
            str.SetVolume(volume); // CO20190205 - more robust
            str.neg_scale = true;
          }
          // DX20181205 - Volume scaling by atomic species - END
        }
      }
      str.SpeciesPutAlphabetic();
    }

    // ---------------------------------------------------------------------------
    // fix title of geometry file
    aurostd::StringSubstInPlace(str.title, label, aurostd::joinWDelimiter(str.species, "") + "/" + label); // vlabel.at(ifound) //use label as we want permutations too //CO20181216
    if (scale_volume_by_species) {
      aurostd::StringSubstInPlace(str.title, " params=1.0", " params=-1");
    } // CO20181216

    // ---------------------------------------------------------------------------
    // if this is a new prototype (i.e., not in library), we should check the
    // symmetry; it is possible that the provided parameters elevate the structure
    // to a higher symmetry
    //
    if (!XHOST.vflag_pflow.flag("PROTO::WEBPAGE")) {
      // this section is skipped when we are generating web pages //MJM20250521 //HE20250524
      if (static_cast<bool>(loaded_prototype)) {
        string updated_label_and_params;
        if (!structureAndLabelConsistent(str, label_anrl, updated_label_and_params)) {
          // if changes symmetry, give the appropriate label
          message << "The structure has a higher symmetry than indicated by the label. To ignore this error use --webpage\n";
          message << "The correct label and parameters for this structure are:" << endl;
          message << updated_label_and_params << endl;
          message << "Please feed this label and set of parameters into the prototype generator.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
      }
    }
    if (number_of_atoms_conventional != str.atoms.size() && XHOST.vflag_pflow.flag("PROTO::STRICT")) {
      message << "The chosen parameters result in overlapping atoms\n";
      message << "Expected " << number_of_atoms_conventional << " atoms but only " << str.atoms.size() << " atoms in final structure after cleaning.";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return str;
  }
} // namespace anrl

#endif // _AFLOW_ANRL_CPP

// Written by Stefano Curtarolo - 2016
// Written by David Hicks (DX) - 2016/2020 (generic prototype generator)
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *           Aflow DAVID HICKS - Duke University 2014-2021                 *
// *                                                                         *
// ***************************************************************************
